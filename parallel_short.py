#!/usr/bin/env python3
"""Run independent short.gp searches in parallel and verify the first result."""

import argparse
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time


ROOT = Path(__file__).resolve().parent
SHORT_GP = ROOT / "short.gp"
VERIFIER = ROOT / "vshortECPP.py"


def positive_int(value):
    value = int(value)
    if value < 1:
        raise argparse.ArgumentTypeError("must be positive")
    return value


def positive_int_list(value):
    try:
        values = [positive_int(item) for item in value.split(",")]
    except (ValueError, argparse.ArgumentTypeError) as error:
        raise argparse.ArgumentTypeError("must be comma-separated positive integers") from error
    if not values:
        raise argparse.ArgumentTypeError("must not be empty")
    return values


def nonnegative_int(value):
    value = int(value)
    if value < 0:
        raise argparse.ArgumentTypeError("must be nonnegative")
    return value


def nonnegative_int_list(value):
    try:
        values = [nonnegative_int(item) for item in value.split(",")]
    except (ValueError, argparse.ArgumentTypeError) as error:
        raise argparse.ArgumentTypeError(
            "must be comma-separated nonnegative integers"
        ) from error
    if not values:
        raise argparse.ArgumentTypeError("must not be empty")
    return values


def parse_args():
    parser = argparse.ArgumentParser(
        description="Search for a short ECPP with parallel PARI/GP workers."
    )
    parser.add_argument("p", type=positive_int, help="probable prime to certify")
    parser.add_argument(
        "-j",
        "--workers",
        type=positive_int,
        default=max(1, os.cpu_count() or 1),
        help="number of GP workers (default: number of logical CPUs)",
    )
    parser.add_argument(
        "--gp-threads",
        type=positive_int,
        default=1,
        help="PARI threads allowed per worker (default: 1)",
    )
    parser.add_argument(
        "--factor-seconds",
        type=positive_int_list,
        default=[20, 20, 20, 2, 5],
        help=(
            "comma-separated per-order factorization limits assigned cyclically "
            "to workers (default: 20,20,20,2,5)"
        ),
    )
    parser.add_argument(
        "--seed",
        type=positive_int,
        default=1,
        help="base random seed; worker i uses seed+i (default: 1)",
    )
    parser.add_argument(
        "-o", "--output", type=Path, help="write the verified certificate to this file"
    )
    return parser.parse_args()


def worker_input(p, factor_seconds, seed, gp_threads):
    return (
        f"default(nbthreads,{gp_threads}); SC_tlim={factor_seconds}; "
        f"setrand({seed}); c=shortcert({p}); "
        'for(i=1,#c,printf("%d%s",c[i],if(i<#c," ","\\n"))); '
        'warning("curves=",SC_curves);\n'
    )


def stop_workers(workers):
    for process, _ in workers:
        if process.poll() is None:
            process.terminate()
    deadline = time.monotonic() + 5
    for process, _ in workers:
        if process.poll() is None:
            try:
                process.wait(timeout=max(0, deadline - time.monotonic()))
            except subprocess.TimeoutExpired:
                process.kill()
    for process, _ in workers:
        process.wait()


def verify(certificate):
    fields = certificate.split()
    if not fields or any(not field.isdigit() for field in fields):
        return False
    result = subprocess.run(
        [sys.executable, str(VERIFIER), *fields],
        cwd=ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    return result.returncode == 0 and result.stdout.strip() == "True"


def main():
    args = parse_args()
    if not SHORT_GP.exists() or not VERIFIER.exists():
        raise SystemExit("short.gp and vshortECPP.py must be beside this script")

    workers = []
    started = time.monotonic()
    try:
        for index in range(args.workers):
            factor_seconds = args.factor_seconds[index % len(args.factor_seconds)]
            log = tempfile.TemporaryFile(mode="w+")
            process = subprocess.Popen(
                ["gp", "-q", str(SHORT_GP)],
                cwd=ROOT,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=log,
                text=True,
            )
            process.stdin.write(
                worker_input(
                    args.p,
                    factor_seconds,
                    args.seed + index,
                    args.gp_threads,
                )
            )
            process.stdin.close()
            workers.append((process, log))
        print(
            f"started {args.workers} workers for p={args.p} "
            f"with seeds {args.seed}..{args.seed + args.workers - 1} "
            f"and factor limits "
            f"{[args.factor_seconds[i % len(args.factor_seconds)] for i in range(args.workers)]}; "
            f"{args.gp_threads} PARI thread(s) each",
            file=sys.stderr,
            flush=True,
        )

        last_report = started
        while True:
            for index, (process, log) in enumerate(workers):
                returncode = process.poll()
                if returncode is None:
                    continue
                output = process.stdout.read().strip()
                if returncode == 0 and verify(output):
                    elapsed = time.monotonic() - started
                    log.seek(0)
                    worker_detail = log.read().strip().splitlines()
                    stop_workers(workers)
                    if args.output:
                        args.output.write_text(output + "\n")
                    print(output)
                    print(
                        f"worker {index} found and verified a certificate "
                        f"in {elapsed:.1f} seconds"
                        + (f" ({worker_detail[-1]})" if worker_detail else ""),
                        file=sys.stderr,
                    )
                    return 0
                log.seek(0)
                detail = log.read().strip()
                raise RuntimeError(
                    f"worker {index} exited without a valid certificate "
                    f"(status {returncode}): {detail or output or 'no output'}"
                )

            now = time.monotonic()
            if now - last_report >= 60:
                print(f"searching: {now - started:.0f} seconds elapsed", file=sys.stderr)
                last_report = now
            time.sleep(0.25)
    except (KeyboardInterrupt, Exception) as error:
        stop_workers(workers)
        if isinstance(error, KeyboardInterrupt):
            print("search interrupted", file=sys.stderr)
            return 130
        raise
    finally:
        for _, log in workers:
            log.close()


if __name__ == "__main__":
    signal.signal(signal.SIGTERM, lambda *_: (_ for _ in ()).throw(KeyboardInterrupt()))
    raise SystemExit(main())
