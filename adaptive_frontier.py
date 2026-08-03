#!/usr/bin/env python3
"""Iterate bounded saved-order factoring while the exact frontier improves."""

import argparse
import math
import subprocess
import sys
from pathlib import Path

from parallel_short import ROOT, load_resume_candidates, verify


DRIVER = ROOT / "parallel_short.py"


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Repeatedly factor the ranked saved-order frontier, reloading exact "
            "residual progress after every bounded pass."
        )
    )
    parser.add_argument("p", type=int, help="probable prime to certify")
    parser.add_argument("--candidate-dir", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("-o", "--output", type=Path, required=True)
    parser.add_argument("-j", "--workers", type=int, default=4)
    parser.add_argument("--frontier-size", type=int, default=40)
    parser.add_argument(
        "--frontier-offset",
        type=int,
        default=0,
        help="skip this many exactly ranked checkpoints before each pass",
    )
    parser.add_argument("--max-rounds", type=int, default=4)
    parser.add_argument(
        "--stale-rounds",
        type=int,
        default=1,
        help="stop after this many consecutive passes without exact progress",
    )
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--factor-seconds", type=int, default=5)
    parser.add_argument("--factor-flags", type=int, default=9)
    parser.add_argument("--prefactor-seconds", type=int, default=2)
    parser.add_argument("--pm1-bound", type=int, default=5_000_000)
    parser.add_argument("--pp1-bound", type=int, default=5_000_000)
    parser.add_argument("--ecm-bound", type=int, default=100_000)
    parser.add_argument("--ecm-curves", type=int, default=8)
    parser.add_argument("--branch-curves", type=int, default=96)
    args = parser.parse_args()
    for name in ("workers", "frontier_size", "max_rounds", "stale_rounds"):
        if getattr(args, name) < 1:
            parser.error(f"--{name.replace('_', '-')} must be positive")
    if args.frontier_offset < 0:
        parser.error("--frontier-offset must be nonnegative")
    return args


def candidate_key(candidate):
    if candidate["kind"] == "level":
        return tuple(
            candidate.get(field, "") for field in ("kind", "A", "x", "o", "q")
        )
    if candidate["kind"] == "cm_screen":
        return tuple(
            candidate.get(field, "") for field in ("kind", "D", "N", "o", "q")
        )
    return (
        candidate["kind"],
        candidate.get("A", ""),
        candidate.get("xden", ""),
        candidate.get("N", ""),
    )


def snapshot(candidates, size):
    return tuple(
        (candidate_key(candidate), candidate.get("residual", ""))
        for candidate in candidates[:size]
    )


def describe(candidate):
    return (
        f"kind={candidate['kind']} "
        f"D={candidate.get('A', candidate.get('D', '-'))} "
        f"smooth={candidate['smooth_bits']} bits "
        f"residual={candidate['residual_bits']} bits "
        f"minimum-q={candidate['minimum_q_bits']} bits "
        f"gap={candidate['cofactor_gap_bits']} bits"
    )


def verified_for(certificate, p):
    fields = certificate.split()
    return bool(fields) and fields[0] == str(p) and verify(certificate)


def main():
    args = parse_args()
    candidate_dir = args.candidate_dir.resolve()
    manifest = args.manifest.resolve()
    output = args.output.resolve()
    if output.exists():
        certificate = output.read_text(encoding="utf-8").strip()
        if verified_for(certificate, args.p):
            print(certificate)
            print(f"using existing verified certificate {output}", file=sys.stderr)
            return 0
        raise SystemExit(f"refusing to overwrite non-certificate output {output}")

    stale = 0
    for round_index in range(args.max_rounds):
        before = load_resume_candidates(candidate_dir, args.p)
        before_frontier = before[
            args.frontier_offset : args.frontier_offset + args.frontier_size
        ]
        if not before_frontier:
            raise SystemExit("frontier offset is beyond the ranked checkpoint queue")
        before_snapshot = snapshot(before_frontier, args.frontier_size)
        print(
            f"frontier round {round_index + 1}/{args.max_rounds}: "
            f"{len(before_snapshot)} checkpoint(s) after offset "
            f"{args.frontier_offset}; {describe(before_frontier[0])}",
            file=sys.stderr,
            flush=True,
        )
        per_worker = math.ceil(len(before_snapshot) / args.workers)
        command = [
            sys.executable,
            str(DRIVER),
            str(args.p),
            "--workers",
            str(args.workers),
            "--seed",
            str(args.seed + round_index * args.workers),
            "--factor-seconds",
            str(args.factor_seconds),
            "--factor-flags",
            str(args.factor_flags),
            "--factor-rounds",
            "2",
            "--prefactor-seconds",
            str(args.prefactor_seconds),
            "--prefactor-flags",
            "8",
            "--pm1-bounds",
            str(args.pm1_bound),
            "--pp1-bounds",
            str(args.pp1_bound),
            "--ecm-bounds",
            str(args.ecm_bound),
            "--ecm-curves",
            str(args.ecm_curves),
            "--ecm-rounds",
            "1",
            "--curve-families",
            "5",
            "--candidate-bits",
            "1",
            "--candidate-dir",
            str(candidate_dir),
            "--resume-candidates",
            str(candidate_dir),
            "--resume-offset",
            str(args.frontier_offset),
            "--resume-top",
            str(len(before_snapshot)),
            "--resume-per-worker",
            str(per_worker),
            "--resume-only",
            "--branch-curves",
            str(args.branch_curves),
            "--manifest",
            str(manifest),
            "--output",
            str(output),
        ]
        result = subprocess.run(command, cwd=ROOT, check=False)
        if result.returncode:
            return result.returncode
        if output.exists():
            certificate = output.read_text(encoding="utf-8").strip()
            if not verified_for(certificate, args.p):
                raise SystemExit(f"driver wrote an invalid certificate to {output}")
            print(certificate)
            return 0

        after = load_resume_candidates(candidate_dir, args.p)
        after_frontier = after[
            args.frontier_offset : args.frontier_offset + args.frontier_size
        ]
        if not after_frontier:
            print("selected frontier slice exhausted", file=sys.stderr, flush=True)
            break
        after_snapshot = snapshot(after_frontier, args.frontier_size)
        if after_snapshot == before_snapshot:
            stale += 1
            print(
                f"frontier made no exact progress ({stale}/{args.stale_rounds} "
                "stale pass(es))",
                file=sys.stderr,
                flush=True,
            )
            if stale >= args.stale_rounds:
                break
        else:
            stale = 0
            print(
                f"frontier improved; new leader: {describe(after_frontier[0])}",
                file=sys.stderr,
            )

    print("adaptive frontier exhausted without a certificate", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
