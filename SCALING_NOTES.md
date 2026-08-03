# Scaling notes

This file records reproducible performance observations from the first search beyond
the 200-digit challenge entry.  Search manifests and candidate files are intentionally
not committed; the commands below are sufficient to reproduce their structure.

## Test system

- Apple M4, 10 physical cores, 16 GB RAM
- PARI/GP 2.17.4 with `seadata`
- GMP-ECM 7.0.6
- macOS, 3 August 2026

The target used below is the least prime greater than `10^210`, the first line of
`targets.txt`.

## Smooth-part extraction

`smoothpart` now caches the square-free product of every prime through the certificate
bound and repeatedly takes gcds with the order.  On ten representative 701-bit orders,
the old prime-by-prime loop took 53 ms and the cached-gcd path took 4 ms, a 13.2x
screening speedup.  Building the 701-bit target's cache took 128 ms; its product was
701,527 bits (about 86 KiB).  The old and new paths returned identical smooth and rough
parts.

The cache is shared by every root checkpoint attempted in one GP process.  It is
therefore especially useful with `--resume-per-worker`, which avoids rebuilding both
the cache and PARI's private factor table between candidate orders.

## Two-pass CM search

Before invoking Cornacchia, the scanner now tests the necessary quadratic-residue
condition.  A representation `4p=t^2+Dv^2` implies that `-D` is a square modulo
`p` (and the corresponding statement holds for the even-discriminant form), so a
Kronecker-symbol rejection is exact.  On the first 50,000 odd/even discriminant
tests at the 210-digit target, this reduced Cornacchia calls to 15,229 and wall
time from 22.65 seconds to 10.03 seconds, a 2.26x speedup, while preserving all
239 represented orders.

The production enumeration used no SEA calls and no factoring:

```sh
read -r p < targets.txt
python3 parallel_short.py "$p" -j 10 --seed 2113000 \
  --cm-bound 100000000 --cm-screen-only \
  --candidate-bits 40 --candidate-dir candidates/210 \
  --manifest search-runs-210.jsonl
```

It completed in 1,264.6 seconds of wall time and 9,634.9 aggregate child CPU seconds.
The run saved 1,095 distinct live orders before later factoring and tombstones changed
the frontier, reached an 86-bit smooth part, and made zero SEA calls.  Equal-width
discriminant partitioning kept the finite scan bounded and avoided the severe tail seen
when factor work and representation density were both assigned by square-root ranges.

The expensive second stage ranked the saved orders and assigned twenty to each worker:

```sh
python3 parallel_short.py "$p" -j 10 --seed 2115000 \
  --factor-seconds 20,20,20,20,20,20,5,5,5,5 \
  --factor-flags 0,0,0,0,0,0,9,9,9,9 \
  --factor-rounds 3,3,3,3,3,3,2,2,2,2 \
  --prefactor-seconds 2 --pm1-bounds 1000000 --pp1-bounds 1000000 \
  --ecm-bounds 0,0,0,0,0,0,0,0,100000,100000 \
  --ecm-curves 0,0,0,0,0,0,0,0,8,8 \
  --resume-candidates candidates/210 --resume-top 200 \
  --resume-per-worker 20 --resume-only --candidate-bits 40 \
  --candidate-dir candidates/210 --curve-families 5 \
  --branch-curves 96 --manifest search-runs-210.jsonl -o cert210.txt
```

This pass exhausted its 200 assigned checkpoints in 564.6 seconds and 2,387.0 child
CPU seconds without a certificate.  It nevertheless demonstrated why checkpoints must
retain partial factor progress: one cofactor fell from 502 to 397 bits, and another
from 512 to 444 bits.  Re-ranking these reduced cofactors made them the first two items
in the next portfolio.  In the fast worker profiles, P-1 and P+1 found factors in many
orders; the two GMP-ECM profiles reported 17 factor-bearing batches in 34 attempts.

The practical workflow is consequently an adaptive frontier, not a fixed seed list:

1. enumerate many exact orders with no factor budget;
2. apply cheap P-1, P+1, and short partial-factor passes to the ranked frontier;
3. persist every smaller residual and re-rank;
4. spend independent ECM curves and longer complete-factor limits only on the new
   frontier; and
5. report both wall time and aggregate core-seconds, plus actual `resume_attempts`.

## Low-smoothness portfolio lane

The 40-bit CM screen is a ranking optimization, not a completeness threshold.  Applying
the verifier's exact smooth/rough split to the published 110--200 digit root levels gives
only 3--29 smooth bits.  A portfolio that factors only the 40-bit frontier therefore
omits the regime that produced every existing large certificate.

A minority lane should instead set `--cm-smooth-bits 0`, use the proven complete
`--factor-seconds 20 --factor-flags 0` profile, and checkpoint from one smooth bit.  Exact
CM orders remove the SEA cost from this lane, so its throughput is limited primarily by
the intentionally bounded factor calls:

```sh
python3 parallel_short.py "$p" -j 1 --seed 2123100 \
  --factor-seconds 20 --factor-flags 0 --factor-rounds 3 \
  --cm-bound 5000000 --cm-smooth-bits 0 --cm-only \
  --candidate-bits 1 --candidate-dir candidates/210 \
  --branch-curves 64 --manifest search-runs-210.jsonl -o cert210.txt
```

A factor-free four-worker scan through `|D|=50,000,000` with
`--cm-smooth-bits 1 --candidate-bits 1 --cm-screen-only` completed in 1,382.4 seconds
wall and 2,567.6 aggregate child CPU seconds.  It wrote 8,850 append-only checkpoints
representing 8,026 distinct live orders before later factoring.  This measured pipeline
is more productive than factoring inline: screen the finite discriminant interval at
full speed, then give the globally ranked low-smoothness orders their bounded factor
passes.

A five-worker continuation over `50,000,001 <= |D| <= 200,000,000` completed in
3,083.5 seconds wall and 7,699.7 aggregate child CPU seconds.  It wrote 9,103
checkpoints representing 8,233 distinct viable orders before downstream factoring.
The wider interval therefore remains practical as a first-stage production scan, while
the explicit child-CPU total makes clear that its wall-clock speed comes from parallel
discriminant coverage rather than reduced aggregate work.

Keep these full-factor lanes at modest discriminants.  `scmontgomerylevel` currently
calls PARI's integer `polclass`; a rare factor win at a very large discriminant can make
curve reconstruction more expensive than the order search.  Computing one CM root
directly modulo `p` is the appropriate next implementation step before extending this
lane to substantially larger discriminants.

The frontier also applies an exact impossibility bound.  Every unresolved prime factor
is larger than `n^2`; for a composite residual `R`, a future child prime is therefore at
most `R/(n^2+1)`.  Once that upper bound times the known smooth part cannot clear the
certificate window, GP writes a tombstone and the Python loader suppresses every older
checkpoint for the same order.  A focused ECM pass triggered this rule by reducing a
397-bit residual to 310 bits: another split would leave at most 291 bits, below that
order's 299-bit minimum.

## Adaptive frontier iteration

A ranked queue should be refreshed after every cheap saved-factor round instead of being
treated as a static batch.  One 20-checkpoint pass at the 210-digit target used 5-second
partial PARI factoring, 2-second prefactoring, one P-1 and P+1 pass at `B1=5,000,000`,
and eight GMP-ECM curves at `B1=100,000`.  It completed in 335.9 seconds wall and 167.3
child CPU seconds.  Two ECM batches found factors: the `D=-7331` residual fell from 505
to 425 bits, and the `D=-15600715` residual fell from 519 to 452 bits.

The first reduction changed the exact maximum complementary cofactor from 177 bits to
97 bits, promoting the order ahead of nearly the entire previous queue.  Continuing the
old static list would not act on this much stronger lead.  The operational loop is thus:

1. run a bounded cheap queue;
2. reload every append-only checkpoint;
3. retain the numerically smallest residual even when its bit length is unchanged;
4. rank by the exact bound `floor(residual/minimum_q)`, not its coarse bit length; and
5. launch deeper independent ECM only on the refreshed frontier.

This separates high-throughput factor discovery from expensive tail work and makes the
portfolio responsive to partial progress without weakening any certificate condition.
The checked-in `adaptive_frontier.py` controller implements this loop and stops after a
full pass leaves the exact ranked snapshot unchanged.  This prevents a long unattended
run from repeatedly applying the same cheap layer after its useful frontier reductions
have been harvested.

A later low-smoothness order at `D=-15226867` reached a 419-bit residual with an exact
maximum complementary factor of 83 bits.  An independent batch of 100 GMP-ECM curves at
`B1=1,000,000` completed in 325.6 seconds wall and 146.8 child CPU seconds without a
factor.  The exact bound therefore ranks the opportunity correctly but is not evidence
that such a factor exists; after one appropriately deep negative batch, fresh orders
again have higher expected value than repeatedly grinding the same residual.

A later 210-digit frontier order reduced from 425 to 332 residual bits, leaving an exact
maximum complementary cofactor of 17,675,219.  A complete prime-divisor scan from
`n^2+1` through that bound took 110 ms and found no divisor, proving that the order could
not contain an admissible child.  The same worker had otherwise begun 300 ECM curves at
`B1=1,000,000`.  Both order paths now perform this exact scan automatically whenever the
bound is at most 20,000,000, converting a potentially long randomized tail into a fast
deterministic tombstone.

For comparison, a one-off audit raised the limit for the apparent `D=-36644`
leader, whose exact complementary bound was 6,593,364,716.  The exhaustive scan
finished in 105.1 seconds wall and 47.4 child CPU seconds, found no divisor, and
proved that order unusable.  The default remains 20,000,000 to keep per-order
latency negligible; a multi-billion scan is worthwhile only as an explicit
audit of an unusually strong frontier item.

A later `D=-14995663` order illustrates the same deterministic finish after a
successful deep factor pass.  GMP-ECM reduced its residual from 455 to 361 bits,
leaving an exact complementary bound of 33,761,519,371 (35 bits).  A one-off scan
through that bound completed in 274.6 seconds wall and 236.4 child CPU seconds,
found no divisor, and tombstoned the order.  The result distinguishes a genuinely
impossible root from another randomized-factor failure and prevents all future
portfolios from replaying it.

## Torsion-conditioned curves

Curve family 5 ports the optimized `X_1(27)` construction from OneShotSEA.  The retained
side has full rational 2-torsion and a point of order four, so its group order is
divisible by 216; the twist remains available for the same point count.  Generation at
701 bits takes a fraction of a second, whereas a PARI SEA call on an unconditioned curve
at this size takes roughly three to four minutes on one core.  A complete two-worker
20-digit smoke search found and verified a certificate with family 5 in 0.3 seconds.

## Smaller CM class invariants

CM reconstruction now selects PARI's Weber class invariant whenever the discriminant
permits it.  For `D = 1 (mod 8)`, invariant 1 supplies Weber `f` and maps a modular root
back through `j=(f^24-16)^3/f^24`; if 3 ramifies, invariant 3 supplies `f^3` and uses the
same formula with its eighth power.  Other discriminants retain `gamma_2` or classical
`j`.  This does not change the exact order or certificate model, but can reduce the
height of the class polynomial substantially before a screened order is reconstructed.

As a small correctness and size check, `D=-2263` produced degree-22 polynomials occupying
752 bytes with Weber `f` versus 1,232 bytes with `gamma_2`.  Modular roots for invariant
1 (`D=-23`) and the 3-ramified invariant 3 fallback (`D=-15`) were converted to `j` and
used to recover Montgomery curves having the prescribed CM orders.  The advantage grows
in coefficient size with `|D|`; the small examples remain overhead-dominated in elapsed
time, so this is a reconstruction-memory improvement rather than a claim that every tiny
class polynomial is faster.

## OneShotSEA boundary at 701 bits

OneShotSEA commit `c0bbf81` was built with its native C++20/GMP path.  Its generic
Weber-compatible curve at seed 2115000, index 0 was counted with one SEA thread.

With the checked-in 77 levels through 401, the run used 210.66 seconds wall and 179.95
seconds user CPU, then correctly stopped at the level limit with about 7.95e22 compatible
traces.  An arbitrary Montgomery `X_1(27)` curve was not eligible for this path because it
had no rational Weber lift.

The pinned archive was then materialized through level 601 (108 authenticated tables,
53 MB).  The second run stopped after level 461 in 557.07 seconds wall and 317.26 seconds
user CPU.  It soundly enumerated the complete remaining set of 67,306 traces, but did not
produce a unique point count.  Thus the current generic Weber path is not a practical
replacement for PARI SEA at 701 bits on this machine.

The useful longer-term route is narrower:

- generate `X_1(27)` curves together with their authenticated Weber coordinate, rather
  than converting an unrelated Montgomery curve after generation;
- feed the exact 216/432 torsion trace prior into SEA;
- apply exact smooth-part rejection to the complete partial trace set before paying for
  more levels; and
- replace finite high-level tables with direct finite-field specialization as the input
  grows.

The CM two-pass search is currently the better practical engine for the 210--300 digit
targets.  The Weber route remains the asymptotically motivated next-round design, but the
measurements above identify the missing integration points instead of assuming a
crossover that has not occurred.
