#!/usr/bin/env python3
"""
Verification of a short ECPP (github.com/AndrewVSutherland/ShortPrimalityProofs)
in quasi-quadratic time.  Montgomery-ladder core follows voneshot.py (Opus 4.8 /
A.V. Sutherland); chain logic per the ShortPrimalityProofs definition, as
revised in August 2026 (the radical-capped format).

A short ECPP for p_0 is a flat integer sequence

    (p_0, A_0, x_0, o_0, A_1, x_1, o_1, ..., A_k, x_k, o_k)

with o_i = m_i * p_{i+1} and p_{k+1} = 1, in which, writing n = ceil(log2 p_0)
(fixed at the top level for the whole chain) and B = floor(n^2 / log2 n):

  - each p_i is odd (p_0 given; later moduli are recovered from the previous
    level as p_{i+1} = the B-rough part of o_i), with p_{i+1} = 1 or
    n^2 < p_{i+1} and p_{i+1}^2 < p_i (the floor is explicit: roughness only
    guarantees p_{i+1} > B);
  - m_i, the B-smooth part of o_i, satisfies the radical cap
        floor(log2 rad(m_i)) < n / log2 n,
    where rad(m) is the product of the distinct primes of m, and, with r_i its
    least prime divisor,
        L_i < o_i < r_i * L_i,
    where q_i = isqrt(p_i) and L_i = q_i + 1 + isqrt(4*q_i) is an upper bound
    on the order of a point on an elliptic curve over F_q for any q <= sqrt(p_i);
  - 0 <= A_i < p_i with gcd(A_i^2 - 4, p_i) = 1, and 0 <= x_i < p_i;
  - x_i is the x-coordinate of a point of order exactly o_i on the Montgomery
    curve B y^2 = x^3 + A_i x^2 + x over Z/p_i for some B -- equivalently on
    E_{A_i} or its quadratic twist; the x-only ladder never needs B or y.

Every certificate of this revised format is also valid under the original
(n^2-smooth, uncapped) definition: the new conditions only restrict.

Why this proves p_0 prime: by induction from the top.  At level i the verifier
establishes ord(P_i) = o_i modulo every prime divisor l of p_i, with every prime
factor of o_i certified (primes <= B by trial division, p_{i+1} by the next
level, 1 trivially).  If some prime l <= sqrt(p_i) divided p_i, then o_i =
ord(P_i mod l) <= #E(F_l) <= l + 1 + floor(2*sqrt(l)) <= L_i < o_i, a
contradiction; hence p_i is prime.  The minimality window o_i < r_i * L_i keeps
the certificate size O(n) (sum of the level sizes is geometric).

Soundness for composite input is inherited from the voneshot.py machinery: the
x-only ladder is valid on the Kummer line of E and its twist so no on-curve test
is needed; [o]P = O must be reached as a genuine (X:0) with gcd(X, p) = 1
(rejecting the degenerate (0:0) collapse that arises from multiplying past the
order modulo a factor of p); each (o/q)P != O leaf requires gcd(Z, p) = 1, i.e.
nonvanishing modulo every prime divisor of p.

Cost (FFT integer multiplication assumed): this is what the August 2026 caps
buy.  The primorial P_B = prod of the primes <= B has Theta(n^2/log n) bits, so
it is built per certificate in O(n^2 log n) bit operations -- no precomputed
tables.  Per level, gcd(P_B mod o_i, o_i) equals rad(m_i) exactly (P_B is
squarefree and complete to B), an integer of fewer than n/log2 n bits whose
prime factors are at most B; the radical cap is read off its bit length, and
factoring it yields the distinct primes of m_i, whose valuations recover m_i
and p_{i+1}.  The cap bounds the order-exactness tree by fewer than n ladder
bits per level, so the elliptic work is O(sum_i b_i M(b_i)) = O(n^2 log n)
with level sizes decaying geometrically.  Total: O(n^2 log n) bit operations,
worst case over certificates, and O(n^2/log n) bits of memory.  (This
implementation factors each radical by direct trial division over the sieved
primes -- asymptotically the coarsest step, but the simplest, and dominant
only far beyond feasible sizes; a remainder tree over the primorial's product
tree, or Pollard--Strassen with search radius B, keeps even that step within
O(n^2 log n).)

usage: python3 vsmallECPP.py p0 A0 x0 o0 [A1 x1 o1 ...]
       python3 vsmallECPP.py --test
Prints True and exits 0 iff the sequence is a valid short ECPP (p_0 is prime).
"""

from math import gcd, isqrt, log2


# --------------------------------------------------------------------------
# Montgomery x-only (X:Z) arithmetic, valid on the Kummer line of E_A and of
# its quadratic twist (verbatim from voneshot.py).
# --------------------------------------------------------------------------
def xdbl(X, Z, A, p):
    XX = X * X % p
    ZZ = Z * Z % p
    XZ = X * Z % p
    X2 = (XX - ZZ) * (XX - ZZ) % p
    Z2 = 4 * XZ % p * ((XX + A * XZ + ZZ) % p) % p
    return X2, Z2


def xadd(X1, Z1, X2, Z2, Xd, Zd, p):
    a = (X1 - Z1) * (X2 + Z2) % p
    b = (X1 + Z1) * (X2 - Z2) % p
    s = (a + b) % p
    d = (a - b) % p
    X3 = Zd * (s * s % p) % p
    Z3 = Xd * (d * d % p) % p
    return X3, Z3


def ladder(k, XP, ZP, A, p):
    if k == 0:
        return (1, 0)
    XP %= p
    ZP %= p
    if k == 1:
        return (XP, ZP)
    Xd, Zd = XP, ZP
    X0, Z0 = XP, ZP
    X1, Z1 = xdbl(XP, ZP, A, p)
    for bit in bin(k)[3:]:
        if bit == '0':
            X1, Z1 = xadd(X0, Z0, X1, Z1, Xd, Zd, p)
            X0, Z0 = xdbl(X0, Z0, A, p)
        else:
            X0, Z0 = xadd(X0, Z0, X1, Z1, Xd, Zd, p)
            X1, Z1 = xdbl(X1, Z1, A, p)
    return X0, Z0


def sieve_primes(limit):
    if limit < 2:
        return []
    is_p = bytearray([1]) * (limit + 1)
    is_p[0] = is_p[1] = 0
    for i in range(2, isqrt(limit) + 1):
        if is_p[i]:
            is_p[i * i::i] = bytearray(len(is_p[i * i::i]))
    return [i for i in range(2, limit + 1) if is_p[i]]


def balanced_product(xs):
    """prod(xs) by power-of-two pairing: O(M(S) log k) bit operations for k
    factors of total size S, versus Theta(S^2/k) for sequential accumulation."""
    if not xs:
        return 1
    if len(xs) <= 8:
        r = xs[0]
        for x in xs[1:]:
            r *= x
        return r
    xs = list(xs)
    while len(xs) > 1:
        nxt = [xs[i] * xs[i + 1] for i in range(0, len(xs) - 1, 2)]
        if len(xs) & 1:
            nxt.append(xs[-1])
        xs = nxt
    return xs[0]


# --------------------------------------------------------------------------
# Order-exactness tree (from voneshot.py): with Q = (o/R)P for R the product of
# the distinct primes of o, each leaf holds (o/q)P, whose Z must be a unit.
# The radical cap keeps this tree logarithmically small.
# --------------------------------------------------------------------------
def check_orders(XQ, ZQ, primes, A, p):
    t = len(primes)
    if t == 0:
        return True
    if t == 1:
        return gcd(ZQ % p, p) == 1
    mid = t // 2
    Lh, Rh = primes[:mid], primes[mid:]
    hL = balanced_product(Lh)
    hR = balanced_product(Rh)
    XL, ZL = ladder(hL, XQ, ZQ, A, p)
    XR, ZR = ladder(hR, XQ, ZQ, A, p)
    return check_orders(XL, ZL, Rh, A, p) and check_orders(XR, ZR, Lh, A, p)


# --------------------------------------------------------------------------
# The short-ECPP verifier.
# --------------------------------------------------------------------------
def verify(seq):
    """True iff seq = (p0, A0, x0, o0, ..., Ak, xk, ok) is a valid short ECPP."""
    if len(seq) < 4 or (len(seq) - 1) % 3 != 0:
        return False
    if any(not isinstance(v, int) or v < 0 for v in seq):
        return False
    p = seq[0]
    if p < 5 or p % 2 == 0:
        return False
    n = p.bit_length()            # = ceil(log2 p0): p0 is odd, never a power of 2
    lg = log2(n)
    B = int(n * n / lg)           # floor(n^2 / log2 n): the smoothness bound
    radlim = n / lg               # require floor(log2 rad(m)) < radlim

    # collect the level orders and pre-screen their sizes so that the work below
    # runs on a certificate-independent budget.  Both bounds are implied by
    # validity, so rejecting on them is sound: a valid chain has at most
    # 1 + log2(n) levels (the moduli at least halve in bit length and stay above
    # n^2 until the last), and every level order is below the Hasse bound of its
    # modulus, hence below p0^2.
    os = [seq[i + 2] for i in range(1, len(seq), 3)]
    if len(os) > n:
        return False
    if any(o < 2 or o.bit_length() > 2 * n + 2 for o in os):
        return False

    primes = sieve_primes(B)
    P = balanced_product(primes)                  # the primorial of B

    for i in range(1, len(seq), 3):
        A, x, o = seq[i], seq[i + 1], seq[i + 2]
        if p < 3 or p % 2 == 0:   # a mid-chain modulus collapsed to 1 (or worse)
            return False
        if not (0 <= A < p) or not (0 <= x < p):
            return False
        if gcd((A * A - 4) % p, p) != 1:          # nonsingular mod every divisor of p
            return False
        if o < 2:
            return False

        # recover rad(m) = gcd(P mod o, o) -- exact, since P is squarefree and
        # complete to B -- then the primes of m, then m and p_next by valuation
        g = gcd(P % o, o)
        if g <= 1:                                # m = 1: r undefined, reject
            return False
        if g.bit_length() - 1 >= radlim:          # the radical cap
            return False
        small = []                                # ascending prime factors of g
        gg = g
        for q in primes:
            if q * q > gg:
                break
            if gg % q == 0:
                small.append(q)
                gg //= q
        if gg > 1:
            small.append(gg)                      # prime (trial division passed sqrt)
        m = 1                                     # m = prod over q | g of q^{v_q(o)}
        for q in small:
            oo = o
            while oo % q == 0:
                oo //= q
                m *= q
        p_next = o // m
        r = small[0]                              # least prime divisor of m (and of o)

        # size window: L < o < r*L, L an upper bound on point orders over any F_q,
        # q <= sqrt(p)
        q_ = isqrt(p)
        L = q_ + 1 + isqrt(4 * q_)
        if not (L < o < r * L):
            return False

        # descent: p_next = 1, or n^2 < p_next (explicit floor) and p_next^2 < p
        if p_next != 1 and (p_next <= n * n or p_next * p_next >= p):
            return False

        # [o]P = O reached as a genuine (X:0) with X a unit mod p
        Xo, Zo = ladder(o, x, 1, A, p)
        if Zo % p != 0 or gcd(Xo % p, p) != 1:
            return False

        # (o/q)P != O (mod every prime divisor of p) for each prime q | o;
        # p_next participates as a factor, its primality certified by level i+1
        divisors = small + ([p_next] if p_next > 1 else [])
        R = balanced_product(divisors)
        XQ, ZQ = ladder(o // R, x, 1, A, p)
        if not check_orders(XQ, ZQ, divisors, A, p):
            return False

        p = p_next
    return p == 1                                 # the chain must terminate exactly


# --------------------------------------------------------------------------
# Self-tests: valid vectors produced by short.gp and by exhaustive search over
# small p, plus tamper cases exercising each rejection path (including the new
# radical cap and smoothness bound).  The CRT split attack is inherited from
# the original test suite and must still be rejected.
# --------------------------------------------------------------------------
_VALID = [
    # single-level certificates found by exhaustive search (see short8all.txt)
    "5 4 3 8",
    "7 0 2 8",
    "11 8 8 8",
    "13 7 2 8",
    "251 180 183 24",
    # the (migrated) 10^10+19 entry of certs.csv
    "10000000019 9322349340 1921958667 116108 15213 3538 243",
]

_INVALID = [
    "251 0 10 63",                 # valid in the ORIGINAL format; the radical cap rejects
    "11 0 3 12",                   # ditto (rad 6: floor(log2 6) = 2 is not < 4/2)
    "251 0 10 126",                # m doubled: breaks the minimality window
    "251 2 10 32",                 # singular curve (A = 2)
    "221 5 2 34",                  # composite p0 (221 = 13*17)
    "3 0 0 6",                     # p = 3 admits no short ECPP at all
    "10000000019 9322349340 1921958667 116108 15213 3538 243 0 1 8",
                                   # trailing level after the chain terminated
    # CRT split attack: p0 = 2098153*2102167; (x0,1) has order exactly 8*525029
    # in E(Z/p0) but order 8 mod one factor and 525029 mod the other
    "4410667997551 1365834658413 107710304518 4200232 199129 175565 880",
]


if __name__ == "__main__":
    import sys
    if len(sys.argv) == 2 and sys.argv[1] == "--test":
        good = True
        for s in _VALID:
            v = verify([int(t) for t in s.split()])
            print(f"valid vector accepted: {v}")
            good = good and v
        for s in _INVALID:
            v = verify([int(t) for t in s.split()])
            print(f"invalid vector rejected: {not v}")
            good = good and not v
        print("ALL TESTS PASSED" if good else "TESTS FAILED")
        sys.exit(0 if good else 1)
    try:
        seq = [int(t) for t in sys.argv[1:]]
    except ValueError:
        print("usage: vsmallECPP.py p0 A0 x0 o0 [A1 x1 o1 ...]")
        sys.exit(2)
    r = verify(seq)
    print(r)
    sys.exit(0 if r else 1)
