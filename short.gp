/* short.gp -- compute a short ECPP certificate for a probable prime p >= 5.
 * Written by Fable 5.  A toy implementation, in the style of oneshot.gp; the goal is
 * clarity, not speed.
 *
 * Format: github.com/AndrewVSutherland/ShortPrimalityProofs.  The certificate is the flat
 * sequence  (p_0, A_0, x_0, o_0, A_1, x_1, o_1, ..., A_k, x_k, o_k)  with o_i = m_i*p_{i+1},
 * p_{k+1} = 1 and n = ceil(log_2 p_0) FIXED for the whole chain: on E_{A_i} : y^2 = x^3 +
 * A_i x^2 + x over F_{p_i} the point with x-coordinate x_i has order exactly o_i, where m_i
 * is B-smooth for B = ceil(n^2/log_2 n) with log_2 rad(m_i) <= B (rad(m) is
 * the product of the distinct primes of m), p_{i+1} is B-rough and either a prime in
 * (B^2, sqrt(p_i)) proven prime by the next level, or -- at the last level only -- 1 or
 * below B^2 (a B-rough integer below B^2 is 1 or prime, so it certifies itself by size),
 * and L_i < o_i < r_i L_i with L_i = (p_i^{1/4}+1)^2 and r_i the least prime divisor of m_i.
 *
 * Method.  For each level we search random curves E_A/F_{p_i}, compute N = #E_A by SEA
 * (ellcard), and look for a divisor o | N inside the narrow window (L_i, r_i L_i) -- note
 * o is only about sqrt(p_i), so most of N is discarded.  A usable o must be of the shape
 * (B-smooth with a small radical) * (one prime), so we trial-divide N up to B to get its
 * smooth part s, and then look for the single large prime among the factors of the rough
 * part N/s.  That last step is the expensive one: the rough part is factored under an
 * `alarm` time budget (SC_tlim seconds) and the curve is abandoned if the budget runs out
 * -- an early abort that keeps the search on curves whose order factors easily.  A level
 * whose large prime is below B^2 (or absent) ends the chain.
 *
 * Usage:
 *     echo 'printshort(nextprime(10^30))' | gp -q short.gp
 *     echo 'SC_tlim = 60; printshort(nextprime(10^100))' | gp -q short.gp
 */

default(parisizemax, 2^32);                              \\ let the stack grow for large ellcard calls

SC_curves = 0;                                           \\ curves tried (all levels)
SC_tlim = 20;                                            \\ seconds allowed per rough-part factorization
                                                         \\ (SC_tlim = 0 disables step (c): then the only
                                                         \\  large prime considered is a prime rough part)
SC_maxcurves = 0;                                        \\ 0 = unlimited; else give up after this many curves

/* B-smooth part s of N and the rough cofactor r = N/s, by trial division over primes <= B */
smoothpart(N, B) = {
  my(s = 1, r = N);
  forprime(q = 2, B, while(r % q == 0, r /= q; s *= q));
  [s, r];
};

scbound(p) = sqrtint(p) + 1 + sqrtint(4 * sqrtint(p));   \\ integer form of L = (p^{1/4}+1)^2

/* radical condition: log_2 rad(m) <= B */
radcap(m, B) = #binary(vecprod(factor(m)[,1])) <= B;

/* A point of order exactly o on E (N = #E, o | N, fo = the primes of o), or 0 if none found. */
scpoint(E, N, o, fo) = {
  my(Q, ok);
  for(t = 1, 64,
    Q = ellmul(E, random(E), N/o);                       \\ order(Q) divides o
    if(#Q == 1, next);                                   \\ Q = O, resample
    ok = 1;
    for(i = 1, #fo, if(#ellmul(E, Q, o/fo[i]) == 1, ok = 0; break));
    if(ok, return(Q))
  );
  0;
};

/* One level: search random curves over F_p for [A, x0, o, p_next].  p_next = 1 ends the chain
 * (including the self-certifying case: a prime q < B^2 stays inside o and needs no recursion). */
sclevel(p, B, B2) = {
  my(L = scbound(p), rt = sqrtint(p), A, E, N, sr, s, R, dv, F, m, o, q, Q);
  while(1,
    if(SC_maxcurves && SC_curves >= SC_maxcurves, return(0));
    SC_curves++;
    A = random(p); E = ellinit([0, A, 0, 1, 0], p);
    if(#E == 0, next);                                   \\ singular (A = +-2 mod p)
    N = ellcard(E);                                      \\ SEA
    sr = smoothpart(N, B); s = sr[1]; R = sr[2];
    dv = divisors(s);

    \\ (a) terminal level: o = m is B-smooth all by itself, radical-capped, in the window
    for(i = 1, #dv, m = dv[i];
      if(m > L && m < factor(m)[1,1] * L && radcap(m, B),
        Q = scpoint(E, N, m, factor(m)[,1]);
        if(Q != 0, return([A, lift(Q[1]), m, 1]))));

    \\ (a') terminal level with a self-certifying prime: the rough part is below B^2, hence
    \\      automatically prime (no primality test needed!), and stays inside o
    if(R > 1 && R < B2,
      for(i = 1, #dv, m = dv[i]; o = m * R;
        if(m > 1 && o > L && o < factor(m)[1,1] * L && radcap(m, B),
          Q = scpoint(E, N, o, concat(factor(m)[,1], [R]~));
          if(Q != 0, return([A, lift(Q[1]), o, 1])))));

    \\ (b) cheap descent: the rough part is itself prime, so q = R needs no factoring.  Then the
    \\     discarded cofactor N/o = s/m divides s, i.e. the whole ~sqrt(p) worth of N that we
    \\     throw away must be smooth, which is very rare for large p -- hence (c).  Costs one
    \\     primality test.
    if(s > L && R > B2 && R < rt && ispseudoprime(R),
      q = R;
      for(i = 1, #dv, m = dv[i]; o = m * q;
        if(m > 1 && o > L && o < factor(m)[1,1] * L && radcap(m, B),
          Q = scpoint(E, N, o, concat(factor(m)[,1], [q]~));
          if(Q != 0, return([A, lift(Q[1]), o, q])))));

    \\ (c) general descent: o = m*q with q any prime factor of the rough part, m | s, m > 1.
    \\     Here the discarded cofactor N/o is unconstrained, which is far likelier -- but q must
    \\     be dug out of R, so we factor R under a time budget and abandon slow curves.
    if(SC_tlim > 0 && R > 1 && s * rt > L,
      F = iferr(alarm(SC_tlim, factor(R)[,1]), e, 0);    \\ early abort on slow factorizations
      if(type(F) == "t_COL",
        for(j = 1, #F, q = F[j];
          if(q < rt,
            for(i = 1, #dv, m = dv[i]; o = m * q;
              if(m > 1 && o > L && o < factor(m)[1,1] * L && radcap(m, B),
                Q = scpoint(E, N, o, concat(factor(m)[,1], [q]~));
                if(Q != 0, return([A, lift(Q[1]), o, if(q < B2, 1, q)]))))))))
  );
};

/* The full chain: returns the flat sequence (p, A_0, x_0, o_0, ..., A_k, x_k, o_k). */
shortcert(p) = {
  if(!ispseudoprime(p), error("short: p is composite"));
  if(p < 5, error("short: need p >= 5"));
  my(n = #binary(p), lg = log(n)/log(2), B = ceil(n^2/lg), B2 = B^2,
     seq = [p], cur = p, lev);
  SC_curves = 0;
  while(cur > 1,
    lev = sclevel(cur, B, B2);
    if(lev == 0, return(0));                             \\ gave up (SC_maxcurves reached)
    seq = concat(seq, [lev[1], lev[2], lev[3]]);
    cur = lev[4]
  );
  seq;
};

printshort(p) = {                                        \\ "p A0 x0 o0 A1 x1 o1 ..."
  my(c = shortcert(p));
  for(i = 1, #c, printf("%d%s", c[i], if(i < #c, " ", "\n")));
};
