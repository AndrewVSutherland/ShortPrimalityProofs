# ShortPrimalityProofs
For the purpose of thie repository, a **short ECPP** is a sequence of integers $(p_0,A_0,x_0,m_0p_1,A_1,x_1,\ldots,m_{k-1}p_k,A_k,x_k,m_kp_{k+1})$ in which
- we have $k \ge 0$, $p_0 \ge 5$, and put $n:=\lceil \log_2 p_i\rceil$,
- the $p_i$ are odd integers that satisfy $n^2 < p_{i+1} < \sqrt{p_i}$ for $0\le i < k$ and $p_{k+1}=1$,
- the $m_i$ are $n^2$-smooth integers, where $n=\lceil\log_2 p_0\rceil$, satisfiying $L_i < m_ip_{i+i} < r_iL_i$, where $L_i = q_i+1+\lfloor 2\sqrt{q_i}\rfloor$ with $q_i=\lfloor\sqrt{p_i}\rfloor$ and $r_i$ is the least prime divisor of $m_i$,
- each $A_i$ is a nonnegative integer less than $p_i$ with $a_i\ne \pm 2\bmod p_i$,
- each $x_i$ is a nonnegative integer less than $p_i$,

such that for each $0\le i \le k$ there exist integers $B_i,y_i\in [0,p_i-1]$ for which $(x_i,y_i)$ is a point of order $m_ip_{i+1}$ on the [Montgomery curve](https://en.wikipedia.org/wiki/Montgomery_curve) $B_iy^2 = x^3 + A_ix^2 +x$.
