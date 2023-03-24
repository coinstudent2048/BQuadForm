# Miscellaneous utilities for bquadform.py
#
# Use this code only for prototyping

# integer square root (floor)
# reference: Algorithm 1.7.1 (p.38) of Cohen -
# "A Course in Computational Algebraic Number theory" (GTM 138)
def isqrt(n):
    if not isinstance(n, int):
        raise TypeError("input is not integer")
    if n < 0:
        raise ValueError("input is negative")
    if n == 0:
        return 0
    # [Initialize]
    x = 1 << (n.bit_length() + 2 >> 1)
    while True:
        # [Newtonian step]
        y = (x + n // x) >> 1
        # [Finished?]
        if y >= x:
            return x
        x = y

# # integer square root (ceiling)
# def isqrt_ceil(n):
#     if n == 0:
#         return 0
#     else:
#         return 1 + isqrt(n - 1)

# Euclidean division: always ensures that
# 0 <= r < |b| regardless of sign of divisor
def divmod_euclid(a, b):
    q, r = divmod(a, b)   # divmod uses floor division
    if r < 0:
        q += 1
        r -= b
    return (q, r)

# extended Euclidean algorithm (assumes a >= 0 & b >= 0)
# reference: Algorithm 1.3.6 (p.16) of Cohen -
# "A Course in Computational Algebraic Number theory" (GTM 138)
def ext_euclid(a, b, u=1, v1=0):
    # [Initialize]
    d = a
    # v to be computed in ext_euclid_front()
    if b == 0:
        return (u, d)
    v3 = b
    # [Finished?]
    while v3 != 0:
        # [Euclidean step]
        q, t3 = divmod(d, v3)
        t1 = u - q * v1
        u = v1
        d = v3
        v1 = t1
        v3 = t3
    # [Finished?] cont. moved to ext_euclid_front()
    return (u, d)

# extended partial Euclidean algorithm
# reference: Sub-algorithm PARTEUCL(a, b) (p. 248) of Cohen -
# "A Course in Computational Algebraic Number theory" (GTM 138)
def part_euclid(d, v3, v, v2, L):
    # [Initialize]
    z = 0
    # [Finished?]
    while abs(v3) > L:
        # [Euclidean step]
        q, t3 = divmod_euclid(d, v3)
        t2 = v - q * v2
        v = v2
        d = v3
        v2 = t2
        v3 = t3
        z += 1
    # [Finished?] cont. moved to part_euclid_front()
    return (v, d, v2, v3, z)

# most significant digit of a, and the value of b in same place
# in base M (assumes a >= b, a >= 0, and b >= 0)
def same_msd(a, b, M):
    while a >= M:
        a //= M
        b //= M
    return a, b

# Lehmer extended (assumes a >= b, a >= 0, and b >= 0)
# reference: Algorithm 1.3.7 (p. 17) of Cohen -
# "A Course in Computational Algebraic Number theory" (GTM 138)
def lehmer(a, b, M):
    # [Initialize]
    u = 1
    v1 = 0
    # [Finished?]
    while b >= M:
        a_hat, b_hat = same_msd(a, b, M)
        A = 1
        B = 0
        C = 0
        D = 1
        # [Test quotient]
        while not (b_hat + C == 0 or b_hat + D == 0):
            q = (a_hat + A) // (b_hat + C)
            if q != ((a_hat + B) // (b_hat + D)):
                break
            # [Euclidean step]
            T = A - q * C
            A = C
            C = T
            T = B - q * D
            B = D
            D = T
            T = a_hat - q * b_hat
            a_hat = b_hat
            b_hat = T
        # [Multi-precision step]
        if B == 0:
            q, t = divmod(a, b)
            a = b
            b = t
            t = u - q * v1
            u = v1
            v1 = t
        else:
            t = A * a + B * b
            r = C * a + D * b
            a = t
            b = r
            t = A * u + B * v1
            r = C * u + D * v1
            u = t
            v1 = r
    return a, b, u, v1

# "frontend" for extended Euclidean algorithm
def ext_euclid_front(a, b, use_lehmer=True, M=(1 << 32)):
    # init: the algorithms assume that a >= 0 & b >= 0
    orig_a = a
    orig_b = b
    if orig_a < 0:
        a = -a
    if orig_b < 0:
        b = -b
    # execute algorithms
    if use_lehmer and a < b:
        at = a
        bt = b
        b, a, u, v1 = lehmer(b, a, M)
        u, d = ext_euclid(b, a, u, v1)
        v = u
        u = (d - bt * v) // at   # [Finished?] cont. of Euclid Extended
    elif use_lehmer:
        at = a
        bt = b
        a, b, u, v1 = lehmer(a, b, M)
        u, d = ext_euclid(a, b, u, v1)
        v = (d - at * u) // bt   # [Finished?] cont. of Euclid Extended
    else:
        u, d = ext_euclid(a, b)
        v = (d - a * u) // b   # [Finished?] cont. of Euclid Extended
    # final: check sign of orig a & b
    if orig_a < 0:
        a = -a
        u = -u
    if b < 0:
        b = -b
        v = -v
    return (u, v, d)

# "frontend" for extended partial Euclidean algorithm
def part_euclid_front(d, v3, v, v2, L):   # TODO: use_lehmer=True, M=(1 << 32)
    # execute algorithms
    v, d, v2, v3, z = part_euclid(d, v3, v, v2, L)
    # final: [Finished?] cont. of PARTEUCL
    if z % 2 == 1:
        v2 = -v2
        v3 = -v3
    return (v, d, v2, v3, z)