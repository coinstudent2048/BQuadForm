# Primitive positive definite binary quadratic forms
#
# Use this code only for prototyping

from bquadform_utils import isqrt, divmod_euclid, ext_euclid, part_euclid

class BQuadForm:
    def __init__(self, a, b, c = None, disc = None):
        if (c == None and disc == None) or (c != None and disc != None):
            raise AttributeError("invalid arguments")
        if not (isinstance(a, int) and isinstance(b, int)):
            raise TypeError("'a' or 'b' is not integer")
        # optional: check if a < 0 to detect non-positive-definiteness

        if c != None:
            if not isinstance(c, int):
                raise TypeError("'c' is not integer")
            # c is given, so compute discriminant
            disc = b * b - ((a * c) << 2)
            if disc >= 0:
                raise ValueError("non-negative discriminant")

        elif disc != None:
            if not isinstance(disc, int):
                raise TypeError("'disc' is not integer")
            if disc >= 0:
                raise ValueError("non-negative discriminant")
            # disc is given, so compute c
            c, c_test = divmod(b * b - disc, a << 2)   # a is assumed positive
            if c_test != 0:
                raise TypeError("'c' is not integer")

        # optional: check if gcd(a, b, c) != 1 to detect non-primitiveness
        self.a = a
        self.b = b
        self.c = c
        self.disc = disc
        self.L = isqrt(isqrt(abs(disc) >> 2))   # 4th root of |disc/4|
        # only work with reduced forms
        self.__reduction()

    # reduction of positive definite forms
    # reference: Algorithm 5.4.2 (p.243) of Cohen -
    # "A Course in Computational Algebraic Number theory" (GTM 138)
    def __reduction(self):
        # [Initialize]: if -a < b <= a goto [Finished?]
        if -self.a < self.b and self.b <= self.a:
            self.__reduction_finished()
            return
        self.__reduction_euclidean()
        self.__reduction_finished()

    # [Euclidean step] in reduction of positive definite forms
    def __reduction_euclidean(self):
        # let b = 2aq + r with 0 <= r < 2a be the Euclidean division of b by 2a
        q, r = divmod_euclid(self.b, self.a << 1)
        # if r > a, set r = r - 2a, q = q + 1 (i.e. we want b = 2aq + r, where -a < r <= a)
        if r > self.a:
            r -= self.a << 1
            q += 1
        # then set c = c - (1/2)(b + r)q, b = r
        self.c -= ((self.b + r) * q) >> 1
        self.b = r

    # [Finished?] in reduction of positive definite forms
    def __reduction_finished(self):
        while not self.isReduced():
            # if a > c set b = -b, exchange a and c and goto [Euclidean step]
            if self.a > self.c:
                self.b = -self.b
                self.a, self.c = self.c, self.a
                self.__reduction_euclidean()
            # otherwise, if a = c and b < 0, set b = -b
            elif self.a == self.c and self.b < 0:
                self.b = -self.b

    # check if the positive definite form is reduced
    # reference: Definition 5.3.2 (p.231) of Cohen -
    # "A Course in Computational Algebraic Number theory" (GTM 138)
    def isReduced(self):
        abs_b = abs(self.b)
        # |b| < a < c
        if abs_b < self.a and self.a < self.c:
            return True
        # |b| = a and b >= 0
        if abs_b == self.a and self.b >= 0:
            return True
        # a = c and b >= 0
        if self.a == self.c and self.b >= 0:
            return True
        return False

    # equality
    def __eq__(self, y):
        if isinstance(y, BQuadForm):
            return self.a == y.a and self.b == y.b and self.c == y.c and self.disc == y.disc
        raise TypeError("syntax error")

    # string representation
    def __repr__(self):
        return str((self.a, self.b, self.c)) + " disc: " + str(self.disc)

    # inverse of primitive quadratic form
    # reference: p.250 of Cohen -
    # "A Course in Computational Algebraic Number theory" (GTM 138)
    def inverse(self):
        return BQuadForm(self.a, -self.b, self.c)

    # identity element := BQuadForm * BQuadForm.inverse()
    def identity(self):
        return self * self.inverse()

    # composition (NUCOMP algorithm)
    # reference: Algorithm 5.4.9 (p.249) of Cohen -
    # "A Course in Computational Algebraic Number theory" (GTM 138)
    def __mul__(self, y):
        if not isinstance(y, BQuadForm):
            raise TypeError("syntax error")
        if self.disc != y.disc:
            raise ValueError("discriminants are different")
        # [Initialize]
        if self.a > y.a:
            a1 = self.a
            b1 = self.b
            c1 = self.c
            a2 = y.a
            b2 = y.b
            c2 = y.c
        else:    # if a1 < a2 exchange f1 and f2
            a1 = y.a
            b1 = y.b
            c1 = y.c
            a2 = self.a
            b2 = self.b
            c2 = self.c
        s = (b1 + b2) >> 1
        n = b2 - s
        # [First Euclidean step]
        u, v ,d = ext_euclid(a2, a1)
        if s % d == 0:   # if d|s
            A = -u * n
            d1 = d
            if d != 1:
                a1 //= d1
                a2 //= d1
                s //= d1
        else:   # if not d|s
            # [Second Euclidean step]
            u1, v1, d1 = ext_euclid(s, d)
            if d1 > 1:
                a1 //= d1
                a2 //= d1
                s //= d1
                d //= d1
            # [Initialization of reduction]
            l = (-u1 * (u * (c1 % d) + v * (c2 % d))) % d
            A = -u * (n // d) + l * (a1 // d)
        # [Partial reduction]
        A %= a1
        A1 = a1 - A
        if A1 < A:
            A = -A1
        # input sequence: d, v3, v, v2, L = bound
        v, d, v2, v3, z = part_euclid(a1, A, 0, 1, self.L)
        if z % 2 == 1:   # [Finished?] cont. of PARTEUCL
            v2 = -v2
            v3 = -v3
        # [Special case]
        if z == 0:
            Q1 = a2 * v3
            Q2 = Q1 + n
            f = Q2 // d
            g = (v3 * s + c2) // d
            a3 = d * a2
            # c3 = v3 * f + g * d1   # not required if disc is provided
            b3 = 2 * Q1 + b2
            return BQuadForm(a3, b3, disc = self.disc)
        # [Final computations]
        b = (a2 * d + n * v) // a1
        Q1 = b * v3
        Q2 = Q1 + n
        f = Q2 // d
        e = (s * d + c2 * v) // a1
        Q3 = e * v2
        Q4 = Q3 - s
        g = Q4 // v
        if d1 > 1:
            v2 *= d1
            v *= d1
        a3 = d * b + e * v
        # c3 = v3 * f + g * v2   # not required if disc is provided
        b3 = Q1 + Q2 + d1 * (Q3 + Q4)
        return BQuadForm(a3, b3, disc = self.disc)

    # squaring (NUDUPL algorithm)
    # reference: Algorithm 5.4.8 (p.248) of Cohen -
    # "A Course in Computational Algebraic Number theory" (GTM 138)
    def square(self):
        a = self.a
        b = self.b
        c = self.c
        # [Euclidean step]
        u, _, d1 = ext_euclid(b, a)
        A = a // d1
        B = b // d1
        C = (- c * u) % A
        C1 = A - C
        if C1 < C:
            C = -C1
        # [Partial reduction]
        # input sequence: d, v3, v, v2, L = bound
        v, d, v2, v3, z = part_euclid(A, C, 0, 1, self.L)
        if z % 2 == 1:   # [Finished?] cont. of PARTEUCL
            v2 = -v2
            v3 = -v3
        # [Special case]
        if z == 0:
            g = (B * v3 + c) // d
            a2 = d * d
            # c2 = v3 * v3 + g * d1   # not required if disc is provided
            b2 = b + 2 * d * v3   # simplified
            return BQuadForm(a2, b2, disc = self.disc)
        # [Final computations]
        e = (c * v + B * d) // A
        g = (e * v2 - B) // v
        b2 = e * v2 + v * g
        if d1 > 1:
            b2 *= d1
            v *= d1
            v2 *= d1
        a2 = d * d + e * v
        # c2 = v3 * v3 + g * v2   # not required if disc is provided
        b2 += 2 * d * v3   # simplified
        return BQuadForm(a2, b2, disc = self.disc)

    # exponentiation using Non-Adjacent Form (NAF) of an integer
    # reference: Algorithm 3.1 (p.27) of Sayles -
    # "Improved Arithmetic in the Ideal Class Group of Imaginary Quadratic Number Fields
    # with an Application to Integer Factoring"
    def __pow__(self, y):
        if not (isinstance(y, int) and y >= 0):
            raise TypeError("syntax error")
        c = 0
        T = self
        R = self.identity()
        # note: instead of using counter i for the while loop, I used right-shifting of y instead
        rsh_y = y
        while rsh_y > 0:
            twobit_test = (rsh_y + c) % 4
            if twobit_test == 1:
                R *= T
                c = 0
            elif twobit_test == 3:
                R *= T.inverse()
                c = 1
            T = T.square()
            rsh_y >>= 1
        if c == 1:
            R *= T
        return R
