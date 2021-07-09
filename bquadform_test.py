# Tests for bquadform.py
# Reference: https://github.com/getamis/alice/blob/master/crypto/binaryquadraticform/binaryquadraticform_test.go

import unittest
from bquadform import BQuadForm

class TestBQuadForm(unittest.TestCase):
    def test_init(self):
        print("test_init: error messages of the following must be \"non-negative discriminant\"")
        with self.assertRaises(ValueError) as test: BQuadForm(0, 0, 5)
        print("test_init: " + str(test.exception))
        with self.assertRaises(ValueError) as test: BQuadForm(1, 10, 10)
        print("test_init: " + str(test.exception))
    def test_isReduced(self):
        # output of __init__ is reduced, but for this test, we force a non-reduced form.
        test = BQuadForm(33, 11, 5)
        test.a = 33
        test.b = 11
        test.c = 5
        self.assertFalse(test.isReduced())
    def test_get(self):
        a = 101
        b = 38
        c = 4898
        test = BQuadForm(a, b, c)
        self.assertEqual(test.a, a)
        self.assertEqual(test.b, b)
        self.assertEqual(test.c, c)
        # discriminant is not tested here, unlike in getamis/alice
    def test_reduction(self):
        test = BQuadForm(33, 11, 5)
        self.assertEqual((test.a, test.b, test.c), (5, -1, 27))
        test = BQuadForm(15, 0, 15)
        self.assertEqual((test.a, test.b, test.c), (15, 0, 15))
        test = BQuadForm(6, 3, 1)
        self.assertEqual((test.a, test.b, test.c), (1, 1, 4))
        test = BQuadForm(1, 2, 3)
        self.assertEqual((test.a, test.b, test.c), (1, 0, 2))
        test = BQuadForm(1, 2, 30)
        self.assertEqual((test.a, test.b, test.c), (1, 0, 29))
        test = BQuadForm(4, 5, 3)
        self.assertEqual((test.a, test.b, test.c), (2, -1, 3))
    def test_compose(self):
        test = BQuadForm(1, 1, 6) * BQuadForm(1, 1, 6)
        self.assertEqual((test.a, test.b, test.c), (1, 1, 6))
        test = BQuadForm(2, -1, 3) * BQuadForm(2, -1, 3)
        self.assertEqual((test.a, test.b, test.c), (2, 1, 3))
        test = BQuadForm(2, -1, 3) * BQuadForm(2, 1, 3)
        self.assertEqual((test.a, test.b, test.c), (1, 1, 6))
        test = BQuadForm(31, 24, 15951) * BQuadForm(31, 24, 15951)
        self.assertEqual((test.a, test.b, test.c), (517, 100, 961))
        test = BQuadForm(142, 130, 3511) * BQuadForm(677, 664, 893)
        self.assertEqual((test.a, test.b, test.c), (591, 564, 971))
    def test_square(self):
        test = BQuadForm(1, 1, 6).square()
        self.assertEqual((test.a, test.b, test.c), (1, 1, 6))
        test = BQuadForm(19, 18, 26022).square()
        self.assertEqual((test.a, test.b, test.c), (361, -286, 1426))
        test = BQuadForm(19, -12, 262).square()
        self.assertEqual((test.a, test.b, test.c), (46, -32, 113))
        test = BQuadForm(31, 24, 15951).square()
        self.assertEqual((test.a, test.b, test.c), (517, 100, 961))
        test = BQuadForm(517, 100, 961).square()
        self.assertEqual((test.a, test.b, test.c), (529, -378, 1002))
        test = BQuadForm(3, -2, 176081).square()
        self.assertEqual((test.a, test.b, test.c), (9, 4, 58694))
        test = BQuadForm(729, 626, 859).square()
        self.assertEqual((test.a, test.b, test.c), (419, -412, 1362))
    def test_exp(self):
        test = BQuadForm(2, 1, 3) ** 6
        self.assertEqual((test.a, test.b, test.c), (1, 1, 6))
        test = BQuadForm(31, 24, 15951) ** 200
        self.assertEqual((test.a, test.b, test.c), (517, -276, 993))
        test = BQuadForm(78, -52, 6781) ** 500
        self.assertEqual((test.a, test.b, test.c), (738, -608, 841))
        test = BQuadForm(101, 38, 4898) ** 508
        self.assertEqual((test.a, test.b, test.c), (66, 54, 7501))
        test = BQuadForm(101, 38, 4898) ** 1
        self.assertEqual((test.a, test.b, test.c), (101, 38, 4898))
        test = BQuadForm(101, 38, 4898) ** 22999971
        self.assertEqual((test.a, test.b, test.c), (101, 38, 4898))

if __name__ == '__main__':
    unittest.main()
