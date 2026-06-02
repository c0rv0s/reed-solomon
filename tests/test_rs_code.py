import itertools
import unittest

from rs_grand_list_decoding.finite_field import smooth_domain
from rs_grand_list_decoding.rs_code import (
    enumerate_rs_codewords,
    fold_codeword,
    hamming_distance,
    interleave_codewords,
)


class RSCodeTests(unittest.TestCase):
    def test_enumerate_rs_codeword_count(self):
        p = 7
        domain = smooth_domain(p, 6)
        codewords = enumerate_rs_codewords(p, domain, k=3)
        self.assertEqual(len(codewords), p**3)
        self.assertEqual(len(set(codewords)), p**3)

    def test_rs_minimum_distance(self):
        p = 7
        domain = smooth_domain(p, 6)
        k = 3
        codewords = enumerate_rs_codewords(p, domain, k=k)
        min_distance = min(
            hamming_distance(a, b) for a, b in itertools.combinations(codewords, 2)
        )
        self.assertEqual(min_distance, len(domain) - k + 1)

    def test_fold_codeword(self):
        self.assertEqual(fold_codeword((1, 2, 3, 4), 2), ((1, 2), (3, 4)))
        with self.assertRaises(ValueError):
            fold_codeword((1, 2, 3), 2)

    def test_interleave_codewords(self):
        self.assertEqual(
            interleave_codewords([(1, 2, 3), (4, 5, 6)]),
            ((1, 4), (2, 5), (3, 6)),
        )
        with self.assertRaises(ValueError):
            interleave_codewords([(1, 2), (3,)])


if __name__ == "__main__":
    unittest.main()
