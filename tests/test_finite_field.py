import unittest

from rs_grand_list_decoding.finite_field import factor, is_prime, primitive_root, smooth_domain


class FiniteFieldTests(unittest.TestCase):
    def test_is_prime(self):
        self.assertFalse(is_prime(1))
        self.assertTrue(is_prime(2))
        self.assertTrue(is_prime(17))
        self.assertFalse(is_prime(21))

    def test_factor(self):
        self.assertEqual(factor(1), {})
        self.assertEqual(factor(360), {2: 3, 3: 2, 5: 1})

    def test_primitive_root_generates_field_units(self):
        p = 17
        g = primitive_root(p)
        self.assertEqual({pow(g, i, p) for i in range(p - 1)}, set(range(1, p)))

    def test_smooth_domain_subgroup(self):
        domain = smooth_domain(17, 8)
        self.assertEqual(len(domain), 8)
        self.assertEqual(len(set(domain)), 8)
        self.assertTrue(all(x != 0 for x in domain))
        self.assertEqual({(a * b) % 17 for a in domain for b in domain}, set(domain))

    def test_smooth_domain_rejects_invalid_size(self):
        with self.assertRaises(ValueError):
            smooth_domain(17, 7)


if __name__ == "__main__":
    unittest.main()
