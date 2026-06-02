import unittest

from rs_grand_list_decoding.exact_search import build_scalar_rs_instance, exact_profile
from rs_grand_list_decoding.finite_field import smooth_domain
from rs_grand_list_decoding.rs_code import encode_rs
from rs_grand_list_decoding.structured_bad_centers import (
    count_nearby_codewords_bruteforce,
    count_nearby_codewords_by_subsets,
    lagrange_interpolate,
    monomial_center,
    poly_degree,
    search_monomial_centers,
    sparse_center,
)


class StructuredBadCenterTests(unittest.TestCase):
    def test_monomial_center(self):
        domain = smooth_domain(5, 4)
        self.assertEqual(monomial_center(5, domain, 3), tuple(pow(x, 3, 5) for x in domain))

    def test_sparse_center(self):
        domain = [1, 2, 3]
        center = sparse_center(7, domain, [(2, 1), (3, 2)])
        self.assertEqual(center, tuple((2 * x + 3 * x * x) % 7 for x in domain))

    def test_lagrange_interpolate_recovers_polynomial(self):
        p = 7
        xs = [1, 2, 3]
        coeffs = [3, 2, 1]
        ys = [sum(coeff * pow(x, i, p) for i, coeff in enumerate(coeffs)) % p for x in xs]
        recovered = lagrange_interpolate(p, xs, ys)
        self.assertEqual(recovered, coeffs)
        self.assertEqual(poly_degree(recovered), 2)

    def test_bruteforce_and_subset_counts_agree(self):
        p = 5
        domain = smooth_domain(p, 4)
        center = monomial_center(p, domain, 3)
        brute = count_nearby_codewords_bruteforce(p, domain, k=2, center=center, radius=2)
        subsets = count_nearby_codewords_by_subsets(p, domain, k=2, center=center, radius=2)
        self.assertEqual(subsets, brute)

    def test_codeword_center_radius_zero_count_is_one(self):
        p = 5
        domain = smooth_domain(p, 4)
        center = encode_rs(p, domain, [2, 3])
        self.assertEqual(
            count_nearby_codewords_by_subsets(p, domain, k=2, center=center, radius=0),
            1,
        )

    def test_structured_search_does_not_exceed_exact_profile(self):
        instance = build_scalar_rs_instance(p=5, n=4, k=2)
        profile = exact_profile(instance["codewords"], instance["alphabet"])
        exact_max = profile[2]["max_list"]
        rows = search_monomial_centers(p=5, n=4, k=2, radius=2)
        self.assertLessEqual(max(row["count"] for row in rows), exact_max)

    def test_search_report_shape(self):
        rows = search_monomial_centers(p=5, n=4, k=2, radius=2, exponent_range=range(2))
        self.assertEqual(
            {
                "p",
                "n",
                "k",
                "rho",
                "radius",
                "relative_radius",
                "center_type",
                "center_parameters",
                "count",
                "method",
            },
            set(rows[0]),
        )


if __name__ == "__main__":
    unittest.main()
