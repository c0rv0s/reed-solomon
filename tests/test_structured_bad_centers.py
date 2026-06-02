import unittest

import math

from rs_grand_list_decoding.exact_search import (
    build_folded_rs_instance,
    build_scalar_rs_instance,
    exact_profile,
)
from rs_grand_list_decoding.finite_field import random_domain, smooth_domain
from rs_grand_list_decoding.rs_code import encode_rs, fold_codeword, hamming_distance
from rs_grand_list_decoding.structured_bad_centers import (
    baseline_metrics,
    count_nearby_codewords_bruteforce,
    count_nearby_codewords_by_subsets,
    count_nearby_folded_codewords_bruteforce,
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
        self.assertTrue(
            {
                "p",
                "n",
                "k",
                "rho",
                "mode",
                "m",
                "N",
                "domain_type",
                "domain_label",
                "radius",
                "relative_radius",
                "agreement_required",
                "subset_count",
                "generic_expected",
                "generic_ratio",
                "boundary_ratio",
                "boundary_case",
                "center_type",
                "center_parameters",
                "count",
                "duplicate_parameter_count",
                "all_parameters_sample",
                "method",
            }.issubset(set(rows[0])),
        )

    def test_duplicate_monomial_exponents_dedupe(self):
        rows = search_monomial_centers(
            p=5,
            n=4,
            k=2,
            radius=2,
            exponent_range=[0, 4],
        )
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["duplicate_parameter_count"], 2)

    def test_baseline_metrics(self):
        metrics = baseline_metrics(p=5, N=4, k=2, radius=1, m=1)
        s = 3
        self.assertEqual(metrics["agreement_required"], s)
        self.assertEqual(metrics["subset_count"], math.comb(4, s))
        self.assertEqual(metrics["generic_expected"], math.comb(4, s) * 5 ** (2 - s))
        self.assertFalse(metrics["boundary_case"])

    def test_boundary_case_marked_at_s_equals_k(self):
        rows = search_monomial_centers(p=5, n=4, k=2, radius=2, exponent_range=[0])
        self.assertTrue(rows[0]["boundary_case"])

    def test_smooth_and_random_domains_are_valid(self):
        smooth = smooth_domain(17, 8)
        random = random_domain(17, 8, seed=1)
        self.assertEqual(len(smooth), 8)
        self.assertEqual(len(random), 8)
        self.assertEqual(len(set(smooth)), 8)
        self.assertEqual(len(set(random)), 8)
        self.assertTrue(all(x != 0 for x in smooth + random))

    def test_folded_count_agrees_with_exact_folded_bruteforce(self):
        p = 5
        n = 4
        k = 2
        m = 2
        domain = smooth_domain(p, n)
        center = monomial_center(p, domain, 3)
        radius = 1
        count = count_nearby_folded_codewords_bruteforce(
            p, domain, k=k, center=center, radius=radius, m=m
        )
        instance = build_folded_rs_instance(p=p, n=n, k=k, m=m)
        folded_center = fold_codeword(center, m)
        exact_count = sum(
            hamming_distance(folded_center, codeword) <= radius
            for codeword in instance["codewords"]
        )
        self.assertEqual(count, exact_count)


if __name__ == "__main__":
    unittest.main()
