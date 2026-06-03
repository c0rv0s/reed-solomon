import unittest

from rs_grand_list_decoding.lower_bound_constructions import (
    construct_coset_union_codewords,
    coset_union_lower_bound_count,
    cosets_of_subgroup,
    lower_bound_grid,
    subgroup_elements,
    subgroup_generator_for_size,
    subgroup_of_subgroup,
    vanishing_polynomial,
    verify_against_counting,
)
from rs_grand_list_decoding.rs_code import encode_rs
from rs_grand_list_decoding.structured_bad_centers import (
    count_nearby_codewords_by_subsets,
    poly_degree,
)


class LowerBoundConstructionTests(unittest.TestCase):
    def test_validated_binomial_row_is_coset_union_construction(self):
        p = 193
        n = 12
        k = 3
        s = 4
        ell = 2
        h_coeffs = [0, 0, 1]
        H = subgroup_elements(p, n)
        center = encode_rs(p, H, [0, 0, 1, 0, 1])

        self.assertEqual(coset_union_lower_bound_count(n, k, s, ell), 15)
        constructed = construct_coset_union_codewords(
            p,
            n,
            k,
            s,
            ell,
            h_coeffs=h_coeffs,
            max_emit=100,
        )
        self.assertEqual(len(constructed), 15)
        self.assertTrue(all(row["degree"] < k for row in constructed))
        self.assertTrue(all(row["selected_agreement_count"] == s for row in constructed))

        exact_count = count_nearby_codewords_by_subsets(
            p,
            H,
            k,
            center,
            radius=n - s,
        )
        self.assertEqual(exact_count, 15)

        verification = verify_against_counting(p, n, k, s, ell, h_coeffs=h_coeffs)
        self.assertEqual(verification["constructed_count"], 15)
        self.assertEqual(verification["exact_count"], 15)
        self.assertTrue(verification["matches_exact"])

    def test_validated_pure_monomial_center_also_matches_exact_count(self):
        verification = verify_against_counting(p=193, n=12, k=3, s=4, ell=2)
        self.assertEqual(verification["lower_bound_count"], 15)
        self.assertEqual(verification["constructed_count"], 15)
        self.assertEqual(verification["exact_count"], 15)

    def test_pair_construction_count(self):
        n = 16
        k = 7
        s = 8
        ell = 2
        self.assertEqual(coset_union_lower_bound_count(n, k, s, ell), 70)
        constructed = construct_coset_union_codewords(
            p=97,
            n=n,
            k=k,
            s=s,
            ell=ell,
            max_emit=100,
        )
        self.assertEqual(len(constructed), 70)
        self.assertTrue(all(row["degree"] < k for row in constructed))

    def test_reject_invalid_divisibility_cases(self):
        self.assertEqual(coset_union_lower_bound_count(n=12, k=3, s=5, ell=2), 0)
        self.assertEqual(coset_union_lower_bound_count(n=12, k=3, s=4, ell=5), 0)
        with self.assertRaises(ValueError):
            subgroup_of_subgroup(p=193, n=12, ell=5)
        with self.assertRaises(ValueError):
            construct_coset_union_codewords(p=193, n=12, k=3, s=5, ell=2)

    def test_distinct_coset_unions_emit_distinct_polynomials(self):
        constructed = construct_coset_union_codewords(
            p=193,
            n=12,
            k=3,
            s=4,
            ell=2,
            h_coeffs=[0, 0, 1],
            max_emit=100,
        )
        polynomials = {tuple(row["polynomial"]) for row in constructed}
        self.assertEqual(len(polynomials), len(constructed))

    def test_vanishing_polynomial_uses_only_ell_periodic_powers(self):
        p = 193
        n = 12
        ell = 2
        H = subgroup_elements(p, n)
        M = subgroup_of_subgroup(p, n, ell)
        cosets = cosets_of_subgroup(p, H, M)
        S = cosets[0] + cosets[1]
        coeffs = vanishing_polynomial(p, S)
        self.assertEqual(poly_degree(coeffs), len(S))
        for exponent, coeff in enumerate(coeffs):
            if exponent % ell != 0:
                self.assertEqual(coeff % p, 0)

    def test_subgroup_generator_has_expected_order(self):
        generator = subgroup_generator_for_size(193, 12)
        self.assertEqual(pow(generator, 12, 193), 1)
        self.assertNotEqual(pow(generator, 6, 193), 1)

    def test_lower_bound_grid_smoke(self):
        rows = lower_bound_grid(
            q_bits_values=(64,),
            eps_bits=128,
            n_values=(2**10,),
            rates=(1 / 2, 1 / 4),
        )
        self.assertTrue(rows)
        self.assertTrue(
            {
                "q_bits",
                "n",
                "k",
                "ell",
                "s",
                "radius",
                "log2_list_lower_bound",
                "beats_budget",
                "below_capacity",
            }.issubset(rows[0])
        )


if __name__ == "__main__":
    unittest.main()
