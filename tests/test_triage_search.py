import unittest

from rs_grand_list_decoding.finite_field import smooth_domain
from rs_grand_list_decoding.structured_bad_centers import (
    count_nearby_codewords_by_subsets,
    monomial_center,
)
from rs_grand_list_decoding.triage_search import (
    build_triage_comparison_rows,
    count_by_sampled_subsets,
    extract_pattern_rows,
    generate_triage_rows,
    pattern_fields,
    sampled_index_subsets,
    stable_seed,
)


class TriageSearchTests(unittest.TestCase):
    def test_sampled_index_subsets_exact_when_budget_suffices(self):
        subsets, exact = sampled_index_subsets(N=5, s=2, sample_budget=20, seed=1)
        self.assertTrue(exact)
        self.assertEqual(len(subsets), 10)

    def test_sampled_index_subsets_deterministic_when_sampled(self):
        a, exact_a = sampled_index_subsets(N=12, s=6, sample_budget=15, seed=7)
        b, exact_b = sampled_index_subsets(N=12, s=6, sample_budget=15, seed=7)
        self.assertFalse(exact_a)
        self.assertFalse(exact_b)
        self.assertEqual(a, b)

    def test_sampled_count_matches_exact_when_full_subset_space_scanned(self):
        p = 5
        domain = smooth_domain(p, 4)
        center = monomial_center(p, domain, 3)
        exact = count_nearby_codewords_by_subsets(p, domain, k=2, center=center, radius=2)
        sampled = count_by_sampled_subsets(
            p,
            domain,
            k=2,
            center=center,
            radius=2,
            sample_budget=100,
            seed=0,
        )
        self.assertTrue(sampled["exact_subset_scan"])
        self.assertEqual(sampled["count_lower_bound"], exact)

    def test_pattern_fields(self):
        fields = pattern_fields("binomial", "a=3;b=11;lambda=2", n=8)
        self.assertEqual(fields["a_mod_n"], 3)
        self.assertEqual(fields["b_mod_n"], 3)
        self.assertEqual(fields["b_minus_a_mod_n"], 0)
        self.assertEqual(fields["lambda"], 2)

    def test_stable_seed(self):
        self.assertEqual(stable_seed("a", 1), stable_seed("a", 1))
        self.assertNotEqual(stable_seed("a", 1), stable_seed("a", 2))

    def test_generate_triage_rows_smoke(self):
        rows = generate_triage_rows(
            primes=(17,),
            max_n=4,
            random_seed_count=2,
            agreement_offsets=(1,),
            exponent_limit=3,
            lambda_limit=1,
            constant_limit=1,
            sample_budget=20,
        )
        self.assertTrue(rows)
        self.assertTrue({"smooth", "random"}.issubset({row["domain_type"] for row in rows}))
        self.assertIn("generic_ratio_raw", rows[0])

    def test_comparison_and_patterns_smoke(self):
        rows = generate_triage_rows(
            primes=(17,),
            max_n=4,
            random_seed_count=2,
            agreement_offsets=(1,),
            exponent_limit=3,
            lambda_limit=1,
            constant_limit=1,
            sample_budget=20,
        )
        comparisons = build_triage_comparison_rows(rows)
        patterns = extract_pattern_rows(rows, min_count=1)
        self.assertTrue(comparisons)
        self.assertTrue(patterns)


if __name__ == "__main__":
    unittest.main()
