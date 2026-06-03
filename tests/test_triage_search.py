import unittest

from rs_grand_list_decoding.finite_field import smooth_domain
from rs_grand_list_decoding.structured_bad_centers import (
    count_nearby_codewords_by_subsets,
    monomial_center,
)
from rs_grand_list_decoding.triage_search import (
    _domain_specs,
    add_low_k_fiber_classifier,
    add_sampling_metrics,
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
        self.assertEqual(fields["pattern_domain_group"], "smooth_coset")
        self.assertEqual(fields["a_mod_n"], 3)
        self.assertEqual(fields["b_mod_n"], 3)
        self.assertEqual(fields["b_minus_a_mod_n"], 0)
        self.assertEqual(fields["lambda"], 2)

    def test_pattern_fields_random_uses_raw_and_field_modulus(self):
        fields = pattern_fields(
            "binomial",
            "a=3;b=19;lambda=2",
            n=8,
            p=17,
            domain_type="random",
        )
        self.assertEqual(fields["pattern_domain_group"], "random")
        self.assertEqual(fields["a_raw"], 3)
        self.assertEqual(fields["b_raw"], 19)
        self.assertEqual(fields["b_mod_p_minus_1"], 3)
        self.assertEqual(fields["b_minus_a_raw"], 16)
        self.assertEqual(fields["b_minus_a_mod_p_minus_1"], 0)
        self.assertNotIn("b_mod_n", fields)

    def test_sampling_metrics(self):
        row = {
            "subsets_checked": 10,
            "subset_space": 100,
            "agreement_required": 3,
            "generic_ratio_raw": 7.0,
            "count": 5,
        }
        add_sampling_metrics(row, p=11, k=2)
        self.assertEqual(row["sample_fraction"], 0.1)
        self.assertEqual(row["sampled_generic_expected"], 10 * 11 ** -1)
        self.assertEqual(row["sampled_generic_ratio"], 5.0)
        self.assertEqual(row["full_generic_ratio_lower_bound"], 7.0)

    def test_low_k_fiber_classifier_explains_monomial_artifact(self):
        row = {
            "k": 1,
            "center_type": "monomial",
            "domain_type": "smooth",
            "center_parameters": "exponent=2",
            "n": 12,
            "agreement_required": 2,
            "count": 6,
        }
        add_low_k_fiber_classifier(row)
        self.assertEqual(row["low_k_fiber_size"], 2)
        self.assertEqual(row["low_k_fiber_predicted_count"], 6)
        self.assertTrue(row["low_k_fiber_explained"])

    def test_domain_specs_dedupe_coset_shift_inside_subgroup(self):
        p = 17
        n = 8
        subgroup_shift = smooth_domain(p, n)[1]
        specs = _domain_specs(p, n, random_seeds=(), coset_shifts=(subgroup_shift,))
        self.assertEqual([(kind, label) for kind, label, _ in specs], [("smooth", "subgroup")])

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
        self.assertIn("sampled_generic_ratio", rows[0])
        self.assertIn("low_k_fiber_explained", rows[0])

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
        smooth_coset_patterns = extract_pattern_rows(rows, min_count=1, domain_group="smooth_coset")
        random_patterns = extract_pattern_rows(rows, min_count=1, domain_group="random")
        self.assertTrue(comparisons)
        self.assertTrue(patterns)
        self.assertTrue(smooth_coset_patterns)
        self.assertTrue(random_patterns)


if __name__ == "__main__":
    unittest.main()
