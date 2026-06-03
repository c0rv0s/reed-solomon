import unittest

from rs_grand_list_decoding.finite_field import smooth_domain
from rs_grand_list_decoding.structured_bad_centers import monomial_center
from rs_grand_list_decoding.triage_validation import (
    center_from_row,
    parse_int_list,
    validate_candidate,
)


class TriageValidationTests(unittest.TestCase):
    def test_parse_int_list(self):
        self.assertEqual(parse_int_list("100, 1000 10000"), [100, 1000, 10000])

    def test_center_from_row_reconstructs_monomial(self):
        row = {"center_type": "monomial", "center_parameters": "exponent=2"}
        domain = smooth_domain(17, 4)
        self.assertEqual(center_from_row(17, domain, row), monomial_center(17, domain, 2))

    def test_validate_candidate_smoke(self):
        row = {
            "p": "17",
            "n": "4",
            "k": "1",
            "radius": "2",
            "center_type": "monomial",
            "center_parameters": "exponent=2",
            "count": "2",
            "generic_ratio": "2.0",
            "generic_ratio_raw": "2.0",
            "boundary_case": "False",
        }
        result = validate_candidate(
            row,
            budgets=(10,),
            random_seed_count=3,
            exact_threshold=1000,
        )
        self.assertEqual(result["subset_space"], 6)
        self.assertTrue(result["exact_subset_scan"])
        self.assertEqual(result["smooth_count_budget_10"], result["smooth_best_count"])
        self.assertEqual(result["random_seed_count"], 3)


if __name__ == "__main__":
    unittest.main()
