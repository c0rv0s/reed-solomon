import unittest

from rs_grand_list_decoding.quotient_lower_bound_search import (
    find_prime_for_domain,
    generate_quotient_candidates,
    is_B_smooth,
    is_probable_prime,
    materialize_candidate,
    minimal_realizing_ell,
    quotient_candidate,
    rate_to_fraction,
    select_materialization_rows,
)


class QuotientLowerBoundSearchTests(unittest.TestCase):
    def test_rate_to_fraction(self):
        self.assertEqual(rate_to_fraction("1/4"), (1, 4))
        self.assertEqual(rate_to_fraction(0.25), (1, 4))

    def test_rho_quarter_M_225_candidate(self):
        rho_num, rho_den = rate_to_fraction("1/4")
        candidate = quotient_candidate(225, rho_num, rho_den)
        self.assertIsNotNone(candidate)
        assert candidate is not None
        self.assertEqual(candidate["r"], 57)
        self.assertAlmostEqual(candidate["radius"], 1 - 57 / 225)
        self.assertGreater(candidate["log2_list_lower_bound"], 170)
        ell = minimal_realizing_ell(225, rho_num, rho_den)
        self.assertEqual(ell, 4)
        self.assertEqual(ell * 225, 900)
        self.assertEqual((rho_num * ell * 225) // rho_den, 225)
        self.assertEqual(candidate["r"] * ell, 228)

    def test_power_of_two_quotient_collapses_for_dyadic_rate(self):
        rho_num, rho_den = rate_to_fraction("1/4")
        self.assertIsNone(quotient_candidate(256, rho_num, rho_den))

    def test_is_B_smooth(self):
        self.assertTrue(is_B_smooth(900, 5))
        self.assertFalse(is_B_smooth(900, 3))

    def test_generated_rows_include_M_225_q256_candidate(self):
        rows = generate_quotient_candidates(
            q_bits_values=(256,),
            eps_bits=128,
            rates=("1/4",),
            M_max=225,
            smoothness_bounds=(16,),
            n_scale_multipliers=(1,),
        )
        matches = [row for row in rows if row["M"] == 225 and row["q_bits"] == 256]
        self.assertEqual(len(matches), 1)
        row = matches[0]
        self.assertEqual(row["r"], 57)
        self.assertEqual(row["ell"], 4)
        self.assertEqual(row["n"], 900)
        self.assertEqual(row["k"], 225)
        self.assertEqual(row["s"], 228)
        self.assertTrue(row["beats_budget"])
        self.assertTrue(row["below_capacity"])
        self.assertGreater(row["log2_margin"], 50)
        self.assertEqual(row["M_smooth_B"], 5)
        self.assertEqual(row["n_smooth_B"], 5)

    def test_prime_search_small_case(self):
        p = find_prime_for_domain(n=12, q_bits=8, max_trials=1000)
        self.assertIsNotNone(p)
        assert p is not None
        self.assertEqual((p - 1) % 12, 0)
        self.assertEqual(p.bit_length(), 8)
        self.assertTrue(is_probable_prime(p))

    def test_materialize_candidate_smoke(self):
        rows = generate_quotient_candidates(
            q_bits_values=(32,),
            eps_bits=16,
            rates=("1/4",),
            M_max=45,
            smoothness_bounds=(16,),
            n_scale_multipliers=(1,),
        )
        row = rows[0]
        materialized = materialize_candidate(row, max_trials=10_000)
        self.assertEqual(materialized["n"], row["n"])
        self.assertEqual(materialized["s"], row["s"])
        self.assertIn("prime_found", materialized)

    def test_materialization_selection_includes_M_225_candidate(self):
        rows = generate_quotient_candidates(
            q_bits_values=(256,),
            eps_bits=128,
            rates=("1/16", "1/4"),
            M_max=225,
            smoothness_bounds=(128,),
            n_scale_multipliers=(1,),
        )
        selected = select_materialization_rows(rows, limit=1)
        self.assertTrue(any(row["M"] == 225 and row["rho_den"] == 4 for row in selected))


if __name__ == "__main__":
    unittest.main()
