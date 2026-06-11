import unittest

from rs_grand_list_decoding.quotient_lower_bound_search import (
    FLAGSHIP_P,
    PROTH_K,
    PROTH_N,
    PROTH_WITNESS,
    flagship_counterexample,
    find_prime_for_domain,
    generate_quotient_candidates,
    is_B_smooth,
    is_probable_prime,
    materialize_candidate,
    materialize_top_candidates,
    minimal_realizing_ell,
    quotient_candidate,
    rate_to_fraction,
    select_materialization_rows,
    smoothness_level,
    verify_proth_certificate,
    verify_flagship_prime,
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
        materialized = materialize_candidate(
            225,
            rho="1/4",
            ell_multiplier=1,
            q_bits=256,
            eps_bits=128,
            p=FLAGSHIP_P,
        )
        self.assertGreater(materialized["log2_list_lower_bound"], 179)
        self.assertGreater(materialized["log2_margin"], 51)
        self.assertTrue(materialized["beats_budget"])
        self.assertTrue(materialized["below_capacity"])

    def test_power_of_two_quotient_collapses_for_dyadic_rate(self):
        rho_num, rho_den = rate_to_fraction("1/4")
        self.assertIsNone(quotient_candidate(224, rho_num, rho_den))

    def test_is_B_smooth(self):
        self.assertTrue(is_B_smooth(900, 5))
        self.assertFalse(is_B_smooth(900, 3))
        self.assertEqual(smoothness_level(900, bounds=(2, 4, 8, 16)), 8)
        self.assertIsNone(smoothness_level(900, bounds=(2, 4)))

    def test_generated_rows_include_M_225_q256_candidate(self):
        rows = generate_quotient_candidates(
            q_bits_values=(256,),
            eps_bits=128,
            rates=("1/4",),
            M_max=225,
            smoothness_bounds=(8, 16),
            ell_multipliers=(1,),
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
        self.assertEqual(row["M_smooth_B"], 8)
        self.assertEqual(row["n_smooth_B"], 8)

    def test_prime_search_small_case(self):
        p = find_prime_for_domain(n=12, q_bits=8, max_trials=1000)
        self.assertIsNotNone(p)
        assert p is not None
        self.assertEqual((p - 1) % 12, 0)
        self.assertEqual(p.bit_length(), 8)
        self.assertTrue(is_probable_prime(p))

    def test_hardcoded_flagship_prime(self):
        verification = verify_flagship_prime()
        self.assertEqual(FLAGSHIP_P, PROTH_K * 2**PROTH_N + 1)
        self.assertEqual(PROTH_K % 2, 1)
        self.assertLess(PROTH_K, 2**PROTH_N)
        self.assertEqual(pow(PROTH_WITNESS, (FLAGSHIP_P - 1) // 2, FLAGSHIP_P), FLAGSHIP_P - 1)
        self.assertTrue(
            verify_proth_certificate(FLAGSHIP_P, PROTH_K, PROTH_N, PROTH_WITNESS)
        )
        self.assertTrue(verification["proth_certificate_verified"])
        self.assertTrue(verification["prime_verified"])
        self.assertEqual(verification["p_mod_n"], 1)
        self.assertEqual(FLAGSHIP_P % 900, 1)
        self.assertEqual(FLAGSHIP_P.bit_length(), 256)

    def test_flagship_counterexample(self):
        row = flagship_counterexample()
        self.assertEqual(row["p"], FLAGSHIP_P)
        self.assertEqual(row["n"], 900)
        self.assertEqual(row["k"], 225)
        self.assertEqual(row["ell"], 4)
        self.assertEqual(row["s"], 228)
        self.assertEqual(row["M"], 225)
        self.assertEqual(row["r"], 57)
        self.assertTrue(row["prime_verified"])
        self.assertEqual(row["p_mod_n"], 1)
        self.assertTrue(row["beats_budget"])
        self.assertTrue(row["below_capacity"])

    def test_materialize_candidate_smoke(self):
        rows = generate_quotient_candidates(
            q_bits_values=(32,),
            eps_bits=16,
            rates=("1/4",),
            M_max=45,
            smoothness_bounds=(16,),
            ell_multipliers=(1,),
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
            ell_multipliers=(1,),
        )
        selected = select_materialization_rows(rows, limit=1)
        self.assertTrue(any(row["M"] == 225 and row["rho_den"] == 4 for row in selected))
        materialized = materialize_top_candidates(rows, limit=1, max_trials=1)
        self.assertTrue(any(row["M"] == 225 and row["n"] == 900 for row in materialized))


if __name__ == "__main__":
    unittest.main()
