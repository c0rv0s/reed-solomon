import math
import unittest

from rs_grand_list_decoding.rs_capacity_threshold import h_q, inverse_h_q, threshold_params


class CapacityThresholdTests(unittest.TestCase):
    def test_entropy_inverse_round_trip(self):
        ln_q = 256 * math.log(2)
        for delta in (0.1, 0.25, 0.5, 0.9):
            target = h_q(delta, ln_q)
            self.assertLess(abs(inverse_h_q(target, ln_q) - delta), 1e-10)

    def test_capacity_values_match_research_notes(self):
        expected = {
            1 / 2: 0.49609,
            1 / 4: 0.74683,
            1 / 8: 0.87288,
            1 / 16: 0.93618,
        }
        for rho, approx in expected.items():
            out = threshold_params(q_bits=256, n=2**20, rho=rho, m=1, mode="folded")
            self.assertLess(abs(float(out["delta_entropy"]) - approx), 5e-5)

    def test_m_improves_folded_threshold_toward_capacity(self):
        for rho in (1 / 2, 1 / 4, 1 / 8, 1 / 16):
            m1 = threshold_params(q_bits=256, n=2**20, rho=rho, m=1, mode="folded")
            m4 = threshold_params(q_bits=256, n=2**20, rho=rho, m=4, mode="folded")
            capacity = 1.0 - rho
            self.assertLess(
                abs(capacity - float(m4["delta_entropy"])),
                abs(capacity - float(m1["delta_entropy"])),
            )

    def test_volume_grid_does_not_exceed_entropy_by_more_than_one_step(self):
        out = threshold_params(q_bits=256, n=2**20, rho=1 / 4, m=4, mode="folded")
        one_step = 1 / int(out["N"])
        self.assertLessEqual(
            float(out["delta_volume_grid"]),
            float(out["delta_entropy"]) + one_step + 1e-12,
        )


if __name__ == "__main__":
    unittest.main()
