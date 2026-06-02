import unittest

from rs_grand_list_decoding.exact_search import (
    build_folded_rs_instance,
    build_interleaved_rs_instance,
    build_scalar_rs_instance,
    exact_profile,
)


class ExactSearchTests(unittest.TestCase):
    def test_repetition_code_profile(self):
        code = [(0, 0, 0), (1, 1, 1)]
        profile = exact_profile(code, [0, 1])
        self.assertEqual([row["max_list"] for row in profile], [1, 1, 2, 2])

    def test_tiny_scalar_rs_profile(self):
        instance = build_scalar_rs_instance(p=5, n=4, k=2)
        self.assertEqual(len(instance["codewords"]), 25)
        profile = exact_profile(instance["codewords"], instance["alphabet"])
        self.assertEqual(profile[0]["max_list"], 1)
        max_lists = [row["max_list"] for row in profile]
        self.assertEqual(max_lists, sorted(max_lists))

    def test_tiny_scalar_rs_min_distance_boundary(self):
        instance = build_scalar_rs_instance(p=5, n=4, k=2)
        profile = exact_profile(instance["codewords"], instance["alphabet"])
        self.assertEqual(profile[2]["max_list"], 6)
        self.assertEqual(profile[3]["max_list"], 17)

    def test_scalar_rs_7_6_3_boundary(self):
        instance = build_scalar_rs_instance(p=7, n=6, k=3)
        profile = exact_profile(
            instance["codewords"],
            instance["alphabet"],
            max_operations=50_000_000,
        )
        self.assertEqual(profile[3]["max_list"], 20)
        self.assertEqual(profile[4]["max_list"], 91)

    def test_exact_profile_refuses_large_received_space(self):
        code = [(0,) * 10, (1,) * 10]
        with self.assertRaises(ValueError):
            exact_profile(code, [0, 1], max_words=100)

    def test_folded_smoke(self):
        instance = build_folded_rs_instance(p=5, n=4, k=2, m=2)
        self.assertEqual(len(instance["codewords"]), 25)
        self.assertEqual(len(instance["codewords"][0]), 2)

    def test_interleaved_smoke(self):
        instance = build_interleaved_rs_instance(p=5, n=4, k=1, m=2)
        self.assertEqual(len(instance["codewords"]), 25)
        self.assertEqual(len(instance["codewords"][0]), 4)


if __name__ == "__main__":
    unittest.main()
