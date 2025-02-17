import unittest
from src.main import scoretype


class TestScore(unittest.TestCase):

    def test_percentile(self):
        data = (1.0, 3.0, 4.0, 5.0, 7.0, 8.0)
        self.assertEqual(
                1.0, scoretype.sorted_percentile(data, 0.0))
        self.assertEqual(
                8.0, scoretype.sorted_percentile(data, 1.0))
        self.assertEqual(
                4.5, scoretype.sorted_percentile(data, 0.5))
        self.assertEqual(
                3.25, scoretype.sorted_percentile(data, 0.25))
        self.assertEqual(
                6.5, scoretype.sorted_percentile(data, 0.75))
