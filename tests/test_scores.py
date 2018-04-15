import unittest


import pairwise_snp_counter as psc


class TestScoreThreshold(unittest.TestCase):

    def test_get_score_threshold_1(self):
        """
        Example 1 from https://en.wikipedia.org/wiki/Percentile#The_nearest-rank_method
        """
        scores = {'contig_1': [35, 15, 50], 'contig_2': [40, 20]}
        self.assertEqual(psc.get_score_threshold(scores, 0.0), 15)
        self.assertEqual(psc.get_score_threshold(scores, 5.0), 15)
        self.assertEqual(psc.get_score_threshold(scores, 30.0), 20)
        self.assertEqual(psc.get_score_threshold(scores, 40.0), 20)
        self.assertEqual(psc.get_score_threshold(scores, 50.0), 35)
        self.assertEqual(psc.get_score_threshold(scores, 100.0), 50)

    def test_get_score_threshold_2(self):
        """
        Example 2 from https://en.wikipedia.org/wiki/Percentile#The_nearest-rank_method
        """
        scores = {'contig_1': [13, 15, 6, 7, 8, 20, 3, 8, 10, 16]}
        self.assertEqual(psc.get_score_threshold(scores, 0.0), 3)
        self.assertEqual(psc.get_score_threshold(scores, 25.0), 7)
        self.assertEqual(psc.get_score_threshold(scores, 50.0), 8)
        self.assertEqual(psc.get_score_threshold(scores, 75.0), 15)
        self.assertEqual(psc.get_score_threshold(scores, 100.0), 20)

    def test_get_score_threshold_3(self):
        """
        Example 3 from https://en.wikipedia.org/wiki/Percentile#The_nearest-rank_method
        """
        scores = {'contig_1': [3, 8], 'contig_2': [6, 9], 'contig_3': [7, 10], 'contig_4': [8, 13],
                  'contig_5': [20], 'contig_6': [16], 'contig_7': [15]}
        self.assertEqual(psc.get_score_threshold(scores, 0.0), 3)
        self.assertEqual(psc.get_score_threshold(scores, 25.0), 7)
        self.assertEqual(psc.get_score_threshold(scores, 50.0), 9)
        self.assertEqual(psc.get_score_threshold(scores, 75.0), 15)
        self.assertEqual(psc.get_score_threshold(scores, 100.0), 20)

