"""
Copyright 2018 Ryan Wick (rrwick@gmail.com), Stephen Watts, Alex Tokolyi
https://github.com/rrwick/Pairwise-SNP-counter

This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version. This program is distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. You should
have received a copy of the GNU General Public License along with this program.  If not, see
<http://www.gnu.org/licenses/>.
"""

import unittest
import pairwise_snp_counter as psc
import pathlib


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

    def test_get_score_threshold_4(self):
        """
        If there are no scores, then the threshold is always zero.
        """
        scores = {}
        self.assertEqual(psc.get_score_threshold(scores, 0.0), 0.0)
        self.assertEqual(psc.get_score_threshold(scores, 5.0), 0.0)
        self.assertEqual(psc.get_score_threshold(scores, 20.0), 0.0)
        self.assertEqual(psc.get_score_threshold(scores, 100.0), 0.0)


class TestMpileupOutputParsing(unittest.TestCase):

    def test_get_base_scores_from_mpileup_output_1(self):
        """
        This mpileup output (mpileup_output_1) comes from an Illumina read assembly.
        """
        mpileup_output_filename = pathlib.Path(__file__).parent / 'mpileup_output_1'
        with open(mpileup_output_filename, 'rt') as mpileup_output_file:
            mpileup_output = mpileup_output_file.read()
        scores = psc.get_base_scores_from_mpileup_output(mpileup_output)

        # Test some individual bases.
        self.assertAlmostEqual(scores['chromosome'][0], 6 / 14)
        self.assertAlmostEqual(scores['chromosome'][1], 5 / 20)
        self.assertAlmostEqual(scores['chromosome'][50], 86 / 88)
        self.assertAlmostEqual(scores['chromosome'][51], 88 / 89)
        self.assertAlmostEqual(scores['chromosome'][52], 89 / 89)

        # These bases contain an indel which covers positions 57 and 58 (58 and 59 in the VCF's
        # 1-based indexing). This means that position 58 should not get a score of 1 as indicated
        # by its VCF line, but rather a score of 89/90 as indicated by the indel VCF line.
        self.assertAlmostEqual(scores['chromosome'][56], 90 / 90)
        self.assertAlmostEqual(scores['chromosome'][57], 89 / 90)
        self.assertAlmostEqual(scores['chromosome'][58], 89 / 90)
        self.assertAlmostEqual(scores['chromosome'][59], 90 / 91)

    def test_get_base_scores_from_mpileup_output_2(self):
        """
        This mpileup output (mpileup_output_2) comes from a Nanopore read assembly, where indels
        are much more common.
        """
        mpileup_output_filename = pathlib.Path(__file__).parent / 'mpileup_output_2'
        with open(mpileup_output_filename, 'rt') as mpileup_output_file:
            mpileup_output = mpileup_output_file.read()
        scores = psc.get_base_scores_from_mpileup_output(mpileup_output)

        # This base is (unusually for a Nanopore assembly) not covered by an indel.
        self.assertAlmostEqual(scores['chromosome'][29], 50 / 57)

        # These bases are covered by an indel as well as the normal substitution line, so should
        # have the worst score of the two.
        self.assertAlmostEqual(scores['chromosome'][30], 46 / 57)
        self.assertAlmostEqual(scores['chromosome'][50], 34 / 59)

        # This base is covered by two indels as well as the normal substitution line, so should
        # have the worst score of the three.
        self.assertAlmostEqual(scores['chromosome'][226], 30 / 60)
