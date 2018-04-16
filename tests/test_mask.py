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
import pathlib
import pairwise_snp_counter as psc


class TestMaskFiles(unittest.TestCase):

    def setUp(self):
        self.assembly_fp = pathlib.Path(__file__).parent / 'temp'
        self.mask_fp = pathlib.Path(__file__).parent / 'temp.mask'
        self.scores = {'contig_1': [0.1, 0.1, 0.2, 0.2, 0.5, 0.0],
                       'contig_2': [0.8, 0.8, 0.0, 0.1, 0.2, 0.3]}

    def tearDown(self):
        if self.mask_fp.is_file():
            self.mask_fp.unlink()

    def test_mask_file_1(self):
        psc.write_mask_file(self.scores, 0.15, self.assembly_fp)
        target_mask = {'contig_1': [0, 1, 5], 'contig_2': [2, 3]}
        self.assertEqual(psc.load_mask_file(self.mask_fp), target_mask)

    def test_mask_file_2(self):
        """
        A threshold of 0.1 should not mask scores of 0.1 - scores must be lower than the threshold
        to be masked.
        """
        psc.write_mask_file(self.scores, 0.1, self.assembly_fp)
        target_mask = {'contig_1': [5], 'contig_2': [2]}
        self.assertEqual(psc.load_mask_file(self.mask_fp), target_mask)

    def test_mask_file_3(self):
        psc.write_mask_file(self.scores, 0.0, self.assembly_fp)
        target_mask = {'contig_1': [], 'contig_2': []}
        self.assertEqual(psc.load_mask_file(self.mask_fp), target_mask)

    def test_mask_file_4(self):
        psc.write_mask_file(self.scores, 1.0, self.assembly_fp)
        target_mask = {'contig_1': [0, 1, 2, 3, 4, 5], 'contig_2': [0, 1, 2, 3, 4, 5]}
        self.assertEqual(psc.load_mask_file(self.mask_fp), target_mask)
