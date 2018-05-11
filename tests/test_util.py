"""
Copyright 2018 Ryan Wick (rrwick@gmail.com), Stephen Watts, Alex Tokolyi
https://github.com/rrwick/Snouter

This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version. This program is distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. You should
have received a copy of the GNU General Public License along with this program.  If not, see
<http://www.gnu.org/licenses/>.
"""

import unittest
import snouter
import pathlib
from . import data_directory


class TestUtilFunctions(unittest.TestCase):

    def test_execute_command_1(self):
        result = snouter.execute_command('echo test', check=False)
        self.assertEqual(result.returncode, 0)
        self.assertEqual(result.stdout, 'test\n')
        self.assertEqual(result.stderr, '')

    def test_execute_command_2(self):
        result = snouter.execute_command('echo test 1>&2', check=False)
        self.assertEqual(result.returncode, 0)
        self.assertEqual(result.stdout, '')
        self.assertEqual(result.stderr, 'test\n')

    def test_execute_command_3(self):
        result = snouter.execute_command('this_is_not_a_valid_command', check=False)
        self.assertNotEqual(result.returncode, 0)

    def test_default_thread_count(self):
        """
        The number of threads returned by default_thread_count depends on the machine's CPU, but it
        should never be more than 8.
        """
        thread_count = snouter.default_thread_count()
        self.assertTrue(1 <= thread_count <= 8)

    def test_bulk_mask_parse(self):
        # Tests a standard input file from data/bulk.tsv
        # These mask inputs can be used for future mask tests
        bulk_fn = data_directory / 'bulk.tsv'

        # obj.threads,exclude,read_type,assembly_fp,read_fps1,read_fps2,read_fps3
        result = snouter.bulk_mask_parse(bulk_fn, 8)

        #illumina	as1.fasta	r1_1.fastq.gz	r1_2.fastq.gz
        mask = result[0]
        self.assertEqual(mask.threads,8)
        self.assertEqual(mask.exclude,None)
        self.assertEqual(mask.read_type,'illumina')
        self.assertEqual(mask.assembly_fp,pathlib.Path('as1.fasta'))
        self.assertEqual(mask.read_fps[0],pathlib.Path('r1_1.fastq.gz'))
        self.assertEqual(mask.read_fps[1],pathlib.Path('r1_2.fastq.gz'))

        #long	as2.fasta	r2.fastq.gz
        mask = result[1]
        self.assertEqual(mask.threads,8)
        self.assertEqual(mask.exclude,None)
        self.assertEqual(mask.read_type,'long')
        self.assertEqual(mask.assembly_fp,pathlib.Path('as2.fasta'))
        self.assertEqual(mask.read_fps[0],pathlib.Path('r2.fastq.gz'))

        #7	illumina	as3.fasta	r3_1.fastq.gz	r3_2.fastq.gz
        mask = result[2]
        self.assertEqual(mask.threads,8)
        self.assertEqual(mask.exclude,7)
        self.assertEqual(mask.read_type,'illumina')
        self.assertEqual(mask.assembly_fp,pathlib.Path('as3.fasta'))
        self.assertEqual(mask.read_fps[0],pathlib.Path('r3_1.fastq.gz'))
        self.assertEqual(mask.read_fps[1],pathlib.Path('r3_2.fastq.gz'))



