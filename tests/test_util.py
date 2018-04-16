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


class TestUtilFunctions(unittest.TestCase):

    def test_execute_command_1(self):
        result = psc.execute_command('echo test', check=False)
        self.assertEqual(result.returncode, 0)
        self.assertEqual(result.stdout, 'test\n')
        self.assertEqual(result.stderr, '')

    def test_execute_command_2(self):
        result = psc.execute_command('echo test 1>&2', check=False)
        self.assertEqual(result.returncode, 0)
        self.assertEqual(result.stdout, '')
        self.assertEqual(result.stderr, 'test\n')

    def test_execute_command_3(self):
        result = psc.execute_command('this_is_not_a_valid_command', check=False)
        self.assertNotEqual(result.returncode, 0)

    def test_default_thread_count(self):
        """
        The number of threads returned by default_thread_count depends on the machine's CPU, but it
        should never be more than 8.
        """
        thread_count = psc.default_thread_count()
        self.assertTrue(1 <= thread_count <= 8)
