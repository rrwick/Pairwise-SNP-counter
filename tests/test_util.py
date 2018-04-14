import unittest


import pairwise_snp_counter as psc


class TestUtilFunctions(unittest.TestCase):

    def test_execute_command(self):
        result = psc.execute_command('echo -n test')
        self.assertEqual(result.returncode, 0)
        self.assertEqual(result.stdout, 'test')
        self.assertEqual(result.stderr, '')
