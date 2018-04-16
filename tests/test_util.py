import unittest


import pairwise_snp_counter as psc


class TestUtilFunctions(unittest.TestCase):

    def test_execute_command(self):
        result = psc.execute_command('echo test', check=False)
        print(result.returncode)
        print(result.stdout)
        print(result.stderr)
        self.assertEqual(result.returncode, 0)
        self.assertEqual(result.stdout, 'test\n')
        self.assertEqual(result.stderr, '')
