import unittest


import pairwise_snp_counter as psc


class TestUtilFunctions(unittest.TestCase):

    def test_execute_command(self):
        result = psc.execute_command('echo test', check=False)
        self.assertEqual(result.returncode, 0)
        self.assertEqual(result.stdout, 'test\n')
        self.assertEqual(result.stderr, '')

    def test_default_thread_count(self):
        """
        The number of threads returned by default_thread_count depends on the machine's CPU, but it
        should never be more than 8.
        """
        thread_count = psc.default_thread_count()
        self.assertTrue(1 <= thread_count <= 8)
