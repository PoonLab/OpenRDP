import unittest
import numpy as np
import os

from openrdp.bootscan import Bootscan
from openrdp.common import TripletGenerator, Triplet, read_fasta
from itertools import combinations


basepath = os.path.dirname(os.path.abspath(__file__))
SHORT_INFILE = os.path.join(basepath, 'short.fasta')
LONG_INFILE = os.path.join(basepath, 'long.fasta')
HIV_INFILE = os.path.join(basepath, 'CRF_07_test.fasta')


class TestBootscan(unittest.TestCase):
    def setUp(self):
        # Set up short test fixtures
        test_settings = {'max_pvalue': '1.0', 'win_size': '6', 'step_size': '1', 'num_replicates': '3',
                         'random_seed': '3', 'cutoff_percentage': '0.7', 'scan': 'distances', 'np': '2'}
        with open(SHORT_INFILE) as short_test:
            names, test_seqs = read_fasta(short_test)
        self.short_align = np.array(list(map(list, test_seqs)))
        self.test_short = Bootscan(self.short_align, settings=test_settings, quiet=True)
        self.short_names = names

        # Set up long test fixtures
        test_settings = {'max_pvalue': '0.05', 'win_size': '50', 'step_size': '5', 'num_replicates': '100',
                         'random_seed': '3', 'cutoff_percentage': '0.7', 'scan': 'distances', 'np': '2'}
        with open(LONG_INFILE) as long_test:
            names, test_seqs = read_fasta(long_test)
        self.long_align = np.array(list(map(list, test_seqs)))
        self.test_long = Bootscan(self.long_align, settings=test_settings, quiet=True)
        self.long_names = names
        
        # Set up HIV test fixtures
        test_settings = {'max_pvalue': '0.05', 'win_size': '200', 'step_size': '20', 'num_replicates': '100',
                         'random_seed': '3', 'cutoff_percentage': '0.7', 'scan': 'distances', 'np': '2',
                         'p_value_calculation': 'binomial', 'model': 'Jukes-Cantor'}
        with open(HIV_INFILE) as hiv_test:
            names, test_seqs = read_fasta(hiv_test)
        self.hiv_align = np.array(list(map(list, test_seqs)))
        self.test_hiv = Bootscan(self.hiv_align, settings=test_settings, quiet=True)
        self.hiv_names = names

    def tearDown(self):
        if os.path.exists(self.test_short.dt_matrix_file):
            os.remove(self.test_short.dt_matrix_file)

        if os.path.exists(self.test_long.dt_matrix_file):
            os.remove(self.test_long.dt_matrix_file)

        if os.path.exists(self.test_hiv.dt_matrix_file):
            os.remove(self.test_hiv.dt_matrix_file)

    def test_set_and_validate_options(self):
        self.assertEqual(6, self.test_short.win_size)
        self.assertEqual(1, self.test_short.step_size)
        self.assertEqual(3, self.test_short.num_replicates)
        self.assertEqual(3, self.test_short.random_seed)
        self.assertEqual(0.7, self.test_short.cutoff)

        self.assertEqual(50, self.test_long.win_size)
        self.assertEqual(5, self.test_long.step_size)
        self.assertEqual(100, self.test_long.num_replicates)
        self.assertEqual(3, self.test_long.random_seed)
        self.assertEqual(0.7, self.test_long.cutoff)

        self.assertEqual(200, self.test_hiv.win_size)
        self.assertEqual(20, self.test_hiv.step_size)
        self.assertEqual(100, self.test_hiv.num_replicates)
        self.assertEqual(3, self.test_hiv.random_seed)
        self.assertEqual(0.7, self.test_hiv.cutoff)

    def test_execute_short(self):
        total_num_trps = sum(1 for _ in combinations(range(self.short_align.shape[0]), 3))
        expected = []   # No breakpoints found
        self.test_short.execute_all(total_combinations=total_num_trps, seq_names=self.short_names)
        result = self.test_short.merge_breakpoints()
        self.assertEqual(expected, result)

    def test_execute_long(self):
        total_num_trps = sum(1 for _ in combinations(range(self.long_align.shape[0]), 3))
        expected = []   # P-value of breakpoints is outside the threshold
        self.test_long.execute_all(total_combinations=total_num_trps, seq_names=self.long_names)
        result = self.test_long.merge_breakpoints()
        self.assertEqual(expected, result)

    def test_execute_hiv(self):
        expected = [('07_BC', ('B', 'C'), 580, 1260, 1.0158213242512536e-09),
                    ('07_BC', ('B', 'C'), 1280, 1340, 0.024030305239701163),
                    ('07_BC', ('B', 'C'), 2060, 2560, 1.7427448207965451e-06),
                    ('07_BC', ('B', 'C'), 2620, 2900, 2.1518293256796478e-07),
                    ('07_BC', ('B', 'C'), 3000, 3160, 0.003489818124092514),
                    ('07_BC', ('B', 'C'), 3260, 5680, 1.5236266249248996e-40),
                    ('07_BC', ('B', 'C'), 8960, 9140, 0.0010069350767726085),
                    ('B', ('07_BC', 'C'), 100, 560, 3.53955031313494e-19),
                    ('B', ('07_BC', 'C'), 1280, 1340, 0.022012231965328217),
                    ('B', ('07_BC', 'C'), 2060, 2560, 1.6836162768997676e-11),
                    ('B', ('07_BC', 'C'), 3000, 3160, 8.361164579231644e-08),
                    ('B', ('07_BC', 'C'), 8960, 9140, 1.5107519378439202e-05),
                    ('C', ('07_BC', 'B'), 580, 1260, 0.04109525585896175),
                    ('C', ('07_BC', 'B'), 2620, 2900, 0.00022445087998712194),
                    ('C', ('07_BC', 'B'), 3260, 5680, 9.28749065950454e-34)]

        total_num_trps = sum(1 for _ in combinations(range(self.hiv_align.shape[0]), 3))
        self.test_hiv.execute_all(total_combinations=total_num_trps, seq_names=self.hiv_names)
        result = self.test_hiv.merge_breakpoints()
        self.assertEqual(expected, result)


if __name__ == '__main__':
    unittest.main()
