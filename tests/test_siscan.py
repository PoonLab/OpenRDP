import os
import unittest
import numpy as np
import random

from openrdp.common import TripletGenerator, read_fasta, merge_breakpoints
from openrdp.siscan import Siscan


class TestSiscan(unittest.TestCase):
    def setUp(self):
        # Set up test example
        test_settings = {'max_pvalue': '1.0', 'win_size': '6', 'step_size': '2', 'strip_gaps': 'False',
                         'pvalue_perm_num': '1100', 'scan_perm_num': '100', 'random_seed': '3'}
        short_seq_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'short.fasta')
        with open(short_seq_path) as small_test:
            names, test_seqs = read_fasta(small_test)
            self.short_align = np.array(list(map(list, test_seqs)))
            self.test_short = Siscan(self.short_align, names, settings=test_settings)

        self.short_triplets = [trp for trp in TripletGenerator(self.short_align, names)]

        # Set up test example 2
        test_settings = {'max_pvalue': '0.8', 'win_size': '50', 'step_size': '20', 'strip_gaps': 'True',
                         'pvalue_perm_num': '1100', 'scan_perm_num': '100', 'random_seed': '3'}
        long_seq_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'long.fasta')
        with open(long_seq_path) as test:
            names, test_seqs = read_fasta(test)
            self.long_align = np.array(list(map(list, test_seqs)))
            self.test_long = Siscan(self.long_align, names, settings=test_settings)

        self.long_triplets = [trp for trp in TripletGenerator(self.long_align, names)]

        # Set up HIV CRF07 test case
        test_settings = {'max_pvalue': '0.05', 'win_size': '200', 'step_size': '20', 'strip_gaps': 'True',
                         'pvalue_perm_num': '1100', 'scan_perm_num': '100', 'random_seed': '3'}
        hiv_seq_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'CRF_07_test.fasta')
        with open(hiv_seq_path) as hiv_test:
            names, crf07_seqs = read_fasta(hiv_test)
            self.hiv_align = np.array(list(map(list, crf07_seqs)))
            self.test_hiv = Siscan(self.hiv_align, names, settings=test_settings)

        self.hiv_triplets = [trp for trp in TripletGenerator(self.hiv_align, names)]

        self.pattern_aln = [np.array(['A','G','C','T','A','C','T','G','A','C','A','C','T','G','C','A','T','G','A','C','G','T','A','C','G','T','A','C','G','T']),
                            np.array(['C','A','C','T','T','G','G','T','A','C','C','T','A','C','C','T','G','A','G','T','G','T','A','C','C','A','T','G','G','T']),
                            np.array(['G','T','G','A','A','G','G','A','A','C','G','T','T','G','G','C','T','G','T','G','G','T','G','T','G','T','G','G','G','T'])]
        self.siscan = Siscan(align = self.pattern_aln,scan_perm_num=10, win_size = 1)

    def test_set_and_validate_options(self):
        self.assertEqual(6, self.test_short.win_size)
        self.assertEqual(False, self.test_short.strip_gaps)
        self.assertEqual(2, self.test_short.step_size)
        self.assertEqual(1100, self.test_short.pvalue_perm_num)
        self.assertEqual(100, self.test_short.scan_perm_num)
        self.assertEqual(3, self.test_short.random_seed)

        self.assertEqual(50, self.test_long.win_size)
        self.assertEqual(True, self.test_long.strip_gaps)
        self.assertEqual(20, self.test_long.step_size)
        self.assertEqual(1100, self.test_long.pvalue_perm_num)
        self.assertEqual(100, self.test_long.scan_perm_num)
        self.assertEqual(3, self.test_long.random_seed)

        self.assertEqual(200, self.test_hiv.win_size)
        self.assertEqual(True, self.test_hiv.strip_gaps)
        self.assertEqual(20, self.test_hiv.step_size)
        self.assertEqual(1100, self.test_hiv.pvalue_perm_num)
        self.assertEqual(100, self.test_hiv.scan_perm_num)
        self.assertEqual(3, self.test_hiv.random_seed)

    def test_all_patterns_nucleotides(self):
        expected = [2, 2, 2, 4, 2, 1, 1, 1, 1, 1, 4, 2, 4, 1, 2]

        d = np.array(['T','C','T','G','C','A','A','G','C','T','A','C','T','G','G','C','G','C','A','C','A','C','A','C','G','T','T','G','G','T'])
        counts = self.siscan.count_patterns(self.pattern_aln + [d])
        self.assertEqual(counts, expected)


    def test_sum_pattern(self):
        expected = [6, 6, 6, 8, 8, 8, 8, 8, 8]

        counts = [2] * 15
        self.assertEqual(self.siscan.sum_pattern_counts(counts),expected)
        
    def test_adjust_nt_f(self):
        expected = 3

        a = np.array(['A','A','A','C','C','C','T','T'])
        b = np.array(['C','C','C','C','A','A','T','T'])
        diff = np.array(['C','C','C','A','A','A','T','T'])

        start, end = 0, 7
        out = self.siscan.adjust_nt_f(a, b, diff, start, end)

        self.assertEqual(out, expected)

    def test_adjust_nt_r(self):
        expected = 3

        a = np.array(['A','A','A','C','C','C','T','T'])
        b = np.array(['C','C','C','C','A','A','T','T'])
        diff = np.array(['C','C','C','A','A','A','T','T'])

        start, end = 7, 0
        out = self.siscan.adjust_nt_r(a, b, diff, start, end)
        
        self.assertEqual(out, expected)

    def test_shuffle_pattern(self):
        """
        test shuffle if we are looking for a specific pattern
        """
        np.random.seed(1)
        expected = -0.06917
        parent1, parent2, recomb = self.siscan.align
        out = np.array(['T','C','T','G','C','A','A','G','C','T','A','C','T','G','G','C','G','C','A','C','A','C','A','C','G','T','T','G','G','T'])
        ps_ind = [True, 1]

        # this will be the null distribution we care about
        # [1, 4, 2, 4, 1, 1, 0, 3, 4, 1]
        zscore = self.siscan.shuffle(parent1, parent2, recomb, out, ps_ind)
        zscore = round(zscore, 5)

        self.assertEqual(zscore, expected)

    def test_shuffle_sum(self):
        """
        test shuffle if we are looking for a specific pattern
        """
        np.random.seed(1)
        expected = 1.90516
        parent1, parent2, recomb = self.siscan.align
        out = np.array(['T','C','T','G','C','A','A','G','C','T','A','C','T','G','G','C','G','C','A','C','A','C','A','C','G','T','T','G','G','T'])
        ps_ind = [False, 2]

        # this will be the null distribution we care about
        # [4, 5, 6, 7, 2, 3, 3, 3, 4, 5]
        zscore = self.siscan.shuffle(parent1, parent2, recomb, out, ps_ind)
        zscore = round(zscore, 5)

        self.assertEqual(zscore, expected)

    def get_ps(self):
        expected = (False, 0) # seq_1_2 ind 2

        self.assertEqual(self.siscan.get_ps(0,0,2), expected)

    def test_find_interval(self):
        # Set up 2D arrays (list of lists of 3 values)
        maj = [
            [100, 100, 100],  # index 0
            [100, 100, 100],  # index 1
            [100, 100, 100],  # index 2
            [100, 100, 100],  # index 3
            [100, 100, 100],  # index 4
            [100, 100, 100],  # index 5
            [100, 100, 100],  # index 6
            [100, 100, 100],  # index 7
            [100, 0, 100],    # index 8
            [100, 100, 100]   # index 9
        ]
        close = [
            [1, 1, 1],  
            [1, 1, 1],
            [1, 1, 1],
            [1, 1, 1],
            [1, 1, 1],
            [1, 1, 1],
            [1, 1, 1],
            [1, 1, 1],
            [1, 1, 1],
            [1, 1, 1]
        ]
        expected = (1,7)
        self.assertEqual(self.siscan.find_interval(1, maj, close, 1), expected)

    def test_find_signal(self):
        close = np.array([
            [5] * 14,
            [5] * 14,
            [5] * 14,
            [5] * 14,
            [5] * 14,
            [5] * 14,
            [5] * 14,
            [5] * 14,
            [5] * 14,
            [5] * 14
        ])  # shape (10, 14)

        maj1 = np.array([
            [1] * 14,
            [1] * 14,
            [1] * 14,
            [7] * 14,
            [7] * 14,
            [7] * 14,
            [7] * 14,
            [7] * 14,
            [7] * 14,
            [7] * 14
        ])  # shape (10, 14)

        maj2 = np.array([
            [1] * 14,
            [7] * 14,
            [7] * 14,
            [7] * 14,
            [7] * 14,
            [7] * 14,
            [7] * 14,
            [7] * 14,
            [7] * 14,
            [1] * 14
        ])  # shape (10, 14)

        ind = 2
        expected1 = [
            (3, 9, 2, (True, 2)),
            (3, 9, 2, (True, 8)),
            (3, 9, 2, (False, 4)),
            (3, 9, 2, (False, 5))
        ]

        expected2 = [
            (1, 8, 2, (True, 4)),
            (1, 8, 2, (True, 9)),
            (1, 8, 2, (False, 2)),
            (1, 8, 2, (False, 6))
        ]

        self.assertEqual(self.siscan.find_signal(close, maj1,maj2,ind), (expected1, expected2))

    # These tests don't really make sense with the refactoring

    # def test_execute_short(self):
    #     expected = [('A', ('B', 'C'), 2, 11, 0.7787397893226587),
    #                 ('A', ('B', 'D'), 3, 11, 0.7524567011843551),
    #                 ('A', ('B', 'E'), 2, 11, 0.7787397893226587),
    #                 ('A', ('C', 'D'), 2, 11, 0.7685248858128534),
    #                 ('A', ('D', 'E'), 2, 11, 0.7260095477693602),
    #                 ('B', ('C', 'E'), 3, 11, 0.7750581068141527),
    #                 ('B', ('D', 'E'), 3, 11, 0.745172605746384),
    #                 ('C', ('A', 'E'), 2, 11, 0.804804577675132),
    #                 ('C', ('B', 'D'), 2, 11, 0.7838806032521947),
    #                 ('C', ('D', 'E'), 2, 11, 0.7465417522802897)]

    #     for trp in self.short_triplets:
    #         self.test_short.execute(trp)
    #     result = merge_breakpoints(self.test_short.raw_results, self.test_short.max_pvalues)
    #     self.assertEqual(expected, result)

    # def test_execute_long(self):
    #     expected = [('Test1 ', ('Test2', 'Test3'), 2, 55, 0.7521364874425968),
    #                 ('Test1 ', ('Test2', 'Test4'), 2, 55, 0.7669294010627822),
    #                 ('Test1 ', ('Test3', 'Test4'), 2, 55, 0.7651446184495413),
    #                 ('Test2', ('Test3', 'Test4'), 2, 55, 0.7651446184495413)]

    #     for trp in self.long_triplets:
    #         self.test_long.execute(trp)
    #     result = merge_breakpoints(self.test_long.raw_results, self.test_long.max_pvalues)
    #     self.assertEqual(expected, result)

    # def test_execute_hiv(self):
    #     expected = []   # Breakpoints have p_values that are too large (above threshold)
    #     for trp in self.hiv_triplets:
    #         self.test_hiv.execute(trp)
    #     result = merge_breakpoints(self.test_hiv.raw_results, self.test_hiv.max_pvalues)
    #     self.assertEqual(expected, result)


if __name__ == '__main__':
    unittest.main()
