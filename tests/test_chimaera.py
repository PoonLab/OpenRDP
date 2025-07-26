import configparser
import os
import unittest

import numpy as np

from openrdp.chimaera import Chimaera
from openrdp.common import TripletGenerator, read_fasta, merge_breakpoints


class TestChimaera(unittest.TestCase):
    def setUp(self):
        # Set up test example
        test_settings = {'max_pvalue': '1.0', 'win_size': '5', 'strip_gaps': 'False',
                         'fixed_win_size': 'True', 'num_var_sites': '5', 'frac_var_sites': '0.1'}

        short_seq_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'short.fasta')
        with open(short_seq_path) as small_test:
            names, test_seqs = read_fasta(small_test)
            self.short_align = np.array(list(map(list, test_seqs)))
            self.test_short = Chimaera(self.short_align, settings=test_settings)

        self.short_triplets = [trp for trp in TripletGenerator(self.short_align, names)]

        # Set up test example 2
        test_settings = {'max_pvalue': '0.05', 'win_size': '40', 'strip_gaps': 'False',
                         'fixed_win_size': 'True', 'num_var_sites': '70', 'frac_var_sites': '0.1'}

        long_seq_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'long.fasta')
        with open(long_seq_path) as test:
            names, test_seqs = read_fasta(test)
            self.long_align = np.array(list(map(list, test_seqs)))
            self.test_long = Chimaera(self.long_align, settings=test_settings)

        self.long_triplets = [trp for trp in TripletGenerator(self.long_align, names)]

        # Set up HIV CRF07 test case
        test_settings = {'max_pvalue': '0.05', 'win_size': '50', 'strip_gaps': 'False',
                         'fixed_win_size': 'True', 'num_var_sites': '60', 'frac_var_sites': '0.1'}

        hiv_seq_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'CRF_07_test.fasta')
        with open(hiv_seq_path) as hiv_test:
            names, crf07_seqs = read_fasta(hiv_test)
            self.hiv_align = np.array(list(map(list, crf07_seqs)))
            self.test_hiv = Chimaera(self.hiv_align, settings=test_settings)

        self.hiv_triplets = [trp for trp in TripletGenerator(self.hiv_align, names)]

    def test_set_and_validate_options(self):
        self.assertEqual(5, self.test_short.win_size)
        self.assertEqual(False, self.test_short.strip_gaps)
        self.assertEqual(True, self.test_short.fixed_win_size)
        self.assertEqual(5, self.test_short.num_var_sites)
        self.assertEqual(0.1, self.test_short.frac_var_sites)

        self.assertEqual(40, self.test_long.win_size)
        self.assertEqual(False, self.test_long.strip_gaps)
        self.assertEqual(True, self.test_long.fixed_win_size)
        self.assertEqual(70, self.test_long.num_var_sites)
        self.assertEqual(0.1, self.test_long.frac_var_sites)

        self.assertEqual(False, self.test_hiv.strip_gaps)
        self.assertEqual(True, self.test_hiv.fixed_win_size)
        self.assertEqual(60, self.test_hiv.num_var_sites)
        self.assertEqual(0.1, self.test_hiv.frac_var_sites)

    def test_refine_nt(self):
        seq = np.array(list('00010000100111110001100110')).astype(int)
        start, end = 5, 18

        result = Chimaera.refine_nt(seq, start, end)
        expected = [3,19]

        self.assertEqual(expected, result)

    def test_compress_recombinant(self):
        new_aln = np.array([['T', 'G', 'T', 'G'],
                            ['T', 'T', 'T', 'T'],
                            ['G', 'G', 'G', 'G']])

        expected = [0, 1, 0, 1]
        result = Chimaera.compress_triplet_aln(new_aln)
        self.assertEqual(expected, result)

    def test_get_window_positions_short(self):
        # Test short
        exp_reg_left = [0, 0]
        exp_reg_right = [1, 0]

        seq = [0, 0, 1, 0, 0, 1, 0, 1]
        res_reg_left, res_reg_right = Chimaera.get_window_positions(seq, 0, 4)

        self.assertEqual(exp_reg_left, res_reg_left)
        self.assertEqual(exp_reg_right, res_reg_right)

    def test_get_window_positions_long(self):
        # Test long
        exp_reg_left = [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        exp_reg_right = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0, 0, 0, 0]

        seq = [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

        res_reg_left, res_reg_right = Chimaera.get_window_positions(seq, 0, 100)

        self.assertEqual(exp_reg_left, res_reg_left)
        self.assertEqual(exp_reg_right, res_reg_right)

    def test_get_window_positions_hiv(self):
        # Test HIV (no gaps stripped, first window position)
        exp_reg_left = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0]
        exp_reg_right = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0]

        seq = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
               1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
               1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
               1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
               1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0,
               1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1,
               1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1,
               1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1,
               1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1,
               1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0,
               0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1,
               1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1,
               1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0,
               0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1,
               1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1,
               0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0,
               1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1,
               0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0,
               1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1,
               1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1,
               1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1,
               1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1,
               1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1,
               0, 0, 0, 0, 0, 0]

        res_reg_left, res_reg_right = Chimaera.get_window_positions(seq, 0, 200)

        self.assertEqual(exp_reg_left, res_reg_left)
        self.assertEqual(exp_reg_right, res_reg_right)

        # Testing internal window position
        exp_reg_left = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1, 1, 1, 1]
        exp_reg_right = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0,
                         0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1,
                         1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1,
                         1, 1, 1, 1]
        res_reg_left, res_reg_right = Chimaera.get_window_positions(seq, 500, 200)
        self.assertEqual(exp_reg_left, res_reg_left)
        self.assertEqual(exp_reg_right, res_reg_right)

    def test_compute_contingency_table_short(self):
        reg_left = [0, 1]
        reg_right = [0, 0, 1]
        half_win_size = 5

        expected = [[1, 4, 5], [1, 4, 5], [2, 8, 10]]
        result = Chimaera.compute_contingency_table(reg_left, reg_right, half_win_size)
        self.assertEqual(expected, result)

        reg_left = [0, 1]
        reg_right = [1, 1, 1]
        expected = [[3, 2, 5], [1, 4, 5], [4, 6, 10]]
        result = Chimaera.compute_contingency_table(reg_left, reg_right, half_win_size)
        self.assertEqual(expected, result)

    def test_compute_contingency_table_long(self):
        reg_left = [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        reg_right = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                     0, 0, 0, 0, 0, 0, 0, 0]
        half_win_size = 50

        expected = [[4, 46, 50], [2, 48, 50], [6, 94, 100]]
        result = Chimaera.compute_contingency_table(reg_left, reg_right, half_win_size)
        self.assertEqual(expected, result)

    def test_compute_contingency_table_hiv(self):
        reg_left = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1]
        reg_right = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0,
                     0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1,
                     1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1,
                     1, 1, 1, 1]
        half_win_size = 100

        expected = [[71, 29, 100], [99, 1, 100], [170, 30, 200]]
        result = Chimaera.compute_contingency_table(reg_left, reg_right, half_win_size)
        self.assertEqual(expected, result)

    # def test_execute_short(self):
    #     expected = [('A', ('B', 'C'), 0, 11, 0.8556951983876534),
    #                 ('A', ('B', 'D'), 1, 9, 0.8556951983876534),
    #                 ('A', ('B', 'E'), 3, 11, 0.8556951983876534),
    #                 ('A', ('C', 'D'), 1, 9, 0.8556951983876534),
    #                 ('A', ('C', 'E'), 4, 12, 0.8556951983876534),
    #                 ('A', ('D', 'E'), 1, 11, 0.8556951983876534),
    #                 ('B', ('C', 'D'), 1, 15, 0.8556951983876534),
    #                 ('B', ('D', 'E'), 1, 10, 0.8556951983876534),
    #                 ('C', ('B', 'D'), 3, 11, 0.8556951983876534),
    #                 ('E', ('A', 'D'), 3, 11, 0.8556951983876534)]
    #     for trp in self.short_triplets:
    #         self.test_short.execute(trp)
    #     result = merge_breakpoints(self.test_short.raw_results)
    #     self.assertEqual(expected, result)

    # def test_execute_long(self):
    #     expected = [('Test1 ', ('Test2', 'Test4'), 99, 142, 0.01174095674176845),
    #                 ('Test1 ', ('Test3', 'Test4'), 192, 235, 0.02047438504938101),
    #                 ('Test2', ('Test1 ', 'Test3'), 243, 286, 0.0018132288986577026),
    #                 ('Test2', ('Test3', 'Test4'), 176, 219, 0.0019834358538684586)]
    #     for trp in self.long_triplets:
    #         self.test_long.execute(trp)
    #     result = merge_breakpoints(self.test_long.raw_results)
    #     self.assertEqual(expected, result)

    # def test_execute_hiv(self):
    #     self.maxDiff = None
    #     expected = [('B', ('07_BC', 'C'), 150, 210, 7.72722639263188e-08),
    #                 ('B', ('07_BC', 'C'), 326, 384, 0.0001118857269439367),
    #                 ('B', ('07_BC', 'C'), 403, 520, 2.6235181749165126e-07),
    #                 ('B', ('07_BC', 'C'), 978, 1031, 0.010181534711170576),
    #                 ('B', ('07_BC', 'C'), 1760, 1813, 0.04154873114149176),
    #                 ('B', ('07_BC', 'C'), 1985, 2038, 0.01735126523666451),
    #                 ('B', ('07_BC', 'C'), 2245, 2298, 0.036092176952611674),
    #                 ('B', ('07_BC', 'C'), 3071, 3124, 3.221473917412952e-05),
    #                 ('B', ('07_BC', 'C'), 3673, 3778, 0.04499874628000287),
    #                 ('B', ('07_BC', 'C'), 3959, 4012, 0.008660655610568911),
    #                 ('B', ('07_BC', 'C'), 4264, 4317, 0.021266969405618293),
    #                 ('B', ('07_BC', 'C'), 4838, 4891, 0.026839498926228344),
    #                 ('B', ('07_BC', 'C'), 5111, 5238, 0.04926030700200055),
    #                 ('B', ('07_BC', 'C'), 5909, 5964, 0.0007511156049633134),
    #                 ('B', ('07_BC', 'C'), 6451, 6504, 0.04926030700200055),
    #                 ('C', ('07_BC', 'B'), 611, 675, 0.0014988494708724491)]
    #     for trp in self.hiv_triplets:
    #         self.test_hiv.execute(trp)
    #     result = merge_breakpoints(self.test_hiv.raw_results)

    #     for index, exp_tuple in enumerate(expected):
    #         self.assertEqual(exp_tuple[0], result[index][0])
    #         self.assertEqual(exp_tuple[1], result[index][1])
    #         self.assertEqual(exp_tuple[2], result[index][2])
    #         self.assertEqual(exp_tuple[3], result[index][3])
    #         self.assertAlmostEqual(exp_tuple[4], result[index][4])


if __name__ == '__main__':
    unittest.main()
