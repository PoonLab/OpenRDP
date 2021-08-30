import configparser
import os
import unittest

import numpy as np

from openrdp import read_fasta
from scripts.chimaera import Chimaera
from scripts.common import generate_triplets, Triplet


class TestChimaera(unittest.TestCase):

    def setUp(self):
        # Set up test example
        config = configparser.ConfigParser()
        short_cfg_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_short.ini')
        config.read(short_cfg_path)
        test_settings = dict(config.items('Chimaera'))

        short_seq_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'short.fasta')
        with open(short_seq_path) as small_test:
            names, test_seqs = read_fasta(small_test)
            self.short_align = np.array(list(map(list, test_seqs)))
            self.test_short = Chimaera(self.short_align, names, settings=test_settings)

        self.short_triplets = []
        for trp in generate_triplets(self.short_align):
            self.short_triplets.append(Triplet(self.short_align, names, trp))

        # Set up test example 2
        config = configparser.ConfigParser()
        long_cfg_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_long.ini')
        config.read(long_cfg_path)
        test_settings = dict(config.items('Chimaera'))

        long_seq_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'long.fasta')
        with open(long_seq_path) as test:
            names, test_seqs = read_fasta(test)
            self.long_align = np.array(list(map(list, test_seqs)))
            self.test_long = Chimaera(self.long_align, names, settings=test_settings)

        self.long_triplets = []
        for trp in generate_triplets(self.long_align):
            self.long_triplets.append(Triplet(self.long_align, names, trp))

        # Set up HIV CRF07 test case
        config = configparser.ConfigParser()
        hiv_cfg_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'default.ini')
        config.read(hiv_cfg_path)
        settings = dict(config.items('Chimaera'))

        hiv_seq_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'CRF_07_test.fasta')
        with open(hiv_seq_path) as hiv_test:
            names, crf07_seqs = read_fasta(hiv_test)
            self.hiv_align = np.array(list(map(list, crf07_seqs)))
            self.test_hiv = Chimaera(self.hiv_align, names, settings=settings)

        self.hiv_triplets = []
        for trp in generate_triplets(self.hiv_align):
            self.hiv_triplets.append(Triplet(self.hiv_align, names, trp))

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

    def test_execute_short(self):
        expected = {('A', 'B', 'C'): [(0, 8, 0.8556951983876534)],
                    ('A', 'B', 'D'): [],
                    ('A', 'B', 'E'): [(3, 11, 0.8556951983876534)],
                    ('A', 'C', 'D'): [(1, 9, 0.8556951983876534)],
                    ('A', 'C', 'E'): [],
                    ('A', 'D', 'E'): [(3, 11, 0.8556951983876534)],
                    ('B', 'C', 'A'): [(0, 8, 0.8556951983876534)],
                    ('B', 'C', 'D'): [(3, 11, 0.8556951983876534)],
                    ('B', 'C', 'E'): [],
                    ('B', 'D', 'A'): [],
                    ('B', 'D', 'E'): [(1, 9, 0.8556951983876534)],
                    ('B', 'E', 'A'): [(3, 11, 0.8556951983876534)],
                    ('C', 'B', 'A'): [(0, 8, 0.8556951983876534)],
                    ('C', 'D', 'A'): [(1, 9, 0.8556951983876534)],
                    ('C', 'D', 'B'): [(3, 11, 0.8556951983876534)],
                    ('C', 'D', 'E'): [],
                    ('C', 'E', 'A'): [],
                    ('C', 'E', 'B'): [],
                    ('D', 'B', 'A'): [],
                    ('D', 'C', 'A'): [(1, 9, 0.8556951983876534)],
                    ('D', 'C', 'B'): [(3, 11, 0.8556951983876534)],
                    ('D', 'E', 'A'): [(3, 11, 0.8556951983876534)],
                    ('D', 'E', 'B'): [(1, 9, 0.8556951983876534)],
                    ('D', 'E', 'C'): [],
                    ('E', 'B', 'A'): [(3, 11, 0.8556951983876534)],
                    ('E', 'C', 'A'): [],
                    ('E', 'C', 'B'): [],
                    ('E', 'D', 'A'): [(3, 11, 0.8556951983876534)],
                    ('E', 'D', 'B'): [(1, 9, 0.8556951983876534)],
                    ('E', 'D', 'C'): []}
        result = self.test_short.execute(self.short_triplets, quiet=True)
        self.assertEqual(expected, result)

    def test_execute_long(self):
        expected = {('Test1 ', 'Test2', 'Test3'): [],
                    ('Test1 ', 'Test2', 'Test4'): [],
                    ('Test1 ', 'Test3', 'Test4'): [],
                    ('Test2', 'Test3', 'Test1 '): [],
                    ('Test2', 'Test3', 'Test4'): [(176, 219, 1.0)],
                    ('Test2', 'Test4', 'Test1 '): [],
                    ('Test3', 'Test2', 'Test1 '): [],
                    ('Test3', 'Test4', 'Test1 '): [],
                    ('Test3', 'Test4', 'Test2'): [(176, 219, 1.0)],
                    ('Test4', 'Test2', 'Test1 '): [],
                    ('Test4', 'Test3', 'Test1 '): [],
                    ('Test4', 'Test3', 'Test2'): [(176, 219, 1.0)]}
        result = self.test_long.execute(self.long_triplets, quiet=True)
        self.assertEqual(expected, result)

    def test_execute_hiv(self):
        self.maxDiff = None
        expected = {('B', 'C', '07_BC'): [(457, 520, 6.081144750529128e-08),
                                          (611, 667, 0.0014988494708724491),
                                          (2998, 3051, 0.003441468681534075),
                                          (3071, 3124, 3.221473917412952e-05),
                                          (3133, 3186, 0.006267877878701045),
                                          (3332, 3385, 1.0),
                                          (3723, 3776, 0.026961195032282365),
                                          (3936, 3989, 1.0),
                                          (4748, 4802, 1.0),
                                          (4838, 4891, 0.026839498926228344),
                                          (5111, 5164, 0.04926030700200058),
                                          (5185, 5238, 0.010181534711170576)],
                    ('C', '07_BC', 'B'): [(457, 520, 6.081144750529128e-08),
                                          (611, 667, 0.0014988494708724491),
                                          (2998, 3051, 0.003441468681534075),
                                          (3071, 3124, 3.221473917412952e-05),
                                          (3133, 3186, 0.006267877878701045),
                                          (3332, 3385, 1.0),
                                          (3723, 3776, 0.026961195032282365),
                                          (3936, 3989, 1.0),
                                          (4748, 4802, 1.0),
                                          (4838, 4891, 0.026839498926228344),
                                          (5111, 5164, 0.04926030700200058),
                                          (5185, 5238, 0.010181534711170576)],
                    ('07_BC', 'C', 'B'): [(457, 520, 6.081144750529128e-08),
                                          (611, 667, 0.0014988494708724491),
                                          (2998, 3051, 0.003441468681534075),
                                          (3071, 3124, 3.221473917412952e-05),
                                          (3133, 3186, 0.006267877878701045),
                                          (3332, 3385, 1.0),
                                          (3723, 3776, 0.026961195032282365),
                                          (3936, 3989, 1.0),
                                          (4748, 4802, 1.0),
                                          (4838, 4891, 0.026839498926228344),
                                          (5111, 5164, 0.04926030700200058),
                                          (5185, 5238, 0.010181534711170576)]}

        result = self.test_hiv.execute(self.hiv_triplets, quiet=True)
        self.assertEqual(expected, result)


if __name__ == '__main__':
    unittest.main()
