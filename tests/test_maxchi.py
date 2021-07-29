import unittest
from scripts.maxchi import MaxChi
from scripts.main import read_fasta
from scripts.common import generate_triplets, Triplet
import configparser
import numpy as np


class TestMaxChi(unittest.TestCase):

    def setUp(self):
        # Set up test example
        config = configparser.ConfigParser()
        config.read('test_short.ini')
        test_settings = dict(config.items('MaxChi'))

        with open('short.fasta') as small_test:
            names, test_seqs = read_fasta(small_test)
            self.short_align = np.array(list(map(list, test_seqs)))
            self.test_short = MaxChi(self.short_align, names, settings=test_settings)

        self.short_triplets = []
        for trp in generate_triplets(self.short_align):
            self.short_triplets.append(Triplet(self.short_align, names, trp))

        # Set up test example 2
        config = configparser.ConfigParser()
        config.read('test_long.ini')
        test_settings = dict(config.items('MaxChi'))

        with open('long.fasta') as test:
            names, test_seqs = read_fasta(test)
            self.long_align = np.array(list(map(list, test_seqs)))
            self.test_long = MaxChi(self.long_align, names, settings=test_settings)

        self.long_triplets = []
        for trp in generate_triplets(self.long_align):
            self.long_triplets.append(Triplet(self.long_align, names, trp))

        # Set up HIV CRF07 test case
        config = configparser.ConfigParser()
        config.read('default.ini')
        settings = dict(config.items('MaxChi'))

        with open('CRF_07_test.fasta') as hiv_test:
            names, crf07_seqs = read_fasta(hiv_test)
            self.hiv_align = np.array(list(map(list, crf07_seqs)))
            self.test_hiv = MaxChi(self.hiv_align, names, settings=settings)

        self.hiv_triplets = []
        for trp in generate_triplets(self.hiv_align):
            self.hiv_triplets.append(Triplet(self.hiv_align, names, trp))

    def test_set_and_validate_options(self):
        self.assertEqual(6, self.test_short.win_size)
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
        self.assertEqual(70, self.test_hiv.num_var_sites)
        self.assertEqual(0.1, self.test_hiv.frac_var_sites)

    def test_get_window_positions(self):
        # Test short
        exp_reg1_right = 'GACGA'
        exp_reg1_left = 'ATGCT'
        exp_reg2_right = 'TGGTG'
        exp_reg2_left = 'AACCT'
        exp_reg1 = 'ATGCTGACGA'
        exp_reg2 = 'AACCTTGGTG'

        seq1 = self.short_triplets[0].sequences[0]
        seq2 = self.short_triplets[0].sequences[1]

        res_reg1_left, res_reg2_left, res_reg1_right, res_reg2_right, res_reg1, res_reg2 = \
            MaxChi.get_window_positions(seq1, seq2, 0, 10)

        self.assertEqual(exp_reg1_right, ''.join(res_reg1_right))
        self.assertEqual(exp_reg1, ''.join(res_reg1))
        self.assertEqual(exp_reg2, ''.join(res_reg2))
        self.assertEqual(exp_reg1_left, ''.join(res_reg1_left))
        self.assertEqual(exp_reg2_right, ''.join(res_reg2_right))
        self.assertEqual(exp_reg2_left, ''.join(res_reg2_left))

        # Test long
        exp_reg1_right = 'GACGA'
        exp_reg1_left = 'ATGCT'
        exp_reg2_right = 'TGGTG'
        exp_reg2_left = 'AACCT'
        exp_reg1 = 'ATGCTGACGA'
        exp_reg2 = 'AACCTTGGTG'

        seq1 = self.long_triplets[1].sequences[1]
        seq2 = self.long_triplets[1].sequences[2]

        res_reg1_left, res_reg2_left, res_reg1_right, res_reg2_right, res_reg1, res_reg2 = \
            MaxChi.get_window_positions(seq1, seq2, 0, 100)

        self.assertEqual(exp_reg1_right, ''.join(res_reg1_right))
        self.assertEqual(exp_reg1, ''.join(res_reg1))
        self.assertEqual(exp_reg2, ''.join(res_reg2))
        self.assertEqual(exp_reg1_left, ''.join(res_reg1_left))
        self.assertEqual(exp_reg2_right, ''.join(res_reg2_right))
        self.assertEqual(exp_reg2_left, ''.join(res_reg2_left))

    def test_execute_short(self):
        expected = {('A', 'B'): [],
                    ('A', 'C'): [],
                    ('A', 'D'): [],
                    ('A', 'E'): [],
                    ('B', 'C'): [],
                    ('B', 'D'): [],
                    ('B', 'E'): [],
                    ('C', 'D'): [],
                    ('C', 'E'): [],
                    ('D', 'E'): []}
        result = self.test_short.execute()
        self.assertEqual(expected, result)

    def test_execute_long(self):
        expected = {('Test1 ', 'Test2'): [(6, 50, 0.026098291541307085)],
                    ('Test1 ', 'Test3'): [],
                    ('Test1 ', 'Test4'): [],
                    ('Test2', 'Test3'): [(2, 47, 0.023868442164574358)],
                    ('Test2', 'Test4'): [],
                    ('Test3', 'Test4'): [(0, 41, 0.01572529975450535)]}
        result = self.test_long.execute()
        self.assertEqual(expected, result)

    def test_execute_hiv(self):
        expected = {('07_BC', 'B'): [(27, 144, 0.035258260485873466),
                                     (140, 252, 0.01869359769155257),
                                     (253, 361, 0.009873332091192726),
                                     (419, 527, 0.005183132922906616),
                                     (562, 671, 0.005183132922906616),
                                     (677, 787, 0.0026997960632601883)],
                    ('07_BC', 'C'): [(35, 145, 0.035258260485873466),
                                     (141, 255, 0.035258260485873466),
                                     (262, 376, 0.01869359769155257),
                                     (460, 567, 0.009873332091192726),
                                     (566, 673, 0.005183132922906616)],
                    ('B', 'C'): [(43, 160, 0.035258260485873466),
                                 (141, 247, 0.01869359769155257)]}
        result = self.test_hiv.execute()
        self.assertEqual(expected, result)
