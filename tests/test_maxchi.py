import configparser
import os
import unittest
from io import StringIO
import numpy as np

from openrdp.common import TripletGenerator, read_fasta, merge_breakpoints
from openrdp.maxchi import MaxChi


class TestMaxChi(unittest.TestCase):

    def setUp(self):
        # Set up test example
        test_settings = {'max_pvalue': '1.0', 'win_size': '6', 'strip_gaps': 'False',
                         'fixed_win_size': 'True', 'num_var_sites': '5', 'frac_var_sites': '0.1'}

        short_seq_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'short.fasta')
        with open(short_seq_path) as small_test:
            names, test_seqs = read_fasta(small_test)
            self.short_align = np.array(list(map(list, test_seqs)))
            self.test_short = MaxChi(self.short_align, settings=test_settings)

        self.short_triplets = [trp for trp in TripletGenerator(self.short_align, names)]

        # Set up test example 2
        test_settings = {'max_pvalue': '0.05', 'win_size': '40', 'strip_gaps': 'False',
                         'fixed_win_size': 'True', 'num_var_sites': '70', 'frac_var_sites': '0.1'}

        long_seq_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'long.fasta')
        with open(long_seq_path) as test:
            names, test_seqs = read_fasta(test)
            self.long_align = np.array(list(map(list, test_seqs)))
            self.test_long = MaxChi(self.long_align, settings=test_settings)

        self.long_triplets = [trp for trp in TripletGenerator(self.long_align, names)]

        # Set up HIV CRF07 test case
        test_settings = {'max_pvalue': '0.05', 'win_size': '100', 'strip_gaps': 'False',
                         'fixed_win_size': 'True', 'num_var_sites': '70', 'frac_var_sites': '0.1'}

        hiv_seq_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'CRF_07_test.fasta')
        with open(hiv_seq_path) as hiv_test:
            names, crf07_seqs = read_fasta(hiv_test)
            self.hiv_align = np.array(list(map(list, crf07_seqs)))
            self.test_hiv = MaxChi(self.hiv_align, settings=test_settings)

        self.hiv_triplets = [trp for trp in TripletGenerator(self.hiv_align, names)]

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

    def test_get_window_positions_short(self):
        # Test short
        exp_reg1_right = 'GACGA'
        exp_reg1_left = 'ATGCT'
        exp_reg2_right = 'TGGTG'
        exp_reg2_left = 'AACCT'

        seq1 = self.short_triplets[0].sequences[0]
        seq2 = self.short_triplets[0].sequences[1]

        res_reg1_left, res_reg2_left, res_reg1_right, res_reg2_right = \
            MaxChi.get_window_positions(seq1, seq2, 5, 0, 10)

        self.assertEqual(exp_reg1_right, ''.join(res_reg1_right))
        self.assertEqual(exp_reg1_left, ''.join(res_reg1_left))
        self.assertEqual(exp_reg2_right, ''.join(res_reg2_right))
        self.assertEqual(exp_reg2_left, ''.join(res_reg2_left))

    def test_get_window_positions_long(self):
        # Test long
        exp_reg1_right = 'AAGGCGCGGGAGTGACTTATTTAGAGCCGTCCGCCAGCCAAATCGGGCAT'
        exp_reg1_left = 'AAAAACCTCTACTCGGACGCGCTGCGCGTTTGAAGTCGCCGCGCGCGATC'
        exp_reg2_right = 'AAGGGGCGGGGGTGACTTATCTGGAGCCGTCCGCCAGCCAAATCAGGCAT'
        exp_reg2_left = 'AAAAACATCGACCCGCACCCGCTGCGCGTTTGAAGTCGCCGCACGTGACC'

        seq1 = self.long_triplets[1].sequences[1]
        seq2 = self.long_triplets[1].sequences[2]

        res_reg1_left, res_reg2_left, res_reg1_right, res_reg2_right = \
            MaxChi.get_window_positions(seq1, seq2, 50, 0, 100)

        self.assertEqual(exp_reg1_right, ''.join(res_reg1_right))
        self.assertEqual(exp_reg1_left, ''.join(res_reg1_left))
        self.assertEqual(exp_reg2_right, ''.join(res_reg2_right))
        self.assertEqual(exp_reg2_left, ''.join(res_reg2_left))

    def test_get_window_positions_hiv(self):
        # Test HIV (no gaps stripped)
        exp_reg1_right = '----------------------------------------------------' \
                         '------------------------------------------------'
        exp_reg1_left = '-----------------------------------------------------' \
                        '-----------------------------------------------'
        exp_reg2_right = 'ACCAGGACCAGGGACCAGATTTCCACTGACTTTTGGGTGGTGCTTCAAGCTAG' \
                         'TACCAGTTGACCCAGGGGAAGTAGAAGAGGCCAACGAAGGAGAAGAC'
        exp_reg2_left = 'GAATTCTGGAAGGGTTAATTTACTCTAAGAAAAGGCAAGAGATCCTTGATTTGT' \
                        'GGGTCTATCACACACAAGGCTACTTCCCTGATTGGCACAACTACAC'

        seq1 = self.hiv_triplets[0].sequences[0]
        seq2 = self.hiv_triplets[0].sequences[2]

        res_reg1_left, res_reg2_left, res_reg1_right, res_reg2_right = \
            MaxChi.get_window_positions(seq1, seq2, 100, 0, 200)

        self.assertEqual(exp_reg1_right, ''.join(res_reg1_right))
        self.assertEqual(exp_reg1_left, ''.join(res_reg1_left))
        self.assertEqual(exp_reg2_right, ''.join(res_reg2_right))
        self.assertEqual(exp_reg2_left, ''.join(res_reg2_left))

    def test_compute_contingency_table_short(self):
        reg1_right = np.array(list('GACGA'))
        reg2_right = np.array(list('TGGTG'))
        reg1_left  = np.array(list('ATGCT'))
        reg2_left  = np.array(list('AACCT'))

        expected = [[3, 2, 5],
                    [0, 5, 5],
                    [3, 7, 10]]
        result = MaxChi.compute_contingency_table(reg1_right, reg2_right, reg1_left, reg2_left)
        self.assertEqual(expected, result.tolist())

    def test_compute_contingency_table_long(self):
        reg1_right = np.array(list('AAGGCGCGGGAGTGACTTATTTAGAGCCGTCCGCCAGCCAAATCGGGCAT'))
        reg1_left  = np.array(list('AAAAACCTCTACTCGGACGCGCTGCGCGTTTGAAGTCGCCGCGCGCGATC'))
        reg2_right = np.array(list('AAGGGGCGGGGGTGACTTATCTGGAGCCGTCCGCCAGCCAAATCAGGCAT'))
        reg2_left  = np.array(list('AAAAACATCGACCCGCACCCGCTGCGCGTTTGAAGTCGCCGCACGTGACC'))

        expected = [[42, 8, 50],
                    [45, 5, 50],
                    [87, 13, 100]]
        
        result = MaxChi.compute_contingency_table(reg1_right, reg2_right, reg1_left, reg2_left)
        self.assertEqual(expected, result.tolist())

    def test_compute_contingency_table_hiv(self):
        reg1_right = np.array(list(
            '----------------------------------------------------'
            '------------------------------------------------'
        ))

        reg1_left = np.array(list(
            '-----------------------------------------------------'
            '-----------------------------------------------'
        ))

        reg2_right = np.array(list(
            'ACCAGGACCAGGGACCAGATTTCCACTGACTTTTGGGTGGTGCTTCAAGCTAG'
            'TACCAGTTGACCCAGGGGAAGTAGAAGAGGCCAACGAAGGAGAAGAC'
        ))

        reg2_left = np.array(list(
            'GAATTCTGGAAGGGTTAATTTACTCTAAGAAAAGGCAAGAGATCCTTGATTTGT'
            'GGGTCTATCACACACAAGGCTACTTCCCTGATTGGCACAACTACAC'
        ))

        expected = [[0, 100, 100], [0, 100, 100], [0, 200, 200]]
        result = MaxChi.compute_contingency_table(reg1_right, reg2_right, reg1_left, reg2_left)
        self.assertEqual(expected, result.tolist())


if __name__ == '__main__':
    unittest.main()
