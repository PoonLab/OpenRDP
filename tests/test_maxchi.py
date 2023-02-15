import configparser
import os
import unittest
from io import StringIO
import numpy as np

from openrdp.common import TripletGenerator, Triplet, read_fasta
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
            self.test_short = MaxChi(self.short_align, names, settings=test_settings)

        self.short_triplets = [trp for trp in TripletGenerator(self.short_align, names)]

        # Set up test example 2
        test_settings = {'max_pvalue': '0.05', 'win_size': '40', 'strip_gaps': 'False',
                         'fixed_win_size': 'True', 'num_var_sites': '70', 'frac_var_sites': '0.1'}

        long_seq_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'long.fasta')
        with open(long_seq_path) as test:
            names, test_seqs = read_fasta(test)
            self.long_align = np.array(list(map(list, test_seqs)))
            self.test_long = MaxChi(self.long_align, names, settings=test_settings)

        self.long_triplets = [trp for trp in TripletGenerator(self.long_align, names)]

        # Set up HIV CRF07 test case
        test_settings = {'max_pvalue': '0.05', 'win_size': '100', 'strip_gaps': 'False',
                         'fixed_win_size': 'True', 'num_var_sites': '70', 'frac_var_sites': '0.1'}

        hiv_seq_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'CRF_07_test.fasta')
        with open(hiv_seq_path) as hiv_test:
            names, crf07_seqs = read_fasta(hiv_test)
            self.hiv_align = np.array(list(map(list, crf07_seqs)))
            self.test_hiv = MaxChi(self.hiv_align, names, settings=test_settings)

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

    def test_get_window_positions_long(self):
        # Test long
        exp_reg1_right = 'AAGGCGCGGGAGTGACTTATTTAGAGCCGTCCGCCAGCCAAATCGGGCAT'
        exp_reg1_left = 'AAAAACCTCTACTCGGACGCGCTGCGCGTTTGAAGTCGCCGCGCGCGATC'
        exp_reg2_right = 'AAGGGGCGGGGGTGACTTATCTGGAGCCGTCCGCCAGCCAAATCAGGCAT'
        exp_reg2_left = 'AAAAACATCGACCCGCACCCGCTGCGCGTTTGAAGTCGCCGCACGTGACC'
        exp_reg1 = 'AAAAACCTCTACTCGGACGCGCTGCGCGTTTGAAGTCGCCGCGCGCGATCAAGGCGCGGGAGT' \
                   'GACTTATTTAGAGCCGTCCGCCAGCCAAATCGGGCAT'
        exp_reg2 = 'AAAAACATCGACCCGCACCCGCTGCGCGTTTGAAGTCGCCGCACGTGACCAAGGGGCGGGGGTG' \
                   'ACTTATCTGGAGCCGTCCGCCAGCCAAATCAGGCAT'

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
        exp_reg1 = '------------------------------------------------------------' \
                   '------------------------------------------------------------' \
                   '------------------------------------------------------------' \
                   '--------------------'
        exp_reg2 = 'GAATTCTGGAAGGGTTAATTTACTCTAAGAAAAGGCAAGAGATCCTTGATTTGTGGGTCT' \
                   'ATCACACACAAGGCTACTTCCCTGATTGGCACAACTACACACCAGGACCAGGGACCAGAT' \
                   'TTCCACTGACTTTTGGGTGGTGCTTCAAGCTAGTACCAGTTGACCCAGGGGAAGTAGAAG' \
                   'AGGCCAACGAAGGAGAAGAC'

        seq1 = self.hiv_triplets[0].sequences[0]
        seq2 = self.hiv_triplets[0].sequences[2]

        res_reg1_left, res_reg2_left, res_reg1_right, res_reg2_right, res_reg1, res_reg2 = \
            MaxChi.get_window_positions(seq1, seq2, 0, 200)

        self.assertEqual(exp_reg1_right, ''.join(res_reg1_right))
        self.assertEqual(exp_reg1, ''.join(res_reg1))
        self.assertEqual(exp_reg2, ''.join(res_reg2))
        self.assertEqual(exp_reg1_left, ''.join(res_reg1_left))
        self.assertEqual(exp_reg2_right, ''.join(res_reg2_right))
        self.assertEqual(exp_reg2_left, ''.join(res_reg2_left))

    def test_compute_contingency_table_short(self):
        reg1_right = 'GACGA'
        reg2_right = 'TGGTG'
        reg1_left = 'ATGCT'
        reg2_left = 'AACCT'
        half_win_size = 5

        expected = [[0, 5, 5],
                    [0, 5, 5],
                    [0, 10, 10]]
        result = MaxChi.compute_contingency_table(reg1_right, reg2_right, reg1_left, reg2_left, half_win_size)
        self.assertEqual(expected, result)

    def test_compute_contingency_table_long(self):
        reg1_right = 'AAGGCGCGGGAGTGACTTATTTAGAGCCGTCCGCCAGCCAAATCGGGCAT'
        reg1_left = 'AAAAACCTCTACTCGGACGCGCTGCGCGTTTGAAGTCGCCGCGCGCGATC'
        reg2_right = 'AAGGGGCGGGGGTGACTTATCTGGAGCCGTCCGCCAGCCAAATCAGGCAT'
        reg2_left = 'AAAAACATCGACCCGCACCCGCTGCGCGTTTGAAGTCGCCGCACGTGACC'
        half_win_size = 50

        expected = [[0, 50, 50],
                    [0, 50, 50],
                    [0, 100, 100]]
        result = MaxChi.compute_contingency_table(reg1_right, reg2_right, reg1_left, reg2_left, half_win_size)
        self.assertEqual(expected, result)

    def test_compute_contingency_table_hiv(self):
        reg1_right = '----------------------------------------------------' \
                     '------------------------------------------------'
        reg1_left = '-----------------------------------------------------' \
                    '-----------------------------------------------'
        reg2_right = 'ACCAGGACCAGGGACCAGATTTCCACTGACTTTTGGGTGGTGCTTCAAGCTAG' \
                     'TACCAGTTGACCCAGGGGAAGTAGAAGAGGCCAACGAAGGAGAAGAC'
        reg2_left = 'GAATTCTGGAAGGGTTAATTTACTCTAAGAAAAGGCAAGAGATCCTTGATTTGT' \
                    'GGGTCTATCACACACAAGGCTACTTCCCTGATTGGCACAACTACAC'
        half_win_size = 100

        expected = [[0, 100, 100], [0, 100, 100], [0, 200, 200]]
        result = MaxChi.compute_contingency_table(reg1_right, reg2_right, reg1_left, reg2_left, half_win_size)
        self.assertEqual(expected, result)

    def test_execute_short(self):
        expected = [('A', ('B', 'C'), 2, 25, 0.8780986177504423),
                    ('A', ('B', 'D'), 4, 23, 0.5578254003710748),
                    ('A', ('B', 'E'), 10, 26, 0.5578254003710748),
                    ('A', ('C', 'D'), 9, 22, 0.5578254003710748),
                    ('A', ('C', 'E'), 4, 21, 0.5578254003710748),
                    ('A', ('D', 'E'), 12, 22, 0.5578254003710748),
                    ('B', ('A', 'E'), 10, 20, 0.5578254003710748),
                    ('B', ('C', 'D'), 2, 23, 0.8780986177504423),
                    ('B', ('C', 'E'), 3, 22, 0.8780986177504423),
                    ('B', ('D', 'E'), 4, 20, 0.5578254003710748),
                    ('C', ('D', 'E'), 10, 21, 0.5578254003710748),
                    ('D', ('A', 'C'), 6, 16, 0.5578254003710748),
                    ('D', ('A', 'E'), 12, 21, 0.9553750807650524),
                    ('D', ('C', 'E'), 3, 13, 0.8780986177504423),
                    ('E', ('A', 'B'), 5, 17, 0.8780986177504423),
                    ('E', ('A', 'D'), 6, 17, 0.5578254003710748)]

        for trp in self.short_triplets:
            self.test_short.execute(trp)
        result = self.test_short.merge_breakpoints()
        self.assertEqual(expected, result)

    def test_execute_long(self):
        expected = [('Test1 ', ('Test2', 'Test3'), 475, 518, 0.04042768199451279),
                    ('Test1 ', ('Test2', 'Test4'), 475, 518, 0.04042768199451279),
                    ('Test1 ', ('Test3', 'Test4'), 439, 482, 0.04042768199451279)]

        for trp in self.long_triplets:
            self.test_long.execute(trp)
        result = self.test_long.merge_breakpoints()
        self.assertEqual(expected, result)

    def test_execute_hiv(self):
        expected = [('B', ('07_BC', 'C'), 431, 584, 1.3857496123286656e-17),
                    ('B', ('07_BC', 'C'), 760, 865, 1.1539403917027959e-05),
                    ('B', ('07_BC', 'C'), 4138, 4430, 0.004935527092627837),
                    ('B', ('07_BC', 'C'), 5043, 5246, 0.03739498509166499),
                    ('B', ('07_BC', 'C'), 7673, 7811, 0.03739106813700015),
                    ('B', ('07_BC', 'C'), 8241, 8392, 0.004821912746418881),
                    ('C', ('07_BC', 'B'), 620, 725, 9.2704597756479e-17)]
        for trp in self.hiv_triplets:
            self.test_hiv.execute(trp)
        result = self.test_hiv.merge_breakpoints()
        self.assertEqual(expected, result)


if __name__ == '__main__':
    unittest.main()
