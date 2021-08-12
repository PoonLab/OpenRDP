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

    def test_execute_short(self):
        expected = {('A', 'B'): [(10, 23, 0.5578254003710748)],
                    ('A', 'C'): [(9, 21, 0.8780986177504423)],
                    ('A', 'D'): [(6, 16, 0.5578254003710748), (12, 22, 0.5578254003710748)],
                    ('A', 'E'): [(7, 17, 0.8780986177504423), (14, 23, 1.0)],
                    ('B', 'C'): [(3, 12, 0.8780986177504423), (13, 22, 0.8780986177504423)],
                    ('B', 'D'): [(4, 15, 0.5578254003710748)],
                    ('B', 'E'): [(5, 15, 0.8780986177504423), (10, 20, 0.5578254003710748)],
                    ('C', 'D'): [(10, 20, 0.5578254003710748)],
                    ('C', 'E'): [(4, 13, 0.5578254003710748), (12, 21, 0.5578254003710748)],
                    ('D', 'E'): [(3, 13, 0.8780986177504423), (12, 21, 1.0)]}
        result = self.test_short.execute(self.short_triplets, False)
        self.assertEqual(expected, result)

    def test_execute_long(self):
        expected = {('Test1 ', 'Test2'): [(475, 518, 0.04042768199451279)],
                    ('Test1 ', 'Test3'): [(439, 482, 0.04042768199451279)],
                    ('Test1 ', 'Test4'): [],
                    ('Test2', 'Test3'): [],
                    ('Test2', 'Test4'): [],
                    ('Test3', 'Test4'): []}
        result = self.test_long.execute(self.long_triplets, False)
        self.assertEqual(expected, result)

    def test_execute_hiv(self):
        expected = {('07_BC', 'B'): [(431, 560, 1.3857496123286656e-17),
                                     (4213, 4317, 0.004935527092627837),
                                     (4327, 4430, 0.0038900938364277763),
                                     (4641, 4744, 1.0),
                                     (4748, 4852, 1.0),
                                     (4933, 5036, 0.000536748557278202),
                                     (5043, 5146, 0.03739498509166499),
                                     (5141, 5246, 0.01765020273817586),
                                     (7673, 7776, 0.026776606487834902),
                                     (8167, 8270, 0.01281262547553272),
                                     (8288, 8392, 0.004821912746418881)],
                    ('07_BC', 'C'): [(620, 725, 9.2704597756479e-17),
                                     (760, 865, 8.287948833018469e-09),
                                     (884, 987, 0.000868758267805495),
                                     (9216, 9320, 1.0),
                                     (9464, 9567, 0.04887860165041302)],
                    ('B', 'C'): [(457, 584, 2.5120278983309063e-18),
                                 (620, 725, 1.609327335603253e-15),
                                 (760, 865, 1.1539403917027959e-05),
                                 (877, 980, 1.0),
                                 (4138, 4241, 0.045376001451184415),
                                 (4240, 4343, 0.0010705880196384548),
                                 (4346, 4450, 1.0),
                                 (4641, 4744, 1.0),
                                 (4748, 4852, 1.0),
                                 (4933, 5036, 0.002924799618437854),
                                 (7708, 7811, 0.03739106813700015),
                                 (8241, 8346, 0.005648215229618876),
                                 (9091, 9194, 0.026776606487834902),
                                 (9465, 9568, 1.0)]}
        result = self.test_hiv.execute(self.hiv_triplets, True)
        self.assertEqual(expected, result)


if __name__ == '__main__':
    unittest.main()
