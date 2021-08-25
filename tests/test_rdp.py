import unittest
from scripts.rdp import RdpMethod
from scripts.main import read_fasta
from scripts.common import generate_triplets, Triplet
import configparser
import numpy as np


class TestRdpMethod(unittest.TestCase):

    def setUp(self):
        # Set up test example
        config = configparser.ConfigParser()
        config.read('test_short.ini')
        test_settings = dict(config.items('RDP'))

        with open('short.fasta') as small_test:
            names, test_seqs = read_fasta(small_test)
            self.short_align = np.array(list(map(list, test_seqs)))
            self.test_short = RdpMethod(self.short_align, names, settings=test_settings)

        self.short_triplets = []
        for trp in generate_triplets(self.short_align):
            self.short_triplets.append(Triplet(self.short_align, names, trp))

        # Set up test example 2
        config = configparser.ConfigParser()
        config.read('test_long.ini')
        test_settings = dict(config.items('RDP'))

        with open('long.fasta') as test:
            names, test_seqs = read_fasta(test)
            self.long_align = np.array(list(map(list, test_seqs)))
            self.test_long = RdpMethod(self.long_align, names, settings=test_settings)

        self.long_triplets = []
        for trp in generate_triplets(self.long_align):
            self.long_triplets.append(Triplet(self.long_align, names, trp))

        # Set up HIV CRF07 test case
        config = configparser.ConfigParser()
        config.read('default.ini')
        settings = dict(config.items('RDP'))

        with open('CRF_07_test.fasta') as hiv_test:
            names, crf07_seqs = read_fasta(hiv_test)
            self.hiv_align = np.array(list(map(list, crf07_seqs)))
            self.test_hiv = RdpMethod(self.hiv_align, names, settings=settings)

        self.hiv_triplets = []
        for trp in generate_triplets(self.hiv_align):
            self.hiv_triplets.append(Triplet(self.hiv_align, names, trp))

    def test_set_and_validate_options(self):
        self.assertEqual(6, self.test_short.win_size)
        self.assertEqual(None, self.test_short.reference)
        self.assertEqual(0, self.test_short.min_id)
        self.assertEqual(100, self.test_short.max_id)

        self.assertEqual(40, self.test_long.win_size)
        self.assertEqual(None, self.test_long.reference)
        self.assertEqual(0, self.test_long.min_id)
        self.assertEqual(100, self.test_long.max_id)

        self.assertEqual(60, self.test_hiv.win_size)
        self.assertEqual(None, self.test_hiv.reference)
        self.assertEqual(0, self.test_hiv.min_id)
        self.assertEqual(100, self.test_hiv.max_id)

    def test_identity_short(self):
        reg_ab = np.array([['C', 'G', 'C', 'G', 'A', 'T',],
                           ['G', 'C', 'C', 'T', 'A', 'G',]])
        reg_bc = np.array([['G', 'C', 'C', 'T', 'A', 'G'],
                           ['C', 'C', 'G', 'G', 'T', 'G']])
        reg_ac = np.array([['C', 'G', 'C', 'G', 'A', 'T'],
                           ['C', 'C', 'G', 'G', 'T', 'G']])

        expected_id = (33.33333333333333, 33.33333333333333, 33.33333333333333)
        result_id = self.test_short.pairwise_identity(reg_ab, reg_bc, reg_ac)
        self.assertEqual(expected_id, result_id)

        expected_trps = [['A', 'B', 'C'], ['A', 'B', 'D'], ['A', 'B', 'E'], ['A', 'C', 'D'], ['A', 'C', 'E'],
                         ['A', 'D', 'E'], ['B', 'C', 'D'], ['B', 'C', 'E'], ['B', 'D', 'E'], ['C', 'D', 'E']]
        result_trps = self.test_short.triplet_identity(self.short_triplets)
        for i, res_trps in enumerate(result_trps):
            self.assertEqual(expected_trps[i], res_trps.names)

    def test_identity_long(self):
        reg_ab = np.array([['C', 'C', 'A', 'T', 'T', 'C', 'G', 'C', 'T', 'T', 'A', 'C', 'G', 'T', 'C', 'G', 'C', 'A',
                            'C', 'T', 'C', 'C', 'A', 'G', 'C', 'C', 'C', 'C', 'G', 'C', 'G', 'T', 'G', 'T', 'A', 'A',
                            'T', 'C', 'G', 'T'],
                           ['T', 'C', 'A', 'C', 'T', 'C', 'A', 'T', 'C', 'C', 'G', 'T', 'G', 'C', 'T', 'A', 'T', 'G',
                            'T', 'C', 'C', 'C', 'G', 'G', 'C', 'C', 'C', 'C', 'C', 'T', 'A', 'C', 'G', 'A', 'G', 'G',
                            'C', 'T', 'A', 'C']])
        reg_bc = np.array([['T', 'C', 'A', 'C', 'T', 'C', 'A', 'T', 'C', 'C', 'G', 'T', 'G', 'C', 'T', 'A', 'T', 'G',
                            'T', 'C', 'C', 'C', 'G', 'G', 'C', 'C', 'C', 'C', 'C', 'T', 'A', 'C', 'G', 'A', 'G', 'G',
                            'C', 'T', 'A', 'C'],
                           ['T', 'T', 'G', 'C', 'C', 'A', 'A', 'T', 'C', 'C', 'A', 'C', 'A', 'T', 'T', 'A', 'T', 'G',
                            'T', 'C', 'T', 'A', 'G', 'A', 'T', 'G', 'A', 'T', 'C', 'T', 'A', 'T', 'A', 'A', 'G', 'A',
                            'C', 'C', 'A', 'C']])
        reg_ac = np.array([['C', 'C', 'A', 'T', 'T', 'C', 'G', 'C', 'T', 'T', 'A', 'C', 'G', 'T', 'C', 'G', 'C', 'A',
                            'C', 'T', 'C', 'C', 'A', 'G', 'C', 'C', 'C', 'C', 'G', 'C', 'G', 'T', 'G', 'T', 'A', 'A',
                            'T', 'C', 'G', 'T'],
                           ['T', 'T', 'G', 'C', 'C', 'A', 'A', 'T', 'C', 'C', 'A', 'C', 'A', 'T', 'T', 'A', 'T', 'G',
                            'T', 'C', 'T', 'A', 'G', 'A', 'T', 'G', 'A', 'T', 'C', 'T', 'A', 'T', 'A', 'A', 'G', 'A',
                            'C', 'C', 'A', 'C']])

        expected_id = (32.5, 52.5, 15.0)
        result_id = self.test_long.pairwise_identity(reg_ab, reg_bc, reg_ac)
        self.assertEqual(expected_id, result_id)

        expected_trps = [['Test1 ', 'Test2', 'Test3'], ['Test1 ', 'Test2', 'Test4'],
                         ['Test1 ', 'Test3', 'Test4'], ['Test2', 'Test3', 'Test4']]
        result_trps = self.test_long.triplet_identity(self.long_triplets)
        for i, res_trps in enumerate(result_trps):
            self.assertEqual(expected_trps[i], res_trps.names)

    def test_no_identity(self):
        reg_ab = np.array([['T', 'G', 'A', 'C', 'G', 'T'],
                           ['G', 'A', 'C', 'T', 'A', 'A']])
        reg_bc = np.array([['G', 'A', 'C', 'T', 'A', 'A'],
                           ['C', 'C', 'G', 'G', 'T', 'G']])
        reg_ac = np.array([['T', 'G', 'A', 'C', 'A', 'T'],
                           ['C', 'C', 'G', 'G', 'T', 'G']])

        expected_id = (0.0, 0.0, 0.0)
        result_id = self.test_short.pairwise_identity(reg_ab, reg_bc, reg_ac)
        self.assertEqual(expected_id, result_id)

        # Set up triplet object
        align = np.array([['T', 'G', 'A', 'C', 'G', 'T'],
                          ['G', 'A', 'C', 'T', 'A', 'A'],
                          ['C', 'C', 'G', 'G', 'T', 'G']])

        trps = []
        for trp in generate_triplets(align):
            trps.append(Triplet(align, ['1', '2', '3'], trp))
        no_id = RdpMethod(align, ['1', '2', '3'])

        expected_trps = []
        result_trps = no_id.triplet_identity(trps)
        for i, res_trps in enumerate(result_trps):
            self.assertEqual(expected_trps[i], res_trps.names)

    def test_execute_short(self):
        expected = {('A', 'B', 'C'): [],
                    ('A', 'B', 'D'): [],
                    ('A', 'B', 'E'): [((11, 14), 0.2682960024379628, 2.682960024379628)],
                    ('A', 'C', 'D'): [((1, 5), 0.083859305437408, 0.8385930543740799)],
                    ('A', 'C', 'E'): [((2, 12), 0.00012352588146716804, 0.0012352588146716803)],
                    ('A', 'D', 'E'): [((3, 13), 0.00016844438381886549, 0.0016844438381886549)],
                    ('B', 'C', 'D'): [((2, 10), 0.026774913658655287, 0.26774913658655286)],
                    ('B', 'C', 'E'): [((6, 7), 1.7355371900826446, 17.355371900826448)],
                    ('B', 'D', 'E'): [],
                    ('C', 'D', 'E'): [((4, 5), 2.0826446280991737, 20.82644628099174)]}
        result = self.test_short.execute(self.short_triplets)
        self.assertEqual(expected, result)

    def test_execute_long(self):
        expected = {('Test1 ', 'Test2', 'Test3'): [((6, 15), 7.663881850891308, 30.65552740356523)],
                    ('Test1 ', 'Test2', 'Test4'): [],
                    ('Test1 ', 'Test3', 'Test4'): [((6, 504), 2.7608273396425555e-05, 0.00011043309358570222)],
                    ('Test2', 'Test3', 'Test4'): [((36, 481), 0.00031154368806080143, 0.0012461747522432057)]}
        result = self.test_long.execute(self.long_triplets)
        self.assertEqual(expected, result)

    def test_execute_hiv(self):
        self.maxDiff = None
        expected = {('B', 'C', '07_BC'): [((431, 6337), 0.0, 0.0),
                                          ((6550, 6647), 0.49389675506775216, 0.49389675506775216),
                                          ((6655, 9748), 6.341013016851143e-50, 6.341013016851143e-50)]}
        result = self.test_hiv.execute(self.hiv_triplets)
        self.assertEqual(expected, result)


if __name__ == '__main__':
    unittest.main()
