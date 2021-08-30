import configparser
import itertools
import os
import unittest

import numpy as np

from openrdp import read_fasta
from scripts.bootscan import Bootscan
from scripts.common import Triplet, generate_triplets


class TestBootscan(unittest.TestCase):

    def setUp(self):
        # Set up test example
        config = configparser.ConfigParser()
        short_cfg_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_short.ini')
        config.read(short_cfg_path)
        test_settings = dict(config.items('Bootscan'))

        short_seq_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'short.fasta')
        with open(short_seq_path) as small_test:
            names, test_seqs = read_fasta(small_test)
            self.short_align = np.array(list(map(list, test_seqs)))
            self.test_short = Bootscan(self.short_align, names, settings=test_settings, quiet=True)

        self.short_triplets = []
        for trp in generate_triplets(self.short_align):
            self.short_triplets.append(Triplet(self.short_align, names, trp))

        # Set up test example 2
        config = configparser.ConfigParser()
        long_cfg_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_long.ini')
        config.read(long_cfg_path)
        test_settings = dict(config.items('Bootscan'))

        long_seq_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'long.fasta')
        with open(long_seq_path) as test:
            names, test_seqs = read_fasta(test)
            self.long_align = np.array(list(map(list, test_seqs)))
            self.test_long = Bootscan(self.long_align, names, settings=test_settings, quiet=True)

        self.long_triplets = []
        for trp in generate_triplets(self.long_align):
            self.long_triplets.append(Triplet(self.long_align, names, trp))

        # Set up HIV CRF07 test case
        config = configparser.ConfigParser()
        hiv_cfg_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'default.ini')
        config.read(hiv_cfg_path)
        settings = dict(config.items('Bootscan'))

        hiv_seq_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'CRF_07_test.fasta')
        with open(hiv_seq_path) as hiv_test:
            names, crf07_seqs = read_fasta(hiv_test)
            self.hiv_align = np.array(list(map(list, crf07_seqs)))
            self.test_hiv = Bootscan(self.hiv_align, names, settings=settings, quiet=True)

        self.hiv_triplets = []
        for trp in generate_triplets(self.hiv_align):
            self.hiv_triplets.append(Triplet(self.hiv_align, names, trp))

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

    def test_percent_diff_short(self):
        # Generate all pairs of sequences (5 sequences)
        # (0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)
        pairs = list(itertools.combinations(range(self.test_short.align.shape[0]), 2))

        expected = [0.7272727272727273, 0.8181818181818182, 0.6818181818181818, 0.7727272727272727, 0.7727272727272727,
                    0.6363636363636364, 0.6818181818181818, 0.6818181818181818, 0.6363636363636364, 0.5909090909090909]
        for i, pair in enumerate(pairs):
            s1 = self.test_short.align[pair[0]]
            s2 = self.test_short.align[pair[1]]
            result = self.test_short.percent_diff(s1, s2)
            self.assertEqual(expected[i], result)

    def test_percent_diff_long(self):
        # Generate all pairs of sequences (4 sequences)
        # (0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)
        pairs = list(itertools.combinations(range(self.test_long.align.shape[0]), 2))

        expected = [0.054637865311308764, 0.15374841168996187, 0.13468869123252858, 0.1207115628970775,
                    0.12198221092757307, 0.08005082592121983]
        for i, pair in enumerate(pairs):
            s1 = self.test_long.align[pair[0]]
            s2 = self.test_long.align[pair[1]]
            result = self.test_long.percent_diff(s1, s2)
            self.assertEqual(expected[i], result)

    def test_percent_diff_hiv(self):
        # Generate all pairs of sequences (3 sequences)
        # (0, 1), (0, 2), (1, 2)
        pairs = list(itertools.combinations(range(self.test_hiv.align.shape[0]), 2))

        expected = [0.1383217410814449, 0.12356954225352113, 0.07229984475493458]
        for i, pair in enumerate(pairs):
            s1 = self.test_hiv.align[pair[0]]
            s2 = self.test_hiv.align[pair[1]]
            result = self.test_hiv.percent_diff(s1, s2)
            self.assertEqual(expected[i], result)

    def test_jc_distance_short(self):
        # Generate all pairs of sequences (5 sequences)
        # (0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)
        pairs = list(itertools.combinations(range(self.test_short.align.shape[0]), 2))

        expected = [2.6223806710998607, 1, 1.7984214545987776, 1, 1, 1.415302236774285, 1.7984214545987776,
                    1.7984214545987776, 1.415302236774285, 1.1629480593083754]
        for i, pair in enumerate(pairs):
            s1 = self.test_short.align[pair[0]]
            s2 = self.test_short.align[pair[1]]
            result = self.test_short.jc_distance(s1, s2)
            self.assertEqual(expected[i], result)

    def test_jc_distance_long(self):
        # Generate all pairs of sequences (4 sequences)
        # (0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)
        pairs = list(itertools.combinations(range(self.test_long.align.shape[0]), 2))

        expected = [0.056730329671987476, 0.1720578753742531, 0.1484586552588878,
                    0.13161261779022151, 0.13312853543617206, 0.08465351745505377]
        for i, pair in enumerate(pairs):
            s1 = self.test_long.align[pair[0]]
            s2 = self.test_long.align[pair[1]]
            result = self.test_long.jc_distance(s1, s2)
            self.assertEqual(expected[i], result)

    def test_jc_distance_hiv(self):
        # Generate all pairs of sequences (3 sequences)
        # (0, 1), (0, 2), (1, 2)
        pairs = list(itertools.combinations(range(self.test_hiv.align.shape[0]), 2))

        expected = [0.15290008723420206, 0.13502657966855233, 0.0760261989839483]
        for i, pair in enumerate(pairs):
            s1 = self.test_hiv.align[pair[0]]
            s2 = self.test_hiv.align[pair[1]]
            result = self.test_hiv.jc_distance(s1, s2)
            self.assertEqual(expected[i], result)

    def test_execute_short(self):
        expected = {('A', 'B', 'C'): [],
                    ('A', 'B', 'D'): [],
                    ('A', 'B', 'E'): [],
                    ('A', 'C', 'D'): [],
                    ('A', 'C', 'E'): [],
                    ('A', 'D', 'E'): [],
                    ('B', 'C', 'D'): [],
                    ('B', 'C', 'E'): [],
                    ('B', 'D', 'E'): [],
                    ('C', 'D', 'E'): []}
        result = self.test_short.execute(self.short_triplets, True)
        self.assertEqual(expected, result)

    def test_execute_long(self):
        expected = {('Test1 ', 'Test2', 'Test3'): [],
                    ('Test1 ', 'Test2', 'Test4'): [],
                    ('Test1 ', 'Test3', 'Test4'): [('Test1 ', 155, 195, 0.0, 0.0),
                                                   ('Test1 ', 550, 560, 0.0, 0.0),
                                                   ('Test3', 155, 195, 0.0, 0.0),
                                                   ('Test3', 550, 560, 24.974927890672568, 99.89971156269027)],
                    ('Test2', 'Test3', 'Test4'): []}
        result = self.test_long.execute(self.long_triplets, True)
        self.assertEqual(expected, result)

    def test_execute_hiv(self):
        expected = {('B', 'C', '07_BC'): [('B', 100, 560, 7.549245243949326e-18, 7.549245243949326e-18),
                                          ('B', 1280, 1340, 3.599366796863919, 3.599366796863919),
                                          ('B', 2060, 2560, 3.303591858532724e-10, 3.303591858532724e-10),
                                          ('B', 3000, 3160, 5.126961605427604e-06, 5.126961605427604e-06),
                                          ('B', 8960, 9140, 0.00082344373678815, 0.00082344373678815),
                                          ('C', 100, 560, 2.467762123557746, 2.467762123557746),
                                          ('C', 580, 1260, 0.5929199341651085, 0.5929199341651085),
                                          ('C', 2620, 2900, 0.007864598512691619, 0.007864598512691619),
                                          ('C', 3260, 5680, 3.765271523156985e-33, 3.765271523156985e-33),
                                          ('07_BC', 580, 1260, 1.4656210312101543e-08, 1.4656210312101543e-08),
                                          ('07_BC', 1280, 1340, 3.929355411778469, 3.929355411778469),
                                          ('07_BC', 2060, 2560, 3.419613887366981e-05, 3.419613887366981e-05),
                                          ('07_BC', 2620, 2900, 7.539856255086794e-06, 7.539856255086794e-06),
                                          ('07_BC', 3000, 3160, 0.21399128509669785, 0.21399128509669785),
                                          ('07_BC', 3260, 5680, 6.176983808734789e-40, 6.176983808734789e-40),
                                          ('07_BC', 8960, 9140, 0.05488355576786701, 0.05488355576786701)]}
        result = self.test_hiv.execute(self.hiv_triplets, True)
        self.assertEqual(expected, result)


if __name__ == '__main__':
    unittest.main()
