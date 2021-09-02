import configparser
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

    def test_execute_short(self):
        expected = []
        result = self.test_short.execute(self.short_triplets, True)
        self.assertEqual(expected, result)

    def test_execute_long(self):
        expected = [('Test1 ', ('Test3', 'Test4'), 155, 195, 0.0),
                    ('Test1 ', ('Test3', 'Test4'), 550, 560, 0.0),
                    ('Test3', ('Test1 ', 'Test4'), 155, 195, 0.0),
                    ('Test3', ('Test1 ', 'Test4'), 550, 560, 99.89971156269027)]
        result = self.test_long.execute(self.long_triplets, True)
        self.assertEqual(expected, result)

    def test_execute_hiv(self):
        expected = [('B', ('07_BC', 'C'), 100, 560, 7.549245243949326e-18),
                    ('B', ('07_BC', 'C'), 1280, 1340, 3.599366796863919),
                    ('B', ('07_BC', 'C'), 3000, 3160, 5.126961605427604e-06),
                    ('B', ('07_BC', 'C'), 8960, 9140, 0.00082344373678815),
                    ('C', ('07_BC', 'B'), 100, 560, 2.467762123557746),
                    ('C', ('07_BC', 'B'), 580, 1260, 0.5929199341651085),
                    ('C', ('07_BC', 'B'), 3260, 5680, 3.765271523156985e-33),
                    ('07_BC', ('B', 'C'), 580, 1260, 1.4656210312101543e-08),
                    ('07_BC', ('B', 'C'), 1280, 1340, 3.929355411778469),
                    ('07_BC', ('B', 'C'), 2060, 2560, 3.419613887366981e-05),
                    ('07_BC', ('B', 'C'), 3000, 3160, 0.21399128509669785),
                    ('07_BC', ('B', 'C'), 3260, 5680, 6.176983808734789e-40),
                    ('07_BC', ('B', 'C'), 8960, 9140, 0.05488355576786701)]
        result = self.test_hiv.execute(self.hiv_triplets, True)
        self.assertEqual(expected, result)


if __name__ == '__main__':
    unittest.main()
