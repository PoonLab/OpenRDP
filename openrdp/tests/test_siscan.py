import configparser
import os
import unittest

import numpy as np

from openrdp.main import read_fasta
from openrdp.scripts.common import generate_triplets, Triplet
from openrdp.scripts.siscan import Siscan


class TestSiscan(unittest.TestCase):

    def setUp(self):
        # Set up test example
        config = configparser.ConfigParser()
        short_cfg_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_short.ini')
        config.read(short_cfg_path)
        test_settings = dict(config.items('Siscan'))

        short_seq_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'short.fasta')
        with open(short_seq_path) as small_test:
            names, test_seqs = read_fasta(small_test)
            self.short_align = np.array(list(map(list, test_seqs)))
            self.test_short = Siscan(self.short_align, names, settings=test_settings)

        self.short_triplets = []
        for trp in generate_triplets(self.short_align):
            self.short_triplets.append(Triplet(self.short_align, names, trp))

        # Set up test example 2
        config = configparser.ConfigParser()
        long_cfg_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_long.ini')
        config.read(long_cfg_path)
        test_settings = dict(config.items('Siscan'))

        long_seq_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'long.fasta')
        with open(long_seq_path) as test:
            names, test_seqs = read_fasta(test)
            self.long_align = np.array(list(map(list, test_seqs)))
            self.test_long = Siscan(self.long_align, names, settings=test_settings)

        self.long_triplets = []
        for trp in generate_triplets(self.long_align):
            self.long_triplets.append(Triplet(self.long_align, names, trp))

        # Set up HIV CRF07 test case
        config = configparser.ConfigParser()
        hiv_cfg_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'default.ini')
        config.read(hiv_cfg_path)
        settings = dict(config.items('Siscan'))

        hiv_seq_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'CRF_07_test.fasta')
        with open(hiv_seq_path) as hiv_test:
            names, crf07_seqs = read_fasta(hiv_test)
            self.hiv_align = np.array(list(map(list, crf07_seqs)))
            self.test_hiv = Siscan(self.hiv_align, names, settings=settings)

        self.hiv_triplets = []
        for trp in generate_triplets(self.hiv_align):
            self.hiv_triplets.append(Triplet(self.hiv_align, names, trp))

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

    def test_execute_short(self):
        expected = [('A', ('B', 'C'), 2, 11, 0.7787397893226586),
                    ('A', ('B', 'D'), 3, 11, 0.7524567011843551),
                    ('A', ('B', 'E'), 2, 11, 0.7787397893226586),
                    ('A', ('C', 'D'), 2, 11, 0.7685248858128534),
                    ('C', ('A', 'E'), 2, 11, 0.804804577675132),
                    ('A', ('D', 'E'), 2, 11, 0.7260095477693602),
                    ('C', ('B', 'D'), 2, 11, 0.7838806032521947),
                    ('B', ('C', 'E'), 3, 11, 0.7750581068141527),
                    ('B', ('D', 'E'), 3, 11, 0.745172605746384),
                    ('C', ('D', 'E'), 2, 11, 0.7465417522802897)]

        for trp in self.short_triplets:
            self.test_short.execute(trp)
        result = self.test_short.merge_breakpoints()
        self.assertEqual(expected, result)

    def test_execute_long(self):
        expected = [('Test1 ', ('Test2', 'Test3'), 2, 55, 0.7521364874425969),
                    ('Test1 ', ('Test2', 'Test4'), 2, 55, 0.7669294010627822),
                    ('Test1 ', ('Test3', 'Test4'), 2, 55, 0.7651446184495414),
                    ('Test2', ('Test3', 'Test4'), 2, 55, 0.7651446184495414)]

        for trp in self.long_triplets:
            self.test_long.execute(trp)
        result = self.test_long.merge_breakpoints()
        self.assertEqual(expected, result)

    def test_execute_hiv(self):
        expected = []   # Breakpoints have p_values that are too large (above threshold)
        for trp in self.hiv_triplets:
            self.test_hiv.execute(trp)
        result = self.test_hiv.merge_breakpoints()
        self.assertEqual(expected, result)


if __name__ == '__main__':
    unittest.main()
