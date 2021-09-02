import configparser
import os
import unittest

import numpy as np

from openrdp import read_fasta
from scripts.common import generate_triplets, Triplet
from scripts.siscan import Siscan


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
        expected = [('A', ('B', 'C'), 2, 14, 0.6936497660301494),
                    ('A', ('B', 'D'), 0, 11, 0.2712988780903053),
                    ('A', ('B', 'E'), 0, 14, 0.017491303540256298),
                    ('A', ('C', 'D'), 0, 14, 0.017491303540256298),
                    ('A', ('C', 'E'), 3, 11, 0.7623133685530908),
                    ('A', ('D', 'E'), 0, 9, 0.16809551850276472),
                    ('B', ('A', 'E'), 3, 11, 0.5823189648769475),
                    ('B', ('C', 'D'), 4, 14, 0.6148083133771495),
                    ('B', ('C', 'E'), 3, 11, 0.7154326277820072),
                    ('B', ('D', 'E'), 1, 13, 0.7164743750877498),
                    ('C', ('A', 'E'), 3, 14, 0.732581251074242),
                    ('C', ('B', 'D'), 1, 8, 0.3121175827892353),
                    ('C', ('B', 'E'), 2, 14, 0.802097907199899),
                    ('C', ('D', 'E'), 2, 13, 0.7465417522802897),
                    ('D', ('A', 'E'), 5, 14, 0.694825913849001),
                    ('E', ('A', 'B'), 3, 11, 0.7154326277820072),
                    ('E', ('A', 'C'), 3, 11, 0.6919841462633497),
                    ('E', ('C', 'D'), 2, 14, 0.7836931309660096)]
        result = self.test_short.execute(self.short_triplets, quiet=True)
        self.assertEqual(expected, result)

    def test_execute_long(self):
        expected = [('Test1 ', ('Test2', 'Test3'), 0, 55, 0.040289150619780056),
                    ('Test1 ', ('Test2', 'Test4'), 0, 58, 0.016337725552167226),
                    ('Test1 ', ('Test3', 'Test4'), 0, 58, 0.013332210495479169),
                    ('Test2', ('Test1 ', 'Test3'), 2, 55, 0.7365813968714297),
                    ('Test2', ('Test1 ', 'Test4'), 2, 58, 0.7657205072921975),
                    ('Test2', ('Test3', 'Test4'), 0, 55, 0.01754957684396985),
                    ('Test3', ('Test1 ', 'Test2'), 1, 55, 0.06817619324278229),
                    ('Test3', ('Test1 ', 'Test4'), 1, 55, 0.06938539363644475),
                    ('Test3', ('Test2', 'Test4'), 1, 55, 0.09766520576061671),
                    ('Test4', ('Test1 ', 'Test2'), 2, 55, 0.7378186976647264),
                    ('Test4', ('Test1 ', 'Test3'), 2, 55, 0.6777401038341523),
                    ('Test4', ('Test2', 'Test3'), 2, 55, 0.723655056240975)]
        result = self.test_long.execute(self.long_triplets, quiet=True)
        self.assertEqual(expected, result)

    def test_execute_hiv(self):
        self.maxDiff = None
        expected = [('B', ('07_BC', 'C'), 5, 208, 0.32316506014686414)]
        result = self.test_hiv.execute(self.hiv_triplets, quiet=True)
        self.assertEqual(expected, result)


if __name__ == '__main__':
    unittest.main()
