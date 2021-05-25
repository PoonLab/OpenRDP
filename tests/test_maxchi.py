import unittest
from scripts.maxchi import MaxChi
from scripts.main import read_fasta
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
            small_aln = np.array(list(map(list, test_seqs)))
            self.test_short = MaxChi(small_aln, names, settings=test_settings)

        # Set up test example 2
        config = configparser.ConfigParser()
        config.read('test_long.ini')
        test_settings = dict(config.items('MaxChi'))

        with open('long.fasta') as test:
            names, test_seqs = read_fasta(test)
            aln = np.array(list(map(list, test_seqs)))
            self.test_long = MaxChi(aln, names, settings=test_settings)

        # Set up HIV CRF07 test case
        config = configparser.ConfigParser()
        config.read('default.ini')
        settings = dict(config.items('MaxChi'))

        with open('CRF_07_test.fasta') as hiv_test:
            names, crf07_seqs = read_fasta(hiv_test)
            hiv_aln = np.array(list(map(list, crf07_seqs)))
            self.test_hiv = MaxChi(hiv_aln, names, settings=settings)

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
                    ('B', 'C'): [(43, 160, 0.035258260485873466), (141, 247, 0.01869359769155257)]}
        result = self.test_hiv.execute()
        self.assertEqual(expected, result)
