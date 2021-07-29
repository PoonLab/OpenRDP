import unittest
from scripts.rdp import RdpMethod
from scripts.main import read_fasta
import configparser
import numpy as np


class TestRdpMethod(unittest.TestCase):

    def setUp(self):
        # Set up test example
        config = configparser.ConfigParser()
        config.read('test_short.ini')
        test_settings = dict(config.items('MaxChi'))

        with open('short.fasta') as small_test:
            names, test_seqs = read_fasta(small_test)
            small_aln = np.array(list(map(list, test_seqs)))
            self.test_short = RdpMethod(small_aln, names, settings=test_settings)

        # Set up test example 2
        config = configparser.ConfigParser()
        config.read('test_long.ini')
        test_settings = dict(config.items('MaxChi'))

        with open('long.fasta') as test:
            names, test_seqs = read_fasta(test)
            aln = np.array(list(map(list, test_seqs)))
            self.test_long = RdpMethod(aln, names, settings=test_settings)

        # Set up HIV CRF07 test case
        config = configparser.ConfigParser()
        config.read('default.ini')
        settings = dict(config.items('MaxChi'))

        with open('CRF_07_test.fasta') as hiv_test:
            names, crf07_seqs = read_fasta(hiv_test)
            hiv_aln = np.array(list(map(list, crf07_seqs)))
            self.test_hiv = RdpMethod(hiv_aln, names, settings=settings)

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