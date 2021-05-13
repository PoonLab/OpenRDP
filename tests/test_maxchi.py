import unittest
from scripts.maxchi import MaxChi
from scripts.main import read_fasta
import configparser
import numpy as np


# 790 1269 C
# 1270 1410 B
# 1411 2080 C
# 2081 2620 B
# 2621 3010 C
# 3011 3310 B
# 3311 5699 C
# 5700 6310 B
# 6311 8796 C
# 8797 9058 B
# 9059 9417 C


class TestMaxChi(unittest.TestCase):

    def setUp(self):
        # Set up test example
        self.test_cfg = configparser.ConfigParser()
        self.test_cfg.read('test_short.ini')
        self.test_settings = self.test_cfg['MaxChi']

        with open('short.fasta') as small_test:
            names, test_seqs = read_fasta(small_test)
            small_aln = np.array(list(map(list, test_seqs)))
            self.test_short = MaxChi(small_aln, names, self.test_settings)

        # Set up test example 2
        self.test_config = configparser.ConfigParser()
        self.test_config.read('test_long.ini')
        self.test2_settings = self.test_config['MaxChi']

        with open('long.fasta') as test:
            names, test_seqs = read_fasta(test)
            aln = np.array(list(map(list, test_seqs)))
            self.test_long = MaxChi(aln, names, self.test2_settings)

        # Set up HIV CRF07 test case
        self.config = configparser.ConfigParser()
        self.config.read('default.ini')
        self.settings = self.config['MaxChi']

        with open('CRF_07_test.fasta') as hiv_test:
            names, crf07_seqs = read_fasta(hiv_test)
            hiv_aln = np.array(list(map(list, crf07_seqs)))
            self.test_hiv = MaxChi(hiv_aln, names, self.settings)

    def test_set_and_validate_options(self):
        self.assertEqual(5, self.test_short.win_size)
        self.assertEqual(False, self.test_short.strip_gaps)
        self.assertEqual(True, self.test_short.fixed_win_size)
        self.assertEqual(5, self.test_short.num_var_sites)
        self.assertEqual(0.1, self.test_short.frac_var_sites)

        self.assertEqual(200, self.test_long.win_size)
        self.assertEqual(False, self.test_long.strip_gaps)
        self.assertEqual(True, self.test_long.fixed_win_size)
        self.assertEqual(70, self.test_long.num_var_sites)
        self.assertEqual(0.1, self.test_long.frac_var_sites)

        self.assertEqual(False, self.test_hiv.strip_gaps)
        self.assertEqual(True, self.test_hiv.fixed_win_size)
        self.assertEqual(70, self.test_hiv.num_var_sites)
        self.assertEqual(0.1, self.test_hiv.frac_var_sites)
