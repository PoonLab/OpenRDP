import configparser
import os
import unittest
import sys
import numpy as np
import os
from pathlib import Path

from openrdp.tests.test_common import SHORT_INFILE, LONG_INFILE, HIV_INFILE
from openrdp.bootscan import Bootscan
from openrdp.common import Triplet, generate_triplets, read_fasta

path = Path(__file__)
root_dir = path.parent.absolute()

class TestBootscan(unittest.TestCase):
    def setUp(self):
        # Set up short test fixtures
        config = configparser.ConfigParser()
        short_cfg_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_short.ini')
        config.read(short_cfg_path)

        test_settings = dict(config.items('Bootscan'))
        print(test_settings)

        small_test = open(SHORT_INFILE)
        names, test_seqs = read_fasta(small_test)
        small_test.close()
        self.short_align = np.array(list(map(list, test_seqs)))
        self.test_short = Bootscan(self.short_align, names, settings=test_settings, quiet=True)

        self.short_triplets = []
        for trp in generate_triplets(self.short_align):
            self.short_triplets.append(Triplet(self.short_align, names, trp))

        # Set up long test fixtures
        config = configparser.ConfigParser()
        long_cfg_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_long.ini')
        config.read(long_cfg_path)

        test_settings = dict(config.items('Bootscan'))
        long_test = open(LONG_INFILE)
        names, test_seqs = read_fasta(long_test)
        self.long_align = np.array(list(map(list, test_seqs)))
        self.test_long = Bootscan(self.long_align, names, settings=test_settings, quiet=True)

        self.long_triplets = []
        for trp in generate_triplets(self.long_align):
            self.long_triplets.append(
                Triplet(self.long_align, names, trp)
            )

        # Set up HIV test fixtures
        config = configparser.ConfigParser()
        hiv_cfg_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'default.ini')
        config.read(hiv_cfg_path)

        test_settings = dict(config.items('Bootscan'))
        hiv_test = open(HIV_INFILE)
        names, test_seqs = read_fasta(hiv_test)
        hiv_test.close()

        self.hiv_align = np.array(list(map(list, test_seqs)))
        self.test_hiv = Bootscan(self.hiv_align, names, settings=test_settings, quiet=True)

        self.hiv_triplets = []
        for trp in generate_triplets(self.hiv_align):
            self.hiv_triplets.append(
                Triplet(self.hiv_align, names, trp)
            )


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
        expected = []   # No breakpoints found
        for trp in self.short_triplets:
            self.test_short.execute(trp)
        result = self.test_short.merge_breakpoints()
        self.assertEqual(expected, result)

    def test_execute_long(self):
        expected = []   # P-value of breakpoints is outside the threshold
        for trp in self.long_triplets:
            self.test_long.execute(trp)
        result = self.test_long.merge_breakpoints()
        self.assertEqual(expected, result)

    def test_execute_hiv(self):
        expected = [('07_BC', ('B', 'C'), 580, 1260, 1.0158213242512536e-09),
                    ('07_BC', ('B', 'C'), 2060, 2560, 1.7427448207965451e-06),
                    ('07_BC', ('B', 'C'), 3000, 3160, 0.003489818124092514),
                    ('07_BC', ('B', 'C'), 8960, 9140, 0.0010069350767726085),
                    ('B', ('07_BC', 'C'), 100, 560, 3.53955031313494e-19),
                    ('B', ('07_BC', 'C'), 2060, 2560, 1.6836162768997676e-11),
                    ('B', ('07_BC', 'C'), 8960, 9140, 1.5107519378439202e-05),
                    ('C', ('07_BC', 'B'), 2620, 2900, 0.00022445087998712194)]

        for trp in self.hiv_triplets:
            self.test_hiv.execute(trp)
        result = self.test_hiv.merge_breakpoints()
        self.assertEqual(expected, result)


if __name__ == '__main__':
    unittest.main()
