import configparser
import os
import random
import unittest

from openrdp.geneconv import GeneConv


class TestGeneConv(unittest.TestCase):
    def setUp(self):
        # Test with default settings
        test_settings = {'indels_as_polymorphisms': 'True', 'mismatch_penalty': '1',
                         'min_len': '1', 'min_poly': '2', 'min_score': '2', 'max_num': '1'}
        self.gc_default = GeneConv(settings=test_settings)

        # Test with modified parameters
        test_settings = {'indels_as_polymorphisms': 'False', 'mismatch_penalty': '1',
                              'min_len': '1', 'min_poly': '2', 'min_score': '2', 'max_num': '1'}
        self.gc_test = GeneConv(settings=test_settings)

    def test_set_and_validate_options(self):
        self.assertEqual(1, self.gc_default.gscale)
        self.assertEqual(False, self.gc_default.ignore_indels)
        self.assertEqual(1, self.gc_default.min_length)
        self.assertEqual(2, self.gc_default.min_poly)
        self.assertEqual(2, self.gc_default.min_score)
        self.assertEqual(1, self.gc_default.max_overlap)

        self.assertEqual(1, self.gc_test.gscale)
        self.assertEqual(True, self.gc_test.ignore_indels)
        self.assertEqual(1, self.gc_test.min_length)
        self.assertEqual(2, self.gc_test.min_poly)
        self.assertEqual(2, self.gc_test.min_score)
        self.assertEqual(1, self.gc_test.max_overlap)

    def test_parse_results(self):
        random.seed(9)
        # Test small example with default settings
        long_seq_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'long.fasta')
        self.gc_default.execute(long_seq_path)
        expected = [('Test2', ['Test3', '-'], (1, 204), '0.00002'),
                    ('Test1', ['Test3', '-'], (151, 195), '0.00210'),
                    ('Test1', ['Test2', '-'], (203, 507), '0.00829'),
                    ('Test1', ['Test2', '-'], (539, 759), '0.15378'),
                    ('Test4', ['-', '-'], (151, 193), '0.02202'),
                    ('Test1', ['-', '-'], (56, 170), '0.02728')]
        res_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'long.frags')
        result = self.gc_default.parse_output(res_path)
        self.assertEqual(expected, result)

        # Test BC_07 with default settings
        hiv_seq_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'CRF_07_test.fasta')
        self.gc_default.execute(hiv_seq_path)
        expected = [('B', ['07_BC', '-'], (8972, 9114), '0.00000'),
                    ('C', ['07_BC', '-'], (4416, 5139), '0.00117'),
                    ('B', ['07_BC', '-'], (5927, 5957), '0.00904'),
                    ('C', ['07_BC', '-'], (8219, 8322), '0.03100'),
                    ('C', ['-', '-'], (8972, 9118), '0.00000'),
                    ('C', ['-', '-'], (2193, 2365), '0.00017'),
                    ('B', ['-', '-'], (4416, 5139), '0.00758'),
                    ('C', ['-', '-'], (5700, 5771), '0.01917'),
                    ('B', ['-', '-'], (1, 1229), '0.02098'),
                    ('C', ['-', '-'], (5927, 5957), '0.06170')]
        res_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'CRF_07_test.frags')
        result = self.gc_default.parse_output(res_path)
        self.assertEqual(expected, result)


if __name__ == '__main__':
    unittest.main()
