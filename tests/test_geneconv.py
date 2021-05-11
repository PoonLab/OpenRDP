import random
import unittest
import configparser
from scripts.geneconv import GeneConv


class TestGeneConv(unittest.TestCase):
    def setUp(self):
        # Test with default settings
        self.config = configparser.ConfigParser()
        self.config.read('default_config.ini')
        self.settings = self.config['Geneconv']
        self.gc_default = GeneConv(self.settings)

        # Test with modified parameters
        self.test_config = configparser.ConfigParser()
        self.test_config.read('test_config.ini')
        self.test_settings = self.test_config['Geneconv']
        self.gc_test = GeneConv(self.test_settings)

    def test_set_and_validate_options(self):
        self.assertEqual(1, self.gc_default.gscale)
        self.assertEqual(False, self.gc_default.ignore_indels)
        self.assertEqual(1, self.gc_default.min_length)
        self.assertEqual(2, self.gc_default.min_poly)
        self.assertEqual(2, self.gc_default.min_score)
        self.assertEqual(1, self.gc_default.max_overlap)

        self.assertEqual(-1, self.gc_test.gscale)
        self.assertEqual(True, self.gc_test.ignore_indels)
        self.assertEqual(1, self.gc_test.min_length)
        self.assertEqual(2, self.gc_test.min_poly)
        self.assertEqual(-2, self.gc_test.min_score)
        self.assertEqual(1, self.gc_test.max_overlap)

    def test_parse_results(self):
        random.seed(9)
        # Test small example with default settings
        self.gc_default.execute('test.fasta')
        expected = {('1', '204'): ['Test2;Test3', '0.0000', '0.00002', 'GI'],
                    ('151', '193'): ['Test4', '0.0085', '0.02202', 'GO'],
                    ('151', '195'): ['Test1;Test3', '0.0004', '0.00210', 'GI'],
                    ('203', '507'): ['Test1;Test2', '0.0016', '0.00829', 'GI'],
                    ('539', '759'): ['Test1;Test2', '0.0240', '0.15378', 'AI'],
                    ('56', '170'): ['Test1', '0.0110', '0.02728', 'GO']}
        result = self.gc_default.parse_output('test.frags')

        # Check the breakpoint locations as p-values can vary
        self.assertEqual(expected.keys(), result.keys())

        # Test small example with modified parameters
        self.gc_test.execute('test.fasta')
        result = self.gc_test.parse_output('test.frags')
        # Check the breakpoint locations as p-values can vary
        self.assertEqual(expected.keys(), result.keys())

        # Test BC_07 with default settings
        self.gc_default.execute('CRF_07_test.fasta')
        expected = {('1', '1229'): ['B', '0.0181', '0.02098', 'GO'],
                    ('2193', '2365'): ['C', '0.0000', '0.00017', 'GO'],
                    ('4416', '5139'): ['B', '0.0055', '0.00758', 'GO'],
                    ('5700', '5771'): ['C', '0.0167', '0.01917', 'GO'],
                    ('5927', '5957'): ['C', '0.0411', '0.06170', 'GO'],
                    ('8219', '8322'): ['C;07_BC', '0.0144', '0.03100', 'GI'],
                    ('8972', '9114'): ['B;07_BC', '0.0000', '0.00000', 'GI'],
                    ('8972', '9118'): ['C', '0.0000', '0.00000', 'GO']}
        result = self.gc_default.parse_output('CRF_07_test.frags')
        # Check the breakpoint locations as p-values can vary
        self.assertEqual(expected.keys(), result.keys())
