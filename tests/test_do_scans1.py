import unittest
from scripts.do_scans import ThreeSeq, GeneConv, MaxChi
import configparser
import numpy as np


class Test3Seq(unittest.TestCase):
    def setUp(self):
        self.threeseq = ThreeSeq('test.fasta')

    def test_parse_results(self):
        expected = [[['Test3', 'Test1', 'Test2'], '0.000000000025', '5.982096e-10', ['202-204 & 742-759', '202-204 & 784-787']],
                    [['Test2', 'Test4', 'Test3'], '0.000000220615', '5.294757e-06', ['181-193 & 742-750', '181-193 & 784-787']]]
        results = self.threeseq.parse_output("test.fasta.3s.rec")
        self.assertEqual(expected, results)


class TestGeneConv(unittest.TestCase):
    def setUp(self):
        self.config = configparser.ConfigParser()
        self.config.read('default_config.ini')
        self.settings = self.config['Geneconv']
        self.geneconv = GeneConv(self.settings)

        self.test_config = configparser.ConfigParser()
        self.test_config.read('test_config.ini')
        self.test_settings = self.test_config['Geneconv']
        self.geneconv2 = GeneConv(self.test_settings)

    def test_set_and_validate_options(self):
        self.assertEqual(1, self.geneconv.gscale)
        self.assertEqual(False, self.geneconv.ignore_indels)
        self.assertEqual(1, self.geneconv.min_length)
        self.assertEqual(2, self.geneconv.min_poly)
        self.assertEqual(2, self.geneconv.min_score)
        self.assertEqual(1, self.geneconv.max_overlap)

        self.assertEqual(-1, self.geneconv2.gscale)
        self.assertEqual(True, self.geneconv2.ignore_indels)
        self.assertEqual(1, self.geneconv2.min_length)
        self.assertEqual(2, self.geneconv2.min_poly)
        self.assertEqual(-2, self.geneconv2.min_score)
        self.assertEqual(1, self.geneconv2.max_overlap)

    def test_parse_results(self):
        expected = [['Test2;Test3', '0.0000', '0.00002', ('1', '204'), 'GI'],
                    ['Test1;Test3', '0.0006', '0.00210', ('151', '195'), 'GI'],
                    ['Test1;Test2', '0.0019', '0.00829', ('203', '507'), 'GI'],
                    ['Test1;Test2', '0.0222', '0.15378', ('539', '759'), 'AI'],
                    ['Test4', '0.0078', '0.02202', ('151', '193'), 'GO'],
                    ['Test1', '0.0109', '0.02728', ('56', '170'), 'GO']]
        result = self.geneconv.parse_output("test.frags")
        self.assertEqual(expected, result)


class TestMaxChi(unittest.TestCase):
    def setUp(self):

        self.config = configparser.ConfigParser()
        self.config.read('default_config.ini')
        self.settings = self.config['MaxChi']
        self.maxchi = MaxChi(ALN, self.settings)

        self.test_config = configparser.ConfigParser()
        self.test_config.read('test_config.ini')
        self.test_settings = self.test_config['MaxChi']
        self.maxchi2 = MaxChi(ALN, self.test_settings)

    def test_set_and_validate_options(self):
        self.assertEqual(200, self.maxchi.win_size)
        self.assertEqual(False, self.maxchi.strip_gaps)
        self.assertEqual(True, self.maxchi.fixed_win_size)
        self.assertEqual(70, self.maxchi.num_var_sites)
        self.assertEqual(0.1, self.maxchi.frac_var_sites)

        self.assertEqual(200, self.maxchi2.win_size)
        self.assertEqual(False, self.maxchi2.strip_gaps)
        self.assertEqual(True, self.maxchi2.fixed_win_size)
        self.assertEqual(70, self.maxchi2.num_var_sites)
        self.assertEqual(0.1, self.maxchi2.frac_var_sites)

    def test_monomorphic_sites(self):
        exp_aln = None
        exp_poly_sites = None
        res_aln, res_poly_sites = self.maxchi.remove_monomorphic_sites(ALN)
        self.assertEqual(exp_aln, res_aln)
        self.assertEqual(exp_poly_sites, res_poly_sites)

