import unittest
from scripts.main import *
from scripts.common import remove_monomorphic_sites
import numpy as np


class TestMain(unittest.TestCase):

    def test_convert_fasta(self):
        with open('short.fasta') as test_handle:
            exp_names = ['A', 'B', 'C', 'D', 'E']
            exp_seqs = ['ATGCTGACGACGTAGCAGGTAA',
                        'AACCTTGGTGCGAAATGCAAGT',
                        'AGCTGACAGCGATGAGCGAATG',
                        'ATGCGACTAGCTAGCTAGAGGC',
                        'ACCGAGCGATATCGATCGATGA']
            res_names, res_seqs = read_fasta(test_handle)
            self.assertEqual(exp_names, res_names)
            self.assertEqual(exp_seqs, res_seqs)

    def test_valid_alignment(self):
        with open('short.fasta') as test_handle:
            names, aln = read_fasta(test_handle)
            expected = True
            result = valid_alignment(aln)
            self.assertEqual(expected, result)

        aln2 = ['ATGCGCGC',
                'TGACACAATGC']
        expected = False
        result = valid_alignment(aln2)
        self.assertEqual(expected, result)

    def test_valid_chars(self):
        aln = ['ZXCXZATC',
               'ATGCGGATGGGG',
               'TGTTCAGA']
        expected = False
        result = valid_chars(aln)
        self.assertEqual(expected, result)

        aln = ['TCGCGACGTCAA',
               'ATGCGGATGGGG']
        expected = True
        result = valid_chars(aln)
        self.assertEqual(expected, result)


class TestCommon(unittest.TestCase):

    def setUp(self):
        short_infile = 'short.fasta'
        with open(short_infile) as short_handle:
            _, aln = read_fasta(short_handle)
        self.short_align = np.array(list(map(list, aln)))

        test_infile = 'long.fasta'
        with open(test_infile) as test_handle:
            _, aln = read_fasta(test_handle)
        self.long_align = np.array(list(map(list, aln)))

        infile = 'CRF_07_test.fasta'
        with open(infile) as in_handle:
            _, aln = read_fasta(in_handle)
        self.hiv_align = np.array(list(map(list, aln)))

    def test_remove_monomorphic_sites_short(self):
        exp_aln = ['TGCTGACGACGTAGCAGGTAA',
                   'ACCTTGGTGCGAAATGCAAGT',
                   'GCTGACAGCGATGAGCGAATG',
                   'TGCGACTAGCTAGCTAGAGGC',
                   'CCGAGCGATATCGATCGATGA']
        exp_poly_sites = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21]
        res_aln, res_poly_sites = remove_monomorphic_sites(self.short_align)

        for pos in range(res_aln.shape[0]):
            self.assertEqual(exp_aln[pos], ''.join(res_aln[pos]))

        self.assertEqual(exp_poly_sites, res_poly_sites)

    def test_remove_monomorphic_sites_long(self):

        exp_aln = ['ACGCATCACTCCCCACCGAGCAGGCGATCCTAACTGATGAGGCTCATATCTTGGCGCTCATAGAATGGATCGCGGTACCGCAACCACTCAACCCGCTCCAGGTTCGCTTCACGTCCGCGACTCCAGCCCCGCGTGTAATCGGTG',
                   'CTTGGCCGCTCGATGTTAGGTGCCGGGCGCTAACTAGCGAGGCTCATATCTGGGCGCTCATAGAATGGAGCGCGGTATCGCAACCACCCAACCCGCTCGGAATTTATTCTGCGTCCGCGACTCCAGCCCCGCGTGTAACCCATG',
                   'CTTGGCTGCCAGATGTTAGGTGTCAAGCGCTAACTGGCCGATCCGCCGCTTTAATCGTGGGGCGGCCCTGACTAACGTACGAGTTGTCGGTTGTAAATTGGGCCAATCCCACATATATCGTCTAGATGATCTATAAGACCTACA',
                   'AGCCCCCATCGGGTGTCGGATGTCGGGCGGGCGTGGGCCGGCTCCCCGCCCTGGCCACGGGGCGGCCCAGACTAACGCCCGCGTCGCCGGTTGCGGATCAGGCTCATCCCGTGCGTATTGTCCCGGCCCCCTACGAGGCTCACA']

        exp_poly_sites = [6, 9, 12, 15, 18, 24, 36, 42, 45, 48, 54, 59, 60, 62, 63, 66, 70, 72, 75, 94, 99, 102, 120,
                          123, 129, 132, 138, 147, 149, 170, 171, 177, 178, 179, 180, 193, 195, 201, 204, 210, 216, 225,
                          226, 243, 246, 249, 258, 267, 273, 279, 297, 299, 312, 322, 330, 347, 348, 351, 354, 361, 363,
                          366, 369, 372, 375, 376, 378, 381, 402, 405, 406, 410, 412, 427, 428, 432, 433, 435, 438, 439,
                          440, 444, 446, 452, 453, 457, 465, 468, 470, 474, 478, 479, 480, 481, 493, 496, 498, 504, 507,
                          510, 511, 512, 513, 515, 516, 518, 522, 525, 531, 534, 537, 552, 558, 566, 567, 582, 588, 589,
                          609, 615, 626, 636, 637, 640, 642, 663, 678, 681, 696, 708, 720, 721, 723, 729, 732, 738, 741,
                          750, 759, 760, 762, 766, 780, 783]

        res_aln, res_poly_sites = remove_monomorphic_sites(self.long_align)

        for pos in range(res_aln.shape[0]):
            self.assertEqual(exp_aln[pos], ''.join(res_aln[pos]))

        self.assertEqual(exp_poly_sites, res_poly_sites)
