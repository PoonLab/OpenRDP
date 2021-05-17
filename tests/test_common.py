import unittest
from scripts.main import read_fasta
from scripts.common import remove_monomorphic_sites, generate_triplets
import numpy as np


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

    def test_generate_triplets_short(self):
        exp_trps = [(0, 1, 2), (0, 1, 3), (0, 1, 4), (0, 2, 3), (0, 2, 4),
                    (0, 3, 4), (1, 2, 3), (1, 2, 4), (1, 3, 4), (2, 3, 4)]
        res_trps = []
        for trp in generate_triplets(self.short_align):
            res_trps.append(trp)
        self.assertEqual(exp_trps, res_trps)

    def test_generate_triplets_long(self):
        exp_trps = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
        res_trps = []
        for trp in generate_triplets(self.long_align):
            res_trps.append(trp)
        self.assertEqual(exp_trps, res_trps)

    def test_generate_triplets_hiv(self):
        exp_trps = [(0, 1, 2)]
        res_trps = []
        for trp in generate_triplets(self.hiv_align):
            res_trps.append(trp)

        self.assertEqual(exp_trps, res_trps)
