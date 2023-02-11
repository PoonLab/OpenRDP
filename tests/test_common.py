import os
import unittest
import itertools

from openrdp import Scanner
from openrdp.common import *
import numpy as np

basepath = os.path.dirname(os.path.abspath(__file__))
SHORT_INFILE = os.path.join(basepath, 'short.fasta')
LONG_INFILE = os.path.join(basepath, 'long.fasta')
HIV_INFILE = os.path.join(basepath, 'CRF_07_test.fasta')


class TestCommon(unittest.TestCase):
    def setUp(self):
        self.scanner = Scanner()
        self.scanner._import_data(SHORT_INFILE)
        self.short_align = self.scanner.alignment
        self.short_names = self.scanner.seq_names

        self.scanner._import_data(LONG_INFILE)
        self.long_align = self.scanner.alignment
        self.long_names = self.scanner.seq_names

        self.scanner._import_data(HIV_INFILE)
        self.hiv_align = self.scanner.alignment
        self.hiv_names = self.scanner.seq_names

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

    def test_reduce_to_unique_seqs(self):
        aln = np.array(['ATGCATTGCGA',
                        'ATGCATTGCGA',
                        'TCGCACGACAA'])
        exp = ['ATGCATTGCGA', 'TCGCACGACAA']
        result = reduce_to_unique_seqs(aln)
        self.assertEqual(sorted(exp), sorted(result))

    def test_calculate_chi2(self):
        c_table = [[338, 363, 701],
                   [125, 156, 281],
                   [463, 519, 982]]
        exp_chi2 = 1.12
        exp_p_value = 0.891
        res_chi2, res_p_value = calculate_chi2(c_table, 500)
        self.assertEqual(exp_chi2, round(res_chi2, 2))
        self.assertEqual(exp_p_value, round(res_p_value, 3))

    def test_percent_diff_short(self):
        # Generate all pairs of sequences (5 sequences)
        # (0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)
        pairs = list(itertools.combinations(range(self.short_align.shape[0]), 2))

        expected = [0.7272727272727273, 0.8181818181818182, 0.6818181818181818, 0.7727272727272727, 0.7727272727272727,
                    0.6363636363636364, 0.6818181818181818, 0.6818181818181818, 0.6363636363636364, 0.5909090909090909]
        for i, pair in enumerate(pairs):
            s1 = self.short_align[pair[0]]
            s2 = self.short_align[pair[1]]
            result = percent_diff(s1, s2)
            self.assertEqual(expected[i], result)

    def test_percent_diff_long(self):
        # Generate all pairs of sequences (4 sequences)
        # (0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)
        pairs = list(itertools.combinations(range(self.long_align.shape[0]), 2))

        expected = [0.054637865311308764, 0.15374841168996187, 0.13468869123252858, 0.1207115628970775,
                    0.12198221092757307, 0.08005082592121983]
        for i, pair in enumerate(pairs):
            s1 = self.long_align[pair[0]]
            s2 = self.long_align[pair[1]]
            result = percent_diff(s1, s2)
            self.assertEqual(expected[i], result)

    def test_percent_diff_hiv(self):
        # Generate all pairs of sequences (3 sequences)
        # (0, 1), (0, 2), (1, 2)
        pairs = list(itertools.combinations(range(self.hiv_align.shape[0]), 2))

        expected = [0.1383217410814449, 0.12356954225352113, 0.07229984475493458]
        for i, pair in enumerate(pairs):
            s1 = self.hiv_align[pair[0]]
            s2 = self.hiv_align[pair[1]]
            result = percent_diff(s1, s2)
            self.assertEqual(expected[i], result)

    def test_jc_distance_short(self):
        # Generate all pairs of sequences (5 sequences)
        # (0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)
        pairs = list(itertools.combinations(range(self.short_align.shape[0]), 2))

        expected = [2.6223806710998607, 1, 1.7984214545987776, 1, 1, 1.415302236774285, 1.7984214545987776,
                    1.7984214545987776, 1.415302236774285, 1.1629480593083754]
        for i, pair in enumerate(pairs):
            s1 = self.short_align[pair[0]]
            s2 = self.short_align[pair[1]]
            result = jc_distance(s1, s2)
            self.assertEqual(expected[i], result)

    def test_jc_distance_long(self):
        # Generate all pairs of sequences (4 sequences)
        # (0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)
        pairs = list(itertools.combinations(range(self.long_align.shape[0]), 2))

        expected = [0.056730329671987476, 0.1720578753742531, 0.1484586552588878,
                    0.13161261779022151, 0.13312853543617206, 0.08465351745505377]
        for i, pair in enumerate(pairs):
            s1 = self.long_align[pair[0]]
            s2 = self.long_align[pair[1]]
            result = jc_distance(s1, s2)
            self.assertEqual(expected[i], result)

    def test_jc_distance_hiv(self):
        # Generate all pairs of sequences (3 sequences)
        # (0, 1), (0, 2), (1, 2)
        pairs = list(itertools.combinations(range(self.hiv_align.shape[0]), 2))

        expected = [0.15290008723420206, 0.13502657966855233, 0.0760261989839483]
        for i, pair in enumerate(pairs):
            s1 = self.hiv_align[pair[0]]
            s2 = self.hiv_align[pair[1]]
            result = jc_distance(s1, s2)
            self.assertEqual(expected[i], result)


class TestTriplet(unittest.TestCase):
    def setUp(self):
        self.scanner = Scanner()
        self.scanner._import_data(SHORT_INFILE)
        self.short_align = self.scanner.alignment
        self.short_names = self.scanner.seq_names

        self.scanner._import_data(LONG_INFILE)
        self.long_align = self.scanner.alignment
        self.long_names = self.scanner.seq_names

        self.scanner._import_data(HIV_INFILE)
        self.hiv_align = self.scanner.alignment
        self.hiv_names = self.scanner.seq_names

        self.short_triplets = []
        for trp in generate_triplets(self.short_align):
            self.short_triplets.append(
                Triplet(self.short_align, self.short_names, trp)
            )

        self.long_triplets = []
        for trp in generate_triplets(self.long_align):
            self.long_triplets.append(
                Triplet(self.long_align, self.long_names, trp)
            )

        self.hiv_triplets = []
        for trp in generate_triplets(self.hiv_align):
            self.hiv_triplets.append(Triplet(self.hiv_align, self.hiv_names, trp))

    def test_get_triplets(self):
        expected = ['ATGCTGACGACGTAGCAGGTAA', 'AACCTTGGTGCGAAATGCAAGT', 'AGCTGACAGCGATGAGCGAATG']
        result = self.short_triplets[0].get_triplets(self.short_align)
        for i, res in enumerate(result):
            self.assertEqual(expected[i], ''.join(res))

    def test_get_trp_names(self):
        expected = [['A', 'B', 'C'],
                    ['A', 'B', 'D'],
                    ['A', 'B', 'E'],
                    ['A', 'C', 'D'],
                    ['A', 'C', 'E'],
                    ['A', 'D', 'E'],
                    ['B', 'C', 'D'],
                    ['B', 'C', 'E'],
                    ['B', 'D', 'E'],
                    ['C', 'D', 'E']]
        for i, trp in enumerate(self.short_triplets):
            result = trp.get_trp_names(['A', 'B', 'C', 'D', 'E'])
            self.assertEqual(expected[i], result)

    def test_get_sequence_name(self):
        expected = 'C'
        result = self.short_triplets[0].get_sequence_name(2)
        self.assertEqual(expected, result)

    def test_remove_monomorphic_sites_short(self):
        # Test with small sequence
        exp_aln = ['TGCTGACGACGTAGCAGGTAA',
                   'ACCTTGGTGCGAAATGCAAGT',
                   'GCTGACAGCGATGAGCGAATG']
        exp_poly_sites = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21]

        trp = self.short_triplets[0]
        res_aln, res_poly_sites = trp.remove_monomorphic_sites()

        for pos in range(res_aln.shape[0]):
            self.assertEqual(exp_aln[pos], ''.join(res_aln[pos]))
        self.assertEqual(exp_poly_sites, res_poly_sites)
        self.assertEqual(len(exp_aln), len(res_aln))

        # Test with long example
        exp_aln = ['ACGCATACTCCCCACCGAGCAGGCATCCTAACTGATGAGCTATATTTGCTCATAGAATGGTCGCGGTACGCAACATCAACCCTCCAGGTCGCTTCACTC'
                   'CGCGACTAGCGTTAATCGGTG',
                   'CTTGGCGCTCGATGTTAGGTGCCGGCGCTAACTAGCGAGCTATATTGGCTCATAGAATGGGCGCGGTATGCAACACCAACCCTCGGAATTATTCTGCTC'
                   'CGCGACTAGCGTTAACCCATG',
                   'AGCCCCATCGGGTGTCGGATGTCGGCGGGCGTGGGCCGCTCCCGCCTCACGGGGCGGCCCGACTAACGCCGCGTGCGGTTGGATCAGGCCATCCCGTCG'
                   'TATTGTCGCTACAGGCTCACA']
        exp_poly_sites = [6, 9, 12, 15, 18, 24, 42, 45, 48, 54, 59, 60, 62, 63, 66, 70, 72, 75, 94, 99, 102, 120, 123,
                          129, 138, 147, 149, 170, 171, 177, 178, 179, 180, 193, 195, 201, 204, 210, 225, 226, 243,
                          249, 258, 267, 273, 297, 299, 347, 348, 351, 354, 361, 363, 366, 369, 372, 375, 376, 378,
                          381, 405, 406, 410, 412, 427, 428, 432, 433, 435, 439, 440, 444, 446, 452, 457, 468, 470,
                          474, 478, 479, 480, 496, 498, 504, 507, 510, 511, 512, 513, 516, 518, 522, 525, 531, 534,
                          537, 552, 566, 567, 582, 588, 589, 609, 615, 626, 636, 642, 720, 721, 723, 729, 738, 741,
                          750, 759, 760, 762, 766, 780, 783]

        trp = self.long_triplets[1]
        res_aln, res_poly_sites = trp.remove_monomorphic_sites()

        for pos in range(res_aln.shape[0]):
            self.assertEqual(exp_aln[pos], ''.join(res_aln[pos]))
        self.assertEqual(exp_poly_sites, res_poly_sites)
        self.assertEqual(len(exp_aln), len(res_aln))

    def test_remove_uninformative_sites(self):
        # Small example
        exp_aln = ['GCTGCGTAGGGT',
                   'CCTTCGAAACAA',
                   'CTGGGATGAGAA']
        exp_infor_sites = [2, 3, 4, 8, 10, 11, 12, 13, 14, 17, 18, 19]
        exp_uninfor_sites = [0, 1, 5, 6, 7, 9, 15, 16, 20, 21]

        trp = self.short_triplets[0]
        res_aln, res_infor_sites, res_uninfor_sites = trp.remove_uninformative_sites()
        for i, seq in enumerate(res_aln):
            self.assertEqual(exp_aln[i], ''.join(seq))
        self.assertEqual(exp_uninfor_sites, res_uninfor_sites)
        self.assertEqual(exp_infor_sites, res_infor_sites)

        # Larger example
        exp_aln = ['ACTACTCCCACCGAGCAGCATCCTAACTGATGAGCTATATTTGCTCATAGAATGGTCGCGGTACGCAACATCAACCCTCCAGGTCGCTTCACTCCGCGA'
                   'CTAGCGTTAATCGGTG',
                   'CGCGCTCGTGTTAGGTGCGGCGCTAACTAGCGAGCTATATTGGCTCATAGAATGGGCGCGGTATGCAACACCAACCCTCGGAATTATTCTGCTCCGCGA'
                   'CTAGCGTTAACCCATG',
                   'ACCATCGGTGTCGGATGCGGCGGGCGTGGGCCGCTCCCGCCTCACGGGGCGGCCCGACTAACGCCGCGTGCGGTTGGATCAGGCCATCCCGTCGTATTG'
                   'TCGCTACAGGCTCACA']

        exp_infor_sites = [6, 15, 24, 42, 45, 48, 54, 59, 62, 63, 66, 70, 72, 75, 94, 99, 102, 123, 129, 138, 147,
                           149, 170, 171, 177, 178, 179, 180, 193, 195, 201, 204, 210, 225, 226, 243, 249, 258, 267,
                           273, 297, 299, 347, 348, 351, 354, 361, 363, 366, 369, 372, 375, 376, 378, 381, 405, 406,
                           410, 412, 427, 428, 432, 433, 435, 439, 440, 444, 446, 452, 457, 468, 470, 474, 478, 479,
                           480, 496, 498, 504, 507, 510, 511, 512, 513, 516, 518, 522, 525, 531, 534, 537, 552, 566,
                           567, 582, 588, 589, 609, 615, 626, 636, 642, 720, 721, 723, 729, 738, 741, 750, 759, 760,
                           762, 766, 780, 783]

        exp_uninfor_sites = [0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 25, 26, 27,
                             28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 43, 44, 46, 47, 49, 50, 51, 52, 53,
                             55, 56, 57, 58, 60, 61, 64, 65, 67, 68, 69, 71, 73, 74, 76, 77, 78, 79, 80, 81, 82, 83, 84,
                             85, 86, 87, 88, 89, 90, 91, 92, 93, 95, 96, 97, 98, 100, 101, 103, 104, 105, 106, 107, 108,
                             109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 124, 125, 126, 127,
                             128, 130, 131, 132, 133, 134, 135, 136, 137, 139, 140, 141, 142, 143, 144, 145, 146, 148,
                             150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167,
                             168, 169, 172, 173, 174, 175, 176, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191,
                             192, 194, 196, 197, 198, 199, 200, 202, 203, 205, 206, 207, 208, 209, 211, 212, 213, 214,
                             215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 227, 228, 229, 230, 231, 232, 233, 234,
                             235, 236, 237, 238, 239, 240, 241, 242, 244, 245, 246, 247, 248, 250, 251, 252, 253, 254,
                             255, 256, 257, 259, 260, 261, 262, 263, 264, 265, 266, 268, 269, 270, 271, 272, 274, 275,
                             276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293,
                             294, 295, 296, 298, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313,
                             314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331,
                             332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 349, 350, 352,
                             353, 355, 356, 357, 358, 359, 360, 362, 364, 365, 367, 368, 370, 371, 373, 374, 377, 379,
                             380, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392,  393, 394, 395, 396, 397, 398,
                             399, 400, 401, 402, 403, 404, 407, 408, 409, 411, 413, 414, 415, 416, 417, 418, 419, 420,
                             421, 422, 423, 424, 425, 426, 429, 430, 431, 434, 436, 437, 438, 441, 442, 443, 445, 447,
                             448, 449, 450, 451, 453, 454, 455, 456, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467,
                             469, 471, 472, 473, 475, 476, 477, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491,
                             492, 493, 494, 495, 497, 499, 500, 501, 502, 503, 505, 506, 508, 509, 514, 515, 517, 519,
                             520, 521, 523, 524, 526, 527, 528, 529, 530, 532, 533, 535, 536, 538, 539, 540, 541, 542,
                             543, 544, 545, 546, 547, 548, 549, 550, 551, 553, 554, 555, 556, 557, 558, 559, 560, 561,
                             562, 563, 564, 565, 568, 569, 570, 571, 572, 573, 574, 575, 576, 577, 578, 579, 580, 581,
                             583, 584, 585, 586, 587, 590, 591, 592, 593, 594, 595, 596, 597, 598, 599, 600, 601, 602,
                             603, 604, 605, 606, 607, 608, 610, 611, 612, 613, 614, 616, 617, 618, 619, 620, 621, 622,
                             623, 624, 625, 627, 628, 629, 630, 631, 632, 633, 634, 635, 637, 638, 639, 640, 641, 643,
                             644, 645, 646, 647, 648, 649, 650, 651, 652, 653, 654, 655, 656, 657, 658, 659, 660, 661,
                             662, 663, 664, 665, 666, 667, 668, 669, 670, 671, 672, 673, 674, 675, 676, 677, 678, 679,
                             680, 681, 682, 683, 684, 685, 686, 687, 688, 689, 690, 691, 692, 693, 694, 695, 696, 697,
                             698, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715,
                             716, 717, 718, 719, 722, 724, 725, 726, 727, 728, 730, 731, 732, 733, 734, 735, 736, 737,
                             739, 740, 742, 743, 744, 745, 746, 747, 748, 749, 751, 752, 753, 754, 755, 756, 757, 758,
                             761, 763, 764, 765, 767, 768, 769, 770, 771, 772, 773, 774, 775, 776, 777, 778, 779, 781,
                             782, 784, 785, 786]

        trp = self.long_triplets[1]
        res_aln, res_infor_sites, res_uninfor_sites = trp.remove_uninformative_sites()

        for i, seq in enumerate(res_aln):
            self.assertEqual(exp_aln[i], ''.join(seq))
        self.assertEqual(exp_uninfor_sites, res_uninfor_sites)
        self.assertEqual(exp_infor_sites, res_infor_sites)

    def testGetWindowSize(self):
        # Test fixed window size
        expected = 3
        result = self.short_triplets[0].get_win_size(offset=0, win_size=3, fixed_win_size=True,
                                                     num_var_sites=3, frac_var_sites=0)
        self.assertEqual(expected, result)

        # Test variable window size
        expected = 16
        result = self.short_triplets[0].get_win_size(offset=0, win_size=8, fixed_win_size=False,
                                                     num_var_sites=3, frac_var_sites=0)
        self.assertEqual(expected, result)

        expected = 3
        result = self.short_triplets[0].get_win_size(offset=0, win_size=14, fixed_win_size=False,
                                                     num_var_sites=3, frac_var_sites=0)
        self.assertEqual(expected, result)

        expected = 22
        result = self.short_triplets[0].get_win_size(offset=0, win_size=14, fixed_win_size=False,
                                                     num_var_sites=0, frac_var_sites=0.05)
        self.assertEqual(expected, result)

        aln = ['TTTTTTTTTTTTTTTTT',
               'TTGCATGCATTTTTTTT',
               'TTTTTTTTTTTTTTTTT']
        align = np.array(list(map(list, aln)))
        triplet = Triplet(align, ['1', '2', '3'], (0, 1, 2))

        expected = 5
        result = triplet.get_win_size(offset=0, win_size=6, fixed_win_size=False,
                                      num_var_sites=0, frac_var_sites=0.20)
        self.assertEqual(expected, result)


if __name__ == '__main__':
    unittest.main()
