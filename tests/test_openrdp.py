import os
import unittest

from openrdp import *

SHORT_SEQ_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'short.fasta')


class TestMain(unittest.TestCase):

    def test_convert_fasta(self):
        with open(SHORT_SEQ_PATH) as test_handle:
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
        with open(SHORT_SEQ_PATH) as test_handle:
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


if __name__ == '__main__':
    unittest.main()
