import os
import unittest

from openrdp.common import read_fasta
from openrdp import Scanner
from io import StringIO


class TestMain(unittest.TestCase):
    def setUp(self):
        self.scanner = Scanner()

    def test_convert_fasta(self):
        test_fasta = ">A\nATGCTGACGACGTAGCAGGTAA\n>B\nAACCTTGGTGCGAAATGCAAGT\n" \
                     ">C\nAGCTGACAGCGATGAGCGAATG\n>D\nATGCGACTAGCTAGCTAGAGGC\n" \
                     ">E\nACCGAGCGATATCGATCGATGA\n"
        exp_names = ['A', 'B', 'C', 'D', 'E']
        exp_seqs = ['ATGCTGACGACGTAGCAGGTAA',
                    'AACCTTGGTGCGAAATGCAAGT',
                    'AGCTGACAGCGATGAGCGAATG',
                    'ATGCGACTAGCTAGCTAGAGGC',
                    'ACCGAGCGATATCGATCGATGA']
        res_names, res_seqs = read_fasta(StringIO(test_fasta))
        self.assertEqual(exp_names, res_names)
        self.assertEqual(exp_seqs, res_seqs)

    def test_valid_alignment(self):
        aln1 = StringIO(">seq1\nATGCGCGC\n>seq2\nATGCTCGC\n")
        self.scanner._import_data(aln1)

        aln2 = StringIO(">seq3\nATGCGCGC\n>seq4\nTGACACAATGC\n")
        with self.assertRaises(SystemExit):
            self.scanner._import_data(aln2)

    def test_valid_chars(self):
        aln3 = StringIO('>seq1\nZXCXZATC\n>seq2\nATCATATC\n>seq3\nTGTTCAGA\n')
        with self.assertRaises(SystemExit):
            self.scanner._import_data(aln3)


if __name__ == '__main__':
    unittest.main()
