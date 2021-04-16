import unittest
from scripts.run_scans import *


class TestRunScans(unittest.TestCase):
    def test_run_scans(self):
        aln = [['>Test1', 'TCGCGACGTCAA'],
               ['>Test2', 'ATGCGGATGGGG'],
               ['>Test3', 'CATGCGACTACG']]