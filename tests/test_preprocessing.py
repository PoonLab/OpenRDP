import unittest
import random
import numpy as np
from scripts.preprocessing import Sequence, Alignment

S1 = Sequence('Test1', 'ATGC----------', 0)
S2 = Sequence('Test2', 'AT--AT--CG--GC', 1)
test_seqs = [S1, S2]


class TestSequence(unittest.TestCase):
    def setUp(self):
        random.seed(0)
        self.s1 = S1
        self.s2 = S2

    def test_gap_at_pos(self):
        result = self.s1.gap_at_pos(1)
        self.assertEqual(False, result)

        result = self.s1.gap_at_pos(4)
        self.assertEqual(True, result)

    def test_find_gaps(self):
        expected = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13]
        result = self.s1.find_gaps()
        self.assertEqual(expected, result)

        expected = [2, 3, 6, 7, 10, 11]
        result = self.s2.find_gaps()
        self.assertEqual(expected, result)

    def test_encode_seq(self):
        expected = np.array([[195], [165], [0], [0], [0], [0], [0]], dtype='uint8')
        result = self.s1.encode_seq()
        self.assertEqual(expected, result)


class TestAlignment(unittest.TestCase):
    def setUp(self):
        self.A1 = Alignment(test_seqs)

    def test_pairwise_distances(self):
        expected = 0.571
        result = round(self.A1.pairwise_distances())
        self.assertEqual(expected, result)
