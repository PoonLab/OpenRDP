import os
import unittest

from openrdp.threeseq import ThreeSeq


class Test3Seq(unittest.TestCase):
    def setUp(self):
        long_seq_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'long.fasta')
        self.long_test = ThreeSeq(long_seq_path)
        crf_seq_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'CRF_07_test.fasta')
        self.crf_test = ThreeSeq(crf_seq_path)

    def test_parse_results_long(self):
        result = self.long_test.execute()
        expected = [('Test3', ('Test1', 'Test2'), '202', '787', '5.982096e-10'),
                    ('Test2', ('Test3', 'Test4'), '181', '787', '5.294757e-06')]
        self.assertEqual(expected, result)

    def test_parse_results_crf(self):
        result = self.crf_test.execute()
        expected = [('C', ('07_BC', 'B'), '1989', '2644', '1.270678e-11')]
        self.assertEqual(expected, result)


if __name__ == '__main__':
    unittest.main()
