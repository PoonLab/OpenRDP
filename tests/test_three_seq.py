import os
import unittest

from scripts.threeseq import ThreeSeq


class Test3Seq(unittest.TestCase):
    def setUp(self):
        long_seq_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'long.fasta')
        self.long_test = ThreeSeq(long_seq_path)
        crf_seq_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'CRF_07_test.fasta')
        self.crf_test = ThreeSeq(crf_seq_path)

    def test_parse_results_long(self):
        self.long_test.execute()
        expected = [[['Test3', 'Test1', 'Test2'], '0.000000000025', '5.982096e-10', ['202-204 & 742-759', '202-204 & 784-787']],
                    [['Test2', 'Test4', 'Test3'], '0.000000220615', '5.294757e-06', ['181-193 & 742-750', '181-193 & 784-787']]]
        long_res_path = os.path.join(os.getcwd(), 'long.fasta.3s.rec')
        print(long_res_path)
        result = self.long_test.parse_output(long_res_path)
        self.assertEqual(expected, result)

    def test_parse_results_crf(self):
        self.crf_test.execute()
        expected = [[['C', 'B', '07_BC'],
                     '0.000000000002',
                     '1.270678e-11',
                     ['1989-2030 & 2617-2617', '1989-2030 & 2636-2644']]]
        crf_res_path = os.path.join(os.getcwd(), 'CRF_07_test.fasta.3s.rec')
        result = self.crf_test.parse_output(crf_res_path)
        self.assertEqual(expected, result)


if __name__ == '__main__':
    unittest.main()
