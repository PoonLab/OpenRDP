import unittest
from scripts.threeseq import ThreeSeq


class Test3Seq(unittest.TestCase):
    def setUp(self):
        self.small_test = ThreeSeq('long.fasta')
        self.crf_07_test = ThreeSeq('CRF_07_test.fasta')

    def test_parse_results(self):
        self.small_test.execute()
        expected = [[['Test3', 'Test1', 'Test2'], '0.000000000025', '5.982096e-10', ['202-204 & 742-759', '202-204 & 784-787']],
                    [['Test2', 'Test4', 'Test3'], '0.000000220615', '5.294757e-06', ['181-193 & 742-750', '181-193 & 784-787']]]
        result = self.small_test.parse_output('long.fasta.3s.rec')
        self.assertEqual(expected, result)

        self.crf_07_test.execute()
        expected = [[['C', 'B', '07_BC'],
                     '0.000000000002',
                     '1.270678e-11',
                     ['1989-2030 & 2617-2617', '1989-2030 & 2636-2644']]]
        result = self.crf_07_test.parse_output('CRF_07_test.fasta.3s.rec')
        self.assertEqual(expected, result)
