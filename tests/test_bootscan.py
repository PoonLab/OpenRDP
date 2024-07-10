import unittest
import numpy as np
import os

from openrdp.bootscan import Bootscan
from openrdp.common import TripletGenerator, read_fasta


basepath = os.path.dirname(os.path.abspath(__file__))
SHORT_INFILE = os.path.join(basepath, 'short.fasta')
LONG_INFILE = os.path.join(basepath, 'long.fasta')
HIV_INFILE = os.path.join(basepath, 'CRF_07_test.fasta')

SHORT_REFFILE = os.path.join(basepath, 'short_ref.fasta')
LONG_REFFILE = os.path.join(basepath, 'long_ref.fasta')
HIV_REFFILE = os.path.join(basepath, 'CRF_07_ref.fasta')


class TestBootscan(unittest.TestCase):
    def setUp(self):        
        # Set up short test fixtures
        test_settings = {'max_pvalue': '1.0', 'win_size': '6', 'step_size': '1', 'num_replicates': '3',
                         'random_seed': '3', 'cutoff_percentage': '0.7', 'scan': 'distances', 'np': '2'}
        with open(SHORT_INFILE) as short_test, open(SHORT_REFFILE) as short_ref:
            names, test_seqs = read_fasta(short_test)
            ref_names, ref_seqs = read_fasta(short_ref)
        self.short_align = np.array(list(map(list, test_seqs)))
        self.short_names = names
        self.short_ref_names = ref_names
        self.short_align_r = np.array(list(map(list, ref_seqs)))
        self.test_short = Bootscan(self.short_align, settings=test_settings, verbose=False,
                                   ref_align=self.short_align_r)
        temp = []
        for i in range(0, self.short_align.shape[1], self.test_short .step_size):
            temp.append(self.test_short.scan(i))
        self.test_short.dt_matrix_file = self.test_short.collate_scanning_phase(temp)

        # Set up long test fixtures
        test_settings = {'max_pvalue': '0.05', 'win_size': '50', 'step_size': '5', 'num_replicates': '100',
                         'random_seed': '3', 'cutoff_percentage': '0.7', 'scan': 'distances', 'np': '2'}
        with open(LONG_INFILE) as long_test, open(LONG_REFFILE) as long_ref:
            names, test_seqs = read_fasta(long_test)
            ref_names, ref_seqs = read_fasta(long_ref)
        self.long_align = np.array(list(map(list, test_seqs)))
        self.long_names = names
        self.long_ref_names = ref_names
        self.long_align_r = np.array(list(map(list, ref_seqs)))
        self.test_long = Bootscan(self.long_align, settings=test_settings, verbose=False,
                                  ref_align=self.long_align_r)
        temp = []
        for i in range(0, self.long_align.shape[1], self.test_long.step_size):
            temp.append(self.test_long.scan(i))
        self.test_long.dt_matrix_file = self.test_long.collate_scanning_phase(temp)

        # Set up HIV test fixtures
        test_settings = {'max_pvalue': '0.05', 'win_size': '200', 'step_size': '20', 'num_replicates': '100',
                         'random_seed': '3', 'cutoff_percentage': '0.7', 'scan': 'distances', 'np': '2',
                         'p_value_calculation': 'binomial', 'model': 'Jukes-Cantor'}
        with open(HIV_INFILE) as hiv_test, open(HIV_REFFILE) as hiv_ref:
            names, test_seqs = read_fasta(hiv_test)
            ref_names, ref_seqs = read_fasta(hiv_ref)
        self.hiv_align = np.array(list(map(list, test_seqs)))
        self.hiv_align_r = np.array(list(map(list, ref_seqs)))
        self.hiv_names = names
        self.hiv_ref_names = ref_names
        self.test_hiv = Bootscan(self.hiv_align, settings=test_settings, verbose=False,
                                 ref_align=self.hiv_align_r)
        temp = []
        for i in range(0, self.hiv_align.shape[1], self.test_hiv.step_size):
            temp.append(self.test_hiv.scan(i))
        self.test_hiv.dt_matrix_file = self.test_hiv.collate_scanning_phase(temp)
        # self.test_hiv_r = Bootscan(self.hiv_align, self.hiv_align_r, settings=test_settings, verbose=False)

    def tearDown(self):
        if os.path.exists(self.test_short.dt_matrix_file):
            os.remove(self.test_short.dt_matrix_file)

        if os.path.exists(self.test_long.dt_matrix_file):
            os.remove(self.test_long.dt_matrix_file)

        if os.path.exists(self.test_hiv.dt_matrix_file):
            os.remove(self.test_hiv.dt_matrix_file)

    def test_set_and_validate_options(self):
        self.assertEqual(6, self.test_short.win_size)
        self.assertEqual(1, self.test_short.step_size)
        self.assertEqual(3, self.test_short.num_replicates)
        self.assertEqual(3, self.test_short.random_seed)
        self.assertEqual(0.7, self.test_short.cutoff)

        self.assertEqual(50, self.test_long.win_size)
        self.assertEqual(5, self.test_long.step_size)
        self.assertEqual(100, self.test_long.num_replicates)
        self.assertEqual(3, self.test_long.random_seed)
        self.assertEqual(0.7, self.test_long.cutoff)

        self.assertEqual(200, self.test_hiv.win_size)
        self.assertEqual(20, self.test_hiv.step_size)
        self.assertEqual(100, self.test_hiv.num_replicates)
        self.assertEqual(3, self.test_hiv.random_seed)
        self.assertEqual(0.7, self.test_hiv.cutoff)

    def test_execute_short(self):
        expected = []   # No breakpoints found

        triplets = TripletGenerator(self.short_align, self.short_names,
                                    ref_align=self.short_align_r,
                                    ref_names=self.short_ref_names)

        temp = []
        for trp_count, triplet in enumerate(triplets):
            temp.append(self.test_short.execute((trp_count, triplet)))

        self.test_short.raw_results = [l for j in temp for l in j]
        result = self.test_short.merge_breakpoints()
        self.assertEqual(expected, result)

    def test_execute_long(self):
        expected = [('Test3', ('X', 'Z'), 730, 765, 1.7349144672495443e-06),
                    ('Test4', ('X', 'Z'), 740, 765, 9.437118236855531e-05),
                    ('X', ('Test3', 'Z'), 730, 765, 9.367208339129421e-07),
                    ('X', ('Test4', 'Z'), 675, 720, 1.5115935397493938e-08),
                    ('X', ('Test4', 'Z'), 740, 765, 5.2804433825752846e-05)]

        triplets = TripletGenerator(self.long_align, self.long_names,
                                    ref_align=self.long_align_r,
                                    ref_names=self.long_ref_names)
        temp = []
        for trp_count, triplet in enumerate(triplets):
            temp.append(self.test_long.execute((trp_count, triplet)))

        self.test_long.raw_results = [l for j in temp for l in j]
        result = self.test_long.merge_breakpoints()
        self.assertEqual(expected, result)

    def test_execute_hiv(self):
        expected = [('07_BC', ('D', 'N'), 100, 120, 0.00016114151306454157),
                    ('07_BC', ('D', 'N'), 3160, 3280, 1.458628445847309e-22),
                    ('07_BC', ('D', 'S'), 2620, 2640, 0.0002501551852276991),
                    ('07_BC', ('N', 'S'), 4360, 4440, 1.1667614687407429e-15),
                    ('07_BC', ('N', 'S'), 4660, 4720, 6.312989816178787e-12),
                    ('07_BC', ('N', 'S'), 9280, 9480, 4.6500240269794826e-38),
                    ('B', ('D', 'S'), 3140, 3160, 0.0001399725672939798),
                    ('B', ('N', 'S'), 100, 380, 2.982401314971437e-57),
                    ('B', ('N', 'S'), 4320, 4440, 5.954071357483718e-25),
                    ('B', ('N', 'S'), 9280, 9320, 8.412429127473028e-09),
                    ('C', ('D', 'N'), 2380, 2420, 1.1791230407234813e-08),
                    ('C', ('D', 'N'), 3100, 3180, 1.3834955760112245e-16),
                    ('C', ('D', 'S'), 100, 520, 2.7934091397449077e-80),
                    ('C', ('N', 'S'), 100, 520, 1.3590882600638936e-87),
                    ('C', ('N', 'S'), 4660, 4720, 3.8945231346265203e-13),
                    ('N', ('07_BC', 'D'), 100, 120, 7.520513397326907e-05),
                    ('N', ('07_BC', 'D'), 2580, 2780, 2.8791612866906974e-40),
                    ('N', ('07_BC', 'D'), 3160, 3280, 1.8860642064252804e-24),
                    ('N', ('07_BC', 'D'), 5720, 5760, 1.0720513377765899e-08),
                    ('N', ('B', 'D'), 2560, 2760, 1.6234709623492212e-41),
                    ('N', ('B', 'D'), 5720, 5780, 5.672929964151832e-13),
                    ('N', ('C', 'D'), 2260, 2320, 2.353586510938943e-13),
                    ('N', ('C', 'D'), 2380, 2420, 3.813185790775636e-09),
                    ('N', ('C', 'D'), 2580, 2780, 8.128006333945204e-43),
                    ('N', ('C', 'D'), 3100, 3180, 1.4562505126128315e-17),
                    ('N', ('C', 'D'), 3380, 3440, 2.3592104323839885e-13),
                    ('N', ('C', 'D'), 5720, 5760, 3.730560757117377e-09),
                    ('S', ('07_BC', 'D'), 200, 260, 3.3155881292016285e-12),
                    ('S', ('07_BC', 'D'), 2440, 2460, 0.0001476233362251771),
                    ('S', ('07_BC', 'D'), 2620, 2640, 0.0001492876792339202),
                    ('S', ('07_BC', 'D'), 2720, 2860, 1.6550144001454064e-27),
                    ('S', ('07_BC', 'D'), 8980, 9020, 1.780712548539637e-08),
                    ('S', ('07_BC', 'N'), 2380, 2420, 1.9009310996513495e-08),
                    ('S', ('07_BC', 'N'), 4360, 4440, 3.614352822874917e-16),
                    ('S', ('07_BC', 'N'), 4660, 4720, 2.6212253934165498e-12),
                    ('S', ('07_BC', 'N'), 7860, 7920, 2.610134585332796e-12),
                    ('S', ('07_BC', 'N'), 9220, 9240, 0.00013777059147130865),
                    ('S', ('07_BC', 'N'), 9280, 9480, 2.48357319198434e-39),
                    ('S', ('B', 'D'), 2420, 2540, 1.2876461299158664e-24),
                    ('S', ('B', 'D'), 2720, 2860, 1.9318586489804064e-28),
                    ('S', ('B', 'D'), 3140, 3160, 7.4288692646716e-05),
                    ('S', ('B', 'D'), 7780, 7840, 1.3261014181654115e-12),
                    ('S', ('B', 'N'), 100, 380, 2.2843135345823474e-56),
                    ('S', ('B', 'N'), 760, 780, 0.00010570179471101775),
                    ('S', ('B', 'N'), 2380, 2420, 1.125108332719097e-08),
                    ('S', ('B', 'N'), 4320, 4440, 1.4247950821532757e-24),
                    ('S', ('B', 'N'), 9280, 9320, 1.1252450518809477e-08),
                    ('S', ('C', 'D'), 100, 520, 2.4486956382114977e-83),
                    ('S', ('C', 'D'), 7760, 7820, 1.571795328991421e-12),
                    ('S', ('C', 'D'), 8960, 9120, 3.3584820710807995e-32),
                    ('S', ('C', 'N'), 100, 520, 1.2769009576196612e-83),
                    ('S', ('C', 'N'), 2400, 2420, 9.241898245445662e-05),
                    ('S', ('C', 'N'), 4660, 4720, 1.4388016130344453e-12),
                    ('S', ('C', 'N'), 4760, 4880, 2.070177779166627e-24),
                    ('S', ('C', 'N'), 6640, 6660, 0.00011279847560181666),
                    ('S', ('C', 'N'), 9080, 9100, 0.00011250132041173)]

        triplets = TripletGenerator(self.hiv_align, self.hiv_names,
                                    ref_align=self.hiv_align_r,
                                    ref_names=self.hiv_ref_names)

        temp = []
        for trp_count, triplet in enumerate(triplets):
            temp.append(self.test_hiv.execute((trp_count, triplet)))

        self.test_hiv.raw_results = [l for j in temp for l in j]
        result = self.test_hiv.merge_breakpoints()
        self.assertEqual(expected, result)


if __name__ == '__main__':
    unittest.main()
