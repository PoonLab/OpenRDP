import unittest
import numpy as np
import os

from openrdp.bootscan import Bootscan
from openrdp.common import TripletGenerator, Triplet, read_fasta
from itertools import combinations, product


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
        self.short_align_r = np.array(list(map(list, ref_seqs)))
        self.test_short = Bootscan(self.short_align, settings=test_settings, verbose=False)
        self.test_short_r = Bootscan(self.short_align, self.short_align_r, settings=test_settings, verbose=False)
        self.short_names = names
        self.short_names_r = ref_names

        # Set up long test fixtures
        test_settings = {'max_pvalue': '0.05', 'win_size': '50', 'step_size': '5', 'num_replicates': '100',
                         'random_seed': '3', 'cutoff_percentage': '0.7', 'scan': 'distances', 'np': '2'}
        with open(LONG_INFILE) as long_test, open(LONG_REFFILE) as long_ref:
            names, test_seqs = read_fasta(long_test)
            ref_names, ref_seqs = read_fasta(long_ref)
        self.long_align = np.array(list(map(list, test_seqs)))
        self.long_align_r = np.array(list(map(list, ref_seqs)))
        self.test_long = Bootscan(self.long_align, settings=test_settings, verbose=False)
        self.test_long_r = Bootscan(self.long_align, self.long_align_r, settings=test_settings, verbose=False)
        self.long_names = names
        self.long_names_r = ref_names
        
        # Set up HIV test fixtures
        test_settings = {'max_pvalue': '0.05', 'win_size': '200', 'step_size': '20', 'num_replicates': '100',
                         'random_seed': '3', 'cutoff_percentage': '0.7', 'scan': 'distances', 'np': '2',
                         'p_value_calculation': 'binomial', 'model': 'Jukes-Cantor'}
        with open(HIV_INFILE) as hiv_test, open(HIV_REFFILE) as hiv_ref:
            names, test_seqs = read_fasta(hiv_test)
            ref_names, ref_seqs = read_fasta(hiv_ref)
        self.hiv_align = np.array(list(map(list, test_seqs)))
        self.hiv_align_r = np.array(list(map(list, ref_seqs)))
        self.test_hiv = Bootscan(self.hiv_align, settings=test_settings, verbose=False)
        self.test_hiv_r = Bootscan(self.hiv_align, self.hiv_align_r, settings=test_settings, verbose=False)
        self.hiv_names = names
        self.hiv_names_r = ref_names

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

        # Repeating above tests using reference sequencess
        self.assertEqual(6, self.test_short_r.win_size)
        self.assertEqual(1, self.test_short_r.step_size)
        self.assertEqual(3, self.test_short_r.num_replicates)
        self.assertEqual(3, self.test_short_r.random_seed)
        self.assertEqual(0.7, self.test_short_r.cutoff)

        self.assertEqual(50, self.test_long_r.win_size)
        self.assertEqual(5, self.test_long_r.step_size)
        self.assertEqual(100, self.test_long_r.num_replicates)
        self.assertEqual(3, self.test_long_r.random_seed)
        self.assertEqual(0.7, self.test_long_r.cutoff)

        self.assertEqual(200, self.test_hiv_r.win_size)
        self.assertEqual(20, self.test_hiv_r.step_size)
        self.assertEqual(100, self.test_hiv_r.num_replicates)
        self.assertEqual(3, self.test_hiv_r.random_seed)
        self.assertEqual(0.7, self.test_hiv_r.cutoff)

    def test_execute_short(self):
        total_num_trps = sum(1 for _ in combinations(range(self.short_align.shape[0]), 3))
        expected = []   # No breakpoints found
        self.test_short.execute_all(total_combinations=total_num_trps, seq_names=self.short_names)
        result = self.test_short.merge_breakpoints()
        self.assertEqual(expected, result)

        total_num_trps = sum(1 for _ in product(
                list(combinations(range(self.short_align.shape[0]), 1)),
                list(combinations(range(self.short_align_r.shape[0]), 2))
            ))
        expected = []   # No breakpoints found
        self.test_short_r.execute_all(total_combinations=total_num_trps, seq_names=self.short_names, ref_names=self.short_names_r)
        result = self.test_short_r.merge_breakpoints()
        self.assertEqual(expected, result)

    def test_execute_long(self):
        total_num_trps = sum(1 for _ in combinations(range(self.long_align.shape[0]), 3))
        expected = []   # P-value of breakpoints is outside the threshold
        self.test_long.execute_all(total_combinations=total_num_trps, seq_names=self.long_names)
        result = self.test_long.merge_breakpoints()
        self.assertEqual(expected, result)

        total_num_trps = sum(1 for _ in product(
                list(combinations(range(self.long_align.shape[0]), 1)),
                list(combinations(range(self.long_align_r.shape[0]), 2))
            ))
        expected = [('Test1 ', ('W', 'X'), 40, 145, 1.3845385451851523e-18),
                    ('Test1 ', ('W', 'Y'), 25, 150, 5.821450054495946e-21),
                    ('Test1 ', ('W', 'Y'), 565, 805, 1.4088341047479353e-39),
                    ('Test1 ', ('W', 'Z'), 715, 805, 1.2772268245687568e-15),
                    ('Test1 ', ('Y', 'Z'), 155, 195, 1.1918191282617655e-06),
                    ('Test1 ', ('Y', 'Z'), 550, 560, 0.03245925992192347),
                    ('Test2', ('W', 'Y'), 740, 780, 2.4859484734458336e-07),
                    ('Test2', ('W', 'Y'), 795, 805, 0.022343466949465467),
                    ('Test2', ('W', 'Z'), 795, 805, 0.020986093366838553),
                    ('Test2', ('X', 'Z'), 565, 730, 3.832653596513257e-27),
                    ('Test3', ('W', 'X'), 670, 695, 4.3080308597911745e-05),
                    ('Test3', ('W', 'Y'), 795, 805, 0.023286139899502414),
                    ('Test3', ('W', 'Z'), 295, 325, 1.2830681870496585e-05),
                    ('Test3', ('X', 'Y'), 25, 100, 5.860167264224008e-12),
                    ('Test3', ('X', 'Y'), 120, 145, 0.00017880944506494796),
                    ('Test3', ('X', 'Z'), 340, 525, 3.5725634411819413e-31),
                    ('Test3', ('Y', 'Z'), 300, 325, 0.00022935844518956437),
                    ('Test4', ('W', 'X'), 670, 695, 6.678516268095026e-05),
                    ('Test4', ('W', 'Y'), 795, 805, 0.025265094739408644),
                    ('Test4', ('W', 'Z'), 295, 325, 1.9769954398658228e-05),
                    ('Test4', ('X', 'Y'), 25, 100, 5.8781054633729135e-12),
                    ('Test4', ('X', 'Y'), 120, 145, 0.00017880944506494796),
                    ('Test4', ('X', 'Z'), 340, 525, 1.6490332941978137e-30),
                    ('Test4', ('Y', 'Z'), 300, 325, 0.000229228126283948),
                    ('W', ('Test1 ', 'X'), 40, 145, 6.84986365070102e-20),
                    ('W', ('Test1 ', 'Y'), 25, 150, 7.170946834782862e-22),
                    ('W', ('Test1 ', 'Y'), 565, 805, 2.5275914936231947e-41),
                    ('W', ('Test1 ', 'Z'), 300, 325, 4.445931155247956e-05),
                    ('W', ('Test1 ', 'Z'), 655, 695, 9.247950777090451e-08),
                    ('W', ('Test1 ', 'Z'), 715, 805, 2.2646247157225583e-16),
                    ('W', ('Test2', 'Y'), 795, 805, 0.019278518923610127),
                    ('W', ('Test2', 'Z'), 795, 805, 0.016925677601307337),
                    ('W', ('Test3', 'X'), 670, 695, 3.0576449465841275e-05),
                    ('W', ('Test3', 'Y'), 25, 30, 0.04064382184176061),
                    ('W', ('Test3', 'Y'), 795, 805, 0.018467083400591427),
                    ('W', ('Test3', 'Z'), 295, 325, 7.567009623762942e-06),
                    ('W', ('Test4', 'X'), 670, 695, 3.824711050357666e-05),
                    ('W', ('Test4', 'Y'), 25, 30, 0.04546103876022602),
                    ('W', ('Test4', 'Y'), 795, 805, 0.022343466949465467),
                    ('W', ('Test4', 'Z'), 295, 325, 9.782608554226772e-06),
                    ('X', ('Test2', 'Z'), 565, 730, 1.7446904773314294e-28),
                    ('X', ('Test3', 'Y'), 25, 100, 2.0361212562943013e-14),
                    ('X', ('Test4', 'Y'), 25, 100, 4.722003369745256e-14),
                    ('Y', ('Test1 ', 'Z'), 155, 195, 6.790480162411059e-07),
                    ('Y', ('Test1 ', 'Z'), 550, 560, 0.02562413201556452),
                    ('Y', ('Test2', 'W'), 740, 780, 1.4925076386415515e-06),
                    ('Y', ('Test3', 'X'), 120, 145, 0.00010680021028385066),
                    ('Y', ('Test3', 'Z'), 300, 325, 0.00016378282522499786),
                    ('Y', ('Test3', 'Z'), 340, 530, 1.6928839943818628e-29),
                    ('Y', ('Test4', 'X'), 120, 145, 8.868817811576461e-05),
                    ('Y', ('Test4', 'Z'), 300, 325, 0.0001412977004683572),
                    ('Y', ('Test4', 'Z'), 340, 530, 5.517866508189488e-30),
                    ('Z', ('Test1 ', 'W'), 300, 325, 0.0001156362054264059),
                    ('Z', ('Test1 ', 'W'), 655, 695, 4.92656153021994e-07),
                    ('Z', ('Test3', 'X'), 340, 525, 3.508004364601102e-30),
                    ('Z', ('Test3', 'Y'), 340, 530, 1.7492318850812004e-31),
                    ('Z', ('Test4', 'X'), 340, 525, 3.50800436460127e-30),
                    ('Z', ('Test4', 'Y'), 340, 530, 1.2099624005873301e-30)]
        self.test_long_r.execute_all(total_combinations=total_num_trps, seq_names=self.long_names, ref_names=self.long_names_r)
        result = self.test_long_r.merge_breakpoints()
        self.assertEqual(expected, result)

    def test_execute_hiv(self):
        expected = [('07_BC', ('B', 'C'), 580, 1260, 1.0158213242512536e-09),
                    ('07_BC', ('B', 'C'), 1280, 1340, 0.024030305239701163),
                    ('07_BC', ('B', 'C'), 2060, 2560, 1.7427448207965451e-06),
                    ('07_BC', ('B', 'C'), 2620, 2900, 2.1518293256796478e-07),
                    ('07_BC', ('B', 'C'), 3000, 3160, 0.003489818124092514),
                    ('07_BC', ('B', 'C'), 3260, 5680, 1.5236266249248996e-40),
                    ('07_BC', ('B', 'C'), 8960, 9140, 0.0010069350767726085),
                    ('B', ('07_BC', 'C'), 100, 560, 3.53955031313494e-19),
                    ('B', ('07_BC', 'C'), 1280, 1340, 0.022012231965328217),
                    ('B', ('07_BC', 'C'), 2060, 2560, 1.6836162768997676e-11),
                    ('B', ('07_BC', 'C'), 3000, 3160, 8.361164579231644e-08),
                    ('B', ('07_BC', 'C'), 8960, 9140, 1.5107519378439202e-05),
                    ('C', ('07_BC', 'B'), 580, 1260, 0.04109525585896175),
                    ('C', ('07_BC', 'B'), 2620, 2900, 0.00022445087998712194),
                    ('C', ('07_BC', 'B'), 3260, 5680, 9.28749065950454e-34)]

        total_num_trps = sum(1 for _ in combinations(range(self.hiv_align.shape[0]), 3))
        self.test_hiv.execute_all(total_combinations=total_num_trps, seq_names=self.hiv_names)
        result = self.test_hiv.merge_breakpoints()
        self.assertEqual(expected, result)

        expected = [('07_BC', ('D', 'N'), 100, 360, 4.881427058185601e-48),
                    ('07_BC', ('D', 'N'), 1280, 1340, 1.0595039034315693e-11),
                    ('07_BC', ('D', 'N'), 2060, 2580, 2.382833012446353e-95),
                    ('07_BC', ('D', 'N'), 2980, 3200, 9.273248781882314e-41),
                    ('07_BC', ('D', 'N'), 5800, 5940, 3.3446094314319793e-26),
                    ('07_BC', ('D', 'N'), 6140, 6400, 4.88142705819169e-48),
                    ('07_BC', ('D', 'N'), 8960, 9140, 1.761639412097178e-33),
                    ('07_BC', ('D', 'S'), 100, 560, 1.448296664276467e-83),
                    ('07_BC', ('D', 'S'), 580, 1260, 3.484473924295404e-123),
                    ('07_BC', ('D', 'S'), 2620, 2900, 3.7685488477306826e-51),
                    ('07_BC', ('N', 'S'), 100, 560, 1.365044950915786e-86),
                    ('07_BC', ('N', 'S'), 580, 640, 6.311292745859333e-12),
                    ('07_BC', ('N', 'S'), 680, 840, 1.3613129217492749e-30),
                    ('07_BC', ('N', 'S'), 1200, 1380, 2.5159261286591272e-34),
                    ('07_BC', ('N', 'S'), 2840, 2880, 2.416580662256423e-08),
                    ('07_BC', ('N', 'S'), 3000, 3180, 2.515989943996236e-34),
                    ('07_BC', ('N', 'S'), 3300, 3560, 2.935565371247222e-49),
                    ('07_BC', ('N', 'S'), 3720, 3740, 0.00015364197645138482),
                    ('07_BC', ('N', 'S'), 4820, 5000, 2.515993974743578e-34),
                    ('07_BC', ('N', 'S'), 5240, 5360, 3.9854030208559366e-23),
                    ('07_BC', ('N', 'S'), 5480, 5500, 0.00018283510051390274),
                    ('07_BC', ('N', 'S'), 7400, 7440, 3.4157026740246995e-08),
                    ('07_BC', ('N', 'S'), 7820, 7840, 0.0001283538103179031),
                    ('07_BC', ('N', 'S'), 7960, 8420, 1.3650449509157897e-86),
                    ('07_BC', ('N', 'S'), 9260, 9280, 0.00018425967304589247),
                    ('07_BC', ('N', 'S'), 9340, 9380, 3.4157882054062814e-08),
                    ('B', ('D', 'N'), 100, 560, 6.722257329590765e-90),
                    ('B', ('D', 'N'), 1280, 1340, 2.206148653469563e-12),
                    ('B', ('D', 'N'), 2060, 2560, 1.1840935474819853e-97),
                    ('B', ('D', 'N'), 3000, 3160, 9.626833968794724e-32),
                    ('B', ('D', 'N'), 8960, 9140, 1.2776707863712195e-35),
                    ('B', ('D', 'S'), 2620, 2940, 2.6131494126158124e-62),
                    ('B', ('D', 'S'), 6060, 6080, 0.0001297342192715822),
                    ('C', ('D', 'N'), 100, 360, 2.918988522279219e-52),
                    ('C', ('D', 'N'), 1280, 1340, 1.0870034646348963e-12),
                    ('C', ('D', 'N'), 2060, 2580, 8.520493993200118e-104),
                    ('C', ('D', 'N'), 2980, 3200, 2.4754728855475345e-44),
                    ('C', ('D', 'N'), 5800, 5940, 1.7797315605509396e-28),
                    ('C', ('D', 'N'), 6140, 6400, 2.918988522274107e-52),
                    ('C', ('D', 'N'), 8960, 9140, 2.0993457026384212e-36),
                    ('C', ('D', 'S'), 100, 560, 7.405867461585028e-88),
                    ('C', ('D', 'S'), 580, 1260, 1.579441542995601e-129),
                    ('C', ('D', 'S'), 2620, 2900, 9.206397024627418e-54),
                    ('C', ('N', 'S'), 100, 560, 7.248090857184098e-96),
                    ('C', ('N', 'S'), 580, 640, 3.894609793204167e-13),
                    ('C', ('N', 'S'), 680, 840, 8.08899389659263e-34),
                    ('C', ('N', 'S'), 1200, 1380, 5.907301042633563e-38),
                    ('C', ('N', 'S'), 2840, 2880, 4.110001158688038e-09),
                    ('C', ('N', 'S'), 3000, 3180, 5.907197658951816e-38),
                    ('C', ('N', 'S'), 3300, 3560, 1.6801334934754224e-54),
                    ('C', ('N', 'S'), 3720, 3740, 5.898112808213309e-05),
                    ('C', ('N', 'S'), 4820, 5000, 5.907319782295677e-38),
                    ('C', ('N', 'S'), 5240, 5360, 1.5167055433906316e-25),
                    ('C', ('N', 'S'), 5480, 5500, 7.203848528220467e-05),
                    ('C', ('N', 'S'), 7400, 7440, 5.332823098813123e-09),
                    ('C', ('N', 'S'), 7820, 7840, 4.8372974027197606e-05),
                    ('C', ('N', 'S'), 7960, 8420, 7.248090857184098e-96),
                    ('C', ('N', 'S'), 9260, 9280, 7.274132246895203e-05),
                    ('C', ('N', 'S'), 9340, 9380, 5.333051421533079e-09),
                    ('D', ('07_BC', 'N'), 2620, 2940, 5.987754444592552e-61),
                    ('D', ('07_BC', 'N'), 6060, 6080, 0.00016663021506877894),
                    ('D', ('07_BC', 'S'), 580, 1260, 2.4953954883393473e-127),
                    ('D', ('07_BC', 'S'), 1280, 1340, 6.736813087332313e-12),
                    ('D', ('07_BC', 'S'), 2060, 2580, 1.534785775305691e-97),
                    ('D', ('07_BC', 'S'), 2620, 2900, 7.403044078511966e-53),
                    ('D', ('07_BC', 'S'), 3000, 3180, 3.072401990173266e-34),
                    ('D', ('07_BC', 'S'), 8960, 9140, 3.072402236920229e-34),
                    ('D', ('B', 'N'), 580, 1260, 1.7967628822515222e-130),
                    ('D', ('B', 'N'), 1280, 1340, 3.5606555379764214e-12),
                    ('D', ('B', 'N'), 2060, 2560, 3.970912067404395e-96),
                    ('D', ('B', 'N'), 2620, 2900, 3.76165144442256e-54),
                    ('D', ('B', 'N'), 3000, 3160, 2.962408059129308e-31),
                    ('D', ('B', 'N'), 8960, 9140, 4.524811559040887e-35),
                    ('D', ('B', 'S'), 1280, 1340, 4.172250902051943e-12),
                    ('D', ('B', 'S'), 2060, 2560, 1.4819449290358866e-95),
                    ('D', ('B', 'S'), 2980, 3180, 1.1703959859207535e-38),
                    ('D', ('B', 'S'), 5800, 5940, 2.8043548047407118e-27),
                    ('D', ('B', 'S'), 6140, 6400, 4.884641262008904e-50),
                    ('D', ('B', 'S'), 8960, 9140, 7.269416347802067e-35),
                    ('D', ('C', 'N'), 2620, 2940, 3.10216544900611e-61),
                    ('D', ('C', 'N'), 6060, 6080, 0.00015986980274114264),
                    ('D', ('C', 'S'), 580, 1260, 7.267844587748626e-129),
                    ('D', ('C', 'S'), 1280, 1340, 4.920227470382949e-12),
                    ('D', ('C', 'S'), 2060, 2580, 1.02721920481641e-98),
                    ('D', ('C', 'S'), 2620, 2900, 1.72603412531948e-53),
                    ('D', ('C', 'S'), 3000, 3180, 1.204925825875451e-34),
                    ('D', ('C', 'S'), 8960, 9140, 1.2049258296788594e-34),
                    ('N', ('07_BC', 'D'), 100, 360, 3.954093584452123e-52),
                    ('N', ('07_BC', 'D'), 1280, 1340, 1.3731514512506501e-12),
                    ('N', ('07_BC', 'D'), 2060, 2580, 1.5634856154806153e-103),
                    ('N', ('07_BC', 'D'), 2620, 2940, 5.430339837751515e-64),
                    ('N', ('07_BC', 'D'), 2980, 3200, 3.200321202732338e-44),
                    ('N', ('07_BC', 'D'), 5800, 5940, 2.0964607554535662e-28),
                    ('N', ('07_BC', 'D'), 6060, 6080, 9.095025851436264e-05),
                    ('N', ('07_BC', 'D'), 6140, 6400, 3.954093594577339e-52),
                    ('N', ('07_BC', 'D'), 8960, 9140, 2.5902411149014285e-36),
                    ('N', ('07_BC', 'S'), 580, 640, 1.438670875187207e-12),
                    ('N', ('07_BC', 'S'), 1200, 1380, 2.978978535620813e-36),
                    ('N', ('07_BC', 'S'), 3000, 3180, 2.9789760307758125e-36),
                    ('N', ('07_BC', 'S'), 3720, 3740, 0.00011250132041173),
                    ('N', ('07_BC', 'S'), 5480, 5500, 0.00011250132041173),
                    ('N', ('07_BC', 'S'), 7400, 7440, 1.2499981322585755e-08),
                    ('N', ('07_BC', 'S'), 7820, 7840, 9.241898245445662e-05),
                    ('N', ('07_BC', 'S'), 9340, 9380, 1.2744815443958386e-08),
                    ('N', ('B', 'D'), 100, 560, 1.527647591619792e-94),
                    ('N', ('B', 'D'), 580, 1260, 2.067819024048125e-139),
                    ('N', ('B', 'D'), 2620, 2900, 7.845539261205802e-58),
                    ('N', ('C', 'D'), 100, 360, 1.918582903612156e-55),
                    ('N', ('C', 'D'), 1280, 1340, 2.360280087690839e-13),
                    ('N', ('C', 'D'), 2060, 2580, 3.680960358033879e-110),
                    ('N', ('C', 'D'), 2620, 2940, 4.528726421739898e-68),
                    ('N', ('C', 'D'), 2980, 3200, 5.023235731684667e-47),
                    ('N', ('C', 'D'), 5800, 5940, 3.4434150250666295e-30),
                    ('N', ('C', 'D'), 6060, 6080, 5.583453128943263e-05),
                    ('N', ('C', 'D'), 6140, 6400, 1.9185829036121357e-55),
                    ('N', ('C', 'D'), 8960, 9140, 1.3151840963243797e-38),
                    ('N', ('C', 'S'), 580, 640, 2.3604538629011573e-13),
                    ('N', ('C', 'S'), 1200, 1380, 1.3151840947381604e-38),
                    ('N', ('C', 'S'), 3000, 3180, 1.3151835237701565e-38),
                    ('N', ('C', 'S'), 3720, 3740, 6.154828049278473e-05),
                    ('N', ('C', 'S'), 5480, 5500, 6.154828049278473e-05),
                    ('N', ('C', 'S'), 7400, 7440, 3.730560757117377e-09),
                    ('N', ('C', 'S'), 7820, 7840, 4.9658051093140255e-05),
                    ('N', ('C', 'S'), 9340, 9380, 3.819235931181965e-09),
                    ('S', ('07_BC', 'D'), 100, 560, 1.0107848486621713e-88),
                    ('S', ('07_BC', 'D'), 1280, 1340, 3.329208729986169e-12),
                    ('S', ('07_BC', 'D'), 2060, 2580, 3.3651584141568535e-100),
                    ('S', ('07_BC', 'D'), 3000, 3180, 3.690100802328052e-35),
                    ('S', ('07_BC', 'D'), 8960, 9140, 3.690049439406398e-35),
                    ('S', ('07_BC', 'N'), 100, 560, 1.6168771209832646e-89),
                    ('S', ('07_BC', 'N'), 680, 840, 1.306255411297981e-31),
                    ('S', ('07_BC', 'N'), 2840, 2880, 1.5123935899049965e-08),
                    ('S', ('07_BC', 'N'), 3300, 3560, 6.510285764706049e-51),
                    ('S', ('07_BC', 'N'), 4820, 5000, 1.7999889072375998e-35),
                    ('S', ('07_BC', 'N'), 5240, 5360, 6.871414984475036e-24),
                    ('S', ('07_BC', 'N'), 7960, 8420, 1.6168771209832646e-89),
                    ('S', ('07_BC', 'N'), 9260, 9280, 0.0001374278243352618),
                    ('S', ('B', 'D'), 1280, 1340, 1.325918664612657e-12),
                    ('S', ('B', 'D'), 2060, 2560, 1.0507307736000933e-99),
                    ('S', ('B', 'D'), 2620, 2940, 4.505619671022191e-64),
                    ('S', ('B', 'D'), 2980, 3180, 2.5621028370654944e-40),
                    ('S', ('B', 'D'), 5800, 5940, 1.9320466407434272e-28),
                    ('S', ('B', 'D'), 6060, 6080, 8.98626878848598e-05),
                    ('S', ('B', 'D'), 6140, 6400, 3.397625871983414e-52),
                    ('S', ('B', 'D'), 8960, 9140, 2.3320111160042457e-36),
                    ('S', ('C', 'D'), 100, 560, 3.320582793279103e-91),
                    ('S', ('C', 'D'), 1280, 1340, 1.579127227806961e-12),
                    ('S', ('C', 'D'), 2060, 2580, 5.243648495653591e-103),
                    ('S', ('C', 'D'), 3000, 3180, 3.937784248901982e-36),
                    ('S', ('C', 'D'), 8960, 9140, 3.9325198503587943e-36),
                    ('S', ('C', 'N'), 100, 560, 1.6274419956304624e-91),
                    ('S', ('C', 'N'), 680, 840, 2.638466169416362e-32),
                    ('S', ('C', 'N'), 2840, 2880, 1.0038822440963653e-08),
                    ('S', ('C', 'N'), 3300, 3560, 4.839093218622394e-52),
                    ('S', ('C', 'N'), 4820, 5000, 2.9722788407446952e-36),
                    ('S', ('C', 'N'), 5240, 5360, 2.0703554611687758e-24),
                    ('S', ('C', 'N'), 7960, 8420, 1.6274419956304624e-91),
                    ('S', ('C', 'N'), 9260, 9280, 0.00011250132041173)]
        total_num_trps = sum(1 for _ in product(
                list(combinations(range(self.hiv_align.shape[0]), 1)),
                list(combinations(range(self.hiv_align_r.shape[0]), 2))
            ))
        self.test_hiv_r.execute_all(total_combinations=total_num_trps, seq_names=self.hiv_names, ref_names=self.hiv_names_r)
        result = self.test_hiv_r.merge_breakpoints()
        self.assertEqual(expected, result)

if __name__ == '__main__':
    unittest.main()
