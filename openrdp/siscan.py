import random
import json
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
from scipy.stats import norm

from .common import identify_recombinant


class Siscan:
    def __init__(self, align, win_size=150, step_size=100, strip_gaps=True, pvalue_perm_num=1100,
                 scan_perm_num=100, random_seed=3, max_pvalue=0.05, settings=None, ref_align=None, verbose=False):
        """
        Constructs a Siscan object
        :param win_size: the size of the sliding window
        :param step_size: the step size of the sliding window
        :param strip_gaps: Strip gaps or keep gaps
        :param pvalue_perm_num: p-value permutation number
        :param scan_perm_num: number of permutations of the scans
        :param random_seed: the random seed
        """
        self.align = align
        if settings:
            self.set_options_from_config(settings)
            self.validate_options(align)
        else: # pragma: no cover
            self.win_size = win_size
            self.step_size = step_size
            self.strip_gaps = strip_gaps
            self.pvalue_perm_num = pvalue_perm_num
            self.scan_perm_num = scan_perm_num
            self.random_seed = random_seed
            self.max_pvalues = max_pvalue

        self.raw_results = []
        self.results = []
        self.name = 'siscan'

    def set_options_from_config(self, settings):
        """
        Set the parameters of Siscan from the  config file
        :param settings: a dictionary of settings
        """
        self.win_size = int(settings['win_size'])
        self.step_size = int(settings['step_size'])
        self.max_pvalues = abs(float(settings['max_pvalue']))

        if settings['strip_gaps'] == 'True':
            self.strip_gaps = True
        else:
            self.strip_gaps = False

        self.pvalue_perm_num = int(settings['pvalue_perm_num'])
        self.scan_perm_num = int(settings['scan_perm_num'])
        self.random_seed = int(settings['random_seed'])

    def validate_options(self, align):
        """
        Check if the options from the config file are valid
        If the options are invalid, the default value will be used instead
        """
        if self.win_size < 0 or self.win_size > align.shape[1]:
            print("Invalid option for 'win_size'.\nUsing default value (200) instead.")
            self.win_size = 200

        if self.step_size < 0 or self.step_size >= self.win_size or self.step_size > align.shape[1]:
            print("Invalid option for 'step_size'.\nUsing default value (20) instead.")
            self.step_size = 20

        if self.pvalue_perm_num <= 0:
            print("Invalid option for 'pvalue_perm_num'.\nUsing default value (1100) instead.")
            self.pvalue_perm_num = 1100

        if self.scan_perm_num <= 0:
            print("Invalid option for 'scan_perm_num'.\nUsing default value (100) instead.")
            self.scan_perm_num = 100

        if self.random_seed <= 0:
            print("Invalid option for 'random_seed'.\nUsing default value (3) instead.")
            self.random_seed = 3

    @staticmethod
    def count_patterns(seq_array):
        a, b, c, d = seq_array

        pat_counts = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

        ab = a == b
        bc = b == c
        ac = a == c
        ad = a == d
        cd = c == d
        bd = b == d

        # 1 ~ 2 ~ 3 ~ 4
        pat_counts[0] = np.sum(np.logical_not(ab) & np.logical_not(ac) & np.logical_not(ad) & np.logical_not(bc) & np.logical_not(bd) & np.logical_not(cd))
        # 1 = 2 ~ 3 ~ 4
        pat_counts[1] = np.sum(ab & np.logical_not(ac) & np.logical_not(ad) & np.logical_not(cd))
        # 1 = 3 ~ 2 ~ 4
        pat_counts[2] = np.sum(ac & np.logical_not(ab) & np.logical_not(ad) & np.logical_not (bd))
        # 1 = 4 ~ 2 ~ 3
        pat_counts[3] = np.sum(ad & np.logical_not(ab) & np.logical_not(ac) & np.logical_not (bc))
        # 2 = 3 ~ 1 ~ 4
        pat_counts[4] = np.sum(bc & np.logical_not(ab) & np.logical_not(bd) & np.logical_not(ad))
        # 2 = 4 ~ 1 ~ 3
        pat_counts[5] = np.sum(bd & np.logical_not(ab) & np.logical_not(bc) & np.logical_not (ac))
        # 3 = 4 ~ 1 ~ 2
        pat_counts[6] = np.sum(cd & np.logical_not(bc) & np.logical_not(ac) & np.logical_not(ab))
        # 1 = 2 ~ 3 = 4
        pat_counts[7] = np.sum(ab & cd & np.logical_not(bc))
        # 1 = 3 ~ 2 = 4
        pat_counts[8] = np.sum(ac & bd & np.logical_not(bc))
        # 1 = 4 ~ 2 = 3
        pat_counts[9] = np.sum(ad & bc & np.logical_not(ab))
        # 1 = 2 = 3 ~ 4
        pat_counts[10] = np.sum(ab & bc & np.logical_not(ad))
        # 1 = 2 = 4 ~ 3
        pat_counts[11] = np.sum(ab & bd & np.logical_not(ac))
        # 1 = 3 = 4 ~ 2
        pat_counts[12] = np.sum(ac & cd & np.logical_not(ab))
        # 2 = 3 = 4 ~ 1
        pat_counts[13] = np.sum(bc & cd & np.logical_not(ab))
        # 1 = 2 = 3 = 4
        pat_counts[14] = np.sum(ab & bc & ad)

        print(pat_counts)

        return pat_counts

    @staticmethod
    def sum_pattern_counts(pat_counts):
        """
        Sum counts where 2 sequences are identical
        :param pat_counts: vector of counts corresponding to each pattern
        :return: vector of summed counts
        """
        sum_pat_counts = [0, 0, 0, 0, 0, 0, 0, 0, 0]

        sum_pat_counts[0] = pat_counts[1] + pat_counts[6] + pat_counts[7]  # 2 + 7 + 8
        sum_pat_counts[1] = pat_counts[2] + pat_counts[5] + pat_counts[8]  # 3 + 6 + 9
        sum_pat_counts[2] = pat_counts[3] + pat_counts[4] + pat_counts[9]  # 4 + 5 + 10
        sum_pat_counts[3] = pat_counts[1] + pat_counts[7] + pat_counts[10] + pat_counts[11]  # 2 + 8 + 11 + 12
        sum_pat_counts[4] = pat_counts[2] + pat_counts[8] + pat_counts[10] + pat_counts[12]  # 3 + 9 + 11 + 13
        sum_pat_counts[5] = pat_counts[3] + pat_counts[9] + pat_counts[11] + pat_counts[12]  # 4 + 10 + 12 + 13
        sum_pat_counts[6] = pat_counts[4] + pat_counts[9] + pat_counts[10] + pat_counts[13]  # 5 + 10 + 11 + 14
        sum_pat_counts[7] = pat_counts[5] + pat_counts[8] + pat_counts[11] + pat_counts[13]  # 6 + 9 + 12 + 14
        sum_pat_counts[8] = pat_counts[6] + pat_counts[7] + pat_counts[12] + pat_counts[13]  # 7 + 8 + 13 + 14

        return sum_pat_counts

    def find_interval(major, minor, ind):
        """
        major zlist
        minor zlist
        ind of list to start looking
        """
        z_cut = norm.ppf(1 - 0.05/len(major)) # zscore corrected to number of windows
        end = int
        while minor[end] > major[end] and minor[end] > z_cut:
            end += 1
        return int, end

    def find_signal(self, major, minor1, minor2):
        """
        major_comb is pattern/sum set of two most related
        minor_combs is least related
        """
        min1_intervals, min2_intervals = [], []
        z_cut = norm.ppf(1 - 0.05/len(major))
        for ind, maj, min1, min2 in enumerate(zip(major, minor1, minor2)):
            if min1 > maj and min1 > z_cut:
                min1_intervals.append(self.find_interval(maj, min1, ind))
            if min2 > maj and min2 > z_cut:
                min2_intervals.append(self.find_interval(maj, min1, ind))
        return min1_intervals, min2_intervals

    def fine_adj(major, minor1, minor2, interval):
        """
        
        """
        
        
        
        

    def execute(self, triplet):
        """
        Do Sister-scanning as described in Gibbs, Armstrong, and Gibbs (2000), using a randomized 4th sequence
        :param triplet: a triplet object
        """
        random.seed(self.random_seed)
        np.random.seed(self.random_seed)

        # Initialize list to map z_values to window positions
        z_pattern_values = [[None for i in range(14)] for i in range(triplet.sequences.shape[1])]
        z_sum_values = [[None for i in range(9)] for i in range(triplet.sequences.shape[1])] 


        # Based on leading edge of the window
        for window in range(0, self.align.shape[1], self.step_size):
            win_end = window + self.win_size

            # shouldn't run if it doesn't get plotted
            if (window + win_end)//2 >= len(z_pattern_values):
                continue

            # Label sequences
            a = triplet.sequences[0][window: win_end]
            b = triplet.sequences[1][window: win_end]
            c = triplet.sequences[2][window: win_end]

            # Create the fourth sequence through horizontal randomization
            selected_seq = random.choice((0, 1, 2))
            d = triplet.sequences[selected_seq][window: win_end]
            np.random.shuffle(d)

            # (1) Count number of positions within a window that conform to each pattern
            seq_array = np.array([a, b, c, d])
            pat_counts = self.count_patterns(seq_array)

            # (2) Sum counts where 2 sequences are identical 
            # RDP5 does 2 and three together along with permutations
            sum_pat_counts = self.sum_pattern_counts(pat_counts)

            # (4) Create 4 vertically randomized sequences (steps 1 and 2 + 3 [as per RDP5]), repeat for 100 times
            p_counts = []
            sum_p_counts = []
            for i in range(self.scan_perm_num):

                # Generate 4 vertically randomized sequences (shuffle the values in the columns)
                a1 = np.apply_along_axis(np.random.permutation, 0, seq_array)

                # Count number of patterns and sum counts
                curr_p = self.count_patterns(a1)
                p_counts.append(curr_p)
                curr_s = self.sum_pattern_counts(curr_p)
                sum_p_counts.append(curr_s)

            # (5) Calculate Z-scores for each pattern and sum of patterns for each window, position should be middle of window            
            pop_mean_pcounts = np.mean(p_counts, axis=0) # array of [mean pat 1, mean pat 2, ...]
            pop_mean_patsum = np.mean(sum_p_counts, axis=0)
            pop_std_pcounts = np.std(p_counts, axis=0)
            pop_std_patsum = np.std(sum_p_counts, axis=0)

            z_pat_counts = [0 for i in pat_counts]
            for num, value in enumerate(pat_counts): # if none it will give error
                z_pat_counts[num] = (value - pop_mean_pcounts[num]) / pop_std_pcounts[num]
                
            z_sum_counts = [0 for i in sum_pat_counts]
            for num, value in enumerate(sum_pat_counts):
                z_sum_counts[num] = (value - pop_mean_patsum[num]) / pop_std_patsum[num]

            # update to total by middle of window
            z_pattern_values[(window + win_end)//2] = z_pat_counts
            z_sum_values[(window + win_end)//2] = z_sum_counts

        seq_1_2 = [z_pat_counts[3], z_pat_counts[9], z_sum_counts[2], z_sum_counts[5]] # p2/8, s1/4, 1==2
        seq_1_3 = [z_pat_counts[4], z_pat_counts[10], z_sum_counts[6], z_sum_counts[7]] # p3/9, s5/6, 1==3
        seq_2_3 = [z_pat_counts[6], z_pat_counts[11], z_sum_counts[4], z_sum_counts[8]] # p5/10, s3/7, 2==3

        # find major combination
        major = None

        #
        if major == 1: # seq1/2 major 
            min1, min2 = self.find_signal(seq_1_2, seq_1_3, seq_2_3)
        if major == 2: # seq2/3 major
            min1, min2 = self.find_signal(seq_2_3, seq_1_2, seq_1_3)
        if major == 3: # seq1/3 major
            min1, min2 = self.find_signal(seq_1_3, seq_1_2, seq_2_3)

        for ind, start, end in enumerate(min2):
            nt_start = start * self.step_size
            nt_end = end * self.step_size + self.win_size + 1
            
            # adjust the start of the window til it matches the pattern
            for ind, seq1, seq2, seq3 in enumerate(zip(a[nt_start:nt_end], b[nt_start:nt_end], c[nt_start:nt_end])):
                if major == 1: # seq1/2 major
                    if seq1 == seq2 and seq1 != seq3:
                        nt_start = nt_start + ind
                if major == 2:
                    if seq1 == seq3 and seq1 != seq2:
                        nt_start = nt_start + ind
                if major == 3:
                    if seq2 == seq3 and seq1 != seq2:
                        nt_start = nt_start + ind
            min2[ind][0] = nt_start # update new index based on NT
            
            nt_start = start * self.step_size - 1
            nt_end = end * self.step_size + self.win_size
            # adjust the end of the interval til it matches the pattern
            for ind, seq1, seq2, seq3 in enumerate(zip(a[nt_end:nt_start:-1], b[nt_end:nt_start:-1], c[nt_end:nt_start:-1])):
                if major == 1: # seq1/2 major
                    if seq1 == seq2 and seq1 != seq3:
                        nt_start = nt_start + ind
                if major == 2:
                    if seq1 == seq3 and seq1 != seq2:
                        nt_start = nt_start + ind
                if major == 3:
                    if seq2 == seq3 and seq1 != seq2:
                        nt_start = nt_start + ind
            min2[ind][1] = nt_start # update new end index based on NT

        for ind, start, end in enumerate(min1):
            nt_start = start * self.step_size
            nt_end = end * self.step_size + self.win_size + 1
            
            # adjust the start of the window til it matches the pattern
            for ind, seq1, seq2, seq3 in enumerate(zip(a[nt_start:nt_end], b[nt_start:nt_end], c[nt_start:nt_end])):
                if major == 1: # seq1/2 major
                    if seq1 == seq2 and seq1 != seq3:
                        nt_start = nt_start + ind
                if major == 2:
                    if seq1 == seq3 and seq1 != seq2:
                        nt_start = nt_start + ind
                if major == 3:
                    if seq2 == seq3 and seq1 != seq2:
                        nt_start = nt_start + ind
            min1[ind][0] = nt_start # update new index based on NT
            
            nt_start = start * self.step_size - 1
            nt_end = end * self.step_size + self.win_size
            # adjust the end of the interval til it matches the pattern
            for ind, seq1, seq2, seq3 in enumerate(zip(a[nt_end:nt_start:-1], b[nt_end:nt_start:-1], c[nt_end:nt_start:-1])):
                if major == 1: # seq1/2 major
                    if seq1 == seq2 and seq1 != seq3:
                        nt_start = nt_start + ind
                if major == 2:
                    if seq1 == seq3 and seq1 != seq2:
                        nt_start = nt_start + ind
                if major == 3:
                    if seq2 == seq3 and seq1 != seq2:
                        nt_start = nt_start + ind
            min1[ind][1] = nt_start # update new end index based on NT

        if major == 1: # The recombinant should be the sequence that's the least releated
            rec_name, parents = triplet.name[0], (triplet.name[1], triplet.name[2])
        if major == 2:
            rec_name, parents = triplet.name[1], (triplet.name[0], triplet.name[2])
        if major == 3:
            rec_name, parents = triplet.name[2], (triplet.name[0], triplet.name[1])

        for start, end in min1:
            self.raw_results.append((rec_name, parents, (start, end), None)) # TODO need to add Z score
        for start, end in min2:
            self.raw_results.append((rec_name, parents, (start, end), None))

        # peaks = find_peaks(sum_pat_zscore, distance=self.win_size)
        # for k, peak in enumerate(peaks[0]):
        #     search_win_size = 1
        #     while peak - search_win_size > 0 \
        #             and peak + search_win_size < len(sum_pat_zscore) - 1 \
        #             and sum_pat_zscore[peak + search_win_size] > 0.3 * sum_pat_zscore[peak] \
        #             and sum_pat_zscore[peak - search_win_size] > 0.3 * sum_pat_zscore[peak]:
        #         search_win_size += 1

        #     if sum_pat_zscore[peak + search_win_size] > sum_pat_zscore[peak - search_win_size]:
        #         aln_pos = (int(peak), int(peak + search_win_size + self.win_size))
        #     else:
        #         aln_pos = (int(peak - search_win_size), int(peak + self.win_size))

        #     rec_name, parents = identify_recombinant(triplet, aln_pos)
        #     if (rec_name, parents, *aln_pos, abs(sum_pat_zscore[peak])) not in self.raw_results:
        #         self.raw_results.append((rec_name, parents, *aln_pos, abs(sum_pat_zscore[peak])))
