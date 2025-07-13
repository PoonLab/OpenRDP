import random

import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
from scipy.stats import norm

from .common import identify_recombinant, find_parent


class Siscan:
    def __init__(self, align, win_size=200, step_size=20, strip_gaps=True, pvalue_perm_num=1000,
                 scan_perm_num=100, random_seed=3, max_pvalue=0.05, settings=None, 
                 ref_align=None, verbose=False):
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

    def find_interval(self, ind, maj, close, ps):
        """
        close, major/major zlist
        maj, major/minor zlist
        ind of list to start looking
        ps int, ind of which p/s on the array of maj is the one we care about
        """
        z_cut = norm.ppf(1-(0.05/(len(close)//self.win_size)))  # zscore corrected to number of windows
        end = ind
        while True: 
            if maj[end][ps] and close[end][ps]:
                if maj[end][ps] > close[end][ps] and maj[end][ps] > z_cut:
                    end += 1
                else:
                    break
            else: # it's none so move forward a NT
                end += 1
            if end == len(maj)-1:
                break
        return ind, end

    def find_signal(self, close, maj1, maj2, signal):
        """
        close is the list of p/s of the two most related
        maj1/maj2 is the list of pattern/sum set of related to non related
        signal, index used in see which set is maj/min

        return list of indexes that matter, ind is the pattern/sum that matters
        """
        # use f1 and f2 to keep the intervals where it is found to be significant
        m1_intervals, m2_intervals = [], []
        z_cut = norm.ppf(1-(0.05/(self.align.shape[1]//self.step_size))) 
        # check per pattern/sum
        for ps in range(4):

            ind, ind2 = 0, 0
            while ind < len(close):
                
                value = maj1[ind][ps]

                # if any pattern for the minor and major is greater than major major
                if value:
                    if value > close[ind][ps] and value > z_cut:
                        # in this case `i` is the index of the pattern of importance
                        ind, new = self.find_interval(ind, maj1, close, ps)
                        ps_fin = self.get_ps(signal, 0, ps) # identifier for later use to optimize shuffling (ok, upon further checking it's super slow still here)
                        m1_intervals.append((ind, new, signal, ps_fin)) # which pattern/sum, and starting index and ending index of the recombinant region
                        ind = new # so we don't keep adding them back and forth
                ind += 1

            while ind2 < len(close):
                # do again with maj2
                value = maj2[ind2][ps]
                if value:
                    if value > close[ind2][ps] and value > z_cut:
                        ind, new = self.find_interval(ind2, maj2, close, ps)
                        ps_fin = self.get_ps(signal, 1, ps)
                        m2_intervals.append((ind2, new, signal, ps_fin))
                        ind2 = new
                ind2 += 1

        return m1_intervals, m2_intervals        

    def adjust_nt_f(self, eq1, eq2, diff, start, end):
        """
        adjust the range of nucleotide positions to match the significant patterns for start
        
        eq1/eq2 are the triplet.sequences that we want to equal one another
        diff is the sequence that we want to be different
        start int, is the start index of the nucleotide sequences
        end int, is the end index of window 
        """
        while start < end:
            if eq1[start] == eq2[start] and eq1[start] != diff[start]:
                return start
            start += 1
        return start # assumes that it's on the same window

    def adjust_nt_r(self, eq1, eq2, diff, start, end):
        """
        same as above, just other way
        """
        while start > end:
            if eq1[start] == eq2[start] and eq1[start] != diff[start]:
                return start
            start -= 1
        return start 

    def shuffle(self, parent1, parent2, recombinant, out, ps_ind):
        """
        recalculate the value of pattern and sum counts for the window by regenerating a null distribution

        parent1, parent2, strings of the window for the sequence that are most related
        recombinant, string of the window of the recombinant
        out, string of outlier sequence
        ps, (true/false, index), index of the pattern/sum we care about. T = pattern, F = sum
        ind, int, second index is which one that matters
        
        return zscore for the final p-value
        """

        #TODO find a way to not calculate all numbers and just find the one you care about
        # ps_ind doesn't work the way you want it to right now
        seqs = [parent1, parent2, recombinant, out]
        seqs = list(map(np.array, seqs))
        ps, ind = ps_ind

        p = self.count_patterns(seqs)
        if ps:
            true_val = p[ind]
        else:
            s = self.sum_pattern_counts(p)
            true_val = s[ind]

        null_dist = []
        for i in range(self.scan_perm_num):
            # Generate 4 vertically randomized sequences (shuffle the values in the columns)
            a1 = np.apply_along_axis(np.random.permutation, 0, seqs)

            curr_p = self.count_patterns(a1)    
            # Count number of patterns and sum counts
            if ps: # it's a pattern we care about
                null_dist.append(curr_p[ind])

            else: # it's a sum we care about
                curr_s = self.sum_pattern_counts(curr_p)
                null_dist.append(curr_s[ind])


        mean = np.mean(null_dist)
        sd = np.std(null_dist, axis=0)

        return (true_val - mean)/sd


    def get_ps(self, signal, ind, lr):
        """
        get a tuple that works as an identifier for if it is a sum or a pattern, and which one

        signal, int, denotes which seq_X_Y was used in self.execute() in the signal_orders list
        lr, int, determines which min_maj combo was used in self.find_signal() 0 should be min1, 1 should be min2
        ind, int, index of position the pattern that denotes how X = Y ~ Z (so which position in seq_X_Y)
        
        return: (True/False denoting Pattern/Sum, index of the position that matter)
        """
        
        seq_1_2 = [(True, 2), (True, 7), (False, 0), (False, 3)]
        seq_1_3 = [(True, 2), (True, 8), (False, 4), (False, 5)]
        seq_2_3 = [(True, 4), (True, 9), (False, 2), (False, 6)]

        signal_orders = [
            (seq_1_2, seq_1_3),  # min == 0 2==3
            (seq_2_3, seq_1_2),  # min == 1 1==3
            (seq_1_3, seq_2_3)   # min == 2 1==2
        ]
        return signal_orders[signal][ind][lr]


    def execute(self, triplet, tree):
        """
        Do Sister-scanning as described in Gibbs, Armstrong, and Gibbs (2000), using a randomized 4th sequence
        updated with notes from Darren on how RDP4/5 works
        :param triplet: a triplet object
        :param tree: node object from upgma
        """
        random.seed(self.random_seed)
        np.random.seed(self.random_seed)

        # Initialize list to map z_values to window positions
        z_pattern_values = [[None for i in range(14)] for j in range(triplet.sequences.shape[1])]
        z_sum_values = [[None for i in range(9)] for j in range(triplet.sequences.shape[1])] 

        # Create the fourth sequence through horizontal randomization
        selected_seq = random.choice((0, 1, 2))
        out = triplet.sequences[selected_seq]
        np.random.shuffle(out)

        # Based on leading edge of the window
        for window in range(0, self.align.shape[1], self.step_size):
            win_end = window + self.win_size
            
            # shouldn't run if it doesn't get plotted
            # TODO check if this is the right way to handle it and not to just put at end
            if (window + win_end)//2 >= len(z_pattern_values):
                continue

            # Label sequences
            a = triplet.sequences[0][window: win_end]
            b = triplet.sequences[1][window: win_end]
            c = triplet.sequences[2][window: win_end]
            d = out[window: win_end]

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
                if pop_std_pcounts[num] == 0:
                    continue
                z_pat_counts[num] = (value - pop_mean_pcounts[num]) / pop_std_pcounts[num]

            z_sum_counts = [0 for i in sum_pat_counts]
            for num, value in enumerate(sum_pat_counts):
                if pop_std_patsum[num] == 0:
                    continue
                z_sum_counts[num] = (value - pop_mean_patsum[num]) / pop_std_patsum[num]

            # update to total by middle of window
            z_pattern_values[(window + win_end)//2] = z_pat_counts
            z_sum_values[(window + win_end)//2] = z_sum_counts

        # extract the patterns and sums according to which patterns are describing each pattern 
        seq_1_2, seq_1_3, seq_2_3 = [], [], []
        for i in range(len(z_sum_values)):
            seq_1_2.append([z_pattern_values[i][2], z_pattern_values[i][7], z_sum_values[i][0], z_sum_values[i][3]]) # p2/8, s1/4, 1==2
            seq_1_3.append([z_pattern_values[i][2], z_pattern_values[i][8], z_sum_values[i][4], z_sum_values[i][5]]) # p3/9, s5/6, 1==3
            seq_2_3.append([z_pattern_values[i][4], z_pattern_values[i][9], z_sum_values[i][2], z_sum_values[i][6]]) # p5/10, s3/7, 2==3

        # DEBUGGING
        # with open('siscan.json', 'w') as file:
        #     json.dump([seq_1_2, seq_2_3, seq_1_3], file)

        # find major combination
        s1, s2, s3 = triplet.names
        names = [s1, s2, s3]
        parent1, parent2 = find_parent(s1, s2, s3, tree)

        for i in triplet.names:
            if not i in (parent1, parent2):
                recombinant = i

        # use this variable as indicator to see which group is major and which is minor        
        signal = names.index(recombinant)

        # functions input order is relative
        signal_orders = [
            (seq_2_3, seq_1_2, seq_1_3),  # min == 0 2==3
            (seq_1_3, seq_2_3, seq_1_2),  # min == 1 1==3
            (seq_1_2, seq_1_3, seq_2_3)   # min == 2 1==2
        ]

        x, y, z = signal_orders[signal]
        min1, min2 = self.find_signal(x, y, z, signal)

        for signals in [min1, min2]:
            for ind, (start, end, signal, lr) in enumerate(signals):
                a,b,c = triplet.sequences

                # TODO make sequence d an outgroup rather than randomized sequences of the three
                # this is a much cleaner way
                sequence_orders = [
                    (b, c, a),
                    (a, c, b),
                    (a, b, c)
                ]

                x, y, z = sequence_orders[signal]
                start = self.adjust_nt_f(x, y, z, start, end)
                end = self.adjust_nt_r(x, y, z, end, start)
                
                # if the window does not continue, don't draw new window and discard
                # TODO: in second thought, I have no idea if this is right or not, need to check
                if end == start:
                    continue 

                p_val = norm.sf(self.shuffle(x, y, z, out, lr))
                
                if p_val < 0.05/(len(range(0, self.align.shape[1], self.step_size))):
                    self.raw_results.append((recombinant, (parent1, parent2), start, end, p_val))

        return
