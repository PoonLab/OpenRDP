import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
from scipy.stats import chisqprob
import json

from .common import calculate_chi2, identify_recombinant


class MaxChi:
    def __init__(self, align, max_pvalue=0.05, win_size=70, strip_gaps=True, fixed_win_size=True,
                 num_var_sites=None, frac_var_sites=None, settings=None, ref_align=None, verbose=False):
        """
        Constructs a MaxChi Object
        :param win_size: Size of the sliding window
        :param strip_gaps: Determines whether gaps are stripped or kept
        :param fixed_win_size: Determines the type of sliding window (fixed or variable)
        :param num_var_sites: The number of variable sites in a variable length sliding window
        :param frac_var_sites: The fraction of variable sites in a variable length sliding window
        """
        if settings:
            self.set_options_from_config(settings)
            self.validate_options(align)

        else: # pragma: no cover
            self.align = align
            self.raw_results = {}
            self.max_pvalues = max_pvalue           
            self.fixed_win_size = win_size
            self.strip_gaps = strip_gaps
            self.fixed_win = fixed_win_size
            self.num_var_sites = num_var_sites
            self.frac_var_sites = frac_var_sites

        # see issue #99
        if len(align[0]) < self.win_size:
            win_size = len(align[0])

        self.raw_results = []
        self.results = []
        self.name = 'maxchi'

    def set_options_from_config(self, settings):
        """
        Set the parameters of MaxChi from the config file
        :param settings: a dictionary of settings
        """
        self.win_size = abs(int(settings['win_size']))
        self.max_pvalues = abs(float(settings['max_pvalue']))

        if settings['strip_gaps'] == 'False':
            self.strip_gaps = False
        else:
            self.strip_gaps = True

        if settings['fixed_win_size'] == 'True':
            self.fixed_win_size = True
        else:
            self.fixed_win_size = False

        self.num_var_sites = abs(int(settings['num_var_sites']))
        self.frac_var_sites = abs(float(settings['frac_var_sites']))

    def validate_options(self, align):
        """
        Check if the options from the config file are valid
        If the options are invalid, the default value will be used instead
        """
        if self.win_size < 0 or self.win_size > align.shape[1]:
            print("Invalid option for 'win_size'.\nUsing default value (200) instead.")
            self.win_size = 200

        if self.num_var_sites < 0:
            print("Invalid option for 'num_var_sites'.\nUsing default value (70) instead.")
            self.num_var_sites = 70

        if self.frac_var_sites > 1 or self.frac_var_sites < 0:
            print("Invalid option for 'frac_var_sites'.\nUsing default value (0.1) instead.")
            self.frac_var_sites = 0.1

    @staticmethod
    def get_window_positions(seq1, seq2, k, win_size):
        """
        Get the left and right half of the window
        :param seq1: the first sequence
        :param seq2: the second sequence
        :param k: the offset for the window
        :param win_size: the size of the window
        :return: the left and right regions on either side of the partition
        """
        half_win_size = int(win_size // 2)
        reg1_left = seq1[k: half_win_size + k]
        reg2_left = seq2[k: half_win_size + k]
        reg1_right = seq1[k + half_win_size: k + win_size]
        reg2_right = seq2[k + half_win_size: k + win_size]

        reg1 = seq1[k: k + win_size]
        reg2 = seq2[k: k + win_size]

        return reg1_left, reg2_left, reg1_right, reg2_right, reg1, reg2

    @staticmethod
    def compute_contingency_table(reg1_right, reg2_right, reg1_left, reg2_left, half_win_size):
        """
        Calculate the number of variable sites on either side of the partition
        :param reg1_right: the right half of the window from the first sequence
        :param reg2_right: the left half of the window from the second sequence
        :param reg1_left: the left half of the window from the first sequence
        :param reg2_left: the left half of the window from the second sequence
        :param half_win_size: half the width of the window
        :return: the contingency table
        """
        # Record the totals for the rows and columns
        c_table = [[0, 0, 0],
                   [0, 0, 0],
                   [0, 0, 0]]

        # Compute contingency table for each window position
        r_matches = np.sum((reg1_right == reg2_right))
        c_table[0][0] = int(r_matches)
        c_table[0][1] = half_win_size - r_matches

        l_matches = np.sum((reg1_left == reg2_left))
        c_table[1][0] = int(l_matches)
        c_table[1][1] = half_win_size - l_matches

        # Sum the rows and columns
        c_table[0][2] = c_table[0][0] + c_table[0][1]
        c_table[1][2] = c_table[1][0] + c_table[1][1]
        c_table[2][0] = c_table[0][0] + c_table[1][0]
        c_table[2][1] = c_table[0][1] + c_table[1][1]
        c_table[2][2] = c_table[0][2] + c_table[1][2]

        return c_table

    def execute(self, triplet):
        """
        Executes the MaxChi algorithm
        :param triplet: a triplet object
        """

        # 1. Sample two sequences
        pairs = ((0, 1), (1, 2), (2, 0))
        chi2_values = [np.zeros(triplet.sequences.shape[1])] * 3

        for pair, (i, j) in enumerate(pairs):
            seq1 = triplet.sequences[i]
            seq2 = triplet.sequences[j]

            # Initialize lists to map chi2 and p-values to window positions

            # Get the size of the first window
            win_size = triplet.get_win_size(0, self.win_size, self.fixed_win_size, self.num_var_sites,
                                            self.frac_var_sites)

            # Slide along the sequences
            k2 = int(win_size // 2) # using these variable names to be consistant with the original paper
            for k in range(triplet.poly_sites_align.shape[1] - win_size):

                reg1_left, reg2_left, reg1_right, reg2_right, reg1, reg2 = self.get_window_positions(seq1, seq2,
                                                                                                     k, win_size)

                c_table = self.compute_contingency_table(reg1_right, reg2_right,
                                                         reg1_left, reg2_left, k2)


                # Compute chi-squared value
                chi2, p_value = calculate_chi2(c_table, self.max_pvalues)
                if chi2 is not None and p_value is not None:
                    # Insert p-values and chi2 values so they correspond to positions in the original alignment
                    chi2_values[pair][triplet.poly_sites[k + k2]] = chi2  # centred window

                # win_size = triplet.get_win_size(k, self.win_size, self.fixed_win_size, self.num_var_sites,
                #                                 self.frac_var_sites)

        with open('chi2.json', 'w') as file:
            json.dump(
                [gaussian_filter1d(x, 1.5).tolist() for x in chi2_values], 
                file
            )

        # Smooth chi2-values
        
        # p_values = gaussian_filter1d(p_values, 1.5)
        # self.plot_chi2_values(chi2_values, p_values)

        # method as according to Darren:
        peaks = find_peaks(chi2_values, distance=self.win_size)

        # 1. find highest peak amongst all three pairs
        best = [0,0,0]
        for pair, values in enumerate(chi2_values):
            for pos, chi in enumerate(values):
                if chi > best[1]:
                    best = [pair, chi, pos]

        # 2. find optimized window of maxChi values by expanding and contracting the window and recalc Chi2 while counting fails

        ##### problem outlined:
        # assuming you have an array of integers, n, with a start index, k, and a left and right value index of the window. You have to expand and contract this window to get the max value. 
        # However, you are limited by your actions. You can either increase one end of the window by 1 position or contract the window by 1 position. If the new window has a lower max value than the previous best, you add a counter to fail. if you reach 100 fails you automatically return the best max window.
        # NOTE BY WILL: i coulnd't find how RDP solved this problem, so i'm going to use a greedy algorithm to expand the window similar to LeetCode 1004
        # NOTE 2: a dp solution is probably the optimal solution. gets out of [100, ... n * fails limit -1, 10001]

        if chisqprob(best[1],1) < 0.05/len(chi2_values[0] * 3): #BF correction

            # using degree of freedom of 1 because of 2x2 chi2 test
            initial_window = 5
            fail_count = 0
                        
            # get pair that mattered and get the windows
            i, j = pairs[best[0]]
            seq1 = triplet.sequences[i]
            seq2 = triplet.sequences[j]

            # initalize variables to expand the window
            left, right = best[2] - initial_window // 2, best[2] + initial_window // 2

            # TODO: check source code, rn the 100 is arbitrary
            fail_count = 0
            while fail_count < 100:
                candidates = []

                # expand left
                if left > 0:
                    el = left - 1
                    c_table = self.compute_contingency_table(
                        seq1[k:right], seq2[k:right],
                        seq1[el:k], seq2[el:k], (right - el)//2
                        )
                    chi2, p_value = calculate_chi2(c_table, self.max_pvalues)
                    candidates.append(('expand_left', chi2, el, right))

                # expand right
                if right < len(seq1) - 1:
                    er = right + 1
                    c_table = self.compute_contingency_table(
                        seq1[k:er], seq2[k:er],
                        seq1[left:k], seq2[left:k], (right - er)//2
                    )
                    chi2, p_value = calculate_chi2(c_table, self.max_pvalues)
                    candidates.append(('expand_right', chi2, left, er))

                # contract left
                if left < k:
                    cl = left + 1
                    c_table = self.compute_contingency_table(
                        seq1[k:right], seq2[k:right],
                        seq1[cl:k], seq2[cl:k], (right - cl)//2
                    )
                    chi2, p_value = calculate_chi2(c_table, self.max_pvalues)
                    candidates.append(('contract_left', chi2, cl, right))

                # contract right
                if right > k:
                    cr = right - 1
                    c_table = self.compute_contingency_table(
                        seq1[k:cr], seq2[k:cr],
                        seq1[left:k], seq2[left:k], (cr - left)//2
                    )
                    chi2, p_value = calculate_chi2(c_table, self.max_pvalues)
                    candidates.append(('contract_right', chi2, left, cr))

                if not candidates:
                    break  # no valid moves

                # Choose the best candidate
                best_move, best_chi2, new_left, new_right = max(candidates, key=lambda x: x[1])

                if best_chi2 > best_score:
                    best_score = best_chi2
                    best_left = new_left
                    best_right = new_right
                    fail_count = 0
                else:
                    fail_count += 1

                # Apply the move explicitly
                if best_move == 'expand_left':
                    left -= 1
                elif best_move == 'expand_right':
                    right += 1
                elif best_move == 'contract_left':
                    left += 1
                elif best_move == 'contract_right':
                    right -= 1

        # 3, determine if left or right side of the new window is the best

        # breakpoint one is the original peak, now find the biggest peak for the second breakpoint
        break2 = [0, 0] # ind, val
        for ind, val in enumerate(chi2_values[best[0]][k+1:right+1]):
            if val > break2[1]:
                break2[1] = val
                break2[0] = ind
        
        # check left side of window now
        for ind, val in enumerate(chi2_values[best[0]][left:k]):
            if val > break2[1]:
                break2[1] = val
                break2[0] = ind
        
        # 4, take higher chi2 value as the p-value
        p_val = min(chisqprob(chi2_values[break2[0]],1), chisqprob(chi2_values[k], 1))

        # now add back to solutions
        # TODO: check if this is the right way
        aln_pos = [min(break2[0], k), max(break2[0])]

        rec_name, parents = identify_recombinant(triplet, aln_pos)

        self.raw_results.append((rec_name, parents, aln_pos, p_val))


        return