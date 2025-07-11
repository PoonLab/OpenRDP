import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
from scipy.stats import chi2 as chisq

from .common import calculate_chi2, identify_recombinant


class Chimaera:
    def __init__(self, align, max_pvalue=0.05, win_size=200, strip_gaps=True, fixed_win_size=True,
                 num_var_sites=None, frac_var_sites=None, settings=None, ref_align=None, verbose=False):
        """
        Constructs a Chimaera Object
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
            self.max_pvalues = max_pvalue
            self.win_size = win_size
            self.strip_gaps = strip_gaps
            self.fixed_win_size = fixed_win_size
            self.num_var_sites = num_var_sites
            self.frac_var_sites = frac_var_sites

        self.align = align
        self.raw_results = []
        self.results = []
        self.name = 'chimaera'

    def set_options_from_config(self, settings):
        """
        Set the parameters of Chimaera from the config file
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

        if not self.fixed_win_size:
            if self.num_var_sites < 0:
                print("Invalid option for 'num_var_sites'.\nUsing default value (60) instead.")
                self.num_var_sites = 60

            if self.frac_var_sites > 1 or self.frac_var_sites < 0:
                print("Invalid option for 'frac_var_sites'.\nUsing default value (0.1) instead.")
                self.frac_var_sites = 0.1

    @staticmethod
    def get_window_positions(comp_seq, k, win_size):
        """
        Get the left and right half of the window
        :param comp_seq: the compressed recombinant sequence
        :param k: the offset for the window
        :param win_size: the size of the window
        :return: the left and right regions on either side of the partition
        """
        half_win_size = int(win_size // 2)
        reg_left = comp_seq[k: half_win_size + k]
        reg_right = comp_seq[k + half_win_size: k + win_size]

        return reg_left, reg_right

    @staticmethod
    def compute_contingency_table(reg_left, reg_right, half_win_size):
        """
        Calculate the number of variable sites on either side of the partition
        :param reg_right: the right half of the window
        :param reg_left: the left half of the window
        :param half_win_size: half the width of the window
        :return: the contingency table
        """

        c_table = [[0, 0, 0],
                   [0, 0, 0],
                   [0, 0, 0]]

        # Compute contingency table for each window position
        count_r_ones = np.count_nonzero(reg_right)
        c_table[0][0] = count_r_ones
        c_table[0][1] = abs(half_win_size - count_r_ones)

        count_l_ones = np.count_nonzero(reg_left)
        c_table[1][0] = count_l_ones
        c_table[1][1] = abs(half_win_size - count_l_ones)

        # Sum the rows and columns
        c_table[0][2] = c_table[0][0] + c_table[0][1]
        c_table[1][2] = c_table[1][0] + c_table[1][1]
        c_table[2][0] = c_table[0][0] + c_table[1][0]
        c_table[2][1] = c_table[0][1] + c_table[1][1]
        c_table[2][2] = c_table[0][2] + c_table[1][2]

        return c_table

    @staticmethod
    def compress_triplet_aln(new_aln):
        """
        Compress the sequences into a string of 0s and 1s where:
            "1" represents a match between the "recombinant" (the first sequence in the alignment) and parent A
            "0" represents a match between the "recombinant" and parent B
        Since the alignment consists only of informative sites, the "recombinant" must match either parent A or parent B
        :param new_aln: alignment containing only informative sites
            The first sequence is the "recombinant", the second is parent A, and the third is parent B
        :return: a bitstring that represents the alignment
        """
        comp_seq = []
        for i in range(new_aln.shape[1]):
            if new_aln[0][i] == new_aln[1][i]:
                comp_seq.append(0)
            elif new_aln[0][i] == new_aln[2][i]:
                comp_seq.append(1)
        return comp_seq

    def refine_nt(self, seq, start, end):
        """
        adjust the range of nucleotide positions to match the significant patterns for start
        
        seq, compressed bit string
        start int, left index breakpoint
        end int, right index breakpoint
        """
        adj = 0
        l, r = False, False
        final = [None,None] # final left and right index
        while start - adj >= 0 or  start + adj < len(seq):
            # prioritize the contraction over expansion
            if start + adj < len(seq):
                if seq[start + adj] == 1:
                    final[0] = start + adj
                    break
            if start - adj < len(seq):
                if seq[start - adj] == 1:
                    final[0] = start - adj
                    break
            adj += 1
        
        if not final[0]: # couldn't fine a refined spot
            final[0] = start

        adj = 0
        while end - adj > start or end + adj < len(seq):
            # prioritize the contraction over expansion
            if start + adj < len(seq):
                if seq[end - adj] == 1:
                    final[1] = end - adj
                    break
            if start - adj < len(seq):
                if seq[end + adj] == 1:
                    final[1] = end + adj
                    break
            adj += 1

        if not final[1]:
            final[1] = end 
        return final

    def execute(self, triplet):
        """
        Executes the Chimaera algorithm
        :param triplets: a triplet object
        """
        # Try every possible combination of "parentals" and "recombinant"
        combos = [(0, 1, 2), (1, 2, 0), (2, 1, 0)]
        for run in combos:

            # Initialize lists to map chi2 and p-values to window positions
            chi2_values = np.zeros(self.align.shape[1])
            p_values = np.ones(self.align.shape[1])  # Map window position to p-values

            # Build new alignment with uninformative sites removed
            # Uninformative sites where neither parental matches recombinant
            new_aln = triplet.sequences[:, triplet.info_sites]

            # Rearrange the new alignment to try every possible combination of "parentals" and "recombinant"
            new_aln = np.array([new_aln[run[0]], new_aln[run[1]], new_aln[run[2]]])

            # Compress into a bitstring
            comp_seq = self.compress_triplet_aln(new_aln)

            # Get the size of the first window
            win_size = triplet.get_win_size(0, self.win_size, self.fixed_win_size, self.num_var_sites,
                                            self.frac_var_sites)

            # Move sliding window along compressed sequence, 1 position at a time
            # Slide along the sequences
            half_win_size = int(win_size // 2)
            for k in range(len(comp_seq) - win_size):
                reg_left, reg_right = self.get_window_positions(comp_seq, k, win_size)

                c_table = self.compute_contingency_table(reg_left, reg_right, half_win_size)

                # Compute chi-squared value
                chi2, p_value = calculate_chi2(c_table)
                if chi2 is not None and p_value is not None:
                    # Insert p-values and chi2 values so they correspond to positions in the original alignment
                    chi2_values[triplet.poly_sites[k + half_win_size]] = chi2  # centred window
                    p_values[triplet.poly_sites[k + half_win_size]] = p_value

                win_size = triplet.get_win_size(k, self.win_size, self.fixed_win_size, self.num_var_sites,
                                                self.frac_var_sites)

            # Smooth chi2-values
            chi2_values = gaussian_filter1d(chi2_values, 1.5)

            # Locate "peaks" in chi2 values as "peaks" represent potential breakpoints
            best = (0,0) # maxCHI2, ind
            for i, value in enumerate(chi2_values):
                if value > best[0]:
                    best = (value, i)

            # start same algorithm in maxchi
            k = best[1]
            if chisq.sf(best[0],1) < 0.05/len(comp_seq) - win_size: #bf correction
                
                initial_window = self.win_size

                # initial variables
                left, right = best[1] - initial_window // 2, best[1] + initial_window // 2

                fail_count, best_score = 0, 0
                while fail_count < 100:
                    print(left, right)
                    candidates = []

                    half_win_size = (right - left) //2

                    # expand left and right at the same time
                    if left > 0 and right < len(chi2_values):
                        el = left - 1
                        er = right + 1
                        reg_left, reg_right = comp_seq[el:k], comp_seq[k:er]
                        c_table = self.compute_contingency_table(reg_left, reg_right, half_win_size)
                        chi2, p_value = calculate_chi2(c_table)
                        if chi2:
                            candidates.append(('expand', chi2))

                    # contract left
                    if left < k:
                        cl = left + 1
                        reg_left, reg_right = comp_seq[cl:k], comp_seq[k:right]
                        c_table = self.compute_contingency_table(reg_left, reg_right, half_win_size)
                        chi2, p_value = calculate_chi2(c_table)
                        if chi2:
                            candidates.append(('contract_left', chi2))

                    # contract right
                    if right > k:
                        cr = right - 1
                        reg_left, reg_right = comp_seq[left:k], comp_seq[k:cr]
                        c_table = self.compute_contingency_table(reg_left, reg_right, half_win_size)
                        chi2, p_value = calculate_chi2(c_table)
                        if chi2:
                            candidates.append(('contract_right', chi2))      

                    if not candidates:
                        break

                    else:
                        best_move, best_chi2 = max(candidates, key=lambda x: x[1])

                        if best_chi2 > best_score:
                            best_score = best_chi2

                        else:
                            fail_count += 1

                        # Apply the move
                        if best_move == 'expand':
                            left -= 1
                            right += 1
                        elif best_move == 'contract_left':
                            left += 1
                        elif best_move == 'contract_right':
                            right -= 1

                left_peak = chi2_values[left]
                right_peak = chi2_values[right]
                midpoint = left - (left - right)//2

                # determine which way is the second breakpoint
                if left_peak > right_peak:
                    primary_breakpoint = left
                else:
                    primary_breakpoint = midpoint

                if primary_breakpoint == left:
                    second_break = midpoint
                else:
                    second_break = right

                # new midpoint for recalculation for chi2
                midpoint = second_break - (second_break - primary_breakpoint)//2

                reg_left, reg_right = comp_seq[primary_breakpoint:midpoint], comp_seq[midpoint:second_break]

                c_table = self.compute_contingency_table(reg_left, reg_right,  (second_break - primary_breakpoint)//2)
                chi2, p_value = calculate_chi2(c_table)

                if chi2:
                    final_chi2 = max(best_score, chi2)
                else:
                    final_chi2 = best_score
                final_p = chisq.sf(final_chi2, df=1)
                
                aln_pos = self.refine_nt(comp_seq, primary_breakpoint, second_break)
                rec_name, parents = identify_recombinant(triplet, aln_pos)
                self.raw_results.append((rec_name, parents, *aln_pos, final_p))


        return
