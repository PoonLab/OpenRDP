import numpy as np
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d
from scripts.common import calculate_chi2


class Chimaera:
    def __init__(self, align, max_pvalue=0.05, win_size=200, strip_gaps=True, fixed_win_size=True,
                 num_var_sites=None, frac_var_sites=None, settings=None):
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
        else:
            self.max_pvalue = max_pvalue
            self.win_size = win_size
            self.strip_gaps = strip_gaps
            self.fixed_win_size = fixed_win_size
            self.num_var_sites = num_var_sites
            self.frac_var_sites = frac_var_sites

        self.align = align
        self.results = {}

    def set_options_from_config(self, settings):
        """
        Set the parameters of Chimaera from the config file
        :param settings: a dictionary of settings
        """
        self.win_size = abs(int(settings['win_size']))
        self.max_pvalue = abs(float(settings['max_pvalue']))

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

    def execute(self, triplets, quiet=False):
        """
        Executes the Chimaera algorithm
        :param triplets: a list of triplet objects, representing the sequence triplets
        :param quiet: report progress
        """
        trp_count = 1
        total_num_trps = len(triplets)
        for triplet in triplets:
            if not quiet:
                print("Scanning triplet {} / {}".format(trp_count, total_num_trps))
            trp_count += 1

            # Try every possible combination of "parentals" and "recombinant"
            combos = [(0, 1, 2), (1, 2, 0), (2, 1, 0)]
            for run in combos:

                rec_name = triplet.get_sequence_name(run[0])
                p1_name = triplet.get_sequence_name(run[1])
                p2_name = triplet.get_sequence_name(run[2])

                # Prepare dictionary to store sequences involved in recombination
                names = tuple([rec_name, p1_name, p2_name])
                self.results[names] = []

                # Initialize lists to map chi2 and p-values to window positions
                chi2_values = np.zeros(self.align.shape[1])
                p_values = np.ones(self.align.shape[1])  # Map window position to p-values

                # Build new alignment with uninformative sites removed
                # Uninformative sites where neither parental matches recombinant
                new_aln = triplet.sequences[:, triplet.info_sites]

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

                    # Using notation from Maynard Smith (1992)
                    n = self.win_size
                    k2 = half_win_size
                    s = np.sum(reg_left == reg_right)
                    r = np.sum(reg_left != reg_right)

                    # Avoid dividing by 0 (both "parental" sequences are identical in the window)
                    if s > 0:
                        cur_val = (float(n) * (k2 * s - n * r) * (k2 * s - n * r)) / (float(k2 * s) * (n - k2) * (n - s))
                    else:
                        cur_val = -1

                    # Compute chi-squared value
                    chi2, p_value = calculate_chi2(c_table, self.max_pvalue)
                    if chi2 is not None and p_value is not None:
                        # Insert p-values and chi2 values so they correspond to positions in the original alignment
                        chi2_values[triplet.poly_sites[k + half_win_size]] = chi2  # centred window
                        p_values[triplet.poly_sites[k + half_win_size]] = p_value

                    win_size = triplet.get_win_size(k, self.win_size, self.fixed_win_size, self.num_var_sites,
                                                    self.frac_var_sites)

                # Smooth chi2-values
                chi2_values = gaussian_filter1d(chi2_values, 1.5)

                # Locate "peaks" in chi2 values as "peaks" represent potential breakpoints
                peaks = find_peaks(chi2_values, distance=self.win_size)
                for k, peak in enumerate(peaks[0]):
                    search_win_size = 1
                    while peak - search_win_size > 0 \
                            and peak + search_win_size < len(chi2_values) - 1 \
                            and chi2_values[peak + search_win_size] > 0.3 * chi2_values[peak] \
                            and chi2_values[peak - search_win_size] > 0.3 * chi2_values[peak]:
                        search_win_size += 1

                    if chi2_values[peak + search_win_size] > chi2_values[peak - search_win_size]:
                        aln_pos = (int(peak), int(peak + search_win_size + self.win_size))
                    else:
                        aln_pos = (int(peak - search_win_size), int(peak + self.win_size))

                    # Check that breakpoint has not already been detected
                    if (*aln_pos, p_values[k]) not in self.results[names]:
                        self.results[names].append((*aln_pos, p_values[peak]))

        return self.results
