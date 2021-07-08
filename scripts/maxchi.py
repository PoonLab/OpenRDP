from scipy.signal import find_peaks, medfilt
from scipy.ndimage import gaussian_filter1d
import numpy as np
from scripts.common import generate_triplets, calculate_chi2
import matplotlib.pyplot as plt


class MaxChi:
    def __init__(self, align, names, max_pvalue=0.05, win_size=200, strip_gaps=True, fixed_win_size=True,
                 num_var_sites=None, frac_var_sites=None, settings=None):
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

        else:
            self.align = align
            self.s_names = names
            self.results = {}
            self.max_pvalue = max_pvalue
            self.fixed_win_size = win_size
            self.strip_gaps = strip_gaps
            self.fixed_win = fixed_win_size
            self.num_var_sites = num_var_sites
            self.frac_var_sites = frac_var_sites

        self.results = {}

    def set_options_from_config(self, settings):
        """
        Set the parameters of MaxChi from the config file
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

        if self.num_var_sites < 0:
            print("Invalid option for 'num_var_sites'.\nUsing default value (70) instead.")
            self.num_var_sites = 70

        if self.frac_var_sites > 1 or self.frac_var_sites < 0:
            print("Invalid option for 'frac_var_sites'.\nUsing default value (0.1) instead.")
            self.frac_var_sites = 0.1




    def execute(self, num_trp):
        """
        Executes the MaxChi algorithm
        :param num_trp: the number of triplets
        """
        trp_count = 1
        for trp_idx in generate_triplets(self.align):
            print("Scanning triplet {} / {}".format(trp_count, num_trp))
            trp_count += 1

            # 1. Select the 3 processed sequences
            seqs = []
            for idx in trp_idx:
                seqs.append(self.new_align[idx])

            # 2. Sample two sequences
            pairs = ((0, 1), (1, 2), (2, 0))
            for i, j in pairs:
                seq1 = seqs[i]
                seq2 = seqs[j]

                # Prepare dictionary to store sequences involved in recombination
                names = tuple(sorted([self.s_names[trp_idx[i]], self.s_names[trp_idx[j]]]))
                self.results[names] = []

                # Initialize lists to map chi2 and p-values to window positions
                chi2_values = np.zeros(self.align.shape[1])
                p_values = np.ones(self.align.shape[1])       # Map window position to p-values

                # Get the size of the first window
                win_size = self.get_win_size(offset=0)

                # Slide along the sequences
                half_win_size = int(win_size // 2)
                for k in range(self.new_align.shape[1] - win_size):

                    # Get the left and right half of the window
                    reg1_left = seq1[k: half_win_size + k]
                    reg2_left = seq2[k: half_win_size + k]
                    reg1_right = seq1[k + half_win_size: k + win_size]
                    reg2_right = seq2[k + half_win_size: k + win_size]

                    reg1 = seq1[k: k + win_size]
                    reg2 = seq2[k: k + win_size]

                    s = np.sum(reg1 != reg2)
                    r = np.sum(reg1_left != reg2_left)

                    c_table = [[0, 0],
                               [0, 0]]

                    # Compute contingency table for each window position
                    r_matches = np.sum((reg1_right == reg2_right))
                    c_table[0][0] = int(r_matches)
                    c_table[0][1] = half_win_size - r_matches

                    l_matches = np.sum((reg1_left == reg2_left))
                    c_table[1][0] = int(l_matches)
                    c_table[1][1] = half_win_size - l_matches

                    n = self.win_size
                    k2 = half_win_size

                    cur_val = (float(n) * (k2 * s - n * r) * (k2 * s - n * r)) / (float(k2 * s) * (n - k2) * (n - s))

                    # Compute chi-squared value
                    chi2, p_value = calculate_chi2(c_table, self.max_pvalue)
                    if chi2 is not None and p_value is not None:
                        # Insert p-values and chi2 values so they correspond to positions in the original alignment
                        chi2_values[self.poly_sites[k + half_win_size]] = cur_val  # centred window
                        p_values[self.poly_sites[k + half_win_size]] = p_value

                    win_size = self.get_win_size(k)

                # Smooth chi2-values
                chi2_values = gaussian_filter1d(chi2_values, 1.5)
                # p_values = gaussian_filter1d(p_values, 1.5)
                # self.plot_chi2_values(chi2_values, p_values)

                peaks = find_peaks(chi2_values, distance=self.win_size)
                for k, peak in enumerate(peaks[0]):
                    search_win_size = 1
                    while peak - search_win_size > 0\
                            and peak + search_win_size < len(chi2_values) - 1\
                            and chi2_values[peak + search_win_size] > 0.3 * chi2_values[peak]\
                            and chi2_values[peak - search_win_size] > 0.3 * chi2_values[peak]:
                        search_win_size += 1

                    if chi2_values[peak + search_win_size] > chi2_values[peak - search_win_size]:
                        aln_pos = (int(peak), int(peak + search_win_size + win_size))
                    else:
                        aln_pos = (int(peak - search_win_size), int(peak + win_size))

                    # Check that breakpoint has not already been detected
                    if (*aln_pos, p_values[k]) not in self.results[names]:
                        self.results[names].append((*aln_pos, p_values[peak]))

        return self.results

    @staticmethod
    def plot_chi2_values(self, chi_values, p_values):
        plt.figure()
        p = -np.log(p_values)
        plt.plot(chi_values)
        plt.plot(p)
        plt.show()
