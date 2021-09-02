import copy

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
from scipy.stats import pearsonr

from scripts.common import calculate_chi2, jc_distance, all_items_equal


class MaxChi:
    def __init__(self, align, max_pvalue=0.05, win_size=200, strip_gaps=True, fixed_win_size=True,
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
            self.raw_results = {}
            self.max_pvalue = max_pvalue
            self.fixed_win_size = win_size
            self.strip_gaps = strip_gaps
            self.fixed_win = fixed_win_size
            self.num_var_sites = num_var_sites
            self.frac_var_sites = frac_var_sites

        self.raw_results = []
        self.results = []

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

    def execute(self, triplets, quiet=False):
        """
        Executes the MaxChi algorithm
        :param triplets: a list of triplet objects, representing the sequence triplets
        :param quiet: report progress
        """
        trp_count = 1
        total_num_trps = len(triplets)
        for triplet in triplets:
            if not quiet:
                print("Scanning triplet {} / {}".format(trp_count, total_num_trps))
            trp_count += 1

            # 1. Sample two sequences
            pairs = ((0, 1), (1, 2), (2, 0))
            for i, j in pairs:
                seq1 = triplet.sequences[i]
                seq2 = triplet.sequences[j]

                # Initialize lists to map chi2 and p-values to window positions
                chi2_values = np.zeros(triplet.sequences.shape[1])
                p_values = np.ones(triplet.sequences.shape[1])  # Map window position to p-values

                # Get the size of the first window
                win_size = triplet.get_win_size(0, self.win_size, self.fixed_win_size, self.num_var_sites,
                                                self.frac_var_sites)

                # Slide along the sequences
                half_win_size = int(win_size // 2)
                for k in range(triplet.poly_sites_align.shape[1] - win_size):

                    reg1_left, reg2_left, reg1_right, reg2_right, reg1, reg2 = self.get_window_positions(seq1, seq2,
                                                                                                           k, win_size)

                    s = np.sum(reg1 != reg2)
                    r = np.sum(reg1_left != reg2_left)

                    c_table = self.compute_contingency_table(reg1_right, reg2_right,
                                                               reg1_left, reg2_left, half_win_size)

                    n = self.win_size
                    k2 = half_win_size

                    if (float(k2 * s) * (n - k2) * (n - s)) == 0:
                        continue
                    else:
                        cur_val = (float(n) * (k2 * s - n * r) * (k2 * s - n * r)) / (float(k2 * s) * (n - k2) * (n - s))

                    # Compute chi-squared value
                    chi2, p_value = calculate_chi2(c_table, self.max_pvalue)
                    if chi2 is not None and p_value is not None:
                        # Insert p-values and chi2 values so they correspond to positions in the original alignment
                        chi2_values[triplet.poly_sites[k + half_win_size]] = cur_val  # centred window
                        p_values[triplet.poly_sites[k + half_win_size]] = p_value

                    win_size = triplet.get_win_size(k, self.win_size, self.fixed_win_size, self.num_var_sites,
                                                    self.frac_var_sites)

                # Smooth chi2-values
                chi2_values = gaussian_filter1d(chi2_values, 1.5)
                # p_values = gaussian_filter1d(p_values, 1.5)
                # self.plot_chi2_values(chi2_values, p_values)

                peaks = find_peaks(chi2_values, distance=self.win_size)
                for k, peak in enumerate(peaks[0]):
                    search_win_size = 1
                    while peak - search_win_size > 0 \
                            and peak + search_win_size < len(chi2_values) - 1 \
                            and chi2_values[peak + search_win_size] > 0.3 * chi2_values[peak] \
                            and chi2_values[peak - search_win_size] > 0.3 * chi2_values[peak]:
                        search_win_size += 1

                    if chi2_values[peak + search_win_size] > chi2_values[peak - search_win_size]:
                        aln_pos = (int(peak), int(peak + search_win_size + win_size))
                    else:
                        aln_pos = (int(peak - search_win_size), int(peak + win_size))

                    rec_name, parents = self.identify_recombinant(triplet, aln_pos)
                    if ((rec_name, parents, *aln_pos)) not in self.raw_results and p_values[peak] != 1.0:
                        self.raw_results.append((rec_name, parents, *aln_pos, p_values[peak]))

                    self.raw_results = sorted(self.raw_results)

        self.results = self.merge_breakpoints()

        return self.results

    @staticmethod
    def plot_chi2_values(chi_values, p_values):
        plt.figure()
        p = -np.log(p_values)
        plt.plot(chi_values)
        plt.plot(p)
        plt.show()

    def identify_recombinant(self, trp, aln_pos):
        """
        Find the most likely recombinant sequence using the PhPr method described in the RDP5 documentation
        and Weiler GF (1998) Phylogenetic profiles: A graphical method for detecting genetic recombinations
        in homologous sequences. Mol Biol Evol 15: 326–335
        :return: name of the recombinant sequence and the names of the parental sequences
        """
        upstream_dists = []
        downstream_dists = []

        # Get possible breakpoint locations
        for i in range(len(trp.names)):
            upstream_dists.append([])
            downstream_dists.append([])
            for j in range(len(trp.names)):
                if i != j:
                    # Calculate pairwise Jukes-Cantor distances for regions upstream and downstream of breakpoint
                    upstream_dists[i].append(jc_distance(trp.sequences[i][0: aln_pos[0]],
                                                         trp.sequences[j][0: aln_pos[0]]))
                    downstream_dists[i].append(jc_distance(trp.sequences[i][aln_pos[1]: trp.sequences.shape[1]],
                                                           trp.sequences[j][aln_pos[1]: trp.sequences.shape[1]]))

        # Calculate Pearson's correlation coefficient for 2 lists
        r_coeff = [0, 0, 0]
        if all_items_equal(upstream_dists[0]) or all_items_equal(downstream_dists[0]):
            r_coeff[0] = float('NaN')
        elif all_items_equal(upstream_dists[1]) or all_items_equal(downstream_dists[1]):
            r_coeff[1] = float('NaN')
        elif all_items_equal(upstream_dists[2]) or all_items_equal(downstream_dists[2]):
            r_coeff[2] = float('NaN')

        else:
            r_coeff[0], _ = pearsonr(upstream_dists[0], downstream_dists[0])
            r_coeff[1], _ = pearsonr(upstream_dists[1], downstream_dists[1])
            r_coeff[2], _ = pearsonr(upstream_dists[2], downstream_dists[2])

        # Most likely recombinant sequence is sequence with lowest coefficient
        trp_names = copy.copy(trp.names)
        rec_name = trp_names.pop(np.argmin(r_coeff))
        p_names = trp_names

        return rec_name, p_names

    def merge_breakpoints(self):
        """
        Merge overlapping breakpoint locations
        :return: list of breakpoint locations where overlapping intervals are merged
        """
        results_dict = {}
        results = []

        # Gather all regions with the same recombinant
        for i, bp in enumerate(self.raw_results):
            rec_name = self.raw_results[i][0]
            parents = tuple(sorted(self.raw_results[i][1]))
            key = (rec_name, parents)
            if key not in results_dict:
                results_dict[key] = []
            results_dict[key].append(self.raw_results[i][2:])

        # Merge any locations that overlap - eg [1, 5] and [3, 7] would become [1, 7]
        for key in results_dict:
            merged_regions = []
            for region in results_dict[key]:
                region = list(region)
                old_regions = list(results_dict[key])
                for region2 in old_regions:
                    start = region[0]
                    end = region[1]
                    start2 = region2[0]
                    end2 = region2[1]
                    if start <= start2 <= end or start <= end2 <= end:
                        region[0] = min(start,start2)
                        region[1] = max(end, end2)
                        results_dict[key].remove(region2)
                merged_regions.append(region)

            # Output the results
            for region in merged_regions:
                rec_name = key[0]
                parents = key[1]
                start = region[0]
                end = region[1]
                p_value = region[2]
                results.append((rec_name, parents, start, end, p_value))

        return results
