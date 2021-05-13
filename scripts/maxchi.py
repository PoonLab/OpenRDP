from scipy.stats import chi2_contingency
from scipy.signal import find_peaks
import numpy as np
from scripts.common import remove_monomorphic_sites


class MaxChi:
    def __init__(self, align, names, settings=None, win_size=200, strip_gaps=True, fixed_win_size=True, num_var_sites=None,
                 frac_var_sites=None):
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
            self.win_size = win_size
            self.strip_gaps = strip_gaps
            self.fixed_win_size = fixed_win_size
            self.num_var_sites = num_var_sites
            self.frac_var_sites = frac_var_sites

        self.s_names = names
        self.align = align
        self.new_align, self.poly_sites = remove_monomorphic_sites(align)
        self.results = []

    def set_options_from_config(self, settings):
        """
        Set the parameters of MaxChi from the config file
        :param settings: a dictionary of settings
        """
        self.win_size = int(settings['win_size'])

        if settings['strip_gaps'] == 'False':
            self.strip_gaps = False
        else:
            self.strip_gaps = True

        if settings['fixed_win_size'] == 'True':
            self.fixed_win_size = True
        else:
            self.fixed_win_size = False

        self.num_var_sites = int(settings['num_var_sites'])
        self.frac_var_sites = float(settings['frac_var_sites'])

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

    def execute(self, triplet):
        """
        Executes the MaxChi algorithm
        :param triplet: a list of three sequences
        """
        # Get the triplet sequences
        trp_seqs = []
        for seq_num in triplet:
            trp_seqs.append(self.new_align[seq_num])

        # 2. Sample two sequences
        pairs = ((0, 1), (1, 2), (2, 0))
        for i, j in pairs:
            seq1 = trp_seqs[i]
            seq2 = trp_seqs[j]

            seq1_name = self.s_names[triplet[i]]
            seq2_name = self.s_names[triplet[j]]

            chi2_values = []
            p_values = []

            # Slide along the sequences
            half_win_size = int(self.win_size // 2)
            for k in range(self.new_align.shape[1] - self.win_size):
                reg1_left = seq1[k: half_win_size + k]
                reg2_left = seq2[k: half_win_size + k]
                reg1_right = seq1[k + half_win_size: k + self.win_size]
                reg2_right = seq2[k + half_win_size: k + self.win_size]

                c_table = [[0, 0],
                           [0, 0]]

                # Compute contingency table for each window position
                r_matches = np.sum((reg1_right == reg2_right))
                c_table[0][0] = int(r_matches)
                c_table[0][1] = half_win_size - r_matches

                l_matches = np.sum((reg1_left == reg2_left))
                c_table[1][0] = int(l_matches)
                c_table[1][1] = half_win_size - l_matches

                if c_table[0][0] == 0 or c_table[0][1] == 0:
                    chi2_values.append(0)
                    p_values.append(0)

                # Compute chi-squared value
                else:
                    chi2, p_value, _, _ = chi2_contingency(c_table)
                    chi2_values.append(chi2)
                    p_values.append(p_value)

            peaks = find_peaks(chi2_values, wlen=self.win_size, distance=self.win_size)
            for i, peak in enumerate(peaks[0]):
                search_win_size = 1
                while peak - search_win_size > 0\
                        and peak + search_win_size < len(chi2_values) - 1\
                        and chi2_values[peak + search_win_size] > 0.5 * chi2_values[peak]\
                        and chi2_values[peak - search_win_size] > 0.5 * chi2_values[peak]:
                    search_win_size += 1

                if chi2_values[peak + search_win_size] > chi2_values[peak - search_win_size]:
                    aln_pos = (int(peak), int(peak + search_win_size + self.win_size))
                else:
                    aln_pos = (int(peak - search_win_size), int(peak + self.win_size))

                self.results.append((seq1_name, seq2_name, *aln_pos, p_values[i]))

        return self.results
