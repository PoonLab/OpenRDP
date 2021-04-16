from scipy.stats import chi2_contingency
import numpy as np
from scipy.signal import find_peaks


class Chimaera:
    def __init__(self, align, s_names, settings=None, win_size=200, strip_gaps=True, fixed_win_size=True, num_var_sites=None,
                 frac_var_sites=None):
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
            self.win_size = win_size
            self.strip_gaps = strip_gaps
            self.fixed_win_size = fixed_win_size
            self.num_var_sites = num_var_sites
            self.frac_var_sites = frac_var_sites

        self.s_names = s_names
        self.new_align, self.poly_sites = self.remove_monomorphic_sites(align)
        self.results = []

    def set_options_from_config(self, settings):
        """
        Set the parameters of Chimaera from the config file
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

        if not self.fixed_win_size:
            if self.num_var_sites < 0:
                print("Invalid option for 'num_var_sites'.\nUsing default value (60) instead.")
                self.num_var_sites = 60

            if self.frac_var_sites > 1 or self.frac_var_sites < 0:
                print("Invalid option for 'frac_var_sites'.\nUsing default value (0.1) instead.")
                self.frac_var_sites = 0.1

    @staticmethod
    def remove_monomorphic_sites(align):
        """
        Remove monomorphic sites
        :param align: n x m numpy array where n is the length of the alignment and m is the number of sequences
        :return: a tuple containing the polymorphic sites and the positions of polymorphic sites
        """
        poly_sites = []
        # Find positions of polymorphic sites
        for i in range(align.shape[1]):
            col = align[:, i]
            if not np.all(col == col[0]):
                poly_sites.append(i)

        # Build "new alignment"
        new_align = align[:, poly_sites]

        return new_align, poly_sites

    def execute(self, triplet):
        """
        Executes the Chimaera algorithm
        :param triplet: a list of three sequences
        """

        # Get the triplet sequences
        trp_seqs = []
        for seq_num in triplet:
            trp_seqs.append(self.new_align[seq_num])

        combos = [(0, 1, 2), (1, 2, 0), (2, 1, 0)]
        for run in combos:
            recombinant = trp_seqs[run[0]]
            parental_1 = trp_seqs[run[1]]
            parental_2 = trp_seqs[run[2]]

            rec_name = self.s_names[triplet[run[0]]]
            p1_name = self.s_names[triplet[run[1]]]
            p2_name = self.s_names[triplet[run[2]]]

            chi2_values = []
            p_values = []

            # Remove sites where neither the parental match the recombinant
            informative_sites = []
            aln = np.array([recombinant, parental_1, parental_2])
            # Find indices for informative sites
            for i in range(aln.shape[1]):
                col = aln[:, i]
                if col[0] == col[1] or col[0] == col[2]:
                    informative_sites.append(i)

            # Build new alignment with uninformative sites removed
            new_aln = aln[:, informative_sites]

            # 2. Compress "recombinant" into bit-strings
            comp_seq = []
            for i in range(new_aln.shape[1]):
                if new_aln[0][i] == new_aln[1][i]:
                    comp_seq.append(0)
                elif new_aln[0][i] == new_aln[2][i]:
                    comp_seq.append(1)

            # Move sliding window along compressed sequence , 1 position at a time
            # Slide along the sequences
            half_win_size = int(self.win_size // 2)
            for k in range(len(comp_seq) - self.win_size):
                reg_left = comp_seq[k: half_win_size + k]
                reg_right = comp_seq[k + half_win_size: k + self.win_size]

                c_table = [[0, 0],
                           [0, 0]]

                # Compute contingency table for each window position
                count_r_ones = np.count_nonzero(reg_right)
                c_table[0][0] = count_r_ones
                c_table[0][1] = half_win_size - count_r_ones

                count_l_ones = np.count_nonzero(reg_left)
                c_table[1][0] = count_l_ones
                c_table[1][1] = half_win_size - count_l_ones

                # Compute chi-squared value
                if c_table[0][0] == 0 or c_table[0][1] == 0:
                    chi2_values.append(0)
                    p_values.append(0)

                else:
                    chi2, p_value, _, _ = chi2_contingency(c_table)
                    chi2_values.append(chi2)
                    p_values.append(p_value)

            peaks = find_peaks(chi2_values, wlen=self.win_size, distance=self.win_size)
            for i, peak in enumerate(peaks[0]):
                search_win_size = 1
                while peak - search_win_size > 0 \
                        and peak + search_win_size < len(chi2_values) - 1 \
                        and chi2_values[peak + search_win_size] > 0.5 * chi2_values[peak] \
                        and chi2_values[peak - search_win_size] > 0.5 * chi2_values[peak]:
                    search_win_size += 1

                if chi2_values[peak + search_win_size] > chi2_values[peak - search_win_size]:
                    aln_pos = (int(peak), int(peak + search_win_size + self.win_size))
                else:
                    aln_pos = (int(peak - search_win_size), int(peak + self.win_size))

                self.results.append((rec_name, p1_name, p2_name, *aln_pos, p_values[i]))

        return
