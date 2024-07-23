import re
from .common import identify_recombinant
import numpy as np
from itertools import combinations


class RdpMethod:
    """
    Executes RDP method
    """
    def __init__(self, align, win_size=30, reference=None, min_id=0, max_id=100,
                 settings=None, ref_align=None, verbose=False, max_pvalues = 100):
        """
        @param align:  Numpy array of sequences as lists
        @param win_size:  window length
        @param reference:  str, user-specified sequence to fix in triplets (NOT IMPLEMENTED)
        @param min_id:  int, minimum sequence identity as a percentage (1-100)
                        FIXME: why not float?
        @param max_id:  int, maximum sequence identity as a percentage (1-100)
        @param settings:  dict, optionally initialize member variables
        @param ref_align:  unused
        @param verbose:  unused
        @param max_pvalues:  unused
        """
        if settings:
            self.set_options_from_config(settings)
            self.validate_options()
        else:
            self.win_size = win_size
            self.reference = reference
            self.min_id = min_id
            self.max_id = max_id

        # TODO check for valid p-val since no innate max_pvalue, output has given values greater than 31 before
        self.max_pvalues = max_pvalues
        self.align = align
        self.raw_results = []
        self.results = []
        self.name = 'rdp'

    def set_options_from_config(self, settings):
        """
        Set the parameters of the RDP method from the config file
        :param settings: a dictionary of settings
        """
        self.win_size = int(settings['window_size'])
        self.reference = settings['reference_sequence']
        self.min_id = int(settings['min_identity'])
        self.max_id = int(settings['max_identity'])

    def validate_options(self):
        """
        Check if the options from the config file are valid
        If the options are invalid, the default value will be used instead
        """
        if self.reference == 'None':
            self.reference = None

        if self.win_size < 0:
            print("Invalid option for 'window_size'.\nUsing default value (30) instead.")
            self.win_size = 30

        if self.min_id < 0 or self.min_id > 100:
            print("Invalid option for 'min_identity'.\nUsing default value (0) instead.")
            self.min_id = 0

        if self.max_id < 0 or self.max_id > 100:
            print("Invalid option for 'max_identity'.\nUsing default value (100) instead.")
            self.min_id = 100

    def triplet_identity(self, triplets):
        """
        Calculate the percent identity of each triplet and validate (UNUSED)
        :param triplets: a list of all triplets
        :return: triplets whose identity is greater than the minimum identity and less than the 
                 maximum identity
        """
        trps = []
        for trp in triplets:
            ab = np.array([trp.sequences[0], trp.sequences[1]])
            bc = np.array([trp.sequences[1], trp.sequences[2]])
            ac = np.array([trp.sequences[0], trp.sequences[2]])
            ab, bc, ac = self.pairwise_identity(ab, bc, ac)

            # Include only triplets whose identity is valid
            if self.min_id < ab < self.max_id and self.min_id < bc < self.max_id and self.min_id < ac < self.max_id:
                trps.append(trp)

        return trps

    @staticmethod
    def pairwise_identity(reg_ab, reg_bc, reg_ac):
        """
        Calculate the pairwise identity of each sequence within the triplet

        :param reg_ab:  matrix of size 2 x sequence_length that contains sequences A and B
        :param reg_bc:  matrix of size 2 x sequence_length that contains sequences B and C
        :param reg_ac:  matrix of size 2 x sequence_length that contains sequences A and C
        :return:  list, binary vector of identity between sequences A and B
                  list, binary vector of identity between sequences B and C
                  list, binary vector of identity between sequences A and C
        """
        # FIXME: all sequences should probably have the same length...
        a_b = [int(reg_ab[0, j] == reg_ab[1, j]) for i in range(reg_ab.shape[1])]
        b_c = [int(reg_bc[0, j] == reg_bc[1, j]) for i in range(reg_bc.shape[1])]
        a_c = [int(reg_ac[0, j] == reg_ac[1, j]) for i in range(reg_ac.shape[1])]
        return a_b, b_c, a_c

    def execute(self, triplet):
        """
        Performs RDP detection method for one triplet of sequences
        :param triplet:  object of class Triplet
        :return: the coordinates of the potential recombinant region and the p_value
        """
        G = sum(1 for _ in combinations(range(self.align.shape[0]), 3))  # Number of triplets

        # Get the three pairs of sequences from alignment of phylogenetically informative sites
        ab = np.array([triplet.info_sites_align[0], triplet.info_sites_align[1]])
        bc = np.array([triplet.info_sites_align[1], triplet.info_sites_align[2]])
        ac = np.array([triplet.info_sites_align[0], triplet.info_sites_align[2]])

        ab_id, bc_id, ac_id = self.pairwise_identity(ab, bc, ac)
        percent_ident = [sum(pair) / len(pair) * 100 for pair in [ab_id, bc_id, ac_id]]

        # Include only triplets whose percent identity is in the valid range
        if min(percent_ident) > self.min_id and max(percent_ident) < self.max_id:
            len_trp = triplet.info_sites_align.shape[1]  # sequence length

            # 2. Sliding window over subsequence and calculate average percent identity at each position
            recombinant_regions = ''  # Recombinant regions denoted by ones
            coord = []
            for i in range(len_trp - self.win_size):
                # Calculate percent identity in each window
                pid_ab = sum(ab_id[:, i: self.win_size + i]) / self.win_size * 100
                pid_bc = sum(bc_id[:, i: self.win_size + i]) / self.win_size * 100
                pid_ac = sum(ac_id[:, i: self.win_size + i]) / self.win_size * 100

                # Identify potential recombinant regions
                if pid_ac > pid_ab or pid_bc > pid_ab:
                    recombinant_regions += "1"
                    coord.append(i)  # FIXME: UNUSED
                else:
                    recombinant_regions += "0"

            # 3. locate runs of 
            recomb_idx = [(m.span()) for m in re.finditer('1+', recombinant_regions)]

            # Convert coordinates from  window-level to alignment-level and record number of windows
            coords = []
            for x, y in recomb_idx:
                coords.append((triplet.info_sites[x], triplet.info_sites[y - 1]))

            for coord in coords:
                n = coord[1] - coord[0]     # Length of putative recombinant region

                if n > 0:
                    # m is the proportion of nts in common between either A or B and C in the recombinant region
                    nts_in_a = triplet.sequences[0][coord[0]: coord[1]]
                    nts_in_c = triplet.sequences[2][coord[0]: coord[1]]
                    m = 0
                    for i in range(n):
                        if nts_in_a[i] == nts_in_c[i]:
                            m += 1

                    # p is the proportion of nts in common between either A or B and C in the entire subsequence
                    id_in_seq = 0
                    for j in range(triplet.sequences.shape[1]):
                        if triplet.sequences[0][j] == triplet.sequences[2][j]:
                            id_in_seq += 1
                    p = id_in_seq / triplet.sequences.shape[1]

                    # Calculate p_value as binomial probability
                    val = 0
                    log_n_fact = np.sum(np.log(np.arange(1, n+1)))  # Convert to log space to prevent integer overflow
                    for i in range(m, n):
                        log_i_fact = np.sum(np.log(np.arange(1, i+1)))
                        log_ni_fact = np.sum(np.log(np.arange(1, n-i+1)))
                        with np.errstate(divide='ignore'):  # Ignore floating point error
                            val += np.math.exp(
                                (log_n_fact - (log_i_fact + log_ni_fact)) + n * np.log(p) + (n-i) * np.log(1-p)
                            )

                    uncorr_pvalue = (len_trp / n) * val
                    corr_p_value = G * uncorr_pvalue

                else:
                    corr_p_value = 'NS'

                if corr_p_value != 'NS' and corr_p_value != 0.0:
                    rec_name, parents = identify_recombinant(triplet, coord)
                    self.raw_results.append((rec_name, parents, *coord, corr_p_value))

            return  # FIXME: why is this still here?
