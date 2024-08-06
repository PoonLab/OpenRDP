import re
from .common import identify_recombinant
import numpy as np
from itertools import combinations
from scipy.stats import binom
import math


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
        FIXME: I think this should return a True/False value
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
        a_b = [int(reg_ab[0, i] == reg_ab[1, i]) for i in range(reg_ab.shape[1])]
        b_c = [int(reg_bc[0, i] == reg_bc[1, i]) for i in range(reg_bc.shape[1])]
        a_c = [int(reg_ac[0, i] == reg_ac[1, i]) for i in range(reg_ac.shape[1])]
        return a_b, b_c, a_c
        
        
    def execute(self, triplet):
        """
        Performs RDP detection method for one triplet of sequences
        :param triplet:  object of class Triplet
        :return: the coordinates of the potential recombinant region and the p_value
        """
        
        # total number of triplets - to adjust for multiple comparisons?
        G = sum(1 for _ in combinations(range(self.align.shape[0]), 3))

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
            recombinant_regions = ''
            for i in range(len_trp - self.win_size):
                # Calculate percent identity in each window
                try:
                    pid_ab = sum(ab_id[i:(self.win_size+i)]) / self.win_size * 100
                except TypeError:
                    print(i)
                    print(self.win_size)
                    raise
                pid_bc = sum(bc_id[i:(self.win_size+i)]) / self.win_size * 100
                pid_ac = sum(ac_id[i:(self.win_size+i)]) / self.win_size * 100

                # Potential recombinant regions where % ident of A-C or B-C is higher than A-B
                recombinant_regions += "1" if pid_ac > pid_ab or pid_bc > pid_ab else "0"

            # 3. locate runs of 1's
            recomb_idx = [m.span() for m in re.finditer('1+', recombinant_regions)]

            # Convert coordinates from window-level to alignment-level and record number of windows
            coords = [(triplet.info_sites[x], triplet.info_sites[y-1]) for x, y in recomb_idx]
            
            for left, right in coords:
                N = right - left  # length of putative recombinant region
                if N <= 0:
                    continue
                
                # retrieve full sequences
                seq_a = triplet.sequences[0]
                seq_b = triplet.sequences[1]
                seq_c = triplet.sequences[2]
                
                M = 0  # number of nts in common between (A or B) and C in this region
                p = 0  # proportion of nts in common between (A or B) and C in the entire sequence
                L = len(seq_a)  # length of the entire sequence
                for i in range(N):
                    if seq_a[i]==seq_c[i] or seq_b[i]==seq_c[i]:
                        p += 1
                        if i >= left and i < right:
                            M += 1  # within putative recombinant region
                p /= L
                
                # Calculate p_value as binomial survival function (1-CDF)
                pvalue = math.exp(math.log(G) + math.log(L) - math.log(N) + 
                                  binom.logsf(M-1, n=N, p=p))  # RDP sums from M to N inclusive

                if pvalue != 0.0:
                    # FIXME: what's wrong with P=0?
                    rec_name, parents = identify_recombinant(triplet, [left, right])
                    self.raw_results.append((rec_name, parents, left, right, pvalue))
