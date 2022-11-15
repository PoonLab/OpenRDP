from itertools import combinations

import numpy as np
from scipy.stats import chi2_contingency
from scipy.stats import pearsonr
import copy


def generate_triplets(align):
    """
    Generate all possible combinations of sequence triplets
    :param align: a numpy character array of 'n' sequences
    :return: indices for every possible combination of sequence triplets
    """
    return combinations(range(align.shape[0]), 3)

def generate_combination(align, n):
    """
    Generate all possible combinations of sequence triplets
    :param align: a numpy character array of 'n' sequences
    :return: indices for every possible combination of sequence triplets
    """
    return combinations(range(align.shape[0]), n)


def calculate_chi2(c_table, max_pvalue):
    """
    Computes the chi-squared value and returns the chi-squared and p-value if the difference is significant
    :param c_table: a 2x2 contingency table
    :param max_pvalue: the p-value threshold
    :return: a tuple containing the chi-squared value and p-value
    """
    # Compute chi-squared value if the expected frequencies are valid
    if (c_table[0][0] > 0 and c_table[0][1] > 0) or (c_table[1][0] > 0 and c_table[1][1] > 0):

        chi2, p_value, _, _ = chi2_contingency(c_table)

        # Record only significant events
        if p_value < max_pvalue:
            return chi2, p_value

    return None, None


def reduce_to_unique_seqs(aln):
    """
    Remove identical sequences
    :param aln: list of aligned sequences
    :return: a list of unique sequences
    """
    return list(set(aln))


def percent_diff(s1, s2):
    s1_valid = (s1 == 'A') | (s1 == 'T') | (s1 == 'G') | (s1 == 'C')
    s2_valid = (s2 == 'A') | (s2 == 'T') | (s2 == 'G') | (s2 == 'C')
    valid = s1_valid & s2_valid
    diffs = np.sum(s1[valid] != s2[valid])
    num_valid = valid.sum()
    return diffs.sum() / num_valid if num_valid else 0


def jc_distance(s1, s2):
    """
    Calculate the pairwise Jukes-Cantor distance between 2 sequences
    :param s1: the first sequence
    :param s2: the second sequence
    :return: the pairwise JC69 distance between 2 sequences
    """
    p_dist = percent_diff(s1, s2)
    if p_dist >= 0.75:
        return 1
    return -0.75 * np.log(1 - (p_dist * 4 / 3)) if p_dist else 0


def all_items_equal(x):
    """
    Check if all items in a list are identical
    :param x: the list
    :return: True if all the items are identical, false otherwise
    """
    return x.count(x[0]) == len(x)


def identify_recombinant(trp, aln_pos):
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


class Triplet:

    def __init__(self, alignment, seq_names, trp_idxs):
        self.idxs = trp_idxs
        self.sequences = self.get_triplets(alignment)
        self.names = self.get_trp_names(seq_names)
        self.poly_sites_align, self.poly_sites = self.remove_monomorphic_sites()
        self.info_sites_align, self.info_sites, self.uninfo_sites = self.remove_uninformative_sites()

    def get_triplets(self, alignment):
        return np.take(alignment, self.idxs, axis=0)

    def get_trp_names(self, seq_names):
        return [seq_names[i] for i in self.idxs]

    def get_sequence_name(self, trp_idx):
        return self.names[trp_idx]

    def get_seq_from_name(self, trp_name):
        for i, name in enumerate(self.names):
            if trp_name == name:
                return self.sequences[i]

    def remove_uninformative_sites(self):
        """
        Remove sites that are all the same or all different
        """
        infor_sites = []
        uninfor_sites = []
        # Find positions of sites that are all the same sites or all sites that are different
        for i in range(self.sequences.shape[1]):
            col = self.sequences[:, i]
            if np.unique(col).shape[0] == 2:
                infor_sites.append(i)
            else:
                uninfor_sites.append(i)

        # Build "new alignment"
        new_aln = self.sequences[:, infor_sites]

        return new_aln, infor_sites, uninfor_sites

    def remove_monomorphic_sites(self):
        """
        Remove monomorphic sites
        :return: a tuple containing the polymorphic sites and the positions of polymorphic sites
        """
        poly_sites = []
        # Find positions of polymorphic sites
        for i in range(self.sequences.shape[1]):
            col = self.sequences[:, i]
            if not np.all(col == col[0]):
                # Record polymorphic sites because it seems shorter than
                poly_sites.append(i)

        # Build "new alignment"
        new_align = self.sequences[:, poly_sites]

        return new_align, poly_sites

    def get_win_size(self, offset, win_size, fixed_win_size, num_var_sites, frac_var_sites):
        """
        Sets the size (can be fixed or variable) of the sliding window
        :param offset: the step size
        """
        # Fixed window length
        if fixed_win_size:
            return win_size

        # User specified number of variable sites
        else:
            if len(self.poly_sites) > 1.5 * win_size:
                return round(0.75 * len(self.poly_sites))

            elif num_var_sites:
                curr = 0
                while curr < num_var_sites - 1:
                    curr += 1
                return self.poly_sites[curr] + offset

            # User specified proportion of variable sites
            else:
                num_var_sites = 1
                frac = 0
                while frac < frac_var_sites:
                    for var_site in self.poly_sites:
                        frac = num_var_sites / var_site
                        num_var_sites += 1
                return round(frac * num_var_sites)
