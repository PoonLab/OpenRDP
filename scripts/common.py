import numpy as np
from itertools import combinations
from scipy.stats import chi2_contingency


def generate_triplets(align):
    """
    Generate all possible combinations of sequence triplets
    :param align: a numpy character array of 'n' sequences
    :return: indices for every possible combination of sequence triplets
    """
    return combinations(range(align.shape[0]), 3)


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


class Triplet:

    def __init__(self, sequences):
        self.sequences = sequences
        self.poly_sites_align, self.poly_sites = self.remove_monomorphic_sites()
        self.info_sites_align, self.info_sites, self.uninfo_sites = self.remove_uninformative_sites()

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
            if len(self.poly_sites) > 1.5 * fixed_win_size:
                return 0.75 * len(self.poly_sites)

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
                return frac * num_var_sites
