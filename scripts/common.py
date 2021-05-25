import numpy as np
from itertools import combinations
from scipy.stats import chi2_contingency


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

            # Record polymorphic sites because it seems shorter than
            poly_sites.append(i)

    # Build "new alignment"
    new_align = align[:, poly_sites]

    return new_align, poly_sites


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
    # Compute chi-squared value if the sexpecetd frequencies are valid
    if (c_table[0][0] > 0 and c_table[0][1] > 0) or (c_table[1][0] > 0 and c_table[1][1] > 0):
        chi2, p_value, _, _ = chi2_contingency(c_table)

        # Record only significant events
        if p_value < max_pvalue:
            return chi2, p_value

    return None, None
