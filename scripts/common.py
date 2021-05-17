import numpy as np
from itertools import combinations


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
