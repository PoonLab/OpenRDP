from itertools import combinations
import numpy as np
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor


class Sequence:
    """
    Represents an aligned sequence
    """
    def __init__(self, seq_name, sequence):
        self.sequence = np.array(list(sequence))
        self.seq_name = seq_name
        self.gap_pos = self.find_gaps()

    def gap_at_pos(self, pos):
        """
        Finds if there is a gap in the aligned sequence at a particular location
        :param pos: the position of interest
        :return: True if there is a gap, False otherwise
        """
        if self.sequence[pos] == '-':
            return True
        else:
            return False

    def find_gaps(self):
        """
        Finds position of gaps in the aligned sequence
        :return: a list of positions where there is a gap
        """
        return [i for i, val in enumerate(self.sequence) if val == '-']


class Alignment:
    """
    Represents an alignment of sequences
    """
    def __init__(self, aligned_seqs, circular=False):
        self.sequences = aligned_seqs    # List of references to sequences
        self.num_seqs = len(aligned_seqs)
        self.start_pos = 0
        self.end_pos = len(aligned_seqs[0].sequence)
        self.site_types = self.make_seq_cat_count()  # Analogous to SeqCatCount in RDP

    def make_seq_cat_count(self):
        """
        Categorize sites based on the presence of a gap and the number of variants
        (NucXX in Module4.bas)
        :return: list of tuples where:
                    - the first element represents if there is a gap
                    - the second element represents the number of nucleotide variants
        """

        # Create a numpy array of aligned sequences
        seqs = np.array([])
        for seq in self.sequences:
            np.append(seqs, np.array(seq.sequence), axis=0)

        # Classify the sites column by column
        unique, counts = np.unique(seqs, return_counts=True, axis=1)
        site_counts = dict(zip(unique, counts))

        sites = []
        for pos, d in site_counts:
            gaps = False
            if d[0] > 0:
                gaps = True
            num_variants = len(unique)
            sites.append((gaps, num_variants))

        return sites

# TODO: accessor functions for sites with different categories


def rdp_method(aln, win_size, automask):
    pass
    # Handle masked sequences and update the number of query and reference sequences

    # Calculate A, T, G, C

    # Write Initial file (Module 2)

    # Create analysis list best reference sequence

    # Set Up Globals

    # Calculate distances between sequences (Module 2)-- MakeCatCountSeq2P

    # UPGMA


    for triplet in combinations(aln, 3):
        print(triplet)

        # Remove unimformative sites
        subtriplet = remove_uninformative_sites(triplet)
        # ...


def remove_uninformative_sites(triplet):
    pass
