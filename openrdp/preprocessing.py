from itertools import combinations
import numpy as np
from scipy.spatial.distance import hamming


class Sequence:
    """
    Represents an aligned sequence
    """
    def __init__(self, seq_name, sequence, num):
        self.sequence = sequence
        self.length = len(self.sequence)
        self.seq_name = seq_name

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

    def encode_seq(self):
        """
        Encodes 2 nucleotides as 1 byte
        :return: the sequence encoded as a uint8 numpy array
        """
        # TODO: cleaner (numpy-less) way to do this?
        # Note: Using 4 bits for ease of implementation
        alpha_bin = {'A': [1, 1, 0, 0], 'T': [0, 0, 1, 1],
                     'G': [1, 0, 1, 0], 'C': [0, 1, 0, 1], '-': [0, 0, 0, 0]}

        encoded_seq = np.array([])

        # Extend length of the sequence to an even number for bit packing
        seq = self.sequence
        if len(self.sequence) % 2 != 0:
            np.concatenate((seq, ['-']), axis=None)

        i = 0
        # Read sequence 2 at a time
        while i < len(seq):
            byte_arr = np.concatenate(alpha_bin[str(seq[i])], alpha_bin[str(seq[i+1])])
            encoded_seq = np.append(byte_arr)
            i += 2

        bit_seq = np.packbits(encoded_seq)
        return bit_seq

    def remove_gaps(self):
        return str(self.sequence).replace('-', '')


class Alignment:
    """
    Represents an alignment of sequences
    """
    def __init__(self, aligned_seqs):
        self.sequences = aligned_seqs    # List of references to sequences
        self.num_seqs = len(aligned_seqs)
        self.start_pos = 0
        self.end_pos = len(aligned_seqs[0].sequence)
        # self.site_types = self.make_seq_cat_count()
        # self.pairwise_dists = self.pairwise_distances()

    def get_sequence(self, idx):
        return self.sequences[idx]

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

    def pairwise_distances(self):
        """
        Calculates the pairwise distance between the sequences using Hamming Distance and
        normalizes the distance based on the length of the sequence.
        The distance will be between 0 and 1 where 0 indicates no matches, and 1 indicates complete matches
        (see FastDistCalcZ from Module4.bas, MakeCompressor)
        :return: matrix of normalized pairwise distances
        """
        pairwise_dists = np.array([self.num_seqs, self.num_seqs])

        # Calculate pairwise hamming distance, normalized by the sequence length
        len_seq = len(self.sequences[1].sequence)
        for seq1 in self.sequences:
            for seq2 in self.sequences:
                if seq1 != seq2:
                    h_dist = hamming(seq1.sequence, seq2.sequence)
                    norm_h_dist = (1 - h_dist) / len_seq
                    # h_dist = seq1.bit_seq ^ seq2.bit_seq
                    # matches = 0
                    # while h_dist > 0:
                    #     matches += 1
                    #     h_dist >>= 1
                    # norm_h_dist = h_dist / len_seq

                    pairwise_dists[seq1.num][seq2.num] = norm_h_dist

        return pairwise_dists
