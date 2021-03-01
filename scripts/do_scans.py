import os
import sys
import logging
import subprocess
import random
import numpy as np
from itertools import combinations
from scipy.stats import chi2_contingency


class RdpMethod:
    """
    Executes command for RDP method
    """
    def __init__(self, win_size=30, reference=None):
        self._win_size = win_size
        self.reference = reference

    def set_options_from_config(self, settings={}):
        pass

    def validate_options(self):
        pass

    def execute(self, align):
        pass


class GeneConv:
    def __init__(self, gscale=1, ignore_indels=False, min_length=1, min_poly=2, min_score=2, max_overlap=1):
        self.gscale = gscale                # Mismatch penalty
        self.ignore_indels = ignore_indels  # Ignore indels or treat indels as one polymorphism (default = False)
        self.min_length = min_length
        self.min_poly = min_poly
        self.min_score = min_score
        self.max_overlap = max_overlap

    def set_options_from_config(self, settings={}):
        pass

    def validate_options(self):
        pass

    def execute(self, data_path):

        # Check if valid options
        if not self.validate_options():
            return None

        # Create config file
        with open("geneconv.cfg", 'w+') as cfg_handle:
            cfg_handle.write('#GCONV_CONFIG\n')
            cfg_handle.write('  -inputpath={}\n'.format(os.path.realpath(data_path.name)))

            if not self.ignore_indels:
                cfg_handle.write('  -Indel_blocs\n')

            cfg_handle.write('  -Gscale={}\n'.format(self.gscale))
            cfg_handle.write('  -Minlength={}\n'.format(self.min_length))
            cfg_handle.write('  -Minnpoly={}\n'.format(self.min_poly))
            cfg_handle.write('  -Minscore={}\n'.format(self.min_score))
            cfg_handle.write('  -Maxoverlap={}\n'.format(self.max_overlap))

            cfg_path = format(os.path.realpath(cfg_handle.name))

        # Path to GENECONV executables
        script_path = os.path.dirname(os.path.abspath(__file__))

        if sys.platform.startswith("win"):
            bin_path = os.path.join(script_path, 'bin/GENECONV/windows_geneconv.exe')
        else:
            bin_path = os.path.join(script_path, 'bin/GENECONV/unix_geneconv.exe')

        if not os.path.isfile(bin_path):
            logging.error("No file exists")

        # Run GENECONV
        if sys.platform.startswith("win"):
            gc_output = subprocess.check_output([bin_path, cfg_path], shell=False, stderr=subprocess.DEVNULL)
        else:
            gc_output = subprocess.check_output([bin_path, cfg_path], shell=False, stderr=subprocess.STDOUT)

        # Remove file
        os.remove(cfg_path)

        return gc_output


class ThreeSeq:
    def __init__(self, p_value_table=None):
        self.p_value_table = p_value_table

    def load_config(self):
        pass

    @staticmethod
    def execute(data_path):
        # Path to 3Seq executables
        script_path = os.path.dirname(os.path.abspath(__file__))

        if sys.platform.startswith("win"):
            bin_path = os.path.join(script_path, 'bin/3Seq/windows_3seq.exe')
        else:
            bin_path = os.path.join(script_path, 'bin/GENECONV/unix_3seq.exe')

        if not os.path.isfile(bin_path):
            logging.error("No file exists")

        # Run 3Seq
        if sys.platform.startswith("win"):
            tseq_output = subprocess.check_output([bin_path, '-f', os.path.realpath(data_path.name)],
                                                  shell=False, stderr=subprocess.DEVNULL)
        else:
            tseq_output = subprocess.check_output([bin_path, '-f', os.path.realpath(data_path.name)],
                                                  shell=False, stderr=subprocess.STDOUT)

        return tseq_output


class Siscan:
    def __init__(self, win_size=200, step_size=20, strip_gaps=True, use_pos=3, pvalue_perm_num=1100, scan_perm_num=100,
                 random_seed=3, fourth_seq_sel='random'):
        self.win_size = win_size
        self.step_size = step_size
        self.strip_gaps = strip_gaps
        self.use_pos = use_pos
        self.pvalue_perm_num = pvalue_perm_num
        self.scan_perm_num = scan_perm_num
        self.random_seed = random_seed
        self.fourth_seq = fourth_seq_sel

    def set_options_from_config(self, settings={}):
        pass

    def validate_options(self):
        pass

    def count_patterns(self, a, b, c, d):
        """
        Count the number of positions within the window that conform to each of the 15 patterns
        :param a: the first sequence in the triplet
        :param b: the second sequence in the triplet
        :param c: the third sequence in the triplet
        :param d: the fourth sequence
        :return: counts for each pattern
        """
        pat_counts = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

        # TODO: simplify logic?

        for pos in range(self.win_size):
            # Pattern 2
            if a[pos] == b[pos] and a[pos] != c[pos] and a[pos] != d[pos]:
                pat_counts[1] += 1
            # Pattern 3
            elif a[pos] == c[pos] and a[pos] != b[pos] and a[pos] != d[pos]:
                pat_counts[2] += 1
            # Pattern 4
            elif a[pos] == d[pos] and a[pos] != b[pos] and a[pos] != c[pos]:
                pat_counts[3] += 1
            # Pattern 5
            elif b[pos] == c[pos] and b[pos] != a[pos] and b[pos] != d[pos]:
                pat_counts[4] += 1
            # Pattern 6
            elif b[pos] == d[pos] and b[pos] != a[pos] and b[pos] != c[pos]:
                pat_counts[5] += 1
            # Pattern 7
            elif c[pos] == d[pos] and c[pos] != b[pos] and c[pos] != a[pos]:
                pat_counts[6] += 1
            # Pattern 8
            elif a[pos] == b[pos] and c[pos] == d[pos] and b[pos] != c[pos]:
                pat_counts[7] += 1
            # Pattern 9
            elif a[pos] == c[pos] and b[pos] == d[pos] and b[pos] != c[pos]:
                pat_counts[8] += 1
            # Pattern 10
            elif a[pos] == d[pos] and b[pos] == c[pos] and a[pos] != b[pos]:
                pat_counts[9] += 1
            # Pattern 11
            elif a[pos] == b[pos] and b[pos] == c[pos] and a[pos] != d[pos]:
                pat_counts[10] += 1
            # Pattern 12
            elif a[pos] == b[pos] and b[pos] == d[pos] and a[pos] != c[pos]:
                pat_counts[11] += 1
            # Pattern 13
            elif a[pos] == c[pos] and c[pos] == d[pos] and a[pos] != b[pos]:
                pat_counts[12] += 1
            # Pattern 14
            elif b[pos] == c[pos] and c[pos] == d[pos] and b[pos] != a[pos]:
                pat_counts[13] += 1
            # Pattern 15
            elif a[pos] == b[pos] and b[pos] == c[pos] and c[pos] == d[pos]:
                pat_counts[14] += 1
            # Pattern 1
            else:
                pat_counts[0] += 1

        return pat_counts

    @staticmethod
    def sum_pattern_counts(pat_counts):
        """
        Sum counts where 2 sequences are identical
        :param pat_counts: vector of counts corresponding to each pattern
        :return: vector of summed counts
        """
        sum_pat_counts = [0, 0, 0, 0, 0, 0, 0, 0, 0]

        sum_pat_counts[3] = pat_counts[1] + pat_counts[7] + pat_counts[10] + pat_counts[11]  # 2 + 8 + 11 + 12
        sum_pat_counts[4] = pat_counts[2] + pat_counts[8] + pat_counts[10] + pat_counts[12]  # 3 + 9 + 11 + 13
        sum_pat_counts[5] = pat_counts[3] + pat_counts[9] + pat_counts[11] + pat_counts[12]  # 4 + 10 + 12 + 13
        sum_pat_counts[6] = pat_counts[4] + pat_counts[9] + pat_counts[10] + pat_counts[13]  # 5 + 10 + 11 + 14
        sum_pat_counts[7] = pat_counts[5] + pat_counts[8] + pat_counts[11] + pat_counts[13]  # 6 + 9 + 12 + 14
        sum_pat_counts[8] = pat_counts[6] + pat_counts[7] + pat_counts[12] + pat_counts[13]  # 7 + 8 + 13 + 14

        return sum_pat_counts

    def execute(self, alignment):
        """
        Do Sister-scanning as described in Gibbs, Armstrong, and Gibbs (2000)

        """
        random.seed(self.random_seed)

        # Remove gaps from sequences
        if self.strip_gaps:
            ungapped = []
            for s in alignment.sequences:
                ungapped.append(s.ungapped)

        # Generate all pairs of triplet sequences
        for triplet in list(combinations(alignment.sequences, 3)):

            # TODO: pad sequences if not divisible by window length
            # Based on leading edge of the window
            for window in range(len(alignment.sequences.length), self.step_size):
                win_start = 0
                win_end = win_start + self.win_size

                # Label sequences
                a = triplet[0][win_start: win_end]
                b = triplet[1][win_start: win_end]
                c = triplet[2][win_start: win_end]

                # Create the fourth sequence through horizontal randomization
                selected_seq = random.choice(('a', 'b', 'c'))
                d = np.array([selected_seq])[win_start: win_end]
                random.shuffle(d)

                # (1) Count number of positions within a window that conform to each pattern
                pat_counts = self.count_patterns(a, b, c, d)

                # (2) Sum counts where 2 sequences are identical
                sum_pat_counts = self.sum_pattern_counts(pat_counts)

                # (3) Sum counts of each kind of informative site for each window
                sum_pat_counts[0] = pat_counts[1] + pat_counts[6] + pat_counts[7]  # 2 + 7 + 8
                sum_pat_counts[1] = pat_counts[2] + pat_counts[5] + pat_counts[8]  # 3 + 6 + 9
                sum_pat_counts[2] = pat_counts[3] + pat_counts[4] + pat_counts[9]  # 4 + 5 + 10

                # (4) Create 4 vertically randomized sequences (steps 1 and 2), repeat for 100 times
                p_counts = []
                sum_p_counts = []
                for i in range(100):

                    seq_array = np.array([a, b, c, d], axis=0)
                    # Generate 4 vertically randomized sequences (shuffle the columns)
                    a1 = seq_array.shuffle(axis=0)
                    b1 = seq_array.shuffle(axis=0)
                    c1 = seq_array.shuffle(axis=0)
                    d1 = seq_array.shuffle(axis=0)

                    # Count number of patterns and sum counts
                    p_counts = self.count_patterns(a1, b1, c1, d1)
                    sum_p_counts = self.sum_pattern_counts(p_counts)

                # (5) Calculate Z-scores for each pattern and sum of patterns for each window
                pop_mean_pcounts = np.mean(p_counts)
                pop_mean_patsum = np.mean(sum_p_counts)
                pop_std_pcounts = np.std(p_counts)
                pop_std_patsum = np.std(sum_p_counts)

                pat_zscore = []
                for val in pat_counts:
                    z = (val - pop_mean_pcounts) / pop_std_pcounts
                    pat_zscore.append(z)

                sum_pat_zscore = []
                for val in sum_pat_counts:
                    z = (val - pop_mean_patsum) / pop_std_patsum
                    sum_pat_zscore.append(z)

        return pat_zscore, sum_pat_zscore


class MaxChi:
    def __init__(self, win_size=200, strip_gaps=True):
        self.win_size = win_size
        self.strip_gaps = strip_gaps

    def set_options_from_config(self, settings={}):
        pass

    def validate_options(self):
        pass

    def execute(self, align, triplets):
        """
        Executes the MaxChi algorithm
        :param align: a n x m numpy array where n is the length of the alignment and m is the number of sequences
        :param triplets: list of all combinations of
        :return:
        """

        # 1. Remove monomorphic sites
        poly_sites = []
        # Find positions of polymorphic sites
        for i in range(len(align)):
            col = align[:, i]
            if not np.all(col == col[0]):
                poly_sites.append(i)

        # Build "new alignment"
        new_align = align[:, poly_sites]

        maxchi_out = []

        # 2. Select 3 sequences
        for a, b, c in triplets:
            trps = [new_align[a], new_align[b], new_align[c]]

            # 3. Sample two sequences
            pairs = [(0, 1), (1, 2), (2, 0)]
            for i, j in pairs:
                seq1 = trps[i]
                seq2 = trps[j]

                # Slide along the sequences
                for k in range(len(new_align)):
                    reg1_r = seq1[k: self.win_size/2]
                    reg2_r = seq2[k: self.win_size/2]
                    reg1_l = seq1[self.win_size/2: self.win_size]
                    reg2_l = seq2[self.win_size/2: self.win_size]

                    c_table = [[0, 0],
                               [0, 0]]

                    # Compute contingency table for each window position
                    left_matches = np.sum((reg1_l == reg2_l))
                    print(left_matches)
                    c_table[0][0] = int(left_matches)
                    c_table[0][1] = int(self.win_size / 2) - left_matches

                    right_matches = np.sum((reg1_r == reg2_r))
                    print(right_matches)
                    c_table[1][0] = int(right_matches)
                    c_table[1][1] = int(self.win_size / 2) - right_matches

                    # Compute chi-squared value
                    chi2, p_value, _, _ = chi2_contingency(c_table)

                    maxchi_out.append((chi2, p_value))

        return maxchi_out


class Chimaera:
    def __init__(self, win_size=200, strip_gaps=True):
        self.win_size = win_size
        self.strip_gaps =strip_gaps

    def set_options_from_config(self, settings={}):
        pass

    def validate_options(self):
        pass

    def execute(self, align, triplets):

        # 1. Remove monomorphic sites
        poly_sites = []
        # Find positions of polymorphic sites
        for i in range(len(align)):
            col = align[:, i]
            if not np.all(col == col[0]):
                poly_sites.append(i)

        # Build "new alignment"
        new_align = align[:, poly_sites]

        chimaera_out = []

        for a, b, c in triplets:
            trps = [new_align[a], new_align[b], new_align[c]]

            combos = [(0, 1, 2), (1, 2, 0), (2, 1, 0)]
            for run in combos:
                recombinant = trps[run[0]]
                parental_1 = trps[run[1]]
                parental_2 = trps[run[2]]

                # Remove sites where neither the parental match the recombinant
                informative_sites = []
                aln = np.array([recombinant, parental_1, parental_2])
                # Find indices for informative sites
                for i in range(len(aln)):
                    col = align[:, i]
                    if len(np.unique(col)) != 3:
                        informative_sites += i

                # Build new alignment with uninformative sites removed
                new_aln = aln[:, informative_sites]

                # 2. Compress "recombinant" into bitstrings
                comp_seq = []
                for i in range(len(new_align)):
                    if new_align[i, 0] == new_align[i, 1]:
                        comp_seq.append(0)
                    elif new_align[i, 0] == new_align[i, 2]:
                        comp_seq.append(1)

                # Move sliding window along compressed sequence , 1 position at a time
                # Slide along the sequences
                for k in range(len(comp_seq)):
                    reg_r = comp_seq[k: self.win_size/2]
                    reg_l = comp_seq[self.win_size/2: self.win_size]

                    c_table = [[0, 0],
                               [0, 0]]

                    # Compute contingency table for each window position
                    count_left_ones = np.count_nonzero(reg_l)
                    c_table[0][0] = count_left_ones
                    c_table[0][1] = self.win_size/2 - count_left_ones

                    count_right_ones = np.count_nonzero(reg_r)
                    c_table[1][0] = count_right_ones
                    c_table[1][1] = self.win_size/2 - count_right_ones

                    # Compute chi-squared value
                    chi2, p_value, _, _ = chi2_contingency(c_table)

                    chimaera_out.append((chi2, p_value))
