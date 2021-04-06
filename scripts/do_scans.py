import glob
import logging
import os
import random
import re
import subprocess
import sys
import copy

import numpy as np
from scipy.signal import find_peaks
from scipy.stats import chi2_contingency


# TODO: multiple comparison correction


class ThreeSeq:
    def __init__(self, in_path):
        self.in_path = os.path.realpath(in_path)
        self.in_name = os.path.basename(in_path)

    def execute(self):
        """
        Execute the 3Seq algorithm.
            Lam HM, Ratmann O, Boni MF.
            Improved algorithmic complexity for the 3SEQ recombination detection algorithm.
            Mol Biol Evol, 35(1):247-251, 2018.
        :return: A list containing the results of the 3Seq analysis
                    Format: [triplets, uncorrected p-value, corrected p-value, breakpoint locations]
        """
        # Clear output files
        out_files = glob.glob('./*.3s.*')
        for f in out_files:
            try:
                os.remove(f)
            except OSError:
                pass

        # Set paths to 3Seq executables
        if sys.platform.startswith("win"):
            bin_path = os.path.abspath('../utils/bin/3Seq/windows_3seq.exe')
        else:
            bin_path = os.path.abspath('../utils/bin/3Seq/unix_3seq.exe')

        if not os.path.isfile(bin_path):
            logging.error("No 3Seq executable file exists.")

        # Run 3Seq
        if sys.platform.startswith("win"):
            try:
                subprocess.check_output([bin_path, "-f", self.in_path, "-d", "-id", self.in_name],
                                        shell=False, input=b"Y\n")  # Respond to prompt
            except subprocess.CalledProcessError as e:
                print(e.output, e.returncode)
        else:
            try:
                subprocess.check_output([bin_path, "-f", self.in_path, "-d", "-id", self.in_name],
                                        shell=False, input=b"Y\n")  # Respond to prompt
            except subprocess.CalledProcessError as e:
                print(e.output, e.returncode)

        # Parse the output of 3Seq
        tseq_out = self.in_name + '.3s.rec'
        out_path = os.path.join(os.getcwd(), tseq_out)
        ts_results = self.parse_output(out_path)

        return ts_results

    def parse_output(self, out_path):
        """
        Parse the output of the 3Seq analysis
        :param out_path: Path to the output file containing information about recombinant sequences
        :return: List of triplets, corrected and uncorrected p-values, and breakpoint locations
        """
        ts_results = []

        # Check that the out file exists
        try:
            with open(out_path) as out_handle:
                out_handle.readline()  # Read first line

                for line in out_handle:
                    line = line.split('\t')
                    line = [l.strip() for l in line]
                    triplet = ([line[0], line[1], line[2]])  # Record the triplets
                    uncorr_p_value = line[6]  # Uncorrected p-value
                    corr_p_value = line[10]  # Dunn-Sidak corrected p-value
                    locations = line[12:]  # Breakpoint locations
                    ts_results.append([triplet, uncorr_p_value, corr_p_value, locations])
        except FileNotFoundError as e:
            print(e)

        return ts_results


class GeneConv:
    def __init__(self, settings=None, gscale=1, ignore_indels=False, min_length=1, min_poly=2, min_score=2,
                 max_overlap=1):
        """
        Constructs a GeneConv object
        :param gscale: mismatch penalty
        :param ignore_indels: Ignore indels or treat indels as one polymorphism (default = False)
        :param min_length: Minimum length of the fragments
        :param min_poly: Minimum number of polymorphic sites
        :param min_score: Minimum pairwise score
        :param max_overlap: Maximum number of overlapping pairs
        """
        if settings is not None:
            self.set_options_from_config(settings)
            self.validate_options()

        else:
            self.gscale = gscale
            self.ignore_indels = ignore_indels
            self.min_length = min_length
            self.min_poly = min_poly
            self.min_score = min_score
            self.max_overlap = max_overlap

    def set_options_from_config(self, settings):
        """
        Set the parameters of GENECONV from the config file
        :param settings: a dictionary of settings
        """
        self.gscale = int(settings['mismatch_penalty'])

        if settings['indels_as_polymorphisms'] == 'True':
            self.ignore_indels = False
        elif settings['indels_as_polymorphisms'] == 'False':
            self.ignore_indels = True
        else:
            self.ignore_indels = None

        self.min_length = int(settings['min_len'])
        self.min_poly = int(settings['min_poly'])
        self.min_score = int(settings['min_score'])
        self.max_overlap = int(settings['max_num'])

    def validate_options(self):
        """
        Check if the options from the config file are valid
        If the options are invalid, the default value will be used instead
        """
        if not isinstance(self.ignore_indels, bool):
            print("Invalid option for 'indels_as_polymorphisms'.\nUsing default value (False) instead.\n")
            self.ignore_indels = False

        if self.min_length <= 0:
            print("Invalid option for 'min_len'.\nUsing default value (1) instead.")
            self.min_length = 1

        if self.min_poly < 1:
            print("Invalid option for 'min_score'.\nUsing default value (2) instead.")
            self.min_poly = 2

        if self.max_overlap < 0:
            print("Invalid option for 'max_num'.\nUsing default value (1) instead.")
            self.min_score = 1

    def execute(self, in_path):
        """
        Execute the GENECONV algorithm
            S. A. Sawyer (1999)
            GENECONV: A computer package for the statistical detection of gene conversion.
            Distributed by the author, Department of Mathematics, Washington University in St. Louis,
            Available at http://www.math.wustl.edu/~sawyer.
        :param in_path: Path to the input alignment file
        :return: A list of results
        """
        # Clear output files
        out_files = glob.glob('../utils/bin/GENECONV/*.frags') + glob.glob('../utils/bin/GENECONV/*.sum')
        for f in out_files:
            try:
                os.remove(f)
            except OSError:
                pass

        # Create config file
        with open("geneconv.cfg", 'w+') as cfg_handle:
            cfg_handle.write('#GCONV_CONFIG\n')
            cfg_handle.write('  -inputpath={}\n'.format(os.path.realpath(in_path)))

            if not self.ignore_indels:
                cfg_handle.write('  -Indel_blocs\n')

            cfg_handle.write('  -Gscale={}\n'.format(self.gscale))
            cfg_handle.write('  -Minlength={}\n'.format(self.min_length))
            cfg_handle.write('  -Minnpoly={}\n'.format(self.min_poly))
            cfg_handle.write('  -Minscore={}\n'.format(self.min_score))
            cfg_handle.write('  -Maxoverlap={}\n'.format(self.max_overlap))

            cfg_path = format(os.path.realpath("geneconv.cfg"))

        # Path to GENECONV executables
        if sys.platform.startswith("win"):
            bin_path = os.path.abspath('../utils/bin/GENECONV/windows_geneconv.exe')
        else:
            bin_path = os.path.abspath('../utils/bin/GENECONV/unix_geneconv.exe')

        if not os.path.isfile(bin_path):
            logging.error("No GENECONV executable file exists")

        # Run GENECONV
        if sys.platform.startswith("win"):
            try:
                subprocess.check_output([bin_path, "-Seqfile={}".format(in_path),
                                         "-Config={}".format(cfg_path), "-nolog"], shell=False)
            except subprocess.CalledProcessError as e:
                print(e.output, e.returncode)
        else:
            try:
                subprocess.check_output([bin_path, "-Seqfile={}".format(in_path),
                                         "-Config={}".format(cfg_path), "-nolog"], shell=False)
            except subprocess.CalledProcessError as e:
                print(e.output, e.returncode)

        # Remove GENECONV config file
        try:
            os.remove(cfg_path)
        except OSError:
            pass

        # Parse the output of 3Seq
        in_name = os.path.basename(in_path).split('.')[0]
        out_name = in_name + '.frags'
        out_path = os.path.join(os.path.dirname(in_path), out_name)
        gc_results = self.parse_output(out_path)

        return gc_results

    def parse_output(self, out_path):
        """
        Parse the output of the GENECONV analysis
        :param out_path: Path to the output file
        :return: List of results
        """
        gc_results = []

        # Check that the out file exists
        try:
            with open(out_path) as out_handle:

                for line in out_handle:
                    if not line.startswith('#'):
                        line = line.strip()
                        line = line.split()
                        seqs = line[1]  # Sequence pairs
                        uncorr_p_value = line[2]  # Uncorrected p_value
                        corr_p_value = line[3]  # Bonferroni Corrected - Karlin-Altschul
                        locations = (line[4], line[5])  # Locations in alignment
                        type = line[0]  # Global inner (GI), global outer (GO), additional inner (AI)
                        gc_results.append([seqs, uncorr_p_value, corr_p_value, locations, type])

        except FileNotFoundError as e:
            print(e)

        return gc_results


class MaxChi:
    def __init__(self, align, settings=None, win_size=200, strip_gaps=True, fixed_win_size=True, num_var_sites=None,
                 frac_var_sites=None):
        """
        Constructs a MaxChi Object
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

        self.align = align
        self.new_align, self.poly_sites = self.remove_monomorphic_sites(align)
        self.results = []

    def set_options_from_config(self, settings):
        """
        Set the parameters of MaxChi from the config file
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
                print("Invalid option for 'num_var_sites'.\nUsing default value (70) instead.")
                self.num_var_sites = 70

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
        Executes the MaxChi algorithm
        :param triplet: a list of three sequences
        """

        # Get the triplet sequences
        trp_seqs = []
        for seq_num in triplet:
            trp_seqs.append(self.new_align[seq_num])

        # 2. Sample two sequences
        pairs = [(0, 1), (1, 2), (2, 0)]
        for i, j in pairs:
            seq1 = trp_seqs[i]
            seq2 = trp_seqs[j]

            # Slide along the sequences
            half_win_size = int(self.win_size // 2)
            for k in range(self.new_align.shape[1] - self.win_size):
                reg1_left = seq1[k: half_win_size + k]
                reg2_left = seq2[k: half_win_size + k]
                reg1_right = seq1[k + half_win_size: k + self.win_size]
                reg2_right = seq2[k + half_win_size: k + self.win_size]

                c_table = [[0, 0],
                           [0, 0]]

                # Compute contingency table for each window position
                r_matches = np.sum((reg1_right == reg2_right))
                c_table[0][0] = int(r_matches)
                c_table[0][1] = half_win_size - r_matches

                l_matches = np.sum((reg1_left == reg2_left))
                c_table[1][0] = int(l_matches)
                c_table[1][1] = half_win_size - l_matches

                # Compute chi-squared value
                chi2, p_value, _, _ = chi2_contingency(c_table)

                self.results.append((triplet, (i, j), chi2, p_value))

        return


class Chimaera:
    def __init__(self, align, settings=None, win_size=200, strip_gaps=True, fixed_win_size=True, num_var_sites=None,
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
                chi2, p_value, _, _ = chi2_contingency(c_table)

                self.results.append((triplet, run, chi2, p_value))

        return


class Siscan:
    def __init__(self, align, settings=None, win_size=200, step_size=20, strip_gaps=True, pvalue_perm_num=1100,
                 scan_perm_num=100, random_seed=3):
        """
        Constructs a Siscan object
        :param win_size: the size of the sliding window
        :param step_size: the step size of the sliding window
        :param strip_gaps: Strip gaps or keep gaps
        :param pvalue_perm_num: p-value permutation number
        :param scan_perm_num: number of permutations of the scans
        :param random_seed: the random seed
        """
        self.align = align
        if settings:
            self.set_options_from_config(settings)
            self.validate_options(align)
        else:
            self.win_size = win_size
            self.step_size = step_size
            self.strip_gaps = strip_gaps
            self.pvalue_perm_num = pvalue_perm_num
            self.scan_perm_num = scan_perm_num
            self.random_seed = random_seed

        self.results = []

    def set_options_from_config(self, settings):
        """
        Set the parameters of Siscan from the  config file
        :param settings: a dictionary of settings
        """
        self.win_size = int(settings['win_size'])
        self.step_size = int(settings['step_size'])

        if settings['strip_gaps'] == 'True':
            self.strip_gaps = True
        else:
            self.strip_gaps = False

        self.pvalue_perm_num = int(settings['pvalue_perm_num'])
        self.scan_perm_num = int(settings['scan_perm_num'])
        self.random_seed = int(settings['random_seed'])

    def validate_options(self, align):
        """
        Check if the options from the config file are valid
        If the options are invalid, the default value will be used instead
        """
        if self.win_size < 0 or self.win_size > align.shape[1]:
            print("Invalid option for 'win_size'.\nUsing default value (200) instead.")
            self.win_size = 200

        if self.step_size < 0 or self.step_size >= self.win_size or self.step_size > align.shape[1]:
            print("Invalid option for 'step_size'.\nUsing default value (20) instead.")
            self.step_size = 20

        if self.pvalue_perm_num <= 0:
            print("Invalid option for 'pvalue_perm_num'.\nUsing default value (1100) instead.")
            self.pvalue_perm_num = 1100

        if self.scan_perm_num <= 0:
            print("Invalid option for 'scan_perm_num'.\nUsing default value (100) instead.")
            self.scan_perm_num = 100

        if self.random_seed <= 0:
            print("Invalid option for 'random_seed'.\nUsing default value (3) instead.")
            self.random_seed = 3

    def count_patterns(self, s_array):
        """
        Count the number of positions within the window that conform to each of the 15 patterns
        :param a: the first sequence in the triplet
        :param b: the second sequence in the triplet
        :param c: the third sequence in the triplet
        :param d: the fourth sequence
        :return: counts for each pattern
        """
        a = s_array[0]
        b = s_array[1]
        c = s_array[2]
        d = s_array[3]

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

    @staticmethod
    def vertical_randomization(seq_array):
        """
        Vertically randomizes a sequence
        Adapted from https://stackoverflow.com/questions/27486677/best-way-to-permute-contents-of-each-column-in-numpy
        :param seq_array: n x m numpy array of sequences (n = number of seqeunces, m = length of seqeunces)
        :return: a vertically randomized version of the sequence
        """
        s_array = copy.deepcopy(seq_array)
        ix_i = np.random.sample(s_array.shape).argsort(axis=0)
        ix_j = np.tile(np.arange(s_array.shape[1]), (s_array.shape[0], 1))
        v_random_seq = s_array[ix_i, ix_j]
        return v_random_seq

    def execute(self, alignment, triplet):
        """
        Do Sister-scanning as described in Gibbs, Armstrong, and Gibbs (2000)
        """
        random.seed(self.random_seed)

        # Get the triplet sequences
        trp_seqs = []
        for seq_num in triplet:
            trp_seqs.append(alignment[seq_num])

        # TODO: pad sequences if not divisible by window length
        # Based on leading edge of the window
        for window in range(0, alignment.shape[1], self.step_size):
            win_start = 0
            win_end = win_start + self.win_size

            # Label sequences
            a = trp_seqs[0][win_start: win_end]
            b = trp_seqs[1][win_start: win_end]
            c = trp_seqs[2][win_start: win_end]

            # Create the fourth sequence through horizontal randomization
            selected_seq = random.choice((0, 1, 2))
            d = trp_seqs[selected_seq][win_start: win_end]
            np.random.shuffle(d)

            # (1) Count number of positions within a window that conform to each pattern
            seq_array = np.array([a, b, c, d])
            pat_counts = self.count_patterns(seq_array)

            # (2) Sum counts where 2 sequences are identical
            sum_pat_counts = self.sum_pattern_counts(pat_counts)

            # (3) Sum counts of each kind of informative site for each window
            sum_pat_counts[0] = pat_counts[1] + pat_counts[6] + pat_counts[7]  # 2 + 7 + 8
            sum_pat_counts[1] = pat_counts[2] + pat_counts[5] + pat_counts[8]  # 3 + 6 + 9
            sum_pat_counts[2] = pat_counts[3] + pat_counts[4] + pat_counts[9]  # 4 + 5 + 10

            # (4) Create 4 vertically randomized sequences (steps 1 and 2), repeat for 100 times
            p_counts = []
            sum_p_counts = []
            for i in range(self.scan_perm_num):

                # Generate 4 vertically randomized sequences (shuffle the values in the columns)
                a1 = self.vertical_randomization(seq_array)

                # Count number of patterns and sum counts
                p_counts = self.count_patterns(a1)
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

            self.results.append((triplet, window, pat_zscore, sum_pat_zscore))

        return


class RdpMethod:
    """
    Executes RDP method
    """
    def __init__(self, align, settings=None, win_size=30, reference=None, min_id=0, max_id=100):
        if settings:
            self.set_options_from_config(settings)
            self.validate_options()

        else:
            self.win_size = win_size
            self.reference = reference
            self.min_id = min_id
            self.max_id = max_id

        self.results = []

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
        if self.win_size < 0:
            print("Invalid option for 'window_size'.\nUsing default value (30) instead.")
            self.win_size = 30

        if self.min_id < 0 or self.min_id > 100:
            print("Invalid option for 'min_identity'.\nUsing default value (0) instead.")
            self.min_id = 0

        if self.max_id < 0 or self.max_id > 100:
            print("Invalid option for 'max_identity'.\nUsing default value (100) instead.")
            self.min_id = 100

    @staticmethod
    def remove_uninformative_sites(triplet):
        """
        Remove sites that are all the same or all different
        :param triplet: n x 3 numpy array where n is the length of the alignment
        :return: a tuple containing the informative sites and the positions of informative sites
        """
        infor_sites = []
        uninfor_sites = []
        # Find positions of sites that are all the same sites or all sites that are different
        for i in range(triplet.shape[1]):
            col = triplet[:, i]
            if np.unique(col).shape[0] == 2:
                infor_sites.append(i)
            else:
                uninfor_sites.append(i)

        # Build "new alignment"
        new_trp = triplet[:, infor_sites]

        return new_trp, infor_sites, uninfor_sites

    def execute(self, alignment, triplet, G):
        """
        Performs RDP detection method for one triplet of sequences
        :param alignment: the Alignment object
        :param triplet: the triplet
        :param G: the number of possible triplets
        :return: the coordinates of the potential recombinant region and the p_value
        """
        # Get the triplet sequences
        trp_seqs = []
        for seq_num in triplet:
            trp_seqs.append(alignment[seq_num])
        trp_seqs = np.array(trp_seqs)

        # Find the informative sites
        new_trp, sites, unsites = self.remove_uninformative_sites(trp_seqs)

        # Get the three pairs of sequences
        ab = np.array([new_trp[0], new_trp[1]])
        bc = np.array([new_trp[1], new_trp[2]])
        ac = np.array([new_trp[0], new_trp[2]])

        len_trp = new_trp.shape[1]

        # 2. Sliding window over subsequence and calculate average percent identity at each position
        recombinant_regions = ''  # Recombinant regions denoted by ones
        coord = []
        for i in range(len_trp - self.win_size):
            reg_ab = ab[:, i: self.win_size + i]
            reg_bc = bc[:, i: self.win_size + i]
            reg_ac = ac[:, i: self.win_size + i]

            # Calculate percent identity in each window
            a_b, b_c, a_c = 0, 0, 0
            for j in range(reg_ab.shape[1]):
                if reg_ab[0, j] == reg_ab[1, j]:
                    a_b += 1
                if reg_bc[0, j] == reg_bc[1, j]:
                    b_c += 1
                if reg_ac[0, j] == reg_ac[1, j]:
                    a_c += 1

            percent_identity_ab = a_b / len_trp * 100
            percent_identity_bc = b_c / len_trp * 100
            percent_identity_ac = a_c / len_trp * 100

            # Identify recombinant regions
            if percent_identity_ac > percent_identity_ab or percent_identity_bc > percent_identity_ab:
                recombinant_regions += "1"
                coord.append(i)
            else:
                recombinant_regions += "0"

        # 3. Record significance of events
        recomb_idx = [(m.span()) for m in re.finditer('1+', recombinant_regions)]

        # Convert coordinates from  window-level to alignment-level and record number of windows
        coords = []
        for x, y in recomb_idx:
            coords.append((sites[x], sites[y - 1]))

        for coord in coords:
            n = coord[1] - coord[0]     # Length of putative recombinant region

            if n > 0:
                # m is the proportion of nts in common between either A or B and C in the recombinant region
                nts_in_a = trp_seqs[0][coord[0]: coord[1]]
                nts_in_c = trp_seqs[2][coord[0]: coord[1]]
                m = 0
                for i in range(n):
                    if nts_in_a[i] == nts_in_c[i]:
                        m += 1

                # p is the proportion of nts in common between either A or B and C in the entire subsequence
                id_in_seq = 0
                for j in range(trp_seqs.shape[1]):
                    if trp_seqs[0][j] == trp_seqs[2][j]:
                        id_in_seq += 1
                p = id_in_seq / trp_seqs.shape[1]

                # Calculate p_value
                val = 0
                for i in range(m, n):
                    val += (np.math.factorial(n) / (np.math.factorial(i) * np.math.factorial(n - i))) * p ** n * (
                            1 - p) ** (n - i)

                uncorr_pvalue = (len_trp / n) * val
                corr_p_value = G * uncorr_pvalue

            else:
                uncorr_pvalue = 'NS'
                corr_p_value = 'NS'

            self.results.append((triplet, coord, uncorr_pvalue, corr_p_value))

        return


class Bootscan:
    def __init__(self, win_size=200, step_size=20, use_distances=True, num_replicates=100, random_seed=3,
                 cutoff=0.7, p_val_calc='binomial', model='JC69'):
        self.win_size = win_size
        self.step_size = step_size
        self.use_distances = use_distances
        self.num_replicates = num_replicates
        self.random_seed = random_seed
        self.cutoff = cutoff
        self.p_val_calc = p_val_calc
        self.model = model

    def set_options_from_config(self, settings={}):
        pass

    def validate_options(self):
        """
        Check if the options from the config file are valid
        If the options are invalid, the default value will be used instead
        """
        pass

    def percent_diff(self, s1, s2):
        alphabet = ['A', 'T', 'G', 'C']
        diffs = 0
        num_valid = 0
        for x, y in zip(s1, s2):
            if x in alphabet and y in alphabet:
                num_valid += 1
                if x != y:
                    diffs += 1
        return diffs / num_valid if num_valid else 0

    def jc_distance(self, s1, s2):
        p_dist = self.percent_diff(s1, s2)
        return 0.75 * np.log(1 - (p_dist * 4 / 3)) if p_dist else 0

    def jc_matrix(self, aln):
        mat = []  # Linearized distances
        for s1 in aln:
            for s2 in aln:
                mat.append(self.jc_distance(s1, s2))
        return mat

    def binomial_p(self, G, l, n, m, p):
        # Calculate p_value
        val = 0
        for i in range(m, n):
            val += (np.math.factorial(n) / (np.math.factorial(i) * np.math.factorial(n - i))) * p ** n * (1 - p) ** (
                    n - i)
        p_value = G * (l / n) * val
        return p_value

    def find_potential_events(self, pair1, pair2):
        # Loop over coordinates for peaks to high bootstrap support that alternates between two pairs
        putative_regions = []
        possible_start = False
        possible_region = False
        possible_end = False
        start = 0
        end = 0

        for val in range(len(pair1)):
            # Find regions of high bootstrap support in one sequence
            if pair1[val] >= self.cutoff and pair1[val + 1] > self.cutoff:
                possible_start = True
                start = val

            if possible_start:
                if pair2[val] < self.cutoff and pair2[val + 1] < self.cutoff:
                    possible_region = True

            if possible_region:
                if pair2[val] > self.cutoff:
                    possible_end = True
                    end = val

            if possible_end:
                putative_regions.append((start, end))

        return putative_regions

    def execute(self, align, triplet, trp_pos_in_mat):
        random.seed(self.random_seed)
        # Scanning phase
        p_values = []
        all_dists = []
        for i in range(len(align), self.step_size):
            window = align[i:i + self.win_size]
            dists = {}

            # Make bootstrap replicates of alignment
            for rep in range(self.num_replicates):
                rep_window = np.random.choice(window, replace=True)  # Shuffle columns
                dist_mat = self.jc_matrix(rep_window)
                dists[align](dist_mat)
            all_dists.append(dists)

            a = align.sequence[triplet[0]]
            b = align.sequence[triplet[1]]
            c = align.sequence[triplet[2]]

            # Look at boostrap support for sequence pairs
            ab_dists = all_dists[trp_pos_in_mat[a]]
            bc_dists = all_dists[trp_pos_in_mat[b]]
            ca_dists = all_dists[trp_pos_in_mat[c]]

            # Find peaks
            ab_max = find_peaks(ab_dists)
            bc_max = find_peaks(bc_dists)
            ca_max = find_peaks(ca_dists)

            # Loop over coordinates for peaks to high bootstrap support that alternates between two pairs
            pairs = [(ab_max, bc_max), (bc_max, ca_max), (ab_max, ca_max)]
            possible_regions = []
            for pair1, pair2 in pairs:
                possible_regions.append(self.find_potential_events(pair1, pair2))

            # Find p-value for regions
            for event in possible_regions:
                n = align[event[0]] - align[event[1]]
                G = align[event[1]]  # num windows involved
                l = event[1] - event[0]

                # m is the proportion of nts in common between either A or B and C in the recombinant region
                m = np.all(align[event[0]][:, l], align[event[1]][:, l]).sum()

                # p is the proportion of nts in common between either A or B and C in the entire sequence
                p = np.all(align[event[0]], align[event[1]]).sum() / align.length

                p_val = self.binomial_p(G, l, n, m, p)

                p_values.append((p_val, event))

        return p_values
