import multiprocessing
import numpy as np
from scipy.spatial.distance import pdist, squareform
import functools
import random
from math import factorial
from common import generate_triplets


class Bootscan:
    def __init__(self, alignment, s_names, settings=None, win_size=200, step_size=20, use_distances=True,
                 num_replicates=100, random_seed=3, cutoff=0.7, model='JC69'):
        if settings:
            self.set_options_from_config(settings)
            self.validate_options(alignment)
        else:
            self.win_size = win_size
            self.step_size = step_size
            self.use_distances = use_distances
            self.num_replicates = num_replicates
            self.random_seed = random_seed
            self.cutoff = cutoff
            self.model = model

        self.names = s_names
        self.align = alignment
        random.seed(self.random_seed)
        print('Starting Scanning Phase of Bootscan/Recscan')
        self.dists = self.do_scanning_phase(alignment)
        print('Finished Scanning Phase of Bootscan/Recscan')

        self.results = {}

    def set_options_from_config(self, settings):
        """
        Set the parameters of Siscan from the  config file
        :param settings: a dictionary of settings
        """
        self.win_size = int(settings['win_size'])
        self.step_size = int(settings['step_size'])
        self.num_replicates = int(settings['num_replicates'])
        self.random_seed = int(settings['random_seed'])
        self.cutoff = float(settings['cutoff_percentage'])

        if settings['scan'] == 'distances':
            self.use_distances = True

    def validate_options(self, alignment):
        """
        Check if the options from the config file are valid
        If the options are invalid, the default value will be used instead
        """
        if self.win_size < 0 or self.win_size > alignment.shape[1]:
            print("Invalid option for 'window_size'.\nUsing default value (200) instead.")
            self.win_size = 200

        if self.step_size < 0 or self.step_size >= self.win_size:
            print("Invalid option for 'step_size'.\nUsing default value (20) instead.")
            self.step_size = 20

        if self.num_replicates < 0:
            print("Invalid option for 'num_replicates'.\nUsing default value (100) instead.")
            self.num_replicates = 100

        if self.random_seed < 0:
            print("Invalid option for 'random_seed'.\nUsing default value (3) instead.")
            self.random_seed = 3

        if self.cutoff <= 0 or self.cutoff > 1:
            print("Invalid option for 'cutoff_percentage'.\nUsing default value (0.7) instead.")
            self.cutoff = 0.7

    def percent_diff(self, s1, s2):
        s1_valid = (s1 == 'A') | (s1 == 'T') | (s1 == 'G') | (s1 == 'C')
        s2_valid = (s2 == 'A') | (s2 == 'T') | (s2 == 'G') | (s2 == 'C')
        valid = s1_valid & s2_valid
        diffs = np.sum(s1[valid] != s2[valid])
        num_valid = valid.sum()
        return diffs.sum() / num_valid if num_valid else 0

    def jc_distance(self, s1, s2):
        """
        Calculate the pairwise Jukes-Cantor distance between 2 sequences
        :param s1: the first sequence
        :param s2: the second sequence
        :return: the pairwise JC69 distance between 2 sequences
        """
        p_dist = self.percent_diff(s1, s2)
        if p_dist >= 0.75:
            return 1
        return -0.75 * np.log(1 - (p_dist * 4 / 3)) if p_dist else 0

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

    def scan(self, i):
        window = self.align[:, i:i + self.win_size]
        # Make bootstrap replicates of alignment
        dists = []
        for rep in range(self.num_replicates):
            # Shuffle columns with replacement
            rep_window = window[:, np.random.randint(0, window.shape[1], window.shape[1])]
            dist_mat = squareform(pdist(rep_window, self.jc_distance))
            dists.append(dist_mat)
        return dists

    def do_scanning_phase(self, align):
        """
        Perform scanning phase of the Bootscan/Recscan algorithm
        :param align: a n x m array of aligned sequences
        """
        scan = functools.partial(self.scan, align)
        with multiprocessing.Pool() as p:
            all_dists = p.map(scan, range(0, align.shape[1], self.step_size))

        return all_dists

    def execute(self):
        """
        Executes the exploratory version of the BOOTSCAN from RDP5 uses the RECSCAN algorithm,
            which does not require that recombinants are known
        """
        G = int(factorial(self.align.shape[0]) / (factorial(3) * factorial((self.align.shape[0]) - 3)))
        trp_count = 1
        num_seqs = self.align.shape[0]
        for trp_idx in generate_triplets(self.align):
            print("Scanning triplet {} / {}".format(trp_count, G))
            trp_count += 1

            # Get the triplet sequences
            trp_seqs = []
            for seq_num in trp_idx:
                trp_seqs.append(self.align[seq_num])
            trp_seqs = np.array(trp_seqs)

            # Detection phase
            pairs = ((0, 1), (1, 2), (0, 2))

            # Look at boostrap support for sequence pairs
            ab_support = []
            bc_support = []
            ac_support = []
            for dists in self.dists:
                supports = []
                for dist_mat in dists:
                    # Access pairwise distances for each pair
                    ab_dist = dist_mat[trp_idx[0], trp_idx[1]]
                    bc_dist = dist_mat[trp_idx[1], trp_idx[2]]
                    ac_dist = dist_mat[trp_idx[0], trp_idx[2]]
                    supports.append(np.argmin([ab_dist, bc_dist, ac_dist]))

                ab_support.append(np.sum(np.equal(supports, 0)) / self.num_replicates)
                bc_support.append(np.sum(np.equal(supports, 1)) / self.num_replicates)
                ac_support.append(np.sum(np.equal(supports, 2)) / self.num_replicates)

            supports = np.array([ab_support, bc_support, ac_support])
            supports_max = np.argmax(supports, axis=0)
            supports_thresh = supports > self.cutoff

            transition_window_locations = []
            for i in range(1, supports.shape[1] - self.step_size):
                if np.any(supports_thresh[:,i]):
                    max1 = np.argmax(supports_thresh[:,i])
                    if not supports_thresh[max1,i+1]:
                        for j in range(1,self.step_size):
                            max2 = supports_max[i+j]
                            if max2 != max1 and supports_thresh[max2, i+j]:
                                transition_window_locations.append(i)
                                transition_window_locations.append(i+j)
                                break

            # Identify areas where the bootstrap support alternates between two different sequence pairs
            transition_window_locations = [0] + transition_window_locations + [supports.shape[1] - 1]
            possible_regions = []
            groupings = ((0, 2), (0, 1), (1, 2))
            trps = (0, 1, 2)
            for rec_pot in range(3):
                for i in range(len(transition_window_locations) - 1):
                    begin = transition_window_locations[i]
                    end = transition_window_locations[i + 1]
                    pair = supports_max[begin]
                    if np.all(supports_thresh[pair, begin:end+1]) and pair in groupings[rec_pot]:
                        region = (rec_pot, (begin * self.step_size + self.win_size // 2, end * self.step_size + self.win_size // 2))
                        possible_regions.append(region)

            # Find p-value for regions
            for recomb_candidate, event in possible_regions:
                n = event[1] - event[0]
                l = self.align.shape[1]

                # m is the proportion of nts in common between either A or B and C in the recombinant region
                recomb_region_cand = trp_seqs[recomb_candidate, event[0]: event[1]]
                other_seqs = trp_seqs[trps[:recomb_candidate] + trps[recomb_candidate+1:], event[0]: event[1]]
                m = np.sum(np.any(recomb_region_cand == other_seqs, axis=0))

                # p is the proportion of nts in common between either A or B and C in the entire sequence
                recomb_region_cand = trp_seqs[recomb_candidate, :]
                other_seqs = trp_seqs[trps[:recomb_candidate] + trps[recomb_candidate + 1:], :]
                p = np.sum(np.any(recomb_region_cand == other_seqs, axis=0)) / l

                if n > 0:
                    # Calculate p_value
                    val = 0
                    log_n_fact = np.sum(np.log(np.arange(1, n + 1)))  # Convert to log space to prevent integer overflow
                    for i in range(m, n):
                        log_i_fact = np.sum(np.log(np.arange(1, i + 1)))
                        log_ni_fact = np.sum(np.log(np.arange(1, n - i + 1)))
                        val += np.math.exp(
                            (log_n_fact - (log_i_fact + log_ni_fact)) + np.log(p ** n) + np.log((1 - p) ** (n - i)))

                    uncorr_pvalue = (l / n) * val
                    corr_p_value = G * uncorr_pvalue

                    rec_name = self.names[trp_idx[recomb_candidate]]

                    try:
                        self.results[names].append((coord, uncorr_pvalue, corr_p_value))
                    except KeyError:
                        self.results[names] = (coord, uncorr_pvalue, corr_p_value)

                    self.results.append((rec_name, *event, uncorr_pvalue, corr_p_value))

                else:
                    return

            return
