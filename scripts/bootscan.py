import copy
import multiprocessing
import random
from itertools import combinations

import numpy as np
from scipy.spatial.distance import pdist, squareform

from scripts.common import jc_distance


class Bootscan:
    def __init__(self, alignment, win_size=200, step_size=20, use_distances=True, num_replicates=100,
                 random_seed=3, cutoff=0.7, model='JC69', quiet=False, max_pvalue=0.05, settings=None,):
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
            self.max_pvalue = max_pvalue

        self.align = alignment
        random.seed(self.random_seed)
        np.random.seed(self.random_seed)
        if not quiet:
            print('Starting Scanning Phase of Bootscan/Recscan')
        self.dists = self.do_scanning_phase(alignment)
        if not quiet:
            print('Finished Scanning Phase of Bootscan/Recscan')

        self.raw_results = []
        self.results = []

    def set_options_from_config(self, settings):
        """
        Set the parameters of Bootscan from the  config file
        :param settings: a dictionary of settings
        """
        self.win_size = int(settings['win_size'])
        self.step_size = int(settings['step_size'])
        self.max_pvalue = abs(float(settings['max_pvalue']))
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
        np.random.seed(self.random_seed)
        random.seed(self.random_seed)
        for rep in range(self.num_replicates):
            # Shuffle columns with replacement
            rep_window = window[:, np.random.randint(0, window.shape[1], window.shape[1])]
            dist_mat = squareform(pdist(rep_window, jc_distance))
            dists.append(dist_mat)
        return dists

    def do_scanning_phase(self, align):
        """
        Perform scanning phase of the Bootscan/Recscan algorithm
        :param align: a n x m array of aligned sequences
        """
        with multiprocessing.Pool() as p:
            all_dists = p.map(self.scan, range(0, align.shape[1], self.step_size))

        return all_dists

    def execute(self, triplet):
        """
        Executes the exploratory version of the BOOTSCAN from RDP5 using the RECSCAN algorithm.
        :param triplet: a triplet object
        """

        # Look at boostrap support for sequence pairs
        ab_support = []
        bc_support = []
        ac_support = []
        for dists in self.dists:
            supports = []
            for dist_mat in dists:
                # Access pairwise distances for each pair
                ab_dist = dist_mat[triplet.idxs[0], triplet.idxs[1]]
                bc_dist = dist_mat[triplet.idxs[1], triplet.idxs[2]]
                ac_dist = dist_mat[triplet.idxs[0], triplet.idxs[2]]
                supports.append(np.argmin([ab_dist, bc_dist, ac_dist]))

            ab_support.append(np.sum(np.equal(supports, 0)) / self.num_replicates)
            bc_support.append(np.sum(np.equal(supports, 1)) / self.num_replicates)
            ac_support.append(np.sum(np.equal(supports, 2)) / self.num_replicates)

        supports = np.array([ab_support, bc_support, ac_support])
        supports_max = np.argmax(supports, axis=0)
        supports_thresh = supports > self.cutoff

        transition_window_locations = []
        for i in range(1, supports.shape[1] - self.step_size):
            if np.any(supports_thresh[:, i]):
                max1 = np.argmax(supports_thresh[:, i])
                if not supports_thresh[max1, i+1]:
                    for j in range(1, self.step_size):
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
            recomb_region_cand = triplet.sequences[recomb_candidate, event[0]: event[1]]
            other_seqs = triplet.sequences[trps[:recomb_candidate] + trps[recomb_candidate+1:], event[0]: event[1]]
            m = np.sum(np.any(recomb_region_cand == other_seqs, axis=0))

            # p is the proportion of nts in common between either A or B and C in the entire sequence
            recomb_region_cand = triplet.sequences[recomb_candidate, :]
            other_seqs = triplet.sequences[trps[:recomb_candidate] + trps[recomb_candidate + 1:], :]
            p = np.sum(np.any(recomb_region_cand == other_seqs, axis=0)) / l

            if n > 0:
                # Calculate p_value
                val = 0
                log_n_fact = np.sum(np.log(np.arange(1, n + 1)))  # Convert to log space to prevent integer overflow
                for i in range(m, n):
                    log_i_fact = np.sum(np.log(np.arange(1, i + 1)))
                    log_ni_fact = np.sum(np.log(np.arange(1, n - i + 1)))
                    if (log_i_fact + log_ni_fact) != 0:
                        val += np.math.exp(
                            (log_n_fact - (log_i_fact + log_ni_fact)) + np.log(p ** n) + np.log((1 - p) ** (n - i)))

                # Get potential recombinant and the parents
                trp_names = copy.copy(triplet.names)
                for i, name in enumerate(trp_names):
                    if i == recomb_candidate:
                        rec_name = trp_names.pop(i)
                        parents = trp_names
                        if val != 0.0:
                            self.raw_results.append((rec_name, parents, *event, val))

        return

    def merge_breakpoints(self):
        """
        Merge overlapping breakpoint locations
        :return: list of breakpoint locations where overlapping intervals are merged
        """
        self.raw_results = sorted(self.raw_results)
        results_dict = {}
        results = []

        # Gather all regions with the same recombinant
        for i, bp in enumerate(self.raw_results):
            rec_name = self.raw_results[i][0]
            parents = tuple(sorted(self.raw_results[i][1]))
            key = (rec_name, parents)
            if key not in results_dict:
                results_dict[key] = []
            results_dict[key].append(self.raw_results[i][2:])

        # Merge any locations that overlap - eg [1, 5] and [3, 7] would become [1, 7]
        for key in results_dict:
            merged_regions = []
            for region in results_dict[key]:
                region = list(region)
                old_regions = list(results_dict[key])
                for region2 in old_regions:
                    start = region[0]
                    end = region[1]
                    start2 = region2[0]
                    end2 = region2[1]
                    if start <= start2 <= end or start <= end2 <= end:
                        region[0] = min(start,start2)
                        region[1] = max(end, end2)
                        results_dict[key].remove(region2)
                merged_regions.append(region)

            # Output the results
            for region in merged_regions:
                rec_name = key[0]
                parents = key[1]
                start = region[0]
                end = region[1]
                p_value = region[2]
                if p_value < self.max_pvalue:
                    results.append((rec_name, parents, start, end, p_value))
        return results
