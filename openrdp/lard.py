import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2 as chisq
from scipy.stats import chi2_contingency
from .common import identify_recombinant

class Lard:
    def __init__(self, align, step_size=20, min_win_size=100, settings=None, ref_align=None, verbose=False):
        """
        Constructs a Lard Object
        :param step_size: the step size made in exhaustive search to determine recombination point
        :param min_win_size: size of the region that undergoes the region for recombination (THIS IS SUPER IMPORTANT THIS DETERMINES THE SMALLEST VALUE THE BREAKPOINT WILL BE, see exhaust_region())
        """
        self.align = align
        self.raw_results = []
        self.step_size = step_size
        self.min_win_size = min_win_size
        self.nt_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        self.name = 'lard'

    def optimize_tree_null(self, seq_lenghts):
        """
        given a tree, find the optimized tree lengths and return the -loglikelihood to be use as null

        param lengths: list of lengths

        return: result.x list, internal node tree lengths
        return: result.fun float, -loglikelihood
        """
        result = minimize(self.log_likelihood, seq_lenghts, method='L-BFGS-B', bounds=[(1e-5, 100)] * 3)

        return result.x, result.fun


    def optimize_tree(self, seq_lenghts, seqs):
        """
        given a tree, find the optimized tree lengths. function to be used in scipy.optimize

        param lengths: list of lengths

        return: result.x list, internal node tree lengths
        """
        result = minimize(lambda x: self.log_likelihood(x,seqs), seq_lenghts, method='L-BFGS-B', bounds=[(1e-5, 100)] * 3)

        return result.x, result.fun

    def jc69_prob(self, i, j, t):
        """
        jc probability of two sequences
        """
        if i == j:
            return 0.25 + 0.75 * np.exp(-4/3 * t)
        else:
            return 0.25 - 0.25 * np.exp(-4/3 * t)
        
    def site_likelihood(self, site, t1, t2, t3):
        """
        JC log prob at one site given three nucleotides
        """
        if not all(base in self.nt_to_idx for base in site):
            return 0
        
        likelihood = 0.0
        for x in range(4):  # internal node state
            p1 = self.jc69_prob(x, self.nt_to_idx[site[0]], t1)
            p2 = self.jc69_prob(x, self.nt_to_idx[site[1]], t2)
            p3 = self.jc69_prob(x, self.nt_to_idx[site[2]], t3)
            likelihood += 0.25 * p1 * p2 * p3  # prior for x is 0.25 under JC69
        return np.log(likelihood)

    def log_likelihood(self, cur_tip_lens, seqs=None):
        """
        log prob for JC across all sites, used in scipy to minimize

        param cur_tip_lens: list of 3 ints, tip lengths from tip of tip to root
        param len_seq: int, NT length of the alignment for a given sequnece 

        """
        if seqs is None:
            seqs = self.align # need to check this

        t1, t2, t3 = cur_tip_lens
        if any(t <= 0.00001 for t in (t1, t2, t3)):
            return 1e10  # penalize small or negative branch lengths

        total = 0.0
        for i in range(len(seqs)):
            site = seqs[:, i]
            try:
                total += self.site_likelihood(site, t1, t2, t3)
            except:
                return 1e10

        return -total
    
    def exhaustive_point(self):
        """
        perform an exhaustive search for the location of breakpoints assuming there is only 1 breakpoint

        return likelihoods, list, value for likelihood of recombination
        reutn best, int, index for lowest log likelihood value
        """
        seq1, seq2, seq3 = self.align # need to check this
        likelihoods = [1e10 for i in range(len(seq1))] # null array
        single_breakpoints = []

        # check over all possible recombination points with the given step size, starting at stepsize
        for breakpoint in range(self.step_size, len(seq1), self.step_size):
            seq1_left, seq1_right = seq1[:breakpoint], seq1[breakpoint:]
            seq2_left, seq2_right = seq2[:breakpoint], seq2[breakpoint:]
            seq3_left, seq3_right = seq3[:breakpoint], seq3[breakpoint:]
            left_tree, right_tree = [1,1,1],[1,1,1] # initialization values [dist node 1, dist2, dist3]
            
            # find likelihood
            _, left_likelihood = self.optimize_tree(left_tree, np.array([seq1_left, seq2_left, seq3_left]))
            _, right_likelihood = self.optimize_tree(right_tree, np.array([seq1_right, seq2_right, seq3_right]))
            likelihoods[breakpoint]= left_likelihood+right_likelihood

            # save likelihood
            single_breakpoints.append((breakpoint, likelihoods[breakpoint]))

        return likelihoods, single_breakpoints

    def exhaustive_region(self, left_likihoods):
        """
        perform exhaustive search for the location of 2 breakpoints.
        this alrogithm works by having a region, defined by left and right index limited by the minimum amount self.min_win_size.
        
        1. the left index is incremented by self.step_size, and the right increases by the same amount.
        2. the right index is continually incremented by self.step_size until it reaches the end of the alignment
        3. once the right index reaches the end, the left increment is increase by self.step_size and the right index is put at the end of the end of the min_window_size 

        This is obviously incredibly slow due to needing to use scipy.optimize n^2 times. An alternative to this is to use a greedy algorithm that constacts the left and expands the right
        """
        seq1, seq2, seq3 = self.align
        likelihoods = []

        # check over all possible recombination points with the given step size, starting at stepsize
        for breakpoint in range(self.step_size, len(seq1) - self.step_size - self.min_win_size, self.step_size):
            left_likelihood = left_likihoods[breakpoint]

            for right_breakpoint in range(breakpoint + self.min_win_size, min(len(seq1) - self.step_size, len(seq1) - self.min_win_size), self.step_size):
                seq1_right = seq1[right_breakpoint:]
                seq2_right = seq2[right_breakpoint:]
                seq3_right = seq3[right_breakpoint:]
                right_tree = [10,10,10] # initialization values [dist node 1, dist2, dist3]
                
                # find likelihood
                _, right_likelihood = self.optimize_tree(right_tree, np.array([seq1_right, seq2_right, seq3_right]))
                curr_likelihood = left_likelihood+right_likelihood
                likelihoods.append((breakpoint, right_breakpoint, curr_likelihood))

        return likelihoods

    def remove_insignificant(self, breakpoints, null):
        """
        takes the input of exhaustive searches, filters out insig ones

        breakpoints, list tuple, breakpoint output from exhuast point/region        
        """
        valid = []
        num_trials = len(breakpoints)

        for breakpoint in breakpoints:
            ratio = 2 * (breakpoint[-1] - null)
            pval = 1 - chisq.cdf(ratio, 3) # DF of 3
            
            if pval < 0.05/num_trials:
                valid.append((pval, *breakpoint)) 
        
        return valid

    def merge_contiguous_regions(self, tuples):
        """
        Merge tuples with incrementally increasing start values (index 1) that are part of a contiguous region.
        
        tuples: list of 4-element tuples (pval, start, end, stat)
        step: expected increment between consecutive starts to consider as contiguous
        
        List of merged tuples.
        """
        if not tuples:
            return []

        merged = []
        current = list(tuples[0])

        for next_t in tuples[1:]:
            # Check if next region is contiguous with current
            if next_t[1] == current[1] + self.step_size:
                current[0] = min(current[0], next_t[0])  # min p-value
                current[2] = next_t[2]  # update end
                current[3] = next_t[3]  # use latest stat (or average if preferred)
                current[1] = current[1]  # keep original start
            else:
                merged.append(tuple(current))
                current = list(next_t)

        merged.append(tuple(current))
        return merged


    def execute(self, triplet):
        """
        run the thing
        """
        # initialize tree
        dists = [1,1,1]

        # find the optimized branch lengths for the null tree and null likelihood
        dists, null = self.optimize_tree_null(dists)
        # assuming only 1 breakpoint location
        left_likelihood, single_breakpoints = self.exhaustive_point()

        # assuming a region for breakpoints (2 breakpoints)
        region_breakpoints = self.exhaustive_region(left_likelihood)

        # remove the ones with high p-value
        single_breakpoints = self.remove_insignificant(single_breakpoints, null)
        region_breakpoints = self.remove_insignificant(region_breakpoints, null)

        # merge same region so it makes identify_recombinant computationally lighter
        single_breakpoints = self.merge_contiguous_regions(single_breakpoints)
        region_breakpoints = self.merge_contiguous_regions(region_breakpoints)

        # format into raw_results acceptable tuples while adding parents
        # (rec_name, parents, *aln_pos, final_p)
        for breakpoint in region_breakpoints:
            recs, parents = identify_recombinant(triplet, breakpoint[1:3])
            self.raw_results.append((recs, parents, breakpoint[1], breakpoint[2], breakpoint[0]))

        for breakpoint in single_breakpoints:
            recs, parents = identify_recombinant(triplet, breakpoint[1])
            self.raw_results.append((recs, parents, breakpoint[1], breakpoint[1], breakpoint[0]))

