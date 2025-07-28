"""
changes made from original LARD scripts

exhaustive point/region search for original breakpoint:
- always use gamma with 4 rate categories and JC69 subst. model
"""


import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import minimize
from scipy.stats import chi2 as chisq
from scipy.stats import chi2_contingency
from common import Node

class Lard:
    def __init__(self, align, step_size=20, min_win_size=100):
        """
        Constructs a MaxChi Object
        :param step_size: the step size made in exhaustive search to determine recombination point
        :param min_win_size: size of the region that undergoes the region for recombination (THIS IS SUPER IMPORTANT THIS DETERMINES THE SMALLEST VALUE THE BREAKPOINT WILL BE, see exhaust_region())
        """
        self.align = align
        self.raw_results = []
        self.step_size = step_size
        self.min_win_size = min_win_size
        self.name = 'lard'

    
    def generate_null_tree(self, triplet):
        """
        generate tree for null hypothesis assuming that there is no recombination
        should 1 internal node with three terminal nodes connected to root node

        param tiplet: triplet object from commmon

        return root: Node object of the one internal node
        """
        seq1, seq2, seq3 = Node(name=triplet.names[0], dist=1), Node(name=triplet.names[1], dist=1), Node(triplet.names[2], dist=1)
        root = Node(name='root', left=seq1, middle=seq2, right=seq3, terminal=False)
        return root


    def optimize_tree(self, seq_lenghts, seqs):
        """
        given a tree, find the optimized tree lengths

        param lengths: list of lengths

        return: list, internal node tree lengths
        """
        result = minimize(self.log_likelihood, (seq_lenghts, seqs), method='L-BFGS-B', bounds=[(1e-5, 100)] * 3)

        # update node length values
        return result.x

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
        likelihood = 0.0
        for x in range(4):  # internal node state
            p1 = self.jc69_prob(x, site[0], t1)
            p2 = self.jc69_prob(x, site[1], t2)
            p3 = self.jc69_prob(x, site[2], t3)
            likelihood += 0.25 * p1 * p2 * p3  # prior for x is 0.25 under JC69
        return np.log(likelihood)

    def log_likelihood(self, cur_tip_lens, len_seq, seqs):
        """
        log prob for JC across all sites, used in scipy to minimize

        param cur_tip_lens: list of 3 ints, tip lengths from tip of tip to root
        param len_seq: int, NT length of the alignment for a given sequnece 

        """
        t1, t2, t3 = cur_tip_lens
        if any(t <= 0.00001 for t in (t1, t2, t3)):
            return 1e10  # penalize small or negative branch lengths
        total = 0.0
        for i in range(len_seq):
            site = seqs[:, i]
            try:
                total += self.site_likelihood(site, t1, t2, t3)
            except:
                return 1e10
        return -total  # minimize negative log-likelihood

    def exhaustive_point(self):
        """
        perform an exhaustive search for the location of breakpoints assuming there is only 1 breakpoint

        return likelihoods, list, value for likelihood of recombination
        reutn best, int, index for lowest log likelihood value
        """
        seq1, seq2, seq3 = self.align
        likelihoods = [1e10 for i in range(len(seq1))] # null array
        single_breakpoints = []

        # check over all possible recombination points with the given step size, starting at stepsize
        for breakpoint in range(self.step_size, len(seq1), self.step_size):
            seq1_left, seq1_right = seq1[:breakpoint], seq1[breakpoint:]
            seq2_left, seq2_right = seq2[:breakpoint], seq2[breakpoint:]
            seq3_left, seq3_right = seq3[:breakpoint], seq3[breakpoint:]
            left_tree, right_tree = [10,10,10],[10,10,10] # initialization values [dist node 1, dist2, dist3]
            
            # find likelihood
            left_likelihood = self.optimize_tree(left_tree, np.array([seq1_left, seq2_left, seq3_left]))
            right_likelihood = self.optimize_tree(right_tree, np.array([seq1_right, seq2_right, seq3_right]))
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
        best = np.inf # track index of best

        # check over all possible recombination points with the given step size, starting at stepsize
        for breakpoint in range(self.step_size, len(seq1) - self.step_size - self.win_size, self.step_size):
            left_likelihood = left_likihoods[breakpoint]

            for right_breakpoint in range(breakpoint + self.min_win_size, len(seq1) - self.step_size - self.win_size, self.step_size):
                seq1_right = seq1[right_breakpoint:]
                seq2_right = seq2[right_breakpoint:]
                seq3_right = seq3[right_breakpoint:]
                right_tree = [10,10,10],[10,10,10] # initialization values [dist node 1, dist2, dist3]
                
                # find likelihood
                right_likelihood = self.optimize_tree(right_tree, np.array([seq1_right, seq2_right, seq3_right]))
                curr_likelihood = left_likelihood+right_likelihood
                best = min(curr_likelihood, best)
                likelihoods.append((breakpoint, right_likelihood, curr_likelihood))

        return likelihoods, best

    def remove_insignificant(self, breakpoints):
        """
        takes the input of exhaustive searches
        """
        return

    def execute(self, triplet):
        """
        run the thing
        """

        # initialize tree
        root = self.generate_null_tree(triplet.sequences)

        # find the optimized branch lengths for the null tree
        root = self.optimize_tree(root, self.align)

        # assuming only 1 breakpoint location
        left_likelihood, single_breakpoints = self.exhaustive_point()

        # assuming a region for breakpoints (2 breakpoints)
        region_breakpoints = self.exhaustive_region(left_likelihood)

        # remove the ones with high p-value
        single_breakpoints = self.remove_insignificant()
        region_breakpoints = self.remove_insignificant()
    
        self.raw_results.append(single_breakpoints)
        self.raw_results.append(region_breakpoints)