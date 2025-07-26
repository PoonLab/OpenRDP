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
    def __init__(self, align):
        """
        Constructs a MaxChi Object
        :param min_regin_step: the step size made in exhaustive search to determine recombination point
        :param 
        """
        self.align = align
        self.raw_results = []
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


    def optimize_null_tree(self, root):
        """
        given the null tree, find the optimized tree lengths

        param root: root node from generate_null_tree

        return root: node object of internal node now with tree lengths
        """
        seq_lenghts = [root.left.dist, root.middle.dist, root.right.dist]
        result = minimize(self.log_likelihood, seq_lenghts, method='L-BFGS-B', bounds=[(1e-5, 100)] * 3)

        # update node length values
        root.left.dist, root.middle.dist, root.right.dist = result.x
        return root
    
    def exhaustive_region(self, root):
        """
        perform exhaustive search for the lcoation of 2 breakpoints 
        """

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

    def exhaustive_point(self, root):
        """
        perform an exhaustive search for the location of breakpoints assuming there is only 1 breakpoint
        """
        seq1, seq2, seq3 = self.align

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
        root = self.optimize_null_tree(root)

        # assuming only 1 breakpoint location
        single_breakpoints = self.exhaustive_point()

        # assuming a region for breakpoints (2 breakpoints)
        region_breakpoints = self.exhaustive_region()

        # remove the ones with high p-value
        single_breakpoints = self.remove_insignificant()
        region_breakpoints = self.remove_insignificant()
    
        self.raw_results.append(single_breakpoints)
        self.raw_results.append(region_breakpoints)