from itertools import combinations, product

import numpy as np
from scipy.stats import chi2_contingency
from scipy.stats import pearsonr
import copy
import sys

class node:
    """
    node class for UPGMA tree

    dist is distance from parent 
    total_dist is distance from current node to the end of the tree
    """
    def __init__(self, name=None, left=None, right=None, dist=0, p_dist=0, terminal=True):
        self.left = left
        self.right = right
        self.name = name
        self.dist = dist
        self.p_dist = p_dist # distance from parent node to this node
        self.terminal = terminal # is terminal branch 

def merge_breakpoints(raw_results, max_pvalue=100):
    """
    took from siscan
    added max_pvalue clause
    """
    results_dict = {}
    results = []

    raw_results = sorted(raw_results)

    # Gather all regions with the same recombinant
    for i, bp in enumerate(raw_results):
        rec_name = raw_results[i][0]
        parents = tuple(sorted(raw_results[i][1]))
        key = (rec_name, parents)
        if key not in results_dict:
            results_dict[key] = []
        results_dict[key].append(raw_results[i][2:])

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
                    region[0] = min(start, start2)
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
            if float(p_value) < max_pvalue:
                    results.append((rec_name, parents, start, end, p_value))

    return results

def read_fasta(handle):
    """
    Converts a FASTA formatted file to a tuple containing a list of headers and sequences
    :param handle: file stream for the FASTA file
    :return: tuple of headers (list) and sequences (list)
    """
    # Verifies file have the correct formatting
    if hasattr(handle, "seek"):
        found = False
        for line in handle:
            if line.startswith('>'):
                found = True
                break
        if not found:
            print(f"Error: Input {handle.name} does not appear to be in a FASTA format.")
            sys.exit(1)
        handle.seek(0)  # Reset pointer to beginning of file

    headers, seqs = [], []
    sequence, h = '', ''
    for line in handle:
        if line.startswith('>'):
            if len(sequence) > 0:
                headers.append(h)
                seqs.append(sequence)
                sequence = ''
            h = line.strip('>\t\n\r')
        else:
            sequence += line.strip('\n\r').upper()

    # Handle the last entry
    seqs.append(sequence)
    headers.append(h)
    return headers, seqs


def calculate_chi2(c_table, max_pvalue):
    """
    Computes the chi-squared value and returns the chi-squared and p-value if the difference is significant
    :param c_table: a 2x2 contingency table
    :param max_pvalue: the p-value threshold
    :return: a tuple containing the chi-squared value and p-value
    """
    # Compute chi-squared value if the expected frequencies are valid
    if (c_table[0][0] > 0 and c_table[0][1] > 0) or (c_table[1][0] > 0 and c_table[1][1] > 0):

        chi2, p_value, _, _ = chi2_contingency(c_table)

        # Record only significant events
        if p_value < max_pvalue:
            return chi2, p_value

    return None, None


def reduce_to_unique_seqs(aln):
    """
    Remove identical sequences
    :param aln: list of aligned sequences
    :return: a list of unique sequences
    """
    return list(set(aln))


def percent_diff(s1, s2):
    s1_valid = (s1 == 'A') | (s1 == 'T') | (s1 == 'G') | (s1 == 'C')
    s2_valid = (s2 == 'A') | (s2 == 'T') | (s2 == 'G') | (s2 == 'C')
    valid = s1_valid & s2_valid
    diffs = np.sum(s1[valid] != s2[valid])
    num_valid = valid.sum()
    return diffs.sum() / num_valid if num_valid else 0


def jc_distance(s1, s2):
    """
    Calculate the pairwise Jukes-Cantor distance between 2 sequences
    :param s1: the first sequence
    :param s2: the second sequence
    :return: the pairwise JC69 distance between 2 sequences
    """
    p_dist = percent_diff(s1, s2)
    if p_dist >= 0.75:
        return 1
    return -0.75 * np.log(1 - (p_dist * 4 / 3)) if p_dist else 0


def find_parent(id1, id2, id3, root):
    """
    determine which of the two are the closest related

    id1, id2, id3: names of the sequences in the triplet
    root: root node of dendrogram

    return: two seq ids
    """
    d1 = find_dist(root, id1, id2)
    d2 = find_dist(root, id1, id3)
    d3 = find_dist(root, id2, id3)

    res = [(d1, id1, id2), (d2, id1, id3), (d3, id2, id3)]

    # return the two ids that are the smallest, these are the major parents
    # should probably clean this up
    return res[list(map(lambda x:x[0], res)).index(min(map(lambda x:x[0], res)))][1:]


def setup_upgma(seqs, names):
    """
    calculate the pairwise distance matrix for upgma and get list of tree nodes
    
    :param seqs: list of nucleotide sequence from import_data
    :param names: list of sequence names from import_data

    :return tree: list, nodes
    :return matrix: pairwise distance matrix where [y][x] is the two sequences indexed by tree
    """

    # get list of of ndoes, all are terminal with names of seuqences
    tree = [node(name=name) for name in names]
    matrix = [[0 for j in names] for i in names] # emtpy distance matrix with size of len(seq)

    for y, seq in enumerate(seqs):
        seq = np.array(seq)
        for x, seq2 in enumerate(seqs[y:]):
            x += y
            seq2 = np.array(seq2)

            matrix[y][x] = sum(seq!=seq2)


    matrix = np.array([np.array(row) for row in matrix])
    matrix += matrix.T
    return tree, matrix


def recalculate_dist(matrix, pos1, pos2):
    """
    Recalculate the distance matrix after merging two clusters.

    Parameters:
    matrix (np.ndarray): Pairwise distance matrix.
    pos1 (int): Index of the first cluster (to keep).
    pos2 (int): Index of the second cluster (to merge and remove).

    Returns:
    np.ndarray: Updated distance matrix.
    """
    # Create a copy of the distance matrix
    new_dist_mat = matrix.copy()

    # Update distances for the merged cluster
    for i in range(len(matrix)):
        if i != pos1 and i != pos2:  # Skip the merged clusters
            dist1 = matrix[pos1, i]
            dist2 = matrix[pos2, i]
            new_dist_mat[pos1, i] = (dist1 + dist2) / 2
            new_dist_mat[i, pos1] = (dist1 + dist2) / 2

    # Remove the row and column corresponding to pos2
    new_dist_mat = np.delete(new_dist_mat, pos2, axis=0)
    new_dist_mat = np.delete(new_dist_mat, pos2, axis=1)

    return new_dist_mat

def upgma(headers, matrix):
    """
    headers, list, ids/nodes
    matrix, numpy 2d array

    returns: root node, 
    """
    # get closest
    (x, y), smallest = find_min(matrix)
    id1, id2 = headers[x], headers[y]
    b_len = smallest / 2 

    # turn closest into new node
    new = node(name = id1.name + ';' + id2.name, left = id1, right = id2, dist = smallest/2, terminal=False)
    id1.p_dist = b_len - id1.dist
    id2.p_dist = b_len - id2.dist
    
    # update nodes
    headers[x] = new
    headers.pop(y)
    
    # create new distance matrix
    new_dist_matrix = recalculate_dist(matrix, x, y)
    
    if len(headers) == 1: # last node
        return headers[0], matrix

    return upgma(headers, new_dist_matrix)
        
    
def find_min(matrix):
    """
    matrix, adj matrix
    return inds, (column, row indexes) of closest
    """
    smallest = np.inf
    ind = None
    for row_ind, row in enumerate(matrix):
        for column_ind, value in enumerate(row):
            if value != 0 and value < smallest:
                ind = (column_ind, row_ind)
                smallest = value
    return ind, smallest

def find_dist(curr, n1, n2):
    
    if curr.terminal: # this is the last node
        if curr.name in (n1, n2):
            return curr.p_dist, True
        
        return 0, False
    
    # recursively traverse the tree
    l1, f1 = find_dist(curr.left, n1, n2)
    l2, f2 = find_dist(curr.right, n1, n2)

    # both terminal nodes are found, don't add distance upwards now
    if f1 and f2:
        return l1 + l2, True
    
    # if found only one, don't add otherside
    if f1:
        return l1 + curr.p_dist, True
    
    # other case
    if f2:
        return l2 + curr.p_dist, True

    # if this node contains no children/itself isn't the node we care about
    return 0, False


def all_items_equal(x):
    """
    Check if all items in a list are identical
    :param x: the list
    :return: True if all the items are identical, false otherwise
    """
    return x.count(x[0]) == len(x)


def identify_recombinant(trp, aln_pos):
    """
    Find the most likely recombinant sequence using the PhPr method described in the RDP5 documentation
    and Weiler GF (1998) Phylogenetic profiles: A graphical method for detecting genetic recombinations
    in homologous sequences. Mol Biol Evol 15: 326–335
    :return: name of the recombinant sequence and the names of the parental sequences
    """
    upstream_dists = []
    downstream_dists = []

    # Get possible breakpoint locations
    for i in range(len(trp.names)):
        upstream_dists.append([])
        downstream_dists.append([])
        for j in range(len(trp.names)):
            if i != j:
                # Calculate pairwise Jukes-Cantor distances for regions upstream and downstream of breakpoint
                upstream_dists[i].append(jc_distance(trp.sequences[i][0: aln_pos[0]],
                                                     trp.sequences[j][0: aln_pos[0]]))
                downstream_dists[i].append(jc_distance(trp.sequences[i][aln_pos[1]: trp.sequences.shape[1]],
                                                       trp.sequences[j][aln_pos[1]: trp.sequences.shape[1]]))

    # Calculate Pearson's correlation coefficient for 2 lists
    r_coeff = [0, 0, 0]
    if all_items_equal(upstream_dists[0]) or all_items_equal(downstream_dists[0]):
        r_coeff[0] = float('NaN')
    elif all_items_equal(upstream_dists[1]) or all_items_equal(downstream_dists[1]):
        r_coeff[1] = float('NaN')
    elif all_items_equal(upstream_dists[2]) or all_items_equal(downstream_dists[2]):
        r_coeff[2] = float('NaN')

    else:
        r_coeff[0], _ = pearsonr(upstream_dists[0], downstream_dists[0])
        r_coeff[1], _ = pearsonr(upstream_dists[1], downstream_dists[1])
        r_coeff[2], _ = pearsonr(upstream_dists[2], downstream_dists[2])

    # Most likely recombinant sequence is sequence with lowest coefficient
    trp_names = copy.copy(trp.names)
    rec_name = trp_names.pop(np.argmin(r_coeff))
    p_names = trp_names

    return rec_name, p_names


class TripletGenerator:
    def __init__(self, alignment, seq_names, ref_align=None, ref_names=None):
        """
        :param alignment:  numpy.array, a character array of 'n' sequences
        :param seq_names:  list, sequence labels / names
        :param ref_align:  numpy.array, a character array of 'n' sequences
        :param ref_names:  list, reference labels / names
        """
        self.alignment = alignment
        self.seq_names = seq_names
        self.ref_align = ref_align
        self.ref_names = ref_names

        if isinstance(self.ref_align, np.ndarray):
            self.combinations = product(
                list(combinations(range(self.alignment.shape[0]), 1)),
                list(combinations(range(self.ref_align.shape[0]), 2))
            )  # Example: ((1,), (1, 2))
        else:
            self.combinations = combinations(
                range(self.alignment.shape[0]),  # number of "rows" (sequences)
            3)  # all combinations of three elements

    def __iter__(self):
        return self

    def __next__(self):
        return self.next()

    def next(self):
        idxs = next(self.combinations)

        if isinstance(self.ref_align, np.ndarray):
            query_idx = idxs[0]
            ref_idxs = idxs[1]
            seqs = np.take(self.alignment, query_idx, axis=0) # one row of query sequences
            refs = np.take(self.ref_align, ref_idxs, axis=0) # two rows of ref sequences
            seqs = np.concatenate((seqs, refs))
            names = [self.seq_names[i] for i in query_idx] + [self.ref_names[i] for i in ref_idxs]
        else:
            seqs = np.take(self.alignment, idxs, axis=0)
            names = [self.seq_names[i] for i in idxs]
            
        return Triplet(seqs, names, idxs)


class Triplet:
    def __init__(self, seqs, names, idxs=None):
        self.sequences = seqs
        self.names = names
        self.idxs = idxs
        self.poly_sites_align, self.poly_sites = self.remove_monomorphic_sites()
        self.info_sites_align, self.info_sites, self.uninfo_sites = self.remove_uninformative_sites()

    def get_sequence_name(self, trp_idx):
        return self.names[trp_idx]

    def get_seq_from_name(self, trp_name):
        for i, name in enumerate(self.names):
            if trp_name == name:
                return self.sequences[i]

    def get_triplets(self):
        return self.sequences

    def get_trp_names(self):
        return self.names

    def remove_uninformative_sites(self):
        """
        Remove sites that are all the same or all different
        """
        infor_sites = []
        uninfor_sites = []
        # Find positions of sites that are all the same sites or all sites that are different
        for i in range(self.sequences.shape[1]):
            col = self.sequences[:, i]
            if np.unique(col).shape[0] == 2:
                infor_sites.append(i)
            else:
                uninfor_sites.append(i)

        # Build "new alignment"
        new_aln = self.sequences[:, infor_sites]

        return new_aln, infor_sites, uninfor_sites

    def remove_monomorphic_sites(self):
        """
        Remove monomorphic sites
        :return: a tuple containing the polymorphic sites and the positions of polymorphic sites
        """
        poly_sites = []
        # Find positions of polymorphic sites
        for i in range(self.sequences.shape[1]):
            col = self.sequences[:, i]
            if not np.all(col == col[0]):
                # Record polymorphic sites because it seems shorter than
                poly_sites.append(i)

        # Build "new alignment"
        new_align = self.sequences[:, poly_sites]

        return new_align, poly_sites

    def get_win_size(self, offset, win_size, fixed_win_size, num_var_sites, frac_var_sites):
        """
        Sets the size (can be fixed or variable) of the sliding window
        :param offset: the step size
        """
        # Fixed window length
        if fixed_win_size:
            return win_size

        # User specified number of variable sites
        else:
            if len(self.poly_sites) > 1.5 * win_size:
                return round(0.75 * len(self.poly_sites))

            elif num_var_sites:
                curr = 0
                while curr < num_var_sites - 1:
                    curr += 1
                return self.poly_sites[curr] + offset

            # User specified proportion of variable sites
            else:
                num_var_sites = 1
                frac = 0
                while frac < frac_var_sites:
                    for var_site in self.poly_sites:
                        frac = num_var_sites / var_site
                        num_var_sites += 1
                return round(frac * num_var_sites)
