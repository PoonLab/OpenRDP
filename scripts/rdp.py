from itertools import combinations


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
