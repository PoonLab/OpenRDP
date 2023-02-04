import sys
import os
import configparser
import numpy as np
from itertools import combinations

from openrdp import __path__ as basepath
from openrdp.bootscan import Bootscan
from openrdp.chimaera import Chimaera
from openrdp.common import generate_triplets, Triplet
from openrdp.geneconv import GeneConv
from openrdp.maxchi import MaxChi
from openrdp.rdp import RdpMethod
from openrdp.siscan import Siscan
from openrdp.threeseq import ThreeSeq


# list of all recombination detection methods
aliases = {
    'geneconv': {'key': "Geneconv", 'method': GeneConv},
    'bootscan': {'key': 'Bootscan', 'method': Bootscan},
    'maxchi': {'key': 'MaxChi', 'method': MaxChi},
    'siscan': {'key': 'Siscan', 'method': Siscan},
    'chimaera': {'key': "Chimaera", 'method': Chimaera},
    'threeseq': {'key': "3Seq", 'method': ThreeSeq},
    'rdp': {'key': "RDP", 'method': RdpMethod}
}
DNA_ALPHABET = ['A', 'T', 'G', 'C', '-', '*']


class ScanResults:
    """ Object to return from Scanner, derived from dict """
    def __init__(self, d):
        self.dict = d

    def write(self, outfile):
        """
        Return CSV-formatted string suitable for writing to a file
        :param outfile:  file to write output
        """
        outfile.write('Method,StartLocation,EndLocation,Recombinant,Parent1,Parent2,Pvalue\n')
        for method, events in self.dict.items():
            for e in events:
                if method == 'geneconv':
                    outfile.write(f'Geneconv,{e[2][0]},{e[2][1]},{e[0]},{e[1][0]},{e[1][1]},{e[3]}\n')
                else:
                    outfile.write(f'{method.title()},{e[2]},{e[3]},{e[0]},{e[1][0]},{e[1][1]},{e[4]}\n')

    def __str__(self):
        """ Print results to console """
        print('\n{:<20} {:<20} {:<20} {:<20} {:<20} {:<20} {:<20}'
              .format('Method', 'StartLocation', 'EndLocation', 'Recombinant',
                      'Parent1', 'Parent2', 'Pvalue'))

        for method, events in self.dict.items():
            key = aliases[method]['key']
            for e in events:
                if method == 'geneconv':
                    print(f"{key:<20} {e[2][0]:<20} {e[2][1]:<20} {e[0]:<20} "
                          f"{e[1][0]:<20} {e[1][1]:<20} {e[3]:<20}")
                else:
                    print(f"{key:<20} {e[2]:<20} {e[3]:<20} {e[0]:<20} "
                          f"{e[1][0]:<20} {e[1][1]:<20} {e[4]:<20}")


class Scanner:
    def __init__(self, names, infile, outfile, cfg=None, methods=None,
                 quiet=False):
        """
        :param names:  list, sequence labels
        :param infile:  str, path to input FASTA
        :param outfile:  str, path to write output CSV
        :param cfg:  str, path to configuration file.  Defaults to None, causing
                     each method to use default settings.
        :param methods:  tuple, names of methods to run
        :param quiet:  bool, if True, suppress console messages
        """
        self.seq_names = names
        self.infile = infile

        self.config = None
        self.cfg_file = cfg
        if self.cfg_file:
            self.config = configparser.ConfigParser()
            self.config.read(self.cfg_file)

        self.methods = methods
        self.quiet = quiet
        self.outfile = outfile

    def print(self, msg):
        """ Implements self.quiet """
        if not self.quiet:
            print(msg)

    def get_config(self):
        """ Return ConfigParser as dict """
        d = {}
        for section in self.config.sections():
            d.update({section: {}})
            for key, val in self.config[section].items():
                if val.isnumeric() or val in ['True', 'False']:
                    val = eval(val)
                d[section].update({key: val})
        return d

    def set_config(self, usr):
        """
        Modify configuration by passing a dict with same structure
        :param usr:  dict, should be same structure as return value of get_config()
        """
        for section in self.config.sections():
            if section not in usr:
                continue
            for key in self.config[section]:
                if key not in usr[section]:
                    continue
                self.config[section][key] = str(usr[section][key])

    def run_scans(self, aln):
        """
        Run the selected recombination detection analyses
        :param aln:  list, sequences imported from FASTA file
        """
        aln = list(set(aln))  # Remove identical sequences

        # Check that sequences are of the same length
        seqlens = [len(s) for s in aln]
        if len(set(seqlens)) > 1:
            print(f"ERROR: Sequences in input {self.infile} are not the same length!")
            sys.exit()

        # check that all sequences are the same length
        seqlen = set([len(s) for s in aln])
        if len(seqlen) > 1:
            print("ERROR: The input sequences are different lengths!  Input must be aligned.")
            sys.exit()

        # Create an m x n array of sequences
        alignment = np.array(list(map(list, aln)))

        # prepare return value
        results = ScanResults(dict([(method, {}) for method in aliases.keys()]))

        # Run 3Seq
        if 'threeseq' in self.methods:
            three_seq = ThreeSeq(self.infile)
            self.print("Starting 3Seq Analysis")
            results.dict['threeseq'] = three_seq.execute()
            self.print("Finished 3Seq Analysis")

        # Run GENECONV
        if 'geneconv' in self.methods:
            # Parse output file if available
            if self.config:
                geneconv = GeneConv(settings=dict(self.config.items('Geneconv')))
            else:
                geneconv = GeneConv()  # default config

            self.print("Starting GENECONV Analysis")
            results.dict['geneconv'] = geneconv.execute(self.infile)
            self.print("Finished GENECONV Analysis")

        # Exit early if 3Seq and Geneconv are the only methods selected
        check = set(aliases.keys()).intersection(self.methods)
        if len(check) == 0:
            return results

        tmethods = []
        for alias, a in aliases.items():
            self.print(f"Setting up {alias} analysis...")
            if self.config:
                settings = dict(self.config.items(a['key']))
                tmethods.append(a['method'](alignment, settings=settings, quiet=self.quiet))
            else:
                tmethods.append(a['method'](alignment, quiet=self.quiet))

        # iterate over all triplets in the alignment
        trp_count = 1
        total_num_trps = sum(1 for _ in combinations(range(alignment.shape[0]), 3))
        for trp in generate_triplets(alignment):
            triplet = Triplet(alignment, self.seq_names, trp)
            self.print("Scanning triplet {} / {}".format(trp_count, total_num_trps))
            trp_count += 1
            for tmethod in tmethods:
                tmethod.execute(triplet)

        # Process results by joining breakpoint locations that overlap
        for tmethod in tmethods:
            results.dict[tmethod.name] = tmethod.merge_breakpoints()

        return results


def valid_alignment(alignment):
    """
    Check that the input alignment is valid
    :param alignment: a list of lists containing the sequence headers and the
                      aligned sequences
    :return True if the alignment is valid, false otherwise
    """
    aln_len = len(alignment[0])
    for seq in alignment:
        if len(seq) != aln_len:
            print("Improper alignment! Not all alignments are the same length.")
            return False
    return True


def valid_chars(alignment):
    """
    Check that the alignment only contains valid characters
    :param alignment: a list of lists containing the sequence headers and the
                      aligned sequences
    :return: True if the alignment contains only valid characters, False
             otherwise
    """
    for s in alignment:
        if not all(pos in DNA_ALPHABET for pos in s):
            print("Alignment contains invalid characters.")
            return False
    return True


def read_fasta(handle):
    """
    Converts a FASTA formatted file to a tuple containing a list of headers and sequences
    :param handle: file stream for the FASTA file
    :return: tuple of headers (list) and sequences (list)
    """
    headers, seqs = [], []
    sequence, h = '', ''

    # Verifies files have the correct formatting
    for i, line in enumerate(handle):
        if line.startswith('>'):
            break
        else:
            print("No header")
            raise NameError

    # Reset pointer to beginning of file
    if hasattr(handle, 'seek'):
        handle.seek(0)

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


def openrdp(infile, outfile, cfg=None, methods=None, quiet=False):
    """
    Main function
    :param infile:  str, path to input FASTA file
    :param outfile:  str, path to write output CSV
    :param methods:  list, names of methods to run
    :param quiet:  bool, if True, suppress console messages
    :return:  dict, results from each method
    """
    if methods is None:
        methods = list(aliases.keys())

    # Check that the OS is valid
    platform = sys.platform
    try:
        platform.startswith("win") or platform.startswith("linux") or sys.platform == "darwin"
    except OSError:
        print("OSError: {} is not supported".format(sys.platform))

    if not (infile.endswith('.fa') or infile.endswith('.fasta')):
        print(f"Expected '.fa' or '.fasta' suffix for input FASTA {infile}, ignoring")

    # import labels and sequences from FASTA file
    with open(infile) as in_handle:
        names, aln = read_fasta(in_handle)

    if not valid_alignment(aln) and not valid_chars(aln):
        sys.exit(1)

    scanner = Scanner(names, infile, outfile, cfg=cfg, methods=methods, quiet=quiet)
    results = scanner.run_scans(aln)
    return results
