import _io
import sys
import os
import configparser
import numpy as np
from collections import OrderedDict
from scipy.special import comb as ncomb

from openrdp import __path__ as basepath
from openrdp.bootscan import Bootscan
from openrdp.chimaera import Chimaera
from openrdp.common import TripletGenerator, Triplet, read_fasta
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
DNA_ALPHABET = ['A', 'T', 'G', 'C', '-', 'N']


class ScanResults:
    """ Object to return from Scanner, derived from dict """
    def __init__(self, d):
        self.dict = d

    def write(self, outfile):
        """
        Return CSV-formatted string suitable for writing to a file
        :param outfile:  file to write output
        """
        outfile.write('Method,Start,End,Recombinant,Parent1,Parent2,Pvalue\n')
        for method, events in self.dict.items():
            for e in events:
                if method == 'geneconv':
                    outfile.write(f'Geneconv,{e[2][0]},{e[2][1]},{e[0]},{e[1][0]},{e[1][1]},{e[3]}\n')
                else:
                    outfile.write(f'{method.title()},{e[2]},{e[3]},{e[0]},{e[1][0]},{e[1][1]},{e[4]}\n')

    def __str__(self):
        """ Print results to console """
        outstr = '\n' + '\t'.join([
            'Method  ', 'Start', 'End', 'Recombinant', 'Parent1', 'Parent2',
            'Pvalue']) + '\n' + '-'*72 + '\n'

        for method, events in self.dict.items():
            key = aliases[method]['key']
            for e in events:
                if method == 'geneconv':
                    outstr += f"{key:<8}\t{e[2][0]}\t{e[2][1]}\t{e[0]:<11}\t" \
                              f"{e[1][0]:<7}\t{e[1][1]:<7}\t{float(e[3]):.2E}\n"
                else:
                    outstr += f"{key:<8}\t{e[2]}\t{e[3]}\t{e[0]:<11}\t" \
                              f"{e[1][0]:<7}\t{e[1][1]:<7}\t{float(e[4]):.2E}\n"
        return outstr

    def __getitem__(self, key):
        events = self.dict.get(key, None)
        if key == 'geneconv':
            return [{'start': e[2][0], 'end': e[2][1], 'recombinant': e[0],
                     'parent1': e[1][0], 'parent2': e[1][1], 'pvalue': float(e[3])}
                    for e in events]
        else:
            return [{'start': e[2], 'end': e[3], 'recombinant': e[0],
                     'parent1': e[1][0], 'parent2': e[1][1], 'pvalue': float(e[4])}
                    for e in events]

    def keys(self):
        return self.dict.keys()


class Scanner:
    def __init__(self, cfg=None, methods=None, quiet=True):
        """
        :param cfg:  str, path to configuration file.  Defaults to None, causing
                     each method to use default settings.
        :param methods:  tuple, names of methods to setup and run
        :param quiet:  bool, if True, suppress console messages
        """
        # Check that the OS is valid
        sp = sys.platform
        if sp not in ['win32', 'cygwin', 'darwin', 'linux']:
            print(f"Error: binaries do not support {sp} - please contact "
                  f"the developers at https://github.com/PoonLab/OpenRDP.")
            sys.exit()

        if methods is None:
            methods = list(aliases.keys())
        self.methods = methods
        self.quiet = quiet

        self.config = configparser.ConfigParser()
        self.cfg_file = cfg
        if self.cfg_file is None:
            # load default configuration from package file
            self.cfg_file = os.path.join(basepath[0], 'default.ini')
        self.print(f"Loading configuration from {self.cfg_file}")
        self.config.read(self.cfg_file)

        self.seq_names = []
        self.alignment = None  # np.array
        self.ref_names = []
        self.ref_align = None  # np.array

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

    def _import_data(self, infile, is_ref=False):
        """
        Import labels and sequences from FASTA file and do some quality control.
        Stores sequences as a character matrix.
        :param infile:  str or File, input FASTA
        :param is_ref:  bool, specifies if infile is a query sequence or a reference
        """
        if type(infile) == str:
            if not os.path.exists(infile):
                print(f"Error: No file found at path {infile}")
                sys.exit(1)
            with open(infile) as handle:
                names, aln = read_fasta(handle)
        elif hasattr(infile, "read"):
            names, aln = read_fasta(infile)
        else:
            print(f"Error: Scanner.run_scans() must be called with a file path or File")
            sys.exit(1)

        # validate alignment
        seqlens = [len(s) for s in aln]
        if len(set(seqlens)) > 1:
            print(f"Error: {infile} does not appear to contain a valid alignment!")
            sys.exit(1)

        # Remove identical sequences
        unique = OrderedDict()  # remembers order that entries were added
        for label, seq in zip(names, aln):
            if seq not in unique:
                # validate sequence
                charset = set(seq)
                invalid = charset.difference(DNA_ALPHABET)
                if invalid:
                    print(f"Alignment contains invalid characters {''.join(invalid)}.")
                    print(f"Sequences can only contain {','.join(DNA_ALPHABET)}.")
                    sys.exit(1)
                unique.update({seq: []})
            unique[seq].append(label)

        new_aln = []
        if is_ref:
            self.ref_names = []
        else:
            self.seq_names = []
        for seq, labels in unique.items():
            if is_ref:
                self.ref_names.append(labels[0])
            else:
                self.seq_names.append(labels[0])
            if len(labels) > 1:
                for label in labels[1:]:
                    self.print(f"{label} is a duplicate of {labels[0]}")
            new_aln.append(seq)

        # Create an m x n array of sequences
        alignment = np.array(list(map(list, new_aln)))
        if is_ref:
            self.ref_align = alignment
        else:
            self.alignment = alignment

    def run_scans(self, infile, ref_file):  # pragma: no cover
        """
        Run the selected recombination detection analyses
        :param infile:  str, path to input FASTA file
        :param ref_file:  str, path to input FASTA reference file
        """
        # prepare return value
        results = ScanResults(dict([(method, {}) for method in aliases.keys()]))

        # Run methods with external binaries
        if 'threeseq' in self.methods:
            three_seq = ThreeSeq(infile)
            self.print("Starting 3Seq Analysis")
            results.dict['threeseq'] = three_seq.execute()
            self.print("Finished 3Seq Analysis")

        if 'geneconv' in self.methods:
            # Parse output file if available
            if self.config:
                geneconv = GeneConv(settings=dict(self.config.items('Geneconv')))
            else:
                geneconv = GeneConv()  # default config

            self.print("Starting GENECONV Analysis")
            results.dict['geneconv'] = geneconv.execute(infile)
            self.print("Finished GENECONV Analysis")

        # Exit early if 3Seq and Geneconv are the only methods selected
        check = set(aliases.keys()).intersection(self.methods)
        if len(check) == 0:
            return results

        # Run internal methods
        self._import_data(infile)  # sets seq_names and alignment
        if ref_file:
            self._import_data(ref_file, True)

        tmethods = {}
        for alias, a in aliases.items():
            if alias in ['threeseq', 'geneconv'] or alias not in self.methods:
                continue

            self.print(f"Setting up {alias} analysis...")
            if self.config:
                settings = dict(self.config.items(a['key']))
                tmethod = a['method'](self.alignment, settings=settings, quiet=self.quiet,
                                      ref_align=self.ref_align if ref_file else None)
            else:
                tmethod = a['method'](self.alignment, quiet=self.quiet,
                                      ref_align=self.ref_align if ref_file else None)
            tmethods.update({alias: tmethod})

        # iterate over all triplets in the alignment
        if ref_file:
            total_num_trps = self.alignment.shape[0] * int(ncomb(self.ref_align.shape[0], 2))
        else:
            total_num_trps = int(ncomb(self.alignment.shape[0], 3))
        if 'bootscan' in tmethods:
            bootscan = tmethods['bootscan']
            bootscan.execute_all(total_combinations=total_num_trps, seq_names=self.seq_names,
                                 ref_names=self.ref_names if ref_file else None)
            if os.path.exists(bootscan.dt_matrix_file):
                os.remove(bootscan.dt_matrix_file)

        triplets = TripletGenerator(self.alignment, self.seq_names,
                                    ref_align=self.ref_align if ref_file else None,
                                    ref_names=self.ref_names if ref_file else None)
        for trp_count, triplet in enumerate(triplets):
            self.print("Scanning triplet {} / {}".format(trp_count + 1, total_num_trps))
            for alias, tmethod in tmethods.items():
                if alias == 'bootscan':
                    continue
                tmethod.execute(triplet)

        # Process results by joining breakpoint locations that overlap
        for alias, tmethod in tmethods.items():
            results.dict[alias] = tmethod.merge_breakpoints()

        return results
