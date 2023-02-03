import sys
import configparser
import numpy as np
from itertools import combinations

from openrdp.bootscan import Bootscan
from openrdp.chimaera import Chimaera
from openrdp.common import generate_triplets, Triplet
from openrdp.geneconv import GeneConv
from openrdp.maxchi import MaxChi
from openrdp.rdp import RdpMethod
from openrdp.siscan import Siscan
from openrdp.threeseq import ThreeSeq


# keywords for recombination detection methods implemented in this package
all_methods = ['geneconv', 'bootscan', 'maxchi', 'siscan', 'chimaera',
               'threeseq', 'rdp']


class Scanner:
    def __init__(self, names, infile, outfile, cfg, methods=None, quiet=False):
        self.seq_names = names
        self.infile = infile
        self.cfg_file = cfg
        self.methods = methods
        self.quiet = quiet
        self.outfile = outfile

    def run_scans(self, aln):
        """
        Run the selected recombination detection analyses
        :param aln:  list, sequences imported from FASTA file
        """
        # Parse config file
        if self.cfg_file:
            config = configparser.ConfigParser()
            config.read(self.cfg_file)
        else:
            config = None

        # Remove identical sequences
        aln = list(set(aln))

        # check that all sequences are the same length
        seqlen = set([len(s) for s in aln])
        if len(seqlen) > 1:
            print("ERROR: The input sequences are different lengths!  Input must be aligned.")
            sys.exit()

        # Create an m x n array of sequences
        alignment = np.array(list(map(list, aln)))
        results = dict([(method, None) for method in methods])

        # Run 3Seq
        if 'threeseq' in self.methods:
            three_seq = ThreeSeq(self.infile)
            if not self.quiet:
                print("Starting 3Seq Analysis")

            threeseq_res = three_seq.execute()
            if not self.quiet:
                print("Finished 3Seq Analysis")

        # Run GENECONV
        if self.geneconv:
            # Parse output file if available
            if config:
                geneconv = GeneConv(settings=dict(config.items('Geneconv')))
            else:
                geneconv = GeneConv()

            if not self.quiet:
                print("Starting GENECONV Analysis")
            geneconv_res = geneconv.execute(self.infile)

            if not self.quiet:
                print("Finished GENECONV Analysis")

        # Exit early if 3Seq and Geneconv are the only methods selected
        if not self.maxchi and not self.chimaera and not self.siscan and not self.rdp and not self.bootscan:
            self.write_output(threeseq_res, geneconv_res, bootscan_res, maxchi_res, chimaera_res, rdp_res, siscan_res)
            return

        maxchi, chimaera, rdp, bootscan, siscan = None, None, None, None, None

        # Setup MaxChi
        if self.maxchi:
            if not self.quiet:
                print("Setting Up MaxChi Analysis")
            if config:
                maxchi = MaxChi(alignment, settings=dict(config.items('MaxChi')))
            else:
                maxchi = MaxChi(alignment)

        # Setup Chimaera
        if self.chimaera:
            if not self.quiet:
                print("Setting Up Chimaera Analysis")
            if config:
                chimaera = Chimaera(alignment, settings=dict(config.items('Chimaera')))
            else:
                chimaera = Chimaera(alignment)

        # Setup RDP
        if self.rdp:
            if not self.quiet:
                print("Setting Up RDP Analysis")
            if config:
                rdp = RdpMethod(alignment, settings=dict(config.items('RDP')))
            else:
                rdp = RdpMethod(alignment)

        # Setup Bootscan
        if self.bootscan:
            if not self.quiet:
                print("Setting Up Bootscan Analysis")
            if config:
                bootscan = Bootscan(alignment, settings=dict(config.items('Bootscan')), quiet=self.quiet)
            else:
                bootscan = Bootscan(alignment, quiet=self.quiet)

        # Setup Siscan
        if self.siscan:
            if not self.quiet:
                print("Setting Up Siscan Analysis")
            if config:
                siscan = Siscan(alignment, settings=dict(config.items('Siscan')))
            else:
                siscan = Siscan(alignment)

        trp_count = 1
        total_num_trps = sum(1 for _ in combinations(range(alignment.shape[0]), 3))
        for trp in generate_triplets(alignment):
            triplet = Triplet(alignment, self.seq_names, trp)
            if not self.quiet:
                print("Scanning triplet {} / {}".format(trp_count, total_num_trps))
            trp_count += 1

            # Run MaxChi
            if self.maxchi:
                maxchi.execute(triplet)
            # Run Chimaera
            if self.chimaera:
                chimaera.execute(triplet)
            # Run RDP Method
            if self.rdp:
                rdp.execute(triplet)
            # Run Bootscan
            if self.bootscan:
                bootscan.execute(triplet)
            # Run Siscan
            if self.siscan:
                siscan.execute(triplet)

        # Process results by joining breakpoint locations that overlap
        if self.maxchi:
            maxchi_res = maxchi.merge_breakpoints()
        if self.chimaera:
            chimaera_res = chimaera.merge_breakpoints()
        if self.rdp:
            rdp_res = rdp.merge_breakpoints()
        if self.bootscan:
            bootscan_res = bootscan.merge_breakpoints()
        if self.siscan:
            siscan_res = siscan.merge_breakpoints()

        self.write_output(threeseq_res, geneconv_res, bootscan_res, maxchi_res, chimaera_res, rdp_res, siscan_res)

        if not self.quiet:
            self.print_output(threeseq_res, geneconv_res, bootscan_res, maxchi_res, chimaera_res, rdp_res, siscan_res)

    def write_output(self, threeseq_res, geneconv_res, bootscan_res, maxchi_res, chimaera_res, rdp_res, siscan_res):
        """
        Write results to a file
        :param threeseq_res: results of 3Seq analysis
        :param geneconv_res: results of Geneconv analysis
        :param bootscan_res: results of Bootscan analysis
        :param maxchi_res: results of MaxChi analysis
        :param chimaera_res: results of Chimaera analysis
        :param rdp_res: results of RDP analysis
        :param siscan_res: results of Siscan analysis
        """
        with open(self.outfile, 'w+') as outfile:
            outfile.write('Method,StartLocation,EndLocation,Recombinant,Parent1,Parent2,Pvalue\n')

            if self.threeseq:
                for event in threeseq_res:
                    outfile.write('3Seq,{},{},{},{},{},{}\n'
                                  .format(event[2], event[3], event[0], event[1][0], event[1][1], event[4]))
            if self.geneconv:
                for event in geneconv_res:
                    outfile.write('Geneconv,{},{},{},{},{},{}\n'
                                  .format(event[2][0], event[2][1], event[0], event[1][0], event[1][1], event[3]))
            if self.bootscan:
                for event in bootscan_res:
                    outfile.write('Bootscan,{},{},{},{},{},{}\n'
                                  .format(event[2], event[3], event[0], event[1][0], event[1][1], event[4]))
            if self.maxchi:
                for event in maxchi_res:
                    outfile.write('MaxChi,{},{},{},{},{},{}\n'
                                  .format(event[2], event[3], event[0], event[1][0], event[1][1], event[4]))
            if self.chimaera:
                for event in chimaera_res:
                    outfile.write('Chimaera,{},{},{},{},{},{}\n'
                                  .format(event[2], event[3], event[0], event[1][0], event[1][1], event[4]))
            if self.rdp:
                for event in rdp_res:
                    outfile.write('RDP,{},{},{},{},{},{}\n'
                                  .format(event[2], event[3], event[0], event[1][0], event[1][1], event[4]))
            if self.siscan:
                for event in siscan_res:
                    outfile.write('Siscan,{},{},{},{},{},{}\n'
                                  .format(event[2], event[3], event[0], event[1][0], event[1][1], event[4]))

    def print_output(self, threeseq_res, geneconv_res, bootscan_res, maxchi_res, chimaera_res, rdp_res, siscan_res):
        """
        Print results to console
        :param threeseq_res: results of 3Seq analysis
        :param geneconv_res: results of Geneconv analysis
        :param bootscan_res: results of Bootscan analysis
        :param maxchi_res: results of MaxChi analysis
        :param chimaera_res: results of Chimaera analysis
        :param rdp_res: results of RDP analysis
        :param siscan_res: results of Siscan analysis
        """
        print('\n{:<20} {:<20} {:<20} {:<20} {:<20} {:<20} {:<20}'
              .format('Method', 'StartLocation', 'EndLocation', 'Recombinant', 'Parent1', 'Parent2', 'Pvalue'))

        if self.threeseq:
            for event in threeseq_res:
                print('{:<20} {:<20} {:<20} {:<20} {:<20} {:<20} {:<20}'
                      .format('3Seq', event[2], event[3], event[0], event[1][0], event[1][1], event[4]))
        if self.geneconv:
            for event in geneconv_res:
                print('{:<20} {:<20} {:<20} {:<20} {:<20} {:<20} {:<20}'
                      .format('Geneconv', event[2][0], event[2][1], event[0], event[1][0], event[1][1], event[3]))
        if self.bootscan:
            for event in bootscan_res:
                print('{:<20} {:<20} {:<20} {:<20} {:<20} {:<20} {:<20}'
                      .format('Bootscan', event[2], event[3], event[0], event[1][0], event[1][1], event[4]))
        if self.maxchi:
            for event in maxchi_res:
                print('{:<20} {:<20} {:<20} {:<20} {:<20} {:<20} {:<20}'
                      .format('MaxChi', event[2], event[3], event[0], event[1][0], event[1][1], event[4]))
        if self.chimaera:
            for event in chimaera_res:
                print('{:<20} {:<20} {:<20} {:<20} {:<20} {:<20} {:<20}'
                      .format('Chimaera', event[2], event[3], event[0], event[1][0], event[1][1], event[4]))
        if self.rdp:
            for event in rdp_res:
                print('{:<20} {:<20} {:<20} {:<20} {:<20} {:<20} {:<20}'
                      .format('RDP', event[2], event[3], event[0], event[1][0], event[1][1], event[4]))
        if self.siscan:
            for event in siscan_res:
                print('{:<20} {:<20} {:<20} {:<20} {:<20} {:<20} {:<20}'
                      .format('Siscan', event[2], event[3], event[0], event[1][0], event[1][1], event[4]))


DNA_ALPHABET = ['A', 'T', 'G', 'C', '-', '*']


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
    :return:  FIXME: right now, nothing!
    """
    # Check that the OS is valid
    platform = sys.platform
    try:
        platform.startswith("win") or platform.startswith("linux") or sys.platform == "darwin"
    except OSError:
        print("OSError: {} is not supported".format(sys.platform))

    if not (infile.endswith('.fa') or infile.endswith('.fasta')):
        print(f"Expected '.fa' or '.fasta' suffix for input FASTA {infile}, ignoring")

    with open(infile) as in_handle:
        names, aln = read_fasta(in_handle)

    if not valid_alignment(aln) and not valid_chars(aln):
        sys.exit(1)

    scanner = Scanner(names, infile, outfile, cfg=cfg, methods=methods, quiet=quiet)
    scanner.run_scans(aln)

