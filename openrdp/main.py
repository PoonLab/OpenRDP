import argparse
import sys
from openrdp import run_scans
from datetime import datetime

DNA_ALPHABET = ['A', 'T', 'G', 'C', '-', '*']


def valid_alignment(alignment):
    """
    Check that the input alignment is valid
    :param alignment: a list of lists containing the sequence headers and the aligned sequences
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
    :param alignment: a list of lists containing the sequence headers and the aligned sequences
    :return: True if the alignment contains only valid characters, False otherwise
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
    :return: tuple of headers and sequences
    """
    result = []
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
    # Check that the OS is valid
    platform = sys.platform
    try:
        platform.startswith("win") or platform.startswith("linux") or sys.platform == "darwin"
    except OSError:
        print("OSError: {} is not supported".format(sys.platform))

    with open(infile) as in_handle:
        if infile.endswith('.fa') or infile.endswith('.fasta'):
            names, aln = read_fasta(in_handle)

    if not valid_alignment(aln) and not valid_chars(aln):
        sys.exit(1)

    scanner = run_scans.Scanner(names, infile, outfile, cfg=cfg, methods=methods, quiet=quiet)
    scanner.run_scans(aln)

