import argparse
import sys
from do_scans import *
from run_scans import Scanner

DNA_ALPHABET = ['A', 'T', 'G', 'C', '-', '*']


def valid_arguments(alignment):
    """
    Check that the input alignment is valid
    :param alignment: a list of lists containing the sequence headers and the aligned sequences
    :return True of the alignment is valid, false otherwise
    """
    aln_len = len(alignment[0][1])
    for pair in alignment:
        if len(pair[1]) != aln_len:
            print("Improper alignment! Not all alignments are the same length.")
            return False
    return True


def valid_chars(alignment):
    for aln in alignment:
        if not all(pos in DNA_ALPHABET for pos in aln.upper()):
            print("Alignment contains invalid characters.")
            return False
    return True


def convert_fasta(handle):
    result = []
    sequence, h = '', ''

    # Verifies files have the correct formatting
    for i, line in enumerate(handle):
        if line.startswith('$'):
            continue
        elif line.startswith('>') or line.startswith('#'):
            break
        else:
            print("No header")
            raise NameError

    # Reset pointer to beginning of file
    if hasattr(handle, 'seek'):
        handle.seek(0)

    for line in handle:
        if line.startswith('$'):
            continue
        elif line.startswith('>') or line.startswith('#'):
            if len(sequence) > 0:
                result.append([h, sequence])
                sequence = ''
            h = line.strip('>#\t\n\r')
        else:
            sequence += line.strip('\n\r').upper()

    result.append([h, sequence])
    return result


def parse_args():
    parser = argparse.ArgumentParser(
        description='An open source implementation of RDP5.\n'
    )

    parser.add_argument('infile',
                        help='File containing sequence alignment (FASTA or CLUSTAL) format')

    parser.add_argument('-cfg',
                        help='Path to file that contains parameters')

    parser.add_argument('-geneconv',
                        help='Perform GeneConv analysis',
                        action='store_true')

    parser.add_argument('-bootscan',
                        help='Perform Bootscan analysis',
                        action='store_true')

    parser.add_argument('-maxchi',
                        help='Perform MaxChi analysis',
                        action='store_true')

    parser.add_argument('-siscan',
                        help='Perform Siscan analysis',
                        action='store_true')

    parser.add_argument('-chimaera',
                        help='Perform Chimaera analysis',
                        action='store_true')

    parser.add_argument('-threeseq',
                        help='Perform 3Seq analysis',
                        action='store_true')

    parser.add_argument('-rdp',
                        help='{Perform RDP analysis',
                        action='store_true')

    return parser.parse_args()


def main():
    args = parse_args()

    with open(args.infile) as in_handle:
        if args.infile.endswith('.fa') or args.infile.endswith('.fasta'):
            aln = convert_fasta(in_handle)

    if not valid_arguments(aln) and not valid_arguments(aln):
        sys.exit(1)

    # Check that the OS is valid
    platform = sys.platform
    try:
        platform.startswith("win") or platform.startswith("linux") or sys.platform == "darwin"
    except OSError:
        print("OSError: {} is not supported".format(sys.platform))

    scanner = Scanner(aln, args)
    results = scanner.run_scans()

    # print(results)


if __name__ == '__main__':
    main()
