import argparse
import sys
from rdp import Sequence
from rdp import Alignment

DNA_ALPHABET = ['A', 'T', 'G', 'C', '-', '*']


def valid_arguments(alignment):
    aln_length = len(alignment[0])
    for aln in alignment:
        if len(aln) != aln_length:
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
    result = {}
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
                result[h] = sequence
                sequence = ''
            h = line.strip('>#\n\r')
        else:
            sequence += line.strip('\n\r').upper()

    result[h] = sequence
    return result


def convert_clustal(handle):
    result = {}
    # Loop over each sequence in the clustal alignment
    for ln in handle[3]:
        # Skip first 3 lines
        if len(ln) > 0:
            # Match sequence label and build up sequence
            if ln[0].isalpha():
                line = ln.split()
                header = line[0].strip()
                seq = line[1].strip()
                result.setdefault(header, {})
                result[header] += seq

    return result


def parse_args():
    parser = argparse.ArgumentParser(
        description='An open source implementation of RDP4.\n'
    )

    parser.add_argument('infile',
                        help='File containing sequence alignment (FASTA or CLUSTAL) format')

    parser.add_argument('-am', '-o',
                        help='Optimize auto-masking for maximum recombination detection')

    parser.add_argument('-winsize',
                        help='The window size used for the RDP analysis')

    # TODO: Add arguments for GENECONV, Bootscan, MaxChi, Chimaera, Siscan, 3Seq

    return parser.parse_args()


def main():
    args = parse_args()

    with open(args.infile) as in_handle:
        if args.infile.endswith('.fa') or args.infile.endswith('.fasta'):
            aln = convert_fasta(in_handle)
        elif args.infile.endswith('.aln'):
            aln = convert_clustal(in_handle)

    if not valid_arguments(aln) and not valid_arguments(aln):
        sys.exit(1)

    # Create sequence objects
    seqs_in_aln = []
    for a in aln:
        seq = Sequence(a, aln[a])
        seqs_in_aln.append(seq)

    # Create alignment
    alignment = Alignment(seqs_in_aln)


if __name__ == '__main__':
    main()