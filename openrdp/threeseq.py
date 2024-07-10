import logging
import os
import subprocess
import sys
from tempfile import NamedTemporaryFile
from openrdp.common import merge_breakpoints

class ThreeSeq:
    def __init__(self, in_path, verbose=False):
        self.in_path = os.path.realpath(in_path)  # input FASTA file
        self.in_name = os.path.basename(in_path)
        self.raw_results = []
        self.results = []
        self.name = 'threeseq'
        self.binaries = {
            'win32': 'windows_3seq.exe',
            'cygwin': 'windows_3seq.exe',
            'darwin': '3seq.macOS',
            'linux': '3seq.Unix'
        }

    def execute(self):
        """
        Execute the 3Seq algorithm.
            Lam HM, Ratmann O, Boni MF.
            Improved algorithmic complexity for the 3SEQ recombination detection algorithm.
            Mol Biol Evol, 35(1):247-251, 2018.
        :return: A list containing the results of the 3Seq analysis
                    Format: [triplets, uncorrected p-value, corrected p-value, breakpoint locations]
        """

        # Set paths to 3Seq executables
        bin_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'bin')
        bin_path = os.path.join(bin_dir, '3Seq', self.binaries[sys.platform])
        if not os.path.isfile(bin_path):
            logging.error("No 3Seq executable file exists.")

        with NamedTemporaryFile(delete=False) as tempf:
            cmd = [
                bin_path,
                "-full", self.in_path,  # process entire sequences
                "-ptable myPvalueTable",  # included with 3Seq
                "-d",  # distinct sequences only
                "-id", tempf.name  # write outputs to temporary file
            ]
            subprocess.check_output(cmd, shell=False, input=b"Y\n")  # Respond to prompt

            # Parse the output of 3Seq
            out_path = tempf.name + '.3s.rec'
            if not os.path.exists(out_path):
                out_path += '.csv'
                if not os.path.exists(out_path):
                    # No recombinant triplets were written to this output file
                    return []
            ts_results = self.parse_output(out_path)

        return ts_results

    def parse_output(self, out_path):
        """
        Parse the output of the 3Seq analysis
        :param out_path: Path to the output file containing information about recombinant sequences
        :return: List of triplets, corrected and uncorrected p-values, and breakpoint locations
        """
        delimiter = '\t'
        if out_path.endswith('.csv'):
            delimiter = ','

        with open(out_path) as out_handle:
            out_handle.readline()  # skip first line
            for line in out_handle:
                values = [item.strip() for item in line.split(delimiter)]
                rec = values[0]
                ps = [values[1], values[2]]
                corr_p_value = values[10]  # Dunn-Sidak corrected p-value

                for loc in values[12:]:  # Breakpoint locations
                    parts = loc.split(' & ')
                    # Take the widest interval 3Seq returns
                    start_pos = parts[0].split('-')
                    end_pos = parts[1].split('-')
                    self.raw_results.append((rec, ps, start_pos[0], end_pos[-1], corr_p_value))

        self.results = merge_breakpoints(self.raw_results) # uses no pvalue cut off, set to default 100
        return self.results

