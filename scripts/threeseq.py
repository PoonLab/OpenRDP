import subprocess
import os
import glob
import sys
import logging


class ThreeSeq:
    def __init__(self, in_path):
        self.in_path = os.path.realpath(in_path)
        self.in_name = os.path.basename(in_path)

    def execute(self):
        """
        Execute the 3Seq algorithm.
            Lam HM, Ratmann O, Boni MF.
            Improved algorithmic complexity for the 3SEQ recombination detection algorithm.
            Mol Biol Evol, 35(1):247-251, 2018.
        :return: A list containing the results of the 3Seq analysis
                    Format: [triplets, uncorrected p-value, corrected p-value, breakpoint locations]
        """
        # Clear output files
        out_files = glob.glob('./*.3s.*')
        for f in out_files:
            try:
                os.remove(f)
            except OSError:
                pass

        # Set paths to 3Seq executables
        if sys.platform.startswith("win"):
            bin_path = os.path.abspath('../bin/3Seq/windows_3seq.exe')
        else:
            bin_path = os.path.abspath('../bin/3Seq/unix_3seq.exe')

        if not os.path.isfile(bin_path):
            logging.error("No 3Seq executable file exists.")

        # Run 3Seq
        if sys.platform.startswith("win"):
            try:
                subprocess.check_output([bin_path, "-f", self.in_path, "-pTable myPvalueTable", "-d", "-id", self.in_name],
                                        shell=False, input=b"Y\n")  # Respond to prompt
            except subprocess.CalledProcessError as e:
                print(e.output, e.returncode)
        else:
            try:
                subprocess.check_output([bin_path, "-f", self.in_path, "-pTable myPvalueTable", "-d", "-id", self.in_name],
                                        shell=False, input=b"Y\n")  # Respond to prompt
            except subprocess.CalledProcessError as e:
                print(e.output, e.returncode)

        # Parse the output of 3Seq
        tseq_out = self.in_name + '.3s.rec'
        out_path = os.path.join(os.getcwd(), tseq_out)
        ts_results = self.parse_output(out_path)

        return ts_results

    def parse_output(self, out_path):
        """
        Parse the output of the 3Seq analysis
        :param out_path: Path to the output file containing information about recombinant sequences
        :return: List of triplets, corrected and uncorrected p-values, and breakpoint locations
        """
        ts_results = []

        # Check that the out file exists
        try:
            with open(out_path) as out_handle:
                out_handle.readline()  # Read first line

                for line in out_handle:
                    line = line.split('\t')
                    line = [l.strip() for l in line]
                    triplet = ([line[0], line[1], line[2]])  # Record the triplets
                    uncorr_p_value = line[6]  # Uncorrected p-value
                    corr_p_value = line[10]  # Dunn-Sidak corrected p-value
                    locations = line[12:]  # Breakpoint locations
                    ts_results.append([triplet, uncorr_p_value, corr_p_value, locations])
        except FileNotFoundError as e:
            print(e)

        return ts_results
