import logging
import os
import subprocess
import sys
from tempfile import NamedTemporaryFile


class ThreeSeq:
    def __init__(self, in_path):
        self.in_path = os.path.realpath(in_path)  # input FASTA file
        self.in_name = os.path.basename(in_path)
        self.raw_results = []
        self.results = []
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
        bin_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, 'bin')
        bin_path = os.path.join(bin_dir, '3Seq', self.binaries[sys.platform])
        if not os.path.isfile(bin_path):
            logging.error("No 3Seq executable file exists.")

        with NamedTemporaryFile(delete=False) as tempf:
            subprocess.check_output([
                bin_path,
                "-full", self.in_path,  # process entire sequences
                "-ptable myPvalueTable",
                "-d",  # distinct sequences only
                "-id", tempf.name  # write outputs to temporary file
            ], shell=False, input=b"Y\n")  # Respond to prompt

            # Parse the output of 3Seq
            out_path = tempf.name + '.3s.rec.csv'
            if not os.path.exists(out_path):
                out_path = tempf.name + '.3s.rec'
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

                loc_line = values[12:]    # Breakpoint locations
                for loc in loc_line:
                    parts = loc.split(' & ')
                    # Take the widest interval 3Seq returns
                    start_pos = parts[0].split('-')
                    end_pos = parts[1].split('-')
                    self.raw_results.append((rec, ps, start_pos[0], end_pos[-1], corr_p_value))

        self.results = self.merge_breakpoints()
        return self.results

    def merge_breakpoints(self):
        """
        Merge overlapping breakpoint locations
        :return: list of breakpoint locations where overlapping intervals are merged
        """
        results_dict = {}
        results = []

        # Gather all regions with the same recombinant
        for i, bp in enumerate(self.raw_results):
            rec_name = self.raw_results[i][0]
            parents = tuple(sorted(self.raw_results[i][1]))
            key = (rec_name, parents)
            if key not in results_dict:
                results_dict[key] = []
            results_dict[key].append(self.raw_results[i][2:])

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
                        region[0] = min(start,start2)
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
                results.append((rec_name, parents, start, end, p_value))

        return results

