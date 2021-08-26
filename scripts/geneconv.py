import subprocess
import os
import sys
import logging
import glob


class GeneConv:
    def __init__(self, settings=None, gscale=1, ignore_indels=False, min_length=1, min_poly=2, min_score=2,
                 max_overlap=1):
        """
        Constructs a GeneConv object
        :param gscale: mismatch penalty
        :param ignore_indels: Ignore indels or treat indels as one polymorphism (default = False)
        :param min_length: Minimum length of the fragments
        :param min_poly: Minimum number of polymorphic sites
        :param min_score: Minimum pairwise score
        :param max_overlap: Maximum number of overlapping pairs
        """
        if settings is not None:
            self.set_options_from_config(settings)
            self.validate_options()

        else:
            self.gscale = gscale
            self.ignore_indels = ignore_indels
            self.min_length = min_length
            self.min_poly = min_poly
            self.min_score = min_score
            self.max_overlap = max_overlap

    def set_options_from_config(self, settings):
        """
        Set the parameters of GENECONV from the config file
        :param settings: a dictionary of settings
        """
        self.gscale = int(settings['mismatch_penalty'])

        if settings['indels_as_polymorphisms'] == 'True':
            self.ignore_indels = False
        elif settings['indels_as_polymorphisms'] == 'False':
            self.ignore_indels = True
        else:
            self.ignore_indels = None

        self.min_length = int(settings['min_len'])
        self.min_poly = int(settings['min_poly'])
        self.min_score = int(settings['min_score'])
        self.max_overlap = int(settings['max_num'])

    def validate_options(self):
        """
        Check if the options from the config file are valid
        If the options are invalid, the default value will be used instead
        """
        if not isinstance(self.ignore_indels, bool):
            print("Invalid option for 'indels_as_polymorphisms'.\nUsing default value (False) instead.\n")
            self.ignore_indels = False

        if self.min_length <= 0:
            print("Invalid option for 'min_len'.\nUsing default value (1) instead.")
            self.min_length = 1

        if self.min_poly < 1:
            print("Invalid option for 'min_score'.\nUsing default value (2) instead.")
            self.min_poly = 2

        if self.max_overlap < 0:
            print("Invalid option for 'max_num'.\nUsing default value (1) instead.")
            self.max_overlap = 1

    def execute(self, in_path):
        """
        Execute the GENECONV algorithm
            S. A. Sawyer (1999)
            GENECONV: A computer package for the statistical detection of gene conversion.
            Distributed by the author, Department of Mathematics, Washington University in St. Louis,
            Available at http://www.math.wustl.edu/~sawyer.
        :param in_path: Path to the input alignment file
        :return: A list of results
        """
        # Clear output files
        out_files = glob.glob('../bin/GENECONV/*.frags') + glob.glob('../bin/GENECONV/*.sum')
        for f in out_files:
            try:
                os.remove(f)
            except OSError:
                pass

        # Create config file
        with open("geneconv.cfg", 'w+') as cfg_handle:
            cfg_handle.write('#GCONV_CONFIG\n')
            cfg_handle.write('  -inputpath={}\n'.format(os.path.realpath(in_path)))

            if not self.ignore_indels:
                cfg_handle.write('  -Indel_blocs\n')

            cfg_handle.write('  -Gscale={}\n'.format(self.gscale))
            cfg_handle.write('  -Minlength={}\n'.format(self.min_length))
            cfg_handle.write('  -Minnpoly={}\n'.format(self.min_poly))
            cfg_handle.write('  -Minscore={}\n'.format(self.min_score))
            cfg_handle.write('  -Maxoverlap={}\n'.format(self.max_overlap))

            cfg_path = format(os.path.realpath("geneconv.cfg"))

        # Path to GENECONV executables
        if sys.platform.startswith("win"):
            bin_path = os.path.abspath('../bin/GENECONV/windows_geneconv.exe')
        else:
            bin_path = os.path.abspath('../bin/GENECONV/unix_geneconv.exe')

        if not os.path.isfile(bin_path):
            logging.error("No GENECONV executable file exists")

        # Run GENECONV
        if sys.platform.startswith("win"):
            try:
                subprocess.check_output([bin_path, "-Seqfile={}".format(in_path),
                                         "-Config={}".format(cfg_path), "-nolog"], shell=False)
            except subprocess.CalledProcessError as e:
                print(e.output, e.returncode)
        else:
            try:
                subprocess.check_output([bin_path, "-Seqfile={}".format(in_path),
                                         "-Config={}".format(cfg_path), "-nolog"], shell=False)
            except subprocess.CalledProcessError as e:
                print(e.output, e.returncode)

        # Remove GENECONV config file
        try:
            os.remove(cfg_path)
        except OSError:
            pass

        # Parse the output of 3Seq
        in_name = os.path.basename(in_path).split('.')[0]
        out_name = in_name + '.frags'
        out_path = os.path.join(os.path.dirname(in_path), out_name)
        gc_results = self.parse_output(out_path)

        return gc_results

    def parse_output(self, out_path):
        """
        Parse the output of the GENECONV analysis
        :param out_path: Path to the output file
        :return: List of results
        """
        gc_results = {}

        # Check that the out file exists
        try:
            with open(out_path) as out_handle:

                for line in out_handle:
                    if not line.startswith('#'):
                        line = line.strip()
                        line = line.split()
                        seqs = line[1]  # Sequence pairs
                        uncorr_p_value = line[2]  # Uncorrected p_value
                        corr_p_value = line[3]  # Bonferroni Corrected - Karlin-Altschul
                        locations = (line[4], line[5])  # Locations in alignment
                        type = line[0]  # Global inner (GI), global outer (GO), additional inner (AI)
                        gc_results[locations] = [seqs, uncorr_p_value, corr_p_value, type]

        except FileNotFoundError as e:
            print(e)

        return gc_results
