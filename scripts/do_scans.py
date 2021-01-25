import os
import sys
import logging
import subprocess


class RdpMethod:
    """
    Executes command for RDP method
    """
    def __init__(self, win_size=30, reference=None):
        self._win_size = win_size
        self.reference = reference

    def load_config(self):
        pass

    def execute(self, data_path):
        pass


class GeneConv:
    def __init__(self, gscale=1, ignore_indels=False, min_length=1, min_poly=2, min_score=2, max_overlap=1):
        self.gscale = gscale                # Mismatch penalty
        self.ignore_indels = ignore_indels  # Ignore indels or treat indels as one polymorphism (default = False)
        self.min_length = min_length
        self.min_poly = min_poly
        self.min_score = min_score
        self.max_overlap = max_overlap

    def set_options_from_config(self, settings={}):
        pass

    def validate_options(self):
        pass

    def execute(self, data_path):

        # Check if valid options
        if not self.validate_options():
            return None

        # Create config file
        with open("geneconv.cfg", 'w+') as cfg_handle:
            cfg_handle.write('#GCONV_CONFIG\n')
            cfg_handle.write('  -inputpath={}\n'.format(os.path.realpath(data_path.name)))

            if not self.ignore_indels:
                cfg_handle.write('  -Indel_blocs\n')

            cfg_handle.write('  -Gscale={}\n'.format(self.gscale))
            cfg_handle.write('  -Minlength={}\n'.format(self.min_length))
            cfg_handle.write('  -Minnpoly={}\n'.format(self.min_poly))
            cfg_handle.write('  -Minscore={}\n'.format(self.min_score))
            cfg_handle.write('  -Maxoverlap={}\n'.format(self.max_overlap))

            cfg_path = format(os.path.realpath(cfg_handle.name))

        # Path to GENECONV executables
        script_path = os.path.dirname(os.path.abspath(__file__))

        if sys.platform.startswith("win"):
            bin_path = os.path.join(script_path, 'bin/GENECONV/windows_geneconv.exe')
        else:
            bin_path = os.path.join(script_path, 'bin/GENECONV/unix_geneconv.exe')

        if not os.path.isfile(bin_path):
            logging.error("No file exists")

        # Run GENECONV
        if sys.platform.startswith("win"):
            gc_output = subprocess.check_output([bin_path, cfg_path], shell=False, stderr=subprocess.DEVNULL)
        else:
            gc_output = subprocess.check_output([bin_path, cfg_path], shell=False, stderr=subprocess.STDOUT)

        # Remove file
        os.remove(cfg_path)

        return gc_output


class ThreeSeq:
    def __init__(self, p_value_table=None):
        self.p_value_table = p_value_table

    def load_config(self):
        pass

    def execute(self, data_path):
        # Path to 3Seq executables
        script_path = os.path.dirname(os.path.abspath(__file__))

        if sys.platform.startswith("win"):
            bin_path = os.path.join(script_path, 'bin/3Seq/windows_3seq.exe')
        else:
            bin_path = os.path.join(script_path, 'bin/GENECONV/unix_3seq.exe')

        if not os.path.isfile(bin_path):
            logging.error("No file exists")

        # Run 3Seq
        if sys.platform.startswith("win"):
            tseq_output = subprocess.check_output([bin_path, '-f', os.path.realpath(data_path.name)],
                                                  shell=False, stderr=subprocess.DEVNULL)
        else:
            tseq_output = subprocess.check_output([bin_path, '-f', os.path.realpath(data_path.name)],
                                                  shell=False, stderr=subprocess.STDOUT)

        return tseq_output
