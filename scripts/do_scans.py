from utils import geneconv

import os
from utils.geneconv import run_geneconv


class RdpMethod:
    """
    Executes command for RDP method
    """

    def __init__(self, aln, win_size, automask):
        self._aln = aln
        self._win_size = win_size
        self._automask = automask

    def do_rdp(self):
        pass


class GeneConvMethod:
    def __init__(self, gscale, ignore_indels, min_length, min_poly, min_score, max_overlap):
        self.gscale = gscale                # Mismatch penalty (default = 1)
        self.ignore_indels = ignore_indels  # Ignore indels or treat indels as one polymorphism (default = False)
        self.min_length = min_length        # Default = 1
        self.min_poly = min_poly            # Default = 2
        self.min_score = min_score          # Default = 2
        self.max_overlap = max_overlap      # Default = 1

    def make_cfg(self):
        with open("geneconv.cfg", 'w+') as cfg_handle:
            cfg_handle.write('#GCONV_CONFIG\n')
            cfg.handle.write('  -inputpath={}\n'.format(os.path.realpath(cfg_handle.name)))

            if not self.ignore_indels:
                cfg_handle.write('  -Indel_blocs\n')

            cfg_handle.write('  -Gscale={}\n'.format(self.gscale))
            cfg_handle.write('  -Minlength={}\n'.format(self.min_length))
            cfg_handle.write('  -Minnpoly={}\n'.format(self.min_poly))
            cfg_handle.write('  -Minscore={}\n'.format(self.min_score))
            cfg_handle.write('  -Maxoverlap={}\n'.format(self.max_overlap))

            return cfg_handle.name

    def do_geneconv(self):
        cfg_path = self.make_cfg()
        raw_output = run_geneconv(cfg_path)


class Scan_Invoker:
    """
    Receives the command for the scan and invokes the appropriate receiver
    """
    def __init__(self, geneconv, bootscan, maxchi, siscan, chimaera, seq3):
        self._method_name = ''

    def set_name(self, method_name):
        self._method_name = method_name

    @staticmethod
    def on_start(method_name):
        print('Starting {}'.format(method_name))

    @staticmethod
    def on_finish(method_name):
        print('Finished {}'.format(method_name))

