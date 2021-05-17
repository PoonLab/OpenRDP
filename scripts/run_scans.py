import configparser
import numpy as np
import json

from scripts.geneconv import GeneConv
from scripts.maxchi import MaxChi
from scripts.threeseq import ThreeSeq
from scripts.chimaera import Chimaera
from scripts.rdp import RdpMethod
from scripts.siscan import Siscan
from scripts.bootscan import Bootscan


class Scanner:
    def __init__(self, aln, names, infile, cfg, run_geneconv=False, run_three_seq=False, run_rdp=False,
                 run_siscan=False, run_maxchi=False, run_chimaera=False, run_bootscan=False):
        self.aln = aln
        self.seq_names = names
        self.infile = infile
        self.cfg_file = cfg
        self.geneconv = run_geneconv
        self.threeseq = run_three_seq
        self.rdp = run_rdp
        self.siscan = run_siscan
        self.maxchi = run_maxchi
        self.chimaera = run_chimaera
        self.bootscan = run_bootscan

    def run_scans(self):
        """
        Run the selected recombination detection analyses
        """
        # Parse config file
        if self.cfg_file:
            config = configparser.ConfigParser()
            config.read(self.cfg_file)
        else:
            config = None

        # Create an m x n array of sequences (n = length, m = number of sequences)
        alignment = np.array(list(map(list, self.aln)))

        with open('results.txt', 'w+') as outfile:
            # Run 3Seq
            if self.threeseq:
                three_seq = ThreeSeq(self.infile)
                print("Staring 3Seq Analysis")
                ts_results = three_seq.execute()
                print("Finished 3Seq Analysis")
                print(ts_results)
                outfile.write('3Seq\n')
                for res in ts_results:
                    outfile.write(json.dumps(res))
                    outfile.write('\n')

            # Run GENECONV
            if self.geneconv:
                if config:
                    geneconv = GeneConv(settings=config['Geneconv'])
                else:
                    geneconv = GeneConv()
                print("Starting GENECONV Analysis")
                gc_results = geneconv.execute(self.infile)

                print("Finished GENECONV Analysis")
                print(gc_results)

                outfile.write('Geneconv\n')
                for res in gc_results:
                    outfile.write(json.dumps(res))
                    outfile.write('\n')

        # Exit early if 3Seq and Geneconv are the only methods selected
        if not self.maxchi and not self.chimaera and not self.siscan and not self.rdp and not self.bootscan:
            return

        # Setup MaxChi
        if self.maxchi:
            print("Starting MaxChi Analysis")
            if config:
                maxchi = MaxChi(alignment, self.seq_names, settings=config['MaxChi'])
            else:
                maxchi = MaxChi(alignment, self.seq_names)

        # Setup Chimaera
        if self.chimaera:
            print("Starting Chimaera Analysis")
            if config:
                chimaera = Chimaera(alignment, self.seq_names, settings=config['Chimaera'])
            else:
                chimaera = Chimaera(alignment, self.seq_names)

        # Setup Siscan
        if self.siscan:
            print("Starting Siscan Analysis")
            if config:
                siscan = Siscan(alignment, self.seq_names, settings=config['SisScan'])
            else:
                siscan = Siscan(alignment, self.seq_names)

        # Setup RDP
        if self.rdp:
            print("Starting RDP Analysis")
            if config:
                rdp = RdpMethod(alignment, self.seq_names, settings=config['RDP'])
            else:
                rdp = RdpMethod(alignment, self.seq_names)

        # Setup Bootscan
        if self.bootscan:
            if config:
                bootscan = Bootscan(alignment, self.seq_names, settings=config['Bootscan'])
            else:
                bootscan = Bootscan(alignment, self.seq_names)

        # Run MaxChi
        if self.maxchi:
            maxchi.execute()

            # Run Chimaera
        if self.chimaera:
            chimaera.execute()

            # Run Siscan
        if self.siscan:
            siscan.execute(alignment)

            # Run RDP Method
        if self.rdp:
            rdp.execute(alignment)

            # Run Bootscan
        if self.bootscan:
            bootscan.execute(alignment)

        with open('results.txt', 'a') as outfile:
            # Report the results
            if self.maxchi:
                print("Finished MaxChi Analysis")
                outfile.write('MaxChi\n')
                for res in maxchi.results:
                    outfile.write(json.dumps(res))
                    outfile.write('\n')

                print(maxchi.results)

            if self.chimaera:
                print("Finished Chimaera Analysis")
                outfile.write('Chimaera\n')
                for res in chimaera.results:
                    outfile.write(json.dumps(res))
                    outfile.write('\n')

                print(chimaera.results)

            if self.siscan:
                print("Finished Siscan Analysis")
                outfile.write('Siscan\n')
                for res in siscan.results:
                    outfile.write(json.dumps(res))
                    outfile.write('\n')
                print(siscan.results)

            if self.rdp:
                print("Finished RDP Analysis")
                outfile.write('RDP\n')
                for res in rdp.results:
                    outfile.write(json.dumps(res))
                    outfile.write('\n')
                print(rdp.results)

            if self.bootscan:
                print("Finished Bootscan Analysis")
                outfile.write('Bootscan\n')
                for res in bootscan.results:
                    outfile.write(json.dumps(res))
                    outfile.write('\n')
                print(bootscan.results)
