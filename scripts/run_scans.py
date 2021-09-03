import configparser

import numpy as np

from scripts.bootscan import Bootscan
from scripts.chimaera import Chimaera
from scripts.common import generate_triplets, Triplet
from scripts.geneconv import GeneConv
from scripts.maxchi import MaxChi
from scripts.rdp import RdpMethod
from scripts.siscan import Siscan
from scripts.threeseq import ThreeSeq


class Scanner:
    def __init__(self, names, infile, cfg, run_geneconv=False, run_three_seq=False, run_rdp=False,
                 run_siscan=False, run_maxchi=False, run_chimaera=False, run_bootscan=False, quiet=False):
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
        self.quiet = quiet

    def run_scans(self, aln):
        """
        Run the selected recombination detection analyses
        """
        # Parse config file
        if self.cfg_file:
            config = configparser.ConfigParser()
            config.read(self.cfg_file)
        else:
            config = None

        # Remove identical sequences
        aln = list(set(aln))

        # Create an m x n array of sequences
        alignment = np.array(list(map(list, aln)))

        threeseq_res, geneconv_res = None, None
        bootscan_res, maxchi_res, chimaera_res, rdp_res, siscan_res = None, None, None, None, None

        # Run 3Seq
        if self.threeseq:
            three_seq = ThreeSeq(self.infile)
            if not self.quiet:
                print("Staring 3Seq Analysis")

            threeseq_res = three_seq.execute()
            if not self.quiet:
                print("Finished 3Seq Analysis")

        # Run GENECONV
        if self.geneconv:
            # Parse output file if available
            if config:
                geneconv = GeneConv(settings=dict(config.items('Geneconv')))
            else:
                geneconv = GeneConv()

            if not self.quiet:
                print("Starting GENECONV Analysis")
            geneconv_res = geneconv.execute(self.infile)

            if not self.quiet:
                print("Finished GENECONV Analysis")

        # Exit early if 3Seq and Geneconv are the only methods selected
        if not self.maxchi and not self.chimaera and not self.siscan and not self.rdp and not self.bootscan:
            self.write_output(threeseq_res, geneconv_res, bootscan_res, maxchi_res, chimaera_res, rdp_res, siscan_res)
            return

        triplets = []
        for trp in generate_triplets(alignment):
            triplets.append(Triplet(alignment, self.seq_names, trp))

        maxchi, chimaera, rdp, bootscan, siscan = None, None, None, None, None

        # Setup MaxChi
        if self.maxchi:
            if not self.quiet:
                print("Starting MaxChi Analysis")
            if config:
                maxchi = MaxChi(alignment, settings=dict(config.items('MaxChi')))
            else:
                maxchi = MaxChi(alignment)

        # Setup Chimaera
        if self.chimaera:
            if not self.quiet:
                print("Starting Chimaera Analysis")
            if config:
                chimaera = Chimaera(alignment, settings=dict(config.items('Chimaera')))
            else:
                chimaera = Chimaera(alignment)

        # Setup RDP
        if self.rdp:
            if not self.quiet:
                print("Starting RDP Analysis")
            if config:
                rdp = RdpMethod(alignment, settings=dict(config.items('RDP')))
            else:
                rdp = RdpMethod(alignment)

        # Setup Bootscan
        if self.bootscan:
            if not self.quiet:
                print("Starting Bootscan Analysis")
            if config:
                bootscan = Bootscan(alignment, settings=dict(config.items('Bootscan')), quiet=self.quiet)
            else:
                bootscan = Bootscan(alignment, quiet=self.quiet)

        # Setup Siscan
        if self.siscan:
            if not self.quiet:
                print("Starting Siscan Analysis")
            if config:
                siscan = Siscan(alignment, settings=dict(config.items('Siscan')))
            else:
                siscan = Siscan(alignment)

        # Run MaxChi
        if self.maxchi:
            maxchi_res = maxchi.execute(triplets, self.quiet)
            if not self.quiet:
                print("Finished MaxChi Analysis")

        # Run Chimaera
        if self.chimaera:
            chimaera_res = chimaera.execute(triplets, self.quiet)
            if not self.quiet:
                print("Finished Chimaera Analysis")

        # Run RDP Method
        if self.rdp:
            rdp_res = rdp.execute(triplets, self.quiet)
            if not self.quiet:
                print("Finished RDP Analysis")

        # Run Bootscan
        if self.bootscan:
            bootscan_res = bootscan.execute(triplets, self.quiet)
            if not self.quiet:
                print("Finished Bootscan Analysis")

        # Run Siscan
        if self.siscan:
            siscan_res = siscan.execute(triplets, self.quiet)
            if not self.quiet:
                print("Finished Siscan Analysis")

        self.write_output(threeseq_res, geneconv_res, bootscan_res, maxchi_res, chimaera_res, rdp_res, siscan_res)

        if not self.quiet:
            self.print_output(threeseq_res, geneconv_res, bootscan_res, maxchi_res, chimaera_res, rdp_res, siscan_res)

    def write_output(self, threeseq_res, geneconv_res, bootscan_res, maxchi_res, chimaera_res, rdp_res, siscan_res):
        """
        Write results to a file
        :param threeseq_res: results of 3Seq analysis
        :param geneconv_res: results of Geneconv analysis
        :param bootscan_res: results of Bootscan analysis
        :param maxchi_res: results of MaxChi analysis
        :param chimaera_res: results of Chimaera analysis
        :param rdp_res: results of RDP analysis
        :param siscan_res: results of Siscan analysis
        """
        with open('results.csv', 'w+') as outfile:
            outfile.write('Method,StartLocation,EndLocation,Recombinant,Parent1,Parent2,Pvalue\n')

            if self.threeseq:
                for event in threeseq_res:
                    outfile.write('3Seq,{},{},{},{},{},{}\n'
                                  .format(event[2], event[3], event[0], event[1][0], event[1][1], event[4]))
            if self.geneconv:
                for event in geneconv_res:
                    outfile.write('Geneconv,{},{},{},{},{},{}\n'
                                  .format(event[2][0], event[2][1], event[0], event[1][0], event[1][1], event[3]))
            if self.bootscan:
                for event in bootscan_res:
                    outfile.write('Bootscan,{},{},{},{},{},{}\n'
                                  .format(event[2], event[3], event[0], event[1][0], event[1][1], event[4]))
            if self.maxchi:
                for event in maxchi_res:
                    outfile.write('MaxChi,{},{},{},{},{},{}\n'
                                  .format(event[2], event[3], event[0], event[1][0], event[1][1], event[4]))
            if self.chimaera:
                for event in chimaera_res:
                    outfile.write('Chimaera,{},{},{},{},{},{}\n'
                                  .format(event[2], event[3], event[0], event[1][0], event[1][1], event[4]))
            if self.rdp:
                for event in rdp_res:
                    outfile.write('RDP,{},{},{},{},{},{}\n'
                                  .format(event[2], event[3], event[0], event[1][0], event[1][1], event[4]))
            if self.siscan:
                for event in siscan_res:
                    outfile.write('Siscan,{},{},{},{},{},{}\n'
                                  .format(event[2], event[3], event[0], event[1][0], event[1][1], event[4]))

    def print_output(self, threeseq_res, geneconv_res, bootscan_res, maxchi_res, chimaera_res, rdp_res, siscan_res):
        """
        Print results to console
        :param threeseq_res: results of 3Seq analysis
        :param geneconv_res: results of Geneconv analysis
        :param bootscan_res: results of Bootscan analysis
        :param maxchi_res: results of MaxChi analysis
        :param chimaera_res: results of Chimaera analysis
        :param rdp_res: results of RDP analysis
        :param siscan_res: results of Siscan analysis
        """
        print('Method,StartLocation,EndLocation,Recombinant,Parent1,Parent2,Pvalue')
        if self.threeseq:
            for event in threeseq_res:
                print('3Seq\t{}\t{}\t{}\t{}\t{}\t{}'
                      .format(event[2], event[3], event[0], event[1][0], event[1][1], event[4]))
        if self.geneconv:
            for event in geneconv_res:
                print('Geneconv\t{}\t{}\t{}\t{}\t{}\t{}'
                      .format(event[2][0], event[2][1], event[0], event[1][0], event[1][1], event[3]))
        if self.bootscan:
            for event in bootscan_res:
                print('Bootscan\t{}\t{}\t{}\t{}\t{}\t{}'
                      .format(event[2], event[3], event[0], event[1][0], event[1][1], event[4]))
        if self.maxchi:
            for event in maxchi_res:
               print('MaxChi\t{}\t{}\t{}\t{}\t{}\t{}'
                     .format(event[2], event[3], event[0], event[1][0], event[1][1], event[4]))
        if self.chimaera:
            for event in chimaera_res:
                print('Chimaera\t{}\t{}\t{}\t{}\t{}\t{}'
                      .format(event[2], event[3], event[0], event[1][0], event[1][1], event[4]))
        if self.rdp:
            for event in rdp_res:
                print('RDP\t{}\t{}\t{}\t{}\t{}\t{}'
                      .format(event[2], event[3], event[0], event[1][0], event[1][1], event[4]))
        if self.siscan:
            for event in siscan_res:
                print('Siscan\t{}\t{}\t{}\t{}\t{}\t{}'
                      .format(event[2], event[3], event[0], event[1][0], event[1][1], event[4]))
