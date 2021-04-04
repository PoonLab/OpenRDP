from preprocessing import *
from do_scans import *
from itertools import combinations
import configparser


class Scanner:
    def __init__(self, aln, args):
        self.infile = args.infile
        self.aln = aln
        self.rdp = True if args.rdp else False
        self.geneconv = True if args.geneconv else False
        self.threeseq = True if args.threeseq else False
        self.maxchi = True if args.maxchi else False
        self.chimaera = True if args.chimaera else False
        self.siscan = True if args.siscan else False
        self.bootscan = True if args.bootscan else False

        if args.cfg:
            self.cfg_file = args.cfg
        else:
            self.cfg_file = None

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

        # Create sequence objects
        seqs_in_aln = []
        num = 0
        for a in self.aln:
            seq = Sequence(a, self.aln, num)
            seqs_in_aln.append(seq)
            num += 1

        # Create alignment object
        alignment = Alignment(seqs_in_aln)

        # Run 3Seq
        if self.threeseq:
            three_seq = ThreeSeq(self.infile)
            print("Staring 3Seq Analysis")
            ts_results = three_seq.execute()
            print("Finished 3Seq Analysis")
            print(ts_results)

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

        # Exit early if 3Seq and Geneconv are the only methods selected
        if not self.maxchi and not self.chimaera and not self.siscan and not self.rdp and not self.bootscan:
            return

        # Setup MaxChi
        if self.maxchi:
            if config:
                maxchi = MaxChi(settings=config['MaxChi'])
            else:
                maxchi = MaxChi()

        # Setup Chimaera
        if self.chimaera:
            if config:
                chimaera = Chimaera(settings=config['Chimaera'])
            else:
                chimaera = Chimaera()

        # Setup Siscan
        if self.siscan:
            if config:
                siscan = Siscan(settings=config['SisScan'])
            else:
                siscan = Siscan()

        # Setup Bootscan
        if self.bootscan:
            if config:
                bootscan = Bootscan(settings=config['Bootscan'])
            else:
                bootscan = Bootscan()

        # Setup RDP
        if self.rdp:
            if config:
                rdp = RdpMethod(settings=config['RDP'])
            else:
                rdp = RdpMethod()

        for triplet in list(combinations(alignment.sequences, 3)):
            # Run MaxChi
            if self.maxchi:
                print("Starting MaxChi Analysis")
                mc_results = maxchi.execute(self.aln, triplet)
                print("Finished MaxChi Analysis")
                print(mc_results)

            # Run Chimaera
            if self.chimaera:
                print("Starting Chimaera Analysis")
                cm_results = chimaera.execute(self.aln, triplet)
                print("Finished Chimaera Analysis")
                print(cm_results)

            # Run Siscan
            if self.siscan:
                print("Starting Siscan Analysis")
                ss_results = siscan.execute(self.aln, triplet)
                print("Finished Siscan Analysis")
                print(ss_results)

            # Run Bootscan
            if self.bootscan:
                print("Starting Bootscan Analysis")
                bs_results = bootscan.execute(self.aln, triplet)
                print("Finished Bootscan Analysis")
                print(bs_results)

            # Run RDP Method
            if self.rdp:
                print("Starting RDP Analysis")
                rdp_results = rdp.execute(self.aln, triplet)
                print("Finished RDP Analysis")
                print(rdp_results)
