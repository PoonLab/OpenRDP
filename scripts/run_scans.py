from preprocessing import *
from do_scans import *
from itertools import combinations


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

    def run_scans(self):
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
            geneconv = GeneConv()
            geneconv.validate_options()
            print("Starting GENECONV Analysis")
            gc_results = geneconv.execute(self.aln)
            print("Finished GENECONV Analysis")
            print(gc_results)

        for triplet in list(combinations(alignment.sequences, 3)):

            # Run MaxChi
            if self.maxchi:
                maxchi = MaxChi()
                maxchi.validate_options()
                print("Starting MaxChi Analysis")
                mc_results = maxchi.execute(self.aln, triplet)
                print("Finished MaxChi Analysis")
                print(mc_results)

            # Run Chimaera
            if self.chimaera:
                chimarea = Chimaera()
                chimarea.validate_options()
                print("Starting Chimaera Analysis")
                cm_results = chimarea.execute(self.aln, triplet)
                print("Finished Chimaera Analysis")
                print(cm_results)

            # Run Siscan
            if self.siscan:
                siscan = Siscan()
                siscan.validate_options()
                print("Starting Siscan Analysis")
                ss_results = siscan.execute(self.aln, triplet)
                print("Finished Siscan Analysis")
                print(ss_results)

            if self.bootscan:
                bootscan = Bootscan()
                bootscan.validate_options()
                print("Starting Bootscan Analysis")
                bs_results = bootscan.execute(self.aln, triplet)
                print("Finished Bootscan Analysis")
                print(bs_results)
