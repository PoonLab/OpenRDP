import os
import unittest
import openrdp
from mpi4py import MPI

basepath = os.path.dirname(os.path.abspath(__file__))
LONG_INFILE = os.path.join(basepath, 'long.fasta')
CONFIG = os.path.join(basepath, 'test_cfg.ini')

class TestMPI(unittest.TestCase):
    def setUp(self):
        try:
            self.comm = MPI.COMM_WORLD
            self.my_rank = self.comm.Get_rank()
            self.size = self.comm.Get_size()
        except ImportError as error:
            print(f"Error importing MPI, {error}")
            exit()

        self.LONG_EXPECTED = \
                {'geneconv': [('Test1', ['-', '-'], (56, 170), '0.02728'),
                              ('Test1', ['Test2', '-'], (203, 507), '0.00829'),
                              ('Test1', ['Test2', '-'], (539, 759), '0.15378'),
                              ('Test1', ['Test3', '-'], (151, 195), '0.00210'),
                              ('Test2', ['Test3', '-'], (1, 204), '0.00002'),
                              ('Test4', ['-', '-'], (151, 193), '0.02202')],
                 'bootscan': [('Test2', ('Test1 ', 'Test3'), 760, 765, 0.06513627245570731)],
                 'threeseq': [('Test2', ('Test3', 'Test4'), '181', '787', '5.294757e-06'),
                              ('Test3', ('Test1', 'Test2'), '202', '787', '5.982096e-10')],
                'maxchi': [],
                'siscan': [],
                'chimaera': [],
                 'rdp': [('Test1 ', ('Test2', 'Test3'), 6, 15, 0.04173127820248898),
                         ('Test1 ', ('Test3', 'Test4'), 6, 504, 9.289462019022438e-63),
                         ('Test2', ('Test3', 'Test4'), 36, 481, 2.0810851150451177e-49)]}


    def test_MPI(self):
        scanner = openrdp.Scanner(cfg=CONFIG)
        results = scanner.run_scans(LONG_INFILE)

        results = results.dict if results else None

        self.comm.Barrier()
        all_results = self.comm.gather(results, root=0)
        results_sorted = {}

        if self.my_rank == 0:
            for method, output in all_results[0].items():
                results_sorted.update({method : sorted(output)})

            self.assertEqual(self.LONG_EXPECTED, results_sorted)

if __name__ == '__main__':
    unittest.main()
