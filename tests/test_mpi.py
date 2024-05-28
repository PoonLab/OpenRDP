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

        # expected results set because mpi shuffles the order
        self.LONG_EXPECTED = {
            "('Test2', ['Test3', '-'], (1, 204), '0.00002')", # genecov
                         "('Test1', ['Test3', '-'], (151, 195), '0.00210')",
                         "('Test1', ['Test2', '-'], (203, 507), '0.00829')",
                         "('Test1', ['Test2', '-'], (539, 759), '0.15378')",
                         "('Test4', ['-', '-'], (151, 193), '0.02202')",
                         "('Test1', ['-', '-'], (56, 170), '0.02728')",
            "('Test2', ('Test1 ', 'Test3'), 760, 765, 0.06513627245570731)", # bootscan
            "('Test1 ', ('Test2', 'Test3'), 475, 518, 0.04042768199451279)", # maxchi
                       "('Test1 ', ('Test2', 'Test4'), 475, 518, 0.04042768199451279)",
                       "('Test1 ', ('Test3', 'Test4'), 439, 482, 0.04042768199451279)",
            "('Test1 ', ('Test2', 'Test3'), 2, 45, 0.7490314187294909)", # siscan
                       "('Test1 ', ('Test2', 'Test4'), 2, 45, 0.7624489264586414)",
                       "('Test1 ', ('Test3', 'Test4'), 2, 45, 0.7639703544536595)",
                       "('Test2', ('Test3', 'Test4'), 2, 45, 0.7639703544536595)",
            "('Test1 ', ('Test3', 'Test4'), 198, 241, 0.02047438504938101)", # chimaera
                         "('Test2', ('Test1 ', 'Test4'), 170, 213, 0.0018132288986577026)",
                         "('Test2', ('Test3', 'Test4'), 176, 219, 0.004701217146256585)",
            "('Test3', ('Test1', 'Test2'), '202', '787', '5.982096e-10')", # threeseq
                         "('Test2', ('Test3', 'Test4'), '181', '787', '5.294757e-06')",
            "('Test1 ', ('Test2', 'Test3'), 6, 15, 31.159735288655444)", # rdp
                    "('Test1 ', ('Test3', 'Test4'), 6, 504, 6.867230672435156e-07)",
                    "('Test2', ('Test3', 'Test4'), 36, 481, 1.666567338773972e-05)"
        }

    def test_MPI(self):
        scanner = openrdp.Scanner(cfg=CONFIG)
        results = scanner.run_scans(LONG_INFILE).dict

        self.comm.Barrier()
        all_results = self.comm.gather(results, root=0)

        if self.my_rank == 0:
            overall = set()
            for result in all_results:
                for method, output in result.items():
                    for line in output:
                        overall.add(str(line))
            self.assertEqual(self.LONG_EXPECTED, overall)

if __name__ == '__main__':
    unittest.main()
