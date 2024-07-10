# OpenRDP

OpenRDP is an open-source re-implementation of the [RDP](http://web.cbio.uct.ac.za/~darren/rdp.html) recombination detection program.
RDP is a widely popular program for detecting recombination in genetic sequences - however, it has only been available as Windows binaries
(These binaries are not compatible with Windows 11).

The developers of RDP, specifically [Darren Martin](http://www.idm.uct.ac.za/Darren_Martin), kindly provided the source code for the most recent version ([RDP5](https://academic.oup.com/ve/article/7/1/veaa087/6020281)) and permission to modify and redistribute this code.
This program runs a suite of methods to detect recombination events by comparing a nucleotide sequence to other sequences in a multiple alignment.

We have re-implemented most of this code in the [Python](https://python.org) programming language, with the exception of some of the standalone third-party binaries ([3Seq](https://mol.ax/software/3seq/) and [GENECONV](https://www.math.wustl.edu/~sawyer/geneconv/index.html)).
[3Seq](https://mol.ax/software/3seq/) is released under the Creative Commons license ([CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/)), which restricts commerical use.
Similarly, [GENECONV](https://www.math.wustl.edu/~sawyer/geneconv/index.html) is licensed for academic use and distribution, but commercial use is restricted.
For your convenience, we package these binaries along with their respective licenses.


## Dependencies 
* Python 3 (tested on version [3.8.10](https://www.python.org/downloads/release/python-3810/))
* `numpy` version [1.17.4](https://numpy.org/devdocs/release/1.17.4-notes.html) or later
* `scipy` version [1.5.0](https://docs.scipy.org/doc/scipy/reference/release.1.5.0.html) or later
* `h5py` version [3.8.0](https://docs.h5py.org/en/stable/build.html) or later
* `mpi4py` optional for parallel processing

## Installation

1. Clone the OpenRDP repository:
```console
git clone https://github.com/PoonLab/OpenRDP
```
If you do not have `git` installed, you can [download a release](https://github.com/PoonLab/OpenRDP/releases).

2. Switch into the package directory (`cd OpenRDP`) and run the installation script using either:
```
python3 -m pip install --user . 
```
or if you want to make `openrdp` available system-wide and you have super-user privileges on the computer:
```
sudo python3 -m pip install .
```
to do a local (single user) install.

## Usage

### Command-line interface (CLI)

```console
art@Wernstrom OpenRDP % openrdp -h
usage: openrdp [-h] [-o OUTFILE] [-c CFG] [-m  [...]] [-q] infile

An open source implementation of RDP5.

positional arguments:
  infile                File containing sequence alignment (FASTA or CLUSTAL) format

options:
  -h, --help            show this help message and exit
  -o OUTFILE, --outfile OUTFILE
                        Path to write CSV output
  -c CFG, --cfg CFG     Path to file that contains parameters. Defaults to
                        default.ini.
  -m  [ ...], --methods  [ ...]
                        Space-delimited list of recombination methods to run.
                        Permitted values are: {geneconv bootscan maxchi siscan
                        chimaera threeseq rdp}. All methods are run if not specified.
  -q, --quiet           Hide progress messages
```


The simplest command will run all seven recombination detection methods under the default settings, and prints the results to the console:
```
openrdp <input FASTA>
```

You can specify a different [configuration file](https://docs.python.org/3/library/configparser.html#supported-ini-file-structure) using the `-c` flag.  For example:
```console
art@Wernstrom OpenRDP % openrdp -c tests/test_cfg.ini tests/test_neisseria.fasta
Loading configuration from openrdp/tests/test_cfg.ini
Starting 3Seq Analysis
Finished 3Seq Analysis
Starting GENECONV Analysis
Finished GENECONV Analysis
Setting up bootscan analysis...
Starting Scanning Phase of Bootscan/Recscan
Finished Scanning Phase of Bootscan/Recscan
Setting up maxchi analysis...
Setting up siscan analysis...
Setting up chimaera analysis...
Setting up rdp analysis...
Scanning triplet 1 / 4
Scanning triplet 0 / 4
Scanning triplet 2 / 4
Scanning triplet 3 / 4
Scanning triplet 1 / 4
Scanning triplet 2 / 4
Scanning triplet 3 / 4
Scanning triplet 4 / 4

Method          Start   End     Recombinant     Parent1 Parent2 Pvalue
------------------------------------------------------------------------
Geneconv        1       204     X64866          X64869  -       2.00E-05
Geneconv        151     195     X64860          X64869  -       2.10E-03
Geneconv        203     507     X64860          X64866  -       8.29E-03
Geneconv        539     759     X64860          X64866  -       1.54E-01
Geneconv        151     193     X64873          -       -       2.20E-02
Geneconv        56      170     X64860          -       -       2.73E-02
Bootscan        760     765     X64873          X64860  X64866  6.51E-02
MaxChi          439     482     X64860          X64866  X64869  4.04E-02
MaxChi          475     518     X64860          X64866  X64873  4.04E-02
MaxChi          475     518     X64860          X64869  X64873  4.04E-02
Siscan          2       45      X64860          X64866  X64869  7.64E-01
Siscan          2       45      X64860          X64866  X64873  7.49E-01
Siscan          2       45      X64860          X64869  X64873  7.64E-01
Siscan          2       45      X64866          X64869  X64873  7.65E-01
Chimaera        198     241     X64860          X64866  X64869  2.05E-02
Chimaera        179     265     X64866          X64869  X64873  4.70E-03
Chimaera        170     213     X64873          X64860  X64869  1.81E-03
3Seq            202     787     X64869          X64860  X64866  5.98E-10
3Seq            181     787     X64866          X64869  X64873  5.29E-06
RDP             6       504     X64860          X64866  X64869  6.87E-07
RDP             6       496     X64860          X64866  X64873  3.05E-02
RDP             6       479     X64860          X64869  X64873  3.40E-02
RDP             6       18      X64866          X64869  X64873  2.36E+01
```

You can also provide a parental reference sequence using the `-r` flag. For example:
```console
(venv) sandeep@Jesry:~/OpenRDP$ openrdp -c tests/test_cfg.ini tests/short.fasta -r tests/short_ref.fasta
Loading configuration from tests/test_cfg.ini
Starting 3Seq Analysis
Finished 3Seq Analysis
Starting GENECONV Analysis
Finished GENECONV Analysis
Setting up bootscan analysis...
Starting Scanning Phase of Bootscan/Recscan
Finished Scanning Phase of Bootscan/Recscan
Setting up maxchi analysis...
Invalid option for 'win_size'.
Using default value (200) instead.
Setting up siscan analysis...
Invalid option for 'win_size'.
Using default value (200) instead.
Setting up chimaera analysis...
Invalid option for 'win_size'.
Using default value (200) instead.
Setting up rdp analysis...
Scanning triplet 1 / 15
Scanning triplet 3 / 15
Scanning triplet 2 / 15
Scanning triplet 4 / 15
Scanning triplet 5 / 15
Scanning triplet 7 / 15
Scanning triplet 6 / 15
Scanning triplet 8 / 15
Scanning triplet 9 / 15
Scanning triplet 10 / 15
Scanning triplet 11 / 15
Scanning triplet 13 / 15
Scanning triplet 12 / 15
Scanning triplet 14 / 15
Scanning triplet 15 / 15
Scanning triplet 1 / 15
Scanning triplet 2 / 15
Scanning triplet 3 / 15
Scanning triplet 4 / 15
Scanning triplet 5 / 15
Scanning triplet 6 / 15
Scanning triplet 7 / 15
Scanning triplet 8 / 15
Scanning triplet 9 / 15
Scanning triplet 10 / 15
Scanning triplet 11 / 15
Scanning triplet 12 / 15
Scanning triplet 13 / 15
Scanning triplet 14 / 15
Scanning triplet 15 / 15

Method          Start   End     Recombinant     Parent1 Parent2 Pvalue
------------------------------------------------------------------------
Siscan          2       205     A               J       R       7.84E-01
Siscan          2       205     A               K       R       7.56E-01
Siscan          2       205     A               J       K       7.77E-01
Siscan          2       205     B               J       R       7.77E-01
Siscan          3       205     B               K       R       7.34E-01
Siscan          2       205     B               J       K       7.57E-01
Siscan          2       205     C               J       R       7.77E-01
Siscan          2       205     C               K       R       7.08E-01
Siscan          2       205     C               J       K       7.62E-01
Siscan          2       205     D               J       R       7.84E-01
Siscan          2       205     D               K       R       7.58E-01
Siscan          3       205     D               J       K       7.57E-01
Siscan          2       205     E               J       R       7.67E-01
Siscan          2       205     E               K       R       7.08E-01
Siscan          2       205     E               J       K       7.60E-01
```
If you wish to run openrdp with `mpi4py`, you can do so with the following commands, where `-n` is the number of cores. Given the excessive text output, messages are mute. If you wish to display them anyways, use `-v`
```console
(myenv) [will@Paphlagon OpenRDP]$ mpirun -n 4 openrdp tests/long.fasta

Method          Start   End     Recombinant     Parent1 Parent2 Pvalue
------------------------------------------------------------------------
Geneconv        1       204     Test2           Test3   -       2.00E-05
Geneconv        151     195     Test1           Test3   -       2.10E-03
Geneconv        203     507     Test1           Test2   -       8.29E-03
Geneconv        539     759     Test1           Test2   -       1.54E-01
Geneconv        151     193     Test4           -       -       2.20E-02
Geneconv        56      170     Test1           -       -       2.73E-02
Bootscan        100     160     Test2           Test3   Test4   2.29E-02
3Seq            181     787     Test2           Test3   Test4   5.29E-06
3Seq            202     787     Test3           Test1   Test2   5.98E-10
RDP             6       433     Test1           Test3   Test4   9.28E-17
RDP             6       412     Test2           Test3   Test4   1.04E-15
```
### Using OpenRDP as a Python module 

OpenRDP can be used as a Python module.  In a typical workflow, we would start by creating an instance of the object class `Scanner`:
```python
>>> from openrdp import Scanner
>>> scanner = Scanner(cfg="tests/test_cfg.ini")
```
Calling an instance of Scanner with no arguments, *i.e.*, `Scanner()` loads the default package configuration file, `default.ini`.
We can also extract the configuration and modify specific settings:
```python
>>> cfg = scanner.get_config()
>>> cfg
{'Geneconv': {'indels_as_polymorphisms': False, 'mismatch_penalty': 1, 'min_len': 1, 'min_poly': 2, 'min_score': 2, 'max_num': 1}, 'MaxChi': {'max_pvalue': '0.05', 'win_size': 40, 'strip_gaps': False, 'fixed_win_size': True, 'num_var_sites': 70, 'frac_var_sites': '0.1'}, 'Chimaera': {'max_pvalue': '0.05', 'win_size': 40, 'strip_gaps': False, 'fixed_win_size': True, 'num_var_sites': 70, 'frac_var_sites': '0.1'}, 'RDP': {'max_pvalue': '0.05', 'reference_sequence': 'None', 'window_size': 40, 'min_identity': 0, 'max_identity': 100}, 'Bootscan': {'max_pvalue': '0.1', 'win_size': 20, 'step_size': 5, 'num_replicates': 100, 'random_seed': 3, 'cutoff_percentage': '0.7', 'scan': 'distances', 'np': 2}, 'Siscan': {'max_pvalue': '0.8', 'win_size': 40, 'step_size': 5, 'strip_gaps': True, 'pvalue_perm_num': 1100, 'scan_perm_num': 100, 'random_seed': 3}}
>>> cfg["Siscan"]["win_size"] = 50
>>> scanner.set_config(cfg)
```

The Scanner object has a `run_scans` method that takes the input file path as its only argument:
```python
>>> results = scanner.run_scans("tests/test_neisseria.fasta")
```

Scanner returns an instance of the object class ScanResults, which has a custom `__str__` attribute:
```python
>>> print(results)
Method  	Start	End	Recombinant	Parent1	Parent2	Pvalue
------------------------------------------------------------------------
Geneconv	1	204	X64866     	X64869 	-      	2.00E-05
Geneconv	151	195	X64860     	X64869 	-      	2.10E-03
Geneconv	203	507	X64860     	X64866 	-      	8.29E-03
Geneconv	539	759	X64860     	X64866 	-      	1.54E-01
Geneconv	151	193	X64873     	-      	-      	2.20E-02
Geneconv	56	170	X64860     	-      	-      	2.73E-02
Bootscan	760	765	X64866     	X64860 	X64869 	6.51E-02
MaxChi  	475	518	X64860     	X64866 	X64869 	4.04E-02
MaxChi  	475	518	X64860     	X64866 	X64873 	4.04E-02
MaxChi  	439	482	X64860     	X64869 	X64873 	4.04E-02
Siscan  	2	55	X64860     	X64866 	X64869 	7.52E-01
Siscan  	2	55	X64860     	X64866 	X64873 	7.67E-01
Siscan  	2	55	X64860     	X64869 	X64873 	7.65E-01
Siscan  	2	55	X64866     	X64869 	X64873 	7.65E-01
Chimaera	198	241	X64860     	X64869 	X64873 	2.05E-02
Chimaera	170	213	X64866     	X64860 	X64873 	1.81E-03
Chimaera	178	221	X64866     	X64869 	X64873 	1.17E-02
3Seq    	202	787	X64869     	X64860 	X64866 	5.98E-10
3Seq    	181	787	X64866     	X64869 	X64873 	5.29E-06
RDP     	6	15	X64860     	X64866 	X64869 	3.11E+01
RDP     	6	504	X64860     	X64869 	X64873 	5.51E-07
RDP     	36	481	X64866     	X64869 	X64873 	1.89E-05
```

This object also has a method to write CSV-formatted results to a file:
```python
>>> with open("results.csv", 'w') as outfile:
...     results.write(outfile)
...
```

A `ScanResults` object also behaves like a dictionary for accessing specific results, which are stored as a list of `dict` objects for each detection method:
```python
>>> list(results.keys())
['geneconv', 'bootscan', 'maxchi', 'siscan', 'chimaera', 'threeseq', 'rdp']
>>> results["geneconv"]
[{'start': 1, 'end': 204, 'recombinant': 'X64866', 'parent1': 'X64869', 'parent2': '-', 'pvalue': 2e-05}, {'start': 151, 'end': 195, 'recombinant': 'X64860', 'parent1': 'X64869', 'parent2': '-', 'pvalue': 0.0021}, {'start': 203, 'end': 507, 'recombinant': 'X64860', 'parent1': 'X64866', 'parent2': '-', 'pvalue': 0.00829}, {'start': 539, 'end': 759, 'recombinant': 'X64860', 'parent1': 'X64866', 'parent2': '-', 'pvalue': 0.15378}, {'start': 151, 'end': 193, 'recombinant': 'X64873', 'parent1': '-', 'parent2': '-', 'pvalue': 0.02202}, {'start': 56, 'end': 170, 'recombinant': 'X64860', 'parent1': '-', 'parent2': '-', 'pvalue': 0.02728}]
```
The actual dictionary can be accessed from `ScanResults.dict`.

### Docker

If you prefer to use Docker to run OpenRDP, a Dockerfile has been included in this repository.

To build the Docker image, navigate to the location of the Dockerfile and run the following command:

`docker build . -t openrdp`

The following command can then be used to run all seven recombination detection methods under the default settings:

`cat <input FASTA> | docker run -i --rm openrdp > <output CSV file>`
