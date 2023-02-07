# OpenRDP

OpenRDP is an open-source re-implementation of the [RDP](http://web.cbio.uct.ac.za/~darren/rdp.html) recombination detection program.
RDP is a widely popular program for detecting recombination in genetic sequences - however, it has only been available as Windows binaries.
The developers of RDP, specifically [Darren Martin](http://www.idm.uct.ac.za/Darren_Martin), kindly provided the source code for the most recent version ([RDP5](https://academic.oup.com/ve/article/7/1/veaa087/6020281)) and permission to modify and redistribute this code.
This program runs a suite of recombination detection methods to filter recombinant nucleotide sequences from a multiple sequence alignment; or to detect recombination breakpoints in a query sequence by comparing it against a reference set.

We have re-implemented most of this code in the [Python](https://python.org) programming language, with the exception of some of the standalone third-party binaries ([3Seq](https://mol.ax/software/3seq/) and [GENECONV](https://www.math.wustl.edu/~sawyer/geneconv/index.html)).
[3Seq](https://mol.ax/software/3seq/) is released under the Creative Commons license ([CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/)), which restricts commerical use.
Similarly, [GENECONV](https://www.math.wustl.edu/~sawyer/geneconv/index.html) is licensed for academic use and distribution, but commercial use is restricted.
For your convenience, we package these binaries along with their respective licenses.


## Dependencies 
* Python 3 (tested on version [3.8.10](https://www.python.org/downloads/release/python-3810/))
* `numpy` version [1.17.4](https://numpy.org/devdocs/release/1.17.4-notes.html) or later
* `scipy` version [1.5.0](https://docs.scipy.org/doc/scipy/reference/release.1.5.0.html) or later
* `h5py` version [3.8.0](https://docs.h5py.org/en/stable/build.html) or later


## Installation

1. Clone the OpenRDP repository:
```console
git clone https://github.com/PoonLab/OpenRDP
```
If you do not have `git` installed, you can [download a release](https://github.com/PoonLab/OpenRDP/releases).

2. Switch into the package directory (`cd OpenRDP`) and run the installation script using either:
```
sudo python3 setup.py install
```
if you have super-user privileges on the computer, or:
```
python3 setup.py install --user
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
Loading configuration from tests/test_cfg.ini
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


### Using OpenRDP as a Python module 

Using the following steps, OpenRDP can be used as a Python module 

```console
kaitlyn@boo-ubuntu:~/OpenRDP$ ipython
Python 3.8.10 (default, Jun  2 2021, 10:49:15) 
Type 'copyright', 'credits' or 'license' for more information
IPython 7.18.1 -- An enhanced Interactive Python. Type '?' for help.

In [1]: from openrdp.main import openrdp

In [2]: openrdp('./openrdp/tests/test_neisseria.fasta', './test.csv', './openrdp/tests/test_cfg.ini')
Starting 3Seq Analysis
Finished 3Seq Analysis
Starting GENECONV Analysis
Finished GENECONV Analysis
Setting Up MaxChi Analysis
Setting Up Chimaera Analysis
Setting Up RDP Analysis
Setting Up Bootscan Analysis
Starting Scanning Phase of Bootscan/Recscan
Finished Scanning Phase of Bootscan/Recscan
Setting Up Siscan Analysis
Scanning triplet 1 / 4
Scanning triplet 2 / 4
Scanning triplet 3 / 4
Scanning triplet 4 / 4

Method               StartLocation        EndLocation          Recombinant          Parent1              Parent2              Pvalue              
3Seq                 202                  787                  X64869               X64860               X64866               5.982096e-10        
3Seq                 181                  787                  X64866               X64869               X64873               5.294757e-06        
Geneconv             1                    204                  X64866               X64869               -                    0.00002             
Geneconv             151                  195                  X64860               X64869               -                    0.00210             
Geneconv             203                  507                  X64860               X64866               -                    0.00829             
Geneconv             539                  759                  X64860               X64866               -                    0.15378             
Geneconv             151                  193                  X64873               -                    -                    0.02202             
Geneconv             56                   170                  X64860               -                    -                    0.02728             
Bootscan             760                  765                  X64873               X64860               X64869               0.06513627245570731 
MaxChi               439                  482                  X64860               X64866               X64869               0.04042768199451279 
MaxChi               475                  518                  X64860               X64869               X64873               0.04042768199451279 
MaxChi               475                  518                  X64866               X64869               X64873               0.04042768199451279 
Chimaera             176                  219                  X64860               X64866               X64873               0.0019834358538684586
Chimaera             192                  235                  X64869               X64860               X64866               0.02047438504938101 
Chimaera             99                   142                  X64869               X64866               X64873               0.01174095674176845 
Chimaera             243                  286                  X64873               X64860               X64869               0.0018132288986577026
RDP                  6                    496                  X64860               X64869               X64873               0.000567470376377077
RDP                  6                    479                  X64866               X64869               X64873               0.0007659780246283629
Siscan               2                    45                   X64860               X64866               X64869               0.7593011319150192  
Siscan               2                    45                   X64860               X64866               X64873               0.7646778896282416  
Siscan               2                    45                   X64860               X64869               X64873               0.7663004734577327  
Siscan               2                    45                   X64866               X64869               X64873               0.7677102603374792  
```

