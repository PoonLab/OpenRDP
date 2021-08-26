# OpenRDP

OpenRDP is an open-source re-implementation of the [RDP](http://web.cbio.uct.ac.za/~darren/rdp.html) recombination detection program.
RDP is a widely popular program for detecting recombination in genetic sequences - however, it has only been available as Windows binaries.
The developers of RDP, specifically [Darren Martin](http://www.idm.uct.ac.za/Darren_Martin), kindly provided the source code for the most recent version ([RDP5](https://academic.oup.com/ve/article/7/1/veaa087/6020281)) and permission to modify and redistribute this code.
We have re-implemented most of this code in the [Python](https://python.org) programming language, with the exception of some of the standalone third-party binaries ([3Seq](https://mol.ax/software/3seq/) and [GENECONV](https://www.math.wustl.edu/~sawyer/geneconv/index.html)).

This program runs a suite of recombination detection methods to filter recombinant nucleotide sequences from a multiple sequence alignment; or to detect recombination breakpoints in a query sequence by comparing it against a reference set.


## Installation

1. Clone the OpenRDP repository:
```console
git clone https://github.com/PoonLab/OpenRDP
```
If you do not have `git` installed, you can [download a release]().

2. Enter the repository directory `OpenRDP` and change directories to the `bin` subfolder:
```console
cd utils/bin
```

Since GENECONV is unlicensed, we require users to download the respective binary executable from [their site](https://www.math.wustl.edu/~sawyer/geneconv/index.html) (see "GENECONV PROGRAM FILES").

For Linux and macOS, you need to download and compile the source code:
```console
wget https://www.math.wustl.edu/~sawyer/geneconv/unix.source.tar.gz --no-check-certificate
tar -zvxf unix.source.tar.gz 
cd unix.source
csh makgconv.csh
```

[3Seq](https://mol.ax/software/3seq/) is released under the Creative Commons license ([CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/)), which restricts commerical use.
