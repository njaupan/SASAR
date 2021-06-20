# SASAR 

SASAR is a meta-assembly tool to reconcile different assemblies without a reference guide. 

## Table of Contents

- [Introduction](#intro)
- [Users' Guide](#uguide)
  - [Installation](#install)
  - [Usage](#Usage)
  - [Example](#example)
- [Limitations](#limit)

## <a name="intro"></a>Introduction

## <a name="install"></a>Installation

- [Minimap2 version 2.20](https://github.com/lh3/minimap2) 
   ```
   git clone https://github.com/lh3/minimap2
   cd minimap2 && make
   ```
- [gfatools version 0.5](https://github.com/lh3/gfatools) 
   ```
   git clone https://github.com/lh3/gfatools
   cd gfatools && make
   ```

- [Minigraph version 0.15](https://github.com/lh3/minigraph)
   ```
   git clone https://github.com/lh3/minigraph
   cd minigraph && make
   ```
- [Python 3 ](https://salishsea-meopar-docs.readthedocs.io/en/latest/work_env/python3_conda_environment.html)
    * biopython
    * pandas
    * pybedtools
   ```
   # Install eihter by conda or pip
   conda install -c conda-forge biopython
   conda install -c bioconda pybedtools
   conda install -c anaconda pandas
   ```
- [Silix version 1.2.11](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-116)
   * A tool for Ultra-fast Sequence Clustering from Similarity Networks 
   * Make sure Boost libraries are also installed (for Ubuntu issue the following commands, requires **root permissions**)
   ```bash
   sudo apt-get install libboost-dev
   sudo apt-get install libboost-program-options-dev
   ```
   * Install SiLiX, **requires root**
   ```
    wget ftp://pbil.univ-lyon1.fr/pub/logiciel/silix/silix-1.2.11.tar.gz
    tar zxvf silix-1.2.11.tar.gz
    cd silix-1.2.11
    ./configure
    make
    make check
    sudo make install
    ```

## <a name="Usage"></a>Usage
```
usage: python SASAR.py in_dir

Long read assembly reconciliation

Positional arguments:
  in_dir             input directory containing all assemblies

Settings:
  -t INT             number of CPU threads for whole genome alignment
  -i INT             minimum identity confidence score [95]
  -c INT             minimum coverage confidence score [95]
  -m INT             minium overlap length [50000]
  --repeat_size INT  repeat size [50000]

Output options:
  -o PATH            output directory for SASAR-assembly [./SASAR_output]

Other:
  -h, --help         show this help message and exit
  -v, --version      show program's version number and exit
```
## Dataset 
### Haploid
- [D. melanogaster ISO1](https://www.ncbi.nlm.nih.gov/sra/SRX3676783) (144Mb |ONT 66x)
- [A. thaliana Col-0](https://www.ebi.ac.uk/ena/browser/view/PRJEB34954) (130Mb |ONT 130x)
- [O. sativa Nipponbare IRGSP1](https://www.ebi.ac.uk/ena/browser/view/PRJEB34954) (380Mb |ONT 34x)
- O. sativa Nipponbare
- [S. pennellii](https://plabipd.de/portal/solanum-pennellii) (0.9Gb |ONT 160x)

### Diploid
- [F. ananassa](https://www.ncbi.nlm.nih.gov/sra/?term=SRR11606867) (0.8Gb	|PB HiFi 36x)

## Docs
## Citation
