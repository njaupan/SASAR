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

- [Minimap2](https://github.com/lh3/minimap2)
- [Minigraph](https://github.com/lh3/minigraph)
- [SiLiX](http://lbbe.univ-lyon1.fr/-SiLiX-?lang=en)
- Python 3 
    - biopython
    - pandas
    - pybedtools 

## <a name="Usage"></a>Usage
`
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
  `
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
