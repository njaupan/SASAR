# SASAR 

SASAR is a meta-assembly tool to reconcile different assemblies without a reference guide. 

## Getting Started
## Install required dependencies

- [Minimap2](https://github.com/lh3/minimap2)
- [Minigraph](https://github.com/lh3/minigraph)
- [SiLiX](http://lbbe.univ-lyon1.fr/-SiLiX-?lang=en)

   * Download version **1.2.11**
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
- Python 3 (with the following auto-installed packages)
    - biopython
    - pandas
    - pybedtools 
```bash
# install with conda
conda create --name SASAR python=3.8
conda install -n SASAR -c anaconda -y biopython
conda install -n SASAR -c bioconda -y minimap2 pybedtools pandas

```
## Dataset 
### Haploid
- [D. melanogaster ISO1](https://www.ncbi.nlm.nih.gov/sra/SRX3676783) (144Mb |ONT 66x)
- [A. thaliana Col-0](https://www.ebi.ac.uk/ena/browser/view/PRJEB34954) (130Mb |ONT 130x)
- [O. sativa Nipponbare IRGSP1](https://www.ebi.ac.uk/ena/browser/view/PRJEB34954) (380Mb |ONT 34x)
- [O. sativa Nipponbare]
- [S. pennellii](https://plabipd.de/portal/solanum-pennellii) (0.9Gb |ONT 160x)

### Diploid
- [F. ananassa](https://www.ncbi.nlm.nih.gov/sra/?term=SRR11606867) (0.8Gb	|PB HiFi 36x)

## Docs
## Citation
