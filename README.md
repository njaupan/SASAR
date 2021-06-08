# SASAR
Packages that required

SiLiX

Download version 1.2.11
Make sure Boost libraries are also installed (for Ubuntu issue the following commands, requires root permissions)
sudo apt-get install libboost-dev
sudo apt-get install libboost-program-options-dev

Install SiLiX, requires root

wget ftp://pbil.univ-lyon1.fr/pub/logiciel/silix/silix-1.2.11.tar.gz
tar zxvf silix-1.2.11.tar.gz
cd silix-1.2.11
./configure
make
make check
sudo make install
