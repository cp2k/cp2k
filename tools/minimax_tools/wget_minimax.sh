#!/bin/bash -e

# See also https://doi.org/10.1007/s00791-018-00308-4

#wget https://www.mis.mpg.de/scicomp/EXP_SUM/1_x/1_xData.zip
wget https://www.cp2k.org/static/downloads/1_xData.zip

echo "7be2e56d83d0cb17683bbc8ab85dae1dc9a9c937e1dc1bad1514857938e687cb  1_xData.zip" | sha256sum --check
unzip -q -d 1_xData 1_xData.zip

#EOF
