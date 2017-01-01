#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import operator
import fileinput
from os import path
from optparse import OptionParser

def main():
    usage = "Write a new kernel parameter file as an unique merge of an old parameter file and a new one called parameters.txt as created by collect.py"
    parser = OptionParser(usage)
    parser.add_option("-p", "--params", metavar="filename.txt",
          default="parameters_P100.txt",
          help="Default: %default")

    (options, args) = parser.parse_args()
    assert(len(args) == 0)
    param_fn = options.params

    list = []
    f = fileinput.input(files=(param_fn, "parameters.txt"))
    for line in f:
       if line.find("Kernel") > -1:
          i1 = line.find("m=") + 2
          i2 = line[i1:].find(",")
          m = int(line[i1:i1+i2])
          i1 = line.find("n=") + 2
          i2 = line[i1:].find(",")
          n = int(line[i1:i1+i2])
          i1 = line.find("k=") + 2
          i2 = line[i1:].find(",")
          k = int(line[i1:i1+i2])
          list.append([m,n,k,line.rstrip("\n")])
    f.close()

    f = open("parameters.new","w")
    f.write("# *****************************************************************************\n")
    f.write("# * CP2K: A general program to perform molecular dynamics simulations         *\n")
    f.write("# * Copyright (C) 2000 - 2017  CP2K developers group                          *\n")
    f.write("# *****************************************************************************\n")
    f.write("\n[\n")

    sorted_list = sorted(list)
    i = 0
    while i < len(sorted_list):
       m1 = sorted_list[i][0]
       n1 = sorted_list[i][1]
       k1 = sorted_list[i][2]
       i1 = sorted_list[i][3].find("#") + 2
       i2 = sorted_list[i][3].find("GFlop") - 1
       pwin = float(sorted_list[i][3][i1:i2])
       iwin = i
       j = i
       while j < len(sorted_list):
          if j == len(sorted_list) - 1:
             break
          j += 1
          m = sorted_list[j][0]
          n = sorted_list[j][1]
          k = sorted_list[j][2]
          i1 = sorted_list[j][3].find("#") + 2
          i2 = sorted_list[j][3].find("GFlop") - 1
          p = float(sorted_list[j][3][i1:i2])
          if m == m1 and n == n1 and k == k1:
             if p > pwin:
                pwin = p
                iwin = j
             i = j
          else:
             break
       f.write(sorted_list[iwin][3]+"\n")
       i += 1

    f.write("]\n\n#EOF\n")
    f.close()

    print("Wrote parameters.new")

#===============================================================================
if(len(sys.argv)==2 and sys.argv[-1]=="--selftest"):
    pass #TODO implement selftest
else:
    main()

#EOF
