#! /usr/bin/env python

import sys
import re
import string
from sys import argv

def convert_into_dp(code):
    c1=code
    floatRe=re.compile(r"([0-9]+\.?[0-9]*|\.[0-9]+)[dD]([-+]?[0-9]+)")
    c1=floatRe.sub(r"\1e\2_dp",c1)
    return c1

if __name__ == '__main__':
    import os.path
    if len(sys.argv)<2:
        print "usage:", sys.argv[0]," file1 [file2 ...] "
    else:
        for fileName in sys.argv[1:]:
            try:
               tmpName="TmpFile"
               infile=open(fileName,'r')
               outfile=open(tmpName,'w')
               while 1:
                   line=infile.readline().replace("\t",8*" ")
                   if not line: break
                   outfile.write(convert_into_dp(line))
               infile.close()
               outfile.close()
               os.rename(tmpName, fileName)
            except:
               print "error for file", fileName

