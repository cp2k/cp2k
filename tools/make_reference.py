#! /usr/bin/env python
# performs test runs of cp2k

import sys, re, os, os.path, commands, time, shutil, unittest, glob
from sys import argv
from os.path import join, basename

def make_ref(filePaths,cp2kRootDir,generic=0):
    if generic:
        destDir=join(cp2kRootDir,"tests","QS","test_outputs","generic")
    else:
        archName=commands.getoutput(join(self.cp2kRoot,"tools",
                                         "get_arch_code"))
        destDir=join(cp2kRootDir,"tests","QS","test_outputs",archName)
    if not os.access(destDir,os.W_OK): os.makedir(destDir)
    for filePath in filePaths:
        shutil.copy(filePath,destDir)

if __name__=="__main__":
    usage=os.path.basename(sys.argv[0])+""" [--[no-]generic] [--cp2kRoot /path/to/cp2k/] outputFile1 [outputFile2 ...]

    makes the given files defaults output for the current architecture
    (default) or if --generic was given for the generic outputs.
    tries to get the cp2kRoot path from the path of the first file."""

    optionRe=re.compile(r"--(no-)(generic)|--cp2kRoot")
    start=1
    generic=0
    cp2kRoot=""
    while start<len(argv):
        m=optionRe.match(argv[start])
        if not m: break
        if (m.groups()[1]):
            start=start+1
            generic=m.groups()[0]
        else:
            start=start+2
            cp2kRoot=argv[start-1]
    if not cp2kRoot:
        cp2kRoot=os.path.normpath(os.path.abspath(
            join(os.path.dirname(argv[start]),"..","..")))
    
    make_ref(filePaths=argv[start:],cp2kRootDir=cp2kRoot,
             generic=generic)
        
