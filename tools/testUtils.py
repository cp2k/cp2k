#! /usr/bin/env python
# prepares cp2k for checkin

import sys, re, os, os.path, commands, time, shutil, unittest, glob
from sys import argv
from os.path import join, basename
import prettify
import instantiateTemplates
import diffEpsilon

class simpleRunTest(unittest.TestCase):
    def __init__(self,cp2kRoot,testName,testsRoot,command="cp2k.sopt"):
        unittest.TestCase.__init__(self)
        self.cp2kRoot=cp2kRoot
        self.testName=testName
        self.command=command
        self.testRoot=join(cp2kRoot,testsRoot,testName+"-"+basename(command))
        self.inputFile=testName+".inp"
        self.outputFile=testName+".out"
    def setUp(self):
        if not os.access(os.path.normpath(join(self.testRoot,"..")),os.W_OK):
            os.mkdir(os.path.normpath(join(self.testRoot,"..")))
        os.mkdir(self.testRoot)
        os.symlink(join(self.cp2kRoot,"tests","QS","BASIS_SET"),
                   join(self.testRoot,"BASIS_SET"))
        os.symlink(join(self.cp2kRoot,"tests","QS","POTENTIAL"),
                   join(self.testRoot,"POTENTIAL"))
        os.symlink(join(self.cp2kRoot,"tests","QS","OT_BASIS"),
                   join(self.testRoot,"OT_BASIS"))
        shutil.copy(join(self.cp2kRoot,"tests","QS",self.inputFile),
                    self.testRoot)
    def tearDown(self):
        os.remove(join(self.testRoot,"BASIS_SET"))
        os.remove(join(self.testRoot,"POTENTIAL"))
        os.remove(join(self.testRoot,"OT_BASIS"))
        for f_path in glob.glob(join(self.testRoot,"RESTART*")):
            os.remove(f_path)
    def runTest(self):
        os.chdir(self.testRoot)
        pipe=os.popen("{ { "+self.command+" "+self.inputFile+
                      "; } 2>&1 ; } >>"+self.outputFile)
        if (pipe.close()):
            self.fail('error, the command returned an error')
        else:
            file1=open(join(self.cp2kRoot,"tests","QS","test_outputs",
                            commands.getoutput(join(self.cp2kRoot,"tools",
                                                    "get_arch_code"))
                            ,self.outputFile))
            file2=open(join(self.testRoot),self.outputFile)
            diffs=open(join(self.testRoot),self.testName+".diff","w")
            if diffEpsilon.compareCp2kOutput(file1,file2,logFile=diffs)!=0:
                self.fail('output was different than expected, see '+
                          diffs.name)
    def id(self):
        return "%s.%s-%s" % (self.__class__, self.testName,
                             basename(self.command))
    def shortDescription(self):
        return "interactive run of the test %s with %s" % (self.testName,
                                                           basename(self.command))

def simpleTests(cp2kRoot,testsRoot=None,command="cp2k.sopt"):
    if not testsRoot: testsRoot="test-"+time.strftime("%y%m%d-%H:%M")
    suite=unittest.TestSuite()
    suite.addTest(simpleRunTest(cp2kRoot=cp2kRoot,testName="H2O",
                                testsRoot=testsRoot,command=command))
    suite.addTest(simpleRunTest(cp2kRoot=cp2kRoot,testName="H2O-force",
                                testsRoot=testsRoot,command=command))
    return suite

if __name__=="__main__":
    cp2kRoot=os.path.normpath(join(os.getcwd(),
                                   os.path.dirname(sys.argv[0]),".."))
    runner=unittest.TextTestRunner()
    suite=simpleTests(cp2kRoot)
    runner.run(suite)
