#! /usr/bin/env python
# performs test runs of cp2k

import sys, re, os, os.path, commands, time, shutil, unittest, glob
from sys import argv
from os.path import join, basename
import prettify
import instantiateTemplates
import diffEpsilon

defaultCp2kRoot=os.path.expanduser("~/cp2k")
defaultTestsRoot=os.path.expanduser("~/cp2k/test-"+time.strftime("%y%m%d-%H:%M"))

class simpleRunTest(unittest.TestCase):
    def __init__(self,testName,cp2kRoot=defaultCp2kRoot,
                 testsRoot=defaultTestsRoot,command="cp2k.sopt"):
        unittest.TestCase.__init__(self)
        self.cp2kRoot=cp2kRoot
        self.testName=testName
        self.archName=commands.getoutput(join(self.cp2kRoot,"tools",
                                              "get_arch_code"))
        self.command=join(cp2kRoot,"exe",self.archName,command)
        self.testRoot=join(cp2kRoot,testsRoot,testName+"-"+basename(command))
        self.inputFile=testName+".inp"
        self.outputFile=testName+"-"+os.path.basename(command)+".out"
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
            file1=open(join(self.testRoot,self.outputFile))
            if os.access(join(self.cp2kRoot,"tests","QS","test_outputs",
                              self.archName,
                              self.outputFile),os.R_OK):
                diffs=open(join(self.testRoot,self.testName+"-"
                                +basename(self.command)+".spec.diff"),"w")
                file2=open(join(self.cp2kRoot,"tests","QS","test_outputs",
                                self.archName,
                                self.outputFile))
                if diffEpsilon.compareCp2kOutput(file1,file2,logFile=diffs)!=0:
                    arch_spec_test="failed"
                else:
                    arch_spec_test="succeded"
            else:
                arch_spec_test=0
            if os.access(join(self.cp2kRoot,"tests","QS","test_outputs",
                              "generic",self.outputFile),os.R_OK):
                diffs=open(join(self.testRoot,self.testName+"-"
                                +basename(self.command)+".gen.diff"),"w")
                file2=open(join(self.cp2kRoot,"tests","QS","test_outputs",
                                "generic",self.outputFile))
                diffVal=diffEpsilon.compareCp2kOutput(file1,file2,logFile=diffs)
            elif not arch_spec_test:
                self.fail('no generic output for this test')
            else:
                diffVal=0
                print 'no generic output for this test'
            if diffVal>1.0e-14:
                if arch_spec_test=="succeded":
                    print "** generic diff to be updated? Too different from validated output"
                    print self.id(),"diff=",str(diffVal),'see',diffs.name
                elif arch_spec_test:
                    self.fail('output was different than expected for this platform (and too different from generic output), see '+
                              diffs.name)
                else:
                    self.fail('output was too different from generic output ('
                              +str(diffVal)+'), see '+diffs.name)
            else:
                if arch_spec_test!="succeded":
                    self.fail('output was different than expected for this platform (but generic output within error", see '+
                              diffs.name)
                print self.id(),"generic_diff=",diffVal
    
    def id(self):
        return "%s.%s-%s" % (self.__class__, self.testName,
                             basename(self.command))
    def shortDescription(self):
        return "interactive run of the test %s with %s" % (self.testName,
                                                           basename(self.command))

def simpleTests(cp2kRoot=defaultCp2kRoot,testsRoot=defaultTestsRoot,
                command="cp2k.sopt"):
    suite=unittest.TestSuite()
    suite.addTest(simpleRunTest(cp2kRoot=cp2kRoot,testName="H2O",
                                testsRoot=testsRoot,command=command))
    suite.addTest(simpleRunTest(cp2kRoot=cp2kRoot,testName="H2O-force",
                                testsRoot=testsRoot,command=command))
    suite.addTest(simpleRunTest(cp2kRoot=cp2kRoot,testName="H2O-32",
                                testsRoot=testsRoot,command=command))
    return suite

simpleTestsSopt=lambda :simpleTests()
simpleTestsSdbg=lambda :simpleTests(command="cp2k.sdbg")

for test in simpleTestsSopt()._tests:
    globals()[test.testName.replace("-","_")]=lambda :test

if __name__=="__main__":

    defaultCp2kRoot=os.path.normpath(join(os.getcwd(),
                                          os.path.dirname(sys.argv[0]),".."))
    defaultTestsRoot="test-"+time.strftime("%y%m%d-%H:%M")
    print "testsRoot:",join(defaultCp2kRoot,defaultTestsRoot)
    os.mkdir(join(defaultCp2kRoot,defaultTestsRoot))

    unittest.TestProgram(defaultTest="simpleTests")
