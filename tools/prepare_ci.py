#! /usr/bin/env python
# prepares cp2k for checkin

import sys, re, os, os.path, commands, time, shutil, unittest
from sys import argv
from os.path import join
import instantiateTemplates
import diffEpsilon
import buildUtils
import testUtils

# read directives
directives={"normalize-use":1,"upcase-keywords":1,"clean":1,
            "replace":0,"synopsis":1,"prettify-cvs":1,"popt":0,"tests":0}
directiveRe=re.compile(r"--(no-)?(normalize-use|upcase-keywords|"+
                       r"replace|synopsis|prettify-cvs|clean|popt|tests)$")
descStr=("usage:"+sys.argv[0]+"""
  [--[no-]normalize-use] [--[no-]upcase-keywords] [--[no-]replace]
  [--help] [--[no-]synopsis] [--[no-]prettify-cvs] [--[no-]clean]
  [--[no-]tests]

  Prepares for checkin the source.
  defaults="""+str(directives)
  )       
if "--help" in sys.argv[1:]:
    print descStr
    sys.exit(0)
for directive in sys.argv[1:]:
    m=directiveRe.match(directive)
    if m:
        directives[m.groups()[1]]=not m.groups()[0]
    else:
        print " ** ERROR **\nUnknown argument",directive
        print descStr
        sys.exit(-1)

cp2kRoot=os.path.abspath(os.path.join(os.path.dirname(argv[0]),".."))
logDirPath=join(cp2kRoot,"test-"+
                commands.getoutput(join(cp2kRoot,"tools","get_arch_code"))+
                "-"+time.strftime("%y%m%d-%H:%M"))
os.mkdir(logDirPath)
mainLog=open(join(logDirPath,"main.log"),'w')
mainLog.write(" ******* "+time.strftime("%y%m%d-%H:%M")+
              " prepare check-in BEGAN *******\n")
mainLog.flush()
outDir= join(cp2kRoot,"src","outDir")

print "logDirectory: "+logDirPath
print "main log: "+mainLog.name

buildUtils.prettifyCp2k(cp2kRoot=cp2kRoot,buildType="sdbg",
                        logDirPath=logDirPath,mainLog=mainLog,
                        directives=directives)

# clean? compile sopt
mainLog.write("====== clean="+str(directives["clean"])+
              " compile cp2k sopt ======\n")
mainLog.flush()
mainLog.write("  compilation logFile in 'cp2kBuildSopt.log'\n")
if not buildUtils.buildCp2k(cp2kRoot,"sopt",
                            join(logDirPath,"cp2kBuildSopt.log"),
                            clean=directives["clean"]):
    mainLog.write("+++ ERROR, build FAILED! +++\n")
else:
    mainLog.write("+++ build SUCESSFULL! +++\n")
mainLog.flush()

# clean? compile popt
if directives["popt"]:
    mainLog.write("====== clean="+directives["clean"]+
                  " compile cp2k popt ======\n")
    mainLog.flush()
    mainLog.write("  compilation logFile in 'cp2kBuildPopt.log'\n")
    if not buildUtils.buildCp2k(cp2kRoot,"popt",
                                join(logDirPath,"cp2kBuildPopt.log"),
                                clean=directives["clean"]):
        mainLog.write("+++ ERROR, build FAILED! +++\n")
    else:
        mainLog.write("+++ build SUCESSFULL! +++\n")

# do tests
if directives["tests"]:
    mainLog.write("====== tests ======\n")
    suite=testUtils.simpleTests(cp2kRoot=cp2kRoot,
                                testsRoot=logDirPath)
    testRunner=unittest.TextTestRunner(stream=mainLog,verbosity=2)
    testRunner.run(suite)

mainLog.write(" ******* "+time.strftime("%y%m%d-%H:%M")+
              " prepare check-in FINISHED *******\n")
mainLog.close()


