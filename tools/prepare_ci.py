#! /usr/bin/env python
# prepares cp2k for checkin

import sys, re, os, os.path, commands, time, shutil
from sys import argv
from os.path import join
import instantiateTemplates
import diffEpsilon
import buildUtils

# read directives
directives={"normalize-use":1,"upcase-keywords":1,"clean":1,
            "replace":0,"synopsis":1,"prettify-cvs":1,"popt":0,"tests":1}
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
    mainLog.write("====== H2O test ======\n")
    mainLog.flush()
    exePath= join(cp2kRoot,"exe",
                  commands.getoutput(join(cp2kRoot,"tools","get_arch_code")),
                  "cp2k.sopt")
    log1Path=join(logDirPath,"H2O.out")
    os.chdir(join(cp2kRoot,"tests","QS"))
    commands.getoutput("{ { "+exePath+" "+
                       join(cp2kRoot,"tests","QS","H2O.inp")+
                       " ; } 2>&1 ; } > "+log1Path)
    logFile=open(os.path.join(logDirPath,"H2O.diffs"),'w')
    file1=open(join(cp2kRoot,"tests","QS","H2O.out"),'r')
    file2=open(log1Path)
    diffVal=diffEpsilon.compareCp2kOutput(file1,file2,
                                          0.0001,logFile)
    file1.close(); file2.close()
    logFile.write("totalDiff="+`diffVal`+"\n")
    if diffVal>0.0001:
        mainLog.write("+++ ERROR, H2O test failed +++\n diff="+`diffVal`+
                      " more info in H2O.diffs\n")
    else:
        mainLog.write("+++ H2O test SUCESSFULL (diff="+`diffVal`+")! +++\n")
    logFile.close()

mainLog.write(" ******* prepare check-in FINISHED *******\n")
mainLog.close()


