#! /usr/bin/env python
# prepares cp2k for checkin

import sys
from sys import argv
import os
from os.path import join
import commands
import addSynopsis # .addSynopsisInDir
import instantiateTemplates # .evaluateInstantiationFile
import time
import diffEpsilon

def buildCp2k(cp2kRoot,buildType="sopt",logFilePath=None,clean=None):
    if not logFilePath:
        logFilePath=join(cp2kRoot,"build"+time.strftime("%y%m%d-%H:%M")+
                "."+buildType+".log")
    buildDir= join(cp2kRoot,"obj",
	commands.getoutput(join(cp2kRoot,"tools","get_arch_code")),buildType)
    if clean and os.access(buildDir,os.W_OK):
        commands.getoutput('rm -rf "'+buildDir+'"')
    os.chdir(join(cp2kRoot,"makefiles"))
    if os.access("/usr/bin/gnumake",os.X_OK): makeCmd="/usr/bin/gnumake"
    else: makeCmd="gmake"
    pipe=os.popen("{ { "+makeCmd+" "+buildType+"; } 2>&1 ; } >>"+logFilePath)
    logFile=open(logFilePath,'a')
    if (pipe.close()):
        logFile.write("\n+++ ERROR, build "+buildType+" FAILED! +++\n")
        logFile.close()
        return None
    else:
        logFile.write("\n+++ build "+buildType+" SUCESSFULLY! +++\n")
        logFile.close()
        return 1


cp2kRoot=os.path.abspath(os.path.join(os.path.dirname(argv[0]),".."))
logDirPath=join(cp2kRoot,"test-"+
                commands.getoutput(join(cp2kRoot,"tools","get_arch_code"))+
                "-"+time.strftime("%y%m%d-%H:%M"))
os.mkdir(logDirPath)
mainLog=open(join(logDirPath,"main.log"),'w')
outDir= join(cp2kRoot,"src","outDir")

print "logDirectory: "+logDirPath
print "main log: "+mainLog.name

mainLog.write("===== cleaning forpar =====\n")
os.chdir(join(cp2kRoot,"tools"))
if os.access("/usr/bin/gnumake",os.X_OK): makeCmd="/usr/bin/gnumake"
else: makeCmd="gmake"
os.popen("{ { "+makeCmd+" clean ; } 2>&1 ; } >/dev/null").close()
os.popen("{ { "+makeCmd+" ; } 2>&1 ; } >/dev/null").close()

mainLog.write("===== instantiating templates =====\n")
mainLog.flush()
import glob
logFile=open(os.path.join(logDirPath,"templateInstantiation.log"),'w')
templateDir=join(logDirPath,"outTemplates")
os.mkdir(templateDir)
instantiationFiles=glob.glob(os.path.join(cp2kRoot,"src","*.instantiation"))
templateInstances=[]
for instantiationFile in instantiationFiles:
    templateInstances.extend(
        instantiateTemplates.evaluateInstantiationFile(instantiationFile,
                                                       logFile,templateDir))
mainLog.write(" template generation logFile in '%s'\n"%
              os.path.basename(logFile.name))
logFile.close()

# build sdbg
mainLog.write("====== dirty compiling cp2k sdbg ======\n")
mainLog.flush()
mainLog.write("  compilation logFile in 'cp2kBuildSdbg.log'\n")
if not buildCp2k(cp2kRoot,"sdbg",join(logDirPath,"cp2kBuildSdbg.log")):
    mainLog.write("+++ ERROR, build FAILED! +++\n")
else:
    mainLog.write("+++ build SUCESSFULL! +++\n")

# add synopsis
mainLog.write("====== updating synopsis ======\n")
mainLog.flush()
logFile=open(os.path.join(logDirPath,"addSynopsis.log"),'w')
buildDir= join(cp2kRoot,"obj",
    commands.getoutput(join(cp2kRoot,"tools","get_arch_code")),
    "sdbg")
outDir=join(logDirPath,"synopsisDir")
if os.access(outDir,os.W_OK): commands.getoutput('rm -rf "'+outDir+'"')
os.mkdir(outDir)

filesToSyn2=(glob.glob(os.path.join(cp2kRoot,"src","cp_*.F"))+
            glob.glob(os.path.join(cp2kRoot,"src","pao_*.F"))+
            glob.glob(os.path.join(templateDir,"*.F")))
filesToSyn=templateInstances
baseNames=map(os.path.basename,filesToSyn)
for fileToSyn in filesToSyn2:
    if not os.path.basename(fileToSyn) in baseNames:
	filesToSyn.append(fileToSyn)
addSynopsis.addSynopsisInDir(buildDir,outDir,filesToSyn,logFile)


# check synopsis diffs
os.chdir(outDir)
filesSyn=glob.glob("*.F")
for fileSyn in filesSyn:
    if not os.path.exists(join(cp2kRoot,"src",fileSyn)):
        os.rename(fileSyn,join(cp2kRoot,"src",fileSyn))
        diffs=''
    else:
        origFilePath=join(cp2kRoot,"src",fileSyn)
        diffs=commands.getoutput("diff "+fileSyn+" "+
                                 origFilePath)
        if diffs!='' and fileSyn in baseNames:
            origFilePath=templateInstances[baseNames.index(fileSyn)]
            diffs=commands.getoutput("diff "+fileSyn+" "+
                                     origFilePath)
    if diffs=='':
        diffOk=None
    else:
        diffOk=1
        import re
        commentOrHeader = re.compile("([0-9]+,*[0-9]*[adc][0-9]+,*[0-9]*$|[><] *|[><] *!| *$|---$)")
        for line in diffs.splitlines():
            if not commentOrHeader.match(line):
                mainLog.write("+ WARNING checking synopsis for "+fileSyn+
                              "  found the unaccetable line.\n  line='"+line+
                              "'  addSynopsis not performed\n")
                diffOk=None
                break
    if diffOk:
        os.rename(fileSyn,join(cp2kRoot,"src",os.path.basename(fileSyn)))
mainLog.write("  addSynopsis logFile in '%s'\n"%os.path.basename(logFile.name))
logFile.close()

# clean? compile
shouldClean=not "-noclean" in sys.argv[1:]
mainLog.write("====== clean="+`shouldClean`+" compile cp2k sopt ======\n")
mainLog.flush()
mainLog.write("  compilation logFile in 'cp2kBuildSopt.log'\n")
if not buildCp2k(cp2kRoot,"sopt",join(logDirPath,"cp2kBuildSopt.log"),
                 clean=shouldClean):
    mainLog.write("+++ ERROR, build FAILED! +++\n")
else:
    mainLog.write("+++ build SUCESSFULL! +++\n")

# do tests
exePath= join(cp2kRoot,"exe",
              commands.getoutput(join(cp2kRoot,"tools","get_arch_code")),
              "cp2k.sopt")
log1Path=join(logDirPath,"Ar.out")
os.chdir(join(cp2kRoot,"tests","QS"))
commands.getoutput("{ { "+exePath+" "+join(cp2kRoot,"tests","QS","Ar.inp")+
                   " ; } 2>&1 ; } > "+log1Path)
logFile=open(os.path.join(logDirPath,"Ar.diffs"),'w')
file1=open(join(cp2kRoot,"tests","QS","Ar.out"),'r')
file2=open(log1Path)
diffVal=diffEpsilon.compareCp2kOutput(file1,file2,
                                      0.0001,logFile)
file1.close(); file2.close()
logFile.write("totalDiff="+`diffVal`+"\n")
if diffVal>0.0001:
    mainLog.write("+++ ERROR, Argon test failed +++\n diff="+`diffVal`+
                  " more info in Argon.diffs\n")
else:
    mainLog.write("+++ Argon test SUCESSFULL (diff="+`diffVal`+")! +++\n")
logFile.close()


# cvs
if "-cvs" in sys.argv[1:]:
    os.chdir(cp2kRoot)
    # cvs -n update
    mainLog.write("====== cvs -n update ======\n")
    mainLog.flush()
    os.popen('{ { cvs -n update; } 2>&1 ; } >"%s"'%
             os.path.join(logDirPath,"cvsUpdate.log"))
    mainLog.write("  cvs update log in 'cvsUpdate.log'\n")

    # cvs diff
    mainLog.write("====== cvs diff ======\n")
    mainLog.flush()
    os.popen('{ { cvs diff; } 2>&1 ; } >"%s"'%
             os.path.join(logDirPath,"cvsDiff.log"))
    mainLog.write("  cvs diff log in 'cvsDiff.log'\n")

mainLog.close()
