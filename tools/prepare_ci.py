#! /usr/bin/env python
# prepares cp2k for checkin

import sys, re, os, os.path, commands, time, shutil
from sys import argv
from os.path import join
import prettify
import instantiateTemplates
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
    pipe=os.popen("{ { "+makeCmd+" VERSION="+buildType+"; } 2>&1 ; } >>"+logFilePath)
    if (pipe.close()):
        logFile=open(logFilePath,'a')
        logFile.write("\n+++ ERROR, build "+buildType+" FAILED! +++\n")
        logFile.close()
        return None
    else:
        logFile=open(logFilePath,'a')
        logFile.write("\n+++ build "+buildType+" SUCESSFULLY! +++\n")
        logFile.close()
        return 1


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
mainLog.write(" ******* prepare check-in BEGAN *******\n")
mainLog.flush()
outDir= join(cp2kRoot,"src","outDir")

print "logDirectory: "+logDirPath
print "main log: "+mainLog.name

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

# prettifyes files
mainLog.write("====== prettifying files ======\n")
mainLog.flush()
logFile=open(os.path.join(logDirPath,"prettify.log"),'w')
buildDir= join(cp2kRoot,"obj",
    commands.getoutput(join(cp2kRoot,"tools","get_arch_code")),
    "sdbg")
outDir=join(logDirPath,"prettifyDir")
if os.access(outDir,os.W_OK): commands.getoutput('rm -rf "'+outDir+'"')
os.mkdir(outDir)

filesToPret2=[] # glob.glob(os.path.join(cp2kRoot,"src","cp_*.F"))
filesToPret=templateInstances
baseNames=map(os.path.basename,filesToPret)
for fileToPret in filesToPret2:
    if not os.path.basename(fileToPret) in baseNames:
	filesToPret.append(fileToPret)

# if requested adds cvs modified files to files to prettify
if directives["prettify-cvs"]:
    mainLog.write("+ adding cvs modified files to the files to prettify\n")
    mainLog.flush()
    os.chdir(os.path.join(cp2kRoot,"src"))
    filesC=commands.getoutput("cvs -n update")
    shouldUpdate=0
    filesToPret2=[]
    conflicts=0
    conflictsDir=os.path.join(outDir,"conflicts")
    fileCRe=re.compile(r"([ACRMUP]?) ([a-zA-Z_\.\-]+\.F)$")
    
    for line in filesC.splitlines():
        m=fileCRe.match(line)
        if m:
            if m.groups()[0] in ["A","M","C"]:
                filesToPret2.append(os.path.join(cp2kRoot,"src",m.groups()[1]))
            if m.groups()[0]=="C":
                conflicts=conflicts+1
                if not os.path.isdir(conflictsDir): os.mkdir(conflictsDir)
                shutil.copyfile(m.groups()[1],
                                os.path.join(conflictsDir,m.groups()[1]))
                logFile.write(" copied "+m.groups()[1]+" to "+
                              conflictsDir+"\n")
            if m.groups()[0] in ["U","P","C"]:
                shouldUpdate=1

    baseNames=map(os.path.basename,filesToPret)
    for fileToPret in filesToPret2:
        if not os.path.basename(fileToPret) in baseNames:
            filesToPret.append(fileToPret)
    
    logFile.write("cvs modified file to prettify:\n"+str(filesToPret2)+"\n")
    if shouldUpdate:
        mainLog.write("++ WARNING cvs update is needed\n")
        mainLog.write("++ there were "+str(conflicts)+" conflicts\n")
        if conflicts:
            mainLog.write("++ consider restoring the original files from\n"+
                          "++   "+conflictsDir+"\n")
        mainLog.flush()

# start prettyfication
logFile.write("Files to prettify:\n"+str(filesToPret)+"\n")
logFile.write("\nStarting prettification\n")
logFile.flush()
mainLog.write("+ starting prettify process\n")
mainLog.flush()

if not directives["synopsis"]: buildDir=None
errors=0
for fileP in filesToPret:
    try:
        logFile.write("\n+ processing file '"+os.path.basename(fileP)+"'\n")
        infile=open(fileP,"r")
        outfile=open(os.path.join(outDir,os.path.basename(fileP)),"w")
        prettify.prettifyFile(infile,outfile,
                              normalize_use=directives["normalize-use"],
                              upcase_keywords=directives["upcase-keywords"],
                              interfaces_dir=buildDir,
                              replace=directives["replace"],
                              logFile=logFile)
        infile.close()
        outfile.close()
    except:
        logFile.write("\n** ERROR prettifying the file "+fileP+"\n")
        import traceback
        logFile.write('-'*60+"\n")
        traceback.print_exc(file=logFile)
        logFile.write('-'*60+"\n")
        logFile.flush()
        errors=errors+1
        if fileP in templateInstances:
            shutil.copyfile(fileP,os.path.join(outDir,os.path.basename(fileP)))
            logFile.write("+ used unprettified template instance for file "+
                          os.path.basename(fileP)+"\n")
os.mkdir(os.path.join(outDir,"orig"))
os.chdir(outDir)
for fileP in os.listdir(outDir):
    try:
        if os.path.isfile(fileP) and fileP[-4:]!=".err":
            logFile.write("moving of file "+fileP)
            origF=os.path.join(cp2kRoot,"src",os.path.basename(fileP))
            if os.path.isfile(origF):
                if commands.getoutput('diff "'+origF+'" "'+fileP+'"'):
                    os.rename(origF,os.path.join(outDir,"orig",
                                                 os.path.basename(fileP)))
                    os.rename(fileP,origF)
                else:
                    logFile.write(" NOT")
            else:
                os.rename(fileP,origF)
            logFile.write(" done\n")
    except:
        logFile.write("\n** ERROR moving the file "+
                      os.path.basename(fileP)+"\n")
        import traceback
        logFile.write('-'*60+"\n")
        traceback.print_exc(file=logFile)
        logFile.write('-'*60+"\n")
        logFile.flush()
        errors=errors+1

logFile.write(" ***** prettification finished *****\n")
logFile.write("+ there were "+str(errors)+" errors\n")
logFile.close()
mainLog.write("+ there were "+str(errors)+" errors while prettifying\n")
mainLog.write("  prettification logFile in '%s'\n"%os.path.basename(logFile.name))
mainLog.flush()
# os.rename returns before the operation is complete, so wait a little 
# (use command.getoutput intread of os.rename?)
time.sleep(3.0)

# clean? compile sopt
mainLog.write("====== clean="+str(directives["clean"])+
              " compile cp2k sopt ======\n")
mainLog.flush()
mainLog.write("  compilation logFile in 'cp2kBuildSopt.log'\n")
if not buildCp2k(cp2kRoot,"sopt",join(logDirPath,"cp2kBuildSopt.log"),
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
    if not buildCp2k(cp2kRoot,"popt",join(logDirPath,"cp2kBuildPopt.log"),
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


