#! /usr/bin/env python
"""Functions that builds cp2k, rebuild templates, prettify,...."""

import sys, re, os, os.path, commands, time, shutil
from os.path import join
import prettify
import instantiateTemplates

defaultDirectives={"normalize-use":1,"upcase-keywords":1,"clean":1,
            "replace":0,"synopsis":1,"prettify-cvs":1}

def buildCp2k(cp2kRoot,buildType="sopt",logFilePath=None,clean=None):
    "builds cp2k"
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

def prettifyCp2k(cp2kRoot,buildType="sdbg",logDirPath=None,mainLog=sys.stdout,
                 directives=defaultDirectives):
    if not logDirPath:
        logDirPath=join(cp2kRoot,"prettify"+time.strftime("%y%m%d-%H:%M"))
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

    # build
    mainLog.write("====== dirty compiling cp2k "+buildType+" ======\n")
    mainLog.flush()
    mainLog.write("  compilation logFile in 'cp2kBuildSdbg.log'\n")
    if not buildCp2k(cp2kRoot,buildType,join(logDirPath,
                                             "cp2kBuild"+buildType+".log")):
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

