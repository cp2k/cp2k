#!/usr/bin/env python

import sys
import re, tempfile
import os, os.path
from formatting import normalizeFortranFile
from formatting import replacer
from formatting import addSynopsis
from sys import argv

operatorsStr=r"\.(?:and|eqv?|false|g[et]|l[et]|n(?:e(?:|qv)|ot)|or|true)\."

keywordsStr="(?:a(?:llocat(?:able|e)|ssign(?:|ment))|c(?:a(?:ll|se)|haracter|lose|o(?:m(?:mon|plex)|nt(?:ains|inue))|ycle)|d(?:ata|eallocate|imension|o(?:|uble))|e(?:lse(?:|if|where)|n(?:d(?:|do|file|if)|try)|quivalence|x(?:it|ternal))|f(?:or(?:all|mat)|unction)|goto|i(?:f|mplicit|n(?:clude|quire|t(?:e(?:ger|nt|rface)|rinsic)))|logical|module|n(?:amelist|one|ullify)|o(?:nly|p(?:en|erator|tional))|p(?:a(?:rameter|use)|ointer|r(?:ecision|i(?:nt|vate)|o(?:cedure|gram))|ublic)|re(?:a[dl]|cursive|sult|turn|wind)|s(?:ave|e(?:lect|quence)|top|ubroutine)|t(?:arget|hen|ype)|use|w(?:h(?:ere|ile)|rite))"

intrinsic_procStr=r"(?:a(?:bs|c(?:har|os)|djust[lr]|i(?:mag|nt)|ll(?:|ocated)|n(?:int|y)|s(?:in|sociated)|tan2?)|b(?:it_size|test)|c(?:eiling|har|mplx|o(?:njg|sh?|unt)|shift)|d(?:ate_and_time|ble|i(?:gits|m)|ot_product|prod)|e(?:oshift|psilon|xp(?:|onent))|f(?:loor|raction)|huge|i(?:a(?:char|nd)|b(?:clr|its|set)|char|eor|n(?:dex|t)|or|shftc?)|kind|l(?:bound|en(?:|_trim)|g[et]|l[et]|og(?:|10|ical))|m(?:a(?:tmul|x(?:|exponent|loc|val))|erge|in(?:|exponent|loc|val)|od(?:|ulo)|vbits)|n(?:earest|int|ot)|p(?:ack|r(?:e(?:cision|sent)|oduct))|r(?:a(?:dix|n(?:dom_(?:number|seed)|ge))|e(?:peat|shape)|rspacing)|s(?:ca(?:le|n)|e(?:lected_(?:int_kind|real_kind)|t_exponent)|hape|i(?:gn|nh?|ze)|p(?:acing|read)|qrt|um|ystem_clock)|t(?:anh?|iny|r(?:ans(?:fer|pose)|im))|u(?:bound|npack)|verify)(?= *\()"

toUpcaseRe=re.compile("(?<![A-Za-z0-9_%#])(?<!% )(?P<toUpcase>"+operatorsStr+
                      "|"+ keywordsStr +"|"+ intrinsic_procStr +
                      ")(?![A-Za-z0-9_%])",flags=re.IGNORECASE)
linePartsRe=re.compile("(?P<commands>[^\"'!]*)(?P<comment>!.*)?"+
                       "(?P<string>(?P<qchar>[\"']).*?(?P=qchar))?")

def upcaseStringKeywords(line):
    """Upcases the fortran keywords, operators and intrinsic routines
    in line"""
    res=""
    start=0
    while start<len(line):
        m=linePartsRe.match(line[start:])
        if not m: raise SyntaxError("Syntax error, open string")
        res=res+toUpcaseRe.sub(lambda match: match.group("toUpcase").upper(),
                               m.group("commands"))
        if m.group("comment"):
            res=res+m.group("comment")
        if m.group("string"):
            res=res+m.group("string")
        start=start+m.end()
    return res

def upcaseKeywords(infile,outfile,logFile=sys.stdout):
    """Writes infile to outfile with all the fortran keywords upcased"""
    lineNr=0
    try:
        while 1:
            line=infile.readline()
            lineNr=lineNr+1
            if not line: break
            outfile.write(upcaseStringKeywords(line))
    except SyntaxError, e:
        e.lineno=lineNr
        e.text=line
        raise

def prettifyFile(infile,normalize_use=1, upcase_keywords=1,
             interfaces_dir=None,replace=None,logFile=sys.stdout):
    """prettifyes the fortran source in infile into a temporary file that is
    returned. It can be the same as infile.
    if normalize_use normalizes the use statements (defaults to true)
    if upcase_keywords upcases the keywords (defaults to true)
    if interfaces_dir is defined (and contains the directory with the
    interfaces) updates the synopsis
    if replace does the replacements contained in replacer.py (defaults
    to false)

    does not close the input file"""
    ifile=infile
    orig_filename=infile.name
    tmpfile=None
    try:
        if replace:
            tmpfile2=os.tmpfile()
            replacer.replaceWords(ifile,tmpfile2,logFile=logFile)
            tmpfile2.seek(0)
            if tmpfile:
                tmpfile.close()
            tmpfile=tmpfile2
            ifile=tmpfile
        if normalize_use:
            tmpfile2=os.tmpfile()
            normalizeFortranFile.rewriteFortranFile(ifile,tmpfile2,logFile,
                                                    orig_filename=orig_filename)
            tmpfile2.seek(0)
            if tmpfile:
                tmpfile.close()
            tmpfile=tmpfile2
            ifile=tmpfile
        if upcase_keywords:
            tmpfile2=os.tmpfile()
            upcaseKeywords(ifile,tmpfile2,logFile)
            tmpfile2.seek(0)
            if tmpfile:
                tmpfile.close()
            tmpfile=tmpfile2
            ifile=tmpfile
        if interfaces_dir:
            fileName=os.path.basename(infile.name)
            fileName=fileName[:fileName.rfind(".")]
            try:
                interfaceFile=open(os.path.join(interfaces_dir,
                                                fileName+".int"),"r")
            except:
                logFile.write("error opening file "+
                              os.path.join(interfaces_dir,
                                           fileName+".int")+"\n")
                logFile.write("skipping addSynopsis step for "+fileName+"\n")
                interfaceFile=None
            if interfaceFile:
                tmpfile2=os.tmpfile()
                addSynopsis.addSynopsisToFile(interfaceFile,ifile,
                                              tmpfile2,logFile=logFile)
                tmpfile2.seek(0)
                if tmpfile:
                    tmpfile.close()
                tmpfile=tmpfile2
                ifile=tmpfile
        return ifile
    except:
        logFile.write("error processing file '"+infile.name+"'\n")
        raise

def prettfyInplace(fileName,bkDir="preprettify",normalize_use=1,
                   upcase_keywords=1, interfaces_dir=None,
                   replace=None,logFile=sys.stdout):
    """Same as prettify, but inplace, replaces only if needed"""
    if not os.path.exists(bkDir):
        os.mkdir(bkDir)
    if not os.path.isdir(bkDir):
        raise Error("bk-dir must be a directory, was "+bkDir)
    infile=open(fileName,'r')
    outfile=prettifyFile(infile, normalize_use,
                         upcase_keywords, interfaces_dir, replace)
    if (infile==outfile):
        return
    infile.seek(0)
    outfile.seek(0)
    same=1
    while 1:
        l1=outfile.readline()
        l2=infile.readline()
        if (l1!=l2):
            same=0
            break
        if not l1:
            break
    if (not same):
        bkName=os.path.join(bkDir,os.path.basename(fileName))
        bName=bkName
        i=0
        while os.path.exists(bkName):
            i+=1
            bkName=bName+"."+str(i)
        infile.seek(0)
        bkFile=file(bkName,"w")
        while 1:
            l1=infile.readline()
            if not l1: break
            bkFile.write(l1)
        bkFile.close()
        outfile.seek(0)
        newFile=file(fileName,'w')
        while 1:
            l1=outfile.readline()
            if not l1: break
            newFile.write(l1)
        newFile.close()
    infile.close()
    outfile.close()

                   
if __name__ == '__main__':
    defaultsDict={'upcase':1,'normalize-use':1,'replace':1,
                  'interface-dir':None,
                  'backup-dir':'preprettify'}
    usageDesc=("usage:\n"+sys.argv[0]+ """
    [--[no-]upcase] [--[no-]normalize-use] [--[no-]replace]
    [--interface-dir=~/cp2k/obj/platform/target] [--help]
    [--backup-dir=bk_dir] file1 [file2 ...]

    replaces file1,... with their prettified version after performing on
    them upcase of the fortran keywords, and normalizion the use statements.
    If the interface direcory is given updates also the synopsis.
    If requested the replacements performed by the replacer.py script
    are also preformed.
    """+str(defaultsDict))
    
    replace=None
    if "--help" in sys.argv:
        print usageDesc
        sys.exit(0)
    args=[]
    for arg in sys.argv[1:]:
        m=re.match(r"--(no-)?(normalize-use|upcase|replace)",arg)
        if m:
            defaultsDict[m.groups()[1]]=not m.groups()[0]
        else:
            m=re.match(r"--(interface-dir|backup-dir)=(.*)",arg)
            if m:
                path=os.path.abspath(os.path.expanduser(m.groups()[1]))
                defaultsDict[m.groups()[0]]=path
            else:
                args.append(arg)
    if len(args)<1:
        print usageDesc
    else:
        bkDir=defaultsDict['backup-dir']
        if not os.path.exists(bkDir):
            # Another parallel running instance might just have created the dir.
            try:
                os.mkdir(bkDir)
            except:
                assert(os.path.exists(bkDir))
        if not os.path.isdir(bkDir):
            print "bk-dir must be a directory"
            print usageDesc
        else:
            failure=0
            for fileName in args:
                if not os.path.isfile(fileName):
                    print "file",fileName,"does not exists!"
                else:
                    try:
                        prettfyInplace(fileName,bkDir,
                                   normalize_use=defaultsDict['normalize-use'],
                                   upcase_keywords=defaultsDict['upcase'],
                                   interfaces_dir=defaultsDict['interface-dir'],
                                   replace=defaultsDict['replace'])
                    except:
                        failure+=1
                        import traceback
                        sys.stdout.write('-'*60+"\n")
                        traceback.print_exc(file=sys.stdout)
                        sys.stdout.write('-'*60+"\n")
                        sys.stdout.write("Processing file '"+fileName+"'\n")
            sys.exit(failure>0)
