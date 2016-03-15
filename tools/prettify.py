#!/usr/bin/env python

import sys
import re, tempfile
import os, os.path
import md5
from formatting import normalizeFortranFile
from formatting import replacer
from formatting import addSynopsis
from formatting import reformatFortranFile
from sys import argv

operatorsStr=r"\.(?:and|eqv?|false|g[et]|l[et]|n(?:e(?:|qv)|ot)|or|true)\."

keywordsStr="(?:a(?:llocat(?:able|e)|ssign(?:|ment))|c(?:a(?:ll|se)|haracter|lose|o(?:m(?:mon|plex)|nt(?:ains|inue))|ycle)|d(?:ata|eallocate|imension|o(?:|uble))|e(?:lse(?:|if|where)|n(?:d(?:|do|file|if)|try)|quivalence|x(?:it|ternal))|f(?:or(?:all|mat)|unction)|goto|i(?:f|mplicit|n(?:clude|quire|t(?:e(?:ger|nt|rface)|rinsic)))|logical|module|n(?:amelist|one|ullify)|o(?:nly|p(?:en|erator|tional))|p(?:a(?:rameter|use)|ointer|r(?:ecision|i(?:nt|vate)|o(?:cedure|gram))|ublic)|re(?:a[dl]|cursive|sult|turn|wind)|s(?:ave|e(?:lect|quence)|top|ubroutine)|t(?:arget|hen|ype)|use|w(?:h(?:ere|ile)|rite))"

intrinsic_procStr=r"(?:a(?:bs|c(?:har|os)|djust[lr]|i(?:mag|nt)|ll(?:|ocated)|n(?:int|y)|s(?:in|sociated)|tan2?)|b(?:it_size|test)|c(?:eiling|har|mplx|o(?:njg|sh?|unt)|shift)|d(?:ate_and_time|ble|i(?:gits|m)|ot_product|prod)|e(?:oshift|psilon|xp(?:|onent))|f(?:loor|raction)|huge|i(?:a(?:char|nd)|b(?:clr|its|set)|char|eor|n(?:dex|t)|or|shftc?)|kind|l(?:bound|en(?:|_trim)|g[et]|l[et]|og(?:|10|ical))|m(?:a(?:tmul|x(?:|exponent|loc|val))|erge|in(?:|exponent|loc|val)|od(?:|ulo)|vbits)|n(?:earest|int|ot)|p(?:ack|r(?:e(?:cision|sent)|oduct))|r(?:a(?:dix|n(?:dom_(?:number|seed)|ge))|e(?:peat|shape)|rspacing)|s(?:ca(?:le|n)|e(?:lected_(?:int_kind|real_kind)|t_exponent)|hape|i(?:gn|nh?|ze)|p(?:acing|read)|qrt|um|ystem_clock)|t(?:anh?|iny|r(?:ans(?:fer|pose)|im))|u(?:bound|npack)|verify)(?= *\()"

ompDir=r"(?:atomic|barrier|c(?:apture|ritical)|do|end|flush|if|master|num_threads|ordered|parallel|read|s(?:ection(?:|s)|ingle)|t(?:ask(?:|wait|yield)|hreadprivate)|update|w(?:orkshare|rite)|!\$omp)"

ompClause=r"(?:a|co(?:llapse|py(?:in|private))|default|fi(?:nal|rstprivate)|i(?:and|eor|or)|lastprivate|m(?:ax|ergeable|in)|n(?:one|owait)|ordered|private|reduction|shared|untied|\.(?:and|eqv|neqv|or)\.)"

ompEnv=r"omp_(?:dynamic|max_active_levels|n(?:ested|um_threads)|proc_bind|s(?:tacksize|chedule)|thread_limit|wait_policy)"

# FIXME: does not correctly match operator '.op.' if it is not separated by whitespaces.
toUpcaseRe=re.compile("(?<![A-Za-z0-9_%#])(?<!% )(?P<toUpcase>"+operatorsStr+
                      "|"+ keywordsStr +"|"+ intrinsic_procStr +
                      ")(?![A-Za-z0-9_%])",flags=re.IGNORECASE)
toUpcaseOMPRe=re.compile("(?<![A-Za-z0-9_%#])(?P<toUpcase>"
                         +ompDir+"|"+ompClause+"|"+ompEnv + 
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

def upcaseOMP(line):
    """Upcases OpenMP stuff."""
    return toUpcaseOMPRe.sub(lambda match: match.group("toUpcase").upper(),line)

def upcaseKeywords(infile,outfile,upcase_omp,logFile=sys.stdout):
    """Writes infile to outfile with all the fortran keywords upcased"""
    while 1:
        line=infile.readline()
        if not line: break
        line = upcaseStringKeywords(line)
        if upcase_omp:
            if normalizeFortranFile.ompDirRe.match(line):
                line=upcaseOMP(line)
        outfile.write(line)

def prettifyFile(infile, normalize_use=1, reformat=0, indent=2, whitespace=2, upcase_keywords=1,
                 upcase_omp=0, interfaces_dir=None,replace=None,logFile=sys.stdout):
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
    max_pretty_iter = 5
    n_pretty_iter = 0

    while True:
        n_pretty_iter += 1
        hash_prev = md5.new()
        hash_prev.update(ifile.read())
        ifile.seek(0)
        try:
            if replace:
                tmpfile2=os.tmpfile()
                replacer.replaceWords(ifile,tmpfile2,logFile=logFile)
                tmpfile2.seek(0)
                if tmpfile:
                    tmpfile.close()
                tmpfile=tmpfile2
                ifile=tmpfile
            if reformat: # reformat needs to be done first
                tmpfile2=os.tmpfile()
                reformatFortranFile.reformat_ffile(ifile,tmpfile2,logFile=logFile,
                                                        indent_size=indent,whitespace=whitespace,
                                                        orig_filename=orig_filename)
                tmpfile2.seek(0)
                if tmpfile:
                    tmpfile.close()
                tmpfile=tmpfile2
                ifile=tmpfile
            if normalize_use:
                tmpfile2=os.tmpfile()
                normalizeFortranFile.rewriteFortranFile(ifile,tmpfile2,indent,logFile,
                                                        orig_filename=orig_filename)
                tmpfile2.seek(0)
                if tmpfile:
                    tmpfile.close()
                tmpfile=tmpfile2
                ifile=tmpfile
            if upcase_keywords:
                tmpfile2=os.tmpfile()
                upcaseKeywords(ifile,tmpfile2,upcase_omp,logFile)
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
            hash_next = md5.new()
            hash_next.update(ifile.read())
            ifile.seek(0)
            if hash_prev.digest() == hash_next.digest():
                return ifile
            elif n_pretty_iter >= max_pretty_iter:
                raise RuntimeError("Prettify did not converge in", max_pretty_iter, "steps.")
        except:
            logFile.write("error processing file '"+infile.name+"'\n")
            raise

def prettfyInplace(fileName,bkDir="preprettify",normalize_use=1,
                   reformat=0,indent=2,whitespace=2,
                   upcase_keywords=1, upcase_omp=0, interfaces_dir=None,
                   replace=None,logFile=sys.stdout):
    """Same as prettify, but inplace, replaces only if needed"""
    if not os.path.exists(bkDir):
        os.mkdir(bkDir)
    if not os.path.isdir(bkDir):
        raise Error("bk-dir must be a directory, was "+bkDir)
    infile=open(fileName,'r')
    outfile=prettifyFile(infile, normalize_use, reformat,
                         indent, whitespace, upcase_keywords, upcase_omp,
                         interfaces_dir, replace)
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

# FIXME: 'use' statements means not only use, but also variable declarations, add this to docu
def main():
    # future defaults
    defaultsDict={'upcase':1,'normalize-use':1,'omp-upcase':1,
                  'reformat':1, 'indent':3, 'whitespace':1,
                  'replace':1, 'interface-dir':None,
                  'backup-dir':'preprettify'}

    # current defaults (FIXME: change to future defaults)
    defaultsDict={'upcase':1,'normalize-use':1,'omp-upcase':0,
                  'reformat':0, 'indent':2, 'whitespace':1,
                  'replace':1, 'interface-dir':None,
                  'backup-dir':'preprettify'}

    usageDesc=("usage:\n"+sys.argv[0]+ """
    [--[no-]upcase] [--[no-]normalize-use] [--[no-]omp-upcase] [--[no-]replace]
    [--[no-]reformat] --indent=3 --whitespace=1 [--interface-dir=~/cp2k/obj/platform/target] [--help]
    [--backup-dir=bk_dir] file1 [file2 ...]

    Auto-format F90 source file1, file2, ...:
    - normalization of use statements (--normalize-use) (FIXME: give more detailed description)
    - upcasing fortran keywords (--upcase)
    - auto-indentation, auto-alignment and whitespace formatting (--reformat).
      Amount of whitespace controlled by --whitespace = 0, 1, 2. 
      For indenting with a relative width of n columns specify --indent=n.
      For manual formatting of specific lines:
      - disable auto-alignment by starting line continuation with an ampersand '&'
      - completely disable reformatting by adding a comment '!&'
      For manual formatting of a code block, use:
      - start a manually formatted block with a '!&<' comment and close it with a '!&>' comment
    - upcasing OMP directives (--omp-upcase)
    - If requested the replacements performed by the replacer.py script are also performed (--replace) (FIXME: ???)
    - If the interface direcory is given updates also the synopsis. (--interface-dir) (FIXME: ???)

    Defaults:
    """+str(defaultsDict))

    replace=None
    if "--help" in sys.argv:
        print(usageDesc)
        sys.exit(0)
    args=[]
    for arg in sys.argv[1:]:
        m=re.match(r"--(no-)?(normalize-use|upcase|omp-upcase|replace|reformat)",arg)
        if m:
            defaultsDict[m.groups()[1]]=not m.groups()[0]
        else:
            m=re.match(r"--(indent|whitespace)=(.*)",arg)
            if m:
                defaultsDict[m.groups()[0]]=int(m.groups()[1])
            else:
                m=re.match(r"--(interface-dir|backup-dir)=(.*)",arg)
                if m:
                    path=os.path.abspath(os.path.expanduser(m.groups()[1]))
                    defaultsDict[m.groups()[0]]=path
                else:
                    if arg.startswith('--'):
                        print('unknown option',arg)
                    else:
                        args.append(arg)
    if len(args)<1:
        print(usageDesc)
    else:
        bkDir=defaultsDict['backup-dir']
        if not os.path.exists(bkDir):
            # Another parallel running instance might just have created the dir.
            try:
                os.mkdir(bkDir)
            except:
                assert(os.path.exists(bkDir))
        if not os.path.isdir(bkDir):
            print("bk-dir must be a directory")
            print(usageDesc)
        else:
            failure=0
            for fileName in args:
                if not os.path.isfile(fileName):
                    print("file",fileName,"does not exists!")
                else:
                    try:
                        prettfyInplace(fileName,bkDir,
                                   normalize_use=defaultsDict['normalize-use'],
                                   reformat=defaultsDict['reformat'],
                                   indent=defaultsDict['indent'],
                                   whitespace=defaultsDict['whitespace'],
                                   upcase_keywords=defaultsDict['upcase'],
                                   upcase_omp=defaultsDict['omp-upcase'],
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


#===============================================================================
if(__name__ == '__main__'):
    if(len(sys.argv)==2 and sys.argv[-1]=="--selftest"):
        pass #TODO implement selftest
    else:
        main()
#EOF
