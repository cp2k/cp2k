#! /usr/bin/env python

import sys
import re
import string
from sys import argv

def parseUse(inFile):
    """Parses the use statements in inFile
    The parsing stops at the first non use statement.
    Returns something like:
    ([{'module':'module1','only':['el1',el2=>el3']},...],
     '! comment1\\n!comment2...\\n',
     'last line (the line that stopped the parsing)')
    """
    useStartRe=re.compile(
        r" *(?P<use>use[^&!]*)(?P<continue>&?) *(?P<comment>!.*)?$",
        flags=re.IGNORECASE)
    commentRe=re.compile(r" *!.*$")
    contLineRe=re.compile(
        r"(?P<contLine>[^&!]*)(?P<continue>&?) *(?P<comment>!.*)?$")
    useParseRe=re.compile(
        r" *use +(?P<module>[a-zA-Z_][a-zA-Z_0-9]*)(?P<only> *, *only *:)? *(?P<imports>.*)$",
        flags=re.IGNORECASE)
    lineNr=0
    comments=""
    modules=[]
    line=""
    while 1:
        line=inFile.readline()
        lineNr=lineNr+1
        if not line: break

        m=useStartRe.match(line)
        if m:
            # read whole use
            compactedUse=m.group('use')
            useComments=""
            if m.group('comment'): useComments=m.group('comment')+'\n'
            while m.group('continue'):
                lineNr=lineNr+1
                m=contLineRe.match(inFile.readline())
                compactedUse=compactedUse+m.group('contLine')
                if m.group('comment'):
                    useComments=useComments+m.group('comment')+'\n'
            # parse use
            m=useParseRe.match(compactedUse)
            if not m:
                raise SyntaxError("could not parse use ending at line "+
                                  str(lineNr)+" (compactedUse="+compactedUse+
                                  ")")
            useAtt={'module':m.group('module')}
            if m.group('only'):
                useAtt['only']=map(string.strip,
                                   string.split(m.group('imports'),','))
            else:
                useAtt['renames']=map(string.strip,
                                      string.split(m.group('imports'),','))
		if useAtt['renames']==[""]: del useAtt['renames']
            if useComments : useAtt['comments']=useComments
            # add use to modules
            modules.append(useAtt)
        elif commentRe.match(line):
            comments=comments+line
        elif line and not line.isspace():
            break
    return (modules,comments,line)

def normalizeModules(modules):
    """Sorts the modules and their export and removes duplicates.
    renames aren't sorted correctly"""
    # orders modules
    modules.sort(lambda x,y:cmp(x['module'],y['module']) )
    for i in range(len(modules)-1,0,-1):
        if modules[i]['module']==modules[i-1]['module']:
            if not (modules[i-1].has_key('only') and
                    modules[i].has_key('only')):
                raise SyntaxError('rejoining of module '+
                                  str(modules[i]['module'])+
                                  ' failed as at least one of the use is not a use ...,only:')
            modules[i-1]['only'].extend(modules[i]['only'])
            del modules[i]
    # orders imports
    for m in modules:
        if m.has_key('only'):
            m['only'].sort()
            for i in range(len(m['only'])-1,0,-1):
                if m['only'][i-1]==m['only'][i]: del m['only'][i]

def writeUseLong(modules,outFile):
    for m in modules:
        if m.has_key('only'):
            outFile.write("  USE "+m['module']+","+
                          string.rjust('ONLY: ',38-len(m['module'])))
            if m['only']: outFile.write(m['only'][0])
            for i in range(1,len(m['only'])):
                outFile.write(",&\n"+string.ljust("",45)+m['only'][i])
        else:
            outFile.write("  USE "+m['module'])
            if m.has_key('renames') and m['renames']:
                outFile.write(","+string.ljust("",38)+
                              m['renames'][0])
                for i in range(1,len(m['renames'])):
                    outFile.write(",&\n"+string.ljust("",45)+m['renames'][i])
        if m.has_key('comments') and m['comments']:
            comments=m['comments'].splitlines()
            outFile.write("&")
            for i in range(0,len(comments)-1):
                outFile.write("\n&"+comments[i])
            outFile.write("\n"+comments[-1])
        outFile.write("\n")

def cleanUse(modules,rest):
    """Removes the unneded modules (the ones that are not used in rest)"""
    exceptions={"cp_a_l":1,"cp_to_string":1,"cp_error_type":1,"cp_assert":1,
                "cp_failure_level":1,"cp_warning_level":1,"cp_note_level":1,
                "cp_fatal_level":1,"cp_logger_type":1,"timeset":1,"timestop":1,
                "dp":1,"cp_error_get_logger":1, "cp_error_message":1}
    rest=rest.lower()
    for i in range(len(modules)-1,-1,-1):
        if modules[i].has_key("only"):
            els=modules[i]['only']
            for j in range(len(els)-1,-1,-1):
                if not exceptions.has_key(els[j].lower()):
                    if rest.find(els[j].lower())==-1:
                        del els[j]
            if len(modules[i]['only'])==0:
                del modules[i]

def rewriteUse(inFile,outFile,logFile=sys.stdout):
    """rewrites the use statements of in file to outFile.
    It sorts them and removes the repetitions."""
    import os.path
    moduleRe=re.compile(r" *module (?P<moduleName>[a-zA-Z_][a-zA-Z_0-9]*) *(?:!.*)?$",
                        flags=re.IGNORECASE)
    while 1:
        line=inFile.readline()
        if not line: break
        outFile.write(line)
        m=moduleRe.match(line)
        if m:
            if (m.group('moduleName')!=os.path.basename(inFile.name)[0:-2]) :
                raise SyntaxError("Module name is different from filename ("+
                                  m.group('moduleName')+
                                  "!="+os.path.basename(inFile.name)[0:-2]+")")
            break
    try:
        (modules,comments,line)=parseUse(inFile)
        rest=line+inFile.read()
        cleanUse(modules,rest)
        outFile.write(comments)
        normalizeModules(modules)
        writeUseLong(modules,outFile)
        outFile.write(rest)
    except:
        import traceback
        logFile.write('-'*60+"\n")
        traceback.print_exc(file=logFile)
        logFile.write('-'*60+"\n")

        logFile.write("Processing file '"+inFile.name+"'\n")
        

if __name__ == '__main__':
    import os.path
    if len(sys.argv)<2:
        print "usage:", sys.argv[0]," out_dir file1 [file2 ...]"
    else:
        outDir=sys.argv[1]
        if not os.path.isdir(outDir):
            print "out_dir must be a directory"
            print "usage:", sys.argv[0]," out_dir file1 [file2 ...]"
        else:
            for fileName in sys.argv[2:]:
                infile=open(fileName,'r')
                outfile=open(os.path.join(outDir,
                                          os.path.basename(fileName)),'w')
                rewriteUse(infile,outfile)
