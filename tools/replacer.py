#! /usr/bin/env python

import re
import sys

repl={
    'routine_name':'routineN',
    'module_name':'moduleN'
    }
specialRepl=None
# { re.compile(r"(.*:: *moduleN) *= *(['\"])[a-zA-Z_0-9]+\2",flags=re.IGNORECASE):r"character(len=*), parameter :: moduleN = '__MODULE_NAME__'" }

def replaceWords(infile,outfile,replacements=repl,
                 specialReplacements=specialRepl,
                 logFile=sys.stdout):
    """Replaces the words in infile writing the output to outfile.
    
    replacements is a dictionary with the words to replace.
    specialReplacements is a dictionary with general regexp replacements.
    """
    lineNr=0
    nonWordRe=re.compile(r"(\W+)")
    
    while 1:
        line= infile.readline()
        lineNr=lineNr+1
        if not line: break
        
        if specialReplacements:
            for subs in specialReplacements.keys():
                line=subs.sub(specialReplacements[subs],line)
                
        tokens=nonWordRe.split(line)
        for token in tokens:
            if replacements.has_key(token):
                outfile.write(replacements[token])
            else:
                outfile.write(token)
        

if __name__ == '__main__':
    if len(sys.argv)<2:
        print "usage:", sys.argv[0]," file_in file_out"
    else:
        infile=open(sys.argv[1],'r')
        outfile=open(sys.argv[2],'w')
        replaceWords(infile,outfile,replacements=repl,
                     specialReplacements=specialRepl)

