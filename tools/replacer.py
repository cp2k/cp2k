#! /usr/bin/env python

import re
import sys

repl={
    'blacs_cholesky_decompose':'cp_fm_cholesky_decompose',
    'blacs_triangular_multiply':'cp_fm_triangular_multiply',
    'blacs_cholesky_invert':'cp_fm_cholesky_invert',
    'blacs_cholesky_reduce':'cp_fm_cholesky_reduce',
    'blacs_cholesky_restore':'cp_fm_cholesky_restore',
    'cp_fm_pool':'cp_fm_pool_types',
    'cp_f_matrix_struct':'cp_fm_struct',
    'cp_full_matrix':'cp_fm_types'
    }

specialRepl={}
# {re.compile(r"%blacs_matrix(\W)",flags=re.IGNORECASE):r"%matrix\1"}

def replaceWords(infile,outfile,replacements,specialReplacements=None,
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

