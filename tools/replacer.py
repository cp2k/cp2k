#! /usr/bin/env python

import re
import sys

repl={
    'allocate_blacs_matrix':'cp_fm_create2',
    'allocate_blacs_matrix_vect':'cp_fm_vect_create2',
    'blacs_add':'cp_fm_add',
    'blacs_gemm':'cp_fm_gemm',
    'blacs_matrix_p_type':'cp_full_matrix_p_type',
    'blacs_matrix_type':'cp_full_matrix_type',
    'blacs_symm':'cp_fm_symm',
    'blacs_syrk':'cp_fm_syrk',
    'deallocate_blacs_matrix':'cp_fm_release',
    'deallocate_blacs_matrix_vect':'cp_fm_vect_dealloc',
    'get_blacs_matrix_info':'cp_fm_get_info',
    'sparse_plus_blacs_blacst':'cp_sm_plus_fm_fm_t',
    'sparse_times_blacs':'cp_sm_fm_multiply',
    }

specialRepl={re.compile(r"%blacs_matrix(\W)",flags=re.IGNORECASE):r"%matrix\1"}

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

