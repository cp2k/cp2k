#! /usr/bin/env python

import re
import sys

repl={
    'force_control':'force_env_types',
    'fragment':'subsys',
    'cp_fragment_types':'cp_subsystem_types',
    'cp_fragment_type':'cp_subsystem_type',
    'cp_fragment_p_type':'cp_subsystem_p_type',
    'fragment_create':'cp_subsys_create',
    'fragment_retain':'cp_subsys_retain',
    'fragment_release':'cp_subsys_release',
    'fragment_get':'cp_subsys_get',
    'fragment_set':'cp_subsys_set',
    'qs_geoopt':'geo_opt',
    'qs_md':'md_run',
    'md_qs_energies':'md_energies',
    }

specialRepl=None
# { re.compile(r"%fragment",flags=re.IGNORECASE):r"%subsys" }

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

