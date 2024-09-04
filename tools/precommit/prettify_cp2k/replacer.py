#!/usr/bin/env python3

import re
import sys

repl = {"routine_name": "routineN", "module_name": "moduleN"}
specialRepl = None
# { re.compile(r"(.*:: *moduleN) *= *(['\"])[a-zA-Z_0-9]+\2",flags=re.IGNORECASE):r"character(len=*), parameter :: moduleN = '__MODULE_NAME__'" }


def replaceWords(infile, outfile, replacements=repl, specialReplacements=specialRepl):
    """Replaces the words in infile writing the output to outfile.

    replacements is a dictionary with the words to replace.
    specialReplacements is a dictionary with general regexp replacements.
    """
    lineNr = 0
    nonWordRe = re.compile(r"(\W+)")

    while 1:
        line = infile.readline()
        lineNr = lineNr + 1
        if not line:
            break

        if specialReplacements:
            for subs in specialReplacements.keys():
                line = subs.sub(specialReplacements[subs], line)

        tokens = nonWordRe.split(line)
        for token in tokens:
            if token in replacements.keys():
                outfile.write(replacements[token])
            else:
                outfile.write(token)


# EOF
