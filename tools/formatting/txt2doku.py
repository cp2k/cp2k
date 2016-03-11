#!/usr/bin/python
# -*- coding: utf-8 -*-

# Converts a readable text file into a DokuWiki page
#
# author: Ole Schuett

import sys
import re
from os.path import normpath

#-------------------------------------------------------------------------------
def main():
    if(len(sys.argv) != 3):
        print("Usage: txt2doku.py <input-file> <output-file>")
        sys.exit(1)
    in_fn, out_fn = sys.argv[1:]
    assert(normpath(in_fn) != normpath(out_fn))

    text = open(in_fn).read()

    # format compiler flags as monospaced
    text = re.sub('(-D_[-_=<>A-Za-z0-9]+)', r"''%%\1%%''", text)

    # Create shiny note-boxes
    text = re.sub(r'(\nNote .*?\n)(?=\n)', note_box, text, flags=re.DOTALL)


    # fix multi-line list items (dokuwiki requires each item to be a single line)
    lines_out = []
    got_item = False
    for line in text.split("\n"):
        # start of item
        if(re.match("^  [-*] [^ ].*$", line)):
            lines_out.append(linkify(line))
            got_item = True

        # item continuation
        elif(got_item and re.match("^    [^ ].*$", line)):
            lines_out[-1] += " "+linkify(line.lstrip())

         # code-line as item continuation
        elif(got_item and re.match("^      [^ ].*$", line)):
            lines_out[-1] += " <code>"+line.strip()+"</code>"

        # normal code line
        elif(re.match("^  [^ ].*$", line)):
            lines_out.append(line)
            got_item = False

        # normal line
        else:
            lines_out.append(linkify(line))
            got_item = False

    text = "\n".join(lines_out)
    text = text.replace("</code> <code>", "\n")

    f = open(out_fn, "w")
    f.write(text)
    f.close()
    print("Wrote "+out_fn)

#-------------------------------------------------------------------------------
def note_box(m):
    txt = m.group(1).strip()
    return("\n<note important>\n"+txt+"\n</note>\n")

#-------------------------------------------------------------------------------
def linkify(line):
    return re.sub(r'(?<= )(cp2k/[a-zA-Z0-9_/.]+)(?=$| )', r"[[ src>\1 ]]", line)

#-------------------------------------------------------------------------------
if(len(sys.argv)==2 and sys.argv[-1]=="--selftest"):
    pass #TODO implement selftest
else:
    main()
#EOF
