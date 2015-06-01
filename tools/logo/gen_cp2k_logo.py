#!/usr/bin/python
# -*- coding: utf-8 -*-

# Generates the CP2K Logo
#
# author: Ole Schuett

import os

#-------------------------------------------------------------------------------
def main():
    gen_povray()

    # run povray
    for res in (100,300, 500, 800):
        cmd = "povray -D +UA +H%d +W%d +Q11 +A +Ocp2k_logo_%d.png logo.pov"%(res,res, res)
        print "Running: "+cmd
        os.system(cmd)

#-------------------------------------------------------------------------------
def gen_povray():
    txt = ["XXXX XXXX      X  X",
           "X    X  X XXXX X X ",
           "X    X  X    X XX  ",
           "X    XXXX    X XX  ",
           "X    X    XXXX X X ",
           "XXXX X    X    X  X",
           "          X        ",
           "          XXXX     "]

    coords = []
    for z in range(3):
        for y, line in enumerate(txt):
            for x, c in enumerate(line):
                if(c=="X"):
                    coords.append((x-9,-y+3, -z))

    output  = '#include "colors.inc"\n'
    output += '#include "textures.inc"\n'
    output += 'global_settings { assumed_gamma 1.0 }\n'

    output += 'light_source { <-100, -70, -300> White }\n'
    output += 'light_source { <40, -20, 50> White }\n'

    output += 'camera {\n'
    output += '  location <0.52, -9.1, -28.6>\n'
    output += '  sky   <0.15, 1, 0>\n'
    output += '  angle 40\n'
    output += '  look_at <0,0,0>\n'
    output += '}\n'

    for c in coords:
        output += "sphere { <%f, %f, %f>, 0.73 \n"%c
        output += "  texture {\n"
        output += "    finish { Shiny ambient 0.3 }\n"
        output += "    pigment { color rgb <1, 0.2, 0> }\n"
        output += "  }\n"
        output += "}\n"

    f = open("logo.pov", "w")
    f.write(output)
    f.close()


#-------------------------------------------------------------------------------
main()
#EOF
