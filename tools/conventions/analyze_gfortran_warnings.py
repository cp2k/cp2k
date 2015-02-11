#!/usr/bin/env python
# -*- coding: utf-8 -*-

# author: Ole Schuett

import re
import sys
from os import path

blas_re = re.compile("[SDCZ]"
                     +"(ROTG|ROTMG|ROT|ROTM|SWAP|SCAL|COPY|AXPY|DOT|DOTU|DOTC" #level 1
                     +"|GEMV|GBMV|HEMV|HBMV|HPMV|SYMV|SBMV|SPMV|TRMV|TBMV|TPMV|TRSV" # level 2
                     +"|TBSV|TPSV|GER|GERU|GERC|HER|HPR|HER2|HPR2|SYR|SPR|SYR2|SPR2"
                     +"|GEMM|SYMM|HEMM|SYRK|HERK|SYR2K|HER2K|TRMM|TRSM"  # level 3
                     +"|LANGE|LARNV|LAMCH|NRM2|CNRM2|ZNRM2|DSCAL)") # aux

lapack_re = re.compile("ILAENV|"
                       +"([SDCZ]" + "(BD|DI|GB|GE|GG|GT|HB|HE|HG|HP|HS|OP"
                       +"|OR|PB|PO|PP|PT|SB|SP|ST|SY|TB|TG|TP|TR|TZ|UN|UP)"
                       +"(BAK|BAL|BRD|CON|EBZ|EDC|EIN|EQR|EGR|EQU|EQZ|ERF|EVC"
                       +"|EXC|GBR|GHR|GLQ|GQL|GQR|GRQ|GST|GTR|HRD|LQF|MBR|MHR"
                       +"|MLQ|MQL|MQR|MRQ|MRZ|MTR|QLF|QPF|QRF|RFS|RQF|RZF|SDC"
                       +"|SEN|SJA|SNA|SQR|SVP|SYL|TRD|TRF|TRI|TRS"
                       +"|SDD|EV|GV|SV|BS2D|BR2D))")


#===============================================================================
def main():
    if(len(sys.argv) < 2):
        print("Usage: analyse_gfortran_warnings.py <warn-file-1> ... <warn-file-N>")
        print("       This tool checks the given stderr output from gfortran for violations of the coding conventions.")
        print("       For generating the warn-files run gfortran with as all warning flags and redirect the output to a file.")
        print('       This can be achieved by putting "FCLOGPIPE = 2>$(notdir $<).warn" in the cp2k arch-file.')
        sys.exit(1)

    files = sys.argv[1:]
    for fn in files:
        check_warnings(fn)

#===============================================================================
def check_warnings(fn):
    content = open(fn).read()

    lines = content.split("\n")
    loc = loc_short = None
    for i, line in enumerate(lines):
        if(len(line)==0): continue
        if(line[0]=="/" and line[-1]==":"):
            loc = line.rsplit(":", 2)[0].strip()
            loc_short = path.basename(loc)

        # we are only looking for warnings
        if(not line.startswith("Warning:")): continue

        # we ignore these warnings
        if("Creating array temporary" in line): continue
        if("quality comparison" in line): continue
        if("Unused" in line): 
		if("error" in line or "routinep" in line):
			continue
        if("defined but not used" in line): continue
        if("Removing call to function" in line): continue
        if("Conversion from" in line): continue
        if("CHARACTER expression" in line and "truncated" in line): continue
        if("POINTER-valued function appears on right-hand side" in line): continue

        # ok this warning we should handle
        if("called with an implicit interface" in line):
            parts = line.split()
            assert(parts[1] == "Procedure")
            routine = parts[2].strip("'").upper()
            if(may_call_implicit(loc, routine)): continue
            print "%s: Routine %s called with an implicit interface."%(loc_short, routine)
        else:
            print "%s: %s"%(loc_short, line) # unkown warning, just output

#===============================================================================
def may_call_implicit(loc, routine):
    if(blas_re.match(routine)):
        return(True)  # BLAS calls are allowed everywhere
    if(lapack_re.match(routine)):
        return(True) # Lapack calls are allowed everywhere

    pkg = path.dirname(loc)
    manifest_fn = pkg+"/PACKAGE"
    manifest = eval(open(manifest_fn).read())
    if(not manifest.has_key("implicit")): return(False)

    return(re.match(manifest["implicit"], routine))

#===============================================================================
main()
#EOF
