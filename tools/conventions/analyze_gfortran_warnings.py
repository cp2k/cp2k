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

warning_re = re.compile(".*[Ww]arning: (.*)")


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
    loc = loc_short = ""
    for i, line in enumerate(lines):
        if(len(line)==0): continue
        if(line[0]=="/" and line[-1]==":"):
            loc = line.rsplit(":", 2)[0].strip()
            loc_short = path.basename(loc)

        if(loc.endswith("include/fftw3.f")): continue # an external file

        m = warning_re.match(line)
        if(not m): continue # we are only looking for warnings
        warning = m.group(1)

        # we ignore these warnings
        if("-Wmaybe-uninitialized" in warning): continue
        if("Creating array temporary" in warning): continue
        if("quality comparison" in warning): continue
        if("Unused" in warning and ("'error'" in warning or "'routinep'" in warning)): continue
        if("defined but not used" in warning): continue
        if("Removing call to function" in warning): continue
        if("Conversion from" in warning): continue
        if("CHARACTER expression" in warning and "truncated" in warning): continue
        if("POINTER-valued function appears on right-hand side" in warning): continue

        # ok this warning we should handle
        if("called with an implicit interface" in warning):
            parts = warning.split()
            assert(parts[0] == "Procedure")
            routine = parts[1].strip("'").upper()
            if(may_call_implicit(loc, routine)): continue
            print "%s: Routine %s called with an implicit interface."%(loc_short, routine)
        else:
            print "%s: %s"%(loc_short, warning) # unknown warning, just output

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
