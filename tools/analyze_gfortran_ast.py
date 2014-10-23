#!/usr/bin/python
# -*- coding: utf-8 -*-

# author: Ole Schuett

import sys
import re

BANNED_STM  = ('GOTO', 'OPEN', 'CLOSE', )

#taken from f77_blas_poison.F
BANNED_CALL  = ('SROTG',  'DROTG', 'CROTG',  'ZROTG', 'SROTMG', 'DROTMG', 'SROT',   'DROT',
                'ZROT',   'CSROT', 'ZDROT',  'SROTM', 'DROTM',  'SSWAP',  'DSWAP',  'CSWAP',
                'ZSWAP',  'SSCAL', 'DSCAL',  'CSCAL', 'ZSCAL',  'CSSCAL', 'ZDSCAL', 'SCOPY',
                'DCOPY',  'CCOPY', 'ZCOPY',  'SAXPY', 'DAXPY',  'CAXPY',  'ZAXPY',  'SDOT',
                'DDOT',   'CDOTU', 'ZDOTU',  'CDOTC', 'ZDOTC',  'SNRM2',  'DNRM2',  'SCNRM2',
                'DZNRM2', 'SASUM', 'SCASUM', 'DASUM', 'DZASUM', 'ISAMAX', 'IDAMAX', 'ICAMAX',
                'IZAMAX')

BANNED_CALL += ('SGEMV', 'DGEMV', 'CGEMV', 'ZGEMV', 'SGBMV', 'DGBMV', 'CGBMV', 'ZGBMV',
                'CHEMV', 'ZHEMV', 'CHBMV', 'ZHBMV', 'CHPMV', 'ZHPMV', 'SSYMV', 'DSYMV',
                'SSBMV', 'DSBMV', 'SSPMV', 'DSPMV', 'STRMV', 'DTRMV', 'CTRMV', 'ZTRMV',
                'STBMV', 'DTBMV', 'CTBMV', 'ZTBMV', 'STPMV', 'DTPMV', 'CTPMV', 'ZTPMV',
                'STRSV', 'DTRSV', 'CTRSV', 'ZTRSV', 'STBSV', 'DTBSV', 'CTBSV', 'ZTBSV',
                'STPSV', 'DTPSV', 'CTPSV', 'ZTPSV', 'SGER',  'DGER',  'CGERU', 'ZGERU',
                'CGERC', 'ZGERC', 'CHER',  'ZHER',  'CHPR',  'ZHPR',  'CHER2', 'ZHER2',
                'CHPR2', 'ZHPR2', 'SSYR',  'DSYR',  'SSPR',  'DSPR',  'SSYR2', 'DSYR2',
                'SSPR2', 'DSPR2')

BANNED_CALL += ('SGEMM',  'DGEMM',  'CGEMM',  'ZGEMM',  'SSYMM',  'DSYMM',  'CSYMM',  'ZSYMM',
                'CHEMM',  'ZHEMM',  'SSYRK',  'DSYRK',  'CSYRK',  'ZSYRK',  'CHERK',  'ZHERK',
                'SSYR2K', 'DSYR2K', 'CSYR2K', 'ZSYR2K', 'CHER2K', 'ZHER2K', 'STRMM',  'DTRMM',
                'CTRMM',  'ZTRMM',  'STRSM',  'DTRSM',  'CTRSM',  'ZTRSM',  'SDSDOT', 'DSDOT',
                'DCABS1', 'LSAME',  'SCABS1')

BANNED_CALL += ('CP_FM_GEMM', )

# precompile regex
re_procedure = re.compile(r"^\s*procedure name = (.*)$")
re_symbol    = re.compile(r"^\s*symtree.* symbol: '([^']+)'.*$")
re_attr      = re.compile(r"^\s*attributes: (.*)$")

#===============================================================================
def main():
    if(len(sys.argv) < 2):
        print("Usage: analyse_gfortran_ast.py [--suppressions=<supp-file>] <ast-file-1> ... <ast-file-N>")
        print("       This tool checks the given ASTs for violations of the coding conventions.\n")
        print("       For generating the abstract syntax tree (ast) run gfortran")
        print('       with "-fdump-fortran-original" and redirect output to file.')
        print('       This can be achieved by putting "FCLOGPIPE = >$(notdir $<).ast" in the cp2k arch-file.')
        sys.exit(1)

    suppress = []
    log_files = sys.argv[1:]
    if sys.argv[1].startswith("--suppressions="):
        content = open(sys.argv[1].split("=")[1]).read()
        lines = [l.strip() for l in content.split("\n")]
        suppress = [l for l in lines if not l.startswith("#")]
        log_files = sys.argv[2:]

    issues = []
    for fn in log_files:
        issues += process_log_file(fn)

    issues = sorted(set(issues))
    issues_shown = [i for i in issues if(i not in suppress)]

    for i in issues_shown:
        print i

    n = len(issues_shown)
    m = len(issues) - n
    print "Summary: Found %d issues (%d suppressed)"%(n, m)
    print "Status: " + ("OK" if n==0 else "FAILED")

#===============================================================================
def process_log_file(fn):
    issues = []
    ast = open(fn).read()
    lines = ast.split("\n")

    curr_symbol = curr_procedure = None
    for line in lines:
        m = re_procedure.match(line)
        if(m):
            curr_procedure = m.group(1)
            continue

        m = re_symbol.match(line)
        if(m):
            curr_symbol = m.group(1)
            continue

        m = re_attr.match(line)
        if(m and ("IMPLICIT-SAVE" in line) and ("PARAMETER" not in line)):
            issues.append(fn+': Symbol "'+curr_symbol+'" in procedure "'+curr_procedure+'" is IMPLICIT-SAVE')
            continue

        m = re_attr.match(line)
        if(m and ("IMPLICIT-TYPE" in line)):
            issues.append(fn+': Symbol "'+curr_symbol+'" in procedure "'+curr_procedure+'" is IMPLICIT-TYPE')
            continue

        tokens = line.upper().strip().split()
        if(tokens and tokens[0] in BANNED_STM):
            issues.append(fn+": Found "+tokens[0]+' statement in procedure "'+curr_procedure+'"')
            continue

        if(tokens and tokens[0]=="CALL" and tokens[1] in BANNED_CALL):
            issues.append(fn+": Found CALL "+tokens[1]+' in procedure "'+curr_procedure+'"')
            continue

        if(tokens and tokens[0]=="!$OMP"):
            if("OMP PARALLEL" in line and "DEFAULT(NONE)" not in line):
                issues.append(fn+': OMP PARALLEL without DEFAULT(NONE) found in "'+curr_procedure+'"')

    return(issues)

#===============================================================================
main()

#EOF
