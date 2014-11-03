#!/usr/bin/python
# -*- coding: utf-8 -*-

# author: Ole Schuett

import sys
import re

BANNED_STM  = ('GOTO', 'OPEN', 'CLOSE', )
BANNED_CALL = ('CP_FM_GEMM', )
USE_EXCEPTIONS = ("OMP_LIB", "OMP_LIB_KINDS", "F77_BLAS", "LAPACK",)

# false positives
UNDEF_EXCEPTIONS  = ('RANDOM_NUMBER', 'RANDOM_SEED', 'GET_COMMAND_ARGUMENT')

#BLAS routines
UNDEF_EXCEPTIONS += ('SROTG', 'DROTG', 'CROTG', 'ZROTG', 'SROTMG', 'DROTMG', 'SROT', 'DROT',
                     'ZROT', 'CSROT', 'ZDROT', 'SROTM', 'DROTM', 'SSWAP', 'DSWAP', 'CSWAP',
                     'ZSWAP', 'SSCAL', 'DSCAL', 'CSCAL', 'ZSCAL', 'CSSCAL', 'ZDSCAL', 'SCOPY',
                     'DCOPY', 'CCOPY', 'ZCOPY', 'SAXPY', 'DAXPY', 'CAXPY', 'ZAXPY', 'SDOT',
                     'DDOT', 'CDOTU', 'ZDOTU', 'CDOTC', 'ZDOTC', 'SNRM2', 'DNRM2', 'SCNRM2',
                     'DZNRM2', 'SASUM', 'SCASUM', 'DASUM', 'DZASUM', 'ISAMAX', 'IDAMAX', 'ICAMAX',
                     'IZAMAX',
                     'SGEMV', 'DGEMV', 'CGEMV', 'ZGEMV', 'SGBMV', 'DGBMV', 'CGBMV', 'ZGBMV',
                     'CHEMV', 'ZHEMV', 'CHBMV', 'ZHBMV', 'CHPMV', 'ZHPMV', 'SSYMV', 'DSYMV',
                     'SSBMV', 'DSBMV', 'SSPMV', 'DSPMV', 'STRMV', 'DTRMV', 'CTRMV', 'ZTRMV',
                     'STBMV', 'DTBMV', 'CTBMV', 'ZTBMV', 'STPMV', 'DTPMV', 'CTPMV', 'ZTPMV',
                     'STRSV', 'DTRSV', 'CTRSV', 'ZTRSV', 'STBSV', 'DTBSV', 'CTBSV', 'ZTBSV',
                     'STPSV', 'DTPSV', 'CTPSV', 'ZTPSV', 'SGER', 'DGER', 'CGERU', 'ZGERU',
                     'CGERC', 'ZGERC', 'CHER', 'ZHER', 'CHPR', 'ZHPR', 'CHER2', 'ZHER2',
                     'CHPR2', 'ZHPR2', 'SSYR', 'DSYR', 'SSPR', 'DSPR', 'SSYR2', 'DSYR2',
                     'SSPR2', 'DSPR2',
                     'SGEMM', 'DGEMM', 'CGEMM', 'ZGEMM', 'SSYMM', 'DSYMM', 'CSYMM', 'ZSYMM',
                     'CHEMM', 'ZHEMM', 'SSYRK', 'DSYRK', 'CSYRK', 'ZSYRK', 'CHERK', 'ZHERK',
                     'SSYR2K', 'DSYR2K', 'CSYR2K', 'ZSYR2K', 'CHER2K', 'ZHER2K', 'STRMM', 'DTRMM',
                     'CTRMM', 'ZTRMM', 'STRSM', 'DTRSM', 'CTRSM', 'ZTRSM', 'SDSDOT', 'DSDOT',
                     'DCABS1', 'LSAME', 'SCABS1')

#LAPACK routines
UNDEF_EXCEPTIONS += ('SGEEV', 'SLARNV', 'SPOTRF', 'SPOTRI',
                     'DGEEV', 'DGEQRF', 'DGESDD', 'DGESV', 'DGETRF', 'DGETRI',
                     'DGETRS', 'DLACPY', 'DLARNV', 'DPOTRF', 'DPOTRI', 'DSYEV', 'DSYEVD',
                     'DSYEVX', 'DSYGST', 'DTRTRI',
                     'CPOTRF', 'CPOTRI', 'CGEEV', 'CLARNV',
                     'ZGEEV', 'ZGETRF', 'ZGETRS', 'ZHEEVD', 'ZLARNV', 'ZPOTRF', 'ZPOTRI', 'ZTRTRI')

# precompile regex
re_symbol    = re.compile(r"^\s*symtree.* symbol: '([^']+)'.*$")
re_use       = re.compile(r" USE-ASSOC\(([^)]+)\)")


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

    public_symbols = set()
    used_symbols = set()
    issues = []
    for fn in log_files:
        issues += process_log_file(fn, public_symbols, used_symbols)

    #unused_public_symbols = public_symbols - used_symbols
    #print unused_public_symbols
    #if(len(unused_public_symbols)>0):
    #    issues.append("Found %d unused public symbols"%len(unused_public_symbols))

    issues = sorted(set(issues))
    issues_shown = [i for i in issues if(i not in suppress)]

    for i in issues_shown:
        print i

    n = len(issues_shown)
    m = len(issues) - n
    print "Summary: Found %d issues (%d suppressed)"%(n, m)
    print "Status: " + ("OK" if n==0 else "FAILED")

#===============================================================================
def process_log_file(fn, public_symbols, used_symbols):
    issues = []
    ast = open(fn).read()
    lines = ast.split("\n")
    module_name = None

    curr_symbol = curr_procedure = curr_symbol_defined = None

    for line in lines:
        line = line.strip()
        tokens = line.split()

        if(line.startswith("procedure name =")):
            curr_procedure = line.split("=")[1].strip()
            if(not module_name):
                module_name = curr_procedure

        elif(line.startswith("symtree: ") or len(line)==0):
            if(curr_symbol):
                if(not curr_symbol_defined and curr_symbol.upper() not in UNDEF_EXCEPTIONS):
                    issues.append(fn+': Symbol "'+curr_symbol+'" in "'+curr_procedure+'" not defined')

            curr_symbol = curr_symbol_defined = None
            if(len(line)==0): continue

            curr_symbol_defined = False
            curr_symbol = re_symbol.match(line).group(1)
            if("from namespace" in line):
                curr_symbol_defined = True

        elif(line.startswith("attributes:")):
            if("USE-ASSOC" in line):
                mod = re_use.search(line).group(1)
                used_symbols.add(mod+"::"+curr_symbol)
                curr_symbol_defined = True
                if("MODULE  USE-ASSOC" in line and mod.upper() not in USE_EXCEPTIONS):
                    issues.append(fn+': Module "'+mod+'" USEd without ONLY clause or not PRIVATE')
            #if(("SAVE" in line) and ("PARAMETER" not in line) and ("PUBLIC" in line)):
            #    issues.append(fn+': Symbol "'+curr_symbol+'" in procedure "'+curr_procedure+'" is PUBLIC-SAVE')
            if(("IMPLICIT-SAVE" in line) and ("PARAMETER" not in line) and ("USE-ASSOC" not in line) and (curr_procedure != module_name)):
                issues.append(fn+': Symbol "'+curr_symbol+'" in procedure "'+curr_procedure+'" is IMPLICIT-SAVE')
            if(("IMPLICIT-TYPE" in line) and ("USE-ASSOC" not in line) and ("FUNCTION" not in line)): #TODO sure about last clause?
                issues.append(fn+': Symbol "'+curr_symbol+'" in procedure "'+curr_procedure+'" is IMPLICIT-TYPE')
            if("INTRINSIC-PROC" in line):
                curr_symbol_defined = True
            if("INTRINSIC" in line):
                curr_symbol_defined = True
            if("INTERNAL-PROC" in line):
                curr_symbol_defined = True
            if("MODULE-PROC" in line):
                curr_symbol_defined = True
            if("EXTERNAL" in line):
                curr_symbol_defined = True
            if("MODULE" in line):
                curr_symbol_defined = True
            if("PROGRAM" in line):
                curr_symbol_defined = True
            if("LABEL" in line):
                curr_symbol_defined = True
            if("PUBLIC" in line):
                public_symbols.add(module_name+"::"+curr_symbol)

        elif(line.startswith("type spec :")):
            if("UNKNOWN" not in line):
                curr_symbol_defined = True
        elif(line.startswith("value:")):
            curr_symbol_defined = True
        elif(line.startswith("result:")):
            curr_symbol_defined = True
        elif(line.startswith("components:")):
            curr_symbol_defined = True
        elif(line.startswith("Generic interfaces:")):
            curr_symbol_defined = True
        elif(line.startswith("Formal arglist:")):
            curr_symbol_defined = True

        elif(line.startswith("!$OMP PARALLEL")):
            if("DEFAULT(NONE)" not in line):
                issues.append(fn+': OMP PARALLEL without DEFAULT(NONE) found in "'+curr_procedure+'"')

        elif(line.startswith("CALL")):
            if(tokens[1] in BANNED_CALL):
                issues.append(fn+": Found CALL "+tokens[1]+' in procedure "'+curr_procedure+'"')

        elif(tokens and tokens[0] in BANNED_STM):
            issues.append(fn+": Found "+tokens[0]+' statement in procedure "'+curr_procedure+'"')

    return(issues)

#===============================================================================
main()

#EOF
