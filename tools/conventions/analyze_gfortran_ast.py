#!/usr/bin/python
# -*- coding: utf-8 -*-

# author: Ole Schuett

import re
import sys
from os import path


BANNED_STM  = ('GOTO', 'FORALL', 'OPEN', 'CLOSE', 'STOP',)
BANNED_CALL = ('cp_fm_gemm', 'm_abort', 'mp_abort',)
USE_EXCEPTIONS = ("omp_lib", "omp_lib_kinds", "lapack",)

assert(all([x.upper()==x for x in BANNED_STM]))
assert(all([x.lower()==x for x in BANNED_CALL]))
assert(all([x.lower()==x for x in USE_EXCEPTIONS]))

# precompile regex
re_symbol    = re.compile(r"^\s*symtree.* symbol: '([^']+)'.*$")
re_use       = re.compile(r" USE-ASSOC\(([^)]+)\)")


#===============================================================================
def main():
    if(len(sys.argv) < 2):
        print("Usage: analyse_gfortran_ast.py <ast-file-1> ... <ast-file-N>")
        print("       This tool checks the given ASTs for violations of the coding conventions.\n")
        print("       For generating the abstract syntax tree (ast) run gfortran")
        print('       with "-fdump-fortran-original" and redirect output to file.')
        print('       This can be achieved by putting "FCLOGPIPE = >$(notdir $<).ast" in the cp2k arch-file.')
        sys.exit(1)

    log_files = sys.argv[1:]
    public_symbols = set()
    used_symbols = set()
    issues = []
    for fn in log_files:
        process_log_file(fn, public_symbols, used_symbols)

    #unused_public_symbols = public_symbols - used_symbols
    #print unused_public_symbols
    #if(len(unused_public_symbols)>0):
    #    issues.append("Found %d unused public symbols"%len(unused_public_symbols))

#===============================================================================
def process_log_file(fn, public_symbols, used_symbols):
    ast = open(fn).read()
    assert(fn.endswith(".ast"))
    loc = path.basename(fn)[:-4]
    lines = ast.split("\n")
    module_name = None

    curr_symbol = curr_procedure = stat_var = stat_stm = None
    skip_until_DT_END = False

    for line in lines:
        line = line.strip()
        tokens = line.split()

        if(skip_until_DT_END):
            # skip TRANSFERs which are part of READ/WRITE statement
            assert(tokens[0] in ("DO", "TRANSFER", "END", "DT_END"))
            if(tokens[0] == "DT_END"):
                skip_until_DT_END = False

        elif(stat_var):
            if(stat_var in line): # Ok, something was done with stat_var
                stat_var = stat_stm = None # reset
            elif(line=="ENDIF"):
                pass # skip, check may happen in outer scope
            elif(stat_stm=="ALLOCATE" and line.split()[0] == "ASSIGN"):
                pass # skip lines, it's part of the ALLOCATE statment
            else:
                print(loc+': Found %s with unchecked STAT in "%s"'%(stat_stm,curr_procedure))
                stat_var = stat_stm = None # reset

        elif(line.startswith("procedure name =")):
            curr_procedure = line.split("=")[1].strip()
            if(not module_name):
                module_name = curr_procedure

        elif(line.startswith("symtree: ") or len(line)==0):
            curr_symbol = None
            if(len(line)==0): continue
            curr_symbol = re_symbol.match(line).group(1)

        elif(line.startswith("attributes:")):
            if("USE-ASSOC" in line):
                mod = re_use.search(line).group(1)
                used_symbols.add(mod+"::"+curr_symbol)
                if("MODULE  USE-ASSOC" in line and mod.lower() not in USE_EXCEPTIONS):
                    print(loc+': Module "'+mod+'" USEd without ONLY clause or not PRIVATE')
            #if(("SAVE" in line) and ("PARAMETER" not in line) and ("PUBLIC" in line)):
            #    print(loc+': Symbol "'+curr_symbol+'" in procedure "'+curr_procedure+'" is PUBLIC-SAVE')
            if(("IMPLICIT-SAVE" in line) and ("PARAMETER" not in line) and ("USE-ASSOC" not in line) and (curr_procedure != module_name)):
                print(loc+': Symbol "'+curr_symbol+'" in procedure "'+curr_procedure+'" is IMPLICIT-SAVE')
            if(("IMPLICIT-TYPE" in line) and ("USE-ASSOC" not in line) and ("FUNCTION" not in line)): #TODO sure about last clause?
                print(loc+': Symbol "'+curr_symbol+'" in procedure "'+curr_procedure+'" is IMPLICIT-TYPE')
            if("THREADPRIVATE" in line):
                print(loc+': Symbol "'+curr_symbol+'" in procedure "'+curr_procedure+'" is THREADPRIVATE')
            if("PUBLIC" in line):
                public_symbols.add(module_name+"::"+curr_symbol)

        elif(line.startswith("!$OMP PARALLEL")):
            if("DEFAULT(NONE)" not in line):
                print(loc+': OMP PARALLEL without DEFAULT(NONE) found in "'+curr_procedure+'"')

        elif(line.startswith("CALL")):
            if(tokens[1].lower() in BANNED_CALL):
                print(loc+": Found CALL "+tokens[1]+' in procedure "'+curr_procedure+'"')
            elif(tokens[1].lower().startswith("_gfortran_arandom_")):
                print(loc+': Found CALL RANDOM_NUMBER in procedure "'+curr_procedure+'"')
            elif(tokens[1].lower().startswith("_gfortran_random_seed_")):
                print(loc+': Found CALL RANDOM_SEED in procedure "'+curr_procedure+'"')

        elif(tokens and tokens[0] in BANNED_STM):
            print(loc+": Found "+tokens[0]+' statement in procedure "'+curr_procedure+'"')

        elif(line.startswith("WRITE")):
            unit = tokens[1].split("=")[1]
            if(unit.isdigit()):
                print(loc+': Found WRITE statement with hardcoded unit in "'+curr_procedure+'"')

        elif(line.startswith("DEALLOCATE") and "STAT=" in line):
            print(loc+': Found DEALLOCATE with STAT argument in "'+curr_procedure+'"')

        elif("STAT=" in line): # catches also IOSTAT
            stat_var = line.split("STAT=",1)[1].split()[0]
            stat_stm = line.split()[0]
            skip_until_DT_END = stat_stm in ("READ", "WRITE",)


    # check for run-away DT_END search
    assert(skip_until_DT_END==False)
#===============================================================================
main()

#EOF
