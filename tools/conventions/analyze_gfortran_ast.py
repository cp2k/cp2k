#!/usr/bin/env python3

# author: Ole Schuett

import argparse
import re
from os import path


BANNED_STM = ("GOTO", "FORALL", "OPEN", "CLOSE", "STOP")
BANNED_CALL = ("cp_fm_gemm", "m_abort", "mp_abort")
USE_EXCEPTIONS = ("omp_lib", "omp_lib_kinds", "lapack")

# precompile regex
re_symbol = re.compile(r"^\s*symtree.* symbol: '([^']+)'.*$")
re_use = re.compile(r" USE-ASSOC\(([^)]+)\)")
re_conv = re.compile(
    r"__(real_[48]_r8|real_4_c8|cmplx1_4_r8_r8)\[\["
)  # ignores integers


def process_log_file(fhandle):
    public_symbols = set()
    used_symbols = set()

    def lprint(*args, **kwargs):
        return print("{}:".format(path.basename(fhandle.name)[:-4]), *args, **kwargs)

    module_name = None

    curr_symbol = curr_procedure = stat_var = stat_stm = None
    skip_until_DT_END = False

    for line in fhandle:
        line = line.strip()
        tokens = line.split()

        if skip_until_DT_END:
            # skip TRANSFERs which are part of READ/WRITE statement
            assert tokens[0] in ("DO", "TRANSFER", "END", "DT_END")
            if tokens[0] == "DT_END":
                skip_until_DT_END = False

        elif stat_var:
            if stat_var in line:  # Ok, something was done with stat_var
                stat_var = stat_stm = None  # reset
            elif line == "ENDIF":
                pass  # skip, check may happen in outer scope
            elif stat_stm == "ALLOCATE" and tokens[0] == "ASSIGN":
                pass  # skip lines, it's part of the ALLOCATE statement
            else:
                lprint(
                    'Found %s with unchecked STAT in "%s"' % (stat_stm, curr_procedure)
                )
                stat_var = stat_stm = None  # reset

        elif line.startswith("procedure name ="):
            curr_procedure = line.split("=")[1].strip()
            if not module_name:
                module_name = curr_procedure

        elif line.startswith("symtree: ") or len(line) == 0:
            curr_symbol = None
            if len(line) == 0:
                continue
            curr_symbol = re_symbol.match(line).group(1)

        elif line.startswith("attributes:"):
            if "USE-ASSOC" in line:
                mod = re_use.search(line).group(1)
                used_symbols.add(mod + "::" + curr_symbol)
                if "MODULE  USE-ASSOC" in line and mod.lower() not in USE_EXCEPTIONS:
                    lprint(
                        'Module "{}" USEd without ONLY clause or not PRIVATE'.format(
                            mod
                        )
                    )
            # if(("SAVE" in line) and ("PARAMETER" not in line) and ("PUBLIC" in line)):
            #    print(loc+': Symbol "'+curr_symbol+'" in procedure "'+curr_procedure+'" is PUBLIC-SAVE')
            if (
                ("IMPLICIT-SAVE" in line)
                and ("PARAMETER" not in line)
                and ("USE-ASSOC" not in line)
                and (curr_procedure != module_name)
            ):
                lprint(
                    'Symbol "'
                    + curr_symbol
                    + '" in procedure "'
                    + curr_procedure
                    + '" is IMPLICIT-SAVE'
                )
            if (
                ("IMPLICIT-TYPE" in line)
                and ("USE-ASSOC" not in line)
                and ("FUNCTION" not in line)
            ):  # TODO sure about last clause?
                lprint(
                    'Symbol "'
                    + curr_symbol
                    + '" in procedure "'
                    + curr_procedure
                    + '" is IMPLICIT-TYPE'
                )
            if "THREADPRIVATE" in line:
                lprint(
                    'Symbol "'
                    + curr_symbol
                    + '" in procedure "'
                    + curr_procedure
                    + '" is THREADPRIVATE'
                )
            if "PUBLIC" in line:
                public_symbols.add(module_name + "::" + curr_symbol)

        elif line.startswith("!$OMP PARALLEL"):
            if "DEFAULT(NONE)" not in line:
                lprint(
                    'OMP PARALLEL without DEFAULT(NONE) found in "'
                    + curr_procedure
                    + '"'
                )

        elif line.startswith("CALL"):
            if tokens[1].lower() in BANNED_CALL:
                lprint(
                    "Found CALL " + tokens[1] + ' in procedure "' + curr_procedure + '"'
                )
            elif tokens[1].lower().startswith("_gfortran_arandom_"):
                lprint('Found CALL RANDOM_NUMBER in procedure "' + curr_procedure + '"')
            elif tokens[1].lower().startswith("_gfortran_random_seed_"):
                lprint('Found CALL RANDOM_SEED in procedure "' + curr_procedure + '"')

        elif tokens and tokens[0] in BANNED_STM:
            lprint(
                "Found "
                + tokens[0]
                + ' statement in procedure "'
                + curr_procedure
                + '"'
            )

        elif line.startswith("WRITE"):
            unit = tokens[1].split("=")[1]
            if unit.isdigit():
                lprint(
                    'Found WRITE statement with hardcoded unit in "'
                    + curr_procedure
                    + '"'
                )

        elif line.startswith("DEALLOCATE") and "STAT=" in line:
            if ":ignore __final_" not in line:  # skip over auto-generated destructors
                lprint(
                    'Found DEALLOCATE with STAT argument in "' + curr_procedure + '"'
                )

        elif "STAT=" in line:  # catches also IOSTAT
            stat_var = line.split("STAT=", 1)[1].split()[0]
            stat_stm = line.split()[0]
            skip_until_DT_END = stat_stm in ("READ", "WRITE")

        elif "_gfortran_float" in line:
            lprint('Found FLOAT in "' + curr_procedure + '"')

        elif re_conv.search(line):
            for m in re_conv.finditer(line):
                args = parse_args(line[m.end() :])
                if not re.match(r"\((kind = )?[48]\)", args[-1]):
                    lprint(
                        'Found lossy conversion %s without KIND argument in "%s"'
                        % (m.group(1), curr_procedure)
                    )

    # check for run-away DT_END search
    assert skip_until_DT_END is False

    return (public_symbols, used_symbols)


def parse_args(line):
    assert line[0] == "("
    parentheses = 1
    args = list()
    for i in range(1, len(line)):
        if line[i] == "(":
            if parentheses == 1:
                a = i  # beginning of argument
            parentheses += 1
        elif line[i] == ")":
            parentheses -= 1
            if parentheses == 1:  # end of argument
                args.append(line[a : i + 1])
            if parentheses == 0:
                return args

    raise Exception("Could not find matching parentheses")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Checks the given ASTs for violations of the coding conventions",
        epilog="""\
For generating the abstract syntax tree (ast) run gfortran
with "-fdump-fortran-original" and redirect output to file.
This can be achieved by putting
    FCLOGPIPE = >$(notdir $<).ast
in the cp2k arch-file.
""",
    )
    parser.add_argument(
        "files",
        metavar="<ast-file>",
        type=str,
        nargs="+",
        help="files containing dumps of the AST",
    )
    args = parser.parse_args()

    for fn in args.files:
        assert fn.endswith(".ast")

        with open(fn) as fhandle:
            process_log_file(fhandle)
