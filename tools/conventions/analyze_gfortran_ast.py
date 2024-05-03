#!/usr/bin/env python3

# author: Ole Schuett

import argparse
import re
from os import path
from typing import TextIO, List


USE_EXCEPTIONS = ("omp_lib", "omp_lib_kinds", "lapack")

# precompile regex
re_symbol = re.compile(r"^\s*symtree.* symbol: '([^']+)'.*$")
re_derived_type = re.compile(r"^\s*type spec\s*:\s*\(DERIVED ([^']+)\).*$")
re_use = re.compile(r" USE-ASSOC\(([^)]+)\)")
re_conv = re.compile(
    r"__(real_[48]_r8|real_4_c8|cmplx1_4_r8_r8)\[\["
)  # ignores integers


# ======================================================================================
def process_log_file(fhandle: TextIO) -> None:
    public_symbols = set()
    used_symbols = set()

    def msg(message: str, conv_num: int) -> None:
        short_filename = path.basename(fhandle.name)[:-4]
        print(f"{short_filename}: {message} https://cp2k.org/conv#c{conv_num:03}")

    module_name = None

    cur_sym = cur_proc = cur_derived_type = cur_value = stat_var = stat_stm = None
    skip_until_DT_END = inside_omp_parallel = False

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
                msg(f'Found {stat_stm} with unchecked STAT in "{cur_proc}"', 1)
                stat_var = stat_stm = None  # reset

        elif line.startswith("procedure name ="):
            assert not inside_omp_parallel
            cur_proc = line.split("=")[1].strip()
            if not module_name:
                module_name = cur_proc

        elif line.startswith("symtree: ") or len(line) == 0:
            # Found new symbole, check default initializers of previous symbole.
            if cur_sym and cur_sym.startswith("__def_init_"):
                assert cur_derived_type
                ignore = cur_derived_type in ["c_ptr", "c_funptr"]
                initializer = cur_value or ""
                is_incomplete = any(x in initializer for x in [" () ", "(() ", " ())"])
                if not ignore and (not initializer or is_incomplete):
                    msg(f"Found type {cur_derived_type} without initializer", 16)

            # Parse new symbole.
            cur_sym = cur_derived_type = cur_value = None
            if len(line) == 0:
                continue
            match = re_symbol.match(line)
            assert match
            cur_sym = match.group(1)

        elif line.startswith("type spec"):
            # Only look for derived types because we want to check their initializers.
            match = re_derived_type.match(line)
            if match:
                cur_derived_type = match.group(1)

        elif line.startswith("value:"):
            # Hold onto values that could be passed to initializes of derived types.
            cur_value = line

        elif line.startswith("attributes:"):
            assert module_name and cur_sym
            is_imported = "USE-ASSOC" in line
            is_param = "PARAMETER" in line
            is_func = "FUNCTION" in line
            is_impl_save = "IMPLICIT-SAVE" in line
            is_impl_type = "IMPLICIT-TYPE" in line
            is_module_name = cur_proc == module_name

            if is_imported:
                match = re_use.search(line)
                assert match
                mod = match.group(1)
                used_symbols.add(mod + "::" + cur_sym)
                if "MODULE  USE-ASSOC" in line and mod.lower() not in USE_EXCEPTIONS:
                    msg(f'Module "{mod}" USEd without ONLY clause or not PRIVATE', 2)

            # if(("SAVE" in line) and ("PARAMETER" not in line) and ("PUBLIC" in line)):
            #    print(loc+': Symbol "'+cur_sym+'" in procedure "'+cur_proc+'" is PUBLIC-SAVE')

            if is_impl_save and not is_param and not is_imported and not is_module_name:
                msg(f'Symbol "{cur_sym}" in procedure "{cur_proc}" is IMPLICIT-SAVE', 3)

            if is_impl_type and not is_imported and not is_func:
                msg(f'Symbol "{cur_sym}" in procedure "{cur_proc}" is IMPLICIT-TYPE', 4)

            if "THREADPRIVATE" in line:
                msg(f'Symbol "{cur_sym}" in procedure "{cur_proc}" is THREADPRIVATE', 5)

            if "PUBLIC" in line:
                public_symbols.add(module_name + "::" + cur_sym)

        elif line.startswith("!$OMP PARALLEL"):
            assert not inside_omp_parallel
            inside_omp_parallel = True
            if "DEFAULT(NONE)" not in line:
                msg(f'OMP PARALLEL without DEFAULT(NONE) found in "{cur_proc}"', 6)

        elif line.startswith("!$OMP END PARALLEL"):
            assert inside_omp_parallel
            inside_omp_parallel = False

        elif line.startswith("ASSOCIATE") and inside_omp_parallel:
            text = f'Found ASSOCIATE statement inside OMP PARALLEL region in procedure "{cur_proc}"'
            msg(text, 8)

        elif line.startswith("CALL"):
            if "NULL()" in line:
                msg(f'Found CALL with NULL() as argument in procedure "{cur_proc}"', 7)
            elif tokens[1].lower() == "cp_fm_gemm":
                msg(f'Found CALL cp_fm_gemm in procedure "{cur_proc}"', 101)
            elif tokens[1].lower() == "m_abort":
                msg(f'Found CALL m_abort in procedure "{cur_proc}"', 102)
            elif tokens[1].lower() == "mp_abort":
                msg(f'Found CALL mp_abort in procedure "{cur_proc}"', 103)
            elif tokens[1].lower().startswith("_gfortran_arandom_"):
                msg(f'Found CALL RANDOM_NUMBER in procedure "{cur_proc}"', 104)
            elif tokens[1].lower().startswith("_gfortran_random_seed_"):
                msg(f'Found CALL RANDOM_SEED in procedure "{cur_proc}"', 105)
            elif tokens[1].lower().startswith("_gfortran_execute_command_line"):
                msg(f'Found CALL EXECUTE_COMMAND_LINE in procedure "{cur_proc}"', 106)

        elif line.startswith("GOTO"):
            msg(f'Found GOTO statement in procedure "{cur_proc}"', 201)
        elif line.startswith("FORALL"):
            msg(f'Found FORALL statement in procedure "{cur_proc}"', 202)
        elif line.startswith("OPEN"):
            msg(f'Found OPEN statement in procedure "{cur_proc}"', 203)
        elif line.startswith("CLOSE"):
            msg(f'Found CLOSE statement in procedure "{cur_proc}"', 204)
        elif line.startswith("STOP"):
            msg(f'Found STOP statement in procedure "{cur_proc}"', 205)

        elif line.startswith("WRITE"):
            unit = tokens[1].split("=")[1]
            if unit.isdigit():
                msg(f'Found WRITE statement with hardcoded unit in "{cur_proc}"', 12)

        elif line.startswith("DEALLOCATE") and "STAT=" in line:
            if ":ignore __final_" not in line:  # skip over auto-generated destructors
                msg(f'Found DEALLOCATE with STAT argument in "{cur_proc}"', 13)

        elif "STAT=" in line:  # catches also IOSTAT
            stat_var = line.split("STAT=", 1)[1].split()[0]
            stat_stm = line.split()[0]
            skip_until_DT_END = stat_stm in ("READ", "WRITE")

        elif "_gfortran_float" in line:
            msg(f'Found FLOAT in "{cur_proc}"', 14)

        elif re_conv.search(line):
            for m in re_conv.finditer(line):
                args = parse_args(line[m.end() :])
                if not re.match(r"\((kind = )?[48]\)", args[-1]):
                    text = f'Found lossy conversion {m.group(1)} without KIND argument in "{cur_proc}"'
                    msg(text, 15)

    # check for run-away DT_END search
    assert skip_until_DT_END is False


# ======================================================================================
def parse_args(line: str) -> List[str]:
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


# ======================================================================================
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

        with open(fn, encoding="utf8") as fhandle:
            process_log_file(fhandle)

# EOF
