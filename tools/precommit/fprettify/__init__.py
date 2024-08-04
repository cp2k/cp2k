#!/usr/bin/env python
# -*- coding: utf-8 -*-
###############################################################################
#    This file is part of fprettify.
#    Copyright (C) 2016-2019 Patrick Seewald, CP2K developers group
#
#    fprettify is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    fprettify is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with fprettify. If not, see <http://www.gnu.org/licenses/>.
###############################################################################

"""
Impose white space conventions and indentation based on scopes / subunits

normalization of white spaces supported for following operators:
- relational operators:
  .EQ. .NE. .LT. .LE. .GT. .GE.
  ==   /=   <    <=    >   >=
- logical operators:
  .AND. .OR. .EQV. .NEQV.
  .NOT.
- bracket delimiters
- commas and semicolons:
- arithmetic operators:
  *  /  **  +  -
- other operators:
  %  - (sign)  = (function argument)
  = (assignment)  => (pointer assignment)

supported criteria for alignment / indentation:
 Fortran lines:
 - if, else, endif
 - do, enddo
 - select case, case, end select
 - select rank, rank, end select
 - subroutine, end subroutine
 - function, end function
 - module, end module
 - program, end program
 - interface, end interface
 - type, end type
 Actual lines (parts of Fortran lines separated by linebreaks):
 - bracket delimiters (.), (/./), and [.]
 - assignments by value = and pointer =>.

LIMITATIONS
- assumes that all subunits are explicitly ended within same file,
  no treatment of #include statements
- can not deal with f77 constructs (files are ignored)

FIXME's
- internal errors should not happen
- wrap regular expression parser. This allows to extend parser by constructs
  that are not regular expressions (and support e.g. forall construct).
- strip whitespaces once and for all and then assume no trailing / leading
  whitespaces
- open files only when needed
"""
import re
import sys
import logging
import os
import io

sys.stdin = io.TextIOWrapper(
    sys.stdin.detach(), encoding='UTF-8', line_buffering=True)
sys.stdout = io.TextIOWrapper(
    sys.stdout.detach(), encoding='UTF-8', line_buffering=True)


from .fparse_utils import (VAR_DECL_RE, OMP_COND_RE,
                           InputStream, CharFilter,
                           FprettifyException, FprettifyParseException, FprettifyInternalException,
                           CPP_RE, NOTFORTRAN_LINE_RE, NOTFORTRAN_FYPP_LINE_RE, FYPP_LINE_RE, RE_FLAGS,
                           STR_OPEN_RE, parser_re, FYPP_WITHOUT_PREPRO_RE)

# recognize fortran files by extension
FORTRAN_EXTENSIONS = [".f", ".for", ".ftn",
                      ".f90", ".f95", ".f03", ".fpp"]
FORTRAN_EXTENSIONS += [_.upper() for _ in FORTRAN_EXTENSIONS]

# constants, mostly regular expressions:
FORMATTER_ERROR_MESSAGE = (" Wrong usage of formatting-specific directives"
                           " '&', '!&', '!&<' or '!&>'.")
LINESPLIT_MESSAGE = ("auto indentation failed due to chars limit, "
                     "line should be split")

EOL_STR = r"\s*;?\s*$"  # end of fortran line
EOL_SC = r"\s*;\s*$"  # whether line is ended with semicolon
SOL_STR = r"^\s*"  # start of fortran line

STATEMENT_LABEL_RE = re.compile(r"^\s*(\d+)", RE_FLAGS)

# regular expressions for parsing statements that start, continue or end a
# subunit:
IF_RE = re.compile(
    SOL_STR + r"(\w+\s*:)?\s*IF\s*\(.*\)\s*THEN" + EOL_STR, RE_FLAGS)
ELSE_RE = re.compile(
    SOL_STR + r"ELSE(\s*IF\s*\(.*\)\s*THEN)?" + EOL_STR, RE_FLAGS)
ENDIF_RE = re.compile(SOL_STR + r"END\s*IF(\s+\w+)?" + EOL_STR, RE_FLAGS)

DO_RE = re.compile(SOL_STR + r"(\w+\s*:)?\s*DO(" + EOL_STR + r"|\s+\w)", RE_FLAGS)
ENDDO_RE = re.compile(SOL_STR + r"END\s*DO(\s+\w+)?" + EOL_STR, RE_FLAGS)

SELCASE_RE = re.compile(
    SOL_STR + r"SELECT\s*(CASE|RANK|TYPE)\s*\(.*\)" + EOL_STR, RE_FLAGS)
CASE_RE = re.compile(
    SOL_STR + r"((CASE|RANK|TYPE\s+IS|CLASS\s+IS)\s*(\(.*\)|DEFAULT)|CLASS\s+DEFAULT)" + EOL_STR, RE_FLAGS)
ENDSEL_RE = re.compile(SOL_STR + r"END\s*SELECT" + EOL_STR, RE_FLAGS)

ASSOCIATE_RE = re.compile(SOL_STR + r"ASSOCIATE\s*\(.*\)" + EOL_STR, RE_FLAGS)
ENDASSOCIATE_RE = re.compile(SOL_STR + r"END\s*ASSOCIATE" + EOL_STR, RE_FLAGS)

BLK_RE = re.compile(SOL_STR + r"(\w+\s*:)?\s*BLOCK" + EOL_STR, RE_FLAGS)
ENDBLK_RE = re.compile(SOL_STR + r"END\s*BLOCK(\s+\w+)?" + EOL_STR, RE_FLAGS)

SUBR_RE = re.compile(
    r"^([^\"']* )?SUBROUTINE\s+\w+\s*(\(.*\))?" + EOL_STR, RE_FLAGS)
ENDSUBR_RE = re.compile(
    SOL_STR + r"END\s*SUBROUTINE(\s+\w+)?" + EOL_STR, RE_FLAGS)

FCT_RE = re.compile(
    r"^([^\"']* )?FUNCTION\s+\w+\s*(\(.*\))?(\s*RESULT\s*\(\w+\))?" + EOL_STR,
    RE_FLAGS)
ENDFCT_RE = re.compile(
    SOL_STR + r"END\s*FUNCTION(\s+\w+)?" + EOL_STR, RE_FLAGS)

MOD_RE = re.compile(SOL_STR + r"MODULE\s+\w+" + EOL_STR, RE_FLAGS)
ENDMOD_RE = re.compile(SOL_STR + r"END\s*MODULE(\s+\w+)?" + EOL_STR, RE_FLAGS)

SMOD_RE = re.compile(SOL_STR + r"SUBMODULE\s*\(\w+\)\s+\w+" + EOL_STR, RE_FLAGS)
ENDSMOD_RE = re.compile(SOL_STR + r"END\s*SUBMODULE(\s+\w+)?" + EOL_STR, RE_FLAGS)

TYPE_RE = re.compile(
    SOL_STR +
    r"TYPE(\s*,\s*(BIND\s*\(\s*C\s*\)|EXTENDS\s*\(.*\)|ABSTRACT|PUBLIC|PRIVATE))*(\s*,\s*)?(\s*::\s*|\s+)\w+" + EOL_STR,
    RE_FLAGS)
ENDTYPE_RE = re.compile(SOL_STR + r"END\s*TYPE(\s+\w+)?" + EOL_STR, RE_FLAGS)

PROG_RE = re.compile(SOL_STR + r"PROGRAM\s+\w+" + EOL_STR, RE_FLAGS)
ENDPROG_RE = re.compile(
    SOL_STR + r"END\s*PROGRAM(\s+\w+)?" + EOL_STR, RE_FLAGS)

INTERFACE_RE = re.compile(
    r"^([^\"']* )?INTERFACE(\s+\w+|\s+(OPERATOR|ASSIGNMENT)\s*\(.*\))?" + EOL_STR, RE_FLAGS)
ENDINTERFACE_RE = re.compile(
    SOL_STR + r"END\s*INTERFACE(\s+\w+|\s+(OPERATOR|ASSIGNMENT)\s*\(.*\))?" + EOL_STR, RE_FLAGS)

CONTAINS_RE = re.compile(SOL_STR + r"CONTAINS" + EOL_STR, RE_FLAGS)

ENUM_RE = re.compile(
    SOL_STR + r"ENUM(\s*,\s*(BIND\s*\(\s*C\s*\)))?((\s*::\s*|\s+)\w+)?" + EOL_STR,
    RE_FLAGS)
ENDENUM_RE = re.compile(SOL_STR + r"END\s*ENUM(\s+\w+)?" + EOL_STR, RE_FLAGS)

ENDANY_RE = re.compile(SOL_STR + r"END" + EOL_STR, RE_FLAGS)

# Regular expressions for where and forall block constructs
FORALL_RE = re.compile(SOL_STR + r"(\w+\s*:)?\s*FORALL\s*\(.*\)" + EOL_STR, RE_FLAGS)
ENDFORALL_RE = re.compile(SOL_STR + r"END\s*FORALL(\s+\w+)?" + EOL_STR, RE_FLAGS)

WHERE_RE = re.compile(SOL_STR + r"(\w+\s*:)?\s*WHERE\s*\(.*\)" + EOL_STR, RE_FLAGS)
ELSEWHERE_RE = re.compile(SOL_STR + r"ELSE\s*WHERE(\(.*\))?(\s*\w+)?" + EOL_STR, RE_FLAGS)
ENDWHERE_RE = re.compile(SOL_STR + r"END\s*WHERE(\s+\w+)?" + EOL_STR, RE_FLAGS)

# Regular expressions for preprocessor directives

FYPP_DEF_RE = re.compile(SOL_STR + r"#:DEF\s+", RE_FLAGS)
FYPP_ENDDEF_RE = re.compile(SOL_STR + r"#:ENDDEF", RE_FLAGS)

FYPP_IF_RE = re.compile(SOL_STR + r"#:IF\s+", RE_FLAGS)
FYPP_ELIF_ELSE_RE = re.compile(SOL_STR + r"#:(ELIF\s+|ELSE)", RE_FLAGS)
FYPP_ENDIF_RE = re.compile(SOL_STR + r"#:ENDIF", RE_FLAGS)

FYPP_FOR_RE = re.compile(SOL_STR + r"#:FOR\s+", RE_FLAGS)
FYPP_ENDFOR_RE = re.compile(SOL_STR + r"#:ENDFOR", RE_FLAGS)

FYPP_BLOCK_RE = re.compile(SOL_STR + r"#:BLOCK\s+", RE_FLAGS)
FYPP_ENDBLOCK_RE = re.compile(SOL_STR + r"#:ENDBLOCK", RE_FLAGS)

FYPP_CALL_RE = re.compile(SOL_STR + r"#:CALL\s+", RE_FLAGS)
FYPP_ENDCALL_RE = re.compile(SOL_STR + r"#:ENDCALL", RE_FLAGS)

FYPP_MUTE_RE = re.compile(SOL_STR + r"#:MUTE", RE_FLAGS)
FYPP_ENDMUTE_RE = re.compile(SOL_STR + r"#:ENDMUTE", RE_FLAGS)

PRIVATE_RE = re.compile(SOL_STR + r"PRIVATE\s*::", RE_FLAGS)
PUBLIC_RE = re.compile(SOL_STR + r"PUBLIC\s*::", RE_FLAGS)

END_RE = re.compile(SOL_STR + r"(END)\s*(IF|DO|SELECT|ASSOCIATE|BLOCK|SUBROUTINE|FUNCTION|MODULE|SUBMODULE|TYPE|PROGRAM|INTERFACE|ENUM|WHERE|FORALL)", RE_FLAGS)

# intrinsic statements with parenthesis notation that are not functions
INTR_STMTS_PAR = (r"(ALLOCATE|DEALLOCATE|"
                  r"OPEN|CLOSE|READ|WRITE|"
                  r"FLUSH|ENDFILE|REWIND|BACKSPACE|INQUIRE|"
                  r"FORALL|WHERE|ASSOCIATE|NULLIFY)")

# regular expressions for parsing linebreaks
LINEBREAK_STR = r"(&)[\s]*(?:!.*)?$"

# regular expressions for parsing operators
# Note: +/- in real literals and sign operator is ignored
PLUSMINUS_RE = re.compile(r"(?<=[\w\)\]])\s*(\+|-)\s*", RE_FLAGS)
# Note: ** or // (or any multiples of * or /) are ignored
#       we also ignore any * or / before a :: because we may be seeing 'real*8'
MULTDIV_RE = re.compile(
    r"(?<=[\w\)\]])\s*((?<!\*)\*(?!\*)|(?<!/)/(?!/))(?=[\s\w\(])(?!.*::)", RE_FLAGS)
REL_OP_RE = re.compile(
    r"(?<!\()\s*(\.(?:EQ|NE|LT|LE|GT|GE)\.|(?:==|\/=|<(?!=)|<=|(?<!=)>(?!=)|>=))\s*(?!\))",
    RE_FLAGS)
LOG_OP_RE = re.compile(r"\s*(\.(?:AND|OR|EQV|NEQV)\.)\s*", RE_FLAGS)
PRINT_RE = re.compile(r"(?:(?<=\bPRINT)|(?<=\bREAD))\s*(\*,?)\s*", RE_FLAGS)

# regular expressions for parsing delimiters
DEL_OPEN_STR = r"(\(\/?|\[)"
DEL_OPEN_RE = re.compile(r"^" + DEL_OPEN_STR, RE_FLAGS)
DEL_CLOSE_STR = r"(\/?\)|\])"
DEL_CLOSE_RE = re.compile(r"^" + DEL_CLOSE_STR, RE_FLAGS)

# empty line regex
EMPTY_RE = re.compile(SOL_STR + r"$", RE_FLAGS)

PREPRO_NEW_SCOPE = [parser_re(FYPP_DEF_RE), parser_re(FYPP_IF_RE), parser_re(FYPP_FOR_RE),
                       parser_re(FYPP_BLOCK_RE), parser_re(FYPP_CALL_RE), parser_re(FYPP_MUTE_RE)]
PREPRO_CONTINUE_SCOPE = [None, parser_re(FYPP_ELIF_ELSE_RE), None, None, None, None]
PREPRO_END_SCOPE = [parser_re(FYPP_ENDDEF_RE), parser_re(FYPP_ENDIF_RE), parser_re(FYPP_ENDFOR_RE),
                       parser_re(FYPP_ENDBLOCK_RE), parser_re(FYPP_ENDCALL_RE),
                       parser_re(FYPP_ENDMUTE_RE)]

class plusminus_parser(parser_re):
    """parser for +/- in addition
    """
    def __init__(self, regex):
        self._re = regex
        self._re_excl = re.compile(r"\b(\d+\.?\d*|\d*\.?\d+)[de]" + EOL_STR, RE_FLAGS)

    def split(self, line):
        partsplit = self._re.split(line)
        partsplit_out = []

        # exclude splits due to '+/-' in real literals
        for n, part in enumerate(partsplit):
            if re.search(r"^(\+|-)$", part):
                if self._re_excl.search(partsplit[n-1]):
                    if n==1: partsplit_out = [partsplit[n-1]]
                    if n + 1 >= len(partsplit) or not partsplit_out:
                        raise FprettifyParseException("non-standard expression involving + or -",'',0)
                    partsplit_out[-1] += part + partsplit[n+1]
                else:
                    if n==1: partsplit_out = [partsplit[n-1]]
                    if n + 1 >= len(partsplit):
                        raise FprettifyParseException("non-standard expression involving + or -",'',0)
                    partsplit_out += [part, partsplit[n+1]]

        if not partsplit_out: partsplit_out = partsplit

        return partsplit_out

# two-sided operators
LR_OPS_RE = [REL_OP_RE, LOG_OP_RE, plusminus_parser(PLUSMINUS_RE), MULTDIV_RE, PRINT_RE]

USE_RE = re.compile(
    SOL_STR + "USE(\s+|(,.+?)?::\s*)\w+?((,.+?=>.+?)+|,\s*only\s*:.+?)?$" + EOL_STR, RE_FLAGS)

# markups to deactivate formatter
NO_ALIGN_RE = re.compile(SOL_STR + r"&\s*[^\s*]+")

class where_parser(parser_re):
    """parser for where / forall construct
    """
    def search(self, line):
        match = self._re.search(line)

        if match:
            level = 0
            for pos, char in CharFilter(line):
                [what_del_open, what_del_close] = get_curr_delim(line, pos)

                if what_del_open:
                    if what_del_open.group() == r'(': level += 1

                if what_del_close and what_del_close.group() == r')':
                    if level == 1:
                        if EMPTY_RE.search(line[pos+1:]):
                            return True
                        else:
                            return False
                    else:
                        level += -1

        return False

forall_parser = where_parser

def build_scope_parser(fypp=True, mod=True):
    parser = {}
    parser['new'] = \
        [parser_re(IF_RE), parser_re(DO_RE), parser_re(SELCASE_RE), parser_re(SUBR_RE),
         parser_re(FCT_RE),
         parser_re(INTERFACE_RE), parser_re(TYPE_RE), parser_re(ENUM_RE), parser_re(ASSOCIATE_RE),
         None, parser_re(BLK_RE), where_parser(WHERE_RE), forall_parser(FORALL_RE)]

    parser['continue'] = \
        [parser_re(ELSE_RE), None, parser_re(CASE_RE), parser_re(CONTAINS_RE),
         parser_re(CONTAINS_RE),
         None, parser_re(CONTAINS_RE), None, None,
         None, None, parser_re(ELSEWHERE_RE), None]

    parser['end'] = \
        [parser_re(ENDIF_RE), parser_re(ENDDO_RE), parser_re(ENDSEL_RE), parser_re(ENDSUBR_RE),
         parser_re(ENDFCT_RE),
         parser_re(ENDINTERFACE_RE), parser_re(ENDTYPE_RE), parser_re(ENDENUM_RE), parser_re(ENDASSOCIATE_RE),
         parser_re(ENDANY_RE,spec=False), parser_re(ENDBLK_RE), parser_re(ENDWHERE_RE), parser_re(ENDFORALL_RE)]

    if mod:
        parser['new'].extend([parser_re(MOD_RE), parser_re(SMOD_RE), parser_re(PROG_RE)])
        parser['continue'].extend([parser_re(CONTAINS_RE), parser_re(CONTAINS_RE), parser_re(CONTAINS_RE)])
        parser['end'].extend([parser_re(ENDMOD_RE), parser_re(ENDSMOD_RE), parser_re(ENDPROG_RE)])

    if fypp:
        parser['new'].extend(PREPRO_NEW_SCOPE)
        parser['continue'].extend(PREPRO_CONTINUE_SCOPE)
        parser['end'].extend(PREPRO_END_SCOPE)

    return parser

# match namelist names
NML_RE = re.compile(r"(/\w+/)", RE_FLAGS)
# find namelists and data statements
NML_STMT_RE = re.compile(SOL_STR + r"NAMELIST.*/.*/", RE_FLAGS)
DATA_STMT_RE = re.compile(SOL_STR + r"DATA\s+\w", RE_FLAGS)

## Regexp for f90 keywords'
F90_KEYWORDS_RE = re.compile(r"\b(" + "|".join((
    "allocatable", "allocate", "assign", "assignment", "backspace",
    "block", "call", "case", "character", "close", "common", "complex",
    "contains", "continue", "cycle", "data", "deallocate",
    "dimension", "do", "double", "else", "elseif", "elsewhere", "end",
    "enddo", "endfile", "endif", "entry", "equivalence", "exit",
    "external", "forall", "format", "function", "goto", "if",
    "implicit", "include", "inquire", "integer", "intent",
    "interface", "intrinsic", "logical", "module", "namelist", "none",
    "nullify", "only", "open", "operator", "optional", "parameter",
    "pause", "pointer", "precision", "print", "private", "procedure",
    "program", "public", "read", "real", "recursive", "result", "return",
    "rewind", "save", "select", "sequence", "stop", "subroutine",
    "target", "then", "type", "use", "where", "while", "write",
    ## F95 keywords.
    "elemental", "pure",
    ## F2003
    "abstract", "associate", "asynchronous", "bind", "class",
    "deferred", "enum", "enumerator", "extends", "extends_type_of",
    "final", "generic", "import", "non_intrinsic", "non_overridable",
    "nopass", "pass", "protected", "same_type_as", "value", "volatile",
    ## F2008.
    "contiguous", "submodule", "concurrent", "codimension",
    "sync all", "sync memory", "critical", "image_index",
    )) + r")\b", RE_FLAGS)

## Regexp whose first part matches F90 intrinsic procedures.
## Add a parenthesis to avoid catching non-procedures.
F90_PROCEDURES_RE = re.compile(r"\b(" + "|".join((
    "abs", "achar", "acos", "adjustl", "adjustr", "aimag", "aint",
    "all", "allocated", "anint", "any", "asin", "associated",
    "atan", "atan2", "bit_size", "btest", "ceiling", "char", "cmplx",
    "conjg", "cos", "cosh", "count", "cshift", "date_and_time", "dble",
    "digits", "dim", "dot_product", "dprod", "eoshift", "epsilon",
    "exp", "exponent", "floor", "fraction", "huge", "iachar", "iand",
    "ibclr", "ibits", "ibset", "ichar", "ieor", "index", "int", "ior",
    "ishft", "ishftc", "kind", "lbound", "len", "len_trim", "lge", "lgt",
    "lle", "llt", "log", "log10", "logical", "matmul", "max",
    "maxexponent", "maxloc", "maxval", "merge", "min", "minexponent",
    "minloc", "minval", "mod", "modulo", "mvbits", "nearest", "nint",
    "not", "pack", "precision", "present", "product", "radix",
    ## Real is taken out here to avoid highlighting declarations.
    "random_number", "random_seed", "range", ## "real"
    "repeat", "reshape", "rrspacing", "scale", "scan",
    "selected_int_kind", "selected_real_kind", "set_exponent",
    "shape", "sign", "sin", "sinh", "size", "spacing", "spread", "sqrt",
    "sum", "system_clock", "tan", "tanh", "tiny", "transfer",
    "transpose", "trim", "ubound", "unpack", "verify",
    ## F95 intrinsic functions.
    "null", "cpu_time",
    ## F2003.
    "move_alloc", "command_argument_count", "get_command",
    "get_command_argument", "get_environment_variable",
    "selected_char_kind", "wait", "flush", "new_line",
    "extends", "extends_type_of", "same_type_as", "bind",
    ## F2003 ieee_arithmetic intrinsic module.
    "ieee_support_underflow_control", "ieee_get_underflow_mode",
    "ieee_set_underflow_mode",
    ## F2003 iso_c_binding intrinsic module.
    "c_loc", "c_funloc", "c_associated", "c_f_pointer",
    "c_f_procpointer",
    ## F2008.
    "bge", "bgt", "ble", "blt", "dshiftl", "dshiftr", "leadz", "popcnt",
    "poppar", "trailz", "maskl", "maskr", "shifta", "shiftl", "shiftr",
    "merge_bits", "iall", "iany", "iparity", "storage_size",
    "bessel_j0", "bessel_j1", "bessel_jn",
    "bessel_y0", "bessel_y1", "bessel_yn",
    "erf", "erfc", "erfc_scaled", "gamma", "hypot", "log_gamma",
    "norm2", "parity", "findloc", "is_contiguous",
    "sync images", "lock", "unlock", "image_index",
    "lcobound", "ucobound", "num_images", "this_image",
    ## F2008 iso_fortran_env module.
    "compiler_options", "compiler_version",
    ## F2008 iso_c_binding module.
    "c_sizeof"

    )) + r")\b", RE_FLAGS)

F90_MODULES_RE = re.compile(r"\b(" + "|".join((
    ## F2003/F2008 module names
    "iso_fortran_env",
    "iso_c_binding",
    "ieee_exceptions",
    "ieee_arithmetic",
    "ieee_features"
    )) + r")\b", RE_FLAGS)

## Regexp matching intrinsic operators
F90_OPERATORS_RE = re.compile(r"(" + "|".join([r"\." + a + r"\." for a in (
    "and", "eq", "eqv", "false", "ge", "gt", "le", "lt", "ne",
    "neqv", "not", "or", "true"
    )]) + r")", RE_FLAGS)

## Regexp for Fortran intrinsic constants
F90_CONSTANTS_RE = re.compile(r"\b(" + "|".join((
    ## F2003 iso_fortran_env constants.
    "input_unit", "output_unit", "error_unit",
    "iostat_end", "iostat_eor",
    "numeric_storage_size", "character_storage_size",
    "file_storage_size",
    ## F2003 iso_c_binding constants.
    "c_int", "c_short", "c_long", "c_long_long", "c_signed_char",
    "c_size_t",
    "c_int8_t", "c_int16_t", "c_int32_t", "c_int64_t",
    "c_int_least8_t", "c_int_least16_t", "c_int_least32_t",
    "c_int_least64_t",
    "c_int_fast8_t", "c_int_fast16_t", "c_int_fast32_t",
    "c_int_fast64_t",
    "c_intmax_t", "c_intptr_t",
    "c_float", "c_double", "c_long_double",
    "c_float_complex", "c_double_complex", "c_long_double_complex",
    "c_bool", "c_char",
    "c_null_char", "c_alert", "c_backspace", "c_form_feed",
    "c_new_line", "c_carriage_return", "c_horizontal_tab",
    "c_vertical_tab",
    "c_ptr", "c_funptr", "c_null_ptr", "c_null_funptr",
    ## F2008 iso_fortran_env constants.
    "character_kinds", "int8", "int16", "int32", "int64",
    "integer_kinds", "iostat_inquire_internal_unit",
    "logical_kinds", "real_kinds", "real32", "real64", "real128",
    "lock_type", "atomic_int_kind", "atomic_logical_kind",
    )) + r")\b", RE_FLAGS)

F90_INT_RE = r"[-+]?[0-9]+"
F90_FLOAT_RE = r"[-+]?([0-9]+\.[0-9]*|\.[0-9]+)"
F90_NUMBER_RE = "(" + F90_INT_RE + "|" + F90_FLOAT_RE + ")"
F90_FLOAT_EXP_RE = F90_NUMBER_RE + r"[eEdD]" + F90_NUMBER_RE
F90_NUMBER_ALL_RE = "(" + F90_NUMBER_RE + "|" + F90_FLOAT_EXP_RE + ")"
F90_NUMBER_ALL_REC = re.compile(F90_NUMBER_ALL_RE, RE_FLAGS)

## F90_CONSTANTS_TYPES_RE = re.compile(r"\b" + F90_NUMBER_ALL_RE + "_(" + "|".join([a + r"\b" for a in (
F90_CONSTANTS_TYPES_RE = re.compile(
    r"(" + F90_NUMBER_ALL_RE + ")*_(" + "|".join((
    ## F2003 iso_fortran_env constants.
    ## F2003 iso_c_binding constants.
    "c_int", "c_short", "c_long", "c_long_long", "c_signed_char",
    "c_size_t",
    "c_int8_t", "c_int16_t", "c_int32_t", "c_int64_t",
    "c_int_least8_t", "c_int_least16_t", "c_int_least32_t",
    "c_int_least64_t",
    "c_int_fast8_t", "c_int_fast16_t", "c_int_fast32_t",
    "c_int_fast64_t",
    "c_intmax_t", "c_intptr_t",
    "c_float", "c_double", "c_long_double",
    "c_float_complex", "c_double_complex", "c_long_double_complex",
    "c_bool", "c_char",
    ## F2008 iso_fortran_env constants.
    "character_kinds", "int8", "int16", "int32", "int64",
    "integer_kinds",
    "logical_kinds", "real_kinds", "real32", "real64", "real128",
    "lock_type", "atomic_int_kind", "atomic_logical_kind",
    )) + r")\b", RE_FLAGS)


class F90Indenter(object):
    """
    Parses encapsulation of subunits / scopes line by line
    and updates the indentation.
    """

    def __init__(self, scope_parser, first_indent, rel_indent, filename):
        # scopes / subunits:
        self._scope_storage = []
        # indents for all fortran lines:
        self._indent_storage = []
        # indents of actual lines of current fortran line
        self._line_indents = []

        self._parser = scope_parser

        self._filename = filename
        self._aligner = F90Aligner(filename)

        # no lines have been processed yet:
        self._initial = True

        # implicit scopes: we define implicit scopes, as many as match
        # first_indent and rel_indent. This allows for, e.g., a properly
        # indented "END FUNCTION" without matching "FUNCTION" statement:
        if rel_indent > 0:
            for n_impl in range(first_indent % rel_indent, first_indent + 1, rel_indent):
                self._indent_storage += [n_impl]

        if not self._indent_storage:
            self._indent_storage = [0]

    def process_lines_of_fline(self, f_line, lines, rel_ind, rel_ind_con,
                               line_nr, indent_fypp=True, manual_lines_indent=None):
        """
        Process all lines that belong to a Fortran line `f_line`.

        Impose a relative indent of `rel_ind` for current Fortran line,
        and `rel_ind_con` for line continuation.
        By default line continuations are auto-aligned by F90Aligner
        :param f_line: fortran line
        :param lines: actual lines belonging to f_line
        :param rel_ind: relative scope indent size for this line
        :rel_ind_con: relative continuation indent size for this line
        :line_nr: line number
        :indent_fypp: whether or not to include fypp preprocessor lines
        :manual_lines_indent: don't use F90Aligner but manually impose
                              indents for continuations
        """

        if (self._initial and
            (PROG_RE.match(f_line) or MOD_RE.match(f_line))):
            self._indent_storage[-1] = 0

        self._line_indents = [0] * len(lines)
        br_indent_list = [0] * len(lines)

        # local variables to avoid self hassle:
        line_indents = self._line_indents

        scopes = self._scope_storage
        indents = self._indent_storage
        filename = self._filename

        # check statements that start new scope
        is_new = False
        valid_new = False

        f_filter = CharFilter(f_line, filter_fypp=not indent_fypp)
        f_line_filtered = f_filter.filter_all()

        for new_n, newre in enumerate(self._parser['new']):
            if newre and newre.search(f_line_filtered) and \
                not self._parser['end'][new_n].search(f_line_filtered):
                what_new = new_n
                is_new = True
                valid_new = True
                scopes.append(what_new)
                log_message("{}: {}".format(what_new, f_line),
                            "debug", filename, line_nr)

        # check statements that continue scope
        is_con = False
        valid_con = False
        for con_n, conre in enumerate(self._parser['continue']):
            if conre and conre.search(f_line_filtered):
                what_con = con_n
                is_con = True
                log_message("{}: {}".format(
                    what_con, f_line), "debug", filename, line_nr)
                if len(scopes) > 0:
                    what = scopes[-1]
                    if what == what_con:
                        valid_con = True

        # check statements that end scope
        is_end = False
        valid_end = False
        for end_n, endre in enumerate(self._parser['end']):
            if endre and endre.search(f_line_filtered):
                what_end = end_n
                is_end = True
                log_message("{}: {}".format(
                    what_end, f_line), "debug", filename, line_nr)
                if len(scopes) > 0:
                    what = scopes.pop()
                    if what == what_end or not self._parser['end'][what_end].spec:
                        valid_end = True
                        log_message("{}: {}".format(
                            what_end, f_line), "debug", filename, line_nr)
                else:
                    valid_end = True

        # fypp preprocessor scopes may be within continuation lines
        if indent_fypp and len(lines) > 1 and not FYPP_LINE_RE.search(f_line_filtered):

            for new_n, newre in enumerate(PREPRO_NEW_SCOPE):
                for l in lines:
                    if(newre and newre.search(l)):
                        is_new = True
                        valid_new = True
                        scopes.append(new_n)

            for end_n, endre in enumerate(PREPRO_END_SCOPE):
                for l in lines:
                    if(endre and endre.search(l)):
                        is_end = True
                        valid_end = True
                        if len(scopes) > 0:
                            what = scopes.pop()

        # deal with line breaks
        if not manual_lines_indent:
            self._aligner.process_lines_of_fline(
                f_line, lines, rel_ind_con, line_nr)
            br_indent_list = self._aligner.get_lines_indent()
        else:
            br_indent_list = manual_lines_indent

        for pos in range(0, len(lines) - 1):
            line_indents[pos + 1] = br_indent_list[pos + 1]

        if is_new and not is_end:
            if not valid_new:
                log_message('invalid scope opening statement',
                            "info", filename, line_nr)

            line_indents = [ind + indents[-1] for ind in line_indents]

            indents.append(rel_ind + indents[-1])

        elif (not is_new) and (is_con or is_end):
            valid = valid_con if is_con else valid_end

            if not valid:
                line_indents = [ind + indents[-1] for ind in line_indents]
                log_message('invalid scope closing statement',
                            "info", filename, line_nr)
            else:
                if len(indents) > 1 or self._initial:
                    line_indents = [ind + indents[-2 + self._initial]
                                    for ind in line_indents]

            if is_end and valid:
                if len(indents) > 1:
                    indents.pop()
                else:
                    indents[-1] = 0

        else:
            line_indents = [ind + indents[-1] for ind in line_indents]

        # we have processed first line:
        self._initial = False

        # reassigning self.* to the updated variables
        self._line_indents = line_indents
        self._scope_storage = scopes
        self._indent_storage = indents

    def get_fline_indent(self):
        """after processing, retrieve the indentation of the full Fortran line."""
        return self._indent_storage[-1]

    def get_lines_indent(self):
        """after processing, retrieve the indents of all line parts."""
        return self._line_indents


class F90Aligner(object):
    """
    Alignment of continuations of a broken line,
    based on the following heuristics:

    if line break in brackets
        We are parsing the level of nesting
        and align to most inner bracket delimiter.

    else if line is an assignment
        alignment to '=' or '=>'.
        note: assignment operator recognized as any '=' that is not
        part of another operator and that is not enclosed in bracket

    else if line is a declaration
        alignment to '::'

    else
        default indent
    """

    def __init__(self, filename):
        self._filename = filename
        self.__init_line(0)

    def __init_line(self, line_nr):
        """initialization before processing new line"""
        self._line_nr = line_nr
        self._line_indents = [0]
        self._level = 0
        self._br_indent_list = [0]

    def process_lines_of_fline(self, f_line, lines, rel_ind, line_nr):
        """
        process all lines that belong to a Fortran line `f_line`,
        `rel_ind` is the relative indentation size.
        """

        self.__init_line(line_nr)

        is_decl = VAR_DECL_RE.search(f_line) or PUBLIC_RE.search(f_line) or PRIVATE_RE.match(f_line)
        is_use = USE_RE.search(f_line)
        for pos, line in enumerate(lines):
            self.__align_line_continuations(
                line, is_decl, is_use, rel_ind, self._line_nr + pos)
            if pos + 1 < len(lines):
                self._line_indents.append(self._br_indent_list[-1])

        if len(self._br_indent_list) > 2 or self._level:
            log_message('unpaired bracket delimiters',
                        "info", self._filename, self._line_nr)

    def get_lines_indent(self):
        """after processing, retrieve the indents of all line parts."""
        return self._line_indents

    def __align_line_continuations(self, line, is_decl, is_use, indent_size, line_nr):
        """align continuation lines."""

        indent_list = self._br_indent_list
        level = self._level
        filename = self._filename

        pos_eq = 0
        pos_ldelim = []
        pos_rdelim = []
        ldelim = []
        rdelim = []

        # find delimiters that are not ended on this line.
        # find proper alignment to most inner delimiter
        # or alignment to assignment operator
        rel_ind = indent_list[-1]  # indentation of prev. line

        end_of_delim = -1

        for pos, char in CharFilter(line):

            what_del_open = None
            what_del_close = None
            if pos > end_of_delim:
                [what_del_open, what_del_close] = get_curr_delim(line, pos)

            if what_del_open:
                what_del_open = what_del_open.group()
                end_of_delim = pos + len(what_del_open) - 1
                level += 1
                indent_list.append(pos + len(what_del_open) + rel_ind)
                pos_ldelim.append(pos)
                ldelim.append(what_del_open)
            if what_del_close:
                what_del_close = what_del_close.group()
                end_of_delim = pos + len(what_del_close) - 1
                if level > 0:
                    level += -1
                    indent_list.pop()
                else:
                    log_message('unpaired bracket delimiters',
                                "info", filename, line_nr)

                if pos_ldelim:
                    pos_ldelim.pop()
                    what_del_open = ldelim.pop()
                    valid = False
                    if what_del_open == r"(":
                        valid = what_del_close == r")"
                    if what_del_open == r"(/":
                        valid = what_del_close == r"/)"
                    if what_del_open == r"[":
                        valid = what_del_close == r"]"
                    if not valid:
                        log_message('unpaired bracket delimiters',
                                    "info", filename, line_nr)

                else:
                    pos_rdelim.append(pos)
                    rdelim.append(what_del_close)
            if char == ',' and not level and pos_eq > 0:
                # a top level comma removes previous alignment position.
                # (see issue #11)
                pos_eq = 0
                indent_list.pop()
            if not level and not is_decl and char == '=' and not REL_OP_RE.search(
                    line[max(0, pos - 1):min(pos + 2, len(line))]):
                        # should only have one assignment per line!
                if pos_eq > 0:
                    raise FprettifyInternalException(
                            "found more than one assignment in the same Fortran line", filename, line_nr)
                is_pointer = line[pos + 1] == '>'
                pos_eq = pos + 1
                # don't align if assignment operator directly before
                # line break
                if not re.search(r"=>?\s*" + LINEBREAK_STR, line,
                                 RE_FLAGS):
                    indent_list.append(
                        pos_eq + 1 + is_pointer + indent_list[-1])
            elif is_decl and line[pos:pos + 2] == '::' and not re.search(r"::\s*" + LINEBREAK_STR, line, RE_FLAGS):
                indent_list.append(pos + 3 + indent_list[-1])
            elif is_use and line[pos] == ':' and not re.search(r":\s*" + LINEBREAK_STR, line, RE_FLAGS):
                indent_list.append(pos + 2 + indent_list[-1])

        # Don't align if delimiter opening directly before line break
        if level and re.search(DEL_OPEN_STR + r"\s*" + LINEBREAK_STR, line,
                               RE_FLAGS):
            if len(indent_list) > 1:
                indent_list[-1] = indent_list[-2]
            else:
                indent_list[-1] = 0

        if not indent_list[-1]:
            indent_list[-1] = indent_size

        self._level = level


def inspect_ffile_format(infile, indent_size, strict_indent, indent_fypp=False, orig_filename=None):
    """
    Determine indentation by inspecting original Fortran file.

    This is mainly for finding aligned blocks of DO/IF statements.
    Also check if it has f77 constructs.
    :param infile: open file
    :param indent_size: the default indent size
    :orig_filename: filename used for messages
    :returns: [ target indent sizes for each line,
                indent of first line (offset) ]
    """
    if not orig_filename:
        orig_filename = infile.name

    num_labels = False
    indents = []
    stream = InputStream(infile, filter_fypp=not indent_fypp, orig_filename=orig_filename)
    prev_offset = 0
    first_indent = -1

    while 1:
        f_line, _, lines = stream.next_fortran_line()
        if not lines:
            break

        f_line, lines, label = preprocess_labels(f_line, lines)

        offset = len(lines[0]) - len(lines[0].lstrip(' '))
        if f_line.strip() and first_indent == -1:
            first_indent = offset
        indents.append(offset - prev_offset)

        # don't impose indentation for blocked do/if constructs:
        if (IF_RE.search(f_line) or DO_RE.search(f_line)):
            if (prev_offset != offset or strict_indent):
                indents[-1] = indent_size
        else:
            indents[-1] = indent_size

        prev_offset = offset

    return indents, first_indent


def replace_relational_single_fline(f_line, cstyle):
    """
    format a single Fortran line - replaces scalar relational
    operators in logical expressions to either Fortran or C-style.
    .lt.  <-->  <
    .le.  <-->  <=
    .gt.  <-->  >
    .ge.  <-->  >=
    .eq.  <-->  ==
    .ne.  <-->  /=
    """

    new_line = f_line

    # only act on lines that do contain a relation
    if REL_OP_RE.search(f_line):
        # check that relation is not inside quotes, a string, or commented
        # (think of underlining a heading with === or things like markup being printed which we do not replace)
        pos_prev = -1
        pos = -1
        line_parts = ['']
        for pos, char in CharFilter(f_line):
            if pos > pos_prev + 1: # skipped string
                line_parts.append(f_line[pos_prev + 1:pos].strip()) # append string
                line_parts.append('')

            line_parts[-1] += char

            pos_prev = pos

        if pos + 1 < len(f_line):
            line_parts.append(f_line[pos + 1:])

        for pos, part in enumerate(line_parts):
            # exclude comments, strings:
            if not STR_OPEN_RE.match(part):
                # also exclude / if we see a namelist and data statement
                if cstyle:
                    part = re.sub(r"\.LT\.", "<   ", part, flags=RE_FLAGS)
                    part = re.sub(r"\.LE\.", "<=  ", part, flags=RE_FLAGS)
                    part = re.sub(r"\.GT\.", ">   ", part, flags=RE_FLAGS)
                    part = re.sub(r"\.GE\.", ">=  ", part, flags=RE_FLAGS)
                    part = re.sub(r"\.EQ\.", "==  ", part, flags=RE_FLAGS)
                    part = re.sub(r"\.NE\.", "/=  ", part, flags=RE_FLAGS)
                else:
                    part = re.sub(r"<=",  ".le.", part, flags=RE_FLAGS)
                    part = re.sub(r"<",   ".lt.", part, flags=RE_FLAGS)
                    part = re.sub(r">=",  ".ge.", part, flags=RE_FLAGS)
                    part = re.sub(r">",   ".gt.", part, flags=RE_FLAGS)
                    part = re.sub(r"==",  ".eq.", part, flags=RE_FLAGS)
                    part = re.sub(r"\/=", ".ne.", part, flags=RE_FLAGS)

            line_parts[pos] = part

        new_line = ''.join(line_parts)

    return new_line


def replace_keywords_single_fline(f_line, case_dict):
    """
    format a single Fortran line - change case of keywords
    """

    new_line = f_line

    # Collect words list
    pos_prev = -1
    pos = -1
    line_parts = ['']
    for pos, char in CharFilter(f_line):
        if pos > pos_prev + 1: # skipped string
            line_parts.append(f_line[pos_prev + 1:pos].strip()) # append string
            line_parts.append('')

        line_parts[-1] += char

        pos_prev = pos

    if pos + 1 < len(f_line):
        line_parts.append(f_line[pos + 1:])

    line_parts = [[a] if STR_OPEN_RE.match(a) else re.split(F90_OPERATORS_RE,a)
                  for a in line_parts]  # problem, split "."
    line_parts = [b for a in line_parts for b in a]

    ## line_parts = [[a] if STR_OPEN_RE.match(a) else re.split('(\W)',a)
    ##               for a in line_parts]  # problem, split "."
    line_parts = [[a] if STR_OPEN_RE.match(a)
                  else re.split('([^a-zA-Z0-9_.])',a)
                  for a in line_parts]
    line_parts = [b for a in line_parts for b in a]

    swapcase = lambda s, a: s if a==0 else (s.lower() if a==1 else s.upper())

    nbparts = len(line_parts)
    for pos, part in enumerate(line_parts):
        # exclude comments, strings:
        if part.strip() and not STR_OPEN_RE.match(part):
            if F90_KEYWORDS_RE.match(part):
                part = swapcase(part, case_dict['keywords'])
            elif F90_MODULES_RE.match(part):
                part = swapcase(part, case_dict['procedures'])
            elif F90_PROCEDURES_RE.match(part):
                ok = False
                for pos2 in range(pos+1, nbparts):
                    part2 = line_parts[pos2]
                    if part2.strip() and not (part2 == '\n' or STR_OPEN_RE.match(part2)):
                        ok = (part2 == '(')
                        break
                if ok:
                    part = swapcase(part, case_dict['procedures'])
            elif F90_OPERATORS_RE.match(part):
                part = swapcase(part, case_dict['operators'])
            elif F90_CONSTANTS_RE.match(part):
                part = swapcase(part, case_dict['constants'])
            elif F90_CONSTANTS_TYPES_RE.match(part):
                part = swapcase(part, case_dict['constants'])
            elif F90_NUMBER_ALL_REC.match(part):
                part = swapcase(part, case_dict['constants'])

            line_parts[pos] = part

    new_line = ''.join(line_parts)

    return new_line


def format_single_fline(f_line, whitespace, whitespace_dict, linebreak_pos,
                        ampersand_sep, scope_parser, format_decl, filename, line_nr,
                        auto_format=True):
    """
    format a single Fortran line - imposes white space formatting
    and inserts linebreaks.
    Takes a logical Fortran line `f_line` as input as well as the positions
    of the linebreaks (`linebreak_pos`), and the number of
    separating whitespace characters before ampersand (`ampersand_sep`).
    `filename` and `line_nr` just for error messages.
    The higher `whitespace`, the more white space characters inserted -
    whitespace = 0, 1, 2, 3 are currently supported.
    whitespace formatting can additionally controlled more fine-grained
    via a dictionary of bools (whitespace_dict)
    auto formatting can be turned off by setting `auto_format` to False.
    """

    # define whether to put whitespaces around operators:
    mapping = {
            'comma': 0,           # 0: comma, semicolon
            'assignments': 1,     # 1: assignment operators
            'relational': 2,      # 2: relational operators
            'logical': 3,         # 3: logical operators
            'plusminus': 4,       # 4: arithm. operators plus and minus
            'multdiv': 5,         # 5: arithm. operators multiply and divide
            'print': 6,           # 6: print / read statements
            'type': 7,            # 7: select type components
            'intrinsics': 8,      # 8: intrinsics
            'decl': 9             # 9: declarations
            }

    if whitespace == 0:
        spacey = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    elif whitespace == 1:
        spacey = [1, 1, 1, 1, 0, 0, 1, 0, 1, 1]
    elif whitespace == 2:
        spacey = [1, 1, 1, 1, 1, 0, 1, 0, 1, 1]
    elif whitespace == 3:
        spacey = [1, 1, 1, 1, 1, 1, 1, 0, 1, 1]
    elif whitespace == 4:
        spacey = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    else:
        raise NotImplementedError("unknown value for whitespace")

    if whitespace_dict:
        # iterate over dictionary and override settings for 'spacey'
        for key, value in mapping.items():
            if whitespace_dict[key] == True:
                spacey[value] = 1
            elif whitespace_dict[key] == False:
                spacey[value] = 0

    line = f_line
    line_orig = line

    if auto_format:

        line = rm_extra_whitespace(line, format_decl)
        line = add_whitespace_charwise(line, spacey, scope_parser, format_decl, filename, line_nr)
        line = add_whitespace_context(line, spacey)

    lines_out = split_reformatted_line(
        line_orig, linebreak_pos, ampersand_sep, line, filename, line_nr)
    return lines_out


def rm_extra_whitespace(line, format_decl):
    """rm all unneeded whitespace chars, except for declarations"""
    line_ftd = ''
    pos_prev = -1
    pos = -1
    for pos, char in CharFilter(line):
        if format_decl:
            is_decl = False
        else:
            is_decl = line[pos:].lstrip().startswith('::') or line[
                :pos].rstrip().endswith('::')

        if pos > pos_prev + 1: # skipped string
            line_ftd = line_ftd + line[pos_prev + 1:pos]

        if char == ' ':
            # remove double spaces:
            if line_ftd and (re.search(r'[\w]', line_ftd[-1]) or is_decl):
                line_ftd = line_ftd + char
        else:
            if (line_ftd and line_ftd[-1] == ' ' and
                    (not re.search(r'[\w]', char) and not is_decl)):
                line_ftd = line_ftd[:-1]  # remove spaces except between words
            line_ftd = line_ftd + char
        pos_prev = pos

    line_ftd = line_ftd + line[pos+1:]
    return line_ftd


def add_whitespace_charwise(line, spacey, scope_parser, format_decl, filename, line_nr):
    """add whitespace character wise (no need for context aware parsing)"""
    line_ftd = line
    pos_eq = []
    end_of_delim = -1
    level = 0
    for pos, char in CharFilter(line):
        # offset w.r.t. unformatted line
        offset = len(line_ftd) - len(line)

        # format delimiters
        what_del_open = None
        what_del_close = None
        if pos > end_of_delim:
            [what_del_open, what_del_close] = get_curr_delim(line, pos)

        if what_del_open or what_del_close:
            sep1 = 0
            sep2 = 0

            if what_del_open:
                delim = what_del_open.group()
            else:
                delim = what_del_close.group()

            lhs = line_ftd[:pos + offset]
            rhs = line_ftd[pos + len(delim) + offset:]

            # format opening delimiters
            if what_del_open:
                level += 1  # new scope
                # add separating whitespace before opening delimiter
                # with some exceptions:
                # FIXME: duplication of regex, better to include them into
                # INTR_STMTS_PAR
                if ((not re.search((r"(" + DEL_OPEN_STR +
                                    r"|[\w\*/=\+\-:])\s*$"),
                                   line[:pos], RE_FLAGS) and
                     not EMPTY_RE.search(line[:pos])) or
                        re.search(SOL_STR + r"(\w+\s*:)?(ELSE)?\s*IF\s*$",
                                  line[:pos], RE_FLAGS) or
                        re.search(SOL_STR + r"(\w+\s*:)?\s*DO\s+WHILE\s*$",
                                  line[:pos], RE_FLAGS) or
                        re.search(SOL_STR + r"(SELECT)?\s*CASE\s*$",
                                  line[:pos], RE_FLAGS) or
                        re.search(SOL_STR + r"(SELECT)?\s*RANK\s*$",
                                  line[:pos], RE_FLAGS) or
                        re.search(SOL_STR + r"SELECT\s*TYPE\s*$",
                                  line[:pos], RE_FLAGS) or
                        re.search(SOL_STR + r"CLASS\s*DEFAULT\s*$",
                                  line[:pos], RE_FLAGS) or
                        re.search(SOL_STR + r"(TYPE|CLASS)\s+IS\s*$",
                                  line[:pos], RE_FLAGS) or
                        re.search(r"(?<!%)\b" + INTR_STMTS_PAR + r"\s*$",
                                  line[:pos], RE_FLAGS)):
                    sep1 = 1 * spacey[8]

            # format closing delimiters
            else:
                if level > 0:
                    level += -1  # close scope
                else:
                    log_message('unpaired bracket delimiters',
                                "info", filename, line_nr)

                # add separating whitespace after closing delimiter
                # with some exceptions:
                if not re.search(r"^\s*(" + DEL_CLOSE_STR + r"|[,%:/\*])",
                                 line[pos + 1:], RE_FLAGS):
                    sep2 = 1
                elif re.search(r"^\s*::", line[pos + 1:], RE_FLAGS):
                    sep2 = len(rhs) - len(rhs.lstrip(' ')) if not format_decl else 1

            # where delimiter token ends
            end_of_delim = pos + len(delim) - 1

            line_ftd = lhs.rstrip(' ') + ' ' * sep1 + \
                delim + ' ' * sep2 + rhs.lstrip(' ')

        # format commas and semicolons
        if char in [',', ';']:
            lhs = line_ftd[:pos + offset]
            rhs = line_ftd[pos + 1 + offset:]
            line_ftd = lhs.rstrip(' ') + char + ' ' * \
                spacey[0] + rhs.lstrip(' ')
            line_ftd = line_ftd.rstrip(' ')

        # format type selector %
        if char == "%":
            lhs = line_ftd[:pos + offset]
            rhs = line_ftd[pos + 1 + offset:]
            line_ftd = lhs.rstrip(' ') \
                    + ' ' * spacey[7] \
                    + char \
                    + ' ' * spacey[7] \
                    + rhs.lstrip(' ')
            line_ftd = line_ftd.rstrip(' ')

        # format '::'
        if format_decl and line[pos:pos+2] == "::":
            lhs = line_ftd[:pos + offset]
            rhs = line_ftd[pos + 2 + offset:]
            line_ftd = lhs.rstrip(' ') \
                    + ' ' * spacey[9] \
                    + '::' + ' ' * spacey[9] \
                    + rhs.lstrip(' ')
            line_ftd = line_ftd.rstrip(' ')

        # format .NOT.
        if re.search(r"^\.NOT\.", line[pos:pos + 5], RE_FLAGS):
            lhs = line_ftd[:pos + offset]
            rhs = line_ftd[pos + 5 + offset:]
            line_ftd = lhs.rstrip(
                ' ') + line[pos:pos + 5] + ' ' * spacey[3] + rhs.lstrip(' ')

        # strip whitespaces from '=' and prepare assignment operator
        # formatting:
        if char == '=' and not REL_OP_RE.search(line[pos - 1:pos + 2]):
            lhs = line_ftd[:pos + offset]
            rhs = line_ftd[pos + 1 + offset:]
            line_ftd = lhs.rstrip(' ') + '=' + rhs.lstrip(' ')
            is_pointer = line[pos + 1] == '>'
            if (not level) or is_pointer:  # remember position of assignment operator
                pos_eq.append(len(lhs.rstrip(' ')))

    line = line_ftd

    for pos in pos_eq:
        offset = len(line_ftd) - len(line)
        is_pointer = line[pos + 1] == '>'
        lhs = line_ftd[:pos + offset]
        rhs = line_ftd[pos + 1 + is_pointer + offset:]
        if is_pointer:
            assign_op = '=>'  # pointer assignment
        else:
            assign_op = '='  # assignment
        line_ftd = (lhs.rstrip(' ') +
                    ' ' * spacey[1] + assign_op +
                    ' ' * spacey[1] + rhs.lstrip(' '))
        # offset w.r.t. unformatted line

    is_end = False
    if END_RE.search(line_ftd):
        for endre in scope_parser['end']:
            if endre and endre.search(line_ftd):
                is_end = True
    if is_end:
        line_ftd = END_RE.sub(r'\1' + ' '*spacey[8] + r'\2', line_ftd)

    if level != 0:
        log_message('unpaired bracket delimiters', "info", filename, line_nr)

    return line_ftd


def add_whitespace_context(line, spacey):
    """
    for context aware whitespace formatting we extract line parts that are
    not comments or strings in order to be able to apply a context aware regex.
    """


    pos_prev = -1
    pos = -1
    line_parts = ['']
    for pos, char in CharFilter(line):
        if pos > pos_prev + 1: # skipped string
            line_parts.append(line[pos_prev + 1:pos].strip()) # append string
            line_parts.append('')

        line_parts[-1] += char

        pos_prev = pos

    if pos + 1 < len(line):
        line_parts.append(line[pos + 1:])

    # format namelists with spaces around /
    if NML_STMT_RE.match(line):
        for pos, part in enumerate(line_parts):
            # exclude comments, strings:
            if not STR_OPEN_RE.match(part):
                partsplit = NML_RE.split(part)
                line_parts[pos] = (' '.join(partsplit))

    # Two-sided operators
    for n_op, lr_re in enumerate(LR_OPS_RE):
        for pos, part in enumerate(line_parts):
            # exclude comments, strings:
            if not STR_OPEN_RE.match(part):
                # also exclude / if we see a namelist and data statement
                if not ( NML_STMT_RE.match(line) or DATA_STMT_RE.match(line) ):
                    partsplit = lr_re.split(part)
                    line_parts[pos] = (' ' * spacey[n_op + 2]).join(partsplit)

    line = ''.join(line_parts)

    for newre in [IF_RE, DO_RE, BLK_RE]:
        if newre.search(line) and re.search(SOL_STR + r"\w+\s*:", line):
            line = ': '.join(_.strip() for _ in line.split(':', 1))

    # format ':' for labels and use only statements
    if USE_RE.search(line):
        line = re.sub(r'(only)\s*:\s*', r'\g<1>:' + ' ' *
                      spacey[0], line, flags=RE_FLAGS)

    return line


def split_reformatted_line(line_orig, linebreak_pos_orig, ampersand_sep, line, filename, line_nr):
    """
    Infer linebreak positions of formatted line from linebreak positions in
    original line and split line.
    """
    # shift line break positions from original to reformatted line
    pos_new = 0
    pos_old = 0
    linebreak_pos_orig.sort(reverse=True)
    linebreak_pos_ftd = []
    while 1:

        if pos_new == len(line) or pos_old == len(line_orig):
            break

        if line[pos_new] != line_orig[pos_old]:
            raise FprettifyInternalException(
                "failed at finding line break position", filename, line_nr)

        if linebreak_pos_orig and pos_old > linebreak_pos_orig[-1]:
            linebreak_pos_orig.pop()
            linebreak_pos_ftd.append(pos_new)
            continue

        pos_new += 1
        while pos_new < len(line) and line[pos_new] == ' ':
            pos_new += 1

        pos_old += 1
        while pos_old < len(line_orig) and line_orig[pos_old] == ' ':
            pos_old += 1

    linebreak_pos_ftd.insert(0, 0)

    # We split line into parts and we insert ampersands at line end, but not
    # for empty lines and comment lines
    lines_split = [(line[l:r].rstrip(' ') +
                       ' ' * ampersand_sep[pos] + '&' * min(1, r - l))
                      for pos, (l, r) in enumerate(zip(linebreak_pos_ftd[0:-1],
                                                       linebreak_pos_ftd[1:]))]

    lines_split.append(line[linebreak_pos_ftd[-1]:])

    return lines_split


def diff(a, b, a_name, b_name):
    # type: (str, str, str, str) -> str

    """Return a unified diff string between strings `a` and `b`."""
    import difflib

    a_lines = [line + "\n" for line in a.splitlines()]
    b_lines = [line + "\n" for line in b.splitlines()]
    return "".join(
        difflib.unified_diff(a_lines, b_lines, fromfile=a_name, tofile=b_name, n=5)
    )

def reformat_inplace(filename, stdout=False, diffonly=False, **kwargs):  # pragma: no cover
    """reformat a file in place."""
    if filename == '-':
        infile = io.StringIO()
        infile.write(sys.stdin.read())
    else:
        infile = io.open(filename, 'r', encoding='utf-8')

    newfile = io.StringIO()
    reformat_ffile(infile, newfile,
                   orig_filename=filename, **kwargs)

    if diffonly:
        infile.seek(0)
        newfile.seek(0)
        diff_contents=diff(infile.read(),newfile.read(),filename,filename)
        sys.stdout.write(diff_contents)
    else:

        if stdout:
            sys.stdout.write(newfile.getvalue())
        else:
            outfile = io.open(filename, 'r', encoding='utf-8')

            # write to outfile only if content has changed

            import hashlib
            hash_new = hashlib.md5()
            hash_new.update(newfile.getvalue().encode('utf-8'))
            hash_old = hashlib.md5()
            hash_old.update(outfile.read().encode('utf-8'))

            outfile.close()

            if hash_new.digest() != hash_old.digest():
                outfile = io.open(filename, 'w', encoding='utf-8')
                outfile.write(newfile.getvalue())

def reformat_ffile(infile, outfile, impose_indent=True, indent_size=3, strict_indent=False, impose_whitespace=True,
                   case_dict={},
                   impose_replacements=False, cstyle=False, whitespace=2, whitespace_dict={}, llength=132,
                   strip_comments=False, format_decl=False, orig_filename=None, indent_fypp=True, indent_mod=True):
    """main method to be invoked for formatting a Fortran file."""

    # note: whitespace formatting and indentation may require different parsing rules
    # (e.g. preprocessor statements may be indented but not whitespace formatted)
    # therefore we invoke reformat_ffile independently for:
    # 1) whitespace formatting
    # 2) indentation

    if not orig_filename:
        orig_filename = infile.name

    # 1) whitespace formatting
    oldfile = infile
    newfile = infile

    if impose_whitespace:
        _impose_indent = False

        newfile = io.StringIO()
        reformat_ffile_combined(oldfile, newfile, _impose_indent, indent_size, strict_indent, impose_whitespace,
                                case_dict,
                                impose_replacements, cstyle, whitespace, whitespace_dict, llength,
                                strip_comments, format_decl, orig_filename, indent_fypp, indent_mod)
        oldfile = newfile

    # 2) indentation
    if impose_indent:

        _impose_whitespace = False
        _impose_replacements = False

        newfile = io.StringIO()
        reformat_ffile_combined(oldfile, newfile, impose_indent, indent_size, strict_indent, _impose_whitespace,
                                case_dict,
                                _impose_replacements, cstyle, whitespace, whitespace_dict, llength,
                                strip_comments, format_decl, orig_filename, indent_fypp, indent_mod)


    outfile.write(newfile.getvalue())


def reformat_ffile_combined(infile, outfile, impose_indent=True, indent_size=3, strict_indent=False, impose_whitespace=True,
                            case_dict={},
                            impose_replacements=False, cstyle=False, whitespace=2, whitespace_dict={}, llength=132,
                            strip_comments=False, format_decl=False, orig_filename=None, indent_fypp=True, indent_mod=True):

    if not orig_filename:
        orig_filename = infile.name

    if not impose_indent:
        indent_fypp = False

    scope_parser = build_scope_parser(fypp=indent_fypp, mod=indent_mod)

    infile.seek(0)
    req_indents, first_indent = inspect_ffile_format(
        infile, indent_size, strict_indent, indent_fypp, orig_filename)
    infile.seek(0)

    # initialization

    # special cases for indentation:
    # indent_special = 0: parse syntax and impose indent
    # indent_special = 1: no indentation
    # indent_special = 2: use indent from previous line
    # indent_special = 3: take indent from input file (leave as is)
    indent_special = 0

    if impose_indent:
        indenter = F90Indenter(scope_parser, first_indent, indent_size, orig_filename)
    else:
        indent_special = 3

    impose_case = not all(v == 0 for v in case_dict.values())

    nfl = 0  # fortran line counter
    use_same_line = False
    stream = InputStream(infile, not indent_fypp, orig_filename=orig_filename)
    skip_blank = False
    in_format_off_block = False

    while 1:
        f_line, comments, lines = stream.next_fortran_line()

        if not lines:
            break

        nfl += 1
        orig_lines = lines

        if indent_special != 3:
            indent = [0] * len(lines)
        else:
            indent = [len(l) - len((l.lstrip(' ')).lstrip('&'))  for l in lines]

        comment_lines = format_comments(lines, comments, strip_comments)

        auto_align, auto_format, in_format_off_block = parse_fprettify_directives(
            lines, comment_lines, in_format_off_block, orig_filename, stream.line_nr)
        f_line, lines, is_omp_conditional = preprocess_omp(
            f_line, lines)
        f_line, lines, label = preprocess_labels(f_line, lines)

        lines, do_format, prev_indent, is_blank, is_special = preprocess_line(
            f_line, lines, comments, orig_filename, stream.line_nr, indent_fypp)

        if is_special[0]:
            indent_special = 3

        if prev_indent and indent_special == 0:
            indent_special = 2

        if is_blank and skip_blank:
            continue
        if (not do_format):
            if indent_special == 2:
                # inherit indent from previous line
                indent[:] = [indenter.get_fline_indent()]*len(indent)
            elif indent_special == 0:
                indent_special = 1
        else:

            if not auto_align:
                manual_lines_indent = get_manual_alignment(lines)
            else:
                manual_lines_indent = []

            lines, pre_ampersand, ampersand_sep = remove_pre_ampersands(
                lines, is_special, orig_filename, stream.line_nr)

            linebreak_pos = get_linebreak_pos(lines, filter_fypp=not indent_fypp)

            f_line = f_line.strip(' ')

            if impose_replacements:
                f_line = replace_relational_single_fline(f_line, cstyle)

            if impose_case:
                f_line = replace_keywords_single_fline(f_line, case_dict)

            if impose_whitespace:
                lines = format_single_fline(
                    f_line, whitespace, whitespace_dict, linebreak_pos, ampersand_sep,
                    scope_parser, format_decl, orig_filename, stream.line_nr, auto_format)

                lines = append_comments(lines, comment_lines, is_special)

            # target indent for next line
            rel_indent = req_indents[nfl] if nfl < len(req_indents) else 0

            if indent_special != 3:
                indenter.process_lines_of_fline(
                    f_line, lines, rel_indent, indent_size,
                    stream.line_nr, indent_fypp, manual_lines_indent)
                indent = indenter.get_lines_indent()

            lines, indent = prepend_ampersands(lines, indent, pre_ampersand)

        if any(is_special):
            for pos, line in enumerate(lines):
                if is_special[pos]:
                    indent[pos] = len(line) - len(line.lstrip(' '))
                    lines[pos] = line.lstrip(' ')

        lines = remove_trailing_whitespace(lines)

        write_formatted_line(outfile, indent, lines, orig_lines, indent_special, llength,
                             use_same_line, is_omp_conditional, label, orig_filename, stream.line_nr)

        do_indent, use_same_line = pass_defaults_to_next_line(f_line)

        if impose_indent:
            if do_indent:
                indent_special = 0
            else:
                indent_special = 1

        # rm subsequent blank lines
        skip_blank = EMPTY_RE.search(
            f_line) and not any(comments) and not is_omp_conditional and not label


def format_comments(lines, comments, strip_comments):
    comments_ftd = []
    for line, comment in zip(lines, comments):
        has_comment = bool(comment.strip())
        if has_comment:
            if strip_comments:
                sep = not comment.strip() == line.strip()
            else:
                line_minus_comment = line.replace(comment,"")
                sep = len(line_minus_comment.rstrip('\n')) - len(line_minus_comment.rstrip())
        else:
            sep = 0

        if line.strip():  # empty lines between linebreaks are ignored
            comments_ftd.append(' ' * sep + comment.strip())
    return comments_ftd


def parse_fprettify_directives(lines, comment_lines, in_format_off_block, filename, line_nr):
    """
    parse formatter directives '!&' and line continuations starting with an
    ampersand.
    """
    auto_align = not any(NO_ALIGN_RE.search(_) for _ in lines)
    auto_format = not (in_format_off_block or any(
        _.lstrip().startswith('!&') for _ in comment_lines))
    if not auto_format:
        auto_align = False
    if (len(lines)) == 1:
        valid_directive = True
        if lines[0].strip().startswith('!&<'):
            if in_format_off_block:
                valid_directive = False
            else:
                in_format_off_block = True
        if lines[0].strip().startswith('!&>'):
            if not in_format_off_block:
                valid_directive = False
            else:
                in_format_off_block = False
        if not valid_directive:
            raise FprettifyParseException(
                FORMATTER_ERROR_MESSAGE, filename, line_nr)

    return [auto_align, auto_format, in_format_off_block]


def preprocess_omp(f_line, lines):
    """convert omp conditional to normal fortran"""

    is_omp_conditional = bool(OMP_COND_RE.search(f_line))
    if is_omp_conditional:
        f_line = OMP_COND_RE.sub('   ', f_line, count=1)
        lines = [OMP_COND_RE.sub('   ', l, count=1) for l in lines]

    return [f_line, lines, is_omp_conditional]

def preprocess_labels(f_line, lines):
    """remove statement labels"""

    match = STATEMENT_LABEL_RE.search(f_line)
    if match:
        label = match.group(1)
    else:
        label = ''

    if label:
        f_line = STATEMENT_LABEL_RE.sub(len(label)*' ', f_line, count=1)
        lines[0] = STATEMENT_LABEL_RE.sub(len(label)*' ', lines[0], count=1)

    return [f_line, lines, label]

def preprocess_line(f_line, lines, comments, filename, line_nr, indent_fypp):
    """preprocess lines: identification and formatting of special cases"""
    is_blank = False
    prev_indent = False
    do_format = False

    # is_special: special directives that should not be treated as Fortran
    # currently supported: fypp preprocessor directives or comments for FORD documentation
    is_special = [False]*len(lines)

    for pos, line in enumerate(lines):
        line_strip = line.lstrip()
        if indent_fypp:
            is_special[pos] = line_strip.startswith('!!') or \
                              (FYPP_LINE_RE.search(line_strip) if pos > 0 else False)
        else:
            is_special[pos] = FYPP_LINE_RE.search(line_strip) or line_strip.startswith('!!')

    # if first line is special, all lines should be special
    if is_special[0]: is_special = [True]*len(lines)

    if EMPTY_RE.search(f_line):  # empty lines including comment lines
        if any(comments):
            if lines[0].startswith(' '):
                # indent comment lines only if they were not indented before.
                prev_indent = True
        else:
            is_blank = True
        lines = [l.strip(' ') if not is_special[n] else l for n, l in enumerate(lines)]
    else:
        do_format = True

    return [lines, do_format, prev_indent, is_blank, is_special]


def pass_defaults_to_next_line(f_line):
    """defaults to be transferred from f_line to next line"""
    if re.search(r";\s*$", f_line, RE_FLAGS):
        # if line ended with semicolon, don't indent next line
        do_indent = False
        use_same_line = True
    else:
        do_indent = True
        use_same_line = False

    return [do_indent, use_same_line]


def remove_trailing_whitespace(lines):
    """remove trailing whitespaces from lines"""
    lines = [re.sub(r"\s+$", '\n', l, RE_FLAGS)
             for l in lines]
    return lines


def prepend_ampersands(lines, indent, pre_ampersand):
    """prepend ampersands and correct indent"""
    for pos, line in enumerate(lines):
        amp_insert = pre_ampersand[pos]
        if amp_insert:
            indent[pos] += -1
            lines[pos] = amp_insert + line.lstrip()

    return [lines, indent]


def append_comments(lines, comment_lines, is_special):
    """append comments to lines"""
    for pos, (line, comment) in enumerate(zip(lines, comment_lines)):
        if pos < len(lines) - 1:
            has_nl = True  # has next line
            if not line.strip() and not is_special[pos]: comment = comment.lstrip()
        else:
            has_nl = not re.search(EOL_SC, line)
        lines[pos] = lines[pos].rstrip(' ') + comment + '\n' * has_nl

    return lines


def get_linebreak_pos(lines, filter_fypp=True):
    """extract linebreak positions in Fortran line from lines"""
    linebreak_pos = []
    if filter_fypp:
        notfortran_re = NOTFORTRAN_LINE_RE
    else:
        notfortran_re = NOTFORTRAN_FYPP_LINE_RE

    for line in lines:
        found = None
        for char_pos, _ in CharFilter(line, filter_strings=False):
            if re.match(LINEBREAK_STR, line[char_pos:], RE_FLAGS):
                found = char_pos
        if found:
            linebreak_pos.append(found)
        elif notfortran_re.search(line.lstrip(' ')):
            linebreak_pos.append(0)

    linebreak_pos = [sum(linebreak_pos[0:_ + 1]) -
                     1 for _ in range(0, len(linebreak_pos))]

    return linebreak_pos


def remove_pre_ampersands(lines, is_special, filename, line_nr):
    """
    remove and return preceding ampersands ('pre_ampersand'). Also return
    number of whitespace characters before ampersand of previous line
    ('ampersand_sep').

    Note: Don't do any whitespace formatting on ampersands if next line starts
    with an ampersand but remember the original number of spaces
    (ampersand_sep). This "special rule" is necessary since ampersands starting
    a line can be used to break literals, so changing the number of whitespaces
    before the ampersand ending the previous line may lead to invalid syntax or
    may change the number of whitespace characters in a string.
    """
    pre_ampersand = []
    ampersand_sep = []

    for pos, line in enumerate(lines):
        match = re.search(SOL_STR + r'(&\s*)', line)
        if match:
            pre_ampersand.append(match.group(1))
            # amount of whitespace before ampersand of previous line:
            m = re.search(r'(\s*)&[\s]*(?:!.*)?$', lines[pos - 1])
            if not m:
                raise FprettifyParseException(
                    "Bad continuation line format", filename, line_nr)
            sep = len(m.group(1))

            ampersand_sep.append(sep)
        else:
            pre_ampersand.append('')
            if pos > 0:
                # use default 1 whitespace character before ampersand
                ampersand_sep.append(1)

    lines = [l.strip(' ').strip('&') if not s else l for l, s in zip(lines, is_special)]
    return [lines, pre_ampersand, ampersand_sep]


def get_manual_alignment(lines):
    """extract manual indents for line continuations from line"""
    manual_lines_indent = [
        len(l) - len(l.lstrip(' ').lstrip('&')) for l in lines]
    manual_lines_indent = [ind - manual_lines_indent[0]
                           for ind in manual_lines_indent]
    return manual_lines_indent


def write_formatted_line(outfile, indent, lines, orig_lines, indent_special, llength, use_same_line, is_omp_conditional, label, filename, line_nr):
    """Write reformatted line to file"""

    for ind, line, orig_line in zip(indent, lines, orig_lines):

        # get actual line length excluding comment:
        line_length = 0
        for line_length, _ in CharFilter(line):
            pass
        line_length += 1

        if indent_special != 1:
            ind_use = ind
        else:
            if use_same_line:
                ind_use = 1
            else:
                ind_use = 0

        if CPP_RE.search(line.lstrip()):
            ind_use = 0

        if label:
            label_use = label + ' '
            label = '' # no label for continuation lines
        else:
            label_use = ''

        if ind_use + line_length <= (llength+1):  # llength (default 132) plus 1 newline char
            outfile.write('!$ ' * is_omp_conditional + label_use +
                          ' ' * (ind_use - 3 * is_omp_conditional - len(label_use) +
                                 len(line) - len(line.lstrip(' '))) +
                          line.lstrip(' '))
        elif line_length <= (llength+1):
            outfile.write('!$ ' * is_omp_conditional + label_use + ' ' *
                          ((llength+1) - 3 * is_omp_conditional - len(label_use) -
                           len(line.lstrip(' '))) + line.lstrip(' '))

            log_message(LINESPLIT_MESSAGE+" (limit: "+str(llength)+")", "warning",
                        filename, line_nr)
        else:
            outfile.write(orig_line)
            log_message(LINESPLIT_MESSAGE+" (limit: "+str(llength)+")", "warning",
                        filename, line_nr)


def get_curr_delim(line, pos):
    """get delimiter token in line starting at pos, if it exists"""
    what_del_open = DEL_OPEN_RE.search(line[pos:pos + 2])
    what_del_close = DEL_CLOSE_RE.search(line[pos:pos + 2])
    return [what_del_open, what_del_close]


def set_fprettify_logger(level):
    """setup custom logger"""
    logger = logging.getLogger('fprettify-logger')
    logger.setLevel(level)
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(level)
    formatter = logging.Formatter(
        '%(levelname)s: File %(ffilename)s, line %(fline)s\n    %(message)s')
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)


def log_exception(e, message):
    """log an exception and a message"""
    log_message(message, "exception", e.filename, e.line_nr)


def log_message(message, level, filename, line_nr):
    """log a message"""

    logger = logging.getLogger('fprettify-logger')
    logger_d = {'ffilename': filename, 'fline': line_nr}
    logger_to_use = getattr(logger, level)
    logger_to_use(message, extra=logger_d)


def run(argv=sys.argv):  # pragma: no cover
    """Command line interface"""

    try:
        import configargparse as argparse
    except ImportError:
        import argparse

    def str2bool(str):
        """helper function to convert strings to bool"""
        if str.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif str.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            return None

    def get_config_file_list(filename):
        """helper function to create list of config files found in parent directories"""
        config_file_list = []
        dir = os.path.dirname(filename)
        while True:
            config_file = os.path.join(dir, '.fprettify.rc')
            if os.path.isfile(config_file):
                config_file_list.insert(0, config_file)
            parent=os.path.dirname(dir)
            if parent == dir:
                break
            dir = parent
        return config_file_list

    arguments = {'prog': argv[0],
                 'description': 'Auto-format modern Fortran source files.',
                 'formatter_class': argparse.ArgumentDefaultsHelpFormatter}

    if argparse.__name__ == "configargparse":
        arguments['args_for_setting_config_path'] = ['-c', '--config-file']
        arguments['description'] = arguments['description'] + " Config files ('.fprettify.rc') in the home (~) directory and any such files located in parent directories of the input file will be used. When the standard input is used, the search is started from the current directory."

    def get_arg_parser(args):
        """helper function to create the parser object"""
        parser = argparse.ArgumentParser(**args)

        parser.add_argument("-i", "--indent", type=int, default=3,
                            help="relative indentation width")
        parser.add_argument("-l", "--line-length", type=int, default=132,
                            help="column after which a line should end, viz. -ffree-line-length-n for GCC")
        parser.add_argument("-w", "--whitespace", type=int,
                            choices=range(0, 5), default=2, help="Presets for the amount of whitespace - "
                                                                 "   0: minimal whitespace"
                                                                 " | 1: operators (except arithmetic), print/read"
                                                                 " | 2: operators, print/read, plus/minus"
                                                                 " | 3: operators, print/read, plus/minus, muliply/divide"
                                                                 " | 4: operators, print/read, plus/minus, muliply/divide, type component selector")
        parser.add_argument("--whitespace-comma", type=str2bool, nargs="?", default="None", const=True,
                            help="boolean, en-/disable whitespace for comma/semicolons")
        parser.add_argument("--whitespace-assignment", type=str2bool, nargs="?", default="None", const=True,
                            help="boolean, en-/disable whitespace for assignments")
        parser.add_argument("--whitespace-decl", type=str2bool, nargs="?", default="None", const=True,
                            help="boolean, en-/disable whitespace for declarations (requires '--enable-decl')")
        parser.add_argument("--whitespace-relational", type=str2bool, nargs="?", default="None", const=True,
                            help="boolean, en-/disable whitespace for relational operators")
        parser.add_argument("--whitespace-logical", type=str2bool, nargs="?", default="None", const=True,
                            help="boolean, en-/disable whitespace for logical operators")
        parser.add_argument("--whitespace-plusminus", type=str2bool, nargs="?", default="None", const=True,
                            help="boolean, en-/disable whitespace for plus/minus arithmetic")
        parser.add_argument("--whitespace-multdiv", type=str2bool, nargs="?", default="None", const=True,
                            help="boolean, en-/disable whitespace for multiply/divide arithmetic")
        parser.add_argument("--whitespace-print", type=str2bool, nargs="?", default="None", const=True,
                            help="boolean, en-/disable whitespace for print/read statements")
        parser.add_argument("--whitespace-type", type=str2bool, nargs="?", default="None", const=True,
                            help="boolean, en-/disable whitespace for select type components")
        parser.add_argument("--whitespace-intrinsics", type=str2bool, nargs="?", default="None", const=True,
                            help="boolean, en-/disable whitespace for intrinsics like if/write/close")
        parser.add_argument("--strict-indent", action='store_true', default=False, help="strictly impose indentation even for nested loops")
        parser.add_argument("--enable-decl", action="store_true", default=False, help="enable whitespace formatting of declarations ('::' operator).")
        parser.add_argument("--disable-indent", action='store_true', default=False, help="don't impose indentation")
        parser.add_argument("--disable-whitespace", action='store_true', default=False, help="don't impose whitespace formatting")
        parser.add_argument("--enable-replacements", action='store_true', default=False, help="replace relational operators (e.g. '.lt.' <--> '<')")
        parser.add_argument("--c-relations", action='store_true', default=False, help="C-style relational operators ('<', '<=', ...)")
        parser.add_argument("--case", nargs=4, default=[0,0,0,0], type=int, help="Enable letter case formatting of intrinsics by specifying which of "
                            "keywords, procedures/modules, operators and constants (in this order) should be lowercased or uppercased - "
                            "   0: do nothing"
                            " | 1: lowercase"
                            " | 2: uppercase")

        parser.add_argument("--strip-comments", action='store_true', default=False, help="strip whitespaces before comments")
        parser.add_argument('--disable-fypp', action='store_true', default=False,
                            help="Disables the indentation of fypp preprocessor blocks.")
        parser.add_argument('--disable-indent-mod', action='store_true', default=False,
                            help="Disables the indentation after module / program.")

        parser.add_argument("-d","--diff", action='store_true', default=False,
                             help="Write file differences to stdout instead of formatting inplace")
        parser.add_argument("-s", "--stdout", action='store_true', default=False,
                            help="Write to stdout instead of formatting inplace")

        group = parser.add_mutually_exclusive_group()
        group.add_argument("-S", "--silent", "--no-report-errors", action='store_true',
                           default=False, help="Don't write any errors or warnings to stderr")
        group.add_argument("-D", "--debug", action='store_true',
                           default=False, help=argparse.SUPPRESS)
        parser.add_argument("path", type=str, nargs='*',
                            help="Paths to files to be formatted inplace. If no paths are given, stdin (-) is used by default. Path can be a directory if --recursive is used.", default=['-'])
        parser.add_argument('-r', '--recursive', action='store_true',
                            default=False, help="Recursively auto-format all Fortran files in subdirectories of specified path; recognized filename extensions: {}". format(", ".join(FORTRAN_EXTENSIONS)))
        parser.add_argument('-e', '--exclude', action='append',
                            default=[], type=str,
                            help="File or directory patterns to be excluded when searching for Fortran files to format")
        parser.add_argument('-f', '--fortran', type=str, action='append', default=[],
                            help="Overrides default fortran extensions recognized by --recursive. Repeat this option to specify more than one extension.")
        parser.add_argument('--version', action='version',
                            version='%(prog)s 0.3.7')
        return parser

    parser = get_arg_parser(arguments)

    args = parser.parse_args(argv[1:])

    def build_ws_dict(args):
        """helper function to build whitespace dictionary"""
        ws_dict = {}
        ws_dict['comma'] = args.whitespace_comma
        ws_dict['assignments'] = args.whitespace_assignment
        ws_dict['decl'] = args.whitespace_decl
        ws_dict['relational'] = args.whitespace_relational
        ws_dict['logical'] = args.whitespace_logical
        ws_dict['plusminus'] = args.whitespace_plusminus
        ws_dict['multdiv'] = args.whitespace_multdiv
        ws_dict['print'] = args.whitespace_print
        ws_dict['type'] = args.whitespace_type
        ws_dict['intrinsics'] = args.whitespace_intrinsics
        return ws_dict

    case_dict = {
            'keywords' : args.case[0],
            'procedures' : args.case[1],
            'operators' : args.case[2],
            'constants' : args.case[3]
            }

    # support legacy input:
    if 'stdin' in args.path and not os.path.isfile('stdin'):
        args.path = ['-' if _ == 'stdin' else _ for _ in args.path]

    for directory in args.path:
        if directory == '-':
            if args.recursive:
                sys.stderr.write("--recursive requires a directory.\n")
                sys.exit()
        else:
            if not os.path.exists(directory):
                sys.stderr.write("directory " + directory +
                                 " does not exist!\n")
                sys.exit()
            if not os.path.isfile(directory) and directory != '-' and not args.recursive:
                sys.stderr.write("file " + directory + " does not exist!\n")
                sys.exit()

        if not args.recursive:
            filenames = [directory]
        else:
            if args.fortran:
                ext = args.fortran
            else:
                ext = FORTRAN_EXTENSIONS
            filenames = []

            from fnmatch import fnmatch

            for dirpath, dirnames, files in os.walk(directory,topdown=True):

                # Prune excluded patterns from list of child directories
                dirnames[:] = [dirname for dirname in dirnames if not any(
                    [fnmatch(dirname,exclude_pattern) or fnmatch(os.path.join(dirpath,dirname),exclude_pattern)
                            for exclude_pattern in args.exclude]
                )]

                for ffile in [os.path.join(dirpath, f) for f in files
                              if any(f.endswith(_) for _ in ext)
                              and not any([
                                  fnmatch(f,exclude_pattern)
                                  for exclude_pattern in args.exclude])]:
                    filenames.append(ffile)

        for filename in filenames:

            # reparse arguments using the file's list of config files
            filearguments = arguments
            if argparse.__name__ == "configargparse":
                filearguments['default_config_files'] = ['~/.fprettify.rc'] + get_config_file_list(os.path.abspath(filename) if filename != '-' else os.getcwd())
            file_argparser = get_arg_parser(filearguments)
            file_args = file_argparser.parse_args(argv[1:])
            ws_dict = build_ws_dict(file_args)

            stdout = file_args.stdout or directory == '-'
            diffonly=file_args.diff

            if file_args.debug:
                level = logging.DEBUG
            elif args.silent:
                level = logging.CRITICAL
            else:
                level = logging.WARNING

            set_fprettify_logger(level)

            try:
                reformat_inplace(filename,
                                 stdout=stdout,
                                 diffonly=diffonly,
                                 impose_indent=not file_args.disable_indent,
                                 indent_size=file_args.indent,
                                 strict_indent=file_args.strict_indent,
                                 impose_whitespace=not file_args.disable_whitespace,
                                 impose_replacements=file_args.enable_replacements,
                                 cstyle=file_args.c_relations,
                                 case_dict=case_dict,
                                 whitespace=file_args.whitespace,
                                 whitespace_dict=ws_dict,
                                 llength=1024 if file_args.line_length == 0 else file_args.line_length,
                                 strip_comments=file_args.strip_comments,
                                 format_decl=file_args.enable_decl,
                                 indent_fypp=not file_args.disable_fypp,
                                 indent_mod=not file_args.disable_indent_mod)
            except FprettifyException as e:
                log_exception(e, "Fatal error occured")
