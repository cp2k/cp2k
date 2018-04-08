#!/usr/bin/env python
# -*- coding: utf-8 -*-
###############################################################################
#    This file is part of fprettify.
#    Copyright (C) 2016-2018 Patrick Seewald, CP2K developers group
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
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import re
import sys
import logging
import os
import io
import argparse

# allow for unicode for stdin / stdout, it's a mess
try:
    # python 3
    sys.stdin = io.TextIOWrapper(
        sys.stdin.detach(), encoding='UTF-8', line_buffering=True)
    sys.stdout = io.TextIOWrapper(
        sys.stdout.detach(), encoding='UTF-8', line_buffering=True)
except AttributeError:
    # python 2
    import codecs
    utf8_reader = codecs.getreader('UTF-8')
    utf8_writer = codecs.getwriter('UTF-8')
    sys.stdin = utf8_reader(sys.stdin)
    sys.stdout = utf8_writer(sys.stdout)

from .fparse_utils import (VAR_DECL_RE, OMP_RE, OMP_DIR_RE,
                           InputStream, CharFilter,
                           FprettifyException, FprettifyParseException, FprettifyInternalException, RE_FLAGS)

# recognize fortran files by extension
FORTRAN_EXTENSIONS = [".f", ".for", ".ftn",
                      ".f90", ".f95", ".f03", ".fpp"]
FORTRAN_EXTENSIONS += [_.upper() for _ in FORTRAN_EXTENSIONS]

# constants, mostly regular expressions:
FORMATTER_ERROR_MESSAGE = (" Wrong usage of formatting-specific directives"
                           " '&', '!&', '!&<' or '!&>'.")
LINESPLIT_MESSAGE = ("auto indentation failed due to 132 chars limit, "
                     "line should be splitted")

EOL_STR = r"\s*;?\s*$"  # end of fortran line
EOL_SC = r"\s*;\s*$"  # whether line is ended with semicolon
SOL_STR = r"^\s*"  # start of fortran line

F77_STYLE = re.compile(r"^\s*\d", RE_FLAGS)

# regular expressions for parsing statements that start, continue or end a
# subunit:
IF_RE = re.compile(
    SOL_STR + r"(\w+\s*:)?\s*IF\s*\(.+\)\s*THEN" + EOL_STR, RE_FLAGS)
ELSE_RE = re.compile(
    SOL_STR + r"ELSE(\s*IF\s*\(.+\)\s*THEN)?" + EOL_STR, RE_FLAGS)
ENDIF_RE = re.compile(SOL_STR + r"END\s*IF(\s+\w+)?" + EOL_STR, RE_FLAGS)

DO_RE = re.compile(SOL_STR + r"(\w+\s*:)?\s*DO(" + EOL_STR + r"|\s)", RE_FLAGS)
ENDDO_RE = re.compile(SOL_STR + r"END\s*DO(\s+\w+)?" + EOL_STR, RE_FLAGS)

SELCASE_RE = re.compile(
    SOL_STR + r"SELECT\s*(CASE|TYPE)\s*\(.+\)" + EOL_STR, RE_FLAGS)
CASE_RE = re.compile(
    SOL_STR + r"((CASE|TYPE\s+IS|CLASS\s+IS)\s*(\(.+\)|DEFAULT)|CLASS\s+DEFAULT)" + EOL_STR, RE_FLAGS)
ENDSEL_RE = re.compile(SOL_STR + r"END\s*SELECT" + EOL_STR, RE_FLAGS)

ASSOCIATE_RE = re.compile(SOL_STR + r"ASSOCIATE\s*\(.+\)" + EOL_STR, RE_FLAGS)
ENDASSOCIATE_RE = re.compile(SOL_STR + r"END\s*ASSOCIATE" + EOL_STR, RE_FLAGS)

BLK_RE = re.compile(SOL_STR + r"(\w+\s*:)?\s*BLOCK" + EOL_STR, RE_FLAGS)
ENDBLK_RE = re.compile(SOL_STR + r"END\s*BLOCK(\s+\w+)?" + EOL_STR, RE_FLAGS)

SUBR_RE = re.compile(
    r"^([^\"'!]* )?SUBROUTINE\s+\w+\s*(\(.*\))?" + EOL_STR, RE_FLAGS)
ENDSUBR_RE = re.compile(
    SOL_STR + r"END\s*SUBROUTINE(\s+\w+)?" + EOL_STR, RE_FLAGS)

FCT_RE = re.compile(
    r"^([^\"'!]* )?FUNCTION\s+\w+\s*(\(.*\))?(\s*RESULT\s*\(\w+\))?" + EOL_STR,
    RE_FLAGS)
ENDFCT_RE = re.compile(
    SOL_STR + r"END\s*FUNCTION(\s+\w+)?" + EOL_STR, RE_FLAGS)

MOD_RE = re.compile(SOL_STR + r"MODULE\s+\w+" + EOL_STR, RE_FLAGS)
ENDMOD_RE = re.compile(SOL_STR + r"END\s*MODULE(\s+\w+)?" + EOL_STR, RE_FLAGS)

TYPE_RE = re.compile(
    SOL_STR +
    r"TYPE(\s*,\s*(BIND\s*\(\s*C\s*\)|EXTENDS\s*\(.*\)|ABSTRACT))?(\s*,\s*(PUBLIC|PRIVATE))?(\s*::\s*|\s+)\w+" + EOL_STR,
    RE_FLAGS)
ENDTYPE_RE = re.compile(SOL_STR + r"END\s*TYPE(\s+\w+)?" + EOL_STR, RE_FLAGS)

PROG_RE = re.compile(SOL_STR + r"PROGRAM\s+\w+" + EOL_STR, RE_FLAGS)
ENDPROG_RE = re.compile(
    SOL_STR + r"END\s*PROGRAM(\s+\w+)?" + EOL_STR, RE_FLAGS)

INTERFACE_RE = re.compile(
    r"^([^\"'!]* )?INTERFACE(\s+\w+|\s+(OPERATOR|ASSIGNMENT)\s*\(.*\))?" + EOL_STR, RE_FLAGS)
ENDINTERFACE_RE = re.compile(
    SOL_STR + r"END\s*INTERFACE(\s+\w+)?" + EOL_STR, RE_FLAGS)

CONTAINS_RE = re.compile(SOL_STR + r"CONTAINS" + EOL_STR, RE_FLAGS)

ENUM_RE = re.compile(
    SOL_STR + r"ENUM(\s*,\s*(BIND\s*\(\s*C\s*\)))?(\s*::\s*|\s+)\w+" + EOL_STR,
    RE_FLAGS)
ENDENUM_RE = re.compile(SOL_STR + r"END\s*ENUM(\s+\w+)?" + EOL_STR, RE_FLAGS)

ENDANY_RE = re.compile(SOL_STR + r"END" + EOL_STR, RE_FLAGS)

PRIVATE_RE = re.compile(SOL_STR + r"PRIVATE\s*::", RE_FLAGS)
PUBLIC_RE = re.compile(SOL_STR + r"PUBLIC\s*::", RE_FLAGS)

# intrinsic statements with parenthesis notation that are not functions
INTR_STMTS_PAR = (r"(ALLOCATE|DEALLOCATE|REWIND|BACKSPACE|INQUIRE|"
                  r"OPEN|CLOSE|READ|WRITE|"
                  r"FORALL|WHERE|ASSOCIATE|NULLIFY)")

# regular expressions for parsing linebreaks
LINEBREAK_STR = r"(&)[\s]*(?:!.*)?$"

# regular expressions for parsing operators
# Note: +/- in real literals and sign operator is ignored
PLUSMINUS_RE = re.compile(
    r"(?<=[\w\)\]])(?<![\d\.][de])\s*(\+|-)\s*", RE_FLAGS)
# Note: ** or // (or any multiples of * or /) are ignored
MULTDIV_RE = re.compile(
    r"(?<=[\w\)\]])\s*((?<!\*)\*(?!\*)|(?<!/)/(?!/))(?=[\s\w\(])", RE_FLAGS)
REL_OP_RE = re.compile(
    r"(?<!\()\s*(\.(?:EQ|NE|LT|LE|GT|GE)\.|(?:==|\/=|<(?!=)|<=|(?<!=)>(?!=)|>=))\s*(?!\))",
    RE_FLAGS)
LOG_OP_RE = re.compile(r"\s*(\.(?:AND|OR|EQV|NEQV)\.)\s*", RE_FLAGS)
PRINT_RE = re.compile(r"(?<=\w)\s*(\*,)\s*", RE_FLAGS)

# regular expressions for parsing delimiters
DEL_OPEN_STR = r"(\(\/?|\[)"
DEL_OPEN_RE = re.compile(r"^" + DEL_OPEN_STR, RE_FLAGS)
DEL_CLOSE_STR = r"(\/?\)|\])"
DEL_CLOSE_RE = re.compile(r"^" + DEL_CLOSE_STR, RE_FLAGS)

# empty line regex
EMPTY_RE = re.compile(SOL_STR + r"([!#].*)?$", RE_FLAGS)

# two-sided operators
LR_OPS_RE = [REL_OP_RE, LOG_OP_RE, PLUSMINUS_RE, MULTDIV_RE, PRINT_RE]

USE_RE = re.compile(
    SOL_STR + "USE(\s+|(,.+?)?::\s*)\w+?((,.+?=>.+?)+|,\s*only\s*:.+?)?$" + EOL_STR, RE_FLAGS)

# markups to deactivate formatter
NO_ALIGN_RE = re.compile(SOL_STR + r"&\s*[^\s*]+")

# combine regex that define subunits
NEW_SCOPE_RE = [IF_RE, DO_RE, SELCASE_RE, SUBR_RE,
                FCT_RE, MOD_RE, PROG_RE, INTERFACE_RE, TYPE_RE, ENUM_RE, ASSOCIATE_RE, None, BLK_RE]
CONTINUE_SCOPE_RE = [ELSE_RE, None, CASE_RE, CONTAINS_RE,
                     CONTAINS_RE, CONTAINS_RE, CONTAINS_RE, None, CONTAINS_RE, None, None, None, None]
END_SCOPE_RE = [ENDIF_RE, ENDDO_RE, ENDSEL_RE, ENDSUBR_RE,
                ENDFCT_RE, ENDMOD_RE, ENDPROG_RE, ENDINTERFACE_RE, ENDTYPE_RE, ENDENUM_RE, ENDASSOCIATE_RE, ENDANY_RE, ENDBLK_RE]

# match namelist names
NML_RE = re.compile(r"(/\w+/)", RE_FLAGS)
# find namelists and data statements
NML_STMT_RE = re.compile(SOL_STR + r"NAMELIST.*/.*/", RE_FLAGS)
DATA_STMT_RE = re.compile(SOL_STR + r"DATA\s+\w", RE_FLAGS)

class F90Indenter(object):
    """
    Parses encapsulation of subunits / scopes line by line
    and updates the indentation.
    """

    def __init__(self, first_indent, rel_indent, filename):
        # scopes / subunits:
        self._scope_storage = []
        # indents for all fortran lines:
        self._indent_storage = []
        # indents of actual lines of current fortran line
        self._line_indents = []

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
                               line_nr, manual_lines_indent=None):
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
        :manual_lines_indent: don't use F90Aligner but manually impose
                              indents for continuations
        """

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

        for new_n, newre in enumerate(NEW_SCOPE_RE):
            if newre and newre.search(f_line) and not END_SCOPE_RE[new_n].search(f_line):
                what_new = new_n
                is_new = True
                valid_new = True
                scopes.append(what_new)
                log_message("{}: {}".format(what_new, f_line),
                            "debug", filename, line_nr)

        # check statements that continue scope
        is_con = False
        valid_con = False
        for con_n, conre in enumerate(CONTINUE_SCOPE_RE):
            if conre and conre.search(f_line):
                what_con = con_n
                is_con = True
                if len(scopes) > 0:
                    what = scopes[-1]
                    if what == what_con:
                        valid_con = True
                        log_message("{}: {}".format(
                            what_con, f_line), "debug", filename, line_nr)

        # check statements that end scope
        is_end = False
        valid_end = False
        for end_n, endre in enumerate(END_SCOPE_RE):
            if endre and endre.search(f_line):
                what_end = end_n
                is_end = True
                if len(scopes) > 0:
                    what = scopes.pop()
                    if what == what_end:
                        valid_end = True
                        log_message("{}: {}".format(
                            what_end, f_line), "debug", filename, line_nr)

        # deal with line breaks
        if not manual_lines_indent:
            self._aligner.process_lines_of_fline(
                f_line, lines, rel_ind_con, line_nr)
            br_indent_list = self._aligner.get_lines_indent()
        else:
            br_indent_list = manual_lines_indent

        for pos in range(0, len(lines) - 1):
            line_indents[pos + 1] = br_indent_list[pos + 1]

        if is_new:
            if not valid_new:
                log_message('invalid scope opening statement',
                            "info", filename, line_nr)

            line_indents = [ind + indents[-1] for ind in line_indents]

            indents.append(rel_ind + indents[-1])

        elif is_con or is_end:
            valid = valid_con if is_con else valid_end

            if not valid:
                log_message('invalid scope closing statement',
                            "info", filename, line_nr)
            try:
                line_indents = [ind + indents[-2 + self._initial]
                                for ind in line_indents]
            except IndexError:
                assert not valid

            if is_end:
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
        for pos, line in enumerate(lines):
            self.__align_line_continuations(
                line, is_decl, rel_ind, self._line_nr + pos)
            if pos + 1 < len(lines):
                self._line_indents.append(self._br_indent_list[-1])

        if len(self._br_indent_list) > 2 or self._level:
            log_message('unpaired bracket delimiters',
                        "info", self._filename, self._line_nr)

    def get_lines_indent(self):
        """after processing, retrieve the indents of all line parts."""
        return self._line_indents

    def __align_line_continuations(self, line, is_decl, indent_size, line_nr):
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

        for pos, char in CharFilter(enumerate(line)):

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


def inspect_ffile_format(infile, indent_size, orig_filename=None):
    """
    Determine indentation by inspecting original Fortran file.

    This is mainly for finding aligned blocks of DO/IF statements.
    Also check if it has f77 constructs.
    :param infile: open file
    :param indent_size: the default indent size, if <= 0,
                        adopt original indents
    :orig_filename: filename used for messages
    :returns: [ target indent sizes for each line,
                indent of first line (offset),
                whether file is sufficiently modern fortran ]
    """
    if not orig_filename:
        orig_filename = infile.name

    adopt = indent_size <= 0

    num_labels = False
    indents = []
    stream = InputStream(infile, orig_filename)
    prev_offset = 0
    first_indent = -1
    f77_line_nr = 0

    while 1:
        f_line, _, lines = stream.next_fortran_line()
        if not lines:
            break

        offset = len(lines[0]) - len(lines[0].lstrip(' '))
        if f_line.strip() and first_indent == -1:
            first_indent = offset
        indents.append(offset - prev_offset)

        # do not adopt indentations but impose fixed rel. ind.
        # but don't impose indentation for blocked do/if constructs:
        if not adopt and (prev_offset != offset or (not IF_RE.search(f_line) and
                                                    not DO_RE.search(f_line))):
            indents[-1] = indent_size
        prev_offset = offset

        if not num_labels and F77_STYLE.search(f_line):
            num_labels = True
            f77_line_nr = stream.line_nr

    modern_fortran = not num_labels

    return indents, first_indent, modern_fortran, f77_line_nr


def format_single_fline(f_line, whitespace, linebreak_pos, ampersand_sep,
                        filename, line_nr, auto_format=True):
    """
    format a single Fortran line - imposes white space formatting
    and inserts linebreaks.
    Takes a logical Fortran line `f_line` as input as well as the positions
    of the linebreaks (`linebreak_pos`), and the number of
    separating whitespace characters before ampersand (`ampersand_sep`).
    `filename` and `line_nr` just for error messages.
    The higher `whitespace`, the more white space characters inserted -
    whitespace = 0, 1, 2, 3 are currently supported.
    auto formatting can be turned off by setting `auto_format` to False.
    """

    # define whether to put whitespaces around operators:
    # 0: comma, semicolon
    # 1: assignment operators
    # 2: relational operators
    # 3: logical operators
    # 4: arithm. operators plus and minus
    # 5: arithm. operators multiply and divide
    # 6: print / read statements

    if whitespace == 0:
        spacey = [0, 0, 0, 0, 0, 0, 0]
    elif whitespace == 1:
        spacey = [1, 1, 1, 1, 0, 0, 1]
    elif whitespace == 2:
        spacey = [1, 1, 1, 1, 1, 0, 1]
    elif whitespace == 3:
        spacey = [1, 1, 1, 1, 1, 1, 1]
    else:
        raise NotImplementedError("unknown value for whitespace")

    line = f_line
    line_orig = line

    if auto_format:
        line = rm_extra_whitespace(line)
        line = add_whitespace_charwise(line, spacey, filename, line_nr)
        line = add_whitespace_context(line, spacey)

    lines_out = split_reformatted_line(
        line_orig, linebreak_pos, ampersand_sep, line, filename, line_nr)
    return lines_out


def rm_extra_whitespace(line):
    """rm all unneeded whitespace chars, except for declarations"""
    line_ftd = ''
    pos_prev = -1
    for pos, char in CharFilter(enumerate(line)):
        is_decl = line[pos:].lstrip().startswith('::') or line[
            :pos].rstrip().endswith('::')

        if char == ' ':
            # remove double spaces:
            if line_ftd and (re.search(r'[\w"]', line_ftd[-1]) or is_decl):
                line_ftd = line_ftd + char
        else:
            if (line_ftd and line_ftd[-1] == ' ' and
                    (not re.search(r"[\w\"']", char) and not is_decl)):
                line_ftd = line_ftd[:-1]  # remove spaces except between words
            line_ftd = line_ftd + line[pos_prev + 1:pos + 1]
        pos_prev = pos
    return line_ftd


def add_whitespace_charwise(line, spacey, filename, line_nr):
    """add whitespace character wise (no need for context aware parsing)"""
    line_ftd = line
    pos_eq = []
    end_of_delim = -1
    level = 0
    for pos, char in CharFilter(enumerate(line)):
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
                        re.search(SOL_STR + r"(SELECT)?\s*CASE\s*",
                                  line[:pos], RE_FLAGS) or
                        re.search(SOL_STR + r"SELECT\s*TYPE\s*",
                                  line[:pos], RE_FLAGS) or
                        re.search(SOL_STR + r"CLASS\s*DEFAULT\s*",
                                  line[:pos], RE_FLAGS) or
                        re.search(SOL_STR + r"(TYPE|CLASS)\s+IS\s*",
                                  line[:pos], RE_FLAGS) or
                        re.search(r"\b" + INTR_STMTS_PAR + r"\s*$",
                                  line[:pos], RE_FLAGS)):
                    sep1 = 1

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
                    sep2 = len(rhs) - len(rhs.lstrip(' '))

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
            if not level:  # remember position of assignment operator
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

    if level != 0:
        log_message('unpaired bracket delimiters', "info", filename, line_nr)

    return line_ftd


def add_whitespace_context(line, spacey):
    """
    for context aware whitespace formatting we extract line parts that are
    not comments or strings in order to be able to apply a context aware regex.
    """

    line_parts = []
    str_end = -1
    instring = ''
    for pos, char in enumerate(line):
        if char in ['"', "'"]:  # skip string
            if not instring:
                str_start = pos
                line_parts.append(line[str_end + 1:str_start])
                instring = char
            elif instring == char:
                str_end = pos
                line_parts.append(line[str_start:str_end + 1])
                instring = ''
        if pos == len(line) - 1:
            line_parts.append(line[str_end + 1:])

    # format namelists with spaces around /
    if NML_STMT_RE.match(line):
        for pos, part in enumerate(line_parts):
            # exclude comments, strings:
            if not re.match(r"['\"!]", part, RE_FLAGS):
                partsplit = NML_RE.split(part)
                line_parts[pos] = (' '.join(partsplit))

    # Two-sided operators
    for n_op, lr_re in enumerate(LR_OPS_RE):
        for pos, part in enumerate(line_parts):
            # exclude comments, strings:
            if not re.search(r"^['\"!]", part, RE_FLAGS):
                # also exclude / if we see a namelist and data statement
                if not ( NML_STMT_RE.match(line) or DATA_STMT_RE.match(line) ):
                    partsplit = lr_re.split(part)
                    line_parts[pos] = (' ' * spacey[n_op + 2]).join(partsplit)

    line = ''.join(line_parts)

    # format ':' for labels and use only statements
    for newre in NEW_SCOPE_RE[0:2]:
        if newre.search(line) and re.search(SOL_STR + r"\w+\s*:", line):
            line = ': '.join(_.strip() for _ in line.split(':', 1))

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
    lines_splitted = [(line[l:r].rstrip(' ') +
                       ' ' * ampersand_sep[pos] + '&' * min(1, r - l))
                      for pos, (l, r) in enumerate(zip(linebreak_pos_ftd[0:-1],
                                                       linebreak_pos_ftd[1:]))]

    lines_splitted.append(line[linebreak_pos_ftd[-1]:])

    return lines_splitted


def reformat_inplace(filename, stdout=False, **kwargs):  # pragma: no cover
    """reformat a file in place."""
    if filename == '-':
        infile = io.StringIO()
        infile.write(sys.stdin.read())
    else:
        infile = io.open(filename, 'r', encoding='utf-8')

    newfile = io.StringIO()
    reformat_ffile(infile=infile, outfile=newfile,
                   orig_filename=filename, **kwargs)

    if stdout:
        sys.stdout.write(newfile.getvalue())
    else:
        outfile = io.open(filename, 'w', encoding='utf-8')
        outfile.write(newfile.getvalue())


def reformat_ffile(infile, outfile, indent_size=3, whitespace=2,
                   orig_filename=None):
    """main method to be invoked for formatting a Fortran file."""

    if not orig_filename:
        orig_filename = infile.name

    infile.seek(0)
    req_indents, first_indent, modern, f77_line_nr = inspect_ffile_format(
        infile, indent_size, orig_filename)
    infile.seek(0)

    if not modern:
        raise FprettifyParseException(
            "fprettify failed because of fixed format or f77 constructs.", orig_filename, f77_line_nr)

    # initialization
    indenter = F90Indenter(first_indent, indent_size, orig_filename)
    nfl = 0  # fortran line counter
    do_indent = True
    use_same_line = False
    stream = InputStream(infile, orig_filename)
    skip_blank = False
    in_format_off_block = False

    while 1:
        f_line, comments, lines = stream.next_fortran_line()

        if not lines:
            break

        indent = [0] * len(lines)
        nfl += 1
        orig_lines = lines

        comment_lines = format_comments(lines, comments)
        auto_align, auto_format, in_format_off_block = parse_fprettify_directives(
            lines, comment_lines, in_format_off_block, orig_filename, stream.line_nr)
        f_line, lines, is_omp, is_omp_conditional = preprocess_omp(
            f_line, lines)

        lines, do_format, use_indent, is_blank = preprocess_line(
            f_line, lines, comments, orig_filename, stream.line_nr)
        if is_blank and skip_blank:
            continue
        if not do_format and do_indent:
            if use_indent:
                assert len(indent) == 1
                # inherit indent from previous line
                indent[0] = indenter.get_fline_indent()
            else:
                do_indent = False
        else:

            if not auto_align:
                manual_lines_indent = get_manual_alignment(lines)
            else:
                manual_lines_indent = []

            lines, pre_ampersand, ampersand_sep = remove_pre_ampersands(
                lines, orig_filename, stream.line_nr)
            linebreak_pos = get_linebreak_pos(lines)

            f_line = f_line.strip(' ')

            lines = format_single_fline(
                f_line, whitespace, linebreak_pos, ampersand_sep,
                orig_filename, stream.line_nr, auto_format)

            lines = append_comments(lines, comment_lines)

            # target indent for next line
            rel_indent = req_indents[nfl] if nfl < len(req_indents) else 0

            indenter.process_lines_of_fline(
                f_line, lines, rel_indent, indent_size,
                stream.line_nr, manual_lines_indent)
            indent = indenter.get_lines_indent()

            lines, indent = prepend_ampersands(lines, indent, pre_ampersand)

        lines = remove_trailing_whitespace(lines)

        write_formatted_line(outfile, indent, lines, orig_lines, do_indent,
                             use_same_line, is_omp_conditional, orig_filename, stream.line_nr)

        do_indent, use_same_line = pass_defaults_to_next_line(f_line)

        # rm subsequent blank lines
        skip_blank = EMPTY_RE.search(
            f_line) and not any(comments) and not is_omp


def format_comments(lines, comments):
    comments_ftd = []
    for line, comment in zip(lines, comments):
        has_comment = bool(comment.strip())
        sep = has_comment and not comment.strip() == line.strip()
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

    is_omp = OMP_RE.search(f_line)
    is_omp_conditional = bool(is_omp and not OMP_DIR_RE.search(f_line))
    if is_omp_conditional:
        f_line = OMP_RE.sub('  ', f_line, count=1)
        lines = [OMP_RE.sub('  ', l, count=1) for l in lines]

    return [f_line, lines, is_omp, is_omp_conditional]


def preprocess_line(f_line, lines, comments, filename, line_nr):
    """preprocess lines: identification and formatting of special cases"""
    is_blank = False
    use_indent = False
    do_format = False

    if OMP_DIR_RE.search(f_line):
        # move '!$OMP' to line start, otherwise don't format omp directives
        lines = ['!$OMP' + (len(l) - len(l.lstrip())) *
                 ' ' + OMP_DIR_RE.sub('', l, count=1) for l in lines]
    elif lines[0].startswith('#'):  # preprocessor macros
        if len(lines) != 1:
            raise FprettifyInternalException(
                "Continuation lines for preprocessor statement", filename, line_nr)
    elif EMPTY_RE.search(f_line):  # empty lines including comment lines
        if len(lines) != 1:
            raise FprettifyInternalException(
                "Continuation lines for comment lines", filename, line_nr)
        if any(comments):
            if not lines[0].startswith('!'):
                # indent comment lines only if they were not indented before.
                use_indent = True
        else:
            is_blank = True
        lines = [l.strip(' ') for l in lines]
    else:
        do_format = True

    return [lines, do_format, use_indent, is_blank]


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
            lines[pos] = amp_insert + line

    return [lines, indent]


def append_comments(lines, comment_lines):
    """append comments to lines"""
    for pos, (line, comment) in enumerate(zip(lines, comment_lines)):
        if pos < len(lines) - 1:
            has_nl = True  # has next line
        else:
            has_nl = not re.search(EOL_SC, line)
        lines[pos] = lines[pos].rstrip(' ') + comment + '\n' * has_nl

    return lines


def get_linebreak_pos(lines):
    """extract linebreak positions in Fortran line from lines"""
    linebreak_pos = []
    for line in lines:
        found = None
        for char_pos, _ in CharFilter(enumerate(line)):
            if re.match(LINEBREAK_STR, line[char_pos:], RE_FLAGS):
                found = char_pos
        if found:
            linebreak_pos.append(found)
        elif line.lstrip(' ')[0] in ['!','#']:
            linebreak_pos.append(0)

    linebreak_pos = [sum(linebreak_pos[0:_ + 1]) -
                     1 for _ in range(0, len(linebreak_pos))]

    return linebreak_pos


def remove_pre_ampersands(lines, filename, line_nr):
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
            try:
                # amount of whitespace before ampersand of previous line:
                sep = len(re.search(r'(\s*)&[\s]*(?:!.*)?$',
                                    lines[pos - 1]).group(1))
            except AttributeError:
                raise FprettifyParseException(
                    "Bad continuation line format", filename, line_nr)

            ampersand_sep.append(sep)
        else:
            pre_ampersand.append('')
            if pos > 0:
                # use default 1 whitespace character before ampersand
                ampersand_sep.append(1)

    lines = [l.strip(' ').strip('&') for l in lines]

    return [lines, pre_ampersand, ampersand_sep]


def get_manual_alignment(lines):
    """extract manual indents for line continuations from line"""
    manual_lines_indent = [
        len(l) - len(l.lstrip(' ').lstrip('&')) for l in lines]
    manual_lines_indent = [ind - manual_lines_indent[0]
                           for ind in manual_lines_indent]
    return manual_lines_indent


def write_formatted_line(outfile, indent, lines, orig_lines, do_indent, use_same_line, is_omp_conditional, filename, line_nr):
    """Write reformatted line to file"""
    for ind, line, orig_line in zip(indent, lines, orig_lines):

        # get actual line length excluding comment:
        line_length = 0
        for line_length, _ in CharFilter(enumerate(line)):
            pass
        line_length += 1

        if do_indent:
            ind_use = ind
        else:
            if use_same_line:
                ind_use = 1
            else:
                ind_use = 0

        if line.lstrip().startswith('#'):
            ind_use = 0

        if ind_use + line_length <= 133:  # 132 plus 1 newline char
            outfile.write('!$' * is_omp_conditional +
                          ' ' * (ind_use - 2 * is_omp_conditional +
                                 len(line) - len(line.lstrip(' '))) +
                          line.lstrip(' '))
        elif line_length <= 133:
            outfile.write('!$' * is_omp_conditional + ' ' *
                          (133 - 2 * is_omp_conditional -
                           len(line.lstrip(' '))) + line.lstrip(' '))

            log_message(LINESPLIT_MESSAGE, "warning",
                        filename, line_nr)
        else:
            outfile.write(orig_line)
            log_message(LINESPLIT_MESSAGE, "warning",
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

    parser = argparse.ArgumentParser(prog=argv[0],
                                     description='Auto-format modern Fortran source files.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--indent", type=int, default=3,
                        help="relative indentation width")
    parser.add_argument("-w", "--whitespace", type=int,
                        choices=range(0, 4), default=2, help="Amount of whitespace")
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
    parser.add_argument('-f', '--fortran', type=str, action='append', default=[],
                        help="Overrides default fortran extensions recognized by --recursive. Repeat this option to specify more than one extension.")
    parser.add_argument('--version', action='version',
                        version='%(prog)s 0.3.2')

    args = parser.parse_args(argv[1:])

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
            for dirpath, _, files in os.walk(directory):
                for ffile in [os.path.join(dirpath, f) for f in files if any(f.endswith(_) for _ in ext)]:
                    filenames.append(ffile)

        for filename in filenames:
            stdout = args.stdout or directory == '-'

            if args.debug:
                level = logging.DEBUG
            elif args.silent:
                level = logging.CRITICAL
            else:
                level = logging.WARNING

            set_fprettify_logger(level)

            try:
                reformat_inplace(filename,
                                 stdout=stdout,
                                 indent_size=args.indent,
                                 whitespace=args.whitespace)
            except FprettifyException as e:
                log_exception(e, "Fatal error occured")
