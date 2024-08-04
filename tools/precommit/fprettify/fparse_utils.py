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

"""This is a collection of Fortran parsing utilities."""

import re
from collections import deque

RE_FLAGS = re.IGNORECASE | re.UNICODE

# FIXME bad ass regex!
VAR_DECL_RE = re.compile(
    r"^ *(?P<type>integer(?: *\* *[0-9]+)?|logical|character(?: *\* *[0-9]+)?|real(?: *\* *[0-9]+)?|complex(?: *\* *[0-9]+)?|type) *(?P<parameters>\((?:[^()]+|\((?:[^()]+|\([^()]*\))*\))*\))? *(?P<attributes>(?: *, *[a-zA-Z_0-9]+(?: *\((?:[^()]+|\((?:[^()]+|\([^()]*\))*\))*\))?)+)? *(?P<dpnt>::)?(?P<vars>[^\n]+)\n?", RE_FLAGS)

OMP_COND_RE = re.compile(r"^\s*(!\$ )", RE_FLAGS)

# supported preprocessors
FYPP_LINE_STR = r"^(#!|#:|\$:|@:)"
FYPP_WITHOUT_PREPRO_STR = r"^(#!|\$:|@:)"
CPP_STR = r"^#[^!:{}]"
COMMENT_LINE_STR = r"^!"
FYPP_OPEN_STR = r"(#{|\${|@{)"
FYPP_CLOSE_STR = r"(}#|}\$|}@)"
NOTFORTRAN_LINE_RE = re.compile(r"("+FYPP_LINE_STR+r"|"+CPP_STR+r"|"+COMMENT_LINE_STR+r")", RE_FLAGS)
NOTFORTRAN_FYPP_LINE_RE = re.compile(r"("+CPP_STR+r"|"+COMMENT_LINE_STR+r")", RE_FLAGS)
FYPP_LINE_RE = re.compile(FYPP_LINE_STR, RE_FLAGS)
FYPP_WITHOUT_PREPRO_RE = re.compile(FYPP_WITHOUT_PREPRO_STR, RE_FLAGS)
FYPP_OPEN_RE = re.compile(FYPP_OPEN_STR, RE_FLAGS)
FYPP_CLOSE_RE = re.compile(FYPP_CLOSE_STR, RE_FLAGS)

STR_OPEN_RE = re.compile(r"("+FYPP_OPEN_STR+r"|"+r"'|\"|!)", RE_FLAGS)
CPP_RE = re.compile(CPP_STR, RE_FLAGS)

class fline_parser(object):
    def __init__(self):
        pass
    def search(self, line):
        pass

class parser_re(fline_parser):
    def __init__(self, regex, spec=True):
        self._re = regex
        self.spec = spec

    def search(self, line):
        return self._re.search(line)

    def split(self, line):
        return self._re.split(line)

class FprettifyException(Exception):
    """Base class for all custom exceptions"""

    def __init__(self, msg, filename, line_nr):
        super(FprettifyException, self).__init__(msg)
        self.filename = filename
        self.line_nr = line_nr


class FprettifyParseException(FprettifyException):
    """Exception for unparseable Fortran code (user's fault)."""

    pass


class FprettifyInternalException(FprettifyException):
    """Exception for potential internal errors (fixme's)."""

    pass


class CharFilter(object):
    """
    An iterator to wrap the iterator returned by `enumerate(string)`
    and ignore comments and characters inside strings
    """

    def __init__(self, string, filter_comments=True, filter_strings=True,
                 filter_fypp=True):
        self._content = string
        self._it = enumerate(self._content)
        self._instring = ''
        self._infypp = False
        self._incomment = ''
        self._instring = ''
        self._filter_comments = filter_comments
        self._filter_strings = filter_strings
        if filter_fypp:
            self._notfortran_re = NOTFORTRAN_LINE_RE
        else:
            self._notfortran_re = NOTFORTRAN_FYPP_LINE_RE


    def update(self, string, filter_comments=True, filter_strings=True,
               filter_fypp=True):
        self._content = string
        self._it = enumerate(self._content)
        self._filter_comments = filter_comments
        self._filter_strings = filter_strings
        if filter_fypp:
            self._notfortran_re = NOTFORTRAN_LINE_RE
        else:
            self._notfortran_re = NOTFORTRAN_FYPP_LINE_RE

    def __iter__(self):
        return self

    def __next__(self):

        pos, char = next(self._it)

        char2 = self._content[pos:pos+2]

        if not self._instring:
            if not self._incomment:
                if FYPP_OPEN_RE.search(char2):
                    self._instring = char2
                    self._infypp = True
                elif (self._notfortran_re.search(char2)):
                    self._incomment = char
                elif char in ['"', "'"]:
                    self._instring = char
        else:
            if self._infypp:
                if FYPP_CLOSE_RE.search(char2):
                    self._instring = ''
                    self._infypp = False
                    if self._filter_strings:
                        self.__next__()
                        return self.__next__()

            elif char in ['"', "'"]:
                if self._instring == char:
                    self._instring = ''
                    if self._filter_strings:
                        return self.__next__()

        if self._filter_comments:
            if self._incomment:
                raise StopIteration

        if self._filter_strings:
            if self._instring:
                return self.__next__()

        return (pos, char)

    def filter_all(self):
        filtered_str = ''
        for pos, char in self:
            filtered_str += char
        return filtered_str

    def instring(self):
        return self._instring

class InputStream(object):
    """Class to read logical Fortran lines from a Fortran file."""

    def __init__(self, infile, filter_fypp=True, orig_filename=None):
        if not orig_filename:
            orig_filename = infile.name
        self.line_buffer = deque([])
        self.infile = infile
        self.line_nr = 0
        self.filename = orig_filename
        self.endpos = deque([])
        self.what_omp = deque([])
        if filter_fypp:
            self.notfortran_re = NOTFORTRAN_LINE_RE
        else:
            self.notfortran_re = NOTFORTRAN_FYPP_LINE_RE

    def next_fortran_line(self):
        """Reads a group of connected lines (connected with &, separated by newline or semicolon)
        returns a touple with the joined line, and a list with the original lines.
        Doesn't support multiline character constants!
        """
        joined_line = ""
        comments = []
        lines = []
        continuation = 0
        fypp_cont = 0
        instring = ''

        string_iter = CharFilter('')
        fypp_cont = 0
        while 1:
            if not self.line_buffer:
                line = self.infile.readline().replace("\t", 8 * " ")
                self.line_nr += 1
                # convert OMP-conditional fortran statements into normal fortran statements
                # but remember to convert them back

                what_omp = OMP_COND_RE.search(line)

                if what_omp:
                    what_omp = what_omp.group(1)
                else:
                    what_omp = ''

                if what_omp:
                    line = line.replace(what_omp, '', 1)
                line_start = 0

                pos = -1

                # multiline string: prepend line continuation with '&'
                if string_iter.instring() and not line.lstrip().startswith('&'):
                    line = '&' + line

                # update instead of CharFilter(line) to account for multiline strings
                string_iter.update(line)
                for pos, char in string_iter:
                    if char == ';' or pos + 1 == len(line):
                        self.endpos.append(pos - line_start)
                        self.line_buffer.append(line[line_start:pos + 1])
                        self.what_omp.append(what_omp)
                        what_omp = ''
                        line_start = pos + 1

                if pos + 1 < len(line):
                   if fypp_cont:
                       self.endpos.append(-1)
                       self.line_buffer.append(line)
                       self.what_omp.append(what_omp)
                   else:
                       for pos_add, char in CharFilter(line[pos+1:], filter_comments=False):
                           char2 = line[pos+1+pos_add:pos+3+pos_add]
                           if self.notfortran_re.search(char2):
                               self.endpos.append(pos + pos_add - line_start)
                               self.line_buffer.append(line[line_start:])
                               self.what_omp.append(what_omp)
                               break

                if not self.line_buffer:
                    self.endpos.append(len(line))
                    self.line_buffer.append(line)
                    self.what_omp.append('')


            line = self.line_buffer.popleft()
            endpos = self.endpos.popleft()
            what_omp = self.what_omp.popleft()

            if not line:
                break

            lines.append(what_omp + line)

            line_core = line[:endpos + 1]

            if self.notfortran_re.search(line[endpos+1:endpos+3]) or fypp_cont:
                line_comments = line[endpos + 1:]
            else:
                line_comments = ''

            if line_core:
                newline = (line_core[-1] == '\n')
            else:
                newline = False

            line_core = line_core.strip()

            if line_core and not NOTFORTRAN_LINE_RE.search(line_core):
                continuation = 0
            if line_core.endswith('&'):
                continuation = 1

            if line_comments:
                if (FYPP_LINE_RE.search(line[endpos+1:endpos+3]) or fypp_cont) and line_comments.strip()[-1] == '&':
                    fypp_cont = 1
                else:
                    fypp_cont = 0

            line_core = line_core.strip('&')

            comments.append(line_comments.rstrip('\n'))
            if joined_line.strip():
                joined_line = joined_line.rstrip(
                    '\n') + line_core + '\n' * newline
            else:
                joined_line = what_omp + line_core + '\n' * newline

            if not (continuation or fypp_cont):
                break

        return (joined_line, comments, lines)
