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

"""Dynamically create tests based on examples in examples/before."""
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

import sys
import os
import unittest
import hashlib
import logging
import io
import re
import difflib
import subprocess
import inspect

sys.stderr = io.TextIOWrapper(
    sys.stderr.detach(), encoding='UTF-8', line_buffering=True)

import fprettify
from fprettify.fparse_utils import FprettifyParseException, FprettifyInternalException


def joinpath(path1, path2):
    return os.path.normpath(os.path.join(path1, path2))

MYPATH = os.path.dirname(os.path.abspath(
    inspect.getfile(inspect.currentframe())))

BEFORE_DIR = joinpath(MYPATH, r'../../fortran_tests/before/')
AFTER_DIR = joinpath(MYPATH, r'../../fortran_tests/after/')
RESULT_DIR = joinpath(MYPATH, r'../../fortran_tests/test_results/')
RESULT_FILE = joinpath(RESULT_DIR, r'expected_results')
FAILED_FILE = joinpath(RESULT_DIR, r'failed_results')

RUNSCRIPT = joinpath(MYPATH, r"../../fprettify.py")

fprettify.set_fprettify_logger(logging.ERROR)


class AlienInvasion(Exception):
    """Should not happen"""
    pass


def eprint(*args, **kwargs):
    """
    Print to stderr - to print output compatible with default unittest output.
    """

    print(*args, file=sys.stderr, flush=True, **kwargs)

class FPrettifyTestCase(unittest.TestCase):
    """
    test class to be recognized by unittest.
    """

    def shortDescription(self):
        """don't print doc string of testmethod"""
        return None

    def setUp(self):
        """
        setUp to be recognized by unittest.
        We have large files to compare, raise the limit
        """
        self.maxDiff = None

    @classmethod
    def setUpClass(cls):
        """
        setUpClass to be recognized by unittest.
        """

        cls.n_success = 0
        cls.n_parsefail = 0
        cls.n_internalfail = 0
        cls.n_unexpectedfail = 0

        eprint("-" * 70)
        eprint("recognized Fortran files")
        eprint(", ".join(fprettify.FORTRAN_EXTENSIONS))
        eprint("-" * 70)
        eprint("Testing with Fortran files in " + BEFORE_DIR)
        eprint("Writing formatted Fortran files to " + AFTER_DIR)
        eprint("Storing expected results in " + RESULT_FILE)
        eprint("Storing failed results in " + FAILED_FILE)
        eprint("-" * 70)

    @classmethod
    def tearDownClass(cls):
        """
        tearDownClass to be recognized by unittest. Used for test summary
        output.
        """
        if cls.n_parsefail + cls.n_internalfail > 0:
            format = "{:<20}{:<6}"
            eprint('\n' + "=" * 70)
            eprint("IGNORED errors: invalid or old Fortran")
            eprint("-" * 70)
            eprint(format.format("parse errors: ", cls.n_parsefail))
            eprint(format.format("internal errors: ", cls.n_internalfail))

    @staticmethod
    def write_result(filename, content, sep_str):  # pragma: no cover
        with io.open(filename, 'a', encoding='utf-8') as outfile:
            outfile.write(sep_str.join(content) + '\n')

    def test_whitespace(self):
        """simple test for whitespace formatting options -w in [0, 1, 2]"""
        instring = "(/-a-b-(a+b-c)/(-c)*d**e,f[1]%v/)"
        outstring_exp = ["(/-a-b-(a+b-c)/(-c)*d**e,f[1]%v/)",
                         "(/-a-b-(a+b-c)/(-c)*d**e, f[1]%v/)",
                         "(/-a - b - (a + b - c)/(-c)*d**e, f[1]%v/)",
                         "(/-a - b - (a + b - c) / (-c) * d**e, f[1]%v/)"]

        outstring = []
        for w, out in zip(range(0, 4), outstring_exp):
            args = ['-w', str(w)]
            self.assert_fprettify_result(args, instring, out)

    def test_type_selector(self):
        """test for whitespace formatting option -w 4"""
        instring = "A%component=func(mytype%a,mytype%abc+mytype%abcd)"
        outstring_exp = "A % component = func(mytype % a, mytype % abc + mytype % abcd)"

        self.assert_fprettify_result(['-w 4'], instring, outstring_exp)

    def test_indent(self):
        """simple test for indent options -i in [0, 3, 4]"""

        indents = [0, 3, 4]

        instring = "iF(teSt)ThEn\nCaLl subr(a,b,&\nc,(/d,&\ne,f/))\nEnD iF"
        outstring_exp = [
            "iF (teSt) ThEn\n" +
            " " * ind + "CaLl subr(a, b, &\n" +
            " " * (10 + ind) + "c, (/d, &\n" +
            " " * (15 + ind) + "e, f/))\nEnD iF"
            for ind in indents
        ]

        for ind, out in zip(indents, outstring_exp):
            args = ['-i', str(ind)]
            self.assert_fprettify_result(args, instring, out)

    def test_nested(self):
        """test correct indentation of nested loops"""
        instring = ("integer :: i,j\ndo i=1,2\ndo j=1,3\n"
                    "print*,i,j,i*j\nend do\nend do")
        outstring_exp_default = ("integer :: i, j\ndo i = 1, 2\ndo j = 1, 3\n"
                                 "   print *, i, j, i*j\nend do\nend do")
        outstring_exp_strict = ("integer :: i, j\ndo i = 1, 2\n   do j = 1, 3\n"
                                "      print *, i, j, i*j\n   end do\nend do")

        self.assert_fprettify_result([], instring, outstring_exp_default)
        self.assert_fprettify_result(['--strict-indent'], instring, outstring_exp_strict)

    def test_reset_indent(self):
        """test of reset indentation at file start"""
        instring = ("integer :: i,j\ndo i=1,2\ndo j=1,3\n"
                    "print*,i,j,i*j\nend do\nend do",
                    "   module a\ninteger :: 1\n")
        outstring = ("integer :: i, j\ndo i = 1, 2\ndo j = 1, 3\n"
                     "   print *, i, j, i*j\nend do\nend do",
                     "module a\n   integer :: 1")

        for ind, out in zip(instring, outstring):
            self.assert_fprettify_result([],ind, out)

    def test_disable(self):
        """test disabling indentation and/or whitespace formatting"""
        instring = ("if(&\nl==111)&\n then\n   do m   =1,  2\n A=&\nB+C\n    end  do;   endif")
        outstring_exp_default = ("if ( &\n   l == 111) &\n   then\n   do m = 1, 2\n"
                                 "      A = &\n         B + C\n   end do; end if")
        outstring_exp_nowhitespace = ("if(&\n   l==111)&\n   then\n   do m   =1,  2\n"
                                      "      A=&\n         B+C\n   end  do; endif")
        outstring_exp_noindent = ("if ( &\nl == 111) &\n then\n   do m = 1, 2\n"
                                  " A = &\nB + C\n    end do;   end if")

        self.assert_fprettify_result([], instring, outstring_exp_default)
        self.assert_fprettify_result(['--disable-whitespace'], instring, outstring_exp_nowhitespace)
        self.assert_fprettify_result(['--disable-indent'], instring, outstring_exp_noindent)
        self.assert_fprettify_result(['--disable-indent', '--disable-whitespace'], instring, instring)

    def test_comments(self):
        """test options related to comments"""
        instring = ("TYPE mytype\n!  c1\n  !c2\n   INTEGER :: a   !  c3\n"
                    "   REAL :: b, &   ! c4\n! c5\n                  ! c6\n"
                    "           d      ! c7\nEND TYPE  ! c8")
        outstring_exp_default = ("TYPE mytype\n!  c1\n   !c2\n   INTEGER :: a   !  c3\n"
                                 "   REAL :: b, &   ! c4\n           ! c5\n           ! c6\n"
                                 "           d      ! c7\nEND TYPE  ! c8")
        outstring_exp_strip = ("TYPE mytype\n!  c1\n   !c2\n   INTEGER :: a !  c3\n"
                               "   REAL :: b, & ! c4\n           ! c5\n           ! c6\n"
                               "           d ! c7\nEND TYPE ! c8")

        self.assert_fprettify_result([], instring, outstring_exp_default)
        self.assert_fprettify_result(['--strip-comments'], instring, outstring_exp_strip)

    def test_directive(self):
        """
        test deactivate directives '!&' (inline) and '!&<', '!&>' (block)
        and manual alignment (continuation line prefixed with '&')
        """

        # manual alignment
        instring = "align_me = [ -1,  10,0,  &\n    &     0,1000 ,  0,&\n            &0 , -1,  1]"
        outstring_exp = "align_me = [-1, 10, 0,  &\n    &     0, 1000, 0,&\n            &0, -1, 1]"
        self.assert_fprettify_result([], instring, outstring_exp)

        # inline deactivate
        instring2 = '\n'.join(_ + ' !&' for _ in instring.splitlines())
        outstring_exp = instring2
        self.assert_fprettify_result([], instring2, outstring_exp)

        # block deactivate
        instring3 = '!&<\n' + instring + '\n!&>'
        outstring_exp = instring3
        self.assert_fprettify_result([], instring3, outstring_exp)

    def assert_fprettify_result(self, args, instring, outstring_exp):
        """
        assert that result of calling fprettify with args on instring gives
        outstring_exp
        """
        args.insert(0, RUNSCRIPT)
        p1 = subprocess.Popen(
            args, stdout=subprocess.PIPE, stdin=subprocess.PIPE)
        outstring = p1.communicate(instring.encode(
            'UTF-8'))[0].decode('UTF-8').rstrip()
        self.assertEqual(outstring_exp.rstrip(), outstring)

    def test_io(self):
        """simple test for io (file inplace, stdin & stdout)"""

        # io and unicode
        outstring = []
        instring = "CALL  alien_invasion( 👽 )"
        outstring_exp = "CALL alien_invasion(👽)"

        alien_file = "alien_invasion.f90"
        if os.path.isfile(alien_file):
            raise AlienInvasion(
                "remove file alien_invasion.f90")  # pragma: no cover

        try:
            with io.open(alien_file, 'w', encoding='utf-8') as infile:
                infile.write(instring)

            # testing stdin --> stdout
            p1 = subprocess.Popen(RUNSCRIPT,
                                  stdout=subprocess.PIPE, stdin=subprocess.PIPE)
            outstring.append(p1.communicate(
                instring.encode('UTF-8'))[0].decode('UTF-8'))

            # testing file --> stdout
            p1 = subprocess.Popen([RUNSCRIPT, alien_file, '--stdout'],
                                  stdout=subprocess.PIPE)
            outstring.append(p1.communicate(
                instring.encode('UTF-8')[0])[0].decode('UTF-8'))

            # testing file --> file (inplace)
            p1 = subprocess.Popen([RUNSCRIPT, alien_file])
            p1.wait()

            with io.open(alien_file, 'r', encoding='utf-8') as infile:
                outstring.append(infile.read())

            for outstr in outstring:
                self.assertEqual(outstring_exp, outstr.strip())
        except:  # pragma: no cover
            if os.path.isfile(alien_file):
                os.remove(alien_file)
            raise
        else:
            os.remove(alien_file)

    def test_multi_alias(self):
        """test for issue #11 (multiple alias and alignment)"""
        instring="use A,only:B=>C,&\nD=>E"
        outstring="use A, only: B => C, &\n             D => E"
        self.assert_fprettify_result([], instring, outstring)

    def test_use(self):
        """test for alignment of use statements"""
        instring1="use A,only:B,C,&\nD,E"
        instring2="use A,only:&\nB,C,D,E"
        outstring1="use A, only: B, C, &\n             D, E"
        outstring2="use A, only: &\n   B, C, D, E"
        self.assert_fprettify_result([], instring1, outstring1)
        self.assert_fprettify_result([], instring2, outstring2)

    def test_wrongkind(self):
        """test whitespacing of deprecated kind definition"""
        instring = ["REAL*8 :: r, f  !  some reals",
                    "REAL * 8 :: r, f  !  some reals",
                    "INTEGER * 4 :: c, i  !  some integers",
                    "INTEGER*4 :: c, i  !  some integers"]
        outstring = ["REAL*8 :: r, f  !  some reals",
                     "REAL*8 :: r, f  !  some reals",
                     "INTEGER*4 :: c, i  !  some integers",
                     "INTEGER*4 :: c, i  !  some integers"]

        for i in range(0, len(instring)):
            self.assert_fprettify_result([], instring[i], outstring[i])

    def test_new_intrinsics(self):
        """test new I/O intrinsics"""
        instring = ["REWIND(12)",
                    "BACKSPACE(13)",
                    "INQUIRE(14)"]
        outstring = ["REWIND (12)",
                     "BACKSPACE (13)",
                     "INQUIRE (14)"]

        for i in range(0, len(instring)):
            self.assert_fprettify_result([], instring[i], outstring[i])

    def test_associate(self):
        """test correct formatting of associate construct"""
        instring = ("associate(a=>b , c  =>d ,e=> f  )\n"
                    "e=a+c\n"
                    "end associate")
        outstring = ("associate (a => b, c => d, e => f)\n"
                    "   e = a + c\n"
                    "end associate")

        self.assert_fprettify_result([], instring, outstring)

    def test_line_length(self):
        """test line length option"""
        instring = ["REAL(KIND=4) :: r,f  !  some reals",
                    "if(   min == max.and.min .eq. thres  )",
                    "INQUIRE(14)"]
        instring_ = "if( min == max.and.min .eq. thres ) one_really_long_function_call_to_hit_the_line_limit(parameter1, parameter2,parameter3,parameter4,parameter5,err) ! this line would be too long"
        outstring = ["REAL(KIND=4) :: r, f  !  some reals",
                     "REAL(KIND=4) :: r,f  !  some reals",
                     "if (min == max .and. min .eq. thres)",
                     "if(   min == max.and.min .eq. thres  )",
                     "INQUIRE (14)",
                     "INQUIRE (14)"]
        outstring_ = ["if( min == max.and.min .eq. thres ) one_really_long_function_call_to_hit_the_line_limit(parameter1, parameter2,parameter3,parameter4,parameter5,err) ! this line would be too long",
                      "if (min == max .and. min .eq. thres) one_really_long_function_call_to_hit_the_line_limit(parameter1, parameter2, parameter3, parameter4, parameter5, err) ! this line would be too long"]

        # test shorter lines first, after all the actual length doesn't matter
        for i in range(0, len(instring)):
            self.assert_fprettify_result(['-S'], instring[i], outstring[2*i])
            self.assert_fprettify_result(['-S', '-l 20'], instring[i], outstring[2*i + 1])
        # now test a long line
        self.assert_fprettify_result(['-S'], instring_, outstring_[0])
        self.assert_fprettify_result(['-S', '-l 0'], instring_, outstring_[1])

    def test_relation_replacement(self):
        """test replacement of relational statements"""
        instring = ["if ( min < max .and. min .lt. thres)",
                    "if (min > max .and. min .gt. thres )",
                    "if (   min == max .and. min .eq. thres  )",
                    "if(min /= max .and. min .ne. thres)",
                    "if(min >= max .and. min .ge. thres )",
                    "if( min <= max .and. min .le. thres)",
                    "'==== heading",
                    "if (vtk%my_rank .eq. 0) write (vtk%filehandle_par, '(\"<DataArray",
                    "'(\"</Collection>\",",
                    "if (abc(1) .lt. -bca .or. &\n qwe .gt. ewq) then"]
        f_outstring = ["if (min .lt. max .and. min .lt. thres)",
                     "if (min .gt. max .and. min .gt. thres)",
                     "if (min .eq. max .and. min .eq. thres)",
                     "if (min .ne. max .and. min .ne. thres)",
                     "if (min .ge. max .and. min .ge. thres)",
                     "if (min .le. max .and. min .le. thres)",
                     "'==== heading",
                     "if (vtk%my_rank .eq. 0) write (vtk%filehandle_par, '(\"<DataArray",
                     "'(\"</Collection>\",",
                     "if (abc(1) .lt. -bca .or. &\n    qwe .gt. ewq) then"]
        c_outstring = ["if (min < max .and. min < thres)",
                     "if (min > max .and. min > thres)",
                     "if (min == max .and. min == thres)",
                     "if (min /= max .and. min /= thres)",
                     "if (min >= max .and. min >= thres)",
                     "if (min <= max .and. min <= thres)",
                     "'==== heading",
                     "if (vtk%my_rank == 0) write (vtk%filehandle_par, '(\"<DataArray",
                      "'(\"</Collection>\",",
                     "if (abc(1) < -bca .or. &\n    qwe > ewq) then"]
        for i in range(0, len(instring)):
            self.assert_fprettify_result(['--enable-replacements', '--c-relations'], instring[i], c_outstring[i])
            self.assert_fprettify_result(['--enable-replacements'], instring[i], f_outstring[i])

    def test_swap_case(self):
        """test replacement of keyword character case"""
        instring = (
            "MODULE exAmple",
            "INTEGER,   PARAMETER :: SELECTED_REAL_KIND = 1*2",
            "INTEGER,   PARAMETER :: dp1 = SELECTED_REAL_KIND ( 15 , 307)",
            'CHARACTER(LEN=*), PARAMETER :: a = "INTEGER,   PARAMETER" // "b"',
            "CHARACTER(LEN=*), PARAMETER :: a = 'INTEGER,   PARAMETER' // 'b'",
            "INTEGER(kind=int64), PARAMETER :: l64 = 2_int64",
            "REAL(kind=real64), PARAMETER :: r64a = 2._real64",
            "REAL(kind=real64), PARAMETER :: r64b = 2.0_real64",
            "REAL(kind=real64), PARAMETER :: r64c = .0_real64",
            "REAL(kind=real64), PARAMETER :: r64a = 2.e3_real64",
            "REAL(kind=real64), PARAMETER :: r64b = 2.0e3_real64",
            "REAL(kind=real64), PARAMETER :: r64c = .0e3_real64",
            "REAL, PARAMETER :: r32 = 2.e3",
            "REAL, PARAMETER :: r32 = 2.0d3",
            "REAL, PARAMETER :: r32 = .2e3",
            "USE ISO_FORTRAN_ENV, ONLY: int64",
            "INTEGER, INTENT(IN) :: r, i, j, k",
            "IF (l.EQ.2) l=MAX  (l64, 2_int64)",
            "PURE SUBROUTINE mypure()"
            )
        outstring = (
            "module exAmple",
            "integer, parameter :: SELECTED_REAL_KIND = 1*2",
            "integer, parameter :: dp1 = selected_real_kind(15, 307)",
            'character(LEN=*), parameter :: a = "INTEGER,   PARAMETER"//"b"',
            "character(LEN=*), parameter :: a = 'INTEGER,   PARAMETER'//'b'",
            "integer(kind=INT64), parameter :: l64 = 2_INT64",
            "real(kind=REAL64), parameter :: r64a = 2._REAL64",
            "real(kind=REAL64), parameter :: r64b = 2.0_REAL64",
            "real(kind=REAL64), parameter :: r64c = .0_REAL64",
            "real(kind=REAL64), parameter :: r64a = 2.E3_REAL64",
            "real(kind=REAL64), parameter :: r64b = 2.0E3_REAL64",
            "real(kind=REAL64), parameter :: r64c = .0E3_REAL64",
            "real, parameter :: r32 = 2.E3",
            "real, parameter :: r32 = 2.0D3",
            "real, parameter :: r32 = .2E3",
            "use iso_fortran_env, only: INT64",
            "integer, intent(IN) :: r, i, j, k",
            "if (l .eq. 2) l = max(l64, 2_INT64)",
            "pure subroutine mypure()"
            )
        for i in range(len(instring)):
            self.assert_fprettify_result(['--case', '1', '1', '1', '2'],
                                         instring[i], outstring[i])

    def test_do(self):
        """test correct parsing of do statement"""
        instring = "do = 1\nb = 2"

        self.assert_fprettify_result([], instring, instring)

    def test_omp(self):
        """test formatting of omp directives"""
        instring = ("PROGRAM test_omp\n"
                    " !$OMP    PARALLEL DO\n"
                    "b=4\n"
                    "!$a=b\n"
                    "!$  a=b\n"
                    "   !$    c=b\n"
                    "!$acc parallel loop\n"
                    "!$OMP END  PARALLEL DO\n"
                    "END PROGRAM")
        outstring = ("PROGRAM test_omp\n"
                     "   !$OMP    PARALLEL DO\n"
                     "   b = 4\n"
                     "!$a=b\n"
                     "!$ a = b\n"
                     "!$ c = b\n"
                     "!$acc parallel loop\n"
                     "!$OMP END  PARALLEL DO\n"
                     "END PROGRAM")

        self.assert_fprettify_result([], instring, outstring)

    def test_ford(self):
        """test formatting of ford comments"""
        instring =  ("   a = b\n"
                     "     !!  ford docu\n"
                     "b=c\n"
                     "  !! ford docu\n"
                     "subroutine test(a,b,&\n"
                     "  !! ford docu\n"
                     "  c, d, e)"
                     )
        outstring = ("   a = b\n"
                     "     !!  ford docu\n"
                     "   b = c\n"
                     "  !! ford docu\n"
                     "   subroutine test(a, b, &\n"
                     "  !! ford docu\n"
                     "                   c, d, e)"
                     )

        self.assert_fprettify_result([], instring, outstring)

    def test_plusminus(self):
        """test corner cases of +/-"""
        instring = "val_1d-1-1.0e-9-2.0d-08+.2e-1-val_2d-3.e-12+4"
        outstring = "val_1d - 1 - 1.0e-9 - 2.0d-08 + .2e-1 - val_2d - 3.e-12 + 4"
        self.assert_fprettify_result([], instring, outstring)

    def test_fypp(self):
        """test formatting of fypp statements"""

        instring = []
        outstring = []

        instring += [
"""
#:if DEBUG>  0
print *, "hola"
if( .not. (${cond}$) ) then
#:if  ASSERT(cond)
print *, "Assert failed!"
#:endif
error stop
end if
#:endif
"""
]

        outstring += [
"""
#:if DEBUG>  0
   print *, "hola"
   if (.not. (${cond}$)) then
      #:if  ASSERT(cond)
         print *, "Assert failed!"
      #:endif
      error stop
   end if
#:endif
"""
]

        instring += [
"""
if  (.not. (${cond}$)) then
   #:for element in list
   print *, "Element is in list!"
 #:endfor
   error stop
end if
"""
]

        outstring += [
"""
if (.not. (${cond}$)) then
   #:for element in list
      print *, "Element is in list!"
   #:endfor
   error stop
end if
"""
]

        instring += [
"""
#:if aa > 1
print  *, "Number is more than 1"
if (condition) then
     #:def something
   print *, "Added Definition!"
   #:enddef
end if
#:endif
"""
]

        outstring += [
"""
#:if aa > 1
   print *, "Number is more than 1"
   if (condition) then
      #:def something
         print *, "Added Definition!"
      #:enddef
   end if
#:endif
"""
]

        instring += [
"""
#:def DEBUG_CODE( code)
  #:if DEBUG > 0
    $:code
  #:endif
#:enddef DEBUG_CODE
"""
]

        outstring += [
"""
#:def DEBUG_CODE( code)
   #:if DEBUG > 0
      $:code
   #:endif
#:enddef DEBUG_CODE
"""
]


        instring += [
"""
#:block DEBUG_CODE
  if (a <b) then
    print *, "DEBUG: a is less than b"
  end if
#:endblock  DEBUG_CODE
"""
]

        outstring += [
"""
#:block DEBUG_CODE
   if (a < b) then
      print *, "DEBUG: a is less than b"
   end if
#:endblock  DEBUG_CODE
"""
]

        instring += [
"""
#:call DEBUG_CODE
  if (a < b) then
    print *, "DEBUG: a is less than b"
  end if
#:endcall DEBUG_CODE
"""
]

        outstring += [
"""
#:call DEBUG_CODE
   if (a < b) then
      print *, "DEBUG: a is less than b"
   end if
#:endcall DEBUG_CODE
"""
]

        instring += [
"""
#:if DEBUG > 0
print *, "hola"
if (.not. (${cond}$)) then
   #:mute
   print *, "Muted"
   #:endmute
   error stop
end if
#:endif
"""
]

        outstring += [
"""
#:if DEBUG > 0
   print *, "hola"
   if (.not. (${cond}$)) then
      #:mute
         print *, "Muted"
      #:endmute
      error stop
   end if
#:endif
"""
]

        instring += [
"""
program try
#:def mydef
a = &
#:if dothat
b + &
#:else
c + &
#:endif
d
#:enddef
end program
"""
]

        outstring += [
"""
program try
   #:def mydef
      a = &
#:if dothat
         b + &
#:else
         c + &
#:endif
         d
   #:enddef
end program
"""
]

        instring += [
"""
#:if worktype
      ${worktype}$, &
#:else
      ${type}$, &
#:endif
         DIMENSION(${arr_exp}$), &
         POINTER :: work
"""
]

        outstring += [
"""
#:if worktype
${worktype}$, &
#:else
   ${type}$, &
#:endif
   DIMENSION(${arr_exp}$), &
   POINTER :: work
"""
]



        for instr, outstr in zip(instring, outstring):
            self.assert_fprettify_result([], instr, outstr)

    def test_mod(self):
        """test indentation of module / program"""
        instring_mod = "module my_module\nintrinsic none\ncontains\nfunction my_func()\nend\nend module"
        instring_prog = "program my_program\nintrinsic none\ncontains\nfunction my_func()\nend\nend program"

        outstring_mod = "module my_module\n   intrinsic none\ncontains\n   function my_func()\n   end\nend module"
        outstring_mod_disable = "module my_module\nintrinsic none\ncontains\nfunction my_func()\nend\nend module"

        outstring_prog = "program my_program\n   intrinsic none\ncontains\n   function my_func()\n   end\nend program"
        outstring_prog_disable = "program my_program\nintrinsic none\ncontains\nfunction my_func()\nend\nend program"

        self.assert_fprettify_result([], instring_mod, outstring_mod)
        self.assert_fprettify_result([], instring_prog, outstring_prog)

        self.assert_fprettify_result(['--disable-indent-mod'], instring_mod, outstring_mod_disable)
        self.assert_fprettify_result(['--disable-indent-mod'], instring_prog, outstring_prog_disable)

    def test_decl(self):
        """test formatting of declarations"""
        instring_1 = "integer    ::     a"
        instring_2 = "integer, dimension(:)    ::     a"
        outstring_1 = "integer :: a"
        outstring_2 = "integer, dimension(:) :: a"
        outstring_2_min = "integer, dimension(:)::a"

        self.assert_fprettify_result([], instring_1, instring_1)
        self.assert_fprettify_result([], instring_2, instring_2)
        self.assert_fprettify_result(['--enable-decl'], instring_1, outstring_1)
        self.assert_fprettify_result(['--enable-decl'], instring_2, outstring_2)
        self.assert_fprettify_result(['--enable-decl', '--whitespace-decl=0'], instring_2, outstring_2_min)

    def test_statement_label(self):
        instring = "1003  FORMAT(2(1x, i4), 5x, '-', 5x, '-', 3x, '-', 5x, '-', 5x, '-', 8x, '-', 3x, &\n    1p, 2(1x, d10.3))"
        outstring = "1003 FORMAT(2(1x, i4), 5x, '-', 5x, '-', 3x, '-', 5x, '-', 5x, '-', 8x, '-', 3x, &\n            1p, 2(1x, d10.3))"
        self.assert_fprettify_result([], instring, outstring)

    def test_multiline_str(self):

        instring = []
        outstring = []

        instring += [
'''
      CHARACTER(len=*), PARAMETER      :: serialized_string = &
         "qtb_rng_gaussian                         1 F T F   0.0000000000000000E+00&
                          12.0                12.0                12.0&
                          12.0                12.0                12.0&
                          12.0                12.0                12.0&
                          12.0                12.0                12.0&
                          12.0                12.0                12.0&
                          12.0                12.0                12.0"
'''
]

        outstring += [
'''
      CHARACTER(len=*), PARAMETER      :: serialized_string = &
         "qtb_rng_gaussian                         1 F T F   0.0000000000000000E+00&
&                          12.0                12.0                12.0&
&                          12.0                12.0                12.0&
&                          12.0                12.0                12.0&
&                          12.0                12.0                12.0&
&                          12.0                12.0                12.0&
&                          12.0                12.0                12.0"
'''
]

        instring += [
'''
      CHARACTER(len=*), PARAMETER      :: serialized_string = &
         "qtb_rng_gaussian                         1 F T F   0.0000000000000000E+00&
                 &         12.0                12.0                12.0&
                 &         12.0                12.0                12.0&
                 &         12.0                12.0                12.0&
                 &         12.0                12.0                12.0&
                 &         12.0                12.0                12.0&
                 &         12.0                12.0                12.0"
'''
]

        outstring += [
'''
      CHARACTER(len=*), PARAMETER      :: serialized_string = &
         "qtb_rng_gaussian                         1 F T F   0.0000000000000000E+00&
                 &         12.0                12.0                12.0&
                 &         12.0                12.0                12.0&
                 &         12.0                12.0                12.0&
                 &         12.0                12.0                12.0&
                 &         12.0                12.0                12.0&
                 &         12.0                12.0                12.0"
'''
]

        for instr, outstr in zip(instring, outstring):
            self.assert_fprettify_result([], instr, outstr)



def addtestmethod(testcase, fpath, ffile):
    """add a test method for each example."""

    def testmethod(testcase):
        """this is the test method invoked for each example."""

        dirpath_before = joinpath(BEFORE_DIR, fpath)
        dirpath_after = joinpath(AFTER_DIR, fpath)
        if not os.path.exists(dirpath_after):
            os.makedirs(dirpath_after)

        example_before = joinpath(dirpath_before, ffile)
        example_after = joinpath(dirpath_after, ffile)

        if os.path.isfile(example_after):
            os.remove(example_after)

        def test_result(path, info):
            return [os.path.relpath(path, BEFORE_DIR), info]

        with io.open(example_before, 'r', encoding='utf-8') as infile:

            outstring = io.StringIO()

            try:
                fprettify.reformat_ffile(infile, outstring)
                m = hashlib.sha256()
                m.update(outstring.getvalue().encode('utf-8'))

                test_info = "checksum"
                test_content = test_result(example_before, m.hexdigest())

                with io.open(example_after, 'w', encoding='utf-8') as outfile:
                    outfile.write(outstring.getvalue())
                FPrettifyTestCase.n_success += 1
            except FprettifyParseException as e:
                test_info = "parse error"
                fprettify.log_exception(e, test_info)
                test_content = test_result(example_before, test_info)
                FPrettifyTestCase.n_parsefail += 1
            except FprettifyInternalException as e:
                test_info = "internal error"
                fprettify.log_exception(e, test_info)
                test_content = test_result(example_before, test_info)
                FPrettifyTestCase.n_internalfail += 1
            except:  # pragma: no cover
                FPrettifyTestCase.n_unexpectedfail += 1
                raise

        after_exists = os.path.isfile(example_after)
        if after_exists:
            with io.open(example_before, 'r', encoding='utf-8') as infile:
                before_content = infile.read()
                before_nosp = re.sub(
                    r'\n{3,}', r'\n\n', before_content.lower().replace(' ', '').replace('\t', ''))

            with io.open(example_after, 'r', encoding='utf-8') as outfile:
                after_content = outfile.read()
                after_nosp = after_content.lower().replace(' ', '')

            testcase.assertMultiLineEqual(before_nosp, after_nosp)

        sep_str = ' : '
        with io.open(RESULT_FILE, 'r', encoding='utf-8') as infile:
            found = False
            for line in infile:
                line_content = line.strip().split(sep_str)
                if line_content[0] == test_content[0]:
                    found = True
                    eprint(test_info, end=" ")
                    msg = '{} (old) != {} (new)'.format(
                        line_content[1], test_content[1])
                    if test_info == "checksum" and after_exists and after_content.count('\n') < 10000:
                        # difflib can not handle large files
                        result = list(difflib.unified_diff(before_content.splitlines(
                            True), after_content.splitlines(True), fromfile=test_content[0], tofile=line_content[0]))
                        msg += '\n' + ''.join(result)
                    try:
                        testcase.assertEqual(
                            line_content[1], test_content[1], msg)
                    except AssertionError:  # pragma: no cover
                        FPrettifyTestCase.write_result(
                            FAILED_FILE, test_content, sep_str)
                        raise
                    break

        if not found:  # pragma: no cover
            eprint(test_info + " new", end=" ")
            FPrettifyTestCase.write_result(RESULT_FILE, test_content, sep_str)

    # not sure why this even works, using "test something" (with a space) as function name...
    # however it gives optimal test output
    testmethod.__name__ = ("test " + joinpath(fpath, ffile))

    setattr(testcase, testmethod.__name__, testmethod)

# make sure all directories exist
if not os.path.exists(BEFORE_DIR):  # pragma: no cover
    os.makedirs(BEFORE_DIR)
if not os.path.exists(AFTER_DIR):  # pragma: no cover
    os.makedirs(AFTER_DIR)
if not os.path.exists(RESULT_DIR):  # pragma: no cover
    os.makedirs(RESULT_DIR)
if not os.path.exists(RESULT_FILE):  # pragma: no cover
    io.open(RESULT_FILE, 'w', encoding='utf-8').close()
if os.path.exists(FAILED_FILE):  # pragma: no cover
    # erase failures from previous testers
    io.open(FAILED_FILE, 'w', encoding='utf-8').close()

# this prepares FPrettifyTestCase class when module is loaded by unittest
for dirpath, _, filenames in os.walk(BEFORE_DIR):
    for example in [f for f in filenames if any(f.endswith(_) for _ in fprettify.FORTRAN_EXTENSIONS)]:
        rel_dirpath = os.path.relpath(dirpath, start=BEFORE_DIR)
        addtestmethod(FPrettifyTestCase, rel_dirpath, example)
