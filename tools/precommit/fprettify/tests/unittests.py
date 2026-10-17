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

from __future__ import absolute_import, division, print_function, unicode_literals

import io
import logging
import os
import subprocess
import sys

sys.stderr = io.TextIOWrapper(
    sys.stderr.detach(), encoding="UTF-8", line_buffering=True
)

import fprettify
from fprettify.tests.test_common import RUNSCRIPT, FprettifyTestCase

fprettify.set_fprettify_logger(logging.ERROR)


class FprettifyUnitTestCase(FprettifyTestCase):
    def assert_fprettify_result(self, args, instring, outstring_exp):
        """
        assert that result of calling fprettify with args on instring gives
        outstring_exp
        """

        parser = fprettify.get_arg_parser()
        args = parser.parse_args(args)
        args = fprettify.process_args(args)

        outfile = io.StringIO()
        infile = io.StringIO(instring)

        fprettify.reformat_ffile(infile, outfile, orig_filename="StringIO", **args)
        outstring = outfile.getvalue()
        self.assertEqual(outstring_exp.rstrip(), outstring.rstrip())

    def test_whitespace(self):
        """simple test for whitespace formatting options -w in [0, 1, 2]"""
        instring = "(/-a-b-(a+b-c)/(-c)*d**e,f[1]%v/)"
        outstring_exp = [
            "(/-a-b-(a+b-c)/(-c)*d**e,f[1]%v/)",
            "(/-a-b-(a+b-c)/(-c)*d**e, f[1]%v/)",
            "(/-a - b - (a + b - c)/(-c)*d**e, f[1]%v/)",
            "(/-a - b - (a + b - c) / (-c) * d**e, f[1]%v/)",
        ]

        outstring = []
        for w, out in zip(range(0, 4), outstring_exp):
            args = ["-w", str(w)]
            self.assert_fprettify_result(args, instring, out)

    def test_type_selector(self):
        """test for whitespace formatting option -w 4"""
        instring = "A%component=func(mytype%a,mytype%abc+mytype%abcd)"
        outstring_exp = "A % component = func(mytype % a, mytype % abc + mytype % abcd)"

        self.assert_fprettify_result(["-w 4"], instring, outstring_exp)

    def test_concat(self):
        """test for concat operator whitespace formatting"""
        instring = "str=a//b//c"
        outstring_w0 = "str=a//b//c"
        outstring_w2 = "str = a//b//c"
        outstring_w4 = "str = a // b // c"
        outstring_explicit = "str = a // b // c"
        instring_in_string = 'msg = "URL: http://example.com"'
        instring_in_comment = "a = b  ! http://example.com"

        self.assert_fprettify_result(["-w", "0"], instring, outstring_w0)
        self.assert_fprettify_result(["-w", "2"], instring, outstring_w2)
        self.assert_fprettify_result(["-w", "4"], instring, outstring_w4)
        self.assert_fprettify_result(
            ["--whitespace-concat"], instring, outstring_explicit
        )
        self.assert_fprettify_result([], instring_in_string, instring_in_string)
        self.assert_fprettify_result([], instring_in_comment, instring_in_comment)

    def test_indent(self):
        """simple test for indent options -i in [0, 3, 4]"""

        indents = [0, 3, 4]

        instring = "iF(teSt)ThEn\nCaLl subr(a,b,&\nc,(/d,&\ne,f/))\nEnD iF"
        outstring_exp = [
            "iF (teSt) ThEn\n"
            + " " * ind
            + "CaLl subr(a, b, &\n"
            + " " * (10 + ind)
            + "c, (/d, &\n"
            + " " * (15 + ind)
            + "e, f/))\nEnD iF"
            for ind in indents
        ]

        for ind, out in zip(indents, outstring_exp):
            args = ["-i", str(ind)]
            self.assert_fprettify_result(args, instring, out)

    def test_nested(self):
        """test correct indentation of nested loops"""
        instring = (
            "integer :: i,j\ndo i=1,2\ndo j=1,3\n" "print*,i,j,i*j\nend do\nend do"
        )
        outstring_exp_default = (
            "integer :: i, j\ndo i = 1, 2\ndo j = 1, 3\n"
            "   print *, i, j, i*j\nend do\nend do"
        )
        outstring_exp_strict = (
            "integer :: i, j\ndo i = 1, 2\n   do j = 1, 3\n"
            "      print *, i, j, i*j\n   end do\nend do"
        )

        self.assert_fprettify_result([], instring, outstring_exp_default)
        self.assert_fprettify_result(
            ["--strict-indent"], instring, outstring_exp_strict
        )

    def test_reset_indent(self):
        """test of reset indentation at file start"""
        instring = (
            "integer :: i,j\ndo i=1,2\ndo j=1,3\n" "print*,i,j,i*j\nend do\nend do",
            "   module a\ninteger :: 1\n",
            "     module a\nend\nend",
        )
        outstring = (
            "integer :: i, j\ndo i = 1, 2\ndo j = 1, 3\n"
            "   print *, i, j, i*j\nend do\nend do",
            "module a\n   integer :: 1",
            "module a\nend\nend",
        )

        for ind, out in zip(instring, outstring):
            self.assert_fprettify_result([], ind, out)

    def test_disable(self):
        """test disabling indentation and/or whitespace formatting"""
        instring = (
            "if(&\nl==111)&\n then\n   do m   =1,  2\n A=&\nB+C\n    end  do;   endif"
        )
        outstring_exp_default = (
            "if ( &\n   l == 111) &\n   then\n   do m = 1, 2\n"
            "      A = &\n         B + C\n   end do; end if"
        )
        outstring_exp_nowhitespace = (
            "if(&\n   l==111)&\n   then\n   do m   =1,  2\n"
            "      A=&\n         B+C\n   end  do; endif"
        )
        outstring_exp_noindent = (
            "if ( &\nl == 111) &\n then\n   do m = 1, 2\n"
            " A = &\nB + C\n    end do;   end if"
        )

        self.assert_fprettify_result([], instring, outstring_exp_default)
        self.assert_fprettify_result(
            ["--disable-whitespace"], instring, outstring_exp_nowhitespace
        )
        self.assert_fprettify_result(
            ["--disable-indent"], instring, outstring_exp_noindent
        )
        self.assert_fprettify_result(
            ["--disable-indent", "--disable-whitespace"], instring, instring
        )

    def test_comments(self):
        """test options related to comments"""
        instring = (
            "TYPE mytype\n!  c1\n  !c2\n   INTEGER :: a   !  c3\n"
            "   REAL :: b, &   ! c4\n! c5\n                  ! c6\n"
            "           d      ! c7\nEND TYPE  ! c8"
        )
        outstring_exp_default = (
            "TYPE mytype\n!  c1\n   !c2\n   INTEGER :: a   !  c3\n"
            "   REAL :: b, &   ! c4\n           ! c5\n           ! c6\n"
            "           d      ! c7\nEND TYPE  ! c8"
        )
        outstring_exp_strip = (
            "TYPE mytype\n!  c1\n   !c2\n   INTEGER :: a !  c3\n"
            "   REAL :: b, & ! c4\n           ! c5\n           ! c6\n"
            "           d ! c7\nEND TYPE ! c8"
        )
        outstring_exp_strip_spacing3 = (
            "TYPE mytype\n!  c1\n   !c2\n   INTEGER :: a   !  c3\n"
            "   REAL :: b, &   ! c4\n           ! c5\n           ! c6\n"
            "           d   ! c7\nEND TYPE   ! c8"
        )

        self.assert_fprettify_result([], instring, outstring_exp_default)
        self.assert_fprettify_result(
            ["--strip-comments"], instring, outstring_exp_strip
        )
        self.assert_fprettify_result(
            ["--strip-comments", "--comment-spacing", "3"],
            instring,
            outstring_exp_strip_spacing3,
        )

    def test_directive(self):
        """
        test deactivate directives '!&' (inline) and '!&<', '!&>' (block)
        and manual alignment (continuation line prefixed with '&')
        """

        # manual alignment
        instring = "align_me = [ -1,  10,0,  &\n    &     0,1000 ,  0,&\n            &0 , -1,  1]"
        outstring_exp = (
            "align_me = [-1, 10, 0,  &\n    &     0, 1000, 0,&\n            &0, -1, 1]"
        )
        self.assert_fprettify_result([], instring, outstring_exp)

        # inline deactivate
        instring2 = "\n".join(_ + " !&" for _ in instring.splitlines())
        outstring_exp = instring2
        self.assert_fprettify_result([], instring2, outstring_exp)

        # block deactivate
        instring3 = "!&<\n" + instring + "\n!&>"
        outstring_exp = instring3
        self.assert_fprettify_result([], instring3, outstring_exp)

    def test_io(self):
        """simple test for io (file inplace, stdin & stdout)"""

        # io and unicode
        outstring = []
        instring = "CALL  alien_invasion( 👽 )"
        outstring_exp = "CALL alien_invasion(👽)"

        alien_file = "alien_invasion.f90"
        if os.path.isfile(alien_file):
            raise Exception("remove file alien_invasion.f90")  # pragma: no cover

        try:
            with io.open(alien_file, "w", encoding="utf-8") as infile:
                infile.write(instring)

            # testing stdin --> stdout
            p1 = subprocess.Popen(
                RUNSCRIPT, stdout=subprocess.PIPE, stdin=subprocess.PIPE
            )
            outstring.append(
                p1.communicate(instring.encode("UTF-8"))[0].decode("UTF-8")
            )

            # testing file --> stdout
            p1 = subprocess.Popen(
                [RUNSCRIPT, alien_file, "--stdout"], stdout=subprocess.PIPE
            )
            outstring.append(
                p1.communicate(instring.encode("UTF-8")[0])[0].decode("UTF-8")
            )

            # testing file --> file (inplace)
            p1 = subprocess.Popen([RUNSCRIPT, alien_file])
            p1.wait()

            with io.open(alien_file, "r", encoding="utf-8") as infile:
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
        instring = "use A,only:B=>C,&\nD=>E"
        outstring = "use A, only: B => C, &\n             D => E"
        self.assert_fprettify_result([], instring, outstring)

    def test_use(self):
        """test for alignment of use statements"""
        instring1 = "use A,only:B,C,&\nD,E"
        instring2 = "use A,only:&\nB,C,D,E"
        outstring1 = "use A, only: B, C, &\n             D, E"
        outstring2 = "use A, only: &\n   B, C, D, E"
        self.assert_fprettify_result([], instring1, outstring1)
        self.assert_fprettify_result([], instring2, outstring2)

    def test_wrongkind(self):
        """test whitespacing of deprecated kind definition"""
        instring = [
            "REAL*8 :: r, f  !  some reals",
            "REAL * 8 :: r, f  !  some reals",
            "INTEGER * 4 :: c, i  !  some integers",
            "INTEGER*4 :: c, i  !  some integers",
        ]
        outstring = [
            "REAL*8 :: r, f  !  some reals",
            "REAL*8 :: r, f  !  some reals",
            "INTEGER*4 :: c, i  !  some integers",
            "INTEGER*4 :: c, i  !  some integers",
        ]

        for i in range(0, len(instring)):
            self.assert_fprettify_result([], instring[i], outstring[i])

    def test_new_intrinsics(self):
        """test new I/O intrinsics"""
        instring = ["REWIND(12)", "BACKSPACE(13)", "INQUIRE(14)"]
        outstring = ["REWIND (12)", "BACKSPACE (13)", "INQUIRE (14)"]

        for i in range(0, len(instring)):
            self.assert_fprettify_result([], instring[i], outstring[i])

    def test_associate(self):
        """test correct formatting of associate construct"""
        instring = "associate(a=>b , c  =>d ,e=> f  )\n" "e=a+c\n" "end associate"
        outstring = (
            "associate (a => b, c => d, e => f)\n" "   e = a + c\n" "end associate"
        )

        self.assert_fprettify_result([], instring, outstring)

    def test_line_length(self):
        """test line length option"""
        instring = [
            "REAL(KIND=4) :: r,f  !  some reals",
            "if(   min == max.and.min .eq. thres  )",
            "INQUIRE(14)",
        ]
        instring_ = "if( min == max.and.min .eq. thres ) one_really_long_function_call_to_hit_the_line_limit(parameter1, parameter2,parameter3,parameter4,parameter5,err) ! this line would be too long"
        outstring = [
            "REAL(KIND=4) :: r, f  !  some reals",
            "REAL(KIND=4) :: r, f\n!  some reals",
            "if (min == max .and. min .eq. thres)",
            "if(   min == max.and.min .eq. thres  )",
            "INQUIRE (14)",
            "INQUIRE (14)",
        ]
        outstring_ = [
            "\n".join(
                [
                    "if (min == max .and. min .eq. thres) one_really_long_function_call_to_hit_the_line_limit(parameter1, parameter2, parameter3, &",
                    "                                                                                         parameter4, parameter5, err)",
                    "! this line would be too long",
                ]
            ),
            "if (min == max .and. min .eq. thres) one_really_long_function_call_to_hit_the_line_limit(parameter1, parameter2, parameter3, parameter4, parameter5, err) ! this line would be too long",
        ]

        # test shorter lines first, after all the actual length doesn't matter
        for i in range(0, len(instring)):
            self.assert_fprettify_result(["-S"], instring[i], outstring[2 * i])
            self.assert_fprettify_result(
                ["-S", "-l 20"], instring[i], outstring[2 * i + 1]
            )
        # now test a long line
        self.assert_fprettify_result(["-S"], instring_, outstring_[0])
        self.assert_fprettify_result(["-S", "-l 0"], instring_, outstring_[1])

    def test_relation_replacement(self):
        """test replacement of relational statements"""
        instring = [
            "if ( min < max .and. min .lt. thres)",
            "if (min > max .and. min .gt. thres )",
            "if (   min == max .and. min .eq. thres  )",
            "if(min /= max .and. min .ne. thres)",
            "if(min >= max .and. min .ge. thres )",
            "if( min <= max .and. min .le. thres)",
            "'==== heading",
            "if (vtk%my_rank .eq. 0) write (vtk%filehandle_par, '(\"<DataArray",
            '\'("</Collection>",',
            "if (abc(1) .lt. -bca .or. &\n qwe .gt. ewq) then",
        ]
        f_outstring = [
            "if (min .lt. max .and. min .lt. thres)",
            "if (min .gt. max .and. min .gt. thres)",
            "if (min .eq. max .and. min .eq. thres)",
            "if (min .ne. max .and. min .ne. thres)",
            "if (min .ge. max .and. min .ge. thres)",
            "if (min .le. max .and. min .le. thres)",
            "'==== heading",
            "if (vtk%my_rank .eq. 0) write (vtk%filehandle_par, '(\"<DataArray",
            '\'("</Collection>",',
            "if (abc(1) .lt. -bca .or. &\n    qwe .gt. ewq) then",
        ]
        c_outstring = [
            "if (min < max .and. min < thres)",
            "if (min > max .and. min > thres)",
            "if (min == max .and. min == thres)",
            "if (min /= max .and. min /= thres)",
            "if (min >= max .and. min >= thres)",
            "if (min <= max .and. min <= thres)",
            "'==== heading",
            "if (vtk%my_rank == 0) write (vtk%filehandle_par, '(\"<DataArray",
            '\'("</Collection>",',
            "if (abc(1) < -bca .or. &\n    qwe > ewq) then",
        ]
        for i in range(0, len(instring)):
            self.assert_fprettify_result(
                ["--enable-replacements", "--c-relations"], instring[i], c_outstring[i]
            )
            self.assert_fprettify_result(
                ["--enable-replacements"], instring[i], f_outstring[i]
            )

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
            "PURE SUBROUTINE mypure()",
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
            "pure subroutine mypure()",
        )
        for i in range(len(instring)):
            self.assert_fprettify_result(
                ["--case", "1", "1", "1", "2"], instring[i], outstring[i]
            )

    def test_do(self):
        """test correct parsing of do statement"""
        instring = "do = 1\nb = 2"

        self.assert_fprettify_result([], instring, instring)

    def test_omp(self):
        """test formatting of omp directives"""
        instring = (
            "PROGRAM test_omp\n"
            " !$OMP    PARALLEL DO\n"
            "b=4\n"
            "!$a=b\n"
            "!$  a=b\n"
            "   !$    c=b\n"
            "!$acc parallel loop\n"
            "!$OMP END  PARALLEL DO\n"
            "END PROGRAM"
        )
        outstring = (
            "PROGRAM test_omp\n"
            "!$OMP    PARALLEL DO\n"
            "   b = 4\n"
            "!$a=b\n"
            "!$ a = b\n"
            "!$ c = b\n"
            "!$acc parallel loop\n"
            "!$OMP END  PARALLEL DO\n"
            "END PROGRAM"
        )

        self.assert_fprettify_result([], instring, outstring)

    def test_ford(self):
        """test formatting of ford comments"""
        instring = (
            "   a = b\n"
            "     !!  ford docu\n"
            "b=c\n"
            "  !! ford docu\n"
            "subroutine test(a,b,&\n"
            "  !! ford docu\n"
            "  c, d, e)"
        )
        outstring = (
            "   a = b\n"
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

        instring += ["""
#:if DEBUG>  0
print *, "hola"
if( .not. (${cond}$) ) then
#:if  ASSERT(cond)
print *, "Assert failed!"
#:endif
error stop
end if
#:endif
"""]

        outstring += ["""
#:if DEBUG>  0
   print *, "hola"
   if (.not. (${cond}$)) then
      #:if  ASSERT(cond)
         print *, "Assert failed!"
      #:endif
      error stop
   end if
#:endif
"""]

        instring += ["""
if  (.not. (${cond}$)) then
   #:for element in list
   print *, "Element is in list!"
 #:endfor
   error stop
end if
"""]

        outstring += ["""
if (.not. (${cond}$)) then
   #:for element in list
      print *, "Element is in list!"
   #:endfor
   error stop
end if
"""]

        instring += ["""
#:if aa > 1
print  *, "Number is more than 1"
if (condition) then
     #:def something
   print *, "Added Definition!"
   #:enddef
end if
#:endif
"""]

        outstring += ["""
#:if aa > 1
   print *, "Number is more than 1"
   if (condition) then
      #:def something
         print *, "Added Definition!"
      #:enddef
   end if
#:endif
"""]

        instring += ["""
#:def DEBUG_CODE( code)
  #:if DEBUG > 0
    $:code
  #:endif
#:enddef DEBUG_CODE
"""]

        outstring += ["""
#:def DEBUG_CODE( code)
   #:if DEBUG > 0
      $:code
   #:endif
#:enddef DEBUG_CODE
"""]

        instring += ["""
#:block DEBUG_CODE
  if (a <b) then
    print *, "DEBUG: a is less than b"
  end if
#:endblock  DEBUG_CODE
"""]

        outstring += ["""
#:block DEBUG_CODE
   if (a < b) then
      print *, "DEBUG: a is less than b"
   end if
#:endblock  DEBUG_CODE
"""]

        instring += ["""
#:call DEBUG_CODE
  if (a < b) then
    print *, "DEBUG: a is less than b"
  end if
#:endcall DEBUG_CODE
"""]

        outstring += ["""
#:call DEBUG_CODE
   if (a < b) then
      print *, "DEBUG: a is less than b"
   end if
#:endcall DEBUG_CODE
"""]

        instring += ["""
#:if DEBUG > 0
print *, "hola"
if (.not. (${cond}$)) then
   #:mute
   print *, "Muted"
   #:endmute
   error stop
end if
#:endif
"""]

        outstring += ["""
#:if DEBUG > 0
   print *, "hola"
   if (.not. (${cond}$)) then
      #:mute
         print *, "Muted"
      #:endmute
      error stop
   end if
#:endif
"""]

        instring += ["""
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
"""]

        outstring += ["""
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
"""]

        instring += ["""
#:if worktype
      ${worktype}$, &
#:else
      ${type}$, &
#:endif
         DIMENSION(${arr_exp}$), &
         POINTER :: work
"""]

        outstring += ["""
#:if worktype
${worktype}$, &
#:else
   ${type}$, &
#:endif
   DIMENSION(${arr_exp}$), &
   POINTER :: work
"""]

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

        self.assert_fprettify_result(
            ["--disable-indent-mod"], instring_mod, outstring_mod_disable
        )
        self.assert_fprettify_result(
            ["--disable-indent-mod"], instring_prog, outstring_prog_disable
        )

    def test_decl(self):
        """test formatting of declarations"""
        instring_1 = "integer    ::     a"
        instring_2 = "integer, dimension(:)    ::     a"
        outstring_1 = "integer :: a"
        outstring_2 = "integer, dimension(:) :: a"
        outstring_2_min = "integer, dimension(:)::a"

        self.assert_fprettify_result([], instring_1, instring_1)
        self.assert_fprettify_result([], instring_2, instring_2)
        self.assert_fprettify_result(["--enable-decl"], instring_1, outstring_1)
        self.assert_fprettify_result(["--enable-decl"], instring_2, outstring_2)
        self.assert_fprettify_result(
            ["--enable-decl", "--whitespace-decl=0"], instring_2, outstring_2_min
        )

    def test_statement_label(self):
        instring = "1003  FORMAT(2(1x, i4), 5x, '-', 5x, '-', 3x, '-', 5x, '-', 5x, '-', 8x, '-', 3x, &\n    1p, 2(1x, d10.3))"
        outstring = "1003  FORMAT(2(1x, i4), 5x, '-', 5x, '-', 3x, '-', 5x, '-', 5x, '-', 8x, '-', 3x, &\n             1p, 2(1x, d10.3))"
        self.assert_fprettify_result([], instring, outstring)

        instring = "print *, 'hello'\n1003  FORMAT(2(1x, i4), 5x, '-', 5x, '-', 3x, '-', 5x, '-', 5x, '-', 8x, '-', 3x, &\n    1p, 2(1x, d10.3))"
        outstring = "print *, 'hello'\n1003 FORMAT(2(1x, i4), 5x, '-', 5x, '-', 3x, '-', 5x, '-', 5x, '-', 8x, '-', 3x, &\n            1p, 2(1x, d10.3))"
        self.assert_fprettify_result([], instring, outstring)

    def test_statement_label_auto_break(self):
        instring = "1003  FORMAT(2(1x, i4), 5x, '-', 5x, '-', 3x, '-', 5x, '-', 5x, '-', 8x, '-', 3x, 1p, 2(1x, d10.3)) ! comment"
        # the comment does not fit within the line length limit at any
        # indentation, so it ends up right-aligned at the limit (idempotent)
        outstring = "1003  FORMAT(2(1x, i4), 5x, '-', 5x, '-', 3x, '-', 5x, '-', 5x, '-', 8x, '-', 3x, 1p, 2(1x, d10.3))\n ! comment"
        self.assert_fprettify_result(["--line-length=10"], instring, outstring)

    def test_multiline_str(self):

        instring = []
        outstring = []

        instring += ["""
      CHARACTER(len=*), PARAMETER      :: serialized_string = &
         "qtb_rng_gaussian                         1 F T F   0.0000000000000000E+00&
                          12.0                12.0                12.0&
                          12.0                12.0                12.0&
                          12.0                12.0                12.0&
                          12.0                12.0                12.0&
                          12.0                12.0                12.0&
                          12.0                12.0                12.0"
"""]

        outstring += ["""
      CHARACTER(len=*), PARAMETER      :: serialized_string = &
         "qtb_rng_gaussian                         1 F T F   0.0000000000000000E+00&
&                          12.0                12.0                12.0&
&                          12.0                12.0                12.0&
&                          12.0                12.0                12.0&
&                          12.0                12.0                12.0&
&                          12.0                12.0                12.0&
&                          12.0                12.0                12.0"
"""]

        instring += ["""
      CHARACTER(len=*), PARAMETER      :: serialized_string = &
         "qtb_rng_gaussian                         1 F T F   0.0000000000000000E+00&
                 &         12.0                12.0                12.0&
                 &         12.0                12.0                12.0&
                 &         12.0                12.0                12.0&
                 &         12.0                12.0                12.0&
                 &         12.0                12.0                12.0&
                 &         12.0                12.0                12.0"
"""]

        outstring += ["""
      CHARACTER(len=*), PARAMETER      :: serialized_string = &
         "qtb_rng_gaussian                         1 F T F   0.0000000000000000E+00&
                 &         12.0                12.0                12.0&
                 &         12.0                12.0                12.0&
                 &         12.0                12.0                12.0&
                 &         12.0                12.0                12.0&
                 &         12.0                12.0                12.0&
                 &         12.0                12.0                12.0"
"""]

        for instr, outstr in zip(instring, outstring):
            self.assert_fprettify_result([], instr, outstr)

    def test_label(self):
        instring = """
MODULE cp_lbfgs
CONTAINS
20000    FORMAT('RUNNING THE L-BFGS-B CODE', /, /,                          &
     &    'it    = iteration number', /,                                    &
     &    'Machine precision =', 1p, d10.3)
2        FORMAT('RUNNING THE L-BFGS-B CODE', /, /,                          &
     &    'it    = iteration number', /,                                    &
     &    'Machine precision =', 1p, d10.3)
20000    FORMAT('RUNNING THE L-BFGS-B CODE', /, /,                          &
          'it    = iteration number', /,                                    &
          'Machine precision =', 1p, d10.3)
2        FORMAT('RUNNING THE L-BFGS-B CODE', /, /,                          &
          'it    = iteration number', /,                                    &
          'Machine precision =', 1p, d10.3)
END MODULE
"""

        outstring = """
MODULE cp_lbfgs
CONTAINS
20000 FORMAT('RUNNING THE L-BFGS-B CODE', /, /,                          &
  &    'it    = iteration number', /,                                    &
  &    'Machine precision =', 1p, d10.3)
2  FORMAT('RUNNING THE L-BFGS-B CODE', /, /,                          &
&    'it    = iteration number', /,                                    &
&    'Machine precision =', 1p, d10.3)
20000 FORMAT('RUNNING THE L-BFGS-B CODE', /, /, &
             'it    = iteration number', /, &
             'Machine precision =', 1p, d10.3)
2  FORMAT('RUNNING THE L-BFGS-B CODE', /, /, &
          'it    = iteration number', /, &
          'Machine precision =', 1p, d10.3)
END MODULE
"""

        self.assert_fprettify_result([], instring, outstring)

    def test_ampersand_string(self):
        """test linebreaks within strings"""
        instring = [
            'write  ( * ,  * )  "a&\nstring"',
            'write  ( * ,  * )"a& \n    string"',
            'write  ( * ,  * )    "a& \n    &  string"',
        ]
        outstring = [
            'write (*, *) "a&\n&string"',
            'write (*, *) "a&\n&    string"',
            'write (*, *) "a&\n    &  string"',
        ]

        for instr, outstr in zip(instring, outstring):
            self.assert_fprettify_result([], instr, outstr)

    def test_first_line_non_code(self):
        """test whether first non-code line gets correctly indented"""
        instr = "  ! a comment\n     module mod\n ! a comment\nend"
        outstr = "! a comment\nmodule mod\n   ! a comment\nend"
        self.assert_fprettify_result([], instr, outstr)

        instr = "  ! a comment\n     program mod\n ! a comment\nend"
        outstr = "! a comment\nprogram mod\n   ! a comment\nend"
        self.assert_fprettify_result([], instr, outstr)

        instr = "  ! a comment\n     function fun()\n ! a comment\nend"
        outstr = "     ! a comment\n     function fun()\n        ! a comment\n     end"
        self.assert_fprettify_result([], instr, outstr)

    def test_indent_preserves_line_length_limit(self):
        """indentation should remain stable when exceeding line length"""
        in_lines = [
            "subroutine demo(tokens, stmt_start)",
            "   type(dummy), intent(in) :: tokens(:)",
            "   integer, intent(in) :: stmt_start",
            "   integer :: i, nesting_level",
            "",
            '   if (tokens(stmt_start)%text == "if") then',
            '      if (tokens(i)%text == "endif") then',
            "         nesting_level = nesting_level - 1",
            '      else if (tokens(i)%text == "end" .and. i + 1 <= size(tokens) .and. &',
            '               tokens(i + 1)%kind == TK_KEYWORD .and. tokens(i + 1)%text == "if") then',
            "         nesting_level = nesting_level - 1",
            "      end if",
            "   end if",
            "",
            '   if (tokens(i)%text == "end") then',
            "      if (i + 1 <= size(tokens) .and. tokens(i + 1)%kind == TK_KEYWORD) then",
            '         if (tokens(i + 1)%text == "do" .and. tokens(stmt_start)%text == "do") then',
            "            nesting_level = nesting_level - 1",
            '         else if (tokens(i + 1)%text == "select" .and. tokens(stmt_start)%text == "select") then',
            "            nesting_level = nesting_level - 1",
            '         else if (tokens(i + 1)%text == "where" .and. tokens(stmt_start)%text == "where") then',
            "            nesting_level = nesting_level - 1",
            "         end if",
            "      end if",
            "   end if",
            "end subroutine demo",
            "",
        ]

        out_lines = [
            "subroutine demo(tokens, stmt_start)",
            "   type(dummy), intent(in) :: tokens(:)",
            "   integer, intent(in) :: stmt_start",
            "   integer :: i, nesting_level",
            "",
            '   if (tokens(stmt_start)%text == "if") then',
            '      if (tokens(i)%text == "endif") then',
            "         nesting_level = nesting_level - 1",
            '      else if (tokens(i)%text == "end" .and. i + 1 <= size(tokens) .and. &',
            '               tokens(i + 1)%kind == TK_KEYWORD .and. tokens(i + 1)%text == "if") then',
            "         nesting_level = nesting_level - 1",
            "      end if",
            "   end if",
            "",
            '   if (tokens(i)%text == "end") then',
            "      if (i + 1 <= size(tokens) .and. tokens(i + 1)%kind == TK_KEYWORD) then",
            '         if (tokens(i + 1)%text == "do" .and. tokens(stmt_start)%text == "do") then',
            "            nesting_level = nesting_level - 1",
            '         else if (tokens(i + 1)%text == "select" .and. tokens(stmt_start)%text == &',
            '                  "select") then',
            "            nesting_level = nesting_level - 1",
            '         else if (tokens(i + 1)%text == "where" .and. tokens(stmt_start)%text == &',
            '                  "where") then',
            "            nesting_level = nesting_level - 1",
            "         end if",
            "      end if",
            "   end if",
            "end subroutine demo",
            "",
        ]

        instring = "\n".join(in_lines)
        outstring_exp = "\n".join(out_lines)

        self.assert_fprettify_result(["-l", "90"], instring, outstring_exp)

    def test_auto_split_long_logical_line(self):
        """automatically split long logical lines that exceed the limit after indentation"""
        instring = (
            "subroutine demo()\n"
            "    integer :: a\n"
            "    if (this_condition_is_lengthy .or. second_lengthy_condition) cycle\n"
            "end subroutine demo"
        )

        outstring_exp = (
            "subroutine demo()\n"
            "    integer :: a\n"
            "    if (this_condition_is_lengthy .or. &\n"
            "        second_lengthy_condition) cycle\n"
            "end subroutine demo"
        )

        self.assert_fprettify_result(["-i", "4", "-l", "68"], instring, outstring_exp)

    def test_auto_split_handles_bang_in_string(self):
        """ensure split logic ignores exclamation marks inside string literals"""
        instring = (
            "subroutine demo(str)\n"
            "    character(len=*), intent(in) :: str\n"
            '    if (str .eq. "This string has a ! bang inside") print *, str//", wow!"\n'
            "end subroutine demo"
        )

        outstring_exp = (
            "subroutine demo(str)\n"
            "    character(len=*), intent(in) :: str\n"
            '    if (str .eq. "This string has a ! bang inside") print *, &\n'
            '        str//", wow!"\n'
            "end subroutine demo"
        )

        self.assert_fprettify_result(["-i", "4", "-l", "72"], instring, outstring_exp)

    def test_auto_split_after_indent_adjustment(self):
        """splitting must also run during the indentation pass to stay idempotent"""
        instring = (
            "program demo\n"
            "    integer :: i\n"
            "    if (.true.) then\n"
            "        if (.true.) then\n"
            '  if (i > 1 .and. this_is_a_pretty_freaking_long_parameter_name .eq. 42) print *, "too long"\n'
            "        end if\n"
            "    end if\n"
            "end program demo\n"
        )

        outstring_exp = (
            "program demo\n"
            "    integer :: i\n"
            "    if (.true.) then\n"
            "        if (.true.) then\n"
            "            if (i > 1 .and. this_is_a_pretty_freaking_long_parameter_name .eq. 42) print &\n"
            '                *, "too long"\n'
            "        end if\n"
            "    end if\n"
            "end program demo\n"
        )

        self.assert_fprettify_result(["-i", "4", "-l", "100"], instring, outstring_exp)

    def test_auto_split_when_whitespace_disabled(self):
        """indent-only runs must still split long logical lines"""
        instring = (
            "program demo\n"
            "    if (.true.) then\n"
            "      if (.true.) then\n"
            '  if (i > 1 .and. identifier_that_is_far_too_long .eq. 42) print *, "oops"\n'
            "      end if\n"
            "    end if\n"
            "end program demo\n"
        )

        outstring_exp = (
            "program demo\n"
            "    if (.true.) then\n"
            "        if (.true.) then\n"
            "            if (i > 1 .and. identifier_that_is_far_too_long .eq. 42) &\n"
            '                print *, "oops"\n'
            "        end if\n"
            "    end if\n"
            "end program demo\n"
        )

        self.assert_fprettify_result(
            ["-i", "4", "-l", "70", "--disable-whitespace"], instring, outstring_exp
        )

    def test_line_length_detaches_inline_comment(self):
        """inline comments should move to their own line when they exceed the limit"""
        instring = (
            "program demo\n"
            "    if (.true.) then\n"
            "        print *, 'prefix '//'and '//'suffix' ! trailing comment\n"
            "    end if\n"
            "end program demo\n"
        )

        outstring_exp = (
            "program demo\n"
            "    if (.true.) then\n"
            "        print *, 'prefix '//'and '//'suffix'\n"
            "        ! trailing comment\n"
            "    end if\n"
            "end program demo\n"
        )

        self.assert_fprettify_result(["-i", "4", "-l", "60"], instring, outstring_exp)

    def test_line_length_comment_then_split(self):
        """detaching the comment must still allow the code line to split further"""
        instring = (
            "program demo\n"
            "    if (.true.) then\n"
            "        if (.true.) then\n"
            "            if (foo_bar_identifier .and. bar_baz_identifier) print *, long_identifier, another_long_identifier ! note\n"
            "        end if\n"
            "    end if\n"
            "end program demo\n"
        )

        outstring_exp = (
            "program demo\n"
            "    if (.true.) then\n"
            "        if (.true.) then\n"
            "            if (foo_bar_identifier .and. bar_baz_identifier) print *, &\n"
            "                long_identifier, another_long_identifier\n"
            "            ! note\n"
            "        end if\n"
            "    end if\n"
            "end program demo\n"
        )

        self.assert_fprettify_result(["-i", "4", "-l", "72"], instring, outstring_exp)
