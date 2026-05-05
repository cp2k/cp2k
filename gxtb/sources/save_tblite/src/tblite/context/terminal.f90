! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later
!
! tblite is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! tblite is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with tblite.  If not, see <https://www.gnu.org/licenses/>.

!> @file tblite/context/terminal.f90
!> Provides a terminal class to safely use ANSI escape sequences

!> Support for ANSI escape sequences to get colorful terminal output
module tblite_context_terminal
   use mctc_env, only : i1
   implicit none
   private

   public :: escape, operator(+), operator(//)

   !> Container for terminal escape code
   type :: color
      !> Style descriptor
      integer(i1) :: style = -1_i1
      !> Background color descriptor
      integer(i1) :: bg = -1_i1
      !> Foreground color descriptor
      integer(i1) :: fg = -1_i1
   end type

   interface operator(+)
      module procedure :: add
   end interface operator(+)

   interface operator(//)
      module procedure :: concat_left
      module procedure :: concat_right
   end interface operator(//)

   type(color), parameter :: &
      reset = color(style=0_i1), &
      bold = color(style=1_i1), &
      dim = color(style=2_i1), &
      italic = color(style=3_i1), &
      underline = color(style=4_i1), &
      blink = color(style=5_i1), &
      reverse = color(style=7_i1), &
      hidden = color(style=8_i1)

   type(color), parameter :: &
      black = color(fg=0_i1), &
      red = color(fg=1_i1), &
      green = color(fg=2_i1), &
      yellow = color(fg=3_i1), &
      blue = color(fg=4_i1), &
      magenta = color(fg=5_i1), &
      cyan = color(fg=6_i1), &
      white = color(fg=7_i1)

   type(color), parameter :: &
      bg_black = color(bg=0_i1), &
      bg_red = color(bg=1_i1), &
      bg_green = color(bg=2_i1), &
      bg_yellow = color(bg=3_i1), &
      bg_blue = color(bg=4_i1), &
      bg_magenta = color(bg=5_i1), &
      bg_cyan = color(bg=6_i1), &
      bg_white = color(bg=7_i1)


   !> Colorizer class for handling colorful output in the terminal
   type, public :: context_terminal

      type(color) :: &
         reset = color(), &
         bold = color(), &
         dim = color(), &
         italic = color(), &
         underline = color(), &
         blink = color(), &
         reverse = color(), &
         hidden = color()

      type(color) :: &
         black = color(), &
         red = color(), &
         green = color(), &
         yellow = color(), &
         blue = color(), &
         magenta = color(), &
         cyan = color(), &
         white = color()

      type(color) :: &
         bg_black = color(), &
         bg_red = color(), &
         bg_green = color(), &
         bg_yellow = color(), &
         bg_blue = color(), &
         bg_magenta = color(), &
         bg_cyan = color(), &
         bg_white = color()

      type(color) :: &
         bold_black = color(), &
         bold_red = color(), &
         bold_green = color(), &
         bold_yellow = color(), &
         bold_blue = color(), &
         bold_magenta = color(), &
         bold_cyan = color(), &
         bold_white = color()
   end type context_terminal

   !> Constructor for the colorizer
   interface context_terminal
      module procedure :: new_terminal
   end interface context_terminal

contains


!> Create a new colorizer object
function new_terminal(use_color) result(new)
   !> Enable color output
   logical, intent(in) :: use_color
   !> New instance of the colorizer
   type(context_terminal) :: new

   if (use_color) then
      new%reset = reset
      new%bold = bold
      new%dim = dim
      new%italic = italic
      new%underline = underline
      new%blink = blink
      new%reverse = reverse
      new%hidden = hidden
      new%black = black
      new%red = red
      new%green = green
      new%yellow = yellow
      new%blue = blue
      new%magenta = magenta
      new%cyan = cyan
      new%white = white
      new%bg_black = bg_black
      new%bg_red = bg_red
      new%bg_green = bg_green
      new%bg_yellow = bg_yellow
      new%bg_blue = bg_blue
      new%bg_magenta = bg_magenta
      new%bg_cyan = bg_cyan
      new%bg_white = bg_white
      new%bold_black = bold + black
      new%bold_red = bold + red
      new%bold_green = bold + green
      new%bold_yellow = bold + yellow
      new%bold_blue = bold + blue
      new%bold_magenta = bold + magenta
      new%bold_cyan = bold + cyan
      new%bold_white = bold + white
   end if
end function new_terminal

!> Add two escape sequences, attributes in the right value override the left value ones.
pure function add(lval, rval) result(code)
   !> First escape code
   type(color), intent(in) :: lval
   !> Second escape code
   type(color), intent(in) :: rval
   !> Combined escape code
   type(color) :: code

   code = color( &
      style=merge(rval%style, lval%style, rval%style >= 0), &
      fg=merge(rval%fg, lval%fg, rval%fg >= 0), &
      bg=merge(rval%bg, lval%bg, rval%bg >= 0))
end function add

!> Concatenate an escape code with a string and turn it into an actual escape sequence
pure function concat_left(lval, code) result(str)
   !> String to add the escape code to
   character(len=*), intent(in) :: lval
   !> Escape sequence
   type(color), intent(in) :: code
   !> Concatenated string
   character(len=:), allocatable :: str

   str = lval // escape(code)
end function concat_left

!> Concatenate an escape code with a string and turn it into an actual escape sequence
pure function concat_right(code, rval) result(str)
   !> String to add the escape code to
   character(len=*), intent(in) :: rval
   !> Escape sequence
   type(color), intent(in) :: code
   !> Concatenated string
   character(len=:), allocatable :: str

   str = escape(code) // rval
end function concat_right

!> Transform a color code into an actual ANSI escape sequence
pure function escape(code) result(str)
   !> Color code to be used
   type(color), intent(in) :: code
   !> ANSI escape sequence representing the color code
   character(len=:), allocatable :: str
   character, parameter :: chars(0:9) = &
      ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]

   if (anycolor(code)) then
      str = achar(27) // "[0"  ! Always reset the style
      if (code%style > 0 .and. code%style < 10) str = str // ";" // chars(code%style)
      if (code%fg >= 0 .and. code%fg < 10) str = str // ";3" // chars(code%fg)
      if (code%bg >= 0 .and. code%bg < 10) str = str // ";4" // chars(code%bg)
      str = str // "m"
   else
      str = ""
   end if
end function escape

!> Check whether the code describes any color or is just a stub
pure function anycolor(code)
   !> Escape sequence
   type(color), intent(in) :: code
   !> Any color / style is active
   logical :: anycolor

   anycolor = code%fg >= 0 .or. code%bg >= 0 .or. code%style >= 0
end function anycolor

end module tblite_context_terminal
