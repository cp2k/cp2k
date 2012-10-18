!! Copyright (C) 2003-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU Lesser General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!! GNU Lesser General Public License for more details.
!!
!! You should have received a copy of the GNU Lesser General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id: libxc.f90 3550 2007-11-19 14:32:49Z marques $
!-------------------------------------------------------------------
MODULE xc_f90_types_m



  INTEGER, PUBLIC, PARAMETER :: xc_f90_kind = SELECTED_REAL_KIND(14)


  TYPE xc_f90_pointer_t
    PRIVATE
    INTEGER, POINTER :: buffer
  END TYPE xc_f90_pointer_t

END MODULE xc_f90_types_m
