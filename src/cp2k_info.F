!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/cp2k_info [1.0] *
!!
!!   NAME
!!     cp2k_info
!!
!!   FUNCTION
!!     some minimal info about CP2K, including its version and license
!!
!!   AUTHOR
!!     Joost VandeVondele
!!
!!   MODIFICATION HISTORY
!!     - created (2007-09, Joost VandeVondele)
!!
!!******
!******************************************************************************

MODULE cp2k_info
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: cp2k_year,cp2k_version,cp2k_home
  PUBLIC :: compile_arch,compile_date,compile_host,compile_lastcvs
  PUBLIC :: print_cp2k_license

  ! the version string of CP2K intended to be adjust after releases and branches
  CHARACTER(LEN=*), PARAMETER :: cp2k_year    = "2007"
  CHARACTER(LEN=*), PARAMETER :: cp2k_version = "CP2K version 2.0.0 (development version)"
  CHARACTER(LEN=*), PARAMETER :: cp2k_home    = "http://cp2k.berlios.de/"

  ! compile time information
#if defined(__COMPILE_ARCH)
  CHARACTER(LEN=*), PARAMETER :: compile_arch = __COMPILE_ARCH
#else
  CHARACTER(LEN=*), PARAMETER :: compile_arch = "unknown: -D__COMPILE_ARCH=?"
#endif

#if defined(__COMPILE_DATE)
  CHARACTER(LEN=*), PARAMETER :: compile_date = __COMPILE_DATE
#else
  CHARACTER(LEN=*), PARAMETER :: compile_date = "unknown: -D__COMPILE_DATE=?"
#endif

#if defined(__COMPILE_HOST)
  CHARACTER(LEN=*), PARAMETER :: compile_host = __COMPILE_HOST
#else
  CHARACTER(LEN=*), PARAMETER :: compile_host = "unknown: -D__COMPILE_HOST=?"
#endif

#if defined(__COMPILE_LASTCVS)
  CHARACTER(LEN=*), PARAMETER :: compile_lastcvs = __COMPILE_LASTCVS
#else
  CHARACTER(LEN=*), PARAMETER :: compile_lastcvs = "/No CVS/Entry in src dir//"
#endif

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp2k_info'

!******************************************************************************

CONTAINS

  SUBROUTINE print_cp2k_license(iunit)
     INTEGER :: iunit 

     write(iunit,'(T2,A)') '!-----------------------------------------------------------------------------!'
     write(iunit,'(T2,A)') '!                                                                             !'
     write(iunit,'(T2,A)') '!   CP2K: A general program to perform molecular dynamics simulations         !'
     write(iunit,'(T2,A)') '!   Copyright (C) 2000, 2001, 2002, 2003  CP2K developers group               !'
     write(iunit,'(T2,A)') '!   Copyright (C) 2004, 2005, 2006, 2007  CP2K developers group               !'
     write(iunit,'(T2,A)') '!                                                                             !'
     write(iunit,'(T2,A)') '!   This program is free software; you can redistribute it and/or modify      !'
     write(iunit,'(T2,A)') '!   it under the terms of the GNU General Public License as published by      !'
     write(iunit,'(T2,A)') '!   the Free Software Foundation; either version 2 of the License, or         !'
     write(iunit,'(T2,A)') '!   (at your option) any later version.                                       !'
     write(iunit,'(T2,A)') '!                                                                             !'
     write(iunit,'(T2,A)') '!   This program is distributed in the hope that it will be useful,           !'
     write(iunit,'(T2,A)') '!   but WITHOUT ANY WARRANTY; without even the implied warranty of            !'
     write(iunit,'(T2,A)') '!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             !'
     write(iunit,'(T2,A)') '!   GNU General Public License for more details.                              !'
     write(iunit,'(T2,A)') '!                                                                             !'
     write(iunit,'(T2,A)') '!   You should have received a copy of the GNU General Public License         !'
     write(iunit,'(T2,A)') '!   along with this program; if not, write to the Free Software               !'
     write(iunit,'(T2,A)') '!   Foundation, Inc., 51 Franklin Street, Fifth Floor,                        !'
     write(iunit,'(T2,A)') '!   Boston, MA  02110-1301, USA.                                              !'
     write(iunit,'(T2,A)') '!                                                                             !'
     write(iunit,'(T2,A)') '!   See also http://www.fsf.org/licensing/licenses/gpl.html                   !'
     write(iunit,'(T2,A)') '!                                                                             !'
     write(iunit,'(T2,A)') '!-----------------------------------------------------------------------------!'
     write(iunit,'(T2,A)') '!   CP2K, including its sources and pointers to the authors                   !'
     write(iunit,'(T2,A)') '!   can be found at  http://cp2k.berlios.de/                                  !'
     write(iunit,'(T2,A)') '!-----------------------------------------------------------------------------!'

  END SUBROUTINE print_cp2k_license

END MODULE

