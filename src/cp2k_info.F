!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief some minimal info about CP2K, including its version and license
!> \par History
!>      - created (2007-09, Joost VandeVondele)
!> \author Joost VandeVondele
! *****************************************************************************
MODULE cp2k_info
  USE input_constants,                 ONLY: id_development_version,&
                                             id_release_version

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: cp2k_year,cp2k_version,cp2k_home,id_cp2k_version
  PUBLIC :: compile_arch,compile_date,compile_host,compile_lastcvs
  PUBLIC :: print_cp2k_license

  ! the version string of CP2K intended to be adjust after releases and branches
  LOGICAL, PARAMETER, PUBLIC  :: enable_color_tags=.FALSE.
#if defined(__RELEASE_VERSION)
  INTEGER, PARAMETER          :: id_cp2k_version = id_release_version     ! (Uncomment for release     version)
#else
  INTEGER, PARAMETER          :: id_cp2k_version = id_development_version ! (Uncomment for development version)
#endif

  CHARACTER(LEN=*), PARAMETER :: version_nr      = "2.0.0"
  CHARACTER(LEN=*), PARAMETER :: cp2k_year       = "2008"
  CHARACTER(LEN=*), PARAMETER :: cp2k_version    = "CP2K version "//TRIM(version_nr)
  CHARACTER(LEN=*), PARAMETER :: cp2k_home       = "http://cp2k.berlios.de/"

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

CONTAINS

! *****************************************************************************
  SUBROUTINE print_cp2k_license(iunit)
    INTEGER                                  :: iunit

     WRITE(iunit,'(T2,A)') '!-----------------------------------------------------------------------------!'
     WRITE(iunit,'(T2,A)') '!                                                                             !'
     WRITE(iunit,'(T2,A)') '!   CP2K: A general program to perform molecular dynamics simulations         !'
     WRITE(iunit,'(T2,A)') '!   Copyright (C) 2000, 2001, 2002, 2003  CP2K developers group               !'
     WRITE(iunit,'(T2,A)') '!   Copyright (C) 2004, 2005, 2006, 2007  CP2K developers group               !'
     WRITE(iunit,'(T2,A)') '!   Copyright (C) 2008  CP2K developers group                                 !'
     WRITE(iunit,'(T2,A)') '!                                                                             !'
     WRITE(iunit,'(T2,A)') '!   This program is free software; you can redistribute it and/or modify      !'
     WRITE(iunit,'(T2,A)') '!   it under the terms of the GNU General Public License as published by      !'
     WRITE(iunit,'(T2,A)') '!   the Free Software Foundation; either version 2 of the License, or         !'
     WRITE(iunit,'(T2,A)') '!   (at your option) any later version.                                       !'
     WRITE(iunit,'(T2,A)') '!                                                                             !'
     WRITE(iunit,'(T2,A)') '!   This program is distributed in the hope that it will be useful,           !'
     WRITE(iunit,'(T2,A)') '!   but WITHOUT ANY WARRANTY; without even the implied warranty of            !'
     WRITE(iunit,'(T2,A)') '!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             !'
     WRITE(iunit,'(T2,A)') '!   GNU General Public License for more details.                              !'
     WRITE(iunit,'(T2,A)') '!                                                                             !'
     WRITE(iunit,'(T2,A)') '!   You should have received a copy of the GNU General Public License         !'
     WRITE(iunit,'(T2,A)') '!   along with this program; if not, write to the Free Software               !'
     WRITE(iunit,'(T2,A)') '!   Foundation, Inc., 51 Franklin Street, Fifth Floor,                        !'
     WRITE(iunit,'(T2,A)') '!   Boston, MA  02110-1301, USA.                                              !'
     WRITE(iunit,'(T2,A)') '!                                                                             !'
     WRITE(iunit,'(T2,A)') '!   See also http://www.fsf.org/licensing/licenses/gpl.html                   !'
     WRITE(iunit,'(T2,A)') '!                                                                             !'
     WRITE(iunit,'(T2,A)') '!-----------------------------------------------------------------------------!'
     WRITE(iunit,'(T2,A)') '!   CP2K, including its sources and pointers to the authors                   !'
     WRITE(iunit,'(T2,A)') '!   can be found at  http://cp2k.berlios.de/                                  !'
     WRITE(iunit,'(T2,A)') '!-----------------------------------------------------------------------------!'

  END SUBROUTINE print_cp2k_license

END MODULE

