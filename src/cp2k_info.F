!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief some minimal info about CP2K, including its version and license
!> \par History
!>      - created (2007-09, Joost VandeVondele)
!>      - moved into this module information related to runtime:pid, user_name,
!>        host_name, cwd, datx  (2009-06, Teodoro Laino)
!> \author Joost VandeVondele
! *****************************************************************************
MODULE cp2k_info
  USE f77_blas
  USE input_constants,                 ONLY: id_development_version,&
                                             id_release_version
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length
  USE machine,                         ONLY: m_datum,&
                                             m_getcwd,&
                                             m_getlog,&
                                             m_getpid,&
                                             m_hostnm

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: cp2k_year,cp2k_version,cp2k_home,id_cp2k_version
  PUBLIC :: compile_arch,compile_date,compile_host,compile_revision
  PUBLIC :: print_cp2k_license, get_runtime_info

  ! the version string of CP2K intended to be adjust after releases and branches
#if defined(__RELEASE_VERSION)
  INTEGER, PARAMETER          :: id_cp2k_version = id_release_version     ! (Uncomment for release     version)
#else
  INTEGER, PARAMETER          :: id_cp2k_version = id_development_version ! (Uncomment for development version)
! INTEGER, PARAMETER          :: id_cp2k_version = 2                      ! (Uncomment for branch      version)
#endif

#if defined(__COMPILE_REVISION)
  CHARACTER(LEN=*), PARAMETER :: compile_revision = __COMPILE_REVISION
#else
  CHARACTER(LEN=*), PARAMETER :: compile_revision = "unknown"
#endif

  CHARACTER(LEN=*), PARAMETER :: version_nr      = "2.6"
  CHARACTER(LEN=*), PARAMETER :: cp2k_year       = "2014"
  CHARACTER(LEN=*), PARAMETER :: cp2k_version    = "CP2K version "//TRIM(version_nr)
  CHARACTER(LEN=*), PARAMETER :: cp2k_home       = "http://www.cp2k.org/"

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

  ! Local runtime informations
  CHARACTER(LEN=26), PUBLIC                        :: r_datx
  CHARACTER(LEN=default_path_length), PUBLIC       :: r_cwd
  CHARACTER(LEN=default_string_length), PUBLIC     :: r_host_name, r_user_name
  INTEGER, PUBLIC                                  :: r_pid

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp2k_info'
CONTAINS

! *****************************************************************************
!> \brief ...
!> \param iunit ...
! *****************************************************************************
  SUBROUTINE print_cp2k_license(iunit)
    INTEGER                                  :: iunit

     WRITE(iunit,'(T2,A)') '!-----------------------------------------------------------------------------!'
     WRITE(iunit,'(T2,A)') '!                                                                             !'
     WRITE(iunit,'(T2,A)') '!   CP2K: A general program to perform molecular dynamics simulations         !'
     WRITE(iunit,'(T2,A)') '!   Copyright (C) 2000, 2001, 2002, 2003  CP2K developers group               !'
     WRITE(iunit,'(T2,A)') '!   Copyright (C) 2004, 2005, 2006, 2007  CP2K developers group               !'
     WRITE(iunit,'(T2,A)') '!   Copyright (C) 2008, 2009, 2010, 2011  CP2K developers group               !'
     WRITE(iunit,'(T2,A)') '!   Copyright (C) 2012, 2013, 2014        CP2K developers group               !'
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
     WRITE(iunit,'(T2,A)') '!   can be found at  http://www.cp2k.org/                                     !'
     WRITE(iunit,'(T2,A)') '!-----------------------------------------------------------------------------!'

  END SUBROUTINE print_cp2k_license

! *****************************************************************************
!> \brief ...
! *****************************************************************************
  SUBROUTINE get_runtime_info

    r_datx = ""
    r_cwd  = ""
    r_host_name = ""
    r_user_name = ""
    r_pid = -1

    CALL m_getpid(r_pid)
    CALL m_getlog(r_user_name)
    CALL m_hostnm(r_host_name)
    CALL m_datum(r_datx)
    CALL m_getcwd(r_cwd)

  END SUBROUTINE get_runtime_info

END MODULE cp2k_info

