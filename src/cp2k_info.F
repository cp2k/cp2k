!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
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
  
  USE input_constants,                 ONLY: id_development_version,&
                                             id_release_version
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length
  USE machine,                         ONLY: m_datum,&
                                             m_getcwd,&
                                             m_getlog,&
                                             m_getpid,&
                                             m_hostnm
  USE string_utilities,                ONLY: integer_to_string

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: cp2k_year,cp2k_version,cp2k_home,id_cp2k_version, cp2k_flags
  PUBLIC :: compile_arch,compile_date,compile_host,compile_revision
  PUBLIC :: print_cp2k_license, get_runtime_info, write_restart_header

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

  CHARACTER(LEN=*), PARAMETER :: version_nr      = "2.7"
  CHARACTER(LEN=*), PARAMETER :: cp2k_year       = "2015"
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
!> \brief list all compile time options that influence the capabilities of cp2k.
!>        All new flags should be added here (and be unique grep-able)
!> \retval flags ...
! *****************************************************************************
  FUNCTION cp2k_flags() RESULT(flags)
    CHARACTER(len=10*default_string_length)  :: flags

    CHARACTER(len=default_string_length)     :: tmp_str

    flags="cp2kflags:"
!$  flags=TRIM(flags)//" omp"

#if defined(__LIBINT)
    flags=TRIM(flags)//" libint"
#endif
#if defined(__FFTW3)
    flags=TRIM(flags)//" fftw3"
#endif
#if defined(__LIBXC2)
    flags=TRIM(flags)//" libxc2"
#endif
#if defined(__LIBPEXSI)
    flags=TRIM(flags)//" pexsi"
#endif
#if defined(__ELPA)
    flags=TRIM(flags)//" elpa"
#endif
#if defined(__parallel)
    flags=TRIM(flags)//" parallel"
#endif
#if defined(__SCALAPACK)
    flags=TRIM(flags)//" scalapack"
#endif

#if defined(__QUIP)
    flags=TRIM(flags)//" quip"
#endif
#if defined(__PLUMED_CP2K)
    flags=TRIM(flags)//" plumed"
#endif

#if defined(__ACC)
    flags=TRIM(flags)//" acc"
#endif
#if defined(__PW_CUDA)
    flags=TRIM(flags)//" pw_cuda"
#endif

#if defined __GEMINI || __SEASTAR || __BLUEGENE || __NET || __HWLOC || __LIBNUMA
    flags=TRIM(flags)//" macharch"
#endif


#if defined __HAS_smm_vec
   flags=TRIM(flags)//" smm_vec"
#endif
#if defined __HAS_smm_snn
   flags=TRIM(flags)//" smm_snn"
#endif
#if defined __HAS_smm_snt
   flags=TRIM(flags)//" smm_snt"
#endif
#if defined __HAS_smm_stn
   flags=TRIM(flags)//" smm_stn"
#endif
#if defined __HAS_smm_stt
   flags=TRIM(flags)//" smm_stt"
#endif
#if defined __HAS_smm_znn
   flags=TRIM(flags)//" smm_znn"
#endif
#if defined __HAS_smm_znt
   flags=TRIM(flags)//" smm_znt"
#endif
#if defined __HAS_smm_ztn
   flags=TRIM(flags)//" smm_ztn"
#endif
#if defined __HAS_smm_ztt
   flags=TRIM(flags)//" smm_ztt"
#endif
#if defined __HAS_smm_cnn
   flags=TRIM(flags)//" smm_cnn"
#endif
#if defined __HAS_smm_cnt
   flags=TRIM(flags)//" smm_cnt"
#endif
#if defined __HAS_smm_ctn
   flags=TRIM(flags)//" smm_ctn"
#endif
#if defined __HAS_smm_ctt
   flags=TRIM(flags)//" smm_ctt"
#endif
#if defined __HAS_smm_dnn
   flags=TRIM(flags)//" smm_dnn"
#endif
#if defined __HAS_smm_dnt
   flags=TRIM(flags)//" smm_dnt"
#endif
#if defined __HAS_smm_dtn
   flags=TRIM(flags)//" smm_dtn"
#endif
#if defined __HAS_smm_dtt
   flags=TRIM(flags)//" smm_dtt"
#endif
IF(INDEX(flags, " smm_") > 0) THEN
   flags=TRIM(flags)//" smm"
ENDIF

#if defined __CRAY_PM_ACCEL_ENERGY
   flags=TRIM(flags)//" cray_pm_accel_energy"
#endif
#if defined __CRAY_PM_ENERGY
   flags=TRIM(flags)//" cray_pm_energy"
#endif
#if defined __CRAY_PM_FAKE_ENERGY
   flags=TRIM(flags)//" cray_pm_fake_energy"
#endif
#if defined __CUBLASDP
   flags=TRIM(flags)//" cublasdp"
#endif
#if defined __CUDA_PROFILING
   flags=TRIM(flags)//" cuda_profiling"
#endif
#if defined __DBCSR_ACC
   flags=TRIM(flags)//" dbcsr_acc"
#endif
#if defined __DEBUG_HWTOPO
   flags=TRIM(flags)//" debug_hwtopo"
#endif
#if defined __HAS_LIBGRID
   flags=TRIM(flags)//" libgrid"
#endif
#if defined __HAS_NO_CUDA_STREAM_PRIORITIES
   flags=TRIM(flags)//" has_no_cuda_stream_priorities"
#endif
#if defined __HAS_NO_MPI_MOD
   flags=TRIM(flags)//" has_no_mpi_mod"
#endif
#if defined __HAS_NO_SHARED_GLIBC
   flags=TRIM(flags)//" has_no_shared_glibc"
#endif
#if defined __LIBDERIV_MAX_AM1
   CALL integer_to_string(__LIBDERIV_MAX_AM1, tmp_str)
   flags=TRIM(flags)//" libderiv_max_am1="//TRIM(tmp_str)
#endif
#if defined __LIBINT_MAX_AM
   CALL integer_to_string(__LIBINT_MAX_AM, tmp_str)
   flags=TRIM(flags)//" libint_max_am="//TRIM(tmp_str)
#endif
#if defined __LINUMA
   flags=TRIM(flags)//" linuma"
#endif
#if defined __MAX_CONTR
   CALL integer_to_string(__MAX_CONTR, tmp_str)
   flags=TRIM(flags)//" max_contr="//TRIM(tmp_str)
#endif
#if defined __NO_IPI_DRIVER
   flags=TRIM(flags)//" no_ipi_driver"
#endif
#if defined __NO_MPI_THREAD_SUPPORT_CHECK
   flags=TRIM(flags)//" no_mpi_thread_support_check"
#endif
#if defined __NO_STATM_ACCESS
   flags=TRIM(flags)//" no_statm_access"
#endif
#if defined __PW_CUDA_NO_HOSTALLOC
   flags=TRIM(flags)//" pw_cuda_no_hostalloc"
#endif
#if defined __STATM_RESIDENT
   flags=TRIM(flags)//" statm_resident"
#endif
#if defined __STATM_TOTAL
   flags=TRIM(flags)//" statm_total"
#endif
#if defined __USE_LEGACY_WEIGHTS
   flags=TRIM(flags)//" use_legacy_weights"
#endif

  END FUNCTION cp2k_flags

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
     WRITE(iunit,'(T2,A)') '!   Copyright (C) 2012, 2013, 2014, 2015  CP2K developers group               !'
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

  END SUBROUTINE


! *****************************************************************************
!> \brief Writes the header for the restart file
!> \param iunit ...
!> \par History
!>      01.2008 [created] - Splitted from write_restart
!> \author Teodoro Laino - University of Zurich - 01.2008
! *****************************************************************************
  SUBROUTINE write_restart_header(iunit)
    INTEGER, INTENT(IN)                      :: iunit

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_restart_header', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=256)                       :: cwd, datx

    CALL m_datum(datx)
    CALL m_getcwd(cwd)

    WRITE (UNIT=iunit,FMT="(T2,A)") "# Version information for this restart file "
    WRITE (UNIT=iunit,FMT="(T2,A)") "# current date "//TRIM(datx)
    WRITE (UNIT=iunit,FMT="(T2,A)") "# current working dir "//TRIM(cwd)

    WRITE (UNIT=iunit,FMT="(T2,A,T31,A50)")&
         "# Program compiled at",&
         ADJUSTR(compile_date(1:MIN(50,LEN(compile_date))))
    WRITE (UNIT=iunit,FMT="(T2,A,T31,A50)")&
         "# Program compiled on",&
         ADJUSTR(compile_host(1:MIN(50,LEN(compile_host))))
    WRITE (UNIT=iunit,FMT="(T2,A,T31,A50)")&
         "# Program compiled for",&
         ADJUSTR(compile_arch(1:MIN(50,LEN(compile_arch))))
    WRITE (UNIT=iunit,FMT="(T2,A,T31,A50)")&
         "# Source code revision number",&
         ADJUSTR(compile_revision)

  END SUBROUTINE write_restart_header


END MODULE cp2k_info

