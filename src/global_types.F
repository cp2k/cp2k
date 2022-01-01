!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Define type storing the global information of a run. Keep the amount of stored data small.
!>        Use it sparingly and try not to send it too deep in your structures.
!> \par History
!>      - print keys, basis_set_file name and potential_file_name added to the
!>        global type (27.02.2001, MK)
!>      - added pp_library_path to type (28.11.2001, JGH)
!>      - Merged with MODULE print_keys (17.01.2002, MK)
!>      - reference counting, create (08.2004, fawzi)
!>      - new (parallel) random number generator (11.03.2006, MK)
!>      - add eps_check_diag, remove id_nr from globenv, and revise initialization (04.05.2021, MK)
!> \author JGH, MK, fawzi
! **************************************************************************************************
MODULE global_types

   USE cp_blacs_env,                    ONLY: BLACS_GRID_SQUARE
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE machine,                         ONLY: m_walltime
   USE parallel_rng_types,              ONLY: rng_stream_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ! Global parameters

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'global_types'

   INTEGER, PARAMETER :: SILENT = 0, &
                         LOW = 1, &
                         MEDIUM = 2, &
                         HIGH = 3, &
                         DEBUG = 4

   ! Public data types

   PUBLIC :: global_environment_type

   ! Public subroutines

   PUBLIC :: globenv_create, &
             globenv_retain, &
             globenv_release

! **************************************************************************************************
!> \brief contains the initially parsed file and the initial parallel environment
!> \param ref_count reference count (see doc/ReferenceCounting.html)
!> \param handle handle with the total time of the computation
!>
!>      Personally I think that all the other attributes should go away
!>      (and maybe add  a logger)[fawzi]
!> \note
!>      This is not but really should be passed as pointer and use reference
!>      counting. Use it accordingly wherever possible.
! **************************************************************************************************
   TYPE global_environment_type
      INTEGER :: ref_count = 0
      TYPE(rng_stream_type), ALLOCATABLE      :: gaussian_rng_stream
      CHARACTER(LEN=default_string_length)    :: diag_library = "ScaLAPACK"
      CHARACTER(LEN=default_string_length)    :: default_fft_library = "FFTSG"
      CHARACTER(LEN=default_path_length)      :: fftw_wisdom_file_name = "/etc/fftw/wisdom"

      INTEGER :: fft_pool_scratch_limit = 0 ! limit number of used FFT scratches
      INTEGER :: fftw_plan_type = 0 ! which kind of planning to use with FFTW
      INTEGER :: idum = 0 ! random number seed
      INTEGER :: prog_name_id = 0 ! index to define the type of program
      INTEGER :: run_type_id = 0 ! index to define the run_tupe
      INTEGER :: blacs_grid_layout = BLACS_GRID_SQUARE ! will store the user preference for the BLACS grid
      INTEGER :: k_elpa = 1 ! optimized kernel for the ELPA diagonalization library
      INTEGER :: elpa_neigvec_min = 0 ! Minimum number of eigenvectors for ELPA usage
      LOGICAL :: elpa_qr = .FALSE. ! allow ELPA to use QR during diagonalization
      LOGICAL :: elpa_print = .FALSE. ! if additional information about ELPA diagonalization should be printed
      LOGICAL :: elpa_qr_unsafe = .FALSE. ! enable potentially unsafe ELPA options
      LOGICAL :: blacs_repeatable = .FALSE. ! will store the user preference for the repeatability of BLACS collectives
      REAL(KIND=dp) :: cp2k_start_time = 0.0_dp
      REAL(KIND=dp) :: cp2k_target_time = HUGE(0.0_dp) ! Maximum run time in seconds
      ! Threshold value for the orthonormality of the eigenvectors after diagonalization
      ! A negative threshold value disables the check which is the default
      REAL(KIND=dp) :: eps_check_diag = -1.0_dp
      INTEGER :: handle = 0
   END TYPE global_environment_type

CONTAINS

! **************************************************************************************************
!> \brief Creates the global environment globenv
!> \param globenv the globenv to create
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE globenv_create(globenv)
      TYPE(global_environment_type), POINTER             :: globenv

      CPASSERT(.NOT. ASSOCIATED(globenv))
      ALLOCATE (globenv)
      globenv%ref_count = 1
      globenv%cp2k_start_time = m_walltime()

   END SUBROUTINE globenv_create

! **************************************************************************************************
!> \brief Retains the global environment globenv
!> \param globenv the global environment to retain
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE globenv_retain(globenv)
      TYPE(global_environment_type), POINTER             :: globenv

      CPASSERT(ASSOCIATED(globenv))
      CPASSERT(globenv%ref_count > 0)
      globenv%ref_count = globenv%ref_count + 1

   END SUBROUTINE globenv_retain

! **************************************************************************************************
!> \brief Releases the global environment globenv
!> \param globenv the global environment to release
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE globenv_release(globenv)
      TYPE(global_environment_type), POINTER             :: globenv

      IF (ASSOCIATED(globenv)) THEN
         CPASSERT(globenv%ref_count > 0)
         globenv%ref_count = globenv%ref_count - 1
         IF (globenv%ref_count == 0) THEN
            IF (ALLOCATED(globenv%gaussian_rng_stream)) &
               DEALLOCATE (globenv%gaussian_rng_stream)
            DEALLOCATE (globenv)
         END IF
      END IF

      NULLIFY (globenv)

   END SUBROUTINE globenv_release

END MODULE global_types
