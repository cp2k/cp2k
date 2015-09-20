!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief represent the global information of a run: the input file,
!>      parallel environment (and partially output).
!>      Use sparingly, try not to send it too deep in your structures.
!> \par History
!>      - print keys, basis_set_file name and potential_file_name added to the
!>        global type (27.02.2001, MK)
!>      - JGH (28.11.2001) : added pp_library_path to type
!>      - Merged with MODULE print_keys (17.01.2002, MK)
!>      -  reference counting, create (08.2004, fawzi)
!>      - new (parallel) random number generator (11.03.06,MK)
!> \author JGH,MK,fawzi
! *****************************************************************************
MODULE global_types

  USE cp_blacs_env,                    ONLY: BLACS_GRID_SQUARE
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE machine,                         ONLY: m_walltime
  USE parallel_rng_types,              ONLY: delete_rng_stream,&
                                             rng_stream_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  ! Global parameters

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'global_types'

  INTEGER, PARAMETER :: SILENT = 0,&
                        LOW    = 1,&
                        MEDIUM = 2,&
                        HIGH   = 3,&
                        DEBUG  = 4

  INTEGER, SAVE :: last_globenv_id=0

  ! Public parameters


  ! Public data types
  PUBLIC :: global_environment_type

  ! Public subroutines
  PUBLIC :: globenv_create,&
            globenv_retain,&
            globenv_release

! *****************************************************************************
!> \brief contains the initally parsed file and the initial parallel environment
!> \param id_nr identification number (unique)
!> \param ref_count reference count (see doc/ReferenceCounting.html)
!> \param handle handle with the total time of the computation
!>
!>      Personally I think that all the other attributes should go away
!>      (and maybe add  a logger)[fawzi]
!> \note
!>      This is not but really should be passed as pointer and use reference
!>      counting. Use it accordingly wherever possible.
! *****************************************************************************
  TYPE global_environment_type
     INTEGER :: id_nr, ref_count
     TYPE(rng_stream_type), POINTER          :: gaussian_rng_stream
     CHARACTER(LEN=default_string_length)    :: diag_library
     CHARACTER(LEN=default_string_length)    :: default_fft_library
     CHARACTER(LEN=default_path_length)      :: fftw_wisdom_file_name

     INTEGER :: fft_pool_scratch_limit !! limit used for fft scratches
     INTEGER :: fftw_plan_type !! which kind of planning to use with fftw
     INTEGER :: idum         !! random number seed
     INTEGER :: prog_name_id !! index to define the type of program
     INTEGER :: run_type_id  !! index to define the run_tupe
     INTEGER :: blacs_grid_layout !! will store the user preference for the blacs grid
     INTEGER :: k_elpa !! optimized kernel for the ELPA diagonalization library
     LOGICAL :: blacs_repeatable !! will store the user preference for the repeatability of blacs collectives
     REAL(KIND=dp) :: cp2k_start_time, cp2k_target_time
     INTEGER :: handle
  END TYPE global_environment_type

CONTAINS

! *****************************************************************************
!> \brief creates a globenv
!> \param globenv the globenv to create
!> \author fawzi
! *****************************************************************************
  SUBROUTINE globenv_create(globenv)
    TYPE(global_environment_type), POINTER   :: globenv

    CHARACTER(len=*), PARAMETER :: routineN = 'globenv_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPASSERT(.NOT.ASSOCIATED(globenv))
    ALLOCATE(globenv,stat=stat)
    CPASSERT(stat==0)
    last_globenv_id=last_globenv_id+1
    globenv%id_nr=last_globenv_id
    globenv%ref_count=1
    globenv%run_type_id=0
    globenv%diag_library="SL"
    globenv%k_elpa=1
    globenv%default_fft_library="FFTSG"
    globenv%fftw_wisdom_file_name="/etc/fftw/wisdom"
    globenv%prog_name_id=0
    globenv%idum=0                 !! random number seed
    globenv%blacs_grid_layout=BLACS_GRID_SQUARE
    globenv%cp2k_start_time=m_walltime()
    NULLIFY (globenv%gaussian_rng_stream)
  END SUBROUTINE globenv_create

! *****************************************************************************
!> \brief retains the global environment
!> \param globenv the global environment to retain
!> \author fawzi
! *****************************************************************************
SUBROUTINE globenv_retain(globenv)
    TYPE(global_environment_type), POINTER   :: globenv

    CHARACTER(len=*), PARAMETER :: routineN = 'globenv_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPASSERT(ASSOCIATED(globenv))
  CPASSERT(globenv%ref_count>0)
  globenv%ref_count=globenv%ref_count+1
END SUBROUTINE globenv_retain

! *****************************************************************************
!> \brief releases the global environment
!> \param globenv the global environment to release
!> \author fawzi
! *****************************************************************************
SUBROUTINE globenv_release(globenv)
    TYPE(global_environment_type), POINTER   :: globenv

    CHARACTER(len=*), PARAMETER :: routineN = 'globenv_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.

  IF (ASSOCIATED(globenv)) THEN
     CPASSERT(globenv%ref_count>0)
     globenv%ref_count=globenv%ref_count-1
     IF (globenv%ref_count==0) THEN
        IF (ASSOCIATED(globenv%gaussian_rng_stream)) THEN
          CALL delete_rng_stream(globenv%gaussian_rng_stream)
        END IF
        DEALLOCATE(globenv,stat=stat)
        CPASSERT(stat==0)
     END IF
  END IF
  NULLIFY(globenv)
END SUBROUTINE globenv_release

END MODULE global_types
