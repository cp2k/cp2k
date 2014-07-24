!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Utilities for Molecular Dynamics
!> \author Teodoro Laino [tlaino] - University of Zurich - 09.2007
! *****************************************************************************
MODULE md_util
  USE f77_blas
  USE input_cp2k_restarts,             ONLY: write_restart
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE md_energies,                     ONLY: md_write_output
  USE md_environment_types,            ONLY: md_environment_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'md_util'

  PUBLIC :: md_output

CONTAINS

! *****************************************************************************
!> \brief collects the part of the MD that, basically, does the output
!> \param md_env ...
!> \param md_section ...
!> \param root_section ...
!> \param forced_io ...
!> \param error ...
!> \par History
!>      03.2006 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE md_output(md_env,md_section,root_section,forced_io,error)
    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(section_vals_type), POINTER         :: md_section, root_section
    LOGICAL, INTENT(IN)                      :: forced_io
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'md_output', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: do_print
    TYPE(section_vals_type), POINTER         :: print_section

    CALL timeset(routineN,handle)
    do_print = .TRUE.
    IF (forced_io) THEN
       print_section => section_vals_get_subs_vals(md_section,"PRINT",error=error)
       CALL section_vals_val_get(print_section,"FORCE_LAST",l_val=do_print,error=error)
    END IF
    IF (do_print) THEN
       ! Dumps all files related to the MD run
       CALL md_write_output(md_env, error)
       CALL write_restart(md_env=md_env,root_section=root_section, error=error)
    END IF
    CALL timestop(handle)

  END SUBROUTINE md_output

END MODULE md_util
