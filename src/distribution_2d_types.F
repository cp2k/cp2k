!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2003 CP2K developers group
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* distribution_2d_types [1.0] *
!!
!!   NAME
!!     distribution_2d_types
!!
!!   FUNCTION
!!     stores a mapping of 2D info (e.g. matrix) on a 
!!     2D processor distribution (i.e. blacs grid)
!!     where cpus in the same blacs row own the same rows of the 2D info 
!!     (and similar for the cols)
!!
!!   NOTES
!!
!!   AUTHOR
!!     Joost VandeVondele (2003-08)
!!
MODULE distribution_2d_types
  
  USE kinds,                           ONLY: int_size,&
                                             wp => dp,&
                                             wp_size => dp_size

  USE cp_para_types,                   ONLY: cp_blacs_env_type
  USE cp_blacs_env,                    ONLY: cp_blacs_env_retain, &
                                             cp_blacs_env_release

  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error,&
                                             cp_unimplemented_error
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_logger_generate_filename,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level

  USE timings,                         ONLY: timeset,&
                                             timestop

  USE termination,                     ONLY: stop_memory

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_distribution_2d_type'

  TYPE distribution_2d_type
    INTEGER, DIMENSION(:), POINTER       :: row_distribution
    INTEGER, DIMENSION(:), POINTER       :: col_distribution
    INTEGER, DIMENSION(:), POINTER       :: local_rows
    INTEGER, DIMENSION(:), POINTER       :: local_cols
    TYPE(cp_blacs_env_type), POINTER     :: blacs_env
    INTEGER                              :: ref_count
  END TYPE distribution_2d_type

  PUBLIC :: distribution_2d_type

  PUBLIC :: distribution_2d_create, &
            distribution_2d_release, &
            distribution_2d_retain

CONTAINS

!!****f* distribution_2d_types/distribution_2d_create [1.0] *
!!
!!   NAME
!!     distribution_2d_create
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     initializes the distribution_2d by retaining info and computing local info
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Joost VandeVondele
!!
!!   MODIFICATION HISTORY
!!
!!*** *********************************************************************
SUBROUTINE distribution_2d_create(distribution_2d, row_distribution, col_distribution, blacs_env,error)
   TYPE(distribution_2d_type), POINTER :: distribution_2d
   INTEGER, DIMENSION(:), POINTER      :: row_distribution, col_distribution
   TYPE(cp_blacs_env_type), POINTER    :: blacs_env
   TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

   CHARACTER(len=*), PARAMETER :: routineN='distribution_2d_create',&
                                  routineP=moduleN//':'//routineN

   LOGICAL :: failure
   INTEGER :: istat

   failure=.FALSE.
   CPPrecondition(ASSOCIATED(row_distribution),cp_failure_level,routineP,error,failure)
   CPPrecondition(ASSOCIATED(col_distribution),cp_failure_level,routineP,error,failure)
   CPPrecondition(ASSOCIATED(blacs_env),cp_failure_level,routineP,error,failure)
   CPPrecondition(.NOT.ASSOCIATED(distribution_2d),cp_failure_level,routineP,error,failure)

   IF (.NOT. failure) THEN
      ALLOCATE(distribution_2d,stat=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      IF (.NOT. failure) THEN
         distribution_2d%ref_count = 0
         NULLIFY(distribution_2d%col_distribution)
         NULLIFY(distribution_2d%row_distribution)
         NULLIFY(distribution_2d%local_rows)
         NULLIFY(distribution_2d%local_cols)
         distribution_2d%blacs_env => blacs_env
         CALL cp_blacs_env_retain(distribution_2d%blacs_env)
         ALLOCATE(distribution_2d%row_distribution(SIZE(row_distribution)),stat=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         ALLOCATE(distribution_2d%col_distribution(SIZE(col_distribution)),stat=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         IF (.NOT. failure) THEN
             distribution_2d%row_distribution=row_distribution
             distribution_2d%col_distribution=col_distribution
         ENDIF ! failure
      ENDIF ! failure
   ENDIF ! failure

   IF(.NOT. failure) distribution_2d%ref_count=1

END SUBROUTINE distribution_2d_create

!!****f* distribution_2d_types/distribution_2d_retain [1.0] *
!!
!!   NAME
!!     distribution_2d_retain
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Joost VandeVondele
!!
!!   MODIFICATION HISTORY
!!
!!*** *********************************************************************
SUBROUTINE distribution_2d_retain(distribution_2d,error)
   TYPE(distribution_2d_type), POINTER :: distribution_2d
   TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

   CHARACTER(len=*), PARAMETER :: routineN='distribution_2d_retain',&
                                  routineP=moduleN//':'//routineN

   IF (ASSOCIATED(distribution_2d)) THEN
      distribution_2d%ref_count=distribution_2d%ref_count+1
   ENDIF
END SUBROUTINE distribution_2d_retain

!!****f* distribution_2d_types/distribution_2d_release [1.0] *
!!
!!   NAME
!!     distribution_2d_release
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!*** *********************************************************************
SUBROUTINE distribution_2d_release(distribution_2d,error)
   TYPE(distribution_2d_type), POINTER :: distribution_2d
   TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

   CHARACTER(len=*), PARAMETER :: routineN='distribution_2d_release',&
                                  routineP=moduleN//':'//routineN
   LOGICAL :: failure
   INTEGER :: istat

   IF (ASSOCIATED(distribution_2d)) THEN
      failure=.FALSE.
      CPPrecondition(distribution_2d%ref_count>0,cp_failure_level,routineP,error,failure)
      IF (.NOT. failure ) THEN
         distribution_2d%ref_count=distribution_2d%ref_count-1
         IF (distribution_2d%ref_count == 0 ) THEN
            CALL cp_blacs_env_release(distribution_2d%blacs_env)
            DEALLOCATE(distribution_2d%col_distribution,stat=istat)
            CPPostcondition(istat==0,cp_warning_level,routineP,error,failure)
            DEALLOCATE(distribution_2d%row_distribution,stat=istat)
            CPPostcondition(istat==0,cp_warning_level,routineP,error,failure)
            DEALLOCATE(distribution_2d%local_rows,stat=istat)
            CPPostcondition(istat==0,cp_warning_level,routineP,error,failure)
            DEALLOCATE(distribution_2d%local_cols,stat=istat)
            CPPostcondition(istat==0,cp_warning_level,routineP,error,failure)
            DEALLOCATE(distribution_2d,stat=istat)
            CPPostcondition(istat==0,cp_warning_level,routineP,error,failure)
         ENDIF
      ENDIF
   ENDIF
END SUBROUTINE distribution_2d_release

END MODULE distribution_2d_types
