!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/greens_fns [1.0] *
!!
!!   NAME
!!     greens_fns
!!
!!   FUNCTION
!!     initializations of the green functions
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     06.2004 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE greens_fns
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE greens_fn,                       ONLY: ANALYTIC0D,&
                                             ANALYTIC1D,&
                                             ANALYTIC2D,&
                                             HOCKNEY0D,&
                                             HOCKNEY1D,&
                                             HOCKNEY2D,&
                                             MT0D,&
                                             MT1D,&
                                             MT2D,&
                                             PERIODIC3D,&
                                             greens_function_type,&
                                             influence_factor
  USE kinds,                           ONLY: dp
  USE mt_util,                         ONLY: MTin_create_screen_fn
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: REALDATA1D,&
                                             RECIPROCALSPACE
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: pw_green_fn_rebuild
  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='greens_fns'

!!***
!****************************************************************************
CONTAINS

!!****** greens_fn/pw_green_fn_rebuild [1.0] *
!!
!!   NAME
!!     pw_green_fn_rebuild
!!
!!   FUNCTION
!!     reinitializes the grid dependent parts
!!
!!   AUTHOR
!!     JGH (12-Mar-2001)
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE pw_green_fn_rebuild ( gftype, pw_pool, p3m, alpha, error )


    TYPE(greens_function_type), &
      INTENT(INOUT),TARGET                   :: gftype
    TYPE(pw_pool_type), POINTER              :: pw_pool
    INTEGER, INTENT(IN), OPTIONAL            :: p3m
    REAL(KIND=dp), INTENT(IN), OPTIONAL          :: alpha
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = "pw_green_fn_rebuild", &
      routineP = moduleN//':'//routineN

!------------------------------------------------------------------------------

    CALL pw_pool_give_back_pw(pw_pool,gftype % p3m_charge, accept_non_compatible=.TRUE.)
    CALL pw_pool_give_back_pw(pw_pool,gftype % influence_function, accept_non_compatible=.TRUE.)
    CALL pw_pool_give_back_pw(pw_pool,gftype % screen_function, accept_non_compatible=.TRUE.)
    SELECT CASE ( gftype % method )
    CASE ( PERIODIC3D, ANALYTIC2D, ANALYTIC1D, ANALYTIC0D, MT2D, MT1D, MT0D )
       CALL pw_pool_create_pw ( pw_pool, gftype % influence_function,&
            use_data=REALDATA1D, in_space=RECIPROCALSPACE )
       IF ( PRESENT ( p3m ) ) THEN
          gftype % p3m = .TRUE.
          gftype % p3m_order = p3m
          IF ( PRESENT ( alpha ) ) THEN
             gftype % p3m_alpha = alpha
          ELSE
             CALL stop_program ( "pw_green_fn_rebuild", "P3M needs alpha PARAMETER")
          END IF
          CALL pw_pool_create_pw ( pw_pool, gftype % p3m_charge, use_data=REALDATA1D, &
               in_space=RECIPROCALSPACE)
          CALL influence_factor ( gftype )
       ELSE
          gftype % p3m = .FALSE.
       END IF
       !
       SELECT CASE( gftype % method )
       CASE(MT0D,MT1D,MT2D)
          CALL MTin_create_screen_fn(gftype%screen_function,pw_pool=pw_pool,method=gftype%method,&
               alpha=gftype%alpha)
       END SELECT
       !
    CASE ( HOCKNEY2D, HOCKNEY1D, HOCKNEY0D )
       CALL stop_program ( "pw_green_fn_rebuild", " Method not implemented ")
    CASE DEFAULT
       CALL stop_program ( "pw_green_fn_rebuild", " Method not implemented ")
    END SELECT

  END SUBROUTINE pw_green_fn_rebuild
!!*****
!******************************************************************************

END MODULE greens_fns
