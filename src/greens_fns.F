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
  USE input_constants,                 ONLY: ANALYTIC0D,&
                                             ANALYTIC1D,&
                                             ANALYTIC2D,&
                                             HOCKNEY0D,&
                                             HOCKNEY1D,&
                                             HOCKNEY2D,&
                                             MT0D,&
                                             MT1D,&
                                             MT2D,&
                                             PERIODIC3D
  USE greens_fn,                       ONLY: greens_function_type,&
                                             influence_factor,&
                                             pw_green_create
  USE kinds,                           ONLY: dp,&
                                             default_string_length
  USE mt_util,                         ONLY: MTin_create_screen_fn
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: REALDATA1D,&
                                             RECIPROCALSPACE
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE input_cp2k,                      ONLY: parsed_cp2k_input

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: pw_green_fn_rebuild,&
            pw_green_setup
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


    TYPE(greens_function_type), pointer      :: gftype
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


!!****** greens_fn/pw_green_setup [1.0] *
!!
!!   NAME
!!     pw_green_setup
!!
!!   FUNCTION
!!     Setup the green function using input parameters
!!
!!   AUTHOR
!!     Teodoro Laino [teo] 03.2005
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE pw_green_setup ( gftype, poisson_section, perd, abc, error )
    IMPLICIT NONE
    ! Arguments
    TYPE(greens_function_type), POINTER       :: gftype
    TYPE(section_vals_type), POINTER          :: poisson_section
    INTEGER, DIMENSION(3), INTENT(IN)         :: perd
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)   :: abc
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                :: error
    ! Local Variables
    CHARACTER(LEN=default_string_length)      :: poisson_solver
    CHARACTER(len=*), PARAMETER :: routineN = "pw_green_setup", &
      routineP = moduleN//':'//routineN
    LOGICAL :: failure
    INTEGER :: dim, i

    failure = .FALSE.
    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.(ASSOCIATED(gftype)),cp_failure_level,routineP,error,failure)
       CALL pw_green_create(gftype, error=error)
       dim = COUNT(perd == 1)
       
       SELECT CASE (dim)
       CASE (0)
          gftype%method = ANALYTIC0D
          gftype%radius = 0.5_dp*MINVAL(abc)
       CASE (1)
          gftype%method = ANALYTIC1D
          gftype%special_dimension = MAXLOC(perd,1)
          gftype%radius = MAXVAL(abc)
          DO i=1,3
             IF (i == gftype%special_dimension) CYCLE
             gftype%radius = MIN(gftype%radius,0.5_dp*abc(i))
          END DO
       CASE (2)
          gftype%method = ANALYTIC2D
          i = MINLOC(perd,1)
          gftype%special_dimension = i
          gftype%slab_size = abc(i)
       CASE (3)
          gftype%method = PERIODIC3D
       CASE DEFAULT
          CALL stop_program(routineP,"Illegal value for perd found")
       END SELECT
       
       IF (parsed_cp2k_input(poisson_section,check_this_section=.TRUE.,error=error)) THEN
          CALL section_vals_val_get(poisson_section,"POISSON_SOLVER",c_val=poisson_solver,error=error)
          SELECT CASE (TRIM(poisson_solver))
          CASE ("PERIODIC")
             gftype%method = PERIODIC3D
             IF (dim /= 3) THEN
                CALL stop_program(routineP,&
                     "Illegal combination of periodicity and Poisson solver (perd)")
             END IF
          CASE ("ANALYTIC")
             SELECT CASE (dim)
             CASE (0)
                gftype%method = ANALYTIC0D
             CASE (1)
                gftype%method = ANALYTIC1D
             CASE (2)
                gftype%method = ANALYTIC2D
             CASE (3)
                gftype%method = PERIODIC3D
             CASE DEFAULT
                CALL stop_program(routineP,"Illegal value for perd found")
             END SELECT
          CASE ("HOCKNEY")
             SELECT CASE (dim )
             CASE (0)
                gftype%method = HOCKNEY0D
             CASE (1)
                gftype%method = HOCKNEY1D
             CASE (2)
                gftype%method = HOCKNEY2D
             CASE (3)
                CALL stop_program(routineP,&
                     "Illegal combination of periodicity and Poisson solver (HOCKNEY)")
             CASE DEFAULT
                CALL stop_program(routineP,"Illegal value for perd found")
             END SELECT
          CASE ("MT")
             SELECT CASE (dim)
             CASE (0)
                gftype%method = MT0D
             CASE (1)
                gftype%method = MT1D
             CASE (2)
                gftype%method = MT2D
             CASE (3)
                CALL stop_program(routineP,"Illegal combination of periodicity and Poisson solver (MT)")
             CASE DEFAULT
                CALL stop_program(routineP,"Illegal value for perd found")
             END SELECT
          CASE DEFAULT
             CALL stop_program(routineP,"An unknown Poisson solver was specified")
          END SELECT
       END IF
    END IF
  END SUBROUTINE pw_green_setup

END MODULE greens_fns
