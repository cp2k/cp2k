!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/green_methods [1.0] *
!!
!!   NAME
!!     green_methods
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
MODULE green_methods
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
                                             PERIODIC3D,&
                                             use_periodic,&
                                             use_analytic,&
                                             use_mt,&
                                             use_hockney,&
                                             use_perd_x,&
                                             use_perd_y,&
                                             use_perd_z,&
                                             use_perd_xy,&
                                             use_perd_xz,&
                                             use_perd_yz,&
                                             use_perd_xyz,&
                                             use_perd_none                                             
  USE green_types,                     ONLY: greens_function_type,&
                                             influence_factor,&
                                             pw_green_create
  USE kinds,                           ONLY: dp,&
                                             default_string_length
  USE mt_util,                         ONLY: MTin_create_screen_fn
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: REALDATA1D,&
                                             RECIPROCALSPACE,&
                                             COMPLEXDATA1D
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE input_cp2k,                      ONLY: parsed_cp2k_input
  USE simulation_cell,                 ONLY: cell_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_grids,                        ONLY: pw_grid_construct,&
                                             pw_grid_find_bounds,&
                                             pw_grid_setup

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: pw_green_fn_rebuild,&
            pw_green_setup,&
            setup_MT_ref_grid
  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='green_methods'

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
  SUBROUTINE pw_green_fn_rebuild ( green, pw_pool, p3m, alpha, cell, para_env, error )
    IMPLICIT NONE
    ! Arguments
    TYPE(greens_function_type), pointer      :: green
    TYPE(pw_pool_type), POINTER              :: pw_pool
    INTEGER, INTENT(IN), OPTIONAL            :: p3m
    REAL(KIND=dp), INTENT(IN), OPTIONAL          :: alpha
    TYPE(cell_type), POINTER, OPTIONAL       :: cell
    TYPE(cp_para_env_type), POINTER, OPTIONAL:: para_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = "pw_green_fn_rebuild", &
      routineP = moduleN//':'//routineN
    LOGICAL :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       
       CALL pw_pool_give_back_pw(pw_pool,green % p3m_charge, accept_non_compatible=.TRUE.)
       CALL pw_pool_give_back_pw(pw_pool,green % influence_function, accept_non_compatible=.TRUE.)
       CALL pw_pool_give_back_pw(pw_pool,green % screen_function, accept_non_compatible=.TRUE.)
       SELECT CASE ( green % method )
       CASE ( PERIODIC3D, ANALYTIC2D, ANALYTIC1D, ANALYTIC0D, MT2D, MT1D, MT0D )
          CALL pw_pool_create_pw ( pw_pool, green % influence_function,&
               use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE )
          IF ( PRESENT ( p3m ) ) THEN
             green % p3m = .TRUE.
             green % p3m_order = p3m
             IF ( PRESENT ( alpha ) ) THEN
                green % p3m_alpha = alpha
             ELSE
                CALL stop_program ( "pw_green_fn_rebuild", "P3M needs alpha PARAMETER")
             END IF
             CALL pw_pool_create_pw ( pw_pool, green % p3m_charge, use_data=REALDATA1D, &
                  in_space=RECIPROCALSPACE)
             CALL influence_factor ( green )
          ELSE
             green % p3m = .FALSE.
          END IF
          !
          SELECT CASE( green % method )
          CASE(MT0D,MT1D,MT2D)
             CPPrecondition(PRESENT(cell),cp_failure_level,routineP,error,failure)
             CPPrecondition(PRESENT(para_env),cp_failure_level,routineP,error,failure)
             CALL MTin_create_screen_fn(green%screen_function,pw_pool=pw_pool,method=green%method,&
                  alpha=green%MT_alpha, rel_cutoff=green%MT_rel_cutoff, cell=cell, &
                  special_dimension=green%special_dimension, slab_size=green%slab_size, &
                  para_env=para_env, super_ref_pw_grid=green%super_ref_pw_grid, error=error)
          END SELECT
          !
       CASE ( HOCKNEY2D, HOCKNEY1D, HOCKNEY0D )
          CALL stop_program ( "pw_green_fn_rebuild", " Method not implemented ")
       CASE DEFAULT
          CALL stop_program ( "pw_green_fn_rebuild", " Method not implemented ")
       END SELECT
    END IF
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
  SUBROUTINE pw_green_setup ( green, poisson_section, abc, error )
    IMPLICIT NONE
    ! Arguments
    TYPE(greens_function_type), POINTER       :: green
    TYPE(section_vals_type), POINTER          :: poisson_section
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)   :: abc
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = "pw_green_setup", &
      routineP = moduleN//':'//routineN
    TYPE(section_vals_type), POINTER          :: mt_section
    LOGICAL :: failure
    INTEGER :: dim, i, my_val, my_per
    INTEGER, DIMENSION(3) :: perd

    failure = .FALSE.
    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.(ASSOCIATED(green)),cp_failure_level,routineP,error,failure)
       CALL pw_green_create(green, error=error)
       CALL section_vals_val_get(poisson_section,"PERIODIC",i_val=my_per,error=error)
       SELECT CASE(my_per)
       CASE(use_perd_x)
          perd = (/1,0,0/)
       CASE(use_perd_y)
          perd = (/0,1,0/)
       CASE(use_perd_z)
          perd = (/0,0,1/)
       CASE(use_perd_xy)
          perd = (/1,1,0/)
       CASE(use_perd_xz)
          perd = (/1,0,1/)
       CASE(use_perd_yz)
          perd = (/0,1,1/)
       CASE(use_perd_xyz)
          perd = (/1,1,1/)
       CASE(use_perd_none)
          perd = (/0,0,0/)
       CASE DEFAULT
          CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT
       dim = COUNT(perd == 1)
       CALL section_vals_val_get(poisson_section,"POISSON_SOLVER",i_val=my_val,error=error)
       SELECT CASE (my_val)
       CASE (use_periodic)
          green%method = PERIODIC3D
          IF (dim /= 3) THEN
             CALL stop_program(routineP,&
                  "Illegal combination of periodicity and Poisson solver (perd)")
          END IF
       CASE (use_analytic)
          SELECT CASE (dim)
          CASE (0)
             green%method = ANALYTIC0D
             green%radius = 0.5_dp*MINVAL(abc)
          CASE (1)
             green%method = ANALYTIC1D
             green%special_dimension = MAXLOC(perd,1)
             green%radius = MAXVAL(abc)
             DO i=1,3
                IF (i == green%special_dimension) CYCLE
                green%radius = MIN(green%radius,0.5_dp*abc(i))
             END DO
          CASE (2)
             green%method = ANALYTIC2D
             i = MINLOC(perd,1)
             green%special_dimension = i
             green%slab_size = abc(i)
          CASE (3)
             green%method = PERIODIC3D
          CASE DEFAULT
             CALL stop_program(routineP,"Illegal value for perd found")
          END SELECT
       CASE (use_hockney)
          SELECT CASE (dim )
          CASE (0)
             green%method = HOCKNEY0D
          CASE (1)
             green%method = HOCKNEY1D
          CASE (2)
             green%method = HOCKNEY2D
          CASE (3)
             CALL stop_program(routineP,&
                  "Illegal combination of periodicity and Poisson solver (HOCKNEY)")
          CASE DEFAULT
             CALL stop_program(routineP,"Illegal value for perd found")
          END SELECT
       CASE (use_mt)
          mt_section => section_vals_get_subs_vals(poisson_section,"MT",error=error)
          CALL section_vals_val_get(mt_section,"REL_CUTOFF",r_val=green%MT_rel_cutoff,error=error)
          CALL section_vals_val_get(mt_section,"ALPHA",r_val=green%MT_alpha,error=error)
          green%MT_alpha=green%MT_alpha/MINVAL(abc)
          SELECT CASE (dim)
          CASE (0)
             green%method = MT0D
             green%radius = 0.5_dp*MINVAL(abc)
          CASE (1)
             green%method = MT1D
             green%special_dimension = MAXLOC(perd,1)
             green%radius = MAXVAL(abc)
             DO i=1,3
                IF (i == green%special_dimension) CYCLE
                green%radius = MIN(green%radius,0.5_dp*abc(i))
             END DO
          CASE (2)
             green%method = MT2D
             i = MINLOC(perd,1)
             green%special_dimension = i
             green%slab_size = abc(i)
          CASE (3)
             CALL stop_program(routineP,"Illegal combination of periodicity and Poisson solver (MT)")
          CASE DEFAULT
             CALL stop_program(routineP,"Illegal value for perd found")
          END SELECT
       CASE DEFAULT
          CALL stop_program(routineP,"An unknown Poisson solver was specified")
       END SELECT

    END IF
  END SUBROUTINE pw_green_setup

!!****** greens_fn/setup_MT_ref_grid [1.0] *
!!
!!   NAME
!!     setup_MT_ref_grid
!!
!!   FUNCTION
!!     Initialize the super-reference grid if Tuckerman is present
!!
!!   AUTHOR
!!     03-2005 Teodoro Laino [teo]
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE setup_MT_ref_grid(green, cutilev, grid_span, cell_ref, para_env, &
          old_pw_grid, input, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(greens_function_type), POINTER       :: green
    TYPE(section_vals_type), POINTER          :: input
    REAL(KIND=dp), INTENT(IN) :: cutilev
    INTEGER, INTENT(IN)       :: grid_span
    TYPE(cell_type), pointer  :: cell_ref
    TYPE(cp_para_env_type), POINTER :: para_env
    TYPE(pw_grid_type), POINTER     :: old_pw_grid
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error    
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = "setup_MT_ref_grid", &
      routineP = moduleN//':'//routineN
    TYPE(section_vals_type), POINTER          :: poisson_section
    INTEGER :: my_val, stat
    INTEGER, DIMENSION(2, 3) :: bounds1, bounds2
    REAL(KIND=dp) :: my_cutilev
    LOGICAL :: failure

    failure = .FALSE.
    NULLIFY(poisson_section)
    IF (.NOT. failure) THEN
       poisson_section => section_vals_get_subs_vals(input,"DFT%POISSON_DFT",error=error)    
       CALL section_vals_val_get(poisson_section,"POISSON_SOLVER",i_val=my_val,error=error)
       SELECT CASE(my_val)
       CASE(use_mt)
          !
          ! Check if grids will be the same... In this case we don't use a super-reference grid
          !
          my_cutilev = cutilev * green%MT_rel_cutoff
          CALL pw_grid_find_bounds ( bounds1, cell_ref % h_inv, &
               ABS(cutilev), fft=.TRUE. )
          CALL pw_grid_find_bounds ( bounds2, cell_ref % h_inv, &
               ABS(my_cutilev), fft=.TRUE. )
          IF (ANY (bounds1 /= bounds2)) THEN 
             ALLOCATE(green%super_ref_pw_grid,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             CALL pw_grid_construct(green%super_ref_pw_grid)
             green%super_ref_pw_grid%para%rs_dims ( 1 ) = para_env%num_pe
             green%super_ref_pw_grid%para%rs_dims ( 2 ) = 1
             green%super_ref_pw_grid%grid_span = grid_span
             CALL pw_grid_setup(cell_ref, green%super_ref_pw_grid, my_cutilev,&
                  pe_group=para_env%group, blocked=.FALSE.)
             old_pw_grid => green%super_ref_pw_grid
          END IF
       CASE DEFAULT
          ! Do nothing 
       END SELECT
    END IF
END SUBROUTINE setup_MT_ref_grid

END MODULE green_methods
