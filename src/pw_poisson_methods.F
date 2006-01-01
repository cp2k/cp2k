!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****h* cp2k/pw_poisson_methods *
!!
!!   NAME
!!     pw_poisson_methods
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     09.2005 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE pw_poisson_methods
  USE kinds,                           ONLY: dp
  USE pw_poisson_types,                ONLY: pw_poisson_type, greens_fn_type,&
                                             pw_green_create,pw_green_release
  USE timings,                         ONLY: timeset, timestop
  USE input_constants,                 ONLY: ANALYTIC0D,&
                                             ANALYTIC1D,&
                                             ANALYTIC2D,&
                                             MT0D,&
                                             MT1D,&
                                             MT2D,&
                                             PERIODIC3D,&
                                             MULTIPOLE0D,&
                                             use_periodic,&
                                             use_analytic,&
                                             use_mt,&
                                             use_none,&
                                             use_multipole,&
                                             use_perd_x,&
                                             use_perd_y,&
                                             use_perd_z,&
                                             use_perd_xy,&
                                             use_perd_xz,&
                                             use_perd_yz,&
                                             use_perd_xyz,&
                                             use_perd_none,&
                                             do_ewald_none,&
                                             do_ewald_spme
  USE fft_tools,                       ONLY: FFT_RADIX_NEXT,&
                                             FFT_RADIX_NEXT_ODD,&
                                             fft_radix_operations
  USE mt_util,                         ONLY: MTin_create_screen_fn
  USE pw_types,                        ONLY: REALDATA1D,&
                                             RECIPROCALSPACE,&
                                             COMPLEXDATA1D,&
                                             pw_type,&
                                             pw_derive,&
                                             pw_transfer,&
                                             pw_p_type,&
                                             pw_release,&
                                             pw_integral_ab,&
                                             pw_copy
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_release,&
                                             section_vals_retain
  USE simulation_cell,                 ONLY: cell_type, cell_release, cell_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE pw_grid_types,                   ONLY: pw_grid_type,&
                                             HALFSPACE
  USE pw_grids,                        ONLY: pw_grid_create,&
                                             pw_grid_find_bounds,&
                                             pw_grid_setup,&
                                             pw_grid_release,&
                                             pw_grid_retain 
  USE mathconstants,                   ONLY: fourpi,twopi
  USE pw_pool_types,                   ONLY: pw_pool_type,&
                                             pw_pool_create,&
                                             pw_pool_p_type,&
                                             pw_pools_dealloc,&
                                             pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pools_copy
                                             
#include "cp_common_uses.h"

  IMPLICIT NONE
  private

  logical, private, parameter :: debug_this_module=.true.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='pw_poisson_methods'

PUBLIC :: pw_poisson_cleanup, pw_poisson_rebuild, pw_poisson_check, &
     pw_poisson_solve, pw_poisson_set
  !!***
  !****************************************************************************
contains

  !!****** pw_poisson_methods/pw_poisson_cleanup [1.0] *
  !!
  !!   NAME
  !!     pw_poisson_cleanup
  !!
  !!   FUNCTION
  !!     removes all the object created from the parameters pw_pools and cell
  !!     and used to solve the poisson equation like the green function and
  !!     all the things allocated in pw_poisson_rebuild
  !!
  !!   MODIFICATION HISTORY
  !!     none
  !!
  !!*** *************************************************************************
  SUBROUTINE pw_poisson_cleanup (poisson_env, error)
    TYPE(pw_poisson_type), pointer :: poisson_env
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

    LOGICAL :: can_give_back, failure
    TYPE(pw_pool_type), POINTER :: pw_pool
    TYPE(greens_fn_type), pointer :: gftype
    CHARACTER(len=*), PARAMETER :: routineN="pw_poisson_cleanup",&
         routineP=moduleN//':'//routineN

    failure=.false.
    CPPrecondition(ASSOCIATED(poisson_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(poisson_env%ref_count>0,cp_failure_level,routineP,error,failure)

    IF (.NOT.failure) THEN
       nullify(pw_pool)
       IF (ASSOCIATED(poisson_env%pw_pools)) THEN
          pw_pool => poisson_env%pw_pools(poisson_env%pw_level)%pool
       END IF
       CALL pw_green_release(poisson_env%green_fft,pw_pool=pw_pool,error=error)
       poisson_env%rebuild=.true.
    END IF
  END SUBROUTINE pw_poisson_cleanup
  !******************************************************************************

  !!****f* pw_poisson_methods/pw_poisson_check *
  !!
  !!   NAME
  !!     pw_poisson_check
  !!
  !!   FUNCTION
  !!     checks if pw_poisson_rebuild has to be called and calls it if needed
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     - poisson_env: the object to be checked
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** **********************************************************************
  SUBROUTINE pw_poisson_check(poisson_env,error)
    TYPE(pw_poisson_type), pointer :: poisson_env
    type(cp_error_type), optional, intent(inout) :: error

    LOGICAL :: failure, rebuild
    character(len=*), parameter :: routineN='pw_poisson_check',&
         routineP=moduleN//':'//routineN
    INTEGER :: i_val, ewald_kind,s_order
    TYPE(section_vals_type), pointer :: ewald_section
    REAL(kind=dp) :: alpha
    TYPE(greens_fn_type), pointer :: green

    failure=.false.

    nullify(ewald_section)
    CPPrecondition(ASSOCIATED(poisson_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(poisson_env%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPrecondition(ASSOCIATED(poisson_env%pw_pools),cp_failure_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(poisson_env%cell),cp_failure_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(poisson_env%parameters),cp_failure_level,routineP,error,failure)
       CPPrecondition(poisson_env%pw_level>=LBOUND(poisson_env%pw_pools,1),cp_failure_level,routineP,error,failure)
       CPPrecondition(poisson_env%pw_level<=UBOUND(poisson_env%pw_pools,1),cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT.failure) THEN
       green => poisson_env%green_fft
       rebuild=poisson_env%rebuild
       call section_vals_val_get(poisson_env%parameters,"POISSON_SOLVER",&
            i_val=i_val,error=error)
       ewald_section => section_vals_get_subs_vals(poisson_env%parameters,&
            "EWALD", can_return_null=.true., error=error)
       rebuild=rebuild.OR.(poisson_env%method/=i_val)&
            .OR..NOT.ASSOCIATED(green)
       poisson_env%method=i_val
       IF (.NOT.rebuild) THEN
          IF (ASSOCIATED(ewald_section)) THEN
             CALL section_vals_val_get(ewald_section,"ewald_type",i_val=ewald_kind,error=error)
             IF (ewald_kind==do_ewald_spme) THEN
                CALL section_vals_val_get(ewald_section,"alpha",r_val=alpha,error=error)
                rebuild=(alpha/=green%p3m_alpha).OR.rebuild
                CALL section_vals_val_get(ewald_section,"o_spline",i_val=s_order,error=error)
                rebuild=(s_order/=green%p3m_order).OR.rebuild
             END IF
          END IF
          SELECT CASE(poisson_env%method)
          CASE(use_analytic)
             SELECT CASE(green%method)
             CASE(ANALYTIC0D,ANALYTIC1D,ANALYTIC2D,PERIODIC3D)
             CASE default
                rebuild=.TRUE.
             END SELECT
          CASE(use_mt)
             SELECT CASE(green%method)
             CASE(MT0D,MT1D,MT2D)
             CASE default
                rebuild=.TRUE.
             END SELECT
             CALL section_vals_val_get(poisson_env%parameters,&
                  "MT%ALPHA",r_val=alpha,error=error)
             rebuild=(alpha/=green%mt_alpha).OR.rebuild
          CASE default
             CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
          END SELECT
       END IF
       IF (rebuild) THEN
          poisson_env%rebuild=.true.
          CALL pw_poisson_cleanup(poisson_env,error=error)
       END IF
    END IF
  END SUBROUTINE pw_poisson_check
  !***************************************************************************

  !!****f* pw_poisson_methods/pw_poisson_rebuild *
  !!
  !!   NAME
  !!     pw_poisson_rebuild
  !!
  !!   FUNCTION
  !!     rebuilds all the internal values needed to use the poisson solver
  !!
  !!   NOTES
  !!     rebuilds if poisson_env%rebuild is true
  !!
  !!   ARGUMENTS
  !!     - poisson_env: the environment to rebuild
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** **********************************************************************
  SUBROUTINE pw_poisson_rebuild(poisson_env,error)
    TYPE(pw_poisson_type), pointer :: poisson_env
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='pw_poisson_rebuild',&
         routineP=moduleN//':'//routineN

    failure=.false.

    CPPrecondition(ASSOCIATED(poisson_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(poisson_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(poisson_env%cell),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(poisson_env%pw_pools),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(poisson_env%parameters),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       IF (poisson_env%rebuild) THEN
          CALL pw_poisson_cleanup(poisson_env,error=error)

          CALL section_vals_val_get(poisson_env%parameters,"POISSON_SOLVER",&
               i_val=poisson_env%method,error=error)
          SELECT CASE (poisson_env%method)
          CASE(use_periodic,use_analytic,use_mt,&
               use_multipole)
             CALL pw_green_create(poisson_env%green_fft,cell=poisson_env%cell,&
                  pw_pool=poisson_env%pw_pools(poisson_env%pw_level)%pool,&
                  poisson_section=poisson_env%parameters,&
                  mt_super_ref_pw_grid=poisson_env%mt_super_ref_pw_grid,&
                  error=error)
          CASE(use_none)
          CASE default
             CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
          END SELECT
          poisson_env%rebuild=.false.
       END IF
    END IF
  END SUBROUTINE pw_poisson_rebuild
  !***************************************************************************

  !!****** pw_poisson_methods/pw_poisson_solve [1.0] *
  !!
  !!   NAME
  !!     pw_poisson_solve
  !!
  !!   FUNCTION
  !!     Solve Poisson equation in a plane wave basis set
  !!     Obtains electrostatic potential and its derivatives with respect to r
  !!     from the density
  !!
  !!   AUTHOR
  !!     apsi
  !!
  !!   MODIFICATION HISTORY
  !!     JGH (13-Mar-2001) : completely revised
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pw_poisson_solve ( poisson_env, density, ehartree, vhartree,&
       dvhartree, h_stress, error )

    TYPE(pw_poisson_type), POINTER :: poisson_env
    TYPE(pw_type), pointer             :: density
    REAL(kind=dp), INTENT(out), optional :: ehartree
    TYPE(pw_type), pointer, &
         OPTIONAL                       :: vhartree
    TYPE(pw_p_type), DIMENSION(3), &
         OPTIONAL                :: dvhartree
    REAL(KIND=dp), DIMENSION(3, 3), &
         INTENT(OUT), OPTIONAL                  :: h_stress
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error

    INTEGER                                  :: alpha, beta, i, ng
    INTEGER, DIMENSION(3)                    :: n
    REAL(KIND=dp)                                :: ffa
    TYPE(pw_grid_type), POINTER              :: pw_grid
    TYPE(pw_p_type)                            :: dvg( 3 )
    TYPE(pw_type), POINTER :: rhog,tmpg
    TYPE(pw_pool_type), pointer :: pw_pool
    logical :: failure
    CHARACTER(len=*), PARAMETER :: routineN="pw_poisson_solve",&
         routineP=moduleN//':'//routineN
    !------------------------------------------------------------------------------

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(poisson_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(poisson_env%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CALL pw_poisson_rebuild(poisson_env,error=error)
       CALL cp_error_check(error,failure)
    END IF
    IF (.NOT.failure) THEN
       NULLIFY (pw_grid,rhog)
       ! point pw
       pw_pool => poisson_env%pw_pools(poisson_env%pw_level)%pool
       pw_grid => pw_pool % pw_grid
       IF ( PRESENT ( vhartree ) ) THEN
          CPPrecondition(ASSOCIATED(vhartree),cp_failure_level,routineP,error,failure)
          CALL cp_assert(pw_pool%pw_grid%id_nr==vhartree%pw_grid%id_nr,cp_assertion_failed,&
               cp_failure_level,routineP,&
               "vhartree has a different grid than the poisson solver",error,failure)
       END IF
       ! density in G space
       CALL pw_pool_create_pw ( pw_pool,rhog, use_data=COMPLEXDATA1D,in_space = RECIPROCALSPACE,&
            error=error)

       ! apply the greens function
       ng = SIZE ( pw_grid % gsq )
       SELECT CASE ( poisson_env%green_fft % method )
       CASE ( PERIODIC3D, ANALYTIC2D, ANALYTIC1D, ANALYTIC0D, MT2D, MT1D, MT0D, MULTIPOLE0D )
          CALL pw_transfer ( density, rhog )
          IF (PRESENT(ehartree).AND..NOT.PRESENT(vhartree)) THEN
             CALL pw_pool_create_pw ( pw_pool,tmpg, use_data=COMPLEXDATA1D,&
                  in_space = RECIPROCALSPACE, error=error)
             CALL pw_copy(rhog,tmpg)
          END IF

          rhog % cc ( 1 : ng )  = rhog % cc ( 1 : ng ) *  &
               poisson_env%green_fft % influence_fn % cc  ( 1 : ng )

          IF ( PRESENT ( vhartree ) ) THEN
             CALL pw_transfer ( rhog, vhartree )
             IF ( PRESENT ( ehartree ) ) THEN
                ehartree = 0.5_dp * pw_integral_ab ( density , vhartree )
             END IF
          ELSE if (present(ehartree)) then
             ehartree=0.5_dp*pw_integral_ab(rhog,tmpg)
             CALL pw_pool_give_back_pw(pw_pool,tmpg,error=error)
          END IF
       CASE DEFAULT
          CALL cp_unimplemented_error(routineP,"unknown poisson method"//&
               cp_to_string(poisson_env%green_fft%method),error)
       END SELECT

       ! do we need to calculate the derivative of the potential?
       IF ( PRESENT ( h_stress ) .OR. PRESENT ( dvhartree ) ) THEN
          IF ( poisson_env%green_fft % method /= PERIODIC3D ) THEN
             CALL cp_unimplemented_error(routineP,&
                  "Pressure calculation works ONLY for fully priodic systems",error)
          END IF
          DO i = 1, 3
             nullify(dvg(i)%pw)
             CALL pw_pool_create_pw (pw_pool, dvg ( i )%pw, use_data=COMPLEXDATA1D,&
                  in_space= RECIPROCALSPACE,error=error)
             n = 0
             n ( i ) = 1
             CALL pw_copy ( rhog, dvg ( i )%pw )
             CALL pw_derive ( dvg ( i )%pw, n )
          END DO
          ! save the derivatives
          IF ( PRESENT ( dvhartree ) ) THEN
             DO i = 1, 3
                CALL pw_transfer ( dvg ( i )%pw, dvhartree ( i ) % pw )
             END DO
          END IF
          ! calculate the contribution to the stress tensor
          ! this is only the contribution from the Greens function
          ! and the volume factor of the plane waves
          IF ( PRESENT ( h_stress ) ) THEN

             ffa = -1.0_dp / fourpi
             h_stress = 0.0_dp
             DO alpha = 1, 3
                h_stress ( alpha, alpha ) = ehartree
                DO beta = alpha, 3
                   h_stress ( alpha, beta ) = h_stress ( alpha, beta ) &
                        + ffa * pw_integral_ab ( dvg ( alpha )%pw, dvg ( beta )%pw )
                   h_stress ( beta, alpha ) = h_stress ( alpha, beta )
                END DO
             END DO

          END IF
          DO i = 1, 3
             CALL pw_pool_give_back_pw ( pw_pool,&
                  dvg ( i )%pw, error=error )
          END DO
       END IF

       CALL pw_pool_give_back_pw (pw_pool,&
            rhog, error=error )
    END IF
  END SUBROUTINE pw_poisson_solve
  !****************************************************************************

  !!****** pw_poisson_methods/pw_poisson_set [1.0] *
  !!
  !!   NAME
  !!     pw_poisson_set
  !!
  !!   FUNCTION
  !!     sets cell, grids and parameters used by the poisson solver
  !!     You should call this at least once (and set everything)
  !!     before using the poisson solver.
  !!     Smart, doesn't set the thing twice to the same value
  !!     Keeps track of the need to rebuild the poisson_env
  !!
  !!   NOTES
  !!     Checks everything at the end. This means that after *each* call to
  !!     this method the poisson env must be fully ready, so the first time
  !!     you have to set everything at once. Change this behaviour?
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** *************************************************************************
  SUBROUTINE pw_poisson_set ( poisson_env, cell, parameters, pw_pools,&
       use_level, mt_super_ref_pw_grid, force_rebuild, error )

    TYPE(pw_poisson_type), pointer :: poisson_env
    TYPE(cell_type), POINTER, optional :: cell
    TYPE(section_vals_type), POINTER, optional :: parameters
    TYPE(pw_pool_p_type), DIMENSION(:), POINTER, optional :: pw_pools
    INTEGER, INTENT(in), optional :: use_level
    TYPE(pw_grid_type), POINTER, OPTIONAL :: mt_super_ref_pw_grid
    LOGICAL, INTENT(in), optional :: force_rebuild
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_pw_poisson_solver', &
         routineP = moduleN//':'//routineN
    TYPE(pw_pool_p_type), DIMENSION(:), POINTER :: tmp_pools
    LOGICAL                                  :: failure, same
    integer :: i

    IF (PRESENT(parameters)) THEN
       CPPrecondition(ASSOCIATED(parameters),cp_failure_level,routineP,error,failure)
       CALL section_vals_retain(parameters,error=error)
       CALL section_vals_release(poisson_env%parameters,error=error)
       poisson_env%parameters => parameters
    END IF
    IF (PRESENT(cell)) THEN
       CPPrecondition(ASSOCIATED(cell),cp_failure_level,routineP,error,failure)
       IF (ASSOCIATED(poisson_env%cell)) THEN
          IF (poisson_env%cell%id_nr/=cell%id_nr) THEN
             CALL pw_poisson_cleanup(poisson_env,error=error)
          END IF
       END IF
       CALL cell_retain(cell,error=error)
       CALL cell_release(poisson_env%cell,error=error)
       poisson_env%cell => cell
       poisson_env%rebuild=.true.
    END IF
    IF (PRESENT(pw_pools)) THEN
       CPPrecondition(ASSOCIATED(pw_pools),cp_failure_level,routineP,error,failure)
       same=.false.
       IF (ASSOCIATED(poisson_env%pw_pools)) THEN
          same=SIZE(poisson_env%pw_pools)==size(pw_pools)
          IF (same) THEN
             DO i =1,SIZE(pw_pools)
                IF (poisson_env%pw_pools(i)%pool%id_nr/=&
                     pw_pools(i)%pool%id_nr) same=.false.
             END DO
          END IF
       END IF
       IF (.NOT.same) THEN
          poisson_env%rebuild=.TRUE.
          CALL pw_pools_copy(pw_pools,tmp_pools,error=error)
          CALL pw_pools_dealloc(poisson_env%pw_pools,error=error)
          poisson_env%pw_pools => tmp_pools
       END IF
    END IF
    IF (PRESENT(use_level)) poisson_env%pw_level=use_level
    IF (PRESENT(mt_super_ref_pw_grid)) THEN
       IF (ASSOCIATED(mt_super_ref_pw_grid)) THEN
          CALL pw_grid_retain(mt_super_ref_pw_grid,error=error)
       END IF
       CALL pw_grid_release(poisson_env%mt_super_ref_pw_grid,error=error)
       poisson_env%mt_super_ref_pw_grid => mt_super_ref_pw_grid
    END IF
    IF (PRESENT(force_rebuild)) THEN
       IF (force_rebuild) poisson_env%rebuild=.true.
    END IF
    CALL pw_poisson_check(poisson_env,error=error)
  END SUBROUTINE pw_poisson_set
  !******************************************************************************

END MODULE pw_poisson_methods
