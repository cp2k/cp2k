!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      09.2005 created [fawzi]
!> \author fawzi
! *****************************************************************************
MODULE pw_poisson_methods
  USE f77_blas
  USE input_constants,                 ONLY: &
       ANALYTIC0D, ANALYTIC1D, ANALYTIC2D, MT0D, MT1D, MT2D, MULTIPOLE0D, &
       PERIODIC3D, WAVELET0D, WAVELET1D, WAVELET2D, WAVELET3D, do_ewald_spme, &
       use_analytic, use_gs_grid, use_mt, use_multipole, use_none, &
       use_periodic, use_rs_grid, use_wavelet
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_release,&
                                             section_vals_retain,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: fourpi
  USE ps_wavelet_types,                ONLY: cp2k_distribution_to_z_slices,&
                                             ps_wavelet_create,&
                                             ps_wavelet_solve,&
                                             ps_wavelet_type,&
                                             z_slices_to_cp2k_distribution
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_grids,                        ONLY: pw_grid_compare,&
                                             pw_grid_release,&
                                             pw_grid_retain
  USE pw_methods,                      ONLY: pw_copy,&
                                             pw_derive,&
                                             pw_integral_ab,&
                                             pw_transfer
  USE pw_poisson_types,                ONLY: greens_fn_type,&
                                             pw_green_create,&
                                             pw_green_release,&
                                             pw_poisson_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_p_type,&
                                             pw_pool_type,&
                                             pw_pools_copy,&
                                             pw_pools_dealloc
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type,&
                                             pw_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pw_poisson_methods'

  PUBLIC :: pw_poisson_cleanup, pw_poisson_rebuild, pw_poisson_check, &
       pw_poisson_solve, pw_poisson_set
CONTAINS

! *****************************************************************************
!> \brief removes all the object created from the parameters pw_pools and cell
!>      and used to solve the poisson equation like the green function and
!>      all the things allocated in pw_poisson_rebuild
!> \par History
!>      none
! *****************************************************************************
  SUBROUTINE pw_poisson_cleanup (poisson_env, error)
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_poisson_cleanup', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(pw_pool_type), POINTER              :: pw_pool

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(poisson_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(poisson_env%ref_count>0,cp_failure_level,routineP,error,failure)

    IF (.NOT.failure) THEN
       NULLIFY(pw_pool)
       IF (ASSOCIATED(poisson_env%pw_pools)) THEN
          pw_pool => poisson_env%pw_pools(poisson_env%pw_level)%pool
       END IF
       CALL pw_green_release(poisson_env%green_fft,pw_pool=pw_pool,error=error)
       poisson_env%rebuild=.TRUE.

    END IF
  END SUBROUTINE pw_poisson_cleanup

! *****************************************************************************
!> \brief checks if pw_poisson_rebuild has to be called and calls it if needed
!> \param poisson_env the object to be checked
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE pw_poisson_check(poisson_env,error)
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_poisson_check', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ewald_kind, i_val, itype_scf, &
                                                s_order
    LOGICAL                                  :: failure, rebuild
    REAL(kind=dp)                            :: alpha
    TYPE(greens_fn_type), POINTER            :: green
    TYPE(ps_wavelet_type), POINTER           :: wavelet
    TYPE(section_vals_type), POINTER         :: ewald_section

    failure=.FALSE.

    NULLIFY(ewald_section)
    CPPrecondition(ASSOCIATED(poisson_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(poisson_env%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPrecondition(ASSOCIATED(poisson_env%pw_pools),cp_failure_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(poisson_env%parameters),cp_failure_level,routineP,error,failure)
       CPPrecondition(poisson_env%pw_level>=LBOUND(poisson_env%pw_pools,1),cp_failure_level,routineP,error,failure)
       CPPrecondition(poisson_env%pw_level<=UBOUND(poisson_env%pw_pools,1),cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT.failure) THEN
       green => poisson_env%green_fft
       wavelet => poisson_env%wavelet
       rebuild=poisson_env%rebuild
       CALL section_vals_val_get(poisson_env%parameters,"POISSON_SOLVER",&
            i_val=i_val,error=error)
       ewald_section => section_vals_get_subs_vals(poisson_env%parameters,&
            "EWALD", can_return_null=.TRUE., error=error)
       rebuild=rebuild.OR.(poisson_env%method/=i_val)&
            .OR..NOT.ASSOCIATED(green)
       poisson_env%method=i_val
       IF(i_val==use_wavelet)THEN
         poisson_env%used_grid=use_rs_grid
       ELSE
         poisson_env%used_grid=use_gs_grid
       END IF
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
             CALL section_vals_val_get(poisson_env%parameters,"MT%ALPHA",r_val=alpha,error=error)
             rebuild=(alpha/=green%mt_alpha).OR.rebuild
          CASE(use_wavelet)
             CALL section_vals_val_get(poisson_env%parameters,"WAVELET%SCF_TYPE",i_val=itype_scf,error=error)
             rebuild=(itype_scf/=wavelet%itype_scf).OR.rebuild
          CASE default
             CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
          END SELECT
       END IF
       IF (rebuild) THEN
          poisson_env%rebuild=.TRUE.
          CALL pw_poisson_cleanup(poisson_env,error=error)
       END IF
    END IF
  END SUBROUTINE pw_poisson_check

! *****************************************************************************
!> \brief rebuilds all the internal values needed to use the poisson solver
!> \param poisson_env the environment to rebuild
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      rebuilds if poisson_env%rebuild is true
!> \author fawzi
! *****************************************************************************
  SUBROUTINE pw_poisson_rebuild(poisson_env,density,error)
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(pw_type), OPTIONAL, POINTER         :: density
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_poisson_rebuild', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(poisson_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(poisson_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(poisson_env%pw_pools),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(poisson_env%parameters),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       IF (poisson_env%rebuild) THEN
          CALL pw_poisson_cleanup(poisson_env,error=error)

          CALL section_vals_val_get(poisson_env%parameters,"POISSON_SOLVER",&
               i_val=poisson_env%method,error=error)
           SELECT CASE (poisson_env%method)
           CASE(use_periodic,use_analytic,use_mt,use_multipole)
              CALL pw_green_create(poisson_env%green_fft,cell_hmat=poisson_env%cell_hmat,&
                   pw_pool=poisson_env%pw_pools(poisson_env%pw_level)%pool,&
                   poisson_section=poisson_env%parameters,&
                   mt_super_ref_pw_grid=poisson_env%mt_super_ref_pw_grid,&
                   error=error)
           CASE(use_wavelet)
              CPPrecondition(ASSOCIATED(density%pw_grid),cp_failure_level,routineP,error,failure)
              CALL ps_wavelet_create(poisson_env%parameters,poisson_env%wavelet,&
                                     density%pw_grid,error)
           CASE(use_none)
           CASE default
              CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
           END SELECT
          poisson_env%rebuild=.FALSE.
       END IF
    END IF
  END SUBROUTINE pw_poisson_rebuild

! *****************************************************************************
!> \brief Solve Poisson equation in a plane wave basis set
!>      Obtains electrostatic potential and its derivatives with respect to r
!>      from the density
!> \par History
!>      JGH (13-Mar-2001) : completely revised
!> \author apsi
! *****************************************************************************
  SUBROUTINE pw_poisson_solve ( poisson_env, density, ehartree, vhartree,&
       dvhartree, h_stress, error )

    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(pw_type), POINTER                   :: density
    REAL(kind=dp), INTENT(out), OPTIONAL     :: ehartree
    TYPE(pw_type), OPTIONAL, POINTER         :: vhartree
    TYPE(pw_p_type), DIMENSION(3), OPTIONAL  :: dvhartree
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(OUT), OPTIONAL                  :: h_stress
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_poisson_solve', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: alpha, beta, handle, i, ig, ng
    INTEGER, DIMENSION(3)                    :: n
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: ffa
    TYPE(pw_grid_type), POINTER              :: pw_grid
    TYPE(pw_p_type)                          :: dvg( 3 )
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(pw_type), POINTER                   :: rhog, rhor, tmpg

    CALL timeset(routineN,handle)

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(poisson_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(poisson_env%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CALL pw_poisson_rebuild(poisson_env,density,error=error)
       CALL cp_error_check(error,failure)
    END IF
    IF (.NOT.failure) THEN
       NULLIFY (pw_grid,rhog)
       ! point pw
       pw_pool => poisson_env%pw_pools(poisson_env%pw_level)%pool
       pw_grid => pw_pool % pw_grid
       IF ( PRESENT ( vhartree ) ) THEN
          CPPrecondition(ASSOCIATED(vhartree),cp_failure_level,routineP,error,failure)
          CALL cp_assert(pw_grid_compare(pw_pool%pw_grid,vhartree%pw_grid),cp_assertion_failed,&
               cp_failure_level,routineP,&
               "vhartree has a different grid than the poisson solver",error,failure)
       END IF
       ! density in G space
       CALL pw_pool_create_pw ( pw_pool,rhog, use_data=COMPLEXDATA1D,in_space = RECIPROCALSPACE,&
            error=error)
       ! apply the greens function
       ng = SIZE ( pw_grid % gsq )
       SELECT CASE(poisson_env%used_grid)
       CASE(use_gs_grid)
          SELECT CASE ( poisson_env%green_fft%method )
          CASE (PERIODIC3D, ANALYTIC2D, ANALYTIC1D, ANALYTIC0D, MT2D, MT1D, MT0D, MULTIPOLE0D)
             CALL pw_transfer ( density, rhog , error=error)
             IF (PRESENT(ehartree).AND..NOT.PRESENT(vhartree)) THEN
                CALL pw_pool_create_pw ( pw_pool,tmpg, use_data=COMPLEXDATA1D,&
                     in_space = RECIPROCALSPACE, error=error)
                CALL pw_copy(rhog,tmpg, error=error)
             END IF

             DO ig=1,ng
                rhog%cc(ig) = rhog%cc(ig) *  poisson_env%green_fft%influence_fn%cc(ig)
             ENDDO

             IF ( PRESENT ( vhartree ) ) THEN
                CALL pw_transfer ( rhog, vhartree , error=error)
                IF ( PRESENT ( ehartree ) ) THEN
                   ehartree = 0.5_dp * pw_integral_ab ( density , vhartree , error=error)
                END IF
             ELSE IF (PRESENT(ehartree)) THEN
                ehartree=0.5_dp*pw_integral_ab(rhog,tmpg, error=error)
                CALL pw_pool_give_back_pw(pw_pool,tmpg,error=error)
             END IF
          CASE DEFAULT
             CALL cp_unimplemented_error(routineP,"unknown poisson method"//&
                  cp_to_string(poisson_env%green_fft%method),error)
          END SELECT
       CASE(use_rs_grid)
          CALL pw_transfer ( density, rhog , error=error)
          ! point pw
          pw_pool => poisson_env%pw_pools(poisson_env%pw_level)%pool
          pw_grid => pw_pool % pw_grid
          IF ( PRESENT ( vhartree ) ) THEN
             CPPrecondition(ASSOCIATED(vhartree),cp_failure_level,routineP,error,failure)
             CALL cp_assert(pw_grid_compare(pw_pool%pw_grid,vhartree%pw_grid),cp_assertion_failed,&
                  cp_failure_level,routineP,&
                  "vhartree has a different grid than the poisson solver",error,failure)
          END IF

          IF (PRESENT(ehartree).AND..NOT.PRESENT(vhartree)) THEN
             CALL pw_pool_create_pw ( pw_pool,tmpg, use_data=COMPLEXDATA1D,&
                  in_space = RECIPROCALSPACE, error=error)
             CALL pw_transfer(density,tmpg, error=error)
          END IF

          ! density in R space
          CALL pw_pool_create_pw ( pw_pool,rhor, use_data=REALDATA3D,in_space = REALSPACE,&
               error=error)
          CALL pw_transfer(density,rhor, error=error)
          CALL cp2k_distribution_to_z_slices(rhor ,poisson_env%wavelet, rhor%pw_grid, error)
          CALL ps_wavelet_solve(poisson_env%wavelet,rhor%pw_grid,ehartree,error)

          CALL z_slices_to_cp2k_distribution(rhor ,poisson_env% wavelet,rhor% pw_grid, error)
          IF ( PRESENT ( vhartree ) ) THEN
             CALL pw_transfer(rhor,vhartree, error=error)
             IF ( PRESENT ( ehartree ) ) THEN
                ehartree = 0.5_dp * pw_integral_ab ( density , vhartree , error=error)
             END IF
          ELSE IF (PRESENT(ehartree)) THEN
             ehartree=0.5_dp*pw_integral_ab(rhog,tmpg, error=error)
             CALL pw_pool_give_back_pw(pw_pool, tmpg, error=error)
          END IF
          CALL pw_pool_give_back_pw (pw_pool, rhor, error=error )
       END SELECT

       ! do we need to calculate the derivative of the potential?
       IF ( PRESENT ( h_stress ) .OR. PRESENT ( dvhartree ) ) THEN
          DO i = 1, 3
             NULLIFY(dvg(i)%pw)
             CALL pw_pool_create_pw (pw_pool, dvg ( i )%pw, use_data=COMPLEXDATA1D,&
                  in_space= RECIPROCALSPACE,error=error)
             n = 0
             n ( i ) = 1
             CALL pw_copy ( rhog, dvg ( i )%pw , error=error)
             CALL pw_derive ( dvg ( i )%pw, n , error=error)
          END DO
          ! save the derivatives
          IF ( PRESENT ( dvhartree ) ) THEN
             DO i = 1, 3
                CALL pw_transfer ( dvg ( i )%pw, dvhartree ( i ) % pw , error=error)
             END DO
          END IF
          ! Calculate the contribution to the stress tensor this is only the contribution from
          ! the Greens FUNCTION and the volume factor of the plane waves
          IF ( PRESENT ( h_stress ) ) THEN
             ffa = -1.0_dp / fourpi
             h_stress = 0.0_dp
             DO alpha = 1, 3
                h_stress ( alpha, alpha ) = ehartree
                DO beta = alpha, 3
                   h_stress ( alpha, beta ) = h_stress ( alpha, beta ) &
                        + ffa * pw_integral_ab ( dvg ( alpha )%pw, dvg ( beta )%pw , error=error)
                   h_stress ( beta, alpha ) = h_stress ( alpha, beta )
                END DO
             END DO
             ! Handle the periodicity cases for the Stress Tensor
             SELECT CASE(poisson_env%used_grid)
             CASE(use_gs_grid)
                ! FFT based Poisson-Solver
                SELECT CASE(poisson_env%green_fft%method)
                CASE(PERIODIC3D)
                   ! Do Nothing
                CASE(ANALYTIC2D, MT2D)
                   ! Zero the 1 non-periodic component
                   alpha = poisson_env%green_fft%special_dimension
                   h_stress(:,alpha) = 0.0_dp
                   h_stress(alpha,:) = 0.0_dp
                   CALL cp_unimplemented_error(routineP,"Stress Tensor not tested for 2D systems.",error)
                CASE(ANALYTIC1D, MT1D)
                   ! Zero the 2 non-periodic components
                   DO alpha = 1, 3
                      DO beta = alpha, 3
                         IF  ((alpha/=poisson_env%green_fft%special_dimension).OR.&
                              ( beta/=poisson_env%green_fft%special_dimension)) THEN
                            h_stress(alpha,beta) = 0.0_dp
                            h_stress(beta,alpha) = 0.0_dp
                         END IF
                      END DO
                   END DO
                   CALL cp_unimplemented_error(routineP,"Stress Tensor not tested for 1D systems.",error)
                CASE(ANALYTIC0D, MT0D, MULTIPOLE0D)
                   ! Zero the full stress tensor
                   h_stress = 0.0_dp
                CASE DEFAULT
                   CALL cp_unimplemented_error(routineP,"unknown poisson method"//&
                        cp_to_string(poisson_env%green_fft%method),error)
                END SELECT
             CASE(use_rs_grid)
                ! Wavelet based Poisson-Solver
                SELECT CASE(poisson_env%wavelet%method)
                CASE(WAVELET3D)
                   ! Do Nothing
                CASE(WAVELET2D)
                   ! Zero the 1 non-periodic component
                   alpha = poisson_env%wavelet%special_dimension
                   h_stress(:,alpha) = 0.0_dp
                   h_stress(alpha,:) = 0.0_dp
                   CALL cp_unimplemented_error(routineP,"Stress Tensor not tested for 2D systems.",error)
                CASE(WAVELET1D)
                   ! Zero the 2 non-periodic components
                   CALL cp_unimplemented_error(routineP,"WAVELET 1D not implemented!", error)
                CASE(WAVELET0D)
                   ! Zero the full stress tensor
                   h_stress = 0.0_dp
                END SELECT
             END SELECT
          END IF
          DO i = 1, 3
             CALL pw_pool_give_back_pw ( pw_pool, dvg ( i )%pw, error=error )
          END DO
       END IF

       CALL pw_pool_give_back_pw (pw_pool, rhog, error=error )
    ELSE
       ! stop on failure ?!
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END IF

    CALL timestop(handle)

  END SUBROUTINE pw_poisson_solve

! *****************************************************************************
!> \brief sets cell, grids and parameters used by the poisson solver
!>      You should call this at least once (and set everything)
!>      before using the poisson solver.
!>      Smart, doesn't set the thing twice to the same value
!>      Keeps track of the need to rebuild the poisson_env
!> \note
!>      Checks everything at the end. This means that after *each* call to
!>      this method the poisson env must be fully ready, so the first time
!>      you have to set everything at once. Change this behaviour?
!> \author fawzi
! *****************************************************************************
  SUBROUTINE pw_poisson_set ( poisson_env, cell_hmat, parameters, pw_pools,&
       use_level, mt_super_ref_pw_grid, force_rebuild, error )

    TYPE(pw_poisson_type), POINTER           :: poisson_env
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN), OPTIONAL                   :: cell_hmat
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: parameters
    TYPE(pw_pool_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: pw_pools
    INTEGER, INTENT(in), OPTIONAL            :: use_level
    TYPE(pw_grid_type), OPTIONAL, POINTER    :: mt_super_ref_pw_grid
    LOGICAL, INTENT(in), OPTIONAL            :: force_rebuild
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_poisson_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure, same
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: tmp_pools

    IF (PRESENT(parameters)) THEN
       CPPrecondition(ASSOCIATED(parameters),cp_failure_level,routineP,error,failure)
       CALL section_vals_retain(parameters,error=error)
       CALL section_vals_release(poisson_env%parameters,error=error)
       poisson_env%parameters => parameters
    END IF
    IF (PRESENT(cell_hmat)) THEN
       IF (ANY(poisson_env%cell_hmat /= cell_hmat)) &
          CALL pw_poisson_cleanup(poisson_env,error=error)
       poisson_env%cell_hmat(:,:) = cell_hmat(:,:)
       poisson_env%rebuild=.TRUE.
    END IF
    IF (PRESENT(pw_pools)) THEN
       CPPrecondition(ASSOCIATED(pw_pools),cp_failure_level,routineP,error,failure)
       same=.FALSE.
       IF (ASSOCIATED(poisson_env%pw_pools)) THEN
          same=SIZE(poisson_env%pw_pools)==SIZE(pw_pools)
          IF (same) THEN
             DO i =1,SIZE(pw_pools)
                IF (poisson_env%pw_pools(i)%pool%id_nr/=&
                     pw_pools(i)%pool%id_nr) same=.FALSE.
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
       IF (force_rebuild) poisson_env%rebuild=.TRUE.
    END IF
    CALL pw_poisson_check(poisson_env,error=error)
  END SUBROUTINE pw_poisson_set

END MODULE pw_poisson_methods
