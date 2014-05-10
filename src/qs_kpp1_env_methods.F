!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief module that builds the second order perturbation kernel
!>      kpp1 = delta_rho|_P delta_rho|_P E drho(P1) drho
!> \par History
!>      07.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE qs_kpp1_env_methods
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_add,&
                                             cp_dbcsr_allocate_matrix_set,&
                                             cp_dbcsr_copy,&
                                             cp_dbcsr_deallocate_matrix_set,&
                                             cp_dbcsr_init,&
                                             cp_dbcsr_p_type,&
                                             cp_dbcsr_scale,&
                                             cp_dbcsr_set
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE f77_blas
  USE hartree_local_methods,           ONLY: Vh_1c_gg_integrals
  USE input_constants,                 ONLY: do_method_gapw,&
                                             do_method_gapw_xc,&
                                             tddfpt_excitations,&
                                             tddfpt_singlet,&
                                             tddfpt_triplet
  USE input_section_types,             ONLY: section_get_ival,&
                                             section_get_rval,&
                                             section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kahan_sum,                       ONLY: accurate_sum
  USE kinds,                           ONLY: dp
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_methods,                      ONLY: pw_axpy,&
                                             pw_copy,&
                                             pw_integrate_function,&
                                             pw_transfer,&
                                             pw_zero
  USE pw_poisson_methods,              ONLY: pw_poisson_solve
  USE pw_poisson_types,                ONLY: pw_poisson_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_p_type,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_create,&
                                             pw_p_type,&
                                             pw_release,&
                                             pw_retain
  USE qs_energy_types,                 ONLY: allocate_qs_energy,&
                                             deallocate_qs_energy,&
                                             qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_gapw_densities,               ONLY: prepare_gapw_den
  USE qs_integrate_potential,          ONLY: integrate_v_rspace
  USE qs_kpp1_env_types,               ONLY: qs_kpp1_env_type
  USE qs_ks_atom,                      ONLY: update_ks_atom
  USE qs_ks_methods,                   ONLY: qs_ks_build_kohn_sham_matrix
  USE qs_p_env_types,                  ONLY: qs_p_env_type
  USE qs_rho0_ggrid,                   ONLY: integrate_vhg0_rspace
  USE qs_rho_types,                    ONLY: qs_rho_get,&
                                             qs_rho_type
  USE qs_vxc_atom,                     ONLY: calculate_xc_2nd_deriv_atom
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc,                              ONLY: xc_calc_2nd_deriv,&
                                             xc_prep_2nd_deriv
  USE xc_derivative_set_types,         ONLY: xc_dset_release
  USE xc_derivatives,                  ONLY: xc_functionals_get_needs
  USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
  USE xc_rho_set_types,                ONLY: xc_rho_set_create,&
                                             xc_rho_set_release,&
                                             xc_rho_set_type,&
                                             xc_rho_set_update
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_kpp1_env_methods'
  INTEGER, PRIVATE, SAVE :: last_kpp1_id_nr=0

  PUBLIC :: kpp1_create, &
            kpp1_calc_k_p_p1, &
            kpp1_calc_k_p_p1_fdiff, &
            kpp1_did_change

CONTAINS

! *****************************************************************************
!> \brief allocates and initializes a kpp1_env
!> \param kpp 1_env: the environement to initialize
!> \param qs_env the qs environment in which this environment lives
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      07.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE kpp1_create(kpp1_env,qs_env, error)
    TYPE(qs_kpp1_env_type), POINTER          :: kpp1_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kpp1_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    ALLOCATE(kpp1_env,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       NULLIFY(kpp1_env%v_rspace, kpp1_env%v_ao, kpp1_env%drho_r,&
            kpp1_env%rho_set, &
            kpp1_env%deriv_set, kpp1_env%spin_pot, kpp1_env%grad_pot,&
            kpp1_env%ndiag_term)
       kpp1_env%ref_count=1
       last_kpp1_id_nr=last_kpp1_id_nr+1
       kpp1_env%id_nr=last_kpp1_id_nr
       kpp1_env%iter=0
       kpp1_env%print_count=0
    END IF
  END SUBROUTINE kpp1_create

! *****************************************************************************
!> \brief calcualtes the k_p_p1 kernel of the pertubation theory
!> \param kpp 1_env: the kpp1 environement in which perform the calculation
!> \param qs_env kpp1's qs_env
!> \param k_p_p 1: the sparse matrix that will contain the kernel k_p_p1
!> \param rho the density where to evaluate the derivatives (i.e. p along
!>        with with its grid representations, that must be valid)
!> \param rho 1: the density that represent the first direction along which
!>        you should evaluate the derivatives
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE kpp1_calc_k_p_p1(kpp1_env, p_env, qs_env, k_p_p1, rho, rho1, rho1_xc, error)

    TYPE(qs_kpp1_env_type), POINTER          :: kpp1_env
    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: k_p_p1
    TYPE(qs_rho_type), POINTER               :: rho, rho1
    TYPE(qs_rho_type), OPTIONAL, POINTER     :: rho1_xc
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kpp1_calc_k_p_p1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: excitations, handle, ispin, &
                                                nspins, output_unit, &
                                                res_etype, stat
    INTEGER, DIMENSION(2, 3)                 :: bo
    LOGICAL                                  :: explicit, failure, gapw, &
                                                gapw_xc, ionode, lsd, &
                                                lsd_singlets
    REAL(KIND=dp)                            :: energy_hartree, &
                                                energy_hartree_1c, fac
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: rho1_tot_gspace, &
                                                v_hartree_gspace, &
                                                v_hartree_rspace
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho1_g_pw, rho1_r, rho1_r_pw, &
                                                tau_pw, v_rspace_new, v_xc
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(section_vals_type), POINTER         :: input, scf_section, &
                                                xc_fun_section, xc_section
    TYPE(xc_rho_cflags_type)                 :: needs
    TYPE(xc_rho_set_type), POINTER           :: rho1_set

    CALL timeset(routineN,handle)
    failure=.FALSE.

    NULLIFY(auxbas_pw_pool, pw_pools, pw_env, v_rspace_new, &
            rho1_r, rho1_g_pw, tau_pw, v_xc, rho1_set,&
            poisson_env, input, scf_section)
    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source

    CPPrecondition(ASSOCIATED(kpp1_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(k_p_p1),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(rho),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(rho1),cp_failure_level,routineP,error,failure)

    IF (.NOT.failure) THEN
       CPPrecondition(kpp1_env%ref_count>0,cp_failure_level,routineP,error,failure)
       CALL kpp1_check_i_alloc(kpp1_env,qs_env=qs_env,error=error)
       CALL get_qs_env(qs_env=qs_env,&
            pw_env=pw_env,&
            input=input,&
            error=error)

       gapw=(section_get_ival(input,"DFT%QS%METHOD",error=error)==do_method_gapw)
       gapw_xc=(section_get_ival(input,"DFT%QS%METHOD",error=error)==do_method_gapw_xc)
       IF(gapw_xc) THEN
          CPPrecondition(ASSOCIATED(rho1_xc),cp_failure_level,routineP,error,failure)
       END IF

       nspins = SIZE(k_p_p1)
       lsd = (nspins==2)

       xc_section => section_vals_get_subs_vals(input,"DFT%XC",error=error)
       scf_section => section_vals_get_subs_vals(input,"DFT%SCF",error=error)
       CALL section_vals_val_get(input,"DFT%EXCITATIONS",&
            i_val=excitations,error=error)
       IF (excitations==tddfpt_excitations) THEN
          xc_section => section_vals_get_subs_vals(input,"DFT%TDDFPT%XC",error=error)
          !FM this check should already had happened and section made explicit, give an error?
          CALL section_vals_get(xc_section,explicit=explicit,error=error)
          IF (.NOT.explicit) THEN
             xc_section => section_vals_get_subs_vals(input,"DFT%XC",error=error)
          END IF
       END IF

       CALL section_vals_val_get(input,"DFT%TDDFPT%LSD_SINGLETS",&
            l_val=lsd_singlets,error=error)
       CALL section_vals_val_get(input,"DFT%TDDFPT%RES_ETYPE",&
            i_val=res_etype,error=error)
    END IF

    IF (.NOT.failure) THEN
       kpp1_env%iter=kpp1_env%iter+1
    END IF

! gets the tmp grids
    IF (.NOT. failure) THEN
       CPPrecondition(ASSOCIATED(pw_env),cp_failure_level,routineP,error,failure)
       CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,&
            pw_pools=pw_pools, poisson_env=poisson_env,error=error)
       ALLOCATE(v_rspace_new(nspins), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT.failure) THEN
       CALL pw_pool_create_pw(auxbas_pw_pool,v_hartree_gspace%pw,&
            use_data = COMPLEXDATA1D,&
            in_space = RECIPROCALSPACE, error=error)
       CALL pw_pool_create_pw(auxbas_pw_pool,v_hartree_rspace%pw,&
            use_data = REALDATA3D,&
            in_space = REALSPACE, error=error)
    END IF

    IF (gapw .OR. gapw_xc) &
       CALL prepare_gapw_den(qs_env,p_env%local_rho_set, do_rho0=(.NOT.gapw_xc), error=error)

! *** calculate the hartree potential on the total density ***
    IF (.NOT. failure) THEN

       CALL pw_pool_create_pw(auxbas_pw_pool, rho1_tot_gspace%pw,&
            use_data = COMPLEXDATA1D,&
            in_space = RECIPROCALSPACE, error=error)

       CALL pw_copy(rho1%rho_g(1)%pw,rho1_tot_gspace%pw,error=error)
       DO ispin=2,nspins
          CALL pw_axpy(rho1%rho_g(ispin)%pw, rho1_tot_gspace%pw, error=error)
       END DO
       IF (gapw) &
            CALL pw_axpy(p_env%local_rho_set%rho0_mpole%rho0_s_gs%pw, rho1_tot_gspace%pw,&
            error=error)

       IF (cp_print_key_should_output(logger%iter_info,scf_section,"PRINT%TOTAL_DENSITIES",&
            error=error)/=0) THEN
          output_unit = cp_print_key_unit_nr(logger,scf_section,"PRINT%TOTAL_DENSITIES",&
               extension=".scfLog",error=error)
          CALL print_densities(kpp1_env, rho1, rho1_tot_gspace, output_unit, error=error)
          CALL cp_print_key_finished_output(output_unit,logger,scf_section,&
               "PRINT%TOTAL_DENSITIES", error=error)
       END IF

       IF (.NOT.(nspins==1 .AND. excitations==tddfpt_excitations .AND. &
                 res_etype /= tddfpt_singlet )) THEN
          CALL pw_poisson_solve(poisson_env,rho1_tot_gspace%pw, &
                                 energy_hartree, &
                                 v_hartree_gspace%pw,error=error)
          CALL pw_transfer(v_hartree_gspace%pw,v_hartree_rspace%pw,error=error)
       END IF

       CALL pw_pool_give_back_pw(auxbas_pw_pool, rho1_tot_gspace%pw,&
            error=error)

! *** calculate the xc potential ***
       IF(gapw_xc) THEN
         CALL qs_rho_get(rho1_xc, rho_r=rho1_r, error=error)
       ELSE
         CALL qs_rho_get(rho1, rho_r=rho1_r, error=error)
       END IF

       IF (nspins == 1 .AND. excitations==tddfpt_excitations .AND. &
           (lsd_singlets .OR. res_etype == tddfpt_triplet)) THEN

          lsd = .TRUE.
          ALLOCATE(rho1_r_pw(2))
          DO ispin=1, 2
             NULLIFY(rho1_r_pw(ispin)%pw)
             CALL pw_create(rho1_r_pw(ispin)%pw, rho1_r(1)%pw%pw_grid, &
                  rho1_r(1)%pw%in_use, rho1_r(1)%pw%in_space,error=error)
             CALL pw_transfer(rho1_r(1)%pw, rho1_r_pw(ispin)%pw,error=error)
          END DO

       ELSE

          ALLOCATE(rho1_r_pw(nspins))
          DO ispin=1, nspins
             rho1_r_pw(ispin)%pw => rho1_r(ispin)%pw
             CALL pw_retain(rho1_r_pw(ispin)%pw,error=error)
          END DO

       END IF

       NULLIFY(tau_pw)

       !------!
       ! rho1 !
       !------!
       bo = rho1_r(1)%pw%pw_grid%bounds_local
       ! create the place where to store the argument for the functionals
       CALL xc_rho_set_create(rho1_set, bo, &
                              rho_cutoff=section_get_rval(xc_section, "DENSITY_CUTOFF",error=error), &
                              drho_cutoff=section_get_rval(xc_section, "GRADIENT_CUTOFF",error=error), &
                              tau_cutoff=section_get_rval(xc_section, "TAU_CUTOFF",error=error), &
                              error=error)

       xc_fun_section => section_vals_get_subs_vals(xc_section,"XC_FUNCTIONAL",&
            error=error)
       needs=xc_functionals_get_needs(xc_fun_section,lsd,.TRUE.,error)

       ! calculate the arguments needed by the functionals
       CALL xc_rho_set_update(rho1_set, rho1_r_pw, rho1_g_pw, tau_pw, needs,&
                             section_get_ival(xc_section,"XC_GRID%XC_DERIV",error=error),&
                             section_get_ival(xc_section,"XC_GRID%XC_SMOOTH_RHO",error=error),&
                             auxbas_pw_pool,  error)

       ALLOCATE(v_xc(nspins),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO ispin=1, nspins
          NULLIFY(v_xc(ispin)%pw)
          CALL pw_pool_create_pw(auxbas_pw_pool,v_xc(ispin)%pw,&
                                 use_data = REALDATA3D,&
                                 in_space = REALSPACE, error=error)
          CALL pw_zero(v_xc(ispin)%pw,error=error)
       END DO

       fac=0._dp
       IF (nspins==1.AND.excitations==tddfpt_excitations) THEN
          IF (lsd_singlets) fac=1.0_dp
          IF (res_etype == tddfpt_triplet) fac=-1.0_dp
       END IF

       CALL  xc_calc_2nd_deriv(v_xc, kpp1_env%deriv_set, kpp1_env%rho_set, &
            rho1_set, auxbas_pw_pool,xc_section=xc_section,&
            tddfpt_fac=fac, error=error)

       DO ispin=1,nspins
          v_rspace_new(ispin)%pw => v_xc(ispin)%pw
       END DO
       DEALLOCATE(v_xc,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       IF (gapw) CALL calculate_xc_2nd_deriv_atom(p_env, qs_env, xc_section, error=error)

       CALL xc_rho_set_release(rho1_set,error=error)

       DO ispin=1,SIZE(rho1_r_pw)
          CALL pw_release(rho1_r_pw(ispin)%pw,error=error)
       END DO
       DEALLOCATE(rho1_r_pw, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       !-------------------------------!
       ! Add both hartree and xc terms !
       !-------------------------------!
       DO ispin=1,nspins

          IF(gapw_xc) THEN
          ! XC and Hartree are integrated separatedly
          ! XC uses the sofft basis set only
            v_rspace_new(ispin)%pw%cr3d = v_rspace_new(ispin)%pw%cr3d * &
                                        v_rspace_new(ispin)%pw%pw_grid%dvol

            IF (excitations==tddfpt_excitations .AND. nspins==1) THEN

                IF (.NOT.(lsd_singlets .OR. &
                       res_etype == tddfpt_triplet)) THEN

                   v_rspace_new(1)%pw%cr3d = 2.0_dp * v_rspace_new(1)%pw%cr3d

                END IF
                ! remove kpp1_env%v_ao and work directly on k_p_p1 ?
                CALL cp_dbcsr_set(kpp1_env%v_ao(ispin)%matrix,0.0_dp,error=error)
                CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin),&
                     p=rho%rho_ao(ispin),&
                     h=kpp1_env%v_ao(ispin),&
                     qs_env=qs_env,&
                     calculate_forces=.FALSE.,gapw=gapw_xc,error=error)

               ! add hartree only for SINGLETS
                IF (res_etype == tddfpt_singlet) THEN
                    v_hartree_rspace%pw%cr3d = v_hartree_rspace%pw%cr3d * &
                                               v_hartree_rspace%pw%pw_grid%dvol
                    v_rspace_new(1)%pw%cr3d  = 2.0_dp * v_hartree_rspace%pw%cr3d

                    CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin),&
                         p=rho%rho_ao(ispin),&
                         h=kpp1_env%v_ao(ispin),&
                         qs_env=qs_env,&
                         calculate_forces=.FALSE.,gapw=gapw,error=error)
                END IF
             ELSE
               ! remove kpp1_env%v_ao and work directly on k_p_p1 ?
               CALL cp_dbcsr_set(kpp1_env%v_ao(ispin)%matrix,0.0_dp,error=error)
               CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin),&
                    p=rho%rho_ao(ispin),&
                    h=kpp1_env%v_ao(ispin),&
                    qs_env=qs_env,&
                    calculate_forces=.FALSE.,gapw=gapw_xc,error=error)

               IF (ispin == 1) THEN
                  v_hartree_rspace%pw%cr3d = v_hartree_rspace%pw%cr3d * &
                                             v_hartree_rspace%pw%pw_grid%dvol
               END IF
               v_rspace_new(ispin)%pw%cr3d  = v_hartree_rspace%pw%cr3d
               CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin),&
                    p=rho%rho_ao(ispin),&
                    h=kpp1_env%v_ao(ispin),&
                    qs_env=qs_env,&
                    calculate_forces=.FALSE.,gapw=gapw,error=error)
             END IF

          ELSE
            v_rspace_new(ispin)%pw%cr3d = v_rspace_new(ispin)%pw%cr3d * &
                                          v_rspace_new(ispin)%pw%pw_grid%dvol

            IF (excitations==tddfpt_excitations .AND. nspins==1) THEN

               IF (.NOT.(lsd_singlets .OR. &
                         res_etype == tddfpt_triplet)) THEN

                  v_rspace_new(1)%pw%cr3d = 2.0_dp * v_rspace_new(1)%pw%cr3d

               END IF

               ! add hartree only for SINGLETS
               IF (res_etype == tddfpt_singlet) THEN
                  v_hartree_rspace%pw%cr3d = v_hartree_rspace%pw%cr3d * &
                                             v_hartree_rspace%pw%pw_grid%dvol
                  v_rspace_new(1)%pw%cr3d  = v_rspace_new(1)%pw%cr3d + &
                                             2.0_dp * v_hartree_rspace%pw%cr3d
               END IF
            ELSE
               IF (ispin == 1) THEN
                  v_hartree_rspace%pw%cr3d = v_hartree_rspace%pw%cr3d * &
                                             v_hartree_rspace%pw%pw_grid%dvol
               END IF
               v_rspace_new(ispin)%pw%cr3d  = v_rspace_new(ispin)%pw%cr3d + &
                                              v_hartree_rspace%pw%cr3d
            END IF

            ! remove kpp1_env%v_ao and work directly on k_p_p1 ?
            CALL cp_dbcsr_set(kpp1_env%v_ao(ispin)%matrix,0.0_dp,error=error)
            CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin),&
                 p=rho%rho_ao(ispin),&
                 h=kpp1_env%v_ao(ispin),&
                 qs_env=qs_env,&
                 calculate_forces=.FALSE.,gapw=gapw,error=error)

          END IF

          CALL cp_dbcsr_copy(k_p_p1(ispin)%matrix,kpp1_env%v_ao(ispin)%matrix,error=error)
       END DO

       IF (gapw) THEN
          IF (.NOT. (excitations==tddfpt_excitations .AND. (nspins == 1 .AND. &
               res_etype == tddfpt_triplet))) THEN
             CALL Vh_1c_gg_integrals(qs_env,energy_hartree_1c, .TRUE., p_env=p_env,error=error)
             CALL integrate_vhg0_rspace(qs_env, v_hartree_rspace, &
                                        .FALSE., .TRUE., p_env=p_env, error=error)
          END IF
!         ***  Add single atom contributions to the KS matrix ***
          CALL update_ks_atom(qs_env,k_p_p1,rho%rho_ao,.FALSE.,.TRUE.,p_env,error)
       END IF

       CALL pw_pool_give_back_pw(auxbas_pw_pool,v_hartree_gspace%pw,&
            error=error)
       CALL pw_pool_give_back_pw(auxbas_pw_pool,v_hartree_rspace%pw,&
            error=error)
       DO ispin=1,nspins
          CALL pw_pool_give_back_pw(auxbas_pw_pool,v_rspace_new(ispin)%pw,&
               error=error)
       END DO
       DEALLOCATE(v_rspace_new, stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF

    CALL timestop(handle)
  END SUBROUTINE kpp1_calc_k_p_p1

! *****************************************************************************
!> \brief calcualtes the k_p_p1 kernel of the pertubation theory with finite
!>      differences
!> \param kpp 1_env: the kpp1 environement in which perform the calculation
!> \param qs_env kpp1's qs_env
!> \param k_p_p 1: the sparse matrix that will contain the kernel k_p_p1
!> \param rho the density where to evaluate the derivatives (i.e. p along
!>        with with its grid representations, that must be valid)
!> \param rho 1: the density that represent the first direction along which
!>        you should evaluate the derivatives
!> \param diff the amount of the finite difference step
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      useful for testing purposes.
!>      rescale my_diff depending on the norm of rho1?
!> \par History
!>      01.2003 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE kpp1_calc_k_p_p1_fdiff(kpp1_env,qs_env,k_p_p1,rho,rho1,&
       diff, error)
    TYPE(qs_kpp1_env_type), POINTER          :: kpp1_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: k_p_p1
    TYPE(qs_rho_type), POINTER               :: rho, rho1
    REAL(KIND=dp), INTENT(in), OPTIONAL      :: diff
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kpp1_calc_k_p_p1_fdiff', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin, nspins
    REAL(KIND=dp)                            :: my_diff
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: ks_2, matrix_s, rho_ao
    TYPE(qs_energy_type), POINTER            :: qs_energy

    NULLIFY(ks_2,matrix_s,qs_energy)
    nspins = SIZE(k_p_p1)
    my_diff=1.0e-6_dp
    IF (PRESENT(diff)) my_diff=diff
    CALL allocate_qs_energy(qs_energy)

    CALL qs_rho_get(rho, rho_ao=rho_ao,error=error)
    CALL get_qs_env(qs_env, &
         matrix_s=matrix_s,error=error)

    ! rho = rho0+h/2*rho1
    my_diff=my_diff/2.0_dp
    DO ispin=1,SIZE(k_p_p1)
       CALL cp_dbcsr_add(rho%rho_ao(ispin)%matrix,rho1%rho_ao(ispin)%matrix,&
            alpha_scalar=1.0_dp,beta_scalar=my_diff,error=error)
       rho%rho_r(ispin)%pw%cr3d=rho%rho_r(ispin)%pw%cr3d+&
            my_diff*rho1%rho_r(ispin)%pw%cr3d
       rho%rho_g(ispin)%pw%cc=rho%rho_g(ispin)%pw%cc+&
            my_diff*rho1%rho_g(ispin)%pw%cc
    END DO

    CALL qs_ks_build_kohn_sham_matrix(ks_env=qs_env%ks_env,&
         qs_env=qs_env,ks_matrix=k_p_p1,rho=rho,energy=qs_energy,&
         calculate_forces=.FALSE.,&
         just_energy=.FALSE.,error=error)

    CALL cp_dbcsr_allocate_matrix_set(ks_2,nspins,error=error)
    DO ispin=1,nspins
       ALLOCATE(ks_2(ispin)%matrix)
       CALL cp_dbcsr_init(ks_2(ispin)%matrix,error=error)
       CALL cp_dbcsr_copy(ks_2(ispin)%matrix,matrix_s(1)%matrix,&
            name="tmp_ks2-"//ADJUSTL(cp_to_string(ispin)),error=error)
    END DO

    ! rho = rho0-h/2*rho1
    my_diff=-2.0_dp*my_diff
    DO ispin=1,nspins
       CALL cp_dbcsr_add(rho%rho_ao(ispin)%matrix,rho1%rho_ao(ispin)%matrix,&
            alpha_scalar=1.0_dp,beta_scalar=my_diff,error=error)
       rho%rho_r(ispin)%pw%cr3d=rho%rho_r(ispin)%pw%cr3d+&
            my_diff*rho1%rho_r(ispin)%pw%cr3d
       rho%rho_g(ispin)%pw%cc=rho%rho_g(ispin)%pw%cc+&
            my_diff*rho1%rho_g(ispin)%pw%cc
    END DO

    CALL qs_ks_build_kohn_sham_matrix(ks_env=qs_env%ks_env,&
         qs_env=qs_env,ks_matrix=ks_2,rho=rho,energy=qs_energy,&
         calculate_forces=.FALSE.,&
         just_energy=.FALSE.,error=error)

    ! rho = rho0
    my_diff=-0.5_dp*my_diff
    DO ispin=1,nspins
       CALL cp_dbcsr_add(rho%rho_ao(ispin)%matrix,rho1%rho_ao(ispin)%matrix,&
            alpha_scalar=1.0_dp,beta_scalar=my_diff,error=error)
       rho%rho_r(ispin)%pw%cr3d=rho%rho_r(ispin)%pw%cr3d+&
            my_diff*rho1%rho_r(ispin)%pw%cr3d
       rho%rho_g(ispin)%pw%cc=rho%rho_g(ispin)%pw%cc+&
            my_diff*rho1%rho_g(ispin)%pw%cc
    END DO

    ! k_p_p1=(H(rho0+h/2 rho1)-H(rho0-h/2 rho1))/h
    DO ispin=1,nspins
       CALL cp_dbcsr_add(k_p_p1(ispin)%matrix,ks_2(ispin)%matrix,&
            alpha_scalar=1.0_dp,beta_scalar=-1.0_dp,error=error)
       CALL cp_dbcsr_scale(k_p_p1(ispin)%matrix,alpha_scalar=0.5_dp/my_diff,error=error)
    END DO

    CALL cp_dbcsr_deallocate_matrix_set(ks_2,error=error)
    CALL deallocate_qs_energy(qs_energy)
  END SUBROUTINE kpp1_calc_k_p_p1_fdiff

! *****************************************************************************
!> \brief checks that the intenal storage is allocated, and allocs it if needed
!> \param kpp 1_env: the environment to check
!> \param qs_env the qs environment this kpp1_env lives in
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      private routine
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE kpp1_check_i_alloc(kpp1_env, qs_env, error)

    TYPE(qs_kpp1_env_type), POINTER          :: kpp1_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kpp1_check_i_alloc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: excitations, ispin, nspins, &
                                                res_etype, stat
    LOGICAL                                  :: explicit, failure, &
                                                lsd_singlets
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: my_rho_r, rho_r
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(section_vals_type), POINTER         :: input, xc_section

! ------------------------------------------------------------------

    failure=.FALSE.
    NULLIFY(pw_env,auxbas_pw_pool,matrix_s,rho,rho_r,input)

    CPPrecondition(ASSOCIATED(kpp1_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(kpp1_env%ref_count>0,cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN

       CALL get_qs_env(qs_env, pw_env=pw_env,&
                       matrix_s=matrix_s, input=input, error=error, rho=rho)

       CALL qs_rho_get(rho, rho_r=rho_r, error=error)
       nspins=SIZE(rho_r)

       CALL pw_env_get(pw_env, auxbas_pw_pool = auxbas_pw_pool, error=error)

       IF (.NOT.ASSOCIATED(kpp1_env%v_rspace)) THEN
          ALLOCATE(kpp1_env%v_rspace(nspins),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          IF (.NOT.failure) THEN
             DO ispin=1,nspins
                CALL pw_pool_create_pw(auxbas_pw_pool, &
                     kpp1_env%v_rspace(ispin)%pw,&
                     use_data=REALDATA3D, in_space=REALSPACE,error=error)
             END DO
          END IF
       END IF

       IF (.NOT.ASSOCIATED(kpp1_env%v_ao)) THEN
          CALL cp_dbcsr_allocate_matrix_set(kpp1_env%v_ao,nspins,error=error)
             DO ispin=1,nspins
                ALLOCATE(kpp1_env%v_ao(ispin)%matrix)
                CALL cp_dbcsr_init(kpp1_env%v_ao(ispin)%matrix,error=error)
                CALL cp_dbcsr_copy(kpp1_env%v_ao(ispin)%matrix,matrix_s(1)%matrix,&
                     name="kpp1%v_ao-"//ADJUSTL(cp_to_string(ispin)),error=error)
             END DO
       END IF

       IF (.not.ASSOCIATED(kpp1_env%deriv_set)) THEN

          CALL section_vals_val_get(input,"DFT%EXCITATIONS",&
               i_val=excitations,error=error)
          CALL section_vals_val_get(input,"DFT%TDDFPT%LSD_SINGLETS",&
               l_val=lsd_singlets,error=error)
          CALL section_vals_val_get(input,"DFT%TDDFPT%RES_ETYPE",&
               i_val=res_etype,error=error)
          IF (nspins==1.AND.(excitations==tddfpt_excitations.and.&
               (lsd_singlets .OR. res_etype == tddfpt_triplet))) THEN
             ALLOCATE(my_rho_r(2),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DO ispin=1,2
                CALL pw_pool_create_pw(auxbas_pw_pool,my_rho_r(ispin)%pw, &
                     use_data=rho_r(1)%pw%in_use, in_space=rho_r(1)%pw%in_space,&
                     error=error)
                my_rho_r(ispin)%pw%cr3d = 0.5_dp * rho_r(1)%pw%cr3d
             END DO
          ELSE
             ALLOCATE(my_rho_r(SIZE(rho_r)),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DO ispin=1,SIZE(rho_r)
                my_rho_r(ispin)%pw => rho_r(ispin)%pw
                CALL pw_retain(my_rho_r(ispin)%pw,error=error)
             END DO
          END IF

          xc_section => section_vals_get_subs_vals(input,"DFT%XC",error=error)
          CALL section_vals_val_get(input,"DFT%EXCITATIONS",&
               i_val=excitations,error=error)
          IF (excitations==tddfpt_excitations) THEN
             xc_section => section_vals_get_subs_vals(input,"DFT%TDDFPT%XC",error=error)
!FM this check should already had happened and section made explicit, give an error?
             CALL section_vals_get(xc_section,explicit=explicit,error=error)
             IF (.NOT.explicit) THEN
                xc_section => section_vals_get_subs_vals(input,"DFT%XC",error=error)
             END IF
          END IF

          CALL xc_prep_2nd_deriv(kpp1_env%deriv_set, kpp1_env%rho_set, &
               my_rho_r, auxbas_pw_pool, &
               xc_section=xc_section, error=error)

          DO ispin=1,SIZE(my_rho_r)
             CALL pw_release(my_rho_r(ispin)%pw,error=error)
          END DO
          DEALLOCATE(my_rho_r,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF

  END SUBROUTINE kpp1_check_i_alloc

! *****************************************************************************
!> \brief function to advise of changes either in the grids or in the the
!>      matrix structure (to update kpp1_env internal cache)
!> \param kpp 1_env: the kpp1_env
!> \param grid_changed true if the grids have changed and have to be rebuilt
!> \param s_struct_changed true if the structure of the s matrix has changed
!> \param psi 0_changed: true if the value of psi0 (qs_env%rho) has changed
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      11.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE kpp1_did_change(kpp1_env, s_struct_changed, grid_changed,&
       psi0_changed, error)
    TYPE(qs_kpp1_env_type), POINTER          :: kpp1_env
    LOGICAL, INTENT(in), OPTIONAL            :: s_struct_changed, &
                                                grid_changed, psi0_changed
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kpp1_did_change', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure, my_psi0_changed, &
                                                my_s_struct_changed

    failure=.FALSE.
    my_s_struct_changed=.FALSE.
    my_psi0_changed=.FALSE.

    IF (PRESENT(s_struct_changed)) my_s_struct_changed=s_struct_changed
    IF (PRESENT(psi0_changed)) my_psi0_changed=psi0_changed

    CPPrecondition(ASSOCIATED(kpp1_env),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPrecondition(kpp1_env%ref_count>0,cp_failure_level,routineP,error,failure)
    END IF
    IF (.not.failure) THEN
       IF (my_s_struct_changed) THEN
          IF (ASSOCIATED(kpp1_env%v_ao)) THEN
             CALL cp_dbcsr_deallocate_matrix_set(kpp1_env%v_ao,error=error)
          END IF
       END IF
       IF (my_s_struct_changed.or.my_psi0_changed) THEN
          IF (ASSOCIATED(kpp1_env%drho_r)) THEN
             DEALLOCATE(kpp1_env%drho_r, stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF (ASSOCIATED(kpp1_env%deriv_set)) THEN
             CALL xc_dset_release(kpp1_env%deriv_set, error=error)
             NULLIFY(kpp1_env%deriv_set)
          END IF
          IF (ASSOCIATED(kpp1_env%spin_pot)) THEN
             DEALLOCATE(kpp1_env%spin_pot, stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF (ASSOCIATED(kpp1_env%grad_pot)) THEN
             DEALLOCATE(kpp1_env%grad_pot, stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF (ASSOCIATED(kpp1_env%ndiag_term)) THEN
             DEALLOCATE(kpp1_env%ndiag_term, stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          CALL xc_rho_set_release(kpp1_env%rho_set,error=error) ! it would be better to pass a pw pool
       END IF
       IF (PRESENT(grid_changed)) THEN
          IF (grid_changed) THEN
             IF (ASSOCIATED(kpp1_env%v_rspace)) THEN
                DO i=1,SIZE(kpp1_env%v_rspace)
                   CALL pw_release(kpp1_env%v_rspace(i)%pw,error=error)
                END DO
                DEALLOCATE(kpp1_env%v_rspace,stat=stat)
                CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             END IF
          END IF
       END IF
    END IF
  END SUBROUTINE kpp1_did_change

! *****************************************************************************
  SUBROUTINE print_densities(kpp1_env, rho1, rho1_tot_gspace, output_unit, error)

    TYPE(qs_kpp1_env_type), POINTER          :: kpp1_env
    TYPE(qs_rho_type), POINTER               :: rho1
    TYPE(pw_p_type), INTENT(IN)              :: rho1_tot_gspace
    INTEGER                                  :: output_unit
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'print_densities', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: total_rho_gspace

    total_rho_gspace = pw_integrate_function(rho1_tot_gspace%pw,isign=-1,error=error)
    IF (output_unit>0) THEN
       WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
            "KPP1 total charge density (r-space):",&
            accurate_sum(rho1%tot_rho_r),&
            "KPP1 total charge density (g-space):",&
            total_rho_gspace
    END IF

  END SUBROUTINE print_densities

END MODULE qs_kpp1_env_methods
