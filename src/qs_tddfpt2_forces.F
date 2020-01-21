!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

MODULE qs_tddfpt2_forces
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_plus_fm_fm_t,&
                                              cp_dbcsr_sm_fm_multiply
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE dbcsr_api,                       ONLY: dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_type
   USE input_section_types,             ONLY: section_get_ival,&
                                              section_get_rval,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_scale,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_poisson_methods,              ONLY: pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_p_type
   USE qs_collocate_density,            ONLY: calculate_rho_elec
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integrate_potential,          ONLY: integrate_v_rspace
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_tddfpt2_types,                ONLY: tddfpt_ground_state_mos
   USE xc,                              ONLY: xc_calc_2nd_deriv,&
                                              xc_prep_2nd_deriv
   USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                              xc_dset_release
   USE xc_derivatives,                  ONLY: xc_functionals_get_needs
   USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
   USE xc_rho_set_types,                ONLY: xc_rho_set_create,&
                                              xc_rho_set_release,&
                                              xc_rho_set_type,&
                                              xc_rho_set_update
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt2_forces'

   PUBLIC :: tddfpt_forces, tddfpt_resvec1, tddfpt_resvec2

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Calculate direct tddft forces
!> \param qs_env ...
!> \param evects   guess vectors distributed across all processors (initialised on exit)
!> \param evals    guessed transition energies (initialised on exit)
!> \param istate ...
!> \par History
!>    * 01.2020 screated [JGH]
! **************************************************************************************************
   SUBROUTINE tddfpt_forces(qs_env, evects, evals, istate)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: evects
      REAL(kind=dp), DIMENSION(:), INTENT(in)            :: evals
      INTEGER, INTENT(in)                                :: istate

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_forces', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, nspins, nstates

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(qs_env))
      nspins = SIZE(evects, 1)
      nstates = SIZE(evals)
      CPASSERT(nstates >= istate)

      CALL timestop(handle)

   END SUBROUTINE tddfpt_forces

! **************************************************************************************************
!> \brief ...
!> \param evect ...
!> \param mos_occ ...
!> \param matrix_hz ...
! **************************************************************************************************
   SUBROUTINE tddfpt_resvec1(evect, mos_occ, matrix_hz)

      TYPE(cp_fm_type), POINTER                          :: evect, mos_occ
      TYPE(dbcsr_type), POINTER                          :: matrix_hz

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_resvec1', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, nao, norb
      TYPE(cp_fm_struct_type), POINTER                   :: fmstruct, fmstruct2
      TYPE(cp_fm_type), POINTER                          :: cxmat, xxmat

      CALL timeset(routineN, handle)
      ! X*X^T
      CALL cp_fm_get_info(mos_occ, nrow_global=nao, ncol_global=norb)
      CALL cp_dbcsr_plus_fm_fm_t(matrix_hz, matrix_v=evect, ncol=norb)
      ! X^T*X
      CALL cp_fm_get_info(evect, matrix_struct=fmstruct)
      NULLIFY (fmstruct2)
      CALL cp_fm_struct_create(fmstruct=fmstruct2, template_fmstruct=fmstruct, &
                               nrow_global=norb, ncol_global=norb)
      CALL cp_fm_create(xxmat, matrix_struct=fmstruct2)
      CALL cp_fm_struct_release(fmstruct2)
      CALL cp_fm_create(cxmat, matrix_struct=fmstruct)
      CALL cp_gemm('T', 'N', norb, norb, nao, 1.0_dp, evect, evect, 0.0_dp, xxmat)
      CALL cp_gemm('N', 'N', nao, norb, norb, 1.0_dp, mos_occ, xxmat, 0.0_dp, cxmat)
      CALL cp_fm_release(xxmat)
      ! C*C^T*XX
      CALL cp_dbcsr_plus_fm_fm_t(matrix_hz, matrix_v=evect, ncol=norb)
      CALL cp_dbcsr_plus_fm_fm_t(matrix_hz, matrix_v=mos_occ, matrix_g=cxmat, &
                                 ncol=norb, alpha=-1.0_dp)
      CALL cp_fm_release(cxmat)
      CALL timestop(handle)

   END SUBROUTINE tddfpt_resvec1

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param matrix_hz ...
!> \param gs_mos ...
!> \param cpmos ...
! **************************************************************************************************
   SUBROUTINE tddfpt_resvec2(qs_env, matrix_hz, gs_mos, cpmos)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_hz
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER                                         :: gs_mos
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: cpmos

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_resvec2', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, norb, nspins
      INTEGER, DIMENSION(2, 3)                           :: bo
      LOGICAL                                            :: lsd
      REAL(KIND=dp)                                      :: thartree, total_rho
      TYPE(cp_fm_type), POINTER                          :: mos
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: rho_tot_gspace, v_hartree_gspace, &
                                                            v_hartree_rspace
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_g, rho_r, tau_r, trho_g, trho_r, v_xc
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: input, xc_fun_section, xc_section
      TYPE(xc_derivative_set_type), POINTER              :: deriv_set
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type), POINTER                     :: rho1_set, rho_set

      CALL timeset(routineN, handle)

      nspins = SIZE(gs_mos)

      NULLIFY (pw_env)
      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env, ks_env=ks_env)
      CPASSERT(ASSOCIATED(pw_env))

      NULLIFY (auxbas_pw_pool, poisson_env)
      ! gets the tmp grids
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                      poisson_env=poisson_env)

      NULLIFY (v_hartree_gspace%pw, rho_tot_gspace%pw, v_hartree_rspace%pw)
      CALL pw_pool_create_pw(auxbas_pw_pool, v_hartree_gspace%pw, &
                             use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, rho_tot_gspace%pw, &
                             use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, v_hartree_rspace%pw, &
                             use_data=REALDATA3D, in_space=REALSPACE)

      ALLOCATE (trho_r(nspins), trho_g(nspins))
      DO ispin = 1, nspins
         NULLIFY (trho_r(ispin)%pw, trho_g(ispin)%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool, trho_r(ispin)%pw, &
                                use_data=REALDATA3D, in_space=REALSPACE)
         CALL pw_pool_create_pw(auxbas_pw_pool, trho_g(ispin)%pw, &
                                use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      ENDDO
      CALL pw_zero(rho_tot_gspace%pw)
      DO ispin = 1, nspins
         CALL calculate_rho_elec(ks_env=ks_env, matrix_p=matrix_hz(ispin)%matrix, &
                                 rho=trho_r(ispin), &
                                 rho_gspace=trho_g(ispin), &
                                 total_rho=total_rho)
         CALL pw_axpy(trho_g(ispin)%pw, rho_tot_gspace%pw)
      END DO
      ! calculate associated hartree potential
      CALL pw_poisson_solve(poisson_env, rho_tot_gspace%pw, thartree, &
                            v_hartree_gspace%pw)
      CALL pw_transfer(v_hartree_gspace%pw, v_hartree_rspace%pw)
      CALL pw_scale(v_hartree_rspace%pw, v_hartree_rspace%pw%pw_grid%dvol)

      ! Fxc*drho term
      CALL get_qs_env(qs_env, rho=rho)
      CALL qs_rho_get(rho, rho_r=rho_r, rho_g=rho_g, tau_r=tau_r)
      ALLOCATE (v_xc(nspins))
      DO ispin = 1, nspins
         NULLIFY (v_xc(ispin)%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool, v_xc(ispin)%pw, &
                                use_data=REALDATA3D, in_space=REALSPACE)
         CALL pw_zero(v_xc(ispin)%pw)
      END DO
      CALL get_qs_env(qs_env, input=input)
      xc_section => section_vals_get_subs_vals(input, "DFT%XC")
      lsd = (nspins == 2)
      NULLIFY (deriv_set, rho_set, rho1_set)
      CALL xc_prep_2nd_deriv(deriv_set, rho_set, rho_r, auxbas_pw_pool, xc_section=xc_section)
      bo = trho_r(1)%pw%pw_grid%bounds_local
      CALL xc_rho_set_create(rho1_set, bo, &
                             rho_cutoff=section_get_rval(xc_section, "DENSITY_CUTOFF"), &
                             drho_cutoff=section_get_rval(xc_section, "GRADIENT_CUTOFF"), &
                             tau_cutoff=section_get_rval(xc_section, "TAU_CUTOFF"))

      xc_fun_section => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
      needs = xc_functionals_get_needs(xc_fun_section, lsd, .TRUE.)

      ! calculate the arguments needed by the functionals
      CALL xc_rho_set_update(rho1_set, trho_r, trho_g, tau_r, needs, &
                             section_get_ival(xc_section, "XC_GRID%XC_DERIV"), &
                             section_get_ival(xc_section, "XC_GRID%XC_SMOOTH_RHO"), &
                             auxbas_pw_pool)
      CALL xc_calc_2nd_deriv(v_xc, deriv_set, rho_set, &
                             rho1_set, auxbas_pw_pool, xc_section=xc_section)
      CALL xc_dset_release(deriv_set)
      CALL xc_rho_set_release(rho_set)
      CALL xc_rho_set_release(rho1_set)

      ! vtot = v_xc(ispin) + v_hartree
      DO ispin = 1, nspins
         CALL dbcsr_set(matrix_hz(ispin)%matrix, 0.0_dp)
         CALL pw_scale(v_xc(ispin)%pw, v_xc(ispin)%pw%pw_grid%dvol)
         CALL pw_axpy(v_hartree_rspace%pw, v_xc(ispin)%pw)
         CALL integrate_v_rspace(qs_env=qs_env, v_rspace=v_xc(ispin), &
                                 hmat=matrix_hz(ispin), &
                                 calculate_forces=.FALSE.)
      END DO

      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_gspace%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_rspace%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_tot_gspace%pw)
      DO ispin = 1, nspins
         CALL pw_pool_give_back_pw(auxbas_pw_pool, trho_r(ispin)%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, trho_g(ispin)%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, v_xc(ispin)%pw)
      END DO
      DEALLOCATE (trho_r, trho_g, v_xc)

      DO ispin = 1, nspins
         mos => gs_mos(ispin)%mos_occ
         CALL cp_fm_get_info(mos, ncol_global=norb)
         CALL cp_dbcsr_sm_fm_multiply(matrix_hz(ispin)%matrix, mos, cpmos(ispin)%matrix, &
                                      norb, alpha=1.0_dp, beta=0.0_dp)
      END DO

      CALL timestop(handle)

   END SUBROUTINE tddfpt_resvec2

END MODULE qs_tddfpt2_forces
