!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE qs_tddfpt2_operators
   USE admm_types,                      ONLY: admm_type
   USE cp_dbcsr_operations,             ONLY: copy_fm_to_dbcsr,&
                                              cp_dbcsr_sm_fm_multiply
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                              cp_fm_scale_and_add
   USE cp_fm_struct,                    ONLY: cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE dbcsr_api,                       ONLY: dbcsr_p_type,&
                                              dbcsr_set
   USE hfx_admm_utils,                  ONLY: tddft_hfx_matrix
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_poisson_methods,              ONLY: pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_p_type,&
                                              pw_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kernel_types,                 ONLY: full_kernel_env_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_tddfpt2_types,                ONLY: tddfpt_ground_state_mos
   USE xc,                              ONLY: xc_calc_2nd_deriv_analytical,&
                                              xc_calc_2nd_deriv_numerical
   USE xc_rho_set_types,                ONLY: xc_rho_set_type,&
                                              xc_rho_set_update
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt2_operators'

   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .FALSE.
   ! number of first derivative components (3: d/dx, d/dy, d/dz)
   INTEGER, PARAMETER, PRIVATE          :: nderivs = 3
   INTEGER, PARAMETER, PRIVATE          :: maxspins = 2

   PUBLIC :: tddfpt_apply_energy_diff, tddfpt_apply_coulomb, tddfpt_apply_xc, tddfpt_apply_hfx

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Apply orbital energy difference term:
!>        Aop_evects(spin,state) += KS(spin) * evects(spin,state) -
!>                                  S * evects(spin,state) * diag(evals_occ(spin))
!> \param Aop_evects  action of TDDFPT operator on trial vectors (modified on exit)
!> \param evects      trial vectors C_{1,i}
!> \param S_evects    S * C_{1,i}
!> \param gs_mos      molecular orbitals optimised for the ground state (only occupied orbital
!>                    energies [component %evals_occ] are needed)
!> \param matrix_ks   Kohn-Sham matrix
!> \par History
!>    * 05.2016 initialise all matrix elements in one go [Sergey Chulkov]
!>    * 03.2017 renamed from tddfpt_init_energy_diff(), altered prototype [Sergey Chulkov]
!> \note Based on the subroutine p_op_l1() which was originally created by
!>       Thomas Chassaing on 08.2002.
! **************************************************************************************************
   SUBROUTINE tddfpt_apply_energy_diff(Aop_evects, evects, S_evects, gs_mos, matrix_ks)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: Aop_evects, evects, S_evects
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(in)       :: matrix_ks

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_apply_energy_diff'

      INTEGER                                            :: handle, ispin, ivect, nactive, nao, &
                                                            nspins, nvects
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct
      TYPE(cp_fm_type), POINTER                          :: hevec

      CALL timeset(routineN, handle)

      nspins = SIZE(evects, 1)
      nvects = SIZE(evects, 2)

      DO ispin = 1, nspins
         CALL cp_fm_get_info(matrix=evects(ispin, 1)%matrix, matrix_struct=matrix_struct, &
                             nrow_global=nao, ncol_global=nactive)
         NULLIFY (hevec)
         CALL cp_fm_create(hevec, matrix_struct)

         DO ivect = 1, nvects
            CALL cp_dbcsr_sm_fm_multiply(matrix_ks(ispin)%matrix, evects(ispin, ivect)%matrix, &
                                         Aop_evects(ispin, ivect)%matrix, ncol=nactive, &
                                         alpha=1.0_dp, beta=1.0_dp)

            IF (ASSOCIATED(gs_mos(ispin)%evals_occ_matrix)) THEN
               ! orbital energy correction: evals_occ_matrix is not a diagonal matrix
               CALL cp_gemm('N', 'N', nao, nactive, nactive, 1.0_dp, &
                            S_evects(ispin, ivect)%matrix, gs_mos(ispin)%evals_occ_matrix, &
                            0.0_dp, hevec)
            ELSE
               CALL cp_fm_to_fm(S_evects(ispin, ivect)%matrix, hevec)
               CALL cp_fm_column_scale(hevec, gs_mos(ispin)%evals_occ)
            END IF

            ! KS * C1 - S * C1 * occupied_orbital_energies
            CALL cp_fm_scale_and_add(1.0_dp, Aop_evects(ispin, ivect)%matrix, -1.0_dp, hevec)
         END DO
         CALL cp_fm_release(hevec)
      END DO

      CALL timestop(handle)

   END SUBROUTINE tddfpt_apply_energy_diff

! **************************************************************************************************
!> \brief Update v_rspace by adding coulomb term.
!> \param A_ia_rspace    action of TDDFPT operator on the trial vector expressed in a plane wave
!>                       representation (modified on exit)
!> \param rho_ia_g       response density in reciprocal space for the given trial vector
!> \param pw_env         plain wave environment
!> \param work_v_gspace  work reciprocal-space grid to store Coulomb potential (modified on exit)
!> \param work_v_rspace  work real-space grid to store Coulomb potential (modified on exit)
!> \par History
!>    * 05.2016 compute all coulomb terms in one go [Sergey Chulkov]
!>    * 03.2017 proceed excited states sequentially; minimise the number of conversions between
!>              DBCSR and FM matrices [Sergey Chulkov]
!>    * 06.2018 return the action expressed in the plane wave representation instead of the one
!>              in the atomic basis set representation
!> \note Based on the subroutine kpp1_calc_k_p_p1() which was originally created by
!>       Mohamed Fawzi on 10.2002.
! **************************************************************************************************
   SUBROUTINE tddfpt_apply_coulomb(A_ia_rspace, rho_ia_g, pw_env, work_v_gspace, work_v_rspace)
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: A_ia_rspace
      TYPE(pw_type), POINTER                             :: rho_ia_g
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), INTENT(inout)                     :: work_v_gspace, work_v_rspace

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_apply_coulomb'

      INTEGER                                            :: handle, ispin, nspins
      REAL(kind=dp)                                      :: alpha, pair_energy
      TYPE(pw_poisson_type), POINTER                     :: poisson_env

      CALL timeset(routineN, handle)

      nspins = SIZE(A_ia_rspace)
      CALL pw_env_get(pw_env, poisson_env=poisson_env)

      IF (nspins > 1) THEN
         alpha = 1.0_dp
      ELSE
         ! spin-restricted case: alpha == 2 due to singlet state.
         ! In case of triplet states alpha == 0, so we should not call this subroutine at all.
         alpha = 2.0_dp
      END IF

      CALL pw_poisson_solve(poisson_env, rho_ia_g, pair_energy, work_v_gspace%pw)
      CALL pw_transfer(work_v_gspace%pw, work_v_rspace%pw)

      ! (i a || j b) = ( i_alpha a_alpha + i_beta a_beta || j_alpha b_alpha + j_beta b_beta) =
      !                tr (Cj_alpha^T * [J_i{alpha}a{alpha}_munu + J_i{beta}a{beta}_munu] * Cb_alpha) +
      !                tr (Cj_beta^T * [J_i{alpha}a{alpha}_munu + J_i{beta}a{beta}_munu] * Cb_beta)
      DO ispin = 1, nspins
         CALL pw_axpy(work_v_rspace%pw, A_ia_rspace(ispin)%pw, alpha)
      END DO

      CALL timestop(handle)

   END SUBROUTINE tddfpt_apply_coulomb

! **************************************************************************************************
!> \brief Driver routine for applying fxc (analyic vs. finite difference for testing
!> \param A_ia_rspace      action of TDDFPT operator on trial vectors expressed in a plane wave
!>                         representation (modified on exit)
!> \param kernel_env       kernel environment
!> \param rho_ia_struct    response density for the given trial vector
!> \param is_rks_triplets  indicates that the triplet excited states calculation using
!>                         spin-unpolarised molecular orbitals has been requested
!> \param pw_env           plain wave environment
!> \param work_v_xc        work real-space grid to store the gradient of the exchange-correlation
!>                         potential with respect to the response density (modified on exit)
!> \param work_v_xc_tau ...
! **************************************************************************************************
   SUBROUTINE tddfpt_apply_xc(A_ia_rspace, kernel_env, rho_ia_struct, is_rks_triplets, pw_env, work_v_xc, work_v_xc_tau)
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: A_ia_rspace
      TYPE(full_kernel_env_type), INTENT(in)             :: kernel_env
      TYPE(qs_rho_type), POINTER                         :: rho_ia_struct
      LOGICAL, INTENT(in)                                :: is_rks_triplets
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: work_v_xc, work_v_xc_tau

      INTEGER                                            :: ispin

      IF (kernel_env%deriv2_analytic) THEN
         CALL tddfpt_apply_xc_analytic(A_ia_rspace, kernel_env, rho_ia_struct, is_rks_triplets, &
                                       pw_env, work_v_xc, work_v_xc_tau)
      ELSE
         CALL tddfpt_apply_xc_fd(A_ia_rspace, kernel_env, rho_ia_struct, is_rks_triplets, &
                                 pw_env, work_v_xc, work_v_xc_tau)
      END IF

      DO ispin = 1, SIZE(work_v_xc)
         ! pw2 = pw2 + alpha * pw1
         CALL pw_axpy(work_v_xc(ispin)%pw, A_ia_rspace(ispin)%pw, kernel_env%alpha)
      END DO

   END SUBROUTINE tddfpt_apply_xc

! **************************************************************************************************
!> \brief Update A_ia_munu by adding exchange-correlation term.
!> \param A_ia_rspace      action of TDDFPT operator on trial vectors expressed in a plane wave
!>                         representation (modified on exit)
!> \param kernel_env       kernel environment
!> \param rho_ia_struct    response density for the given trial vector
!> \param is_rks_triplets  indicates that the triplet excited states calculation using
!>                         spin-unpolarised molecular orbitals has been requested
!> \param pw_env           plain wave environment
!> \param work_v_xc        work real-space grid to store the gradient of the exchange-correlation
!>                         potential with respect to the response density (modified on exit)
!> \param work_v_xc_tau ...
!> \par History
!>    * 05.2016 compute all kernel terms in one go [Sergey Chulkov]
!>    * 03.2017 proceed excited states sequentially; minimise the number of conversions between
!>              DBCSR and FM matrices [Sergey Chulkov]
!>    * 06.2018 return the action expressed in the plane wave representation instead of the one
!>              in the atomic basis set representation
!> \note Based on the subroutine kpp1_calc_k_p_p1() which was originally created by
!>       Mohamed Fawzi on 10.2002.
! **************************************************************************************************
   SUBROUTINE tddfpt_apply_xc_analytic(A_ia_rspace, kernel_env, rho_ia_struct, is_rks_triplets, pw_env, work_v_xc, work_v_xc_tau)
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: A_ia_rspace
      TYPE(full_kernel_env_type), INTENT(in)             :: kernel_env
      TYPE(qs_rho_type), POINTER                         :: rho_ia_struct
      LOGICAL, INTENT(in)                                :: is_rks_triplets
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: work_v_xc, work_v_xc_tau

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_apply_xc_analytic'

      INTEGER                                            :: handle, ispin, nspins
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_ia_g, rho_ia_g2, rho_ia_r, &
                                                            rho_ia_r2, tau_ia_r, tau_ia_r2
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool

      CALL timeset(routineN, handle)

      nspins = SIZE(A_ia_rspace)
      CALL qs_rho_get(rho_ia_struct, rho_g=rho_ia_g, rho_r=rho_ia_r, tau_r=tau_ia_r)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      IF (debug_this_module) THEN
         CPASSERT(SIZE(rho_ia_g) == nspins)
         CPASSERT(SIZE(rho_ia_r) == nspins)
         CPASSERT((.NOT. ASSOCIATED(tau_ia_r)) .OR. SIZE(tau_ia_r) == nspins)
         CPASSERT((.NOT. is_rks_triplets) .OR. nspins == 1)
      END IF

      NULLIFY (tau_ia_r2)
      IF (is_rks_triplets) THEN
         ALLOCATE (rho_ia_r2(2))
         ALLOCATE (rho_ia_g2(2))
         rho_ia_r2(1)%pw => rho_ia_r(1)%pw
         rho_ia_r2(2)%pw => rho_ia_r(1)%pw
         rho_ia_g2(1)%pw => rho_ia_g(1)%pw
         rho_ia_g2(2)%pw => rho_ia_g(1)%pw

         IF (ASSOCIATED(tau_ia_r)) THEN
            ALLOCATE (tau_ia_r2(2))
            tau_ia_r2(1)%pw => tau_ia_r(1)%pw
            tau_ia_r2(2)%pw => tau_ia_r(1)%pw
         END IF
      ELSE
         ALLOCATE (rho_ia_r2(nspins))
         ALLOCATE (rho_ia_g2(nspins))
         DO ispin = 1, nspins
            rho_ia_r2(ispin)%pw => rho_ia_r(ispin)%pw
            rho_ia_g2(ispin)%pw => rho_ia_g(ispin)%pw
         END DO

         IF (ASSOCIATED(tau_ia_r)) THEN
            ALLOCATE (tau_ia_r2(nspins))
            DO ispin = 1, nspins
               tau_ia_r2(ispin)%pw => tau_ia_r(ispin)%pw
            END DO
         END IF
      END IF

      DO ispin = 1, nspins
         CALL pw_zero(work_v_xc(ispin)%pw)
         IF (ASSOCIATED(work_v_xc_tau)) CALL pw_zero(work_v_xc_tau(ispin)%pw)
      END DO

      CALL xc_rho_set_update(rho_set=kernel_env%xc_rho1_set, rho_r=rho_ia_r2, rho_g=rho_ia_g2, tau=tau_ia_r2, &
                             needs=kernel_env%xc_rho1_cflags, xc_deriv_method_id=kernel_env%deriv_method_id, &
                             xc_rho_smooth_id=kernel_env%rho_smooth_id, pw_pool=auxbas_pw_pool)

      CALL xc_calc_2nd_deriv_analytical(v_xc=work_v_xc, v_xc_tau=work_v_xc_tau, deriv_set=kernel_env%xc_deriv_set, &
                                        rho_set=kernel_env%xc_rho_set, &
                                        rho1_set=kernel_env%xc_rho1_set, pw_pool=auxbas_pw_pool, &
                                        xc_section=kernel_env%xc_section, gapw=.FALSE., tddfpt_fac=kernel_env%beta)

      DEALLOCATE (rho_ia_g2, rho_ia_r2)

      CALL timestop(handle)

   END SUBROUTINE tddfpt_apply_xc_analytic

! **************************************************************************************************
!> \brief Update A_ia_munu by adding exchange-correlation term using finite difference methods.
!> \param A_ia_rspace      action of TDDFPT operator on trial vectors expressed in a plane wave
!>                         representation (modified on exit)
!> \param kernel_env       kernel environment
!> \param rho_ia_struct    response density for the given trial vector
!> \param is_rks_triplets  indicates that the triplet excited states calculation using
!>                         spin-unpolarised molecular orbitals has been requested
!> \param pw_env           plain wave environment
!> \param work_v_xc        work real-space grid to store the gradient of the exchange-correlation
!>                         potential with respect to the response density (modified on exit)
!> \param work_v_xc_tau ...
! **************************************************************************************************
   SUBROUTINE tddfpt_apply_xc_fd(A_ia_rspace, kernel_env, rho_ia_struct, is_rks_triplets, pw_env, work_v_xc, work_v_xc_tau)
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: A_ia_rspace
      TYPE(full_kernel_env_type), INTENT(in)             :: kernel_env
      TYPE(qs_rho_type), POINTER                         :: rho_ia_struct
      LOGICAL, INTENT(in)                                :: is_rks_triplets
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: work_v_xc, work_v_xc_tau

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_apply_xc_fd'

      INTEGER                                            :: handle, ispin, nspins
      LOGICAL                                            :: lsd, singlet, triplet
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho1_g, rho1_r, tau1_r
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(xc_rho_set_type), POINTER                     :: rho_set

      CALL timeset(routineN, handle)

      nspins = SIZE(A_ia_rspace)
      CALL qs_rho_get(rho_ia_struct, rho_r=rho1_r, rho_g=rho1_g, tau_r=tau1_r)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
      DO ispin = 1, nspins
         CALL pw_zero(work_v_xc(ispin)%pw)
      END DO
      rho_set => kernel_env%xc_rho_set

      singlet = .FALSE.
      triplet = .FALSE.
      lsd = .FALSE.
      IF (nspins == 1 .AND. .NOT. is_rks_triplets) THEN
         singlet = .TRUE.
      ELSE IF (nspins == 1 .AND. is_rks_triplets) THEN
         triplet = .TRUE.
      ELSE IF (nspins == 2) THEN
         lsd = .TRUE.
      ELSE
         CPABORT("illegal options")
      END IF

      IF (ASSOCIATED(tau1_r)) THEN
      DO ispin = 1, nspins
         CALL pw_zero(work_v_xc_tau(ispin)%pw)
      END DO
      END IF

      CALL xc_calc_2nd_deriv_numerical(work_v_xc, work_v_xc_tau, rho_set, rho1_r, rho1_g, tau1_r, &
                                       auxbas_pw_pool, kernel_env%xc_section, &
                                       is_rks_triplets)

      CALL timestop(handle)

   END SUBROUTINE tddfpt_apply_xc_fd

! **************************************************************************************************
!> \brief Update action of TDDFPT operator on trial vectors by adding exact-exchange term.
!> \param Aop_evects      action of TDDFPT operator on trial vectors (modified on exit)
!> \param evects          trial vectors
!> \param gs_mos          molecular orbitals optimised for the ground state (only occupied
!>                        molecular orbitals [component %mos_occ] are needed)
!> \param do_admm         perform auxiliary density matrix method calculations
!> \param qs_env          Quickstep environment
!> \param work_rho_ia_ao  work sparse matrix with shape [nao x nao] distributed globally
!>                        to store response density (modified on exit)
!> \param work_hmat       work sparse matrix with shape [nao x nao] distributed globally
!>                        (modified on exit)
!> \param wfm_rho_orb     work dense matrix with shape [nao x nao] distributed globally
!>                        (modified on exit)
!> \par History
!>    * 05.2016 compute all exact-exchange terms in one go [Sergey Chulkov]
!>    * 03.2017 code related to ADMM correction is now moved to tddfpt_apply_admm_correction()
!>              in order to compute this correction within parallel groups [Sergey Chulkov]
!> \note Based on the subroutine kpp1_calc_k_p_p1() which was originally created by
!>       Mohamed Fawzi on 10.2002.
! **************************************************************************************************
   SUBROUTINE tddfpt_apply_hfx(Aop_evects, evects, gs_mos, do_admm, qs_env, &
                               work_rho_ia_ao, work_hmat, wfm_rho_orb)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: Aop_evects, evects
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      LOGICAL, INTENT(in)                                :: do_admm
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(INOUT)    :: work_rho_ia_ao
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(INOUT), &
         TARGET                                          :: work_hmat
      TYPE(cp_fm_type), INTENT(IN), POINTER              :: wfm_rho_orb

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'tddfpt_apply_hfx'

      INTEGER                                            :: handle, ispin, ivect, nao, nao_aux, &
                                                            nspins, nvects
      INTEGER, DIMENSION(maxspins)                       :: nactive
      LOGICAL                                            :: do_hfx
      REAL(kind=dp)                                      :: alpha
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(section_vals_type), POINTER                   :: hfx_section, input

      CALL timeset(routineN, handle)

      ! Check for hfx section
      CALL get_qs_env(qs_env, input=input)
      hfx_section => section_vals_get_subs_vals(input, "DFT%XC%HF")
      CALL section_vals_get(hfx_section, explicit=do_hfx)

      IF (do_hfx) THEN
         nspins = SIZE(evects, 1)
         nvects = SIZE(evects, 2)

         IF (nspins > 1) THEN
            alpha = 1.0_dp
         ELSE
            alpha = 2.0_dp
         END IF

         CALL cp_fm_get_info(gs_mos(1)%mos_occ, nrow_global=nao)
         DO ispin = 1, nspins
            CALL cp_fm_get_info(evects(ispin, 1)%matrix, ncol_global=nactive(ispin))
         END DO

         IF (do_admm) THEN
            CALL get_qs_env(qs_env, admm_env=admm_env)
            CALL cp_fm_get_info(admm_env%A, nrow_global=nao_aux)
         END IF

         ! some stuff from qs_ks_build_kohn_sham_matrix
         ! TO DO: add SIC support
         DO ivect = 1, nvects
            DO ispin = 1, nspins
               CALL cp_gemm('N', 'T', nao, nao, nactive(ispin), 1.0_dp, gs_mos(ispin)%mos_occ, &
                            evects(ispin, ivect)%matrix, 0.0_dp, wfm_rho_orb)
               CALL cp_gemm('N', 'T', nao, nao, nactive(ispin), 1.0_dp, evects(ispin, ivect)%matrix, &
                            gs_mos(ispin)%mos_occ, 1.0_dp, wfm_rho_orb)

               CALL dbcsr_set(work_hmat(ispin)%matrix, 0.0_dp)
               IF (do_admm) THEN
                  CALL cp_gemm('N', 'N', nao_aux, nao, nao, 1.0_dp, admm_env%A, &
                               wfm_rho_orb, 0.0_dp, admm_env%work_aux_orb)
                  CALL cp_gemm('N', 'T', nao_aux, nao_aux, nao, 1.0_dp, admm_env%A, admm_env%work_aux_orb, &
                               0.0_dp, admm_env%work_aux_aux)
                  CALL copy_fm_to_dbcsr(admm_env%work_aux_aux, work_rho_ia_ao(ispin)%matrix, keep_sparsity=.TRUE.)
               ELSE
                  CALL copy_fm_to_dbcsr(wfm_rho_orb, work_rho_ia_ao(ispin)%matrix, keep_sparsity=.TRUE.)
               END IF
            END DO

            CALL tddft_hfx_matrix(work_hmat, work_rho_ia_ao, qs_env)

            IF (do_admm) THEN
               DO ispin = 1, nspins
                  CALL cp_dbcsr_sm_fm_multiply(work_hmat(ispin)%matrix, admm_env%A, admm_env%work_aux_orb, &
                                               ncol=nao, alpha=1.0_dp, beta=0.0_dp)

                  CALL cp_gemm('T', 'N', nao, nao, nao_aux, 1.0_dp, admm_env%A, &
                               admm_env%work_aux_orb, 0.0_dp, wfm_rho_orb)

                  CALL cp_gemm('N', 'N', nao, nactive(ispin), nao, alpha, wfm_rho_orb, &
                               gs_mos(ispin)%mos_occ, 1.0_dp, Aop_evects(ispin, ivect)%matrix)
               END DO
            ELSE
               DO ispin = 1, nspins
                  CALL cp_dbcsr_sm_fm_multiply(work_hmat(ispin)%matrix, gs_mos(ispin)%mos_occ, &
                                               Aop_evects(ispin, ivect)%matrix, ncol=nactive(ispin), &
                                               alpha=alpha, beta=1.0_dp)
               END DO
            END IF
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE tddfpt_apply_hfx

END MODULE qs_tddfpt2_operators
