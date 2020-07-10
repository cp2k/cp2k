!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to calculate 2nd order kernels from a given response density in ao basis
!>      linear response scf
!> \par History
!>      created 08-2020 [Frederick Stein], Code by M. Iannuzzi
!> \author Frederick Stein
! **************************************************************************************************
MODULE qs_2nd_kernel_ao
   USE admm_types,                      ONLY: admm_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              cp_dbcsr_plus_fm_fm_t,&
                                              cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE dbcsr_api,                       ONLY: dbcsr_add,&
                                              dbcsr_copy,&
                                              dbcsr_create,&
                                              dbcsr_deallocate_matrix,&
                                              dbcsr_p_type,&
                                              dbcsr_release,&
                                              dbcsr_set,&
                                              dbcsr_type
   USE hfx_admm_utils,                  ONLY: tddft_hfx_matrix
   USE input_constants,                 ONLY: do_admm_aux_exch_func_none,&
                                              do_admm_basis_projection,&
                                              do_admm_exch_scaling_none,&
                                              do_admm_purify_none
   USE input_section_types,             ONLY: section_get_ival,&
                                              section_get_rval,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_zero
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
   USE qs_linres_types,                 ONLY: linres_control_type
   USE qs_p_env_types,                  ONLY: qs_p_env_type
   USE task_list_types,                 ONLY: task_list_type
   USE xc,                              ONLY: xc_calc_2nd_deriv
   USE xc_derivatives,                  ONLY: xc_functionals_get_needs
   USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
   USE xc_rho_set_types,                ONLY: xc_rho_set_create,&
                                              xc_rho_set_release,&
                                              xc_rho_set_type,&
                                              xc_rho_set_update
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ! *** Public subroutines ***
   PUBLIC :: build_dm_response
   PUBLIC :: apply_hfx_ao
   PUBLIC :: apply_xc_admm_ao
   PUBLIC :: admm_aux_reponse_density

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_2nd_kernel_ao'

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief This routine builds response density in dbcsr format
!> \param c0 coefficients of unperturbed system (not changed)
!> \param c1 coefficients of response (not changed)
!> \param dm response density matrix
! **************************************************************************************************
   SUBROUTINE build_dm_response(c0, c1, dm)
      !
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(IN)       :: c0, c1
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(INOUT)    :: dm

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_dm_response', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ispin, ncol, nspins

      nspins = SIZE(dm, 1)

      DO ispin = 1, nspins
         CALL dbcsr_set(dm(ispin)%matrix, 0.0_dp)
         CALL cp_fm_get_info(c0(ispin)%matrix, ncol_global=ncol)
         CALL cp_dbcsr_plus_fm_fm_t(dm(ispin)%matrix, &
                                    matrix_v=c0(ispin)%matrix, &
                                    matrix_g=c1(ispin)%matrix, &
                                    ncol=ncol, alpha=1.0_dp, &
                                    keep_sparsity=.TRUE.)
         CALL cp_dbcsr_plus_fm_fm_t(dm(ispin)%matrix, &
                                    matrix_v=c1(ispin)%matrix, &
                                    matrix_g=c0(ispin)%matrix, &
                                    ncol=ncol, alpha=1.0_dp, &
                                    keep_sparsity=.TRUE.)
      ENDDO

   END SUBROUTINE build_dm_response

! **************************************************************************************************
!> \brief THis routine applies the Hartree-Fock Exchange kernel to a perturbation density matrix considering ADMM
!> \param qs_env the Quickstep environment
!> \param rho1_ao perturbation density matrix in AO basis (does not change)
!> \param rho1_ao_admm ...
!> \param h1_mat the resulting kernel matrix in AO basis (on output, it will have added the hfx kernel)
!> \param recalc_integrals whether the integrals are to be recalculated (default: no)
! **************************************************************************************************
   SUBROUTINE apply_hfx_ao(qs_env, rho1_ao, rho1_ao_admm, h1_mat, recalc_integrals)
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN), &
         TARGET                                          :: rho1_ao, rho1_ao_admm, h1_mat
      LOGICAL, INTENT(IN), OPTIONAL                      :: recalc_integrals

      CHARACTER(LEN=*), PARAMETER :: routineN = 'apply_hfx_ao', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, nao, nao_aux, nspins
      LOGICAL                                            :: my_recalc_integrals
      REAL(KIND=dp)                                      :: alpha
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, rho1, work_hmat
      TYPE(dbcsr_type)                                   :: work_hmat2
      TYPE(dft_control_type), POINTER                    :: dft_control

      CALL timeset(routineN, handle)

      my_recalc_integrals = .FALSE.
      IF (PRESENT(recalc_integrals)) my_recalc_integrals = recalc_integrals

      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control)

      IF (dft_control%do_admm) THEN
         IF (dft_control%admm_control%purification_method /= do_admm_purify_none) THEN
            CPABORT("ADMM: Linear Response needs purification_method=none")
         END IF
         IF (dft_control%admm_control%scaling_model /= do_admm_exch_scaling_none) THEN
            CPABORT("ADMM: Linear Response needs scaling_model=none")
         END IF
         IF (dft_control%admm_control%method /= do_admm_basis_projection) THEN
            CPABORT("ADMM: Linear Response needs admm_method=basis_projection")
         END IF
         !
      END IF

      nspins = SIZE(rho1_ao)

      DO ispin = 1, nspins
         CPASSERT(ASSOCIATED(rho1_ao(ispin)%matrix))
         CPASSERT(ASSOCIATED(h1_mat(ispin)%matrix))
      END DO

      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control)

      nspins = dft_control%nspins

      IF (dft_control%do_admm) THEN
         rho1 => rho1_ao_admm
      ELSE
         rho1 => rho1_ao
      END IF

      NULLIFY (work_hmat)
      CALL dbcsr_allocate_matrix_set(work_hmat, nspins)
      DO ispin = 1, nspins
         ALLOCATE (work_hmat(ispin)%matrix)
         CALL dbcsr_create(work_hmat(ispin)%matrix, template=rho1(ispin)%matrix)
         CALL dbcsr_copy(work_hmat(ispin)%matrix, rho1(ispin)%matrix)
         CALL dbcsr_set(work_hmat(ispin)%matrix, 0.0_dp)
      END DO

      ! Calculate kernel
      CALL tddft_hfx_matrix(work_hmat, rho1, qs_env, .FALSE., my_recalc_integrals)

      alpha = 2.0_dp
      IF (nspins == 2) alpha = 1.0_dp

      ! Get admm kernel in AO basis and copy final result
      IF (dft_control%do_admm) THEN
         ! Create working matrix to allow adding hfx kernel to original kernel
         CALL get_qs_env(qs_env=qs_env, admm_env=admm_env, matrix_s=matrix_s)
         CALL dbcsr_create(work_hmat2, template=matrix_s(1)%matrix)
         CALL dbcsr_copy(work_hmat2, matrix_s(1)%matrix)
         CALL dbcsr_set(work_hmat2, 0.0_dp)

         CALL cp_fm_get_info(admm_env%A, nrow_global=nao_aux, ncol_global=nao)

         DO ispin = 1, nspins
            ! tv = K * A_admm
            CALL cp_dbcsr_sm_fm_multiply(work_hmat(ispin)%matrix, admm_env%A, admm_env%work_aux_orb, &
                                         ncol=nao, alpha=1.0_dp, beta=0.0_dp)
            ! fm_ao_ao = A_admm^T * tv
            CALL cp_gemm('T', 'N', nao, nao, nao_aux, 1.0_dp, admm_env%A, &
                         admm_env%work_aux_orb, 0.0_dp, admm_env%work_orb_orb)

            ! Copy kernel matrix to sparse matrix and add it to the actual kernel
            CALL dbcsr_set(work_hmat2, 0.0_dp)
            CALL copy_fm_to_dbcsr(admm_env%work_orb_orb, work_hmat2, keep_sparsity=.TRUE.)
            CALL dbcsr_add(h1_mat(ispin)%matrix, work_hmat2, 1.0_dp, alpha)
         END DO
         CALL dbcsr_release(work_hmat2)
      ELSE
         DO ispin = 1, nspins
            CALL dbcsr_add(h1_mat(ispin)%matrix, work_hmat(ispin)%matrix, 1.0_dp, alpha)
         END DO
      END IF

      CALL dbcsr_deallocate_matrix_set(work_hmat)

      CALL timestop(handle)

   END SUBROUTINE apply_hfx_ao

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param p_env perturbation environment
! **************************************************************************************************
   SUBROUTINE apply_xc_admm_ao(qs_env, p_env)
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      TYPE(qs_p_env_type), INTENT(IN), POINTER           :: p_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'apply_xc_admm_ao'

      INTEGER                                            :: handle, ispin, nao, nao_aux, nspins
      INTEGER, DIMENSION(2, 3)                           :: bo
      LOGICAL                                            :: lsd
      REAL(KIND=dp)                                      :: alpha, fac, tot_rho_aux
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(cp_fm_type), POINTER                          :: tc0
      TYPE(dbcsr_p_type)                                 :: xcmat
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dbcsr_type)                                   :: work_hmat
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(linres_control_type), POINTER                 :: linres_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho1_aux_g, rho1_aux_r, tau_pw, v_xc
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(section_vals_type), POINTER                   :: xc_fun_section, xc_section
      TYPE(task_list_type), POINTER                      :: task_list_aux_fit
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type), POINTER                     :: rho1_set

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control)

      IF (dft_control%admm_control%aux_exch_func == do_admm_aux_exch_func_none) THEN
         ! nothing to do
      ELSE
         CALL get_qs_env(qs_env=qs_env, linres_control=linres_control)
         CPASSERT(.NOT. dft_control%qs_control%gapw)
         CPASSERT(.NOT. dft_control%qs_control%gapw_xc)
         CPASSERT(.NOT. dft_control%qs_control%lrigpw)
         CPASSERT(.NOT. linres_control%lr_triplet)

         nspins = dft_control%nspins

         ! AUX basis contribution
         CALL get_qs_env(qs_env=qs_env, pw_env=pw_env)
         CPASSERT(ASSOCIATED(pw_env))
         CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
         ALLOCATE (v_xc(nspins))
         DO ispin = 1, nspins
            NULLIFY (v_xc(ispin)%pw)
            CALL pw_pool_create_pw(auxbas_pw_pool, v_xc(ispin)%pw, &
                                   use_data=REALDATA3D, in_space=REALSPACE)
            CALL pw_zero(v_xc(ispin)%pw)
         END DO
         NULLIFY (tau_pw)
         ! calculate the xc potential
         lsd = (nspins == 2)
         CALL get_qs_env(qs_env=qs_env, ks_env=ks_env, matrix_s_aux_fit=matrix_s, task_list_aux_fit=task_list_aux_fit)
         ALLOCATE (xcmat%matrix)
         CALL dbcsr_create(xcmat%matrix, template=matrix_s(1)%matrix)
         ALLOCATE (rho1_aux_r(nspins), rho1_aux_g(nspins))
         DO ispin = 1, nspins
            NULLIFY (rho1_aux_r(ispin)%pw, rho1_aux_g(ispin)%pw)
            CALL pw_pool_create_pw(auxbas_pw_pool, rho1_aux_r(ispin)%pw, &
                                   use_data=REALDATA3D, in_space=REALSPACE)
            CALL pw_pool_create_pw(auxbas_pw_pool, rho1_aux_g(ispin)%pw, &
                                   in_space=RECIPROCALSPACE, use_data=COMPLEXDATA1D)
         END DO

         DO ispin = 1, nspins
            CALL calculate_rho_elec(matrix_p=p_env%p1_admm(ispin)%matrix, &
                                    rho=rho1_aux_r(ispin), rho_gspace=rho1_aux_g(ispin), &
                                    total_rho=tot_rho_aux, ks_env=ks_env, &
                                    basis_type="AUX_FIT", &
                                    task_list_external=task_list_aux_fit)
         END DO

         CALL get_qs_env(qs_env, admm_env=admm_env)
         xc_section => admm_env%xc_section_aux
         bo = rho1_aux_r(1)%pw%pw_grid%bounds_local
         ! create the place where to store the argument for the functionals
         NULLIFY (rho1_set)
         CALL xc_rho_set_create(rho1_set, bo, &
                                rho_cutoff=section_get_rval(xc_section, "DENSITY_CUTOFF"), &
                                drho_cutoff=section_get_rval(xc_section, "GRADIENT_CUTOFF"), &
                                tau_cutoff=section_get_rval(xc_section, "TAU_CUTOFF"))

         xc_fun_section => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
         needs = xc_functionals_get_needs(xc_fun_section, lsd, .TRUE.)

         ! calculate the arguments needed by the functionals
         CALL xc_rho_set_update(rho1_set, rho1_aux_r, rho1_aux_g, tau_pw, needs, &
                                section_get_ival(xc_section, "XC_GRID%XC_DERIV"), &
                                section_get_ival(xc_section, "XC_GRID%XC_SMOOTH_RHO"), &
                                auxbas_pw_pool)

         fac = 0._dp
         CALL xc_calc_2nd_deriv(v_xc, p_env%kpp1_env%deriv_set_admm, &
                                p_env%kpp1_env%rho_set_admm, &
                                rho1_set, auxbas_pw_pool, xc_section=xc_section, &
                                tddfpt_fac=fac)
         CALL xc_rho_set_release(rho1_set)

         tc0 => admm_env%work_aux_orb
         CALL dbcsr_create(work_hmat, template=p_env%kpp1(1)%matrix)
         CALL dbcsr_copy(work_hmat, p_env%kpp1(1)%matrix)

         alpha = 1.0_dp
         IF (nspins == 1) alpha = 2.0_dp

         CALL cp_fm_get_info(admm_env%A, nrow_global=nao_aux, ncol_global=nao)
         DO ispin = 1, nspins
            v_xc(ispin)%pw%cr3d = v_xc(ispin)%pw%cr3d*v_xc(ispin)%pw%pw_grid%dvol
            CALL dbcsr_copy(xcmat%matrix, matrix_s(1)%matrix)
            CALL dbcsr_set(xcmat%matrix, 0.0_dp)
            CALL integrate_v_rspace(v_rspace=v_xc(ispin), hmat=xcmat, qs_env=qs_env, &
                                    calculate_forces=.FALSE., basis_type="AUX_FIT", &
                                    task_list_external=task_list_aux_fit)
            CALL cp_dbcsr_sm_fm_multiply(xcmat%matrix, admm_env%A, tc0, &
                                         ncol=nao, alpha=1.0_dp, beta=0.0_dp)
            CALL cp_gemm('T', 'N', nao, nao, nao_aux, 1.0_dp, admm_env%A, &
                         tc0, 0.0_dp, admm_env%work_orb_orb)
            CALL copy_fm_to_dbcsr(admm_env%work_orb_orb, work_hmat, keep_sparsity=.TRUE.)
            CALL dbcsr_add(p_env%kpp1(ispin)%matrix, work_hmat, 1.0_dp, alpha)
         END DO

         CALL dbcsr_release(work_hmat)

         DO ispin = 1, nspins
            CALL pw_pool_give_back_pw(auxbas_pw_pool, v_xc(ispin)%pw)
            CALL pw_pool_give_back_pw(auxbas_pw_pool, rho1_aux_r(ispin)%pw)
            CALL pw_pool_give_back_pw(auxbas_pw_pool, rho1_aux_g(ispin)%pw)
         END DO
         DEALLOCATE (v_xc, rho1_aux_r, rho1_aux_g)
         CALL dbcsr_deallocate_matrix(xcmat%matrix)

      END IF

      CALL timestop(handle)

   END SUBROUTINE apply_xc_admm_ao

! **************************************************************************************************
!> \brief Calculate ADMM auxiliary response density
!> \param qs_env ...
!> \param dm ...
!> \param dm_admm ...
! **************************************************************************************************
   SUBROUTINE admm_aux_reponse_density(qs_env, dm, dm_admm)
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN)       :: dm
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(INOUT)    :: dm_admm

      CHARACTER(LEN=*), PARAMETER :: routineN = 'admm_aux_reponse_density'

      INTEGER                                            :: handle, ispin, nao, nao_aux, ncol, nspins
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(dft_control_type), POINTER                    :: dft_control

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, admm_env=admm_env, dft_control=dft_control)

      nspins = dft_control%nspins

      CPASSERT(ASSOCIATED(admm_env%A))
      CPASSERT(ASSOCIATED(admm_env%work_orb_orb))
      CPASSERT(ASSOCIATED(admm_env%work_aux_orb))
      CPASSERT(ASSOCIATED(admm_env%work_aux_aux))
      CALL cp_fm_get_info(admm_env%A, nrow_global=nao_aux, ncol_global=nao)

      ! P1 -> AUX BASIS
      CALL cp_fm_get_info(admm_env%work_orb_orb, nrow_global=nao, ncol_global=ncol)
      DO ispin = 1, nspins
         CALL copy_dbcsr_to_fm(dm(ispin)%matrix, admm_env%work_orb_orb)
         CALL cp_gemm('N', 'N', nao_aux, ncol, nao, 1.0_dp, admm_env%A, &
                      admm_env%work_orb_orb, 0.0_dp, admm_env%work_aux_orb)
         CALL cp_gemm('N', 'T', nao_aux, nao_aux, nao, 1.0_dp, admm_env%A, &
                      admm_env%work_aux_orb, 0.0_dp, admm_env%work_aux_aux)
         CALL copy_fm_to_dbcsr(admm_env%work_aux_aux, dm_admm(ispin)%matrix, keep_sparsity=.TRUE.)
      ENDDO

      CALL timestop(handle)

   END SUBROUTINE admm_aux_reponse_density

END MODULE qs_2nd_kernel_ao
