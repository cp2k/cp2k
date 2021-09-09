!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief linres kernel functions
!> \par History
!>      created from qs_linres_methods
!> \author JGH
! **************************************************************************************************
MODULE qs_linres_kernel
   USE admm_types,                      ONLY: admm_type
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_plus_fm_fm_t,&
                                              cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_copy,&
                                              dbcsr_create,&
                                              dbcsr_deallocate_matrix,&
                                              dbcsr_p_type,&
                                              dbcsr_set
   USE hartree_local_methods,           ONLY: Vh_1c_gg_integrals
   USE hfx_energy_potential,            ONLY: integrate_four_center
   USE hfx_ri,                          ONLY: hfx_ri_update_ks
   USE hfx_types,                       ONLY: hfx_type
   USE input_constants,                 ONLY: do_admm_aux_exch_func_none,&
                                              do_admm_basis_projection,&
                                              do_admm_exch_scaling_none,&
                                              do_admm_purify_none,&
                                              kg_tnadd_embed
   USE input_cp2k_hfx,                  ONLY: ri_mo
   USE input_section_types,             ONLY: section_get_ival,&
                                              section_get_lval,&
                                              section_get_rval,&
                                              section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kg_correction,                   ONLY: kg_ekin_subset
   USE kg_environment_types,            ONLY: kg_environment_type
   USE kinds,                           ONLY: dp
   USE lri_environment_types,           ONLY: lri_density_type,&
                                              lri_environment_type,&
                                              lri_kind_type
   USE lri_ks_methods,                  ONLY: calculate_lri_ks_matrix
   USE message_passing,                 ONLY: mp_sum
   USE mulliken,                        ONLY: ao_charges
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_copy,&
                                              pw_transfer
   USE pw_poisson_methods,              ONLY: pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_p_type,&
                                              pw_release,&
                                              pw_retain
   USE qs_collocate_density,            ONLY: calculate_rho_elec
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_fxc,                          ONLY: qs_fxc_analytic,&
                                              qs_fxc_fdiff
   USE qs_gapw_densities,               ONLY: prepare_gapw_den
   USE qs_integrate_potential,          ONLY: integrate_v_rspace,&
                                              integrate_v_rspace_diagonal,&
                                              integrate_v_rspace_one_center
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_kpp1_env_types,               ONLY: qs_kpp1_env_type
   USE qs_ks_atom,                      ONLY: update_ks_atom
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_linres_types,                 ONLY: linres_control_type
   USE qs_p_env_types,                  ONLY: qs_p_env_type
   USE qs_rho0_ggrid,                   ONLY: integrate_vhg0_rspace
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_vxc_atom,                     ONLY: calculate_xc_2nd_deriv_atom
   USE task_list_types,                 ONLY: task_list_type
   USE xc,                              ONLY: xc_calc_2nd_deriv,&
                                              xc_prep_2nd_deriv
   USE xc_derivatives,                  ONLY: xc_functionals_get_needs
   USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
   USE xc_rho_set_types,                ONLY: xc_rho_set_create,&
                                              xc_rho_set_release,&
                                              xc_rho_set_type,&
                                              xc_rho_set_update
   USE xtb_ehess,                       ONLY: xtb_coulomb_hessian
   USE xtb_types,                       ONLY: get_xtb_atom_param,&
                                              xtb_atom_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ! *** Public subroutines ***
   PUBLIC :: apply_op_2

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_linres_kernel'

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param p_env ...
!> \param c0 ...
!> \param v ...
!> \param Av ...
! **************************************************************************************************
   SUBROUTINE apply_op_2(qs_env, p_env, c0, v, Av)
      !
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_p_env_type), POINTER                       :: p_env
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: c0, v, Av

      TYPE(dft_control_type), POINTER                    :: dft_control

      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control)
      IF (dft_control%qs_control%semi_empirical) THEN
         CPABORT("Linear response not available with SE methods")
      ELSEIF (dft_control%qs_control%dftb) THEN
         CPABORT("Linear response not available with DFTB")
      ELSEIF (dft_control%qs_control%xtb) THEN
         CALL apply_op_2_xtb(qs_env, p_env, c0, v, Av)
      ELSE
         CALL apply_op_2_dft(qs_env, p_env, c0, v, Av)
         CALL apply_hfx(qs_env, p_env, c0, v, Av)
         CALL apply_xc_admm(qs_env, p_env, c0, v, Av)
      END IF

   END SUBROUTINE apply_op_2

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param p_env ...
!> \param c0 ...
!> \param v ...
!> \param Av ...
! **************************************************************************************************
   SUBROUTINE apply_op_2_dft(qs_env, p_env, c0, v, Av)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_p_env_type), POINTER                       :: p_env
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: c0, v, Av

      CHARACTER(len=*), PARAMETER                        :: routineN = 'apply_op_2_dft'

      INTEGER                                            :: handle, ikind, ispin, ncol, nkind, ns, &
                                                            nspins
      LOGICAL                                            :: deriv2_analytic, gapw, gapw_xc, &
                                                            lr_triplet, lrigpw
      REAL(KIND=dp)                                      :: alpha, ekin_mol, energy_hartree, &
                                                            energy_hartree_1c
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: k1mat, matrix_s, rho1_ao, rho_ao
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ksmat, psmat
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kg_environment_type), POINTER                 :: kg_env
      TYPE(linres_control_type), POINTER                 :: linres_control
      TYPE(lri_density_type), POINTER                    :: lri_density
      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(lri_kind_type), DIMENSION(:), POINTER         :: lri_v_int
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: rho1_tot_gspace, v_hartree_gspace, &
                                                            v_hartree_rspace
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho1_g, rho1_r, rho_r, tau, &
                                                            v_rspace_new, v_xc
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_kpp1_env_type), POINTER                    :: kpp1_env
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho, rho1, rho1_xc, rho1a, rho_aux
      TYPE(section_vals_type), POINTER                   :: input, xc_section, xc_section_aux

      CALL timeset(routineN, handle)

      NULLIFY (auxbas_pw_pool, pw_env, v_rspace_new, para_env, rho1_r, &
               v_xc, rho1_ao, rho_ao, poisson_env, input, rho, dft_control, &
               logger, rho1_g)
      logger => cp_get_default_logger()

      energy_hartree = 0.0_dp
      energy_hartree_1c = 0.0_dp

      CPASSERT(ASSOCIATED(c0))
      CPASSERT(ASSOCIATED(v))
      CPASSERT(ASSOCIATED(Av))

      CPASSERT(ASSOCIATED(p_env%kpp1))
      CPASSERT(ASSOCIATED(p_env%kpp1_env))
      kpp1_env => p_env%kpp1_env
      CPASSERT(kpp1_env%ref_count > 0)

      CALL get_qs_env(qs_env=qs_env, &
                      ks_env=ks_env, &
                      pw_env=pw_env, &
                      input=input, &
                      admm_env=admm_env, &
                      para_env=para_env, &
                      rho=rho, &
                      linres_control=linres_control, &
                      dft_control=dft_control)

      gapw = dft_control%qs_control%gapw
      gapw_xc = dft_control%qs_control%gapw_xc
      lr_triplet = linres_control%lr_triplet

      rho1 => p_env%rho1
      rho1_xc => p_env%rho1_xc
      CPASSERT(ASSOCIATED(rho1))
      IF (gapw_xc) THEN
         CPASSERT(ASSOCIATED(rho1_xc))
      END IF

      CALL qs_rho_get(rho, rho_ao=rho_ao, rho_r=rho_r)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      nspins = SIZE(p_env%kpp1)
      lrigpw = dft_control%qs_control%lrigpw
      IF (lrigpw) THEN
         CALL get_qs_env(qs_env, &
                         lri_env=lri_env, &
                         lri_density=lri_density, &
                         atomic_kind_set=atomic_kind_set)
      ENDIF

      IF (.NOT. ASSOCIATED(kpp1_env%v_rspace)) THEN
         ALLOCATE (kpp1_env%v_rspace(nspins))
         DO ispin = 1, nspins
            CALL pw_pool_create_pw(auxbas_pw_pool, &
                                   kpp1_env%v_rspace(ispin)%pw, &
                                   use_data=REALDATA3D, in_space=REALSPACE)
         END DO
      END IF

      IF (.NOT. ASSOCIATED(kpp1_env%v_ao)) THEN
         CALL get_qs_env(qs_env, matrix_s=matrix_s)
         CALL dbcsr_allocate_matrix_set(kpp1_env%v_ao, nspins)
         DO ispin = 1, nspins
            ALLOCATE (kpp1_env%v_ao(ispin)%matrix)
            CALL dbcsr_copy(kpp1_env%v_ao(ispin)%matrix, matrix_s(1)%matrix, &
                            name="kpp1%v_ao-"//ADJUSTL(cp_to_string(ispin)))
         END DO
      END IF

      IF (dft_control%do_admm) THEN
         xc_section => admm_env%xc_section_primary
      ELSE
         xc_section => section_vals_get_subs_vals(input, "DFT%XC")
      END IF

      kpp1_env%iter = kpp1_env%iter + 1

      ! gets the tmp grids
      CPASSERT(ASSOCIATED(pw_env))
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                      poisson_env=poisson_env)
      ALLOCATE (v_rspace_new(nspins))
      CALL pw_pool_create_pw(auxbas_pw_pool, v_hartree_gspace%pw, &
                             use_data=COMPLEXDATA1D, &
                             in_space=RECIPROCALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, v_hartree_rspace%pw, &
                             use_data=REALDATA3D, &
                             in_space=REALSPACE)

      IF (gapw .OR. gapw_xc) &
         CALL prepare_gapw_den(qs_env, p_env%local_rho_set, do_rho0=(.NOT. gapw_xc))

      ! *** calculate the hartree potential on the total density ***
      CALL pw_pool_create_pw(auxbas_pw_pool, rho1_tot_gspace%pw, &
                             use_data=COMPLEXDATA1D, &
                             in_space=RECIPROCALSPACE)

      CALL qs_rho_get(rho1, rho_g=rho1_g)
      CALL pw_copy(rho1_g(1)%pw, rho1_tot_gspace%pw)
      DO ispin = 2, nspins
         CALL pw_axpy(rho1_g(ispin)%pw, rho1_tot_gspace%pw)
      END DO
      IF (gapw) &
         CALL pw_axpy(p_env%local_rho_set%rho0_mpole%rho0_s_gs%pw, rho1_tot_gspace%pw)

      IF (.NOT. (nspins == 1 .AND. lr_triplet)) THEN
         CALL pw_poisson_solve(poisson_env, rho1_tot_gspace%pw, &
                               energy_hartree, &
                               v_hartree_gspace%pw)
         CALL pw_transfer(v_hartree_gspace%pw, v_hartree_rspace%pw)
      ENDIF

      CALL pw_pool_give_back_pw(auxbas_pw_pool, rho1_tot_gspace%pw)

      ! *** calculate the xc potential ***
      NULLIFY (rho1a)
      IF (gapw_xc) THEN
         rho1a => rho1_xc
      ELSE
         rho1a => rho1
      END IF

      deriv2_analytic = section_get_lval(xc_section, "2ND_DERIV_ANALYTICAL")
      IF (deriv2_analytic) THEN
         CALL qs_rho_get(rho1a, rho_r=rho1_r)
         CALL qs_fxc_analytic(rho, rho1_r, xc_section, auxbas_pw_pool, lr_triplet, v_xc)
         IF (gapw) CALL calculate_xc_2nd_deriv_atom(p_env%local_rho_set, qs_env, xc_section, para_env, &
                                                    do_tddft=.FALSE., do_triplet=lr_triplet)
      ELSE
         NULLIFY (tau)
         CALL qs_fxc_fdiff(ks_env, rho, rho1a, xc_section, 6, lr_triplet, v_xc, tau)
         CPASSERT(.NOT. gapw)
      END IF

      DO ispin = 1, nspins
         v_rspace_new(ispin)%pw => v_xc(ispin)%pw
      END DO
      DEALLOCATE (v_xc)

      ! ADMM Correction
      IF (dft_control%do_admm) THEN
         IF (dft_control%admm_control%aux_exch_func /= do_admm_aux_exch_func_none) THEN
            IF (.NOT. ASSOCIATED(kpp1_env%deriv_set_admm)) THEN
               CPASSERT(.NOT. lr_triplet)
               xc_section_aux => admm_env%xc_section_aux
               CALL get_qs_env(qs_env=qs_env, rho_aux_fit=rho_aux)
               CALL qs_rho_get(rho_aux, rho_r=rho_r)
               CALL xc_prep_2nd_deriv(kpp1_env%deriv_set_admm, kpp1_env%rho_set_admm, &
                                      rho_r, auxbas_pw_pool, &
                                      xc_section=xc_section_aux)
            END IF
         END IF
      END IF

      !-------------------------------!
      ! Add both hartree and xc terms !
      !-------------------------------!
      DO ispin = 1, nspins

         IF (gapw_xc) THEN
            ! XC and Hartree are integrated separatedly
            ! XC uses the sofft basis set only
            v_rspace_new(ispin)%pw%cr3d = v_rspace_new(ispin)%pw%cr3d* &
                                          v_rspace_new(ispin)%pw%pw_grid%dvol

            IF (nspins == 1) THEN

               IF (.NOT. (lr_triplet)) THEN
                  v_rspace_new(1)%pw%cr3d = 2.0_dp*v_rspace_new(1)%pw%cr3d
               END IF
               CALL qs_rho_get(rho1, rho_ao=rho1_ao)
               ! remove kpp1_env%v_ao and work directly on k_p_p1 ?
               CALL dbcsr_set(kpp1_env%v_ao(ispin)%matrix, 0.0_dp)
               CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin), &
                                       pmat=rho1_ao(ispin), &
                                       hmat=kpp1_env%v_ao(ispin), &
                                       qs_env=qs_env, &
                                       calculate_forces=.FALSE., gapw=gapw_xc)

               ! add hartree only for SINGLETS
               IF (.NOT. lr_triplet) THEN
                  v_hartree_rspace%pw%cr3d = v_hartree_rspace%pw%cr3d* &
                                             v_hartree_rspace%pw%pw_grid%dvol
                  v_rspace_new(1)%pw%cr3d = 2.0_dp*v_hartree_rspace%pw%cr3d

                  CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin), &
                                          pmat=rho_ao(ispin), &
                                          hmat=kpp1_env%v_ao(ispin), &
                                          qs_env=qs_env, &
                                          calculate_forces=.FALSE., gapw=gapw)
               END IF
            ELSE
               ! remove kpp1_env%v_ao and work directly on k_p_p1 ?
               CALL dbcsr_set(kpp1_env%v_ao(ispin)%matrix, 0.0_dp)
               CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin), &
                                       pmat=rho_ao(ispin), &
                                       hmat=kpp1_env%v_ao(ispin), &
                                       qs_env=qs_env, &
                                       calculate_forces=.FALSE., gapw=gapw_xc)

               IF (ispin == 1) THEN
                  v_hartree_rspace%pw%cr3d = v_hartree_rspace%pw%cr3d* &
                                             v_hartree_rspace%pw%pw_grid%dvol
               END IF
               v_rspace_new(ispin)%pw%cr3d = v_hartree_rspace%pw%cr3d
               CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin), &
                                       pmat=rho_ao(ispin), &
                                       hmat=kpp1_env%v_ao(ispin), &
                                       qs_env=qs_env, &
                                       calculate_forces=.FALSE., gapw=gapw)
            END IF

         ELSE

            v_rspace_new(ispin)%pw%cr3d = v_rspace_new(ispin)%pw%cr3d* &
                                          v_rspace_new(ispin)%pw%pw_grid%dvol

            IF (nspins == 1) THEN
               IF (.NOT. (lr_triplet)) THEN
                  v_rspace_new(1)%pw%cr3d = 2.0_dp*v_rspace_new(1)%pw%cr3d
               ENDIF
               ! add hartree only for SINGLETS
               !IF (res_etype == tddfpt_singlet) THEN
               IF (.NOT. lr_triplet) THEN
                  v_hartree_rspace%pw%cr3d = v_hartree_rspace%pw%cr3d* &
                                             v_hartree_rspace%pw%pw_grid%dvol
                  v_rspace_new(1)%pw%cr3d = v_rspace_new(1)%pw%cr3d+ &
                                            2.0_dp*v_hartree_rspace%pw%cr3d
               END IF
            ELSE
               IF (ispin == 1) THEN
                  v_hartree_rspace%pw%cr3d = v_hartree_rspace%pw%cr3d* &
                                             v_hartree_rspace%pw%pw_grid%dvol
               END IF
               v_rspace_new(ispin)%pw%cr3d = v_rspace_new(ispin)%pw%cr3d+ &
                                             v_hartree_rspace%pw%cr3d
            END IF

            ! remove kpp1_env%v_ao and work directly on k_p_p1 ?
            CALL dbcsr_set(kpp1_env%v_ao(ispin)%matrix, 0.0_dp)

            IF (lrigpw) THEN
               lri_v_int => lri_density%lri_coefs(ispin)%lri_kinds
               CALL get_qs_env(qs_env, nkind=nkind)
               DO ikind = 1, nkind
                  lri_v_int(ikind)%v_int = 0.0_dp
               END DO
               CALL integrate_v_rspace_one_center(v_rspace_new(ispin), qs_env, &
                                                  lri_v_int, .FALSE., "LRI_AUX")
               DO ikind = 1, nkind
                  CALL mp_sum(lri_v_int(ikind)%v_int, para_env%group)
               END DO
               ALLOCATE (k1mat(1))
               k1mat(1)%matrix => kpp1_env%v_ao(ispin)%matrix
               IF (lri_env%exact_1c_terms) THEN
                  CALL integrate_v_rspace_diagonal(v_rspace_new(ispin), k1mat(1)%matrix, &
                                                   rho_ao(ispin)%matrix, qs_env, .FALSE., "ORB")
               END IF
               CALL calculate_lri_ks_matrix(lri_env, lri_v_int, k1mat, atomic_kind_set)
               DEALLOCATE (k1mat)
            ELSE
               CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin), &
                                       pmat=rho_ao(ispin), &
                                       hmat=kpp1_env%v_ao(ispin), &
                                       qs_env=qs_env, &
                                       calculate_forces=.FALSE., gapw=gapw)
            END IF

         END IF

         CALL dbcsr_copy(p_env%kpp1(ispin)%matrix, kpp1_env%v_ao(ispin)%matrix)
      END DO

      IF (gapw) THEN

         IF (.NOT. ((nspins == 1 .AND. lr_triplet))) THEN
            CALL Vh_1c_gg_integrals(qs_env, energy_hartree_1c, &
                                    p_env%hartree_local%ecoul_1c, &
                                    p_env%local_rho_set, &
                                    para_env, tddft=.TRUE.)

            CALL integrate_vhg0_rspace(qs_env, v_hartree_rspace, para_env, &
                                       calculate_forces=.FALSE., &
                                       local_rho_set=p_env%local_rho_set)
         END IF

         ! ***  Add single atom contributions to the KS matrix ***
         ! remap pointer
         ns = SIZE(p_env%kpp1)
         ksmat(1:ns, 1:1) => p_env%kpp1(1:ns)
         ns = SIZE(rho_ao)
         psmat(1:ns, 1:1) => rho_ao(1:ns)

         CALL update_ks_atom(qs_env, ksmat, psmat, forces=.FALSE., tddft=.TRUE., &
                             rho_atom_external=p_env%local_rho_set%rho_atom_set)

      END IF

      ! KG embedding, contribution of kinetic energy functional to kernel
      IF (dft_control%qs_control%do_kg .AND. .NOT. (lr_triplet .OR. gapw .OR. gapw_xc)) THEN
         IF (qs_env%kg_env%tnadd_method == kg_tnadd_embed) THEN

            CALL qs_rho_get(rho1, rho_ao=rho1_ao)
            alpha = 1.0_dp

            ekin_mol = 0.0_dp
            CALL get_qs_env(qs_env, kg_env=kg_env)
            CALL kg_ekin_subset(qs_env=qs_env, &
                                ks_matrix=p_env%kpp1, &
                                ekin_mol=ekin_mol, &
                                calc_force=.FALSE., &
                                do_kernel=.TRUE., &
                                pmat_ext=rho1_ao)
         END IF
      END IF

      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_gspace%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_rspace%pw)
      DO ispin = 1, nspins
         CALL pw_pool_give_back_pw(auxbas_pw_pool, v_rspace_new(ispin)%pw)
      END DO
      DEALLOCATE (v_rspace_new)

      DO ispin = 1, nspins
         CALL cp_fm_get_info(c0(ispin)%matrix, ncol_global=ncol)
         CALL cp_dbcsr_sm_fm_multiply(p_env%kpp1(ispin)%matrix, &
                                      c0(ispin)%matrix, &
                                      Av(ispin)%matrix, &
                                      ncol=ncol, alpha=1.0_dp, beta=1.0_dp)
      ENDDO

      CALL timestop(handle)

   END SUBROUTINE apply_op_2_dft

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param p_env ...
!> \param c0 ...
!> \param v ...
!> \param Av ...
! **************************************************************************************************
   SUBROUTINE apply_op_2_xtb(qs_env, p_env, c0, v, Av)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_p_env_type), POINTER                       :: p_env
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: c0, v, Av

      CHARACTER(len=*), PARAMETER                        :: routineN = 'apply_op_2_xtb'

      INTEGER                                            :: atom_a, handle, iatom, ikind, is, ispin, &
                                                            na, natom, natorb, ncol, nkind, ns, &
                                                            nsgf, nspins
      INTEGER, DIMENSION(25)                             :: lao
      INTEGER, DIMENSION(5)                              :: occ
      LOGICAL                                            :: lr_triplet
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: mcharge, mcharge1
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: aocg, aocg1, charges, charges1
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ao
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_p, matrix_p1, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(linres_control_type), POINTER                 :: linres_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kpp1_env_type), POINTER                    :: kpp1_env
      TYPE(qs_rho_type), POINTER                         :: rho, rho1
      TYPE(xtb_atom_type), POINTER                       :: xtb_kind

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(c0))
      CPASSERT(ASSOCIATED(v))
      CPASSERT(ASSOCIATED(Av))
      MARK_USED(v)

      CPASSERT(ASSOCIATED(p_env%kpp1_env))
      CPASSERT(ASSOCIATED(p_env%kpp1))
      kpp1_env => p_env%kpp1_env

      rho1 => p_env%rho1
      CPASSERT(ASSOCIATED(rho1))

      CPASSERT(kpp1_env%ref_count > 0)

      CALL get_qs_env(qs_env=qs_env, &
                      pw_env=pw_env, &
                      para_env=para_env, &
                      rho=rho, &
                      linres_control=linres_control, &
                      dft_control=dft_control)

      CALL qs_rho_get(rho, rho_ao=rho_ao)

      lr_triplet = linres_control%lr_triplet
      CPASSERT(.NOT. lr_triplet)

      nspins = SIZE(p_env%kpp1)
      kpp1_env%iter = kpp1_env%iter + 1

      DO ispin = 1, nspins
         CALL dbcsr_set(p_env%kpp1(ispin)%matrix, 0.0_dp)
      ENDDO

      IF (dft_control%qs_control%xtb_control%coulomb_interaction) THEN
         ! Mulliken charges
         CALL get_qs_env(qs_env, particle_set=particle_set, matrix_s_kp=matrix_s)
         natom = SIZE(particle_set)
         CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
         CALL qs_rho_get(rho1, rho_ao_kp=matrix_p1)
         ALLOCATE (mcharge(natom), charges(natom, 5))
         ALLOCATE (mcharge1(natom), charges1(natom, 5))
         charges = 0.0_dp
         charges1 = 0.0_dp
         CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set)
         nkind = SIZE(atomic_kind_set)
         CALL get_qs_kind_set(qs_kind_set, maxsgf=nsgf)
         ALLOCATE (aocg(nsgf, natom))
         aocg = 0.0_dp
         ALLOCATE (aocg1(nsgf, natom))
         aocg1 = 0.0_dp
         CALL ao_charges(matrix_p, matrix_s, aocg, para_env)
         CALL ao_charges(matrix_p1, matrix_s, aocg1, para_env)
         DO ikind = 1, nkind
            CALL get_atomic_kind(atomic_kind_set(ikind), natom=na)
            CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_kind)
            CALL get_xtb_atom_param(xtb_kind, natorb=natorb, lao=lao, occupation=occ)
            DO iatom = 1, na
               atom_a = atomic_kind_set(ikind)%atom_list(iatom)
               charges(atom_a, :) = REAL(occ(:), KIND=dp)
               DO is = 1, natorb
                  ns = lao(is) + 1
                  charges(atom_a, ns) = charges(atom_a, ns) - aocg(is, atom_a)
                  charges1(atom_a, ns) = charges1(atom_a, ns) - aocg1(is, atom_a)
               END DO
            END DO
         END DO
         DEALLOCATE (aocg, aocg1)
         DO iatom = 1, natom
            mcharge(iatom) = SUM(charges(iatom, :))
            mcharge1(iatom) = SUM(charges1(iatom, :))
         END DO
         ! Coulomb Kernel
         CALL xtb_coulomb_hessian(qs_env, p_env%kpp1, charges1, mcharge1, mcharge)
         !
         DEALLOCATE (charges, mcharge, charges1, mcharge1)
      END IF

      DO ispin = 1, nspins
         CALL cp_fm_get_info(c0(ispin)%matrix, ncol_global=ncol)
         CALL cp_dbcsr_sm_fm_multiply(p_env%kpp1(ispin)%matrix, &
                                      c0(ispin)%matrix, &
                                      Av(ispin)%matrix, &
                                      ncol=ncol, alpha=1.0_dp, beta=1.0_dp)
      ENDDO

      CALL timestop(handle)

   END SUBROUTINE apply_op_2_xtb

! **************************************************************************************************
!> \brief Update action of TDDFPT operator on trial vectors by adding exact-exchange term.
!> \param qs_env ...
!> \param p_env ...
!> \param c0 ...
!> \param v ...
!> \param Av ...
!> \par History
!>    * 11.2019 adapted from tddfpt_apply_hfx
! **************************************************************************************************
   SUBROUTINE apply_hfx(qs_env, p_env, c0, v, Av)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_p_env_type), POINTER                       :: p_env
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: c0, v, Av

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'apply_hfx'

      INTEGER                                            :: handle, ispin, nao, nao_aux, ncol, nspins
      LOGICAL                                            :: do_hfx
      REAL(KIND=dp)                                      :: alpha
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(cp_fm_type), POINTER                          :: tc0, tv
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, rho1_ao, work
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(section_vals_type), POINTER                   :: hfx_section, input

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()

      CPASSERT(ASSOCIATED(c0))
      CPASSERT(ASSOCIATED(v))
      CPASSERT(ASSOCIATED(Av))

      CALL get_qs_env(qs_env=qs_env, &
                      input=input, &
                      matrix_s=matrix_s, &
                      dft_control=dft_control)
      nspins = dft_control%nspins

      hfx_section => section_vals_get_subs_vals(input, "DFT%XC%HF")
      CALL section_vals_get(hfx_section, explicit=do_hfx)

      IF (do_hfx) THEN

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
            CALL get_qs_env(qs_env, admm_env=admm_env)
            CPASSERT(ASSOCIATED(admm_env%A))
            CPASSERT(ASSOCIATED(admm_env%work_aux_orb))
            CPASSERT(ASSOCIATED(admm_env%work_aux_orb2))
            CALL cp_fm_get_info(admm_env%A, nrow_global=nao_aux)
            CALL get_qs_env(qs_env=qs_env, matrix_s_aux_fit=matrix_s)
            NULLIFY (work, rho1_ao)
            CALL dbcsr_allocate_matrix_set(work, nspins)
            CALL dbcsr_allocate_matrix_set(rho1_ao, nspins)
            DO ispin = 1, nspins
               ALLOCATE (work(ispin)%matrix, rho1_ao(ispin)%matrix)
               CALL dbcsr_create(work(ispin)%matrix, template=matrix_s(1)%matrix)
               CALL dbcsr_copy(work(ispin)%matrix, matrix_s(1)%matrix)
               CALL dbcsr_set(work(ispin)%matrix, 0.0_dp)
               CALL dbcsr_create(rho1_ao(ispin)%matrix, template=matrix_s(1)%matrix)
               CALL dbcsr_copy(rho1_ao(ispin)%matrix, matrix_s(1)%matrix)
               CALL dbcsr_set(rho1_ao(ispin)%matrix, 0.0_dp)
            END DO
            ! P1 -> AUX BASIS
            DO ispin = 1, nspins
               CALL cp_fm_get_info(c0(ispin)%matrix, nrow_global=nao, ncol_global=ncol)
               tv => admm_env%work_aux_orb
               tc0 => admm_env%work_aux_orb2
               CALL cp_gemm('N', 'N', nao_aux, ncol, nao, 1.0_dp, admm_env%A, &
                            v(ispin)%matrix, 0.0_dp, tv)
               CALL cp_gemm('N', 'N', nao_aux, ncol, nao, 1.0_dp, admm_env%A, &
                            c0(ispin)%matrix, 0.0_dp, tc0)
               CALL cp_dbcsr_plus_fm_fm_t(rho1_ao(ispin)%matrix, matrix_v=tv, matrix_g=tc0, &
                                          ncol=ncol, alpha=1.0_dp)
               CALL cp_dbcsr_plus_fm_fm_t(rho1_ao(ispin)%matrix, matrix_v=tc0, matrix_g=tv, &
                                          ncol=ncol, alpha=1.0_dp)
            ENDDO
         ELSE
            NULLIFY (work, rho1_ao)
            CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s)
            CALL dbcsr_allocate_matrix_set(work, nspins)
            DO ispin = 1, nspins
               ALLOCATE (work(ispin)%matrix)
               CALL dbcsr_create(work(ispin)%matrix, template=matrix_s(1)%matrix)
               CALL dbcsr_copy(work(ispin)%matrix, matrix_s(1)%matrix)
               CALL dbcsr_set(work(ispin)%matrix, 0.0_dp)
            END DO
            rho1_ao => p_env%p1
         END IF

         CALL hfx_matrix(work, rho1_ao, qs_env, hfx_section)

         alpha = 2.0_dp
         IF (nspins == 2) alpha = 1.0_dp
         IF (dft_control%do_admm) THEN
            DO ispin = 1, nspins
               CALL cp_fm_get_info(c0(ispin)%matrix, nrow_global=nao, ncol_global=ncol)
               CALL cp_gemm('N', 'N', nao_aux, ncol, nao, 1.0_dp, admm_env%A, &
                            c0(ispin)%matrix, 0.0_dp, tc0)
               CALL cp_dbcsr_sm_fm_multiply(work(ispin)%matrix, tc0, tv, &
                                            ncol=ncol, alpha=alpha, beta=0.0_dp)
               CALL cp_gemm('T', 'N', nao, ncol, nao_aux, 1.0_dp, admm_env%A, &
                            tv, 1.0_dp, Av(ispin)%matrix)
            END DO
            CALL dbcsr_deallocate_matrix_set(rho1_ao)
            CALL dbcsr_deallocate_matrix_set(work)
         ELSE
            DO ispin = 1, nspins
               CALL cp_fm_get_info(c0(ispin)%matrix, ncol_global=ncol)
               CALL cp_dbcsr_sm_fm_multiply(work(ispin)%matrix, c0(ispin)%matrix, Av(ispin)%matrix, &
                                            ncol=ncol, alpha=alpha, beta=1.0_dp)
            END DO
            CALL dbcsr_deallocate_matrix_set(work)
         END IF

      END IF

      CALL timestop(handle)

   END SUBROUTINE apply_hfx

! **************************************************************************************************
!> \brief Add the hfx contributions to the Hamiltonian
!>
!> \param matrix_ks ...
!> \param rho_ao ...
!> \param qs_env ...
!> \param hfx_sections ...
!> \note
!>     Simplified version of subroutine hfx_ks_matrix()
! **************************************************************************************************
   SUBROUTINE hfx_matrix(matrix_ks, rho_ao, qs_env, hfx_sections)
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, rho_ao
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: hfx_sections

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'hfx_matrix'

      INTEGER                                            :: handle, irep, ispin, mspin, n_rep_hf, &
                                                            nspins
      LOGICAL                                            :: distribute_fock_matrix, &
                                                            hfx_treat_lsd_in_core, &
                                                            s_mstruct_changed
      REAL(KIND=dp)                                      :: eh1
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks_kp, rho_ao_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(hfx_type), DIMENSION(:, :), POINTER           :: x_data

      CALL timeset(routineN, handle)

      NULLIFY (dft_control, para_env, matrix_ks_kp, rho_ao_kp)

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      para_env=para_env, &
                      s_mstruct_changed=s_mstruct_changed, &
                      x_data=x_data)

      CPASSERT(dft_control%nimages == 1)
      nspins = dft_control%nspins

      CALL section_vals_get(hfx_sections, n_repetition=n_rep_hf)
      CALL section_vals_val_get(hfx_sections, "TREAT_LSD_IN_CORE", l_val=hfx_treat_lsd_in_core, &
                                i_rep_section=1)

      CALL section_vals_get(hfx_sections, n_repetition=n_rep_hf)
      distribute_fock_matrix = .TRUE.

      mspin = 1
      IF (hfx_treat_lsd_in_core) mspin = nspins

      matrix_ks_kp(1:nspins, 1:1) => matrix_ks(1:nspins)
      rho_ao_kp(1:nspins, 1:1) => rho_ao(1:nspins)

      DO irep = 1, n_rep_hf
         eh1 = 0.0_dp

         IF (x_data(irep, 1)%do_hfx_ri) THEN
            IF (x_data(irep, 1)%ri_data%flavor == ri_mo) THEN
               CPABORT("NYI with RI_FLAVOR MO")
            END IF
            CALL hfx_ri_update_ks(qs_env, x_data(irep, 1)%ri_data, matrix_ks_kp, eh1, &
                                  rho_ao=rho_ao_kp, geometry_did_change=s_mstruct_changed, &
                                  nspins=nspins, hf_fraction=x_data(irep, 1)%general_parameter%fraction)

         ELSE

            DO ispin = 1, mspin
               CALL integrate_four_center(qs_env, x_data, matrix_ks_kp, eh1, rho_ao_kp, hfx_sections, para_env, &
                                          s_mstruct_changed, irep, distribute_fock_matrix, ispin=ispin)
            END DO

         END IF
      END DO

      CALL timestop(handle)

   END SUBROUTINE hfx_matrix

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param p_env ...
!> \param c0 ...
!> \param v ...
!> \param Av ...
! **************************************************************************************************
   SUBROUTINE apply_xc_admm(qs_env, p_env, c0, v, Av)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_p_env_type), POINTER                       :: p_env
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: c0, v, Av

      CHARACTER(len=*), PARAMETER                        :: routineN = 'apply_xc_admm'

      INTEGER                                            :: handle, ispin, nao, nao_aux, ncol, nspins
      INTEGER, DIMENSION(2, 3)                           :: bo
      LOGICAL                                            :: lsd
      REAL(KIND=dp)                                      :: fac
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(cp_fm_type), POINTER                          :: tc0, tc1
      TYPE(dbcsr_p_type)                                 :: xcmat
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(linres_control_type), POINTER                 :: linres_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho1_aux_g, rho1_aux_r, tau_pw, v_xc
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(section_vals_type), POINTER                   :: xc_fun_section, xc_section
      TYPE(task_list_type), POINTER                      :: task_list_aux_fit
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type), POINTER                     :: rho1_set

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(c0))
      CPASSERT(ASSOCIATED(v))
      CPASSERT(ASSOCIATED(Av))

      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, task_list_aux_fit=task_list_aux_fit)

      IF (dft_control%do_admm) THEN
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
            NULLIFY (tau_pw)
            ! calculate the xc potential
            lsd = (nspins == 2)
            CALL get_qs_env(qs_env=qs_env, matrix_s_aux_fit=matrix_s)
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

            CALL admm_aux_reponse_density(qs_env, c0, v, rho1_aux_r, rho1_aux_g)
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
            CALL xc_calc_2nd_deriv(v_xc, p_env%kpp1_env%deriv_set_admm, p_env%kpp1_env%rho_set_admm, &
                                   rho1_aux_r, rho1_aux_g, auxbas_pw_pool, gapw=.FALSE., &
                                   xc_section=xc_section)
            CALL xc_rho_set_release(rho1_set)

            tc0 => admm_env%work_aux_orb
            tc1 => admm_env%work_aux_orb2
            CALL cp_fm_get_info(admm_env%A, nrow_global=nao_aux)
            DO ispin = 1, nspins
               v_xc(ispin)%pw%cr3d = v_xc(ispin)%pw%cr3d*v_xc(ispin)%pw%pw_grid%dvol
               IF (nspins == 1) THEN
                  v_xc(ispin)%pw%cr3d = 2.0_dp*v_xc(ispin)%pw%cr3d
               END IF
               CALL dbcsr_copy(xcmat%matrix, matrix_s(1)%matrix)
               CALL dbcsr_set(xcmat%matrix, 0.0_dp)
               CALL integrate_v_rspace(v_rspace=v_xc(ispin), hmat=xcmat, qs_env=qs_env, &
                                       calculate_forces=.FALSE., basis_type="AUX_FIT", &
                                       task_list_external=task_list_aux_fit)
               CALL cp_fm_get_info(c0(ispin)%matrix, nrow_global=nao, ncol_global=ncol)
               CALL cp_gemm('N', 'N', nao_aux, ncol, nao, 1.0_dp, admm_env%A, &
                            c0(ispin)%matrix, 0.0_dp, tc0)
               CALL cp_dbcsr_sm_fm_multiply(xcmat%matrix, tc0, tc1, &
                                            ncol=ncol, alpha=1.0_dp, beta=0.0_dp)
               CALL cp_gemm('T', 'N', nao, ncol, nao_aux, 1.0_dp, admm_env%A, &
                            tc1, 1.0_dp, Av(ispin)%matrix)
            END DO

            DO ispin = 1, nspins
               CALL pw_pool_give_back_pw(auxbas_pw_pool, v_xc(ispin)%pw)
               CALL pw_pool_give_back_pw(auxbas_pw_pool, rho1_aux_r(ispin)%pw)
               CALL pw_pool_give_back_pw(auxbas_pw_pool, rho1_aux_g(ispin)%pw)
            END DO
            DEALLOCATE (v_xc, rho1_aux_r, rho1_aux_g)
            CALL dbcsr_deallocate_matrix(xcmat%matrix)

         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE apply_xc_admm

! **************************************************************************************************
!> \brief Calculate ADMM auxiliary response density
!> \param qs_env ...
!> \param c0 ...
!> \param c1 ...
!> \param rho1_aux_r ...
!> \param rho1_aux_g ...
! **************************************************************************************************
   SUBROUTINE admm_aux_reponse_density(qs_env, c0, c1, rho1_aux_r, rho1_aux_g)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: c0, c1
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho1_aux_r, rho1_aux_g

      CHARACTER(LEN=*), PARAMETER :: routineN = 'admm_aux_reponse_density'

      INTEGER                                            :: handle, ispin, nao, nao_aux, ncol, nspins
      REAL(KIND=dp)                                      :: tot_rho_aux
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(cp_fm_type), POINTER                          :: tc0, tc1
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, rho1_ao
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(task_list_type), POINTER                      :: task_list_aux_fit

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, &
                      ks_env=ks_env, &
                      matrix_s=matrix_s, &
                      dft_control=dft_control)
      nspins = dft_control%nspins

      CALL get_qs_env(qs_env, admm_env=admm_env, &
                      task_list_aux_fit=task_list_aux_fit)
      CPASSERT(ASSOCIATED(admm_env%A))
      CPASSERT(ASSOCIATED(admm_env%work_aux_orb))
      CPASSERT(ASSOCIATED(admm_env%work_aux_orb2))
      CALL cp_fm_get_info(admm_env%A, nrow_global=nao_aux)
      CALL get_qs_env(qs_env=qs_env, matrix_s_aux_fit=matrix_s)
      NULLIFY (rho1_ao)
      CALL dbcsr_allocate_matrix_set(rho1_ao, nspins)
      DO ispin = 1, nspins
         ALLOCATE (rho1_ao(ispin)%matrix)
         CALL dbcsr_create(rho1_ao(ispin)%matrix, template=matrix_s(1)%matrix)
         CALL dbcsr_copy(rho1_ao(ispin)%matrix, matrix_s(1)%matrix)
         CALL dbcsr_set(rho1_ao(ispin)%matrix, 0.0_dp)
      END DO
      ! P1 -> AUX BASIS
      DO ispin = 1, nspins
         CALL cp_fm_get_info(c0(ispin)%matrix, nrow_global=nao, ncol_global=ncol)
         tc0 => admm_env%work_aux_orb
         tc1 => admm_env%work_aux_orb2
         CALL cp_gemm('N', 'N', nao_aux, ncol, nao, 1.0_dp, admm_env%A, &
                      c1(ispin)%matrix, 0.0_dp, tc1)
         CALL cp_gemm('N', 'N', nao_aux, ncol, nao, 1.0_dp, admm_env%A, &
                      c0(ispin)%matrix, 0.0_dp, tc0)
         CALL cp_dbcsr_plus_fm_fm_t(rho1_ao(ispin)%matrix, matrix_v=tc1, matrix_g=tc0, &
                                    ncol=ncol, alpha=1.0_dp)
         CALL cp_dbcsr_plus_fm_fm_t(rho1_ao(ispin)%matrix, matrix_v=tc0, matrix_g=tc1, &
                                    ncol=ncol, alpha=1.0_dp)
      ENDDO

      DO ispin = 1, nspins
         CALL calculate_rho_elec(matrix_p=rho1_ao(ispin)%matrix, &
                                 rho=rho1_aux_r(ispin), rho_gspace=rho1_aux_g(ispin), &
                                 total_rho=tot_rho_aux, ks_env=ks_env, &
                                 basis_type="AUX_FIT", &
                                 task_list_external=task_list_aux_fit)
      END DO

      CALL dbcsr_deallocate_matrix_set(rho1_ao)

      CALL timestop(handle)

   END SUBROUTINE admm_aux_reponse_density

! **************************************************************************************************
!> \brief ...
!> \param kpp1_env ...
!> \param qs_env ...
!> \param lr_triplet ...
! **************************************************************************************************
   SUBROUTINE kpp1_check_i_alloc(kpp1_env, qs_env, lr_triplet)

      TYPE(qs_kpp1_env_type), POINTER                    :: kpp1_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: lr_triplet

      INTEGER                                            :: ispin, nspins
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: my_rho_r, rho_r
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: input, xc_section

      NULLIFY (pw_env, auxbas_pw_pool, matrix_s, rho, rho_r, input)

      CPASSERT(ASSOCIATED(kpp1_env))
      CPASSERT(kpp1_env%ref_count > 0)

      CALL get_qs_env(qs_env, pw_env=pw_env, &
                      matrix_s=matrix_s, input=input, rho=rho)

      CALL qs_rho_get(rho, rho_r=rho_r)
      nspins = SIZE(rho_r)

      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      IF (.NOT. ASSOCIATED(kpp1_env%v_rspace)) THEN
         ALLOCATE (kpp1_env%v_rspace(nspins))
         DO ispin = 1, nspins
            CALL pw_pool_create_pw(auxbas_pw_pool, &
                                   kpp1_env%v_rspace(ispin)%pw, &
                                   use_data=REALDATA3D, in_space=REALSPACE)
         END DO
      END IF

      IF (.NOT. ASSOCIATED(kpp1_env%v_ao)) THEN
         CALL dbcsr_allocate_matrix_set(kpp1_env%v_ao, nspins)
         DO ispin = 1, nspins
            ALLOCATE (kpp1_env%v_ao(ispin)%matrix)
            CALL dbcsr_copy(kpp1_env%v_ao(ispin)%matrix, matrix_s(1)%matrix, &
                            name="kpp1%v_ao-"//ADJUSTL(cp_to_string(ispin)))
         END DO
      END IF

      CALL get_qs_env(qs_env, admm_env=admm_env, dft_control=dft_control)

      IF (.NOT. ASSOCIATED(kpp1_env%deriv_set)) THEN
         IF (nspins == 1 .AND. lr_triplet) THEN
            ALLOCATE (my_rho_r(2))
            DO ispin = 1, 2
               CALL pw_pool_create_pw(auxbas_pw_pool, my_rho_r(ispin)%pw, &
                                      use_data=rho_r(1)%pw%in_use, in_space=rho_r(1)%pw%in_space)
               my_rho_r(ispin)%pw%cr3d = 0.5_dp*rho_r(1)%pw%cr3d
            END DO
         ELSE
            ALLOCATE (my_rho_r(SIZE(rho_r)))
            DO ispin = 1, SIZE(rho_r)
               my_rho_r(ispin)%pw => rho_r(ispin)%pw
               CALL pw_retain(my_rho_r(ispin)%pw)
            END DO
         END IF

         IF (dft_control%do_admm) THEN
            xc_section => admm_env%xc_section_primary
         ELSE
            xc_section => section_vals_get_subs_vals(input, "DFT%XC")
         END IF

         CALL xc_prep_2nd_deriv(kpp1_env%deriv_set, kpp1_env%rho_set, &
                                my_rho_r, auxbas_pw_pool, &
                                xc_section=xc_section)

         DO ispin = 1, SIZE(my_rho_r)
            CALL pw_release(my_rho_r(ispin)%pw)
         ENDDO
         DEALLOCATE (my_rho_r)
      ENDIF

      ! ADMM Correction
      IF (dft_control%do_admm) THEN
         IF (dft_control%admm_control%aux_exch_func /= do_admm_aux_exch_func_none) THEN
            IF (.NOT. ASSOCIATED(kpp1_env%deriv_set_admm)) THEN
               CPASSERT(.NOT. lr_triplet)
               xc_section => admm_env%xc_section_aux
               CALL get_qs_env(qs_env=qs_env, rho_aux_fit=rho)
               CALL qs_rho_get(rho, rho_r=rho_r)
               CALL xc_prep_2nd_deriv(kpp1_env%deriv_set_admm, kpp1_env%rho_set_admm, &
                                      rho_r, auxbas_pw_pool, &
                                      xc_section=xc_section)
            END IF
         END IF
      END IF

   END SUBROUTINE kpp1_check_i_alloc

END MODULE qs_linres_kernel
