!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utility subroutine for qs energy calculation
!> \par History
!>      none
!> \author MK (29.10.2002)
! **************************************************************************************************
MODULE qs_matrix_w
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_gemm
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_release,&
                                              cp_cfm_to_fm,&
                                              cp_cfm_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type,&
                                              dbcsr_set
   USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                              fm_pool_get_el_struct
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE dft_plus_u,                      ONLY: lowdin_kp_apply_mo,&
                                              plus_u
   USE input_constants,                 ONLY: plus_u_lowdin
   USE kinds,                           ONLY: dp
   USE kpoint_methods,                  ONLY: kpoint_density_matrices,&
                                              kpoint_density_transform,&
                                              kpoint_ot_energy_weighted_density
   USE kpoint_types,                    ONLY: kpoint_env_type,&
                                              kpoint_type
   USE qs_density_matrices,             ONLY: calculate_w_matrix,&
                                              calculate_w_matrix_ot
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kpoint_operators,             ONLY: kpoint_operator_context_create,&
                                              kpoint_operator_context_release,&
                                              kpoint_operator_context_type,&
                                              kpoint_operator_get,&
                                              kpoint_operator_prepare,&
                                              kpoint_slot_ks
   USE qs_matrix_pools,                 ONLY: mpools_get
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE scf_control_types,               ONLY: scf_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_matrix_w'

   PUBLIC :: compute_matrix_w

CONTAINS

! **************************************************************************************************
!> \brief Refactoring of qs_energies_scf. Moves computation of matrix_w
!>        into separate subroutine
!> \param qs_env ...
!> \param calc_forces ...
!> \param include_lowdin_forces include the Lowdin overlap derivative before the W(k) transform
!> \par History
!>      05.2013 created [Florian Schiffmann]
! **************************************************************************************************

   SUBROUTINE compute_matrix_w(qs_env, calc_forces, include_lowdin_forces)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: calc_forces
      LOGICAL, INTENT(IN), OPTIONAL                      :: include_lowdin_forces

      CHARACTER(len=*), PARAMETER                        :: routineN = 'compute_matrix_w'

      INTEGER                                            :: handle, is, ispin, nao, nspin
      LOGICAL                                            :: do_kpoints, has_unit_metric
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s, matrix_w, &
                                                            mo_derivs, rho_ao
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mo_set_type), POINTER                         :: mo_set
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(scf_control_type), POINTER                    :: scf_control

      CALL timeset(routineN, handle)

      ! if calculate forces, time to compute the w matrix
      CALL get_qs_env(qs_env, has_unit_metric=has_unit_metric)

      IF (calc_forces .AND. .NOT. has_unit_metric) THEN
         CALL get_qs_env(qs_env, do_kpoints=do_kpoints)

         IF (do_kpoints) THEN
            BLOCK
               TYPE(cp_fm_type), DIMENSION(2)                   :: fmwork
               TYPE(cp_fm_struct_type), POINTER                   :: ao_ao_fmstruct
               TYPE(cp_fm_type), POINTER                          :: mo_coeff
               TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks_kp, matrix_s_kp, &
                                                                     matrix_w_kp
               TYPE(kpoint_type), POINTER                         :: kpoints
               TYPE(neighbor_list_set_p_type), DIMENSION(:), &
                  POINTER                                         :: sab_nl

               CALL get_qs_env(qs_env, &
                               matrix_w_kp=matrix_w_kp, &
                               matrix_ks_kp=matrix_ks_kp, &
                               matrix_s_kp=matrix_s_kp, &
                               sab_orb=sab_nl, &
                               mos=mos, &
                               dft_control=dft_control, &
                               kpoints=kpoints, &
                               scf_control=scf_control)

               CALL get_mo_set(mos(1), mo_coeff=mo_coeff, nao=nao)
               CALL cp_fm_struct_create(fmstruct=ao_ao_fmstruct, nrow_global=nao, ncol_global=nao, &
                                        template_fmstruct=mo_coeff%matrix_struct)

               DO is = 1, SIZE(fmwork)
                  CALL cp_fm_create(fmwork(is), matrix_struct=ao_ao_fmstruct)
               END DO
               CALL cp_fm_struct_release(ao_ao_fmstruct)

               ! energy weighted density matrices in k-space
               IF (scf_control%use_ot) THEN
                  IF (ALLOCATED(kpoints%lowdin_v)) CALL plus_u(qs_env, matrix_h=matrix_ks_kp)
                  CALL kpoint_ot_energy_weighted_matrices(kpoints, matrix_ks_kp, fmwork(1))
               ELSE
                  CALL kpoint_density_matrices(kpoints, energy_weighted=.TRUE.)
               END IF
               ! Property callers need the ordinary energy-weighted density. Force callers
               ! also need -dE_U/dS; transform both together, once, on the overlap support.
               IF (PRESENT(include_lowdin_forces)) THEN
                  IF (include_lowdin_forces .AND. dft_control%dft_plus_u .AND. &
                      dft_control%plus_u_method_id == plus_u_lowdin) THEN
                     CALL plus_u(qs_env, matrix_w=matrix_w_kp)
                  END IF
               END IF
               ! energy weighted density matrices in real space
               CALL kpoint_density_transform(kpoints, matrix_w_kp, .TRUE., &
                                             matrix_s_kp(1, 1)%matrix, sab_nl, fmwork)

               DO is = 1, SIZE(fmwork)
                  CALL cp_fm_release(fmwork(is))
               END DO

            END BLOCK
         ELSE

            NULLIFY (dft_control, rho_ao)
            CALL get_qs_env(qs_env, &
                            matrix_w=matrix_w, &
                            matrix_ks=matrix_ks, &
                            matrix_s=matrix_s, &
                            mo_derivs=mo_derivs, &
                            scf_control=scf_control, &
                            mos=mos, &
                            rho=rho, &
                            dft_control=dft_control)

            CALL qs_rho_get(rho, rho_ao=rho_ao)

            nspin = SIZE(mos)
            DO ispin = 1, nspin
               mo_set => mos(ispin)
               IF (dft_control%roks) THEN
                  IF (scf_control%use_ot) THEN
                     IF (ispin > 1) THEN
                        ! not very elegant, indeed ...
                        CALL dbcsr_set(matrix_w(ispin)%matrix, 0.0_dp)
                     ELSE
                        CALL calculate_w_matrix_ot(mo_set, mo_derivs(ispin)%matrix, &
                                                   matrix_w(ispin)%matrix, matrix_s(1)%matrix)
                     END IF
                  ELSE
                     CALL calculate_w_matrix(mo_set=mo_set, &
                                             matrix_ks=matrix_ks(ispin)%matrix, &
                                             matrix_p=rho_ao(ispin)%matrix, &
                                             matrix_w=matrix_w(ispin)%matrix)
                  END IF
               ELSE
                  IF (scf_control%use_ot) THEN
                     CALL calculate_w_matrix_ot(mo_set, mo_derivs(ispin)%matrix, &
                                                matrix_w(ispin)%matrix, matrix_s(1)%matrix)
                  ELSE
                     CALL calculate_w_matrix(mo_set, matrix_w(ispin)%matrix)
                  END IF
               END IF
            END DO

         END IF

      END IF

      CALL timestop(handle)

   END SUBROUTINE compute_matrix_w

! **************************************************************************************************
!> \brief Build OT energy-weighted densities using Hamiltonians on the k-point groups.
!> \param kpoint K-point environment
!> \param matrix_ks current real-space Kohn-Sham matrices
!> \param ao_ao_fm square full-matrix template on the global parallel environment
! **************************************************************************************************
   SUBROUTINE kpoint_ot_energy_weighted_matrices(kpoint, matrix_ks, ao_ao_fm)

      TYPE(kpoint_type), POINTER                         :: kpoint
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks
      TYPE(cp_fm_type), INTENT(IN)                       :: ao_ao_fm

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_ot_energy_weighted_matrices'
      COMPLEX(KIND=dp), PARAMETER                        :: zone = (1.0_dp, 0.0_dp), &
                                                            zzero = (0.0_dp, 0.0_dp)

      INTEGER                                            :: handle, ikpgr, ispin, kplocal, kpmax, &
                                                            nao, nmo, nspin
      REAL(KIND=dp), DIMENSION(:), POINTER               :: occupation
      TYPE(cp_cfm_type)                                  :: h_k, hc, wmat
      TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER     :: ao_ao_fm_pools
      TYPE(cp_fm_struct_type), POINTER                   :: h_k_struct
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(kpoint_operator_context_type)                 :: op_ctx

      CALL timeset(routineN, handle)
      NULLIFY (ao_ao_fm_pools, h_k_struct, occupation)
      CPASSERT(.NOT. kpoint%use_real_wfn)
      CPASSERT(ASSOCIATED(matrix_ks))

      nspin = SIZE(matrix_ks, 1)
      kplocal = SIZE(kpoint%kp_env)
      kpmax = MAXVAL(kpoint%kp_dist(2, :) - kpoint%kp_dist(1, :) + 1)
      CALL kpoint_operator_context_create(op_ctx, kpoint, ao_ao_fm, nspin)
      ! Capture the current Hamiltonian after SCF, independently of OT's iteration caches.
      CALL kpoint_operator_prepare(op_ctx, matrix_ks, kpoint_slot_ks)
      CALL mpools_get(kpoint%mpools, ao_ao_fm_pools=ao_ao_fm_pools)
      h_k_struct => fm_pool_get_el_struct(ao_ao_fm_pools(1)%pool)
      CALL cp_cfm_create(h_k, h_k_struct)

      DO ikpgr = 1, kpmax
         DO ispin = 1, nspin
            ! Every group calls the collective service; the AO/AO pool layout keeps H(k)
            ! valid even on inactive rounds, when fallback assembly may communicate globally.
            CALL kpoint_operator_get(op_ctx, ikpgr, ispin, matrix_ks, &
                                     cmat=h_k, &
                                     matrix_row=ispin, slot=kpoint_slot_ks)
            IF (ikpgr <= kplocal) THEN
               kp => kpoint%kp_env(ikpgr)%kpoint_env
               CALL get_mo_set(kp%mos(ispin), nmo=nmo, &
                               occupation_numbers=occupation)
               CALL cp_cfm_get_info(kp%mos(ispin)%cmo_coeff, nrow_global=nao)
               CPASSERT(nmo >= 1)
               CPASSERT(SIZE(occupation) >= nmo)

               ! Smearing occupations are k-point resolved, whereas mo_set%homo is not.
               CALL cp_cfm_create(hc, kp%mos(ispin)%cmo_coeff%matrix_struct, set_zero=.TRUE.)
               CALL cp_cfm_gemm("N", "N", nao, nmo, nao, zone, h_k, &
                                kp%mos(ispin)%cmo_coeff, zzero, hc)
               IF (ALLOCATED(kpoint%lowdin_v)) THEN
                  CALL lowdin_kp_apply_mo(kp, kpoint%lowdin_v(:, ispin), &
                                          kp%mos(ispin)%cmo_coeff, hc)
               END IF

               CALL cp_cfm_create(wmat, kp%wmat(1, ispin)%matrix_struct)
               CALL kpoint_ot_energy_weighted_density( &
                  kp%mos(ispin)%cmo_coeff, hc, occupation, wmat)
               ! wmat remains split-FM storage; keep that conversion at its boundary.
               CALL cp_cfm_to_fm(wmat, kp%wmat(1, ispin), kp%wmat(2, ispin))

               CALL cp_cfm_release(wmat)
               CALL cp_cfm_release(hc)
            END IF
         END DO
      END DO

      CALL cp_cfm_release(h_k)
      CALL kpoint_operator_context_release(op_ctx)
      CALL timestop(handle)

   END SUBROUTINE kpoint_ot_energy_weighted_matrices

END MODULE qs_matrix_w
