!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Methods for preparing and committing k-point orbital states to the QS environment.
! **************************************************************************************************
MODULE qs_kpoint_state
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_gemm,&
                                              cp_cfm_scale
   USE cp_cfm_diag,                     ONLY: cp_cfm_geeig,&
                                              cp_cfm_heevd
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_release,&
                                              cp_cfm_to_cfm,&
                                              cp_cfm_to_fm,&
                                              cp_cfm_type,&
                                              cp_fm_to_cfm
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type,&
                                              cp_fm_vectorsnorm
   USE kinds,                           ONLY: dp
   USE kpoint_methods,                  ONLY: kpoint_density_matrices,&
                                              kpoint_density_transform,&
                                              kpoint_set_mo_occupation
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_env_type,&
                                              kpoint_type
   USE mathconstants,                   ONLY: z_one,&
                                              z_zero
   USE message_passing,                 ONLY: mp_para_env_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_ks_types,                     ONLY: qs_ks_did_change
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_rho_methods,                  ONLY: qs_rho_update_rho
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE scf_control_types,               ONLY: scf_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_kpoint_state'

   PUBLIC :: qs_kpoint_copy_spin_mos, qs_kpoint_mos_initialized, &
             qs_kpoint_state_canonicalize_fixed, qs_kpoint_state_commit, &
             qs_kpoint_state_prepare_fixed_density

CONTAINS

! **************************************************************************************************
!> \brief Reconstructs fixed-rank complex k-point orbitals from a density kernel.
!>        The density and overlap must first be transferred to kp%ot_hmat and kp%ot_smat. The
!>        largest natural-orbital subspace gives a metric-defined fixed-rank projector without
!>        replacing the guess by the ground-state projector of its initial Hamiltonian.
!> \param kpoints K-point environment containing the transferred density and overlap matrices.
! **************************************************************************************************
   SUBROUTINE qs_kpoint_state_prepare_fixed_density(kpoints)
      TYPE(kpoint_type), POINTER                         :: kpoints

      INTEGER                                            :: ikpoint_local, ispin, nao, nmo, nocc
      LOGICAL                                            :: use_real_wfn
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: density_eigenvalues
      TYPE(cp_cfm_type)                                  :: density_metric, density_work, matrix_p, &
                                                            matrix_s, matrix_s_work, projected
      TYPE(cp_fm_struct_type), POINTER                   :: active_struct, square_struct
      TYPE(cp_fm_type)                                   :: active_im, active_re
      TYPE(cp_fm_type), POINTER                          :: mo_coeff_im, mo_coeff_re
      TYPE(kpoint_env_type), POINTER                     :: kp

      NULLIFY (active_struct, kp, mo_coeff_im, mo_coeff_re, square_struct)

      CPASSERT(ASSOCIATED(kpoints))
      CPASSERT(ASSOCIATED(kpoints%kp_env))
      CALL get_kpoint_info(kpoints, use_real_wfn=use_real_wfn)
      IF (use_real_wfn) THEN
         CPABORT("Fixed-density k-point OT initialization requires complex wavefunctions.")
      END IF

      DO ikpoint_local = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(ikpoint_local)%kpoint_env
         CPASSERT(ASSOCIATED(kp))
         CPASSERT(ASSOCIATED(kp%mos))
         CPASSERT(ASSOCIATED(kp%ot_hmat))
         CPASSERT(ASSOCIATED(kp%ot_smat))
         CPASSERT(SIZE(kp%ot_hmat, 1) >= 2)
         CPASSERT(SIZE(kp%ot_hmat, 2) >= SIZE(kp%mos, 2))
         CPASSERT(SIZE(kp%ot_smat) >= 2)

         CALL cp_fm_get_info(kp%ot_hmat(1, 1), matrix_struct=square_struct, nrow_global=nao)
         CALL cp_cfm_create(matrix_p, square_struct, name="natural-orbital density P(k)")
         CALL cp_cfm_create(matrix_s, square_struct, name="natural-orbital overlap S(k)")
         CALL cp_fm_to_cfm(kp%ot_smat(1), kp%ot_smat(2), matrix_s)

         DO ispin = 1, SIZE(kp%mos, 2)
            CALL get_mo_set(kp%mos(1, ispin), homo=nocc, mo_coeff=mo_coeff_re, nmo=nmo)
            CALL get_mo_set(kp%mos(2, ispin), mo_coeff=mo_coeff_im)
            CPASSERT(ASSOCIATED(mo_coeff_re))
            CPASSERT(ASSOCIATED(mo_coeff_im))
            CPASSERT(nocc >= 1 .AND. nocc <= nmo)

            CALL cp_fm_struct_create(active_struct, template_fmstruct=mo_coeff_re%matrix_struct, &
                                     ncol_global=nocc)
            CALL cp_fm_to_cfm(kp%ot_hmat(1, ispin), kp%ot_hmat(2, ispin), matrix_p)
            CALL cp_cfm_create(projected, active_struct, name="fixed-density natural orbitals")
            CALL cp_cfm_create(density_metric, square_struct, name="natural-orbital metric")
            CALL cp_cfm_create(density_work, square_struct, name="natural-orbital work")
            CALL cp_cfm_create(matrix_s_work, square_struct, name="natural-orbital overlap work")

            ! The largest eigenvectors of S*P*S*C = S*C*n are the natural orbitals of
            ! S**(1/2)*P*S**(1/2). Negating the left-hand side lets the generalized
            ! eigensolver request only its lowest nocc eigenpairs when supported.
            CALL cp_cfm_gemm('N', 'N', nao, nao, nao, z_one, matrix_p, matrix_s, &
                             z_zero, density_work)
            CALL cp_cfm_gemm('N', 'N', nao, nao, nao, z_one, matrix_s, density_work, &
                             z_zero, density_metric)
            CALL cp_cfm_scale(-z_one, density_metric)
            CALL cp_cfm_to_cfm(matrix_s, matrix_s_work)
            ALLOCATE (density_eigenvalues(nocc))
            CALL cp_cfm_geeig(density_metric, matrix_s_work, projected, &
                              density_eigenvalues, density_work, lowest_subset=.TRUE.)
            DEALLOCATE (density_eigenvalues)

            CALL cp_fm_create(active_re, active_struct, name="fixed-density orbitals real")
            CALL cp_fm_create(active_im, active_struct, name="fixed-density orbitals imaginary")
            CALL cp_cfm_to_fm(projected, active_re, active_im)
            CALL cp_fm_set_all(mo_coeff_re, 0.0_dp)
            CALL cp_fm_set_all(mo_coeff_im, 0.0_dp)
            CALL cp_fm_to_fm(active_re, mo_coeff_re, nocc)
            CALL cp_fm_to_fm(active_im, mo_coeff_im, nocc)

            CALL cp_fm_release(active_im)
            CALL cp_fm_release(active_re)
            CALL cp_cfm_release(matrix_s_work)
            CALL cp_cfm_release(density_work)
            CALL cp_cfm_release(density_metric)
            CALL cp_cfm_release(projected)
            CALL cp_fm_struct_release(active_struct)
         END DO

         CALL cp_cfm_release(matrix_s)
         CALL cp_cfm_release(matrix_p)
      END DO

   END SUBROUTINE qs_kpoint_state_prepare_fixed_density

! **************************************************************************************************
!> \brief Rotates a uniformly occupied k-point subspace into Ritz states of the current
!>        Hamiltonian and stores their energy labels. The Hamiltonian must first be transferred
!>        to kp%ot_hmat. Since the rotation is internal to a uniformly occupied subspace, it does
!>        not change the density represented by the orbitals.
!> \param kpoints K-point environment containing the orbitals and transferred Hamiltonians.
! **************************************************************************************************
   SUBROUTINE qs_kpoint_state_canonicalize_fixed(kpoints)
      TYPE(kpoint_type), POINTER                         :: kpoints

      INTEGER                                            :: ikpoint_local, ispin, nao, nmo, nocc
      LOGICAL                                            :: use_real_wfn
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: ritz_values
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues_im, eigenvalues_re
      TYPE(cp_cfm_type)                                  :: coeff, eigenvectors, hblock, hc, &
                                                            matrix_h, rotated
      TYPE(cp_fm_struct_type), POINTER                   :: active_struct, square_struct
      TYPE(cp_fm_type)                                   :: active_im, active_re
      TYPE(cp_fm_type), POINTER                          :: mo_coeff_im, mo_coeff_re
      TYPE(kpoint_env_type), POINTER                     :: kp

      NULLIFY (active_struct, eigenvalues_im, eigenvalues_re, kp, mo_coeff_im, mo_coeff_re, &
               square_struct)

      CPASSERT(ASSOCIATED(kpoints))
      CPASSERT(ASSOCIATED(kpoints%kp_env))
      CALL get_kpoint_info(kpoints, use_real_wfn=use_real_wfn)
      IF (use_real_wfn) THEN
         CPABORT("Fixed-density k-point OT canonicalization requires complex wavefunctions.")
      END IF

      DO ikpoint_local = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(ikpoint_local)%kpoint_env
         CPASSERT(ASSOCIATED(kp))
         CPASSERT(ASSOCIATED(kp%mos))
         CPASSERT(ASSOCIATED(kp%ot_hmat))
         CPASSERT(SIZE(kp%ot_hmat, 1) >= 2)
         CPASSERT(SIZE(kp%ot_hmat, 2) >= SIZE(kp%mos, 2))

         CALL cp_fm_get_info(kp%ot_hmat(1, 1), matrix_struct=square_struct, nrow_global=nao)
         CALL cp_cfm_create(matrix_h, square_struct, name="fixed-density Ritz Hamiltonian")

         DO ispin = 1, SIZE(kp%mos, 2)
            CALL get_mo_set(kp%mos(1, ispin), homo=nocc, mo_coeff=mo_coeff_re, nmo=nmo, &
                            eigenvalues=eigenvalues_re)
            CALL get_mo_set(kp%mos(2, ispin), mo_coeff=mo_coeff_im, &
                            eigenvalues=eigenvalues_im)
            CPASSERT(ASSOCIATED(mo_coeff_re))
            CPASSERT(ASSOCIATED(mo_coeff_im))
            CPASSERT(ASSOCIATED(eigenvalues_re))
            CPASSERT(ASSOCIATED(eigenvalues_im))
            CPASSERT(nocc >= 1 .AND. nocc <= nmo)
            CPASSERT(SIZE(eigenvalues_re) >= nocc)
            CPASSERT(SIZE(eigenvalues_im) >= nocc)

            CALL cp_fm_struct_create(active_struct, template_fmstruct=mo_coeff_re%matrix_struct, &
                                     ncol_global=nocc)
            CALL cp_fm_create(active_re, active_struct, name="fixed-density Ritz orbitals real")
            CALL cp_fm_create(active_im, active_struct, name="fixed-density Ritz orbitals imaginary")
            CALL cp_fm_to_fm(mo_coeff_re, active_re, nocc)
            CALL cp_fm_to_fm(mo_coeff_im, active_im, nocc)

            CALL cp_cfm_create(coeff, active_struct, name="fixed-density Ritz orbitals")
            CALL cp_cfm_create(hc, active_struct, name="fixed-density Ritz HC")
            CALL cp_cfm_create(rotated, active_struct, name="fixed-density Ritz rotated orbitals")
            CALL cp_cfm_create(hblock, active_struct, nrow=nocc, ncol=nocc, &
                               name="fixed-density projected Hamiltonian")
            CALL cp_cfm_create(eigenvectors, active_struct, nrow=nocc, ncol=nocc, &
                               name="fixed-density Ritz eigenvectors")
            CALL cp_fm_to_cfm(kp%ot_hmat(1, ispin), kp%ot_hmat(2, ispin), matrix_h)
            CALL cp_fm_to_cfm(active_re, active_im, coeff)

            CALL cp_cfm_gemm('N', 'N', nao, nocc, nao, z_one, matrix_h, coeff, z_zero, hc)
            CALL cp_cfm_gemm('C', 'N', nocc, nocc, nao, z_one, coeff, hc, z_zero, hblock)
            ALLOCATE (ritz_values(nocc))
            CALL cp_cfm_heevd(hblock, eigenvectors, ritz_values)
            CALL cp_cfm_gemm('N', 'N', nao, nocc, nocc, z_one, coeff, eigenvectors, &
                             z_zero, rotated)
            CALL cp_cfm_to_fm(rotated, active_re, active_im)
            CALL cp_fm_to_fm(active_re, mo_coeff_re, nocc)
            CALL cp_fm_to_fm(active_im, mo_coeff_im, nocc)
            eigenvalues_re(1:nocc) = ritz_values
            eigenvalues_im(1:nocc) = ritz_values
            DEALLOCATE (ritz_values)

            CALL cp_cfm_release(eigenvectors)
            CALL cp_cfm_release(hblock)
            CALL cp_cfm_release(rotated)
            CALL cp_cfm_release(hc)
            CALL cp_cfm_release(coeff)
            CALL cp_fm_release(active_im)
            CALL cp_fm_release(active_re)
            CALL cp_fm_struct_release(active_struct)
         END DO

         CALL cp_cfm_release(matrix_h)
      END DO

   END SUBROUTINE qs_kpoint_state_canonicalize_fixed

! **************************************************************************************************
!> \brief Rebuilds the physical density from the current k-point orbitals and occupations.
!> \param qs_env QS environment containing the orbital state to commit.
!> \param update_occupations Recompute occupations from the current MO energies and SCF settings.
!> \param separate_spin_occupations Fill restricted energy channels independently.
!> \param fixed_occupations Restore the fixed occupied rank required by OT.
! **************************************************************************************************
   SUBROUTINE qs_kpoint_state_commit(qs_env, update_occupations, &
                                     separate_spin_occupations, fixed_occupations)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: update_occupations
      LOGICAL, INTENT(IN), OPTIONAL                      :: separate_spin_occupations, &
                                                            fixed_occupations

      INTEGER                                            :: ispin
      LOGICAL                                            :: fixed, separate_spins
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s_kp, rho_ao_kp
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mo_set_type), DIMENSION(:, :), POINTER        :: mos
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control

      NULLIFY (energy, kpoints, matrix_s_kp, mos, rho_ao_kp, sab_nl, rho, scf_env, scf_control)

      fixed = .FALSE.
      separate_spins = .FALSE.
      IF (PRESENT(fixed_occupations)) fixed = fixed_occupations
      IF (PRESENT(separate_spin_occupations)) separate_spins = separate_spin_occupations
      IF (fixed .AND. .NOT. update_occupations) THEN
         CPABORT("Fixed occupations can only be requested while updating occupations.")
      END IF

      CALL get_qs_env(qs_env, kpoints=kpoints, matrix_s_kp=matrix_s_kp, &
                      energy=energy, rho=rho, scf_env=scf_env, scf_control=scf_control)
      CALL get_kpoint_info(kpoints, sab_nl=sab_nl)

      IF (update_occupations) THEN
         CALL kpoint_set_mo_occupation( &
            kpoints, scf_control%smear, separate_spin_occupations=separate_spins)
         IF (fixed) CALL qs_kpoint_set_fixed_occupations(kpoints)
      END IF
      CALL kpoint_density_matrices(kpoints)
      CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)
      CALL kpoint_density_transform(kpoints, rho_ao_kp, .FALSE., &
                                    matrix_s_kp(1, 1)%matrix, sab_nl, scf_env%scf_work1, &
                                    overlap_rs=matrix_s_kp)
      CALL qs_rho_update_rho(rho, qs_env=qs_env)
      CALL qs_ks_did_change(qs_env%ks_env, rho_changed=.TRUE.)

      energy%kTS = 0.0_dp
      energy%efermi = 0.0_dp
      mos => kpoints%kp_env(1)%kpoint_env%mos
      CPASSERT(ASSOCIATED(mos))
      DO ispin = 1, SIZE(mos, 2)
         energy%kTS = energy%kTS + mos(1, ispin)%kTS
         energy%efermi = energy%efermi + mos(1, ispin)%mu
      END DO
      energy%efermi = energy%efermi/REAL(SIZE(mos, 2), KIND=dp)

   END SUBROUTINE qs_kpoint_state_commit

! **************************************************************************************************
!> \brief Checks whether every local k-point channel contains nonzero occupied orbitals.
!> \param kpoints K-point environment.
!> \param require_full_space Require every available orbital instead of the occupied block.
!> \return True if every channel contains a nonzero real or imaginary coefficient.
! **************************************************************************************************
   FUNCTION qs_kpoint_mos_initialized(kpoints, require_full_space) RESULT(initialized)
      TYPE(kpoint_type), POINTER                         :: kpoints
      LOGICAL, INTENT(IN), OPTIONAL                      :: require_full_space
      LOGICAL                                            :: initialized

      INTEGER                                            :: icomponent, ikpoint_local, ispin, &
                                                            nactive, nmo, nocc
      LOGICAL                                            :: full_space
      REAL(KIND=dp)                                      :: coefficient_norm
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: column_norm, component_norm
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(mp_para_env_type), POINTER                    :: para_env

      NULLIFY (kp, mo_coeff, para_env)
      CPASSERT(ASSOCIATED(kpoints))
      CPASSERT(ASSOCIATED(kpoints%kp_env))

      full_space = .FALSE.
      IF (PRESENT(require_full_space)) full_space = require_full_space
      coefficient_norm = HUGE(1.0_dp)
      DO ikpoint_local = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(ikpoint_local)%kpoint_env
         CPASSERT(ASSOCIATED(kp))
         CPASSERT(ASSOCIATED(kp%mos))
         DO ispin = 1, SIZE(kp%mos, 2)
            CALL get_mo_set(kp%mos(1, ispin), homo=nocc, mo_coeff=mo_coeff, nmo=nmo)
            CPASSERT(ASSOCIATED(mo_coeff))
            nactive = MERGE(nmo, nocc, full_space)
            IF (nactive < 1 .OR. nactive > nmo) THEN
               coefficient_norm = 0.0_dp
               CYCLE
            END IF
            ALLOCATE (column_norm(nmo), component_norm(nmo))
            column_norm(:) = 0.0_dp
            DO icomponent = 1, SIZE(kp%mos, 1)
               CALL get_mo_set(kp%mos(icomponent, ispin), mo_coeff=mo_coeff)
               CPASSERT(ASSOCIATED(mo_coeff))
               CALL cp_fm_vectorsnorm(mo_coeff, component_norm)
               column_norm(:) = column_norm(:) + component_norm(:)**2
            END DO
            coefficient_norm = MIN(coefficient_norm, SQRT(MINVAL(column_norm(1:nactive))))
            DEALLOCATE (column_norm, component_norm)
         END DO
      END DO

      CALL get_kpoint_info(kpoints, para_env=para_env)
      CPASSERT(ASSOCIATED(para_env))
      CALL para_env%min(coefficient_norm)
      initialized = coefficient_norm > 100.0_dp*EPSILON(1.0_dp)

   END FUNCTION qs_kpoint_mos_initialized

! **************************************************************************************************
!> \brief Copies the first k-point spin channel into the remaining channels.
!> \param kpoints K-point environment.
!> \param nspin Number of spin channels to populate.
! **************************************************************************************************
   SUBROUTINE qs_kpoint_copy_spin_mos(kpoints, nspin)
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER, INTENT(IN)                                :: nspin

      INTEGER                                            :: icomponent, ikpoint_local, ispin
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues_source, eigenvalues_target
      TYPE(cp_fm_type), POINTER                          :: mo_coeff_source, mo_coeff_target
      TYPE(kpoint_env_type), POINTER                     :: kp

      NULLIFY (eigenvalues_source, eigenvalues_target, kp, mo_coeff_source, mo_coeff_target)
      CPASSERT(ASSOCIATED(kpoints))
      CPASSERT(ASSOCIATED(kpoints%kp_env))
      CPASSERT(nspin >= 1)

      DO ikpoint_local = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(ikpoint_local)%kpoint_env
         CPASSERT(ASSOCIATED(kp))
         CPASSERT(ASSOCIATED(kp%mos))
         CPASSERT(SIZE(kp%mos, 2) >= nspin)
         DO icomponent = 1, SIZE(kp%mos, 1)
            CALL get_mo_set(kp%mos(icomponent, 1), mo_coeff=mo_coeff_source, &
                            eigenvalues=eigenvalues_source)
            CPASSERT(ASSOCIATED(mo_coeff_source))
            CPASSERT(ASSOCIATED(eigenvalues_source))
            DO ispin = 2, nspin
               CALL get_mo_set(kp%mos(icomponent, ispin), mo_coeff=mo_coeff_target, &
                               eigenvalues=eigenvalues_target)
               CPASSERT(ASSOCIATED(mo_coeff_target))
               CPASSERT(ASSOCIATED(eigenvalues_target))
               CALL cp_fm_to_fm(mo_coeff_source, mo_coeff_target)
               eigenvalues_target(:) = eigenvalues_source(:)
            END DO
         END DO
      END DO

   END SUBROUTINE qs_kpoint_copy_spin_mos

! **************************************************************************************************
!> \brief Restores the fixed occupied rank required by k-point OT.
!> \param kpoints K-point environment.
! **************************************************************************************************
   SUBROUTINE qs_kpoint_set_fixed_occupations(kpoints)
      TYPE(kpoint_type), POINTER                         :: kpoints

      INTEGER                                            :: icomponent, ikpoint_local, ispin, &
                                                            nelectron, nmo, nocc
      REAL(KIND=dp)                                      :: maxocc
      REAL(KIND=dp), DIMENSION(:), POINTER               :: occupation_numbers
      TYPE(kpoint_env_type), POINTER                     :: kp

      NULLIFY (kp, occupation_numbers)
      CPASSERT(ASSOCIATED(kpoints))
      CPASSERT(ASSOCIATED(kpoints%kp_env))

      DO ikpoint_local = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(ikpoint_local)%kpoint_env
         CPASSERT(ASSOCIATED(kp))
         CPASSERT(ASSOCIATED(kp%mos))
         DO ispin = 1, SIZE(kp%mos, 2)
            DO icomponent = 1, SIZE(kp%mos, 1)
               CALL get_mo_set(kp%mos(icomponent, ispin), homo=nocc, maxocc=maxocc, &
                               nelectron=nelectron, nmo=nmo, occupation_numbers=occupation_numbers)
               CPASSERT(ASSOCIATED(occupation_numbers))
               IF (nocc < 1 .OR. nocc > nmo) THEN
                  CPABORT("Fixed-occupation k-point OT requires a nonempty occupied subspace.")
               END IF
               IF (ABS(maxocc*REAL(nocc, KIND=dp) - REAL(nelectron, KIND=dp)) > &
                   100.0_dp*EPSILON(1.0_dp)) THEN
                  CPABORT("Fixed-occupation k-point OT requires uniform integer occupations.")
               END IF
               occupation_numbers(:) = 0.0_dp
               occupation_numbers(1:nocc) = maxocc
               kp%mos(icomponent, ispin)%kTS = 0.0_dp
            END DO
         END DO
      END DO

   END SUBROUTINE qs_kpoint_set_fixed_occupations

END MODULE qs_kpoint_state
