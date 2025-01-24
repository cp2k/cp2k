!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for an external energy correction on top of a Kohn-Sham calculation
!> \par History
!>       10.2024 created
!> \author JGH
! **************************************************************************************************
MODULE ec_external
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cell_types,                      ONLY: cell_type
   USE core_ppl,                        ONLY: build_core_ppl
   USE core_ppnl,                       ONLY: build_core_ppnl
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_add,&
                                              dbcsr_copy,&
                                              dbcsr_create,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_type,&
                                              dbcsr_type_symmetric
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_plus_fm_fm_t,&
                                              cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_diag,&
                                              cp_fm_get_info,&
                                              cp_fm_maxabsval,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE ec_env_types,                    ONLY: energy_correction_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE mathlib,                         ONLY: det_3x3
   USE message_passing,                 ONLY: mp_para_env_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: pascal
   USE qs_core_energies,                ONLY: calculate_ptrace
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_kinetic,                      ONLY: build_kinetic_matrix
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_overlap,                      ONLY: build_overlap_matrix
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE virial_methods,                  ONLY: one_third_sum_diag
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ec_external'

   PUBLIC :: ec_ext_energy

CONTAINS

! **************************************************************************************************
!> \brief External energy method
!> \param qs_env ...
!> \param ec_env ...
!> \param calculate_forces ...
!> \par History
!>       10.2024 created
!> \author JGH
! **************************************************************************************************
   SUBROUTINE ec_ext_energy(qs_env, ec_env, calculate_forces)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(energy_correction_type), POINTER              :: ec_env
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(len=*), PARAMETER                        :: routineN = 'ec_ext_energy'

      INTEGER                                            :: handle, ispin, nocc, nspins, unit_nr
      REAL(KIND=dp)                                      :: focc
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: cpmos, mo_occ, mo_ref
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, dft_control=dft_control)
      nspins = dft_control%nspins

      logger => cp_get_default_logger()
      IF (logger%para_env%is_source()) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      END IF

      CALL get_qs_env(qs_env, mos=mos)
      ALLOCATE (cpmos(nspins), mo_ref(nspins), mo_occ(nspins))
      DO ispin = 1, nspins
         CALL get_mo_set(mo_set=mos(ispin), mo_coeff=mo_coeff, homo=nocc)
         NULLIFY (fm_struct)
         CALL cp_fm_struct_create(fm_struct, ncol_global=nocc, &
                                  template_fmstruct=mo_coeff%matrix_struct)
         CALL cp_fm_create(cpmos(ispin), fm_struct)
         CALL cp_fm_set_all(cpmos(ispin), 0.0_dp)
         CALL cp_fm_create(mo_ref(ispin), fm_struct)
         CALL cp_fm_set_all(mo_ref(ispin), 0.0_dp)
         CALL cp_fm_create(mo_occ(ispin), fm_struct)
         CALL cp_fm_to_fm(mo_coeff, mo_occ(ispin), nocc)
         CALL cp_fm_struct_release(fm_struct)
      END DO

      CALL cp_fm_release(ec_env%mo_occ)
      CALL cp_fm_release(ec_env%cpmos)
      IF (ec_env%debug_external) THEN
         !
         ec_env%mo_occ => mo_ref
         CALL ec_ext_debug(qs_env, ec_env, calculate_forces, unit_nr)
         !
         IF (calculate_forces) THEN
            focc = 2.0_dp
            IF (nspins == 1) focc = 4.0_dp
            DO ispin = 1, nspins
               CALL get_mo_set(mo_set=mos(ispin), mo_coeff=mo_coeff, homo=nocc)
               CALL cp_dbcsr_sm_fm_multiply(ec_env%matrix_h(1, 1)%matrix, ec_env%mo_occ(ispin), &
                                            cpmos(ispin), nocc, &
                                            alpha=focc, beta=0.0_dp)
            END DO
         END IF
         ec_env%cpmos => cpmos
      ELSE
         ! get external information
         CALL ec_ext_interface(qs_env, ec_env, mo_ref, cpmos, calculate_forces, unit_nr)
         ec_env%mo_occ => mo_ref
         ec_env%cpmos => cpmos
      END IF

      IF (calculate_forces) THEN
         ! check for orbital rotations
         CALL get_qs_env(qs_env, matrix_s=matrix_s)
         DO ispin = 1, nspins
            CALL align_vectors(ec_env%cpmos(ispin), ec_env%mo_occ(ispin), mo_occ(ispin), &
                               matrix_s(1)%matrix, unit_nr)
         END DO
         ! set up matrices for response
         CALL ec_ext_setup(qs_env, ec_env, .TRUE., unit_nr)
         ! orthogonality force
         CALL matrix_r_forces(qs_env, ec_env%cpmos, ec_env%mo_occ, &
                              ec_env%matrix_w(1, 1)%matrix, unit_nr)
      ELSE
         CALL ec_ext_setup(qs_env, ec_env, .FALSE., unit_nr)
      END IF

      CALL cp_fm_release(mo_occ)

      CALL timestop(handle)

   END SUBROUTINE ec_ext_energy

! **************************************************************************************************

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param ec_env ...
!> \param mo_ref ...
!> \param cpmos ...
!> \param calculate_forces ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE ec_ext_interface(qs_env, ec_env, mo_ref, cpmos, calculate_forces, unit_nr)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(energy_correction_type), POINTER              :: ec_env
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: mo_ref, cpmos
      LOGICAL, INTENT(IN)                                :: calculate_forces
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'ec_ext_interface'

      INTEGER                                            :: handle, ispin, nao, nocc(2), nspins
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      MARK_USED(qs_env)
      MARK_USED(ec_env)
      MARK_USED(mo_ref)
      MARK_USED(cpmos)
      MARK_USED(calculate_forces)
      MARK_USED(unit_nr)

      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, para_env=para_env)

      nspins = dft_control%nspins
      CALL get_qs_env(qs_env, mos=mos)
      nocc = 0
      DO ispin = 1, nspins
         CALL get_mo_set(mo_set=mos(ispin), nao=nao, homo=nocc(ispin))
      END DO

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A)') "Read EXTERNAL Response from file: "//TRIM(ec_env%exresp_fn)
      END IF

      CALL timestop(handle)

   END SUBROUTINE ec_ext_interface

! **************************************************************************************************

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param ec_env ...
!> \param calculate_forces ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE ec_ext_debug(qs_env, ec_env, calculate_forces, unit_nr)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(energy_correction_type), POINTER              :: ec_env
      LOGICAL, INTENT(IN)                                :: calculate_forces
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'ec_ext_debug'

      CHARACTER(LEN=default_string_length)               :: headline
      INTEGER                                            :: handle, ispin, nocc, nspins
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_h, matrix_p, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(qs_rho_type), POINTER                         :: rho

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      sab_orb=sab_orb, &
                      rho=rho, &
                      matrix_s_kp=matrix_s, &
                      matrix_h_kp=matrix_h)

      nspins = dft_control%nspins
      CALL get_qs_env(qs_env, mos=mos)
      DO ispin = 1, nspins
         CALL get_mo_set(mo_set=mos(ispin), mo_coeff=mo_coeff, homo=nocc)
         CALL cp_fm_to_fm(mo_coeff, ec_env%mo_occ(ispin), nocc)
      END DO

      ! Core Hamiltonian matrix
      IF (ASSOCIATED(ec_env%matrix_h)) CALL dbcsr_deallocate_matrix_set(ec_env%matrix_h)
      CALL dbcsr_allocate_matrix_set(ec_env%matrix_h, 1, 1)
      headline = "CORE HAMILTONIAN MATRIX"
      ALLOCATE (ec_env%matrix_h(1, 1)%matrix)
      CALL dbcsr_create(ec_env%matrix_h(1, 1)%matrix, name=TRIM(headline), &
                        template=matrix_h(1, 1)%matrix, matrix_type=dbcsr_type_symmetric)
      CALL cp_dbcsr_alloc_block_from_nbl(ec_env%matrix_h(1, 1)%matrix, sab_orb)
      CALL dbcsr_copy(ec_env%matrix_h(1, 1)%matrix, matrix_h(1, 1)%matrix)

      ! Get density matrix of reference calculation
      CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
      ! Use Core energy as model energy
      CALL calculate_ptrace(ec_env%matrix_h, matrix_p, ec_env%ex, nspins)

      IF (calculate_forces) THEN
         ! force of model energy
         CALL ec_debug_force(qs_env, matrix_p, unit_nr)
      END IF

      CALL timestop(handle)

   END SUBROUTINE ec_ext_debug

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param matrix_p ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE ec_debug_force(qs_env, matrix_p, unit_nr)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_p
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'ec_debug_force'

      INTEGER                                            :: handle, iounit, nder, nimages
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: calculate_forces, debug_forces, &
                                                            debug_stress, use_virial
      REAL(KIND=dp)                                      :: eps_ppnl, fconv
      REAL(KIND=dp), DIMENSION(3)                        :: fodeb
      REAL(KIND=dp), DIMENSION(3, 3)                     :: stdeb, sttot
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: scrm
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb, sac_ppl, sap_ppnl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      debug_forces = .TRUE.
      debug_stress = .TRUE.
      iounit = unit_nr

      calculate_forces = .TRUE.

      ! no k-points possible
      NULLIFY (cell, dft_control, force, ks_env, para_env, virial)
      CALL get_qs_env(qs_env=qs_env, &
                      cell=cell, &
                      dft_control=dft_control, &
                      force=force, &
                      ks_env=ks_env, &
                      para_env=para_env, &
                      virial=virial)
      nimages = dft_control%nimages
      IF (nimages /= 1) THEN
         CPABORT("K-points not implemented")
      END IF

      ! check for virial
      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)

      fconv = 1.0E-9_dp*pascal/cell%deth
      IF (debug_stress .AND. use_virial) THEN
         sttot = virial%pv_virial
      END IF

      ! check for GAPW/GAPW_XC
      IF (dft_control%qs_control%gapw .OR. dft_control%qs_control%gapw_xc) THEN
         CPABORT("GAPW not implemented")
      END IF

      ! get neighbor lists
      NULLIFY (sab_orb, sac_ppl, sap_ppnl)
      CALL get_qs_env(qs_env=qs_env, &
                      sab_orb=sab_orb, sac_ppl=sac_ppl, sap_ppnl=sap_ppnl)

      ! initialize src matrix
      NULLIFY (scrm)
      CALL dbcsr_allocate_matrix_set(scrm, 1, 1)
      ALLOCATE (scrm(1, 1)%matrix)
      CALL dbcsr_create(scrm(1, 1)%matrix, template=matrix_p(1, 1)%matrix)
      CALL cp_dbcsr_alloc_block_from_nbl(scrm(1, 1)%matrix, sab_orb)

      nder = 1
      IF (SIZE(matrix_p, 1) == 2) THEN
         CALL dbcsr_add(matrix_p(1, 1)%matrix, matrix_p(2, 1)%matrix, &
                        alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
      END IF

      ! kinetic energy
      IF (debug_forces) fodeb(1:3) = force(1)%kinetic(1:3, 1)
      IF (debug_stress .AND. use_virial) stdeb = virial%pv_ekinetic
      CALL build_kinetic_matrix(ks_env, matrixkp_t=scrm, &
                                matrix_name="KINETIC ENERGY MATRIX", &
                                basis_type="ORB", &
                                sab_nl=sab_orb, calculate_forces=.TRUE., &
                                matrixkp_p=matrix_p)
      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%kinetic(1:3, 1) - fodeb(1:3)
         CALL para_env%sum(fodeb)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pout*dT    ", fodeb
      END IF
      IF (debug_stress .AND. use_virial) THEN
         stdeb = fconv*(virial%pv_ekinetic - stdeb)
         CALL para_env%sum(stdeb)
         IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
            'STRESS| Pout*dT', one_third_sum_diag(stdeb), det_3x3(stdeb)
      END IF
      IF (SIZE(matrix_p, 1) == 2) THEN
         CALL dbcsr_add(matrix_p(1, 1)%matrix, matrix_p(2, 1)%matrix, &
                        alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
      END IF

      ! compute the ppl contribution to the core hamiltonian
      NULLIFY (atomic_kind_set, particle_set, qs_kind_set)
      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set, particle_set=particle_set, &
                      atomic_kind_set=atomic_kind_set)

      IF (ASSOCIATED(sac_ppl)) THEN
         IF (calculate_forces .AND. debug_forces) fodeb(1:3) = force(1)%gth_ppl(1:3, 1)
         IF (debug_stress .AND. use_virial) stdeb = virial%pv_ppl
         CALL build_core_ppl(scrm, matrix_p, force, &
                             virial, calculate_forces, use_virial, nder, &
                             qs_kind_set, atomic_kind_set, particle_set, sab_orb, sac_ppl, &
                             nimages, cell_to_index, "ORB")
         IF (calculate_forces .AND. debug_forces) THEN
            fodeb(1:3) = force(1)%gth_ppl(1:3, 1) - fodeb(1:3)
            CALL para_env%sum(fodeb)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pout*dH_PPL ", fodeb
         END IF
         IF (debug_stress .AND. use_virial) THEN
            stdeb = fconv*(virial%pv_ppl - stdeb)
            CALL para_env%sum(stdeb)
            IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
               'STRESS| Pout*dH_PPL', one_third_sum_diag(stdeb), det_3x3(stdeb)
         END IF
      END IF

      ! compute the ppnl contribution to the core hamiltonian ***
      eps_ppnl = dft_control%qs_control%eps_ppnl
      IF (ASSOCIATED(sap_ppnl)) THEN
         IF (calculate_forces .AND. debug_forces) fodeb(1:3) = force(1)%gth_ppnl(1:3, 1)
         IF (debug_stress .AND. use_virial) stdeb = virial%pv_ppnl
         CALL build_core_ppnl(scrm, matrix_p, force, &
                              virial, calculate_forces, use_virial, nder, &
                              qs_kind_set, atomic_kind_set, particle_set, &
                              sab_orb, sap_ppnl, eps_ppnl, &
                              nimages, cell_to_index, "ORB")
         IF (calculate_forces .AND. debug_forces) THEN
            fodeb(1:3) = force(1)%gth_ppnl(1:3, 1) - fodeb(1:3)
            CALL para_env%sum(fodeb)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pout*dH_PPNL", fodeb
         END IF
         IF (debug_stress .AND. use_virial) THEN
            stdeb = fconv*(virial%pv_ppnl - stdeb)
            CALL para_env%sum(stdeb)
            IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
               'STRESS| Pout*dH_PPNL', one_third_sum_diag(stdeb), det_3x3(stdeb)
         END IF
      END IF

      IF (debug_stress .AND. use_virial) THEN
         stdeb = fconv*(virial%pv_virial - sttot)
         CALL para_env%sum(stdeb)
         IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
            'STRESS| Stress Pout*dHcore   ', one_third_sum_diag(stdeb), det_3x3(stdeb)
         IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") ' '
      END IF

      ! delete scr matrix
      CALL dbcsr_deallocate_matrix_set(scrm)

      CALL timestop(handle)

   END SUBROUTINE ec_debug_force

! **************************************************************************************************

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param ec_env ...
!> \param calc_forces ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE ec_ext_setup(qs_env, ec_env, calc_forces, unit_nr)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(energy_correction_type), POINTER              :: ec_env
      LOGICAL, INTENT(IN)                                :: calc_forces
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'ec_ext_setup'

      CHARACTER(LEN=default_string_length)               :: headline
      INTEGER                                            :: handle, ispin, nao, nocc, nspins
      REAL(KIND=dp)                                      :: a_max, c_max
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct, mat_struct
      TYPE(cp_fm_type)                                   :: emat, ksmo, smo
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_h, matrix_ks, matrix_p, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(qs_rho_type), POINTER                         :: rho

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      sab_orb=sab_orb, &
                      rho=rho, &
                      matrix_s_kp=matrix_s, &
                      matrix_ks_kp=matrix_ks, &
                      matrix_h_kp=matrix_h)

      nspins = dft_control%nspins

      ! KS Hamiltonian matrix
      IF (ASSOCIATED(ec_env%matrix_ks)) CALL dbcsr_deallocate_matrix_set(ec_env%matrix_ks)
      CALL dbcsr_allocate_matrix_set(ec_env%matrix_ks, nspins, 1)
      headline = "HAMILTONIAN MATRIX"
      DO ispin = 1, nspins
         ALLOCATE (ec_env%matrix_ks(ispin, 1)%matrix)
         CALL dbcsr_create(ec_env%matrix_ks(ispin, 1)%matrix, name=TRIM(headline), &
                           template=matrix_ks(ispin, 1)%matrix, matrix_type=dbcsr_type_symmetric)
         CALL cp_dbcsr_alloc_block_from_nbl(ec_env%matrix_ks(ispin, 1)%matrix, sab_orb)
         CALL dbcsr_copy(ec_env%matrix_ks(ispin, 1)%matrix, matrix_ks(ispin, 1)%matrix)
      END DO

      ! Overlap matrix
      IF (ASSOCIATED(ec_env%matrix_s)) CALL dbcsr_deallocate_matrix_set(ec_env%matrix_s)
      CALL dbcsr_allocate_matrix_set(ec_env%matrix_s, 1, 1)
      headline = "OVERLAP MATRIX"
      ALLOCATE (ec_env%matrix_s(1, 1)%matrix)
      CALL dbcsr_create(ec_env%matrix_s(1, 1)%matrix, name=TRIM(headline), &
                        template=matrix_s(1, 1)%matrix, matrix_type=dbcsr_type_symmetric)
      CALL cp_dbcsr_alloc_block_from_nbl(ec_env%matrix_s(1, 1)%matrix, sab_orb)
      CALL dbcsr_copy(ec_env%matrix_s(1, 1)%matrix, matrix_s(1, 1)%matrix)

      ! density matrix
      ! Get density matrix of reference calculation
      CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
      IF (ASSOCIATED(ec_env%matrix_p)) CALL dbcsr_deallocate_matrix_set(ec_env%matrix_p)
      CALL dbcsr_allocate_matrix_set(ec_env%matrix_p, nspins, 1)
      headline = "DENSITY MATRIX"
      DO ispin = 1, nspins
         ALLOCATE (ec_env%matrix_p(ispin, 1)%matrix)
         CALL dbcsr_create(ec_env%matrix_p(ispin, 1)%matrix, name=TRIM(headline), &
                           template=matrix_p(ispin, 1)%matrix, matrix_type=dbcsr_type_symmetric)
         CALL cp_dbcsr_alloc_block_from_nbl(ec_env%matrix_p(ispin, 1)%matrix, sab_orb)
         CALL dbcsr_copy(ec_env%matrix_p(ispin, 1)%matrix, matrix_p(ispin, 1)%matrix)
      END DO

      IF (calc_forces) THEN
         ! energy weighted density matrix
         ! for security, we recalculate W here (stored in qs_env)
         IF (ASSOCIATED(ec_env%matrix_w)) CALL dbcsr_deallocate_matrix_set(ec_env%matrix_w)
         CALL dbcsr_allocate_matrix_set(ec_env%matrix_w, nspins, 1)
         headline = "ENERGY WEIGHTED DENSITY MATRIX"
         DO ispin = 1, nspins
            ALLOCATE (ec_env%matrix_w(ispin, 1)%matrix)
            CALL dbcsr_create(ec_env%matrix_w(ispin, 1)%matrix, name=TRIM(headline), &
                              template=matrix_p(ispin, 1)%matrix, matrix_type=dbcsr_type_symmetric)
            CALL cp_dbcsr_alloc_block_from_nbl(ec_env%matrix_w(ispin, 1)%matrix, sab_orb)
            CALL dbcsr_set(ec_env%matrix_w(ispin, 1)%matrix, 0.0_dp)
         END DO

         ! hz matrix
         IF (ASSOCIATED(ec_env%matrix_hz)) CALL dbcsr_deallocate_matrix_set(ec_env%matrix_hz)
         CALL dbcsr_allocate_matrix_set(ec_env%matrix_hz, nspins)
         headline = "Hz MATRIX"
         DO ispin = 1, nspins
            ALLOCATE (ec_env%matrix_hz(ispin)%matrix)
            CALL dbcsr_create(ec_env%matrix_hz(ispin)%matrix, name=TRIM(headline), &
                              template=matrix_s(1, 1)%matrix, matrix_type=dbcsr_type_symmetric)
            CALL cp_dbcsr_alloc_block_from_nbl(ec_env%matrix_hz(ispin)%matrix, sab_orb)
            CALL dbcsr_set(ec_env%matrix_hz(ispin)%matrix, 0.0_dp)
         END DO

         ! Test for consistency of orbitals and KS matrix
         DO ispin = 1, nspins
            mat_struct => ec_env%mo_occ(ispin)%matrix_struct
            CALL cp_fm_create(ksmo, mat_struct)
            CALL cp_fm_get_info(ksmo, nrow_global=nao, ncol_global=nocc)
            CALL cp_dbcsr_sm_fm_multiply(ec_env%matrix_ks(ispin, 1)%matrix, ec_env%mo_occ(ispin), &
                                         ksmo, nocc, alpha=1.0_dp, beta=0.0_dp)
            CALL cp_fm_create(smo, mat_struct)
            CALL cp_dbcsr_sm_fm_multiply(ec_env%matrix_s(1, 1)%matrix, ec_env%mo_occ(ispin), &
                                         smo, nocc, alpha=1.0_dp, beta=0.0_dp)
            CALL cp_fm_struct_create(fm_struct, ncol_global=nocc, nrow_global=nocc, &
                                     template_fmstruct=mat_struct)
            CALL cp_fm_create(emat, fm_struct)
            CALL parallel_gemm('T', 'N', nocc, nocc, nao, 1.0_dp, ec_env%mo_occ(ispin), ksmo, 0.0_dp, emat)
            CALL parallel_gemm('N', 'N', nao, nocc, nocc, -1.0_dp, smo, emat, 1.0_dp, ksmo)
            CALL cp_fm_maxabsval(ksmo, a_max)
            CALL cp_fm_struct_release(fm_struct)
            CALL cp_fm_release(smo)
            CALL cp_fm_release(ksmo)
            CALL cp_fm_release(emat)
            CALL cp_fm_maxabsval(ec_env%mo_occ(ispin), c_max)
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, "(T3,A,T50,I2,T61,F20.12)") "External:: Max value of MO coeficients", ispin, c_max
               WRITE (unit_nr, "(T3,A,T50,I2,T61,F20.12)") "External:: Max value of MO gradients", ispin, a_max
            END IF
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE ec_ext_setup

! **************************************************************************************************
!> \brief ...
!> \param cpmos ...
!> \param mo_ref ...
!> \param mo_occ ...
!> \param matrix_s ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE align_vectors(cpmos, mo_ref, mo_occ, matrix_s, unit_nr)
      TYPE(cp_fm_type), INTENT(IN)                       :: cpmos, mo_ref, mo_occ
      TYPE(dbcsr_type), POINTER                          :: matrix_s
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'align_vectors'

      INTEGER                                            :: handle, i, nao, nocc, scg
      REAL(KIND=dp)                                      :: a_max
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: diag
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct, mat_struct
      TYPE(cp_fm_type)                                   :: emat, smo

      CALL timeset(routineN, handle)

      mat_struct => cpmos%matrix_struct
      CALL cp_fm_create(smo, mat_struct)
      CALL cp_fm_get_info(smo, nrow_global=nao, ncol_global=nocc)
      CALL cp_dbcsr_sm_fm_multiply(matrix_s, mo_occ, smo, nocc, alpha=1.0_dp, beta=0.0_dp)
      CALL cp_fm_struct_create(fm_struct, ncol_global=nocc, nrow_global=nocc, &
                               template_fmstruct=mat_struct)
      CALL cp_fm_create(emat, fm_struct)
      CALL parallel_gemm('T', 'N', nocc, nocc, nao, 1.0_dp, mo_ref, smo, 0.0_dp, emat)
      CALL parallel_gemm('N', 'N', nao, nocc, nocc, 1.0_dp, cpmos, emat, 0.0_dp, smo)
      CALL cp_fm_to_fm(smo, cpmos)
      CALL cp_fm_to_fm(mo_occ, mo_ref)
      !
      ALLOCATE (diag(nocc))
      CALL cp_fm_get_diag(emat, diag)
      a_max = nocc - SUM(diag)
      scg = 0
      DO i = 1, nocc
         IF (ABS(diag(i) + 1.0_dp) < 0.001) scg = scg + 1
      END DO
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, "(T3,A,T61,F20.8)") "External:: Orbital rotation index", a_max
         WRITE (unit_nr, "(T3,A,T71,I10)") "External:: Number of orbital phase changes", scg
      END IF

      DEALLOCATE (diag)
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_release(smo)
      CALL cp_fm_release(emat)

      CALL timestop(handle)

   END SUBROUTINE align_vectors

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param matrix_w ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE matrix_w_forces(qs_env, matrix_w, unit_nr)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_w
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'matrix_w_forces'

      INTEGER                                            :: handle, iounit, nder, nimages
      LOGICAL                                            :: debug_forces, debug_stress, use_virial
      REAL(KIND=dp)                                      :: fconv
      REAL(KIND=dp), DIMENSION(3)                        :: fodeb
      REAL(KIND=dp), DIMENSION(3, 3)                     :: stdeb, sttot
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: scrm
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      debug_forces = .TRUE.
      debug_stress = .TRUE.

      iounit = unit_nr

      ! no k-points possible
      CALL get_qs_env(qs_env=qs_env, &
                      cell=cell, &
                      dft_control=dft_control, &
                      force=force, &
                      ks_env=ks_env, &
                      sab_orb=sab_orb, &
                      para_env=para_env, &
                      virial=virial)
      nimages = dft_control%nimages
      IF (nimages /= 1) THEN
         CPABORT("K-points not implemented")
      END IF

      ! check for virial
      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)

      fconv = 1.0E-9_dp*pascal/cell%deth
      IF (debug_stress .AND. use_virial) THEN
         sttot = virial%pv_virial
      END IF

      ! initialize src matrix
      NULLIFY (scrm)
      CALL dbcsr_allocate_matrix_set(scrm, 1, 1)
      ALLOCATE (scrm(1, 1)%matrix)
      CALL dbcsr_create(scrm(1, 1)%matrix, template=matrix_w(1, 1)%matrix)
      CALL cp_dbcsr_alloc_block_from_nbl(scrm(1, 1)%matrix, sab_orb)

      nder = 1
      IF (SIZE(matrix_w, 1) == 2) THEN
         CALL dbcsr_add(matrix_w(1, 1)%matrix, matrix_w(2, 1)%matrix, &
                        alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
      END IF

      ! Overlap and kinetic energy matrices
      IF (debug_forces) fodeb(1:3) = force(1)%overlap(1:3, 1)
      IF (debug_stress .AND. use_virial) stdeb = virial%pv_overlap
      CALL build_overlap_matrix(ks_env, matrixkp_s=scrm, &
                                matrix_name="OVERLAP MATRIX", &
                                basis_type_a="ORB", &
                                basis_type_b="ORB", &
                                sab_nl=sab_orb, calculate_forces=.TRUE., &
                                matrixkp_p=matrix_w)

      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%overlap(1:3, 1) - fodeb(1:3)
         CALL para_env%sum(fodeb)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Wout*dS    ", fodeb
      END IF
      IF (debug_stress .AND. use_virial) THEN
         stdeb = fconv*(virial%pv_overlap - stdeb)
         CALL para_env%sum(stdeb)
         IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
            'STRESS| Wout*dS', one_third_sum_diag(stdeb), det_3x3(stdeb)
      END IF
      IF (SIZE(matrix_w, 1) == 2) THEN
         CALL dbcsr_add(matrix_w(1, 1)%matrix, matrix_w(2, 1)%matrix, &
                        alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
      END IF

      ! delete scrm matrix
      CALL dbcsr_deallocate_matrix_set(scrm)

      CALL timestop(handle)

   END SUBROUTINE matrix_w_forces

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param cpmos ...
!> \param mo_occ ...
!> \param matrix_r ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE matrix_r_forces(qs_env, cpmos, mo_occ, matrix_r, unit_nr)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: cpmos, mo_occ
      TYPE(dbcsr_type), POINTER                          :: matrix_r
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'matrix_r_forces'

      INTEGER                                            :: handle, iounit, ispin, nao, nocc, nspins
      LOGICAL                                            :: debug_forces, debug_stress, use_virial
      REAL(KIND=dp)                                      :: fconv, focc
      REAL(KIND=dp), DIMENSION(3)                        :: fodeb
      REAL(KIND=dp), DIMENSION(3, 3)                     :: stdeb
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct, mat_struct
      TYPE(cp_fm_type)                                   :: chcmat, rcvec
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: scrm
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      debug_forces = .TRUE.
      debug_stress = .TRUE.
      iounit = unit_nr

      nspins = SIZE(mo_occ)
      focc = 1.0_dp
      IF (nspins == 1) focc = 2.0_dp
      focc = 0.25_dp*focc

      CALL dbcsr_set(matrix_r, 0.0_dp)
      DO ispin = 1, nspins
         CALL cp_fm_get_info(cpmos(ispin), matrix_struct=fm_struct, nrow_global=nao, ncol_global=nocc)
         CALL cp_fm_create(rcvec, fm_struct)
         CALL cp_fm_struct_create(mat_struct, nrow_global=nocc, ncol_global=nocc, template_fmstruct=fm_struct)
         CALL cp_fm_create(chcmat, mat_struct)
         CALL parallel_gemm("T", "N", nocc, nocc, nao, 1.0_dp, mo_occ(ispin), cpmos(ispin), 0.0_dp, chcmat)
         CALL parallel_gemm("N", "N", nao, nocc, nocc, 1.0_dp, mo_occ(ispin), chcmat, 0.0_dp, rcvec)
         CALL cp_dbcsr_plus_fm_fm_t(matrix_r, matrix_v=rcvec, matrix_g=mo_occ(ispin), ncol=nocc, alpha=focc)
         CALL cp_fm_struct_release(mat_struct)
         CALL cp_fm_release(rcvec)
         CALL cp_fm_release(chcmat)
      END DO

      CALL get_qs_env(qs_env=qs_env, &
                      cell=cell, &
                      dft_control=dft_control, &
                      force=force, &
                      ks_env=ks_env, &
                      sab_orb=sab_orb, &
                      para_env=para_env, &
                      virial=virial)
      ! check for virial
      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
      fconv = 1.0E-9_dp*pascal/cell%deth

      IF (debug_forces) fodeb(1:3) = force(1)%overlap(1:3, 1)
      IF (debug_stress .AND. use_virial) stdeb = virial%pv_overlap
      NULLIFY (scrm)
      CALL build_overlap_matrix(ks_env, matrix_s=scrm, &
                                matrix_name="OVERLAP MATRIX", &
                                basis_type_a="ORB", basis_type_b="ORB", &
                                sab_nl=sab_orb, calculate_forces=.TRUE., &
                                matrix_p=matrix_r)
      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%overlap(1:3, 1) - fodeb(1:3)
         CALL para_env%sum(fodeb)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Wz*dS ", fodeb
      END IF
      IF (debug_stress .AND. use_virial) THEN
         stdeb = fconv*(virial%pv_overlap - stdeb)
         CALL para_env%sum(stdeb)
         IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
            'STRESS| Wz   ', one_third_sum_diag(stdeb), det_3x3(stdeb)
      END IF
      CALL dbcsr_deallocate_matrix_set(scrm)

      CALL timestop(handle)

   END SUBROUTINE matrix_r_forces

END MODULE ec_external
