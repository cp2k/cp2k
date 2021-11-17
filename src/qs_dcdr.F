!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculate the derivatives of the MO coefficients wrt nuclear coordinates
!> \author Sandra Luber, Edward Ditler
! **************************************************************************************************

MODULE qs_dcdr

!#include "./common/cp_common_uses.f90"
   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_array_utils,                  ONLY: cp_2d_r_p_type
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm,&
                                              cp_fm_scale,&
                                              cp_fm_scale_and_add,&
                                              cp_fm_trace
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_diag,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE dbcsr_api,                       ONLY: dbcsr_add,&
                                              dbcsr_copy,&
                                              dbcsr_desymmetrize,&
                                              dbcsr_p_type,&
                                              dbcsr_set
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE molecule_types,                  ONLY: molecule_of_atom,&
                                              molecule_type
   USE particle_types,                  ONLY: particle_type
   USE qs_dcdr_ao,                      ONLY: apply_op_constant_term,&
                                              core_dR,&
                                              d_core_charge_density_dR,&
                                              d_vhxc_dR,&
                                              hr_mult_by_delta_1d,&
                                              vhxc_R_perturbed_basis_functions
   USE qs_dcdr_utils,                   ONLY: dcdr_read_restart,&
                                              dcdr_write_restart,&
                                              multiply_localization,&
                                              shift_wannier_into_cell
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_linres_methods,               ONLY: linres_solver
   USE qs_linres_types,                 ONLY: dcdr_env_type,&
                                              linres_control_type
   USE qs_moments,                      ONLY: build_local_moment_matrix,&
                                              dipole_deriv_ao
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_p_env_types,                  ONLY: qs_p_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: prepare_per_atom, dcdr_response_dR, dcdr_build_op_dR, apt_dR, apt_dR_localization

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_dcdr'

CONTAINS

! **************************************************************************************************
!> \brief Prepare the environment for a choice of lambda
!> \param dcdr_env ...
!> \param qs_env ...
!> \author Edward Ditler
! **************************************************************************************************
   SUBROUTINE prepare_per_atom(dcdr_env, qs_env)
      TYPE(dcdr_env_type)                                :: dcdr_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'prepare_per_atom', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, ispin, j, natom
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_all
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      NULLIFY (sab_all, qs_kind_set, particle_set)
      CALL get_qs_env(qs_env=qs_env, &
                      sab_all=sab_all, &
                      qs_kind_set=qs_kind_set, &
                      particle_set=particle_set)

      natom = SIZE(particle_set)
      IF (dcdr_env%distributed_origin) dcdr_env%ref_point(:) = particle_set(dcdr_env%lambda)%r(:)

      dcdr_env%delta_basis_function = 0._dp
      dcdr_env%delta_basis_function(:, dcdr_env%lambda) = 1._dp

      ! S matrix
      ! S1 = - < da/dr | b > * delta_a - < a | db/dr > * delta_b

      ! matrix_s(2:4) are anti-symmetric matrices and contain derivatives wrt. to < a |
      !               = < da/dR | b > = - < da/dr | b > = < a | db/dr >
      ! matrix_s1(2:4) = d/dR < a | b >
      !                and it's built as
      !                = - matrix_s      * delta_b  +  matrix_s      * delta_a
      !                = - < da/dR | b > * delta_b  +  < da/dR | b > * delta_a
      !                = + < da/dr | b > * delta_b  -  < da/dr | b > * delta_a
      !                = - < a | db/dr > * delta_b  -  < da/dr | b > * delta_a

      DO i = 1, 3
         ! S matrix
         CALL dbcsr_set(dcdr_env%matrix_nosym_temp(i)%matrix, 0._dp)
         CALL dbcsr_desymmetrize(dcdr_env%matrix_s(1 + i)%matrix, dcdr_env%matrix_s1(1 + i)%matrix)
         CALL dbcsr_desymmetrize(dcdr_env%matrix_s(1 + i)%matrix, dcdr_env%matrix_nosym_temp(i)%matrix)

         CALL hr_mult_by_delta_1d(dcdr_env%matrix_s1(1 + i)%matrix, qs_kind_set, "ORB", &
                                  sab_all, dcdr_env%lambda, direction_Or=.TRUE.)
         CALL hr_mult_by_delta_1d(dcdr_env%matrix_nosym_temp(i)%matrix, qs_kind_set, "ORB", &
                                  sab_all, dcdr_env%lambda, direction_Or=.FALSE.)

         CALL dbcsr_add(dcdr_env%matrix_s1(1 + i)%matrix, dcdr_env%matrix_nosym_temp(i)%matrix, -1._dp, +1._dp)
         CALL dbcsr_set(dcdr_env%matrix_nosym_temp(i)%matrix, 0._dp)

         ! T matrix
         CALL dbcsr_set(dcdr_env%matrix_nosym_temp(i)%matrix, 0._dp)
         CALL dbcsr_desymmetrize(dcdr_env%matrix_t(1 + i)%matrix, dcdr_env%matrix_t1(1 + i)%matrix)
         CALL dbcsr_desymmetrize(dcdr_env%matrix_t(1 + i)%matrix, dcdr_env%matrix_nosym_temp(i)%matrix)

         CALL hr_mult_by_delta_1d(dcdr_env%matrix_t1(1 + i)%matrix, qs_kind_set, "ORB", &
                                  sab_all, dcdr_env%lambda, direction_Or=.TRUE.)
         CALL hr_mult_by_delta_1d(dcdr_env%matrix_nosym_temp(i)%matrix, qs_kind_set, "ORB", &
                                  sab_all, dcdr_env%lambda, direction_Or=.FALSE.)

         CALL dbcsr_add(dcdr_env%matrix_t1(1 + i)%matrix, dcdr_env%matrix_nosym_temp(i)%matrix, -1._dp, +1._dp)
         CALL dbcsr_set(dcdr_env%matrix_nosym_temp(i)%matrix, 0._dp)
      END DO

      ! Operator:
      ispin = 1
      DO i = 1, 3
         CALL dbcsr_set(dcdr_env%matrix_ppnl_1(i)%matrix, 0.0_dp)
         CALL dbcsr_set(dcdr_env%matrix_hc(i)%matrix, 0.0_dp)
         CALL dbcsr_set(dcdr_env%matrix_vhxc_perturbed_basis(ispin, i)%matrix, 0.0_dp)
         CALL dbcsr_set(dcdr_env%matrix_vhxc_perturbed_basis(ispin, i + 3)%matrix, 0.0_dp)
         CALL dbcsr_set(dcdr_env%matrix_d_vhxc_dR(i)%matrix, 0.0_dp)
         CALL dbcsr_set(dcdr_env%matrix_core_charge_1(i)%matrix, 0.0_dp)
      END DO

      CALL core_dR(qs_env, dcdr_env)  ! dcdr_env%matrix_ppnl_1, hc
      CALL d_vhxc_dR(qs_env, dcdr_env)  ! dcdr_env%matrix_d_vhxc_dR
      CALL d_core_charge_density_dR(qs_env, dcdr_env)  ! dcdr_env%matrix_core_charge_1
      CALL vhxc_R_perturbed_basis_functions(qs_env, dcdr_env)  ! dcdr_env%matrix_vhxc_perturbed_basis

      ! APT:
      DO i = 1, 3
         DO j = 1, 3
            CALL dbcsr_set(dcdr_env%matrix_difdip(i, j)%matrix, 0._dp)
         END DO
      END DO

      CALL dipole_deriv_ao(qs_env, dcdr_env%matrix_difdip, dcdr_env%delta_basis_function, 1, dcdr_env%ref_point)

      CALL timestop(handle)
   END SUBROUTINE prepare_per_atom

! **************************************************************************************************
!> \brief Build the operator for the position perturbation
!> \param dcdr_env ...
!> \param qs_env ...
!> \authors SL, ED
! **************************************************************************************************
   SUBROUTINE dcdr_build_op_dR(dcdr_env, qs_env)

      TYPE(dcdr_env_type)                                :: dcdr_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dcdr_build_op_dR', &
         routineP = moduleN//':'//routineN
      REAL(KIND=dp), PARAMETER                           :: one = 1.0_dp, zero = 0.0_dp

      INTEGER                                            :: handle, nao, nmo
      TYPE(cp_fm_type), POINTER                          :: buf
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: opdr_sym

      CALL timeset(routineN, handle)

      nao = dcdr_env%nao
      nmo = dcdr_env%nmo

      ! allocate matrix for the sum of the perturbation terms of the operator (dbcsr matrix)
      NULLIFY (opdr_sym)
      CALL dbcsr_allocate_matrix_set(opdr_sym, 1)
      ALLOCATE (opdr_sym(1)%matrix)
      CALL dbcsr_copy(opdr_sym(1)%matrix, dcdr_env%matrix_s1(1)%matrix)  ! symmetric
      CALL dbcsr_set(opdr_sym(1)%matrix, 0.0_dp)

      CALL apply_op_constant_term(qs_env, dcdr_env)  ! dcdr_env%matrix_apply_op_constant
      ! Hartree and Exchange-Correlation contributions
      CALL dbcsr_add(opdr_sym(1)%matrix, dcdr_env%matrix_core_charge_1(dcdr_env%beta)%matrix, zero, one)
      CALL dbcsr_add(opdr_sym(1)%matrix, dcdr_env%matrix_d_vhxc_dR(dcdr_env%beta)%matrix, one, one)
      CALL dbcsr_add(opdr_sym(1)%matrix, dcdr_env%matrix_vhxc_perturbed_basis(1, dcdr_env%beta)%matrix, one, one)

      ! Core Hamiltonian contributions
      CALL dbcsr_add(opdr_sym(1)%matrix, dcdr_env%matrix_hc(dcdr_env%beta)%matrix, one, one)
      CALL dbcsr_add(opdr_sym(1)%matrix, dcdr_env%matrix_ppnl_1(dcdr_env%beta)%matrix, one, one)
      CALL dbcsr_add(opdr_sym(1)%matrix, dcdr_env%matrix_apply_op_constant(1)%matrix, one, one)

      CALL dbcsr_desymmetrize(opdr_sym(1)%matrix, dcdr_env%hamiltonian1(1)%matrix)
      CALL dbcsr_add(dcdr_env%hamiltonian1(1)%matrix, dcdr_env%matrix_t1(dcdr_env%beta + 1)%matrix, one, one)

      CALL cp_dbcsr_sm_fm_multiply(dcdr_env%hamiltonian1(1)%matrix, dcdr_env%mo_coeff(1)%matrix, &
                                   dcdr_env%op_dR(1)%matrix, ncol=nmo)

      ! The overlap derivative terms for the Sternheimer equation
      ! buf = mo * (-mo * matrix_ks * mo)
      CALL cp_fm_create(buf, dcdr_env%likemos_fm_struct)
      CALL cp_gemm('N', 'N', nao, nmo, nmo, &
                   -1.0_dp, dcdr_env%mo_coeff(1)%matrix, dcdr_env%chc(1)%matrix, &
                   0.0_dp, buf)

      CALL cp_dbcsr_sm_fm_multiply(dcdr_env%matrix_s1(dcdr_env%beta + 1)%matrix, buf, dcdr_env%op_dR(1)%matrix, &
                                   nmo, alpha=1.0_dp, beta=1.0_dp)
      CALL cp_fm_release(buf)

      ! SL multiply by -1 for response solver (H-S<H> C + dR_coupled= - (op_dR)
      CALL cp_fm_scale(-1.0_dp, dcdr_env%op_dR(1)%matrix)

      CALL dbcsr_deallocate_matrix_set(opdr_sym)

      CALL timestop(handle)
   END SUBROUTINE dcdr_build_op_dR

! **************************************************************************************************
!> \brief Get the dC/dR by solving the Sternheimer equation, using the op_dR matrix
!> \param dcdr_env ...
!> \param p_env ...
!> \param qs_env ...
!> \authors SL, ED
! **************************************************************************************************
   SUBROUTINE dcdr_response_dR(dcdr_env, p_env, qs_env)

      TYPE(dcdr_env_type)                                :: dcdr_env
      TYPE(qs_p_env_type), POINTER                       :: p_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dcdr_response_dR', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, output_unit
      LOGICAL                                            :: should_stop
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: h1_psi0, psi0_order, psi1
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(linres_control_type), POINTER                 :: linres_control
      TYPE(section_vals_type), POINTER                   :: lr_section

      CALL timeset(routineN, handle)
      NULLIFY (linres_control, lr_section, logger, psi1, h1_psi0)

      CALL get_qs_env(qs_env=qs_env, &
                      linres_control=linres_control)

      logger => cp_get_default_logger()
      lr_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%LINRES")

      output_unit = cp_print_key_unit_nr(logger, lr_section, "PRINT%PROGRAM_RUN_INFO", &
                                         extension=".linresLog")
      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(T10,A,/)") &
            "*** Self consistent optimization of the response wavefunction ***"
      END IF

      ! allocate the vectors
      ALLOCATE (psi0_order(1), psi1(1), h1_psi0(1))
      NULLIFY (psi1(1)%matrix, h1_psi0(1)%matrix)

      psi0_order(1)%matrix => dcdr_env%mo_coeff(1)%matrix
      CALL cp_fm_create(psi1(1)%matrix, dcdr_env%likemos_fm_struct)
      CALL cp_fm_create(h1_psi0(1)%matrix, dcdr_env%likemos_fm_struct)

      ! Restart
      IF (linres_control%linres_restart) THEN
         CALL dcdr_read_restart(qs_env, lr_section, psi1, dcdr_env%lambda, dcdr_env%beta, "dCdR")
      ELSE
         CALL cp_fm_set_all(psi1(1)%matrix, 0.0_dp)
      END IF

      IF (output_unit > 0) THEN
         WRITE (output_unit, "(T10,A,I4,A)") &
            "Response to the perturbation operator referring to atom ", dcdr_env%lambda, &
            " displaced in "//ACHAR(dcdr_env%beta + 119)
      END IF
      CALL cp_fm_set_all(dcdr_env%dCR(1)%matrix, 0.0_dp)
      CALL cp_fm_to_fm(dcdr_env%op_dR(1)%matrix, h1_psi0(1)%matrix)

      linres_control%lr_triplet = .FALSE. ! we do singlet response
      linres_control%do_kernel = .TRUE.
      linres_control%converged = .FALSE.

      ! Position perturbation to get dCR
      ! (H0-E0) psi1 = (H1-E1) psi0
      ! psi1 = the perturbed wavefunction
      ! h1_psi0 = (H1-E1-S1*\varepsilon)
      ! psi0_order = the unperturbed wavefunction
      CALL linres_solver(p_env, qs_env, psi1, h1_psi0, psi0_order, &
                         output_unit, should_stop)
      CALL cp_fm_to_fm(psi1(1)%matrix, dcdr_env%dCR(1)%matrix)

      ! Write the new result to the restart file
      IF (linres_control%linres_restart) THEN
         CALL dcdr_write_restart(qs_env, lr_section, psi1, dcdr_env%lambda, dcdr_env%beta, "dCdR")
      END IF

      ! clean up
      CALL cp_fm_release(psi1(1)%matrix)
      CALL cp_fm_release(h1_psi0(1)%matrix)
      DEALLOCATE (psi1, h1_psi0, psi0_order)
      CALL cp_print_key_finished_output(output_unit, logger, lr_section, &
                                        "PRINT%PROGRAM_RUN_INFO")

      CALL timestop(handle)

   END SUBROUTINE dcdr_response_dR

! **************************************************************************************************
!> \brief Calculate atomic polar tensor
!> \param qs_env ...
!> \param dcdr_env ...
!> \author Edward Ditler
! **************************************************************************************************
   SUBROUTINE apt_dR(qs_env, dcdr_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dcdr_env_type)                                :: dcdr_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'apt_dR', routineP = moduleN//':'//routineN

      INTEGER                                            :: alpha, handle, ikind, nao, nmo
      LOGICAL                                            :: ghost
      REAL(dp)                                           :: apt_basis_derivative, &
                                                            apt_coeff_derivative, charge, f_spin
      REAL(dp), DIMENSION(:, :, :), POINTER              :: apt_el, apt_nuc
      TYPE(cp_fm_p_type), POINTER                        :: overlap1_MO, tmp_fm_like_mos
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      apt_basis_derivative = 0._dp
      apt_coeff_derivative = 0._dp

      CALL timeset(routineN, handle)

      NULLIFY (qs_kind_set, particle_set)
      CALL get_qs_env(qs_env=qs_env, &
                      qs_kind_set=qs_kind_set, &
                      particle_set=particle_set)

      nao = dcdr_env%nao
      nmo = dcdr_env%nmo
      mo_coeff => dcdr_env%mo_coeff(1)%matrix
      apt_el => dcdr_env%apt_el_dcdr
      apt_nuc => dcdr_env%apt_nuc_dcdr

      f_spin = 2._dp

      ALLOCATE (tmp_fm_like_mos)
      ALLOCATE (overlap1_MO)
      CALL cp_fm_create(tmp_fm_like_mos%matrix, dcdr_env%likemos_fm_struct)
      CALL cp_fm_create(overlap1_MO%matrix, dcdr_env%momo_fm_struct)

      ! Compute S^(1,R)_(ij)
      CALL cp_fm_set_all(tmp_fm_like_mos%matrix, 0.0_dp)
      CALL cp_fm_scale_and_add(0._dp, dcdr_env%dCR_prime(1)%matrix, 1._dp, dcdr_env%dCR(1)%matrix)
      CALL cp_dbcsr_sm_fm_multiply(dcdr_env%matrix_s1(dcdr_env%beta + 1)%matrix, mo_coeff, &
                                   tmp_fm_like_mos%matrix, ncol=nmo)
      CALL cp_fm_gemm("T", "N", nmo, nmo, nao, &
                      1.0_dp, mo_coeff, tmp_fm_like_mos%matrix, &
                      0.0_dp, overlap1_MO%matrix)

      !   C^1 <- -dCR - 0.5 * mo_coeff @ S1_ij
      !    We get the negative of the coefficients out of the linres solver
      !    And apply the constant correction due to the overlap derivative.
      CALL cp_fm_gemm("N", "N", nao, nmo, nmo, &
                      -0.5_dp, mo_coeff, overlap1_MO%matrix, &
                      -1.0_dp, dcdr_env%dCR_prime(1)%matrix)

      DO alpha = 1, 3
         ! FIRST CONTRIBUTION: dCR * moments * mo
         CALL cp_fm_set_all(tmp_fm_like_mos%matrix, 0._dp)
         CALL dbcsr_desymmetrize(dcdr_env%matrix_s1(1)%matrix, dcdr_env%matrix_nosym_temp(1)%matrix)
         CALL dbcsr_desymmetrize(dcdr_env%moments(alpha)%matrix, dcdr_env%matrix_nosym_temp(2)%matrix)
         CALL dbcsr_add(dcdr_env%matrix_nosym_temp(1)%matrix, dcdr_env%matrix_nosym_temp(2)%matrix, &
                        -dcdr_env%ref_point(alpha), 1._dp)
         CALL cp_dbcsr_sm_fm_multiply(dcdr_env%matrix_nosym_temp(1)%matrix, dcdr_env%dCR_prime(1)%matrix, &
                                      tmp_fm_like_mos%matrix, ncol=nmo)
         CALL cp_fm_trace(mo_coeff, tmp_fm_like_mos%matrix, apt_coeff_derivative)

         apt_coeff_derivative = (-2._dp)*f_spin*apt_coeff_derivative
         apt_el(dcdr_env%beta, alpha, dcdr_env%lambda) &
            = apt_el(dcdr_env%beta, alpha, dcdr_env%lambda) + apt_coeff_derivative
      END DO ! alpha

      DO alpha = 1, 3
         ! SECOND CONTRIBUTION: We assemble all combinations of r_i, d(chi)/d(idir)
         ! difdip contains derivatives with respect to atom dcdr_env%lambda
         ! difdip(alpha, beta): < a | r_alpha | db/dR_beta >
         ! Multiply by the MO coefficients
         CALL cp_fm_set_all(tmp_fm_like_mos%matrix, 0.0_dp)
         CALL cp_dbcsr_sm_fm_multiply(dcdr_env%matrix_difdip(alpha, dcdr_env%beta)%matrix, mo_coeff, &
                                      tmp_fm_like_mos%matrix, ncol=nmo)
         CALL cp_fm_trace(mo_coeff, tmp_fm_like_mos%matrix, apt_basis_derivative)

         ! The negative sign is due to dipole_deriv_ao computing the derivatives with respect to nuclear coordinates.
         apt_basis_derivative = -f_spin*apt_basis_derivative
         apt_el(dcdr_env%beta, alpha, dcdr_env%lambda) = &
            apt_el(dcdr_env%beta, alpha, dcdr_env%lambda) + apt_basis_derivative

      END DO ! alpha

      ! Finally the nuclear contribution: nuclear charge * Kronecker_delta_{dcdr_env%beta,i}
      CALL get_atomic_kind(particle_set(dcdr_env%lambda)%atomic_kind, kind_number=ikind)
      CALL get_qs_kind(qs_kind_set(ikind), core_charge=charge, ghost=ghost)
      IF (.NOT. ghost) THEN
         apt_nuc(dcdr_env%beta, dcdr_env%beta, dcdr_env%lambda) = &
            apt_nuc(dcdr_env%beta, dcdr_env%beta, dcdr_env%lambda) + charge
      END IF

      ! And deallocate all the things!
      CALL cp_fm_release(tmp_fm_like_mos%matrix)
      CALL cp_fm_release(overlap1_MO%matrix)
      DEALLOCATE (tmp_fm_like_mos)
      DEALLOCATE (overlap1_MO)

      CALL timestop(handle)
   END SUBROUTINE apt_dR

! **************************************************************************************************
!> \brief Calculate atomic polar tensor using the localized dipole operator
!> \param qs_env ...
!> \param dcdr_env ...
!> \author Edward Ditler
! **************************************************************************************************
   SUBROUTINE apt_dR_localization(qs_env, dcdr_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dcdr_env_type)                                :: dcdr_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'apt_dR_localization', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: alpha, handle, i, icenter, ikind, &
                                                            map_atom, map_molecule, nao, natom, &
                                                            nmo, nsubset
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: mapping_atom_molecule, &
                                                            mapping_wannier_atom
      LOGICAL                                            :: ghost
      REAL(dp)                                           :: apt_basis_derivative, &
                                                            apt_coeff_derivative, charge, f_spin, &
                                                            smallest_r, this_factor, tmp_aptcontr, &
                                                            tmp_r
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: diagonal_elements
      REAL(dp), DIMENSION(3)                             :: distance, r_shifted
      REAL(dp), DIMENSION(:, :, :), POINTER              :: apt_el, apt_nuc
      REAL(dp), DIMENSION(:, :, :, :), POINTER           :: apt_center, apt_subset
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_2d_r_p_type), DIMENSION(:), POINTER        :: centers_set
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: tmp_fm_like_mos, tmp_fm_momo
      TYPE(cp_fm_p_type), POINTER                        :: tmp_fm
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      NULLIFY (qs_kind_set, particle_set, molecule_set, cell)

      CALL get_qs_env(qs_env=qs_env, &
                      qs_kind_set=qs_kind_set, &
                      particle_set=particle_set, &
                      molecule_set=molecule_set, &
                      cell=cell)

      nsubset = SIZE(molecule_set)
      natom = SIZE(particle_set)
      apt_el => dcdr_env%apt_el_dcdr
      apt_nuc => dcdr_env%apt_nuc_dcdr
      apt_subset => dcdr_env%apt_el_dcdr_per_subset
      apt_center => dcdr_env%apt_el_dcdr_per_center

      ! Map wannier functions to atoms
      centers_set => dcdr_env%centers_set
      ALLOCATE (mapping_wannier_atom(dcdr_env%nbr_center(1)))
      DO icenter = 1, dcdr_env%nbr_center(1)
         ! For every center we check which atom is closest
         CALL shift_wannier_into_cell(r=centers_set(1)%array(1:3, icenter), &
                                      cell=cell, &
                                      r_shifted=r_shifted)

         smallest_r = HUGE(0._dp)
         DO i = 1, natom
            distance = pbc(r_shifted, particle_set(i)%r(1:3), cell)
            tmp_r = SUM(distance**2)
            IF (tmp_r < smallest_r) THEN
               mapping_wannier_atom(icenter) = i
               smallest_r = tmp_r
            END IF
         END DO
      END DO

      ! Map atoms to molecules
      ALLOCATE (mapping_atom_molecule(natom))
      CALL molecule_of_atom(molecule_set, atom_to_mol=mapping_atom_molecule)
      IF (dcdr_env%lambda == 1 .AND. dcdr_env%beta == 1) THEN
         ! PRINT *, 'Mapping of atoms to subsets'
         DO icenter = 1, dcdr_env%nbr_center(1)
            map_atom = mapping_wannier_atom(icenter)
            map_molecule = mapping_atom_molecule(map_atom)
            ! PRINT *, 'Center', icenter, 'map_atom', map_atom
         END DO
      END IF

      nao = dcdr_env%nao
      nmo = dcdr_env%nmo
      mo_coeff => dcdr_env%mo_coeff(1)%matrix
      f_spin = 2._dp

      ALLOCATE (tmp_fm, tmp_fm_momo(3), tmp_fm_like_mos(3))
      ALLOCATE (diagonal_elements(nmo))

      ! Allocate temporary matrices
      CALL cp_fm_create(tmp_fm%matrix, dcdr_env%likemos_fm_struct)
      DO i = 1, 3
         CALL cp_fm_create(tmp_fm_momo(i)%matrix, dcdr_env%momo_fm_struct)
         CALL cp_fm_create(tmp_fm_like_mos(i)%matrix, dcdr_env%likemos_fm_struct)
      END DO

      ! Build the full coefficient derivatives.
      CALL cp_fm_set_all(tmp_fm%matrix, 0.0_dp)
      CALL cp_fm_scale_and_add(0._dp, dcdr_env%dCR_prime(1)%matrix, 1._dp, dcdr_env%dCR(1)%matrix)
      CALL cp_dbcsr_sm_fm_multiply(dcdr_env%matrix_s1(dcdr_env%beta + 1)%matrix, mo_coeff, &
                                   tmp_fm%matrix, ncol=nmo)
      CALL cp_fm_gemm("T", "N", nmo, nmo, nao, &
                      1.0_dp, mo_coeff, tmp_fm%matrix, &
                      0.0_dp, tmp_fm_momo(1)%matrix)

      !   C^1 <- -dCR - 0.5 * mo_coeff @ S1_ij
      CALL cp_fm_gemm("N", "N", nao, nmo, nmo, &
                      -0.5_dp, mo_coeff, tmp_fm_momo(1)%matrix, &
                      -1.0_dp, dcdr_env%dCR_prime(1)%matrix)

      ! FIRST CONTRIBUTION: dCR * moments * mo
      this_factor = -2._dp*f_spin
      DO alpha = 1, 3
         DO icenter = 1, dcdr_env%nbr_center(1)
            CALL build_local_moment_matrix(qs_env, dcdr_env%moments, 1, &
                                           ref_point=centers_set(1)%array(1:3, icenter))
            CALL multiply_localization(ao_matrix=dcdr_env%moments(alpha)%matrix, &
                                       mo_coeff=dcdr_env%dCR_prime(1)%matrix, work=tmp_fm%matrix, nmo=nmo, &
                                       icenter=icenter, &
                                       res=tmp_fm_like_mos(alpha)%matrix)
            CALL dbcsr_set(dcdr_env%moments(alpha)%matrix, 0.0_dp)
         END DO

         CALL cp_fm_gemm("T", "N", nmo, nmo, nao, &
                         1.0_dp, mo_coeff, tmp_fm_like_mos(alpha)%matrix, &
                         0.0_dp, tmp_fm_momo(alpha)%matrix)
         CALL cp_fm_get_diag(tmp_fm_momo(alpha)%matrix, diagonal_elements)

         DO icenter = 1, dcdr_env%nbr_center(1)
            map_atom = mapping_wannier_atom(icenter)
            map_molecule = mapping_atom_molecule(map_atom)
            tmp_aptcontr = this_factor*diagonal_elements(icenter)

            apt_subset(dcdr_env%beta, alpha, dcdr_env%lambda, map_molecule) &
               = apt_subset(dcdr_env%beta, alpha, dcdr_env%lambda, map_molecule) + tmp_aptcontr

            apt_center(dcdr_env%beta, alpha, dcdr_env%lambda, icenter) &
               = apt_center(dcdr_env%beta, alpha, dcdr_env%lambda, icenter) + tmp_aptcontr
         END DO

         apt_coeff_derivative = this_factor*SUM(diagonal_elements)
         apt_el(dcdr_env%beta, alpha, dcdr_env%lambda) &
            = apt_el(dcdr_env%beta, alpha, dcdr_env%lambda) + apt_coeff_derivative
      END DO

      ! SECOND CONTRIBUTION: We assemble all combinations of r_i, dphi/d(idir)
      ! build part with AOs differentiated with respect to nuclear coordinates
      ! difdip contains derivatives with respect to atom dcdr_env%lambda
      ! difdip(alpha, beta): < a | r_alpha | d b/dR_beta >
      this_factor = -f_spin
      DO alpha = 1, 3
         DO icenter = 1, dcdr_env%nbr_center(1)
            ! Build the AO matrix with the right wannier center as reference point
            CALL dbcsr_set(dcdr_env%matrix_difdip(1, dcdr_env%beta)%matrix, 0._dp)
            CALL dbcsr_set(dcdr_env%matrix_difdip(2, dcdr_env%beta)%matrix, 0._dp)
            CALL dbcsr_set(dcdr_env%matrix_difdip(3, dcdr_env%beta)%matrix, 0._dp)
            CALL dipole_deriv_ao(qs_env, dcdr_env%matrix_difdip, dcdr_env%delta_basis_function, &
                                 1, centers_set(1)%array(1:3, icenter))
            CALL multiply_localization(ao_matrix=dcdr_env%matrix_difdip(alpha, dcdr_env%beta)%matrix, &
                                       mo_coeff=mo_coeff, work=tmp_fm%matrix, nmo=nmo, &
                                       icenter=icenter, &
                                       res=tmp_fm_like_mos(alpha)%matrix)
         END DO ! icenter

         CALL cp_fm_gemm("T", "N", nmo, nmo, nao, &
                         1.0_dp, mo_coeff, tmp_fm_like_mos(alpha)%matrix, &
                         0.0_dp, tmp_fm_momo(alpha)%matrix)
         CALL cp_fm_get_diag(tmp_fm_momo(alpha)%matrix, diagonal_elements)

         DO icenter = 1, dcdr_env%nbr_center(1)
            map_atom = mapping_wannier_atom(icenter)
            map_molecule = mapping_atom_molecule(map_atom)
            tmp_aptcontr = this_factor*diagonal_elements(icenter)

            apt_subset(dcdr_env%beta, alpha, dcdr_env%lambda, map_molecule) &
               = apt_subset(dcdr_env%beta, alpha, dcdr_env%lambda, map_molecule) + tmp_aptcontr

            apt_center(dcdr_env%beta, alpha, dcdr_env%lambda, icenter) &
               = apt_center(dcdr_env%beta, alpha, dcdr_env%lambda, icenter) + tmp_aptcontr
         END DO

         ! The negative sign is due to dipole_deriv_ao computing the derivatives with respect to nuclear coordinates.
         apt_basis_derivative = this_factor*SUM(diagonal_elements)

         apt_el(dcdr_env%beta, alpha, dcdr_env%lambda) &
            = apt_el(dcdr_env%beta, alpha, dcdr_env%lambda) + apt_basis_derivative

      END DO  ! alpha

      ! Finally the nuclear contribution: nuclear charge * Kronecker_delta_{dcdr_env%beta,i}
      CALL get_atomic_kind(particle_set(dcdr_env%lambda)%atomic_kind, kind_number=ikind)
      CALL get_qs_kind(qs_kind_set(ikind), core_charge=charge, ghost=ghost)
      IF (.NOT. ghost) THEN  ! Those come from the pseudopotential, right?
         apt_nuc(dcdr_env%beta, dcdr_env%beta, dcdr_env%lambda) = &
            apt_nuc(dcdr_env%beta, dcdr_env%beta, dcdr_env%lambda) + charge

         map_molecule = mapping_atom_molecule(dcdr_env%lambda)
         apt_subset(dcdr_env%beta, dcdr_env%beta, dcdr_env%lambda, map_molecule) &
            = apt_subset(dcdr_env%beta, dcdr_env%beta, dcdr_env%lambda, map_molecule) + charge
      END IF

      ! And deallocate all the things!
      CALL cp_fm_release(tmp_fm%matrix)
      DO i = 1, 3
         CALL cp_fm_release(tmp_fm_like_mos(i)%matrix)
         CALL cp_fm_release(tmp_fm_momo(i)%matrix)
      END DO
      DEALLOCATE (tmp_fm, tmp_fm_like_mos, tmp_fm_momo)
      DEALLOCATE (diagonal_elements)

      CALL timestop(handle)
   END SUBROUTINE apt_dR_localization

END MODULE qs_dcdr
