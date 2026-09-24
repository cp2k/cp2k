!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief   Add the DFT+U contribution to the Hamiltonian matrix
!> \details The implemented methods refers to:\n
!>          S. L. Dudarev, D. Nguyen Manh, and A. P. Sutton,
!>          Philos. Mag. B \b 75, 613 (1997)\n
!>          S. L. Dudarev et al.,
!>          Phys. Rev. B \b 57, 1505 (1998)
!> \author  Matthias Krack (MK)
!> \date    14.01.2008
!> \version 1.0
! **************************************************************************************************
MODULE dft_plus_u
   USE OMP_LIB,                         ONLY: omp_get_max_threads,&
                                              omp_get_thread_num
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE bibliography,                    ONLY: Dudarev1997,&
                                              Dudarev1998,&
                                              cite_reference
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_column_scale,&
                                              cp_cfm_transpose
   USE cp_cfm_diag,                     ONLY: cp_cfm_heevd
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_release,&
                                              cp_cfm_to_cfm,&
                                              cp_cfm_type,&
                                              cp_fm_to_cfm
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_copy, dbcsr_create, dbcsr_deallocate_matrix, dbcsr_desymmetrize, dbcsr_finalize, &
        dbcsr_get_block_p, dbcsr_init_p, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_multiply, &
        dbcsr_p_type, dbcsr_put_block, dbcsr_release, dbcsr_set, dbcsr_type
   USE cp_dbcsr_contrib,                ONLY: dbcsr_get_block_diag,&
                                              dbcsr_get_diag,&
                                              dbcsr_set_diag
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              cp_dbcsr_plus_fm_fm_t,&
                                              cp_dbcsr_sm_fm_multiply
   USE cp_dbcsr_output,                 ONLY: cp_dbcsr_write_sparse_matrix,&
                                              write_fm_with_basis_info
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add,&
                                              cp_fm_schur_product,&
                                              cp_fm_transpose
   USE cp_fm_diag,                      ONLY: choose_eigv_solver
   USE cp_fm_struct,                    ONLY: cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_set_submatrix,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr,&
                                              low_print_level
   USE input_constants,                 ONLY: diag_update_method_adiis,&
                                              plus_u_lowdin,&
                                              plus_u_mulliken,&
                                              plus_u_mulliken_charges,&
                                              plus_u_tensorial
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE kpoint_types,                    ONLY: kpoint_env_type,&
                                              kpoint_type
   USE local_gemm_api,                  ONLY: LOCAL_GEMM_PU_GPU,&
                                              local_gemm_ctxt_type
   USE mathconstants,                   ONLY: z_one,&
                                              z_zero
   USE mathlib,                         ONLY: invert_matrix,&
                                              jacobi
   USE message_passing,                 ONLY: mp_para_env_type
   USE orbital_symbols,                 ONLY: sgf_symbol
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: evolt
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type,&
                                              set_qs_kind
   USE qs_kpoint_operators,             ONLY: kpoint_operator_context_create,&
                                              kpoint_operator_context_release,&
                                              kpoint_operator_context_type,&
                                              kpoint_operator_get
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_scf_types,                    ONLY: qs_scf_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dft_plus_u'

   PUBLIC :: plus_u
   PUBLIC :: lowdin_kp_init, lowdin_kp_initial_density, lowdin_kp_population, lowdin_kp_apply
   PUBLIC :: lowdin_kp_forces
   PUBLIC :: lowdin_kp_apply_mo, lowdin_kp_kernel
   PUBLIC :: lowdin_kp_prepare

CONTAINS
! **************************************************************************************************
!> \brief         Add the DFT+U contribution to the Hamiltonian matrix.\n
!>                Wrapper routine for all "+U" methods
!> \param[in]     qs_env Quickstep environment
!> \param[in,out] matrix_h Hamiltonian matrices for each spin
!> \param[in,out] matrix_w Energy weighted density matrices for each spin
!> \date          14.01.2008
!> \author        Matthias Krack (MK)
!> \version       1.0
! **************************************************************************************************
   SUBROUTINE plus_u(qs_env, matrix_h, matrix_w)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), OPTIONAL, &
         POINTER                                         :: matrix_h, matrix_w

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'plus_u'

      INTEGER                                            :: handle, output_unit, print_level
      LOGICAL                                            :: orthonormal_basis, should_output
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(section_vals_type), POINTER                   :: input

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(qs_env))

      NULLIFY (input, dft_control)

      logger => cp_get_default_logger()

      CALL get_qs_env(qs_env=qs_env, &
                      input=input, &
                      dft_control=dft_control)

      CALL cite_reference(Dudarev1997)
      CALL cite_reference(Dudarev1998)

      ! Later we could save here some time, if the method in use has this property
      ! which then has to be figured out here

      orthonormal_basis = .FALSE.

      ! Setup print control

      print_level = logger%iter_info%print_level
      should_output = (BTEST(cp_print_key_should_output(logger%iter_info, input, &
                                                        "DFT%PRINT%PLUS_U"), cp_p_file) .AND. &
                       (.NOT. PRESENT(matrix_w)))
      output_unit = cp_print_key_unit_nr(logger, input, "DFT%PRINT%PLUS_U", &
                                         extension=".plus_u", &
                                         ignore_should_output=should_output, &
                                         log_filename=.FALSE.)

      ! Select DFT+U method

      SELECT CASE (dft_control%plus_u_method_id)
      CASE (plus_u_lowdin)
         IF (orthonormal_basis) THEN
            ! For an orthonormal basis the Lowdin method and the Mulliken method
            ! are equivalent
            CALL mulliken(qs_env, orthonormal_basis, matrix_h, &
                          should_output, output_unit, print_level)
         ELSE
            CALL lowdin(qs_env, matrix_h, matrix_w, &
                        should_output, output_unit, print_level)
         END IF
      CASE (plus_u_mulliken)
         CALL mulliken(qs_env, orthonormal_basis, matrix_h, &
                       should_output, output_unit, print_level)
      CASE (plus_u_mulliken_charges)
         CALL mulliken_charges(qs_env, orthonormal_basis, matrix_h, matrix_w, &
                               should_output, output_unit, print_level)
      CASE (plus_u_tensorial)
         CALL tensorial(qs_env, matrix_h, matrix_w, &
                        should_output, output_unit, print_level)
      CASE DEFAULT
         CPABORT("Invalid DFT+U method requested")
      END SELECT

      CALL cp_print_key_finished_output(output_unit, logger, input, "DFT%PRINT%PLUS_U", &
                                        ignore_should_output=should_output)

      CALL timestop(handle)

   END SUBROUTINE plus_u

! **************************************************************************************************
!> \brief         Add a DFT+U contribution to the Hamiltonian matrix\n
!>                using a method based on Lowdin charges
!>                \f[Q = S^{1/2} P S^{1/2}\f]
!>                where \b P and \b S are the density and the
!>                overlap matrix, respectively.
!> \param[in]     qs_env Quickstep environment
!> \param[in,out] matrix_h Hamiltonian matrices for each spin
!> \param[in,out] matrix_w Energy weighted density matrices for each spin
!> \param should_output ...
!> \param output_unit ...
!> \param print_level ...
!> \date          02.07.2008
!> \par
!>  \f{eqnarray*}{
!>   E^{\rm DFT+U} & = & E^{\rm DFT} +  E^{\rm U}
!>                 & = & E^{\rm DFT} +  \frac{1}{2}(U - J)\sum_\mu (q_\mu - q_\mu^2)\\[1ex]
!>   V_{\mu\nu}^{\rm DFT+U} & = & V_{\mu\nu}^{\rm DFT} + V_{\mu\nu}^{\rm U}\\\
!>                          & = & \frac{\partial E^{\rm DFT}}
!>                                     {\partial P_{\mu\nu}} +
!>                                \frac{\partial E^{\rm U}}
!>                                     {\partial P_{\mu\nu}}\\\
!>                          & = & H_{\mu\nu} +
!>                                \frac{\partial E^{\rm U}}{\partial q_\mu}
!>                                \frac{\partial q_\mu}{\partial P_{\mu\nu}}\\\
!>  \f}
!> \author        Matthias Krack (MK)
!> \version       1.0
! **************************************************************************************************
   SUBROUTINE lowdin(qs_env, matrix_h, matrix_w, should_output, output_unit, &
                     print_level)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), OPTIONAL, &
         POINTER                                         :: matrix_h, matrix_w
      LOGICAL, INTENT(IN)                                :: should_output
      INTEGER, INTENT(IN)                                :: output_unit, print_level

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'lowdin'

      CHARACTER(LEN=10)                                  :: spin_info
      CHARACTER(LEN=6), ALLOCATABLE, DIMENSION(:)        :: symbol
      CHARACTER(LEN=default_string_length)               :: atomic_kind_name
      INTEGER :: atom_a, handle, i, i0, iatom, ikind, iorb, isb, iset, isgf, ishell, ispin, j, &
         jsb, jsgf, lu, m, max_scf, n, natom_of_kind, nkind, norb, nsb, nsbsize, nset, nsgf, &
         nsgf_kind, nspin
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: shell_sgf
      INTEGER, DIMENSION(1)                              :: iloc
      INTEGER, DIMENSION(:), POINTER                     :: atom_list, nshell, orbitals
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf, l, last_sgf
      LOGICAL                                            :: debug, dft_plus_u_atom, do_kpoints, &
                                                            found, if_plus_j, just_energy, smear
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: orb_occ
      REAL(KIND=dp) :: eps_scf, eps_u_ramping, fspin, hund_j, occ, perturbation_strength, sij, &
         trq, trq2, trqxq_hund_j, u_minus_j, u_minus_j_target, u_ramping
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigval, q_eigval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: orbq, q_eigvec, q_matrix, q_opp_matrix, &
                                                            q_work, slam
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: local_data
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: q_block, q_opp_block, v_block, vhxc_block
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_fm_struct_type), POINTER                   :: fmstruct
      TYPE(cp_fm_type)                                   :: fm_sev, fm_work1, fm_work2, slambda
      TYPE(cp_fm_type), POINTER                          :: fm_s_half
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_vhxc
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_p, matrix_s
      TYPE(dbcsr_type)                                   :: sm_q, sm_q_opp, sm_v
      TYPE(dbcsr_type), POINTER                          :: sm_h, sm_p, sm_p_opp, sm_s, sm_vhxc, sm_w
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_scf_env_type), POINTER                     :: scf_env

      CALL timeset(routineN, handle)

      debug = .FALSE. ! Set to .TRUE. to print debug information

      NULLIFY (sm_h, sm_p, sm_p_opp, sm_s, sm_w, vhxc_block, sm_vhxc, matrix_vhxc)

      smear = .FALSE.
      max_scf = -1
      eps_scf = 1.0E30_dp
      if_plus_j = .FALSE.

      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      dft_control=dft_control, &
                      do_kpoints=do_kpoints, &
                      kpoints=kpoints, &
                      energy=energy, &
                      matrix_s_kp=matrix_s, &
                      matrix_vhxc=matrix_vhxc, &
                      rho=rho, &
                      scf_env=scf_env, &
                      para_env=para_env)

      CALL qs_rho_get(rho, rho_ao_kp=matrix_p) ! Density matrices in sparse format

      energy%dft_plus_u = 0.0_dp

      nspin = dft_control%nspins
      IF (dft_control%mtlr_dft_with_perturbation) THEN
         IF (.NOT. ASSOCIATED(matrix_vhxc)) THEN
            CPABORT("MTLR requires the projected Hxc matrix for every spin channel.")
         END IF
         IF (SIZE(matrix_vhxc) /= nspin) THEN
            CPABORT("The number of projected Hxc matrices does not match the spin channels.")
         END IF
      END IF

      IF (nspin == 2) THEN
         fspin = 1.0_dp
      ELSE
         fspin = 0.5_dp
      END IF

      ! Get the number of basis functions and atomic kinds.
      CALL get_qs_kind_set(qs_kind_set, nsgf=nsgf)

      nkind = SIZE(atomic_kind_set)

      IF (PRESENT(matrix_h) .OR. PRESENT(matrix_w)) THEN
         just_energy = .FALSE.
      ELSE
         just_energy = .TRUE.
      END IF

      ! Check if unlike spin +J correction is needed and reject an
      ! incompatible enforced occupation setup before any expensive work.
      DO ikind = 1, nkind
         NULLIFY (orbitals)
         CALL get_qs_kind(qs_kind_set(ikind), &
                          hund_j=hund_j, &
                          orbitals=orbitals)
         IF (hund_j /= 0.0_dp) THEN
            if_plus_j = .TRUE.
            IF (ASSOCIATED(orbitals)) THEN
               CALL cp_abort(__LOCATION__, "ENFORCE_OCCUPATION is incompatible with a nonzero Hund J "// &
                             "in the Lowdin DFT+U+J method.")
            END IF
         END IF
      END DO

      NULLIFY (orbitals)

      IF (.NOT. do_kpoints) THEN
         ! Retrieve S^(1/2) from the SCF environment
         fm_s_half => scf_env%s_half
         CPASSERT(ASSOCIATED(fm_s_half))
         ! work matrices
         CALL cp_fm_get_info(fm_s_half, matrix_struct=fmstruct)
         CALL cp_fm_create(matrix=fm_work1, matrix_struct=fmstruct, &
                           name="FULL WORK MATRIX 1")
         CALL cp_fm_create(matrix=fm_work2, matrix_struct=fmstruct, &
                           name="FULL WORK MATRIX 2")
      END IF

      ! Calculate S eigenvectors and Lambda matrix for forces
      ! See sTDA forces (get_lowdin_mo_coefficients in qs_tddfpt2_stda_utils
      ! A. Hehn et al JCTC 2022, 18, 4186
      IF (PRESENT(matrix_w) .AND. .NOT. do_kpoints) THEN
         CALL cp_fm_create(matrix=fm_sev, matrix_struct=fmstruct)
         CALL cp_fm_create(matrix=slambda, matrix_struct=fmstruct)
         ALLOCATE (eigval(nsgf), slam(nsgf, 1))
         sm_s => matrix_s(1, 1)%matrix
         CALL copy_dbcsr_to_fm(sm_s, fm_work1)
         CALL choose_eigv_solver(fm_work1, fm_sev, eigval)
         !
         DO i = 1, nsgf
            IF (eigval(i) > 0._dp) THEN
               slam(i, 1) = SQRT(eigval(i))
            ELSE
               CPABORT("S matrix not positive definit")
            END IF
         END DO
         DO i = 1, nsgf
            CALL cp_fm_set_submatrix(slambda, slam, 1, i, nsgf, 1, 1.0_dp, 0.0_dp)
         END DO
         DO i = 1, nsgf
            CALL cp_fm_set_submatrix(slambda, slam, i, 1, 1, nsgf, 1.0_dp, 1.0_dp, .TRUE.)
         END DO
         CALL cp_fm_get_info(slambda, local_data=local_data)
         DO i = 1, SIZE(local_data, 2)
            DO j = 1, SIZE(local_data, 1)
               sij = local_data(j, i)
               IF (sij > 0.0_dp) sij = 1.0_dp/sij
               local_data(j, i) = sij
            END DO
         END DO
         DEALLOCATE (eigval, slam)
      END IF

      ! Calculate S^(1/2)*P*S^(1/2)
      IF (do_kpoints) THEN
         ! OT and ADIIS evaluate the current trial density, including energy-only line searches.
         ! Density mixing instead owns separate accepted coordinates.
         IF (qs_env%scf_control%use_ot .OR. &
             qs_env%scf_control%diagonalization%update_method == diag_update_method_adiis .OR. PRESENT(matrix_w)) THEN
            CALL lowdin_kp_population(kpoints, alpha=1.0_dp)
         END IF
         ALLOCATE (orbq(nsgf, nspin))
         orbq = 0.0_dp
         orbq(kpoints%lowdin_orbitals, :) = kpoints%lowdin_q
      END IF

      ! Create local block diagonal matrices
      sm_s => matrix_s(1, 1)%matrix
      CALL dbcsr_get_block_diag(sm_s, sm_q)
      IF (if_plus_j .AND. nspin == 2) THEN
         CALL dbcsr_get_block_diag(sm_s, sm_q_opp)
      END IF
      CALL dbcsr_get_block_diag(sm_s, sm_v)

      ! Loop over all spins
      DO ispin = 1, nspin

         CALL dbcsr_set(sm_q, 0.0_dp)
         IF (if_plus_j .AND. nspin == 2) THEN
            CALL dbcsr_set(sm_q_opp, 0.0_dp)
         END IF
         CALL dbcsr_set(sm_v, 0.0_dp)

         IF (dft_control%mtlr_dft_with_perturbation) THEN
            IF (ispin == 1) perturbation_strength = 0.9_dp*dft_control%perturbation_strength
            IF (ispin == 2) perturbation_strength = 1.1_dp*dft_control%perturbation_strength
            IF (.NOT. ASSOCIATED(matrix_vhxc(ispin)%matrix)) THEN
               CPABORT("MTLR projected Hxc matrix is not initialized.")
            END IF
            sm_vhxc => matrix_vhxc(ispin)%matrix
         END IF

         IF (do_kpoints) THEN
            CALL dbcsr_set_diag(sm_q, orbq(:, ispin))
            IF (if_plus_j .AND. nspin == 2) CALL dbcsr_set_diag(sm_q_opp, orbq(:, 3 - ispin))
         ELSE
            ! Calculate S^(1/2)*P*S^(1/2) as a full matrix (Lowdin)
            sm_p => matrix_p(ispin, 1)%matrix
            CALL cp_dbcsr_sm_fm_multiply(sm_p, fm_s_half, fm_work1, nsgf)
            CALL parallel_gemm(transa="N", &
                               transb="N", &
                               m=nsgf, &
                               n=nsgf, &
                               k=nsgf, &
                               alpha=1.0_dp, &
                               matrix_a=fm_s_half, &
                               matrix_b=fm_work1, &
                               beta=0.0_dp, &
                               matrix_c=fm_work2)
            IF (debug) THEN
               CALL cp_dbcsr_write_sparse_matrix(sm_p, 4, 6, qs_env, para_env, &
                                                 output_unit=output_unit)
               CALL write_fm_with_basis_info(fm_s_half, 4, 6, qs_env, para_env, &
                                             output_unit=output_unit)
               CALL write_fm_with_basis_info(fm_work2, 4, 6, qs_env, para_env, &
                                             output_unit=output_unit)
            END IF ! debug
            ! Copy occupation matrix to sparse matrix format, finally we are only
            ! interested in the diagonal (atomic) blocks, i.e. the previous full
            ! matrix product is not the most efficient choice, anyway.
            CALL copy_fm_to_dbcsr(fm_work2, sm_q, keep_sparsity=.TRUE.)

            ! If nspin == 1 and one would like to add "unlike" spin +J corrections,
            ! sm_p is already enough and sm_p_opp is not needed
            IF (if_plus_j .AND. nspin == 2) THEN
               ! Calculate S^(1/2)*P*S^(1/2) for the other spin channel if
               ! unlike spin +J correction is activated.
               sm_p_opp => matrix_p(3 - ispin, 1)%matrix
               CALL cp_dbcsr_sm_fm_multiply(sm_p_opp, fm_s_half, fm_work1, nsgf)
               CALL parallel_gemm(transa="N", &
                                  transb="N", &
                                  m=nsgf, &
                                  n=nsgf, &
                                  k=nsgf, &
                                  alpha=1.0_dp, &
                                  matrix_a=fm_s_half, &
                                  matrix_b=fm_work1, &
                                  beta=0.0_dp, &
                                  matrix_c=fm_work2)
               CALL copy_fm_to_dbcsr(fm_work2, sm_q_opp, keep_sparsity=.TRUE.)
            END IF
         END IF

         ! E[DFT+U] = E[DFT] + E[U]
         !          = E[DFT] + (U - J)*(Tr(q) - Tr(q*q))/2

         ! V(i,j)[DFT+U] = V(i,j)[DFT] + V(i,j)[U]
         !               = dE[DFT]/dP(i,j) + dE[U]/dP(i,j)
         !               = dE[DFT]/dP(i,j) + (dE(U)/dq)*(dq/dP(i,j))

         ! Loop over all atomic kinds
         DO ikind = 1, nkind

            ! Load the required atomic kind data
            CALL get_atomic_kind(atomic_kind_set(ikind), &
                                 atom_list=atom_list, &
                                 name=atomic_kind_name, &
                                 natom=natom_of_kind)

            CALL get_qs_kind(qs_kind_set(ikind), &
                             dft_plus_u_atom=dft_plus_u_atom, &
                             l_of_dft_plus_u=lu, &
                             nsgf=nsgf_kind, &
                             basis_set=orb_basis_set, &
                             u_minus_j=u_minus_j, &
                             hund_j=hund_j, &
                             u_minus_j_target=u_minus_j_target, &
                             u_ramping=u_ramping, &
                             eps_u_ramping=eps_u_ramping, &
                             orbitals=orbitals, &
                             eps_scf=eps_scf, &
                             max_scf=max_scf, &
                             smear=smear)

            ! Check, if the atoms of this atomic kind need a DFT+U correction
            IF (.NOT. ASSOCIATED(orb_basis_set)) CYCLE
            IF (.NOT. dft_plus_u_atom) CYCLE
            IF (lu < 0) CYCLE

            ! Apply U ramping if requested
            IF ((ispin == 1) .AND. (u_ramping > 0.0_dp) .AND. .NOT. do_kpoints) THEN
               IF (qs_env%scf_env%iter_delta <= eps_u_ramping) THEN
                  u_minus_j = MIN(u_minus_j + u_ramping, u_minus_j_target)
                  CALL set_qs_kind(qs_kind_set(ikind), u_minus_j=u_minus_j)
               END IF
               IF (should_output .AND. (output_unit > 0)) THEN
                  WRITE (UNIT=output_unit, FMT="(T3,A,3X,A,F0.3,A)") &
                     "Kind name: "//TRIM(ADJUSTL(atomic_kind_name)), &
                     "U(eff) = ", u_minus_j*evolt, " eV"
               END IF
            END IF

            IF (u_minus_j == 0.0_dp .AND. hund_j == 0.0_dp) CYCLE

            ! Load the required Gaussian basis set data
            CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                   first_sgf=first_sgf, &
                                   l=l, &
                                   last_sgf=last_sgf, &
                                   nset=nset, &
                                   nshell=nshell)

            ! Select this kind's shell once; occupations and potentials use the same AO order.
            ALLOCATE (shell_sgf(nsgf_kind))
            n = 0
            DO iset = 1, nset
               DO ishell = 1, nshell(iset)
                  IF (l(ishell, iset) /= lu) CYCLE
                  DO isgf = first_sgf(ishell, iset), last_sgf(ishell, iset)
                     n = n + 1
                     shell_sgf(n) = isgf
                  END DO
               END DO
            END DO
            nsbsize = 2*lu + 1
            nsb = n/nsbsize

            ALLOCATE (q_matrix(n, n))
            q_matrix(:, :) = 0.0_dp

            IF (hund_j /= 0.0_dp .AND. nspin == 2) THEN
               ALLOCATE (q_opp_matrix(n, n))
               q_opp_matrix(:, :) = 0.0_dp
            END IF

            IF (dft_control%mtlr_dft_with_perturbation .AND. &
                dft_control%mtlr_ikind == ikind) THEN
               dft_control%vhxc(ispin) = 0.0_dp
               dft_control%trq(ispin) = 0.0_dp
            END IF

            ! Print headline if requested
            IF (should_output .AND. (print_level > low_print_level)) THEN
               IF (output_unit > 0) THEN
                  ALLOCATE (symbol(nsbsize))
                  DO m = -lu, lu
                     symbol(lu + m + 1) = sgf_symbol(0, lu, m)
                  END DO
                  IF (nspin > 1) THEN
                     WRITE (UNIT=spin_info, FMT="(A8,I2)") " of spin", ispin
                  ELSE
                     spin_info = ""
                  END IF
                  WRITE (UNIT=output_unit, FMT="(/,T3,A,I0,A,/,/,T5,A,10(2X,A6))") &
                     "DFT+U occupations"//TRIM(spin_info)//" for the atoms of atomic kind ", ikind, &
                     ": "//TRIM(atomic_kind_name), &
                     "Atom   Shell  ", (ADJUSTR(symbol(i)), i=1, nsbsize), " Trace"
                  DEALLOCATE (symbol)
               END IF
            END IF

            ! Loop over all atoms of the current atomic kind
            DO iatom = 1, natom_of_kind
               atom_a = atom_list(iatom)
               q_matrix(:, :) = 0.0_dp
               IF (ALLOCATED(q_opp_matrix)) THEN
                  q_opp_matrix(:, :) = 0.0_dp
               END IF

               ! Get diagonal block
               CALL dbcsr_get_block_p(matrix=sm_q, &
                                      row=atom_a, &
                                      col=atom_a, &
                                      block=q_block, &
                                      found=found)
               IF (hund_j /= 0.0_dp .AND. nspin == 2) THEN
                  CALL dbcsr_get_block_p(matrix=sm_q_opp, &
                                         row=atom_a, &
                                         col=atom_a, &
                                         block=q_opp_block, &
                                         found=found)
               END IF
               IF (dft_control%mtlr_dft_with_perturbation .AND. &
                   dft_control%mtlr_ikind == ikind .AND. &
                   qs_kind_set(ikind)%dft_plus_u%lr_atom == atom_a) THEN
                  CALL dbcsr_get_block_p(matrix=sm_vhxc, &
                                         row=atom_a, &
                                         col=atom_a, &
                                         block=vhxc_block, &
                                         found=found)
               END IF

               IF (ASSOCIATED(q_block)) THEN
                  ! Calculate energy contribution to E(U)
                  DO i = 1, n
                     isgf = shell_sgf(i)
                     q_matrix(i, i) = q_block(isgf, isgf)
                     IF (hund_j /= 0.0_dp .AND. nspin == 2 .AND. ASSOCIATED(q_opp_block)) THEN
                        q_opp_matrix(i, i) = q_opp_block(isgf, isgf)
                     END IF
                     IF (ASSOCIATED(vhxc_block) .AND. dft_control%mtlr_dft_with_perturbation .AND. &
                         dft_control%mtlr_ikind == ikind .AND. qs_kind_set(ikind)%dft_plus_u%lr_atom == atom_a) THEN
                        dft_control%vhxc(ispin) = dft_control%vhxc(ispin) + vhxc_block(isgf, isgf)
                        dft_control%trq(ispin) = dft_control%trq(ispin) + q_block(isgf, isgf)
                     END IF
                  END DO

                  ! Perform the requested manipulations of the (initial) orbital occupations
                  IF (ASSOCIATED(orbitals)) THEN
                     IF ((qs_env%scf_env%iter_delta >= eps_scf) .OR. &
                         ((qs_env%scf_env%outer_scf%iter_count == 0) .AND. &
                          (qs_env%scf_env%iter_count <= max_scf))) THEN
                        ALLOCATE (orb_occ(nsbsize))
                        ALLOCATE (q_eigval(n))
                        q_eigval(:) = 0.0_dp
                        ALLOCATE (q_eigvec(n, n))
                        q_eigvec(:, :) = 0.0_dp
                        norb = SIZE(orbitals)
                        CALL jacobi(q_matrix, q_eigval, q_eigvec)
                        q_matrix(:, :) = 0.0_dp
                        DO isb = 1, nsb
                           trq = 0.0_dp
                           DO i = (isb - 1)*nsbsize + 1, isb*nsbsize
                              trq = trq + q_eigval(i)
                           END DO
                           IF (smear) THEN
                              occ = trq/REAL(norb, KIND=dp)
                           ELSE
                              occ = 1.0_dp/fspin
                           END IF
                           orb_occ(:) = .FALSE.
                           iloc = MAXLOC(q_eigvec(:, isb*nsbsize))
                           jsb = INT((iloc(1) - 1)/nsbsize) + 1
                           i = 0
                           i0 = (jsb - 1)*nsbsize + 1
                           iorb = -1000
                           DO j = i0, jsb*nsbsize
                              i = i + 1
                              IF (i > norb) THEN
                                 DO m = -lu, lu
                                    IF (.NOT. orb_occ(lu + m + 1)) THEN
                                       iorb = i0 + lu + m
                                       orb_occ(lu + m + 1) = .TRUE.
                                    END IF
                                 END DO
                              ELSE
                                 iorb = i0 + lu + orbitals(i)
                                 orb_occ(lu + orbitals(i) + 1) = .TRUE.
                              END IF
                              CPASSERT(iorb /= -1000)
                              iloc = MAXLOC(q_eigvec(iorb, :))
                              q_eigval(iloc(1)) = MIN(occ, trq)
                              q_matrix(:, iloc(1)) = q_eigval(iloc(1))*q_eigvec(:, iloc(1)) ! backtransform left
                              trq = trq - q_eigval(iloc(1))
                           END DO
                        END DO
                        q_matrix(:, :) = MATMUL(q_matrix, TRANSPOSE(q_eigvec)) ! backtransform right
                        DEALLOCATE (orb_occ)
                        DEALLOCATE (q_eigval)
                        DEALLOCATE (q_eigvec)
                     END IF
                  END IF ! orbitals associated

                  trq = fspin*SUM([(q_matrix(i, i), i=1, n)])
                  trq2 = fspin*fspin*SUM(q_matrix*TRANSPOSE(q_matrix))
                  trqxq_hund_j = 0.0_dp
                  IF (hund_j /= 0.0_dp) THEN
                     IF (nspin == 2) THEN
                        trqxq_hund_j = SUM(q_matrix*TRANSPOSE(q_opp_matrix))
                     ELSE
                        trqxq_hund_j = SUM(q_matrix*TRANSPOSE(q_matrix))
                     END IF
                  END IF

                  ! Calculate energy contribution to E(U)
                  energy%dft_plus_u = energy%dft_plus_u + 0.5_dp*u_minus_j*(trq - trq2)/fspin
                  IF (hund_j /= 0.0_dp) THEN
                     trqxq_hund_j = fspin*fspin*trqxq_hund_j
                     energy%dft_plus_u = energy%dft_plus_u + 0.5_dp*hund_j*trqxq_hund_j/fspin
                  END IF

                  ! Calculate potential V(U) = dE(U)/dq
                  IF (.NOT. just_energy) THEN
                     CALL dbcsr_get_block_p(matrix=sm_v, &
                                            row=atom_a, &
                                            col=atom_a, &
                                            block=v_block, &
                                            found=found)
                     CPASSERT(ASSOCIATED(v_block))

                     DO j = 1, n
                        jsgf = shell_sgf(j)
                        DO i = 1, n
                           isgf = shell_sgf(i)
                           IF (i == j) THEN
                              v_block(isgf, isgf) = u_minus_j*(0.5_dp - fspin*q_matrix(i, i))
                              IF (dft_control%mtlr_dft_with_perturbation .AND. &
                                  dft_control%mtlr_ikind == ikind .AND. &
                                  qs_kind_set(ikind)%dft_plus_u%lr_atom == atom_a) THEN
                                 v_block(isgf, isgf) = v_block(isgf, isgf) + perturbation_strength
                              END IF
                           ELSE
                              CPASSERT(ABS(q_matrix(j, i)) < 1.0E-14_dp)
                              IF (hund_j /= 0.0_dp .AND. nspin == 2) THEN
                                 CPASSERT(ABS(q_opp_matrix(j, i)) < 1.0E-14_dp)
                              END IF
                              v_block(isgf, jsgf) = -u_minus_j*fspin*q_matrix(j, i)
                           END IF
                           IF (hund_j /= 0.0_dp) THEN
                              IF (nspin == 2) THEN
                                 v_block(isgf, jsgf) = v_block(isgf, jsgf) + hund_j*fspin*q_opp_matrix(j, i)
                              ELSE
                                 v_block(isgf, jsgf) = v_block(isgf, jsgf) + hund_j*fspin*q_matrix(j, i)
                              END IF
                           END IF
                        END DO
                     END DO
                  END IF ! not just energy

               END IF ! q_block associated

               ! Consider print requests
               IF (should_output .AND. (print_level > low_print_level)) THEN
                  CALL para_env%sum(q_matrix)
                  IF (output_unit > 0) THEN
                     ALLOCATE (q_work(nsb, nsbsize))
                     q_work(:, :) = 0.0_dp
                     DO isb = 1, nsb
                        j = 0
                        DO i = (isb - 1)*nsbsize + 1, isb*nsbsize
                           j = j + 1
                           q_work(isb, j) = q_matrix(i, i)
                        END DO
                     END DO
                     DO isb = 1, nsb
                        WRITE (UNIT=output_unit, FMT="(T3,I6,2X,I6,2X,10F8.3)") &
                           atom_a, isb, q_work(isb, :), SUM(q_work(isb, :))
                     END DO
                     WRITE (UNIT=output_unit, FMT="(T12,A,2X,10F8.3)") &
                        "Total", (SUM(q_work(:, i)), i=1, nsbsize), SUM(q_work)
                     WRITE (UNIT=output_unit, FMT="(A)") ""
                     DEALLOCATE (q_work)
                     IF (debug) THEN
                        ! Print the DFT+U occupation matrix
                        WRITE (UNIT=output_unit, FMT="(T9,70I10)") (i, i=1, n)
                        DO i = 1, n
                           WRITE (UNIT=output_unit, FMT="(T3,I6,70F10.6)") i, q_matrix(i, :)
                        END DO
                        ! Print the eigenvalues and eigenvectors of the occupation matrix
                        ALLOCATE (q_eigval(n))
                        q_eigval(:) = 0.0_dp
                        ALLOCATE (q_eigvec(n, n))
                        q_eigvec(:, :) = 0.0_dp
                        CALL jacobi(q_matrix, q_eigval, q_eigvec)
                        WRITE (UNIT=output_unit, FMT="(/,T9,70I10)") (i, i=1, n)
                        WRITE (UNIT=output_unit, FMT="(T9,71F10.6)") (q_eigval(i), i=1, n), &
                           SUM(q_eigval(1:n))
                        DO i = 1, n
                           WRITE (UNIT=output_unit, FMT="(T3,I6,70F10.6)") i, q_eigvec(i, :)
                        END DO
                        DEALLOCATE (q_eigval)
                        DEALLOCATE (q_eigvec)
                     END IF ! debug
                  END IF
                  IF (debug) THEN
                     ! Print the full atomic occupation matrix block
                     ALLOCATE (q_work(nsgf_kind, nsgf_kind))
                     q_work(:, :) = 0.0_dp
                     IF (ASSOCIATED(q_block)) q_work(:, :) = q_block(:, :)
                     CALL para_env%sum(q_work)
                     IF (output_unit > 0) THEN
                        norb = SIZE(q_work, 1)
                        WRITE (UNIT=output_unit, FMT="(/,T9,200I10)") (i, i=1, norb)
                        DO i = 1, norb
                           WRITE (UNIT=output_unit, FMT="(T3,I6,200F10.6)") i, q_work(i, :)
                        END DO
                        ALLOCATE (q_eigval(norb))
                        q_eigval(:) = 0.0_dp
                        ALLOCATE (q_eigvec(norb, norb))
                        q_eigvec(:, :) = 0.0_dp
                        CALL jacobi(q_work, q_eigval, q_eigvec)
                        WRITE (UNIT=output_unit, FMT="(/,T9,200I10)") (i, i=1, norb)
                        WRITE (UNIT=output_unit, FMT="(T9,201F10.6)") (q_eigval(i), i=1, norb), &
                           SUM(q_eigval(1:norb))
                        DO i = 1, norb
                           WRITE (UNIT=output_unit, FMT="(T3,I6,200F10.6)") i, q_eigvec(i, :)
                        END DO
                        DEALLOCATE (q_eigval)
                        DEALLOCATE (q_eigvec)
                     END IF
                     DEALLOCATE (q_work)
                  END IF ! debug
               END IF ! should output

            END DO ! next atom "iatom" of atomic kind "ikind"

            IF (dft_control%mtlr_dft_with_perturbation .AND. &
                dft_control%mtlr_ikind == ikind) THEN
               dft_control%vhxc(ispin) = dft_control%vhxc(ispin)/REAL(n, dp)
               CALL para_env%sum(dft_control%vhxc(ispin))
               CALL para_env%sum(dft_control%trq(ispin))
            END IF

            DEALLOCATE (shell_sgf, q_matrix)

            IF (ALLOCATED(q_opp_matrix)) THEN
               DEALLOCATE (q_opp_matrix)
            END IF

         END DO ! next atomic kind "ikind"

         ! Add V(i,j)[U] to V(i,j)[DFT]
         IF (PRESENT(matrix_h) .OR. (PRESENT(matrix_w) .AND. do_kpoints)) THEN
            IF (do_kpoints) THEN
               BLOCK
                  REAL(dp) :: diagonal(nsgf)
                  CALL dbcsr_get_diag(sm_v, diagonal)
                  CALL para_env%sum(diagonal)
                  kpoints%lowdin_v(:, ispin) = diagonal(kpoints%lowdin_orbitals)
               END BLOCK
            ELSE
               sm_h => matrix_h(ispin, 1)%matrix
               CALL cp_dbcsr_sm_fm_multiply(sm_v, fm_s_half, fm_work1, nsgf)
               CALL cp_fm_transpose(fm_work1, fm_work2)
               CALL cp_dbcsr_plus_fm_fm_t(sm_h, fm_s_half, fm_work2, nsgf)
            END IF
         END IF ! An update of the Hamiltonian matrix is requested

         ! Calculate the contribution (non-Pulay part) to the derivatives
         ! w.r.t. the nuclear positions
         IF (PRESENT(matrix_w) .AND. .NOT. do_kpoints) THEN

            sm_p => matrix_p(ispin, 1)%matrix
            sm_w => matrix_w(ispin, 1)%matrix

            CALL cp_dbcsr_sm_fm_multiply(sm_v, fm_s_half, fm_work1, nsgf)
            CALL cp_fm_transpose(fm_work1, fm_work2)
            CALL cp_dbcsr_sm_fm_multiply(sm_p, fm_work2, fm_work1, nsgf)
            CALL parallel_gemm('N', 'N', nsgf, nsgf, nsgf, 1.0_dp, fm_work1, fm_sev, 0.0_dp, fm_work2)
            CALL parallel_gemm('T', 'N', nsgf, nsgf, nsgf, 1.0_dp, fm_sev, fm_work2, 0.0_dp, fm_work1)
            CALL cp_fm_schur_product(fm_work1, slambda, fm_work2)
            CALL cp_fm_transpose(fm_work2, fm_work1)
            ! Symmetrize before projecting onto the stored triangular sparse blocks.
            CALL cp_fm_scale_and_add(alpha=1.0_dp, matrix_a=fm_work2, matrix_b=fm_work1)
            CALL parallel_gemm('N', 'N', nsgf, nsgf, nsgf, -1.0_dp, fm_sev, fm_work2, 0.0_dp, fm_work1)
            CALL cp_dbcsr_plus_fm_fm_t(sm_w, fm_work1, fm_sev, nsgf)

         END IF ! W matrix update requested

      END DO ! next spin "ispin"

      IF (PRESENT(matrix_w)) THEN
         IF (do_kpoints) CALL lowdin_kp_forces(kpoints)
         CALL cp_fm_release(matrix=fm_sev)
         CALL cp_fm_release(matrix=slambda)
      END IF

      ! Collect the energy contributions from all processes

      CALL para_env%sum(energy%dft_plus_u)

      IF (energy%dft_plus_u < 0.0_dp) THEN
         CALL cp_warn(__LOCATION__, &
                      "DFT+U energy contribution is negative possibly due "// &
                      "to unphysical Lowdin charges!")
      END IF

      ! Release (local) full matrices
      NULLIFY (fm_s_half)
      CALL cp_fm_release(matrix=fm_work1)
      CALL cp_fm_release(matrix=fm_work2)

      ! Release (local) sparse matrices
      CALL dbcsr_release(sm_q)
      CALL dbcsr_release(sm_v)

      IF (if_plus_j .AND. nspin == 2) THEN
         CALL dbcsr_release(sm_q_opp)
      END IF

      CALL timestop(handle)

   END SUBROUTINE lowdin

! **************************************************************************************************
!> \brief Select the AO columns used by the existing orbital-diagonal Lowdin functional.
!> \param qs_env Quickstep environment, after S(k)**(1/2) has been prepared.
! **************************************************************************************************
   SUBROUTINE lowdin_kp_init(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      INTEGER                                            :: i, ia, ik, ikind, iset, ish, lu, m, nao, &
                                                            nset, nspin
      INTEGER, ALLOCATABLE                               :: first_atom(:)
      INTEGER, POINTER                                   :: atoms(:), first(:, :), l(:, :), &
                                                            last(:, :), nshell(:), orbitals(:)
      LOGICAL                                            :: active
      LOGICAL, ALLOCATABLE                               :: selected(:)
      REAL(KIND=dp)                                      :: j, u
      TYPE(atomic_kind_type), POINTER                    :: kinds_set(:)
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_type), POINTER                  :: basis
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(particle_type), POINTER                       :: particles(:)
      TYPE(qs_kind_type), POINTER                        :: qs_kinds(:)

      CALL get_qs_env(qs_env, kpoints=kpoints, atomic_kind_set=kinds_set, qs_kind_set=qs_kinds, &
                      particle_set=particles, dft_control=dft_control)
      IF (dft_control%mtlr_dft_with_perturbation) THEN
         CPABORT("MTLR with Lowdin DFT+U and k-points is not implemented")
      END IF
      DO ik = 1, SIZE(kpoints%kp_sym)
         IF (kpoints%kp_sym(ik)%kpoint_sym%apply_symmetry) THEN
            CPABORT("Lowdin DFT+U requires a full k-point grid or time-reversal-only reduction")
         END IF
      END DO
      CALL get_qs_kind_set(qs_kinds, nsgf=nao)
      nspin = dft_control%nspins
      ALLOCATE (first_atom(SIZE(particles)), selected(nao))
      CALL get_particle_set(particles, qs_kinds, first_sgf=first_atom)
      selected = .FALSE.
      DO ikind = 1, SIZE(qs_kinds)
         NULLIFY (orbitals)
         CALL get_qs_kind(qs_kinds(ikind), dft_plus_u_atom=active, l_of_dft_plus_u=lu, &
                          basis_set=basis, u_minus_j_target=u, hund_j=j, orbitals=orbitals)
         IF (.NOT. active) CYCLE
         IF (ASSOCIATED(orbitals)) THEN
            CPABORT("Lowdin DFT+U with k-points does not support ENFORCE_OCCUPATION")
         END IF
         IF (u == 0.0_dp .AND. j == 0.0_dp) CYCLE
         CALL get_atomic_kind(kinds_set(ikind), atom_list=atoms)
         CALL get_gto_basis_set(basis, first_sgf=first, last_sgf=last, l=l, nset=nset, nshell=nshell)
         DO ia = 1, SIZE(atoms)
            DO iset = 1, nset
               DO ish = 1, nshell(iset)
                  IF (l(ish, iset) /= lu) CYCLE
                  selected(first_atom(atoms(ia)) + first(ish, iset) - 1: &
                           first_atom(atoms(ia)) + last(ish, iset) - 1) = .TRUE.
               END DO
            END DO
         END DO
      END DO
      IF (ALLOCATED(kpoints%lowdin_orbitals)) THEN
         DEALLOCATE (kpoints%lowdin_orbitals, kpoints%lowdin_q, kpoints%lowdin_q_new, kpoints%lowdin_v)
      END IF
      m = COUNT(selected)
      ALLOCATE (kpoints%lowdin_orbitals(m))
      kpoints%lowdin_orbitals(:) = PACK([(i, i=1, nao)], selected)
      ALLOCATE (kpoints%lowdin_q(m, nspin), kpoints%lowdin_q_new(m, nspin), kpoints%lowdin_v(m, nspin))
      kpoints%lowdin_q(:, :) = 0.0_dp
      kpoints%lowdin_q_new(:, :) = 0.0_dp
      kpoints%lowdin_v(:, :) = 0.0_dp
      kpoints%lowdin_density_ready = .FALSE.
      kpoints%lowdin_population_ready = .FALSE.
      DO ik = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(ik)%kpoint_env
         CALL lowdin_kp_prepare(kp, kpoints%lowdin_orbitals, kpoints%use_real_wfn)
      END DO
   END SUBROUTINE lowdin_kp_init

! **************************************************************************************************
!> \brief Prepare selected Lowdin projectors and their reusable product workspace.
!> \param kp One k point; SCF already provides its overlap square root.
!> \param orbitals Selected global AO indices.
!> \param use_real_wfn Whether to use the real overlap square root.
!> \param smat Optional complex overlap for post-SCF reconstruction; remains unchanged.
! **************************************************************************************************
   SUBROUTINE lowdin_kp_prepare(kp, orbitals, use_real_wfn, smat)
      TYPE(kpoint_env_type), POINTER                     :: kp
      INTEGER, INTENT(IN)                                :: orbitals(:)
      LOGICAL, INTENT(IN)                                :: use_real_wfn
      TYPE(cp_cfm_type), INTENT(IN), OPTIONAL            :: smat

      INTEGER                                            :: i, m, n
      REAL(dp), ALLOCATABLE                              :: eigenvalues(:)
      TYPE(cp_cfm_type)                                  :: evec, work

      CALL cp_fm_release(kp%lowdin_b)
      CALL cp_fm_release(kp%lowdin_work)
      CALL cp_cfm_release(kp%lowdin_cb)
      CALL cp_cfm_release(kp%lowdin_cwork)
      CALL kp%lowdin_gemm%destroy()
      m = SIZE(orbitals)
      IF (m == 0) RETURN
      IF (PRESENT(smat)) THEN
         CPASSERT(.NOT. use_real_wfn)
         CALL cp_cfm_get_info(smat, nrow_global=n)
         ALLOCATE (eigenvalues(n))
         CALL cp_cfm_create(evec, smat%matrix_struct)
         CALL cp_cfm_create(work, smat%matrix_struct)
         CALL cp_cfm_release(kp%cshalf)
         CALL cp_cfm_create(kp%cshalf, smat%matrix_struct)
         CALL cp_cfm_to_cfm(smat, work)
         CALL cp_cfm_heevd(work, evec, eigenvalues)
         CPASSERT(ALL(eigenvalues >= 0.0_dp))
         CALL cp_cfm_to_cfm(evec, work)
         CALL cp_cfm_column_scale(work, CMPLX(SQRT(eigenvalues), 0.0_dp, dp))
         CALL parallel_gemm('N', 'C', n, n, n, z_one, work, evec, z_zero, kp%cshalf)
         CALL cp_cfm_release(evec)
         CALL cp_cfm_release(work)
      END IF
      CALL kp%lowdin_gemm%create(LOCAL_GEMM_PU_GPU, timing=.FALSE.)
      IF (use_real_wfn) THEN
         CALL cp_fm_create(kp%lowdin_b, kp%shalf%matrix_struct, ncol=m)
         CALL cp_fm_create(kp%lowdin_work, kp%lowdin_b%matrix_struct)
         DO i = 1, m
            CALL cp_fm_to_fm(kp%shalf, kp%lowdin_b, 1, orbitals(i), i)
         END DO
         ! The real driver also uses complex storage during its ADIIS/CDIIS handover.
         CALL cp_cfm_create(kp%lowdin_cb, kp%lowdin_b%matrix_struct)
         CALL cp_fm_to_cfm(kp%lowdin_b, mtarget=kp%lowdin_cb)
      ELSE
         CALL cp_cfm_create(kp%lowdin_cb, kp%cshalf%matrix_struct, ncol=m)
         DO i = 1, m
            CALL cp_cfm_to_cfm(kp%cshalf, kp%lowdin_cb, 1, orbitals(i), i)
         END DO
      END IF
      CALL cp_cfm_create(kp%lowdin_cwork, kp%lowdin_cb%matrix_struct)
      IF (PRESENT(smat)) CALL cp_cfm_release(kp%cshalf)
   END SUBROUTINE lowdin_kp_prepare

! **************************************************************************************************
!> \brief Complete the initial density's Lowdin coordinates, including atomic and density restarts.
!> \param qs_env Quickstep environment with an initialized real-space density.
!> \param fmwork Global AO matrix workspace from the SCF environment being initialized.
! **************************************************************************************************
   SUBROUTINE lowdin_kp_initial_density(qs_env, fmwork)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), INTENT(IN)                       :: fmwork

      INTEGER                                            :: ik, ispin, kplocal, kpmax
      TYPE(cp_cfm_type)                                  :: cpmat
      TYPE(cp_fm_type)                                   :: pmat
      TYPE(dbcsr_p_type), POINTER                        :: density(:, :)
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(kpoint_operator_context_type)                 :: ctx
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(qs_rho_type), POINTER                         :: rho

      CALL get_qs_env(qs_env, kpoints=kpoints, rho=rho)
      IF (.NOT. ALLOCATED(kpoints%lowdin_orbitals)) RETURN
      IF (.NOT. kpoints%lowdin_density_ready .AND. SIZE(kpoints%lowdin_orbitals) > 0) THEN
         ! The initial real-space guess defines a finite Fourier series. Subsequent
         ! accepted densities use their full Bloch matrices, never this sparse reconstruction.
         CALL qs_rho_get(rho, rho_ao_kp=density)
         CALL kpoint_operator_context_create(ctx, kpoints, fmwork, SIZE(density, 1))
         kp => kpoints%kp_env(1)%kpoint_env
         IF (kpoints%use_real_wfn) THEN
            CALL cp_fm_create(pmat, kp%shalf%matrix_struct)
         ELSE
            CALL cp_cfm_create(cpmat, kp%cshalf%matrix_struct)
         END IF
         kplocal = SIZE(kpoints%kp_env)
         kpmax = MAXVAL(kpoints%kp_dist(2, :) - kpoints%kp_dist(1, :) + 1)
         kpoints%lowdin_q_new(:, :) = 0.0_dp
         DO ik = 1, kpmax
            DO ispin = 1, SIZE(density, 1)
               IF (kpoints%use_real_wfn) THEN
                  CALL kpoint_operator_get(ctx, ik, ispin, density, fm_re=pmat, matrix_row=ispin)
               ELSE
                  CALL kpoint_operator_get(ctx, ik, ispin, density, cmat=cpmat, matrix_row=ispin)
               END IF
               IF (ik > kplocal) CYCLE
               kp => kpoints%kp_env(ik)%kpoint_env
               IF (kpoints%use_real_wfn) THEN
                  CALL lowdin_kp_accumulate(kp, kpoints%lowdin_q_new(:, ispin), pmat=pmat)
               ELSE
                  CALL lowdin_kp_accumulate(kp, kpoints%lowdin_q_new(:, ispin), cpmat=cpmat)
               END IF
            END DO
         END DO
         CALL kpoints%para_env%sum(kpoints%lowdin_q_new)
         CALL kpoint_operator_context_release(ctx)
         CALL cp_fm_release(pmat)
         CALL cp_cfm_release(cpmat)
         kpoints%lowdin_population_ready = .TRUE.
      END IF
      CALL lowdin_kp_population(kpoints, alpha=1.0_dp)
   END SUBROUTINE lowdin_kp_initial_density

! **************************************************************************************************
!> \brief Cache the current occupations, optionally measure their residual and accept or mix them.
!> \param kpoints Optional k-point state; absent or inactive Lowdin state is a no-op.
!> \param delta Include the unmixed occupation residual in this convergence measure.
!> \param alpha Weight of the current occupations in the accepted state; absent leaves it unchanged.
! **************************************************************************************************
   SUBROUTINE lowdin_kp_population(kpoints, delta, alpha)
      TYPE(kpoint_type), OPTIONAL, POINTER               :: kpoints
      REAL(KIND=dp), INTENT(INOUT), OPTIONAL             :: delta
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: alpha

      IF (.NOT. PRESENT(kpoints)) RETURN
      IF (.NOT. ASSOCIATED(kpoints)) RETURN
      IF (.NOT. ALLOCATED(kpoints%lowdin_q)) RETURN
      IF (.NOT. kpoints%lowdin_population_ready) THEN
         IF (SIZE(kpoints%lowdin_q) > 0) THEN
            CPASSERT(kpoints%lowdin_density_ready)
            CALL lowdin_kp_project(kpoints, kpoints%lowdin_q_new)
         END IF
         kpoints%lowdin_population_ready = .TRUE.
      END IF
      IF (PRESENT(delta)) delta = MAX(delta, MAXVAL(ABS(kpoints%lowdin_q_new - kpoints%lowdin_q)))
      IF (PRESENT(alpha)) THEN
         IF (alpha == 1.0_dp) THEN
            kpoints%lowdin_q(:, :) = kpoints%lowdin_q_new
         ELSE IF (alpha /= 0.0_dp) THEN
            kpoints%lowdin_q(:, :) = (1.0_dp - alpha)*kpoints%lowdin_q + alpha*kpoints%lowdin_q_new
         END IF
      END IF
   END SUBROUTINE lowdin_kp_population

! **************************************************************************************************
!> \brief Project and reduce a Bloch density with one accumulator and matrix workspace per worker.
!> \param kpoints Selected overlap columns and parallel distribution.
!> \param q Weighted occupations; does not change the ground-state occupation caches.
!> \param density Optional split-complex response; absent uses the current Bloch density.
! **************************************************************************************************
   SUBROUTINE lowdin_kp_project(kpoints, q, density)
      TYPE(kpoint_type), POINTER                         :: kpoints
      REAL(KIND=dp), CONTIGUOUS, INTENT(OUT)             :: q(:, :)
      TYPE(cp_fm_type), INTENT(IN), OPTIONAL             :: density(:, :, :)

      INTEGER                                            :: handle, ik, ispin, kplocal, m, nspin, &
                                                            nworkers, thread
      TYPE(cp_cfm_type), ALLOCATABLE                     :: work(:)
      TYPE(kpoint_env_type), POINTER                     :: kp

      m = SIZE(q, 1)
      IF (m == 0) RETURN
      CALL timeset("lowdin_kp_population", handle)
      nspin = SIZE(q, 2)
      kplocal = SIZE(kpoints%kp_env)
      kp => kpoints%kp_env(1)%kpoint_env
      nworkers = 1
!$    IF (kpoints%para_env_kp%num_pe == 1) nworkers = MIN(omp_get_max_threads(), kplocal)
      ALLOCATE (work(nworkers))
      q(:, :) = 0.0_dp
      IF (PRESENT(density) .OR. .NOT. kpoints%use_real_wfn) THEN
         DO thread = 1, nworkers
            CALL cp_cfm_create(work(thread), kp%pmat(1, 1)%matrix_struct)
         END DO
      END IF
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) NUM_THREADS(nworkers) IF(nworkers > 1) &
!$OMP SHARED(kpoints,density,work,nspin,kplocal,nworkers) PRIVATE(ik,ispin,kp,thread) REDUCTION(+:q)
      DO ik = 1, kplocal
         thread = 1
!$       thread = omp_get_thread_num() + 1
         kp => kpoints%kp_env(ik)%kpoint_env
         DO ispin = 1, nspin
            IF (PRESENT(density)) THEN
               CALL cp_fm_to_cfm(density(ik, 1, ispin), density(ik, 2, ispin), work(thread))
               CALL lowdin_kp_accumulate(kp, q(:, ispin), cpmat=work(thread))
            ELSE IF (kpoints%use_real_wfn) THEN
               CALL lowdin_kp_accumulate(kp, q(:, ispin), pmat=kp%pmat(1, ispin))
            ELSE
               CALL cp_fm_to_cfm(kp%pmat(1, ispin), kp%pmat(2, ispin), work(thread))
               CALL lowdin_kp_accumulate(kp, q(:, ispin), cpmat=work(thread))
            END IF
         END DO
      END DO
!$OMP END PARALLEL DO
      ! Each entry is a local tile contraction: reduce once across groups AND group ranks.
      CALL kpoints%para_env%sum(q)
      DO thread = 1, nworkers
         CALL cp_cfm_release(work(thread))
      END DO
      CALL timestop(handle)
   END SUBROUTINE lowdin_kp_project

! **************************************************************************************************
!> \brief Accumulate weighted local entries of diag(B**H P B), without a full transformed density.
!> \param kp One k point, with a private reusable product workspace.
!> \param q Accumulated local contributions; the caller zeroes and reduces the worker buffers.
!> \param pmat Real Bloch density.
!> \param cpmat Complex Bloch density (exclusive with pmat).
! **************************************************************************************************
   SUBROUTINE lowdin_kp_accumulate(kp, q, pmat, cpmat)
      TYPE(kpoint_env_type), POINTER                     :: kp
      REAL(KIND=dp), INTENT(INOUT)                       :: q(:)
      TYPE(cp_fm_type), INTENT(IN), OPTIONAL             :: pmat
      TYPE(cp_cfm_type), INTENT(IN), OPTIONAL            :: cpmat

      INTEGER                                            :: j, m, n, nc, nr
      INTEGER, POINTER                                   :: columns(:)
      LOGICAL                                            :: local

      m = SIZE(q)
      CPASSERT(PRESENT(pmat) .NEQV. PRESENT(cpmat))
      IF (PRESENT(pmat)) THEN
         CALL cp_fm_get_info(kp%lowdin_b, nrow_global=n, nrow_local=nr, ncol_local=nc, col_indices=columns)
         local = PRODUCT(kp%lowdin_b%matrix_struct%context%num_pe) == 1
         IF (local) THEN
            CALL kp%lowdin_gemm%gemm('N', 'N', n, m, n, 1.0_dp, pmat%local_data, SIZE(pmat%local_data, 1), &
                                     kp%lowdin_b%local_data, SIZE(kp%lowdin_b%local_data, 1), 0.0_dp, &
                                     kp%lowdin_work%local_data, SIZE(kp%lowdin_work%local_data, 1))
         ELSE
            CALL parallel_gemm('N', 'N', n, m, n, 1.0_dp, pmat, kp%lowdin_b, 0.0_dp, kp%lowdin_work)
         END IF
         DO j = 1, nc
            q(columns(j)) = q(columns(j)) + kp%wkp* &
                            SUM(kp%lowdin_b%local_data(:nr, j)*kp%lowdin_work%local_data(:nr, j))
         END DO
      ELSE
         CALL cp_cfm_get_info(kp%lowdin_cb, nrow_global=n, nrow_local=nr, ncol_local=nc, col_indices=columns)
         CALL lowdin_kp_gemm('N', 'N', n, m, n, cpmat, kp%lowdin_cb, z_zero, kp%lowdin_cwork, kp%lowdin_gemm)
         DO j = 1, nc
            q(columns(j)) = q(columns(j)) + kp%wkp*REAL(SUM(CONJG(kp%lowdin_cb%local_data(:nr, j))* &
                                                            kp%lowdin_cwork%local_data(:nr, j)), dp)
         END DO
      END IF
   END SUBROUTINE lowdin_kp_accumulate

! **************************************************************************************************
!> \brief Add B diag(v) B**H before any solver, DIIS residual, or history storage uses H(k).
!> \param kp One k point, with its private workspace.
!> \param v Selected orbital potentials for the same raw/mixed state as the supplied Hamiltonian.
!> \param rmat Real Hamiltonian.
!> \param cmat Complex Hamiltonian (exclusive with rmat).
! **************************************************************************************************
   SUBROUTINE lowdin_kp_apply(kp, v, rmat, cmat)
      TYPE(kpoint_env_type), POINTER                     :: kp
      REAL(KIND=dp), INTENT(IN)                          :: v(:)
      TYPE(cp_fm_type), INTENT(INOUT), OPTIONAL          :: rmat
      TYPE(cp_cfm_type), INTENT(INOUT), OPTIONAL         :: cmat

      INTEGER                                            :: j, m, n, nc, nr
      INTEGER, POINTER                                   :: columns(:)
      LOGICAL                                            :: local

      m = SIZE(v)
      IF (m == 0) RETURN
      IF (ALL(v == 0.0_dp)) RETURN
      CPASSERT(PRESENT(rmat) .NEQV. PRESENT(cmat))
      IF (PRESENT(rmat)) THEN
         CALL cp_fm_get_info(kp%lowdin_b, nrow_global=n, nrow_local=nr, ncol_local=nc, col_indices=columns)
         DO j = 1, nc
            kp%lowdin_work%local_data(:nr, j) = v(columns(j))*kp%lowdin_b%local_data(:nr, j)
         END DO
         local = PRODUCT(kp%lowdin_b%matrix_struct%context%num_pe) == 1
         IF (local) THEN
            CALL kp%lowdin_gemm%gemm('N', 'T', n, n, m, 1.0_dp, &
                                     kp%lowdin_work%local_data, SIZE(kp%lowdin_work%local_data, 1), &
                                     kp%lowdin_b%local_data, SIZE(kp%lowdin_b%local_data, 1), 1.0_dp, &
                                     rmat%local_data, SIZE(rmat%local_data, 1))
         ELSE
            CALL parallel_gemm('N', 'T', n, n, m, 1.0_dp, kp%lowdin_work, kp%lowdin_b, 1.0_dp, rmat)
         END IF
      ELSE
         CALL cp_cfm_get_info(kp%lowdin_cb, nrow_global=n, nrow_local=nr, ncol_local=nc, col_indices=columns)
         DO j = 1, nc
            kp%lowdin_cwork%local_data(:nr, j) = v(columns(j))*kp%lowdin_cb%local_data(:nr, j)
         END DO
         CALL lowdin_kp_gemm('N', 'C', n, n, m, kp%lowdin_cwork, kp%lowdin_cb, z_one, cmat, kp%lowdin_gemm)
      END IF
   END SUBROUTINE lowdin_kp_apply

! **************************************************************************************************
!> \brief Add B diag(v) B**H C to the split-complex orbital gradient, without forming an AO matrix.
!> \param kp K point owning the selected overlap columns and a private GEMM context.
!> \param v Orbital potential, without k weights or occupations (applied by the caller).
!> \param coeff_re Real part of the active orbitals.
!> \param coeff_im Imaginary part of the active orbitals.
!> \param hc_re Real part of H*C, updated in place.
!> \param hc_im Imaginary part of H*C, updated in place.
! **************************************************************************************************
   SUBROUTINE lowdin_kp_apply_mo(kp, v, coeff_re, coeff_im, hc_re, hc_im)
      TYPE(kpoint_env_type), POINTER                     :: kp
      REAL(KIND=dp), INTENT(IN)                          :: v(:)
      TYPE(cp_fm_type), INTENT(IN)                       :: coeff_re, coeff_im
      TYPE(cp_fm_type), INTENT(INOUT)                    :: hc_re, hc_im

      INTEGER                                            :: i, m, n, nc, nmo, nr
      INTEGER, POINTER                                   :: rows(:)
      TYPE(cp_cfm_type)                                  :: coeff, projected

      m = SIZE(v)
      IF (m == 0 .OR. ALL(v == 0.0_dp)) RETURN
      CALL cp_fm_get_info(coeff_re, nrow_global=n, ncol_global=nmo)
      CALL cp_cfm_create(coeff, coeff_re%matrix_struct)
      CALL cp_cfm_create(projected, coeff_re%matrix_struct, nrow=m)
      CALL cp_fm_to_cfm(coeff_re, coeff_im, coeff)
      CALL lowdin_kp_gemm('C', 'N', m, nmo, n, kp%lowdin_cb, coeff, z_zero, projected, kp%lowdin_gemm)
      CALL cp_cfm_get_info(projected, nrow_local=nr, ncol_local=nc, row_indices=rows)
      DO i = 1, nr
         projected%local_data(i, :nc) = v(rows(i))*projected%local_data(i, :nc)
      END DO
      CALL lowdin_kp_gemm('N', 'N', n, nmo, m, kp%lowdin_cb, projected, z_zero, coeff, kp%lowdin_gemm)
      CALL cp_fm_get_info(coeff_re, nrow_local=nr, ncol_local=nc)
      hc_re%local_data(:nr, :nc) = hc_re%local_data(:nr, :nc) + REAL(coeff%local_data(:nr, :nc), dp)
      hc_im%local_data(:nr, :nc) = hc_im%local_data(:nr, :nc) + AIMAG(coeff%local_data(:nr, :nc))
      CALL cp_cfm_release(projected)
      CALL cp_cfm_release(coeff)
   END SUBROUTINE lowdin_kp_apply_mo

! **************************************************************************************************
!> \brief Project a Bloch density response and return its orbital potential response, including +J.
!> \param qs_env Quickstep environment with selected Lowdin projectors and current U/J.
!> \param density Unweighted split-complex Bloch density response (local k point, component, spin).
!> \param response Globally summed orbital potential response; ground-state caches are preserved.
! **************************************************************************************************
   SUBROUTINE lowdin_kp_kernel(qs_env, density, response)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), INTENT(IN)                       :: density(:, :, :)
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT)            :: response(:, :)

      INTEGER                                            :: i, iatom, ikind, m, nspin
      INTEGER, ALLOCATABLE                               :: last(:)
      REAL(KIND=dp)                                      :: j, u
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(particle_type), POINTER                       :: particles(:)
      TYPE(qs_kind_type), POINTER                        :: qs_kinds(:)

      CALL get_qs_env(qs_env, kpoints=kpoints, particle_set=particles, qs_kind_set=qs_kinds)
      m = SIZE(kpoints%lowdin_orbitals)
      nspin = SIZE(density, 3)
      ALLOCATE (response(m, nspin))
      response(:, :) = 0.0_dp
      IF (m == 0) RETURN
      CALL lowdin_kp_project(kpoints, response, density)
      ALLOCATE (last(SIZE(particles)))
      CALL get_particle_set(particles, qs_kinds, last_sgf=last)
      iatom = 1
      DO i = 1, m
         DO WHILE (kpoints%lowdin_orbitals(i) > last(iatom))
            iatom = iatom + 1
         END DO
         CALL get_atomic_kind(particles(iatom)%atomic_kind, kind_number=ikind)
         CALL get_qs_kind(qs_kinds(ikind), u_minus_j=u, hund_j=j)
         response(i, :) = 0.5_dp*REAL(nspin, dp)*(j*response(i, nspin:1:-1) - u*response(i, :))
      END DO
   END SUBROUTINE lowdin_kp_kernel

! **************************************************************************************************
!> \brief Add -dE_U/dS(k) to W(k), before its ordinary inverse Fourier transform.
!>        If A=S**(1/2), A*dA+dA*A=dS. In the overlap eigenbasis its adjoint is
!>        G_ij=(U**H*(P*A*V+V*A*P)*U)_ij/(sqrt(s_i)+sqrt(s_j)).
!>        W <- W-U*G*U**H; the existing overlap derivatives give forces and virial.
!>        k-point weights belong to the Fourier transform and are not applied here.
!> \param kpoints Final Bloch densities, orbital potentials and cached overlap spectra.
! **************************************************************************************************
   SUBROUTINE lowdin_kp_forces(kpoints)
      TYPE(kpoint_type), POINTER                         :: kpoints

      INTEGER                                            :: handle, ik, ispin, iw, nworkers, thread
      TYPE(cp_cfm_type), ALLOCATABLE                     :: work(:, :)
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(local_gemm_ctxt_type), ALLOCATABLE            :: gemm(:)

      IF (SIZE(kpoints%lowdin_orbitals) == 0) RETURN
      IF (ALL(kpoints%lowdin_v == 0.0_dp)) RETURN
      CALL timeset("lowdin_kp_forces", handle)
      nworkers = 1
!$    IF (kpoints%para_env_kp%num_pe == 1) &
!$       nworkers = MIN(omp_get_max_threads(), SIZE(kpoints%kp_env)*SIZE(kpoints%lowdin_v, 2))
      ALLOCATE (work(3, nworkers), gemm(nworkers))
      kp => kpoints%kp_env(1)%kpoint_env
      DO thread = 1, nworkers
         DO iw = 1, 2
            CALL cp_cfm_create(work(iw, thread), kp%pmat(1, 1)%matrix_struct)
         END DO
         CALL cp_cfm_create(work(3, thread), kp%lowdin_cb%matrix_struct)
         CALL gemm(thread)%create(LOCAL_GEMM_PU_GPU, timing=.FALSE.)
      END DO
!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(NONE) SCHEDULE(DYNAMIC, 1) NUM_THREADS(nworkers) IF(nworkers > 1) &
!$OMP SHARED(kpoints,work,gemm,nworkers) PRIVATE(ik,ispin,thread)
      DO ispin = 1, SIZE(kpoints%lowdin_v, 2)
         DO ik = 1, SIZE(kpoints%kp_env)
            thread = 1
!$          thread = omp_get_thread_num() + 1
            CALL lowdin_kp_force_matrix(kpoints%kp_env(ik)%kpoint_env, kpoints%lowdin_orbitals, ispin, &
                                        kpoints%lowdin_v(:, ispin), work(:, thread), gemm(thread))
         END DO
      END DO
!$OMP END PARALLEL DO
      DO thread = 1, nworkers
         CALL gemm(thread)%destroy()
         DO iw = 1, 3
            CALL cp_cfm_release(work(iw, thread))
         END DO
      END DO
      CALL timestop(handle)
   END SUBROUTINE lowdin_kp_forces

! **************************************************************************************************
!> \brief One spin/k-point's Hermitian overlap derivative, factoring through the selected AO columns.
!> \param kp K-point state; W is updated, P and the overlap spectrum are preserved.
!> \param orbitals Selected AO columns.
!> \param ispin Spin channel to update.
!> \param v Derivatives of the common U/J functional with respect to this spin's q.
!> \param work Two square matrices and selected-column scratch, private to this worker.
!> \param gemm Private local GEMM context; concurrent spins must not share the k-point context.
! **************************************************************************************************
   SUBROUTINE lowdin_kp_force_matrix(kp, orbitals, ispin, v, work, gemm)
      TYPE(kpoint_env_type), POINTER                     :: kp
      INTEGER, INTENT(IN)                                :: orbitals(:), ispin
      REAL(KIND=dp), INTENT(IN)                          :: v(:)
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: work(3)
      TYPE(local_gemm_ctxt_type), INTENT(INOUT)          :: gemm

      COMPLEX(KIND=dp)                                   :: value
      INTEGER                                            :: i, j, m, n, nc, nr
      INTEGER, POINTER                                   :: columns(:), rows(:)
      LOGICAL                                            :: local, real_only

      IF (ALL(v == 0.0_dp)) RETURN
      CPASSERT(ALLOCATED(kp%lowdin_sqrt))
      CPASSERT(ASSOCIATED(kp%lowdin_evec%local_data))
      CPASSERT(ASSOCIATED(kp%wmat))
      CALL cp_cfm_get_info(work(1), nrow_global=n, nrow_local=nr, ncol_local=nc, &
                           row_indices=rows, col_indices=columns)
      m = SIZE(orbitals)
      local = PRODUCT(work(1)%matrix_struct%context%num_pe) == 1
      real_only = SIZE(kp%pmat, 1) == 1
      IF (real_only) THEN
         CALL cp_fm_to_cfm(kp%pmat(1, ispin), mtarget=work(1))
      ELSE
         CALL cp_fm_to_cfm(kp%pmat(1, ispin), kp%pmat(2, ispin), work(1))
      END IF
      ! U**H * P * A * V * U = (U**H * P * B) * diag(v) * U(orbitals,:).
      ! Keep these products thin: only the final backtransform needs square GEMMs.
      CALL lowdin_kp_gemm('N', 'N', n, m, n, work(1), kp%lowdin_cb, z_zero, work(2), gemm)
      CALL lowdin_kp_gemm('C', 'N', n, m, n, kp%lowdin_evec, work(2), z_zero, work(1), gemm)
      DO j = 1, nc
         IF (columns(j) <= m) work(1)%local_data(:nr, j) = v(columns(j))*work(1)%local_data(:nr, j)
      END DO
      IF (local) THEN
         DO j = 1, m
            work(3)%local_data(:n, j) = CONJG(kp%lowdin_evec%local_data(orbitals(j), :n))
         END DO
      ELSE
         CALL cp_cfm_transpose(kp%lowdin_evec, 'C', work(2))
         DO j = 1, m
            CALL cp_cfm_to_cfm(work(2), work(3), 1, orbitals(j), j)
         END DO
      END IF
      CALL lowdin_kp_gemm('N', 'C', n, n, m, work(1), work(3), z_zero, work(2), gemm)
      ! Eigenvalue sums remain regular for equal or nearly equal overlap eigenvalues.
      DO j = 1, nc
         DO i = 1, nr
            work(2)%local_data(i, j) = work(2)%local_data(i, j)/ &
                                       (kp%lowdin_sqrt(rows(i)) + kp%lowdin_sqrt(columns(j)))
         END DO
      END DO
      CALL lowdin_kp_gemm('N', 'N', n, n, n, kp%lowdin_evec, work(2), z_zero, work(1), gemm)
      CALL lowdin_kp_gemm('N', 'C', n, n, n, work(1), kp%lowdin_evec, z_zero, work(2), gemm)
      IF (.NOT. local) CALL cp_cfm_transpose(work(2), 'C', work(1))
      DO j = 1, nc
         DO i = 1, nr
            IF (local) THEN
               value = work(2)%local_data(i, j) + CONJG(work(2)%local_data(j, i))
            ELSE
               value = work(2)%local_data(i, j) + work(1)%local_data(i, j)
            END IF
            kp%wmat(1, ispin)%local_data(i, j) = kp%wmat(1, ispin)%local_data(i, j) - REAL(value, dp)
            IF (.NOT. real_only) THEN
               kp%wmat(2, ispin)%local_data(i, j) = kp%wmat(2, ispin)%local_data(i, j) - AIMAG(value)
            END IF
         END DO
      END DO
   END SUBROUTINE lowdin_kp_force_matrix

! **************************************************************************************************
!> \brief C = op(A)*op(B) + beta*C, with a private local context or collective parallel GEMM.
!> \param ta Left transpose flag.
!> \param tb Right transpose flag.
!> \param nrow Output rows.
!> \param ncol Output columns.
!> \param inner Contraction dimension; may use only the selected columns of a square workspace.
!> \param a Left operand.
!> \param b Right operand.
!> \param beta Scale of the existing output.
!> \param c Output, distinct from both inputs.
!> \param gemm Context owned by the current k-point or force worker.
! **************************************************************************************************
   SUBROUTINE lowdin_kp_gemm(ta, tb, nrow, ncol, inner, a, b, beta, c, gemm)
      CHARACTER, INTENT(IN)                              :: ta, tb
      INTEGER, INTENT(IN)                                :: nrow, ncol, inner
      TYPE(cp_cfm_type), INTENT(IN)                      :: a, b
      COMPLEX(KIND=dp), INTENT(IN)                       :: beta
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: c
      TYPE(local_gemm_ctxt_type), INTENT(INOUT)          :: gemm

      IF (PRODUCT(c%matrix_struct%context%num_pe) == 1) THEN
         CALL gemm%gemm(ta, tb, nrow, ncol, inner, z_one, &
                        a%local_data, SIZE(a%local_data, 1), b%local_data, SIZE(b%local_data, 1), &
                        beta, c%local_data, SIZE(c%local_data, 1))
      ELSE
         CALL parallel_gemm(ta, tb, nrow, ncol, inner, z_one, a, b, beta, c)
      END IF
   END SUBROUTINE lowdin_kp_gemm

! **************************************************************************************************
!> \brief         Add a DFT+U contribution to the Hamiltonian matrix\n
!>                using a method based on the Mulliken population analysis
!>                \f[q_{\mu\nu} = \frac{1}{2} (P_{\mu\nu} S_{\nu\mu} +
!>                                             S_{\mu\nu} P_{\nu\mu})\f]
!>                where \b P and \b S are the density and the
!>                overlap matrix, respectively.
!> \param[in]     qs_env Quickstep environment
!> \param orthonormal_basis ...
!> \param[in,out] matrix_h Hamiltonian matrices for each spin
!> \param should_output ...
!> \param output_unit ...
!> \param print_level ...
!> \date          03.07.2008
!> \par
!>  \f{eqnarray*}{
!>   E^{\rm DFT+U} & = & E^{\rm DFT} + E^{\rm U}\\\
!>                 & = & E^{\rm DFT} + \frac{1}{2}\sum_A(U_A - J_A)(Tr(q_A) - Tr(q^2_A))\\[1ex]
!>   V_{\mu\nu}^{\rm DFT+U} & = & V_{\mu\nu}^{\rm DFT} + V_{\mu\nu}^{\rm U}\\\
!>                          & = & \frac{\partial E^{\rm DFT}}
!>                                     {\partial P_{\mu\nu}} +
!>                                \frac{\partial E^{\rm U}}
!>                                     {\partial P_{\mu\nu}}\\\
!>                          & = & H_{\mu\nu} + \sum_A
!>                                \frac{\partial E^{\rm U}}{\partial q_A}
!>                                \frac{\partial q_A}{\partial P_{\mu\nu}}\\\
!>  \f}
!> \author        Matthias Krack (MK)
!> \version       1.0
! **************************************************************************************************
   SUBROUTINE mulliken(qs_env, orthonormal_basis, matrix_h, should_output, &
                       output_unit, print_level)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: orthonormal_basis
      TYPE(dbcsr_p_type), DIMENSION(:, :), OPTIONAL, &
         POINTER                                         :: matrix_h
      LOGICAL, INTENT(IN)                                :: should_output
      INTEGER, INTENT(IN)                                :: output_unit, print_level

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'mulliken'

      CHARACTER(LEN=10)                                  :: spin_info
      CHARACTER(LEN=6), ALLOCATABLE, DIMENSION(:)        :: symbol
      CHARACTER(LEN=default_string_length)               :: atomic_kind_name
      INTEGER :: atom_a, handle, i, i0, iatom, ic, ikind, iorb, isb, iset, isgf, ishell, ispin, j, &
         jsb, jset, jsgf, jshell, lu, m, max_scf, n, natom, natom_of_kind, nimg, nkind, norb, nsb, &
         nsbsize, nset, nsgf_kind, nspin
      INTEGER, DIMENSION(1)                              :: iloc
      INTEGER, DIMENSION(:), POINTER                     :: atom_list, nshell, orbitals
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf, l, last_sgf
      LOGICAL                                            :: debug, dft_plus_u_atom, found, &
                                                            just_energy, occupation_enforced, smear
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: is_plus_u_kind, orb_occ
      REAL(KIND=dp)                                      :: eps_scf, eps_u_ramping, fspin, occ, trq, &
                                                            trq2, u_minus_j, u_minus_j_target, &
                                                            u_ramping
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: q_eigval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: q_eigvec, q_matrix, q_work
      REAL(KIND=dp), DIMENSION(:), POINTER               :: nelec
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: h_block, p_block, q_block, s_block, &
                                                            v_block
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atomic_kind_type), POINTER                    :: kind_a
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_p, matrix_s
      TYPE(dbcsr_type), POINTER                          :: sm_h, sm_p, sm_q, sm_s, sm_v
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho

      CALL timeset(routineN, handle)

      debug = .FALSE. ! Set to .TRUE. to print debug information

      NULLIFY (atom_list)
      NULLIFY (atomic_kind_set)
      NULLIFY (qs_kind_set)
      NULLIFY (dft_control)
      NULLIFY (energy)
      NULLIFY (first_sgf)
      NULLIFY (h_block)
      NULLIFY (matrix_p)
      NULLIFY (matrix_s)
      NULLIFY (l)
      NULLIFY (last_sgf)
      NULLIFY (nelec)
      NULLIFY (nshell)
      NULLIFY (orb_basis_set)
      NULLIFY (p_block)
      NULLIFY (particle_set)
      NULLIFY (q_block)
      NULLIFY (rho)
      NULLIFY (s_block)
      NULLIFY (orbitals)
      NULLIFY (sm_h)
      NULLIFY (sm_p)
      NULLIFY (sm_q)
      NULLIFY (sm_s)
      NULLIFY (sm_v)
      NULLIFY (v_block)
      NULLIFY (para_env)

      smear = .FALSE.
      max_scf = -1
      eps_scf = 1.0E30_dp
      occupation_enforced = .FALSE.

      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      dft_control=dft_control, &
                      energy=energy, &
                      particle_set=particle_set, &
                      rho=rho, &
                      para_env=para_env)

      CPASSERT(ASSOCIATED(atomic_kind_set))
      CPASSERT(ASSOCIATED(dft_control))
      CPASSERT(ASSOCIATED(energy))
      CPASSERT(ASSOCIATED(particle_set))
      CPASSERT(ASSOCIATED(rho))

      IF (orthonormal_basis) THEN
         NULLIFY (sm_s)
      ELSE
         ! Get overlap matrix in sparse format
         CALL get_qs_env(qs_env=qs_env, &
                         matrix_s_kp=matrix_s)
         CPASSERT(ASSOCIATED(matrix_s))
      END IF
      nimg = dft_control%nimages

      ! Get density matrices in sparse format

      CALL qs_rho_get(rho, rho_ao_kp=matrix_p)

      energy%dft_plus_u = 0.0_dp

      nspin = dft_control%nspins

      IF (nspin == 2) THEN
         fspin = 1.0_dp
      ELSE
         fspin = 0.5_dp
      END IF

      ! Get the total number of atoms, contracted spherical Gaussian basis
      ! functions, and atomic kinds

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                               natom=natom)

      nkind = SIZE(atomic_kind_set)

      ALLOCATE (is_plus_u_kind(nkind))
      is_plus_u_kind(:) = .FALSE.

      IF (PRESENT(matrix_h)) THEN
         just_energy = .FALSE.
      ELSE
         just_energy = .TRUE.
      END IF

      ! Loop over all spins
      DO ispin = 1, nspin

         ! Loop over cell images
         DO ic = 1, nimg
            IF (.NOT. orthonormal_basis) THEN
               sm_s => matrix_s(1, ic)%matrix
            END IF

            IF (PRESENT(matrix_h)) THEN
               ! Hamiltonian matrix for spin ispin in sparse format
               sm_h => matrix_h(ispin, ic)%matrix
            ELSE
               NULLIFY (sm_h)
            END IF

            ! Get density matrix for spin ispin in sparse format

            sm_p => matrix_p(ispin, ic)%matrix

            IF (.NOT. ASSOCIATED(sm_q)) THEN
               ALLOCATE (sm_q)
               CALL dbcsr_get_block_diag(sm_p, sm_q)
            END IF
            CALL dbcsr_set(sm_q, 0.0_dp)

            IF (.NOT. ASSOCIATED(sm_v)) THEN
               ALLOCATE (sm_v)
               CALL dbcsr_get_block_diag(sm_p, sm_v)
            END IF
            CALL dbcsr_set(sm_v, 0.0_dp)

            DO iatom = 1, natom

               CALL dbcsr_get_block_p(matrix=sm_p, &
                                      row=iatom, &
                                      col=iatom, &
                                      block=p_block, &
                                      found=found)

               IF (.NOT. ASSOCIATED(p_block)) CYCLE

               CALL dbcsr_get_block_p(matrix=sm_q, &
                                      row=iatom, &
                                      col=iatom, &
                                      block=q_block, &
                                      found=found)
               CPASSERT(ASSOCIATED(q_block))

               IF (orthonormal_basis) THEN
                  ! S is the unit matrix
                  DO isgf = 1, SIZE(q_block, 1)
                     q_block(isgf, isgf) = p_block(isgf, isgf)
                  END DO
               ELSE
                  CALL dbcsr_get_block_p(matrix=sm_s, &
                                         row=iatom, &
                                         col=iatom, &
                                         block=s_block, &
                                         found=found)
                  CPASSERT(ASSOCIATED(s_block))
                  ! Exploit that P and S are symmetric
                  DO jsgf = 1, SIZE(p_block, 2)
                     DO isgf = 1, SIZE(p_block, 1)
                        q_block(isgf, jsgf) = p_block(isgf, jsgf)*s_block(isgf, jsgf)
                     END DO
                  END DO
               END IF ! orthonormal basis set

            END DO ! next atom "iatom"

            ! E[DFT+U] = E[DFT] + E[U]
            !          = E[DFT] + (U - J)*(Tr(q) - Tr(q*q))/2

            ! V(i,j)[DFT+U] = V(i,j)[DFT] + V(i,j)[U]
            !               = dE[DFT]/dP(i,j) + dE[U]/dP(i,j)
            !               = dE[DFT]/dP(i,j) + (dE(U)/dq)*(dq/dP(i,j))

            ! Loop over all atomic kinds

            DO ikind = 1, nkind

               ! Load the required atomic kind data

               CALL get_atomic_kind(atomic_kind_set(ikind), &
                                    atom_list=atom_list, &
                                    name=atomic_kind_name, &
                                    natom=natom_of_kind)

               CALL get_qs_kind(qs_kind_set(ikind), &
                                dft_plus_u_atom=dft_plus_u_atom, &
                                l_of_dft_plus_u=lu, &
                                nsgf=nsgf_kind, &
                                basis_set=orb_basis_set, &
                                u_minus_j=u_minus_j, &
                                u_minus_j_target=u_minus_j_target, &
                                u_ramping=u_ramping, &
                                eps_u_ramping=eps_u_ramping, &
                                nelec=nelec, &
                                orbitals=orbitals, &
                                eps_scf=eps_scf, &
                                max_scf=max_scf, &
                                smear=smear)

               ! Check, if the atoms of this atomic kind need a DFT+U correction

               IF (.NOT. ASSOCIATED(orb_basis_set)) CYCLE
               IF (.NOT. dft_plus_u_atom) CYCLE
               IF (lu < 0) CYCLE

               ! Apply U ramping if requested

               IF ((ispin == 1) .AND. (u_ramping > 0.0_dp)) THEN
                  IF (qs_env%scf_env%iter_delta <= eps_u_ramping) THEN
                     u_minus_j = MIN(u_minus_j + u_ramping, u_minus_j_target)
                     CALL set_qs_kind(qs_kind_set(ikind), u_minus_j=u_minus_j)
                  END IF
                  IF (should_output .AND. (output_unit > 0)) THEN
                     WRITE (UNIT=output_unit, FMT="(T3,A,3X,A,F0.3,A)") &
                        "Kind name: "//TRIM(ADJUSTL(atomic_kind_name)), &
                        "U(eff) = ", u_minus_j*evolt, " eV"
                  END IF
               END IF

               IF (u_minus_j == 0.0_dp) CYCLE

               is_plus_u_kind(ikind) = .TRUE.

               ! Load the required Gaussian basis set data

               CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                      first_sgf=first_sgf, &
                                      l=l, &
                                      last_sgf=last_sgf, &
                                      nset=nset, &
                                      nshell=nshell)

               ! Count the relevant shell blocks of this atomic kind

               nsb = 0
               DO iset = 1, nset
                  DO ishell = 1, nshell(iset)
                     IF (l(ishell, iset) == lu) nsb = nsb + 1
                  END DO
               END DO

               nsbsize = (2*lu + 1)
               n = nsb*nsbsize

               ALLOCATE (q_matrix(n, n))
               q_matrix(:, :) = 0.0_dp

               ! Print headline if requested

               IF (should_output .AND. (print_level > low_print_level)) THEN
                  IF (output_unit > 0) THEN
                     ALLOCATE (symbol(nsbsize))
                     DO m = -lu, lu
                        symbol(lu + m + 1) = sgf_symbol(0, lu, m)
                     END DO
                     IF (nspin > 1) THEN
                        WRITE (UNIT=spin_info, FMT="(A8,I2)") " of spin", ispin
                     ELSE
                        spin_info = ""
                     END IF
                     WRITE (UNIT=output_unit, FMT="(/,T3,A,I0,A,/,/,T5,A,10(2X,A6))") &
                        "DFT+U occupations"//TRIM(spin_info)//" for the atoms of atomic kind ", ikind, &
                        ": "//TRIM(atomic_kind_name), &
                        "Atom   Shell  ", (ADJUSTR(symbol(i)), i=1, nsbsize), " Trace"
                     DEALLOCATE (symbol)
                  END IF
               END IF

               ! Loop over all atoms of the current atomic kind

               DO iatom = 1, natom_of_kind

                  atom_a = atom_list(iatom)

                  q_matrix(:, :) = 0.0_dp

                  ! Get diagonal block

                  CALL dbcsr_get_block_p(matrix=sm_q, &
                                         row=atom_a, &
                                         col=atom_a, &
                                         block=q_block, &
                                         found=found)

                  ! Calculate energy contribution to E(U)

                  IF (ASSOCIATED(q_block)) THEN

                     i = 0
                     DO iset = 1, nset
                        DO ishell = 1, nshell(iset)
                           IF (l(ishell, iset) /= lu) CYCLE
                           DO isgf = first_sgf(ishell, iset), last_sgf(ishell, iset)
                              i = i + 1
                              j = 0
                              DO jset = 1, nset
                                 DO jshell = 1, nshell(jset)
                                    IF (l(jshell, jset) /= lu) CYCLE
                                    DO jsgf = first_sgf(jshell, jset), last_sgf(jshell, jset)
                                       j = j + 1
                                       q_matrix(i, j) = q_block(isgf, jsgf)
                                    END DO ! next contracted spherical Gaussian function "jsgf"
                                 END DO ! next shell "jshell"
                              END DO ! next shell set "jset"
                           END DO ! next contracted spherical Gaussian function "isgf"
                        END DO ! next shell "ishell"
                     END DO ! next shell set "iset"

                     ! Perform the requested manipulations of the (initial) orbital occupations

                     IF (ASSOCIATED(orbitals)) THEN
                        IF ((qs_env%scf_env%iter_delta >= eps_scf) .OR. &
                            ((qs_env%scf_env%outer_scf%iter_count == 0) .AND. &
                             (qs_env%scf_env%iter_count <= max_scf))) THEN
                           ALLOCATE (orb_occ(nsbsize))
                           ALLOCATE (q_eigval(n))
                           q_eigval(:) = 0.0_dp
                           ALLOCATE (q_eigvec(n, n))
                           q_eigvec(:, :) = 0.0_dp
                           norb = SIZE(orbitals)
                           CALL jacobi(q_matrix, q_eigval, q_eigvec)
                           q_matrix(:, :) = 0.0_dp
                           IF (nelec(ispin) >= 0.5_dp) THEN
                              trq = nelec(ispin)/SUM(q_eigval(1:n))
                              q_eigval(1:n) = trq*q_eigval(1:n)
                           END IF
                           DO isb = 1, nsb
                              trq = 0.0_dp
                              DO i = (isb - 1)*nsbsize + 1, isb*nsbsize
                                 trq = trq + q_eigval(i)
                              END DO
                              IF (smear) THEN
                                 occ = trq/REAL(norb, KIND=dp)
                              ELSE
                                 occ = 1.0_dp/fspin
                              END IF
                              orb_occ(:) = .FALSE.
                              iloc = MAXLOC(q_eigvec(:, isb*nsbsize))
                              jsb = INT((iloc(1) - 1)/nsbsize) + 1
                              i = 0
                              i0 = (jsb - 1)*nsbsize + 1
                              iorb = -1000
                              DO j = i0, jsb*nsbsize
                                 i = i + 1
                                 IF (i > norb) THEN
                                    DO m = -lu, lu
                                       IF (.NOT. orb_occ(lu + m + 1)) THEN
                                          iorb = i0 + lu + m
                                          orb_occ(lu + m + 1) = .TRUE.
                                       END IF
                                    END DO
                                 ELSE
                                    iorb = i0 + lu + orbitals(i)
                                    orb_occ(lu + orbitals(i) + 1) = .TRUE.
                                 END IF
                                 CPASSERT(iorb /= -1000)
                                 iloc = MAXLOC(q_eigvec(iorb, :))
                                 q_eigval(iloc(1)) = MIN(occ, trq)
                                 q_matrix(:, iloc(1)) = q_eigval(iloc(1))*q_eigvec(:, iloc(1)) ! backtransform left
                                 trq = trq - q_eigval(iloc(1))
                              END DO
                           END DO
                           q_matrix(:, :) = MATMUL(q_matrix, TRANSPOSE(q_eigvec)) ! backtransform right
                           DEALLOCATE (orb_occ)
                           DEALLOCATE (q_eigval)
                           DEALLOCATE (q_eigvec)
                           occupation_enforced = .TRUE.
                        END IF
                     END IF ! orbitals associated

                     trq = 0.0_dp
                     trq2 = 0.0_dp

                     DO i = 1, n
                        trq = trq + q_matrix(i, i)
                        DO j = 1, n
                           trq2 = trq2 + q_matrix(i, j)*q_matrix(j, i)
                        END DO
                     END DO

                     trq = fspin*trq
                     trq2 = fspin*fspin*trq2

                     ! Calculate energy contribution to E(U)

                     energy%dft_plus_u = energy%dft_plus_u + 0.5_dp*u_minus_j*(trq - trq2)/fspin

                     ! Calculate potential V(U) = dE(U)/dq

                     IF (.NOT. just_energy) THEN

                        CALL dbcsr_get_block_p(matrix=sm_v, &
                                               row=atom_a, &
                                               col=atom_a, &
                                               block=v_block, &
                                               found=found)
                        CPASSERT(ASSOCIATED(v_block))

                        i = 0
                        DO iset = 1, nset
                           DO ishell = 1, nshell(iset)
                              IF (l(ishell, iset) /= lu) CYCLE
                              DO isgf = first_sgf(ishell, iset), last_sgf(ishell, iset)
                                 i = i + 1
                                 j = 0
                                 DO jset = 1, nset
                                    DO jshell = 1, nshell(jset)
                                       IF (l(jshell, jset) /= lu) CYCLE
                                       DO jsgf = first_sgf(jshell, jset), last_sgf(jshell, jset)
                                          j = j + 1
                                          IF (isgf == jsgf) THEN
                                             v_block(isgf, isgf) = u_minus_j*(0.5_dp - fspin*q_matrix(i, i))
                                          ELSE
                                             v_block(isgf, jsgf) = -u_minus_j*fspin*q_matrix(j, i)
                                          END IF
                                       END DO ! next contracted spherical Gaussian function "jsgf"
                                    END DO ! next shell "jshell"
                                 END DO ! next shell set "jset"
                              END DO ! next contracted spherical Gaussian function "isgf"
                           END DO ! next shell "ishell"
                        END DO ! next shell set "iset"

                     END IF ! not just energy

                  END IF ! q_block associated

                  ! Consider print requests

                  IF (should_output .AND. (print_level > low_print_level)) THEN
                     CALL para_env%sum(q_matrix)
                     IF (output_unit > 0) THEN
                        ALLOCATE (q_work(nsb, nsbsize))
                        q_work(:, :) = 0.0_dp
                        DO isb = 1, nsb
                           j = 0
                           DO i = (isb - 1)*nsbsize + 1, isb*nsbsize
                              j = j + 1
                              q_work(isb, j) = q_matrix(i, i)
                           END DO
                        END DO
                        DO isb = 1, nsb
                           WRITE (UNIT=output_unit, FMT="(T3,I6,2X,I6,2X,10F8.3)") &
                              atom_a, isb, q_work(isb, :), SUM(q_work(isb, :))
                        END DO
                        WRITE (UNIT=output_unit, FMT="(T12,A,2X,10F8.3)") &
                           "Total", (SUM(q_work(:, i)), i=1, nsbsize), SUM(q_work)
                        WRITE (UNIT=output_unit, FMT="(A)") ""
                        DEALLOCATE (q_work)
                        IF (debug) THEN
                           ! Print the DFT+U occupation matrix
                           WRITE (UNIT=output_unit, FMT="(T9,70I10)") (i, i=1, n)
                           DO i = 1, n
                              WRITE (UNIT=output_unit, FMT="(T3,I6,70F10.6)") i, q_matrix(i, :)
                           END DO
                           ! Print the eigenvalues and eigenvectors of the occupation matrix
                           ALLOCATE (q_eigval(n))
                           q_eigval(:) = 0.0_dp
                           ALLOCATE (q_eigvec(n, n))
                           q_eigvec(:, :) = 0.0_dp
                           CALL jacobi(q_matrix, q_eigval, q_eigvec)
                           WRITE (UNIT=output_unit, FMT="(/,T9,70I10)") (i, i=1, n)
                           WRITE (UNIT=output_unit, FMT="(T9,71F10.6)") (q_eigval(i), i=1, n), &
                              SUM(q_eigval(1:n))
                           DO i = 1, n
                              WRITE (UNIT=output_unit, FMT="(T3,I6,70F10.6)") i, q_eigvec(i, :)
                           END DO
                           DEALLOCATE (q_eigval)
                           DEALLOCATE (q_eigvec)
                        END IF ! debug
                     END IF
                     IF (debug) THEN
                        ! Print the full atomic occupation matrix block
                        ALLOCATE (q_work(nsgf_kind, nsgf_kind))
                        q_work(:, :) = 0.0_dp
                        IF (ASSOCIATED(q_block)) q_work(:, :) = q_block(:, :)
                        CALL para_env%sum(q_work)
                        IF (output_unit > 0) THEN
                           norb = SIZE(q_work, 1)
                           WRITE (UNIT=output_unit, FMT="(/,T9,200I10)") (i, i=1, norb)
                           DO i = 1, norb
                              WRITE (UNIT=output_unit, FMT="(T3,I6,200F10.6)") i, q_work(i, :)
                           END DO
                           ALLOCATE (q_eigval(norb))
                           q_eigval(:) = 0.0_dp
                           ALLOCATE (q_eigvec(norb, norb))
                           q_eigvec(:, :) = 0.0_dp
                           CALL jacobi(q_work, q_eigval, q_eigvec)
                           WRITE (UNIT=output_unit, FMT="(/,T9,200I10)") (i, i=1, norb)
                           WRITE (UNIT=output_unit, FMT="(T9,201F10.6)") (q_eigval(i), i=1, norb), &
                              SUM(q_eigval(1:norb))
                           DO i = 1, norb
                              WRITE (UNIT=output_unit, FMT="(T3,I6,200F10.6)") i, q_eigvec(i, :)
                           END DO
                           DEALLOCATE (q_eigval)
                           DEALLOCATE (q_eigvec)
                        END IF
                        DEALLOCATE (q_work)
                     END IF ! debug
                  END IF ! should output

               END DO ! next atom "iatom" of atomic kind "ikind"

               IF (ALLOCATED(q_matrix)) THEN
                  DEALLOCATE (q_matrix)
               END IF

            END DO ! next atomic kind "ikind"

            ! Add V(i,j)[U] to V(i,j)[DFT]

            IF (ASSOCIATED(sm_h)) THEN

               DO ikind = 1, nkind

                  IF (.NOT. is_plus_u_kind(ikind)) CYCLE

                  kind_a => atomic_kind_set(ikind)

                  CALL get_atomic_kind(atomic_kind=kind_a, &
                                       atom_list=atom_list, &
                                       natom=natom_of_kind)

                  DO iatom = 1, natom_of_kind

                     atom_a = atom_list(iatom)

                     CALL dbcsr_get_block_p(matrix=sm_h, &
                                            row=atom_a, &
                                            col=atom_a, &
                                            block=h_block, &
                                            found=found)

                     IF (.NOT. ASSOCIATED(h_block)) CYCLE

                     CALL dbcsr_get_block_p(matrix=sm_v, &
                                            row=atom_a, &
                                            col=atom_a, &
                                            block=v_block, &
                                            found=found)
                     CPASSERT(ASSOCIATED(v_block))

                     IF (orthonormal_basis) THEN
                        DO isgf = 1, SIZE(h_block, 1)
                           h_block(isgf, isgf) = h_block(isgf, isgf) + v_block(isgf, isgf)
                        END DO
                     ELSE
                        CALL dbcsr_get_block_p(matrix=sm_s, &
                                               row=atom_a, &
                                               col=atom_a, &
                                               block=s_block, &
                                               found=found)
                        CPASSERT(ASSOCIATED(s_block))
                        DO jsgf = 1, SIZE(h_block, 2)
                           DO isgf = 1, SIZE(h_block, 1)
                              h_block(isgf, jsgf) = h_block(isgf, jsgf) + v_block(isgf, jsgf)*s_block(isgf, jsgf)
                           END DO
                        END DO
                     END IF ! orthonormal basis set

                  END DO ! next atom "iatom" of atomic kind "ikind"

               END DO ! Next atomic kind "ikind"

            END IF ! An update of the Hamiltonian matrix is requested

         END DO ! next cell image

      END DO ! next spin "ispin"

      ! Collect the energy contributions from all processes

      CALL para_env%sum(energy%dft_plus_u)

      IF (energy%dft_plus_u < 0.0_dp) THEN
         IF (.NOT. occupation_enforced) THEN
            CALL cp_warn(__LOCATION__, &
                         "DFT+U energy contribution is negative possibly due "// &
                         "to unphysical Mulliken charges!")
         END IF
      END IF

      CALL dbcsr_deallocate_matrix(sm_q)
      CALL dbcsr_deallocate_matrix(sm_v)

      CALL timestop(handle)

   END SUBROUTINE mulliken

! **************************************************************************************************
!> \brief         Add a DFT+U contribution to the Hamiltonian matrix\n
!>                using a method based on Mulliken charges
!>                \f[q_\mu = \sum_\nu \frac{1}{2}(P_{\mu\nu} S_{\nu\mu} +
!>                                                S_{\mu\nu} P_{\nu\mu})
!>                         = \sum_\nu P_{\mu\nu} S_{\nu\mu}\f]
!>                where \b P and \b S are the density and the
!>                overlap matrix, respectively.
!> \param[in]     qs_env Quickstep environment
!> \param orthonormal_basis ...
!> \param[in,out] matrix_h Hamiltonian matrices for each spin
!> \param[in,out] matrix_w Energy weighted density matrices for each spin
!> \param should_output ...
!> \param output_unit ...
!> \param print_level ...
!> \date          11.01.2008
!> \par
!>  \f{eqnarray*}{
!>   E^{\rm DFT+U} & = & E^{\rm DFT} +  E^{\rm U}\\\
!>                 & = & E^{\rm DFT} +  \frac{1}{2}(U - J)\sum_\mu (q_\mu - q_\mu^2)\\[1ex]
!>   V_{\mu\nu}^{\rm DFT+U} & = & V_{\mu\nu}^{\rm DFT} + V_{\mu\nu}^{\rm U}\\\
!>                          & = & \frac{\partial E^{\rm DFT}}
!>                                     {\partial P_{\mu\nu}} +
!>                                \frac{\partial E^{\rm U}}
!>                                     {\partial P_{\mu\nu}}\\\
!>                          & = & H_{\mu\nu} +
!>                                \frac{\partial E^{\rm U}}{\partial q_\mu}
!>                                \frac{\partial q_\mu}{\partial P_{\mu\nu}}\\\
!>                          & = & H_{\mu\nu} +
!>                                \frac{1}{2}(U - J)(1 - q_\mu - q_\nu) S_{\mu\nu}\\\
!>  \f}
!> \author        Matthias Krack (MK)
!> \version       1.0
!> \note          The use of any full matrices was avoided. Thus no ScaLAPACK
!>                calls are performed
! **************************************************************************************************
   SUBROUTINE mulliken_charges(qs_env, orthonormal_basis, matrix_h, matrix_w, &
                               should_output, output_unit, print_level)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: orthonormal_basis
      TYPE(dbcsr_p_type), DIMENSION(:, :), OPTIONAL, &
         POINTER                                         :: matrix_h, matrix_w
      LOGICAL, INTENT(IN)                                :: should_output
      INTEGER, INTENT(IN)                                :: output_unit, print_level

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'mulliken_charges'

      CHARACTER(LEN=10)                                  :: spin_info
      CHARACTER(LEN=6), ALLOCATABLE, DIMENSION(:)        :: symbol
      CHARACTER(LEN=default_string_length)               :: atomic_kind_name
      INTEGER :: atom_a, handle, i, iatom, ic, ikind, isb, iset, isgf, ishell, ispin, jatom, jsgf, &
         lu, m, natom, natom_of_kind, nimg, nkind, nsb, nset, nsgf, nspin, sgf
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: first_sgf_atom
      INTEGER, DIMENSION(:), POINTER                     :: atom_list, nshell
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf, l, last_sgf
      LOGICAL                                            :: dft_plus_u_atom, found, just_energy
      REAL(KIND=dp)                                      :: eps_u_ramping, fspin, q, u_minus_j, &
                                                            u_minus_j_target, u_ramping, v
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: dEdq, trps
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: q_ii
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: h_block, p_block, s_block, w_block
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_p, matrix_s
      TYPE(dbcsr_type), POINTER                          :: sm_h, sm_p, sm_s, sm_w
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho

      CALL timeset(routineN, handle)

      NULLIFY (atom_list)
      NULLIFY (atomic_kind_set)
      NULLIFY (qs_kind_set)
      NULLIFY (dft_control)
      NULLIFY (energy)
      NULLIFY (first_sgf)
      NULLIFY (h_block)
      NULLIFY (matrix_p)
      NULLIFY (matrix_s)
      NULLIFY (l)
      NULLIFY (last_sgf)
      NULLIFY (nshell)
      NULLIFY (orb_basis_set)
      NULLIFY (p_block)
      NULLIFY (particle_set)
      NULLIFY (rho)
      NULLIFY (s_block)
      NULLIFY (sm_h)
      NULLIFY (sm_p)
      NULLIFY (sm_s)
      NULLIFY (w_block)
      NULLIFY (para_env)

      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      dft_control=dft_control, &
                      energy=energy, &
                      particle_set=particle_set, &
                      rho=rho, &
                      para_env=para_env)

      CPASSERT(ASSOCIATED(atomic_kind_set))
      CPASSERT(ASSOCIATED(dft_control))
      CPASSERT(ASSOCIATED(energy))
      CPASSERT(ASSOCIATED(particle_set))
      CPASSERT(ASSOCIATED(rho))

      IF (orthonormal_basis) THEN
         NULLIFY (sm_s)
      ELSE
         ! Get overlap matrix in sparse format
         CALL get_qs_env(qs_env=qs_env, &
                         matrix_s_kp=matrix_s)
         CPASSERT(ASSOCIATED(matrix_s))
      END IF

      ! Get density matrices in sparse format

      CALL qs_rho_get(rho, rho_ao_kp=matrix_p)

      energy%dft_plus_u = 0.0_dp

      nspin = dft_control%nspins
      nimg = dft_control%nimages

      IF (nspin == 2) THEN
         fspin = 1.0_dp
      ELSE
         fspin = 0.5_dp
      END IF

      ! Get the total number of atoms, contracted spherical Gaussian basis
      ! functions, and atomic kinds

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, natom=natom)
      CALL get_qs_kind_set(qs_kind_set, nsgf=nsgf)

      nkind = SIZE(atomic_kind_set)

      ALLOCATE (first_sgf_atom(natom))
      first_sgf_atom(:) = 0

      CALL get_particle_set(particle_set, qs_kind_set, &
                            first_sgf=first_sgf_atom)

      ALLOCATE (trps(nsgf))
      trps(:) = 0.0_dp

      IF (PRESENT(matrix_h) .OR. PRESENT(matrix_w)) THEN
         ALLOCATE (dEdq(nsgf))
         just_energy = .FALSE.
      ELSE
         just_energy = .TRUE.
      END IF

      ! Loop over all spins

      DO ispin = 1, nspin

         IF (.NOT. just_energy) dEdq(:) = 0.0_dp

         ! Calculate Trace(P*S) assuming symmetric matrices

         trps(:) = 0.0_dp

         DO ic = 1, nimg
            IF (orthonormal_basis) THEN
               NULLIFY (sm_s)
            ELSE
               sm_s => matrix_s(1, ic)%matrix
            END IF
            sm_p => matrix_p(ispin, ic)%matrix ! Density matrix for spin ispin in sparse format

            CALL dbcsr_iterator_start(iter, sm_p)

            DO WHILE (dbcsr_iterator_blocks_left(iter))

               CALL dbcsr_iterator_next_block(iter, iatom, jatom, p_block)

               IF (orthonormal_basis) THEN

                  IF (iatom /= jatom) CYCLE

                  IF (ASSOCIATED(p_block)) THEN
                     sgf = first_sgf_atom(iatom)
                     DO isgf = 1, SIZE(p_block, 1)
                        trps(sgf) = trps(sgf) + p_block(isgf, isgf)
                        sgf = sgf + 1
                     END DO
                  END IF

               ELSE

                  CALL dbcsr_get_block_p(matrix=sm_s, &
                                         row=iatom, &
                                         col=jatom, &
                                         block=s_block, &
                                         found=found)
                  CPASSERT(ASSOCIATED(s_block))

                  sgf = first_sgf_atom(jatom)
                  DO jsgf = 1, SIZE(p_block, 2)
                     DO isgf = 1, SIZE(p_block, 1)
                        trps(sgf) = trps(sgf) + p_block(isgf, jsgf)*s_block(isgf, jsgf)
                     END DO
                     sgf = sgf + 1
                  END DO

                  IF (iatom /= jatom) THEN
                     sgf = first_sgf_atom(iatom)
                     DO isgf = 1, SIZE(p_block, 1)
                        DO jsgf = 1, SIZE(p_block, 2)
                           trps(sgf) = trps(sgf) + p_block(isgf, jsgf)*s_block(isgf, jsgf)
                        END DO
                        sgf = sgf + 1
                     END DO
                  END IF

               END IF ! orthonormal basis set

            END DO ! next atom "iatom"

            CALL dbcsr_iterator_stop(iter)

         END DO ! cell images

         CALL para_env%sum(trps)

         ! q <- Trace(PS)

         ! E[DFT+U] = E[DFT] + E[U]
         !          = E[DFT] + (U - J)*(q - q**2))/2

         ! V(i,j)[DFT+U] = V(i,j)[DFT] + V(i,j)[U]
         !               = dE[DFT]/dP(i,j) + dE[U]/dP(i,j)
         !               = dE[DFT]/dP(i,j) + (dE(U)/dq)*(dq/dP(i,j))

         ! Loop over all atomic kinds

         DO ikind = 1, nkind

            ! Load the required atomic kind data
            CALL get_atomic_kind(atomic_kind_set(ikind), &
                                 atom_list=atom_list, &
                                 name=atomic_kind_name, &
                                 natom=natom_of_kind)

            CALL get_qs_kind(qs_kind_set(ikind), &
                             dft_plus_u_atom=dft_plus_u_atom, &
                             l_of_dft_plus_u=lu, &
                             basis_set=orb_basis_set, &
                             u_minus_j=u_minus_j, &
                             u_minus_j_target=u_minus_j_target, &
                             u_ramping=u_ramping, &
                             eps_u_ramping=eps_u_ramping)

            ! Check, if this atom needs a DFT+U correction

            IF (.NOT. ASSOCIATED(orb_basis_set)) CYCLE
            IF (.NOT. dft_plus_u_atom) CYCLE
            IF (lu < 0) CYCLE

            ! Apply U ramping if requested

            IF ((ispin == 1) .AND. (u_ramping > 0.0_dp)) THEN
               IF (qs_env%scf_env%iter_delta <= eps_u_ramping) THEN
                  u_minus_j = MIN(u_minus_j + u_ramping, u_minus_j_target)
                  CALL set_qs_kind(qs_kind_set(ikind), u_minus_j=u_minus_j)
               END IF
               IF (should_output .AND. (output_unit > 0)) THEN
                  WRITE (UNIT=output_unit, FMT="(T3,A,3X,A,F0.3,A)") &
                     "Kind name: "//TRIM(ADJUSTL(atomic_kind_name)), &
                     "U(eff) = ", u_minus_j*evolt, " eV"
               END IF
            END IF

            IF (u_minus_j == 0.0_dp) CYCLE

            ! Load the required Gaussian basis set data

            CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                   first_sgf=first_sgf, &
                                   l=l, &
                                   last_sgf=last_sgf, &
                                   nset=nset, &
                                   nshell=nshell)

            ! Count the relevant shell blocks of this atomic kind

            nsb = 0
            DO iset = 1, nset
               DO ishell = 1, nshell(iset)
                  IF (l(ishell, iset) == lu) nsb = nsb + 1
               END DO
            END DO

            ALLOCATE (q_ii(nsb, 2*lu + 1))

            ! Print headline if requested

            IF (should_output .AND. (print_level > low_print_level)) THEN
               IF (output_unit > 0) THEN
                  ALLOCATE (symbol(2*lu + 1))
                  DO m = -lu, lu
                     symbol(lu + m + 1) = sgf_symbol(0, lu, m)
                  END DO
                  IF (nspin > 1) THEN
                     WRITE (UNIT=spin_info, FMT="(A8,I2)") " of spin", ispin
                  ELSE
                     spin_info = ""
                  END IF
                  WRITE (UNIT=output_unit, FMT="(/,T3,A,I0,A,/,/,T5,A,10(2X,A6))") &
                     "DFT+U occupations"//TRIM(spin_info)//" for the atoms of atomic kind ", ikind, &
                     ": "//TRIM(atomic_kind_name), &
                     "Atom   Shell  ", (ADJUSTR(symbol(i)), i=1, 2*lu + 1), " Trace"
                  DEALLOCATE (symbol)
               END IF
            END IF

            ! Loop over all atoms of the current atomic kind

            DO iatom = 1, natom_of_kind

               atom_a = atom_list(iatom)

               q_ii(:, :) = 0.0_dp

               ! Get diagonal block

               CALL dbcsr_get_block_p(matrix=sm_p, &
                                      row=atom_a, &
                                      col=atom_a, &
                                      block=p_block, &
                                      found=found)

               ! Calculate E(U) and dE(U)/dq

               IF (ASSOCIATED(p_block)) THEN

                  sgf = first_sgf_atom(atom_a)

                  isb = 0
                  DO iset = 1, nset
                     DO ishell = 1, nshell(iset)
                        IF (l(ishell, iset) == lu) THEN
                           isb = isb + 1
                           i = 0
                           DO isgf = first_sgf(ishell, iset), last_sgf(ishell, iset)
                              q = fspin*trps(sgf)
                              i = i + 1
                              q_ii(isb, i) = q
                              energy%dft_plus_u = energy%dft_plus_u + &
                                                  0.5_dp*u_minus_j*(q - q**2)/fspin
                              IF (.NOT. just_energy) THEN
                                 dEdq(sgf) = dEdq(sgf) + u_minus_j*(0.5_dp - q)
                              END IF
                              sgf = sgf + 1
                           END DO ! next contracted spherical Gaussian function "isgf"
                        ELSE
                           sgf = sgf + last_sgf(ishell, iset) - first_sgf(ishell, iset) + 1
                        END IF ! angular momentum requested for DFT+U correction
                     END DO ! next shell "ishell"
                  END DO ! next shell set "iset"

               END IF ! this process is the owner of the sparse matrix block?

               ! Consider print requests

               IF (should_output .AND. (print_level > low_print_level)) THEN
                  CALL para_env%sum(q_ii)
                  IF (output_unit > 0) THEN
                     DO isb = 1, nsb
                        WRITE (UNIT=output_unit, FMT="(T3,I6,2X,I6,2X,10F8.3)") &
                           atom_a, isb, q_ii(isb, :), SUM(q_ii(isb, :))
                     END DO
                     WRITE (UNIT=output_unit, FMT="(T12,A,2X,10F8.3)") &
                        "Total", (SUM(q_ii(:, i)), i=1, 2*lu + 1), SUM(q_ii)
                     WRITE (UNIT=output_unit, FMT="(A)") ""
                  END IF
               END IF ! should output

            END DO ! next atom "iatom" of atomic kind "ikind"

            IF (ALLOCATED(q_ii)) THEN
               DEALLOCATE (q_ii)
            END IF

         END DO ! next atomic kind "ikind"

         IF (.NOT. just_energy) THEN
            CALL para_env%sum(dEdq)
         END IF

         ! Add V(i,j)[U] to V(i,j)[DFT]

         IF (PRESENT(matrix_h)) THEN

            DO ic = 1, nimg
               IF (orthonormal_basis) THEN
                  NULLIFY (sm_s)
               ELSE
                  sm_s => matrix_s(1, ic)%matrix
               END IF
               sm_h => matrix_h(ispin, ic)%matrix

               CALL dbcsr_iterator_start(iter, sm_h)

               DO WHILE (dbcsr_iterator_blocks_left(iter))

                  CALL dbcsr_iterator_next_block(iter, iatom, jatom, h_block)

                  IF (orthonormal_basis) THEN

                     IF (iatom /= jatom) CYCLE

                     IF (ASSOCIATED(h_block)) THEN
                        sgf = first_sgf_atom(iatom)
                        DO isgf = 1, SIZE(h_block, 1)
                           h_block(isgf, isgf) = h_block(isgf, isgf) + dEdq(sgf)
                           sgf = sgf + 1
                        END DO
                     END IF

                  ELSE

                     ! Request katom just to check for consistent sparse matrix pattern

                     CALL dbcsr_get_block_p(matrix=sm_s, &
                                            row=iatom, &
                                            col=jatom, &
                                            block=s_block, &
                                            found=found)
                     CPASSERT(ASSOCIATED(s_block))

                     ! Consider the symmetric form 1/2*(P*S + S*P) for the calculation

                     sgf = first_sgf_atom(iatom)

                     DO isgf = 1, SIZE(h_block, 1)
                        IF (dEdq(sgf) /= 0.0_dp) THEN
                           v = 0.5_dp*dEdq(sgf)
                           DO jsgf = 1, SIZE(h_block, 2)
                              h_block(isgf, jsgf) = h_block(isgf, jsgf) + v*s_block(isgf, jsgf)
                           END DO
                        END IF
                        sgf = sgf + 1
                     END DO

                     sgf = first_sgf_atom(jatom)

                     DO jsgf = 1, SIZE(h_block, 2)
                        IF (dEdq(sgf) /= 0.0_dp) THEN
                           v = 0.5_dp*dEdq(sgf)
                           DO isgf = 1, SIZE(h_block, 1)
                              h_block(isgf, jsgf) = h_block(isgf, jsgf) + v*s_block(isgf, jsgf)
                           END DO
                        END IF
                        sgf = sgf + 1
                     END DO

                  END IF ! orthonormal basis set

               END DO ! Next atom "iatom"

               CALL dbcsr_iterator_stop(iter)

            END DO

         END IF ! An update of the Hamiltonian matrix is requested

         ! Calculate the contribution (non-Pulay part) to the derivatives
         ! w.r.t. the nuclear positions, which requires an update of the
         ! energy weighted density W.

         IF (PRESENT(matrix_w) .AND. (.NOT. orthonormal_basis)) THEN

            DO ic = 1, nimg
               sm_s => matrix_s(1, ic)%matrix
               sm_p => matrix_p(ispin, ic)%matrix
               sm_w => matrix_w(ispin, ic)%matrix

               CALL dbcsr_iterator_start(iter, sm_p)

               DO WHILE (dbcsr_iterator_blocks_left(iter))

                  CALL dbcsr_iterator_next_block(iter, iatom, jatom, p_block)

                  ! Skip the diagonal blocks of the W matrix

                  IF (iatom == jatom) CYCLE

                  ! Request katom just to check for consistent sparse matrix patterns

                  CALL dbcsr_get_block_p(matrix=sm_w, &
                                         row=iatom, &
                                         col=jatom, &
                                         block=w_block, &
                                         found=found)
                  CPASSERT(ASSOCIATED(w_block))

                  ! Consider the symmetric form 1/2*(P*S + S*P) for the calculation

                  sgf = first_sgf_atom(iatom)

                  DO isgf = 1, SIZE(w_block, 1)
                     IF (dEdq(sgf) /= 0.0_dp) THEN
                        v = -0.5_dp*dEdq(sgf)
                        DO jsgf = 1, SIZE(w_block, 2)
                           w_block(isgf, jsgf) = w_block(isgf, jsgf) + v*p_block(isgf, jsgf)
                        END DO
                     END IF
                     sgf = sgf + 1
                  END DO

                  sgf = first_sgf_atom(jatom)

                  DO jsgf = 1, SIZE(w_block, 2)
                     IF (dEdq(sgf) /= 0.0_dp) THEN
                        v = -0.5_dp*dEdq(sgf)
                        DO isgf = 1, SIZE(w_block, 1)
                           w_block(isgf, jsgf) = w_block(isgf, jsgf) + v*p_block(isgf, jsgf)
                        END DO
                     END IF
                     sgf = sgf + 1
                  END DO

               END DO ! next block node "jatom"

               CALL dbcsr_iterator_stop(iter)

            END DO

         END IF ! W matrix update requested

      END DO ! next spin "ispin"

      ! Collect the energy contributions from all processes

      CALL para_env%sum(energy%dft_plus_u)

      IF (energy%dft_plus_u < 0.0_dp) THEN
         CALL cp_warn(__LOCATION__, &
                      "DFT+U energy contribution is negative possibly due "// &
                      "to unphysical Mulliken charges!")
      END IF

      ! Release local work storage

      IF (ALLOCATED(first_sgf_atom)) THEN
         DEALLOCATE (first_sgf_atom)
      END IF

      IF (ALLOCATED(trps)) THEN
         DEALLOCATE (trps)
      END IF

      IF (ALLOCATED(dEdq)) THEN
         DEALLOCATE (dEdq)
      END IF

      CALL timestop(handle)

   END SUBROUTINE mulliken_charges

! **************************************************************************************************
!> \brief Calculate the tensorial DFT+U+J energy contribution, Hamiltonian correction,
!>        and optional force contribution.
!> \param[in]     qs_env         Quickstep environment
!> \param[in,out] matrix_h       Hamiltonian matrices
!> \param[in,out] matrix_w       Energy-weighted density matrices for each spin channel
!> \param[in]     should_output  Whether detailed output should be printed
!> \param[in]     output_unit     ...
!> \param[in]     print_level     ...
!> \date          14.03.2026
!> \author        Ziwei Chai (ZC)
! **************************************************************************************************
   SUBROUTINE tensorial(qs_env, matrix_h, matrix_w, should_output, output_unit, print_level)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), OPTIONAL, &
         POINTER                                         :: matrix_h, matrix_w
      LOGICAL, INTENT(IN)                                :: should_output
      INTEGER, INTENT(IN)                                :: output_unit, print_level

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'tensorial'

      CHARACTER(LEN=10)                                  :: spin_info
      CHARACTER(LEN=6), ALLOCATABLE                      :: symbol(:)
      CHARACTER(LEN=default_string_length)               :: atomic_kind_name
      INTEGER :: atom_a, atom_aa, handle, i, i0, iatom, ikind, iorb, isb, iset, isgf, ishell, &
         ispin, j, jatom, jj, jkind, jsb, jset, jsgf, jshell, katom, lu, m, matom_of_kind, &
         max_scf, n, natom, natom_of_kind, nkind, norb, nsb, nsbsize, nset, nsgf, nsgf_kind, &
         nsgf_kind2, nspin, timevalues(8)
      INTEGER, ALLOCATABLE                               :: first_sgf_atom(:)
      INTEGER, DIMENSION(1)                              :: iloc
      INTEGER, DIMENSION(:), POINTER                     :: atom_list, atom_list2, nshell, orbitals
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf, l, last_sgf
      LOGICAL                                            :: debug, dft_plus_u_atom, found, &
                                                            just_energy, smear
      LOGICAL, ALLOCATABLE                               :: orb_occ(:)
      REAL(KIND=dp) :: eps_scf, eps_u_ramping, fspin, hund_j, occ, perturbation_strength, &
         some_real, trq, trq2, trqxq_hund_j, u_minus_j, u_minus_j_target, u_ramping
      REAL(KIND=dp), ALLOCATABLE :: a(:), q_a_matrix(:, :), q_b_matrix(:, :), q_eigval(:), &
         q_eigvec(:, :), q_work(:, :), s_inv_matrix(:, :), s_matrix(:, :), tmp_matrix1(:, :), &
         tmp_matrix2(:, :), tmp_matrix3(:, :), v_matrix(:, :), vhxc_matrix(:, :)

      TYPE :: block3d
         REAL(KIND=dp), ALLOCATABLE                       :: matrix(:, :, :)
      END TYPE block3d
      TYPE :: block4d
         REAL(KIND=dp), ALLOCATABLE                       :: matrix(:, :, :, :)
      END TYPE block4d
      TYPE(block3d), ALLOCATABLE                          :: a_matrix(:, :), d_matrix(:, :), &
                                                             e_matrix(:, :), i_matrix(:, :), &
                                                             j_matrix(:, :)
      TYPE(block4d), ALLOCATABLE                          :: dd_matrix(:, :), ee_matrix(:, :), &
                                                             ff_matrix(:, :)
      REAL(KIND=dp), DIMENSION(:), POINTER                :: ao_coef
      REAL(KIND=dp), DIMENSION(:, :), POINTER             :: a_block, d_block, e_block, f_block, &
                                                             i_block, j_block, q_a_block, q_b_block, &
                                                             s_block, tmp_block, v_block, vhxc_block

      TYPE(atomic_kind_type), DIMENSION(:), POINTER       :: atomic_kind_set
      TYPE(mp_para_env_type), POINTER                     :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER           :: matrix_p, matrix_s, matrix_vhxc
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER        :: matrixkp_s
      TYPE(dbcsr_type), POINTER                           :: sm_a, sm_b, sm_c, sm_d, sm_e, sm_f, &
                                                             sm_h, sm_i, sm_j, sm_k, sm_q_a, sm_q_b, &
                                                             sm_s, sm_v, sm_w, sm_mid, sm_vhxc, sm_p
      TYPE(dbcsr_iterator_type)                           :: iter
      TYPE(dft_control_type), POINTER                     :: dft_control
      TYPE(gto_basis_set_type), POINTER                   :: orb_basis_set
      TYPE(particle_type), DIMENSION(:), POINTER          :: particle_set
      TYPE(qs_energy_type), POINTER                       :: energy
      TYPE(qs_force_type), DIMENSION(:), POINTER          :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER           :: qs_kind_set
      TYPE(qs_rho_type), POINTER                          :: rho
      TYPE(qs_scf_env_type), POINTER                      :: scf_env

      CALL timeset(routineN, handle)

      debug = .FALSE.
      smear = .FALSE.
      max_scf = -1
      eps_scf = 1.0E30_dp

      NULLIFY (atom_list)
      NULLIFY (atom_list2)
      NULLIFY (atomic_kind_set)
      NULLIFY (qs_kind_set)
      NULLIFY (dft_control)
      NULLIFY (energy)
      NULLIFY (force)
      NULLIFY (first_sgf)
      NULLIFY (matrix_p)
      NULLIFY (matrix_s)
      NULLIFY (matrix_vhxc)
      NULLIFY (l)
      NULLIFY (last_sgf)
      NULLIFY (nshell)
      NULLIFY (orb_basis_set)
      NULLIFY (orbitals)
      NULLIFY (ao_coef)
      NULLIFY (particle_set)
      NULLIFY (q_a_block)
      NULLIFY (q_b_block)
      NULLIFY (vhxc_block)
      NULLIFY (rho)
      NULLIFY (scf_env)
      NULLIFY (sm_h)
      NULLIFY (sm_p)
      NULLIFY (sm_q_a)
      NULLIFY (sm_q_b)
      NULLIFY (sm_s)
      NULLIFY (sm_v)
      NULLIFY (sm_a)
      NULLIFY (sm_b)
      NULLIFY (sm_c)
      NULLIFY (sm_d)
      NULLIFY (sm_e)
      NULLIFY (sm_f)
      NULLIFY (sm_i)
      NULLIFY (sm_j)
      NULLIFY (sm_k)
      NULLIFY (sm_mid)
      NULLIFY (sm_vhxc)
      NULLIFY (v_block)
      NULLIFY (a_block)
      NULLIFY (d_block)
      NULLIFY (e_block)
      NULLIFY (f_block)
      NULLIFY (i_block)
      NULLIFY (j_block)
      NULLIFY (para_env)
      NULLIFY (s_block)
      NULLIFY (tmp_block)

      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      dft_control=dft_control, &
                      energy=energy, &
                      matrix_s=matrix_s, &
                      matrix_vhxc=matrix_vhxc, &
                      matrix_s_kp=matrixkp_s, &
                      particle_set=particle_set, &
                      rho=rho, &
                      scf_env=scf_env, &
                      para_env=para_env)

      CPASSERT(ASSOCIATED(atomic_kind_set))
      CPASSERT(ASSOCIATED(dft_control))
      CPASSERT(ASSOCIATED(energy))
      CPASSERT(ASSOCIATED(matrix_s))
      CPASSERT(ASSOCIATED(particle_set))
      CPASSERT(ASSOCIATED(rho))

      sm_s => matrix_s(1)%matrix
      CALL qs_rho_get(rho, rho_ao=matrix_p)

      energy%dft_plus_u = 0.0_dp

      nspin = dft_control%nspins
      IF (dft_control%mtlr_dft_with_perturbation) THEN
         IF (.NOT. ASSOCIATED(matrix_vhxc)) THEN
            CPABORT("MTLR requires the projected Hxc matrix for every spin channel.")
         END IF
         IF (SIZE(matrix_vhxc) /= nspin) THEN
            CPABORT("The number of projected Hxc matrices does not match the spin channels.")
         END IF
      END IF

      IF (nspin == 2) THEN
         fspin = 1.0_dp
      ELSE
         fspin = 0.5_dp
      END IF

      ! Get the total number of atoms, contracted spherical Gaussian basis
      ! functions, and atomic kinds

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, natom=natom)
      CALL get_qs_kind_set(qs_kind_set, nsgf=nsgf)

      nkind = SIZE(atomic_kind_set)

      ALLOCATE (first_sgf_atom(natom))
      first_sgf_atom(:) = 0

      CALL get_particle_set(particle_set, qs_kind_set, first_sgf=first_sgf_atom)

      IF (PRESENT(matrix_h) .OR. PRESENT(matrix_w)) THEN
         just_energy = .FALSE.
      ELSE
         just_energy = .TRUE.
      END IF

      CALL dbcsr_init_p(sm_mid)
      CALL dbcsr_create(sm_mid, template=sm_s, matrix_type="N")

      ! Create local block diagonal matrices

      ALLOCATE (sm_q_a)
      CALL dbcsr_get_block_diag(sm_s, sm_q_a)

      ALLOCATE (sm_q_b)
      CALL dbcsr_get_block_diag(sm_s, sm_q_b)

      CALL dbcsr_set(sm_q_a, 0.0_dp)
      CALL dbcsr_set(sm_q_b, 0.0_dp)

      ! Loop over all spins for calculating sm_q_a and sm_q_b

      DO ispin = 1, nspin

         sm_p => matrix_p(ispin)%matrix

         IF (nspin == 2) THEN
            IF (ispin == 1) THEN
               CALL dbcsr_multiply("N", "N", 1.0_dp, sm_s, sm_p, 0.0_dp, sm_mid)
               CALL dbcsr_multiply("N", "N", 1.0_dp, sm_mid, sm_s, 0.0_dp, sm_q_a)
            ELSE
               CALL dbcsr_multiply("N", "N", 1.0_dp, sm_s, sm_p, 0.0_dp, sm_mid)
               CALL dbcsr_multiply("N", "N", 1.0_dp, sm_mid, sm_s, 0.0_dp, sm_q_b)
            END IF
         ELSE
            CALL dbcsr_multiply("N", "N", 1.0_dp, sm_s, sm_p, 0.0_dp, sm_mid)
            CALL dbcsr_multiply("N", "N", 1.0_dp, sm_mid, sm_s, 0.0_dp, sm_q_a)
            CALL dbcsr_copy(sm_q_b, sm_q_a)
         END IF

      END DO

      ALLOCATE (sm_v)
      CALL dbcsr_get_block_diag(sm_s, sm_v)

      IF (PRESENT(matrix_w)) THEN

         CALL dbcsr_init_p(sm_b)
         CALL dbcsr_create(sm_b, template=matrix_s(1)%matrix, matrix_type="N")

         CALL dbcsr_init_p(sm_d)
         CALL dbcsr_create(sm_d, template=matrix_s(1)%matrix, matrix_type="N")

         CALL dbcsr_init_p(sm_c)
         CALL dbcsr_create(sm_c, template=matrix_s(1)%matrix, matrix_type="N")

         CALL dbcsr_init_p(sm_e)
         CALL dbcsr_create(sm_e, template=matrix_s(1)%matrix, matrix_type="N")

         CALL dbcsr_init_p(sm_f)
         CALL dbcsr_create(sm_f, template=matrix_s(1)%matrix, matrix_type="N")

         CALL dbcsr_init_p(sm_i)
         CALL dbcsr_create(sm_i, template=matrix_s(1)%matrix, matrix_type="N")

         CALL dbcsr_init_p(sm_j)
         CALL dbcsr_create(sm_j, template=matrix_s(1)%matrix, matrix_type="N")

         CALL dbcsr_init_p(sm_k)
         CALL dbcsr_create(sm_k, template=matrix_s(1)%matrix, matrix_type="N")

      END IF

      ALLOCATE (sm_a)
      CALL dbcsr_get_block_diag(sm_s, sm_a)

      ! Loop over all spins

      DO ispin = 1, nspin

         IF (PRESENT(matrix_h)) THEN
            sm_h => matrix_h(ispin, 1)%matrix
         ELSE
            NULLIFY (sm_h)
         END IF

         IF (PRESENT(matrix_w)) THEN
            sm_w => matrix_w(ispin, 1)%matrix
         ELSE
            NULLIFY (sm_w)
         END IF

         CALL dbcsr_set(sm_v, 0.0_dp)
         CALL dbcsr_set(sm_a, 0.0_dp)

         IF (dft_control%mtlr_dft_with_perturbation) THEN
            IF (ispin == 1) perturbation_strength = 0.9_dp*dft_control%perturbation_strength
            IF (ispin == 2) perturbation_strength = 1.1_dp*dft_control%perturbation_strength
            IF (.NOT. ASSOCIATED(matrix_vhxc(ispin)%matrix)) THEN
               CPABORT("MTLR projected Hxc matrix is not initialized.")
            END IF
            sm_vhxc => matrix_vhxc(ispin)%matrix
         END IF

         ! Loop over all atomic kinds

         DO ikind = 1, nkind

            CALL get_atomic_kind(atomic_kind_set(ikind), &
                                 atom_list=atom_list, &
                                 name=atomic_kind_name, &
                                 natom=natom_of_kind)

            CALL get_qs_kind(qs_kind_set(ikind), &
                             dft_plus_u_atom=dft_plus_u_atom, &
                             l_of_dft_plus_u=lu, &
                             nsgf=nsgf_kind, &
                             basis_set=orb_basis_set, &
                             u_minus_j=u_minus_j, &
                             hund_j=hund_j, &
                             ao_coef=ao_coef, &
                             u_minus_j_target=u_minus_j_target, &
                             u_ramping=u_ramping, &
                             eps_u_ramping=eps_u_ramping, &
                             orbitals=orbitals, &
                             eps_scf=eps_scf, &
                             max_scf=max_scf, &
                             smear=smear)

            IF (.NOT. ASSOCIATED(orb_basis_set)) CYCLE
            IF (.NOT. dft_plus_u_atom) CYCLE
            IF (lu < 0) CYCLE
            IF (u_minus_j == 0.0_dp .AND. hund_j == 0.0_dp) CYCLE
            IF (.NOT. ASSOCIATED(ao_coef)) THEN
               CPABORT("Tensorial DFT+U requires associated ao_coef")
            END IF

            ALLOCATE (a(SIZE(ao_coef)))
            a(:) = ao_coef(:)

            ! Apply U ramping if requested

            IF ((ispin == 1) .AND. (u_ramping > 0.0_dp)) THEN
               IF (qs_env%scf_env%iter_delta <= eps_u_ramping) THEN
                  u_minus_j = MIN(u_minus_j + u_ramping, u_minus_j_target)
                  CALL set_qs_kind(qs_kind_set(ikind), u_minus_j=u_minus_j)
               END IF

               IF (should_output .AND. (output_unit > 0)) THEN
                  WRITE (UNIT=output_unit, FMT="(T3,A,3X,A,F0.3,A)") &
                     "Kind name: "//TRIM(ADJUSTL(atomic_kind_name)), &
                     "U(eff) = ", u_minus_j*evolt, " eV"
               END IF
            END IF

            CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                   first_sgf=first_sgf, &
                                   l=l, &
                                   last_sgf=last_sgf, &
                                   nset=nset, &
                                   nshell=nshell)

            nsb = 1
            nsbsize = (2*lu + 1)
            n = nsb*nsbsize

            ALLOCATE (q_a_matrix(n, n))
            ALLOCATE (q_b_matrix(n, n))
            ALLOCATE (s_matrix(n, n))
            ALLOCATE (s_inv_matrix(n, n))
            ALLOCATE (v_matrix(n, n))

            q_a_matrix(:, :) = 0.0_dp
            q_b_matrix(:, :) = 0.0_dp
            s_matrix(:, :) = 0.0_dp
            s_inv_matrix(:, :) = 0.0_dp
            v_matrix(:, :) = 0.0_dp

            ! Print headline if requested

            IF (should_output .AND. (print_level > low_print_level)) THEN
               IF (output_unit > 0) THEN
                  ALLOCATE (symbol(nsbsize))
                  DO m = -lu, lu
                     symbol(lu + m + 1) = sgf_symbol(0, lu, m)
                  END DO

                  IF (nspin > 1) THEN
                     WRITE (UNIT=spin_info, FMT="(A8,I2)") " of spin", ispin
                  ELSE
                     spin_info = ""
                  END IF

                  WRITE (UNIT=output_unit, FMT="(/,T3,A,I0,A,/,/,T5,A,10(2X,A6))") &
                     "DFT+U occupations"//TRIM(spin_info)//" for the atoms of atomic kind ", ikind, &
                     ": "//TRIM(atomic_kind_name), &
                     "Atom   Shell  ", (ADJUSTR(symbol(i)), i=1, nsbsize), " Trace"

                  DEALLOCATE (symbol)
               END IF
            END IF

            IF (dft_control%mtlr_dft_with_perturbation .AND. &
                dft_control%mtlr_ikind == ikind) THEN
               dft_control%vhxc(ispin) = 0.0_dp
               dft_control%trq(ispin) = 0.0_dp
               ALLOCATE (vhxc_matrix(n, n))
               vhxc_matrix(:, :) = 0.0_dp
            END IF

            ! Loop over all atoms of the current atomic kind

            DO iatom = 1, natom_of_kind

               atom_a = atom_list(iatom)

               q_a_matrix(:, :) = 0.0_dp
               q_b_matrix(:, :) = 0.0_dp
               s_matrix(:, :) = 0.0_dp
               s_inv_matrix(:, :) = 0.0_dp
               v_matrix(:, :) = 0.0_dp

               CALL dbcsr_get_block_p(matrix=sm_q_a, &
                                      row=atom_a, &
                                      col=atom_a, &
                                      block=q_a_block, &
                                      found=found)

               CALL dbcsr_get_block_p(matrix=sm_q_b, &
                                      row=atom_a, &
                                      col=atom_a, &
                                      block=q_b_block, &
                                      found=found)

               CALL dbcsr_get_block_p(matrix=sm_s, &
                                      row=atom_a, &
                                      col=atom_a, &
                                      block=s_block, &
                                      found=found)

               IF (dft_control%mtlr_dft_with_perturbation .AND. &
                   dft_control%mtlr_ikind == ikind .AND. &
                   qs_kind_set(ikind)%dft_plus_u%lr_atom == atom_a) THEN
                  CALL dbcsr_get_block_p(matrix=sm_vhxc, &
                                         row=atom_a, &
                                         col=atom_a, &
                                         block=vhxc_block, &
                                         found=found)
               END IF

               IF (ASSOCIATED(q_a_block) .AND. ASSOCIATED(q_b_block) .AND. ASSOCIATED(s_block)) THEN

                  isb = 0
                  DO iset = 1, nset
                     DO ishell = 1, nshell(iset)
                        IF (l(ishell, iset) /= lu) CYCLE
                        isb = isb + 1
                        i = 0
                        DO isgf = first_sgf(ishell, iset), last_sgf(ishell, iset)
                           i = i + 1
                           jsb = 0
                           DO jset = 1, nset
                              DO jshell = 1, nshell(jset)
                                 IF (l(jshell, jset) /= lu) CYCLE
                                 jsb = jsb + 1
                                 j = 0
                                 DO jsgf = first_sgf(jshell, jset), last_sgf(jshell, jset)
                                    j = j + 1
                                    q_a_matrix(i, j) = q_a_matrix(i, j) + q_a_block(isgf, jsgf)*a(isb)*a(jsb)
                                    q_b_matrix(i, j) = q_b_matrix(i, j) + q_b_block(isgf, jsgf)*a(isb)*a(jsb)
                                    s_matrix(i, j) = s_matrix(i, j) + s_block(isgf, jsgf)*a(isb)*a(jsb)
                                    IF (ASSOCIATED(vhxc_block) .AND. &
                                        dft_control%mtlr_dft_with_perturbation .AND. &
                                        dft_control%mtlr_ikind == ikind .AND. &
                                        qs_kind_set(ikind)%dft_plus_u%lr_atom == atom_a) THEN
                                       vhxc_matrix(i, j) = vhxc_matrix(i, j) + vhxc_block(isgf, jsgf)*a(isb)*a(jsb)
                                    END IF
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO

                  CALL invert_matrix(s_matrix, s_inv_matrix, some_real)

                  q_a_matrix(:, :) = MATMUL(q_a_matrix, s_inv_matrix)
                  q_b_matrix(:, :) = MATMUL(q_b_matrix, s_inv_matrix)

                  IF (ASSOCIATED(vhxc_block) .AND. &
                      dft_control%mtlr_dft_with_perturbation .AND. &
                      dft_control%mtlr_ikind == ikind .AND. &
                      qs_kind_set(ikind)%dft_plus_u%lr_atom == atom_a) THEN
                     vhxc_matrix(:, :) = MATMUL(vhxc_matrix, s_inv_matrix)
                     DO i = 1, n
                        dft_control%vhxc(ispin) = dft_control%vhxc(ispin) + vhxc_matrix(i, i)
                        IF (ispin == 1) dft_control%trq(1) = dft_control%trq(1) + q_a_matrix(i, i)
                        IF (ispin == 2) dft_control%trq(2) = dft_control%trq(2) + q_b_matrix(i, i)
                     END DO
                  END IF

                  ! Perform the requested manipulations of the (initial) orbital occupations

                  IF (ASSOCIATED(orbitals)) THEN
                     IF ((qs_env%scf_env%iter_delta >= eps_scf) .OR. &
                         ((qs_env%scf_env%outer_scf%iter_count == 0) .AND. &
                          (qs_env%scf_env%iter_count <= max_scf))) THEN

                        ALLOCATE (orb_occ(nsbsize))
                        ALLOCATE (q_eigval(n))
                        ALLOCATE (q_eigvec(n, n))

                        q_eigval(:) = 0.0_dp
                        q_eigvec(:, :) = 0.0_dp
                        norb = SIZE(orbitals)

                        IF (ispin == 1) THEN
                           CALL jacobi(q_a_matrix, q_eigval, q_eigvec)
                           q_a_matrix(:, :) = 0.0_dp
                        ELSE
                           CALL jacobi(q_b_matrix, q_eigval, q_eigvec)
                           q_b_matrix(:, :) = 0.0_dp
                        END IF

                        DO isb = 1, nsb
                           trq = 0.0_dp
                           DO i = (isb - 1)*nsbsize + 1, isb*nsbsize
                              trq = trq + q_eigval(i)
                           END DO

                           IF (smear) THEN
                              occ = trq/REAL(norb, KIND=dp)
                           ELSE
                              occ = 1.0_dp/fspin
                           END IF

                           orb_occ(:) = .FALSE.
                           iloc = MAXLOC(q_eigvec(:, isb*nsbsize))
                           jsb = INT((iloc(1) - 1)/nsbsize) + 1
                           i = 0
                           i0 = (jsb - 1)*nsbsize + 1
                           iorb = -1000

                           DO j = i0, jsb*nsbsize
                              i = i + 1
                              IF (i > norb) THEN
                                 DO m = -lu, lu
                                    IF (.NOT. orb_occ(lu + m + 1)) THEN
                                       iorb = i0 + lu + m
                                       orb_occ(lu + m + 1) = .TRUE.
                                    END IF
                                 END DO
                              ELSE
                                 iorb = i0 + lu + orbitals(i)
                                 orb_occ(lu + orbitals(i) + 1) = .TRUE.
                              END IF

                              CPASSERT(iorb /= -1000)

                              iloc = MAXLOC(q_eigvec(iorb, :))
                              q_eigval(iloc(1)) = MIN(occ, trq)

                              IF (ispin == 1) THEN
                                 q_a_matrix(:, iloc(1)) = q_eigval(iloc(1))*q_eigvec(:, iloc(1))
                              ELSE
                                 q_b_matrix(:, iloc(1)) = q_eigval(iloc(1))*q_eigvec(:, iloc(1))
                              END IF

                              trq = trq - q_eigval(iloc(1))
                           END DO
                        END DO

                        IF (ispin == 1) THEN
                           q_a_matrix(:, :) = MATMUL(q_a_matrix, TRANSPOSE(q_eigvec))
                        ELSE
                           q_b_matrix(:, :) = MATMUL(q_b_matrix, TRANSPOSE(q_eigvec))
                        END IF

                        DEALLOCATE (orb_occ)
                        DEALLOCATE (q_eigval)
                        DEALLOCATE (q_eigvec)
                     END IF
                  END IF

                  trq = 0.0_dp
                  trq2 = 0.0_dp
                  trqxq_hund_j = 0.0_dp

                  IF (ispin == 1) THEN
                     DO i = 1, n
                        trq = trq + q_a_matrix(i, i)
                        DO j = 1, n
                           trq2 = trq2 + q_a_matrix(i, j)*q_a_matrix(j, i)
                           trqxq_hund_j = trqxq_hund_j + q_a_matrix(i, j)*q_b_matrix(j, i)

                           IF (i == j) THEN
                              v_matrix(i, i) = u_minus_j*(0.5_dp - fspin*q_a_matrix(i, i))
                              v_matrix(i, i) = v_matrix(i, i) + hund_j*fspin*q_b_matrix(i, i)
                              IF (dft_control%mtlr_dft_with_perturbation .AND. &
                                  dft_control%mtlr_ikind == ikind .AND. &
                                  qs_kind_set(ikind)%dft_plus_u%lr_atom == atom_a) THEN
                                 v_matrix(i, i) = v_matrix(i, i) + perturbation_strength
                              END IF
                           ELSE
                              v_matrix(i, j) = -u_minus_j*fspin*q_a_matrix(j, i)
                              v_matrix(i, j) = v_matrix(i, j) + hund_j*fspin*q_b_matrix(j, i)
                           END IF
                        END DO
                     END DO
                  ELSE
                     DO i = 1, n
                        trq = trq + q_b_matrix(i, i)
                        DO j = 1, n
                           trq2 = trq2 + q_b_matrix(i, j)*q_b_matrix(j, i)
                           trqxq_hund_j = trqxq_hund_j + q_b_matrix(i, j)*q_a_matrix(j, i)

                           IF (i == j) THEN
                              v_matrix(i, i) = u_minus_j*(0.5_dp - fspin*q_b_matrix(i, i))
                              v_matrix(i, i) = v_matrix(i, i) + hund_j*fspin*q_a_matrix(i, i)
                              IF (dft_control%mtlr_dft_with_perturbation .AND. &
                                  dft_control%mtlr_ikind == ikind .AND. &
                                  qs_kind_set(ikind)%dft_plus_u%lr_atom == atom_a) THEN
                                 v_matrix(i, i) = v_matrix(i, i) + perturbation_strength
                              END IF
                           ELSE
                              v_matrix(i, j) = -u_minus_j*fspin*q_b_matrix(j, i)
                              v_matrix(i, j) = v_matrix(i, j) + hund_j*fspin*q_a_matrix(j, i)
                           END IF
                        END DO
                     END DO
                  END IF

                  v_matrix(:, :) = MATMUL(s_inv_matrix, v_matrix)

                  trq = fspin*trq
                  trq2 = fspin*fspin*trq2
                  trqxq_hund_j = fspin*fspin*trqxq_hund_j

                  energy%dft_plus_u = energy%dft_plus_u + 0.5_dp*u_minus_j*(trq - trq2)/fspin
                  energy%dft_plus_u = energy%dft_plus_u + 0.5_dp*hund_j*trqxq_hund_j/fspin

                  ! Calculate potential V(U) = dE(U)/dq

                  IF (.NOT. just_energy) THEN

                     CALL dbcsr_get_block_p(matrix=sm_v, &
                                            row=atom_a, &
                                            col=atom_a, &
                                            block=v_block, &
                                            found=found)
                     CPASSERT(ASSOCIATED(v_block))

                     CALL dbcsr_get_block_p(matrix=sm_a, &
                                            row=atom_a, &
                                            col=atom_a, &
                                            block=a_block, &
                                            found=found)
                     CPASSERT(ASSOCIATED(a_block))

                     isb = 0
                     DO iset = 1, nset
                        DO ishell = 1, nshell(iset)
                           IF (l(ishell, iset) /= lu) CYCLE
                           isb = isb + 1
                           i = 0
                           DO isgf = first_sgf(ishell, iset), last_sgf(ishell, iset)
                              i = i + 1
                              jsb = 0
                              DO jset = 1, nset
                                 DO jshell = 1, nshell(jset)
                                    IF (l(jshell, jset) /= lu) CYCLE
                                    jsb = jsb + 1
                                    j = 0
                                    DO jsgf = first_sgf(jshell, jset), last_sgf(jshell, jset)
                                       j = j + 1
                                       v_block(isgf, jsgf) = v_matrix(i, j)*a(isb)*a(jsb)
                                       a_block(isgf, jsgf) = s_inv_matrix(i, j)*a(isb)*a(jsb)
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO

                  END IF

               END IF

               ! Consider print requests

               IF (should_output .AND. (print_level > low_print_level)) THEN
                  IF (ispin == 1) THEN
                     CALL para_env%sum(q_a_matrix)
                  ELSE
                     CALL para_env%sum(q_b_matrix)
                  END IF

                  IF (output_unit > 0) THEN
                     ALLOCATE (q_work(nsb, nsbsize))
                     q_work(:, :) = 0.0_dp

                     DO isb = 1, nsb
                        j = 0
                        DO i = (isb - 1)*nsbsize + 1, isb*nsbsize
                           j = j + 1
                           IF (ispin == 1) THEN
                              q_work(isb, j) = q_a_matrix(i, i)
                           ELSE
                              q_work(isb, j) = q_b_matrix(i, i)
                           END IF
                        END DO
                     END DO

                     DO isb = 1, nsb
                        WRITE (UNIT=output_unit, FMT="(T3,I6,2X,I6,2X,10F8.3)") &
                           atom_a, isb, q_work(isb, :), SUM(q_work(isb, :))
                     END DO

                     WRITE (UNIT=output_unit, FMT="(T12,A,2X,10F8.3)") &
                        "Total", (SUM(q_work(:, i)), i=1, nsbsize), SUM(q_work)
                     WRITE (UNIT=output_unit, FMT="(A)") ""

                     DEALLOCATE (q_work)

                     IF (debug) THEN
                        WRITE (UNIT=output_unit, FMT="(T9,70I10)") (i, i=1, n)
                        DO i = 1, n
                           IF (ispin == 1) THEN
                              WRITE (UNIT=output_unit, FMT="(T3,I6,70F10.6)") i, q_a_matrix(i, :)
                           ELSE
                              WRITE (UNIT=output_unit, FMT="(T3,I6,70F10.6)") i, q_b_matrix(i, :)
                           END IF
                        END DO

                        ALLOCATE (q_eigval(n))
                        ALLOCATE (q_eigvec(n, n))
                        q_eigval(:) = 0.0_dp
                        q_eigvec(:, :) = 0.0_dp

                        IF (ispin == 1) THEN
                           CALL jacobi(q_a_matrix, q_eigval, q_eigvec)
                        ELSE
                           CALL jacobi(q_b_matrix, q_eigval, q_eigvec)
                        END IF

                        WRITE (UNIT=output_unit, FMT="(/,T9,70I10)") (i, i=1, n)
                        WRITE (UNIT=output_unit, FMT="(T9,71F10.6)") (q_eigval(i), i=1, n), SUM(q_eigval(1:n))

                        DO i = 1, n
                           WRITE (UNIT=output_unit, FMT="(T3,I6,70F10.6)") i, q_eigvec(i, :)
                        END DO

                        DEALLOCATE (q_eigval)
                        DEALLOCATE (q_eigvec)
                     END IF
                  END IF

                  IF (debug) THEN
                     ALLOCATE (q_work(nsgf_kind, nsgf_kind))
                     q_work(:, :) = 0.0_dp

                     IF (ispin == 1) THEN
                        IF (ASSOCIATED(q_a_block)) q_work(:, :) = q_a_block(:, :)
                     ELSE
                        IF (ASSOCIATED(q_b_block)) q_work(:, :) = q_b_block(:, :)
                     END IF

                     CALL para_env%sum(q_work)

                     IF (output_unit > 0) THEN
                        norb = SIZE(q_work, 1)

                        WRITE (UNIT=output_unit, FMT="(/,T9,200I10)") (i, i=1, norb)
                        DO i = 1, norb
                           WRITE (UNIT=output_unit, FMT="(T3,I6,200F10.6)") i, q_work(i, :)
                        END DO

                        ALLOCATE (q_eigval(norb))
                        ALLOCATE (q_eigvec(norb, norb))
                        q_eigval(:) = 0.0_dp
                        q_eigvec(:, :) = 0.0_dp

                        CALL jacobi(q_work, q_eigval, q_eigvec)

                        WRITE (UNIT=output_unit, FMT="(/,T9,200I10)") (i, i=1, norb)
                        WRITE (UNIT=output_unit, FMT="(T9,201F10.6)") (q_eigval(i), i=1, norb), SUM(q_eigval(1:norb))

                        DO i = 1, norb
                           WRITE (UNIT=output_unit, FMT="(T3,I6,200F10.6)") i, q_eigvec(i, :)
                        END DO

                        DEALLOCATE (q_eigval)
                        DEALLOCATE (q_eigvec)
                     END IF

                     DEALLOCATE (q_work)
                  END IF
               END IF

            END DO

            IF (ALLOCATED(q_a_matrix)) DEALLOCATE (q_a_matrix)
            IF (ALLOCATED(q_b_matrix)) DEALLOCATE (q_b_matrix)
            IF (ALLOCATED(s_matrix)) DEALLOCATE (s_matrix)
            IF (ALLOCATED(s_inv_matrix)) DEALLOCATE (s_inv_matrix)
            IF (ALLOCATED(v_matrix)) DEALLOCATE (v_matrix)

            DEALLOCATE (a)

            IF (dft_control%mtlr_dft_with_perturbation .AND. &
                dft_control%mtlr_ikind == ikind) THEN
               dft_control%vhxc(ispin) = dft_control%vhxc(ispin)/REAL(n, dp)
               CALL para_env%sum(dft_control%vhxc(ispin))
               CALL para_env%sum(dft_control%trq(ispin))
               DEALLOCATE (vhxc_matrix)
            END IF

         END DO

         ! Add V(i,j)[U] to V(i,j)[DFT]

         IF (ASSOCIATED(sm_h)) THEN
            CALL dbcsr_multiply("N", "N", 1.0_dp, sm_s, sm_v, 0.0_dp, sm_mid)
            CALL dbcsr_multiply("N", "N", 1.0_dp, sm_mid, sm_s, 1.0_dp, sm_h)
         END IF

      END DO

      CALL dbcsr_deallocate_matrix(sm_mid)

      IF (PRESENT(matrix_w)) THEN

         CALL get_qs_env(qs_env=qs_env, force=force)
         DO ikind = 1, nkind
            force(ikind)%tensorial_u(:, :) = 0.0_dp
         END DO
         IF (nspin == 2) THEN
            CALL dbcsr_multiply("N", "N", 1.0_dp, sm_s, matrix_p(1)%matrix, 0.0_dp, sm_b)
            CALL dbcsr_multiply("N", "N", 1.0_dp, sm_s, matrix_p(2)%matrix, 0.0_dp, sm_c)
         ELSE
            CALL dbcsr_multiply("N", "N", 0.5_dp, sm_s, matrix_p(1)%matrix, 0.0_dp, sm_b)
            CALL dbcsr_multiply("N", "N", 0.5_dp, sm_s, matrix_p(1)%matrix, 0.0_dp, sm_c)
         END IF
         CALL dbcsr_multiply("N", "N", 1.0_dp, sm_b, sm_s, 0.0_dp, sm_i)
         CALL dbcsr_multiply("N", "N", 1.0_dp, sm_c, sm_s, 0.0_dp, sm_j)

         DO jj = 1, 3
            CALL dbcsr_desymmetrize(matrixkp_s(jj + 1, 1)%matrix, sm_k)

            CALL dbcsr_set(sm_f, 0.0_dp)

            CALL dbcsr_iterator_start(iter, sm_k)
            DO WHILE (dbcsr_iterator_blocks_left(iter))
               NULLIFY (tmp_block)
               CALL dbcsr_iterator_next_block(iter, jatom, katom, tmp_block)
               IF (jatom == katom) CYCLE
               CALL dbcsr_put_block(matrix=sm_f, row=jatom, col=katom, block=tmp_block)
            END DO
            CALL dbcsr_iterator_stop(iter)
            CALL dbcsr_finalize(sm_f)

            CALL dbcsr_multiply("N", "N", 1.0_dp, sm_b, sm_f, 0.0_dp, sm_d)
            CALL dbcsr_multiply("N", "N", 1.0_dp, sm_c, sm_f, 0.0_dp, sm_e)

            ALLOCATE (a_matrix(nkind, nkind))
            ALLOCATE (d_matrix(nkind, nkind))
            ALLOCATE (e_matrix(nkind, nkind))
            ALLOCATE (i_matrix(nkind, nkind))
            ALLOCATE (j_matrix(nkind, nkind))
            ALLOCATE (dd_matrix(nkind, nkind))
            ALLOCATE (ee_matrix(nkind, nkind))
            ALLOCATE (ff_matrix(nkind, nkind))

            DO ikind = 1, nkind

               CALL get_atomic_kind(atomic_kind_set(ikind), &
                                    atom_list=atom_list, &
                                    natom=natom_of_kind)

               CALL get_qs_kind(qs_kind_set(ikind), &
                                u_minus_j=u_minus_j, &
                                l_of_dft_plus_u=lu, &
                                basis_set=orb_basis_set, &
                                hund_j=hund_j, &
                                nsgf=nsgf_kind, &
                                dft_plus_u_atom=dft_plus_u_atom)

               IF (.NOT. dft_plus_u_atom) CYCLE
               IF (.NOT. ASSOCIATED(orb_basis_set)) CYCLE
               IF (lu < 0) CYCLE
               IF (u_minus_j == 0.0_dp .AND. hund_j == 0.0_dp) CYCLE

               CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                      first_sgf=first_sgf, &
                                      l=l, &
                                      last_sgf=last_sgf, &
                                      nset=nset, &
                                      nshell=nshell)

               nsb = 0
               DO iset = 1, nset
                  DO ishell = 1, nshell(iset)
                     IF (l(ishell, iset) == lu) nsb = nsb + 1
                  END DO
               END DO

               nsbsize = (2*lu + 1)
               n = nsb*nsbsize

               ALLOCATE (a_matrix(ikind, ikind)%matrix(n, n, natom_of_kind))
               ALLOCATE (d_matrix(ikind, ikind)%matrix(n, n, natom_of_kind))
               ALLOCATE (e_matrix(ikind, ikind)%matrix(n, n, natom_of_kind))
               ALLOCATE (i_matrix(ikind, ikind)%matrix(n, n, natom_of_kind))
               ALLOCATE (j_matrix(ikind, ikind)%matrix(n, n, natom_of_kind))

               a_matrix(ikind, ikind)%matrix(:, :, :) = 0.0_dp
               d_matrix(ikind, ikind)%matrix(:, :, :) = 0.0_dp
               e_matrix(ikind, ikind)%matrix(:, :, :) = 0.0_dp
               i_matrix(ikind, ikind)%matrix(:, :, :) = 0.0_dp
               j_matrix(ikind, ikind)%matrix(:, :, :) = 0.0_dp

               DO iatom = 1, natom_of_kind

                  atom_aa = atom_list(iatom)

                  NULLIFY (a_block)
                  CALL dbcsr_get_block_p(matrix=sm_a, &
                                         row=atom_aa, &
                                         col=atom_aa, &
                                         block=a_block, &
                                         found=found)

                  NULLIFY (d_block)
                  CALL dbcsr_get_block_p(matrix=sm_d, &
                                         row=atom_aa, &
                                         col=atom_aa, &
                                         block=d_block, &
                                         found=found)

                  NULLIFY (e_block)
                  CALL dbcsr_get_block_p(matrix=sm_e, &
                                         row=atom_aa, &
                                         col=atom_aa, &
                                         block=e_block, &
                                         found=found)

                  NULLIFY (i_block)
                  CALL dbcsr_get_block_p(matrix=sm_i, &
                                         row=atom_aa, &
                                         col=atom_aa, &
                                         block=i_block, &
                                         found=found)

                  NULLIFY (j_block)
                  CALL dbcsr_get_block_p(matrix=sm_j, &
                                         row=atom_aa, &
                                         col=atom_aa, &
                                         block=j_block, &
                                         found=found)

                  IF (ASSOCIATED(a_block) .AND. ASSOCIATED(d_block) .AND. ASSOCIATED(e_block) &
                      .AND. ASSOCIATED(i_block) .AND. ASSOCIATED(j_block)) THEN

                     i = 0
                     DO iset = 1, nset
                        DO ishell = 1, nshell(iset)
                           IF (l(ishell, iset) /= lu) CYCLE
                           DO isgf = first_sgf(ishell, iset), last_sgf(ishell, iset)
                              i = i + 1
                              j = 0
                              DO jset = 1, nset
                                 DO jshell = 1, nshell(jset)
                                    IF (l(jshell, jset) /= lu) CYCLE
                                    DO jsgf = first_sgf(jshell, jset), last_sgf(jshell, jset)
                                       j = j + 1
                                       a_matrix(ikind, ikind)%matrix(i, j, iatom) = a_block(isgf, jsgf)
                                       d_matrix(ikind, ikind)%matrix(i, j, iatom) = -d_block(isgf, jsgf)
                                       e_matrix(ikind, ikind)%matrix(i, j, iatom) = -e_block(isgf, jsgf)
                                       i_matrix(ikind, ikind)%matrix(i, j, iatom) = i_block(isgf, jsgf)
                                       j_matrix(ikind, ikind)%matrix(i, j, iatom) = j_block(isgf, jsgf)
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO

                  END IF

                  DO jkind = 1, nkind

                     CALL get_atomic_kind(atomic_kind_set(jkind), &
                                          atom_list=atom_list2, &
                                          natom=matom_of_kind)

                     CALL get_qs_kind(qs_kind_set(jkind), nsgf=nsgf_kind2)

                     IF (.NOT. ALLOCATED(dd_matrix(ikind, jkind)%matrix)) THEN
                        ALLOCATE (dd_matrix(ikind, jkind)%matrix(n, nsgf_kind2, natom_of_kind, matom_of_kind))
                        dd_matrix(ikind, jkind)%matrix(:, :, :, :) = 0.0_dp
                     END IF

                     IF (.NOT. ALLOCATED(ee_matrix(ikind, jkind)%matrix)) THEN
                        ALLOCATE (ee_matrix(ikind, jkind)%matrix(n, nsgf_kind2, natom_of_kind, matom_of_kind))
                        ee_matrix(ikind, jkind)%matrix(:, :, :, :) = 0.0_dp
                     END IF

                     IF (.NOT. ALLOCATED(ff_matrix(ikind, jkind)%matrix)) THEN
                        ALLOCATE (ff_matrix(ikind, jkind)%matrix(nsgf_kind2, n, natom_of_kind, matom_of_kind))
                        ff_matrix(ikind, jkind)%matrix(:, :, :, :) = 0.0_dp
                     END IF

                     DO jatom = 1, matom_of_kind

                        atom_a = atom_list2(jatom)

                        NULLIFY (d_block)
                        CALL dbcsr_get_block_p(matrix=sm_b, &
                                               row=atom_aa, &
                                               col=atom_a, &
                                               block=d_block, &
                                               found=found)

                        NULLIFY (e_block)
                        CALL dbcsr_get_block_p(matrix=sm_c, &
                                               row=atom_aa, &
                                               col=atom_a, &
                                               block=e_block, &
                                               found=found)

                        IF (ASSOCIATED(d_block) .AND. ASSOCIATED(e_block)) THEN
                           i = 0
                           DO iset = 1, nset
                              DO ishell = 1, nshell(iset)
                                 IF (l(ishell, iset) /= lu) CYCLE
                                 DO isgf = first_sgf(ishell, iset), last_sgf(ishell, iset)
                                    i = i + 1
                                    dd_matrix(ikind, jkind)%matrix(i, :, iatom, jatom) = d_block(isgf, :)
                                    ee_matrix(ikind, jkind)%matrix(i, :, iatom, jatom) = e_block(isgf, :)
                                 END DO
                              END DO
                           END DO
                        END IF

                        NULLIFY (f_block)
                        CALL dbcsr_get_block_p(matrix=sm_f, &
                                               row=atom_a, &
                                               col=atom_aa, &
                                               block=f_block, &
                                               found=found)

                        IF (ASSOCIATED(f_block)) THEN
                           i = 0
                           DO iset = 1, nset
                              DO ishell = 1, nshell(iset)
                                 IF (l(ishell, iset) /= lu) CYCLE
                                 DO isgf = first_sgf(ishell, iset), last_sgf(ishell, iset)
                                    i = i + 1
                                    ff_matrix(ikind, jkind)%matrix(:, i, iatom, jatom) = f_block(:, isgf)
                                 END DO
                              END DO
                           END DO
                        END IF

                     END DO
                  END DO
               END DO
            END DO

            DO ikind = 1, nkind

               CALL get_atomic_kind(atomic_kind_set(ikind), &
                                    atom_list=atom_list, &
                                    natom=natom_of_kind)

               CALL get_qs_kind(qs_kind_set(ikind), &
                                u_minus_j=u_minus_j, &
                                l_of_dft_plus_u=lu, &
                                basis_set=orb_basis_set, &
                                hund_j=hund_j, &
                                nsgf=nsgf_kind, &
                                dft_plus_u_atom=dft_plus_u_atom)

               IF (.NOT. dft_plus_u_atom) CYCLE
               IF (.NOT. ASSOCIATED(orb_basis_set)) CYCLE
               IF (lu < 0) CYCLE
               IF (u_minus_j == 0.0_dp .AND. hund_j == 0.0_dp) CYCLE

               CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                      first_sgf=first_sgf, &
                                      l=l, &
                                      last_sgf=last_sgf, &
                                      nset=nset, &
                                      nshell=nshell)

               CALL para_env%sum(a_matrix(ikind, ikind)%matrix(:, :, :))
               CALL para_env%sum(d_matrix(ikind, ikind)%matrix(:, :, :))
               CALL para_env%sum(e_matrix(ikind, ikind)%matrix(:, :, :))
               CALL para_env%sum(i_matrix(ikind, ikind)%matrix(:, :, :))
               CALL para_env%sum(j_matrix(ikind, ikind)%matrix(:, :, :))

               DO jkind = 1, nkind
                  CALL para_env%sum(dd_matrix(ikind, jkind)%matrix(:, :, :, :))
                  CALL para_env%sum(ee_matrix(ikind, jkind)%matrix(:, :, :, :))
                  CALL para_env%sum(ff_matrix(ikind, jkind)%matrix(:, :, :, :))
               END DO

            END DO

            ! Force calculation

            DO ikind = 1, nkind

               CALL get_atomic_kind(atomic_kind_set(ikind), &
                                    atom_list=atom_list, &
                                    natom=natom_of_kind)

               CALL get_qs_kind(qs_kind_set(ikind), &
                                u_minus_j=u_minus_j, &
                                l_of_dft_plus_u=lu, &
                                basis_set=orb_basis_set, &
                                hund_j=hund_j, &
                                nsgf=nsgf_kind, &
                                dft_plus_u_atom=dft_plus_u_atom)

               IF (.NOT. dft_plus_u_atom) CYCLE
               IF (.NOT. ASSOCIATED(orb_basis_set)) CYCLE
               IF (lu < 0) CYCLE
               IF (u_minus_j == 0.0_dp .AND. hund_j == 0.0_dp) CYCLE

               CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                      first_sgf=first_sgf, &
                                      l=l, &
                                      last_sgf=last_sgf, &
                                      nset=nset, &
                                      nshell=nshell)

               nsb = 0
               DO iset = 1, nset
                  DO ishell = 1, nshell(iset)
                     IF (l(ishell, iset) == lu) nsb = nsb + 1
                  END DO
               END DO

               nsbsize = (2*lu + 1)
               n = nsb*nsbsize

               ALLOCATE (tmp_matrix1(n, n))
               ALLOCATE (tmp_matrix2(n, n))
               ALLOCATE (tmp_matrix3(n, n))

               DO iatom = 1, natom_of_kind

                  atom_aa = atom_list(iatom)

                  tmp_matrix1(:, :) = MATMUL(d_matrix(ikind, ikind)%matrix(:, :, iatom), &
                                             a_matrix(ikind, ikind)%matrix(:, :, iatom))
                  tmp_matrix2(:, :) = MATMUL(e_matrix(ikind, ikind)%matrix(:, :, iatom), &
                                             a_matrix(ikind, ikind)%matrix(:, :, iatom))

                  i_matrix(ikind, ikind)%matrix(:, :, iatom) = MATMUL(i_matrix(ikind, ikind)%matrix(:, :, iatom), &
                                                                      a_matrix(ikind, ikind)%matrix(:, :, iatom))
                  j_matrix(ikind, ikind)%matrix(:, :, iatom) = MATMUL(j_matrix(ikind, ikind)%matrix(:, :, iatom), &
                                                                      a_matrix(ikind, ikind)%matrix(:, :, iatom))

                  tmp_matrix3(:, :) = 0.5_dp*u_minus_j*tmp_matrix1 - &
                                      u_minus_j*MATMUL(tmp_matrix1, i_matrix(ikind, ikind)%matrix(:, :, iatom)) &
                                      + hund_j*MATMUL(tmp_matrix2, i_matrix(ikind, ikind)%matrix(:, :, iatom))
                  tmp_matrix3(:, :) = tmp_matrix3 + 0.5_dp*u_minus_j*tmp_matrix2 - &
                                      u_minus_j*MATMUL(tmp_matrix2, j_matrix(ikind, ikind)%matrix(:, :, iatom)) &
                                      + hund_j*MATMUL(tmp_matrix1, j_matrix(ikind, ikind)%matrix(:, :, iatom))

                  DO i = 1, n
                     force(ikind)%tensorial_u(jj, iatom) = force(ikind)%tensorial_u(jj, iatom) + tmp_matrix3(i, i)
                  END DO

                  DO jkind = 1, nkind

                     CALL get_atomic_kind(atomic_kind_set(jkind), &
                                          atom_list=atom_list2, &
                                          natom=matom_of_kind)

                     CALL get_qs_kind(qs_kind_set(jkind), nsgf=nsgf_kind2)

                     DO jatom = 1, matom_of_kind

                        atom_a = atom_list2(jatom)

                        tmp_matrix1(:, :) = MATMUL(dd_matrix(ikind, jkind)%matrix(:, :, iatom, jatom), &
                                                   ff_matrix(ikind, jkind)%matrix(:, :, iatom, jatom))
                        tmp_matrix1(:, :) = MATMUL(tmp_matrix1, a_matrix(ikind, ikind)%matrix(:, :, iatom))

                        tmp_matrix2(:, :) = MATMUL(ee_matrix(ikind, jkind)%matrix(:, :, iatom, jatom), &
                                                   ff_matrix(ikind, jkind)%matrix(:, :, iatom, jatom))
                        tmp_matrix2(:, :) = MATMUL(tmp_matrix2, a_matrix(ikind, ikind)%matrix(:, :, iatom))

                        tmp_matrix3(:, :) = 0.5_dp*u_minus_j*tmp_matrix1 - &
                                            u_minus_j*MATMUL(tmp_matrix1, i_matrix(ikind, ikind)%matrix(:, :, iatom)) &
                                            + hund_j*MATMUL(tmp_matrix2, i_matrix(ikind, ikind)%matrix(:, :, iatom))
                        tmp_matrix3(:, :) = tmp_matrix3 + 0.5_dp*u_minus_j*tmp_matrix2 - &
                                            u_minus_j*MATMUL(tmp_matrix2, j_matrix(ikind, ikind)%matrix(:, :, iatom)) &
                                            + hund_j*MATMUL(tmp_matrix1, j_matrix(ikind, ikind)%matrix(:, :, iatom))

                        IF (atom_a /= atom_aa) THEN
                           DO i = 1, n
                              force(jkind)%tensorial_u(jj, jatom) = force(jkind)%tensorial_u(jj, jatom) &
                                                                    + tmp_matrix3(i, i)
                           END DO
                        END IF

                     END DO
                  END DO
               END DO

               DEALLOCATE (tmp_matrix1)
               DEALLOCATE (tmp_matrix2)
               DEALLOCATE (tmp_matrix3)

            END DO

            DO ikind = 1, nkind

               IF (ALLOCATED(a_matrix(ikind, ikind)%matrix)) DEALLOCATE (a_matrix(ikind, ikind)%matrix)
               IF (ALLOCATED(d_matrix(ikind, ikind)%matrix)) DEALLOCATE (d_matrix(ikind, ikind)%matrix)
               IF (ALLOCATED(e_matrix(ikind, ikind)%matrix)) DEALLOCATE (e_matrix(ikind, ikind)%matrix)
               IF (ALLOCATED(i_matrix(ikind, ikind)%matrix)) DEALLOCATE (i_matrix(ikind, ikind)%matrix)
               IF (ALLOCATED(j_matrix(ikind, ikind)%matrix)) DEALLOCATE (j_matrix(ikind, ikind)%matrix)

               DO jkind = 1, nkind
                  IF (ALLOCATED(dd_matrix(ikind, jkind)%matrix)) DEALLOCATE (dd_matrix(ikind, jkind)%matrix)
                  IF (ALLOCATED(ee_matrix(ikind, jkind)%matrix)) DEALLOCATE (ee_matrix(ikind, jkind)%matrix)
                  IF (ALLOCATED(ff_matrix(ikind, jkind)%matrix)) DEALLOCATE (ff_matrix(ikind, jkind)%matrix)
               END DO

            END DO

            DEALLOCATE (a_matrix)
            DEALLOCATE (d_matrix)
            DEALLOCATE (e_matrix)
            DEALLOCATE (i_matrix)
            DEALLOCATE (j_matrix)
            DEALLOCATE (dd_matrix)
            DEALLOCATE (ee_matrix)
            DEALLOCATE (ff_matrix)

         END DO

         DO ikind = 1, nkind
            force(ikind)%tensorial_u(:, :) = force(ikind)%tensorial_u(:, :)*2.0_dp
         END DO

         CALL DATE_AND_TIME(values=timevalues)

      END IF

      ! Collect the energy contributions from all processes

      CALL para_env%sum(energy%dft_plus_u)

      IF (energy%dft_plus_u < 0.0_dp) THEN
         CALL cp_warn(__LOCATION__, &
                      "DFT+U energy contibution is negative possibly due "// &
                      "to unphysical Lowdin charges. Check your input, "// &
                      "if this warning persists or try a different method!")
      END IF

      ! Release local sparse matrices

      CALL dbcsr_deallocate_matrix(sm_q_a)
      CALL dbcsr_deallocate_matrix(sm_q_b)
      CALL dbcsr_deallocate_matrix(sm_v)
      CALL dbcsr_deallocate_matrix(sm_a)

      IF (PRESENT(matrix_w)) THEN
         CALL dbcsr_deallocate_matrix(sm_b)
         CALL dbcsr_deallocate_matrix(sm_c)
         CALL dbcsr_deallocate_matrix(sm_d)
         CALL dbcsr_deallocate_matrix(sm_e)
         CALL dbcsr_deallocate_matrix(sm_f)
         CALL dbcsr_deallocate_matrix(sm_k)
         CALL dbcsr_deallocate_matrix(sm_i)
         CALL dbcsr_deallocate_matrix(sm_j)
      END IF

      CALL timestop(handle)

   END SUBROUTINE tensorial

END MODULE dft_plus_u
