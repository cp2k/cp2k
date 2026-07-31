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
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE bibliography,                    ONLY: Dudarev1997,&
                                              Dudarev1998,&
                                              cite_reference
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_copy, dbcsr_create, dbcsr_deallocate_matrix, dbcsr_desymmetrize, dbcsr_finalize, &
        dbcsr_get_block_p, dbcsr_init_p, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_multiply, &
        dbcsr_p_type, dbcsr_put_block, dbcsr_release, dbcsr_set, dbcsr_type
   USE cp_dbcsr_contrib,                ONLY: dbcsr_get_block_diag
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
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr,&
                                              low_print_level
   USE input_constants,                 ONLY: plus_u_lowdin,&
                                              plus_u_mulliken,&
                                              plus_u_mulliken_charges,&
                                              plus_u_tensorial
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE kpoint_methods,                  ONLY: lowdin_kp_trans
   USE kpoint_types,                    ONLY: kpoint_type
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
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_scf_types,                    ONLY: qs_scf_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dft_plus_u'

   PUBLIC :: plus_u

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
         jsb, jset, jsgf, jshell, lu, m, max_scf, n, natom, natom_of_kind, nimg, nkind, norb, nsb, &
         nsbsize, nset, nsgf, nsgf_kind, nspin
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: first_sgf_atom
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
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: fm_wmat
      TYPE(cp_fm_type), POINTER                          :: fm_s_half
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_vhxc
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_p, matrix_s
      TYPE(dbcsr_type)                                   :: sm_q, sm_q_opp, sm_v
      TYPE(dbcsr_type), POINTER                          :: sm_h, sm_p, sm_p_opp, sm_s, sm_vhxc, sm_w
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
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
                      particle_set=particle_set, &
                      rho=rho, &
                      scf_env=scf_env, &
                      para_env=para_env)

      CALL qs_rho_get(rho, rho_ao_kp=matrix_p) ! Density matrices in sparse format

      energy%dft_plus_u = 0.0_dp

      nspin = dft_control%nspins
      nimg = dft_control%nimages
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

      CALL get_particle_set(particle_set, qs_kind_set, &
                            first_sgf=first_sgf_atom)

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

      IF (do_kpoints) THEN
         fm_wmat => scf_env%scf_work1
         fmstruct => fm_wmat(1)%matrix_struct
      ELSE
         ! Retrieve S^(1/2) from the SCF environment
         fm_s_half => scf_env%s_half
         CPASSERT(ASSOCIATED(fm_s_half))
         ! work matrices
         CALL cp_fm_get_info(fm_s_half, matrix_struct=fmstruct)
      END IF
      CALL cp_fm_create(matrix=fm_work1, matrix_struct=fmstruct, &
                        name="FULL WORK MATRIX 1")
      CALL cp_fm_create(matrix=fm_work2, matrix_struct=fmstruct, &
                        name="FULL WORK MATRIX 2")

      ! Calculate S eigenvectors and Lambda matrix for forces
      ! See sTDA forces (get_lowdin_mo_coefficients in qs_tddfpt2_stda_utils
      ! A. Hehn et al JCTC 2022, 18, 4186
      IF (PRESENT(matrix_w)) THEN
         IF (do_kpoints) THEN
            CPABORT("Lowdin forces with k-points NYA in DFT+U")
         END IF
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
         CPABORT("Lowdin option with k-points NYA in DFT+U")
         ALLOCATE (orbq(nsgf, nspin))
         CALL lowdin_kp_trans(kpoints, orbq)
         DEALLOCATE (orbq)
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
            CPABORT("Lowdin option with k-points NYA in DFT+U")
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

            IF (u_minus_j == 0.0_dp .AND. hund_j == 0.0_dp) CYCLE

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
                                    IF (isgf == jsgf) q_matrix(i, j) = q_block(isgf, jsgf)
                                    IF (hund_j /= 0.0_dp .AND. nspin == 2 .AND. ASSOCIATED(q_opp_block)) THEN
                                       IF (isgf == jsgf) q_opp_matrix(i, j) = q_opp_block(isgf, jsgf)
                                    END IF
                                    IF (ASSOCIATED(vhxc_block) .AND. &
                                        dft_control%mtlr_dft_with_perturbation .AND. &
                                        dft_control%mtlr_ikind == ikind .AND. &
                                        qs_kind_set(ikind)%dft_plus_u%lr_atom == atom_a .AND. &
                                        isgf == jsgf) THEN
                                       dft_control%vhxc(ispin) = dft_control%vhxc(ispin) + vhxc_block(isgf, isgf)
                                       dft_control%trq(ispin) = dft_control%trq(ispin) + q_block(isgf, isgf)
                                    END IF
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

                  trq = 0.0_dp
                  trq2 = 0.0_dp
                  trqxq_hund_j = 0.0_dp
                  DO i = 1, n
                     trq = trq + q_matrix(i, i)
                     DO j = 1, n
                        trq2 = trq2 + q_matrix(i, j)*q_matrix(j, i)
                        IF (hund_j /= 0.0_dp) THEN
                           IF (nspin == 2) THEN
                              trqxq_hund_j = trqxq_hund_j + q_matrix(i, j)*q_opp_matrix(j, i)
                           ELSE
                              trqxq_hund_j = trqxq_hund_j + q_matrix(i, j)*q_matrix(j, i)
                           END IF
                        END IF
                     END DO
                  END DO
                  trq = fspin*trq
                  trq2 = fspin*fspin*trq2

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

            IF (dft_control%mtlr_dft_with_perturbation .AND. &
                dft_control%mtlr_ikind == ikind) THEN
               dft_control%vhxc(ispin) = dft_control%vhxc(ispin)/REAL(n, dp)
               CALL para_env%sum(dft_control%vhxc(ispin))
               CALL para_env%sum(dft_control%trq(ispin))
            END IF

            IF (ALLOCATED(q_matrix)) THEN
               DEALLOCATE (q_matrix)
            END IF

            IF (ALLOCATED(q_opp_matrix)) THEN
               DEALLOCATE (q_opp_matrix)
            END IF

         END DO ! next atomic kind "ikind"

         ! Add V(i,j)[U] to V(i,j)[DFT]
         IF (PRESENT(matrix_h)) THEN
            IF (do_kpoints) THEN
               CPABORT("Lowdin option with k-points NYA in DFT+U")
            ELSE
               sm_h => matrix_h(ispin, 1)%matrix
               CALL cp_dbcsr_sm_fm_multiply(sm_v, fm_s_half, fm_work1, nsgf)
               CALL cp_fm_transpose(fm_work1, fm_work2)
               CALL cp_dbcsr_plus_fm_fm_t(sm_h, fm_s_half, fm_work2, nsgf)
            END IF
         END IF ! An update of the Hamiltonian matrix is requested

         ! Calculate the contribution (non-Pulay part) to the derivatives
         ! w.r.t. the nuclear positions
         IF (PRESENT(matrix_w)) THEN

            sm_p => matrix_p(ispin, 1)%matrix
            sm_w => matrix_w(ispin, 1)%matrix

            CALL cp_dbcsr_sm_fm_multiply(sm_v, fm_s_half, fm_work1, nsgf)
            CALL cp_fm_transpose(fm_work1, fm_work2)
            CALL cp_dbcsr_sm_fm_multiply(sm_p, fm_work2, fm_work1, nsgf)
            CALL parallel_gemm('N', 'N', nsgf, nsgf, nsgf, 1.0_dp, fm_work1, fm_sev, 0.0_dp, fm_work2)
            CALL parallel_gemm('T', 'N', nsgf, nsgf, nsgf, 1.0_dp, fm_sev, fm_work2, 0.0_dp, fm_work1)
            CALL cp_fm_schur_product(fm_work1, slambda, fm_work2)
            CALL cp_fm_transpose(fm_work2, fm_work1)
            CALL cp_fm_scale_and_add(alpha=1.0_dp, matrix_a=fm_work1, matrix_b=fm_work2)
            CALL parallel_gemm('N', 'N', nsgf, nsgf, nsgf, -2.0_dp, fm_sev, fm_work2, 0.0_dp, fm_work1)
            CALL cp_dbcsr_plus_fm_fm_t(sm_w, fm_work1, fm_sev, nsgf)

         END IF ! W matrix update requested

      END DO ! next spin "ispin"

      IF (PRESENT(matrix_w)) THEN
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
