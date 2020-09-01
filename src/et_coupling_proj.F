!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief calculates the electron transfer coupling elements by projection-operator approach
!>        Kondov et al. J.Phys.Chem.C 2007, 111, 11970-11981
!> \author Z. Futera (02.2017)
! **************************************************************************************************
MODULE et_coupling_proj

   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE bibliography,                    ONLY: Futera2017,&
                                              cite_reference
   USE cell_types,                      ONLY: cell_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              cp_dbcsr_sm_fm_multiply
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale
   USE cp_fm_diag,                      ONLY: choose_eigv_solver,&
                                              cp_fm_power
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_equivalent,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_to_fm_submat,&
                                              cp_fm_type,&
                                              cp_fm_vectorssum
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE input_section_types,             ONLY: section_get_ivals,&
                                              section_get_lval,&
                                              section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE kpoint_types,                    ONLY: kpoint_type
   USE orbital_pointers,                ONLY: nso
   USE particle_list_types,             ONLY: particle_list_type
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: evolt
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_p_type,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_p_type
   USE qs_collocate_density,            ONLY: calculate_wavefunction
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_mo_methods,                   ONLY: make_mo_eig
   USE qs_mo_occupation,                ONLY: set_mo_occupation
   USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                              deallocate_mo_set,&
                                              mo_set_p_type,&
                                              mo_set_type
   USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                              qs_subsys_type
   USE scf_control_types,               ONLY: scf_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'et_coupling_proj'

   ! Electronic-coupling calculation data structure
   !
   ! n_atoms      - number of atoms in the blocks
   ! n_blocks     - number of atomic blocks (donor,acceptor,bridge,...)
   ! fermi        - system Fermi level (alpha/beta spin component)
   ! m_transf     - transformation matrix for basis-set orthogonalization (S^{-1/2})
   ! m_transf_inv - inversion transformation matrix
   ! block        - atomic data blocks
   TYPE et_cpl
      INTEGER                                            :: n_atoms
      INTEGER                                            :: n_blocks
      REAL(KIND=dp), DIMENSION(:), POINTER               :: fermi
      TYPE(cp_fm_type), POINTER                          :: m_transf
      TYPE(cp_fm_type), POINTER                          :: m_transf_inv
      TYPE(et_cpl_block), DIMENSION(:), POINTER          :: block
   END TYPE et_cpl

   ! Electronic-coupling data block
   !
   ! n_atoms     - number of atoms
   ! n_electrons - number of electrons
   ! n_ao        - number of AO basis functions
   ! atom        - list of atoms
   ! mo          - electronic states
   ! hab         - electronic-coupling elements
   TYPE et_cpl_block
      INTEGER                                            :: n_atoms
      INTEGER                                            :: n_electrons
      INTEGER                                            :: n_ao
      TYPE(et_cpl_atom), DIMENSION(:), POINTER           :: atom
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mo
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER        :: hab
   END TYPE et_cpl_block

   ! Electronic-coupling block-atom data
   ! id     - atom ID
   ! n_ao   - number of AO basis functions
   ! ao_pos - position of atom in array of AO functions
   TYPE et_cpl_atom
      INTEGER                                            :: id
      INTEGER                                            :: n_ao
      INTEGER                                            :: ao_pos
   END TYPE et_cpl_atom

   PUBLIC :: calc_et_coupling_proj

CONTAINS

! **************************************************************************************************
!> \brief Release memory allocate for electronic coupling data structures
!> \param ec electronic coupling data structure
!> \author Z. Futera (02.2017)
! **************************************************************************************************
   SUBROUTINE release_ec_data(ec)

      ! Routine arguments
      TYPE(et_cpl), POINTER                              :: ec

      INTEGER                                            :: i, j, k

! Routine name for debug purposes

      IF (ASSOCIATED(ec)) THEN

         IF (ASSOCIATED(ec%fermi)) &
            DEALLOCATE (ec%fermi)
         IF (ASSOCIATED(ec%m_transf)) &
            CALL cp_fm_release(matrix=ec%m_transf)
         IF (ASSOCIATED(ec%m_transf_inv)) &
            CALL cp_fm_release(matrix=ec%m_transf_inv)

         IF (ASSOCIATED(ec%block)) THEN

            DO i = 1, SIZE(ec%block)
               IF (ASSOCIATED(ec%block(i)%atom)) &
                  DEALLOCATE (ec%block(i)%atom)
               IF (ASSOCIATED(ec%block(i)%mo)) THEN
                  DO j = 1, SIZE(ec%block(i)%mo)
                     IF (ASSOCIATED(ec%block(i)%mo(j)%mo_set)) &
                        CALL deallocate_mo_set(ec%block(i)%mo(j)%mo_set)
                  END DO
                  DEALLOCATE (ec%block(i)%mo)
               END IF
               IF (ASSOCIATED(ec%block(i)%hab)) THEN
                  DO j = 1, SIZE(ec%block(i)%hab, 1)
                     DO k = 1, SIZE(ec%block(i)%hab, 2)
                        IF (ASSOCIATED(ec%block(i)%hab(j, k)%matrix)) &
                           CALL cp_fm_release(matrix=ec%block(i)%hab(j, k)%matrix)
                     END DO
                  END DO
                  DEALLOCATE (ec%block(i)%hab)
               END IF
            END DO

            DEALLOCATE (ec%block)

         END IF

         DEALLOCATE (ec)

      END IF

   END SUBROUTINE release_ec_data

! **************************************************************************************************
!> \brief check the electronic-coupling input section and set the atomic block data
!> \param qs_env QuickStep environment containing all system data
!> \param et_proj_sec the electronic-coupling input section
!> \param ec electronic coupling data structure
!> \author Z. Futera (02.2017)
! **************************************************************************************************
   SUBROUTINE set_block_data(qs_env, et_proj_sec, ec)

      ! Routine arguments
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: et_proj_sec
      TYPE(et_cpl), POINTER                              :: ec

      INTEGER                                            :: i, j, k, l, n, n_ao, n_atoms, n_set
      INTEGER, DIMENSION(:), POINTER                     :: atom_id, atom_nf, atom_ps, n_shell, t
      INTEGER, DIMENSION(:, :), POINTER                  :: ang_mom_id
      LOGICAL                                            :: found
      TYPE(gto_basis_set_type), POINTER                  :: ao_basis_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: block_sec

! Routine name for debug purposes

      NULLIFY (ao_basis_set)
      NULLIFY (particle_set)
      NULLIFY (qs_kind_set)
      NULLIFY (n_shell)
      NULLIFY (ang_mom_id)
      NULLIFY (atom_nf)
      NULLIFY (atom_id)
      NULLIFY (block_sec)

      ! Initialization
      ec%n_atoms = 0
      ec%n_blocks = 0
      NULLIFY (ec%fermi)
      NULLIFY (ec%m_transf)
      NULLIFY (ec%m_transf_inv)
      NULLIFY (ec%block)

      ! Number of atoms / atom types
      CALL get_qs_env(qs_env, particle_set=particle_set, qs_kind_set=qs_kind_set, natom=n_atoms)
      ! Number of AO basis functions
      CALL get_qs_kind_set(qs_kind_set, nsgf=n_ao)

      ! Number of AO functions per atom
      ALLOCATE (atom_nf(n_atoms))
      CPASSERT(ASSOCIATED(atom_nf))

      atom_nf = 0
      DO i = 1, n_atoms
         CALL get_atomic_kind(particle_set(i)%atomic_kind, kind_number=j)
         CALL get_qs_kind(qs_kind_set(j), basis_set=ao_basis_set)
         IF (.NOT. ASSOCIATED(ao_basis_set)) &
            CPABORT('Unsupported basis set type. ')
         CALL get_gto_basis_set(gto_basis_set=ao_basis_set, &
                                nset=n_set, nshell=n_shell, l=ang_mom_id)
         DO j = 1, n_set
            DO k = 1, n_shell(j)
               atom_nf(i) = atom_nf(i) + nso(ang_mom_id(k, j))
            END DO
         END DO
      END DO

      ! Sanity check
      n = 0
      DO i = 1, n_atoms
         n = n + atom_nf(i)
      END DO
      CPASSERT(n == n_ao)

      ! Atom position in AO array
      ALLOCATE (atom_ps(n_atoms))
      CPASSERT(ASSOCIATED(atom_ps))
      atom_ps = 1
      DO i = 1, n_atoms - 1
         atom_ps(i + 1) = atom_ps(i) + atom_nf(i)
      END DO

      ! Number of blocks
      block_sec => section_vals_get_subs_vals(et_proj_sec, 'BLOCK')
      CALL section_vals_get(block_sec, n_repetition=ec%n_blocks)
      ALLOCATE (ec%block(ec%n_blocks))
      CPASSERT(ASSOCIATED(ec%block))

      ! Block data
      ALLOCATE (t(n_atoms))
      CPASSERT(ASSOCIATED(t))

      ec%n_atoms = 0
      DO i = 1, ec%n_blocks

         ! Initialization
         ec%block(i)%n_atoms = 0
         ec%block(i)%n_electrons = 0
         ec%block(i)%n_ao = 0
         NULLIFY (ec%block(i)%atom)
         NULLIFY (ec%block(i)%mo)
         NULLIFY (ec%block(i)%hab)

         ! Number of electrons
         CALL section_vals_val_get(block_sec, i_rep_section=i, &
                                   keyword_name='NELECTRON', i_val=ec%block(i)%n_electrons)

         ! User-defined atom array
         CALL section_vals_val_get(block_sec, i_rep_section=i, &
                                   keyword_name='ATOMS', i_vals=atom_id)

         ! Count unique atoms
         DO j = 1, SIZE(atom_id)
            ! Check atom ID validity
            IF (atom_id(j) < 1 .OR. atom_id(j) > n_atoms) &
               CPABORT('invalid fragment atom ID ('//TRIM(ADJUSTL(cp_to_string(atom_id(j))))//')')
            ! Check if the atom is not in previously-defined blocks
            found = .FALSE.
            DO k = 1, i - 1
               DO l = 1, ec%block(k)%n_atoms
                  IF (ec%block(k)%atom(l)%id == atom_id(j)) THEN
                     CPWARN('multiple definition of atom'//TRIM(ADJUSTL(cp_to_string(atom_id(j)))))
                     found = .TRUE.
                     EXIT
                  END IF
               END DO
            END DO
            ! Check if the atom is not in already defined in the present block
            IF (.NOT. found) THEN
               DO k = 1, ec%block(i)%n_atoms
                  IF (t(k) == atom_id(j)) THEN
                     CPWARN('multiple definition of atom'//TRIM(ADJUSTL(cp_to_string(atom_id(j)))))
                     found = .TRUE.
                     EXIT
                  END IF
               END DO
            END IF
            ! Save the atom
            IF (.NOT. found) THEN
               ec%block(i)%n_atoms = ec%block(i)%n_atoms + 1
               t(ec%block(i)%n_atoms) = atom_id(j)
            END IF
         END DO

         ! Memory allocation
         ALLOCATE (ec%block(i)%atom(ec%block(i)%n_atoms))
         CPASSERT(ASSOCIATED(ec%block(i)%atom))

         ! Save atom IDs and number of AOs
         DO j = 1, ec%block(i)%n_atoms
            ec%block(i)%atom(j)%id = t(j)
            ec%block(i)%atom(j)%n_ao = atom_nf(ec%block(i)%atom(j)%id)
            ec%block(i)%atom(j)%ao_pos = atom_ps(ec%block(i)%atom(j)%id)
            ec%block(i)%n_ao = ec%block(i)%n_ao + ec%block(i)%atom(j)%n_ao
         END DO

         ec%n_atoms = ec%n_atoms + ec%block(i)%n_atoms
      END DO

      ! Clean memory
      IF (ASSOCIATED(atom_nf)) &
         DEALLOCATE (atom_nf)
      IF (ASSOCIATED(atom_ps)) &
         DEALLOCATE (atom_ps)
      IF (ASSOCIATED(t)) &
         DEALLOCATE (t)

   END SUBROUTINE set_block_data

! **************************************************************************************************
!> \brief check the electronic-coupling input section and set the atomic block data
!> \param ec electronic coupling data structure
!> \param fa system Fermi level (alpha spin)
!> \param fb system Fermi level (beta spin)
!> \author Z. Futera (02.2017)
! **************************************************************************************************
   SUBROUTINE set_fermi(ec, fa, fb)

      ! Routine arguments
      TYPE(et_cpl), POINTER                              :: ec
      REAL(KIND=dp)                                      :: fa
      REAL(KIND=dp), OPTIONAL                            :: fb

! Routine name for debug purposes

      NULLIFY (ec%fermi)

      IF (PRESENT(fb)) THEN

         ALLOCATE (ec%fermi(2))
         CPASSERT(ASSOCIATED(ec%fermi))
         ec%fermi(1) = fa
         ec%fermi(2) = fb

      ELSE

         ALLOCATE (ec%fermi(1))
         CPASSERT(ASSOCIATED(ec%fermi))
         ec%fermi(1) = fa

      END IF

   END SUBROUTINE set_fermi

! **************************************************************************************************
!> \brief copy full-matrix elements to real array on local core
!> \param mat the full matrix
!> \param arr the array for the matrix elements
!> \author Z. Futera (02.2017)
! **************************************************************************************************
   SUBROUTINE get_fm_matrix_array(mat, arr)

      IMPLICIT NONE

      ! Routine arguments
      TYPE(cp_fm_type), POINTER                          :: mat
      REAL(KIND=dp), DIMENSION(:, :), POINTER             :: arr

      ! Routine name for debug purposes

      ! Local variables
      INTEGER                                            :: i, j

#if defined(__SCALAPACK)

      INTEGER                                            :: c_row, c_col, c_row_f, c_col_f
      INTEGER                                            :: i_p_row, i_p_col
      INTEGER                                            :: n_p_rows, n_p_cols

      n_p_rows = mat%matrix_struct%context%num_pe(1)
      n_p_cols = mat%matrix_struct%context%num_pe(2)

      c_row = mat%matrix_struct%context%mepos(1)
      c_col = mat%matrix_struct%context%mepos(2)

      DO i = 1, mat%matrix_struct%nrow_global
         DO j = 1, mat%matrix_struct%ncol_global

            CALL infog2l(i, j, mat%matrix_struct%descriptor, &
                         n_p_rows, n_p_cols, &
                         c_row, c_col, &
                         i_p_row, i_p_col, &
                         c_row_f, c_col_f)

            IF ((c_row_f == c_row) .AND. (c_col_f == c_col)) THEN
               arr(i, j) = mat%local_data(i_p_row, i_p_col)
               CALL dgebs2d(mat%matrix_struct%context%group, 'All', ' ', 1, 1, arr(i, j), 1)
            ELSE
               CALL dgebr2d(mat%matrix_struct%context%group, 'All', ' ', 1, 1, arr(i, j), 1, c_row_f, c_col_f)
            END IF

         END DO
      END DO

#else

      DO i = 1, mat%matrix_struct%nrow_global
         DO j = 1, mat%matrix_struct%ncol_global
            arr(i, j) = mat%local_data(i, j)
         END DO
      END DO

#endif

   END SUBROUTINE get_fm_matrix_array

! **************************************************************************************************
!> \brief reorder Hamiltonian matrix according to defined atomic blocks
!> \param ec electronic coupling data structure
!> \param mat_h the Hamiltonian matrix
!> \param mat_w working matrix of the same dimension
!> \author Z. Futera (02.2017)
! **************************************************************************************************
   SUBROUTINE reorder_hamiltonian_matrix(ec, mat_h, mat_w)

      IMPLICIT NONE

      ! Routine arguments
      TYPE(et_cpl), POINTER                              :: ec
      TYPE(cp_fm_type), POINTER                          :: mat_h
      TYPE(cp_fm_type), POINTER                          :: mat_w

      ! Routine name for debug purposes

      ! Local variables
      INTEGER                                            :: ir, ic, jr, jc, kr, kc, nr, nc
#if defined(__SCALAPACK)
      INTEGER                                            :: iw_p_row, iw_p_col, ih_p_row, ih_p_col
      INTEGER                                            :: nw_p_rows, nw_p_cols, nh_p_rows, nh_p_cols
      INTEGER                                            :: cw_row, cw_col, ch_row, ch_col
      INTEGER                                            :: cw_row_f, cw_col_f, ch_row_f, ch_col_f
      REAL(KIND=dp)                                      :: xh

      IF (.NOT. cp_fm_struct_equivalent(mat_h%matrix_struct, mat_w%matrix_struct)) &
         CPABORT('cannot reorder Hamiltonian, working-matrix structure is not equivalent')

      ! number of processors
      nw_p_rows = mat_w%matrix_struct%context%num_pe(1)
      nw_p_cols = mat_w%matrix_struct%context%num_pe(2)
      nh_p_rows = mat_h%matrix_struct%context%num_pe(1)
      nh_p_cols = mat_h%matrix_struct%context%num_pe(2)

      ! position of processors
      cw_row = mat_w%matrix_struct%context%mepos(1)
      cw_col = mat_w%matrix_struct%context%mepos(2)
      ch_row = mat_h%matrix_struct%context%mepos(1)
      ch_col = mat_h%matrix_struct%context%mepos(2)
#endif

      ! Matrix-element reordering
      nr = 1
      ! Rows
      DO ir = 1, ec%n_blocks
         DO jr = 1, ec%block(ir)%n_atoms
            DO kr = 1, ec%block(ir)%atom(jr)%n_ao
               ! Columns
               nc = 1
               DO ic = 1, ec%n_blocks
                  DO jc = 1, ec%block(ic)%n_atoms
                     DO kc = 1, ec%block(ic)%atom(jc)%n_ao
#if defined(__SCALAPACK)
                        CALL infog2l(nr, nc, mat_w%matrix_struct%descriptor, &
                                     nw_p_rows, nw_p_cols, cw_row, cw_col, &
                                     iw_p_row, iw_p_col, cw_row_f, cw_col_f)
                        CALL infog2l( &
                           ec%block(ir)%atom(jr)%ao_pos + kr - 1, &
                           ec%block(ic)%atom(jc)%ao_pos + kc - 1, &
                           mat_h%matrix_struct%descriptor, &
                           nh_p_rows, nh_p_cols, ch_row, ch_col, &
                           ih_p_row, ih_p_col, ch_row_f, ch_col_f)
                        ! Local H element
                        IF ((ch_row_f == ch_row) .AND. (ch_col_f == ch_col)) THEN
                           xh = mat_h%local_data(ih_p_row, ih_p_col)
                           CALL dgebs2d(mat_h%matrix_struct%context%group, 'All', ' ', 1, 1, xh, 1)
                           IF ((cw_row_f == cw_row) .AND. (cw_col_f == cw_col)) THEN
                              mat_w%local_data(iw_p_row, iw_p_col) = xh
                           END IF
                           ! Remote H element
                        ELSE
                           CALL dgebr2d(mat_h%matrix_struct%context%group, &
                                        'All', ' ', 1, 1, xh, 1, ch_row_f, ch_col_f)
                           IF ((cw_row_f == cw_row) .AND. (cw_col_f == cw_col)) THEN
                              mat_w%local_data(iw_p_row, iw_p_col) = xh
                           END IF
                        END IF
#else
                        mat_w%local_data(nr, nc) = mat_h%local_data( &
                                                   ec%block(ir)%atom(jr)%ao_pos + kr - 1, &
                                                   ec%block(ic)%atom(jc)%ao_pos + kc - 1)
#endif
                        nc = nc + 1
                     END DO
                  END DO
               END DO
               nr = nr + 1
            END DO
         END DO
      END DO

      ! Copy the reordered matrix to original data array
      CALL cp_fm_to_fm(mat_w, mat_h)

   END SUBROUTINE reorder_hamiltonian_matrix

! **************************************************************************************************
!> \brief calculated transformation matrix for basis-set orthogonalization (S^{-1/2})
!> \param qs_env QuickStep environment containing all system data
!> \param mat_t storage for the transformation matrix
!> \param mat_i storage for the inversion transformation matrix
!> \param mat_w working matrix of the same dimension
!> \author Z. Futera (02.2017)
! **************************************************************************************************
   SUBROUTINE get_s_half_inv_matrix(qs_env, mat_t, mat_i, mat_w)

      ! Routine arguments
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), POINTER                          :: mat_t, mat_i, mat_w

      INTEGER                                            :: n_deps
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_s
      TYPE(scf_control_type), POINTER                    :: scf_cntrl

! Routine name for debug purposes

      NULLIFY (mat_s)
      NULLIFY (scf_cntrl)

      CALL get_qs_env(qs_env, matrix_s=mat_s)
      CALL copy_dbcsr_to_fm(mat_s(1)%matrix, mat_t)
      CALL copy_dbcsr_to_fm(mat_s(1)%matrix, mat_i)

      ! Transformation S -> S^{-1/2}
      CALL get_qs_env(qs_env, scf_control=scf_cntrl)
      CALL cp_fm_power(mat_t, mat_w, -0.5_dp, scf_cntrl%eps_eigval, n_deps)
      CALL cp_fm_power(mat_i, mat_w, +0.5_dp, scf_cntrl%eps_eigval, n_deps)
      ! Sanity check
      IF (n_deps /= 0) THEN
         CALL cp_warn(__LOCATION__, &
                      "Overlap matrix exhibits linear dependencies. At least some "// &
                      "eigenvalues have been quenched.")
      END IF

   END SUBROUTINE get_s_half_inv_matrix

! **************************************************************************************************
!> \brief transform KS hamiltonian to orthogonalized block-separated basis set
!> \param qs_env QuickStep environment containing all system data
!> \param ec electronic coupling data structure
!> \param fm_s full-matrix structure used for allocation of KS matrices
!> \param mat_t storage for pointers to the transformed KS matrices
!> \param mat_w working matrix of the same dimension
!> \param n_ao total number of AO basis functions
!> \param n_spins number of spin components
!> \author Z. Futera (02.2017)
! **************************************************************************************************
   SUBROUTINE get_block_hamiltonian(qs_env, ec, fm_s, mat_t, mat_w, n_ao, n_spins)

      ! Routine arguments
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(et_cpl), POINTER                              :: ec
      TYPE(cp_fm_struct_type), POINTER                   :: fm_s
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: mat_t
      TYPE(cp_fm_type), POINTER                          :: mat_w
      INTEGER                                            :: n_ao, n_spins

      INTEGER                                            :: i
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_h

! Routine name for debug purposes

      NULLIFY (mat_h)

      ! Memory allocation
      ALLOCATE (mat_t(n_spins))
      CPASSERT(ASSOCIATED(mat_t))

      ! KS Hamiltonian
      CALL get_qs_env(qs_env, matrix_ks=mat_h)
      ! Transformation matrix
      CALL cp_fm_create(matrix=ec%m_transf, matrix_struct=fm_s, &
                        name='S^(-1/2) TRANSFORMATION MATRIX')
      CALL cp_fm_create(matrix=ec%m_transf_inv, matrix_struct=fm_s, &
                        name='S^(+1/2) TRANSFORMATION MATRIX')
      CALL get_s_half_inv_matrix(qs_env, ec%m_transf, ec%m_transf_inv, mat_w)

      DO i = 1, n_spins

         ! Full-matrix format
         CALL cp_fm_create(matrix=mat_t(i)%matrix, matrix_struct=fm_s, &
                           name='KS HAMILTONIAN IN SEPARATED ORTHOGONALIZED BASIS SET')
         CALL copy_dbcsr_to_fm(mat_h(i)%matrix, mat_t(i)%matrix)

         ! Transform KS Hamiltonian to the orthogonalized AO basis set
         CALL cp_gemm("N", "N", n_ao, n_ao, n_ao, 1.0_dp, ec%m_transf, mat_t(i)%matrix, 0.0_dp, mat_w)
         CALL cp_gemm("N", "N", n_ao, n_ao, n_ao, 1.0_dp, mat_w, ec%m_transf, 0.0_dp, mat_t(i)%matrix)

         ! Reorder KS Hamiltonain elements to defined block structure
         CALL reorder_hamiltonian_matrix(ec, mat_t(i)%matrix, mat_w)

      END DO

   END SUBROUTINE get_block_hamiltonian

! **************************************************************************************************
!> \brief Diagonalize diagonal blocks of the KS hamiltonian in separated orthogonalized basis set
!> \param qs_env QuickStep environment containing all system data
!> \param ec electronic coupling data structure
!> \param mat_h Hamiltonian in separated orthogonalized basis set
!> \author Z. Futera (02.2017)
! **************************************************************************************************
   SUBROUTINE hamiltonian_block_diag(qs_env, ec, mat_h)

      ! Routine arguments
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(et_cpl), POINTER                              :: ec
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: mat_h

      INTEGER                                            :: i, j, k, l, n_spins, spin
      REAL(KIND=dp), DIMENSION(:), POINTER               :: vec_e
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: dat
      TYPE(cp_fm_struct_type), POINTER                   :: fm_s
      TYPE(cp_fm_type), POINTER                          :: mat_u
      TYPE(cp_para_env_type), POINTER                    :: para_env

! Routine name for debug purposes

      NULLIFY (vec_e)
      NULLIFY (blacs_env)
      NULLIFY (para_env)
      NULLIFY (fm_s)
      NULLIFY (dat)
      NULLIFY (mat_u)

      ! Parallel environment
      CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env)

      ! Storage for block sub-matrices
      ALLOCATE (dat(ec%n_blocks))
      CPASSERT(ASSOCIATED(dat))

      ! Storage for electronic states and couplings
      n_spins = SIZE(mat_h)
      DO i = 1, ec%n_blocks
         ALLOCATE (ec%block(i)%mo(n_spins))
         CPASSERT(ASSOCIATED(ec%block(i)%mo))
         ALLOCATE (ec%block(i)%hab(n_spins, ec%n_blocks))
         CPASSERT(ASSOCIATED(ec%block(i)%hab))
         DO j = 1, n_spins
            NULLIFY (ec%block(i)%mo(j)%mo_set)
            DO k = 1, ec%n_blocks
               NULLIFY (ec%block(i)%hab(j, k)%matrix)
            END DO
         END DO
      END DO

      ! Spin components
      DO spin = 1, n_spins

         ! Diagonal blocks
         j = 1
         DO i = 1, ec%n_blocks

            ! Memory allocation
            CALL cp_fm_struct_create(fmstruct=fm_s, para_env=para_env, context=blacs_env, &
                                     nrow_global=ec%block(i)%n_ao, ncol_global=ec%block(i)%n_ao)
            CALL cp_fm_create(matrix=dat(i)%matrix, matrix_struct=fm_s, &
                              name='H_KS DIAGONAL BLOCK')

            ALLOCATE (vec_e(ec%block(i)%n_ao))
            CPASSERT(ASSOCIATED(vec_e))

            ! Copy block data
            CALL cp_fm_to_fm_submat(mat_h(spin)%matrix, &
                                    dat(i)%matrix, ec%block(i)%n_ao, &
                                    ec%block(i)%n_ao, j, j, 1, 1)

            ! Diagonalization
            CALL cp_fm_create(matrix=mat_u, matrix_struct=fm_s, name='UNITARY MATRIX')
            CALL choose_eigv_solver(dat(i)%matrix, mat_u, vec_e)
            CALL cp_fm_to_fm(mat_u, dat(i)%matrix)

            ! Save state energies / vectors
            CALL create_block_mo_set(qs_env, ec, i, spin, mat_u, vec_e)

            ! Clean memory
            CALL cp_fm_struct_release(fmstruct=fm_s)
            CALL cp_fm_release(matrix=mat_u)
            DEALLOCATE (vec_e)

            ! Off-set for next block
            j = j + ec%block(i)%n_ao

         END DO

         ! Off-diagonal blocks
         k = 1
         DO i = 1, ec%n_blocks
            l = 1
            DO j = 1, ec%n_blocks
               IF (i /= j) THEN

                  ! Memory allocation
                  CALL cp_fm_struct_create(fmstruct=fm_s, para_env=para_env, context=blacs_env, &
                                           nrow_global=ec%block(i)%n_ao, ncol_global=ec%block(j)%n_ao)
                  CALL cp_fm_create(matrix=ec%block(i)%hab(spin, j)%matrix, matrix_struct=fm_s, &
                                    name='H_KS OFF-DIAGONAL BLOCK')

                  ! Copy block data
                  CALL cp_fm_to_fm_submat(mat_h(spin)%matrix, &
                                          ec%block(i)%hab(spin, j)%matrix, ec%block(i)%n_ao, &
                                          ec%block(j)%n_ao, k, l, 1, 1)

                  ! Transformation
                  CALL cp_fm_create(matrix=mat_u, matrix_struct=fm_s, name='FULL WORK MATRIX')
                  CALL cp_gemm("T", "N", ec%block(i)%n_ao, ec%block(j)%n_ao, ec%block(i)%n_ao, &
                               1.0_dp, dat(i)%matrix, ec%block(i)%hab(spin, j)%matrix, 0.0_dp, mat_u)
                  CALL cp_gemm("N", "N", ec%block(i)%n_ao, ec%block(j)%n_ao, ec%block(j)%n_ao, &
                               1.0_dp, mat_u, dat(j)%matrix, 0.0_dp, ec%block(i)%hab(spin, j)%matrix)

                  ! Clean memory
                  CALL cp_fm_struct_release(fmstruct=fm_s)
                  CALL cp_fm_release(matrix=mat_u)

               END IF
               ! Off-set for next block
               l = l + ec%block(j)%n_ao
            END DO
            ! Off-set for next block
            k = k + ec%block(i)%n_ao
         END DO

         ! Clean memory
         IF (ASSOCIATED(dat)) THEN
            DO i = 1, SIZE(dat)
               IF (ASSOCIATED(dat(i)%matrix)) &
                  CALL cp_fm_release(matrix=dat(i)%matrix)
            END DO
         END IF
      END DO

      ! Clean memory
      IF (ASSOCIATED(dat)) &
         DEALLOCATE (dat)

   END SUBROUTINE hamiltonian_block_diag

! **************************************************************************************************
!> \brief Return sum of selected squared MO coefficients
!> \param blk_at list of atoms in the block
!> \param mo array of MO sets
!> \param id state index
!> \param atom list of atoms for MO coefficient summing
!> \return ...
!> \author Z. Futera (02.2017)
! **************************************************************************************************
   FUNCTION get_mo_c2_sum(blk_at, mo, id, atom) RESULT(c2)

      IMPLICIT NONE

      ! Routine arguments
      TYPE(et_cpl_atom), DIMENSION(:), POINTER           :: blk_at
      TYPE(cp_fm_type), POINTER                          :: mo
      INTEGER                                            :: id
      INTEGER, DIMENSION(:), POINTER                     :: atom

      ! Returning value
      REAL(KIND=dp)                                      :: c2

      ! Routine name for debug purposes

      ! Local variables
      LOGICAL                                            :: found
      INTEGER                                            :: i, j, k
      REAL(KIND=dp)                                      :: c
#if defined(__SCALAPACK)
      INTEGER                                            :: c_row, c_col, c_row_f, c_col_f
      INTEGER                                            :: i_p_row, i_p_col
      INTEGER                                            :: n_p_rows, n_p_cols

      ! number of processors
      n_p_rows = mo%matrix_struct%context%num_pe(1)
      n_p_cols = mo%matrix_struct%context%num_pe(2)

      ! position of processors
      c_row = mo%matrix_struct%context%mepos(1)
      c_col = mo%matrix_struct%context%mepos(2)
#endif

      ! initialization
      c2 = 0.0d0

      ! selected atoms
      DO i = 1, SIZE(atom)

         ! find atomic function offset
         found = .FALSE.
         DO j = 1, SIZE(blk_at)
            IF (blk_at(j)%id == atom(i)) THEN
               found = .TRUE.
               EXIT
            END IF
         END DO

         IF (.NOT. found) &
            CPABORT('MO-fraction atom ID not defined in the block')

         ! sum MO coefficients from the atom
         DO k = 1, blk_at(j)%n_ao
#if defined(__SCALAPACK)
            c = 0.0_dp
            CALL infog2l(blk_at(j)%ao_pos + k - 1, id, mo%matrix_struct%descriptor, &
                         n_p_rows, n_p_cols, c_row, c_col, &
                         i_p_row, i_p_col, c_row_f, c_col_f)
            ! local element
            IF ((c_row_f == c_row) .AND. (c_col_f == c_col)) THEN
               c = mo%local_data(i_p_row, i_p_col)
               CALL dgebs2d(mo%matrix_struct%context%group, 'All', ' ', 1, 1, c, 1)
               ! remote element
            ELSE
               CALL dgebr2d(mo%matrix_struct%context%group, 'All', ' ', 1, 1, c, 1, c_row_f, c_col_f)
            END IF
#else
            c = mo%local_data(blk_at(j)%ao_pos + k - 1, id)
#endif
            c2 = c2 + c*c
         END DO

      END DO

   END FUNCTION get_mo_c2_sum

! **************************************************************************************************
!> \brief Print out specific MO coefficients
!> \param output_unit unit number of the open output stream
!> \param qs_env QuickStep environment containing all system data
!> \param ec electronic coupling data structure
!> \param blk atomic-block ID
!> \param n_spins number of spin components
!> \author Z. Futera (02.2017)
! **************************************************************************************************
   SUBROUTINE print_mo_coeff(output_unit, qs_env, ec, blk, n_spins)

      ! Routine arguments
      INTEGER                                            :: output_unit
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(et_cpl), POINTER                              :: ec
      INTEGER                                            :: blk, n_spins

      INTEGER                                            :: j, k, l, m, n, n_ao, n_mo
      INTEGER, DIMENSION(:), POINTER                     :: list_at, list_mo
      REAL(KIND=dp)                                      :: c1, c2
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: mat_w
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: block_sec, print_sec

! Routine name for debug purposes

      NULLIFY (block_sec)
      NULLIFY (print_sec)
      NULLIFY (qs_kind_set)

      ! Atomic block data
      block_sec => section_vals_get_subs_vals(qs_env%input, &
                                              'PROPERTIES%ET_COUPLING%PROJECTION%BLOCK')

      print_sec => section_vals_get_subs_vals(block_sec, 'PRINT', i_rep_section=blk)

      ! List of atoms
      CALL section_vals_val_get(print_sec, keyword_name='MO_COEFF_ATOM', n_rep_val=n)

      IF (n > 0) THEN

         IF (output_unit > 0) &
            WRITE (output_unit, '(/,T3,A/)') 'Block state fractions:'

         ! Number of AO functions
         CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set)
         CALL get_qs_kind_set(qs_kind_set, nsgf=n_ao)

         ! MOs in orthonormal basis set
         ALLOCATE (mat_w(n_spins))
         CPASSERT(ASSOCIATED(mat_w))
         DO j = 1, n_spins
            n_mo = ec%block(blk)%n_ao
            CALL cp_fm_create(matrix=mat_w(j)%matrix, &
                              matrix_struct=ec%block(blk)%mo(j)%mo_set%mo_coeff%matrix_struct, &
                              name='BLOCK MOs IN ORTHONORMAL BASIS SET')
            CALL cp_gemm("N", "N", n_ao, n_mo, n_ao, 1.0_dp, ec%m_transf_inv, &
                         ec%block(blk)%mo(j)%mo_set%mo_coeff, 0.0_dp, mat_w(j)%matrix)
         END DO

         DO j = 1, n
            NULLIFY (list_at)
            CALL section_vals_val_get(print_sec, keyword_name='MO_COEFF_ATOM', &
                                      i_rep_val=j, i_vals=list_at)
            IF (ASSOCIATED(list_at)) THEN

               ! List of states
               CALL section_vals_val_get(print_sec, keyword_name='MO_COEFF_ATOM_STATE', n_rep_val=m)

               IF (m > 0) THEN

                  DO k = 1, m
                     NULLIFY (list_mo)
                     CALL section_vals_val_get(print_sec, keyword_name='MO_COEFF_ATOM_STATE', &
                                               i_rep_val=k, i_vals=list_mo)
                     IF (ASSOCIATED(list_mo)) THEN

                        IF (j > 1) THEN
                           IF (output_unit > 0) &
                              WRITE (output_unit, *)
                        END IF

                        DO l = 1, SIZE(list_mo)

                           IF (n_spins > 1) THEN
                              c1 = get_mo_c2_sum(ec%block(blk)%atom, mat_w(1)%matrix, &
                                                 list_mo(l), list_at)
                              c2 = get_mo_c2_sum(ec%block(blk)%atom, mat_w(2)%matrix, &
                                                 list_mo(l), list_at)
                              IF (output_unit > 0) &
                                 WRITE (output_unit, '(I5,A,I5,2F20.10)') j, ' /', list_mo(l), c1, c2
                           ELSE
                              c1 = get_mo_c2_sum(ec%block(blk)%atom, mat_w(1)%matrix, &
                                                 list_mo(l), list_at)
                              IF (output_unit > 0) &
                                 WRITE (output_unit, '(I5,A,I5,F20.10)') j, ' /', list_mo(l), c1
                           END IF

                        END DO

                     END IF
                  END DO

               END IF

            END IF
         END DO

         ! Clean memory
         DO j = 1, n_spins
            CALL cp_fm_release(matrix=mat_w(j)%matrix)
         END DO
         DEALLOCATE (mat_w)

      END IF

   END SUBROUTINE print_mo_coeff

! **************************************************************************************************
!> \brief Print out electronic states (MOs)
!> \param output_unit unit number of the open output stream
!> \param mo array of MO sets
!> \param n_spins number of spin components
!> \param label output label
!> \param mx_mo_a maximum number of alpha states to print out
!> \param mx_mo_b maximum number of beta states to print out
!> \param fermi print out Fermi level and number of electrons
!> \author Z. Futera (02.2017)
! **************************************************************************************************
   SUBROUTINE print_states(output_unit, mo, n_spins, label, mx_mo_a, mx_mo_b, fermi)

      ! Routine arguments
      INTEGER                                            :: output_unit
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mo
      INTEGER                                            :: n_spins
      CHARACTER(LEN=*)                                   :: label
      INTEGER, OPTIONAL                                  :: mx_mo_a, mx_mo_b
      LOGICAL, OPTIONAL                                  :: fermi

      INTEGER                                            :: i, mx_a, mx_b, n
      LOGICAL                                            :: prnt_fm

! Routine name for debug purposes

      prnt_fm = .FALSE.
      IF (PRESENT(fermi)) &
         prnt_fm = fermi

      IF (output_unit > 0) THEN

         WRITE (output_unit, '(/,T3,A/)') 'State energies ('//TRIM(ADJUSTL(label))//'):'

         ! Spin-polarized calculation
         IF (n_spins > 1) THEN

            mx_a = mo(1)%mo_set%nmo
            IF (PRESENT(mx_mo_a)) &
               mx_a = MIN(mo(1)%mo_set%nmo, mx_mo_a)
            mx_b = mo(2)%mo_set%nmo
            IF (PRESENT(mx_mo_b)) &
               mx_b = MIN(mo(2)%mo_set%nmo, mx_mo_b)
            n = MAX(mx_a, mx_b)

            DO i = 1, n
               WRITE (output_unit, '(T3,I10)', ADVANCE='no') i
               IF (i <= mx_a) THEN
                  WRITE (output_unit, '(2F12.4)', ADVANCE='no') &
                     mo(1)%mo_set%occupation_numbers(i), mo(1)%mo_set%eigenvalues(i)
               ELSE
                  WRITE (output_unit, '(A)', ADVANCE='no') '                        '
               END IF
               WRITE (output_unit, '(A)', ADVANCE='no') '     '
               IF (i <= mx_b) THEN
                  WRITE (output_unit, '(2F12.4)') &
                     mo(2)%mo_set%occupation_numbers(i), mo(2)%mo_set%eigenvalues(i)
               ELSE
                  WRITE (output_unit, *)
               END IF
            END DO

            IF (prnt_fm) THEN
               WRITE (output_unit, '(/,T3,I10,F24.4,I10,F19.4)') &
                  mo(1)%mo_set%nelectron, mo(1)%mo_set%mu, &
                  mo(2)%mo_set%nelectron, mo(2)%mo_set%mu
            END IF

            ! Spin-restricted calculation
         ELSE

            mx_a = mo(1)%mo_set%nmo
            IF (PRESENT(mx_mo_a)) &
               mx_a = MIN(mo(1)%mo_set%nmo, mx_mo_a)

            DO i = 1, mx_a
               WRITE (output_unit, '(T3,I10,2F12.4)') &
                  i, mo(1)%mo_set%occupation_numbers(i), mo(1)%mo_set%eigenvalues(i)
            END DO

            IF (prnt_fm) THEN
               WRITE (output_unit, '(/,T3,I10,F24.4)') &
                  mo(1)%mo_set%nelectron, mo(1)%mo_set%mu
            END IF

         END IF

      END IF

   END SUBROUTINE print_states

! **************************************************************************************************
!> \brief Print out donor-acceptor state couplings
!> \param ec_sec ...
!> \param output_unit unit number of the open output stream
!> \param logger ...
!> \param ec electronic coupling data structure
!> \param mo ...
!> \author Z. Futera (02.2017)
! **************************************************************************************************
   SUBROUTINE print_couplings(ec_sec, output_unit, logger, ec, mo)

      IMPLICIT NONE

      ! Routine arguments
      TYPE(section_vals_type), POINTER                   :: ec_sec
      INTEGER                                            :: output_unit
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(et_cpl), POINTER                              :: ec
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mo

      ! Routine name for debug purposes

      ! Local variables
      CHARACTER(LEN=default_path_length)                 :: filename, my_pos, title
      INTEGER                                            :: i, j, k, l, nspins, unit_nr
      INTEGER                                            :: n_states(2)
#if defined(__SCALAPACK)
      INTEGER                                            :: nr, nc
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: w1, w2
#endif
      LOGICAL                                            :: append
      TYPE(section_vals_type), POINTER                   :: print_key

      n_states = 0
      DO i = 1, SIZE(mo)
         n_states(i) = mo(i)%mo_set%nmo
      END DO
      nspins = 1
      IF (n_states(2) > 0) nspins = 2

      print_key => section_vals_get_subs_vals(section_vals=ec_sec, &
                                              subsection_name="PRINT%COUPLINGS")

      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), &
                cp_p_file)) THEN

         my_pos = "REWIND"
         append = section_get_lval(print_key, "APPEND")
         IF (append) THEN
            my_pos = "APPEND"
         END IF

         IF (output_unit > 0) &
            WRITE (output_unit, '(/,T3,A/)') 'Printing coupling elements to output files'

#if defined(__SCALAPACK)

         DO i = 1, ec%n_blocks
            DO j = i + 1, ec%n_blocks

               nr = ec%block(i)%hab(1, j)%matrix%matrix_struct%nrow_global
               nc = ec%block(i)%hab(1, j)%matrix%matrix_struct%ncol_global

               ALLOCATE (w1(nr, nc))
               CPASSERT(ASSOCIATED(w1))
               CALL get_fm_matrix_array(ec%block(i)%hab(1, j)%matrix, w1)
               IF (nspins > 1) THEN
                  ALLOCATE (w2(nr, nc))
                  CPASSERT(ASSOCIATED(w2))
                  CALL get_fm_matrix_array(ec%block(i)%hab(2, j)%matrix, w2)
               END IF

               IF (output_unit > 0) THEN

                  WRITE (filename, '(a5,I1.1,a1,I1.1)') "ET_BL_", i, "-", j
                  unit_nr = cp_print_key_unit_nr(logger, ec_sec, "PRINT%COUPLINGS", extension=".elcoup", &
                                                 middle_name=TRIM(filename), file_position=my_pos, log_filename=.FALSE.)

                  WRITE (title, *) 'Coupling elements [meV] between blocks:', i, j

                  WRITE (unit_nr, *) TRIM(title)
                  IF (nspins > 1) THEN
                     WRITE (unit_nr, '(T3,A8,T13,A8,T28,A,A)') 'State A', 'State B', 'Coupling spin 1', '  Coupling spin 2'
                  ELSE
                     WRITE (unit_nr, '(T3,A8,T13,A8,T28,A)') 'State A', 'State B', 'Coupling'
                  END IF

                  DO k = 1, MIN(ec%block(i)%n_ao, n_states(1))
                     DO l = 1, MIN(ec%block(j)%n_ao, n_states(1))

                        IF (nspins > 1) THEN

                           WRITE (unit_nr, '(T3,I5,T13,I5,T22,E20.6)', ADVANCE='no') &
                              k, l, w1(k, l)*evolt*1000.0_dp
                           IF ((k <= n_states(2)) .AND. (l <= n_states(2))) THEN
                              WRITE (unit_nr, '(E20.6)') &
                                 w2(k, l)*evolt*1000.0_dp
                           ELSE
                              WRITE (unit_nr, *)
                           END IF

                        ELSE

                           WRITE (unit_nr, '(T3,I5,T13,I5,T22,E20.6)') &
                              k, l, w1(k, l)*evolt*1000.0_dp
                        END IF

                     END DO
                     WRITE (unit_nr, *)
                  END DO
                  CALL cp_print_key_finished_output(unit_nr, logger, ec_sec, "PRINT%COUPLINGS")

               END IF

               IF (ASSOCIATED(w1)) &
                  DEALLOCATE (w1)
               IF (ASSOCIATED(w2)) &
                  DEALLOCATE (w2)

            END DO
         END DO

#else

         IF (output_unit > 0) THEN

            DO i = 1, ec%n_blocks
               DO j = i + 1, ec%n_blocks

                  WRITE (filename, '(a5,I1.1,a1,I1.1)') "ET_BL_", i, "-", j
                  unit_nr = cp_print_key_unit_nr(logger, ec_sec, "PRINT%COUPLINGS", extension=".elcoup", &
                                                 middle_name=TRIM(filename), file_position=my_pos, log_filename=.FALSE.)

                  WRITE (title, *) 'Coupling elements [meV] between blocks:', i, j

                  WRITE (unit_nr, *) TRIM(title)
                  IF (nspins > 1) THEN
                     WRITE (unit_nr, '(T3,A8,T13,A8,T28,A,A)') 'State A', 'State B', 'Coupling spin 1', '  Coupling spin 2'
                  ELSE
                     WRITE (unit_nr, '(T3,A8,T13,A8,T28,A)') 'State A', 'State B', 'Coupling'
                  END IF

                  DO k = 1, MIN(ec%block(i)%n_ao, n_states(1))
                     DO l = 1, MIN(ec%block(j)%n_ao, n_states(1))

                        IF (nspins > 1) THEN

                           WRITE (unit_nr, '(T3,I5,T13,I5,T22,E20.6)', ADVANCE='no') &
                              k, l, &
                              ec%block(i)%hab(1, j)%matrix%local_data(k, l)*evolt*1000.0_dp
                           IF ((k <= n_states(2)) .AND. (l <= n_states(2))) THEN
                              WRITE (unit_nr, '(E20.6)') &
                                 ec%block(i)%hab(2, j)%matrix%local_data(k, l)*evolt*1000.0_dp
                           ELSE
                              WRITE (unit_nr, *)
                           END IF

                        ELSE

                           WRITE (unit_nr, '(T3,I5,T13,I5,T22,E20.6)') &
                              k, l, &
                              ec%block(i)%hab(1, j)%matrix%local_data(k, l)*evolt*1000.0_dp

                        END IF

                     END DO
                     WRITE (unit_nr, *)
                  END DO

                  CALL cp_print_key_finished_output(unit_nr, logger, ec_sec, "PRINT%COUPLINGS")
               END DO
            END DO

         END IF

#endif

      END IF
   END SUBROUTINE print_couplings

! **************************************************************************************************
!> \brief Normalize set of MO vectors
!> \param qs_env QuickStep environment containing all system data
!> \param mo storage for the MO data set
!> \param n_ao number of AO basis functions
!> \param n_mo number of block states
!> \author Z. Futera (02.2017)
! **************************************************************************************************
   SUBROUTINE normalize_mo_vectors(qs_env, mo, n_ao, n_mo)

      ! Routine arguments
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mo_set_type), POINTER                         :: mo
      INTEGER                                            :: n_ao, n_mo

      REAL(KIND=dp), DIMENSION(:), POINTER               :: vec_t
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_s
      TYPE(cp_fm_type), POINTER                          :: mat_sc, mat_t
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_s

! Routine name for debug purposes

      ! Initialization
      NULLIFY (blacs_env)
      NULLIFY (para_env)
      NULLIFY (fm_s)
      NULLIFY (mat_s)
      NULLIFY (mat_sc)
      NULLIFY (mat_t)
      NULLIFY (vec_t)

      ! Overlap matrix
      CALL get_qs_env(qs_env, matrix_s=mat_s)

      ! Calculate S*C product
      CALL cp_fm_create(matrix=mat_sc, matrix_struct=mo%mo_coeff%matrix_struct, &
                        name='S*C PRODUCT MATRIX')
      CALL cp_dbcsr_sm_fm_multiply(mat_s(1)%matrix, mo%mo_coeff, mat_sc, n_mo)

      ! Calculate C^T*S*C
      CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env)
      CALL cp_fm_struct_create(fmstruct=fm_s, para_env=para_env, context=blacs_env, &
                               nrow_global=n_mo, ncol_global=n_mo)
      CALL cp_fm_create(matrix=mat_t, matrix_struct=fm_s, &
                        name='C^T*S*C OVERLAP PRODUCT MATRIX')
      CALL cp_gemm('T', 'N', n_mo, n_mo, n_ao, 1.0_dp, mo%mo_coeff, mat_sc, 0.0_dp, mat_t)

      ! Normalization
      ALLOCATE (vec_t(n_mo))
      CPASSERT(ASSOCIATED(vec_t))
      CALL cp_fm_vectorssum(mat_t, vec_t)
      vec_t = 1.0_dp/DSQRT(vec_t)
      CALL cp_fm_column_scale(mo%mo_coeff, vec_t)

      ! Clean memory
      CALL cp_fm_struct_release(fmstruct=fm_s)
      CALL cp_fm_release(matrix=mat_sc)
      CALL cp_fm_release(matrix=mat_t)
      IF (ASSOCIATED(vec_t)) &
         DEALLOCATE (vec_t)

   END SUBROUTINE normalize_mo_vectors

! **************************************************************************************************
!> \brief Transform block MO coefficients to original non-orthogonal basis set and save them
!> \param qs_env QuickStep environment containing all system data
!> \param ec electronic coupling data structure
!> \param id block ID
!> \param mo storage for the MO data set
!> \param mat_u matrix of the block states
!> \param n_ao number of AO basis functions
!> \param n_mo number of block states
!> \author Z. Futera (02.2017)
! **************************************************************************************************
   SUBROUTINE set_mo_coefficients(qs_env, ec, id, mo, mat_u, n_ao, n_mo)

      IMPLICIT NONE

      ! Routine arguments
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(et_cpl), POINTER                              :: ec
      INTEGER                                            :: id
      TYPE(mo_set_type), POINTER                         :: mo
      TYPE(cp_fm_type), POINTER                          :: mat_u
      INTEGER                                            :: n_ao
      INTEGER                                            :: n_mo

      ! Routine name for debug purposes

      ! Local variables
      INTEGER                                            :: ir, ic, jr, jc, nr, nc
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_type), POINTER                          :: mat_w
      TYPE(cp_para_env_type), POINTER                    :: para_env
#if defined(__SCALAPACK)
      INTEGER                                            :: is_p_row, is_p_col, iu_p_row, iu_p_col
      INTEGER                                            :: ns_p_rows, ns_p_cols, nu_p_rows, nu_p_cols
      INTEGER                                            :: cs_row, cs_col, cu_row, cu_col
      INTEGER                                            :: cs_row_f, cs_col_f, cu_row_f, cu_col_f
      REAL(KIND=dp)                                      :: xu

      ! number of processors
      ns_p_rows = mo%mo_coeff%matrix_struct%context%num_pe(1)
      ns_p_cols = mo%mo_coeff%matrix_struct%context%num_pe(2)
      nu_p_rows = mat_u%matrix_struct%context%num_pe(1)
      nu_p_cols = mat_u%matrix_struct%context%num_pe(2)

      ! position of processors
      cs_row = mo%mo_coeff%matrix_struct%context%mepos(1)
      cs_col = mo%mo_coeff%matrix_struct%context%mepos(2)
      cu_row = mat_u%matrix_struct%context%mepos(1)
      cu_col = mat_u%matrix_struct%context%mepos(2)
#endif

      NULLIFY (blacs_env)
      NULLIFY (para_env)
      NULLIFY (mat_w)

      ! Working matrix
      CALL cp_fm_create(matrix=mat_w, matrix_struct=mo%mo_coeff%matrix_struct, &
                        name='BLOCK MO-TRANSFORMATION WORKING MATRIX')
      CALL cp_fm_set_all(mat_w, 0.0_dp)

      ! Matrix-element reordering
      nr = 1
      ! Rows
      DO ir = 1, ec%block(id)%n_atoms
         DO jr = 1, ec%block(id)%atom(ir)%n_ao
            ! Columns
            nc = 1
            DO ic = 1, ec%block(id)%n_atoms
               DO jc = 1, ec%block(id)%atom(ic)%n_ao
#if defined(__SCALAPACK)
                  CALL infog2l(ec%block(id)%atom(ir)%ao_pos + jr - 1, nc, &
                               mat_w%matrix_struct%descriptor, &
                               ns_p_rows, ns_p_cols, cs_row, cs_col, &
                               is_p_row, is_p_col, cs_row_f, cs_col_f)
                  CALL infog2l(nr, nc, mat_u%matrix_struct%descriptor, &
                               nu_p_rows, nu_p_cols, cu_row, cu_col, &
                               iu_p_row, iu_p_col, cu_row_f, cu_col_f)
                  ! Local U element
                  IF ((cu_row_f == cu_row) .AND. (cu_col_f == cu_col)) THEN
                     xu = mat_u%local_data(iu_p_row, iu_p_col)
                     CALL dgebs2d(mat_u%matrix_struct%context%group, 'All', ' ', 1, 1, xu, 1)
                     IF ((cs_row_f == cs_row) .AND. (cs_col_f == cs_col)) THEN
                        mat_w%local_data(is_p_row, is_p_col) = xu
                     END IF
                     ! Remote U element
                  ELSE
                     CALL dgebr2d(mat_u%matrix_struct%context%group, &
                                  'All', ' ', 1, 1, xu, 1, cu_row_f, cu_col_f)
                     IF ((cs_row_f == cs_row) .AND. (cs_col_f == cs_col)) THEN
                        mat_w%local_data(is_p_row, is_p_col) = xu
                     END IF
                  END IF
#else
                  mat_w%local_data(ec%block(id)%atom(ir)%ao_pos + jr - 1, nc) = &
                     mat_u%local_data(nr, nc)
#endif
                  nc = nc + 1
               END DO
            END DO
            nr = nr + 1
         END DO
      END DO

      ! Transformation to original non-orthogonal basis set
      CALL cp_gemm("N", "N", n_ao, n_mo, n_ao, 1.0_dp, ec%m_transf, mat_w, 0.0_dp, mo%mo_coeff)
      CALL normalize_mo_vectors(qs_env, mo, n_ao, n_mo)

      ! Clean memory
      CALL cp_fm_release(matrix=mat_w)

   END SUBROUTINE set_mo_coefficients

! **************************************************************************************************
!> \brief Creates MO set corresponding to one atomic data block
!> \param qs_env QuickStep environment containing all system data
!> \param ec electronic coupling data structure
!> \param id block ID
!> \param spin spin component
!> \param mat_u matrix of the block states
!> \param vec_e array of the block eigenvalues
!> \author Z. Futera (02.2017)
! **************************************************************************************************
   SUBROUTINE create_block_mo_set(qs_env, ec, id, spin, mat_u, vec_e)

      ! Routine arguments
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(et_cpl), POINTER                              :: ec
      INTEGER                                            :: id, spin
      TYPE(cp_fm_type), POINTER                          :: mat_u
      REAL(KIND=dp), DIMENSION(:), POINTER               :: vec_e

      INTEGER                                            :: n_ao, n_el, n_mo
      REAL(KIND=dp)                                      :: mx_occ
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_s
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_cntrl
      TYPE(mo_set_type), POINTER                         :: mo
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(scf_control_type), POINTER                    :: scf_cntrl

! Routine name for debug purposes

      NULLIFY (blacs_env)
      NULLIFY (dft_cntrl)
      NULLIFY (para_env)
      NULLIFY (qs_kind_set)
      NULLIFY (fm_s)
      NULLIFY (scf_cntrl)
      NULLIFY (mo)

      ! Number of basis functions
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set)
      CALL get_qs_kind_set(qs_kind_set, nsgf=n_ao)

      ! Number of states
      n_mo = mat_u%matrix_struct%nrow_global
      IF (n_mo /= mat_u%matrix_struct%ncol_global) &
         CPABORT('block state matrix is not square')
      IF (n_mo /= SIZE(vec_e)) &
         CPABORT('inconsistent number of states / energies')

      ! Maximal occupancy
      CALL get_qs_env(qs_env, dft_control=dft_cntrl)
      mx_occ = 2.0_dp
      IF (dft_cntrl%nspins > 1) &
         mx_occ = 1.0_dp

      ! Number of electrons
      n_el = ec%block(id)%n_electrons
      IF (dft_cntrl%nspins > 1) THEN
         n_el = n_el/2
         IF (MOD(ec%block(id)%n_electrons, 2) == 1) THEN
            IF (spin == 1) &
               n_el = n_el + 1
         END IF
      END IF

      ! Memory allocation
      NULLIFY (ec%block(id)%mo(spin)%mo_set)
      CALL allocate_mo_set(ec%block(id)%mo(spin)%mo_set, n_ao, n_mo, n_el, REAL(n_el, dp), mx_occ, 0.0_dp)
      mo => ec%block(id)%mo(spin)%mo_set

      ! State energies
      ALLOCATE (mo%eigenvalues(n_mo))
      CPASSERT(ASSOCIATED(mo%eigenvalues))
      mo%eigenvalues = vec_e

      ! States coefficients
      CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env)
      CALL cp_fm_struct_create(fmstruct=fm_s, para_env=para_env, context=blacs_env, &
                               nrow_global=n_ao, ncol_global=n_mo)
      CALL cp_fm_create(matrix=mo%mo_coeff, matrix_struct=fm_s, name='BLOCK STATES')

      ! Transform MO coefficients to original non-orthogonal basis set
      CALL set_mo_coefficients(qs_env, ec, id, mo, mat_u, n_ao, n_mo)

      ! Occupancies
      ALLOCATE (mo%occupation_numbers(n_mo))
      CPASSERT(ASSOCIATED(mo%occupation_numbers))
      mo%occupation_numbers = 0.0_dp

      IF (n_el > 0) THEN
         CALL get_qs_env(qs_env, scf_control=scf_cntrl)
         CALL set_mo_occupation(mo_set=mo, smear=scf_cntrl%smear)
      END IF

      ! Clean memory
      CALL cp_fm_struct_release(fmstruct=fm_s)

   END SUBROUTINE create_block_mo_set

! **************************************************************************************************
!> \brief save given electronic state to cube files
!> \param qs_env QuickStep environment containing all system data
!> \param logger output logger
!> \param input input-file block print setting section
!> \param mo electronic states data
!> \param ib block ID
!> \param im state ID
!> \param is spin ID
!> \author Z. Futera (02.2017)
! **************************************************************************************************
   SUBROUTINE save_mo_cube(qs_env, logger, input, mo, ib, im, is)

      ! Routine arguments
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: input
      TYPE(mo_set_type), POINTER                         :: mo
      INTEGER                                            :: ib, im, is

      CHARACTER(LEN=default_path_length)                 :: filename
      CHARACTER(LEN=default_string_length)               :: title
      INTEGER                                            :: unit_nr
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: wf_g, wf_r
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_subsys_type), POINTER                      :: subsys

! Routine name for debug purposes

      ! Initialization
      NULLIFY (particles)
      NULLIFY (subsys)

      NULLIFY (pw_env)
      NULLIFY (pw_pools)
      NULLIFY (auxbas_pw_pool)

      NULLIFY (atomic_kind_set)
      NULLIFY (cell)
      NULLIFY (dft_control)
      NULLIFY (particle_set)
      NULLIFY (qs_kind_set)

      ! Name of the cube file
      WRITE (filename, '(A4,I1.1,A1,I5.5,A1,I1.1)') 'BWF_', ib, '_', im, '_', is
      ! Open the file
      unit_nr = cp_print_key_unit_nr(logger, input, 'MO_CUBES', extension='.cube', &
                                     middle_name=TRIM(filename), file_position='REWIND', log_filename=.FALSE.)
      ! Title of the file
      WRITE (title, *) 'WAVEFUNCTION ', im, ' block ', ib, ' spin ', is

      ! List of all atoms
      CALL get_qs_env(qs_env, subsys=subsys)
      CALL qs_subsys_get(subsys, particles=particles)

      ! Grids for wavefunction
      CALL get_qs_env(qs_env, pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, pw_pools=pw_pools)
      CALL pw_pool_create_pw(auxbas_pw_pool, wf_r%pw, &
                             use_data=REALDATA3D, in_space=REALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, wf_g%pw, &
                             use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)

      ! Calculate the grid values
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set, &
                      cell=cell, dft_control=dft_control, particle_set=particle_set)
      CALL calculate_wavefunction(mo%mo_coeff, im, wf_r, wf_g, atomic_kind_set, &
                                  qs_kind_set, cell, dft_control, particle_set, pw_env)
      CALL cp_pw_to_cube(wf_r%pw, unit_nr, title, particles=particles, &
                         stride=section_get_ivals(input, 'MO_CUBES%STRIDE'))

      ! Close file
      CALL cp_print_key_finished_output(unit_nr, logger, input, 'MO_CUBES')

      ! Clean memory
      CALL pw_pool_give_back_pw(auxbas_pw_pool, wf_r%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, wf_g%pw)

   END SUBROUTINE save_mo_cube

! **************************************************************************************************
!> \brief save specified electronic states to cube files
!> \param qs_env QuickStep environment containing all system data
!> \param ec electronic coupling data structure
!> \param n_spins number of spin states
!> \author Z. Futera (02.2017)
! **************************************************************************************************
   SUBROUTINE save_el_states(qs_env, ec, n_spins)

      ! Routine arguments
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(et_cpl), POINTER                              :: ec
      INTEGER                                            :: n_spins

      INTEGER                                            :: i, j, k, l, n
      INTEGER, DIMENSION(:), POINTER                     :: list
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(mo_set_type), POINTER                         :: mo
      TYPE(section_vals_type), POINTER                   :: block_sec, mo_sec, print_sec

! Routine name for debug purposes

      NULLIFY (logger)
      NULLIFY (block_sec)
      NULLIFY (print_sec)
      NULLIFY (mo_sec)

      ! Output logger
      logger => cp_get_default_logger()
      block_sec => section_vals_get_subs_vals(qs_env%input, &
                                              'PROPERTIES%ET_COUPLING%PROJECTION%BLOCK')

      ! Print states of all blocks
      DO i = 1, ec%n_blocks

         print_sec => section_vals_get_subs_vals(block_sec, 'PRINT', i_rep_section=i)

         ! Check if the print input section is active
         IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                              print_sec, 'MO_CUBES'), cp_p_file)) THEN

            mo_sec => section_vals_get_subs_vals(print_sec, 'MO_CUBES')

            ! Spin states
            DO j = 1, n_spins

               mo => ec%block(i)%mo(j)%mo_set

               CALL section_vals_val_get(mo_sec, keyword_name='MO_LIST', n_rep_val=n)

               ! List of specific MOs
               IF (n > 0) THEN

                  DO k = 1, n
                     NULLIFY (list)
                     CALL section_vals_val_get(mo_sec, keyword_name='MO_LIST', &
                                               i_rep_val=k, i_vals=list)
                     IF (ASSOCIATED(list)) THEN
                        DO l = 1, SIZE(list)
                           CALL save_mo_cube(qs_env, logger, print_sec, mo, i, list(l), j)
                        END DO
                     END IF
                  END DO

                  ! Frontier MOs
               ELSE

                  ! Occupied states
                  CALL section_vals_val_get(mo_sec, keyword_name='NHOMO', i_val=n)

                  IF (n > 0) THEN
                     DO k = MAX(1, mo%homo - n + 1), mo%homo
                        CALL save_mo_cube(qs_env, logger, print_sec, mo, i, k, j)
                     END DO
                  END IF

                  ! Unoccupied states
                  CALL section_vals_val_get(mo_sec, keyword_name='NLUMO', i_val=n)

                  IF (n > 0) THEN
                     DO k = mo%lfomo, MIN(mo%lfomo + n - 1, mo%nmo)
                        CALL save_mo_cube(qs_env, logger, print_sec, mo, i, k, j)
                     END DO
                  END IF

               END IF

            END DO

         END IF

      END DO

   END SUBROUTINE save_el_states

! **************************************************************************************************
!> \brief calculates the electron transfer coupling elements by projection-operator approach
!>        Kondov et al. J.Phys.Chem.C 2007, 111, 11970-11981
!> \param qs_env QuickStep environment containing all system data
!> \author Z. Futera (02.2017)
! **************************************************************************************************
   SUBROUTINE calc_et_coupling_proj(qs_env)

      ! Routine arguments
      TYPE(qs_environment_type), POINTER                 :: qs_env

      INTEGER                                            :: i, j, k, n_ao, n_atoms, output_unit
      LOGICAL                                            :: do_kp, master
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: mat_h
      TYPE(cp_fm_struct_type), POINTER                   :: fm_s
      TYPE(cp_fm_type), POINTER                          :: mat_w
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: ks, mo_der
      TYPE(dft_control_type), POINTER                    :: dft_cntrl
      TYPE(et_cpl), POINTER                              :: ec
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mo
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), POINTER                   :: et_proj_sec

! Routine name for debug purposes

      ! Pointer initialization
      NULLIFY (logger)

      NULLIFY (blacs_env)
      NULLIFY (para_env)
      NULLIFY (dft_cntrl)
      NULLIFY (kpoints)
      NULLIFY (qs_kind_set)
      NULLIFY (et_proj_sec)

      NULLIFY (fm_s)
      NULLIFY (mat_h)
      NULLIFY (mat_w)
      NULLIFY (ks, mo_der)

      NULLIFY (ec)

      ! Reference
      CALL cite_reference(Futera2017)

      ! Stream for output to LOG file
      logger => cp_get_default_logger()

      et_proj_sec => section_vals_get_subs_vals(qs_env%input, 'PROPERTIES%ET_COUPLING%PROJECTION')

      output_unit = cp_print_key_unit_nr(logger, et_proj_sec, &
                                         'PROGRAM_RUN_INFO', extension='.log')

      ! Parallel calculation - master thread
      master = .FALSE.
      IF (output_unit > 0) &
         master = .TRUE.

      ! Header
      IF (master) THEN
         WRITE (output_unit, '(/,T2,A)') &
            '!-----------------------------------------------------------------------------!'
         WRITE (output_unit, '(T17,A)') &
            'Electronic coupling - Projection-operator method'
      END IF

      ! Main data structure
      ALLOCATE (ec)
      CPASSERT(ASSOCIATED(ec))
      CALL set_block_data(qs_env, et_proj_sec, ec)

      ! Number of atoms and AO functions
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, natom=n_atoms)
      CALL get_qs_kind_set(qs_kind_set, nsgf=n_ao)

      ! Print out info about system partitioning
      IF (master) THEN

         WRITE (output_unit, '(/,T3,A,I10)') &
            'Number of atoms                    = ', n_atoms
         WRITE (output_unit, '(T3,A,I10)') &
            'Number of fragments                = ', ec%n_blocks
         WRITE (output_unit, '(T3,A,I10)') &
            'Number of fragment atoms           = ', ec%n_atoms
         WRITE (output_unit, '(T3,A,I10)') &
            'Number of unassigned atoms         = ', n_atoms - ec%n_atoms
         WRITE (output_unit, '(T3,A,I10)') &
            'Number of AO basis functions       = ', n_ao

         DO i = 1, ec%n_blocks

            WRITE (output_unit, '(/,T3,A,I0,A)') &
               'Block ', i, ':'
            WRITE (output_unit, '(T3,A,I10)') &
               'Number of block atoms              = ', ec%block(i)%n_atoms
            WRITE (output_unit, '(T3,A,I10)') &
               'Number of block electrons          = ', ec%block(i)%n_electrons
            WRITE (output_unit, '(T3,A,I10)') &
               'Number of block AO functions       = ', ec%block(i)%n_ao

            IF (ec%block(i)%n_atoms < 10) THEN

               WRITE (output_unit, '(T3,A,10I6)') &
                  'Block atom IDs                     =     ', &
                  (ec%block(i)%atom(j)%id, j=1, ec%block(i)%n_atoms)

            ELSE

               WRITE (output_unit, '(T3,A)') 'Block atom IDs                     ='
               DO j = 1, ec%block(i)%n_atoms/10
                  WRITE (output_unit, '(T3,A,10I6)') '      ', &
                     (ec%block(i)%atom((j - 1)*10 + k)%id, k=1, 10)
               END DO
               IF (MOD(ec%block(i)%n_atoms, 10) /= 0) THEN
                  WRITE (output_unit, '(T3,A,10I6)') '      ', &
                     (ec%block(i)%atom(k + 10*(ec%block(i)%n_atoms/10))%id, &
                      k=1, MOD(ec%block(i)%n_atoms, 10))
               END IF

            END IF

         END DO

      END IF

      ! Full matrix data structure
      CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env)
      CALL cp_fm_struct_create(fmstruct=fm_s, para_env=para_env, context=blacs_env, &
                               nrow_global=n_ao, ncol_global=n_ao)
      CALL cp_fm_create(matrix=mat_w, matrix_struct=fm_s, name='FULL WORK MATRIX')

      ! Spin polarization / K-point sampling
      CALL get_qs_env(qs_env, dft_control=dft_cntrl, do_kpoints=do_kp)
      CALL get_qs_env(qs_env, mos=mo, matrix_ks=ks, mo_derivs=mo_der, scf_control=scf_control)
      CALL make_mo_eig(mo, dft_cntrl%nspins, ks, scf_control, mo_der)

      IF (do_kp) THEN
         CPABORT('ET_COUPLING not implemented with kpoints')
      ELSE
         !  no K-points
         IF (master) &
            WRITE (output_unit, '(T3,A)') 'No K-point sampling (Gamma point only)'
      END IF

      IF (dft_cntrl%nspins == 2) THEN

         IF (master) &
            WRITE (output_unit, '(/,T3,A)') 'Spin-polarized calculation'

         !<--- Open shell / No K-points ------------------------------------------------>!

         ! State eneries of the whole system
         IF (mo(1)%mo_set%nao /= mo(2)%mo_set%nao) &
            CPABORT('different number of alpha/beta AO basis functions')
         IF (master) THEN
            WRITE (output_unit, '(/,T3,A,I10)') &
               'Number of AO basis functions       = ', mo(1)%mo_set%nao
            WRITE (output_unit, '(T3,A,I10)') &
               'Number of alpha states             = ', mo(1)%mo_set%nmo
            WRITE (output_unit, '(T3,A,I10)') &
               'Number of beta states              = ', mo(2)%mo_set%nmo
         END IF
         CALL print_states(output_unit, mo, dft_cntrl%nspins, 'the whole system', fermi=.TRUE.)
         CALL set_fermi(ec, mo(1)%mo_set%mu, mo(2)%mo_set%mu)

         ! KS Hamiltonian
         CALL get_block_hamiltonian(qs_env, ec, fm_s, mat_h, mat_w, n_ao, dft_cntrl%nspins)

         ! Block diagonization
         CALL hamiltonian_block_diag(qs_env, ec, mat_h)

         ! Print out energies and couplings
         DO i = 1, ec%n_blocks
            IF (output_unit > 0) THEN
               CALL print_states(output_unit, ec%block(i)%mo, dft_cntrl%nspins, &
                                 'block '//TRIM(ADJUSTL(cp_to_string(i)))//' states', &
                                 mx_mo_a=mo(1)%mo_set%nmo, mx_mo_b=mo(2)%mo_set%nmo, fermi=.TRUE.)
            END IF
            CALL print_mo_coeff(output_unit, qs_env, ec, i, dft_cntrl%nspins)
         END DO

         CALL print_couplings(et_proj_sec, output_unit, logger, ec, mo)

      ELSE

         IF (master) &
            WRITE (output_unit, '(/,T3,A)') 'Spin-restricted calculation'

         !<--- Close shell / No K-points ----------------------------------------------->!

         ! State eneries of the whole system
         IF (master) THEN
            WRITE (output_unit, '(/,T3,A,I10)') &
               'Number of AO basis functions       = ', mo(1)%mo_set%nao
            WRITE (output_unit, '(T3,A,I10)') &
               'Number of states                   = ', mo(1)%mo_set%nmo
         END IF
         CALL print_states(output_unit, mo, dft_cntrl%nspins, 'the whole system', fermi=.TRUE.)
         CALL set_fermi(ec, mo(1)%mo_set%mu)

         ! KS Hamiltonian
         CALL get_block_hamiltonian(qs_env, ec, fm_s, mat_h, mat_w, n_ao, dft_cntrl%nspins)

         ! Block diagonization
         CALL hamiltonian_block_diag(qs_env, ec, mat_h)

         ! Print out energies and couplings
         DO i = 1, ec%n_blocks
            IF (output_unit > 0) THEN
               CALL print_states(output_unit, ec%block(i)%mo, dft_cntrl%nspins, &
                                 'block '//TRIM(ADJUSTL(cp_to_string(i)))//' states', &
                                 mx_mo_a=mo(1)%mo_set%nmo, fermi=.TRUE.)
            END IF
            CALL print_mo_coeff(output_unit, qs_env, ec, i, dft_cntrl%nspins)
         END DO

         CALL print_couplings(et_proj_sec, output_unit, logger, ec, mo)

      END IF

      ! Save electronic states
      CALL save_el_states(qs_env, ec, dft_cntrl%nspins)

      ! Footer
      IF (master) WRITE (output_unit, '(/,T2,A)') &
         '!-----------------------------------------------------------------------------!'

      ! Clean memory
      CALL cp_fm_struct_release(fmstruct=fm_s)
      CALL cp_fm_release(matrix=mat_w)
      IF (ASSOCIATED(mat_h)) THEN
         DO i = 1, SIZE(mat_h)
            IF (ASSOCIATED(mat_h(i)%matrix)) &
               CALL cp_fm_release(matrix=mat_h(i)%matrix)
         END DO
         DEALLOCATE (mat_h)
      END IF
      CALL release_ec_data(ec)

      ! Close output stream
      CALL cp_print_key_finished_output(output_unit, logger, et_proj_sec, 'PROGRAM_RUN_INFO')

   END SUBROUTINE calc_et_coupling_proj

END MODULE et_coupling_proj
