!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Determine active space Hamiltonian
!> \par History
!>      04.2016 created [JGH]
!> \author JGH
! **************************************************************************************************
MODULE qs_active_space_methods
   USE admm_types, ONLY: admm_type, &
                         get_admm_env, &
                         admm_env_release
   USE atomic_kind_types, ONLY: atomic_kind_type
   USE basis_set_types, ONLY: allocate_sto_basis_set, &
                              create_gto_from_sto_basis, &
                              deallocate_sto_basis_set, &
                              gto_basis_set_type, &
                              init_orb_basis_set, &
                              set_sto_basis_set, &
                              srules, &
                              sto_basis_set_type
   USE cell_types, ONLY: cell_type, use_perd_none, use_perd_xyz
   USE cell_methods, ONLY: init_cell, set_cell_param, write_cell_low
   USE cp_blacs_env, ONLY: cp_blacs_env_type, cp_blacs_env_create, cp_blacs_env_release, BLACS_GRID_SQUARE
   USE cp_control_types, ONLY: dft_control_type, qs_control_type
   USE cp_dbcsr_operations, ONLY: cp_dbcsr_plus_fm_fm_t, &
                                  cp_dbcsr_sm_fm_multiply, &
                                  dbcsr_allocate_matrix_set, &
                                  cp_dbcsr_m_by_n_from_template, copy_dbcsr_to_fm
   USE cp_dbcsr_output, ONLY: cp_dbcsr_write_sparse_matrix
   USE cp_files, ONLY: close_file, &
                       file_exists, &
                       open_file
   USE cp_fm_basic_linalg, ONLY: cp_fm_column_scale
   USE cp_fm_struct, ONLY: cp_fm_struct_create, &
                           cp_fm_struct_release, &
                           cp_fm_struct_type
   USE cp_fm_types, ONLY: &
      cp_fm_create, cp_fm_get_element, cp_fm_get_info, cp_fm_init_random, cp_fm_release, &
      cp_fm_set_all, cp_fm_set_element, cp_fm_to_fm, cp_fm_type, cp_fm_write_formatted
   USE cp_log_handling, ONLY: cp_get_default_logger, &
                              cp_logger_get_default_io_unit, &
                              cp_logger_type
   USE cp_output_handling, ONLY: &
      cp_p_file, cp_print_key_finished_output, cp_print_key_should_output, cp_print_key_unit_nr, &
      debug_print_level, high_print_level, low_print_level, medium_print_level, &
      silent_print_level
   USE cp_realspace_grid_cube, ONLY: cp_pw_to_cube
   USE cp_dbcsr_api, ONLY: &
      dbcsr_copy, dbcsr_csr_create, dbcsr_csr_type, dbcsr_p_type, dbcsr_type, dbcsr_release, &
      dbcsr_type_no_symmetry, dbcsr_create, dbcsr_set, dbcsr_multiply, dbcsr_iterator_next_block, &
      dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_blocks_left, &
      dbcsr_iterator_type, dbcsr_type_symmetric, dbcsr_get_occupation, dbcsr_get_info
   USE erf_complex, ONLY: erfz_fast
   USE group_dist_types, ONLY: get_group_dist, release_group_dist, group_dist_d1_type
   USE input_constants, ONLY: &
      casci_canonical, eri_method_full_gpw, eri_method_gpw_ht, eri_operator_coulomb, &
      eri_operator_erf, eri_operator_erfc, eri_operator_gaussian, eri_operator_yukawa, &
      eri_operator_trunc, eri_operator_lr_trunc, &
      manual_selection, mao_projection, no_solver, qiskit_solver, wannier_projection, &
      eri_poisson_analytic, eri_poisson_periodic, eri_poisson_mt, high_spin_roks
   USE input_section_types, ONLY: section_vals_get, section_vals_get_subs_vals, &
                                  section_vals_set_subs_vals, section_vals_type, &
                                  section_vals_val_get, &
                                  section_vals_val_set
   USE ISO_C_BINDING, ONLY: c_null_char
   USE kinds, ONLY: default_path_length, &
                    default_string_length, &
                    dp, &
                    int_8
   USE hfx_types, ONLY: hfx_create, hfx_release
   USE machine, ONLY: m_walltime, m_flush
   USE mathlib, ONLY: diamat_all
   USE mathconstants, ONLY: fourpi, twopi, pi, rootpi
   USE memory_utilities, ONLY: reallocate
   USE message_passing, ONLY: mp_comm_type, &
                              mp_para_env_type, &
                              mp_para_env_release
   USE mp2_gpw, ONLY: create_mat_munu, grep_rows_in_subgroups, build_dbcsr_from_rows
   USE mt_util, ONLY: MT0D
   USE parallel_gemm_api, ONLY: parallel_gemm
   USE particle_list_types, ONLY: particle_list_type
   USE particle_types, ONLY: particle_type
   USE periodic_table, ONLY: ptable
   USE physcon, ONLY: angstrom, bohr
   USE preconditioner_types, ONLY: preconditioner_type
   USE pw_env_methods, ONLY: pw_env_create, &
                             pw_env_rebuild
   USE pw_env_types, ONLY: pw_env_get, &
                           pw_env_release, &
                           pw_env_type
   USE pw_methods, ONLY: pw_integrate_function, &
                         pw_multiply, &
                         pw_multiply_with, &
                         pw_transfer, &
                         pw_zero, pw_integral_ab, pw_scale, &
                         pw_gauss_damp, pw_compl_gauss_damp
   USE pw_poisson_methods, ONLY: pw_poisson_rebuild, &
                                 pw_poisson_solve
   USE pw_poisson_types, ONLY: ANALYTIC0D, &
                               PERIODIC3D, &
                               greens_fn_type, &
                               pw_poisson_analytic, &
                               pw_poisson_periodic, &
                               pw_poisson_type
   USE pw_pool_types, ONLY: &
      pw_pool_type
   USE pw_types, ONLY: &
      pw_c1d_gs_type, &
      pw_r3d_rs_type
   USE qcschema, ONLY: qcschema_env_create, &
                       qcschema_env_release, &
                       qcschema_to_hdf5, &
                       qcschema_type
   USE qs_active_space_types, ONLY: active_space_type, &
                                    create_active_space_type, &
                                    csr_idx_from_combined, &
                                    csr_idx_to_combined, &
                                    eri_type, &
                                    eri_type_eri_element_func, &
                                    get_irange_csr
   USE qs_active_space_utils, ONLY: eri_to_array, &
                                    subspace_matrix_to_array
   USE qs_collocate_density, ONLY: calculate_wavefunction
   USE qs_density_matrices, ONLY: calculate_density_matrix
   USE qs_energy_types, ONLY: qs_energy_type
   USE qs_environment_types, ONLY: get_qs_env, &
                                   qs_environment_type, &
                                   set_qs_env
   USE qs_integrate_potential, ONLY: integrate_v_rspace
   USE qs_kind_types, ONLY: qs_kind_type
   USE qs_ks_methods, ONLY: qs_ks_update_qs_env, qs_ks_build_kohn_sham_matrix, &
                            evaluate_core_matrix_traces
   USE qs_ks_types, ONLY: qs_ks_did_change, &
                          qs_ks_env_type, set_ks_env
   USE qs_mo_io, ONLY: write_mo_set_to_output_unit
   USE qs_mo_methods, ONLY: calculate_subspace_eigenvalues
   USE qs_mo_types, ONLY: allocate_mo_set, &
                          get_mo_set, &
                          init_mo_set, &
                          mo_set_type
   USE qs_neighbor_list_types, ONLY: neighbor_list_set_p_type, release_neighbor_list_sets
   USE qs_ot_eigensolver, ONLY: ot_eigensolver
   USE qs_rho_methods, ONLY: qs_rho_update_rho
   USE qs_rho_types, ONLY: qs_rho_get, &
                           qs_rho_type
   USE qs_subsys_types, ONLY: qs_subsys_get, &
                              qs_subsys_type
   USE qs_scf_post_scf, ONLY: qs_scf_compute_properties
   USE scf_control_types, ONLY: scf_control_type
#ifndef __NO_SOCKETS
   USE sockets_interface, ONLY: accept_socket, &
                                close_socket, &
                                listen_socket, &
                                open_bind_socket, &
                                readbuffer, &
                                remove_socket_file, &
                                writebuffer
#endif
   USE task_list_methods, ONLY: generate_qs_task_list
   USE task_list_types, ONLY: allocate_task_list, &
                              deallocate_task_list, &
                              task_list_type
   USE util, ONLY: get_limit
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_active_space_methods'

   PUBLIC :: active_space_main

   TYPE, EXTENDS(eri_type_eri_element_func) :: eri_fcidump_print
      INTEGER :: unit_nr = -1, bra_start = -1, ket_start = -1
   CONTAINS
      PROCEDURE :: func => eri_fcidump_print_func
   END TYPE eri_fcidump_print

   TYPE, EXTENDS(eri_type_eri_element_func) :: eri_fcidump_checksum
      INTEGER :: bra_start = 0, ket_start = 0
      REAL(KIND=dp) :: checksum = 0.0_dp
   CONTAINS
      PROCEDURE, PASS :: set => eri_fcidump_set
      PROCEDURE :: func => eri_fcidump_checksum_func
   END TYPE eri_fcidump_checksum

CONTAINS

! **************************************************************************************************
!> \brief Sets the starting indices of the bra and ket.
!> \param this object reference
!> \param bra_start starting index of the bra
!> \param ket_start starting index of the ket
! **************************************************************************************************
   SUBROUTINE eri_fcidump_set(this, bra_start, ket_start)
      CLASS(eri_fcidump_checksum) :: this
      INTEGER, INTENT(IN) :: bra_start, ket_start
      this%bra_start = bra_start
      this%ket_start = ket_start
   END SUBROUTINE eri_fcidump_set

! **************************************************************************************************
!> \brief Main method for determining the active space Hamiltonian
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE active_space_main(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'active_space_main'

      CHARACTER(len=10)                                  :: cshell, lnam(5)
      CHARACTER(len=default_path_length)                 :: qcschema_filename
      CHARACTER(LEN=default_string_length)               :: basis_type
      INTEGER :: as_solver, eri_method, eri_operator, eri_print, group_size, handle, i, iatom, &
         ishell, isp, ispin, iw, j, jm, m, max_orb_ind, mselect, n1, n2, nao, natom, nel, &
         nelec_active, nelec_inactive, nelec_total, nmo, nmo_active, nmo_available, nmo_inactive, &
         nmo_inactive_remaining, nmo_occ, nmo_virtual, nn1, nn2, nrow_global, nspins
      INTEGER, DIMENSION(5)                              :: nshell
      INTEGER, DIMENSION(:), POINTER                     :: invals
      LOGICAL                                            :: do_ddapc, do_kpoints, ex_omega, &
                                                            ex_operator, ex_perd, ex_rcut, &
                                                            explicit, stop_after_print, store_wfn
      REAL(KIND=dp) :: alpha, eri_eps_filter, eri_eps_grid, eri_eps_int, eri_gpw_cutoff, &
         eri_op_omega, eri_rcut, eri_rel_cutoff, fel, focc, maxocc, nze_percentage
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: eigenvalues
      REAL(KIND=dp), DIMENSION(:), POINTER               :: evals_virtual
      TYPE(active_space_type), POINTER                   :: active_space_env
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: context
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp
      TYPE(cp_fm_type)                                   :: fm_dummy, mo_virtual
      TYPE(cp_fm_type), POINTER                          :: fm_target_active, fm_target_inactive, &
                                                            fmat, mo_coeff, mo_ref, mo_target
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_csr_type), POINTER                      :: eri_mat
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: ks_matrix, rho_ao, s_matrix
      TYPE(dbcsr_type), POINTER                          :: denmat
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mo_set_type), POINTER                         :: mo_set, mo_set_active, mo_set_inactive
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(preconditioner_type), POINTER                 :: local_preconditioner
      TYPE(qcschema_type)                                :: qcschema_env
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), POINTER                   :: adiabatic_rescaling, as_input, &
                                                            hfx_section, input, loc_print, &
                                                            loc_section, print_orb, xc_section

      !--------------------------------------------------------------------------------------------!

      CALL get_qs_env(qs_env, input=input)
      as_input => section_vals_get_subs_vals(input, "DFT%ACTIVE_SPACE")
      CALL section_vals_get(as_input, explicit=explicit)
      IF (.NOT. explicit) RETURN
      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      iw = cp_logger_get_default_io_unit(logger)

      IF (iw > 0) THEN
         WRITE (iw, '(/,T2,A)') &
            '!-----------------------------------------------------------------------------!'
         WRITE (iw, '(T26,A)') "Active Space Embedding Module"
         WRITE (iw, '(T2,A)') &
            '!-----------------------------------------------------------------------------!'
      END IF

      ! k-points?
      CALL get_qs_env(qs_env, do_kpoints=do_kpoints, dft_control=dft_control)
      IF (do_kpoints) THEN
         CALL cp_abort(__LOCATION__, "k-points not supported in active space module")
      END IF

      ! adiabatic rescaling?
      adiabatic_rescaling => section_vals_get_subs_vals(input, "DFT%XC%ADIABATIC_RESCALING")
      CALL section_vals_get(adiabatic_rescaling, explicit=explicit)
      IF (explicit) THEN
         CALL cp_abort(__LOCATION__, "Adiabatic rescaling not supported in active space module")
      END IF

      ! Setup the possible usage of DDAPC charges
      do_ddapc = dft_control%qs_control%ddapc_restraint .OR. &
                 qs_env%cp_ddapc_ewald%do_decoupling .OR. &
                 qs_env%cp_ddapc_ewald%do_qmmm_periodic_decpl .OR. &
                 qs_env%cp_ddapc_ewald%do_solvation
      IF (do_ddapc) THEN
         CALL cp_abort(__LOCATION__, "DDAPC charges are not supported in the active space module")
      END IF
      IF (dft_control%do_sccs) THEN
         CALL cp_abort(__LOCATION__, "SCCS is not supported in the active space module")
      END IF
      IF (dft_control%correct_surf_dip) THEN
         IF (dft_control%surf_dip_correct_switch) THEN
            CALL cp_abort(__LOCATION__, "Surface dipole correction not supported in the AS module")
         END IF
      END IF
      IF (dft_control%smeagol_control%smeagol_enabled) THEN
         CALL cp_abort(__LOCATION__, "SMEAGOL is not supported in the active space module")
      END IF
      IF (dft_control%qs_control%do_kg) THEN
         CALL cp_abort(__LOCATION__, "KG correction not supported in the active space module")
      END IF

      NULLIFY (active_space_env)
      CALL create_active_space_type(active_space_env)
      active_space_env%energy_total = 0.0_dp
      active_space_env%energy_ref = 0.0_dp
      active_space_env%energy_inactive = 0.0_dp
      active_space_env%energy_active = 0.0_dp

      ! input options

      ! figure out what needs to be printed/stored
      IF (BTEST(cp_print_key_should_output(logger%iter_info, as_input, "FCIDUMP"), cp_p_file)) THEN
         active_space_env%fcidump = .TRUE.
      END IF

      CALL section_vals_val_get(as_input, "QCSCHEMA", c_val=qcschema_filename, explicit=explicit)
      IF (explicit) THEN
         active_space_env%qcschema = .TRUE.
         active_space_env%qcschema_filename = qcschema_filename
      END IF

      CALL section_vals_val_get(as_input, "ACTIVE_ELECTRONS", i_val=nelec_active)
      CALL get_qs_env(qs_env, nelectron_total=nelec_total)

      IF (nelec_active <= 0) CPABORT("Specify a positive number of active electrons.")
      IF (nelec_active > nelec_total) CPABORT("More active electrons than total electrons.")

      nelec_inactive = nelec_total - nelec_active
      IF (MOD(nelec_inactive, 2) /= 0) THEN
         CPABORT("The remaining number of inactive electrons has to be even.")
      END IF

      CALL section_vals_val_get(as_input, "ALPHA", r_val=alpha)
      IF (alpha < 0.0 .OR. alpha > 1.0) CPABORT("Specify a damping factor between 0 and 1.")
      active_space_env%alpha = alpha

      IF (iw > 0) THEN
         WRITE (iw, '(T3,A,T70,I10)') "Total number of electrons", nelec_total
         WRITE (iw, '(T3,A,T70,I10)') "Number of inactive electrons", nelec_inactive
         WRITE (iw, '(T3,A,T70,I10)') "Number of active electrons", nelec_active
      END IF

      CALL get_qs_env(qs_env, dft_control=dft_control)
      nspins = dft_control%nspins

      active_space_env%nelec_active = nelec_active
      active_space_env%nelec_inactive = nelec_inactive
      active_space_env%nelec_total = nelec_total
      active_space_env%nspins = nspins
      active_space_env%multiplicity = dft_control%multiplicity

      ! define the active/inactive space orbitals
      CALL section_vals_val_get(as_input, "ACTIVE_ORBITALS", explicit=explicit, i_val=nmo_active)
      IF (.NOT. explicit) THEN
         CALL cp_abort(__LOCATION__, "Number of Active Orbitals has to be specified.")
      END IF
      active_space_env%nmo_active = nmo_active
      ! this is safe because nelec_inactive is always even
      nmo_inactive = nelec_inactive/2
      active_space_env%nmo_inactive = nmo_inactive

      CALL section_vals_val_get(as_input, "ORBITAL_SELECTION", i_val=mselect)
      IF (iw > 0) THEN
         SELECT CASE (mselect)
         CASE DEFAULT
            CPABORT("Unknown orbital selection method")
         CASE (casci_canonical)
            WRITE (iw, '(/,T3,A)') &
               "Active space orbitals selected using energy ordered canonical orbitals"
         CASE (wannier_projection)
            WRITE (iw, '(/,T3,A)') &
               "Active space orbitals selected using projected Wannier orbitals"
         CASE (mao_projection)
            WRITE (iw, '(/,T3,A)') &
               "Active space orbitals selected using modified atomic orbitals (MAO)"
         CASE (manual_selection)
            WRITE (iw, '(/,T3,A)') &
               "Active space orbitals selected manually"
         END SELECT

         WRITE (iw, '(T3,A,T70,I10)') "Number of inactive orbitals", nmo_inactive
         WRITE (iw, '(T3,A,T70,I10)') "Number of active orbitals", nmo_active
      END IF

      ! get projection spaces
      CALL section_vals_val_get(as_input, "SUBSPACE_ATOM", i_val=iatom, explicit=explicit)
      IF (explicit) THEN
         CALL get_qs_env(qs_env, natom=natom)
         IF (iatom <= 0 .OR. iatom > natom) THEN
            IF (iw > 0) THEN
               WRITE (iw, '(/,T3,A,I3)') "ERROR: SUBSPACE_ATOM number is not valid", iatom
            END IF
            CPABORT("Select a valid SUBSPACE_ATOM")
         END IF
      END IF
      CALL section_vals_val_get(as_input, "SUBSPACE_SHELL", c_val=cshell, explicit=explicit)
      nshell = 0
      lnam = ""
      IF (explicit) THEN
         cshell = ADJUSTL(cshell)
         n1 = 1
         DO i = 1, 5
            ishell = i
            IF (cshell(n1:n1) == " ") THEN
               ishell = ishell - 1
               EXIT
            END IF
            READ (cshell(n1:), "(I1,A1)") nshell(i), lnam(i)
            n1 = n1 + 2
         END DO
      END IF

      ! generate orbitals
      SELECT CASE (mselect)
      CASE DEFAULT
         CPABORT("Unknown orbital selection method")
      CASE (casci_canonical)
         CALL get_qs_env(qs_env, mos=mos)

         ! total number of occupied orbitals, i.e. inactive plus active MOs
         nmo_occ = nmo_inactive + nmo_active

         ! set inactive orbital indices, these are trivially 1...nmo_inactive
         ALLOCATE (active_space_env%inactive_orbitals(nmo_inactive, nspins))
         DO ispin = 1, nspins
            DO i = 1, nmo_inactive
               active_space_env%inactive_orbitals(i, ispin) = i
            END DO
         END DO

         ! set active orbital indices, these are shifted by nmo_inactive
         ALLOCATE (active_space_env%active_orbitals(nmo_active, nspins))
         DO ispin = 1, nspins
            DO i = 1, nmo_active
               active_space_env%active_orbitals(i, ispin) = nmo_inactive + i
            END DO
         END DO

         ! allocate and initialize inactive and active mo coefficients.
         ! These are stored in a data structure for the full occupied space:
         ! for inactive mos, the active subset is set to zero, vice versa for the active mos
         ! TODO: allocate data structures only for the eaxct number MOs
         maxocc = 2.0_dp
         IF (nspins > 1) maxocc = 1.0_dp
         ALLOCATE (active_space_env%mos_active(nspins))
         ALLOCATE (active_space_env%mos_inactive(nspins))
         DO ispin = 1, nspins
            CALL get_mo_set(mos(ispin), mo_coeff=mo_ref, nao=nao)
            CALL cp_fm_get_info(mo_ref, context=context, para_env=para_env, nrow_global=nrow_global)
            ! the right number of active electrons per spin channel is initialized further down
            CALL allocate_mo_set(active_space_env%mos_active(ispin), nao, nmo_occ, 0, 0.0_dp, maxocc, 0.0_dp)
            CALL cp_fm_struct_create(fm_struct_tmp, para_env=para_env, context=context, &
                                     nrow_global=nrow_global, ncol_global=nmo_occ)
            CALL init_mo_set(active_space_env%mos_active(ispin), fm_struct=fm_struct_tmp, name="Active Space MO")
            CALL cp_fm_struct_release(fm_struct_tmp)
            IF (nspins == 2) THEN
               nel = nelec_inactive/2
            ELSE
               nel = nelec_inactive
            END IF
            CALL allocate_mo_set(active_space_env%mos_inactive(ispin), nao, nmo_occ, nel, &
                                 REAL(nel, KIND=dp), maxocc, 0.0_dp)
            CALL cp_fm_struct_create(fm_struct_tmp, para_env=para_env, context=context, &
                                     nrow_global=nrow_global, ncol_global=nmo_occ)
            CALL init_mo_set(active_space_env%mos_inactive(ispin), fm_struct=fm_struct_tmp, name="Inactive Space MO")
            CALL cp_fm_struct_release(fm_struct_tmp)
         END DO

         ! create canonical orbitals
         CALL get_qs_env(qs_env, scf_control=scf_control)
         IF (dft_control%roks .AND. scf_control%roks_scheme /= high_spin_roks) THEN
            CPABORT("Unclear how we define MOs in the general restricted case ... stopping")
         ELSE
            IF (dft_control%do_admm) THEN
               IF (dft_control%do_admm_mo) THEN
                  CPABORT("ADMM currently possible only with purification none_dm")
               END IF
            END IF

            ALLOCATE (eigenvalues(nmo_occ, nspins))
            eigenvalues = 0.0_dp
            CALL get_qs_env(qs_env, matrix_ks=ks_matrix, matrix_s=s_matrix, scf_control=scf_control)

            ! calculate virtual MOs and copy inactive and active orbitals
            IF (iw > 0) THEN
               WRITE (iw, '(/,T3,A)') "Calculating virtual MOs..."
            END IF
            DO ispin = 1, nspins
               ! nmo_available is the number of MOs available from the SCF calculation:
               ! this is at least the number of occupied orbitals in the SCF, plus
               ! any number of added MOs (virtuals) requested in the SCF section
               CALL get_mo_set(mos(ispin), mo_coeff=mo_ref, nmo=nmo_available)

               ! calculate how many extra MOs we still have to compute
               nmo_virtual = nmo_occ - nmo_available
               nmo_virtual = MAX(nmo_virtual, 0)

               NULLIFY (evals_virtual)
               ALLOCATE (evals_virtual(nmo_virtual))

               CALL cp_fm_get_info(mo_ref, context=context, para_env=para_env, &
                                   nrow_global=nrow_global)

               CALL cp_fm_struct_create(fm_struct_tmp, para_env=para_env, context=context, &
                                        nrow_global=nrow_global, ncol_global=nmo_virtual)
               CALL cp_fm_create(mo_virtual, fm_struct_tmp, name="virtual")
               CALL cp_fm_struct_release(fm_struct_tmp)
               CALL cp_fm_init_random(mo_virtual, nmo_virtual)

               NULLIFY (local_preconditioner)

               ! compute missing virtual MOs
               CALL ot_eigensolver(matrix_h=ks_matrix(ispin)%matrix, matrix_s=s_matrix(1)%matrix, &
                                   matrix_c_fm=mo_virtual, matrix_orthogonal_space_fm=mo_ref, &
                                   eps_gradient=scf_control%eps_lumos, &
                                   preconditioner=local_preconditioner, &
                                   iter_max=scf_control%max_iter_lumos, &
                                   size_ortho_space=nmo_available)

               ! get the eigenvalues
               CALL calculate_subspace_eigenvalues(mo_virtual, ks_matrix(ispin)%matrix, evals_virtual)

               ! we need to send the copy of MOs to preserve the sign
               CALL cp_fm_create(fm_dummy, mo_ref%matrix_struct)
               CALL cp_fm_to_fm(mo_ref, fm_dummy)
               CALL calculate_subspace_eigenvalues(fm_dummy, ks_matrix(ispin)%matrix, &
                                                   evals_arg=eigenvalues(:, ispin), do_rotation=.TRUE.)

               ! copy inactive orbitals
               mo_set => active_space_env%mos_inactive(ispin)
               CALL get_mo_set(mo_set, mo_coeff=mo_target)
               DO i = 1, SIZE(active_space_env%inactive_orbitals, 1)
                  m = active_space_env%inactive_orbitals(i, ispin)
                  CALL cp_fm_to_fm(mo_ref, mo_target, 1, m, m)
                  mo_set%eigenvalues(m) = eigenvalues(m, ispin)
                  IF (nspins > 1) THEN
                     mo_set%occupation_numbers(m) = 1.0
                  ELSE
                     mo_set%occupation_numbers(m) = 2.0
                  END IF
               END DO

               ! copy active orbitals
               mo_set => active_space_env%mos_active(ispin)
               CALL get_mo_set(mo_set, mo_coeff=mo_target)
               ! for mult > 1, put the polarized electrons in the alpha channel
               IF (nspins == 2) THEN
                  IF (ispin == 1) THEN
                     nel = (nelec_active + active_space_env%multiplicity - 1)/2
                  ELSE
                     nel = (nelec_active - active_space_env%multiplicity + 1)/2
                  END IF
               ELSE
                  nel = nelec_active
               END IF
               mo_set%nelectron = nel
               mo_set%n_el_f = REAL(nel, KIND=dp)
               DO i = 1, nmo_active
                  m = active_space_env%active_orbitals(i, ispin)
                  IF (m > nmo_available) THEN
                     CALL cp_fm_to_fm(mo_virtual, mo_target, 1, m - nmo_available, m)
                     eigenvalues(m, ispin) = evals_virtual(m - nmo_available)
                     mo_set%occupation_numbers(m) = 0.0
                  ELSE
                     CALL cp_fm_to_fm(mo_ref, mo_target, 1, m, m)
                     mo_set%occupation_numbers(m) = mos(ispin)%occupation_numbers(m)
                  END IF
                  mo_set%eigenvalues(m) = eigenvalues(m, ispin)
               END DO
               ! Release
               DEALLOCATE (evals_virtual)
               CALL cp_fm_release(fm_dummy)
               CALL cp_fm_release(mo_virtual)
            END DO

            IF (iw > 0) THEN
               DO ispin = 1, nspins
                  WRITE (iw, '(/,T3,A,I3,T66,A)') "Canonical Orbital Selection for spin", ispin, &
                     "[atomic units]"
                  DO i = 1, nmo_inactive, 4
                     jm = MIN(3, nmo_inactive - i)
                     WRITE (iw, '(T3,4(F14.6,A5))') (eigenvalues(i + j, ispin), " [I]", j=0, jm)
                  END DO
                  DO i = nmo_inactive + 1, nmo_inactive + nmo_active, 4
                     jm = MIN(3, nmo_inactive + nmo_active - i)
                     WRITE (iw, '(T3,4(F14.6,A5))') (eigenvalues(i + j, ispin), " [A]", j=0, jm)
                  END DO
                  WRITE (iw, '(/,T3,A,I3)') "Active Orbital Indices for spin", ispin
                  DO i = 1, SIZE(active_space_env%active_orbitals, 1), 4
                     jm = MIN(3, SIZE(active_space_env%active_orbitals, 1) - i)
                     WRITE (iw, '(T3,4(I4))') (active_space_env%active_orbitals(i + j, ispin), j=0, jm)
                  END DO
               END DO
            END IF
            DEALLOCATE (eigenvalues)
         END IF

      CASE (manual_selection)
         ! create canonical orbitals
         IF (dft_control%roks) THEN
            CPABORT("Unclear how we define MOs in the restricted case ... stopping")
         ELSE
            IF (dft_control%do_admm) THEN
               ! For admm_mo, the auxiliary density is computed from the MOs, which never change
               ! in the rs-dft embedding, therefore the energy is wrong as the LR HFX never changes.
               ! For admm_dm, the auxiliary density is computed from the density matrix, which is
               ! updated at each iteration and therefore works.
               IF (dft_control%do_admm_mo) THEN
                  CPABORT("ADMM currently possible only with purification none_dm")
               END IF
            END IF

            CALL section_vals_val_get(as_input, "ACTIVE_ORBITAL_INDICES", explicit=explicit, i_vals=invals)
            IF (.NOT. explicit) THEN
               CALL cp_abort(__LOCATION__, "Manual orbital selection requires to explicitly "// &
                             "set the active orbital indices via ACTIVE_ORBITAL_INDICES")
            END IF

            IF (nspins == 1) THEN
               CPASSERT(SIZE(invals) == nmo_active)
            ELSE
               CPASSERT(SIZE(invals) == 2*nmo_active)
            END IF
            ALLOCATE (active_space_env%inactive_orbitals(nmo_inactive, nspins))
            ALLOCATE (active_space_env%active_orbitals(nmo_active, nspins))

            DO ispin = 1, nspins
               DO i = 1, nmo_active
                  active_space_env%active_orbitals(i, ispin) = invals(i + (ispin - 1)*nmo_active)
               END DO
            END DO

            CALL get_qs_env(qs_env, mos=mos)

            ! include MOs up to the largest index in the list
            max_orb_ind = MAXVAL(invals)
            maxocc = 2.0_dp
            IF (nspins > 1) maxocc = 1.0_dp
            ALLOCATE (active_space_env%mos_active(nspins))
            ALLOCATE (active_space_env%mos_inactive(nspins))
            DO ispin = 1, nspins
               ! init active orbitals
               CALL get_mo_set(mos(ispin), mo_coeff=mo_ref, nao=nao)
               CALL cp_fm_get_info(mo_ref, context=context, para_env=para_env, nrow_global=nrow_global)
               CALL allocate_mo_set(active_space_env%mos_active(ispin), nao, max_orb_ind, 0, 0.0_dp, maxocc, 0.0_dp)
               CALL cp_fm_struct_create(fm_struct_tmp, para_env=para_env, context=context, &
                                        nrow_global=nrow_global, ncol_global=max_orb_ind)
               CALL init_mo_set(active_space_env%mos_active(ispin), fm_struct=fm_struct_tmp, name="Active Space MO")
               CALL cp_fm_struct_release(fm_struct_tmp)

               ! init inactive orbitals
               IF (nspins == 2) THEN
                  nel = nelec_inactive/2
               ELSE
                  nel = nelec_inactive
               END IF
               CALL allocate_mo_set(active_space_env%mos_inactive(ispin), nao, max_orb_ind, nel, REAL(nel, KIND=dp), maxocc, 0.0_dp)
               CALL cp_fm_struct_create(fm_struct_tmp, para_env=para_env, context=context, &
                                        nrow_global=nrow_global, ncol_global=max_orb_ind)
               CALL init_mo_set(active_space_env%mos_inactive(ispin), fm_struct=fm_struct_tmp, name="Inactive Space MO")
               ! small hack: set the correct inactive occupations down below
               active_space_env%mos_inactive(ispin)%occupation_numbers = 0.0_dp
               CALL cp_fm_struct_release(fm_struct_tmp)
            END DO

            ALLOCATE (eigenvalues(max_orb_ind, nspins))
            eigenvalues = 0.0_dp
            CALL get_qs_env(qs_env, matrix_ks=ks_matrix, matrix_s=s_matrix, scf_control=scf_control)

            ! calculate virtual MOs and copy inactive and active orbitals
            IF (iw > 0) THEN
               WRITE (iw, '(/,T3,A)') "Calculating virtual MOs..."
            END IF
            DO ispin = 1, nspins
               CALL get_mo_set(mos(ispin), mo_coeff=mo_ref, nmo=nmo_available)
               nmo_virtual = max_orb_ind - nmo_available
               nmo_virtual = MAX(nmo_virtual, 0)

               NULLIFY (evals_virtual)
               ALLOCATE (evals_virtual(nmo_virtual))

               CALL cp_fm_get_info(mo_ref, context=context, para_env=para_env, &
                                   nrow_global=nrow_global)

               CALL cp_fm_struct_create(fm_struct_tmp, para_env=para_env, context=context, &
                                        nrow_global=nrow_global, ncol_global=nmo_virtual)
               CALL cp_fm_create(mo_virtual, fm_struct_tmp, name="virtual")
               CALL cp_fm_struct_release(fm_struct_tmp)
               CALL cp_fm_init_random(mo_virtual, nmo_virtual)

               NULLIFY (local_preconditioner)

               CALL ot_eigensolver(matrix_h=ks_matrix(ispin)%matrix, matrix_s=s_matrix(1)%matrix, &
                                   matrix_c_fm=mo_virtual, matrix_orthogonal_space_fm=mo_ref, &
                                   eps_gradient=scf_control%eps_lumos, &
                                   preconditioner=local_preconditioner, &
                                   iter_max=scf_control%max_iter_lumos, &
                                   size_ortho_space=nmo_available)

               CALL calculate_subspace_eigenvalues(mo_virtual, ks_matrix(ispin)%matrix, &
                                                   evals_virtual)

               ! We need to send the copy of MOs to preserve the sign
               CALL cp_fm_create(fm_dummy, mo_ref%matrix_struct)
               CALL cp_fm_to_fm(mo_ref, fm_dummy)

               CALL calculate_subspace_eigenvalues(fm_dummy, ks_matrix(ispin)%matrix, &
                                                   evals_arg=eigenvalues(:, ispin), do_rotation=.TRUE.)

               mo_set_active => active_space_env%mos_active(ispin)
               CALL get_mo_set(mo_set_active, mo_coeff=fm_target_active)
               mo_set_inactive => active_space_env%mos_inactive(ispin)
               CALL get_mo_set(mo_set_inactive, mo_coeff=fm_target_inactive)

               ! copy orbitals
               nmo_inactive_remaining = nmo_inactive
               DO i = 1, max_orb_ind
                  ! case for i being an active orbital
                  IF (ANY(active_space_env%active_orbitals(:, ispin) == i)) THEN
                     IF (i > nmo_available) THEN
                        CALL cp_fm_to_fm(mo_virtual, fm_target_active, 1, i - nmo_available, i)
                        eigenvalues(i, ispin) = evals_virtual(i - nmo_available)
                        mo_set_active%occupation_numbers(i) = 0.0
                     ELSE
                        CALL cp_fm_to_fm(fm_dummy, fm_target_active, 1, i, i)
                        mo_set_active%occupation_numbers(i) = mos(ispin)%occupation_numbers(i)
                     END IF
                     mo_set_active%eigenvalues(i) = eigenvalues(i, ispin)
                     ! if it was not an active orbital, check whether it is an inactive orbital
                  ELSEIF (nmo_inactive_remaining > 0) THEN
                     CALL cp_fm_to_fm(fm_dummy, fm_target_inactive, 1, i, i)
                     ! store on the fly the mapping of inactive orbitals
                     active_space_env%inactive_orbitals(nmo_inactive - nmo_inactive_remaining + 1, ispin) = i
                     mo_set_inactive%eigenvalues(i) = eigenvalues(i, ispin)
                     mo_set_inactive%occupation_numbers(i) = mos(ispin)%occupation_numbers(i)
                     ! hack: set homo and lumo manually
                     IF (nmo_inactive_remaining == 1) THEN
                        mo_set_inactive%homo = i
                        mo_set_inactive%lfomo = i + 1
                     END IF
                     nmo_inactive_remaining = nmo_inactive_remaining - 1
                  ELSE
                     CYCLE
                  END IF
               END DO

               ! Release
               DEALLOCATE (evals_virtual)
               CALL cp_fm_release(fm_dummy)
               CALL cp_fm_release(mo_virtual)
            END DO

            IF (iw > 0) THEN
               DO ispin = 1, nspins
                  WRITE (iw, '(/,T3,A,I3,T66,A)') "Orbital Energies and Selection for spin", ispin, "[atomic units]"

                  DO i = 1, max_orb_ind, 4
                     jm = MIN(3, max_orb_ind - i)
                     WRITE (iw, '(T4)', advance="no")
                     DO j = 0, jm
                        IF (ANY(active_space_env%active_orbitals(:, ispin) == i + j)) THEN
                           WRITE (iw, '(T3,F12.6,A5)', advance="no") eigenvalues(i + j, ispin), " [A]"
                        ELSEIF (ANY(active_space_env%inactive_orbitals(:, ispin) == i + j)) THEN
                           WRITE (iw, '(T3,F12.6,A5)', advance="no") eigenvalues(i + j, ispin), " [I]"
                        ELSE
                           WRITE (iw, '(T3,F12.6,A5)', advance="no") eigenvalues(i + j, ispin), " [V]"
                        END IF
                     END DO
                     WRITE (iw, *)
                  END DO
                  WRITE (iw, '(/,T3,A,I3)') "Active Orbital Indices for spin", ispin
                  DO i = 1, SIZE(active_space_env%active_orbitals, 1), 4
                     jm = MIN(3, SIZE(active_space_env%active_orbitals, 1) - i)
                     WRITE (iw, '(T3,4(I4))') (active_space_env%active_orbitals(i + j, ispin), j=0, jm)
                  END DO
               END DO
            END IF
            DEALLOCATE (eigenvalues)
         END IF

      CASE (wannier_projection)
         NULLIFY (loc_section, loc_print)
         loc_section => section_vals_get_subs_vals(as_input, "LOCALIZE")
         CPASSERT(ASSOCIATED(loc_section))
         loc_print => section_vals_get_subs_vals(as_input, "LOCALIZE%PRINT")
         !
         CPABORT("not yet available")
         !
      CASE (mao_projection)
         !
         CPABORT("not yet available")
         !
      END SELECT

      ! Print orbitals on Cube files
      print_orb => section_vals_get_subs_vals(as_input, "PRINT_ORBITAL_CUBES")
      CALL section_vals_get(print_orb, explicit=explicit)
      CALL section_vals_val_get(print_orb, "STOP_AFTER_CUBES", l_val=stop_after_print)
      IF (explicit) THEN
         !
         CALL print_orbital_cubes(print_orb, qs_env, active_space_env%mos_active)
         !
         IF (stop_after_print) THEN

            IF (iw > 0) THEN
               WRITE (iw, '(/,T2,A)') &
                  '!----------------- Early End of Active Space Interface -----------------------!'
            END IF

            CALL timestop(handle)

            RETURN
         END IF
      END IF

      ! calculate inactive density matrix
      CALL get_qs_env(qs_env, rho=rho)
      CALL qs_rho_get(rho, rho_ao=rho_ao)
      CPASSERT(ASSOCIATED(rho_ao))
      CALL dbcsr_allocate_matrix_set(active_space_env%pmat_inactive, nspins)
      DO ispin = 1, nspins
         ALLOCATE (denmat)
         CALL dbcsr_copy(denmat, rho_ao(ispin)%matrix)
         mo_set => active_space_env%mos_inactive(ispin)
         CALL calculate_density_matrix(mo_set, denmat)
         active_space_env%pmat_inactive(ispin)%matrix => denmat
      END DO

      ! read in ERI parameters
      CALL section_vals_val_get(as_input, "ERI%METHOD", i_val=eri_method)
      active_space_env%eri%method = eri_method
      CALL section_vals_val_get(as_input, "ERI%OPERATOR", i_val=eri_operator, explicit=ex_operator)
      active_space_env%eri%operator = eri_operator
      CALL section_vals_val_get(as_input, "ERI%OMEGA", r_val=eri_op_omega, explicit=ex_omega)
      active_space_env%eri%omega = eri_op_omega
      CALL section_vals_val_get(as_input, "ERI%CUTOFF_RADIUS", r_val=eri_rcut, explicit=ex_rcut)
      active_space_env%eri%cutoff_radius = eri_rcut  ! this is already converted to bohr!
      CALL section_vals_val_get(as_input, "ERI%PERIODICITY", i_vals=invals, explicit=ex_perd)
      CALL section_vals_val_get(as_input, "ERI%EPS_INTEGRAL", r_val=eri_eps_int)
      active_space_env%eri%eps_integral = eri_eps_int
      ! if eri periodicity is explicitly set, we use it, otherwise we use the cell periodicity
      IF (ex_perd) THEN
         IF (SIZE(invals) == 1) THEN
            active_space_env%eri%periodicity(1:3) = invals(1)
         ELSE
            active_space_env%eri%periodicity(1:3) = invals(1:3)
         END IF
      ELSE
         CALL get_qs_env(qs_env, cell=cell)
         active_space_env%eri%periodicity(1:3) = cell%perd(1:3)
      END IF
      IF (iw > 0) THEN
         WRITE (iw, '(/,T3,A)') "Calculation of Electron Repulsion Integrals"

         SELECT CASE (eri_method)
         CASE (eri_method_full_gpw)
            WRITE (iw, '(T3,A,T50,A)') "Integration method", "GPW Fourier transform over MOs"
         CASE (eri_method_gpw_ht)
            WRITE (iw, '(T3,A,T44,A)') "Integration method", "Half transformed integrals from GPW"
         CASE DEFAULT
            CPABORT("Unknown ERI method")
         END SELECT

         SELECT CASE (eri_operator)
         CASE (eri_operator_coulomb)
            WRITE (iw, '(T3,A,T73,A)') "ERI operator", "Coulomb"

         CASE (eri_operator_yukawa)
            WRITE (iw, '(T3,A,T74,A)') "ERI operator", "Yukawa"
            IF (.NOT. ex_omega) CALL cp_abort(__LOCATION__, &
                                              "Yukawa operator requires OMEGA to be explicitly set")
            WRITE (iw, '(T3,A,T66,F14.3)') "ERI operator parameter OMEGA", eri_op_omega

         CASE (eri_operator_erf)
            WRITE (iw, '(T3,A,T63,A)') "ERI operator", "Longrange Coulomb"
            IF (.NOT. ex_omega) CALL cp_abort(__LOCATION__, &
                                              "Longrange operator requires OMEGA to be explicitly set")
            WRITE (iw, '(T3,A,T66,F14.3)') "ERI operator parameter OMEGA", eri_op_omega

         CASE (eri_operator_erfc)
            WRITE (iw, '(T3,A,T62,A)') "ERI operator", "Shortrange Coulomb"
            IF (.NOT. ex_omega) CALL cp_abort(__LOCATION__, &
                                              "Shortrange operator requires OMEGA to be explicitly set")
            WRITE (iw, '(T3,A,T66,F14.3)') "ERI operator parameter OMEGA", eri_op_omega

         CASE (eri_operator_trunc)
            WRITE (iw, '(T3,A,T63,A)') "ERI operator", "Truncated Coulomb"
            IF (.NOT. ex_rcut) CALL cp_abort(__LOCATION__, &
                                             "Cutoff radius not specified for trunc. Coulomb operator")
            WRITE (iw, '(T3,A,T66,F14.3)') "ERI operator cutoff radius (au)", eri_rcut

         CASE (eri_operator_lr_trunc)
            WRITE (iw, '(T3,A,T53,A)') "ERI operator", "Longrange truncated Coulomb"
            IF (.NOT. ex_rcut) CALL cp_abort(__LOCATION__, &
                                             "Cutoff radius not specified for trunc. longrange operator")
            WRITE (iw, '(T3,A,T66,F14.3)') "ERI operator cutoff radius (au)", eri_rcut
            IF (.NOT. ex_omega) CALL cp_abort(__LOCATION__, &
                                              "LR truncated operator requires OMEGA to be explicitly set")
            WRITE (iw, '(T3,A,T66,F14.3)') "ERI operator parameter OMEGA", eri_op_omega
            IF (eri_op_omega < 0.01_dp) THEN
               CPABORT("LR truncated operator requires OMEGA >= 0.01 to be stable")
            END IF

         CASE DEFAULT
            CPABORT("Unknown ERI operator")

         END SELECT

         WRITE (iw, '(T3,A,T68,E12.4)') "Accuracy of ERIs", eri_eps_int
         WRITE (iw, '(T3,A,T71,3I3)') "Periodicity", active_space_env%eri%periodicity(1:3)

         ! TODO: should be moved after ERI calculation, as it depends on screening
         IF (nspins < 2) THEN
            WRITE (iw, '(T3,A,T68,I12)') "Total Number of ERI", (nmo_active**4)/8
         ELSE
            WRITE (iw, '(T3,A,T68,I12)') "Total Number of ERI (aa|aa)", (nmo_active**4)/8
            WRITE (iw, '(T3,A,T68,I12)') "Total Number of ERI (bb|bb)", (nmo_active**4)/8
            WRITE (iw, '(T3,A,T68,I12)') "Total Number of ERI (aa|bb)", (nmo_active**4)/4
         END IF
      END IF

      ! allocate container for integrals (CSR matrix)
      CALL get_qs_env(qs_env, para_env=para_env)
      m = (nspins*(nspins + 1))/2
      ! With ROHF/ROKS, we need ERIs from only a single set of orbitals
      IF (dft_control%roks) m = 1
      ALLOCATE (active_space_env%eri%eri(m))
      DO i = 1, m
         CALL get_mo_set(active_space_env%mos_active(1), nmo=nmo)
         ALLOCATE (active_space_env%eri%eri(i)%csr_mat)
         eri_mat => active_space_env%eri%eri(i)%csr_mat
         IF (i == 1) THEN
            n1 = nmo
            n2 = nmo
         ELSEIF (i == 2) THEN
            n1 = nmo
            n2 = nmo
         ELSE
            n1 = nmo
            n2 = nmo
         END IF
         nn1 = (n1*(n1 + 1))/2
         nn2 = (n2*(n2 + 1))/2
         CALL dbcsr_csr_create(eri_mat, nn1, nn2, 0_int_8, 0, 0, para_env%get_handle())
         active_space_env%eri%norb = nmo
      END DO

      SELECT CASE (eri_method)
      CASE (eri_method_full_gpw, eri_method_gpw_ht)
         CALL section_vals_val_get(as_input, "ERI_GPW%EPS_GRID", r_val=eri_eps_grid)
         active_space_env%eri%eri_gpw%eps_grid = eri_eps_grid
         CALL section_vals_val_get(as_input, "ERI_GPW%EPS_FILTER", r_val=eri_eps_filter)
         active_space_env%eri%eri_gpw%eps_filter = eri_eps_filter
         CALL section_vals_val_get(as_input, "ERI_GPW%CUTOFF", r_val=eri_gpw_cutoff)
         active_space_env%eri%eri_gpw%cutoff = eri_gpw_cutoff
         CALL section_vals_val_get(as_input, "ERI_GPW%REL_CUTOFF", r_val=eri_rel_cutoff)
         active_space_env%eri%eri_gpw%rel_cutoff = eri_rel_cutoff
         CALL section_vals_val_get(as_input, "ERI_GPW%PRINT_LEVEL", i_val=eri_print)
         active_space_env%eri%eri_gpw%print_level = eri_print
         CALL section_vals_val_get(as_input, "ERI_GPW%STORE_WFN", l_val=store_wfn)
         active_space_env%eri%eri_gpw%store_wfn = store_wfn
         CALL section_vals_val_get(as_input, "ERI_GPW%GROUP_SIZE", i_val=group_size)
         active_space_env%eri%eri_gpw%group_size = group_size
         ! Always redo Poisson solver for now
         active_space_env%eri%eri_gpw%redo_poisson = .TRUE.
         ! active_space_env%eri%eri_gpw%redo_poisson = (ex_operator .OR. ex_perd)
         IF (iw > 0) THEN
            WRITE (iw, '(/,T2,A,T71,F10.1)') "ERI_GPW| Energy cutoff [Ry]", eri_gpw_cutoff
            WRITE (iw, '(T2,A,T71,F10.1)') "ERI_GPW| Relative energy cutoff [Ry]", eri_rel_cutoff
         END IF
         !
         CALL calculate_eri_gpw(active_space_env%mos_active, active_space_env%active_orbitals, active_space_env%eri, qs_env, iw, &
                                dft_control%roks)
         !
      CASE DEFAULT
         CPABORT("Unknown ERI method")
      END SELECT
      IF (iw > 0) THEN
         DO isp = 1, SIZE(active_space_env%eri%eri)
            eri_mat => active_space_env%eri%eri(isp)%csr_mat
            nze_percentage = 100.0_dp*(REAL(eri_mat%nze_total, KIND=dp) &
                                       /REAL(eri_mat%nrows_total, KIND=dp))/REAL(eri_mat%ncols_total, KIND=dp)
            WRITE (iw, '(/,T2,A,I2,T30,A,T68,I12)') "ERI_GPW| Spinmatrix:", isp, &
               "Number of  CSR non-zero elements:", eri_mat%nze_total
            WRITE (iw, '(T2,A,I2,T30,A,T68,F12.4)') "ERI_GPW| Spinmatrix:", isp, &
               "Percentage CSR non-zero elements:", nze_percentage
            WRITE (iw, '(T2,A,I2,T30,A,T68,I12)') "ERI_GPW| Spinmatrix:", isp, &
               "nrows_total", eri_mat%nrows_total
            WRITE (iw, '(T2,A,I2,T30,A,T68,I12)') "ERI_GPW| Spinmatrix:", isp, &
               "ncols_total", eri_mat%ncols_total
            WRITE (iw, '(T2,A,I2,T30,A,T68,I12)') "ERI_GPW| Spinmatrix:", isp, &
               "nrows_local", eri_mat%nrows_local
         END DO
         CALL m_flush(iw)
      END IF
      CALL para_env%sync()

      ! set the reference active space density matrix
      nspins = active_space_env%nspins
      ALLOCATE (active_space_env%p_active(nspins))
      DO isp = 1, nspins
         mo_set => active_space_env%mos_active(isp)
         CALL get_mo_set(mo_set, mo_coeff=mo_coeff, nmo=nmo)
         CALL create_subspace_matrix(mo_coeff, active_space_env%p_active(isp), nmo)
      END DO
      SELECT CASE (mselect)
      CASE DEFAULT
         CPABORT("Unknown orbital selection method")
      CASE (casci_canonical, manual_selection)
         focc = 2.0_dp
         IF (nspins == 2) focc = 1.0_dp
         DO isp = 1, nspins
            fmat => active_space_env%p_active(isp)
            CALL cp_fm_set_all(fmat, alpha=0.0_dp)
            IF (nspins == 2) THEN
               IF (isp == 1) THEN
                  nel = (active_space_env%nelec_active + active_space_env%multiplicity - 1)/2
               ELSE
                  nel = (active_space_env%nelec_active - active_space_env%multiplicity + 1)/2
               END IF
            ELSE
               nel = active_space_env%nelec_active
            END IF
            DO i = 1, nmo_active
               m = active_space_env%active_orbitals(i, isp)
               fel = MIN(focc, REAL(nel, KIND=dp))
               CALL cp_fm_set_element(fmat, m, m, fel)
               nel = nel - NINT(fel)
               nel = MAX(nel, 0)
            END DO
         END DO
      CASE (wannier_projection)
         CPABORT("NOT IMPLEMENTED")
      CASE (mao_projection)
         CPABORT("NOT IMPLEMENTED")
      END SELECT

      ! compute alpha-beta overlap matrix in case of spin-polarized calculation
      CALL calculate_spin_pol_overlap(active_space_env%mos_active, qs_env, active_space_env)

      ! figure out if we have a new xc section for the AS
      xc_section => section_vals_get_subs_vals(input, "DFT%ACTIVE_SPACE%XC")
      explicit = .FALSE.
      IF (ASSOCIATED(xc_section)) CALL section_vals_get(xc_section, explicit=explicit)

      ! rebuild KS matrix if needed
      IF (explicit) THEN
         ! release the hfx data if it was part of the SCF functional
         IF (ASSOCIATED(qs_env%x_data)) CALL hfx_release(qs_env%x_data)
         ! also release the admm environment in case we are using admm
         IF (ASSOCIATED(qs_env%admm_env)) CALL admm_env_release(qs_env%admm_env)

         CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set, &
                         particle_set=particle_set, cell=cell, ks_env=ks_env)
         IF (dft_control%do_admm) THEN
            basis_type = 'AUX_FIT'
         ELSE
            basis_type = 'ORB'
         END IF
         hfx_section => section_vals_get_subs_vals(xc_section, "HF")
         CALL hfx_create(qs_env%x_data, para_env, hfx_section, atomic_kind_set, &
                         qs_kind_set, particle_set, dft_control, cell, orb_basis=basis_type, &
                         nelectron_total=nelec_total)

         qs_env%requires_matrix_vxc = .TRUE.  ! needs to be set only once

         ! a bit of a hack: this forces a new re-init of HFX
         CALL set_ks_env(ks_env, s_mstruct_changed=.TRUE.)
         CALL qs_ks_build_kohn_sham_matrix(qs_env, calculate_forces=.FALSE., &
                                           just_energy=.FALSE., &
                                           ext_xc_section=xc_section)
         ! we need to reset it to false
         CALL set_ks_env(ks_env, s_mstruct_changed=.FALSE.)
      ELSE
         xc_section => section_vals_get_subs_vals(input, "DFT%XC")
      END IF
      ! set the xc_section
      active_space_env%xc_section => xc_section

      CALL get_qs_env(qs_env, energy=energy)
      ! transform KS/Fock, Vxc and Hcore to AS MO basis
      CALL calculate_operators(active_space_env%mos_active, qs_env, active_space_env)
      ! set the reference energy in the active space
      active_space_env%energy_ref = energy%total
      ! calculate inactive energy and embedding potential
      CALL subspace_fock_matrix(active_space_env, dft_control%roks)

      ! associate the active space environment with the qs environment
      CALL set_qs_env(qs_env, active_space=active_space_env)

      ! Perform the embedding calculation only if qiskit is specified
      CALL section_vals_val_get(as_input, "AS_SOLVER", i_val=as_solver)
      SELECT CASE (as_solver)
      CASE (no_solver)
         IF (iw > 0) THEN
            WRITE (iw, '(/,T3,A)') "No active space solver specified, skipping embedding calculation"
            CALL m_flush(iw)
         END IF
         CALL para_env%sync()
      CASE (qiskit_solver)
         CALL rsdft_embedding(qs_env, active_space_env, as_input)
         CALL qs_scf_compute_properties(qs_env, wf_type="MC-DFT", do_mp2=.FALSE.)
      CASE DEFAULT
         CPABORT("Unknown active space solver")
      END SELECT

      ! Output a FCIDUMP file if requested
      IF (active_space_env%fcidump) CALL fcidump(active_space_env, as_input, dft_control%roks)

      ! Output a QCSchema file if requested
      IF (active_space_env%qcschema) THEN
         CALL qcschema_env_create(qcschema_env, qs_env)
         CALL qcschema_to_hdf5(qcschema_env, active_space_env%qcschema_filename)
         CALL qcschema_env_release(qcschema_env)
      END IF

      IF (iw > 0) THEN
         WRITE (iw, '(/,T2,A)') &
            '!-------------------- End of Active Space Interface --------------------------!'
         CALL m_flush(iw)
      END IF
      CALL para_env%sync()

      CALL timestop(handle)

   END SUBROUTINE active_space_main

! **************************************************************************************************
!> \brief computes the alpha-beta overlap within the active subspace
!> \param mos the molecular orbital set within the active subspace
!> \param qs_env ...
!> \param active_space_env ...
!> \par History
!>      04.2016 created [JGH]
! **************************************************************************************************
   SUBROUTINE calculate_spin_pol_overlap(mos, qs_env, active_space_env)

      TYPE(mo_set_type), DIMENSION(:), INTENT(IN)        :: mos
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(active_space_type), POINTER                   :: active_space_env

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_spin_pol_overlap'

      INTEGER                                            :: handle, nmo, nspins
      TYPE(cp_fm_type), POINTER                          :: mo_coeff_a, mo_coeff_b
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: s_matrix

      CALL timeset(routineN, handle)

      nspins = active_space_env%nspins

      ! overlap in AO
      IF (nspins > 1) THEN
         CALL get_qs_env(qs_env, matrix_s=s_matrix)
         ALLOCATE (active_space_env%sab_sub(1))

         CALL get_mo_set(mo_set=mos(1), mo_coeff=mo_coeff_a, nmo=nmo)
         CALL get_mo_set(mo_set=mos(2), mo_coeff=mo_coeff_b, nmo=nmo)
         CALL subspace_operator(mo_coeff_a, nmo, s_matrix(1)%matrix, active_space_env%sab_sub(1), mo_coeff_b)
      END IF

      CALL timestop(handle)

   END SUBROUTINE calculate_spin_pol_overlap

! **************************************************************************************************
!> \brief computes the one-electron operators in the subspace of the provided orbital set
!> \param mos the molecular orbital set within the active subspace
!> \param qs_env ...
!> \param active_space_env ...
!> \par History
!>      04.2016 created [JGH]
! **************************************************************************************************
   SUBROUTINE calculate_operators(mos, qs_env, active_space_env)

      TYPE(mo_set_type), DIMENSION(:), INTENT(IN)        :: mos
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(active_space_type), POINTER                   :: active_space_env

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_operators'

      INTEGER                                            :: handle, ispin, nmo, nspins
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: h_matrix, ks_matrix

      CALL timeset(routineN, handle)

      nspins = active_space_env%nspins

      ! Kohn-Sham / Fock operator
      CALL cp_fm_release(active_space_env%ks_sub)
      CALL get_qs_env(qs_env, matrix_ks_kp=ks_matrix)
      ALLOCATE (active_space_env%ks_sub(nspins))
      DO ispin = 1, nspins
         CALL get_mo_set(mo_set=mos(ispin), mo_coeff=mo_coeff, nmo=nmo)
         CALL subspace_operator(mo_coeff, nmo, ks_matrix(ispin, 1)%matrix, active_space_env%ks_sub(ispin))
      END DO

      ! Core Hamiltonian
      CALL cp_fm_release(active_space_env%h_sub)

      NULLIFY (h_matrix)
      CALL get_qs_env(qs_env=qs_env, matrix_h_kp=h_matrix)
      ALLOCATE (active_space_env%h_sub(nspins))
      DO ispin = 1, nspins
         CALL get_mo_set(mo_set=mos(ispin), mo_coeff=mo_coeff, nmo=nmo)
         CALL subspace_operator(mo_coeff, nmo, h_matrix(1, 1)%matrix, active_space_env%h_sub(ispin))
      END DO

      CALL timestop(handle)

   END SUBROUTINE calculate_operators

! **************************************************************************************************
!> \brief computes a one-electron operator in the subspace of the provided orbital set
!> \param mo_coeff the orbital coefficient matrix
!> \param nmo the number of subspace orbitals
!> \param op_matrix operator matrix in AO basis
!> \param op_sub operator in orbital basis
!> \param mo_coeff_b the beta orbital coefficients
!> \par History
!>      04.2016 created [JGH]
! **************************************************************************************************
   SUBROUTINE subspace_operator(mo_coeff, nmo, op_matrix, op_sub, mo_coeff_b)

      TYPE(cp_fm_type), INTENT(IN)                       :: mo_coeff
      INTEGER, INTENT(IN)                                :: nmo
      TYPE(dbcsr_type), POINTER                          :: op_matrix
      TYPE(cp_fm_type), INTENT(INOUT)                    :: op_sub
      TYPE(cp_fm_type), INTENT(IN), OPTIONAL             :: mo_coeff_b

      CHARACTER(len=*), PARAMETER                        :: routineN = 'subspace_operator'

      INTEGER                                            :: handle, ncol, nrow
      TYPE(cp_fm_type)                                   :: vectors

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(matrix=mo_coeff, ncol_global=ncol, nrow_global=nrow)
      CPASSERT(nmo <= ncol)

      IF (nmo > 0) THEN
         CALL cp_fm_create(vectors, mo_coeff%matrix_struct, "vectors")
         CALL create_subspace_matrix(mo_coeff, op_sub, nmo)

         IF (PRESENT(mo_coeff_b)) THEN
            ! if beta orbitals are present, compute the cross alpha_beta term
            CALL cp_dbcsr_sm_fm_multiply(op_matrix, mo_coeff_b, vectors, nmo)
         ELSE
            ! otherwise the same spin, whatever that is
            CALL cp_dbcsr_sm_fm_multiply(op_matrix, mo_coeff, vectors, nmo)
         END IF

         CALL parallel_gemm('T', 'N', nmo, nmo, nrow, 1.0_dp, mo_coeff, vectors, 0.0_dp, op_sub)
         CALL cp_fm_release(vectors)
      END IF

      CALL timestop(handle)

   END SUBROUTINE subspace_operator

! **************************************************************************************************
!> \brief creates a matrix of subspace size
!> \param orbitals the orbital coefficient matrix
!> \param op_sub operator in orbital basis
!> \param n the number of orbitals
!> \par History
!>      04.2016 created [JGH]
! **************************************************************************************************
   SUBROUTINE create_subspace_matrix(orbitals, op_sub, n)

      TYPE(cp_fm_type), INTENT(IN)                       :: orbitals
      TYPE(cp_fm_type), INTENT(OUT)                      :: op_sub
      INTEGER, INTENT(IN)                                :: n

      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct

      IF (n > 0) THEN

         NULLIFY (fm_struct)
         CALL cp_fm_struct_create(fm_struct, nrow_global=n, ncol_global=n, &
                                  para_env=orbitals%matrix_struct%para_env, &
                                  context=orbitals%matrix_struct%context)
         CALL cp_fm_create(op_sub, fm_struct, name="Subspace operator")
         CALL cp_fm_struct_release(fm_struct)

      END IF

   END SUBROUTINE create_subspace_matrix

! **************************************************************************************************
!> \brief computes the electron repulsion integrals using the GPW technology
!> \param mos the molecular orbital set within the active subspace
!> \param orbitals ...
!> \param eri_env ...
!> \param qs_env ...
!> \param iw ...
!> \param restricted ...
!> \par History
!>      04.2016 created [JGH]
! **************************************************************************************************
   SUBROUTINE calculate_eri_gpw(mos, orbitals, eri_env, qs_env, iw, restricted)

      TYPE(mo_set_type), DIMENSION(:), INTENT(IN)        :: mos
      INTEGER, DIMENSION(:, :), POINTER                  :: orbitals
      TYPE(eri_type)                                     :: eri_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: iw
      LOGICAL, INTENT(IN)                                :: restricted

      CHARACTER(len=*), PARAMETER                        :: routineN = 'calculate_eri_gpw'

      INTEGER :: col_local, color, handle, i1, i2, i3, i4, i_multigrid, icount2, intcount, &
         irange(2), isp, isp1, isp2, ispin, iwa1, iwa12, iwa2, iwb1, iwb12, iwb2, iwbs, iwbt, &
         iwfn, n_multigrid, ncol_global, ncol_local, nmm, nmo, nmo1, nmo2, nrow_global, &
         nrow_local, nspins, number_of_subgroups, nx, row_local, stored_integrals
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: eri_index
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: print1, print2, &
                                                            skip_load_balance_distributed
      REAL(KIND=dp)                                      :: dvol, erint, pair_int, &
                                                            progression_factor, rc, rsize, t1, t2
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eri
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env, blacs_env_sub
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_matrix_pq_rnu, fm_matrix_pq_rs, &
                                                            fm_mo_coeff_as
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dbcsr_p_type)                                 :: mat_munu
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:)        :: matrix_pq_rnu, mo_coeff_as
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb_sub
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_c1d_gs_type)                               :: pot_g, rho_g
      TYPE(pw_env_type), POINTER                         :: pw_env_sub
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_r3d_rs_type)                               :: rho_r, wfn_r
      TYPE(pw_r3d_rs_type), ALLOCATABLE, &
         DIMENSION(:, :), TARGET                         :: wfn_a
      TYPE(pw_r3d_rs_type), POINTER                      :: wfn1, wfn2, wfn3, wfn4
      TYPE(qs_control_type), POINTER                     :: qs_control, qs_control_old
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(task_list_type), POINTER                      :: task_list_sub

      CALL timeset(routineN, handle)

      IF (iw > 0) t1 = m_walltime()

      ! print levels
      SELECT CASE (eri_env%eri_gpw%print_level)
      CASE (silent_print_level)
         print1 = .FALSE.
         print2 = .FALSE.
      CASE (low_print_level)
         print1 = .FALSE.
         print2 = .FALSE.
      CASE (medium_print_level)
         print1 = .TRUE.
         print2 = .FALSE.
      CASE (high_print_level)
         print1 = .TRUE.
         print2 = .TRUE.
      CASE (debug_print_level)
         print1 = .TRUE.
         print2 = .TRUE.
      CASE DEFAULT
         ! do nothing
      END SELECT

      ! Check the input group
      CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env)
      IF (eri_env%eri_gpw%group_size < 1) eri_env%eri_gpw%group_size = para_env%num_pe
      IF (MOD(para_env%num_pe, eri_env%eri_gpw%group_size) /= 0) &
         CPABORT("Group size must be a divisor of the total number of processes!")
      ! Create a new para_env or reuse the old one
      IF (eri_env%eri_gpw%group_size == para_env%num_pe) THEN
         eri_env%para_env_sub => para_env
         CALL eri_env%para_env_sub%retain()
         blacs_env_sub => blacs_env
         CALL blacs_env_sub%retain()
         number_of_subgroups = 1
         color = 0
      ELSE
         number_of_subgroups = para_env%num_pe/eri_env%eri_gpw%group_size
         color = para_env%mepos/eri_env%eri_gpw%group_size
         ALLOCATE (eri_env%para_env_sub)
         CALL eri_env%para_env_sub%from_split(para_env, color)
         NULLIFY (blacs_env_sub)
         CALL cp_blacs_env_create(blacs_env_sub, eri_env%para_env_sub, BLACS_GRID_SQUARE, .TRUE.)
      END IF
      CALL eri_env%comm_exchange%from_split(para_env, eri_env%para_env_sub%mepos)

      ! This should be done differently! Copied from MP2 code
      CALL get_qs_env(qs_env, dft_control=dft_control)
      ALLOCATE (qs_control)
      qs_control_old => dft_control%qs_control
      qs_control = qs_control_old
      dft_control%qs_control => qs_control
      progression_factor = qs_control%progression_factor
      n_multigrid = SIZE(qs_control%e_cutoff)
      nspins = SIZE(mos)
      ! In case of ROHF/ROKS, we assume the orbital coefficients in both spin channels to be the same
      ! and save operations by calculating ERIs from only one spin channel
      IF (restricted) nspins = 1
      ! Allocate new cutoffs (just in private qs_control, not in qs_control_old)
      ALLOCATE (qs_control%e_cutoff(n_multigrid))

      qs_control%cutoff = eri_env%eri_gpw%cutoff*0.5_dp
      qs_control%e_cutoff(1) = qs_control%cutoff
      DO i_multigrid = 2, n_multigrid
         qs_control%e_cutoff(i_multigrid) = qs_control%e_cutoff(i_multigrid - 1) &
                                            /progression_factor
      END DO
      qs_control%relative_cutoff = eri_env%eri_gpw%rel_cutoff*0.5_dp

      ! For now, we will distribute neighbor lists etc. within the global communicator
      CALL get_qs_env(qs_env, ks_env=ks_env)
      CALL create_mat_munu(mat_munu, qs_env, eri_env%eri_gpw%eps_grid, blacs_env_sub, sab_orb_sub=sab_orb_sub, &
                           do_alloc_blocks_from_nbl=.TRUE., dbcsr_sym_type=dbcsr_type_symmetric)
      CALL dbcsr_set(mat_munu%matrix, 0.0_dp)

      ! Generate the appropriate pw_env
      NULLIFY (pw_env_sub)
      CALL pw_env_create(pw_env_sub)
      CALL pw_env_rebuild(pw_env_sub, qs_env, external_para_env=eri_env%para_env_sub)
      CALL pw_env_get(pw_env_sub, auxbas_pw_pool=auxbas_pw_pool, poisson_env=poisson_env)

      ! TODO: maybe we can let `pw_env_rebuild` do what we manually overwrite here?
      IF (eri_env%eri_gpw%redo_poisson) THEN
         ! We need to rebuild the Poisson solver on the fly
         IF (SUM(eri_env%periodicity) /= 0) THEN
            poisson_env%parameters%solver = pw_poisson_periodic
         ELSE
            poisson_env%parameters%solver = pw_poisson_analytic
         END IF
         poisson_env%parameters%periodic = eri_env%periodicity

         ! Rebuilds the poisson green (influence) function according
         ! to the poisson solver and parameters set so far.
         ! Also sets the variable poisson_env%rebuild to .FALSE.
         CALL pw_poisson_rebuild(poisson_env)

         ! set the cutoff radius for the Greens function in case we use ANALYTIC Poisson solver
         CALL get_qs_env(qs_env, cell=cell)
         rc = cell%hmat(1, 1)
         DO iwa1 = 1, 3
            ! TODO: I think this is not the largest possible radius inscribed in the cell
            rc = MIN(rc, 0.5_dp*cell%hmat(iwa1, iwa1))
         END DO
         poisson_env%green_fft%radius = rc

         ! Overwrite the Greens function with the one we want
         CALL pw_eri_green_create(poisson_env%green_fft, eri_env)

         IF (iw > 0) THEN
            CALL get_qs_env(qs_env, cell=cell)
            IF (SUM(cell%perd) /= SUM(eri_env%periodicity)) THEN
               IF (SUM(eri_env%periodicity) /= 0) THEN
                  WRITE (UNIT=iw, FMT="(/,T2,A,T51,A30)") &
                     "ERI_GPW| Switching Poisson solver to", "PERIODIC"
               ELSE
                  WRITE (UNIT=iw, FMT="(/,T2,A,T51,A30)") &
                     "ERI_GPW| Switching Poisson solver to", "ANALYTIC"
               END IF
            END IF
            ! print out the Greens function to check it matches the Poisson solver
            SELECT CASE (poisson_env%green_fft%method)
            CASE (PERIODIC3D)
               WRITE (UNIT=iw, FMT="(T2,A,T51,A30)") &
                  "ERI_GPW| Poisson Greens function", "PERIODIC"
            CASE (ANALYTIC0D)
               WRITE (UNIT=iw, FMT="(T2,A,T51,A30)") &
                  "ERI_GPW| Poisson Greens function", "ANALYTIC"
               WRITE (UNIT=iw, FMT="(T2,A,T71,F10.4)") "ERI_GPW| Poisson cutoff radius", &
                  poisson_env%green_fft%radius*angstrom
            CASE DEFAULT
               CPABORT("Wrong Greens function setup")
            END SELECT
         END IF
      END IF

      ALLOCATE (mo_coeff_as(nspins), fm_mo_coeff_as(nspins))
      DO ispin = 1, nspins
         BLOCK
            REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: C
            INTEGER :: nmo
            TYPE(group_dist_d1_type) :: gd_array
            TYPE(cp_fm_type), POINTER :: mo_coeff
            CALL get_mo_set(mos(ispin), mo_coeff=mo_coeff, nmo=nmo)
            CALL grep_rows_in_subgroups(para_env, eri_env%para_env_sub, mo_coeff, gd_array, C)

            CALL build_dbcsr_from_rows(eri_env%para_env_sub, mo_coeff_as(ispin), &
                                       C(:, :nmo), mat_munu%matrix, gd_array, eri_env%eri_gpw%eps_filter)
            CALL release_group_dist(gd_array)
            DEALLOCATE (C)
         END BLOCK

         CALL dbcsr_get_info(mo_coeff_as(ispin), nfullrows_total=nrow_global, nfullcols_total=ncol_global)

         NULLIFY (fm_struct)
         CALL cp_fm_struct_create(fm_struct, context=blacs_env_sub, para_env=eri_env%para_env_sub, &
                                  nrow_global=nrow_global, ncol_global=ncol_global)
         CALL cp_fm_create(fm_mo_coeff_as(ispin), fm_struct)
         CALL cp_fm_struct_release(fm_struct)

         CALL copy_dbcsr_to_fm(mo_coeff_as(ispin), fm_mo_coeff_as(ispin))
      END DO

      IF (eri_env%method == eri_method_gpw_ht) THEN
         ! We need a task list
         NULLIFY (task_list_sub)
         skip_load_balance_distributed = dft_control%qs_control%skip_load_balance_distributed
         CALL allocate_task_list(task_list_sub)
         CALL generate_qs_task_list(ks_env, task_list_sub, basis_type="ORB", &
                                    reorder_rs_grid_ranks=.TRUE., &
                                    skip_load_balance_distributed=skip_load_balance_distributed, &
                                    pw_env_external=pw_env_sub, sab_orb_external=sab_orb_sub)

         ! Create sparse matrices carrying the matrix products, Code borrowed from the MP2 GPW method
         ! Create equal distributions for them (no sparsity present)
         ! We use the routines from mp2 suggesting that one may replicate the grids later for better performance
         ALLOCATE (matrix_pq_rnu(nspins), fm_matrix_pq_rnu(nspins), fm_matrix_pq_rs(nspins))
         DO ispin = 1, nspins
            CALL dbcsr_create(matrix_pq_rnu(ispin), template=mo_coeff_as(ispin))
            CALL dbcsr_set(matrix_pq_rnu(ispin), 0.0_dp)

            CALL dbcsr_get_info(matrix_pq_rnu(ispin), nfullrows_total=nrow_global, nfullcols_total=ncol_global)

            NULLIFY (fm_struct)
            CALL cp_fm_struct_create(fm_struct, context=blacs_env_sub, para_env=eri_env%para_env_sub, &
                                     nrow_global=nrow_global, ncol_global=ncol_global)
            CALL cp_fm_create(fm_matrix_pq_rnu(ispin), fm_struct)
            CALL cp_fm_struct_release(fm_struct)

            NULLIFY (fm_struct)
            CALL cp_fm_struct_create(fm_struct, context=blacs_env_sub, para_env=eri_env%para_env_sub, &
                                     nrow_global=ncol_global, ncol_global=ncol_global)
            CALL cp_fm_create(fm_matrix_pq_rs(ispin), fm_struct)
            CALL cp_fm_struct_release(fm_struct)
         END DO

         ! Copy the active space of the MOs into DBCSR matrices
      END IF

      CALL auxbas_pw_pool%create_pw(wfn_r)
      CALL auxbas_pw_pool%create_pw(rho_g)
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, cell=cell, &
                      particle_set=particle_set, atomic_kind_set=atomic_kind_set)

      ! pre-calculate wavefunctions on reals space grid
      nspins = SIZE(mos)
      ! In case of ROHF/ROKS, we assume the orbital coefficients in both spin channels to be the same
      ! and save operations by calculating ERIs from only one spin channel
      IF (restricted) nspins = 1
      IF (eri_env%eri_gpw%store_wfn) THEN
         ! pre-calculate wavefunctions on reals space grid
         rsize = 0.0_dp
         nmo = 0
         DO ispin = 1, nspins
            CALL get_mo_set(mo_set=mos(ispin), nmo=nx)
            nmo = MAX(nmo, nx)
            rsize = REAL(SIZE(wfn_r%array), KIND=dp)*nx
         END DO
         IF (print1 .AND. iw > 0) THEN
            rsize = rsize*8._dp/1000000._dp
            WRITE (iw, "(T2,'ERI_GPW|',' Store active orbitals on real space grid ',T66,F12.3,' MB')") rsize
         END IF
         ALLOCATE (wfn_a(nmo, nspins))
         DO ispin = 1, nspins
            CALL get_mo_set(mo_set=mos(ispin), mo_coeff=mo_coeff, nmo=nmo)
            DO i1 = 1, SIZE(orbitals, 1)
               iwfn = orbitals(i1, ispin)
               CALL auxbas_pw_pool%create_pw(wfn_a(iwfn, ispin))
               CALL calculate_wavefunction(mo_coeff, iwfn, wfn_a(iwfn, ispin), rho_g, atomic_kind_set, &
                                           qs_kind_set, cell, dft_control, particle_set, pw_env_sub)
               IF (print2 .AND. iw > 0) THEN
                  WRITE (iw, "(T2,'ERI_GPW|',' Orbital stored ',I4,'  Spin ',i1)") iwfn, ispin
               END IF
            END DO
         END DO
      ELSE
         ! Even if we do not store all WFNs, we still need containers for the functions to store
         ALLOCATE (wfn1, wfn2)
         CALL auxbas_pw_pool%create_pw(wfn1)
         CALL auxbas_pw_pool%create_pw(wfn2)
         IF (eri_env%method /= eri_method_gpw_ht) THEN
            ALLOCATE (wfn3, wfn4)
            CALL auxbas_pw_pool%create_pw(wfn3)
            CALL auxbas_pw_pool%create_pw(wfn4)
         END IF
      END IF

      ! get some of the grids ready
      CALL auxbas_pw_pool%create_pw(rho_r)
      CALL auxbas_pw_pool%create_pw(pot_g)

      ! run the FFT once, to set up buffers and to take into account the memory
      CALL pw_zero(rho_r)
      CALL pw_transfer(rho_r, rho_g)
      dvol = rho_r%pw_grid%dvol

      IF (iw > 0) THEN
         CALL m_flush(iw)
      END IF
      ! calculate the integrals
      stored_integrals = 0
      DO isp1 = 1, nspins
         CALL get_mo_set(mo_set=mos(isp1), nmo=nmo1)
         nmm = (nmo1*(nmo1 + 1))/2
         irange = get_irange_csr(nmm, eri_env%comm_exchange)
         DO i1 = 1, SIZE(orbitals, 1)
            iwa1 = orbitals(i1, isp1)
            IF (eri_env%eri_gpw%store_wfn) THEN
               wfn1 => wfn_a(iwa1, isp1)
            ELSE
               CALL calculate_wavefunction(fm_mo_coeff_as(isp1), iwa1, wfn1, rho_g, atomic_kind_set, &
                                           qs_kind_set, cell, dft_control, particle_set, pw_env_sub)
            END IF
            DO i2 = i1, SIZE(orbitals, 1)
               iwa2 = orbitals(i2, isp1)
               iwa12 = csr_idx_to_combined(iwa1, iwa2, nmo1)
               ! Skip calculation directly if the pair is not part of our subgroup
               IF (iwa12 < irange(1) .OR. iwa12 > irange(2)) CYCLE
               iwa12 = iwa12 - irange(1) + 1
               IF (eri_env%eri_gpw%store_wfn) THEN
                  wfn2 => wfn_a(iwa2, isp1)
               ELSE
                  CALL calculate_wavefunction(fm_mo_coeff_as(isp1), iwa2, wfn2, rho_g, atomic_kind_set, &
                                              qs_kind_set, cell, dft_control, particle_set, pw_env_sub)
               END IF
               ! calculate charge distribution and potential
               CALL pw_zero(rho_r)
               CALL pw_multiply(rho_r, wfn1, wfn2)
               CALL pw_transfer(rho_r, rho_g)
               CALL pw_poisson_solve(poisson_env, rho_g, pair_int, pot_g)

               ! screening using pair_int
               IF (pair_int < eri_env%eps_integral) CYCLE
               CALL pw_transfer(pot_g, rho_r)
               !
               IF (eri_env%method == eri_method_gpw_ht) THEN
                  CALL pw_scale(rho_r, dvol)
                  DO isp2 = isp1, nspins
                     CALL get_mo_set(mo_set=mos(isp2), nmo=nmo2)
                     nx = (nmo2*(nmo2 + 1))/2
                     ALLOCATE (eri(nx), eri_index(nx))
                     CALL dbcsr_set(mat_munu%matrix, 0.0_dp)
                     CALL integrate_v_rspace(rho_r, hmat=mat_munu, qs_env=qs_env, &
                                             calculate_forces=.FALSE., compute_tau=.FALSE., gapw=.FALSE., &
                                             pw_env_external=pw_env_sub, task_list_external=task_list_sub)

                     CALL dbcsr_multiply("N", "N", 1.0_dp, mat_munu%matrix, mo_coeff_as(isp2), &
                                         0.0_dp, matrix_pq_rnu(isp2), filter_eps=eri_env%eri_gpw%eps_filter)
                     CALL copy_dbcsr_to_fm(matrix_pq_rnu(isp2), fm_matrix_pq_rnu(isp2))

                     CALL cp_fm_get_info(fm_matrix_pq_rnu(isp2), ncol_global=ncol_global, nrow_global=nrow_global)

                     CALL parallel_gemm("T", "N", ncol_global, ncol_global, nrow_global, 0.5_dp, &
                                        fm_matrix_pq_rnu(isp2), fm_mo_coeff_as(isp2), &
                                        0.0_dp, fm_matrix_pq_rs(isp2))
                     CALL parallel_gemm("T", "N", ncol_global, ncol_global, nrow_global, 0.5_dp, &
                                        fm_mo_coeff_as(isp2), fm_matrix_pq_rnu(isp2), &
                                        1.0_dp, fm_matrix_pq_rs(isp2))

                     CALL cp_fm_get_info(fm_matrix_pq_rs(isp2), ncol_local=ncol_local, nrow_local=nrow_local, &
                                         col_indices=col_indices, row_indices=row_indices)

                     icount2 = 0
                     DO col_local = 1, ncol_local
                        iwb2 = orbitals(col_indices(col_local), isp2)
                        DO row_local = 1, nrow_local
                           iwb1 = orbitals(row_indices(row_local), isp2)

                           IF (iwb1 <= iwb2) THEN
                              iwb12 = csr_idx_to_combined(iwb1, iwb2, nmo2)
                              erint = fm_matrix_pq_rs(isp2)%local_data(row_local, col_local)
                              IF (ABS(erint) > eri_env%eps_integral .AND. (irange(1) - 1 + iwa12 <= iwb12 .OR. isp1 /= isp2)) THEN
                                 icount2 = icount2 + 1
                                 eri(icount2) = erint
                                 eri_index(icount2) = iwb12
                              END IF
                           END IF
                        END DO
                     END DO
                     stored_integrals = stored_integrals + icount2
                     !
                     isp = (isp1 - 1)*isp2 + (isp2 - isp1 + 1)
                     CALL update_csr_matrix(eri_env%eri(isp)%csr_mat, icount2, eri, eri_index, iwa12)
                     !
                     DEALLOCATE (eri, eri_index)
                  END DO
               ELSEIF (eri_env%method == eri_method_full_gpw) THEN
                  DO isp2 = isp1, nspins
                     CALL get_mo_set(mo_set=mos(isp2), nmo=nmo2)
                     nx = (nmo2*(nmo2 + 1))/2
                     ALLOCATE (eri(nx), eri_index(nx))
                     icount2 = 0
                     iwbs = 1
                     IF (isp1 == isp2) iwbs = i1
                     isp = (isp1 - 1)*isp2 + (isp2 - isp1 + 1)
                     DO i3 = iwbs, SIZE(orbitals, 1)
                        iwb1 = orbitals(i3, isp2)
                        IF (eri_env%eri_gpw%store_wfn) THEN
                           wfn3 => wfn_a(iwb1, isp2)
                        ELSE
                           CALL calculate_wavefunction(fm_mo_coeff_as(isp1), iwb1, wfn3, rho_g, atomic_kind_set, &
                                                       qs_kind_set, cell, dft_control, particle_set, pw_env_sub)
                        END IF
                        CALL pw_zero(wfn_r)
                        CALL pw_multiply(wfn_r, rho_r, wfn3)
                        iwbt = i3
                        IF (isp1 == isp2 .AND. i1 == i3) iwbt = i2
                        DO i4 = iwbt, SIZE(orbitals, 1)
                           iwb2 = orbitals(i4, isp2)
                           IF (eri_env%eri_gpw%store_wfn) THEN
                              wfn4 => wfn_a(iwb2, isp2)
                           ELSE
                              CALL calculate_wavefunction(fm_mo_coeff_as(isp1), iwb2, wfn4, rho_g, atomic_kind_set, &
                                                          qs_kind_set, cell, dft_control, particle_set, pw_env_sub)
                           END IF
                           ! We reduce the amount of communication by collecting the local sums first and sum globally later
                           erint = pw_integral_ab(wfn_r, wfn4, local_only=.TRUE.)
                           icount2 = icount2 + 1
                           eri(icount2) = erint
                           eri_index(icount2) = csr_idx_to_combined(iwb1, iwb2, nmo2)
                        END DO
                     END DO
                     ! Now, we sum the integrals globally
                     CALL eri_env%para_env_sub%sum(eri)
                     ! and we reorder the integrals to prevent storing too small integrals
                     intcount = 0
                     icount2 = 0
                     iwbs = 1
                     IF (isp1 == isp2) iwbs = i1
                     isp = (isp1 - 1)*isp2 + (isp2 - isp1 + 1)
                     DO i3 = iwbs, SIZE(orbitals, 1)
                        iwb1 = orbitals(i3, isp2)
                        iwbt = i3
                        IF (isp1 == isp2 .AND. i1 == i3) iwbt = i2
                        DO i4 = iwbt, SIZE(orbitals, 1)
                           iwb2 = orbitals(i4, isp2)
                           intcount = intcount + 1
                           erint = eri(intcount)
                           IF (ABS(erint) > eri_env%eps_integral) THEN
                              IF (MOD(intcount, eri_env%para_env_sub%num_pe) == eri_env%para_env_sub%mepos) THEN
                                 icount2 = icount2 + 1
                                 eri(icount2) = erint
                                 eri_index(icount2) = eri_index(intcount)
                              END IF
                           END IF
                        END DO
                     END DO
                     stored_integrals = stored_integrals + icount2
                     !
                     CALL update_csr_matrix(eri_env%eri(isp)%csr_mat, icount2, eri, eri_index, iwa12)
                     !
                     DEALLOCATE (eri, eri_index)
                  END DO
               ELSE
                  CPABORT("Unknown option")
               END IF
            END DO
         END DO
      END DO

      IF (print1 .AND. iw > 0) THEN
         WRITE (iw, "(T2,'ERI_GPW|',' Number of Integrals stored locally',T71,I10)") stored_integrals
      END IF

      IF (eri_env%eri_gpw%store_wfn) THEN
         DO ispin = 1, nspins
            DO i1 = 1, SIZE(orbitals, 1)
               iwfn = orbitals(i1, ispin)
               CALL wfn_a(iwfn, ispin)%release()
            END DO
         END DO
         DEALLOCATE (wfn_a)
      ELSE
         CALL wfn1%release()
         CALL wfn2%release()
         DEALLOCATE (wfn1, wfn2)
         IF (eri_env%method /= eri_method_gpw_ht) THEN
            CALL wfn3%release()
            CALL wfn4%release()
            DEALLOCATE (wfn3, wfn4)
         END IF
      END IF
      CALL auxbas_pw_pool%give_back_pw(wfn_r)
      CALL auxbas_pw_pool%give_back_pw(rho_g)
      CALL auxbas_pw_pool%give_back_pw(rho_r)
      CALL auxbas_pw_pool%give_back_pw(pot_g)

      IF (eri_env%method == eri_method_gpw_ht) THEN
         DO ispin = 1, nspins
            CALL dbcsr_release(mo_coeff_as(ispin))
            CALL dbcsr_release(matrix_pq_rnu(ispin))
            CALL cp_fm_release(fm_matrix_pq_rnu(ispin))
            CALL cp_fm_release(fm_matrix_pq_rs(ispin))
         END DO
         DEALLOCATE (matrix_pq_rnu, fm_matrix_pq_rnu, fm_matrix_pq_rs)
         CALL deallocate_task_list(task_list_sub)
      END IF
      DO ispin = 1, nspins
         CALL dbcsr_release(mo_coeff_as(ispin))
         CALL cp_fm_release(fm_mo_coeff_as(ispin))
      END DO
      DEALLOCATE (mo_coeff_as, fm_mo_coeff_as)
      CALL release_neighbor_list_sets(sab_orb_sub)
      CALL cp_blacs_env_release(blacs_env_sub)
      CALL dbcsr_release(mat_munu%matrix)
      DEALLOCATE (mat_munu%matrix)
      CALL pw_env_release(pw_env_sub)
      ! Return to the old qs_control
      dft_control%qs_control => qs_control_old
      DEALLOCATE (qs_control%e_cutoff)
      DEALLOCATE (qs_control)

      ! print out progress
      IF (iw > 0) THEN
         t2 = m_walltime()
         WRITE (iw, '(/,T2,A,T66,F14.2)') "ERI_GPW| ERI calculation took (sec)", t2 - t1
         CALL m_flush(iw)
      END IF

      CALL timestop(handle)

   END SUBROUTINE calculate_eri_gpw

! **************************************************************************************************
!> \brief Sets the Green's function for the ERI calculation. Here we deal with the G=0 case!
!> \param green ...
!> \param eri_env ...
!> \par History
!>      04.2016 created [JGH]
!>      08.2025 added support for the LR truncation [SB]
! **************************************************************************************************
   SUBROUTINE pw_eri_green_create(green, eri_env)

      TYPE(greens_fn_type), INTENT(INOUT)                :: green
      TYPE(eri_type)                                     :: eri_env

      COMPLEX(KIND=dp)                                   :: erf_fac_p, z_p
      INTEGER                                            :: ig
      REAL(KIND=dp)                                      :: cossin_fac, ea, erfcos_fac, exp_prefac, &
                                                            g, G0, g2, g3d, ga, Ginf, omega, &
                                                            omega2, Rc, Rc2

      ! initialize influence function
      ASSOCIATE (gf => green%influence_fn, grid => green%influence_fn%pw_grid)
         SELECT CASE (green%method)
         CASE (PERIODIC3D)

            SELECT CASE (eri_env%operator)
            CASE (eri_operator_coulomb)
               DO ig = grid%first_gne0, grid%ngpts_cut_local
                  g2 = grid%gsq(ig)
                  gf%array(ig) = fourpi/g2
               END DO
               IF (grid%have_g0) gf%array(1) = 0.0_dp

            CASE (eri_operator_yukawa)
               CALL cp_warn(__LOCATION__, "Yukawa operator has not been tested")
               omega2 = eri_env%omega**2
               DO ig = grid%first_gne0, grid%ngpts_cut_local
                  g2 = grid%gsq(ig)
                  gf%array(ig) = fourpi/(omega2 + g2)
               END DO
               IF (grid%have_g0) gf%array(1) = fourpi/omega2

            CASE (eri_operator_erf)
               omega2 = eri_env%omega**2
               DO ig = grid%first_gne0, grid%ngpts_cut_local
                  g2 = grid%gsq(ig)
                  gf%array(ig) = fourpi/g2*EXP(-0.25_dp*g2/omega2)
               END DO
               IF (grid%have_g0) gf%array(1) = 0.0_dp

            CASE (eri_operator_erfc)
               omega2 = eri_env%omega**2
               DO ig = grid%first_gne0, grid%ngpts_cut_local
                  g2 = grid%gsq(ig)
                  gf%array(ig) = fourpi/g2*(1.0_dp - EXP(-0.25_dp*g2/omega2))
               END DO
               IF (grid%have_g0) gf%array(1) = pi/omega2

            CASE (eri_operator_trunc)
               Rc = eri_env%cutoff_radius
               DO ig = grid%first_gne0, grid%ngpts_cut_local
                  g2 = grid%gsq(ig)
                  g = SQRT(g2)
                  ! Taylor expansion around zero
                  IF (g*Rc >= 0.005_dp) THEN
                     gf%array(ig) = fourpi/g2*(1.0_dp - COS(g*Rc))
                  ELSE
                     gf%array(ig) = fourpi/g2*(g*Rc)**2/2.0_dp*(1.0_dp - (g*Rc)**2/12.0_dp)
                  END IF
               END DO
               IF (grid%have_g0) gf%array(1) = twopi*Rc**2

            CASE (eri_operator_lr_trunc)
               omega = eri_env%omega
               omega2 = omega**2
               Rc = eri_env%cutoff_radius
               Rc2 = Rc**2
               G0 = 0.001_dp ! threshold for the G=0 case
               Ginf = 20.0_dp  ! threshold for the Taylor exapnsion arounf G=∞
               DO ig = grid%first_gne0, grid%ngpts_cut_local
                  g2 = grid%gsq(ig)
                  g = SQRT(g2)
                  IF (g <= 2.0_dp*G0) THEN
                     gf%array(ig) = -pi/omega2*erf(omega*Rc) &
                                    + twopi*Rc2*erf(omega*Rc) &
                                    + 2*rootpi*Rc*EXP(-omega2*Rc2)/omega
                  ELSE IF (g >= 2.0_dp*Ginf*omega) THEN
                     ! exponential prefactor
                     exp_prefac = EXP(-omega2*Rc2)/(rootpi*(omega2*Rc2 + 0.25_dp*g2/omega2))
                     ! cos sin factor
                     cossin_fac = omega*Rc*COS(g*Rc) - 0.5_dp*g/omega*SIN(g*Rc)
                     ! real erf term with cosine
                     erfcos_fac = ERF(omega*Rc)*COS(g*Rc)
                     ! Combine terms
                     gf%array(ig) = fourpi/g2*(-exp_prefac*cossin_fac - erfcos_fac)
                  ELSE
                     ! exponential prefactor
                     exp_prefac = twopi/g2*EXP(-0.25_dp*g2/omega2)
                     ! Compute complex arguments for erf
                     z_p = CMPLX(omega*Rc, 0.5_dp*g/omega, kind=dp)
                     ! Evaluate complex error functions
                     erf_fac_p = 2.0_dp*REAL(erfz_fast(z_p))
                     ! Real erf term with cosine
                     erfcos_fac = fourpi/g2*ERF(omega*Rc)*COS(g*Rc)
                     ! Combine terms
                     gf%array(ig) = exp_prefac*erf_fac_p - erfcos_fac
                  END IF
               END DO
               IF (grid%have_g0) THEN
                  gf%array(1) = -pi/omega2*ERF(omega*Rc) &
                                + twopi*Rc2*ERF(omega*Rc) &
                                + 2*rootpi*Rc*EXP(-omega2*Rc2)/omega
               END IF

            CASE DEFAULT
               CPABORT("Please specify a valid operator for the periodic Poisson solver")
            END SELECT

            ! The analytic Poisson solver simply limits the domain of integration
            ! of the Fourier transform to a sphere of radius Rc, rather than integrating
            ! over all space (-∞,∞)
         CASE (ANALYTIC0D)

            SELECT CASE (eri_env%operator)
               ! This is identical to the truncated Coulomb operator integrated
               ! over all space, when the truncation radius is equal to the radius of
               ! the Poisson solver
            CASE (eri_operator_coulomb, eri_operator_trunc)
               IF (eri_env%operator == eri_operator_coulomb) THEN
                  Rc = green%radius
               ELSE
                  Rc = eri_env%cutoff_radius
               END IF
               DO ig = grid%first_gne0, grid%ngpts_cut_local
                  g2 = grid%gsq(ig)
                  g = SQRT(g2)
                  ! Taylor expansion around zero
                  IF (g*Rc >= 0.005_dp) THEN
                     gf%array(ig) = fourpi/g2*(1.0_dp - COS(g*Rc))
                  ELSE
                     gf%array(ig) = fourpi/g2*(g*Rc)**2/2.0_dp*(1.0_dp - (g*Rc)**2/12.0_dp)
                  END IF
               END DO
               IF (grid%have_g0) gf%array(1) = twopi*Rc**2

               ! Not tested
            CASE (eri_operator_yukawa)
               CALL cp_warn(__LOCATION__, "Yukawa operator has not been tested")
               Rc = green%radius
               omega = eri_env%omega
               ea = EXP(-omega*Rc)
               DO ig = grid%first_gne0, grid%ngpts_cut_local
                  g2 = grid%gsq(ig)
                  g = SQRT(g2)
                  g3d = fourpi/(omega**2 + g2)
                  gf%array(ig) = g3d*(1.0_dp - ea*(COS(g*Rc) + omega/g*SIN(g*Rc)))
               END DO
               IF (grid%have_g0) gf%array(1) = fourpi/(omega**2)*(1.0_dp - ea*(1.0_dp + omega*Rc))

               ! Long-range Coulomb
               ! TODO: this should be equivalent to LR truncated Coulomb from above!
            CASE (eri_operator_erf, eri_operator_lr_trunc)
               IF (eri_env%operator == eri_operator_erf) THEN
                  Rc = green%radius
               ELSE
                  Rc = eri_env%cutoff_radius
               END IF
               omega2 = eri_env%omega**2
               DO ig = grid%first_gne0, grid%ngpts_cut_local
                  g2 = grid%gsq(ig)
                  g = SQRT(g2)
                  ga = -0.25_dp*g2/omega2
                  gf%array(ig) = fourpi/g2*EXP(ga)*(1.0_dp - COS(g*Rc))
               END DO
               IF (grid%have_g0) gf%array(1) = twopi*Rc**2

               ! Short-range Coulomb
               ! TODO: this should actually be properly derived and see whether it is correct
            CASE (eri_operator_erfc)
               CALL cp_warn(__LOCATION__, &
                            "Short-range Coulomb operator may be incorrect with ANALYTIC0D Poisson solver")
               Rc = green%radius
               omega2 = eri_env%omega**2
               DO ig = grid%first_gne0, grid%ngpts_cut_local
                  g2 = grid%gsq(ig)
                  g = SQRT(g2)
                  ga = -0.25_dp*g2/omega2
                  gf%array(ig) = fourpi/g2*(1.0_dp - EXP(ga))*(1.0_dp - COS(g*Rc))
               END DO
               IF (grid%have_g0) gf%array(1) = pi/omega2

            CASE DEFAULT
               CPABORT("Unsupported operator")
            END SELECT

         CASE DEFAULT
            CPABORT("Unsupported Poisson solver")
         END SELECT
      END ASSOCIATE

   END SUBROUTINE pw_eri_green_create

! **************************************************************************************************
!> \brief Adds data for a new row to the csr matrix
!> \param csr_mat ...
!> \param nnz ...
!> \param rdat ...
!> \param rind ...
!> \param irow ...
!> \par History
!>      04.2016 created [JGH]
! **************************************************************************************************
   SUBROUTINE update_csr_matrix(csr_mat, nnz, rdat, rind, irow)

      TYPE(dbcsr_csr_type), INTENT(INOUT)                :: csr_mat
      INTEGER, INTENT(IN)                                :: nnz
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rdat
      INTEGER, DIMENSION(:), INTENT(IN)                  :: rind
      INTEGER, INTENT(IN)                                :: irow

      INTEGER                                            :: k, nrow, nze, nze_new

      IF (irow /= 0) THEN
         nze = csr_mat%nze_local
         nze_new = nze + nnz
         ! values
         CALL reallocate(csr_mat%nzval_local%r_dp, 1, nze_new)
         csr_mat%nzval_local%r_dp(nze + 1:nze_new) = rdat(1:nnz)
         ! col indices
         CALL reallocate(csr_mat%colind_local, 1, nze_new)
         csr_mat%colind_local(nze + 1:nze_new) = rind(1:nnz)
         ! rows
         nrow = csr_mat%nrows_local
         CALL reallocate(csr_mat%rowptr_local, 1, irow + 1)
         csr_mat%rowptr_local(nrow + 1:irow) = nze + 1
         csr_mat%rowptr_local(irow + 1) = nze_new + 1
         ! nzerow
         CALL reallocate(csr_mat%nzerow_local, 1, irow)
         DO k = nrow + 1, irow
            csr_mat%nzerow_local(k) = csr_mat%rowptr_local(k + 1) - csr_mat%rowptr_local(k)
         END DO
         csr_mat%nrows_local = irow
         csr_mat%nze_local = csr_mat%nze_local + nnz
      END IF
      csr_mat%nze_total = csr_mat%nze_total + nnz
      csr_mat%has_indices = .TRUE.

   END SUBROUTINE update_csr_matrix

! **************************************************************************************************
!> \brief Computes and prints the active orbitals on Cube Files
!> \param input ...
!> \param qs_env the qs_env in which the qs_env lives
!> \param mos ...
! **************************************************************************************************
   SUBROUTINE print_orbital_cubes(input, qs_env, mos)
      TYPE(section_vals_type), POINTER                   :: input
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mo_set_type), DIMENSION(:), INTENT(IN)        :: mos

      CHARACTER(LEN=default_path_length)                 :: filebody, filename, title
      INTEGER                                            :: i, imo, isp, nmo, str(3), unit_nr
      INTEGER, DIMENSION(:), POINTER                     :: alist, blist, istride
      LOGICAL                                            :: do_mo, explicit_a, explicit_b
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_c1d_gs_type)                               :: wf_g
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_r3d_rs_type)                               :: wf_r
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_subsys_type), POINTER                      :: subsys
      TYPE(section_vals_type), POINTER                   :: dft_section, scf_input

      CALL section_vals_val_get(input, "FILENAME", c_val=filebody)
      CALL section_vals_val_get(input, "STRIDE", i_vals=istride)
      IF (SIZE(istride) == 1) THEN
         str(1:3) = istride(1)
      ELSEIF (SIZE(istride) == 3) THEN
         str(1:3) = istride(1:3)
      ELSE
         CPABORT("STRIDE arguments inconsistent")
      END IF
      CALL section_vals_val_get(input, "ALIST", i_vals=alist, explicit=explicit_a)
      CALL section_vals_val_get(input, "BLIST", i_vals=blist, explicit=explicit_b)

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      para_env=para_env, &
                      subsys=subsys, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      cell=cell, &
                      particle_set=particle_set, &
                      pw_env=pw_env, &
                      input=scf_input)

      CALL qs_subsys_get(subsys, particles=particles)
      !
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
      CALL auxbas_pw_pool%create_pw(wf_r)
      CALL auxbas_pw_pool%create_pw(wf_g)
      !
      dft_section => section_vals_get_subs_vals(scf_input, "DFT")
      !
      DO isp = 1, SIZE(mos)
         CALL get_mo_set(mo_set=mos(isp), mo_coeff=mo_coeff, nmo=nmo)

         IF (SIZE(mos) > 1) THEN
            SELECT CASE (isp)
            CASE (1)
               CALL write_mo_set_to_output_unit(mos(isp), qs_kind_set, particle_set, &
                                                dft_section, 4, 0, final_mos=.TRUE., spin="ALPHA")
            CASE (2)
               CALL write_mo_set_to_output_unit(mos(isp), qs_kind_set, particle_set, &
                                                dft_section, 4, 0, final_mos=.TRUE., spin="BETA")
            CASE DEFAULT
               CPABORT("Invalid spin")
            END SELECT
         ELSE
            CALL write_mo_set_to_output_unit(mos(isp), qs_kind_set, particle_set, &
                                             dft_section, 4, 0, final_mos=.TRUE.)
         END IF

         DO imo = 1, nmo
            IF (isp == 1 .AND. explicit_a) THEN
               IF (alist(1) == -1) THEN
                  do_mo = .TRUE.
               ELSE
                  do_mo = .FALSE.
                  DO i = 1, SIZE(alist)
                     IF (imo == alist(i)) do_mo = .TRUE.
                  END DO
               END IF
            ELSE IF (isp == 2 .AND. explicit_b) THEN
               IF (blist(1) == -1) THEN
                  do_mo = .TRUE.
               ELSE
                  do_mo = .FALSE.
                  DO i = 1, SIZE(blist)
                     IF (imo == blist(i)) do_mo = .TRUE.
                  END DO
               END IF
            ELSE
               do_mo = .TRUE.
            END IF
            IF (.NOT. do_mo) CYCLE
            CALL calculate_wavefunction(mo_coeff, imo, wf_r, wf_g, atomic_kind_set, &
                                        qs_kind_set, cell, dft_control, particle_set, pw_env)
            IF (para_env%is_source()) THEN
               WRITE (filename, '(A,A1,I4.4,A1,I1.1,A)') TRIM(filebody), "_", imo, "_", isp, ".cube"
               CALL open_file(filename, unit_number=unit_nr, file_status="UNKNOWN", file_action="WRITE")
               WRITE (title, *) "Active Orbital ", imo, " spin ", isp
            ELSE
               unit_nr = -1
            END IF
            CALL cp_pw_to_cube(wf_r, unit_nr, title, particles=particles, stride=istride)
            IF (para_env%is_source()) THEN
               CALL close_file(unit_nr)
            END IF
         END DO
      END DO

      CALL auxbas_pw_pool%give_back_pw(wf_r)
      CALL auxbas_pw_pool%give_back_pw(wf_g)

   END SUBROUTINE print_orbital_cubes

! **************************************************************************************************
!> \brief Writes a FCIDUMP file
!> \param active_space_env ...
!> \param as_input ...
!> \param restricted ...
!> \par History
!>      04.2016 created [JGH]
! **************************************************************************************************
   SUBROUTINE fcidump(active_space_env, as_input, restricted)

      TYPE(active_space_type), POINTER                   :: active_space_env
      TYPE(section_vals_type), POINTER                   :: as_input
      LOGICAL, INTENT(IN)                                :: restricted

      INTEGER                                            :: i, i1, i2, i3, i4, isym, iw, m1, m2, &
                                                            nmo, norb, nspins
      REAL(KIND=dp)                                      :: checksum, esub
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: fmat
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(eri_fcidump_checksum)                         :: eri_checksum

      checksum = 0.0_dp

      logger => cp_get_default_logger()
      iw = cp_print_key_unit_nr(logger, as_input, "FCIDUMP", &
                                extension=".fcidump", file_status="REPLACE", file_action="WRITE", file_form="FORMATTED")
      !
      nspins = active_space_env%nspins
      norb = SIZE(active_space_env%active_orbitals, 1)
      IF (nspins == 1 .OR. restricted) THEN
         ! Closed shell or restricted open-shell
         ASSOCIATE (ms2 => active_space_env%multiplicity, &
                    nelec => active_space_env%nelec_active)

            IF (iw > 0) THEN
               WRITE (iw, "(A,A,I4,A,I4,A,I2,A)") "&FCI", " NORB=", norb, ",NELEC=", nelec, ",MS2=", ms2, ","
               isym = 1
               WRITE (iw, "(A,1000(I1,','))") "  ORBSYM=", (isym, i=1, norb)
               isym = 0
               WRITE (iw, "(A,I1,A)") "  ISYM=", isym, ","
               IF (restricted) WRITE (iw, "(A,I1,A)") "  UHF=", 0, ","
               WRITE (iw, "(A)") " /"
            END IF
            !
            ! Print integrals: ERI
            CALL active_space_env%eri%eri_foreach(1, active_space_env%active_orbitals, &
                                                  eri_fcidump_print(iw, 1, 1), 1, 1)
            CALL eri_checksum%set(1, 1)
            CALL active_space_env%eri%eri_foreach(1, active_space_env%active_orbitals, eri_checksum, 1, 1)

            ! Print integrals: Fij
            ! replicate Fock matrix
            nmo = active_space_env%eri%norb
            ALLOCATE (fmat(nmo, nmo))
            CALL replicate_and_symmetrize_matrix(nmo, active_space_env%fock_sub(1), fmat)
            IF (iw > 0) THEN
               i3 = 0; i4 = 0
               DO m1 = 1, SIZE(active_space_env%active_orbitals, 1)
                  i1 = active_space_env%active_orbitals(m1, 1)
                  DO m2 = m1, SIZE(active_space_env%active_orbitals, 1)
                     i2 = active_space_env%active_orbitals(m2, 1)
                     checksum = checksum + ABS(fmat(i1, i2))
                     WRITE (iw, "(ES23.16,4I4)") fmat(i1, i2), m1, m2, i3, i4
                  END DO
               END DO
            END IF
            DEALLOCATE (fmat)
            ! Print energy
            esub = active_space_env%energy_inactive
            i1 = 0; i2 = 0; i3 = 0; i4 = 0
            checksum = checksum + ABS(esub)
            IF (iw > 0) WRITE (iw, "(ES23.16,4I4)") esub, i1, i2, i3, i4
         END ASSOCIATE

      ELSE
         ASSOCIATE (ms2 => active_space_env%multiplicity, &
                    nelec => active_space_env%nelec_active)

            IF (iw > 0) THEN
               WRITE (iw, "(A,A,I4,A,I4,A,I2,A)") "&FCI", " NORB=", norb, ",NELEC=", nelec, ",MS2=", ms2, ","
               isym = 1
               WRITE (iw, "(A,1000(I1,','))") "  ORBSYM=", (isym, i=1, norb)
               isym = 0
               WRITE (iw, "(A,I1,A)") "  ISYM=", isym, ","
               WRITE (iw, "(A,I1,A)") "  UHF=", 1, ","
               WRITE (iw, "(A)") " /"
            END IF
            !
            ! Print integrals: ERI
            ! alpha-alpha
            CALL active_space_env%eri%eri_foreach(1, active_space_env%active_orbitals, &
                                                  eri_fcidump_print(iw, 1, 1), 1, 1)
            CALL eri_checksum%set(1, 1)
            CALL active_space_env%eri%eri_foreach(1, active_space_env%active_orbitals, eri_checksum, 1, 1)
            ! alpha-beta
            CALL active_space_env%eri%eri_foreach(2, active_space_env%active_orbitals, &
                                                  eri_fcidump_print(iw, 1, norb + 1), 1, 2)
            CALL eri_checksum%set(1, norb + 1)
            CALL active_space_env%eri%eri_foreach(2, active_space_env%active_orbitals, eri_checksum, 1, 2)
            ! beta-beta
            CALL active_space_env%eri%eri_foreach(3, active_space_env%active_orbitals, &
                                                  eri_fcidump_print(iw, norb + 1, norb + 1), 2, 2)
            CALL eri_checksum%set(norb + 1, norb + 1)
            CALL active_space_env%eri%eri_foreach(3, active_space_env%active_orbitals, eri_checksum, 2, 2)
            ! Print integrals: Fij
            ! alpha
            nmo = active_space_env%eri%norb
            ALLOCATE (fmat(nmo, nmo))
            CALL replicate_and_symmetrize_matrix(nmo, active_space_env%fock_sub(1), fmat)
            IF (iw > 0) THEN
               i3 = 0; i4 = 0
               DO m1 = 1, norb
                  i1 = active_space_env%active_orbitals(m1, 1)
                  DO m2 = m1, norb
                     i2 = active_space_env%active_orbitals(m2, 1)
                     checksum = checksum + ABS(fmat(i1, i2))
                     WRITE (iw, "(ES23.16,4I4)") fmat(i1, i2), m1, m2, i3, i4
                  END DO
               END DO
            END IF
            DEALLOCATE (fmat)
            ! beta
            ALLOCATE (fmat(nmo, nmo))
            CALL replicate_and_symmetrize_matrix(nmo, active_space_env%fock_sub(2), fmat)
            IF (iw > 0) THEN
               i3 = 0; i4 = 0
               DO m1 = 1, SIZE(active_space_env%active_orbitals, 1)
                  i1 = active_space_env%active_orbitals(m1, 2)
                  DO m2 = m1, SIZE(active_space_env%active_orbitals, 1)
                     i2 = active_space_env%active_orbitals(m2, 2)
                     checksum = checksum + ABS(fmat(i1, i2))
                     WRITE (iw, "(ES23.16,4I4)") fmat(i1, i2), m1 + norb, m2 + norb, i3, i4
                  END DO
               END DO
            END IF
            DEALLOCATE (fmat)
            ! Print energy
            esub = active_space_env%energy_inactive
            i1 = 0; i2 = 0; i3 = 0; i4 = 0
            checksum = checksum + ABS(esub)
            IF (iw > 0) WRITE (iw, "(ES23.16,4I4)") esub, i1, i2, i3, i4
         END ASSOCIATE
      END IF
      !
      CALL cp_print_key_finished_output(iw, logger, as_input, "FCIDUMP")

      !>>
      iw = cp_logger_get_default_io_unit(logger)
      IF (iw > 0) WRITE (iw, '(T4,A,T66,F12.8)') "FCIDUMP| Checksum:", eri_checksum%checksum + checksum
      !<<

   END SUBROUTINE fcidump

! **************************************************************************************************
!> \brief replicate and symmetrize a matrix
!> \param norb the number of orbitals
!> \param distributed_matrix ...
!> \param replicated_matrix ...
! **************************************************************************************************
   SUBROUTINE replicate_and_symmetrize_matrix(norb, distributed_matrix, replicated_matrix)
      INTEGER, INTENT(IN)                                :: norb
      TYPE(cp_fm_type), INTENT(IN)                       :: distributed_matrix
      REAL(dp), DIMENSION(:, :), INTENT(INOUT)           :: replicated_matrix

      INTEGER                                            :: i1, i2
      REAL(dp)                                           :: mval

      replicated_matrix(:, :) = 0.0_dp
      DO i1 = 1, norb
         DO i2 = i1, norb
            CALL cp_fm_get_element(distributed_matrix, i1, i2, mval)
            replicated_matrix(i1, i2) = mval
            replicated_matrix(i2, i1) = mval
         END DO
      END DO
   END SUBROUTINE replicate_and_symmetrize_matrix

! **************************************************************************************************
!> \brief Calculates active space Fock matrix and inactive energy
!> \param active_space_env ...
!> \param restricted ...
!> \par History
!>      06.2016 created [JGH]
! **************************************************************************************************
   SUBROUTINE subspace_fock_matrix(active_space_env, restricted)

      TYPE(active_space_type), POINTER                   :: active_space_env
      LOGICAL, INTENT(IN)                                :: restricted

      INTEGER                                            :: i1, i2, is, norb, nspins
      REAL(KIND=dp)                                      :: eeri, eref, esub, mval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: ks_a_mat, ks_a_ref, ks_b_mat, ks_b_ref, &
                                                            ks_mat, ks_ref, p_a_mat, p_b_mat, p_mat
      TYPE(cp_fm_type), POINTER                          :: matrix, mo_coef
      TYPE(dbcsr_csr_type), POINTER                      :: eri, eri_aa, eri_ab, eri_bb

      eref = active_space_env%energy_ref
      nspins = active_space_env%nspins

      IF (nspins == 1) THEN
         CALL get_mo_set(active_space_env%mos_active(1), nmo=norb, mo_coeff=mo_coef)
         !
         ! Loop over ERI, calculate subspace HF energy and Fock matrix
         !
         ! replicate KS, Core, and P matrices
         ALLOCATE (ks_mat(norb, norb), ks_ref(norb, norb), p_mat(norb, norb))
         ks_ref = 0.0_dp

         ! ks_mat contains the KS/Fock matrix (of full density) projected onto the AS MO subspace (f_ref in eq. 19)
         CALL replicate_and_symmetrize_matrix(norb, active_space_env%ks_sub(1), ks_mat)
         CALL replicate_and_symmetrize_matrix(norb, active_space_env%p_active(1), p_mat)

         ! compute ks_ref = V_H[rho^A] + V_HFX[rho^A]
         eri => active_space_env%eri%eri(1)%csr_mat
         CALL build_subspace_fock_matrix(active_space_env%active_orbitals, eri, p_mat, ks_ref, &
                                         active_space_env%eri%comm_exchange)

         ! compute eeri = E_H[rho^A] + E_HFX[rho^A] as
         ! eeri = 1/2 * (SUM_pq (V_H[rho^A] + V_HFX[rho^A])_pq * D^A_pq)
         eeri = 0.5_dp*SUM(ks_ref*p_mat)

         ! now calculate the inactive energy acoording to eq. 19, that is
         ! esub = E^I = E_ref - f_ref .* D^A + E_H[rho^A] + E_HFX[rho^A]
         ! where f^ref = ks_mat, which is the KS/Fock matrix in MO basis, transformed previously
         ! and is equal to ks_mat = h^0 + V_core + V_H[rho] + V_HFX[rho]
         esub = eref - SUM(ks_mat(1:norb, 1:norb)*p_mat(1:norb, 1:norb)) + eeri

         ! reuse ks_mat to store f^I = f^ref - (V_H[rho^A] + V_HFX[rho^A]) according to eq. 20
         ks_mat(1:norb, 1:norb) = ks_mat(1:norb, 1:norb) - ks_ref(1:norb, 1:norb)
         ! this is now the embedding potential for the AS calculation!

         active_space_env%energy_inactive = esub

         CALL cp_fm_release(active_space_env%fock_sub)
         ALLOCATE (active_space_env%fock_sub(nspins))
         DO is = 1, nspins
            matrix => active_space_env%ks_sub(is)
            CALL cp_fm_create(active_space_env%fock_sub(is), matrix%matrix_struct, &
                              name="Active Fock operator")
         END DO
         matrix => active_space_env%fock_sub(1)
         DO i1 = 1, norb
            DO i2 = 1, norb
               mval = ks_mat(i1, i2)
               CALL cp_fm_set_element(matrix, i1, i2, mval)
            END DO
         END DO
      ELSE

         CALL get_mo_set(active_space_env%mos_active(1), nmo=norb)
         !
         ! Loop over ERI, calculate subspace HF energy and Fock matrix
         !
         ! replicate KS, Core, and P matrices
         ALLOCATE (ks_a_mat(norb, norb), ks_b_mat(norb, norb), &
              &    ks_a_ref(norb, norb), ks_b_ref(norb, norb), &
              &     p_a_mat(norb, norb), p_b_mat(norb, norb))
         ks_a_ref(:, :) = 0.0_dp; ks_b_ref(:, :) = 0.0_dp

         CALL replicate_and_symmetrize_matrix(norb, active_space_env%p_active(1), p_a_mat)
         CALL replicate_and_symmetrize_matrix(norb, active_space_env%p_active(2), p_b_mat)
         CALL replicate_and_symmetrize_matrix(norb, active_space_env%ks_sub(1), ks_a_mat)
         CALL replicate_and_symmetrize_matrix(norb, active_space_env%ks_sub(2), ks_b_mat)
         !
         !
         IF (restricted) THEN
            ! In the restricted case, we use the same ERIs for each spin channel
            eri_aa => active_space_env%eri%eri(1)%csr_mat
            CALL build_subspace_spin_fock_matrix(active_space_env%active_orbitals, eri_aa, eri_aa, p_a_mat, p_b_mat, ks_a_ref, &
                                                 tr_mixed_eri=.FALSE., comm_exchange=active_space_env%eri%comm_exchange)
            CALL build_subspace_spin_fock_matrix(active_space_env%active_orbitals, eri_aa, eri_aa, p_b_mat, p_a_mat, ks_b_ref, &
                                                 tr_mixed_eri=.TRUE., comm_exchange=active_space_env%eri%comm_exchange)
         ELSE
            eri_aa => active_space_env%eri%eri(1)%csr_mat
            eri_ab => active_space_env%eri%eri(2)%csr_mat
            eri_bb => active_space_env%eri%eri(3)%csr_mat
            CALL build_subspace_spin_fock_matrix(active_space_env%active_orbitals, eri_aa, eri_ab, p_a_mat, p_b_mat, ks_a_ref, &
                                                 tr_mixed_eri=.FALSE., comm_exchange=active_space_env%eri%comm_exchange)
            CALL build_subspace_spin_fock_matrix(active_space_env%active_orbitals, eri_bb, eri_ab, p_b_mat, p_a_mat, ks_b_ref, &
                                                 tr_mixed_eri=.TRUE., comm_exchange=active_space_env%eri%comm_exchange)
         END IF
         !
         ! calculate energy
         eeri = 0.0_dp
         eeri = 0.5_dp*(SUM(ks_a_ref*p_a_mat) + SUM(ks_b_ref*p_b_mat))
         esub = eref - SUM(ks_a_mat*p_a_mat) - SUM(ks_b_mat*p_b_mat) + eeri
         ks_a_mat(:, :) = ks_a_mat(:, :) - ks_a_ref(:, :)
         ks_b_mat(:, :) = ks_b_mat(:, :) - ks_b_ref(:, :)
         !
         active_space_env%energy_inactive = esub
         !
         CALL cp_fm_release(active_space_env%fock_sub)
         ALLOCATE (active_space_env%fock_sub(nspins))
         DO is = 1, nspins
            matrix => active_space_env%ks_sub(is)
            CALL cp_fm_create(active_space_env%fock_sub(is), matrix%matrix_struct, &
                              name="Active Fock operator")
         END DO

         matrix => active_space_env%fock_sub(1)
         DO i1 = 1, norb
            DO i2 = 1, norb
               mval = ks_a_mat(i1, i2)
               CALL cp_fm_set_element(matrix, i1, i2, mval)
            END DO
         END DO
         matrix => active_space_env%fock_sub(2)
         DO i1 = 1, norb
            DO i2 = 1, norb
               mval = ks_b_mat(i1, i2)
               CALL cp_fm_set_element(matrix, i1, i2, mval)
            END DO
         END DO

      END IF

   END SUBROUTINE subspace_fock_matrix

! **************************************************************************************************
!> \brief build subspace fockian
!> \param active_orbitals the active orbital indices
!> \param eri two electon integrals in MO
!> \param p_mat density matrix
!> \param ks_ref fockian matrix
!> \param comm_exchange ...
! **************************************************************************************************
   SUBROUTINE build_subspace_fock_matrix(active_orbitals, eri, p_mat, ks_ref, comm_exchange)
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: active_orbitals
      TYPE(dbcsr_csr_type), INTENT(IN)                   :: eri
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: p_mat
      REAL(dp), DIMENSION(:, :), INTENT(INOUT)           :: ks_ref
      TYPE(mp_comm_type), INTENT(IN)                     :: comm_exchange

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_subspace_fock_matrix'

      INTEGER                                            :: handle, i1, i12, i12l, i2, i3, i34, &
                                                            i34l, i4, irptr, m1, m2, nindex, &
                                                            nmo_total, norb
      INTEGER, DIMENSION(2)                              :: irange
      REAL(dp)                                           :: erint
      TYPE(mp_comm_type)                                 :: mp_group

      CALL timeset(routineN, handle)

      ! Nothing to do
      norb = SIZE(active_orbitals, 1)
      nmo_total = SIZE(p_mat, 1)
      nindex = (nmo_total*(nmo_total + 1))/2
      CALL mp_group%set_handle(eri%mp_group%get_handle())
      irange = get_irange_csr(nindex, comm_exchange)
      DO m1 = 1, norb
         i1 = active_orbitals(m1, 1)
         DO m2 = m1, norb
            i2 = active_orbitals(m2, 1)
            i12 = csr_idx_to_combined(i1, i2, nmo_total)
            IF (i12 >= irange(1) .AND. i12 <= irange(2)) THEN
               i12l = i12 - irange(1) + 1
               irptr = eri%rowptr_local(i12l) - 1
               DO i34l = 1, eri%nzerow_local(i12l)
                  i34 = eri%colind_local(irptr + i34l)
                  CALL csr_idx_from_combined(i34, nmo_total, i3, i4)
                  erint = eri%nzval_local%r_dp(irptr + i34l)
                  ! Coulomb
                  ks_ref(i1, i2) = ks_ref(i1, i2) + erint*p_mat(i3, i4)
                  IF (i3 /= i4) THEN
                     ks_ref(i1, i2) = ks_ref(i1, i2) + erint*p_mat(i3, i4)
                  END IF
                  IF (i12 /= i34) THEN
                     ks_ref(i3, i4) = ks_ref(i3, i4) + erint*p_mat(i1, i2)
                     IF (i1 /= i2) THEN
                        ks_ref(i3, i4) = ks_ref(i3, i4) + erint*p_mat(i1, i2)
                     END IF
                  END IF
                  ! Exchange
                  erint = -0.5_dp*erint
                  ks_ref(i1, i3) = ks_ref(i1, i3) + erint*p_mat(i2, i4)
                  IF (i1 /= i2) THEN
                     ks_ref(i2, i3) = ks_ref(i2, i3) + erint*p_mat(i1, i4)
                  END IF
                  IF (i3 /= i4) THEN
                     ks_ref(i1, i4) = ks_ref(i1, i4) + erint*p_mat(i2, i3)
                  END IF
                  IF (i1 /= i2 .AND. i3 /= i4) THEN
                     ks_ref(i2, i4) = ks_ref(i2, i4) + erint*p_mat(i1, i3)
                  END IF
               END DO
            END IF
         END DO
      END DO
      !
      DO m1 = 1, norb
         i1 = active_orbitals(m1, 1)
         DO m2 = m1, norb
            i2 = active_orbitals(m2, 1)
            ks_ref(i2, i1) = ks_ref(i1, i2)
         END DO
      END DO
      CALL mp_group%sum(ks_ref)

      CALL timestop(handle)

   END SUBROUTINE build_subspace_fock_matrix

! **************************************************************************************************
!> \brief build subspace fockian for unrestricted spins
!> \param active_orbitals the active orbital indices
!> \param eri_aa two electon integrals in MO with parallel spins
!> \param eri_ab two electon integrals in MO with anti-parallel spins
!> \param p_a_mat density matrix for up-spin
!> \param p_b_mat density matrix for down-spin
!> \param ks_a_ref fockian matrix for up-spin
!> \param tr_mixed_eri boolean to indicate Coulomb interaction alignment
!> \param comm_exchange ...
! **************************************************************************************************
   SUBROUTINE build_subspace_spin_fock_matrix(active_orbitals, eri_aa, eri_ab, p_a_mat, p_b_mat, ks_a_ref, tr_mixed_eri, &
                                              comm_exchange)
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: active_orbitals
      TYPE(dbcsr_csr_type), INTENT(IN)                   :: eri_aa, eri_ab
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: p_a_mat, p_b_mat
      REAL(dp), DIMENSION(:, :), INTENT(INOUT)           :: ks_a_ref
      LOGICAL, INTENT(IN)                                :: tr_mixed_eri
      TYPE(mp_comm_type), INTENT(IN)                     :: comm_exchange

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_subspace_spin_fock_matrix'

      INTEGER                                            :: handle, i1, i12, i12l, i2, i3, i34, &
                                                            i34l, i4, irptr, m1, m2, nindex, &
                                                            nmo_total, norb, spin1, spin2
      INTEGER, DIMENSION(2)                              :: irange
      REAL(dp)                                           :: erint
      TYPE(mp_comm_type)                                 :: mp_group

      CALL timeset(routineN, handle)

      norb = SIZE(active_orbitals, 1)
      nmo_total = SIZE(p_a_mat, 1)
      nindex = (nmo_total*(nmo_total + 1))/2
      irange = get_irange_csr(nindex, comm_exchange)
      IF (tr_mixed_eri) THEN
         spin1 = 2
         spin2 = 1
      ELSE
         spin1 = 1
         spin2 = 2
      END IF
      DO m1 = 1, norb
         i1 = active_orbitals(m1, spin1)
         DO m2 = m1, norb
            i2 = active_orbitals(m2, spin1)
            i12 = csr_idx_to_combined(i1, i2, nmo_total)
            IF (i12 >= irange(1) .AND. i12 <= irange(2)) THEN
               i12l = i12 - irange(1) + 1
               irptr = eri_aa%rowptr_local(i12l) - 1
               DO i34l = 1, eri_aa%nzerow_local(i12l)
                  i34 = eri_aa%colind_local(irptr + i34l)
                  CALL csr_idx_from_combined(i34, nmo_total, i3, i4)
                  erint = eri_aa%nzval_local%r_dp(irptr + i34l)
                  ! Coulomb
                  !F_ij += (ij|kl)*d_kl
                  ks_a_ref(i1, i2) = ks_a_ref(i1, i2) + erint*p_a_mat(i3, i4)
                  IF (i12 /= i34) THEN
                     !F_kl += (ij|kl)*d_ij
                     ks_a_ref(i3, i4) = ks_a_ref(i3, i4) + erint*p_a_mat(i1, i2)
                  END IF
                  ! Exchange
                  erint = -1.0_dp*erint
                  !F_ik -= (ij|kl)*d_jl
                  ks_a_ref(i1, i3) = ks_a_ref(i1, i3) + erint*p_a_mat(i2, i4)
                  IF (i1 /= i2) THEN
                     !F_jk -= (ij|kl)*d_il
                     ks_a_ref(i2, i3) = ks_a_ref(i2, i3) + erint*p_a_mat(i1, i4)
                  END IF
                  IF (i3 /= i4) THEN
                     !F_il -= (ij|kl)*d_jk
                     ks_a_ref(i1, i4) = ks_a_ref(i1, i4) + erint*p_a_mat(i2, i3)
                  END IF
                  IF (i1 /= i2 .AND. i3 /= i4) THEN
                     !F_jl -= (ij|kl)*d_ik
                     ks_a_ref(i2, i4) = ks_a_ref(i2, i4) + erint*p_a_mat(i1, i3)
                  END IF
               END DO
            END IF
         END DO
      END DO
      !

      irange = get_irange_csr(nindex, comm_exchange)
      DO m1 = 1, norb
         i1 = active_orbitals(m1, 1)
         DO m2 = m1, norb
            i2 = active_orbitals(m2, 1)
            i12 = csr_idx_to_combined(i1, i2, nmo_total)
            IF (i12 >= irange(1) .AND. i12 <= irange(2)) THEN
               i12l = i12 - irange(1) + 1
               irptr = eri_ab%rowptr_local(i12l) - 1
               DO i34l = 1, eri_ab%nzerow_local(i12l)
                  i34 = eri_ab%colind_local(irptr + i34l)
                  CALL csr_idx_from_combined(i34, nmo_total, i3, i4)
                  erint = eri_ab%nzval_local%r_dp(irptr + i34l)
                  ! Coulomb
                  IF (tr_mixed_eri) THEN
                     !F_kl += (kl beta|ij alpha )*d_alpha_ij
                     ks_a_ref(i3, i4) = ks_a_ref(i3, i4) + erint*p_b_mat(i1, i2)
                  ELSE
                     !F_ij += (ij alpha|kl beta )*d_beta_kl
                     ks_a_ref(i1, i2) = ks_a_ref(i1, i2) + erint*p_b_mat(i3, i4)
                  END IF
               END DO
            END IF
         END DO
      END DO
      !
      DO m1 = 1, norb
         i1 = active_orbitals(m1, spin1)
         DO m2 = m1, norb
            i2 = active_orbitals(m2, spin1)
            ks_a_ref(i2, i1) = ks_a_ref(i1, i2)
         END DO
      END DO
      CALL mp_group%set_handle(eri_aa%mp_group%get_handle())
      CALL mp_group%sum(ks_a_ref)

      CALL timestop(handle)

   END SUBROUTINE build_subspace_spin_fock_matrix

! **************************************************************************************************
!> \brief Creates a local basis
!> \param pro_basis_set ...
!> \param zval ...
!> \param ishell ...
!> \param nshell ...
!> \param lnam ...
!> \par History
!>      05.2016 created [JGH]
! **************************************************************************************************
   SUBROUTINE create_pro_basis(pro_basis_set, zval, ishell, nshell, lnam)
      TYPE(gto_basis_set_type), POINTER                  :: pro_basis_set
      INTEGER, INTENT(IN)                                :: zval, ishell
      INTEGER, DIMENSION(:), INTENT(IN)                  :: nshell
      CHARACTER(len=*), DIMENSION(:), INTENT(IN)         :: lnam

      CHARACTER(len=6), DIMENSION(:), POINTER            :: sym
      INTEGER                                            :: i, l, nj
      INTEGER, DIMENSION(4, 7)                           :: ne
      INTEGER, DIMENSION(:), POINTER                     :: lq, nq
      REAL(KIND=dp), DIMENSION(:), POINTER               :: zet
      TYPE(sto_basis_set_type), POINTER                  :: sto_basis_set

      CPASSERT(.NOT. ASSOCIATED(pro_basis_set))
      NULLIFY (sto_basis_set)

      ! electronic configuration
      ne = 0
      DO l = 1, 4 !lq(1)+1
         nj = 2*(l - 1) + 1
         DO i = l, 7 ! nq(1)
            ne(l, i) = ptable(zval)%e_conv(l - 1) - 2*nj*(i - l)
            ne(l, i) = MAX(ne(l, i), 0)
            ne(l, i) = MIN(ne(l, i), 2*nj)
         END DO
      END DO
      ALLOCATE (nq(ishell), lq(ishell), zet(ishell), sym(ishell))
      DO i = 1, ishell
         nq(i) = nshell(i)
         SELECT CASE (lnam(i))
         CASE ('S', 's')
            lq(i) = 0
         CASE ('P', 'p')
            lq(i) = 1
         CASE ('D', 'd')
            lq(i) = 2
         CASE ('F', 'f')
            lq(i) = 3
         CASE DEFAULT
            CPABORT("Wrong l QN")
         END SELECT
         sym(i) = lnam(i)
         zet(i) = srules(zval, ne, nq(1), lq(1))
      END DO
      CALL allocate_sto_basis_set(sto_basis_set)
      CALL set_sto_basis_set(sto_basis_set, nshell=1, nq=nq, lq=lq, zet=zet, symbol=sym)
      CALL create_gto_from_sto_basis(sto_basis_set, pro_basis_set, 6)
      pro_basis_set%norm_type = 2
      CALL init_orb_basis_set(pro_basis_set)
      CALL deallocate_sto_basis_set(sto_basis_set)

   END SUBROUTINE create_pro_basis

! **************************************************************************************************
!> \brief Update the density matrix in AO basis with the active density contribution
!> \param active_space_env the active space environment
!> \param rho_ao the density matrix in AO basis
! **************************************************************************************************
   SUBROUTINE update_density_ao(active_space_env, rho_ao)
      TYPE(active_space_type), POINTER                   :: active_space_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ao

      INTEGER                                            :: ispin, nao, nmo, nspins
      TYPE(cp_fm_type)                                   :: R, U
      TYPE(cp_fm_type), POINTER                          :: C_active, p_active_mo
      TYPE(dbcsr_type), POINTER                          :: p_inactive_ao
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos_active

      ! Transform the AS density matrix P_MO to the atomic orbital basis,
      ! this is simply C * P_MO * C^T
      nspins = active_space_env%nspins
      mos_active => active_space_env%mos_active
      DO ispin = 1, nspins
         ! size of p_inactive_ao is (nao x nao)
         p_inactive_ao => active_space_env%pmat_inactive(ispin)%matrix

         ! copy p_inactive_ao to rho_ao
         CALL dbcsr_copy(rho_ao(ispin)%matrix, p_inactive_ao)

         ! size of p_active_mo is (nmo x nmo)
         p_active_mo => active_space_env%p_active(ispin)

         ! calculate R = p_mo
         CALL cp_fm_create(R, p_active_mo%matrix_struct)
         CALL cp_fm_to_fm(p_active_mo, R)

         ! calculate U = C * p_mo
         CALL get_mo_set(mos_active(ispin), mo_coeff=C_active, nao=nao, nmo=nmo)
         CALL cp_fm_create(U, C_active%matrix_struct)
         CALL parallel_gemm("N", "N", nao, nmo, nmo, 1.0_dp, C_active, R, 0.0_dp, U)

         CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=rho_ao(ispin)%matrix, &
                                    matrix_v=U, matrix_g=C_active, ncol=nmo, alpha=1.0_dp)

         CALL cp_fm_release(R)
         CALL cp_fm_release(U)
      END DO

   END SUBROUTINE update_density_ao

! **************************************************************************************************
!> \brief Print each value on the master node
!> \param this object reference
!> \param i i-index
!> \param j j-index
!> \param k k-index
!> \param l l-index
!> \param val value of the integral at (i,j,k.l)
!> \return always true to dump all integrals
! **************************************************************************************************
   LOGICAL FUNCTION eri_fcidump_print_func(this, i, j, k, l, val) RESULT(cont)
      CLASS(eri_fcidump_print), INTENT(inout) :: this
      INTEGER, INTENT(in)                     :: i, j, k, l
      REAL(KIND=dp), INTENT(in)               :: val

      ! write to the actual file only on the master
      IF (this%unit_nr > 0) THEN
         WRITE (this%unit_nr, "(ES23.16,4I4)") val, i + this%bra_start - 1, j + this%bra_start - 1, &
              &                                     k + this%ket_start - 1, l + this%ket_start - 1
      END IF

      cont = .TRUE.
   END FUNCTION eri_fcidump_print_func

! **************************************************************************************************
!> \brief checksum each value on the master node
!> \param this object reference
!> \param i i-index
!> \param j j-index
!> \param k k-index
!> \param l l-index
!> \param val value of the integral at (i,j,k.l)
!> \return always true to dump all integrals
! **************************************************************************************************
   LOGICAL FUNCTION eri_fcidump_checksum_func(this, i, j, k, l, val) RESULT(cont)
      CLASS(eri_fcidump_checksum), INTENT(inout) :: this
      INTEGER, INTENT(in)                     :: i, j, k, l
      REAL(KIND=dp), INTENT(in)               :: val
      MARK_USED(i)
      MARK_USED(j)
      MARK_USED(k)
      MARK_USED(l)

      this%checksum = this%checksum + ABS(val)

      cont = .TRUE.
   END FUNCTION eri_fcidump_checksum_func

! **************************************************************************************************
!> \brief Update active space density matrix from a fortran array
!> \param p_act_mo density matrix in active space MO basis
!> \param active_space_env active space environment
!> \param ispin spin index
!> \author Vladimir Rybkin
! **************************************************************************************************
   SUBROUTINE update_active_density(p_act_mo, active_space_env, ispin)
      REAL(KIND=dp), DIMENSION(:)                        :: p_act_mo
      TYPE(active_space_type), POINTER                   :: active_space_env
      INTEGER                                            :: ispin

      INTEGER                                            :: i1, i2, m1, m2, nmo_active
      REAL(KIND=dp)                                      :: alpha, pij_new, pij_old
      TYPE(cp_fm_type), POINTER                          :: p_active

      p_active => active_space_env%p_active(ispin)
      nmo_active = active_space_env%nmo_active
      alpha = active_space_env%alpha

      DO i1 = 1, nmo_active
         m1 = active_space_env%active_orbitals(i1, ispin)
         DO i2 = 1, nmo_active
            m2 = active_space_env%active_orbitals(i2, ispin)
            CALL cp_fm_get_element(p_active, m1, m2, pij_old)
            pij_new = p_act_mo(i2 + (i1 - 1)*nmo_active)
            pij_new = alpha*pij_new + (1.0_dp - alpha)*pij_old
            CALL cp_fm_set_element(p_active, m1, m2, pij_new)
         END DO
      END DO

   END SUBROUTINE update_active_density

! **************************************************************************************************
!> \brief Compute and print the AS rdm and the natural orbitals occupation numbers
!> \param active_space_env active space environment
!> \param iw output unit
!> \author Stefano Battaglia
! **************************************************************************************************
   SUBROUTINE print_pmat_noon(active_space_env, iw)
      TYPE(active_space_type), POINTER                   :: active_space_env
      INTEGER                                            :: iw

      INTEGER                                            :: i1, i2, ii, ispin, jm, m1, m2, &
                                                            nmo_active, nspins
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: noon, pmat
      TYPE(cp_fm_type), POINTER                          :: p_active

      nspins = active_space_env%nspins
      nmo_active = active_space_env%nmo_active

      ALLOCATE (noon(nmo_active, nspins))
      ALLOCATE (pmat(nmo_active, nmo_active))

      DO ispin = 1, nspins
         p_active => active_space_env%p_active(ispin)
         noon(:, ispin) = 0.0_dp
         pmat = 0.0_dp

         DO i1 = 1, nmo_active
            m1 = active_space_env%active_orbitals(i1, ispin)
            DO i2 = 1, nmo_active
               m2 = active_space_env%active_orbitals(i2, ispin)
               CALL cp_fm_get_element(p_active, m1, m2, pmat(i1, i2))
            END DO
         END DO

         IF (iw > 0) THEN
            WRITE (iw, '(/,T3,A,I2,A)') "Active space density matrix for spin ", ispin
            DO i1 = 1, nmo_active
               DO ii = 1, nmo_active, 8
                  jm = MIN(7, nmo_active - ii)
                  WRITE (iw, '(T3,6(F9.4))') (pmat(i1, ii + i2), i2=0, jm)
               END DO
            END DO
         END IF

         ! diagonalize the density matrix
         CALL diamat_all(pmat, noon(:, ispin))

         IF (iw > 0) THEN
            WRITE (iw, '(/,T3,A,I2,A)') "Natural orbitals occupation numbers for spin ", ispin
            DO i1 = 1, nmo_active, 8
               jm = MIN(7, nmo_active - i1)
               ! noons are stored in ascending order, so reverse-print them
               WRITE (iw, '(T3,6(F9.4))') (noon(nmo_active - i1 - i2 + 1, ispin), i2=0, jm)
            END DO
         END IF

      END DO

      DEALLOCATE (noon)
      DEALLOCATE (pmat)

   END SUBROUTINE print_pmat_noon

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param active_space_env ...
!> \param as_input ...
! **************************************************************************************************
   SUBROUTINE rsdft_embedding(qs_env, active_space_env, as_input)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(active_space_type), POINTER                   :: active_space_env
      TYPE(section_vals_type), POINTER                   :: as_input

      CHARACTER(len=*), PARAMETER                        :: routineN = 'rsdft_embedding'
      INTEGER                                            :: handle

#ifdef __NO_SOCKETS
      CALL timeset(routineN, handle)
      CPABORT("CP2K was compiled with the __NO_SOCKETS option!")
      MARK_USED(qs_env)
      MARK_USED(active_space_env)
      MARK_USED(as_input)
#else

      INTEGER                                            :: iw, client_fd, socket_fd, iter, max_iter
      LOGICAL                                            :: converged, do_scf_embedding, ionode
      REAL(KIND=dp)                                      :: delta_E, energy_corr, energy_new, &
                                                            energy_old, energy_scf, eps_iter, t1, t2
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ao
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos_active
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(dft_control_type), POINTER                    :: dft_control

      CALL timeset(routineN, handle)

      t1 = m_walltime()

      logger => cp_get_default_logger()
      iw = cp_logger_get_default_io_unit(logger)

      CALL get_qs_env(qs_env, para_env=para_env, dft_control=dft_control)
      ionode = para_env%is_source()

      ! get info from the input
      CALL section_vals_val_get(as_input, "SCF_EMBEDDING", l_val=do_scf_embedding)
      active_space_env%do_scf_embedding = do_scf_embedding
      CALL section_vals_val_get(as_input, "MAX_ITER", i_val=max_iter)
      IF (max_iter < 0) CPABORT("Specify a non-negative number of max iterations.")
      CALL section_vals_val_get(as_input, "EPS_ITER", r_val=eps_iter)
      IF (eps_iter < 0.0) CPABORT("Specify a non-negative convergence threshold.")

      ! create the socket and wait for the client to connect
      CALL initialize_socket(socket_fd, client_fd, as_input, ionode)
      CALL para_env%sync()

      ! send two-electron integrals to the client
      CALL send_eri_to_client(client_fd, active_space_env, para_env)

      ! get pointer to density in ao basis
      CALL get_qs_env(qs_env, rho=rho, energy=energy, ks_env=ks_env)
      CALL qs_rho_get(rho, rho_ao=rho_ao)

      IF (iw > 0) THEN
         WRITE (UNIT=iw, FMT="(/,T2,A,/)") &
            "RANGE-SEPARATED DFT EMBEDDING SELF-CONSISTENT OPTIMIZATION"

         WRITE (iw, '(T3,A,T68,I12)') "Max. iterations", max_iter
         WRITE (iw, '(T3,A,T68,E12.4)') "Conv. threshold", eps_iter
         WRITE (iw, '(T3,A,T66,F14.2)') "Density damping", active_space_env%alpha

         WRITE (UNIT=iw, FMT="(/,T3,A,T11,A,T21,A,T34,A,T55,A,T75,A,/,T3,A)") &
            "Iter", "Update", "Time", "Corr. energy", "Total energy", "Change", REPEAT("-", 78)
      END IF
      ! CALL cp_add_iter_level(logger%iter_info, "QS_SCF")

      iter = 0
      converged = .FALSE.
      ! store the scf energy
      energy_scf = active_space_env%energy_ref
      energy_new = energy_scf
      mos_active => active_space_env%mos_active
      ! CALL set_qs_env(qs_env, active_space=active_space_env)

      ! start the self-consistent embedding loop
      DO WHILE (iter < max_iter)
         iter = iter + 1

         ! send V_emb and E_ina to the active space solver and update
         ! the active space environment with the new active energy and density
         CALL send_fock_to_client(client_fd, active_space_env, para_env)

         ! update energies
         energy_old = energy_new
         energy_new = active_space_env%energy_total
         energy_corr = energy_new - energy_scf
         delta_E = energy_new - energy_old

         ! get timer
         t2 = t1
         t1 = m_walltime()
         ! print out progress
         IF ((iw > 0)) THEN
            WRITE (UNIT=iw, &
                   FMT="(T3,I4,T11,A,T19,F6.1,T28,F18.10,T49,F18.10,T70,ES11.2)") &
               iter, 'P_Mix', t1 - t2, energy_corr, energy_new, delta_E
            CALL m_flush(iw)
         END IF

         ! update total density in AO basis with the AS contribution
         CALL update_density_ao(active_space_env, rho_ao) ! rho_ao is updated

         ! calculate F_ks in AO basis (which contains Vxc) with the new density
         CALL qs_rho_update_rho(rho, qs_env=qs_env) ! updates rho_r and rho_g using rho_ao
         CALL qs_ks_did_change(qs_env%ks_env, rho_changed=.TRUE.) ! set flags about the change
         ! Re-evaluate the traces between the density matrix and the core Hamiltonians
         CALL evaluate_core_matrix_traces(qs_env)
         ! the ks matrix will be rebuilt so this is fine now
         ! CALL set_ks_env(qs_env%ks_env, potential_changed=.FALSE.)
         CALL qs_ks_build_kohn_sham_matrix(qs_env, calculate_forces=.FALSE., &
                                           just_energy=.FALSE., &
                                           ext_xc_section=active_space_env%xc_section)

         ! update the reference energy
         active_space_env%energy_ref = energy%total

         ! transform KS/Fock, Vxc and Hcore from AO to MO basis
         CALL calculate_operators(mos_active, qs_env, active_space_env)

         ! calculate the new inactive energy and embedding potential
         CALL subspace_fock_matrix(active_space_env, dft_control%roks)

         ! check if it is a one-shot correction
         IF (.NOT. active_space_env%do_scf_embedding) THEN
            IF (iw > 0) THEN
               WRITE (UNIT=iw, FMT="(/,T3,A,I5,A)") &
                  "*** one-shot embedding correction finished ***"
            END IF
            converged = .TRUE.
            EXIT
            ! check for convergence
         ELSEIF (ABS(delta_E) <= eps_iter) THEN
            IF (iw > 0) THEN
               WRITE (UNIT=iw, FMT="(/,T3,A,I5,A)") &
                  "*** rs-DFT embedding run converged in ", iter, " iteration(s) ***"
            END IF
            converged = .TRUE.
            EXIT
         END IF
      END DO

      IF (.NOT. converged) THEN
         IF (iw > 0) THEN
            WRITE (UNIT=iw, FMT="(/,T3,A,I5,A)") &
               "*** rs-DFT embedding did not converged after ", iter, " iteration(s) ***"
         END IF
      END IF

      ! update qs total energy to the final rs-DFT energy
      energy%total = active_space_env%energy_total

      ! print final energy contributions
      IF (iw > 0) THEN
         WRITE (UNIT=iw, FMT="(/,T3,A)") &
            "Final energy contributions:"
         WRITE (UNIT=iw, FMT="(T6,A,T56,F20.10)") &
            "Inactive energy:", active_space_env%energy_inactive
         WRITE (UNIT=iw, FMT="(T6,A,T56,F20.10)") &
            "Active energy:", active_space_env%energy_active
         WRITE (UNIT=iw, FMT="(T6,A,T56,F20.10)") &
            "Correlation energy:", energy_corr
         WRITE (UNIT=iw, FMT="(T6,A,T56,F20.10)") &
            "Total rs-DFT energy:", active_space_env%energy_total
      END IF

      ! print the AS rdm and the natural orbital occupation numbers
      CALL print_pmat_noon(active_space_env, iw)

      CALL finalize_socket(socket_fd, client_fd, as_input, ionode)
      CALL para_env%sync()
#endif

      CALL timestop(handle)

   END SUBROUTINE rsdft_embedding

#ifndef __NO_SOCKETS
! **************************************************************************************************
!> \brief Creates the socket, spawns the client and connects to it
!> \param socket_fd the socket file descriptor
!> \param client_fd the client file descriptor
!> \param as_input active space inpute section
!> \param ionode logical flag indicating if the process is the master
! **************************************************************************************************
   SUBROUTINE initialize_socket(socket_fd, client_fd, as_input, ionode)
      INTEGER, INTENT(OUT)                               :: socket_fd, client_fd
      TYPE(section_vals_type), INTENT(IN), POINTER       :: as_input
      LOGICAL, INTENT(IN)                                :: ionode

      CHARACTER(len=*), PARAMETER                        :: routineN = 'initialize_socket'
      INTEGER, PARAMETER                                 :: backlog = 10

      CHARACTER(len=default_path_length)                 :: hostname
      INTEGER                                            :: handle, iw, port, protocol
      LOGICAL                                            :: inet
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      iw = cp_logger_get_default_io_unit(logger)

      ! protocol == 0 for UNIX, protocol > 0 for INET
      CALL section_vals_val_get(as_input, "SOCKET%INET", l_val=inet)
      IF (inet) THEN
         protocol = 1
      ELSE
         protocol = 0
      END IF
      CALL section_vals_val_get(as_input, "SOCKET%HOST", c_val=hostname)
      CALL section_vals_val_get(as_input, "SOCKET%PORT", i_val=port)

      IF (ionode) THEN
         CALL open_bind_socket(socket_fd, protocol, port, TRIM(hostname)//C_NULL_CHAR)
         WRITE (iw, '(/,T2,A,A)') "@SERVER: Created socket with address ", TRIM(hostname)
         CALL listen_socket(socket_fd, backlog)

         ! wait until a connetion request arrives
         WRITE (iw, '(T2,A)') "@SERVER: Waiting for requests..."
         CALL accept_socket(socket_fd, client_fd)
         WRITE (iw, '(T2,A,I2)') "@SERVER: Accepted socket with fd ", client_fd
      END IF

      CALL timestop(handle)

   END SUBROUTINE initialize_socket

! **************************************************************************************************
!> \brief Closes the connection to the socket and deletes the file
!> \param socket_fd the socket file descriptor
!> \param client_fd the client file descriptor
!> \param as_input active space inpute section
!> \param ionode logical flag indicating if the process is the master
! **************************************************************************************************
   SUBROUTINE finalize_socket(socket_fd, client_fd, as_input, ionode)
      INTEGER, INTENT(IN)                                :: socket_fd, client_fd
      TYPE(section_vals_type), INTENT(IN), POINTER       :: as_input
      LOGICAL, INTENT(IN)                                :: ionode

      CHARACTER(len=*), PARAMETER                        :: routineN = 'finalize_socket'
      INTEGER, PARAMETER                                 :: header_len = 12

      CHARACTER(len=default_path_length)                 :: hostname
      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL section_vals_val_get(as_input, "SOCKET%HOST", c_val=hostname)

      IF (ionode) THEN
         ! signal the client to quit
         CALL writebuffer(client_fd, "QUIT        ", header_len)
         ! close the connection
         CALL close_socket(client_fd)
         CALL close_socket(socket_fd)

         ! delete the socket file
         IF (file_exists(TRIM(hostname))) THEN
            CALL remove_socket_file(TRIM(hostname)//C_NULL_CHAR)
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE finalize_socket

! **************************************************************************************************
!> \brief Sends the two-electron integrals to the client vie the socket
!> \param client_fd the client file descriptor
!> \param active_space_env active space environment
!> \param para_env parallel environment
! **************************************************************************************************
   SUBROUTINE send_eri_to_client(client_fd, active_space_env, para_env)
      INTEGER, INTENT(IN)                                :: client_fd
      TYPE(active_space_type), INTENT(IN), POINTER       :: active_space_env
      TYPE(mp_para_env_type), INTENT(IN), POINTER        :: para_env

      CHARACTER(len=*), PARAMETER :: routineN = 'send_eri_to_client'
      INTEGER, PARAMETER                                 :: header_len = 12

      CHARACTER(len=default_string_length)               :: header
      INTEGER                                            :: handle, iw
      LOGICAL                                            :: ionode
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eri_aa, eri_ab, eri_bb, s_ab
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      iw = cp_logger_get_default_io_unit(logger)
      ionode = para_env%is_source()

      ALLOCATE (eri_aa(active_space_env%nmo_active**4))
      CALL eri_to_array(active_space_env%eri, eri_aa, active_space_env%active_orbitals, 1, 1)
      IF (active_space_env%nspins == 2) THEN
         ALLOCATE (eri_ab(active_space_env%nmo_active**4))
         CALL eri_to_array(active_space_env%eri, eri_ab, active_space_env%active_orbitals, 1, 2)
         ALLOCATE (eri_bb(active_space_env%nmo_active**4))
         CALL eri_to_array(active_space_env%eri, eri_bb, active_space_env%active_orbitals, 2, 2)
         ! get the overlap_ab matrix into Fortran array
         ALLOCATE (s_ab(active_space_env%nmo_active**2))
         ASSOCIATE (act_indices_a => active_space_env%active_orbitals(:, 1), &
                    act_indices_b => active_space_env%active_orbitals(:, 2))
            CALL subspace_matrix_to_array(active_space_env%sab_sub(1), s_ab, act_indices_a, act_indices_b)
         END ASSOCIATE
      END IF

      ! ask the status of the client
      IF (ionode) CALL writebuffer(client_fd, "STATUS      ", header_len)
      DO
         header = ""
         CALL para_env%sync()
         IF (ionode) THEN
            ! IF (iw > 0) WRITE(iw, *) "@SERVER: Waiting for messages..."
            CALL readbuffer(client_fd, header, header_len)
         END IF
         CALL para_env%bcast(header, para_env%source)

         ! IF (iw > 0) WRITE(iw, *) "@SERVER: Message from client: ", TRIM(header)

         IF (TRIM(header) == "READY") THEN
            ! if the client is ready, send the data
            CALL para_env%sync()
            IF (ionode) THEN
               CALL writebuffer(client_fd, "TWOBODY     ", header_len)
               CALL writebuffer(client_fd, active_space_env%nspins)
               CALL writebuffer(client_fd, active_space_env%nmo_active)
               CALL writebuffer(client_fd, active_space_env%nelec_active)
               CALL writebuffer(client_fd, active_space_env%multiplicity)
               ! send the alpha component
               CALL writebuffer(client_fd, eri_aa, SIZE(eri_aa))
               ! send the beta part for unrestricted calculations
               IF (active_space_env%nspins == 2) THEN
                  CALL writebuffer(client_fd, eri_ab, SIZE(eri_ab))
                  CALL writebuffer(client_fd, eri_bb, SIZE(eri_bb))
                  CALL writebuffer(client_fd, s_ab, SIZE(s_ab))
               END IF
            END IF
         ELSE IF (TRIM(header) == "RECEIVED") THEN
            EXIT
         END IF
      END DO

      DEALLOCATE (eri_aa)
      IF (active_space_env%nspins == 2) THEN
         DEALLOCATE (eri_ab)
         DEALLOCATE (eri_bb)
         DEALLOCATE (s_ab)
      END IF

      CALL para_env%sync()

      CALL timestop(handle)

   END SUBROUTINE send_eri_to_client

! **************************************************************************************************
!> \brief Sends the one-electron embedding potential and the inactive energy to the client
!> \param client_fd the client file descriptor
!> \param active_space_env active space environment
!> \param para_env parallel environment
! **************************************************************************************************
   SUBROUTINE send_fock_to_client(client_fd, active_space_env, para_env)
      INTEGER, INTENT(IN)                                :: client_fd
      TYPE(active_space_type), INTENT(IN), POINTER       :: active_space_env
      TYPE(mp_para_env_type), INTENT(IN), POINTER        :: para_env

      CHARACTER(len=*), PARAMETER :: routineN = 'send_fock_to_client'
      INTEGER, PARAMETER                                 :: header_len = 12

      CHARACTER(len=default_string_length)               :: header
      INTEGER                                            :: handle, iw
      LOGICAL                                            :: debug, ionode
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: fock_a, fock_b, p_act_mo_a, p_act_mo_b
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)

      ! Set to .TRUE. to activate debug output
      debug = .FALSE.

      logger => cp_get_default_logger()
      iw = cp_logger_get_default_io_unit(logger)
      ionode = para_env%is_source()

      ALLOCATE (p_act_mo_a(active_space_env%nmo_active**2))
      ALLOCATE (fock_a(active_space_env%nmo_active**2))
      IF (active_space_env%nspins == 2) THEN
         ALLOCATE (p_act_mo_b(active_space_env%nmo_active**2))
         ALLOCATE (fock_b(active_space_env%nmo_active**2))
      END IF

      ! get the fock matrix into Fortran arrays
      ASSOCIATE (act_indices => active_space_env%active_orbitals(:, 1))
         CALL subspace_matrix_to_array(active_space_env%fock_sub(1), fock_a, act_indices, act_indices)
      END ASSOCIATE

      IF (active_space_env%nspins == 2) THEN
         ASSOCIATE (act_indices => active_space_env%active_orbitals(:, 2))
            CALL subspace_matrix_to_array(active_space_env%fock_sub(2), fock_b, act_indices, act_indices)
         END ASSOCIATE
      END IF

      ! ask the status of the client
      IF (ionode) CALL writebuffer(client_fd, "STATUS      ", header_len)
      DO
         header = ""

         CALL para_env%sync()
         IF (ionode) THEN
            IF (debug .AND. iw > 0) WRITE (iw, *) "@SERVER: Waiting for messages..."
            CALL readbuffer(client_fd, header, header_len)
         END IF
         CALL para_env%bcast(header, para_env%source)

         IF (debug .AND. iw > 0) WRITE (iw, *) "@SERVER: Message from client: ", TRIM(header)

         IF (TRIM(header) == "READY") THEN
            ! if the client is ready, send the data
            CALL para_env%sync()
            IF (ionode) THEN
               CALL writebuffer(client_fd, "ONEBODY     ", header_len)
               CALL writebuffer(client_fd, active_space_env%energy_inactive)
               ! send the alpha component
               CALL writebuffer(client_fd, fock_a, SIZE(fock_a))
               ! send the beta part for unrestricted calculations
               IF (active_space_env%nspins == 2) THEN
                  CALL writebuffer(client_fd, fock_b, SIZE(fock_b))
               END IF
            END IF

         ELSE IF (TRIM(header) == "HAVEDATA") THEN
            ! qiskit has data to transfer, let them know we want it and wait for it
            CALL para_env%sync()
            IF (ionode) THEN
               IF (debug .AND. iw > 0) WRITE (iw, *) "@SERVER: Qiskit has data to transfer"
               CALL writebuffer(client_fd, "GETDENSITY  ", header_len)

               ! read the active energy and density
               CALL readbuffer(client_fd, active_space_env%energy_active)
               CALL readbuffer(client_fd, p_act_mo_a, SIZE(p_act_mo_a))
               IF (active_space_env%nspins == 2) THEN
                  CALL readbuffer(client_fd, p_act_mo_b, SIZE(p_act_mo_b))
               END IF
            END IF

            ! broadcast the data to all processors
            CALL para_env%bcast(active_space_env%energy_active, para_env%source)
            CALL para_env%bcast(p_act_mo_a, para_env%source)
            IF (active_space_env%nspins == 2) THEN
               CALL para_env%bcast(p_act_mo_b, para_env%source)
            END IF

            ! update total and reference energies in active space enviornment
            active_space_env%energy_total = active_space_env%energy_inactive + active_space_env%energy_active

            ! update the active density matrix in the active space environment
            CALL update_active_density(p_act_mo_a, active_space_env, 1)
            IF (active_space_env%nspins == 2) THEN
               CALL update_active_density(p_act_mo_b, active_space_env, 2)
            END IF

            ! the non-iterative part is done, we can continue
            EXIT
         END IF

      END DO

      DEALLOCATE (p_act_mo_a)
      DEALLOCATE (fock_a)
      IF (active_space_env%nspins == 2) THEN
         DEALLOCATE (p_act_mo_b)
         DEALLOCATE (fock_b)
      END IF

      CALL para_env%sync()

      CALL timestop(handle)

   END SUBROUTINE send_fock_to_client
#endif

END MODULE qs_active_space_methods
