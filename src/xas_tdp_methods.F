!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Methods for X-Ray absorption spectroscopy (XAS) using TDDFPT
!> \author AB (11.2017)
! **************************************************************************************************

MODULE xas_tdp_methods
   USE admm_types,                      ONLY: admm_type
   USE admm_utils,                      ONLY: admm_correct_for_eigenvalues,&
                                              admm_uncorrect_for_eigenvalues
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE basis_set_types,                 ONLY: &
        allocate_sto_basis_set, create_gto_from_sto_basis, deallocate_gto_basis_set, &
        deallocate_sto_basis_set, get_gto_basis_set, gto_basis_set_type, init_orb_basis_set, &
        set_sto_basis_set, srules, sto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              cp_dbcsr_sm_fm_multiply
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_fm_diag,                      ONLY: cp_fm_power,&
                                              cp_fm_syevd
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: &
        cp_fm_create, cp_fm_get_diag, cp_fm_get_info, cp_fm_get_submatrix, cp_fm_p_type, &
        cp_fm_read_unformatted, cp_fm_release, cp_fm_set_all, cp_fm_to_fm, cp_fm_to_fm_submat, &
        cp_fm_type, cp_fm_write_unformatted
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_generate_filename,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr,&
                                              debug_print_level
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_add_on_diag,&
                                              dbcsr_copy,&
                                              dbcsr_create,&
                                              dbcsr_finalize,&
                                              dbcsr_get_info,&
                                              dbcsr_multiply,&
                                              dbcsr_p_type,&
                                              dbcsr_set
   USE hfx_libint_wrapper,              ONLY: cp_libint_static_init
   USE input_constants,                 ONLY: &
        do_loc_none, do_potential_coulomb, do_potential_short, do_potential_truncated, &
        op_loc_berry, state_loc_list, tddfpt_singlet, tddfpt_spin_cons, tddfpt_spin_flip, &
        tddfpt_triplet, xas_1s_type, xas_2p_type, xas_2s_type, xas_dip_len, xas_dip_vel, &
        xas_not_excited, xas_tdp_by_index, xas_tdp_by_kind, xas_tdp_diag_syevr
   USE input_cp2k_dft,                  ONLY: create_localize_section
   USE input_section_types,             ONLY: &
        section_release, section_type, section_vals_create, section_vals_duplicate, &
        section_vals_get, section_vals_get_subs_vals, section_vals_get_subs_vals2, &
        section_vals_release, section_vals_type, section_vals_val_get, section_vals_val_set
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE machine,                         ONLY: m_flush
   USE mathlib,                         ONLY: get_diag
   USE memory_utilities,                ONLY: reallocate
   USE message_passing,                 ONLY: mp_bcast
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE periodic_table,                  ONLY: ptable
   USE physcon,                         ONLY: a_fine,&
                                              angstrom,&
                                              evolt
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_interactions,                 ONLY: init_interaction_radii_orb_basis
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_loc_methods,                  ONLY: centers_spreads_berry,&
                                              qs_loc_driver,&
                                              qs_print_cubes
   USE qs_loc_types,                    ONLY: get_qs_loc_env,&
                                              localized_wfn_control_create,&
                                              localized_wfn_control_type,&
                                              qs_loc_env_create,&
                                              qs_loc_env_new_type,&
                                              qs_loc_env_release
   USE qs_loc_utils,                    ONLY: qs_loc_control_init,&
                                              qs_loc_env_init,&
                                              set_loc_centers
   USE qs_mo_methods,                   ONLY: calculate_subspace_eigenvalues
   USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                              deallocate_mo_set,&
                                              get_mo_set,&
                                              init_mo_set,&
                                              mo_set_p_type,&
                                              mo_set_type
   USE qs_operators_ao,                 ONLY: p_xyz_ao,&
                                              rRc_xyz_ao
   USE qs_pdos,                         ONLY: calculate_projected_dos
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_scf_types,                    ONLY: ot_method_nr
   USE util,                            ONLY: locate,&
                                              sort_unique
   USE xas_methods,                     ONLY: calc_stogto_overlap
   USE xas_tdp_atom,                    ONLY: init_xas_atom_env,&
                                              integrate_fxc_atoms,&
                                              integrate_soc_atoms
   USE xas_tdp_types,                   ONLY: &
        donor_state_create, donor_state_type, free_ds_memory, free_exat_memory, &
        read_xas_tdp_control, set_donor_state, set_xas_tdp_env, xas_atom_env_create, &
        xas_atom_env_release, xas_atom_env_type, xas_tdp_control_create, xas_tdp_control_type, &
        xas_tdp_env_create, xas_tdp_env_release, xas_tdp_env_type
   USE xas_tdp_utils,                   ONLY: compute_ri_coulomb2_int,&
                                              compute_ri_exchange2_int,&
                                              compute_ri_o3c_coulomb,&
                                              compute_ri_o3c_exchange,&
                                              include_os_soc,&
                                              include_rcs_soc,&
                                              setup_xas_tdp_prob,&
                                              solve_xas_tdp_prob
   USE xc_derivatives,                  ONLY: xc_functional_get_info
   USE xc_libxc,                        ONLY: libxc_version_info
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xas_tdp_methods'

   PUBLIC :: xas_tdp

CONTAINS

! **************************************************************************************************
!> \brief Driver for XAS TDDFT calculations.
!> \param qs_env the inherited qs_environment
!> \author AB
!> \note Empty for now...
! **************************************************************************************************
   SUBROUTINE xas_tdp(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'xas_tdp', routineP = moduleN//':'//routineN

      CHARACTER(default_string_length)                   :: rst_filename
      INTEGER                                            :: handle, n_rep, output_unit
      LOGICAL                                            :: do_restart
      TYPE(section_vals_type), POINTER                   :: xas_tdp_section

      CALL timeset(routineN, handle)

!  Logger initialization and XAS TDP banner printing
      NULLIFY (xas_tdp_section)
      xas_tdp_section => section_vals_get_subs_vals(qs_env%input, "DFT%XAS_TDP")
      output_unit = cp_logger_get_default_io_unit()

      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T3,A,/,T3,A,/,T3,A,/,T3,A,/)") &
            "!===========================================================================!", &
            "!                              XAS TDP                                      !", &
            "!    Starting TDDFPT driven X-rays absorption spectroscopy calculations     !", &
            "!===========================================================================!"
      END IF

!  Check whether this is a restart calculation, i.e. is a restart file is provided
      CALL section_vals_val_get(xas_tdp_section, "RESTART_FROM_FILE", n_rep_val=n_rep)

      IF (n_rep < 1) THEN
         do_restart = .FALSE.
      ELSE
         CALL section_vals_val_get(xas_tdp_section, "RESTART_FROM_FILE", c_val=rst_filename)
         do_restart = .TRUE.
      END IF

!  Restart the calculation if needed
      IF (do_restart) THEN

         IF (output_unit > 0) THEN
            WRITE (UNIT=output_unit, FMT="(/,T3,A)") &
               "# This is a RESTART calculation for PDOS and/or CUBE printing"
         END IF

         CALL restart_calculation(rst_filename, xas_tdp_section, qs_env)

!  or run the core XAS_TDP routine if not
      ELSE
         CALL xas_tdp_core(xas_tdp_section, qs_env)
      END IF

      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T3,A,/,T3,A,/,T3,A,/)") &
            "!===========================================================================!", &
            "!     End of TDDFPT driven X-rays absorption spectroscopy calculations      !", &
            "!===========================================================================!"
      END IF

      CALL timestop(handle)

   END SUBROUTINE xas_tdp

! **************************************************************************************************
!> \brief The core workflow of the XAS_TDP method
!> \param xas_tdp_section the input values for XAS_TDP
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE xas_tdp_core(xas_tdp_section, qs_env)

      TYPE(section_vals_type), POINTER                   :: xas_tdp_section
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'xas_tdp_core', routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_string_length)               :: kind_name
      INTEGER                                            :: current_state_index, iat, iatom, ikind, &
                                                            ispin, istate, output_unit, tmp_index
      INTEGER, DIMENSION(:), POINTER                     :: atoms_of_kind
      LOGICAL                                            :: do_os
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(donor_state_type), POINTER                    :: current_state
      TYPE(gto_basis_set_type), POINTER                  :: tmp_basis
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(xas_atom_env_type), POINTER                   :: xas_atom_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env

      NULLIFY (xas_tdp_env, xas_tdp_control, atomic_kind_set, atoms_of_kind, current_state)
      NULLIFY (xas_atom_env, dft_control, matrix_ks, admm_env, qs_kind_set, tmp_basis)

!  Initialization
      output_unit = cp_logger_get_default_io_unit()

      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T3,A)") &
            "# Create and initialize the XAS_TDP environment"
      END IF
      CALL get_qs_env(qs_env, dft_control=dft_control)
      CALL xas_tdp_init(xas_tdp_env, xas_tdp_control, qs_env)
      CALL print_info(output_unit, xas_tdp_control, qs_env)

      IF (output_unit > 0) THEN
         IF (xas_tdp_control%check_only) THEN
            CPWARN("This is a CHECK_ONLY run for donor MOs verification")
         END IF
      END IF

!  Localization of the core orbitals if requested (used for better identification of donor states)
      IF (xas_tdp_control%do_loc) THEN
         IF (output_unit > 0) THEN
            WRITE (UNIT=output_unit, FMT="(/,T3,A,/)") &
               "# Localizing core orbitals for better identification"
         END IF
!        closed shell or ROKS => myspin=1
         IF (xas_tdp_control%do_uks) THEN
            CALL qs_loc_driver(qs_env, xas_tdp_env%qs_loc_env, xas_tdp_control%print_loc_subsection)
         ELSE
            CALL qs_loc_driver(qs_env, xas_tdp_env%qs_loc_env, &
                               xas_tdp_control%print_loc_subsection, myspin=1)
         END IF
      END IF

!  Find the MO centers
      CALL find_mo_centers(xas_tdp_env, xas_tdp_control, qs_env)

!  Assign lowest energy orbitals to excited atoms
      CALL assign_mos_to_ex_atoms(xas_tdp_env, xas_tdp_control, qs_env)

!  Once assigned, diagonalize the MOs wrt the KS matrix in the subspace associated to each atom
      IF (xas_tdp_control%do_loc) THEN
         IF (output_unit > 0) THEN
            WRITE (UNIT=output_unit, FMT="(/,T3,A,/,T5,A)") &
               "# Diagonalize localized MOs wrt the KS matrix in the subspace of each excited", &
               "atom for better donor state identification."
         END IF
         CALL diagonalize_assigned_mo_subset(xas_tdp_env, xas_tdp_control, qs_env)
         ! update MO centers
         CALL find_mo_centers(xas_tdp_env, xas_tdp_control, qs_env)
      END IF

      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T3,A,I4,A,/)") &
            "# Assign the relevant subset of the ", xas_tdp_control%n_search, &
            "  lowest energy MOs to excited atoms"
      END IF
      CALL write_mos_to_ex_atoms_association(xas_tdp_env, xas_tdp_control, qs_env)

!  If CHECK_ONLY run, check the donor MOs
      IF (xas_tdp_control%check_only) CALL print_checks(xas_tdp_env, xas_tdp_control, qs_env)

!  If not simply exact exchange, setup a xas_atom_env and compute the xc integrals on the atomic grids
!  Also needed if SOC is included. Done before looping on atoms as it's all done at once
      IF ((xas_tdp_control%do_xc .OR. xas_tdp_control%do_soc) .AND. .NOT. xas_tdp_control%check_only) THEN

         IF (output_unit > 0 .AND. xas_tdp_control%do_xc) THEN
            WRITE (UNIT=output_unit, FMT="(/,T3,A,I4,A)") &
               "# Integrating the xc kernel on the atomic grids ..."
            CALL m_flush(output_unit)
         END IF

         CALL xas_atom_env_create(xas_atom_env)
         CALL init_xas_atom_env(xas_atom_env, xas_tdp_env, xas_tdp_control, qs_env)
         do_os = xas_tdp_control%do_uks .OR. xas_tdp_control%do_roks

         IF (xas_tdp_control%do_xc) THEN
            CALL integrate_fxc_atoms(xas_tdp_env%ri_fxc, xas_atom_env, xas_tdp_control, qs_env)
         END IF

         IF (xas_tdp_control%do_soc) THEN
            CALL integrate_soc_atoms(xas_tdp_env%orb_soc, xas_atom_env, qs_env)
         END IF

         CALL xas_atom_env_release(xas_atom_env)
      END IF

!  Compute the 3-center Coulomb integrals for the whole system
      IF ((.NOT. xas_tdp_control%check_only) .AND. &
          (xas_tdp_control%do_xc .OR. xas_tdp_control%do_coulomb)) THEN
         IF (output_unit > 0) THEN
            WRITE (UNIT=output_unit, FMT="(/,T3,A,I4,A)") &
               "# Computing the RI 3-center Coulomb integrals ..."
            CALL m_flush(output_unit)
         END IF
         CALL compute_ri_o3c_coulomb(xas_tdp_env, qs_env)

      END IF

!  Loop over donor states for calculation
      CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set)
      current_state_index = 1

!     Loop over atomic kinds
      DO ikind = 1, SIZE(atomic_kind_set)

         IF (xas_tdp_control%check_only) EXIT

         CALL get_atomic_kind(atomic_kind=atomic_kind_set(ikind), name=kind_name, &
                              atom_list=atoms_of_kind)

         IF (.NOT. ANY(xas_tdp_env%ex_kind_indices == ikind)) CYCLE
         ! compute the RI coulomb2 inverse for this kind, and RI exchange2 if needed
         CALL compute_ri_coulomb2_int(ikind, xas_tdp_env, xas_tdp_control, qs_env)
         IF (xas_tdp_control%do_hfx) THEN
            CALL compute_ri_exchange2_int(ikind, xas_tdp_env, xas_tdp_control, qs_env)
         END IF

!        Loop over atoms of kind
         DO iat = 1, SIZE(atoms_of_kind)
            iatom = atoms_of_kind(iat)

            IF (.NOT. ANY(xas_tdp_env%ex_atom_indices == iatom)) CYCLE
            tmp_index = locate(xas_tdp_env%ex_atom_indices, iatom)

            !if dipole in length rep, compute the dipole in the AO basis for this atom
            !if quadrupole is required, compute it there too (in length rep)
            IF (xas_tdp_control%dipole_form == xas_dip_len .OR. xas_tdp_control%do_quad) THEN
               CALL compute_lenrep_multipole(iatom, xas_tdp_env, xas_tdp_control, qs_env)
            END IF

            !compute the three center RI exchange integrals if so required
            IF (xas_tdp_control%do_hfx) THEN
               IF (output_unit > 0) THEN
                  WRITE (UNIT=output_unit, FMT="(/,T3,A,I4,A1,A)") &
                     "# Computing the RI 3-center Exchange integrals for excited atom (idx/kind): ", &
                     iatom, "/", TRIM(kind_name)
                  CALL m_flush(output_unit)
               END IF
               CALL compute_ri_o3c_exchange(iatom, xas_tdp_env, xas_tdp_control, qs_env)
            END IF

!           Loop over states of excited atom of kind
            DO istate = 1, SIZE(xas_tdp_env%state_types, 1)

               IF (xas_tdp_env%state_types(istate, tmp_index) == xas_not_excited) CYCLE

               current_state => xas_tdp_env%donor_states(current_state_index)
               CALL set_donor_state(current_state, at_index=iatom, &
                                    at_symbol=kind_name, kind_index=ikind, &
                                    state_type=xas_tdp_env%state_types(istate, tmp_index))

!              Initial write for the donor state of interest
               IF (output_unit > 0) THEN
                  WRITE (UNIT=output_unit, FMT="(/,T3,A,A2,A,I4,A,A,/)") &
                     "# Start of calculations for donor state of type ", &
                     xas_tdp_env%state_type_char(current_state%state_type), " for atom", &
                     current_state%at_index, " of kind ", TRIM(current_state%at_symbol)
                  CALL m_flush(output_unit)
               END IF

!              Assign best fitting MO(s) to current core donnor state
               CALL assign_mos_to_donor_state(current_state, xas_tdp_env, xas_tdp_control, qs_env)

!              Perform MO restricted Mulliken pop analysis for verification
               CALL perform_mulliken_on_donor_state(current_state, qs_env)

!              Do main calculations here
               CALL setup_xas_tdp_prob(current_state, qs_env, xas_tdp_env, xas_tdp_control)

!              Warn against possible failure of SYEVR if this is the solver of choice
               IF (output_unit > 0 .AND. xas_tdp_control%diag_type == xas_tdp_diag_syevr) THEN
                  WRITE (UNIT=output_unit, FMT="(/T5,A,/,T5,A,/,T5,A)") &
                     "Diagonalizing using ScaLAPACK's PDSYEVR routine. This routine is not always", &
                     "stable and might lead to a crash. If it so happens, switch back to the", &
                     "default DIAG_TYPE STD keyword."
                  CALL m_flush(output_unit)
               END IF

               IF (xas_tdp_control%do_spin_cons) THEN
                  CALL solve_xas_tdp_prob(current_state, xas_tdp_control, xas_tdp_env, qs_env, &
                                          ex_type=tddfpt_spin_cons)
                  CALL compute_dipole_fosc(current_state, xas_tdp_control, xas_tdp_env)
                  IF (xas_tdp_control%do_quad) CALL compute_quadrupole_fosc(current_state, &
                                                                            xas_tdp_control, xas_tdp_env)
                  CALL xas_tdp_post(tddfpt_spin_cons, current_state, xas_tdp_env, &
                                    xas_tdp_section, qs_env)
                  CALL write_donor_state_restart(tddfpt_spin_cons, current_state, xas_tdp_section)
               END IF

               IF (xas_tdp_control%do_spin_flip) THEN
                  CALL solve_xas_tdp_prob(current_state, xas_tdp_control, xas_tdp_env, qs_env, &
                                          ex_type=tddfpt_spin_flip)
                  !no dipole in spin-flip (spin-forbidden)
                  CALL xas_tdp_post(tddfpt_spin_flip, current_state, xas_tdp_env, &
                                    xas_tdp_section, qs_env)
                  CALL write_donor_state_restart(tddfpt_spin_flip, current_state, xas_tdp_section)
               END IF

               IF (xas_tdp_control%do_singlet) THEN
                  CALL solve_xas_tdp_prob(current_state, xas_tdp_control, xas_tdp_env, qs_env, &
                                          ex_type=tddfpt_singlet)
                  CALL compute_dipole_fosc(current_state, xas_tdp_control, xas_tdp_env)
                  IF (xas_tdp_control%do_quad) CALL compute_quadrupole_fosc(current_state, &
                                                                            xas_tdp_control, xas_tdp_env)
                  CALL xas_tdp_post(tddfpt_singlet, current_state, xas_tdp_env, &
                                    xas_tdp_section, qs_env)
                  CALL write_donor_state_restart(tddfpt_singlet, current_state, xas_tdp_section)
               END IF

               IF (xas_tdp_control%do_triplet) THEN
                  CALL solve_xas_tdp_prob(current_state, xas_tdp_control, xas_tdp_env, qs_env, &
                                          ex_type=tddfpt_triplet)
                  !no dipole for triplets by construction
                  CALL xas_tdp_post(tddfpt_triplet, current_state, xas_tdp_env, &
                                    xas_tdp_section, qs_env)
                  CALL write_donor_state_restart(tddfpt_triplet, current_state, xas_tdp_section)
               END IF

!              Include the SOC if required
               IF (xas_tdp_control%do_soc) THEN
                  IF (xas_tdp_control%do_singlet .AND. xas_tdp_control%do_triplet) THEN
                     CALL include_rcs_soc(current_state, xas_tdp_env, xas_tdp_control, qs_env)
                  END IF
                  IF (xas_tdp_control%do_spin_cons .AND. xas_tdp_control%do_spin_flip) THEN
                     CALL include_os_soc(current_state, xas_tdp_env, xas_tdp_control, qs_env)
                  END IF
               END IF

!              Print the requested properties
               CALL print_xas_tdp_to_file(current_state, xas_tdp_env, xas_tdp_control, xas_tdp_section)

!              Free some unneeded attributes of current_state
               CALL free_ds_memory(current_state)

               current_state_index = current_state_index+1
               NULLIFY (current_state)

            END DO ! state type

            CALL free_exat_memory(xas_tdp_env, iatom)
         END DO ! atom of kind
      END DO ! kind

!  Return to ususal KS matrix
      IF (dft_control%do_admm) THEN
         CALL get_qs_env(qs_env, matrix_ks=matrix_ks, admm_env=admm_env)
         DO ispin = 1, dft_control%nspins
            CALL admm_uncorrect_for_eigenvalues(ispin, admm_env, matrix_ks(ispin)%matrix)
         END DO
      END IF

!  Return to initial basis set radii
      IF (xas_tdp_control%eps_pgf > 0.0_dp) THEN
         CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set)
         DO ikind = 1, SIZE(atomic_kind_set)
            CALL get_qs_kind(qs_kind_set(ikind), basis_set=tmp_basis, basis_type="ORB")
            CALL init_interaction_radii_orb_basis(tmp_basis, eps_pgf_orb=dft_control%qs_control%eps_pgf_orb)
         END DO
      END IF

!  Clean-up
      CALL xas_tdp_env_release(xas_tdp_env)

   END SUBROUTINE xas_tdp_core

! **************************************************************************************************
!> \brief Overall control and  environment types initialization
!> \param xas_tdp_env the environment type to initialize
!> \param xas_tdp_control the control type to initialize
!> \param qs_env the inherited qs environement type
! **************************************************************************************************
   SUBROUTINE xas_tdp_init(xas_tdp_env, xas_tdp_control, qs_env)

      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'xas_tdp_init', routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_string_length)               :: kind_name
      INTEGER                                            :: at_ind, i, ispin, j, k, kind_ind, &
                                                            n_donor_states, n_kinds, nao, &
                                                            nat_of_kind, natom, nex_atoms, &
                                                            nex_kinds, nmatch, nspins
      INTEGER, DIMENSION(2)                              :: homo, n_mo, n_moloc
      INTEGER, DIMENSION(:), POINTER                     :: ind_of_kind
      LOGICAL                                            :: do_os, do_uks, unique
      REAL(dp)                                           :: fact
      REAL(dp), DIMENSION(:), POINTER                    :: mo_evals
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: at_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s, rho_ao
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_type), POINTER                  :: tmp_basis
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_loc_env_new_type), POINTER                 :: qs_loc_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_type), POINTER                        :: dummy_section
      TYPE(section_vals_type), POINTER                   :: loc_section, xas_tdp_section

      NULLIFY (xas_tdp_section, at_kind_set, ind_of_kind, dft_control, qs_kind_set, tmp_basis)
      NULLIFY (qs_loc_env, loc_section, mos, particle_set, rho, rho_ao, mo_evals, cell)
      NULLIFY (mo_coeff, matrix_ks, admm_env, dummy_section)

!  XAS TDP control type initialization
      xas_tdp_section => section_vals_get_subs_vals(qs_env%input, "DFT%XAS_TDP")

      CALL get_qs_env(qs_env, dft_control=dft_control)
      CALL xas_tdp_control_create(dft_control%xas_tdp_control)
      CALL read_xas_tdp_control(dft_control%xas_tdp_control, xas_tdp_section)

      xas_tdp_control => dft_control%xas_tdp_control

!  Check the qs_env for a LSD/ROKS calculation
      IF (dft_control%uks) xas_tdp_control%do_uks = .TRUE.
      IF (dft_control%roks) xas_tdp_control%do_roks = .TRUE.
      do_uks = xas_tdp_control%do_uks
      do_os = do_uks .OR. xas_tdp_control%do_roks

!  XAS TDP environment type initialization
      CALL xas_tdp_env_create(xas_tdp_env)

!  Retrieving the excited atoms indices and correspondig state types
      IF (xas_tdp_control%define_excited == xas_tdp_by_index) THEN

!        simply copy indices from xas_tdp_control
         nex_atoms = SIZE(xas_tdp_control%list_ex_atoms)
         CALL set_xas_tdp_env(xas_tdp_env, nex_atoms=nex_atoms)
         ALLOCATE (xas_tdp_env%ex_atom_indices(nex_atoms))
         ALLOCATE (xas_tdp_env%state_types(SIZE(xas_tdp_control%state_types, 1), nex_atoms))
         xas_tdp_env%ex_atom_indices = xas_tdp_control%list_ex_atoms
         xas_tdp_env%state_types = xas_tdp_control%state_types

!        Test that these indices are within the range of available atoms
         CALL get_qs_env(qs_env=qs_env, natom=natom)
         IF (ANY(xas_tdp_env%ex_atom_indices > natom)) THEN
            CPABORT("Invalid index for the ATOM_LIST keyword.")
         END IF

!        Check atom kinds and fill corresponding array
         ALLOCATE (xas_tdp_env%ex_kind_indices(nex_atoms))
         xas_tdp_env%ex_kind_indices = 0
         k = 0
         CALL get_qs_env(qs_env, particle_set=particle_set)
         DO i = 1, nex_atoms
            at_ind = xas_tdp_env%ex_atom_indices(i)
            CALL get_atomic_kind(particle_set(at_ind)%atomic_kind, kind_number=j)
            IF (ALL(ABS(xas_tdp_env%ex_kind_indices-j) .NE. 0)) THEN
               k = k+1
               xas_tdp_env%ex_kind_indices(k) = j
            END IF
         END DO
         nex_kinds = k
         CALL set_xas_tdp_env(xas_tdp_env, nex_kinds=nex_kinds)
         CALL reallocate(xas_tdp_env%ex_kind_indices, 1, nex_kinds)

      ELSE IF (xas_tdp_control%define_excited == xas_tdp_by_kind) THEN

!        need to find out which atom of which kind is excited
         CALL get_qs_env(qs_env=qs_env, atomic_kind_set=at_kind_set)
         n_kinds = SIZE(at_kind_set)
         nex_atoms = 0

         nex_kinds = SIZE(xas_tdp_control%list_ex_kinds)
         ALLOCATE (xas_tdp_env%ex_kind_indices(nex_kinds))
         k = 0

         DO i = 1, n_kinds
            CALL get_atomic_kind(atomic_kind=at_kind_set(i), name=kind_name, &
                                 natom=nat_of_kind, kind_number=kind_ind)
            IF (ANY(xas_tdp_control%list_ex_kinds == kind_name)) THEN
               nex_atoms = nex_atoms+nat_of_kind
               k = k+1
               xas_tdp_env%ex_kind_indices(k) = kind_ind
            END IF
         END DO

         ALLOCATE (xas_tdp_env%ex_atom_indices(nex_atoms))
         ALLOCATE (xas_tdp_env%state_types(SIZE(xas_tdp_control%state_types, 1), nex_atoms))
         nex_atoms = 0
         nmatch = 0

         DO i = 1, n_kinds
            CALL get_atomic_kind(atomic_kind=at_kind_set(i), name=kind_name, &
                                 natom=nat_of_kind, atom_list=ind_of_kind)
            DO j = 1, nex_kinds
               IF (xas_tdp_control%list_ex_kinds(j) == kind_name) THEN
                  xas_tdp_env%ex_atom_indices(nex_atoms+1:nex_atoms+nat_of_kind) = ind_of_kind
                  DO k = 1, SIZE(xas_tdp_control%state_types, 1)
                     xas_tdp_env%state_types(k, nex_atoms+1:nex_atoms+nat_of_kind) = &
                        xas_tdp_control%state_types(k, j)
                  END DO
                  nex_atoms = nex_atoms+nat_of_kind
                  nmatch = nmatch+1
               END IF
            END DO
         END DO

         CALL set_xas_tdp_env(xas_tdp_env, nex_atoms=nex_atoms, nex_kinds=nex_kinds)

!        Verifying that the input was valid
         IF (nmatch .NE. SIZE(xas_tdp_control%list_ex_kinds)) THEN
            CPABORT("Invalid kind(s) for the KIND_LIST keyword.")
         END IF

      END IF

!  Sort the excited atoms indices (for convinience and use of locate function)
      CALL sort_unique(xas_tdp_env%ex_atom_indices, unique)
      IF (.NOT. unique) THEN
         CPABORT("Excited atoms not uniquely defined.")
      END IF

!  Check for periodicity
      CALL get_qs_env(qs_env, cell=cell)
      IF (ALL(cell%perd == 0)) THEN
         xas_tdp_control%is_periodic = .FALSE.
      ELSE IF (ALL(cell%perd == 1)) THEN
         xas_tdp_control%is_periodic = .TRUE.
      ELSE
         CPABORT("XAS TDP only implemented for full PBCs or non-PBCs")
      END IF

!  Allocating memory for the array of donor states
      n_donor_states = COUNT(xas_tdp_env%state_types /= xas_not_excited)
      ALLOCATE (xas_tdp_env%donor_states(n_donor_states))
      DO i = 1, n_donor_states
         CALL donor_state_create(xas_tdp_env%donor_states(i))
      END DO

!  In case of ADMM, for the whole duration of the XAS_TDP, we need the total KS matrix
      IF (dft_control%do_admm) THEN
         CALL get_qs_env(qs_env, admm_env=admm_env, matrix_ks=matrix_ks)

         DO ispin = 1, dft_control%nspins
            CALL admm_correct_for_eigenvalues(ispin, admm_env, matrix_ks(ispin)%matrix)
         END DO
      END IF

!  In case of externally imposed EPS_PGF_XAS, need to update ORB and RI_XAS interaction radii
      IF (xas_tdp_control%eps_pgf > 0.0_dp) THEN
         CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set)

         DO i = 1, SIZE(qs_kind_set)
            CALL get_qs_kind(qs_kind_set(i), basis_set=tmp_basis, basis_type="ORB")
            CALL init_interaction_radii_orb_basis(tmp_basis, eps_pgf_orb=xas_tdp_control%eps_pgf)
            CALL get_qs_kind(qs_kind_set(i), basis_set=tmp_basis, basis_type="RI_XAS")
            CALL init_interaction_radii_orb_basis(tmp_basis, eps_pgf_orb=xas_tdp_control%eps_pgf)
         END DO
      END IF

!  In case of ground state OT optimization, compute the MO eigenvalues and get canonical MOs
      IF (qs_env%scf_env%method == ot_method_nr) THEN

         CALL get_qs_env(qs_env, mos=mos, matrix_ks=matrix_ks)
         nspins = 1; IF (do_uks) nspins = 2

         DO ispin = 1, nspins
            CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff, eigenvalues=mo_evals)
            CALL calculate_subspace_eigenvalues(mo_coeff, matrix_ks(ispin)%matrix, evals_arg=mo_evals)
         END DO
      END IF

!  Initializing the qs_loc_env from the LOCALIZE subsection of XAS_TDP (largely inpired by MI's XAS)
!  We create the LOCALIZE subsection here, since it is completely overwritten anyways
      CALL create_localize_section(dummy_section)
      CALL section_vals_create(xas_tdp_control%loc_subsection, dummy_section)
      CALL section_vals_val_set(xas_tdp_control%loc_subsection, "_SECTION_PARAMETERS_", &
                                l_val=xas_tdp_control%do_loc)
      CALL section_release(dummy_section)
      xas_tdp_control%print_loc_subsection => section_vals_get_subs_vals( &
                                              xas_tdp_control%loc_subsection, "PRINT")

      CALL qs_loc_env_create(xas_tdp_env%qs_loc_env)
      qs_loc_env => xas_tdp_env%qs_loc_env
      loc_section => xas_tdp_control%loc_subsection
!     getting the number of MOs
      CALL get_qs_env(qs_env, mos=mos)
      CALL get_mo_set(mos(1)%mo_set, nmo=n_mo(1), homo=homo(1), nao=nao)
      n_mo(2) = n_mo(1)
      homo(2) = homo(1)
      nspins = 1
      IF (do_os) CALL get_mo_set(mos(2)%mo_set, nmo=n_mo(2), homo=homo(2))
      IF (do_uks) nspins = 2 !in roks, same MOs for both spins

      ! by default, all (doubly occupied) homo are localized
      IF (xas_tdp_control%n_search < 0 .OR. xas_tdp_control%n_search > MINVAL(homo)) &
         xas_tdp_control%n_search = MINVAL(homo)
      CALL qs_loc_control_init(qs_loc_env, loc_section, do_homo=.TRUE., do_xas=.TRUE., &
                               nloc_xas=xas_tdp_control%n_search, spin_xas=1)

      ! do_xas argument above only prepares spin-alpha localization
      IF (do_uks) THEN
         qs_loc_env%localized_wfn_control%nloc_states(2) = xas_tdp_control%n_search
         qs_loc_env%localized_wfn_control%lu_bound_states(1, 2) = 1
         qs_loc_env%localized_wfn_control%lu_bound_states(2, 2) = xas_tdp_control%n_search
      END IF

!     final qs_loc_env initialization. Impose Berry operator
      qs_loc_env%localized_wfn_control%operator_type = op_loc_berry
      IF (.NOT. xas_tdp_control%do_loc) THEN
         qs_loc_env%localized_wfn_control%localization_method = do_loc_none
      ELSE
         n_moloc = qs_loc_env%localized_wfn_control%nloc_states
         CALL set_loc_centers(qs_loc_env%localized_wfn_control, n_moloc, nspins)
         IF (do_uks) THEN
            CALL qs_loc_env_init(qs_loc_env, qs_loc_env%localized_wfn_control, &
                                 qs_env, do_localize=.TRUE.)
         ELSE
            CALL qs_loc_env_init(qs_loc_env, qs_loc_env%localized_wfn_control, &
                                 qs_env, do_localize=.TRUE., myspin=1)
         END IF
      END IF

!  Allocating memory for the array of excited atoms MOs. Worst case senario, all searched MOs are
!  associated to the same atom
      ALLOCATE (xas_tdp_env%mos_of_ex_atoms(xas_tdp_control%n_search, nex_atoms, nspins))

!  Compute the projector on the unoccupied, unperturbed ground state: Q = 1 - SP, sor each spin
      IF (do_os) nspins = 2
      CALL get_qs_env(qs_env, rho=rho, matrix_s=matrix_s)
      CALL qs_rho_get(rho, rho_ao=rho_ao)

      ALLOCATE (xas_tdp_env%q_projector(nspins))
      ALLOCATE (xas_tdp_env%q_projector(1)%matrix)
      CALL dbcsr_create(xas_tdp_env%q_projector(1)%matrix, name="Q PROJECTOR ALPHA", &
                        template=matrix_s(1)%matrix, matrix_type="N")
      IF (do_os) THEN
         ALLOCATE (xas_tdp_env%q_projector(2)%matrix)
         CALL dbcsr_create(xas_tdp_env%q_projector(2)%matrix, name="Q PROJECTOR BETA", &
                           template=matrix_s(1)%matrix, matrix_type="N")
      END IF

!     In the case of spin-restricted calculations, rho_ao includes double occupency => 0.5 prefactor
      fact = -0.5_dp; IF (do_os) fact = -1.0_dp
      CALL dbcsr_multiply('N', 'N', fact, matrix_s(1)%matrix, rho_ao(1)%matrix, 0.0_dp, &
                          xas_tdp_env%q_projector(1)%matrix)
      CALL dbcsr_add_on_diag(xas_tdp_env%q_projector(1)%matrix, 1.0_dp)
      CALL dbcsr_finalize(xas_tdp_env%q_projector(1)%matrix)

      IF (do_os) THEN
         CALL dbcsr_multiply('N', 'N', fact, matrix_s(1)%matrix, rho_ao(2)%matrix, 0.0_dp, &
                             xas_tdp_env%q_projector(2)%matrix)
         CALL dbcsr_add_on_diag(xas_tdp_env%q_projector(2)%matrix, 1.0_dp)
         CALL dbcsr_finalize(xas_tdp_env%q_projector(2)%matrix)
      END IF

!  Create the structure for the dipole in the AO basis
      ALLOCATE (xas_tdp_env%dipmat(3))
      DO i = 1, 3
         ALLOCATE (xas_tdp_env%dipmat(i)%matrix)
         CALL dbcsr_copy(xas_tdp_env%dipmat(i)%matrix, matrix_s(1)%matrix, name="XAS TDP dipole matrix")
         CALL dbcsr_set(xas_tdp_env%dipmat(i)%matrix, 0.0_dp)
      END DO

!  Create the structure for the electric quadrupole in the AO basis, if required
      IF (xas_tdp_control%do_quad) THEN
         ALLOCATE (xas_tdp_env%quadmat(6))
         DO i = 1, 6
            ALLOCATE (xas_tdp_env%quadmat(i)%matrix)
            CALL dbcsr_copy(xas_tdp_env%quadmat(i)%matrix, matrix_s(1)%matrix, name="XAS TDP quadrupole matrix")
            CALL dbcsr_set(xas_tdp_env%quadmat(i)%matrix, 0.0_dp)
         END DO
      END IF

!     Precompute it in the velocity representation, if so chosen
      IF (xas_tdp_control%dipole_form == xas_dip_vel) THEN
         !enforce minimum image to avoid any PBCs related issues. Ok because very localized densities
         CALL p_xyz_ao(xas_tdp_env%dipmat, qs_env, minimum_image=.TRUE.)
      END IF

!  Allocate SOC in AO basis matrices
      IF (xas_tdp_control%do_soc) THEN
         ALLOCATE (xas_tdp_env%orb_soc(3))
         DO i = 1, 3
            ALLOCATE (xas_tdp_env%orb_soc(i)%matrix)
         END DO
      END IF

!  Check that everything is allowed
      CALL safety_check(xas_tdp_control)

!  Initialize libint for the 3-center integrals
      CALL cp_libint_static_init()

   END SUBROUTINE xas_tdp_init

! **************************************************************************************************
!> \brief Checks for forbidden keywords combinations
!> \param xas_tdp_control ...
! **************************************************************************************************
   SUBROUTINE safety_check(xas_tdp_control)

      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control

      CHARACTER(len=*), PARAMETER :: routineN = 'safety_check', routineP = moduleN//':'//routineN

      !PB only available without exact exchange
      IF (xas_tdp_control%is_periodic .AND. xas_tdp_control%do_hfx &
          .AND. xas_tdp_control%x_op == do_potential_coulomb) THEN
         CPABORT("XAS TDP with Coulomb operator for exact exchange only supports non-periodic BCs")
      END IF

      !open-shell/closed-shell tests
      IF (xas_tdp_control%do_roks .OR. xas_tdp_control%do_uks) THEN

         IF (.NOT. (xas_tdp_control%do_spin_cons .OR. xas_tdp_control%do_spin_flip)) THEN
            CPABORT("Need spin-conserving and/or spin-flip excitations for open-shell systems")
         END IF

         IF (xas_tdp_control%do_singlet .OR. xas_tdp_control%do_triplet) THEN
            CPABORT("Singlet/triplet excitations only for restricted closed-shell systems")
         END IF

         IF (xas_tdp_control%do_soc .AND. .NOT. &
             (xas_tdp_control%do_spin_flip .AND. xas_tdp_control%do_spin_cons)) THEN

            CPABORT("Both spin-conserving and spin-flip excitations are required for SOC")
         END IF
      ELSE

         IF (.NOT. (xas_tdp_control%do_singlet .OR. xas_tdp_control%do_triplet)) THEN
            CPABORT("Need singlet and/or triplet excitations for closed-shell systems")
         END IF

         IF (xas_tdp_control%do_spin_cons .OR. xas_tdp_control%do_spin_flip) THEN
            CPABORT("Spin-conserving/spin-flip excitations only for open-shell systems")
         END IF

         IF (xas_tdp_control%do_soc .AND. .NOT. &
             (xas_tdp_control%do_singlet .AND. xas_tdp_control%do_triplet)) THEN

            CPABORT("Both singlet and triplet excitations are needed for SOC")
         END IF
      END IF

      !Look for eigenvalues within energy range and SOC not compatible
      IF (xas_tdp_control%e_range .GT. 0.0_dp .AND. xas_tdp_control%do_soc) THEN
         CPABORT("Use N_EXCITED instead of ENERGY_RANGE for SOC calculations")
      END IF

      !TDA, full-TDDFT and diagonalization
      IF (.NOT. xas_tdp_control%tamm_dancoff) THEN

         IF (xas_tdp_control%do_spin_flip) THEN
            CPABORT("Spin-flip kernel only implemented for Tamm-Dancoff approximation")
         END IF

      END IF

   END SUBROUTINE safety_check

! **************************************************************************************************
!> \brief Prints some basic info about the chosen parameters
!> \param ou the output unis
!> \param xas_tdp_control ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE print_info(ou, xas_tdp_control, qs_env)

      INTEGER, INTENT(IN)                                :: ou
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'print_info', routineP = moduleN//':'//routineN

      CHARACTER(LEN=2*default_string_length)             :: shortform
      CHARACTER(LEN=20)                                  :: tmpStr
      CHARACTER(LEN=20*default_string_length)            :: reference
      INTEGER                                            :: i, il, irep, nrep
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(section_vals_type), POINTER                   :: input, libxc_fun, xc_fun, xc_functionals

      NULLIFY (input, xc_functionals, libxc_fun, xc_fun, dft_control)

      IF (ou .LE. 0) RETURN

      CALL get_qs_env(qs_env, input=input, dft_control=dft_control)

      !Reference calculation
      IF (xas_tdp_control%do_uks) THEN
         WRITE (UNIT=ou, FMT="(/,T3,A)") &
            "XAS_TDP| Reference calculation: Unrestricted Kohn-Sham"
      ELSE IF (xas_tdp_control%do_roks) THEN
         WRITE (UNIT=ou, FMT="(/,T3,A)") &
            "XAS_TDP| Reference calculation: Restricted Open-Shell Kohn-Sham"
      ELSE
         WRITE (UNIT=ou, FMT="(/,T3,A)") &
            "XAS_TDP| Reference calculation: Restricted Closed-Shell Kohn-Sham"
      END IF

      !TDA
      IF (xas_tdp_control%tamm_dancoff) THEN
         WRITE (UNIT=ou, FMT="(T3,A)") &
            "XAS_TDP| Tamm-Dancoff Approximation (TDA): On"
      ELSE
         WRITE (UNIT=ou, FMT="(T3,A)") &
            "XAS_TDP| Tamm-Dancoff Approximation (TDA): Off"
      END IF

      !Dipole form
      IF (xas_tdp_control%dipole_form == xas_dip_vel) THEN
         WRITE (UNIT=ou, FMT="(T3,A)") &
            "XAS_TDP| Transition Dipole Representation: VELOCITY"
      ELSE
         WRITE (UNIT=ou, FMT="(T3,A)") &
            "XAS_TDP| Transition Dipole Representation: LENGTH"
      END IF

      !Quadrupole
      IF (xas_tdp_control%do_quad) THEN
         WRITE (UNIT=ou, FMT="(T3,A)") &
            "XAS_TDP| Transition Quadrupole: On"
      END IF

      !EPS_PGF
      IF (xas_tdp_control%eps_pgf > 0.0_dp) THEN
         WRITE (UNIT=ou, FMT="(T3,A,ES7.1)") &
            "XAS_TDP| EPS_PGF_XAS: ", xas_tdp_control%eps_pgf
      ELSE
         WRITE (UNIT=ou, FMT="(T3,A,ES7.1,A)") &
            "XAS_TDP| EPS_PGF_XAS: ", dft_control%qs_control%eps_pgf_orb, " (= EPS_PGF_ORB)"
      END IF

      !Grid info
      IF (xas_tdp_control%do_soc .OR. xas_tdp_control%do_xc) THEN
         WRITE (UNIT=ou, FMT="(T3,A)") &
            "XAS_TDP| Radial Grid(s) Info: Kind,  na,  nr"
         DO i = 1, SIZE(xas_tdp_control%grid_info, 1)
            WRITE (UNIT=ou, FMT="(T3,A,A6,A,A,A,A)") &
               "                            ", TRIM(xas_tdp_control%grid_info(i, 1)), ", ", &
               TRIM(xas_tdp_control%grid_info(i, 2)), ", ", TRIM(xas_tdp_control%grid_info(i, 3))
         END DO
      END IF

      !No kernel
      IF (.NOT. xas_tdp_control%do_coulomb) THEN
         WRITE (UNIT=ou, FMT="(/,T3,A)") &
            "XAS_TDP| No kernel (standard DFT)"
      END IF

      !XC kernel
      IF (xas_tdp_control%do_xc) THEN

         WRITE (UNIT=ou, FMT="(/,T3,A,F5.2,A)") &
            "XAS_TDP| RI Region's Radius: ", xas_tdp_control%ri_radius*angstrom, " Ang"

         !parallelization scheme
         WRITE (UNIT=ou, FMT="(T3,A,I4,A)") &
            "XAS_TDP| Excited atoms distributed on batches of ", xas_tdp_control%batch_size, " procs"

         CALL libxc_version_info(tmpStr)
         WRITE (UNIT=ou, FMT="(T3,A,A,/,T3,A,/)") &
            "XAS_TDP| XC Kernel Functional(s) from LIBXC Vers. ", TRIM(tmpStr(1:5)), &
            "         (Marques, Oliveira, Burnus, CPC 183, 2272 (2012)):"
         xc_functionals => section_vals_get_subs_vals(input, "DFT%XAS_TDP%KERNEL%XC_FUNCTIONAL")
         xc_fun => section_vals_get_subs_vals2(xc_functionals, i_section=1)
         CALL section_vals_get(xc_fun, n_repetition=nrep)
         DO irep = 1, nrep
            NULLIFY (libxc_fun)
            CALL section_vals_duplicate(xc_fun, libxc_fun, i_rep_start=irep, i_rep_end=irep)
            IF (.NOT. ASSOCIATED(libxc_fun)) EXIT
            CALL xc_functional_get_info(libxc_fun, lsd=.TRUE., reference=reference, shortform=shortform)
            CALL section_vals_release(libxc_fun)
            WRITE (ou, FMT="('  FUNCTIONAL| ',A,':')") TRIM(shortform)
            DO il = 1, LEN_TRIM(reference), 67
               WRITE (ou, FMT="('  FUNCTIONAL| ',a67)") reference(il:)
            END DO
         END DO
      END IF

      !HFX kernel
      IF (xas_tdp_control%do_hfx) THEN
         WRITE (UNIT=ou, FMT="(/,T3,A,/,/,T3,A,F5.3)") &
            "XAS_TDP| Exat Exchange Kernel: On ", &
            "EXACT_EXCHANGE| Scale: ", xas_tdp_control%sx
         IF (xas_tdp_control%x_op == do_potential_coulomb) THEN
            WRITE (UNIT=ou, FMT="(T3,A)") &
               "EXACT_EXCHANGE| Operator : Coulomb"
         ELSE IF (xas_tdp_control%x_op == do_potential_truncated) THEN
            WRITE (UNIT=ou, FMT="(T3,A,/,T3,A,F5.2,A,/,T3,A,A)") &
               "EXACT_EXCHANGE| Operator: Truncated Coulomb", &
               "EXACT_EXCHANGE| Range: ", xas_tdp_control%x_range*angstrom, ", (Ang)", &
               "EXACT_EXCHANGE| T_C_G_DATA: ", xas_tdp_control%x_t_c_filename
         ELSE IF (xas_tdp_control%x_op == do_potential_short) THEN
            WRITE (UNIT=ou, FMT="(T3,A,/,T3,A,F5.2,A,/,T3,A,F5.2,A,/,T3,A,ES7.1)") &
               "EXACT_EXCHANGE| Operator: Short Range", &
               "EXACT_EXCHANGE| Omega: ", xas_tdp_control%x_omega, ", (1/a0)", &
               "EXACT_EXCHANGE| Effective Range: ", xas_tdp_control%x_range*angstrom, ", (Ang)", &
               "EXACT_EXCHANGE| EPS_RANGE: ", xas_tdp_control%eps_range
         END IF
         IF (xas_tdp_control%eps_screen > 1.0E-16) THEN
            WRITE (UNIT=ou, FMT="(T3,A,ES7.1)") &
               "EXACT_EXCHANGE| EPS_SCREENING: ", xas_tdp_control%eps_screen
         END IF
      END IF

   END SUBROUTINE print_info

! **************************************************************************************************
!> \brief Assosciate (possibly localized) lowest energy  MOs to each excited atoms. The procedure
!>        looks for MOs "centered" on the excited atoms by comparing distances. It
!>        then fills the mos_of_ex_atoms arrays of the xas_tdp_env. Only the xas_tdp_control%n_search
!>        lowest energy MOs are considered. Largely inspired by MI's implementation of XAS
!>        It is assumed that the Berry phase is used to compute centers.
!> \param xas_tdp_env ...
!> \param xas_tdp_control ...
!> \param qs_env ...
!> \note Whether localization took place or not, the procedure is the same as centers are stored in
!>       xas_tdp_env%qs_loc_env%localized_wfn_control%centers_set
!>       Assumes that find_mo_centers has been run previously
! **************************************************************************************************
   SUBROUTINE assign_mos_to_ex_atoms(xas_tdp_env, xas_tdp_control, qs_env)

      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'assign_mos_to_ex_atoms', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: at_index, iat, iat_memo, imo, ispin, &
                                                            n_atoms, n_search, nex_atoms, nspins
      INTEGER, DIMENSION(3)                              :: perd_init
      INTEGER, DIMENSION(:, :, :), POINTER               :: mos_of_ex_atoms
      REAL(dp)                                           :: dist, dist_min
      REAL(dp), DIMENSION(3)                             :: at_pos, r_ac, wfn_center
      TYPE(cell_type), POINTER                           :: cell
      TYPE(localized_wfn_control_type), POINTER          :: localized_wfn_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      NULLIFY (localized_wfn_control, mos_of_ex_atoms, cell, particle_set)

!  Initialization. mos_of_ex_atoms filled with -1, meaning no assigned state
      mos_of_ex_atoms => xas_tdp_env%mos_of_ex_atoms
      mos_of_ex_atoms(:, :, :) = -1
      n_search = xas_tdp_control%n_search
      nex_atoms = xas_tdp_env%nex_atoms
      localized_wfn_control => xas_tdp_env%qs_loc_env%localized_wfn_control
      CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, cell=cell)
      n_atoms = SIZE(particle_set)
      nspins = 1; IF (xas_tdp_control%do_uks) nspins = 2

!     Temporarly impose periodic BCs because of Berry's phase operator used for localization
      perd_init = cell%perd
      cell%perd = 1

!  Loop over n_search lowest energy MOs and all atoms, for each spin
      DO ispin = 1, nspins
         DO imo = 1, n_search
!           retrieve MO wave function center coordinates.
            wfn_center(1:3) = localized_wfn_control%centers_set(ispin)%array(1:3, imo)
            iat_memo = 0

!           a large enough value to avoid bad surprises
            dist_min = 10000.0_dp
            DO iat = 1, n_atoms
               at_pos = particle_set(iat)%r
               r_ac = pbc(at_pos, wfn_center, cell)
               dist = NORM2(r_ac)

!              keep memory of which atom is the closest to the wave function center
               IF (dist < dist_min) THEN
                  iat_memo = iat
                  dist_min = dist
               END IF
            END DO

!           Verify that the closest atom is actually excited and assign the MO if so
            IF (ANY(xas_tdp_env%ex_atom_indices == iat_memo)) THEN
               at_index = locate(xas_tdp_env%ex_atom_indices, iat_memo)
               mos_of_ex_atoms(imo, at_index, ispin) = 1
            END IF
         END DO !imo
      END DO !ispin

!  Go back to initial BCs
      cell%perd = perd_init

   END SUBROUTINE assign_mos_to_ex_atoms

! **************************************************************************************************
!> \brief Re-initialize the qs_loc_env to the current MOs.
!> \param qs_loc_env the env to re-initialize
!> \param n_loc_states the number of states to include
!> \param do_uks in cas of spin unrestricted calculation, initialize for both spins
!> \param qs_env ...
!> \note  Useful when one needs to make use of qs_loc features and it is either with canonical MOs
!>        or the localized MOs have been modified. do_localize is overwritten.
!>        Same loc range for both spins
! **************************************************************************************************
   SUBROUTINE reinit_qs_loc_env(qs_loc_env, n_loc_states, do_uks, qs_env)

      TYPE(qs_loc_env_new_type), POINTER                 :: qs_loc_env
      INTEGER, INTENT(IN)                                :: n_loc_states
      LOGICAL, INTENT(IN)                                :: do_uks
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'reinit_qs_loc_env', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, nspins
      TYPE(localized_wfn_control_type), POINTER          :: loc_wfn_control

!  First, release the old env
      CALL qs_loc_env_release(qs_loc_env)

!  Re-create it
      CALL qs_loc_env_create(qs_loc_env)
      CALL localized_wfn_control_create(qs_loc_env%localized_wfn_control)
      loc_wfn_control => qs_loc_env%localized_wfn_control

!  Initialize it
      loc_wfn_control%localization_method = do_loc_none
      loc_wfn_control%operator_type = op_loc_berry
      loc_wfn_control%nloc_states(:) = n_loc_states
      loc_wfn_control%eps_occ = 0.0_dp
      loc_wfn_control%lu_bound_states(1, :) = 1
      loc_wfn_control%lu_bound_states(2, :) = n_loc_states
      loc_wfn_control%set_of_states = state_loc_list
      loc_wfn_control%do_homo = .TRUE.
      ALLOCATE (loc_wfn_control%loc_states(n_loc_states, 2))
      DO i = 1, n_loc_states
         loc_wfn_control%loc_states(i, :) = i
      END DO

      nspins = 1; IF (do_uks) nspins = 2
      CALL set_loc_centers(loc_wfn_control, loc_wfn_control%nloc_states, nspins=nspins)
      ! need to set do_localize=.TRUE. because otherwise no routine works
      IF (do_uks) THEN
         CALL qs_loc_env_init(qs_loc_env, loc_wfn_control, qs_env, do_localize=.TRUE.)
      ELSE
         CALL qs_loc_env_init(qs_loc_env, loc_wfn_control, qs_env, myspin=1, do_localize=.TRUE.)
      END IF

   END SUBROUTINE reinit_qs_loc_env

! *************************************************************************************************
!> \brief Diagonalize the subset of previously localized MOs that are associated to each excited
!>        atoms. Updates the MO coeffs accordingly.
!> \param xas_tdp_env ...
!> \param xas_tdp_control ...
!> \param qs_env ...
!> \note  Needed because after localization, the MOs loose their identity (1s, 2s , 2p, etc)
! **************************************************************************************************
   SUBROUTINE diagonalize_assigned_mo_subset(xas_tdp_env, xas_tdp_control, qs_env)

      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'diagonalize_assigned_mo_subset', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, iat, ilmo, ispin, nao, nlmo, nspins
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: evals
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: ks_struct, lmo_struct
      TYPE(cp_fm_type), POINTER                          :: evecs, ks_fm, lmo_fm, mo_coeff, work
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos

      NULLIFY (mos, mo_coeff, matrix_ks, para_env, blacs_env, lmo_struct, lmo_fm, ks_struct, ks_fm)
      NULLIFY (evecs, work)

      ! Get what we need from qs_env
      CALL get_qs_env(qs_env, mos=mos, matrix_ks=matrix_ks, para_env=para_env, blacs_env=blacs_env)

      nspins = 1; IF (xas_tdp_control%do_uks) nspins = 2

      ! Loop over the excited atoms and spin
      DO ispin = 1, nspins
         DO iat = 1, xas_tdp_env%nex_atoms

            ! get the MOs
            CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff, nao=nao)

            ! count how many MOs are associated to this atom and create a fm/struct
            nlmo = COUNT(xas_tdp_env%mos_of_ex_atoms(:, iat, ispin) == 1)
            CALL cp_fm_struct_create(lmo_struct, nrow_global=nao, ncol_global=nlmo, &
                                     para_env=para_env, context=blacs_env)
            CALL cp_fm_create(lmo_fm, lmo_struct)
            CALL cp_fm_create(work, lmo_struct)

            CALL cp_fm_struct_create(ks_struct, nrow_global=nlmo, ncol_global=nlmo, &
                                     para_env=para_env, context=blacs_env)
            CALL cp_fm_create(ks_fm, ks_struct)
            CALL cp_fm_create(evecs, ks_struct)

            ! Loop over the localized MOs associated to this atom
            i = 0
            DO ilmo = 1, xas_tdp_control%n_search
               IF (xas_tdp_env%mos_of_ex_atoms(ilmo, iat, ispin) == -1) CYCLE

               i = i+1
               ! put the coeff in our atom-restricted lmo_fm
               CALL cp_fm_to_fm_submat(mo_coeff, lmo_fm, nrow=nao, ncol=1, s_firstrow=1, &
                                       s_firstcol=ilmo, t_firstrow=1, t_firstcol=i)

            END DO !ilmo

            ! Computing the KS matrix in the subset of MOs
            CALL cp_dbcsr_sm_fm_multiply(matrix_ks(ispin)%matrix, lmo_fm, work, ncol=nlmo)
            CALL cp_gemm('T', 'N', nlmo, nlmo, nao, 1.0_dp, lmo_fm, work, 0.0_dp, ks_fm)

            ! Diagonalizing the KS matrix in the subset of MOs
            ALLOCATE (evals(nlmo))
            CALL cp_fm_syevd(ks_fm, evecs, evals)
            DEALLOCATE (evals)

            ! Express the MOs in the basis that diagonalizes KS
            CALL cp_gemm('N', 'N', nao, nlmo, nlmo, 1.0_dp, lmo_fm, evecs, 0.0_dp, work)

            ! Replacing the new MOs back in the MO coeffs
            i = 0
            DO ilmo = 1, xas_tdp_control%n_search
               IF (xas_tdp_env%mos_of_ex_atoms(ilmo, iat, ispin) == -1) CYCLE

               i = i+1
               CALL cp_fm_to_fm_submat(work, mo_coeff, nrow=nao, ncol=1, s_firstrow=1, &
                                       s_firstcol=i, t_firstrow=1, t_firstcol=ilmo)

            END DO

            ! Excited atom clean-up
            CALL cp_fm_release(lmo_fm)
            CALL cp_fm_release(work)
            CALL cp_fm_struct_release(lmo_struct)
            CALL cp_fm_release(ks_fm)
            CALL cp_fm_release(evecs)
            CALL cp_fm_struct_release(ks_struct)
         END DO !iat
      END DO !ispin

   END SUBROUTINE diagonalize_assigned_mo_subset

! **************************************************************************************************
!> \brief Assign core MO(s) to a given donor_state, taking the type (1S, 2S, etc) into account.
!>        The projection on a representative Slater-type orbital basis is used as a indicator.
!>        It is assumed that MOs are already assigned to excited atoms based on their center
!> \param donor_state the donor_state to which a MO must be assigned
!> \param xas_tdp_env ...
!> \param xas_tdp_control ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE assign_mos_to_donor_state(donor_state, xas_tdp_env, xas_tdp_control, qs_env)

      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'assign_mos_to_donor_state', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: at_index, i, iat, imo, ispin, l, my_mo, &
                                                            n_search, n_states, nao, ndo_so, nj, &
                                                            nsgf_kind, nsgf_sto, nspins, &
                                                            output_unit, zval
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: my_mos
      INTEGER, DIMENSION(2)                              :: next_best_overlap_ind
      INTEGER, DIMENSION(4, 7)                           :: ne
      INTEGER, DIMENSION(:), POINTER                     :: first_sgf, lq, nq
      INTEGER, DIMENSION(:, :, :), POINTER               :: mos_of_ex_atoms
      REAL(dp)                                           :: zeff
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: diag, overlap, sto_overlap
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: max_overlap
      REAL(dp), DIMENSION(2)                             :: next_best_overlap
      REAL(dp), DIMENSION(:), POINTER                    :: mo_evals, zeta
      REAL(dp), DIMENSION(:, :), POINTER                 :: overlap_matrix, tmp_coeff
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: eval_mat_struct, gs_struct
      TYPE(cp_fm_type), POINTER                          :: eval_mat, gs_coeffs, mo_coeff, work_mat
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks
      TYPE(gto_basis_set_type), POINTER                  :: kind_basis_set, sto_to_gto_basis_set
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(sto_basis_set_type), POINTER                  :: sto_basis_set

      NULLIFY (sto_basis_set, sto_to_gto_basis_set, qs_kind_set, kind_basis_set, lq, nq, zeta)
      NULLIFY (overlap_matrix, mos, mo_coeff, mos_of_ex_atoms, tmp_coeff, first_sgf, particle_set)
      NULLIFY (mo_evals, matrix_ks, para_env, blacs_env, work_mat)
      NULLIFY (eval_mat_struct, eval_mat, gs_struct, gs_coeffs)

      output_unit = cp_logger_get_default_io_unit()

      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set, mos=mos, particle_set=particle_set, &
                      matrix_ks=matrix_ks, para_env=para_env, blacs_env=blacs_env)

      nspins = 1; IF (xas_tdp_control%do_uks) nspins = 2

!  Contruction of a STO that fits the type of orbital we look for
      ALLOCATE (zeta(1))
      ALLOCATE (lq(1))
      ALLOCATE (nq(1))
!     Retrieving quantum numbers
      IF (donor_state%state_type == xas_1s_type) THEN
         nq(1) = 1
         lq(1) = 0
         n_states = 1
      ELSE IF (donor_state%state_type == xas_2s_type) THEN
         nq(1) = 2
         lq(1) = 0
         n_states = 1
      ELSE IF (donor_state%state_type == xas_2p_type) THEN
         nq(1) = 2
         lq(1) = 1
         n_states = 3
      ELSE
         CPABORT("Procedure for required type not implemented")
      END IF
      ALLOCATE (my_mos(n_states, nspins))
      ALLOCATE (max_overlap(n_states, nspins))

!     Getting the atomic number
      CALL get_qs_kind(qs_kind_set(donor_state%kind_index), zeff=zeff)
      zval = INT(zeff)

!     Electronic configuration (copied from MI's XAS)
      ne = 0
      DO l = 1, 4
         nj = 2*(l-1)+1
         DO i = l, 7
            ne(l, i) = ptable(zval)%e_conv(l-1)-2*nj*(i-l)
            ne(l, i) = MAX(ne(l, i), 0)
            ne(l, i) = MIN(ne(l, i), 2*nj)
         END DO
      END DO

!     computing zeta with the Slater sum rules
      zeta(1) = srules(zval, ne, nq(1), lq(1))

!     Allocating memory and initiate STO
      CALL allocate_sto_basis_set(sto_basis_set)
      CALL set_sto_basis_set(sto_basis_set, nshell=1, nq=nq, lq=lq, zet=zeta)

!     Some clean-up
      DEALLOCATE (nq, lq, zeta)

!  Expanding the STO into (normalized) GTOs for later calculations, use standard 3 gaussians
      CALL create_gto_from_sto_basis(sto_basis_set=sto_basis_set, &
                                     gto_basis_set=sto_to_gto_basis_set, &
                                     ngauss=3)
      sto_to_gto_basis_set%norm_type = 2
      CALL init_orb_basis_set(sto_to_gto_basis_set)

!  Retrieving the atomic kind related GTO in which MOs are expanded
      CALL get_qs_kind(qs_kind_set(donor_state%kind_index), basis_set=kind_basis_set)

!  Allocating and computing the overlap between the two basis (they share the same center)
      CALL get_gto_basis_set(gto_basis_set=kind_basis_set, nsgf=nsgf_kind)
      CALL get_gto_basis_set(gto_basis_set=sto_to_gto_basis_set, nsgf=nsgf_sto)
      ALLOCATE (overlap_matrix(nsgf_sto, nsgf_kind))

!     Making use of MI's subroutine
      CALL calc_stogto_overlap(sto_to_gto_basis_set, kind_basis_set, overlap_matrix)

!     Some clean-up
      CALL deallocate_sto_basis_set(sto_basis_set)
      CALL deallocate_gto_basis_set(sto_to_gto_basis_set)

!  Looping over the potential donor states to compute overlap with STO basis
      mos_of_ex_atoms => xas_tdp_env%mos_of_ex_atoms
      n_search = xas_tdp_control%n_search
      at_index = donor_state%at_index
      iat = locate(xas_tdp_env%ex_atom_indices, at_index)
      ALLOCATE (first_sgf(SIZE(particle_set))) !probably do not need that
      CALL get_particle_set(particle_set=particle_set, qs_kind_set=qs_kind_set, first_sgf=first_sgf)
      ALLOCATE (tmp_coeff(nsgf_kind, 1))
      ALLOCATE (sto_overlap(nsgf_kind))
      ALLOCATE (overlap(n_search))

      next_best_overlap = 0.0_dp
      max_overlap = 0.0_dp

      DO ispin = 1, nspins

         CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff, nao=nao)
         overlap = 0.0_dp

         my_mo = 0
         DO imo = 1, n_search
            IF (mos_of_ex_atoms(imo, iat, ispin) > 0) THEN

               sto_overlap = 0.0_dp
               tmp_coeff = 0.0_dp

!              Getting the relevant coefficients for the candidate state
               CALL cp_fm_get_submatrix(fm=mo_coeff, target_m=tmp_coeff, start_row=first_sgf(at_index), &
                                        start_col=imo, n_rows=nsgf_kind, n_cols=1, transpose=.FALSE.)

!              Computing the product overlap_matrix*coeffs
               CALL dgemm('N', 'N', nsgf_sto, 1, nsgf_kind, 1.0_dp, overlap_matrix, nsgf_sto, &
                          tmp_coeff, nsgf_kind, 0.0_dp, sto_overlap, nsgf_sto)

!              Each element of column vector sto_overlap is the overlap of a basis element of the
!              generated STO basis with the kind specific orbital basis. Take the sum of the absolute
!              values so that rotation (of the px, py, pz for example) does not hinder our search
               overlap(imo) = SUM(ABS(sto_overlap))

            END IF
         END DO

!     Finding the best overlap(s)
         DO i = 1, n_states
            my_mo = MAXLOC(overlap, 1)
            my_mos(i, ispin) = my_mo
            max_overlap(i, ispin) = MAXVAL(overlap, 1)
            overlap(my_mo) = 0.0_dp
         END DO
!        Getting the next best overlap (for validation purposes)
         next_best_overlap(ispin) = MAXVAL(overlap, 1)
         next_best_overlap_ind(ispin) = MAXLOC(overlap, 1)

      END DO !ispin

!     Some clean-up
      DEALLOCATE (overlap_matrix, tmp_coeff)

!  Dealing with the result
      IF (ALL(my_mos > 0) .AND. ALL(my_mos <= n_search)) THEN
!        Assigning the MO indices to the donor_state
         ALLOCATE (donor_state%mo_indices(n_states, nspins))
         donor_state%mo_indices = my_mos
         donor_state%ndo_mo = n_states

!        Storing the MOs in the donor_state, as vectors column: first columns alpha spin, then beta
         CALL cp_fm_struct_create(gs_struct, nrow_global=nao, ncol_global=n_states*nspins, &
                                  para_env=para_env, context=blacs_env)
         CALL cp_fm_create(donor_state%gs_coeffs, gs_struct)

         DO ispin = 1, nspins
            CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff)
            DO i = 1, n_states
               CALL cp_fm_to_fm_submat(msource=mo_coeff, mtarget=donor_state%gs_coeffs, nrow=nao, &
                                       ncol=1, s_firstrow=1, s_firstcol=my_mos(i, ispin), &
                                       t_firstrow=1, t_firstcol=(ispin-1)*n_states+i)
            END DO
         END DO
         gs_coeffs => donor_state%gs_coeffs

         !Keep the subset of the coeffs centered on the excited atom as global array (used a lot)
         ALLOCATE (donor_state%contract_coeffs(nsgf_kind, n_states*nspins))
         CALL cp_fm_get_submatrix(gs_coeffs, donor_state%contract_coeffs, start_row=first_sgf(at_index), &
                                  start_col=1, n_rows=nsgf_kind, n_cols=n_states*nspins)

!     Assigning corresponding energy eigenvalues and writing some info in standard input file

         !standard eigenvalues as gotten from the KS diagonalization in the ground state
         IF (.NOT. xas_tdp_control%do_loc .AND. .NOT. xas_tdp_control%do_roks) THEN
            IF (output_unit > 0) THEN
               WRITE (UNIT=output_unit, FMT="(T5,A,/,T5,A,/,T5,A)") &
                  "The following canonical MO(s) have been associated with the donor state(s)", &
                  "based on the overlap with the components of a minimal STO basis: ", &
                  "                                         Spin   MO index     overlap(sum)"
            END IF

            ALLOCATE (donor_state%energy_evals(n_states, nspins))
            donor_state%energy_evals = 0.0_dp

!           Canonical MO, no change in eigenvalues, only diagonal elements
            DO ispin = 1, nspins
               CALL get_mo_set(mos(ispin)%mo_set, eigenvalues=mo_evals)
               DO i = 1, n_states
                  donor_state%energy_evals(i, ispin) = mo_evals(my_mos(i, ispin))

                  IF (output_unit > 0) THEN
                     WRITE (UNIT=output_unit, FMT="(T46,I4,I11,F17.5)") &
                        ispin, my_mos(i, ispin), max_overlap(i, ispin)
                  END IF
               END DO
            END DO

            !either localization of MOs or ROKS, in both cases the MO eigenvalues from the KS
            !digonalization mat have changed
         ELSE
            IF (output_unit > 0) THEN
               WRITE (UNIT=output_unit, FMT="(T5,A,/,T5,A,/,T5,A)") &
                  "The following localized MO(s) have been associated with the donor state(s)", &
                  "based on the overlap with the components of a minimal STO basis: ", &
                  "                                         Spin   MO index     overlap(sum)"
            END IF

!           Loop over the donor states  and print
            DO ispin = 1, nspins
               DO i = 1, n_states

!                 Print info
                  IF (output_unit > 0) THEN
                     WRITE (UNIT=output_unit, FMT="(T46,I4,I11,F17.5)") &
                        ispin, my_mos(i, ispin), max_overlap(i, ispin)
                  END IF
               END DO
            END DO

!           MO have been rotated or non-physical ROKS MO eigrenvalues:
!           => need epsilon_ij = <psi_i|F|psi_j> = sum_{pq} c_{qi}c_{pj} F_{pq}
!           Note: only have digonal elements by construction
            ndo_so = nspins*n_states
            CALL cp_fm_create(work_mat, gs_struct)
            CALL cp_fm_struct_create(eval_mat_struct, nrow_global=ndo_so, ncol_global=ndo_so, &
                                     para_env=para_env, context=blacs_env)
            CALL cp_fm_create(eval_mat, eval_mat_struct)
            ALLOCATE (diag(ndo_so))

            IF (.NOT. xas_tdp_control%do_roks) THEN

               ALLOCATE (donor_state%energy_evals(n_states, nspins))
               donor_state%energy_evals = 0.0_dp

!              Compute gs_coeff^T * matrix_ks * gs_coeff to get the epsilon_ij matrix
               DO ispin = 1, nspins
                  CALL cp_dbcsr_sm_fm_multiply(matrix_ks(ispin)%matrix, gs_coeffs, work_mat, ncol=ndo_so)
                  CALL cp_gemm('T', 'N', ndo_so, ndo_so, nao, 1.0_dp, gs_coeffs, work_mat, 0.0_dp, eval_mat)

!                 Put the epsilon_ii into the donor_state. No off-diagonal element because of subset diag
                  CALL cp_fm_get_diag(eval_mat, diag)
                  donor_state%energy_evals(:, ispin) = diag((ispin-1)*n_states+1:ispin*n_states)

               END DO

            ELSE
               ! If ROKS, slightly different procedure => 2 KS matrices but one type of MOs
               ALLOCATE (donor_state%energy_evals(n_states, 2))
               donor_state%energy_evals = 0.0_dp

!              Compute gs_coeff^T * matrix_ks * gs_coeff to get the epsilon_ij matrix
               DO ispin = 1, 2
                  CALL cp_dbcsr_sm_fm_multiply(matrix_ks(ispin)%matrix, gs_coeffs, work_mat, ncol=ndo_so)
                  CALL cp_gemm('T', 'N', ndo_so, ndo_so, nao, 1.0_dp, gs_coeffs, work_mat, 0.0_dp, eval_mat)

                  CALL cp_fm_get_diag(eval_mat, diag)
                  donor_state%energy_evals(:, ispin) = diag(:)

               END DO

               DEALLOCATE (diag)
            END IF

!           Clean-up
            CALL cp_fm_release(work_mat)
            CALL cp_fm_release(eval_mat)
            CALL cp_fm_struct_release(eval_mat_struct)

         END IF ! do_localize and/or ROKS

!        Clean-up
         CALL cp_fm_struct_release(gs_struct)
         DEALLOCATE (first_sgf)

         IF (output_unit > 0) WRITE (UNIT=output_unit, FMT="(T5,A)") " "

         DO ispin = 1, nspins
            IF (output_unit > 0) THEN
               WRITE (UNIT=output_unit, FMT="(T5,A,I1,A,F7.5,A,I4)") &
                  "The next best overlap for spin ", ispin, " is ", next_best_overlap(ispin), &
                  " for MO with index ", next_best_overlap_ind(ispin)
            END IF
         END DO
         IF (output_unit > 0) WRITE (UNIT=output_unit, FMT="(T5,A)") " "

      ELSE
         CPABORT("A core donor state could not be assigned MO(s). Increasing NSEARCH might help.")
      END IF

   END SUBROUTINE assign_mos_to_donor_state

! **************************************************************************************************
!> \brief Compute the centers and spreads of (core) MOs using the Berry phase operator
!> \param xas_tdp_env ...
!> \param xas_tdp_control ...
!> \param qs_env ...
!> \note xas_tdp_env%qs_loc_env is used and modified. OK since no localization done after this
!>       subroutine is used.
! **************************************************************************************************
   SUBROUTINE find_mo_centers(xas_tdp_env, xas_tdp_control, qs_env)

      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'find_mo_centers', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: dim_op, i, ispin, j, n_centers, nao, &
                                                            nspins
      REAL(dp), DIMENSION(6)                             :: weights
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: moloc_coeff
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: zij_fm_set
      TYPE(cp_fm_struct_type), POINTER                   :: tmp_fm_struct
      TYPE(cp_fm_type), POINTER                          :: opvec, vectors
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: op_sm_set
      TYPE(qs_loc_env_new_type), POINTER                 :: qs_loc_env
      TYPE(section_vals_type), POINTER                   :: print_loc_section, prog_run_info

      NULLIFY (qs_loc_env, zij_fm_set, cell, print_loc_section, op_sm_set, moloc_coeff, vectors)
      NULLIFY (opvec, tmp_fm_struct, para_env, blacs_env, prog_run_info)

!  Initialization
      print_loc_section => xas_tdp_control%print_loc_subsection
      n_centers = xas_tdp_control%n_search
      CALL get_qs_env(qs_env=qs_env, para_env=para_env, blacs_env=blacs_env, cell=cell)

!  Set print option to debug to keep clean output file
      prog_run_info => section_vals_get_subs_vals(print_loc_section, "PROGRAM_RUN_INFO")
      CALL section_vals_val_set(prog_run_info, keyword_name="_SECTION_PARAMETERS_", &
                                i_val=debug_print_level)

!  Re-initialize the qs_loc_env to get the current MOs. Use force_loc because needed for centers
      CALL reinit_qs_loc_env(xas_tdp_env%qs_loc_env, n_centers, xas_tdp_control%do_uks, qs_env)
      qs_loc_env => xas_tdp_env%qs_loc_env

!  Get what we need from the qs_lovc_env
      CALL get_qs_loc_env(qs_loc_env=qs_loc_env, weights=weights, op_sm_set=op_sm_set, &
                          moloc_coeff=moloc_coeff)

!  Prepare for zij
      vectors => moloc_coeff(1)%matrix
      CALL cp_fm_get_info(vectors, nrow_global=nao)
      CALL cp_fm_create(opvec, vectors%matrix_struct)

      CALL cp_fm_struct_create(tmp_fm_struct, para_env=para_env, context=blacs_env, &
                               ncol_global=n_centers, nrow_global=n_centers)

      IF (cell%orthorhombic) THEN
         dim_op = 3
      ELSE
         dim_op = 6
      END IF
      ALLOCATE (zij_fm_set(2, dim_op))
      DO i = 1, dim_op
         DO j = 1, 2
            CALL cp_fm_create(zij_fm_set(j, i)%matrix, tmp_fm_struct)
         END DO
      END DO

      !  If spin-unrestricted, need to go spin by spin
      nspins = 1; IF (xas_tdp_control%do_uks) nspins = 2

      DO ispin = 1, nspins
!     zij computation, copied from qs_loc_methods:optimize_loc_berry
         vectors => moloc_coeff(ispin)%matrix
         DO i = 1, dim_op
            DO j = 1, 2
               CALL cp_fm_set_all(zij_fm_set(j, i)%matrix, 0.0_dp)
               CALL cp_dbcsr_sm_fm_multiply(op_sm_set(j, i)%matrix, vectors, opvec, ncol=n_centers)
               CALL cp_gemm("T", "N", n_centers, n_centers, nao, 1.0_dp, vectors, opvec, 0.0_dp, &
                            zij_fm_set(j, i)%matrix)
            END DO
         END DO

!     Compute centers (and spread)
         CALL centers_spreads_berry(qs_loc_env=qs_loc_env, zij=zij_fm_set, nmoloc=n_centers, &
                                    cell=cell, weights=weights, ispin=ispin, &
                                    print_loc_section=print_loc_section, only_initial_out=.TRUE.)
      END DO !ispins

!  Clean-up
      CALL cp_fm_release(opvec)
      CALL cp_fm_struct_release(tmp_fm_struct)
      DO i = 1, dim_op
         DO j = 1, 2
            CALL cp_fm_release(zij_fm_set(j, i)%matrix)
         END DO
      END DO
      DEALLOCATE (zij_fm_set)

!  Make sure we leave with the correct do_loc value
      qs_loc_env%do_localize = xas_tdp_control%do_loc

   END SUBROUTINE find_mo_centers

! **************************************************************************************************
!> \brief Prints the MO to donor_state assocaition with overlap and Mulliken population analysis
!> \param xas_tdp_env ...
!> \param xas_tdp_control ...
!> \param qs_env ...
!> \note  Called only in case of CHECK_ONLY run
! **************************************************************************************************
   SUBROUTINE print_checks(xas_tdp_env, xas_tdp_control, qs_env)

      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'print_checks', routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_string_length)               :: kind_name
      INTEGER                                            :: current_state_index, iat, iatom, ikind, &
                                                            istate, output_unit, tmp_index
      INTEGER, DIMENSION(:), POINTER                     :: atoms_of_kind
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(donor_state_type), POINTER                    :: current_state

      NULLIFY (atomic_kind_set, atoms_of_kind, current_state)

      output_unit = cp_logger_get_default_io_unit()

      IF (output_unit > 0) THEN
         WRITE (output_unit, "(/,T3,A,/,T3,A,/,T3,A)") &
            "# Check the donor states for their quality. They need to have a well defined type ", &
            "  (1s, 2s, etc) which is indicated by the overlap. They also need to be localized, ", &
            "  for which the Mulliken population analysis is one indicator (must be close to 1.0)"
      END IF

!  Loop over the donor states (as in the main xas_tdp loop)
      CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set)
      current_state_index = 1

      !loop over atomic kinds
      DO ikind = 1, SIZE(atomic_kind_set)

         CALL get_atomic_kind(atomic_kind=atomic_kind_set(ikind), name=kind_name, &
                              atom_list=atoms_of_kind)

         IF (.NOT. ANY(xas_tdp_env%ex_kind_indices == ikind)) CYCLE

         !loop over atoms of kind
         DO iat = 1, SIZE(atoms_of_kind)
            iatom = atoms_of_kind(iat)

            IF (.NOT. ANY(xas_tdp_env%ex_atom_indices == iatom)) CYCLE
            tmp_index = locate(xas_tdp_env%ex_atom_indices, iatom)

            !loop over states of excited atom
            DO istate = 1, SIZE(xas_tdp_env%state_types, 1)

               IF (xas_tdp_env%state_types(istate, tmp_index) == xas_not_excited) CYCLE

               current_state => xas_tdp_env%donor_states(current_state_index)
               CALL set_donor_state(current_state, at_index=iatom, &
                                    at_symbol=kind_name, kind_index=ikind, &
                                    state_type=xas_tdp_env%state_types(istate, tmp_index))

               IF (output_unit > 0) THEN
                  WRITE (output_unit, "(/,T4,A,A2,A,I4,A,A,A)") &
                     "-Donor state of type ", xas_tdp_env%state_type_char(current_state%state_type), &
                     " for atom", current_state%at_index, " of kind ", TRIM(current_state%at_symbol), ":"
               END IF

               !Assign the MOs and perform Mulliken
               CALL assign_mos_to_donor_state(current_state, xas_tdp_env, xas_tdp_control, qs_env)
               CALL perform_mulliken_on_donor_state(current_state, qs_env)

               current_state_index = current_state_index+1
               NULLIFY (current_state)

            END DO !istate
         END DO !iat
      END DO !ikind

      IF (output_unit > 0) THEN
         WRITE (output_unit, "(/,T5,A)") &
            "Use LOCALIZE and/or increase N_SEARCH for better results, if so required."
      END IF

   END SUBROUTINE print_checks

! **************************************************************************************************
!> \brief Computes the required multipole moment in the lenght representation for a given atom
!> \param iatom index of the given atom
!> \param xas_tdp_env ...
!> \param xas_tdp_control ...
!> \param qs_env ...
!> \note Assumes that wither dipole or quadrupole in length rep is required
! **************************************************************************************************
   SUBROUTINE compute_lenrep_multipole(iatom, xas_tdp_env, xas_tdp_control, qs_env)

      INTEGER, INTENT(IN)                                :: iatom
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'compute_lenrep_multipole', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, order
      REAL(dp), DIMENSION(3)                             :: rc
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: work
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      NULLIFY (work, particle_set)

      CALL get_qs_env(qs_env, particle_set=particle_set)
      rc = particle_set(iatom)%r

      ALLOCATE (work(9))
      IF (xas_tdp_control%dipole_form == xas_dip_len) THEN
         DO i = 1, 3
            CALL dbcsr_set(xas_tdp_env%dipmat(i)%matrix, 0.0_dp)
            work(i)%matrix => xas_tdp_env%dipmat(i)%matrix
         END DO
         order = 1
      END IF
      IF (xas_tdp_control%do_quad) THEN
         DO i = 1, 6
            CALL dbcsr_set(xas_tdp_env%quadmat(i)%matrix, 0.0_dp)
            work(3+i)%matrix => xas_tdp_env%quadmat(i)%matrix
         END DO
         order = 2
         IF (xas_tdp_control%dipole_form == xas_dip_vel) order = -2
      END IF

      !enforce minimum image to avoid PBCs related issues, ok because localized densities
      CALL rRc_xyz_ao(work, qs_env, rc, order=order, minimum_image=.TRUE.)
      DEALLOCATE (work)

   END SUBROUTINE compute_lenrep_multipole

! **************************************************************************************************
!> \brief Computes the oscillator strength based on the dipole moment (velocity or length rep) for
!>        all available excitation energies and store the results in the donor_state. There is no
!>        triplet dipole in the spin-restricted ground state.
!> \param donor_state the donor state which is excited
!> \param xas_tdp_control ...
!> \param xas_tdp_env ...
!> \note The oscillator strength is a scalar: osc_str = 2/(3*omega)*(dipole_v)^2 in the velocity rep
!>       or : osc_str = 2/3*omega*(dipole_r)^2 in the length representation
!>       The formulae for the dipoles come from the trace of the dipole operator with the transition
!>       densities, i.e. what we get from solving the xas_tdp problem. Same procedure with or wo TDA
! **************************************************************************************************
   SUBROUTINE compute_dipole_fosc(donor_state, xas_tdp_control, xas_tdp_env)

      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env

      CHARACTER(len=*), PARAMETER :: routineN = 'compute_dipole_fosc', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iosc, j, nao, ndo_mo, ndo_so, &
                                                            ngs, nosc, nspins
      LOGICAL                                            :: do_sc, do_sg
      REAL(dp)                                           :: pref
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: tot_contr
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: dip_block
      REAL(dp), DIMENSION(:), POINTER                    :: lr_evals, osc_str
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: col_struct, mat_struct
      TYPE(cp_fm_type), POINTER                          :: col_work, lr_coeffs, mat_work
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: dipmat

      NULLIFY (dipmat, col_work, mat_work, col_struct, mat_struct, para_env, blacs_env, lr_coeffs)
      NULLIFY (lr_evals, osc_str)

      CALL timeset(routineN, handle)

!  Initialization
      do_sc = xas_tdp_control%do_spin_cons
      do_sg = xas_tdp_control%do_singlet
      IF (do_sc) THEN
         nspins = 2
         lr_evals => donor_state%sc_evals
         lr_coeffs => donor_state%sc_coeffs
      ELSE IF (do_sg) THEN
         nspins = 1
         lr_evals => donor_state%sg_evals
         lr_coeffs => donor_state%sg_coeffs
      ELSE
         CPABORT("Dipole oscilaltor strength only for singlets and spin-conserving excitations.")
      END IF
      ndo_mo = donor_state%ndo_mo
      ndo_so = ndo_mo*nspins
      ngs = ndo_so; IF (xas_tdp_control%do_roks) ngs = ndo_mo !in ROKS, same gs coeffs
      nosc = SIZE(lr_evals)
      ALLOCATE (donor_state%osc_str(nosc))
      osc_str => donor_state%osc_str
      osc_str = 0.0_dp
      dipmat => xas_tdp_env%dipmat

      ! do some work matrix initialization
      CALL cp_fm_get_info(donor_state%gs_coeffs, matrix_struct=col_struct, para_env=para_env, &
                          context=blacs_env, nrow_global=nao)
      CALL cp_fm_struct_create(mat_struct, para_env=para_env, context=blacs_env, &
                               nrow_global=ndo_so*nosc, ncol_global=ngs)
      CALL cp_fm_create(mat_work, mat_struct)
      CALL cp_fm_create(col_work, col_struct)

      ALLOCATE (tot_contr(ndo_mo), dip_block(ndo_so, ngs))
      pref = 2.0_dp; IF (do_sc) pref = 1.0_dp !because of singlet definition u = 1/sqrt(2)(c_a+c_b)

!  Looping over cartesian coord
      DO j = 1, 3

         !Compute dip*gs_coeffs
         CALL cp_dbcsr_sm_fm_multiply(dipmat(j)%matrix, donor_state%gs_coeffs, col_work, ncol=ngs)
         !compute lr_coeffs*dip*gs_coeffs
         CALL cp_gemm('T', 'N', ndo_so*nosc, ngs, nao, 1.0_dp, lr_coeffs, col_work, 0.0_dp, mat_work)

         !Loop over the excited states
         DO iosc = 1, nosc

            tot_contr = 0.0_dp
            CALL cp_fm_get_submatrix(fm=mat_work, target_m=dip_block, start_row=(iosc-1)*ndo_so+1, &
                                     start_col=1, n_rows=ndo_so, n_cols=ngs)
            IF (do_sg) THEN
               tot_contr(:) = get_diag(dip_block)
            ELSE IF (do_sc .AND. xas_tdp_control%do_uks) THEN
               tot_contr(:) = get_diag(dip_block(1:ndo_mo, 1:ndo_mo)) !alpha
               tot_contr(:) = tot_contr(:)+get_diag(dip_block(ndo_mo+1:ndo_so, ndo_mo+1:ndo_so)) !beta
            ELSE
               !roks
               tot_contr(:) = get_diag(dip_block(1:ndo_mo, :)) !alpha
               tot_contr(:) = tot_contr(:)+get_diag(dip_block(ndo_mo+1:ndo_so, :)) !beta
            END IF

            osc_str(iosc) = osc_str(iosc)+SUM(tot_contr)**2

         END DO !iosc
      END DO !j

      !compute the prefactor
      IF (xas_tdp_control%dipole_form == xas_dip_len) THEN
         osc_str(:) = pref*2.0_dp/3.0_dp*lr_evals(:)*osc_str(:)
      ELSE
         osc_str(:) = pref*2.0_dp/3.0_dp/lr_evals(:)*osc_str(:)
      END IF

      !clean-up
      CALL cp_fm_release(mat_work)
      CALL cp_fm_release(col_work)
      CALL cp_fm_struct_release(mat_struct)

      CALL timestop(handle)

   END SUBROUTINE compute_dipole_fosc

! **************************************************************************************************
!> \brief Computes the oscillator strength due to the electric quadrupole moment and store it in
!>        the donor_state (for singlet or spin-conserving)
!> \param donor_state the donor state which is excited
!> \param xas_tdp_control ...
!> \param xas_tdp_env ...
!> \note Formula: 1/20*a_fine^2*omega^3 * sum_ab (sum_i r_ia*r_ib - 1/3*ri^2*delta_ab)
! **************************************************************************************************
   SUBROUTINE compute_quadrupole_fosc(donor_state, xas_tdp_control, xas_tdp_env)

      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env

      CHARACTER(len=*), PARAMETER :: routineN = 'compute_quadrupole_fosc', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iosc, j, nao, ndo_mo, ndo_so, &
                                                            ngs, nosc, nspins
      LOGICAL                                            :: do_sc, do_sg
      REAL(dp)                                           :: pref
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: tot_contr, trace
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: quad_block
      REAL(dp), DIMENSION(:), POINTER                    :: lr_evals, osc_str
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: col_struct, mat_struct
      TYPE(cp_fm_type), POINTER                          :: col_work, lr_coeffs, mat_work
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: quadmat

      NULLIFY (lr_evals, osc_str, mat_work, lr_coeffs, col_work, col_struct, mat_struct, para_env)
      NULLIFY (blacs_env)

      CALL timeset(routineN, handle)

      ! Initialization
      do_sc = xas_tdp_control%do_spin_cons
      do_sg = xas_tdp_control%do_singlet
      IF (do_sc) THEN
         nspins = 2
         lr_evals => donor_state%sc_evals
         lr_coeffs => donor_state%sc_coeffs
      ELSE IF (do_sg) THEN
         nspins = 1
         lr_evals => donor_state%sg_evals
         lr_coeffs => donor_state%sg_coeffs
      ELSE
         CPABORT("Quadrupole oscillator strengths only for singlet and spin-conserving excitations")
      END IF
      ndo_mo = donor_state%ndo_mo
      ndo_so = ndo_mo*nspins
      ngs = ndo_so; IF (xas_tdp_control%do_roks) ngs = ndo_mo !only alpha do_mo in ROKS
      nosc = SIZE(lr_evals)
      ALLOCATE (donor_state%quad_osc_str(nosc))
      osc_str => donor_state%quad_osc_str
      osc_str = 0.0_dp
      quadmat => xas_tdp_env%quadmat

      !work matrices init
      CALL cp_fm_get_info(donor_state%gs_coeffs, matrix_struct=col_struct, para_env=para_env, &
                          context=blacs_env, nrow_global=nao)
      CALL cp_fm_struct_create(mat_struct, para_env=para_env, context=blacs_env, &
                               nrow_global=ndo_so*nosc, ncol_global=ngs)
      CALL cp_fm_create(mat_work, mat_struct)
      CALL cp_fm_create(col_work, col_struct)

      ALLOCATE (quad_block(ndo_so, ngs), tot_contr(ndo_mo))
      pref = 2.0_dp; IF (do_sc) pref = 1.0_dp !because of singlet definition u = 1/sqrt(2)*...
      ALLOCATE (trace(nosc))
      trace = 0.0_dp

      !Loop over the cartesioan coord :x2, xy, xz, y2, yz, z2
      DO j = 1, 6

         !Compute quad*gs_coeffs
         CALL cp_dbcsr_sm_fm_multiply(quadmat(j)%matrix, donor_state%gs_coeffs, col_work, ncol=ngs)
         !compute lr_coeffs*quadmat*gs_coeffs
         CALL cp_gemm('T', 'N', ndo_so*nosc, ngs, nao, 1.0_dp, lr_coeffs, col_work, 0.0_dp, mat_work)

         !Loop over the excited states
         DO iosc = 1, nosc

            tot_contr = 0.0_dp
            CALL cp_fm_get_submatrix(fm=mat_work, target_m=quad_block, start_row=(iosc-1)*ndo_so+1, &
                                     start_col=1, n_rows=ndo_so, n_cols=ngs)

            IF (do_sg) THEN
               tot_contr(:) = get_diag(quad_block)
            ELSE IF (do_sc .AND. xas_tdp_control%do_uks) THEN
               tot_contr(:) = get_diag(quad_block(1:ndo_mo, 1:ndo_mo)) !alpha
               tot_contr(:) = tot_contr(:)+get_diag(quad_block(ndo_mo+1:ndo_so, ndo_mo+1:ndo_so)) !beta
            ELSE
               !roks
               tot_contr(:) = get_diag(quad_block(1:ndo_mo, :)) !alpha
               tot_contr(:) = tot_contr(:)+get_diag(quad_block(ndo_mo+1:ndo_so, :)) !beta
            END IF

            !if x2, y2, or z2 direction, need to update the trace (for later)
            IF (j == 1 .OR. j == 4 .OR. j == 6) THEN
               osc_str(iosc) = osc_str(iosc)+SUM(tot_contr)**2
               trace(iosc) = trace(iosc)+SUM(tot_contr)

               !if xy, xz or yz, need to count twice the contribution (for yx, zx and zy)
            ELSE
               osc_str(iosc) = osc_str(iosc)+2.0_dp*SUM(tot_contr)**2
            END IF

         END DO !iosc
      END DO !j

      !compute the prefactor, and remove 1/3*trace^2
      osc_str(:) = pref*1._dp/20._dp*a_fine**2*lr_evals(:)**3*(osc_str(:)-1._dp/3._dp*trace(:)**2)

      !clean-up
      CALL cp_fm_release(mat_work)
      CALL cp_fm_release(col_work)
      CALL cp_fm_struct_release(mat_struct)

      CALL timestop(handle)

   END SUBROUTINE compute_quadrupole_fosc

! **************************************************************************************************
!> \brief Writes the core MOs to excited atoms associations in the main output file
!> \param xas_tdp_env ...
!> \param xas_tdp_control ...
!> \param qs_env ...
!> \note Look at alpha spin MOs, as we are dealing with core states and alpha/beta MOs are the same
! **************************************************************************************************
   SUBROUTINE write_mos_to_ex_atoms_association(xas_tdp_env, xas_tdp_control, qs_env)

      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'write_mos_to_ex_atoms_association', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_string_length)               :: kind_name
      INTEGER                                            :: at_index, imo, ispin, nmo, nspins, &
                                                            output_unit, tmp_index
      INTEGER, DIMENSION(3)                              :: perd_init
      INTEGER, DIMENSION(:), POINTER                     :: ex_atom_indices
      INTEGER, DIMENSION(:, :, :), POINTER               :: mos_of_ex_atoms
      REAL(dp)                                           :: dist, mo_spread
      REAL(dp), DIMENSION(3)                             :: at_pos, r_ac, wfn_center
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      NULLIFY (cell, particle_set, mos_of_ex_atoms, ex_atom_indices)

      output_unit = cp_logger_get_default_io_unit()

      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T3,A,/,T3,A,/,T3,A)") &
            "                  Associated    Associated        Distance to   MO spread (Ang^2)", &
            "Spin  MO index    atom index     atom kind    MO center (Ang)   -w_i ln(|z_ij|^2)", &
            "---------------------------------------------------------------------------------"
      END IF

!  Initialization
      nspins = 1; IF (xas_tdp_control%do_uks) nspins = 2
      mos_of_ex_atoms => xas_tdp_env%mos_of_ex_atoms
      ex_atom_indices => xas_tdp_env%ex_atom_indices
      nmo = xas_tdp_control%n_search
      CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, cell=cell)

!     because the use of Berry's phase operator implies PBCs
      perd_init = cell%perd
      cell%perd = 1

!  Retrieving all the info for each MO and spin
      DO imo = 1, nmo
         DO ispin = 1, nspins

!           each Mo is associated to at most one atom (only 1 in array of -1)
            IF (ANY(mos_of_ex_atoms(imo, :, ispin) == 1)) THEN
               tmp_index = MAXLOC(mos_of_ex_atoms(imo, :, ispin), 1)
               at_index = ex_atom_indices(tmp_index)
               kind_name = particle_set(at_index)%atomic_kind%name

               at_pos = particle_set(at_index)%r
               wfn_center = xas_tdp_env%qs_loc_env%localized_wfn_control%centers_set(ispin)%array(1:3, imo)
               r_ac = pbc(at_pos, wfn_center, cell)
               dist = NORM2(r_ac)
!              convert distance from a.u. to Angstrom
               dist = dist*angstrom

               mo_spread = xas_tdp_env%qs_loc_env%localized_wfn_control%centers_set(ispin)%array(4, imo)
               mo_spread = mo_spread*angstrom*angstrom

               IF (output_unit > 0) THEN
                  WRITE (UNIT=output_unit, FMT="(T3,I4,I10,I14,A14,ES19.3,ES20.3)") &
                     ispin, imo, at_index, TRIM(kind_name), dist, mo_spread
               END IF

            END IF
         END DO !ispin
      END DO !imo

      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(T3,A,/)") &
            "---------------------------------------------------------------------------------"
      END IF

!  Go back to initial BCs
      cell%perd = perd_init

   END SUBROUTINE write_mos_to_ex_atoms_association

! **************************************************************************************************
!> \brief Performs Mulliken population analysis for the MO(s) of a donor_state_type so that user
!>        can verify it is indeed a core state
!> \param donor_state ...
!> \param qs_env ...
!> \note This is a specific case of Mulliken analysis. In general one computes sum_i (SP)_ii, where
!>       i labels the basis function centered on the atom of interest. For a specific MO with index
!>       j, one need to compute sum_{ik} c_{ij} S_{ik} c_{kj}, k = 1,nao
! **************************************************************************************************
   SUBROUTINE perform_mulliken_on_donor_state(donor_state, qs_env)
      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'perform_mulliken_on_donor_state', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: at_index, i, ispin, nao, natom, ndo_mo, &
                                                            ndo_so, nsgf, nspins, output_unit
      INTEGER, DIMENSION(:), POINTER                     :: first_sgf, last_sgf
      INTEGER, DIMENSION(:, :), POINTER                  :: mo_indices
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: mul_pop, pop_mat
      REAL(dp), DIMENSION(:, :), POINTER                 :: work_array
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: col_vect_struct
      TYPE(cp_fm_type), POINTER                          :: gs_coeffs, work_vect
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (mo_indices, qs_kind_set, particle_set, first_sgf, work_array)
      NULLIFY (matrix_s, para_env, blacs_env, col_vect_struct, work_vect, last_sgf)

!  Initialization
      at_index = donor_state%at_index
      mo_indices => donor_state%mo_indices
      ndo_mo = donor_state%ndo_mo
      gs_coeffs => donor_state%gs_coeffs
      output_unit = cp_logger_get_default_io_unit()
      nspins = 1; IF (SIZE(mo_indices, 2) == 2) nspins = 2
      ndo_so = ndo_mo*nspins
      ALLOCATE (mul_pop(ndo_mo, nspins))
      mul_pop = 0.0_dp

      CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, qs_kind_set=qs_kind_set, &
                      para_env=para_env, blacs_env=blacs_env, matrix_s=matrix_s)
      CALL cp_fm_get_info(gs_coeffs, nrow_global=nao, matrix_struct=col_vect_struct)

      natom = SIZE(particle_set, 1)
      ALLOCATE (first_sgf(natom))
      ALLOCATE (last_sgf(natom))

      CALL get_particle_set(particle_set, qs_kind_set, first_sgf=first_sgf, last_sgf=last_sgf)
      nsgf = last_sgf(at_index)-first_sgf(at_index)+1

      CALL cp_fm_create(work_vect, col_vect_struct)

!  Take the product of S*coeffs
      CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, gs_coeffs, work_vect, ncol=ndo_so)

!  Only consider the product coeffs^T * S * coeffs on the atom of interest
      ALLOCATE (work_array(nsgf, ndo_so))
      ALLOCATE (pop_mat(ndo_so, ndo_so))

      CALL cp_fm_get_submatrix(fm=work_vect, target_m=work_array, start_row=first_sgf(at_index), &
                               start_col=1, n_rows=nsgf, n_cols=ndo_so, transpose=.FALSE.)

      CALL dgemm('T', 'N', ndo_so, ndo_so, nsgf, 1.0_dp, donor_state%contract_coeffs, nsgf, &
                 work_array, nsgf, 0.0_dp, pop_mat, ndo_so)

!  The Mullikan population for the MOs in on the diagonal.
      DO ispin = 1, nspins
         DO i = 1, ndo_mo
            mul_pop(i, ispin) = pop_mat((ispin-1)*ndo_mo+i, (ispin-1)*ndo_mo+i)
         END DO
      END DO

!  Printing in main output file
      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(T5,A,/,T5,A)") &
            "Mulliken population analysis retricted to the associated MO(s) yields: ", &
            "                                              Spin  MO index     charge"
         DO ispin = 1, nspins
            DO i = 1, ndo_mo
               WRITE (UNIT=output_unit, FMT="(T51,I4,I10,F11.3)") &
                  ispin, mo_indices(i, ispin), mul_pop(i, ispin)
            END DO
         END DO
      END IF

!  Clean-up
      DEALLOCATE (first_sgf, last_sgf, work_array)
      CALL cp_fm_release(work_vect)

   END SUBROUTINE perform_mulliken_on_donor_state

! **************************************************************************************************
!> \brief write the PDOS wrt the LR-orbitals for the current donor_state and/or the CUBES files
!> \param ex_type the excitation type: singlet, triplet, spin-conserving, etc
!> \param donor_state ...
!> \param xas_tdp_env ...
!> \param xas_tdp_section ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE xas_tdp_post(ex_type, donor_state, xas_tdp_env, xas_tdp_section, qs_env)

      INTEGER, INTENT(IN)                                :: ex_type
      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(section_vals_type), POINTER                   :: xas_tdp_section
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'xas_tdp_post', routineP = moduleN//':'//routineN

      CHARACTER(len=default_string_length)               :: domo, domon, excite, pos, xas_mittle
      INTEGER                                            :: handle, ic, ido_mo, imo, irep, ispin, &
                                                            n_dependent, n_rep, nao, ncubes, &
                                                            ndo_mo, ndo_so, nlumo, nmo, nspins, &
                                                            output_unit
      INTEGER, DIMENSION(:), POINTER                     :: bounds, list, state_list
      LOGICAL                                            :: append_cube, do_cubes, do_pdos
      REAL(dp), DIMENSION(:), POINTER                    :: lr_evals
      REAL(dp), DIMENSION(:, :), POINTER                 :: centers
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct, mo_struct
      TYPE(cp_fm_type), POINTER                          :: lr_coeffs, mo_coeff, work_fm
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(mo_set_type), POINTER                         :: mo_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: print_key

      NULLIFY (atomic_kind_set, particle_set, qs_kind_set, mo_set, lr_evals, lr_coeffs, mo_coeff)
      NULLIFY (mo_struct, para_env, blacs_env, fm_struct, matrix_s, work_fm, print_key, logger)
      NULLIFY (bounds, state_list, list)

      !Tests on what to do
      logger => cp_get_default_logger()
      do_pdos = .FALSE.; do_cubes = .FALSE.

      IF (BTEST(cp_print_key_should_output(logger%iter_info, xas_tdp_section, &
                                           "PRINT%PDOS"), cp_p_file)) do_pdos = .TRUE.

      IF (BTEST(cp_print_key_should_output(logger%iter_info, xas_tdp_section, &
                                           "PRINT%CUBES"), cp_p_file)) do_cubes = .TRUE.

      IF (.NOT. (do_pdos .OR. do_cubes)) RETURN

      CALL timeset(routineN, handle)

      !Initialization
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, particle_set=particle_set, &
                      qs_kind_set=qs_kind_set, para_env=para_env, blacs_env=blacs_env, &
                      matrix_s=matrix_s)

      SELECT CASE (ex_type)
      CASE (tddfpt_spin_cons)
         lr_evals => donor_state%sc_evals
         lr_coeffs => donor_state%sc_coeffs
         nspins = 2
         excite = "spincons"
      CASE (tddfpt_spin_flip)
         lr_evals => donor_state%sf_evals
         lr_coeffs => donor_state%sf_coeffs
         nspins = 2
         excite = "spinflip"
      CASE (tddfpt_singlet)
         lr_evals => donor_state%sg_evals
         lr_coeffs => donor_state%sg_coeffs
         nspins = 1
         excite = "singlet"
      CASE (tddfpt_triplet)
         lr_evals => donor_state%tp_evals
         lr_coeffs => donor_state%tp_coeffs
         nspins = 1
         excite = "triplet"
      END SELECT

      SELECT CASE (donor_state%state_type)
      CASE (xas_1s_type)
         domo = "1s"
      CASE (xas_2s_type)
         domo = "2s"
      CASE (xas_2p_type)
         domo = "2p"
      END SELECT

      ndo_mo = donor_state%ndo_mo
      ndo_so = ndo_mo*nspins
      nmo = SIZE(lr_evals)
      CALL dbcsr_get_info(matrix_s(1)%matrix, nfullrows_total=nao)

      CALL cp_fm_struct_create(mo_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=nao, ncol_global=nmo)
      CALL cp_fm_create(mo_coeff, mo_struct)

      !PDOS related stuff
      IF (do_pdos) THEN

         !If S^0.5 not yet stored, compute it once and for all
         IF (.NOT. ASSOCIATED(xas_tdp_env%matrix_shalf) .AND. do_pdos) THEN
            CALL cp_fm_struct_create(fm_struct, context=blacs_env, para_env=para_env, &
                                     nrow_global=nao, ncol_global=nao)
            CALL cp_fm_create(xas_tdp_env%matrix_shalf, fm_struct)
            CALL cp_fm_create(work_fm, fm_struct)

            CALL copy_dbcsr_to_fm(matrix_s(1)%matrix, xas_tdp_env%matrix_shalf)
            CALL cp_fm_power(xas_tdp_env%matrix_shalf, work_fm, 0.5_dp, EPSILON(0.0_dp), n_dependent)

            CALL cp_fm_release(work_fm)
            CALL cp_fm_struct_release(fm_struct)
         END IF

         !Giving some PDOS info
         output_unit = cp_logger_get_default_io_unit()
         IF (output_unit > 0) THEN
            WRITE (UNIT=output_unit, FMT="(/,T5,A,/,T5,A,/,T5,A)") &
               "Computing the PDOS of linear-response orbitals for spectral features analysis", &
               "Note: using standard PDOS routines => ignore mentions of KS states and MO ", &
               "      occupation numbers. Eigenvalues in *.pdos files are excitations energies."
         END IF

         !Check on NLUMO
         CALL section_vals_val_get(xas_tdp_section, "PRINT%PDOS%NLUMO", i_val=nlumo)
         IF (nlumo .NE. 0) THEN
            CPWARN("NLUMO is irrelevant for XAS_TDP PDOS. It was overwritten to 0.")
         END IF
         CALL section_vals_val_set(xas_tdp_section, "PRINT%PDOS%NLUMO", i_val=0)
      END IF

      !CUBES related stuff
      IF (do_cubes) THEN

         print_key => section_vals_get_subs_vals(xas_tdp_section, "PRINT%CUBES")

         CALL section_vals_val_get(print_key, "CUBES_LU_BOUNDS", i_vals=bounds)
         ncubes = bounds(2)-bounds(1)+1
         IF (ncubes > 0) THEN
            ALLOCATE (state_list(ncubes))
            DO ic = 1, ncubes
               state_list(ic) = bounds(1)+ic-1
            END DO
         END IF

         IF (.NOT. ASSOCIATED(state_list)) THEN
            CALL section_vals_val_get(print_key, "CUBES_LIST", n_rep_val=n_rep)

            ncubes = 0
            DO irep = 1, n_rep
               NULLIFY (list)
               CALL section_vals_val_get(print_key, "CUBES_LIST", i_rep_val=irep, i_vals=list)
               IF (ASSOCIATED(list)) THEN
                  CALL reallocate(state_list, 1, ncubes+SIZE(list))
                  DO ic = 1, SIZE(list)
                     state_list(ncubes+ic) = list(ic)
                  END DO
                  ncubes = ncubes+SIZE(list)
               END IF
            END DO
         END IF

         IF (.NOT. ASSOCIATED(state_list)) THEN
            ncubes = 1
            ALLOCATE (state_list(1))
            state_list(1) = 1
         END IF

         CALL section_vals_val_get(print_key, "APPEND", l_val=append_cube)
         pos = "REWIND"
         IF (append_cube) pos = "APPEND"

         ALLOCATE (centers(6, ncubes))
         centers = 0.0_dp

      END IF

      !Loop over MOs and spin, one PDOS/CUBE for each
      DO ido_mo = 1, ndo_mo
         DO ispin = 1, nspins

            !need to create a mo set for the LR-orbitals
            CALL allocate_mo_set(mo_set, nao=nao, nmo=nmo, nelectron=nmo, n_el_f=REAL(nmo, dp), &
                                 maxocc=1.0_dp, flexible_electron_count=0.0_dp)
            CALL init_mo_set(mo_set, fm_ref=mo_coeff, name="PDOS XAS_TDP MOs")
            mo_set%eigenvalues(:) = lr_evals(:)

            !get the actual coeff => most common case: closed-shell K-edge, can directly take lr_coeffs
            IF (nspins == 1 .AND. ndo_mo == 1) THEN
               CALL cp_fm_to_fm(lr_coeffs, mo_set%mo_coeff)
            ELSE
               DO imo = 1, nmo
                  CALL cp_fm_to_fm_submat(msource=lr_coeffs, mtarget=mo_set%mo_coeff, &
                                          nrow=nao, ncol=1, s_firstrow=1, &
                                          s_firstcol=(imo-1)*ndo_so+(ispin-1)*ndo_mo+ido_mo, &
                                          t_firstrow=1, t_firstcol=imo)
               END DO
            END IF

            !naming the output
            domon = domo
            IF (donor_state%state_type == xas_2p_type) domon = TRIM(domo)//TRIM(ADJUSTL(cp_to_string(ido_mo)))
            xas_mittle = 'xasat'//TRIM(ADJUSTL(cp_to_string(donor_state%at_index)))//'_'// &
                         TRIM(domon)//'_'//TRIM(excite)

            IF (do_pdos) THEN
               CALL calculate_projected_dos(mo_set, atomic_kind_set, qs_kind_set, particle_set, &
                                            qs_env, xas_tdp_section, ispin, xas_mittle, &
                                            external_matrix_shalf=xas_tdp_env%matrix_shalf)
            END IF

            IF (do_cubes) THEN
               CALL qs_print_cubes(qs_env, mo_set%mo_coeff, ncubes, state_list, centers, &
                                   print_key=print_key, root=xas_mittle, ispin=ispin, &
                                   file_position=pos)
            END IF

            !clean-up
            CALL deallocate_mo_set(mo_set)

         END DO
      END DO

      !clean-up
      CALL cp_fm_release(mo_coeff)
      CALL cp_fm_struct_release(mo_struct)
      IF (do_cubes) DEALLOCATE (centers, state_list)

      CALL timestop(handle)

   END SUBROUTINE xas_tdp_post

! **************************************************************************************************
!> \brief Prints the excitation energies and the oscillator strengths for a given donor_state in a file
!> \param donor_state the donor_state to print
!> \param xas_tdp_env ...
!> \param xas_tdp_control ...
!> \param xas_tdp_section ...
! **************************************************************************************************
   SUBROUTINE print_xas_tdp_to_file(donor_state, xas_tdp_env, xas_tdp_control, xas_tdp_section)

      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(section_vals_type), POINTER                   :: xas_tdp_section

      CHARACTER(len=*), PARAMETER :: routineN = 'print_xas_tdp_to_file', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, output_unit, xas_tdp_unit
      TYPE(cp_logger_type), POINTER                      :: logger

      NULLIFY (logger)
      logger => cp_get_default_logger()

      xas_tdp_unit = cp_print_key_unit_nr(logger, xas_tdp_section, "PRINT%SPECTRUM", &
                                          extension=".spectrum", file_position="APPEND", &
                                          file_action="WRITE", file_form="FORMATTED")

      output_unit = cp_logger_get_default_io_unit()

      IF (output_unit > 0) THEN
         WRITE (output_unit, FMT="(/,T5,A)") &
            "Calculations done: "
      END IF

      IF (xas_tdp_control%do_spin_cons) THEN
         IF (xas_tdp_unit > 0) THEN

!           Printing the general donor state information
            WRITE (xas_tdp_unit, FMT="(A,/,A,A,A/,A,I5,A,I5,A,A,/,A)") &
               "==================================================================================", &
               "XAS TDP open-shell spin-conserving (no SOC) excitations for DONOR STATE: ", &
               xas_tdp_env%state_type_char(donor_state%state_type), ",", &
               "from EXCITED ATOM: ", donor_state%at_index, ", of KIND (index/symbol): ", &
               donor_state%kind_index, "/", TRIM(donor_state%at_symbol), &
               "=================================================================================="

!           Simply dump the excitation energies/ oscillator strength as they come

            IF (xas_tdp_control%do_quad) THEN
               WRITE (xas_tdp_unit, FMT="(T3,A)") &
                  " Index     Excitation energy (eV)    fosc dipole (a.u.)   fosc quadrupole (a.u.)"
               DO i = 1, SIZE(donor_state%sc_evals)
                  WRITE (xas_tdp_unit, FMT="(T3,I6,F27.6,F22.6,F25.6)") &
                     i, donor_state%sc_evals(i)*evolt, donor_state%osc_str(i), &
                     donor_state%quad_osc_str(i)
               END DO
            ELSE
               WRITE (xas_tdp_unit, FMT="(T3,A)") &
                  " Index     Excitation energy (eV)    fosc dipole (a.u.)"
               DO i = 1, SIZE(donor_state%sc_evals)
                  WRITE (xas_tdp_unit, FMT="(T3,I6,F27.6,F22.6)") &
                     i, donor_state%sc_evals(i)*evolt, donor_state%osc_str(i)
               END DO
            END IF

            WRITE (xas_tdp_unit, FMT="(A,/)") " "
         END IF !xas_tdp_unit > 0

         IF (output_unit > 0) THEN
            WRITE (output_unit, FMT="(T5,A,F25.6)") &
               "First spin-conserving XAS excitation energy (eV): ", donor_state%sc_evals(1)*evolt
         END IF

      END IF ! do_spin_cons

      IF (xas_tdp_control%do_spin_flip) THEN
         IF (xas_tdp_unit > 0) THEN

!           Printing the general donor state information
            WRITE (xas_tdp_unit, FMT="(A,/,A,A,A/,A,I5,A,I5,A,A,/,A)") &
               "==================================================================================", &
               "XAS TDP open-shell spin-flip (no SOC) excitations for DONOR STATE: ", &
               xas_tdp_env%state_type_char(donor_state%state_type), ",", &
               "from EXCITED ATOM: ", donor_state%at_index, ", of KIND (index/symbol): ", &
               donor_state%kind_index, "/", TRIM(donor_state%at_symbol), &
               "==========================================================--======================"

!           Simply dump the excitation energies/ oscillator strength as they come

            IF (xas_tdp_control%do_quad) THEN
               WRITE (xas_tdp_unit, FMT="(T3,A)") &
                  " Index     Excitation energy (eV)    fosc dipole (a.u.)   fosc quadrupole (a.u.)"
               DO i = 1, SIZE(donor_state%sf_evals)
                  WRITE (xas_tdp_unit, FMT="(T3,I6,F27.6,F22.6,F25.6)") &
                     i, donor_state%sf_evals(i)*evolt, 0.0_dp, 0.0_dp !spin-forbidden !
               END DO
            ELSE
               WRITE (xas_tdp_unit, FMT="(T3,A)") &
                  " Index     Excitation energy (eV)    fosc dipole (a.u.)"
               DO i = 1, SIZE(donor_state%sf_evals)
                  WRITE (xas_tdp_unit, FMT="(T3,I6,F27.6,F22.6)") &
                     i, donor_state%sf_evals(i)*evolt, 0.0_dp
               END DO
            END IF

            WRITE (xas_tdp_unit, FMT="(A,/)") " "
         END IF !xas_tdp_unit

         IF (output_unit > 0) THEN
            WRITE (output_unit, FMT="(T5,A,F25.6)") &
               "First spin-flip XAS excitation energy (eV): ", donor_state%sf_evals(1)*evolt
         END IF
      END IF ! do_spin_flip

      IF (xas_tdp_control%do_singlet) THEN
         IF (xas_tdp_unit > 0) THEN

!           Printing the general donor state information
            WRITE (xas_tdp_unit, FMT="(A,/,A,A,A/,A,I5,A,I5,A,A,/,A)") &
               "==================================================================================", &
               "XAS TDP singlet excitations (no SOC) for DONOR STATE: ", &
               xas_tdp_env%state_type_char(donor_state%state_type), ",", &
               "from EXCITED ATOM: ", donor_state%at_index, ", of KIND (index/symbol): ", &
               donor_state%kind_index, "/", TRIM(donor_state%at_symbol), &
               "=================================================================================="

!           Simply dump the excitation energies/ oscillator strength as they come

            IF (xas_tdp_control%do_quad) THEN
               WRITE (xas_tdp_unit, FMT="(T3,A)") &
                  " Index     Excitation energy (eV)    fosc dipole (a.u.)   fosc quadrupole (a.u.)"
               DO i = 1, SIZE(donor_state%sg_evals)
                  WRITE (xas_tdp_unit, FMT="(T3,I6,F27.6,F22.6,F25.6)") &
                     i, donor_state%sg_evals(i)*evolt, donor_state%osc_str(i), &
                     donor_state%quad_osc_str(i)
               END DO
            ELSE
               WRITE (xas_tdp_unit, FMT="(T3,A)") &
                  " Index     Excitation energy (eV)    fosc dipole (a.u.)"
               DO i = 1, SIZE(donor_state%sg_evals)
                  WRITE (xas_tdp_unit, FMT="(T3,I6,F27.6,F22.6)") &
                     i, donor_state%sg_evals(i)*evolt, donor_state%osc_str(i)
               END DO
            END IF

            WRITE (xas_tdp_unit, FMT="(A,/)") " "
         END IF !xas_tdp_unit

         IF (output_unit > 0) THEN
            WRITE (output_unit, FMT="(T5,A,F25.6)") &
               "First singlet XAS excitation energy (eV): ", donor_state%sg_evals(1)*evolt
         END IF
      END IF ! do_singlet

      IF (xas_tdp_control%do_triplet) THEN
         IF (xas_tdp_unit > 0) THEN

!           Printing the general donor state information
            WRITE (xas_tdp_unit, FMT="(A,/,A,A,A/,A,I5,A,I5,A,A,/,A)") &
               "==================================================================================", &
               "XAS TDP triplet excitations (no SOC) for DONOR STATE: ", &
               xas_tdp_env%state_type_char(donor_state%state_type), ",", &
               "from EXCITED ATOM: ", donor_state%at_index, ", of KIND (index/symbol): ", &
               donor_state%kind_index, "/", TRIM(donor_state%at_symbol), &
               "=================================================================================="

!           Simply dump the excitation energies/ oscillator strength as they come

            IF (xas_tdp_control%do_quad) THEN
               WRITE (xas_tdp_unit, FMT="(T3,A)") &
                  " Index     Excitation energy (eV)    fosc dipole (a.u.)   fosc quadrupole (a.u.)"
               DO i = 1, SIZE(donor_state%tp_evals)
                  WRITE (xas_tdp_unit, FMT="(T3,I6,F27.6,F22.6,F25.6)") &
                     i, donor_state%tp_evals(i)*evolt, 0.0_dp, 0.0_dp !spin-forbidden !
               END DO
            ELSE
               WRITE (xas_tdp_unit, FMT="(T3,A)") &
                  " Index     Excitation energy (eV)    fosc dipole (a.u.)"
               DO i = 1, SIZE(donor_state%tp_evals)
                  WRITE (xas_tdp_unit, FMT="(T3,I6,F27.6,F22.6)") &
                     i, donor_state%tp_evals(i)*evolt, 0.0_dp
               END DO
            END IF

            WRITE (xas_tdp_unit, FMT="(A,/)") " "
         END IF !xas_tdp_unit

         IF (output_unit > 0) THEN
            WRITE (output_unit, FMT="(T5,A,F25.6)") &
               "First triplet XAS excitation energy (eV): ", donor_state%tp_evals(1)*evolt
         END IF
      END IF ! do_triplet

      IF (xas_tdp_control%do_soc) THEN
         IF (xas_tdp_unit > 0) THEN

!           Printing the general donor state information
            WRITE (xas_tdp_unit, FMT="(A,/,A,A,A/,A,I5,A,I5,A,A,/,A)") &
               "==================================================================================", &
               "XAS TDP  excitations after spin-orbit coupling for DONOR STATE: ", &
               xas_tdp_env%state_type_char(donor_state%state_type), ",", &
               "from EXCITED ATOM: ", donor_state%at_index, ", of KIND (index/symbol): ", &
               donor_state%kind_index, "/", TRIM(donor_state%at_symbol), &
               "=================================================================================="

!           Simply dump the excitation energies/ oscillator strength as they come
            IF (xas_tdp_control%do_quad) THEN
               WRITE (xas_tdp_unit, FMT="(T3,A)") &
                  " Index     Excitation energy (eV)    fosc dipole (a.u.)   fosc quadrupole (a.u.)"
               DO i = 1, SIZE(donor_state%soc_evals)
                  WRITE (xas_tdp_unit, FMT="(T3,I6,F27.6,F22.6,F25.6)") &
                     i, donor_state%soc_evals(i)*evolt, donor_state%soc_osc_str(i), &
                     donor_state%soc_quad_osc_str(i)
               END DO
            ELSE
               WRITE (xas_tdp_unit, FMT="(T3,A)") &
                  " Index     Excitation energy (eV)    fosc dipole (a.u.)"
               DO i = 1, SIZE(donor_state%soc_evals)
                  WRITE (xas_tdp_unit, FMT="(T3,I6,F27.6,F22.6)") &
                     i, donor_state%soc_evals(i)*evolt, donor_state%soc_osc_str(i)
               END DO
            END IF

            WRITE (xas_tdp_unit, FMT="(A,/)") " "
         END IF !xas_tdp_unit

         IF (output_unit > 0) THEN
            WRITE (output_unit, FMT="(T5,A,F25.6)") &
               "First SOC XAS excitation energy (eV): ", donor_state%soc_evals(1)*evolt
         END IF
      END IF !do_soc

      CALL cp_print_key_finished_output(xas_tdp_unit, logger, xas_tdp_section, "PRINT%SPECTRUM")

   END SUBROUTINE print_xas_tdp_to_file

! **************************************************************************************************
!> \brief Prints the donor_state and excitation_type info into a RESTART file for cheap PDOS and/or
!>        CUBE printing without expensive computation
!> \param ex_type singlet, triplet, etc.
!> \param donor_state ...
!> \param xas_tdp_section ...
! **************************************************************************************************
   SUBROUTINE write_donor_state_restart(ex_type, donor_state, xas_tdp_section)

      INTEGER, INTENT(IN)                                :: ex_type
      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(section_vals_type), POINTER                   :: xas_tdp_section

      CHARACTER(len=*), PARAMETER :: routineN = 'write_donor_state_restart', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=default_path_length)                 :: filename
      CHARACTER(len=default_string_length)               :: domo, excite, my_middle
      INTEGER                                            :: ex_atom, handle, nao, ndo_mo, nex, &
                                                            nspins, output_unit, rst_unit, &
                                                            state_type
      LOGICAL                                            :: do_print
      REAL(dp), DIMENSION(:), POINTER                    :: lr_evals
      TYPE(cp_fm_type), POINTER                          :: lr_coeffs
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: print_key

      NULLIFY (logger, lr_coeffs, lr_evals, print_key)

      !Initialization
      logger => cp_get_default_logger()
      do_print = .FALSE.
      IF (BTEST(cp_print_key_should_output(logger%iter_info, xas_tdp_section, &
                                           "PRINT%RESTART", used_print_key=print_key), cp_p_file)) do_print = .TRUE.

      IF (.NOT. do_print) RETURN

      CALL timeset(routineN, handle)

      output_unit = cp_logger_get_default_io_unit()

      !Get general info
      SELECT CASE (ex_type)
      CASE (tddfpt_spin_cons)
         lr_evals => donor_state%sc_evals
         lr_coeffs => donor_state%sc_coeffs
         excite = "spincons"
         nspins = 2
      CASE (tddfpt_spin_flip)
         lr_evals => donor_state%sf_evals
         lr_coeffs => donor_state%sf_coeffs
         excite = "spinflip"
         nspins = 2
      CASE (tddfpt_singlet)
         lr_evals => donor_state%sg_evals
         lr_coeffs => donor_state%sg_coeffs
         excite = "singlet"
         nspins = 1
      CASE (tddfpt_triplet)
         lr_evals => donor_state%tp_evals
         lr_coeffs => donor_state%tp_coeffs
         excite = "triplet"
         nspins = 1
      END SELECT

      SELECT CASE (donor_state%state_type)
      CASE (xas_1s_type)
         domo = "1s"
      CASE (xas_2s_type)
         domo = "2s"
      CASE (xas_2p_type)
         domo = "2p"
      END SELECT

      ndo_mo = donor_state%ndo_mo
      nex = SIZE(lr_evals)
      CALL cp_fm_get_info(lr_coeffs, nrow_global=nao)
      state_type = donor_state%state_type
      ex_atom = donor_state%at_index

      !Opening restart file
      rst_unit = -1
      my_middle = 'xasat'//TRIM(ADJUSTL(cp_to_string(ex_atom)))//'_'//TRIM(domo)//'_'//TRIM(excite)
      rst_unit = cp_print_key_unit_nr(logger, xas_tdp_section, "PRINT%RESTART", extension=".rst", &
                                      file_status="REPLACE", file_action="WRITE", &
                                      file_form="UNFORMATTED", middle_name=TRIM(my_middle))

      filename = cp_print_key_generate_filename(logger, print_key, middle_name=TRIM(my_middle), &
                                                extension=".rst", my_local=.FALSE.)

      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/T5,A,A,A)") &
            "Linear-response orbitals and excitation energies are written in: ", &
            '"', TRIM(filename), '"'
      END IF

      !Writing
      IF (rst_unit > 0) THEN
         WRITE (rst_unit) ex_atom, state_type, ndo_mo, ex_type
         WRITE (rst_unit) nao, nex, nspins
         WRITE (rst_unit) lr_evals(:)
      END IF
      CALL cp_fm_write_unformatted(lr_coeffs, rst_unit)

      !closing
      CALL cp_print_key_finished_output(rst_unit, logger, xas_tdp_section, "PRINT%RESTART")

      CALL timestop(handle)

   END SUBROUTINE write_donor_state_restart

! **************************************************************************************************
!> \brief Reads donor_state info from a restart file
!> \param donor_state the pre-allocated donor_state
!> \param ex_type the excitations stored in this specific file
!> \param filename the restart file to read from
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE read_donor_state_restart(donor_state, ex_type, filename, qs_env)

      TYPE(donor_state_type), POINTER                    :: donor_state
      INTEGER, INTENT(OUT)                               :: ex_type
      CHARACTER(len=*), INTENT(IN)                       :: filename
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'read_donor_state_restart', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: group, handle, nao, nex, nspins, &
                                                            output_unit, read_params(7), rst_unit, &
                                                            source
      LOGICAL                                            :: file_exists
      REAL(dp), DIMENSION(:), POINTER                    :: lr_evals
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: lr_coeffs
      TYPE(cp_para_env_type), POINTER                    :: para_env

      NULLIFY (lr_evals, lr_coeffs, para_env, fm_struct, blacs_env)

      CALL timeset(routineN, handle)

      output_unit = cp_logger_get_default_io_unit()
      CPASSERT(ASSOCIATED(donor_state))
      CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env)
      group = para_env%group
      source = para_env%source

      file_exists = .FALSE.
      rst_unit = -1

      IF (para_env%ionode) THEN

         INQUIRE (FILE=filename, EXIST=file_exists)
         IF (.NOT. file_exists) CPABORT("Trying to read non-existing XAS_TDP restart file")

         CALL open_file(file_name=TRIM(filename), file_action="READ", file_form="UNFORMATTED", &
                        file_position="REWIND", file_status="OLD", unit_number=rst_unit)
      END IF

      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/,T5,A,A,A)") &
            "Reading linear-response orbitals and excitation energies from file: ", &
            '"', filename, '"'
      END IF

      !read general params
      IF (rst_unit > 0) THEN
         READ (rst_unit) read_params(1:4)
         READ (rst_unit) read_params(5:7)
      END IF
      CALL mp_bcast(read_params, source, group)
      donor_state%at_index = read_params(1)
      donor_state%state_type = read_params(2)
      donor_state%ndo_mo = read_params(3)
      ex_type = read_params(4)
      nao = read_params(5)
      nex = read_params(6)
      nspins = read_params(7)

      !read evals
      ALLOCATE (lr_evals(nex))
      IF (rst_unit > 0) READ (rst_unit) lr_evals(1:nex)
      CALL mp_bcast(lr_evals, source, group)

      !read evecs
      CALL cp_fm_struct_create(fm_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=nao, ncol_global=nex*donor_state%ndo_mo*nspins)
      CALL cp_fm_create(lr_coeffs, fm_struct)
      CALL cp_fm_read_unformatted(lr_coeffs, rst_unit)
      CALL cp_fm_struct_release(fm_struct)

      !closing file
      IF (para_env%ionode) THEN
         CALL close_file(unit_number=rst_unit)
      END IF

      !case study on excitation type
      SELECT CASE (ex_type)
      CASE (tddfpt_spin_cons)
         donor_state%sc_evals => lr_evals
         donor_state%sc_coeffs => lr_coeffs
      CASE (tddfpt_spin_flip)
         donor_state%sf_evals => lr_evals
         donor_state%sf_coeffs => lr_coeffs
      CASE (tddfpt_singlet)
         donor_state%sg_evals => lr_evals
         donor_state%sg_coeffs => lr_coeffs
      CASE (tddfpt_triplet)
         donor_state%tp_evals => lr_evals
         donor_state%tp_coeffs => lr_coeffs
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE read_donor_state_restart

! **************************************************************************************************
!> \brief Checks whether this is a restart calculation and runs it if so
!> \param rst_filename the file to read for restart
!> \param xas_tdp_section ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE restart_calculation(rst_filename, xas_tdp_section, qs_env)

      CHARACTER(len=*), INTENT(IN)                       :: rst_filename
      TYPE(section_vals_type), POINTER                   :: xas_tdp_section
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'restart_calculation', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ex_type
      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env

      NULLIFY (xas_tdp_env, donor_state)

      !create a donor_state that we fill with the information we read
      ALLOCATE (donor_state)
      CALL donor_state_create(donor_state)
      CALL read_donor_state_restart(donor_state, ex_type, rst_filename, qs_env)

      !create a dummy xas_tdp_env and compute the post XAS_TDP stuff
      CALL xas_tdp_env_create(xas_tdp_env)
      CALL xas_tdp_post(ex_type, donor_state, xas_tdp_env, xas_tdp_section, qs_env)

      !clean-up
      CALL xas_tdp_env_release(xas_tdp_env)
      CALL free_ds_memory(donor_state)
      DEALLOCATE (donor_state)

   END SUBROUTINE restart_calculation

END MODULE xas_tdp_methods
