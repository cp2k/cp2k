!, unique------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Methods for X-Ray absorption spectroscopy (XAS) using TDDFPT
!> \author AB (11.2017)
! **************************************************************************************************

MODULE xas_tdp_methods
   USE arnoldi_api,                     ONLY: arnoldi_data_type, arnoldi_ev,&
                                              deallocate_arnoldi_data, get_selected_ritz_val, &
                                              setup_arnoldi_data
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE basis_set_types,                 ONLY: &
        allocate_sto_basis_set, create_gto_from_sto_basis, deallocate_gto_basis_set, &
        deallocate_sto_basis_set, get_gto_basis_set, gto_basis_set_type, init_orb_basis_set, &
        set_sto_basis_set, srules, sto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_array_utils,                  ONLY: cp_1d_r_p_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply, copy_dbcsr_to_fm, &
                                              copy_fm_to_dbcsr
   USE cp_fm_diag,                      ONLY: choose_eigv_solver, cp_fm_geeig
   USE cp_fm_struct,                    ONLY: cp_fm_struct_type, cp_fm_struct_create,&
                                              cp_fm_struct_release
   USE cp_fm_types,                     ONLY: cp_fm_type, cp_fm_get_submatrix, cp_fm_create,&
                                              cp_fm_release, cp_fm_set_all, cp_fm_p_type,&
                                              cp_fm_get_info, cp_fm_to_fm_submat, & 
                                              cp_fm_set_submatrix, cp_fm_get_diag, cp_fm_to_fm
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_log_handling,                 ONLY: cp_get_default_logger, cp_logger_type, &
                                              cp_logger_get_default_io_unit
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_output_handling,              ONLY: cp_print_key_unit_nr, debug_print_level, &
                                              cp_print_key_finished_output
   USE dbcsr_api,                       ONLY: dbcsr_p_type, dbcsr_type, dbcsr_release, dbcsr_print, &
                                              dbcsr_create, dbcsr_multiply, dbcsr_add_on_diag, &
                                              dbcsr_finalize, dbcsr_copy, dbcsr_set, dbcsr_add, &
                                              dbcsr_get_info
   USE hfx_libint_wrapper,              ONLY: cp_libint_static_init
   USE input_constants,                 ONLY: xas_tdp_by_index, xas_tdp_by_kind, xas_not_excited,&
                                              do_loc_none, xas_1s_type, xas_2s_type, xas_2p_type,&
                                              op_loc_berry, tddfpt_singlet, tddfpt_spin_cons, &
                                              tddfpt_triplet, xas_tdp_uks, xas_tdp_roks, &
                                              general_roks, high_spin_roks, state_loc_list, &
                                              tddfpt_spin_flip, do_hfx_potential_coulomb, &
                                              xas_dip_vel, xas_dip_len, do_hfx_potential_truncated,&
                                              do_hfX_potential_short
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_get_subs_vals,&
                                              section_vals_val_set, &
                                              section_vals_duplicate, &
                                              section_vals_release, &
                                              section_vals_get, &
                                              section_vals_val_get, &
                                              section_vals_get_subs_vals2
   USE kinds,                           ONLY: dp, default_string_length
   USE machine,                         ONLY: m_flush
   USE memory_utilities,                ONLY: reallocate
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE periodic_table,                  ONLY: ptable
   USE physcon,                         ONLY: angstrom,&
                                              evolt
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_loc_methods,                  ONLY: centers_spreads_berry,&
                                              qs_loc_driver
   USE qs_loc_types,                    ONLY: get_qs_loc_env,&
                                              localized_wfn_control_create,&
                                              localized_wfn_control_type,&
                                              qs_loc_env_create,&
                                              qs_loc_env_new_type,&
                                              qs_loc_env_release
   USE qs_loc_utils,                    ONLY: qs_loc_control_init,&
                                              qs_loc_env_init,&
                                              set_loc_centers,&
                                              set_loc_wfn_lists
   USE qs_mo_methods,                   ONLY: calculate_subspace_eigenvalues
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE qs_operators_ao,                 ONLY: p_xyz_ao, rRc_xyz_ao
   USE qs_rho_types,                    ONLY: qs_rho_get, qs_rho_type
   USE qs_scf_methods,                  ONLY: combine_ks_matrices
   USE qs_scf_types,                    ONLY: ot_method_nr
   USE scf_control_types,               ONLY: scf_control_type
   USE util,                            ONLY: locate, sort_unique
   USE xas_methods,                     ONLY: calc_stogto_overlap
   USE xas_tdp_atom,                    ONLY: init_xas_atom_env,&
                                              integrate_fxc_atoms,&
                                              integrate_soc_atoms
   USE xas_tdp_types,                   ONLY: xas_tdp_env_type,&
                                              xas_tdp_env_create,&
                                              xas_tdp_env_release,&
                                              set_xas_tdp_env,&
                                              xas_tdp_control_type, &
                                              read_xas_tdp_control, &
                                              xas_tdp_control_create, & 
                                              donor_state_type,&
                                              set_donor_state, &
                                              free_ds_memory, &
                                              xas_atom_env_type, &
                                              xas_atom_env_create, &
                                              xas_atom_env_release, &
                                              donor_state_create, &
                                              free_exat_memory
   USE xas_tdp_utils,                   ONLY: setup_xas_tdp_prob, &
                                              solve_xas_tdp_prob, &
                                              compute_ri_coulomb2_int, &
                                              compute_ri_exchange2_int, &
                                              include_rcs_soc, &
                                              include_os_soc, &
                                              compute_ri_o3c_coulomb, &
                                              compute_ri_o3c_exchange, &
                                              get_opt_3c_dist2d
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

      INTEGER                                            :: handle, ikind, iat, istate, tmp_index,&
                                                            current_state_index, output_unit, &
                                                            iatom
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      CHARACTER(len=2)                                   :: symbol_of_kind
      INTEGER :: current_state_index, current_state_index, handle, handle, iat, iat, iatom, ikind, &
         ikind, istate, istate, output_unit, output_unit, tmp_index, tmp_index
      INTEGER, DIMENSION(:), POINTER                     :: atoms_of_kind
      TYPE(donor_state_type), POINTER                    :: current_state
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: dft_section, xas_tdp_section
      TYPE(xas_atom_env_type), POINTER                   :: xas_atom_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      LOGICAL                                            :: do_os
      REAL(dp)                                           :: rc(3)
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      NULLIFY(xas_tdp_env, xas_tdp_control, atomic_kind_set, atoms_of_kind, current_state)
      NULLIFY(xas_atom_env, dft_control, particle_set)

      CALL timeset(routineN, handle)

!  Logger initialization and XAS TDP banner printing
      NULLIFY (logger, dft_section, xas_tdp_section)
      logger => cp_get_default_logger()
      dft_section => section_vals_get_subs_vals(qs_env%input, "DFT")
      xas_tdp_section => section_vals_get_subs_vals(dft_section, "XAS_TDP")
      output_unit = cp_print_key_unit_nr(logger, xas_tdp_section, "PRINT%PROGRAM_RUN_INFO", &
                                         extension=".Log")

      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T3,A,/,T3,A,/,T3,A,/,T3,A,/)") &
            "!===========================================================================!", &
            "!                              XAS TDP                                      !", &
            "!    Starting TDDFPT driven X-rays absorption spectroscopy calculations     !", &
            "!===========================================================================!"
      END IF

!  Initialization
      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T3,A)") &
            "# Create and initialize the XAS_TDP environment"
      END IF
      CALL get_qs_env(qs_env, dft_control=dft_control)
      CALL xas_tdp_init(xas_tdp_env, xas_tdp_control, qs_env)
      CALL print_info(output_unit, xas_tdp_control, qs_env)

      IF (output_unit > 0) THEN 
         IF (xas_tdp_control%check_only) THEN
            CPWARN("This is a CHECK_ONLY run for donor MOs (and possibly guess LUMOs) verification")
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
            WRITE(UNIT=output_unit, FMT="(/,T3,A,/,T5,A)") &
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

!  If CHECK_ONLY run, check the donor MOs and the guess LUMOs
      IF (xas_tdp_control%check_only) CALL print_checks(xas_tdp_env, xas_tdp_control, qs_env)

!  If not simply exact exchange, setup a xas_atom_env and compute the xc integrals on the atomic grids
!  Also needed if SOC is included. Done before looping on atoms as it's all done at once
      IF ((xas_tdp_control%do_xc .OR. xas_tdp_control%do_soc) .AND. .NOT. xas_tdp_control%check_only) THEN

         IF (output_unit > 0) THEN
            WRITE(UNIT=output_unit, FMT="(/,T3,A,I4,A)") &
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
      IF (.NOT. xas_tdp_control%check_only) THEN
         IF (output_unit > 0) THEN
            WRITE(UNIT=output_unit, FMT="(/,T3,A,I4,A)") &
               "# Computing the RI 3-center Coulomb integrals ..."
            CALL m_flush(output_unit)
         END IF
         CALL compute_ri_o3c_coulomb(xas_tdp_env, qs_env)
 
      END IF

!  Loop over donor states for calculation
      CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set)      
      current_state_index = 1

!     Loop over atomic kinds
      DO ikind = 1,SIZE(atomic_kind_set)

         IF (xas_tdp_control%check_only) EXIT
         
         CALL get_atomic_kind(atomic_kind=atomic_kind_set(ikind), element_symbol=symbol_of_kind,&
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
            IF (xas_tdp_control%dipole_form == xas_dip_len) THEN
               CALL get_qs_env(qs_env, particle_set=particle_set)
               rc = particle_set(iatom)%r
               CALL dbcsr_set(xas_tdp_env%dipmat(1)%matrix, 0.0_dp)
               CALL dbcsr_set(xas_tdp_env%dipmat(2)%matrix, 0.0_dp)
               CALL dbcsr_set(xas_tdp_env%dipmat(3)%matrix, 0.0_dp)
               CALL rRc_xyz_ao(xas_tdp_env%dipmat, qs_env, rc, order=1)
            END IF

            !compute the three center RI exchange integrals if so required
            IF (xas_tdp_control%do_hfx) THEN
               IF (output_unit > 0) THEN
                  WRITE(UNIT=output_unit, FMT="(/,T3,A,I4,A1,A2)") &
                     "# Computing the RI 3-center Exchange integrals for excited atom (idx/kind): ", &
                     iatom, "/", symbol_of_kind
                  CALL m_flush(output_unit)
               END IF
               CALL compute_ri_o3c_exchange(iatom, xas_tdp_env, xas_tdp_control, qs_env)
            END IF

!           Loop over states of excited atom of kind
            DO istate = 1, SIZE(xas_tdp_env%state_types, 1)

               IF (xas_tdp_env%state_types(istate, tmp_index) == xas_not_excited) CYCLE

               current_state => xas_tdp_env%donor_states(current_state_index)
               CALL set_donor_state(current_state, at_index=iatom, &
                                    at_symbol=symbol_of_kind, kind_index=ikind, &
                                    state_type=xas_tdp_env%state_types(istate, tmp_index))

!              Initial write for the donor state of interest
               IF (output_unit > 0) THEN
                  WRITE(UNIT=output_unit, FMT="(/,T3,A,A2,A,I4,A,A2,/)") &
                     "# Start of calculations for donor state of type ", &
                     xas_tdp_env%state_type_char(current_state%state_type), " for atom",&
                     current_state%at_index, " of kind ", current_state%at_symbol
                  CALL m_flush(output_unit)
               END IF

!              Assign best fitting MO(s) to current core donnor state
               CALL assign_mos_to_donor_state(current_state, xas_tdp_env, xas_tdp_control, qs_env)

!              Perform MO restricted Mulliken pop analysis for verification
               CALL perform_mulliken_on_donor_state(current_state,qs_env)

!              Do main calculations here
               CALL setup_xas_tdp_prob(current_state, qs_env, xas_tdp_env, xas_tdp_control)

               IF (xas_tdp_control%do_spin_cons) THEN
                  CALL solve_xas_tdp_prob(current_state, xas_tdp_control, xas_tdp_env, qs_env, &
                                          ex_type=tddfpt_spin_cons)
                  CALL compute_dipole_fosc(current_state, xas_tdp_control, xas_tdp_env)
               END IF
               IF (xas_tdp_control%do_spin_flip) THEN
                  CALL solve_xas_tdp_prob(current_state, xas_tdp_control, xas_tdp_env, qs_env, &
                                          ex_type=tddfpt_spin_flip)
                  !no dipole in spin-flip (spin-forbidden)
               END IF
               IF (xas_tdp_control%do_singlet) THEN
                  CALL solve_xas_tdp_prob(current_state, xas_tdp_control, xas_tdp_env, qs_env, &
                                          ex_type=tddfpt_singlet)
                  CALL compute_dipole_fosc(current_state, xas_tdp_control, xas_tdp_env)
               END IF
               IF (xas_tdp_control%do_triplet) THEN
                  CALL solve_xas_tdp_prob(current_state, xas_tdp_control, xas_tdp_env, qs_env, &
                                          ex_type=tddfpt_triplet)
                  !no dipole for triplets by construction
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

               CALL print_xas_tdp_to_file(current_state, xas_tdp_env, xas_tdp_control, xas_tdp_section)

!              Free some unneeded attributes of current_state
               CALL free_ds_memory(current_state)

               current_state_index = current_state_index + 1
               NULLIFY(current_state)
                  
            END DO ! state type

            CALL free_exat_memory(xas_tdp_env, iatom)
         END DO ! atom of kind 
      END DO ! kind

!  Clean-up
      CALL xas_tdp_env_release(xas_tdp_env)
      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T3,A,/,T3,A,/,T3,A,/)") &
            "!===========================================================================!", &
            "!     End of TDDFPT driven X-rays absorption spectroscopy calculations      !", &
            "!===========================================================================!"
      END IF

      CALL cp_print_key_finished_output(output_unit, logger, xas_tdp_section, &
                                        "PRINT%PROGRAM_RUN_INFO")

      CALL timestop(handle)

   END SUBROUTINE xas_tdp

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
      
      TYPE(section_vals_type), POINTER                   :: dft_section, xas_tdp_section,&
                                                            loc_section
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: at_kind_set
      INTEGER                                            :: nex_atoms, n_kinds, i, nat_of_kind,&
                                                            j, nex_kinds, k, n_donor_states,&
                                                            nspins, natom, nmatch, kind_ind, &
                                                            at_ind, nao, ispin
      INTEGER, DIMENSION(2)                              :: n_mo, n_moloc, homo
      INTEGER, DIMENSION(:), POINTER                     :: ind_of_kind
      CHARACTER(len=2)                                   :: symbol
      TYPE(qs_loc_env_new_type), POINTER                 :: qs_loc_env
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      REAL(dp), DIMENSION(:), POINTER                    :: mo_evals
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, rho_ao, matrix_ks
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(scf_control_type), POINTER                    :: scf_control
      REAL(dp)                                           :: fact
      LOGICAL                                            :: do_uks, do_os, unique
      TYPE(cell_type), POINTER                           :: cell

      NULLIFY(dft_section, xas_tdp_section, at_kind_set, ind_of_kind, qs_kind_set, dft_control)
      NULLIFY(qs_loc_env, loc_section, mos, particle_set, rho, rho_ao, mo_evals, cell)
      NULLIFY(mo_coeff, scf_control, matrix_ks)

!  XAS TDP control type initialization
      dft_section => section_vals_get_subs_vals(qs_env%input, "DFT")
      xas_tdp_section => section_vals_get_subs_vals(dft_section, "XAS_TDP")

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
            CALL get_atomic_kind(atomic_kind=at_kind_set(i), element_symbol=symbol, &
                                 natom=nat_of_kind, kind_number=kind_ind)
            IF (ANY(xas_tdp_control%list_ex_kinds == symbol)) THEN
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
            CALL get_atomic_kind(atomic_kind=at_kind_set(i), element_symbol=symbol, &
                                 natom=nat_of_kind, atom_list=ind_of_kind)
            DO j = 1, nex_kinds
               IF (xas_tdp_control%list_ex_kinds(j) == symbol) THEN
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

!  In case of ground state OT optimization, compute the MO eigenvalues => we need them. The arrays
!  are allocated anyways so the only cost is that of a MM
      IF (qs_env%scf_env%method == ot_method_nr) THEN

         CALL get_qs_env(qs_env, mos=mos, matrix_ks=matrix_ks)
         nspins = 1; IF (do_uks) nspins = 2

         DO ispin = 1, nspins
            CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff, eigenvalues=mo_evals)
            CALL calculate_subspace_eigenvalues(mo_coeff, matrix_ks(ispin)%matrix, evals_arg=mo_evals)
         END DO
      END IF

!  Initializing the qs_loc_env from the LOCALIZE subsection of XAS_TDP (largely inpired by MI's XAS)
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
         qs_loc_env%localized_wfn_control%lu_bound_states(1,2) = 1
         qs_loc_env%localized_wfn_control%lu_bound_states(2,2) = xas_tdp_control%n_search
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
      ALLOCATE(xas_tdp_env%mos_of_ex_atoms(xas_tdp_control%n_search, nex_atoms, nspins))

!  Compute the projector on the unoccupied, unperturbed ground state: Q = 1 - SP, sor each spin
      IF (do_os) nspins = 2
      CALL get_qs_env(qs_env, rho=rho, matrix_s=matrix_s)
      CALL qs_rho_get(rho, rho_ao=rho_ao)

      ALLOCATE(xas_tdp_env%q_projector(nspins))
      ALLOCATE(xas_tdp_env%q_projector(1)%matrix)
      CALL dbcsr_create(xas_tdp_env%q_projector(1)%matrix, name="Q PROJECTOR ALPHA", &
                        template=matrix_s(1)%matrix, matrix_type="N") 
      IF (do_os) THEN
         ALLOCATE(xas_tdp_env%q_projector(2)%matrix)
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
      ALLOCATE(xas_tdp_env%dipmat(3))
      DO i =1,3
         ALLOCATE(xas_tdp_env%dipmat(i)%matrix)          
         CALL dbcsr_copy(xas_tdp_env%dipmat(i)%matrix, matrix_s(1)%matrix, name="XAS TDP dipole matrix")
         CALL dbcsr_set(xas_tdp_env%dipmat(i)%matrix, 0.0_dp)
      END DO

!     Precompute it in the velocity representation, if so chosen
      IF (xas_tdp_control%dipole_form == xas_dip_vel) CALL p_xyz_ao(xas_tdp_env%dipmat, qs_env)

!  Allocate SOC in AO basis matrices
      IF (xas_tdp_control%do_soc) THEN
         ALLOCATE (xas_tdp_env%orb_soc(3))
         DO i = 1, 3
            ALLOCATE (xas_tdp_env%orb_soc(i)%matrix)
         END DO
      END IF

!  Computing the first LUMOs for guesses in case of OT eigensolver
      IF (.NOT. xas_tdp_control%full_diag) THEN
         CALL make_lumo_guess(xas_tdp_env, xas_tdp_control, qs_env)        
      END IF

!  Check that everything is allowed
      CALL safety_check(xas_tdp_control)

!  Initialize libint for the 3-center integrals
      CALL cp_libint_static_init()

!  Compute the optimized distribution_2d for the 3-center inegrals
      CALL get_opt_3c_dist2d(xas_tdp_env, qs_env)

   END SUBROUTINE xas_tdp_init

!> *************************************************************************************************
!> \brief Checks for forbidden keywords combinations
!> \param xas_tdp_control ...
!> *************************************************************************************************
   SUBROUTINE safety_check(xas_tdp_control)

      TYPE(xas_tdp_control_type), POINTER             :: xas_tdp_control

      CHARACTER(len=*), PARAMETER :: routineN = "safety_check", routineP = moduleN//":"//routineN

      !PB only available without exact exchange
      IF (xas_tdp_control%is_periodic .AND. xas_tdp_control%do_hfx &
          .AND. xas_tdp_control%x_op == do_hfx_potential_coulomb) THEN
         CPABORT("XAS TDP with Coulomb operator for exact exchange only supports non-periodic BCs")
      END IF

      !open-shell/closed-shell tests
      IF (xas_tdp_control%do_roks .OR. xas_tdp_control%do_uks) THEN

         IF (xas_tdp_control%do_singlet .OR. xas_tdp_control%do_triplet) THEN
            CPABORT("Singlet/triplet excitations only for restricted closed-shell systems")
         END IF

         IF (xas_tdp_control%do_soc .AND. .NOT. &
            (xas_tdp_control%do_spin_flip .AND. xas_tdp_control%do_spin_cons)) THEN

            CPABORT("Both spin-conserving and spin-flip excitations are required for SOC")
         END IF
      ELSE

         IF (xas_tdp_control%do_spin_cons .OR. xas_tdp_control%do_spin_flip) THEN
            CPABORT("Spin-conserving/spin-flip excitations only for open-shell systems")
         END IF

         IF (xas_tdp_control%do_soc .AND. .NOT. &
            (xas_tdp_control%do_singlet .AND. xas_tdp_control%do_triplet)) THEN

            CPABORT("Both singlet and triplet excitations are needed for SOC")
         END IF
      END IF

      !TDA, full-TDDFT and diagonalization
      IF (.NOT. xas_tdp_control%tamm_dancoff) THEN

         IF (xas_tdp_control%do_spin_flip) THEN
            CPABORT("Spin-flip kernel only implemented for Tamm-Dancoff approximation")
         END IF

         IF (.NOT. xas_tdp_control%full_diag) THEN
            CPABORT("Non-TDA (full TDDFT) calculation only available with stendard fm diagonalization")
         END IF

      END IF

   END SUBROUTINE safety_check

!> *************************************************************************************************
!> \brief Prints some basic info about the chosen parameters
!> \param ou the output unis
!> \param xas_tdp_control ...
!> \param qs_env ...
!> *************************************************************************************************
   SUBROUTINE print_info(ou, xas_tdp_control, qs_env)

      INTEGER, INTENT(IN)                             :: ou
      TYPE(xas_tdp_control_type), POINTER             :: xas_tdp_control
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = "print_info", routineP = moduleN//":"//routineN

      INTEGER                                         :: i, irep, nrep, il
      TYPE(section_vals_type), POINTER                :: input, xc_functionals, libxc_fun, xc_fun
      CHARACTER(LEN=20)                               :: tmpStr
      CHARACTER(LEN=10*default_string_length)         :: reference
      CHARACTER(LEN=2*default_string_length)          :: shortform

      NULLIFY(input, xc_functionals, libxc_fun, xc_fun)

      IF (ou .LE. 0) RETURN

      CALL get_qs_env(qs_env, input=input)

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
         WRITE(UNIT=ou, FMT="(T3,A)") &
            "XAS_TDP| Transition Dipole Representation: VELOCITY"
      ELSE
         WRITE(UNIT=ou, FMT="(T3,A)") &
            "XAS_TDP| Transition Dipole Representation: LENGTH"
      END IF

      !Grid info
      IF (xas_tdp_control%do_soc .OR. xas_tdp_control%do_xc) THEN
         WRITE(UNIT=ou, FMT="(T3,A)") &
            "XAS_TDP| Radial Grid(s) Info: Symbol,  na,  nr" 
         DO i = 1, SIZE(xas_tdp_control%grid_info,1)
            WRITE(UNIT=ou, FMT="(T3,A,A6,A,A,A,A)") &
               "                              ", TRIM(xas_tdp_control%grid_info(i,1)), ", ", &
               TRIM(xas_tdp_control%grid_info(i,2)), ", ", TRIM(xas_tdp_control%grid_info(i,3))
         END DO
      END IF

      !XC kernel
      IF (xas_tdp_control%do_xc) THEN

         WRITE(UNIT=ou, FMT="(/,T3,A,F5.2,A)") &
            "XAS_TDP| RI Region's Radius: ", xas_tdp_control%ri_radius, " Ang"

         CALL libxc_version_info(tmpStr)
         WRITE(UNIT=ou, FMT="(T3,A,A,/,T3,A,/)") &
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
         WRITE(UNIT=ou, FMT="(/,T3,A,/,/,T3,A,F5.3)") &
            "XAS_TDP| Exat Exchange Kernel: On ", &
            "EXACT_EXCHANGE| Scale: ", xas_tdp_control%sx
         IF (xas_tdp_control%x_op == do_hfx_potential_coulomb) THEN
            WRITE(UNIT=ou, FMT="(T3,A)") &
               "EXACT_EXCHANGE| Operator : Coulomb"
         ELSE IF (xas_tdp_control%x_op == do_hfx_potential_truncated) THEN
            WRITE(UNIT=ou, FMT="(T3,A,/,T3,A,F5.2,A,/,T3,A,A)") &
               "EXACT_EXCHANGE| Operator: Truncated Coulomb", &
               "EXACT_EXCHANGE| Range: ", xas_tdp_control%x_range*angstrom, ", (Ang)", &
               "EXACT_EXCHANGE| T_C_G_DATA: ", xas_tdp_control%x_t_c_filename
         ELSE IF (xas_tdp_control%x_op == do_hfx_potential_short) THEN
            WRITE(UNIT=ou, FMT="(T3,A,/,T3,A,F5.2,A,/,T3,A,F5.2,A,/,T3,A,ES7.1)") &
               "EXACT_EXCHANGE| Operator: Short Range", &
               "EXACT_EXCHANGE| Omega: ", xas_tdp_control%x_omega, ", (1/a0)", &
               "EXACT_EXCHANGE| Effective Range: ", xas_tdp_control%x_range*angstrom, ", (Ang)", &
               "EXACT_EXCHANGE| EPS_RANGE: ", xas_tdp_control%eps_range
         END IF
      END IF

   END SUBROUTINE print_info

!> *************************************************************************************************
!> \brief Computed the first N_EXCITED LUMOs for the OT eigensolver guesses
!> \param xas_tdp_env ...
!> \param xas_tdp_control ...
!> \param qs_env ...
!> \note Uses stendard diagonalization. Do not use the stendard make_lumo subroutine as it uses
!>       the OT eigensolver and there is no guarantee that it will converge fast
!> *************************************************************************************************
   SUBROUTINE make_lumo_guess(xas_tdp_env, xas_tdp_control, qs_env)

      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = "make_lumo_guess", routineP = moduleN//":"//routineN

      INTEGER                                            :: handle, nspins, ispin, nevals(2), &
                                                            nelec_spin(2), nao, minel, nocc(2)
      LOGICAL                                            :: do_os
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s 
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: lumo_struct, fm_struct
      TYPE(cp_fm_type), POINTER                          :: amatrix, bmatrix, work_fm, evecs
      REAL(dp), DIMENSION(:), ALLOCATABLE                :: evals

      NULLIFY(matrix_ks, matrix_s, para_env, blacs_env)
      NULLIFY(lumo_struct, amatrix, bmatrix, work_fm, fm_struct, evecs)

      CALL timeset(routineN, handle)

      do_os = xas_tdp_control%do_uks .OR. xas_tdp_control%do_roks
      nspins = 1; IF (do_os) nspins = 2 
      ALLOCATE(xas_tdp_env%lumo_guess(nspins))
      IF (xas_tdp_control%check_only) ALLOCATE(xas_tdp_env%lumo_evals(nspins))
      CALL get_qs_env(qs_env, matrix_ks=matrix_ks, matrix_s=matrix_s, nelectron_spin=nelec_spin, &
                      para_env=para_env, blacs_env=blacs_env)
      CALL dbcsr_get_info(matrix_s(1)%matrix, nfullrows_total=nao)
      IF (do_os) THEN
         nevals = nao - nelec_spin
         nocc = nelec_spin
      ELSE
         nevals = nao -nelec_spin(1)/2
         nocc = nelec_spin(1)/2
      END IF
      IF (xas_tdp_control%n_excited > 0 .AND. ANY(xas_tdp_control%n_excited < nevals)) THEN
         nevals = xas_tdp_control%n_excited/nspins !total # of evals is n_excited

         !Add some padding to make sure we have enough on both spins
         !Less electron for a given spin => more evals
         IF (do_os) THEN
            minel = MINLOC(nelec_spin, 1)
            nevals(minel) = nevals(minel) + ABS(nelec_spin(2)-nelec_spin(1))
         END IF
      END IF

      DO ispin = 1,nspins

         !Going through fm to diagonalize
         CALL cp_fm_struct_create(fm_struct, para_env=para_env, context=blacs_env, &
                                  nrow_global=nao, ncol_global=nao)
         CALL cp_fm_create(amatrix, fm_struct)
         CALL cp_fm_create(bmatrix, fm_struct)
         CALL cp_fm_create(evecs, fm_struct)
         CALL cp_fm_create(work_fm, fm_struct)
         ALLOCATE(evals(nocc(ispin)+nevals(ispin)))

         CALL copy_dbcsr_to_fm(matrix_ks(ispin)%matrix, amatrix)
         CALL copy_dbcsr_to_fm(matrix_s(1)%matrix, bmatrix)

         !The actual diagonalization through Cholesky decomposition
         CALL cp_fm_geeig(amatrix, bmatrix, evecs, evals, work_fm)

         !Storing results
         IF (xas_tdp_control%check_only) THEN
            ALLOCATE(xas_tdp_env%lumo_evals(ispin)%array(nevals(ispin)))
            xas_tdp_env%lumo_evals(ispin)%array(:) = evals(nocc(ispin)+1:nocc(ispin)+nevals(ispin))
         END IF

         CALL cp_fm_struct_create(lumo_struct, para_env=para_env, context=blacs_env, &
                                  nrow_global=nao, ncol_global=nevals(ispin))
         CALL cp_fm_create(xas_tdp_env%lumo_guess(ispin)%matrix, lumo_struct)

         CALL cp_fm_to_fm_submat(evecs, xas_tdp_env%lumo_guess(ispin)%matrix, nrow=nao, &
                                 ncol=nevals(ispin), s_firstrow=1, s_firstcol=nocc(ispin)+1, &
                                 t_firstrow=1, t_firstcol=1)

         !clean-up
         CALL cp_fm_release(amatrix)
         CALL cp_fm_release(bmatrix)
         CALL cp_fm_release(evecs)
         CALL cp_fm_release(work_fm)
         CALL cp_fm_struct_release(fm_struct)
         CALL cp_fm_struct_release(lumo_struct)
         DEALLOCATE(evals)
      END DO

      CALL timestop(handle)

   END SUBROUTINE make_lumo_guess

!> *************************************************************************************************
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
!> *************************************************************************************************
   SUBROUTINE assign_mos_to_ex_atoms(xas_tdp_env, xas_tdp_control, qs_env)

      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'assign_mos_to_ex_atoms',&
                                     routineP = moduleN//':'//routineN
      
      INTEGER, DIMENSION(:,:,:), POINTER                 :: mos_of_ex_atoms
      TYPE(localized_wfn_control_type), POINTER          :: localized_wfn_control
      INTEGER                                            :: iat, imo, n_search, nex_atoms, ispin, &
                                                            iat_memo, n_atoms, at_index, nspins
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      INTEGER, DIMENSION(3)                              :: perd_init
      INTEGER, DIMENSION(:, :), POINTER                  :: mos_of_ex_atoms
      REAL(dp)                                           :: dist, dist_min
      REAL(dp), DIMENSION(3)                             :: at_pos, r_ac, wfn_center
      TYPE(cell_type), POINTER                           :: cell
      TYPE(localized_wfn_control_type), POINTER          :: localized_wfn_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      NULLIFY (localized_wfn_control, mos_of_ex_atoms, cell, particle_set)

!  Initialization. mos_of_ex_atoms filled with -1, meaning no assigned state
      mos_of_ex_atoms => xas_tdp_env%mos_of_ex_atoms
      mos_of_ex_atoms(:,:,:) = -1
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
         DO imo = 1,n_search
!           retrieve MO wave function center coordinates. 
            wfn_center(1:3) = localized_wfn_control%centers_set(ispin)%array(1:3, imo) 
            iat_memo = 0

!           a large enough value to avoid bad surprises
            dist_min = 10000.0_dp
            DO iat = 1,n_atoms
               at_pos = particle_set(iat)%r
               r_ac = pbc(at_pos, wfn_center, cell)
               dist = SQRT(r_ac(1)*r_ac(1) + r_ac(2)*r_ac(2) + r_ac(3)*r_ac(3))
         
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

!> *************************************************************************************************
!> \brief Re-initialize the qs_loc_env to the current MOs.
!> \param qs_loc_env the env to re-initialize
!> \param n_loc_states the number of states to include
!> \param do_uks in cas of spin unrestricted calculation, initialize for both spins
!> \param qs_env ...
!> \note  Useful when one needs to make use of qs_loc features and it is either with canonical MOs
!>        or the localized MOs have been modified. do_localize is overwritten. 
!>        Same loc range for both spins
!> *************************************************************************************************
   SUBROUTINE reinit_qs_loc_env(qs_loc_env, n_loc_states, do_uks, qs_env)

      TYPE(qs_loc_env_new_type), POINTER              :: qs_loc_env
      INTEGER, INTENT(IN)                             :: n_loc_states
      LOGICAL, INTENT(IN)                             :: do_uks
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'reinit_qs_loc_env', &
         routineP = moduleN//':'//routineN

      INTEGER                                         :: i, nspins
      TYPE(localized_wfn_control_type), POINTER       :: loc_wfn_control

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
      loc_wfn_control%lu_bound_states(1,:) = 1
      loc_wfn_control%lu_bound_states(2,:) = n_loc_states
      loc_wfn_control%set_of_states = state_loc_list
      loc_wfn_control%do_homo = .TRUE.
      ALLOCATE(loc_wfn_control%loc_states(n_loc_states,2))
      DO i = 1,n_loc_states
         loc_wfn_control%loc_states(i,:) = i
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

!> *************************************************************************************************
!> \brief Diagonalize the subset of previously localized MOs that are associated to each excited
!>        atoms. Updates the MO coeffs accordingly.
!> \param xas_tdp_env ...
!> \param xas_tdp_control ...
!> \param qs_env ...
!> \note  Needed because after localization, the MOs loose their identity (1s, 2s , 2p, etc)
!> *************************************************************************************************
   SUBROUTINE diagonalize_assigned_mo_subset(xas_tdp_env, xas_tdp_control, qs_env)

      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'diagonalize_assigned_mo_subset', &
         routineP = moduleN//':'//routineN

      INTEGER                                         :: iat, nlmo, nao, ilmo, i, ispin, nspins
      TYPE(mo_set_p_type), DIMENSION(:), POINTER      :: mos
      TYPE(cp_fm_type), POINTER                       :: mo_coeff
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_ks
      TYPE(cp_para_env_type), POINTER                 :: para_env
      TYPE(cp_blacs_env_type), POINTER                :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                :: lmo_struct, ks_struct
      TYPE(cp_fm_type), POINTER                       :: lmo_fm, ks_fm, work, evecs
      REAL(dp), DIMENSION(:), ALLOCATABLE             :: evals

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
            nlmo = COUNT(xas_tdp_env%mos_of_ex_atoms(:,iat,ispin) == 1) 
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
            DO ilmo = 1,xas_tdp_control%n_search
               IF (xas_tdp_env%mos_of_ex_atoms(ilmo,iat,ispin) == -1) CYCLE 
            
               i = i+1
               ! put the coeff in our atom-restricted lmo_fm
               CALL cp_fm_to_fm_submat(mo_coeff, lmo_fm, nrow=nao, ncol=1, s_firstrow=1, &
                                       s_firstcol=ilmo, t_firstrow=1, t_firstcol=i)

            END DO !ilmo

         ! Computing the KS matrix in the subset of MOs
            CALL cp_dbcsr_sm_fm_multiply(matrix_ks(ispin)%matrix, lmo_fm, work, ncol=nlmo)
            CALL cp_gemm('T', 'N', nlmo, nlmo, nao, 1.0_dp, lmo_fm, work, 0.0_dp, ks_fm)
         
         ! Diagonalizing the KS matrix in the subset of MOs
            ALLOCATE(evals(nlmo))
            CALL choose_eigv_solver(ks_fm, evecs, evals)
            DEALLOCATE(evals)

         ! Express the MOs in the basis that diagonalizes KS
            CALL cp_gemm('N', 'N', nao, nlmo, nlmo, 1.0_dp, lmo_fm, evecs, 0.0_dp, work)

         ! Replacing the new MOs back in the MO coeffs
            i = 0
            DO ilmo = 1,xas_tdp_control%n_search
               IF (xas_tdp_env%mos_of_ex_atoms(ilmo,iat,ispin) == -1) CYCLE

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

!> *************************************************************************************************
!> \brief Assign core MO(s) to a given donor_state, taking the type (1S, 2S, etc) into account.
!>        The projection on a representative Slater-type orbital basis is used as a indicator.
!>        It is assumed that MOs are already assigned to excited atoms based on their center
!> \param donor_state the donor_state to which a MO must be assigned
!> \param xas_tdp_env ...
!> \param xas_tdp_control ...
!> \param qs_env ...
!> *************************************************************************************************
   SUBROUTINE assign_mos_to_donor_state(donor_state, xas_tdp_env, xas_tdp_control, qs_env)

      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(xas_tdp_env_type), POINTER                 :: xas_tdp_env        
      TYPE(xas_tdp_control_type), POINTER             :: xas_tdp_control
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN ='assign_mos_to_donor_state', &
                                                      routineP = moduleN//':'//routineN

      INTEGER                                         :: zval, i, l, nj, nsgf_sto, nsgf_kind, &
                                                         n_search, imo, iat, nao, at_index, my_mo,&
                                                         n_states, output_unit, ispin, nspins, ndo_so
      REAL(dp), DIMENSION(:), POINTER                 :: zeta
      INTEGER, DIMENSION(:), POINTER                  :: lq, nq, first_sgf
      TYPE(sto_basis_set_type), POINTER               :: sto_basis_set
      TYPE(gto_basis_set_type), POINTER               :: sto_to_gto_basis_set, kind_basis_set
      INTEGER, DIMENSION(4,7)                         :: ne
      INTEGER, DIMENSION(:,:), ALLOCATABLE            :: my_mos
      REAL(dp)                                        :: zeff 
      REAL(dp), DIMENSION(2)                          :: next_best_overlap
      INTEGER, DIMENSION(2)                           :: next_best_overlap_ind 
      TYPE(qs_kind_type), DIMENSION(:), POINTER       :: qs_kind_set
      TYPE(mo_set_p_type), DIMENSION(:), POINTER      :: mos
      TYPE(cp_fm_type), POINTER                       :: mo_coeff, work_mat, eval_mat, gs_coeffs
      INTEGER, DIMENSION(:,:,:), POINTER              :: mos_of_ex_atoms
      REAL(dp), DIMENSION(:,:), POINTER               :: tmp_coeff, overlap_matrix
      TYPE(particle_type), DIMENSION(:), POINTER      :: particle_set
      REAL(dp), DIMENSION(:), POINTER                 :: mo_evals
      REAL(dp), DIMENSION(:), ALLOCATABLE             :: overlap, sto_overlap, diag
      REAL(dp), DIMENSION(:,:), ALLOCATABLE           :: max_overlap
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_ks
      TYPE(cp_fm_struct_type), POINTER                :: eval_mat_struct, gs_struct 
      TYPE(cp_para_env_type), POINTER                 :: para_env                                    
      TYPE(cp_blacs_env_type), POINTER                :: blacs_env 
      
      NULLIFY(sto_basis_set, sto_to_gto_basis_set, qs_kind_set, kind_basis_set, lq, nq, zeta)
      NULLIFY(overlap_matrix, mos, mo_coeff, mos_of_ex_atoms, tmp_coeff, first_sgf, particle_set)
      NULLIFY(mo_evals, matrix_ks, para_env, blacs_env, work_mat)
      NULLIFY(eval_mat_struct, eval_mat, gs_struct, gs_coeffs)
      
      output_unit = cp_logger_get_default_io_unit()

      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set, mos=mos, particle_set=particle_set,&
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
      ALLOCATE(my_mos(n_states, nspins))
      ALLOCATE(max_overlap(n_states, nspins))

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
      ALLOCATE(first_sgf(SIZE(particle_set))) !probably do not need that
      CALL get_particle_set(particle_set=particle_set, qs_kind_set=qs_kind_set, first_sgf=first_sgf)
      ALLOCATE(tmp_coeff(nsgf_kind,1))
      ALLOCATE(sto_overlap(nsgf_kind))
      ALLOCATE(overlap(n_search))

      next_best_overlap = 0.0_dp
      max_overlap = 0.0_dp

      DO ispin = 1, nspins

         CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff, nao=nao)
         overlap = 0.0_dp

         my_mo = 0
         DO imo = 1,n_search
            IF(mos_of_ex_atoms(imo,iat,ispin) > 0) THEN

               sto_overlap = 0.0_dp
               tmp_coeff = 0.0_dp

!              Getting the relevant coefficients for the candidate state
               CALL cp_fm_get_submatrix(fm=mo_coeff, target_m=tmp_coeff, start_row=first_sgf(at_index),&
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
         DO i = 1,n_states
            my_mo = MAXLOC(overlap,1)
            my_mos(i,ispin) = my_mo 
            max_overlap(i,ispin) = MAXVAL(overlap,1)
            overlap(my_mo) = 0.0_dp
         END DO
!        Getting the next best overlap (for validation purposes)
         next_best_overlap(ispin) = MAXVAL(overlap,1) 
         next_best_overlap_ind(ispin) = MAXLOC(overlap,1)

      END DO !ispin

!     Some clean-up
      DEALLOCATE(overlap_matrix, tmp_coeff) 

!  Dealing with the result
      IF (ALL(my_mos > 0) .AND. ALL(my_mos <= n_search)) THEN
!        Assigning the MO indices to the donor_state
         ALLOCATE(donor_state%mo_indices(n_states, nspins))
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
                                       ncol=1, s_firstrow=1, s_firstcol=my_mos(i,ispin), &
                                       t_firstrow=1, t_firstcol=(ispin-1)*n_states+i)
            END DO
         END DO
         gs_coeffs => donor_state%gs_coeffs

         !Keep the subset of the coeffs centered on the excited atom as global array (used a lot)
         ALLOCATE(donor_state%contract_coeffs(nsgf_kind, n_states*nspins))
         CALL cp_fm_get_submatrix(gs_coeffs, donor_state%contract_coeffs, start_row=first_sgf(at_index),&
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

            ALLOCATE(donor_state%energy_evals(n_states, nspins)) 
            donor_state%energy_evals = 0.0_dp

!           Canonical MO, no change in eigenvalues, only diagonal elements
            DO ispin = 1, nspins
               CALL get_mo_set(mos(ispin)%mo_set, eigenvalues=mo_evals)
               DO i = 1,n_states 
                  donor_state%energy_evals(i,ispin) = mo_evals(my_mos(i,ispin))
               
                  IF (output_unit > 0 ) THEN
                     WRITE (UNIT=output_unit, FMT="(T46,I4,I11,F17.5)") &
                        ispin, my_mos(i,ispin), max_overlap(i,ispin)
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
            DO ispin = 1,nspins
               DO i = 1,n_states

!                 Print info
                  IF (output_unit > 0 ) THEN
                     WRITE (UNIT=output_unit, FMT="(T46,I4,I11,F17.5)") &
                        ispin, my_mos(i,ispin), max_overlap(i,ispin) 
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
            ALLOCATE(diag(ndo_so))

            IF (.NOT. xas_tdp_control%do_roks) THEN

               ALLOCATE(donor_state%energy_evals(n_states, nspins)) 
               donor_state%energy_evals = 0.0_dp

!              Compute gs_coeff^T * matrix_ks * gs_coeff to get the epsilon_ij matrix
               DO ispin = 1, nspins
                  CALL cp_dbcsr_sm_fm_multiply(matrix_ks(ispin)%matrix, gs_coeffs, work_mat, ncol=ndo_so)
                  CALL cp_gemm('T', 'N', ndo_so, ndo_so, nao, 1.0_dp, gs_coeffs, work_mat, 0.0_dp, eval_mat)

!                 Put the epsilon_ii into the donor_state. No off-diagonal element because of subset diag
                  CALL cp_fm_get_diag(eval_mat, diag)
                  donor_state%energy_evals(:,ispin) = diag((ispin-1)*n_states+1:ispin*n_states)
                  
               END DO

            ELSE
            ! If ROKS, slightly different procedure => 2 KS matrices but one type of MOs
               ALLOCATE(donor_state%energy_evals(n_states, 2)) 
               donor_state%energy_evals = 0.0_dp

!              Compute gs_coeff^T * matrix_ks * gs_coeff to get the epsilon_ij matrix
               DO ispin = 1, 2
                  CALL cp_dbcsr_sm_fm_multiply(matrix_ks(ispin)%matrix, gs_coeffs, work_mat, ncol=ndo_so)
                  CALL cp_gemm('T', 'N', ndo_so, ndo_so, nao, 1.0_dp, gs_coeffs, work_mat, 0.0_dp, eval_mat)

                  CALL cp_fm_get_diag(eval_mat, diag)
                  donor_state%energy_evals(:,ispin) = diag(:)

               END DO

               DEALLOCATE(diag)
            END IF


!           Clean-up
            CALL cp_fm_release(work_mat)
            CALL cp_fm_release(eval_mat)
            CALL cp_fm_struct_release(eval_mat_struct)

         END IF ! do_localize and/or ROKS

!        Clean-up
         CALL cp_fm_struct_release(gs_struct) 
         DEALLOCATE(first_sgf)

         IF (output_unit > 0) WRITE(UNIT=output_unit, FMT="(T5,A)") " "

         DO ispin = 1, nspins
            IF (output_unit > 0) THEN
               WRITE(UNIT=output_unit, FMT="(T5,A,I1,A,F7.5,A,I4)") &
                  "The next best overlap for spin ", ispin, " is ", next_best_overlap(ispin), &
                  " for MO with index ", next_best_overlap_ind(ispin)
            END IF
         END DO
         IF (output_unit > 0) WRITE(UNIT=output_unit, FMT="(T5,A)") " "

      ELSE
         CPABORT("A core donor state could not be assigned MO(s). Increasing NSEARCH might help.")
      END IF

   END SUBROUTINE assign_mos_to_donor_state

!> *************************************************************************************************
!> \brief Compute the centers and spreads of (core) MOs using the Berry phase operator
!> \param centers the coordinates of the MO centers and spreads. First 3 elements of each coloum
!>        contain the MO centers x,y,z coordinates, last 2 contain the spread for 2 functionals
!> \param xas_tdp_env ...
!> \param xas_tdp_control ...
!> \param qs_env ...
!> \note xas_tdp_env%qs_loc_env is used and modified. OK since no localization done after this
!>       subroutine is used.
!> *************************************************************************************************
   SUBROUTINE find_mo_centers(xas_tdp_env, xas_tdp_control, qs_env)

      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(qs_environment_type), POINTER                 :: qs_env

      INTEGER                                         :: n_centers, nao, i, j, dim_op, nspins, ispin
      TYPE(qs_loc_env_new_type), POINTER              :: qs_loc_env
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER    :: zij_fm_set 
      TYPE(cell_type), POINTER                        :: cell 
      REAL(dp), DIMENSION(6)                          :: weights 
      TYPE(section_vals_type), POINTER                :: print_loc_section, prog_run_info
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER    :: op_sm_set 
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER       :: moloc_coeff
      TYPE(cp_fm_type), POINTER                       :: vectors, opvec
      TYPE(cp_fm_struct_type), POINTER                :: tmp_fm_struct  
      TYPE(cp_para_env_type), POINTER                 :: para_env                                    
      TYPE(cp_blacs_env_type), POINTER                :: blacs_env

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
      ALLOCATE(zij_fm_set(2,dim_op))
      DO i = 1,dim_op
         DO j = 1,2
            CALL cp_fm_create(zij_fm_set(j, i)%matrix, tmp_fm_struct)
         END DO
      END DO

 !  If spin-unrestricted, need to go spin by spin
      nspins = 1; IF (xas_tdp_control%do_uks) nspins = 2

      DO ispin = 1,nspins                
!     zij computation, copied from qs_loc_methods:optimize_loc_berry
         vectors => moloc_coeff(ispin)%matrix
         DO i = 1,dim_op
            DO j = 1,2 
               CALL cp_fm_set_all(zij_fm_set(j, i)%matrix, 0.0_dp)
               CALL cp_dbcsr_sm_fm_multiply(op_sm_set(j, i)%matrix, vectors, opvec, ncol=n_centers) 
               CALL cp_gemm("T", "N", n_centers, n_centers, nao, 1.0_dp, vectors, opvec, 0.0_dp, & 
                            zij_fm_set(j, i)%matrix)
            END DO
         END DO

!     Compute centers (and spread)
         CALL centers_spreads_berry(qs_loc_env=qs_loc_env, zij=zij_fm_set, nmoloc=n_centers,&
                                    cell=cell, weights=weights, ispin=ispin ,&
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
      DEALLOCATE(zij_fm_set)

!  Make sure we leave with the correct do_loc value
      qs_loc_env%do_localize = xas_tdp_control%do_loc

   END SUBROUTINE find_mo_centers

! **************************************************************************************************
!> \brief Prints the MO to donor_state assocaition with overlap and Mulliken population analysis, 
!>        as well as the gaps in the LUMO eigenvalues in cas of OT solver. 
!> \param xas_tdp_env ...
!> \param xas_tdp_control ...
!> \param qs_env ...
!> \note  Called only in case of CHECK_ONLY run
! **************************************************************************************************
   SUBROUTINE print_checks(xas_tdp_env, xas_tdp_control, qs_env)   

      TYPE(xas_tdp_env_type), POINTER                 :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER             :: xaS_tdp_control
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = "print_checks", routineP = moduleN//":"//routineN

      INTEGER                                         :: current_state_index, ikind, iat, iatom, &
                                                         tmp_index, istate, output_unit
      TYPE(atomic_kind_type), DIMENSION(:), POINTER   :: atomic_kind_set
      INTEGER, DIMENSION(:), POINTER                  :: atoms_of_kind
      CHARACTER(len=2)                                :: symbol_of_kind
      TYPE(donor_state_type), POINTER                 :: current_state

      NULLIFY(atomic_kind_set, atoms_of_kind, current_state)

      output_unit = cp_logger_get_default_io_unit()

      IF (output_unit > 0) THEN
         WRITE(output_unit, "(/,T3,A,/,T3,A,/,T3,A)") &
            "# Check the donor states for their quality. They need to have a well defined type ", &
            "  (1s, 2s, etc) which is indicated by the overlap. They also need to be localized, ", &
            "  for which the Mulliken population analysis is one indicator (must be close to 1.0)"
      END IF

!  Loop over the donor states (as in the main xas_tdp loop)
      CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set)
      current_state_index = 1

      !loop over atomic kinds
      DO ikind = 1,SIZE(atomic_kind_set)

         CALL get_atomic_kind(atomic_kind=atomic_kind_set(ikind), element_symbol=symbol_of_kind,&
                              atom_list=atoms_of_kind)

         IF (.NOT. ANY(xas_tdp_env%ex_kind_indices == ikind)) CYCLE

         !loop over atoms of kind
         DO iat = 1,SIZE(atoms_of_kind)
            iatom = atoms_of_kind(iat)

            IF (.NOT. ANY(xas_tdp_env%ex_atom_indices == iatom)) CYCLE
            tmp_index = locate(xas_tdp_env%ex_atom_indices, iatom)

            !loop over states of excited atom
            DO istate = 1,SIZE(xas_tdp_env%state_types,1)

               IF (xas_tdp_env%state_types(istate, tmp_index) == xas_not_excited) CYCLE

               current_state => xas_tdp_env%donor_states(current_state_index)
               CALL set_donor_state(current_state, at_index=iatom,&
                                    at_symbol=symbol_of_kind, kind_index=ikind,&
                                    state_type=xas_tdp_env%state_types(istate,tmp_index))

               IF (output_unit > 0) THEN
                  WRITE(output_unit, "(/,T4,A,A2,A,I4,A,A2,A)") &
                     "-Donor state of type ", xas_tdp_env%state_type_char(current_state%state_type),&
                     " for atom", current_state%at_index, " of kind ", current_state%at_symbol, ":"
               END IF

               !Assign the MOs and perform Mulliken
               CALL assign_mos_to_donor_state(current_state, xas_tdp_env, xas_tdp_control, qs_env)
               CALL perform_mulliken_on_donor_state(current_state,qs_env)

               current_state_index = current_state_index + 1
               NULLIFY(current_state)

            END DO !istate
         END DO !iat
      END DO !ikind

      IF (output_unit > 0) THEN
         WRITE(output_unit, "(/,T5,A)") &
            "Use LOCALIZE and/or increase N_SEARCH for better results, if so required."
      END IF

!  If OT eigensolver, then also analyse the LUMOs
      IF (.NOT. xas_tdp_control%full_diag) CALL check_lumo_guess(xas_tdp_control, xas_tdp_env, qs_env)

   END SUBROUTINE print_checks

! **************************************************************************************************
!> \brief Analyze the LUMOs from the ground state calculations to be used as guess for the LR
!>        orbitals. Mainly look for gaps in the spctrum such that OT converges optimally
!> \param xas_tdp_control ...
!> \param xas_tdp_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE check_lumo_guess(xas_tdp_control, xas_tdp_env, qs_env)

      TYPE(xas_tdp_control_type), POINTER             :: xas_tdp_control
      TYPE(xas_tdp_env_type), POINTER                 :: xas_tdp_env
      TYPE(qs_environment_type), POINTER              :: qs_env
      
      CHARACTER(len=*), PARAMETER :: routineN = "check_lumo_guess", routineP = moduleN//":"//routineN

      INTEGER                                         :: nlumo(2), imo, output_unit, ilumo, nspins, &
                                                         offset, nelec_spin(2), ispin
      LOGICAL                                         :: do_os
      REAL(dp)                                        :: min_gap
      REAL(dp), DIMENSION(:), ALLOCATABLE             :: ref
      TYPE(cp_1d_r_p_type), DIMENSION(:), POINTER     :: lumo_evals
      TYPE(mo_set_p_type), DIMENSION(:), POINTER      :: mos
      REAL(dp), DIMENSION(:), POINTER                 :: mo_evals
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_ks, matrix_s, matrices
      TYPE(arnoldi_data_type)                         :: my_arnoldi

      NULLIFY(lumo_evals, mos, mo_evals, matrix_ks, matrix_s, matrices)

      output_unit = cp_logger_get_default_io_unit()

!  Get the LUMOs and their energies
      do_os = xas_tdp_control%do_uks .OR. xas_tdp_control%do_roks
      lumo_evals => xas_tdp_env%lumo_evals
      nspins = SIZE(lumo_evals)
      nlumo(1) = SIZE(lumo_evals(1)%array)
      IF (do_os) nlumo(2) = SIZE(lumo_evals(2)%array)
      CALL get_qs_env(qs_env, mos=mos, nelectron_spin=nelec_spin, matrix_ks=matrix_ks, &
                      matrix_s=matrix_s)

      IF (output_unit > 0) THEN
         WRITE(output_unit, "(/,T3,A,/,T3,A,/,T3,A,/,T3,A)") &
            "# Check the LUMOs eigenvalues for gaps, as the OT iterative eigensolver", &
            "  takes LUMOs as guesses and converges well only in presence of a gap.", &
            "  The LUMO energies below have the (spin) corresponding lowest MO ", &
            "  eigenvalue as reference (plus a shift to make it more readable)."
      END IF

!  Set the minimal gap to an arbitrary 5.0 eV = 0.18 Ha
      min_gap = 5.0_dp/evolt

!  Get the energy reference
      ALLOCATE(ref(nspins))
      ref = 0.0_dp
      CALL get_mo_set(mos(1)%mo_set, eigenvalues=mo_evals)
      ref(1) = mo_evals(1)
      IF (nspins == 2) THEN
         CALL get_mo_set(mos(2)%mo_set, eigenvalues=mo_evals)
         ref(2) = mo_evals(1)
      END IF
      ref = ref - MINVAL(ref, 1)
      
      !if ROKS, the KS eigenvalues are not physical => need MO^T * KS * MO
      IF (xas_tdp_control%do_roks) THEN
         ALLOCATE(matrices(2))
         matrices(2)%matrix => matrix_s(1)%matrix
         DO ispin = 1, 2
            matrices(1)%matrix => matrix_ks(ispin)%matrix

            CALL setup_arnoldi_data(my_arnoldi, matrices, max_iter=50, threshold=1.0E-4_dp, &
                                    selection_crit=3, nval_request=1, nrestarts=8, &
                                    generalized_ev=.TRUE., iram=.FALSE.)
            CALL arnoldi_ev(matrices, my_arnoldi)
            ref(ispin) = REAL(get_selected_ritz_val(my_arnoldi, 1), dp)
            CALL deallocate_arnoldi_data(my_arnoldi)

         END DO
         DEALLOCATE(matrices)
         ref = ref - MINVAL(ref, 1)
      END IF
      
!  Two different approaches depending on calculation type
!  If spin-restricted
      IF (.NOT. xas_tdp_control%do_uks .AND. .NOT. xas_tdp_control%do_roks) THEN
         DO imo = 1, nlumo(1)-1
            IF (min_gap <= lumo_evals(1)%array(imo+1)-lumo_evals(1)%array(imo)) THEN
               IF (output_unit > 0) THEN
                   WRITE(output_unit, "(/,T5,A22,I5,/,T8,A)") &
                        "Found gap after LUMO # ", imo, &
                        "LUMO index       Indicative energy (eV)"
               END IF
               DO ilumo = MAX(imo-1, 1), MIN(imo+1,nlumo(1))
                  IF (output_unit > 0) THEN
                        WRITE(output_unit, "(T8,I10,F29.3)") &
                           ilumo, (lumo_evals(1)%array(ilumo)-ref(1))*evolt
                  END IF
               END DO !ilumo
               IF (output_unit > 0) THEN
                  WRITE(output_unit, "(T5,A,I5)") &
                     "Corresponding N_EXCITED value: ", imo
               END IF
            END IF !gap
         END DO !imo
         IF (output_unit > 0) THEN
            WRITE(output_unit, "(/,T5,A,/,T5,A)") &
               "Use one of the above value as N_EXCITED for optimal OT convergence.", &
               "Increase current N_EXCITED and restart a CHECK_ONLY run for more options."
         END IF
      END IF !spin-unrestricted

!  If open-shell, ignore the singly occupied MOs, consider the two LUMO eigenvalue vector with the
!  appropriate reference (also with a shift)
      IF (do_os) THEN
         !assume more alpha electrons
         offset = nelec_spin(1)-nelec_spin(2)
         DO imo = 1, MINVAL(nlumo)-offset-1
            !check that both spin chanels have gaps
            IF ((min_gap <= lumo_evals(1)%array(imo+1)-lumo_evals(1)%array(imo)) .AND. &
                (min_gap <= lumo_evals(2)%array(imo+offset+1)-lumo_evals(2)%array(imo+offset))) THEN
               
               IF (output_unit > 0) THEN
                  WRITE(output_unit, "(/,T5,A,I5,A,I5,A,/,T8,A)") &
                     "Found gaps after LUMO # ", imo, " for spin 1     and LUMO #", imo+offset, &
                     " for spin 2", &
                     "LUMO index  Indicative energy (eV),       LUMO index  Indicative energy (eV)"
               END IF
               DO ilumo = MAX(imo-1, 1), MIN(imo+1,MINVAL(nlumo))
                  IF (output_unit > 0) THEN
                     WRITE(output_unit, "(T8,I10,F24.3,I18,F24.3)") &
                        !different ref 
                        ilumo, (lumo_evals(1)%array(ilumo)-ref(1))*evolt, &
                        ilumo+offset, (lumo_evals(2)%array(ilumo+offset)-ref(2))*evolt
                  END IF
               END DO !ilumo
               IF (output_unit > 0) THEN
                  WRITE(output_unit, "(T5,A,I5)") &
                     "Corresponding N_EXCITED value: ", 2*imo+offset
               END IF
            END IF !gap
         END DO !imo
         IF (output_unit > 0) THEN
            WRITE(output_unit, "(/,T5,A,/,T5,A)") &
               "Use one of the above value as N_EXCITED (or -1) for optimal OT convergence.", &
               "Increase current N_EXCITED and restart a CHECK_ONLY run for more options."
         END IF
      END IF

   END SUBROUTINE check_lumo_guess

! **************************************************************************************************
!> \brief Computes the oscillator strength based on the dipole moment (velocity rep) for
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

      INTEGER                                         :: nosc, iosc, ndo_mo, j, nao, nspins, ispin,&
                                                         ndo_so, is, handle
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: dipmat
      TYPE(cp_fm_type), POINTER                       :: mat_work, lr_coeffs
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER       :: col_work
      TYPE(cp_fm_struct_type), POINTER                :: col_struct, mat_struct
      REAL(dp), DIMENSION(:), ALLOCATABLE             :: diag, tot_contr
      REAL(dp), DIMENSION(:), POINTER                 :: lr_evals, osc_str
      TYPE(cp_para_env_type), POINTER                 :: para_env
      TYPE(cp_blacs_env_type), POINTER                :: blacs_env
      REAL(dp)                                        :: dip_mom, pref
      LOGICAL                                         :: do_sg, do_sc

      NULLIFY(dipmat, col_work, mat_work, col_struct, mat_struct, para_env, blacs_env, lr_coeffs)
      NULLIFY(lr_evals, osc_str)

      CALL timeset(routineN, handle)

!  Initialization
      do_sc = xas_tdp_control%do_spin_cons
      do_sg = xas_tdp_control%do_singlet
      IF (do_sc) THEN
         nspins = 2
         lr_evals => donor_state%sc_evals
         osc_str => donor_state%sc_osc_str
      ELSE IF (do_sg) THEN
         nspins = 1
         lr_evals => donor_state%sg_evals
         osc_str => donor_state%sg_osc_str
      ELSE
         CPABORT("Dipole oscilaltor strength only for singlets and spin-conserving excitations.")
      END IF
      ndo_mo = donor_state%ndo_mo
      ndo_so = ndo_mo 
      IF (xas_tdp_control%do_uks) ndo_so = ndo_mo*2 !both spins explicitely in the ground state
      nosc = SIZE(lr_evals)

      dipmat => xas_tdp_env%dipmat
      !compute the prefactor
      IF (xas_tdp_control%dipole_form == xas_dip_len) THEN
         osc_str(:) = 2.0_dp/3.0_dp*lr_evals(:)
      ELSE
         osc_str(:) = 2.0_dp/3.0_dp/lr_evals(:)
      END IF

      ! do some work matrix initialization
      CALL cp_fm_get_info(donor_state%gs_coeffs, matrix_struct=col_struct, para_env=para_env, &
                          context=blacs_env, nrow_global=nao)
      CALL cp_fm_struct_create(mat_struct, para_env=para_env, context=blacs_env, &
                               nrow_global=ndo_mo, ncol_global=ndo_mo)
      CALL cp_fm_create(mat_work, mat_struct)
      ALLOCATE(col_work(3))
      DO j = 1,3
         CALL cp_fm_create(col_work(j)%matrix, col_struct)
         ! get the product dipmat * c^0 for both spins (first columns for alpha)
         CALL cp_dbcsr_sm_fm_multiply(dipmat(j)%matrix, donor_state%gs_coeffs, col_work(j)%matrix, &
                                      ncol=ndo_so)
      END DO

      ALLOCATE(diag(ndo_mo), tot_contr(ndo_mo))
      tot_contr = 0.0_dp
      pref = 2.0_dp; IF (do_sc) pref = 1.0_dp !because of singlet definition u = 1/sqrt(2)(c_a+c_b)
    
!  Looping over the excited singlets
      DO iosc = 1, nosc

         dip_mom = 0.0_dp
         IF (do_sc) THEN
            lr_coeffs => donor_state%sc_coeffs(iosc)%matrix
         ELSE
            lr_coeffs => donor_state%sg_coeffs(iosc)%matrix
         END IF

         !looping over the cartesian coordinates
         DO j = 1, 3

            tot_contr = 0.0_dp
            DO ispin = 1,nspins
               is = ispin
               IF (.NOT. xas_tdp_control%do_uks) is = 1 !only have alpha MOs

               !computing c^i * dipmat(j) * c^0
               CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, lr_coeffs, col_work(j)%matrix, &
                            0.0_dp, mat_work, a_first_col=(ispin-1)*ndo_mo+1, & 
                            b_first_col=(is-1)*ndo_mo+1)

               CALL cp_fm_get_diag(mat_work, diag)
               tot_contr(:) = tot_contr(:) + diag(:)

            END DO !ispin

            !get the sum of the dipole moments. 
            dip_mom = dip_mom + pref*SUM(tot_contr)**2

         END DO !j

         !Taking the excitation energy into account
         osc_str(iosc) = osc_str(iosc)*dip_mom

      END DO !iosc

      !clean-up
      CALL cp_fm_release(mat_work)
      CALL cp_fm_struct_release(mat_struct)
      DO j = 1,3
         CALL cp_fm_release(col_work(j)%matrix)
      END DO
      DEALLOCATE(diag, col_work)

      CALL timestop(handle)

   END SUBROUTINE compute_dipole_fosc

!> *************************************************************************************************
!> \brief Writes the core MOs to excited atoms associations in the main output file
!> \param xas_tdp_env ...
!> \param qs_env ...
!> \note Look at alpha spin MOs, as we are dealing with core states and alpha/beta MOs are the same
!> *************************************************************************************************
   SUBROUTINE write_mos_to_ex_atoms_association(xas_tdp_env, xas_tdp_control, qs_env) 

      TYPE(xas_tdp_env_type), POINTER                 :: xas_tdp_env 
      TYPE(xas_tdp_control_type), POINTER             :: xas_tdp_control
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'write_mos_to_ex_atoms_association', &
         routineP = moduleN//':'//routineN

      INTEGER                                         :: imo, nmo, at_index, tmp_index, &
                                                         output_unit, ispin, nspins
      REAL(dp), DIMENSION(3)                          :: wfn_center, at_pos, r_ac
      INTEGER, DIMENSION(3)                           :: perd_init
      TYPE(cell_type), POINTER                        :: cell
      TYPE(particle_type), DIMENSION(:), POINTER      :: particle_set
      INTEGER, DIMENSION(:,:,:), POINTER              :: mos_of_ex_atoms
      INTEGER, DIMENSION(:), POINTER                  :: ex_atom_indices
      REAL(dp)                                        :: dist, mo_spread
      CHARACTER(LEN=2)                                :: at_symbol

      NULLIFY (cell, particle_set, mos_of_ex_atoms, ex_atom_indices)

      output_unit = cp_logger_get_default_io_unit()

      IF (output_unit > 0) THEN
         WRITE(UNIT=output_unit, FMT="(/,T3,A,/,T3,A,/,T3,A)") &
            "                  Associated    Associated        Distance to   MO spread (Ang^2)",&
            "Spin  MO index    atom index     atom kind    MO center (Ang)   -w_i ln(|z_ij|^2)",&
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
      DO imo = 1,nmo
         DO ispin = 1,nspins

!           each Mo is associated to at most one atom (only 1 in array of -1)
            IF (ANY(mos_of_ex_atoms(imo,:,ispin) == 1)) THEN
               tmp_index = MAXLOC(mos_of_ex_atoms(imo,:,ispin),1) 
               at_index = ex_atom_indices(tmp_index)
               at_symbol = particle_set(at_index)%atomic_kind%element_symbol

               at_pos = particle_set(at_index)%r
               wfn_center = xas_tdp_env%qs_loc_env%localized_wfn_control%centers_set(ispin)%array(1:3,imo)
               r_ac = pbc(at_pos, wfn_center, cell)
               dist = SQRT(r_ac(1)*r_ac(1) + r_ac(2)*r_ac(2) + r_ac(3)*r_ac(3))
!              convert distance from a.u. to Angstrom
               dist = dist*angstrom         
         
               mo_spread = xas_tdp_env%qs_loc_env%localized_wfn_control%centers_set(ispin)%array(4,imo)
               mo_spread = mo_spread*angstrom*angstrom
         
               IF (output_unit > 0 ) THEN
                  WRITE(UNIT=output_unit, FMT="(T3,I4,I10,I14,A14,ES19.3,ES20.3)") &
                     ispin, imo, at_index, at_symbol, dist, mo_spread
               END IF

            END IF
         END DO !ispin
      END DO !imo

      IF (output_unit > 0 ) THEN
         WRITE(UNIT=output_unit, FMT="(T3,A,/)") &
            "---------------------------------------------------------------------------------"
      END IF

!  Go back to initial BCs
      cell%perd = perd_init

   END SUBROUTINE write_mos_to_ex_atoms_association

!> *************************************************************************************************
!> \brief Performs Mulliken population analysis for the MO(s) of a donor_state_type so that user
!>        can verify it is indeed a core state
!> \param donor_state ...
!> \param qs_env ...
!> \note This is a specific case of Mulliken analysis. In general one computes sum_i (SP)_ii, where
!>       i labels the basis function centered on the atom of interest. For a specific MO with index
!>       j, one need to compute sum_{ik} c_{ij} S_{ik} c_{kj}, k = 1,nao
!> *************************************************************************************************
   SUBROUTINE perform_mulliken_on_donor_state(donor_state, qs_env)
      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(qs_environment_type), POINTER              :: qs_env
      
      CHARACTER(len=*), PARAMETER :: routineN = "perfomr_mulliken_on_donor_state", &
                                     routineP = moduleN//":"//routineN

      INTEGER                                         :: at_index, i, nao, natom, ndo_mo, nsgf, &
                                                         output_unit, nspins, ispin, ndo_so
      INTEGER, DIMENSION(:), POINTER                  :: first_sgf, last_sgf
      INTEGER, DIMENSION(:,:), POINTER                :: mo_indices
      TYPE(qs_kind_type), DIMENSION(:), POINTER       :: qs_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER      :: particle_set
      TYPE(cp_fm_type), POINTER                       :: work_vect, gs_coeffs
      REAL(dp), DIMENSION(:,:), POINTER               :: work_array
      REAL(dp), DIMENSION(:,:), ALLOCATABLE           :: pop_mat, mul_pop
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_s
      TYPE(cp_para_env_type), POINTER                 :: para_env
      TYPE(cp_blacs_env_type), POINTER                :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                :: col_vect_struct

      NULLIFY(mo_indices, qs_kind_set, particle_set, first_sgf,work_array)
      NULLIFY(matrix_s, para_env, blacs_env, col_vect_struct, work_vect, last_sgf)

!  Initialization
      at_index = donor_state%at_index
      mo_indices => donor_state%mo_indices
      ndo_mo = donor_state%ndo_mo
      gs_coeffs => donor_state%gs_coeffs
      output_unit = cp_logger_get_default_io_unit()
      nspins = 1; IF (SIZE(mo_indices,2) == 2) nspins = 2
      ndo_so = ndo_mo*nspins
      ALLOCATE(mul_pop(ndo_mo, nspins))
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
      ALLOCATE(work_array(nsgf,ndo_so))
      ALLOCATE(pop_mat(ndo_so,ndo_so))

      CALL cp_fm_get_submatrix(fm=work_vect, target_m=work_array, start_row=first_sgf(at_index), &
                               start_col=1, n_rows=nsgf, n_cols=ndo_so, transpose=.FALSE.)

      CALL dgemm('T', 'N', ndo_so, ndo_so, nsgf, 1.0_dp, donor_state%contract_coeffs, nsgf, &
                 work_array, nsgf, 0.0_dp, pop_mat, ndo_so)

!  The Mullikan population for the MOs in on the diagonal.
      DO ispin = 1,nspins
         DO i = 1,ndo_mo
            mul_pop(i,ispin) = pop_mat((ispin-1)*ndo_mo+i,(ispin-1)*ndo_mo+i)
         END DO
      END DO

!  Printing in main output file
      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(T5,A,/,T5,A)") &
            "Mulliken population analysis retricted to the associated MO(s) yields: ", &
            "                                              Spin  MO index     charge"
         DO ispin = 1, nspins
            DO i = 1,ndo_mo
               WRITE (UNIT=output_unit, FMT="(T51,I4,I10,F11.3)") &
                  ispin, mo_indices(i,ispin), mul_pop(i,ispin)
            END DO
         END DO
      END IF

!  Clean-up
      DEALLOCATE(first_sgf, last_sgf, work_array)
      CALL cp_fm_release(work_vect)

   END SUBROUTINE perform_mulliken_on_donor_state

! **************************************************************************************************
!> \brief Prints the excitation energies and the oscillator strengths for a given donor_state in a file
!> \param donor_state the donor_state to print
! **************************************************************************************************
   SUBROUTINE print_xas_tdp_to_file(donor_state, xas_tdp_env, xas_tdp_control, xas_tdp_section)

      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(section_vals_type), POINTER                   :: xas_tdp_section

      CHARACTER(len=*), PARAMETER :: routineN = 'print_xas_tdp_to_file', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, xas_tdp_unit
      TYPE(cp_logger_type), POINTER                      :: logger

      NULLIFY (logger)
      logger => cp_get_default_logger()

      xas_tdp_unit = cp_print_key_unit_nr(logger, xas_tdp_section, "PRINT%SPECTRUM", &
                                          extension=".spectrum", file_position="APPEND", &
                                          file_action="WRITE", file_form="FORMATTED")

      IF (xas_tdp_unit > 0) THEN

         IF (xas_tdp_control%do_spin_cons) THEN

!           Printing the general donor state information
            WRITE(xas_tdp_unit, FMT="(A,/,A,A,A/,A,I5,A,I5,A,A,/,A)") &
               "================================================================================", &
               "XAS TDP open-shell spin-conserving (no SOC) excitations for DONOR STATE: ", &
                xas_tdp_env%state_type_char(donor_state%state_type), ",",&
               "from EXCITED ATOM: ", donor_state%at_index, ", of KIND (index/symbol): ", &
               donor_state%kind_index, "/", donor_state%at_symbol,&
               "================================================================================"

!           Simply dump the excitation energies/ oscillator strength as they come
            WRITE(xas_tdp_unit, FMT="(T3,A)") &
               " Index            Excitation energy (eV)            Oscillator strength (a.u.)"

            DO i =1,SIZE(donor_state%sc_evals)
               WRITE(xas_tdp_unit, FMT="(T3,I6,F34.6,F38.6)") &
                  i, donor_state%sc_evals(i)*evolt, donor_state%sc_osc_str(i)
            END DO

            WRITE(xas_tdp_unit, FMT="(A,/)") " "
         END IF ! do_spin_cons

         IF (xas_tdp_control%do_spin_flip) THEN

!           Printing the general donor state information
            WRITE(xas_tdp_unit, FMT="(A,/,A,A,A/,A,I5,A,I5,A,A,/,A)") &
               "================================================================================", &
               "XAS TDP open-shell spin-flip (no SOC) excitations for DONOR STATE: ", &
                xas_tdp_env%state_type_char(donor_state%state_type), ",",&
               "from EXCITED ATOM: ", donor_state%at_index, ", of KIND (index/symbol): ", &
               donor_state%kind_index, "/", donor_state%at_symbol,&
               "================================================================================"

!           Simply dump the excitation energies/ oscillator strength as they come
            WRITE(xas_tdp_unit, FMT="(T3,A)") &
               " Index            Excitation energy (eV)            Oscillator strength (a.u.)"

            DO i =1,SIZE(donor_state%sf_evals)
               WRITE(xas_tdp_unit, FMT="(T3,I6,F34.6,F38.6)") &
                  i, donor_state%sf_evals(i)*evolt, 0.0_dp !spin-forbidden !
            END DO

            WRITE(xas_tdp_unit, FMT="(A,/)") " "
         END IF ! do_spin_flip

         IF (xas_tdp_control%do_singlet) THEN

!           Printing the general donor state information
            WRITE(xas_tdp_unit, FMT="(A,/,A,A,A/,A,I5,A,I5,A,A,/,A)") &
               "================================================================================", &
               "XAS TDP singlet excitations (no SOC) for DONOR STATE: ", &
               xas_tdp_env%state_type_char(donor_state%state_type), ",", &
               "from EXCITED ATOM: ", donor_state%at_index, ", of KIND (index/symbol): ", &
               donor_state%kind_index, "/", donor_state%at_symbol,&
               "================================================================================"

!           Simply dump the excitation energies/ oscillator strength as they come
            WRITE(xas_tdp_unit, FMT="(T3,A)") &
               " Index            Excitation energy (eV)            Oscillator strength (a.u.)"

            DO i =1,SIZE(donor_state%sg_evals)
               WRITE(xas_tdp_unit, FMT="(T3,I6,F34.6,F38.6)") &
                  i, donor_state%sg_evals(i)*evolt, donor_state%sg_osc_str(i)
            END DO

            WRITE (xas_tdp_unit, FMT="(A,/)") " "
         END IF ! do_singlet

         IF (xas_tdp_control%do_triplet) THEN

!           Printing the general donor state information
            WRITE(xas_tdp_unit, FMT="(A,/,A,A,A/,A,I5,A,I5,A,A,/,A)") &
               "================================================================================", &
               "XAS TDP triplet excitations (no SOC) for DONOR STATE: ", &
               xas_tdp_env%state_type_char(donor_state%state_type), ",", &
               "from EXCITED ATOM: ", donor_state%at_index, ", of KIND (index/symbol): ", &
               donor_state%kind_index, "/", donor_state%at_symbol,&
               "================================================================================"

!           Simply dump the excitation energies/ oscillator strength as they come
            WRITE(xas_tdp_unit, FMT="(T3,A)") &
               " Index            Excitation energy (eV)            Oscillator strength (a.u.)"

            DO i =1,SIZE(donor_state%tp_evals)
               WRITE(xas_tdp_unit, FMT="(T3,I6,F34.6,F38.6)") &
                  i, donor_state%tp_evals(i)*evolt, 0.0_dp !triplets from closed-shell have none
            END DO

            WRITE (xas_tdp_unit, FMT="(A,/)") " "
         END IF ! do_triplet

         IF (xas_tdp_control%do_soc) THEN

!           Printing the general donor state information
            WRITE(xas_tdp_unit, FMT="(A,/,A,A,A/,A,I5,A,I5,A,A,/,A)") &
               "================================================================================", &
               "XAS TDP  excitations after spin-orbit coupling for DONOR STATE: ", &
               xas_tdp_env%state_type_char(donor_state%state_type), ",", &
               "from EXCITED ATOM: ", donor_state%at_index, ", of KIND (index/symbol): ", &
               donor_state%kind_index, "/", donor_state%at_symbol,&
               "================================================================================"

!           Simply dump the excitation energies/ oscillator strength as they come
            WRITE(xas_tdp_unit, FMT="(T3,A)") &
               " Index            Excitation energy (eV)            Oscillator strength (a.u.)"

            DO i =1,SIZE(donor_state%soc_evals)
               WRITE(xas_tdp_unit, FMT="(T3,I6,F34.6,F38.6)") &
                  i, donor_state%soc_evals(i)*evolt, donor_state%soc_osc_str(i)
            END DO

            WRITE (xas_tdp_unit, FMT="(A,/)") " "

         END IF !do_soc

      END IF

      CALL cp_print_key_finished_output(xas_tdp_unit, logger, xas_tdp_section, "PRINT%SPECTRUM")

   END SUBROUTINE print_xas_tdp_to_file
END MODULE xas_tdp_methods
