!------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Methods for X-Ray absorption spectroscopy (XAS) using TDDFPT
!> \author AB (11.2017)
! **************************************************************************************************

MODULE xas_tdp_methods

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
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply
   USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm
   USE cp_fm_struct,                    ONLY: cp_fm_struct_type, cp_fm_struct_create,&
                                              cp_fm_struct_release
   USE cp_fm_types,                     ONLY: cp_fm_type, cp_fm_get_submatrix, cp_fm_create,&
                                              cp_fm_release, cp_fm_set_all, cp_fm_p_type,&
                                              cp_fm_get_info, cp_fm_to_fm_submat, & 
                                              cp_fm_set_submatrix
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_unit_nr,&
                                              debug_print_level
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_output_handling,              ONLY: cp_print_key_unit_nr, debug_print_level, &
                                              cp_print_key_finished_output
   USE dbcsr_api,                       ONLY: dbcsr_p_type, dbcsr_type, dbcsr_release, dbcsr_print, &
                                              dbcsr_create, dbcsr_multiply, dbcsr_add_on_diag, &
                                              dbcsr_finalize, dbcsr_copy, dbcsr_set, dbcsr_add
   USE input_constants,                 ONLY: xas_tdp_by_index, xas_tdp_by_kind, xas_not_excited,&
                                              do_loc_none, xas_1s_type, xas_2s_type, xas_2p_type,&
                                              op_loc_berry, xas_dip_len, xas_dip_vel, xas_tdp_hfx
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_get_subs_vals,&
                                              section_vals_val_set
   USE kinds,                           ONLY: dp
   USE memory_utilities,                ONLY: reallocate
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE periodic_table,                  ONLY: ptable
   USE physcon,                         ONLY: angstrom, evolt
   USE pw_types,                        ONLY: pw_p_type, pw_release
   USE qs_environment_types,            ONLY: qs_environment_type,&
                                              get_qs_env
   USE qs_kind_types,                   ONLY: get_qs_kind,&   
                                              qs_kind_type,&
                                              get_qs_kind_set
   USE qs_loc_methods,                  ONLY: qs_loc_driver,&
                                              centers_spreads_berry
   USE qs_loc_types,                    ONLY: localized_wfn_control_type,&
                                              qs_loc_env_create,&
                                              qs_loc_env_new_type
   USE qs_loc_utils,                    ONLY: qs_loc_control_init,&
                                              qs_loc_env_init,&
                                              set_loc_centers,&
                                              set_loc_wfn_lists
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE qs_operators_ao,                 ONLY: p_xyz_ao, rRc_xyz_ao       
   USE qs_rho_types,                    ONLY: qs_rho_get, qs_rho_type
   USE util,                            ONLY: locate
   USE xas_methods,                     ONLY: calc_stogto_overlap
   USE xas_tdp_atom,                    ONLY: xas_atom_env_type, init_xas_atom_env, &
                                              release_xas_atom_env, integrate_xas_atoms
   USE xas_tdp_types,                   ONLY: xas_tdp_env_type,&
                                              xas_tdp_env_create,&
                                              xas_tdp_env_release,&
                                              set_xas_tdp_env,&
                                              xas_tdp_control_type, &
                                              read_xas_tdp_control, &
                                              xas_tdp_control_create, & 
                                              donor_state_type,&
                                              set_donor_state, &
                                              free_ds_memory
   USE xas_tdp_utils,                   ONLY: setup_xas_tdp_full_prob, &
                                              solve_xas_tdp_full_prob, &
                                              compute_dipole_fosc, &
                                              compute_ri_coulomb3_int, &
                                              compute_ri_coulomb2_int
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
                                                            current_state_index, output_unit, iatom
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set       
      CHARACTER(len=2)                                   :: symbol_of_kind
      CHARACTER(len=2), DIMENSION(3)                     :: state_type_char
      INTEGER                                            :: current_state_index, handle, iat, ikind, &
                                                            istate, output_unit, tmp_index
      INTEGER, DIMENSION(:), POINTER                     :: atoms_of_kind
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type), POINTER                          :: matrix_lhs, matrix_tdp
      TYPE(donor_state_type), POINTER                    :: current_state
      TYPE(section_vals_type), POINTER                   :: dft_section, xas_tdp_section
      TYPE(cp_para_env_type), POINTER                    :: para_env
      REAL(dp), DIMENSION(3)                             :: rc
      TYPE(cell_type), POINTER                           :: cell
      TYPE(xas_atom_env_type), POINTER                   :: xas_atom_env
      TYPE(dft_control_type), POINTER                    :: dft_control

      NULLIFY(xas_tdp_env, xas_tdp_control, atomic_kind_set, atoms_of_kind, current_state, para_env)
      NULLIFY(particle_set, cell, xas_atom_env, dft_control)

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

!  Safety checks
      IF (dft_control%lsd) THEN
         CPABORT("XAS TDP only implemented for closed shell ground state calculations")
      END IF
!     Check on periodicity: only deal with non-periodic systems for now
      CALL get_qs_env(qs_env, cell=cell)
      IF (.NOT. ALL(cell%perd .EQ. 0)) THEN
         CPABORT("XAS TDP calculation only supports non-periodic BCs (for now)")
      END IF
      IF (xas_tdp_env%qs_loc_env%do_localize .AND. &
          xas_tdp_env%qs_loc_env%localized_wfn_control%operator_type /= op_loc_berry) THEN
         CPABORT("Please use the Berry phase operator for MO localization.")
      END IF

!  Localization of the core orbitals if requested (used for better identification of donor states)
      IF (xas_tdp_env%qs_loc_env%do_localize) THEN
         IF (output_unit > 0) THEN
            WRITE (UNIT=output_unit, FMT="(/,T3,A,/)") &
               "# Localizing core orbitals for better identification"
         END IF
!        closed shell => myspin=1
         CALL qs_loc_driver(qs_env, xas_tdp_env%qs_loc_env, xas_tdp_control%print_loc_subsection, &
                            myspin=1)

!  If no active localization, need to find centers and spreads of MO
      ELSE IF (.NOT. xas_tdp_env%qs_loc_env%do_localize) THEN
         CALL find_mo_centers(xas_tdp_env, xas_tdp_control, qs_env)
      END IF

!  Assign lowest energy orbitals to excited atoms
      CALL assign_mos_to_ex_atoms(xas_tdp_env, xas_tdp_control, qs_env)

      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T3,A,I4,A,/)") &
            "# Assign the relevant subset of the ", xas_tdp_control%n_search, &
            "  lowest energy MOs to excited atoms"
      END IF
      CALL write_mos_to_ex_atoms_association(xas_tdp_env, qs_env, output_unit)

!  If not simply exact exchange, setup a xas_atom_env and compute the xc integrals on the atomic grids
!  Done before looping on atoms as it's all done at once
      IF (.NOT. xas_tdp_control%fxc_func == xas_tdp_hfx) THEN
         CALL init_xas_atom_env(xas_atom_env, xas_tdp_env%ex_atom_indices, &
                                xas_tdp_env%ex_kind_indices, xas_tdp_control%grid_info, qs_env) 

         ALLOCATE(xas_tdp_env%ri_atom_int%matrix)

         CALL integrate_xas_atoms(xas_tdp_env%ri_atom_int, xas_tdp_control%do_singlet, &
                                  xas_atom_env, qs_env)

         CALL release_xas_atom_env(xas_atom_env)
      END IF

!  Loop over donor states for calculation
!  TODO: I suspect that for parallel run with mpirun, all the processors do not go at the same speed
!        causing some pointer to be dealloacted before they are called by the other processors
!        Get bad pointer memory leak

      CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, para_env=para_env)      
      current_state_index = 1

!     Loop over atomic kinds
      DO ikind = 1, SIZE(atomic_kind_set)

         CALL get_atomic_kind(atomic_kind=atomic_kind_set(ikind), element_symbol=symbol_of_kind, &
                              atom_list=atoms_of_kind)

!        Loop over atoms of kind
         DO iat = 1,SIZE(atoms_of_kind)
            iatom = atoms_of_kind(iat)

            IF (.NOT. ANY(xas_tdp_env%ex_atom_indices == iatom)) CYCLE
            tmp_index = locate(xas_tdp_env%ex_atom_indices, iatom)

!           If length representation for the dipole, compute it here in the AO basis, for this atom
            IF (xas_tdp_control%dipole_form == xas_dip_len) THEN
               CALL get_qs_env(qs_env, particle_set=particle_set)
               rc = particle_set(iatom)%r    
               CALL dbcsr_set(xas_tdp_env%dipmat(1)%matrix, 0.0_dp)
               CALL dbcsr_set(xas_tdp_env%dipmat(2)%matrix, 0.0_dp)
               CALL dbcsr_set(xas_tdp_env%dipmat(3)%matrix, 0.0_dp)
               CALL rRc_xyz_ao(xas_tdp_env%dipmat, qs_env, rc, order=1)
            END IF

!           For this excited atom, also compute the RI coulomb integrals
            CALL compute_ri_coulomb3_int(iatom, ikind, xas_tdp_env, qs_env)

!           Loop over states of excited atom of kind
            DO istate = 1,SIZE(xas_tdp_env%state_types,1)

               IF (xas_tdp_env%state_types(istate, tmp_index) == xas_not_excited) CYCLE
                    
               current_state => xas_tdp_env%donor_states(current_state_index)
               CALL set_donor_state(current_state, at_index=iatom,&
                                    at_symbol=symbol_of_kind, kind_index=ikind,&
                                    state_type=xas_tdp_env%state_types(istate,tmp_index))

!              Initial write for the donor state of interest
               IF (output_unit > 0) THEN
                  WRITE(UNIT=output_unit, FMT="(/,T3,A,A2,A,I4,A,A2,/)") &
                     "# Start of calculations for donor state(s) of type ", &
                     xas_tdp_env%state_type_char(current_state%state_type), " for atom",&
                     current_state%at_index, " of kind ", current_state%at_symbol
                           
               END IF

!              Assign best fitting MO(s) to current core donnor state
               CALL assign_mo_to_donor_state(current_state, xas_tdp_env, xas_tdp_control,&
                                             qs_env, output_unit)

!              Perform MO restricted Mulliken pop analysis for verification
               CALL perform_mulliken_on_donor_state(current_state,qs_env,output_unit)

               IF (xas_tdp_control%check_donor_states) CYCLE
!              Do main calculations here

               IF (xas_tdp_control%full_diag) THEN
                  CALL setup_xas_tdp_full_prob(current_state, qs_env, xas_tdp_env, xas_tdp_control)
                  CALL solve_xas_tdp_full_prob(current_state, xas_tdp_control, qs_env)
                  CALL compute_dipole_fosc(current_state, xas_tdp_control, xas_tdp_env)
                           
               ELSE IF (.NOT. xas_tdp_control%full_diag) THEN
                           
                  CPABORT("Iterative solver not implemented yet")

               END IF !full_diag

               CALL print_xas_tdp_to_file(current_state, xas_tdp_env, xas_tdp_section)

                  

!              Free some unneeded attributes of current_state
               CALL free_ds_memory(current_state)

               current_state_index = current_state_index + 1
               NULLIFY(current_state)

            END DO ! state type
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
                                                            nspins, homo, natom, nmatch, kind_ind, &
                                                            at_ind
      INTEGER, DIMENSION(2)                              :: n_mo, n_moloc
      INTEGER, DIMENSION(:), POINTER                     :: ind_of_kind
      CHARACTER(len=2)                                   :: symbol
      TYPE(qs_loc_env_new_type), POINTER                 :: qs_loc_env
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, rho_ao
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(dft_control_type), POINTER                    :: dft_control

      NULLIFY(dft_section, xas_tdp_section, at_kind_set, ind_of_kind, qs_kind_set, dft_control)
      NULLIFY(qs_loc_env, loc_section, mos, particle_set, rho, rho_ao)

!  XAS TDP control type initialization
      dft_section => section_vals_get_subs_vals(qs_env%input, "DFT")
      xas_tdp_section => section_vals_get_subs_vals(dft_section, "XAS_TDP")

      CALL get_qs_env(qs_env, dft_control=dft_control)
      CALL xas_tdp_control_create(dft_control%xas_tdp_control)
      CALL read_xas_tdp_control(dft_control%xas_tdp_control, xas_tdp_section)

      xas_tdp_control => dft_control%xas_tdp_control

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
         ALLOCATE(xas_tdp_env%ex_kind_indices(nex_atoms))
         xas_tdp_env%ex_kind_indices = 0
         k = 0
         CALL get_qs_env(qs_env, particle_set=particle_set)
         DO i =1,nex_atoms
            at_ind = xas_tdp_env%ex_atom_indices(i)
            CALL get_atomic_kind(particle_set(at_ind)%atomic_kind, kind_number=j)
            IF (ALL(ABS(xas_tdp_env%ex_kind_indices-j) .NE. 0)) THEN
               k = k +1
               xas_tdp_env%ex_kind_indices(k) = j
            END IF
         END DO
         nex_kinds = k
         CALL set_xas_tdp_env(xas_tdp_env, nex_kinds=nex_kinds)
         CALL reallocate(xas_tdp_env%ex_kind_indices, 1 , nex_kinds)
         

      ELSE IF (xas_tdp_control%define_excited == xas_tdp_by_kind) THEN

!        need to find out which atom of which kind is excited
         CALL get_qs_env(qs_env=qs_env, atomic_kind_set=at_kind_set)
         n_kinds = SIZE(at_kind_set)
         nex_atoms = 0

         nex_kinds = SIZE(xas_tdp_control%list_ex_kinds)
         ALLOCATE(xas_tdp_env%ex_kind_indices(nex_kinds))
         k = 0

         DO i = 1,n_kinds
            CALL get_atomic_kind(atomic_kind=at_kind_set(i), element_symbol=symbol,&
                                 natom=nat_of_kind, kind_number=kind_ind)
            IF (ANY(xas_tdp_control%list_ex_kinds == symbol)) THEN
               nex_atoms = nex_atoms + nat_of_kind
               k=k+1
               xas_tdp_env%ex_kind_indices(k) = kind_ind
            END IF
         END DO

         ALLOCATE(xas_tdp_env%ex_atom_indices(nex_atoms))
         ALLOCATE(xas_tdp_env%state_types(SIZE(xas_tdp_control%state_types,1),nex_atoms))
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

!  Allocating memory for the array of donor states
      n_donor_states = COUNT(xas_tdp_env%state_types /= xas_not_excited)
      ALLOCATE (xas_tdp_env%donor_states(n_donor_states))

!  Initializing the qs_loc_env from the LOCALIZE subsection of XAS_TDP (largely inpired by MI's XAS)
!  Closed shell ground state calculation is assumed (spin-unpolarized)
      CALL qs_loc_env_create(xas_tdp_env%qs_loc_env)
      qs_loc_env => xas_tdp_env%qs_loc_env
      loc_section => xas_tdp_control%loc_subsection
!     getting the number of (doubly occupied)  MOs
      CALL get_qs_env(qs_env, mos=mos)
      CALL get_mo_set(mos(1)%mo_set, nmo=n_mo(1), homo=homo)
      n_mo(2) = n_mo(1)
      nspins = 1
!     by default, all homo are localized
      IF (xas_tdp_control%n_search < 0) xas_tdp_control%n_search = homo
      CALL qs_loc_control_init(qs_loc_env, loc_section, do_homo=.TRUE., &
                               do_xas=.TRUE., nloc_xas=xas_tdp_control%n_search)
!     final qs_loc_env initialization
      IF (.NOT. qs_loc_env%do_localize) THEN
         qs_loc_env%localized_wfn_control%localization_method = do_loc_none
      ELSE IF (qs_loc_env%do_localize) THEN
         n_moloc = qs_loc_env%localized_wfn_control%nloc_states
         CALL set_loc_wfn_lists(qs_loc_env%localized_wfn_control, n_moloc, n_mo, nspins, my_spin=1)  
         CALL set_loc_centers(qs_loc_env%localized_wfn_control, n_moloc, nspins)  
!        closed shell => myspin=1
         CALL qs_loc_env_init(qs_loc_env, qs_loc_env%localized_wfn_control, &
                              qs_env, myspin=1, do_localize=qs_loc_env%do_localize)
      END IF

!  Allocating memory for the array of excited atoms MOs. Worst case senario, all searched MOs are
!  associated to the same atom
      ALLOCATE(xas_tdp_env%mos_of_ex_atoms(xas_tdp_control%n_search, nex_atoms))

!  Precomputing for RI => compute and invert the two-center two-electron repulsion matrix for the RI
!  aux basis, for each excited kind

      CALL compute_ri_coulomb2_int(xas_tdp_env, qs_env)

!  Compute the projector on the unoccupied, unperturbed ground state: Q = 1 - SP
      CALL get_qs_env(qs_env, rho=rho, matrix_s=matrix_s)
      CALL qs_rho_get(rho, rho_ao=rho_ao)

      ALLOCATE(xas_tdp_env%q_projector)
      CALL dbcsr_create(xas_tdp_env%q_projector, name="Q PROJECTOR", template=matrix_s(1)%matrix, &
                        matrix_type="N") 
!     TODO: Add a filter_eps in the multiplication ? Cuz many almost zero blocks
!     Spin-restricted ground state: occupation number 2 included in density matrix => hence factor 0.5
      CALL dbcsr_multiply('N', 'N', -0.5_dp, matrix_s(1)%matrix, rho_ao(1)%matrix, 0.0_dp, &
                          xas_tdp_env%q_projector)
      CALL dbcsr_add_on_diag(xas_tdp_env%q_projector, 1.0_dp)
      CALL dbcsr_finalize(xas_tdp_env%q_projector)

!  Create the structure for the dipole in the AO basis 
      ALLOCATE(xas_tdp_env%dipmat(3))
      DO i =1,3
         ALLOCATE(xas_tdp_env%dipmat(i)%matrix)          
         CALL dbcsr_copy(xas_tdp_env%dipmat(i)%matrix, matrix_s(1)%matrix, name="XAS TDP dipole matrix")
         CALL dbcsr_set(xas_tdp_env%dipmat(i)%matrix, 0.0_dp)
      END DO

!     Precompute it if velocity represnetation. If length, need to do it for each excited atom 
      IF (xas_tdp_control%dipole_form == xas_dip_vel) CALL p_xyz_ao(xas_tdp_env%dipmat, qs_env)      

   END SUBROUTINE xas_tdp_init

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
!> *************************************************************************************************
   SUBROUTINE assign_mos_to_ex_atoms(xas_tdp_env, xas_tdp_control, qs_env)

      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'assign_mos_to_ex_atoms', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: at_index, iat, iat_memo, imo, n_atoms, &
                                                            n_search, nex_atoms
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
      mos_of_ex_atoms(:, :) = -1
      nex_atoms = xas_tdp_env%nex_atoms
      n_search = xas_tdp_control%n_search
      localized_wfn_control => xas_tdp_env%qs_loc_env%localized_wfn_control
      CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, cell=cell)
      n_atoms = SIZE(particle_set)

!     Temporarly impose periodic BCs in case of Berry's phase operator used for localization
      perd_init = cell%perd
      cell%perd = 1

!  Loop over n_search lowest energy MOs and all atoms
      DO imo = 1, n_search
!        retrieve MO wave function center coordinates. Closed-shell => centers_set(spin=1)
         wfn_center(1:3) = localized_wfn_control%centers_set(1)%array(1:3, imo)
         iat_memo = 0

!        a large enough value to avoid surprises
         dist_min = 10000.0_dp
         DO iat = 1, n_atoms
            at_pos = particle_set(iat)%r
            r_ac = pbc(at_pos, wfn_center, cell)
            dist = SQRT(r_ac(1)*r_ac(1)+r_ac(2)*r_ac(2)+r_ac(3)*r_ac(3))

!           keep memory of which atom is the closest to the wave function center
            IF (dist < dist_min) THEN
               iat_memo = iat
               dist_min = dist
            END IF
         END DO

!        Verify that the closest atom is actually excited and assign the MO if so
         IF (ANY(xas_tdp_env%ex_atom_indices == iat_memo)) THEN
            at_index = locate(xas_tdp_env%ex_atom_indices, iat_memo)
            mos_of_ex_atoms(imo, at_index) = 1
         END IF
      END DO

!  Go back to initial BCs
      cell%perd = perd_init

   END SUBROUTINE assign_mos_to_ex_atoms

!> *************************************************************************************************
!> \brief Assign core MO(s) to a given donor_state, taking the type (1S, 2S, etc) into account.
!>        The projection on a representative Slater-type orbital basis is used as a indicator.
!>        It is assumed that MOs are already assigned to excited atoms based on their center.
!> \param donor_state the donor_state to which a MO must be assigned
!> \param xas_tdp_env ...
!> \param xas_tdp_control ...
!> \param qs_env ...
!> \param output_unit ...
!> *************************************************************************************************
   SUBROUTINE assign_mo_to_donor_state(donor_state, xas_tdp_env, xas_tdp_control, qs_env, &
                                       output_unit)

      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(xas_tdp_env_type), POINTER                 :: xas_tdp_env        
      TYPE(xas_tdp_control_type), POINTER             :: xas_tdp_control
      TYPE(qs_environment_type), POINTER              :: qs_env
      INTEGER, INTENT(IN)                             :: output_unit

      CHARACTER(len=*), PARAMETER :: routineN ='assign_mo_to_donor_state', &
                                                      routineP = moduleN//':'//routineN

      REAL(dp), DIMENSION(:), POINTER                 :: zeta
      INTEGER, DIMENSION(:), POINTER                  :: lq, nq, first_sgf
      TYPE(sto_basis_set_type), POINTER               :: sto_basis_set
      TYPE(gto_basis_set_type), POINTER               :: sto_to_gto_basis_set, kind_basis_set
      INTEGER, DIMENSION(4,7)                         :: ne
      INTEGER                                         :: zval, i, l, nj, nsgf_sto, nsgf_kind, j,&
                                                         n_search, imo, iat, nao, at_index, my_mo,&
                                                         n_states, next_best_overlap_ind
      INTEGER, DIMENSION(:), ALLOCATABLE              :: my_mos
      REAL(dp)                                        :: zeff, next_best_overlap
      TYPE(qs_kind_type), DIMENSION(:), POINTER       :: qs_kind_set
      TYPE(mo_set_p_type), DIMENSION(:), POINTER      :: mos
      TYPE(cp_fm_type), POINTER                       :: mo_coeff, work_mat, eval_mat, gs_coeffs
      INTEGER, DIMENSION(:,:), POINTER                :: mos_of_ex_atoms
      REAL(dp), DIMENSION(:,:), POINTER               :: tmp_coeff, overlap_matrix
      TYPE(particle_type), DIMENSION(:), POINTER      :: particle_set
      REAL(dp), DIMENSION(:), POINTER                 :: mo_evals
      REAL(dp), DIMENSION(:), ALLOCATABLE             :: max_overlap, overlap, sto_overlap
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_ks
      TYPE(cp_fm_struct_type), POINTER                :: eval_mat_struct, gs_struct 
      TYPE(cp_para_env_type), POINTER                 :: para_env                                    
      TYPE(cp_blacs_env_type), POINTER                :: blacs_env 
      
      NULLIFY(sto_basis_set, sto_to_gto_basis_set, qs_kind_set, kind_basis_set, lq, nq, zeta)
      NULLIFY(overlap_matrix, mos, mo_coeff, mos_of_ex_atoms, tmp_coeff, first_sgf, particle_set)
      NULLIFY(mo_evals, matrix_ks, para_env, blacs_env, work_mat)
      NULLIFY(eval_mat_struct, eval_mat, gs_struct, gs_coeffs)
      
      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set, mos=mos, particle_set=particle_set,&
                      matrix_ks=matrix_ks, para_env=para_env, blacs_env=blacs_env)

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
      ALLOCATE (my_mos(n_states))
      ALLOCATE (max_overlap(n_states))

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
      CALL create_gto_from_sto_basis(sto_basis_set=sto_basis_set,&
                                     gto_basis_set=sto_to_gto_basis_set,&
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
      CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff, nao=nao, eigenvalues=mo_evals)
      ALLOCATE (first_sgf(SIZE(particle_set)))
      CALL get_particle_set(particle_set=particle_set, qs_kind_set=qs_kind_set, first_sgf=first_sgf)
      ALLOCATE(tmp_coeff(nsgf_kind,1))
      ALLOCATE(sto_overlap(nsgf_kind))

      next_best_overlap = 0.0_dp
      max_overlap(:) = 0.0_dp
      ALLOCATE (overlap(n_search))
      overlap(:) = 0.0_dp
      my_mo = 0
      DO imo = 1,n_search
         IF(mos_of_ex_atoms(imo,iat) > 0) THEN
            
            sto_overlap = 0.0_dp
            tmp_coeff = 0.0_dp

!           Getting the relevant coefficients for the candidate state
            CALL cp_fm_get_submatrix(fm=mo_coeff, target_m=tmp_coeff, start_row=first_sgf(at_index),&
                                     start_col=imo, n_rows=nsgf_kind, n_cols=1, transpose=.FALSE.)
                                  
!           Computing the product overlap_matrix*coeffs
            CALL dgemm('N', 'N', nsgf_sto, 1, nsgf_kind, 1.0_dp, overlap_matrix, nsgf_sto, &
                       tmp_coeff, nsgf_kind, 0.0_dp, sto_overlap, nsgf_sto)

!           Each element of column vector sto_overlap is the overlap of a basis element of the 
!           generated STO basis with the kind specific orbital basis. Take the sum of the absolute
!           values so that rotation (of the px, py, pz for example) does not hinder our search
            overlap(imo) = SUM(ABS(sto_overlap))

         END IF
      END DO

!  Finding the best overlap(s)
      DO i = 1, n_states
         my_mo = MAXLOC(overlap, 1)
         my_mos(i) = my_mo
         max_overlap(i) = MAXVAL(overlap, 1)
         overlap(my_mo) = 0.0_dp
      END DO
!     Getting the next best overlap (for validation purposes)
      next_best_overlap = MAXVAL(overlap, 1)
      next_best_overlap_ind = MAXLOC(overlap, 1)

!     Some clean-up
      DEALLOCATE(first_sgf, overlap_matrix, tmp_coeff) 

!  Dealing with the result
      IF (ALL(my_mos > 0) .AND. ALL(my_mos <= n_search)) THEN
!        Assigning the MO indices to the donor_state
         ALLOCATE (donor_state%mo_indices(n_states))
         donor_state%mo_indices = my_mos

!        Storing the MOs in the donor_state, as vectors column
         CALL cp_fm_struct_create(gs_struct, nrow_global=nao, ncol_global=n_states, &
                                  para_env=para_env, context=blacs_env) 
         CALL cp_fm_create(donor_state%gs_coeffs, gs_struct)

         DO i = 1, n_states
            CALL cp_fm_to_fm_submat(msource=mo_coeff, mtarget=donor_state%gs_coeffs, nrow=nao, &
                                    ncol=1, s_firstrow=1, s_firstcol=my_mos(i), t_firstrow=1, &           
                                    t_firstcol=i)
         END DO
         gs_coeffs => donor_state%gs_coeffs

!        Assigning corresponding energy eigenvalues and writing some info in standard input file
         ALLOCATE (donor_state%energy_evals(n_states, n_states))
         donor_state%energy_evals = 0.0_dp
         IF (.NOT. xas_tdp_env%qs_loc_env%do_localize) THEN
            IF (output_unit > 0) THEN
               WRITE (UNIT=output_unit, FMT="(T5,A,/,T5,A,/,T5,A)") &
                  "The following canonical MO(s) have been associated with the donor state(s)", &
                  "based on the overlap with the components of a minimal STO basis: ", &
                  "                                                MO index      overlap(sum)"
            END IF

!           Canonical MO, no change in eigenvalues, only diagonal elements
            DO i = 1, n_states
               donor_state%energy_evals(i, i) = mo_evals(my_mos(i))

               IF (output_unit > 0) THEN
                  WRITE (UNIT=output_unit, FMT="(T53,I8,F18.5)") &
                     my_mos(i), max_overlap(i)
               END IF
            END DO

         ELSE IF (xas_tdp_env%qs_loc_env%do_localize) THEN
            IF (output_unit > 0) THEN
               WRITE (UNIT=output_unit, FMT="(T5,A,/,T5,A,/,T5,A)") &
                  "The following localized MO(s) have been associated with the donor state(s)", &
                  "based on the overlap with the components of a minimal STO basis: ", &
                  "                                                MO index      overlap(sum)"
            END IF

!           MO have been rotated, need epsilon_ij = <psi_i|F|psi_j> = sum_{pq} c_{qi}c_{pj} F_{pq}
            CALL cp_fm_create(work_mat, gs_struct)
            CALL cp_fm_struct_create(eval_mat_struct, nrow_global=n_states, ncol_global=n_states, &
                                     para_env=para_env, context=blacs_env)
            CALL cp_fm_create(eval_mat, eval_mat_struct)

!           Loop over the donor states  and print
            DO j = 1,n_states

!              Print info
               IF (output_unit > 0 ) THEN
                  WRITE (UNIT=output_unit, FMT="(T53,I8,F18.5)") &
                     my_mos(j), max_overlap(j) 
               END IF

            END DO

!           Compute gs_coeff^T * matrix_ks * gs_coeff to get the epsilon_ij matrix
            CALL cp_dbcsr_sm_fm_multiply(matrix_ks(1)%matrix, gs_coeffs, work_mat, ncol=n_states)
            CALL cp_fm_gemm('T', 'N', n_states, n_states, nao, 1.0_dp, gs_coeffs, work_mat, 0.0_dp,&
                            eval_mat)

!           Put the epsilon_ij into the donor_state
            CALL cp_fm_get_submatrix(fm=eval_mat,target_m=donor_state%energy_evals)

!           Clean-up
            CALL cp_fm_release(work_mat)
            CALL cp_fm_release(eval_mat)
            CALL cp_fm_struct_release(eval_mat_struct)

         END IF ! do_localize

!        Clean-up
         CALL cp_fm_struct_release(gs_struct) 
         
         IF (output_unit > 0) THEN
            WRITE (UNIT=output_unit, FMT="(T5,A,F7.5,A,I4/)") &
               "The next best overlap is ", next_best_overlap, " for MO with index ", &
               next_best_overlap_ind
         END IF

      ELSE
         CPABORT("A core donor state could not be assigned MO(s). Increasing NSEARCH might help.")
      END IF

   END SUBROUTINE assign_mo_to_donor_state

!> *************************************************************************************************
!> \brief Compute the centers and spreads of (core) MOs using the Berry phase operator
!> \param centers the coordinates of the MO centers and spreads. First 3 elements of each coloum
!>        contain the MO centers x,y,z coordinates, last 2 contain the spread for 2 functionals
!> \param xas_tdp_env ...
!> \param xas_tdp_control ...
!> \param qs_env ...
!> \note xas_tdp_env%qs_loc_env is used and modified. OK since no localization done when this
!>       subroutine is used
!> *************************************************************************************************
   SUBROUTINE find_mo_centers(xas_tdp_env, xas_tdp_control, qs_env)

      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'find_mo_centers', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: dim_op, i, j, n_centers, nao
      REAL(dp), DIMENSION(6)                             :: weights
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: moloc_coeff
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: zij_fm_set
      TYPE(cp_fm_struct_type), POINTER                   :: tmp_fm_struct
      TYPE(cp_fm_type), POINTER                          :: opvec, vectors
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: op_sm_set
      TYPE(localized_wfn_control_type), POINTER          :: loc_wfn_control
      TYPE(qs_loc_env_new_type), POINTER                 :: qs_loc_env
      TYPE(section_vals_type), POINTER                   :: print_loc_section, prog_run_info

      NULLIFY (qs_loc_env, zij_fm_set, cell, print_loc_section, op_sm_set, moloc_coeff, vectors)
      NULLIFY (opvec, loc_wfn_control, tmp_fm_struct, para_env, blacs_env, prog_run_info)

!  Initialization
      qs_loc_env => xas_tdp_env%qs_loc_env
      loc_wfn_control => qs_loc_env%localized_wfn_control
      print_loc_section => xas_tdp_control%print_loc_subsection
      n_centers = xas_tdp_control%n_search
      CALL get_qs_env(qs_env=qs_env, para_env=para_env, blacs_env=blacs_env, cell=cell)

!  Set print option to debug to keep clean output file
      prog_run_info => section_vals_get_subs_vals(print_loc_Section, "PROGRAM_RUN_INFO")
      CALL section_vals_val_set(prog_run_info, keyword_name="_SECTION_PARAMETERS_", &
                                i_val=debug_print_level)

!  Assume qs_loc_env not initialized yet. Impose values we need
      loc_wfn_control%localization_method = do_loc_none
      loc_wfn_control%operator_type = 0 !op_loc_berry
      loc_wfn_control%nloc_states(:) = n_centers
      loc_wfn_control%eps_occ = 0.0_dp
      loc_wfn_control%lu_bound_states(1, :) = 1
      loc_wfn_control%lu_bound_states(2, :) = n_centers
      loc_wfn_control%set_of_states = 2 !state_loc_list
      loc_wfn_control%do_homo = .FALSE.
      ALLOCATE (loc_wfn_control%loc_states(n_centers, 2))
      DO i = 1, n_centers
         loc_wfn_control%loc_states(i, :) = i
      END DO
      CALL set_loc_centers(loc_wfn_control, loc_wfn_control%nloc_states, nspins=1)
      CALL qs_loc_env_init(qs_loc_env, loc_wfn_control, qs_env, myspin=1, do_localize=.TRUE.)

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

!  zij computation, copied from qs_loc_methods:optimize_loc_berry
      DO i = 1, dim_op
         DO j = 1, 2
            CALL cp_fm_create(zij_fm_set(j, i)%matrix, tmp_fm_struct)
            CALL cp_fm_set_all(zij_fm_set(j, i)%matrix, 0.0_dp)
            CALL cp_dbcsr_sm_fm_multiply(op_sm_set(j, i)%matrix, vectors, opvec, ncol=n_centers)
            CALL cp_gemm("T", "N", n_centers, n_centers, nao, 1.0_dp, vectors, opvec, 0.0_dp, &
                         zij_fm_set(j, i)%matrix)
         END DO
      END DO

!  Compute centers (and spread)
      CALL centers_spreads_berry(qs_loc_env=qs_loc_env, zij=zij_fm_set, nmoloc=n_centers, &
                                 cell=cell, weights=weights, ispin=1, &
                                 print_loc_section=print_loc_section, only_initial_out=.TRUE.)

!  Clean-up
      CALL cp_fm_release(opvec)
      CALL cp_fm_struct_release(tmp_fm_struct)
      DO i = 1, dim_op
         DO j = 1, 2
            CALL cp_fm_release(zij_fm_set(j, i)%matrix)
         END DO
      END DO
      DEALLOCATE (zij_fm_set)
      qs_loc_env%do_localize = .FALSE.

   END SUBROUTINE find_mo_centers

!> *************************************************************************************************
!> \brief Writes the core MOs to excited atoms associations in the main output file
!> \param xas_tdp_env ...
!> \param output_unit ...
!> *************************************************************************************************
   SUBROUTINE write_mos_to_ex_atoms_association(xas_tdp_env, qs_env, output_unit)

      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: output_unit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'write_mos_to_ex_atoms_association', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=2)                                   :: at_symbol
      INTEGER                                            :: at_index, imo, nmo, tmp_index
      INTEGER, DIMENSION(3)                              :: perd_init
      INTEGER, DIMENSION(:), POINTER                     :: ex_atom_indices
      INTEGER, DIMENSION(:, :), POINTER                  :: mos_of_ex_atoms
      REAL(dp)                                           :: dist, mo_spread
      REAL(dp), DIMENSION(3)                             :: at_pos, r_ac, wfn_center
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      NULLIFY (cell, particle_set, mos_of_ex_atoms, ex_atom_indices)

      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T3,A,/,T3,A,/,T3,A)") &
            "              Associated    Associated        Distance to   MO spread (Ang^2)", &
            "MO index      atom index     atom kind    MO center (Ang)   -w_i ln(|z_ij|^2)", &
            "-----------------------------------------------------------------------------"
      END IF

!  Initialization
      mos_of_ex_atoms => xas_tdp_env%mos_of_ex_atoms
      ex_atom_indices => xas_tdp_env%ex_atom_indices
      nmo = SIZE(mos_of_ex_atoms, 1)
      CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, cell=cell)

!     because the use of Berry's phase operator implies PBCs
      perd_init = cell%perd
      cell%perd = 1

!  Retrieving all the info for each MO
      DO imo = 1, nmo

!        each Mo is associated to at most one atom (only 1 in array of -1)
         IF (ANY(mos_of_ex_atoms(imo, :) == 1)) THEN
            tmp_index = MAXLOC(mos_of_ex_atoms(imo, :), 1)
            at_index = ex_atom_indices(tmp_index)
            at_symbol = particle_set(at_index)%atomic_kind%element_symbol

            at_pos = particle_set(at_index)%r
            wfn_center = xas_tdp_env%qs_loc_env%localized_wfn_control%centers_set(1)%array(1:3, imo)
            r_ac = pbc(at_pos, wfn_center, cell)
            dist = SQRT(r_ac(1)*r_ac(1)+r_ac(2)*r_ac(2)+r_ac(3)*r_ac(3))
!           convert distance from a.u. to Angstrom
            dist = dist*angstrom

            mo_spread = xas_tdp_env%qs_loc_env%localized_wfn_control%centers_set(1)%array(4, imo)
            mo_spread = mo_spread*angstrom*angstrom

            IF (output_unit > 0) THEN
               WRITE (UNIT=output_unit, FMT="(T3,I8,I16,A14,ES19.3,ES20.3)") &
                  imo, at_index, at_symbol, dist, mo_spread
            END IF

         END IF

      END DO

      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(T3,A,/)") &
            "-----------------------------------------------------------------------------"
      END IF

!  Go back to initial BCs
      cell%perd = perd_init

   END SUBROUTINE write_mos_to_ex_atoms_association

!> *************************************************************************************************
!> \brief Performs Mulliken population analysis for the MO(s) of a donor_state_type so that user
!>        can verify it is indeed a core state
!> \param donor_state ...
!> \param qs_env ...
!> \param output_unit ...
!> \note This is a specific case of Mulliken analysis. In general one computes sum_i (SP)_ii, where
!>       i labels the basis function centered on the atom of interest. For a specific MO with index
!>       j, one need to compute sum_{ik} c_{ij} S_{ik} c_{kj}, k = 1,nao
!> *************************************************************************************************
   SUBROUTINE perform_mulliken_on_donor_state(donor_state, qs_env, output_unit)
      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: output_unit

      CHARACTER(len=*), PARAMETER :: routineN = 'perform_mulliken_on_donor_state', &
         routineP = moduleN//':'//routineN

      INTEGER                                         :: at_index, i, nao, natom, ndo_mo, nsgf
      INTEGER, DIMENSION(:), POINTER                  :: mo_indices, first_sgf, last_sgf
      TYPE(mo_set_p_type), DIMENSION(:), POINTER      :: mos
      TYPE(qs_kind_type), DIMENSION(:), POINTER       :: qs_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER      :: particle_set
      TYPE(cp_fm_type), POINTER                       :: work_vect, gs_coeffs
      REAL(dp), DIMENSION(:,:), POINTER               :: tmp_coeff, work_array
      REAL(dp), DIMENSION(:), ALLOCATABLE             :: mul_pop
      REAL(dp), DIMENSION(:,:), ALLOCATABLE           :: pop_mat
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_s
      TYPE(cp_para_env_type), POINTER                 :: para_env
      TYPE(cp_blacs_env_type), POINTER                :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                :: col_vect_struct

      NULLIFY(mos, mo_indices, qs_kind_set, particle_set, tmp_coeff, first_sgf,work_array)
      NULLIFY(matrix_s, para_env, blacs_env, col_vect_struct, work_vect, last_sgf)

!  Initialization
      at_index = donor_state%at_index
      mo_indices => donor_state%mo_indices
      ndo_mo = SIZE(donor_state%mo_indices,1)
      ALLOCATE(mul_pop(ndo_mo))
      mul_pop = 0.0_dp
      gs_coeffs => donor_state%gs_coeffs

      CALL get_qs_env(qs_env=qs_env, mos=mos, particle_set=particle_set, qs_kind_set=qs_kind_set, &
                      para_env=para_env, blacs_env=blacs_env, matrix_s=matrix_s)
      CALL cp_fm_get_info(gs_coeffs, nrow_global=nao, matrix_struct=col_vect_struct)

      natom = SIZE(particle_set, 1)
      ALLOCATE (first_sgf(natom))
      ALLOCATE (last_sgf(natom))

      CALL get_particle_set(particle_set, qs_kind_set, first_sgf=first_sgf, last_sgf=last_sgf)
      nsgf = last_sgf(at_index)-first_sgf(at_index)+1

      CALL cp_fm_create(work_vect, col_vect_struct)

!  Take the product of S*coeffs
      CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, gs_coeffs, work_vect, ncol=ndo_mo) 

!  Only consider the product coeffs^T * S * coeffs on the atom of interest
      ALLOCATE(tmp_coeff(ndo_mo,nsgf))
      ALLOCATE(work_array(nsgf,ndo_mo))
      ALLOCATE(pop_mat(ndo_mo,ndo_mo))

      CALL cp_fm_get_submatrix(fm=gs_coeffs, target_m=tmp_coeff, start_row=first_sgf(at_index), &
                               start_col=1, n_rows=nsgf, n_cols=ndo_mo, transpose=.TRUE.)

      CALL cp_fm_get_submatrix(fm=work_vect, target_m=work_array, start_row=first_sgf(at_index), &
                               start_col=1, n_rows=nsgf, n_cols=ndo_mo, transpose=.FALSE.)

      CALL dgemm('N', 'N', ndo_mo, ndo_mo, nsgf, 1.0_dp, tmp_coeff, ndo_mo, work_array, nsgf, 0.0_dp, &
                 pop_mat, ndo_mo)

!  The Mullikan population for the MOs in on the diagonal. Factor of 2 because doubly occupied
      DO i = 1,ndo_mo
         mul_pop(i) = 2*pop_mat(i,i)
      END DO

!  Printing in main output file
      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(T5,A,/,T5,A)") &
            "Mulliken population analysis retricted to the associated MO(s) yields: ", &
            "                                                     MO index      charge"
         DO i = 1, SIZE(mo_indices, 1)
            WRITE (UNIT=output_unit, FMT="(T58,I8,F12.3)") &
               mo_indices(i), mul_pop(i)
         END DO
      END IF

!  Clean-up
      DEALLOCATE(first_sgf, last_sgf, tmp_coeff, work_array)
      CALL cp_fm_release(work_vect)

   END SUBROUTINE perform_mulliken_on_donor_state

! **************************************************************************************************
!> \brief Prints the excitation energies and the oscillator strengths for a given donor_state in a file 
!> \param donor_state the donor_state to print
! **************************************************************************************************
   SUBROUTINE print_xas_tdp_to_file(donor_state, xas_tdp_env, xas_tdp_section)

      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(xas_tdp_env_type), POINTER                 :: xas_tdp_env
      TYPE(section_vals_type), POINTER                :: xas_tdp_section   

      CHARACTER(len=*), PARAMETER :: routineN = "print_xas_tdp_to_file", &
                                     routineP = moduleN//":"//routineN

      TYPE(cp_logger_type), POINTER                   :: logger
      INTEGER                                         :: xas_tdp_unit, i

      NULLIFY(logger)
      logger => cp_get_default_logger()

      xas_tdp_unit = cp_print_key_unit_nr(logger, xas_tdp_section, "PRINT%SPECTRUM", &
                                      extension=".spectrum", file_position="APPEND", &
                                      file_action="WRITE", file_form="FORMATTED")

      IF (xas_tdp_unit > 0) THEN

!        Printing the general donor state information
         WRITE(xas_tdp_unit, FMT="(A,/,A,A,A/,A,I5,A,I5,A,A,/,A)") &
            "================================================================", &
            "XAS TDP spectral data for DONOR STATE: ", xas_tdp_env%state_type_char(donor_state%state_type), ",",&
            "from EXCITED ATOM: ", donor_state%at_index, ", of KIND (index/symbol): ", donor_state%kind_index,&
            "/", donor_state%at_symbol,&
            "================================================================"

!        Simply dump the excitation energies/ oscillator strength as they come
!        Carefull analysis is the user's responsibilty 
         WRITE(xas_tdp_unit, FMT="(T3,A)") &
            "Excitation energy (eV)            Oscillator strength (a.u.)"

         DO i =1,SIZE(donor_state%lr_evals)
            WRITE(xas_tdp_unit, FMT="(T3,F22.6,F38.6)") &
               donor_state%lr_evals(i)*evolt, donor_state%osc_str(i)
               !donor_state%lr_evals(i), donor_state%osc_str(i)
         END DO

         WRITE(xas_tdp_unit, FMT="(A,/)") " "
            
      END IF


      CALL cp_print_key_finished_output(xas_tdp_unit, logger, xas_tdp_section, "PRINT%SPECTRUM")

   END SUBROUTINE print_xas_tdp_to_file
END MODULE xas_tdp_methods
