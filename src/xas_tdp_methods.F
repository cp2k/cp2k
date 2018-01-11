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
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: &
        cp_fm_create, cp_fm_get_info, cp_fm_get_submatrix, cp_fm_p_type, cp_fm_release, &
        cp_fm_set_all, cp_fm_set_submatrix, cp_fm_to_fm_submat, cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_unit_nr,&
                                              debug_print_level
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_p_type,&
                                              dbcsr_release,&
                                              dbcsr_type
   USE input_constants,                 ONLY: do_loc_none,&
                                              op_loc_berry,&
                                              xas_1s_type,&
                                              xas_2p_type,&
                                              xas_2s_type,&
                                              xas_not_excited,&
                                              xas_tdp_by_index,&
                                              xas_tdp_by_kind
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_set
   USE kinds,                           ONLY: dp
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE periodic_table,                  ONLY: ptable
   USE physcon,                         ONLY: angstrom
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_loc_methods,                  ONLY: centers_spreads_berry,&
                                              qs_loc_driver
   USE qs_loc_types,                    ONLY: get_qs_loc_env,&
                                              localized_wfn_control_type,&
                                              qs_loc_env_create,&
                                              qs_loc_env_new_type
   USE qs_loc_utils,                    ONLY: qs_loc_control_init,&
                                              qs_loc_env_init,&
                                              set_loc_centers,&
                                              set_loc_wfn_lists
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE xas_methods,                     ONLY: calc_stogto_overlap
   USE xas_tdp_types,                   ONLY: &
        donor_state_type, read_xas_tdp_control, set_donor_state, set_xas_tdp_env, &
        xas_tdp_control_create, xas_tdp_control_type, xas_tdp_env_create, xas_tdp_env_release, &
        xas_tdp_env_type
   USE xas_tdp_utils,                   ONLY: compute_oscillator_strength,&
                                              setup_xas_tdp_geeig_prob,&
                                              solve_xas_tdp_geeig_prob
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xas_tdp_methods'

   PUBLIC :: xas_tdp

CONTAINS

! **************************************************************************************************
!> \brief Driver for XAS TDDFT calculations.
!> \param qs_env the inherited qs_environment
!> \param dft_control the set of input DFT keywords and values
!> \author AB
!> \note Empty for now...
! **************************************************************************************************
   SUBROUTINE xas_tdp(qs_env, dft_control)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dft_control_type), POINTER                    :: dft_control

      CHARACTER(len=*), PARAMETER :: routineN = 'xas_tdp', routineP = moduleN//':'//routineN

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
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env

      NULLIFY (xas_tdp_env, xas_tdp_control, atomic_kind_set, atoms_of_kind, current_state)
      NULLIFY (matrix_tdp, matrix_lhs)

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
      CALL xas_tdp_init(xas_tdp_env, xas_tdp_control, qs_env, dft_control)

!     Building an array for better writing of donor state info. Indices correspond input constants
      state_type_char(1) = "1s"
      state_type_char(2) = "2s"
      state_type_char(3) = "2p"

!  Safety checks
      IF (dft_control%lsd) THEN
         CPABORT("XAS TDP only implemented for closed shell ground state calculations")
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

!  Loop over donor states for calculation

      CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set)
      current_state_index = 1

!     Loop over atomic kinds
      DO ikind = 1, SIZE(atomic_kind_set)

         CALL get_atomic_kind(atomic_kind=atomic_kind_set(ikind), element_symbol=symbol_of_kind, &
                              atom_list=atoms_of_kind)

!        Loop over atoms of kind
         DO iat = 1, SIZE(atoms_of_kind)

            IF (ANY(xas_tdp_env%ex_atom_indices == atoms_of_kind(iat))) THEN
               tmp_index = MINLOC(ABS(xas_tdp_env%ex_atom_indices-atoms_of_kind(iat)), 1)

!              Loop over states of excited atom of kind
               DO istate = 1, SIZE(xas_tdp_env%state_types, 1)

                  IF (xas_tdp_env%state_types(istate, tmp_index) /= xas_not_excited) THEN

                     current_state => xas_tdp_env%donor_states(current_state_index)
                     CALL set_donor_state(current_state, &
                                          at_index=xas_tdp_env%ex_atom_indices(tmp_index), &
                                          at_symbol=symbol_of_kind, kind_index=ikind, &
                                          state_type=xas_tdp_env%state_types(istate, tmp_index))

!                    Initial write for the donor state of interest
                     IF (output_unit > 0) THEN
                        WRITE (UNIT=output_unit, FMT="(/,T3,A,A2,A,I4,A,A2,/)") &
                           "# Start of calculations for donor state(s) of type ", &
                           state_type_char(current_state%state_type), " for atom", &
                           current_state%at_index, " of kind ", current_state%at_symbol

                     END IF

!                    Assign best fitting MO(s) to current core donnor state
                     CALL assign_mo_to_donor_state(current_state, xas_tdp_env, xas_tdp_control, &
                                                   qs_env, output_unit)

!                    Perform MO restricted Mulliken pop analysis for verification
                     CALL perform_mulliken_on_donor_state(current_state, qs_env, output_unit)

                     IF (.NOT. xas_tdp_control%check_donor_states) THEN
!                       Do main calculations here
                        PRINT *, "Should perform calculations"

!                       Testing
                        ALLOCATE (matrix_tdp)
                        ALLOCATE (matrix_lhs)
                        CALL setup_xas_tdp_geeig_prob(matrix_tdp, matrix_lhs, current_state, &
                                                      qs_env, xas_tdp_control)
                        CALL solve_xas_tdp_geeig_prob(matrix_tdp, matrix_lhs, current_state, &
                                                      qs_env, xas_tdp_control)
                        CALL compute_oscillator_strength(current_state, qs_env, xas_tdp_control)

                        PRINT *, "evals: "
                        PRINT *, current_state%lr_evals

                        CALL dbcsr_release(matrix_tdp)
                        CALL dbcsr_release(matrix_lhs)
                        DEALLOCATE (matrix_tdp)
                        DEALLOCATE (matrix_lhs)
                     END IF

                     current_state_index = current_state_index+1
                     NULLIFY (current_state)
                  END IF
               END DO ! state type
            END IF ! if excited
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

      CALL timestop(handle)

   END SUBROUTINE xas_tdp

! **************************************************************************************************
!> \brief Overall control and  environment types initialization
!> \param xas_tdp_env the environment type to initialize
!> \param xas_tdp_control the control type to initialize
!> \param qs_env the inherited qs environement type
!> \param dft_control the type in which xas_tdp_control lives
! **************************************************************************************************
   SUBROUTINE xas_tdp_init(xas_tdp_env, xas_tdp_control, qs_env, dft_control)

      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dft_control_type), POINTER                    :: dft_control

      CHARACTER(len=*), PARAMETER :: routineN = 'xas_tdp_init', routineP = moduleN//':'//routineN

      CHARACTER(len=2)                                   :: symbol
      INTEGER                                            :: homo, i, j, k, n_donor_states, n_kinds, &
                                                            nat_of_kind, natom, nex_atoms, &
                                                            nex_kinds, nmatch, nspins
      INTEGER, DIMENSION(2)                              :: n_mo, n_moloc
      INTEGER, DIMENSION(:), POINTER                     :: ind_of_kind
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: at_kind_set
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(qs_loc_env_new_type), POINTER                 :: qs_loc_env
      TYPE(section_vals_type), POINTER                   :: dft_section, loc_section, xas_tdp_section

      NULLIFY (dft_section, xas_tdp_section, at_kind_set, ind_of_kind)
      NULLIFY (qs_loc_env, loc_section, mos)

!  XAS TDP control type initialization
      dft_section => section_vals_get_subs_vals(qs_env%input, "DFT")
      xas_tdp_section => section_vals_get_subs_vals(dft_section, "XAS_TDP")

      CALL xas_tdp_control_create(dft_control%xas_tdp_control)
      CALL read_xas_tdp_control(dft_control%xas_tdp_control, xas_tdp_section)

      xas_tdp_control => dft_control%xas_tdp_control

!  XAS TDP environment type initialization
      CALL xas_tdp_env_create(xas_tdp_env)

!  Retrieving the excited atoms indices and correspondig state types
      IF (xas_tdp_control%define_excited == xas_tdp_by_index) THEN

!        simply copy from xas_tdp_control
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

      ELSE IF (xas_tdp_control%define_excited == xas_tdp_by_kind) THEN

!        need to find out which atom of which kind is excited
         CALL get_qs_env(qs_env=qs_env, atomic_kind_set=at_kind_set)
         n_kinds = SIZE(at_kind_set)
         nex_atoms = 0

         DO i = 1, n_kinds
            CALL get_atomic_kind(atomic_kind=at_kind_set(i), element_symbol=symbol, &
                                 natom=nat_of_kind)
            IF (ANY(xas_tdp_control%list_ex_kinds == symbol)) nex_atoms = nex_atoms+nat_of_kind
         END DO

         ALLOCATE (xas_tdp_env%ex_atom_indices(nex_atoms))
         ALLOCATE (xas_tdp_env%state_types(SIZE(xas_tdp_control%state_types, 1), nex_atoms))
         nex_kinds = SIZE(xas_tdp_control%list_ex_kinds)
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

         CALL set_xas_tdp_env(xas_tdp_env, nex_atoms=nex_atoms)

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
         CALL set_loc_wfn_lists(qs_loc_env%localized_wfn_control, n_moloc, n_mo, nspins)
         CALL set_loc_centers(qs_loc_env%localized_wfn_control, n_moloc, nspins)
!        closed shell => myspin=1
         CALL qs_loc_env_init(qs_loc_env, qs_loc_env%localized_wfn_control, &
                              qs_env, myspin=1, do_localize=qs_loc_env%do_localize)
      END IF

!  Allocating memory for the array of excited atoms MOs. Worst case senario, all searched MOs are
!  associated to the same atom
      ALLOCATE (xas_tdp_env%mos_of_ex_atoms(xas_tdp_control%n_search, nex_atoms))

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
            at_index = MINLOC(ABS(xas_tdp_env%ex_atom_indices-iat_memo), 1)
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

      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: output_unit

      CHARACTER(len=*), PARAMETER :: routineN = 'assign_mo_to_donor_state', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: at_index, i, iat, imo, isgf, j, k, l, &
                                                            my_mo, n_search, n_states, nao, &
                                                            next_best_overlap_ind, nj, nsgf_kind, &
                                                            nsgf_sto, zval
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: my_mos
      INTEGER, DIMENSION(4, 7)                           :: ne
      INTEGER, DIMENSION(:), POINTER                     :: first_sgf, lq, nq
      INTEGER, DIMENSION(:, :), POINTER                  :: mos_of_ex_atoms
      REAL(dp)                                           :: component, next_best_overlap, zeff
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: max_overlap, overlap
      REAL(dp), DIMENSION(:), POINTER                    :: mo_evals, zeta
      REAL(dp), DIMENSION(:, :), POINTER                 :: overlap_matrix, tmp_coeff, work_array
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: col_vect_struct
      TYPE(cp_fm_type), POINTER                          :: col_coeff_j, mo_coeff, work_vector
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks
      TYPE(gto_basis_set_type), POINTER                  :: kind_basis_set, sto_to_gto_basis_set
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(sto_basis_set_type), POINTER                  :: sto_basis_set

      NULLIFY (sto_basis_set, sto_to_gto_basis_set, qs_kind_set, kind_basis_set, lq, nq, zeta)
      NULLIFY (overlap_matrix, mos, mo_coeff, mos_of_ex_atoms, tmp_coeff, first_sgf, particle_set)
      NULLIFY (mo_evals, matrix_ks, col_vect_struct, para_env, blacs_env, col_coeff_j, work_vector)
      NULLIFY (work_array)

      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set, mos=mos, particle_set=particle_set, &
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

!     Electronic configureation (copied from MI's XAS)
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

!  Expanding the STO into (normalized) GTOs for later calculations
      CALL create_gto_from_sto_basis(sto_basis_set=sto_basis_set, &
                                     gto_basis_set=sto_to_gto_basis_set, &
                                     ngauss=xas_tdp_control%n_gauss)
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
      iat = MINLOC(ABS(xas_tdp_env%ex_atom_indices-at_index), 1)
      CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff, nao=nao, eigenvalues=mo_evals)
      ALLOCATE (first_sgf(SIZE(particle_set)))
      CALL get_particle_set(particle_set=particle_set, qs_kind_set=qs_kind_set, first_sgf=first_sgf)
      ALLOCATE (tmp_coeff(nao, 1))

      next_best_overlap = 0.0_dp
      max_overlap(:) = 0.0_dp
      ALLOCATE (overlap(n_search))
      overlap(:) = 0.0_dp
      my_mo = 0
      DO imo = 1, n_search
         IF (mos_of_ex_atoms(imo, iat) > 0) THEN
!           Getting the relevant coefficients for the candidate state
            CALL cp_fm_get_submatrix(fm=mo_coeff, target_m=tmp_coeff, start_row=1, &
                                     start_col=imo, n_rows=nao, n_cols=1, transpose=.FALSE.)

!           Computing the overlap
            DO i = 1, nsgf_sto
               component = 0.0_dp
!              Here we go line by line and take the absolute value each time, i.e. we add the
!              absolute values of the overlap of the MO with each component of the STO basis
!              This allows us to discriminate between state types even if the basis is rotated
               DO j = 1, nsgf_kind
                  isgf = first_sgf(at_index)+j-1
                  component = component+overlap_matrix(i, j)*tmp_coeff(isgf, 1)
               END DO
               overlap(imo) = overlap(imo)+ABS(component)
            END DO

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
      DEALLOCATE (first_sgf, overlap_matrix)

!  Dealing with the result
      IF (ALL(my_mos > 0) .AND. ALL(my_mos <= n_search)) THEN
!        Assigning the MO indices to the donor_state
         ALLOCATE (donor_state%mo_indices(n_states))
         donor_state%mo_indices = my_mos

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
            CALL cp_fm_struct_create(col_vect_struct, nrow_global=nao, ncol_global=1, &
                                     para_env=para_env, context=blacs_env)
            CALL cp_fm_create(col_coeff_j, col_vect_struct)
            CALL cp_fm_create(work_vector, col_vect_struct)
            ALLOCATE (work_array(nao, 1))

            DO j = 1, n_states
               IF (output_unit > 0) THEN
                  WRITE (UNIT=output_unit, FMT="(T53,I8,F18.5)") &
                     my_mos(j), max_overlap(j)
               END IF

!              Retrieving mo_coeff for psi_j
               CALL cp_fm_to_fm_submat(msource=mo_coeff, mtarget=col_coeff_j, nrow=nao, ncol=1, &
                                       s_firstrow=1, s_firstcol=my_mos(j), t_firstrow=1, &
                                       t_firstcol=1)

!              Performing matrix product F|psi_j>
               CALL cp_dbcsr_sm_fm_multiply(matrix_ks(1)%matrix, col_coeff_j, work_vector, ncol=1)
               CALL cp_fm_get_submatrix(work_vector, work_array)

               DO i = 1, n_states
!                 Retrieving mo_coeff for psi_i
                  CALL cp_fm_get_submatrix(fm=mo_coeff, target_m=tmp_coeff, start_row=1, &
                                           start_col=my_mos(i), n_rows=nao, n_cols=1, &
                                           transpose=.FALSE.)

!                 Do the final dot product
                  DO k = 1, nao
                     donor_state%energy_evals(i, j) = donor_state%energy_evals(i, j)+ &
                                                      tmp_coeff(k, 1)*work_array(k, 1)
                  END DO
               END DO
            END DO

            DEALLOCATE (work_array)
         END IF

         IF (output_unit > 0) THEN
            WRITE (UNIT=output_unit, FMT="(T5,A,F7.5,A,I4/)") &
               "The next best overlap is ", next_best_overlap, " for MO with index ", &
               next_best_overlap_ind
         END IF

      ELSE
         CPABORT("A core donor state could not be assigned MO(s). Increasing NSEARCH might help.")
      END IF

!  Clean-up
      CALL cp_fm_struct_release(col_vect_struct)
      CALL cp_fm_release(col_coeff_j)
      CALL cp_fm_release(work_vector)
      DEALLOCATE (tmp_coeff)

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

      INTEGER                                            :: at_index, i, imo, j, nao, natom
      INTEGER, DIMENSION(:), POINTER                     :: first_sgf, last_sgf, mo_indices
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: mul_pop
      REAL(dp), DIMENSION(:, :), POINTER                 :: tmp_coeff, work_array
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: col_vect_struct
      TYPE(cp_fm_type), POINTER                          :: col_coeff, mo_coeff, work_vect
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (mos, mo_indices, qs_kind_set, particle_set, mo_coeff, tmp_coeff, first_sgf, work_array)
      NULLIFY (matrix_s, para_env, blacs_env, col_vect_struct, col_coeff, work_vect, last_sgf)

!  Initialization
      at_index = donor_state%at_index
      mo_indices => donor_state%mo_indices
      ALLOCATE (mul_pop(SIZE(mo_indices, 1)))
      mul_pop = 0.0_dp

      CALL get_qs_env(qs_env=qs_env, mos=mos, particle_set=particle_set, qs_kind_set=qs_kind_set, &
                      para_env=para_env, blacs_env=blacs_env, matrix_s=matrix_s)
      CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff, nao=nao)
      ALLOCATE (tmp_coeff(nao, 1))
      ALLOCATE (work_array(nao, 1))

      natom = SIZE(particle_set, 1)
      ALLOCATE (first_sgf(natom))
      ALLOCATE (last_sgf(natom))

      CALL get_particle_set(particle_set, qs_kind_set, first_sgf=first_sgf, last_sgf=last_sgf)

      CALL cp_fm_struct_create(col_vect_struct, nrow_global=nao, ncol_global=1, para_env=para_env, &
                               context=blacs_env)
      CALL cp_fm_create(col_coeff, col_vect_struct)
      CALL cp_fm_create(work_vect, col_vect_struct)

!  Loop over the MOs assigned to the donor_state
      DO i = 1, SIZE(mo_indices, 1)
         imo = mo_indices(i)
         tmp_coeff = 0.0_dp
         CALL cp_fm_get_submatrix(fm=mo_coeff, target_m=tmp_coeff, start_row=1, start_col=imo, &
                                  n_rows=nao, n_cols=1, transpose=.FALSE.)

         CALL cp_fm_set_submatrix(col_coeff, tmp_coeff)

!        Do the matrix product S*c_j, where c_j is the column vector of the relevant MO coefficients
         CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, col_coeff, work_vect, ncol=1)
         CALL cp_fm_get_submatrix(work_vect, work_array)

!        We restrict oureselves to the basis function centered on the atom of interest
         DO j = first_sgf(at_index), last_sgf(at_index)
            mul_pop(i) = mul_pop(i)+tmp_coeff(j, 1)*work_array(j, 1)
         END DO

!        Factor 2 because doubly occupied MOs
         mul_pop(i) = 2*mul_pop(i)
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
      DEALLOCATE (first_sgf, last_sgf, tmp_coeff, work_array)
      CALL cp_fm_struct_release(col_vect_struct)
      CALL cp_fm_release(col_coeff)
      CALL cp_fm_release(work_vect)

   END SUBROUTINE perform_mulliken_on_donor_state
END MODULE xas_tdp_methods
