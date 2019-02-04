!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

!> *************************************************************************************************
!> \brief Define XAS TDP control type and associated create, release, etc subroutines, as well as
!>        XAS TDP environment type and associated set, get, etc subroutines
!> \author AB (11.2017)
!> *************************************************************************************************
MODULE xas_tdp_types
   USE cp_array_utils,                  ONLY: cp_2d_r_p_type, cp_1d_r_p_type , cp_1d_i_p_type, &
                                              cp_3d_r_p_type
   USE cp_fm_types,                     ONLY: cp_fm_type, cp_fm_release, cp_fm_p_type
   USE dbcsr_api,                       ONLY: dbcsr_type, dbcsr_release, dbcsr_p_type, &
                                              dbcsr_release_p
   USE input_constants,                 ONLY: xas_dip_vel,&
                                              xas_tdp_by_index,&
                                              xas_tdp_by_kind, &
                                              tddfpt_singlet, &
                                              tddfpt_triplet, &
                                              ot_mini_cg, &
                                              ot_mini_diis
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get,&
                                              section_vals_get_subs_vals
   USE kinds,                           ONLY: dp ,default_string_length
   USE memory_utilities,                ONLY: reallocate
   USE qs_grid_atom,                    ONLY: grid_atom_type, deallocate_grid_atom
   USE qs_harmonics_atom,               ONLY: harmonics_atom_type, deallocate_harmonics_atom
   USE qs_loc_types,                    ONLY: qs_loc_env_new_type,&
                                              qs_loc_env_release
   USE qs_o3c_types,                    ONLY: o3c_container_type, release_o3c_container
   USE qs_ot_types,                     ONLY: qs_ot_settings_type, qs_ot_settings_init
                                     
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! **************************************************************************************************
!> \brief Type containing control information for TDP XAS calculations
!> \param define_excited whether excited atoms are chosen by kind or index
!> \param dipole_form whether the dipole moment is computed in the length or velocity representation
!> \param n_search # of lowest energy MOs to search for donor orbitals
!> \param check_only whether a check run for donor MOs and LUMOs guesses is conducted
!> \param do_hfx whether exact exchange is included
!> \param do_xc wheter xc functional(s) is(are) included
!> \param sx the scaling applied to exact exchange
!> \param tamm_dancoff whether the calculations should be done in the Tamm-Dancoff approximation
!> \param list_ex_atoms list of excited atom indices, kept empty if define_excited=by_kind
!> \param list_ex_kinds list of excited atom kinds, kept empty if define_excited=by_index
!> \param do_loc whether the core MOs should be localized
!> \param do_lsd whether the calculation is spin-unrestricted
!> \param do_roks whether the calculation is restricted open-shell
!> \param do_singlet whether singlet excitations should be computed
!> \param do_triplet whether triplet excitations should be computed
!> \param do_soc whether spin-orbit coupling should be included
!> \param n_excited the number of excited states to compute
!> \param state_types columns correspond to the states to excite for each atom kind/index
!>                    the number of rows is the number of times the keyword is repeated
!> \param grid_info the information about the atomic grids used for the xc kernel integrals
!> \param full_diag whether the full direct diagonalization of the problem should take place
!> \param ot_settings the settings for the OT diagonalization if so required
!> \param max_iter the maximum allowed number of iteration for OT diagonalization
!> \param eps_iter the required accuracy for OT diagonalization
! **************************************************************************************************
   TYPE xas_tdp_control_type
      INTEGER                                 :: define_excited   
      INTEGER                                 :: dipole_form
      INTEGER                                 :: n_search
      INTEGER                                 :: n_excited
      INTEGER                                 :: max_iter
      REAL(dp)                                :: sx
      REAL(dp)                                :: eps_iter
      LOGICAL                                 :: do_hfx
      LOGICAL                                 :: do_xc
      LOGICAL                                 :: check_only
      LOGICAL                                 :: tamm_dancoff
      LOGICAL                                 :: do_loc
      LOGICAL                                 :: do_lsd
      LOGICAL                                 :: do_roks
      LOGICAL                                 :: do_soc
      LOGICAL                                 :: do_singlet
      LOGICAL                                 :: do_triplet
      LOGICAL                                 :: full_diag
      INTEGER, DIMENSION(:), POINTER          :: list_ex_atoms
      CHARACTER(len=default_string_length),&
                     DIMENSION(:), POINTER    :: list_ex_kinds
      INTEGER, DIMENSION(:,:), POINTER        :: state_types 
      TYPE(section_vals_type), POINTER        :: loc_subsection
      TYPE(section_vals_type), POINTER        :: print_loc_subsection
      CHARACTER(len=default_string_length),&
                     DIMENSION(:,:), POINTER  :: grid_info
      TYPE(qs_ot_settings_type), POINTER      :: ot_settings

   END TYPE xas_tdp_control_type
   
!> *************************************************************************************************
!> \brief Type containing informations such as inputs and results for TDP XAS calculations
!> \param state_type_char an array containing the general donor state types as char (1s, 2s, 2p, ...)
!> \param nex_atoms number of excited atoms
!> \param nex_kinds number of excited kinds
!> \param ex_atom_indices array containing the indices of the excited atoms
!> \param ex_kind_indices array containing the indices of the excited kinds
!> \param state_types columns correspond to the different donor states of each excited atom
!> \param qs_loc_env the environemnt type dealing with the possible localization of donor orbitals
!> \param mos_of_ex_atoms links lowest energy MOs to excited atoms. Elements of value 1 mark the
!>        association between the MO irow and the excited atom icolumn. The third index is for spin
!> \param ri_inv_coul the inverse coulomb RI integral (P|Q)^-1, updated for each excited kind
!>        excited kinds, based on basis functions of the RI auxiliary basis used for Coulomb and
!>        exchange kernels. 
!> \param q_projector the projector on the unperturbed, unoccupied ground state as a dbcsr matrix, 
!>        fpr each spin
!> \param dipmat the dbcsr matrices containing the dipole in x,y,z directions evaluated on the 
!>        contracted spherical gaussians. It can either be in the length or the velocity 
!>        representation. For length representation, it has to be computed once with the origin on 
!>        each excited atom
!> \param roks_ks_mat the KS combined ROKS matrix, if needed
!> \param ri_o3c_coul the container for the RI 3-center Coulomb integrals (updated for each ex atom)
!> \param ri_o3c_ex the container for the RI 3-center exchange integrals (updated for each ex atom)
!> \param lb_o3c_coul the lower bond of a,b present in the 3-center RI Coulomb integrals (ab|P)
!> \param ub_o3c_coul the upper bond of a,b present in the 3-center RI Coulomb integrals (ab|P)
!> \param ri_fxc the array of xc integrals of type (P|fxc|Q), for alpha-alpha, alpha-beta and beta-beta
!> \param orb_soc the matrix where the SOC is evaluated wrt the orbital basis set, for x,y,z
!> \param lumo_guess the coefficients for the LUMOs as guees for the OT eigensolver
!> \param lumo_evals the eigenvalues of the guess LUMOs
!> *************************************************************************************************
   TYPE xas_tdp_env_type
      CHARACTER(len=2), DIMENSION(3)          :: state_type_char
      INTEGER                                 :: nex_atoms
      INTEGER                                 :: nex_kinds
      INTEGER, DIMENSION(:), POINTER          :: ex_atom_indices
      INTEGER, DIMENSION(:), POINTER          :: ex_kind_indices
      INTEGER, DIMENSION(:,:), POINTER        :: state_types
      TYPE(o3c_container_type), POINTER       :: ri_o3c_coul
      TYPE(o3c_container_type), POINTER       :: ri_o3c_ex
      INTEGER                                 :: lb_o3c_coul
      INTEGER                                 :: ub_o3c_coul
      TYPE(cp_fm_p_type), DIMENSION(:), &
                                    POINTER   :: lumo_guess
      TYPE(cp_1d_r_p_type), DIMENSION(:), &
                                    POINTER   :: lumo_evals
      TYPE(donor_state_type), DIMENSION(:),&
                                    POINTER   :: donor_states
      INTEGER, DIMENSION(:,:,:), POINTER      :: mos_of_ex_atoms
      TYPE(qs_loc_env_new_type), POINTER      :: qs_loc_env
      REAL(dp), DIMENSION(:,:), POINTER       :: ri_inv_coul
      TYPE(dbcsr_type), POINTER               :: roks_ks_mat
      TYPE(dbcsr_p_type), DIMENSION(:), &
                                    POINTER   :: q_projector
      TYPE(dbcsr_p_type), DIMENSION(:), &
                                    POINTER   :: dipmat
      TYPE(cp_2d_r_p_type), DIMENSION(:,:), &
                                    POINTER   :: ri_fxc
      TYPE(dbcsr_p_type), DIMENSION(:), &
                                    POINTER   :: orb_soc
   END TYPE xas_tdp_env_type

!> *************************************************************************************************
!> \brief Type containing informations about a single donor state
!> \param at_index the index of the atom to which the state belongs
!> \param kind_index the index of the atomic kind to which the state belongs
!> \param ndo_mo the number of donor MOs per spin
!> \param at_symbol the chemical symbol of the atom to which the state belongs
!> \param state_type whether this is a 1s, 2s, etc state
!> \param energy_evals the energy eigenvalue of the donor state, for each spin
!> \param mo_indices indices of associated MOs. Greater than 1 when not a s-type state.
!> \param os_coeffs solutions of the linear-response TDDFT equation for open-shell calculations
!> \param sg_coeffs solutions of the linear-response TDDFT singlet equations
!> \param tp_coeffs solutions of the linear-response TDDFT triplet equations
!> \param gs_coeffs the ground state MO coefficients 
!> \param os_evals open-shell excitation energies
!> \param sg_evals singlet excitation energies => the eigenvalues of the linear response equation
!> \param tp_evals triplet excitation energies => the eigenvalues of the linear response equation
!> \param soc_evals excitation energies after inclusion of SOC
!> \param os_osc_str oscilaltor strenghts for open-shell excitations
!> \param sg_osc_str oscillator strengths for the different excitation energies of singlets
!> \param tp_osc_str oscillator strengths for the different excitation energies of triplets
!> \param soc_osc_str oscillator strengths after the inclusion of SOC
!> \param os_matrix_tdp the dbcsr matrix to be diagonalized for open-shell calculations
!> \param sg_matrix_tdp the dbcsr matrix to be diagonalized to solve the problem for singlets
!> \param tp_matrix_tdp the dbcsr matrix to be diagonalized to solve the problem for triplets
!> \param metric the metric of the linear response problem M*c = omega*S*c and its inverse
!> \param matrix_aux the auxiliary matrix (A-D+E)^1/2 used to make the problem Hermitian
!> *************************************************************************************************
   TYPE donor_state_type
      INTEGER                                 :: at_index
      INTEGER                                 :: kind_index
      INTEGER                                 :: ndo_mo
      CHARACTER(len=2)                        :: at_symbol
      INTEGER                                 :: state_type
      REAL(dp), DIMENSION(:,:), POINTER       :: energy_evals
      INTEGER, DIMENSION(:,:), POINTER        :: mo_indices
      TYPE(cp_fm_p_type), DIMENSION(:), &
                                    POINTER   :: os_coeffs
      TYPE(cp_fm_p_type), DIMENSION(:), &
                                    POINTER   :: sg_coeffs
      TYPE(cp_fm_p_type), DIMENSION(:), &
                                    POINTER   :: tp_coeffs
      TYPE(cp_fm_type), POINTER               :: gs_coeffs
      REAL(dp), DIMENSION(:), POINTER         :: os_evals
      REAL(dp), DIMENSION(:), POINTER         :: sg_evals
      REAL(dp), DIMENSION(:), POINTER         :: tp_evals
      REAL(dp), DIMENSION(:), POINTER         :: soc_evals
      REAL(dp), DIMENSION(:), POINTER         :: os_osc_str
      REAL(dp), DIMENSION(:), POINTER         :: sg_osc_str
      REAL(dp), DIMENSION(:), POINTER         :: tp_osc_str
      REAL(dp), DIMENSION(:), POINTER         :: soc_osc_str
      TYPE(dbcsr_type), POINTER               :: os_matrix_tdp
      TYPE(dbcsr_type), POINTER               :: sg_matrix_tdp
      TYPE(dbcsr_type), POINTER               :: tp_matrix_tdp
      TYPE(dbcsr_p_type), DIMENSION(:), &
                                    POINTER   :: metric
      TYPE(dbcsr_type), POINTER               :: matrix_aux

   END TYPE donor_state_type

!  Some helper types
   TYPE grid_atom_p_type
      TYPE(grid_atom_type), POINTER                   :: grid_atom
   END TYPE grid_atom_p_type
   
   TYPE harmonics_atom_p_type
      TYPE(harmonics_atom_type), POINTER              :: harmonics_atom
   END TYPE harmonics_atom_p_type

! **************************************************************************************************
!> \brief a environment type that contains all the info needed for XAS_TDP atomic grid calculations
!> \param excited_atoms the atoms for which RI xc-kernel calculations must be done 
!> \param excited_kinds the kinds for which RI xc-kernel calculations must be done 
!> \param proc_of_exat the mapping of the excited atoms on the processors
!> \param grid_atom_set the set of atomic grid for each kind
!> \param ri_dcoeff the expansion coefficients to express the spin-density in the RI basis for each atom
!> \param exat_neighbors the neighbors of each excited atom
!> \param kind_3c_overlap contains the (abP) overlap for each kind, where a,b and P share the same center
!> \param ri_sphi_so contains the coefficient for direct contraction from so to sgf, for the ri basis
!> \param orb_sphi_so contains the coefficient for direct contraction from so to sgf, for the orb basis
!> *************************************************************************************************
   TYPE xas_atom_env_type
      INTEGER, DIMENSION(:), POINTER                  :: excited_atoms
      INTEGER, DIMENSION(:), POINTER                  :: excited_kinds
      INTEGER, DIMENSION(:), POINTER                  :: proc_of_exat
      TYPE(grid_atom_p_type), DIMENSION(:), POINTER   :: grid_atom_set
      TYPE(harmonics_atom_p_type), DIMENSION(:), & 
                                             POINTER  :: harmonics_atom_set
      TYPE(cp_1d_r_p_type), DIMENSION(:,:), POINTER   :: ri_dcoeff
      TYPE(cp_3d_r_p_type), DIMENSION(:), POINTER     :: kind_3c_overlap
      TYPE(cp_1d_i_p_type), DIMENSION(:), POINTER     :: exat_neighbors
      TYPE(cp_2d_r_p_type), DIMENSION(:), POINTER     :: ri_sphi_so
      TYPE(cp_2d_r_p_type), DIMENSION(:), POINTER     :: orb_sphi_so
   END TYPE xas_atom_env_type

   CHARACTER(len=*), PARAMETER, PRIVATE  :: moduleN = 'xas_tdp_types'

! *** Public data types ***
   PUBLIC :: xas_tdp_env_type, donor_state_type, xas_tdp_control_type, xas_atom_env_type

! *** Public subroutines ***
   PUBLIC :: set_donor_state, free_ds_memory, &
             xas_tdp_env_create, xas_tdp_env_release, set_xas_tdp_env, &
             xas_tdp_control_create, xas_tdp_control_release, read_xas_tdp_control, &
             xas_atom_env_create, xas_atom_env_release, donor_state_create, free_exat_memory

CONTAINS

! **************************************************************************************************
!> \brief Creates and initializes the xas_tdp_control_type
!> \param xas_tdp_control the type to initialize
! **************************************************************************************************
   SUBROUTINE xas_tdp_control_create(xas_tdp_control)
      
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control

      CHARACTER(len=*), PARAMETER                        :: routineN = 'xas_tdp_control_create', &
                                                            routineP = moduleN//':'//routineN
      
      CPASSERT(.NOT. ASSOCIATED(xas_tdp_control))
      ALLOCATE (xas_tdp_control)
      
      xas_tdp_control%dipole_form = xas_dip_vel
      xas_tdp_control%define_excited = xas_tdp_by_index
      xas_tdp_control%n_search = -1
      xas_tdp_control%do_hfx = .TRUE.
      xas_tdp_control%do_xc = .FALSE.
      xas_tdp_control%sx = 1.0_dp
      xas_tdp_control%check_only = .FALSE.
      xas_tdp_control%tamm_dancoff = .FALSE.
      xas_tdp_control%do_loc = .FALSE.
      xas_tdp_control%do_lsd = .FALSE.
      xas_tdp_control%do_roks = .FALSE.
      xas_tdp_control%do_soc = .FALSE.
      xas_tdp_control%do_singlet = .TRUE.
      xas_tdp_control%do_triplet = .FALSE.
      xas_tdp_control%full_diag = .TRUE.
      xas_tdp_control%n_excited = -1
      xas_tdp_control%eps_iter = 1.0E-6_dp
      xas_tdp_control%max_iter = 1000
      NULLIFY (xas_tdp_control%state_types)
      NULLIFY (xas_tdp_control%list_ex_atoms)
      NULLIFY (xas_tdp_control%list_ex_kinds)
      NULLIFY (xas_tdp_control%loc_subsection)
      NULLIFY (xas_tdp_control%print_loc_subsection)
      NULLIFY (xas_tdp_control%grid_info)
      NULLIFY (xas_tdp_control%ot_settings)

   END SUBROUTINE xas_tdp_control_create

! **************************************************************************************************
!> \brief Releases the xas_tdp_control_type
!> \param xas_tdp_control the type to release
! **************************************************************************************************
   SUBROUTINE xas_tdp_control_release(xas_tdp_control)
      
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control

      CHARACTER(len=*), PARAMETER                        :: routineN = 'xas_tdp_control_release', &
                                                            routineP = moduleN//':'//routineN
      
      IF (ASSOCIATED(xas_tdp_control)) THEN
         IF (ASSOCIATED(xas_tdp_control%list_ex_atoms)) THEN
            DEALLOCATE (xas_tdp_control%list_ex_atoms)
         END IF
         IF (ASSOCIATED(xas_tdp_control%list_ex_kinds)) THEN
            DEALLOCATE (xas_tdp_control%list_ex_kinds)
         END IF
         IF (ASSOCIATED(xas_tdp_control%state_types)) THEN
            DEALLOCATE(xas_tdp_control%state_types)
         END IF
         IF (ASSOCIATED(xas_tdp_control%grid_info)) THEN
            DEALLOCATE(xas_tdp_control%grid_info)
         END IF
         IF (ASSOCIATED(xas_tdp_control%ot_settings)) THEN
            DEALLOCATE(xas_tdp_control%ot_settings)
         END IF
         DEALLOCATE (xas_tdp_control)
      END IF

   END SUBROUTINE xas_tdp_control_release

! **************************************************************************************************
!> \brief Reads the inputs and stores in xas_tdp_control_type
!> \param xas_tdp_control the type where inputs are stored
!> \param xas_tdp_section the section from which input are read
! **************************************************************************************************
   SUBROUTINE read_xas_tdp_control(xas_tdp_control, xas_tdp_section)

      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(section_vals_type), POINTER                   :: xas_tdp_section
      
      CHARACTER(len=*), PARAMETER                        :: routineN = 'read_xas_tdp_control', &
                                                            routineP = moduleN//':'//routineN

      INTEGER                                            :: irep, nrep, nexc, excitation, ot_method
      CHARACTER(len=default_string_length),&
                     DIMENSION(:), POINTER               :: k_list
      INTEGER, DIMENSION(:), POINTER                     :: a_list, t_list

      NULLIFY (k_list, a_list, t_list)
      
      CALL section_vals_val_get(xas_tdp_section, "DEFINE_EXCITED", &
                                i_val=xas_tdp_control%define_excited)

      CALL section_vals_val_get(xas_tdp_section, "DIPOLE_FORM", i_val=xas_tdp_control%dipole_form)

      CALL section_vals_val_get(xas_tdp_section, "CHECK_ONLY", &
                                l_val=xas_tdp_control%check_only)

      CALL section_vals_val_get(xas_tdp_section, "TAMM_DANCOFF", &
                                l_val=xas_tdp_control%tamm_dancoff)
                           
      CALL section_vals_val_get(xas_tdp_section, "SPIN_ORBIT_COUPLING", &
                                l_val=xas_tdp_control%do_soc)

      CALL section_vals_val_get(xas_tdp_section, "EXCITATIONS", i_val=excitation)
      IF (excitation == tddfpt_singlet) THEN
         xas_tdp_control%do_singlet = .TRUE.
         xas_tdp_control%do_triplet = .FALSE.
      ELSE IF (excitation == tddfpt_triplet) THEN
         xas_tdp_control%do_singlet = .FALSE.
         xas_tdp_control%do_triplet = .TRUE.
      ELSE
         xas_tdp_control%do_singlet = .TRUE.
         xas_tdp_control%do_triplet = .TRUE.
      END IF

      IF (.NOT. ASSOCIATED(xas_tdp_control%list_ex_kinds)) THEN
         IF (xas_tdp_control%define_excited .EQ. xas_tdp_by_index) THEN

            ALLOCATE (xas_tdp_control%list_ex_kinds(0))

         ELSE IF (xas_tdp_control%define_excited .EQ. xas_tdp_by_kind) THEN                   

            CALL section_vals_val_get(xas_tdp_section, "KIND_LIST", c_vals=k_list)
         
            IF (ASSOCIATED(k_list)) THEN
               nexc = SIZE(k_list)
               ALLOCATE (xas_tdp_control%list_ex_kinds(nexc))
               xas_tdp_control%list_ex_kinds = k_list
            END IF

         END IF
      END IF

      IF (.NOT. ASSOCIATED(xas_tdp_control%list_ex_atoms)) THEN
         IF (xas_tdp_control%define_excited .EQ. xas_tdp_by_kind) THEN
            
            ALLOCATE (xas_tdp_control%list_ex_atoms(0))

         ELSE IF (xas_tdp_control%define_excited .EQ. xas_tdp_by_index) THEN

            CALL section_vals_val_get(xas_tdp_section, "ATOM_LIST", i_vals=a_list)

            IF (ASSOCIATED(a_list)) THEN
               nexc = SIZE(a_list)
               CALL reallocate(xas_tdp_control%list_ex_atoms, 1, nexc)
               xas_tdp_control%list_ex_atoms = a_list
            END IF

         END IF
      END IF

      CALL section_vals_val_get(xas_tdp_section, "STATE_TYPES", n_rep_val=nrep)

      IF(.NOT. ASSOCIATED(xas_tdp_control%state_types)) THEN
         ALLOCATE (xas_tdp_control%state_types(nrep,nexc))
         DO irep = 1,nrep
            CALL section_vals_val_get(xas_tdp_section, "STATE_TYPES", i_rep_val=irep, i_vals=t_list)
            IF (SIZE(t_list) .NE. nexc) THEN
               CPABORT("The STATE_TYPES keywords do not have the correct number of entries.")
            END IF
            xas_tdp_control%state_types(irep,:) = t_list
         END DO
      END IF
      
      CALL section_vals_val_get(xas_tdp_section, "N_SEARCH", i_val=xas_tdp_control%n_search)

      CALL section_vals_val_get(xas_tdp_section, "GRID", n_rep_val=nrep)

      IF (.NOT. ASSOCIATED(xas_tdp_control%grid_info)) THEN
         ALLOCATE(xas_tdp_control%grid_info(nrep,3))
         DO irep = 1, nrep
            CALL section_vals_val_get(xas_tdp_section, "GRID", i_rep_val=irep, c_vals=k_list)
            IF (SIZE(k_list) .NE. 3) CPABORT("The GRID keyword needs three values")
            xas_tdp_control%grid_info(irep,:) = k_list
         END DO
      END IF

      xas_tdp_control%loc_subsection => section_vals_get_subs_vals(xas_tdp_section, "LOCALIZE")       
      xas_tdp_control%print_loc_subsection => section_vals_get_subs_vals(&
                                                         xas_tdp_control%loc_subsection, "PRINT")                      
      CALL section_vals_val_get(xas_tdp_control%loc_subsection, "_SECTION_PARAMETERS_", &
                                l_val=xas_tdp_control%do_loc)

      !Kernel analysis
      CALL section_vals_val_get(xas_tdp_section, "KERNEL%XC_FUNCTIONAL%_SECTION_PARAMETERS_", &
                                l_val=xas_tdp_control%do_xc) 
      CALL section_vals_val_get(xas_tdp_section, "KERNEL%EXACT_EXCHANGE%_SECTION_PARAMETERS_", &
                                l_val=xas_tdp_control%do_hfx)

      IF (xas_tdp_control%do_hfx) THEN
         CALL section_vals_val_get(xas_tdp_section, "KERNEL%EXACT_EXCHANGE%SCALE", &
                                   r_val=xas_tdp_control%sx)
      END IF

      IF ((.NOT. xas_tdp_control%do_xc) .AND. (.NOT. xas_tdp_control%do_hfx)) THEN
         CPABORT("You must at least define a XC functional OR exact exchange.")
      END IF

      !Diagonalization analysis
      CALL section_vals_val_get(xas_tdp_section, "DIAGONALIZATION%N_EXCITED", &
                                i_val=xas_tdp_control%n_excited)
      CALL section_vals_val_get(xas_tdp_section, "DIAGONALIZATION%FULL_DIAG", &
                                l_val=xas_tdp_control%full_diag)

      ALLOCATE(xas_tdp_control%ot_settings)
      CALL qs_ot_settings_init(xas_tdp_control%ot_settings)
      CALL section_vals_val_get(xas_tdp_section, "DIAGONALIZATION%OT_DIAG%MINIMIZER", i_val=ot_method)
      SELECT CASE (ot_method)
      CASE (ot_mini_cg)
         xas_tdp_control%ot_settings%ot_method = "CG"
      CASE (ot_mini_diis)
         xas_tdp_control%ot_settings%ot_method = "DIIS"
      END SELECT

      CALL section_vals_val_get(xas_tdp_section, "DIAGONALIZATION%OT_DIAG%MAX_ITER", &
                                i_val=xas_tdp_control%max_iter)
      CALL section_vals_val_get(xas_tdp_section, "DIAGONALIZATION%OT_DIAG%EPS_ITER", &
                                r_val=xas_tdp_control%eps_iter)
      
   END SUBROUTINE read_xas_tdp_control

!> *************************************************************************************************
!> \brief Creates a TDP XAS environment type
!> \param xas_tdp_env the type to create
!> *************************************************************************************************
   SUBROUTINE xas_tdp_env_create(xas_tdp_env)

      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env

      CHARACTER(len=*), PARAMETER :: routineN = 'xas_tdp_env_create', &
         routineP = moduleN//':'//routineN

      ALLOCATE (xas_tdp_env)

      xas_tdp_env%nex_atoms = 1
      xas_tdp_env%nex_kinds = 1

      NULLIFY (xas_tdp_env%ex_atom_indices)
      NULLIFY (xas_tdp_env%ex_kind_indices)
      NULLIFY (xas_tdp_env%state_types)
      NULLIFY (xas_tdp_env%donor_states)
      NULLIFY (xas_tdp_env%qs_loc_env)
      NULLIFY (xas_tdp_env%mos_of_ex_atoms)
      NULLIFY (xas_tdp_env%ri_inv_coul)
      NULLIFY (xas_tdp_env%lumo_guess)
      NULLIFY (xas_tdp_env%lumo_evals)
      NULLIFY (xas_tdp_env%roks_ks_mat)
      NULLIFY (xas_tdp_env%q_projector)
      NULLIFY (xas_tdp_env%dipmat)
      NULLIFY (xas_tdp_env%ri_o3c_coul)
      NULLIFY (xas_tdp_env%ri_o3c_ex)
      NULLIFY (xas_tdp_env%ri_fxc)
      NULLIFY (xas_tdp_env%orb_soc)

!     Putting the state types as char manually
      xas_tdp_env%state_type_char(1) = "1s"
      xas_tdp_env%state_type_char(2) = "2s"
      xas_tdp_env%state_type_char(3) = "2p"

   END SUBROUTINE xas_tdp_env_create

!> *************************************************************************************************
!> \brief Releases the TDP XAS environment type
!> \param xas_tdp_env the type to release
!> *************************************************************************************************
   SUBROUTINE xas_tdp_env_release(xas_tdp_env)

      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env

      CHARACTER(len=*), PARAMETER :: routineN = 'xas_tdp_env_release', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, j

      IF (ASSOCIATED(xas_tdp_env)) THEN
         IF (ASSOCIATED(xas_tdp_env%ex_atom_indices)) THEN
            DEALLOCATE (xas_tdp_env%ex_atom_indices)
         END IF
         IF (ASSOCIATED(xas_tdp_env%ex_kind_indices)) THEN
            DEALLOCATE(xas_tdp_env%ex_kind_indices)
         END IF

         IF (ASSOCIATED(xas_tdp_env%state_types)) THEN
            DEALLOCATE (xas_tdp_env%state_types)
         END IF
         IF (ASSOCIATED(xas_tdp_env%donor_states)) THEN
            CALL deallocate_donor_state_set(xas_tdp_env%donor_states)
         END IF
         IF (ASSOCIATED(xas_tdp_env%qs_loc_env)) THEN
            CALL qs_loc_env_release(xas_tdp_env%qs_loc_env)
         END IF
         IF (ASSOCIATED(xas_tdp_env%mos_of_ex_atoms)) THEN
            DEALLOCATE(xas_tdp_env%mos_of_ex_atoms)
         END IF
         IF (ASSOCIATED(xas_tdp_env%ri_inv_coul)) THEN
            DEALLOCATE(xas_tdp_env%ri_inv_coul)
         END IF
         IF (ASSOCIATED(xas_tdp_env%lumo_guess)) THEN
            DO i = 1,SIZE(xas_tdp_env%lumo_guess)
               CALL cp_fm_release(xas_tdp_env%lumo_guess(i)%matrix)
            END DO
            DEALLOCATE(xas_tdp_env%lumo_guess)
         END IF
         IF (ASSOCIATED(xas_tdp_env%lumo_evals)) THEN
            DO i = 1,SIZE(xas_tdp_env%lumo_evals)
               DEALLOCATE(xas_tdp_env%lumo_evals(i)%array)
            END DO
            DEALLOCATE(xas_tdp_env%lumo_evals)
         END IF
         IF (ASSOCIATED(xas_tdp_env%roks_ks_mat)) THEN
            CALL dbcsr_release_p(xas_tdp_env%roks_ks_mat)
         END IF
         IF (ASSOCIATED(xas_tdp_env%q_projector)) THEN
            DO i=1,SIZE(xas_tdp_env%q_projector)
               CALL dbcsr_release_p(xas_tdp_env%q_projector(i)%matrix)
            END DO
            DEALLOCATE(xas_tdp_env%q_projector)
         END IF
         IF (ASSOCIATED(xas_tdp_env%dipmat)) THEN
            DO i = 1,SIZE(xas_tdp_env%dipmat)
               CALL dbcsr_release_p(xas_tdp_env%dipmat(i)%matrix)
            END DO
            DEALLOCATE(xas_tdp_env%dipmat)
         END IF
         IF (ASSOCIATED(xas_tdp_env%ri_o3c_coul)) THEN
            CALL release_o3c_container(xas_tdp_env%ri_o3c_coul)
            DEALLOCATE(xas_tdp_env%ri_o3c_coul)
         END IF
         IF (ASSOCIATED(xas_tdp_env%ri_o3c_ex)) THEN
            CALL release_o3c_container(xas_tdp_env%ri_o3c_ex)
            DEALLOCATE(xas_tdp_env%ri_o3c_ex)
         END IF
         IF (ASSOCIATED(xas_tdp_env%ri_fxc)) THEN
            DO i =1,SIZE(xas_tdp_env%ri_fxc, 1)
               DO j = 1,SIZE(xas_tdp_env%ri_fxc, 2)
                  IF (ASSOCIATED(xas_tdp_env%ri_fxc(i,j)%array)) THEN
                     DEALLOCATE(xas_tdp_env%ri_fxc(i,j)%array)
                  END IF
               END DO
            END DO
            DEALLOCATE(xas_tdp_env%ri_fxc)
         END IF
         IF (ASSOCIATED(xas_tdp_env%orb_soc)) THEN
            DO i = 1,SIZE(xas_tdp_env%orb_soc)
               CALL dbcsr_release(xas_tdp_env%orb_soc(i)%matrix)
               DEALLOCATE(xas_tdp_env%orb_soc(i)%matrix)
            END DO
            DEALLOCATE(xas_tdp_env%orb_soc)
         END IF
         DEALLOCATE (xas_tdp_env)
      END IF
   END SUBROUTINE xas_tdp_env_release

!> *************************************************************************************************
!> \brief Sets values of selected variables within the TDP XAS environment type
!> *************************************************************************************************
   SUBROUTINE set_xas_tdp_env(xas_tdp_env, nex_atoms, nex_kinds)

      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      INTEGER, INTENT(IN), OPTIONAL                      :: nex_atoms
      INTEGER, INTENT(IN), OPTIONAL                      :: nex_kinds

      CHARACTER(len=*), PARAMETER :: routineN = 'set_xas_tdp_env', &
         routineP = moduleN//':'//routineN

      CPASSERT(ASSOCIATED(xas_tdp_env))

      IF (PRESENT(nex_atoms)) xas_tdp_env%nex_atoms = nex_atoms
      IF (PRESENT(nex_kinds)) xas_tdp_env%nex_kinds = nex_kinds
      
   END SUBROUTINE set_xas_tdp_env

!> *************************************************************************************************
!> \brief Creates a donor_state
!> \param donor_state ...
!> *************************************************************************************************
   SUBROUTINE donor_state_create(donor_state)

      TYPE(donor_state_type), INTENT(INOUT)              :: donor_state

      CHARACTER(len=*), PARAMETER                        :: routineN = 'donor_state_create',&
                                                            routineP = moduleN//':'//routineN

      NULLIFY (donor_state%energy_evals)
      NULLIFY (donor_state%mo_indices)
      NULLIFY (donor_state%os_coeffs)
      NULLIFY (donor_state%sg_coeffs)
      NULLIFY (donor_state%tp_coeffs)
      NULLIFY (donor_state%gs_coeffs)
      NULLIFY (donor_state%os_evals)
      NULLIFY (donor_state%sg_evals)
      NULLIFY (donor_state%tp_evals)
      NULLIFY (donor_state%soc_evals)
      NULLIFY (donor_state%soc_osc_str)
      NULLIFY (donor_state%os_osc_str)
      NULLIFY (donor_state%sg_osc_str)
      NULLIFY (donor_state%tp_osc_str)
      NULLIFY (donor_state%os_matrix_tdp)
      NULLIFY (donor_state%sg_matrix_tdp)
      NULLIFY (donor_state%tp_matrix_tdp)
      NULLIFY (donor_state%metric)
      NULLIFY (donor_state%matrix_aux)

   END SUBROUTINE donor_state_create

!> *************************************************************************************************
!> \brief sets specified values of the donor state type
!> \param donor_state the type which values should be set
!> \param at_index ...
!> \param at_symbol ...
!> \param state_type ...
!> *************************************************************************************************
   SUBROUTINE set_donor_state(donor_state, at_index, at_symbol, kind_index, state_type)

      TYPE(donor_state_type), POINTER                    :: donor_state
      INTEGER, INTENT(IN), OPTIONAL                      :: at_index
      CHARACTER(len=2), INTENT(IN), OPTIONAL             :: at_symbol
      INTEGER, INTENT(IN), OPTIONAL                      :: kind_index, state_type

      CHARACTER(len=*), PARAMETER :: routineN = 'set_donor_state', &
         routineP = moduleN//':'//routineN

      CPASSERT(ASSOCIATED(donor_state))

      IF (PRESENT(at_index)) donor_state%at_index = at_index
      IF (PRESENT(kind_index)) donor_state%kind_index = kind_index
      IF (PRESENT(state_type)) donor_state%state_type = state_type
      IF (PRESENT(at_symbol)) donor_state%at_symbol = at_symbol

   END SUBROUTINE set_donor_state

!> *************************************************************************************************
!> \brief Deallocate a set of donor states
!> \param donor_state_set the set of donor states to deallocate
!> *************************************************************************************************
   SUBROUTINE deallocate_donor_state_set(donor_state_set)

      TYPE(donor_state_type), DIMENSION(:), POINTER      :: donor_state_set

      CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_donor_state_set', &
                                     routineP = moduleN//':'//routineN

      INTEGER                                            :: i, j

      IF (ASSOCIATED(donor_state_set)) THEN
         DO i = 1, SIZE(donor_state_set)

            IF (ASSOCIATED(donor_state_set(i)%os_coeffs)) THEN
               DO j = 1,SIZE(donor_state_set(i)%os_coeffs)
                  IF(ASSOCIATED(donor_state_set(i)%os_coeffs(j)%matrix)) THEN
                     CALL cp_fm_release(donor_state_set(i)%os_coeffs(j)%matrix)
                  END IF
               END DO
               DEALLOCATE(donor_state_set(i)%os_coeffs)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%sg_coeffs)) THEN
               DO j = 1,SIZE(donor_state_set(i)%sg_coeffs)
                  IF(ASSOCIATED(donor_state_set(i)%sg_coeffs(j)%matrix)) THEN
                     CALL cp_fm_release(donor_state_set(i)%sg_coeffs(j)%matrix)
                  END IF
               END DO
               DEALLOCATE(donor_state_set(i)%sg_coeffs)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%tp_coeffs)) THEN
               DO j = 1,SIZE(donor_state_set(i)%tp_coeffs)
                  IF(ASSOCIATED(donor_state_set(i)%tp_coeffs(j)%matrix)) THEN
                     CALL cp_fm_release(donor_state_set(i)%tp_coeffs(j)%matrix)
                  END IF
               END DO
               DEALLOCATE(donor_state_set(i)%tp_coeffs)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%gs_coeffs)) THEN
               CALL cp_fm_release(donor_state_set(i)%gs_coeffs)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%os_evals)) THEN
               DEALLOCATE(donor_state_set(i)%os_evals)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%sg_evals)) THEN
               DEALLOCATE(donor_state_set(i)%sg_evals)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%tp_evals)) THEN
               DEALLOCATE(donor_state_set(i)%tp_evals)
            END IF

             IF (ASSOCIATED(donor_state_set(i)%soc_evals)) THEN
               DEALLOCATE(donor_state_set(i)%soc_evals)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%os_osc_str)) THEN
               DEALLOCATE(donor_state_set(i)%os_osc_str)
            END IF
            
            IF (ASSOCIATED(donor_state_set(i)%sg_osc_str)) THEN
               DEALLOCATE(donor_state_set(i)%sg_osc_str)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%tp_osc_str)) THEN
               DEALLOCATE(donor_state_set(i)%tp_osc_str)
            END IF 

            IF (ASSOCIATED(donor_state_set(i)%soc_osc_str)) THEN
               DEALLOCATE(donor_state_set(i)%soc_osc_str)
            END IF

            If (ASSOCIATED(donor_state_set(i)%energy_evals)) THEN
               DEALLOCATE(donor_state_set(i)%energy_evals)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%mo_indices)) THEN
               DEALLOCATE (donor_state_set(i)%mo_indices)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%os_matrix_tdp)) THEN
               CALL dbcsr_release(donor_state_set(i)%os_matrix_tdp)
               DEALLOCATE(donor_state_set(i)%os_matrix_tdp)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%sg_matrix_tdp)) THEN
               CALL dbcsr_release(donor_state_set(i)%sg_matrix_tdp)
               DEALLOCATE(donor_state_set(i)%sg_matrix_tdp)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%tp_matrix_tdp)) THEN
               CALL dbcsr_release(donor_state_set(i)%tp_matrix_tdp)
               DEALLOCATE(donor_state_set(i)%tp_matrix_tdp)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%metric)) THEN
               DO j = 1,SIZE(donor_state_set(i)%metric)
                  IF (ASSOCIATED(donor_state_set(i)%metric(j)%matrix)) THEN  
                     CALL dbcsr_release(donor_state_set(i)%metric(j)%matrix)      
                     DEALLOCATE(donor_state_set(i)%metric(j)%matrix)           
                  END IF
               END DO
               DEALLOCATE(donor_state_set(i)%metric)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%matrix_aux)) THEN
               CALL dbcsr_release(donor_state_set(i)%matrix_aux)      
               DEALLOCATE(donor_state_set(i)%matrix_aux)
            END IF

         END DO
         DEALLOCATE (donor_state_set)
      END IF

   END SUBROUTINE deallocate_donor_state_set

!> *************************************************************************************************
!> \brief Deallocate a donor_state's heavy attributes
!> \param donor_state ...
!> *************************************************************************************************
   SUBROUTINE free_ds_memory(donor_state)
      
      TYPE(donor_state_type), POINTER                    :: donor_state

      CHARACTER(len=*), PARAMETER :: routineN = "free_ds_memory", routineP = moduleN//":"//routineN

      INTEGER                                            :: i

      IF (ASSOCIATED(donor_state%os_evals)) DEALLOCATE(donor_state%os_evals)
      IF (ASSOCIATED(donor_state%sg_evals)) DEALLOCATE(donor_state%sg_evals)
      IF (ASSOCIATED(donor_state%tp_evals)) DEALLOCATE(donor_state%tp_evals)
      IF (ASSOCIATED(donor_state%soc_evals)) DEALLOCATE(donor_state%soc_evals)
      IF (ASSOCIATED(donor_state%os_osc_str)) DEALLOCATE(donor_state%os_osc_str)
      IF (ASSOCIATED(donor_state%sg_osc_str)) DEALLOCATE(donor_state%sg_osc_str)
      IF (ASSOCIATED(donor_state%tp_osc_str)) DEALLOCATE(donor_state%tp_osc_str)
      IF (ASSOCIATED(donor_state%soc_osc_str)) DEALLOCATE(donor_state%soc_osc_str)
      IF (ASSOCIATED(donor_state%gs_coeffs)) CALL cp_fm_release(donor_state%gs_coeffs)

      IF (ASSOCIATED(donor_state%os_coeffs)) THEN
         DO i = 1,SIZE(donor_state%os_coeffs)
            CALL cp_fm_release(donor_state%os_coeffs(i)%matrix)
         END DO
         DEALLOCATE(donor_state%os_coeffs)
      END IF

      IF (ASSOCIATED(donor_state%sg_coeffs)) THEN
         DO i = 1,SIZE(donor_state%sg_coeffs)
            CALL cp_fm_release(donor_state%sg_coeffs(i)%matrix)
         END DO
         DEALLOCATE(donor_state%sg_coeffs)
      END IF   

      IF (ASSOCIATED(donor_state%tp_coeffs)) THEN
         DO i = 1,SIZE(donor_state%tp_coeffs)
            CALL cp_fm_release(donor_state%tp_coeffs(i)%matrix)
         END DO
         DEALLOCATE(donor_state%tp_coeffs)
      END IF

      IF (ASSOCIATED(donor_state%os_matrix_tdp)) THEN
         CALL dbcsr_release(donor_state%os_matrix_tdp)
         DEALLOCATE(donor_state%os_matrix_tdp)
      END IF

      IF (ASSOCIATED(donor_state%sg_matrix_tdp)) THEN
         CALL dbcsr_release(donor_state%sg_matrix_tdp)
         DEALLOCATE(donor_state%sg_matrix_tdp)
      END IF

      IF (ASSOCIATED(donor_state%tp_matrix_tdp)) THEN
         CALL dbcsr_release(donor_state%tp_matrix_tdp)
         DEALLOCATE(donor_state%tp_matrix_tdp)
      END IF

      IF (ASSOCIATED(donor_state%metric)) THEN
         DO i = 1,SIZE(donor_state%metric)
            IF (ASSOCIATED(donor_state%metric(i)%matrix)) THEN   
               CALL dbcsr_release(donor_state%metric(i)%matrix) 
               DEALLOCATE(donor_state%metric(i)%matrix)  
            END IF
         END DO
         DEALLOCATE(donor_state%metric)
      END IF

      IF (ASSOCIATED(donor_state%matrix_aux)) THEN
         CALL dbcsr_release(donor_state%matrix_aux) 
         DEALLOCATE(donor_state%matrix_aux)
      END IF

   END SUBROUTINE free_ds_memory


!> *************************************************************************************************
!> \brief Creates a xas_atom_env type 
!> \param xas_atom_env ...
!> *************************************************************************************************
   SUBROUTINE xas_atom_env_create(xas_atom_env)

      TYPE(xas_atom_env_type), POINTER                   :: xas_atom_env

      CHARACTER(len=*), PARAMETER :: routineN = "xas_atom_env_create", &
                                     routineP = moduleN//":"//routineN

      ALLOCATE(xas_atom_env)
      
      NULLIFY(xas_atom_env%excited_atoms)
      NULLIFY(xas_atom_env%excited_kinds)
      NULLIFY(xas_atom_env%proc_of_exat)
      NULLIFY(xas_atom_env%grid_atom_set)
      NULLIFY(xas_atom_env%harmonics_atom_set)
      NULLIFY(xas_atom_env%exat_neighbors)
      NULLIFY(xas_atom_env%ri_dcoeff)
      NULLIFY(xas_atom_env%kind_3c_overlap)
      NULLIFY(xas_atom_env%ri_sphi_so)
      NULLIFY(xas_atom_env%orb_sphi_so)

   END SUBROUTINE xas_atom_env_create

! **************************************************************************************************
!> \brief Releases the xas_atom_env type 
!> \param xas_atom_env the type to release
! **************************************************************************************************
   SUBROUTINE xas_atom_env_release(xas_atom_env)

      TYPE(xas_atom_env_type), POINTER                   :: xas_atom_env

      CHARACTER(len=*), PARAMETER :: routineN = "xas_atom_env_release", &
                                     routineP = moduleN//":"//routineN

      INTEGER                                            :: i, j

      IF (ASSOCIATED(xas_atom_env%grid_atom_set)) THEN 
         DO i = 1,SIZE(xas_atom_env%grid_atom_set)
            IF (ASSOCIATED(xas_atom_env%grid_atom_set(i)%grid_atom)) THEN
               CALL deallocate_grid_atom(xas_atom_env%grid_atom_set(i)%grid_atom)
            END IF
         END DO
         DEALLOCATE(xas_atom_env%grid_atom_set)
      END IF

      IF (ASSOCIATED(xas_atom_env%harmonics_atom_set)) THEN
         DO i = 1,SIZE(xas_atom_env%harmonics_atom_set)
            IF (ASSOCIATED(xas_atom_env%harmonics_atom_set(i)%harmonics_atom)) THEN 
               CALL deallocate_harmonics_atom(xas_atom_env%harmonics_atom_set(i)%harmonics_atom)
            END IF
         END DO
         DEALLOCATE(xas_atom_env%harmonics_atom_set)
      END IF

      ! Note that excited_atoms and excited_kinds are not deallocated because they point to other 
      ! ressources, namely xas_tdp_env

      IF (ASSOCIATED(xas_atom_env%proc_of_exat)) THEN
         DEALLOCATE(xas_atom_env%proc_of_exat)
      END IF

      IF (ASSOCIATED(xas_atom_env%exat_neighbors)) THEN
         DO i = 1,SIZE(xas_atom_env%exat_neighbors)
            IF (ASSOCIATED(xas_atom_env%exat_neighbors(i)%array)) THEN
               DEALLOCATE(xas_atom_env%exat_neighbors(i)%array)
            END IF
         END DO
         DEALLOCATE(xas_atom_env%exat_neighbors)
      END IF

      IF (ASSOCIATED(xas_atom_env%kind_3c_overlap)) THEN
         DO i = 1,SIZE(xas_atom_env%kind_3c_overlap)
            IF (ASSOCIATED(xas_atom_env%kind_3c_overlap(i)%array)) THEN
               DEALLOCATE(xas_atom_env%kind_3c_overlap(i)%array)
            END IF
         END DO
         DEALLOCATE(xas_atom_env%kind_3c_overlap)
      END IF

      IF (ASSOCIATED(xas_atom_env%ri_dcoeff)) THEN
         DO i = 1,SIZE(xas_atom_env%ri_dcoeff,1)
            DO j = 1,SIZE(xas_atom_env%ri_dcoeff,2)
               IF (ASSOCIATED(xas_atom_env%ri_dcoeff(i,j)%array)) THEN
                  DEALLOCATE(xas_atom_env%ri_dcoeff(i,j)%array)
               END IF
            END DO
         END DO
         DEALLOCATE(xas_atom_env%ri_dcoeff)
      END IF

      IF (ASSOCIATED(xas_atom_env%ri_sphi_so)) THEN
         DO i = 1,SIZE(xas_atom_env%ri_sphi_so)
            IF (ASSOCIATED(xas_atom_env%ri_sphi_so(i)%array)) THEN 
               DEALLOCATE(xas_atom_env%ri_sphi_so(i)%array)
            END IF
         END DO
         DEALLOCATE(xas_atom_env%ri_sphi_so)
      END IF

      IF (ASSOCIATED(xas_atom_env%orb_sphi_so)) THEN
         DO i = 1,SIZE(xas_atom_env%orb_sphi_so)
            IF (ASSOCIATED(xas_atom_env%orb_sphi_so(i)%array)) THEN 
               DEALLOCATE(xas_atom_env%orb_sphi_so(i)%array)
            END IF
         END DO
         DEALLOCATE(xas_atom_env%orb_sphi_so)
      END IF

      DEALLOCATE(xas_atom_env)

   END SUBROUTINE xas_atom_env_release

! **************************************************************************************************
!> \brief Releases the memory heavy attribute of xas_tdp_env that are specific to the current
!>        excited atom
!> \param xas_tdp_env ...
!> \param atom the index of the current excited atom
! **************************************************************************************************
   SUBROUTINE free_exat_memory(xas_tdp_env, atom)

      TYPE(xas_tdp_env_type), POINTER                 :: xas_tdp_env
      INTEGER, INTENT(IN)                             :: atom

      CHARACTER(len=*), PARAMETER :: routineN = "free_exat_memory", routineP = moduleN//":"//routineN

      INTEGER                                         :: i

      IF (ASSOCIATED(xas_tdp_env%ri_o3c_ex)) THEN
         CALL release_o3c_container(xas_tdp_env%ri_o3c_ex)
         DEALLOCATE(xas_tdp_env%ri_o3c_ex)
      END IF

      IF (ASSOCIATED(xas_tdp_env%ri_fxc)) THEN
         DO i = 1, SIZE(xas_tdp_env%ri_fxc,2)
            IF (ASSOCIATED(xas_tdp_env%ri_fxc(atom,i)%array)) THEN
               DEALLOCATE(xas_tdp_env%ri_fxc(atom,i)%array) 
            END IF
         END DO
      END IF

   END SUBROUTINE free_exat_memory

END MODULE xas_tdp_types
