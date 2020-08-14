!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

!> *************************************************************************************************
!> \brief Define XAS TDP control type and associated create, release, etc subroutines, as well as
!>        XAS TDP environment type and associated set, get, etc subroutines
!> \author AB (11.2017)
!> *************************************************************************************************
MODULE xas_tdp_types
   USE cp_array_utils,                  ONLY: cp_1d_i_p_type,&
                                              cp_1d_r_p_type,&
                                              cp_2d_i_p_type,&
                                              cp_2d_r_p_type,&
                                              cp_3d_r_p_type
   USE cp_files,                        ONLY: file_exists
   USE cp_fm_types,                     ONLY: cp_fm_release,&
                                              cp_fm_type
   USE cp_para_env,                     ONLY: cp_para_env_release
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_distribution_release,&
                                              dbcsr_distribution_type,&
                                              dbcsr_p_type,&
                                              dbcsr_release,&
                                              dbcsr_release_p,&
                                              dbcsr_type
   USE dbcsr_tensor_api,                ONLY: dbcsr_t_destroy,&
                                              dbcsr_t_type
   USE distribution_2d_types,           ONLY: distribution_2d_release,&
                                              distribution_2d_type
   USE input_constants,                 ONLY: &
        do_potential_coulomb, do_potential_short, do_potential_truncated, tddfpt_singlet, &
        tddfpt_spin_cons, tddfpt_spin_flip, tddfpt_triplet, xas_dip_vel, xas_tdp_by_index, &
        xas_tdp_by_kind
   USE input_section_types,             ONLY: section_vals_release,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE libint_2c_3c,                    ONLY: libint_potential_type
   USE libint_wrapper,                  ONLY: cp_libint_static_cleanup
   USE mathlib,                         ONLY: erfc_cutoff
   USE memory_utilities,                ONLY: reallocate
   USE physcon,                         ONLY: bohr,&
                                              evolt
   USE qs_grid_atom,                    ONLY: deallocate_grid_atom,&
                                              grid_atom_type
   USE qs_harmonics_atom,               ONLY: deallocate_harmonics_atom,&
                                              harmonics_atom_type
   USE qs_loc_types,                    ONLY: qs_loc_env_new_type,&
                                              qs_loc_env_release
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! **************************************************************************************************
!> \brief Type containing control information for TDP XAS calculations
!> \param define_excited whether excited atoms are chosen by kind or index
!> \param dipole_form whether the dipole moment is computed in the length or velocity representation
!> \param n_search # of lowest energy MOs to search for donor orbitals
!> \param check_only whether a check run for donor MOs is conducted
!> \param do_hfx whether exact exchange is included
!> \param do_xc wheter xc functional(s) is(are) included (libxc)
!> \param do_coulomb whether the coulomb kernel is computed, .FALSE. if no xc nor hfx => normal dft
!> \param sx the scaling applied to exact exchange
!> \param x_potential the potential used for exact exchange (incl. cutoff, t_c_file, omega)
!> \param ri_m_potential the potential used for exact exchange RI metric
!> \param do_ri_metric whether a metric is used fir the RI
!> \param eps_range the threshold to determine the effective range of the short range operator
!> \param eps_pgf the threshold to determine the extent of all pgf in the method
!> \param eps_filter threshold for dbcsr operations
!> \param ri_radius the radius of the sphere defining the neighbors in the RI projection of the dens
!> \param tamm_dancoff whether the calculations should be done in the Tamm-Dancoff approximation
!> \param do_quad whether the electric quadrupole transition moments should be computed
!> \param list_ex_atoms list of excited atom indices, kept empty if define_excited=by_kind
!> \param list_ex_kinds list of excited atom kinds, kept empty if define_excited=by_index
!> \param do_loc whether the core MOs should be localized
!> \param do_uks whether the calculation is spin-unrestricted
!> \param do_roks whether the calculation is restricted open-shell
!> \param do_singlet whether singlet excitations should be computed
!> \param do_triplet whether triplet excitations should be computed
!> \param do_spin_cons whether spin-conserving excitation (for open-shell) should be computed
!> \param do_spin_flip whether spin-flip excitation (for open-shell) should be computed
!> \param do_soc whether spin-orbit coupling should be included
!> \param n_excited the number of excited states to compute
!> \param e_range the energy range where to look for eigenvalues
!> \param state_types columns correspond to the states to excite for each atom kind/index
!>                    the number of rows is the number of times the keyword is repeated
!> \param grid_info the information about the atomic grids used for the xc kernel integrals
!> \param is_periodic self-explanatory
! **************************************************************************************************
   TYPE xas_tdp_control_type
      INTEGER                                 :: define_excited
      INTEGER                                 :: dipole_form
      INTEGER                                 :: n_search
      INTEGER                                 :: n_excited
      REAL(dp)                                :: e_range
      REAL(dp)                                :: sx
      REAL(dp)                                :: eps_range
      REAL(dp)                                :: eps_screen
      REAL(dp)                                :: eps_pgf
      REAL(dp)                                :: eps_filter
      TYPE(libint_potential_type)             :: x_potential
      TYPE(libint_potential_type)             :: ri_m_potential
      REAL(dp)                                :: ri_radius
      LOGICAL                                 :: do_hfx
      LOGICAL                                 :: do_xc
      LOGICAL                                 :: do_coulomb
      LOGICAL                                 :: do_ri_metric
      LOGICAL                                 :: check_only
      LOGICAL                                 :: tamm_dancoff
      LOGICAL                                 :: do_quad
      LOGICAL                                 :: do_loc
      LOGICAL                                 :: do_uks
      LOGICAL                                 :: do_roks
      LOGICAL                                 :: do_soc
      LOGICAL                                 :: do_singlet
      LOGICAL                                 :: do_triplet
      LOGICAL                                 :: do_spin_cons
      LOGICAL                                 :: do_spin_flip
      LOGICAL                                 :: is_periodic
      INTEGER, DIMENSION(:), POINTER          :: list_ex_atoms
      CHARACTER(len=default_string_length), &
         DIMENSION(:), POINTER    :: list_ex_kinds
      INTEGER, DIMENSION(:, :), POINTER        :: state_types
      TYPE(section_vals_type), POINTER        :: loc_subsection
      TYPE(section_vals_type), POINTER        :: print_loc_subsection
      CHARACTER(len=default_string_length), &
         DIMENSION(:, :), POINTER  :: grid_info

   END TYPE xas_tdp_control_type

!> *************************************************************************************************
!> \brief Type containing informations such as inputs and results for TDP XAS calculations
!> \param state_type_char an array containing the general donor state types as char (1s, 2s, 2p, ...)
!> \param nex_atoms number of excited atoms
!> \param nex_kinds number of excited kinds
!> \param ex_atom_indices array containing the indices of the excited atoms
!> \param ex_kind_indices array containing the indices of the excited kinds
!> \param state_types columns correspond to the different donor states of each excited atom
!> \param qs_loc_env the environment type dealing with the possible localization of donor orbitals
!> \param mos_of_ex_atoms links lowest energy MOs to excited atoms. Elements of value 1 mark the
!>        association between the MO irow and the excited atom icolumn. The third index is for spin
!> \param ri_inv_coul the inverse coulomb RI integral (P|Q)^-1, updated for each excited kind
!>        based on basis functions of the RI_XAS basis for that kind
!> \param ri_inv_ex the inverse exchange RI integral (P|Q)^-1, updated for each excited kind
!>        based on basis functions of the RI_XAS basis for that kind, and with the exchange operator
!>        Optionally, if a RI metric is present, contains M^-1 (P|Q) M^-1
!> \param q_projector the projector on the unperturbed, unoccupied ground state as a dbcsr matrix,
!>        for each spin
!> \param dipmat the dbcsr matrices containing the dipole in x,y,z directions evaluated on the
!>        contracted spherical gaussians. It can either be in the length or the velocity
!>        representation. For length representation, it has to be computed once with the origin on
!>        each excited atom
!> \param quadmat the dbcsr matrices containing the electric quadrupole in x2, xy, xz, y2, yz and z2
!>        directions in the AO basis. It is always in the length representation with the origin
!>        set to the current excited atom
!> \param ri_3c_coul the tensor containing the RI 3-cetner Coulomb integrals (computed once)
!> \param ri_3c_ex the tensor containing the RI 3-center exchange integrals (computed for each ex atom)
!> \param opt_dist2d_coul an optimized distribution_2d for localized Coulomb 3-center integrals
!> \param opt_dist2d_ex an optimized distribution_2d for localized exchange 3-center integrals
!> \param ri_fxc the array of xc integrals of type (P|fxc|Q), for alpha-alpha, alpha-beta and beta-beta
!> \param fxc_avail a boolean telling whwther fxc is availavle on all procs
!> \param orb_soc the matrix where the SOC is evaluated wrt the orbital basis set, for x,y,z
!> \param matrix_shalf the SQRT of the orbital overlap matrix, stored for PDOS use
!> *************************************************************************************************
   TYPE xas_tdp_env_type
      CHARACTER(len=2), DIMENSION(3)          :: state_type_char
      INTEGER                                 :: nex_atoms
      INTEGER                                 :: nex_kinds
      INTEGER, DIMENSION(:), POINTER          :: ex_atom_indices
      INTEGER, DIMENSION(:), POINTER          :: ex_kind_indices
      INTEGER, DIMENSION(:, :), POINTER       :: state_types
      TYPE(dbcsr_t_type), POINTER             :: ri_3c_coul
      TYPE(dbcsr_t_type), POINTER             :: ri_3c_ex
      TYPE(donor_state_type), DIMENSION(:), &
         POINTER   :: donor_states
      INTEGER, DIMENSION(:, :, :), POINTER      :: mos_of_ex_atoms
      TYPE(qs_loc_env_new_type), POINTER      :: qs_loc_env
      REAL(dp), DIMENSION(:, :), POINTER       :: ri_inv_coul
      REAL(dp), DIMENSION(:, :), POINTER       :: ri_inv_ex
      TYPE(distribution_2d_type), POINTER     :: opt_dist2d_coul
      TYPE(distribution_2d_type), POINTER     :: opt_dist2d_ex
      TYPE(dbcsr_p_type), DIMENSION(:), &
         POINTER   :: q_projector
      TYPE(dbcsr_p_type), DIMENSION(:), &
         POINTER   :: dipmat
      TYPE(dbcsr_p_type), DIMENSION(:), &
         POINTER   :: quadmat
      TYPE(cp_2d_r_p_type), DIMENSION(:, :), &
         POINTER   :: ri_fxc
      LOGICAL                                 :: fxc_avail
      TYPE(dbcsr_p_type), DIMENSION(:), &
         POINTER   :: orb_soc
      TYPE(cp_fm_type), POINTER               :: matrix_shalf
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
!> \param sc_coeffs solutions of the linear-response TDDFT equation for spin-conserving open-shell
!> \param sf_coeffs solutions of the linear-response TDDFT equation for spin-flip open-shell
!> \param sg_coeffs solutions of the linear-response TDDFT singlet equations
!> \param tp_coeffs solutions of the linear-response TDDFT triplet equations
!> \param gs_coeffs the ground state MO coefficients
!> \param contract_coeffs the subset of gs_coeffs centered on excited atom, used for RI contraction
!> \param sc_evals open-shell spin-conserving excitation energies
!> \param sf_evals open-shell spin-flip excitation energies
!> \param sg_evals singlet excitation energies => the eigenvalues of the linear response equation
!> \param tp_evals triplet excitation energies => the eigenvalues of the linear response equation
!> \param soc_evals excitation energies after inclusion of SOC
!> \param osc_str dipole oscilaltor strengths
!> \param soc_osc_str dipole oscillator strengths after the inclusion of SOC
!> \param quad_osc_str quadrupole oscilaltor strengths
!> \param soc_quad_osc_str quadrupole oscillator strengths after the inclusion of SOC
!> \param sc_matrix_tdp the dbcsr matrix to be diagonalized for open-shell spin-conserving calculations
!> \param sf_matrix_tdp the dbcsr matrix to be diagonalized for open-shell spin-flip calculations
!> \param sg_matrix_tdp the dbcsr matrix to be diagonalized to solve the problem for singlets
!> \param tp_matrix_tdp the dbcsr matrix to be diagonalized to solve the problem for triplets
!> \param metric the metric of the linear response problem M*c = omega*S*c and its inverse
!> \param matrix_aux the auxiliary matrix (A-D+E)^1/2 used to make the problem Hermitian
!> \param blk_size the col/row block size of the dbcsr matrices
!> \param dbcsr_dist the distribution of the dbcsr matrices
!> *************************************************************************************************
   TYPE donor_state_type
      INTEGER                                 :: at_index
      INTEGER                                 :: kind_index
      INTEGER                                 :: ndo_mo
      CHARACTER(LEN=default_string_length)    :: at_symbol
      INTEGER                                 :: state_type
      INTEGER, DIMENSION(:), POINTER          :: blk_size
      REAL(dp), DIMENSION(:, :), POINTER      :: energy_evals
      INTEGER, DIMENSION(:, :), POINTER       :: mo_indices
      TYPE(cp_fm_type), POINTER               :: sc_coeffs
      TYPE(cp_fm_type), POINTER               :: sf_coeffs
      TYPE(cp_fm_type), POINTER               :: sg_coeffs
      TYPE(cp_fm_type), POINTER               :: tp_coeffs
      TYPE(cp_fm_type), POINTER               :: gs_coeffs
      REAL(dp), DIMENSION(:, :), POINTER       :: contract_coeffs
      REAL(dp), DIMENSION(:), POINTER         :: sc_evals
      REAL(dp), DIMENSION(:), POINTER         :: sf_evals
      REAL(dp), DIMENSION(:), POINTER         :: sg_evals
      REAL(dp), DIMENSION(:), POINTER         :: tp_evals
      REAL(dp), DIMENSION(:), POINTER         :: soc_evals
      REAL(dp), DIMENSION(:), POINTER         :: osc_str
      REAL(dp), DIMENSION(:), POINTER         :: soc_osc_str
      REAL(dp), DIMENSION(:), POINTER         :: quad_osc_str
      REAL(dp), DIMENSION(:), POINTER         :: soc_quad_osc_str
      TYPE(dbcsr_type), POINTER               :: sc_matrix_tdp
      TYPE(dbcsr_type), POINTER               :: sf_matrix_tdp
      TYPE(dbcsr_type), POINTER               :: sg_matrix_tdp
      TYPE(dbcsr_type), POINTER               :: tp_matrix_tdp
      TYPE(dbcsr_p_type), DIMENSION(:), &
         POINTER   :: metric
      TYPE(dbcsr_type), POINTER               :: matrix_aux
      TYPE(dbcsr_distribution_type), POINTER  :: dbcsr_dist

   END TYPE donor_state_type

!  Some helper types for xas_tdp_atom
   TYPE grid_atom_p_type
      TYPE(grid_atom_type), POINTER                   :: grid_atom
   END TYPE grid_atom_p_type

   TYPE harmonics_atom_p_type
      TYPE(harmonics_atom_type), POINTER              :: harmonics_atom
   END TYPE harmonics_atom_p_type

   TYPE batch_info_type
      TYPE(cp_para_env_type), POINTER             :: para_env
      INTEGER                                     :: batch_size
      INTEGER                                     :: nbatch
      INTEGER                                     :: ibatch
      INTEGER                                     :: ipe
      INTEGER, DIMENSION(:), ALLOCATABLE          :: nso_proc
      INTEGER, DIMENSION(:, :), ALLOCATABLE       :: so_bo
      TYPE(cp_2d_i_p_type), POINTER, DIMENSION(:) :: so_proc_info
   END TYPE batch_info_type

! **************************************************************************************************
!> \brief a environment type that contains all the info needed for XAS_TDP atomic grid calculations
!> \param ri_radius defines the neighbors in the RI projection of the density
!> \param nspins ...
!> \param excited_atoms the atoms for which RI xc-kernel calculations must be done
!> \param excited_kinds the kinds for which RI xc-kernel calculations must be done
!> \param grid_atom_set the set of atomic grid for each kind
!> \param ri_dcoeff the expansion coefficients to express the density in the RI basis for each atom
!> \param exat_neighbors the neighbors of each excited atom
!> \param ri_sphi_so contains the coefficient for direct contraction from so to sgf, for the ri basis
!> \param orb_sphi_so contains the coefficient for direct contraction from so to sgf, for the orb basis
!> \param ga the angular part of the spherical gaussians on the grid of excited kinds
!> \param gr the radial part of the spherical gaussians on the grid of excited kinds
!> \param dgr1 first radial part of the gradient of the RI spherical gaussians
!> \param dgr2 second radial part of the gradient of the RI spherical gaussians
!> \param dga1 first angular part of the gradient of the RI spherical gaussians
!> \param dga2 second angular part of the gradient of the RI spherical gaussians
!> *************************************************************************************************
   TYPE xas_atom_env_type
      INTEGER                                         :: nspins
      REAL(dp)                                        :: ri_radius
      INTEGER, DIMENSION(:), POINTER                  :: excited_atoms
      INTEGER, DIMENSION(:), POINTER                  :: excited_kinds
      INTEGER, DIMENSION(:), POINTER                  :: proc_of_exat
      TYPE(grid_atom_p_type), DIMENSION(:), POINTER   :: grid_atom_set
      TYPE(harmonics_atom_p_type), DIMENSION(:), &
         POINTER  :: harmonics_atom_set
      TYPE(cp_1d_r_p_type), DIMENSION(:, :, :), POINTER :: ri_dcoeff
      TYPE(cp_2d_r_p_type), DIMENSION(:), POINTER     :: ri_sphi_so
      TYPE(cp_2d_r_p_type), DIMENSION(:), POINTER     :: orb_sphi_so
      TYPE(cp_1d_i_p_type), DIMENSION(:), POINTER     :: exat_neighbors
      TYPE(cp_2d_r_p_type), DIMENSION(:), POINTER     :: ga, gr, dgr1, dgr2
      TYPE(cp_3d_r_p_type), DIMENSION(:), POINTER     :: dga1, dga2
   END TYPE xas_atom_env_type

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xas_tdp_types'

! *** Public data types ***
   PUBLIC :: xas_tdp_env_type, donor_state_type, xas_tdp_control_type, xas_atom_env_type, &
             batch_info_type

! *** Public subroutines ***
   PUBLIC :: set_donor_state, free_ds_memory, release_batch_info, &
             xas_tdp_env_create, xas_tdp_env_release, set_xas_tdp_env, &
             xas_tdp_control_create, xas_tdp_control_release, read_xas_tdp_control, &
             xas_atom_env_create, xas_atom_env_release, donor_state_create, free_exat_memory, &
             get_proc_batch_sizes

CONTAINS

! **************************************************************************************************
!> \brief Creates and initializes the xas_tdp_control_type
!> \param xas_tdp_control the type to initialize
! **************************************************************************************************
   SUBROUTINE xas_tdp_control_create(xas_tdp_control)

      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control

      CHARACTER(len=*), PARAMETER :: routineN = 'xas_tdp_control_create'

      CPASSERT(.NOT. ASSOCIATED(xas_tdp_control))
      ALLOCATE (xas_tdp_control)

      xas_tdp_control%define_excited = xas_tdp_by_index
      xas_tdp_control%n_search = -1
      xas_tdp_control%dipole_form = xas_dip_vel
      xas_tdp_control%do_hfx = .FALSE.
      xas_tdp_control%do_xc = .FALSE.
      xas_tdp_control%do_coulomb = .TRUE.
      xas_tdp_control%do_ri_metric = .FALSE.
      xas_tdp_control%sx = 1.0_dp
      xas_tdp_control%eps_range = 1.0E-6_dp
      xas_tdp_control%eps_screen = 1.0E-10_dp
      xas_tdp_control%eps_pgf = -1.0_dp
      xas_tdp_control%eps_filter = 1.0E-10_dp
      xas_tdp_control%ri_radius = 0.0_dp
      xas_tdp_control%x_potential%potential_type = do_potential_coulomb
      xas_tdp_control%x_potential%cutoff_radius = 0.0_dp
      xas_tdp_control%x_potential%omega = 0.0_dp
      xas_tdp_control%x_potential%filename = " "
      xas_tdp_control%ri_m_potential%potential_type = do_potential_coulomb
      xas_tdp_control%ri_m_potential%cutoff_radius = 0.0_dp
      xas_tdp_control%ri_m_potential%omega = 0.0_dp
      xas_tdp_control%ri_m_potential%filename = " "
      xas_tdp_control%check_only = .FALSE.
      xas_tdp_control%tamm_dancoff = .FALSE.
      xas_tdp_control%do_quad = .FALSE.
      xas_tdp_control%do_loc = .FALSE.
      xas_tdp_control%do_uks = .FALSE.
      xas_tdp_control%do_roks = .FALSE.
      xas_tdp_control%do_soc = .FALSE.
      xas_tdp_control%do_singlet = .FALSE.
      xas_tdp_control%do_triplet = .FALSE.
      xas_tdp_control%do_spin_cons = .FALSE.
      xas_tdp_control%do_spin_flip = .FALSE.
      xas_tdp_control%is_periodic = .FALSE.
      xas_tdp_control%n_excited = -1
      xas_tdp_control%e_range = -1.0_dp
      NULLIFY (xas_tdp_control%state_types)
      NULLIFY (xas_tdp_control%list_ex_atoms)
      NULLIFY (xas_tdp_control%list_ex_kinds)
      NULLIFY (xas_tdp_control%loc_subsection)
      NULLIFY (xas_tdp_control%print_loc_subsection)
      NULLIFY (xas_tdp_control%grid_info)

   END SUBROUTINE xas_tdp_control_create

! **************************************************************************************************
!> \brief Releases the xas_tdp_control_type
!> \param xas_tdp_control the type to release
! **************************************************************************************************
   SUBROUTINE xas_tdp_control_release(xas_tdp_control)

      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control

      CHARACTER(len=*), PARAMETER :: routineN = 'xas_tdp_control_release'

      IF (ASSOCIATED(xas_tdp_control)) THEN
         IF (ASSOCIATED(xas_tdp_control%list_ex_atoms)) THEN
            DEALLOCATE (xas_tdp_control%list_ex_atoms)
         END IF
         IF (ASSOCIATED(xas_tdp_control%list_ex_kinds)) THEN
            DEALLOCATE (xas_tdp_control%list_ex_kinds)
         END IF
         IF (ASSOCIATED(xas_tdp_control%state_types)) THEN
            DEALLOCATE (xas_tdp_control%state_types)
         END IF
         IF (ASSOCIATED(xas_tdp_control%grid_info)) THEN
            DEALLOCATE (xas_tdp_control%grid_info)
         END IF
         IF (ASSOCIATED(xas_tdp_control%loc_subsection)) THEN
            !recursive, print_loc_subsection removed too
            CALL section_vals_release(xas_tdp_control%loc_subsection)
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

      CHARACTER(len=*), PARAMETER :: routineN = 'read_xas_tdp_control'

      CHARACTER(len=default_string_length), &
         DIMENSION(:), POINTER                           :: k_list
      INTEGER                                            :: excitation, irep, nexc, nrep
      INTEGER, DIMENSION(:), POINTER                     :: a_list, t_list

      NULLIFY (k_list, a_list, t_list)

!  Deal with the lone keywords

      CALL section_vals_val_get(xas_tdp_section, "CHECK_ONLY", &
                                l_val=xas_tdp_control%check_only)

      CALL section_vals_val_get(xas_tdp_section, "TAMM_DANCOFF", &
                                l_val=xas_tdp_control%tamm_dancoff)

      CALL section_vals_val_get(xas_tdp_section, "SPIN_ORBIT_COUPLING", &
                                l_val=xas_tdp_control%do_soc)

      CALL section_vals_val_get(xas_tdp_section, "DIPOLE_FORM", i_val=xas_tdp_control%dipole_form)

      CALL section_vals_val_get(xas_tdp_section, "QUADRUPOLE", l_val=xas_tdp_control%do_quad)

      CALL section_vals_val_get(xas_tdp_section, "EPS_PGF_XAS", n_rep_val=nrep)
      IF (nrep > 0) CALL section_vals_val_get(xas_tdp_section, "EPS_PGF_XAS", r_val=xas_tdp_control%eps_pgf)

      CALL section_vals_val_get(xas_tdp_section, "EPS_FILTER", r_val=xas_tdp_control%eps_filter)

      CALL section_vals_val_get(xas_tdp_section, "GRID", n_rep_val=nrep)

      IF (.NOT. ASSOCIATED(xas_tdp_control%grid_info)) THEN
         ALLOCATE (xas_tdp_control%grid_info(nrep, 3))
         DO irep = 1, nrep
            CALL section_vals_val_get(xas_tdp_section, "GRID", i_rep_val=irep, c_vals=k_list)
            IF (SIZE(k_list) .NE. 3) CPABORT("The GRID keyword needs three values")
            xas_tdp_control%grid_info(irep, :) = k_list
         END DO
      END IF

      CALL section_vals_val_get(xas_tdp_section, "EXCITATIONS", n_rep_val=nrep)
      DO irep = 1, nrep
         CALL section_vals_val_get(xas_tdp_section, "EXCITATIONS", i_rep_val=irep, i_val=excitation)
         IF (excitation == tddfpt_singlet) xas_tdp_control%do_singlet = .TRUE.
         IF (excitation == tddfpt_triplet) xas_tdp_control%do_triplet = .TRUE.
         IF (excitation == tddfpt_spin_cons) xas_tdp_control%do_spin_cons = .TRUE.
         IF (excitation == tddfpt_spin_flip) xas_tdp_control%do_spin_flip = .TRUE.
      END DO

      CALL section_vals_val_get(xas_tdp_section, "N_EXCITED", &
                                i_val=xas_tdp_control%n_excited)
      CALL section_vals_val_get(xas_tdp_section, "ENERGY_RANGE", &
                                r_val=xas_tdp_control%e_range)
      !store the range in Hartree, not eV
      xas_tdp_control%e_range = xas_tdp_control%e_range/evolt

!  Deal with the DONOR_STATES subsection

      CALL section_vals_val_get(xas_tdp_section, "DONOR_STATES%DEFINE_EXCITED", &
                                i_val=xas_tdp_control%define_excited)

      IF (.NOT. ASSOCIATED(xas_tdp_control%list_ex_kinds)) THEN
         IF (xas_tdp_control%define_excited .EQ. xas_tdp_by_index) THEN

            ALLOCATE (xas_tdp_control%list_ex_kinds(0))

         ELSE IF (xas_tdp_control%define_excited .EQ. xas_tdp_by_kind) THEN

            CALL section_vals_val_get(xas_tdp_section, "DONOR_STATES%KIND_LIST", c_vals=k_list)

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

            CALL section_vals_val_get(xas_tdp_section, "DONOR_STATES%ATOM_LIST", i_vals=a_list)

            IF (ASSOCIATED(a_list)) THEN
               nexc = SIZE(a_list)
               CALL reallocate(xas_tdp_control%list_ex_atoms, 1, nexc)
               xas_tdp_control%list_ex_atoms = a_list
            END IF

         END IF
      END IF

      CALL section_vals_val_get(xas_tdp_section, "DONOR_STATES%STATE_TYPES", n_rep_val=nrep)

      IF (.NOT. ASSOCIATED(xas_tdp_control%state_types)) THEN
         ALLOCATE (xas_tdp_control%state_types(nrep, nexc))
         DO irep = 1, nrep
            CALL section_vals_val_get(xas_tdp_section, "DONOR_STATES%STATE_TYPES", i_rep_val=irep, i_vals=t_list)
            IF (SIZE(t_list) .NE. nexc) THEN
               CPABORT("The STATE_TYPES keywords do not have the correct number of entries.")
            END IF
            xas_tdp_control%state_types(irep, :) = t_list
         END DO
      END IF
      IF (ALL(xas_tdp_control%state_types == 0)) CPABORT("Please specify STATE_TYPES")

      CALL section_vals_val_get(xas_tdp_section, "DONOR_STATES%N_SEARCH", i_val=xas_tdp_control%n_search)

      CALL section_vals_val_get(xas_tdp_section, "DONOR_STATES%LOCALIZE", l_val=xas_tdp_control%do_loc)

!  Deal with the KERNEL subsection
      CALL section_vals_val_get(xas_tdp_section, "KERNEL%XC_FUNCTIONAL%_SECTION_PARAMETERS_", &
                                l_val=xas_tdp_control%do_xc)
      CALL section_vals_val_get(xas_tdp_section, "KERNEL%EXACT_EXCHANGE%_SECTION_PARAMETERS_", &
                                l_val=xas_tdp_control%do_hfx)

      CALL section_vals_val_get(xas_tdp_section, "KERNEL%RI_REGION", r_val=xas_tdp_control%ri_radius)
      xas_tdp_control%ri_radius = bohr*xas_tdp_control%ri_radius

      IF (xas_tdp_control%do_hfx) THEN
         !The main exact echange potential and related params
         CALL section_vals_val_get(xas_tdp_section, "KERNEL%EXACT_EXCHANGE%SCALE", &
                                   r_val=xas_tdp_control%sx)
         CALL section_vals_val_get(xas_tdp_section, "KERNEL%EXACT_EXCHANGE%POTENTIAL_TYPE", &
                                   i_val=xas_tdp_control%x_potential%potential_type)
         !truncated Coulomb
         IF (xas_tdp_control%x_potential%potential_type == do_potential_truncated) THEN
            CALL section_vals_val_get(xas_tdp_section, "KERNEL%EXACT_EXCHANGE%T_C_G_DATA", &
                                      c_val=xas_tdp_control%x_potential%filename)
            IF (.NOT. file_exists(xas_tdp_control%x_potential%filename)) THEN
               CPABORT("Could not find provided T_C_G_DATA file.")
            END IF
            CALL section_vals_val_get(xas_tdp_section, "KERNEL%EXACT_EXCHANGE%CUTOFF_RADIUS", &
                                      r_val=xas_tdp_control%x_potential%cutoff_radius)
            !store the range in bohrs
            xas_tdp_control%x_potential%cutoff_radius = bohr*xas_tdp_control%x_potential%cutoff_radius
         END IF

         !short range erfc
         IF (xas_tdp_control%x_potential%potential_type == do_potential_short) THEN
            CALL section_vals_val_get(xas_tdp_section, "KERNEL%EXACT_EXCHANGE%OMEGA", &
                                      r_val=xas_tdp_control%x_potential%omega)
            CALL section_vals_val_get(xas_tdp_section, "KERNEL%EXACT_EXCHANGE%EPS_RANGE", &
                                      r_val=xas_tdp_control%eps_range)
            !get the effective range (omega in 1/a0, range in a0)
            CALL erfc_cutoff(xas_tdp_control%eps_range, xas_tdp_control%x_potential%omega, &
                             xas_tdp_control%x_potential%cutoff_radius)

         END IF

         CALL section_vals_val_get(xas_tdp_section, "KERNEL%EXACT_EXCHANGE%EPS_SCREENING", &
                                   r_val=xas_tdp_control%eps_screen)
         !The RI metric stuff
         CALL section_vals_val_get(xas_tdp_section, "KERNEL%EXACT_EXCHANGE%RI_METRIC%_SECTION_PARAMETERS_", &
                                   l_val=xas_tdp_control%do_ri_metric)
         IF (xas_tdp_control%do_ri_metric) THEN

            CALL section_vals_val_get(xas_tdp_section, "KERNEL%EXACT_EXCHANGE%RI_METRIC%POTENTIAL_TYPE", &
                                      i_val=xas_tdp_control%ri_m_potential%potential_type)

            !truncated Coulomb
            IF (xas_tdp_control%ri_m_potential%potential_type == do_potential_truncated) THEN
               CALL section_vals_val_get(xas_tdp_section, "KERNEL%EXACT_EXCHANGE%RI_METRIC%T_C_G_DATA", &
                                         c_val=xas_tdp_control%ri_m_potential%filename)
               IF (.NOT. file_exists(xas_tdp_control%ri_m_potential%filename)) THEN
                  CPABORT("Could not find provided T_C_G_DATA file.")
               END IF
               CALL section_vals_val_get(xas_tdp_section, "KERNEL%EXACT_EXCHANGE%RI_METRIC%CUTOFF_RADIUS", &
                                         r_val=xas_tdp_control%ri_m_potential%cutoff_radius)
               !store the range in bohrs
               xas_tdp_control%ri_m_potential%cutoff_radius = bohr*xas_tdp_control%ri_m_potential%cutoff_radius
            END IF

            !short range erfc
            IF (xas_tdp_control%ri_m_potential%potential_type == do_potential_short) THEN
               CALL section_vals_val_get(xas_tdp_section, "KERNEL%EXACT_EXCHANGE%RI_METRIC%OMEGA", &
                                         r_val=xas_tdp_control%ri_m_potential%omega)
               !get the effective range (omega in 1/a0, range in a0)
               CALL erfc_cutoff(xas_tdp_control%eps_range, xas_tdp_control%ri_m_potential%omega, &
                                xas_tdp_control%ri_m_potential%cutoff_radius)

            END IF
         ELSE
            !No defined metric, V-approximation, set all ri_m_potential params to those of x_pot
            xas_tdp_control%ri_m_potential = xas_tdp_control%x_potential

         END IF

      END IF

      IF ((.NOT. xas_tdp_control%do_xc) .AND. (.NOT. xas_tdp_control%do_hfx)) THEN
         !then no coulomb either and go full DFT
         xas_tdp_control%do_coulomb = .FALSE.
      END IF

   END SUBROUTINE read_xas_tdp_control

! **************************************************************************************************
!> \brief Creates a TDP XAS environment type
!> \param xas_tdp_env the type to create
! **************************************************************************************************
   SUBROUTINE xas_tdp_env_create(xas_tdp_env)

      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env

      CHARACTER(len=*), PARAMETER :: routineN = 'xas_tdp_env_create'

      ALLOCATE (xas_tdp_env)

      xas_tdp_env%nex_atoms = 1
      xas_tdp_env%nex_kinds = 1
      xas_tdp_env%fxc_avail = .FALSE.

      NULLIFY (xas_tdp_env%ex_atom_indices)
      NULLIFY (xas_tdp_env%ex_kind_indices)
      NULLIFY (xas_tdp_env%state_types)
      NULLIFY (xas_tdp_env%donor_states)
      NULLIFY (xas_tdp_env%qs_loc_env)
      NULLIFY (xas_tdp_env%mos_of_ex_atoms)
      NULLIFY (xas_tdp_env%ri_inv_coul)
      NULLIFY (xas_tdp_env%ri_inv_ex)
      NULLIFY (xas_tdp_env%opt_dist2d_coul)
      NULLIFY (xas_tdp_env%opt_dist2d_ex)
      NULLIFY (xas_tdp_env%q_projector)
      NULLIFY (xas_tdp_env%dipmat)
      NULLIFY (xas_tdp_env%quadmat)
      NULLIFY (xas_tdp_env%ri_3c_coul)
      NULLIFY (xas_tdp_env%ri_3c_ex)
      NULLIFY (xas_tdp_env%ri_fxc)
      NULLIFY (xas_tdp_env%orb_soc)
      NULLIFY (xas_tdp_env%matrix_shalf)

!     Putting the state types as char manually
      xas_tdp_env%state_type_char(1) = "1s"
      xas_tdp_env%state_type_char(2) = "2s"
      xas_tdp_env%state_type_char(3) = "2p"

   END SUBROUTINE xas_tdp_env_create

! **************************************************************************************************
!> \brief Releases the TDP XAS environment type
!> \param xas_tdp_env the type to release
! **************************************************************************************************
   SUBROUTINE xas_tdp_env_release(xas_tdp_env)

      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env

      CHARACTER(len=*), PARAMETER :: routineN = 'xas_tdp_env_release'

      INTEGER                                            :: i, j

      IF (ASSOCIATED(xas_tdp_env)) THEN
         IF (ASSOCIATED(xas_tdp_env%ex_atom_indices)) THEN
            DEALLOCATE (xas_tdp_env%ex_atom_indices)
         END IF
         IF (ASSOCIATED(xas_tdp_env%ex_kind_indices)) THEN
            DEALLOCATE (xas_tdp_env%ex_kind_indices)
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
            DEALLOCATE (xas_tdp_env%mos_of_ex_atoms)
         END IF
         IF (ASSOCIATED(xas_tdp_env%ri_inv_coul)) THEN
            DEALLOCATE (xas_tdp_env%ri_inv_coul)
         END IF
         IF (ASSOCIATED(xas_tdp_env%ri_inv_ex)) THEN
            DEALLOCATE (xas_tdp_env%ri_inv_ex)
         END IF
         IF (ASSOCIATED(xas_tdp_env%opt_dist2d_coul)) THEN
            CALL distribution_2d_release(xas_tdp_env%opt_dist2d_coul)
         END IF
         IF (ASSOCIATED(xas_tdp_env%opt_dist2d_ex)) THEN
            CALL distribution_2d_release(xas_tdp_env%opt_dist2d_ex)
         END IF
         IF (ASSOCIATED(xas_tdp_env%q_projector)) THEN
            DO i = 1, SIZE(xas_tdp_env%q_projector)
               CALL dbcsr_release_p(xas_tdp_env%q_projector(i)%matrix)
            END DO
            DEALLOCATE (xas_tdp_env%q_projector)
         END IF
         IF (ASSOCIATED(xas_tdp_env%dipmat)) THEN
            DO i = 1, SIZE(xas_tdp_env%dipmat)
               CALL dbcsr_release_p(xas_tdp_env%dipmat(i)%matrix)
            END DO
            DEALLOCATE (xas_tdp_env%dipmat)
         END IF
         IF (ASSOCIATED(xas_tdp_env%quadmat)) THEN
            DO i = 1, SIZE(xas_tdp_env%quadmat)
               CALL dbcsr_release_p(xas_tdp_env%quadmat(i)%matrix)
            END DO
            DEALLOCATE (xas_tdp_env%quadmat)
         END IF
         IF (ASSOCIATED(xas_tdp_env%ri_3c_coul)) THEN
            CALL dbcsr_t_destroy(xas_tdp_env%ri_3c_coul)
            DEALLOCATE (xas_tdp_env%ri_3c_coul)
         END IF
         IF (ASSOCIATED(xas_tdp_env%ri_3c_ex)) THEN
            CALL dbcsr_t_destroy(xas_tdp_env%ri_3c_ex)
            DEALLOCATE (xas_tdp_env%ri_3c_ex)
         END IF
         IF (ASSOCIATED(xas_tdp_env%ri_fxc)) THEN
            DO i = 1, SIZE(xas_tdp_env%ri_fxc, 1)
               DO j = 1, SIZE(xas_tdp_env%ri_fxc, 2)
                  IF (ASSOCIATED(xas_tdp_env%ri_fxc(i, j)%array)) THEN
                     DEALLOCATE (xas_tdp_env%ri_fxc(i, j)%array)
                  END IF
               END DO
            END DO
            DEALLOCATE (xas_tdp_env%ri_fxc)
         END IF
         IF (ASSOCIATED(xas_tdp_env%orb_soc)) THEN
            DO i = 1, SIZE(xas_tdp_env%orb_soc)
               CALL dbcsr_release(xas_tdp_env%orb_soc(i)%matrix)
               DEALLOCATE (xas_tdp_env%orb_soc(i)%matrix)
            END DO
            DEALLOCATE (xas_tdp_env%orb_soc)
         END IF
         IF (ASSOCIATED(xas_tdp_env%matrix_shalf)) THEN
            CALL cp_fm_release(xas_tdp_env%matrix_shalf)
         END IF
         DEALLOCATE (xas_tdp_env)
      END IF
   END SUBROUTINE xas_tdp_env_release

! **************************************************************************************************
!> \brief Sets values of selected variables within the TDP XAS environment type
!> \param xas_tdp_env ...
!> \param nex_atoms ...
!> \param nex_kinds ...
! **************************************************************************************************
   SUBROUTINE set_xas_tdp_env(xas_tdp_env, nex_atoms, nex_kinds)

      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      INTEGER, INTENT(IN), OPTIONAL                      :: nex_atoms, nex_kinds

      CHARACTER(len=*), PARAMETER                        :: routineN = 'set_xas_tdp_env'

      CPASSERT(ASSOCIATED(xas_tdp_env))

      IF (PRESENT(nex_atoms)) xas_tdp_env%nex_atoms = nex_atoms
      IF (PRESENT(nex_kinds)) xas_tdp_env%nex_kinds = nex_kinds

   END SUBROUTINE set_xas_tdp_env

! **************************************************************************************************
!> \brief Creates a donor_state
!> \param donor_state ...
! **************************************************************************************************
   SUBROUTINE donor_state_create(donor_state)

      TYPE(donor_state_type), INTENT(INOUT)              :: donor_state

      CHARACTER(len=*), PARAMETER :: routineN = 'donor_state_create'

      NULLIFY (donor_state%energy_evals)
      NULLIFY (donor_state%mo_indices)
      NULLIFY (donor_state%sc_coeffs)
      NULLIFY (donor_state%sf_coeffs)
      NULLIFY (donor_state%sg_coeffs)
      NULLIFY (donor_state%tp_coeffs)
      NULLIFY (donor_state%gs_coeffs)
      NULLIFY (donor_state%contract_coeffs)
      NULLIFY (donor_state%sc_evals)
      NULLIFY (donor_state%sf_evals)
      NULLIFY (donor_state%sg_evals)
      NULLIFY (donor_state%tp_evals)
      NULLIFY (donor_state%soc_evals)
      NULLIFY (donor_state%soc_osc_str)
      NULLIFY (donor_state%osc_str)
      NULLIFY (donor_state%soc_quad_osc_str)
      NULLIFY (donor_state%quad_osc_str)
      NULLIFY (donor_state%sc_matrix_tdp)
      NULLIFY (donor_state%sf_matrix_tdp)
      NULLIFY (donor_state%sg_matrix_tdp)
      NULLIFY (donor_state%tp_matrix_tdp)
      NULLIFY (donor_state%metric)
      NULLIFY (donor_state%matrix_aux)
      NULLIFY (donor_state%blk_size)
      NULLIFY (donor_state%dbcsr_dist)

   END SUBROUTINE donor_state_create

! **************************************************************************************************
!> \brief sets specified values of the donor state type
!> \param donor_state the type which values should be set
!> \param at_index ...
!> \param at_symbol ...
!> \param kind_index ...
!> \param state_type ...
! **************************************************************************************************
   SUBROUTINE set_donor_state(donor_state, at_index, at_symbol, kind_index, state_type)

      TYPE(donor_state_type), POINTER                    :: donor_state
      INTEGER, INTENT(IN), OPTIONAL                      :: at_index
      CHARACTER(LEN=default_string_length), OPTIONAL     :: at_symbol
      INTEGER, INTENT(IN), OPTIONAL                      :: kind_index, state_type

      CHARACTER(len=*), PARAMETER                        :: routineN = 'set_donor_state'

      CPASSERT(ASSOCIATED(donor_state))

      IF (PRESENT(at_index)) donor_state%at_index = at_index
      IF (PRESENT(kind_index)) donor_state%kind_index = kind_index
      IF (PRESENT(state_type)) donor_state%state_type = state_type
      IF (PRESENT(at_symbol)) donor_state%at_symbol = at_symbol

   END SUBROUTINE set_donor_state

! **************************************************************************************************
!> \brief Deallocate a set of donor states
!> \param donor_state_set the set of donor states to deallocate
! **************************************************************************************************
   SUBROUTINE deallocate_donor_state_set(donor_state_set)
      TYPE(donor_state_type), DIMENSION(:), POINTER      :: donor_state_set

      CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_donor_state_set'

      INTEGER                                            :: i, j

      IF (ASSOCIATED(donor_state_set)) THEN
         DO i = 1, SIZE(donor_state_set)

            IF (ASSOCIATED(donor_state_set(i)%sc_coeffs)) THEN
               CALL cp_fm_release(donor_state_set(i)%sc_coeffs)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%sf_coeffs)) THEN
               CALL cp_fm_release(donor_state_set(i)%sf_coeffs)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%sg_coeffs)) THEN
               CALL cp_fm_release(donor_state_set(i)%sg_coeffs)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%tp_coeffs)) THEN
               CALL cp_fm_release(donor_state_set(i)%tp_coeffs)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%gs_coeffs)) THEN
               CALL cp_fm_release(donor_state_set(i)%gs_coeffs)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%contract_coeffs)) THEN
               DEALLOCATE (donor_state_set(i)%contract_coeffs)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%sc_evals)) THEN
               DEALLOCATE (donor_state_set(i)%sc_evals)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%sf_evals)) THEN
               DEALLOCATE (donor_state_set(i)%sf_evals)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%sg_evals)) THEN
               DEALLOCATE (donor_state_set(i)%sg_evals)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%tp_evals)) THEN
               DEALLOCATE (donor_state_set(i)%tp_evals)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%soc_evals)) THEN
               DEALLOCATE (donor_state_set(i)%soc_evals)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%osc_str)) THEN
               DEALLOCATE (donor_state_set(i)%osc_str)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%soc_osc_str)) THEN
               DEALLOCATE (donor_state_set(i)%soc_osc_str)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%quad_osc_str)) THEN
               DEALLOCATE (donor_state_set(i)%quad_osc_str)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%soc_quad_osc_str)) THEN
               DEALLOCATE (donor_state_set(i)%soc_quad_osc_str)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%energy_evals)) THEN
               DEALLOCATE (donor_state_set(i)%energy_evals)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%mo_indices)) THEN
               DEALLOCATE (donor_state_set(i)%mo_indices)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%sc_matrix_tdp)) THEN
               CALL dbcsr_release(donor_state_set(i)%sc_matrix_tdp)
               DEALLOCATE (donor_state_set(i)%sc_matrix_tdp)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%sf_matrix_tdp)) THEN
               CALL dbcsr_release(donor_state_set(i)%sf_matrix_tdp)
               DEALLOCATE (donor_state_set(i)%sf_matrix_tdp)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%sg_matrix_tdp)) THEN
               CALL dbcsr_release(donor_state_set(i)%sg_matrix_tdp)
               DEALLOCATE (donor_state_set(i)%sg_matrix_tdp)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%tp_matrix_tdp)) THEN
               CALL dbcsr_release(donor_state_set(i)%tp_matrix_tdp)
               DEALLOCATE (donor_state_set(i)%tp_matrix_tdp)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%metric)) THEN
               DO j = 1, SIZE(donor_state_set(i)%metric)
                  IF (ASSOCIATED(donor_state_set(i)%metric(j)%matrix)) THEN
                     CALL dbcsr_release(donor_state_set(i)%metric(j)%matrix)
                     DEALLOCATE (donor_state_set(i)%metric(j)%matrix)
                  END IF
               END DO
               DEALLOCATE (donor_state_set(i)%metric)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%matrix_aux)) THEN
               CALL dbcsr_release(donor_state_set(i)%matrix_aux)
               DEALLOCATE (donor_state_set(i)%matrix_aux)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%blk_size)) THEN
               DEALLOCATE (donor_state_set(i)%blk_size)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%dbcsr_dist)) THEN
               CALL dbcsr_distribution_release(donor_state_set(i)%dbcsr_dist)
               DEALLOCATE (donor_state_set(i)%dbcsr_dist)
            END IF
         END DO
         DEALLOCATE (donor_state_set)
      END IF

   END SUBROUTINE deallocate_donor_state_set

! **************************************************************************************************
!> \brief Deallocate a donor_state's heavy attributes
!> \param donor_state ...
! **************************************************************************************************
   SUBROUTINE free_ds_memory(donor_state)

      TYPE(donor_state_type), POINTER                    :: donor_state

      CHARACTER(len=*), PARAMETER                        :: routineN = 'free_ds_memory'

      INTEGER                                            :: i

      IF (ASSOCIATED(donor_state%sc_evals)) DEALLOCATE (donor_state%sc_evals)
      IF (ASSOCIATED(donor_state%contract_coeffs)) DEALLOCATE (donor_state%contract_coeffs)
      IF (ASSOCIATED(donor_state%sf_evals)) DEALLOCATE (donor_state%sf_evals)
      IF (ASSOCIATED(donor_state%sg_evals)) DEALLOCATE (donor_state%sg_evals)
      IF (ASSOCIATED(donor_state%tp_evals)) DEALLOCATE (donor_state%tp_evals)
      IF (ASSOCIATED(donor_state%soc_evals)) DEALLOCATE (donor_state%soc_evals)
      IF (ASSOCIATED(donor_state%osc_str)) DEALLOCATE (donor_state%osc_str)
      IF (ASSOCIATED(donor_state%soc_osc_str)) DEALLOCATE (donor_state%soc_osc_str)
      IF (ASSOCIATED(donor_state%quad_osc_str)) DEALLOCATE (donor_state%quad_osc_str)
      IF (ASSOCIATED(donor_state%soc_quad_osc_str)) DEALLOCATE (donor_state%soc_quad_osc_str)
      IF (ASSOCIATED(donor_state%gs_coeffs)) CALL cp_fm_release(donor_state%gs_coeffs)
      IF (ASSOCIATED(donor_state%blk_size)) DEALLOCATE (donor_state%blk_size)

      IF (ASSOCIATED(donor_state%sc_coeffs)) THEN
         CALL cp_fm_release(donor_state%sc_coeffs)
      END IF

      IF (ASSOCIATED(donor_state%sf_coeffs)) THEN
         CALL cp_fm_release(donor_state%sf_coeffs)
      END IF

      IF (ASSOCIATED(donor_state%sg_coeffs)) THEN
         CALL cp_fm_release(donor_state%sg_coeffs)
      END IF

      IF (ASSOCIATED(donor_state%tp_coeffs)) THEN
         CALL cp_fm_release(donor_state%tp_coeffs)
      END IF

      IF (ASSOCIATED(donor_state%sc_matrix_tdp)) THEN
         CALL dbcsr_release(donor_state%sc_matrix_tdp)
         DEALLOCATE (donor_state%sc_matrix_tdp)
      END IF

      IF (ASSOCIATED(donor_state%sf_matrix_tdp)) THEN
         CALL dbcsr_release(donor_state%sf_matrix_tdp)
         DEALLOCATE (donor_state%sf_matrix_tdp)
      END IF

      IF (ASSOCIATED(donor_state%sg_matrix_tdp)) THEN
         CALL dbcsr_release(donor_state%sg_matrix_tdp)
         DEALLOCATE (donor_state%sg_matrix_tdp)
      END IF

      IF (ASSOCIATED(donor_state%tp_matrix_tdp)) THEN
         CALL dbcsr_release(donor_state%tp_matrix_tdp)
         DEALLOCATE (donor_state%tp_matrix_tdp)
      END IF

      IF (ASSOCIATED(donor_state%metric)) THEN
         DO i = 1, SIZE(donor_state%metric)
            IF (ASSOCIATED(donor_state%metric(i)%matrix)) THEN
               CALL dbcsr_release(donor_state%metric(i)%matrix)
               DEALLOCATE (donor_state%metric(i)%matrix)
            END IF
         END DO
         DEALLOCATE (donor_state%metric)
      END IF

      IF (ASSOCIATED(donor_state%matrix_aux)) THEN
         CALL dbcsr_release(donor_state%matrix_aux)
         DEALLOCATE (donor_state%matrix_aux)
      END IF

      IF (ASSOCIATED(donor_state%dbcsr_dist)) THEN
         CALL dbcsr_distribution_release(donor_state%dbcsr_dist)
         DEALLOCATE (donor_state%dbcsr_dist)
      END IF

   END SUBROUTINE free_ds_memory

! **************************************************************************************************
!> \brief Creates a xas_atom_env type
!> \param xas_atom_env ...
! **************************************************************************************************
   SUBROUTINE xas_atom_env_create(xas_atom_env)

      TYPE(xas_atom_env_type), POINTER                   :: xas_atom_env

      CHARACTER(len=*), PARAMETER :: routineN = 'xas_atom_env_create'

      ALLOCATE (xas_atom_env)

      xas_atom_env%nspins = 1
      xas_atom_env%ri_radius = 0.0_dp
      NULLIFY (xas_atom_env%excited_atoms)
      NULLIFY (xas_atom_env%excited_kinds)
      NULLIFY (xas_atom_env%grid_atom_set)
      NULLIFY (xas_atom_env%harmonics_atom_set)
      NULLIFY (xas_atom_env%ri_dcoeff)
      NULLIFY (xas_atom_env%ri_sphi_so)
      NULLIFY (xas_atom_env%orb_sphi_so)
      NULLIFY (xas_atom_env%exat_neighbors)
      NULLIFY (xas_atom_env%gr)
      NULLIFY (xas_atom_env%ga)
      NULLIFY (xas_atom_env%dgr1)
      NULLIFY (xas_atom_env%dgr2)
      NULLIFY (xas_atom_env%dga1)
      NULLIFY (xas_atom_env%dga2)

   END SUBROUTINE xas_atom_env_create

! **************************************************************************************************
!> \brief Releases the xas_atom_env type
!> \param xas_atom_env the type to release
! **************************************************************************************************
   SUBROUTINE xas_atom_env_release(xas_atom_env)

      TYPE(xas_atom_env_type), POINTER                   :: xas_atom_env

      CHARACTER(len=*), PARAMETER :: routineN = 'xas_atom_env_release'

      INTEGER                                            :: i, j, k

      IF (ASSOCIATED(xas_atom_env%grid_atom_set)) THEN
         DO i = 1, SIZE(xas_atom_env%grid_atom_set)
            IF (ASSOCIATED(xas_atom_env%grid_atom_set(i)%grid_atom)) THEN
               CALL deallocate_grid_atom(xas_atom_env%grid_atom_set(i)%grid_atom)
            END IF
         END DO
         DEALLOCATE (xas_atom_env%grid_atom_set)
      END IF

      IF (ASSOCIATED(xas_atom_env%harmonics_atom_set)) THEN
         DO i = 1, SIZE(xas_atom_env%harmonics_atom_set)
            IF (ASSOCIATED(xas_atom_env%harmonics_atom_set(i)%harmonics_atom)) THEN
               CALL deallocate_harmonics_atom(xas_atom_env%harmonics_atom_set(i)%harmonics_atom)
            END IF
         END DO
         DEALLOCATE (xas_atom_env%harmonics_atom_set)
      END IF

      ! Note that excited_atoms and excited_kinds are not deallocated because they point to other
      ! ressources, namely xas_tdp_env.

      IF (ASSOCIATED(xas_atom_env%ri_dcoeff)) THEN
         DO i = 1, SIZE(xas_atom_env%ri_dcoeff, 1)
            DO j = 1, SIZE(xas_atom_env%ri_dcoeff, 2)
               DO k = 1, SIZE(xas_atom_env%ri_dcoeff, 3)
                  IF (ASSOCIATED(xas_atom_env%ri_dcoeff(i, j, k)%array)) THEN
                     DEALLOCATE (xas_atom_env%ri_dcoeff(i, j, k)%array)
                  END IF
               END DO
            END DO
         END DO
         DEALLOCATE (xas_atom_env%ri_dcoeff)
      END IF

      IF (ASSOCIATED(xas_atom_env%ri_sphi_so)) THEN
         DO i = 1, SIZE(xas_atom_env%ri_sphi_so)
            IF (ASSOCIATED(xas_atom_env%ri_sphi_so(i)%array)) THEN
               DEALLOCATE (xas_atom_env%ri_sphi_so(i)%array)
            END IF
         END DO
         DEALLOCATE (xas_atom_env%ri_sphi_so)
      END IF

      IF (ASSOCIATED(xas_atom_env%exat_neighbors)) THEN
         DO i = 1, SIZE(xas_atom_env%exat_neighbors)
            IF (ASSOCIATED(xas_atom_env%exat_neighbors(i)%array)) THEN
               DEALLOCATE (xas_atom_env%exat_neighbors(i)%array)
            END IF
         END DO
         DEALLOCATE (xas_atom_env%exat_neighbors)
      END IF

      IF (ASSOCIATED(xas_atom_env%gr)) THEN
         DO i = 1, SIZE(xas_atom_env%gr)
            IF (ASSOCIATED(xas_atom_env%gr(i)%array)) THEN
               DEALLOCATE (xas_atom_env%gr(i)%array)
            END IF
         END DO
         DEALLOCATE (xas_atom_env%gr)
      END IF

      IF (ASSOCIATED(xas_atom_env%ga)) THEN
         DO i = 1, SIZE(xas_atom_env%ga)
            IF (ASSOCIATED(xas_atom_env%ga(i)%array)) THEN
               DEALLOCATE (xas_atom_env%ga(i)%array)
            END IF
         END DO
         DEALLOCATE (xas_atom_env%ga)
      END IF

      IF (ASSOCIATED(xas_atom_env%dgr1)) THEN
         DO i = 1, SIZE(xas_atom_env%dgr1)
            IF (ASSOCIATED(xas_atom_env%dgr1(i)%array)) THEN
               DEALLOCATE (xas_atom_env%dgr1(i)%array)
            END IF
         END DO
         DEALLOCATE (xas_atom_env%dgr1)
      END IF

      IF (ASSOCIATED(xas_atom_env%dgr2)) THEN
         DO i = 1, SIZE(xas_atom_env%dgr2)
            IF (ASSOCIATED(xas_atom_env%dgr2(i)%array)) THEN
               DEALLOCATE (xas_atom_env%dgr2(i)%array)
            END IF
         END DO
         DEALLOCATE (xas_atom_env%dgr2)
      END IF

      IF (ASSOCIATED(xas_atom_env%dga1)) THEN
         DO i = 1, SIZE(xas_atom_env%dga1)
            IF (ASSOCIATED(xas_atom_env%dga1(i)%array)) THEN
               DEALLOCATE (xas_atom_env%dga1(i)%array)
            END IF
         END DO
         DEALLOCATE (xas_atom_env%dga1)
      END IF

      IF (ASSOCIATED(xas_atom_env%dga2)) THEN
         DO i = 1, SIZE(xas_atom_env%dga2)
            IF (ASSOCIATED(xas_atom_env%dga2(i)%array)) THEN
               DEALLOCATE (xas_atom_env%dga2(i)%array)
            END IF
         END DO
         DEALLOCATE (xas_atom_env%dga2)
      END IF

      IF (ASSOCIATED(xas_atom_env%orb_sphi_so)) THEN
         DO i = 1, SIZE(xas_atom_env%orb_sphi_so)
            IF (ASSOCIATED(xas_atom_env%orb_sphi_so(i)%array)) THEN
               DEALLOCATE (xas_atom_env%orb_sphi_so(i)%array)
            END IF
         END DO
         DEALLOCATE (xas_atom_env%orb_sphi_so)
      END IF

      !Clean-up libint
      CALL cp_libint_static_cleanup()

      DEALLOCATE (xas_atom_env)

   END SUBROUTINE xas_atom_env_release

! **************************************************************************************************
!> \brief Releases the memory heavy attribute of xas_tdp_env that are specific to the current
!>        excited atom
!> \param xas_tdp_env ...
!> \param atom the index of the current excited atom
! **************************************************************************************************
   SUBROUTINE free_exat_memory(xas_tdp_env, atom)

      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      INTEGER, INTENT(IN)                                :: atom

      CHARACTER(len=*), PARAMETER                        :: routineN = 'free_exat_memory'

      INTEGER                                            :: i

      IF (ASSOCIATED(xas_tdp_env%ri_3c_ex)) THEN
         CALL dbcsr_t_destroy(xas_tdp_env%ri_3c_ex)
         DEALLOCATE (xas_tdp_env%ri_3c_ex)
      END IF

      IF (ASSOCIATED(xas_tdp_env%ri_fxc)) THEN
         DO i = 1, SIZE(xas_tdp_env%ri_fxc, 2)
            IF (ASSOCIATED(xas_tdp_env%ri_fxc(atom, i)%array)) THEN
               DEALLOCATE (xas_tdp_env%ri_fxc(atom, i)%array)
            END IF
         END DO
      END IF

      !because we compute the exchange integrals for each excited atom
      IF (ASSOCIATED(xas_tdp_env%opt_dist2d_ex)) THEN
         CALL distribution_2d_release(xas_tdp_env%opt_dist2d_ex)
      END IF

      xas_tdp_env%fxc_avail = .FALSE.

   END SUBROUTINE free_exat_memory

! **************************************************************************************************
!> \brief Releases a batch_info type
!> \param batch_info ...
! **************************************************************************************************
   SUBROUTINE release_batch_info(batch_info)

      TYPE(batch_info_type)                              :: batch_info

      INTEGER                                            :: i

      CALL cp_para_env_release(batch_info%para_env)

      IF (ASSOCIATED(batch_info%so_proc_info)) THEN
         DO i = 1, SIZE(batch_info%so_proc_info)
            IF (ASSOCIATED(batch_info%so_proc_info(i)%array)) THEN
               DEALLOCATE (batch_info%so_proc_info(i)%array)
            END IF
         END DO
         DEALLOCATE (batch_info%so_proc_info)
      END IF

   END SUBROUTINE release_batch_info

! **************************************************************************************************
!> \brief Uses heuristics to determine a good batching of the processros for fxc integration
!> \param batch_size ...
!> \param nbatch ...
!> \param nex_atom ...
!> \param nprocs ...
!> \note It is here and not in xas_tdp_atom because of circular dependencies issues
! **************************************************************************************************
   SUBROUTINE get_proc_batch_sizes(batch_size, nbatch, nex_atom, nprocs)

      INTEGER, INTENT(OUT)                               :: batch_size, nbatch
      INTEGER, INTENT(IN)                                :: nex_atom, nprocs

      INTEGER                                            :: rest, test_size

      !We have essentially 2 cases nex_atom >= nprocs or nex_atom < nprocs

      IF (nex_atom >= nprocs) THEN

         !If nex_atom >= nprocs, we look from batch size (starting from 1, ending with 4) that yields
         !the best indicative load balance, i.e. the best spread of excited atom per batch
         rest = 100000
         DO test_size = 1, MIN(nprocs, 4)
            nbatch = nprocs/test_size
            IF (MODULO(nex_atom, nbatch) < rest) THEN
               rest = MODULO(nex_atom, nbatch)
               batch_size = test_size
            END IF
         END DO
         nbatch = nprocs/batch_size

      ELSE

         !If nex_atom < nprocs, simply devide processors in nex_atom batches
         nbatch = nex_atom
         batch_size = nprocs/nbatch

      END IF

      !Note: because of possible odd numbers of MPI ranks / excited atoms, a couple of procs can
      !      be excluded from the batching (max 4)

   END SUBROUTINE get_proc_batch_sizes

END MODULE xas_tdp_types
