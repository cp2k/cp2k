!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Defines control structures, which contain the parameters and the
!>      settings for the DFT-based calculations.
! **************************************************************************************************
MODULE cp_control_types
   USE cp_fm_types,                     ONLY: cp_fm_p_type,&
                                              cp_fm_release
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE hirshfeld_types,                 ONLY: hirshfeld_type,&
                                              release_hirshfeld_type
   USE input_constants,                 ONLY: becke_cutoff_global,&
                                              becke_none_conf,&
                                              cdft_combined_all,&
                                              cdft_density_constraint,&
                                              do_full_density,&
                                              outer_scf_none,&
                                              radius_single
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE pw_types,                        ONLY: pw_p_type
   USE xas_control,                     ONLY: xas_control_release,&
                                              xas_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! **************************************************************************************************
! \brief Control parameters for pw grids
! **************************************************************************************************
   TYPE pw_grid_option
      LOGICAL                              :: spherical
      LOGICAL                              :: fullspace
      INTEGER, DIMENSION(2)                :: distribution_layout
      INTEGER                              :: blocked
   END TYPE pw_grid_option

! **************************************************************************************************
! \brief Control parameters for REAL_TIME_PROPAGATION calculations
! **************************************************************************************************
   TYPE rtp_control_type
      LOGICAL                              :: converged
      REAL(KIND=dp)                        :: eps_ener
      INTEGER                              :: max_iter
      INTEGER                              :: mat_exp
      INTEGER                              :: propagator
      LOGICAL                              :: fixed_ions
      INTEGER                              :: initial_wfn
      REAL(dp)                             :: eps_exp
      LOGICAL                              :: initial_step
      LOGICAL                              :: hfx_redistribute
      INTEGER                              :: aspc_order
      INTEGER                              :: sc_check_start
      LOGICAL                              :: apply_delta_pulse
      LOGICAL                              :: periodic
      LOGICAL                              :: linear_scaling
      LOGICAL                              :: write_restart
      INTEGER                              :: mcweeny_max_iter
      INTEGER                              :: acc_ref
      REAL(dp)                             :: mcweeny_eps
      INTEGER, DIMENSION(3)                :: delta_pulse_direction
      REAL(KIND=dp)                        :: delta_pulse_scale
   END TYPE rtp_control_type
! **************************************************************************************************
! \brief Control parameters for DFTB calculations
! **************************************************************************************************
   TYPE dftb_control_type
      LOGICAL                              :: self_consistent
      LOGICAL                              :: orthogonal_basis
      LOGICAL                              :: dispersion
      INTEGER                              :: dispersion_type
      LOGICAL                              :: dftb3_diagonal
      LOGICAL                              :: hb_sr_damp
      REAL(KIND=dp)                        :: hb_sr_para
      REAL(KIND=dp)                        :: eps_disp
      REAL(KIND=dp)                        :: epscn
      REAL(KIND=dp)                        :: rcdisp
      REAL(KIND=dp), DIMENSION(3)          :: sd3
      LOGICAL                              :: do_ewald
      CHARACTER(LEN=default_path_length)   :: sk_file_path
      CHARACTER(LEN=default_path_length)   :: sk_file_list
      CHARACTER(LEN=default_string_length), &
         DIMENSION(:, :), POINTER          :: sk_pair_list
      CHARACTER(LEN=default_path_length)   :: uff_force_field
      CHARACTER(LEN=default_path_length)   :: dispersion_parameter_file
   END TYPE dftb_control_type

! **************************************************************************************************
! \brief Control parameters for SCPTB calculations
! **************************************************************************************************
   TYPE scptb_control_type
      LOGICAL                              :: dispersion
      LOGICAL                              :: do_ewald
      LOGICAL                              :: do_scc
      LOGICAL                              :: do_scp
      INTEGER                              :: sto_ng
      REAL(KIND=dp)                        :: rcdisp
      REAL(KIND=dp)                        :: epscn
      REAL(KIND=dp), DIMENSION(3)          :: sd3
      CHARACTER(LEN=default_path_length)   :: parameter_file
      CHARACTER(LEN=default_path_length)   :: dispersion_parameter_file
      REAL(KIND=dp)                        :: epspair
   END TYPE scptb_control_type

! **************************************************************************************************
! \brief Control parameters for semi empirical calculations
! **************************************************************************************************
   TYPE semi_empirical_control_type
      LOGICAL                              :: orthogonal_basis
      LOGICAL                              :: analytical_gradients
      LOGICAL                              :: force_kdsod_EX
      LOGICAL                              :: do_ewald, do_ewald_r3, do_ewald_gks
      INTEGER                              :: integral_screening, periodic_type
      INTEGER                              :: max_multipole
      INTEGER                              :: ga_ncells
      REAL(KIND=dp)                        :: delta
      ! Dispersion pair potential
      LOGICAL                              :: dispersion
      REAL(KIND=dp)                        :: rcdisp
      REAL(KIND=dp)                        :: epscn
      REAL(KIND=dp), DIMENSION(3)          :: sd3
      CHARACTER(LEN=default_path_length)   :: dispersion_parameter_file
      ! Parameters controlling the evaluation of the integrals
      REAL(KIND=dp)                        :: cutoff_lrc, taper_lrc, range_lrc
      REAL(KIND=dp)                        :: cutoff_cou, taper_cou, range_cou
      REAL(KIND=dp)                        :: cutoff_exc, taper_exc, range_exc
      REAL(KIND=dp)                        :: taper_scr, range_scr
   END TYPE semi_empirical_control_type

! **************************************************************************************************
! \brief Control parameters for GAPW method within QUICKSTEP ***
! **************************************************************************************************
   TYPE gapw_control_type
      REAL(KIND=dp)                        :: eps_fit, &
                                              eps_iso, &
                                              eps_Vrho0, &
                                              eps_svd, &
                                              eps_cpc
      INTEGER                              :: ladd_rho0, &
                                              lmax_rho0, &
                                              lmax_sphere, &
                                              quadrature
      LOGICAL                              :: lrho1_eq_lrho0
      LOGICAL                              :: alpha0_hard_from_input, &
                                              force_paw, &
                                              non_paw_atoms, &
                                              nopaw_as_gpw
      REAL(KIND=dp)                        :: alpha0_hard
      REAL(KIND=dp)                        :: max_rad_local
   END TYPE gapw_control_type
! **************************************************************************************************
! \brief parameters for calculations involving a time dependent electric field
! **************************************************************************************************
   TYPE efield_type
      REAL(KIND=dp)                        :: actual_time
      REAL(KIND=dp), DIMENSION(:), POINTER :: polarisation
      INTEGER                              :: envelop_id
      REAL(KIND=dp), DIMENSION(:), POINTER :: envelop_r_vars
      INTEGER, DIMENSION(:), POINTER       :: envelop_i_vars
      REAL(KIND=dp)                        :: strength
      REAL(KIND=dp)                        :: phase_offset
      REAL(KIND=dp)                        :: wavelength
   END TYPE efield_type

   TYPE efield_p_type
      TYPE(efield_type), POINTER           :: efield
   END TYPE efield_p_type
! **************************************************************************************************
! \brief parameters for calculations involving a time dependent electric field
! **************************************************************************************************
   TYPE period_efield_type
      LOGICAL                              :: displacement_field
      REAL(KIND=dp), DIMENSION(3)          :: polarisation
      REAL(KIND=dp), DIMENSION(3)          :: d_filter
      REAL(KIND=dp)                        :: strength
   END TYPE period_efield_type

! **************************************************************************************************
! \brief some parameters useful for mulliken_restraints
! **************************************************************************************************
   TYPE mulliken_restraint_type
      INTEGER                              :: ref_count
      REAL(KIND=dp)                        :: strength
      REAL(KIND=dp)                        :: TARGET
      INTEGER                              :: natoms
      INTEGER, POINTER, DIMENSION(:)       :: atoms
   END TYPE mulliken_restraint_type

! **************************************************************************************************
! \brief some parameters useful for ddapc_restraints
! **************************************************************************************************
   TYPE ddapc_restraint_type
      INTEGER                              :: ref_count
      REAL(KIND=dp)                        :: strength
      REAL(KIND=dp)                        :: TARGET
      REAL(KIND=dp)                        :: ddapc_order_p
      INTEGER                              :: functional_form
      INTEGER                              :: natoms
      INTEGER, POINTER, DIMENSION(:)       :: atoms
      REAL(KIND=dp), POINTER, DIMENSION(:) :: coeff
      INTEGER                              :: density_type
   END TYPE ddapc_restraint_type

! **************************************************************************************************
! \brief provides a vector of pointers to ddapc_restraint_type
! **************************************************************************************************
   TYPE ddapc_restraint_p_type
      TYPE(ddapc_restraint_type), POINTER  :: ddapc_restraint_control
   END TYPE ddapc_restraint_p_type

! **************************************************************************************************
!> \brief some parameters useful for becke_restraints
!> \param fragment_a_fname filename of cube file holding the isolated electron density of fragment a
!> \param fragment_b_fname filename of cube file holding the isolated electron density of fragment b
!> \param ref_count the ref count
!> \param rglobal global cutoff to use for building the constraint
!> \param target target values of the constraints
!> \param strength Lagrangian multipliers of the constraints
!> \param becke_order_p integrated values of the constraints
!> \param natoms number of constraint atoms
!> \param confine_method the confine method id
!> \param cavity_shape the confinement cavity shape id
!> \param constraint_type the constraint type to use
!> \param cutoff_type the cutoff type to use for building the constraint
!> \param combined_type for combined density+spin constraint, determines which atoms to apply the spin
!>                      constraint to
!> \param atoms list of constraint atoms
!> \param becke_pot the Becke real space potential
!> \param cavity the Gaussian confinement cavity: the constraint is nonzero outside this cavity
!> \param combined_weight the Becke real space potential for the spin constraint when using combined
!>                        density+spin constraint
!> \param charge atomic Becke real space potentials needed to calculate atomic Becke charges
!> \param fragments container for isolated fragment densities read from cube files
!> \param need_pot logical which determines if the Becke potential needs to be built
!> \param save_pot logical which determines if the Becke potential should be saved until forces
!>                 have been evaluated
!> \param in_memory logical which determines if the gradients of the Becke potential should be
!>                  computed simultaneously with the potential instead of separately
!> \param adjust logical which determines if the Becke potential is adjusted with atomic radii
!> \param atomic_charges logical which determines if atomic Becke charges should be computed
!> \param confine logical which determines if static confinement is active
!> \param dynamic_confine logical which determines if dynamic confinement is active
!> \param cavity_confine logical which determines if cavity confinement is active
!> \param should_skip logical which determines is grid points should be skipped if all constraint
!>                    atoms are found to reside beyond the cutoff distance from it
!> \param print_cavity logical to print the Gaussian confinement cavity
!> \param external_control logical which determines if the constraint has already been built
!>                         in a mixed_env that holds multiple CDFT states
!> \param first_iteration a flag to mark the first iteration for printing of additional data
!> \param fragment_density use isolated fragment densities as a reference for the constraint
!> \param fragments_integrated logical to determine if the fragment densities have been integrated
!> \param use_bohr decides whether to use angstrom or bohr units for the confinement cavity radius
!> \param coeff determines how to sum up the constraint atoms to form the value of the constraint
!> \param cutoffs element specific cutoffs
!> \param cutoffs_tmp same as cutoffs but a temporary read during parsing of this type
!> \param charges_fragment the values of the integrated isolated fragment densities
!> \param radii_tmp temporary list of element specific atomic radii used to adjust the Becke cells
!> \param radii permanent copy of radii_tmp
!> \param confine_dir for dynamic/mixed/static confinement determines the direction of confinement
!> \param confine_bounds_int grid point indices outside which the constraint vanishes along confine_dir
!> \param confine_bounds same as confine_bounds_int, but in distance units
!> \param dynamic_radius the radius parameter used for creating the dynamic confinement bounds
!> \param rcavity an optional global radius parameter used to define the Gaussian confinement cavity
!> \param eps_cavity threshold used screen small values of the Gaussian cavity density
!> \param cavity_env the structure used to build the Gaussian cavity
!> \param aij pairwise parameters used to adjust the Becke cell boundaries built from atomic radii
!> \param skip_list a list of grid points skipped my keyword should_skip
!> \param cavity_mat a compacted version of cavity
!> \param combined_mat a compated version of combined_weight
!> \param gradients atomic gradients of the weight function
! **************************************************************************************************
   TYPE becke_restraint_type
      CHARACTER(LEN=default_path_length)   :: fragment_a_fname, &
                                              fragment_b_fname
      INTEGER                              :: ref_count
      REAL(KIND=dp)                        :: rglobal
      REAL(KIND=dp), DIMENSION(:), POINTER :: TARGET, strength, &
                                              becke_order_p
      INTEGER                              :: natoms
      INTEGER                              :: confine_method, cavity_shape, cutoff_type, &
                                              constraint_type, combined_type
      INTEGER, POINTER, DIMENSION(:)       :: atoms
      TYPE(pw_p_type)                      :: becke_pot, cavity, combined_weight
      TYPE(pw_p_type), POINTER, &
         DIMENSION(:)                      :: charge, fragments
      LOGICAL                              :: need_pot, save_pot, in_memory, &
                                              adjust, atomic_charges, confine, &
                                              dynamic_confine, cavity_confine, &
                                              should_skip, print_cavity, &
                                              external_control, first_iteration, &
                                              fragment_density, fragments_integrated, &
                                              use_bohr, transfer_pot
      REAL(KIND=dp), POINTER, DIMENSION(:) :: coeff, cutoffs, cutoffs_tmp, &
                                              charges_fragment, radii_tmp, &
                                              radii
      INTEGER                              :: confine_dir, confine_bounds_int(2)
      REAL(KIND=dp)                        :: confine_bounds(2), dynamic_radius, &
                                              rcavity, eps_cavity
      TYPE(hirshfeld_type), POINTER        :: cavity_env
      REAL(KIND=dp), POINTER, &
         DIMENSION(:, :)                   :: aij
      LOGICAL, POINTER, DIMENSION(:, :, :) :: skip_list
      REAL(KIND=dp), POINTER, &
         DIMENSION(:, :, :)                :: cavity_mat, combined_mat
      REAL(KIND=dp), POINTER, &
         DIMENSION(:, :, :, :)             :: gradients
   END TYPE becke_restraint_type

! **************************************************************************************************
! \brief control parameters for hirshfeld constraint
!> \param hirshfeld_env holds the information needed to build the Gaussian Hirshfeld weight
!> \param natoms the number of constraint atoms
!> \param constraint_type the constraint type to use
!> \param combined_type for combined density+spin constraint, determines which atoms to apply the spin
!>                      constraint to
!> \param atoms list of constraint atoms
!> \param coeff determines how to sum up the constraint atoms to form the value of the constraint
! **************************************************************************************************
   TYPE hirshfeld_constraint_type
      TYPE(hirshfeld_type), POINTER        :: hirshfeld_env
      INTEGER                              :: natoms, constraint_type, &
                                              combined_type
      INTEGER, POINTER, DIMENSION(:)       :: atoms
      REAL(KIND=dp), POINTER, DIMENSION(:) :: coeff
   END TYPE hirshfeld_constraint_type

! **************************************************************************************************
! \brief control parameters for CDFT with OT
!        TODO: make cdft_control parent section for Becke restraint type
!> \param ref_count the ref count
!> \param total_steps counter to keep track of the total number of SCF steps
!> \param type type of CDFT constraint
!> \param precond_freq preconditioner can be used if SCF converged in less than precond_freq steps
!> \param nreused determines how many times the current OT preconditioner has been reused
!> \param max_reuse the same preconditioner can be used a maximum of max_reuse times
!> \param purge_freq determines how large nbad_conv can grow before purging the wfn/constraint history
!> \param nbad_conv a running counter keeping track of the number of CDFT SCF loops when the first
!>                  CDFT SCF iteration required more than 1 outer SCF loop. Reset when convergence is
!>                  smooth
!> \param purge_offset purging is only allowed when more than purge_offset steps have passed since
!>                     last purge
!> \param istep a counter to keep track of how many steps have passed since the last purge
!> \param ienergy a counter tracking the total number of CDFT energy evaluations
!> \param constraint_type the constraint type to use
!> \param combined_type for combined density+spin constraint, determines which atoms to apply the spin
!>                      constraint to
!> \param constraints holds information about the CDFT SCF loop
!> \param constraint_control the outer_scf_control_type for the constraint
!> \param ot_control the outer_scf_control_type for OT where data is stashed when outside the OT
!>                   outer loop
!> \param hirshfeld_control control parameters for a Hirshfeld constraint
!> \param strength Lagrangian multipliers of the constraints
!> \param target target values of the constraints
!> \param value integrated values of the constraints
!> \param weight the constraint potential
!> \param need_pot logical which determines if the constraint potential needs to be built
!> \param save_pot logical which determines if the constraint potential should be saved until forces
!>                 have been evaluated
!> \param do_et logical which determines if a ET coupling calculation was requested
!> \param reuse_precond logical which determines if a preconditioner can be reused
!> \param purge_history logical which determines if the wfn/constraint history can be purged
!> \param should_purge logical which determines if purging should take place after this CDFT SCF loop
!> \param calculate_metric logical which determines if the ET coupling reliablity metric is computed
!> \param mo_coeff save the MO coeffs (for do_et)
!> \param wmat matrix representation of the weight function (for do_et)
!> \param matrix_s save the overlap matrix (for do_et)
!> \param matrix_p save the density matrix (for calculate_metric)
! **************************************************************************************************
   ! Copied from qs_scf_types to avoid circular dependency
   TYPE qs_outer_scf_type
      INTEGER                              :: iter_count
      LOGICAL                              :: deallocate_jacobian
      REAL(KIND=dp), DIMENSION(:), POINTER :: energy
      REAL(KIND=dp), DIMENSION(:, :), &
         POINTER                           :: variables
      REAL(KIND=dp), DIMENSION(:, :), &
         POINTER                           :: gradient
      REAL(KIND=dp), DIMENSION(:, :), &
         POINTER                           :: inv_jacobian
      INTEGER, DIMENSION(:), POINTER       :: count
   END TYPE qs_outer_scf_type

   ! Copied from scf_control_types to avoid circular dependency
   TYPE outer_scf_control_type
      LOGICAL                              :: have_scf
      LOGICAL                              :: build_jacobian
      LOGICAL                              :: broyden_update
      INTEGER                              :: max_scf, &
                                              ijacobian
      REAL(KIND=dp)                        :: eps_scf, step_size, &
                                              jacobian_step, &
                                              newton_step
      INTEGER                              :: TYPE
      INTEGER                              :: optimizer
      INTEGER                              :: diis_buffer_length
      INTEGER                              :: extrapolation_order
      INTEGER                              :: bisect_trust_count
      INTEGER                              :: jacobian_type
      INTEGER                              :: broyden_type
      INTEGER                              :: jacobian_freq(2)
   END TYPE outer_scf_control_type

   TYPE cdft_control_type
      INTEGER                              :: ref_count, total_steps, TYPE, &
                                              precond_freq, nreused, max_reuse, &
                                              purge_freq, nbad_conv, purge_offset, &
                                              istep, constraint_type, combined_type, &
                                              ienergy
      TYPE(qs_outer_scf_type)              :: constraint
      TYPE(outer_scf_control_type)         :: constraint_control, ot_control
      TYPE(hirshfeld_constraint_type), &
         POINTER                           :: hirshfeld_control
      REAL(KIND=dp), DIMENSION(:), POINTER :: strength, TARGET, value
      TYPE(pw_p_type)                      :: weight
      LOGICAL                              :: need_pot, save_pot, do_et, &
                                              reuse_precond, purge_history, &
                                              should_purge, calculate_metric
      TYPE(cp_fm_p_type), DIMENSION(:), &
         POINTER                           :: mo_coeff
      TYPE(dbcsr_p_type)                   :: wmat, matrix_s
      TYPE(dbcsr_p_type), DIMENSION(:), &
         POINTER                           :: matrix_p
   END TYPE cdft_control_type

! **************************************************************************************************
! \brief some parameters useful for s2_restraints
! **************************************************************************************************
   TYPE s2_restraint_type
      INTEGER                              :: ref_count
      REAL(KIND=dp)                        :: strength
      REAL(KIND=dp)                        :: TARGET
      REAL(KIND=dp)                        :: s2_order_p
      INTEGER                              :: functional_form
   END TYPE s2_restraint_type

! **************************************************************************************************
! \brief some parameters useful for auxiliary density matrix method
! **************************************************************************************************
   TYPE admm_block_type
      INTEGER, DIMENSION(:), ALLOCATABLE   :: list
   END TYPE admm_block_type

   TYPE admm_control_type
      REAL(KIND=dp)                        :: eps_filter
      INTEGER                              :: purification_method
      INTEGER                              :: method
      LOGICAL                              :: charge_constrain
      INTEGER                              :: scaling_model
      INTEGER                              :: aux_exch_func
      LOGICAL                              :: aux_exch_func_param
      REAL(KIND=dp), DIMENSION(3)          :: aux_x_param
      TYPE(admm_block_type), DIMENSION(:), &
         ALLOCATABLE                       :: blocks
   END TYPE admm_control_type

! **************************************************************************************************
! \brief Control parameters for a QUICKSTEP and KIM-GORDON calculation ***
!        eps_pgf_orb: Cutoff value for the interaction of the primitive
!                     Gaussian-type functions (primitive basis functions).
! **************************************************************************************************
   TYPE qs_control_type
      CHARACTER(LEN=10)                    :: method
      INTEGER                              :: method_id
      REAL(KIND=dp)                        :: eps_core_charge, &
                                              eps_kg_orb, &
                                              eps_pgf_orb, &
                                              eps_ppl, &
                                              eps_ppnl, &
                                              eps_rho_gspace, &
                                              eps_rho_rspace, &
                                              eps_filter_matrix, &
                                              eps_gvg_rspace, &
                                              progression_factor, &
                                              relative_cutoff
      LOGICAL                              :: do_almo_scf
      LOGICAL                              :: do_ls_scf
      LOGICAL                              :: do_kg
      LOGICAL                              :: commensurate_mgrids
      LOGICAL                              :: realspace_mgrids
      LOGICAL                              :: map_consistent
      LOGICAL                              :: gapw, gapw_xc, gpw, pao
      LOGICAL                              :: lrigpw
      LOGICAL                              :: lri_optbas
      LOGICAL                              :: ofgpw
      LOGICAL                              :: dftb
      LOGICAL                              :: scptb
      LOGICAL                              :: semi_empirical
      LOGICAL                              :: mulliken_restraint
      LOGICAL                              :: ddapc_restraint
      LOGICAL                              :: ddapc_restraint_is_spin
      LOGICAL                              :: ddapc_explicit_potential
      LOGICAL                              :: becke_restraint
      LOGICAL                              :: cdft
      LOGICAL                              :: et_coupling_calc
      LOGICAL                              :: s2_restraint
      INTEGER                              :: do_ppl_method
      INTEGER                              :: wf_interpolation_method_nr
      INTEGER                              :: wf_extrapolation_order
      REAL(KIND=dp)                        :: cutoff
      REAL(KIND=dp), DIMENSION(:), POINTER :: e_cutoff
      TYPE(mulliken_restraint_type), &
         POINTER                           :: mulliken_restraint_control
      TYPE(ddapc_restraint_p_type), &
         DIMENSION(:), POINTER             :: ddapc_restraint_control
      TYPE(becke_restraint_type), POINTER  :: becke_control
      TYPE(cdft_control_type), POINTER     :: cdft_control
      TYPE(s2_restraint_type), POINTER     :: s2_restraint_control
      TYPE(dftb_control_type), POINTER     :: dftb_control
      TYPE(scptb_control_type), POINTER    :: scptb_control
      TYPE(semi_empirical_control_type), &
         POINTER                           :: se_control
      TYPE(gapw_control_type), POINTER     :: gapw_control
      TYPE(pw_grid_option)                 :: pw_grid_opt
      LOGICAL                              :: skip_load_balance_distributed
   END TYPE qs_control_type

! **************************************************************************************************
! \brief Control parameters for the SCCS models
! **************************************************************************************************
   TYPE sccs_control_type
      LOGICAL                              :: sccs_activated
      INTEGER                              :: derivative_method, &
                                              max_iter, &
                                              method_id, &
                                              ref_count
      REAL(KIND=dp)                        :: alpha_solvent, &
                                              beta, &
                                              beta_solvent, &
                                              delta_rho, &
                                              eps_sccs, &
                                              eps_scf, &
                                              epsilon_solvent, &
                                              gamma_solvent, &
                                              mixing, &
                                              rho_zero, &
                                              rho_max, &
                                              rho_min
   END TYPE sccs_control_type

! **************************************************************************************************
! \brief Control parameters for a TIME-DEPENDENT PERTURBATION calculation
! \par  ATTRIBUTES
!   - n_ev       : number of eigenvalues to calculate
!   - n_reortho  : how many time to reorthogonalize (in the lanczos algorithm)
!   - do_kernel  : wether to evaluate the kernel (this is a debugging option)
!   - res_etype : { SINGLET | TRIPLET } which excitations
!                  to calculate
!   - lumos_eigenvalues : holds the eigenvalues of the lumos (if calculated in QS)
!
! \par  NOTES
!   The lumos are helpfull in choosing a initial vector for the TDDFPT
!   calculation, since they can be used to construct the solutions of the
!   TDDFPT operator without the perturbation kernel.
! **************************************************************************************************
   TYPE tddfpt_control_type
      TYPE(cp_fm_p_type), DIMENSION(:), &
         POINTER                           :: lumos
      REAL(KIND=dp)                        :: tolerance
      INTEGER                              :: n_ev
      INTEGER                              :: max_kv
      INTEGER                              :: n_restarts
      INTEGER                              :: n_reortho
      LOGICAL                              :: do_kernel
      LOGICAL                              :: lsd_singlets
      LOGICAL                              :: invert_S
      LOGICAL                              :: precond
      LOGICAL                              :: drho_by_collocation
      LOGICAL                              :: use_kinetic_energy_density
      INTEGER                              :: res_etype
      INTEGER                              :: diag_method
      INTEGER                              :: oe_corr
      INTEGER                              :: sic_method_id
      INTEGER                              :: sic_list_id
      REAL(KIND=dp)                        :: sic_scaling_a, sic_scaling_b
      REAL(KIND=dp), DIMENSION(:, :), &
         POINTER                           :: lumos_eigenvalues
   END TYPE tddfpt_control_type

! **************************************************************************************************
! \brief Control parameters for a Time-Dependent DFT calculation
! **************************************************************************************************
   TYPE tddfpt2_control_type
      !> compute TDDFPT excitation energies and oscillator strengths
      LOGICAL                              :: enabled
      !> number of excited states to converge
      INTEGER                              :: nstates
      !> number of additional excited states in the active space
      INTEGER                              :: added_states
      !> maximal number of iterations to be performed
      INTEGER                              :: niters
      !> maximal number of Krylov space vectors
      INTEGER                              :: nkvs
      !> minimal number of MPI processes to be used per excited state
      INTEGER                              :: nprocs
      !> maximal number of orthogonalisation iterations
      INTEGER                              :: nreortho
      !> target accuracy
      REAL(kind=dp)                        :: conv
      !> energy threshold which controls when excited states are considered to be degenerate
      REAL(kind=dp)                        :: degenerate_eps
      !> the smallest excitation amplitude to print
      REAL(kind=dp)                        :: min_excitation_amplitude
      !> threshold which controls when two wave functions considered to be orthogonal:
      !> maxabs(Ci^T * S * Cj) <= orthogonal_eps
      REAL(kind=dp)                        :: orthogonal_eps
      !> should all degenerate excited states be automatically added to the active space
      LOGICAL                              :: add_degenerate
      !> read guess wave functions from restart file if exists
      LOGICAL                              :: is_restart
      !> compute triplet excited states using spin-unpolarised molecular orbitals
      LOGICAL                              :: rks_triplets
      !
      ! DIPOLE_MOMENTS subsection
      !
      ! form of the dipole operator used to compute oscillator strengths
      INTEGER                              :: dipole_form
      !> type of the reference point used for calculation of electrostatic dipole moments
      INTEGER                              :: dipole_reference
      !> user-defined reference point
      REAL(kind=dp), DIMENSION(:), POINTER :: dipole_ref_point
      !
      ! MGRID subsection
      !
      !> number of plain-wave grids
      INTEGER                              :: mgrid_ngrids
      !> create commensurate grids (progression factor and cutoff values of sub-grids will be ignored)
      LOGICAL                              :: mgrid_commensurate_mgrids
      !> signals that MGRID section has been explicitly given. Other mgrid_* variables
      !> are not initialised when it is equal to .FALSE. as in this case the default
      !> set of plain-wave grids will be used
      LOGICAL                              :: mgrid_is_explicit
      !> same as qs_control%realspace_mgrids
      LOGICAL                              :: mgrid_realspace_mgrids
      !> do not perform load balancing
      LOGICAL                              :: mgrid_skip_load_balance
      !> cutoff value at the finest grid level
      REAL(kind=dp)                        :: mgrid_cutoff
      !> cutoff at the next grid level will be smaller then the cutoff
      !> at the current grid by this number of times
      REAL(kind=dp)                        :: mgrid_progression_factor
      !> cutoff that determines to which grid a particular Gaussian function will be mapped
      REAL(kind=dp)                        :: mgrid_relative_cutoff
      !> manually provided the list of cutoff values for each grid level
      !> (when it is null(), the cutoff values will be assigned automatically)
      REAL(kind=dp), DIMENSION(:), POINTER :: mgrid_e_cutoff
   END TYPE tddfpt2_control_type

! **************************************************************************************************
! \brief Control parameters for a DFT calculation
! **************************************************************************************************
   TYPE dft_control_type
      TYPE(admm_control_type), POINTER     :: admm_control
      TYPE(period_efield_type), POINTER    :: period_efield
      TYPE(qs_control_type), POINTER       :: qs_control
      TYPE(rtp_control_type), POINTER      :: rtp_control
      TYPE(sccs_control_type), POINTER     :: sccs_control
      TYPE(tddfpt_control_type), POINTER   :: tddfpt_control
      TYPE(tddfpt2_control_type), POINTER  :: tddfpt2_control
      TYPE(xas_control_type), POINTER      :: xas_control
      TYPE(efield_p_type), POINTER, &
         DIMENSION(:)                      :: efield_fields
      INTEGER                              :: nspins, &
                                              charge, &
                                              multiplicity, &
                                              sic_method_id, &
                                              ref_count, &
                                              id_nr, &
                                              plus_u_method_id, &
                                              dir_surf_dip, &
                                              nimages = 1
      INTEGER                              :: sic_list_id
      REAL(KIND=dp)                        :: relax_multiplicity, &
                                              sic_scaling_a, &
                                              sic_scaling_b
      LOGICAL                              :: do_tddfpt_calculation, &
                                              do_xas_calculation, &
                                              drho_by_collocation, &
                                              use_kinetic_energy_density, &
                                              restricted, &
                                              roks, &
                                              uks, &
                                              lsd, &
                                              dft_plus_u, &
                                              apply_efield, &
                                              apply_efield_field, &
                                              apply_period_efield, &
                                              apply_external_potential, &
                                              eval_external_potential, &
                                              do_admm, &
                                              do_admm_dm, &
                                              do_admm_mo, &
                                              smear, &
                                              low_spin_roks, &
                                              apply_external_density, &
                                              read_external_density, &
                                              apply_external_vxc, &
                                              read_external_vxc, &
                                              correct_surf_dip, &
                                              do_sccs
   END TYPE dft_control_type

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_control_types'
   INTEGER, SAVE :: last_dft_control_id = 0

   ! Public data types

   PUBLIC :: dft_control_type, &
             qs_control_type, &
             gapw_control_type, &
             tddfpt_control_type, &
             tddfpt2_control_type, &
             efield_type, &
             mulliken_restraint_type, &
             ddapc_restraint_type, &
             dftb_control_type, &
             scptb_control_type, &
             becke_restraint_type, &
             semi_empirical_control_type, &
             s2_restraint_type, &
             admm_control_type, &
             rtp_control_type, &
             sccs_control_type, &
             cdft_control_type

   ! Public subroutines

   PUBLIC :: dft_control_retain, &
             dft_control_release, &
             dft_control_create, &
             tddfpt_control_create, &
             admm_control_create, &
             ddapc_control_create, &
             sccs_control_create, &
             becke_control_create, &
             becke_control_release

CONTAINS

! **************************************************************************************************
!> \brief create  the mulliken_restraint_type
!> \param mulliken_restraint_control ...
!> \par History
!>      02.2005 created [Joost VandeVondele]
! **************************************************************************************************
   SUBROUTINE mulliken_control_create(mulliken_restraint_control)
      TYPE(mulliken_restraint_type), POINTER             :: mulliken_restraint_control

      CHARACTER(len=*), PARAMETER :: routineN = 'mulliken_control_create', &
         routineP = moduleN//':'//routineN

      CPASSERT(.NOT. ASSOCIATED(mulliken_restraint_control))
      ALLOCATE (mulliken_restraint_control)

      mulliken_restraint_control%ref_count = 1
      mulliken_restraint_control%strength = 0.1_dp
      mulliken_restraint_control%target = 1.0_dp
      mulliken_restraint_control%natoms = 0
      NULLIFY (mulliken_restraint_control%atoms)
   END SUBROUTINE mulliken_control_create

! **************************************************************************************************
!> \brief release the mulliken_restraint_type
!> \param mulliken_restraint_control ...
!> \par History
!>      02.2005 created [Joost VandeVondele]
! **************************************************************************************************
   SUBROUTINE mulliken_control_release(mulliken_restraint_control)
      TYPE(mulliken_restraint_type), POINTER             :: mulliken_restraint_control

      CHARACTER(len=*), PARAMETER :: routineN = 'mulliken_control_release', &
         routineP = moduleN//':'//routineN

      CPASSERT(ASSOCIATED(mulliken_restraint_control))
      CPASSERT(mulliken_restraint_control%ref_count > 0)
      mulliken_restraint_control%ref_count = mulliken_restraint_control%ref_count-1
      IF (mulliken_restraint_control%ref_count == 0) THEN
         IF (ASSOCIATED(mulliken_restraint_control%atoms)) &
            DEALLOCATE (mulliken_restraint_control%atoms)
         mulliken_restraint_control%ref_count = 0
         mulliken_restraint_control%strength = 0.0_dp
         mulliken_restraint_control%target = 0.0_dp
         mulliken_restraint_control%natoms = 0
         DEALLOCATE (mulliken_restraint_control)
      ENDIF
   END SUBROUTINE mulliken_control_release

! **************************************************************************************************
!> \brief retain the mulliken_restraint_type
!> \param mulliken_restraint_control ...
!> \par History
!>      02.2005 created [Joost VandeVondele]
! **************************************************************************************************
   SUBROUTINE mulliken_control_retain(mulliken_restraint_control)
      TYPE(mulliken_restraint_type), POINTER             :: mulliken_restraint_control

      CHARACTER(len=*), PARAMETER :: routineN = 'mulliken_control_retain', &
         routineP = moduleN//':'//routineN

      CPASSERT(ASSOCIATED(mulliken_restraint_control))

      mulliken_restraint_control%ref_count = mulliken_restraint_control%ref_count+1
   END SUBROUTINE mulliken_control_retain

! **************************************************************************************************
!> \brief create the ddapc_restraint_type
!> \param ddapc_restraint_control ...
!> \par History
!>      02.2006 created [Joost VandeVondele]
! **************************************************************************************************
   SUBROUTINE ddapc_control_create(ddapc_restraint_control)
      TYPE(ddapc_restraint_type), POINTER                :: ddapc_restraint_control

      CHARACTER(len=*), PARAMETER :: routineN = 'ddapc_control_create', &
         routineP = moduleN//':'//routineN

      CPASSERT(.NOT. ASSOCIATED(ddapc_restraint_control))
      ALLOCATE (ddapc_restraint_control)

      ddapc_restraint_control%density_type = do_full_density
      ddapc_restraint_control%ref_count = 1
      ddapc_restraint_control%strength = 0.1_dp
      ddapc_restraint_control%ddapc_order_p = 0.0_dp
      ddapc_restraint_control%functional_form = -1
      ddapc_restraint_control%target = 1.0_dp
      ddapc_restraint_control%natoms = 0
      NULLIFY (ddapc_restraint_control%atoms)
      NULLIFY (ddapc_restraint_control%coeff)

   END SUBROUTINE ddapc_control_create

! **************************************************************************************************
!> \brief release the ddapc_restraint_type
!> \param ddapc_restraint_control ...
!> \par History
!>      02.2006 created [Joost VandeVondele]
! **************************************************************************************************
   SUBROUTINE ddapc_control_release(ddapc_restraint_control)
      TYPE(ddapc_restraint_type), POINTER                :: ddapc_restraint_control

      CHARACTER(len=*), PARAMETER :: routineN = 'ddapc_control_release', &
         routineP = moduleN//':'//routineN

      CPASSERT(ASSOCIATED(ddapc_restraint_control))
      CPASSERT(ddapc_restraint_control%ref_count > 0)
      ddapc_restraint_control%ref_count = ddapc_restraint_control%ref_count-1
      IF (ddapc_restraint_control%ref_count == 0) THEN
         IF (ASSOCIATED(ddapc_restraint_control%atoms)) &
            DEALLOCATE (ddapc_restraint_control%atoms)
         IF (ASSOCIATED(ddapc_restraint_control%coeff)) &
            DEALLOCATE (ddapc_restraint_control%coeff)
         ddapc_restraint_control%ref_count = 0
         ddapc_restraint_control%strength = 0.0_dp
         ddapc_restraint_control%target = 0.0_dp
         ddapc_restraint_control%natoms = 0
         DEALLOCATE (ddapc_restraint_control)
      ENDIF
   END SUBROUTINE ddapc_control_release

! **************************************************************************************************
!> \brief retain the ddapc_restraint_type
!> \param ddapc_restraint_control ...
!> \par History
!>      02.2006 created [Joost VandeVondele]
! **************************************************************************************************
   SUBROUTINE ddapc_control_retain(ddapc_restraint_control)
      TYPE(ddapc_restraint_type), POINTER                :: ddapc_restraint_control

      CHARACTER(len=*), PARAMETER :: routineN = 'ddapc_control_retain', &
         routineP = moduleN//':'//routineN

      CPASSERT(ASSOCIATED(ddapc_restraint_control))

      ddapc_restraint_control%ref_count = ddapc_restraint_control%ref_count+1
   END SUBROUTINE ddapc_control_retain

! **************************************************************************************************
!> \brief create the becke_restraint_type
!> \param becke_control the structure to create
!> \par History
!>      02.2007 created [Florian Schiffmann]
! **************************************************************************************************
   SUBROUTINE becke_control_create(becke_control)
      TYPE(becke_restraint_type), POINTER                :: becke_control

      CHARACTER(len=*), PARAMETER :: routineN = 'becke_control_create', &
         routineP = moduleN//':'//routineN

      CPASSERT(.NOT. ASSOCIATED(becke_control))
      ALLOCATE (becke_control)

      becke_control%ref_count = 1
      becke_control%natoms = 0
      becke_control%need_pot = .TRUE.
      becke_control%save_pot = .FALSE.
      becke_control%transfer_pot = .FALSE.
      becke_control%adjust = .FALSE.
      becke_control%atomic_charges = .FALSE.
      becke_control%confine_method = becke_none_conf
      becke_control%cutoff_type = becke_cutoff_global
      becke_control%confine = .FALSE.
      becke_control%dynamic_confine = .FALSE.
      becke_control%cavity_confine = .FALSE.
      becke_control%should_skip = .FALSE.
      becke_control%print_cavity = .FALSE.
      becke_control%in_memory = .FALSE.
      becke_control%first_iteration = .TRUE.
      becke_control%fragment_density = .FALSE.
      becke_control%fragments_integrated = .FALSE.
      becke_control%use_bohr = .FALSE.
      becke_control%confine_dir = 3
      becke_control%confine_bounds = 0.0_dp
      becke_control%confine_bounds_int = 0
      becke_control%dynamic_radius = 6.0_dp
      becke_control%rcavity = 3.0_dp
      becke_control%rglobal = 6.0_dp
      becke_control%eps_cavity = 1.0e-5_dp
      becke_control%cavity_shape = radius_single
      becke_control%external_control = .FALSE.
      becke_control%constraint_type = cdft_density_constraint
      becke_control%combined_type = cdft_combined_all
      NULLIFY (becke_control%becke_pot%pw)
      NULLIFY (becke_control%cavity%pw)
      NULLIFY (becke_control%combined_weight%pw)
      NULLIFY (becke_control%strength)
      NULLIFY (becke_control%target)
      NULLIFY (becke_control%becke_order_p)
      NULLIFY (becke_control%atoms)
      NULLIFY (becke_control%coeff)
      NULLIFY (becke_control%charge)
      NULLIFY (becke_control%aij)
      NULLIFY (becke_control%gradients)
      NULLIFY (becke_control%cavity_mat)
      NULLIFY (becke_control%combined_mat)
      NULLIFY (becke_control%cavity_env)
      NULLIFY (becke_control%skip_list)
      NULLIFY (becke_control%cutoffs)
      NULLIFY (becke_control%cutoffs_tmp)
      NULLIFY (becke_control%charges_fragment)
      NULLIFY (becke_control%fragments)
      NULLIFY (becke_control%radii)
      NULLIFY (becke_control%radii_tmp)
   END SUBROUTINE becke_control_create

! **************************************************************************************************
!> \brief release the becke_restraint_type
!> \param becke_control the structure to release
!> \par History
!>      02.2007 created [Florian Schiffmann]
! **************************************************************************************************
   SUBROUTINE becke_control_release(becke_control)
      TYPE(becke_restraint_type), POINTER                :: becke_control

      CHARACTER(len=*), PARAMETER :: routineN = 'becke_control_release', &
         routineP = moduleN//':'//routineN

      CPASSERT(ASSOCIATED(becke_control))
      CPASSERT(becke_control%ref_count > 0)
      becke_control%ref_count = becke_control%ref_count-1
      IF (becke_control%ref_count == 0) THEN
         IF (ASSOCIATED(becke_control%atoms)) &
            DEALLOCATE (becke_control%atoms)
         IF (ASSOCIATED(becke_control%cutoffs)) &
            DEALLOCATE (becke_control%cutoffs)
         IF (ASSOCIATED(becke_control%cutoffs_tmp)) &
            DEALLOCATE (becke_control%cutoffs_tmp)
         IF (ASSOCIATED(becke_control%radii_tmp)) &
            DEALLOCATE (becke_control%radii_tmp)
         IF (ASSOCIATED(becke_control%radii)) &
            DEALLOCATE (becke_control%radii)
         IF (ASSOCIATED(becke_control%charges_fragment)) &
            DEALLOCATE (becke_control%charges_fragment)
         IF (ASSOCIATED(becke_control%fragments)) &
            DEALLOCATE (becke_control%fragments)
         IF (ASSOCIATED(becke_control%coeff)) &
            DEALLOCATE (becke_control%coeff)
         IF (ASSOCIATED(becke_control%charge)) &
            DEALLOCATE (becke_control%charge)
         IF (ASSOCIATED(becke_control%aij)) &
            DEALLOCATE (becke_control%aij)
         IF (ASSOCIATED(becke_control%gradients)) &
            DEALLOCATE (becke_control%gradients)
         IF (ASSOCIATED(becke_control%cavity_mat)) &
            DEALLOCATE (becke_control%cavity_mat)
         IF (ASSOCIATED(becke_control%skip_list)) &
            DEALLOCATE (becke_control%skip_list)
         IF (ASSOCIATED(becke_control%target)) &
            DEALLOCATE (becke_control%target)
         IF (ASSOCIATED(becke_control%strength)) &
            DEALLOCATE (becke_control%strength)
         IF (ASSOCIATED(becke_control%becke_order_p)) &
            DEALLOCATE (becke_control%becke_order_p)
         IF (ASSOCIATED(becke_control%combined_mat)) &
            DEALLOCATE (becke_control%combined_mat)
         IF (becke_control%cavity_confine) &
            CALL release_hirshfeld_type(becke_control%cavity_env)
         becke_control%ref_count = 0
         becke_control%natoms = 0
         DEALLOCATE (becke_control)
      ENDIF
   END SUBROUTINE becke_control_release

! **************************************************************************************************
!> \brief retain the becke_restraint_type
!> \param becke_control the structure to retain
!> \par History
!>      02.2007 created [Florian Schiffmann]
! **************************************************************************************************
   SUBROUTINE becke_control_retain(becke_control)
      TYPE(becke_restraint_type), POINTER                :: becke_control

      CHARACTER(len=*), PARAMETER :: routineN = 'becke_control_retain', &
         routineP = moduleN//':'//routineN

      CPASSERT(ASSOCIATED(becke_control))
      becke_control%ref_count = becke_control%ref_count+1
   END SUBROUTINE becke_control_retain
! *****************************************************************************
!> \brief create the cdft_control_type
!> \param cdft_control the structure to create
!> \par History
!>      12.2015 created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE cdft_control_create(cdft_control)
      TYPE(cdft_control_type), POINTER                   :: cdft_control

      CHARACTER(len=*), PARAMETER :: routineN = 'cdft_control_create', &
         routineP = moduleN//':'//routineN

      CPASSERT(.NOT. ASSOCIATED(cdft_control))
      ALLOCATE (cdft_control)
      cdft_control%ref_count = 1
      cdft_control%total_steps = 0
      NULLIFY (cdft_control%strength)
      NULLIFY (cdft_control%target)
      NULLIFY (cdft_control%value)
      cdft_control%type = outer_scf_none
      cdft_control%need_pot = .TRUE.
      cdft_control%save_pot = .FALSE.
      cdft_control%do_et = .FALSE.
      cdft_control%reuse_precond = .FALSE.
      cdft_control%nreused = 0
      cdft_control%precond_freq = 0
      cdft_control%max_reuse = 0
      cdft_control%should_purge = .FALSE.
      cdft_control%purge_history = .FALSE.
      cdft_control%calculate_metric = .FALSE.
      cdft_control%purge_freq = 0
      cdft_control%nbad_conv = 0
      cdft_control%purge_offset = 0
      cdft_control%istep = 0
      cdft_control%constraint_type = cdft_density_constraint
      cdft_control%combined_type = cdft_combined_all
      cdft_control%ienergy = 0
      NULLIFY (cdft_control%hirshfeld_control)
      NULLIFY (cdft_control%weight%pw)
      NULLIFY (cdft_control%wmat%matrix)
      NULLIFY (cdft_control%matrix_s%matrix)
      NULLIFY (cdft_control%mo_coeff)
      NULLIFY (cdft_control%matrix_p)
      ! Outer SCF default settings
      cdft_control%ot_control%have_scf = .FALSE.
      cdft_control%ot_control%max_scf = 0
      cdft_control%ot_control%eps_scf = 0.0_dp
      cdft_control%ot_control%step_size = 0.0_dp
      cdft_control%ot_control%type = -1
      cdft_control%ot_control%optimizer = -1
      cdft_control%ot_control%diis_buffer_length = -1
      cdft_control%ot_control%jacobian_type = -1
      cdft_control%ot_control%jacobian_step = 0.0_dp
      cdft_control%ot_control%newton_step = 1.0_dp
      cdft_control%constraint_control%have_scf = .FALSE.
      cdft_control%constraint_control%max_scf = 0
      cdft_control%constraint_control%eps_scf = 0.0_dp
      cdft_control%constraint_control%step_size = 0.0_dp
      cdft_control%constraint_control%type = -1
      cdft_control%constraint_control%optimizer = -1
      cdft_control%constraint_control%diis_buffer_length = -1
      cdft_control%constraint_control%jacobian_type = -1
      cdft_control%constraint_control%broyden_type = -1
      cdft_control%constraint_control%jacobian_freq = 1
      cdft_control%constraint_control%jacobian_step = 0.0_dp
      cdft_control%constraint_control%newton_step = 1.0_dp
      cdft_control%constraint_control%ijacobian = 0
      cdft_control%constraint_control%build_jacobian = .FALSE.
      cdft_control%constraint_control%broyden_update = .FALSE.
      cdft_control%constraint%iter_count = 0
      NULLIFY (cdft_control%constraint%variables)
      NULLIFY (cdft_control%constraint%gradient)
      NULLIFY (cdft_control%constraint%energy)
      NULLIFY (cdft_control%constraint%count)
      NULLIFY (cdft_control%constraint%inv_jacobian)
      cdft_control%constraint%deallocate_jacobian = .TRUE.
   END SUBROUTINE cdft_control_create

! *****************************************************************************
!> \brief release the cdft_control_type
!> \param cdft_control the structure to release
!> \par History
!>      12.2015 created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE cdft_control_release(cdft_control)
      TYPE(cdft_control_type), POINTER                   :: cdft_control

      CHARACTER(len=*), PARAMETER :: routineN = 'cdft_control_release', &
         routineP = moduleN//':'//routineN

      CPASSERT(ASSOCIATED(cdft_control))
      CPASSERT(cdft_control%ref_count > 0)
      cdft_control%ref_count = cdft_control%ref_count-1
      IF (cdft_control%ref_count == 0) THEN
         IF (ASSOCIATED(cdft_control%strength)) &
            DEALLOCATE (cdft_control%strength)
         IF (ASSOCIATED(cdft_control%target)) &
            DEALLOCATE (cdft_control%target)
         IF (ASSOCIATED(cdft_control%value)) &
            DEALLOCATE (cdft_control%value)
         cdft_control%type = outer_scf_none
         NULLIFY (cdft_control%weight%pw)
         IF (ASSOCIATED(cdft_control%hirshfeld_control)) THEN
            CALL release_hirshfeld_type(cdft_control%hirshfeld_control%hirshfeld_env)
            cdft_control%hirshfeld_control%natoms = 0
            IF (ASSOCIATED(cdft_control%hirshfeld_control%atoms)) &
               DEALLOCATE (cdft_control%hirshfeld_control%atoms)
            IF (ASSOCIATED(cdft_control%hirshfeld_control%coeff)) &
               DEALLOCATE (cdft_control%hirshfeld_control%coeff)
            DEALLOCATE (cdft_control%hirshfeld_control)
         END IF
         IF (ASSOCIATED(cdft_control%constraint%variables)) &
            DEALLOCATE (cdft_control%constraint%variables)
         IF (ASSOCIATED(cdft_control%constraint%count)) &
            DEALLOCATE (cdft_control%constraint%count)
         IF (ASSOCIATED(cdft_control%constraint%gradient)) &
            DEALLOCATE (cdft_control%constraint%gradient)
         IF (ASSOCIATED(cdft_control%constraint%energy)) &
            DEALLOCATE (cdft_control%constraint%energy)
         IF (ASSOCIATED(cdft_control%constraint%inv_jacobian)) &
            DEALLOCATE (cdft_control%constraint%inv_jacobian)
         DEALLOCATE (cdft_control)
      END IF
   END SUBROUTINE cdft_control_release

! *****************************************************************************
!> \brief retain the cdft_control_type
!> \param cdft_control the structure to retain
!> \par History
!>      created 12.2015 [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE cdft_control_retain(cdft_control)
      TYPE(cdft_control_type), POINTER                   :: cdft_control

      CHARACTER(len=*), PARAMETER :: routineN = 'cdft_control_retain', &
         routineP = moduleN//':'//routineN

      CPASSERT(ASSOCIATED(cdft_control))
      cdft_control%ref_count = cdft_control%ref_count+1
   END SUBROUTINE cdft_control_retain

! **************************************************************************************************
!> \brief create the s2_restraint_type
!> \param s2_restraint_control ...
!> \par History
!>      03.2006 created [Joost VandeVondele]
! **************************************************************************************************
   SUBROUTINE s2_control_create(s2_restraint_control)
      TYPE(s2_restraint_type), POINTER                   :: s2_restraint_control

      CHARACTER(len=*), PARAMETER :: routineN = 's2_control_create', &
         routineP = moduleN//':'//routineN

      CPASSERT(.NOT. ASSOCIATED(s2_restraint_control))
      ALLOCATE (s2_restraint_control)

      s2_restraint_control%ref_count = 1
      s2_restraint_control%strength = 0.1_dp
      s2_restraint_control%s2_order_p = 0.0_dp
      s2_restraint_control%functional_form = -1
      s2_restraint_control%target = 1.0_dp
   END SUBROUTINE s2_control_create

! **************************************************************************************************
!> \brief release the s2_restraint_type
!> \param s2_restraint_control ...
!> \par History
!>      03.2006 created [Joost VandeVondele]
! **************************************************************************************************
   SUBROUTINE s2_control_release(s2_restraint_control)
      TYPE(s2_restraint_type), POINTER                   :: s2_restraint_control

      CHARACTER(len=*), PARAMETER :: routineN = 's2_control_release', &
         routineP = moduleN//':'//routineN

      CPASSERT(ASSOCIATED(s2_restraint_control))
      CPASSERT(s2_restraint_control%ref_count > 0)
      s2_restraint_control%ref_count = s2_restraint_control%ref_count-1
      IF (s2_restraint_control%ref_count == 0) THEN
         s2_restraint_control%ref_count = 0
         s2_restraint_control%strength = 0.0_dp
         s2_restraint_control%target = 0.0_dp
         DEALLOCATE (s2_restraint_control)
      ENDIF
   END SUBROUTINE s2_control_release

! **************************************************************************************************
!> \brief retain the s2_restraint_type
!> \param s2_restraint_control ...
!> \par History
!>      03.2006 created [Joost VandeVondele]
! **************************************************************************************************
   SUBROUTINE s2_control_retain(s2_restraint_control)
      TYPE(s2_restraint_type), POINTER                   :: s2_restraint_control

      CHARACTER(len=*), PARAMETER :: routineN = 's2_control_retain', &
         routineP = moduleN//':'//routineN

      CPASSERT(ASSOCIATED(s2_restraint_control))
      s2_restraint_control%ref_count = s2_restraint_control%ref_count+1
   END SUBROUTINE s2_control_retain

! **************************************************************************************************
!> \brief allocates and perform a very basic initialization
!> \param dft_control the object to create
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE dft_control_create(dft_control)
      TYPE(dft_control_type), POINTER                    :: dft_control

      CHARACTER(len=*), PARAMETER :: routineN = 'dft_control_create', &
         routineP = moduleN//':'//routineN

      CPASSERT(.NOT. ASSOCIATED(dft_control))
      ALLOCATE (dft_control)
      dft_control%ref_count = 1
      last_dft_control_id = last_dft_control_id+1
      dft_control%id_nr = last_dft_control_id
      NULLIFY (dft_control%xas_control)
      NULLIFY (dft_control%qs_control)
      NULLIFY (dft_control%tddfpt_control)
      NULLIFY (dft_control%tddfpt2_control)
      NULLIFY (dft_control%efield_fields)
      NULLIFY (dft_control%period_efield)
      NULLIFY (dft_control%admm_control)
      NULLIFY (dft_control%rtp_control)
      NULLIFY (dft_control%sccs_control)
      dft_control%do_sccs = .FALSE.
      CALL qs_control_create(dft_control%qs_control)
      CALL tddfpt2_control_create(dft_control%tddfpt2_control)
   END SUBROUTINE dft_control_create

! **************************************************************************************************
!> \brief ...
!> \param dft_control ...
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE dft_control_retain(dft_control)
      TYPE(dft_control_type), POINTER                    :: dft_control

      CHARACTER(len=*), PARAMETER :: routineN = 'dft_control_retain', &
         routineP = moduleN//':'//routineN

      CPASSERT(ASSOCIATED(dft_control))
      CPASSERT(dft_control%ref_count > 0)
      dft_control%ref_count = dft_control%ref_count+1
   END SUBROUTINE dft_control_retain

! **************************************************************************************************
!> \brief ...
!> \param dft_control ...
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE dft_control_release(dft_control)
      TYPE(dft_control_type), POINTER                    :: dft_control

      CHARACTER(len=*), PARAMETER :: routineN = 'dft_control_release', &
         routineP = moduleN//':'//routineN

      IF (ASSOCIATED(dft_control)) THEN
         CPASSERT(dft_control%ref_count > 0)
         dft_control%ref_count = dft_control%ref_count-1
         IF (dft_control%ref_count == 0) THEN
            CALL qs_control_release(dft_control%qs_control)
            CALL tddfpt_control_release(dft_control%tddfpt_control)
            CALL tddfpt2_control_release(dft_control%tddfpt2_control)
            CALL xas_control_release(dft_control%xas_control)
            CALL admm_control_release(dft_control%admm_control)
            CALL efield_fields_release(dft_control%efield_fields)
            CALL sccs_control_release(dft_control%sccs_control)
            IF (ASSOCIATED(dft_control%period_efield)) THEN
               DEALLOCATE (dft_control%period_efield)
            END IF
            IF (ASSOCIATED(dft_control%rtp_control)) THEN
               DEALLOCATE (dft_control%rtp_control)
            END IF
            DEALLOCATE (dft_control)
         END IF
      END IF

   END SUBROUTINE dft_control_release

! **************************************************************************************************
!> \brief ...
!> \param gapw_control ...
! **************************************************************************************************
   SUBROUTINE gapw_control_create(gapw_control)
      TYPE(gapw_control_type), POINTER                   :: gapw_control

      CHARACTER(len=*), PARAMETER :: routineN = 'gapw_control_create', &
         routineP = moduleN//':'//routineN

      CPASSERT(.NOT. ASSOCIATED(gapw_control))
      ALLOCATE (gapw_control)
   END SUBROUTINE gapw_control_create

! **************************************************************************************************
!> \brief ...
!> \param qs_control ...
! **************************************************************************************************
   SUBROUTINE qs_control_create(qs_control)
      TYPE(qs_control_type), POINTER                     :: qs_control

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_control_create', &
         routineP = moduleN//':'//routineN

      CPASSERT(.NOT. ASSOCIATED(qs_control))
      ALLOCATE (qs_control)

      NULLIFY (qs_control%e_cutoff)
      NULLIFY (qs_control%gapw_control)
      NULLIFY (qs_control%mulliken_restraint_control)
      NULLIFY (qs_control%ddapc_restraint_control)
      NULLIFY (qs_control%s2_restraint_control)
      NULLIFY (qs_control%se_control)
      NULLIFY (qs_control%dftb_control)
      NULLIFY (qs_control%scptb_control)
      NULLIFY (qs_control%becke_control)
      NULLIFY (qs_control%cdft_control)
      NULLIFY (qs_control%ddapc_restraint_control)

      CALL mulliken_control_create(qs_control%mulliken_restraint_control)
      CALL becke_control_create(qs_control%becke_control)
      CALL s2_control_create(qs_control%s2_restraint_control)
      CALL gapw_control_create(qs_control%gapw_control)
      CALL se_control_create(qs_control%se_control)
      CALL dftb_control_create(qs_control%dftb_control)
      CALL scptb_control_create(qs_control%scptb_control)
      CALL cdft_control_create(qs_control%cdft_control)
   END SUBROUTINE qs_control_create

! **************************************************************************************************
!> \brief ...
!> \param qs_control ...
! **************************************************************************************************
   SUBROUTINE qs_control_release(qs_control)
      TYPE(qs_control_type), POINTER                     :: qs_control

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_control_release', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i

      IF (ASSOCIATED(qs_control)) THEN
         CALL mulliken_control_release(qs_control%mulliken_restraint_control)
         CALL s2_control_release(qs_control%s2_restraint_control)
         CALL se_control_release(qs_control%se_control)
         CALL dftb_control_release(qs_control%dftb_control)
         CALL scptb_control_release(qs_control%scptb_control)
         CALL becke_control_release(qs_control%becke_control)
         CALL cdft_control_release(qs_control%cdft_control)

         IF (ASSOCIATED(qs_control%e_cutoff)) THEN
            DEALLOCATE (qs_control%e_cutoff)
         END IF
         IF (ASSOCIATED(qs_control%gapw_control)) THEN
            DEALLOCATE (qs_control%gapw_control)
         END IF
         IF (ASSOCIATED(qs_control%ddapc_restraint_control)) THEN
            DO i = 1, SIZE(qs_control%ddapc_restraint_control)
               CALL ddapc_control_release(qs_control%ddapc_restraint_control(i)%ddapc_restraint_control)
            END DO
            DEALLOCATE (qs_control%ddapc_restraint_control)
         END IF
         DEALLOCATE (qs_control)
      END IF
   END SUBROUTINE qs_control_release

! **************************************************************************************************
!> \brief ...
!> \param tddfpt_control ...
! **************************************************************************************************
   SUBROUTINE tddfpt_control_create(tddfpt_control)
      TYPE(tddfpt_control_type), POINTER                 :: tddfpt_control

      CHARACTER(len=*), PARAMETER :: routineN = 'tddfpt_control_create', &
         routineP = moduleN//':'//routineN

      CPASSERT(.NOT. ASSOCIATED(tddfpt_control))
      ALLOCATE (tddfpt_control)
      NULLIFY (tddfpt_control%lumos)
      NULLIFY (tddfpt_control%lumos_eigenvalues)

   END SUBROUTINE tddfpt_control_create

! **************************************************************************************************
!> \brief ...
!> \param tddfpt_control ...
! **************************************************************************************************
   SUBROUTINE tddfpt_control_release(tddfpt_control)
      TYPE(tddfpt_control_type), POINTER                 :: tddfpt_control

      CHARACTER(len=*), PARAMETER :: routineN = 'tddfpt_control_release', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ispin
      LOGICAL                                            :: dummy

      IF (ASSOCIATED(tddfpt_control)) THEN
         IF (ASSOCIATED(tddfpt_control%lumos)) THEN
            DO ispin = 1, SIZE(tddfpt_control%lumos)
               CALL cp_fm_release(tddfpt_control%lumos(ispin)%matrix)
               !MK the following line just avoids a crash of TDDFT runs using
               !MK the sdbg version compiled with the NAG compiler when
               !MK tddfpt_control%lumos is deallocated. This is most likely a
               !MK compiler bug and thus the line might become obsolete
               dummy = ASSOCIATED(tddfpt_control%lumos(ispin)%matrix)
            END DO
            DEALLOCATE (tddfpt_control%lumos)
         END IF
         IF (ASSOCIATED(tddfpt_control%lumos_eigenvalues)) THEN
            DEALLOCATE (tddfpt_control%lumos_eigenvalues)
         END IF
         DEALLOCATE (tddfpt_control)
      END IF
   END SUBROUTINE tddfpt_control_release

! **************************************************************************************************
!> \brief allocate control options for Time-Dependent Density Functional Theory calculation
!> \param tddfpt_control an object to create
!> \par History
!>    * 05.2016 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt2_control_create(tddfpt_control)
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control

      CHARACTER(len=*), PARAMETER :: routineN = 'tddfpt2_control_create', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CPASSERT(.NOT. ASSOCIATED(tddfpt_control))
      CALL timeset(routineN, handle)

      ALLOCATE (tddfpt_control)

      CALL timestop(handle)
   END SUBROUTINE tddfpt2_control_create

! **************************************************************************************************
!> \brief release memory allocated for TDDFT control options
!> \param tddfpt_control an object to release
!> \par History
!>    * 05.2016 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt2_control_release(tddfpt_control)
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control

      CHARACTER(len=*), PARAMETER :: routineN = 'tddfpt2_control_release', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (ASSOCIATED(tddfpt_control)) THEN
         DEALLOCATE (tddfpt_control)
      END IF

      CALL timestop(handle)
   END SUBROUTINE tddfpt2_control_release

! **************************************************************************************************
!> \brief ...
!> \param efield_fields ...
! **************************************************************************************************
   SUBROUTINE efield_fields_release(efield_fields)
      TYPE(efield_p_type), DIMENSION(:), POINTER         :: efield_fields

      CHARACTER(len=*), PARAMETER :: routineN = 'efield_fields_release', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i

      IF (ASSOCIATED(efield_fields)) THEN
         DO i = 1, SIZE(efield_fields)
            IF (ASSOCIATED(efield_fields(i)%efield)) THEN
               IF (ASSOCIATED(efield_fields(i)%efield%envelop_r_vars)) THEN
                  DEALLOCATE (efield_fields(i)%efield%envelop_r_vars)
               END IF
               IF (ASSOCIATED(efield_fields(i)%efield%envelop_i_vars)) THEN
                  DEALLOCATE (efield_fields(i)%efield%envelop_i_vars)
               END IF
               IF (ASSOCIATED(efield_fields(i)%efield%polarisation)) THEN
                  DEALLOCATE (efield_fields(i)%efield%polarisation)
               END IF
               DEALLOCATE (efield_fields(i)%efield)
            END IF
         END DO
         DEALLOCATE (efield_fields)
      END IF
   END SUBROUTINE efield_fields_release

! **************************************************************************************************
!> \brief ...
!> \param dftb_control ...
! **************************************************************************************************
   SUBROUTINE dftb_control_create(dftb_control)
      TYPE(dftb_control_type), POINTER                   :: dftb_control

      CHARACTER(len=*), PARAMETER :: routineN = 'dftb_control_create', &
         routineP = moduleN//':'//routineN

      CPASSERT(.NOT. ASSOCIATED(dftb_control))
      ALLOCATE (dftb_control)

      NULLIFY (dftb_control%sk_pair_list)
   END SUBROUTINE dftb_control_create

! **************************************************************************************************
!> \brief ...
!> \param dftb_control ...
! **************************************************************************************************
   SUBROUTINE dftb_control_release(dftb_control)
      TYPE(dftb_control_type), POINTER                   :: dftb_control

      CHARACTER(len=*), PARAMETER :: routineN = 'dftb_control_release', &
         routineP = moduleN//':'//routineN

      IF (ASSOCIATED(dftb_control)) THEN
         IF (ASSOCIATED(dftb_control%sk_pair_list)) THEN
            DEALLOCATE (dftb_control%sk_pair_list)
         END IF
         DEALLOCATE (dftb_control)
      END IF
   END SUBROUTINE dftb_control_release

! **************************************************************************************************
!> \brief ...
!> \param scptb_control ...
! **************************************************************************************************
   SUBROUTINE scptb_control_create(scptb_control)
      TYPE(scptb_control_type), POINTER                  :: scptb_control

      CHARACTER(len=*), PARAMETER :: routineN = 'scptb_control_create', &
         routineP = moduleN//':'//routineN

      CPASSERT(.NOT. ASSOCIATED(scptb_control))
      ALLOCATE (scptb_control)

   END SUBROUTINE scptb_control_create

! **************************************************************************************************
!> \brief ...
!> \param scptb_control ...
! **************************************************************************************************
   SUBROUTINE scptb_control_release(scptb_control)
      TYPE(scptb_control_type), POINTER                  :: scptb_control

      CHARACTER(len=*), PARAMETER :: routineN = 'scptb_control_release', &
         routineP = moduleN//':'//routineN

      IF (ASSOCIATED(scptb_control)) THEN
         DEALLOCATE (scptb_control)
      END IF
   END SUBROUTINE scptb_control_release

! **************************************************************************************************
!> \brief ...
!> \param se_control ...
! **************************************************************************************************
   SUBROUTINE se_control_create(se_control)
      TYPE(semi_empirical_control_type), POINTER         :: se_control

      CHARACTER(len=*), PARAMETER :: routineN = 'se_control_create', &
         routineP = moduleN//':'//routineN

      CPASSERT(.NOT. ASSOCIATED(se_control))
      ALLOCATE (se_control)
   END SUBROUTINE se_control_create

! **************************************************************************************************
!> \brief ...
!> \param se_control ...
! **************************************************************************************************
   SUBROUTINE se_control_release(se_control)
      TYPE(semi_empirical_control_type), POINTER         :: se_control

      CHARACTER(len=*), PARAMETER :: routineN = 'se_control_release', &
         routineP = moduleN//':'//routineN

      IF (ASSOCIATED(se_control)) THEN
         DEALLOCATE (se_control)
      END IF
   END SUBROUTINE se_control_release

! **************************************************************************************************
!> \brief ...
!> \param admm_control ...
! **************************************************************************************************
   SUBROUTINE admm_control_create(admm_control)
      TYPE(admm_control_type), POINTER                   :: admm_control

      CHARACTER(len=*), PARAMETER :: routineN = 'admm_control_create', &
         routineP = moduleN//':'//routineN

      CPASSERT(.NOT. ASSOCIATED(admm_control))
      ALLOCATE (admm_control)

   END SUBROUTINE admm_control_create

! **************************************************************************************************
!> \brief ...
!> \param admm_control ...
! **************************************************************************************************
   SUBROUTINE admm_control_release(admm_control)
      TYPE(admm_control_type), POINTER                   :: admm_control

      CHARACTER(len=*), PARAMETER :: routineN = 'admm_control_release', &
         routineP = moduleN//':'//routineN

      IF (ASSOCIATED(admm_control)) THEN
         DEALLOCATE (admm_control)
      END IF
   END SUBROUTINE admm_control_release

! **************************************************************************************************
!> \brief   Create sccs_control_type
!> \param sccs_control ...
!> \par     History
!>          - Creation (11.10.2013,MK)
!> \author  Matthias Krack (MK)
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE sccs_control_create(sccs_control)
      TYPE(sccs_control_type), POINTER                   :: sccs_control

      CHARACTER(LEN=*), PARAMETER :: routineN = 'sccs_control_create', &
         routineP = moduleN//':'//routineN

      CPASSERT(.NOT. ASSOCIATED(sccs_control))
      ALLOCATE (sccs_control)

      sccs_control%ref_count = 1
      sccs_control%derivative_method = 0
      sccs_control%max_iter = 0
      sccs_control%method_id = 0
      sccs_control%alpha_solvent = 0.0_dp
      sccs_control%beta = 0.0_dp
      sccs_control%beta_solvent = 0.0_dp
      sccs_control%delta_rho = 0.0_dp
      sccs_control%eps_sccs = 0.0_dp
      sccs_control%eps_scf = 0.0_dp
      sccs_control%epsilon_solvent = 0.0_dp
      sccs_control%gamma_solvent = 0.0_dp
      sccs_control%mixing = 0.0_dp
      sccs_control%rho_max = 0.0_dp
      sccs_control%rho_min = 0.0_dp
      sccs_control%rho_zero = 0.0_dp
      sccs_control%sccs_activated = .FALSE.

   END SUBROUTINE sccs_control_create

! **************************************************************************************************
!> \brief   Release sccs_control_type
!> \param sccs_control ...
!> \par     History
!>          - Creation (11.10.2013,MK)
!> \author  Matthias Krack (MK)
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE sccs_control_release(sccs_control)
      TYPE(sccs_control_type), POINTER                   :: sccs_control

      CHARACTER(LEN=*), PARAMETER :: routineN = 'sccs_control_release', &
         routineP = moduleN//':'//routineN

      IF (ASSOCIATED(sccs_control)) THEN
         CPASSERT(sccs_control%ref_count > 0)
         sccs_control%ref_count = sccs_control%ref_count-1
         IF (sccs_control%ref_count == 0) THEN
            DEALLOCATE (sccs_control)
         END IF
      END IF

   END SUBROUTINE sccs_control_release

END MODULE cp_control_types
