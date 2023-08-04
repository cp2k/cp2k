!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Defines control structures, which contain the parameters and the
!>      settings for the DFT-based calculations.
! **************************************************************************************************
MODULE cp_control_types
   USE cp_fm_types,                     ONLY: cp_fm_release,&
                                              cp_fm_type
   USE input_constants,                 ONLY: do_full_density
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE pair_potential_types,            ONLY: pair_potential_p_release,&
                                              pair_potential_p_type
   USE qs_cdft_types,                   ONLY: cdft_control_create,&
                                              cdft_control_release,&
                                              cdft_control_type
   USE xas_control,                     ONLY: xas_control_release,&
                                              xas_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! **************************************************************************************************
! \brief Control parameters for pw grids
! **************************************************************************************************
   TYPE pw_grid_option
      LOGICAL                              :: spherical = .FALSE.
      LOGICAL                              :: fullspace = .FALSE.
      INTEGER, DIMENSION(2)                :: distribution_layout = 0
      INTEGER                              :: blocked = 0
   END TYPE pw_grid_option
! **************************************************************************************************
! \brief parameters for EMD/RTP calculations involving MO projections
! **************************************************************************************************
   TYPE proj_mo_type
      INTEGER, DIMENSION(:), ALLOCATABLE         :: ref_mo_index
      INTEGER                                    :: ref_mo_spin = 1
      LOGICAL                                    :: sum_on_all_ref = .FALSE.
      INTEGER, DIMENSION(:), ALLOCATABLE         :: td_mo_index
      INTEGER                                    :: td_mo_spin = 1
      LOGICAL                                    :: sum_on_all_td = .FALSE.
      CHARACTER(LEN=default_path_length)         :: ref_mo_file_name = ""
      TYPE(cp_fm_type), DIMENSION(:), &
         ALLOCATABLE                       :: mo_ref
   END TYPE proj_mo_type

   TYPE proj_mo_p_type
      TYPE(proj_mo_type), POINTER                :: proj_mo => NULL()
   END TYPE proj_mo_p_type
! **************************************************************************************************
! \brief Control parameters for REAL_TIME_PROPAGATION calculations
! **************************************************************************************************
   TYPE rtp_control_type
      LOGICAL                              :: converged = .FALSE.
      REAL(KIND=dp)                        :: eps_ener = 0.0_dp
      INTEGER                              :: max_iter = 0
      INTEGER                              :: mat_exp = 0
      INTEGER                              :: propagator = 0
      LOGICAL                              :: fixed_ions = .FALSE.
      INTEGER                              :: initial_wfn = 0
      REAL(dp)                             :: eps_exp = 0.0_dp
      LOGICAL                              :: initial_step = .FALSE.
      LOGICAL                              :: hfx_redistribute = .FALSE.
      INTEGER                              :: aspc_order = 0
      INTEGER                              :: sc_check_start = 0
      LOGICAL                              :: apply_delta_pulse = .FALSE.
      LOGICAL                              :: apply_delta_pulse_mag = .FALSE.
      LOGICAL                              :: periodic = .FALSE.
      LOGICAL                              :: linear_scaling = .FALSE.
      LOGICAL                              :: write_restart = .FALSE.
      INTEGER                              :: mcweeny_max_iter = 0
      INTEGER                              :: acc_ref = 0
      REAL(dp)                             :: mcweeny_eps = 0.0_dp
      INTEGER, DIMENSION(3)                :: delta_pulse_direction = 0
      REAL(KIND=dp)                        :: delta_pulse_scale = 0.0_dp
      LOGICAL                              :: velocity_gauge = .FALSE.
      REAL(KIND=dp), DIMENSION(3)          :: vec_pot = 0.0_dp
      LOGICAL                              :: nl_gauge_transform = .FALSE.
      LOGICAL                              :: is_proj_mo = .FALSE.
      TYPE(proj_mo_p_type), DIMENSION(:), &
         POINTER                :: proj_mo_list => NULL()
   END TYPE rtp_control_type
! **************************************************************************************************
! \brief Control parameters for DFTB calculations
! **************************************************************************************************
   TYPE dftb_control_type
      LOGICAL                              :: self_consistent = .FALSE.
      LOGICAL                              :: orthogonal_basis = .FALSE.
      LOGICAL                              :: dispersion = .FALSE.
      INTEGER                              :: dispersion_type = 0
      LOGICAL                              :: dftb3_diagonal = .FALSE.
      LOGICAL                              :: hb_sr_damp = .FALSE.
      REAL(KIND=dp)                        :: hb_sr_para = 0.0_dp
      REAL(KIND=dp)                        :: eps_disp = 0.0_dp
      REAL(KIND=dp)                        :: epscn = 0.0_dp
      REAL(KIND=dp)                        :: exp_pre = 0.0_dp
      REAL(KIND=dp)                        :: scaling = 0.0_dp
      REAL(KIND=dp)                        :: rcdisp = 0.0_dp
      REAL(KIND=dp), DIMENSION(3)          :: sd3 = 0.0_dp
      REAL(KIND=dp), DIMENSION(4)          :: sd3bj = 0.0_dp
      LOGICAL                              :: do_ewald = .FALSE.
      CHARACTER(LEN=default_path_length)   :: sk_file_path = ""
      CHARACTER(LEN=default_path_length)   :: sk_file_list = ""
      CHARACTER(LEN=default_string_length), &
         DIMENSION(:, :), POINTER          :: sk_pair_list => NULL()
      CHARACTER(LEN=default_path_length)   :: uff_force_field = ""
      CHARACTER(LEN=default_path_length)   :: dispersion_parameter_file = ""
   END TYPE dftb_control_type

! **************************************************************************************************
! \brief Control parameters for xTB calculations
! **************************************************************************************************
   TYPE xtb_control_type
      !
      LOGICAL                              :: do_ewald = .FALSE.
      !
      INTEGER                              :: sto_ng = 0
      INTEGER                              :: h_sto_ng = 0
      !
      CHARACTER(LEN=default_path_length)   :: parameter_file_path = ""
      CHARACTER(LEN=default_path_length)   :: parameter_file_name = ""
      !
      CHARACTER(LEN=default_path_length)   :: dispersion_parameter_file = ""
      REAL(KIND=dp)                        :: epscn = 0.0_dp
      REAL(KIND=dp)                        :: rcdisp = 0.0_dp
      REAL(KIND=dp)                        :: s6 = 0.0_dp, s8 = 0.0_dp
      REAL(KIND=dp)                        :: a1 = 0.0_dp, a2 = 0.0_dp
      !
      REAL(KIND=dp)                        :: ks = 0.0_dp, kp = 0.0_dp, kd = 0.0_dp, ksp = 0.0_dp, k2sh = 0.0_dp
      REAL(KIND=dp)                        :: kg = 0.0_dp, kf = 0.0_dp
      REAL(KIND=dp)                        :: kcns = 0.0_dp, kcnp = 0.0_dp, kcnd = 0.0_dp
      REAL(KIND=dp)                        :: ken = 0.0_dp
      REAL(KIND=dp)                        :: kxr = 0.0_dp, kx2 = 0.0_dp
      !
      LOGICAL                              :: xb_interaction = .FALSE.
      LOGICAL                              :: do_nonbonded = .FALSE.
      LOGICAL                              :: coulomb_interaction = .FALSE.
      LOGICAL                              :: coulomb_lr = .FALSE.
      LOGICAL                              :: tb3_interaction = .FALSE.
      LOGICAL                              :: check_atomic_charges = .FALSE.
      LOGICAL                              :: old_coulomb_damping = .FALSE.
      !
      REAL(KIND=dp)                        :: xb_radius = 0.0_dp
      REAL(KIND=dp)                        :: coulomb_sr_cut = 0.0_dp
      REAL(KIND=dp)                        :: coulomb_sr_eps = 0.0_dp
      !
      CHARACTER(LEN=default_string_length), &
         DIMENSION(:, :), POINTER          :: kab_param => NULL()
      INTEGER, DIMENSION(:, :), POINTER    :: kab_types => NULL()
      INTEGER                              :: kab_nval = 0
      REAL, DIMENSION(:), POINTER          :: kab_vals => NULL()
      !
      TYPE(pair_potential_p_type), POINTER :: nonbonded => NULL()
   END TYPE xtb_control_type

! **************************************************************************************************
! \brief Control parameters for semi empirical calculations
! **************************************************************************************************
   TYPE semi_empirical_control_type
      LOGICAL                              :: orthogonal_basis = .FALSE.
      LOGICAL                              :: analytical_gradients = .FALSE.
      LOGICAL                              :: force_kdsod_EX = .FALSE.
      LOGICAL                              :: do_ewald = .FALSE., do_ewald_r3 = .FALSE., do_ewald_gks = .FALSE.
      INTEGER                              :: integral_screening = 0, periodic_type = 0
      INTEGER                              :: max_multipole = 0
      INTEGER                              :: ga_ncells = 0
      REAL(KIND=dp)                        :: delta = 0.0_dp
      ! Dispersion pair potential
      LOGICAL                              :: dispersion = .FALSE.
      REAL(KIND=dp)                        :: rcdisp = 0.0_dp
      REAL(KIND=dp)                        :: epscn = 0.0_dp
      REAL(KIND=dp), DIMENSION(3)          :: sd3 = 0.0_dp
      CHARACTER(LEN=default_path_length)   :: dispersion_parameter_file = ""
      ! Parameters controlling the evaluation of the integrals
      REAL(KIND=dp)                        :: cutoff_lrc = 0.0_dp, taper_lrc = 0.0_dp, range_lrc = 0.0_dp
      REAL(KIND=dp)                        :: cutoff_cou = 0.0_dp, taper_cou = 0.0_dp, range_cou = 0.0_dp
      REAL(KIND=dp)                        :: cutoff_exc = 0.0_dp, taper_exc = 0.0_dp, range_exc = 0.0_dp
      REAL(KIND=dp)                        :: taper_scr = 0.0_dp, range_scr = 0.0_dp
   END TYPE semi_empirical_control_type

! **************************************************************************************************
! \brief Control parameters for GAPW method within QUICKSTEP ***
! **************************************************************************************************
   TYPE gapw_control_type
      INTEGER                              :: basis_1c = 0
      REAL(KIND=dp)                        :: eps_fit = 0.0_dp, &
                                              eps_iso = 0.0_dp, &
                                              eps_Vrho0 = 0.0_dp, &
                                              eps_svd = 0.0_dp, &
                                              eps_cpc = 0.0_dp
      INTEGER                              :: ladd_rho0 = 0, &
                                              lmax_rho0 = 0, &
                                              lmax_sphere = 0, &
                                              quadrature = 0
      LOGICAL                              :: lrho1_eq_lrho0 = .FALSE.
      LOGICAL                              :: alpha0_hard_from_input = .FALSE., &
                                              force_paw = .FALSE., &
                                              non_paw_atoms = .FALSE., &
                                              nopaw_as_gpw = .FALSE.
      REAL(KIND=dp)                        :: alpha0_hard = 0.0_dp
      REAL(KIND=dp)                        :: max_rad_local = 0.0_dp
   END TYPE gapw_control_type
! **************************************************************************************************
! \brief parameters for calculations involving a time dependent electric field
! **************************************************************************************************
   TYPE efield_type
      REAL(KIND=dp)                        :: actual_time = 0.0_dp
      REAL(KIND=dp), DIMENSION(:), POINTER :: polarisation => NULL()
      INTEGER                              :: envelop_id = 0
      REAL(KIND=dp), DIMENSION(:), POINTER :: envelop_r_vars => NULL()
      INTEGER, DIMENSION(:), POINTER       :: envelop_i_vars => NULL()
      REAL(KIND=dp)                        :: strength = 0.0_dp
      REAL(KIND=dp)                        :: phase_offset = 0.0_dp
      REAL(KIND=dp)                        :: wavelength = 0.0_dp
      REAL(KIND=dp), DIMENSION(3)          :: vec_pot_initial = 0.0_dp
   END TYPE efield_type

   TYPE efield_p_type
      TYPE(efield_type), POINTER           :: efield => NULL()
   END TYPE efield_p_type
! **************************************************************************************************
! \brief parameters for calculations involving a time dependent electric field
! **************************************************************************************************
   TYPE period_efield_type
      LOGICAL                              :: displacement_field = .FALSE.
      REAL(KIND=dp), DIMENSION(3)          :: polarisation = 0.0_dp
      REAL(KIND=dp), DIMENSION(3)          :: d_filter = 0.0_dp
      REAL(KIND=dp)                        :: strength = 0.0_dp
   END TYPE period_efield_type

! **************************************************************************************************
! \brief some parameters useful for mulliken_restraints
! **************************************************************************************************
   TYPE mulliken_restraint_type
      REAL(KIND=dp)                        :: strength = 0.0_dp
      REAL(KIND=dp)                        :: TARGET = 0.0_dp
      INTEGER                              :: natoms = 0
      INTEGER, POINTER, DIMENSION(:)       :: atoms => NULL()
   END TYPE mulliken_restraint_type

! **************************************************************************************************
! \brief some parameters useful for ddapc_restraints
! **************************************************************************************************
   TYPE ddapc_restraint_type
      INTEGER                              :: ref_count = 0
      REAL(KIND=dp)                        :: strength = 0.0_dp
      REAL(KIND=dp)                        :: TARGET = 0.0_dp
      REAL(KIND=dp)                        :: ddapc_order_p = 0.0_dp
      INTEGER                              :: functional_form = 0
      INTEGER                              :: natoms = 0
      INTEGER, POINTER, DIMENSION(:)       :: atoms => NULL()
      REAL(KIND=dp), POINTER, DIMENSION(:) :: coeff => NULL()
      INTEGER                              :: density_type = 0
   END TYPE ddapc_restraint_type

! **************************************************************************************************
! \brief some parameters useful for s2_restraints
! **************************************************************************************************
   TYPE s2_restraint_type
      REAL(KIND=dp)                        :: strength = 0.0_dp
      REAL(KIND=dp)                        :: TARGET = 0.0_dp
      REAL(KIND=dp)                        :: s2_order_p = 0.0_dp
      INTEGER                              :: functional_form = 0
   END TYPE s2_restraint_type

! **************************************************************************************************
! \brief some parameters useful for auxiliary density matrix method
! **************************************************************************************************
   TYPE admm_block_type
      INTEGER, DIMENSION(:), ALLOCATABLE   :: list
   END TYPE admm_block_type

   TYPE admm_control_type
      REAL(KIND=dp)                        :: eps_filter = 0.0_dp
      INTEGER                              :: purification_method = 0
      INTEGER                              :: method = 0
      LOGICAL                              :: charge_constrain = .FALSE.
      INTEGER                              :: scaling_model = 0
      INTEGER                              :: aux_exch_func = 0
      LOGICAL                              :: aux_exch_func_param = .FALSE.
      REAL(KIND=dp), DIMENSION(3)          :: aux_x_param = 0.0_dp
      TYPE(admm_block_type), DIMENSION(:), &
         ALLOCATABLE                       :: blocks
   END TYPE admm_control_type

! **************************************************************************************************
! \brief Parameters for external potential
! **************************************************************************************************
   TYPE expot_control_type
      LOGICAL                              :: read_from_cube = .FALSE.
      LOGICAL                              :: maxwell_solver = .FALSE.
      LOGICAL                              :: static = .FALSE.
      REAL(KIND=dp)                        :: scaling_factor = 0.0_dp
   END TYPE expot_control_type

! **************************************************************************************************
! \brief Parameters useful for Maxwell equation evaluation of external potential
! **************************************************************************************************
   TYPE maxwell_control_type
      LOGICAL                              :: log_test = .FALSE.
      INTEGER                              :: int_test = 0
      REAL(KIND=dp)                        :: real_test = 0.0_dp
   END TYPE maxwell_control_type

! **************************************************************************************************
! \brief Control parameters for a QUICKSTEP and KIM-GORDON calculation ***
!        eps_pgf_orb: Cutoff value for the interaction of the primitive
!                     Gaussian-type functions (primitive basis functions).
! **************************************************************************************************
   TYPE qs_control_type
      INTEGER                              :: method_id = 0
      REAL(KIND=dp)                        :: eps_core_charge = 0.0_dp, &
                                              eps_kg_orb = 0.0_dp, &
                                              eps_pgf_orb = 0.0_dp, &
                                              eps_ppl = 0.0_dp, &
                                              eps_ppnl = 0.0_dp, &
                                              eps_rho_gspace = 0.0_dp, &
                                              eps_rho_rspace = 0.0_dp, &
                                              eps_filter_matrix = 0.0_dp, &
                                              eps_gvg_rspace = 0.0_dp, &
                                              progression_factor = 0.0_dp, &
                                              relative_cutoff = 0.0_dp
      LOGICAL                              :: do_almo_scf = .FALSE.
      LOGICAL                              :: do_ls_scf = .FALSE.
      LOGICAL                              :: do_kg = .FALSE.
      LOGICAL                              :: commensurate_mgrids = .FALSE.
      LOGICAL                              :: realspace_mgrids = .FALSE.
      LOGICAL                              :: gapw = .FALSE., gapw_xc = .FALSE., gpw = .FALSE., pao = .FALSE.
      LOGICAL                              :: lrigpw = .FALSE., rigpw = .FALSE.
      LOGICAL                              :: lri_optbas = .FALSE.
      LOGICAL                              :: ofgpw = .FALSE.
      LOGICAL                              :: dftb = .FALSE.
      LOGICAL                              :: xtb = .FALSE.
      LOGICAL                              :: semi_empirical = .FALSE.
      LOGICAL                              :: mulliken_restraint = .FALSE.
      LOGICAL                              :: ddapc_restraint = .FALSE.
      LOGICAL                              :: ddapc_restraint_is_spin = .FALSE.
      LOGICAL                              :: ddapc_explicit_potential = .FALSE.
      LOGICAL                              :: cdft = .FALSE.
      LOGICAL                              :: et_coupling_calc = .FALSE.
      LOGICAL                              :: s2_restraint = .FALSE.
      INTEGER                              :: do_ppl_method = 0
      INTEGER                              :: wf_interpolation_method_nr = 0
      INTEGER                              :: wf_extrapolation_order = 0
      INTEGER                              :: periodicity = 0
      REAL(KIND=dp)                        :: pairlist_radius = 0.0_dp
      REAL(KIND=dp)                        :: cutoff = 0.0_dp
      REAL(KIND=dp), DIMENSION(:), POINTER :: e_cutoff => NULL()
      TYPE(mulliken_restraint_type), &
         POINTER                           :: mulliken_restraint_control => NULL()
      TYPE(ddapc_restraint_type), &
         DIMENSION(:), POINTER             :: ddapc_restraint_control => NULL()
      TYPE(cdft_control_type), POINTER     :: cdft_control => NULL()
      TYPE(s2_restraint_type), POINTER     :: s2_restraint_control => NULL()
      TYPE(dftb_control_type), POINTER     :: dftb_control => NULL()
      TYPE(xtb_control_type), POINTER      :: xtb_control => NULL()
      TYPE(semi_empirical_control_type), &
         POINTER                           :: se_control => NULL()
      TYPE(gapw_control_type), POINTER     :: gapw_control => NULL()
      TYPE(pw_grid_option)                 :: pw_grid_opt = pw_grid_option()
      LOGICAL                              :: skip_load_balance_distributed = .FALSE.
      ! Types of subsystems for embedding
      LOGICAL                              :: ref_embed_subsys = .FALSE.
      LOGICAL                              :: cluster_embed_subsys = .FALSE.
      LOGICAL                              :: high_level_embed_subsys = .FALSE.
      LOGICAL                              :: dfet_embedded = .FALSE.
      LOGICAL                              :: dmfet_embedded = .FALSE.
   END TYPE qs_control_type

! **************************************************************************************************
! \brief Control parameters for the SCCS models
! **************************************************************************************************
   TYPE sccs_control_type
      LOGICAL                              :: sccs_activated = .FALSE.
      INTEGER                              :: derivative_method = 0, &
                                              max_iter = 0, &
                                              method_id = 0
      REAL(KIND=dp)                        :: alpha_solvent = 0.0_dp, &
                                              beta = 0.0_dp, &
                                              beta_solvent = 0.0_dp, &
                                              delta_rho = 0.0_dp, &
                                              eps_sccs = 0.0_dp, &
                                              eps_scf = 0.0_dp, &
                                              epsilon_solvent = 0.0_dp, &
                                              gamma_solvent = 0.0_dp, &
                                              mixing = 0.0_dp, &
                                              rho_zero = 0.0_dp, &
                                              rho_max = 0.0_dp, &
                                              rho_min = 0.0_dp
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
!   The lumos are helpful in choosing a initial vector for the TDDFPT
!   calculation, since they can be used to construct the solutions of the
!   TDDFPT operator without the perturbation kernel.
! **************************************************************************************************
   TYPE tddfpt_control_type
      TYPE(cp_fm_type), DIMENSION(:), &
         POINTER                           :: lumos => NULL()
      REAL(KIND=dp)                        :: tolerance = 0.0_dp
      INTEGER                              :: n_ev = 0
      INTEGER                              :: max_kv = 0
      INTEGER                              :: n_restarts = 0
      INTEGER                              :: n_reortho = 0
      LOGICAL                              :: do_kernel = .FALSE.
      LOGICAL                              :: lsd_singlets = .FALSE.
      LOGICAL                              :: invert_S = .FALSE.
      LOGICAL                              :: precond = .FALSE.
      LOGICAL                              :: drho_by_collocation = .FALSE.
      LOGICAL                              :: use_kinetic_energy_density = .FALSE.
      INTEGER                              :: res_etype = 0
      INTEGER                              :: diag_method = 0
      INTEGER                              :: oe_corr = 0
      INTEGER                              :: sic_method_id = 0
      INTEGER                              :: sic_list_id = 0
      REAL(KIND=dp)                        :: sic_scaling_a = 0.0_dp, sic_scaling_b = 0.0_dp
      REAL(KIND=dp), DIMENSION(:, :), &
         POINTER                           :: lumos_eigenvalues => NULL()
   END TYPE tddfpt_control_type

! **************************************************************************************************
! \brief Control parameters for simplified Tamm Dancoff approximation (sTDA)
! \par  ATTRIBUTES
! \par  NOTES
! **************************************************************************************************
   TYPE stda_control_type
      LOGICAL                                :: do_ewald = .FALSE.
      LOGICAL                                :: do_exchange = .FALSE.
      REAL(KIND=dp)                          :: hfx_fraction = 0.0_dp
      REAL(KIND=dp)                          :: eps_td_filter = 0.0_dp
      REAL(KIND=dp)                          :: mn_alpha = 0.0_dp
      REAL(KIND=dp)                          :: mn_beta = 0.0_dp
      REAL(KIND=dp)                          :: coulomb_sr_cut = 0.0_dp
      REAL(KIND=dp)                          :: coulomb_sr_eps = 0.0_dp
   END TYPE stda_control_type

! **************************************************************************************************
! \brief Control parameters for a Time-Dependent DFT calculation.
! **************************************************************************************************
   TYPE tddfpt2_control_type
      !> compute TDDFPT excitation energies and oscillator strengths
      LOGICAL                              :: enabled = .FALSE.
      !> number of excited states to converge
      INTEGER                              :: nstates = 0
      !> maximal number of iterations to be performed
      INTEGER                              :: niters = 0
      !> maximal number of Krylov space vectors
      INTEGER                              :: nkvs = 0
      !> number of unoccupied (virtual) molecular orbitals to consider
      INTEGER                              :: nlumo = 0
      !> minimal number of MPI processes to be used per excited state
      INTEGER                              :: nprocs = 0
      !> type of kernel function/approximation to use
      INTEGER                              :: kernel = 0
      !> for full kernel, do we have HFX/ADMM
      LOGICAL                              :: do_hfx = .FALSE.
      LOGICAL                              :: do_admm = .FALSE.
      !> for full kernel, do we have long-range HFX and/or Kxc potential
      LOGICAL                              :: do_hfxlr = .FALSE.
      LOGICAL                              :: do_exck = .FALSE.
      !> options used in sTDA calculation (Kernel)
      TYPE(stda_control_type)              :: stda_control = stda_control_type()
      !> algorithm to correct orbital energies
      INTEGER                              :: oe_corr = 0
      !> eigenvalue shifts
      REAL(KIND=dp)                        :: ev_shift = 0.0_dp, eos_shift = 0.0_dp
      !> target accuracy
      REAL(kind=dp)                        :: conv = 0.0_dp
      !> the smallest excitation amplitude to print
      REAL(kind=dp)                        :: min_excitation_amplitude = 0.0_dp
      !> threshold which controls when two wave functions considered to be orthogonal:
      !> maxabs(Ci^T * S * Cj) <= orthogonal_eps
      REAL(kind=dp)                        :: orthogonal_eps = 0.0_dp
      !> read guess wave functions from restart file if exists
      LOGICAL                              :: is_restart = .FALSE.
      !> compute triplet excited states using spin-unpolarised molecular orbitals
      LOGICAL                              :: rks_triplets = .FALSE.
      !> local resolution of identity for Coulomb contribution
      LOGICAL                              :: do_lrigpw = .FALSE.
      ! automatic generation of auxiliary basis for LRI-TDDFT
      INTEGER                              :: auto_basis_p_lri_aux = 1
      !> use symmetric definition of ADMM Kernel correction
      LOGICAL                              :: admm_symm = .FALSE.
      !> Use/Ignore possible ADMM Kernel XC correction
      LOGICAL                              :: admm_xc_correction = .FALSE.
      !
      ! DIPOLE_MOMENTS subsection
      !
      ! form of the dipole operator used to compute oscillator strengths
      INTEGER                              :: dipole_form = 0
      !> type of the reference point used for calculation of electrostatic dipole moments
      INTEGER                              :: dipole_reference = 0
      !> user-defined reference point
      REAL(kind=dp), DIMENSION(:), POINTER :: dipole_ref_point => NULL()
      !
      ! SOC subsection
      LOGICAL                              :: do_soc = .FALSE.
      !
      ! MGRID subsection
      !
      !> number of plain-wave grids
      INTEGER                              :: mgrid_ngrids = 0
      !> create commensurate grids (progression factor and cutoff values of sub-grids will be ignored)
      LOGICAL                              :: mgrid_commensurate_mgrids = .FALSE.
      !> signals that MGRID section has been explicitly given. Other mgrid_* variables
      !> are not initialised when it is equal to .FALSE. as in this case the default
      !> set of plain-wave grids will be used
      LOGICAL                              :: mgrid_is_explicit = .FALSE.
      !> same as qs_control%realspace_mgrids
      LOGICAL                              :: mgrid_realspace_mgrids = .FALSE.
      !> do not perform load balancing
      LOGICAL                              :: mgrid_skip_load_balance = .FALSE.
      !> cutoff value at the finest grid level
      REAL(kind=dp)                        :: mgrid_cutoff = 0.0_dp
      !> cutoff at the next grid level will be smaller then the cutoff
      !> at the current grid by this number of times
      REAL(kind=dp)                        :: mgrid_progression_factor = 0.0_dp
      !> cutoff that determines to which grid a particular Gaussian function will be mapped
      REAL(kind=dp)                        :: mgrid_relative_cutoff = 0.0_dp
      !> manually provided the list of cutoff values for each grid level
      !> (when it is null(), the cutoff values will be assigned automatically)
      REAL(kind=dp), DIMENSION(:), POINTER :: mgrid_e_cutoff => NULL()
   END TYPE tddfpt2_control_type

! **************************************************************************************************
! \brief Control parameters for a DFT calculation
! \par History
!      10.2019 added variables related to surface dipole correction [Soumya Ghosh]
! **************************************************************************************************
   TYPE dft_control_type
      TYPE(admm_control_type), POINTER     :: admm_control => NULL()
      TYPE(period_efield_type), POINTER    :: period_efield => NULL()
      TYPE(qs_control_type), POINTER       :: qs_control => NULL()
      TYPE(rtp_control_type), POINTER      :: rtp_control => NULL()
      TYPE(sccs_control_type), POINTER     :: sccs_control => NULL()
      TYPE(tddfpt_control_type), POINTER   :: tddfpt_control => NULL()
      TYPE(tddfpt2_control_type), POINTER  :: tddfpt2_control => NULL()
      TYPE(xas_control_type), POINTER      :: xas_control => NULL()
      TYPE(expot_control_type), POINTER    :: expot_control => NULL()
      TYPE(maxwell_control_type), POINTER  :: maxwell_control => NULL()
      TYPE(efield_p_type), POINTER, &
         DIMENSION(:)                      :: efield_fields => NULL()
      INTEGER                              :: nspins = 0, &
                                              charge = 0, &
                                              multiplicity = 0, &
                                              sic_method_id = 0, &
                                              plus_u_method_id = 0, &
                                              dir_surf_dip = 0, &
                                              nimages = 1
      INTEGER                              :: sic_list_id = 0
      INTEGER                              :: auto_basis_ri_aux = 1, &
                                              auto_basis_aux_fit = 1, &
                                              auto_basis_lri_aux = 1, &
                                              auto_basis_p_lri_aux = 1, &
                                              auto_basis_ri_hxc = 1, &
                                              auto_basis_ri_xas = 1, &
                                              auto_basis_ri_hfx = 1
      REAL(KIND=dp)                        :: relax_multiplicity = 0.0_dp, &
                                              sic_scaling_a = 0.0_dp, &
                                              sic_scaling_b = 0.0_dp, &
                                              pos_dir_surf_dip = 0.0_dp
      LOGICAL                              :: do_tddfpt_calculation = .FALSE., &
                                              do_xas_calculation = .FALSE., &
                                              do_xas_tdp_calculation = .FALSE., &
                                              drho_by_collocation = .FALSE., &
                                              use_kinetic_energy_density = .FALSE., &
                                              restricted = .FALSE., &
                                              roks = .FALSE., &
                                              uks = .FALSE., &
                                              lsd = .FALSE., &
                                              dft_plus_u = .FALSE., &
                                              apply_efield = .FALSE., &
                                              apply_efield_field = .FALSE., &
                                              apply_vector_potential = .FALSE., &
                                              apply_period_efield = .FALSE., &
                                              apply_external_potential = .FALSE., &
                                              eval_external_potential = .FALSE., &
                                              do_admm = .FALSE., &
                                              do_admm_dm = .FALSE., &
                                              do_admm_mo = .FALSE., &
                                              smear = .FALSE., &
                                              low_spin_roks = .FALSE., &
                                              apply_external_density = .FALSE., &
                                              read_external_density = .FALSE., &
                                              apply_external_vxc = .FALSE., &
                                              read_external_vxc = .FALSE., &
                                              correct_surf_dip = .FALSE., &
                                              surf_dip_correct_switch = .FALSE., &
                                              switch_surf_dip = .FALSE., &
                                              correct_el_density_dip = .FALSE., &
                                              do_sccs = .FALSE., &
                                              apply_embed_pot = .FALSE., &
                                              apply_dmfet_pot = .FALSE.
   END TYPE dft_control_type

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_control_types'

   ! Public data types

   PUBLIC :: dft_control_type, &
             qs_control_type, &
             gapw_control_type, &
             tddfpt_control_type, &
             tddfpt2_control_type, &
             proj_mo_type, &
             efield_type, &
             mulliken_restraint_type, &
             ddapc_restraint_type, &
             dftb_control_type, &
             xtb_control_type, &
             semi_empirical_control_type, &
             s2_restraint_type, &
             admm_control_type, &
             maxwell_control_type, &
             expot_control_type, &
             rtp_control_type, &
             sccs_control_type, &
             stda_control_type

   ! Public subroutines

   PUBLIC :: dft_control_release, &
             dft_control_create, &
             tddfpt_control_create, &
             admm_control_create, &
             maxwell_control_create, &
             expot_control_create, &
             ddapc_control_create

CONTAINS

! **************************************************************************************************
!> \brief create  the mulliken_restraint_type
!> \param mulliken_restraint_control ...
!> \par History
!>      02.2005 created [Joost VandeVondele]
! **************************************************************************************************
   SUBROUTINE mulliken_control_create(mulliken_restraint_control)
      TYPE(mulliken_restraint_type), INTENT(OUT)         :: mulliken_restraint_control

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
      TYPE(mulliken_restraint_type), INTENT(INOUT)       :: mulliken_restraint_control

      IF (ASSOCIATED(mulliken_restraint_control%atoms)) &
         DEALLOCATE (mulliken_restraint_control%atoms)
      mulliken_restraint_control%strength = 0.0_dp
      mulliken_restraint_control%target = 0.0_dp
      mulliken_restraint_control%natoms = 0
   END SUBROUTINE mulliken_control_release

! **************************************************************************************************
!> \brief create the ddapc_restraint_type
!> \param ddapc_restraint_control ...
!> \par History
!>      02.2006 created [Joost VandeVondele]
! **************************************************************************************************
   SUBROUTINE ddapc_control_create(ddapc_restraint_control)
      TYPE(ddapc_restraint_type), INTENT(OUT)            :: ddapc_restraint_control

      ddapc_restraint_control%density_type = do_full_density
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
      TYPE(ddapc_restraint_type), INTENT(INOUT)          :: ddapc_restraint_control

      IF (ASSOCIATED(ddapc_restraint_control%atoms)) &
         DEALLOCATE (ddapc_restraint_control%atoms)
      IF (ASSOCIATED(ddapc_restraint_control%coeff)) &
         DEALLOCATE (ddapc_restraint_control%coeff)
      ddapc_restraint_control%strength = 0.0_dp
      ddapc_restraint_control%target = 0.0_dp
      ddapc_restraint_control%natoms = 0
   END SUBROUTINE ddapc_control_release

! **************************************************************************************************
!> \brief create the s2_restraint_type
!> \param s2_restraint_control ...
!> \par History
!>      03.2006 created [Joost VandeVondele]
! **************************************************************************************************
   SUBROUTINE s2_control_create(s2_restraint_control)
      TYPE(s2_restraint_type), INTENT(OUT)               :: s2_restraint_control

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
      TYPE(s2_restraint_type), INTENT(INOUT)             :: s2_restraint_control

      s2_restraint_control%strength = 0.0_dp
      s2_restraint_control%target = 0.0_dp
   END SUBROUTINE s2_control_release

! **************************************************************************************************
!> \brief allocates and perform a very basic initialization
!> \param dft_control the object to create
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE dft_control_create(dft_control)
      TYPE(dft_control_type), INTENT(OUT)                :: dft_control

      NULLIFY (dft_control%xas_control)
      NULLIFY (dft_control%qs_control)
      NULLIFY (dft_control%tddfpt_control)
      NULLIFY (dft_control%tddfpt2_control)
      NULLIFY (dft_control%efield_fields)
      NULLIFY (dft_control%period_efield)
      NULLIFY (dft_control%admm_control)
      NULLIFY (dft_control%expot_control)
      NULLIFY (dft_control%maxwell_control)
      NULLIFY (dft_control%rtp_control)
      NULLIFY (dft_control%sccs_control)
      dft_control%do_sccs = .FALSE.
      dft_control%apply_embed_pot = .FALSE.
      dft_control%apply_dmfet_pot = .FALSE.
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
   SUBROUTINE dft_control_release(dft_control)
      TYPE(dft_control_type), INTENT(INOUT)              :: dft_control

      CALL qs_control_release(dft_control%qs_control)
      CALL tddfpt_control_release(dft_control%tddfpt_control)
      CALL tddfpt2_control_release(dft_control%tddfpt2_control)
      IF (ASSOCIATED(dft_control%xas_control)) THEN
         CALL xas_control_release(dft_control%xas_control)
         DEALLOCATE (dft_control%xas_control)
      END IF
      CALL admm_control_release(dft_control%admm_control)
      CALL expot_control_release(dft_control%expot_control)
      CALL maxwell_control_release(dft_control%maxwell_control)
      CALL efield_fields_release(dft_control%efield_fields)
      IF (ASSOCIATED(dft_control%sccs_control)) DEALLOCATE (dft_control%sccs_control)
      IF (ASSOCIATED(dft_control%period_efield)) THEN
         DEALLOCATE (dft_control%period_efield)
      END IF
      IF (ASSOCIATED(dft_control%rtp_control)) THEN
         CALL proj_mo_list_release(dft_control%rtp_control%proj_mo_list)
         DEALLOCATE (dft_control%rtp_control)
      END IF

   END SUBROUTINE dft_control_release

! **************************************************************************************************
!> \brief ...
!> \param qs_control ...
! **************************************************************************************************
   SUBROUTINE qs_control_create(qs_control)
      TYPE(qs_control_type), POINTER                     :: qs_control

      CPASSERT(.NOT. ASSOCIATED(qs_control))
      ALLOCATE (qs_control)

      NULLIFY (qs_control%e_cutoff)
      NULLIFY (qs_control%gapw_control)
      NULLIFY (qs_control%mulliken_restraint_control)
      NULLIFY (qs_control%ddapc_restraint_control)
      NULLIFY (qs_control%s2_restraint_control)
      NULLIFY (qs_control%se_control)
      NULLIFY (qs_control%dftb_control)
      NULLIFY (qs_control%xtb_control)
      NULLIFY (qs_control%cdft_control)
      NULLIFY (qs_control%ddapc_restraint_control)

      ALLOCATE (qs_control%mulliken_restraint_control)
      CALL mulliken_control_create(qs_control%mulliken_restraint_control)
      ALLOCATE (qs_control%s2_restraint_control)
      CALL s2_control_create(qs_control%s2_restraint_control)
      ALLOCATE (qs_control%gapw_control)
      CALL se_control_create(qs_control%se_control)
      CALL dftb_control_create(qs_control%dftb_control)
      CALL xtb_control_create(qs_control%xtb_control)
      ALLOCATE (qs_control%cdft_control)
      CALL cdft_control_create(qs_control%cdft_control)
   END SUBROUTINE qs_control_create

! **************************************************************************************************
!> \brief ...
!> \param qs_control ...
! **************************************************************************************************
   SUBROUTINE qs_control_release(qs_control)
      TYPE(qs_control_type), POINTER                     :: qs_control

      INTEGER                                            :: i

      IF (ASSOCIATED(qs_control)) THEN
         CALL mulliken_control_release(qs_control%mulliken_restraint_control)
         DEALLOCATE (qs_control%mulliken_restraint_control)
         CALL s2_control_release(qs_control%s2_restraint_control)
         DEALLOCATE (qs_control%s2_restraint_control)
         CALL se_control_release(qs_control%se_control)
         CALL dftb_control_release(qs_control%dftb_control)
         CALL xtb_control_release(qs_control%xtb_control)
         IF (ASSOCIATED(qs_control%cdft_control)) THEN
            CALL cdft_control_release(qs_control%cdft_control)
            DEALLOCATE (qs_control%cdft_control)
         END IF

         IF (ASSOCIATED(qs_control%e_cutoff)) THEN
            DEALLOCATE (qs_control%e_cutoff)
         END IF
         IF (ASSOCIATED(qs_control%gapw_control)) THEN
            DEALLOCATE (qs_control%gapw_control)
         END IF
         IF (ASSOCIATED(qs_control%ddapc_restraint_control)) THEN
            DO i = 1, SIZE(qs_control%ddapc_restraint_control)
               CALL ddapc_control_release(qs_control%ddapc_restraint_control(i))
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

      IF (ASSOCIATED(tddfpt_control)) THEN
         CALL cp_fm_release(tddfpt_control%lumos)
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

      CHARACTER(len=*), PARAMETER :: routineN = 'tddfpt2_control_create'

      INTEGER                                            :: handle

      CPASSERT(.NOT. ASSOCIATED(tddfpt_control))
      CALL timeset(routineN, handle)

      ALLOCATE (tddfpt_control)
      tddfpt_control%do_soc = .FALSE.

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

      CHARACTER(len=*), PARAMETER :: routineN = 'tddfpt2_control_release'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (ASSOCIATED(tddfpt_control)) THEN
         DEALLOCATE (tddfpt_control)
      END IF

      CALL timestop(handle)
   END SUBROUTINE tddfpt2_control_release

! **************************************************************************************************
!> \brief ...
!> \param proj_mo_list ...
! **************************************************************************************************
   SUBROUTINE proj_mo_list_release(proj_mo_list)
      TYPE(proj_mo_p_type), DIMENSION(:), POINTER        :: proj_mo_list

      INTEGER                                            :: i, mo_ref_nbr

      IF (ASSOCIATED(proj_mo_list)) THEN
         DO i = 1, SIZE(proj_mo_list)
            IF (ASSOCIATED(proj_mo_list(i)%proj_mo)) THEN
               IF (ALLOCATED(proj_mo_list(i)%proj_mo%ref_mo_index)) &
                  DEALLOCATE (proj_mo_list(i)%proj_mo%ref_mo_index)
               IF (ALLOCATED(proj_mo_list(i)%proj_mo%mo_ref)) THEN
                  DO mo_ref_nbr = 1, SIZE(proj_mo_list(i)%proj_mo%mo_ref)
                     CALL cp_fm_release(proj_mo_list(i)%proj_mo%mo_ref(mo_ref_nbr))
                  END DO
                  DEALLOCATE (proj_mo_list(i)%proj_mo%mo_ref)
               END IF
               IF (ALLOCATED(proj_mo_list(i)%proj_mo%td_mo_index)) &
                  DEALLOCATE (proj_mo_list(i)%proj_mo%td_mo_index)
               DEALLOCATE (proj_mo_list(i)%proj_mo)
            END IF
         END DO
         DEALLOCATE (proj_mo_list)
      END IF
   END SUBROUTINE proj_mo_list_release

! **************************************************************************************************
!> \brief ...
!> \param efield_fields ...
! **************************************************************************************************
   SUBROUTINE efield_fields_release(efield_fields)
      TYPE(efield_p_type), DIMENSION(:), POINTER         :: efield_fields

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
               IF (ASSOCIATED(efield_fields(i)%efield%polarisation)) &
                  DEALLOCATE (efield_fields(i)%efield%polarisation)
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

      IF (ASSOCIATED(dftb_control)) THEN
         IF (ASSOCIATED(dftb_control%sk_pair_list)) THEN
            DEALLOCATE (dftb_control%sk_pair_list)
         END IF
         DEALLOCATE (dftb_control)
      END IF
   END SUBROUTINE dftb_control_release

! **************************************************************************************************
!> \brief ...
!> \param xtb_control ...
! **************************************************************************************************
   SUBROUTINE xtb_control_create(xtb_control)
      TYPE(xtb_control_type), POINTER                    :: xtb_control

      CPASSERT(.NOT. ASSOCIATED(xtb_control))
      ALLOCATE (xtb_control)

      NULLIFY (xtb_control%kab_param)
      NULLIFY (xtb_control%kab_vals)
      NULLIFY (xtb_control%kab_types)
      NULLIFY (xtb_control%nonbonded)

   END SUBROUTINE xtb_control_create

! **************************************************************************************************
!> \brief ...
!> \param xtb_control ...
! **************************************************************************************************
   SUBROUTINE xtb_control_release(xtb_control)
      TYPE(xtb_control_type), POINTER                    :: xtb_control

      IF (ASSOCIATED(xtb_control)) THEN
         IF (ASSOCIATED(xtb_control%kab_param)) THEN
            DEALLOCATE (xtb_control%kab_param)
         END IF
         IF (ASSOCIATED(xtb_control%kab_vals)) THEN
            DEALLOCATE (xtb_control%kab_vals)
         END IF
         IF (ASSOCIATED(xtb_control%kab_types)) THEN
            DEALLOCATE (xtb_control%kab_types)
         END IF
         IF (ASSOCIATED(xtb_control%nonbonded)) THEN
            CALL pair_potential_p_release(xtb_control%nonbonded)
         END IF
         DEALLOCATE (xtb_control)
      END IF
   END SUBROUTINE xtb_control_release

! **************************************************************************************************
!> \brief ...
!> \param se_control ...
! **************************************************************************************************
   SUBROUTINE se_control_create(se_control)
      TYPE(semi_empirical_control_type), POINTER         :: se_control

      CPASSERT(.NOT. ASSOCIATED(se_control))
      ALLOCATE (se_control)
   END SUBROUTINE se_control_create

! **************************************************************************************************
!> \brief ...
!> \param se_control ...
! **************************************************************************************************
   SUBROUTINE se_control_release(se_control)
      TYPE(semi_empirical_control_type), POINTER         :: se_control

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

      CPASSERT(.NOT. ASSOCIATED(admm_control))
      ALLOCATE (admm_control)

   END SUBROUTINE admm_control_create

! **************************************************************************************************
!> \brief ...
!> \param admm_control ...
! **************************************************************************************************
   SUBROUTINE admm_control_release(admm_control)
      TYPE(admm_control_type), POINTER                   :: admm_control

      IF (ASSOCIATED(admm_control)) THEN
         DEALLOCATE (admm_control)
      END IF
   END SUBROUTINE admm_control_release

! **************************************************************************************************
!> \brief ...
!> \param expot_control ...
! **************************************************************************************************
   SUBROUTINE expot_control_create(expot_control)
      TYPE(expot_control_type), POINTER                  :: expot_control

      CPASSERT(.NOT. ASSOCIATED(expot_control))
      ALLOCATE (expot_control)
      expot_control%read_from_cube = .FALSE.
      expot_control%maxwell_solver = .FALSE.
      expot_control%static = .TRUE.
      expot_control%scaling_factor = 1.0_dp

   END SUBROUTINE expot_control_create

! **************************************************************************************************
!> \brief ...
!> \param expot_control ...
! **************************************************************************************************
   SUBROUTINE expot_control_release(expot_control)
      TYPE(expot_control_type), POINTER                  :: expot_control

      IF (ASSOCIATED(expot_control)) THEN
         DEALLOCATE (expot_control)
      END IF

   END SUBROUTINE expot_control_release

! **************************************************************************************************
!> \brief ...
!> \param maxwell_control ...
! **************************************************************************************************
   SUBROUTINE maxwell_control_create(maxwell_control)
      TYPE(maxwell_control_type), POINTER                :: maxwell_control

      CPASSERT(.NOT. ASSOCIATED(maxwell_control))
      ALLOCATE (maxwell_control)

   END SUBROUTINE maxwell_control_create

! **************************************************************************************************
!> \brief ...
!> \param maxwell_control ...
! **************************************************************************************************
   SUBROUTINE maxwell_control_release(maxwell_control)
      TYPE(maxwell_control_type), POINTER                :: maxwell_control

      IF (ASSOCIATED(maxwell_control)) THEN
         DEALLOCATE (maxwell_control)
      END IF

   END SUBROUTINE maxwell_control_release

END MODULE cp_control_types
