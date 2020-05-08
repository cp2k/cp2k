!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Defines control structures, which contain the parameters and the
!>      settings for the DFT-based calculations.
! **************************************************************************************************
MODULE cp_control_types
   USE cp_fm_types,                     ONLY: cp_fm_p_type,&
                                              cp_fm_release
   USE input_constants,                 ONLY: do_full_density
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
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
! \brief Control parameters for xTB calculations
! **************************************************************************************************
   TYPE xtb_control_type
      !
      LOGICAL                              :: do_ewald
      !
      INTEGER                              :: sto_ng
      INTEGER                              :: h_sto_ng
      !
      CHARACTER(LEN=default_path_length)   :: parameter_file_path
      CHARACTER(LEN=default_path_length)   :: parameter_file_name
      !
      CHARACTER(LEN=default_path_length)   :: dispersion_parameter_file
      REAL(KIND=dp)                        :: epscn
      REAL(KIND=dp)                        :: rcdisp
      REAL(KIND=dp)                        :: s6, s8
      REAL(KIND=dp)                        :: a1, a2
      !
      REAL(KIND=dp)                        :: ks, kp, kd, ksp, k2sh
      REAL(KIND=dp)                        :: kg, kf
      REAL(KIND=dp)                        :: kcns, kcnp, kcnd
      REAL(KIND=dp)                        :: ken
      REAL(KIND=dp)                        :: kxr, kx2
      !
      LOGICAL                              :: xb_interaction
      LOGICAL                              :: coulomb_interaction
      LOGICAL                              :: tb3_interaction
      LOGICAL                              :: check_atomic_charges
      !
      REAL(KIND=dp)                        :: xb_radius
      !
      CHARACTER(LEN=default_string_length), &
         DIMENSION(:, :), POINTER          :: kab_param
      INTEGER, DIMENSION(:, :), POINTER    :: kab_types
      INTEGER                              :: kab_nval
      REAL, DIMENSION(:), POINTER          :: kab_vals
   END TYPE xtb_control_type

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
      LOGICAL                              :: gapw, gapw_xc, gpw, pao
      LOGICAL                              :: lrigpw, rigpw
      LOGICAL                              :: lri_optbas
      LOGICAL                              :: ofgpw
      LOGICAL                              :: dftb
      LOGICAL                              :: xtb
      LOGICAL                              :: semi_empirical
      LOGICAL                              :: mulliken_restraint
      LOGICAL                              :: ddapc_restraint
      LOGICAL                              :: ddapc_restraint_is_spin
      LOGICAL                              :: ddapc_explicit_potential
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
      TYPE(cdft_control_type), POINTER     :: cdft_control
      TYPE(s2_restraint_type), POINTER     :: s2_restraint_control
      TYPE(dftb_control_type), POINTER     :: dftb_control
      TYPE(xtb_control_type), POINTER      :: xtb_control
      TYPE(semi_empirical_control_type), &
         POINTER                           :: se_control
      TYPE(gapw_control_type), POINTER     :: gapw_control
      TYPE(pw_grid_option)                 :: pw_grid_opt
      LOGICAL                              :: skip_load_balance_distributed
      ! Types of subsystems for embedding
      LOGICAL                              :: ref_embed_subsys
      LOGICAL                              :: cluster_embed_subsys
      LOGICAL                              :: high_level_embed_subsys
      LOGICAL                              :: dfet_embedded
      LOGICAL                              :: dmfet_embedded
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
!   The lumos are helpful in choosing a initial vector for the TDDFPT
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
! \brief Control parameters for simplified Tamm Dancoff approximation (sTDA)
! \par  ATTRIBUTES
! \par  NOTES
! **************************************************************************************************
   TYPE stda_control_type
      LOGICAL                                :: do_ewald
      REAL(KIND=dp)                          :: hfx_fraction
      REAL(KIND=dp)                          :: eps_td_filter
   END TYPE stda_control_type

! **************************************************************************************************
! \brief Control parameters for a Time-Dependent DFT calculation.
! **************************************************************************************************
   TYPE tddfpt2_control_type
      !> compute TDDFPT excitation energies and oscillator strengths
      LOGICAL                              :: enabled
      !> number of excited states to converge
      INTEGER                              :: nstates
      !> maximal number of iterations to be performed
      INTEGER                              :: niters
      !> maximal number of Krylov space vectors
      INTEGER                              :: nkvs
      !> number of unoccupied (virtual) molecular orbitals to consider
      INTEGER                              :: nlumo
      !> minimal number of MPI processes to be used per excited state
      INTEGER                              :: nprocs
      !> type of kernel function/approximation to use
      INTEGER                              :: kernel
      !> options used in sTDA calculation (Kernel)
      TYPE(stda_control_type)              :: stda_control
      !> algorithm to correct orbital energies
      INTEGER                              :: oe_corr
      !> eigenvalue shifts
      REAL(KIND=dp)                        :: ev_shift, eos_shift
      !> target accuracy
      REAL(kind=dp)                        :: conv
      !> the smallest excitation amplitude to print
      REAL(kind=dp)                        :: min_excitation_amplitude
      !> threshold which controls when two wave functions considered to be orthogonal:
      !> maxabs(Ci^T * S * Cj) <= orthogonal_eps
      REAL(kind=dp)                        :: orthogonal_eps
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
! \par History
!      10.2019 added variables related to surface dipole correction [Soumya Ghosh]
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
      INTEGER                              :: auto_basis_ri_aux = 1, &
                                              auto_basis_aux_fit = 1, &
                                              auto_basis_lri_aux = 1, &
                                              auto_basis_ri_hxc = 1, &
                                              auto_basis_ri_xas = 1, &
                                              auto_basis_ri_hfx = 1
      REAL(KIND=dp)                        :: relax_multiplicity, &
                                              sic_scaling_a, &
                                              sic_scaling_b, &
                                              pos_dir_surf_dip
      LOGICAL                              :: do_tddfpt_calculation, &
                                              do_xas_calculation, &
                                              do_xas_tdp_calculation, &
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
                                              surf_dip_correct_switch, &
                                              switch_surf_dip, &
                                              correct_el_density_dip, &
                                              do_sccs, &
                                              apply_embed_pot, &
                                              apply_dmfet_pot
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
             xtb_control_type, &
             semi_empirical_control_type, &
             s2_restraint_type, &
             admm_control_type, &
             rtp_control_type, &
             sccs_control_type, &
             stda_control_type

   ! Public subroutines

   PUBLIC :: dft_control_retain, &
             dft_control_release, &
             dft_control_create, &
             tddfpt_control_create, &
             admm_control_create, &
             ddapc_control_create, &
             sccs_control_create

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
      mulliken_restraint_control%ref_count = mulliken_restraint_control%ref_count - 1
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

      mulliken_restraint_control%ref_count = mulliken_restraint_control%ref_count + 1
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
      ddapc_restraint_control%ref_count = ddapc_restraint_control%ref_count - 1
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

      ddapc_restraint_control%ref_count = ddapc_restraint_control%ref_count + 1
   END SUBROUTINE ddapc_control_retain

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
      s2_restraint_control%ref_count = s2_restraint_control%ref_count - 1
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
      s2_restraint_control%ref_count = s2_restraint_control%ref_count + 1
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
      last_dft_control_id = last_dft_control_id + 1
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
   SUBROUTINE dft_control_retain(dft_control)
      TYPE(dft_control_type), POINTER                    :: dft_control

      CHARACTER(len=*), PARAMETER :: routineN = 'dft_control_retain', &
         routineP = moduleN//':'//routineN

      CPASSERT(ASSOCIATED(dft_control))
      CPASSERT(dft_control%ref_count > 0)
      dft_control%ref_count = dft_control%ref_count + 1
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
         dft_control%ref_count = dft_control%ref_count - 1
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
      NULLIFY (qs_control%xtb_control)
      NULLIFY (qs_control%cdft_control)
      NULLIFY (qs_control%ddapc_restraint_control)

      CALL mulliken_control_create(qs_control%mulliken_restraint_control)
      CALL s2_control_create(qs_control%s2_restraint_control)
      CALL gapw_control_create(qs_control%gapw_control)
      CALL se_control_create(qs_control%se_control)
      CALL dftb_control_create(qs_control%dftb_control)
      CALL xtb_control_create(qs_control%xtb_control)
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
         CALL xtb_control_release(qs_control%xtb_control)
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
!> \param xtb_control ...
! **************************************************************************************************
   SUBROUTINE xtb_control_create(xtb_control)
      TYPE(xtb_control_type), POINTER                    :: xtb_control

      CHARACTER(len=*), PARAMETER :: routineN = 'xtb_control_create', &
         routineP = moduleN//':'//routineN

      CPASSERT(.NOT. ASSOCIATED(xtb_control))
      ALLOCATE (xtb_control)

      NULLIFY (xtb_control%kab_param)
      NULLIFY (xtb_control%kab_vals)
      NULLIFY (xtb_control%kab_types)

   END SUBROUTINE xtb_control_create

! **************************************************************************************************
!> \brief ...
!> \param xtb_control ...
! **************************************************************************************************
   SUBROUTINE xtb_control_release(xtb_control)
      TYPE(xtb_control_type), POINTER                    :: xtb_control

      CHARACTER(len=*), PARAMETER :: routineN = 'xtb_control_release', &
         routineP = moduleN//':'//routineN

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
         DEALLOCATE (xtb_control)
      END IF
   END SUBROUTINE xtb_control_release

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
         sccs_control%ref_count = sccs_control%ref_count - 1
         IF (sccs_control%ref_count == 0) THEN
            DEALLOCATE (sccs_control)
         END IF
      END IF

   END SUBROUTINE sccs_control_release

END MODULE cp_control_types
