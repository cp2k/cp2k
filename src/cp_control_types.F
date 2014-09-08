!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Defines control structures, which contain the parameters and the
!>      settings for the DFT-based calculations.
! *****************************************************************************
MODULE cp_control_types
  USE cp_fm_types,                     ONLY: cp_fm_p_type,&
                                             cp_fm_release
  USE f77_blas
  USE input_constants,                 ONLY: do_full_density
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE pw_types,                        ONLY: pw_p_type
  USE xas_control,                     ONLY: xas_control_create,&
                                             xas_control_release,&
                                             xas_control_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *****************************************************************************
! \brief Control parameters for pw grids
! *****************************************************************************
  TYPE pw_grid_option
     LOGICAL :: spherical
     LOGICAL :: fullspace
     INTEGER, DIMENSION(2) :: distribution_layout
     INTEGER :: blocked
  END TYPE pw_grid_option

! *****************************************************************************
! \brief Control parameters for REAL_TIME_PROPAGATION calculations
! *****************************************************************************
  TYPE rtp_control_type
     LOGICAL                              :: converged
     REAL(KIND = dp)                      :: eps_ener
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
     LOGICAL                              :: calc_idempotency
     LOGICAL                              :: orthonormal
     INTEGER                              :: mcweeny_max_iter
     INTEGER                              :: acc_ref
     REAL(dp)                             :: mcweeny_eps
     INTEGER, DIMENSION(3)                :: delta_pulse_direction
     REAL(KIND=dp)                        :: delta_pulse_scale
  END TYPE rtp_control_type
! *****************************************************************************
! \brief Control parameters for DFTB calculations
! *****************************************************************************
  TYPE dftb_control_type
     LOGICAL                              :: self_consistent
     LOGICAL                              :: orthogonal_basis
     LOGICAL                              :: dispersion
     INTEGER                              :: dispersion_type
     LOGICAL                              :: dftb3_diagonal
     LOGICAL                              :: hb_sr_damp
     REAL(KIND = dp)                      :: hb_sr_para
     REAL(KIND = dp)                      :: eps_disp
     REAL(KIND = dp)                      :: epscn
     REAL(KIND = dp)                      :: rcdisp
     REAL(KIND = dp), DIMENSION(3)        :: sd3
     LOGICAL                              :: do_ewald
     CHARACTER(LEN=default_path_length)   :: sk_file_path
     CHARACTER(LEN=default_path_length)   :: sk_file_list
     CHARACTER(LEN=default_string_length),&
          DIMENSION(:,:), POINTER         :: sk_pair_list
     CHARACTER(LEN=default_path_length)   :: uff_force_field
     CHARACTER(LEN=default_path_length)   :: dispersion_parameter_file
  END TYPE dftb_control_type

! *****************************************************************************
! \brief Control parameters for SCPTB calculations
! *****************************************************************************
  TYPE scptb_control_type
     LOGICAL                              :: dispersion
     LOGICAL                              :: do_ewald
     LOGICAL                              :: do_scc
     LOGICAL                              :: do_scp
     INTEGER                              :: sto_ng
     REAL(KIND = dp)                      :: rcdisp
     REAL(KIND = dp)                      :: epscn
     REAL(KIND = dp), DIMENSION(3)        :: sd3
     CHARACTER(LEN=default_path_length)   :: parameter_file
     CHARACTER(LEN=default_path_length)   :: dispersion_parameter_file
     REAL(KIND = dp)                      :: epspair
  END TYPE scptb_control_type

! *****************************************************************************
! \brief Control parameters for semi empirical calculations
! *****************************************************************************
  TYPE semi_empirical_control_type
     LOGICAL               :: orthogonal_basis
     LOGICAL               :: analytical_gradients
     LOGICAL               :: force_kdsod_EX
     LOGICAL               :: do_ewald, do_ewald_r3, do_ewald_gks
     LOGICAL               :: scp
     INTEGER               :: integral_screening, periodic_type
     INTEGER               :: max_multipole
     INTEGER               :: ga_ncells
     REAL(KIND = dp)       :: delta
     ! Dispersion pair potential
     LOGICAL               :: dispersion
     REAL(KIND = dp)                      :: rcdisp
     REAL(KIND = dp)                      :: epscn
     REAL(KIND = dp), DIMENSION(3)        :: sd3
     CHARACTER(LEN=default_path_length)   :: dispersion_parameter_file
     ! Parameters controlling the evaluation of the integrals
     REAL(KIND = dp)       :: cutoff_lrc, taper_lrc, range_lrc
     REAL(KIND = dp)       :: cutoff_cou, taper_cou, range_cou
     REAL(KIND = dp)       :: cutoff_exc, taper_exc, range_exc
     REAL(KIND = dp)       :: taper_scr , range_scr
  END TYPE semi_empirical_control_type

! *****************************************************************************
! \brief Control parameters for GAPW method within QUICKSTEP ***
! *****************************************************************************
  TYPE gapw_control_type
    REAL(KIND = dp)                 :: eps_fit,&
                                       eps_iso,&
                                       eps_Vrho0,&
                                       eps_svd,&
                                       eps_cpc
    INTEGER                         :: ladd_rho0,&
                                       lmax_rho0,&
                                       lmax_sphere,&
                                       quadrature
    LOGICAL                         :: lrho1_eq_lrho0
    LOGICAL                         :: alpha0_hard_from_input,&
                                       force_paw,&
                                       non_paw_atoms,&
                                       nopaw_as_gpw
    REAL(KIND = dp)                 :: alpha0_hard
    REAL(KIND = dp)                 :: max_rad_local
  END TYPE gapw_control_type
! *****************************************************************************
! \brief parameters for calculations involving a time dependent electric field
! *****************************************************************************
  TYPE efield_type
     REAL(KIND = dp)                 :: actual_time
     REAL(KIND = dp), DIMENSION(:),POINTER :: polarisation
     INTEGER                         :: envelop_id
     REAL(kind=dp),DIMENSION(:),&
                             POINTER :: envelop_r_vars
     INTEGER,DIMENSION(:),&
                             POINTER :: envelop_i_vars
     REAL(kind=dp)                   :: strength
     REAL(KIND=dp)                   :: phase_offset
     REAL(KIND=dp)                   :: wavelength
  END TYPE efield_type

  TYPE efield_p_type
     TYPE(efield_type),POINTER       :: efield
  END TYPE efield_p_type
! *****************************************************************************
! \brief parameters for calculations involving a time dependent electric field
! *****************************************************************************
  TYPE period_efield_type
     REAL(KIND = dp), DIMENSION(:),POINTER :: polarisation
     REAL(kind=dp)                   :: strength
  END TYPE period_efield_type

! *****************************************************************************
! \brief some parameters useful for mulliken_restraints
! *****************************************************************************
  TYPE mulliken_restraint_type
    INTEGER                         :: ref_count
    REAL(KIND = dp)                 :: strength
    REAL(KIND = dp)                 :: TARGET
    INTEGER                         :: natoms
    INTEGER, POINTER, DIMENSION(:)  :: atoms
  END TYPE mulliken_restraint_type

! *****************************************************************************
! \brief some parameters useful for ddapc_restraints
! *****************************************************************************
  TYPE ddapc_restraint_type
    INTEGER                         :: ref_count
    REAL(KIND = dp)                 :: strength
    REAL(KIND = dp)                 :: TARGET
    REAL(KIND = dp)                 :: ddapc_order_p
    INTEGER                         :: functional_form
    INTEGER                         :: natoms
    INTEGER, POINTER, DIMENSION(:)  :: atoms
    REAL(KIND=dp), POINTER, DIMENSION(:)  :: coeff
    INTEGER                         :: density_type
  END TYPE ddapc_restraint_type

! *****************************************************************************
! \brief provides a vector of pointers to ddapc_restraint_type
! *****************************************************************************
  TYPE ddapc_restraint_p_type
     TYPE(ddapc_restraint_type),POINTER:: ddapc_restraint_control
  END TYPE ddapc_restraint_p_type

! *****************************************************************************
! \brief some parameters useful for becke_restraints
! *****************************************************************************
  TYPE becke_restraint_type
    INTEGER                         :: ref_count
    REAL(KIND = dp)                 :: strength
    REAL(KIND = dp)                 :: TARGET
    REAL(KIND = dp)                 :: becke_order_p
    INTEGER                         :: functional_form
    INTEGER                         :: natoms
    INTEGER, POINTER, DIMENSION(:)  :: atoms
    TYPE(pw_p_type)                 :: becke_pot
    LOGICAL                         :: need_pot
    REAL(KIND=dp), POINTER, DIMENSION(:)  :: coeff
    INTEGER                         :: density_type
 END TYPE becke_restraint_type

! *****************************************************************************
! \brief some parameters useful for s2_restraints
! *****************************************************************************
  TYPE s2_restraint_type
    INTEGER                         :: ref_count
    REAL(KIND = dp)                 :: strength
    REAL(KIND = dp)                 :: TARGET
    REAL(KIND = dp)                 :: s2_order_p
    INTEGER                         :: functional_form
  END TYPE s2_restraint_type

! *****************************************************************************
! \brief some parameters useful for auxiliary density matrix method
! *****************************************************************************
  TYPE admm_block_type
    INTEGER, DIMENSION(:), ALLOCATABLE :: list
  END TYPE admm_block_type

  TYPE admm_control_type
    INTEGER                         :: purification_method
    INTEGER                         :: method
    TYPE(admm_block_type), DIMENSION(:), &
                        ALLOCATABLE :: blocks
  END TYPE admm_control_type

! *****************************************************************************
! \brief Control parameters for a QUICKSTEP and KIM-GORDON calculation ***
!        eps_pgf_orb: Cutoff value for the interaction of the primitive
!                     Gaussian-type functions (primitive basis functions).
! *****************************************************************************
  TYPE qs_control_type
    CHARACTER(LEN=10)                             :: method
    INTEGER                                       :: method_id
    REAL(KIND = dp)                               :: eps_core_charge,&
                                                     eps_kg_orb,&
                                                     eps_pgf_orb,&
                                                     eps_ppl,&
                                                     eps_ppnl,&
                                                     eps_rho_gspace,&
                                                     eps_rho_rspace,&
                                                     eps_filter_matrix,&
                                                     eps_gvg_rspace,&
                                                     progression_factor,&
                                                     relative_cutoff
    LOGICAL                                       :: do_almo_scf
    LOGICAL                                       :: do_ls_scf
    LOGICAL                                       :: do_kg
    LOGICAL                                       :: commensurate_mgrids
    LOGICAL                                       :: realspace_mgrids
    LOGICAL                                       :: map_consistent
    LOGICAL                                       :: gapw,gapw_xc,gpw,pao
    LOGICAL                                       :: lrigpw
    LOGICAL                                       :: lri_optbas
    LOGICAL                                       :: ofgpw
    LOGICAL                                       :: dftb
    LOGICAL                                       :: scptb
    LOGICAL                                       :: semi_empirical
    LOGICAL                                       :: mulliken_restraint
    LOGICAL                                       :: ddapc_restraint
    LOGICAL                                       :: ddapc_restraint_is_spin
    LOGICAL                                       :: ddapc_explicit_potential
    LOGICAL                                       :: becke_restraint
    LOGICAL                                       :: et_coupling_calc
    LOGICAL                                       :: s2_restraint
    INTEGER                                       :: do_ppl_method
    INTEGER                                       :: wf_interpolation_method_nr
    INTEGER                                       :: wf_extrapolation_order
    REAL(KIND = dp)                               :: cutoff
    REAL(KIND = dp), DIMENSION(:), POINTER        :: e_cutoff
    TYPE (mulliken_restraint_type), POINTER       :: mulliken_restraint_control
    TYPE (ddapc_restraint_p_type), DIMENSION(:),POINTER :: ddapc_restraint_control
    TYPE (becke_restraint_type), POINTER          :: becke_control
    TYPE (s2_restraint_type), POINTER             :: s2_restraint_control
    TYPE (dftb_control_type),POINTER              :: dftb_control
    TYPE (scptb_control_type),POINTER             :: scptb_control
    TYPE (semi_empirical_control_type),POINTER    :: se_control
    TYPE (gapw_control_type), POINTER             :: gapw_control
    TYPE (pw_grid_option)                         :: pw_grid_opt
    LOGICAL                                       :: check_bcsr_code
    INTEGER                                       :: bcsr_code
    LOGICAL                                       :: skip_load_balance_distributed
  END TYPE qs_control_type

! *****************************************************************************
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
! *****************************************************************************
  TYPE tddfpt_control_type
     TYPE(cp_fm_p_type), DIMENSION(:), POINTER :: lumos
     REAL(KIND=dp)                             :: tolerance
     INTEGER                                   :: n_ev
     INTEGER                                   :: max_kv
     INTEGER                                   :: n_restarts
     INTEGER                                   :: n_reortho
     LOGICAL                                   :: do_kernel
     LOGICAL                                   :: lsd_singlets
     LOGICAL                                   :: invert_S
     LOGICAL                                   :: precond
     LOGICAL                                   :: drho_by_collocation
     LOGICAL                                   :: use_kinetic_energy_density
     INTEGER                                   :: res_etype
     INTEGER                                   :: diag_method
     INTEGER                                   :: oe_corr
     INTEGER                                   :: sic_method_id
     INTEGER                                   :: sic_list_id
     REAL(KIND=dp)                             :: sic_scaling_a,sic_scaling_b
     REAL(KIND = dp), DIMENSION(:,:), POINTER  :: lumos_eigenvalues
  END TYPE tddfpt_control_type

! *****************************************************************************
! \brief Control parameters for a DFT calculation
! *****************************************************************************
  TYPE dft_control_type
     TYPE(qs_control_type), POINTER    :: qs_control
     TYPE(tddfpt_control_type),POINTER :: tddfpt_control
     TYPE(xas_control_type), POINTER   :: xas_control
     TYPE(rtp_control_type),POINTER    :: rtp_control
     TYPE(admm_control_type), POINTER &
                                       :: admm_control
     TYPE(efield_p_type),POINTER,&
                        DIMENSION(:)   :: efield_fields
     TYPE(period_efield_type),POINTER  :: period_efield
     INTEGER                           :: nspins,&
                                          charge,&
                                          multiplicity,&
                                          sic_method_id,&
                                          ref_count,&
                                          id_nr,&
                                          plus_u_method_id,&
                                          dir_surf_dip
     INTEGER                           :: sic_list_id
     REAL(KIND=dp)                     :: relax_multiplicity,&
                                          sic_scaling_a,&
                                          sic_scaling_b
     LOGICAL                           :: do_tddfpt_calculation, &
                                          do_xas_calculation, &
                                          drho_by_collocation, &
                                          use_kinetic_energy_density, &
                                          restricted,&
                                          roks,&
                                          uks,&
                                          lsd,&
                                          dft_plus_u,&
                                          apply_efield_field,&
                                          apply_period_efield,&
                                          apply_external_potential,&
                                          eval_external_potential,&
                                          do_admm,&
                                          smear,&
                                          low_spin_roks,&
                                          apply_external_density,&
                                          read_external_density,&
                                          apply_external_vxc,&
                                          read_external_vxc,&
                                          correct_surf_dip

  END TYPE dft_control_type

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_control_types'
  INTEGER, SAVE :: last_dft_control_id=0

! *** Public data types ***

  PUBLIC :: dft_control_type,&
            qs_control_type,&
            gapw_control_type,&
            tddfpt_control_type,&
            efield_type,&
            mulliken_restraint_type,&
            ddapc_restraint_type,&
            dftb_control_type,&
            scptb_control_type,&
            becke_restraint_type,&
            semi_empirical_control_type,&
            s2_restraint_type,&
            admm_control_type,&
            rtp_control_type

! *** Public subroutines ***

  PUBLIC :: dft_control_retain,&
            dft_control_release,&
            dft_control_create,&
            xas_control_create,&
            tddfpt_control_create,&
            admm_control_create,&
            ddapc_control_create

CONTAINS

! *****************************************************************************
!> \brief create  the mulliken_restraint_type
!> \param mulliken_restraint_control ...
!> \param error ...
!> \par History
!>      02.2005 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE mulliken_control_create(mulliken_restraint_control,error)
    TYPE(mulliken_restraint_type), POINTER   :: mulliken_restraint_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mulliken_control_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(mulliken_restraint_control),cp_failure_level,routineP,error,failure)
    ALLOCATE(mulliken_restraint_control,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    mulliken_restraint_control%ref_count=1
    mulliken_restraint_control%strength=0.1_dp
    mulliken_restraint_control%target=1.0_dp
    mulliken_restraint_control%natoms=0
    NULLIFY(mulliken_restraint_control%atoms)
  END SUBROUTINE mulliken_control_create

! *****************************************************************************
!> \brief release the mulliken_restraint_type
!> \param mulliken_restraint_control ...
!> \param error ...
!> \par History
!>      02.2005 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE mulliken_control_release(mulliken_restraint_control,error)
    TYPE(mulliken_restraint_type), POINTER   :: mulliken_restraint_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mulliken_control_release', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(mulliken_restraint_control),cp_failure_level,routineP,error,failure)
    CPPrecondition(mulliken_restraint_control%ref_count>0,cp_failure_level,routineP,error,failure)
    mulliken_restraint_control%ref_count=mulliken_restraint_control%ref_count-1
    IF (mulliken_restraint_control%ref_count==0) THEN
       IF (ASSOCIATED(mulliken_restraint_control%atoms)) &
            DEALLOCATE(mulliken_restraint_control%atoms)
       mulliken_restraint_control%ref_count=0
       mulliken_restraint_control%strength=0.0_dp
       mulliken_restraint_control%target=0.0_dp
       mulliken_restraint_control%natoms=0
       DEALLOCATE(mulliken_restraint_control)
    ENDIF
  END SUBROUTINE mulliken_control_release

! *****************************************************************************
!> \brief retain the mulliken_restraint_type
!> \param mulliken_restraint_control ...
!> \param error ...
!> \par History
!>      02.2005 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE mulliken_control_retain(mulliken_restraint_control,error)
    TYPE(mulliken_restraint_type), POINTER   :: mulliken_restraint_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mulliken_control_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(mulliken_restraint_control),cp_failure_level,routineP,error,failure)

    mulliken_restraint_control%ref_count=mulliken_restraint_control%ref_count+1
  END SUBROUTINE mulliken_control_retain

! *****************************************************************************
!> \brief create the ddapc_restraint_type
!> \param ddapc_restraint_control ...
!> \param error ...
!> \par History
!>      02.2006 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE ddapc_control_create(ddapc_restraint_control,error)
    TYPE(ddapc_restraint_type), POINTER      :: ddapc_restraint_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ddapc_control_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(ddapc_restraint_control),cp_failure_level,routineP,error,failure)
    ALLOCATE(ddapc_restraint_control,stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    ddapc_restraint_control%density_type=do_full_density
    ddapc_restraint_control%ref_count=1
    ddapc_restraint_control%strength=0.1_dp
    ddapc_restraint_control%ddapc_order_p=0.0_dp
    ddapc_restraint_control%functional_form=-1
    ddapc_restraint_control%target=1.0_dp
    ddapc_restraint_control%natoms=0
    NULLIFY(ddapc_restraint_control%atoms)
    NULLIFY(ddapc_restraint_control%coeff)

  END SUBROUTINE ddapc_control_create

! *****************************************************************************
!> \brief release the ddapc_restraint_type
!> \param ddapc_restraint_control ...
!> \param error ...
!> \par History
!>      02.2006 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE ddapc_control_release(ddapc_restraint_control,error)
    TYPE(ddapc_restraint_type), POINTER      :: ddapc_restraint_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ddapc_control_release', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(ddapc_restraint_control),cp_failure_level,routineP,error,failure)
    CPPrecondition(ddapc_restraint_control%ref_count>0,cp_failure_level,routineP,error,failure)
    ddapc_restraint_control%ref_count=ddapc_restraint_control%ref_count-1
    IF (ddapc_restraint_control%ref_count==0) THEN
       IF (ASSOCIATED(ddapc_restraint_control%atoms)) &
            DEALLOCATE(ddapc_restraint_control%atoms)
       IF (ASSOCIATED(ddapc_restraint_control%coeff)) &
            DEALLOCATE(ddapc_restraint_control%coeff)
       ddapc_restraint_control%ref_count=0
       ddapc_restraint_control%strength=0.0_dp
       ddapc_restraint_control%target=0.0_dp
       ddapc_restraint_control%natoms=0
       DEALLOCATE(ddapc_restraint_control)
    ENDIF
  END SUBROUTINE ddapc_control_release

! *****************************************************************************
!> \brief retain the ddapc_restraint_type
!> \param ddapc_restraint_control ...
!> \param error ...
!> \par History
!>      02.2006 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE ddapc_control_retain(ddapc_restraint_control,error)
    TYPE(ddapc_restraint_type), POINTER      :: ddapc_restraint_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ddapc_control_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(ddapc_restraint_control),cp_failure_level,routineP,error,failure)

    ddapc_restraint_control%ref_count=ddapc_restraint_control%ref_count+1
  END SUBROUTINE ddapc_control_retain

! *****************************************************************************
!> \brief create the becke_restraint_type
!> \param becke_control ...
!> \param error ...
!> \par History
!>      02.2007 created [Florian Schiffmann]
! *****************************************************************************
  SUBROUTINE becke_control_create(becke_control,error)
    TYPE(becke_restraint_type), POINTER      :: becke_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'becke_control_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(becke_control),cp_failure_level,routineP,error,failure)
    ALLOCATE(becke_control,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)


    becke_control%ref_count=1
    becke_control%density_type=do_full_density
    becke_control%strength=0.1_dp
    becke_control%becke_order_p=0.0_dp
    becke_control%functional_form=-1
    becke_control%target=1.0_dp
    becke_control%natoms=0
    becke_control%need_pot=.TRUE.
    NULLIFY(becke_control%atoms)
    NULLIFY(becke_control%coeff)
  END SUBROUTINE becke_control_create

! *****************************************************************************
!> \brief release the becke_restraint_type
!> \param becke_control ...
!> \param error ...
!> \par History
!>      02.2007 created [Florian Schiffmann]
! *****************************************************************************
  SUBROUTINE becke_control_release(becke_control,error)
    TYPE(becke_restraint_type), POINTER      :: becke_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'becke_control_release', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(becke_control),cp_failure_level,routineP,error,failure)
    CPPrecondition(becke_control%ref_count>0,cp_failure_level,routineP,error,failure)
    becke_control%ref_count=becke_control%ref_count-1
    IF (becke_control%ref_count==0) THEN
       IF (ASSOCIATED(becke_control%atoms)) &
            DEALLOCATE(becke_control%atoms)
       IF (ASSOCIATED(becke_control%coeff)) &
            DEALLOCATE(becke_control%coeff)
       becke_control%ref_count=0
       becke_control%strength=0.0_dp
       becke_control%target=0.0_dp
       becke_control%natoms=0
       DEALLOCATE(becke_control)
    ENDIF
  END SUBROUTINE becke_control_release

! *****************************************************************************
!> \brief retain the becke_restraint_type
!> \param becke_control ...
!> \param error ...
!> \par History
!>      02.2007 created [Florian Schiffmann]
! *****************************************************************************
  SUBROUTINE becke_control_retain(becke_control,error)
    TYPE(becke_restraint_type), POINTER      :: becke_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'becke_control_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(becke_control),cp_failure_level,routineP,error,failure)
    becke_control%ref_count=becke_control%ref_count+1
  END SUBROUTINE becke_control_retain

! *****************************************************************************
!> \brief create the s2_restraint_type
!> \param s2_restraint_control ...
!> \param error ...
!> \par History
!>      03.2006 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE s2_control_create(s2_restraint_control,error)
    TYPE(s2_restraint_type), POINTER         :: s2_restraint_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 's2_control_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(s2_restraint_control),cp_failure_level,routineP,error,failure)
    ALLOCATE(s2_restraint_control,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    s2_restraint_control%ref_count=1
    s2_restraint_control%strength=0.1_dp
    s2_restraint_control%s2_order_p=0.0_dp
    s2_restraint_control%functional_form=-1
    s2_restraint_control%target=1.0_dp
  END SUBROUTINE s2_control_create

! *****************************************************************************
!> \brief release the s2_restraint_type
!> \param s2_restraint_control ...
!> \param error ...
!> \par History
!>      03.2006 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE s2_control_release(s2_restraint_control,error)
    TYPE(s2_restraint_type), POINTER         :: s2_restraint_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 's2_control_release', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(s2_restraint_control),cp_failure_level,routineP,error,failure)
    CPPrecondition(s2_restraint_control%ref_count>0,cp_failure_level,routineP,error,failure)
    s2_restraint_control%ref_count=s2_restraint_control%ref_count-1
    IF (s2_restraint_control%ref_count==0) THEN
       s2_restraint_control%ref_count=0
       s2_restraint_control%strength=0.0_dp
       s2_restraint_control%target=0.0_dp
       DEALLOCATE(s2_restraint_control)
    ENDIF
  END SUBROUTINE s2_control_release

! *****************************************************************************
!> \brief retain the s2_restraint_type
!> \param s2_restraint_control ...
!> \param error ...
!> \par History
!>      03.2006 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE s2_control_retain(s2_restraint_control,error)
    TYPE(s2_restraint_type), POINTER         :: s2_restraint_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 's2_control_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(s2_restraint_control),cp_failure_level,routineP,error,failure)
    s2_restraint_control%ref_count=s2_restraint_control%ref_count+1
  END SUBROUTINE s2_control_retain

! *****************************************************************************
!> \brief allocates and perform a very basic initialization
!> \param dft_control the object to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE dft_control_create(dft_control, error)
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dft_control_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(dft_control),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       ALLOCATE (dft_control,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT.failure) THEN
       dft_control%ref_count=1
       last_dft_control_id=last_dft_control_id+1
       dft_control%id_nr=last_dft_control_id
       NULLIFY(dft_control%xas_control)
       NULLIFY(dft_control%qs_control)
       NULLIFY(dft_control%tddfpt_control)
       NULLIFY(dft_control%efield_fields)
       NULLIFY(dft_control%period_efield)
       NULLIFY(dft_control%admm_control)
       NULLIFY(dft_control%rtp_control)
       CALL qs_control_create(dft_control%qs_control, error=error)
    END IF
  END SUBROUTINE dft_control_create

! *****************************************************************************
!> \brief ...
!> \param dft_control ...
!> \param error ...
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE dft_control_retain(dft_control,error)
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dft_control_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(dft_control),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPreconditionNoFail(dft_control%ref_count>0,cp_failure_level,routineP,error)
       dft_control%ref_count=dft_control%ref_count+1
    END IF
  END SUBROUTINE dft_control_retain

! *****************************************************************************
!> \brief ...
!> \param dft_control ...
!> \param error ...
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE dft_control_release(dft_control,error)
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dft_control_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    IF (ASSOCIATED(dft_control)) THEN
       CPPreconditionNoFail(dft_control%ref_count>0,cp_failure_level,routineP,error)
       dft_control%ref_count=dft_control%ref_count-1
       IF (dft_control%ref_count==0) THEN
          CALL qs_control_release(dft_control%qs_control, error=error)
          CALL tddfpt_control_release(dft_control%tddfpt_control, error=error)
          CALL xas_control_release(dft_control%xas_control, error=error)
          CALL admm_control_release(dft_control%admm_control, error=error)
          CALL efield_fields_release(dft_control%efield_fields, error=error)
          IF (ASSOCIATED(dft_control%period_efield))THEN
             DEALLOCATE(dft_control%period_efield%polarisation,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             DEALLOCATE(dft_control%period_efield,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF (ASSOCIATED(dft_control%rtp_control)) THEN
             DEALLOCATE(dft_control%rtp_control,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          ENDIF
          DEALLOCATE(dft_control, stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF
  END SUBROUTINE dft_control_release

! *****************************************************************************
!> \brief ...
!> \param gapw_control ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE gapw_control_create(gapw_control, error)
    TYPE(gapw_control_type), POINTER         :: gapw_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'gapw_control_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(.NOT.ASSOCIATED(gapw_control),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       ALLOCATE (gapw_control,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
  END SUBROUTINE gapw_control_create

! *****************************************************************************
!> \brief ...
!> \param qs_control ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE qs_control_create(qs_control, error)
    TYPE(qs_control_type), POINTER           :: qs_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_control_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(.NOT.ASSOCIATED(qs_control),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       ALLOCATE (qs_control,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       NULLIFY(qs_control%e_cutoff)
       NULLIFY(qs_control%gapw_control)
       NULLIFY(qs_control%mulliken_restraint_control)
       NULLIFY(qs_control%ddapc_restraint_control)
       NULLIFY(qs_control%s2_restraint_control)
       NULLIFY(qs_control%se_control)
       NULLIFY(qs_control%dftb_control)
       NULLIFY(qs_control%scptb_control)
       NULLIFY(qs_control%becke_control)
       NULLIFY(qs_control%ddapc_restraint_control)

       CALL mulliken_control_create(qs_control%mulliken_restraint_control,error=error)
       CALL becke_control_create(qs_control%becke_control,error=error)
       CALL s2_control_create(qs_control%s2_restraint_control,error=error)
       CALL gapw_control_create(qs_control%gapw_control, error=error)
       CALL se_control_create(qs_control%se_control,error=error)
       CALL dftb_control_create(qs_control%dftb_control,error=error)
       CALL scptb_control_create(qs_control%scptb_control,error=error)
    END IF
  END SUBROUTINE qs_control_create

! *****************************************************************************
!> \brief ...
!> \param qs_control ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE qs_control_release(qs_control, error)
    TYPE(qs_control_type), POINTER           :: qs_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_control_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(qs_control)) THEN
       CALL mulliken_control_release(qs_control%mulliken_restraint_control,error=error)
       CALL s2_control_release(qs_control%s2_restraint_control,error=error)
       CALL se_control_release(qs_control%se_control,error=error)
       CALL dftb_control_release(qs_control%dftb_control,error=error)
       CALL scptb_control_release(qs_control%scptb_control,error=error)
       CALL becke_control_release(qs_control%becke_control,error=error)

       IF (ASSOCIATED(qs_control%e_cutoff)) THEN
          DEALLOCATE(qs_control%e_cutoff,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED(qs_control%gapw_control))THEN
          DEALLOCATE(qs_control%gapw_control,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(qs_control%ddapc_restraint_control))THEN
          DO i = 1 , SIZE(qs_control%ddapc_restraint_control)
             CALL ddapc_control_release(qs_control%ddapc_restraint_control(i)%ddapc_restraint_control,&
                       error=error)
          END DO
          DEALLOCATE(qs_control%ddapc_restraint_control)
       END IF
       DEALLOCATE(qs_control,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
  END SUBROUTINE qs_control_release

! *****************************************************************************
!> \brief ...
!> \param tddfpt_control ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE tddfpt_control_create(tddfpt_control, error)
    TYPE(tddfpt_control_type), POINTER       :: tddfpt_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tddfpt_control_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(.NOT.ASSOCIATED(tddfpt_control),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       ALLOCATE (tddfpt_control,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       NULLIFY(tddfpt_control%lumos)
       NULLIFY(tddfpt_control%lumos_eigenvalues)
    END IF

  END SUBROUTINE tddfpt_control_create

! *****************************************************************************
!> \brief ...
!> \param tddfpt_control ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE tddfpt_control_release(tddfpt_control, error)
    TYPE(tddfpt_control_type), POINTER       :: tddfpt_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tddfpt_control_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin, stat
    LOGICAL                                  :: dummy, failure

    failure = .FALSE.
    IF (ASSOCIATED(tddfpt_control)) THEN
       IF (ASSOCIATED(tddfpt_control%lumos)) THEN
          DO ispin=1, SIZE(tddfpt_control%lumos)
             CALL cp_fm_release(tddfpt_control%lumos(ispin)%matrix,error=error)
             !MK the following line just avoids a crash of TDDFT runs using
             !MK the sdbg version compiled with the NAG compiler when
             !MK tddfpt_control%lumos is deallocated. This is most likely a
             !MK compiler bug and thus the line might become obsolete
             dummy = ASSOCIATED(tddfpt_control%lumos(ispin)%matrix)
          END DO
          DEALLOCATE(tddfpt_control%lumos,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED(tddfpt_control%lumos_eigenvalues)) THEN
          DEALLOCATE(tddfpt_control%lumos_eigenvalues,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       DEALLOCATE(tddfpt_control,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
  END SUBROUTINE tddfpt_control_release


! *****************************************************************************
!> \brief ...
!> \param efield_fields ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE efield_fields_release(efield_fields, error)
    TYPE(efield_p_type), DIMENSION(:), &
      POINTER                                :: efield_fields
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'efield_fields_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(efield_fields)) THEN
       DO i=1,SIZE(efield_fields)
          IF(ASSOCIATED(efield_fields(i)%efield))THEN
             IF(ASSOCIATED(efield_fields(i)%efield%envelop_r_vars))THEN
                DEALLOCATE(efield_fields(i)%efield%envelop_r_vars,stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
             IF(ASSOCIATED(efield_fields(i)%efield%envelop_i_vars))THEN
                DEALLOCATE(efield_fields(i)%efield%envelop_i_vars,stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
             IF(ASSOCIATED(efield_fields(i)%efield%polarisation))THEN
                DEALLOCATE(efield_fields(i)%efield%polarisation,stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
             DEALLOCATE(efield_fields(i)%efield,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
       END DO
       DEALLOCATE(efield_fields,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
  END SUBROUTINE efield_fields_release

! *****************************************************************************
!> \brief ...
!> \param dftb_control ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dftb_control_create(dftb_control,error)
    TYPE(dftb_control_type), POINTER         :: dftb_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dftb_control_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(dftb_control),cp_failure_level,routineP,error,failure)
    ALLOCATE(dftb_control,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    NULLIFY(dftb_control%sk_pair_list)
  END SUBROUTINE dftb_control_create

! *****************************************************************************
!> \brief ...
!> \param dftb_control ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dftb_control_release(dftb_control,error)
    TYPE(dftb_control_type), POINTER         :: dftb_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dftb_control_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(dftb_control)) THEN
       IF (ASSOCIATED(dftb_control%sk_pair_list)) THEN
          DEALLOCATE(dftb_control%sk_pair_list,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
       DEALLOCATE(dftb_control,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF
  END SUBROUTINE dftb_control_release

! *****************************************************************************
!> \brief ...
!> \param scptb_control ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE scptb_control_create(scptb_control,error)
    TYPE(scptb_control_type), POINTER        :: scptb_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scptb_control_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(scptb_control),cp_failure_level,routineP,error,failure)
    ALLOCATE(scptb_control,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE scptb_control_create

! *****************************************************************************
!> \brief ...
!> \param scptb_control ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE scptb_control_release(scptb_control,error)
    TYPE(scptb_control_type), POINTER        :: scptb_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scptb_control_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(scptb_control)) THEN
       DEALLOCATE(scptb_control,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF
  END SUBROUTINE scptb_control_release

! *****************************************************************************
!> \brief ...
!> \param se_control ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE se_control_create(se_control,error)
    TYPE(semi_empirical_control_type), &
      POINTER                                :: se_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'se_control_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(se_control),cp_failure_level,routineP,error,failure)
    ALLOCATE(se_control,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    se_control%scp = .FALSE. ! Flag is needed even in non-semiempirical methods
  END SUBROUTINE se_control_create

! *****************************************************************************
!> \brief ...
!> \param se_control ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE se_control_release(se_control,error)
    TYPE(semi_empirical_control_type), &
      POINTER                                :: se_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'se_control_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(se_control)) THEN
       DEALLOCATE(se_control,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF
  END SUBROUTINE se_control_release

! *****************************************************************************
!> \brief ...
!> \param admm_control ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE admm_control_create(admm_control,error)
    TYPE(admm_control_type), POINTER         :: admm_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'admm_control_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(admm_control),cp_failure_level,routineP,error,failure)
    ALLOCATE(admm_control,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE admm_control_create

! *****************************************************************************
!> \brief ...
!> \param admm_control ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE admm_control_release(admm_control,error)
    TYPE(admm_control_type), POINTER         :: admm_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'admm_control_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(admm_control)) THEN
       DEALLOCATE(admm_control,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF
  END SUBROUTINE admm_control_release

END MODULE cp_control_types

