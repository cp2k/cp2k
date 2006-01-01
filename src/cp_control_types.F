!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2005  CP2K developers group                          !
!-----------------------------------------------------------------------------!


!!****s* cp2k/cp_control_types [1.0] *
!!
!!   NAME
!!     cp_control_types
!!
!!   FUNCTION
!!     Defines control structures, which contain the parameters and the
!!     settings for the calculations.
!!
!!   SOURCE
!******************************************************************************

MODULE cp_control_types
  USE cp_array_r_utils,                ONLY: cp_2d_r_p_type
  USE cp_fm_types,                     ONLY: cp_fm_p_type,&
                                             cp_fm_release
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE input_constants, only: xc_rho_no_smooth, xc_rho_nn50, xc_rho_nn10, &
       xc_rho_spline2_smooth,xc_rho_spline3_smooth,xc_deriv_pw,&
       xc_deriv_spline3, xc_deriv_spline2,&
       xc_deriv_nn50_smooth, xc_deriv_nn10_smooth, xc_deriv_spline2_smooth,&
       xc_deriv_spline3_smooth,sic_none, sic_mauri_us, sic_ad, sic_mauri_spz,&
       no_excitations,&
       do_loc_min,&
       do_loc_jacobi,&
       op_loc_berry,&
       state_loc_all,&
       state_loc_fromitoj,&
       state_loc_list,&
       oe_none,&
       oe_saop,&
       oe_sic,&
       tddfpt_davidson,&
       tddfpt_excitations,&
       tddfpt_lanczos,&
       tddfpt_singlet,&
       tddfpt_triplet
  USE qs_loc_control,                  ONLY: localized_wfn_control_create, &
                                             localized_wfn_control_release,&
                                             localized_wfn_control_retain,&
                                             localized_wfn_control_type,&
                                             set_loc_wfn_lists,&
                                             set_loc_centers
  USE string_utilities,                ONLY: uppercase,&
                                             xstring
  USE memory_utilities,                ONLY: reallocate
  USE xas_control,                     ONLY: xas_control_create,&
                                             xas_control_release,&
                                             xas_control_type
  USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  TYPE xc_grid_type
     LOGICAL                   :: pw_default
!FM     TYPE ( pw_grid_type )     :: pw_grid
     INTEGER                   :: pw_mesh ( 3 )
     REAL (KIND = dp)               :: pw_cutoff
  END TYPE xc_grid_type

! *** Control parameters for pw grids ***

  TYPE pw_grid_option
     LOGICAL :: spherical
     LOGICAL :: fullspace
  END TYPE pw_grid_option

! *** Control parameters for semi empirical calculations
  TYPE semi_empirical_control_type
     LOGICAL               :: orthogonal_basis
     LOGICAL               :: ewald
     REAL(KIND = dp)       :: delta
     REAL(KIND = dp)       :: alpha_gauss
     REAL(KIND = dp)       :: eps_gauss
     REAL(KIND = dp)       :: rc_exchange
     REAL(KIND = dp)       :: rc_coulomb
  END TYPE semi_empirical_control_type

! *** Control parameters for GAPW method within QUICKSTEP ***

  TYPE gapw_control_type
    REAL(KIND = dp)                 :: eps_fit,eps_iso,eps_Vrho0,eps_svd
    INTEGER                         :: ladd_rho0,lmax_rho0,&
                                       lmax_sphere, quadrature
    CHARACTER(LEN=30)               :: quadrature_name
    LOGICAL                         :: lrho1_eq_lrho0
    LOGICAL                         :: alpha0_hard_from_input, &
                                       alpha0_soft_from_input
    REAL(KIND = dp)                 :: alpha0_hard, alpha0_soft
    REAL(KIND = dp)                 :: max_rad_local
  END TYPE gapw_control_type

! *** some parameters useful for mulliken_restraints
  TYPE mulliken_restraint_type
    INTEGER                         :: ref_count
    REAL(KIND = dp)                 :: strength
    REAL(KIND = dp)                 :: TARGET
    INTEGER                         :: natoms
    INTEGER, POINTER, DIMENSION(:)  :: atoms
  END TYPE mulliken_restraint_type


! *** Control parameters for a KIM-GORDON calculation ***
! *** Control parameters for a QUICKSTEP calculation ***

! eps_pgf_orb: Cutoff value for the interaction of the primitive Gaussian-type
!              functions (primitive basis functions).

  TYPE qs_control_type
!MK PRIVATE would be nice
    CHARACTER(LEN=10)               :: method
    INTEGER                         :: method_id
    REAL(KIND = dp)                 :: eps_core_charge,&
                                       eps_pgf_orb,&
                                       eps_ppl,&
                                       eps_ppnl,&
                                       eps_rho_gspace,&
                                       eps_rho_rspace,&
                                       eps_gvg_rspace,&
                                       progression_factor,&
                                       relative_cutoff
    LOGICAL                         :: commensurate_mgrids
    LOGICAL                         :: realspace_mgrids
    LOGICAL                         :: map_paa
    LOGICAL                         :: map_consistent
    LOGICAL                         :: gapw,gapw_xc,gpw,pao
    LOGICAL                         :: semi_empirical
    LOGICAL                         :: polarization
    LOGICAL                         :: mulliken_restraint
    LOGICAL                         :: localized_wfn
    TYPE (mulliken_restraint_type), POINTER :: mulliken_restraint_control
    TYPE(localized_wfn_control_type), POINTER :: localized_wfn_control
    TYPE (semi_empirical_control_type)      :: se_control
    REAL(KIND = dp)                         :: cutoff
    REAL(KIND = dp), DIMENSION(:), POINTER  :: e_cutoff
    TYPE (gapw_control_type), POINTER       :: gapw_control
    TYPE ( pw_grid_option )         :: pw_grid_opt
    INTEGER                         :: pw_grid_opt_id
    INTEGER                         :: rs_distributed_id
    INTEGER                         :: wf_interpolation_method_nr
    INTEGER                         :: wf_extrapolation_order
  END TYPE qs_control_type

! *** Control parameters for a WAVE calculation ***

  TYPE wave_control_type
     REAL(KIND = dp) :: e_cutoff_wf
     REAL(KIND = dp) :: e_cutoff_dual
  END TYPE wave_control_type

  ! *** Control parameters for a TIME-DEPENDENT PERTURBATION calculation ***

!!****s* cp_control_types/tddfpt_control_type
!!
!! NAME
!!   tddfpt_control_type
!!
!! FUNCTION
!!   A type that holds controling information for a time-dependent
!!   perturbation calculation.
!!
!! ATTRIBUTES
!!   - n_ev       : number of eigenvalues to calculate
!!   - n_reortho  : how many time to reorthogonalize (in the lanczos algorithm)
!!   - do_kernel  : wether to evaluate the kernel (this is a debugging option)
!!   - res_etype : { SINGLET | TRIPLET } which excitations
!!                  to calculate
!!   - lumos_eigenvalues : holds the eigenvalues of the lumos (if calculated in QS)
!!
!! NOTES  
!!   The lumos are helpfull in choosing a initial vector for the TDDFPT 
!!   calculation, since they can be used to construct the solutions of the
!!   TDDFPT operator without the perturbation kernel.
!!  
!!***
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
     LOGICAL                                   :: use_kinetic_energy_density
     INTEGER                                   :: res_etype
     INTEGER                                   :: diag_method
     INTEGER                                   :: oe_corr
     INTEGER                                   :: sic_method_id
     REAL(KIND=dp)                             :: sic_scaling_a,sic_scaling_b
     REAL(KIND = dp), DIMENSION(:,:), POINTER  :: lumos_eigenvalues
  END TYPE tddfpt_control_type

  ! *** Control parameters for a DFT calculation ***

  TYPE dft_control_type
     TYPE(qs_control_type)           :: qs_control
     TYPE(tddfpt_control_type)       :: tddfpt_control
     TYPE(wave_control_type)         :: wave_control
     TYPE(xas_control_type), POINTER :: xas_control
     TYPE(xc_grid_type)              :: integration_grid
     INTEGER                         :: nspins,&
                                        charge,&
                                        multiplicity,&
                                        sic_method_id,&
                                        ref_count,&
                                        id_nr,&
                                        xas_estate
     REAL(KIND=dp)                   :: sic_scaling_a,sic_scaling_b
     LOGICAL                         :: forces,&
                                        do_tddfpt_calculation, &
                                        do_xas_calculation, &
                                        derive_function_rho, &
                                        use_kinetic_energy_density, & 
                                        restricted
     !FM add xc_section ?
  END TYPE dft_control_type

  CHARACTER(LEN=*), PARAMETER :: moduleN = "cp_control_types"

  INTEGER, SAVE :: last_dft_control_id=0

! *** Public data types ***

  PUBLIC :: dft_control_type, xc_grid_type,&
            qs_control_type,&
            gapw_control_type,&
            wave_control_type, tddfpt_control_type,&
            pw_grid_option, mulliken_restraint_type,&
            semi_empirical_control_type


! *** Public subroutines ***

  PUBLIC :: dft_control_retain,&
            dft_control_release,&
            dft_control_create,&
            gapw_control_create,&
            set_loc_wfn_lists,&
            set_loc_centers,&
            mulliken_control_create,&
            mulliken_control_retain,&
            mulliken_control_release

!!***
! *****************************************************************************

CONTAINS

!!****f* cp_control_types/mulliken_control_create/retain/release *
!!
!!   NAME
!!      mulliken_control_create/retain/release
!!   FUNCTION
!!      create retain releasethe mulliken_restraint_type
!!   NOTES
!!      
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     02.2005 created [Joost VandeVondele]
!!
!!   SOURCE
!!*** **********************************************************************
SUBROUTINE mulliken_control_create(mulliken_restraint_control,error)
   TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
   TYPE(mulliken_restraint_type), POINTER :: mulliken_restraint_control
   CHARACTER(len=*), PARAMETER :: routineN = 'mulliken_control_create', &
      routineP = moduleN//':'//routineN

   INTEGER :: stat
   LOGICAL :: failure 
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

SUBROUTINE mulliken_control_release(mulliken_restraint_control,error)
   TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
   TYPE(mulliken_restraint_type), POINTER :: mulliken_restraint_control
   CHARACTER(len=*), PARAMETER :: routineN = 'mulliken_control_release', &
      routineP = moduleN//':'//routineN
   LOGICAL :: failure 
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

SUBROUTINE mulliken_control_retain(mulliken_restraint_control,error)
   TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
   TYPE(mulliken_restraint_type), POINTER :: mulliken_restraint_control
   CHARACTER(len=*), PARAMETER :: routineN = 'mulliken_control_retain', &
      routineP = moduleN//':'//routineN
   LOGICAL :: failure 
   failure=.FALSE.
   CPPrecondition(ASSOCIATED(mulliken_restraint_control),cp_failure_level,routineP,error,failure)  

   mulliken_restraint_control%ref_count=mulliken_restraint_control%ref_count+1
END SUBROUTINE mulliken_control_retain

 
! *****************************************************************************

!!****f* cp_control_types/dft_control_create [1.0] *
!!
!!   NAME
!!     dft_control_create
!!
!!   FUNCTION
!!     allocates and perform a very basic initialization
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - dft_control: the object to create
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE dft_control_create(dft_control, error)
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
     CALL qs_control_init(dft_control%qs_control, error=error)
     CALL tddfpt_control_init(dft_control%tddfpt_control, error=error)
     CALL wave_control_init(dft_control%wave_control, error=error)
     CALL xas_control_create(dft_control%xas_control,error=error)
     CALL xc_grid_init(dft_control%integration_grid, error=error)
  END IF
END SUBROUTINE dft_control_create
!***************************************************************************

!!****f* cp_control_types/dft_control_retain [1.0] *
!!
!!   NAME
!!     dft_control_retain
!!
!!   FUNCTION
!!     retains the given dft_control
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - dft_control: the dft_control to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE dft_control_retain(dft_control,error)
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
!***************************************************************************

!!****f* cp_control_types/dft_control_release [1.0] *
!!
!!   NAME
!!     dft_control_release
!!
!!   FUNCTION
!!     releases the given dft_control
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - dft_control: the dft_control to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE dft_control_release(dft_control,error)
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dft_control_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  IF (ASSOCIATED(dft_control)) THEN
     CPPreconditionNoFail(dft_control%ref_count>0,cp_failure_level,routineP,error)
     dft_control%ref_count=dft_control%ref_count-1
     IF (dft_control%ref_count==0) THEN
        CALL qs_control_dealloc_ref(dft_control%qs_control, error=error)
        CALL tddfpt_control_dealloc_ref(dft_control%tddfpt_control, error=error)
        CALL wave_control_dealloc_ref(dft_control%wave_control, error=error)
        CALL xas_control_release(dft_control%xas_control, error=error)
        CALL xc_grid_dealloc_ref(dft_control%integration_grid, error=error)
        DEALLOCATE(dft_control, stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(dft_control)
END SUBROUTINE dft_control_release
!***************************************************************************

SUBROUTINE gapw_control_create(gapw_control, error)
    TYPE(gapw_control_type), POINTER         :: gapw_control
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
!***************************************************************************

SUBROUTINE qs_control_init(qs_control, error)
    TYPE(qs_control_type), INTENT(out)       :: qs_control
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

  NULLIFY(qs_control%e_cutoff)
  NULLIFY(qs_control%gapw_control)
  NULLIFY(qs_control%mulliken_restraint_control)
  NULLIFY(qs_control%localized_wfn_control)

  CALL mulliken_control_create(qs_control%mulliken_restraint_control)

  CALL localized_wfn_control_create(qs_control%localized_wfn_control, error=error)

END SUBROUTINE qs_control_init

SUBROUTINE qs_control_dealloc_ref(qs_control, error)
    TYPE(qs_control_type), INTENT(inout)     :: qs_control
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = "qs_control_dealloc_ref", &
      routineP = moduleN//"/"//routineN

    INTEGER                                  :: stat

    CALL mulliken_control_release(qs_control%mulliken_restraint_control)
    CALL localized_wfn_control_release(qs_control% localized_wfn_control, error=error)

    IF (ASSOCIATED(qs_control%e_cutoff)) THEN
       DEALLOCATE(qs_control%e_cutoff,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF
    IF (ASSOCIATED(qs_control%gapw_control))THEN
       DEALLOCATE(qs_control%gapw_control,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF

END SUBROUTINE qs_control_dealloc_ref

SUBROUTINE tddfpt_control_init(tddfpt_control, error)
    TYPE(tddfpt_control_type), INTENT(out)   :: tddfpt_control
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    NULLIFY(tddfpt_control%lumos,tddfpt_control%lumos_eigenvalues)
END SUBROUTINE tddfpt_control_init

SUBROUTINE tddfpt_control_dealloc_ref(tddfpt_control, error)
    TYPE(tddfpt_control_type), INTENT(inout) :: tddfpt_control
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = "tddfpt_control_dealloc_ref", &
      routineP = moduleN//"/"//routineN

    INTEGER                                  :: ispin, stat

  IF (ASSOCIATED(tddfpt_control%lumos)) THEN
     DO ispin=1, SIZE(tddfpt_control%lumos)
        CALL cp_fm_release(tddfpt_control%lumos(ispin)%matrix)
     END DO
     DEALLOCATE(tddfpt_control%lumos,stat=stat)
     CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
  END IF
  IF (ASSOCIATED(tddfpt_control%lumos_eigenvalues)) THEN
     DEALLOCATE(tddfpt_control%lumos_eigenvalues,stat=stat)
     CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
  END IF
END SUBROUTINE tddfpt_control_dealloc_ref

SUBROUTINE wave_control_init(wave_control, error)
    TYPE(wave_control_type), INTENT(out)     :: wave_control
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

END SUBROUTINE wave_control_init

SUBROUTINE wave_control_dealloc_ref(wave_control, error)
    TYPE(wave_control_type), INTENT(inout)   :: wave_control
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

END SUBROUTINE wave_control_dealloc_ref

SUBROUTINE xc_grid_init(xc_grid, error)

    TYPE(xc_grid_type), INTENT(out)          :: xc_grid
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

END SUBROUTINE xc_grid_init

SUBROUTINE xc_grid_dealloc_ref(xc_grid, error)
    TYPE(xc_grid_type), INTENT(inout)        :: xc_grid
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

END SUBROUTINE xc_grid_dealloc_ref


END MODULE cp_control_types
