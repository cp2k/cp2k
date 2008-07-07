!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Utilities to set up the control types
! *****************************************************************************
MODULE cp_control_utils
  USE bibliography,                    ONLY: &
       Dewar1977, Dewar1985, Elstner1998, Krack2000, Lippert1997, &
       Lippert1999, Porezag1995, Repasky2002, Rocha2006, Seifert1996, &
       Stewart1989, Stewart2007, Thiel1992, VandeVondele2005a, &
       VandeVondele2005b, Zhechkov2005, cite_reference
  USE cp_control_types,                ONLY: &
       dft_control_create, dft_control_type, finite_field_control_create, &
       qs_control_type, scp_control_create, tddfpt_control_create, &
       tddfpt_control_type, wfn_fitting_control_create, xas_control_create
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE f77_blas
  USE input_constants
  USE input_cp2k_check,                ONLY: xc_functionals_expand
  USE input_cp2k_dft,                  ONLY: create_dft_section,&
                                             create_qs_section
  USE input_enumeration_types,         ONLY: enum_i2c,&
                                             enumeration_type
  USE input_keyword_types,             ONLY: keyword_get,&
                                             keyword_type
  USE input_section_types,             ONLY: &
       section_get_keyword, section_release, section_type, section_vals_get, &
       section_vals_get_subs_vals, section_vals_get_subs_vals2, &
       section_vals_type, section_vals_val_get, section_vals_val_set
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE semi_empirical_int_arrays,       ONLY: init_se_intd_array
  USE termination,                     ONLY: stop_program
  USE xas_control,                     ONLY: read_xas_control
  USE xc_derivatives,                  ONLY: xc_functional_get_info,&
                                             xc_functionals_get_needs
  USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_control_utils'

! *** Public subroutines ***

  PUBLIC :: read_dft_control,&
            read_mgrid_section,&
            read_qs_section,&
            read_tddfpt_control,&
            write_dft_control,&
            write_qs_control,&
            read_ddapc_section,&
            read_becke_section

CONTAINS

! *****************************************************************************
  SUBROUTINE read_dft_control(dft_control,dft_section, error)
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(section_vals_type), POINTER         :: dft_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_dft_control', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_path_length)       :: basis_set_file_name, &
                                                potential_file_name
    INTEGER                                  :: excitations, method_id
    LOGICAL                                  :: failure, is_present, l_param, &
                                                not_SE, was_present
    REAL(KIND=dp)                            :: density_cut, gradient_cut, &
                                                tau_cut
    REAL(KIND=dp), DIMENSION(:), POINTER     :: rvec
    TYPE(section_vals_type), POINTER         :: scf_section, scp_section, &
                                                tmp_section, xc_fun_section, &
                                                xc_section
    TYPE(xc_rho_cflags_type)                 :: needs

    failure=.FALSE.
    was_present = .FALSE.
    NULLIFY(tmp_section,xc_fun_section,xc_section)
    CALL dft_control_create(dft_control, error=error)
    ! determine wheather this is a semiempirical or DFTB run 
    ! --> (no XC section needs to be provided)
    not_SE = .TRUE.
    CALL section_vals_val_get(dft_section,"QS%METHOD",i_val=method_id,error=error)
    SELECT CASE (method_id)
    CASE (do_method_dftb,do_method_mndo,do_method_am1,do_method_pm3,&
          do_method_pm6,do_method_pdg,do_method_rm1,do_method_mndod)
       not_SE = .FALSE.
    END SELECT
    ! Check for XC section and XC_FUNCTIONAL section
    xc_section => section_vals_get_subs_vals(dft_section,"XC",error=error)
    CALL section_vals_get(xc_section,explicit=is_present,error=error)
    IF (.NOT.is_present.AND.not_SE) THEN
       CALL stop_program(routineN,"XC section missing.")
    END IF
    IF (is_present) THEN
       CALL section_vals_val_get(xc_section,"density_cutoff", r_val=density_cut,error=error)
       CALL section_vals_val_get(xc_section,"gradient_cutoff", r_val=gradient_cut,error=error)
       CALL section_vals_val_get(xc_section,"tau_cutoff", r_val=tau_cut,error=error)
       ! Perform numerical stability checks and possibly correct the issues
       CALL cp_assert(density_cut>EPSILON(0.0_dp)*100.0_dp,cp_warning_level,cp_assertion_failed,routineP,&
            "DENSITY_CUTOFF lower than 100*EPSILON, where EPSILON is the machine precision. "//&
            "This may lead to numerical problems. Setting up shake_tol to 100*EPSILON! "//&
CPSourceFileRef,&
            only_ionode=.TRUE.)
       density_cut = MAX(EPSILON(0.0_dp)*100.0_dp,density_cut)
       CALL cp_assert(gradient_cut>EPSILON(0.0_dp)*100.0_dp,cp_warning_level,cp_assertion_failed,routineP,&
            "GRADIENT_CUTOFF lower than 100*EPSILON, where EPSILON is the machine precision. "//&
            "This may lead to numerical problems. Setting up shake_tol to 100*EPSILON! "//&
CPSourceFileRef,&
            only_ionode=.TRUE.)
       gradient_cut = MAX(EPSILON(0.0_dp)*100.0_dp,gradient_cut)
       CALL cp_assert(tau_cut>EPSILON(0.0_dp)*100.0_dp,cp_warning_level,cp_assertion_failed,routineP,&
            "TAU_CUTOFF lower than 100*EPSILON, where EPSILON is the machine precision. "//&
            "This may lead to numerical problems. Setting up shake_tol to 100*EPSILON! "//&
CPSourceFileRef,&
            only_ionode=.TRUE.)
       tau_cut = MAX(EPSILON(0.0_dp)*100.0_dp,tau_cut)
       CALL section_vals_val_set(xc_section,"density_cutoff", r_val=density_cut,error=error)
       CALL section_vals_val_set(xc_section,"gradient_cutoff", r_val=gradient_cut,error=error)
       CALL section_vals_val_set(xc_section,"tau_cutoff", r_val=tau_cut,error=error)
    END IF
    xc_fun_section => section_vals_get_subs_vals(xc_section,"XC_FUNCTIONAL",error=error)
    CALL section_vals_get(xc_fun_section,explicit=is_present,error=error)
    IF (.NOT.is_present.AND.not_SE) THEN
       CALL stop_program(routineN,"XC_FUNCTIONAL section missing.")
    END IF
    scf_section => section_vals_get_subs_vals(dft_section,"SCF",error=error)
    CALL section_vals_val_get(dft_section,"UKS",l_val=dft_control%uks,error=error)
    CALL section_vals_val_get(dft_section,"ROKS",l_val=dft_control%roks,error=error)
    IF (dft_control%uks.OR.dft_control%roks) THEN
       dft_control%nspins = 2
    ELSE
       dft_control%nspins = 1
    END IF

    dft_control%lsd = (dft_control%nspins > 1)
    needs = xc_functionals_get_needs(xc_fun_section,&
                                    lsd=dft_control%lsd,&
                                    add_basic_components=.TRUE.,&
                                    error=error)
    dft_control%use_kinetic_energy_density = (needs%tau_spin.OR.needs%tau)

    !! check if we do wavefunction fitting
    tmp_section => section_vals_get_subs_vals(dft_section,"WFN_FITTING",error=error)
    CALL section_vals_get(tmp_section,explicit=is_present,error=error)
    dft_control%do_wfn_fitting = is_present
    IF(is_present) THEN
      CALL section_vals_val_get(dft_section,"WFN_FITTING%METHOD",i_val=method_id,error=error)
      CALL wfn_fitting_control_create(dft_control%wfn_fitting_control, error)
      dft_control%wfn_fitting_control%method_id = method_id
    END IF

    CALL section_vals_val_get(dft_section,"SCP%_SECTION_PARAMETERS_",l_val=dft_control%scp,error=error)
    IF (dft_control%scp) THEN
       CALL scp_control_create(dft_control%scp_control, error=error)
       scp_section => section_vals_get_subs_vals(dft_section,"SCP", error=error)
       CALL section_vals_val_get(scp_section,"DISPERSION",l_val=l_param,error=error)
       dft_control%scp_control%dispersion = l_param
    END IF

    ! Set restricted to true, if both OT and ROKS are requested
    !MK in principle dft_control%restricted could be dropped completely like the
    !MK input key by using only dft_control%roks now
    CALL section_vals_val_get(scf_section,"OT%_SECTION_PARAMETERS_",l_val=l_param,error=error)
    dft_control%restricted = (dft_control%roks.AND.l_param)

    CALL section_vals_val_get(dft_section,"CHARGE",i_val=dft_control%charge,error=error)
    CALL section_vals_val_get(dft_section,"MULTIPLICITY",i_val=dft_control%multiplicity,error=error)

    dft_control%sic_method_id = sic_none
    dft_control%sic_scaling_a = 1.0_dp
    dft_control%sic_scaling_b = 1.0_dp

    dft_control%xas_estate    = -1

    ! DFT+U
    dft_control%dft_plus_u    = .FALSE.
    CALL section_vals_val_get(dft_section,"PLUS_U_METHOD",i_val=method_id,error=error)
    dft_control%plus_u_method_id = method_id

    CALL section_vals_val_get(dft_section,"BASIS_SET_FILE_NAME",&
         c_val=basis_set_file_name , error=error)
    CALL section_vals_val_get(dft_section,"POTENTIAL_FILE_NAME",&
         c_val=potential_file_name , error=error)

!   *** Read the input section ***
    tmp_section => section_vals_get_subs_vals(dft_section,"sic",error=error)
    CALL section_vals_val_get(tmp_section,"SIC_METHOD",&
         i_val=dft_control%sic_method_id,error=error)
    CALL section_vals_val_get(tmp_section,"ORBITAL_SET",&
         i_val=dft_control%sic_list_id,error=error)
    CALL section_vals_val_get(tmp_section,"SIC_SCALING_A",&
         r_val=dft_control%sic_scaling_a,error=error)
    CALL section_vals_val_get(tmp_section,"SIC_SCALING_B",&
         r_val=dft_control%sic_scaling_b,error=error)

!   *** Determine if this is a  TDDFPT  run *** 
    CALL section_vals_val_get(dft_section,"EXCITATIONS",i_val=excitations,error=error)
    dft_control%do_tddfpt_calculation = (excitations==tddfpt_excitations)
    IF (dft_control%do_tddfpt_calculation) THEN
       CALL tddfpt_control_create(dft_control%tddfpt_control, error=error)
    END IF

!   *** Read the input section ***
    tmp_section => section_vals_get_subs_vals(dft_section,"XAS",error=error)
    CALL section_vals_get(tmp_section,explicit=dft_control%do_xas_calculation,error=error)
    IF (dft_control%do_xas_calculation) THEN
       ! Override with section parameter
       CALL section_vals_val_get(tmp_section,"_SECTION_PARAMETERS_",&
            l_val=dft_control%do_xas_calculation,error=error)
       IF (dft_control%do_xas_calculation) THEN
          CALL xas_control_create(dft_control%xas_control,error=error)       
          CALL read_xas_control(dft_control%xas_control,dft_section,error=error)
       END IF
    END IF

!   *** Read the finite field input section ***
    tmp_section => section_vals_get_subs_vals(dft_section,"FIELD",error=error)
    CALL section_vals_get(tmp_section,explicit=dft_control%apply_finite_field,error=error)
    IF (dft_control%apply_finite_field) THEN
       CALL finite_field_control_create(dft_control%finite_field_control, error=error)
       CALL section_vals_val_get(tmp_section,"ELECTRICAL_FIELD_STRENGTH",&
            r_val=dft_control%finite_field_control%efield_strength,error=error)
       CALL section_vals_val_get(tmp_section,"ELECTRICAL_FIELD_POLARISATION",&
            r_vals=rvec,error=error)
       dft_control%finite_field_control%efield_polarisation(1:3)=rvec(1:3)
       CALL section_vals_val_get(tmp_section,"BERRY_PHASE",&
            l_val=dft_control%finite_field_control%berry_phase,error=error)
    END IF
  END SUBROUTINE read_dft_control

! *****************************************************************************
  SUBROUTINE read_mgrid_section(qs_control,dft_section,error)

    TYPE(qs_control_type), INTENT(INOUT)     :: qs_control
    TYPE(section_vals_type), POINTER         :: dft_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_mgrid_section', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: igrid_level, ngrid_level, stat
    LOGICAL                                  :: failure, multigrid_set
    REAL(dp)                                 :: cutoff
    REAL(dp), DIMENSION(:), POINTER          :: cutofflist
    TYPE(section_vals_type), POINTER         :: mgrid_section

    failure=.FALSE.
    NULLIFY(mgrid_section, cutofflist)
    mgrid_section => section_vals_get_subs_vals(dft_section,"MGRID",error=error)

    CALL section_vals_val_get(mgrid_section,"NGRIDS",i_val=ngrid_level,error=error)
    CALL section_vals_val_get(mgrid_section,"MULTIGRID_SET",l_val=multigrid_set,error=error)
    IF ( qs_control%semi_empirical .OR. qs_control%dftb ) THEN
       ngrid_level = 1
       multigrid_set = .FALSE.
    END IF

    ALLOCATE(qs_control%e_cutoff(ngrid_level),STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    CALL section_vals_val_get(mgrid_section,"CUTOFF",r_val=cutoff,error=error)
    qs_control%cutoff = MAX(0.01_dp,cutoff)
    CALL section_vals_val_get(mgrid_section,"PROGRESSION_FACTOR",&
                 r_val=qs_control%progression_factor,error=error)
    CALL section_vals_val_get(mgrid_section,"COMMENSURATE",&
                 l_val=qs_control%commensurate_mgrids,error=error)
    CALL section_vals_val_get(mgrid_section,"REALSPACE",&
                 l_val=qs_control%realspace_mgrids,error=error)
    CALL section_vals_val_get(mgrid_section,"REL_CUTOFF",&
                 r_val=qs_control%relative_cutoff,error=error)

    IF ( qs_control%semi_empirical .OR. qs_control%dftb ) THEN
      qs_control%cutoff = 1._dp
    END IF

    IF(multigrid_set) THEN
      ! Read the values from input
      IF (qs_control%commensurate_mgrids) THEN
          CALL stop_program(routineN,&
               "Do not specify cutoffs for the commensurate grids (NYI)")
      ENDIF

      CALL section_vals_val_get(mgrid_section,"MULTIGRID_CUTOFF",&
           r_vals=cutofflist,error=error)
      IF(ASSOCIATED(cutofflist)) THEN
         CALL cp_assert(SIZE(cutofflist,1) == ngrid_level,&
              cp_failure_level,cp_assertion_failed,routineN,&
              "Inconsistent values for number of multi grids",error,failure)
         IF (.NOT.failure) THEN
            DO igrid_level = 1,ngrid_level
               qs_control%e_cutoff(igrid_level) = cutofflist(igrid_level)*0.5_dp
            ENDDO
         END IF
      END IF
      ! set cutoff to smallest value in multgrid available with >= cutoff
      DO igrid_level=ngrid_level,1,-1
        IF ( qs_control%cutoff <= qs_control%e_cutoff(igrid_level) ) THEN
          qs_control%cutoff = qs_control%e_cutoff(igrid_level)
          EXIT
        END IF
        ! set largest grid value to cutoff
        IF ( igrid_level == 1) THEN
          qs_control%cutoff = qs_control%e_cutoff(1)
        END IF
      END DO
    ELSE
      IF (qs_control%commensurate_mgrids) qs_control%progression_factor = 4.0_dp
      qs_control%e_cutoff(1) = qs_control%cutoff
      DO igrid_level=2,ngrid_level
        qs_control%e_cutoff(igrid_level) = qs_control%e_cutoff(igrid_level-1) &
            / qs_control%progression_factor
      END DO
    END IF
    ! check that multigrids are ordered
    DO igrid_level=2,ngrid_level
      IF ( qs_control%e_cutoff(igrid_level) > &
           qs_control%e_cutoff(igrid_level-1) ) THEN
         CALL stop_program (routineN,"Multi grids not ordered")
      END IF
    END DO
  END SUBROUTINE read_mgrid_section

! *****************************************************************************
  SUBROUTINE read_qs_section(qs_control,qs_section,error)

    TYPE(qs_control_type), INTENT(INOUT)     :: qs_control
    TYPE(section_vals_type), POINTER         :: qs_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_qs_section', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length), &
      DIMENSION(:), POINTER                  :: clist
    INTEGER                                  :: itmp, j, jj, k, n_rep, n_var, &
                                                stat
    INTEGER, DIMENSION(:), POINTER           :: tmplist
    LOGICAL                                  :: failure, was_present
    REAL(dp)                                 :: tmpsqrt, value
    TYPE(enumeration_type), POINTER          :: enum
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: section
    TYPE(section_vals_type), POINTER :: becke_constraint_section, &
      ddapc_restraint_section, dftb_parameter, dftb_section, mull_section, &
      s2_restraint_section, se_section

    failure=.FALSE.
    was_present = .FALSE.
    NULLIFY(mull_section,se_section,dftb_section)

    mull_section => section_vals_get_subs_vals(qs_section,"MULLIKEN_RESTRAINT", error=error)
    ddapc_restraint_section => section_vals_get_subs_vals(qs_section,"DDAPC_RESTRAINT", error=error)
    s2_restraint_section => section_vals_get_subs_vals(qs_section,"S2_RESTRAINT", error=error)
    becke_constraint_section=> section_vals_get_subs_vals(qs_section,"BECKE_CONSTRAINT", error=error)
    se_section => section_vals_get_subs_vals(qs_section,"SE", error=error)
    dftb_section => section_vals_get_subs_vals(qs_section,"DFTB", error=error)
    dftb_parameter => section_vals_get_subs_vals(dftb_section,"PARAMETER", error=error)

    ! Setup all defaults values and overwrite input parameters
    CALL section_vals_val_get(qs_section,"EPS_DEFAULT",r_val=value,error=error)
    tmpsqrt=SQRT(value) ! a trick to work around a NAG 5.1 optimizer bug
    qs_control%eps_core_charge = value/100.0_dp
    qs_control%eps_pgf_orb = tmpsqrt
    qs_control%eps_ppl = value/100.0_dp
    qs_control%eps_ppnl = value/100.0_dp
    qs_control%eps_rho_gspace = value
    qs_control%eps_rho_rspace = value
    qs_control%eps_gvg_rspace = tmpsqrt
    qs_control%eps_kg_gpw = tmpsqrt*1000._dp
    !
    CALL section_vals_val_get(qs_section,"EPS_CORE_CHARGE",n_rep_val=n_rep,error=error)
    IF (n_rep /=0) THEN
       CALL section_vals_val_get(qs_section,"EPS_CORE_CHARGE",r_val=qs_control%eps_core_charge,error=error)
    END IF
    CALL section_vals_val_get(qs_section,"EPS_GVG_RSPACE",n_rep_val=n_rep,error=error)
    IF (n_rep /=0) THEN
       CALL section_vals_val_get(qs_section,"EPS_GVG_RSPACE",r_val=qs_control%eps_gvg_rspace,error=error)
    END IF
    CALL section_vals_val_get(qs_section,"EPS_PGF_ORB",n_rep_val=n_rep,error=error)
    IF (n_rep /=0) THEN
       CALL section_vals_val_get(qs_section,"EPS_PGF_ORB",r_val=qs_control%eps_pgf_orb,error=error)
    END IF
    CALL section_vals_val_get(qs_section,"EPS_PPL",n_rep_val=n_rep,error=error)
    IF (n_rep /=0) THEN
       CALL section_vals_val_get(qs_section,"EPS_PPL",r_val=qs_control%eps_ppl,error=error)
    END IF
    CALL section_vals_val_get(qs_section,"EPS_PPNL",n_rep_val=n_rep,error=error)
    IF (n_rep /=0) THEN
       CALL section_vals_val_get(qs_section,"EPS_PPNL",r_val=qs_control%eps_ppnl,error=error)
    END IF
    CALL section_vals_val_get(qs_section,"EPS_RHO",n_rep_val=n_rep,error=error)
    IF (n_rep /=0) THEN
       CALL section_vals_val_get(qs_section,"EPS_RHO",r_val=qs_control%eps_rho_gspace,error=error)
       qs_control%eps_rho_rspace = qs_control%eps_rho_gspace
    END IF
    CALL section_vals_val_get(qs_section,"EPS_RHO_RSPACE",n_rep_val=n_rep,error=error)
    IF (n_rep /=0) THEN
       CALL section_vals_val_get(qs_section,"EPS_RHO_RSPACE",r_val=qs_control%eps_rho_rspace,error=error)
    END IF
    CALL section_vals_val_get(qs_section,"EPS_RHO_GSPACE",n_rep_val=n_rep,error=error)
    IF (n_rep /=0) THEN
       CALL section_vals_val_get(qs_section,"EPS_RHO_GSPACE",r_val=qs_control%eps_rho_gspace,error=error)
    END IF
    CALL section_vals_val_get(qs_section,"EPS_KG_GPW",n_rep_val=n_rep,error=error)
    IF (n_rep /=0) THEN
       CALL section_vals_val_get(qs_section,"EPS_KG_GPW",r_val=qs_control%eps_kg_gpw,error=error)
    END IF

    CALL section_vals_val_get(qs_section,"EPS_3C_REDUCE",r_val=qs_control%gapw_control%eps_3c_reduce,error=error)
    CALL section_vals_val_get(qs_section,"EPSFIT",r_val=qs_control%gapw_control%eps_fit,error=error)
    CALL section_vals_val_get(qs_section,"EPSISO",r_val=qs_control%gapw_control%eps_iso,error=error)
    CALL section_vals_val_get(qs_section,"EPSSVD",r_val=qs_control%gapw_control%eps_svd,error=error)
    CALL section_vals_val_get(qs_section,"EPSRHO0",r_val=qs_control%gapw_control%eps_Vrho0,error=error)
    CALL section_vals_val_get(qs_section,"ALPHA0_HARD",r_val=qs_control%gapw_control%alpha0_hard,error=error)
    qs_control%gapw_control%lrho1_eq_lrho0 = .FALSE.
    qs_control%gapw_control%alpha0_hard_from_input = .FALSE.
    IF(qs_control%gapw_control%alpha0_hard/=0.0_dp) qs_control%gapw_control%alpha0_hard_from_input = .TRUE.
    CALL section_vals_val_get(qs_section,"ALPHA0_SOFT",r_val=qs_control%gapw_control%alpha0_soft,error=error)
    qs_control%gapw_control%alpha0_soft_from_input = .FALSE.
    IF(qs_control%gapw_control%alpha0_soft/=0.0_dp) qs_control%gapw_control%alpha0_soft_from_input = .TRUE.
    CALL section_vals_val_get(qs_section,"FULL_GAPW",l_val=qs_control%gapw_control%full_gapw,error=error)
    CALL section_vals_val_get(qs_section,"MAX_RAD_LOCAL",r_val=qs_control%gapw_control%max_rad_local,error=error)

    ! Logicals
    CALL section_vals_val_get(qs_section,"MAP_PAA",l_val=qs_control%map_paa,error=error)
    CALL section_vals_val_get(qs_section,"MAP_CONSISTENT",l_val=qs_control%map_consistent,error=error)

    ! Integers gapw
    CALL section_vals_val_get(qs_section,"LMAXN1",i_val=qs_control%gapw_control%lmax_sphere,error=error)
    CALL section_vals_val_get(qs_section,"LMAXN0",i_val=qs_control%gapw_control%lmax_rho0,error=error)
    CALL section_vals_val_get(qs_section,"LADDN0",i_val=qs_control%gapw_control%ladd_rho0,error=error)
    CALL section_vals_val_get(qs_section,"QUADRATURE",i_val=qs_control%gapw_control%quadrature,error=error)
    qs_control%gapw_control%quadrature_name=quadrature_name(qs_control%gapw_control%quadrature)

    ! Integers grids
    CALL section_vals_val_get(qs_section,"PW_GRID",i_val=itmp,error=error)
    SELECT CASE (itmp)
    CASE(do_pwgrid_spherical)
       qs_control%pw_grid_opt%spherical = .TRUE.
       qs_control%pw_grid_opt%fullspace = .FALSE.
    CASE (do_pwgrid_ns_fullspace)
       qs_control%pw_grid_opt%spherical = .FALSE.
       qs_control%pw_grid_opt%fullspace = .TRUE.
    CASE (do_pwgrid_ns_halfspace)
       qs_control%pw_grid_opt%spherical = .FALSE.
       qs_control%pw_grid_opt%fullspace = .FALSE.
    END SELECT

    CALL section_vals_val_get(qs_section,"PW_GRID_LAYOUT",i_vals=tmplist,error=error)
    qs_control%pw_grid_opt%distribution_layout=tmplist

    !Integers extrapolation
    CALL section_vals_val_get(qs_section,"EXTRAPOLATION",i_val=qs_control%wf_interpolation_method_nr,error=error)
    CALL section_vals_val_get(qs_section,"EXTRAPOLATION_ORDER",i_val=qs_control%wf_extrapolation_order,error=error)

    !Method
    CALL section_vals_val_get(qs_section,"METHOD",i_val=qs_control%method_id,error=error)
    NULLIFY(section,keyword,enum)
    CALL create_qs_section(section,error=error)
    keyword => section_get_keyword(section,"METHOD",error=error)
    CALL keyword_get(keyword,enum=enum,error=error)
    qs_control%method = enum_i2c(enum,qs_control%method_id,error=error)
    CALL section_release(section,error=error)
    qs_control%gapw           = .FALSE.
    qs_control%gapw_xc        = .FALSE.
    qs_control%gpw            = .FALSE.
    qs_control%pao            = .FALSE.
    qs_control%polarization   = .FALSE.
    qs_control%dftb           = .FALSE.
    qs_control%semi_empirical = .FALSE.
    SELECT CASE (qs_control%method_id)
    CASE (do_method_gapw)
       CALL cite_reference(Lippert1999)
       CALL cite_reference(Krack2000)
       qs_control%gapw = .TRUE.
    CASE (do_method_gapw_xc)
       qs_control%gapw_xc = .TRUE.
    CASE (do_method_gpw)
       CALL cite_reference(Lippert1997)
       CALL cite_reference(VandeVondele2005a)
       qs_control%gpw = .TRUE.
    CASE (do_method_dftb)
       qs_control%dftb = .TRUE.
       CALL cite_reference(Porezag1995)
       CALL cite_reference(Seifert1996)
    CASE (do_method_mndo)
       CALL cite_reference(Dewar1977)
       qs_control%semi_empirical = .TRUE.
    CASE (do_method_am1)
       CALL cite_reference(Dewar1985)
       qs_control%semi_empirical = .TRUE.
    CASE (do_method_pm3)
       CALL cite_reference(Stewart1989)
       qs_control%semi_empirical = .TRUE.
    CASE (do_method_pm6)
       CALL cite_reference(Stewart2007)
       qs_control%semi_empirical = .TRUE.
    CASE (do_method_pdg)
       CALL cite_reference(Repasky2002)
       qs_control%semi_empirical = .TRUE.
    CASE (do_method_rm1)
       CALL cite_reference(Rocha2006)
       qs_control%semi_empirical = .TRUE.
    CASE (do_method_mndod)
       CALL cite_reference(Dewar1977)
       CALL cite_reference(Thiel1992)
       qs_control%semi_empirical = .TRUE.
    END SELECT
    
    IF(qs_control%method_id == do_method_kg_pol) qs_control%polarization = .TRUE.

    CALL section_vals_get(mull_section,explicit=qs_control%mulliken_restraint,error=error)

    IF (qs_control%mulliken_restraint) THEN
       CALL section_vals_val_get(mull_section,"STRENGTH",r_val=qs_control%mulliken_restraint_control%strength,&
            error=error)
       CALL section_vals_val_get(mull_section,"TARGET",r_val=qs_control%mulliken_restraint_control%target,error=error)
       CALL section_vals_val_get(mull_section,"ATOMS",n_rep_val=n_rep,error=error)
       jj = 0
       DO k = 1,n_rep
          CALL section_vals_val_get(mull_section,"ATOMS",i_rep_val=k,i_vals=tmplist, error=error)
          jj=jj+SIZE(tmplist)
       END DO
       qs_control%mulliken_restraint_control%natoms = jj
       IF (qs_control%mulliken_restraint_control%natoms<1) &
            CALL stop_program(routineN, "Need at least 1 atom to use mulliken contraints")
       ALLOCATE(qs_control%mulliken_restraint_control%atoms(qs_control%mulliken_restraint_control%natoms))
       jj = 0
       DO k = 1,n_rep
          CALL section_vals_val_get(mull_section,"ATOMS",i_rep_val=k,i_vals=tmplist, error=error)
          DO  j = 1,SIZE(tmplist)
             jj = jj+1
             qs_control%mulliken_restraint_control%atoms (jj) = tmplist(j)
          END DO
       END DO
    ENDIF

    CALL section_vals_get(ddapc_restraint_section,explicit=qs_control%ddapc_restraint,error=error)

    IF (qs_control%ddapc_restraint) THEN
       CALL  read_ddapc_section(qs_control,ddapc_restraint_section,error)
    ENDIF

    CALL section_vals_get(s2_restraint_section,explicit=qs_control%s2_restraint,error=error)
    IF (qs_control%s2_restraint) THEN
       CALL section_vals_val_get(s2_restraint_section,"STRENGTH", &
            r_val=qs_control%s2_restraint_control%strength,error=error)
       CALL section_vals_val_get(s2_restraint_section,"TARGET", &
            r_val=qs_control%s2_restraint_control%target,error=error)
       CALL section_vals_val_get(s2_restraint_section,"FUNCTIONAL_FORM", &
            i_val=qs_control%s2_restraint_control%functional_form,error=error)
    ENDIF

    CALL section_vals_get(becke_constraint_section,explicit=qs_control%becke_constraint,error=error)

    IF (qs_control%becke_constraint) THEN
       CALL  read_becke_section(qs_control,becke_constraint_section,error)
    ENDIF

    ! Semi-empirical code
    IF (qs_control%semi_empirical) THEN
       CALL init_se_intd_array(error)
       CALL section_vals_val_get(se_section,"ANALYTICAL_GRADIENTS",&
            l_val=qs_control%se_control%analytical_gradients,error=error)
       CALL section_vals_val_get(se_section,"ORTHOGONAL_BASIS",&
            l_val=qs_control%se_control%orthogonal_basis,error=error)
       CALL section_vals_val_get(se_section,"DO_EWALD",&
            l_val=qs_control%se_control%do_ewald,error=error)
       CALL section_vals_val_get(se_section,"DELTA",&
            r_val=qs_control%se_control%delta,error=error)
       CALL section_vals_val_get(se_section,"RC_INTERACTION",&
            r_val=qs_control%se_control%rc_interaction,error=error)
       CALL section_vals_val_get(se_section,"RC_COULOMB",&
            r_val=qs_control%se_control%rc_coulomb,error=error)
       CALL section_vals_val_get(se_section,"RC_RANGE",&
            r_val=qs_control%se_control%rc_range,error=error)
       IF(  qs_control %method_id == do_method_mndo .OR. &
            qs_control %method_id == do_method_am1  .OR. &
            qs_control %method_id == do_method_mndod.OR. &
            qs_control %method_id == do_method_pdg  .OR. &
            qs_control %method_id == do_method_pm3  .OR. &
            qs_control %method_id == do_method_pm6  .OR. &
            qs_control %method_id == do_method_rm1  ) THEN
          qs_control%se_control%orthogonal_basis=.TRUE.
       END IF
    END IF

    ! DFTB code
    IF (qs_control%dftb) THEN
       CALL section_vals_val_get(dftb_section,"ORTHOGONAL_BASIS",&
            l_val=qs_control%dftb_control%orthogonal_basis,error=error)
       CALL section_vals_val_get(dftb_section,"SELF_CONSISTENT",&
            l_val=qs_control%dftb_control%self_consistent,error=error)
       CALL section_vals_val_get(dftb_section,"DISPERSION",&
            l_val=qs_control%dftb_control%dispersion,error=error)
       CALL section_vals_val_get(dftb_section,"EPS_DISP",&
            r_val=qs_control%dftb_control%eps_disp,error=error)
       CALL section_vals_val_get(dftb_section,"DO_EWALD",&
            l_val=qs_control%dftb_control%do_ewald,error=error)
       CALL section_vals_val_get(dftb_parameter,"PARAM_FILE_PATH",&
            c_val=qs_control%dftb_control%sk_file_path,error=error)
       CALL section_vals_val_get(dftb_parameter,"PARAM_FILE_NAME",&
            c_val=qs_control%dftb_control%sk_file_list,error=error)
       CALL section_vals_val_get(dftb_parameter,"SK_FILE",n_rep_val=n_var,&
            error=error)
       ALLOCATE(qs_control%dftb_control%sk_pair_list(3,n_var),STAT=stat)
       CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
       DO k = 1, n_var
          CALL section_vals_val_get(dftb_parameter,"SK_FILE",i_rep_val=k,&
               c_vals=clist,error=error)
          qs_control%dftb_control%sk_pair_list(1:3,k) = clist(1:3)
       END DO
       CALL section_vals_val_get(dftb_parameter,"UFF_FORCE_FIELD",&
            c_val=qs_control%dftb_control%uff_force_field,error=error)
       
       IF (qs_control%dftb) THEN
          IF (qs_control%dftb_control%dispersion)      CALL cite_reference(Zhechkov2005)
          IF (qs_control%dftb_control%self_consistent) CALL cite_reference(Elstner1998)
       ENDIF
    END IF
  END SUBROUTINE read_qs_section

! *****************************************************************************
  SUBROUTINE read_tddfpt_control(t_control,nspins,dft_section,error)
    TYPE(tddfpt_control_type)                :: t_control
    INTEGER, INTENT(IN)                      :: nspins
    TYPE(section_vals_type), POINTER         :: dft_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_tddfpt_control', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure, kenergy_den
    TYPE(section_vals_type), POINTER         :: sic_section, t_section

    failure=.FALSE.
    kenergy_den = .FALSE.
    NULLIFY(sic_section, t_section)
    t_section => section_vals_get_subs_vals(dft_section,"TDDFPT",error=error)

    CALL section_vals_val_get(t_section,"CONVERGENCE", r_val=t_control%tolerance,error=error)
    CALL section_vals_val_get(t_section,"NEV",i_val=t_control%n_ev,error=error)
    CALL section_vals_val_get(t_section,"MAX_KV",i_val=t_control%max_kv,error=error)
    CALL section_vals_val_get(t_section,"RESTARTS",i_val=t_control%n_restarts,error=error)
    CALL section_vals_val_get(t_section,"NREORTHO",i_val=t_control%n_reortho,error=error)
    CALL section_vals_val_get(t_section,"RES_ETYPE",i_val=t_control%res_etype,error=error)
    CALL section_vals_val_get(t_section,"DIAG_METHOD",i_val=t_control%diag_method,error=error)
    CALL section_vals_val_get(t_section,"KERNEL",l_val=t_control%do_kernel,error=error)
    CALL section_vals_val_get(t_section,"LSD_SINGLETS",l_val=t_control%lsd_singlets,error=error)
    CALL section_vals_val_get(t_section,"INVERT_S",l_val=t_control%invert_S,error=error)
    CALL section_vals_val_get(t_section,"PRECOND",l_val=t_control%precond,error=error)
    CALL section_vals_val_get(t_section,"OE_CORR",i_val=t_control%oe_corr,error=error)

    t_control%use_kinetic_energy_density = .FALSE.
    sic_section=>section_vals_get_subs_vals(t_section,"SIC",error=error)
    CALL section_vals_val_get(sic_section,"SIC_METHOD",i_val=t_control%sic_method_id,error=error)
    CALL section_vals_val_get(sic_section,"ORBITAL_SET",i_val=t_control%sic_list_id,error=error)
    CALL section_vals_val_get(sic_section,"SIC_SCALING_A",r_val=t_control%sic_scaling_a,error=error)
    CALL section_vals_val_get(sic_section,"SIC_SCALING_B",r_val=t_control%sic_scaling_b,error=error)

  END SUBROUTINE read_tddfpt_control

! *****************************************************************************
!> \brief Write the DFT control parameters to the output unit. 
! *****************************************************************************
  SUBROUTINE write_dft_control(dft_control, dft_section, error)
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(section_vals_type), POINTER         :: dft_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'write_dft_control', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=160)                       :: reference
    CHARACTER(len=20)                        :: tmpStr
    INTEGER                                  :: ifun, il, myfun, output_unit
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: density_cut, &
                                                density_smooth_cut_range, &
                                                gradient_cut, tau_cut
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(enumeration_type), POINTER          :: enum
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: section
    TYPE(section_vals_type), POINTER         :: xc_fun, xc_fun_section, &
                                                xc_section

    failure = .FALSE.
    IF (dft_control%qs_control%semi_empirical) RETURN
    IF (dft_control%qs_control%dftb) RETURN
    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    output_unit = cp_print_key_unit_nr(logger,dft_section,&
         "PRINT%DFT_CONTROL_PARAMETERS",extension=".Log",error=error)
    IF (output_unit>0) THEN
       xc_section => section_vals_get_subs_vals(dft_section,&
            "XC",error=error)
       IF (dft_control%uks) THEN
          WRITE (UNIT=output_unit,FMT="(/,T2,A,T78,A)")&
               "DFT| Spin unrestricted (spin-polarized) Kohn-Sham calculation","UKS"
       ELSE IF (dft_control%roks) THEN
          WRITE (UNIT=output_unit,FMT="(/,T2,A,T77,A)")&
               "DFT| Spin restricted open Kohn-Sham calculation","ROKS"
       ELSE
          WRITE (UNIT=output_unit,FMT="(/,T2,A,T78,A)")&
               "DFT| Spin restricted Kohn-Sham (RKS) calculation","RKS"
       END IF

       WRITE (UNIT=output_unit,FMT="(T2,A,T76,I5)")&
            "DFT| Multiplicity",dft_control%multiplicity
       WRITE (UNIT=output_unit,FMT="(T2,A,T76,I5)")&
            "DFT| Number of spin states",dft_control%nspins

       IF (dft_control%apply_finite_field) THEN
          IF (dft_control%finite_field_control%berry_phase) THEN
             WRITE ( output_unit, '( A, T51, A )' ) " DFT| Finite electrical field ",&
                  "Polarisation type: Berry phase"
          ELSE
             WRITE ( output_unit, '( A )' ) " DFT| Finite electrical field "
          END IF
          WRITE ( output_unit, '( A, T61, F20.5 )' ) " DFT| Field strength [a.u.]",&
               dft_control%finite_field_control%efield_strength
          WRITE ( output_unit, '( A, T51, 3F10.5 )' ) " DFT| Field polarisation",&
               dft_control%finite_field_control%efield_polarisation
          CALL stop_program(routineN,moduleN,__LINE__,&
               "Finite field not operational")
       END IF

       WRITE (UNIT=output_unit,FMT="(T2,A,T76,I5)")&
         "DFT| Charge",dft_control%charge

       IF (dft_control%sic_method_id.NE.sic_none) CALL cite_reference(VandeVondele2005b)
       SELECT CASE ( dft_control%sic_method_id )
       CASE ( sic_none )
         tmpstr = "NO"
       CASE ( sic_mauri_spz )
         tmpstr = "SPZ/MAURI SIC"
       CASE ( sic_mauri_us )
         tmpstr = "US/MAURI SIC"
       CASE ( sic_ad )
         tmpstr = "AD SIC"
       CASE ( sic_eo )
         tmpstr = "Explicit Orbital SIC"
       CASE DEFAULT
         ! fix throughout the cp2k for this option
         CALL stop_program(routineN,moduleN,__LINE__,&
                           "SIC option unknown")
       END SELECT

       WRITE (UNIT=output_unit,FMT="(T2,A,T61,A20)")&
         "DFT| Self-interaction correction (SIC)",ADJUSTR(TRIM(tmpstr))

       IF (dft_control%sic_method_id /= sic_none) THEN
         WRITE (UNIT=output_unit,FMT="(T2,A,T66,ES15.6)")&
           "DFT| SIC scaling parameter a",dft_control%sic_scaling_a,&
           "DFT| SIC scaling parameter b",dft_control%sic_scaling_b
       END IF

       IF (dft_control%sic_method_id == sic_eo) THEN
         IF (dft_control%sic_list_id==sic_list_all) THEN
           WRITE (UNIT=output_unit,FMT="(T2,A,T66,A)")&
             "DFT| SIC orbitals","ALL"
         ENDIF
         IF (dft_control%sic_list_id==sic_list_unpaired) THEN
           WRITE (UNIT=output_unit,FMT="(T2,A,T66,A)")&
             "DFT| SIC orbitals","UNPAIRED"
         ENDIF
       END IF

       CALL section_vals_val_get(xc_section,"density_cutoff", r_val=density_cut,error=error)
       CALL section_vals_val_get(xc_section,"gradient_cutoff", r_val=gradient_cut,error=error)
       CALL section_vals_val_get(xc_section,"tau_cutoff", r_val=tau_cut,error=error)
       CALL section_vals_val_get(xc_section,"density_smooth_cutoff_range",r_val=density_smooth_cut_range,error=error)
       
       WRITE (UNIT=output_unit,FMT="(T2,A,T66,ES15.6)")&
            "DFT| Cutoffs: density ",density_cut,&
            "DFT|          gradient",gradient_cut,&
            "DFT|          tau     ",tau_cut,&
            "DFT|          cutoff_smoothing_range",density_smooth_cut_range
       CALL section_vals_val_get(xc_section,"XC_GRID%XC_SMOOTH_RHO",&
            c_val=tmpStr,error=error)
       WRITE ( output_unit, '( A, T61, A )' ) &
            " DFT| XC density smoothing ",ADJUSTR(tmpStr)
       CALL section_vals_val_get(xc_section,"XC_GRID%XC_DERIV",&
            c_val=tmpStr,error=error)
       WRITE ( output_unit, '( A, T61, A )' ) &
            " DFT| XC derivatives ",ADJUSTR(tmpStr)
       IF (dft_control%dft_plus_u) THEN
         NULLIFY (enum,keyword,section)
         CALL create_dft_section(section,error=error)
         keyword => section_get_keyword(section,"PLUS_U_METHOD",error=error)
         CALL keyword_get(keyword,enum=enum,error=error)
         WRITE (UNIT=output_unit,FMT="(T2,A,T41,A40)")&
           "DFT+U| Method",ADJUSTR(TRIM(enum_i2c(enum,dft_control%plus_u_method_id,error=error)))
         WRITE (UNIT=output_unit,FMT="(T2,A)")&
           "DFT+U| Check atomic kind information for details"
         CALL section_release(section,error=error)
       END IF

       xc_fun_section => section_vals_get_subs_vals(xc_section,&
            "XC_FUNCTIONAL",error=error)
       CALL section_vals_val_get(xc_fun_section,"_SECTION_PARAMETERS_",i_val=myfun,error=error)
       IF (myfun/=xc_none) THEN
          CALL section_vals_val_get(xc_section,"FUNCTIONAL_ROUTINE",&
               c_val=tmpStr,error=error)
          WRITE (output_unit, fmt="(' FUNCTIONAL| ROUTINE=',a)")&
               TRIM(tmpStr)
          CALL xc_functionals_expand(xc_fun_section,xc_section,error=error)
          ifun=0
          DO
             ifun=ifun+1
             xc_fun => section_vals_get_subs_vals2(xc_fun_section,i_section=ifun,&
                  error=error)
             IF (.NOT.ASSOCIATED(xc_fun)) EXIT
             CALL xc_functional_get_info(xc_fun,&
                                         lsd=dft_control%lsd,&
                                         reference=reference,&
                                         error=error)
             WRITE (output_unit, fmt="(' FUNCTIONAL| ',a,':')")&
                  TRIM(xc_fun%section%name)
             il = LEN_TRIM(reference)
             IF (il < 67) THEN
               WRITE (output_unit, fmt="(' FUNCTIONAL| ',a)") reference(1:il)
             ELSE IF (il < 134) THEN
               WRITE (output_unit, fmt="(' FUNCTIONAL| ',a)") reference(1:67)
               WRITE (output_unit, fmt="(' FUNCTIONAL| ',a)") reference(68:il)
             ELSE
               WRITE (output_unit, fmt="(' FUNCTIONAL| ',a)") reference(1:67)
               WRITE (output_unit, fmt="(' FUNCTIONAL| ',a)") reference(68:134)
               WRITE (output_unit, fmt="(' FUNCTIONAL| ',a)") reference(134:il)
             END IF
          END DO
       ELSE
          WRITE (output_unit, fmt="(' FUNCTIONAL| NO EXCHAGE-CORRELATION FUNCTIONAL USED.')")
       END IF
    END IF
    CALL cp_print_key_finished_output(output_unit,logger,dft_section,&
         "PRINT%DFT_CONTROL_PARAMETERS",error=error)

  END SUBROUTINE write_dft_control

! *****************************************************************************
!> \brief Purpose: Write the QS control parameters to the output unit.
! *****************************************************************************
  SUBROUTINE write_qs_control(qs_control,dft_section,error)
    TYPE(qs_control_type), INTENT(IN)        :: qs_control
    TYPE(section_vals_type), POINTER         :: dft_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'write_qs_control', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: igrid_level, ngrid_level, &
                                                output_unit
    TYPE(cp_logger_type), POINTER            :: logger

    IF (qs_control%semi_empirical) RETURN
    IF (qs_control%dftb) RETURN
    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    output_unit = cp_print_key_unit_nr(logger, dft_section,&
         "PRINT%DFT_CONTROL_PARAMETERS",extension=".Log",error=error)
    IF (output_unit>0) THEN
       ngrid_level = SIZE(qs_control%e_cutoff)
       WRITE (UNIT=output_unit,FMT="(/,T2,A,T71,A)")&
            "QS| Method:",ADJUSTR(qs_control%method)
       IF ( qs_control% pw_grid_opt%spherical ) THEN
          WRITE (UNIT=output_unit,FMT="(T2,A,T61,A)")&
               "QS| Density plane wave grid type"," SPHERICAL HALFSPACE"
       ELSE IF ( qs_control% pw_grid_opt%fullspace ) THEN
          WRITE (UNIT=output_unit,FMT="(T2,A,T57,A)")&
               "QS| Density plane wave grid type"," NON-SPHERICAL FULLSPACE"
       ELSE
          WRITE (UNIT=output_unit,FMT="(T2,A,T57,A)")&
               "QS| Density plane wave grid type"," NON-SPHERICAL HALFSPACE"
       END IF
       WRITE (UNIT=output_unit,FMT="(T2,A,T71,I10)")&
            "QS| Number of grid levels:",SIZE(qs_control%e_cutoff)
       IF (ngrid_level == 1) THEN
          WRITE (UNIT=output_unit,FMT="(T2,A,T71,F10.1)")&
               "QS| Density cutoff [a.u.]:",qs_control%e_cutoff(1)
       ELSE
          WRITE (UNIT=output_unit,FMT="(T2,A,T71,F10.1)")&
               "QS| Density cutoff [a.u.]:",qs_control%cutoff
          IF (qs_control%commensurate_mgrids) &
               WRITE (UNIT=output_unit,FMT="(T2,A)") "QS| Using commensurate multigrids"
          WRITE (UNIT=output_unit,FMT="(T2,A,T71,F10.1)")&
               "QS| Multi grid cutoff [a.u.]: 1) grid level",qs_control%e_cutoff(1)
          WRITE (UNIT=output_unit,FMT="(T2,A,I3,A,T71,F10.1)")&
               ("QS|                         ",igrid_level,") grid level",&
               qs_control%e_cutoff(igrid_level),&
               igrid_level=2,SIZE(qs_control%e_cutoff))
       END IF
       IF (qs_control%pao) THEN
          WRITE (UNIT=output_unit,FMT="(T2,A)") "QS| PAO active"
       END IF
       WRITE (UNIT=output_unit,FMT="(T2,A,T71,F10.1)")&
            "QS| Grid level progression factor:",qs_control%progression_factor
       WRITE (UNIT=output_unit,FMT="(T2,A,T71,F10.1)")&
            "QS| Relative density cutoff [a.u.]:",qs_control%relative_cutoff
       IF (qs_control%map_consistent) THEN
          WRITE (UNIT=output_unit,FMT="(T2,A)")&
               "QS| Consistent realspace mapping and integration "
       ENDIF
       WRITE (UNIT=output_unit,FMT="(T2,A,T73,ES8.1)")&
            "QS| Interaction thresholds: eps_pgf_orb:",&
            qs_control%eps_pgf_orb,&
            "QS|                         eps_core_charge:",&
            qs_control%eps_core_charge,&
            "QS|                         eps_rho_gspace:",&
            qs_control%eps_rho_gspace,&
            "QS|                         eps_rho_rspace:",&
            qs_control%eps_rho_rspace,&
            "QS|                         eps_gvg_rspace:",&
            qs_control%eps_gvg_rspace,&
            "QS|                         eps_ppl:",&
            qs_control%eps_ppl,&
            "QS|                         eps_ppnl:",&
            qs_control%eps_ppnl
       IF (qs_control%gapw) THEN
          WRITE (UNIT=output_unit,FMT="(T2,A,T73,ES8.1)")&
               "QS| GAPW|                   eps_fit:",&
               qs_control%gapw_control%eps_fit,&
               "QS| GAPW|                   eps_iso:",&
               qs_control%gapw_control%eps_iso,&
               "QS| GAPW|                   eps_svd:",&
               qs_control%gapw_control%eps_svd
          WRITE (UNIT=output_unit,FMT="(T2,A,T55,A30)")&
               "QS| GAPW|   atom-r-grid: quadrature:",&
               qs_control%gapw_control%quadrature_name
          WRITE (UNIT=output_unit,FMT="(T2,A,T71,I10)")&
               "QS| GAPW|      atom-s-grid:  max l :",&
               qs_control%gapw_control%lmax_sphere ,&
               "QS| GAPW|      max_l_rho0 :",&
               qs_control%gapw_control%lmax_rho0
       END IF
       IF (qs_control%gapw_xc) THEN
          WRITE (UNIT=output_unit,FMT="(T2,A,T73,ES8.1)")&
               "QS| GAPW_XC|                eps_fit:",&
               qs_control%gapw_control%eps_fit,&
               "QS| GAPW_XC|                eps_iso:",&
               qs_control%gapw_control%eps_iso,&
               "QS| GAPW_XC|                eps_svd:",&
               qs_control%gapw_control%eps_svd
          WRITE (UNIT=output_unit,FMT="(T2,A,T55,A30)")&
               "QS| GAPW_XC|atom-r-grid: quadrature:",&
               qs_control%gapw_control%quadrature_name
          WRITE (UNIT=output_unit,FMT="(T2,A,T71,I10)")&
               "QS| GAPW_XC|   atom-s-grid:  max l :",&
               qs_control%gapw_control%lmax_sphere
       END IF
       IF (qs_control%mulliken_restraint) THEN
          WRITE (UNIT=output_unit,FMT="(T2,A,T73,ES8.1)")&
               "QS| Mulliken restraint target", qs_control%mulliken_restraint_control%target
          WRITE (UNIT=output_unit,FMT="(T2,A,T73,ES8.1)")&
               "QS| Mulliken restraint strength", qs_control%mulliken_restraint_control%strength
          WRITE (UNIT=output_unit,FMT="(T2,A,T73,I8)")&
               "QS| Mulliken restraint atoms: ", qs_control%mulliken_restraint_control%natoms
          WRITE (UNIT=output_unit,FMT="(5I8)") qs_control%mulliken_restraint_control%atoms
       END IF
       IF (qs_control%ddapc_restraint) THEN
          WRITE (UNIT=output_unit,FMT="(T2,A,T73,ES8.1)")&
               "QS| ddapc restraint target", qs_control%ddapc_restraint_control%target
          WRITE (UNIT=output_unit,FMT="(T2,A,T73,ES8.1)")&
               "QS| ddapc restraint strength", qs_control%ddapc_restraint_control%strength
          WRITE (UNIT=output_unit,FMT="(T2,A,T73,I8)")&
               "QS| ddapc restraint atoms: ", qs_control%ddapc_restraint_control%natoms
          WRITE (UNIT=output_unit,FMT="(5I8)") qs_control%ddapc_restraint_control%atoms
          WRITE (UNIT=output_unit,FMT="(T2,A)") "Coefficients:"
          WRITE (UNIT=output_unit,FMT="(5F6.2)") qs_control%ddapc_restraint_control%coeff
          SELECT CASE(qs_control%ddapc_restraint_control%functional_form)
          CASE (do_ddapc_restraint)
             WRITE (UNIT=output_unit,FMT="(T2,A,T61,A20)")&
                  "QS| ddapc restraint functional form :","RESTRAINT"
          CASE (do_ddapc_constraint)
             WRITE (UNIT=output_unit,FMT="(T2,A,T61,A20)")&
                  "QS| ddapc restraint functional form :","CONSTRAINT"
          CASE DEFAULT
             CALL stop_program(routineN,"unknown ddapc restraint")
          END SELECT
       END IF
       IF (qs_control%s2_restraint) THEN
          WRITE (UNIT=output_unit,FMT="(T2,A,T73,ES8.1)")&
               "QS| s2 restraint target", qs_control%s2_restraint_control%target
          WRITE (UNIT=output_unit,FMT="(T2,A,T73,ES8.1)")&
               "QS| s2 restraint strength", qs_control%s2_restraint_control%strength
          SELECT CASE(qs_control%s2_restraint_control%functional_form)
          CASE (do_s2_restraint)
             WRITE (UNIT=output_unit,FMT="(T2,A,T61,A20)")&
                  "QS| s2 restraint functional form :","RESTRAINT"
             CALL stop_program(routineN,"Not yet implemented")
          CASE (do_s2_constraint)
             WRITE (UNIT=output_unit,FMT="(T2,A,T61,A20)")&
                  "QS| s2 restraint functional form :","CONSTRAINT"
          CASE DEFAULT
             CALL stop_program(routineN,"unknown ddapc restraint")
          END SELECT
       END IF
    END IF
    CALL cp_print_key_finished_output(output_unit,logger,dft_section,&
         "PRINT%DFT_CONTROL_PARAMETERS",error=error)
  END SUBROUTINE write_qs_control

! *****************************************************************************
!> \brief returns a string that describes the smoothing of rho
!> \param xc_rho_smooth_id the id that represent the smoothing
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      05.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  FUNCTION xc_get_rho_smooth_label(xc_rho_smooth_id,error) RESULT(res)
    INTEGER, INTENT(in)                      :: xc_rho_smooth_id
    TYPE(cp_error_type), INTENT(inout)       :: error
    CHARACTER(len=10)                        :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_get_rho_smooth_label', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    SELECT CASE(xc_rho_smooth_id)
    CASE (xc_rho_no_smooth)
       res="NONE"
    CASE (xc_rho_spline2_smooth)
       res="SPLINE2"
    CASE (xc_rho_spline3_smooth)
       res="SPLINE3"
    CASE (xc_rho_nn10)
       res="NN10"
    CASE (xc_rho_nn50)
       res="NN50"
    CASE default
       WRITE (res,"('UKN',i6)") xc_rho_smooth_id
    END SELECT
  END FUNCTION xc_get_rho_smooth_label

! *****************************************************************************
!> \brief returns a string that describes the derivative used in the xc
!>      calculation
!> \param xc_deriv_method_id the id that represent the derivative method
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      06.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  FUNCTION xc_get_deriv_method_label(xc_deriv_method_id,error) RESULT(res)
    INTEGER, INTENT(in)                      :: xc_deriv_method_id
    TYPE(cp_error_type), INTENT(inout)       :: error
    CHARACTER(len=20)                        :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_get_deriv_method_label', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    SELECT CASE(xc_deriv_method_id)
    CASE (xc_deriv_pw)
       res="PW"
    CASE (xc_deriv_spline2)
       res="SPLINE2"
    CASE (xc_deriv_spline3)
       res="SPLINE3"
    CASE (xc_deriv_spline2_smooth)
       res="SPLINE2_SMOOTH"
    CASE (xc_deriv_spline3_smooth)
       res="SPLINE3_SMOOTH"
    CASE (xc_deriv_nn10_smooth)
       res="NN10_SMOOTH"
    CASE (xc_deriv_nn50_smooth)
       res="NN50_SMOOTH"
    CASE default
       WRITE (res,"('UKN',i6)") xc_deriv_method_id
    END SELECT
  END FUNCTION xc_get_deriv_method_label

! *****************************************************************************
!> \brief returns a string that describes the functional routine used in the xc
!>      calculation
!> \param xc_functional_routine_id the id that represent the functional
!>        routine
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Fawzi Mohamed
! *****************************************************************************
  FUNCTION xc_get_routine_label(xc_functional_routine_id,error) RESULT(res)
    INTEGER, INTENT(in)                      :: xc_functional_routine_id
    TYPE(cp_error_type), INTENT(inout)       :: error
    CHARACTER(len=20)                        :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_get_routine_label', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    SELECT CASE(xc_functional_routine_id)
    CASE(xc_old_f_routine)
       res="OLD"
    CASE(xc_new_f_routine)
       res="NEW"
    CASE(xc_test_lsd_f_routine)
       res="TEST_LSD"
    CASE(xc_debug_new_routine)
       res="DEBUG_NEW"
    CASE default
       WRITE (res,"('UKN',i6)") xc_functional_routine_id
    END SELECT
  END FUNCTION xc_get_routine_label

! *****************************************************************************
!> \brief reads the input parameters needed for ddapc.
!> \note
!>      either reads DFT%QS%DDAPC_RESTRAINT or PROPERTIES%ET_coupling
!> \author fschiff
! *****************************************************************************
  SUBROUTINE read_ddapc_section(qs_control,ddapc_restraint_section,error)

    TYPE(qs_control_type), INTENT(INOUT)     :: qs_control
    TYPE(section_vals_type), POINTER         :: ddapc_restraint_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_ddapc_section', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: j, jj, k, n_rep
    INTEGER, DIMENSION(:), POINTER           :: tmplist
    REAL(KIND=dp), DIMENSION(:), POINTER     :: rtmplist

    CALL section_vals_val_get(ddapc_restraint_section,"STRENGTH", &
         r_val=qs_control%ddapc_restraint_control%strength,error=error)
    CALL section_vals_val_get(ddapc_restraint_section,"TARGET", &
         r_val=qs_control%ddapc_restraint_control%target,error=error)
    CALL section_vals_val_get(ddapc_restraint_section,"FUNCTIONAL_FORM", &
         i_val=qs_control%ddapc_restraint_control%functional_form,error=error)
    CALL section_vals_val_get(ddapc_restraint_section,"ATOMS", &
         n_rep_val=n_rep,error=error)
    jj = 0
    DO k = 1,n_rep
       CALL section_vals_val_get(ddapc_restraint_section,"ATOMS",i_rep_val=k,i_vals=tmplist, error=error)
       DO  j = 1,SIZE(tmplist)
          jj = jj+1
       END DO
    END DO
    IF (jj<1) CALL stop_program(routineN, "Need at least 1 atom to use ddapc contraints")
    qs_control%ddapc_restraint_control%natoms=jj
    IF(ASSOCIATED(qs_control%ddapc_restraint_control%atoms))&
         DEALLOCATE(qs_control%ddapc_restraint_control%atoms)
    ALLOCATE(qs_control%ddapc_restraint_control%atoms(qs_control%ddapc_restraint_control%natoms))
    jj = 0
    DO k = 1,n_rep
       CALL section_vals_val_get(ddapc_restraint_section,"ATOMS",i_rep_val=k,i_vals=tmplist, error=error)
       DO  j = 1,SIZE(tmplist)
          jj = jj+1
          qs_control%ddapc_restraint_control%atoms (jj) = tmplist(j)
       END DO
    END DO

    IF(ASSOCIATED(qs_control%ddapc_restraint_control%coeff))&
         DEALLOCATE(qs_control%ddapc_restraint_control%coeff)
    ALLOCATE(qs_control%ddapc_restraint_control%coeff(qs_control%ddapc_restraint_control%natoms))
    qs_control%ddapc_restraint_control%coeff=1.0_dp

    CALL section_vals_val_get(ddapc_restraint_section,"COEFF", &
         n_rep_val=n_rep,error=error)
    jj = 0
    DO k = 1,n_rep
       CALL section_vals_val_get(ddapc_restraint_section,"COEFF",i_rep_val=k,r_vals=rtmplist, error=error)
       DO  j = 1,SIZE(rtmplist)
          jj = jj+1
          IF (jj>qs_control%ddapc_restraint_control%natoms) &
               CALL stop_program(routineN, "Need the same number of coeff as there are atoms ")
          qs_control%ddapc_restraint_control%coeff (jj) = rtmplist(j)
       END DO
    END DO
    IF (jj<qs_control%ddapc_restraint_control%natoms .AND. jj.NE.0) &
         CALL stop_program(routineN, "Need no or the same number of coeff as there are atoms ")

  END SUBROUTINE read_ddapc_section

! *****************************************************************************
!> \brief reads the input parameters needed for evaluating a becke weight population constraint
!> \author fschiff
! *****************************************************************************
  SUBROUTINE read_becke_section(qs_control,becke_section,error)

    TYPE(qs_control_type), INTENT(INOUT)     :: qs_control
    TYPE(section_vals_type), POINTER         :: becke_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_becke_section', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: j, jj, k, n_rep
    INTEGER, DIMENSION(:), POINTER           :: tmplist
    REAL(KIND=dp), DIMENSION(:), POINTER     :: rtmplist

    CALL section_vals_val_get(becke_section,"STRENGTH", &
         r_val=qs_control%becke_control%strength,error=error)
    CALL section_vals_val_get(becke_section,"TARGET", &
         r_val=qs_control%becke_control%target,error=error)
    CALL section_vals_val_get(becke_section,"ATOMS", &
         n_rep_val=n_rep,error=error)
    jj = 0
    DO k = 1,n_rep
       CALL section_vals_val_get(becke_section,"ATOMS",i_rep_val=k,i_vals=tmplist, error=error)
       DO  j = 1,SIZE(tmplist)
          jj = jj+1
       END DO
    END DO
    IF (jj<1) CALL stop_program(routineN, "Need at least 1 atom to use ddapc contraints")
    qs_control%becke_control%natoms=jj
    IF(ASSOCIATED(qs_control%becke_control%atoms))DEALLOCATE(qs_control%becke_control%atoms)
    ALLOCATE(qs_control%becke_control%atoms(qs_control%becke_control%natoms))
    jj = 0
    DO k = 1,n_rep
       CALL section_vals_val_get(becke_section,"ATOMS",i_rep_val=k,i_vals=tmplist, error=error)
       DO  j = 1,SIZE(tmplist)
          jj = jj+1
          qs_control%becke_control%atoms (jj) = tmplist(j)
       END DO
    END DO

    IF(ASSOCIATED(qs_control%becke_control%coeff))&
         DEALLOCATE(qs_control%becke_control%coeff)
    ALLOCATE(qs_control%becke_control%coeff(qs_control%becke_control%natoms))
    qs_control%becke_control%coeff=1.0_dp

    CALL section_vals_val_get(becke_section,"COEFF", &
         n_rep_val=n_rep,error=error)
    jj = 0
    DO k = 1,n_rep
       CALL section_vals_val_get(becke_section,"COEFF",i_rep_val=k,r_vals=rtmplist, error=error)
       DO  j = 1,SIZE(rtmplist)
          jj = jj+1
          IF (jj>qs_control%becke_control%natoms) &
               CALL stop_program(routineN, "Need the same number of coeff as there are atoms ")
          qs_control%becke_control%coeff (jj) = rtmplist(j)
       END DO
    END DO
    IF (jj<qs_control%becke_control%natoms .AND. jj.NE.0) &
         CALL stop_program(routineN, "Need no or the same number of coeff as there are atoms ")
  END SUBROUTINE read_becke_section

END MODULE cp_control_utils
