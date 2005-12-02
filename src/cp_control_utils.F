!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2004  CP2K developers group                          !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****s* cp2k/cp_control_utils [1.0] *
!!
!!   NAME
!!     cp_control_utils
!!
!!   FUNCTION
!!     Utilities to set up the control types
!!
!!   SOURCE
!******************************************************************************

MODULE cp_control_utils
  USE cp_control_types,                ONLY: dft_control_create,&
                                             dft_control_type,&
                                             gapw_control_create,&
                                             qs_control_type,&
                                             tddfpt_control_type,&
                                             wave_control_type
  USE input_constants,                 ONLY: do_loc_min,&
                                             do_loc_max,&
                                             do_loc_jacobi,&
                                             do_loc_crazy,&
                                             do_loc_direct,&
                                             do_loc_none,&
                                             op_loc_berry,&
                                             op_loc_boys,&
                                             op_loc_pipek,&
                                             state_loc_all,&
                                             state_loc_upton,&
                                             state_loc_fromitoj,&
                                             state_loc_list,&
                                             state_loc_none,&
                                             oe_none,&
                                             oe_saop,&
                                             oe_sic,&
                                             sic_ad,&
                                             sic_mauri_spz,&
                                             sic_mauri_us,&
                                             sic_none,&
                                             tddfpt_davidson,&
                                             tddfpt_excitations,&
                                             tddfpt_lanczos,&
                                             tddfpt_singlet,&
                                             tddfpt_triplet,&
                                             run_type_names
  USE input_constants,                 ONLY: xc_debug_new_routine,&
                                             xc_deriv_fd,&
                                             xc_deriv_nn10_smooth,&
                                             xc_deriv_nn50_smooth,&
                                             xc_deriv_pw,&
                                             xc_deriv_spline2,&
                                             xc_deriv_spline2_smooth,&
                                             xc_deriv_spline3,&
                                             xc_deriv_spline3_smooth,&
                                             xc_new_f_routine,&
                                             xc_old_f_routine,&
                                             xc_rho_nn10,&
                                             xc_rho_nn50,&
                                             xc_rho_no_smooth,&
                                             xc_rho_spline2_smooth,&
                                             xc_rho_spline3_smooth,&
                                             xc_test_lsd_f_routine
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_check,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_reset,&
                                             cp_error_type
  USE cp_keywords,                     ONLY: find_keyword,&
                                             keyword_type,&
                                             list_keywords
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: electronic_spectra_run,&
                                             energy_force_run,&
                                             energy_run,&
                                             geo_opt_run,&
                                             mon_car_run,&
                                             hmc_run,&
                                             ldmc_run,&
                                             mol_dyn_run,&
                                             qs_debug_run,&
                                             bsse_run
  USE input_constants,                 ONLY: quadrature_name,&
                                             do_pwgrid_spherical,&
                                             do_pwgrid_ns_fullspace,&
                                             do_pwgrid_ns_halfspace,&
                                             do_rsgrid_distributed,&
                                             do_rsgrid_replicated,&
                                             elstr_method_name,&
                                             do_method_gapw,& 
                                             do_method_gapw_xc,& 
                                             do_method_gpw,& 
                                             do_method_kg_pol,&
                                             do_method_kg_nopol,&
                                             do_method_kg_gpw,&
                                             do_method_eht,& 
                                             do_method_mndo,& 
                                             do_method_am1,& 
                                             do_method_pm3 
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_get_subs_vals2,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp,&
                                             dp_size, default_string_length
  USE memory_utilities,                ONLY: reallocate
  USE qs_loc_control,                  ONLY: read_loc_section
  USE qs_wf_history_types,             ONLY: wfi_linear_p_method_nr,&
                                             wfi_linear_ps_method_nr,&
                                             wfi_linear_wf_method_nr,&
                                             wfi_ps_method_nr,&
                                             wfi_use_guess_method_nr,&
                                             wfi_use_prev_p_method_nr,&
                                             wfi_use_prev_rho_r_method_nr,&
                                             wfi_use_prev_wf_method_nr, &
                                             wfi_frozen_method_nr, &
                                             wfi_aspc
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE xas_control,                     ONLY: read_xas_control,&
                                             xas_control_type
  USE xc_derivatives,                  ONLY: xc_functional_get_info,&
       xc_functionals_get_needs
  USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_setall,&
                                             xc_rho_cflags_type
  USE input_cp2k_dft, only: xc_functionals_expand

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: moduleN = "cp_control_utils"

! *** Public subroutines ***

  PUBLIC :: read_dft_control,&
            read_mgrid_section,&
            read_qs_section,&
            read_tddfpt_control,&
            write_dft_control,&
            write_qs_control

!!***
! *****************************************************************************

CONTAINS

  SUBROUTINE read_dft_control(dft_control,globenv, error)
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE read_dft_control (MODULE dft_types)"

    INTEGER                                  :: excitations, i, run_id
    LOGICAL                                  :: failure, l_param, xc_fun_set
    LOGICAL                                  :: was_present
    TYPE(cp_error_type)                      :: sub_error
    TYPE(section_vals_type), POINTER         :: dft_section, tmp_section, &
                                                xc_fun_section, &
                                                xc_grid_section, xc_section
    CHARACTER(len=default_string_length) :: basis_set_file_name, potential_file_name
    REAL(kind=dp) :: density_cut, gradient_cut, tau_cut
    TYPE(xc_rho_cflags_type) :: needs
         

!   ---------------------------------------------------------------------------

    failure=.FALSE.
    was_present = .FALSE.
    NULLIFY(dft_section,tmp_section,xc_fun_section,xc_grid_section,xc_section)

    CALL dft_control_create(dft_control, error=error)
    dft_section => section_vals_get_subs_vals(globenv%input_file,"FORCE_EVAL",error=error)
    dft_section => section_vals_get_subs_vals(dft_section,"DFT",error=error)
    xc_section => section_vals_get_subs_vals(dft_section,"XC",error=error)
    xc_grid_section => section_vals_get_subs_vals(xc_section,"XC_GRID",error=error)
    xc_fun_section => section_vals_get_subs_vals(xc_section,"XC_FUNCTIONAL",&
         error=error)
    
    
    CALL section_vals_val_get(dft_section,"LSD",l_val=l_param,error=error)
    IF (l_param) THEN
       dft_control % nspins = 2
    ELSE
       dft_control % nspins = 1
    END IF
    needs=xc_functionals_get_needs(xc_fun_section,lsd=l_param,&
         add_basic_components=.TRUE.,error=error)
    dft_control%use_kinetic_energy_density=(needs%tau_spin.or.needs%tau)



    CALL section_vals_val_get(dft_section,"RESTRICTED",l_val=l_param,error=error)
    dft_control % restricted = l_param
    
    dft_control % do_tddfpt_calculation = .FALSE.
    CALL section_vals_val_get(dft_section,"CHARGE",i_val=dft_control % charge,&
         error=error)
    CALL section_vals_val_get(dft_section,"MULTIPLICITY",&
         i_val=dft_control % multiplicity,error=error)

    dft_control % sic_method_id = sic_none
    dft_control % sic_scaling_a = 1.0_dp
    dft_control % sic_scaling_b = 1.0_dp

    dft_control % xas_estate    = -1

    dft_control % derive_function_rho = .FALSE.          !FM to rm
    dft_control % integration_grid % pw_default = .TRUE. !FM to rm
    dft_control % integration_grid % pw_cutoff  = -1._dp !FM to rm
    dft_control % integration_grid % pw_mesh    = -1_dp  !FM to rm

    CALL section_vals_val_get(dft_section,"BASIS_SET_FILE_NAME",&
         c_val=basis_set_file_name , error=error)
    CALL section_vals_val_get(dft_section,"POTENTIAL_FILE_NAME",&
         c_val=potential_file_name , error=error)

    tmp_section => section_vals_get_subs_vals(globenv%input_file,"GLOBAL",&
         error=error)
    CALL section_vals_val_get(tmp_section,"RUN_TYPE",i_val=run_id,error=error)
    SELECT CASE (run_id)!FM to rm, ugly!!
    CASE (energy_run,mon_car_run,energy_force_run,electronic_spectra_run,bsse_run)
      dft_control%forces = .FALSE.
    CASE (mol_dyn_run,hmc_run,ldmc_run,geo_opt_run,qs_debug_run)
      dft_control%forces = .TRUE.
    CASE DEFAULT
      CALL stop_program(routine,"unknown run type "//TRIM(run_type_names(globenv%run_type_id)))
    END SELECT
 
!   *** Read the input section ***
    tmp_section => section_vals_get_subs_vals(dft_section,"sic",&
         error=error)
    CALL section_vals_val_get(tmp_section,"SIC_METHOD",&
         i_val=dft_control%sic_method_id,error=error)
    CALL section_vals_val_get(tmp_section,"SIC_SCALING_A",&
         r_val=dft_control%sic_scaling_a,error=error)
    CALL section_vals_val_get(tmp_section,"SIC_SCALING_B",&
         r_val=dft_control%sic_scaling_b,error=error)
    

    CALL section_vals_val_get(dft_section,"EXCITATIONS",&
         i_val=excitations,error=error)
    dft_control % do_tddfpt_calculation = (excitations==tddfpt_excitations)

    xc_section => section_vals_get_subs_vals(globenv%input_file,&
                   "FORCE_EVAL%DFT%XC", error=error)
!FM    CALL read_xc_control(xc_section, dft_control%xc_control, dft_control%nspins,&
!FM         dft_control%use_kinetic_energy_density, error=error)

    
!   *** Read the input section ***
    tmp_section => section_vals_get_subs_vals(dft_section,"xas",&
         error=error)
    call section_vals_get(tmp_section,explicit=was_present,error=error)

    if (was_present) then

      dft_control % do_xas_calculation = .TRUE.

      CALL read_xas_control(dft_control%xas_control,globenv,error=error)

    else

      dft_control % do_xas_calculation = .FALSE.

    end if

  END SUBROUTINE read_dft_control

! *****************************************************************************
  SUBROUTINE read_mgrid_section(qs_control,globenv,error)

    TYPE(qs_control_type), INTENT(INOUT)     :: qs_control
    TYPE(global_environment_type), &
      INTENT(INOUT)                          :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE read_mgrid_control "

    LOGICAL                                  :: failure,multigrid_set
    INTEGER                                  :: j, igrid_level, istat, jj, k, n_rep, ngrid_level
    INTEGER, DIMENSION(:), POINTER           :: tmplist
    REAL(dp)                                 :: cutoff
    REAL(dp), DIMENSION(:), POINTER          :: cutofflist
    TYPE(section_vals_type), POINTER         :: interp_section, mgrid_section

!   ---------------------------------------------------------------------------

    failure=.FALSE.
    NULLIFY(interp_section,mgrid_section,cutofflist)
    mgrid_section => section_vals_get_subs_vals(globenv%input_file,"FORCE_EVAL",error=error)
    mgrid_section => section_vals_get_subs_vals(mgrid_section,"DFT",error=error)
    mgrid_section => section_vals_get_subs_vals(mgrid_section,"MGRID",error=error)


    CALL section_vals_val_get(mgrid_section,"NGRIDS",i_val=ngrid_level,error=error)
    CALL section_vals_val_get(mgrid_section,"MULTIGRID_SET",l_val=multigrid_set,error=error)

    ALLOCATE(qs_control%e_cutoff(ngrid_level),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"qs_control%e_cutoff",&
                                             ngrid_level*dp_size)
    CALL section_vals_val_get(mgrid_section,"CUTOFF",r_val=cutoff,error=error)
    qs_control%cutoff = MAX(0.01_dp,cutoff) 
    CALL section_vals_val_get(mgrid_section,"PROGRESSION_FACTOR",r_val=qs_control%progression_factor,error=error)
    CALL section_vals_val_get(mgrid_section,"COMMENSURATE",l_val=qs_control%commensurate_mgrids,error=error)
    CALL section_vals_val_get(mgrid_section,"REALSPACE",l_val=qs_control%realspace_mgrids,error=error)
    CALL section_vals_val_get(mgrid_section,"REL_CUTOFF",r_val=qs_control%relative_cutoff)
    
    IF(multigrid_set) THEN
      ! Read the values from input
      IF (qs_control%commensurate_mgrids) THEN
          CALL stop_program(routine,&
               "Do not specify cutoffs for the commensurate grids (NYI)")
      ENDIF

      CALL section_vals_val_get(mgrid_section,"MULTIGRID_CUTOFF",&
           r_vals=cutofflist,error=error)
      IF(ASSOCIATED(cutofflist)) THEN
         CALL cp_assert(SIZE(cutofflist,1) == ngrid_level,&
              cp_failure_level,cp_assertion_failed,routine,&
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
         CALL stop_program (routine,"Multi grids not ordered")
      END IF
    END DO

    interp_section => section_vals_get_subs_vals(mgrid_section,"INTERPOLATOR",error=error)
    ! I don't kow what this is, help is welcome


  END SUBROUTINE read_mgrid_section

  SUBROUTINE read_qs_section(qs_control,globenv,error)

    TYPE(qs_control_type), INTENT(INOUT)     :: qs_control
    TYPE(global_environment_type), &
      INTENT(INOUT)                          :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE read_qs_section "

    LOGICAL                                  :: failure, was_present
    INTEGER                                  :: j, jj, k, n_rep
    INTEGER, DIMENSION(:), POINTER           :: tmplist
    REAL(dp) :: value
    TYPE(section_vals_type), POINTER         :: qs_section,loc_section,&
                                                mull_section,& 
                                                se_section
    
!   ---------------------------------------------------------------------------

    failure=.FALSE.
    was_present = .FALSE.
    NULLIFY(qs_section,loc_section,mull_section,se_section)

    qs_section => section_vals_get_subs_vals(globenv%input_file,"FORCE_EVAL",error=error)
    qs_section => section_vals_get_subs_vals(qs_section,"DFT",error=error)
    qs_section => section_vals_get_subs_vals(qs_section,"QS",error=error)
    loc_section => section_vals_get_subs_vals(qs_section,"LOCALIZE",error=error)
    mull_section => section_vals_get_subs_vals(qs_section,"MULLIKEN", error=error)
    se_section => section_vals_get_subs_vals(qs_section,"SE", error=error)

    ! Setup all defaults values and overwrite input parameters
    CALL section_vals_val_get(qs_section,"EPS_CORE_CHARGE",r_val=qs_control%eps_core_charge,error=error)
    CALL section_vals_val_get(qs_section,"EPS_GVG_RSPACE",r_val=qs_control%eps_gvg_rspace,error=error)
    CALL section_vals_val_get(qs_section,"EPS_PGF_ORB",r_val=qs_control%eps_pgf_orb,error=error)
    qs_control%eps_pgf_aux = qs_control%eps_pgf_orb
    CALL section_vals_val_get(qs_section,"EPS_PPL",r_val=qs_control%eps_ppl,error=error)
    CALL section_vals_val_get(qs_section,"EPS_PPNL",r_val=qs_control%eps_ppnl,error=error)
    CALL section_vals_val_get(qs_section,"EPS_RHO_RSPACE",r_val=qs_control%eps_rho_rspace,error=error)
    qs_control%eps_rho_gspace = qs_control%eps_rho_rspace
    !
    CALL section_vals_val_get(qs_section,"EPS_DEFAULT",n_rep_val=n_rep,error=error)
    IF (n_rep /=0) THEN
       CALL section_vals_val_get(qs_section,"EPS_DEFAULT",r_val=value,error=error)
       qs_control%eps_core_charge = value/100.0_dp
       qs_control%eps_pgf_orb = SQRT(value)
       qs_control%eps_pgf_aux = SQRT(value)
       qs_control%eps_ppl = value/100.0_dp
       qs_control%eps_ppnl = value/100.0_dp
       qs_control%eps_rho_gspace = value
       qs_control%eps_rho_rspace = value
       qs_control%eps_gvg_rspace = SQRT(value)       
    END IF
    !
    CALL section_vals_val_get(qs_section,"EPS_RHO",n_rep_val=n_rep,error=error)
    IF (n_rep /=0) THEN
       CALL section_vals_val_get(qs_section,"EPS_RHO",r_val=qs_control%eps_rho_gspace,error=error)
       qs_control%eps_rho_rspace = qs_control%eps_rho_gspace
    END IF

    ! Gapw
    CALL gapw_control_create(qs_control%gapw_control)

    CALL section_vals_val_get(qs_section,"EPSFIT",r_val=qs_control%gapw_control%eps_fit,error=error)
    CALL section_vals_val_get(qs_section,"EPSISO",r_val=qs_control%gapw_control%eps_iso,error=error)
    CALL section_vals_val_get(qs_section,"EPSRHO0",r_val=qs_control%gapw_control%eps_Vrho0,error=error)
    CALL section_vals_val_get(qs_section,"ALPHA0_HARD",r_val=qs_control%gapw_control%alpha0_hard,error=error)
    qs_control%gapw_control%lrho1_eq_lrho0 = .FALSE.
    qs_control%gapw_control%alpha0_hard_from_input = .FALSE.
    IF(qs_control%gapw_control%alpha0_hard/=0.0_dp) qs_control%gapw_control%alpha0_hard_from_input = .TRUE.
    CALL section_vals_val_get(qs_section,"ALPHA0_SOFT",r_val=qs_control%gapw_control%alpha0_soft,error=error)
    qs_control%gapw_control%alpha0_soft_from_input = .FALSE.
    IF(qs_control%gapw_control%alpha0_soft/=0.0_dp) qs_control%gapw_control%alpha0_soft_from_input = .TRUE.
    CALL section_vals_val_get(qs_section,"MAX_RAD_LOCAL",r_val=qs_control%gapw_control%max_rad_local,error=error)

    ! Logicals
    CALL section_vals_val_get(qs_section,"MAP_PAA",l_val=qs_control%map_paa,error=error)
    CALL section_vals_val_get(qs_section,"MAP_CONSISTENT",l_val=qs_control%map_consistent,error=error)
    CALL section_vals_val_get(qs_section,"MULLIKEN_RESTRAINT",l_val=qs_control%mulliken_restraint,error=error)
    CALL section_vals_val_get(qs_section,"LOCALIZE_WFN",l_val=qs_control%localized_wfn,error=error)

    ! Integers gapw
    CALL section_vals_val_get(qs_section,"LMAXN1",i_val=qs_control%gapw_control%lmax_sphere,error=error)
    CALL section_vals_val_get(qs_section,"LMAXN0",i_val=qs_control%gapw_control%lmax_rho0,error=error)
    CALL section_vals_val_get(qs_section,"LADDN0",i_val=qs_control%gapw_control%ladd_rho0,error=error)
    CALL section_vals_val_get(qs_section,"QUADRATURE",i_val=qs_control%gapw_control%quadrature,error=error)
    qs_control%gapw_control%quadrature_name=quadrature_name(qs_control%gapw_control%quadrature)
    
    ! Integers grids
    CALL section_vals_val_get(qs_section,"PW_GRID",i_val=qs_control%pw_grid_opt_id,error=error)
    SELECT CASE (qs_control%pw_grid_opt_id)
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

    CALL section_vals_val_get(qs_section,"RS_GRID",i_val=qs_control%rs_distributed_id,error=error)

    !Integers extrapolation
    CALL section_vals_val_get(qs_section,"EXTRAPOLATION",i_val=qs_control%wf_interpolation_method_nr,error=error)
    CALL section_vals_val_get(qs_section,"EXTRAPOLATION_ORDER",i_val=qs_control%wf_extrapolation_order,error=error)
    
    !Method
    CALL section_vals_val_get(qs_section,"METHOD",i_val=qs_control%method_id,error=error)
    qs_control %method = elstr_method_name(qs_control %method_id)
    qs_control%gapw = .FALSE.
    qs_control%gapw_xc = .FALSE.
    qs_control%gpw = .FALSE.
    qs_control%pao = .FALSE.
    qs_control%polarization = .FALSE.
    qs_control%semi_empirical = .FALSE.
    IF(qs_control %method_id==do_method_gapw) qs_control%gapw = .TRUE.
    IF(qs_control %method_id==do_method_gapw_xc) qs_control%gapw_xc = .TRUE.
    IF(qs_control %method_id==do_method_gpw) qs_control%gpw = .TRUE.
    IF(qs_control %method_id >= do_method_eht .AND. qs_control %method_id <= do_method_pm3)&
                              qs_control%semi_empirical = .TRUE.
    
    IF(qs_control %method_id == do_method_kg_pol) qs_control%polarization = .TRUE.
    
    IF(qs_control%localized_wfn) THEN
      CALL read_loc_section(qs_control%localized_wfn_control,loc_section,error=error)
    END IF

!    CALL read_mulliken_section()
    CALL section_vals_val_get(mull_section,"STRENGTH",r_val=qs_control%mulliken_restraint_control%strength,error=error)
    CALL section_vals_val_get(mull_section,"TARGET",r_val=qs_control%mulliken_restraint_control%target,error=error)
    CALL section_vals_val_get(mull_section,"NATOMS",i_val=qs_control%mulliken_restraint_control%natoms,error=error)
    IF(qs_control%mulliken_restraint_control%natoms .GT. 0) THEN
      ALLOCATE(qs_control%mulliken_restraint_control%atoms(qs_control%mulliken_restraint_control%natoms))
      CALL section_vals_val_get(mull_section,"ATOMS",n_rep_val=n_rep,error=error)
      jj = 0
      do k = 1,n_rep
       CALL section_vals_val_get(mull_section,"ATOMS",i_rep_val=k,i_vals=tmplist, error=error)
       DO  j = 1,SIZE(tmplist) 
         jj = jj+1
         qs_control%mulliken_restraint_control%atoms (jj) = tmplist(j) 
       END DO 
      end do
    ELSE
      IF(ASSOCIATED(qs_control%mulliken_restraint_control%atoms)) &
          DEALLOCATE(qs_control%mulliken_restraint_control%atoms)
      NULLIFY(qs_control%mulliken_restraint_control%atoms)
    END IF

!    CALL read_se_control_section()
    CALL section_vals_val_get(se_section,"orthogonal_basis",l_val=qs_control%se_control%orthogonal_basis,error=error)
    CALL section_vals_val_get(se_section,"EWALD",l_val=qs_control%se_control%ewald,error=error)
    CALL section_vals_val_get(se_section,"DELTA",r_val=qs_control%se_control%delta,error=error)
    CALL section_vals_val_get(se_section,"ALPHA_GAUSS",r_val=qs_control%se_control%alpha_gauss,error=error)
    CALL section_vals_val_get(se_section,"EPS_GAUSS",r_val=qs_control%se_control%eps_gauss,error=error)
    CALL section_vals_val_get(se_section,"RC_EXCHANGE",r_val=qs_control%se_control%rc_exchange,error=error)
    CALL section_vals_val_get(se_section,"RC_COULOMB",r_val=qs_control%se_control%rc_coulomb,error=error)

    ! do not print dft section for semi empirical methods
    IF(qs_control%semi_empirical) THEN
      globenv%print%dft_control_parameters=.FALSE.
    END IF


  END SUBROUTINE read_qs_section

! *****************************************************************************

  SUBROUTINE read_tddfpt_control(t_control,nspins,globenv, error)
    TYPE(tddfpt_control_type)                :: t_control
    INTEGER, INTENT(IN)                      :: nspins
    TYPE(global_environment_type)            :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = "read_tddfpt_control", &
      routineP = moduleN//"/"//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure, kenergy_den, was_present
    TYPE(section_vals_type), POINTER         :: t_section, sic_section, xc_section
    TYPE(section_vals_type), POINTER         :: tmp_section, xc_grid, xc_fun


    failure=.FALSE.
    kenergy_den = .FALSE.
    nullify(sic_section, t_section, xc_section, xc_grid, xc_fun)
    t_section => section_vals_get_subs_vals(globenv%input_file,"FORCE_EVAL",error=error)
    t_section => section_vals_get_subs_vals(t_section,"DFT",error=error)
    t_section => section_vals_get_subs_vals(t_section,"TDDFPT",error=error)
    
    CALL section_vals_val_get(t_section,"CONVERGENCE", &
         r_val=t_control%tolerance,error=error)

    CALL section_vals_val_get(t_section,"NEV", &
         i_val=t_control%n_ev,error=error)

    CALL section_vals_val_get(t_section,"MAX_KV", &
         i_val=t_control%max_kv,error=error)

    CALL section_vals_val_get(t_section,"RESTARTS", &
         i_val=t_control%n_restarts,error=error)

    CALL section_vals_val_get(t_section,"NREORTHO", &
         i_val=t_control%n_reortho,error=error)

    CALL section_vals_val_get(t_section,"RES_ETYPE", &
         i_val=t_control%res_etype,error=error)

    CALL section_vals_val_get(t_section,"DIAG_METHOD", &
         i_val=t_control%diag_method,error=error)

    CALL section_vals_val_get(t_section,"KERNEL", &
         l_val=t_control%do_kernel,error=error)

    CALL section_vals_val_get(t_section,"LSD_SINGLETS", &
         l_val=t_control%lsd_singlets,error=error)

    CALL section_vals_val_get(t_section,"INVERT_S", &
         l_val=t_control%invert_S,error=error)

    CALL section_vals_val_get(t_section,"PRECOND", &
         l_val=t_control%precond,error=error)

    CALL section_vals_val_get(t_section,"OE_CORR", &
         i_val=t_control%oe_corr,error=error)

    t_control%use_kinetic_energy_density = .FALSE.
    xc_section=>section_vals_get_subs_vals(t_section,"XC",error=error)
    call section_vals_get(xc_section,explicit=was_present,error=error)
    if (was_present) then
!FM       CALL read_xc_control(xc_section, t_control%xc_control, nspins,&
!FM                            t_control%use_kinetic_energy_density, error=error)
    end if

    sic_section=>section_vals_get_subs_vals(t_section,"SIC",error=error)
    CALL section_vals_val_get(sic_section,"SIC_METHOD",&
         i_val=t_control%sic_method_id,error=error)
    CALL section_vals_val_get(sic_section,"SIC_SCALING_A",&
         r_val=t_control%sic_scaling_a,error=error)
    CALL section_vals_val_get(sic_section,"SIC_SCALING_B",&
         r_val=t_control%sic_scaling_b,error=error)

  END SUBROUTINE read_tddfpt_control

!   ***************************************************************************

  SUBROUTINE write_dft_control(dft_control,globenv, error)

    !   Purpose: Write the DFT control parameters to the output unit.

    !   ***************************************************************************

    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER              :: routineN = "write_dft_control"

    CHARACTER(len=160)                       :: reference
    CHARACTER(len=20)                        :: tmpStr
    INTEGER                                  :: i, output_unit, ifun
    LOGICAL                                  :: failure, lsd
    TYPE(section_vals_type), POINTER         :: xc_section,xc_fun_section,&
         xc_fun
    REAL(kind=dp) :: density_cut, gradient_cut, tau_cut, density_smooth_cut_range

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%ionode) RETURN
    failure=.FALSE.

    IF (globenv%print%dft_control_parameters) THEN

       output_unit = globenv%scr
       xc_section => section_vals_get_subs_vals(globenv%input_file,&
            "FORCE_EVAL%DFT%XC",error=error)
       CALL section_vals_val_get(globenv%input_file,"FORCE_EVAL%DFT%LSD",&
            l_val=lsd,error=error)
       IF (dft_control%forces) THEN
          WRITE (UNIT=output_unit,FMT="(/,T2,A,T78,A)")&
               "DFT| Force calculation","YES"
       ELSE
          WRITE (UNIT=output_unit,FMT="(/,T2,A,T78,A)")&
               "DFT| Force calculation"," NO"
       END IF

       IF ( lsd ) THEN
          WRITE ( output_unit, '( A, T60, A )' ) " DFT| ","spin polarized system"
       ELSE
          WRITE ( output_unit, '( A, T60, A, T77, I4 )' ) " DFT| ","sets of spin states:",&
               dft_control % nspins
       END IF

       IF ( dft_control % restricted ) THEN
          WRITE ( output_unit, '( A, T60, A )' ) " DFT| ","      spin restricted"
       ENDIF
       
       SELECT CASE ( dft_control%sic_method_id )
       CASE ( sic_none )
          WRITE ( output_unit, '( A, T60, A )' ) " DFT| ","               NO SIC"
       CASE ( sic_mauri_spz )
          WRITE ( output_unit, '( A, T60, A )' ) " DFT| ","        SPZ/MAURI SIC"
       CASE ( sic_mauri_us )
          WRITE ( output_unit, '( A, T60, A )' ) " DFT| ","         US/MAURI SIC"
       CASE ( sic_ad )
          WRITE ( output_unit, '( A, T60, A )' ) " DFT| ","               AD SIC"
       CASE DEFAULT
            ! fix throughout the cp2k for this option
            CALL stop_program(routineN,"SIC option unknown")
       END SELECT

       IF ( dft_control%sic_method_id .NE. sic_none) THEN
          WRITE ( output_unit, '( A, T50, A, F17.12 )' ) " DFT| ","SIC SCALING A",dft_control%sic_scaling_a
          WRITE ( output_unit, '( A, T50, A, F17.12 )' ) " DFT| ","SIC SCALING B",dft_control%sic_scaling_b
       ENDIF

       IF ( dft_control % charge /= 0 ) THEN
          WRITE ( output_unit, '( A, T60, A, T77, I4 )' ) " DFT| ","charge:",&
               dft_control % charge
       END IF

       CALL section_vals_val_get(xc_section,"density_cutoff",&
            r_val=density_cut,error=error)
       CALL section_vals_val_get(xc_section,"gradient_cutoff",&
            r_val=gradient_cut,error=error)
       CALL section_vals_val_get(xc_section,"tau_cutoff",&
            r_val=tau_cut,error=error)
       CALL section_vals_val_get(xc_section,"density_smooth_cutoff_range",&
            r_val=density_smooth_cut_range,error=error)
       WRITE (UNIT=output_unit,FMT="(T2,A,T67,E14.5)")&
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
       IF ( dft_control % derive_function_rho )  &
            WRITE ( output_unit, '( A, T79, A )' ) &
            " DFT| Density function derivative ","ON"

       IF ( dft_control % integration_grid % pw_default ) THEN
          WRITE ( output_unit, '( A, T74, A )' ) " DFT| XC mesh (PW density)", "default"
       ELSE
          IF ( dft_control % integration_grid % pw_mesh ( 1 ) > 0 ) &
               WRITE ( output_unit, '( A, T51, 3I10 )' ) &
               " DFT| XC mesh (PW density)", dft_control % integration_grid % pw_mesh
          IF ( dft_control % integration_grid % pw_cutoff > 0 ) &
               WRITE ( output_unit, '( A, T67, F10.1, "a.u." )' ) &
               " DFT| XC mesh (PW density)", dft_control % integration_grid % pw_cutoff
       END IF

       IF (.NOT.failure) THEN
          xc_fun_section => section_vals_get_subs_vals(xc_section,&
               "XC_FUNCTIONAL",error=error)
          CALL section_vals_val_get(xc_section,"FUNCTIONAL_ROUTINE",&
               c_val=tmpStr,error=error)
          WRITE (output_unit, fmt="(' FUNCTIONAL| ROUTINE=',a)")&
               TRIM(tmpStr)
          CALL xc_functionals_expand(xc_fun_section,error=error)
          ifun=0
          DO
             ifun=ifun+1
             xc_fun => section_vals_get_subs_vals2(xc_fun_section,i_section=ifun,&
                  error=error)
             IF (.NOT.ASSOCIATED(xc_fun)) EXIT
             CALL xc_functional_get_info(xc_fun,lsd=lsd,&
                  reference=reference, error=error)
             WRITE (output_unit, fmt="(' FUNCTIONAL| ',a,':')")&
                  TRIM(xc_fun%section%name)
             WRITE (output_unit, fmt="(' FUNCTIONAL| ',a)")&
                  TRIM(reference)
          END DO
       END IF
    END IF
  END SUBROUTINE write_dft_control

! *****************************************************************************

  SUBROUTINE write_qs_control(qs_control,globenv)

!   Purpose: Write the DFT control parameters to the output unit.

!   ***************************************************************************

    TYPE(qs_control_type), INTENT(IN)        :: qs_control
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    INTEGER                                  :: igrid_level, ngrid_level, &
                                                nstep,output_unit
    

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%ionode) RETURN

    IF (globenv%print%dft_control_parameters) THEN
      ngrid_level = SIZE(qs_control%e_cutoff)
      output_unit = globenv%scr
      WRITE (UNIT=output_unit,FMT="(/,T2,A,T71,A)")&
        "QS| Method:",ADJUSTR(qs_control%method)
      IF ( qs_control% pw_grid_opt % spherical ) THEN
         WRITE (UNIT=output_unit,FMT="(T2,A,T61,A)")&
            "QS| Density plane wave grid type"," SPHERICAL HALFSPACE"
      ELSE IF ( qs_control% pw_grid_opt % fullspace ) THEN
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
           qs_control%gapw_control%eps_iso
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
           qs_control%gapw_control%eps_iso
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
      IF (qs_control%localized_wfn) THEN
         nstep = qs_control%localized_wfn_control%each_step
         WRITE (UNIT=output_unit,FMT="(T2,A)")&
           "QS| Localization: The spread relative to a set of orbitals is computed"
         WRITE (UNIT=output_unit,FMT="(T2,A,T55,I5,A)")&
           "                  at the end of a scf loop every ",nstep," iteration steps"
         SELECT CASE ( qs_control%localized_wfn_control%set_of_states )
         CASE(state_loc_all)
           WRITE (UNIT=output_unit,FMT="(T2,A)")& 
           "                Set:  All the occupied orbitals"
         CASE(state_loc_upton)
           WRITE (UNIT=output_unit,FMT="(T2,A,T73,I8)")& 
           "                Set:  All the orbital with index lower than ",&
            qs_control%localized_wfn_control%nloc_states(1)
         CASE(state_loc_fromitoj)
           WRITE (UNIT=output_unit,FMT="(T2,A,T73,2I8)")& 
           "                Set:  All the orbital with index between ",&
            qs_control%localized_wfn_control%lu_bound_states(1,1),&
            qs_control%localized_wfn_control%lu_bound_states(2,1)
         CASE(state_loc_list)
           WRITE (UNIT=output_unit,FMT="(T2,A)")& 
           "                Set:  All the orbital with index in the list:"
           WRITE (UNIT=output_unit,FMT=*) qs_control%localized_wfn_control%loc_states
         CASE DEFAULT
           WRITE (UNIT=output_unit,FMT="(T2,A)")& 
           "                Set:  None"
         END SELECT
         SELECT CASE (qs_control%localized_wfn_control%operator_type )
         CASE(op_loc_berry) 
           WRITE (UNIT=output_unit,FMT="(T2,A)")&
           "                Spread: defined by the Berry phase operator "
         CASE(op_loc_boys) 
           WRITE (UNIT=output_unit,FMT="(T2,A)")&
           "                Spread: defined by the Boys operator"
         CASE DEFAULT
            WRITE (UNIT=output_unit,FMT="(T2,A)")&
           "                Spread: defined by the Pipek operator"
         END SELECT
         
         SELECT CASE (qs_control%localized_wfn_control%localization_method)
         CASE(do_loc_jacobi)
            WRITE (UNIT=output_unit,FMT="(T2,A)")&
             "                Optimize  the spread by the jacobi algorithm"
         CASE(do_loc_crazy)
            WRITE (UNIT=output_unit,FMT="(T2,A)")&
             "                Optimize  the spread by the crazy angle algorithm"
         CASE(do_loc_direct)
            WRITE (UNIT=output_unit,FMT="(T2,A)")&
             "                Optimize  the spread by the direct algorithm"
         CASE(do_loc_none)
            WRITE (UNIT=output_unit,FMT="(T2,A)")&
             "                No minimization or maximization is performed "
         END SELECT
      END IF
    END IF

  END SUBROUTINE write_qs_control

! *****************************************************************************

!!****f* cp_control_types/xc_get_rho_smooth_label [1.0] *
!!
!!   NAME
!!     xc_get_rho_smooth_label
!!
!!   FUNCTION
!!     returns a string that describes the smoothing of rho
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - xc_rho_smooth_id: the id that represent the smoothing
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     05.2002 created [fawzi]
!!
!!*** **********************************************************************
FUNCTION xc_get_rho_smooth_label(xc_rho_smooth_id,error) RESULT(res)
    INTEGER, INTENT(in)                      :: xc_rho_smooth_id
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
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
!***************************************************************************

!!****f* cp_control_types/xc_get_deriv_method_label [1.0] *
!!
!!   NAME
!!     xc_get_deriv_method_label
!!
!!   FUNCTION
!!     returns a string that describes the derivative used in the xc 
!!     calculation
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - xc_deriv_method_id: the id that represent the derivative method
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     06.2002 created [fawzi]
!!
!!*** **********************************************************************
FUNCTION xc_get_deriv_method_label(xc_deriv_method_id,error) RESULT(res)
    INTEGER, INTENT(in)                      :: xc_deriv_method_id
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
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
  CASE (xc_deriv_fd)
     res="FD"
  CASE default
     WRITE (res,"('UKN',i6)") xc_deriv_method_id
  END SELECT
END FUNCTION xc_get_deriv_method_label
!***************************************************************************

!!****f* cp_control_types/xc_get_routine_label [1.0] *
!!
!!   NAME
!!     xc_get_routine_label
!!
!!   FUNCTION
!!     returns a string that describes the functional routine used in the xc 
!!     calculation
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - xc_functional_routine_id: the id that represent the functional 
!!       routine
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!*** **********************************************************************
FUNCTION xc_get_routine_label(xc_functional_routine_id,error) RESULT(res)
    INTEGER, INTENT(in)                      :: xc_functional_routine_id
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
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
!***************************************************************************

END MODULE cp_control_utils
