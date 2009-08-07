!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief routines that build the Kohn-Sham matrix (i.e calculate the coulomb
!>      and xc parts
!> \par History
!>      05.2002 moved from qs_scf (see there the history) [fawzi]
!>      JGH [30.08.02] multi-grid arrays independent from density and potential
!>      10.2002 introduced pools, uses updated rho as input,
!>              removed most temporary variables, renamed may vars,
!>              began conversion to LSD [fawzi]
!>      10.2004 moved calculate_w_matrix here [Joost VandeVondele]
!>              introduced energy derivative wrt MOs [Joost VandeVondele]
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE qs_ks_methods

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE cell_types,                      ONLY: cell_type
  USE cp_control_types,                ONLY: becke_restraint_type,&
                                             dft_control_type
  USE cp_ddapc,                        ONLY: cp_ddapc_apply_CD,&
                                             qs_ks_ddapc
  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                             cp_fm_gemm,&
                                             cp_fm_scale_and_add,&
                                             cp_fm_symm,&
                                             cp_fm_transpose,&
                                             cp_fm_upper_to_full
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: copy_fm_to_sm,&
                                             copy_sm_to_fm,&
                                             cp_sm_fm_multiply,&
                                             cp_sm_plus_fm_fm_t
  USE dft_plus_u,                      ONLY: plus_u
  USE efield_utils,                    ONLY: efield_potential
  USE ep_qs_types,                     ONLY: ep_qs_type
  USE f77_blas
  USE harris_env_types,                ONLY: harris_env_type
  USE hartree_local_2centers,          ONLY: integrate_vhg0_g0
  USE hartree_local_3centers,          ONLY: integrate_a_vhg0_b
  USE hartree_local_methods,           ONLY: Vh_1c_gg_integrals
  USE hfx_communication,               ONLY: scale_and_add_fock_to_ks_matrix
  USE hfx_derivatives,                 ONLY: derivatives_four_center
  USE hfx_energy_potential,            ONLY: integrate_four_center
  USE hfx_ri_methods,                  ONLY: hfx_ri_energy_potential
  USE input_constants,                 ONLY: &
       do_adiabatic_hybrid_mcy3, do_adiabatic_model_pade, do_ppl_grid, &
       sic_ad, sic_eo, sic_list_all, sic_list_unpaired, sic_mauri_spz, &
       sic_mauri_us, sic_none, use_aux_fit_basis_set, use_orb_basis_set, &
       xc_none
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kahan_sum,                       ONLY: accurate_sum
  USE kg_environment_types,            ONLY: get_kg_env,&
                                             kg_environment_type
  USE kg_gpw_correction,               ONLY: kg_gpw_ekin_mol
  USE kg_gpw_fm_mol_types,             ONLY: kg_fm_p_type
  USE kg_gpw_pw_env_types,             ONLY: kg_sub_pw_env_type
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_grids,                        ONLY: pw_grid_compare
  USE pw_methods,                      ONLY: pw_axpy,&
                                             pw_copy,&
                                             pw_integrate_function,&
                                             pw_scale,&
                                             pw_transfer,&
                                             pw_zero
  USE pw_poisson_methods,              ONLY: pw_poisson_solve
  USE pw_poisson_types,                ONLY: pw_poisson_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_p_type,&
                                             pw_pool_retain,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_create,&
                                             pw_p_type,&
                                             pw_release,&
                                             pw_type
  USE qs_charges_types,                ONLY: qs_charges_type
  USE qs_collocate_density,            ONLY: calculate_rho_elec,&
                                             calculate_wavefunction
  USE qs_core_energies,                ONLY: calculate_ecore
  USE qs_dftb_matrices,                ONLY: build_dftb_ks_matrix
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_gapw_densities,               ONLY: prepare_gapw_den
  USE qs_integrate_potential,          ONLY: integrate_ppl_rspace,&
                                             integrate_v_core_rspace,&
                                             integrate_v_rspace
  USE qs_ks_apply_restraints,          ONLY: qs_ks_becke_restraint,&
                                             qs_ks_mulliken_restraint,&
                                             qs_ks_s2_restraint
  USE qs_ks_atom,                      ONLY: update_ks_atom
  USE qs_ks_qmmm_methods,              ONLY: qmmm_calculate_energy,&
                                             qmmm_modify_hartree_pot
  USE qs_ks_scp_methods,               ONLY: scp_calculate_gpw_energy,&
                                             scp_calculate_qmmm_energy,&
                                             scp_modify_hartree_pot
  USE qs_ks_scp_types,                 ONLY: qs_ks_scp_env_type
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_mo_methods,                   ONLY: calculate_density_matrix
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             mo_set_type
  USE qs_rho0_ggrid,                   ONLY: integrate_vhg0_rspace
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE qs_vxc_atom,                     ONLY: calculate_vxc_atom
  USE scp_energy,                      ONLY: scp_qs_energies
  USE scp_energy_types,                ONLY: scp_energy_type
  USE scp_environment_types,           ONLY: get_scp_env,&
                                             scp_environment_type
  USE scp_hartree_1center,             ONLY: integrate_a_vhscp_b
  USE se_fock_matrix,                  ONLY: build_se_fock_matrix
  USE sparse_matrix_types,             ONLY: cp_sm_scale_and_add,&
                                             deallocate_matrix,&
                                             real_matrix_p_type,&
                                             real_matrix_type,&
                                             replicate_matrix_structure,&
                                             set_matrix,&
                                             transfer_matrix
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_types,                    ONLY: virial_type
  USE wfn_fitting_methods,             ONLY: fit_mo_coeffs,&
                                             merge_mo_derivs
  USE wfn_fitting_types,               ONLY: wfn_fitting_type
  USE xc,                              ONLY: xc_exc_calc,&
                                             xc_vxc_pw_create1
  USE xc_adiabatic_methods,            ONLY: rescale_MCY3_pade
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  INTERFACE calculate_w_matrix
    MODULE PROCEDURE calculate_w_matrix_1,&
                     calculate_w_matrix_2,&
                     calculate_w_matrix_roks
                     
  END INTERFACE

  LOGICAL, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ks_methods'
  INTEGER, SAVE, PRIVATE :: last_ks_id_nr=0

  PUBLIC :: qs_ks_create, qs_ks_did_change, qs_vxc_create, calc_rho_tot_gspace, &
       qs_ks_update_qs_env, qs_ks_build_kohn_sham_matrix, calculate_w_matrix
!***

CONTAINS

! *****************************************************************************
!> \brief allocates and initializes the given ks_env.
!> \param ks_env the ks env to be initialized
!> \param qs_env the qs environment
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      05.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE qs_ks_create(ks_env, qs_env, error)
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ks_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: auxbas_grid, handle, nspins, &
                                                stat
    LOGICAL                                  :: failure
    TYPE(cell_type), POINTER                 :: cell
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool

    CALL timeset(routineN,handle)

    failure=.FALSE.
    CPPreconditionNoFail(.NOT.ASSOCIATED(ks_env),cp_failure_level,routineP,error)
    ALLOCATE(ks_env, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    IF (.not.failure) THEN
       NULLIFY( pw_env, auxbas_pw_pool,&
            pw_env, cell,pw_pools)
       CALL get_qs_env(qs_env=qs_env,&
            dft_control=dft_control,&
            pw_env=pw_env, cell=cell,error=error)
       CALL pw_env_get(pw_env,auxbas_pw_pool=auxbas_pw_pool,pw_pools=pw_pools,&
            auxbas_grid=auxbas_grid,error=error)

       nspins=dft_control%nspins

       ks_env%s_mstruct_changed=.TRUE.
       ks_env%rho_changed=.TRUE.
       ks_env%potential_changed=.TRUE.
       ks_env%forces_up_to_date=.FALSE.
       ks_env%n_evals=0
       ks_env%ref_count=1
       last_ks_id_nr=last_ks_id_nr+1
       ks_env%id_nr=last_ks_id_nr

    END IF
    IF (.NOT.failure) THEN
       ks_env%auxbas_pw_pool=>auxbas_pw_pool
       CALL pw_pool_retain(ks_env%auxbas_pw_pool,error=error)
       CALL pw_pool_create_pw(auxbas_pw_pool,ks_env%v_hartree_rspace%pw,&
            use_data=REALDATA3D, in_space=REALSPACE,error=error)
    END IF

    CALL timestop(handle)

  END SUBROUTINE qs_ks_create

! *****************************************************************************
!> \brief routine where the real calculations are made: the
!>      KS matrix is calculated
!> \param ks_env the ks_env that old all the temporary objects that
!>        the calculation of the KS matrix needs
!> \param qs_env the qs_env to update
!> \param calculate_forces if true calculate the quantities needed
!>        to calculate the forces. Defaults to false.
!> \param just_energy if true updates the energies but not the
!>        ks matrix. Defaults to false
!> \param ks the KS matrix that is created
!> \param rho the density, must be up to date
!> \param rho_xc the soft density, only for gapw_xc , must be up to date
!> \param energy the place where energies are stored
!> \param qs_charges the place where the grids charges are stored
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      make rho, energy and qs_charges optional, defaulting
!>      to qs_env components?
!> \par History
!>      06.2002 moved from qs_scf to qs_ks_methods, use of ks_env
!>              new did_change scheme [fawzi]
!>      10.2002 introduced pools, uses updated rho as input, LSD [fawzi]
!>      10.2004 build_kohn_sham matrix now also computes the derivatives
!>              of the total energy wrt to the MO coefs, if instructed to
!>              do so. This appears useful for orbital dependent functionals
!>              where the KS matrix alone (however this might be defined)
!>               does not contain the info to construct this derivative.
!> \author Matthias Krack
! *****************************************************************************
  SUBROUTINE qs_ks_build_kohn_sham_matrix(ks_env,qs_env,ks_matrix,&
                                          rho,energy,calculate_forces,&
                                          just_energy,rho_xc,kg_gpw,&
                                          kg_fm_set,kg_sub_pw_env,&
                                          error)

    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(qs_energy_type), POINTER            :: energy
    LOGICAL, INTENT(in)                      :: calculate_forces, just_energy
    TYPE(qs_rho_type), OPTIONAL, POINTER     :: rho_xc
    LOGICAL, INTENT(IN), OPTIONAL            :: kg_gpw
    TYPE(kg_fm_p_type), OPTIONAL, POINTER    :: kg_fm_set
    TYPE(kg_sub_pw_env_type), OPTIONAL, &
      POINTER                                :: kg_sub_pw_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_ks_build_kohn_sham_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER :: adiabatic_functional, adiabatic_model, handle, ikind, irep, &
      ispin, k, k2, n, n_rep_hf, nspins, output_unit, stat
    LOGICAL :: ddapc_restraint_is_spin, do_adiabatic_rescaling, do_ddapc, &
      do_ep, do_hfx, do_hfx_ri, do_ppl, explicit_potential, failure, gapw, &
      gapw_xc, hfx_treat_lsd_in_core, ionode, my_kg_gpw, scp_dft, &
      uniform_occupation, use_virial
    REAL(dp) :: adiabatic_lambda, adiabatic_omega, scale_dDFA, scale_ddW0, &
      scale_dEx1, scale_dEx2, total_energy_xc
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: hf_energy
    REAL(KIND=dp)                            :: ddapc_order_p, ecore_ppl, &
                                                mulliken_order_p, s2_order_p
    REAL(KIND=dp), DIMENSION(3, 3)           :: h_stress
    REAL(KIND=dp), DIMENSION(:), POINTER     :: occupation_numbers
    TYPE(becke_restraint_type), POINTER      :: becke
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mo_derivs, mo_derivs_aux_fit
    TYPE(cp_fm_type), POINTER                :: hc_tmp, mo_coeff, &
                                                mo_coeff_aux_fit
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(ep_qs_type), POINTER                :: ep_qs_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array, mos, mos_aux_fit
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type) :: ep_pot_r_coeff, ep_rho_g, ep_rho_r, rho_tot_gspace, &
      v_efield_rspace, v_hartree_gspace, v_sic_rspace, v_spin_ddapc_rest_r
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: v_rspace_new, &
                                                v_rspace_new_aux_fit, &
                                                v_tau_rspace, &
                                                v_tau_rspace_aux_fit
    TYPE(pw_p_type), POINTER                 :: vppl_rspace
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(pw_type), POINTER                   :: ep_pot_g, ep_pot_r
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_ks_scp_env_type), POINTER        :: ks_scp_env
    TYPE(qs_rho_type), POINTER               :: my_rho_xc
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_dv, matrix_h, &
                                                matrix_p_aux_fit, matrix_s, &
                                                my_rho
    TYPE(scp_energy_type), POINTER           :: scp_energy
    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(section_vals_type), POINTER :: adiabatic_rescaling_section, &
      hfx_ri_section, hfx_section, hfx_sections, input, scf_section
    TYPE(virial_type), POINTER               :: virial
    TYPE(wfn_fitting_type), POINTER          :: wfn_fitting_env

    failure=.FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(matrix_h,matrix_s,dft_control,auxbas_pw_pool, pw_pools,&
         pw_env, cell, logger, v_rspace_new, v_tau_rspace,&
         input, poisson_env, scf_section,my_rho,becke, v_rspace_new_aux_fit,&
         v_tau_rspace_aux_fit)

    CPPrecondition(ASSOCIATED(ks_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ks_env%ref_count>0,cp_failure_level,routineP,error,failure)

    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    ! Check for a KG_GPW calculation: KE for total rho and for molecular rho
    my_kg_gpw = .FALSE.
    IF(PRESENT(kg_gpw)) my_kg_gpw = kg_gpw

    NULLIFY (matrix_dv)

    CALL get_qs_env(qs_env=qs_env,&
         dft_control=dft_control,&
         matrix_h=matrix_h,&
         matrix_s=matrix_s,&
         pw_env=pw_env,&
         cell=cell,&
         para_env=para_env,&
         input=input,&
         virial=virial,&
         matrix_dv=matrix_dv,&
         ep_qs_env=ep_qs_env,&
         error=error)
! Define logicals for SCP
    scp_dft =  dft_control % scp 

    use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)
    do_ep=ASSOCIATED(ep_qs_env)
    IF (do_ep) do_ep=ep_qs_env%ep_active.AND.calculate_forces.and.&
         ASSOCIATED(ep_qs_env%dH_coeffs)

    hfx_sections => section_vals_get_subs_vals(input,"DFT%XC%HF",error=error)
    CALL section_vals_get(hfx_sections,explicit=do_hfx,error=error)
    IF( do_hfx ) THEN
      CALL section_vals_val_get(hfx_sections, "TREAT_LSD_IN_CORE", l_val=hfx_treat_lsd_in_core,&
                                i_rep_section=1,error=error) 
      NULLIFY(hfx_ri_section)
      hfx_ri_section => section_vals_get_subs_vals(hfx_sections,"HFX_RI",error=error)
      CALL section_vals_get(hfx_ri_section,explicit=do_hfx_ri,error=error)
    END IF
    adiabatic_rescaling_section => section_vals_get_subs_vals(input,"DFT%XC%ADIABATIC_RESCALING",error=error)
    CALL section_vals_get(adiabatic_rescaling_section,explicit=do_adiabatic_rescaling,error=error)

    scf_section => section_vals_get_subs_vals(input,"DFT%SCF",error=error)
    nspins=dft_control%nspins
    CPPrecondition(ASSOCIATED(matrix_h),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(matrix_s),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(rho),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(pw_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(SIZE(ks_matrix)>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(rho%rho_r_valid,cp_failure_level,routineP,error,failure)
    CPPrecondition(rho%rho_g_valid,cp_failure_level,routineP,error,failure)

    ! Setup the possible usage of DDAPC charges
    do_ddapc = qs_env%dft_control%qs_control%ddapc_restraint.OR.&
               qs_env%cp_ddapc_ewald%do_decoupling.OR.&
               qs_env%cp_ddapc_ewald%do_qmmm_periodic.OR.&
               qs_env%cp_ddapc_ewald%do_solvation

    ! Check for GAPW method : additional terms for local densities
    gapw    = dft_control%qs_control%gapw
    gapw_xc = dft_control%qs_control%gapw_xc
    IF(gapw_xc .AND. gapw) CALL cp_assert( .FALSE. , cp_failure_level,cp_assertion_failed,routineP,&
                           " GAPW and GAPW_XC are not compatible",error,failure)
    IF((gapw_xc .OR. gapw) .AND. use_virial) &
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="Virial for GAPW not debugged!!", &
            error=error, error_level=cp_warning_level)

    do_ppl = dft_control%qs_control%do_ppl_method == do_ppl_grid
    IF ( do_ppl ) THEN
      CPPrecondition(.NOT.gapw,cp_failure_level,routineP,error,failure)
      CALL get_qs_env ( qs_env = qs_env, vppl = vppl_rspace, error = error )
    END IF

    IF(gapw_xc) THEN
      CPPrecondition(ASSOCIATED(rho_xc),cp_failure_level,routineP,error,failure)
      CPPrecondition(rho_xc%rho_r_valid,cp_failure_level,routineP,error,failure)
      CPPrecondition(rho_xc%rho_g_valid,cp_failure_level,routineP,error,failure)
      my_rho_xc => rho_xc
    ELSE
      NULLIFY(my_rho_xc)
    END IF

    ! gets the tmp grids
    IF (.NOT. failure) THEN
       CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,&
                       pw_pools=pw_pools, poisson_env=poisson_env,error=error)
    END IF

    ! ***  Prepare densities for gapw ***
    IF(gapw .OR. gapw_xc) THEN
       CALL prepare_gapw_den(qs_env,do_rho0=(.NOT.gapw_xc),error=error)
    ENDIF

    IF (.NOT. failure) THEN

       ! *** calculate the hartree potential on the pw density ***
       ! *** Hartree contributions ***

       IF (.NOT.failure) THEN
         CALL pw_pool_create_pw(auxbas_pw_pool,&
                                 v_hartree_gspace%pw, &
                                 use_data=COMPLEXDATA1D,&
                                 in_space=RECIPROCALSPACE,&
                                 error=error)
         CALL pw_pool_create_pw(auxbas_pw_pool,&
                                 rho_tot_gspace%pw,&
                                 use_data=COMPLEXDATA1D,&
                                 in_space=RECIPROCALSPACE,&
                                 error=error)
       END IF

       IF (BTEST(cp_print_key_should_output(logger%iter_info,scf_section,&
            "PRINT%DETAILED_ENERGY",error=error),cp_p_file).AND.(.NOT.gapw).AND.(.NOT.gapw_xc)) THEN
          CALL pw_zero(rho_tot_gspace%pw, error=error)
          CALL calc_rho_tot_gspace(rho_tot_gspace, qs_env, rho, skip_nuclear_density=.TRUE.,error=error)
          CALL pw_poisson_solve(poisson_env,rho_tot_gspace%pw,energy%e_hartree,&
               v_hartree_gspace%pw,error=error)
          CALL pw_zero(rho_tot_gspace%pw, error=error)
          CALL pw_zero(v_hartree_gspace%pw, error=error)
       END IF

       !CJM get the total density in g-space [ions + electrons]
       CALL calc_rho_tot_gspace(rho_tot_gspace, qs_env, rho,error=error)

       output_unit = cp_print_key_unit_nr(logger,scf_section,"PRINT%TOTAL_DENSITIES",&
            extension=".scfLog",error=error)
       CALL print_densities(qs_env, rho_tot_gspace, rho, my_rho_xc, output_unit, error)
       CALL cp_print_key_finished_output(output_unit,logger,scf_section,&
            "PRINT%TOTAL_DENSITIES", error=error)

       ! CJM getting the Hartree energy and Hartree potential.  Also getting the stress tensor
       ! CJM from the Hartree term if needed.  No nuclear force information here
       IF (use_virial .AND. calculate_forces) THEN
         h_stress(:,:) = 0.0_dp
         CALL pw_poisson_solve(poisson_env,rho_tot_gspace%pw,energy%hartree,&
                             v_hartree_gspace%pw,h_stress=h_stress,error=error)
         virial%pv_virial = virial%pv_virial + h_stress/REAL(para_env%num_pe,dp)
       ELSE
         CALL pw_poisson_solve(poisson_env,rho_tot_gspace%pw,energy%hartree,&
                             v_hartree_gspace%pw,error=error)
       END IF
       IF ( scp_dft ) THEN
          ! Compute the contribution rho_scp * v_hartree_gspace.  v_hartree_gpace is
          ! potential due to the ions + electrons.  
          ! Get the SCP environment
          CALL get_qs_env ( qs_env = qs_env, scp_env = scp_env, error = error )
          CALL scp_calculate_gpw_energy ( scp_env, v_hartree_gspace, error )
       END IF 

       ! In case decouple periodic images and/or apply restraints to charges
       IF (do_ddapc) THEN
          CALL qs_ks_ddapc(qs_env, auxbas_pw_pool, rho_tot_gspace, v_hartree_gspace,&
               v_spin_ddapc_rest_r, energy, calculate_forces, ks_env, ks_matrix,    &
               just_energy, ddapc_restraint_is_spin, explicit_potential, error)
       ELSE
          explicit_potential      = .FALSE.
          ddapc_restraint_is_spin = .FALSE.
          IF (.NOT. just_energy) THEN
             CALL pw_transfer(v_hartree_gspace%pw, ks_env%v_hartree_rspace%pw, error=error)
             CALL pw_scale(ks_env%v_hartree_rspace%pw, ks_env%v_hartree_rspace%pw%pw_grid%dvol,&
                  error=error)
          END IF
       END IF
       CALL pw_pool_give_back_pw(auxbas_pw_pool,v_hartree_gspace%pw,error=error)
       CALL pw_pool_give_back_pw(auxbas_pw_pool,rho_tot_gspace%pw,error=error)

       IF(dft_control%apply_efield_field)THEN
          CALL pw_pool_create_pw(auxbas_pw_pool,v_efield_rspace%pw,&
            use_data=REALDATA3D, in_space=REALSPACE,error=error)
          CALL efield_potential(qs_env,v_efield_rspace,calculate_forces,error)
          CALL pw_scale(v_efield_rspace%pw, v_efield_rspace%pw%pw_grid%dvol, error=error)
       END IF


       ! SIC
       CALL calc_v_sic_rspace(v_sic_rspace,energy,qs_env,dft_control,rho,poisson_env,&
            just_energy,calculate_forces,auxbas_pw_pool,error=error)

       IF ( gapw ) CALL Vh_1c_gg_integrals(qs_env,energy%hartree_1c,error=error)

       !***** Check if becke potential is needed to constrain charges *****
       CALL qs_ks_becke_restraint(qs_env,auxbas_pw_pool,calculate_forces,matrix_s,becke,error)
       
       ! *** Add the QM/MM potential ***
       IF (qs_env%qmmm) THEN
          output_unit=cp_print_key_unit_nr(logger,input,"QMMM%PRINT%PROGRAM_RUN_INFO",&
               extension=".qmmmLog",error=error)
          IF (output_unit>0) &
               WRITE (UNIT=output_unit,FMT="(T3,A)")&
               "Adding QM/MM electrostatic potential to the Kohn-Sham potential.."
           CALL qmmm_calculate_energy (qs_env=qs_env,&
                                       rho=rho%rho_r,&
                                       v_qmmm=qs_env%ks_qmmm_env%v_qmmm_rspace,&
                                       qmmm_energy=energy%qmmm_el,&
             ! also take into account the scaling (in addition to the volume element)
                                       error=error)
           CALL cp_print_key_finished_output(output_unit,logger,input,&
                "QMMM%PRINT%PROGRAM_RUN_INFO", error=error)

           IF (.NOT.just_energy) THEN
              CALL qmmm_modify_hartree_pot(v_hartree=ks_env%v_hartree_rspace,&
                   v_qmmm=qs_env%ks_qmmm_env%v_qmmm_rspace,&
                   error=error)
           END IF
       END IF
! ****SCP        *** Add the polarization potential ***
       IF ( scp_dft ) THEN
! Get the KS_SCP environment
          CALL get_scp_env ( scp_env = scp_env, ks_scp_env = ks_scp_env, error = error )
          IF (output_unit>0) &
               WRITE (UNIT=output_unit,FMT="(T3,A)")&
               "Adding SCP  potential to the Kohn-Sham potential.."
          IF (qs_env % qmmm ) THEN
            CALL scp_calculate_qmmm_energy ( scp_env, qs_env%ks_qmmm_env%v_qmmm_rspace, error )
          END IF
          IF (.NOT.just_energy) THEN
! Needs up-to-date KS_SCP_ENV
            CALL scp_modify_hartree_pot(v_hartree=ks_env%v_hartree_rspace,&
                                        v_scp=ks_scp_env%v_scp_rspace,&
                                        error=error)
          END IF
       END IF
! ****SCP   

       !! calculate the density matrix for the fitted mo_coeffs
       IF( dft_control%do_wfn_fitting ) THEN       
         NULLIFY(mo_coeff_aux_fit)
         CALL get_qs_env(qs_env,mos_aux_fit=mos_aux_fit,error=error)
         CALL get_qs_env(qs_env,mos=mos,wfn_fitting_env=wfn_fitting_env,error=error)
         CALL fit_mo_coeffs(wfn_fitting_env, para_env, &
                            qs_env%matrix_s_aux_fit, qs_env%matrix_s_aux_fit_vs_orb,&
                            mos, mos_aux_fit, ks_env%s_mstruct_changed, error)
         CALL set_qs_env(qs_env,wfn_fitting_env=wfn_fitting_env,error=error)

         DO ispin=1,nspins
           matrix_p_aux_fit =>qs_env%rho_aux_fit%rho_ao
           CALL calculate_density_matrix(mos_aux_fit(ispin)%mo_set, &
                                         matrix_p_aux_fit(ispin)%matrix,error=error)
           CALL calculate_rho_elec(matrix_p=qs_env%rho_aux_fit%rho_ao(ispin)%matrix,&
                rho=qs_env%rho_aux_fit%rho_r(ispin),&
                rho_gspace=qs_env%rho_aux_fit%rho_g(ispin),&
                total_rho=qs_env%rho_aux_fit%tot_rho_r(ispin),&
                qs_env=qs_env,soft_valid=.FALSE.,&
                basis_set_id=use_aux_fit_basis_set, error=error)
         END DO
         qs_env%rho_aux_fit%rho_r_valid=.TRUE.
         qs_env%rho_aux_fit%rho_g_valid=.TRUE.
       END IF

       ! only activate stress calculation if 
       IF (use_virial .AND. calculate_forces) virial%pv_calculate = .TRUE.

       IF( .NOT. do_adiabatic_rescaling ) THEN
         ! *** calculate the xc potential on the pw density ***
         ! *** associates v_rspace_new if the xc potential needs to be computed.
! CJM getting the XC potential
         !! If we do wavefunction fitting, we need the vxc_potential in the auxiliary basis set
         IF( dft_control%do_wfn_fitting ) THEN
           CALL qs_vxc_create(qs_env=qs_env,&
              vxc_rho=v_rspace_new_aux_fit, vxc_tau=v_tau_rspace_aux_fit, exc=energy%exc_aux_fit, &
              just_energy=just_energy, gapw_xc=gapw_xc, basis_set_id=use_aux_fit_basis_set,&
              error=error)
         END IF
         CALL qs_vxc_create(qs_env=qs_env,&
              vxc_rho=v_rspace_new, vxc_tau=v_tau_rspace, exc=energy%exc, &
              just_energy=just_energy, gapw_xc=gapw_xc, error=error)
         IF (gapw .OR. gapw_xc) THEN
            CALL calculate_vxc_atom(qs_env,just_energy,error)
         END IF
       ELSE
         !! If we perform adiabatic rescaling, the xc potential has to be scaled by the xc- and
         !! HFX-energy. Thus, let us first calculate the energy 
         CALL qs_vxc_create(qs_env=qs_env,&
              vxc_rho=v_rspace_new, vxc_tau=v_tau_rspace, exc=energy%exc, &
              just_energy=.TRUE., gapw_xc=gapw_xc, error=error)
         IF (gapw .OR. gapw_xc) THEN
           CALL calculate_vxc_atom(qs_env,.TRUE.,error=error)
         END IF

       END IF

       IF (use_virial .AND. calculate_forces) THEN
         virial%pv_virial = virial%pv_virial - virial%pv_xc 
         virial%pv_calculate = .FALSE.
       ENDIF
      
       ! *** Initialize the auxiliary ks matrix to zero if required
       IF( dft_control%do_wfn_fitting ) THEN
         DO ispin = 1,nspins
           CALL set_matrix(qs_env%matrix_ks_aux_fit(ispin)%matrix,0.0_dp) 
         END DO
       END IF
       ! *** Add Hartree-Fock contribution if required ***
       IF ( do_hfx ) THEN
         DO ispin = 1,nspins
           CALL set_matrix(ks_matrix(ispin)%matrix,0.0_dp)
         END DO
         CALL section_vals_get(hfx_sections,n_repetition=n_rep_hf,error=error)
         !! set hf exchange energy to zero
         energy%ex = 0.0_dp

         IF ( do_hfx_ri ) THEN
           IF( calculate_forces ) THEN
             !! initalize force array to zero
             CALL get_qs_env(qs_env=qs_env, force=force, error=error)
             DO ikind = 1,SIZE(force)
               force(ikind)%hfx_ri(:,:) = 0.0_dp
             END DO
           END IF
!deb
           CALL hfx_ri_energy_potential (qs_env,calculate_forces,error)
!deb
         ELSE
           IF( calculate_forces ) THEN
             !! initalize force array to zero
             CALL get_qs_env(qs_env=qs_env, force=force, error=error)
             DO ikind = 1,SIZE(force)
               force(ikind)%fock_4c(:,:) = 0.0_dp
             END DO
           END IF
           ALLOCATE(hf_energy(n_rep_hf), STAT=stat) 
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           DO irep = 1,n_rep_hf
             IF( .NOT. do_adiabatic_rescaling) THEN
               IF( dft_control%do_wfn_fitting ) THEN
                 IF ( .NOT. hfx_treat_lsd_in_core ) THEN 
                   CALL integrate_four_center(qs_env, qs_env%matrix_ks_aux_fit,energy,&
                                              qs_env%rho_aux_fit,hfx_sections,&
                                              para_env,ks_env%s_mstruct_changed,irep,.TRUE.,&
                                              ispin=1, error=error)
                 ELSE
                   DO ispin = 1,nspins
                     CALL integrate_four_center(qs_env, qs_env%matrix_ks_aux_fit,energy,&
                                                qs_env%rho_aux_fit,hfx_sections,&
                                                para_env,ks_env%s_mstruct_changed,irep,.TRUE.,&
                                                ispin=ispin, error=error)
                   END DO
                 END IF
               ELSE
                 IF ( .NOT. hfx_treat_lsd_in_core ) THEN
                   CALL integrate_four_center(qs_env, ks_matrix ,energy,rho,hfx_sections,&
                                              para_env,ks_env%s_mstruct_changed,irep,.TRUE.,&
                                              ispin=1, error=error)
                 ELSE
                   DO ispin = 1,nspins
                     CALL integrate_four_center(qs_env, ks_matrix ,energy,rho,hfx_sections,&
                                                para_env,ks_env%s_mstruct_changed,irep,.TRUE.,&
                                                ispin=ispin, error=error)
                   END DO
                 END IF
               END IF
               IF( calculate_forces ) THEN
                 CPPrecondition(.NOT.use_virial,cp_failure_level,routineP,error,failure)
                 CALL derivatives_four_center(qs_env, ks_matrix,energy,rho,hfx_section,&
                                              para_env, irep, error=error)
               END IF
             ELSE
               IF( .NOT. hfx_treat_lsd_in_core ) THEN
                 !! we calculate everything but we do not add the potential in a first step
                 CALL integrate_four_center(qs_env, ks_matrix,energy,rho,hfx_sections,&
                                            para_env,ks_env%s_mstruct_changed,irep,.FALSE.,&
                                            ispin=1, error=error)
               ELSE
                 CALL cp_unimplemented_error(fromWhere=routineP, &
                       message="HFX_TREAT_LSD_IN_CORE not implemented for adiabatically rescaled hybrids", &
                       error=error, error_level=cp_failure_level)
               END IF
               hf_energy(irep) = energy%ex
               !! If required, the calculation of the forces will be done later
             END IF
             CALL pw_hfx(qs_env,energy,hfx_sections,poisson_env,auxbas_pw_pool,irep,error)
           END DO

         END IF 

         ! *** Add Core-Hamiltonian-Matrix ***
         DO ispin=1,nspins
            CALL cp_sm_scale_and_add(ks_matrix(ispin)%matrix,matrix_b=matrix_h(1)%matrix,error=error)
         END DO
       END IF 

       !! If we perform adiabatic rescaling we are now able to rescale the xc-potential
       IF( do_adiabatic_rescaling ) THEN
         CALL section_vals_val_get(adiabatic_rescaling_section, "FUNCTIONAL_TYPE", i_val=adiabatic_functional,error=error) 
         CALL section_vals_val_get(adiabatic_rescaling_section, "FUNCTIONAL_MODEL", i_val=adiabatic_model,error=error)
         CALL section_vals_val_get(adiabatic_rescaling_section, "LAMBDA", r_val=adiabatic_lambda,error=error)
         CALL section_vals_val_get(adiabatic_rescaling_section, "OMEGA", r_val=adiabatic_omega,error=error)
         SELECT CASE(adiabatic_functional)
         CASE (do_adiabatic_hybrid_mcy3)
           SELECT CASE(adiabatic_model)
           CASE(do_adiabatic_model_pade)
             CALL cp_assert( n_rep_hf == 2 , cp_failure_level,cp_assertion_failed,routineP,&
                            " For this kind of adiababatic hybrid functional 2 HF sections have to be provided. "//&
                            " Please check your input file.",&
                             error,failure)
             CALL rescale_MCY3_pade(qs_env, hf_energy, energy, adiabatic_lambda, &
                                    adiabatic_omega, scale_dEx1, scale_ddW0, scale_dDFA,&
                                    scale_dEx2, total_energy_xc, error)
              
             !! Scale and add Fock matrix to KS matrix
             IF(do_hfx) THEN
               CALL scale_and_add_fock_to_ks_matrix(para_env, qs_env, ks_matrix, 1 ,&
                                                    scale_dEx1, error)
               CALL scale_and_add_fock_to_ks_matrix(para_env, qs_env, ks_matrix, 2 ,&
                                                    scale_dEx2, error)
             END IF
             IF( calculate_forces ) THEN
               CPPrecondition(.NOT.use_virial,cp_failure_level,routineP,error,failure)
               !! we also have to scale the forces!!!!
               CALL derivatives_four_center(qs_env, ks_matrix,energy,rho,hfx_section,&
                                            para_env, 1, adiabatic_rescale_factor=scale_dEx1, error=error)
               CALL derivatives_four_center(qs_env, ks_matrix,energy,rho,hfx_section,&
                                            para_env, 2, adiabatic_rescale_factor=scale_dEx2, error=error)
             END IF

             !! Calculate vxc and rescale it
             CALL qs_vxc_create(qs_env=qs_env,&
                  vxc_rho=v_rspace_new, vxc_tau=v_tau_rspace, exc=energy%exc, &
                  just_energy=just_energy, gapw_xc=gapw_xc, adiabatic_rescale_factor=scale_dDFA, error=error)
             !! Calcualte vxc and rescale it
             IF (gapw .OR. gapw_xc) THEN
                CALL calculate_vxc_atom(qs_env,just_energy,adiabatic_rescale_factor=scale_dDFA,error=error) 
             END IF
             !! Hack for the total energy expression
             energy%ex = 0.0_dp 
             energy%exc1 = 0.0_dp
             energy%exc = total_energy_xc 

           END SELECT
         END SELECT            
       END IF

       IF(do_ppl .AND. calculate_forces)THEN
          CPPrecondition(.NOT.gapw,cp_failure_level,routineP,error,failure)
          DO ispin=1,nspins
            CALL integrate_ppl_rspace(rho%rho_r(ispin),qs_env,error=error)
          END DO
       END IF

       ! ***  Single atom contributions ***
       IF (.NOT. just_energy) THEN
!          IF (calculate_forces .AND. .NOT. gapw) THEN
          IF (calculate_forces ) THEN
             IF (do_ep) THEN
                CALL pw_pool_create_pw(auxbas_pw_pool,ep_rho_r%pw,&
                     use_data=REALDATA3D,in_space=REALSPACE,error=error)
                CALL pw_pool_create_pw(auxbas_pw_pool,ep_rho_g%pw,&
                     use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,error=error)
                CPAssert(nspins==1,cp_failure_level,routineP,error,failure)
                CALL calculate_rho_elec(matrix_p=ep_qs_env%dH_coeffs(1)%matrix,&
                     rho=ep_rho_r,rho_gspace=ep_rho_g, total_rho=ep_qs_env%tot_rho,&
                     qs_env=qs_env, error=error)

                CALL pw_pool_give_back_pw(auxbas_pw_pool,ep_rho_r%pw,error=error)
                CALL pw_pool_create_pw(auxbas_pw_pool,ep_pot_g,&
                     use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,error=error)

                CALL pw_poisson_solve(poisson_env,ep_rho_g%pw,ep_qs_env%core_energy,&
                     ep_pot_g,error=error)

                CALL pw_pool_give_back_pw(auxbas_pw_pool,ep_rho_g%pw,error=error)
                CALL pw_pool_create_pw(auxbas_pw_pool,ep_pot_r,&
                     use_data=REALDATA3D,in_space=REALSPACE,error=error)

                CALL pw_transfer(ep_pot_g,ep_pot_r, error=error)

                CALL pw_pool_give_back_pw(auxbas_pw_pool,ep_pot_g,error=error)

                ep_pot_r_coeff%pw => ep_pot_r
                CALL integrate_v_core_rspace(ep_pot_r_coeff, qs_env,error=error)

                CALL pw_pool_give_back_pw(auxbas_pw_pool,ep_pot_r,error=error)
             ELSE
! CJM getting nuclear force contribution from the core charge density
                CALL integrate_v_core_rspace(ks_env%v_hartree_rspace, qs_env,error=error)
             END IF
          END IF

          IF( .NOT.do_hfx ) THEN
            ! Initialize the Kohn-Sham matrix with the core Hamiltonian matrix
           
            DO ispin=1,nspins
             CALL transfer_matrix(matrix_h(1)%matrix,ks_matrix(ispin)%matrix,error=error)
            END DO
          END IF

          IF (ASSOCIATED(v_rspace_new)) THEN
             IF(gapw_xc) THEN
               ! not implemented (or at least not tested)
               CPPrecondition(dft_control%sic_method_id==sic_none,cp_failure_level,routineP,error,failure)

               DO ispin=1,nspins
                 !Only the xc potential, because it has to be integrated with the soft basis
                  v_rspace_new(ispin)%pw%cr3d  =&
                       v_rspace_new(ispin)%pw%pw_grid%dvol * &
                       v_rspace_new(ispin)%pw%cr3d

                 ! add the xc  part due to v_rspace soft
                 CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin),&
                      p=rho%rho_ao(ispin),h=ks_matrix(ispin),&
                      qs_env=qs_env, &
                      calculate_forces=calculate_forces,&
                      gapw=gapw_xc,error=error)

                 ! Now the Hartree potential to be integrated with the full basis
                 v_rspace_new(ispin)%pw%cr3d  =&
                      ks_env%v_hartree_rspace%pw%cr3d

                 ! add the hartree part due to v_rspace
                 CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin),&
                      p=rho%rho_ao(ispin),h=ks_matrix(ispin),&
                      qs_env=qs_env, &
                      calculate_forces=calculate_forces,&
                      gapw=gapw,error=error)

                 CALL pw_pool_give_back_pw(auxbas_pw_pool,v_rspace_new(ispin)%pw,&
                     error=error)

               END DO  ! ispin
            ELSE
!CJM add v_hartree + v_xc = v_rspace_new
               DO ispin=1,nspins
                  v_rspace_new(ispin)%pw%cr3d  =&
                       v_rspace_new(ispin)%pw%pw_grid%dvol * &
                       v_rspace_new(ispin)%pw%cr3d + &
                       ks_env%v_hartree_rspace%pw%cr3d
                  IF(explicit_potential)THEN
                     IF (ddapc_restraint_is_spin) THEN
                        IF (ispin==1) THEN
                           v_rspace_new(ispin)%pw%cr3d=v_rspace_new(ispin)%pw%cr3d &
                                +v_spin_ddapc_rest_r%pw%cr3d
                        ELSE
                           v_rspace_new(ispin)%pw%cr3d=v_rspace_new(ispin)%pw%cr3d &
                                -v_spin_ddapc_rest_r%pw%cr3d
                        ENDIF
                     ELSE
                        v_rspace_new(ispin)%pw%cr3d=v_rspace_new(ispin)%pw%cr3d &
                             +v_spin_ddapc_rest_r%pw%cr3d
                     END IF
                  END IF
                  IF(qs_env%dft_control%qs_control%becke_restraint)THEN
                     v_rspace_new(ispin)%pw%cr3d=v_rspace_new(ispin)%pw%cr3d &
                          +becke%becke_pot%pw%cr3d*qs_env%dft_control%qs_control%becke_control%strength
                  END IF
                  ! the efield contribution
                  IF(dft_control%apply_efield_field)THEN
                      v_rspace_new(ispin)%pw%cr3d= v_rspace_new(ispin)%pw%cr3d+&
                                                   v_efield_rspace%pw%cr3d
                  END IF
                  IF(do_ppl)THEN
                    CPPrecondition(.NOT.gapw,cp_failure_level,routineP,error,failure)
                    v_rspace_new(ispin)%pw%cr3d  = v_rspace_new(ispin)%pw%cr3d + &
                      vppl_rspace%pw%cr3d*vppl_rspace%pw%pw_grid%dvol
                  END IF
                  ! the electrostatic sic contribution
                  SELECT CASE (dft_control%sic_method_id)
                  CASE (sic_none)
                     !
                  CASE (sic_mauri_us,sic_mauri_spz )
                     IF (ispin==1) THEN
                        v_rspace_new(ispin)%pw%cr3d=v_rspace_new(ispin)%pw%cr3d &
                                                 -v_sic_rspace%pw%cr3d
                     ELSE
                        v_rspace_new(ispin)%pw%cr3d=v_rspace_new(ispin)%pw%cr3d &
                                                 +v_sic_rspace%pw%cr3d
                     ENDIF
                  CASE ( sic_ad )
                    v_rspace_new(ispin)%pw%cr3d=v_rspace_new(ispin)%pw%cr3d-v_sic_rspace%pw%cr3d
                  CASE ( sic_eo )
                     ! NOTHING TO BE DONE
                  END SELECT

                  ! add the part due to v_rspace
                  IF (do_ep) THEN
                     my_rho => ep_qs_env%dH_coeffs
!CJM Pointing my_rho to the density matrix rho_ao
                  ELSE
                     my_rho => rho%rho_ao
                  END IF
!CJM compute matrix elements of V_HXC including nuclear force information through matrix_dv
                  IF (ASSOCIATED(matrix_dv)) THEN
                     CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin),&
                          p=my_rho(ispin),h=ks_matrix(ispin),&
                          qs_env=qs_env, &
                          calculate_forces=calculate_forces,&
                          gapw=gapw,matrix_dv=matrix_dv(3*(ispin-1)+1:3*ispin),&
                          error=error)
                  ELSE
                     CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin),&
                          p=my_rho(ispin),h=ks_matrix(ispin),&
                          qs_env=qs_env, &
                          calculate_forces=calculate_forces,&
                          gapw=gapw, error=error)
                  END IF
                  CALL pw_pool_give_back_pw(auxbas_pw_pool,v_rspace_new(ispin)%pw,&
                        error=error)
               END DO
             END IF

             SELECT CASE (dft_control%sic_method_id)
             CASE (sic_none)
             CASE (sic_mauri_us,sic_mauri_spz, sic_ad )
                CALL pw_pool_give_back_pw(auxbas_pw_pool,v_sic_rspace%pw,error=error)
             END SELECT
             DEALLOCATE(v_rspace_new,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

          ELSE
             ! not implemented (or at least not tested)
             CPPrecondition(dft_control%sic_method_id==sic_none,cp_failure_level,routineP,error,failure)
             CPPrecondition(.NOT.ddapc_restraint_is_spin,cp_failure_level,routineP,error,failure)
             DO ispin=1,nspins
                ! the efield contribution
                IF(dft_control%apply_efield_field)THEN
                   ks_env%v_hartree_rspace%pw%cr3d= ks_env%v_hartree_rspace%pw%cr3d +&
                        v_efield_rspace%pw%cr3d
                END IF
                ! add only v_hartree
                 IF (do_ep) THEN
                    my_rho => ep_qs_env%dH_coeffs
                 ELSE
                    my_rho => rho%rho_ao
                 END IF
                 IF (ASSOCIATED(matrix_dv)) THEN
                    CALL integrate_v_rspace(v_rspace=ks_env%v_hartree_rspace,&
                         p=my_rho(ispin),h=ks_matrix(ispin),&
                         qs_env=qs_env, &
                         calculate_forces=calculate_forces,&
                         gapw=gapw,matrix_dv=matrix_dv(3*(ispin-1)+1:3*ispin),&
                         error=error)
                 ELSE
                    CALL integrate_v_rspace(v_rspace=ks_env%v_hartree_rspace,&
                         p=my_rho(ispin),h=ks_matrix(ispin),&
                         qs_env=qs_env, &
                         calculate_forces=calculate_forces,&
                         gapw=gapw, error=error)

                 END IF
             END DO
          END IF ! ASSOCIATED(v_rspace_new)

          IF (ASSOCIATED(v_tau_rspace)) THEN
             DO ispin=1,nspins
                 v_tau_rspace(ispin)%pw%cr3d =&
                      v_tau_rspace(ispin)%pw%pw_grid%dvol*&
                      v_tau_rspace(ispin)%pw%cr3d

                 IF (ASSOCIATED(matrix_dv)) THEN
                    CALL integrate_v_rspace(v_rspace=v_tau_rspace(ispin),&
                         p=rho%rho_ao(ispin),h=ks_matrix(ispin),&
                         qs_env=qs_env,&
                         calculate_forces=calculate_forces,compute_tau=.TRUE., &
                         gapw=gapw,matrix_dv=matrix_dv(3*(ispin-1)+1:3*ispin),&
                         error=error)
                 ELSE
                    CALL integrate_v_rspace(v_rspace=v_tau_rspace(ispin),&
                         p=rho%rho_ao(ispin),h=ks_matrix(ispin),&
                         qs_env=qs_env,&
                         calculate_forces=calculate_forces,compute_tau=.TRUE., &
                         gapw=gapw,&
                         error=error)
                 END IF
                 CALL pw_pool_give_back_pw(auxbas_pw_pool,v_tau_rspace(ispin)%pw,&
                      error=error)

              END DO
              DEALLOCATE(v_tau_rspace, stat=stat)
              CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ENDIF

          ! Add contributions from wfn_fitting if requested
          IF (dft_control%do_wfn_fitting) THEN
            IF( ASSOCIATED(v_rspace_new_aux_fit)) THEN
              DO ispin=1,nspins
                ! Calculate the xc potential
                v_rspace_new_aux_fit(ispin)%pw%cr3d  =&
                   v_rspace_new_aux_fit(ispin)%pw%pw_grid%dvol * &
                   v_rspace_new_aux_fit(ispin)%pw%cr3d
                ! Add potential to ks_matrix aux_fit
                CALL integrate_v_rspace(v_rspace=v_rspace_new_aux_fit(ispin),&
                                        p=qs_env%rho_aux_fit%rho_ao(ispin),&
                                        h=qs_env%matrix_ks_aux_fit(ispin),&
                                        qs_env=qs_env, &
                                        calculate_forces=calculate_forces,&
                                        gapw=gapw_xc,&
                                        basis_set_id=use_aux_fit_basis_set,&
                                        error=error)
                 CALL pw_pool_give_back_pw(auxbas_pw_pool,v_rspace_new_aux_fit(ispin)%pw,&
                     error=error)
              END DO
              DEALLOCATE(v_rspace_new_aux_fit,stat=stat)
              CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
            END IF
            ! Clean up v_tau_rspace_aux_fit, which is actually not needed
            IF( ASSOCIATED(v_tau_rspace_aux_fit)) THEN
              DO ispin=1,nspins
                CALL pw_pool_give_back_pw(auxbas_pw_pool,v_tau_rspace_aux_fit(ispin)%pw,&
                       error=error)
              END DO
              DEALLOCATE(v_tau_rspace_aux_fit, stat=stat)
              CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
            END IF  
          END IF

          
       END IF  ! .NOT. just energy

       IF (explicit_potential) THEN
           CALL pw_pool_give_back_pw(auxbas_pw_pool,v_spin_ddapc_rest_r%pw,&
                                        error=error)
       ENDIF
       IF (calculate_forces.AND.qs_env%dft_control%qs_control%becke_restraint)&
           CALL pw_pool_give_back_pw(auxbas_pw_pool,becke%becke_pot%pw,&
                                        error=error)
       IF(dft_control%apply_efield_field)&
            CALL pw_pool_give_back_pw(auxbas_pw_pool,v_efield_rspace%pw,&
                                        error=error)   
        
       IF(my_kg_gpw) THEN
         CPPostcondition(kg_fm_set%ref_count>=1,cp_failure_level,routineP,error,failure)
         CPPostcondition(kg_sub_pw_env%ref_count>=1,cp_failure_level,routineP,error,failure)
       ! Calculate the KE for each molecule independently
         energy%kg_gpw_ekin_mol = 0.0_dp
         CALL kg_gpw_ekin_mol(qs_env,&
                              kg_sub_pw_env%molbox_env_set,&
                              kg_fm_set%kg_fm_mol_set,&
                              ks_global=ks_matrix,&
                              p_global=rho%rho_ao,ekin_mol=energy%kg_gpw_ekin_mol,&
                              calculate_forces=calculate_forces,&
                              just_energy=just_energy,error=error)
       ! Change sign : it is the correction due to the double counting
         energy%kg_gpw_ekin_mol = -energy%kg_gpw_ekin_mol
       END IF

       IF(gapw) THEN
!      ***  Integrals of the Hartree potential with g0_soft ***
          CALL integrate_vhg0_rspace(qs_env,ks_env%v_hartree_rspace, &
                                     calculate_forces, error=error)

!      ***  Energy, Integrals and Forces coming from the 2-centers terms
!      ***  Inizialize the energy
          energy%hartree_2c = 0.0_dp
          CALL integrate_vhg0_g0(qs_env,energy%hartree_2c,.FALSE.,calculate_forces,error=error)

!      ***  Energy, Integrals and Forces coming from the 3-centers terms
!      *** Initialize energy
          energy%hartree_3c = 0.0_dp
          CALL integrate_a_vhg0_b(qs_env,ks_matrix,rho%rho_ao,&
                       energy%hartree_3c,.FALSE.,calculate_forces,error=error)

       END IF

! **** SCP
!    Single atom contributions to KS matrix due to self-correction of the SCP with 
!    the orbitals. 
       IF (scp_dft) THEN
          CALL integrate_a_vhscp_b ( qs_env, ks_matrix, rho%rho_ao, &
                                     just_energy, error=error )
          CALL scp_qs_energies ( qs_env, calculate_forces, just_energy, error=error )
          
       ENDIF
! **** SCP
       IF(gapw .OR. gapw_xc) THEN
!      ***  Single atom contributions in the KS matrix ***
          CALL update_ks_atom(qs_env,ks_matrix,rho%rho_ao,calculate_forces,error=error)
       ENDIF
!

       !Calculation of Mulliken restraint, if requested
       CALL qs_ks_mulliken_restraint(energy,dft_control,just_energy,para_env,&
            ks_matrix, matrix_s,rho,mulliken_order_p,error)
 
       ! Add DFT+U contribution, if requested
       IF (dft_control%dft_plus_u) THEN
         IF (just_energy) THEN
           CALL plus_u(qs_env=qs_env,error=error)
         ELSE
           CALL plus_u(qs_env=qs_env,matrix_h=ks_matrix,error=error)
         END IF
       ELSE
         energy%dft_plus_u = 0.0_dp
       END IF

       ! prepare to add electronic entropic contribution is needed
       IF (dft_control%smear) THEN
          energy%kTS = 0.0_dp
          energy%efermi = 0.0_dp
          CALL get_qs_env(qs_env,mos=mo_array,error=error)
          DO ispin=1,SIZE(mo_array)
             energy%kTS = energy%kTS + mo_array(ispin)%mo_set%kTS
             energy%efermi = energy%efermi + mo_array(ispin)%mo_set%mu
          ENDDO          
          energy%efermi = energy%efermi /REAL(SIZE(mo_array))
       ENDIF

       ! **** right now we can compute the orbital derivative here, as it depends currently only on the available
       ! **** Kohn-Sham matrix. This might change in the future, in which case more pieces might need to be assembled
       ! **** from this routine, notice that this part of the calculation in not linear scaling
       ! **** right now this operation is only non-trivial because of occupation numbers and the restricted keyword
       IF (qs_env%requires_mo_derivs .AND. .NOT. just_energy.AND..NOT.qs_env%run_rtp) THEN
          CALL get_qs_env(qs_env,mo_derivs=mo_derivs,mos=mo_array,error=error)
          DO ispin=1,SIZE(mo_derivs)
             CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,mo_coeff=mo_coeff, &
                             occupation_numbers=occupation_numbers )
             CALL cp_fm_get_info(mo_coeff,nrow_global=n,ncol_global=k,error=error)
             CALL cp_sm_fm_multiply(ks_matrix(ispin)%matrix,mo_coeff,&
                                          mo_derivs(ispin)%matrix, k,error=error)

             IF( dft_control%do_wfn_fitting) THEN
               CALL get_qs_env(qs_env,mo_derivs=mo_derivs,mos=mos,error=error)
               CALL get_qs_env(qs_env,mos_aux_fit=mos_aux_fit,error=error)
               CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff)
               CALL get_mo_set(mo_set=mos_aux_fit(ispin)%mo_set,mo_coeff=mo_coeff_aux_fit)
               CALL get_qs_env(qs_env,mo_derivs_aux_fit=mo_derivs_aux_fit,error=error)
               CALL merge_mo_derivs(ispin, qs_env%wfn_fitting_env, mo_coeff, mo_coeff_aux_fit, mo_derivs, &
                                    mo_derivs_aux_fit, qs_env%matrix_ks_aux_fit, error)
             END IF

             IF (dft_control%restricted) THEN
                 ! only the first mo_set are actual variables, but we still need both
                 CPPrecondition(ispin==1, cp_failure_level, routineP, error, failure)
                 CPPrecondition(SIZE(mo_array)==2, cp_failure_level, routineP, error, failure)
                 ! use a temporary array with the same size as the first spin for the second spin
                 CALL cp_fm_create(hc_tmp,mo_derivs(1)%matrix%matrix_struct,error=error)
                 CALL cp_fm_set_all(hc_tmp,0.0_dp,error=error)
                 CALL cp_fm_get_info(hc_tmp,ncol_global=k2,error=error)

                 ! uniform_occupation is needed for this case, otherwise we can no reconstruct things in ot, since we 
                 ! irreversibly sum 
                 CALL get_mo_set(mo_set=mo_array(1)%mo_set, uniform_occupation=uniform_occupation)
                 CPPrecondition(uniform_occupation, cp_failure_level, routineP, error, failure)
                 CALL get_mo_set(mo_set=mo_array(2)%mo_set,mo_coeff=mo_coeff, uniform_occupation=uniform_occupation)
                 CPPrecondition(uniform_occupation, cp_failure_level, routineP, error, failure)

                 CALL cp_fm_get_info(mo_coeff,ncol_global=k,error=error)
                 CALL cp_sm_fm_multiply(ks_matrix(2)%matrix,mo_coeff,&
                                        hc_tmp, k,error=error)

                 ! sum the corresponding columns (relies on the final columns of hc_tmp being zero)
                 CALL cp_fm_scale_and_add(1.0_dp, mo_derivs(1)%matrix, 1.0_dp, hc_tmp,error=error)
                 CALL cp_fm_release(hc_tmp,error=error)
             ENDIF

          ENDDO
       ENDIF

       ! deal with low spin roks
        CALL low_spin_roks(energy,qs_env,dft_control,just_energy,&
                                  calculate_forces,auxbas_pw_pool,error)

       ! deal with sic on explicit orbitals
       CALL sic_explicit_orbitals(energy,qs_env,dft_control,poisson_env,just_energy,&
                                  calculate_forces,auxbas_pw_pool,error)

       ! adds s2_restraint energy and orbital derivatives
       CALL qs_ks_s2_restraint(dft_control,qs_env,matrix_s,&
            mo_array,energy,calculate_forces, just_energy,error)

       IF ( do_ppl ) THEN
         ! update core energy for grid based local pseudopotential
         ecore_ppl = 0._dp
         DO ispin=1,nspins
           ecore_ppl = ecore_ppl+SUM(vppl_rspace%pw%cr3d*rho%rho_r(ispin)%pw%cr3d)*vppl_rspace%pw%pw_grid%dvol
         END DO
         CALL mp_sum(ecore_ppl,para_env%group)
         energy%core = energy%core+ecore_ppl
       END IF

       ! sum all energy terms to obtain the total energy
       energy%total = energy%core_overlap + energy%core_self + &
            energy%core + energy%hartree + &
            energy%hartree_1c + energy%hartree_2c + energy%hartree_3c +&
            energy%exc + energy%exc1 + energy%ex + energy%dispersion +&
            energy%qmmm_el + energy%mulliken + &
            energy%kg_gpw_ekin_mol + SUM(energy%ddapc_restraint) + energy%s2_restraint +&
            energy%becke + energy%dft_plus_u + energy%kTS + energy%efield +energy%efield_core

       IF( dft_control%do_wfn_fitting ) THEN
         energy%total = energy%total + energy%exc_aux_fit
       END IF 
! **** SCP
       IF (scp_dft) THEN
         CALL get_scp_env ( scp_env, energy = scp_energy, error = error )
         energy % total = energy % total + scp_energy % e_scp_total
         energy % total = energy % total + scp_energy % e_scp_disp
       ENDIF 
! **** SCP

       output_unit=cp_print_key_unit_nr(logger,input,"DFT%SCF%PRINT%DETAILED_ENERGY",&
            extension=".scfLog",error=error)
       IF (output_unit>0) THEN
          IF(my_kg_gpw) THEN
            WRITE (UNIT=output_unit,FMT="(/,(T3,A,T61,F20.10))")&
               "Core Hamiltonian energy of separated molecules:",energy%core,&
               "Hartree energy of the total density:           ",energy%hartree,&
               "XC and Kin. energy of the total density:       ",energy%exc,&
               "Kin. energy of the molecular densities:        ",energy%kg_gpw_ekin_mol
          ELSE
            WRITE (UNIT=output_unit,FMT="(/,(T3,A,T61,F20.10))")&
               "Core Hamiltonian energy:                       ",energy%core,&
               "Hartree energy:                                ",energy%hartree,&
               "Exchange-correlation energy:                   ",energy%exc
          END IF
          IF (energy%e_hartree /= 0.0_dp)&
               WRITE (UNIT=output_unit,FMT="(T3,A,T61,F20.10)")&
               "Coulomb (electron-electron) energy:            ",energy%e_hartree
          IF (energy%dispersion/= 0.0_dp)&
               WRITE (UNIT=output_unit,FMT="(T3,A,T61,F20.10)")&
               "Dispersion energy:                             ",energy%dispersion
          IF(gapw) THEN
             WRITE (UNIT=output_unit,FMT="(/,(T3,A,T61,F20.10))")&
                "GAPW| Exc from hard and soft atomic rho1:      ",energy%exc1,&
                "GAPW| local Eh = 1 center integrals:           ",energy%hartree_1c,&
                "GAPW| local Eh = 2 center integrals:           ",energy%hartree_2c,&
                "GAPW| local Eh = 3 centers integrals:          ",energy%hartree_3c
          END IF
          IF(gapw_xc) THEN
             WRITE (UNIT=output_unit,FMT="(/,(T3,A,T61,F20.10))")&
                "GAPW| Exc from hard and soft atomic rho1:      ",energy%exc1
          END IF
          IF (dft_control%dft_plus_u) THEN
            WRITE (UNIT=output_unit,FMT="(T3,A,T61,F20.10)")&
              "DFT+U energy:",energy%dft_plus_u
          END IF
          IF  (qs_env%qmmm) THEN
             WRITE (UNIT=output_unit,FMT="(T3,A,T61,F20.10)")&
               "QM/MM Electrostatic energy:                    ",energy%qmmm_el
          END IF
          IF (dft_control%qs_control%mulliken_restraint) THEN
              WRITE (UNIT=output_unit,FMT="(T3,A,T41,2F20.10)")&
                "Mulliken restraint (order_p,energy) : ",mulliken_order_p,energy%mulliken
          ENDIF
          IF (dft_control%qs_control%ddapc_restraint) THEN
             DO n=1,SIZE( dft_control%qs_control%ddapc_restraint_control)
                ddapc_order_p = &
                  dft_control%qs_control%ddapc_restraint_control(n)%ddapc_restraint_control%ddapc_order_p
                WRITE (UNIT=output_unit,FMT="(T3,A,T41,2F20.10)")&
                     "DDAPC restraint (order_p,energy) : ",ddapc_order_p,energy%ddapc_restraint(n)
             END DO
          ENDIF
          IF (dft_control%qs_control%s2_restraint) THEN
             s2_order_p = dft_control%qs_control%s2_restraint_control%s2_order_p
              WRITE (UNIT=output_unit,FMT="(T3,A,T41,2F20.10)")&
                "S2 restraint (order_p,energy) : ",s2_order_p,energy%s2_restraint
          ENDIF

       END IF
       CALL cp_print_key_finished_output(output_unit,logger,input,&
            "DFT%SCF%PRINT%DETAILED_ENERGY", error=error)

    END IF

    CALL timestop(handle)

  END SUBROUTINE qs_ks_build_kohn_sham_matrix

! *****************************************************************************
!> \brief computes the Hartree-Fock energy brute force in a pw basis
!> \note
!>      only computes the HFX energy, no derivatives as yet
!> \par History
!>      12.2007 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE pw_hfx(qs_env,energy,hfx_section,poisson_env,auxbas_pw_pool,irep,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(section_vals_type), POINTER         :: hfx_section
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    INTEGER                                  :: irep
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(*), PARAMETER :: routineN = 'pw_hfx', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blocksize, handle, iloc, &
                                                iorb, iorb_block, ispin, iw, &
                                                jloc, jorb, jorb_block, norb
    LOGICAL                                  :: do_pw_hfx, failure
    REAL(KIND=dp)                            :: exchange_energy, fraction, &
                                                pair_energy, scaling
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: pot_g, rho_g, rho_r
    TYPE(pw_p_type), ALLOCATABLE, &
      DIMENSION(:)                           :: rho_i, rho_j

    CALL timeset(routineN,handle)
    failure = .FALSE.
    logger => cp_error_get_logger(error)

    CALL section_vals_val_get(hfx_section, "PW_HFX", l_val=do_pw_hfx, i_rep_section=irep, error=error)

    IF (do_pw_hfx) THEN
       CALL section_vals_val_get(hfx_section, "FRACTION", r_val=fraction, error=error)
       CALL section_vals_val_get(hfx_section, "PW_HFX_BLOCKSIZE", i_val=blocksize, error=error)
   
       CALL get_qs_env(qs_env,mos=mo_array,pw_env=pw_env, dft_control=dft_control, &
                       cell=cell, particle_set=particle_set, &
                       atomic_kind_set=atomic_kind_set, error=error)

       ! limit the blocksize by the number of orbitals
       CALL get_mo_set(mo_set=mo_array(1)%mo_set,mo_coeff=mo_coeff)
       CALL cp_fm_get_info(mo_coeff,ncol_global=norb,error=error)
       blocksize=MAX(1,MIN(blocksize,norb))

       CALL pw_pool_create_pw(auxbas_pw_pool,rho_r%pw,&
                               use_data=REALDATA3D,&
                               in_space=REALSPACE,error=error)
       CALL pw_pool_create_pw(auxbas_pw_pool,rho_g%pw,&
                               use_data=COMPLEXDATA1D,&
                               in_space=RECIPROCALSPACE,error=error)
       CALL pw_pool_create_pw(auxbas_pw_pool,pot_g%pw,&
                               use_data=COMPLEXDATA1D,&
                               in_space=RECIPROCALSPACE,error=error)

       ALLOCATE(rho_i(blocksize)) 
       ALLOCATE(rho_j(blocksize)) 
  
       DO iorb_block=1,blocksize 
          NULLIFY(rho_i(iorb_block)%pw)
          CALL pw_create(rho_i(iorb_block)%pw,rho_r%pw%pw_grid,&
                                  use_data=REALDATA3D,&
                                  in_space=REALSPACE,error=error)
          NULLIFY(rho_j(iorb_block)%pw)
          CALL pw_create(rho_j(iorb_block)%pw,rho_r%pw%pw_grid,&
                                  use_data=REALDATA3D,&
                                  in_space=REALSPACE,error=error)
       ENDDO

       exchange_energy = 0.0_dp
   
       DO ispin=1,SIZE(mo_array)
          CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,mo_coeff=mo_coeff)
          CALL cp_fm_get_info(mo_coeff,ncol_global=norb,error=error)

          DO iorb_block=1,norb,blocksize

             DO iorb=iorb_block,MIN(iorb_block+blocksize-1,norb)

                 iloc=iorb-iorb_block+1
                 CALL calculate_wavefunction(mo_coeff,iorb,rho_i(iloc),rho_g, &
                          atomic_kind_set,cell,dft_control,particle_set, &
                          pw_env,error=error)

             ENDDO
   
             DO jorb_block=iorb_block,norb,blocksize
  
                DO jorb=jorb_block,MIN(jorb_block+blocksize-1,norb)

                   jloc=jorb-jorb_block+1
                   CALL calculate_wavefunction(mo_coeff,jorb,rho_j(jloc),rho_g, &
                            atomic_kind_set,cell,dft_control,particle_set, &
                            pw_env,error=error)

                ENDDO
   
                DO iorb=iorb_block,MIN(iorb_block+blocksize-1,norb)
                   iloc=iorb-iorb_block+1
                   DO jorb=jorb_block,MIN(jorb_block+blocksize-1,norb)
                      jloc=jorb-jorb_block+1 
                      IF (jorb<iorb) CYCLE

                      ! compute the pair density
                      rho_r%pw%cr3d = rho_i(iloc)%pw%cr3d * rho_j(jloc)%pw%cr3d
         
                      ! go the g-space and compute hartree energy
                      CALL pw_transfer(rho_r%pw, rho_g%pw, error=error)
                      CALL pw_poisson_solve(poisson_env,rho_g%pw, pair_energy, pot_g%pw,error=error)
         
                      ! sum up to the full energy
                      scaling=fraction
                      IF (SIZE(mo_array)==1) scaling=scaling*2.0_dp
                      IF (iorb/=jorb) scaling=scaling*2.0_dp
         
                      exchange_energy=exchange_energy - scaling * pair_energy

                   ENDDO
                ENDDO
        
             ENDDO
          ENDDO
       ENDDO
       
       DO iorb_block=1,blocksize 
          CALL pw_release(rho_i(iorb_block)%pw,error=error)
          CALL pw_release(rho_j(iorb_block)%pw,error=error)
       ENDDO

       CALL pw_pool_give_back_pw(auxbas_pw_pool,rho_r%pw,error=error)
       CALL pw_pool_give_back_pw(auxbas_pw_pool,rho_g%pw,error=error)
       CALL pw_pool_give_back_pw(auxbas_pw_pool,pot_g%pw,error=error)
   
       iw = cp_print_key_unit_nr(logger,hfx_section,"HF_INFO",&
              extension=".scfLog",error=error)
   
       IF (iw>0) THEN
           WRITE (UNIT=iw,FMT="((T3,A,T61,F20.10))")&
                    "HF_PW_HFX| PW exchange energy:",exchange_energy
           WRITE (UNIT=iw,FMT="((T3,A,T61,F20.10),/)")&
                    "HF_PW_HFX| Gaussian exchange energy:",energy%ex
       ENDIF
   
       CALL cp_print_key_finished_output(iw,logger,hfx_section,&
                                         "HF_INFO", error=error)
   
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE pw_hfx

! *****************************************************************************
!> \brief do ROKS calculations yielding low spin states
! *****************************************************************************
  SUBROUTINE low_spin_roks(energy,qs_env,dft_control,just_energy,&
                               calculate_forces,auxbas_pw_pool,error)

    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(dft_control_type), POINTER          :: dft_control
    LOGICAL, INTENT(IN)                      :: just_energy, calculate_forces
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(*), PARAMETER :: routineN = 'low_spin_roks', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, iterm, k, &
                                                k_alpha, k_beta, n_rep, &
                                                Nelectron, Nspin, Nterms, stat
    INTEGER, DIMENSION(:), POINTER           :: ivec
    INTEGER, DIMENSION(:, :, :), POINTER     :: occupations
    LOGICAL                                  :: failure, in_range, &
                                                uniform_occupation
    REAL(KIND=dp)                            :: exc, total_rho
    REAL(KIND=dp), DIMENSION(:), POINTER     :: energy_scaling, rvec, scaling
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mo_derivs
    TYPE(cp_fm_type), POINTER                :: fm_deriv, fm_scaled, mo_coeff
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: work_v_rspace
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_g, rho_r, tau, vxc, &
                                                vxc_tau
    TYPE(pw_pool_type), POINTER              :: xc_pw_pool
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_h, matrix_p
    TYPE(section_vals_type), POINTER         :: input, low_spin_roks_section, &
                                                xc_section

    IF (.NOT. dft_control%low_spin_roks) RETURN
    failure=.FALSE.

    CALL timeset(routineN,handle)

    CALL get_qs_env(qs_env,mo_derivs=mo_derivs,mos=mo_array,rho=rho,pw_env=pw_env, &
                    input=input, cell=cell, error=error)
    xc_section => section_vals_get_subs_vals(input,"DFT%XC",error=error)

    ! some assumptions need to be checked
    ! we have two spins
    CPPrecondition(SIZE(mo_array,1)==2,cp_failure_level,routineP,error,failure)
    Nspin=2
    ! we want uniform occupations
    CALL get_mo_set(mo_set=mo_array(1)%mo_set, uniform_occupation=uniform_occupation)
    CPPrecondition(uniform_occupation, cp_failure_level, routineP, error, failure)
    CALL get_mo_set(mo_set=mo_array(2)%mo_set,mo_coeff=mo_coeff, uniform_occupation=uniform_occupation)
    CPPrecondition(uniform_occupation, cp_failure_level, routineP, error, failure)

    ! basic info
    CALL get_mo_set(mo_set=mo_array(1)%mo_set,mo_coeff=mo_coeff)
    CALL cp_fm_get_info(mo_coeff,ncol_global=k_alpha,error=error)
    CALL get_mo_set(mo_set=mo_array(2)%mo_set,mo_coeff=mo_coeff)
    CALL cp_fm_get_info(mo_coeff,ncol_global=k_beta,error=error)

    ! read the input
    low_spin_roks_section => section_vals_get_subs_vals(input,"DFT%LOW_SPIN_ROKS", error=error)

    CALL section_vals_val_get(low_spin_roks_section,"ENERGY_SCALING",r_vals=rvec,error=error)
    Nterms=SIZE(rvec)
    ALLOCATE(energy_scaling(Nterms))
    energy_scaling=rvec !? just wondering, should this add up to 1, in which case we should cpp?

    CALL section_vals_val_get(low_spin_roks_section,"SPIN_CONFIGURATION",n_rep_val=n_rep,error=error)
    CPPostcondition(n_rep==Nterms, cp_failure_level, routineP, error, failure)
    CALL section_vals_val_get(low_spin_roks_section,"SPIN_CONFIGURATION",i_rep_val=1,i_vals=ivec,error=error)
    Nelectron=SIZE(ivec)
    CPPostcondition(Nelectron==k_alpha-k_beta, cp_failure_level, routineP, error, failure)
    ALLOCATE(occupations(2,Nelectron,Nterms))
    occupations=0
    DO iterm=1,Nterms
       CALL section_vals_val_get(low_spin_roks_section,"SPIN_CONFIGURATION",i_rep_val=iterm,i_vals=ivec,error=error)
       CPPostcondition(Nelectron==SIZE(ivec), cp_failure_level, routineP, error, failure)
       in_range=ALL(ivec>=1) .AND. ALL(ivec<=2)
       CPPostcondition(in_range, cp_failure_level, routineP, error, failure)
       DO k=1,Nelectron
          occupations(ivec(k),k,iterm)=1 
       ENDDO
    ENDDO

    ! set up general data structures
    ! density matrices, kohn-sham matrices
    ALLOCATE(matrix_p(Nspin),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO ispin=1,Nspin
       NULLIFY(matrix_p(ispin)%matrix)
       CALL replicate_matrix_structure(source=rho%rho_ao(1)%matrix,&
                                       TARGET=matrix_p(ispin)%matrix,&
                                       target_name="density matrix low spin roks",error=error)
    ENDDO

    ALLOCATE(matrix_h(Nspin),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO ispin=1,Nspin
       NULLIFY(matrix_h(ispin)%matrix)
       CALL replicate_matrix_structure(source=rho%rho_ao(1)%matrix,&
                                       TARGET=matrix_h(ispin)%matrix,&
                                       target_name="KS matrix low spin roks",error=error)
    ENDDO

    ! grids in real and g space for rho and vxc
    ! tau functionals are not supported
    NULLIFY(tau,vxc_tau,vxc)
    CALL pw_env_get(pw_env,xc_pw_pool=xc_pw_pool, error=error)

    ALLOCATE(rho_r(Nspin),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(rho_g(Nspin),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO ispin=1,Nspin
       CALL pw_pool_create_pw(auxbas_pw_pool,rho_r(ispin)%pw,&
                            use_data=REALDATA3D,&
                            in_space=REALSPACE,error=error)
       CALL pw_pool_create_pw(auxbas_pw_pool,rho_g(ispin)%pw,&
                            use_data=COMPLEXDATA1D,&
                            in_space=RECIPROCALSPACE,error=error)
    ENDDO
    CALL pw_pool_create_pw(auxbas_pw_pool,work_v_rspace%pw,&
                            use_data=REALDATA3D,&
                            in_space=REALSPACE,error=error)

    ! get mo matrices needed to construct the density matrices
    ! we will base all on the alpha spin matrix, obviously possible in ROKS
    CALL get_mo_set(mo_set=mo_array(1)%mo_set,mo_coeff=mo_coeff)
    CALL cp_fm_create(fm_scaled,mo_coeff%matrix_struct,error=error)
    CALL cp_fm_create(fm_deriv,mo_coeff%matrix_struct,error=error)
    ALLOCATE(scaling(k_alpha),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ! for each term, add it with the given scaling factor to the energy, and compute the required derivatives
    DO iterm=1,Nterms

       DO ispin=1,Nspin
          ! compute the proper density matrices with the required occupations
          CALL set_matrix(matrix_p(ispin)%matrix,0.0_dp)
          scaling=1.0_dp
          scaling(k_alpha-Nelectron+1:k_alpha)=occupations(ispin,:,iterm)
          CALL cp_fm_to_fm(mo_coeff,fm_scaled,error=error)
          CALL cp_fm_column_scale(fm_scaled,scaling)
          CALL cp_sm_plus_fm_fm_t(sparse_matrix=matrix_p(ispin)%matrix,&
                               matrix_v=mo_coeff, matrix_g=fm_scaled,&
                               ncol=k_alpha,error=error)
          ! compute the densities on the grid
          CALL calculate_rho_elec(matrix_p=matrix_p(ispin)%matrix,&
                rho=rho_r(ispin),rho_gspace=rho_g(ispin), total_rho=total_rho,&
                qs_env=qs_env, error=error)
       ENDDO

       ! compute the exchange energies / potential if needed
       IF (just_energy) THEN
           exc=xc_exc_calc(rho_r=rho_r, rho_g=rho_g, tau=tau, xc_section=xc_section,&
                   cell=cell, pw_pool=xc_pw_pool, error=error)
       ELSE
           CALL xc_vxc_pw_create1(vxc_rho=vxc, rho_r=rho_r,&
                   rho_g=rho_g, tau=tau, vxc_tau=vxc_tau, exc=exc, xc_section=xc_section, cell=cell, &
                   pw_pool=xc_pw_pool, error=error)
       END IF

       energy%exc = energy%exc + energy_scaling(iterm) * exc

       ! add the corresponding derivatives to the MO derivatives
       IF (.NOT. just_energy) THEN
           ! get the potential in matrix form
           DO ispin=1,Nspin
              ! use a work_v_rspace
              work_v_rspace%pw%cr3d = (energy_scaling(iterm) * vxc(ispin) % pw %pw_grid%dvol)* &
                                      vxc(ispin) % pw % cr3d 
              ! zero first ?!
              CALL set_matrix(matrix_h(ispin)%matrix,0.0_dp)
              CALL integrate_v_rspace(v_rspace=work_v_rspace,p=matrix_p(ispin),h=matrix_h(ispin),&
                                      qs_env=qs_env,calculate_forces=calculate_forces,error=error) 
              CALL pw_pool_give_back_pw(auxbas_pw_pool,vxc(ispin)%pw,error=error)
           ENDDO
           DEALLOCATE(vxc)

           ! add this to the mo_derivs, again based on the alpha mo_coeff
           DO ispin=1,Nspin
              CALL cp_sm_fm_multiply(matrix_h(ispin)%matrix,mo_coeff,fm_deriv, k_alpha, error=error)
              scaling=1.0_dp
              scaling(k_alpha-Nelectron+1:k_alpha)=occupations(ispin,:,iterm)
              CALL cp_fm_column_scale(fm_deriv,scaling)
              CALL cp_fm_scale_and_add(1.0_dp,mo_derivs(1)%matrix, &
                                       1.0_dp,fm_deriv,error=error)
           ENDDO

       ENDIF

    ENDDO

    ! release allocated memory
    DO ispin=1,Nspin
       CALL pw_pool_give_back_pw(auxbas_pw_pool,rho_r(ispin)%pw,error=error)
       CALL pw_pool_give_back_pw(auxbas_pw_pool,rho_g(ispin)%pw,error=error)
       CALL deallocate_matrix(matrix_p(ispin)%matrix,error=error)
       CALL deallocate_matrix(matrix_h(ispin)%matrix,error=error)
    ENDDO
    DEALLOCATE(rho_r,rho_g,matrix_p,matrix_h)
    CALL pw_pool_give_back_pw(auxbas_pw_pool,work_v_rspace%pw,error=error)

    CALL cp_fm_release(fm_deriv,error=error)
    CALL cp_fm_release(fm_scaled,error=error)

    DEALLOCATE(occupations)
    DEALLOCATE(energy_scaling)
    DEALLOCATE(scaling)

    CALL timestop(handle)

  END SUBROUTINE low_spin_roks
! *****************************************************************************
!> \brief do sic calculations on explicit orbitals
! *****************************************************************************
  SUBROUTINE sic_explicit_orbitals(energy,qs_env,dft_control,poisson_env,just_energy,&
                               calculate_forces,auxbas_pw_pool,error)

    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    LOGICAL, INTENT(IN)                      :: just_energy, calculate_forces
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(*), PARAMETER :: routineN = 'sic_explicit_orbitals', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, Iorb, k_alpha, &
                                                k_beta, Norb, stat
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: sic_orbital_list
    LOGICAL                                  :: failure, uniform_occupation
    REAL(KIND=dp)                            :: ener, exc, total_rho
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mo_derivs, mo_derivs_local
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: matrix_hv, matrix_v, mo_coeff
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: orb_rho_g, orb_rho_r, tmp_g, &
                                                tmp_r, work_v_gspace, &
                                                work_v_rspace
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_g, rho_r, tau, vxc, &
                                                vxc_tau
    TYPE(pw_pool_type), POINTER              :: xc_pw_pool
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(real_matrix_p_type)                 :: orb_density_matrix_p, orb_h_p
    TYPE(real_matrix_type), POINTER          :: orb_density_matrix, orb_h
    TYPE(section_vals_type), POINTER         :: input, xc_section

    IF (dft_control%sic_method_id .NE. sic_eo) RETURN

    CALL timeset(routineN,handle)

    failure=.FALSE.

    NULLIFY(tau,vxc_tau)

    ! generate the lists of orbitals that need sic treatment
    CALL get_qs_env(qs_env,mo_derivs=mo_derivs,mos=mo_array,rho=rho,pw_env=pw_env, &
                    input=input, cell=cell, error=error)
    xc_section => section_vals_get_subs_vals(input,"DFT%XC",error=error)

    CALL pw_env_get(pw_env,xc_pw_pool=xc_pw_pool, error=error)

    ! we have two spins
    CPPrecondition(SIZE(mo_array,1)==2,cp_failure_level,routineP,error,failure)
    ! we want uniform occupations
    CALL get_mo_set(mo_set=mo_array(1)%mo_set, uniform_occupation=uniform_occupation)
    CPPrecondition(uniform_occupation, cp_failure_level, routineP, error, failure)
    CALL get_mo_set(mo_set=mo_array(2)%mo_set,mo_coeff=mo_coeff, uniform_occupation=uniform_occupation)
    CPPrecondition(uniform_occupation, cp_failure_level, routineP, error, failure)

    k_alpha=0 ; k_beta=0
    SELECT CASE(dft_control%sic_list_id)
    CASE(sic_list_all)

      CALL get_mo_set(mo_set=mo_array(1)%mo_set,mo_coeff=mo_coeff)
      CALL cp_fm_get_info(mo_coeff,ncol_global=k_alpha,error=error)

      IF (SIZE(mo_array,1)>1) THEN
          CALL get_mo_set(mo_set=mo_array(2)%mo_set,mo_coeff=mo_coeff)
          CALL cp_fm_get_info(mo_coeff,ncol_global=k_beta,error=error)
      ENDIF

      Norb=k_alpha + k_beta
      ALLOCATE(sic_orbital_list(3,Norb),stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      iorb=0
      DO i=1,k_alpha
         iorb=iorb+1 
         sic_orbital_list(1,iorb)=1
         sic_orbital_list(2,iorb)=i
         sic_orbital_list(3,iorb)=1
      ENDDO
      DO i=1,k_beta
         iorb=iorb+1 
         sic_orbital_list(1,iorb)=2
         sic_orbital_list(2,iorb)=i
         IF (SIZE(mo_derivs,1)==1) THEN
             sic_orbital_list(3,iorb)=1
         ELSE
             sic_orbital_list(3,iorb)=2
         ENDIF
      ENDDO

    CASE(sic_list_unpaired)
      ! we have two spins
      CPPrecondition(SIZE(mo_array,1)==2,cp_failure_level,routineP,error,failure)
      ! we have them restricted 
      CPPrecondition(SIZE(mo_derivs,1)==1,cp_failure_level,routineP,error,failure)
      CPPrecondition(dft_control%restricted,cp_failure_level,routineP,error,failure)

      CALL get_mo_set(mo_set=mo_array(1)%mo_set,mo_coeff=mo_coeff)
      CALL cp_fm_get_info(mo_coeff,ncol_global=k_alpha,error=error)

      CALL get_mo_set(mo_set=mo_array(2)%mo_set,mo_coeff=mo_coeff)
      CALL cp_fm_get_info(mo_coeff,ncol_global=k_beta,error=error)

      Norb=k_alpha-k_beta
      ALLOCATE(sic_orbital_list(3,Norb),stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      
      iorb=0
      DO i=k_beta+1,k_alpha
         iorb=iorb+1 
         sic_orbital_list(1,iorb)=1
         sic_orbital_list(2,iorb)=i
         ! we are guaranteed to be restricted
         sic_orbital_list(3,iorb)=1
      ENDDO

    CASE DEFAULT
      CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT

    ! data needed for each of the orbs
    CALL pw_pool_create_pw(auxbas_pw_pool,orb_rho_r%pw,&
                            use_data=REALDATA3D,&
                            in_space=REALSPACE,error=error)
    CALL pw_pool_create_pw(auxbas_pw_pool,tmp_r%pw,&
                            use_data=REALDATA3D,&
                            in_space=REALSPACE,error=error)
    CALL pw_pool_create_pw(auxbas_pw_pool,orb_rho_g%pw,&
                            use_data=COMPLEXDATA1D,&
                            in_space=RECIPROCALSPACE,error=error)
    CALL pw_pool_create_pw(auxbas_pw_pool,tmp_g%pw,&
                            use_data=COMPLEXDATA1D,&
                            in_space=RECIPROCALSPACE,error=error)
    CALL pw_pool_create_pw(auxbas_pw_pool,work_v_gspace%pw,&
                            use_data = COMPLEXDATA1D,&
                            in_space = RECIPROCALSPACE, error=error)
    CALL pw_pool_create_pw(auxbas_pw_pool,work_v_rspace%pw,&
                            use_data=REALDATA3D,&
                            in_space=REALSPACE,error=error)

    NULLIFY(orb_density_matrix)
    CALL replicate_matrix_structure(source=rho%rho_ao(1)%matrix,&
                                    TARGET=orb_density_matrix,&
                                    target_name="orb_density_matrix",error=error)
    orb_density_matrix_p%matrix=>orb_density_matrix

    NULLIFY(orb_h)
    CALL replicate_matrix_structure(source=rho%rho_ao(1)%matrix,&
                                    TARGET=orb_h,&
                                    target_name="orb_density_matrix",error=error)
    orb_h_p%matrix=>orb_h

    CALL get_mo_set(mo_set=mo_array(1)%mo_set,mo_coeff=mo_coeff)

    CALL cp_fm_struct_create(fm_struct_tmp, ncol_global=1, &
                             template_fmstruct=mo_coeff%matrix_struct, error=error)
    CALL cp_fm_create(matrix_v,fm_struct_tmp, name="matrix_v",error=error)
    CALL cp_fm_create(matrix_hv,fm_struct_tmp, name="matrix_hv",error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)

    ALLOCATE(mo_derivs_local(SIZE(mo_derivs,1)),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO I=1,SIZE(mo_derivs,1)
       CALL cp_fm_create(mo_derivs_local(I)%matrix,mo_derivs(I)%matrix%matrix_struct,error=error)
    ENDDO

    ALLOCATE(rho_r(2),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    rho_r(1)%pw=>orb_rho_r%pw
    rho_r(2)%pw=>tmp_r%pw
    CALL pw_zero(tmp_r%pw, error=error)
    
    ALLOCATE(rho_g(2),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    rho_g(1)%pw=>orb_rho_g%pw
    rho_g(2)%pw=>tmp_g%pw
    CALL pw_zero(tmp_g%pw, error=error)

    NULLIFY(vxc)
    ! ALLOCATE(vxc(2),stat=stat)
    ! CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ! CALL pw_pool_create_pw(xc_pw_pool,vxc(1)%pw,&
    !         in_space=REALSPACE, use_data=REALDATA3D,error=error)
    ! CALL pw_pool_create_pw(xc_pw_pool,vxc(2)%pw,&
    !         in_space=REALSPACE, use_data=REALDATA3D,error=error)

    ! now apply to SIC correction to each selected orbital
    DO iorb=1,Norb
       ! extract the proper orbital from the mo_coeff
       CALL get_mo_set(mo_set=mo_array(sic_orbital_list(1,iorb))%mo_set,mo_coeff=mo_coeff)
       CALL cp_fm_to_fm(mo_coeff,matrix_v,1,sic_orbital_list(2,iorb),1)

       ! construct the density matrix and the corresponding density
       CALL set_matrix(orb_density_matrix,0.0_dp)
       CALL cp_sm_plus_fm_fm_t(orb_density_matrix,matrix_v=matrix_v,ncol=1,&
                                  alpha=1.0_dp,error=error)

       CALL calculate_rho_elec(matrix_p=orb_density_matrix,&
                rho=orb_rho_r,rho_gspace=orb_rho_g, total_rho=total_rho,&
                qs_env=qs_env, error=error)

       ! write(6,*) 'Orbital ',sic_orbital_list(1,iorb),sic_orbital_list(2,iorb)
       ! write(6,*) 'Total orbital rho= ',total_rho

       ! compute the energy functional for this orbital and its derivative

       CALL pw_poisson_solve(poisson_env,orb_rho_g%pw, ener, work_v_gspace%pw,error=error)
       energy%hartree=energy%hartree - dft_control % sic_scaling_a * ener
       IF (.NOT. just_energy) THEN
            CALL pw_transfer(work_v_gspace%pw, work_v_rspace%pw, error=error)
            CALL pw_scale(work_v_rspace%pw, - dft_control % sic_scaling_a * work_v_rspace%pw%pw_grid%dvol,&
                 error=error)
            CALL set_matrix(orb_h,0.0_dp)
       ENDIF

       IF (just_energy) THEN
           exc=xc_exc_calc(rho_r=rho_r, rho_g=rho_g, tau=tau, xc_section=xc_section,&
                   cell=cell, pw_pool=xc_pw_pool, error=error)
       ELSE
           CALL xc_vxc_pw_create1(vxc_rho=vxc, rho_r=rho_r,&
                   rho_g=rho_g, tau=tau, vxc_tau=vxc_tau, exc=exc, xc_section=xc_section, cell=cell, &
                   pw_pool=xc_pw_pool, error=error)
           ! add to the existing work_v_rspace
           work_v_rspace%pw%cr3d = work_v_rspace%pw%cr3d - &
                   dft_control % sic_scaling_b * vxc(1) % pw %pw_grid%dvol *  vxc(1) % pw % cr3d 
       END IF
       energy%exc = energy%exc - dft_control % sic_scaling_b * exc

       IF (.NOT. just_energy) THEN
           ! note, orb_h (which is being pointed to with orb_h_p) is zeroed above
           CALL integrate_v_rspace(v_rspace=work_v_rspace,p=orb_density_matrix_p,h=orb_h_p,&
                                   qs_env=qs_env,calculate_forces=calculate_forces,error=error) 

           ! add this to the mo_derivs
           CALL cp_sm_fm_multiply(orb_h,matrix_v,matrix_hv, 1, error=error)
           ! silly trick, copy to an array of the right size and add to mo_derivs
           CALL cp_fm_set_all(mo_derivs_local(sic_orbital_list(3,iorb))%matrix,0.0_dp,error=error)
           CALL cp_fm_to_fm(matrix_hv,mo_derivs_local(sic_orbital_list(3,iorb))%matrix,1,1,sic_orbital_list(2,iorb))
           CALL cp_fm_scale_and_add(1.0_dp,mo_derivs(sic_orbital_list(3,iorb))%matrix, &
                                    1.0_dp,mo_derivs_local(sic_orbital_list(3,iorb))%matrix,error=error)

           ! need to deallocate vxc
           CALL pw_pool_give_back_pw(xc_pw_pool,vxc(1)%pw,error=error)
           CALL pw_pool_give_back_pw(xc_pw_pool,vxc(2)%pw,error=error)
           DEALLOCATE(vxc)

       ENDIF

    ENDDO

    CALL pw_pool_give_back_pw(auxbas_pw_pool,orb_rho_r%pw,error=error)
    CALL pw_pool_give_back_pw(auxbas_pw_pool,tmp_r%pw,error=error)
    CALL pw_pool_give_back_pw(auxbas_pw_pool,orb_rho_g%pw,error=error)
    CALL pw_pool_give_back_pw(auxbas_pw_pool,tmp_g%pw,error=error)
    CALL pw_pool_give_back_pw(auxbas_pw_pool,work_v_gspace%pw,error=error)
    CALL pw_pool_give_back_pw(auxbas_pw_pool,work_v_rspace%pw,error=error)

    CALL deallocate_matrix(orb_density_matrix,error=error)
    CALL deallocate_matrix(orb_h,error=error)
    CALL cp_fm_release(matrix_v,error)
    CALL cp_fm_release(matrix_hv,error)
    DO I=1,SIZE(mo_derivs,1)
       CALL cp_fm_release(mo_derivs_local(I)%matrix,error=error)
    ENDDO
    DEALLOCATE(mo_derivs_local)
    DEALLOCATE(rho_r)
    DEALLOCATE(rho_g)

    CALL timestop(handle)

  END SUBROUTINE sic_explicit_orbitals

! *****************************************************************************
!> \brief do sic calculations on the spin density
! *****************************************************************************
  SUBROUTINE calc_v_sic_rspace(v_sic_rspace,energy,&
                               qs_env,dft_control,rho,poisson_env,just_energy,&
                               calculate_forces,auxbas_pw_pool,error)

    TYPE(pw_p_type), INTENT(INOUT)           :: v_sic_rspace
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    LOGICAL, INTENT(IN)                      :: just_energy, calculate_forces
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(*), PARAMETER :: routineN = 'calc_v_sic_rspace', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, nelec, nelec_a, nelec_b, &
                                                nforce
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: ener, full_scaling, scaling
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: store_forces
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(pw_p_type)                          :: work_rho, work_v
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force

    failure = .FALSE.
    NULLIFY(mo_array)

    IF (dft_control%sic_method_id == sic_none) RETURN
    IF (dft_control%sic_method_id == sic_eo) RETURN

    CALL cp_assert(.NOT. dft_control%qs_control%gapw, cp_failure_level,cp_assertion_failed,routineP,&
                   "sic and GAPW not yet compatible",error,failure)

    ! OK, right now we like two spins to do sic, could be relaxed for AD
    CPPrecondition(dft_control%nspins == 2,cp_failure_level,routineP,error,failure)

    CALL pw_pool_create_pw(auxbas_pw_pool, work_rho%pw, &
                            use_data = COMPLEXDATA1D,&
                            in_space = RECIPROCALSPACE, error=error)
    CALL pw_pool_create_pw(auxbas_pw_pool, work_v%pw,&
                            use_data = COMPLEXDATA1D,&
                            in_space = RECIPROCALSPACE, error=error)

    ! Hartree sic corrections
    SELECT CASE ( dft_control%sic_method_id )
    CASE ( sic_mauri_us, sic_mauri_spz )
       CALL pw_copy(rho%rho_g(1)%pw,work_rho%pw, error=error)
       CALL pw_axpy(rho%rho_g(2)%pw,work_rho%pw,alpha=-1._dp, error=error)
       CALL pw_poisson_solve(poisson_env,work_rho%pw, ener, work_v%pw,error=error)
    CASE ( sic_ad )
       ! find out how many elecs we have
       CALL get_qs_env(qs_env,mos=mo_array,error=error)
       CALL get_mo_set(mo_set=mo_array(1)%mo_set,nelectron=nelec_a)
       CALL get_mo_set(mo_set=mo_array(2)%mo_set,nelectron=nelec_b)
       nelec = nelec_a + nelec_b
       CALL pw_copy(rho%rho_g(1)%pw,work_rho%pw, error=error)
       CALL pw_axpy(rho%rho_g(2)%pw,work_rho%pw, error=error)
       scaling = 1.0_dp / REAL(nelec,KIND=dp)
       CALL pw_scale(work_rho%pw,scaling, error=error)
       CALL pw_poisson_solve(poisson_env,work_rho%pw, ener, work_v%pw,error=error)
    CASE DEFAULT
       CALL cp_assert( .FALSE., cp_failure_level,cp_assertion_failed,routineP,&
            "Unknown sic method id",error,failure)
    END SELECT

    ! Correct for  DDAP charges (if any)
    ! storing whatever force might be there from previous decoupling
    IF (calculate_forces) THEN
       CALL get_qs_env(qs_env=qs_env,force=force,error=error)
       nforce=0
       DO i=1,SIZE(force)
          nforce=nforce+SIZE(force(i)%ch_pulay,2)
       ENDDO
       ALLOCATE(store_forces(3,nforce))
       nforce=0
       DO i=1,SIZE(force)
          store_forces(1:3,nforce+1:nforce+SIZE(force(i)%ch_pulay,2))=force(i)%ch_pulay(:,:)
          force(i)%ch_pulay(:,:)=0.0_dp
          nforce=nforce+SIZE(force(i)%ch_pulay,2)
       ENDDO
    ENDIF

    CALL cp_ddapc_apply_CD(qs_env,&
                           work_rho,&
                           ener,&
                           v_hartree_gspace=work_v,&
                           calculate_forces=calculate_forces,&
                           Itype_of_density="SPIN",&
                           error=error)

    SELECT CASE ( dft_control%sic_method_id )
    CASE ( sic_mauri_us, sic_mauri_spz )
       full_scaling= - dft_control % sic_scaling_a 
    CASE ( sic_ad )
       full_scaling= - dft_control % sic_scaling_a * nelec
    CASE DEFAULT
       CALL cp_assert( .FALSE., cp_failure_level,cp_assertion_failed,routineP,&
            "Unknown sic method id",error,failure)
    END SELECT
    energy%hartree=energy%hartree + full_scaling * ener

    ! add scaled forces, restoring the old
    IF (calculate_forces) THEN
       nforce=0
       DO i=1,SIZE(force)
          force(i)%ch_pulay(:,:)=force(i)%ch_pulay(:,:)*full_scaling + store_forces(1:3,nforce+1:nforce+SIZE(force(i)%ch_pulay,2))
          nforce=nforce+SIZE(force(i)%ch_pulay,2)
       ENDDO
    ENDIF

    IF (.NOT. just_energy) THEN
       CALL pw_pool_create_pw(auxbas_pw_pool,v_sic_rspace%pw,&
                               use_data=REALDATA3D, in_space=REALSPACE,error=error)
       CALL pw_transfer(work_v%pw, v_sic_rspace%pw, error=error)
       ! also take into account the scaling (in addition to the volume element)
       CALL pw_scale(v_sic_rspace%pw, &
            dft_control % sic_scaling_a * v_sic_rspace%pw%pw_grid%dvol, error=error )
    ENDIF

    CALL pw_pool_give_back_pw(auxbas_pw_pool,work_rho%pw,error=error)
    CALL pw_pool_give_back_pw(auxbas_pw_pool,work_v%pw,error=error)

  END SUBROUTINE calc_v_sic_rspace

! *****************************************************************************
  SUBROUTINE calc_rho_tot_gspace(rho_tot_gspace, qs_env, rho, skip_nuclear_density, error)
    TYPE(pw_p_type), INTENT(INOUT)           :: rho_tot_gspace
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_rho_type), POINTER               :: rho
    LOGICAL, INTENT(IN), OPTIONAL            :: skip_nuclear_density
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(*), PARAMETER :: routineN = 'calc_rho_tot_gspace', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin
    LOGICAL                                  :: failure, my_skip
    TYPE(pw_p_type), POINTER                 :: rho0_s_gs, rho_core
    TYPE(qs_charges_type), POINTER           :: qs_charges

    NULLIFY (rho_core, rho0_s_gs)
    my_skip = .FALSE.
    IF (PRESENT(skip_nuclear_density)) my_skip=skip_nuclear_density
    CALL get_qs_env(qs_env=qs_env,&
                    rho_core=rho_core,&
                    rho0_s_gs=rho0_s_gs,&
                    qs_charges=qs_charges,error=error)

    IF (.NOT.my_skip) THEN
       IF(qs_env%dft_control%qs_control%gapw ) THEN
           IF( qs_env%dft_control%qs_control%gapw_control%nopaw_as_gpw) THEN
              CPPrecondition(ASSOCIATED(rho0_s_gs),cp_failure_level,routineP,error,failure)
              CALL pw_copy(rho0_s_gs%pw,rho_tot_gspace%pw, error=error)
              CALL pw_axpy(rho_core%pw,rho_tot_gspace%pw, error=error)
              CALL pw_axpy(rho%rho_g(1)%pw,rho_tot_gspace%pw, error=error)
           ELSE
              CPPrecondition(ASSOCIATED(rho0_s_gs),cp_failure_level,routineP,error,failure)
              CALL pw_copy(rho0_s_gs%pw,rho_tot_gspace%pw, error=error)
              CALL pw_axpy(rho%rho_g(1)%pw,rho_tot_gspace%pw, error=error)
           END IF
       ELSE
          CALL pw_copy(rho_core%pw,rho_tot_gspace%pw, error=error)
          CALL pw_axpy(rho%rho_g(1)%pw,rho_tot_gspace%pw, error=error)
       END IF
    ELSE
       CALL pw_axpy(rho%rho_g(1)%pw, rho_tot_gspace%pw, error=error)
    END IF
    DO ispin=2, qs_env%dft_control%nspins
       CALL pw_axpy(rho%rho_g(ispin)%pw, rho_tot_gspace%pw, error=error)
    END DO
    IF (.NOT.my_skip) &
         qs_charges%total_rho_gspace = pw_integrate_function(rho_tot_gspace%pw,isign=-1, error=error)

  END SUBROUTINE calc_rho_tot_gspace

! *****************************************************************************
  SUBROUTINE print_densities(qs_env, rho_tot_gspace, rho, rho_xc, output_unit, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(pw_p_type), INTENT(IN)              :: rho_tot_gspace
    TYPE(qs_rho_type), POINTER               :: rho, rho_xc
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    INTEGER                                  :: ispin, n_electrons
    REAL(dp)                                 :: tot1_h, tot1_s, tot_rho_r
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(qs_charges_type), POINTER           :: qs_charges

    NULLIFY(qs_charges, atomic_kind_set, cell)
    CALL get_qs_env(qs_env=qs_env, &
         atomic_kind_set=atomic_kind_set, &
         cell=cell,qs_charges=qs_charges,error=error)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
         nelectron=n_electrons)
    n_electrons = n_electrons - qs_env%dft_control%charge

    tot_rho_r = accurate_sum(rho%tot_rho_r)
    IF(output_unit>0) THEN
      WRITE (UNIT=output_unit,FMT="(/,(T3,A,T41,2F20.10))")&
           "Total electronic density (r-space): ",&
           tot_rho_r,&
           tot_rho_r + &
           REAL(n_electrons,dp),&
           "Total core charge density (r-space):",&
           qs_charges%total_rho_core_rspace,&
           qs_charges%total_rho_core_rspace - REAL(n_electrons,dp)
    END IF
    IF(qs_env%dft_control%qs_control%gapw ) THEN
       tot1_h =  qs_charges%total_rho1_hard(1)
       tot1_s =  qs_charges%total_rho1_soft(1)
       DO ispin=2,qs_env%dft_control%nspins
          tot1_h = tot1_h + qs_charges%total_rho1_hard(ispin)
          tot1_s = tot1_s + qs_charges%total_rho1_soft(ispin)
       END DO
       IF(output_unit>0) THEN
         WRITE (UNIT=output_unit,FMT="(/,(T3,A,T41,2F20.10))")&
              "Hard and soft densities (Lebedev):",&
              tot1_h, tot1_s
         WRITE (UNIT=output_unit,FMT="(T3,A,T61,F20.10)")&
              "Total Rho_soft + Rho1_hard - Rho1_soft (r-space): ",&
              tot_rho_r+ tot1_h - tot1_s ,&
              "Total charge density (r-space):      ",&
              tot_rho_r+ tot1_h - tot1_s &
              + qs_charges%total_rho_core_rspace,&
              "Total Rho_soft + Rho0_soft (g-space):",&
              qs_charges%total_rho_gspace
       END IF
       qs_charges%background=tot_rho_r+ tot1_h - tot1_s+&
                                    qs_charges%total_rho_core_rspace
    ELSE IF( qs_env%dft_control%qs_control%gapw_xc) THEN
       tot1_h =  qs_charges%total_rho1_hard(1)
       tot1_s =  qs_charges%total_rho1_soft(1)
       DO ispin=2,qs_env%dft_control%nspins
          tot1_h = tot1_h + qs_charges%total_rho1_hard(ispin)
          tot1_s = tot1_s + qs_charges%total_rho1_soft(ispin)
       END DO
       IF(output_unit>0) THEN
         WRITE (UNIT=output_unit,FMT="(/,(T3,A,T41,2F20.10))")&
              "Hard and soft densities (Lebedev):",&
              tot1_h, tot1_s
         WRITE (UNIT=output_unit,FMT="(T3,A,T61,F20.10)")&
              "Total Rho_soft + Rho1_hard - Rho1_soft (r-space): ",&
              accurate_sum(rho_xc%tot_rho_r)+ tot1_h - tot1_s
       END IF
       qs_charges%background=tot_rho_r+ &
                                    qs_charges%total_rho_core_rspace
    ELSE
       IF(output_unit>0) THEN
         WRITE (UNIT=output_unit,FMT="(T3,A,T61,F20.10)")&
              "Total charge density (r-space):     ",&
              tot_rho_r+&
              qs_charges%total_rho_core_rspace,&
              "Total charge density (g-space):     ",&
              qs_charges%total_rho_gspace
       END IF
       qs_charges%background=tot_rho_r+ &
                                    qs_charges%total_rho_core_rspace
    END IF
    qs_charges%background=qs_charges%background/cell%deth

  END SUBROUTINE print_densities

! *****************************************************************************
!> \brief updates the Kohn Sham matrix of the given qs_env (facility method)
!> \param ks_env the ks_env that old all the temporary objects that
!>        the calculation of the KS matrix needs
!> \param qs_env the qs_env to update
!> \param calculate_forces if true calculate the quantities needed
!>        to calculate the forces. Defaults to false.
!> \param just_energy if true updates the energies but not the
!>        ks matrix. Defaults to false
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      4.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE qs_ks_update_qs_env(ks_env,qs_env, kg_env,calculate_forces,&
                               just_energy,error)
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(kg_environment_type), OPTIONAL, &
      POINTER                                :: kg_env
    LOGICAL, INTENT(IN), OPTIONAL            :: calculate_forces, just_energy
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_ks_update_qs_env', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin
    LOGICAL                                  :: c_forces, did_update, &
                                                energy_only, failure
    REAL(KIND=dp)                            :: ecore
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(kg_fm_p_type), POINTER              :: kg_fm_set
    TYPE(kg_sub_pw_env_type), POINTER        :: kg_sub_pw_env
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_rho_type), POINTER               :: rho, rho_xc
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_h, matrix_ks

  c_forces=.FALSE.
  energy_only = .FALSE.
  IF (PRESENT(just_energy)) energy_only=just_energy
  IF (PRESENT(calculate_forces)) c_forces=calculate_forces

  IF (c_forces) THEN
     CALL timeset(routineN//' (forces)',handle)
  ELSE
     CALL timeset(routineN,handle)
  ENDIF

  failure=.FALSE.
  did_update=.FALSE.
  NULLIFY(matrix_ks,matrix_h,energy, logger, rho, rho_xc, dft_control, para_env)
  logger => cp_error_get_logger(error)

  CPPrecondition(ASSOCIATED(ks_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(ks_env%ref_count>0,cp_failure_level,routineP,error,failure)

  IF (.NOT. failure) THEN
     IF (ks_env%rho_changed .OR. ks_env%s_mstruct_changed .OR. &
          ks_env%potential_changed .OR. &
          (c_forces.AND..NOT.ks_env%forces_up_to_date)) THEN
        did_update=.TRUE.
        CALL get_qs_env(qs_env,matrix_ks=matrix_ks,rho=rho,&
             rho_xc=rho_xc,energy=energy,&
             dft_control=dft_control,&
             matrix_h=matrix_h,para_env=para_env,error=error)

        ! update ecore
        energy%core=0.0_dp
        DO ispin=1,dft_control%nspins
          CALL calculate_ecore(h=matrix_h(1)%matrix,&
                               p=rho%rho_ao(ispin)%matrix,&
                               ecore=ecore,&
                               para_env=para_env)
          energy%core=energy%core+ecore
        END DO
        ! the ks matrix will be rebuilt so this is fine now
        ks_env%potential_changed=.FALSE.
        IF ( dft_control%qs_control%semi_empirical ) THEN
          CALL build_se_fock_matrix(ks_env,qs_env=qs_env,&
               ks_matrix=matrix_ks,rho=rho,energy=energy,&
               calculate_forces=c_forces,just_energy=energy_only,&
               error=error)
        ELSEIF ( dft_control%qs_control%dftb ) THEN
          CALL build_dftb_ks_matrix(ks_env,qs_env=qs_env,&
               ks_matrix=matrix_ks,rho=rho,energy=energy,&
               calculate_forces=c_forces,just_energy=energy_only,&
               error=error)
        ELSEIF (PRESENT(kg_env)) THEN
          CALL get_qs_env(qs_env=qs_env, kg_sub_pw_env=kg_sub_pw_env,error=error)
          CALL get_kg_env(kg_env=kg_env, kg_fm_set=kg_fm_set,error=error)
          CPPostcondition(kg_fm_set%ref_count>=1,cp_failure_level,routineP,error,failure)
          CPPostcondition(kg_sub_pw_env%ref_count>=1,cp_failure_level,routineP,error,failure)
          CALL qs_ks_build_kohn_sham_matrix(ks_env,qs_env=qs_env,&
               ks_matrix=matrix_ks, rho=rho, energy=energy,&
               calculate_forces=c_forces,just_energy=energy_only,&
               kg_gpw = .TRUE.,kg_sub_pw_env=kg_sub_pw_env, kg_fm_set=kg_fm_set, &
               error=error)
        ELSEIF (dft_control%qs_control%gapw_xc) THEN
          CALL qs_ks_build_kohn_sham_matrix(ks_env,qs_env=qs_env,&
               ks_matrix=matrix_ks, rho=rho, energy=energy,&
               calculate_forces=c_forces,just_energy=energy_only,&
               rho_xc=rho_xc,error=error)
        ELSE
          CALL qs_ks_build_kohn_sham_matrix(ks_env,qs_env=qs_env,&
               ks_matrix=matrix_ks, rho=rho, energy=energy,&
               calculate_forces=c_forces,just_energy=energy_only,&
               error=error)
        END IF
        IF(.NOT.energy_only) THEN
           ks_env%rho_changed=.FALSE.
           ks_env%s_mstruct_changed=.FALSE.
           ks_env%forces_up_to_date=ks_env%forces_up_to_date.or.c_forces
        END IF
     END IF
  END IF

  IF (cp_debug.AND.debug_this_module) THEN
     IF (did_update) THEN
        CALL cp_log(logger,cp_note_level,routineP,&
             "did update")
     ELSE
        CALL cp_log(logger,cp_note_level,routineP,&
             "did NOT update")
     END IF
  END IF

  CALL timestop(handle)

END SUBROUTINE qs_ks_update_qs_env

! *****************************************************************************
!> \brief tells that some of the things relevant to the ks calculation
!>      did change. has to be called when changes happen otherwise
!>      the calculation will give wrong results.
!> \param ks_env the environement that is informed about the changes
!> \param s_mstruct_changed if true it means that the structure of the
!>        overlap matrix has changed
!>        (atoms have moved)
!> \param rho_changed if true it means that the density has changed
!> \param full_reset if true everything has changed
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      4.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE qs_ks_did_change(ks_env,s_mstruct_changed,rho_changed,potential_changed,full_reset,error)
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    LOGICAL, INTENT(in), OPTIONAL            :: s_mstruct_changed, &
                                                rho_changed, &
                                                potential_changed, full_reset
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ks_did_change', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)
    

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(ks_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ks_env%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(rho_changed)) THEN
          IF (rho_changed) ks_env%rho_changed=.TRUE.
       END IF
       IF (PRESENT(potential_changed)) THEN
          IF (potential_changed) ks_env%potential_changed=.TRUE.
       END IF
       IF (PRESENT(s_mstruct_changed)) THEN
          IF (s_mstruct_changed) THEN
             ks_env%s_mstruct_changed=.TRUE.
          END IF
       END IF
       IF (PRESENT(full_reset)) THEN
          IF (full_reset) THEN
             ks_env%potential_changed=.TRUE.
             ks_env%s_mstruct_changed=.TRUE.
          END IF
       END IF
       IF (ks_env%s_mstruct_changed.OR.ks_env%potential_changed.or.&
            ks_env%rho_changed) ks_env%forces_up_to_date=.FALSE.
    END IF

    CALL  timestop(handle)

END SUBROUTINE qs_ks_did_change

! *****************************************************************************
!> \brief calculates and allocates the xc potential, already reducing it to
!>      the dependence on rho and the one on tau
!> \param vxc_rho will contain the v_xc part that depend on rho
!>        (if one of the choosen xc functionals has it it is allocated and you
!>        are responsible for it)
!> \param vxc_tau will contain the kinetic (tau) part of v_xc
!>        (if one of the choosen xc functionals has it it is allocated and you
!>        are responsible for it)
!> \param qs_env qs_environment to get all the needed things (rho,...)
!> \param calculate_forces if the forces should be calculated (so that you
!>        might be able to add forces to the atoms is you do strange stuff in
!>        the xc
!> \param just_energy if true calculates just the energy, and does not
!>        allocate v_*_rspace
!> \param harris If true uses the rho structure of the Harris environment
!> \param nsc_force Uses rho_diff in case the non-self-consistent force should be
!>                   calculated. Works only in combination with the Harris functional!
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      - 05.2002 modified to use the mp_allgather function each pe
!>        computes only part of the grid and this is broadcasted to all
!>        instead of summed.
!>        This scales significantly better (e.g. factor 3 on 12 cpus
!>        32 H2O) [Joost VdV]
!>      - moved to qs_ks_methods [fawzi]
!>      - sic alterations [Joost VandeVondele]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE qs_vxc_create(vxc_rho,vxc_tau, exc, qs_env,&
     just_energy, harris, nsc_force, gapw_xc, epr_xc, adiabatic_rescale_factor, &
     basis_set_id, error)
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: vxc_rho, vxc_tau
    REAL(KIND=dp)                            :: exc
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(in), OPTIONAL            :: just_energy, harris, &
                                                nsc_force, gapw_xc, epr_xc
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: adiabatic_rescale_factor
    INTEGER, INTENT(IN), OPTIONAL            :: basis_set_id
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_vxc_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, &
                                                my_basis_set_id, myfun, &
                                                nelec_spin(2), stat
    LOGICAL :: do_adiabatic_rescaling, failure, harris_flag, my_epr_xc, &
      my_gapw_xc, my_just_energy, nsc_flag, sic_scaling_b_zero, uf_grid
    REAL(KIND=dp)                            :: exc_m, &
                                                my_adiabatic_rescale_factor, &
                                                my_scaling, nelec_s_inv
    TYPE(cell_type), POINTER                 :: cell
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(harris_env_type), POINTER           :: harris_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: my_vxc_rho, my_vxc_tau, &
                                                rho_g, rho_m_gspace, &
                                                rho_m_rspace, rho_r, tau
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool, xc_pw_pool
    TYPE(pw_type), POINTER                   :: tmp_g, tmp_g2, tmp_pw
    TYPE(qs_rho_type), POINTER               :: rho_struct
    TYPE(section_vals_type), POINTER         :: input, xc_section
    TYPE(virial_type), POINTER               :: virial

  CALL timeset( routineN ,handle)

  failure=.FALSE.
  CPPrecondition(.NOT.ASSOCIATED(vxc_rho),cp_failure_level,routineP,error,failure)
  CPPrecondition(.NOT.ASSOCIATED(vxc_tau),cp_failure_level,routineP,error,failure)
  NULLIFY(dft_control,rho_struct,pw_env,auxbas_pw_pool,xc_pw_pool,harris_env,&
       cell, my_vxc_rho, tmp_pw,tmp_g,tmp_g2,&
       my_vxc_tau, rho_g, rho_r, tau, rho_m_rspace,rho_m_gspace, input)

  my_just_energy=.FALSE.
  IF (PRESENT(just_energy)) my_just_energy=just_energy
  my_gapw_xc=.FALSE.
  IF(PRESENT(gapw_xc)) my_gapw_xc=gapw_xc
  my_epr_xc=.FALSE.
  IF(PRESENT(epr_xc)) my_epr_xc=epr_xc
  IF (PRESENT(harris)) THEN
     harris_flag = harris
  ELSE
     harris_flag = .FALSE.
  END IF
  IF (PRESENT(nsc_force)) THEN
     nsc_flag = nsc_force
  ELSE
     nsc_flag = .FALSE.
  END IF

  my_adiabatic_rescale_factor = 1.0_dp
  do_adiabatic_rescaling = .FALSE.
  IF( PRESENT(adiabatic_rescale_factor)) THEN
    my_adiabatic_rescale_factor = adiabatic_rescale_factor
    do_adiabatic_rescaling = .TRUE.
  END IF

  IF( PRESENT(basis_set_id) ) THEN
    my_basis_set_id = basis_set_id
  ELSE
    my_basis_set_id = use_orb_basis_set
  END IF


  CALL get_qs_env(qs_env=qs_env, dft_control=dft_control,input=input, &
       pw_env=pw_env, cell=cell,virial=virial, error=error)
  IF(my_epr_xc) THEN
     xc_section => section_vals_get_subs_vals(input,&
                   "PROPERTIES%LINRES%EPR%PRINT%G_TENSOR%XC",error=error)
  ELSE
     xc_section => section_vals_get_subs_vals(input,"DFT%XC",error=error)
  END IF
  CALL section_vals_val_get(xc_section,"XC_FUNCTIONAL%_SECTION_PARAMETERS_",&
       i_val=myfun,error=error)
  IF (myfun/=xc_none) THEN
     ! would be better to pass these as arguments
     IF (harris_flag) THEN
        CALL get_qs_env(qs_env=qs_env,  harris_env=harris_env, error=error)
        IF (nsc_flag) THEN
           rho_struct => harris_env%rho_diff
        ELSE
           rho_struct => harris_env%rho
        END IF
     ELSEIF(my_gapw_xc) THEN
        CALL get_qs_env(qs_env=qs_env, rho_xc=rho_struct, error=error)
     ELSE
        SELECT CASE (my_basis_set_id)
        CASE (use_orb_basis_set)
          CALL get_qs_env(qs_env=qs_env, rho=rho_struct, error=error)
        CASE (use_aux_fit_basis_set)
          CALL get_qs_env(qs_env=qs_env, rho_aux_fit=rho_struct, error=error)
        END SELECT
     END IF

     ! test if the real space density is available
     CPPrecondition(ASSOCIATED(rho_struct),cp_failure_level,routineP,error,failure)
     CPPrecondition(rho_struct%ref_count>0,cp_failure_level,routineP,error,failure)
     CPPrecondition(rho_struct%rho_r_valid,cp_failure_level,routineP,error,failure)
     CALL cp_assert( dft_control%nspins == 1 .OR. dft_control%nspins == 2,&
          cp_failure_level,cp_assertion_failed,routineP,&
          "nspins must be 1 or 2",error,failure)
     ! there are some options related to SIC here.
     ! Normal DFT computes E(rho_alpha,rho_beta) (or its variant E(2*rho_alpha) for non-LSD)
     ! SIC can             E(rho_alpha,rho_beta)-b*(E(rho_alpha,rho_beta)-E(rho_beta,rho_beta))
     ! or compute          E(rho_alpha,rho_beta)-b*E(rho_alpha-rho_beta,0)

     ! my_scaling is the scaling needed of the standard E(rho_alpha,rho_beta) term
     my_scaling=1.0_dp
     SELECT CASE (dft_control%sic_method_id)
     CASE ( sic_none )
        ! all fine
     CASE ( sic_mauri_spz, sic_ad )
        ! no idea yet what to do here in that case
        CPPrecondition(.NOT.rho_struct%tau_r_valid ,cp_failure_level,routineP,error,failure)
     CASE ( sic_mauri_us )
        my_scaling=1.0_dp-dft_control % sic_scaling_b
        ! no idea yet what to do here in that case
        CPPrecondition(.NOT.rho_struct%tau_r_valid ,cp_failure_level,routineP,error,failure)
     CASE ( sic_eo )
        ! NOTHING TO BE DONE
     CASE DEFAULT
        ! this case has not yet been treated here
        CALL cp_assert(.FALSE., cp_failure_level,cp_assertion_failed,routineP,"NYI",error,failure)
     END SELECT

     IF (dft_control % sic_scaling_b .EQ. 0.0_dp) THEN
        sic_scaling_b_zero = .TRUE.
     ELSE
        sic_scaling_b_zero = .FALSE.
     ENDIF

     IF ( .NOT. failure ) THEN
        CALL pw_env_get(pw_env,xc_pw_pool=xc_pw_pool,auxbas_pw_pool=auxbas_pw_pool,&
             error=error)
        uf_grid = .NOT. pw_grid_compare(auxbas_pw_pool%pw_grid,xc_pw_pool%pw_grid)

        ALLOCATE(rho_r(dft_control%nspins),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        IF (.not.uf_grid) THEN
           DO ispin=1,dft_control%nspins
              rho_r(ispin)%pw => rho_struct%rho_r(ispin)%pw
           END DO

           IF (rho_struct%tau_r_valid) THEN
              ALLOCATE(tau(dft_control%nspins),stat=stat)
              CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
              DO ispin=1,dft_control%nspins
                 tau(ispin)%pw => rho_struct%tau_r(ispin)%pw
              END DO
           END IF

           ! for gradient corrected functional the density in g space might
           ! be useful so if we have it, we pass it in
           IF ( rho_struct%rho_g_valid ) THEN
              ALLOCATE(rho_g(dft_control%nspins),stat=stat)
              CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
              DO ispin=1,dft_control%nspins
                 rho_g(ispin)%pw => rho_struct%rho_g(ispin)%pw
              END DO
           END IF
        ELSE
           CPPrecondition(rho_struct%rho_g_valid,cp_failure_level,routineP,error,failure)
           ALLOCATE(rho_g(dft_control%nspins),stat=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           DO ispin=1,dft_control%nspins
              CALL pw_pool_create_pw(xc_pw_pool,rho_g(ispin)%pw,&
                   in_space=RECIPROCALSPACE, use_data=COMPLEXDATA1D,error=error)
              CALL pw_transfer(rho_struct%rho_g(ispin)%pw,rho_g(ispin)%pw, error=error)
           END DO
           DO ispin=1,dft_control%nspins
              CALL pw_pool_create_pw(xc_pw_pool,rho_r(ispin)%pw,&
                   in_space=REALSPACE, use_data=REALDATA3D,error=error)
              CALL pw_transfer(rho_g(ispin)%pw,rho_r(ispin)%pw, error=error)
           END DO
           CPPrecondition(.NOT.rho_struct%tau_r_valid,cp_failure_level,routineP,error,failure)
           IF (rho_struct%tau_r_valid) THEN
              ALLOCATE(tau(dft_control%nspins),stat=stat)
              CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
              DO ispin=1,dft_control%nspins
                 CALL pw_pool_create_pw(xc_pw_pool,tau(ispin)%pw,&
                      in_space=REALSPACE, use_data=REALDATA3D,error=error)

                 CALL pw_pool_create_pw(xc_pw_pool,tmp_g,&
                      in_space=RECIPROCALSPACE,use_data=COMPLEXDATA1D,error=error)
                 CALL pw_pool_create_pw(auxbas_pw_pool,tmp_g2,&
                      in_space=RECIPROCALSPACE,use_data=COMPLEXDATA1D,error=error)
                 CALL pw_transfer(tau(ispin)%pw,tmp_g, error=error)
                 CALL pw_transfer(tmp_g,tmp_g2, error=error)
                 CALL pw_transfer(tmp_g2,tmp_pw, error=error)
                 CALL pw_pool_give_back_pw(auxbas_pw_pool,tmp_g2,error=error)
                 CALL pw_pool_give_back_pw(xc_pw_pool,tmp_g,error=error)
              END DO
           END IF
        END IF
        !
        ! here the rho_r, rho_g, tau is what it should be
        ! we get back the right my_vxc_rho and my_vxc_tau as required
        !
        IF (my_just_energy) THEN
           exc=xc_exc_calc(rho_r=rho_r,tau=tau,&
                rho_g=rho_g, xc_section=xc_section,&
                cell=cell, pw_pool=xc_pw_pool,&
                error=error)

        ELSE
           CALL xc_vxc_pw_create1(vxc_rho=my_vxc_rho,vxc_tau=my_vxc_tau, rho_r=rho_r,&
                rho_g=rho_g,tau=tau,exc=exc,&
                xc_section=xc_section,&
                cell=cell, pw_pool=xc_pw_pool,&
                error=error,&
                virial=virial)
        END IF

        !! Apply rescaling to the potential if requested
        IF(.NOT. my_just_energy) THEN
          IF(do_adiabatic_rescaling) THEN
            IF( ASSOCIATED(my_vxc_rho)) THEN
              DO ispin=1,SIZE(my_vxc_rho)
                my_vxc_rho(ispin)%pw%cr3d=my_vxc_rho(ispin)%pw%cr3d*my_adiabatic_rescale_factor
              END DO
            END IF
          END IF
        END IF

        IF (my_scaling .NE. 1.0_dp) THEN
           exc=exc * my_scaling
           IF (ASSOCIATED(my_vxc_rho)) THEN
              DO ispin=1,SIZE(my_vxc_rho)
                 my_vxc_rho(ispin)%pw%cr3d=my_vxc_rho(ispin)%pw%cr3d*my_scaling
              ENDDO
           ENDIF
           IF (ASSOCIATED(my_vxc_tau)) THEN
              DO ispin=1,SIZE(my_vxc_tau)
                 my_vxc_tau(ispin)%pw%cr3d=my_vxc_tau(ispin)%pw%cr3d*my_scaling
              ENDDO
           ENDIF
        ENDIF

        ! we have pw data for the xc, qs_ks requests coeff structure, here we transfer
        ! pw -> coeff
        IF (ASSOCIATED(my_vxc_rho)) THEN
           ALLOCATE(vxc_rho(dft_control%nspins),stat=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           DO ispin=1,dft_control%nspins
              vxc_rho(ispin)%pw => my_vxc_rho(ispin)%pw
           END DO
           DEALLOCATE(my_vxc_rho,stat=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        END IF
        IF (ASSOCIATED(my_vxc_tau)) THEN
           ALLOCATE(vxc_tau(dft_control%nspins),stat=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           DO ispin=1,dft_control%nspins
              vxc_tau(ispin)%pw => my_vxc_tau(ispin)%pw
           END DO
           DEALLOCATE(my_vxc_tau,stat=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        END IF

        ! compute again the xc but now for Exc(m,o) and the opposite sign
        IF (dft_control%sic_method_id .EQ. sic_mauri_spz .AND. .NOT. sic_scaling_b_zero) THEN
           ALLOCATE(rho_m_rspace(2),rho_m_gspace(2))
           CALL pw_pool_create_pw(xc_pw_pool, rho_m_gspace(1)%pw,&
                use_data = COMPLEXDATA1D,&
                in_space = RECIPROCALSPACE, error=error)
           CALL pw_pool_create_pw(xc_pw_pool, rho_m_rspace(1)%pw,&
                use_data = REALDATA3D,&
                in_space = REALSPACE, error=error)
           CALL pw_copy(rho_struct%rho_r(1)%pw,rho_m_rspace(1)%pw, error=error)
           CALL pw_axpy(rho_struct%rho_r(2)%pw,rho_m_rspace(1)%pw,alpha=-1._dp, error=error)
           CALL pw_copy(rho_struct%rho_g(1)%pw,rho_m_gspace(1)%pw, error=error)
           CALL pw_axpy(rho_struct%rho_g(2)%pw,rho_m_gspace(1)%pw,alpha=-1._dp, error=error)
           ! bit sad, these will be just zero...
           CALL pw_pool_create_pw(xc_pw_pool, rho_m_gspace(2)%pw,&
                use_data = COMPLEXDATA1D,&
                in_space = RECIPROCALSPACE, error=error)
           CALL pw_pool_create_pw(xc_pw_pool, rho_m_rspace(2)%pw,&
                use_data = REALDATA3D,&
                in_space = REALSPACE, error=error)
           CALL pw_zero(rho_m_rspace(2)%pw, error=error)
           CALL pw_zero(rho_m_gspace(2)%pw, error=error)

           rho_g(1)%pw => rho_m_gspace(1)%pw
           rho_g(2)%pw => rho_m_gspace(2)%pw
           rho_r(1)%pw => rho_m_rspace(1)%pw
           rho_r(2)%pw => rho_m_rspace(2)%pw

           IF (my_just_energy) THEN
              exc_m=xc_exc_calc(rho_r=rho_r,tau=tau,&
                   rho_g=rho_g, xc_section=xc_section,&
                   cell=cell, pw_pool=xc_pw_pool,&
                   error=error)
           ELSE
              CALL xc_vxc_pw_create1(vxc_rho=my_vxc_rho,vxc_tau=my_vxc_tau, rho_r=rho_r,&
                   rho_g=rho_g,tau=tau,exc=exc_m,&
                   xc_section=xc_section,&
                   cell=cell, pw_pool=xc_pw_pool,&
                   error=error)
           END IF

           exc = exc - dft_control % sic_scaling_b * exc_m

           ! and take care of the potential only vxc_rho is taken into account
           IF (.NOT. my_just_energy) THEN
              vxc_rho(1)%pw%cr3d=vxc_rho(1)%pw%cr3d-dft_control % sic_scaling_b *&
                   my_vxc_rho(1)%pw%cr3d
              vxc_rho(2)%pw%cr3d=vxc_rho(2)%pw%cr3d+dft_control % sic_scaling_b *&
                   my_vxc_rho(1)%pw%cr3d ! 1=m
              CALL pw_release(my_vxc_rho(1)%pw,error=error)
              CALL pw_release(my_vxc_rho(2)%pw,error=error)
              DEALLOCATE(my_vxc_rho,stat=stat)
              CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           ENDIF

           DO ispin=1,2
              CALL pw_pool_give_back_pw(xc_pw_pool,rho_m_rspace(ispin)%pw,&
                   error=error)
              CALL pw_pool_give_back_pw(xc_pw_pool,rho_m_gspace(ispin)%pw,&
                   error=error)
           ENDDO
           DEALLOCATE(rho_m_rspace)
           DEALLOCATE(rho_m_gspace)

        ENDIF

        ! now we have - sum_s N_s * Exc(rho_s/N_s,0)
        IF ( dft_control % sic_method_id .EQ. sic_ad  .AND. .NOT. sic_scaling_b_zero ) THEN

           ! find out how many elecs we have
           CALL get_qs_env(qs_env,mos=mo_array,error=error)
           CALL get_mo_set(mo_set=mo_array(1)%mo_set,nelectron=nelec_spin(1))
           CALL get_mo_set(mo_set=mo_array(2)%mo_set,nelectron=nelec_spin(2))

           ALLOCATE(rho_m_rspace(2),rho_m_gspace(2))
           DO ispin=1,2
              CALL pw_pool_create_pw(xc_pw_pool, rho_m_gspace(ispin)%pw,&
                   use_data = COMPLEXDATA1D,&
                   in_space = RECIPROCALSPACE, error=error)
              CALL pw_pool_create_pw(xc_pw_pool, rho_m_rspace(ispin)%pw,&
                   use_data = REALDATA3D,&
                   in_space = REALSPACE, error=error)
           ENDDO

           rho_g(1)%pw => rho_m_gspace(1)%pw
           rho_g(2)%pw => rho_m_gspace(2)%pw
           rho_r(1)%pw => rho_m_rspace(1)%pw
           rho_r(2)%pw => rho_m_rspace(2)%pw

           DO ispin=1,2
              IF (nelec_spin(ispin) .GT. 0 ) THEN
                 nelec_s_inv=1.0_dp/REAL(nelec_spin(ispin),KIND=dp)
              ELSE
                 ! does it matter if there are no electrons with this spin (H) ?
                 nelec_s_inv=0.0_dp
              ENDIF
              CALL pw_copy(rho_struct%rho_r(ispin)%pw,rho_m_rspace(1)%pw, error=error)
              CALL pw_copy(rho_struct%rho_g(ispin)%pw,rho_m_gspace(1)%pw, error=error)
              CALL pw_scale(rho_m_rspace(1)%pw,nelec_s_inv, error=error)
              CALL pw_scale(rho_m_gspace(1)%pw,nelec_s_inv, error=error)
              CALL pw_zero(rho_m_rspace(2)%pw, error=error)
              CALL pw_zero(rho_m_gspace(2)%pw, error=error)

              IF (my_just_energy) THEN
                 exc_m=xc_exc_calc(rho_r=rho_r,tau=tau,&
                      rho_g=rho_g, xc_section=xc_section,&
                      cell=cell, pw_pool=xc_pw_pool,&
                      error=error)
              ELSE
                 CALL xc_vxc_pw_create1(vxc_rho=my_vxc_rho,vxc_tau=my_vxc_tau, rho_r=rho_r,&
                      rho_g=rho_g,tau=tau,exc=exc_m,&
                      xc_section=xc_section,&
                      cell=cell, pw_pool=xc_pw_pool,&
                      error=error)
              END IF

              exc = exc - dft_control % sic_scaling_b * nelec_spin(ispin) * exc_m

              ! and take care of the potential only vxc_rho is taken into account
              IF (.NOT. my_just_energy) THEN
                 vxc_rho(ispin)%pw%cr3d=vxc_rho(ispin)%pw%cr3d-dft_control % sic_scaling_b *&
                      my_vxc_rho(1)%pw%cr3d
                 CALL pw_release(my_vxc_rho(1)%pw,error=error)
                 CALL pw_release(my_vxc_rho(2)%pw,error=error)
                 DEALLOCATE(my_vxc_rho,stat=stat)
                 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
              ENDIF
           ENDDO

           DO ispin=1,2
              CALL pw_pool_give_back_pw(xc_pw_pool,rho_m_rspace(ispin)%pw,&
                   error=error)
              CALL pw_pool_give_back_pw(xc_pw_pool,rho_m_gspace(ispin)%pw,&
                   error=error)
           ENDDO
           DEALLOCATE(rho_m_rspace)
           DEALLOCATE(rho_m_gspace)

        ENDIF

        ! compute again the xc but now for Exc(n_down,n_down)
        IF (dft_control%sic_method_id .EQ. sic_mauri_us .AND. .NOT. sic_scaling_b_zero ) THEN
           rho_r(1)%pw => rho_struct%rho_r(2)%pw
           rho_r(2)%pw => rho_struct%rho_r(2)%pw
           IF ( rho_struct%rho_g_valid ) THEN
              rho_g(1)%pw => rho_struct%rho_g(2)%pw
              rho_g(2)%pw => rho_struct%rho_g(2)%pw
           ENDIF

           IF (my_just_energy) THEN
              exc_m=xc_exc_calc(rho_r=rho_r,tau=tau,&
                   rho_g=rho_g, xc_section=xc_section,&
                   cell=cell, pw_pool=xc_pw_pool,&
                   error=error)
           ELSE
              CALL xc_vxc_pw_create1(vxc_rho=my_vxc_rho,vxc_tau=my_vxc_tau, rho_r=rho_r,&
                   rho_g=rho_g,tau=tau,exc=exc_m,&
                   xc_section=xc_section,&
                   cell=cell, pw_pool=xc_pw_pool,&
                   error=error)
           END IF

           exc = exc + dft_control % sic_scaling_b * exc_m

           ! and take care of the potential
           IF (.NOT. my_just_energy) THEN
              ! both go to minority spin
              vxc_rho(2)%pw%cr3d = vxc_rho(2)%pw%cr3d + &
                   2.0_dp * dft_control % sic_scaling_b * my_vxc_rho(1)%pw%cr3d
              CALL pw_release(my_vxc_rho(1)%pw,error=error)
              CALL pw_release(my_vxc_rho(2)%pw,error=error)
              DEALLOCATE(my_vxc_rho)
              CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           ENDIF

        ENDIF

        !
        ! cleanups
        !
        IF (uf_grid) THEN
           DO ispin=1,SIZE(rho_r)
              CALL pw_pool_give_back_pw(xc_pw_pool,rho_r(ispin)%pw,error=error)
           END DO
           IF (ASSOCIATED(vxc_rho)) THEN
              DO ispin=1,SIZE(vxc_rho)
                 CALL pw_pool_create_pw(auxbas_pw_pool,tmp_pw,&
                      in_space=REALSPACE,use_data=REALDATA3D,error=error)

                 CALL pw_pool_create_pw(xc_pw_pool,tmp_g,&
                      in_space=RECIPROCALSPACE,use_data=COMPLEXDATA1D,error=error)
                 CALL pw_pool_create_pw(auxbas_pw_pool,tmp_g2,&
                      in_space=RECIPROCALSPACE,use_data=COMPLEXDATA1D,error=error)
                 CALL pw_transfer(vxc_rho(ispin)%pw,tmp_g, error=error)
                 CALL pw_transfer(tmp_g,tmp_g2, error=error)
                 CALL pw_transfer(tmp_g2,tmp_pw, error=error)
                 CALL pw_pool_give_back_pw(auxbas_pw_pool,tmp_g2,error=error)
                 CALL pw_pool_give_back_pw(xc_pw_pool,tmp_g,error=error)
                 !FM              CALL pw_zero(tmp_pw,error=error)
                 !FM              CALL pw_restrict_s3(vxc_rho(ispin)%pw,tmp_pw,&
                 !FM                   auxbas_pw_pool,param_section=interp_section,error=error)
                 CALL pw_pool_give_back_pw(xc_pw_pool,vxc_rho(ispin)%pw,error=error)
                 vxc_rho(ispin)%pw => tmp_pw
                 NULLIFY(tmp_pw)
              END DO
           END IF
           IF (ASSOCIATED(vxc_tau)) THEN
              DO ispin=1,SIZE(vxc_tau)
                 CALL pw_pool_create_pw(auxbas_pw_pool,tmp_pw,&
                      in_space=REALSPACE,use_data=REALDATA3D,error=error)

                 CALL pw_pool_create_pw(xc_pw_pool,tmp_g,&
                      in_space=RECIPROCALSPACE,use_data=COMPLEXDATA1D,error=error)
                 CALL pw_pool_create_pw(auxbas_pw_pool,tmp_g2,&
                      in_space=RECIPROCALSPACE,use_data=COMPLEXDATA1D,error=error)
                 CALL pw_transfer(vxc_tau(ispin)%pw,tmp_g, error=error)
                 CALL pw_transfer(tmp_g,tmp_g2, error=error)
                 CALL pw_transfer(tmp_g2,tmp_pw, error=error)
                 CALL pw_pool_give_back_pw(auxbas_pw_pool,tmp_g2,error=error)
                 CALL pw_pool_give_back_pw(xc_pw_pool,tmp_g,error=error)
                 !FM              CALL pw_zero(tmp_pw,error=error)
                 !FM              CALL pw_restrict_s3(vxc_rho(ispin)%pw,tmp_pw,&
                 !FM                   auxbas_pw_pool,param_section=interp_section,error=error)
                 CALL pw_pool_give_back_pw(xc_pw_pool,vxc_tau(ispin)%pw,error=error)
                 vxc_tau(ispin)%pw => tmp_pw
                 NULLIFY(tmp_pw)
              END DO
           END IF

        END IF
        DEALLOCATE(rho_r,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        IF (ASSOCIATED(rho_g)) THEN
           IF (uf_grid) THEN
              DO ispin=1,SIZE(rho_g)
                 CALL pw_pool_give_back_pw(xc_pw_pool,rho_g(ispin)%pw,error=error)
              END DO
           END IF
           DEALLOCATE(rho_g,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(tau)) THEN
           IF (uf_grid) THEN
              DO ispin=1,SIZE(tau)
                 CALL pw_pool_give_back_pw(xc_pw_pool,tau(ispin)%pw,error=error)
              END DO
           END IF
           DEALLOCATE(tau,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF

     END IF
  END IF
  CALL timestop(handle)

END SUBROUTINE qs_vxc_create

! *****************************************************************************
!> \brief Calculate the W matrix from the MO eigenvectors, MO eigenvalues,
!>       and the MO occupation numbers. Only works if they are eigenstates
!> \param mo_set type containing the full matrix of the MO and the eigenvalues
!> \param w_matrix sparse matrix
!>        error
!> \par History
!>         Creation (03.03.03,MK)
!>         Modification that computes it as a full block, several times (e.g. 20)
!>               faster at the cost of some additional memory
!> \author MK
! *****************************************************************************
  SUBROUTINE calculate_w_matrix_1(mo_set,w_matrix,error)

    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(real_matrix_type), POINTER          :: w_matrix
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_w_matrix_1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, imo, istat
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: eigocc
    TYPE(cp_fm_type), POINTER                :: weighted_vectors

    CALL timeset(routineN,handle)
    failure=.FALSE.
    NULLIFY(weighted_vectors)

    CALL set_matrix(w_matrix,0.0_dp)
    CALL cp_fm_create(weighted_vectors,mo_set%mo_coeff%matrix_struct,"weighted_vectors",error=error)
    CALL cp_fm_to_fm(mo_set%mo_coeff,weighted_vectors,error=error)

    ! scale every column with the occupation
    ALLOCATE(eigocc(mo_set%homo),stat=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DO imo=1,mo_set%homo
       eigocc(imo) = mo_set%eigenvalues(imo)*mo_set%occupation_numbers(imo)
    ENDDO
    CALL cp_fm_column_scale(weighted_vectors,eigocc)
    DEALLOCATE(eigocc)

    CALL cp_sm_plus_fm_fm_t(sparse_matrix=w_matrix,&
                            matrix_v=mo_set%mo_coeff,&
                            matrix_g=weighted_vectors,&
                            ncol=mo_set%homo,error=error)

    CALL cp_fm_release(weighted_vectors,error=error)
    CALL timestop(handle)

  END SUBROUTINE calculate_w_matrix_1

! *****************************************************************************
!> \brief Calculate the W matrix from the MO coefs, MO derivs
!>        could overwrite the mo_derivs for increased memory efficiency
!> \param mo_set type containing the full matrix of the MO coefs
!>        mo_deriv:
!> \param w_matrix sparse matrix
!> \param s_matrix sparse matrix for the overlap
!>        error
!> \par History
!>         Creation (JV)
!> \author MK
! *****************************************************************************
  SUBROUTINE calculate_w_matrix_2(mo_set,mo_deriv,w_matrix,s_matrix,error)

    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_fm_type), POINTER                :: mo_deriv
    TYPE(real_matrix_type), POINTER          :: w_matrix, s_matrix
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_w_matrix_2', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: check_gradient = .FALSE., &
                                                do_symm = .FALSE.

    INTEGER                                  :: handle, ncol_block, &
                                                ncol_global, nrow_block, &
                                                nrow_global, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: occupation_numbers, &
                                                scaling_factor
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: gradient, h_block, h_block_t, &
                                                weighted_vectors

    CALL timeset(routineN,handle)
    failure=.FALSE.
    NULLIFY(weighted_vectors,h_block,fm_struct_tmp)

    CALL cp_fm_get_info(matrix=mo_set%mo_coeff,&
                        ncol_global=ncol_global,&
                        nrow_global=nrow_global,&
                        nrow_block=nrow_block,&
                        ncol_block=ncol_block,error=error)

    CALL cp_fm_create(weighted_vectors,mo_set%mo_coeff%matrix_struct,"weighted_vectors",error=error)
    CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ncol_global,  ncol_global=ncol_global, &
                                            para_env=mo_set%mo_coeff%matrix_struct%para_env, &
                                            context=mo_set%mo_coeff%matrix_struct%context,error=error)
    CALL cp_fm_create(h_block,fm_struct_tmp, name="h block",error=error)
    IF (do_symm) CALL cp_fm_create(h_block_t,fm_struct_tmp, name="h block t",error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)
  
    CALL get_mo_set(mo_set=mo_set, occupation_numbers=occupation_numbers )
    ALLOCATE(scaling_factor(SIZE(occupation_numbers)),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    scaling_factor=2.0_dp*occupation_numbers
    CALL cp_fm_to_fm(mo_deriv,weighted_vectors,error)
    CALL cp_fm_column_scale(weighted_vectors,scaling_factor)
    DEALLOCATE(scaling_factor,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ! the convention seems to require the half here, the factor of two is presumably taken care of
    ! internally in qs_core_hamiltonian
    CALL cp_fm_gemm('T','N',ncol_global,ncol_global,nrow_global,0.5_dp, &
                    mo_set%mo_coeff,weighted_vectors,0.0_dp,h_block,error=error)

    IF (do_symm) THEN
       ! at the minimum things are anyway symmetric, but numerically it might not be the case
       ! needs some investigation to find out if using this is better
       CALL cp_fm_transpose(h_block,h_block_t,error=error)
       CALL cp_fm_scale_and_add(0.5_dp,h_block,0.5_dp,h_block_t,error=error)
    ENDIF

    ! this could overwrite the mo_derivs to save the weighted_vectors
    CALL cp_fm_gemm('N','N',nrow_global,ncol_global,ncol_global,1.0_dp, &
                    mo_set%mo_coeff,h_block,0.0_dp,weighted_vectors,error=error)

    CALL set_matrix(w_matrix,0.0_dp)
    CALL cp_sm_plus_fm_fm_t(sparse_matrix=w_matrix,&
                            matrix_v=mo_set%mo_coeff,&
                            matrix_g=weighted_vectors,&
                            ncol=mo_set%homo,error=error)

    IF (check_gradient) THEN
       CALL cp_fm_create(gradient,mo_set%mo_coeff%matrix_struct,"gradient",error=error)
       CALL cp_sm_fm_multiply(s_matrix,weighted_vectors,&
                                          gradient, ncol_global,error=error)

       ALLOCATE(scaling_factor(SIZE(occupation_numbers)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       scaling_factor=2.0_dp*occupation_numbers
       CALL cp_fm_to_fm(mo_deriv,weighted_vectors,error)
       CALL cp_fm_column_scale(weighted_vectors,scaling_factor)
       DEALLOCATE(scaling_factor,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       
       WRITE(6,*) " maxabs difference ", MAXVAL(ABS(weighted_vectors%local_data-2.0_dp*gradient%local_data))
       CALL cp_fm_release(gradient,error=error)
    ENDIF

    IF (do_symm) CALL cp_fm_release(h_block_t,error=error)
    CALL cp_fm_release(weighted_vectors,error=error)
    CALL cp_fm_release(h_block,error=error)
    CALL timestop(handle)

  END SUBROUTINE calculate_w_matrix_2

! *****************************************************************************
  SUBROUTINE calculate_w_matrix_roks(mo_set,matrix_ks,matrix_p,matrix_w,error)

    ! Calculate the energy-weighted density matrix W if ROKS is active.
    ! The W matrix is returned in matrix_w.

    ! - Creation (04.05.06,MK)

    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(real_matrix_type), POINTER          :: matrix_ks, matrix_p, matrix_w
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_w_matrix_roks', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nao
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct
    TYPE(cp_fm_type), POINTER                :: c, ks, p, work
    TYPE(cp_para_env_type), POINTER          :: para_env

! -------------------------------------------------------------------------

    CALL timeset(routineN,handle)

    NULLIFY (c)
    NULLIFY (context)
    NULLIFY (fm_struct)
    NULLIFY (ks)
    NULLIFY (p)
    NULLIFY (para_env)
    NULLIFY (work)

    CALL get_mo_set(mo_set=mo_set,mo_coeff=c)
    CALL cp_fm_get_info(c,context=context,nrow_global=nao,para_env=para_env,&
                        error=error)
    CALL cp_fm_struct_create(fm_struct,context=context,nrow_global=nao,&
                             ncol_global=nao,para_env=para_env,error=error)
    CALL cp_fm_create(ks,fm_struct,name="Kohn-Sham matrix",error=error)
    CALL cp_fm_create(p,fm_struct,name="Density matrix",error=error)
    CALL cp_fm_create(work,fm_struct,name="Work matrix",error=error)
    CALL cp_fm_struct_release(fm_struct,error=error)
    CALL copy_sm_to_fm(matrix_ks,ks,error=error)
    CALL copy_sm_to_fm(matrix_p,p,error=error)
    CALL cp_fm_upper_to_full(p,work,error)
    CALL cp_fm_symm("L","U",nao,nao,1.0_dp,ks,p,0.0_dp,work,error=error)
    CALL cp_fm_gemm("T","N",nao,nao,nao,1.0_dp,p,work,0.0_dp,ks,error=error)
    CALL set_matrix(matrix_w,0.0_dp)
    CALL copy_fm_to_sm(ks,matrix_w)
    CALL cp_fm_release(work,error=error)
    CALL cp_fm_release(p,error=error)
    CALL cp_fm_release(ks,error=error)

    CALL timestop(handle)

  END SUBROUTINE calculate_w_matrix_roks

END MODULE qs_ks_methods
