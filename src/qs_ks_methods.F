!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
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
  USE admm_methods,                    ONLY: admm_merge_ks_matrix,&
                                             admm_merge_mo_derivs
  USE admm_types,                      ONLY: admm_type
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cell_types,                      ONLY: cell_type
  USE cp_control_types,                ONLY: becke_restraint_type,&
                                             dft_control_type
  USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_col_block_sizes, cp_dbcsr_copy, cp_dbcsr_create, &
       cp_dbcsr_distribution, cp_dbcsr_filter, cp_dbcsr_get_occupation, &
       cp_dbcsr_init, cp_dbcsr_multiply, cp_dbcsr_name, cp_dbcsr_p_type, &
       cp_dbcsr_row_block_sizes, cp_dbcsr_set, cp_dbcsr_type
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                             copy_fm_to_dbcsr,&
                                             cp_dbcsr_allocate_matrix_set,&
                                             cp_dbcsr_plus_fm_fm_t,&
                                             cp_dbcsr_sm_fm_multiply
  USE cp_ddapc,                        ONLY: qs_ks_ddapc
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
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_should_output
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE dbcsr_api,                       ONLY: array_i1d_obj,&
                                             dbcsr_distribution_obj,&
                                             dbcsr_type_symmetric
  USE dft_plus_u,                      ONLY: plus_u
  USE efield_utils,                    ONLY: efield_potential
  USE ep_qs_types,                     ONLY: ep_qs_type
  USE f77_blas
  USE hartree_local_methods,           ONLY: Vh_1c_gg_integrals
  USE hfx_admm_utils,                  ONLY: admm_density_matrix,&
                                             hfx_ks_matrix
  USE input_constants,                 ONLY: do_ppl_grid
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kahan_sum,                       ONLY: accurate_sum
  USE kg_correction,                   ONLY: kg_ekin_subset
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE lri_environment_types,           ONLY: lri_density_type,&
                                             lri_environment_type,&
                                             lri_kind_type
  USE lri_forces,                      ONLY: calculate_lri_forces
  USE lri_ks_methods,                  ONLY: calculate_lri_ks_matrix
  USE message_passing,                 ONLY: mp_sum
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
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
                                             pw_p_type,&
                                             pw_release
  USE qmmm_image_charge,               ONLY: add_image_pot_to_hartree_pot,&
                                             calculate_image_pot,&
                                             integrate_potential_devga_rspace
  USE qs_charges_types,                ONLY: qs_charges_type
  USE qs_core_energies,                ONLY: calculate_ecore
  USE qs_dftb_matrices,                ONLY: build_dftb_ks_matrix
  USE qs_efield_berry,                 ONLY: qs_efield_berry_phase
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_gapw_densities,               ONLY: prepare_gapw_den
  USE qs_integrate_potential,          ONLY: integrate_ppl_rspace,&
                                             integrate_rho_nlcc,&
                                             integrate_v_core_rspace,&
                                             integrate_v_rspace_one_center
  USE qs_ks_apply_restraints,          ONLY: qs_ks_becke_restraint,&
                                             qs_ks_mulliken_restraint,&
                                             qs_ks_s2_restraint
  USE qs_ks_atom,                      ONLY: update_ks_atom
  USE qs_ks_qmmm_methods,              ONLY: qmmm_calculate_energy,&
                                             qmmm_modify_hartree_pot
  USE qs_ks_types,                     ONLY: qs_ks_env_type,&
                                             qs_ks_release
  USE qs_ks_utils,                     ONLY: &
       calc_v_sic_rspace, calculate_zmp_potential, compute_matrix_vxc, &
       ep_v_core, low_spin_roks, print_densities, print_detailed_energy, &
       sic_explicit_orbitals, sum_up_and_integrate
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             mo_set_type
  USE qs_rho0_ggrid,                   ONLY: integrate_vhg0_rspace
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE qs_vxc,                          ONLY: qs_vxc_create
  USE qs_vxc_atom,                     ONLY: calculate_vxc_atom
  USE rtp_admm_methods,                ONLY: rtp_admm_density_matrix,&
                                             rtp_admm_merge_ks_matrix
  USE scptb_ks_matrix,                 ONLY: build_scptb_ks_matrix
  USE se_fock_matrix_control,          ONLY: build_se_fock_matrix_control
  USE surface_dipole,                  ONLY: calc_dipsurf_potential
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_types,                    ONLY: virial_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  INTERFACE calculate_w_matrix
    MODULE PROCEDURE calculate_w_matrix_1,&
                     calculate_w_matrix_roks

  END INTERFACE

  LOGICAL, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ks_methods'

  PUBLIC :: qs_ks_create, qs_ks_did_change,  calc_rho_tot_gspace, &
       qs_ks_update_qs_env, qs_ks_build_kohn_sham_matrix, calculate_w_matrix, calculate_w_matrix_ot,&
       qs_ks_allocate_basics
!***

  TYPE qs_potential_type
    TYPE(pw_p_type), POINTER :: v_efield_rspace, v_hartree_gspace, &
          vppl_rspace, v_sic_rspace, v_spin_ddapc_rest_r
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: v_rspace_new, &
                                                v_rspace_new_aux_fit, &
                                                v_tau_rspace, &
                                                v_tau_rspace_aux_fit
  END TYPE qs_potential_type

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
                                          just_energy,print_active,rho_xc,error)

    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(qs_energy_type), POINTER            :: energy
    LOGICAL, INTENT(in)                      :: calculate_forces, just_energy
    LOGICAL, INTENT(IN), OPTIONAL            :: print_active
    TYPE(qs_rho_type), OPTIONAL, POINTER     :: rho_xc
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_ks_build_kohn_sham_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, nspins
    LOGICAL :: do_adiabatic_rescaling, do_ddapc, do_ep, do_hfx, do_ppl, &
      failure, gapw, gapw_xc, hfx_treat_lsd_in_core, just_energy_xc, &
      my_print, use_virial
    REAL(KIND=dp)                            :: ecore_ppl, edisp, ee_ener, &
                                                ekin_mol, mulliken_order_p
    REAL(KIND=dp), DIMENSION(3, 3)           :: h_stress
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(becke_restraint_type), POINTER      :: becke
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_h, matrix_p, &
                                                matrix_p_aux_fit, matrix_s, &
                                                mo_derivs, my_rho
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_aux_fit
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(ep_qs_type), POINTER                :: ep_qs_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array, mos_aux_fit
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type) :: rho_tot_gspace, v_efield_rspace, v_hartree_gspace, &
      v_sic_rspace, v_spin_ddapc_rest_r
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: v_rspace_new, &
                                                v_rspace_new_aux_fit, &
                                                v_tau_rspace, &
                                                v_tau_rspace_aux_fit
    TYPE(pw_p_type), POINTER                 :: vppl_rspace
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_rho_type), POINTER               :: my_rho_xc, rho_struct
    TYPE(section_vals_type), POINTER         :: adiabatic_rescaling_section, &
                                                hfx_sections, input, &
                                                scf_section, xc_section
    TYPE(virial_type), POINTER               :: virial

    failure=.FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(admm_env, becke, cell, dft_control, ep_qs_env, logger, &
         mo_array, mo_coeff, mo_coeff_aux_fit, mo_derivs, mos_aux_fit, my_rho, &
         my_rho_xc, rho_struct, para_env, pw_env, &
         virial, vppl_rspace)
    NULLIFY(adiabatic_rescaling_section,hfx_sections,input, scf_section,xc_section)
    NULLIFY(matrix_h,matrix_p,matrix_p_aux_fit,matrix_s, auxbas_pw_pool, &
         poisson_env, pw_pools, v_rspace_new, v_rspace_new_aux_fit, &
         v_tau_rspace,v_tau_rspace_aux_fit)

    CPPrecondition(ASSOCIATED(ks_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ks_env%ref_count>0,cp_failure_level,routineP,error,failure)

    logger => cp_error_get_logger(error)
    my_print = .TRUE.
    IF(PRESENT(print_active)) my_print = print_active

    CALL get_qs_env(qs_env=qs_env,&
         dft_control=dft_control,&
         matrix_h=matrix_h,&
         matrix_s=matrix_s,&
         pw_env=pw_env,&
         cell=cell,&
         para_env=para_env,&
         input=input,&
         virial=virial,&
         ep_qs_env=ep_qs_env,&
         error=error)

    use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)
    do_ep=ASSOCIATED(ep_qs_env)
    IF (do_ep) do_ep=ep_qs_env%ep_active.AND.calculate_forces.and.&
         ASSOCIATED(ep_qs_env%dH_coeffs)

    hfx_sections => section_vals_get_subs_vals(input,"DFT%XC%HF",error=error)
    CALL section_vals_get(hfx_sections,explicit=do_hfx,error=error)
    IF( do_hfx ) THEN
      CALL section_vals_val_get(hfx_sections, "TREAT_LSD_IN_CORE", l_val=hfx_treat_lsd_in_core,&
                                i_rep_section=1,error=error)
    END IF
    adiabatic_rescaling_section => section_vals_get_subs_vals(input,"DFT%XC%ADIABATIC_RESCALING",error=error)
    CALL section_vals_get(adiabatic_rescaling_section,explicit=do_adiabatic_rescaling,error=error)
    just_energy_xc=just_energy
    IF(do_adiabatic_rescaling) THEN
     !! If we perform adiabatic rescaling, the xc potential has to be scaled by the xc- and
     !! HFX-energy. Thus, let us first calculate the energy
     just_energy_xc=.TRUE.
    END IF

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
               qs_env%cp_ddapc_ewald%do_qmmm_periodic_decpl.OR.&
               qs_env%cp_ddapc_ewald%do_solvation

    ! Check for GAPW method : additional terms for local densities
    gapw    = dft_control%qs_control%gapw
    gapw_xc = dft_control%qs_control%gapw_xc
    IF(gapw_xc .AND. gapw) CALL cp_assert( .FALSE. , cp_failure_level,cp_assertion_failed,routineP,&
                           " GAPW and GAPW_XC are not compatible",error,failure)
!    IF((gapw_xc .OR. gapw) .AND. use_virial) &
!       CALL cp_unimplemented_error(fromWhere=routineP, &
!            message="Virial for GAPW not debugged!!", &
!            error=error, error_level=cp_warning_level)

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

       scf_section => section_vals_get_subs_vals(input,"DFT%SCF",error=error)
       IF (BTEST(cp_print_key_should_output(logger%iter_info,scf_section,&
            "PRINT%DETAILED_ENERGY",error=error),cp_p_file).AND.(.NOT.gapw).AND.(.NOT.gapw_xc)) THEN
          CALL pw_zero(rho_tot_gspace%pw, error=error)
          CALL calc_rho_tot_gspace(rho_tot_gspace, qs_env, rho, skip_nuclear_density=.TRUE.,error=error)
          CALL pw_poisson_solve(poisson_env,rho_tot_gspace%pw,energy%e_hartree,&
               v_hartree_gspace%pw,error=error)
          CALL pw_zero(rho_tot_gspace%pw, error=error)
          CALL pw_zero(v_hartree_gspace%pw, error=error)
       END IF

       ! Get the total density in g-space [ions + electrons]
       CALL calc_rho_tot_gspace(rho_tot_gspace, qs_env, rho,error=error)

       IF(my_print) THEN
         CALL print_densities(qs_env, rho, my_rho_xc, error=error)
       END IF

       ! Getting the Hartree energy and Hartree potential.  Also getting the stress tensor
       ! from the Hartree term if needed.  No nuclear force information here
       IF (use_virial .AND. calculate_forces) THEN
         h_stress(:,:) = 0.0_dp
         CALL pw_poisson_solve(poisson_env,rho_tot_gspace%pw,energy%hartree,&
                             v_hartree_gspace%pw,h_stress=h_stress,error=error)
         virial%pv_virial = virial%pv_virial + h_stress/REAL(para_env%num_pe,dp)
       ELSE
         CALL pw_poisson_solve(poisson_env,rho_tot_gspace%pw,energy%hartree,&
                             v_hartree_gspace%pw,error=error)
       END IF

       ! In case decouple periodic images and/or apply restraints to charges
       IF (do_ddapc) THEN
          CALL qs_ks_ddapc(qs_env, auxbas_pw_pool, rho_tot_gspace, v_hartree_gspace,&
               v_spin_ddapc_rest_r, energy, calculate_forces, ks_env, ks_matrix,    &
               just_energy, error)
       ELSE
          dft_control%qs_control%ddapc_explicit_potential      = .FALSE.
          dft_control%qs_control%ddapc_restraint_is_spin = .FALSE.
          IF (.NOT. just_energy) THEN
             CALL pw_transfer(v_hartree_gspace%pw, ks_env%v_hartree_rspace%pw, error=error)
             CALL pw_scale(ks_env%v_hartree_rspace%pw, ks_env%v_hartree_rspace%pw%pw_grid%dvol,&
                  error=error)
          END IF
       END IF
       CALL pw_pool_give_back_pw(auxbas_pw_pool,v_hartree_gspace%pw,error=error)

       IF(dft_control%apply_efield_field)THEN
          CALL pw_pool_create_pw(auxbas_pw_pool,v_efield_rspace%pw,&
            use_data=REALDATA3D, in_space=REALSPACE,error=error)
          CALL efield_potential(qs_env,v_efield_rspace,error)
          CALL pw_scale(v_efield_rspace%pw, v_efield_rspace%pw%pw_grid%dvol, error=error)
       END IF

       IF(dft_control%correct_surf_dip) THEN
         CALL calc_dipsurf_potential(qs_env,ks_env,energy,error=error)
         energy%hartree = energy%hartree + energy%surf_dipole
       END IF

       ! SIC
       CALL calc_v_sic_rspace(v_sic_rspace,energy,qs_env,dft_control,rho,poisson_env,&
            just_energy,calculate_forces,auxbas_pw_pool,error=error)

       IF (gapw) CALL Vh_1c_gg_integrals(qs_env,energy%hartree_1c,error=error)

       ! Check if becke potential is needed to constrain charges
       CALL qs_ks_becke_restraint(qs_env,auxbas_pw_pool,calculate_forces,matrix_s,becke,error)

       ! Adds the External Potential if requested
       IF(dft_control%apply_external_potential) THEN
          ! Compute the energy due to the external potential
          ee_ener=0.0_dp
          DO ispin=1,nspins
             ee_ener=ee_ener+accurate_sum(qs_env%vee%pw%cr3d*rho%rho_r(ispin)%pw%cr3d)*&
                  qs_env%vee%pw%pw_grid%dvol
          END DO
          CALL mp_sum(ee_ener,qs_env%para_env%group)
          ! the sign accounts for the charge of the electrons
          energy%ee=-ee_ener
       END IF

       ! Adds the QM/MM potential
       IF (qs_env%qmmm) THEN
          CALL qmmm_calculate_energy (qs_env=qs_env,&
                                      rho=rho%rho_r,&
                                      v_qmmm=qs_env%ks_qmmm_env%v_qmmm_rspace,&
                                      qmmm_energy=energy%qmmm_el,&
                                      error=error)
             IF (qs_env%qmmm_env_qm%image_charge) THEN
               CALL calculate_image_pot(v_hartree_rspace=ks_env%v_hartree_rspace,&
                                       rho_hartree_gspace=rho_tot_gspace,&
                                       energy=energy,&
                                       qmmm_env=qs_env%qmmm_env_qm,&
                                       qs_env=qs_env,error=error)
               IF (.NOT.just_energy) THEN
                  CALL add_image_pot_to_hartree_pot(v_hartree=ks_env%v_hartree_rspace,&
                                                 v_metal=qs_env%ks_qmmm_env%v_metal_rspace,&
                                                 qs_env=qs_env,error=error)
                  IF(calculate_forces) THEN
                    CALL  integrate_potential_devga_rspace(&
                                  potential=ks_env%v_hartree_rspace,coeff=qs_env%image_coeff,&
                                  forces=qs_env%qmmm_env_qm%image_charge_pot%image_forcesMM,&
                                  qmmm_env=qs_env%qmmm_env_qm,qs_env=qs_env,error=error)
                  ENDIF
               ENDIF
                  CALL pw_release(qs_env%ks_qmmm_env%v_metal_rspace%pw,error=error)
             END IF
          IF (.NOT.just_energy) THEN
             CALL qmmm_modify_hartree_pot(v_hartree=ks_env%v_hartree_rspace,&
                  v_qmmm=qs_env%ks_qmmm_env%v_qmmm_rspace,scale=1.0_dp,&
                  error=error)
          END IF
       END IF
       CALL pw_pool_give_back_pw(auxbas_pw_pool,rho_tot_gspace%pw,error=error)

       ! calculate the density matrix for the fitted mo_coeffs
       IF( dft_control%do_admm ) THEN
          IF(.NOT.qs_env%run_rtp)THEN
             CALL admm_density_matrix(qs_env,ks_env,dft_control,input,error=error)
          ELSE
             CALL rtp_admm_density_matrix(qs_env,ks_env,dft_control,input,error=error)
          END IF
       END IF

       ! only activate stress calculation if
       IF (use_virial .AND. calculate_forces) virial%pv_calculate = .TRUE.

       ! *** calculate the xc potential on the pw density ***
       ! *** associates v_rspace_new if the xc potential needs to be computed.
       ! If we do wavefunction fitting, we need the vxc_potential in the auxiliary basis set
       IF( dft_control%do_admm ) THEN
         CALL get_qs_env(qs_env, admm_env=admm_env,error=error)
         xc_section => admm_env%xc_section_aux
   
         IF (gapw_xc) THEN
           CALL get_qs_env(qs_env=qs_env, rho_xc=rho_struct, error=error)
         ELSE
           CALL get_qs_env(qs_env=qs_env, rho_aux_fit=rho_struct, error=error)
         END IF

         ! here we ignore a possible vdW section in admm_env%xc_section_aux
         CALL qs_vxc_create(qs_env=qs_env, rho_struct=rho_struct, xc_section=xc_section, &
              vxc_rho=v_rspace_new_aux_fit, vxc_tau=v_tau_rspace_aux_fit, exc=energy%exc_aux_fit, &
              just_energy=just_energy_xc, error=error)

         NULLIFY(rho_struct)

         IF (use_virial .AND. calculate_forces) THEN
           virial%pv_virial = virial%pv_virial - virial%pv_xc
           ! ** virial%pv_xc will be zeroed in the xc routines
         END IF
         xc_section => admm_env%xc_section_primary
       ELSE
         xc_section => section_vals_get_subs_vals(input,"DFT%XC",error=error)
       END IF

       IF (gapw_xc) THEN
         CALL get_qs_env(qs_env=qs_env, rho_xc=rho_struct, error=error) 
       ELSE
         CALL get_qs_env(qs_env=qs_env, rho=rho_struct, error=error) 
       END IF

       ! zmp 
       IF (qs_env%dft_control%apply_external_density .OR. qs_env%dft_control%apply_external_vxc) THEN
         energy%exc = 0.
         CALL calculate_zmp_potential(qs_env, v_rspace_new, rho, exc=energy%exc,error=error)
       ELSE
         CALL qs_vxc_create(qs_env=qs_env, rho_struct=rho_struct, xc_section=xc_section, &
              vxc_rho=v_rspace_new, vxc_tau=v_tau_rspace, exc=energy%exc, &
              edisp=edisp, dispersion_env=qs_env%dispersion_env, &
              just_energy=just_energy_xc, error=error)
         IF(edisp /= 0.0_dp) energy%dispersion = edisp
         IF (qs_env%requires_matrix_vxc) THEN
            CALL compute_matrix_vxc(qs_env=qs_env,v_rspace=v_rspace_new, matrix_vxc=qs_env%matrix_vxc, error=error)
         ENDIF

         IF (gapw .OR. gapw_xc) THEN
             CALL calculate_vxc_atom(qs_env,just_energy_xc,error)
         ENDIF
       ENDIF

       NULLIFY(rho_struct)
       IF (use_virial .AND. calculate_forces) THEN
          virial%pv_virial = virial%pv_virial - virial%pv_xc
       ENDIF

       ! *** Add Hartree-Fock contribution if required ***
       IF ( do_hfx ) THEN
         CALL hfx_ks_matrix(ks_env,qs_env,ks_matrix,rho,energy,calculate_forces,&
              just_energy,v_rspace_new,v_tau_rspace,error=error)
!!       Adiabatic rescaling  only if do_hfx; right?????
       END IF !do_hfx

       IF(do_ppl .AND. calculate_forces)THEN
          CPPrecondition(.NOT.gapw,cp_failure_level,routineP,error,failure)
          DO ispin=1,nspins
             CALL integrate_ppl_rspace(rho%rho_r(ispin),qs_env,error=error)
          END DO
       END IF

       IF(ASSOCIATED(qs_env%rho_nlcc) .AND. calculate_forces)THEN
          CPPrecondition(.NOT.(gapw .OR. gapw_xc),cp_failure_level,routineP,error,failure)
          DO ispin=1,nspins
             CALL integrate_rho_nlcc(v_rspace_new(ispin),qs_env,error=error)
             IF(dft_control%do_admm)CALL integrate_rho_nlcc(v_rspace_new_aux_fit(ispin),qs_env,error=error)
          END DO
       ENDIF

       ! calculate KG correction
       IF (qs_env%dft_control%qs_control%do_kg.AND.just_energy) THEN

         CALL kg_ekin_subset(qs_env, ks_matrix, gapw, gapw_xc, v_rspace_new, v_tau_rspace, ekin_mol, &
                             calculate_forces, error)

         ! substract kg corr from the total energy
         energy%exc = energy%exc - ekin_mol

       END IF  

       ! ***  Single atom contributions ***
       IF (.NOT. just_energy) THEN
          IF (calculate_forces ) THEN
             IF (do_ep) THEN
                ! Getting  nuclear force contribution from ep potential
                CALL ep_v_core(qs_env,ep_qs_env,nspins,error=error)
             ELSE
                ! Getting nuclear force contribution from the core charge density
                CALL integrate_v_core_rspace(ks_env%v_hartree_rspace, qs_env,error=error)
             END IF
          END IF

          IF (.NOT.do_hfx) THEN
             ! Initialize the Kohn-Sham matrix with the core Hamiltonian matrix
             DO ispin=1,nspins
                CALL cp_dbcsr_copy(ks_matrix(ispin)%matrix,matrix_h(1)%matrix,&
                                   name=cp_dbcsr_name(ks_matrix(ispin)%matrix),error=error)
             END DO
          END IF

          ! sum up potentials and integrate
          IF (do_ep) THEN
             my_rho => ep_qs_env%dH_coeffs
          ELSE
             ! Pointing my_rho to the density matrix rho_ao
             my_rho => rho%rho_ao
          END IF
          CALL  sum_up_and_integrate(ks_env,qs_env,ks_matrix,rho,my_rho,vppl_rspace,&
                v_rspace_new,v_rspace_new_aux_fit,v_tau_rspace,v_tau_rspace_aux_fit,&
                v_efield_rspace,v_sic_rspace,v_spin_ddapc_rest_r,&
                becke,calculate_forces,error=error)

          IF (qs_env%dft_control%qs_control%do_kg) THEN
       
            CALL kg_ekin_subset(qs_env, ks_matrix, gapw, gapw_xc, v_rspace_new, v_tau_rspace, ekin_mol, &
                                calculate_forces, error)

            ! substract kg corr from the total energy
            energy%exc = energy%exc - ekin_mol

          END IF

       END IF  ! .NOT. just energy

       IF (dft_control%qs_control%ddapc_explicit_potential) THEN
           CALL pw_pool_give_back_pw(auxbas_pw_pool,v_spin_ddapc_rest_r%pw,&
                                        error=error)
       ENDIF
       IF (calculate_forces.AND.qs_env%dft_control%qs_control%becke_restraint)&
           CALL pw_pool_give_back_pw(auxbas_pw_pool,becke%becke_pot%pw,&
                                        error=error)
       IF(dft_control%apply_efield_field)THEN
            CALL pw_pool_give_back_pw(auxbas_pw_pool,v_efield_rspace%pw,&
                                        error=error)
       END IF
       IF(gapw) THEN
          ! Integrals of the Hartree potential with g0_soft
          CALL integrate_vhg0_rspace(qs_env,ks_env%v_hartree_rspace, &
                                     calculate_forces, error=error)
       END IF

       IF(gapw .OR. gapw_xc) THEN
          ! Single atom contributions in the KS matrix ***
          CALL update_ks_atom(qs_env,ks_matrix,rho%rho_ao,calculate_forces,error=error)
       ENDIF

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

       ! At this point the ks matrix should be up to date, filter it if requested
       DO ispin=1,nspins
          IF(.FALSE.)WRITE(*,*) 'before',cp_dbcsr_get_occupation(ks_matrix(ispin)%matrix)
          CALL cp_dbcsr_filter(ks_matrix(ispin)%matrix,&
               qs_env%dft_control%qs_control%eps_filter_matrix,&
               error=error)
          IF(.FALSE.)WRITE(*,*) 'after',cp_dbcsr_get_occupation(ks_matrix(ispin)%matrix)
       ENDDO

       ! prepare to add electronic entropic contribution if needed
       IF (dft_control%smear) THEN
          energy%kTS = 0.0_dp
          energy%efermi = 0.0_dp
          CALL get_qs_env(qs_env,mos=mo_array,error=error)
          DO ispin=1,SIZE(mo_array)
             energy%kTS = energy%kTS + mo_array(ispin)%mo_set%kTS
             energy%efermi = energy%efermi + mo_array(ispin)%mo_set%mu
          ENDDO
          energy%efermi = energy%efermi /REAL(SIZE(mo_array),KIND=dp)
       ENDIF

       !** merge the auxiliary KS matrix and the primary one
       IF( dft_control%do_admm) THEN
         CALL get_qs_env(qs_env,mos=mo_array,error=error)
         DO ispin=1,nspins
           matrix_p =>qs_env%rho%rho_ao
           CALL get_qs_env(qs_env,mos_aux_fit=mos_aux_fit,error=error)
           CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,mo_coeff=mo_coeff)
           CALL get_mo_set(mo_set=mos_aux_fit(ispin)%mo_set,mo_coeff=mo_coeff_aux_fit)
           IF(.NOT.qs_env%run_rtp)THEN
              CALL admm_merge_ks_matrix(ispin, qs_env%admm_env, mo_array(ispin)%mo_set, mo_coeff,&
                   mo_coeff_aux_fit, qs_env%matrix_ks, qs_env%matrix_ks_aux_fit, qs_env%matrix_s,&
                   matrix_p_aux_fit, matrix_p, error)
           ELSE
              CALL rtp_admm_merge_ks_matrix(ispin,qs_env%rtp, qs_env%admm_env, mo_array(ispin)%mo_set, mo_coeff,&
                   mo_coeff_aux_fit, qs_env%matrix_ks, qs_env%matrix_ks_aux_fit, qs_env%matrix_s,&
                   matrix_p_aux_fit, matrix_p, error)

           END IF
         END DO
       END IF

       ! Right now we can compute the orbital derivative here, as it depends currently only on the available
       ! Kohn-Sham matrix. This might change in the future, in which case more pieces might need to be assembled
       ! from this routine, notice that this part of the calculation in not linear scaling
       ! right now this operation is only non-trivial because of occupation numbers and the restricted keyword
       IF (qs_env%requires_mo_derivs .AND. .NOT. just_energy.AND..NOT.qs_env%run_rtp) THEN
           CALL get_qs_env(qs_env,mo_derivs=mo_derivs,error=error)
           CALL calc_mo_derivatives(qs_env,ks_matrix,mo_derivs,error)
       ENDIF

       ! deal with low spin roks
        CALL low_spin_roks(energy,qs_env,dft_control,just_energy,&
                                  calculate_forces,auxbas_pw_pool,error)

       ! deal with sic on explicit orbitals
       CALL sic_explicit_orbitals(energy,qs_env,dft_control,poisson_env,just_energy,&
                                  calculate_forces,auxbas_pw_pool,error)
       CALL qs_efield_berry_phase(ks_env,qs_env,calculate_forces,error)
       ! adds s2_restraint energy and orbital derivatives
       CALL qs_ks_s2_restraint(dft_control,qs_env,matrix_s,&
            mo_array,energy,calculate_forces, just_energy,error)

       IF ( do_ppl ) THEN
         ! update core energy for grid based local pseudopotential
         ecore_ppl = 0._dp
         DO ispin=1,nspins
           ecore_ppl = ecore_ppl+&
              SUM(vppl_rspace%pw%cr3d*rho%rho_r(ispin)%pw%cr3d)*vppl_rspace%pw%pw_grid%dvol
         END DO
         CALL mp_sum(ecore_ppl,para_env%group)
         energy%core = energy%core+ecore_ppl
       END IF

       ! sum all energy terms to obtain the total energy
       energy%total = energy%core_overlap + energy%core_self + &
            energy%core + energy%hartree + &
            energy%hartree_1c + &
            energy%exc + energy%exc1 + energy%ex + energy%dispersion +&
            energy%qmmm_el + energy%mulliken + &
            SUM(energy%ddapc_restraint) + energy%s2_restraint +&
            energy%becke + energy%dft_plus_u + energy%kTS + &
            energy%efield +energy%efield_core+&
            energy%ee + energy%ee_core+&
            energy%exc_aux_fit+&
            energy%image_charge
       
       !Print detailed energy
       IF(my_print) THEN
         CALL print_detailed_energy(qs_env,dft_control,input,energy,mulliken_order_p,error=error)
       END IF ! my_print

    END IF ! failure

    CALL timestop(handle)

  END SUBROUTINE qs_ks_build_kohn_sham_matrix

! *****************************************************************************
!> \brief calculation of the LRIGPW KS matrix
!> \param ks_env the ks_env that holds all the temporary objects that
!>        the calculation of the LRI KS matrix needs
!> \param qs_env the qs_env to update
!> \param ks the KS matrix that is created
!> \param rho the density, must be up to date
!> \param energy the place where energies are stored
!> \param just_energy if true updates the energies but not the
!>        ks matrix. Defaults to false
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note  might be integrate in qs_ks_build_kohn_sham_matrix later
!> \par History
!>      01.2014 created
!> \author Dorothea Golze
! *****************************************************************************
  SUBROUTINE lri_build_ks_matrix(ks_env,qs_env,ks_matrix,rho,energy,&
             calculate_forces,just_energy,print_active,error)

    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(qs_energy_type), POINTER            :: energy
    LOGICAL, INTENT(in)                      :: calculate_forces, just_energy
    LOGICAL, INTENT(in), OPTIONAL            :: print_active
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'lri_build_ks_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, nspins, stat
    LOGICAL                                  :: failure, just_energy_xc, &
                                                my_print
    REAL(KIND=dp)                            :: edisp
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_h, matrix_p, mo_derivs
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(lri_density_type), POINTER          :: lri_density
    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(lri_kind_type), DIMENSION(:), &
      POINTER                                :: lri_v_int
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: rho_tot_gspace, &
                                                v_hartree_gspace
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: v_rspace_new, v_tau_rspace
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_rho_type), POINTER               :: my_rho_xc
    TYPE(section_vals_type), POINTER         :: input, xc_section

    failure=.FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(dft_control, input,  mo_derivs, my_rho_xc, para_env,&
            pw_env, xc_section)
    NULLIFY(auxbas_pw_pool, matrix_h, matrix_p, poisson_env, pw_pools, v_rspace_new, &
            v_tau_rspace)
    NULLIFY(atomic_kind_set, lri_env, lri_density, lri_v_int)

    CPPrecondition(ASSOCIATED(ks_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ks_env%ref_count>0,cp_failure_level,routineP,error,failure)

    my_print = .TRUE.
    IF(PRESENT(print_active)) my_print = print_active

    CALL get_qs_env(qs_env=qs_env,&
         atomic_kind_set=atomic_kind_set,&
         dft_control=dft_control,&
         matrix_h=matrix_h,&
         lri_env=lri_env,&
         lri_density=lri_density,&
         pw_env=pw_env,&
         para_env=para_env,&
         input=input,&
         error=error)

    nspins=dft_control%nspins
    CPPrecondition(SIZE(ks_matrix)>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(pw_env),cp_failure_level,routineP,error,failure)

    ! gets the tmp grids
    IF (.NOT. failure) THEN
       CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,&
                       pw_pools=pw_pools, poisson_env=poisson_env,error=error)
    END IF

    IF (.NOT. failure) THEN

       ! *** calculate the hartree potential on the LRI pw density ***
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

       CALL calc_rho_tot_gspace(rho_tot_gspace, qs_env, rho, error=error)

       IF(my_print) THEN
         CALL print_densities(qs_env, rho, my_rho_xc, error=error)
       END IF

       CALL pw_poisson_solve(poisson_env,rho_tot_gspace%pw,energy%hartree,&
                             v_hartree_gspace%pw,error=error)

       IF (.NOT. just_energy) THEN
          CALL pw_transfer(v_hartree_gspace%pw, ks_env%v_hartree_rspace%pw, error=error)
          CALL pw_scale(ks_env%v_hartree_rspace%pw, ks_env%v_hartree_rspace%pw%pw_grid%dvol,&
               error=error)
       END IF
       CALL pw_pool_give_back_pw(auxbas_pw_pool,v_hartree_gspace%pw,error=error)
       CALL pw_pool_give_back_pw(auxbas_pw_pool,rho_tot_gspace%pw,error=error)

       ! core contribution
       IF (.NOT. just_energy) THEN
             ! Initialize the Kohn-Sham matrix with the core Hamiltonian matrix
             DO ispin=1,nspins
                CALL cp_dbcsr_copy(ks_matrix(ispin)%matrix,matrix_h(1)%matrix,&
                                   name=cp_dbcsr_name(ks_matrix(ispin)%matrix),error=error)
             END DO
       END IF

      ! exchange-correlation
      xc_section => section_vals_get_subs_vals(input,"DFT%XC",error=error)
      just_energy_xc=just_energy

      CALL qs_vxc_create(qs_env=qs_env, rho_struct=rho, xc_section=xc_section, &
           vxc_rho=v_rspace_new, vxc_tau=v_tau_rspace, exc=energy%exc, &
           edisp=edisp, dispersion_env=qs_env%dispersion_env, &
           just_energy=just_energy_xc, error=error)
      IF(edisp /= 0.0_dp) energy%dispersion = edisp

      ! sum v_hartree + v_xc
      IF (ASSOCIATED(v_rspace_new)) THEN
        DO ispin=1,nspins
         v_rspace_new(ispin)%pw%cr3d  =&
                  v_rspace_new(ispin)%pw%pw_grid%dvol * v_rspace_new(ispin)%pw%cr3d + &
                  ks_env%v_hartree_rspace%pw%cr3d 
        ENDDO
      ENDIF
      ! NO TAU YET
 
      !do integration and build matrix
      IF (.NOT. just_energy) THEN

       DO ispin=1,nspins
          lri_v_int => lri_density%lri_coefs(ispin)%lri_kinds
          CALL integrate_v_rspace_one_center(v_rspace_new(ispin),qs_env,&
               lri_v_int,calculate_forces,error)
          CALL pw_pool_give_back_pw(auxbas_pw_pool,v_rspace_new(ispin)%pw,&
               error=error)
          CALL calculate_lri_ks_matrix(lri_env, lri_v_int, ks_matrix(ispin)%matrix,&
               atomic_kind_set, para_env, error)
       ENDDO
       DEALLOCATE(v_rspace_new,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

       IF (calculate_forces ) THEN
           ! Getting nuclear force contribution from the core charge density
            CALL integrate_v_core_rspace(ks_env%v_hartree_rspace, qs_env,error=error)
           ! Getting electronic force contribution from fitted LRI charge density
            matrix_p =>qs_env%rho%rho_ao
            CALL calculate_lri_forces(lri_env,lri_density,qs_env,matrix_p,atomic_kind_set,error)
       END IF

      ENDIF
      
      ! Right now we can compute the orbital derivative here, as it depends currently only on the available
      ! Kohn-Sham matrix. This might change in the future, in which case more pieces might need to be assembled
      ! from this routine, notice that this part of the calculation in not linear scaling
      ! right now this operation is only non-trivial because of occupation numbers and the restricted keyword
      IF (qs_env%requires_mo_derivs .AND. .NOT. just_energy.AND..NOT.qs_env%run_rtp) THEN
         CALL get_qs_env(qs_env,mo_derivs=mo_derivs,error=error)
         CALL calc_mo_derivatives(qs_env,ks_matrix,mo_derivs,error)
      ENDIF

      energy%total = energy%core_overlap + energy%core_self + &
            energy%core + energy%hartree + &
            energy%hartree_1c + &
            energy%exc + energy%exc1 + energy%ex + energy%dispersion +&
            energy%qmmm_el + energy%mulliken + &
            SUM(energy%ddapc_restraint) + energy%s2_restraint +&
            energy%becke + energy%dft_plus_u + energy%kTS + &
            energy%efield +energy%efield_core+&
            energy%ee + energy%ee_core+&
            energy%exc_aux_fit+&
            energy%image_charge


    END IF ! failure

    CALL timestop(handle)

  END SUBROUTINE lri_build_ks_matrix

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
!> \brief compute MO derivatives
!> \param ks_env the ks_env that old all the temporary objects that
!>        the calculation of the KS matrix needs
!> \param qs_env the qs_env to update
!> \par History
!>      01.2014 created, transferred from qs_ks_build_kohn_sham_matrix in
!>      separate subroutine
!> \author Dorothea Golze
! *****************************************************************************
  SUBROUTINE calc_mo_derivatives(qs_env,ks_matrix,mo_derivs,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix, mo_derivs
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_mo_derivatives', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin
    LOGICAL                                  :: failure, uniform_occupation
    REAL(KIND=dp), DIMENSION(:), POINTER     :: occupation_numbers
    TYPE(cp_dbcsr_type), POINTER             :: mo_coeff_b
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mo_derivs_aux_fit, &
                                                mo_derivs_tmp
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_aux_fit
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array, mos_aux_fit

    failure=.FALSE.
    NULLIFY(dft_control, mo_array, mo_coeff, mo_coeff_aux_fit, mo_coeff_b, &
            mo_derivs_aux_fit, mo_derivs_tmp, mos_aux_fit,&
            occupation_numbers)


    CALL get_qs_env(qs_env,dft_control=dft_control,mos=mo_array,error=error)

    IF(dft_control%do_admm) THEN !fm->dbcsr
       NULLIFY(mo_derivs_tmp)!fm->dbcsr
       ALLOCATE(mo_derivs_tmp(SIZE(mo_derivs)))
       DO ispin=1,SIZE(mo_derivs)!fm->dbcsr
          CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,mo_coeff=mo_coeff)!fm->dbcsr
          NULLIFY(mo_derivs_tmp(ispin)%matrix)
          CALL cp_fm_create(mo_derivs_tmp(ispin)%matrix,mo_coeff%matrix_struct,error=error)!fm->dbcsr
       ENDDO!fm->dbcsr
    ENDIF!fm->dbcsr


    DO ispin=1,SIZE(mo_derivs)

       CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,mo_coeff=mo_coeff,&
            mo_coeff_b=mo_coeff_b, occupation_numbers=occupation_numbers )
       CALL cp_dbcsr_multiply('n','n',1.0_dp,ks_matrix(ispin)%matrix,mo_coeff_b,&
             0.0_dp,mo_derivs(ispin)%matrix, error=error)

       IF( dft_control%do_admm) THEN
         CALL get_qs_env(qs_env, mos_aux_fit=mos_aux_fit,&
              mo_derivs_aux_fit=mo_derivs_aux_fit,error=error)
         CALL get_mo_set(mo_set=mos_aux_fit(ispin)%mo_set,mo_coeff=mo_coeff_aux_fit)

         CALL copy_dbcsr_to_fm(mo_derivs(ispin)%matrix,&
              mo_derivs_tmp(ispin)%matrix,error=error)!fm->dbcsr
         CALL admm_merge_mo_derivs(ispin, qs_env%admm_env, mo_array(ispin)%mo_set, &
              mo_coeff, mo_coeff_aux_fit, mo_derivs_tmp, mo_derivs_aux_fit, &
              qs_env%matrix_ks_aux_fit, error=error)
         CALL copy_fm_to_dbcsr(mo_derivs_tmp(ispin)%matrix,mo_derivs(ispin)%matrix,error=error)!fm->dbcsr
       END IF

       IF (dft_control%restricted) THEN
           ! only the first mo_set are actual variables, but we still need both
           CPPrecondition(ispin==1, cp_failure_level, routineP, error, failure)
           CPPrecondition(SIZE(mo_array)==2, cp_failure_level, routineP, error, failure)
           ! use a temporary array with the same size as the first spin for the second spin

           ! uniform_occupation is needed for this case, otherwise we can no
           ! reconstruct things in ot, since we irreversibly sum
           CALL get_mo_set(mo_set=mo_array(1)%mo_set, uniform_occupation=uniform_occupation)
           CPPrecondition(uniform_occupation, cp_failure_level, routineP, error, failure)
           CALL get_mo_set(mo_set=mo_array(2)%mo_set,&
                uniform_occupation=uniform_occupation)
           CPPrecondition(uniform_occupation, cp_failure_level, routineP, error, failure)

           CALL get_mo_set(mo_set=mo_array(2)%mo_set,mo_coeff_b=mo_coeff_b)
           CALL cp_dbcsr_multiply('n','n',1.0_dp,ks_matrix(2)%matrix,mo_coeff_b,&
                      1.0_dp, mo_derivs(1)%matrix, error=error)
       ENDIF

    ENDDO

    IF(dft_control%do_admm) THEN !fm->dbcsr
       DO ispin=1,SIZE(mo_derivs)!fm->dbcsr
          CALL cp_fm_release(mo_derivs_tmp(ispin)%matrix,error=error)!fm->dbcsr
       ENDDO!fm->dbcsr
       DEALLOCATE(mo_derivs_tmp)!fm->dbcsr
    ENDIF!fm->dbcsr

  END SUBROUTINE calc_mo_derivatives

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
  SUBROUTINE qs_ks_update_qs_env(ks_env,qs_env,calculate_forces,&
       just_energy,print_active,error)
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN), OPTIONAL            :: calculate_forces, &
                                                just_energy, print_active
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_ks_update_qs_env', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin
    LOGICAL                                  :: c_forces, did_update, &
                                                energy_only, failure
    REAL(KIND=dp)                            :: ecore, ekin
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_h, matrix_ks, matrix_t
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_rho_type), POINTER               :: rho, rho_xc

    c_forces=.FALSE.
    energy_only = .FALSE.
    IF (PRESENT(just_energy)) energy_only=just_energy
    IF (PRESENT(calculate_forces)) c_forces=calculate_forces

    IF (c_forces) THEN
       CALL timeset(routineN//'_forces',handle)
    ELSE
       CALL timeset(routineN,handle)
    ENDIF

    failure=.FALSE.
    did_update=.FALSE.
    NULLIFY(matrix_ks,matrix_t,matrix_h,energy, logger, rho, rho_xc, dft_control, para_env)
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
               dft_control=dft_control,kinetic=matrix_t, &
               matrix_h=matrix_h,para_env=para_env,error=error)
          ! update ecore
          energy%core=0.0_dp
          DO ispin=1,dft_control%nspins
             CALL calculate_ecore(h=matrix_h(1)%matrix,&
                  p=rho%rho_ao(ispin)%matrix,&
                  ecore=ecore,&
                  error=error)
             energy%core=energy%core+ecore
          END DO

          ! update ekin, but only when doing scptb, otherwise segfault.
          ! who knows, it might segfault anyway, this is a kindof ugly solution
!          IF (dft_control%qs_control%scptb) THEN
          IF(ASSOCIATED(matrix_t)) THEN
              energy%kinetic=0.0_dp
              DO ispin=1,dft_control%nspins
                 CALL calculate_ecore(h=matrix_t(1)%matrix,&
                      p=rho%rho_ao(ispin)%matrix,&
                      ecore=ekin,&
                      error=error)
                 !abusing the ecore method to calculate the same trace.
                 energy%kinetic=energy%kinetic+ekin
              END DO
          END IF

          ! the ks matrix will be rebuilt so this is fine now
          ks_env%potential_changed=.FALSE.
          IF ( dft_control%qs_control%semi_empirical ) THEN
             CALL build_se_fock_matrix_control(ks_env,qs_env=qs_env,&
                  ks_matrix=matrix_ks,rho=rho,energy=energy,&
                  calculate_forces=c_forces,just_energy=energy_only,&
                  error=error)
          ELSEIF ( dft_control%qs_control%dftb ) THEN
             CALL build_dftb_ks_matrix(ks_env,qs_env=qs_env,&
                  ks_matrix=matrix_ks,rho=rho,energy=energy,&
                  calculate_forces=c_forces,just_energy=energy_only,&
                  error=error)
          ELSEIF ( dft_control%qs_control%scptb ) THEN
             CALL build_scptb_ks_matrix(ks_env,qs_env=qs_env,&
                  ks_matrix=matrix_ks,rho=rho,energy=energy,&
                  calculate_forces=c_forces,just_energy=energy_only,&
                  error=error)
          ELSEIF (dft_control%qs_control%gapw_xc) THEN
             CALL qs_ks_build_kohn_sham_matrix(ks_env,qs_env=qs_env,&
                  ks_matrix=matrix_ks, rho=rho, energy=energy,&
                  calculate_forces=c_forces,just_energy=energy_only,&
                  print_active=print_active,rho_xc=rho_xc,error=error)
          ELSEIF (dft_control%qs_control%lrigpw) THEN
            CALL lri_build_ks_matrix(ks_env,qs_env,ks_matrix=matrix_ks,& 
                 rho=rho,energy=energy,calculate_forces=c_forces, &
                 just_energy=energy_only,print_active=print_active,&
                 error=error)
          ELSE
             CALL qs_ks_build_kohn_sham_matrix(ks_env,qs_env=qs_env,&
                  ks_matrix=matrix_ks, rho=rho, energy=energy,&
                  calculate_forces=c_forces,just_energy=energy_only,&
                  print_active=print_active,error=error)
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
    TYPE(cp_dbcsr_type), POINTER             :: w_matrix
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

    CALL cp_dbcsr_set(w_matrix,0.0_dp,error=error)
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

    CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=w_matrix,&
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
  SUBROUTINE calculate_w_matrix_ot(mo_set,mo_deriv,w_matrix,s_matrix,error)

    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_dbcsr_type), POINTER             :: mo_deriv, w_matrix, s_matrix
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_w_matrix_ot', &
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
    CALL copy_dbcsr_to_fm(mo_deriv,weighted_vectors,error=error)
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

    CALL cp_dbcsr_set(w_matrix,0.0_dp,error=error)
    CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=w_matrix,&
                            matrix_v=mo_set%mo_coeff,&
                            matrix_g=weighted_vectors,&
                            ncol=mo_set%homo,error=error)

    IF (check_gradient) THEN
       CALL cp_fm_create(gradient,mo_set%mo_coeff%matrix_struct,"gradient",error=error)
       CALL cp_dbcsr_sm_fm_multiply(s_matrix,weighted_vectors,&
                                 gradient, ncol_global,error=error)

       ALLOCATE(scaling_factor(SIZE(occupation_numbers)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       scaling_factor=2.0_dp*occupation_numbers
       CALL copy_dbcsr_to_fm(mo_deriv,weighted_vectors,error=error)
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

  END SUBROUTINE calculate_w_matrix_ot

! *****************************************************************************
!> \brief Calculate the energy-weighted density matrix W if ROKS is active.
!>        The W matrix is returned in matrix_w.
!> \author 04.05.06,MK
! *****************************************************************************
  SUBROUTINE calculate_w_matrix_roks(mo_set,matrix_ks,matrix_p,matrix_w,error)
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_dbcsr_type), POINTER             :: matrix_ks, matrix_p, matrix_w
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_w_matrix_roks', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nao
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct
    TYPE(cp_fm_type), POINTER                :: c, ks, p, work
    TYPE(cp_para_env_type), POINTER          :: para_env

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
    CALL copy_dbcsr_to_fm(matrix_ks,ks,error=error)
    CALL copy_dbcsr_to_fm(matrix_p,p,error=error)
    CALL cp_fm_upper_to_full(p,work,error)
    CALL cp_fm_symm("L","U",nao,nao,1.0_dp,ks,p,0.0_dp,work,error=error)
    CALL cp_fm_gemm("T","N",nao,nao,nao,1.0_dp,p,work,0.0_dp,ks,error=error)
    CALL cp_dbcsr_set(matrix_w,0.0_dp,error=error)
    CALL copy_fm_to_dbcsr(ks,matrix_w,keep_sparsity=.TRUE., error=error)

    CALL cp_fm_release(work,error=error)
    CALL cp_fm_release(p,error=error)
    CALL cp_fm_release(ks,error=error)

    CALL timestop(handle)

  END SUBROUTINE calculate_w_matrix_roks

! *****************************************************************************
!> \brief Allocate ks_matrix and ks_env if necessary
!> \par History
!>    refactoring 04.03.2011 [MI]
!> \author
! *****************************************************************************

  SUBROUTINE qs_ks_allocate_basics(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_ks_allocate_basics', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: headline
    INTEGER                                  :: ispin
    TYPE(array_i1d_obj)                      :: col_blk_size
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_ks_aux_fit, &
                                                matrix_s
    TYPE(dbcsr_distribution_obj)             :: dist
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_ks_env_type), POINTER            :: ks_env

    NULLIFY(dft_control, matrix_s, matrix_ks_aux_fit, matrix_ks, ks_env)

    CALL get_qs_env(qs_env=qs_env,&
         dft_control=dft_control,&
         matrix_ks=matrix_ks,&
         matrix_ks_aux_fit=matrix_ks_aux_fit,&
         matrix_s=matrix_s,&
         ks_env=ks_env,&
         error=error)

    IF (.NOT.ASSOCIATED(matrix_ks)) THEN
       CALL cp_dbcsr_allocate_matrix_set(matrix_ks, dft_control%nspins, error)
       DO ispin=1,dft_control%nspins
          ALLOCATE(matrix_ks(ispin)%matrix)
          CALL cp_dbcsr_init(matrix_ks(ispin)%matrix,error=error)
          IF (dft_control%nspins > 1) THEN
             IF (ispin == 1) THEN
                headline="KOHN-SHAM MATRIX FOR ALPHA SPIN"
             ELSE
                headline="KOHN-SHAM MATRIX FOR BETA SPIN"
             END IF
          ELSE
             headline="KOHN-SHAM MATRIX"
          END IF
          CALL cp_dbcsr_create(matrix=matrix_ks(ispin)%matrix,&
               name=TRIM(headline),&
               dist=cp_dbcsr_distribution(matrix_s(1)%matrix), matrix_type=dbcsr_type_symmetric,&
               row_blk_size=cp_dbcsr_row_block_sizes(matrix_s(1)%matrix),&
               col_blk_size=cp_dbcsr_col_block_sizes(matrix_s(1)%matrix),&
               nblks=0, nze=0, error=error)
          CALL cp_dbcsr_alloc_block_from_nbl(matrix_ks(ispin)%matrix,qs_env%sab_orb,error=error)
          CALL cp_dbcsr_set(matrix_ks(ispin)%matrix,0.0_dp,error=error)
       ENDDO

       CALL set_qs_env(qs_env=qs_env,&
            matrix_ks=matrix_ks,error=error)

    END IF

!   *** Allocate matrix_ks_aux_fit if requested and put it in the QS environment ***

    IF( dft_control%do_admm) THEN
      IF (.NOT.ASSOCIATED(matrix_ks_aux_fit)) THEN
         CALL cp_dbcsr_allocate_matrix_set(matrix_ks_aux_fit, dft_control%nspins, error)
         DO ispin=1,dft_control%nspins
            ALLOCATE(matrix_ks_aux_fit(ispin)%matrix)
            CALL cp_dbcsr_init(matrix_ks_aux_fit(ispin)%matrix,error=error)

            CALL cp_dbcsr_create(matrix=matrix_ks_aux_fit(ispin)%matrix, &
                 name="KOHN-SHAM_MATRIX for ADMM", &
                 dist=cp_dbcsr_distribution(qs_env%matrix_s_aux_fit(1)%matrix), matrix_type=dbcsr_type_symmetric, &
                 row_blk_size=cp_dbcsr_row_block_sizes(qs_env%matrix_s_aux_fit(1)%matrix),&
                 col_blk_size=cp_dbcsr_col_block_sizes(qs_env%matrix_s_aux_fit(1)%matrix),&
                 nblks=0, nze=0, error=error)
            CALL cp_dbcsr_alloc_block_from_nbl(matrix_ks_aux_fit(ispin)%matrix,qs_env%sab_aux_fit,&
                 error=error)
            CALL cp_dbcsr_set(matrix_ks_aux_fit(ispin)%matrix,0.0_dp,error=error)
         ENDDO

         CALL set_qs_env(qs_env=qs_env,&
              matrix_ks_aux_fit=matrix_ks_aux_fit,error=error)
      END IF
    END IF

!   *** allocate the ks env **
    IF (.NOT.ASSOCIATED(ks_env)) THEN
       CALL qs_ks_create(ks_env,qs_env=qs_env,error=error)
       CALL set_qs_env(qs_env, ks_env=ks_env,error=error)
       CALL qs_ks_release(ks_env,error=error)
    END IF


  END SUBROUTINE qs_ks_allocate_basics

END MODULE qs_ks_methods
