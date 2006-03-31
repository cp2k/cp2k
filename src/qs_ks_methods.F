!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002,2003,2004  CP2K developers group                       !
!-----------------------------------------------------------------------------!


!!****h* cp2k/qs_ks_methods [1.0] *
!!
!!   NAME
!!     qs_ks_methods
!!
!!   FUNCTION
!!     routines that build the Kohn-Sham matrix (i.e calculate the coulomb
!!     and xc parts
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     05.2002 moved from qs_scf (see there the history) [fawzi]
!!     JGH [30.08.02] multi-grid arrays independent from density and potential
!!     10.2002 introduced pools, uses updated rho as input,
!!             removed most temporary variables, renamed may vars,
!!             began conversion to LSD [fawzi]
!!     10.2004 moved calculate_w_matrix here [Joost VandeVondele]
!!             introduced energy derivative wrt MOs [Joost VandeVondele]
!!
!!   SOURCE
!****************************************************************************
MODULE qs_ks_methods
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE cell_types,                      ONLY: cell_type
  USE coefficient_types,               ONLY: &
       coeff_add, coeff_copy, coeff_init, coeff_scale, coeff_subtract, &
       coeff_sumup, coeff_transform_space, coeff_type, coeff_zero
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_ddapc,                        ONLY: cp_ddapc_apply
  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                             cp_fm_gemm,&
                                             cp_fm_scale_and_add,&
                                             cp_fm_transpose
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
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply,&
                                             cp_sm_plus_fm_fm_t
  USE harris_env_types,                ONLY: harris_env_type
  USE hartree_local_2centers,          ONLY: integrate_vhg0_g0
  USE hartree_local_3centers,          ONLY: integrate_a_vhg0_b
  USE hartree_local_methods,           ONLY: Vh_1c_gg_integrals
  USE input_constants,                 ONLY: sic_ad,&
                                             sic_mauri_spz,&
                                             sic_mauri_us,&
                                             sic_none
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kahan_sum,                       ONLY: accurate_sum
  USE kg_environment_types,            ONLY: get_kg_env,&
                                             kg_environment_type
  USE kg_gpw_correction,               ONLY: kg_gpw_ekin_mol
  USE kg_gpw_fm_mol_types,             ONLY: kg_fm_p_type
  USE kg_gpw_pw_env_types,             ONLY: kg_sub_pw_env_type
  USE kinds,                           ONLY: dp
  USE mulliken,                        ONLY: mulliken_restraint
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_poisson_methods,              ONLY: pw_poisson_solve
  USE pw_poisson_types,                ONLY: pw_poisson_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_coeff,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_init_coeff,&
                                             pw_pool_p_type,&
                                             pw_pool_retain,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type,&
                                             pw_release,&
                                             pw_transfer,&
                                             pw_type
  USE qs_charges_types,                ONLY: qs_charges_type
  USE qs_collocate_density,            ONLY: calculate_total_rho
  USE qs_core_energies,                ONLY: calculate_ecore
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_gapw_densities,               ONLY: prepare_gapw_den
  USE qs_integrate_potential,          ONLY: integrate_v_core_rspace,&
                                             integrate_v_rspace
  USE qs_ks_atom,                      ONLY: update_ks_atom
  USE qs_ks_qmmm_methods,              ONLY: qmmm_calculate_energy,&
                                             qmmm_modify_hartree_pot
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             mo_set_type
  USE qs_rho0_ggrid,                   ONLY: integrate_vhg0_rspace
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE qs_semi_empirical_hamiltonian,   ONLY: build_se_fock_matrix
  USE qs_vxc_atom,                     ONLY: calculate_vxc_atom
  USE s_square_methods,                ONLY: s2_restraint
  USE sparse_matrix_types,             ONLY: real_matrix_p_type,&
                                             real_matrix_type,&
                                             set_matrix,&
                                             transfer_matrix
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc,                              ONLY: xc_exc_calc,&
                                             xc_vxc_pw_create1
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  INTERFACE calculate_w_matrix
    MODULE PROCEDURE calculate_w_matrix_1,&
                     calculate_w_matrix_2
  END INTERFACE

  LOGICAL, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ks_methods'
  INTEGER, SAVE, PRIVATE :: last_ks_id_nr=0

  PUBLIC :: qs_ks_create, qs_ks_did_change, qs_vxc_create, calc_rho_tot_gspace, &
       qs_ks_update_qs_env, qs_ks_build_kohn_sham_matrix, calculate_w_matrix
!***
!****************************************************************************

CONTAINS

!!****f* qs_ks_methods/qs_ks_create [1.0] *
!!
!!   NAME
!!     qs_ks_create
!!
!!   FUNCTION
!!     allocates and initializes the given ks_env.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - ks_env: the ks env to be initialized
!!     - qs_env: the qs environment
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
  SUBROUTINE qs_ks_create(ks_env, qs_env, error)
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.
    CPPreconditionNoFail(.NOT.ASSOCIATED(ks_env),cp_failure_level,routineP,error)
    ALLOCATE(ks_env, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    IF (.not.failure) THEN
       NULLIFY( pw_env, auxbas_pw_pool,&
            pw_env, cell,pw_pools)
       CALL get_qs_env(qs_env=qs_env,&
            dft_control=dft_control,&
            pw_env=pw_env, cell=cell)
       CALL pw_env_get(pw_env,auxbas_pw_pool=auxbas_pw_pool,pw_pools=pw_pools,&
            auxbas_grid=auxbas_grid)
       
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
       CALL pw_pool_retain(ks_env%auxbas_pw_pool)
       CALL pw_pool_init_coeff(auxbas_pw_pool,ks_env%v_hartree_rspace,&
            use_data=REALDATA3D, in_space=REALSPACE)
    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE qs_ks_create
!***************************************************************************

!!****f* qs_ks_methods/qs_ks_build_kohn_sham_matrix [1.0] *
!!
!!   NAME
!!     qs_ks_build_kohn_sham_matrix
!!
!!   FUNCTION
!!     routine where the real calculations are made: the
!!     KS matrix is calculated
!!
!!   NOTES
!!     make rho, energy and qs_charges optional, defaulting
!!     to qs_env components?
!!
!!   INPUTS
!!     - ks_env: the ks_env that old all the temporary objects that
!!       the calculation of the KS matrix needs
!!     - qs_env: the qs_env to update
!!     - calculate_forces: if true calculate the quantities needed
!!       to calculate the forces. Defaults to false.
!!     - just_energy: if true updates the energies but not the
!!       ks matrix. Defaults to false
!!     - ks: the KS matrix that is created
!!     - rho: the density, must be up to date
!!     - rho_xc: the soft density, only for gapw_xc , must be up to date
!!     - energy: the place where energies are stored
!!     - qs_charges: the place where the grids charges are stored
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Matthias Krack
!!
!!   MODIFICATION HISTORY
!!     06.2002 moved from qs_scf to qs_ks_methods, use of ks_env
!!             new did_change scheme [fawzi]
!!     10.2002 introduced pools, uses updated rho as input, LSD [fawzi]
!!     10.2004 build_kohn_sham matrix now also computes the derivatives
!!             of the total energy wrt to the MO coefs, if instructed to
!!             do so. This appears useful for orbital dependent functionals
!!             where the KS matrix alone (however this might be defined)
!!              does not contain the info to construct this derivative.
!!
!!*** *********************************************************************
  SUBROUTINE qs_ks_build_kohn_sham_matrix(ks_env,qs_env,ks_matrix,&
       rho,energy,calculate_forces,just_energy,rho_xc,kg_gpw,kg_fm_set,kg_sub_pw_env, error)

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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ks_build_kohn_sham_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, k, k2, n, &
                                                nspins, output_unit, stat
    LOGICAL                                  :: failure, gapw, gapw_xc, &
                                                ionode, my_kg_gpw
    REAL(KIND=dp)                            :: ddapc_order_p, &
                                                mulliken_order_p, s2_order_p
    REAL(KIND=dp), DIMENSION(:), POINTER     :: occupation_numbers, &
                                                scaling_factor
    TYPE(cell_type), POINTER                 :: cell
    TYPE(coeff_type)                         :: rho_tot_gspace, &
                                                v_hartree_gspace, v_sic_rspace
    TYPE(coeff_type), DIMENSION(:), POINTER  :: v_rspace_new, v_tau_rspace
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mo_derivs
    TYPE(cp_fm_type), POINTER                :: hc_tmp, mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_rho_type), POINTER               :: my_rho_xc
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_h, matrix_s
    TYPE(section_vals_type), POINTER         :: input, scf_section

! shorter for output

    failure=.FALSE.
    NULLIFY(matrix_h,matrix_s,dft_control,auxbas_pw_pool, pw_pools,&
         pw_env, cell, logger, v_rspace_new, v_tau_rspace,&
         input, poisson_env, scf_section)

    !   -----------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

    CPPrecondition(ASSOCIATED(ks_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ks_env%ref_count>0,cp_failure_level,routineP,error,failure)

    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    ! Check for a KG_GPW calculation: KE for total rho and for molecular rho
    my_kg_gpw = .FALSE.
    IF(PRESENT(kg_gpw)) my_kg_gpw = kg_gpw

    CALL get_qs_env(qs_env=qs_env,&
         dft_control=dft_control,&
         matrix_h=matrix_h,&
         matrix_s=matrix_s,&
         pw_env=pw_env,&
         cell=cell,&
         para_env=para_env,&
         input=input,&
         error=error)
    
    scf_section => section_vals_get_subs_vals(input,"DFT%SCF",error=error)
    nspins=dft_control%nspins
    CPPrecondition(ASSOCIATED(matrix_h),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(matrix_s),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(rho),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(pw_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(SIZE(ks_matrix)>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(rho%rho_r_valid,cp_failure_level,routineP,error,failure)
    CPPrecondition(rho%rho_g_valid,cp_failure_level,routineP,error,failure)

    ! Check for GAPW method : additional terms for local densities
    gapw = dft_control%qs_control%gapw
    gapw_xc = dft_control%qs_control%gapw_xc
    IF(gapw_xc .AND. gapw) CALL cp_assert( .FALSE. , cp_failure_level,cp_assertion_failed,routineP,&
                           " GAPW and GAPW_XC are not compatible",error,failure)

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
           CALL pw_pool_init_coeff(auxbas_pw_pool, v_hartree_gspace, &
                                   use_data = COMPLEXDATA1D,&
                                   in_space = RECIPROCALSPACE, error=error)
           CALL pw_pool_init_coeff(auxbas_pw_pool, rho_tot_gspace,&
                                   use_data = COMPLEXDATA1D,&
                                   in_space = RECIPROCALSPACE, error=error)
       END IF

       CALL calc_rho_tot_gspace(rho_tot_gspace, qs_env, rho)

       output_unit = cp_print_key_unit_nr(logger,scf_section,"PRINT%TOTAL_DENSITIES",&
            extension=".scfLog",error=error)       
       IF (output_unit>0) THEN
          CALL print_densities(qs_env, rho_tot_gspace, rho, my_rho_xc, output_unit, error)
       END IF
       CALL cp_print_key_finished_output(output_unit,logger,scf_section,&
            "PRINT%TOTAL_DENSITIES", error=error)
       CALL pw_poisson_solve(poisson_env,rho_tot_gspace%pw,energy%hartree,&
            v_hartree_gspace%pw)
      
       ! In case decouple periodic images and/or apply restraints to charges
       CALL cp_ddapc_apply(qs_env,&
                           rho_tot_gspace,&
                           energy%hartree,&
                           energy%ddapc_restraint,&
                           v_hartree_gspace,&
                           calculate_forces=calculate_forces,&
                           error=error)

       IF (.NOT. just_energy) THEN
          CALL coeff_transform_space(v_hartree_gspace,&
               ks_env%v_hartree_rspace)
          CALL coeff_scale(ks_env%v_hartree_rspace, &
               ks_env%v_hartree_rspace%pw%pw_grid%dvol)
       ENDIF

       CALL pw_pool_give_back_coeff(auxbas_pw_pool,rho_tot_gspace,&
            error=error)
       CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_hartree_gspace,&
            error=error)

       CALL calc_v_sic_rspace(v_sic_rspace,energy,&
                              qs_env,dft_control,rho,poisson_env,just_energy,&
                              auxbas_pw_pool)

       IF(gapw) CALL Vh_1c_gg_integrals(qs_env,energy%hartree_1c)

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

       ! *** calculate the xc potential on the pw density ***       
       ! *** associates v_rspace_new if the xc potential needs to be computed.

       CALL qs_vxc_create(qs_env=qs_env,&
            vxc_rho=v_rspace_new, vxc_tau=v_tau_rspace, exc=energy%exc, &
            just_energy=just_energy, gapw_xc=gapw_xc, error=error)

       IF (gapw .OR. gapw_xc) THEN
          CALL calculate_vxc_atom(qs_env,just_energy,error)
       END IF

       ! ***  Single atom contributions ***
       IF (.NOT. just_energy) THEN
          IF (calculate_forces .AND. .NOT. gapw) THEN
             CALL integrate_v_core_rspace(ks_env%v_hartree_rspace, qs_env)
          END IF

          DO ispin=1,nspins
             ! copy the core ham into the ks matrix
             CALL transfer_matrix(matrix_h(1)%matrix,ks_matrix(ispin)%matrix)
          END DO

          IF (ASSOCIATED(v_rspace_new)) THEN
             IF(gapw_xc) THEN
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

                 !Now the Hartree potential to be integrated with the full basis
                 v_rspace_new(ispin)%pw%cr3d  =&
                      ks_env%v_hartree_rspace%pw%cr3d

                  ! add the hartree part due to v_rspace
                 CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin),&
                      p=rho%rho_ao(ispin),h=ks_matrix(ispin),&
                      qs_env=qs_env, &
                      calculate_forces=calculate_forces,&
                      gapw=gapw,error=error)

                 CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_rspace_new(ispin),&
                     error=error)

               END DO  ! ispin
             ELSE
               DO ispin=1,nspins
                  v_rspace_new(ispin)%pw%cr3d  =&
                       v_rspace_new(ispin)%pw%pw_grid%dvol * &
                       v_rspace_new(ispin)%pw%cr3d + &
                       ks_env%v_hartree_rspace%pw%cr3d
 
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
                 END SELECT 

                  ! add the part due to v_rspace
                  CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin),&
                       p=rho%rho_ao(ispin),h=ks_matrix(ispin),&
                       qs_env=qs_env, &
                       calculate_forces=calculate_forces,&
                       gapw=gapw,error=error)

                  CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_rspace_new(ispin),&
                       error=error)
               END DO
             END IF
             SELECT CASE (dft_control%sic_method_id) 
             CASE (sic_none)
             CASE (sic_mauri_us,sic_mauri_spz, sic_ad )
                CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_sic_rspace,error=error)
             END SELECT 
             DEALLOCATE(v_rspace_new,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          ELSE
             IF (dft_control%sic_method_id .NE. sic_none) THEN
                STOP "Nope"
             ENDIF
             DO ispin=1,nspins
                ! add only v_hartree
                CALL integrate_v_rspace(v_rspace=ks_env%v_hartree_rspace,&
                     p=rho%rho_ao(ispin),h=ks_matrix(ispin), qs_env=qs_env,&
                     calculate_forces=calculate_forces,&
                     gapw=gapw,error=error)
             END DO
          END IF

          IF (ASSOCIATED(v_tau_rspace)) THEN
             DO ispin=1,nspins
                 v_tau_rspace(ispin)%pw%cr3d =&
                      v_tau_rspace(ispin)%pw%pw_grid%dvol*&
                      v_tau_rspace(ispin)%pw%cr3d 

                 CALL integrate_v_rspace(v_rspace=v_tau_rspace(ispin),&
                      p=rho%rho_ao(ispin),h=ks_matrix(ispin),&
                      qs_env=qs_env,&
                      calculate_forces=calculate_forces,compute_tau=.TRUE., &
                      gapw=gapw,error=error)

              CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_tau_rspace(ispin),&
                   error=error)                 

              END DO
              DEALLOCATE(v_tau_rspace, stat=stat)
              CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ENDIF
       END IF  ! .NOT. just energy
  
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
          CALL integrate_vhg0_g0(qs_env,energy%hartree_2c,.FALSE.,calculate_forces)

!      ***  Energy, Integrals and Forces coming from the 3-centers terms
!      *** Initialize energy
          energy%hartree_3c = 0.0_dp
          CALL integrate_a_vhg0_b(qs_env,ks_matrix,rho%rho_ao,&
                       energy%hartree_3c,.FALSE.,calculate_forces,error=error)

       END IF

       IF(gapw .OR. gapw_xc) THEN
!      ***  Single atom contributions in the KS matrix ***
          CALL update_ks_atom(qs_env,ks_matrix,rho%rho_ao,&
               calculate_forces,error=error)
       ENDIF
!      
  
       energy%mulliken = 0.0_dp
       IF (dft_control%qs_control%mulliken_restraint) THEN
           IF (just_energy) THEN
               CALL mulliken_restraint(dft_control%qs_control%mulliken_restraint_control, &
                             para_env,matrix_s(1)%matrix, rho%rho_ao,energy=energy%mulliken, &
                             order_p=mulliken_order_p)
           ELSE
               CALL mulliken_restraint(dft_control%qs_control%mulliken_restraint_control,&
                             para_env,matrix_s(1)%matrix,rho%rho_ao,energy=energy%mulliken, &
                             ks_matrix=ks_matrix, order_p=mulliken_order_p)
           ENDIF
       ENDIF

       ! **** right now we can compute the orbital derivative here, as it depends currently only on the available
       ! **** Kohn-Sham matrix. This might change in the future, in which case more pieces might need to be assembled 
       ! **** from this routine, notice that this part of the calculation in not linear scaling
       ! **** right now this operation is only non-trivial because of occupation numbers and the restricted keyword
       IF (qs_env%requires_mo_derivs .AND. .NOT. just_energy) THEN
          CALL get_qs_env(qs_env,mo_derivs=mo_derivs,mos=mo_array)
          DO ispin=1,SIZE(mo_derivs)          
             CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,mo_coeff=mo_coeff, &
                             occupation_numbers=occupation_numbers )
             CALL cp_fm_get_info(mo_coeff,nrow_global=n,ncol_global=k)
             CALL cp_sm_fm_multiply(ks_matrix(ispin)%matrix,mo_coeff,&
                                          mo_derivs(ispin)%matrix, k)
             ALLOCATE(scaling_factor(SIZE(occupation_numbers)))
             scaling_factor=2.0_dp*occupation_numbers
             CALL cp_fm_column_scale(mo_derivs(ispin)%matrix,scaling_factor)
             DEALLOCATE(scaling_factor)

             IF (dft_control%restricted) THEN 
                 ! only the first mo_set are actual variables, but we still need both
                 CPPrecondition(ispin==1, cp_failure_level, routineP, error, failure)  
                 CPPrecondition(SIZE(mo_array)==2, cp_failure_level, routineP, error, failure)  
                 ! use a temporary array with the same size as the first spin for the second spin
                 CALL cp_fm_create(hc_tmp,mo_derivs(1)%matrix%matrix_struct)
                 CALL cp_fm_set_all(hc_tmp,0.0_dp)
                 CALL cp_fm_get_info(hc_tmp,ncol_global=k2)

                 CALL get_mo_set(mo_set=mo_array(2)%mo_set,mo_coeff=mo_coeff, &
                                 occupation_numbers=occupation_numbers )
                 CALL cp_fm_get_info(mo_coeff,ncol_global=k)
                 CALL cp_sm_fm_multiply(ks_matrix(2)%matrix,mo_coeff,&
                                        hc_tmp, k)
                 ALLOCATE(scaling_factor(k2))
                 scaling_factor=0.0_dp
                 scaling_factor(1:k)=2.0_dp*occupation_numbers(1:k)
                 CALL cp_fm_column_scale(hc_tmp,scaling_factor)
                 DEALLOCATE(scaling_factor)

                 ! sum the corresponding columns (relies on the final columns of hc_tmp being zero)
                 CALL cp_fm_scale_and_add(1.0_dp, mo_derivs(1)%matrix, 1.0_dp, hc_tmp)
                 CALL cp_fm_release(hc_tmp)
             ENDIF

          ENDDO
       ENDIF

       ! adds s2_restraint energy and orbital derivatives
       IF (dft_control%qs_control%s2_restraint) THEN
           CPPrecondition(dft_control%nspins == 2,cp_failure_level,routineP,error,failure)
           CPPrecondition(qs_env%requires_mo_derivs,cp_failure_level,routineP,error,failure)
           ! forces are not implemented (not difficult, but ... )
           CPPrecondition(.NOT. calculate_forces, cp_failure_level,routineP,error,failure) 
           CALL get_qs_env(qs_env,mo_derivs=mo_derivs,mos=mo_array)
           CALL s2_restraint(mo_array, matrix_s, mo_derivs, energy%s2_restraint, &
                             dft_control%qs_control%s2_restraint_control, just_energy, error) 
       ELSE
          energy%s2_restraint = 0.0_dp
       ENDIF


       ! sum all energy terms to obtain the total energy

       energy%total = energy%core_overlap + energy%core_self + &
            energy%core + energy%hartree + &
            energy%hartree_1c + energy%hartree_2c + energy%hartree_3c +&
            energy%exc + energy%exc1 + &
            energy%qmmm_el + energy%mulliken + &
            energy%kg_gpw_ekin_mol + energy%ddapc_restraint + energy%s2_restraint

       output_unit=cp_print_key_unit_nr(logger,input,"DFT%SCF%PRINT%DETAILED_ENERGY",&
            extension=".scfLog",error=error)
       IF (output_unit>0) THEN
          IF(my_kg_gpw) THEN
            WRITE (UNIT=output_unit,FMT="(/,(T3,A,T60,F20.10))")&
               "Core Hamiltonian energy of separated molecules:",energy%core,&
               "Hartree energy of the total density:           ",energy%hartree,&
               "XC and Kin. energy of the total density:       ",energy%exc,&
               "Kin. energy of the molecular densities:        ",energy%kg_gpw_ekin_mol
          ELSE
            WRITE (UNIT=output_unit,FMT="(/,(T3,A,T60,F20.10))")&
               "Core Hamiltonian energy:                       ",energy%core,&
               "Hartree energy:                                ",energy%hartree,&
               "Exchange-correlation energy:                   ",energy%exc
          END IF
          IF(gapw) THEN
             WRITE (UNIT=output_unit,FMT="(/,(T3,A,T60,F20.10))")&
                "GAPW| Exc from hard and soft atomic rho1:      ",energy%exc1,&
                "GAPW| local Eh = 1 center integrals:           ",energy%hartree_1c,&
                "GAPW| local Eh = 2 center integrals:           ",energy%hartree_2c,&
                "GAPW| local Eh = 3 centers integrals:          ",energy%hartree_3c
          END IF
          IF(gapw_xc) THEN
             WRITE (UNIT=output_unit,FMT="(/,(T3,A,T60,F20.10))")&
                "GAPW| Exc from hard and soft atomic rho1:      ",energy%exc1
          END IF
          IF  (qs_env%qmmm) THEN
             WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
               "QM/MM Electrostatic energy:                    ",energy%qmmm_el
          END IF
          IF (dft_control%qs_control%mulliken_restraint) THEN
              WRITE (UNIT=output_unit,FMT="(T3,A,T40,2F20.10)")&
                "Mulliken restraint (order_p,energy) : ",mulliken_order_p,energy%mulliken
          ENDIF
          IF (dft_control%qs_control%ddapc_restraint) THEN
             ddapc_order_p = dft_control%qs_control%ddapc_restraint_control%ddapc_order_p
              WRITE (UNIT=output_unit,FMT="(T3,A,T40,2F20.10)")&
                "DDAPC restraint (order_p,energy) : ",ddapc_order_p,energy%ddapc_restraint
          ENDIF
          IF (dft_control%qs_control%s2_restraint) THEN
             s2_order_p = dft_control%qs_control%s2_restraint_control%s2_order_p
              WRITE (UNIT=output_unit,FMT="(T3,A,T40,2F20.10)")&
                "S2 restraint (order_p,energy) : ",s2_order_p,energy%s2_restraint
          ENDIF

       END IF
       CALL cp_print_key_finished_output(output_unit,logger,input,&
            "DFT%SCF%PRINT%DETAILED_ENERGY", error=error)

    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE qs_ks_build_kohn_sham_matrix

! *****************************************************************************

  SUBROUTINE calc_v_sic_rspace(v_sic_rspace,energy,&
                               qs_env,dft_control,rho,poisson_env,just_energy,&
                               auxbas_pw_pool,error)

    TYPE(coeff_type), INTENT(INOUT)          :: v_sic_rspace
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    LOGICAL, INTENT(IN)                      :: just_energy
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(*), PARAMETER :: routineN = 'calc_v_sic_rspace', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nelec, nelec_a, nelec_b
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: ener, scaling
    TYPE(coeff_type)                         :: work_rho, work_v
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array

    failure = .FALSE.
    NULLIFY(mo_array)

    IF (dft_control%sic_method_id == sic_none) RETURN
    CALL cp_assert(.NOT. dft_control%qs_control%gapw, cp_failure_level,cp_assertion_failed,routineP,&
                   "sic and GAPW not yet compatible",error,failure)

    ! OK, right now we like two spins to do sic, could be relaxed for AD 
    CPPrecondition(dft_control%nspins == 2,cp_failure_level,routineP,error,failure)

    CALL pw_pool_init_coeff(auxbas_pw_pool, work_rho, &
                            use_data = COMPLEXDATA1D,&
                            in_space = RECIPROCALSPACE, error=error)
    CALL pw_pool_init_coeff(auxbas_pw_pool, work_v,&
                            use_data = COMPLEXDATA1D,&
                            in_space = RECIPROCALSPACE, error=error)
    
    ! Hartree sic corrections
    SELECT CASE ( dft_control%sic_method_id )
    CASE ( sic_mauri_us, sic_mauri_spz )
       CALL coeff_subtract(rho%rho_g(1),rho%rho_g(2),work_rho)
       CALL pw_poisson_solve(poisson_env,work_rho%pw, ener, work_v%pw)
       energy%hartree=energy%hartree - dft_control % sic_scaling_a * ener
    CASE ( sic_ad )
       ! find out how many elecs we have
       CALL get_qs_env(qs_env,mos=mo_array)
       CALL get_mo_set(mo_set=mo_array(1)%mo_set,nelectron=nelec_a)
       CALL get_mo_set(mo_set=mo_array(2)%mo_set,nelectron=nelec_b)
       nelec = nelec_a + nelec_b
       CALL coeff_add(rho%rho_g(1),rho%rho_g(2),work_rho)
       scaling = 1.0_dp / REAL(nelec,KIND=dp)
       CALL coeff_scale(work_rho,scaling)
       CALL pw_poisson_solve(poisson_env,work_rho%pw, ener, work_v%pw)
       energy%hartree=energy%hartree - dft_control % sic_scaling_a * nelec * ener
    CASE DEFAULT
       CALL cp_assert( .FALSE., cp_failure_level,cp_assertion_failed,routineP,&
            "Unknown sic method id",error,failure)
    END SELECT

    IF (.NOT. just_energy) THEN
       CALL pw_pool_init_coeff(auxbas_pw_pool,v_sic_rspace,&
                               use_data=REALDATA3D, in_space=REALSPACE)
       CALL coeff_transform_space(work_v, v_sic_rspace)
       ! also take into account the scaling (in addition to the volume element)
       CALL coeff_scale(v_sic_rspace, &
            dft_control % sic_scaling_a * v_sic_rspace%pw%pw_grid%dvol )
    ENDIF

    CALL pw_pool_give_back_coeff(auxbas_pw_pool,work_rho,error=error)
    CALL pw_pool_give_back_coeff(auxbas_pw_pool,work_v,error=error)

  END SUBROUTINE calc_v_sic_rspace

! *****************************************************************************

  SUBROUTINE calc_rho_tot_gspace(rho_tot_gspace, qs_env, rho, error)
    TYPE(coeff_type), INTENT(INOUT)          :: rho_tot_gspace
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(*), PARAMETER :: routineN = 'calc_rho_tot_gspace', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin
    LOGICAL                                  :: failure
    TYPE(coeff_type), POINTER                :: rho0_s_gs, rho_core
    TYPE(qs_charges_type), POINTER           :: qs_charges

    NULLIFY (rho_core, rho0_s_gs)
    CALL get_qs_env(qs_env=qs_env,&
                    rho_core=rho_core,&
                    rho0_s_gs=rho0_s_gs,&
                    qs_charges=qs_charges)
    
    IF(qs_env%dft_control%qs_control%gapw) THEN
       CPPrecondition(ASSOCIATED(rho0_s_gs), cp_failure_level, routineP, error, failure)
       CALL coeff_add(rho0_s_gs, rho%rho_g(1), rho_tot_gspace)
    ELSE
       CALL coeff_add(rho_core, rho%rho_g(1), rho_tot_gspace)
    END IF
    DO ispin=2, qs_env%dft_control%nspins
       CALL coeff_sumup(rho%rho_g(ispin), rho_tot_gspace)
    END DO
    qs_charges%total_rho_gspace = calculate_total_rho(rho_tot_gspace)

  END SUBROUTINE calc_rho_tot_gspace

! *****************************************************************************

  SUBROUTINE print_densities(qs_env, rho_tot_gspace, rho, rho_xc, output_unit, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(coeff_type), INTENT(IN)             :: rho_tot_gspace
    TYPE(qs_rho_type), POINTER               :: rho, rho_xc
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    INTEGER                                  :: ispin, n_electrons
    REAL(dp)                                 :: tot1_h, tot1_s
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(qs_charges_type), POINTER           :: qs_charges

    NULLIFY(qs_charges, atomic_kind_set)
    CALL get_qs_env(qs_env=qs_env, &
         atomic_kind_set=atomic_kind_set, &
         qs_charges=qs_charges)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
         nelectron=n_electrons)
    n_electrons = n_electrons - qs_env%dft_control%charge
    
    WRITE (UNIT=output_unit,FMT="(/,(T3,A,T40,2F20.10))")&
         "Total electronic density (r-space): ",&
         accurate_sum(rho%tot_rho_r),&
         accurate_sum(rho%tot_rho_r) + &
         REAL(n_electrons,dp),&
         "Total core charge density (r-space):",&
         qs_charges%total_rho_core_rspace,&
         qs_charges%total_rho_core_rspace - REAL(n_electrons,dp)
    IF(qs_env%dft_control%qs_control%gapw ) THEN
       tot1_h =  qs_charges%total_rho1_hard(1)
       tot1_s =  qs_charges%total_rho1_soft(1)
       DO ispin=2,qs_env%dft_control%nspins
          tot1_h = tot1_h + qs_charges%total_rho1_hard(ispin)
          tot1_s = tot1_s + qs_charges%total_rho1_soft(ispin)
       END DO
       WRITE (UNIT=output_unit,FMT="(/,(T3,A,T40,2F20.10))")&
            "Hard and soft densities (Lebedev):",&
            tot1_h, tot1_s
       WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
            "Total Rho_soft + Rho1_hard - Rho1_soft (r-space): ",&
            accurate_sum(rho%tot_rho_r)+ tot1_h - tot1_s ,&
            "Total charge density (r-space):      ",&
            accurate_sum(rho%tot_rho_r)+ tot1_h - tot1_s &
            + qs_charges%total_rho_core_rspace,&
            "Total Rho_soft + Rho0_soft (g-space):",&
            qs_charges%total_rho_gspace
    ELSE IF( qs_env%dft_control%qs_control%gapw_xc) THEN
       tot1_h =  qs_charges%total_rho1_hard(1)
       tot1_s =  qs_charges%total_rho1_soft(1)
       DO ispin=2,qs_env%dft_control%nspins
          tot1_h = tot1_h + qs_charges%total_rho1_hard(ispin)
          tot1_s = tot1_s + qs_charges%total_rho1_soft(ispin)
       END DO
       WRITE (UNIT=output_unit,FMT="(/,(T3,A,T40,2F20.10))")&
            "Hard and soft densities (Lebedev):",&
            tot1_h, tot1_s
       WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
            "Total Rho_soft + Rho1_hard - Rho1_soft (r-space): ",&
            accurate_sum(rho_xc%tot_rho_r)+ tot1_h - tot1_s 
    ELSE  
       WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
            "Total charge density (r-space):     ",&
            accurate_sum(rho%tot_rho_r)+&
            qs_charges%total_rho_core_rspace,&
            "Total charge density (g-space):     ",&
            qs_charges%total_rho_gspace
    END IF
           
  END SUBROUTINE print_densities

! *****************************************************************************

!!****f* qs_ks_methods/qs_ks_update_qs_env [1.0] *
!!
!!   NAME
!!     qs_ks_update_qs_env
!!
!!   FUNCTION
!!     updates the Kohn Sham matrix of the given qs_env (facility method)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - ks_env: the ks_env that old all the temporary objects that
!!       the calculation of the KS matrix needs
!!     - qs_env: the qs_env to update
!!     - calculate_forces: if true calculate the quantities needed
!!       to calculate the forces. Defaults to false.
!!     - just_energy: if true updates the energies but not the
!!       ks matrix. Defaults to false
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE qs_ks_update_qs_env(ks_env,qs_env, kg_env, &
     calculate_forces, just_energy, error)
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(kg_environment_type), OPTIONAL, &
      POINTER                                :: kg_env
    LOGICAL, INTENT(in), OPTIONAL            :: calculate_forces, just_energy
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ks_update_qs_env', &
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
     CALL timeset(routineN//" (forces)","I","",handle)
  ELSE 
     CALL timeset(routineN,"I","",handle)
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
        ! the ks matrix we be rebuild so this is fine now
        ks_env%potential_changed=.FALSE.
        IF ( dft_control%qs_control%semi_empirical ) THEN
          CALL build_se_fock_matrix(ks_env,qs_env=qs_env,&
               ks_matrix=matrix_ks,rho=rho,energy=energy,&
               calculate_forces=c_forces,just_energy=energy_only,&
               error=error)
        ELSEIF (PRESENT(kg_env)) THEN
          CALL get_qs_env(qs_env=qs_env, kg_sub_pw_env=kg_sub_pw_env)
          CALL get_kg_env(kg_env=kg_env, kg_fm_set=kg_fm_set)
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
           ks_env%forces_up_to_date=.NOT.c_forces
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

  CALL timestop(0.0_dp,handle)

END SUBROUTINE qs_ks_update_qs_env
!***************************************************************************

!!****f* qs_ks_methods/qs_ks_did_change [1.0] *
!!
!!   NAME
!!     qs_ks_did_change
!!
!!   FUNCTION
!!     tells that some of the things relevant to the ks calculation
!!     did change. has to be called when changes happen otherwise
!!     the calculation will give wrong results.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - ks_env: the environement that is informed about the changes
!!     - s_mstruct_changed: if true it means that the structure of the
!!       overlap matrix has changed
!!       (atoms have moved)
!!     - rho_changed: if true it means that the density has changed
!!     - full_reset: if true everything has changed
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE qs_ks_did_change(ks_env,s_mstruct_changed,rho_changed,potential_changed,full_reset,error)
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    LOGICAL, INTENT(in), OPTIONAL            :: s_mstruct_changed, &
                                                rho_changed, &
                                                potential_changed, full_reset
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ks_did_change', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

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
  END IF
END SUBROUTINE qs_ks_did_change
!***************************************************************************

!!****f* qs_ks_methods/qs_vxc_create [1.0] *
!!
!!   NAME
!!     qs_vxc_create
!!
!!   FUNCTION
!!     calculates and allocates the xc potential, already reducing it to 
!!     the dependence on rho and the one on tau
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - vxc_rho: will contain the v_xc part that depend on rho
!!       (if one of the choosen xc functionals has it it is allocated and you
!!       are responsible for it)
!!     - vxc_tau: will contain the kinetic (tau) part of v_xc
!!       (if one of the choosen xc functionals has it it is allocated and you
!!       are responsible for it)
!!     - qs_env: qs_environment to get all the needed things (rho,...)
!!     - calculate_forces: if the forces should be calculated (so that you
!!       might be able to add forces to the atoms is you do strange stuff in
!!       the xc
!!     - just_energy: if true calculates just the energy, and does not 
!!       allocate v_*_rspace
!!     - harris: If true uses the rho structure of the Harris environment
!!     - nsc_force: Uses rho_diff in case the non-self-consistent force should be
!!                  calculated. Works only in combination with the Harris functional!
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     - 05.2002 modified to use the mp_allgather function each pe
!!       computes only part of the grid and this is broadcasted to all 
!!       instead of summed. 
!!       This scales significantly better (e.g. factor 3 on 12 cpus
!!       32 H2O) [Joost VdV]
!!     - moved to qs_ks_methods [fawzi]
!!     - sic alterations [Joost VandeVondele]
!!
!!*** *********************************************************************
SUBROUTINE qs_vxc_create(vxc_rho,vxc_tau, exc, qs_env,&
      just_energy, harris, nsc_force, gapw_xc, error)
    TYPE(coeff_type), DIMENSION(:), POINTER  :: vxc_rho, vxc_tau
    REAL(KIND=dp)                            :: exc
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(in), OPTIONAL            :: just_energy, harris, &
                                                nsc_force, gapw_xc
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_vxc_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, nelec_spin(2), &
                                                stat
    LOGICAL                                  :: failure, harris_flag, &
                                                my_gapw_xc, my_just_energy, &
                                                nsc_flag, sic_scaling_b_zero, &
                                                uf_grid
    REAL(KIND=dp)                            :: exc_m, my_scaling, nelec_s_inv
    TYPE(cell_type), POINTER                 :: cell
    TYPE(coeff_type), DIMENSION(:), POINTER  :: rho_m_gspace, rho_m_rspace
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(harris_env_type), POINTER           :: harris_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: my_vxc_rho, my_vxc_tau, &
                                                rho_g, rho_r, tau
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool, xc_pw_pool
    TYPE(pw_type), POINTER                   :: tmp_g, tmp_g2, tmp_pw
    TYPE(qs_rho_type), POINTER               :: rho_struct
    TYPE(section_vals_type), POINTER         :: input, xc_section

!   ---------------------------------------------------------------------------

  CALL timeset( routineN ,"I","",handle)

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

  CALL get_qs_env(qs_env=qs_env, dft_control=dft_control,input=input, &
       pw_env=pw_env, cell=cell, error=error)
  xc_section => section_vals_get_subs_vals(input,"DFT%XC",error=error)
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
    CALL get_qs_env(qs_env=qs_env, rho=rho_struct, error=error)
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
     uf_grid=(auxbas_pw_pool%pw_grid%id_nr/=xc_pw_pool%pw_grid%id_nr)
     
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
!FM        DO ispin=1,dft_control%nspins
!FM           CALL pw_pool_create_pw(xc_pw_pool,rho_r(ispin)%pw,&
!FM                in_space=REALSPACE, use_data=REALDATA3D,error=error)
!FM           call pw_zero(rho_r(ispin)%pw)
!FM           CALL pw_prolongate_s3(rho_struct%rho_r(ispin)%pw,rho_r(ispin)%pw,&
!FM                auxbas_pw_pool,param_section=interp_section,error=error)
!FM        END DO
!FM        
!FM        IF (rho_struct%tau_r_valid) THEN
!FM           ALLOCATE(tau(dft_control%nspins),stat=stat)
!FM           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!FM           DO ispin=1,dft_control%nspins
!FM              CALL pw_pool_create_pw(xc_pw_pool,tau(ispin)%pw,&
!FM                   in_space=REALSPACE, use_data=REALDATA3D,error=error)
!FM              CALL pw_zero(rho_r(ispin)%pw)
!FM              CALL pw_prolongate_s3(rho_struct%tau_r(ispin)%pw,tau(ispin)%pw,&
!FM                   auxbas_pw_pool,param_section=interp_section,error=error)
!FM           END DO
!FM        END IF
!FM
!FM        ! for gradient corrected functional the density in g space might
!FM        ! be useful so if we have it, we pass it in
!FM        IF ( rho_struct%rho_g_valid ) THEN
!FM           ALLOCATE(rho_g(dft_control%nspins),stat=stat)
!FM           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!FM           DO ispin=1,dft_control%nspins
!FM              CALL pw_pool_create_pw(xc_pw_pool,rho_g(ispin)%pw,&
!FM                   in_space=REALSPACE, use_data=REALDATA3D,error=error)
!FM              CALL pw_transfer(rho_struct%rho_g(ispin)%pw,rho_g(ispin)%pw)
!FM           END DO
!FM        END IF
        CPPrecondition(rho_struct%rho_g_valid,cp_failure_level,routineP,error,failure)
        ALLOCATE(rho_g(dft_control%nspins),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        DO ispin=1,dft_control%nspins
           CALL pw_pool_create_pw(xc_pw_pool,rho_g(ispin)%pw,&
                in_space=RECIPROCALSPACE, use_data=COMPLEXDATA1D,error=error)
           CALL pw_transfer(rho_struct%rho_g(ispin)%pw,rho_g(ispin)%pw)
        END DO
        DO ispin=1,dft_control%nspins
           CALL pw_pool_create_pw(xc_pw_pool,rho_r(ispin)%pw,&
                in_space=REALSPACE, use_data=REALDATA3D,error=error)
           CALL pw_transfer(rho_g(ispin)%pw,rho_r(ispin)%pw)
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
              CALL pw_transfer(tau(ispin)%pw,tmp_g)
              CALL pw_transfer(tmp_g,tmp_g2)
              CALL pw_transfer(tmp_g2,tmp_pw)
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
             error=error)
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
           CALL coeff_init(vxc_rho(ispin),pw=my_vxc_rho(ispin)%pw,&
                error=error)
        END DO
        DEALLOCATE(my_vxc_rho,stat=stat)
     END IF
     IF (ASSOCIATED(my_vxc_tau)) THEN
        ALLOCATE(vxc_tau(dft_control%nspins),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        DO ispin=1,dft_control%nspins
           CALL coeff_init(vxc_tau(ispin),pw=my_vxc_tau(ispin)%pw,&
                error=error)
        END DO
        DEALLOCATE(my_vxc_tau,stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     END IF

     ! compute again the xc but now for Exc(m,o) and the opposite sign
     IF (dft_control%sic_method_id .EQ. sic_mauri_spz .AND. .NOT. sic_scaling_b_zero) THEN
        ALLOCATE(rho_m_rspace(2),rho_m_gspace(2))
        CALL pw_pool_init_coeff(xc_pw_pool, rho_m_gspace(1),&
                                   use_data = COMPLEXDATA1D,&
                                   in_space = RECIPROCALSPACE, error=error)
        CALL pw_pool_init_coeff(xc_pw_pool, rho_m_rspace(1),&
                                   use_data = REALDATA3D,&
                                   in_space = REALSPACE, error=error)
        CALL coeff_subtract(rho_struct%rho_r(1),rho_struct%rho_r(2),rho_m_rspace(1))
        CALL coeff_subtract(rho_struct%rho_g(1),rho_struct%rho_g(2),rho_m_gspace(1))
        ! bit sad, these will be just zero...
        CALL pw_pool_init_coeff(xc_pw_pool, rho_m_gspace(2),&
                                   use_data = COMPLEXDATA1D,&
                                   in_space = RECIPROCALSPACE, error=error)
        CALL pw_pool_init_coeff(xc_pw_pool, rho_m_rspace(2),&
                                   use_data = REALDATA3D,&
                                   in_space = REALSPACE, error=error)
        CALL coeff_zero(rho_m_rspace(2))
        CALL coeff_zero(rho_m_gspace(2))

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
           vxc_rho(1)%pw%cr3d=vxc_rho(1)%pw%cr3d-dft_control % sic_scaling_b * my_vxc_rho(1)%pw%cr3d
           vxc_rho(2)%pw%cr3d=vxc_rho(2)%pw%cr3d+dft_control % sic_scaling_b * my_vxc_rho(1)%pw%cr3d ! 1=m
           CALL pw_release(my_vxc_rho(1)%pw) 
           CALL pw_release(my_vxc_rho(2)%pw) 
           DEALLOCATE(my_vxc_rho,stat=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        ENDIF

        DO ispin=1,2
           CALL pw_pool_give_back_coeff(xc_pw_pool,rho_m_rspace(ispin),&
                                        error=error)                 
           CALL pw_pool_give_back_coeff(xc_pw_pool,rho_m_gspace(ispin),&
                                        error=error)                 
        ENDDO
        DEALLOCATE(rho_m_rspace)
        DEALLOCATE(rho_m_gspace)
        
     ENDIF

     ! now we have - sum_s N_s * Exc(rho_s/N_s,0)
     IF ( dft_control % sic_method_id .EQ. sic_ad  .AND. .NOT. sic_scaling_b_zero ) THEN

        ! find out how many elecs we have
        CALL get_qs_env(qs_env,mos=mo_array)
        CALL get_mo_set(mo_set=mo_array(1)%mo_set,nelectron=nelec_spin(1))
        CALL get_mo_set(mo_set=mo_array(2)%mo_set,nelectron=nelec_spin(2))

        ALLOCATE(rho_m_rspace(2),rho_m_gspace(2))
        DO ispin=1,2
           CALL pw_pool_init_coeff(xc_pw_pool, rho_m_gspace(ispin),&
                                   use_data = COMPLEXDATA1D,&
                                   in_space = RECIPROCALSPACE, error=error)
           CALL pw_pool_init_coeff(xc_pw_pool, rho_m_rspace(ispin),&
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
           CALL coeff_copy(rho_struct%rho_r(ispin),rho_m_rspace(1))
           CALL coeff_copy(rho_struct%rho_g(ispin),rho_m_gspace(1))
           CALL coeff_scale(rho_m_rspace(1),nelec_s_inv)
           CALL coeff_scale(rho_m_gspace(1),nelec_s_inv)
           CALL coeff_zero(rho_m_rspace(2))
           CALL coeff_zero(rho_m_gspace(2))

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
              vxc_rho(ispin)%pw%cr3d=vxc_rho(ispin)%pw%cr3d-dft_control % sic_scaling_b * my_vxc_rho(1)%pw%cr3d
              CALL pw_release(my_vxc_rho(1)%pw)
              CALL pw_release(my_vxc_rho(2)%pw)
              DEALLOCATE(my_vxc_rho,stat=stat)
              CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           ENDIF
        ENDDO

        DO ispin=1,2
           CALL pw_pool_give_back_coeff(xc_pw_pool,rho_m_rspace(ispin),&
                                        error=error)
           CALL pw_pool_give_back_coeff(xc_pw_pool,rho_m_gspace(ispin),&
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
           CALL pw_release(my_vxc_rho(1)%pw) 
           CALL pw_release(my_vxc_rho(2)%pw) 
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
              CALL pw_transfer(vxc_rho(ispin)%pw,tmp_g)
              CALL pw_transfer(tmp_g,tmp_g2)
              CALL pw_transfer(tmp_g2,tmp_pw)
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
              CALL pw_transfer(vxc_tau(ispin)%pw,tmp_g)
              CALL pw_transfer(tmp_g,tmp_g2)
              CALL pw_transfer(tmp_g2,tmp_pw)
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

  CALL timestop(0.0_dp,handle)

END SUBROUTINE qs_vxc_create

! *****************************************************************************
!!****f* qs_ks_methods/calculate_w_matrix_1 [1.0] *
!!
!!    NAME
!!      calculate_w_matrix_1
!!
!!    FUNCTION
!!      Calculate the W matrix from the MO eigenvectors, MO eigenvalues,
!!      and the MO occupation numbers. Only works if they are eigenstates
!!
!!    NOTE
!!      -
!!
!!    INPUTS
!!       mo_set : type containing the full matrix of the MO and the eigenvalues
!!       w_matrix : sparse matrix 
!!       error
!!
!!    AUTHOR
!!       MK
!!
!!    MODIFICATION HISTORY
!!        Creation (03.03.03,MK)
!!        Modification that computes it as a full block, several times (e.g. 20) 
!!              faster at the cost of some additional memory
!!
!!*** *************************************************************************
  SUBROUTINE calculate_w_matrix_1(mo_set,w_matrix,error)

    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(real_matrix_type), POINTER          :: w_matrix
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_w_matrix_1', &
      routine = "SUBROUTINE "//routineN//" (MODULE "//moduleN//")"

    INTEGER                                  :: handle, imo, istat
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: eigocc
    TYPE(cp_fm_type), POINTER                :: weighted_vectors

!   ---------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.
    NULLIFY(weighted_vectors)

    CALL set_matrix(w_matrix,0.0_dp)
    CALL cp_fm_create(weighted_vectors,mo_set%mo_coeff%matrix_struct,"weighted_vectors")
    CALL cp_fm_to_fm(mo_set%mo_coeff,weighted_vectors)

    ! scale every column with the occupation
    ALLOCATE(eigocc(mo_set%homo),stat=istat)
    CPPostcondition(istat==0,cp_failure_level,routineN,error,failure)

    DO imo=1,mo_set%homo
       eigocc(imo) = mo_set%eigenvalues(imo)*mo_set%occupation_numbers(imo)
    ENDDO
    CALL cp_fm_column_scale(weighted_vectors,eigocc)
    DEALLOCATE(eigocc)

    CALL cp_sm_plus_fm_fm_t(sparse_matrix=w_matrix,&
                            matrix_v=mo_set%mo_coeff,&
                            matrix_g=weighted_vectors,&
                            ncol=mo_set%homo)

    CALL cp_fm_release(weighted_vectors)
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE calculate_w_matrix_1

! *****************************************************************************
!!****f* qs_ks_methods/calculate_w_matrix_2 [1.0] *
!!
!!    NAME
!!      calculate_w_matrix_2
!!
!!    FUNCTION
!!      Calculate the W matrix from the MO coefs, MO derivs 
!!       could overwrite the mo_derivs for increased memory efficiency
!!
!!    NOTE
!!      -
!!
!!    INPUTS
!!       mo_set : type containing the full matrix of the MO coefs
!!       mo_deriv:
!!       w_matrix : sparse matrix 
!!       s_matrix: sparse matrix for the overlap
!!       error
!!
!!    AUTHOR
!!       MK
!!
!!    MODIFICATION HISTORY
!!        Creation (JV)
!!
!!*** *************************************************************************
  
  SUBROUTINE calculate_w_matrix_2(mo_set,mo_deriv,w_matrix,s_matrix,error)

    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_fm_type), POINTER                :: mo_deriv
    TYPE(real_matrix_type), POINTER          :: w_matrix, s_matrix
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    LOGICAL, PARAMETER                       :: check_gradient = .FALSE., &
                                                do_symm = .FALSE.

    INTEGER                                  :: handle, ncol_block, &
                                                ncol_global, nrow_block, &
                                                nrow_global
    LOGICAL                                  :: failure
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: gradient, h_block, h_block_t, &
                                                weighted_vectors

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_w_matrix_2","I","",handle)
    failure=.FALSE.
    NULLIFY(weighted_vectors,h_block,fm_struct_tmp)

    CALL cp_fm_get_info(matrix=mo_set%mo_coeff,&
                        ncol_global=ncol_global,&
                        nrow_global=nrow_global,&
                        nrow_block=nrow_block,&
                        ncol_block=ncol_block)

    CALL cp_fm_create(weighted_vectors,mo_set%mo_coeff%matrix_struct,"weighted_vectors")
    CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ncol_global,  ncol_global=ncol_global, &
                                            para_env=mo_set%mo_coeff%matrix_struct%para_env, &
                                            context=mo_set%mo_coeff%matrix_struct%context)
    CALL cp_fm_create(h_block,fm_struct_tmp, name="h block")
    IF (do_symm) CALL cp_fm_create(h_block_t,fm_struct_tmp, name="h block t")
    CALL cp_fm_struct_release(fm_struct_tmp)

    ! the convention seems to require the half here, the factor of two is presumably taken care of 
    ! internally in qs_core_hamiltonian
    CALL cp_fm_gemm('T','N',ncol_global,ncol_global,nrow_global,0.5_dp, &
                    mo_set%mo_coeff,mo_deriv,0.0_dp,h_block)

    IF (do_symm) THEN
       ! at the minimum things are anyway symmetric, but numerically it might not be the case
       ! needs some investigation to find out if using this is better
       CALL cp_fm_transpose(h_block,h_block_t)
       CALL cp_fm_scale_and_add(0.5_dp,h_block,0.5_dp,h_block_t)
    ENDIF

    ! this could overwrite the mo_derivs to save the weighted_vectors
    CALL cp_fm_gemm('N','N',nrow_global,ncol_global,ncol_global,1.0_dp, &
                    mo_set%mo_coeff,h_block,0.0_dp,weighted_vectors)

    CALL set_matrix(w_matrix,0.0_dp)
    CALL cp_sm_plus_fm_fm_t(sparse_matrix=w_matrix,&
                            matrix_v=mo_set%mo_coeff,&
                            matrix_g=weighted_vectors,&
                            ncol=mo_set%homo)

    IF (check_gradient) THEN
       CALL cp_fm_create(gradient,mo_set%mo_coeff%matrix_struct,"gradient")
       CALL cp_sm_fm_multiply(s_matrix,weighted_vectors,&
                                          gradient, ncol_global)
       WRITE(6,*) " maxabs difference ", MAXVAL(ABS(mo_deriv%local_data-2.0_dp*gradient%local_data))
       CALL cp_fm_release(gradient)
    ENDIF

    IF (do_symm) CALL cp_fm_release(h_block_t)
    CALL cp_fm_release(weighted_vectors)
    CALL cp_fm_release(h_block)
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE calculate_w_matrix_2

END MODULE qs_ks_methods

! *****************************************************************************
