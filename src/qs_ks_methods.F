!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

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
!!
!!   SOURCE
!****************************************************************************
MODULE qs_ks_methods
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE coefficient_types,               ONLY: coeff_add,&
                                             coeff_allocate,&
                                             coeff_copy,&
                                             coeff_deallocate,&
                                             coeff_reduce,&
                                             coeff_sumup,&
                                             coeff_transform_space,&
                                             coeff_type,&
                                             coeff_zero
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error,&
                                             cp_unimplemented_error
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cube_utils,                      ONLY: cube_info_type
  USE dft_types,                       ONLY: dft_control_type
  USE gaussian_gridlevels,             ONLY: gridlevel_info_type
  USE greens_fn,                       ONLY: pw_green_fn_init
  USE hartree,                         ONLY: calculate_hartree
  USE kinds,                           ONLY: wp=>dp
  USE l_utils,                         ONLY: l_info_type
  USE message_passing,                 ONLY: mp_allgather,&
                                             mp_max,&
                                             mp_range,&
                                             mp_sum,&
                                             mp_sync
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_give_back_coeff,&
                                             pw_pool_init_coeff,&
                                             pw_pool_p_type,&
                                             pw_pool_type,&
                                             pw_pools_give_back_coeffs,&
                                             pw_pools_init_coeffs
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             COMPLEXDATA3D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_type
  USE pws,                             ONLY: cleanup_pw_poisson_solver,&
                                             init_pw_poisson_solver
  USE qs_charges_types,                ONLY: qs_charges_type
  USE qs_collocate_density,            ONLY: calculate_total_rho
  USE qs_core_energies,                ONLY: calculate_ecore,&
                                             calculate_ecore_overlap,&
                                             calculate_ecore_self
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_integrate_potential,          ONLY: integrate_v_core_rspace,&
                                             integrate_v_rspace
  USE qs_ks_types,                     ONLY: qs_ks_env_p_type,&
                                             qs_ks_env_type,&
                                             qs_ks_release,&
                                             qs_ks_retain
  USE qs_rho_methods,                  ONLY: qs_rho_update_rho
  USE qs_rho_types,                    ONLY: qs_rho_did_change,&
                                             qs_rho_type
  USE simulation_cell,                 ONLY: cell_type,&
                                             get_cell
  USE sparse_matrix_types,             ONLY: add_matrices,&
                                             deallocate_matrix,&
                                             first_block_node,&
                                             get_block_node,&
                                             get_matrix_info,&
                                             next_block_node,&
                                             real_block_node_type,&
                                             real_matrix_p_type,&
                                             replicate_matrix_structure,&
                                             set_matrix,&
                                             symmetrise_diagonal_blocks,&
                                             transfer_matrix
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc,                              ONLY: xc_calculate_pw2
  IMPLICIT NONE

  PRIVATE

  LOGICAL, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER :: moduleN='qs_ks_methods'
  INTEGER, SAVE, PRIVATE :: last_ks_id_nr=0

  PUBLIC :: qs_ks_create, qs_ks_did_change,&
       qs_ks_update_qs_env, qs_ks_build_kohn_sham_matrix
!***
!****************************************************************************

CONTAINS
!!****f* qs_ks_methods/qs_ks_create [1.0] *
!!
!!   NAME
!!     qs_ks_create
!!
!!   SYNOPSIS
!!     Subroutine qs_ks_create(ks_env, qs_env, error)
!!       Type(qs_ks_env_type), Pointer:: ks_env
!!       Type(qs_environment_type), Intent (IN), Target:: qs_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine qs_ks_create
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
    TYPE(qs_ks_env_type),POINTER :: ks_env
    TYPE(qs_environment_type), INTENT(in),TARGET :: qs_env
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    INTEGER, PARAMETER :: maxgridpoints = 100
    CHARACTER(len=*), PARAMETER :: routineN='qs_ks_create',&
         routineP=moduleN//':'//routineN
    LOGICAL :: failure
    INTEGER  :: stat, ispin, nspins, handle
    TYPE(dft_control_type), POINTER :: dft_control
    TYPE(pw_pool_type), POINTER :: auxbas_pw_pool
    TYPE(pw_env_type), POINTER :: pw_env
    TYPE(cell_type), POINTER :: cell

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.
    ALLOCATE(ks_env, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    IF (.not.failure) THEN
       NULLIFY( ks_env%v_rspace, ks_env%v, pw_env, auxbas_pw_pool,&
            pw_env, cell)
       CALL get_qs_env(qs_env=qs_env,&
            dft_control=dft_control,&
            pw_env=pw_env, cell=cell)
       CALL pw_env_get(pw_env,auxbas_pw_pool=auxbas_pw_pool)

       nspins=dft_control%nspins

       ks_env%s_mstruct_changed=.TRUE.
       ks_env%rho_changed=.TRUE.
       ks_env%forces_up_to_date=.FALSE.
       ks_env%n_evals=0
       ks_env%ref_count=1
       last_ks_id_nr=last_ks_id_nr+1
       ks_env%id_nr=last_ks_id_nr

       ALLOCATE(ks_env%v_rspace(nspins),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT.failure) THEN
       DO ispin=1,nspins
          CALL pw_pool_init_coeff(auxbas_pw_pool,ks_env%v_rspace(ispin),&
               use_data=REALDATA3D, in_space=REALSPACE,error=error)
       END DO
       CALL pw_pool_init_coeff(auxbas_pw_pool,ks_env%v_hartree_rspace,&
            use_data=REALDATA3D, in_space=REALSPACE)
    END IF

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE qs_ks_create
!***************************************************************************

!!****f* qs_ks_methods/qs_ks_build_kohn_sham_matrix [1.0] *
!!
!!   NAME
!!     qs_ks_build_kohn_sham_matrix
!!
!!   SYNOPSIS
!!     Subroutine qs_ks_build_kohn_sham_matrix(ks_env, qs_env, ks, rho,&
!!         energy, rebuild, calculate_forces, just_energy, error)
!!       Type(qs_ks_env_type), Pointer:: ks_env
!!       Type(qs_environment_type), Intent (INOUT):: qs_env
!!       Type(real_matrix_p_type), Dimension(:), Intent (INOUT):: ks
!!       Type(qs_rho_type), Pointer:: rho
!!       Type(qs_energy_type), Intent (INOUT), Target:: energy
!!       Logical, Intent (IN):: calculate_forces, rebuild, just_energy
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine qs_ks_build_kohn_sham_matrix
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
!!     - rebuild: if true builds again the potential
!!     - calculate_forces: if true calculate the quantities needed
!!       to calculate the forces. Defaults to false.
!!     - just_energy: if true updates the energies but not the
!!       ks matrix. Defaults to false
!!     - ks: the KS matrix that is created
!!     - rho: the density, must be up to date
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
!!
!!*** *********************************************************************
  SUBROUTINE qs_ks_build_kohn_sham_matrix(ks_env,qs_env,ks,&
       rho,energy,rebuild,calculate_forces,just_energy,error)
    TYPE(qs_ks_env_type), POINTER :: ks_env
    TYPE(qs_environment_type), INTENT(INOUT)  :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), INTENT(inout) :: ks
    TYPE(qs_rho_type), POINTER :: rho
    TYPE(qs_energy_type), INTENT(inout), TARGET :: energy
    LOGICAL, INTENT(in) :: calculate_forces,rebuild, just_energy
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

    !   *** Local variables ***

    REAL(wp)          :: total_rho_elec_rspace
    CHARACTER(len=*), PARAMETER :: routineN='build_kohn_sham',&  ! shorter for output
         routineP=moduleN//':'//routineN
    INTEGER :: handle,output_unit, nspins, ispin, stat, n_electrons,handle2
    LOGICAL :: ionode, failure,r_build
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: h,s
    TYPE(dft_control_type), POINTER :: dft_control
    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(pw_pool_p_type), DIMENSION(:), POINTER :: pw_pools
    TYPE(pw_pool_type), POINTER :: auxbas_pw_pool
    TYPE(pw_env_type), POINTER :: pw_env
    TYPE(coeff_type), POINTER :: rho_core
    TYPE(cell_type), POINTER :: cell
    TYPE(cp_logger_type), POINTER :: logger
    TYPE(qs_charges_type), POINTER :: qs_charges

    TYPE(coeff_type), DIMENSION(:), ALLOCATABLE   :: v_rspace_new
    TYPE(coeff_type)                              :: v_hartree_gspace,&
                                                     rho_tot_gspace
    TYPE(coeff_type), DIMENSION(:), POINTER       :: mgrid_gspace, &
                                                     mgrid_rspace

    failure=.FALSE.
    r_build=rebuild
    NULLIFY(h,s,dft_control,atomic_kind_set,auxbas_pw_pool, pw_pools,&
         pw_env,rho_core,cell, logger, mgrid_rspace, mgrid_gspace)
    !   -----------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)
    CALL timeset(routineN//"_energy","I","",handle2)

    CPPrecondition(ASSOCIATED(ks_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ks_env%ref_count>0,cp_failure_level,routineP,error,failure)

    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    output_unit = cp_logger_get_default_unit_nr(logger,local=.FALSE.)

    CALL get_qs_env(qs_env=qs_env,&
                    dft_control=dft_control,&
                    atomic_kind_set=atomic_kind_set,&
                    h=h,&
                    s=s,&
                    pw_env=pw_env,&
                    cell=cell,&
                    rho_core=rho_core,&
                    qs_charges=qs_charges)

    nspins=dft_control%nspins
    CPPrecondition(ASSOCIATED(h),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(s),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(rho),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(pw_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(SIZE(ks)>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(rho%rho_r_valid,cp_failure_level,routineP,error,failure)
    CPPrecondition(rho%rho_g_valid,cp_failure_level,routineP,error,failure)
    
    ! gets the tmp grids
    IF (.NOT. failure) THEN
       CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,&
            pw_pools=pw_pools)
       ALLOCATE(v_rspace_new(nspins), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(mgrid_rspace(SIZE(pw_pools)), mgrid_gspace(SIZE(pw_pools)),&
            stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT.failure) THEN
       CALL pw_pool_init_coeff(auxbas_pw_pool,v_hartree_gspace,&
               use_data = COMPLEXDATA1D,&
               in_space = RECIPROCALSPACE, error=error)
       CALL pw_pools_init_coeffs(pw_pools,mgrid_gspace,&
            use_data = COMPLEXDATA1D,&
            in_space = RECIPROCALSPACE, error=error)
       CALL pw_pools_init_coeffs(pw_pools,mgrid_rspace,&
            use_data = REALDATA3D,&
            in_space = REALSPACE, error=error)
    END IF

    IF (.NOT.failure) THEN
       IF (ks_env%s_mstruct_changed) THEN
          IF (ASSOCIATED(ks_env%v)) THEN
             DO ispin=1,SIZE(ks_env%v)
                IF (ASSOCIATED(ks_env%v(ispin)%matrix)) THEN
                   CALL deallocate_matrix(ks_env%v(ispin)%matrix)
                END IF
             END DO
             DEALLOCATE(ks_env%v,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          ALLOCATE(ks_env%v(nspins),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          IF (.NOT. failure) THEN
             DO ispin=1,nspins
                NULLIFY(ks_env%v(ispin)%matrix)
                CALL replicate_matrix_structure(h(1)%matrix,&
                     ks_env%v(ispin)%matrix,&
                     "POTENTIAL MATRIX "//ADJUSTL(cp_to_string(ispin)))
             END DO
          END IF
          r_build=.TRUE.
          ks_env%s_mstruct_changed=.FALSE.
       END IF
    END IF

! *** calculate the hartree potential on the total density ***
    IF (.NOT. failure) THEN
       CALL pw_pool_init_coeff(auxbas_pw_pool,rho_tot_gspace,&
            use_data = COMPLEXDATA1D,&
            in_space = RECIPROCALSPACE, error=error)

       CALL coeff_add(rho_core,rho%rho_g(1),rho_tot_gspace)
       DO ispin=2,nspins
          CALL coeff_sumup(rho%rho_g(ispin),rho_tot_gspace)
       END DO
       qs_charges%total_rho_gspace = calculate_total_rho(rho_tot_gspace)

       IF (ionode.AND.logger%print_keys%total_densities) THEN
          CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
               nelectron=n_electrons)
          n_electrons=n_electrons-dft_control%charge

          WRITE (UNIT=output_unit,FMT="(/,(T3,A,T40,2F20.10))")&
               "Total electronic density (r-space): ",&
               SUM(rho%tot_rho_r),&
               SUM(rho%tot_rho_r) + &
               REAL(n_electrons,wp),&
               "Total core charge density (r-space):",&
               qs_charges%total_rho_core_rspace,&
               qs_charges%total_rho_core_rspace - REAL(n_electrons,wp)
          WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
               "Total charge density (r-space):     ",&
               SUM(rho%tot_rho_r)+&
               qs_charges%total_rho_core_rspace,&
               "Total charge density (g-space):     ",&
               qs_charges%total_rho_gspace
       END IF

       CALL init_pw_poisson_solver(cell%green,error=error)
       CALL calculate_hartree(rho_tot_gspace,energy%hartree,&
            v_hartree_gspace)
       CALL cleanup_pw_poisson_solver(cell%green,error=error)

       CALL pw_pool_give_back_coeff(auxbas_pw_pool,rho_tot_gspace,&
            error=error)

! *** calculate the xc potential ***       

       DO ispin=1,nspins
          CALL pw_pool_init_coeff(auxbas_pw_pool,v_rspace_new(ispin),&
               use_data = REALDATA3D,&
               in_space = REALSPACE, error=error)
          CALL coeff_zero(v_rspace_new(ispin))
       END DO

       CALL calculate_xc_potential(qs_env=qs_env,&
            v_rspace_new=v_rspace_new,&
            error=error)

       energy%total = energy%core_overlap + energy%core_self + &
            energy%core + energy%hartree + energy%exc

       IF (ionode.AND.logger%print_keys%scf_energies) THEN
          WRITE (UNIT=output_unit,FMT="(/,(T3,A,T60,F20.10))")&
            "Core Hamiltonian energy:                       ",energy%core,&
            "Hartree energy:                                ",energy%hartree,&
            "Exchange-correlation energy:                   ",energy%exc
       END IF
       CALL timestop(0.0_wp,handle2)

       IF (.NOT. just_energy) THEN
          CALL timeset(routineN//"_pot","I","",handle2)

          CALL coeff_transform_space(v_hartree_gspace,&
               ks_env%v_hartree_rspace)
          ks_env%v_hartree_rspace%pw%cr3d =&
               ks_env%v_hartree_rspace%pw%pw_grid%dvol*&
               ks_env%v_hartree_rspace%pw%cr3d

          IF (calculate_forces) THEN
             CALL integrate_v_core_rspace(ks_env%v_hartree_rspace,&
                  pw_env%cube_info(1),&
                  pw_env%l_info,qs_env)
          END IF

          DO ispin=1,nspins
             v_rspace_new(ispin)%pw%cr3d =&
                  v_rspace_new(ispin)%pw%pw_grid%dvol*&
                  v_rspace_new(ispin)%pw%cr3d +&
                  ks_env%v_hartree_rspace%pw%cr3d
             
             IF (r_build) CALL set_matrix(ks_env%v(ispin)%matrix,0.0_wp)
             CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin),&
                  mgrid_rspace=mgrid_rspace,&
                  mgrid_gspace=mgrid_gspace,&
                  auxbas_grid=pw_env%auxbas_grid,&
                  gridlevel_info=pw_env%gridlevel_info,&
                  v_rspace_old=ks_env%v_rspace(ispin),&
                  cube_info=pw_env%cube_info,&
                  l_info=pw_env%l_info,&
                  rebuild=r_build,p=rho%rho_ao(ispin),h=ks_env%v(ispin),&
                  qs_env=qs_env,&
                  calculate_forces=calculate_forces,error=error)
             CALL add_matrices(ks(ispin)%matrix,1.0_wp,h(1)%matrix,1.0_wp,&
                  ks_env%v(ispin)%matrix)
          END DO
          CALL timestop(0.0_wp,handle2)
       END IF

       CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_hartree_gspace,&
            error=error)
       DO ispin=1,nspins
          CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_rspace_new(ispin),&
               error=error)
       END DO
       CALL pw_pools_give_back_coeffs(pw_pools,mgrid_gspace,&
            error=error)
       CALL pw_pools_give_back_coeffs(pw_pools,mgrid_rspace,&
            error=error)
       DEALLOCATE(v_rspace_new, stat=stat) ! just to be sure...
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       DEALLOCATE(mgrid_rspace, mgrid_gspace, stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

    END IF

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE qs_ks_build_kohn_sham_matrix
! *****************************************************************************

!!****f* qs_ks_methods/qs_ks_update_qs_env [1.0] *
!!
!!   NAME
!!     qs_ks_update_qs_env
!!
!!   SYNOPSIS
!!     Subroutine qs_ks_update_qs_env(ks_env, qs_env, rebuild,&
!!         calculate_forces, just_energy, error)
!!       Type(qs_ks_env_type), Pointer:: ks_env
!!       Type(qs_environment_type), Intent (INOUT), Target:: qs_env
!!       Logical, Intent (IN), Optional:: rebuild, calculate_forces,&
!!         just_energy
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine qs_ks_update_qs_env
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
!!     - rebuild: if true builds again th
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
SUBROUTINE qs_ks_update_qs_env(ks_env,qs_env,&
     rebuild,calculate_forces, just_energy, error)
  TYPE(qs_ks_env_type), POINTER :: ks_env
  TYPE(qs_environment_type), INTENT(inout), TARGET :: qs_env
  LOGICAL, INTENT(in), OPTIONAL :: rebuild, calculate_forces, just_energy
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure, r_build, c_forces, energy_only, did_update
  CHARACTER(len=*), PARAMETER :: routineN='qs_ks_update_qs_env',&
        routineP=moduleN//':'//routineN
  INTEGER :: handle,ispin
  REAL(kind=wp) :: ecore
  TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: ks,p, h, s
  TYPE(qs_energy_type), POINTER :: energy
  TYPE(cp_logger_type), POINTER :: logger
  TYPE(qs_rho_type), POINTER :: rho
  TYPE(dft_control_type), POINTER :: dft_control
  TYPE(cp_para_env_type), POINTER :: para_env

  CALL timeset(routineN,"I","",handle)

  failure=.FALSE.
  did_update=.FALSE.
  NULLIFY(ks,h,s,energy, logger, rho, dft_control, para_env)
  logger => cp_error_get_logger(error)
  r_build=.FALSE.; c_forces=.FALSE.; energy_only = .FALSE.
  IF (PRESENT(just_energy)) energy_only=just_energy
  IF (PRESENT(calculate_forces)) c_forces=calculate_forces
  IF (PRESENT(rebuild)) r_build=rebuild

  CPPrecondition(ASSOCIATED(ks_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(ks_env%ref_count>0,cp_failure_level,routineP,error,failure)

  IF (.NOT. failure) THEN
     IF (ks_env%rho_changed .OR. ks_env%s_mstruct_changed .OR.&
          (c_forces.AND..NOT.ks_env%forces_up_to_date)) THEN
        did_update=.TRUE.
        CALL get_qs_env(qs_env,k=ks,rho=rho,energy=energy,&
             dft_control=dft_control,&
             h=h,para_env=para_env,s=s,error=error)

        ! update ecore
        CALL cp_assert(rho%rho_ao_valid,cp_warning_level,&
             cp_assertion_failed, routineP,&
             "rho_ao  is not valid, thus energy%core (and total) might "//&
             "be wrong "//&
             CPSourceFileRef,&
             error=error)
        energy%core=0.0_wp
        DO ispin=1,dft_control%nspins
          CALL calculate_ecore(h=h(1)%matrix,&
                               p=rho%rho_ao(ispin)%matrix,&
                               ecore=ecore,&
                               para_env=para_env)
          energy%core=energy%core+ecore
        END DO

        CALL qs_ks_build_kohn_sham_matrix(ks_env,qs_env=qs_env,&
             ks=ks,&
             rho=rho,energy=energy,rebuild=r_build,&
             calculate_forces=c_forces,just_energy=energy_only,&
             error=error)
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

  CALL timestop(0.0_wp,handle)

END SUBROUTINE qs_ks_update_qs_env
!***************************************************************************

!!****f* qs_ks_methods/qs_ks_did_change [1.0] *
!!
!!   NAME
!!     qs_ks_did_change
!!
!!   SYNOPSIS
!!     Subroutine qs_ks_did_change(ks_env, qs_env, s_mstruct_changed,&
!!         rho_changed, full_reset, error)
!!       Type(qs_ks_env_type), Pointer:: ks_env
!!       Type(qs_environment_type), Intent (INOUT), Target:: qs_env
!!       Logical, Intent (IN), Optional:: s_mstruct_changed, rho_changed,&
!!         full_reset
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine qs_ks_did_change
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
!!     - qs_env: the qs_environement
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
SUBROUTINE qs_ks_did_change(ks_env,qs_env,s_mstruct_changed,&
     rho_changed,full_reset,error)
  TYPE(qs_ks_env_type), POINTER :: ks_env
  TYPE(qs_environment_type), INTENT(inout), TARGET :: qs_env
  LOGICAL, INTENT(in), OPTIONAL :: s_mstruct_changed, rho_changed, full_reset
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='qs_ks_did_change',&
        routineP=moduleN//':'//routineN
  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(ks_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(ks_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     IF (PRESENT(rho_changed)) THEN
        IF (rho_changed) ks_env%rho_changed=.TRUE.
     END IF
     IF (PRESENT(s_mstruct_changed)) THEN
        IF (s_mstruct_changed) THEN
           ks_env%s_mstruct_changed=.TRUE.
        END IF
     END IF
     IF (PRESENT(full_reset)) THEN
        IF (full_reset) THEN
           ks_env%s_mstruct_changed=.TRUE.
        END IF
     END IF
  END IF
END SUBROUTINE qs_ks_did_change
!***************************************************************************

!!****f* qs_ks_methods/calculate_xc_potential [1.0] *
!!
!!   NAME
!!     calculate_xc_potential
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
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
!!
!!*** *********************************************************************
  SUBROUTINE calculate_xc_potential(v_rspace_new,qs_env,error)
    IMPLICIT NONE
    TYPE(qs_environment_type), INTENT(in), TARGET :: qs_env
    TYPE(coeff_type), DIMENSION(:), INTENT(inout) :: v_rspace_new
    TYPE(cp_error_type), INTENT(inout), OPTIONAL  :: error

!   *** Local variables ***

    INTEGER :: handle

    INTEGER, DIMENSION(3) :: ng

    REAL(wp), DIMENSION(:,:,:), POINTER :: drho,rho,vxc
    TYPE(dft_control_type), POINTER :: dft_control
    TYPE(qs_rho_type), POINTER :: rho_struct
    TYPE(qs_energy_type), POINTER :: energy
    TYPE(pw_env_type), POINTER :: pw_env
    TYPE(pw_pool_type), POINTER :: auxbas_pw_pool
    TYPE(cell_type), POINTER :: cell

    LOGICAL :: failure
    INTEGER :: i
    CHARACTER(len=*), PARAMETER :: routineN='calculate_xc_potential',&
         routineP=moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    CALL timeset( routineN ,"I","",handle)

    failure=.FALSE.
    NULLIFY(drho,rho,vxc,dft_control,rho_struct,energy,pw_env,auxbas_pw_pool,&
         cell)

    CALL get_qs_env(qs_env=qs_env,&
                    dft_control=dft_control,&
                    rho=rho_struct,&
                    pw_env=pw_env,&
                    energy=energy,&
                    cell=cell)

    ! test if the real space density is available
    CPPrecondition(ASSOCIATED(rho_struct),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_struct%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_struct%rho_r_valid,cp_failure_level,routineP,error,failure)
    IF ( .NOT. failure ) THEN

      ! for gradient corrected functional we need the density in g space
      IF ( ANY(qs_env%dft_control%gradient_functionals) ) THEN
         IF ( .NOT. rho_struct%rho_g_valid ) THEN
            DO i=1,SIZE(rho_struct%rho_r)
               CALL coeff_transform_space(rho_struct%rho_r(i),rho_struct%rho_g(i))
            END DO
         END IF
      END IF

      IF ( dft_control%nspins == 1 .OR. dft_control%nspins == 2) THEN
         CALL pw_env_get(pw_env,auxbas_pw_pool=auxbas_pw_pool,error=error)
         CALL xc_calculate_pw2(rho_r=rho_struct%rho_r,&
              rho_g=rho_struct%rho_g, vxc=v_rspace_new,exc=energy%exc,&
              dft_control=dft_control,&
              cell=cell, pw_pool=auxbas_pw_pool,error=error)
      ELSE
         CALL stop_program(routineP,"nspins > 2")
      END IF

   END IF

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE calculate_xc_potential

! *****************************************************************************

END MODULE qs_ks_methods

! *****************************************************************************
