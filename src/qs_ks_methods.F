!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002,2003,2004  CP2K developers group                       !
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
                                             get_atomic_kind_set
  USE coefficient_types,               ONLY: coeff_add,&
                                             coeff_init,&
                                             coeff_sumup,&
                                             coeff_transform_space,&
                                             coeff_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsystem_types,              ONLY: cp_subsystem_p_type
  USE hartree,                         ONLY: calculate_hartree
  USE hartree_local_methods,           ONLY: Vh_1c_gg_integrals
  USE kinds,                           ONLY: dp
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_give_back_coeff,&
                                             pw_pool_init_coeff,&
                                             pw_pool_p_type,&
                                             pw_pool_retain,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type
  USE pws,                             ONLY: cleanup_pw_poisson_solver,&
                                             init_pw_poisson_solver
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
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_ks_qmmm_types,                ONLY: qs_ks_qmmm_env_type
  USE qs_rho0_ggrid,                   ONLY: integrate_vhg0_rspace
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE qs_vxc_atom,                     ONLY: calculate_vxc_atom
  USE simulation_cell,                 ONLY: cell_type
  USE sparse_matrix_types,             ONLY: real_matrix_p_type,&
                                             transfer_matrix
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc,                              ONLY: xc_exc_calc,&
                                             xc_vxc_pw_create1
  USE qmmm_methods,                    ONLY: qmmm_calculate_energy,&
                                             qmmm_modify_hartree_pot
  USE qs_ks_qmmm_methods,              ONLY: qs_ks_qmmm_create
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
    TYPE(qs_environment_type), INTENT(in), &
      TARGET                                 :: qs_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ks_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nspins, stat
    LOGICAL                                  :: failure
    TYPE(cell_type), POINTER                 :: cell
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.
    CPPreconditionNoFail(.NOT.ASSOCIATED(ks_env),cp_failure_level,routineP,error)
    ALLOCATE(ks_env, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    IF (.not.failure) THEN
       NULLIFY( pw_env, auxbas_pw_pool,&
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
       rho,energy,calculate_forces,just_energy,error)
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_environment_type), INTENT(INOUT) :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      INTENT(inout)                          :: ks
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(qs_energy_type), INTENT(inout), &
      TARGET                                 :: energy
    LOGICAL, INTENT(in)                      :: calculate_forces, just_energy
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_kohn_sham', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, n_electrons, &
                                                nspins, output_unit, stat
    LOGICAL                                  :: failure, gapw, ionode
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(coeff_type)                         :: rho_tot_gspace, &
                                                v_hartree_gspace
    TYPE(coeff_type), DIMENSION(:), POINTER  :: v_rspace_new, v_tau_rspace
    TYPE(coeff_type), POINTER                :: rho_core,&
                                                rho0_s_gs
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_charges_type), POINTER           :: qs_charges
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: h, s
    REAL(dp)                                 :: tot1_h, tot1_s 

! shorter for output

    failure=.FALSE.
    NULLIFY(h,s,dft_control,atomic_kind_set,auxbas_pw_pool, pw_pools,&
         pw_env,rho_core,rho0_s_gs,cell, logger, v_rspace_new, v_tau_rspace)

    !   -----------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

    CPPrecondition(ASSOCIATED(ks_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ks_env%ref_count>0,cp_failure_level,routineP,error,failure)

    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    IF (ionode) output_unit = cp_logger_get_default_unit_nr(logger)

    CALL get_qs_env(qs_env=qs_env,&
         dft_control=dft_control,&
         atomic_kind_set=atomic_kind_set,&
         h=h,&
         s=s,&
         pw_env=pw_env,&
         cell=cell,&
         rho_core=rho_core,&
         rho0_s_gs=rho0_s_gs,&
         qs_charges=qs_charges)

    nspins=dft_control%nspins
    CPPrecondition(ASSOCIATED(h),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(s),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(rho),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(pw_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(SIZE(ks)>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(rho%rho_r_valid,cp_failure_level,routineP,error,failure)
    CPPrecondition(rho%rho_g_valid,cp_failure_level,routineP,error,failure)

    ! gapw method 
    gapw = dft_control%qs_control%gapw

    ! gets the tmp grids
    IF (.NOT. failure) THEN
       CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,&
            pw_pools=pw_pools)
    END IF
    IF (.NOT.failure) THEN
       CALL pw_pool_init_coeff(auxbas_pw_pool,v_hartree_gspace,&
            use_data = COMPLEXDATA1D,&
            in_space = RECIPROCALSPACE, error=error)
    END IF

! ***  Prepare densities for gapw ***
    IF(gapw) THEN
      CALL prepare_gapw_den(qs_env,error)
    ENDIF

! *** calculate the hartree potential on the pw density ***
    IF (.NOT. failure) THEN
       CALL pw_pool_init_coeff(auxbas_pw_pool,rho_tot_gspace,&
            use_data = COMPLEXDATA1D,&
            in_space = RECIPROCALSPACE, error=error)

       IF(gapw) THEN
         CPPrecondition(ASSOCIATED(rho0_s_gs),cp_failure_level,routineP,error,failure)
         CALL coeff_add(rho0_s_gs,rho%rho_g(1),rho_tot_gspace)
       ELSE
         CALL coeff_add(rho_core,rho%rho_g(1),rho_tot_gspace)
       END IF
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
               REAL(n_electrons,dp),&
               "Total core charge density (r-space):",&
               qs_charges%total_rho_core_rspace,&
               qs_charges%total_rho_core_rspace - REAL(n_electrons,dp)
          IF(gapw) THEN
            tot1_h =  qs_charges%total_rho1_hard(1)
            tot1_s =  qs_charges%total_rho1_soft(1)
            DO ispin=2,nspins
              tot1_h = tot1_h + qs_charges%total_rho1_hard(ispin)
              tot1_s = tot1_s + qs_charges%total_rho1_soft(ispin)
            END DO
            WRITE (UNIT=output_unit,FMT="(/,(T3,A,T44,2F18.10))")&
                "Hard and soft local densities (lebedev):",&
                tot1_h, tot1_s
            WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
               "Total Rho_soft + Rho1_hard - Rho1_soft (r-space): ",&
               SUM(rho%tot_rho_r)+ tot1_h - tot1_s ,&
               "Total charge density (r-space):      ",&
               SUM(rho%tot_rho_r)+ tot1_h - tot1_s &
               + qs_charges%total_rho_core_rspace,&
               "Total Rho_soft + Rho0_soft (g-space):",&
               qs_charges%total_rho_gspace
          ELSE  
            WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
               "Total charge density (r-space):     ",&
               SUM(rho%tot_rho_r)+&
               qs_charges%total_rho_core_rspace,&
               "Total charge density (g-space):     ",&
               qs_charges%total_rho_gspace
          END IF

       END IF

       CALL init_pw_poisson_solver(cell%green,error=error)
       CALL calculate_hartree(rho_tot_gspace,energy%hartree,&
            v_hartree_gspace)
       CALL cleanup_pw_poisson_solver(cell%green,error=error)
       CALL pw_pool_give_back_coeff(auxbas_pw_pool,rho_tot_gspace,&
            error=error)
       IF (.NOT. just_energy) THEN
          CALL coeff_transform_space(v_hartree_gspace,&
               ks_env%v_hartree_rspace)
          ks_env%v_hartree_rspace%pw%cr3d =&
               ks_env%v_hartree_rspace%pw%pw_grid%dvol*&
               ks_env%v_hartree_rspace%pw%cr3d
       ENDIF
       CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_hartree_gspace,&
            error=error)

!      WRITE(*,*)  energy%hartree

       ! *** Add the QM/MM potential ***
       If (qs_env%qmmm) Then
          WRITE (UNIT=output_unit,FMT="(T3,A)")&
               "Adding QM/MM Electrostatic Potential to the Kohn-Sham Matrix.."
           CALL qmmm_calculate_energy (qs_env=qs_env,&
                                       rho=rho%rho_r,&
                                       v_qmmm=qs_env%ks_qmmm_env%v_qmmm_rspace,&
                                       qmmm_energy=energy%qmmm_el,&
                                       error=error)
           CALL qmmm_modify_hartree_pot(v_hartree=ks_env%v_hartree_rspace,&
                                        v_qmmm=qs_env%ks_qmmm_env%v_qmmm_rspace,&
                                        error=error)
       End If

! *** calculate the xc potential on the pw density ***       

       CALL qs_vxc_create(qs_env=qs_env,&
            vxc_rho=v_rspace_new, vxc_tau=v_tau_rspace, &
            calculate_forces=calculate_forces, &
            just_energy=just_energy, error=error)

! ***  Single atom contributions ***
       IF(gapw) THEN
         CALL calculate_vxc_atom(qs_env,error)
         CALL Vh_1c_gg_integrals(qs_env,energy)
       ENDIF
!       WRITE(*,*) energy%hartree_1c
!      stop 'hartree'
       energy%total = energy%core_overlap + energy%core_self + &
            energy%core + energy%hartree + energy%hartree_1c +&
            energy%exc + energy%exc1 + &
            energy%qmmm_el

       IF (ionode.AND.logger%print_keys%scf_energies) THEN
          WRITE (UNIT=output_unit,FMT="(/,(T3,A,T60,F20.10))")&
               "Core Hamiltonian energy:                       ",energy%core,&
               "Hartree energy:                                ",energy%hartree,&
               "Exchange-correlation energy:                   ",energy%exc
          IF(gapw) THEN
             WRITE (UNIT=output_unit,FMT="(/,(T3,A,T60,F20.10))")&
                "GAPW| Exc from hard and soft atomic rho1:      ",energy%exc1,&
                "GAPW| local Eh = 1 center integrals:           ",energy%hartree_1c
          END IF
          WRITE(*,*)'QM/MM enabled:', qs_env%qmmm
          If  (qs_env%qmmm) Then
             WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
               "QM/MM Electrostatic energy:                    ",energy%qmmm_el
          End If
       END IF

       IF (.NOT. just_energy) THEN
          IF (calculate_forces .AND. .NOT. gapw) THEN
             CALL integrate_v_core_rspace(ks_env%v_hartree_rspace, qs_env)
          END IF

          DO ispin=1,nspins
             ! copy the core ham into the ks matrix
             CALL transfer_matrix(h(1)%matrix,ks(ispin)%matrix)
          END DO
          IF (ASSOCIATED(v_rspace_new)) THEN
             DO ispin=1,nspins
                v_rspace_new(ispin)%pw%cr3d =&
                     v_rspace_new(ispin)%pw%pw_grid%dvol*&
                     v_rspace_new(ispin)%pw%cr3d +&
                     ks_env%v_hartree_rspace%pw%cr3d

                ! add the part due to v_rspace
                CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin),&
                     p=rho%rho_ao(ispin),h=ks(ispin),&
                     qs_env=qs_env, &
                     calculate_forces=calculate_forces,&
                     gapw=gapw,error=error)

                CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_rspace_new(ispin),&
                     error=error)
             END DO
             DEALLOCATE(v_rspace_new,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             
          ELSE
             DO ispin=1,nspins
                ! add oly v_hartree
                CALL integrate_v_rspace(v_rspace=ks_env%v_hartree_rspace,&
                     p=rho%rho_ao(ispin),h=ks(ispin), qs_env=qs_env,&
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
                      p=rho%rho_ao(ispin),h=ks(ispin),&
                      qs_env=qs_env,&
                      calculate_forces=calculate_forces,compute_tau=.TRUE., &
                      gapw=gapw,error=error)

              CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_tau_rspace(ispin),&
                   error=error)                 

              END DO
              DEALLOCATE(v_tau_rspace, stat=stat)
              CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ENDIF

          IF(gapw) THEN
!         ***  Integrals of the Hartree potential with g0_soft ***
             CALL integrate_vhg0_rspace(qs_env,ks_env%v_hartree_rspace, &
                                        calculate_forces, error=error)

!         ***  Single atom contributions in the KS matrix ***
             DO ispin = 1,nspins
                CALL update_ks_atom(qs_env,ks(ispin),rho%rho_ao(ispin),&
                     ispin,calculate_forces,error)
             ENDDO
          ENDIF

       END IF

    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE qs_ks_build_kohn_sham_matrix
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
SUBROUTINE qs_ks_update_qs_env(ks_env,qs_env,&
     calculate_forces, just_energy, error)
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_environment_type), &
      INTENT(inout), TARGET                  :: qs_env
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
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: h, ks

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
  NULLIFY(ks,h,energy, logger, rho, dft_control, para_env)
  logger => cp_error_get_logger(error)

  CPPrecondition(ASSOCIATED(ks_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(ks_env%ref_count>0,cp_failure_level,routineP,error,failure)

  IF (.NOT. failure) THEN
     IF (ks_env%rho_changed .OR. ks_env%s_mstruct_changed .OR.&
          (c_forces.AND..NOT.ks_env%forces_up_to_date)) THEN
        did_update=.TRUE.
        CALL get_qs_env(qs_env,k=ks,rho=rho,&
             energy=energy,&
             dft_control=dft_control,&
             h=h,para_env=para_env,error=error)

        ! update ecore
        energy%core=0.0_dp
        DO ispin=1,dft_control%nspins
          CALL calculate_ecore(h=h(1)%matrix,&
                               p=rho%rho_ao(ispin)%matrix,&
                               ecore=ecore,&
                               para_env=para_env)
          energy%core=energy%core+ecore
        END DO

        CALL qs_ks_build_kohn_sham_matrix(ks_env,qs_env=qs_env,&
             ks=ks,&
             rho=rho,&
             energy=energy,&
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
SUBROUTINE qs_ks_did_change(ks_env,s_mstruct_changed,rho_changed,full_reset,error)
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    LOGICAL, INTENT(in), OPTIONAL            :: s_mstruct_changed, &
                                                rho_changed, full_reset
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
!!
!!*** *********************************************************************
SUBROUTINE qs_vxc_create(vxc_rho,vxc_tau,qs_env,&
     calculate_forces, just_energy, error)
    TYPE(coeff_type), DIMENSION(:), POINTER  :: vxc_rho, vxc_tau
    TYPE(qs_environment_type), INTENT(in), &
      TARGET                                 :: qs_env
    LOGICAL, INTENT(in), OPTIONAL            :: calculate_forces, just_energy
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_vxc_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, stat
    LOGICAL                                  :: failure, my_just_energy
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsystem_p_type), DIMENSION(1)  :: subsys
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: my_vxc_rho, my_vxc_tau, &
                                                rho_g, rho_r, tau
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_rho_type), POINTER               :: rho_struct

!   ---------------------------------------------------------------------------

  CALL timeset( routineN ,"I","",handle)

  failure=.FALSE.
  CPPrecondition(.NOT.ASSOCIATED(vxc_rho),cp_failure_level,routineP,error,failure)
  CPPrecondition(.NOT.ASSOCIATED(vxc_tau),cp_failure_level,routineP,error,failure)
  NULLIFY(dft_control,rho_struct,energy,pw_env,auxbas_pw_pool,&
       cell, subsys(1)%subsys,my_vxc_rho, my_vxc_tau, rho_g, rho_r, tau)
  my_just_energy=.FALSE.
  IF (PRESENT(just_energy)) my_just_energy=just_energy

  CALL get_qs_env(qs_env=qs_env,&
       dft_control=dft_control,&
       rho=rho_struct,&
       pw_env=pw_env,&
       energy=energy,&
       subsys=subsys(1)%subsys,&
       cell=cell)

  ! test if the real space density is available
  CPPrecondition(ASSOCIATED(rho_struct),cp_failure_level,routineP,error,failure)
  CPPrecondition(rho_struct%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(rho_struct%rho_r_valid,cp_failure_level,routineP,error,failure)
  CALL cp_assert( dft_control%nspins == 1 .OR. dft_control%nspins == 2,&
       cp_failure_level,cp_assertion_failed,routineP,&
       "nspins must be 1 or 2",error,failure)

  IF ( .NOT. failure ) THEN
     ALLOCATE(rho_r(dft_control%nspins),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
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

     CALL pw_env_get(pw_env,auxbas_pw_pool=auxbas_pw_pool,error=error)

     IF (my_just_energy) THEN
        energy%exc=xc_exc_calc(rho_r=rho_r,tau=tau,&
             rho_g=rho_g, dft_control=dft_control,&
             cell=cell, pw_pool=auxbas_pw_pool,atoms=subsys,&
             error=error)
     ELSE
        CALL xc_vxc_pw_create1(vxc_rho=my_vxc_rho,vxc_tau=my_vxc_tau, rho_r=rho_r,&
             rho_g=rho_g,tau=tau,exc=energy%exc,&
             dft_control=dft_control,&
             cell=cell, pw_pool=auxbas_pw_pool,atoms=subsys,&
             calculate_forces=calculate_forces,error=error)
     END IF
     DEALLOCATE(rho_r,stat=stat)
     CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     IF (ASSOCIATED(rho_g)) THEN
        DEALLOCATE(rho_g,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF

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
     END IF

  END IF

  CALL timestop(0.0_dp,handle)

END SUBROUTINE qs_vxc_create

! *****************************************************************************

END MODULE qs_ks_methods

! *****************************************************************************
