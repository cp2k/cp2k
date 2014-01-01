!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation of the Fock matrix for SE methods
!> \author JGH and TLAINO
!> \par History
!>      Teodoro Laino (04.2008) [tlaino] - University of Zurich : d-orbitals
!>      Teodoro Laino (09.2008) [tlaino] - University of Zurich : Speed-up
!>      Teodoro Laino (09.2008) [tlaino] - University of Zurich : Periodic SE
!>      Teodoro Laino (05.2009) [tlaino] - Split and module reorganization
! *****************************************************************************
MODULE se_fock_matrix_mpi
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE atprop_types,                    ONLY: atprop_array_init
  USE cp_control_types,                ONLY: dft_control_type,&
                                             semi_empirical_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_add,&
                                             cp_dbcsr_copy,&
                                             cp_dbcsr_get_info,&
                                             cp_dbcsr_multiply,&
                                             cp_dbcsr_trace
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_p_type,&
                                             cp_dbcsr_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE input_constants,                 ONLY: do_method_am1,&
                                             do_method_mndo,&
                                             do_method_mndod,&
                                             do_method_pdg,&
                                             do_method_pm3,&
                                             do_method_pm6,&
                                             do_method_pnnl,&
                                             do_method_rm1
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: particle_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE scp_delta_fock_matrix,           ONLY: build_delta_fock_matrix
  USE scp_energy_types,                ONLY: scp_energy_type
  USE scp_environment_types,           ONLY: get_scp_env,&
                                             scp_environment_type
  USE se_fock_matrix_coulomb_mpi,      ONLY: build_fock_matrix_coul_lr_r3,&
                                             build_fock_matrix_coulomb,&
                                             build_fock_matrix_coulomb_lr
  USE se_fock_matrix_dbg,              ONLY: dbg_energy_coulomb_lr
  USE se_fock_matrix_exchange,         ONLY: build_fock_matrix_exchange
  USE semi_empirical_store_int_types,  ONLY: semi_empirical_si_finalize,&
                                             semi_empirical_si_initialize,&
                                             semi_empirical_si_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'se_fock_matrix_mpi'
  LOGICAL, PARAMETER, PRIVATE          :: debug_this_module       = .FALSE.
  LOGICAL, PARAMETER, PRIVATE          :: debug_energy_coulomb_lr = .FALSE.

  PUBLIC :: build_se_fock_matrix

CONTAINS

! *****************************************************************************
!> \brief Construction of the Fock matrix for NDDO methods
!> \author JGH
!> \par History
!>         - Teodoro Laino [tlaino] (05.2009) - Split and module reorganization
! *****************************************************************************
  SUBROUTINE build_se_fock_matrix (ks_env,qs_env,ks_matrix,rho,energy,&
       calculate_forces,just_energy,error)

    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(qs_energy_type), POINTER            :: energy
    LOGICAL, INTENT(in)                      :: calculate_forces, just_energy
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_se_fock_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, natom, &
                                                ncol_global, nspins, &
                                                output_unit
    LOGICAL                                  :: failure, scp_nddo
    REAL(KIND=dp)                            :: ecoul, qmmm_el
    REAL(KIND=dp), DIMENSION(:), POINTER     :: occupation_numbers
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_h, matrix_p, mo_derivs
    TYPE(cp_dbcsr_type), POINTER             :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(scp_energy_type), POINTER           :: energy_scp
    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(section_vals_type), POINTER         :: scf_section
    TYPE(semi_empirical_control_type), &
      POINTER                                :: se_control
    TYPE(semi_empirical_si_type), POINTER    :: scp_store_int_env, &
                                                store_int_env

    failure=.FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(matrix_h, dft_control, logger, scf_section, store_int_env, se_control)
    NULLIFY(atomic_kind_set, scp_env, energy_scp, scp_store_int_env)
    logger => cp_error_get_logger(error)
    CPPrecondition(ASSOCIATED(ks_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ks_env%ref_count>0,cp_failure_level,routineP,error,failure)

    CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, matrix_h=matrix_h,&
         para_env=para_env, se_store_int_env=store_int_env, &
         atomic_kind_set=atomic_kind_set, &
         error=error)

    SELECT CASE (dft_control%qs_control%method_id)
    CASE DEFAULT
       ! Abort if the parameterization is an unknown one..
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="Fock Matrix not available for the chosen parameterization! ", &
            error=error, error_level=cp_failure_level)

    CASE (do_method_am1,do_method_rm1,do_method_mndo,do_method_pdg,&
         do_method_pm3,do_method_pm6,do_method_mndod, do_method_pnnl)

       ! Check for properly allocation of Matrixes
       nspins=dft_control%nspins
       CPPrecondition(((nspins>=1).AND.(nspins<=2)),cp_failure_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(matrix_h),cp_failure_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(rho),cp_failure_level,routineP,error,failure)
       CPPrecondition(SIZE(ks_matrix)>0,cp_failure_level,routineP,error,failure)

       se_control  => dft_control%qs_control%se_control
       scp_nddo    =  se_control%scp
       scf_section => section_vals_get_subs_vals(qs_env%input,"DFT%SCF",error=error)
       matrix_p => rho%rho_ao

       IF ( scp_nddo ) THEN
         CALL get_qs_env(qs_env=qs_env, scp_env=scp_env, error=error)
         CALL get_scp_env(scp_env=scp_env, se_store_int_env=scp_store_int_env, &
                          energy=energy_scp, error=error )
       END IF

       energy%qmmm_el = 0.0_dp
       energy%total = 0.0_dp

       DO ispin=1,nspins
          ! Copy the core matrix into the fock matrix
          CALL cp_dbcsr_copy(ks_matrix(ispin)%matrix,matrix_h(1)%matrix,error=error)
       END DO

!      WRITE ( *, * ) 'KS_ENV%s_mstruct_changed', ks_env%s_mstruct_changed
       IF(qs_env%atprop%energy) THEN
         CALL get_qs_env(qs_env=qs_env,particle_set=particle_set,error=error)
         natom = SIZE(particle_set)
         CALL atprop_array_init(qs_env%atprop%atecoul,natom,error)
       END IF

       ! Compute Exchange and Coulomb terms
       CALL semi_empirical_si_initialize(store_int_env, ks_env%s_mstruct_changed, error)
       CALL build_fock_matrix_exchange(qs_env,ks_matrix,matrix_p,calculate_forces,&
            store_int_env,error)
       CALL build_fock_matrix_coulomb(qs_env,ks_matrix,matrix_p,energy,calculate_forces,&
            store_int_env,error)

! Add the SCP energies and additional contribution to fock matrix if requested
       IF ( scp_nddo ) THEN
!         WRITE ( *, * ) ' SCP KS_ENV%s_mstruct_changed', ks_env%s_mstruct_changed
          CALL semi_empirical_si_initialize(scp_store_int_env, ks_env%s_mstruct_changed, error)
          CALL build_delta_fock_matrix (qs_env, ks_matrix,matrix_p,energy,&
               calculate_forces,just_energy,scp_store_int_env,error)
          CALL semi_empirical_si_finalize(scp_store_int_env, ks_env%s_mstruct_changed, error)
       END IF

       ! Debug statements for Long-Range
       IF (debug_energy_coulomb_lr.AND.se_control%do_ewald) THEN
          CALL dbg_energy_coulomb_lr(energy, ks_matrix, nspins, qs_env, matrix_p,&
                                     calculate_forces, store_int_env, error)
       END IF

       ! Long Range Electrostatic
       IF (se_control%do_ewald) THEN
          ! Evaluate Coulomb Long-Range
          CALL build_fock_matrix_coulomb_lr(qs_env,ks_matrix,matrix_p,energy,calculate_forces,&
               store_int_env, error)

          ! Possibly handle the slowly convergent term 1/R^3
          IF (se_control%do_ewald_r3) THEN
            CALL build_fock_matrix_coul_lr_r3(qs_env,ks_matrix,matrix_p,energy,&
                 calculate_forces,error)
          END IF
       END IF
       CALL semi_empirical_si_finalize(store_int_env, ks_env%s_mstruct_changed, error)

       IF(qs_env%atprop%energy) THEN
         qs_env%atprop%atecoul = 0.5_dp*qs_env%atprop%atecoul
       END IF

       ! Compute the Hartree energy
       ! NOTE: If we are performing SCP-NDDO, ks_matrix contains coulomb piece from SCP.
       DO ispin=1,nspins
          CALL cp_dbcsr_trace(ks_matrix(ispin)%matrix,matrix_p(ispin)%matrix,trace=ecoul,error=error)
          energy%hartree = energy%hartree + ecoul
       END DO
!      WRITE ( *, * ) 'AFTER Hartree',  ecoul, energy%hartree


!       CALL build_fock_matrix_ph(qs_env,ks_matrix,error)
       ! QM/MM
       IF (qs_env%qmmm) THEN
          DO ispin = 1, nspins
             ! If QM/MM sumup the 1el Hamiltonian
             CALL cp_dbcsr_add(ks_matrix(ispin)%matrix,qs_env%ks_qmmm_env%matrix_h(1)%matrix,&
                  1.0_dp,1.0_dp,error=error)
             ! Compute QM/MM Energy
             CALL cp_dbcsr_trace(qs_env%ks_qmmm_env%matrix_h(1)%matrix,&
                  matrix_p(ispin)%matrix,trace=qmmm_el,error=error)
             energy%qmmm_el = energy%qmmm_el + qmmm_el
          END DO
       END IF

!      WRITE ( *, * ) ' before TOTAL', energy%total
       ! Collect all the energy terms
       energy%mulliken = 0.0_dp
       energy%exc      = 0.0_dp
       energy%total    = energy%total+energy%core +&
                         energy%core_overlap +&
                         0.5_dp*energy%hartree + &
                         energy%qmmm_el +&
                         energy%mulliken
!      WRITE ( *, * ) ' AFTER TOTAL', energy%total


       output_unit=cp_print_key_unit_nr(logger,scf_section,"PRINT%DETAILED_ENERGY",&
            extension=".scfLog",error=error)

       IF (output_unit>0) THEN
          WRITE (UNIT=output_unit,FMT="(/,(T3,A,T60,F20.10))")&
               "Core Hamiltonian energy:                       ",energy%core,&
               "Two-electron integral energy:                  ",energy%hartree
          IF (scp_nddo) THEN
             WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
                 "SCP electrostatic energy:             ",energy_scp%e_scp_hartree
             WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
                 "SCP kernel energy:                    ",energy_scp%e_scp_kernel
             IF ( ABS(energy_scp%e_scp_disp) > 0._dp ) &
                 WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
                 "SCP dispersion energy:                ",energy_scp%e_scp_disp
          END IF
          IF (qs_env%qmmm) THEN
             WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
                  "QM/MM Electrostatic energy:                    ",energy%qmmm_el
          END IF
       END IF

       CALL cp_print_key_finished_output(output_unit,logger,scf_section,&
            "PRINT%DETAILED_ENERGY", error=error)

       ! Here we compute dE/dC if needed. Assumes dE/dC is H_{ks}C (plus occupation numbers)
       IF (qs_env%requires_mo_derivs .AND. .NOT. just_energy) THEN
          CALL get_qs_env(qs_env,mo_derivs=mo_derivs,mos=mo_array,error=error)
          DO ispin=1,SIZE(mo_derivs)
             CALL get_mo_set(mo_set=mo_array(ispin)%mo_set, &
                  mo_coeff_b=mo_coeff, occupation_numbers=occupation_numbers )
             IF(.NOT.mo_array(ispin)%mo_set%use_mo_coeff_b) THEN
                CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
             ENDIF
             CALL cp_dbcsr_get_info(mo_coeff,nfullcols_total=ncol_global)
             CALL cp_dbcsr_multiply('n','n',1.0_dp,ks_matrix(ispin)%matrix,mo_coeff,&
                  0.0_dp,mo_derivs(ispin)%matrix, error=error)
          ENDDO
       ENDIF

    END SELECT

    CALL timestop(handle)

  END SUBROUTINE build_se_fock_matrix

END MODULE se_fock_matrix_mpi

