!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
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
MODULE se_fock_matrix
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cp_control_types,                ONLY: dft_control_type,&
                                             semi_empirical_control_type
  USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set,&
                                             dbcsr_deallocate_matrix,&
                                             dbcsr_deallocate_matrix_set,&
                                             dbcsr_from_sm,&
                                             dbcsr_obj_p_type,&
                                             sm_from_dbcsr
  USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE dbcsr_operations,                ONLY: dbcsr_add,&
                                             dbcsr_copy,&
                                             dbcsr_trace
  USE dbcsr_types,                     ONLY: dbcsr_obj
  USE distribution_2d_types,           ONLY: distribution_2d_type
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
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE scp_energy_types,                ONLY: scp_energy_type
  USE scp_environment_types,           ONLY: get_scp_env,&
                                             scp_environment_type
  USE scp_nddo_utils,                  ONLY: evaluate_scp_force
  USE se_fock_matrix_coulomb,          ONLY: build_fock_matrix_coul_lr_r3,&
                                             build_fock_matrix_coulomb,&
                                             build_fock_matrix_coulomb_lr
  USE se_fock_matrix_dbg,              ONLY: dbg_energy_coulomb_lr
  USE se_fock_matrix_exchange,         ONLY: build_fock_matrix_exchange
  USE semi_empirical_store_int_types,  ONLY: semi_empirical_si_finalize,&
                                             semi_empirical_si_initialize,&
                                             semi_empirical_si_type
  USE sparse_matrix_types,             ONLY: deallocate_matrix,&
                                             real_matrix_p_type,&
                                             real_matrix_type,&
                                             set_matrix
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'se_fock_matrix'
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
  SUBROUTINE build_se_fock_matrix (ks_env,qs_env,ks_matrix_sm,rho,energy,&
       calculate_forces,just_energy,error)

    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix_sm
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(qs_energy_type), POINTER            :: energy
    LOGICAL, INTENT(in)                      :: calculate_forces, just_energy
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_se_fock_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, ncol_global, &
                                                nspins, output_unit
    LOGICAL                                  :: failure, scp_nddo
    REAL(dp), POINTER                        :: ascp_diag( : )
    REAL(KIND=dp)                            :: ecoul, ecoul_scp, qmmm_el
    REAL(KIND=dp), DIMENSION(:), POINTER     :: occupation_numbers
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mo_derivs
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dbcsr_obj), POINTER                 :: ks_scp, matrix_tmp, pscp
    TYPE(dbcsr_obj_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix, matrix_h, matrix_p
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_h_sm, matrix_p_sm
    TYPE(real_matrix_type), POINTER          :: ks_scp_sm, pscp_sm
    TYPE(scp_energy_type), POINTER           :: energy_scp
    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(section_vals_type), POINTER         :: scf_section
    TYPE(semi_empirical_control_type), &
      POINTER                                :: se_control
    TYPE(semi_empirical_si_type), POINTER    :: store_int_env

    failure=.FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(matrix_h, dft_control, logger, scf_section, store_int_env, se_control, matrix_tmp)
    NULLIFY(scp_env, pscp, ks_scp_sm, ascp_diag, energy_scp,atomic_kind_set, pscp, ks_scp)
    logger => cp_error_get_logger(error)
    CPPrecondition(ASSOCIATED(ks_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ks_env%ref_count>0,cp_failure_level,routineP,error,failure)

    CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, matrix_h=matrix_h_sm,&
         para_env=para_env, se_store_int_env=store_int_env, scp_env=scp_env,&
         atomic_kind_set=atomic_kind_set, distribution_2d=distribution_2d,&
         error=error)

    NULLIFY(ks_matrix,matrix_h)!sm->dbcsr
    CALL dbcsr_allocate_matrix_set( ks_matrix, SIZE(ks_matrix_sm), error )!sm->dbcsr
    CALL dbcsr_allocate_matrix_set( matrix_h, SIZE(matrix_h_sm), error )!sm->dbcsr
    DO ispin=1,SIZE(ks_matrix)!sm->dbcsr
       ALLOCATE(ks_matrix(ispin)%matrix)!sm->dbcsr
       CALL dbcsr_from_sm(ks_matrix(ispin)%matrix, ks_matrix_sm(ispin)%matrix, error)!sm->dbcsr
    ENDDO!sm->dbcsr
    DO ispin=1,SIZE(matrix_h)!sm->dbcsr
       ALLOCATE(matrix_h(ispin)%matrix)!sm->dbcsr
       CALL dbcsr_from_sm(matrix_h(ispin)%matrix, matrix_h_sm(ispin)%matrix, error)!sm->dbcsr
    ENDDO!sm->dbcsr

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
       matrix_p_sm => rho%rho_ao

       NULLIFY(matrix_p)
       CALL dbcsr_allocate_matrix_set( matrix_p, SIZE(matrix_p_sm), error )!sm->dbcsr
       DO ispin=1,SIZE(matrix_p)!sm->dbcsr
          ALLOCATE(matrix_p(ispin)%matrix)!sm->dbcsr
          CALL dbcsr_from_sm(matrix_p(ispin)%matrix, matrix_p_sm(ispin)%matrix, error)!sm->dbcsr
       ENDDO!sm->dbcsr

       energy%qmmm_el = 0.0_dp

       CALL dbcsr_deallocate_matrix_set( ks_matrix, error )!sm->dbcsr !needed a soft copy!
       CALL dbcsr_allocate_matrix_set( ks_matrix, nspins, error )!sm->dbcsr
       DO ispin=1,nspins
          ALLOCATE(ks_matrix(ispin)%matrix)!sm->dbcsr
          ! Copy the core matrix into the fock matrix
          CALL dbcsr_copy(ks_matrix(ispin)%matrix,matrix_h(1)%matrix,error=error)
       END DO

       ! Zero the ks_scp matrix.  Will get info from coulomb routines
       IF (scp_nddo) THEN
          CALL get_scp_env(scp_env, ks_scp=ks_scp_sm, error=error)
          CALL set_matrix(ks_scp_sm,0.0_dp)
       END IF



       ! Compute Exchange and Coulomb terms
       CALL semi_empirical_si_initialize(store_int_env, ks_env%s_mstruct_changed, error)
       CALL build_fock_matrix_exchange(qs_env,ks_matrix,matrix_p,calculate_forces,&
            store_int_env,error)
       CALL build_fock_matrix_coulomb(qs_env,ks_matrix,matrix_p,energy,calculate_forces,&
            store_int_env,error)

       ! Debug statements for Long-Range
       IF (debug_energy_coulomb_lr.AND.se_control%do_ewald) THEN
          CALL dbg_energy_coulomb_lr(energy, ks_matrix, nspins, qs_env, matrix_p,&
                                     calculate_forces, para_env, store_int_env, error)
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



       ! Compute the Hartree energy
       ! NOTE: If we are performing SCP-NDDO, ks_matrix contains coulomb piece from SCP.
       DO ispin=1,nspins
          CALL dbcsr_trace(ks_matrix(ispin)%matrix,matrix_p(ispin)%matrix,trace=ecoul,error=error)
          energy%hartree = energy%hartree + ecoul
       END DO
       ! Compute the SCP Hartree energy and polarization energy and update ks_scp
       IF (scp_nddo) THEN
          CALL get_scp_env ( scp_env, pscp=pscp_sm, ks_scp=ks_scp_sm, energy=energy_scp, &
                             ascp_diag= ascp_diag, error = error )
          ALLOCATE(pscp, ks_scp)!sm->dbcsr
          CALL dbcsr_from_sm(pscp, pscp_sm, error)!sm->dbcsr
          CALL dbcsr_from_sm(ks_scp, ks_scp_sm, error)!sm->dbcsr
          ! Get energy
          CALL dbcsr_trace(ks_scp,pscp,trace=ecoul_scp,local_sum=.TRUE.,error=error)
          energy_scp%e_scp_hartree = ecoul_scp*0.5_dp
          ! evalulate kernel and forces on the SCP coeffs
          CALL evaluate_scp_force ( atomic_kind_set, pscp, ks_scp, ascp_diag, para_env,&
                                    just_energy, energy_scp%e_scp_kernel, error )
          CALL deallocate_matrix(scp_env%ks_scp, error=error)!sm->dbcsr
          CALL sm_from_dbcsr(scp_env%ks_scp, ks_scp, distribution_2d, error=error )!sm->dbcsr
          CALL get_scp_env ( scp_env, ks_scp=ks_scp_sm, error=error )!sm->dbcsr
          CALL dbcsr_deallocate_matrix( pscp, error=error )!sm->dbcsr
          CALL dbcsr_deallocate_matrix( ks_scp, error=error )!sm->dbcsr

          ! NOTE: If needed, at this point, forces on P_SCP are up-to-date and stored 
          !       in KS_SCP as a*ks_scp
       END IF

       ! QM/MM
       IF (qs_env%qmmm) THEN
          DO ispin = 1, nspins
             ! If QM/MM sumup the 1el Hamiltonian
             ALLOCATE(matrix_tmp)!sm->dbcsr
             CALL dbcsr_from_sm(matrix_tmp, qs_env%ks_qmmm_env%matrix_h(1)%matrix, error)!sm->dbcsr
             CALL dbcsr_add(ks_matrix(ispin)%matrix,matrix_tmp,error=error)
             ! Compute QM/MM Energy
             CALL dbcsr_trace(matrix_tmp,matrix_p(ispin)%matrix,trace=qmmm_el,error=error)
             CALL dbcsr_deallocate_matrix(matrix_tmp,error=error)!sm->dbcsr
             energy%qmmm_el = energy%qmmm_el + qmmm_el
          END DO
       END IF
10 CONTINUE
       ! Collect all the energy terms
       energy%mulliken = 0.0_dp
       energy%exc      = 0.0_dp
       energy%total    = energy%core +&
                         energy%core_overlap +&
                         0.5_dp*energy%hartree + &
                         energy%qmmm_el +&
                         energy%mulliken

       ! Modify energy for SCP
       IF (scp_nddo) THEN 
          energy%total = energy%total + energy_scp%e_scp_hartree + energy_scp%e_scp_kernel + &
                                        energy_scp%e_scp_disp
       END IF

       output_unit=cp_print_key_unit_nr(logger,scf_section,"PRINT%DETAILED_ENERGY",&
            extension=".scfLog",error=error)

       IF (output_unit>0) THEN
          WRITE (UNIT=output_unit,FMT="(/,(T3,A,T60,F20.10))")&
               "Core Hamiltonian energy:                       ",energy%core,&
               "Two-electron integral energy:                  ",energy%hartree
          IF (qs_env%qmmm) THEN
             WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
                  "QM/MM Electrostatic energy:                    ",energy%qmmm_el
          END IF
          IF (scp_nddo) THEN
             WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
                  "SCP electrostatic energy:             ",energy_scp%e_scp_hartree
             WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
                  "SCP kernel energy:                    ",energy_scp%e_scp_kernel
             IF ( ABS(energy_scp%e_scp_disp) > 0._dp ) &
               WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
                  "SCP dispersion energy:                ",energy_scp%e_scp_disp
          END IF
       END IF

       CALL cp_print_key_finished_output(output_unit,logger,scf_section,&
            "PRINT%DETAILED_ENERGY", error=error)

       DO ispin=1,SIZE(ks_matrix)!sm->dbcsr
          CALL sm_from_dbcsr ( ks_matrix_sm(ispin)%matrix, ks_matrix(ispin)%matrix, &!sm->dbcsr
               distribution_2d, error=error )!sm->dbcsr
       ENDDO!sm->dbcsr

       ! Here we compute dE/dC if needed. Assumes dE/dC is H_{ks}C (plus occupation numbers)
       IF (qs_env%requires_mo_derivs .AND. .NOT. just_energy) THEN
          CALL get_qs_env(qs_env,mo_derivs=mo_derivs,mos=mo_array,error=error)
          DO ispin=1,SIZE(mo_derivs)
             CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,mo_coeff=mo_coeff, &
                             occupation_numbers=occupation_numbers )
             CALL cp_fm_get_info(mo_coeff,ncol_global=ncol_global,error=error)
             CALL cp_sm_fm_multiply(ks_matrix_sm(ispin)%matrix,mo_coeff,&
                                    mo_derivs(ispin)%matrix, ncol_global,error=error)
          ENDDO
       ENDIF

    END SELECT

    CALL dbcsr_deallocate_matrix_set ( matrix_h, error=error )
    CALL dbcsr_deallocate_matrix_set ( ks_matrix, error=error )
    CALL dbcsr_deallocate_matrix_set ( matrix_p, error=error )
    !CALL dbcsr_deallocate_matrix ( matrix_tmp, error=error )
    !CALL dbcsr_deallocate_matrix ( ks_scp, error=error )
    !CALL dbcsr_deallocate_matrix ( pscp, error=error )

    CALL timestop(handle)

  END SUBROUTINE build_se_fock_matrix
    
END MODULE se_fock_matrix

