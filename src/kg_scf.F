!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2004 CP2K developers group
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****** cp2k/kg_scf [1.0] *
!!
!!   NAME
!!     kg_scf
!!
!!   FUNCTION
!!     Routines for the Kim-Gordon type SCF run.
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************
 
MODULE kg_scf

!  USE atomic_kind_types,               ONLY: atomic_kind_type,&
!                                             get_atomic_kind,&
!                                             get_atomic_kind_set
!  USE basis_set_types,                 ONLY: get_gto_basis_set,&
!                                             gto_basis_set_type
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
!  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm,&
!                                             cp_fm_symm,&
!                                             cp_fm_upper_to_full
!  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
!                                             cp_fm_cholesky_reduce,&
!                                             cp_fm_cholesky_restore
!  USE cp_fm_diag,                      ONLY: cp_fm_block_jacobi_classic,&
!                                             cp_fm_power,&
!                                             cp_fm_syevd,&
!                                             cp_fm_syevx
!  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
!                                             fm_pool_get_el_struct
!  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
!                                             cp_fm_struct_get,&
!                                             cp_fm_struct_release,&
!                                             cp_fm_struct_type
!  USE cp_fm_types,                     ONLY: cp_fm_create,&
!                                             cp_fm_get_element,&
!                                             cp_fm_init_random,&
!                                             cp_fm_p_type,&
!                                             cp_fm_release,&
!                                             cp_fm_set_element,&
!                                             cp_fm_to_fm,&
!                                             cp_fm_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
!  USE cp_para_types,                   ONLY: cp_para_env_type
!  USE cp_sm_fm_interactions,           ONLY: copy_sm_to_fm,&
!                                             cp_sm_fm_multiply
  USE cp_sm_pool_types,                ONLY: cp_sm_pool_p_type,&
                                             sm_pools_create_matrix_vect
!                                             sm_pools_flush_cache,&
!                                             sm_pools_give_back_matrix_vect
!  USE external_potential_types,        ONLY: all_potential_type,&
!                                             get_potential,&
!                                             gth_potential_type
  USE global_types,                    ONLY: global_environment_type
  USE kg_environment_types,            ONLY: kg_environment_type,&
                                              get_kg_env,&
                                              set_kg_env
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             int_size
  USE machine,                         ONLY: m_flush,&
                                             m_walltime
!  USE message_passing,                 ONLY: mp_bcast,&
!                                             mp_max,&
!                                             mp_sum
!  USE pao_types,                       ONLY: pao_env_type
!  USE particle_types,                  ONLY: get_particle_set,&
!                                             particle_type
!  USE preconditioner,                  ONLY: destroy_preconditioner,&
!                                             init_preconditioner,&
!                                             make_preconditioner
!  USE pw_env_types,                    ONLY: pw_env_type
!  USE qs_charges_types,                ONLY: qs_charges_type
  USE qs_core_energies,                ONLY: calculate_ecore_overlap,&
                                             calculate_ecore_self
!  USE qs_diis,                         ONLY: qs_diis_b_clear,&
!                                             qs_diis_b_create,&
!                                             qs_diis_b_step
!  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_ks_methods,                   ONLY: qs_ks_create,&
!                                             qs_ks_did_change,&
                                             qs_ks_update_qs_env
  USE qs_ks_types,                     ONLY: qs_ks_env_type,&
                                             qs_ks_release
  USE qs_matrix_pools,                 ONLY: mpools_get
!  USE qs_mo_methods,                   ONLY: calculate_density_matrix,&
!                                             calculate_orthonormality,&
!                                             make_basis
!  USE qs_mo_types,                     ONLY: correct_mo_eigenvalues,&
!                                             get_mo_set,&
!                                             init_mo_set,&
!                                             mo_set_p_type,&
!                                             mo_set_type,&
!                                             read_mo_set,&
!                                             set_mo_occupation,&
!                                             write_mo_set
  USE qs_overlap,                      ONLY: write_fm_with_basis_info,&
                                             write_sparse_matrix
!  USE qs_rho_atom_methods,             ONLY: calculate_rho_atom_coeff
!  USE qs_rho_methods,                  ONLY: qs_rho_update_rho
  USE qs_rho_types,                    ONLY: qs_rho_type
!  USE qs_scf_post,                     ONLY: scf_post_calculation
   USE kg_scf_types,                    ONLY: kg_scf_env_type,&
                                              kg_scf_env_release,&
                                              kg_scf_env_create
!  USE qs_wf_history_methods,           ONLY: wfi_extrapolate,&
!                                             wfi_get_method_label,&
!                                             wfi_update
!  USE qs_wf_history_types,             ONLY: wfi_use_guess_method_nr
  USE scf_control_types,               ONLY: scf_control_type
!  USE simulation_cell,                 ONLY: cell_type
  USE sparse_matrix_types,             ONLY: real_matrix_p_type
!  USE sparse_matrix_types,             ONLY: cp_sm_scale_and_add,&
!                                             deallocate_matrix_set,&
!                                             first_block_node,&
!                                             get_block_node,&
!                                             get_matrix_diagonal,&
!                                             get_matrix_info,&
!                                             next_block_node,&
!                                             real_block_node_type,&
!                                             real_matrix_p_type,&
!                                             real_matrix_type,&
!                                             set_matrix_diagonal,&
!                                             transfer_matrix
!  USE termination,                     ONLY: stop_memory,&
!                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: moduleN = "kg_scf"

  LOGICAL, PARAMETER, PRIVATE :: debug_this_module=.FALSE.

! *** Public subroutines ***

  PUBLIC :: kg_qs_scf

! *****************************************************************************

CONTAINS

! *****************************************************************************

!!****f* kg_scf/kg_qs_scf [1.0] *
!!
!!   NAME
!!     kg_qs_scf
!!
!!   AUTHOR
!!     fawzi/jgh
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************
  SUBROUTINE kg_qs_scf(kg_env,globenv,error)
    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_qs_scf', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
     TYPE(kg_environment_type), POINTER       :: my_kg_env
     TYPE(kg_scf_env_type), POINTER           :: scf_env

     NULLIFY(scf_env,my_kg_env)
     my_kg_env => kg_env
     failure=.FALSE.
     CPPrecondition(associated(kg_env),cp_failure_level,routineP,error,failure)
     NULLIFY(scf_env)
     IF (.NOT. failure) THEN
        CALL get_kg_env(kg_env,scf_env=scf_env,error=error)
        IF (.not.ASSOCIATED(scf_env)) THEN
           CALL kg_scf_env_create(scf_env, error=error)
           CALL set_kg_env(kg_env,scf_env=scf_env,error=error)
           CALL kg_scf_env_release(scf_env)
           CALL get_kg_env(kg_env,scf_env=scf_env,error=error)
        END IF
        CALL scf_env_do_scf(scf_env,my_kg_env,globenv,error=error)
     END IF
  END SUBROUTINE kg_qs_scf
!***************************************************************************

!!****f* kg_scf/scf_env_do_scf [1.0] *
!!
!!   NAME
!!     scf_env_do_scf
!!
!!   AUTHOR
!!     jgh (adapted from QS)
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************
  SUBROUTINE scf_env_do_scf(scf_env,kg_env,globenv,error)

    TYPE(kg_scf_env_type), POINTER           :: scf_env
    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "scf_env_do_scf", &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: t1, t2
    INTEGER                                  :: handle, handle2
    INTEGER                                  :: output_unit, ispin
    LOGICAL                                  :: failure
    LOGICAL                                  :: energy_only
    LOGICAL                                  :: ionode
!    LOGICAL :: diis_step, do_level_shift,&
!      should_stop, use_cholesky, use_jacobi
!    TYPE(atomic_kind_type), DIMENSION(:), &
!      POINTER                                :: atomic_kind_set
!    TYPE(cell_type), POINTER                 :: cell
!    TYPE(cp_fm_p_type), DIMENSION(:), &
!      POINTER                                :: scf_work1
!    TYPE(cp_fm_type), POINTER                :: ortho, scf_work2
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(qs_ks_env_type), POINTER            :: ks_env
!    TYPE(particle_type), DIMENSION(:), &
!      POINTER                                :: particle_set
!    TYPE(pw_env_type), POINTER               :: pw_env
!    TYPE(qs_charges_type), POINTER           :: qs_charges
!    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: p,ks

!   ---------------------------------------------------------------------------

    NULLIFY(dft_control,scf_control,ks_env,qs_env,p,ks)
    NULLIFY(rho)

    logger => cp_error_get_logger(error)
    t1 = m_walltime()
    failure=.FALSE.
 
    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(kg_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(kg_env%ref_count>0,cp_failure_level,routineP,error,failure)
     
     CALL get_kg_env(kg_env=kg_env,&
          sub_qs_env=qs_env,&
          dft_control=dft_control,&
          scf_control=scf_control)

    ! Quick return, if no SCF iteration is requested
    IF (scf_control%max_scf < 1) RETURN

    CALL write_checkpoint_information("entering "//routineN,globenv)

    CALL timeset(routineN,"I"," ",handle)

    ionode = logger%para_env%source==logger%para_env%mepos
    IF (ionode) output_unit = cp_logger_get_default_unit_nr(logger)
    energy_only = .FALSE.

    IF (ionode.AND.logger%print_keys%scf) THEN
       WRITE (UNIT=output_unit,FMT="(/,/,T2,A)")&
            "SCF WAVEFUNCTION OPTIMIZATION"
    END IF

    CALL init_scf_run(scf_env=scf_env,kg_env=kg_env,&
         globenv=globenv,error=error)

    CALL get_qs_env(qs_env=qs_env,ks_env=ks_env)
!
!    scf_work1 => scf_env%scf_work1
!    scf_work2 => scf_env%scf_work2
!    ortho => scf_env%ortho
!
    CALL get_qs_env(qs_env,matrix_ks=ks,rho=rho)
    p => rho%rho_ao

    scf_env%iter_count = 0
!    diis_step = .FALSE.
!    use_jacobi = .FALSE.

    IF (ionode.AND.logger%print_keys%scf) THEN
       WRITE (UNIT=output_unit,&
            FMT="(/,T3,A,T9,A,T36,A,T49,A,T68,A,/,T3,A)")&
            "Step","Update method","Time","Convergence","Total energy",&
            REPEAT("-",77)
    END IF
    
    IF (.NOT. logger%print_keys%timecumul) t1 = m_walltime()

!   *** SCF loop ***

    scf_loop: DO
       CALL timeset("scf_iter","I"," ",handle2)

       IF (ionode) CALL m_flush(output_unit)

       scf_env%iter_count = scf_env%iter_count + 1

! ** here qs_env%rho%rho_r and qs_env%rho%rho_g should be up to date

       CALL qs_ks_update_qs_env(ks_env,qs_env=qs_env,&
            error=error,&
            calculate_forces=.FALSE.,just_energy=energy_only)

       DO ispin=1,dft_control%nspins

          IF (logger%print_keys%density_matrix.AND.&
               logger%print_keys%each_scf_step) THEN
             CALL write_sparse_matrix(p(ispin)%matrix,4,6,qs_env,globenv)
          END IF

          IF (logger%print_keys%kohn_sham_matrix.AND.logger%print_keys%each_scf_step) THEN
             CALL write_sparse_matrix(ks(ispin)%matrix,4,6,qs_env,globenv)
          END IF

       ENDDO

!      IF (logger%print_keys%mo_orthonormality) THEN
!          CALL calculate_orthonormality(orthonormality,c,s(1)%matrix)
!         IF (ionode) THEN
!            WRITE(output_unit,'(T8,A,T60,E20.4)')  &
!                   " Maximum deviation from MO S-orthonormality",orthonormality
!         ENDIF
!      ENDIF

       scf_env%p_mix = scf_control%p_mix
       scf_env%iter_param = 0.0_dp

       SELECT CASE (scf_env%method)
       CASE DEFAULT
!
!          DO ispin=1,dft_control%nspins 
!             CALL copy_sm_to_fm(ks(ispin)%matrix,scf_work1(ispin)%matrix)
!          ENDDO
!
!          IF (scf_env%iter_count > 1) THEN
!             CALL qs_diis_b_step(scf_env%scf_diis_buffer,c,scf_work1,&
!                  scf_work2,scf_env%iter_delta,diis_error,diis_step,&
!                  !!!!,error=error)
!          END IF
!
!          do_level_shift = ((scf_control%level_shift /= 0.0_dp).AND.&
!               ((scf_control%density_guess == "CORE").OR.(scf_env%iter_count > 1)))
!
!          IF (diis_step) THEN
!             scf_env%p_mix = 1.0_dp
!             scf_env%iter_param = diis_error
!             IF (use_jacobi) THEN
!                scf_env%iter_method = "DIIS/Jacobi"
!             ELSE
!                scf_env%iter_method = "DIIS/Diag."
!             END IF
!          ELSE
!             IF (use_jacobi) THEN
!                scf_env%iter_method = "Mixing/Jacobi"
!             ELSE
!                scf_env%iter_method = "Mixing/Diag."
!             END IF
!          END IF
!
!          IF ((scf_env%iter_count > 1).AND.(scf_env%iter_delta < scf_control%eps_jacobi)) THEN
!             use_jacobi = .TRUE.
!          ELSE
!             use_cholesky = scf_control%use_cholesky
!             use_jacobi = .FALSE.
!          END IF
!
!          scf_env%iter_delta=0.0_dp
!
!          DO ispin=1,dft_control%nspins
!
!             CALL eigensolver(scf_work1(ispin)%matrix,c(ispin)%mo_set,ortho,scf_work2,&
!                  do_level_shift,scf_control%level_shift,&
!                  use_cholesky=use_cholesky,&
!                  work_syevx=scf_control%work_syevx,&
!                  use_jacobi=use_jacobi,&
!                  jacobi_threshold=scf_control%jacobi_threshold,&
!                  smear=scf_control%smear, &
!                  error=error)
!
!             IF (logger%print_keys%each_scf_step) THEN
!                CALL write_mo_set(c(ispin)%mo_set,atomic_kind_set,particle_set,4,6,globenv)
!             END IF
!
!             CALL calculate_density_matrix(c(ispin)%mo_set,&
!                  scf_env%p_mix_new(ispin)%matrix,error=error)
!
!          ENDDO

       END SELECT

!       IF (scf_env%mixing_method.EQ.ao_mixing_nr) THEN
!          CALL scf_env_density_mixing(scf_env,qs_env=qs_env,error=error)
!       ENDIF
!
!       t2 = m_walltime()
!
!       IF (ionode.AND.logger%print_keys%scf.and.scf_env%print_iter_line) THEN
!          WRITE (UNIT=output_unit,&
!               FMT="(T2,I5,2X,A,T22,E10.2,T32,F8.2,T40,2F20.10)")&
!               scf_env%iter_count,TRIM(scf_env%iter_method),&
!               scf_env%iter_param,t2 - t1,scf_env%iter_delta,energy%total
!       END IF
!
!! ** convergece check
!       IF (qs_env%para_env%source==qs_env%para_env%mepos) &
!            INQUIRE(file="EXIT",exist=should_stop)
!       CALL mp_bcast(should_stop,qs_env%para_env%source,qs_env%para_env%group)
!       IF (scf_env%iter_delta < scf_control%eps_scf) THEN
!          IF (ionode.AND.logger%print_keys%scf) THEN
!             WRITE(UNIT=output_unit,FMT="(/,T3,A,I5,A/)")&
!                  "*** SCF run converged in ",scf_env%iter_count," steps ***"
!          END IF
!          CALL timestop(0.0_dp,handle2)
!          EXIT scf_loop
!       ELSE IF (should_stop.OR.&
!            scf_env%iter_count == scf_control%max_scf) THEN
!          IF (ionode.AND.logger%print_keys%scf) THEN
!             WRITE(UNIT=output_unit,FMT="(/,T3,A,/)")&
!                  "*** SCF run NOT converged ***"
!          END IF
!          CALL timestop(0.0_dp,handle2)
!          EXIT scf_loop
!       END IF
!
!!   *** Write restart file ***
!       IF (scf_control%write_restart_each>0.AND.&
!            MODULO(scf_env%iter_count,scf_control%write_restart_each)==0) THEN
!         ! CALL write_mo_set(c,atomic_kind_set,particle_set,globenv,&
!         !      id_nr=qs_env%id_nr)
!          CALL write_mo_set(c,atomic_kind_set,particle_set,globenv,&
!               id_nr=1)
!       END IF
!
!       IF (.NOT. logger%print_keys%timecumul) t1 = m_walltime()
!
!!   *** mixing methods have the new density matrix in p_mix_new
!       IF (scf_env%mixing_method.EQ.ao_mixing_nr) THEN
!           DO ispin=1,dft_control%nspins
!              CALL transfer_matrix(scf_env%p_mix_new(ispin)%matrix,&
!                                       rho%rho_ao(ispin)%matrix)
!           END DO
!       ENDIF
!
!! ** update qs_env%rho
!       CALL qs_rho_update_rho(rho, qs_env=qs_env, error=error)
!       CALL qs_ks_did_change(ks_env,rho_changed=.TRUE.,error=error)      

       CALL timestop(0.0_dp,handle2)

    END DO scf_loop

!    IF (ionode.AND.logger%print_keys%scf) THEN
!       WRITE (UNIT=output_unit,FMT="(/,(T3,A,T40,2F20.10))")&
!            "Total electronic density (r-space): ",&
!            SUM(rho%tot_rho_r),&
!            SUM(rho%tot_rho_r)+ REAL(scf_env%nelectron,dp),&
!            "Total core charge density (r-space):",&
!            qs_charges%total_rho_core_rspace,&
!            qs_charges%total_rho_core_rspace - REAL(scf_env%nelectron+dft_control%charge,dp)
!       WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
!            "Total charge density (r-space):     ",&
!            SUM(rho%tot_rho_r)+&
!            qs_charges%total_rho_core_rspace,&
!            "Total charge density (g-space):     ",qs_charges%total_rho_gspace
!       WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
!            "Overlap energy of the core charge distribution:",energy%core_overlap,&
!            "Self energy of the core charge distribution:   ",energy%core_self,&
!            "Core Hamiltonian energy:                       ",energy%core,&
!            "Hartree energy:                                ",energy%hartree,&
!            "Exchange-correlation energy:                   ",energy%exc,&
!            "Total energy:                                  ",energy%total
!       CALL m_flush(output_unit)
!    END IF
!
!!   *** Write restart file ***
!    IF (scf_control%write_restart_each>=0) THEN
!       CALL write_mo_set(c,atomic_kind_set,particle_set,globenv,&
!            id_nr=1)
!    END IF
!
!!   *** add the converged wavefunction to the wavefunction history
!!   should possibly be moved down after de-mixing
!    IF (ASSOCIATED(qs_env%wf_history)) THEN
!        CALL wfi_update(qs_env%wf_history,qs_env=qs_env,dt=1.0_dp,&
!               error=error)
!    END IF
!
!!   *** mixing methods need to undo mixing of the density matrix (restore original density) ***
!    IF (scf_env%mixing_method.EQ.ao_mixing_nr) THEN
!       CALL scf_env_density_mixing(scf_env,qs_env=qs_env,invert=.TRUE.,&
!                                       error=error)
!       DO ispin=1,dft_control%nspins
!          CALL transfer_matrix(scf_env%p_mix_new(ispin)%matrix,&
!                                     rho%rho_ao(ispin)%matrix)
!       END DO
!       DO ispin=1,dft_control%nspins
!          CALL transfer_matrix(scf_env%p_mix_new(ispin)%matrix,&
!                               rho%rho_ao(ispin)%matrix)
!       END DO
!    ENDIF
!
!!   *** update rspace rho since the mo changed
!!   *** this might not always be needed (i.e. no post calculation / no forces )
!    CALL qs_rho_update_rho(rho, qs_env=qs_env, error=error)
!
!    CALL qs_ks_did_change(ks_env,rho_changed=.TRUE.,error=error)
!
!! *** print out info
!    IF (.NOT.should_stop) THEN
!       CALL scf_post_calculation(scf_env=scf_env, qs_env=qs_env, &
!            globenv=globenv, error=error)
!    END IF
!
!! *** cleanup
!    CALL scf_env_cleanup(scf_env,qs_env=qs_env,error=error)
!
    CALL timestop(0.0_dp,handle)

    CALL write_checkpoint_information("leaving "//routineN,globenv)

  END SUBROUTINE scf_env_do_scf
!**************************************************************************

!!!****f* qs_scf/scf_env_density_mixing [1.0] *
!!!
!!!   NAME
!!!     scf_env_density_mixing
!!!
!!!   FUNCTION
!!!     perform (if requested) a density mixing
!!!
!!!   NOTES
!!!     -
!!!
!!!   ARGUMENTS
!!!     - error: variable to control error logging, stopping,... 
!!!       see module cp_error_handling 
!!!
!!!   AUTHOR
!!!     fawzi
!!!
!!!   MODIFICATION HISTORY
!!!     02.2003 created [fawzi]
!!!
!!!*** **********************************************************************
!  SUBROUTINE scf_env_density_mixing(scf_env,qs_env,invert,error)
!    TYPE(qs_scf_env_type), POINTER           :: scf_env
!    TYPE(qs_environment_type), POINTER       :: qs_env
!    LOGICAL, INTENT(in), OPTIONAL            :: invert
!    TYPE(cp_error_type), INTENT(inout), &
!      OPTIONAL                               :: error
!
!    CHARACTER(len=*), PARAMETER :: routineN = 'scf_env_density_mixing', &
!      routineP = moduleN//':'//routineN
!
!    INTEGER                                  :: ispin
!    LOGICAL                                  :: failure, my_invert
!    REAL(KIND=dp)                            :: tmp
!
!    failure=.FALSE.
!    my_invert=.FALSE.
!    IF (PRESENT(invert)) my_invert=invert
!
!    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
!    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
!    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
!    CPPrecondition(qs_env%ref_count>0,cp_failure_level,routineP,error,failure)
!    IF (.NOT. failure) THEN
!       SELECT CASE (scf_env%mixing_method)
!       CASE(no_mixing_nr)
!          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
!               routineP,"don't call mixing in no_mixing_nr:",error,failure)
!       CASE(ao_mixing_nr)
!          CPPrecondition(ASSOCIATED(scf_env%p_mix_new),cp_failure_level,routineP,error,failure)
!!FM        scf_env%iter_delta=0.0_dp
!          DO ispin=1,SIZE(scf_env%p_mix_new)
!             IF (my_invert) THEN
!                CPPrecondition(scf_env%p_mix/=0.0_dp,cp_failure_level,routineP,error,failure)
!                IF (scf_env%p_mix/=1.0_dp) THEN
!                   CALL cp_sm_scale_and_add(matrix_a=scf_env%p_mix_new(ispin)%matrix,&
!                        alpha=1.0_dp/scf_env%p_mix,&
!                        matrix_b=qs_env%rho%rho_ao(ispin)%matrix,&
!                        beta=(scf_env%p_mix-1.0_dp)/scf_env%p_mix,&
!                        error=error)
!                END IF
!             ELSE
!                CALL cp_sm_mix(m1=scf_env%p_mix_new(ispin)%matrix,&
!                     m2=qs_env%rho%rho_ao(ispin)%matrix,&
!                     p_mix=scf_env%p_mix,&
!                     delta=tmp,&
!                     para_env=qs_env%para_env,&
!                     error=error)
!                scf_env%iter_delta=MAX(scf_env%iter_delta,tmp)
!             END IF
!          END DO
!       CASE default
!          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
!               routineP,"unknown density mixing method:"//&
!               cp_to_string(scf_env%mixing_method),error,failure)
!       END SELECT
!    END IF
!  END SUBROUTINE scf_env_density_mixing
!!***************************************************************************

  SUBROUTINE init_scf_run(scf_env,kg_env,globenv,error)

! ***************************************************************************

    TYPE(kg_scf_env_type), POINTER           :: scf_env
    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(global_environment_type), &
      INTENT(in)                             :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "init_scf_run", &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, output_unit
!    INTEGER                                  :: homo, ispin, nao, &
!                                                ndep, nelectron_spin, &
!                                                stat
    LOGICAL                                  :: failure
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    LOGICAL                                  :: ionode
!    TYPE(atomic_kind_type), DIMENSION(:), &
!      POINTER                                :: atomic_kind_set
!    TYPE(cp_fm_type), POINTER                :: mo_coeff, ortho, &
!                                                scf_work2
!    TYPE(mo_set_p_type), DIMENSION(:), &
!      POINTER                                :: c
!    TYPE(real_matrix_p_type), DIMENSION(:), &
!      POINTER                                :: h, k, s
!    TYPE(scf_control_type), POINTER          :: scf_control
!
!   ---------------------------------------------------------------------------

    CALL timeset("init_scf_run","I"," ",handle)

    NULLIFY(dft_control)
    NULLIFY(qs_env)

    failure=.FALSE.
    logger => cp_error_get_logger(error)
    
    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(kg_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(kg_env%ref_count>0,cp_failure_level,routineP,error,failure)
    
    CALL get_kg_env(kg_env=kg_env,&
                    dft_control=dft_control)

    IF ( dft_control % qs_control % method == "KG_GPW" ) THEN

      CALL get_kg_env(kg_env=kg_env,&
                      sub_qs_env=qs_env)
      CALL scf_env_check_i_alloc(scf_env=scf_env, qs_env=qs_env,&
           error=error)

      CALL get_qs_env(qs_env=qs_env,&
           dft_control=dft_control,&
           error=error)
!         scf_control=scf_control,&
!         atomic_kind_set=atomic_kind_set,&
!         c=c,k=k,&
!         h=h,s=s, &
!    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
!         gth_potential_present=gth_potential_present)
!
      ionode = logger%para_env%source == logger%para_env%mepos
      IF (ionode) output_unit = cp_logger_get_default_unit_nr(logger)
!
!    scf_work2 => scf_env%scf_work2
!    ortho => scf_env%ortho
!
      ! update ecore
      CALL calculate_ecore_self(qs_env)
      CALL calculate_ecore_overlap(qs_env,globenv,dft_control%forces)
      IF (ionode.AND.logger%print_keys%scf_energies) THEN
         WRITE (UNIT=output_unit,FMT="(/,(T3,A,T60,F20.10))")&
              "Overlap energy of the core charge distribution:",&
              qs_env%energy%core_overlap,&
              "Self energy of the core charge distribution:   ",&
              qs_env%energy%core_self
      END IF
!
!! updates the total number of electrons
!    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
!         nelectron=scf_env%nelectron)
!    scf_env%nelectron = scf_env%nelectron - dft_control%charge
!
!! calc ortho matrix
!    ndep = 0
!    IF (scf_env%method /= ot_method_nr) THEN
!       IF (scf_control%use_cholesky ) THEN
!          CALL copy_sm_to_fm(s(1)%matrix,ortho)
!          CALL cp_fm_cholesky_decompose(ortho)
!       ELSE
!          CALL calculate_ortho_matrix(ortho,s(1)%matrix,scf_work2,ndep,&
!               scf_control%eps_eigval,scf_control%work_syevx)
!          IF (ionode.and.logger%print_keys%scf) THEN
!             WRITE (UNIT=output_unit,FMT="(T3,A,I10)")&
!                  "Number of independent orbital functions:",nao - ndep
!          END IF
!       ENDIF
!    ENDIF
!    IF (logger%print_keys%ortho_matrix) THEN
!       CALL write_fm_with_basis_info(ortho,4,6,qs_env,globenv)
!    END IF
!
!! initializes rho and the mos
!    CALL scf_env_initial_rho_setup(scf_env,qs_env=qs_env,&
!         globenv=globenv,error=error)
!
!! *** method dependent initializations ***
!    SELECT CASE (scf_env%method)
!    CASE (general_diag_method_nr)
!       IF (.NOT.ASSOCIATED(scf_env%scf_diis_buffer)) THEN       
!          CALL qs_diis_b_create(scf_env%scf_diis_buffer,&
!               nbuffer=scf_control%max_diis,&
!               error=error)
!       END IF
!       CALL qs_diis_b_clear(scf_env%scf_diis_buffer,error=error)
!    CASE DEFAULT
!    END SELECT

   ELSE

     CPErrorMessage(-101,routineP,"Illegal Method",error)

   END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE init_scf_run

! *****************************************************************************
!!****f* kg_scf/scf_env_check_i_alloc [1.0] *
!!
!!   NAME
!!     scf_env_check_i_alloc
!!
!!   AUTHOR
!!     fawzi/jgh
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************
  SUBROUTINE scf_env_check_i_alloc(scf_env,qs_env,error)
    TYPE(kg_scf_env_type), POINTER           :: scf_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "scf_env_check_i_alloc", &
      routineP = moduleN//":"//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(cp_sm_pool_p_type), DIMENSION(:), &
      POINTER                                :: S_sm_pools
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: k, s

!   ---------------------------------------------------------------------------

    CALL timeset(routineN,"I"," ",handle)

    NULLIFY(k,ks_env,S_sm_pools)
    failure=.FALSE.

    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)

    CALL get_qs_env(qs_env=qs_env,&
         matrix_ks=k,&
         ks_env=ks_env)
    CALL mpools_get(qs_env%mpools,&
         S_sm_pools=S_sm_pools,error=error)

!! also allocate the auxilliary matrices needed for the scf calculation
!! -> definition of kg_scf_type

    ! Allocate k and put it in the QS environment
    IF (.NOT.ASSOCIATED(k)) THEN
       CALL sm_pools_create_matrix_vect(S_sm_pools,k,&
            name="SCF"//TRIM(ADJUSTL(cp_to_string(scf_env%id_nr)))//&
            "KOHN-SHAM_MATRIX-",error=error)
       CALL set_qs_env(qs_env=qs_env,matrix_ks=k)
    END IF

    ! allocate the ks env
    IF (.NOT.ASSOCIATED(ks_env)) THEN
       CALL qs_ks_create(ks_env,qs_env=qs_env)
       CALL set_qs_env(qs_env, ks_env=ks_env)
       CALL qs_ks_release(ks_env)
    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE scf_env_check_i_alloc

!***************************************************************************
!  SUBROUTINE calculate_ortho_matrix(ortho,s,work,ndep,eps_eigval,&
!       work_syevx)
!
!!   Purpose: Calculate the orthogonalization matrix (S**(-1/2))
!
!!   History: - Creation (01.05.2001, Matthias Krack)
!
!!   ***************************************************************************
!
!    TYPE(cp_fm_type), POINTER                :: ortho
!    TYPE(real_matrix_type), POINTER          :: s
!    TYPE(cp_fm_type), POINTER                :: work
!    INTEGER, INTENT(OUT)                     :: ndep
!    REAL(KIND=dp), INTENT(in)                :: eps_eigval, work_syevx
!
!    INTEGER                                  :: handle
!
!!   ---------------------------------------------------------------------------
!
!    CALL timeset("calculate_ortho_matrix","I"," ",handle)
!
!    CALL copy_sm_to_fm(s,ortho)
!    CALL cp_fm_power(ortho,work,-0.5_dp,eps_eigval,ndep,work_syevx)
!    CALL cp_fm_upper_to_full(ortho,work)
!
!    CALL timestop(0.0_dp,handle)
!
!  END SUBROUTINE calculate_ortho_matrix
!
!! *****************************************************************************
!
!  SUBROUTINE eigensolver(ks,mo_set,ortho,work,do_level_shift,level_shift,&
!       use_cholesky, work_syevx, use_jacobi,jacobi_threshold,smear,&
!       error)
!
!!   Purpose: Diagonalise the Kohn-Sham matrix to get a new set of MO eigen-
!!            vectors and MO eigenvalues.
!!   ks will be modified
!
!!   History: - Creation (01.05.2001, Matthias Krack)
!
!!   ***************************************************************************
!
!    TYPE(cp_fm_type), POINTER                :: ks
!    TYPE(mo_set_type), POINTER               :: mo_set
!    TYPE(cp_fm_type), POINTER                :: ortho, work
!    LOGICAL, INTENT(IN)                      :: do_level_shift
!    REAL(KIND=dp), INTENT(IN)                :: level_shift
!    LOGICAL, INTENT(INOUT)                   :: use_cholesky
!    REAL(KIND=dp), INTENT(IN)                :: work_syevx
!    LOGICAL, INTENT(IN)                      :: use_jacobi
!    REAL(KIND=dp), INTENT(IN)                :: jacobi_threshold, smear
!    TYPE(cp_error_type), INTENT(inout), &
!      OPTIONAL                               :: error
!
!    INTEGER                                  :: handle, homo, imo, nao, nmo
!    REAL(KIND=dp)                            :: alpha
!    REAL(KIND=dp), DIMENSION(:), POINTER     :: mo_eigenvalues
!    TYPE(cp_fm_type), POINTER                :: mo_coeff
!
!!   ---------------------------------------------------------------------------
!
!    CALL timeset("eigensolver","I"," ",handle)
!
!!   *** Diagonalise the Kohn-Sham matrix ***
!
!    CALL get_mo_set(mo_set=mo_set,&
!         nao=nao,&
!         nmo=nmo,&
!         homo=homo, &
!         eigenvalues=mo_eigenvalues,&
!         mo_coeff=mo_coeff)
!
!    IF (do_level_shift) THEN
!       IF (use_cholesky) THEN
!          CALL stop_program("eigensolver","level shift not implemented")
!       ENDIF
!    ENDIF
!
!    IF (use_cholesky) THEN
!       ! CALL cp_fm_write(ks,6,long_description=.true.)
!       CALL cp_fm_cholesky_reduce(ks,ortho)
!
!       IF (use_jacobi) THEN
!          CALL cp_fm_syevd(ks,work,mo_eigenvalues)
!          use_cholesky = .FALSE.
!          CALL cp_fm_cholesky_restore(work,nmo,ortho,mo_coeff,"SOLVE")
!          CALL cp_fm_to_fm(mo_coeff,ortho)
!       ELSE IF (smear /= 0.0_dp) THEN
!          CALL cp_fm_syevd(ks,work,mo_eigenvalues)
!          CALL cp_fm_cholesky_restore(work,nmo,ortho,mo_coeff,"SOLVE")
!       ELSE
!          CALL cp_fm_syevx(ks,work,mo_eigenvalues,homo,work_syevx)
!          CALL cp_fm_cholesky_restore(work,homo,ortho,mo_coeff,"SOLVE")
!       END IF
!
!    ELSE IF (use_jacobi) THEN
!
!       CALL cp_fm_symm("L","U",nao,homo,1.0_dp,ks,mo_coeff,0.0_dp,work)
!       CALL cp_fm_gemm("T","N",homo,nao-homo,nao,1.0_dp,work,mo_coeff,&
!            0.0_dp,ks,b_first_col=homo+1)
!
!!     epsi =1e-11
!!     Anz_max=4
!!     N=ks%matrix_struct%nrow_global
!
!!     *** Scale the elements with the inverse ***
!!     *** difference of the old eigenvalues   ***
!
!!     CALL blacs_block_jacobi_scaled(ks,work,mo_eigenvalues,&
!!                                    jacobi_threshold,homo+1)
!
!!     *** Klassisches Jacobi-Block-Verfrahren ***
!!     *** Rotiere Elemente nur, wenn sie groesser als thresh sind ***
!
!       CALL cp_fm_block_jacobi_classic(ks,mo_coeff,mo_eigenvalues,&
!            jacobi_threshold,homo+1)
!
!    ELSE
!
!       CALL cp_fm_symm("L","U",nao,nao,1.0_dp,ks,ortho,0.0_dp,work)
!       CALL cp_fm_gemm("T","N",nao,nao,nao,1.0_dp,ortho,work,0.0_dp,ks)
!
!       IF (do_level_shift) THEN
!          DO imo=homo+1,nmo
!             ! now unnecessary broadcast of element
!             CALL cp_fm_get_element(ks,imo,imo,alpha)
!             alpha=alpha+level_shift
!             CALL cp_fm_set_element(ks,imo,imo,alpha)
!          END DO
!       END IF
!
!       CALL cp_fm_syevx(ks,work,mo_eigenvalues,nmo,work_syevx)
!
!       CALL cp_fm_gemm("N","N",nao,nmo,nao,1.0_dp,ortho,work,0.0_dp,&
!            mo_coeff)
!
!       IF (do_level_shift) THEN
!
!          CALL correct_mo_eigenvalues(mo_set,level_shift)
!
!!     *** Use last MO set as orthogonalization matrix ***
!
!          CALL cp_fm_to_fm(mo_coeff,ortho)
!
!       END IF
!
!    END IF
!
!    IF (smear /= 0.0_dp) CALL set_mo_occupation(mo_set,smear)
!
!    CALL timestop(0.0_dp,handle)
!
!  END SUBROUTINE eigensolver
!
!! *****************************************************************************
!
!!!****f* qs_scf/cp_sm_mix [1.0] *
!!!
!!!   NAME
!!!     cp_sm_mix
!!!
!!!   FUNCTION
!!!     Perform a mixing of the given matrixes into the first matrix
!!!     m1 = m2 + p_mix (m1-m2)
!!!
!!!   NOTES
!!!     if you what to store the result in m2 swap m1 and m2 an use 
!!!     (1-pmix) as pmix
!!!     para_env should be removed (embedded in matrix)
!!!
!!!   ARGUMENTS
!!!     - m1: first (new) matrix, is modified
!!!     - m2: the second (old) matrix
!!!     - p_mix:how much m1 is conserved (0: none, 1: all)
!!!     - delta: maximum norm of m1-m2
!!!     - error: variable to control error logging, stopping,... 
!!!       see module cp_error_handling 
!!!
!!!   AUTHOR
!!!     fawzi
!!!
!!!   MODIFICATION HISTORY
!!!     02.2003 rewamped [fawzi]
!!!
!!!*** *********************************************************************
!  SUBROUTINE cp_sm_mix(m1,m2,p_mix,delta,para_env,error)
!
!    TYPE(real_matrix_type), POINTER          :: m1, m2
!    REAL(KIND=dp), INTENT(IN)                :: p_mix
!    REAL(KIND=dp), INTENT(OUT)               :: delta
!    TYPE(cp_para_env_type), POINTER          :: para_env
!    TYPE(cp_error_type), INTENT(inout), &
!      OPTIONAL                               :: error
!
!    CHARACTER(len=*), PARAMETER :: routineN = "cp_sm_mix", &
!      routineP = moduleN//':'//routineN
!
!    INTEGER                                  :: handle, i, iblock_col, &
!                                                iblock_row, j, nblock_row
!    REAL(KIND=dp)                            :: r_mix
!    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: p_new_block, p_old_block
!    TYPE(real_block_node_type), POINTER      :: block_node
!
!!   ---------------------------------------------------------------------------
!
!    CALL timeset(routineN,"I"," ",handle)
!
!    delta = 0.0_dp
!
!    r_mix = 1.0_dp - p_mix
!
!    CALL get_matrix_info(matrix=m2,nblock_row=nblock_row)
!
!    DO iblock_row=1,nblock_row
!
!       block_node => first_block_node(matrix=m1,&
!            block_row=iblock_row)
!
!       DO WHILE (ASSOCIATED(block_node))
!
!          CALL get_block_node(block_node=block_node,&
!               block_col=iblock_col,&
!               BLOCK=p_new_block)
!
!          CALL get_block_node(matrix=m2,&
!               block_row=iblock_row,&
!               block_col=iblock_col,&
!               BLOCK=p_old_block)
!
!          DO j=1,SIZE(p_new_block,2)
!             DO i=1,SIZE(p_new_block,1)
!                p_new_block(i,j) = p_new_block(i,j) - p_old_block(i,j)
!                delta = MAX(delta,ABS(p_new_block(i,j)))
!                p_new_block(i,j) = p_old_block(i,j) + p_mix*p_new_block(i,j)
!             END DO
!          END DO
!
!          block_node => next_block_node(block_node)
!
!       END DO
!
!    END DO
!
!    CALL mp_max(delta,para_env%group)
!
!    CALL timestop(0.0_dp,handle)
!
!  END SUBROUTINE cp_sm_mix
!
!! *****************************************************************************
!  SUBROUTINE calculate_first_density_matrix(scf_env,qs_env,globenv,error)
!
!    TYPE(qs_scf_env_type), POINTER           :: scf_env
!    TYPE(qs_environment_type), POINTER       :: qs_env
!    TYPE(global_environment_type), &
!      INTENT(IN)                             :: globenv
!    TYPE(cp_error_type), INTENT(inout), &
!      OPTIONAL                               :: error
!
!    CHARACTER(LEN=*), PARAMETER :: &
!      routineN = "calculate_first_density_matrix", &
!      routineP = moduleN//":"//routineN
!
!    INTEGER :: atom_a, group, handle, homo, iatom, ikind, iset, isgf, isgfa, &
!      ishell, ispin, istat, la, maxl, nao, natom, ncount, nelectron, nmo, &
!      nset, nspin, qs_env_id
!    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_sgf
!    INTEGER, DIMENSION(:), POINTER           :: atom_list, elec_conf, nshell
!    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, l, last_sgfa
!    REAL(KIND=dp)                            :: maxocc, nelec, paa, scale, &
!                                                trps1, trps2
!    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: econf, pdiag, sdiag
!    TYPE(all_potential_type), POINTER        :: all_potential
!    TYPE(atomic_kind_type), DIMENSION(:), &
!      POINTER                                :: atomic_kind_set
!    TYPE(atomic_kind_type), POINTER          :: atomic_kind
!    TYPE(cp_fm_p_type), DIMENSION(:), &
!      POINTER                                :: work1
!    TYPE(cp_fm_type), POINTER                :: mo_coeff, ortho, sv, &
!                                                work2
!    TYPE(dft_control_type), POINTER          :: dft_control
!    TYPE(gth_potential_type), POINTER        :: gth_potential
!    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
!    TYPE(mo_set_p_type), DIMENSION(:), &
!      POINTER                                :: mo_array
!    TYPE(particle_type), DIMENSION(:), &
!      POINTER                                :: particle_set
!    TYPE(real_matrix_p_type), DIMENSION(:), &
!      POINTER                                :: p_rmpv, s_sparse
!    TYPE(scf_control_type), POINTER          :: scf_control
!
!!   ---------------------------------------------------------------------------
!
!    NULLIFY(all_potential,atomic_kind, mo_coeff, sv,&
!         gth_potential, orb_basis_set, atomic_kind_set, particle_set,&
!         ortho,work2,work1,mo_array,s_sparse,p_rmpv,scf_control, &
!         dft_control)
!
!    CALL timeset("calculate_first_density_matrix","I"," ",handle)
!
!    CALL get_qs_env(qs_env,atomic_kind_set=atomic_kind_set,&
!         particle_set=particle_set,c=mo_array, s=s_sparse,&
!         scf_control=scf_control, id_nr=qs_env_id, dft_control=dft_control,&
!         error=error)
!    nspin=dft_control%nspins
!
!    p_rmpv => qs_env%rho%rho_ao
!    work1 => scf_env%scf_work1
!    work2 => scf_env%scf_work2
!    ortho => scf_env%ortho
!
!    IF (scf_control%use_ot .AND. .NOT.  &
!         (scf_control%density_guess == "RANDOM" .OR.  &
!         (scf_control%density_guess == "RESTART" .AND. &
!         scf_control%level_shift == 0.0_dp))) THEN
!       CALL stop_program("calculate_first_density_matrix","OT needs GUESS RESTART : other options NYI")
!    ENDIF
!
!    IF (scf_control%density_guess == "RESTART") THEN
!
!       CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set)
!       ! CALL read_mo_set(mo_array,atomic_kind_set,particle_set,globenv,&
!       !      id_nr=qs_env_id)
!       CALL read_mo_set(mo_array,atomic_kind_set,particle_set,globenv,&
!            id_nr=1)
!
!       DO ispin=1,nspin
!          IF (scf_control%level_shift /= 0.0_dp) THEN
!             CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
!                  mo_coeff=mo_coeff)
!             CALL cp_fm_to_fm(mo_coeff,ortho)
!          END IF
!
!          CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
!               mo_coeff=mo_coeff, homo=homo)
!! ortho so that one can restart for different positions (basis sets?)
!          CALL cp_fm_create(sv,mo_coeff%matrix_struct,"SV")
!          CALL cp_sm_fm_multiply(s_sparse(1)%matrix,mo_coeff,sv,homo)
!          CALL make_basis(mo_coeff,homo, &
!               ortho=sv,otype="SV")
!          CALL cp_fm_release(sv)
!          CALL calculate_density_matrix(mo_array(ispin)%mo_set,&
!               p_rmpv(ispin)%matrix,error=error)
!       ENDDO
!
!    ELSE IF (scf_control%density_guess == "RANDOM") THEN
!
!       DO ispin=1,nspin
!          CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
!               mo_coeff=mo_coeff, nmo=nmo)
!          CALL cp_fm_init_random(mo_coeff,nmo)
!! ortho so that one can restart for different positions (basis sets?)
!          CALL cp_fm_create(sv,mo_coeff%matrix_struct,"SV")
!          CALL cp_sm_fm_multiply(s_sparse(1)%matrix,mo_coeff,sv,nmo)
!          CALL make_basis(mo_coeff,nmo, &
!               ortho=sv,otype="SV")
!          CALL cp_fm_release(sv)
!          CALL calculate_density_matrix(mo_array(ispin)%mo_set,&
!               p_rmpv(ispin)%matrix,error=error)
!       ENDDO
!
!    ELSE IF (scf_control%density_guess == "CORE") THEN
!
!! work1(1) should contain the core hamiltonian
!! we need a copy the core matrix for every spin (the dumbest solution,
!! we should copy the mos)..
!       IF (dft_control%nspins.eq.2) THEN
!          IF (.NOT. scf_control%use_ot) THEN
!             CALL cp_fm_to_fm(work1(1)%matrix,&
!                  work1(2)%matrix)
!          ENDIF
!       ENDIF
!       DO ispin=1,nspin
!          CALL eigensolver(work1(ispin)%matrix,mo_array(ispin)%mo_set,ortho,work2, &
!               .FALSE.,0.0_dp,&
!               use_cholesky=scf_control%use_cholesky, &
!               work_syevx=scf_control%work_syevx,&
!               use_jacobi=.FALSE.,&
!               jacobi_threshold=scf_control%jacobi_threshold,&
!               smear=0.0_dp,&
!               error=error)
!          IF (scf_control%level_shift /= 0.0_dp) THEN
!             CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
!                  mo_coeff=mo_coeff)
!             CALL cp_fm_to_fm(mo_coeff,ortho)
!          END IF
!          CALL calculate_density_matrix(mo_array(ispin)%mo_set,&
!               p_rmpv(ispin)%matrix,error=error)
!       ENDDO
!
!    ELSE IF (scf_control%density_guess == "ATOMIC") THEN
!
!      group = qs_env%para_env%group
!
!      natom = SIZE(particle_set)
!      ALLOCATE (first_sgf(natom),STAT=istat)
!      IF (istat /= 0) CALL stop_memory(routineP,"first_sgf",natom*int_size)
!      CALL get_particle_set(particle_set=particle_set,&
!                            first_sgf=first_sgf)
!
!      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,maxlgto=maxl)
!      ALLOCATE (econf(0:maxl),STAT=istat)
!      IF (istat /= 0) CALL stop_memory(routineP,"econf",(maxl + 1)*int_size)
!
!      CALL get_mo_set(mo_array(1)%mo_set,nao=nao)
!      ALLOCATE (pdiag(nao),STAT=istat)
!      IF (istat /= 0) CALL stop_memory(routineP,"pdiag",nao*dp_size)
!      pdiag(:) = 0.0_dp
!
!      ALLOCATE (sdiag(nao),STAT=istat)
!      IF (istat /= 0) CALL stop_memory(routineP,"sdiag",nao*dp_size)
!      CALL get_matrix_diagonal(s_sparse(1)%matrix,sdiag)
!      CALL mp_sum(sdiag,group)
!
!      DO ispin=1,nspin
!
!        CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
!                        maxocc=maxocc,&
!                        nelectron=nelectron)
!
!        ncount = 0
!        trps1 = 0.0_dp
!        trps2 = 0.0_dp
!        pdiag(:) = 0.0_dp
!
!        DO ikind=1,SIZE(atomic_kind_set)
!
!          atomic_kind => atomic_kind_set(ikind)
!
!          CALL get_atomic_kind(atomic_kind=atomic_kind,&
!                               natom=natom,&
!                               atom_list=atom_list,&
!                               all_potential=all_potential,&
!                               gth_potential=gth_potential,&
!                               orb_basis_set=orb_basis_set)
!
!          IF (ASSOCIATED(all_potential)) THEN
!            CALL get_potential(potential=all_potential,elec_conf=elec_conf)
!          ELSE IF (ASSOCIATED(gth_potential)) THEN
!            CALL get_potential(potential=gth_potential,elec_conf=elec_conf)
!          ELSE
!            CYCLE
!          END IF
!
!          maxl = SIZE(elec_conf) - 1
!          econf(:) = 0.0_dp
!          econf(0:maxl) = 0.5_dp*maxocc*REAL(elec_conf(0:maxl),dp)
!
!          CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
!                                 nset=nset,&
!                                 nshell=nshell,&
!                                 l=l,&
!                                 first_sgf=first_sgfa,&
!                                 last_sgf=last_sgfa)
!
!          DO iset=1,nset
!            DO ishell=1,nshell(iset)
!              la = l(ishell,iset)
!              nelec = maxocc*REAL(2*la + 1,dp)
!              IF (econf(la) > 0.0_dp) THEN
!                IF (econf(la) >= nelec) THEN
!                  paa = maxocc
!                  econf(la) = econf(la) - nelec
!                ELSE
!                  paa = maxocc*econf(la)/nelec
!                  econf(la) = 0.0_dp
!                  ncount = ncount + NINT(nelec/maxocc)
!                END IF
!                DO isgfa=first_sgfa(ishell,iset),last_sgfa(ishell,iset)
!                  DO iatom=1,natom
!                    atom_a = atom_list(iatom)
!                    isgf = first_sgf(atom_a) + isgfa - 1
!                    pdiag(isgf) = paa
!                    IF (paa == maxocc) THEN
!                      trps1 = trps1 + paa*sdiag(isgf)
!                    ELSE
!                      trps2 = trps2 + paa*sdiag(isgf)
!                    END IF
!                  END DO
!                END DO
!              END IF
!            END DO
!          END DO
!
!        END DO
!
!        IF (trps2 == 0.0_dp) THEN
!          DO isgf=1,nao
!            IF (sdiag(isgf) > 0.0_dp) pdiag(isgf) = pdiag(isgf)/sdiag(isgf)
!          END DO
!        ELSE
!          scale = (REAL(nelectron,dp) - trps1)/trps2
!          DO isgf=1,nao 
!            IF (pdiag(isgf) < maxocc) pdiag(isgf) = scale*pdiag(isgf)
!          END DO
!        END IF
!
!        CALL set_matrix_diagonal(p_rmpv(ispin)%matrix,pdiag)
!
!      END DO
!
!      DEALLOCATE (econf,STAT=istat)
!      IF (istat /= 0) CALL stop_memory(routineP,"econf")
!
!      DEALLOCATE (first_sgf,STAT=istat)
!      IF (istat /= 0) CALL stop_memory(routineP,"first_sgf")
!
!      DEALLOCATE (pdiag,STAT=istat)
!      IF (istat /= 0) CALL stop_memory(routineP,"pdiag")
!
!      DEALLOCATE (sdiag,STAT=istat)
!      IF (istat /= 0) CALL stop_memory(routineP,"sdiag")
!
!    ELSE
!
!      CALL stop_program(routineP,&
!                        "An invalid keyword for the initial density "//&
!                        "guess was specified")
!
!    END IF
!
!    CALL timestop(0.0_dp,handle)
!
!  END SUBROUTINE calculate_first_density_matrix
!
!!!****f* qs_scf/scf_env_initial_rho_setup [1.0] *
!!!
!!!   NAME
!!!     scf_env_initial_rho_setup
!!!
!!!   FUNCTION
!!!     Initializes rho and the mos, so that an scf cycle can start
!!!
!!!   NOTES
!!!     -
!!!
!!!   ARGUMENTS
!!!     - scf_env: the scf env in which to do the scf
!!!     - qs_env: the qs env the scf_env lives in
!!!     - error: variable to control error logging, stopping,... 
!!!       see module cp_error_handling 
!!!
!!!   AUTHOR
!!!     fawzi
!!!
!!!   MODIFICATION HISTORY
!!!     02.2003 created [fawzi]
!!!
!!!*** **********************************************************************
!  SUBROUTINE scf_env_initial_rho_setup(scf_env, qs_env, globenv, error)
!    TYPE(qs_scf_env_type), POINTER           :: scf_env
!    TYPE(qs_environment_type), POINTER       :: qs_env
!    TYPE(global_environment_type), &
!      INTENT(in), TARGET                     :: globenv
!    TYPE(cp_error_type), INTENT(inout), &
!      OPTIONAL                               :: error
!
!    CHARACTER(len=*), PARAMETER :: routineN = 'scf_env_initial_rho_setup', &
!      routineP = moduleN//':'//routineN
!
!    INTEGER                                  :: extrapolation_method_nr, &
!                                                handle, ispin, nmo
!    LOGICAL                                  :: failure, ionode, orthogonal_wf
!    TYPE(cp_fm_type), POINTER                :: mo_coeff
!    TYPE(cp_logger_type), POINTER            :: logger
!
!    CALL timeset(routineN,"I"," ",handle)
!    failure=.FALSE.
!    NULLIFY(mo_coeff)
!    logger => cp_error_get_logger(error)
!    ionode = logger%para_env%source==logger%para_env%mepos
!    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
!    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
!    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
!    IF (.NOT. failure) THEN
!       extrapolation_method_nr=wfi_use_guess_method_nr
!       IF (ASSOCIATED(qs_env%wf_history)) THEN
!          CALL wfi_extrapolate(globenv,qs_env%wf_history, &
!               qs_env=qs_env, dt=1.0_dp, &
!               extrapolation_method_nr=extrapolation_method_nr,&
!               orthogonal_wf=orthogonal_wf, error=error)
!          ! wfi_use_guess_method_nr the wavefunctions are not yet initialized
!          IF (.NOT.orthogonal_wf.AND.scf_env%method==ot_method_nr.AND. &
!                                .NOT.extrapolation_method_nr==wfi_use_guess_method_nr) THEN
!             DO ispin=1,SIZE(qs_env%c)
!                CALL get_mo_set(qs_env%c(ispin)%mo_set, &
!                     mo_coeff=mo_coeff, nmo=nmo)
!                CALL qs_env_reorthogonalize_vectors(qs_env, &
!                     v_matrix=mo_coeff, ispin=ispin, n_col=nmo,&
!                     error=error)
!             END DO
!          END IF
!       END IF
!       IF (ionode) &
!       WRITE (UNIT=cp_logger_get_default_unit_nr(logger),FMT="(/,T3,A)")&
!         "Extrapolation method: "//&
!         TRIM(wfi_get_method_label(extrapolation_method_nr,error=error))
!       IF (extrapolation_method_nr==wfi_use_guess_method_nr) THEN
!          CALL calculate_first_density_matrix(scf_env=scf_env,qs_env=qs_env,&
!               globenv=globenv, error=error)
!          CALL qs_rho_update_rho(qs_env%rho,qs_env=qs_env, error=error)
!          CALL qs_ks_did_change(qs_env%ks_env,rho_changed=.TRUE.,&
!               error=error)
!       END IF
!    END IF
!
!    CALL timestop(0.0_dp,handle)
!  END SUBROUTINE scf_env_initial_rho_setup
!!***************************************************************************
!
!!!****f* qs_scf/scf_env_cleanup [1.0] *
!!!
!!!   NAME
!!!     scf_env_cleanup
!!!
!!!   FUNCTION
!!!     perform cleanup operations (like releasing temporary storage)
!!!     at the end of the scf
!!!
!!!   NOTES
!!!     -
!!!
!!!   ARGUMENTS
!!!     - error: variable to control error logging, stopping,... 
!!!       see module cp_error_handling 
!!!
!!!   AUTHOR
!!!     fawzi
!!!
!!!   MODIFICATION HISTORY
!!!     02.2003 created [fawzi]
!!!
!!!*** **********************************************************************
!  SUBROUTINE scf_env_cleanup(scf_env,qs_env,error)
!    TYPE(qs_scf_env_type), POINTER           :: scf_env
!    TYPE(qs_environment_type), POINTER       :: qs_env
!    TYPE(cp_error_type), INTENT(inout), &
!      OPTIONAL                               :: error
!
!    CHARACTER(len=*), PARAMETER :: routineN = 'scf_env_cleanup', &
!      routineP = moduleN//':'//routineN
!
!    INTEGER                                  :: handle, ispin, stat
!    LOGICAL                                  :: failure
!    TYPE(cp_sm_pool_p_type), DIMENSION(:), &
!      POINTER                                :: S_sm_pools
!
!    CALL timeset(routineN,"I"," ",handle)
!
!    failure=.FALSE.
!    NULLIFY(S_sm_pools)
!
!    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
!    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
!    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
!    IF (.NOT. failure) THEN
!       CALL mpools_get(qs_env%mpools,S_sm_pools=S_sm_pools,error=error)
!
!!   *** Release SCF work storage ***
!
!       IF (ASSOCIATED(scf_env%scf_work1)) THEN
!          DO ispin=1,SIZE(scf_env%scf_work1)
!             CALL cp_fm_release(scf_env%scf_work1(ispin)%matrix)
!          ENDDO
!          DEALLOCATE(scf_env%scf_work1)
!       ENDIF
!       IF (ASSOCIATED(scf_env%scf_work2)) CALL cp_fm_release(scf_env%scf_work2)
!       IF (ASSOCIATED(scf_env%ortho)) CALL cp_fm_release(scf_env%ortho)
!
!       IF (ASSOCIATED(scf_env%p_mix_new)) THEN
!          CALL sm_pools_give_back_matrix_vect(S_sm_pools,scf_env%p_mix_new,&
!               error=error)
!          CALL sm_pools_flush_cache(S_sm_pools,error=error)
!       END IF
!
!! *** method dependent cleanup
!       SELECT CASE(scf_env%method)
!       CASE(ot_method_nr)
!          CALL destroy_preconditioner(scf_env%ot_preconditioner)
!          DEALLOCATE(scf_env%ot_preconditioner,stat=stat)
!          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
!          DO ispin=1,SIZE(scf_env%qs_ot_env)
!             CALL ot_scf_destroy(scf_env%qs_ot_env(ispin))
!          ENDDO
!          DEALLOCATE(scf_env%qs_ot_env,stat=stat)
!          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
!       CASE(pao_mehod_nr)
!!    call pao_opt_dealloc(pao_optimizer,pao_obj_f,qs_env,globenv,error)
!       CASE(general_diag_method_nr)
!          !
!       CASE default
!          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
!               routineP,"unknown scf method method:"//&
!               cp_to_string(scf_env%method),error,failure)
!       END SELECT
!
!    END IF
!    CALL timestop(0.0_dp,handle)
!  END SUBROUTINE scf_env_cleanup
!!***************************************************************************
!
!!!****f* qs_scf/scf_env_did_change [1.0] *
!!!
!!!   NAME
!!!     scf_env_did_change
!!!
!!!   SYNOPSIS
!!!     Subroutine scf_env_did_change(scf_env, qs_env, s_mstruct_changed,&
!!!         error)
!!!       Type(qs_scf_env_type), Pointer:: scf_env
!!!       Type(qs_environment_type), Pointer:: qs_env
!!!       Logical, Intent (IN):: s_mstruct_changed
!!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!!     End Subroutine scf_env_did_change
!!!
!!!   FUNCTION
!!!     function to be called to inform the scf_env about changes
!!!
!!!   NOTES
!!!     -
!!!
!!!   ARGUMENTS
!!!     - scf_env: the scf env to inform
!!!     - qs_env: the qs_env in which the scf_env lives
!!!     - s_mstruct_changed: treu if the structure of s changed
!!!     - error: variable to control error logging, stopping,... 
!!!       see module cp_error_handling 
!!!
!!!   AUTHOR
!!!     fawzi
!!!
!!!   MODIFICATION HISTORY
!!!     03.2003 created [fawzi]
!!!
!!!*** **********************************************************************
!SUBROUTINE scf_env_did_change(scf_env,error)
!    TYPE(qs_scf_env_type), POINTER           :: scf_env
!    TYPE(cp_error_type), INTENT(inout), &
!      OPTIONAL                               :: error
!
!    CHARACTER(len=*), PARAMETER :: routineN = 'scf_env_did_change', &
!      routineP = moduleN//':'//routineN
!
!    LOGICAL                                  :: failure
!
!  failure=.FALSE.
!  
!  CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
!  CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
!  IF (.NOT. failure) THEN
!     IF (ASSOCIATED(scf_env%p_mix_new)) THEN
!        CALL deallocate_matrix_set(scf_env%p_mix_new)
!     END IF
!  END IF
!END SUBROUTINE scf_env_did_change
!
!***************************************************************************

END MODULE kg_scf
