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

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
!                                             get_atomic_kind,&
                                             get_atomic_kind_set
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
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_sm_pool_types,                ONLY: cp_sm_pool_p_type,&
                                             sm_pools_create_matrix_vect,&
                                             sm_pools_flush_cache,&
                                             sm_pools_give_back_matrix_vect
!  USE external_potential_types,        ONLY: all_potential_type,&
!                                             get_potential,&
!                                             gth_potential_type
  USE global_types,                    ONLY: global_environment_type
  USE kg_diis_methods,                 ONLY: kg_diis_b_clear,&
                                             kg_diis_b_create, &
                                             kg_diis_step_A,&
                                             kg_diis_step_B
  USE kg_gpw_fm_mol_methods,           ONLY: calculate_mol_density_matrix,&
                                             fm_mol_eigensolver,&
                                             fm_mol_orthonormality,&
                                             mol_make_basis,&
                                             multiply_sparse_mol_mo
  USE kg_gpw_fm_mol_types,             ONLY: fm_mol_blocks_type,&
                                             get_kg_fm_mol_set,&
                                             get_fm_mol_block,&
                                             get_mol_mo_set,&
                                             kg_fm_mol_set_type,& 
                                             kg_fm_p_type,&
                                             mol_mo_set_p_type,&
                                             mol_mo_set_type
  USE kg_gpw_fm_mol_utils,             ONLY: copy_sparse2mol_block,&
                                             fm_mol_power,&
                                             mol_mo_random
  USE kg_environment_types,            ONLY: kg_environment_type,&
                                             get_kg_env,&
                                             set_kg_env
  USE kg_scf_types,                    ONLY: kg_scf_env_type,&
                                             kg_scf_env_release
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             int_size
  USE machine,                         ONLY: m_flush,&
                                             m_walltime
  USE mathlib,                         ONLY: diamat_all
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_max,&
                                             mp_sum
  USE particle_types,                  ONLY: get_particle_set,&
                                             particle_type
  USE preconditioner,                  ONLY: destroy_preconditioner,&
                                             init_preconditioner,&
                                             make_preconditioner
!  USE pw_env_types,                    ONLY: pw_env_type
  USE qs_charges_types,                ONLY: qs_charges_type
  USE qs_core_energies,                ONLY: calculate_ecore_overlap,&
                                             calculate_ecore_self
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_ks_methods,                   ONLY: qs_ks_create,&
                                             qs_ks_did_change,&
                                             qs_ks_update_qs_env
  USE qs_ks_types,                     ONLY: qs_ks_env_type,&
                                             qs_ks_release
  USE qs_matrix_pools,                 ONLY: mpools_get
  USE qs_mo_types,                     ONLY: correct_mo_eigenvalues,&
                                             set_mo_occupation
  USE qs_overlap,                      ONLY: write_fm_with_basis_info,&
                                             write_sparse_matrix
  USE qs_ot,                           ONLY: qs_ot_new_preconditioner, &
                                             qs_ot_type
  USE qs_ot_scf,                       ONLY: ot_scf_destroy,&
                                             ot_scf_init,&
                                             ot_scf_mini,&
                                             ot_scf_read_input
  USE qs_rho_methods,                  ONLY: qs_rho_update_rho
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE qs_scf,                          ONLY: scf_env_density_mixing
  USE qs_scf_post,                     ONLY: scf_post_calculation
  USE qs_scf_types,                    ONLY: ao_mixing_nr,&
                                             general_diag_method_nr,&
                                             special_diag_method_nr,&
                                             no_mixing_nr,&
                                             ot_method_nr
  USE qs_wf_history_methods,           ONLY: wfi_get_method_label
                                             !wfi_extrapolate,&
!                                             wfi_update
  USE qs_wf_history_types,             ONLY: wfi_use_guess_method_nr
  USE scf_control_types,               ONLY: scf_control_type
!  USE simulation_cell,                 ONLY: cell_type
  USE sparse_matrix_types,             ONLY: deallocate_matrix_set,&
                                             get_block_node,&
                                             real_matrix_p_type,&
                                             real_matrix_type
!  USE sparse_matrix_types,             ONLY: 
!                                             first_block_node,&
!                                             get_matrix_diagonal,&
!                                             get_matrix_info,&
!                                             next_block_node,&
!                                             real_block_node_type,&
!                                             set_matrix_diagonal,&
!                                             transfer_matrix
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: moduleN = "kg_scf"

  LOGICAL, PARAMETER, PRIVATE :: debug_this_module=.FALSE.

! *** Global variables ***

  INTEGER, SAVE, PRIVATE :: last_scf_env_id=0

! *** Public subroutines ***

  PUBLIC :: kg_qs_scf, kg_scf_env_did_change

! *****************************************************************************

CONTAINS

!!****f* kg_scf/kg_scf_env_create [1.0] *
!!
!!   NAME
!!     kg_scf_env_create
!!
!!   AUTHOR
!!     fawzi/jgh
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************

  SUBROUTINE kg_scf_env_create(scf_env,kg_env, error)

    TYPE(kg_scf_env_type), POINTER           :: scf_env
    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_scf_env_create', &
      routineP = moduleN//':'//routineN

    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(scf_control_type), POINTER          :: scf_control
    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    NULLIFY(dft_control,scf_control )

    CALL get_kg_env(kg_env=kg_env,&
                    dft_control=dft_control,&
                    scf_control=scf_control)

    ALLOCATE(scf_env, stat=stat)
    CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    IF (.NOT. failure) THEN
       scf_env%ref_count=1
       scf_env%print_count=0
       last_scf_env_id=last_scf_env_id+1
       scf_env%id_nr=last_scf_env_id
       scf_env%print_count=0
       scf_env%iter_count=0
! sets the method
       scf_env%method=general_diag_method_nr ! default with diagonalisation
       IF (dft_control%qs_control%semi_empirical) THEN
!         IF (dft_control%qs_control%se_control%orthogonal_basis) &
!             scf_env%method=special_diag_method_nr
          CALL stop_program(routineN,moduleN,__LINE__,&
               "KG_GPW with semi_empirical not implemented")
       END IF
       IF (scf_control%use_ot) scf_env%method=ot_method_nr

!      I still don t know whether KG_GPW and OT can run together
       SELECT CASE(scf_env%method)
       CASE(ot_method_nr)
          scf_env%mixing_method=no_mixing_nr
       CASE DEFAULT
          scf_env%mixing_method=ao_mixing_nr
       END SELECT

       scf_env%iter_param=0.0_dp
       scf_env%iter_delta=0.0_dp
       scf_env%p_mix=scf_control%p_mix
       scf_env%iter_method=""
       scf_env%print_iter_line=.TRUE.

       NULLIFY( scf_env%p_mix_new, scf_env%ot_preconditioner, scf_env%qs_ot_env,&
            scf_env%scf_diis_buffer)

    END IF

  END SUBROUTINE kg_scf_env_create

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
           CALL kg_scf_env_create(scf_env,kg_env, error=error)
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

    TYPE(kg_scf_env_type), POINTER             :: scf_env
    TYPE(kg_environment_type), POINTER         :: kg_env
    TYPE(global_environment_type), &
      INTENT(IN)                               :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                 :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "scf_env_do_scf", &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                              :: t1, t2
    INTEGER                                    :: handle, handle2
    INTEGER                                    :: output_unit, ispin, istat
    INTEGER                                    :: ib, imol, &
                                                  jb, nat_mol, &
                                                  nb, nb1, nelectron_global,&
                                                  nmolecule_global, nmol_local,&
                                                  nmolecule_kind
    LOGICAL                                    :: energy_only
    LOGICAL                                    :: failure
    LOGICAL                                    :: ionode, id_equal
    LOGICAL                                    :: diis_step, do_level_shift,&
                                                  should_stop, use_cholesky, use_jacobi
!    TYPE(atomic_kind_type), DIMENSION(:), &
!      POINTER                                :: atomic_kind_set
!    TYPE(cell_type), POINTER                 :: cell
    TYPE(qs_environment_type), POINTER         :: qs_env
    TYPE(cp_logger_type), POINTER              :: logger
    TYPE(dft_control_type), POINTER            :: dft_control
    TYPE(scf_control_type), POINTER            :: scf_control
    TYPE(qs_ks_env_type), POINTER              :: ks_env
!    TYPE(particle_type), DIMENSION(:), &
!      POINTER                                :: particle_set
!    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(qs_charges_type), POINTER             :: qs_charges
    TYPE(qs_energy_type), POINTER              :: energy
    TYPE(qs_rho_type), POINTER                 :: rho
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                  :: matrix_h,&
                                                  matrix_ks,&
                                                  matrix_p,&
                                                  matrix_s
    TYPE(kg_fm_p_type), POINTER                :: kg_fm_set
    REAL(dp)                                   :: diis_error
    REAL(dp), DIMENSION(:), POINTER            :: diis_ev_global
    REAL(dp), DIMENSION(:,:), POINTER          :: diis_m1_global,diis_m2_global
    REAL(dp), DIMENSION(:), ALLOCATABLE        :: orthonormality

!   ---------------------------------------------------------------------------

    NULLIFY(dft_control,scf_control,ks_env,qs_env)
    NULLIFY(rho, matrix_h, matrix_ks, matrix_p, matrix_s, kg_fm_set)
    NULLIFY(logger)

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

    CALL init_kgscf_run(scf_env=scf_env,kg_env=kg_env,&
                        globenv=globenv,error=error)

    CALL get_qs_env(qs_env=qs_env,&
                    energy=energy,&
                    ks_env=ks_env,&
                    matrix_h=matrix_h,&
                    matrix_ks=matrix_ks,&
                    matrix_s=matrix_s,&
                    qs_charges=qs_charges,&
                    rho=rho)

    CALL get_kg_env(kg_env=kg_env, kg_fm_set=kg_fm_set)
    nmolecule_kind = SIZE(kg_fm_set%kg_fm_mol_set,1)
    nelectron_global = kg_fm_set%nelectron_global
    nmolecule_global = kg_fm_set%nmolecule_global


    matrix_p => rho%rho_ao
 
    ! do some assertions here on these matrices having the same structure, 
    ! as is currently required
    DO ispin=1,SIZE(matrix_s)
      id_equal=(matrix_s(ispin)%matrix%sparsity_id==matrix_s(1)%matrix%sparsity_id)
      CPPrecondition(id_equal,cp_failure_level,routineP,error,failure)
    ENDDO
    DO ispin=1,SIZE(matrix_p)
      id_equal=(matrix_p(ispin)%matrix%sparsity_id==matrix_s(1)%matrix%sparsity_id)
      CPPrecondition(id_equal,cp_failure_level,routineP,error,failure)
    ENDDO
    DO ispin=1,SIZE(matrix_h)
      id_equal=(matrix_h(ispin)%matrix%sparsity_id==matrix_s(1)%matrix%sparsity_id)
      CPPrecondition(id_equal,cp_failure_level,routineP,error,failure)
    ENDDO
    DO ispin=1,SIZE(matrix_ks)
      id_equal=(matrix_ks(ispin)%matrix%sparsity_id==matrix_s(1)%matrix%sparsity_id)
      CPPrecondition(id_equal,cp_failure_level,routineP,error,failure)
    ENDDO
    ! end sparsity check


    scf_env%iter_count = 0
    diis_step = .FALSE.
    use_jacobi = .FALSE.

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
                                calculate_forces=.FALSE.,&
                                just_energy=energy_only)

       DO ispin=1,dft_control%nspins

          IF (logger%print_keys%density_matrix.AND.&
               logger%print_keys%each_scf_step) THEN
             CALL write_sparse_matrix(matrix_p(ispin)%matrix,4,6,qs_env,globenv)
          END IF

          IF (logger%print_keys%kohn_sham_matrix.AND.&
                         logger%print_keys%each_scf_step) THEN
             CALL write_sparse_matrix(matrix_ks(ispin)%matrix,4,6,qs_env,globenv)
          END IF

       ENDDO

       IF (logger%print_keys%mo_orthonormality) THEN
         ALLOCATE(orthonormality(nmolecule_global), STAT =istat)
         IF (istat /= 0) CALL stop_memory(routineP,"orthonormality",nmolecule_global*dp_size)
 
         orthonormality = 0.0_dp
         CALL fm_mol_orthonormality(orthonormality,kg_fm_set%kg_fm_mol_set,&
                                    matrix_s(1)%matrix,dft_control%nspins)
         CALL mp_sum(orthonormality,qs_env%para_env%group)

         IF (ionode) THEN
           DO imol = 1, nmolecule_global
            WRITE(output_unit,'(T8,A,T60,I4,E20.4)')  &
               " Maximum deviation from MOS-orthonormality for molecule",&
               imol,orthonormality(imol)
           END DO
         ENDIF
         DEALLOCATE(orthonormality, STAT =istat)
         IF (istat /= 0) CALL stop_memory(routineP,"deall. orthonormality")
       ENDIF

       scf_env%p_mix = scf_control%p_mix
       scf_env%iter_param = 0.0_dp
       IF (scf_env%mixing_method==ao_mixing_nr) scf_env%iter_param = scf_env%p_mix

 write(*,*) scf_env%mixing_method, ao_mixing_nr
 write(*,*)  scf_env%method, general_diag_method_nr,special_diag_method_nr,ot_method_nr
 STOP 'method'
       SELECT CASE (scf_env%method)
         CASE DEFAULT
            CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP,"unknown scf method method:"//&
               cp_to_string(scf_env%method),error,failure)
         CASE(general_diag_method_nr) ! diagonalisation (default)
  
           IF (scf_env%iter_count > 1) THEN
             ib = MODULO(scf_env%scf_diis_buffer%ncall,scf_env%scf_diis_buffer%nbuffer) + 1
             scf_env%scf_diis_buffer%ncall = scf_env%scf_diis_buffer%ncall + 1
             nb = MIN(scf_env%scf_diis_buffer%ncall, scf_env%scf_diis_buffer%nbuffer)
             nb1 = nb + 1
             ! Allocate the eigenvalue vector for the diis matrix of the global system
             ALLOCATE(diis_ev_global(nb1), STAT = istat)
             IF (istat /= 0) CALL stop_memory(routineP,"diis_ev_global",nb1*dp_size)

             ! Allocate the diis matrix for the global system
             ALLOCATE(diis_m1_global(nb1,nb1), STAT = istat)
             IF (istat /= 0) CALL stop_memory(routineP,"diis_m1_global",nb1*nb1*dp_size)
             ALLOCATE(diis_m2_global(nb1,nb1), STAT = istat)
             IF (istat /= 0) CALL stop_memory(routineP,"diis_m2_global",nb1*nb1*dp_size)
             ! Initialize the diis error
             diis_error = 0.0_dp

             CALL kg_diis_step_A(scf_env%scf_diis_buffer,ib,nb,&
                                 kg_fm_set%kg_fm_mol_set,&
                                 diis_m1_global,&
                                 dft_control%nspins,diis_error,&
                                 matrix_ks,&
                                 overlap=matrix_s(1)%matrix,&
                                 error=error)

             ! Global sum to collect all the contribution to the diis matrix
             CALL mp_sum(diis_m1_global,qs_env%para_env%group)
             ! Global max to get the maximum of the diis error over all the proc.
             CALL mp_max(diis_error,qs_env%para_env%group)

             !   *** Check, if a DIIS step is appropiate ***
             diis_step = ((scf_env%scf_diis_buffer%ncall > 1).AND.&
                         (scf_env%iter_delta < scf_control%eps_diis))

             IF(diis_error < scf_control%eps_diis) THEN

               diis_m1_global(1:nb,nb1) = -1.0_dp
               diis_m1_global(nb1,1:nb) = -1.0_dp
               diis_m1_global(nb1,nb1)  =  0.0_dp
               !     *** Solve the linear DIIS equation system ***
               CALL diamat_all(diis_m1_global(1:nb1,1:nb1),diis_ev_global(1:nb1))
               diis_m2_global(1:nb1,1:nb1) = diis_m1_global(1:nb1,1:nb1)
               DO jb=1,nb1
                 IF (ABS(diis_ev_global(jb)) < 1.0E-12_dp) THEN
                   diis_m2_global(1:nb1,jb) = 0.0_dp
                 ELSE
                   diis_m2_global(1:nb1,jb) = diis_m2_global(1:nb1,jb)/diis_ev_global(jb)
                 END IF
               END DO

               diis_ev_global(1:nb) = MATMUL(diis_m2_global(1:nb,1:nb1),&
                                             diis_m1_global(nb1,1:nb1))

             ELSE
               diis_step = .FALSE.
             END IF

             ! Allocate the diis matrix for the global system
             DEALLOCATE(diis_m1_global, STAT = istat)
             IF (istat /= 0) CALL stop_memory(routineP,"deall. diis_m1_global")
             DEALLOCATE(diis_m2_global, STAT = istat)
             IF (istat /= 0) CALL stop_memory(routineP,"deall. diis_m2_global")

           END IF  ! scf_env%iter_count > 1

           ! Control Variables
           do_level_shift = ((scf_control%level_shift /= 0.0_dp).AND.&
               ((scf_control%density_guess == "CORE").OR.(scf_env%iter_count > 1)))

           IF (diis_step) THEN
               scf_env%p_mix = 1.0_dp
               scf_env%iter_param = diis_error
               IF (use_jacobi) THEN
                  scf_env%iter_method = "DIIS/Jacobi"
               ELSE
                  scf_env%iter_method = "DIIS/Diag."
               END IF
           ELSE
               IF (use_jacobi) THEN
                  scf_env%iter_method = "Mixing/Jacobi"
               ELSE
                  scf_env%iter_method = "Mixing/Diag."
               END IF
           END IF
           scf_env%iter_delta=0.0_dp

           CALL kg_diis_step_B(scf_env%scf_diis_buffer,ib,nb,&
                               kg_fm_set%kg_fm_mol_set,&
                               scf_control,diis_ev_global,&
                               dft_control%nspins,&
                               matrix_ks,matrix_p,&
                               do_iter = (scf_env%iter_count > 1),&
                               diis_step = diis_step,&
                               do_level_shift=do_level_shift,&
                               use_jacobi=use_jacobi,&
                               error=error)

           IF(scf_env%iter_count > 1) THEN
             DEALLOCATE(diis_ev_global, STAT = istat)
             IF (istat /= 0) CALL stop_memory(routineP,"deall. diis_ev_global")
           END IF
         CASE(special_diag_method_nr)
            ! booooooooooooo
            CALL stop_program(routineN,moduleN,__LINE__,&
                     "Special_diag_method_nr not implemented with KG_GPW yet")
         CASE(ot_method_nr) ! orbital transforms
            ! booooooooooooo
            CALL stop_program(routineN,moduleN,__LINE__,&
                     "Ot_method_nr not implemented with KG_GPW yet")
       END SELECT

       IF (scf_env%mixing_method.EQ.ao_mixing_nr) THEN
          CALL scf_env_density_mixing(scf_env%p_mix_new,scf_env%mixing_method,&
                                    scf_env%p_mix,scf_env%iter_delta,&
                                    qs_env=qs_env,error=error)
       ENDIF

       t2 = m_walltime()

       IF (ionode.AND.logger%print_keys%scf.and.scf_env%print_iter_line) THEN
          WRITE (UNIT=output_unit,&
               FMT="(T2,I5,2X,A,T22,E10.2,T32,F8.2,T40,2F20.10)")&
               scf_env%iter_count,TRIM(scf_env%iter_method),&
               scf_env%iter_param,t2 - t1,scf_env%iter_delta,energy%total
       END IF

! ** convergece check
       IF (qs_env%para_env%source==qs_env%para_env%mepos) &
            INQUIRE(file="EXIT",exist=should_stop)
       CALL mp_bcast(should_stop,qs_env%para_env%source,qs_env%para_env%group)
       IF (scf_env%iter_delta < scf_control%eps_scf) THEN
          IF (ionode.AND.logger%print_keys%scf) THEN
             WRITE(UNIT=output_unit,FMT="(/,T3,A,I5,A/)")&
                  "*** SCF run converged in ",scf_env%iter_count," steps ***"
          END IF
          CALL timestop(0.0_dp,handle2)
          EXIT scf_loop
       ELSE IF (should_stop.OR.&
            scf_env%iter_count == scf_control%max_scf) THEN
          IF (ionode.AND.logger%print_keys%scf) THEN
             WRITE(UNIT=output_unit,FMT="(/,T3,A,/)")&
                  "*** SCF run NOT converged ***"
          END IF
          CALL timestop(0.0_dp,handle2)
          EXIT scf_loop
       END IF

!!   *** Write restart file ***
!       IF (scf_control%write_restart_each>0.AND.&
!            MODULO(scf_env%iter_count,scf_control%write_restart_each)==0) THEN
!         ! CALL write_mo_set(c,atomic_kind_set,particle_set,globenv,&
!         !      id_nr=qs_env%id_nr)
!          CALL write_mo_set(c,atomic_kind_set,particle_set,globenv,&
!               id_nr=1)
!       END IF
!
       IF (.NOT. logger%print_keys%timecumul) t1 = m_walltime()
!
!!   *** mixing methods have the new density matrix in p_mix_new
!       IF (scf_env%mixing_method.EQ.ao_mixing_nr) THEN
!           DO ispin=1,dft_control%nspins
!              CALL transfer_matrix(scf_env%p_mix_new(ispin)%matrix,&
!                                       rho%rho_ao(ispin)%matrix)
!           END DO
!       ENDIF
!
! ** update qs_env%rho
       CALL qs_rho_update_rho(rho, qs_env=qs_env, error=error)
       CALL qs_ks_did_change(ks_env,rho_changed=.TRUE.,error=error)      

       CALL timestop(0.0_dp,handle2)

    END DO scf_loop

    IF (ionode.AND.logger%print_keys%scf) THEN
       WRITE (UNIT=output_unit,FMT="(/,(T3,A,T40,2F20.10))")&
            "Total electronic density (r-space): ",&
            SUM(rho%tot_rho_r),&
            SUM(rho%tot_rho_r)+ REAL(nelectron_global,dp),&
            "Total core charge density (r-space):",&
            qs_charges%total_rho_core_rspace,&
            qs_charges%total_rho_core_rspace - REAL(nelectron_global+dft_control%charge,dp)
       WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
            "Total charge density (r-space):     ",&
            SUM(rho%tot_rho_r)+&
            qs_charges%total_rho_core_rspace,&
            "Total charge density (g-space):     ",qs_charges%total_rho_gspace
       WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
            "Overlap energy of the core charge distribution:",energy%core_overlap,&
            "Self energy of the core charge distribution:   ",energy%core_self,&
            "Core Hamiltonian energy:                       ",energy%core,&
            "Hartree energy:                                ",energy%hartree,&
            "Exchange-correlation energy:                   ",energy%exc,&
            "Total energy:                                  ",energy%total
       CALL m_flush(output_unit)
    END IF

!   *** Write restart file ***
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
!       CALL scf_env_density_mixing(scf_env%p_mix_new,scf_env%mixing_method,&
!                                    scf_env%p_mix,scf_env%iter_delta,&
!                                    qs_env=qs_env,invert=.TRUE.,error=error)
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
!   *** update rspace rho since the mo changed
!   *** this might not always be needed (i.e. no post calculation / no forces )
    CALL qs_rho_update_rho(rho, qs_env=qs_env, error=error)

    CALL qs_ks_did_change(ks_env,rho_changed=.TRUE.,error=error)

! *** print out info
!    IF (.NOT.should_stop) THEN
!       CALL scf_post_calculation(scf_env=scf_env, qs_env=qs_env, &
!            globenv=globenv, error=error)
!    END IF

! *** cleanup
    CALL kg_scf_env_cleanup(scf_env,qs_env=qs_env,error=error)

    CALL timestop(0.0_dp,handle)

    CALL write_checkpoint_information("leaving "//routineN,globenv)

  END SUBROUTINE scf_env_do_scf
!**************************************************************************

  SUBROUTINE init_kgscf_run(scf_env,kg_env,globenv,error)

!   Purpose : Initialize the kg_env%scf_env 
!              - allocate the ks_matrix and the ks environment
!                which are contained in kg_env%sub_qs_env
!              - allocate p_mix_new  if the ao_mixing is used (in kg_env%scf_env)
!              - calculate the self and overlap energies (in kg_env%sub_qs_env) 
!              - from the overlap matrix extract the molecular blocks and
!                construct the corresponding ortho matrix
!              - Initialize the separated blocks of the density matrix
!                calculated for each molecule separatedly, and
!                copy them in the corresponding blocks of the sparse matrix matrix_p

! ***************************************************************************

    TYPE(kg_scf_env_type), POINTER           :: scf_env
    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(global_environment_type), &
      INTENT(in)                             :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "init_kgscf_run", &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ndep, output_unit
!    INTEGER                                  :: homo, ispin, nao, &
!                                                ndep, nelectron_spin, &
    LOGICAL                                  :: failure
    LOGICAL                                  :: do_ortho,ionode

    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_sm_pool_p_type), DIMENSION(:), &
      POINTER                                :: S_sm_pools
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(kg_fm_p_type), POINTER              :: kg_fm
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(scf_control_type), POINTER          :: scf_control

!    TYPE(atomic_kind_type), DIMENSION(:), &
!      POINTER                                :: atomic_kind_set
!   ---------------------------------------------------------------------------

    CALL timeset("init_kgscf_run","I"," ",handle)

    NULLIFY(dft_control,kg_fm,scf_control)
    NULLIFY(qs_env,ks_env,matrix_ks)

    failure=.FALSE.
    do_ortho = .TRUE.
    logger => cp_error_get_logger(error)
    
    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(kg_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(kg_env%ref_count>0,cp_failure_level,routineP,error,failure)
    
    CALL get_kg_env(kg_env=kg_env,&
                    dft_control=dft_control,&
                    kg_fm_set=kg_fm)

    IF ( dft_control % qs_control % method == "KG_GPW" ) THEN

      CALL get_kg_env(kg_env=kg_env,&
                      sub_qs_env=qs_env)

      CALL get_qs_env(qs_env=qs_env,&
           dft_control=dft_control,&
           scf_control=scf_control,&
           matrix_ks=matrix_ks,&
           matrix_s=matrix_s,&
           ks_env=ks_env,&
           error=error)

      CALL mpools_get(qs_env%mpools, S_sm_pools=S_sm_pools,&
                       error=error)

!   *** Allocate matrix_ks and put it in the QS environment ***
      IF (.not.ASSOCIATED(matrix_ks)) THEN
         CALL sm_pools_create_matrix_vect(S_sm_pools,matrix_ks,&
              name="SCF"//TRIM(ADJUSTL(cp_to_string(scf_env%id_nr)))//&
              "KOHN-SHAM_MATRIX",&
              error=error)
         CALL set_qs_env(qs_env=qs_env,&
              matrix_ks=matrix_ks)
      END IF

!   *** allocate p_mix_new ***
      IF (scf_env%mixing_method.EQ.ao_mixing_nr) THEN
         IF (.not.ASSOCIATED(scf_env%p_mix_new)) THEN
            CALL sm_pools_create_matrix_vect(S_sm_pools,scf_env%p_mix_new,&
                 name="SCF"//TRIM(ADJUSTL(cp_to_string(scf_env%id_nr)))//&
                 "DENSITY",&
                 error=error)
         END IF
      END IF

!   *** allocate the ks env **
      IF (.not.ASSOCIATED(ks_env)) THEN
         CALL qs_ks_create(ks_env,qs_env=qs_env)
         CALL set_qs_env(qs_env, ks_env=ks_env)
         CALL qs_ks_release(ks_env)
      END IF

!         atomic_kind_set=atomic_kind_set,&
!         c=c,k=k,&
!         h=h,s=s, &
!    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
!         gth_potential_present=gth_potential_present)
!
      ionode = logger%para_env%source == logger%para_env%mepos
      IF (ionode) output_unit = cp_logger_get_default_unit_nr(logger)
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

!! calc ortho matrix
      ndep = 0
      do_ortho = (scf_env%method /= ot_method_nr)
      CALL calculate_ortho_per_molecule(kg_fm%kg_fm_mol_set,matrix_s(1)%matrix,&
                                        ndep,scf_control,do_ortho,error=error)
!          IF (ionode.and.logger%print_keys%scf) THEN
!             WRITE (UNIT=output_unit,FMT="(T3,A,I10)")&
!                  "Number of independent orbital functions:",nao - ndep
!          END IF
!    IF (logger%print_keys%ortho_matrix) THEN
!       CALL write_fm_with_basis_info(ortho,4,6,qs_env,globenv)
!    END IF
!
!   Initializes rho and the mosv
    CALL kgscf_initial_rho(kg_fm%kg_fm_mol_set,qs_env=qs_env,&
                           globenv=globenv,error=error)

! *** method dependent initializations ***
    SELECT CASE (scf_env%method)
      CASE DEFAULT
        CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                 routineP,"unknown scf method method:"//&
                 cp_to_string(scf_env%method),error,failure)
      CASE (general_diag_method_nr)
         IF (.NOT.ASSOCIATED(scf_env%scf_diis_buffer)) THEN       
           CALL kg_diis_b_create(scf_env%scf_diis_buffer,&
                                 kg_fm%kg_fm_mol_set,&
                                 nbuffer=scf_control%max_diis,&
                                 nspins = dft_control%nspins,&
                                 error=error)
         END IF
         CALL kg_diis_b_clear(scf_env%scf_diis_buffer,error=error)

      CASE (special_diag_method_nr)
          CALL stop_program(routineN,moduleN,__LINE__,&
               "Special_diag_method_nr not implemented with KG_GPW yet")
      CASE (ot_method_nr)
          CALL stop_program(routineN,moduleN,__LINE__,&
               "Ot_method_nr not implemented with KG_GPW yet")

    END SELECT

   ELSE

     CPErrorMessage(-101,routineP,"Illegal Method",error)

   END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE init_kgscf_run

!***************************************************************************

  SUBROUTINE calculate_ortho_per_molecule(kg_fm_mol_set,s,ndep,scf_control,&
                                          do_orthogonalization,error)
!
!   Purpose: Calculate the orthogonalization matrix for each molecule
!            separatedly.
!            For a given molecule the corresponding portion of the overlap 
!            matrix is copied in the (n_ao x n_ao) ortho matrix , 
!            where n_ao is restricted to the number of  atomic orbitals 
!            centered on the atoms of the molecule. This matrix is then passed 
!            to the cholesky or the canonical orthogonalization routines. 
!            Each molecule, and cosequently each ortho matrix, is assigned 
!            to one single processor, therefore the orthogonalization operations 
!            are not distributed further. 
!            (I suppose omp can add further parallelization) 
!
!   History: - Creation (24.11.2004, MI)
!
!   ***************************************************************************

    TYPE(kg_fm_mol_set_type), DIMENSION(:), &
      POINTER                                :: kg_fm_mol_set
    TYPE(real_matrix_type), POINTER          :: s
    TYPE(scf_control_type), POINTER          :: scf_control
    INTEGER, INTENT(OUT)                     :: ndep
    LOGICAL,INTENT(IN)                       :: do_orthogonalization
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "calculate_ortho_per_molecule", &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: col_first, col_last, &
                                                handle, iat, iatom, icol,&
                                                icol_ortho, imol, imol_kind, &
                                                info, irow, irow_ortho, jat, &
                                                jatom, ld_ortho, ld_work, natom, &
                                                nmol, nmol_kind, nc_ortho, nc_work,&
                                                row_first, row_last
    REAL(KIND=dp)                            :: eps_eigval
    TYPE(fm_mol_blocks_type), DIMENSION(:),&
      POINTER                                :: fm_mol_blocks
    TYPE(fm_mol_blocks_type), POINTER        :: mol_block
    TYPE(kg_fm_mol_set_type), POINTER        :: fm_mol_set 

    REAL(dp), DIMENSION(:,:),  POINTER        :: ortho
    REAL(dp), DIMENSION(:,:),  POINTER        :: work
    REAL(dp), DIMENSION(:,:),  POINTER        :: s_block
    INTEGER, DIMENSION(:), POINTER            :: i_atom
    INTEGER, DIMENSION(:), POINTER            :: ifirst_ao
    INTEGER, DIMENSION(:), POINTER            :: ilast_ao

    LOGICAL :: failure
    INTEGER :: i,j
!
!   ---------------------------------------------------------------------------
!
    CALL timeset("calculate_ortho_per_molecule","I"," ",handle)

    failure = .FALSE.
    nmol_kind = SIZE(kg_fm_mol_set,1)

    DO imol_kind = 1,nmol_kind
      NULLIFY(fm_mol_set,fm_mol_blocks)
      fm_mol_set => kg_fm_mol_set(imol_kind)
      CALL get_kg_fm_mol_set(kg_fm_mol_set=fm_mol_set,&
                             nmolecule_local=nmol,&
                             natom=natom,&
                             fm_mol_blocks=fm_mol_blocks)

      DO imol = 1,nmol
        NULLIFY(mol_block)
        mol_block => fm_mol_blocks(imol)
        NULLIFY(ortho,work,i_atom,ifirst_ao,ilast_ao)
        CALL get_fm_mol_block(fm_mol_block = mol_block,&
                              index_atom = i_atom,&
                              ortho = ortho,&
                              work = work,&
                              ifirst_ao = ifirst_ao,&
                              ilast_ao = ilast_ao )

        ld_ortho = SIZE(ortho,1)
        nc_ortho = SIZE(ortho,2)
        ld_work = SIZE(work,1)
        nc_work = SIZE(work,2)
        ! Construct ortho from the overlap matrix
        DO iat = 1,natom
           iatom = i_atom(iat)
           row_first = ifirst_ao(iat)
           row_last  = ilast_ao(iat)
           DO jat = iat , natom
             jatom = i_atom(jat)
             col_first = ifirst_ao(jat)
             col_last  = ilast_ao(jat)

             NULLIFY(s_block)
             IF (iatom <= jatom) THEN
               CALL get_block_node(matrix=s,&
                                   block_row=iatom,&
                                   block_col=jatom,&
                                   BLOCK=s_block)
               icol =  1
               DO icol_ortho = col_first, col_last
                  irow = 1
                  DO irow_ortho = row_first,   row_last
                     ortho(irow_ortho,icol_ortho) = s_block(irow,icol)
                     irow = irow + 1
                  END DO  ! irow_ortho 
                  icol = icol + 1
               END DO   ! icol_ortho
             ELSE
               CALL get_block_node(matrix=s,&
                                   block_row=jatom,&
                                   block_col=iatom,&
                                   BLOCK=s_block)

               icol =  1
               DO icol_ortho = col_first, col_last
                  irow = 1
                  DO irow_ortho = row_first,   row_last
                     ortho(irow_ortho,icol_ortho) = s_block(icol,irow)
                     irow = irow + 1
                  END DO  ! irow_ortho 
                  icol = icol + 1
               END DO   ! icol_ortho
             END IF
           END DO  ! jat
        END DO  ! iat

!  write(*,*) 'molecule ' ,imol
!  do i = 1,ld_ortho
!    write(*,'(I4,23f10.5)') i,ortho(i,1:23) 
!  end do

        IF (do_orthogonalization) THEN
          IF (scf_control%use_cholesky ) THEN

            CALL dpotrf('U',ld_ortho,ortho(1,1),ld_ortho,info)
            CPPostcondition(info==0,cp_failure_level,routineP,error,failure)

          ELSE

            CALL fm_mol_power(ortho,work,ld_work,nc_work,-0.5_dp,scf_control%eps_eigval,ndep)
            ! Transform the upper triangular matrix in a full matrix
            DO irow = 1,ld_ortho 
              DO icol = irow+1,nc_ortho
                ortho(icol,irow) = ortho(irow,icol)
              END DO
            END DO
!           IF (ionode.and.logger%print_keys%scf) THEN
!             WRITE (UNIT=output_unit,FMT="(T3,A,I10)")&
!                  "Number of independent orbital functions:",nao - ndep
!           END IF
          ENDIF
        ENDIF

!  write(*,*) 'after ort '
!  do i = 1,ld_ortho
!    write(*,'(I4,23f10.5)') i,ortho(i,1:23) 
!  end do
      END DO  ! imol
    END DO  ! imol_kind

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE calculate_ortho_per_molecule
!
! *****************************************************************************
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
  SUBROUTINE kgscf_first_density_matrix(kg_fm_mol_set,qs_env,globenv,error)

    TYPE(kg_fm_mol_set_type), DIMENSION(:), &
      POINTER                                :: kg_fm_mol_set
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "kgscf_first_density_matrix", &
      routineP = moduleN//":"//routineN

    INTEGER :: handle, imol, imolecule_kind, ispin, istat, n_ao_kind, nat_mol,&
               nmo, nmo_kind(2), nmo_max, nmol_local, nmolecule_kind, nspin ,&
               qs_env_id , i

!    INTEGER :: atom_a, group, handle, homo, iatom, ikind, iset, isgf, isgfa, &
!      ishell, ispin, istat, la, maxl, nao, natom, ncount, nelectron, nmo, &
!      nset, nspin, 
!    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_sgf
!    INTEGER, DIMENSION(:), POINTER           :: atom_list, elec_conf, nshell
!    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, l, last_sgfa
!    REAL(KIND=dp)                            :: maxocc, nelec, paa, scale, &
!                                                trps1, trps2
!    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: econf, pdiag, sdiag
!    TYPE(all_potential_type), POINTER        :: all_potential
!    TYPE(atomic_kind_type), POINTER          :: atomic_kind
!    TYPE(gth_potential_type), POINTER        :: gth_potential
!    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    REAL(dp), DIMENSION(:,:), POINTER         :: mo_coeff
    REAL(dp), DIMENSION(:,:), POINTER         :: aux_matrix
    REAL(dp), DIMENSION(:,:), POINTER         :: ortho
    REAL(dp), DIMENSION(:,:), POINTER         :: work
    REAL(dp), DIMENSION(:,:), POINTER         :: s_mo
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(fm_mol_blocks_type), DIMENSION(:), &
      POINTER                                :: fm_mol_blocks
    TYPE(fm_mol_blocks_type), POINTER        :: mol_block
    TYPE(kg_fm_mol_set_type), POINTER        :: fm_mol_set
    TYPE(mol_mo_set_p_type), DIMENSION(:),&
      POINTER                                :: mos
    TYPE(mol_mo_set_type), POINTER           :: mo_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_h,&
                                                matrix_p,&
                                                matrix_s
    TYPE(scf_control_type), POINTER          :: scf_control
!
!   ---------------------------------------------------------------------------
!
    NULLIFY( aux_matrix, atomic_kind_set, dft_control, &
            para_env, particle_set, &
            matrix_h, matrix_s,matrix_p, scf_control)

    CALL timeset(routineN,"I"," ",handle)

    CALL get_qs_env(qs_env,atomic_kind_set=atomic_kind_set,&
                    para_env=para_env, particle_set=particle_set, &
                    matrix_h=matrix_h,&
                    matrix_s=matrix_s,&
                    scf_control=scf_control, id_nr=qs_env_id,&
                    dft_control=dft_control)

    nspin=dft_control%nspins
    nmolecule_kind = SIZE(kg_fm_mol_set)

    matrix_p => qs_env%rho%rho_ao

    IF (scf_control%use_ot .AND. .NOT.  &
         (scf_control%density_guess == "RANDOM" .OR.  &
         (scf_control%density_guess == "RESTART" .AND. &
         scf_control%level_shift == 0.0_dp))) THEN
       CALL stop_program("calculate_first_density_matrix",&
                         "OT needs GUESS RESTART : other options NYI")
    ENDIF

    IF (scf_control%density_guess == "RESTART") THEN
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
    ELSE IF (scf_control%density_guess == "RANDOM") THEN

      DO imolecule_kind = 1,nmolecule_kind
        NULLIFY(fm_mol_set,fm_mol_blocks)
        fm_mol_set => kg_fm_mol_set(imolecule_kind)
        CALL get_kg_fm_mol_set(kg_fm_mol_set=fm_mol_set, &
                              natom = nat_mol, &
                              n_ao = n_ao_kind , n_mo = nmo_kind, &
                              fm_mol_blocks = fm_mol_blocks,&
                              nmolecule_local = nmol_local)

        nmo_max = MAX(nmo_kind(1),nmo_kind(2))
        ALLOCATE(s_mo(n_ao_kind,nmo_max), STAT = istat)
        IF (istat /= 0) CALL stop_memory(routineP,"s_mo",n_ao_kind*nmo_max*dp_size)

        DO imol = 1, nmol_local
          NULLIFY(mol_block,mos)
          mol_block => fm_mol_blocks(imol)
          CALL get_fm_mol_block(fm_mol_block = mol_block,&
                                mos = mos)
          DO ispin=1,nspin
            NULLIFY(mo_set,mo_coeff)
            mo_set => mos(ispin)%mo_set
            CALL get_mol_mo_set(mo_set=mo_set, mo = mo_coeff,&
                                nmo=nmo)
            ! Initialize random coefficients
            CALL  mol_mo_random(mo_coeff,para_env%mepos,para_env%num_pe) 
            ! Multiply Overlap matric and Coefficient matrix
            CALL  multiply_sparse_mol_mo(matrix_s(ispin)%matrix,mol_block,&
                                         nat_mol,mo_coeff,nmo,s_mo)
            ! orthogonalize the molecular orbitals 
            CALL mol_make_basis(mo_coeff, nmo, matrix_ortho=s_mo, otype="SV")

            CALL calculate_mol_density_matrix(mo_set,mol_block,&
                                       nat_mol,matrix_p(ispin)%matrix,error=error)
          ENDDO  ! ispin
        END DO  ! imol

        DEALLOCATE(s_mo,STAT = istat)
        IF (istat /= 0) CALL stop_memory(routineP,"s_mo")

      END DO  ! imolecule_kind

    ELSE IF (scf_control%density_guess == "CORE") THEN

      DO imolecule_kind = 1,nmolecule_kind
        NULLIFY(fm_mol_set,fm_mol_blocks)
        fm_mol_set => kg_fm_mol_set(imolecule_kind)
        CALL get_kg_fm_mol_set(kg_fm_mol_set=fm_mol_set, &
                              natom = nat_mol, &
                              n_ao = n_ao_kind , n_mo = nmo_kind, &
                              fm_mol_blocks = fm_mol_blocks,&
                              nmolecule_local = nmol_local)

        ALLOCATE(aux_matrix(n_ao_kind,n_ao_kind), STAT = istat)
        IF (istat /= 0) CALL stop_memory(routineP,"aux_matrix",&
                                         n_ao_kind*n_ao_kind*dp_size)

        DO imol = 1, nmol_local
          NULLIFY(mol_block,mos,ortho,work)
          mol_block => fm_mol_blocks(imol)
          CALL get_fm_mol_block(fm_mol_block = mol_block,&
                                ortho = ortho,&
                                work = work,&
                                mos = mos)

          ! Copy the core hamiltonian of this molecule in the aux_matrix
          CALL copy_sparse2mol_block(matrix_h(1)%matrix, mol_block,&
                               aux_matrix, nat_mol, n_ao_kind, n_ao_kind, error)
  
! dbg
!     write(*,*) 'matrix_h',n_ao_kind,nat_mol
!     do i =1,n_ao_kind
!       write(*,'(I4,23f10.5)') i,aux_matrix(i,1:n_ao_kind) 
!     end do
! stop 'h_core'
! dbg
          DO ispin=1,nspin
            NULLIFY(mo_set,mo_coeff)
            mo_set => mos(ispin)%mo_set

            ! Diagonalize the core hamiltonian and copy the eigenvectors in mo_set 
            CALL fm_mol_eigensolver(aux_matrix,mo_set,ortho,work,n_ao_kind,&
                                   .FALSE.,0.0_dp,&
                                    use_cholesky=scf_control%use_cholesky, &
                                    use_jacobi=.FALSE.,&
                                    jacobi_threshold=scf_control%jacobi_threshold,&
                                    smear=0.0_dp, error=error)

            IF (scf_control%level_shift /= 0.0_dp) THEN
              CALL get_mol_mo_set(mo_set=mo_set,&
                                  mo=mo_coeff)
              CALL DCOPY(n_ao_kind*n_ao_kind,mo_coeff(1,1),1,ortho(1,1),1)
            END IF

            CALL calculate_mol_density_matrix(mo_set,mol_block,&
                             nat_mol,matrix_p(ispin)%matrix,error=error)
          ENDDO  ! ispin
        END DO  ! imol

        DEALLOCATE(aux_matrix,STAT = istat)
        IF (istat /= 0) CALL stop_memory(routineP,"aux_matrix")

      END DO  ! imolecule_kind

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
    ELSE

      CALL stop_program(routineP,&
                        "An invalid keyword for the initial density "//&
                        "guess was specified")

    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE kgscf_first_density_matrix

!!****f* kg_scf/kgscf_initial_rho [1.0] *
!!
!!   NAME
!!     kgscf_initial_rho
!!
!!   FUNCTION
!!     Initializes rho and the mos for the kg_gpw method. 
!!     Each molecule contributes a set of mos. Each set
!!     is generated (or read), orthogonalized and opportunely
!!     transformed into the corresponding portion of the density matrix
!!     independently from the others.
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - kg_fm_mol_set: the set of mos and ortho matrixes, one per each molecule 
!!     - qs_env: the qs env that lives in kg_env
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!      MI
!!
!!   MODIFICATION HISTORY
!!     11.2004 created 
!!
!!*************************************************************************

  SUBROUTINE kgscf_initial_rho(kg_fm_mol_set, qs_env, globenv, error)

    TYPE(kg_fm_mol_set_type), DIMENSION(:), &
      POINTER                                :: kg_fm_mol_set
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), &
      INTENT(in), TARGET                     :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kgscf_initial_rho_setup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: extrapolation_method_nr, &
                                                handle, ispin, nmo
    LOGICAL                                  :: failure, ionode, orthogonal_wf
    TYPE(cp_logger_type), POINTER            :: logger

!    -------------------------------------------------------------------------

    CALL timeset(routineN,"I"," ",handle)

    failure=.FALSE.

    logger => cp_error_get_logger(error)
    ionode = logger%para_env%source==logger%para_env%mepos
    CPPrecondition(ASSOCIATED(kg_fm_mol_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN
       extrapolation_method_nr=wfi_use_guess_method_nr
!       IF (ASSOCIATED(qs_env%wf_history)) THEN
!          CALL wfi_extrapolate(globenv,qs_env%wf_history, &
!               qs_env=qs_env, dt=1.0_dp, &
!               extrapolation_method_nr=extrapolation_method_nr,&
!               orthogonal_wf=orthogonal_wf, error=error)
!          ! wfi_use_guess_method_nr the wavefunctions are not yet initialized
!          IF (.NOT.orthogonal_wf.AND.scf_env%method==ot_method_nr.AND. &
!                          .NOT.extrapolation_method_nr==wfi_use_guess_method_nr) THEN
!             DO ispin=1,SIZE(qs_env%c)
!                CALL get_mo_set(qs_env%c(ispin)%mo_set, &
!                     mo_coeff=mo_coeff, nmo=nmo)
!                CALL qs_env_reorthogonalize_vectors(qs_env, &
!                     v_matrix=mo_coeff, ispin=ispin, n_col=nmo,&
!                     error=error)
!             END DO
!          END IF
!       END IF
       IF (ionode) &
       WRITE (UNIT=cp_logger_get_default_unit_nr(logger),FMT="(/,T3,A)")&
         "Extrapolation method: "//&
         TRIM(wfi_get_method_label(extrapolation_method_nr,error=error))
       IF (extrapolation_method_nr==wfi_use_guess_method_nr) THEN
          CALL kgscf_first_density_matrix(kg_fm_mol_set,qs_env=qs_env,&
               globenv=globenv, error=error)
          CALL qs_rho_update_rho(qs_env%rho,qs_env=qs_env, error=error)
          CALL qs_ks_did_change(qs_env%ks_env,rho_changed=.TRUE.,&
               error=error)
       END IF
    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE kgscf_initial_rho
!!***************************************************************************
!
!!****f* kg_scf/kg_scf_env_cleanup [1.0] *
!!
!!   NAME
!!     kg_scf_env_cleanup
!!
!!   FUNCTION
!!     perform cleanup operations (like releasing temporary storage)
!!     at the end of the kg_scf
!!
!!   ARGUMENTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   MODIFICATION HISTORY
!!     11.2004 created 
!!
!!*** **********************************************************************
  SUBROUTINE kg_scf_env_cleanup(scf_env,qs_env,error)

    TYPE(kg_scf_env_type), POINTER           :: scf_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
!
    CHARACTER(len=*), PARAMETER :: routineN = 'kg_scf_env_cleanup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, stat
    LOGICAL                                  :: failure
    TYPE(cp_sm_pool_p_type), DIMENSION(:), &
      POINTER                                :: S_sm_pools

!    -------------------------------------------------------------------------

    CALL timeset(routineN,"I"," ",handle)

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN
       CALL mpools_get(qs_env%mpools,S_sm_pools=S_sm_pools,error=error)

!!   *** Release SCF work storage ***
!
       IF (ASSOCIATED(scf_env%p_mix_new)) THEN
          CALL sm_pools_give_back_matrix_vect(S_sm_pools,scf_env%p_mix_new,&
               error=error)
          CALL sm_pools_flush_cache(S_sm_pools,error=error)
       END IF

! *** method dependent cleanup
       SELECT CASE(scf_env%method)
       CASE(ot_method_nr)
          CALL destroy_preconditioner(scf_env%ot_preconditioner)
          DEALLOCATE(scf_env%ot_preconditioner,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          DO ispin=1,SIZE(scf_env%qs_ot_env)
             CALL ot_scf_destroy(scf_env%qs_ot_env(ispin))
          ENDDO
          DEALLOCATE(scf_env%qs_ot_env,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       CASE(general_diag_method_nr)
          !
       CASE(special_diag_method_nr)
          !
       CASE default
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP,"unknown scf method method:"//&
               cp_to_string(scf_env%method),error,failure)
       END SELECT

    END IF
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE kg_scf_env_cleanup
!***************************************************************************
!
!!****f* kg_scf/kg_scf_env_did_change [1.0] *
!!
!!   NAME
!!     kg_scf_env_did_change
!!
!!   SYNOPSIS
!!     Subroutine kg_scf_env_did_change(scf_env, error)
!!       Type(kg_scf_env_type), Pointer:: scf_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine scf_env_did_change
!!
!   FUNCTION
!!     function to be called to inform the scf_env about changes
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - scf_env: the scf env to inform
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     11.2004 created [MI]
!!
!!*** **********************************************************************
SUBROUTINE kg_scf_env_did_change(scf_env,error)
    TYPE(kg_scf_env_type), POINTER           :: scf_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_scf_env_did_change', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

!    -------------------------------------------------------------------------

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     IF (ASSOCIATED(scf_env%p_mix_new)) THEN
        CALL deallocate_matrix_set(scf_env%p_mix_new)
     END IF
  END IF
END SUBROUTINE kg_scf_env_did_change

!***************************************************************************

END MODULE kg_scf
