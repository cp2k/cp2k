!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2005  CP2K developers group
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****** cp2k/qs_scf [1.0] *
!!
!!   NAME
!!     qs_scf
!!
!!   FUNCTION
!!     Routines for the Quickstep SCF run.
!!
!!   AUTHOR
!!     Matthias Krack (30.04.2001)
!!
!!   MODIFICATION HISTORY
!!     - Joost VandeVondele (02.2002)
!!          added code for: incremental (pab and gvg) update
!!                           initialisation (init_cube, l_info)
!!     - Joost VandeVondele (02.2002)
!!          called the poisson code of the classical part
!!          this takes into account the spherical cutoff and allows for
!!          isolated systems
!!     - Joost VandeVondele (02.2002)
!!          added multiple grid feature
!!          changed to spherical cutoff consistently (?)
!!          therefore removed the gradient correct functionals
!!     - updated with the new QS data structures (10.04.02,MK)
!!     - copy_matrix replaced by transfer_matrix (11.04.02,MK)
!!     - nrebuild_rho and nrebuild_gvg unified (12.04.02,MK)
!!     - set_mo_occupation for smearing of the MO occupation numbers
!!       (17.04.02,MK)
!!     - MO level shifting added (22.04.02,MK)
!!     - Usage of TYPE mo_set_p_type
!!     - Joost VandeVondele (05.2002)
!!           added cholesky based diagonalisation
!!     - 05.2002 added pao method [fawzi]
!!     - parallel FFT (JGH 22.05.2002)
!!     - 06.2002 moved KS matrix construction to qs_build_KS_matrix.F [fawzi]
!!     - started to include more LSD (01.2003,Joost VandeVondele)
!!     - 02.2003 scf_env [fawzi]
!!     - got rid of nrebuild (01.2004, Joost VandeVondele)
!!     - 10.2004 removed pao [fawzi]
!!
!!   SOURCE
!******************************************************************************

MODULE qs_scf

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_files,                        ONLY: close_file,open_file
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm,&
                                             cp_fm_symm,&
                                             cp_fm_upper_to_full
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                             cp_fm_cholesky_reduce,&
                                             cp_fm_cholesky_restore
  USE cp_fm_diag,                      ONLY: cp_fm_block_jacobi_classic,&
                                             cp_fm_power,&
                                             cp_fm_syevd,&
                                             cp_fm_syevx
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                             fm_pool_get_el_struct
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_get,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_element,&
                                             cp_fm_init_random,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_set_element,&
                                             cp_fm_to_fm,&
                                             cp_fm_type,&
                                             cp_fm_write
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: copy_sm_to_fm,&
                                             cp_sm_fm_multiply
  USE cp_sm_pool_types,                ONLY: cp_sm_pool_p_type,&
                                             sm_pools_create_matrix_vect,&
                                             sm_pools_flush_cache,&
                                             sm_pools_give_back_matrix_vect
  USE external_potential_types,        ONLY: all_potential_type,&
                                             get_potential,&
                                             gth_potential_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             int_size
  USE input_constants,                 ONLY: ot_precond_none,&
                                             ot_precond_full_all,&
                                             xas_tp_hh, &
                                             xas_tp_fh, &
                                             xas_dscf
  USE machine,                         ONLY: m_flush,&
                                             m_walltime
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_max,&
                                             mp_sum
  USE particle_types,                  ONLY: get_particle_set,&
                                             particle_type
  USE preconditioner,                  ONLY: destroy_preconditioner,&
                                             init_preconditioner,&
                                             make_preconditioner
  USE pw_env_types,                    ONLY: pw_env_type
  USE physcon,                         ONLY: kcalmol,evolt
  USE qs_charges_types,                ONLY: qs_charges_type
  USE qs_collocate_density,            ONLY: collocate_atomic_charge_density
  USE qs_core_energies,                ONLY: calculate_ecore_overlap,&
                                             calculate_ecore_self
  USE qs_diis,                         ONLY: qs_diis_b_clear,&
                                             qs_diis_b_create,&
                                             qs_diis_b_step
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_env_reorthogonalize_vectors,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_ks_methods,                   ONLY: qs_ks_create,&
                                             qs_ks_did_change,&
                                             qs_ks_update_qs_env,&
                                             qs_ks_build_kohn_sham_matrix
  USE qs_ks_types,                     ONLY: qs_ks_env_type,&
                                             qs_ks_release
  USE qs_loc_methods,                  ONLY: loc_initialize_rho
  USE qs_matrix_pools,                 ONLY: mpools_get
  USE qs_mo_methods,                   ONLY: calculate_density_matrix,&
                                             calculate_orthonormality,&
                                             make_basis,&
                                             calculate_subspace_eigenvalues
  USE qs_mo_types,                     ONLY: correct_mo_eigenvalues,&
                                             get_mo_set,& 
                                             init_mo_set,&
                                             mo_set_p_type,&
                                             mo_set_type,&
                                             read_mo_set,&
                                             mo_set_restrict,&
                                             set_mo_occupation,&
                                             set_mo_set,&
                                             write_mo_set
  USE qs_ot,                           ONLY: qs_ot_new_preconditioner, &
                                             qs_ot_type
  USE qs_ot_scf,                       ONLY: ot_scf_destroy,&
                                             ot_scf_init,&
                                             ot_scf_mini,&
                                             ot_scf_read_input
  USE qs_overlap,                      ONLY: write_fm_with_basis_info,&
                                             write_sparse_matrix
  USE qs_rho_atom_methods,             ONLY: calculate_rho_atom_coeff
  USE qs_rho_methods,                  ONLY: qs_rho_update_rho, &
                                             duplicate_rho_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE qs_scf_post,                     ONLY: scf_post_calculation
  USE qs_scf_types,                    ONLY: ao_mixing_nr,&
                                             special_diag_method_nr,&
                                             general_diag_method_nr,&
                                             no_mixing_nr,&
                                             ot_method_nr,&
                                             qs_scf_env_type,&
                                             scf_env_release
  USE qs_wf_history_methods,           ONLY: wfi_extrapolate,&
                                             wfi_get_method_label,&
                                             wfi_update
  USE qs_wf_history_types,             ONLY: wfi_use_guess_method_nr, &
                                             wfi_frozen_method_nr
  USE scf_control_types,               ONLY: scf_control_type
  USE sparse_matrix_types,             ONLY: copy_matrix,&
                                             cp_sm_scale_and_add,&
                                             deallocate_matrix_set,&
                                             first_block_node,&
                                             get_block_node,&
                                             get_matrix_diagonal,&
                                             get_matrix_info,&
                                             next_block_node,&
                                             real_block_node_type,&
                                             real_matrix_p_type,&
                                             real_matrix_type,&
                                             set_matrix_diagonal,&
                                             transfer_matrix
  USE termination,                     ONLY: external_control,& 
                                             stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE harris_functional,               ONLY: harris_energy_correction, &
                                             harris_eigenvalue_trace_KS_Pmix, &
                                             harris_eigenvalue_calculation, &
                                             harris_postprocessing
  USE harris_env_types,                ONLY: harris_env_type

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: moduleN = "qs_scf"

  LOGICAL, PARAMETER, PRIVATE :: debug_this_module=.FALSE.

! *** Global variables ***

  INTEGER, SAVE, PRIVATE :: last_scf_env_id=0

! *** Public subroutines ***

  PUBLIC :: eigensolver, simple_eigensolver, scf,&
            scf_env_cleanup, scf_env_density_mixing,&
            scf_env_create, scf_env_do_scf, &
            scf_env_did_change, init_scf_run

! *****************************************************************************

CONTAINS

! *****************************************************************************

!!****f* qs_scf/scf_env_create [1.0] *
!!
!!   NAME
!!     scf_env_create
!!
!!   SYNOPSIS
!!     Subroutine scf_env_create(scf_env, qs_env, error)
!!       Type(qs_scf_env_type), Pointer:: scf_env
!!       Type(qs_environment_type), Intent (IN), Target:: qs_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine scf_env_create
!!
!!   FUNCTION
!!     allocates and initialize an scf_env 
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - scf_env: the scf env to initialize
!!     - qs_env: the qs_env this scf lives in
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE scf_env_create(scf_env, qs_env, error)
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scf_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(scf_control_type), POINTER          :: scf_control

    failure=.FALSE.
    NULLIFY(scf_control,dft_control)

    ALLOCATE(scf_env, stat=stat)
    CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    IF (.NOT.failure) THEN
       CALL get_qs_env(qs_env=qs_env,&
            scf_control=scf_control,&
            dft_control=dft_control)

       scf_env%ref_count=1
       scf_env%print_count=0
       last_scf_env_id=last_scf_env_id+1
       scf_env%id_nr=last_scf_env_id
       scf_env%print_count=0
       scf_env%iter_count=0
       scf_env%nelectron=0

! sets the method
       scf_env%method=general_diag_method_nr ! default with diagonalisation
       IF (dft_control%qs_control%semi_empirical) THEN
         IF (dft_control%qs_control%se_control%orthogonal_basis) &
           scf_env%method=special_diag_method_nr
       END IF
       IF (scf_control%use_ot) scf_env%method=ot_method_nr

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

       NULLIFY(scf_env%scf_work1, scf_env%scf_work2, scf_env%ortho,&
            scf_env%p_mix_new, scf_env%ot_preconditioner, scf_env%qs_ot_env,&
            scf_env%scf_diis_buffer)

    END IF
  END SUBROUTINE scf_env_create
!***************************************************************************

!!****f* qs_scf/scf [1.0] *
!!
!!   NAME
!!     scf
!!
!!   SYNOPSIS
!!     Subroutine scf(qs_env, globenv, error)
!!       Type(qs_environment_type), Intent (INOUT), Target:: qs_env
!!       Type(global_environment_type), Intent (IN):: globenv
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine scf
!!
!!   FUNCTION
!!     perform an scf procedure in the given qs_env
!!
!!   NOTES
!!     globenv should be removed
!!
!!   ARGUMENTS
!!     - qs_env: the qs_environment where to perform the scf procedure
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 introduced scf_env, moved real work to scf_env_do_scf [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE scf(qs_env,globenv,calc_forces,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), &
      POINTER                                :: globenv
    LOGICAL, OPTIONAL, INTENT(IN)            :: calc_forces
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scf', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(qs_environment_type), POINTER       :: my_qs_env
    TYPE(qs_scf_env_type), POINTER           :: scf_env

    NULLIFY(scf_env,my_qs_env)
    my_qs_env => qs_env
    failure=.FALSE.
    ! CPPrecondition(associated(qs_env),cp_failure_level,routineP,error,failure)
    NULLIFY(scf_env)
    IF (.NOT.failure) THEN
       CALL get_qs_env(qs_env,scf_env=scf_env,error=error)
       IF (.NOT.ASSOCIATED(scf_env)) THEN
          CALL scf_env_create(scf_env, qs_env=qs_env, error=error)
          CALL set_qs_env(qs_env,scf_env=scf_env,error=error)
          CALL scf_env_release(scf_env)
          CALL get_qs_env(qs_env,scf_env=scf_env,error=error)
       END IF
       IF (PRESENT(calc_forces)) THEN
         CALL scf_env_do_scf(scf_env=scf_env, qs_env=my_qs_env, globenv=globenv, &
                             calc_forces=calc_forces, error=error)
       ELSE
         CALL scf_env_do_scf(scf_env=scf_env, qs_env=my_qs_env, &
                             globenv=globenv, error=error)
       END IF
    END IF
  END SUBROUTINE scf
!***************************************************************************

!!****f* qs_scf/scf_env_do_scf [1.0] *
!!
!!   NAME
!!     scf_env_do_scf
!!
!!   SYNOPSIS
!!     Subroutine scf_env_do_scf(scf_env, qs_env, globenv, error)
!!       Type(qs_scf_env_type), Pointer:: scf_env
!!       Type(global_environment_type), Intent (IN):: globenv
!!       Type(qs_environment_type), Pointer:: qs_env
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine scf_env_do_scf
!!
!!   FUNCTION
!!     perform an scf loop
!!
!!   NOTES
!!     globenv should be removed
!!
!!   ARGUMENTS
!!     - scf_env: the scf_env where to perform the scf procedure
!!     - qs_env: the qs_env, the scf_env lives in
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Matthias Krack
!!
!!   MODIFICATION HISTORY
!!     long history, see cvs and qs_scf module history
!!     02.2003 introduced scf_env [fawzi]
!!     09.2005 Frozen density approximation [TdK]
!!
!!*** **********************************************************************
  SUBROUTINE scf_env_do_scf(scf_env,qs_env,globenv,calc_forces,error)

    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), &
      POINTER                                :: globenv
    LOGICAL, INTENT(IN), OPTIONAL            :: calc_forces
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "scf_env_do_scf", &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, handle2, ispin, &
                                                output_unit
    LOGICAL :: diis_step, do_level_shift, energy_only, failure, gapw, gapw_xc, ionode, &
      should_stop, use_cholesky, use_jacobi, id_equal
    REAL(KIND=dp)                            :: diis_error, orthonormality, &
                                                t1, t2, tot1_h, tot1_s
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: scf_work1
    TYPE(cp_fm_type), POINTER                :: ortho, scf_work2
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(qs_charges_type), POINTER           :: qs_charges
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_rho_type), POINTER               :: rho,rho_xc
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_h, matrix_ks, matrix_p, matrix_s
    TYPE(scf_control_type), POINTER          :: scf_control
    LOGICAL                                  :: has_unit_metric
    TYPE(real_matrix_type), POINTER          :: orthogonality_metric
    LOGICAL                                  :: harris_flag, &
                                                EII_necessary
    TYPE(harris_env_type), POINTER           :: harris_env

!   ---------------------------------------------------------------------------

    NULLIFY(dft_control,matrix_h,matrix_s,matrix_ks,matrix_p,rho,rho_xc,energy, &
            scf_control,logger,pw_env, qs_charges,ks_env,mos,atomic_kind_set, &
            particle_set,scf_work1,ortho,scf_work2,harris_env)

    logger => cp_error_get_logger(error)
    t1 = m_walltime()
    failure=.FALSE.

    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(qs_env%ref_count>0,cp_failure_level,routineP,error,failure)
    
    CALL get_qs_env(qs_env=qs_env,&
         atomic_kind_set=atomic_kind_set,&
         dft_control=dft_control,&
         scf_control=scf_control)

!   *** Quick return, if no SCF iteration is requested ***

    IF (scf_control%max_scf < 1) RETURN

    CALL write_checkpoint_information("entering "//routineN,globenv)

    CALL timeset(routineN,"I"," ",handle)

    ionode = logger%para_env%source==logger%para_env%mepos
    IF (ionode) output_unit = cp_logger_get_default_unit_nr(logger)
    energy_only = .FALSE.
    EII_necessary = .TRUE.

    ! some special checks to eliminate the first problems with restricted
    IF (dft_control%restricted) THEN
       ! OT only right now
       IF (scf_env%method .NE. ot_method_nr) THEN 
          IF (ionode) write(output_unit,*) "OT only for restricted"
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       ENDIF
       CPPrecondition(dft_control%nspins.EQ.2,cp_failure_level,routineP,error,failure)
    ENDIF

    IF (ionode.AND.logger%print_keys%scf) THEN
       WRITE (UNIT=output_unit,FMT="(/,/,T2,A)")&
            "SCF WAVEFUNCTION OPTIMIZATION"
    END IF

!   True the GAPW method is used
    gapw = dft_control%qs_control%gapw
    gapw_xc = dft_control%qs_control%gapw_xc

    harris_flag = qs_env%use_harris

    IF (harris_flag) THEN
      CALL get_qs_env(qs_env=qs_env, harris_env=harris_env, error=error)
    END IF

    IF (PRESENT(calc_forces)) THEN
      CALL init_scf_run(scf_env=scf_env, qs_env=qs_env, globenv=globenv, &
                        calc_forces=calc_forces, error=error)
    ELSE
      CALL init_scf_run(scf_env=scf_env, qs_env=qs_env, &
                        globenv=globenv, error=error)
    END IF

    CALL get_qs_env(qs_env=qs_env,&
         matrix_h=matrix_h,&
         matrix_s=matrix_s,energy=energy,&
         particle_set=particle_set,&
         pw_env=pw_env,&
         qs_charges=qs_charges,&
         ks_env=ks_env, &
         error=error)

     has_unit_metric=.FALSE.
     IF (dft_control%qs_control%semi_empirical) THEN
         IF (dft_control%qs_control%se_control%orthogonal_basis) THEN
             has_unit_metric=.TRUE.
         ENDIF
     END IF


    scf_work1 => scf_env%scf_work1
    scf_work2 => scf_env%scf_work2
    ortho => scf_env%ortho

    CALL get_qs_env(qs_env,matrix_ks=matrix_ks,rho=rho,rho_xc=rho_xc,mos=mos)
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
    
    IF (.NOT.logger%print_keys%timecumul) t1 = m_walltime()

!   *** SCF loop ***

    scf_loop: DO
       CALL timeset("scf_iter","I"," ",handle2)

       IF (ionode) CALL m_flush(output_unit)

       scf_env%iter_count = scf_env%iter_count + 1

! ** here qs_env%rho%rho_r and qs_env%rho%rho_g should be up to date
! ** should also be the case for restricted calculations

       CALL qs_ks_update_qs_env(ks_env,qs_env=qs_env,&
            error=error,&
            calculate_forces=.FALSE.,just_energy=energy_only)

       ! Frozen density approximation
       IF (qs_env%wf_history%interpolation_method_nr==wfi_frozen_method_nr) THEN
         IF (.NOT. ASSOCIATED(qs_env%wf_history%past_states(1)%snapshot)) THEN
           IF (ASSOCIATED(qs_env%wf_history)) THEN
             CALL wfi_update(qs_env%wf_history,qs_env=qs_env,dt=1.0_dp,&
                             error=error)
           END IF
           CALL duplicate_rho_type(rho_input=rho, &
                                   rho_output=qs_env%wf_history%past_states(1)%snapshot%rho_frozen, &
                                   qs_env=qs_env, error=error)

           !CALL mpools_get(qs_env%mpools, S_sm_pools=S_sm_pools, error=error)
           !CALL sm_pools_create_matrix_vect(S_sm_pools, qs_env%wf_history%past_states(1)%snapshot%rho_ao, &
           !                                 error=error)
           !DO ispin = 1,dft_control%nspins
           !  CALL copy_matrix(rho%rho_ao(ispin)%matrix, &
           !       qs_env%wf_history%past_states(1)%snapshot%rho_ao(ispin)%matrix)
           !END DO

         END IF
       END IF

!FM      CALL pw_env_flush_cache(pw_env,error=error)

       DO ispin=1,dft_control%nspins

          IF (logger%print_keys%density_matrix.AND.&
               logger%print_keys%each_scf_step) THEN
             CALL write_sparse_matrix(matrix_p(ispin)%matrix,4,6,qs_env,globenv)
          END IF

          IF (logger%print_keys%kohn_sham_matrix.AND.logger%print_keys%each_scf_step) THEN
             IF (dft_control%qs_control%semi_empirical) THEN
               CALL write_sparse_matrix(matrix_ks(ispin)%matrix,4,6,qs_env,globenv,&
                    scale=evolt)
             ELSE
               CALL write_sparse_matrix(matrix_ks(ispin)%matrix,4,6,qs_env,globenv)
             END IF
          END IF

       ENDDO

       IF (logger%print_keys%mo_orthonormality) THEN
          IF(scf_env%method==special_diag_method_nr) THEN
            CALL calculate_orthonormality(orthonormality,mos)
          ELSE
            CALL calculate_orthonormality(orthonormality,mos,matrix_s(1)%matrix)
          END IF
          IF (ionode) THEN
             WRITE(output_unit,'(T8,A,T60,E20.4)')  &
                    " Maximum deviation from MO S-orthonormality",orthonormality
          ENDIF
       ENDIF

       scf_env%p_mix = scf_control%p_mix
       scf_env%iter_param = 0.0_dp
       IF (scf_env%mixing_method==ao_mixing_nr) scf_env%iter_param = scf_env%p_mix

       SELECT CASE (scf_env%method)
         CASE DEFAULT
            CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                 routineP,"unknown scf method method:"//&
                 cp_to_string(scf_env%method),error,failure)

         CASE(general_diag_method_nr) ! diagonalisation (default)

            DO ispin=1,dft_control%nspins 
             CALL copy_sm_to_fm(matrix_ks(ispin)%matrix,scf_work1(ispin)%matrix)
            ENDDO

            IF (scf_env%iter_count > 1) THEN
               CALL qs_diis_b_step(scf_env%scf_diis_buffer,mos,scf_work1,&
                    scf_work2,scf_env%iter_delta,diis_error,diis_step,&
                    scf_control%eps_diis,s_matrix=matrix_s,error=error)
            END IF

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

            IF ((scf_env%iter_count > 1).AND.(scf_env%iter_delta < scf_control%eps_jacobi)) THEN
               use_jacobi = .TRUE.
            ELSE
               use_cholesky = scf_control%use_cholesky
               use_jacobi = .FALSE.
            END IF

            scf_env%iter_delta=0.0_dp

            DO ispin=1,dft_control%nspins

               CALL eigensolver(scf_work1(ispin)%matrix, mos(ispin)%mo_set,ortho,scf_work2,&
                    do_level_shift,scf_control%level_shift,&
                    use_cholesky=use_cholesky,&
                    work_syevx=scf_control%work_syevx,&
                    use_jacobi=use_jacobi,&
                    jacobi_threshold=scf_control%jacobi_threshold,&
                    smear=scf_control%smear, &
                    error=error)

               IF (logger%print_keys%each_scf_step) THEN
                  CALL write_mo_set(mos(ispin)%mo_set,atomic_kind_set,particle_set,4,6,globenv)
               END IF

               CALL calculate_density_matrix(mos(ispin)%mo_set,&
                    scf_env%p_mix_new(ispin)%matrix,error=error)

            ENDDO

         CASE(special_diag_method_nr)
           ! diagonalization for orthogonal basis sets
            DO ispin=1,dft_control%nspins
              CALL copy_sm_to_fm(matrix_ks(ispin)%matrix,scf_work1(ispin)%matrix)
            ENDDO
            IF (scf_env%iter_count > 1) THEN
               CALL qs_diis_b_step(scf_env%scf_diis_buffer,mos,scf_work1,&
                    scf_work2,scf_env%iter_delta,diis_error,diis_step,&
                    scf_control%eps_diis,error=error)
            END IF
            do_level_shift = ((scf_control%level_shift /= 0.0_dp).AND.&
               ((scf_control%density_guess == "CORE").OR.(scf_env%iter_count > 1)))
            IF (diis_step) THEN
               scf_env%p_mix = 1.0_dp
               scf_env%iter_param = diis_error
               scf_env%iter_method = "DIIS/Diag."
            ELSE
               scf_env%iter_method = "Mixing/Diag."
            END IF
            scf_env%iter_delta=0.0_dp
            DO ispin=1,dft_control%nspins
               CALL simple_eigensolver(scf_work1(ispin)%matrix,&
                    mos(ispin)%mo_set,do_level_shift,scf_control%level_shift,&
                    work_syevx=scf_control%work_syevx,&
                    smear=scf_control%smear,error=error)
               IF (logger%print_keys%each_scf_step) THEN
                  CALL write_mo_set(mos(ispin)%mo_set,atomic_kind_set,&
                       particle_set,4,6,globenv)
               END IF
               CALL calculate_density_matrix(mos(ispin)%mo_set,&
                    scf_env%p_mix_new(ispin)%matrix,error=error)
            ENDDO

         CASE(ot_method_nr) ! orbital transforms

! in case of LSD the first spin qs_ot_env will drive the minimization
! in the case of a restricted calculation, it will make sure the spin orbitals are equal
            IF (has_unit_metric) THEN
               NULLIFY(orthogonality_metric)
            ELSE
               orthogonality_metric=>matrix_s(1)%matrix
            ENDIF

            CALL ot_scf_mini(mos,qs_env%mo_derivs,orthogonality_metric, &
                             energy%total, energy_only,scf_env%iter_delta, &
                             scf_env%qs_ot_env)

            DO ispin=1,dft_control%nspins
               CALL calculate_density_matrix(mos(ispin)%mo_set,&
                     rho%rho_ao(ispin)%matrix,error=error)
            ENDDO

            scf_env%iter_method=scf_env%qs_ot_env(1)%OT_METHOD_FULL
            scf_env%iter_param=scf_env%qs_ot_env(1)%ds_min

       END SELECT

       ! ** calculation of the harris energy correction *** !
       IF (harris_flag) THEN
         CALL harris_energy_correction(qs_env, harris_env, &
                                       EII_necessary=EII_necessary, &
                                       globenv=globenv, fast=.TRUE., error=error)
         IF (scf_env%method .NE. ot_method_nr) THEN
           CALL harris_eigenvalue_trace_KS_Pmix(scf_env, qs_env, harris_env, error=error)
         ELSE
           CALL harris_eigenvalue_calculation(qs_env=qs_env, harris_env=harris_env, error=error)
         END IF
         !energy%total = harris_env%harris_energy%Eharris
       END IF

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

       IF (harris_flag) THEN
          CALL harris_postprocessing(harris_env, error=error)
       END IF

! ** convergence check
       CALL external_control(should_stop,"SCF",error)
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

!   *** Write restart file ***
       IF (scf_control%write_restart_each>0.AND.&
            MODULO(scf_env%iter_count,scf_control%write_restart_each)==0) THEN
          CALL write_mo_set(mos,atomic_kind_set,particle_set,globenv,&
               id_nr=1)
       END IF

       IF (.NOT.logger%print_keys%timecumul) t1 = m_walltime()

!   *** mixing methods have the new density matrix in p_mix_new
       IF (scf_env%mixing_method.EQ.ao_mixing_nr) THEN
           DO ispin=1,dft_control%nspins
              CALL transfer_matrix(scf_env%p_mix_new(ispin)%matrix,&
                                       rho%rho_ao(ispin)%matrix)
           END DO
       ENDIF

! ** update qs_env%rho
       CALL qs_rho_update_rho(rho, qs_env=qs_env, gapw=gapw, error=error)
       IF(gapw ) THEN
         CALL calculate_rho_atom_coeff(qs_env,rho%rho_ao)
       ENDIF
       IF(gapw_xc ) THEN
         CALL calculate_rho_atom_coeff(qs_env,rho%rho_ao)
         CALL qs_rho_update_rho(rho_xc, qs_env=qs_env, gapw=gapw_xc, error=error)
       ENDIF

       CALL qs_ks_did_change(ks_env,rho_changed=.TRUE.,error=error)

       CALL timestop(0.0_dp,handle2)

    END DO scf_loop

    IF (ionode.AND.logger%print_keys%scf) THEN
       WRITE (UNIT=output_unit,FMT="(/,(T3,A,T40,2F20.10))")&
            "Total electronic density (r-space): ",&
            SUM(rho%tot_rho_r),&
            SUM(rho%tot_rho_r)+ REAL(scf_env%nelectron,dp),&
            "Total core charge density (r-space):",&
            qs_charges%total_rho_core_rspace,&
            qs_charges%total_rho_core_rspace - REAL(scf_env%nelectron+dft_control%charge,dp)
       IF(gapw) THEN
          tot1_h =  qs_charges%total_rho1_hard(1)
          tot1_s =  qs_charges%total_rho1_soft(1)
          DO ispin=2,dft_control%nspins
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
            "Total charge density (g-space):     ",qs_charges%total_rho_gspace
       END IF
       IF (dft_control%qs_control%semi_empirical) THEN
         WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
            "Core-core repulsion energy [eV]:               ",energy%core_overlap*evolt,&
            "Core Hamiltonian energy [eV]:                  ",energy%core*evolt,&
            "Two-electron integral energy [eV]:             ",energy%hartree*evolt,&
            "Electronic energy [eV]:                        ",&
                                              (energy%core+0.5_dp*energy%hartree)*evolt
       ELSE
         WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
            "Overlap energy of the core charge distribution:",energy%core_overlap,&
            "Self energy of the core charge distribution:   ",energy%core_self,&
            "Core Hamiltonian energy:                       ",energy%core,&
            "Hartree energy:                                ",energy%hartree,&
            "Exchange-correlation energy:                   ",energy%exc
         IF(gapw) THEN
           WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
              "GAPW| Exc from hard and soft atomic rho1:      ",energy%exc1,&
              "GAPW| local Eh = 1 center integrals:           ",energy%hartree_1c,&
              "GAPW| local Eh = 2 center integrals:           ",energy%hartree_2c,&
              "GAPW| local Eh = 3 centers integrals:          ",energy%hartree_3c
         END IF
         IF(gapw_xc) THEN
           WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
              "GAPW_XC| Exc from hard and soft atomic rho1:      ",energy%exc1
         END IF
       END IF
       IF  (qs_env%qmmm) THEN
          WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
               "QM/MM Electrostatic energy:                    ",energy%qmmm_el
       END IF
       IF (dft_control%qs_control%mulliken_restraint) THEN
          WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
                "Mulliken restraint energy: ",energy%mulliken
       END IF
       IF (dft_control%qs_control%semi_empirical) THEN
         WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
            "Total energy [eV]:                             ",energy%total*evolt
         WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
            "Atomic reference energy [eV]:                  ",energy%core_self*evolt,&
            "Heat of formation [kcal/mol]:                  ",&
            (energy%total+energy%core_self)*kcalmol
       ELSE
         WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
            "Total energy:                                  ",energy%total
       END IF
       CALL m_flush(output_unit)
    END IF

    IF (harris_flag) THEN
      energy%total = harris_env%harris_energy%Eharris
    END IF

!   *** Write restart file ***
    IF (scf_control%write_restart_each>=0) THEN
       CALL write_mo_set(mos,atomic_kind_set,particle_set,globenv,&
            id_nr=1)
    END IF

!   *** add the converged wavefunction to the wavefunction history
!   should possibly be moved down after de-mixing
    IF (ASSOCIATED(qs_env%wf_history)) THEN
        CALL wfi_update(qs_env%wf_history,qs_env=qs_env,dt=1.0_dp,&
               error=error)
    END IF

!   *** mixing methods need to undo mixing of the density matrix (restore original density) ***
    IF (scf_env%mixing_method.EQ.ao_mixing_nr) THEN
        CALL scf_env_density_mixing(scf_env%p_mix_new,scf_env%mixing_method,&
                                    scf_env%p_mix,scf_env%iter_delta,&
                                    qs_env=qs_env,invert=.TRUE.,error=error)
       DO ispin=1,dft_control%nspins
         CALL transfer_matrix(scf_env%p_mix_new(ispin)%matrix,&
                              rho%rho_ao(ispin)%matrix)
       END DO
    ENDIF

!   *** update rspace rho since the mo changed
!   *** this might not always be needed (i.e. no post calculation / no forces )
    CALL qs_rho_update_rho(rho, qs_env=qs_env, gapw = gapw, error=error)
    IF(gapw) THEN
       CALL calculate_rho_atom_coeff(qs_env,rho%rho_ao)
    END IF
    IF(gapw_xc ) THEN
       CALL calculate_rho_atom_coeff(qs_env,rho%rho_ao)
       CALL qs_rho_update_rho(rho_xc, qs_env=qs_env, gapw=gapw_xc, error=error)
    ENDIF

    CALL qs_ks_did_change(ks_env,rho_changed=.TRUE.,error=error)

! *** print out info
    IF ((.NOT.should_stop).AND.&
        (.NOT.dft_control%qs_control%semi_empirical)) THEN
       CALL scf_post_calculation(scf_env=scf_env, qs_env=qs_env, &
            globenv=globenv, error=error)
    END IF

! *** cleanup
    CALL scf_env_cleanup(scf_env,qs_env=qs_env,error=error)

    CALL timestop(0.0_dp,handle)

    CALL write_checkpoint_information("leaving "//routineN,globenv)

  END SUBROUTINE scf_env_do_scf

!!****f* qs_scf/scf_env_density_mixing [1.0] *
!!
!!   NAME
!!     scf_env_density_mixing
!!
!!   FUNCTION
!!     perform (if requested) a density mixing
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE scf_env_density_mixing(p_mix_new,mixing_method,p_mix,iter_delta,&
                                    qs_env,invert,error)

    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: p_mix_new 
    INTEGER, INTENT(IN)                      :: mixing_method
    REAL(KIND = dp), INTENT(IN)              :: p_mix
    REAL(KIND = dp), INTENT(INOUT)           :: iter_delta
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(in), OPTIONAL            :: invert
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scf_env_density_mixing', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin
    LOGICAL                                  :: failure, my_invert
    REAL(KIND=dp)                            :: tmp

    failure=.FALSE.
    my_invert=.FALSE.
    IF (PRESENT(invert)) my_invert=invert

    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(qs_env%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       SELECT CASE (mixing_method)
       CASE(no_mixing_nr)
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP,"don't call mixing in no_mixing_nr:",error,failure)
       CASE(ao_mixing_nr)
          CPPrecondition(ASSOCIATED(p_mix_new),cp_failure_level,routineP,error,failure)
!FM       iter_delta=0.0_dp
          DO ispin=1,SIZE(p_mix_new)
             IF (my_invert) THEN
                CPPrecondition(p_mix/=0.0_dp,cp_failure_level,routineP,error,failure)
                IF (p_mix/=1.0_dp) THEN
                   CALL cp_sm_scale_and_add(matrix_a=p_mix_new(ispin)%matrix,&
                        alpha=1.0_dp/p_mix,&
                        matrix_b=qs_env%rho%rho_ao(ispin)%matrix,&
                        beta=(p_mix-1.0_dp)/p_mix,&
                        error=error)
                END IF
             ELSE
                CALL cp_sm_mix(m1=p_mix_new(ispin)%matrix,&
                     m2=qs_env%rho%rho_ao(ispin)%matrix,&
                     p_mix=p_mix,&
                     delta=tmp,&
                     para_env=qs_env%para_env,&
                     error=error)
                iter_delta=MAX(iter_delta,tmp)
             END IF
          END DO
       CASE default
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP,"unknown density mixing method:"//&
               cp_to_string(mixing_method),error,failure)
       END SELECT
    END IF
  END SUBROUTINE scf_env_density_mixing
!***************************************************************************

! *****************************************************************************

  SUBROUTINE init_scf_run(scf_env,qs_env,globenv,calc_forces,error)

!   Purpose: Initialise a SCF run.

!   History: - Creation (30.04.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), &
      POINTER                                :: globenv
    LOGICAL, INTENT(IN), OPTIONAL            :: calc_forces
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "init_scf_run", &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, homo, ispin, nao, &
                                                ndep, nelectron_spin, &
                                                output_unit, &
                                                stat, number_of_OT_envs
    LOGICAL                                  :: failure, &
                                                gth_potential_present, ionode,&
                                                my_transition_potential
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_fm_type), POINTER                :: mo_coeff, ortho, &
                                                scf_work2, e_vectors
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_h, matrix_ks, matrix_s
    TYPE(scf_control_type), POINTER          :: scf_control
    LOGICAL                                  :: has_unit_metric,do_rotation
    TYPE(real_matrix_type), POINTER          :: orthogonality_metric

!dbg 
    INTEGER :: i
!   ---------------------------------------------------------------------------

    CALL timeset("init_scf_run","I"," ",handle)

    NULLIFY(scf_control,atomic_kind_set,matrix_h,matrix_s,matrix_ks,dft_control,mos,&
         ortho,scf_work2, mo_coeff)
    failure=.FALSE.
    my_transition_potential = .FALSE.
    logger => cp_error_get_logger(error)
    
    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(qs_env%ref_count>0,cp_failure_level,routineP,error,failure)

    CALL scf_env_check_i_alloc(scf_env=scf_env, qs_env=qs_env,&
         error=error)

    CALL get_qs_env(qs_env=qs_env,&
         scf_control=scf_control,&
         dft_control=dft_control,&
         atomic_kind_set=atomic_kind_set,&
         mos=mos,matrix_ks=matrix_ks,&
         matrix_h=matrix_h,matrix_s=matrix_s, error=error)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         gth_potential_present=gth_potential_present)

    has_unit_metric=.FALSE.
    IF (dft_control%qs_control%semi_empirical) THEN
        IF (dft_control%qs_control%se_control%orthogonal_basis) THEN
            has_unit_metric=.TRUE.
        ENDIF
    END IF

    ionode = logger%para_env%source == logger%para_env%mepos
    IF (ionode) output_unit = cp_logger_get_default_unit_nr(logger)

    scf_work2 => scf_env%scf_work2
    ortho => scf_env%ortho

! update ecore
    IF (.NOT.dft_control%qs_control%semi_empirical) THEN
      CALL calculate_ecore_self(qs_env)
      IF (PRESENT(calc_forces)) THEN
        CALL calculate_ecore_overlap(qs_env, globenv, &
                                     calculate_forces=calc_forces)
      ELSE
        CALL calculate_ecore_overlap(qs_env, globenv, &
                                     calculate_forces=dft_control%forces)
      END IF
      IF (ionode.AND.logger%print_keys%scf_energies) THEN
         WRITE (UNIT=output_unit,FMT="(/,(T3,A,T60,F20.10))")&
              "Overlap energy of the core charge distribution:",&
              qs_env%energy%core_overlap,&
              "Self energy of the core charge distribution:   ",&
              qs_env%energy%core_self
      END IF
    END IF

! updates the total number of electrons
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         nelectron=scf_env%nelectron)
    scf_env%nelectron = scf_env%nelectron - dft_control%charge
! print occupation numbers
    IF (ionode.AND.logger%print_keys%scf) THEN
       DO ispin=1,dft_control%nspins
          CALL get_mo_set(mo_set=mos(ispin)%mo_set,&
               homo=homo, &
               nelectron=nelectron_spin, &
               nao=nao)
          IF (dft_control%nspins > 1) THEN
            WRITE (UNIT=output_unit,FMT="(/,T3,A,I2)") "Spin",ispin
          END IF
          WRITE (UNIT=output_unit,FMT="(/,(T3,A,T31,I10))")&
            "Number of electrons:",nelectron_spin,&
            "Number of occupied orbitals:",homo
       ENDDO
       WRITE (UNIT=output_unit,FMT="((T3,A,T31,I10))")&
         "Number of orbital functions:",nao
    END IF

! calc ortho matrix
    ndep = 0
    IF (scf_env%method /= ot_method_nr .AND. &
        scf_env%method /= special_diag_method_nr ) THEN
       IF (scf_control%use_cholesky ) THEN
          CALL copy_sm_to_fm(matrix_s(1)%matrix,ortho)
          CALL cp_fm_cholesky_decompose(ortho)
       ELSE
          CALL calculate_ortho_matrix(ortho,matrix_s(1)%matrix,scf_work2,ndep,&
               scf_control%eps_eigval,scf_control%work_syevx)
          IF (ionode.and.logger%print_keys%scf) THEN
             WRITE (UNIT=output_unit,FMT="(T3,A,I10)")&
                  "Number of independent orbital functions:",nao - ndep
          END IF
       ENDIF
    ENDIF
    IF (logger%print_keys%ortho_matrix) THEN
       CALL write_fm_with_basis_info(ortho,4,6,qs_env,globenv)
    END IF

    my_transition_potential = (dft_control%do_xas_calculation .AND. dft_control%xas_estate>0)
    ! initializes rho and the mos
    IF( my_transition_potential) THEN
      ! initialize the density with the localized mos
        CALL loc_initialize_rho(qs_env,error=error)  
    ELSE
      CALL scf_env_initial_rho_setup(scf_env,qs_env=qs_env,&
           globenv=globenv,error=error)
    END IF

! *** method dependent initializations ***
    SELECT CASE (scf_env%method)
    CASE DEFAULT
      CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP,"unknown scf method method:"//&
               cp_to_string(scf_env%method),error,failure)
    CASE (general_diag_method_nr,special_diag_method_nr)
       IF (.NOT.ASSOCIATED(scf_env%scf_diis_buffer)) THEN       
          CALL qs_diis_b_create(scf_env%scf_diis_buffer,&
               nbuffer=scf_control%max_diis,&
               error=error)
       END IF
       CALL qs_diis_b_clear(scf_env%scf_diis_buffer,error=error)
    CASE (ot_method_nr)
       IF (.NOT.ASSOCIATED(scf_env%ot_preconditioner)) THEN
          ALLOCATE(scf_env%ot_preconditioner, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          CALL init_preconditioner(scf_env%ot_preconditioner,&
               para_env=qs_env%para_env,&
               blacs_env=qs_env%blacs_env)
       END IF
       IF (.NOT.ASSOCIATED(scf_env%qs_ot_env)) THEN

          ! restricted calculations require just on set of OT orbitals
          number_of_OT_envs=dft_control%nspins
          IF (dft_control%restricted) number_of_OT_envs=1

          ALLOCATE(scf_env%qs_ot_env(number_of_OT_envs),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

          ! XXX Joost XXX should disentangle reading input from this part
          CALL ot_scf_read_input(scf_env%qs_ot_env,globenv)

          ! keep a note that we are restricted
          IF (dft_control%restricted) THEN
              scf_env%qs_ot_env(1)%restricted=.TRUE.
              ! requires rotation
              CPPrecondition(scf_env%qs_ot_env(1)%settings%do_rotation,cp_failure_level,routineP,error,failure) 
          ELSE 
              scf_env%qs_ot_env(:)%restricted=.FALSE.
          ENDIF


          ! might need the KS matrix to init properly
          CALL qs_ks_update_qs_env(qs_env%ks_env,qs_env=qs_env,&
               error=error, calculate_forces=.FALSE.,just_energy=.FALSE.)

          IF (scf_env%qs_ot_env(1)%settings%preconditioner_type /= ot_precond_none) THEN 
             !
             ! preconditioning still needs to be done correctly with has_unit_metric
             ! notice that a big part of the preconditioning (S^-1) is fine anyhow
             !
             CPPrecondition(.NOT. has_unit_metric,cp_failure_level,routineP,error,failure) 
             ! constuct the preconditioner based on the first spin set. 
             ! this might or might not be reasonable. Presumably it is.
             CALL get_mo_set(mo_set=mos(1)%mo_set,&
                  mo_coeff=mo_coeff)

             IF (scf_env%qs_ot_env(1)%settings%preconditioner_type == ot_precond_full_all) THEN
                ! this will rotate the MOs to be eigen states, which is not compatible with rotation
                do_rotation=scf_env%qs_ot_env(1)%settings%do_rotation
                CPPrecondition(.NOT.do_rotation,cp_failure_level,routineP,error,failure)
                ! this still needs to be fixed, LSD, requires two separate preconditioners 
                CPPrecondition(dft_control%nspins==1,cp_failure_level,routineP,error,failure)
                ! here we also rotate the derivate, not only the c0
                CALL calculate_subspace_eigenvalues(mo_coeff,matrix_ks(1)%matrix,&
                            para_env=mo_coeff%matrix_struct%para_env, &
                            do_rotation = .TRUE., co_rotate=qs_env%mo_derivs(1)%matrix)
             ENDIF

             IF(scf_env%qs_ot_env(1)%settings%nprec ==1 ) THEN
               CALL make_preconditioner(scf_env%ot_preconditioner, &
                  scf_env%qs_ot_env(1)%settings%preconditioner_type, &
                  matrix_ks(1)%matrix,matrix_s(1)%matrix,&
                  qs_env%kinetic(1)%matrix, &
                  mo_coeff,scf_env%qs_ot_env(1)%settings%energy_gap) 
             ELSE
               CALL make_preconditioner(scf_env%ot_preconditioner, &
                  scf_env%qs_ot_env(1)%settings%preconditioner_type, &
                  matrix_ks(1)%matrix,matrix_s(1)%matrix,&
                  qs_env%kinetic(1)%matrix, &
                  mo_coeff,scf_env%qs_ot_env(1)%settings%energy_gap,&
                  scf_env%qs_ot_env(1)%settings%nprec)
             END IF

          ENDIF

          IF (has_unit_metric) THEN
             NULLIFY(orthogonality_metric)
          ELSE
             orthogonality_metric=>matrix_s(1)%matrix
          ENDIF

          CALL ot_scf_init(mos,orthogonality_metric,scf_env%qs_ot_env)
          
          IF (scf_env%qs_ot_env(1)%settings%preconditioner_type /= ot_precond_none) THEN
             DO ispin=1,SIZE(scf_env%qs_ot_env)
                CALL qs_ot_new_preconditioner(scf_env%qs_ot_env(ispin),&
                     scf_env%ot_preconditioner,scf_env%qs_ot_env(1)%settings%nprec) 
             ENDDO
          ENDIF

       ENDIF
    END SELECT

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE init_scf_run

! *****************************************************************************

!!****f* qs_scf/scf_env_check_i_alloc [1.0] *
!!
!!   NAME
!!     scf_env_check_i_alloc
!!
!!   FUNCTION
!!     checks the allocation status of the needed matrixes, and if necessary
!!     allocate them
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - scf_env: the scf_env to be checked
!!     - qs_env: the qs_env, the scf_env lives in
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE scf_env_check_i_alloc(scf_env,qs_env,error)
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "scf_env_check_i_alloc", &
      routineP = moduleN//":"//routineN

    INTEGER                                  :: handle, homo, ispin, nao, &
                                                nrow_block, stat
    LOGICAL                                  :: failure, gth_potential_present
    LOGICAL                                  :: my_transition_potential,uniform_occupation
    REAL(KIND=dp), DIMENSION(:), POINTER     :: occupation_numbers
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: ao_mo_fm_pools
    TYPE(cp_fm_struct_type), POINTER         :: ao_ao_fmstruct, ao_mo_fmstruct
    TYPE(cp_sm_pool_p_type), DIMENSION(:), &
      POINTER                                :: S_sm_pools
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(cp_fm_p_type), POINTER, DIMENSION(:):: mo_derivs
    TYPE(cp_fm_type), POINTER                :: mo_coeff

!   ---------------------------------------------------------------------------

    CALL timeset(routineN,"I"," ",handle)

    NULLIFY(matrix_ks, ao_mo_fm_pools, matrix_s, ao_mo_fmstruct, ao_ao_fmstruct,&
            dft_control, mos, occupation_numbers, S_sm_pools, ks_env)
    NULLIFY(atomic_kind_set, mo_derivs, mo_coeff)

    failure=.FALSE.

    my_transition_potential = .FALSE.
    uniform_occupation = .TRUE.

    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)

    CALL get_qs_env(qs_env=qs_env,&
         dft_control=dft_control,&
         mos=mos,&
         matrix_ks=matrix_ks,&
         ks_env=ks_env,&
         atomic_kind_set=atomic_kind_set,&
         matrix_s=matrix_s)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         gth_potential_present=gth_potential_present)
    CALL mpools_get(qs_env%mpools, ao_mo_fm_pools=ao_mo_fm_pools,&
         S_sm_pools=S_sm_pools,&
         error=error)

    ! very first tests for xas

    my_transition_potential = (dft_control%do_xas_calculation .AND. dft_control%xas_estate>0)
!   *** finish initialization of the MOs ***
    CPPrecondition(ASSOCIATED(mos),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       DO ispin=1,SIZE(mos)
          CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff,homo=homo)
          IF (.NOT.ASSOCIATED(mo_coeff)) THEN
             CALL init_mo_set(mos(ispin)%mo_set,&
                  ao_mo_fm_pools(ispin)%pool,&
                  name="qs_env"//TRIM(ADJUSTL(cp_to_string(qs_env%id_nr)))//&
                  "%mo"//TRIM(ADJUSTL(cp_to_string(ispin))),&
                  error=error)
          END IF
          IF(my_transition_potential .AND. ispin==1) THEN
            CALL get_mo_set(mos(ispin)%mo_set,&
                 occupation_numbers=occupation_numbers)
            IF(dft_control%xas_control%xas_method == xas_tp_hh) THEN
               occupation_numbers(dft_control%xas_estate) = 0.5_dp
               uniform_occupation = .FALSE.
            ELSE IF(dft_control%xas_control%xas_method == xas_tp_fh) THEN
               occupation_numbers(dft_control%xas_estate) = 0.0_dp
               uniform_occupation = .FALSE.    
            ELSE IF(dft_control%xas_control%xas_method == xas_dscf) THEN
               occupation_numbers(dft_control%xas_estate) = 0.0_dp
               uniform_occupation = .FALSE.
               occupation_numbers(homo) = 1.0_dp
            END IF
            CALL set_mo_set(mos(ispin)%mo_set,&
                 uniform_occupation=uniform_occupation)
          END IF
       END DO
    END IF

!   *** get the mo_derivs OK if needed ***
    IF (qs_env%requires_mo_derivs) THEN
     CALL get_qs_env(qs_env,mo_derivs=mo_derivs)
     IF (.NOT.ASSOCIATED(mo_derivs)) THEN
       IF (dft_control%restricted) THEN ! right now, there might be more mos than needed derivs
          allocate(mo_derivs(1))
          CALL get_mo_set(mos(1)%mo_set,mo_coeff=mo_coeff)
          CALL cp_fm_create(mo_derivs(1)%matrix,mo_coeff%matrix_struct)
       ELSE
          allocate(mo_derivs(dft_control%nspins))
          DO ispin=1,dft_control%nspins
             CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff)
             CALL cp_fm_create(mo_derivs(ispin)%matrix,mo_coeff%matrix_struct)
          ENDDO
       ENDIF
       CALL set_qs_env(qs_env,mo_derivs=mo_derivs)
     ENDIF
    ELSE
       ! nothing should be done
    ENDIF


!   *** Allocate the distributed SCF matrices ***

    IF ((.NOT.ASSOCIATED(scf_env%scf_work1)).OR.&
        (.NOT.ASSOCIATED(scf_env%scf_work2)).OR.&
        (.NOT.ASSOCIATED(scf_env%ortho))) THEN

       ao_mo_fmstruct => fm_pool_get_el_struct(ao_mo_fm_pools(1)%pool,&
            error=error)
       CALL cp_fm_struct_get(ao_mo_fmstruct, nrow_block=nrow_block,&
            error=error)
       CALL get_mo_set(mos(1)%mo_set,nao=nao)
       CALL cp_fm_struct_create(fmstruct=ao_ao_fmstruct,&
            nrow_block=nrow_block,ncol_block=nrow_block,&
            nrow_global=nao, ncol_global=nao,&
            template_fmstruct=ao_mo_fmstruct, error=error)


       IF (.NOT.ASSOCIATED(scf_env%scf_work1)) THEN
          IF (scf_env%method/=ot_method_nr) THEN
             ALLOCATE(scf_env%scf_work1(dft_control%nspins), stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(scf_env%scf_work1)) THEN
             DO ispin=1,SIZE(scf_env%scf_work1)
                NULLIFY(scf_env%scf_work1(ispin)%matrix)
                CALL cp_fm_create(scf_env%scf_work1(ispin)%matrix,&
                     matrix_struct=ao_ao_fmstruct,&
                     name="SCF"//TRIM(ADJUSTL(cp_to_string(scf_env%id_nr)))//&
                     "WORK_MATRIX-1-"//TRIM(ADJUSTL(cp_to_string(ispin))),&
                     error=error)
             ENDDO
          ENDIF
       END IF
       IF ((.NOT.ASSOCIATED(scf_env%ortho)).AND.&
           ((scf_env%method /= ot_method_nr).AND.&
            (scf_env%method /= special_diag_method_nr))) THEN
          CALL cp_fm_create(scf_env%ortho,&
               matrix_struct=ao_ao_fmstruct,&
               name="SCF"//TRIM(ADJUSTL(cp_to_string(scf_env%id_nr)))//&
               "ORTHO_MATRIX",&
               error=error)
       END IF
       IF ((.NOT.ASSOCIATED(scf_env%scf_work2)).AND.&
           (scf_env%method /= ot_method_nr)) THEN
          CALL cp_fm_create(scf_env%scf_work2,&
               matrix_struct=ao_ao_fmstruct,&
               name="SCF"//TRIM(ADJUSTL(cp_to_string(scf_env%id_nr)))//&
               "WORK_MATRIX-2",&
               error=error)
       END IF
       CALL cp_fm_struct_release(ao_ao_fmstruct,error=error)

    END IF

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
       IF (.NOT.ASSOCIATED(scf_env%p_mix_new)) THEN
          CALL sm_pools_create_matrix_vect(S_sm_pools,scf_env%p_mix_new,&
               name="SCF"//TRIM(ADJUSTL(cp_to_string(scf_env%id_nr)))//&
               "DENSITY",&
               error=error)
       END IF
    END IF

!   *** allocate the ks env **
    IF (.NOT.ASSOCIATED(ks_env)) THEN
       CALL qs_ks_create(ks_env,qs_env=qs_env)
       CALL set_qs_env(qs_env, ks_env=ks_env)
       CALL qs_ks_release(ks_env)
    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE scf_env_check_i_alloc
!***************************************************************************

  SUBROUTINE calculate_ortho_matrix(ortho,matrix_s,work,ndep,eps_eigval,&
       work_syevx)

!   Purpose: Calculate the orthogonalization matrix (S**(-1/2))

!   History: - Creation (01.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(cp_fm_type), POINTER                :: ortho
    TYPE(real_matrix_type), POINTER          :: matrix_s
    TYPE(cp_fm_type), POINTER                :: work
    INTEGER, INTENT(OUT)                     :: ndep
    REAL(KIND=dp), INTENT(in)                :: eps_eigval, work_syevx

    INTEGER                                  :: handle

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_ortho_matrix","I"," ",handle)

    CALL copy_sm_to_fm(matrix_s,ortho)
    CALL cp_fm_power(ortho,work,-0.5_dp,eps_eigval,ndep,work_syevx)
    CALL cp_fm_upper_to_full(ortho,work)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE calculate_ortho_matrix

! *****************************************************************************

  SUBROUTINE eigensolver(matrix_ks,mo_set,ortho,work,do_level_shift,level_shift,&
       use_cholesky, work_syevx, use_jacobi,jacobi_threshold,smear,&
       error)

!   Purpose: Diagonalise the Kohn-Sham matrix to get a new set of MO eigen-
!            vectors and MO eigenvalues.
!   ks will be modified

!   History: - Creation (01.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(cp_fm_type), POINTER                :: matrix_ks
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_fm_type), POINTER                :: ortho, work
    LOGICAL, INTENT(IN)                      :: do_level_shift
    REAL(KIND=dp), INTENT(IN)                :: level_shift
    LOGICAL, INTENT(INOUT)                   :: use_cholesky
    REAL(KIND=dp), INTENT(IN)                :: work_syevx
    LOGICAL, INTENT(IN)                      :: use_jacobi
    REAL(KIND=dp), INTENT(IN)                :: jacobi_threshold, smear
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    INTEGER                                  :: handle, homo, imo, lfomo, &
                                                nao, nelectron, nmo
    LOGICAL                                  :: mo_uocc
    REAL(KIND=dp)                            :: alpha, maxocc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mo_eigenvalues
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mo_occupation
    TYPE(cp_fm_type), POINTER                :: mo_coeff

!   ---------------------------------------------------------------------------

    CALL timeset("eigensolver","I"," ",handle)

!   *** Diagonalise the Kohn-Sham matrix ***

    CALL get_mo_set(mo_set=mo_set,&
         nao=nao,&
         nmo=nmo,&
         homo=homo, &
         nelectron=nelectron,&
         maxocc=maxocc,&
         eigenvalues=mo_eigenvalues,&
         occupation_numbers=mo_occupation,&
         mo_coeff=mo_coeff)

    IF (do_level_shift) THEN
       IF (use_cholesky) THEN
          CALL stop_program("eigensolver","level shift not implemented")
       ENDIF
    ENDIF

    IF (use_cholesky) THEN
       CALL cp_fm_cholesky_reduce(matrix_ks,ortho)

       IF (use_jacobi) THEN
          CALL cp_fm_syevd(matrix_ks,work,mo_eigenvalues)
          use_cholesky = .FALSE.
          CALL cp_fm_cholesky_restore(work,nmo,ortho,mo_coeff,"SOLVE")
          CALL cp_fm_to_fm(mo_coeff,ortho)
       ELSE IF (nmo==nao) THEN
          CALL cp_fm_syevd(matrix_ks,work,mo_eigenvalues)
          CALL cp_fm_cholesky_restore(work,nmo,ortho,mo_coeff,"SOLVE")
       ELSE
          CALL cp_fm_syevx(matrix_ks,work,mo_eigenvalues,nmo,work_syevx)
          CALL cp_fm_cholesky_restore(work,nmo,ortho,mo_coeff,"SOLVE")
       END IF

    ELSE IF (use_jacobi) THEN

       CALL cp_fm_symm("L","U",nao,homo,1.0_dp,matrix_ks,mo_coeff,0.0_dp,work)
       CALL cp_fm_gemm("T","N",homo,nao-homo,nao,1.0_dp,work,mo_coeff,&
            0.0_dp,matrix_ks,b_first_col=homo+1,c_first_col=homo+1)

!     epsi =1e-11
!     Anz_max=4
!     N=matrix_ks%matrix_struct%nrow_global

!     *** Scale the elements with the inverse ***
!     *** difference of the old eigenvalues   ***

!     CALL blacs_block_jacobi_scaled(matrix_ks,work,mo_eigenvalues,&
!                                    jacobi_threshold,homo+1)

!     *** Klassisches Jacobi-Block-Verfahren ***
!     *** Rotiere Elemente nur, wenn sie groesser als thresh sind ***

       CALL cp_fm_block_jacobi_classic(matrix_ks,mo_coeff,mo_eigenvalues,&
            jacobi_threshold,homo+1)

    ELSE ! full S^-1 has been computed

       CALL cp_fm_symm("L","U",nao,nao,1.0_dp,matrix_ks,ortho,0.0_dp,work)
       CALL cp_fm_gemm("T","N",nao,nao,nao,1.0_dp,ortho,work,0.0_dp,matrix_ks)

       IF (do_level_shift) THEN
          DO imo=homo+1,nmo
             ! now unnecessary broadcast of element
             CALL cp_fm_get_element(matrix_ks,imo,imo,alpha)
             alpha=alpha+level_shift
             CALL cp_fm_set_element(matrix_ks,imo,imo,alpha)
          END DO
       END IF

       CALL cp_fm_syevx(matrix_ks,work,mo_eigenvalues,nmo,work_syevx)

       CALL cp_fm_gemm("N","N",nao,nmo,nao,1.0_dp,ortho,work,0.0_dp,&
            mo_coeff)

       IF (do_level_shift) THEN

          CALL correct_mo_eigenvalues(mo_set,level_shift)

!     *** Use last MO set as orthogonalization matrix ***

          CALL cp_fm_to_fm(mo_coeff,ortho)

       END IF

    END IF

    mo_uocc = .TRUE.
    IF (smear /= 0.0_dp) THEN
      CALL set_mo_occupation(mo_occupation,mo_eigenvalues,&
                             homo,lfomo,maxocc,nelectron,&
                             mo_uocc,&
                             smear=smear)
      CALL set_mo_set(mo_set=mo_set,&
                      homo=homo,lfomo=lfomo,&
                      uniform_occupation=mo_uocc)
    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE eigensolver

! *****************************************************************************
  SUBROUTINE simple_eigensolver(matrix_ks,mo_set,do_level_shift,level_shift,&
             work_syevx,smear,error)

    TYPE(cp_fm_type), POINTER                :: matrix_ks
    TYPE(mo_set_type), POINTER               :: mo_set
    LOGICAL, INTENT(IN)                      :: do_level_shift
    REAL(KIND=dp), INTENT(IN)                :: level_shift
    REAL(KIND=dp), INTENT(IN)                :: work_syevx
    REAL(KIND=dp), INTENT(IN)                :: smear
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    INTEGER                                  :: handle, homo, imo, lfomo, &
                                                nao, nelectron, nmo
    LOGICAL                                  :: mo_uocc
    REAL(KIND=dp)                            :: alpha, maxocc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mo_eigenvalues
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mo_occupation
    TYPE(cp_fm_type), POINTER                :: mo_coeff

!   ---------------------------------------------------------------------------

    CALL timeset("eigensolver","I"," ",handle)

    NULLIFY(mo_eigenvalues,mo_occupation,mo_coeff)

    CALL get_mo_set(mo_set=mo_set,&
         nao=nao,&
         nmo=nmo,&
         homo=homo, &
         nelectron=nelectron,&
         maxocc=maxocc,&
         eigenvalues=mo_eigenvalues,&
         occupation_numbers=mo_occupation,&
         mo_coeff=mo_coeff)

    IF (do_level_shift) THEN
      DO imo=homo+1,nmo
         ! now unnecessary broadcast of element
         CALL cp_fm_get_element(matrix_ks,imo,imo,alpha)
         alpha=alpha+level_shift
         CALL cp_fm_set_element(matrix_ks,imo,imo,alpha)
      END DO
    END IF

    CALL cp_fm_syevx(matrix_ks,mo_coeff,mo_eigenvalues,nmo,work_syevx)

    IF (do_level_shift) THEN
      CALL correct_mo_eigenvalues(mo_set,level_shift)
    END IF

    mo_uocc = .TRUE.
    IF (smear /= 0.0_dp) THEN
      CALL set_mo_occupation(mo_occupation,mo_eigenvalues,&
                             homo,lfomo,maxocc,nelectron,&
                             mo_uocc, smear=smear)
      CALL set_mo_set(mo_set=mo_set,&
                      homo=homo,lfomo=lfomo,&
                      uniform_occupation=mo_uocc)
    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE simple_eigensolver

! *****************************************************************************

!!****f* qs_scf/cp_sm_mix [1.0] *
!!
!!   NAME
!!     cp_sm_mix
!!
!!   FUNCTION
!!     Perform a mixing of the given matrixes into the first matrix
!!     m1 = m2 + p_mix (m1-m2)
!!
!!   NOTES
!!     if you what to store the result in m2 swap m1 and m2 an use 
!!     (1-pmix) as pmix
!!     para_env should be removed (embedded in matrix)
!!
!!   ARGUMENTS
!!     - m1: first (new) matrix, is modified
!!     - m2: the second (old) matrix
!!     - p_mix:how much m1 is conserved (0: none, 1: all)
!!     - delta: maximum norm of m1-m2
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 rewamped [fawzi]
!!
!!*** *********************************************************************
  SUBROUTINE cp_sm_mix(m1,m2,p_mix,delta,para_env,error)

    TYPE(real_matrix_type), POINTER          :: m1, m2
    REAL(KIND=dp), INTENT(IN)                :: p_mix
    REAL(KIND=dp), INTENT(OUT)               :: delta
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = "cp_sm_mix", &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iblock_col_m1, iblock_col_m2, &
                                                iblock_row, j, nblock_row
    REAL(KIND=dp)                            :: r_mix
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: p_new_block, p_old_block
    TYPE(real_block_node_type), POINTER      :: block_node_m1, block_node_m2
    LOGICAL :: failure

!   ---------------------------------------------------------------------------

    CALL timeset(routineN,"I"," ",handle)
    failure = .FALSE.
    delta = 0.0_dp

    r_mix = 1.0_dp - p_mix

    CALL get_matrix_info(matrix=m2,nblock_row=nblock_row)

    DO iblock_row=1,nblock_row

       block_node_m1 => first_block_node(matrix=m1,&
            block_row=iblock_row)
       block_node_m2 => first_block_node(matrix=m2,&
            block_row=iblock_row)

       DO WHILE (ASSOCIATED(block_node_m1))

          CALL get_block_node(block_node=block_node_m1,&
               block_col=iblock_col_m1,&
               BLOCK=p_new_block)

          CALL get_block_node(block_node=block_node_m2,&
               block_col=iblock_col_m2,&
               BLOCK=p_old_block)
          CPPostcondition(iblock_col_m1==iblock_col_m2,cp_failure_level,routineP,error,failure)

          DO j=1,SIZE(p_new_block,2)
             DO i=1,SIZE(p_new_block,1)
                p_new_block(i,j) = p_new_block(i,j) - p_old_block(i,j)
                delta = MAX(delta,ABS(p_new_block(i,j)))
                p_new_block(i,j) = p_old_block(i,j) + p_mix*p_new_block(i,j)
             END DO
          END DO

          block_node_m1 => next_block_node(block_node_m1)
          block_node_m2 => next_block_node(block_node_m2)

       END DO

    END DO

    CALL mp_max(delta,para_env%group)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE cp_sm_mix

! *****************************************************************************
  SUBROUTINE calculate_first_density_matrix(scf_env,qs_env,globenv,error)

    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = "calculate_first_density_matrix", &
      routineP = moduleN//":"//routineN

    INTEGER :: atom_a, group, handle, homo, iatom, ikind, iset, isgf, isgfa, &
      ishell, ispin, istat, la, maxl, maxll, nao, natom, ncount, nelectron, &
      nmo, nset, nspin, qs_env_id
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_sgf
    INTEGER, DIMENSION(:), POINTER           :: atom_list, elec_conf, nshell
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, l, last_sgfa
    REAL(KIND=dp)                            :: maxocc, nelec, paa, scale, &
                                                trps1, trps2, total_rho
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: econf, pdiag, sdiag
    TYPE(all_potential_type), POINTER        :: all_potential
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: work1
    TYPE(cp_fm_type), POINTER                :: mo_coeff, ortho, sv, &
                                                work2
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gth_potential_type), POINTER        :: gth_potential
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: p_rmpv, s_sparse
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: h_core_sparse
    TYPE(scf_control_type), POINTER          :: scf_control
    LOGICAL                                  :: ortho_basis, &
                                                id_equal, &
                                                failure

!   ---------------------------------------------------------------------------

    NULLIFY(all_potential,atomic_kind, mo_coeff, sv,&
         gth_potential, orb_basis_set, atomic_kind_set, particle_set,&
         ortho,work2,work1,mo_array,s_sparse,p_rmpv,scf_control, &
         dft_control,h_core_sparse)

    CALL timeset("calculate_first_density_matrix","I"," ",handle)

    CALL get_qs_env(qs_env,atomic_kind_set=atomic_kind_set,&
         particle_set=particle_set,mos=mo_array, matrix_s=s_sparse,&
         matrix_h=h_core_sparse,&
         scf_control=scf_control, id_nr=qs_env_id, dft_control=dft_control,&
         error=error)

    nspin=dft_control%nspins

    p_rmpv => qs_env%rho%rho_ao
    work1 => scf_env%scf_work1
    work2 => scf_env%scf_work2
    ortho => scf_env%ortho

    ortho_basis = .FALSE.
    IF ( dft_control%qs_control%semi_empirical ) THEN
       IF ( dft_control%qs_control%se_control%orthogonal_basis ) &
          ortho_basis = .TRUE.
    END IF

    IF (scf_control%use_ot.AND.&
        (.NOT.((scf_control%density_guess == "RANDOM").OR.&
               (scf_control%density_guess == "ATOMIC").OR.&
               ((scf_control%density_guess == "RESTART").AND.&
                (scf_control%level_shift == 0.0_dp))))) THEN
       CALL stop_program("calculate_first_density_matrix",&
            "OT needs GUESS ATOMIC / RESTART : other options NYI")
    END IF

    IF (scf_control%density_guess == "RESTART") THEN

       CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set)
       CALL read_mo_set(mo_array,atomic_kind_set,particle_set,globenv, id_nr=1)

       DO ispin=1,nspin
          IF (scf_control%level_shift /= 0.0_dp) THEN
             CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
                  mo_coeff=mo_coeff)
             CALL cp_fm_to_fm(mo_coeff,ortho)
          END IF

          ! make all nmo vectors present orthonormal
          CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
               mo_coeff=mo_coeff, nmo=nmo, homo=homo)
          IF(ortho_basis) THEN
            CALL make_basis(mo_coeff,nmo)
          ELSE
            ! ortho so that one can restart for different positions (basis sets?)
            CALL cp_fm_create(sv,mo_coeff%matrix_struct,"SV")
            CALL cp_sm_fm_multiply(s_sparse(1)%matrix,mo_coeff,sv,nmo)
            CALL make_basis(mo_coeff,nmo, &
                 ortho=sv,otype="SV")
            CALL cp_fm_release(sv)
          ENDIF
          ! only alpha spin is kept for restricted
          IF (dft_control%restricted) EXIT
       ENDDO
       IF (dft_control%restricted) CALL mo_set_restrict(mo_array)
       DO ispin=1,nspin
          CALL calculate_density_matrix(mo_array(ispin)%mo_set,&
               p_rmpv(ispin)%matrix,error=error)
       ENDDO

    ELSE IF (scf_control%density_guess == "RANDOM") THEN

       DO ispin=1,nspin
          CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
               mo_coeff=mo_coeff, nmo=nmo)
          CALL cp_fm_init_random(mo_coeff,nmo)
          IF(ortho_basis) THEN
            CALL make_basis(mo_coeff,nmo)
          ELSE
            ! ortho so that one can restart for different positions (basis sets?)
            CALL cp_fm_create(sv,mo_coeff%matrix_struct,"SV")
            CALL cp_sm_fm_multiply(s_sparse(1)%matrix,mo_coeff,sv,nmo)
            CALL make_basis(mo_coeff,nmo, &
                 ortho=sv,otype="SV")
            CALL cp_fm_release(sv)
          ENDIF
          ! only alpha spin is kept for restricted
          IF (dft_control%restricted) EXIT
       ENDDO
       IF (dft_control%restricted) CALL mo_set_restrict(mo_array)
       DO ispin=1,nspin
          CALL calculate_density_matrix(mo_array(ispin)%mo_set,&
               p_rmpv(ispin)%matrix,error=error)
       ENDDO

    ELSE IF (scf_control%density_guess == "CORE") THEN
 
       IF (dft_control%restricted) CALL stop_program("calculate_first_density_matrix","Option not yet implemented")

       ! work1(1) should contain the core hamiltonian
       ! we need a copy the core matrix for every spin (the dumbest solution,
       ! we should copy the mos)..
       
       CALL copy_sm_to_fm(h_core_sparse(1)%matrix,work1(1)%matrix)

       IF (dft_control%nspins.eq.2) THEN
          IF (.NOT.scf_control%use_ot) THEN
             CALL cp_fm_to_fm(work1(1)%matrix,&
                  work1(2)%matrix)
          ENDIF
       ENDIF
       DO ispin=1,nspin
          IF(ortho_basis) THEN
            CALL simple_eigensolver(work1(ispin)%matrix,mo_array(ispin)%mo_set,&
                 .FALSE.,0.0_dp,work_syevx=scf_control%work_syevx,&
                 smear=0.0_dp,error=error)
          ELSE
            CALL eigensolver(work1(ispin)%matrix,mo_array(ispin)%mo_set,ortho,work2, &
               .FALSE.,0.0_dp,&
               use_cholesky=scf_control%use_cholesky, &
               work_syevx=scf_control%work_syevx,&
               use_jacobi=.FALSE.,&
               jacobi_threshold=scf_control%jacobi_threshold,&
               smear=0.0_dp,&
               error=error)
          END IF
          IF (scf_control%level_shift /= 0.0_dp) THEN
             CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
                  mo_coeff=mo_coeff)
             CALL cp_fm_to_fm(mo_coeff,ortho)
          END IF
          CALL calculate_density_matrix(mo_array(ispin)%mo_set,&
               p_rmpv(ispin)%matrix,error=error)
       ENDDO

    ELSE IF (scf_control%density_guess == "ATOMIC") THEN

      IF (dft_control%restricted) CALL stop_program("calculate_first_density_matrix","Option not yet implemented")

      group = qs_env%para_env%group

      natom = SIZE(particle_set)
      ALLOCATE (first_sgf(natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineP,"first_sgf",natom*int_size)
      CALL get_particle_set(particle_set=particle_set,&
                            first_sgf=first_sgf)

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,maxlgto=maxl)
      ALLOCATE (econf(0:maxl),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineP,"econf",(maxl + 1)*int_size)

      CALL get_mo_set(mo_array(1)%mo_set,nao=nao)
      ALLOCATE (pdiag(nao),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineP,"pdiag",nao*dp_size)
      pdiag(:) = 0.0_dp

      ALLOCATE (sdiag(nao),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineP,"sdiag",nao*dp_size)
      IF (ortho_basis) THEN
        sdiag=1._dp
      ELSE
        CALL get_matrix_diagonal(s_sparse(1)%matrix,sdiag)
        CALL mp_sum(sdiag,group)
      END IF

      DO ispin=1,nspin

        CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
                        maxocc=maxocc,&
                        nelectron=nelectron)

        ncount = 0
        trps1 = 0.0_dp
        trps2 = 0.0_dp
        pdiag(:) = 0.0_dp

        IF (nelectron /= 0) THEN 

           DO ikind=1,SIZE(atomic_kind_set)

              atomic_kind => atomic_kind_set(ikind)
   
              CALL get_atomic_kind(atomic_kind=atomic_kind,&
                                   natom=natom,&
                                   atom_list=atom_list,&
                                   all_potential=all_potential,&
                                   gth_potential=gth_potential,&
                                   orb_basis_set=orb_basis_set)

              IF (ASSOCIATED(all_potential)) THEN
                 CALL get_potential(potential=all_potential,elec_conf=elec_conf)
              ELSE IF (ASSOCIATED(gth_potential)) THEN
                 CALL get_potential(potential=gth_potential,elec_conf=elec_conf)
              ELSE
                 CYCLE
              END IF

              maxll = MIN(SIZE(elec_conf) - 1,maxl)
              econf(:) = 0.0_dp
              econf(0:maxll) = 0.5_dp*maxocc*REAL(elec_conf(0:maxll),dp)

              CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                                     nset=nset,&
                                     nshell=nshell,&
                                     l=l,&
                                     first_sgf=first_sgfa,&
                                     last_sgf=last_sgfa)

              DO iset=1,nset
                 DO ishell=1,nshell(iset)
                    la = l(ishell,iset)
                    nelec = maxocc*REAL(2*la + 1,dp)
                    IF (econf(la) > 0.0_dp) THEN
                       IF (econf(la) >= nelec) THEN
                          paa = maxocc
                          econf(la) = econf(la) - nelec
                       ELSE
                          paa = maxocc*econf(la)/nelec
                          econf(la) = 0.0_dp
                          ncount = ncount + NINT(nelec/maxocc)
                       END IF
                       DO isgfa=first_sgfa(ishell,iset),last_sgfa(ishell,iset)
                          DO iatom=1,natom
                             atom_a = atom_list(iatom)
                             isgf = first_sgf(atom_a) + isgfa - 1
                             pdiag(isgf) = paa
                             IF (paa == maxocc) THEN
                                trps1 = trps1 + paa*sdiag(isgf)
                             ELSE
                                trps2 = trps2 + paa*sdiag(isgf)
                             END IF
                          END DO
                       END DO
                    END IF
                 END DO
              END DO

           END DO

           IF (trps2 == 0.0_dp) THEN
              DO isgf=1,nao
                 IF (sdiag(isgf) > 0.0_dp) pdiag(isgf) = pdiag(isgf)/sdiag(isgf)
              END DO
           ELSE
              scale = (REAL(nelectron,dp) - trps1)/trps2
              DO isgf=1,nao 
                 IF (pdiag(isgf) < maxocc) pdiag(isgf) = scale*pdiag(isgf)
              END DO
           END IF
        END IF

        CALL set_matrix_diagonal(p_rmpv(ispin)%matrix,pdiag)

        ! do we require an initial wavefunction ?
        ! assume p is a projector (it is not ..) and generate an initial  
        ! wavefunction and the corresponding dm
        IF (scf_control%use_ot) THEN
          CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
               mo_coeff=mo_coeff, nmo=nmo)
          CALL cp_fm_init_random(mo_coeff,nmo)
          CALL cp_fm_create(sv,mo_coeff%matrix_struct,"SV")

          ! multiply times PS
          IF (ortho_basis) THEN
             CALL cp_fm_to_fm(mo_coeff,sv)
          ELSE
             CALL cp_sm_fm_multiply(s_sparse(1)%matrix,mo_coeff,sv,nmo)
          ENDIF
          ! here we could easily multiply with the diag that we actually have replicated already 
          CALL cp_sm_fm_multiply(p_rmpv(ispin)%matrix,sv,mo_coeff,nmo)

          ! and ortho the result
          IF (ortho_basis) THEN
             CALL make_basis(mo_coeff,nmo)
          ELSE
             CALL cp_sm_fm_multiply(s_sparse(1)%matrix,mo_coeff,sv,nmo)
             CALL make_basis(mo_coeff,nmo, ortho=sv,otype="SV")
          ENDIF
          CALL cp_fm_release(sv)

          CALL calculate_density_matrix(mo_array(ispin)%mo_set,&
               p_rmpv(ispin)%matrix,error=error)
        ENDIF

      END DO

      DEALLOCATE (econf,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineP,"econf")

      DEALLOCATE (first_sgf,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineP,"first_sgf")

      DEALLOCATE (pdiag,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineP,"pdiag")

      DEALLOCATE (sdiag,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineP,"sdiag")

    ELSE IF (scf_control%density_guess == "DENSITIES") THEN

      ! Collocation of the density into the PW-grid
      CALL collocate_atomic_charge_density(total_rho=total_rho, qs_env=qs_env, error=error)

      ! do some assertions here on these matrices having the same structure, 
      ! as is currently required
      DO ispin=1,SIZE(qs_env%matrix_s)
        id_equal=(qs_env%matrix_s(ispin)%matrix%sparsity_id==qs_env%matrix_s(1)%matrix%sparsity_id)
        CPPrecondition(id_equal,cp_failure_level,routineP,error,failure)
      ENDDO
      ! matrix_p = rho%rho_ao
      DO ispin=1,SIZE(qs_env%rho%rho_ao)
        id_equal=(qs_env%rho%rho_ao(ispin)%matrix%sparsity_id==qs_env%matrix_s(1)%matrix%sparsity_id)
        CPPrecondition(id_equal,cp_failure_level,routineP,error,failure)
      ENDDO
      DO ispin=1,SIZE(qs_env%matrix_h)
        id_equal=(qs_env%matrix_h(ispin)%matrix%sparsity_id==qs_env%matrix_s(1)%matrix%sparsity_id)
        CPPrecondition(id_equal,cp_failure_level,routineP,error,failure)
      ENDDO
      DO ispin=1,SIZE(qs_env%matrix_ks)
        id_equal=(qs_env%matrix_ks(ispin)%matrix%sparsity_id==qs_env%matrix_s(1)%matrix%sparsity_id)
        CPPrecondition(id_equal,cp_failure_level,routineP,error,failure)
      ENDDO
      ! end sparsity check

      qs_env%scf_env%iter_count = 1

      ! qs_ks_build_kohn_sham_matrix (Integration)
      !CALL qs_ks_build_kohn_sham_matrix(ks_env=qs_env%ks_env, qs_env=qs_env, &
      !                                  ks_matrix=qs_env%matrix_ks, rho=qs_env%rho, &
      !                                  energy=qs_env%energy, calculate_forces=.TRUE., &
      !                                  just_energy=.FALSE., kg_gpw=.FALSE., error=error)
      CALL qs_ks_did_change(ks_env=qs_env%ks_env, rho_changed=.TRUE., error=error)

      CALL qs_ks_update_qs_env(ks_env=qs_env%ks_env, qs_env=qs_env, &
                               calculate_forces=.TRUE., just_energy=.FALSE., error=error)

      ! diagonalization
      DO ispin = 1,qs_env%dft_control%nspins
        CALL copy_sm_to_fm(qs_env%matrix_ks(ispin)%matrix, qs_env%scf_env%scf_work1(ispin)%matrix)
      END DO

      qs_env%scf_env%iter_method = "Mixing/Diag"
      qs_env%scf_env%iter_delta = 0.0_dp

      DO ispin = 1,qs_env%dft_control%nspins
        CALL eigensolver(matrix_ks=qs_env%scf_env%scf_work1(ispin)%matrix, &
                         mo_set=qs_env%mos(ispin)%mo_set, ortho=qs_env%scf_env%ortho, &
                         work=qs_env%scf_env%scf_work2, do_level_shift=.FALSE., &
                         level_shift=scf_control%level_shift, &
                         use_cholesky=scf_control%use_cholesky, work_syevx=scf_control%work_syevx, &
                         use_jacobi=.FALSE., jacobi_threshold=scf_control%jacobi_threshold, &
                         smear=scf_control%smear, error=error)

        ! calculate_density_matrix
        CALL calculate_density_matrix(qs_env%mos(ispin)%mo_set, qs_env%rho%rho_ao(ispin)%matrix, &
                                      error=error)
      END DO

    ELSE

      CALL stop_program(routineP,&
                        "An invalid keyword for the initial density "//&
                        "guess was specified")

    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE calculate_first_density_matrix

!!****f* qs_scf/scf_env_initial_rho_setup [1.0] *
!!
!!   NAME
!!     scf_env_initial_rho_setup
!!
!!   FUNCTION
!!     Initializes rho and the mos, so that an scf cycle can start
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - scf_env: the scf env in which to do the scf
!!     - qs_env: the qs env the scf_env lives in
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE scf_env_initial_rho_setup(scf_env, qs_env, globenv, error)
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), &
      POINTER                                :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scf_env_initial_rho_setup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: extrapolation_method_nr, &
                                                handle, ispin, nmo
    LOGICAL                                  :: failure, gapw, gapw_xc, ionode, orthogonal_wf
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,"I"," ",handle)
    failure=.FALSE.
    NULLIFY(mo_coeff)
    gapw = qs_env%dft_control%qs_control%gapw
    gapw_xc = qs_env%dft_control%qs_control%gapw_xc
    logger => cp_error_get_logger(error)
    ionode = logger%para_env%source==logger%para_env%mepos
    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)

    IF (.NOT.failure) THEN
       extrapolation_method_nr=wfi_use_guess_method_nr
       IF (ASSOCIATED(qs_env%wf_history)) THEN
          CALL wfi_extrapolate(globenv,qs_env%wf_history, &
               qs_env=qs_env, dt=1.0_dp, &
               extrapolation_method_nr=extrapolation_method_nr,&
               orthogonal_wf=orthogonal_wf, gapw=gapw, gapw_xc=gapw_xc, error=error)
          ! wfi_use_guess_method_nr the wavefunctions are not yet initialized
          IF ((.NOT.orthogonal_wf).AND.&
              (scf_env%method == ot_method_nr).AND.&
              (.NOT.(extrapolation_method_nr == wfi_use_guess_method_nr))) THEN
             DO ispin=1,SIZE(qs_env%mos)
                CALL get_mo_set(qs_env%mos(ispin)%mo_set, &
                     mo_coeff=mo_coeff, nmo=nmo)
                CALL qs_env_reorthogonalize_vectors(qs_env, &
                     v_matrix=mo_coeff, ispin=ispin, n_col=nmo,&
                     error=error)
             END DO
          END IF
       END IF
       IF (ionode.AND.logger%print_keys%scf) &
       WRITE (UNIT=cp_logger_get_default_unit_nr(logger),FMT="(/,T3,A)")&
         "Extrapolation method: "//&
         TRIM(wfi_get_method_label(extrapolation_method_nr,error=error))
       IF (extrapolation_method_nr==wfi_use_guess_method_nr) THEN
          CALL calculate_first_density_matrix(scf_env=scf_env,qs_env=qs_env,&
               globenv=globenv, error=error)
          IF (.NOT.(qs_env%scf_control%density_guess=="DENSITIES")) THEN
            CALL qs_rho_update_rho(qs_env%rho,qs_env=qs_env, gapw=gapw, error=error)
            IF(gapw) THEN
               CALL calculate_rho_atom_coeff(qs_env,qs_env%rho%rho_ao)
            ENDIF
            IF(gapw_xc) THEN
               CALL qs_rho_update_rho(qs_env%rho_xc,qs_env=qs_env, gapw=gapw_xc, error=error)
               CALL calculate_rho_atom_coeff(qs_env,qs_env%rho%rho_ao)
            END IF
            CALL qs_ks_did_change(qs_env%ks_env,rho_changed=.TRUE.,&
                 error=error)
          END IF
       END IF
    END IF

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE scf_env_initial_rho_setup
!***************************************************************************

!!****f* qs_scf/scf_env_cleanup [1.0] *
!!
!!   NAME
!!     scf_env_cleanup
!!
!!   FUNCTION
!!     perform cleanup operations (like releasing temporary storage)
!!     at the end of the scf
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE scf_env_cleanup(scf_env,qs_env,error)
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scf_env_cleanup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, stat
    LOGICAL                                  :: failure
    TYPE(cp_sm_pool_p_type), DIMENSION(:), &
      POINTER                                :: S_sm_pools

    CALL timeset(routineN,"I"," ",handle)

    failure=.FALSE.
    NULLIFY(S_sm_pools)

    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CALL mpools_get(qs_env%mpools,S_sm_pools=S_sm_pools,error=error)

!   *** Release SCF work storage ***

       IF (ASSOCIATED(scf_env%scf_work1)) THEN
          DO ispin=1,SIZE(scf_env%scf_work1)
             CALL cp_fm_release(scf_env%scf_work1(ispin)%matrix)
          ENDDO
          DEALLOCATE(scf_env%scf_work1)
       ENDIF
       IF (ASSOCIATED(scf_env%scf_work2)) CALL cp_fm_release(scf_env%scf_work2)
       IF (ASSOCIATED(scf_env%ortho)) CALL cp_fm_release(scf_env%ortho)

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
       CASE DEFAULT
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP,"unknown scf method method:"//&
               cp_to_string(scf_env%method),error,failure)
       END SELECT

    END IF
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE scf_env_cleanup
!***************************************************************************

!!****f* qs_scf/scf_env_did_change [1.0] *
!!
!!   NAME
!!     scf_env_did_change
!!
!!   SYNOPSIS
!!     Subroutine scf_env_did_change(scf_env, qs_env, s_mstruct_changed,&
!!         error)
!!       Type(qs_scf_env_type), Pointer:: scf_env
!!       Type(qs_environment_type), Pointer:: qs_env
!!       Logical, Intent (IN):: s_mstruct_changed
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine scf_env_did_change
!!
!!   FUNCTION
!!     function to be called to inform the scf_env about changes
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - scf_env: the scf env to inform
!!     - qs_env: the qs_env in which the scf_env lives
!!     - s_mstruct_changed: treu if the structure of s changed
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     03.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE scf_env_did_change(scf_env,error)
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scf_env_did_change', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     IF (ASSOCIATED(scf_env%p_mix_new)) THEN
        CALL deallocate_matrix_set(scf_env%p_mix_new)
     END IF
  END IF
END SUBROUTINE scf_env_did_change

!***************************************************************************

END MODULE qs_scf
