!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2002  CP2K developers group
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
!!
!!   SOURCE
!******************************************************************************

MODULE qs_scf
! *****************************************************************************
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE coefficient_types,               ONLY: coeff_allocate,&
                                             coeff_copy,&
                                             coeff_reduce,&
                                             coeff_sumup,&
                                             coeff_transform_space,&
                                             coeff_type,&
                                             coeff_zero
  USE cp_array_utils,                  ONLY: cp_1d_r_output
  USE cp_block_matrix,                 ONLY: cp_block_matrix_p_type,&
                                             cp_block_matrix_type,&
                                             cp_sp_multiplies_blacs
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                             cp_fm_pool_type,&
                                             fm_pool_get_mstruct
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_get,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create
  USE cp_lbfgs_optimizer_gerd_f,       ONLY: cp_lbfgs_opt_gerd_f_type,&
                                             cp_opt_gerd_f_dealloc_ref,&
                                             cp_opt_gerd_f_init,&
                                             cp_opt_gerd_f_next
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_log,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_matrix_utils,                 ONLY: cp_sm_output
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_sm_pool_types,                ONLY: cp_sm_pool_p_type,&
                                             cp_sm_pool_type,&
                                             sm_pool_create_matrix,&
                                             sm_pool_give_back_matrix,&
                                             sm_pools_create_matrix_vect,&
                                             sm_pools_flush_cache,&
                                             sm_pools_give_back_matrix_vect
  USE cube_utils,                      ONLY: cube_info_type,&
                                             destroy_cube_info,&
                                             init_cube_info
  USE dft_types,                       ONLY: dft_control_type
  USE external_potential_types,        ONLY: all_potential_type,&
                                             get_potential,&
                                             gth_potential_type
  USE global_types,                    ONLY: DEBUG,&
                                             global_environment_type
  USE hartree,                         ONLY: calculate_hartree
  USE kinds,                           ONLY: int_size,&
                                             wp => dp
  USE l_utils,                         ONLY: destroy_l_info,&
                                             init_l_info,&
                                             l_info_type
  USE machine,                         ONLY: m_flush,&
                                             m_walltime
  USE mathconstants,                   ONLY: fourpi,&
                                             twopi
  USE message_passing,                 ONLY: mp_allgather,&
                                             mp_max,&
                                             mp_range,&
                                             mp_sum,&
                                             mp_sync
  USE mulliken,                        ONLY: mulliken_charges
  USE pao_glob_angles_methods,         ONLY: pao_g_ang_get,&
                                             pao_g_ang_set
  USE pao_l_angles_methods,            ONLY: pao_guarantee_owned_NUi
  USE pao_obj_function,                ONLY: pao_gerd_f_dealloc_ref,&
                                             pao_gerd_f_init,&
                                             pao_gerd_functional_type
  USE pao_proj_methods,                ONLY: sm_matrix_transf_p_to_m,&
                                             sm_matrix_transf_to_f
  USE pao_qs_env_methods,              ONLY: pao_qs_env_did_change
  USE pao_types,                       ONLY: pao_env_get,&
                                             pao_env_type,&
                                             pao_glob_angles_type
  USE particle_types,                  ONLY: particle_type
  USE preconditioner,                  ONLY: destroy_preconditioner,&
                                             init_preconditioner,&
                                             make_preconditioner,&
                                             preconditioner_type
  USE pw_env_types,                    ONLY: pw_env_flush_cache,&
                                             pw_env_get,&
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
  USE qs_arpack,                       ONLY: arpack_diag_sy
  USE qs_arpack_full,                  ONLY: op_symv_type,&
                                             setup_op_symv
  USE qs_blacs,                        ONLY: blacs_add_to_element,&
                                             blacs_block_jacobi_classic,&
                                             blacs_get_element,&
                                             blacs_make_basis,&
                                             blacs_syevd,&
                                             blacs_syevx,&
                                             copy_blacs_to_blacs_matrix,&
                                             copy_full_matrix_to_real_matrix,&
                                             copy_real_matrix_to_full_matrix,&
                                             cp_fm_cholesky_decompose,&
                                             cp_fm_cholesky_reduce,&
                                             cp_fm_cholesky_restore,&
                                             cp_fm_create2,&
                                             cp_fm_gemm,&
                                             cp_fm_get_info,&
                                             cp_fm_init_random,&
                                             cp_fm_release,&
                                             cp_fm_scale_and_add,&
                                             cp_fm_set_all,&
                                             cp_fm_set_element,&
                                             cp_fm_symm,&
                                             cp_fm_trace,&
                                             cp_full_matrix_p_type,&
                                             cp_full_matrix_type,&
                                             cp_sm_fm_multiply,&
                                             power_blacs_matrix,&
                                             replicate_blacs_matrix,&
                                             symmetrise_blacs_matrix,&
                                             w_blacs_matrix => write_blacs_matrix
  USE qs_charges_types,                ONLY: qs_charges_type
  USE qs_collocate_density,            ONLY: calculate_rho_core,&
                                             calculate_total_abs_rho,&
                                             calculate_total_rho,&
                                             calculate_wavefunction
  USE qs_core_energies,                ONLY: calculate_ecore,&
                                             calculate_ecore_overlap,&
                                             calculate_ecore_self
  USE qs_core_hamiltonian,             ONLY: build_core_hamiltonian_matrix
  USE qs_diis,                         ONLY: qs_diis_b_clear,&
                                             qs_diis_b_create,&
                                             qs_diis_b_step
  USE qs_diis_types,                   ONLY: qs_diis_b_release,&
                                             qs_diis_buffer_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_ks_methods,                   ONLY: qs_ks_create,&
                                             qs_ks_did_change,&
                                             qs_ks_update_qs_env
  USE qs_ks_types,                     ONLY: qs_ks_env_p_type,&
                                             qs_ks_env_type,&
                                             qs_ks_release,&
                                             qs_ks_retain
  USE qs_localization,                 ONLY: qs_loc_compute_set,&
                                             qs_loc_destroy,&
                                             qs_loc_env_type,&
                                             qs_loc_init
  USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                             calculate_density_matrix,&
                                             calculate_subspace_eigenvalues,&
                                             correct_mo_eigenvalues,&
                                             get_mo_set,&
                                             init_mo_set,&
                                             mo_set_p_type,&
                                             mo_set_type,&
                                             read_mo_set,&
                                             set_mo_occupation,&
                                             write_mo_set
  USE qs_ot,                           ONLY: qs_ot_new_preconditioner,&
                                             qs_ot_type
  USE qs_ot_eigensolver,               ONLY: ot_eigensolver
  USE qs_ot_scf,                       ONLY: ot_scf_destroy,&
                                             ot_scf_init,&
                                             ot_scf_mini,&
                                             ot_scf_read_input
  USE qs_overlap,                      ONLY: write_blacs_matrix,&
                                             write_sparse_matrix
  USE qs_parser,                       ONLY: close_file,&
                                             open_file
  USE qs_rho_methods,                  ONLY: qs_rho_update_rho
  USE qs_rho_types,                    ONLY: qs_rho_did_change,&
                                             qs_rho_get,&
                                             qs_rho_type
  USE qs_scf_types,                    ONLY: qs_scf_env_type,&
                                             scf_env_release,&
                                             scf_env_retain
  USE qs_wf_history_methods,           ONLY: wfi_extrapolate,&
                                             wfi_update
  USE qs_wf_history_types,             ONLY: qs_wf_history_type,&
                                             qs_wf_snapshot_p_type,&
                                             qs_wf_snapshot_type
  USE realspace_grid_types,            ONLY: rs_pw_to_cube
  USE scf_control_types,               ONLY: scf_control_type
  USE simulation_cell,                 ONLY: cell_type,&
                                             get_cell
  USE sparse_matrix_types,             ONLY: add_all_real_matrix_blocks,&
                                             add_matrices,&
                                             allocate_matrix,&
                                             allocate_matrix_set,&
                                             cp_sm_scale_and_add,&
                                             deallocate_matrix,&
                                             deallocate_matrix_set,&
                                             first_block_node,&
                                             get_block_node,&
                                             get_matrix_info,&
                                             next_block_node,&
                                             real_block_node_type,&
                                             real_matrix_p_type,&
                                             real_matrix_type,&
                                             replicate_matrix,&
                                             replicate_matrix_structure,&
                                             set_matrix,&
                                             symmetrise_diagonal_blocks,&
                                             transfer_matrix,&
                                             write_blocks_maxabsval
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop,&
                                             trace_debug
  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: moduleN = "qs_scf"
  INTEGER, PARAMETER :: pao_mehod_nr=0, diagonalisation_method_nr=1,&
       ot_method_nr=2
  INTEGER, PARAMETER :: no_mixing_nr=0, ao_mixing_nr=1

  LOGICAL, PARAMETER, PRIVATE :: debug_this_module=.FALSE.

! *** Global variables ***

  INTEGER, SAVE, PRIVATE :: last_scf_env_id=0

! *** Public subroutines ***

  PUBLIC :: scf, scf_env_create, scf_env_do_scf, &
       scf_env_did_change

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
    TYPE(qs_scf_env_type), POINTER :: scf_env
    TYPE(qs_environment_type), INTENT(in), TARGET :: qs_env
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='scf_env_create',&
         routineP=moduleN//':'//routineN
    TYPE(scf_control_type), POINTER :: scf_control
    TYPE(dft_control_type), POINTER   :: dft_control
    TYPE(pao_env_type), POINTER :: pao_env
    INTEGER :: stat

    failure=.FALSE.
    NULLIFY(scf_control,dft_control,pao_env)

    ALLOCATE(scf_env, stat=stat)
    CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    IF (.NOT. failure) THEN
       CALL get_qs_env(qs_env=qs_env,&
            scf_control=scf_control,&
            dft_control=dft_control,&
            pao_env=pao_env)

       scf_env%ref_count=1
       scf_env%print_count=0
       last_scf_env_id=last_scf_env_id+1
       scf_env%id_nr=last_scf_env_id
       scf_env%print_count=0
       scf_env%iter_count=0
       scf_env%nelectron=0

! sets the method
       scf_env%method=diagonalisation_method_nr ! default with diagonalisation
       IF (dft_control%qs_control%pao.AND.ASSOCIATED(pao_env)) THEN
          scf_env%method=pao_mehod_nr
       END IF
       IF (scf_control%use_ot) scf_env%method=ot_method_nr

       SELECT CASE(scf_env%method)
       CASE(ot_method_nr)
          scf_env%mixing_method=no_mixing_nr
       CASE default
          scf_env%mixing_method=ao_mixing_nr
       END SELECT

       scf_env%iter_param=0.0_wp
       scf_env%iter_delta=0.0_wp
       scf_env%p_mix=scf_control%p_mix
       scf_env%iter_method=""
       scf_env%print_iter_line=.TRUE.

       NULLIFY(scf_env%scf_work1, scf_env%scf_work2, scf_env%ortho,&
            scf_env%p_new, scf_env%ot_preconditioner, scf_env%qs_ot_env,&
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
  SUBROUTINE scf(qs_env,globenv,error)
    TYPE(qs_environment_type), INTENT(INOUT), TARGET :: qs_env
    TYPE(global_environment_type), INTENT(IN)  :: globenv
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='scf',&
         routineP=moduleN//':'//routineN
    TYPE(qs_scf_env_type), POINTER :: scf_env
    TYPE(qs_environment_type), POINTER :: my_qs_env


    failure=.FALSE.
    NULLIFY(scf_env,my_qs_env)

!CPPrecondition(associated(qs_env),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL get_qs_env(qs_env,scf_env=scf_env,error=error)
       IF (.not.ASSOCIATED(scf_env)) THEN
          CALL scf_env_create(scf_env, qs_env=qs_env, error=error)
          CALL set_qs_env(qs_env,scf_env=scf_env,error=error)
       END IF
       my_qs_env => qs_env
       CALL scf_env_do_scf(scf_env,my_qs_env,globenv,error=error)
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
!!
!!*** **********************************************************************
  SUBROUTINE scf_env_do_scf(scf_env,qs_env,globenv,error)

    TYPE(qs_scf_env_type), POINTER :: scf_env
    TYPE(global_environment_type), INTENT(IN)  :: globenv
    TYPE(qs_environment_type), POINTER   :: qs_env
    TYPE(cp_error_type),INTENT(inout),OPTIONAL :: error

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routineN = "scf_env_do_scf"
    CHARACTER(LEN=*), PARAMETER :: routineP =&
         moduleN//':'//routineN

!   *** Local variables ***

    TYPE(dft_control_type), POINTER   :: dft_control
    REAL(wp) :: diis_error,t1,t2
    INTEGER  :: handle,handle2,output_unit,ispin
    LOGICAL  :: diis_step,do_level_shift,ionode,&
         use_jacobi,use_cholesky, failure,energy_only, rebuild
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: h,s,&
         ks, p
    TYPE(qs_rho_type), POINTER :: rho
    TYPE(qs_energy_type), POINTER :: energy
    TYPE(scf_control_type), POINTER   :: scf_control
    TYPE(cp_logger_type), POINTER :: logger
    TYPE(pw_env_type), POINTER :: pw_env
    TYPE(qs_charges_type), POINTER :: qs_charges
    TYPE(qs_ks_env_type), POINTER :: ks_env
    TYPE(mo_set_p_type), DIMENSION(:), POINTER    :: c
    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), POINTER    :: particle_set
    TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER :: scf_work1
    TYPE(cp_full_matrix_type), POINTER  :: ortho,scf_work2

!   ---------------------------------------------------------------------------

    NULLIFY(dft_control,h,s,ks,p,rho,energy,scf_control,logger,&
         pw_env, qs_charges,ks_env,c,atomic_kind_set,particle_set,&
         scf_work1,ortho,scf_work2)
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

    CALL timeset(routineN,"I","",handle)

    ionode = logger%para_env%source==logger%para_env%mepos
    output_unit = cp_logger_get_default_unit_nr(logger)
    energy_only = .FALSE.

    IF (ionode.AND.logger%print_keys%scf) THEN
       WRITE (UNIT=output_unit,FMT="(/,/,T2,A)")&
            "SCF WAVEFUNCTION OPTIMIZATION"
    END IF


    CALL init_scf_run(scf_env=scf_env,qs_env=qs_env,&
         globenv=globenv,&
         error=error)

    CALL get_qs_env(qs_env=qs_env,&
         h=h,&
         s=s,energy=energy,&
         particle_set=particle_set,&
         pw_env=pw_env,&
         qs_charges=qs_charges,&
         ks_env=ks_env)

    scf_work1 => scf_env%scf_work1
    scf_work2 => scf_env%scf_work2
    ortho => scf_env%ortho

    CALL get_qs_env(qs_env,k=ks,rho=rho,c=c)
    p => rho%rho_ao

    scf_env%iter_count = 0
    diis_step = .FALSE.
    rebuild = .TRUE.
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

       CALL timeset("scf_iter","I","",handle2)

       IF (ionode) CALL m_flush(output_unit)

       scf_env%iter_count = scf_env%iter_count + 1

! ** here qs_env%rho%rho_r and qs_env%rho%rho_g should be up to date

       CALL qs_ks_update_qs_env(ks_env,qs_env=qs_env,&
            rebuild=rebuild,error=error,&
            calculate_forces=.FALSE.,just_energy=energy_only)

!FM      CALL pw_env_flush_cache(pw_env,error=error)

       DO ispin=1,dft_control%nspins

          IF (logger%print_keys%density_matrix.AND.&
               logger%print_keys%each_scf_step) THEN
             CALL write_sparse_matrix(p(ispin)%matrix,4,6,qs_env,globenv)
          END IF

          IF (logger%print_keys%kohn_sham_matrix.AND.logger%print_keys%each_scf_step) THEN
             CALL write_sparse_matrix(ks(ispin)%matrix,4,6,qs_env,globenv)
          END IF

       ENDDO

       rebuild = (MODULO(scf_env%iter_count,scf_control%nrebuild) == 0)

       ispin=1

       scf_env%p_mix = scf_control%p_mix
       scf_env%iter_param = 0.0_wp
       IF (scf_env%mixing_method==ao_mixing_nr) scf_env%iter_param = scf_env%p_mix

       SELECT CASE (scf_env%method)
       CASE (pao_mehod_nr) ! pao
!FM         call pao_scf_iterate(pao_scf=scf_env%pao_scf, scf_env=scf_env,&
!FM              qs_env=qs_env, globenv=globenv, error=error)
!FM
       CASE(diagonalisation_method_nr) ! diagonalisation (default)

          DO ispin=1,dft_control%nspins 
             CALL copy_real_matrix_to_full_matrix(ks(ispin)%matrix,scf_work1(ispin)%matrix)
          ENDDO

          IF (scf_env%iter_count > 1) THEN
             CALL qs_diis_b_step(scf_env%scf_diis_buffer,c,scf_work1,&
                  scf_work2,scf_env%iter_delta,diis_error,diis_step,&
                  qs_env,error=error)
          END IF

          scf_env%iter_delta=0.0_wp
          do_level_shift = ((scf_control%level_shift /= 0.0_wp).AND.&
               ((scf_control%density_guess == "CORE").OR.(scf_env%iter_count > 1)))

          IF (diis_step) THEN
             scf_env%p_mix = 1.0_wp
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

          DO ispin=1,dft_control%nspins

             CALL eigensolver(scf_work1(ispin)%matrix,c(ispin)%mo_set,ortho,scf_work2,&
                  do_level_shift,scf_control%level_shift,globenv,&
                  use_arpack=scf_control%use_arpack,&
                  use_cholesky=use_cholesky,&
                  work_syevx=scf_control%work_syevx,&
                  use_jacobi=use_jacobi,&
                  eps_jacobi=scf_control%eps_jacobi,&
                  jacobi_threshold=scf_control%jacobi_threshold,&
                  smear=scf_control%smear)

             IF (logger%print_keys%each_scf_step) THEN
                CALL write_mo_set(c(ispin)%mo_set,atomic_kind_set,particle_set,4,6,globenv)
             END IF

             CALL calculate_density_matrix(c(ispin)%mo_set,&
                  scf_env%p_new(ispin)%matrix,error=error)

          ENDDO

       CASE(ot_method_nr) ! orbital transforms

! in case of LSD the first spin qs_ot_env will drive the minimization

          CALL ot_scf_mini(c,ks,s(1)%matrix,energy%total, &
               energy_only,scf_env%iter_delta,scf_env%qs_ot_env,globenv)

          DO ispin=1,dft_control%nspins
             CALL calculate_density_matrix(c(ispin)%mo_set,&
                  scf_env%p_new(ispin)%matrix,error=error)
          ENDDO

          scf_env%iter_method=scf_env%qs_ot_env(1)%OT_METHOD_FULL
          scf_env%iter_param=scf_env%qs_ot_env(1)%ds_min

       END SELECT

! ** at this point scf_env%p_new should be valid **

       CALL scf_env_density_mixing(scf_env,qs_env=qs_env,error=error)

       t2 = m_walltime()

       IF (ionode.AND.logger%print_keys%scf.and.scf_env%print_iter_line) THEN
          WRITE (UNIT=output_unit,&
               FMT="(T2,I5,2X,A,T22,E10.2,T32,F8.2,T40,2F20.10)")&
               scf_env%iter_count,TRIM(scf_env%iter_method),&
               scf_env%iter_param,t2 - t1,scf_env%iter_delta,energy%total
       END IF

! ** convergece check
       IF (scf_env%iter_delta < scf_control%eps_scf) THEN
          IF (ionode.AND.logger%print_keys%scf) THEN
             WRITE(UNIT=output_unit,FMT="(/,T3,A,/)")&
                  "*** SCF run converged ***"
          END IF
          CALL timestop(0.0_wp,handle2)
          EXIT scf_loop
       ELSE IF (scf_env%iter_count == scf_control%max_scf) THEN
          IF (ionode.AND.logger%print_keys%scf) THEN
             WRITE(UNIT=output_unit,FMT="(/,T3,A,/)")&
                  "*** SCF run NOT converged ***"
          END IF
          CALL timestop(0.0_wp,handle2)
          EXIT scf_loop
       END IF

       IF (.NOT. logger%print_keys%timecumul) t1 = m_walltime()

! ** update qs_env%rho
       CALL qs_rho_update_rho(rho, rho_ao_new=scf_env%p_new, qs_env=qs_env,&
            rebuild=rebuild, error=error)
       CALL qs_ks_did_change(ks_env,qs_env,rho_changed=.TRUE.,error=error)      

       CALL timestop(0.0_wp,handle2)

    END DO scf_loop

!   *** Write restart file ***
 
     CALL write_mo_set(c,atomic_kind_set,particle_set,globenv,&
          id_nr=qs_env%id_nr)

    IF (ionode.AND.logger%print_keys%scf) THEN
       WRITE (UNIT=output_unit,FMT="(/,(T3,A,T40,2F20.10))")&
            "Total electronic density (r-space): ",&
            SUM(rho%tot_rho_r),&
            SUM(rho%tot_rho_r)+ REAL(scf_env%nelectron,wp),&
            "Total core charge density (r-space):",&
            qs_charges%total_rho_core_rspace,&
            qs_charges%total_rho_core_rspace - REAL(scf_env%nelectron+dft_control%charge,wp)
       WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
            "Total charge density (r-space):     ",&
            SUM(rho%tot_rho_r)+&
            qs_charges%total_rho_core_rspace,&
            "Total charge density (g-space):     ",qs_charges%total_rho_gspace
       WRITE (UNIT=output_unit,FMT="(/,(T3,A,T60,F20.10))")&
            "Overlap energy of the core charge distribution:",energy%core_overlap,&
            "Self energy of the core charge distribution:   ",energy%core_self,&
            "Core Hamiltonian energy:                       ",energy%core,&
            "Hartree energy:                                ",energy%hartree,&
            "Exchange-correlation energy:                   ",energy%exc,&
            "Total energy:                                  ",energy%total
    END IF

!   *** Undo mixing of the density matrix (restore original density) ***

    CALL scf_env_density_mixing(scf_env,qs_env=qs_env,invert=.TRUE.,&
         error=error)

! ** update qs_env%rho
    CALL qs_rho_update_rho(rho, rho_ao_new=scf_env%p_new, qs_env=qs_env,&
         rebuild=rebuild, error=error)
    CALL qs_ks_did_change(ks_env,qs_env,rho_changed=.TRUE.,error=error)

! *** print out info
    CALL scf_env_post_calc(scf_env=scf_env, qs_env=qs_env, &
         globenv=globenv, error=error)

! *** cleanup
    CALL scf_env_cleanup(scf_env,qs_env=qs_env,error=error)

    CALL timestop(0.0_wp,handle)

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
  SUBROUTINE scf_env_density_mixing(scf_env,qs_env,invert,error)
    TYPE(qs_scf_env_type), POINTER :: scf_env
    TYPE(qs_environment_type), POINTER :: qs_env
    LOGICAL, INTENT(in), OPTIONAL :: invert
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='scf_env_density_mixing',&
         routineP=moduleN//':'//routineN
    REAL(kind=wp) :: tmp
    LOGICAL :: my_invert
    INTEGER :: ispin

    failure=.FALSE.
    my_invert=.FALSE.
    IF (PRESENT(invert)) my_invert=invert

    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(qs_env%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       SELECT CASE (scf_env%mixing_method)
       CASE(no_mixing_nr)
          ;
       CASE(ao_mixing_nr)
!FM        scf_env%iter_delta=0.0_wp
          DO ispin=1,SIZE(scf_env%p_new)
             IF (my_invert) THEN
                CPPrecondition(scf_env%p_mix/=0.0_wp,cp_failure_level,routineP,error,failure)
                IF (scf_env%p_mix/=1.0) THEN
                   CALL cp_sm_scale_and_add(matrix_a=scf_env%p_new(ispin)%matrix,&
                        alpha=1.0_wp/scf_env%p_mix,&
                        matrix_b=qs_env%rho%rho_ao(ispin)%matrix,&
                        beta=(scf_env%p_mix-1.0_wp)/scf_env%p_mix,&
                        error=error)
                END IF
             ELSE
                CALL cp_sm_mix(m1=scf_env%p_new(ispin)%matrix,&
                     m2=qs_env%rho%rho_ao(ispin)%matrix,&
                     p_mix=scf_env%p_mix,&
                     delta=tmp,&
                     para_env=qs_env%para_env,&
                     error=error)
                scf_env%iter_delta=MAX(scf_env%iter_delta,tmp)
             END IF
          END DO
       CASE default
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP,"unknown density mixing method:"//&
               cp_to_string(scf_env%mixing_method),error,failure)
       END SELECT
    END IF
  END SUBROUTINE scf_env_density_mixing
!***************************************************************************

! *****************************************************************************

  SUBROUTINE init_scf_run(scf_env,qs_env,globenv,error)

!   Purpose: Initialise a SCF run.

!   History: - Creation (30.04.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(qs_scf_env_type), POINTER :: scf_env
    TYPE(global_environment_type), INTENT(in) :: globenv
    TYPE(qs_environment_type), POINTER        :: qs_env
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routineN = "init_scf_run"
    CHARACTER(LEN=*), PARAMETER :: routineP =&
         moduleN//':'//routineN

!   *** Local variables ***

    INTEGER :: handle,homo,nao,ndep,ispin,nelectron_spin,&
         output_unit, stat
    LOGICAL :: ionode, gth_potential_present, failure
    TYPE(scf_control_type), POINTER :: scf_control
    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: h,s,k
    TYPE(dft_control_type), POINTER   :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), POINTER    :: c
    TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER :: scf_work1
    TYPE(cp_full_matrix_type), POINTER  :: ortho,scf_work2
    TYPE(cp_logger_type), POINTER :: logger
    TYPE(cp_full_matrix_type), POINTER :: mo_eigenvectors

!   ---------------------------------------------------------------------------

    CALL timeset("init_scf_run","I","",handle)

    NULLIFY(scf_control,atomic_kind_set,h,s,k,dft_control,c,scf_work1,&
         ortho,scf_work2, mo_eigenvectors)
    failure=.FALSE.
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
         c=c,k=k,&
         h=h,s=s, error=error)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         gth_potential_present=gth_potential_present)

    ionode = logger%para_env%source == logger%para_env%mepos
    output_unit = cp_logger_get_default_unit_nr(logger,local=.FALSE.)

    scf_work1 => scf_env%scf_work1
    scf_work2 => scf_env%scf_work2
    ortho => scf_env%ortho

!   *** Redistribute the core Hamiltonian matrix in ***
!   *** order to eliminate redundant atomic blocks  ***

    CALL copy_real_matrix_to_full_matrix(h(1)%matrix,scf_work1(1)%matrix)
    IF ((qs_env%para_env%num_pe > 1).AND.gth_potential_present) THEN
       CALL sm_pool_give_back_matrix(qs_env%H_redundant_sm_pools(1)%pool,&
            qs_env%h(1)%matrix,error=error)
       CALL sm_pool_create_matrix(qs_env%S_sm_pools(1)%pool,h(1)%matrix,&
            name="qs_env"//TRIM(ADJUSTL(cp_to_string(qs_env%id_nr)))//"%h-1",&
            error=error)
       CALL copy_full_matrix_to_real_matrix(scf_work1(1)%matrix,&
            qs_env%h(1)%matrix)
    END IF

    IF (scf_env%method==ot_method_nr) THEN 
       CALL cp_fm_release(scf_work1(1)%matrix)
       DEALLOCATE(scf_work1)
       NULLIFY(scf_work1, scf_env%scf_work1)
    ENDIF

!   *** done redistribute ***

! update ecore
    CALL calculate_ecore_self(qs_env)
    CALL calculate_ecore_overlap(qs_env,globenv,.FALSE.)
    IF (ionode.AND.logger%print_keys%scf_energies) THEN
       WRITE (UNIT=output_unit,FMT="(/,(T3,A,T60,F20.10))")&
            "Overlap energy of the core charge distribution:",&
            qs_env%energy%core_overlap,&
            "Self energy of the core charge distribution:   ",&
            qs_env%energy%core_self
    END IF

! updates the total number of electrons
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         nelectron=scf_env%nelectron)
    scf_env%nelectron = scf_env%nelectron - dft_control%charge
! print occupation numbers
    IF (ionode.AND.logger%print_keys%scf) THEN
       DO ispin=1,dft_control%nspins
          CALL get_mo_set(mo_set=c(ispin)%mo_set,&
               homo=homo, &
               nelectron=nelectron_spin, &
               nao=nao)
          WRITE (UNIT=output_unit,fmt="(/,A,I2)") " For spin ",ispin
          WRITE (UNIT=output_unit,FMT="(/,(T3,A,I10))")&
               "Number of electrons:                    ",nelectron_spin,&
               "Number of occupied orbitals:            ",homo
       ENDDO
       WRITE(UNIT=output_unit,FMT="(/,(T3,A,I10))")&
            "Number of orbital functions:            ",nao
    END IF

! calc ortho matrix
    ndep = 0
    IF (scf_env%method /= ot_method_nr) THEN
       IF (scf_control%use_cholesky .OR. scf_control%use_arpack) THEN
          CALL copy_real_matrix_to_full_matrix(s(1)%matrix,ortho)
          CALL cp_fm_cholesky_decompose(ortho)
       ELSE
          CALL calculate_ortho_matrix(ortho,s(1)%matrix,scf_work2,ndep,&
               scf_control%eps_eigval,scf_control%work_syevx)
          IF (ionode.and.logger%print_keys%scf) THEN
             WRITE (UNIT=output_unit,FMT="(T3,A,I10)")&
                  "Number of independent orbital functions:",nao - ndep
          END IF
       ENDIF
    ENDIF
    IF (logger%print_keys%ortho_matrix) THEN
       CALL write_blacs_matrix(ortho,4,6,qs_env,globenv)
    END IF

! initializes rho and the mos
    CALL scf_env_initial_rho_setup(scf_env,qs_env=qs_env,&
         globenv=globenv,error=error)


! *** method dependent initializations ***
    SELECT CASE (scf_env%method)
    CASE (pao_mehod_nr)
!FM       if (.not.associated(scf_env%pao_scf)) then
!FM          call pao_scf_env_create(scf_env%pao_scf,scf_env=scf_env,qs_env=qs_env,&
!FM               globenv=globenv, error=error)
!FM       end if
    CASE (diagonalisation_method_nr)
       IF (.NOT.ASSOCIATED(scf_env%scf_diis_buffer)) THEN       
          CALL qs_diis_b_create(scf_env%scf_diis_buffer,&
               nbuffer=scf_control%max_diis,&
               error=error)
       END IF
       CALL qs_diis_b_clear(scf_env%scf_diis_buffer,error=error)
    CASE (ot_method_nr)
       IF (.NOT. ASSOCIATED(scf_env%ot_preconditioner)) THEN
          ALLOCATE(scf_env%ot_preconditioner, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          CALL init_preconditioner(scf_env%ot_preconditioner,&
               para_env=qs_env%para_env,&
               blacs_env=qs_env%blacs_env)
       END IF
       IF (.NOT.ASSOCIATED(scf_env%qs_ot_env)) THEN
          ALLOCATE(scf_env%qs_ot_env(dft_control%nspins),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          CALL ot_scf_read_input(scf_env%qs_ot_env,globenv)

          ! might need the KS matrix to init properly
          CALL qs_ks_update_qs_env(qs_env%ks_env,qs_env=qs_env,&
               rebuild=.TRUE.,error=error,&
               calculate_forces=.FALSE.,just_energy=.FALSE.)
       
          IF (scf_env%qs_ot_env(1)%settings%preconditioner_type .NE. "NONE") THEN
             ! constuct the preconditioner based on the first spin set. this might or might not be reasonable. Presumably it is.
             CALL get_mo_set(mo_set=c(1)%mo_set,&
                  eigenvectors=mo_eigenvectors)
             CALL make_preconditioner(scf_env%ot_preconditioner, &
                  scf_env%qs_ot_env(1)%settings%preconditioner_type, &
                  k(1)%matrix,s(1)%matrix,&
                  qs_env%kinetic(1)%matrix, &
                  mo_eigenvectors,0.2_wp) 
          ENDIF
          
          CALL ot_scf_init(c,s(1)%matrix,scf_env%qs_ot_env,globenv)
          
          IF (scf_env%qs_ot_env(1)%settings%preconditioner_type .NE. "NONE") THEN
             DO ispin=1,dft_control%nspins
                CALL qs_ot_new_preconditioner(scf_env%qs_ot_env(ispin),&
                     scf_env%ot_preconditioner) 
             ENDDO
          ENDIF
       ENDIF
    END SELECT

    CALL timestop(0.0_wp,handle)

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
    TYPE(qs_scf_env_type), POINTER :: scf_env
    TYPE(qs_environment_type), POINTER        :: qs_env
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "scf_env_check_i_alloc",&
         routineP = moduleN//":"//routineN
    INTEGER :: stat,nao,ispin,handle,nrow_block
    LOGICAL :: failure
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: k
    TYPE(cp_fm_pool_p_type),DIMENSION(:), POINTER :: ao_mo_fm_pools
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: s
    TYPE(cp_fm_struct_type), POINTER :: ao_mo_fmstruct, ao_ao_fmstruct
    TYPE(cp_full_matrix_type), POINTER :: eigenvectors
    TYPE(dft_control_type), POINTER   :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), POINTER    :: c
    TYPE(cp_sm_pool_p_type), DIMENSION(:), POINTER :: S_sm_pools
    TYPE(qs_ks_env_type), POINTER :: ks_env

!   ---------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

    NULLIFY(k,ao_mo_fm_pools, s, ao_mo_fmstruct, ao_ao_fmstruct,&
         eigenvectors, dft_control, c, S_sm_pools, ks_env)
    failure=.FALSE.

    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)

    CALL get_qs_env(qs_env=qs_env,&
         dft_control=dft_control,&
         c=c,&
         k=k,&
         ks_env=ks_env,&
         ao_mo_fm_pools=ao_mo_fm_pools,&
         S_sm_pools=S_sm_pools,&
         s=s)

!   *** finish initialization of the MOs ***
    CPPrecondition(ASSOCIATED(c),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       DO ispin=1,SIZE(c)
          CALL get_mo_set(c(ispin)%mo_set,eigenvectors=eigenvectors)
          IF (.NOT.ASSOCIATED(eigenvectors)) THEN
             CALL init_mo_set(c(ispin)%mo_set,&
                  ao_mo_fm_pools(ispin)%pool,&
                  name="qs_env"//TRIM(ADJUSTL(cp_to_string(qs_env%id_nr)))//&
                  "%mo"//TRIM(ADJUSTL(cp_to_string(ispin))),&
                  error=error)
          END IF
       END DO
    END IF

!   *** Allocate the distributed SCF matrices ***

    IF (.not.ASSOCIATED(scf_env%scf_work1).OR. &
         .not.ASSOCIATED(scf_env%scf_work2).OR.&
         .not.ASSOCIATED(scf_env%ortho)) THEN

       ao_mo_fmstruct => fm_pool_get_mstruct(ao_mo_fm_pools(1)%pool,&
            error=error)
       CALL cp_fm_struct_get(ao_mo_fmstruct, nrow_block=nrow_block,&
            error=error)
       CALL get_mo_set(c(1)%mo_set,nao=nao)
       CALL cp_fm_struct_create(fmstruct=ao_ao_fmstruct,&
            nrow_block=nrow_block,ncol_block=nrow_block,&
            nrow_global=nao, ncol_global=nao,&
            template_fmstruct=ao_mo_fmstruct, error=error)


       IF (.not.ASSOCIATED(scf_env%scf_work1)) THEN
          IF (scf_env%method/=ot_method_nr) THEN
             ALLOCATE(scf_env%scf_work1(dft_control%nspins), stat=stat)
          ELSE
             ALLOCATE(scf_env%scf_work1(1), stat=stat)
          END IF
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO ispin=1,SIZE(scf_env%scf_work1)
             NULLIFY(scf_env%scf_work1(ispin)%matrix)
             CALL cp_fm_create(scf_env%scf_work1(ispin)%matrix,&
                  matrix_struct=ao_ao_fmstruct,&
                  name="SCF"//TRIM(ADJUSTL(cp_to_string(scf_env%id_nr)))//&
                  "WORK_MATRIX-1-"//TRIM(ADJUSTL(cp_to_string(ispin))),&
                  error=error)
          ENDDO
       END IF
       IF (.not.ASSOCIATED(scf_env%ortho).and.scf_env%method/=ot_method_nr) THEN
          CALL cp_fm_create(scf_env%ortho,&
               matrix_struct=ao_ao_fmstruct,&
               name="SCF"//TRIM(ADJUSTL(cp_to_string(scf_env%id_nr)))//&
               "ORTHO_MATRIX",&
               error=error)
       END IF
       IF (.not.ASSOCIATED(scf_env%scf_work2).AND.&
            scf_env%method/=ot_method_nr) THEN
          CALL cp_fm_create(scf_env%scf_work2,&
               matrix_struct=ao_ao_fmstruct,&
               name="SCF"//TRIM(ADJUSTL(cp_to_string(scf_env%id_nr)))//&
               "WORK_MATRIX-2",&
               error=error)
       END IF
       CALL cp_fm_struct_release(ao_ao_fmstruct,error=error)

    END IF

!   *** Allocate k and put it in the QS environment ***

    IF (.not.ASSOCIATED(k)) THEN
       CALL sm_pools_create_matrix_vect(S_sm_pools,k,&
            name="SCF"//TRIM(ADJUSTL(cp_to_string(scf_env%id_nr)))//&
            "KOHN-SHAM_MATRIX-"//TRIM(ADJUSTL(cp_to_string(ispin))),&
            error=error)
       CALL set_qs_env(qs_env=qs_env,&
            k=k)
    END IF

!   *** allocate p_new ***
    IF (.not.ASSOCIATED(scf_env%p_new)) THEN
       CALL sm_pools_create_matrix_vect(S_sm_pools,scf_env%p_new,&
            name="SCF"//TRIM(ADJUSTL(cp_to_string(scf_env%id_nr)))//&
            "DENSITY",&
            error=error)
    END IF

!   *** allocate the ks env **
    IF (.not.ASSOCIATED(ks_env)) THEN
       CALL qs_ks_create(ks_env,qs_env=qs_env)
       CALL set_qs_env(qs_env, ks_env=ks_env)
    END IF

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE scf_env_check_i_alloc
!***************************************************************************

  SUBROUTINE calculate_ortho_matrix(ortho,s,work,ndep,eps_eigval,&
       work_syevx)

!   Purpose: Calculate the orthogonalization matrix (S**(-1/2))

!   History: - Creation (01.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(cp_full_matrix_type), POINTER        :: ortho,work
    TYPE(real_matrix_type), POINTER           :: s
    INTEGER, INTENT(OUT)                      :: ndep
    REAL(kind=wp), INTENT(in)                 :: eps_eigval,work_syevx

!   *** Local variables ***

    INTEGER :: handle,nao

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_ortho_matrix","I","",handle)

    CALL copy_real_matrix_to_full_matrix(s,ortho)
    CALL power_blacs_matrix(ortho,work,-0.5_wp,eps_eigval,ndep,work_syevx)
    CALL symmetrise_blacs_matrix(ortho,work)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE calculate_ortho_matrix

! *****************************************************************************

  SUBROUTINE eigensolver(ks,mo_set,ortho,work,do_level_shift,level_shift,&
       globenv,use_cholesky, use_arpack,work_syevx,&
       use_jacobi,eps_jacobi,jacobi_threshold,smear)

!   Purpose: Diagonalise the Kohn-Sham matrix to get a new set of MO eigen-
!            vectors and MO eigenvalues.
!   ks will be modified

!   History: - Creation (01.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(cp_full_matrix_type), POINTER        :: ks,ortho,work
    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(mo_set_type), POINTER                :: mo_set
    LOGICAL, INTENT(INOUT)                    :: use_cholesky
    LOGICAL, INTENT(IN)                       :: do_level_shift,use_arpack,&
         use_jacobi
    REAL(wp), INTENT(IN)                      :: level_shift,work_syevx,&
         eps_jacobi,jacobi_threshold,&
         smear

!   *** Local variables ***
    INTEGER :: handle,imo,nao,nmo,homo,ncv
    TYPE(cp_full_matrix_type), POINTER :: mo_eigenvectors
    REAL(wp), DIMENSION(:), POINTER    :: mo_eigenvalues
    TYPE(op_symv_type)                 :: arpack_op

!   ---------------------------------------------------------------------------

    CALL timeset("eigensolver","I","",handle)

!   *** Diagonalise the Kohn-Sham matrix ***

    CALL get_mo_set(mo_set=mo_set,&
         nao=nao,&
         nmo=nmo,&
         homo=homo, &
         eigenvalues=mo_eigenvalues,&
         eigenvectors=mo_eigenvectors)

    IF (do_level_shift) THEN
       IF (use_cholesky .OR. use_arpack) THEN
          CALL stop_program("eigensolver","level shift not implemented")
       ENDIF
    ENDIF

    IF (use_cholesky) THEN

       CALL cp_fm_cholesky_reduce(ks,ortho)

       IF (use_jacobi) THEN
          CALL blacs_syevd(ks,work,mo_eigenvalues)
          use_cholesky = .FALSE.
          CALL cp_fm_cholesky_restore(work,nmo,ortho,mo_eigenvectors,"SOLVE")
          CALL copy_blacs_to_blacs_matrix(mo_eigenvectors,ortho)
       ELSE IF (smear /= 0.0_wp) THEN
          CALL blacs_syevd(ks,work,mo_eigenvalues)
          CALL cp_fm_cholesky_restore(work,nmo,ortho,mo_eigenvectors,"SOLVE")
       ELSE
          CALL blacs_syevx(ks,work,mo_eigenvalues,homo,work_syevx)
          CALL cp_fm_cholesky_restore(work,homo,ortho,mo_eigenvectors,"SOLVE")
       END IF

    ELSE IF (use_arpack) THEN ! not really recommended, since rather slow, except for huge basis sets

!TC       ncv = nmo + 64
       ncv = nmo*2
       CALL cp_fm_cholesky_reduce(ks,ortho)

       CALL setup_op_symv(arpack_op, ks, work, .FALSE.)
       CALL arpack_diag_sy(work,mo_eigenvalues,nmo,ncv,.FALSE.,globenv, &
            op_symv=arpack_op)

       CALL cp_fm_cholesky_restore(work,nmo,ortho,mo_eigenvectors,"SOLVE")

    ELSE IF (use_jacobi) THEN

       CALL cp_fm_symm("L","U",nao,homo,1.0_wp,ks,mo_eigenvectors,0.0_wp,work)
       CALL cp_fm_gemm("T","N",homo,nao-homo,nao,1.0_wp,work,mo_eigenvectors,&
            0.0_wp,ks,b_first_col=homo+1)

!     epsi =1e-11
!     Anz_max=4
!     N=ks%matrix_struct%nrow_global

!     *** Scale the elements with the inverse ***
!     *** difference of the old eigenvalues   ***

!     CALL blacs_block_jacobi_scaled(ks,work,mo_eigenvalues,&
!                                    jacobi_threshold,homo+1)

!     *** Klassisches Jacobi-Block-Verfrahren ***
!     *** Rotiere Elemente nur, wenn sie groesser als thresh sind ***

       CALL blacs_block_jacobi_classic(ks,mo_eigenvectors,mo_eigenvalues,&
            jacobi_threshold,homo+1)

    ELSE

       CALL cp_fm_symm("L","U",nao,nao,1.0_wp,ks,ortho,0.0_wp,work)
       CALL cp_fm_gemm("T","N",nao,nao,nao,1.0_wp,ortho,work,0.0_wp,ks)

       IF (do_level_shift) THEN
          DO imo=homo+1,nmo
             CALL blacs_add_to_element(ks,imo,imo,level_shift)
          END DO
       END IF

       CALL blacs_syevd(ks,work,mo_eigenvalues)

       CALL cp_fm_gemm("N","N",nao,nmo,nao,1.0_wp,ortho,work,0.0_wp,&
            mo_eigenvectors)

       IF (do_level_shift) THEN

          CALL correct_mo_eigenvalues(mo_set,level_shift)

!     *** Use last MO set as orthogonalization matrix ***

          CALL copy_blacs_to_blacs_matrix(mo_eigenvectors,ortho)

       END IF

    END IF

    IF (smear /= 0.0_wp) CALL set_mo_occupation(mo_set,smear)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE eigensolver

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

    TYPE(real_matrix_type), POINTER :: m1,m2
    REAL(wp), INTENT(IN) :: p_mix
    REAL(wp), INTENT(OUT) :: delta
    TYPE(cp_para_env_type), POINTER :: para_env
    TYPE(cp_error_type), INTENT(inout), OPTIONAL:: error

!   *** Local variables ***

    TYPE(real_block_node_type), POINTER :: block_node

    CHARACTER(len=*), PARAMETER :: routineN="cp_sm_mix",&
         routineP=moduleN//':'//routineN
    REAL(wp) :: r_mix
    INTEGER  :: handle,i,iblock_col,iblock_row,j,nblock_row

    REAL(wp), DIMENSION(:,:), POINTER :: p_new_block,p_old_block

!   ---------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

    delta = 0.0_wp

    r_mix = 1.0_wp - p_mix

    CALL get_matrix_info(matrix=m2,nblock_row=nblock_row)

    DO iblock_row=1,nblock_row

       block_node => first_block_node(matrix=m1,&
            block_row=iblock_row)

       DO WHILE (ASSOCIATED(block_node))

          CALL get_block_node(block_node=block_node,&
               block_col=iblock_col,&
               BLOCK=p_new_block)

          CALL get_block_node(matrix=m2,&
               block_row=iblock_row,&
               block_col=iblock_col,&
               BLOCK=p_old_block)

          DO j=1,SIZE(p_new_block,2)
             DO i=1,SIZE(p_new_block,1)
                p_new_block(i,j) = p_new_block(i,j) - p_old_block(i,j)
                delta = MAX(delta,ABS(p_new_block(i,j)))
                p_new_block(i,j) = p_old_block(i,j) + p_mix*p_new_block(i,j)
             END DO
          END DO

          block_node => next_block_node(block_node)

       END DO

    END DO

    CALL mp_max(delta,para_env%group)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE cp_sm_mix

! *****************************************************************************
  SUBROUTINE calculate_first_density_matrix(scf_env,qs_env,globenv,error)

    TYPE(qs_scf_env_type), POINTER :: scf_env
    TYPE(qs_environment_type), POINTER :: qs_env
    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routineN =&
         "calculate_first_density_matrix"
    CHARACTER(LEN=*), PARAMETER :: routineP =&
         moduleN//":"//routineN

!   *** Local variables ***

    TYPE(all_potential_type), POINTER :: all_potential
    TYPE(atomic_kind_type), POINTER   :: atomic_kind
    TYPE(cp_full_matrix_type), POINTER  :: mo_eigenvectors,sv
    TYPE(gth_potential_type), POINTER :: gth_potential
    TYPE(gto_basis_set_type), POINTER :: orb_basis_set
    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), POINTER    :: particle_set

    TYPE(cp_full_matrix_type), POINTER        :: ortho,work2
    TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER        :: work1
    TYPE(mo_set_p_type),  DIMENSION(:), POINTER               :: mo_array
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER   :: s_sparse
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER   :: p_rmpv
    TYPE(scf_control_type), POINTER           :: scf_control
    TYPE(dft_control_type), POINTER :: dft_control

    REAL(wp) :: paa,saa,sum
    INTEGER  :: atom_a,first_sgfa,handle,homo,iatom,ikind,iset,isgf,isgfa,&
         ishell,istat,la,last_sgfa,maxl,natom,nelec,nmo,nset,&
         qs_env_id,restart_unit

    INTEGER, DIMENSION(:), ALLOCATABLE :: econf

    INTEGER, DIMENSION(:), POINTER   :: atom_list,elec_conf,nshell
    INTEGER, DIMENSION(:,:), POINTER :: first_sgf,l,last_sgf
    INTEGER :: nspin,ispin,nelectron

!   ---------------------------------------------------------------------------

    NULLIFY(all_potential,atomic_kind, mo_eigenvectors, sv,&
         gth_potential, orb_basis_set, atomic_kind_set, particle_set,&
         ortho,work2,work1,mo_array,s_sparse,p_rmpv,scf_control, &
         dft_control)

    CALL timeset("calculate_first_density_matrix","I","",handle)

    CALL get_qs_env(qs_env,atomic_kind_set=atomic_kind_set,&
         particle_set=particle_set,c=mo_array, s=s_sparse,&
         scf_control=scf_control, id_nr=qs_env_id, dft_control=dft_control,&
         error=error)
    nspin=dft_control%nspins

    p_rmpv => scf_env%p_new
    work1 => scf_env%scf_work1
    work2 => scf_env%scf_work2
    ortho => scf_env%ortho

    IF (scf_control%use_ot .AND. .NOT.  &
         (scf_control%density_guess == "RANDOM" .OR.  &
         (scf_control%density_guess == "RESTART" .AND. &
         scf_control%level_shift == 0.0_wp))) THEN
       CALL stop_program("calculate_first_density_matrix","OT needs GUESS RESTART : other options NYI")
    ENDIF

    IF (scf_control%density_guess == "RESTART") THEN

       CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set)
       CALL read_mo_set(mo_array,atomic_kind_set,particle_set,globenv,&
            id_nr=qs_env_id)

       DO ispin=1,nspin
          IF (scf_control%level_shift /= 0.0_wp) THEN
             CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
                  eigenvectors=mo_eigenvectors)
             CALL copy_blacs_to_blacs_matrix(mo_eigenvectors,ortho)
          END IF

          CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
               eigenvectors=mo_eigenvectors, homo=homo)
! ortho so that one can restart for different positions (basis sets?)
          CALL replicate_blacs_matrix(mo_eigenvectors,sv,"SV")
          CALL cp_sm_fm_multiply(s_sparse(1)%matrix,mo_eigenvectors,sv,homo,&
               para_env=qs_env%para_env)
          CALL blacs_make_basis(mo_eigenvectors,homo, &
               ortho=sv,otype="SV")
          CALL cp_fm_release(sv)
          CALL calculate_density_matrix(mo_array(ispin)%mo_set,&
               p_rmpv(ispin)%matrix,error=error)
       ENDDO

    ELSE IF (scf_control%density_guess == "RANDOM") THEN

       DO ispin=1,nspin
          CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
               eigenvectors=mo_eigenvectors, nmo=nmo)
          CALL cp_fm_init_random(mo_eigenvectors,nmo)
! ortho so that one can restart for different positions (basis sets?)
          CALL replicate_blacs_matrix(mo_eigenvectors,sv,"SV")
          CALL cp_sm_fm_multiply(s_sparse(1)%matrix,mo_eigenvectors,sv,nmo,&
               para_env=qs_env%para_env)
          CALL blacs_make_basis(mo_eigenvectors,nmo, &
               ortho=sv,otype="SV")
          CALL cp_fm_release(sv)
          CALL calculate_density_matrix(mo_array(ispin)%mo_set,&
               p_rmpv(ispin)%matrix,error=error)
       ENDDO

    ELSE IF (scf_control%density_guess == "CORE") THEN

! work1(1) should contain the core hamiltonian
! we need a copy the core matrix for every spin (the dumbest solution,
! we should copy the mos)..
       IF (dft_control%nspins.eq.2) THEN
          IF (.NOT. scf_control%use_ot) THEN
             CALL copy_blacs_to_blacs_matrix(work1(1)%matrix,&
                  work1(2)%matrix)
          ENDIF
       ENDIF
       DO ispin=1,nspin
          CALL eigensolver(work1(ispin)%matrix,mo_array(ispin)%mo_set,ortho,work2, &
               .FALSE.,0.0_wp,globenv,&
               use_arpack=scf_control%use_arpack, &
               use_cholesky=scf_control%use_cholesky, &
               work_syevx=scf_control%work_syevx,&
               use_jacobi=.FALSE.,&
               eps_jacobi=scf_control%eps_jacobi,&
               jacobi_threshold=scf_control%jacobi_threshold,&
               smear=0.0_wp)
          IF (scf_control%level_shift /= 0.0_wp) THEN
             CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
                  eigenvectors=mo_eigenvectors)
             CALL copy_blacs_to_blacs_matrix(mo_eigenvectors,ortho)
          END IF
          CALL calculate_density_matrix(mo_array(ispin)%mo_set,&
               p_rmpv(ispin)%matrix,error=error)
       ENDDO

    ELSE IF (scf_control%density_guess == "ATOMIC") THEN

       DO ispin=1,nspin
          CALL get_mo_set(mo_array(ispin)%mo_set,nelectron=nelectron)
          CALL cp_fm_set_all(work1(ispin)%matrix,0.0_wp)

          sum = 0.0_wp

          CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,maxl=maxl)

          ALLOCATE (econf(0:maxl),STAT=istat)
          IF (istat /= 0) CALL stop_memory(routineP,"econf",(maxl + 1)*int_size)

          CALL copy_real_matrix_to_full_matrix(s_sparse(1)%matrix,work2)

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

             maxl = SIZE(elec_conf) - 1
             econf(:) = 0
             econf(0:maxl) = elec_conf(0:maxl)

             CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                  nset=nset,&
                  nshell=nshell,&
                  l=l,&
                  first_sgf=first_sgf,&
                  last_sgf=last_sgf)

             DO iset=1,nset
                DO ishell=1,nshell(iset)
                   la = l(ishell,iset)
                   nelec = 4*la + 2
                   IF (econf(la) > 0) THEN
                      IF (econf(la) >= nelec) THEN
                         paa = 2.0_wp
                         econf(la) = econf(la) - nelec
                      ELSE
                         paa = REAL(2*econf(la),wp)/REAL(nelec,wp)
                         econf(la) = 0
                      END IF
                      first_sgfa = first_sgf(ishell,iset)
                      last_sgfa = last_sgf(ishell,iset)
                      DO isgfa=first_sgfa,last_sgfa
                         DO iatom=1,natom
                            atom_a = atom_list(iatom)
                            isgf = particle_set(atom_a)%first_sgf + isgfa - 1
                            CALL cp_fm_set_element(work1(ispin)%matrix,isgf,isgf,paa)
                            CALL blacs_get_element(work2,isgf,isgf,saa)
                            sum = sum + paa*saa
                         END DO
                      END DO
                   END IF
                END DO
             END DO

          END DO

          DEALLOCATE (econf,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routineP,"econf")

          CALL mp_sum(sum,qs_env%para_env%group)

          saa = REAL(nelectron,wp)/sum

          DO ikind=1,SIZE(atomic_kind_set)

             atomic_kind => atomic_kind_set(ikind)

             CALL get_atomic_kind(atomic_kind=atomic_kind,&
                  natom=natom,&
                  atom_list=atom_list,&
                  orb_basis_set=orb_basis_set)

             CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                  nset=nset,&
                  nshell=nshell,&
                  l=l,&
                  first_sgf=first_sgf,&
                  last_sgf=last_sgf)

             DO iset=1,nset
                DO ishell=1,nshell(iset)
                   first_sgfa = first_sgf(ishell,iset)
                   last_sgfa = last_sgf(ishell,iset)
                   DO isgfa=first_sgfa,last_sgfa
                      DO iatom=1,natom
                         atom_a = atom_list(iatom)
                         isgf = particle_set(atom_a)%first_sgf + isgfa - 1
                         CALL blacs_get_element(work1(ispin)%matrix,isgf,isgf,paa)
                         paa = paa*saa
                         CALL cp_fm_set_element(work1(ispin)%matrix,isgf,isgf,paa)
                      END DO
                   END DO
                END DO
             END DO
          END DO

          CALL copy_full_matrix_to_real_matrix(work1(ispin)%matrix,p_rmpv(ispin)%matrix)

       ENDDO

    ELSE

       CALL stop_program(routineP,&
            "An invalid keyword for the initial density "//&
            "guess was specified")

    END IF

    CALL timestop(0.0_wp,handle)

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
    TYPE(qs_scf_env_type), POINTER :: scf_env
    TYPE(qs_environment_type), POINTER :: qs_env
    TYPE(global_environment_type), INTENT(in), TARGET :: globenv
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure, extrapolated
    CHARACTER(len=*), PARAMETER :: routineN='scf_env_initial_rho_setup',&
         routineP=moduleN//':'//routineN
    INTEGER :: handle

    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.
    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       extrapolated=.FALSE.
       IF (ASSOCIATED(qs_env%wf_history)) THEN
          CALL wfi_extrapolate(qs_env%wf_history, &
               qs_env=qs_env, dt=1.0_wp, extrapolated=extrapolated,&
               error=error)
       END IF
       IF (.NOT.extrapolated) THEN
          CALL calculate_first_density_matrix(scf_env=scf_env,qs_env=qs_env,&
               globenv=globenv, error=error)
! update qs_env%rho
          CALL qs_rho_update_rho(qs_env%rho, rho_ao_new=scf_env%p_new,&
               qs_env=qs_env,&
               rebuild=.TRUE., error=error)
          CALL qs_ks_did_change(qs_env%ks_env,qs_env,rho_changed=.TRUE.,&
               error=error)

       END IF
    END IF
    CALL timestop(0.0_wp,handle)
  END SUBROUTINE scf_env_initial_rho_setup
!***************************************************************************

!!****f* qs_scf/scf_env_post_calc [1.0] *
!!
!!   NAME
!!     scf_env_post_calc
!!
!!   FUNCTION
!!     writes out the requested data and does various calculaion at the
!!     end of an scf run
!!
!!   NOTES
!!     globenv should be eliminated
!!
!!   ARGUMENTS
!!     - scf_env: the scf_env whose info should be written out
!!     - qs_env: the qs_env in which the qs_env lives
!!     - globenv: the global environement (should be removed)
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
  SUBROUTINE scf_env_post_calc(scf_env,qs_env,globenv,error)
    TYPE(qs_scf_env_type), POINTER :: scf_env
    TYPE(qs_environment_type), POINTER :: qs_env
    TYPE(global_environment_type), INTENT(in), TARGET :: globenv
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure, ionode
    CHARACTER(len=*), PARAMETER :: routineN='scf_env_post_calc',&
         routineP=moduleN//':'//routineN
    TYPE(dft_control_type), POINTER :: dft_control
    TYPE(pw_env_type), POINTER :: pw_env
    TYPE(pw_pool_type), POINTER :: auxbas_pw_pool
    TYPE(pw_pool_p_type), DIMENSION(:), POINTER :: pw_pools
    TYPE(coeff_type) :: wf_r, wf_g
    REAL(kind=wp) :: total_abs_spin_dens
    TYPE(coeff_type), DIMENSION(:), POINTER :: mgrid_gspace, mgrid_rspace
    TYPE(mo_set_p_type), DIMENSION(:), POINTER :: c
    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), POINTER :: particle_set
    TYPE(qs_rho_type), POINTER :: rho
    TYPE(cp_logger_type), POINTER :: logger
    TYPE(cp_full_matrix_type), POINTER :: mo_eigenvectors
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: ks_rmpv, s
    TYPE(qs_loc_env_type), POINTER :: qs_loc_env
    TYPE(scf_control_type), POINTER :: scf_control
    TYPE(cp_full_matrix_type), POINTER :: unoccupied_orbs
    TYPE(real_matrix_type), POINTER :: pfull
    INTEGER :: ispin,n, output_unit, stat
    CHARACTER(LEN=6)   :: extension
    CHARACTER(LEN=200) :: file_name
    INTEGER            :: file_unit,ivector,nao,nlumos, homo
    CHARACTER(LEN=80)  :: filename

    failure=.FALSE.
    NULLIFY(dft_control,pw_env,auxbas_pw_pool,pw_pools,mgrid_gspace,&
         mgrid_rspace,c,atomic_kind_set,particle_set,rho, mo_eigenvectors,&
         ks_rmpv,qs_loc_env, scf_control,unoccupied_orbs,pfull)
    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    output_unit= cp_logger_get_default_unit_nr(logger,local=.FALSE.)

    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL get_qs_env(qs_env,dft_control=dft_control,&
            c=c,atomic_kind_set=atomic_kind_set,particle_set=particle_set,&
            rho=rho,k=ks_rmpv,scf_control=scf_control,s=s,error=error)

       IF (ASSOCIATED(qs_env%wf_history)) THEN
          CALL wfi_update(qs_env%wf_history,qs_env=qs_env,dt=1.0_wp,&
               error=error)
       END IF

!   *** Write last wavefunction to screen ??????? ***

       DO ispin=1,dft_control%nspins
          CALL write_mo_set(c(ispin)%mo_set,atomic_kind_set,particle_set,4,6,globenv)
       ENDDO


!   *** Integrated absolute spin density ***
       IF (dft_control%nspins.eq.2) THEN
          CALL get_qs_env(qs_env=qs_env,pw_env=pw_env)
          CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,&
               pw_pools=pw_pools)
          CALL pw_pool_init_coeff(auxbas_pw_pool,wf_r,&
               use_data = REALDATA3D,&
               in_space = REALSPACE, error=error)
          CALL coeff_copy(rho%rho_r(1),wf_r)
          CALL coeff_reduce(rho%rho_r(2),wf_r)
          total_abs_spin_dens=calculate_total_abs_rho(wf_r)
          IF (globenv%ionode) WRITE(UNIT=output_unit,fmt='(/,(T3,A,T60,F20.10))') &
               "Integrated absolute spin density : ",total_abs_spin_dens
          CALL pw_pool_give_back_coeff(auxbas_pw_pool,wf_r, error=error)
       ENDIF

       IF (globenv%print%mulliken_analysis)  &
            CALL mulliken_charges(c,s(1)%matrix,globenv%scr,globenv%ionode)

!   *** write additional information about the converged system
       IF (globenv%print%e_density_cube) THEN ! write the electron density
          IF (dft_control%nspins.ne.1) THEN
             CALL get_qs_env(qs_env=qs_env,pw_env=pw_env)
             CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,&
                  pw_pools=pw_pools)
             CALL pw_pool_init_coeff(auxbas_pw_pool,wf_r,&
                  use_data = REALDATA3D,&
                  in_space = REALSPACE, error=error)
             CALL coeff_copy(rho%rho_r(1),wf_r)
             CALL coeff_sumup(rho%rho_r(2),wf_r)
             IF (globenv%ionode) OPEN(UNIT=21,FILE="ELECTRON_DENSITY.cube")
             CALL rs_pw_to_cube(wf_r%pw,21,globenv%ionode,"SPIN 1 + SPIN 2 ELECTRON DENSITY")
             IF (globenv%ionode) CLOSE(21)
             CALL coeff_copy(rho%rho_r(1),wf_r)
             CALL coeff_reduce(rho%rho_r(2),wf_r)
             IF (globenv%ionode) OPEN(UNIT=21,FILE="SPIN_DENSITY.cube")
             CALL rs_pw_to_cube(wf_r%pw,21,globenv%ionode,"SPIN 1 DENSITY - SPIN 2 DENSITY")
             IF (globenv%ionode) CLOSE(21)
             CALL pw_pool_give_back_coeff(auxbas_pw_pool,wf_r, error=error)
          ELSE
             IF (globenv%ionode) OPEN(UNIT=21,FILE="ELECTRON_DENSITY.cube")
             CALL rs_pw_to_cube(rho%rho_r(1)%pw,21,globenv%ionode,"ELECTRON DENSITY")
             IF (globenv%ionode) CLOSE(21)
          ENDIF
       ENDIF

       IF (globenv%print%v_hartree_cube) THEN
! write the hartree potential
          IF (globenv%ionode) OPEN(UNIT=21,FILE="V_HARTREE.cube")
          CALL rs_pw_to_cube(qs_env%ks_env%v_hartree_rspace%pw,&
               21,globenv%ionode,"HARTREE POTENTIAL")
          IF (globenv%ionode) CLOSE(21)
       ENDIF

       IF ((globenv%print%homos .OR. globenv%print%lumos) .AND. globenv%print%cubes) THEN
! gets the tmp grids
          CALL get_qs_env(qs_env=qs_env,pw_env=pw_env)
          CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,&
               pw_pools=pw_pools)
          CALL pw_pool_init_coeff(auxbas_pw_pool,wf_r,&
               use_data = REALDATA3D,&
               in_space = REALSPACE, error=error)
          CALL pw_pool_init_coeff(auxbas_pw_pool,wf_g,&
               use_data = COMPLEXDATA1D,&
               in_space = RECIPROCALSPACE, error=error)
          ALLOCATE(mgrid_rspace(SIZE(pw_pools)), mgrid_gspace(SIZE(pw_pools)),&
               stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          CALL pw_pools_init_coeffs(pw_pools,mgrid_gspace,&
               use_data = COMPLEXDATA1D,&
               in_space = RECIPROCALSPACE, error=error)
          CALL pw_pools_init_coeffs(pw_pools,mgrid_rspace,&
               use_data = REALDATA3D,&
               in_space = REALSPACE, error=error)
       END IF

       IF (globenv%print%homos) THEN
          DO ispin=1,dft_control%nspins
! write eigenvalues and a set of cubes
             CALL get_mo_set(mo_set=c(ispin)%mo_set,eigenvectors=mo_eigenvectors,homo=homo)
             IF (globenv%ionode) WRITE(globenv%scr,*) " "
             IF (globenv%ionode) WRITE(globenv%scr,*) " Eigenvalues of the occupied subspace spin ",ispin
             IF (globenv%ionode) WRITE(globenv%scr,*) " ---------------------------------------------"
             CALL calculate_subspace_eigenvalues(mo_eigenvectors,&
                  ks_rmpv(ispin)%matrix, &
                  para_env=globenv%para_env,ctxt=globenv%blacs_env,scr=globenv%scr, &
                  ionode=globenv%ionode)

             IF (globenv%print%cubes) THEN
                DO ivector=MAX(1,homo-globenv%print%nhomos+1),homo
                   CALL calculate_wavefunction(mo_eigenvectors,ivector,wf_r, &
                        wf_g,mgrid_rspace,&
                        mgrid_gspace,&
                        pw_env%gridlevel_info,pw_env%cube_info, &
                        pw_env%l_info,&
                        qs_env,globenv)

                   WRITE(filename,'(a13,I5.5,a1,I1.1,a5)') "WAVEFUNCTION_",ivector,"_",ispin,".cube"
                   IF (globenv%ionode) OPEN(UNIT=21,FILE=filename)
                   WRITE(filename,*) "WAVEFUNCTION ",ivector," spin ",ispin," i.e. HOMO - ",ivector-homo
                   CALL rs_pw_to_cube(wf_r%pw,21,globenv%ionode,filename)
                   IF (globenv%ionode) CLOSE(21)
                ENDDO
             END IF
          ENDDO
       ENDIF

!   *** Get localization info *** ! notice, computes this for whatever is in mo_eigenvectors
!   in particular, will do it for the eigenstates if you switch on the HOMOS keyword
       IF (globenv%print%localization)THEN
          NULLIFY(qs_loc_env)
          ALLOCATE(qs_loc_env)
          CALL qs_loc_init(qs_loc_env,qs_env)
       ENDIF

! by default we do it for the eigenstates
       IF (globenv%print%localization) THEN
          DO ispin=1,dft_control%nspins
! just get eigenstates
             CALL get_mo_set(mo_set=c(ispin)%mo_set,eigenvectors=mo_eigenvectors)
             IF (.NOT.globenv%print%homos) THEN
                IF (globenv%ionode) WRITE(globenv%scr,*) " "
                IF (globenv%ionode) WRITE(globenv%scr,*) " Eigenvalues of the occupied subspace spin ",ispin
                IF (globenv%ionode) WRITE(globenv%scr,*) " ---------------------------------------------"
                CALL calculate_subspace_eigenvalues(mo_eigenvectors,ks_rmpv(ispin)%matrix, &
                     para_env=globenv%para_env,ctxt=globenv%blacs_env,scr=globenv%scr, &
                     ionode=globenv%ionode)
             ENDIF
! do localization computation
             WRITE(filename,'(a13,I1.1)') "HOMO_CENTERS_",ispin
             IF (globenv%ionode) OPEN(UNIT=21,FILE=filename)
             CALL qs_loc_compute_set(qs_loc_env,mo_eigenvectors,globenv%ionode,21)
             IF (globenv%ionode) CLOSE(21)
          ENDDO
       ENDIF


       IF (globenv%print%lumos) THEN
          DO ispin=1,dft_control%nspins
! write eigenvalues and a set of cubes
             IF (globenv%ionode) WRITE(globenv%scr,*) " "
             IF (globenv%ionode) WRITE(globenv%scr,*) " Lowest Eigenvalues of the unoccupied subspace spin ",ispin
             IF (globenv%ionode) WRITE(globenv%scr,*) " -----------------------------------------------------"
             CALL get_mo_set(mo_set=c(ispin)%mo_set,eigenvectors=mo_eigenvectors,homo=homo,nao=nao)
             CALL cp_fm_get_info(mo_eigenvectors, nrow_global=n)
             nlumos=MAX(1,MIN(globenv%print%nlumos,nao-homo))
             CALL cp_fm_create2(new_matrix=unoccupied_orbs, &
                  nrow_block=scf_control%nrow_block, &
                  ncol_block=scf_control%ncol_block, &
                  nrow_global=n, &
                  ncol_global=nlumos, &
                  name="unoccupied_orbs", &
                  globenv=globenv)
             CALL cp_fm_init_random(unoccupied_orbs,nlumos)

             CALL ot_eigensolver(matrix_h=ks_rmpv(ispin)%matrix,matrix_s=s(1)%matrix, &
                  matrix_c=unoccupied_orbs, &
                  matrix_orthogonal_space=mo_eigenvectors, &
                  eps_gradient=scf_control%eps_lumos, &
                  preconditioner=scf_env%ot_preconditioner, &
                  iter_max=scf_control%max_iter_lumos,globenv=globenv)

             CALL calculate_subspace_eigenvalues(unoccupied_orbs,ks_rmpv(ispin)%matrix,&
                  para_env=globenv%para_env,ctxt=globenv%blacs_env,scr=globenv%scr, &
                  ionode=globenv%ionode)

             IF (globenv%print%localization) THEN
                WRITE(filename,'(a13,I1.1)') "LUMO_CENTERS_",ispin
                IF (globenv%ionode) OPEN(UNIT=21,FILE=filename)
                CALL qs_loc_compute_set(qs_loc_env,unoccupied_orbs,globenv%ionode,21)
                IF (globenv%ionode) CLOSE(21)
             ENDIF

             IF (globenv%print%cubes) THEN
                DO ivector=1,nlumos
                   CALL calculate_wavefunction(mo_vectors=unoccupied_orbs,&
                        ivector=ivector,rho=wf_r, &
                        rho_gspace=wf_g,&
                        mgrid_rspace=mgrid_rspace,&
                        mgrid_gspace=mgrid_gspace,&
                        gridlevel_info=pw_env%gridlevel_info,&
                        cube_info=pw_env%cube_info, &
                        l_info=pw_env%l_info,&
                        qs_env=qs_env,globenv=globenv)

                   WRITE(filename,'(a13,I5.5,a1,I1.1,a5)') "WAVEFUNCTION_",homo+ivector,"_",ispin,".cube"
                   IF (globenv%ionode) OPEN(UNIT=21,FILE=filename)
                   WRITE(filename,*) "WAVEFUNCTION ",homo+ivector," spin ",ispin," i.e. LUMO + ",ivector-1
                   CALL rs_pw_to_cube(wf_r%pw,21,globenv%ionode,filename)
                   IF (globenv%ionode) CLOSE(21)

                ENDDO
             ENDIF
          ENDDO
       ENDIF

       IF ((globenv%print%lumos .OR. globenv%print%homos) .AND. globenv%print%cubes) THEN
! put back the grids
          CALL pw_pool_give_back_coeff(auxbas_pw_pool,wf_r,&
               error=error)
          CALL pw_pool_give_back_coeff(auxbas_pw_pool,wf_g,&
               error=error)
          CALL pw_pools_give_back_coeffs(pw_pools,mgrid_gspace,&
               error=error)
          CALL pw_pools_give_back_coeffs(pw_pools,mgrid_rspace,&
               error=error)
          DEALLOCATE(mgrid_rspace, mgrid_gspace, stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

       END IF
       IF (globenv%print%localization) THEN
          CALL qs_loc_destroy(qs_loc_env)
          DEALLOCATE(qs_loc_env) 
       ENDIF

       IF (globenv%print%density_matrix) THEN
          CALL write_sparse_matrix(rho%rho_ao(1)%matrix,4,6,qs_env,globenv)
       END IF

       IF (Globenv%print%density_matrix_magnitude) THEN
          IF (Globenv%num_pe > 1) THEN
             WRITE (Unit=Extension,Fmt="(I6)") Globenv%mepos
             File_name = "P_magnitude"//"."//ADJUSTL(Extension)
          ELSE
             File_name = "P_magnitude"
          END IF
          CALL Open_file(File_name=File_name,&
               File_action="Write",&
               File_form="Formatted",&
               File_status="Replace",&
               Unit_number=File_unit)

          NULLIFY(Pfull) 
          CALL Replicate_matrix(S(1)%matrix,Pfull,"Full Density Matrix",.FALSE.)
          CALL Add_all_real_matrix_blocks(Pfull,Globenv)
          CALL Calculate_density_matrix(c(1)%mo_set,Pfull,error=error)
          CALL Write_blocks_maxabsval(Pfull,File_unit)
          CALL Deallocate_matrix(Pfull)
          CALL Close_file(Unit_number=File_unit)
       ENDIF

       IF (Globenv%print%kohn_sham_matrix) THEN
          CALL qs_ks_update_qs_env(qs_env%ks_env,qs_env=qs_env,&
               rebuild=.TRUE.,error=error,&
               calculate_forces=.FALSE.,just_energy=.FALSE.)
          CALL Write_sparse_matrix(ks_rmpv(1)%matrix,4,6,Qs_env,Globenv)
       END IF

    END IF
  END SUBROUTINE scf_env_post_calc
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
    TYPE(qs_scf_env_type), POINTER :: scf_env
    TYPE(qs_environment_type), POINTER :: qs_env
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='scf_env_cleanup',&
         routineP=moduleN//':'//routineN
    INTEGER :: stat, ispin,handle
    TYPE(cp_sm_pool_p_type), DIMENSION(:), POINTER :: S_sm_pools

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.
    NULLIFY(S_sm_pools)

    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL get_qs_env(qs_env,S_sm_pools=S_sm_pools,error=error)

!   *** Release SCF work storage ***

       IF (ASSOCIATED(scf_env%scf_work1)) THEN
          DO ispin=1,SIZE(scf_env%scf_work1)
             CALL cp_fm_release(scf_env%scf_work1(ispin)%matrix)
          ENDDO
          DEALLOCATE(scf_env%scf_work1)
       ENDIF
       IF (ASSOCIATED(scf_env%scf_work2)) CALL cp_fm_release(scf_env%scf_work2)
       IF (ASSOCIATED(scf_env%ortho)) CALL cp_fm_release(scf_env%ortho)

       IF (ASSOCIATED(scf_env%p_new)) THEN
          CALL sm_pools_give_back_matrix_vect(S_sm_pools,scf_env%p_new,&
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
       CASE(pao_mehod_nr)
!    call pao_opt_dealloc(pao_optimizer,pao_obj_f,qs_env,globenv,error)
       CASE(diagonalisation_method_nr)
          ;
       CASE default
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP,"unknown scf method method:"//&
               cp_to_string(scf_env%method),error,failure)
       END SELECT

    END IF
    CALL timestop(0.0_wp,handle)
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
SUBROUTINE scf_env_did_change(scf_env,qs_env,s_mstruct_changed,error)
  TYPE(qs_scf_env_type), POINTER :: scf_env
  TYPE(qs_environment_type), POINTER :: qs_env
  LOGICAL, INTENT(in) :: s_mstruct_changed
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='scf_env_did_change',&
        routineP=moduleN//':'//routineN

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     IF (ASSOCIATED(scf_env%p_new)) THEN
        CALL deallocate_matrix_set(scf_env%p_new)
     END IF
  END IF
END SUBROUTINE scf_env_did_change
!***************************************************************************

END MODULE qs_scf
