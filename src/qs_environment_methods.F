!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/qs_environment_methods [1.0] *
!!
!!   NAME
!!     qs_environment_methods
!!
!!   FUNCTION
!!     qs_environement methods that use many other modules
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     09.2002 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE qs_environment_methods
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE cp_b_matrix_structure,           ONLY: cp_b_matrix_struct_type,&
                                             cp_bmstruct_create,&
                                             cp_bmstruct_release
  USE cp_block_matrix,                 ONLY: cp_block_matrix_p_type,&
                                             cp_block_matrix_type,&
                                             cp_sp_array_create,&
                                             cp_sp_array_dealloc,&
                                             cp_sp_create,&
                                             cp_sp_multiplies_blacs,&
                                             cp_sp_set_to
  USE cp_error_handling,               ONLY: cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error,&
                                             cp_unimplemented_error
  USE cp_f_matrix_struct,              ONLY: cp_f_matrix_struct_type,&
                                             cp_fmstruct_create,&
                                             cp_fmstruct_release,&
                                             cp_fmstruct_retain
  USE cp_fm_pool,                      ONLY: cp_fm_pool_type,&
                                             fm_pool_create,&
                                             fm_pool_release,&
                                             fm_pool_retain
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_matrix_utils,                 ONLY: cp_sm_set_dist_sparsity_from_sm
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE dft_types,                       ONLY: dft_control_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl,&
                                             wp=>dp
  USE pw_grid_types,                   ONLY: HALFSPACE,&
                                             pw_grid_type
  USE pw_grids,                        ONLY: pw_grid_construct,&
                                             pw_grid_setup
  USE pw_pool,                         ONLY: pw_pool_create,&
                                             pw_pool_p_type,&
                                             pw_pool_release,&
                                             pw_pool_retain,&
                                             pw_pool_type,&
                                             pw_pools_copy,&
                                             pw_pools_dealloc
  USE qs_blacs,                        ONLY: cp_fm_add,&
                                             cp_fm_create2,&
                                             cp_fm_gemm,&
                                             cp_fm_get_info,&
                                             cp_fm_release,&
                                             cp_fm_vect_create2,&
                                             cp_fm_vect_dealloc,&
                                             cp_full_matrix_p_type,&
                                             cp_full_matrix_type,&
                                             cp_sm_fm_multiply,&
                                             optimal_blacs_col_block_size,&
                                             optimal_blacs_row_block_size
  USE qs_build_ks_matrix,              ONLY: qs_ks_did_change,&
                                             qs_ks_env_type,&
                                             qs_ks_update_qs_env
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_mo_types,                     ONLY: get_mo_set
  USE scf_control_types,               ONLY: scf_c_release,&
                                             scf_c_retain,&
                                             scf_control_type
  USE simulation_cell,                 ONLY: cell_type
  USE sparse_matrix_types,             ONLY: get_matrix_info,&
                                             real_matrix_p_type,&
                                             real_matrix_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qs_environment_methods'

  PUBLIC :: qs_env_update_s_mstruct, qs_env_rebuild_pw_pools,&
    qs_env_rebuild_aomo_pool
!***
!****************************************************************************
CONTAINS
!!****f* qs_environment_types/qs_env_update_s_mstruct [1.0] *
!!
!!   NAME
!!     qs_env_update_s_mstruct
!!
!!   FUNCTION
!!     updates the s_mstruct to reflect the new overlap structure
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - qs_env: the environment to update
!!     - para_env: the mpi environment
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE qs_env_update_s_mstruct(qs_env,para_env,error)
  TYPE(qs_environment_type), INTENT(inout)     :: qs_env
  TYPE(cp_para_env_type), POINTER              :: para_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='qs_env_update_s_mstruct',&
        routineP=moduleN//':'//routineN
  INTEGER :: n_atoms, stat
  INTEGER, DIMENSION(:), ALLOCATABLE :: block_sizes
  INTEGER, DIMENSION(:), POINTER :: first_row, first_col, last_row, last_col
  failure=.FALSE.
  
  IF (ASSOCIATED(qs_env%s_mstruct)) CALL cp_bmstruct_release(qs_env%s_mstruct)
  n_atoms=SIZE(qs_env%particle_set)
  ALLOCATE(block_sizes(n_atoms),stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  block_sizes=qs_env%particle_set(:)%last_sgf-qs_env%particle_set(:)%first_sgf+1
  CALL cp_bmstruct_create(qs_env%s_mstruct,&
       symmetric=.TRUE.,&
       para_env=para_env,&
       b_row_sizes=block_sizes,&
       b_col_sizes=block_sizes,&
       error=error)
  CALL get_matrix_info(qs_env%s(1)%matrix,first_row=first_row,&
       last_row=last_row, first_col=first_col,&
       last_col=last_col)

  CALL cp_assert(ALL(qs_env%s_mstruct%b_dims%block_begins_at_row(1:n_atoms)==&
       first_row),cp_failure_level,cp_assertion_failed,routineP,&
       "s has unexpected size, "//&
       CPSourceFileRef,&
       error)
  CALL cp_assert(ALL(qs_env%s_mstruct%b_dims%block_begins_at_col(1:n_atoms)==&
       first_col),cp_failure_level,cp_assertion_failed,routineP,&
       "s has unexpected size, "//&
       CPSourceFileRef,&
       error)
  CALL cp_assert(qs_env%s_mstruct%b_dims%block_begins_at_col(n_atoms+1)==&
       last_col(n_atoms)+1,cp_failure_level,cp_assertion_failed,routineP,&
       "s has unexpected size, "//&
       CPSourceFileRef,&
       error)
  CALL cp_assert(qs_env%s_mstruct%b_dims%block_begins_at_row(n_atoms+1)==&
       last_row(n_atoms)+1,cp_failure_level,cp_assertion_failed,routineP,&
       "s has unexpected size, "//&
       CPSourceFileRef,&
       error)

  ! at the moment read from s, maybe in the future s will take its structure
  ! from this
  CALL cp_sm_set_dist_sparsity_from_sm(qs_env%s_mstruct,&
       para_env=para_env,&
       sample_matrix=qs_env%s(1)%matrix,&
       error=error)
END SUBROUTINE qs_env_update_s_mstruct
!***************************************************************************

!!****f* qs_environment_types/qs_env_rebuild_pw_pools [1.0] *
!!
!!   NAME
!!     qs_env_rebuild_pw_pools
!!
!!   FUNCTION
!!     rebuilds the pw pools
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - qs_env: the qs environment where to build the pools
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE qs_env_rebuild_pw_pools(qs_env, error)
  TYPE(qs_environment_type), INTENT(inout), TARGET :: qs_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout)     :: error

  LOGICAL :: failure, pw_grid_information
  CHARACTER(len=*), PARAMETER :: routineN='qs_env_rebuild_pw_pools',&
       routineP=moduleN//':'//routineN
  REAL(wp), DIMENSION(:), POINTER :: cutoff
  INTEGER :: handle,i,igrid_level,stat, ngrid_level
  TYPE(cell_type), POINTER          :: cell
  TYPE(pw_grid_type), POINTER :: pw_grid, old_pw_grid
  TYPE(dft_control_type), POINTER :: dft_control
  TYPE(pw_pool_p_type),DIMENSION(:), POINTER :: pw_pools
  TYPE(cp_para_env_type), POINTER :: para_env
  TYPE(cp_logger_type), POINTER :: logger
  failure=.FALSE.
  NULLIFY(cutoff,cell,pw_grid,old_pw_grid,dft_control,pw_pools,para_env)

  CALL timeset(routineN//', '//moduleN,"I","",handle)
  CALL get_qs_env(qs_env, para_env=para_env, dft_control=dft_control,&
       cell=cell)
  cutoff => dft_control%qs_control%e_cutoff
  ngrid_level=SIZE(cutoff)
  ALLOCATE(pw_pools(ngrid_level),stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  logger=>cp_error_get_logger(error)
  pw_grid_information=logger%print_keys%pw_grid_information

  DO igrid_level=1,ngrid_level
     ALLOCATE(pw_grid,stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     CALL pw_grid_construct(pw_grid)
     pw_grid%grid_span = HALFSPACE
     pw_grid%para%rs_dims ( 1 ) = para_env%num_pe
     pw_grid%para%rs_dims ( 2 ) = 1

     IF (igrid_level == 1) THEN
        IF (pw_grid_information) THEN
           CALL pw_grid_setup(cell,pw_grid,cutoff(igrid_level),&
                pe_group=para_env%group,&
                info= cp_logger_get_default_unit_nr(logger,&
                local=.FALSE.))
        ELSE
           CALL pw_grid_setup(cell,pw_grid,cutoff(igrid_level),&
                pe_group=para_env%group)
        END IF
     ELSE
        IF (pw_grid_information) THEN
           CALL pw_grid_setup(cell,pw_grid,cutoff(igrid_level),&
                pe_group=para_env%group,&
                info=cp_logger_get_default_unit_nr(logger),&
                ref_grid=old_pw_grid)
        ELSE
           CALL pw_grid_setup(cell,pw_grid,cutoff(igrid_level),&
                pe_group=para_env%group,&
                ref_grid=old_pw_grid)
        END IF
     END IF
     old_pw_grid => pw_grid

     NULLIFY(pw_pools(igrid_level)%pool)
     CALL pw_pool_create(pw_pools(igrid_level)%pool,&
          pw_grid_ptr=pw_grid,error=error)
  END DO

  CALL set_qs_env(qs_env, pw_pools=pw_pools)
  CALL pw_pools_dealloc(pw_pools,error=error)

  CALL timestop(0.0_wp,handle)
END SUBROUTINE qs_env_rebuild_pw_pools
!***************************************************************************

!!****f* qs_environment_types/qs_env_rebuild_fm_pools [1.0] *
!!
!!   NAME
!!     qs_env_rebuild_fm_pools
!!
!!   FUNCTION
!!     rebuilds the pools of the full matrixes
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - qs_env: the environement where the pools should be rebuilt
!!     - fmstruct: the matrix structure of the ao x mo matrixes
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE qs_env_rebuild_aomo_pool(qs_env, fmstruct, error)
  TYPE(qs_environment_type), INTENT(inout), TARGET :: qs_env
  TYPE(cp_f_matrix_struct_type), POINTER           :: fmstruct
  TYPE(cp_error_type), OPTIONAL, INTENT(inout)     :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='qs_env_rebuild_aomo_pool',&
        routineP=moduleN//':'//routineN
  INTEGER :: handle
  TYPE(cp_fm_pool_type), POINTER          :: aomo_pool
  failure=.FALSE.
  NULLIFY(aomo_pool)
  
  CALL timeset(routineN//', '//moduleN,"I","",handle)
  
  CALL fm_pool_create(aomo_pool,fmstruct,error=error)
  CALL set_qs_env(qs_env, ao_mo_fm_pool=aomo_pool, error=error)
  CALL fm_pool_release(aomo_pool,error=error)

  CALL timestop(0.0_wp,handle)
END SUBROUTINE qs_env_rebuild_aomo_pool
!***************************************************************************

END MODULE qs_environment_methods
