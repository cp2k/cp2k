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
  USE coefficient_types,               ONLY: coeff_deallocate,&
                                             coeff_type
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
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                             cp_fm_pool_type,&
                                             fm_pool_create,&
                                             fm_pool_release,&
                                             fm_pool_retain,&
                                             fm_pools_copy,&
                                             fm_pools_dealloc
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_retain,&
                                             cp_fm_struct_type
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
  USE particle_types,                  ONLY: particle_type
  USE pw_env_methods,                  ONLY: pw_env_create,&
                                             pw_env_rebuild
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_release,&
                                             pw_env_retain,&
                                             pw_env_type
  USE pw_grid_types,                   ONLY: HALFSPACE,&
                                             pw_grid_type
  USE pw_grids,                        ONLY: pw_grid_construct,&
                                             pw_grid_setup
  USE pw_pool_types,                   ONLY: pw_pool_init_coeff,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             COMPLEXDATA3D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_type
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
  USE qs_charges_types,                ONLY: qs_charges_create,&
                                             qs_charges_release,&
                                             qs_charges_retain,&
                                             qs_charges_type
  USE qs_collocate_density,            ONLY: calculate_rho_core
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_rho_methods,                  ONLY: qs_rho_rebuild
  USE qs_rho_types,                    ONLY: qs_rho_release,&
                                             qs_rho_type
  USE scf_control_types,               ONLY: scf_c_release,&
                                             scf_c_retain,&
                                             scf_control_type
  USE simulation_cell,                 ONLY: cell_type
  USE sparse_matrix_types,             ONLY: allocate_matrix_set,&
                                             get_matrix_info,&
                                             real_matrix_p_type,&
                                             real_matrix_type,&
                                             replicate_matrix_structure
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qs_environment_methods'

  PUBLIC :: qs_env_update_s_mstruct, qs_env_rebuild_pw_env,&
       qs_env_rebuild_aomo_pools, qs_env_rebuild_rho, qs_env_setup
!***
!****************************************************************************
CONTAINS

!!****f* qs_environment_methods/qs_env_setup [1.0] *
!!
!!   NAME
!!     qs_env_setup
!!
!!   SYNOPSIS
!!     Subroutine qs_env_setup(qs_env, error)
!!       Type(qs_environment_type), Intent (INOUT):: qs_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine qs_env_setup
!!
!!   FUNCTION
!!     initializes various components of the qs_env, that need only
!!     atomic_kind_set, cell, dft_control, scf_control, c(i)%nmo, 
!!     c(i)%nao, and particle_set to be initialized.
!!     The previous components of qs_env must be valid.
!!     Initializes pools, charges and pw_env.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - qs_env: the qs_env to set up
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     10.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE qs_env_setup(qs_env,error)
  TYPE(qs_environment_type), INTENT(inout)                :: qs_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout)            :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='qs_env_setup',&
        routineP=moduleN//':'//routineN
  TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
  TYPE(cell_type), POINTER                      :: cell
  TYPE(dft_control_type), POINTER               :: dft_control
  TYPE(scf_control_type), POINTER               :: scf_control
  TYPE(particle_type), DIMENSION(:), POINTER    :: particle_set
  TYPE(qs_charges_type), POINTER                :: qs_charges
  TYPE(cp_fm_struct_type), POINTER        :: fmstruct
  INTEGER :: n_mo, n_ao
  failure=.FALSE.
  
  CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set,&
       cell=cell, dft_control=dft_control, scf_control=scf_control,&
       particle_set=particle_set, error=error)
  CPPrecondition(ASSOCIATED(atomic_kind_set),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(cell),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(dft_control),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(scf_control),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(particle_set),cp_failure_level,routineP,error,failure)

  IF (.NOT. failure) THEN
     ! allocate qs_charges
     CALL qs_charges_create(qs_charges,nspins=dft_control%nspins,error=error)
     CALL set_qs_env(qs_env, qs_charges=qs_charges,error=error)
     CALL qs_charges_release(qs_charges,error=error)
     
     ! set up pw_env
     CALL qs_env_rebuild_pw_env(qs_env, error=error)

     ! rebuilds ao_mo_fm_pools
     CALL qs_env_rebuild_aomo_pools(qs_env, error=error)

  END IF
END SUBROUTINE qs_env_setup
!***************************************************************************


!!****f* qs_environment_types/qs_env_update_s_mstruct [1.0] *
!!
!!   NAME
!!     qs_env_update_s_mstruct
!!
!!   SYNOPSIS
!!     Subroutine qs_env_update_s_mstruct(qs_env, error)
!!       Type(qs_environment_type), Intent (INOUT):: qs_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine qs_env_update_s_mstruct
!!
!!   FUNCTION
!!     updates the s_mstruct to reflect the new overlap structure,
!!     and also updates rho_core distribution.
!!     Should be called after the atoms have moved and the new overlap
!!     has been calculated.
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
SUBROUTINE qs_env_update_s_mstruct(qs_env,error)
  TYPE(qs_environment_type), INTENT(inout)     :: qs_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='qs_env_update_s_mstruct',&
        routineP=moduleN//':'//routineN
  INTEGER :: n_atoms, stat
  INTEGER, DIMENSION(:), ALLOCATABLE :: block_sizes
  INTEGER, DIMENSION(:), POINTER :: first_row, first_col, last_row, last_col
  TYPE(cp_para_env_type), POINTER :: para_env

  failure=.FALSE.
  NULLIFY(first_row, first_col, last_row, last_col, para_env)
  
  IF (ASSOCIATED(qs_env%s_mstruct)) CALL cp_bmstruct_release(qs_env%s_mstruct)
  CALL get_qs_env(qs_env, para_env=para_env, error=error)
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
!FM  
  ! *** updates rho core ***
  CPPrecondition(ASSOCIATED(qs_env%rho_core),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(qs_env%pw_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(qs_env%pw_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL calculate_rho_core(qs_env%rho_core,qs_env%pw_env%cube_info(1),&
          qs_env%pw_env%l_info,&
          qs_env%qs_charges%total_rho_core_rspace,qs_env,error=error)
  END IF
  
END SUBROUTINE qs_env_update_s_mstruct
!***************************************************************************

!!****f* qs_environment_types/qs_env_rebuild_aomo_pools [1.0] *
!!
!!   NAME
!!     qs_env_rebuild_aomo_pools
!!
!!   SYNOPSIS
!!     Subroutine qs_env_rebuild_aomo_pools(qs_env, error)
!!       Type(qs_environment_type), Intent (INOUT), Target:: qs_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine qs_env_rebuild_aomo_pools
!!
!!   FUNCTION
!!     rebuilds the pools of the (ao x mo) full matrixes
!!     mo (qs_env%c) must contain up to date nmo and nao
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
SUBROUTINE qs_env_rebuild_aomo_pools(qs_env, error)
  TYPE(qs_environment_type), INTENT(inout), TARGET :: qs_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout)     :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='qs_env_rebuild_aomo_pools',&
        routineP=moduleN//':'//routineN
  INTEGER :: handle, nmo,nao, i, max_nmo, min_nmo, stat
  TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER :: aomo_pools
  TYPE(cp_para_env_type), POINTER                :: para_env
  TYPE(cp_blacs_env_type), POINTER               :: blacs_env
  TYPE(mo_set_p_type), DIMENSION(:), POINTER     :: mos
  TYPE(cp_fm_struct_type), POINTER         :: fmstruct
  TYPE(scf_control_type), POINTER                :: scf_control
  failure=.FALSE.
  NULLIFY(aomo_pools,para_env,blacs_env,mos, fmstruct)
  
  CALL timeset(routineN//', '//moduleN,"I","",handle)
  
  ! frees the memory before allocating the new pools
  CALL set_qs_env(qs_env, ao_mo_fm_pools=aomo_pools,error=error)

  CALL get_qs_env(qs_env,para_env=para_env,blacs_env=blacs_env,&
       c=mos,scf_control=scf_control, error=error)
  
  CPPrecondition(ASSOCIATED(mos),cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     ALLOCATE(aomo_pools(SIZE(mos)), stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF (.NOT. failure) THEN

     CALL get_mo_set(mos(1)%mo_set,nao=nao,nmo=min_nmo)
     max_nmo=min_nmo
     DO i=2,SIZE(mos)
        CALL get_mo_set(mos(i)%mo_set,nmo=nmo)
        IF (max_nmo<nmo) THEN
           CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                routineP,&
                "the mo with the most orbitals must be the first "//&
                CPSourceFileRef,&
                error=error,&
                failure=failure)
        END IF
        min_nmo=MIN(min_nmo,nmo)
     END DO
  END IF
  IF (.NOT.failure) THEN
     IF (max_nmo==min_nmo) THEN
        CALL cp_fm_struct_create(fmstruct, nrow_global=nao,&
             ncol_global=max_nmo, para_env=para_env,&
             context=blacs_env,&
             nrow_block=scf_control%nrow_block,&
             ncol_block=scf_control%ncol_block)
        CALL fm_pool_create(aomo_pools(1)%pool,fmstruct,error=error)
        CALL cp_fm_struct_release(fmstruct,error=error)
        DO i=2,SIZE(mos)
           aomo_pools(i)%pool => aomo_pools(1)%pool
           CALL fm_pool_retain(aomo_pools(1)%pool,error=error)
        END DO
     ELSE
        DO i=1,SIZE(mos)
           NULLIFY(aomo_pools(i)%pool)
           CALL get_mo_set(mos(i)%mo_set,nmo=nmo,nao=nao) 
           CALL cp_fm_struct_create(fmstruct, nrow_global=nao,&
                ncol_global=nmo, para_env=para_env,&
                context=blacs_env,&
                nrow_block=scf_control%nrow_block,&
                ncol_block=scf_control%ncol_block)
           CALL fm_pool_create(aomo_pools(i)%pool,fmstruct,error=error)
           CALL cp_fm_struct_release(fmstruct,error=error)
        END DO
     END IF
     
     CALL set_qs_env(qs_env, ao_mo_fm_pools=aomo_pools,error=error)
     CALL fm_pools_dealloc(aomo_pools,error=error)
  END IF
  
  CALL timestop(0.0_wp,handle)
END SUBROUTINE qs_env_rebuild_aomo_pools
!***************************************************************************

!!****f* qs_environment_methods/qs_env_rebuild_pw_env [1.0] *
!!
!!   NAME
!!     qs_env_rebuild_pw_env
!!
!!   SYNOPSIS
!!     Subroutine qs_env_rebuild_pw_env(qs_env, error)
!!       Type(qs_environment_type), Intent (INOUT), Target:: qs_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine qs_env_rebuild_pw_env
!!
!!   FUNCTION
!!     rebuilds the pw_env in the given qs_env, allocating it if necessary
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - qs_env: the qs_env whose pw_env has to be rebuilt
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     10.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE qs_env_rebuild_pw_env(qs_env,error)
  TYPE(qs_environment_type), INTENT(inout), TARGET :: qs_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='qs_env_rebuild_pw_env',&
        routineP=moduleN//':'//routineN
  TYPE(pw_env_type), POINTER  :: new_pw_env
  TYPE(coeff_type), POINTER   :: rho_core
  TYPE(pw_pool_type), POINTER :: auxbas_pw_pool
  INTEGER                     :: stat
  failure=.FALSE.
  
  ! rebuild pw_env
  CALL get_qs_env(qs_env, pw_env=new_pw_env)
  IF (.NOT.ASSOCIATED(new_pw_env)) THEN
     CALL pw_env_create(new_pw_env,qs_env=qs_env,error=error)
     CALL set_qs_env(qs_env,pw_env=new_pw_env)
     CALL pw_env_release(new_pw_env,error=error)
  ELSE
     CALL pw_env_rebuild(new_pw_env,qs_env=qs_env,error=error)
  END IF

  ! reallocate rho_core
  CALL get_qs_env(qs_env, pw_env=new_pw_env, rho_core=rho_core)
  CPPrecondition(ASSOCIATED(new_pw_env),cp_failure_level,routineP,error,failure)
  IF (ASSOCIATED(rho_core)) THEN
     CALL coeff_deallocate(rho_core)
     DEALLOCATE(rho_core,stat=stat)
     CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
  ELSE
     ALLOCATE(rho_core,stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF (.NOT.failure) THEN
     CALL pw_env_get(new_pw_env, auxbas_pw_pool=auxbas_pw_pool,error=error)
     CALL pw_pool_init_coeff(pool=auxbas_pw_pool, coeff=rho_core, &
          use_data=COMPLEXDATA1D, error=error)
     rho_core%pw%in_space=RECIPROCALSPACE
     CALL set_qs_env(qs_env, rho_core=rho_core,error=error)
  END IF

END SUBROUTINE qs_env_rebuild_pw_env
!***************************************************************************

!!****f* qs_environment_methods/qs_env_rebuild_rho [1.0] *
!!
!!   NAME
!!     qs_env_rebuild_rho
!!
!!   SYNOPSIS
!!     Subroutine qs_env_rebuild_rho(qs_env, rho_ao_with_h_mstruct,&
!!         rebuild_ao, rebuild_grids, error)
!!       Type(qs_environment_type), Intent (INOUT), Target:: qs_env
!!       Logical, Intent (IN), Optional:: rho_ao_with_h_mstruct,&
!!         rebuild_grids, rebuild_ao
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine qs_env_rebuild_rho
!!
!!   FUNCTION
!!     rebuilds the rho structure, making sure that everything is allocated
!!     and has the right size
!!
!!   NOTES
!!     needs updated  pw pools, s_mstruct and h.
!!     The use of p to keep the structure of h (needed for the forces)
!!     is ugly and should be removed.
!!     If necessary rho is created from scratch.
!!
!!   INPUTS
!!     - qs_env: the environment to which rho belongs
!!     - rho: the rho type to rebuild (defaults to qs_env%rho)
!!     - rho_ao_with_h_mstruct: if rho%rho_ao should have the 
!!       redundant structure of h, or just that of s.
!!       If rho is not given defaults to true, otherwise defaults to false.
!!     - rebuild_ao: if it is necessary to rebuild rho_ao. Defaults to true.
!!     - rebuild_grids: if it in necessary to rebuild rho_r and rho_g.
!!       Defaults to false.
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     10.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE qs_env_rebuild_rho(qs_env, rho_ao_with_h_mstruct,&
     rebuild_ao, rebuild_grids, error)
  TYPE(qs_environment_type), INTENT(inout), TARGET :: qs_env
  LOGICAL, INTENT(in), OPTIONAL                    :: rho_ao_with_h_mstruct,&
       rebuild_grids,rebuild_ao
  TYPE(cp_error_type), OPTIONAL, INTENT(inout)     :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='qs_env_rebuild_rho',&
        routineP=moduleN//':'//routineN

  CALL qs_rho_rebuild(qs_env%rho,qs_env=qs_env,&
       rho_ao_with_h_mstruct=rho_ao_with_h_mstruct,&
       rebuild_ao=rebuild_ao, rebuild_grids=rebuild_grids, error=error)
END SUBROUTINE qs_env_rebuild_rho
!***************************************************************************


!!****f* qs_environment_methods/qs_env_update_rho_core [1.0] *
!!
!!   NAME
!!     qs_env_update_rho_core
!!
!!   FUNCTION
!!     Updates rho_core (the grid with the core charge)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - qs_env: the environement whose rho core has to be updated
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     10.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE qs_env_update_rho_core(qs_env,error)
  TYPE(qs_environment_type), INTENT(inout), TARGET :: qs_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='qs_env_update_rho_core',&
        routineP=moduleN//':'//routineN
  TYPE(pw_env_type), POINTER :: pw_env
  TYPE(coeff_type), POINTER  :: rho_core
  TYPE(qs_charges_type), POINTER :: qs_charges
  failure=.FALSE.

  CALL get_qs_env(qs_env, pw_env=pw_env, rho_core=rho_core,&
       qs_charges=qs_charges,error=error)
  CPPrecondition(ASSOCIATED(pw_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(rho_core),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(qs_charges),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL calculate_rho_core(rho_core,pw_env%cube_info(1),&
          pw_env%l_info,&
          qs_charges%total_rho_core_rspace,qs_env,error)
  END IF
END SUBROUTINE qs_env_update_rho_core
!***************************************************************************

END MODULE qs_environment_methods
