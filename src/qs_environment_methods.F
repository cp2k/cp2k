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
!!     - local atom distribution (25.06.2003,MK)
!!
!!   SOURCE
!****************************************************************************

MODULE qs_environment_methods

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE coefficient_types,               ONLY: coeff_deallocate,&
                                             coeff_type
  USE cp_array_i_utils,                ONLY: cp_1d_i_p_type
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
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm,&
                                             cp_fm_scale_and_add
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                             cp_fm_release,&
                                             cp_fm_p_type,&
                                             cp_fm_type
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
  USE cp_local_lists_types,            ONLY: cp_ll_build_global_dist,&
                                             cp_ll_create,&
                                             cp_ll_destroy_global_dist,&
                                             cp_ll_release,&
                                             cp_ll_retain,&
                                             cp_local_lists_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_matrix_utils,                 ONLY: sm_struct_add_sm_local_blocks
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE cp_sm_pool_types,                ONLY: cp_sm_pool_p_type,&
                                             cp_sm_pool_type,&
                                             sm_pool_create,&
                                             sm_pool_get_mstruct,&
                                             sm_pool_release,&
                                             sm_pool_retain,&
                                             sm_pools_dealloc,&
                                             sm_pools_flush_cache
  USE cp_sm_struct,                    ONLY: cp_sm_struct_type,&
                                             sm_struct_clear_local_blocks,&
                                             sm_struct_create,&
                                             sm_struct_release,&
                                             sm_struct_retain,&
                                             sm_struct_write
  USE distribution_2d_types,           ONLY: distribution_2d_type, &
                                             distribution_2d_create
  USE dft_types,                       ONLY: dft_control_type
  USE global_types,                    ONLY: LOW,global_environment_type
  USE kinds,                           ONLY: int_size,&
                                             wp => dp,&
                                             wp_size => dp_size
  USE machine,                         ONLY: m_flush
  USE message_passing,                 ONLY: mp_sum,&
                                             mp_sync
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
  USE qs_charges_types,                ONLY: qs_charges_create,&
                                             qs_charges_release,&
                                             qs_charges_retain,&
                                             qs_charges_type
  USE qs_collocate_density,            ONLY: calculate_rho_core
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_ks_methods,                   ONLY: qs_ks_did_change
  USE qs_matrix_pools,                 ONLY: mpools_rebuild_fm_pools,&
                                             mpools_rebuild_s_sm_pools,&
                                             qs_matrix_pools_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_rho_methods,                  ONLY: qs_rho_rebuild
  USE qs_rho_types,                    ONLY: qs_rho_release,&
                                             qs_rho_type
  USE qs_scf,                          ONLY: scf_env_did_change
  USE scf_control_types,               ONLY: scf_c_release,&
                                             scf_c_retain,&
                                             scf_control_type
  USE simulation_cell,                 ONLY: cell_type
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qs_environment_methods'

  PUBLIC :: qs_env_update_s_mstruct,&
            qs_env_rebuild_pw_env,&
            qs_env_rebuild_rho,&
            qs_env_setup
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
SUBROUTINE qs_env_setup(qs_env,error,globenv)

  TYPE(global_environment_type), INTENT(IN)               :: globenv
  TYPE(qs_environment_type), INTENT(inout)                :: qs_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout)            :: error

  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='qs_env_setup',&
        routineP=moduleN//':'//routineN
  TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
  TYPE(dft_control_type), POINTER               :: dft_control
  TYPE(scf_control_type), POINTER               :: scf_control
  TYPE(qs_charges_type), POINTER                :: qs_charges
  TYPE(cp_para_env_type), POINTER               :: para_env
  TYPE(distribution_2d_type), POINTER           :: distribution_2d
  INTEGER :: n_mo, n_ao

  failure=.FALSE.
  NULLIFY(atomic_kind_set, dft_control, scf_control, qs_charges, para_env)

  CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set,&
       dft_control=dft_control, scf_control=scf_control,&
       para_env=para_env, error=error)
  CPPrecondition(ASSOCIATED(atomic_kind_set),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(dft_control),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(scf_control),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     ! allocate qs_charges
     CALL qs_charges_create(qs_charges,nspins=dft_control%nspins,error=error)
     CALL set_qs_env(qs_env, qs_charges=qs_charges,error=error)
     CALL qs_charges_release(qs_charges,error=error)

     ! set up pw_env
     CALL qs_env_rebuild_pw_env(qs_env, error=error)

     ! rebuilds fm_pools
     CALL mpools_rebuild_fm_pools(qs_env%mpools,mos=qs_env%c,&
          blacs_env=qs_env%blacs_env, para_env=para_env,&
          nrow_block=scf_control%nrow_block,&
          ncol_block=scf_control%ncol_block,&
          error=error)

     ! create local_particles
     CALL rr_local_particles_create(atomic_kind_set,qs_env%local_rows,&
                                    qs_env%local_cols,qs_env%distribution_2d,&
                                    qs_env%blacs_env,globenv)
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
!!       Type(qs_environment_type), Pointer:: qs_env
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
  TYPE(qs_environment_type), POINTER           :: qs_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='qs_env_update_s_mstruct',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat, handle
  TYPE(cp_para_env_type), POINTER :: para_env
  TYPE(dft_control_type), POINTER :: dft_control
  TYPE(cp_logger_type), POINTER :: logger

  CALL timeset(routineN,"I","",handle)
  failure=.FALSE.
  NULLIFY(para_env, dft_control)
  logger => cp_error_get_logger(error)

  CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)

  CALL get_qs_env(qs_env, para_env=para_env, dft_control=dft_control,&
       error=error)

  ! *** rebuilds the S_sm_pools ***
  CALL mpools_rebuild_s_sm_pools(qs_env%mpools,s=qs_env%s,&
       nspins=dft_control%nspins,particle_set=qs_env%particle_set,&
       para_env=para_env, error=error)

  ! *** updates rho core ***
  CPPrecondition(ASSOCIATED(qs_env%rho_core),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(qs_env%pw_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(qs_env%pw_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL calculate_rho_core(qs_env%rho_core,qs_env%pw_env%cube_info(1),&
          qs_env%pw_env%l_info,&
          qs_env%qs_charges%total_rho_core_rspace,qs_env,error=error)
  END IF

  ! *** tell ks_env ***
  IF (ASSOCIATED(qs_env%ks_env)) THEN
     CALL qs_ks_did_change(qs_env%ks_env,qs_env=qs_env,s_mstruct_changed=.TRUE.,&
          error=error)
  END IF

  !   *** Updates rho structure ***
  CALL qs_env_rebuild_rho(qs_env=qs_env)

  ! *** tell scf_env ***
  IF (ASSOCIATED(qs_env%scf_env)) THEN
     CALL scf_env_did_change(qs_env%scf_env,qs_env=qs_env,&
          s_mstruct_changed=.TRUE.,&
          error=error)
  END IF

  CALL timestop(0.0_wp,handle)

END SUBROUTINE qs_env_update_s_mstruct
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
!!     - qs_env: the environment in which rho should be rebuilt
!!     - rho_ao_with_h_mstruct: if rho%rho_ao should have the
!!       redundant structure of h, or just that of s.
!!       Defaults to false.
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

!!****f* qs_environment_methods/rr_local_particles_create [1.0] *
!!
!!   NAME
!!     rr_local_particles_create
!!
!!   FUNCTION
!!     distributes the particles
!!
!!   NOTES
!!     - presumably not the right thing, as they are distributed over a row/a col of cpus
!!       i.e. global lists based on this info are wrong
!!
!!   ARGUMENTS
!!     - local_rows: local atomic block rows (equivalent to local particles)
!!     - local_cols: local atomic block columns
!!     - atomic_kind_set: the atomic kind set to distribute
!!     - para_env: the parallel environement at the basis of the distribution
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

  SUBROUTINE rr_local_particles_create(atomic_kind_set,local_rows,local_cols,&
                                       distribution_2d,blacs_env,globenv)

    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(cp_blacs_env_type), POINTER              :: blacs_env
    TYPE(distribution_2d_type), POINTER           :: distribution_2d
    TYPE(cp_local_lists_type), POINTER            :: local_rows,local_cols
    TYPE(global_environment_type), INTENT(IN)     :: globenv

!   *** Local types ***

    TYPE local_atom_type
      INTEGER, DIMENSION(:), POINTER :: list_local_col,list_local_row
    END TYPE local_atom_type

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "rr_local_particles_create"
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE "//routine_name//" (MODULE "//moduleN//")"

!   *** Local variables ***

    TYPE(atomic_kind_type), POINTER   :: atomic_kind
    TYPE(gto_basis_set_type), POINTER :: orb_basis_set

    INTEGER :: group,i,iatom,iatom_local,ipcol,iprow,ikind,istat,j,&
               mypcol,mype,myprow,natom,nkind,npcol,npe,nprow,nsgf,&
               output_unit
    LOGICAL :: ionode

    TYPE(local_atom_type), DIMENSION(:), ALLOCATABLE :: local_atom
    INTEGER, DIMENSION(:), ALLOCATABLE         :: natom_local_col,&
                                                  natom_local_row,&
                                                  ncol,ncol_old,&
                                                  nparticle,&
                                                  nrow,nrow_old
    INTEGER, DIMENSION(:), POINTER             :: atom_list
    INTEGER, DIMENSION(:), POINTER             :: row_distribution,col_distribution
    INTEGER                                    :: natom_total

!   ---------------------------------------------------------------------------

    group = globenv%group
    ionode = globenv%ionode
    mype = globenv%mepos
    npe = globenv%num_pe
    output_unit = globenv%scr

    myprow = blacs_env%mepos(1) + 1
    mypcol = blacs_env%mepos(2) + 1

    nprow = blacs_env%num_pe(1)
    npcol = blacs_env%num_pe(2)

    ALLOCATE (nrow(nprow),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"nrow",nprow*int_size)
    nrow(:) = 0

    ALLOCATE (nrow_old(nprow),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"nrow_old",nprow*int_size)

    ALLOCATE (ncol(npcol),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"ncol",npcol*int_size)
    ncol(:) = 0

    ALLOCATE (ncol_old(npcol),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"ncol_old",npcol*int_size)

    nkind = SIZE(atomic_kind_set)

    ALLOCATE (local_atom(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"local_atom",nkind*int_size)

    ALLOCATE (natom_local_row(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"natom_local_row",nkind*int_size)

    ALLOCATE (natom_local_col(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"natom_local_col",nkind*int_size)

    natom_total = 0
    DO ikind=1,nkind
      atomic_kind => atomic_kind_set(ikind)
      NULLIFY (atom_list)
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           atom_list=atom_list)
      natom_total = natom_total + SIZE(atom_list)
    ENDDO
    ALLOCATE(row_distribution(natom_total)) 
    ALLOCATE(col_distribution(natom_total)) 

    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)

      NULLIFY (atom_list)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           atom_list=atom_list,&
                           orb_basis_set=orb_basis_set)

      natom = SIZE(atom_list)

      IF (ASSOCIATED(orb_basis_set)) THEN
        CALL get_gto_basis_set(gto_basis_set=orb_basis_set,nsgf=nsgf)
      ELSE
        nsgf = 1
      END IF

!     *** Save current estimated work load ***

      nrow_old(:) = nrow(:)
      ncol_old(:) = ncol(:)

!     *** Get the number of local particles for each process row ***

      natom_local_row(ikind) = 0

      DO iatom=1,natom
        i = 1
        DO j=1,nprow
          IF (nrow(j) < nrow(i)) i = j
        END DO
        nrow(i) = nrow(i) + nsgf
        IF (i == myprow) natom_local_row(ikind) = natom_local_row(ikind) + 1
      END DO

!     *** Get the number of local particles for each process column ***

      natom_local_col(ikind) = 0

      DO iatom=1,natom
        i = 1
        DO j=1,npcol
          IF (ncol(j) < ncol(i)) i = j
        END DO
        ncol(i) = ncol(i) + nsgf
        IF (i == mypcol) natom_local_col(ikind) = natom_local_col(ikind) + 1
      END DO

!     *** Reset work load ***

      nrow(:) = nrow_old(:)
      ncol(:) = ncol_old(:)

!     *** Distribute the particles ***

!     *** Atomic block rows ***

      IF (natom_local_row(ikind) > 0) THEN
        ALLOCATE (local_atom(ikind)%list_local_row(natom_local_row(ikind)),STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"local_atom(ikind)%list_local_row",&
                                         natom_local_row(ikind)*int_size)
      ELSE
        NULLIFY (local_atom(ikind)%list_local_row)
      END IF

      iatom_local = 0

      DO iatom=1,natom
        i = 1
        DO j=1,nprow
          IF (nrow(j) < nrow(i)) i = j
        END DO
        nrow(i) = nrow(i) + nsgf
        row_distribution(atom_list(iatom)) = i - 1 ! generate blacs rows
        IF (i == myprow) THEN
          iatom_local = iatom_local + 1
          local_atom(ikind)%list_local_row(iatom_local) = atom_list(iatom)
        END IF
      END DO

!     *** Atomic block columns ***

      IF (natom_local_col(ikind) > 0) THEN
        ALLOCATE (local_atom(ikind)%list_local_col(natom_local_col(ikind)),STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"local_atom(ikind)%list_local_col",&
                                         natom_local_col(ikind)*int_size)
      ELSE
        NULLIFY (local_atom(ikind)%list_local_col)
      END IF

      iatom_local = 0

      DO iatom=1,natom
        i = 1
        DO j=1,npcol
          IF (ncol(j) < ncol(i)) i = j
        END DO
        ncol(i) = ncol(i) + nsgf
        col_distribution(atom_list(iatom)) = i - 1 ! generate blacs cols
        IF (i == mypcol) THEN
          iatom_local = iatom_local + 1
          local_atom(ikind)%list_local_col(iatom_local) = atom_list(iatom)
        END IF
      END DO

    END DO

!   *** generate the distribution_2d structure
    CALL distribution_2d_create(distribution_2d, row_distribution, &
                                col_distribution, blacs_env)

!   *** Create the local particle structure for block rows and columns ***

    CALL cp_ll_create(local_rows,&
                      n_el=natom_local_row,&
                      para_env=globenv%para_env)

    CALL cp_ll_create(local_cols,&
                      n_el=natom_local_col,&
                      para_env=globenv%para_env)

!   *** Store the generated distribution ***

    DO ikind=1,nkind
      local_rows%list(ikind)%array(:) = local_atom(ikind)%list_local_row(:)
      local_cols%list(ikind)%array(:) = local_atom(ikind)%list_local_col(:)
    END DO

!   *** Print distribution, if requested ***

    IF (globenv%print%distribution) THEN

!     *** Print row distribution ***

      ALLOCATE (nparticle(nprow),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"nparticle",nprow*int_size)
      nparticle(:) = 0

      IF (mypcol == 1) nparticle(myprow) = SUM(natom_local_row)

      CALL mp_sum(nparticle,group)

      IF (ionode) THEN
        WRITE (UNIT=output_unit,&
               FMT="(/,/,T2,A,/,/,T3,A,/,/,(T4,I10,T26,I10,T53,I10))")&
          "DISTRIBUTION OF THE PARTICLES (ROWS)",&
          "Process row   Number of particles      Number of matrix rows",&
          (iprow-1,nparticle(iprow),nrow(iprow),iprow=1,nprow)
        WRITE (UNIT=output_unit,FMT="(/,T11,A3,T26,I10,T53,I10)")&
          "Sum",SUM(nparticle),SUM(nrow)
      END IF

      CALL m_flush(output_unit)
      CALL mp_sync(group)

      IF (globenv%print%level > LOW) THEN
        DO iprow=1,nprow
          IF ((iprow == myprow).AND.(mypcol == 1)) THEN
            IF (iprow == 1) THEN
              WRITE (UNIT=output_unit,FMT="(/,T3,A)")&
                "Process row   Kind   Local particles (global indices)"
            END IF
            DO ikind=1,nkind
              IF (ikind == 1) THEN
                WRITE (UNIT=output_unit,FMT="(/,T4,I10,2X,I5,(T24,10I6))")&
                  iprow-1,ikind,(local_atom(ikind)%list_local_row(iatom_local),&
                                 iatom_local=1,natom_local_row(ikind))
              ELSE
                WRITE (UNIT=output_unit,FMT="(T16,I5,(T24,10I6))")&
                  ikind,(local_atom(ikind)%list_local_row(iatom_local),&
                         iatom_local=1,natom_local_row(ikind))
              END IF
            END DO
          END IF
          CALL m_flush(output_unit)
          CALL mp_sync(group)
        END DO
      END IF

      DEALLOCATE (nparticle,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"nparticle")

!     *** Print column distribution ***

      ALLOCATE (nparticle(npcol),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"nparticle",npcol*int_size)
      nparticle(:) = 0

      IF (myprow == 1) nparticle(mypcol) = SUM(natom_local_col)

      CALL mp_sum(nparticle,group)

      IF (ionode) THEN
        WRITE (UNIT=output_unit,&
               FMT="(/,/,T2,A,/,/,T3,A,/,/,(T4,I10,T26,I10,T53,I10))")&
          "DISTRIBUTION OF THE PARTICLES (COLUMNS)",&
          "Process col   Number of particles   Number of matrix columns",&
          (ipcol-1,nparticle(ipcol),ncol(ipcol),ipcol=1,npcol)
        WRITE (UNIT=output_unit,FMT="(/,T11,A3,T26,I10,T53,I10)")&
          "Sum",SUM(nparticle),SUM(ncol)
      END IF

      CALL m_flush(output_unit)
      CALL mp_sync(group)

      IF (globenv%print%level > LOW) THEN
        DO ipcol=1,npcol
          IF ((myprow == 1).AND.(ipcol == mypcol)) THEN
            IF (ipcol == 1) THEN
              WRITE (UNIT=output_unit,FMT="(/,T3,A)")&
                "Process col   Kind   Local particles (global indices)"
            END IF
            DO ikind=1,nkind
              IF (ikind == 1) THEN
                WRITE (UNIT=output_unit,FMT="(/,T4,I10,2X,I5,(T24,10I6))")&
                  ipcol-1,ikind,(local_atom(ikind)%list_local_col(iatom_local),&
                                 iatom_local=1,natom_local_col(ikind))
              ELSE
                WRITE (UNIT=output_unit,FMT="(T16,I5,(T24,10I6))")&
                  ikind,(local_atom(ikind)%list_local_col(iatom_local),&
                         iatom_local=1,natom_local_col(ikind))
              END IF
            END DO
          END IF
          CALL m_flush(output_unit)
          CALL mp_sync(group)
        END DO
      END IF

      DEALLOCATE (nparticle,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"nparticle")

    END IF

!   *** Release work storage ***

    DEALLOCATE (nrow,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"nrow")

    DEALLOCATE (nrow_old,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"nrow_old")

    DEALLOCATE (ncol,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"ncol")

    DEALLOCATE (ncol_old,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"ncol_old")

    DO ikind=1,nkind
      IF (ASSOCIATED(local_atom(ikind)%list_local_row)) THEN
        DEALLOCATE (local_atom(ikind)%list_local_row,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"local_atom(ikind)%list_local_row")
      END IF
      IF (ASSOCIATED(local_atom(ikind)%list_local_col)) THEN
        DEALLOCATE (local_atom(ikind)%list_local_col,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"local_atom(ikind)%list_local_col")
      END IF
    END DO

    DEALLOCATE (local_atom,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"atom")

    DEALLOCATE (natom_local_row,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"natom_local_row")

    DEALLOCATE (natom_local_col,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"natom_local_col")

    DEALLOCATE (row_distribution,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"row_distribution")

    DEALLOCATE (col_distribution,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"col_distribution")

  END SUBROUTINE rr_local_particles_create

! *****************************************************************************

END MODULE qs_environment_methods
