!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!


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

  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE coefficient_types,               ONLY: coeff_deallocate,&
                                             coeff_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_ddapc_types,                  ONLY: cp_ddapc_release
  USE cp_ddapc_util,                   ONLY: cp_ddapc_init
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_2d_types,           ONLY: distribution_2d_release,&
                                             distribution_2d_type
  USE distribution_methods,            ONLY: distribute_molecules_2d
  USE kinds,                           ONLY: dp
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE particle_types,                  ONLY: particle_type
  USE pw_env_methods,                  ONLY: pw_env_create,&
                                             pw_env_rebuild
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_release,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_init_coeff,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             RECIPROCALSPACE
  USE qs_charges_types,                ONLY: qs_charges_create,&
                                             qs_charges_release,&
                                             qs_charges_type
  USE qs_collocate_density,            ONLY: calculate_rho_core
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_ks_methods,                   ONLY: qs_ks_did_change
  USE qs_matrix_pools,                 ONLY: mpools_rebuild_fm_pools,&
                                             mpools_rebuild_s_sm_pools
  USE qs_outer_scf,                    ONLY: outer_loop_variables_count
  USE qs_rho_methods,                  ONLY: qs_rho_rebuild
  USE qs_scf_types,                    ONLY: scf_env_did_change
  USE scf_control_types,               ONLY: scf_control_type
  USE sparse_matrix_types,             ONLY: deallocate_matrix_set,&
                                             real_matrix_p_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_environment_methods'

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
SUBROUTINE qs_env_setup(qs_env,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_env_setup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nhistory, nvariables, &
                                                stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: outer_scf_history
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_charges_type), POINTER           :: qs_charges
    TYPE(scf_control_type), POINTER          :: scf_control

    CALL timeset(routineN,handle)

    failure = .FALSE.

    NULLIFY(atomic_kind_set, dft_control, scf_control, qs_charges, para_env,&
         distribution_2d, molecule_kind_set, molecule_set, particle_set)

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    dft_control=dft_control,&
                    molecule_kind_set=molecule_kind_set,&
                    molecule_set=molecule_set,&
                    particle_set=particle_set,&
                    scf_control=scf_control,&
                    para_env=para_env,&
                    error=error)
    CPPrecondition(ASSOCIATED(atomic_kind_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(dft_control),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(scf_control),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       ! allocate qs_charges
       CALL qs_charges_create(qs_charges,nspins=dft_control%nspins,error=error)
       CALL set_qs_env(qs_env, qs_charges=qs_charges,error=error)
       CALL qs_charges_release(qs_charges,error=error)

       ! outer scf setup
       IF (scf_control%outer_scf%have_scf) THEN
          nvariables=outer_loop_variables_count(scf_control,error=error)
          nhistory=scf_control%outer_scf%extrapolation_order
          ALLOCATE(outer_scf_history(nvariables,nhistory),stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          CALL set_qs_env(qs_env,outer_scf_history=outer_scf_history,error=error)
          CALL set_qs_env(qs_env,outer_scf_ihistory=0,error=error)
       ENDIF

       ! set up pw_env
       CALL qs_env_rebuild_pw_env(qs_env, error=error)

       ! rebuilds fm_pools

       CALL mpools_rebuild_fm_pools(qs_env%mpools,mos=qs_env%mos,&
              blacs_env=qs_env%blacs_env, para_env=para_env,&
              nrow_block=scf_control%nrow_block,&
              ncol_block=scf_control%ncol_block,&
              error=error)

       ! create 2d distribution

       CALL distribute_molecules_2d(particle_kind_set=atomic_kind_set,&
                                    particle_set=particle_set,&
                                    molecule_kind_set=molecule_kind_set,&
                                    molecule_set=molecule_set,&
                                    distribution_2d=distribution_2d,&
                                    blacs_env=qs_env%blacs_env,&
                                    force_env_section=qs_env%input, error=error)

       CALL set_qs_env(qs_env,distribution_2d=distribution_2d,error=error)
       CALL distribution_2d_release(distribution_2d, error=error)
    END IF

    CALL timestop(handle)

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
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_env_update_s_mstruct', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks

    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.
    NULLIFY(para_env, dft_control,particle_set)

    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)

    CALL get_qs_env(qs_env, para_env=para_env, dft_control=dft_control,&
         particle_set=particle_set, matrix_ks=matrix_ks, error=error)

    ! *** rebuilds the S_sm_pools ***
    CALL mpools_rebuild_s_sm_pools(qs_env%mpools,s=qs_env%matrix_s,&
         nspins=dft_control%nspins,particle_set=particle_set,&
         para_env=para_env, error=error)

    ! *** deallocate matrices that will have the wrong structure ***
    CALL deallocate_matrix_set(matrix_ks,error=error)
    CALL set_qs_env(qs_env,matrix_ks=matrix_ks,error=error)

    IF(dft_control%qs_control%gapw) THEN
      qs_env%qs_charges%total_rho_core_rspace=qs_env%local_rho_set%rhoz_tot
      IF (ASSOCIATED(qs_env%rho_core)) THEN
          CALL coeff_deallocate(qs_env%rho_core,error=error)
          DEALLOCATE(qs_env%rho_core)
      ENDIF
    ELSE IF(dft_control%qs_control%semi_empirical) THEN
      !??
    ELSE IF(dft_control%qs_control%dftb) THEN
      !??
    ELSE
      ! *** updates rho core ***
      CPPrecondition(ASSOCIATED(qs_env%rho_core),cp_failure_level,routineP,error,failure)
      CPPrecondition(ASSOCIATED(qs_env%pw_env),cp_failure_level,routineP,error,failure)
      CPPrecondition(qs_env%pw_env%ref_count>0,cp_failure_level,routineP,error,failure)
      IF (.NOT. failure) THEN
         CALL calculate_rho_core(qs_env%rho_core, &
              qs_env%qs_charges%total_rho_core_rspace,qs_env,error=error)
      END IF
    END IF
    
    ! *** environment for ddapc ***
    IF (ASSOCIATED(qs_env%cp_ddapc_env)) THEN
       CALL cp_ddapc_release(qs_env%cp_ddapc_env,error)
    END IF
    CALL cp_ddapc_init(qs_env, error)

    ! *** tell ks_env ***
    IF (ASSOCIATED(qs_env%ks_env)) THEN
       CALL qs_ks_did_change(qs_env%ks_env,s_mstruct_changed=.TRUE.,&
            error=error)
    END IF

    !   *** Updates rho structure ***
    CALL qs_env_rebuild_rho(qs_env=qs_env,error=error)

    ! *** tell scf_env ***
    IF (ASSOCIATED(qs_env%scf_env)) THEN
       CALL scf_env_did_change(qs_env%scf_env,&
            error=error)
    END IF

    CALL timestop(0.0_dp,handle)

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
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_env_rebuild_pw_env', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure
    TYPE(coeff_type), POINTER                :: rho_core
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_env_type), POINTER               :: new_pw_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool

    CALL timeset(routineN,handle)
    failure=.FALSE.
    ! rebuild pw_env
    NULLIFY(dft_control)
    CALL get_qs_env(qs_env, pw_env=new_pw_env, dft_control=dft_control,error=error)

    IF (.NOT.ASSOCIATED(new_pw_env)) THEN
       CALL pw_env_create(new_pw_env,qs_env=qs_env,error=error)
       CALL set_qs_env(qs_env,pw_env=new_pw_env,error=error)
       CALL pw_env_release(new_pw_env,error=error)
    ELSE
       CALL pw_env_rebuild(new_pw_env,qs_env=qs_env,error=error)
    END IF

    ! reallocate rho_core
    CALL get_qs_env(qs_env, pw_env=new_pw_env, rho_core=rho_core,error=error)
    CPPrecondition(ASSOCIATED(new_pw_env),cp_failure_level,routineP,error,failure)
    IF(dft_control%qs_control%gapw) THEN
      IF (ASSOCIATED(rho_core)) THEN
         CALL coeff_deallocate(rho_core,error=error)
         DEALLOCATE(rho_core,stat=stat)
         CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      END IF
    ELSE IF (dft_control%qs_control%semi_empirical) THEN
    ELSE IF (dft_control%qs_control%dftb) THEN
    ELSE
      IF (ASSOCIATED(rho_core)) THEN
         CALL coeff_deallocate(rho_core,error=error)
         DEALLOCATE(rho_core,stat=stat)
         CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      ENDIF
      ALLOCATE(rho_core,stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      IF (.NOT.failure) THEN
         CALL pw_env_get(new_pw_env, auxbas_pw_pool=auxbas_pw_pool,error=error)
         CALL pw_pool_init_coeff(pool=auxbas_pw_pool, coeff=rho_core, &
              use_data=COMPLEXDATA1D, error=error)
         rho_core%pw%in_space=RECIPROCALSPACE
         CALL set_qs_env(qs_env, rho_core=rho_core,error=error)
      END IF
    END IF

    CALL timestop(handle)

END SUBROUTINE qs_env_rebuild_pw_env
!***************************************************************************

!!****f* qs_environment_methods/qs_env_rebuild_rho [1.0] *
!!
!!   NAME
!!     qs_env_rebuild_rho
!!
!!   SYNOPSIS
!!     Subroutine qs_env_rebuild_rho(qs_env, &
!!         rebuild_ao, rebuild_grids, error)
!!       Type(qs_environment_type), Intent (INOUT), Target:: qs_env
!!       Logical, Intent (IN), Optional:: &
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
SUBROUTINE qs_env_rebuild_rho(qs_env, rebuild_ao, rebuild_grids, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(in), OPTIONAL            :: rebuild_ao, rebuild_grids
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_env_rebuild_rho', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: gapw_xc
    TYPE(dft_control_type), POINTER          :: dft_control

    CALL timeset(routineN,handle)

    CALL get_qs_env(qs_env=qs_env, dft_control=dft_control,error=error)
    gapw_xc=dft_control%qs_control%gapw_xc
    CALL qs_rho_rebuild(qs_env%rho,qs_env=qs_env,&
         rebuild_ao=rebuild_ao, rebuild_grids=rebuild_grids, error=error)
    IF(gapw_xc) THEN
      CALL qs_rho_rebuild(qs_env%rho_xc,qs_env=qs_env,&
         rebuild_ao=rebuild_ao, rebuild_grids=rebuild_grids, gapw_xc=gapw_xc, error=error)
      qs_env%rho_xc%rho_ao => qs_env%rho%rho_ao
    END IF

    CALL timestop(handle)

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
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_env_update_rho_core', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    TYPE(coeff_type), POINTER                :: rho_core
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(qs_charges_type), POINTER           :: qs_charges

    CALL timeset(routineN,handle)
    failure=.FALSE.
    CALL get_qs_env(qs_env, dft_control=dft_control,&
         pw_env=pw_env, rho_core=rho_core,&
         qs_charges=qs_charges,error=error)
    CPPrecondition(ASSOCIATED(pw_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qs_charges),cp_failure_level,routineP,error,failure)
    IF(.NOT. dft_control%qs_control%gapw) THEN
      CPPrecondition(ASSOCIATED(rho_core),cp_failure_level,routineP,error,failure)
      IF (.NOT. failure) THEN
         CALL calculate_rho_core(rho_core, &
            qs_charges%total_rho_core_rspace,qs_env,error)
      END IF
    ELSE
      qs_charges%total_rho_core_rspace=qs_env%local_rho_set%rhoz_tot
    END IF
    CALL timestop(handle)

END SUBROUTINE qs_env_update_rho_core
!***************************************************************************

END MODULE qs_environment_methods
