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
                                             get_atomic_kind
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE coefficient_types,               ONLY: coeff_deallocate,&
                                             coeff_type
  USE cp_array_i_utils,                ONLY: cp_1d_i_p_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE distribution_2d_types,           ONLY: distribution_2d_create,&
                                             distribution_2d_release,&
                                             distribution_2d_type,&
                                             distribution_2d_write
  USE distribution_methods,            ONLY: distribute_molecules_2d
  USE global_types,                    ONLY: LOW
  USE kinds,                           ONLY: dp,&
                                             int_size
  USE machine,                         ONLY: m_flush
  USE message_passing,                 ONLY: mp_sum
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
  USE qs_rho_methods,                  ONLY: qs_rho_rebuild
  USE qs_scf,                          ONLY: scf_env_did_change
  USE scf_control_types,               ONLY: scf_control_type
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
SUBROUTINE qs_env_setup(qs_env,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_env_setup', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
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

  failure=.FALSE.
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
                                  error=error)

     CALL set_qs_env(qs_env,distribution_2d=distribution_2d,error=error)
     CALL distribution_2d_release(distribution_2d, error=error)
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
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_env_update_s_mstruct', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

  CALL timeset(routineN,"I","",handle)
  failure=.FALSE.
  NULLIFY(para_env, dft_control,particle_set)

  CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)

  CALL get_qs_env(qs_env, para_env=para_env, dft_control=dft_control,&
       particle_set=particle_set, error=error)

  ! *** rebuilds the S_sm_pools ***
  CALL mpools_rebuild_s_sm_pools(qs_env%mpools,s=qs_env%matrix_s,&
       nspins=dft_control%nspins,particle_set=particle_set,&
       para_env=para_env, error=error)

  IF(.NOT. dft_control%qs_control%gapw) THEN
  ! *** updates rho core ***
    CPPrecondition(ASSOCIATED(qs_env%rho_core),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qs_env%pw_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(qs_env%pw_env%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL calculate_rho_core(qs_env%rho_core, &
            qs_env%qs_charges%total_rho_core_rspace,qs_env,error=error)
    END IF
  ELSE
    qs_env%qs_charges%total_rho_core_rspace=qs_env%local_rho_set%rhoz_tot
  END IF

  ! *** tell ks_env ***
  IF (ASSOCIATED(qs_env%ks_env)) THEN
     CALL qs_ks_did_change(qs_env%ks_env,s_mstruct_changed=.TRUE.,&
          error=error)
  END IF

  !   *** Updates rho structure ***
  CALL qs_env_rebuild_rho(qs_env=qs_env)

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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_env_rebuild_pw_env', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure
    TYPE(coeff_type), POINTER                :: rho_core
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_env_type), POINTER               :: new_pw_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool

  failure=.FALSE.

  ! rebuild pw_env
  NULLIFY(dft_control)
  CALL get_qs_env(qs_env, pw_env=new_pw_env, dft_control=dft_control)
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
  IF(dft_control%qs_control%gapw) THEN
    IF (ASSOCIATED(rho_core)) THEN
       CALL coeff_deallocate(rho_core)
       DEALLOCATE(rho_core,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF
  ELSE
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
  END IF

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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_env_rebuild_rho', &
      routineP = moduleN//':'//routineN

  CALL qs_rho_rebuild(qs_env%rho,qs_env=qs_env,&
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
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_env_update_rho_core', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(coeff_type), POINTER                :: rho_core
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(qs_charges_type), POINTER           :: qs_charges

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
END SUBROUTINE qs_env_update_rho_core
!***************************************************************************

!!****f* qs_environment_methods/qs_distribution_2d_create [1.0] *
!!
!!   NAME
!!     qs_distribution_2d_create
!!
!!   FUNCTION
!!     distributes the particle pairs creating a 2d distribution
!!     optimally suited for quickstep
!!
!!   NOTES
!!
!!   ARGUMENTS
!!     - distribution_2d: the distribution that will be created by this method
!!     - atomic_kind_set: the atomic kind set to distribute
!!     - blacs_env: the parallel environement at the basis of the distribution
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!     07.2003 local_rows & cols blocksize optimizations [matthias]
!!     09.2003 cleanup of distribution_2d [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE qs_distribution_2d_create(distribution_2d,atomic_kind_set,&
       blacs_env, error)

    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "qs_distribution_2d_create", &
      routineP = "SUBROUTINE "//routineN//" (MODULE "//moduleN//")"

    INTEGER :: group, i, iatom, iatom_local, ikind, ipcol, iprow, j, mypcol, &
      myprow, natom, natom_total, nkind, npcol, nprow, nsgf, output_unit, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: natom_local_col, &
                                                natom_local_row, ncol, &
                                                ncol_old, nparticle, nrow, &
                                                nrow_old
    INTEGER, DIMENSION(:), POINTER           :: atom_list, col_distribution, &
                                                row_distribution
    LOGICAL                                  :: failure, ionode
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_1d_i_p_type), DIMENSION(:), &
      POINTER                                :: local_cols, local_rows
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set

!   ---------------------------------------------------------------------------

    NULLIFY(atomic_kind, orb_basis_set, local_rows, local_cols, atom_list, row_distribution,col_distribution)
    logger => cp_error_get_logger(error)
    output_unit = cp_logger_get_default_unit_nr(logger)
    ionode= (logger%para_env%mepos == logger%para_env%source)

    group = blacs_env%para_env%group

    myprow = blacs_env%mepos(1) + 1
    mypcol = blacs_env%mepos(2) + 1
    nprow = blacs_env%num_pe(1)
    npcol = blacs_env%num_pe(2)
    nkind = SIZE(atomic_kind_set)

    ALLOCATE (nrow(nprow),nrow_old(nprow),ncol(npcol),ncol_old(npcol),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    nrow=0
    ncol=0
    ALLOCATE (local_rows(nkind), local_cols(nkind),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (natom_local_row(nkind),natom_local_col(nkind),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    natom_total = 0
    DO ikind=1,nkind
      atomic_kind => atomic_kind_set(ikind)
      NULLIFY (atom_list)
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           atom_list=atom_list)
      natom_total = natom_total + SIZE(atom_list)
    ENDDO
    ALLOCATE(row_distribution(natom_total), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(col_distribution(natom_total), stat=stat) 
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

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

      ALLOCATE (local_rows(ikind)%array(natom_local_row(ikind)),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

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
          local_rows(ikind)%array(iatom_local) = atom_list(iatom)
        END IF
      END DO

!     *** Atomic block columns ***

      ALLOCATE (local_cols(ikind)%array(natom_local_col(ikind)),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

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
          local_cols(ikind)%array(iatom_local) = atom_list(iatom)
        END IF
      END DO

    END DO

!   *** generate the distribution_2d structure
    CALL distribution_2d_create(distribution_2d, &
         row_distribution_ptr=row_distribution, &
         col_distribution_ptr=col_distribution,&
         local_cols_ptr=local_cols,local_rows_ptr=local_rows,&
         blacs_env=blacs_env, error=error)
    NULLIFY(local_rows, local_cols, row_distribution, col_distribution)

!   *** Print distribution, if requested ***
    
    IF (logger%print_keys%distribution) THEN
!     *** Print row distribution ***

      ALLOCATE (nparticle(nprow),STAT=stat)
      IF (stat /= 0) CALL stop_memory(routineP,"nparticle",nprow*int_size)
      nparticle(:) = 0

      IF (mypcol == 1) nparticle(myprow) = SUM(distribution_2d%n_local_rows)

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

      DEALLOCATE (nparticle,STAT=stat)
      IF (stat /= 0) CALL stop_memory(routineP,"nparticle")

!     *** Print column distribution ***

      ALLOCATE (nparticle(npcol),STAT=stat)
      IF (stat /= 0) CALL stop_memory(routineP,"nparticle",npcol*int_size)
      nparticle(:) = 0

      IF (myprow == 1) nparticle(mypcol) = SUM(distribution_2d%n_local_cols)

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

      DEALLOCATE (nparticle,STAT=stat)
      IF (stat /= 0) CALL stop_memory(routineP,"nparticle")

      IF (logger%print_keys%level>LOW) THEN
         CALL distribution_2d_write(distribution_2d, &
              unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
              local=.TRUE.,long_description=.TRUE.,&
              error=error)
      END IF

    END IF

!   *** Release work storage ***
    DEALLOCATE (nrow,STAT=stat)
    IF (stat /= 0) CALL stop_memory(routineP,"nrow")
    DEALLOCATE (nrow_old,STAT=stat)
    IF (stat /= 0) CALL stop_memory(routineP,"nrow_old")
    DEALLOCATE (ncol,STAT=stat)
    IF (stat /= 0) CALL stop_memory(routineP,"ncol")
    DEALLOCATE (ncol_old,STAT=stat)
    IF (stat /= 0) CALL stop_memory(routineP,"ncol_old")
    DEALLOCATE (natom_local_row,STAT=stat)
    IF (stat /= 0) CALL stop_memory(routineP,"natom_local_row")
    DEALLOCATE (natom_local_col,STAT=stat)
    IF (stat /= 0) CALL stop_memory(routineP,"natom_local_col")

  END SUBROUTINE qs_distribution_2d_create

! *****************************************************************************

END MODULE qs_environment_methods
