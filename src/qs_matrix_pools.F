!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/qs_matrix_pools [1.0] *
!!
!!   NAME
!!     qs_matrix_pools
!!
!!   FUNCTION
!!     wrapper for the pools of matrixes
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     05.2003 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE qs_matrix_pools
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                             cp_fm_pool_type,&
                                             fm_pool_create,&
                                             fm_pool_get_el_struct,&
                                             fm_pool_release,&
                                             fm_pool_retain,&
                                             fm_pools_copy,&
                                             fm_pools_dealloc
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_get,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_matrix_utils,                 ONLY: sm_struct_add_sm_local_blocks
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE cp_sm_pool_types,                ONLY: cp_sm_pool_p_type,&
                                             sm_pool_create,&
                                             sm_pool_retain,&
                                             sm_pools_dealloc
  USE cp_sm_struct,                    ONLY: cp_sm_struct_type,&
                                             sm_struct_create,&
                                             sm_struct_release
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: get_particle_set,&
                                             particle_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE sparse_matrix_types,             ONLY: get_matrix_info,&
                                             real_matrix_p_type
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qs_matrix_pools'
  INTEGER, SAVE, PRIVATE :: last_mpools_id=0

  PUBLIC :: qs_matrix_pools_type
  PUBLIC :: mpools_retain, mpools_release, mpools_set, mpools_get,&
       mpools_create, mpools_rebuild_fm_pools, &
       mpools_rebuild_s_sm_pools

!!***
!****************************************************************************

!!****s* m/qs_matrix_pools_type [1.0] *
!!
!!   NAME
!!     qs_matrix_pools_type
!!
!!   FUNCTION
!!     container for the pools of matrixes used by qs
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - id_nr: identification number
!!     - ref_count: reference count (see doc/ReferenceCounting.html)
!!     - S_sm_pools: pools of sparse matrixes with the structure of S
!!     - ao_mo_fm_pools: pools with (ao x mo) full matrixes (same order as
!!       c).
!!     - ao_ao_fm_pools: pools with (ao x ao) full matrixes (same order as
!!       c).
!!     - mo_mo_fm_pools: pools with (mo x mo) full matrixes (same
!!       order as c).
!!
!!     - maxao_maxao_fm_pools: pool of matrixes big enough to accomodate any
!!       aoxao matrix (useful for temp matrixes)
!!     - maxao_maxmo_fm_pools: pool of matrixes big enough to accomodate any
!!       aoxmo matrix (useful for temp matrixes)
!!     - maxmo_maxmo_fm_pools: pool of matrixes big enough to accomodate any
!!       moxmo matrix (useful for temp matrixes)
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     04.2003 created [fawzi]
!!
!!   SOURCE
!***************************************************************************
  TYPE qs_matrix_pools_type
     INTEGER :: id_nr, ref_count
     TYPE(cp_sm_pool_p_type), DIMENSION(:), POINTER        :: S_sm_pools
     TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER        :: ao_mo_fm_pools,&
          ao_ao_fm_pools,mo_mo_fm_pools
  END TYPE qs_matrix_pools_type
!!***
!****************************************************************************

CONTAINS

!!****f* qs_matrix_pools/mpools_retain [1.0] *
!!
!!   NAME
!!     mpools_retain
!!
!!   SYNOPSIS
!!     Subroutine mpools_retain(mpools, error)
!!       Type(qs_matrix_pools_type), Pointer:: mpools
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine mpools_retain
!!
!!   FUNCTION
!!     retains the given qs_matrix_pools_type
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - mpools: the matrix pools type to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     04.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE mpools_retain(mpools, error)
    TYPE(qs_matrix_pools_type), POINTER      :: mpools
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mpools_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(mpools),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPreconditionNoFail(mpools%ref_count>0,cp_failure_level,routineP,error)
     mpools%ref_count=mpools%ref_count+1
  END IF
END SUBROUTINE mpools_retain
!***************************************************************************

!!****f* qs_matrix_pools/mpools_release [1.0] *
!!
!!   NAME
!!     mpools_release
!!
!!   SYNOPSIS
!!     Subroutine mpools_release(mpools, error)
!!       Type(qs_matrix_pools_type), Pointer:: mpools
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine mpools_release
!!
!!   FUNCTION
!!     releases the given mpools
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - mpools: the matrix pools type to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     04.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE mpools_release(mpools,error)
    TYPE(qs_matrix_pools_type), POINTER      :: mpools
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mpools_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  IF (ASSOCIATED(mpools)) THEN
     CPPreconditionNoFail(mpools%ref_count>0,cp_failure_level,routineP,error)
     mpools%ref_count=mpools%ref_count-1
     IF (mpools%ref_count==0) THEN
        CALL sm_pools_dealloc(mpools%S_sm_pools,error=error)
        CALL fm_pools_dealloc(mpools%ao_mo_fm_pools, error=error)
        CALL fm_pools_dealloc(mpools%ao_ao_fm_pools, error=error)
        CALL fm_pools_dealloc(mpools%mo_mo_fm_pools, error=error)
        DEALLOCATE(mpools, stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(mpools)
END SUBROUTINE mpools_release
!***************************************************************************

!!****f* qs_matrix_pools/mpools_set [1.0] *
!!
!!   NAME
!!     mpools_set
!!
!!   SYNOPSIS
!!     Subroutine mpools_set(mpools, ao_mo_fm_pools, ao_ao_fm_pools,&
!!         mo_mo_fm_pools, s_sm_pools, error)
!!       Type(qs_matrix_pools_type), Pointer:: mpools
!!       Type(cp_sm_pool_p_type), Dimension(:), Optional, Pointer::&
!!         s_sm_pools
!!       Type(cp_fm_pool_p_type), Dimension(:), Optional, Pointer::&
!!         ao_mo_fm_pools, ao_ao_fm_pools, mo_mo_fm_pools
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine mpools_set
!!
!!   FUNCTION
!!     sets various attributes of the mpools
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
!!     04.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE mpools_set(mpools,ao_mo_fm_pools,ao_ao_fm_pools,&
     mo_mo_fm_pools,&
     S_sm_pools,error)
    TYPE(qs_matrix_pools_type), POINTER      :: mpools
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: ao_mo_fm_pools, &
                                                ao_ao_fm_pools, mo_mo_fm_pools
    TYPE(cp_sm_pool_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: S_sm_pools
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mpools_set', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: new_fm_pools

  failure=.FALSE.
  NULLIFY(new_fm_pools )
  
  CPPrecondition(ASSOCIATED(mpools),cp_failure_level,routineP,error,failure)
  CPPrecondition(mpools%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
    IF (PRESENT(ao_mo_fm_pools)) THEN
      IF (ASSOCIATED(ao_mo_fm_pools)) THEN
        CALL fm_pools_copy(ao_mo_fm_pools,new_fm_pools,error=error)
      END IF
      CALL fm_pools_dealloc(mpools%ao_mo_fm_pools,error=error)
      mpools%ao_mo_fm_pools => new_fm_pools
    END IF
    IF (PRESENT(ao_ao_fm_pools)) THEN
      IF (ASSOCIATED(ao_ao_fm_pools)) THEN
        CALL fm_pools_copy(ao_ao_fm_pools,new_fm_pools,error=error)
      END IF
      CALL fm_pools_dealloc(mpools%ao_ao_fm_pools,error=error)
      mpools%ao_ao_fm_pools => new_fm_pools
    END IF
    IF (PRESENT(mo_mo_fm_pools)) THEN
      IF (ASSOCIATED(mo_mo_fm_pools)) THEN
        CALL fm_pools_copy(mo_mo_fm_pools,new_fm_pools,error=error)
      END IF
      CALL fm_pools_dealloc(mpools%mo_mo_fm_pools,error=error)
      mpools%mo_mo_fm_pools => new_fm_pools
    END IF
 END IF
END SUBROUTINE mpools_set
!***************************************************************************

!!****f* qs_matrix_pools/mpools_get [1.0] *
!!
!!   NAME
!!     mpools_get
!!
!!   SYNOPSIS
!!     Subroutine mpools_get(mpools, ao_mo_fm_pools, ao_ao_fm_pools,&
!!         mo_mo_fm_pools, maxao_maxmo_fm_pool, maxao_maxao_fm_pool,&
!!         maxmo_maxmo_fm_pool, s_sm_pools, error)
!!       Type(qs_matrix_pools_type), Pointer:: mpools
!!       Type(cp_sm_pool_p_type), Dimension(:), Optional, Pointer::&
!!         s_sm_pools
!!       Type(cp_fm_pool_p_type), Dimension(:), Optional, Pointer::&
!!         ao_mo_fm_pools, ao_ao_fm_pools, mo_mo_fm_pools
!!       Type(cp_fm_pool_type), Optional, Pointer:: maxao_maxmo_fm_pool,&
!!         maxao_maxao_fm_pool, maxmo_maxmo_fm_pool
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine mpools_get
!!
!!   FUNCTION
!!     returns various attributes of the mpools (notably the pools
!!     contained in it)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - mpools: the matrix pools object you want info about
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!     see qs_matrix_pools_type attributes
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     04.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE mpools_get(mpools,ao_mo_fm_pools,ao_ao_fm_pools,&
     mo_mo_fm_pools, maxao_maxmo_fm_pool,&
     maxao_maxao_fm_pool, maxmo_maxmo_fm_pool,&
     S_sm_pools,error)
    TYPE(qs_matrix_pools_type), POINTER      :: mpools
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: ao_mo_fm_pools, &
                                                ao_ao_fm_pools, mo_mo_fm_pools
    TYPE(cp_fm_pool_type), OPTIONAL, POINTER :: maxao_maxmo_fm_pool, &
                                                maxao_maxao_fm_pool, &
                                                maxmo_maxmo_fm_pool
    TYPE(cp_sm_pool_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: S_sm_pools
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mpools_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  IF (.NOT. failure) THEN
         IF (PRESENT(S_sm_pools)) S_sm_pools => mpools%S_sm_pools
    IF (PRESENT(ao_mo_fm_pools)) ao_mo_fm_pools => mpools%ao_mo_fm_pools
    IF (PRESENT(maxao_maxmo_fm_pool)) THEN
       IF (ASSOCIATED(mpools%ao_mo_fm_pools)) THEN
          maxao_maxmo_fm_pool => mpools%ao_mo_fm_pools(1)%pool
       ELSE
          NULLIFY(maxao_maxmo_fm_pool) ! raise an error?
       END IF
    END IF
    IF (PRESENT(ao_ao_fm_pools)) ao_ao_fm_pools => mpools%ao_ao_fm_pools
    IF (PRESENT(maxao_maxao_fm_pool)) THEN
       IF (ASSOCIATED(mpools%ao_ao_fm_pools)) THEN
          maxao_maxao_fm_pool => mpools%ao_ao_fm_pools(1)%pool
       ELSE
          NULLIFY(maxao_maxao_fm_pool) ! raise an error?
       END IF
    END IF
    IF (PRESENT(mo_mo_fm_pools)) mo_mo_fm_pools => mpools%mo_mo_fm_pools
    IF (PRESENT(maxmo_maxmo_fm_pool)) THEN
       IF (ASSOCIATED(mpools%mo_mo_fm_pools)) THEN
          maxmo_maxmo_fm_pool => mpools%mo_mo_fm_pools(1)%pool
       ELSE
          NULLIFY(maxmo_maxmo_fm_pool) ! raise an error?
       END IF
    END IF
 END IF
END SUBROUTINE mpools_get
!***************************************************************************

!!****f* qs_matrix_pools/mpools_create [1.0] *
!!
!!   NAME
!!     mpools_create
!!
!!   SYNOPSIS
!!     Subroutine mpools_create(mpools, error)
!!       Type(qs_matrix_pools_type), Pointer:: mpools
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine mpools_create
!!
!!   FUNCTION
!!     creates a mpools
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - mpools: the mpools to create
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     04.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE mpools_create(mpools,error)
    TYPE(qs_matrix_pools_type), POINTER      :: mpools
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mpools_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  ALLOCATE(mpools, stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     NULLIFY(mpools%ao_ao_fm_pools, mpools%ao_mo_fm_pools, &
          mpools%mo_mo_fm_pools, mpools%S_sm_pools)
     mpools%ref_count=1
     last_mpools_id=last_mpools_id+1
     mpools%id_nr=last_mpools_id
  END IF
END SUBROUTINE mpools_create
!***************************************************************************

!!****f* qs_matrix_pools/mpools_rebuild_fm_pools [1.0] *
!!
!!   NAME
!!     mpools_rebuild_fm_pools
!!
!!   SYNOPSIS
!!     Subroutine mpools_rebuild_fm_pools(mpools, mos, blacs_env, para_env,&
!!         nrow_block, ncol_block, error)
!!       Type(qs_matrix_pools_type), Pointer:: mpools
!!       Type(mo_set_p_type), Dimension(:), Pointer:: mos
!!       Type(cp_blacs_env_type), Pointer:: blacs_env
!!       Type(cp_para_env_type), Pointer:: para_env
!!       Integer, Intent (IN):: nrow_block, ncol_block
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine mpools_rebuild_fm_pools
!!
!!   FUNCTION
!!     rebuilds the pools of the (ao x mo, ao x ao , mo x mo) full matrixes
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - mpools: the environement where the pools should be rebuilt
!!     - mos: the molecular orbitals (qs_env%c), must contain up to 
!!       date nmo and nao
!!     - blacs_env: the blacs environment of the full matrixes
!!     - para_env: the parallel environment of the matrixes
!!     - nrow_block: number of rows of a block for for the block cyclic
!!       distribution of the elements in parallel
!!     - ncol_block: number of cols of a block for for the block cyclic
!!       distribution of the elements in parallel
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
SUBROUTINE mpools_rebuild_fm_pools(mpools, mos, blacs_env, para_env,&
     nrow_block,ncol_block,error)
    TYPE(qs_matrix_pools_type), POINTER      :: mpools
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(in)                      :: nrow_block, ncol_block
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mpools_rebuild_fm_pools', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, max_nmo, &
                                                min_nmo, nao, ncg, nmo, nrg, &
                                                nspins, stat
    LOGICAL                                  :: failure, should_rebuild
    TYPE(cp_fm_pool_type), POINTER           :: p_att
    TYPE(cp_fm_struct_type), POINTER         :: fmstruct

  CALL timeset(routineN,"I","",handle)

  failure=.FALSE.
  NULLIFY(fmstruct, p_att)

  CPPrecondition(ASSOCIATED(mos),cp_failure_level,routineP,error,failure)
  IF (.NOT.ASSOCIATED(mpools)) THEN
     CALL mpools_create(mpools,error=error)
  END IF
  IF (.NOT.failure) THEN
     nspins=SIZE(mos)

     IF (ASSOCIATED(mpools%ao_mo_fm_pools)) THEN
        IF (nspins/=SIZE(mpools%ao_mo_fm_pools)) THEN
           CALL fm_pools_dealloc(mpools%ao_mo_fm_pools,error=error)
        END IF
     END IF
     IF (.NOT.ASSOCIATED(mpools%ao_mo_fm_pools)) THEN
        ALLOCATE(mpools%ao_mo_fm_pools(nspins), stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        DO ispin=1,nspins
           NULLIFY(mpools%ao_mo_fm_pools(ispin)%pool)
        END DO
     END IF

     IF (ASSOCIATED(mpools%ao_ao_fm_pools)) THEN
        IF (nspins/=SIZE(mpools%ao_ao_fm_pools)) THEN
           CALL fm_pools_dealloc(mpools%ao_ao_fm_pools,error=error)
        END IF
     END IF
     IF (.NOT.ASSOCIATED(mpools%ao_ao_fm_pools)) THEN
        ALLOCATE(mpools%ao_ao_fm_pools(nspins), stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        DO ispin=1,nspins
           NULLIFY(mpools%ao_ao_fm_pools(ispin)%pool)
        END DO
     END IF

     IF (ASSOCIATED(mpools%mo_mo_fm_pools)) THEN
        IF (nspins/=SIZE(mpools%mo_mo_fm_pools)) THEN
           CALL fm_pools_dealloc(mpools%mo_mo_fm_pools,error=error)
        END IF
     END IF
     IF (.NOT.ASSOCIATED(mpools%mo_mo_fm_pools)) THEN
        ALLOCATE(mpools%mo_mo_fm_pools(nspins), stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        DO ispin=1,nspins
           NULLIFY(mpools%mo_mo_fm_pools(ispin)%pool)
        END DO
     END IF
  END IF

  IF (.NOT. failure) THEN
     CALL get_mo_set(mos(1)%mo_set,nao=nao,nmo=min_nmo)
     max_nmo=min_nmo
     DO ispin=2,SIZE(mos)
        CALL get_mo_set(mos(ispin)%mo_set,nmo=nmo)
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
     ! aoao pools
     should_rebuild=.FALSE.
     DO ispin=1,nspins
        p_att => mpools%ao_ao_fm_pools(ispin)%pool
        should_rebuild=should_rebuild.or.&
             .NOT.ASSOCIATED(p_att)
        IF (.NOT.should_rebuild) THEN
           fmstruct => fm_pool_get_el_struct(mpools%ao_ao_fm_pools(ispin)%pool,&
                error=error)
           CALL cp_fm_struct_get(fmstruct, nrow_global=nrg,&
                ncol_global=ncg,error=error)
           CALL get_mo_set(mos(1)%mo_set,nao=nao,nmo=nmo)
           should_rebuild = nao/=nrg.OR.nao/=ncg
        END IF
     END DO
     IF (should_rebuild) THEN
        DO ispin=1,nspins
           CALL fm_pool_release(mpools%ao_ao_fm_pools(ispin)%pool,error=error)
        END DO

        CALL cp_fm_struct_create(fmstruct, nrow_global=nao,&
             ncol_global=nao, para_env=para_env,&
             context=blacs_env,&
             nrow_block=nrow_block,&
             ncol_block=ncol_block)
        CALL fm_pool_create(mpools%ao_ao_fm_pools(1)%pool,fmstruct,error=error)
        CALL cp_fm_struct_release(fmstruct,error=error)
        DO ispin=2,SIZE(mos)
           mpools%ao_ao_fm_pools(ispin)%pool => mpools%ao_ao_fm_pools(1)%pool
           CALL fm_pool_retain(mpools%ao_ao_fm_pools(1)%pool,error=error)
        END DO
     END IF

     ! aomo pools
     should_rebuild=.FALSE.
     DO ispin=1,nspins
        p_att => mpools%ao_mo_fm_pools(ispin)%pool
        should_rebuild=should_rebuild.or.&
             .NOT.ASSOCIATED(p_att)
        IF (.NOT.should_rebuild) THEN
           fmstruct => fm_pool_get_el_struct(mpools%ao_mo_fm_pools(ispin)&
                %pool,error=error)
           CALL cp_fm_struct_get(fmstruct, nrow_global=nrg,&
                ncol_global=ncg,error=error)
           CALL get_mo_set(mos(1)%mo_set,nao=nao,nmo=nmo)
           should_rebuild = nao/=nrg.OR.nmo/=ncg
        END IF
     END DO
     IF (should_rebuild) THEN
        DO ispin=1,nspins
           CALL fm_pool_release(mpools%ao_mo_fm_pools(ispin)%pool,error=error)
        END DO

        IF (max_nmo==min_nmo) THEN
           CALL cp_fm_struct_create(fmstruct, nrow_global=nao,&
                ncol_global=max_nmo, para_env=para_env,&
                context=blacs_env,&
                nrow_block=nrow_block,&
                ncol_block=ncol_block)
           CALL fm_pool_create(mpools%ao_mo_fm_pools(1)%pool,fmstruct,error=error)
           CALL cp_fm_struct_release(fmstruct,error=error)
           DO ispin=2,SIZE(mos)
              mpools%ao_mo_fm_pools(ispin)%pool => mpools%ao_mo_fm_pools(1)%pool
              CALL fm_pool_retain(mpools%ao_mo_fm_pools(1)%pool,error=error)
           END DO
        ELSE
           DO ispin=1,SIZE(mos)
              CALL get_mo_set(mos(ispin)%mo_set,nmo=nmo,nao=nao) 
              CALL cp_fm_struct_create(fmstruct, nrow_global=nao,&
                   ncol_global=nmo, para_env=para_env,&
                   context=blacs_env,&
                   nrow_block=nrow_block,&
                   ncol_block=ncol_block)
              CALL fm_pool_create(mpools%ao_mo_fm_pools(ispin)%pool,&
                   fmstruct,error=error)
              CALL cp_fm_struct_release(fmstruct,error=error)
           END DO
        END IF
     END IF

     ! momo pools
     should_rebuild=.FALSE.
     DO ispin=1,nspins
        p_att => mpools%mo_mo_fm_pools(ispin)%pool
        should_rebuild=should_rebuild.OR.&
             .NOT.ASSOCIATED(p_att)
        IF (.NOT.should_rebuild) THEN
           fmstruct => fm_pool_get_el_struct(p_att,error=error)
           CALL cp_fm_struct_get(fmstruct, nrow_global=nrg,&
                ncol_global=ncg,error=error)
           CALL get_mo_set(mos(1)%mo_set,nao=nao,nmo=nmo)
           should_rebuild = nmo/=nrg.OR.nmo/=ncg
        END IF
     END DO
     IF (should_rebuild) THEN
        DO ispin=1,nspins
           CALL fm_pool_release(mpools%mo_mo_fm_pools(ispin)%pool,error=error)
        END DO

        IF (max_nmo==min_nmo) THEN
           CALL cp_fm_struct_create(fmstruct, nrow_global=max_nmo,&
                ncol_global=max_nmo, para_env=para_env,&
                context=blacs_env,&
                nrow_block=nrow_block,&
                ncol_block=ncol_block)
           CALL fm_pool_create(mpools%mo_mo_fm_pools(1)%pool,&
                fmstruct,error=error)
           CALL cp_fm_struct_release(fmstruct,error=error)
           DO ispin=2,SIZE(mos)
              mpools%mo_mo_fm_pools(ispin)%pool => mpools%mo_mo_fm_pools(1)%pool
              CALL fm_pool_retain(mpools%mo_mo_fm_pools(1)%pool,error=error)
           END DO
        ELSE
           DO ispin=1,SIZE(mos)
              NULLIFY(mpools%mo_mo_fm_pools(ispin)%pool)
              CALL get_mo_set(mos(ispin)%mo_set,nmo=nmo,nao=nao) 
              CALL cp_fm_struct_create(fmstruct, nrow_global=nmo,&
                   ncol_global=nmo, para_env=para_env,&
                   context=blacs_env,&
                   nrow_block=nrow_block,&
                   ncol_block=ncol_block)
              CALL fm_pool_create(mpools%mo_mo_fm_pools(ispin)%pool,&
                   fmstruct,error=error)
              CALL cp_fm_struct_release(fmstruct,error=error)
           END DO
        END IF
     END IF
  END IF

  CALL timestop(0.0_dp,handle)
END SUBROUTINE mpools_rebuild_fm_pools
!***************************************************************************


!!****f* qs_matrix_pools/mpools_rebuild_s_sm_pools [1.0] *
!!
!!   NAME
!!     mpools_rebuild_s_sm_pools
!!
!!   SYNOPSIS
!!     Subroutine mpools_rebuild_s_sm_pools(mpools, s, nspins, particle_set,&
!!         para_env, error)
!!       Type(qs_matrix_pools_type), Pointer:: mpools
!!       Type(real_matrix_p_type), Dimension(:), Pointer:: s
!!       Integer, Intent (IN):: nspins
!!       Type(particle_type), Pointer, Dimension(:):: particle_set
!!       Type(cp_para_env_type), Pointer:: para_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine mpools_rebuild_s_sm_pools
!!
!!   FUNCTION
!!     rebuilds mpools%S_sm_pools
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - mpools: the pools to rebuild
!!     - s: a matrix vector with the structure of the s matrix
!!     - nspins: the number of spins
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     04.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE mpools_rebuild_s_sm_pools(mpools, s, nspins, particle_set, &
     para_env, error)
    TYPE(qs_matrix_pools_type), POINTER      :: mpools
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: s
    INTEGER, INTENT(in)                      :: nspins
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mpools_rebuild_s_sm_pools', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, n_atoms, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: block_begins_at
    INTEGER, DIMENSION(:), POINTER           :: first_col, first_row, &
                                                last_col, last_row
    LOGICAL                                  :: failure
    TYPE(cp_sm_struct_type), POINTER         :: s_mstruct

  CALL timeset(routineN,"I","",handle)
  failure=.FALSE.

  NULLIFY(first_row, first_col, last_row, last_col, s_mstruct)
  CPPrecondition(ASSOCIATED(particle_set),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(para_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(s),cp_failure_level,routineP,error,failure)
  IF (.NOT.ASSOCIATED(mpools)) THEN
     CALL mpools_create(mpools,error=error)
  END IF
  IF (.NOT. failure) THEN
     IF (.NOT.ASSOCIATED(mpools)) THEN
        CALL mpools_create(mpools,error=error)
     END IF
     CALL sm_pools_dealloc(mpools%S_sm_pools,error=error) ! update instead?

     n_atoms = SIZE(particle_set)
     ALLOCATE (block_begins_at(n_atoms+1),STAT=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     CALL get_particle_set(particle_set=particle_set,first_sgf=block_begins_at)

     CALL sm_struct_create(s_mstruct,&
          block_begins_at_row=block_begins_at,&
          block_begins_at_col=block_begins_at,&
          para_env=para_env, symmetry="symmetric",&
          distribution_2d=s(1)%matrix%distribution_2d,&
          sparsity_id=s(1)%matrix%sparsity_id, &
          error=error)

     DEALLOCATE(block_begins_at,stat=stat)
     CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

     CALL get_matrix_info(s(1)%matrix,first_row=first_row,&
          last_row=last_row, first_col=first_col,&
          last_col=last_col)

     CALL cp_assert(ALL(s_mstruct%block_begins_at_row(1:n_atoms)==&
          first_row),cp_failure_level,cp_assertion_failed,routineP,&
          "s has unexpected size, "//&
          CPSourceFileRef,&
          error)
     CALL cp_assert(ALL(s_mstruct%block_begins_at_col(1:n_atoms)==&
          first_col),cp_failure_level,cp_assertion_failed,routineP,&
          "s has unexpected size, "//&
          CPSourceFileRef,&
          error)
     CALL cp_assert(s_mstruct%block_begins_at_col(n_atoms+1)==&
          last_col(n_atoms)+1,cp_failure_level,cp_assertion_failed,routineP,&
          "s has unexpected size, "//&
          CPSourceFileRef,&
          error)
     CALL cp_assert(s_mstruct%block_begins_at_row(n_atoms+1)==&
          last_row(n_atoms)+1,cp_failure_level,cp_assertion_failed,routineP,&
          "s has unexpected size, "//&
          CPSourceFileRef,&
          error)

     ! at the moment read from s, maybe in the future s will take its structure
     ! from this
     CALL sm_struct_add_sm_local_blocks(s_mstruct,&
          matrix=s(1)%matrix,&
          error=error)

     ALLOCATE(mpools%S_sm_pools(nspins),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

     NULLIFY(mpools%S_sm_pools(1)%pool)
     CALL sm_pool_create(mpools%S_sm_pools(1)%pool,&
          matrix_struct=s_mstruct,&
          error=error)
     CALL sm_struct_release(s_mstruct,error=error)

     DO ispin=2,nspins
        mpools%S_sm_pools(ispin)%pool => mpools%S_sm_pools(1)%pool
        CALL sm_pool_retain(mpools%S_sm_pools(ispin)%pool,error=error)
     END DO

  END IF
  CALL timestop(0.0_dp,handle)

END SUBROUTINE mpools_rebuild_s_sm_pools
!***************************************************************************

END MODULE qs_matrix_pools


