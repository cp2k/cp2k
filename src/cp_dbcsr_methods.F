!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Methods for the CP2K DBCSR type
!> \author  Urban Borstnik
!> \date    2010-08-11
!> \version 0.1
!>
!> <b>Modification history:</b>
!> - Created 2010-08-11 [UB]
! *****************************************************************************
MODULE cp_dbcsr_methods
  USE array_types,                     ONLY: array_i1d_obj
  USE cluster_methods,                 ONLY: cluster_allocate,&
                                             cluster_deallocate
  USE cluster_types,                   ONLY: cluster_type
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_mapping_type,&
                                             cp_dbcsr_type

  !$ USE OMP_LIB
#include "cp_common_uses.h"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dbcsr_methods'

  PRIVATE

!  PUBLIC :: cp_dbcsr_get_mapping, cp_dbcsr_get_dbcsr

  PUBLIC :: cp_dbcsr_col_block_sizes, cp_dbcsr_row_block_sizes
  PUBLIC :: cp_dbcsr_col_block_offsets, cp_dbcsr_row_block_offsets
  PUBLIC :: cp_dbcsr_mapping_create, cp_dbcsr_mapping_hold,&
            cp_dbcsr_mapping_release
  PUBLIC :: cp_dbcsr_mapping_defined, cp_dbcsr_has_mapping,&
            cp_dbcsr_mapping_activate


CONTAINS

!  PURE FUNCTION cp_dbcsr_get_dbcsr_matrix (matrix) RESULT (dbcsr_matrix)
!    type(cp_dbcsr_type), intent(IN) :: matrix
!    type(dbcsr_obj) :: dbcsr_matrix
!    dbcsr_matrix = matrix%matrix
!  END FUNCTION cp_dbcsr_get_dbcsr
!
!  PURE FUNCTION cp_dbcsr_get_mapping (matrix) RESULT (mapping)
!    type(cp_dbcsr_type), intent(IN) :: matrix
!    type(mapping_type) :: mapping
!    mapping = matrix%mapping
!  END FUNCTION cp_dbcsr_get_mapping
!
!  PURE FUNCTION cp_dbcsr_get_dbcsr_matrix (matrix) RESULT (dbcsr_matrix)
!    type(cp_dbcsr_type), intent(IN) :: matrix
!    type(dbcsr_obj) :: dbcsr_matrix
!    dbcsr_matrix = matrix%matrix
!  END FUNCTION cp_dbcsr_get_dbcsr
!
!  PURE FUNCTION cp_dbcsr_get_dbcsr_dist (matrix) RESULT (distribution)
!    type(cp_dbcsr_type), intent(IN) :: matrix
!    type(dbcsr_distribution_obj) :: distribution
!    distribution = matrix%distribution
!  END FUNCTION cp_dbcsr_get_dbcsr


  FUNCTION cp_dbcsr_row_block_sizes (matrix) RESULT (row_blk_sizes)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    TYPE(array_i1d_obj)                      :: row_blk_sizes

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_row_block_sizes', &
      routineP = moduleN//':'//routineN

    row_blk_sizes = matrix%cp_row_sizes
  END FUNCTION cp_dbcsr_row_block_sizes

  FUNCTION cp_dbcsr_col_block_sizes (matrix) RESULT (col_blk_sizes)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    TYPE(array_i1d_obj)                      :: col_blk_sizes

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_col_block_sizes', &
      routineP = moduleN//':'//routineN

    col_blk_sizes = matrix%cp_col_sizes
  END FUNCTION cp_dbcsr_col_block_sizes

  FUNCTION cp_dbcsr_col_block_offsets (matrix) RESULT (col_blk_offsets)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    TYPE(array_i1d_obj)                      :: col_blk_offsets

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_col_block_offsets', &
      routineP = moduleN//':'//routineN

    col_blk_offsets = matrix%cp_col_offsets
  END FUNCTION cp_dbcsr_col_block_offsets

  FUNCTION cp_dbcsr_row_block_offsets (matrix) RESULT (row_blk_offsets)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    TYPE(array_i1d_obj)                      :: row_blk_offsets

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_row_block_offsets', &
      routineP = moduleN//':'//routineN

    row_blk_offsets = matrix%cp_row_offsets
  END FUNCTION cp_dbcsr_row_block_offsets


  elemental FUNCTION cp_dbcsr_has_mapping (map) RESULT (has_mapping)
    TYPE(cp_dbcsr_mapping_type), INTENT(in)  :: map
    LOGICAL                                  :: has_mapping

    has_mapping = map%defined .AND. map%active
  END FUNCTION cp_dbcsr_has_mapping

  elemental FUNCTION cp_dbcsr_mapping_defined (map) RESULT (mapping_defined)
    TYPE(cp_dbcsr_mapping_type), INTENT(in)  :: map
    LOGICAL                                  :: mapping_defined

    mapping_defined = map%defined
  END FUNCTION cp_dbcsr_mapping_defined


  SUBROUTINE cp_dbcsr_mapping_create (map, clusters, error)
    TYPE(cp_dbcsr_mapping_type), POINTER     :: map
    TYPE(cluster_type), INTENT(IN), OPTIONAL :: clusters
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_mapping_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nclusters, nparticles, stat
    LOGICAL                                  :: failure

!   ---------------------------------------------------------------------------

    failure = .TRUE.
    CALL cp_assert (.NOT. ASSOCIATED (map), cp_warning_level,&
         cp_wrong_args_error, routineN,&
         "Map is associated or not null.", error=error)
    ALLOCATE (map, stat=stat)
    CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
    map%refcount = 1
    IF (PRESENT (clusters)) THEN
       nparticles = clusters%nparticles
       nclusters = clusters%nclusters
       map%defined = .TRUE.
       map%active = .FALSE.
       map%clusters = clusters
       NULLIFY (map%clusters%cluster_ptr)
       NULLIFY (map%clusters%particle_idx)
       NULLIFY (map%clusters%particle_to_cluster)
       CALL cluster_allocate (map%clusters, nclusters,&
            nparticles, error=error)
       CALL cp_assert (SIZE (clusters%cluster_ptr)-1 .GE. clusters%nclusters,&
            cp_fatal_level, cp_internal_error, routineN,&
            "Cluster index size does not match number of clusters.",&
            error=error)
       map%clusters%cluster_ptr(1:nclusters+1)&
            = clusters%cluster_ptr(1:nclusters+1)
       map%clusters%particle_idx(1:nparticles) = clusters%particle_idx(1:nparticles)
       map%clusters%particle_to_cluster(1:nparticles) = clusters%particle_to_cluster(1:nparticles)
       !@@@
       !write(*,*)routineN//" The mapping arrays"
       !write(*,'(10(1X,I5)))')map%clusters%cluster_ptr(1:nclusters+1)
       !write(*,'(10(1X,I5)))')map%clusters%particle_idx(1:nparticles)
       !write(*,'(10(1X,I5)))')map%clusters%particle_to_cluster(1:nparticles)
    ELSE
       map%defined = .FALSE.
       map%active = .FALSE.
    ENDIF
  END SUBROUTINE cp_dbcsr_mapping_create

  SUBROUTINE cp_dbcsr_mapping_create_idem (map, nparticles, error)
    TYPE(cp_dbcsr_mapping_type), POINTER     :: map
    INTEGER, INTENT(IN)                      :: nparticles
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_mapping_create_idem', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

!   ---------------------------------------------------------------------------

    failure = .TRUE.
    CALL cp_assert (.NOT. ASSOCIATED (map), cp_warning_level,&
         cp_wrong_args_error, routineN,&
         "Map is associated or not null.", error=error)
    ALLOCATE (map, stat=stat)
    CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
    map%refcount = 1
    map%defined = .TRUE.
    map%active = .FALSE.
    NULLIFY (map%clusters%cluster_ptr)
    NULLIFY (map%clusters%particle_idx)
    NULLIFY (map%clusters%particle_to_cluster)
    CALL cluster_allocate (map%clusters, nparticles,&
         nparticles, error=error)
    map%clusters%cluster_ptr(1:nparticles) = (/(i,i=1,nparticles+1)/)
    map%clusters%particle_idx(1:nparticles) = (/(i,i=1,nparticles+1)/)
    map%clusters%particle_to_cluster(1:nparticles) = (/(i,i=1,nparticles+1)/)
  END SUBROUTINE cp_dbcsr_mapping_create_idem


  SUBROUTINE cp_dbcsr_mapping_release (map, error)
    TYPE(cp_dbcsr_mapping_type), POINTER     :: map
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_mapping_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

!   ---------------------------------------------------------------------------
!CALL cp_assert (ASSOCIATED (map), cp_warning_level, cp_wrong_args_error,&
!     routineN, "Map does not exist.", error=error)

    IF (ASSOCIATED (map)) THEN
       map%refcount = map%refcount - 1
       IF (map%refcount .LE. 0) THEN
          IF (map%defined) THEN
             CALL cluster_deallocate (map%clusters, error=error)
             map%defined = .FALSE.
          ENDIF
          failure = .TRUE.
          DEALLOCATE (map, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          NULLIFY (map)
       ENDIF
    ENDIF
  END SUBROUTINE cp_dbcsr_mapping_release

  SUBROUTINE cp_dbcsr_mapping_hold (map, error)
    TYPE(cp_dbcsr_mapping_type), POINTER     :: map
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_mapping_hold', &
      routineP = moduleN//':'//routineN

    CALL cp_assert (ASSOCIATED (map), cp_fatal_level, cp_wrong_args_error,&
         routineN, "Map does not exist.", error=error)
    map%refcount = map%refcount + 1
  END SUBROUTINE cp_dbcsr_mapping_hold

  SUBROUTINE cp_dbcsr_mapping_activate (old_map, new_map, activation, error)
    TYPE(cp_dbcsr_mapping_type), POINTER     :: old_map, new_map
    LOGICAL, INTENT(IN)                      :: activation
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_mapping_activate', &
      routineP = moduleN//':'//routineN

    CALL cp_assert (ASSOCIATED (old_map), cp_fatal_level, cp_wrong_args_error,&
         routineN, "Map does not exist.", error=error)
    NULLIFY (new_map)
    CALL cp_dbcsr_mapping_create (new_map, old_map%clusters, error=error)
    new_map%active = activation
  END SUBROUTINE cp_dbcsr_mapping_activate

END MODULE cp_dbcsr_methods
