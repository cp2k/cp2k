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
  USE dbcsr_methods,                   ONLY: dbcsr_col_block_offsets,&
                                             dbcsr_col_block_sizes,&
                                             dbcsr_row_block_offsets,&
                                             dbcsr_row_block_sizes

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
  PUBLIC :: cp_dbcsr_has_mapping


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

    row_blk_sizes = dbcsr_row_block_sizes (matrix%matrix)
  END FUNCTION cp_dbcsr_row_block_sizes

  FUNCTION cp_dbcsr_col_block_sizes (matrix) RESULT (col_blk_sizes)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    TYPE(array_i1d_obj)                      :: col_blk_sizes

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_col_block_sizes', &
      routineP = moduleN//':'//routineN

    col_blk_sizes = dbcsr_col_block_sizes (matrix%matrix)
  END FUNCTION cp_dbcsr_col_block_sizes

  FUNCTION cp_dbcsr_col_block_offsets (matrix) RESULT (col_blk_offsets)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    TYPE(array_i1d_obj)                      :: col_blk_offsets

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_col_block_offsets', &
      routineP = moduleN//':'//routineN

    col_blk_offsets = dbcsr_col_block_offsets (matrix%matrix)
  END FUNCTION cp_dbcsr_col_block_offsets

  FUNCTION cp_dbcsr_row_block_offsets (matrix) RESULT (row_blk_offsets)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    TYPE(array_i1d_obj)                      :: row_blk_offsets

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_row_block_offsets', &
      routineP = moduleN//':'//routineN

    row_blk_offsets = dbcsr_row_block_offsets (matrix%matrix)
  END FUNCTION cp_dbcsr_row_block_offsets


  elemental FUNCTION cp_dbcsr_has_mapping (map) RESULT (has_mapping)
    TYPE(cp_dbcsr_mapping_type), INTENT(in)  :: map
    LOGICAL                                  :: has_mapping

    has_mapping = map%defined .AND. map%active
  END FUNCTION cp_dbcsr_has_mapping


  SUBROUTINE cp_dbcsr_mapping_create (map, clusters, error)
    TYPE(cp_dbcsr_mapping_type), POINTER     :: map
    TYPE(cluster_type), INTENT(IN), OPTIONAL :: clusters
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_mapping_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .TRUE.
    CALL cp_assert (.NOT. ASSOCIATED (map), cp_warning_level,&
         cp_wrong_args_error, routineN,&
         "Map is associated or not null.", error=error)
    ALLOCATE (map, stat=stat)
    CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
    map%refcount = 1
    IF (PRESENT (clusters)) THEN
       map%defined = .TRUE.
       map%active = .FALSE.
       map%clusters = clusters
       NULLIFY (map%clusters%cluster_ptr)
       NULLIFY (map%clusters%particle_idx)
       NULLIFY (map%clusters%particle_to_cluster)
       CALL cluster_allocate (map%clusters, clusters%nclusters,&
            clusters%nparticles, error=error)
       map%clusters%cluster_ptr(:) = clusters%cluster_ptr(:)
       map%clusters%particle_idx(:) = clusters%particle_idx(:)
       !STOP "check that p2c is actually allocated now"
       !map%clusters%particle_to_cluster(:) = clusters%particle_to_particle(:)
    ELSE
       map%defined = .FALSE.
       map%active = .FALSE.
    ENDIF
  END SUBROUTINE cp_dbcsr_mapping_create

  SUBROUTINE cp_dbcsr_mapping_release (map, error)
    TYPE(cp_dbcsr_mapping_type), POINTER     :: map
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_mapping_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    !CALL cp_assert (ASSOCIATED (map), cp_warning_level, cp_wrong_args_error,&
    !     routineN, "Map does not exist.", error=error)
    IF (ASSOCIATED (map)) THEN
       map%refcount = map%refcount - 1
       IF (map%refcount .LE. 0) THEN
          IF (map%defined) THEN
             CALL cluster_deallocate (map%clusters, error=error)
             !STOP "check that p2c is actually deallocated"
             map%defined = .FALSE.
          ENDIF
       ENDIF
       failure = .TRUE.
       DEALLOCATE (map, stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       NULLIFY (map)
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



END MODULE cp_dbcsr_methods
