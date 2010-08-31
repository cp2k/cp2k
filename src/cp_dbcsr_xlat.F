!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Utility for clustering translations and mapping.
!> \author  Urban Borstnik
!> \date    2010-08-11
!> \version 0.1
!>
!> <b>Modification history:</b>
!> - Created 2010-08-11 [UB]
! *****************************************************************************
MODULE cp_dbcsr_xlat
  USE array_types,                     ONLY: array_data,&
                                             array_hold,&
                                             array_i1d_obj,&
                                             array_new,&
                                             array_nullify,&
                                             array_release
  USE cluster_buffers,                 ONLY: cluster_buffer_flush,&
                                             cluster_buffer_set_ptr
  USE cluster_methods,                 ONLY: cluster_get_nclusters,&
                                             cluster_get_nparticles,&
                                             cluster_member_particles,&
                                             cluster_part2cluster
  USE cp_dbcsr_methods,                ONLY: cp_dbcsr_col_block_offsets,&
                                             cp_dbcsr_col_block_sizes,&
                                             cp_dbcsr_has_mapping,&
                                             cp_dbcsr_mapping_defined,&
                                             cp_dbcsr_row_block_offsets,&
                                             cp_dbcsr_row_block_sizes
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_mapping_type,&
                                             cp_dbcsr_type
  USE dbcsr_block_buffers,             ONLY: dbcsr_buffers_mark_dirty,&
                                             dbcsr_buffers_resize
  USE dbcsr_error_handling,            ONLY: dbcsr_error_type
  USE dbcsr_methods,                   ONLY: dbcsr_col_block_sizes,&
                                             dbcsr_distribution_col_dist,&
                                             dbcsr_distribution_init,&
                                             dbcsr_distribution_mp,&
                                             dbcsr_distribution_new,&
                                             dbcsr_distribution_row_dist,&
                                             dbcsr_row_block_sizes
  USE dbcsr_ptr_util,                  ONLY: dbcsr_ptr_remapping,&
                                             pointer_d_rank_remap2
  USE dbcsr_types,                     ONLY: dbcsr_distribution_obj
  USE dbcsr_util,                      ONLY: swap
  USE kinds,                           ONLY: int_size,&
                                             real_8
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop

  !$ USE OMP_LIB
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dbcsr_xlat'

  LOGICAL, PRIVATE, PARAMETER :: careful_mod = .TRUE.

  PUBLIC :: set_block_pointer_to_cluster
  PUBLIC :: get4blksizes
  PUBLIC :: blksizes_blk_to_cluster, dist_blk_to_cluster
  PUBLIC :: blksizes_blk_to_cluster4

CONTAINS

! *****************************************************************************
!> \brief Sets a block pointer to a region of a clustered block.
!> \par Notes
!>      This routine does all the necessary lookups of blocks sizes
!>      and offsets and then calls another routine for the actual
!>      pointer assignments.
!> \par History
!>      - created 2010-08-10
!> \author UB
! *****************************************************************************
  SUBROUTINE set_block_pointer_to_cluster(matrix,&
       high_row, high_col, high_block, high_tr, contiguous_pointers, read_only,&
       low_block, low_row, low_col, low_tr,&
       high_row_size, high_col_size, low_col_size, low_row_size,&
       error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    INTEGER, INTENT(IN)                      :: high_row, high_col
    REAL(kind=real_8), DIMENSION(:, :), &
      POINTER                                :: high_block
    LOGICAL, INTENT(IN)                      :: high_tr, contiguous_pointers, &
                                                read_only
    REAL(kind=real_8), DIMENSION(:), POINTER :: low_block
    INTEGER, INTENT(IN), OPTIONAL            :: low_row, low_col
    LOGICAL, INTENT(IN)                      :: low_tr
    INTEGER, INTENT(IN), OPTIONAL            :: high_row_size, high_col_size, &
                                                low_col_size, low_row_size
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set_block_pointer_to_cluster', &
      routineP = moduleN//':'//routineN

    INTEGER :: actual_cs, actual_rs, col_offset, handle, high_cs, high_rs, &
      low_c, low_cs, low_r, low_rs, row_offset
    LOGICAL                                  :: failure
    REAL(kind=real_8), DIMENSION(:, :), &
      POINTER                                :: low_block_2d
    TYPE(dbcsr_error_type)                   :: dbcsr_error

!   ---------------------------------------------------------------------------

    IF (careful_mod) CALL timeset (routineN, handle)
    CALL get4blksizes (matrix, high_row, high_col, low_row, low_col,&
       high_row_size, high_col_size, low_row_size, low_col_size,&
       high_rs, high_cs, low_r, low_c, low_rs, low_cs, row_offset, col_offset,&
       error=error)
    !
    CALL cluster_buffer_flush (matrix%buffers, error=error)
    actual_rs = high_rs
    actual_cs = high_cs
    IF (high_tr) CALL swap (actual_rs, actual_cs)
    IF (.NOT. contiguous_pointers .AND. dbcsr_ptr_remapping) THEN
       ! point high pointer into cluster block
       !high_block => low_block(low_row_offset, low_col_offset)
       CALL pointer_d_rank_remap2 (low_block_2d,&
            actual_rs, actual_cs, low_block)
       high_block => low_block_2d(row_offset:row_offset+high_rs-1,&
                     col_offset:col_offset+high_cs-1)
    ELSE
       ! Resize the buffers to correct row size
       CALL dbcsr_buffers_resize (matrix%buffers,&
            row_size=actual_rs, col_size=actual_cs,&
            minimum=.NOT.contiguous_pointers, error=dbcsr_error)
       CALL cluster_buffer_set_ptr (high_block, high_row, high_col,&
            high_rs, high_cs, high_tr, &
            low_rs, low_cs, low_tr, row_offset, col_offset,&
            matrix%buffers, error=error)
       IF (read_only) THEN
          CALL dbcsr_buffers_mark_dirty (matrix%buffers, dirty=.FALSE.,&
               error=dbcsr_error)
       ENDIF
    ENDIF
    IF (careful_mod) CALL timestop (handle)
  END SUBROUTINE set_block_pointer_to_cluster


  SUBROUTINE get4blksizes (matrix, high_row, high_col, low_row, low_col,&
       high_row_size, high_col_size,&
       low_row_size, low_col_size,&
       high_rs, high_cs, low_r, low_c, low_rs, low_cs, row_off, col_off, error)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER, INTENT(IN)                      :: high_row, high_col
    INTEGER, INTENT(IN), OPTIONAL            :: low_row, low_col, &
                                                high_row_size, high_col_size, &
                                                low_row_size, low_col_size
    INTEGER, INTENT(OUT), OPTIONAL           :: high_rs, high_cs
    INTEGER, INTENT(OUT)                     :: low_r, low_c, low_rs, low_cs, &
                                                row_off, col_off
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get4blksizes', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    INTEGER, DIMENSION(:), POINTER           :: cbs, rbs

!   ---------------------------------------------------------------------------

    IF (careful_mod) CALL timeset (routineN, handle)
    IF (PRESENT (high_rs)) THEN
       IF (.NOT. PRESENT (high_row_size)) THEN
          rbs => array_data (cp_dbcsr_row_block_sizes (matrix))
          high_rs = rbs(high_row)
       ELSE
          high_rs = high_row_size
       ENDIF
    ENDIF
    IF (PRESENT (high_cs)) THEN
       IF (.NOT. PRESENT (high_col_size)) THEN
          cbs => array_data (cp_dbcsr_col_block_sizes (matrix))
          high_cs = cbs(high_col)
       ELSE
          high_cs = high_col_size
       ENDIF
    ENDIF
    !
    IF (.NOT. PRESENT (low_row)) THEN
       IF (careful_mod) THEN
          CALL cp_assert (ASSOCIATED (matrix%row_mapping), cp_fatal_level,&
               cp_internal_error, routineN, "Row mapping nonexistent.",&
               error=error)
       ENDIF
       low_r = cp_part2cluster (matrix%row_mapping, high_row, error)
    ELSE
       low_r = low_row
    ENDIF
    IF (.NOT. PRESENT (low_col)) THEN
       IF (careful_mod) THEN
          CALL cp_assert (ASSOCIATED (matrix%col_mapping), cp_fatal_level,&
               cp_internal_error, routineN, "Column mapping nonexistent.",&
               error=error)
       ENDIF
       low_c = cp_part2cluster (matrix%col_mapping, high_col, error)
    ELSE
       low_c = low_col
    ENDIF
    !
    IF (.NOT. PRESENT (low_row_size)) THEN
       rbs => array_data (dbcsr_row_block_sizes (matrix%matrix))
       low_rs = rbs(low_r)
    ELSE
       low_rs = low_row_size
    ENDIF
    IF (.NOT. PRESENT (low_col_size)) THEN
       cbs => array_data (dbcsr_col_block_sizes (matrix%matrix))
       low_cs = cbs(low_c)
    ELSE
       low_cs = low_col_size
    ENDIF
    !
    row_off = get_offset_in_cluster (matrix%row_mapping, high_row,&
         cp_dbcsr_row_block_offsets (matrix), error=error)
    !
    col_off = get_offset_in_cluster (matrix%col_mapping, high_col,&
         cp_dbcsr_col_block_offsets (matrix), error=error) 
    IF (careful_mod) CALL timestop (handle)
  END SUBROUTINE get4blksizes

  FUNCTION get_offset_in_cluster (mapping, particle, particle_offsets, error)&
       RESULT (cluster_offset)
    TYPE(cp_dbcsr_mapping_type), INTENT(IN)  :: mapping
    INTEGER, INTENT(IN)                      :: particle
    TYPE(array_i1d_obj), INTENT(IN)          :: particle_offsets
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    INTEGER                                  :: cluster_offset

    CHARACTER(len=*), PARAMETER :: routineN = 'get_offset_in_cluster', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: cluster, handle, i, p
    INTEGER, DIMENSION(:), POINTER           :: cluster_members, p_offsets

!   ---------------------------------------------------------------------------

    IF (careful_mod)  CALL timeset (routineN, handle)
    ! Checks offset for all particles in the same cluster and returns
    ! their minimum.
    IF (cp_dbcsr_has_mapping (mapping)) THEN
       cluster = cp_part2cluster (mapping, particle, error)
       cluster_members =>&
            cluster_member_particles (mapping%clusters, cluster, error)
       p_offsets => array_data (particle_offsets)
       cluster_offset = 1
       DO i = 1, SIZE (cluster_members)
          p = cluster_members(i)
          IF (p .EQ. particle) THEN
             EXIT
          ELSE
             cluster_offset = cluster_offset + &
                  p_offsets(p+1) - p_offsets(p)
          ENDIF
       ENDDO
    ELSE
       cluster_offset = 1
    ENDIF
    IF (careful_mod) CALL timestop (handle)
  END FUNCTION get_offset_in_cluster


  FUNCTION cp_part2cluster (mapping, particle, error) RESULT (cluster)
    TYPE(cp_dbcsr_mapping_type), INTENT(IN)  :: mapping
    INTEGER, INTENT(IN)                      :: particle
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    INTEGER                                  :: cluster

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_part2cluster', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

!   ---------------------------------------------------------------------------

    IF (careful_mod) CALL timeset (routineN, handle)
    IF (cp_dbcsr_has_mapping (mapping)) THEN
       IF (careful_mod) THEN
          CALL cp_assert (particle .GE. 1,&
               cp_fatal_level, cp_wrong_args_error, routineN,&
               "Particle number less than 1.", error=error)
          CALL cp_assert (particle .LE.&
               cluster_get_nparticles (mapping%clusters),&
               cp_fatal_level, cp_wrong_args_error, routineN,&
               "Particle greater than declared number of particles.",&
               error=error)
       ENDIF
       cluster = cluster_part2cluster (mapping%clusters, particle, error)
    ELSE
       cluster = particle
    ENDIF
    IF (careful_mod) CALL timestop (handle)
  END FUNCTION cp_part2cluster


  SUBROUTINE blksizes_blk_to_cluster(block_sizes,&
       mapping, cluster_sizes, error)
    TYPE(array_i1d_obj), INTENT(IN)          :: block_sizes
    TYPE(cp_dbcsr_mapping_type), INTENT(IN)  :: mapping
    TYPE(array_i1d_obj), INTENT(OUT)         :: cluster_sizes
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'blksizes_blk_to_cluster', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: cluster, handle, ncl, np, &
                                                particle, stat
    INTEGER, DIMENSION(:), POINTER           :: blk_s, cl_s

!   ---------------------------------------------------------------------------

    CALL timeset (routineN, handle)
    IF (cp_dbcsr_mapping_defined (mapping)) THEN
       np = cluster_get_nparticles (mapping%clusters)
       ncl = cluster_get_nclusters (mapping%clusters)
       ALLOCATE (cl_s (ncl), stat=stat)
       IF (stat /= 0) CALL stop_memory (routineN, "cl_s", int_size*ncl)
       CALL array_new (cluster_sizes, cl_s, gift=.TRUE.)
       blk_s => array_data (block_sizes)
       cl_s(:) = 0
       DO particle = 1, np
          cluster = cluster_part2cluster (mapping%clusters, particle, error)
          cl_s(cluster) = cl_s(cluster) + blk_s(particle)
       ENDDO
       !forall (cluster = 1 : ncl)
       !   !diff = mod(old, 4)
       !   !cl_s = old + mod(4-diff,4)
       !   cl_s(cluster) = cl_s(cluster) + MOD(4-MOD(cl_s(cluster),4), 4)
       !end forall
    ELSE
       cluster_sizes = block_sizes
       CALL array_hold (cluster_sizes)
    ENDIF
    CALL timestop (handle)
  END SUBROUTINE blksizes_blk_to_cluster
  SUBROUTINE blksizes_blk_to_cluster4(block_sizes,&
       mapping, cluster_sizes, error)
    TYPE(array_i1d_obj), INTENT(IN)          :: block_sizes
    TYPE(cp_dbcsr_mapping_type), INTENT(IN)  :: mapping
    TYPE(array_i1d_obj), INTENT(OUT)         :: cluster_sizes
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'blksizes_blk_to_cluster4', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: cluster, handle, ncl, np, &
                                                particle, stat
    INTEGER, DIMENSION(:), POINTER           :: blk_s, cl_s

!   ---------------------------------------------------------------------------

    CALL timeset (routineN, handle)
    IF (cp_dbcsr_mapping_defined (mapping)) THEN
       np = cluster_get_nparticles (mapping%clusters)
       ncl = cluster_get_nclusters (mapping%clusters)
       ALLOCATE (cl_s (ncl), stat=stat)
       IF (stat /= 0) CALL stop_memory (routineN, "cl_s", int_size*ncl)
       CALL array_new (cluster_sizes, cl_s, gift=.TRUE.)
       blk_s => array_data (block_sizes)
       cl_s(:) = 0
       DO particle = 1, np
          cluster = cluster_part2cluster (mapping%clusters, particle, error)
          cl_s(cluster) = cl_s(cluster) + blk_s(particle)
       ENDDO
       FORALL (cluster = 1 : ncl)
          !diff = mod(old, 4)
          !cl_s = old + mod(4-diff,4)
          cl_s(cluster) = cl_s(cluster) + MOD(4-MOD(cl_s(cluster),4), 4)
       END FORALL
    ELSE
       cluster_sizes = block_sizes
       CALL array_hold (cluster_sizes)
    ENDIF
    CALL timestop (handle)
  END SUBROUTINE blksizes_blk_to_cluster4

  SUBROUTINE dist_blk_to_cluster (block_dist, row_mapping, col_mapping,&
       cluster_dist, error)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: block_dist
    TYPE(cp_dbcsr_mapping_type), INTENT(IN)  :: row_mapping, col_mapping
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT)                            :: cluster_dist
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dist_blk_to_cluster', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: cl, handle, ncl_col, ncl_row, &
                                                np_col, np_row, stat
    INTEGER, DIMENSION(:), POINTER           :: blk_col_dist, blk_row_dist, &
                                                cl_col_dist, cl_row_dist, &
                                                members
    TYPE(array_i1d_obj)                      :: ccd, crd

!   ---------------------------------------------------------------------------

    CALL timeset (routineN, handle)
    IF (cp_dbcsr_mapping_defined (row_mapping)) THEN
       np_row = cluster_get_nparticles (row_mapping%clusters)
       ncl_row = cluster_get_nclusters (row_mapping%clusters)
       ALLOCATE (cl_row_dist (ncl_row), stat=stat)
       IF (stat /= 0) CALL stop_memory (routineN, "cl_row_dist", int_size*ncl_row)
       blk_row_dist => array_data(dbcsr_distribution_row_dist (block_dist))
       !
       DO cl = 1, ncl_row
          members => cluster_member_particles (row_mapping%clusters, cl, error)
          IF (SIZE (members) .GE. 1) THEN
             cl_row_dist(cl) = blk_row_dist (members(1))
          ELSE
             cl_row_dist(cl) = 0
          ENDIF
       ENDDO
       CALL array_nullify (crd)
       CALL array_new (crd, cl_row_dist, gift=.TRUE.)
    ELSE
       crd = dbcsr_distribution_row_dist (block_dist)
       CALL array_hold (crd)
    ENDIF
    IF (cp_dbcsr_mapping_defined (col_mapping)) THEN
       np_col = cluster_get_nparticles (col_mapping%clusters)
       ncl_col = cluster_get_nclusters (col_mapping%clusters)
       ALLOCATE (cl_col_dist (ncl_col), stat=stat)
       IF (stat /= 0) CALL stop_memory (routineN, "cl_col_dist", int_size*ncl_col)
       blk_col_dist => array_data(dbcsr_distribution_col_dist (block_dist))
       DO cl = 1, ncl_col
          members => cluster_member_particles (col_mapping%clusters, cl, error)
          IF (SIZE (members) .GE. 1) THEN
             cl_col_dist(cl) = blk_col_dist (members(1))
          ELSE
             cl_col_dist(cl) = 0
          ENDIF
       ENDDO
       CALL array_nullify (ccd)
       CALL array_new (ccd, cl_col_dist, gift=.TRUE.)
    ELSE
       ccd = dbcsr_distribution_col_dist (block_dist)
       CALL array_hold (ccd)
    ENDIF
    CALL dbcsr_distribution_init (cluster_dist)
    CALL dbcsr_distribution_new (cluster_dist,&
         dbcsr_distribution_mp (block_dist), crd, ccd)
    CALL array_release (crd)
    CALL array_release (ccd)
    CALL timestop (handle)
  END SUBROUTINE dist_blk_to_cluster

END MODULE cp_dbcsr_xlat
