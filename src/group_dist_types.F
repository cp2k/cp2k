!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Types to describe group distributions
!> \par History
!>       2019.03 created [Frederick Stein]
!> \author Frederick Stein
! **************************************************************************************************
MODULE group_dist_types
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE message_passing,                 ONLY: mp_sum
   USE util,                            ONLY: get_limit
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'group_dist_types'

   PUBLIC :: group_dist_d0_type, group_dist_d1_type, create_group_dist, get_group_dist, release_group_dist

   ! Type to represent start, end and size of one group (usually the own group)
   ! The underscore is added to prevent confusions with fortran intrinsics
   TYPE group_dist_d0_type
      INTEGER :: starts, ends, sizes
   END TYPE group_dist_d0_type

   ! Class to represent 1d-sets of starts, ends and sizes of all groups
   ! Better for exchanging information than using arrays of group_dist_d0_type
   TYPE group_dist_d1_type
      INTEGER, ALLOCATABLE, DIMENSION(:) :: starts, ends, sizes
   END TYPE group_dist_d1_type
   
   INTERFACE create_group_dist
      MODULE PROCEDURE create_group_dist_d1_i1, &
                       create_group_dist_d1_i3, &
                       create_group_dist_d1_gd
   END INTERFACE create_group_dist
   
   INTERFACE get_group_dist
      MODULE PROCEDURE get_group_dist_d1_gd, &
                       get_group_dist_d1_i !, &
!                        get_group_dist_d1_ia
   END INTERFACE

CONTAINS

   SUBROUTINE create_group_dist_d1_i1(this, ngroups, dimen)
      TYPE(group_dist_d1_type), INTENT(INOUT) :: this
      INTEGER, INTENT(IN)                   :: ngroups, dimen

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_group_dist_d1_i1', routineP = moduleN//':'//routineN
   
      INTEGER :: handle, iproc
      INTEGER, DIMENSION(2) :: itmp
      
      CALL timeset(routineN, handle)
      
      ALLOCATE(this%starts(0:ngroups-1))
      this%starts = 0
      ALLOCATE(this%ends(0:ngroups-1))
      this%ends = 0
      ALLOCATE(this%sizes(0:ngroups-1))
      this%sizes = 0

      DO iproc = 0, ngroups-1
         itmp = get_limit(dimen, ngroups, iproc)
         this%starts(iproc) = itmp(1)
         this%ends(iproc) = itmp(2)
         this%sizes(iproc) = itmp(2)-itmp(1)+1
      END DO
      
      CALL timestop(handle)
      
   END SUBROUTINE create_group_dist_d1_i1

   SUBROUTINE create_group_dist_d1_i3(this, ngroups, starts, ends, sizes, para_env)
      TYPE(group_dist_d1_type), INTENT(INOUT) :: this
      INTEGER, INTENT(IN)                   :: ngroups, starts, ends, sizes
      TYPE(cp_para_env_type), POINTER       :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_group_dist_d1_i3', routineP = moduleN//':'//routineN
      
      INTEGER :: handle
      
      CALL timeset(routineN, handle)
      
      ALLOCATE(this%starts(0:ngroups-1))
      this%starts = 0
      ALLOCATE(this%ends(0:ngroups-1))
      this%ends = 0
      ALLOCATE(this%sizes(0:ngroups-1))
      this%sizes = 0

      this%starts(para_env%num_pe) = starts
      this%ends(para_env%num_pe) = ends
      this%sizes(para_env%num_pe) = sizes
      
      CALL mp_sum(this%starts, para_env%group)
      CALL mp_sum(this%ends, para_env%group)
      CALL mp_sum(this%sizes, para_env%group)
      
      CALL timestop(handle)
      
   END SUBROUTINE create_group_dist_d1_i3

   SUBROUTINE create_group_dist_d1_gd(this, dimen, group_dist_ext, para_env, offset)
      TYPE(group_dist_d1_type), INTENT(INOUT) :: this
      INTEGER, INTENT(IN)                   :: dimen
      TYPE(group_dist_d0_type), INTENT(IN)  :: group_dist_ext
      TYPE(cp_para_env_type), POINTER       :: para_env
      INTEGER, INTENT(IN), OPTIONAL         :: offset

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_group_dist_d1_gd', routineP = moduleN//':'//routineN
      
      INTEGER :: handle, my_offset
      
      CALL timeset(routineN, handle)
      
      my_offset = 0
      IF (PRESENT(offset)) my_offset = offset
      
      ALLOCATE(this%starts(my_offset:dimen+my_offset-1))
      this%starts = 0
      ALLOCATE(this%ends(my_offset:dimen+my_offset-1))
      this%ends = 0
      ALLOCATE(this%sizes(my_offset:dimen+my_offset-1))
      this%sizes = 0

      this%starts(para_env%num_pe) = group_dist_ext%starts
      this%ends(para_env%num_pe) = group_dist_ext%ends
      this%sizes(para_env%num_pe) = group_dist_ext%sizes
      
      CALL mp_sum(this%starts, para_env%group)
      CALL mp_sum(this%ends, para_env%group)
      CALL mp_sum(this%sizes, para_env%group)
      
      CALL timestop(handle)
      
   END SUBROUTINE create_group_dist_d1_gd

   PURE SUBROUTINE get_group_dist_d1_gd(this, pos, group_dist_ext)
      TYPE(group_dist_d1_type), INTENT(INOUT) :: this
      INTEGER, INTENT(IN)                   :: pos
      TYPE(group_dist_d0_type), INTENT(OUT) :: group_dist_ext

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_group_dist_d1_gd', routineP = moduleN//':'//routineN

      group_dist_ext%starts = this%starts(pos)
      group_dist_ext%ends = this%ends(pos)
      group_dist_ext%sizes = this%sizes(pos)
      
   END SUBROUTINE get_group_dist_d1_gd

   PURE SUBROUTINE get_group_dist_d1_i(this, pos, starts, ends, sizes)
      TYPE(group_dist_d1_type), INTENT(INOUT) :: this
      INTEGER, INTENT(IN)                   :: pos
      INTEGER, INTENT(OUT), OPTIONAL        :: starts, ends, sizes

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_group_dist_d1_i', routineP = moduleN//':'//routineN

      IF (PRESENT(starts)) starts = this%starts(pos)
      IF (PRESENT(ends)) ends = this%ends(pos)
      IF (PRESENT(sizes)) sizes = this%sizes(pos)
      
   END SUBROUTINE get_group_dist_d1_i

!    PURE SUBROUTINE get_group_dist_d1_ia(this, start_a, end_a, size_a)
!       TYPE(group_dist_d1_type), INTENT(INOUT) :: this
!       INTEGER, DIMENSION(:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: start_a, end_a, size_a
! 
!       CHARACTER(LEN=*), PARAMETER :: routineN = 'get_group_dist_d1_ia', routineP = moduleN//':'//routineN
! 
!       IF (PRESENT(start_a)) THEN
!          IF (ALLOCATED(start_a)) THEN
!             IF (SIZE(start_a) /= SIZE(this%start_)) THEN
!                DEALLOCATE(start_a)
!             END IF
!          END IF
!          IF (.NOT. ALLOCATED(start_a)) ALLOCATE(start_a(SIZE(this%start_)))
!          start_a(:) = this%start_(:)
!       END IF
!       
!       IF (PRESENT(end_a)) THEN
!          IF (ALLOCATED(end_a)) THEN
!             IF (SIZE(end_a) /= SIZE(this%end_)) THEN
!                DEALLOCATE(end_a)
!             END IF
!          END IF
!          IF (.NOT. ALLOCATED(end_a)) ALLOCATE(end_a(SIZE(this%end_)))
!          end_a(:) = this%end_(:)
!       END IF
!       
!       IF (PRESENT(size_a)) THEN
!          IF (ALLOCATED(size_a)) THEN
!             IF (SIZE(size_a) /= SIZE(this%size_)) THEN
!                DEALLOCATE(size_a)
!             END IF
!          END IF
!          IF (.NOT. ALLOCATED(size_a)) ALLOCATE(size_a(SIZE(this%size_)))
!          size_a(:) = this%size_(:)
!       END IF
!       
!    END SUBROUTINE get_group_dist_d1_ia

   PURE SUBROUTINE release_group_dist(this)
      CLASS(group_dist_d1_type), INTENT(INOUT) :: this

      CHARACTER(LEN=*), PARAMETER :: routineN = 'release_group_dist', routineP = moduleN//':'//routineN

      DEALLOCATE(this%starts)
      DEALLOCATE(this%ends)
      DEALLOCATE(this%sizes)
      
   END SUBROUTINE release_group_dist

END MODULE group_dist_types
