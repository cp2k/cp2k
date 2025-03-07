!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Types to describe group distributions
!> \par History
!>       2019.03 created [Frederick Stein]
!> \author Frederick Stein
! **************************************************************************************************
MODULE group_dist_types
   USE message_passing,                 ONLY: mp_comm_type
   USE util,                            ONLY: get_limit
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'group_dist_types'

   PUBLIC :: group_dist_d0_type, group_dist_d1_type, &
             create_group_dist, get_group_dist, release_group_dist, maxsize, group_dist_proc

   ! Type to represent start, end and size of one group (usually the own group)
   TYPE group_dist_d0_type
      INTEGER :: starts = -1, ends = -1, sizes = -1
   END TYPE group_dist_d0_type

   ! Class to represent 1d-sets of starts, ends and sizes of all groups
   ! Better for exchanging information than using arrays of group_dist_d0_type
   TYPE group_dist_d1_type
      INTEGER, ALLOCATABLE, DIMENSION(:) :: starts, ends, sizes
   END TYPE group_dist_d1_type

   INTERFACE create_group_dist
      MODULE PROCEDURE create_group_dist_d0, &
         create_group_dist_d1_i1, &
         create_group_dist_d1_i3, &
         create_group_dist_d1_gd, &
         create_group_dist_d1_0
   END INTERFACE create_group_dist

   INTERFACE get_group_dist
      MODULE PROCEDURE get_group_dist_d0, &
         get_group_dist_d1, &
         get_group_dist_gd1
   END INTERFACE

   INTERFACE release_group_dist
      MODULE PROCEDURE release_group_dist_d1
   END INTERFACE release_group_dist

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param this ...
!> \param ngroups ...
!> \param dimen ...
!> \param pos ...
! **************************************************************************************************
   PURE SUBROUTINE create_group_dist_d0(this, ngroups, dimen, pos)
      TYPE(group_dist_d0_type), INTENT(INOUT)            :: this
      INTEGER, INTENT(IN)                                :: ngroups, dimen, pos

      INTEGER, DIMENSION(2)                              :: itmp

      itmp = get_limit(dimen, ngroups, pos)
      this%starts = itmp(1)
      this%ends = itmp(2)
      this%sizes = itmp(2) - itmp(1) + 1

   END SUBROUTINE create_group_dist_d0

! **************************************************************************************************
!> \brief ...
!> \param this ...
!> \param ngroups ...
!> \param dimen ...
! **************************************************************************************************
   PURE SUBROUTINE create_group_dist_d1_i1(this, ngroups, dimen)
      TYPE(group_dist_d1_type), INTENT(INOUT)            :: this
      INTEGER, INTENT(IN)                                :: ngroups, dimen

      INTEGER                                            :: iproc
      INTEGER, DIMENSION(2)                              :: itmp

      ALLOCATE (this%starts(0:ngroups - 1))
      this%starts = 0
      ALLOCATE (this%ends(0:ngroups - 1))
      this%ends = 0
      ALLOCATE (this%sizes(0:ngroups - 1))
      this%sizes = 0

      DO iproc = 0, ngroups - 1
         itmp = get_limit(dimen, ngroups, iproc)
         this%starts(iproc) = itmp(1)
         this%ends(iproc) = itmp(2)
         this%sizes(iproc) = itmp(2) - itmp(1) + 1
      END DO

   END SUBROUTINE create_group_dist_d1_i1

! **************************************************************************************************
!> \brief ...
!> \param this ...
!> \param ngroups ...
! **************************************************************************************************
   PURE SUBROUTINE create_group_dist_d1_0(this, ngroups)
      TYPE(group_dist_d1_type), INTENT(INOUT)            :: this
      INTEGER, INTENT(IN)                                :: ngroups

      ALLOCATE (this%starts(0:ngroups - 1))
      this%starts = 0
      ALLOCATE (this%ends(0:ngroups - 1))
      this%ends = 0
      ALLOCATE (this%sizes(0:ngroups - 1))
      this%sizes = 0

   END SUBROUTINE create_group_dist_d1_0

! **************************************************************************************************
!> \brief ...
!> \param this ...
!> \param starts ...
!> \param ends ...
!> \param sizes ...
!> \param comm ...
! **************************************************************************************************
   SUBROUTINE create_group_dist_d1_i3(this, starts, ends, sizes, comm)
      TYPE(group_dist_d1_type), INTENT(INOUT)            :: this
      INTEGER, INTENT(IN)                                :: starts, ends, sizes

      CLASS(mp_comm_type), INTENT(IN)                    :: comm

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_group_dist_d1_i3'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ALLOCATE (this%starts(0:comm%num_pe - 1))
      ALLOCATE (this%ends(0:comm%num_pe - 1))
      ALLOCATE (this%sizes(0:comm%num_pe - 1))

      CALL comm%allgather(starts, this%starts)
      CALL comm%allgather(ends, this%ends)
      CALL comm%allgather(sizes, this%sizes)

      CALL timestop(handle)

   END SUBROUTINE create_group_dist_d1_i3

! **************************************************************************************************
!> \brief ...
!> \param this ...
!> \param group_dist_ext ...
!> \param comm ...
! **************************************************************************************************
   SUBROUTINE create_group_dist_d1_gd(this, group_dist_ext, comm)
      TYPE(group_dist_d1_type), INTENT(INOUT)            :: this
      TYPE(group_dist_d0_type), INTENT(IN)               :: group_dist_ext

      CLASS(mp_comm_type), INTENT(IN)                    :: comm

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_group_dist_d1_gd'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ALLOCATE (this%starts(0:comm%num_pe - 1))
      ALLOCATE (this%ends(0:comm%num_pe - 1))
      ALLOCATE (this%sizes(0:comm%num_pe - 1))

      CALL comm%allgather(group_dist_ext%starts, this%starts)
      CALL comm%allgather(group_dist_ext%ends, this%ends)
      CALL comm%allgather(group_dist_ext%sizes, this%sizes)

      CALL timestop(handle)

   END SUBROUTINE create_group_dist_d1_gd

! **************************************************************************************************
!> \brief ...
!> \param this ...
!> \param starts ...
!> \param ends ...
!> \param sizes ...
! **************************************************************************************************
   PURE SUBROUTINE get_group_dist_d0(this, starts, ends, sizes)
      TYPE(group_dist_d0_type), INTENT(IN)               :: this
      INTEGER, INTENT(OUT), OPTIONAL                     :: starts, ends, sizes

      IF (PRESENT(starts)) starts = this%starts
      IF (PRESENT(ends)) ends = this%ends
      IF (PRESENT(sizes)) sizes = this%sizes

   END SUBROUTINE get_group_dist_d0

! **************************************************************************************************
!> \brief ...
!> \param this ...
!> \param pos ...
!> \param starts ...
!> \param ends ...
!> \param sizes ...
!> \param group_dist_ext ...
! **************************************************************************************************
   PURE SUBROUTINE get_group_dist_d1(this, pos, starts, ends, sizes, group_dist_ext)
      TYPE(group_dist_d1_type), INTENT(IN)               :: this
      INTEGER, INTENT(IN)                                :: pos
      INTEGER, INTENT(OUT), OPTIONAL                     :: starts, ends, sizes
      TYPE(group_dist_d0_type), INTENT(OUT), OPTIONAL    :: group_dist_ext

      IF (PRESENT(starts)) starts = this%starts(pos)
      IF (PRESENT(ends)) ends = this%ends(pos)
      IF (PRESENT(sizes)) sizes = this%sizes(pos)

      IF (PRESENT(group_dist_ext)) THEN
         group_dist_ext%starts = this%starts(pos)
         group_dist_ext%ends = this%ends(pos)
         group_dist_ext%sizes = this%sizes(pos)
      END IF

   END SUBROUTINE get_group_dist_d1

! **************************************************************************************************
!> \brief ...
!> \param this ...
!> \param pos ...
!> \param group_dist_ext ...
!> \param pos_ext ...
! **************************************************************************************************
   PURE SUBROUTINE get_group_dist_gd1(this, pos, group_dist_ext, pos_ext)
      TYPE(group_dist_d1_type), INTENT(IN)               :: this
      INTEGER, INTENT(IN)                                :: pos
      TYPE(group_dist_d1_type), INTENT(INOUT)            :: group_dist_ext
      INTEGER, INTENT(IN)                                :: pos_ext

      group_dist_ext%starts(pos_ext) = this%starts(pos)
      group_dist_ext%ends(pos_ext) = this%ends(pos)
      group_dist_ext%sizes(pos_ext) = this%sizes(pos)

   END SUBROUTINE get_group_dist_gd1

! **************************************************************************************************
!> \brief ...
!> \param this ...
! **************************************************************************************************
   PURE SUBROUTINE release_group_dist_d1(this)
      TYPE(group_dist_d1_type), INTENT(INOUT)            :: this

      IF (ALLOCATED(this%starts)) DEALLOCATE (this%starts)
      IF (ALLOCATED(this%ends)) DEALLOCATE (this%ends)
      IF (ALLOCATED(this%sizes)) DEALLOCATE (this%sizes)

   END SUBROUTINE release_group_dist_d1

! **************************************************************************************************
!> \brief ...
!> \param this ...
!> \return ...
! **************************************************************************************************
   ELEMENTAL FUNCTION maxsize(this) RESULT(res)
      TYPE(group_dist_d1_type), INTENT(IN)               :: this
      INTEGER                                            :: res

      res = MAXVAL(this%sizes)

   END FUNCTION maxsize

! **************************************************************************************************
!> \brief ...
!> \param this ...
!> \param pos ...
!> \return ...
! **************************************************************************************************
   ELEMENTAL FUNCTION group_dist_proc(this, pos) RESULT(proc)
      TYPE(group_dist_d1_type), INTENT(IN)               :: this
      INTEGER, INTENT(IN)                                :: pos
      INTEGER                                            :: proc

      INTEGER                                            :: p

      proc = -1
      DO p = 0, SIZE(this%sizes) - 1
         IF (pos <= this%ends(p) .AND. pos >= this%starts(p)) THEN
            proc = p
            RETURN
         END IF
      END DO

   END FUNCTION group_dist_proc

END MODULE group_dist_types
