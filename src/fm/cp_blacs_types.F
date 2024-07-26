!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief wrappers for the actual blacs calls.
!>      all functionality needed in the code should actually be provide by cp_blacs_env
!>      these functions should be private members of that module
!> \note
!>      http://www.netlib.org/blacs/BLACS/QRef.html
!> \par History
!>      12.2003 created [Joost]
!> \author Joost VandeVondele
! **************************************************************************************************
MODULE cp_blacs_types

#if defined(__DLAF)
   USE cp_dlaf_utils_api, ONLY: cp_dlaf_create_grid, &
                                cp_dlaf_free_grid
#endif
   USE kinds, ONLY: dp
   USE message_passing, ONLY: mp_comm_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_blacs_types'
   PRIVATE

   PUBLIC :: cp_blacs_type

   TYPE cp_blacs_type
      PRIVATE
#if defined(__parallel)
      INTEGER :: context_handle = -1
#endif
      INTEGER, DIMENSION(2), PUBLIC :: mepos = -1, num_pe = -1
   CONTAINS
      PROCEDURE, PUBLIC, PASS(this), NON_OVERRIDABLE :: gridinit => cp_blacs_gridinit
      PROCEDURE, PUBLIC, PASS(this), NON_OVERRIDABLE :: gridexit => cp_blacs_gridexit
      PROCEDURE, PRIVATE, PASS(this), NON_OVERRIDABLE :: gridinfo => cp_blacs_gridinfo
      PROCEDURE, PUBLIC, PASS(this), NON_OVERRIDABLE :: set => cp_blacs_set
      PROCEDURE, PUBLIC, PASS(this), NON_OVERRIDABLE :: zgebs2d => cp_blacs_zgebs2d
      PROCEDURE, PUBLIC, PASS(this), NON_OVERRIDABLE :: dgebs2d => cp_blacs_dgebs2d
      PROCEDURE, PUBLIC, PASS(this), NON_OVERRIDABLE :: zgebr2d => cp_blacs_zgebr2d
      PROCEDURE, PUBLIC, PASS(this), NON_OVERRIDABLE :: dgebr2d => cp_blacs_dgebr2d
      PROCEDURE, PUBLIC, PASS(this), NON_OVERRIDABLE :: get_handle => cp_blacs_get_handle

      PROCEDURE, PRIVATE, PASS(this), NON_OVERRIDABLE :: cp_context_is_equal
      GENERIC, PUBLIC :: OPERATOR(==) => cp_context_is_equal

      PROCEDURE, PRIVATE, PASS(this), NON_OVERRIDABLE :: cp_context_is_not_equal
      GENERIC, PUBLIC :: OPERATOR(/=) => cp_context_is_not_equal

      PROCEDURE, PUBLIC, PASS(this), NON_OVERRIDABLE :: interconnect => cp_blacs_interconnect
   END TYPE

!***
CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param this ...
!> \param comm ...
!> \param order ...
!> \param nprow ...
!> \param npcol ...
! **************************************************************************************************
   SUBROUTINE cp_blacs_gridinit(this, comm, order, nprow, npcol)
      CLASS(cp_blacs_type), INTENT(OUT) :: this
      CLASS(mp_comm_type), INTENT(IN) :: comm
      CHARACTER(len=1), INTENT(IN):: order
      INTEGER, INTENT(IN)    :: nprow, npcol
#if defined(__parallel)
      INTEGER :: context_handle
      context_handle = comm%get_handle()
      CALL blacs_gridinit(context_handle, order, nprow, npcol)
      this%context_handle = context_handle
#if defined(__DLAF)
      CALL cp_dlaf_create_grid(context_handle)
#endif
#else
      MARK_USED(this)
      MARK_USED(comm)
      MARK_USED(order)
      MARK_USED(nprow)
      MARK_USED(npcol)
#endif
      CALL this%gridinfo()
   END SUBROUTINE cp_blacs_gridinit

! **************************************************************************************************
!> \brief ...
!> \param this ...
! **************************************************************************************************
   SUBROUTINE cp_blacs_gridexit(this)
      CLASS(cp_blacs_type), INTENT(IN) :: this
#if defined(__parallel)
      CALL blacs_gridexit(this%context_handle)
#if defined(__DLAF)
      CALL cp_dlaf_free_grid(this%context_handle)
#endif
#else
      MARK_USED(this)
#endif
   END SUBROUTINE cp_blacs_gridexit

! **************************************************************************************************
!> \brief ...
!> \param this ...
! **************************************************************************************************
   SUBROUTINE cp_blacs_gridinfo(this)
      CLASS(cp_blacs_type), INTENT(INOUT)  :: this
#if defined(__parallel)
      CALL blacs_gridinfo(this%context_handle, this%num_pe(1), this%num_pe(2), this%mepos(1), this%mepos(2))
#else
      MARK_USED(this)
      this%num_pe = 1
      this%mepos = 0
#endif
   END SUBROUTINE cp_blacs_gridinfo

! **************************************************************************************************
!> \brief ...
!> \param this ...
!> \param what :
!>     WHAT = 0 : Handle indicating default system context;  ! DO NOT USE (i.e. use para_env)
!>     WHAT = 1 : The BLACS message ID range;
!>     WHAT = 2 : The BLACS debug level the library was compiled with;
!>     WHAT = 10: Handle indicating the system context used to define the BLACS context whose handle is ICONTXT;
!>     WHAT = 11: Number of rings multiring topology is presently using;
!>     WHAT = 12: Number of branches general tree topology is presently using.
!>     WHAT = 15: If non-zero, makes topology choice for repeatable collectives
!> \param val ...
! **************************************************************************************************
   SUBROUTINE cp_blacs_set(this, what, val)
      CLASS(cp_blacs_type), INTENT(IN) :: this
      INTEGER, INTENT(IN)  :: what, val
#if defined(__parallel)
      CALL blacs_set(this%context_handle, what, val)
#else
      MARK_USED(this)
      MARK_USED(what)
      MARK_USED(val)
#endif
   END SUBROUTINE cp_blacs_set

! **************************************************************************************************
!> \brief ...
!> \param this ...
!> \param SCOPE ...
!> \param TOP ...
!> \param M ...
!> \param N ...
!> \param A ...
!> \param LDA ...
! **************************************************************************************************
   SUBROUTINE cp_blacs_zgebs2d(this, SCOPE, TOP, M, N, A, LDA)
      CLASS(cp_blacs_type), INTENT(IN)     :: this
      CHARACTER(len=1), INTENT(IN) :: SCOPE, TOP
      INTEGER, INTENT(IN)     :: M, N, LDA
      COMPLEX(KIND=dp)            :: A
#if defined(__parallel)
      CALL zgebs2d(this%context_handle, SCOPE, TOP, M, N, A, LDA)
#else
      MARK_USED(this)
      MARK_USED(SCOPE)
      MARK_USED(TOP)
      MARK_USED(M)
      MARK_USED(N)
      MARK_USED(A)
      MARK_USED(LDA)
#endif
   END SUBROUTINE
! **************************************************************************************************
!> \brief ...
!> \param this ...
!> \param SCOPE ...
!> \param TOP ...
!> \param M ...
!> \param N ...
!> \param A ...
!> \param LDA ...
!> \param RSRC ...
!> \param CSRC ...
! **************************************************************************************************
   SUBROUTINE cp_blacs_zgebr2d(this, SCOPE, TOP, M, N, A, LDA, RSRC, CSRC)
      CLASS(cp_blacs_type), INTENT(IN)     :: this
      CHARACTER(len=1), INTENT(IN) :: SCOPE, TOP
      INTEGER, INTENT(IN)     :: M, N, LDA
      INTEGER, INTENT(IN)     :: RSRC, CSRC
      COMPLEX(KIND=dp)            :: A
#if defined(__parallel)
      CALL zgebr2d(this%context_handle, SCOPE, TOP, M, N, A, LDA, RSRC, CSRC)
#else
      MARK_USED(this)
      MARK_USED(SCOPE)
      MARK_USED(TOP)
      MARK_USED(M)
      MARK_USED(N)
      MARK_USED(A)
      MARK_USED(LDA)
      MARK_USED(RSRC)
      MARK_USED(CSRC)
#endif
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param this ...
!> \param SCOPE ...
!> \param TOP ...
!> \param M ...
!> \param N ...
!> \param A ...
!> \param LDA ...
! **************************************************************************************************
   SUBROUTINE cp_blacs_dgebs2d(this, SCOPE, TOP, M, N, A, LDA)
      CLASS(cp_blacs_type), INTENT(IN)     :: this
      CHARACTER(len=1), INTENT(IN) :: SCOPE, TOP
      INTEGER, INTENT(IN)     :: M, N, LDA
      REAL(KIND=dp)               :: A
#if defined(__parallel)
      CALL dgebs2d(this%context_handle, SCOPE, TOP, M, N, A, LDA)
#else
      MARK_USED(this)
      MARK_USED(SCOPE)
      MARK_USED(TOP)
      MARK_USED(M)
      MARK_USED(N)
      MARK_USED(A)
      MARK_USED(LDA)
#endif
   END SUBROUTINE
! **************************************************************************************************
!> \brief ...
!> \param this ...
!> \param SCOPE ...
!> \param TOP ...
!> \param M ...
!> \param N ...
!> \param A ...
!> \param LDA ...
!> \param RSRC ...
!> \param CSRC ...
! **************************************************************************************************
   SUBROUTINE cp_blacs_dgebr2d(this, SCOPE, TOP, M, N, A, LDA, RSRC, CSRC)
      CLASS(cp_blacs_type), INTENT(IN)     :: this
      CHARACTER(len=1), INTENT(IN) :: SCOPE, TOP
      INTEGER, INTENT(IN)     :: M, N, LDA
      INTEGER, INTENT(IN)     :: RSRC, CSRC
      REAL(KIND=dp)               :: A
#if defined(__parallel)
      CALL dgebr2d(this%context_handle, SCOPE, TOP, M, N, A, LDA, RSRC, CSRC)
#else
      MARK_USED(this)
      MARK_USED(SCOPE)
      MARK_USED(TOP)
      MARK_USED(M)
      MARK_USED(N)
      MARK_USED(A)
      MARK_USED(LDA)
      MARK_USED(RSRC)
      MARK_USED(CSRC)
#endif
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param this ...
!> \return ...
! **************************************************************************************************
   ELEMENTAL INTEGER FUNCTION cp_blacs_get_handle(this)
      CLASS(cp_blacs_type), INTENT(IN) :: this
#if defined(__parallel)
      cp_blacs_get_handle = this%context_handle
#else
      MARK_USED(this)
      cp_blacs_get_handle = -1
#endif
   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param this ...
!> \param other ...
!> \return ...
! **************************************************************************************************
   ELEMENTAL LOGICAL FUNCTION cp_context_is_equal(this, other)
      CLASS(cp_blacs_type), INTENT(IN) :: this, other
#if defined(__parallel)
      cp_context_is_equal = (this%context_handle == other%context_handle)
#else
      MARK_USED(this)
      MARK_USED(other)
      cp_context_is_equal = .TRUE.
#endif
   END FUNCTION cp_context_is_equal

! **************************************************************************************************
!> \brief ...
!> \param this ...
!> \param other ...
!> \return ...
! **************************************************************************************************
   ELEMENTAL LOGICAL FUNCTION cp_context_is_not_equal(this, other)
      CLASS(cp_blacs_type), INTENT(IN) :: this, other
#if defined(__parallel)
      cp_context_is_not_equal = (this%context_handle /= other%context_handle)
#else
      MARK_USED(this)
      MARK_USED(other)
      cp_context_is_not_equal = .FALSE.
#endif
   END FUNCTION cp_context_is_not_equal

! **************************************************************************************************
!> \brief ...
!> \param this ...
!> \param comm_super ...
!> \return ...
! **************************************************************************************************
   TYPE(mp_comm_type) FUNCTION cp_blacs_interconnect(this, comm_super)
      CLASS(cp_blacs_type), INTENT(IN) :: this
      CLASS(mp_comm_type), INTENT(IN) :: comm_super

      INTEGER :: blacs_coord

! We enumerate the processes within the process grid in a linear fashion
      blacs_coord = this%mepos(1)*this%num_pe(2) + this%mepos(2)

      CALL cp_blacs_interconnect%from_split(comm_super, blacs_coord)

   END FUNCTION cp_blacs_interconnect

END MODULE cp_blacs_types
