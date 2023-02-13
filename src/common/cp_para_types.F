!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief type to store parallelization informations (at the moment assumes 1d
!>      position and uses mpi).
!>      type def in this separate module to avoid circular ref.
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
MODULE cp_para_types
   USE message_passing,                 ONLY: mp_cart_type,&
                                              mp_comm_type

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_para_types'

   PUBLIC :: cp_para_env_type, cp_para_env_p_type, cp_para_cart_type

! **************************************************************************************************
!> \brief stores all the informations relevant to an mpi environment
!> \param mepos rank of the actual processor
!> \param num_pe number of processors in the communicator
!> \param source rank of a special processor (for example the one for i-o,
!>        or the master
!> \param owns_group if it owns the group (and thus should free it when
!>        this object is deallocated)
!> \param ref_count the reference count, when it is zero this object gets
!>        deallocated
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   TYPE, EXTENDS(mp_comm_type) :: cp_para_env_type
      LOGICAL :: owns_group = .TRUE.
      INTEGER :: mepos = -1, source = -1, num_pe = -1, ref_count = -1
   CONTAINS
      PROCEDURE, PRIVATE, PASS(comm), NON_OVERRIDABLE :: mp_comm_init => cp_para_env_init
      PROCEDURE, PRIVATE, PASS(comm_new), NON_OVERRIDABLE :: mp_comm_assign => cp_para_env_assign
      PROCEDURE, PRIVATE, PASS(comm), NON_OVERRIDABLE :: mp_comm_free => cp_para_env_free
      PROCEDURE, PUBLIC, PASS(para_env), NON_OVERRIDABLE :: retain => cp_para_env_retain
      PROCEDURE, PUBLIC, PASS(para_env), NON_OVERRIDABLE :: is_source => cp_para_env_is_source
   END TYPE cp_para_env_type

! **************************************************************************************************
!> \brief represent a pointer to a para env (to build arrays)
!> \param para_env the pointer to the para_env
!> \par History
!>      07.2003 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   TYPE cp_para_env_p_type
      TYPE(cp_para_env_type), POINTER :: para_env => NULL()
   END TYPE cp_para_env_p_type

! **************************************************************************************************
!> \brief represent a multidimensional parallel environment
!> \param mepos the position of the actual processor
!> \param num_pe number of processors in the group in each dimension
!> \param source id of a special processor (for example the one for i-o,
!>        or the master
!> \param owns_group if it owns the group (and thus should free it when
!>        this object is deallocated)
!> \param ref_count the reference count, when it is zero this object gets
!>        deallocated
!> \note
!>      not yet implemented for mpi
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   TYPE, EXTENDS(mp_cart_type) :: cp_para_cart_type
      LOGICAL :: owns_group = .TRUE.
      INTEGER :: rank = -1, ntask = -1
      INTEGER, DIMENSION(:), POINTER :: mepos => NULL(), source => NULL(), num_pe => NULL()
      LOGICAL, DIMENSION(:), POINTER :: periodic => NULL()
      INTEGER :: ref_count = -1
   CONTAINS
      PROCEDURE, PRIVATE, PASS(comm), NON_OVERRIDABLE :: mp_comm_free => cp_para_cart_free
      PROCEDURE, PRIVATE, PASS(comm), NON_OVERRIDABLE :: mp_comm_init => cp_para_cart_init
      PROCEDURE, PUBLIC, PASS(cart), NON_OVERRIDABLE :: retain => cp_para_cart_retain
   END TYPE cp_para_cart_type

CONTAINS
! **************************************************************************************************
!> \brief ...
!> \param comm ...
! **************************************************************************************************
   ELEMENTAL IMPURE SUBROUTINE cp_para_env_init(comm)
      CLASS(cp_para_env_type), INTENT(INOUT)                  :: comm

      comm%source = 0
      comm%ref_count = 1
      CALL comm%mp_comm_type%get_size(comm%num_pe)
      CALL comm%mp_comm_type%get_rank(comm%mepos)
   END SUBROUTINE cp_para_env_init

! **************************************************************************************************
!> \brief ...
!> \param comm_new ...
!> \param comm_old ...
! **************************************************************************************************
   ELEMENTAL IMPURE SUBROUTINE cp_para_env_assign(comm_new, comm_old)
      CLASS(cp_para_env_type), INTENT(OUT) :: comm_new
      CLASS(mp_comm_type), INTENT(IN) :: comm_old

      comm_new%mp_comm_type = comm_old
      comm_new%owns_group = .FALSE.

      CALL comm_new%init()

   END SUBROUTINE cp_para_env_assign

! **************************************************************************************************
!> \brief ...
!> \param para_env ...
!> \return ...
! **************************************************************************************************
   ELEMENTAL LOGICAL FUNCTION cp_para_env_is_source(para_env)
      CLASS(cp_para_env_type), INTENT(IN) :: para_env

      cp_para_env_is_source = para_env%source == para_env%mepos

   END FUNCTION cp_para_env_is_source

! **************************************************************************************************
!> \brief ...
!> \param para_env ...
! **************************************************************************************************
   ELEMENTAL SUBROUTINE cp_para_env_retain(para_env)
      CLASS(cp_para_env_type), INTENT(INOUT) :: para_env

      para_env%ref_count = para_env%ref_count + 1

   END SUBROUTINE cp_para_env_retain

! **************************************************************************************************
!> \brief ...
!> \param comm ...
! **************************************************************************************************
   SUBROUTINE cp_para_env_free(comm)
      CLASS(cp_para_env_type), INTENT(INOUT) :: comm

      IF (comm%owns_group) CALL comm%mp_comm_type%free()

   END SUBROUTINE cp_para_env_free

! **************************************************************************************************
!> \brief intializes a cart (multidimensional parallel environment)
!> \param comm ...
!> \author fawzi
! **************************************************************************************************
   ELEMENTAL IMPURE SUBROUTINE cp_para_cart_init(comm)
      CLASS(cp_para_cart_type), INTENT(INOUT)             :: comm

      INTEGER :: ndims

      comm%owns_group = .TRUE.
      ndims = comm%get_ndims()

      ALLOCATE (comm%source(ndims), comm%periodic(ndims), comm%mepos(ndims), &
                comm%num_pe(ndims))

      comm%source = 0
      comm%mepos = 0
      comm%periodic = .FALSE.
      comm%ref_count = 1
      comm%ntask = 1
      CALL comm%get_info_cart(comm%num_pe, comm%mepos, &
                              comm%periodic)
      CALL comm%get_size(comm%ntask)
      CALL comm%get_rank(comm%rank)
   END SUBROUTINE cp_para_cart_init

! **************************************************************************************************
!> \brief ...
!> \param cart ...
! **************************************************************************************************
   ELEMENTAL SUBROUTINE cp_para_cart_retain(cart)
      CLASS(cp_para_cart_type), INTENT(INOUT) :: cart

      cart%ref_count = cart%ref_count + 1

   END SUBROUTINE cp_para_cart_retain

! **************************************************************************************************
!> \brief ...
!> \param comm ...
! **************************************************************************************************
   SUBROUTINE cp_para_cart_free(comm)
      CLASS(cp_para_cart_type), INTENT(INOUT) :: comm

      IF (comm%owns_group) CALL comm%mp_cart_type%free()

      DEALLOCATE (comm%source, comm%periodic, comm%mepos, comm%num_pe)

   END SUBROUTINE cp_para_cart_free

END MODULE cp_para_types
