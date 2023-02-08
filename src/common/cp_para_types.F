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
                                              mp_comm_null,&
                                              mp_comm_type

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_para_types'

   PUBLIC :: cp_para_env_type, cp_para_env_p_type, cp_para_cart_type

! **************************************************************************************************
!> \brief stores all the informations relevant to an mpi environment
!> \param mepos rank of the actual processor
!> \param group id of the actual mpi communicator (I know that it is confusing,
!>        but this is not the mpi group)
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
   TYPE cp_para_env_type
      LOGICAL :: owns_group, ionode
      INTEGER :: mepos, source, num_pe, ref_count
      TYPE(mp_comm_type) :: group = mp_comm_null
   CONTAINS
      PROCEDURE, PUBLIC, PASS(para_env) :: init => cp_para_env_init
      PROCEDURE, PRIVATE, PASS(para_env) :: cp_para_env_assign
      GENERIC, PUBLIC :: ASSIGNMENT(=) => cp_para_env_assign
      PROCEDURE, PUBLIC, PASS(para_env) :: free => cp_para_env_free
   END TYPE cp_para_env_type

! **************************************************************************************************
!> \brief represent a pointer to a para env (to build arrays)
!> \param para_env the pointer to the para_env
!> \par History
!>      07.2003 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   TYPE cp_para_env_p_type
      TYPE(cp_para_env_type), POINTER :: para_env
   END TYPE cp_para_env_p_type

! **************************************************************************************************
!> \brief represent a multidimensional parallel environment
!> \param mepos the position of the actual processor
!> \param group id of the actual group (context, communicator)
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
      LOGICAL :: owns_group
      INTEGER :: rank, ntask
      INTEGER, DIMENSION(:), POINTER :: mepos, source, num_pe
      LOGICAL, DIMENSION(:), POINTER :: periodic
      INTEGER :: ref_count
   CONTAINS
      PROCEDURE, PUBLIC, PASS(comm) :: free => cp_para_cart_free
      PROCEDURE, PUBLIC, PASS(comm) :: mp_comm_init => cp_para_cart_init
   END TYPE cp_para_cart_type

CONTAINS
! **************************************************************************************************
!> \brief ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE cp_para_env_init(para_env)
      CLASS(cp_para_env_type), INTENT(INOUT)                  :: para_env

      para_env%source = 0
      para_env%ref_count = 1
      CALL para_env%group%get_size(para_env%num_pe)
      CALL para_env%group%get_rank(para_env%mepos)
      para_env%ionode = para_env%mepos == para_env%source
   END SUBROUTINE cp_para_env_init

! **************************************************************************************************
!> \brief ...
!> \param para_env ...
!> \param comm ...
! **************************************************************************************************
   SUBROUTINE cp_para_env_assign(para_env, comm)
      CLASS(cp_para_env_type), INTENT(INOUT) :: para_env
      CLASS(mp_comm_type), INTENT(IN) :: comm

      para_env%group = comm
      para_env%owns_group = .FALSE.

      CALL para_env%init()

   END SUBROUTINE cp_para_env_assign

! **************************************************************************************************
!> \brief ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE cp_para_env_free(para_env)
      CLASS(cp_para_env_type), INTENT(INOUT) :: para_env

      IF (para_env%owns_group) CALL para_env%group%free()

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
!> \param comm ...
! **************************************************************************************************
   SUBROUTINE cp_para_cart_free(comm)
      CLASS(cp_para_cart_type), INTENT(INOUT) :: comm

      IF (comm%owns_group) CALL comm%mp_cart_type%free()

      DEALLOCATE (comm%source, comm%periodic, comm%mepos, comm%num_pe)

   END SUBROUTINE cp_para_cart_free

END MODULE cp_para_types
