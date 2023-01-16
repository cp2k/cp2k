!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief type to store parallelization informations (at the moment assumes 1d
!>      position and uses mpi)
!> \par History
!>      07.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
MODULE cp_para_env
   USE cp_para_types,                   ONLY: cp_para_cart_type,&
                                              cp_para_env_type
   USE message_passing,                 ONLY: mp_comm_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_para_env'

   PUBLIC :: cp_para_env_retain, cp_para_env_release, cp_para_env_create, cp_para_env_split
   PUBLIC :: cp_cart_create, cp_cart_release
!***
CONTAINS

! **************************************************************************************************
!> \brief creates a new para environment
!> \param para_env the new parallel environment
!> \param group the id of the actual mpi_group
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE cp_para_env_create(para_env, group)
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CLASS(mp_comm_type), INTENT(in)                     :: group

      CPASSERT(.NOT. ASSOCIATED(para_env))
      ALLOCATE (para_env)
      para_env%group = group
      para_env%owns_group = .TRUE.
      CALL para_env%init()
   END SUBROUTINE cp_para_env_create

! **************************************************************************************************
!> \brief retains the para object (to be called when you want to keep a
!>      shared copy of this object)
!> \param para_env the new group
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE cp_para_env_retain(para_env)
      TYPE(cp_para_env_type), INTENT(IN), POINTER        :: para_env

      CPASSERT(ASSOCIATED(para_env))
      CPASSERT(para_env%ref_count > 0)
      para_env%ref_count = para_env%ref_count + 1
   END SUBROUTINE cp_para_env_retain

! **************************************************************************************************
!> \brief splits the communicator of the template para_env according to the provided color
!> \param para_env para_env with new sub communicator
!> \param para_env_template para_env to be split
!> \param color all processes with same color will belong to the same sub communicator
! **************************************************************************************************
   SUBROUTINE cp_para_env_split(para_env, para_env_template, color)
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(cp_para_env_type), INTENT(IN)                 :: para_env_template
      INTEGER, INTENT(IN)                                :: color

      TYPE(mp_comm_type)                                 :: comm

      CPASSERT(para_env_template%ref_count > 0)
      CALL comm%from_split(para_env_template%group, color)
      NULLIFY (para_env)
      CALL cp_para_env_create(para_env, comm)

   END SUBROUTINE cp_para_env_split

! **************************************************************************************************
!> \brief releases the para object (to be called when you don't want anymore
!>      the shared copy of this object)
!> \param para_env the new group
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      to avoid circular dependencies cp_log_handling has a private copy
!>      of this method (see cp_log_handling:my_cp_para_env_release)!
! **************************************************************************************************
   SUBROUTINE cp_para_env_release(para_env)
      TYPE(cp_para_env_type), POINTER                    :: para_env

      IF (ASSOCIATED(para_env)) THEN
         CPASSERT(para_env%ref_count > 0)
         para_env%ref_count = para_env%ref_count - 1
         IF (para_env%ref_count < 1) THEN
            CALL para_env%free()
            DEALLOCATE (para_env)
         END IF
      END IF
      NULLIFY (para_env)
   END SUBROUTINE cp_para_env_release

! **************************************************************************************************
!> \brief creates a cart (multidimensional parallel environment)
!> \param cart the cart environment to create
!> \param group the mpi communicator
!> \param ndims the number of dimensions of the cart
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE cp_cart_create(cart, group, ndims)
      TYPE(cp_para_cart_type), POINTER                   :: cart

      CLASS(mp_comm_type), INTENT(in)                     :: group
      INTEGER, INTENT(in)                                :: ndims

      CPASSERT(.NOT. ASSOCIATED(cart))
      ALLOCATE (cart)
      cart%ndims = ndims
      cart%group = group

      CALL cart%init()

   END SUBROUTINE cp_cart_create

! **************************************************************************************************
!> \brief releases the given cart
!> \param cart the cart to release
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE cp_cart_release(cart)
      TYPE(cp_para_cart_type), POINTER                   :: cart

      IF (ASSOCIATED(cart)) THEN
         CPASSERT(cart%ref_count > 0)
         cart%ref_count = cart%ref_count - 1
         IF (cart%ref_count == 0) THEN
            CALL cart%free()
            DEALLOCATE (cart)
         END IF
      END IF
      NULLIFY (cart)
   END SUBROUTINE cp_cart_release

END MODULE cp_para_env
