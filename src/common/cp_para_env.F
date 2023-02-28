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

   PUBLIC :: cp_para_env_release, cp_para_env_create
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
      para_env%mp_comm_type = group
      CALL para_env%init()
   END SUBROUTINE cp_para_env_create

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
         CPASSERT(para_env%is_valid())
         CALL para_env%free()
         IF (.NOT. para_env%is_valid()) DEALLOCATE (para_env)
      END IF
      NULLIFY (para_env)
   END SUBROUTINE cp_para_env_release

! **************************************************************************************************
!> \brief creates a cart (multidimensional parallel environment)
!> \param cart the cart environment to create
!> \param group the mpi communicator
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE cp_cart_create(cart, group)
      TYPE(cp_para_cart_type), POINTER                   :: cart

      CLASS(mp_comm_type), INTENT(in)                     :: group

      CPASSERT(.NOT. ASSOCIATED(cart))
      ALLOCATE (cart)
      cart%mp_comm_type = group
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
         CPASSERT(cart%is_valid())
         CALL cart%free()
         IF (.NOT. cart%is_valid()) DEALLOCATE (cart)
      END IF
      NULLIFY (cart)
   END SUBROUTINE cp_cart_release

END MODULE cp_para_env
