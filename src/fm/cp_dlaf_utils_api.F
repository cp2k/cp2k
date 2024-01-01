!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE cp_dlaf_utils_api
   USE ISO_C_BINDING,                   ONLY: C_CHAR,&
                                              C_INT,&
                                              C_LOC,&
                                              C_NULL_CHAR,&
                                              C_PTR
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dlaf_utils_api'

   PUBLIC :: cp_dlaf_initialize, cp_dlaf_finalize
   PUBLIC :: cp_dlaf_create_grid, cp_dlaf_free_grid

CONTAINS

! **************************************************************************************************
!> \brief Initialize DLA-Future and pika runtime
!> \author Rocco Meli
!> \author Mikael Simberg
!> \author Mathieu Taillefumier
! **************************************************************************************************
   SUBROUTINE cp_dlaf_initialize()
      CHARACTER(len=*), PARAMETER :: routineN = 'cp_dlaf_initialize'
      INTEGER, PARAMETER                                 :: dlaf_argc = 1, pika_argc = 1

      CHARACTER(len=5, kind=C_CHAR), ALLOCATABLE, TARGET :: dlaf_argv(:), pika_argv(:)
      INTEGER                                            :: handle
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:)             :: dlaf_argv_ptr, pika_argv_ptr
      INTERFACE
         SUBROUTINE dlaf_init(pika_argc, pika_argv, dlaf_argc, dlaf_argv) BIND(C, name='dlaf_initialize')
            IMPORT :: C_PTR, C_INT
            TYPE(c_ptr), DIMENSION(*) :: pika_argv
            TYPE(c_ptr), DIMENSION(*) :: dlaf_argv
            INTEGER(kind=c_int), value :: pika_argc
            INTEGER(kind=c_int), value :: dlaf_argc
         END SUBROUTINE dlaf_init
      END INTERFACE

      CALL timeset(routineN, handle)

#if defined(__DLAF)
      ALLOCATE (pika_argv(pika_argc))
      pika_argv(1) = "dlaf"//C_NULL_CHAR
      ALLOCATE (dlaf_argv(dlaf_argc))
      dlaf_argv(1) = "dlaf"//C_NULL_CHAR

      ALLOCATE (pika_argv_ptr(pika_argc))
      pika_argv_ptr(1) = c_loc(pika_argv(1))
      ALLOCATE (dlaf_argv_ptr(dlaf_argc))
      dlaf_argv_ptr(1) = c_loc(dlaf_argv(1))

      CALL dlaf_init(pika_argc, pika_argv_ptr, dlaf_argc, dlaf_argv_ptr)
#else
      MARK_USED(pika_argv)
      MARK_USED(pika_argv_ptr)
      MARK_USED(pika_argc)
      MARK_USED(dlaf_argv)
      MARK_USED(dlaf_argv_ptr)
      MARK_USED(dlaf_argc)
#endif

      CALL timestop(handle)
   END SUBROUTINE cp_dlaf_initialize

! **************************************************************************************************
!> \brief Finalize DLA-Future and pika runtime
!> \author Rocco Meli
!> \author Mikael Simberg
!> \author Mathieu Taillefumier
! **************************************************************************************************
   SUBROUTINE cp_dlaf_finalize()
      CHARACTER(len=*), PARAMETER                        :: routineN = 'cp_dlaf_finalize'

      INTEGER                                            :: handle
      INTERFACE
         SUBROUTINE dlaf_finalize_aux() &
            BIND(C, name='dlaf_finalize')
         END SUBROUTINE dlaf_finalize_aux
      END INTERFACE

      CALL timeset(routineN, handle)

#if defined(__DLAF)
      CALL dlaf_finalize_aux()
#endif

      CALL timestop(handle)
   END SUBROUTINE cp_dlaf_finalize

! **************************************************************************************************
!> \brief Create DLA-Future grid from BLACS context
!> \param blacs_context ...
!> \author Rocco Meli
!> \author Mikael Simberg
!> \author Mathieu Taillefumier
! **************************************************************************************************
   SUBROUTINE cp_dlaf_create_grid(blacs_context)
      INTEGER, INTENT(IN)                                :: blacs_context

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_dlaf_create_grid'

      INTEGER                                            :: handle
      INTERFACE
         SUBROUTINE dlaf_create_grid(blacs_contxt) &
            BIND(C, name='dlaf_create_grid_from_blacs')
            IMPORT :: C_INT
            INTEGER(KIND=C_INT), VALUE :: blacs_contxt
         END SUBROUTINE
      END INTERFACE

      CALL timeset(routineN, handle)

#if defined(__DLAF)
      CALL dlaf_create_grid(blacs_context)
#else
      MARK_USED(blacs_context)
      CPABORT("CP2K compiled without the DLA-Future library.")
#endif

      CALL timestop(handle)
   END SUBROUTINE cp_dlaf_create_grid

! **************************************************************************************************
!> \brief Free DLA-Future grid corresponding to BLACS context
!> \param blacs_context ...
!> \author Rocco Meli
!> \author Mikael Simberg
!> \author Mathieu Taillefumier
! **************************************************************************************************
   SUBROUTINE cp_dlaf_free_grid(blacs_context)
      INTEGER, INTENT(IN)                                :: blacs_context

      CHARACTER(len=*), PARAMETER                        :: routineN = 'cp_dlaf_free_grid'

      INTEGER                                            :: handle
      INTERFACE
         SUBROUTINE dlaf_free_grid(blacs_contxt) &
            BIND(C, name='dlaf_free_grid')
            IMPORT :: C_INT
            INTEGER(KIND=C_INT), VALUE :: blacs_contxt
         END SUBROUTINE
      END INTERFACE

      CALL timeset(routineN, handle)

#if defined(__DLAF)
      CALL dlaf_free_grid(blacs_context)
#else
      MARK_USED(blacs_context)
      CPABORT("CP2K compiled without the DLA-Future library.")
#endif

      CALL timestop(handle)
   END SUBROUTINE cp_dlaf_free_grid

END MODULE cp_dlaf_utils_api
