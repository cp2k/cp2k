!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE cp_dlaf_utils
   USE ISO_C_BINDING,                   ONLY: C_INT
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dlaf_utils'

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

      INTEGER                                            :: handle
      INTERFACE
         SUBROUTINE dlaf_init_aux() &
            BIND(C, name='dlaf_init')
         END SUBROUTINE dlaf_init_aux
      END INTERFACE

      CALL timeset(routineN, handle)

#if defined(__DLAF)
      CALL dlaf_init_aux()
#else
      CPABORT("CP2K compiled without the DLA-Future library.")
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
#else
      CPABORT("CP2K compiled without the DLA-Future library.")
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

END MODULE cp_dlaf_utils
