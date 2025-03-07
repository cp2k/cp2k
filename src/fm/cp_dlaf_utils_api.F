!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE cp_dlaf_utils_api

#include "../base/base_uses.f90"

#if defined(__DLAF)
   USE dlaf_fortran, ONLY: dlaf_create_grid_from_blacs, &
                           dlaf_finalize, &
                           dlaf_free_grid, &
                           dlaf_initialize
#endif

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

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)
#if defined(__DLAF)
      CALL dlaf_initialize()
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

      CALL timeset(routineN, handle)
#if defined(__DLAF)
      CALL dlaf_finalize()
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

      CALL timeset(routineN, handle)
#if defined(__DLAF)
      CALL dlaf_create_grid_from_blacs(blacs_context)
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
