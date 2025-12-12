!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Wrapper module for MiMiC Communication Library (MCL) routines
!> \par History
!>       05.2025 Created [AA]
!> \author Andrej Antalik
! **************************************************************************************************

MODULE mcl_api

   USE kinds, ONLY: int_4, int_8, sp, dp
#if defined(__MIMIC)
   USE mcl, ONLY: mcl_finalize, mcl_get_program_id, mcl_receive, mcl_send
#endif

#include "../../base/base_uses.f90"

   IMPLICIT NONE

   #:include 'mcl_wrapper.fypp'

   PRIVATE

   PUBLIC :: mcl_finalize, mcl_get_program_id, mcl_receive, mcl_send

#if !defined(__MIMIC)

   #:for direction in mcl_direction
      INTERFACE mcl_${direction}$
         #:for kind_l in mcl_kind_labels
            #:for dim_l in mcl_dimension_labels
               MODULE PROCEDURE mcl_${direction}$_${kind_l}$_${dim_l}$
            #:endfor
         #:endfor
      END INTERFACE
   #:endfor

CONTAINS

   SUBROUTINE mcl_finalize()
      CPABORT("This CP2K executable has not been linked against MiMiC.")
   END SUBROUTINE mcl_finalize

   SUBROUTINE mcl_get_program_id(id)
      INTEGER :: id
      CPABORT("This CP2K executable has not been linked against MiMiC.")
   END SUBROUTINE mcl_get_program_id

   #:for direction in mcl_direction
      #:for kind_v, kind_l in zip(mcl_kinds, mcl_kind_labels)
         #:for dim_v, dim_l in zip(mcl_dimensions, mcl_dimension_labels)
            SUBROUTINE mcl_${direction}$_${kind_l}$_${dim_l}$ (buffer, length, tag, remote_id, err)
               ${kind_v}$${dim_v}$, TARGET :: buffer
               INTEGER :: length
               INTEGER :: tag
               INTEGER :: remote_id
               INTEGER, OPTIONAL, INTENT(OUT) :: err
               CPABORT("This CP2K executable has not been linked against MiMiC.")
            END SUBROUTINE mcl_${direction}$_${kind_l}$_${dim_l}$
         #:endfor
      #:endfor
   #:endfor

#endif

END MODULE mcl_api

