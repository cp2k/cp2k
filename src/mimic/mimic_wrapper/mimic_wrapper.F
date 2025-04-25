!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Wrapper module for MiMiC and MiMiC Communication Library (MCL) routines
!> \par History
!>       03.2023 Created [AA]
!> \author Andrej Antalik
! **************************************************************************************************

MODULE mimic_wrapper

   USE mimic_send_recv

#if defined(__MIMIC)
   USE mcl,                             ONLY: mcl_finalize
#endif

#include "../../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: mimic_finalize, mimic_receive, mimic_send

   INTERFACE mimic_receive
      MODULE PROCEDURE recv_int, recv_int_array1d, &
                       recv_dp,  recv_dp_array1d 
   END INTERFACE

   INTERFACE mimic_send
      MODULE PROCEDURE send_int, send_int_array1d, &
                       send_dp,  send_dp_array1d, &
                       send_logical 
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief Destroys the endpoint
! **************************************************************************************************
   SUBROUTINE mimic_finalize()

#if defined(__MIMIC)
      CALL mcl_finalize()
#else
      CPABORT("This CP2K executable has not been linked against MiMiC.")
#endif

   END SUBROUTINE mimic_finalize

END MODULE mimic_wrapper

