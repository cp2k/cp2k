!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Wrapper module for MiMiC Communication Library (MCL) request labels
!> \par History
!>       05.2025 Created [AA]
!> \author Andrej Antalik
! **************************************************************************************************

MODULE mcl_requests

   #:include 'mcl_wrapper.fypp'

#if defined(__MIMIC)
   #:for label in mcl_implemented_labels
      USE mcl, ONLY: ${label}$
   #:endfor
#endif

   IMPLICIT NONE

   PRIVATE

   #:for label in mcl_implemented_labels
      PUBLIC :: ${label}$
   #:endfor

#if !defined(__MIMIC)
   #:set req_value = -1
   #:for label in mcl_implemented_labels
      #:set req_value = req_value - 1
      INTEGER, PARAMETER :: ${label}$ = ${req_value}$
   #:endfor
#endif

END MODULE mcl_requests

