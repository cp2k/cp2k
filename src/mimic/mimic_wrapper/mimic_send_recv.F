!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Wrapper module for the send and receive routines of a MiMiC Communication
!>        Library (MCL) client
!> \par History
!>       03.2023 Created [AA]
!> \author Andrej Antalik
! **************************************************************************************************

MODULE mimic_send_recv

   USE kinds,                           ONLY: dp
#if defined(__MIMIC)
   USE mcl,                             ONLY: mcl_receive, &
                                              mcl_send
#endif

#include "../../base/base_uses.f90"

   IMPLICIT NONE

CONTAINS

! **************************************************************************************************
!> \brief Receive an integer from the MCL server/client
!> \param buf buffer to put received data to
!> \param number of elements in the buffer
!> \param tag tag of the messeage (MiMiC command)
!> \param source id of the source of the data
! **************************************************************************************************
   SUBROUTINE recv_int(buf, count, tag, source)
      INTEGER, INTENT(OUT)                             :: buf
      INTEGER, INTENT(IN)                              :: count, tag, source

#if defined(__MIMIC)
      CALL mcl_receive(buf, count, tag, source)
#else
      CPABORT("This CP2K executable has not been linked against MiMiC.")
#endif

   END SUBROUTINE recv_int


! **************************************************************************************************
!> \brief Receive an integer 1D array from the MCL server/client
!> \param buf buffer to put received data to
!> \param number of elements in the buffer
!> \param tag tag of the messeage (MiMiC command)
!> \param source id of the source of the data
! **************************************************************************************************
   SUBROUTINE recv_int_array1d(buf, count, tag, source)
      INTEGER, DIMENSION(:), INTENT(OUT)               :: buf
      INTEGER, INTENT(IN)                              :: count, tag, source

#if defined(__MIMIC)
      CALL mcl_receive(buf, count, tag, source)
#else
      CPABORT("This CP2K executable has not been linked against MiMiC.")
#endif

   END SUBROUTINE recv_int_array1d


! **************************************************************************************************
!> \brief Receive a double precission real from the MCL server/client
!> \param buf buffer to put received data to
!> \param number of elements in the buffer
!> \param tag tag of the messeage (MiMiC command)
!> \param source id of the source of the data
! **************************************************************************************************
   SUBROUTINE recv_dp(buf, count, tag, source)
      REAL(dp), INTENT(OUT)                            :: buf
      INTEGER, INTENT(IN)                              :: count, tag, source

#if defined(__MIMIC)
      CALL mcl_receive(buf, count, tag, source)
#else
      CPABORT("This CP2K executable has not been linked against MiMiC.")
#endif

   END SUBROUTINE recv_dp


! **************************************************************************************************
!> \brief Receive a double precission real 1D array from the MCL server/client
!> \param buf buffer to put received data to
!> \param number of elements in the buffer
!> \param tag tag of the messeage (MiMiC command)
!> \param source id of the source of the data
! **************************************************************************************************
   SUBROUTINE recv_dp_array1d(buf, count, tag, source)
      REAL(dp), DIMENSION(:), INTENT(OUT)              :: buf
      INTEGER, INTENT(IN)                              :: count, tag, source

#if defined(__MIMIC)
      CALL mcl_receive(buf, count, tag, source)
#else
      CPABORT("This CP2K executable has not been linked against MiMiC.")
#endif

   END SUBROUTINE recv_dp_array1d


! **************************************************************************************************
!> \brief Send an integer to the MCL server/client
!> \param buf buffer to send data from
!> \param number of elements in the buffer
!> \param tag of the messeage (MiMiC command)
!> \param dest id of the client to send data to
! **************************************************************************************************
   SUBROUTINE send_int(buf, count, tag, dest)
      INTEGER, INTENT(IN)                              :: buf
      INTEGER, INTENT(IN)                              :: count, tag, dest

#if defined(__MIMIC)
      CALL mcl_send(buf, count, tag, dest)
#else
      CPABORT("This CP2K executable has not been linked against MiMiC.")
#endif

   END SUBROUTINE send_int


! **************************************************************************************************
!> \brief Send an integer 1D array to the MCL server/client
!> \param buf buffer to send data from
!> \param number of elements in the buffer
!> \param tag of the messeage (MiMiC command)
!> \param dest id of the client to send data to
! **************************************************************************************************
   SUBROUTINE send_int_array1d(buf, count, tag, dest)
      INTEGER, DIMENSION(:), INTENT(IN)                :: buf
      INTEGER, INTENT(IN)                              :: count, tag, dest

#if defined(__MIMIC)
      CALL mcl_send(buf, count, tag, dest)
#else
      CPABORT("This CP2K executable has not been linked against MiMiC.")
#endif

   END SUBROUTINE send_int_array1d


! **************************************************************************************************
!> \brief Send a double precission real to the MCL server/client
!> \param buf buffer to send data from
!> \param number of elements in the buffer
!> \param tag of the messeage (MiMiC command)
!> \param dest id of the client to send data to
! **************************************************************************************************
   SUBROUTINE send_dp(buf, count, tag, dest)
      REAL(dp), INTENT(IN)                             :: buf
      INTEGER, INTENT(IN)                              :: count, tag, dest

#if defined(__MIMIC)
      CALL mcl_send(buf, count, tag, dest)
#else
      CPABORT("This CP2K executable has not been linked against MiMiC.")
#endif

   END SUBROUTINE send_dp


! **************************************************************************************************
!> \brief Send a double precission real 1D array to the MCL server/client
!> \param buf buffer to send data from
!> \param number of elements in the buffer
!> \param tag of the messeage (MiMiC command)
!> \param dest id of the client to send data to
! **************************************************************************************************
   SUBROUTINE send_dp_array1d(buf, count, tag, dest)
      REAL(dp), DIMENSION(:), INTENT(IN)               :: buf
      INTEGER, INTENT(IN)                              :: count, tag, dest

#if defined(__MIMIC)
      CALL mcl_send(buf, count, tag, dest)
#else
      CPABORT("This CP2K executable has not been linked against MiMiC.")
#endif

   END SUBROUTINE send_dp_array1d

! **************************************************************************************************
!> \brief Send a logical to the MCL server/client
!> \param buf buffer to send data from
!> \param number of elements in the buffer
!> \param tag of the messeage (MiMiC command)
!> \param dest id of the client to send data to
! **************************************************************************************************
   SUBROUTINE send_logical(val, count, tag, dest)
      LOGICAL, INTENT(IN)                              :: val
      INTEGER, INTENT(IN)                              :: count, tag, dest

      INTEGER                                          :: buf

#if defined(__MIMIC)
      buf = merge(1, 0, val)
      CALL mcl_send(buf, count, tag, dest)
#else
      CPABORT("This CP2K executable has not been linked against MiMiC.")
#endif

   END SUBROUTINE send_logical


END MODULE mimic_send_recv

