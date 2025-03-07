!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Implements UNIX and INET sockets
!> \par History
!>      08.2023 moved here and expanded for AS module by S. Battaglia
!>      03.2012 created by MC in ipi_driver.F
!> \author M. Ceriotti
! **************************************************************************************************

MODULE sockets_interface

   USE ISO_C_BINDING,                   ONLY: C_CHAR,&
                                              C_DOUBLE,&
                                              C_INT,&
                                              C_LOC,&
                                              C_PTR
   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'sockets_interface'

#ifndef __NO_SOCKETS
   PUBLIC :: writebuffer, readbuffer, open_connect_socket, &
             uwait, open_bind_socket, listen_socket, &
             accept_socket, close_socket, remove_socket_file

   INTERFACE writebuffer
      MODULE PROCEDURE writebuffer_s, &
         writebuffer_d, writebuffer_dv, &
         writebuffer_i

   END INTERFACE

   INTERFACE readbuffer
      MODULE PROCEDURE readbuffer_s, &
         readbuffer_dv, readbuffer_d, &
         readbuffer_i

   END INTERFACE

   INTERFACE
      SUBROUTINE uwait(sec) BIND(C, NAME="uwait")
         USE ISO_C_BINDING, ONLY: C_DOUBLE
      REAL(C_DOUBLE)                                     :: sec

      END SUBROUTINE
   END INTERFACE

   INTERFACE
      SUBROUTINE open_connect_socket(psockfd, inet, port, host) BIND(C)
         IMPORT
         INTEGER(KIND=C_INT)                      :: psockfd, inet, port
         CHARACTER(KIND=C_CHAR), DIMENSION(*)     :: host

      END SUBROUTINE open_connect_socket

      SUBROUTINE open_bind_socket(psockfd, inet, port, host) BIND(C)
         IMPORT
         INTEGER(KIND=C_INT)                      :: psockfd, inet, port
         CHARACTER(KIND=C_CHAR), DIMENSION(*)     :: host

      END SUBROUTINE open_bind_socket

      SUBROUTINE listen_socket(psockfd, backlog) BIND(C)
         IMPORT
         INTEGER(KIND=C_INT)                      :: psockfd, backlog

      END SUBROUTINE listen_socket

      SUBROUTINE accept_socket(psockfd, pclientfd) BIND(C)
         IMPORT
         INTEGER(KIND=C_INT)                      :: psockfd, pclientfd

      END SUBROUTINE accept_socket

      SUBROUTINE close_socket(psockfd) BIND(C)
         IMPORT
         INTEGER(KIND=C_INT)                      :: psockfd

      END SUBROUTINE close_socket

      SUBROUTINE remove_socket_file(host) BIND(C)
         IMPORT
         CHARACTER(KIND=C_CHAR), DIMENSION(*)     :: host

      END SUBROUTINE remove_socket_file

      SUBROUTINE writebuffer_csocket(psockfd, pdata, plen) BIND(C, name="writebuffer")
         IMPORT
         INTEGER(KIND=C_INT)                      :: psockfd
         TYPE(C_PTR), VALUE                       :: pdata
         INTEGER(KIND=C_INT)                      :: plen

      END SUBROUTINE writebuffer_csocket

      SUBROUTINE readbuffer_csocket(psockfd, pdata, plen) BIND(C, name="readbuffer")
         IMPORT
         INTEGER(KIND=C_INT)                      :: psockfd
         TYPE(C_PTR), VALUE                       :: pdata
         INTEGER(KIND=C_INT)                      :: plen

      END SUBROUTINE readbuffer_csocket
   END INTERFACE
#endif

CONTAINS

#ifndef __NO_SOCKETS
! **************************************************************************************************
!> \brief ...
!> \param psockfd ...
!> \param fdata ...
! **************************************************************************************************
   SUBROUTINE writebuffer_d(psockfd, fdata)
      INTEGER, INTENT(IN)                                :: psockfd
      REAL(KIND=dp), INTENT(IN)                          :: fdata

      CHARACTER(len=*), PARAMETER                        :: routineN = 'writebuffer_d'

      INTEGER                                            :: handle
      REAL(KIND=C_DOUBLE), TARGET                        :: cdata

      CALL timeset(routineN, handle)

      cdata = fdata
      CALL writebuffer_csocket(psockfd, c_loc(cdata), 8)

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param psockfd ...
!> \param fdata ...
! **************************************************************************************************
   SUBROUTINE writebuffer_i(psockfd, fdata)
      INTEGER, INTENT(IN)                                :: psockfd, fdata

      CHARACTER(len=*), PARAMETER                        :: routineN = 'writebuffer_i'

      INTEGER                                            :: handle
      INTEGER(KIND=C_INT), TARGET                        :: cdata

      CALL timeset(routineN, handle)

      cdata = fdata
      CALL writebuffer_csocket(psockfd, c_loc(cdata), 4)

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param psockfd ...
!> \param fstring ...
!> \param plen ...
! **************************************************************************************************
   SUBROUTINE writebuffer_s(psockfd, fstring, plen)
      INTEGER, INTENT(IN)                                :: psockfd
      CHARACTER(LEN=*), INTENT(IN)                       :: fstring
      INTEGER, INTENT(IN)                                :: plen

      CHARACTER(len=*), PARAMETER                        :: routineN = 'writebuffer_s'

      INTEGER                                            :: handle, i
      CHARACTER(LEN=1, KIND=C_CHAR), TARGET              :: cstring(plen)

      CALL timeset(routineN, handle)

      DO i = 1, plen
         cstring(i) = fstring(i:i)
      END DO
      CALL writebuffer_csocket(psockfd, c_loc(cstring(1)), plen)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param psockfd ...
!> \param fdata ...
!> \param plen ...
! **************************************************************************************************
   SUBROUTINE writebuffer_dv(psockfd, fdata, plen)
      INTEGER, INTENT(IN)                                :: psockfd, plen
      REAL(KIND=dp), INTENT(IN), TARGET                  :: fdata(plen)

      CHARACTER(len=*), PARAMETER                        :: routineN = 'writebuffer_dv'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL writebuffer_csocket(psockfd, c_loc(fdata(1)), 8*plen)

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param psockfd ...
!> \param fdata ...
! **************************************************************************************************
   SUBROUTINE readbuffer_d(psockfd, fdata)
      INTEGER, INTENT(IN)                                :: psockfd
      REAL(KIND=dp), INTENT(OUT)                         :: fdata

      CHARACTER(len=*), PARAMETER                        :: routineN = 'readbuffer_d'

      INTEGER                                            :: handle
      REAL(KIND=C_DOUBLE), TARGET                        :: cdata

      CALL timeset(routineN, handle)

      CALL readbuffer_csocket(psockfd, c_loc(cdata), 8)
      fdata = cdata

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param psockfd ...
!> \param fdata ...
! **************************************************************************************************
   SUBROUTINE readbuffer_i(psockfd, fdata)
      INTEGER, INTENT(IN)                                :: psockfd
      INTEGER, INTENT(OUT)                               :: fdata

      CHARACTER(len=*), PARAMETER                        :: routineN = 'readbuffer_i'

      INTEGER                                            :: handle
      INTEGER(KIND=C_INT), TARGET                        :: cdata

      CALL timeset(routineN, handle)

      CALL readbuffer_csocket(psockfd, c_loc(cdata), 4)
      fdata = cdata

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param psockfd ...
!> \param fstring ...
!> \param plen ...
! **************************************************************************************************
   SUBROUTINE readbuffer_s(psockfd, fstring, plen)
      INTEGER, INTENT(IN)                                :: psockfd
      CHARACTER(LEN=*), INTENT(OUT)                      :: fstring
      INTEGER, INTENT(IN)                                :: plen

      CHARACTER(len=*), PARAMETER                        :: routineN = 'readbuffer_s'

      INTEGER                                            :: handle, i
      CHARACTER(LEN=1, KIND=C_CHAR), TARGET              :: cstring(plen)

      CALL timeset(routineN, handle)

      CALL readbuffer_csocket(psockfd, c_loc(cstring(1)), plen)
      fstring = ""
      DO i = 1, plen
         fstring(i:i) = cstring(i)
      END DO

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param psockfd ...
!> \param fdata ...
!> \param plen ...
! **************************************************************************************************
   SUBROUTINE readbuffer_dv(psockfd, fdata, plen)
      INTEGER, INTENT(IN)                                :: psockfd, plen
      REAL(KIND=dp), INTENT(OUT), TARGET                 :: fdata(plen)

      CHARACTER(len=*), PARAMETER                        :: routineN = 'readbuffer_dv'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL readbuffer_csocket(psockfd, c_loc(fdata(1)), 8*plen)

      CALL timestop(handle)

   END SUBROUTINE
#endif

END MODULE sockets_interface
