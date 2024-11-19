!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Interface to ACE C wrapper.
! **************************************************************************************************

MODULE ace_wrapper
   USE ISO_C_BINDING,                   ONLY: C_CHAR,&
                                              C_DOUBLE,&
                                              C_INT,&
                                              C_NULL_CHAR,&
                                              C_NULL_PTR,&
                                              C_PTR
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: ace_model_type, ace_model_compute, ace_model_initialize, ace_model_release

   TYPE ace_model_type
      TYPE(C_PTR)                                 :: c_ptr = C_NULL_PTR
      CHARACTER(LEN=2), DIMENSION(:), ALLOCATABLE :: symbolc
   END TYPE ace_model_type

   INTERFACE ace_interface
! **************************************************************************************************
!> \brief ...
!> \param natomc ...
!> \param nghostc ...
!> \param neic ...
!> \param neiatc ...
!> \param originc ...
!> \param nlistc ...
!> \param attypec ...
!> \param atposc ...
!> \param forcec ...
!> \param virialc ...
!> \param energyc ...
!> \param acedata_ptr ...
! **************************************************************************************************
      SUBROUTINE AcePotCompute(natomc, nghostc, neic, neiatc, originc, nlistc, &
                               attypec, atposc, forcec, virialc, energyc, &
                               acedata_ptr) BIND(C, NAME="AcePotCompute")

         IMPORT :: C_PTR, C_INT, C_DOUBLE
         INTEGER(C_INT), VALUE                                           :: natomc
         INTEGER(C_INT), VALUE                                           :: nghostc
         INTEGER(C_INT), VALUE                                           :: neic
         INTEGER(C_INT), DIMENSION(0:natomc)                             :: neiatc
         INTEGER(C_INT), DIMENSION(1:nghostc)                            :: originc
         INTEGER(C_INT), DIMENSION(1:neic)                               :: nlistc
         INTEGER(C_INT), DIMENSION(1:natomc + nghostc)                   :: attypec
         REAL(C_DOUBLE), DIMENSION(1:3*(natomc + nghostc))               :: atposc
         REAL(C_DOUBLE), DIMENSION(1:3*natomc)                           :: forcec
         REAL(C_DOUBLE), DIMENSION(1:6)                                  :: virialc
         REAL(C_DOUBLE), DIMENSION(1:natomc)                             :: energyc
         TYPE(C_PTR)                                                     :: acedata_ptr
      END SUBROUTINE
! **************************************************************************************************
!> \brief ...
!> \param ntypec ...
!> \param symbolc ...
!> \param nlen ...
!> \param cname ...
!> \param rcutc ...
!> \param acedata_ptr ...
! **************************************************************************************************
      SUBROUTINE AcePotInitialize(ntypec, symbolc, nlen, cname, rcutc, &
                                  acedata_ptr) BIND(C, NAME="AcePotInitialize")
         USE ISO_C_BINDING, ONLY: C_CHAR, C_INT, C_DOUBLE, C_PTR
      INTEGER(C_INT), VALUE                              :: ntypec
      CHARACTER(KIND=C_CHAR, LEN=1), DIMENSION(*)        :: symbolc
      INTEGER(C_INT), VALUE                              :: nlen
      CHARACTER(KIND=C_CHAR, LEN=1), DIMENSION(*)        :: cname
      REAL(C_DOUBLE), DIMENSION(1:ntypec, 1:ntypec)      :: rcutc
      TYPE(C_PTR)                                        :: acedata_ptr

      END SUBROUTINE AcePotInitialize
! **************************************************************************************************
!> \brief ...
!> \param acedata_ptr ...
! **************************************************************************************************
      SUBROUTINE AcePotFinalize(acedata_ptr) BIND(C, NAME="AcePotFinalize")
         USE ISO_C_BINDING, ONLY: C_PTR
      TYPE(C_PTR)                                        :: acedata_ptr

      END SUBROUTINE AcePotFinalize
   END INTERFACE ace_interface

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param natomc ...
!> \param nghostc ...
!> \param neic ...
!> \param neiatc ...
!> \param originc ...
!> \param nlistc ...
!> \param attypec ...
!> \param atposc ...
!> \param forcec ...
!> \param virialc ...
!> \param energyc ...
!> \param model ...
! **************************************************************************************************
   SUBROUTINE ace_model_compute(natomc, nghostc, neic, neiatc, originc, nlistc, &
                                attypec, atposc, forcec, virialc, energyc, model)

      INTEGER                                            :: natomc, nghostc, neic
      INTEGER, DIMENSION(0:natomc), INTENT(IN)           :: neiatc
      INTEGER, DIMENSION(1:nghostc), INTENT(IN)          :: originc
      INTEGER, DIMENSION(1:neic), INTENT(IN)             :: nlistc
      INTEGER, DIMENSION(1:natomc + nghostc), INTENT(IN) :: attypec
      REAL(KIND=8), DIMENSION(1:3*(natomc + nghostc)), &
         INTENT(IN)                                      :: atposc
      REAL(KIND=8), DIMENSION(1:3*natomc), INTENT(INOUT) :: forcec
      REAL(KIND=8), DIMENSION(1:6), INTENT(INOUT)        :: virialc
      REAL(KIND=8), DIMENSION(1:natomc), INTENT(INOUT)   :: energyc
      TYPE(ace_model_type), INTENT(IN)                   :: model

#if defined(__ACE)

      CALL AcePotCompute(natomc, nghostc, neic, neiatc, originc, nlistc, attypec, &
                         atposc, forcec, virialc, energyc, model%c_ptr)
#else
      MARK_USED(natomc)
      MARK_USED(nghostc)
      MARK_USED(neic)
      MARK_USED(neiatc)
      MARK_USED(originc)
      MARK_USED(nlistc)
      MARK_USED(attypec)
      MARK_USED(atposc)
      MARK_USED(forcec)
      MARK_USED(virialc)
      MARK_USED(energyc)
      MARK_USED(model)
      CPABORT("CP2K was compiled without ACE library.")
#endif
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param ntypec ...
!> \param symbolc ...
!> \param fname ...
!> \param rcutc ...
!> \param model ...
! **************************************************************************************************
   SUBROUTINE ace_model_initialize(ntypec, symbolc, fname, rcutc, model)
      USE ISO_C_BINDING, ONLY: C_CHAR, C_NULL_CHAR
      INTEGER, INTENT(IN)                                :: ntypec
      CHARACTER(KIND=C_CHAR, LEN=2), &
         DIMENSION(1:ntypec), INTENT(IN)                 :: symbolc
      CHARACTER(KIND=C_CHAR, LEN=*), INTENT(IN)          :: fname
      REAL(KIND=8), DIMENSION(1:ntypec, 1:ntypec), &
         INTENT(OUT)                                     :: rcutc
      TYPE(ace_model_type), INTENT(OUT)                  :: model

#if defined(__ACE)

      CHARACTER(KIND=C_CHAR), DIMENSION(LEN(fname) + 1)  :: cname
      INTEGER                                            :: i, nlen

      nlen = LEN(fname)
      DO i = 1, nlen
         cname(i) = fname(i:i)
      END DO
      nlen = nlen + 1
      cname(nlen) = C_NULL_CHAR
      CALL AcePotInitialize(ntypec, symbolc, nlen, cname, rcutc, model%c_ptr)
      ALLOCATE (model%symbolc(ntypec))
      DO i = 1, ntypec
         model%symbolc(i) = symbolc(i)
      END DO
#else
      MARK_USED(ntypec)
      MARK_USED(symbolc)
      MARK_USED(fname)
      MARK_USED(rcutc)
      MARK_USED(model)
      CPABORT("CP2K was compiled without ACE library.")
#endif
   END SUBROUTINE ace_model_initialize

!**************************************************************************************************
!> \brief Releases an ACE model and all its ressources.
!> \param model Pointer to the ACE model.
! **************************************************************************************************
   SUBROUTINE ace_model_release(model)
      USE ISO_C_BINDING, ONLY: C_NULL_PTR
      TYPE(ace_model_type), INTENT(INOUT)                :: model

#if defined(__ACE)

      CALL AcePotFinalize(model%c_ptr)
      model%c_ptr = C_NULL_PTR
      IF (ALLOCATED(model%symbolc)) DEALLOCATE (model%symbolc)
#else
      MARK_USED(model)
      CPABORT("CP2K was compiled without ACE library.")
#endif

   END SUBROUTINE ace_model_release

END MODULE ace_wrapper

