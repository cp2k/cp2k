!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \par History
!>      none
!> \author JHU (9.2022)
! **************************************************************************************************
MODULE gapw_1c_basis_set

   USE basis_set_types,                 ONLY: allocate_gto_basis_set,&
                                              copy_gto_basis_set,&
                                              get_gto_basis_set,&
                                              gto_basis_set_type
   USE kinds,                           ONLY: dp
#include "base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters (only in this module)

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gapw_1c_basis_set'

   INTEGER, PARAMETER :: max_name_length = 60

! *** Public subroutines ***

   PUBLIC :: create_1c_basis

CONTAINS

! **************************************************************************************************
!> \brief   create the one center basis from the orbital basis
!> \param orb_basis ...
!> \param soft_basis ...
!> \param gapw_1c_basis ...
!> \param basis_1c_level ...
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE create_1c_basis(orb_basis, soft_basis, gapw_1c_basis, basis_1c_level)

      TYPE(gto_basis_set_type), POINTER                  :: orb_basis, soft_basis, gapw_1c_basis
      INTEGER, INTENT(IN)                                :: basis_1c_level

      INTEGER                                            :: i, ipgf, iset, l, maxl, maxlo, maxls, &
                                                            mp, nseto, nsets
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nph, nps
      INTEGER, DIMENSION(:), POINTER                     :: lmaxo, lmaxs, lmino, lmins, npgfo, npgfs
      REAL(KIND=dp)                                      :: dex, zms
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: zmaxs
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: zexh, zexs
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: zeto, zets

      CPASSERT(.NOT. ASSOCIATED(gapw_1c_basis))

      CALL allocate_gto_basis_set(gapw_1c_basis)
      ! get information on orbital basis and soft basis
      CALL get_gto_basis_set(gto_basis_set=orb_basis, maxl=maxlo, nset=nseto, &
                             lmin=lmino, lmax=lmaxo, npgf=npgfo, zet=zeto)
      CALL get_gto_basis_set(gto_basis_set=soft_basis, maxl=maxls, nset=nsets, &
                             lmin=lmins, lmax=lmaxs, npgf=npgfs, zet=zets)
      ! determine max soft exponent per l-qn
      maxl = MAX(maxls, maxlo)
      ALLOCATE (zmaxs(0:maxl), nps(0:maxl), nph(0:maxl))
      zmaxs = 0.0_dp
      DO iset = 1, nsets
         zms = MAXVAL(zets(:, iset))
         DO l = lmins(iset), lmaxs(iset)
            zmaxs(l) = MAX(zmaxs(l), zms)
         END DO
      END DO
      ! create pools of exponents (hard,soft) for each l-qn
      nps = 0
      DO iset = 1, nsets
         DO l = lmins(iset), lmaxs(iset)
            nps(l) = nps(l) + npgfs(iset)
         END DO
      END DO
      mp = MAXVAL(nps)
      ALLOCATE (zexs(1:mp, 0:maxl))
      zexs = 0.0_dp
      nph = 0
      DO iset = 1, nseto
         DO l = lmino(iset), lmaxo(iset)
            nph(l) = nph(l) + npgfo(iset)
         END DO
      END DO
      mp = MAXVAL(nph)
      ALLOCATE (zexh(1:mp, 0:maxl))
      zexh = 0.0_dp
      nps = 0
      nph = 0
      DO iset = 1, nseto
         DO ipgf = 1, npgfo(iset)
            DO l = lmino(iset), lmaxo(iset)
               IF (zeto(ipgf, iset) > zmaxs(l)) THEN
                  DO i = 1, nph(l)
                     dex = ABS(zeto(ipgf, iset) - zexh(i, l))/zeto(ipgf, iset)
                     IF (dex < 1.E-3_dp) EXIT
                  END DO
                  IF (dex > 1.E-3_dp) THEN
                     nph(l) = nph(l) + 1
                     zexh(nph(l), l) = zeto(ipgf, iset)
                  END IF
               ELSE
                  DO i = 1, nps(l)
                     dex = ABS(zeto(ipgf, iset) - zexs(i, l))/zeto(ipgf, iset)
                     IF (dex < 1.E-3_dp) EXIT
                  END DO
                  IF (dex > 1.E-3_dp) THEN
                     nph(l) = nph(l) + 1
                     zexs(nph(l), l) = zeto(ipgf, iset)
                  END IF
               END IF
            END DO
         END DO
      END DO

      SELECT CASE (basis_1c_level)
      CASE (1)
      CASE (2)
      CASE (3)
      CASE (4)
      CASE DEFAULT
      END SELECT

      CPWARN("Option not implemented, use Orb Basis for GAPW 1C")
      CALL copy_gto_basis_set(orb_basis, gapw_1c_basis)

      DEALLOCATE (zmaxs, zexs, zexh, nph, nps)

   END SUBROUTINE create_1c_basis

END MODULE gapw_1c_basis_set
