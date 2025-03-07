!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
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
                                              combine_basis_sets,&
                                              copy_gto_basis_set,&
                                              create_primitive_basis_set,&
                                              deallocate_gto_basis_set,&
                                              get_gto_basis_set,&
                                              gto_basis_set_type
   USE kinds,                           ONLY: dp
   USE orbital_transformation_matrices, ONLY: init_spherical_harmonics
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

      INTEGER                                            :: i, ipgf, iset, j, l, lbas, maxl, maxlo, &
                                                            maxls, mp, n1, n2, nn, nseto, nsets
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nps
      INTEGER, DIMENSION(:), POINTER                     :: lmaxo, lmaxs, lmino, lmins, npgfo, npgfs
      REAL(KIND=dp)                                      :: fmin, fr1, fr2, fz, rr, xz, yz, zmall, &
                                                            zms, zz
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: z1, z2, zmaxs
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: zeta, zexs
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: zeto, zets
      TYPE(gto_basis_set_type), POINTER                  :: ext_basis, p_basis

      CPASSERT(.NOT. ASSOCIATED(gapw_1c_basis))

      IF (basis_1c_level == 0) THEN
         ! we use the orbital basis set
         CALL copy_gto_basis_set(orb_basis, gapw_1c_basis)
      ELSE
         CALL copy_gto_basis_set(orb_basis, gapw_1c_basis)
         NULLIFY (ext_basis)
         CALL allocate_gto_basis_set(ext_basis)
         ! get information on orbital basis and soft basis
         CALL get_gto_basis_set(gto_basis_set=orb_basis, maxl=maxlo, nset=nseto, &
                                lmin=lmino, lmax=lmaxo, npgf=npgfo, zet=zeto)
         CALL get_gto_basis_set(gto_basis_set=soft_basis, maxl=maxls, nset=nsets, &
                                lmin=lmins, lmax=lmaxs, npgf=npgfs, zet=zets)
         ! determine max soft exponent per l-qn
         maxl = MAX(maxls, maxlo)
         ALLOCATE (zmaxs(0:maxl), nps(0:maxl))
         zmaxs = 0.0_dp
         DO iset = 1, nsets
            zms = MAXVAL(zets(:, iset))
            DO l = lmins(iset), lmaxs(iset)
               zmaxs(l) = MAX(zmaxs(l), zms)
            END DO
         END DO
         zmall = MAXVAL(zmaxs)
         ! in case of missing soft basis!
         zmall = MAX(zmall, 0.20_dp)
         ! create pools of exponents for each l-qn
         nps = 0
         DO iset = 1, nsets
            DO l = lmins(iset), lmaxs(iset)
               nps(l) = nps(l) + npgfs(iset)
            END DO
         END DO
         mp = MAXVAL(nps)
         ALLOCATE (zexs(1:mp, 0:maxl))
         zexs = 0.0_dp
         nps = 0
         DO iset = 1, nsets
            DO ipgf = 1, npgfs(iset)
               DO l = lmins(iset), lmaxs(iset)
                  nps(l) = nps(l) + 1
                  zexs(nps(l), l) = zets(ipgf, iset)
               END DO
            END DO
         END DO

         SELECT CASE (basis_1c_level)
         CASE (1)
            lbas = maxl
            fr1 = 2.50_dp
            fr2 = 2.50_dp
         CASE (2)
            lbas = maxl
            fr1 = 2.00_dp
            fr2 = 2.50_dp
         CASE (3)
            lbas = maxl + 1
            fr1 = 1.75_dp
            fr2 = 2.50_dp
         CASE (4)
            lbas = maxl + 2
            fr1 = 1.50_dp
            fr2 = 2.50_dp
         CASE DEFAULT
            CPABORT("unknown case")
         END SELECT
         lbas = MIN(lbas, 5)
         !
         CALL init_spherical_harmonics(lbas, 0)
         !
         rr = LOG(zmall/0.05_dp)/LOG(fr1)
         n1 = INT(rr) + 1
         rr = LOG(zmall/0.05_dp)/LOG(fr2)
         n2 = INT(rr) + 1
         ALLOCATE (z1(n1), z2(n2))
         z1 = 0.0_dp
         zz = zmall*SQRT(fr1)
         DO i = 1, n1
            z1(i) = zz/(fr1**(i - 1))
         END DO
         z2 = 0.0_dp
         zz = zmall
         DO i = 1, n2
            z2(i) = zz/(fr2**(i - 1))
         END DO
         ALLOCATE (zeta(MAX(n1, n2), lbas + 1))
         zeta = 0.0_dp
         !
         ext_basis%nset = lbas + 1
         ALLOCATE (ext_basis%lmin(lbas + 1), ext_basis%lmax(lbas + 1))
         ALLOCATE (ext_basis%npgf(lbas + 1))
         DO l = 0, lbas
            ext_basis%lmin(l + 1) = l
            ext_basis%lmax(l + 1) = l
            IF (l <= maxl) THEN
               fmin = 10.0_dp
               nn = 0
               DO i = 1, n1
                  xz = z1(i)
                  DO j = 1, nps(l)
                     yz = zexs(j, l)
                     fz = MAX(xz, yz)/MIN(xz, yz)
                     fmin = MIN(fmin, fz)
                  END DO
                  IF (fmin > fr1**0.25) THEN
                     nn = nn + 1
                     zeta(nn, l + 1) = xz
                  END IF
               END DO
               CPASSERT(nn > 0)
               ext_basis%npgf(l + 1) = nn
            ELSE
               ext_basis%npgf(l + 1) = n2
               zeta(1:n2, l + 1) = z2(1:n2)
            END IF
         END DO
         nn = MAXVAL(ext_basis%npgf)
         ALLOCATE (ext_basis%zet(nn, lbas + 1))
         DO i = 1, lbas + 1
            nn = ext_basis%npgf(i)
            ext_basis%zet(1:nn, i) = zeta(1:nn, i)
         END DO
         ext_basis%name = "extbas"
         ext_basis%kind_radius = orb_basis%kind_radius
         ext_basis%short_kind_radius = orb_basis%short_kind_radius
         ext_basis%norm_type = orb_basis%norm_type

         NULLIFY (p_basis)
         CALL create_primitive_basis_set(ext_basis, p_basis)
         CALL combine_basis_sets(gapw_1c_basis, p_basis)

         CALL deallocate_gto_basis_set(ext_basis)
         CALL deallocate_gto_basis_set(p_basis)
         DEALLOCATE (zmaxs, zexs, nps, z1, z2, zeta)
      END IF

   END SUBROUTINE create_1c_basis

END MODULE gapw_1c_basis_set
