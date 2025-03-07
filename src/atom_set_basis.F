!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
MODULE atom_set_basis
   USE ai_onecenter,                    ONLY: sg_overlap
   USE atom_types,                      ONLY: CGTO_BASIS,&
                                              GTO_BASIS,&
                                              atom_basis_type,&
                                              lmat
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE input_constants,                 ONLY: do_gapw_log
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: dfac,&
                                              twopi
   USE qs_grid_atom,                    ONLY: allocate_grid_atom,&
                                              create_grid_atom,&
                                              grid_atom_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atom_set_basis'

   INTEGER, PARAMETER                                 :: nua = 40, nup = 20
   REAL(KIND=dp), DIMENSION(nua), PARAMETER :: ugbs = (/0.007299_dp, 0.013705_dp, 0.025733_dp, &
                                        0.048316_dp, 0.090718_dp, 0.170333_dp, 0.319819_dp, 0.600496_dp, 1.127497_dp, 2.117000_dp, &
                                                 3.974902_dp, 7.463317_dp, 14.013204_dp, 26.311339_dp, 49.402449_dp, 92.758561_dp, &
                                                      174.164456_dp, 327.013024_dp, 614.003114_dp, 1152.858743_dp, 2164.619772_dp, &
                                                4064.312984_dp, 7631.197056_dp, 14328.416324_dp, 26903.186074_dp, 50513.706789_dp, &
                                         94845.070265_dp, 178082.107320_dp, 334368.848683_dp, 627814.487663_dp, 1178791.123851_dp, &
                                                      2213310.684886_dp, 4155735.557141_dp, 7802853.046713_dp, 14650719.428954_dp, &
                                                  27508345.793637_dp, 51649961.080194_dp, 96978513.342764_dp, 182087882.613702_dp, &
                                                        341890134.751331_dp/)

   PUBLIC :: set_kind_basis_atomic

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param basis ...
!> \param orb_basis_set ...
!> \param has_pp ...
!> \param agrid ...
!> \param cp2k_norm ...
! **************************************************************************************************
   SUBROUTINE set_kind_basis_atomic(basis, orb_basis_set, has_pp, agrid, cp2k_norm)
      TYPE(atom_basis_type), INTENT(INOUT)               :: basis
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      LOGICAL, INTENT(IN)                                :: has_pp
      TYPE(grid_atom_type), OPTIONAL                     :: agrid
      LOGICAL, INTENT(IN), OPTIONAL                      :: cp2k_norm

      INTEGER                                            :: i, ii, ipgf, j, k, l, m, ngp, nj, nr, &
                                                            ns, nset, nsgf, nu, quadtype
      INTEGER, DIMENSION(:), POINTER                     :: lmax, lmin, npgf, nshell
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf, last_sgf, ls
      LOGICAL                                            :: has_basis, set_norm
      REAL(KIND=dp)                                      :: al, an, cn, ear, en, rk
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: zet
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: gcc
      TYPE(grid_atom_type), POINTER                      :: grid

      IF (ASSOCIATED(orb_basis_set)) THEN
         has_basis = .TRUE.
      ELSE
         has_basis = .FALSE.
      END IF

      IF (PRESENT(cp2k_norm)) THEN
         set_norm = cp2k_norm
      ELSE
         set_norm = .FALSE.
      END IF

      NULLIFY (grid)
      IF (PRESENT(agrid)) THEN
         ngp = agrid%nr
         quadtype = agrid%quadrature
      ELSE
         ngp = 400
         quadtype = do_gapw_log
      END IF
      CALL allocate_grid_atom(grid)
      CALL create_grid_atom(grid, ngp, 1, 1, 0, quadtype)
      grid%nr = ngp
      basis%grid => grid

      NULLIFY (basis%am, basis%cm, basis%as, basis%ns, basis%bf, basis%dbf, basis%ddbf)

      IF (has_basis) THEN
         ! fill in the basis data structures
         basis%basis_type = CGTO_BASIS
         basis%eps_eig = 1.e-12_dp
         CALL get_gto_basis_set(orb_basis_set, &
                                nset=nset, nshell=nshell, npgf=npgf, lmin=lmin, lmax=lmax, &
                                l=ls, nsgf=nsgf, zet=zet, gcc=gcc, &
                                first_sgf=first_sgf, last_sgf=last_sgf)
         basis%nprim = 0
         basis%nbas = 0
         DO i = 1, nset
            DO j = lmin(i), MIN(lmax(i), lmat)
               basis%nprim(j) = basis%nprim(j) + npgf(i)
            END DO
            DO j = 1, nshell(i)
               l = ls(j, i)
               IF (l <= lmat) THEN
                  basis%nbas(l) = basis%nbas(l) + 1
                  k = basis%nbas(l)
               END IF
            END DO
         END DO

         nj = MAXVAL(basis%nprim)
         ns = MAXVAL(basis%nbas)
         ALLOCATE (basis%am(nj, 0:lmat))
         basis%am = 0._dp
         ALLOCATE (basis%cm(nj, ns, 0:lmat))
         basis%cm = 0._dp
         DO j = 0, lmat
            nj = 0
            ns = 0
            cn = 2.0_dp**(j + 2)/SQRT(dfac(2*j + 1))/twopi**0.25_dp
            en = (2*j + 3)*0.25_dp
            DO i = 1, nset
               IF (j >= lmin(i) .AND. j <= lmax(i)) THEN
                  DO ipgf = 1, npgf(i)
                     basis%am(nj + ipgf, j) = zet(ipgf, i)
                  END DO
                  DO ii = 1, nshell(i)
                     IF (ls(ii, i) == j) THEN
                        ns = ns + 1
                        IF (set_norm) THEN
                           DO ipgf = 1, npgf(i)
                              an = cn*zet(ipgf, i)**en
                              basis%cm(nj + ipgf, ns, j) = an*gcc(ipgf, ii, i)
                           END DO
                        ELSE
                           DO ipgf = 1, npgf(i)
                              basis%cm(nj + ipgf, ns, j) = gcc(ipgf, ii, i)
                           END DO
                        END IF
                     END IF
                  END DO
                  nj = nj + npgf(i)
               END IF
            END DO
         END DO
         ! Normalization
         IF (set_norm) THEN
            CALL normalize_basis_cp2k(basis)
         END IF
      ELSE
         ! use default basis
         IF (has_pp) THEN
            nu = nup
         ELSE
            nu = nua
         END IF
         basis%geometrical = .FALSE.
         basis%aval = 0._dp
         basis%cval = 0._dp
         basis%start = 0
         basis%eps_eig = 1.e-12_dp

         basis%basis_type = GTO_BASIS
         basis%nbas = nu
         basis%nprim = nu
         ALLOCATE (basis%am(nu, 0:lmat))
         DO i = 0, lmat
            basis%am(1:nu, i) = ugbs(1:nu)
         END DO
      END IF

      ! initialize basis function on a radial grid
      nr = basis%grid%nr
      m = MAXVAL(basis%nbas)
      ALLOCATE (basis%bf(nr, m, 0:lmat))
      ALLOCATE (basis%dbf(nr, m, 0:lmat))
      ALLOCATE (basis%ddbf(nr, m, 0:lmat))

      basis%bf = 0._dp
      basis%dbf = 0._dp
      basis%ddbf = 0._dp
      DO l = 0, lmat
         DO i = 1, basis%nprim(l)
            al = basis%am(i, l)
            IF (basis%basis_type == GTO_BASIS) THEN
               DO k = 1, nr
                  rk = basis%grid%rad(k)
                  ear = EXP(-al*basis%grid%rad(k)**2)
                  basis%bf(k, i, l) = rk**l*ear
                  basis%dbf(k, i, l) = (REAL(l, dp)*rk**(l - 1) - 2._dp*al*rk**(l + 1))*ear
                  basis%ddbf(k, i, l) = (REAL(l*(l - 1), dp)*rk**(l - 2) - &
                                         2._dp*al*REAL(2*l + 1, dp)*rk**(l) + 4._dp*al*rk**(l + 2))*ear
               END DO
            ELSEIF (basis%basis_type == CGTO_BASIS) THEN
               DO k = 1, nr
                  rk = basis%grid%rad(k)
                  ear = EXP(-al*basis%grid%rad(k)**2)
                  DO j = 1, basis%nbas(l)
                     basis%bf(k, j, l) = basis%bf(k, j, l) + rk**l*ear*basis%cm(i, j, l)
                     basis%dbf(k, j, l) = basis%dbf(k, j, l) &
                                          + (REAL(l, dp)*rk**(l - 1) - 2._dp*al*rk**(l + 1))*ear*basis%cm(i, j, l)
                     basis%ddbf(k, j, l) = basis%ddbf(k, j, l) + &
                                           (REAL(l*(l - 1), dp)*rk**(l - 2) - 2._dp*al*REAL(2*l + 1, dp)*rk**(l) + &
                                            4._dp*al*rk**(l + 2))*ear*basis%cm(i, j, l)
                  END DO
               END DO
            ELSE
               CPABORT('Atom basis type?')
            END IF
         END DO
      END DO

   END SUBROUTINE set_kind_basis_atomic

! **************************************************************************************************
!> \brief ...
!> \param basis ...
! **************************************************************************************************
   SUBROUTINE normalize_basis_cp2k(basis)
      TYPE(atom_basis_type), INTENT(INOUT)               :: basis

      INTEGER                                            :: ii, l, n, np
      REAL(KIND=dp)                                      :: fnorm
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: smat

      DO l = 0, lmat
         n = basis%nbas(l)
         np = basis%nprim(l)
         IF (n > 0) THEN
            ALLOCATE (smat(np, np))
            CALL sg_overlap(smat, l, basis%am(1:np, l), basis%am(1:np, l))
            DO ii = 1, basis%nbas(l)
               fnorm = DOT_PRODUCT(basis%cm(1:np, ii, l), MATMUL(smat, basis%cm(1:np, ii, l)))
               fnorm = 1._dp/SQRT(fnorm)
               basis%cm(1:np, ii, l) = fnorm*basis%cm(1:np, ii, l)
            END DO
            DEALLOCATE (smat)
         END IF
      END DO

   END SUBROUTINE normalize_basis_cp2k

! **************************************************************************************************

END MODULE atom_set_basis
