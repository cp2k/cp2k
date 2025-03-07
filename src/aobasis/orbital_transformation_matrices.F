!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of the spherical harmonics and the corresponding orbital
!>        transformation matrices.
!> \par Literature
!>      H. B. Schlegel, M. J. Frisch, Int. J. Quantum Chem. 54, 83 (1995)
!> \par History
!>      - restructured and cleaned (20.05.2004,MK)
!> \author Matthias Krack (08.06.2000)
! **************************************************************************************************
MODULE orbital_transformation_matrices

   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: dfac,&
                                              fac,&
                                              pi
   USE mathlib,                         ONLY: binomial
   USE orbital_pointers,                ONLY: co,&
                                              nco,&
                                              ncoset,&
                                              nso,&
                                              nsoset
   USE orbital_symbols,                 ONLY: cgf_symbol,&
                                              sgf_symbol

!$ USE OMP_LIB, ONLY: omp_get_level

#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'orbital_transformation_matrices'

   TYPE orbtramat_type
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: c2s => NULL(), slm => NULL(), &
                                                 slm_inv => NULL(), s2c => NULL()
   END TYPE orbtramat_type

   TYPE(orbtramat_type), DIMENSION(:), POINTER :: orbtramat => NULL()

   TYPE orbrotmat_type
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: mat => NULL()
   END TYPE orbrotmat_type

   INTEGER, SAVE :: current_maxl = -1

   PUBLIC :: orbtramat
   PUBLIC :: orbrotmat_type, calculate_rotmat, release_rotmat
   PUBLIC :: deallocate_spherical_harmonics, init_spherical_harmonics

CONTAINS

! **************************************************************************************************
!> \brief  Allocate and initialize the orbital transformation matrices for
!>         the transformation and for the backtransformation of orbitals
!>         from the Cartesian representation to the spherical representation.
!> \param maxl ...
!> \date    20.05.2004
!> \author  MK
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE create_spherical_harmonics(maxl)

      INTEGER, INTENT(IN)                                :: maxl

      INTEGER                                            :: expo, i, ic, ic1, ic2, is, j, k, l, lx, &
                                                            lx1, lx2, ly, ly1, ly2, lz, lz1, lz2, &
                                                            m, ma, nc, ns
      REAL(KIND=dp)                                      :: s, s1, s2

!$    IF (omp_get_level() > 0) &
!$       CPABORT("create_spherical_harmonics is not thread-safe")

      IF (current_maxl > -1) THEN
         CALL cp_abort(__LOCATION__, &
                       "Spherical harmonics are already allocated. "// &
                       "Use the init routine for an update")
      END IF

      IF (maxl < 0) THEN
         CALL cp_abort(__LOCATION__, &
                       "A negative maximum angular momentum quantum "// &
                       "number is invalid")
      END IF

      ALLOCATE (orbtramat(0:maxl))
      nc = ncoset(maxl)
      ns = nsoset(maxl)

      DO l = 0, maxl
         nc = nco(l)
         ns = nso(l)
         ALLOCATE (orbtramat(l)%c2s(ns, nc))
         orbtramat(l)%c2s = 0.0_dp
         ALLOCATE (orbtramat(l)%s2c(ns, nc))
         orbtramat(l)%s2c = 0.0_dp
         ALLOCATE (orbtramat(l)%slm(ns, nc))
         orbtramat(l)%slm = 0.0_dp
         ALLOCATE (orbtramat(l)%slm_inv(ns, nc))
         orbtramat(l)%slm_inv = 0.0_dp
      END DO

      ! Loop over all angular momentum quantum numbers l

      DO l = 0, maxl

         ! Build the orbital transformation matrix for the transformation
         ! from Cartesian to spherical orbitals (c2s, formula 15)

         DO lx = 0, l
            DO ly = 0, l - lx
               lz = l - lx - ly
               ic = co(lx, ly, lz)
               DO m = -l, l
                  is = l + m + 1
                  ma = ABS(m)
                  j = lx + ly - ma
                  IF ((j >= 0) .AND. (MODULO(j, 2) == 0)) THEN
                     j = j/2
                     s1 = 0.0_dp
                     DO i = 0, (l - ma)/2
                        s2 = 0.0_dp
                        DO k = 0, j
                           IF (((m < 0) .AND. (MODULO(ABS(ma - lx), 2) == 1)) .OR. &
                               ((m > 0) .AND. (MODULO(ABS(ma - lx), 2) == 0))) THEN
                              expo = (ma - lx + 2*k)/2
                              s = (-1.0_dp)**expo*SQRT(2.0_dp)
                           ELSE IF ((m == 0) .AND. (MODULO(lx, 2) == 0)) THEN
                              expo = k - lx/2
                              s = (-1.0_dp)**expo
                           ELSE
                              s = 0.0_dp
                           END IF
                           s2 = s2 + binomial(j, k)*binomial(ma, lx - 2*k)*s
                        END DO
                        s1 = s1 + binomial(l, i)*binomial(i, j)* &
                             (-1.0_dp)**i*fac(2*l - 2*i)/fac(l - ma - 2*i)*s2
                     END DO
                     orbtramat(l)%c2s(is, ic) = &
                        SQRT((fac(2*lx)*fac(2*ly)*fac(2*lz)*fac(l)*fac(l - ma))/ &
                             (fac(lx)*fac(ly)*fac(lz)*fac(2*l)*fac(l + ma)))*s1/ &
                        (2.0_dp**l*fac(l))
                  ELSE
                     orbtramat(l)%c2s(is, ic) = 0.0_dp
                  END IF
               END DO
            END DO
         END DO

         ! Build the corresponding transformation matrix for the transformation from
         ! spherical to Cartesian orbitals (s2c = s*TRANSPOSE(c2s), formulas 18 and 19)

         DO lx1 = 0, l
            DO ly1 = 0, l - lx1
               lz1 = l - lx1 - ly1
               ic1 = co(lx1, ly1, lz1)
               s1 = SQRT((fac(lx1)*fac(ly1)*fac(lz1))/ &
                         (fac(2*lx1)*fac(2*ly1)*fac(2*lz1)))
               DO lx2 = 0, l
                  DO ly2 = 0, l - lx2
                     lz2 = l - lx2 - ly2
                     ic2 = co(lx2, ly2, lz2)
                     lx = lx1 + lx2
                     ly = ly1 + ly2
                     lz = lz1 + lz2
                     IF ((MODULO(lx, 2) == 0) .AND. &
                         (MODULO(ly, 2) == 0) .AND. &
                         (MODULO(lz, 2) == 0)) THEN
                        s2 = SQRT((fac(lx2)*fac(ly2)*fac(lz2))/ &
                                  (fac(2*lx2)*fac(2*ly2)*fac(2*lz2)))
                        s = fac(lx)*fac(ly)*fac(lz)*s1*s2/ &
                            (fac(lx/2)*fac(ly/2)*fac(lz/2))
                        DO is = 1, nso(l)
                           orbtramat(l)%s2c(is, ic1) = orbtramat(l)%s2c(is, ic1) + &
                                                       s*orbtramat(l)%c2s(is, ic2)
                        END DO
                     END IF
                  END DO
               END DO
            END DO
         END DO

         ! Build up the real spherical harmonics

         s = SQRT(0.25_dp*dfac(2*l + 1)/pi)

         DO lx = 0, l
            DO ly = 0, l - lx
               lz = l - lx - ly
               ic = co(lx, ly, lz)
               DO m = -l, l
                  is = l + m + 1
                  s1 = SQRT(dfac(2*lx - 1)*dfac(2*ly - 1)*dfac(2*lz - 1))
                  orbtramat(l)%slm(is, ic) = s*orbtramat(l)%c2s(is, ic)/s1
                  !MK s2 = (-1.0_dp)**m*s ! alternative S(lm) definition
                  !MK orbtramat(l)%slm(is, ic) = s2*orbtramat(l)%c2s(is,ic)/s1
                  orbtramat(l)%slm_inv(is, ic) = s1*orbtramat(l)%s2c(is, ic)/s
               END DO
            END DO
         END DO

      END DO

      ! Save initialization status

      current_maxl = maxl

   END SUBROUTINE create_spherical_harmonics

! **************************************************************************************************
!> \brief   Deallocate the orbital transformation matrices.
!> \date    20.05.2004
!> \author  MK
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE deallocate_spherical_harmonics()

      INTEGER                                            :: l

!$    IF (omp_get_level() > 0) &
!$       CPABORT("deallocate_spherical_harmonics is not thread-safe")

      IF (current_maxl > -1) THEN
         DO l = 0, SIZE(orbtramat, 1) - 1
            DEALLOCATE (orbtramat(l)%c2s)
            DEALLOCATE (orbtramat(l)%s2c)
            DEALLOCATE (orbtramat(l)%slm)
            DEALLOCATE (orbtramat(l)%slm_inv)
         END DO
         DEALLOCATE (orbtramat)
         current_maxl = -1
      END IF

   END SUBROUTINE deallocate_spherical_harmonics

! **************************************************************************************************
!> \brief  Initialize or update the orbital transformation matrices.
!> \param maxl ...
!> \param output_unit ...
!> \date     09.07.1999
!> \par Variables
!>      - maxl   : Maximum angular momentum quantum number
!> \author MK
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE init_spherical_harmonics(maxl, output_unit)

      INTEGER, INTENT(IN)                                :: maxl
      INTEGER                                            :: output_unit

      CHARACTER(LEN=78)                                  :: headline
      INTEGER                                            :: l, nc, ns

!$    IF (omp_get_level() > 0) &
!$       CPABORT("init_spherical_harmonics is not thread-safe")

      IF (maxl < 0) THEN
         CALL cp_abort(__LOCATION__, &
                       "A negative maximum angular momentum quantum "// &
                       "number is invalid")
      END IF

      IF (maxl > current_maxl) THEN

         CALL deallocate_spherical_harmonics()
         CALL create_spherical_harmonics(maxl)

         ! Print the spherical harmonics and the orbital transformation matrices

         IF (output_unit > 0) THEN

            DO l = 0, maxl

               nc = nco(l)
               ns = nso(l)

               headline = "CARTESIAN ORBITAL TO SPHERICAL ORBITAL "// &
                          "TRANSFORMATION MATRIX"
               CALL write_matrix(orbtramat(l)%c2s, l, output_unit, headline)

               headline = "SPHERICAL ORBITAL TO CARTESIAN ORBITAL "// &
                          "TRANSFORMATION MATRIX"
               CALL write_matrix(orbtramat(l)%s2c, l, output_unit, headline)

               headline = "SPHERICAL HARMONICS"
               CALL write_matrix(orbtramat(l)%slm, l, output_unit, headline)

               headline = "INVERSE SPHERICAL HARMONICS"
               CALL write_matrix(orbtramat(l)%slm_inv, l, output_unit, headline)

            END DO

            WRITE (UNIT=output_unit, FMT="(A)") ""

         END IF

      END IF

   END SUBROUTINE init_spherical_harmonics

! **************************************************************************************************
!> \brief   Calculate rotation matrices for spherical harmonics up to value l
!>          Joseph Ivanic and Klaus Ruedenberg
!>          Rotation Matrices for Real Spherical Harmonics. Direct Determination by Recursion
!>          J. Phys. Chem. 1996, 100, 6342-6347
!> \param orbrotmat ...
!> \param rotmat ...
!> \param lval ...
! **************************************************************************************************
   SUBROUTINE calculate_rotmat(orbrotmat, rotmat, lval)
      TYPE(orbrotmat_type), DIMENSION(:), POINTER        :: orbrotmat
      REAL(KIND=dp), DIMENSION(3, 3)                     :: rotmat
      INTEGER, INTENT(IN)                                :: lval

      INTEGER                                            :: l, m1, m2, ns
      REAL(KIND=dp)                                      :: s3, u, uf, v, vf, w, wf
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: r
      REAL(KIND=dp), DIMENSION(-1:1, -1:1)               :: r1
      REAL(KIND=dp), DIMENSION(-2:2, -2:2)               :: r2
      REAL(KIND=dp), DIMENSION(3, 3)                     :: t

      CALL release_rotmat(orbrotmat)

      ALLOCATE (orbrotmat(0:lval))
      DO l = 0, lval
         ns = nso(l)
         ALLOCATE (orbrotmat(l)%mat(ns, ns))
      END DO

      IF (lval >= 0) THEN
         orbrotmat(0)%mat = 1.0_dp
      END IF
      IF (lval >= 1) THEN
         t(1, 1:3) = rotmat(2, 1:3)
         t(2, 1:3) = rotmat(3, 1:3)
         t(3, 1:3) = rotmat(1, 1:3)
         r1(-1:1, -1) = t(1:3, 2)
         r1(-1:1, 0) = t(1:3, 3)
         r1(-1:1, 1) = t(1:3, 1)
         orbrotmat(1)%mat(1:3, 1:3) = r1(-1:1, -1:1)
      END IF
      IF (lval >= 2) THEN
         s3 = SQRT(3.0_dp)
         ! Table 4
         r2(0, 0) = (3.0_dp*r1(0, 0)**2 - 1.0_dp)*0.5_dp
         r2(0, 1) = s3*r1(0, 1)*r1(0, 0)
         r2(0, -1) = s3*r1(0, -1)*r1(0, 0)
         r2(0, -2) = s3*(r1(0, 1)**2 - r1(0, -1)**2)*0.5_dp
         r2(0, -2) = s3*r1(0, 1)*r1(0, -1)
         !
         r2(1, 0) = s3*r1(1, 0)*r1(0, 0)
         r2(1, 1) = r1(1, 1)*r1(0, 0) + r1(1, 0)*r1(0, 1)
         r2(1, -1) = r1(1, -1)*r1(0, 0) + r1(1, 0)*r1(0, -1)
         r2(1, 2) = r1(1, 1)*r1(0, 1) - r1(1, -1)*r1(0, -1)
         r2(1, -2) = r1(1, 1)*r1(0, -1) + r1(1, -1)*r1(0, 1)
         !
         r2(-1, 0) = s3*r1(-1, 0)*r1(0, 0)
         r2(-1, 1) = r1(-1, 1)*r1(0, 0) + r1(-1, 0)*r1(0, 1)
         r2(-1, -1) = r1(-1, -1)*r1(0, 0) + r1(-1, 0)*r1(0, -1)
         r2(-1, 2) = r1(-1, 1)*r1(0, 1) - r1(-1, -1)*r1(0, -1)
         r2(-1, -2) = r1(-1, 1)*r1(0, -1) + r1(-1, -1)*r1(0, 1)
         !
         r2(2, 0) = s3*(r1(1, 0)**2 - r1(-1, 0)**2)*0.5_dp
         r2(2, 1) = r1(1, 1)*r1(1, 0) - r1(-1, 1)*r1(-1, 0)
         r2(2, -1) = r1(1, -1)*r1(1, 0) - r1(-1, -1)*r1(-1, 0)
         r2(2, 2) = (r1(1, 1)**2 - r1(1, -1)**2 - r1(-1, 1)**2 + r1(-1, -1)**2)*0.5_dp
         r2(2, -2) = r1(1, 1)*r1(1, -1) - r1(-1, 1)*r1(-1, -1)
         !
         r2(-2, 0) = s3*r1(1, 0)*r1(-1, 0)
         r2(2, 1) = r1(1, 1)*r1(-1, 0) + r1(1, 0)*r1(-1, 1)
         r2(2, -1) = r1(1, -1)*r1(-1, 0) + r1(1, 0)*r1(-1, -1)
         r2(2, 2) = r1(1, 1)*r1(-1, 1) - r1(1, -1)*r1(-1, -1)
         r2(2, -2) = r1(1, 1)*r1(-1, -1) + r1(1, -1)*r1(-1, 1)
         !
         orbrotmat(2)%mat(1:5, 1:5) = r2(-2:2, -2:2)
      END IF
      IF (lval >= 3) THEN
         ! General recursion
         ALLOCATE (r(0:lval, -lval:lval, -lval:lval))
         r = 0.0_dp
         r(0, 0, 0) = 1.0_dp
         r(1, -1:1, -1:1) = r1(-1:1, -1:1)
         r(2, -2:2, -2:2) = r2(-2:2, -2:2)
         DO l = 3, lval
            DO m1 = -l, l
               DO m2 = -l, l
                  u = u_func(l, m1, m2)
                  v = v_func(l, m1, m2)
                  w = w_func(l, m1, m2)
                  CALL r_recursion(l, m1, m2, r1, r, lval, uf, vf, wf)
                  r(l, m1, m2) = u*uf + v*vf + w*wf
               END DO
            END DO
         END DO
         DO l = 3, lval
            ns = nso(l)
            orbrotmat(l)%mat(1:ns, 1:ns) = r(l, -l:l, -l:l)
         END DO
         DEALLOCATE (r)
      END IF

   END SUBROUTINE calculate_rotmat

! **************************************************************************************************
!> \brief ...
!> \param l ...
!> \param ma ...
!> \param mb ...
!> \return ...
! **************************************************************************************************
   FUNCTION u_func(l, ma, mb) RESULT(u)
      INTEGER                                            :: l, ma, mb
      REAL(KIND=dp)                                      :: u

      IF (ABS(mb) == l) THEN
         u = REAL((l + ma)*(l - ma), KIND=dp)/REAL(2*l*(2*l - 1), KIND=dp)
         u = SQRT(u)
      ELSE IF (ABS(mb) < l) THEN
         u = REAL((l + ma)*(l - ma), KIND=dp)/REAL((l + mb)*(l - mb), KIND=dp)
         u = SQRT(u)
      ELSE
         CPABORT("Illegal m value")
      END IF
   END FUNCTION u_func

! **************************************************************************************************
!> \brief ...
!> \param l ...
!> \param ma ...
!> \param mb ...
!> \return ...
! **************************************************************************************************
   FUNCTION v_func(l, ma, mb) RESULT(v)
      INTEGER                                            :: l, ma, mb
      REAL(KIND=dp)                                      :: v

      INTEGER                                            :: a1, a2, dm0

      dm0 = 0
      IF (ma == 0) dm0 = 1
      IF (ABS(mb) == l) THEN
         a1 = (1 + dm0)*(l + ABS(ma) - 1)*(l + ABS(ma))
         a2 = 2*l*(2*l - 1)
      ELSE IF (ABS(mb) < l) THEN
         a1 = (1 + dm0)*(l + ABS(ma) - 1)*(l + ABS(ma))
         a2 = (l + mb)*(l - mb)
      ELSE
         CPABORT("Illegal m value")
      END IF
      v = 0.5_dp*SQRT(REAL(a1, KIND=dp)/REAL(a2, KIND=dp))*REAL(1 - 2*dm0, KIND=dp)
   END FUNCTION v_func

! **************************************************************************************************
!> \brief ...
!> \param l ...
!> \param ma ...
!> \param mb ...
!> \return ...
! **************************************************************************************************
   FUNCTION w_func(l, ma, mb) RESULT(w)
      INTEGER                                            :: l, ma, mb
      REAL(KIND=dp)                                      :: w

      INTEGER                                            :: a1, a2, dm0

      dm0 = 0
      IF (ma == 0) dm0 = 1
      IF (ABS(mb) == l) THEN
         a1 = (l - ABS(ma) - 1)*(l - ABS(ma))
         a2 = 2*l*(2*l - 1)
      ELSE IF (ABS(mb) < l) THEN
         a1 = (l - ABS(ma) - 1)*(l - ABS(ma))
         a2 = (l + mb)*(l - mb)
      ELSE
         CPABORT("Illegal m value")
      END IF
      w = -0.5_dp*SQRT(REAL(a1, KIND=dp)/REAL(a2, KIND=dp))*REAL(1 - dm0, KIND=dp)
   END FUNCTION w_func

! **************************************************************************************************
!> \brief ...
!> \param i ...
!> \param l ...
!> \param mu ...
!> \param m ...
!> \param r1 ...
!> \param r ...
!> \param lr ...
!> \return ...
! **************************************************************************************************
   FUNCTION p_func(i, l, mu, m, r1, r, lr) RESULT(p)
      INTEGER                                            :: i, l, mu, m
      REAL(KIND=dp), DIMENSION(-1:1, -1:1)               :: r1
      INTEGER                                            :: lr
      REAL(KIND=dp), DIMENSION(0:lr, -lr:lr, -lr:lr)     :: r
      REAL(KIND=dp)                                      :: p

      IF (m == l) THEN
         p = r1(i, 1)*r(l - 1, mu, m) - r1(i, -1)*r(l - 1, mu, -m)
      ELSE IF (m == -l) THEN
         p = r1(i, 1)*r(l - 1, mu, m) + r1(i, -1)*r(l - 1, mu, -m)
      ELSE IF (ABS(m) < l) THEN
         p = r1(i, 0)*r(l - 1, mu, m)
      ELSE
         CPABORT("Illegal m value")
      END IF
   END FUNCTION p_func

! **************************************************************************************************
!> \brief ...
!> \param l ...
!> \param ma ...
!> \param mb ...
!> \param r1 ...
!> \param r ...
!> \param lr ...
!> \param u ...
!> \param v ...
!> \param w ...
! **************************************************************************************************
   SUBROUTINE r_recursion(l, ma, mb, r1, r, lr, u, v, w)
      INTEGER                                            :: l, ma, mb
      REAL(KIND=dp), DIMENSION(-1:1, -1:1)               :: r1
      INTEGER                                            :: lr
      REAL(KIND=dp), DIMENSION(0:lr, -lr:lr, -lr:lr)     :: r
      REAL(KIND=dp)                                      :: u, v, w

      REAL(KIND=dp)                                      :: dd

      IF (ma == 0) THEN
         u = p_func(0, l, 0, mb, r1, r, lr)
         v = p_func(1, l, 1, mb, r1, r, lr) + p_func(-1, l, -1, mb, r1, r, lr)
         w = 0.0_dp
      ELSE IF (ma > 0) THEN
         dd = 1.0_dp
         IF (ma == 1) dd = 1.0_dp
         u = p_func(0, l, ma, mb, r1, r, lr)
         v = p_func(1, l, ma - 1, mb, r1, r, lr)*SQRT(1.0_dp + dd) - p_func(-1, l, -ma + 1, mb, r1, r, lr)*(1.0_dp - dd)
         w = p_func(1, l, ma + 1, mb, r1, r, lr) + p_func(-1, l, -1 - ma, mb, r1, r, lr)
      ELSE IF (ma < 0) THEN
         dd = 1.0_dp
         IF (ma == -1) dd = 1.0_dp
         u = p_func(0, l, ma, mb, r1, r, lr)
         v = p_func(1, l, ma + 1, mb, r1, r, lr)*(1.0_dp - dd) + p_func(-1, l, -ma - 1, mb, r1, r, lr)*SQRT(1.0_dp + dd)
         w = p_func(1, l, ma - 1, mb, r1, r, lr) - p_func(-1, l, -ma + 1, mb, r1, r, lr)
      END IF
   END SUBROUTINE

! **************************************************************************************************
!> \brief   Release rotation matrices
!> \param orbrotmat ...
! **************************************************************************************************
   SUBROUTINE release_rotmat(orbrotmat)
      TYPE(orbrotmat_type), DIMENSION(:), POINTER        :: orbrotmat

      INTEGER                                            :: i

      IF (ASSOCIATED(orbrotmat)) THEN
         DO i = LBOUND(orbrotmat, 1), UBOUND(orbrotmat, 1)
            IF (ASSOCIATED(orbrotmat(i)%mat)) DEALLOCATE (orbrotmat(i)%mat)
         END DO
         DEALLOCATE (orbrotmat)
      END IF

   END SUBROUTINE release_rotmat

! **************************************************************************************************
!> \brief   Print a matrix to the logical unit "lunit".
!> \param matrix ...
!> \param l ...
!> \param lunit ...
!> \param headline ...
!> \date    07.06.2000
!> \par Variables
!>      - matrix  : Matrix to be printed.
!>      - l       : Angular momentum quantum number
!>      - lunit   : Logical unit number.
!>      - headline: Headline of the matrix.
!> \author  MK
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE write_matrix(matrix, l, lunit, headline)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: matrix
      INTEGER, INTENT(IN)                                :: l, lunit
      CHARACTER(LEN=*), INTENT(IN)                       :: headline

      CHARACTER(12)                                      :: symbol
      CHARACTER(LEN=78)                                  :: string
      INTEGER                                            :: from, i, ic, is, jc, lx, ly, lz, m, nc, &
                                                            to

      ! Write headline

      WRITE (UNIT=lunit, FMT="(/,/,T2,A)") TRIM(headline)

      ! Write the matrix in the defined format

      nc = nco(l)

      DO ic = 1, nc, 3
         from = ic
         to = MIN(nc, from + 2)
         i = 1
         DO lx = l, 0, -1
            DO ly = l - lx, 0, -1
               lz = l - lx - ly
               jc = co(lx, ly, lz)
               IF ((jc >= from) .AND. (jc <= to)) THEN
                  symbol = cgf_symbol(1, (/lx, ly, lz/))
                  WRITE (UNIT=string(i:), FMT="(A18)") TRIM(symbol(3:12))
                  i = i + 18
               END IF
            END DO
         END DO
         WRITE (UNIT=lunit, FMT="(/,T8,A)") TRIM(string)
         symbol = ""
         DO m = -l, l
            is = l + m + 1
            symbol = sgf_symbol(1, l, m)
            WRITE (UNIT=lunit, FMT="(T4,A4,3(1X,F17.12))") &
               symbol(3:6), (matrix(is, jc), jc=from, to)
         END DO
      END DO

   END SUBROUTINE write_matrix

END MODULE orbital_transformation_matrices
