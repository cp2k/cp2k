!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utility functions for CNEO-DFT
!>      (see J. Chem. Theory Comput. 2025, 21, 16, 7865–7877)
!> \par History
!>      08.2025 created [zc62]
!> \author Zehua Chen
! **************************************************************************************************
MODULE qs_cneo_utils
   USE ao_util,                         ONLY: trace_r_AxB
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE kinds,                           ONLY: dp
   USE memory_utilities,                ONLY: reallocate
   USE qs_harmonics_atom,               ONLY: get_none0_cg_list,&
                                              harmonics_atom_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ! *** Global parameters ***
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_cneo_utils'

   ! *** Public subroutines ***
   PUBLIC :: atom_solve_cneo, cneo_gather, cneo_scatter, create_harmonics_atom_cneo, &
             get_maxl_CG_cneo

CONTAINS

! **************************************************************************************************
!> \brief Mostly copied from qs_harmonics_atom::create_harmonics_atom
!> \param harmonics ...
!> \param my_CG ...
!> \param llmax ...
!> \param maxs ...
!> \param max_s_harm ...
! **************************************************************************************************
   SUBROUTINE create_harmonics_atom_cneo(harmonics, my_CG, llmax, maxs, max_s_harm)

      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: my_CG
      INTEGER, INTENT(IN)                                :: llmax, maxs, max_s_harm

      INTEGER                                            :: i, is

      CPASSERT(ASSOCIATED(harmonics))

      harmonics%max_s_harm = max_s_harm
      harmonics%llmax = llmax

      NULLIFY (harmonics%my_CG, harmonics%my_CG_dxyz, harmonics%my_CG_dxyz_asym)
      CALL reallocate(harmonics%my_CG, 1, maxs, 1, maxs, 1, max_s_harm)

      DO i = 1, max_s_harm
         DO is = 1, maxs
            harmonics%my_CG(1:maxs, is, i) = my_CG(1:maxs, is, i)
         END DO
      END DO

   END SUBROUTINE create_harmonics_atom_cneo

! **************************************************************************************************
!> \brief Mostly copied from qs_harmonics_atom::get_maxl_CG
!> \param harmonics ...
!> \param orb_basis ...
!> \param llmax ...
!> \param max_s_harm ...
! **************************************************************************************************
   SUBROUTINE get_maxl_CG_cneo(harmonics, orb_basis, llmax, max_s_harm)

      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis
      INTEGER, INTENT(IN)                                :: llmax, max_s_harm

      INTEGER                                            :: is1, is2, itmp, max_iso_not0, nset
      INTEGER, DIMENSION(:), POINTER                     :: lmax, lmin

      CPASSERT(ASSOCIATED(harmonics))

      CALL get_gto_basis_set(gto_basis_set=orb_basis, lmax=lmax, lmin=lmin, nset=nset)

      !   *** Assign indices for the non null CG coefficients ***
      max_iso_not0 = 0
      DO is1 = 1, nset
         DO is2 = 1, nset
            CALL get_none0_cg_list(harmonics%my_CG, &
                                   lmin(is1), lmax(is1), lmin(is2), lmax(is2), &
                                   max_s_harm, llmax, max_iso_not0=itmp)
            max_iso_not0 = MAX(max_iso_not0, itmp)
         END DO ! is2
      END DO ! is1
      harmonics%max_iso_not0 = max_iso_not0

   END SUBROUTINE get_maxl_CG_cneo

! **************************************************************************************************
!> \brief Mostly copied from atom_utils::atom_solve
!> \param hmat ...
!> \param f ...
!> \param umat ...
!> \param orb ...
!> \param ener ...
!> \param pmat ...
!> \param r ...
!> \param dist ...
!> \param nb ...
!> \param nv ...
! **************************************************************************************************
   SUBROUTINE atom_solve_cneo(hmat, f, umat, orb, ener, pmat, r, dist, nb, nv)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: hmat
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: f
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: umat
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: orb
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: ener
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: pmat
      REAL(KIND=dp), DIMENSION(3), INTENT(INOUT)         :: r
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: dist
      INTEGER, INTENT(IN)                                :: nb, nv

      INTEGER                                            :: info, lwork, m, n
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: w, work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: a, b, h_fx

      CPASSERT(nb >= nv)

      orb = 0._dp
      n = nb
      m = nv
      IF (n > 0 .AND. m > 0) THEN
         lwork = MAX(n*n, n + 100)
         ALLOCATE (a(m, m), b(n, m), w(m), work(lwork))
         IF (DOT_PRODUCT(f, f) /= 0.0_dp) THEN
            ALLOCATE (h_fx(n, n))
            h_fx(1:n, 1:n) = hmat(1:n, 1:n) + f(1)*dist(1:n, 1:n, 1) + &
                             f(2)*dist(1:n, 1:n, 2) + f(3)*dist(1:n, 1:n, 3)
            CALL dgemm("N", "N", n, m, n, 1.0_dp, h_fx, n, umat, n, 0.0_dp, b, n)
            DEALLOCATE (h_fx)
         ELSE
            CALL dgemm("N", "N", n, m, n, 1.0_dp, hmat, n, umat, n, 0.0_dp, b, n)
         END IF
         CALL dgemm("T", "N", m, m, n, 1.0_dp, umat, n, b, n, 0.0_dp, a, m)
         CALL dsyev("V", "U", m, a, m, w, work, lwork, info)
         CALL dgemm("N", "N", n, m, m, 1.0_dp, umat, n, a, m, 0.0_dp, b, n)

         m = MIN(m, SIZE(orb, 2))
         orb(1:n, 1:m) = b(1:n, 1:m)
         ener(1:m) = w(1:m)

         DEALLOCATE (a, b, w, work)

         ! calculate the density matrix using the orbital with the lowest orbital energy
         pmat = 0.0_dp
         CALL dger(n, n, 1.0_dp, orb(:, 1), 1, orb(:, 1), 1, pmat, n)
         ! calculate the expectation position (basis center as the origin)
         r = [trace_r_AxB(dist(1:n, 1:n, 1), n, pmat, n, n, n), &
              trace_r_AxB(dist(1:n, 1:n, 2), n, pmat, n, n, n), &
              trace_r_AxB(dist(1:n, 1:n, 3), n, pmat, n, n, n)]
      END IF

   END SUBROUTINE atom_solve_cneo

! **************************************************************************************************
!> \brief Mostly copied from qs_oce_methods::prj_gather
!> \param ain ...
!> \param aout ...
!> \param nbas ...
!> \param n2oindex ...
! **************************************************************************************************
   SUBROUTINE cneo_gather(ain, aout, nbas, n2oindex)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: ain
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: aout
      INTEGER, INTENT(IN)                                :: nbas
      INTEGER, DIMENSION(:), POINTER                     :: n2oindex

      INTEGER                                            :: i, ip, j, jp

      DO i = 1, nbas
         ip = n2oindex(i)
         DO j = 1, nbas
            jp = n2oindex(j)
            aout(j, i) = ain(jp, ip)
         END DO
      END DO

   END SUBROUTINE cneo_gather

! **************************************************************************************************
!> \brief Mostly copied from qs_oce_methods::prj_scatter
!> \param ain ...
!> \param aout ...
!> \param nbas ...
!> \param n2oindex ...
! **************************************************************************************************
   SUBROUTINE cneo_scatter(ain, aout, nbas, n2oindex)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: ain
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: aout
      INTEGER, INTENT(IN)                                :: nbas
      INTEGER, DIMENSION(:), POINTER                     :: n2oindex

      INTEGER                                            :: i, ip, j, jp

      DO i = 1, nbas
         ip = n2oindex(i)
         DO j = 1, nbas
            jp = n2oindex(j)
            aout(jp, ip) = aout(jp, ip) + ain(j, i)
         END DO
      END DO

   END SUBROUTINE cneo_scatter

END MODULE qs_cneo_utils
