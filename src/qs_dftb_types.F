!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Definition of the DFTB parameter types.
!> \author JGH (24.02.2007)
! **************************************************************************************************
MODULE qs_dftb_types

   USE kinds,                           ONLY: default_string_length,&
                                              dp
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_dftb_types'

! **************************************************************************************************
   TYPE qs_dftb_atom_type
      ! PRIVATE
      CHARACTER(LEN=default_string_length)   :: typ = ""
      CHARACTER(LEN=default_string_length)   :: name = ""
      LOGICAL                                :: defined = .FALSE.
      INTEGER                                :: z = -1 !atomic number
      REAL(KIND=dp)                          :: zeff = -1.0_dp !effective core charge
      INTEGER                                :: natorb = -1 !number of orbitals
      INTEGER                                :: lmax = -1 !max angular momentum
      REAL(KIND=dp), DIMENSION(0:3)          :: skself = -1.0_dp !orbital energy
      REAL(KIND=dp), DIMENSION(0:3)          :: occupation = -1.0_dp !free atom occupation
      REAL(KIND=dp), DIMENSION(0:3)          :: eta = -1.0_dp !orbital hardness
      REAL(KIND=dp)                          :: energy = -1.0_dp !free atom energy
      REAL(KIND=dp)                          :: cutoff = -1.0_dp !cutoff radius for f matrix
      REAL(KIND=dp)                          :: xi = -1.0_dp, di = -1.0_dp !London parameter
      REAL(KIND=dp)                          :: rcdisp = -1.0_dp !cutoff radius for vdW
      REAL(KIND=dp)                          :: dudq = -1.0_dp !DFTB3 hardness derivative
   END TYPE qs_dftb_atom_type

! **************************************************************************************************
   TYPE qs_dftb_pairpot_type
      REAL(KIND=dp)                          :: dgrd = -1.0_dp ! grid spacing
      INTEGER                                :: ngrd = -1 ! number of grid points
      INTEGER                                :: ngrdcut = -1 ! grid cutoff
      INTEGER                                :: llm = -1 ! number of interactions (l,l,m)
      INTEGER                                :: n_urpoly = -1 ! order of polynomial
      REAL(KIND=dp)                          :: urep_cut = -1.0_dp ! cutoff for repulsive pot.
      REAL(KIND=dp), DIMENSION(10)           :: urep = -1.0_dp ! coefficients for repulsive
      ! potential in polynomial form
      INTEGER                                :: spdim = -1 ! number of points for
      ! spline representation
      REAL(KIND=dp)                          :: s_cut = -1.0_dp ! left-hand cutoff
      REAL(KIND=dp), DIMENSION(3)            :: srep = -1.0_dp ! coefficients for extrapolation
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: spxr => NULL() ! spline points
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: scoeff => NULL() ! spline coefficients
      REAL(KIND=dp), DIMENSION(2)            :: surr = -1.0_dp ! coefficients for last point
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: fmat => NULL() ! Slater-Koster table (Hamiltonian)
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: smat => NULL() ! Slater-Koster table (overlap)
      ! van der Waals parameter
      REAL(KIND=dp)                          :: xij = -1.0_dp, dij = -1.0_dp ! standard LJ parameters
      REAL(KIND=dp)                          :: x0ij = -1.0_dp ! Evdw(x0) = 0
      REAL(KIND=dp)                          :: a = -1.0_dp, b = -1.0_dp, c = -1.0_dp ! Short range polynomial coeffs
   END TYPE qs_dftb_pairpot_type

! *** Public data types ***

   PUBLIC :: qs_dftb_atom_type, qs_dftb_pairpot_type, &
             qs_dftb_pairpot_init, qs_dftb_pairpot_create, qs_dftb_pairpot_release

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param pairpot ...
! **************************************************************************************************
   SUBROUTINE qs_dftb_pairpot_init(pairpot)
      TYPE(qs_dftb_pairpot_type), DIMENSION(:, :), &
         POINTER                                         :: pairpot

      INTEGER                                            :: i, j

      IF (ASSOCIATED(pairpot)) THEN
         DO i = 1, SIZE(pairpot, 1)
            DO j = 1, SIZE(pairpot, 2)
               NULLIFY (pairpot(i, j)%spxr, pairpot(i, j)%scoeff, &
                        pairpot(i, j)%smat, pairpot(i, j)%fmat)
            END DO
         END DO
      END IF

   END SUBROUTINE qs_dftb_pairpot_init

! **************************************************************************************************
!> \brief ...
!> \param pairpot ...
!> \param ngrd ...
!> \param llm ...
!> \param spdim ...
! **************************************************************************************************
   SUBROUTINE qs_dftb_pairpot_create(pairpot, ngrd, llm, spdim)
      TYPE(qs_dftb_pairpot_type)                         :: pairpot
      INTEGER, INTENT(IN)                                :: ngrd, llm, spdim

      pairpot%ngrd = ngrd
      pairpot%spdim = spdim
      pairpot%llm = llm

      IF (spdim > 0) THEN
         ALLOCATE (pairpot%spxr(spdim, 2))

         ALLOCATE (pairpot%scoeff(spdim, 4))
      END IF

      ALLOCATE (pairpot%fmat(ngrd, llm))

      ALLOCATE (pairpot%smat(ngrd, llm))

   END SUBROUTINE qs_dftb_pairpot_create

! **************************************************************************************************
!> \brief ...
!> \param pairpot ...
! **************************************************************************************************
   SUBROUTINE qs_dftb_pairpot_release(pairpot)
      TYPE(qs_dftb_pairpot_type), DIMENSION(:, :), &
         POINTER                                         :: pairpot

      INTEGER                                            :: i, j, n1, n2

      IF (ASSOCIATED(pairpot)) THEN
         n1 = SIZE(pairpot, 1)
         n2 = SIZE(pairpot, 2)
         DO i = 1, n1
            DO j = 1, n2
               IF (ASSOCIATED(pairpot(i, j)%spxr)) THEN
                  DEALLOCATE (pairpot(i, j)%spxr)
               END IF
               IF (ASSOCIATED(pairpot(i, j)%scoeff)) THEN
                  DEALLOCATE (pairpot(i, j)%scoeff)
               END IF
               IF (ASSOCIATED(pairpot(i, j)%smat)) THEN
                  DEALLOCATE (pairpot(i, j)%smat)
               END IF
               IF (ASSOCIATED(pairpot(i, j)%fmat)) THEN
                  DEALLOCATE (pairpot(i, j)%fmat)
               END IF
            END DO
         END DO
         DEALLOCATE (pairpot)
      END IF

   END SUBROUTINE qs_dftb_pairpot_release

END MODULE qs_dftb_types

