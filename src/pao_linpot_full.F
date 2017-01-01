!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Full parametrization of Fock matrix, ie. the identity parametrization.
!> \author Ole Schuett
! **************************************************************************************************
MODULE pao_linpot_full
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE kinds,                           ONLY: dp
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_linpot_full'

   PUBLIC :: linpot_full_count_terms, linpot_full_calc_terms

CONTAINS

! **************************************************************************************************
!> \brief Count number of terms for given atomic kind
!> \param qs_env ...
!> \param ikind ...
!> \param nterms ...
! **************************************************************************************************
   SUBROUTINE linpot_full_count_terms(qs_env, ikind, nterms)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: ikind
      INTEGER, INTENT(OUT)                               :: nterms

      INTEGER                                            :: n
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set)
      CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis_set)
      n = basis_set%nsgf
      nterms = n+n*(n-1)/2

   END SUBROUTINE linpot_full_count_terms

! **************************************************************************************************
!> \brief Builds potential terms
!> \param V_blocks ...
! **************************************************************************************************
   SUBROUTINE linpot_full_calc_terms(V_blocks)
      REAL(dp), DIMENSION(:, :, :), INTENT(OUT)          :: V_blocks

      INTEGER                                            :: i, j, kterm, n, nterms

      N = SIZE(V_blocks, 1)
      CPASSERT(SIZE(V_blocks, 2) == N)
      nterms = SIZE(V_blocks, 3)

      V_blocks = 0.0_dp
      kterm = 0
      DO i = 1, N
      DO j = i, N
         kterm = kterm+1
         V_blocks(i, j, kterm) = 1.0_dp
         V_blocks(j, i, kterm) = 1.0_dp
      ENDDO
      ENDDO

      CPASSERT(kterm == nterms)
   END SUBROUTINE linpot_full_calc_terms

END MODULE pao_linpot_full
