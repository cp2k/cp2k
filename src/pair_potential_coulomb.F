!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

MODULE pair_potential_coulomb

  USE erf_fn,                          ONLY: erf,&
                                             erfc
  USE input_constants,                 ONLY: do_ewald_none
  USE kinds,                           ONLY: dp
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pair_potential_coulomb'

  PUBLIC :: potential_coulomb

CONTAINS

! *****************************************************************************
!> \brief Evaluates the electrostatic energy and force
!> \author Toon.Verstraelen@gmail.com
! *****************************************************************************
  FUNCTION potential_coulomb(r2, fscalar, qfac, ewald_type, alpha, mm_radius, &
      interaction_cutoff)

    REAL(KIND=dp), INTENT(IN)                :: r2
    REAL(KIND=dp), INTENT(INOUT)             :: fscalar
    REAL(KIND=dp), INTENT(IN)                :: qfac
    INTEGER, INTENT(IN)                      :: ewald_type
    REAL(KIND=dp), INTENT(IN)                :: alpha, mm_radius, &
                                                interaction_cutoff
    REAL(KIND=dp)                            :: potential_coulomb

    REAL(KIND=dp), PARAMETER :: two_over_sqrt_pi = 1.128379167095512573

    REAL(KIND=dp)                            :: r, x1, x2

    r = SQRT(r2)
    IF (mm_radius > 0.0_dp) THEN
       IF (ewald_type == do_ewald_none) THEN
          x2 = r/mm_radius
          potential_coulomb = erf(x2)/r
          fscalar = fscalar + qfac*(potential_coulomb - &
            two_over_sqrt_pi*EXP(-x2*x2)/mm_radius)/r2
       ELSE
          x1 = alpha*r
          x2 = r/mm_radius
          potential_coulomb = (erf(x2) - erf(x1))/r
          fscalar = fscalar + qfac*(potential_coulomb + &
            two_over_sqrt_pi*(EXP(-x1*x1)*alpha-EXP(-x2*x2)/mm_radius))/r2
       END IF
    ELSE
       IF (ewald_type == do_ewald_none) THEN
          potential_coulomb = 1.0/r
          fscalar = fscalar + qfac*potential_coulomb/r2
       ELSE
          x1 = alpha*r
          potential_coulomb = erfc(x1)/r
          fscalar = fscalar + qfac*(potential_coulomb + &
            two_over_sqrt_pi*EXP(-x1*x1)*alpha)/r2
       END IF
    END IF
    potential_coulomb = qfac*(potential_coulomb - interaction_cutoff)
  END FUNCTION potential_coulomb

END MODULE pair_potential_coulomb

