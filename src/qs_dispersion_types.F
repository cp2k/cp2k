!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Definition of disperson types for DFT calculations
!> \author JGH (20.10.2008)
! *****************************************************************************
MODULE qs_dispersion_types

  USE f77_blas
  USE kinds,                           ONLY: default_string_length,&
                                             dp
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_dispersion_types'

! *****************************************************************************
  INTEGER, PARAMETER                       :: grimme_pp = 1

  TYPE qs_dispersion_type
    INTEGER                                :: TYPE
    CHARACTER (LEN=default_string_length)  :: pp_type
    REAL(KIND=dp)                          :: scaling
    REAL(KIND=dp)                          :: rc_disp
    REAL(KIND=dp)                          :: exp_pre
  END TYPE qs_dispersion_type

  TYPE qs_atom_dispersion_type
    INTEGER                                :: TYPE
    LOGICAL                                :: defined
    REAL(KIND=dp)                          :: vdw_radii  !van der Waals radii
    REAL(KIND=dp)                          :: c6         !c6 coefficients
  END TYPE qs_atom_dispersion_type

! *****************************************************************************

  PUBLIC :: qs_atom_dispersion_type, qs_dispersion_type
  PUBLIC :: grimme_pp

! *****************************************************************************

END MODULE qs_dispersion_types

