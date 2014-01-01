!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Sets the typo for the gaussian treatment of the qm/mm interaction
!> \par History
!>      06.2004 created
!> \author Teodoro Laino
! *****************************************************************************
MODULE qmmm_gaussian_types
  USE f77_blas
  USE kinds,                           ONLY: dp
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

! *****************************************************************************
!> \par History
!>      06.2004 created [tlaino]
!> \author Laino Teodoro
! *****************************************************************************
  TYPE qmmm_gaussian_type
     INTEGER                               :: number_of_gaussians
     REAL(KIND=dp)                         :: Elp_Radius, Elp_Radius_Corr
     INTEGER, DIMENSION(:), POINTER        :: grid_level
     REAL(KIND=dp), DIMENSION(:), POINTER  :: Ak, Gk
  END TYPE qmmm_gaussian_type

! *****************************************************************************
!> \brief represent a pointer to a qmmm_gaussian_type, to be able to create arrays
!>      of pointers
!> \param qmmm_gaussian_type the pointer to the qmmm_gaussian_type
!> \par History
!>      06.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  TYPE qmmm_gaussian_p_type
     TYPE(qmmm_gaussian_type), POINTER :: pgf
  END TYPE qmmm_gaussian_p_type

  PUBLIC  :: qmmm_gaussian_type,&
             qmmm_gaussian_p_type

END MODULE qmmm_gaussian_types

