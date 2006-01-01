!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


!!****h* cp2k/qmmm_gaussian_types [1.0] *
!!
!!   NAME
!!     qmmm_gaussian_types
!!
!!   FUNCTION
!!     Sets the typo for the gaussian treatment of the qm/mm interaction
!!
!!   NOTES
!!     
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     06.2004 created
!!
!!   SOURCE
!****************************************************************************
MODULE qmmm_gaussian_types
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
!!****s* qmmm_gaussian_types/qmmm_gaussian_type [1.0] *
!!
!!   NAME
!!     qmmm_gaussian_types
!!
!!   FUNCTION
!!     
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     -
!!
!!   AUTHOR
!!     Laino Teodoro
!!
!!   MODIFICATION HISTORY
!!     06.2004 created [tlaino]
!!
!!   SOURCE
  !***************************************************************************

  TYPE qmmm_gaussian_type
     INTEGER                               :: number_of_gaussians
     REAL(KIND=dp)                         :: Elp_Radius, Elp_Radius_Corr
     INTEGER, DIMENSION(:), pointer        :: grid_level
     REAL(KIND=dp), DIMENSION(:), POINTER  :: Ak, Gk
  END TYPE qmmm_gaussian_type

!!****s* qmmm_gaussian_types/qmmm_gaussian_p_type [1.0] *
!!
!!   NAME
!!     qmmm_gaussian_p_type
!!
!!   FUNCTION
!!     represent a pointer to a qmmm_gaussian_type, to be able to create arrays
!!     of pointers
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - qmmm_gaussian_type: the pointer to the qmmm_gaussian_type
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     06.2004 created [tlaino]
!!
!!   SOURCE
  !***************************************************************************

  TYPE qmmm_gaussian_p_type
     TYPE(qmmm_gaussian_type), POINTER :: pgf
  END TYPE qmmm_gaussian_p_type

  PUBLIC  :: qmmm_gaussian_type,&
             qmmm_gaussian_p_type

END MODULE qmmm_gaussian_types

