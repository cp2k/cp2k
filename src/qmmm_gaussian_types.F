#include "cp_prep_globals.h"
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
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dbl

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

  TYPE Info_Pot_Type
     REAL(KIND=dbl)  :: Rmax, Rmin
     REAL(KIND=dbl)  :: MaxErr, ChiSq     
  END TYPE Info_Pot_Type

  TYPE LongRange_Type
     REAL(KIND=dbl), DIMENSION(:), POINTER :: Pot, DPot
     REAL(KIND=dbl)  :: Rmax, Rmin, dx
     INTEGER         :: npts
  END TYPE LongRange_Type

  TYPE qmmm_gaussian_type
     INTEGER :: number_of_gaussians
     REAL(KIND=dbl)  :: Elp_Radius
     TYPE(Info_Pot_Type)  :: Info
     TYPE(LongRange_Type) :: LR
     REAL(KIND=dbl)  :: A0 
     REAL(KIND=dbl), DIMENSION(:), POINTER :: Ak, Gk
     INTEGER, DIMENSION(:), pointer :: grid_level
     INTEGER,        DIMENSION(:), POINTER :: mm_atom_index
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
             qmmm_gaussian_p_type,&
             Info_Pot_type,&
             LongRange_Type

END MODULE qmmm_gaussian_types

