!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/pol_thk [1.0] *
!!
!!   NAME
!!     pol_thk
!!
!!   FUNCTION
!!     Calculation of the first order contribution of the Hohenbeg-Kohn
!!     functional t
!!
!!     E (hk) =   c_i* thk_i
!!
!!   AUTHOR
!!     gloria (06 Jul 2001)     
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE pol_thk

! *****************************************************************************

  USE ao_types,                        ONLY: ao_type
  USE empirical_parameters,            ONLY: empirical_parameter_type
  USE kinds,                           ONLY: dbl

  IMPLICIT NONE

  PRIVATE
 
  PUBLIC :: force_thk

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE force_thk(ao, empparm, fthk, energy)

!   ***************************************************************************

    TYPE(ao_type), INTENT(INOUT)             :: ao
    TYPE(empirical_parameter_type), &
      DIMENSION(:), INTENT(IN)               :: empparm
    REAL(dbl), DIMENSION(:), INTENT(INOUT)   :: fthk
    REAL(dbl), INTENT(OUT)                   :: energy

    INTEGER                                  :: icoef, ikind, iset, ncoef
    REAL(dbl)                                :: thk

  ncoef = size(ao%cr)
  fthk(:) = 0.0_dbl

  DO icoef = 1, ncoef
     ikind = ao % coef_to_basis (icoef)
     iset = ao % coef_to_set (icoef)
     thk = empparm (ikind) % hk_param (iset)
     fthk ( icoef ) = - thk  
     energy = energy + ao % cr (icoef) * thk
  END DO

 END SUBROUTINE force_thk

! *****************************************************************************

END MODULE pol_thk

! *****************************************************************************
