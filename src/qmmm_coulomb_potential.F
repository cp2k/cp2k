!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002,2003,2004  CP2K developers group                       !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/qmmm_coulomb_potential [1.0] *
!!
!!   NAME
!!     qmmm_coulomb_potential
!!
!!   FUNCTION
!!     A class of general function to treat the electrostatic 
!!     Coulomb potential into a  QM/MM electrostatic coupling
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!
!!
!!   SOURCE
!****************************************************************************

MODULE qmmm_coulomb_potential
  USE kinds,                           ONLY: dp,&
                                             dbl
  Implicit None
  
  PRIVATE
  PUBLIC :: Vpot
  
  INTERFACE Vpot
     MODULE PROCEDURE Pot1
  END INTERFACE
  
CONTAINS
!!****f* qmmm_coulomb_potential/Pot1 [1.0] *
!!
!!   NAME
!!     Pot1
!!
!!   FUNCTION
!!     Computes the Coulomb Potential using the smooth function as
!!     Optimized in the article:
!!
!!     A. Laio, J. VandeVondele and U. Rothlisberger 
!!     J. Chem. Phys. 116(16), 6941-6947 (2002)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - r0 : core radius
!!     - r  : distance
!!     - Nexp : exponent of the smooth function
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     05.2004 created [tlaino]
!!
!!*** **********************************************************************
  FUNCTION Pot1(Nexp,r,r0) RESULT(vsmooth)
    REAL(dbl), INTENT(IN)   :: r, r0
    INTEGER,   INTENT(IN)   :: Nexp
    REAL(dbl) :: vsmooth
    
    vsmooth = (r**Nexp - r0**Nexp)/(r**(Nexp+1) - r0**(Nexp+1))
    
  End FUNCTION Pot1
End MODULE qmmm_coulomb_potential
