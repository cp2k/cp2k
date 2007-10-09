!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!***** cp2k/csvr_system_types [1.0] *
!!
!!   NAME
!!     csvr_system_types
!!
!!   FUNCTION
!!     Type for the canonical sampling through velocity rescaling     
!!
!!   AUTHOR
!!     Teodoro Laino - 09.2007 University of Zurich [tlaino]
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE csvr_system_types
  USE kinds,                           ONLY: dp

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: csvr_system_type

  TYPE csvr_system_type
     INTEGER       :: region
     REAL(KIND=dp) :: tau_csvr
     REAL(KIND=dp) :: thermostat_energy
  END TYPE csvr_system_type

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'csvr_system_types'

!!*****
!******************************************************************************
END MODULE csvr_system_types
!******************************************************************************
