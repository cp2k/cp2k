!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/kg_main [1.0] *
!!
!!   NAME
!!     kg_main
!!
!!   FUNCTION
!!     perform kim-gordon moecular dynamics 
!!   AUTHOR
!!     gt SEPT-23-2002
!!   SOURCE
!******************************************************************************

MODULE kg_main 

  USE global_types,                    ONLY: global_environment_type
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: kg
  
!!*****
!-----------------------------------------------------------------------------!

CONTAINS

!-----------------------------------------------------------------------------!
! KG KG KG                                                                    !
!-----------------------------------------------------------------------------!
!!****** kg_main/kg [1.0] *
!!
!!   NAME
!!     kg
!!
!!   SYNOPSIS
!!     Subroutine kg(globenv)
!!       Implicit None
!!       Type(global_environment_type), Intent (INOUT):: globenv
!!     End Subroutine kg
!!
!!   FUNCTION
!!     Controls program flow for Kim-Gordon type MD
!!
!!   AUTHOR
!!     GT-CJM
!!
!!   MODIFICATION HISTORY
!!
!!   USED BY
!!     cp2k
!!
!!*** *************************************************************************

SUBROUTINE kg ( globenv )

  IMPLICIT NONE

! Argument
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv



END SUBROUTINE kg

!******************************************************************************

END MODULE kg_main
!******************************************************************************

!******************************************************************************
