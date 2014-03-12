!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief A common interface for passing a callback into the md_run loop.
!> \par History
!> \author Ole
! *****************************************************************************
MODULE mdctrl_types

  USE glbopt_types,                    ONLY: glbopt_mdctrl_data_type
#include "./common/cp_common_uses.f90"

 IMPLICIT NONE
 PRIVATE

 TYPE mdctrl_type
    TYPE(glbopt_mdctrl_data_type), POINTER                 :: glbopt => Null()
    !... and possible more in the future
 END TYPE mdctrl_type


 PUBLIC :: mdctrl_type

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mdctrl_types' 

END MODULE mdctrl_types

