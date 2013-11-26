!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2013  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief A common interface for passing a callback into the md_run loop.
!> \par History
!> \author Ole
! *****************************************************************************
MODULE mdctrl_types


#include "cp_common_uses.h"
 IMPLICIT NONE
 PRIVATE

 TYPE mdctrl_type
    ! Nothing here, yet.
    !... and possible more in the future
 END TYPE mdctrl_type


 PUBLIC :: mdctrl_type

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mdctrl_types' 

END MODULE mdctrl_types

