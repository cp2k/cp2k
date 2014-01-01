!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief performs global geometry optimization
!> \par History
!> \author Ole
! *****************************************************************************
MODULE glbopt_types

  USE kinds,                           ONLY: dp
#include "cp_common_uses.h"

 IMPLICIT NONE
 PRIVATE


 TYPE glbopt_mdctrl_data_type
   INTEGER                                    :: md_bump_counter
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE   :: epot_history
   INTEGER                                    :: output_unit
   INTEGER                                    :: itimes
   INTEGER                                    :: bump_steps_upwards
   INTEGER                                    :: bump_steps_downwards
   INTEGER                                    :: md_bumps_max
 END TYPE glbopt_mdctrl_data_type



 PUBLIC :: glbopt_mdctrl_data_type

END MODULE glbopt_types


