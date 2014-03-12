!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Common types used by the global optimization schemes.
!> \author Ole Schuett
! *****************************************************************************
MODULE glbopt_types

  USE kinds,                           ONLY: dp
#include "./common/cp_common_uses.f90"

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


