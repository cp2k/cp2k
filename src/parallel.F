!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

MODULE parallel
  
  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY : dbl
  USE message_passing, ONLY : mp_start, mp_environ, mp_end
  
  PRIVATE
  PUBLIC :: start_parallel, end_parallel
  
CONTAINS

!******************************************************************************

SUBROUTINE start_parallel ( globenv )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( global_environment_type ) :: globenv
  
!------------------------------------------------------------------------------
  
  CALL mp_start()
  CALL mp_environ(globenv % num_pe,globenv % mepos,globenv % group)
  globenv % source = 0
  globenv % ionode = (globenv % mepos==globenv % source)
  
END SUBROUTINE start_parallel

!******************************************************************************

SUBROUTINE end_parallel()
  
  IMPLICIT NONE
  
!------------------------------------------------------------------------------
  
  CALL mp_end()
  
END SUBROUTINE end_parallel

!******************************************************************************

END MODULE parallel
