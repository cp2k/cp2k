!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/parallel [1.0] *
!!
!!   NAME
!!     parallel
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE parallel

  USE global_types, ONLY : global_environment_type, print_key_type
  USE kinds, ONLY : dbl
  USE message_passing, ONLY : mp_start, mp_world_setup, mp_end
  USE cp_para_env, ONLY: cp_para_env_create
  use cp_log_handling, only: cp_init_default_logger

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: start_parallel, end_parallel, cp2k_is_parallel

#ifdef __parallel
  logical, parameter :: cp2k_is_parallel=.true.
#else
  logical, parameter :: cp2k_is_parallel=.false.
#endif

!******************************************************************************

CONTAINS

!******************************************************************************

SUBROUTINE start_parallel ( globenv )

! Arguments
  TYPE ( global_environment_type ), TARGET :: globenv

!------------------------------------------------------------------------------

  CALL mp_start ( )
  CALL mp_world_setup ( globenv % num_pe, globenv % mepos, globenv % group )
  globenv % source = 0
  globenv % ionode = (globenv % mepos==globenv % source)
  nullify(globenv%para_env, globenv%blacs_env)
  call cp_para_env_create (globenv%para_env, group=globenv%group,&
      owns_group=.false.,source=globenv%source)
  globenv%print%ref_count=1 ! ugly
  call cp_init_default_logger(para_env=globenv%para_env,&
       print_keys=globenv%print)

END SUBROUTINE start_parallel

!******************************************************************************

SUBROUTINE end_parallel ( globenv )

  TYPE ( global_environment_type ), TARGET :: globenv

!------------------------------------------------------------------------------

  CALL mp_end ( globenv%scr )

END SUBROUTINE end_parallel

!******************************************************************************

END MODULE parallel

!******************************************************************************
!!*****
