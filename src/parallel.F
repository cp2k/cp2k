!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2003  CP2K developers group                            !
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

  USE cp_para_env,                     ONLY: cp_para_env_create
  USE global_types,                    ONLY: global_environment_type
  USE message_passing,                 ONLY: mp_end,&
                                             mp_start,&
                                             mp_world_setup

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

    TYPE(global_environment_type), TARGET    :: globenv

!------------------------------------------------------------------------------

  CALL mp_start ( )
  CALL mp_world_setup ( globenv % num_pe, globenv % mepos, globenv % group )
  globenv % source = 0
  globenv % ionode = (globenv % mepos==globenv % source)
  nullify(globenv%para_env)
  call cp_para_env_create (globenv%para_env, group=globenv%group,&
      owns_group=.false.,source=globenv%source)

END SUBROUTINE start_parallel

!******************************************************************************

SUBROUTINE end_parallel ( globenv )

    TYPE(global_environment_type), TARGET    :: globenv

!------------------------------------------------------------------------------

  CALL mp_end ( globenv%scr )

END SUBROUTINE end_parallel

!******************************************************************************

END MODULE parallel

!******************************************************************************
!!*****
