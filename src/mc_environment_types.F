!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2003  CP2K developers group                            !
!-----------------------------------------------------------------------------!
!!****s* cp2k/mc_environment_types [1.0] *
!!
!!   NAME
!!     mc_environment_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     MJM Oct. 15-2003
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE mc_environment_types
  USE cp_para_types,                   ONLY: cp_blacs_env_type
  USE force_env_types,                 ONLY: force_env_release,&
                                             force_env_retain,&
                                             force_env_type
  USE mc_control,                      ONLY: mc_simulation_parameters_type

  IMPLICIT NONE

  PRIVATE

 
  TYPE mc_environment_type
     TYPE ( mc_simulation_parameters_type ), POINTER :: mc_par
     TYPE ( force_env_type ), POINTER :: force_env  
     TYPE ( cp_blacs_env_type ), POINTER :: blacs_env  
  END TYPE mc_environment_type

! *** Public subroutines and data types ***
  PUBLIC :: mc_environment_type, set_mc_env, &
            get_mc_env, init_mc_env, destroy_mc_env

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE get_mc_env ( mc_env, mc_par, force_env, blacs_env)

    IMPLICIT NONE
!   Purpose: Get the integrator environment...works differently than
!            set_mc_env

!   ***************************************************************************

    TYPE ( mc_environment_type ), INTENT ( IN ) :: mc_env
    TYPE ( mc_simulation_parameters_type ), POINTER, OPTIONAL :: mc_par
    TYPE ( force_env_type ), POINTER, OPTIONAL :: force_env
    TYPE ( cp_blacs_env_type ), POINTER, OPTIONAL :: blacs_env
!   ---------------------------------------------------------------------------

    IF ( PRESENT ( mc_par ) ) mc_par => mc_env % mc_par
    IF ( PRESENT ( force_env ) ) force_env => mc_env % force_env
    IF ( PRESENT ( blacs_env ) ) blacs_env => mc_env % blacs_env

  END SUBROUTINE get_mc_env

! *****************************************************************************

  SUBROUTINE init_mc_env ( mc_env )

!   Purpose: Initialise the integrator environment.
!   retain the para_env for this environment (should be used for parallel
!   communications)

!   ***************************************************************************

    TYPE(mc_environment_type), INTENT(OUT) :: mc_env

!   ---------------------------------------------------------------------------
    NULLIFY ( mc_env % mc_par )
    NULLIFY ( mc_env % force_env )
    NULLIFY ( mc_env % blacs_env )

  END SUBROUTINE init_mc_env

  SUBROUTINE destroy_mc_env ( mc_env )

!   Purpose: destroy the integrator environment.
!   release the para_env for this environment 
!   ***************************************************************************

    TYPE(mc_environment_type), INTENT(OUT) :: mc_env

!   ---------------------------------------------------------------------------
    NULLIFY ( mc_env % mc_par )
    NULLIFY ( mc_env % force_env )
    NULLIFY ( mc_env % blacs_env )

  END SUBROUTINE destroy_mc_env

! *****************************************************************************

  SUBROUTINE set_mc_env ( mc_env, mc_par, force_env, blacs_env )

!   Purpose: Set the integrator environment to the correct program.


!   ***************************************************************************
    IMPLICIT NONE
    TYPE ( mc_environment_type ), INTENT ( OUT ) :: mc_env
    TYPE ( mc_simulation_parameters_type ), POINTER, OPTIONAL :: mc_par
    TYPE ( force_env_type ), POINTER, OPTIONAL :: force_env
    TYPE ( cp_blacs_env_type ), POINTER, OPTIONAL :: blacs_env

    IF ( PRESENT ( mc_par ) ) mc_env % mc_par => mc_par
    IF ( PRESENT ( blacs_env ) ) mc_env % blacs_env => blacs_env
    IF ( PRESENT ( force_env ) ) THEN
       CALL force_env_retain( force_env ) ! accept null pointers?
       CALL force_env_release( mc_env % force_env )
       mc_env % force_env => force_env
    END IF
    
!   ---------------------------------------------------------------------------

  END SUBROUTINE set_mc_env

! *****************************************************************************

! *****************************************************************************

!***************************************************************************

END MODULE mc_environment_types

!******************************************************************************
