!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2003  CP2K developers group                            !
!-----------------------------------------------------------------------------!
!!****s* cp2k/md_environment_types [1.0] *
!!
!!   NAME
!!     md_environment_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM SEPT-12-02
!!
!!   MODIFICATION HISTORY
!!     give the md env its own para_env Joost VandeVondele 07.2003
!!
!!   SOURCE
!******************************************************************************

MODULE md_environment_types
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE extended_system_dynamics,        ONLY: lnhc_parameters_type,&
                                             npt_info_type
  USE force_env_types,                   ONLY: force_env_release,&
                                             force_env_retain,&
                                             force_env_type
  USE kinds,                           ONLY: dbl
  USE md,                              ONLY: simulation_parameters_type,&
                                             virial_type
  USE simulation_cell,                 ONLY: cell_type
  USE structure_types,                 ONLY: structure_type
  IMPLICIT NONE

  PRIVATE

  ! para_env is the parallel enviroment of the MD, i.e. the systems that are dealt with by the integrator
  ! e.g in the PIMD this could be parent of every bead

  TYPE md_environment_type
     TYPE(cp_para_env_type), POINTER :: para_env 
     INTEGER, POINTER :: itimes
     REAL ( dbl ), POINTER :: constant
     TYPE ( cell_type ), POINTER :: cell
     TYPE ( lnhc_parameters_type ), POINTER :: nhc_part ( : )
     TYPE ( lnhc_parameters_type ), POINTER :: nhc_coef 
     TYPE ( lnhc_parameters_type ), POINTER :: nhc_baro ( : )
     TYPE ( npt_info_type ), POINTER :: npt ( :, : )
     TYPE ( simulation_parameters_type ), POINTER :: simpar
     TYPE ( structure_type ), POINTER :: struc ( : )
     TYPE ( virial_type ), POINTER :: virial
     TYPE ( force_env_type ), POINTER :: force_env  !? should this also be an array for PIMD
  END TYPE md_environment_type

! *** Public subroutines and data types ***
  PUBLIC :: md_environment_type, set_md_env,   &
            virial_type, get_md_env, init_md_env, destroy_md_env, &
            zero_virial

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE get_md_env ( md_env, itimes, constant, cell, simpar, struc,  &
                          virial, nhc_part, nhc_coef, nhc_baro, npt, &
                          force_env, para_env )

    IMPLICIT NONE
!   Purpose: Get the integrator environment

!   ***************************************************************************

    TYPE ( md_environment_type ), INTENT ( IN ) :: md_env
    TYPE ( virial_type ), POINTER, OPTIONAL :: virial
    TYPE ( cell_type ), POINTER, OPTIONAL :: cell
    TYPE ( simulation_parameters_type ), POINTER, OPTIONAL :: simpar
    TYPE ( structure_type ), POINTER, OPTIONAL :: struc ( : )
    TYPE ( lnhc_parameters_type ), POINTER, OPTIONAL :: nhc_part ( : )
    TYPE ( lnhc_parameters_type ), POINTER, OPTIONAL :: nhc_coef 
    TYPE ( lnhc_parameters_type ), POINTER, OPTIONAL :: nhc_baro ( : )
    TYPE ( npt_info_type ), POINTER, OPTIONAL :: npt ( :, : )
    TYPE ( force_env_type ), POINTER, OPTIONAL :: force_env
    INTEGER, POINTER, OPTIONAL :: itimes
    REAL ( dbl ), POINTER, OPTIONAL :: constant
    TYPE(cp_para_env_type), POINTER, OPTIONAL :: para_env 
!   ---------------------------------------------------------------------------

    IF ( PRESENT ( itimes ) ) itimes => md_env % itimes
    IF ( PRESENT ( constant ) ) constant => md_env % constant
    IF ( PRESENT ( cell ) ) cell => md_env % cell 
    IF ( PRESENT ( simpar ) ) simpar => md_env % simpar
    IF ( PRESENT ( struc ) ) struc => md_env % struc
    IF ( PRESENT ( nhc_part ) ) nhc_part => md_env % nhc_part 
    IF ( PRESENT ( nhc_coef ) ) nhc_coef => md_env % nhc_coef
    IF ( PRESENT ( nhc_baro ) ) nhc_baro => md_env % nhc_baro
    IF ( PRESENT ( npt ) ) npt => md_env % npt
    IF ( PRESENT ( virial ) ) virial => md_env % virial
    IF ( PRESENT ( force_env ) ) force_env => md_env % force_env
    IF ( PRESENT ( para_env ) ) para_env => md_env % para_env

  END SUBROUTINE get_md_env

! *****************************************************************************

  SUBROUTINE init_md_env ( md_env, para_env )

!   Purpose: Initialise the integrator environment.
!   retain the para_env for this environment (should be used for parallel
!   communications)

!   ***************************************************************************

    TYPE(md_environment_type), INTENT(OUT) :: md_env
    TYPE(cp_para_env_type), POINTER :: para_env 

!   ---------------------------------------------------------------------------
    NULLIFY ( md_env % itimes )
    NULLIFY ( md_env % constant )
    NULLIFY ( md_env % cell )
    NULLIFY ( md_env % simpar )
    NULLIFY ( md_env % struc )
    NULLIFY ( md_env % nhc_part )
    NULLIFY ( md_env % nhc_coef )
    NULLIFY ( md_env % nhc_baro )
    NULLIFY ( md_env % npt )
    NULLIFY ( md_env % virial )
    NULLIFY ( md_env % force_env )
    md_env % para_env => para_env
    CALL cp_para_env_retain(md_env % para_env)

  END SUBROUTINE init_md_env

  SUBROUTINE destroy_md_env ( md_env )

!   Purpose: destroy the integrator environment.
!   release the para_env for this environment 
!   ***************************************************************************

    TYPE(md_environment_type), INTENT(OUT) :: md_env

!   ---------------------------------------------------------------------------
    CALL cp_para_env_release(md_env % para_env)
    NULLIFY ( md_env % itimes )
    NULLIFY ( md_env % constant )
    NULLIFY ( md_env % cell )
    NULLIFY ( md_env % simpar )
    NULLIFY ( md_env % struc )
    NULLIFY ( md_env % nhc_part )
    NULLIFY ( md_env % nhc_coef )
    NULLIFY ( md_env % nhc_baro )
    NULLIFY ( md_env % npt )
    NULLIFY ( md_env % virial )
    NULLIFY ( md_env % force_env )

  END SUBROUTINE destroy_md_env

! *****************************************************************************

  SUBROUTINE set_md_env ( md_env, itimes, constant, cell, simpar, struc,  &
                          virial, nhc_part, nhc_coef, nhc_baro, npt,&
                          force_env)

!   Purpose: Set the integrator environment to the correct program.


!   ***************************************************************************
    IMPLICIT NONE
    TYPE ( md_environment_type ), INTENT ( OUT ) :: md_env
    TYPE ( cell_type ), POINTER , OPTIONAL :: cell
    TYPE ( virial_type ), POINTER, OPTIONAL :: virial
    TYPE ( simulation_parameters_type ), POINTER, OPTIONAL :: simpar
    TYPE ( structure_type ), POINTER, OPTIONAL  :: struc ( : )
    TYPE ( lnhc_parameters_type ), POINTER, OPTIONAL :: nhc_part ( : )
    TYPE ( lnhc_parameters_type ), POINTER, OPTIONAL :: nhc_coef 
    TYPE ( lnhc_parameters_type ), POINTER, OPTIONAL :: nhc_baro ( : )
    TYPE ( npt_info_type ), POINTER, OPTIONAL :: npt ( :, : )
    TYPE ( force_env_type ), POINTER, OPTIONAL :: force_env
    INTEGER, POINTER, OPTIONAL :: itimes
    REAL ( dbl ), POINTER, OPTIONAL :: constant

    IF ( PRESENT ( cell ) ) md_env % cell => cell
    IF ( PRESENT ( struc ) ) md_env % struc => struc
    IF ( PRESENT ( simpar ) ) md_env % simpar => simpar
    IF ( PRESENT ( nhc_part ) ) md_env % nhc_part => nhc_part
    IF ( PRESENT ( nhc_coef ) ) md_env % nhc_coef => nhc_coef
    IF ( PRESENT ( nhc_baro ) ) md_env % nhc_baro => nhc_baro
    IF ( PRESENT ( npt ) ) md_env % npt => npt
    IF ( PRESENT ( virial ) ) md_env % virial => virial
    IF ( PRESENT ( itimes ) ) md_env % itimes => itimes
    IF ( PRESENT ( constant ) ) md_env % constant => constant
    IF ( PRESENT ( force_env ) ) THEN
       CALL force_env_retain( force_env ) ! accept null pointers?
       CALL force_env_release( md_env % force_env )
       md_env % force_env => force_env
    END IF
    
!   ---------------------------------------------------------------------------

  END SUBROUTINE set_md_env

! *****************************************************************************

  SUBROUTINE zero_virial ( virial )

!   Purpose:  zero the virials

!   ***************************************************************************

    TYPE (virial_type), INTENT(OUT) :: virial

    virial % pv_total = 0._dbl 
    virial % pv_kinetic = 0._dbl 
    virial % pv_virial = 0._dbl 
    virial % pv_constraint = 0._dbl 
    virial % pv_availability = .TRUE.

!   ---------------------------------------------------------------------------

  END SUBROUTINE  zero_virial

! *****************************************************************************

!***************************************************************************

END MODULE md_environment_types

!******************************************************************************
