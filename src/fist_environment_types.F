!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/fist_environment_types [1.0] *
!!
!!   NAME
!!     fist_environment_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM SEPT-01-02 
!!
!!   MODIFICATION HISTORY
!!   gt sept-23-02 added atomic_kind_set to replica_environment_type
!!   to allow use of kind_based neighbor list
!!
!!   SOURCE
!******************************************************************************

MODULE fist_environment_types
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cp_local_lists_types,            ONLY: cp_local_lists_type
  USE ewald_parameters_types,          ONLY: ewald_parameters_type
  USE kinds,                           ONLY: dbl
  USE linklist_types,                  ONLY: linklist_internal_data_type
  USE md,                              ONLY: virial_type
  USE molecule_types,                  ONLY: molecule_structure_type,&
                                             particle_node_type
  USE pair_potential,                  ONLY: potentialparm_type
  USE particle_types,                  ONLY: particle_type
  USE simulation_cell,                 ONLY: cell_type
  PRIVATE

  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='fist_environment_types'

  TYPE pimd_environment_type
     REAL ( dbl ), DIMENSION ( : , : ), POINTER :: u
     REAL ( dbl ) :: lbd
     REAL ( dbl ) :: slbd
     REAL ( dbl ) :: meff
     REAL ( dbl ) :: ismeff
     REAL ( dbl ), DIMENSION ( 2 , 2 ) :: harmonic_prop
  END TYPE pimd_environment_type

  TYPE fist_thermodynamic_type
     TYPE ( virial_type ), POINTER ::  virial
     REAL ( dbl ) :: kin, pot, e_gspace, e_self, e_neut, e_bonded
  END TYPE fist_thermodynamic_type

  TYPE replica_environment_type
    TYPE ( atomic_kind_type ), POINTER :: atomic_kind_set ( : )
    TYPE ( cell_type ), POINTER :: box
    TYPE ( cell_type ), POINTER :: box_ref
    TYPE ( linklist_internal_data_type), POINTER :: ll_data ( : )
    TYPE ( molecule_structure_type ), POINTER :: molecule ( : )
    TYPE ( particle_type ), POINTER :: part ( : )
    TYPE ( particle_node_type ), POINTER :: pnode ( : )
    TYPE ( pimd_environment_type ), POINTER :: pimd_env
    TYPE ( cp_local_lists_type ), POINTER :: local_particles
  END TYPE replica_environment_type

  TYPE fist_environment_type
    LOGICAL, POINTER :: pimd
    TYPE ( ewald_parameters_type ), POINTER :: ewald_param
    TYPE ( potentialparm_type ), POINTER :: potparm ( :, : )
    TYPE ( replica_environment_type ), POINTER :: rep_env ( : )
    TYPE ( fist_thermodynamic_type ), POINTER :: thermo
  END TYPE fist_environment_type

! *** Public data types ***

  PUBLIC :: fist_environment_type,        & 
            replica_environment_type,     & 
            fist_thermodynamic_type,      &
            pimd_environment_type

! *** Public subroutines ***

  PUBLIC :: get_fist_env,&
            init_fist_env,&
            set_fist_env,&
            replica_env_init

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE get_fist_env( fist_env, pimd, rep_env, ewald_param, potparm, thermo )
    IMPLICIT NONE

!   Purpose: Get the FIST environment.

!   ***************************************************************************

    TYPE ( fist_environment_type ), INTENT(IN)              :: fist_env
    LOGICAL, OPTIONAL, POINTER                              :: pimd
    TYPE ( ewald_parameters_type ), OPTIONAL, POINTER       :: ewald_param
    TYPE ( potentialparm_type ), OPTIONAL, POINTER          :: potparm ( :, : )
    TYPE ( fist_thermodynamic_type ), OPTIONAL, POINTER     :: thermo
    TYPE ( replica_environment_type ), OPTIONAL, POINTER    :: rep_env ( : )

!   ---------------------------------------------------------------------------

    IF ( PRESENT ( pimd ) ) pimd => fist_env % pimd
    IF ( PRESENT ( ewald_param ) ) ewald_param => fist_env % ewald_param
    IF ( PRESENT ( potparm ) ) potparm => fist_env % potparm
    IF ( PRESENT ( thermo ) ) thermo => fist_env % thermo
    IF ( PRESENT ( rep_env ) ) rep_env => fist_env % rep_env


  END SUBROUTINE get_fist_env

! *****************************************************************************

  SUBROUTINE init_fist_env ( fist_env )

!   Purpose: Initialise the FIST environment.

!   ***************************************************************************

    TYPE(fist_environment_type), INTENT(OUT) :: fist_env

!   ---------------------------------------------------------------------------
    NULLIFY (fist_env%pimd)
    NULLIFY (fist_env%ewald_param)
    NULLIFY (fist_env%potparm)
    NULLIFY (fist_env%thermo)
    NULLIFY (fist_env%rep_env)

  END SUBROUTINE init_fist_env

! *****************************************************************************

  SUBROUTINE set_fist_env( fist_env, pimd, rep_env, ewald_param, potparm, thermo )

    IMPLICIT NONE
!   Purpose: Set the FIST environment.

    TYPE ( fist_environment_type ), INTENT(out)             :: fist_env
    LOGICAL, OPTIONAL, POINTER                              :: pimd
    TYPE ( ewald_parameters_type ), OPTIONAL, POINTER       :: ewald_param
    TYPE ( potentialparm_type ), OPTIONAL, POINTER          :: potparm ( :, : )
    TYPE ( fist_thermodynamic_type ), OPTIONAL, POINTER     :: thermo
    TYPE ( replica_environment_type ), OPTIONAL, POINTER    :: rep_env ( : )
  
    IF ( PRESENT ( pimd ) ) fist_env % pimd => pimd
    IF ( PRESENT ( ewald_param ) ) fist_env % ewald_param => ewald_param
    IF ( PRESENT ( potparm ) ) fist_env % potparm => potparm
    IF ( PRESENT ( thermo ) ) fist_env % thermo => thermo
    IF ( PRESENT ( rep_env ) ) fist_env % rep_env => rep_env


!***************************************************************************
  END SUBROUTINE set_fist_env
!***************************************************************************

!!****f* fist_environment_type/replica_env_init [1.0] *
!!
!!   NAME
!!     replica_env_init
!!
!!   SYNOPSIS
!!     Subroutine replica_env_init(replica_env)
!!       Type(replica_environment_type), Intent (INOUT):: replica_env
!!     End Subroutine replica_env_init
!!
!!   FUNCTION
!!     initializes a replica environment
!!
!!   MODIFICATION HISTORY
!!     06.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE replica_env_init(replica_env)
  TYPE(replica_environment_type), INTENT(inout) :: replica_env
  
  NULLIFY(replica_env%atomic_kind_set, replica_env%box_ref, &
       replica_env%ll_data, replica_env%molecule, replica_env%part, &
       replica_env%pnode, replica_env%pimd_env, &
       replica_env%local_particles)
END SUBROUTINE replica_env_init
!***************************************************************************

END MODULE fist_environment_types
