!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief initialize scp environment
!> \author CJM
! *****************************************************************************
MODULE scp_environment
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
! USE fist_environment_types,          ONLY: fist_environment_type
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp,&
                                             int_size
  USE machine,                         ONLY: m_flush
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE scp_coeff_types,                 ONLY: aux_coeff_set_create,&
                                             aux_coeff_set_initialize,&
                                             aux_coeff_set_release,&
                                             aux_coeff_set_type
  USE scp_energy_types,                ONLY: allocate_scp_energy,&
                                             scp_energy_type
  USE scp_environment_methods,         ONLY: scp_qs_env_setup
  USE scp_environment_types,           ONLY: scp_environment_type,&
                                             set_scp_env
  USE scp_force_types,                 ONLY: allocate_scp_force,&
                                             scp_force_type
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: scp_init

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'scp_environment'

!-----------------------------------------------------------------------------!

CONTAINS

!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Initializes the SCP environment using the QS environment
! *****************************************************************************
  SUBROUTINE scp_init ( scp_env, qs_env, error )
!  SUBROUTINE scp_init ( scp_env, qs_env, fist_env, error )

    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(qs_environment_type), OPTIONAL, &
      POINTER                                :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scp_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, istat, nkind, &
                                                output_unit
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: natom_of_kind
    LOGICAL                                  :: ionode
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(aux_coeff_set_type), POINTER        :: aux_coeff_set
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(scp_energy_type), POINTER           :: energy
    TYPE(scp_force_type), POINTER            :: force( : )
    TYPE(section_vals_type), POINTER         :: coeff_section, input

!    TYPE(fist_environment_type), OPTIONAL, &
!                                 POINTER     :: fist_env

    CALL timeset(routineN,"I",'',handle)

    ! nullifying  pointers
!    NULLIFY ( scp_control )
    NULLIFY ( aux_coeff_set )
    NULLIFY ( energy, force )
    NULLIFY ( para_env )
    NULLIFY ( input )
    NULLIFY ( atomic_kind_set )
    NULLIFY ( local_particles )

!    IF ( PRESENT ( qs_env ) ) THEN
      CALL get_qs_env (qs_env, para_env = para_env, input = input, &
                       atomic_kind_set = atomic_kind_set, &
                       local_particles = local_particles, error = error )
!    ELSEIF ( PRESENT ( fist_env ) ) THEN
!      CALL get_fist_env (fist_env, para_env = para_env, input = input, &
!                       atomic_kind_set = atomic_kind_set, &
!                       local_particles = local_particles, error = error )
!    ENDIF

    ionode = para_env%ionode
    logger => cp_error_get_logger(error)
    IF(ionode) output_unit = cp_logger_get_default_unit_nr(logger)

!    scp_section =>  section_vals_get_subs_vals( input, "SCP", error = error )
! Set the input section of the SCP
    CALL set_scp_env(scp_env,input=input,&
         error=error)

!    ! Read/write the input section with the dft and qs control parameters
!    CALL read_scp_section ( scp_control, error )

    ! Allocate and initialize coefs pointer
    CALL aux_coeff_set_create(aux_coeff_set=aux_coeff_set,&
                               atomic_kind_set=atomic_kind_set,&
                               distribution=local_particles,&
                               error=error)
!    coeff_section => section_vals_get_subs_vals(input,"SCP%AUX_COEFF",error=error)
    CALL aux_coeff_set_initialize (aux_coeff_set,coeff_section,para_env, &
                                   error=error)
    CALL set_scp_env ( scp_env=scp_env,&
                       aux_coeff_set = aux_coeff_set, &
                       error=error)
    CALL aux_coeff_set_release(aux_coeff_set, error=error)

   ! Allocate the data structure for energies ***
    CALL allocate_scp_energy ( energy )

   ! Allocate the data structure for energies ***
    nkind = SIZE(atomic_kind_set)
    ALLOCATE (natom_of_kind(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "natom_of_kind",nkind*int_size)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                           natom_of_kind=natom_of_kind)
    CALL allocate_scp_force(force,natom_of_kind)
    DEALLOCATE (natom_of_kind,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "natom_of_kind")

    ! set the scp_env
    CALL set_scp_env ( scp_env=scp_env,&
                       energy=energy,&
                       force=force,&
                       error=error)

   !Sets up pw, grids, densities
    IF ( PRESENT ( qs_env ) )  CALL scp_qs_env_setup(scp_env,qs_env,error)
!    IF ( PRESENT ( fist_env ) ) CALL  scp_fist_env_setup(scp_env,fist_env,error)

    IF ( ionode ) CALL m_flush(output_unit)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE scp_init

END MODULE scp_environment

