!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2006  CP2K developers group                            !
!-----------------------------------------------------------------------------!

!!****h* cp2k/tbmd_environment_methods [1.0] *
!!
!!   NAME
!!     tbmd_environment
!!
!!   FUNCTION
!!     Methods on the tight-binding environment
!!
!!   NOTES
!!     -
!!     
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!   SOURCE
!****************************************************************************
MODULE tbmd_environment_methods
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_2d_types,           ONLY: distribution_2d_type
  USE distribution_methods,            ONLY: distribute_molecules_2d
  USE kinds,                           ONLY: dp
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE particle_types,                  ONLY: particle_type
  USE qs_matrix_pools,                 ONLY: mpools_rebuild_fm_pools
  USE qs_mo_types,                     ONLY: mo_set_p_type,&
                                             mo_set_type
  USE tbmd_environment_types,          ONLY: tbmd_env_get,&
                                             tbmd_env_set,&
                                             tbmd_environment_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tbmd_environment_methods'

! *** Public subroutines ***

  PUBLIC :: tbmd_env_setup

!!***
! *****************************************************************************

CONTAINS

!!****f* tbmd_environment_methods/tbmd_env_setup [1.0] *
!!
!!   NAME
!!     tbmd_env_setup
!!
!!   FUNCTION
!!     Initialize the tbmd environment
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - tbmd_env: The tbmd environment to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE tbmd_env_setup(tbmd_env, error)
    TYPE(tbmd_environment_type), POINTER     :: tbmd_env
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_env_setup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos_ptr
    TYPE(mo_set_type), POINTER               :: mos
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

!   ------------------------------------------------------------------------

    CALL timeset(routineN, "I", "", handle)

    failure = .FALSE.

    NULLIFY(distribution_2d, atomic_kind_set, molecule_kind_set, molecule_set, &
            particle_set, para_env, mos, mos_ptr)
  
    CPPrecondition(ASSOCIATED(tbmd_env), cp_failure_level, routineP, error, failure)
 
    IF (.NOT. failure) THEN

      CALL tbmd_env_get(tbmd_env=tbmd_env, atomic_kind_set=atomic_kind_set, &
                        particle_set=particle_set, para_env=para_env, mos=mos, &
                        molecule_kind_set=molecule_kind_set, molecule_set=molecule_set, &
                        error=error)
      CPPrecondition(ASSOCIATED(atomic_kind_set), cp_failure_level, routineP, error, failure)

      ALLOCATE(mos_ptr(1), STAT=stat)
      CPPostcondition(stat==0, cp_failure_level, routineN, error, failure)
      mos_ptr(1)%mo_set => mos

      CALL mpools_rebuild_fm_pools(mpools=tbmd_env%matrix_pools, mos=mos_ptr, &
                                   blacs_env=tbmd_env%blacs_env, para_env=para_env, &
                                   nrow_block=32, ncol_block=32, error=error)

      DEALLOCATE(mos_ptr, STAT=stat)
      CPPostcondition(stat==0, cp_failure_level, routineN, error, failure)

      CALL distribute_molecules_2d(particle_kind_set=atomic_kind_set, &
                                   particle_set=particle_set, &
                                   molecule_kind_set=molecule_kind_set, &
                                   molecule_set=molecule_set, &
                                   distribution_2d=distribution_2d, &
                                   blacs_env=tbmd_env%blacs_env, &
                                   force_env_section=tbmd_env%force_env_input, &
                                   error=error)

      CALL tbmd_env_set(tbmd_env=tbmd_env, distribution_2d=distribution_2d, &
                        error=error)

    END IF

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE tbmd_env_setup
!***************************************************************************

END MODULE tbmd_environment_methods
