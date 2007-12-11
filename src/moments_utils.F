!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/moments_utils [1.0] *
!!
!!   NAME
!!     moments_utils
!!
!!   FUNCTION
!!     Calculates the moment integrals <a|r^m|b>
!!
!!   AUTHOR
!!     JGH (20.07.2006)
!!
!!   MODIFICATION HISTORY
!!     12.2007 [tlaino] - Splitting common routines to QS and FIST
!!
!!   SOURCE
!******************************************************************************

MODULE moments_utils

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE fist_environment_types,          ONLY: fist_environment_type,&
                                             get_fist_env
  USE input_constants,                 ONLY: use_mom_ref_coac,&
                                             use_mom_ref_com,&
                                             use_mom_ref_nondef,&
                                             use_mom_ref_user,&
                                             use_mom_ref_zero
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'moments_utils'

! *** Public subroutines ***

  PUBLIC :: get_reference_point

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE get_reference_point(rpoint,drpoint,qs_env,fist_env,reference,ref_point,error)
    REAL(dp), DIMENSION(3), INTENT(OUT)      :: rpoint
    REAL(dp), DIMENSION(3), INTENT(OUT), &
      OPTIONAL                               :: drpoint
    TYPE(qs_environment_type), OPTIONAL, &
      POINTER                                :: qs_env
    TYPE(fist_environment_type), OPTIONAL, &
      POINTER                                :: fist_env
    INTEGER, INTENT(IN)                      :: reference
    REAL(KIND=dp), DIMENSION(:), POINTER     :: ref_point
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_reference_point', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ia, iatom, ikind
    REAL(dp)                                 :: charge, mass, mtot, ztot
    REAL(dp), DIMENSION(3)                   :: ria
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    NULLIFY ( cell, particle_set, local_particles, para_env )
    IF (PRESENT(qs_env)) THEN
       CALL get_qs_env ( qs_env, cell=cell, particle_set=particle_set, &
            local_particles=local_particles, para_env=para_env, error=error)
    END IF
    IF (PRESENT(fist_env)) THEN
       CALL get_fist_env( fist_env, cell=cell, particle_set=particle_set, &
            local_particles=local_particles, para_env=para_env, error=error)
    END IF
    IF (PRESENT(drpoint)) drpoint = 0.0_dp
    SELECT CASE (reference)
    CASE DEFAULT
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="Type of reference point not implemented", &
            error=error, error_level=cp_failure_level)
    CASE (use_mom_ref_com)
       rpoint = 0._dp
       mtot = 0._dp
       DO ikind = 1,SIZE(local_particles%n_el)
          DO ia = 1,local_particles%n_el(ikind)
             iatom = local_particles%list(ikind)%array(ia)
             ria = particle_set(iatom)%r
             ria = pbc(ria,cell)
             atomic_kind => particle_set(iatom)%atomic_kind
             CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
             rpoint(:) = rpoint(:) + mass*ria(:)
             IF (PRESENT(drpoint)) drpoint = drpoint + mass*particle_set(iatom)%v
             mtot = mtot + mass
          END DO
       END DO
       CALL mp_sum(rpoint,para_env%group)
       CALL mp_sum(mtot,para_env%group)
       IF (ABS(mtot) > 0._dp) THEN 
          rpoint(:) = rpoint(:)/mtot
          IF (PRESENT(drpoint)) drpoint = drpoint / mtot
       END IF
    CASE (use_mom_ref_coac)
       rpoint = 0._dp
       ztot = 0._dp
       DO ikind = 1,SIZE(local_particles%n_el)
          DO ia = 1,local_particles%n_el(ikind)
             iatom = local_particles%list(ikind)%array(ia)
             ria = particle_set(iatom)%r
             ria = pbc(ria,cell)
             atomic_kind => particle_set(iatom)%atomic_kind
             CALL get_atomic_kind(atomic_kind=atomic_kind,core_charge=charge)
             rpoint(:) = rpoint(:) + charge*ria(:)
             IF (PRESENT(drpoint)) drpoint = drpoint + charge*particle_set(iatom)%v
             ztot = ztot + charge
          END DO
       END DO
       CALL mp_sum(rpoint,para_env%group)
       CALL mp_sum(ztot,para_env%group)
       IF (ABS(ztot) > 0._dp) THEN 
          rpoint(:) = rpoint(:)/ztot
          IF (PRESENT(drpoint)) drpoint = drpoint / ztot
       END IF
    CASE (use_mom_ref_user)
       rpoint = ref_point
    CASE (use_mom_ref_zero)
       rpoint = 0._dp
    CASE (use_mom_ref_nondef)
       rpoint = 0._dp
    END SELECT

  END SUBROUTINE get_reference_point

END MODULE moments_utils

! *****************************************************************************
