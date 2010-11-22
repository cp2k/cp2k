!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Collections of functions working on replicas
!>        Module split from other band modules.
!> \author Teodoro Laino 06.2009
! *****************************************************************************
MODULE band_replica_methods
  USE band_replica_utils,              ONLY: neb_replica_distance
  USE f77_blas
  USE input_constants,                 ONLY: do_band_cartesian
  USE kinds,                           ONLY: dp
  USE mathlib,                         ONLY: matvec_3x3
  USE neb_types,                       ONLY: neb_var_type
  USE particle_types,                  ONLY: particle_type
  USE rmsd,                            ONLY: rmsd3
  USE spline_env_replica_methods,      ONLY: eval_spline_replica,&
                                             reparametrize_spline_replica
  USE spline_env_replica_types,        ONLY: create_spline_env_replica,&
                                             release_spline_env_replica,&
                                             spline_env_replica_type
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'band_replica_methods'
  LOGICAL, PARAMETER, PRIVATE :: debug_this_module = .FALSE.

  PUBLIC :: reorient_images,&
            reparametrize_images

CONTAINS

! *****************************************************************************
!> \brief Reorient iteratively all images of the NEB chain in order to
!>      have always the smaller RMSD between two following images
!> \author Teodoro Laino 09.2006
! *****************************************************************************
  SUBROUTINE reorient_images(rotate_frames, particle_set, coords, vels, iw,&
       distances, number_of_replica, error)
    LOGICAL, INTENT(IN)                      :: rotate_frames
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particle_set
    TYPE(neb_var_type), POINTER              :: coords, vels
    INTEGER, INTENT(IN)                      :: iw
    REAL(KIND=dp), DIMENSION(:), OPTIONAL    :: distances
    INTEGER, INTENT(IN)                      :: number_of_replica
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'reorient_images', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, k, kind
    LOGICAL                                  :: check, failure
    REAL(KIND=dp), DIMENSION(3)              :: tmp
    REAL(KIND=dp), DIMENSION(3, 3)           :: rot

    failure = .FALSE.
    IF (.NOT.failure) THEN
       rot      = 0.0_dp
       rot(1,1) = 1.0_dp
       rot(2,2) = 1.0_dp
       rot(3,3) = 1.0_dp
       DO i = 2, number_of_replica
          ! The rotation of the replica is enabled exclusively when working in
          ! cartesian coordinates
          IF (rotate_frames.AND.(coords%in_use==do_band_cartesian)) THEN
             CALL rmsd3(particle_set, coords%xyz(:,i), coords%xyz(:,i-1), iw,&
                  rotate=.TRUE., rot=rot, error=error)
             ! Rotate velocities
             DO k = 1, SIZE(vels%xyz,1)/3
                kind = (k-1)*3
                tmp  = vels%xyz(kind+1:kind+3,i)
                CALL matvec_3x3(vels%xyz(kind+1:kind+3,i),TRANSPOSE(rot),tmp)
             END DO
          END IF
          IF (PRESENT(distances)) THEN
             check = SIZE(distances)==(number_of_replica-1)
             CPPostcondition(check,cp_failure_level,routineP,error,failure)
             CALL neb_replica_distance(coords=coords, i0=i-1, i=i, distance=distances(i-1),&
                  iw=iw, rotate=.FALSE., error=error)
          END IF
       END DO
    END IF
  END SUBROUTINE reorient_images

! *****************************************************************************
!> \brief Reparametrization of the replica for String Method with splines
!> \author Teodoro Laino - Rodolphe Vuilleumier 09.2008
!> \par History
!>      06.2009    Teodoro Laino [tlaino] - Generalization to order n B-splines
! *****************************************************************************
  SUBROUTINE reparametrize_images(reparametrize_frames, spline_order, smoothing,&
       coords, sline, distances, error)

    LOGICAL, INTENT(IN)                      :: reparametrize_frames
    INTEGER, INTENT(IN)                      :: spline_order
    REAL(KIND=dp), INTENT(IN)                :: smoothing
    TYPE(neb_var_type), POINTER              :: coords, sline
    REAL(KIND=dp), DIMENSION(:)              :: distances
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'reparametrize_images', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, nreplica, nxyz, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: tmp_coords
    TYPE(spline_env_replica_type), POINTER   :: spline_env_replica

    failure = .FALSE.
    IF (reparametrize_frames) THEN
       NULLIFY(spline_env_replica)
       nxyz     = SIZE(coords%wrk,1)
       nreplica = SIZE(coords%wrk,2)
       ALLOCATE(tmp_coords(nxyz,nreplica),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       tmp_coords = coords%wrk
       ! Perform Smoothing
       DO i = 2, nreplica-1
          coords%wrk(:,i) = tmp_coords(:,i)*(1.0_dp-2.0_dp*smoothing) + &
               tmp_coords(:,i-1)*smoothing + tmp_coords(:,i+1)*smoothing
       END DO
       sline%wrk = coords%wrk - tmp_coords + sline%wrk
       tmp_coords =  coords%wrk

       ! Spline and reparametrize
       CALL create_spline_env_replica(spline_env_replica, spline_order, nreplica,&
                                      nxyz, error)
       CALL eval_spline_replica(spline_env_replica, coords, nreplica, nxyz, error)
       CALL reparametrize_spline_replica(spline_env_replica, coords, distances,&
                                         nreplica, nxyz, error)
       CALL release_spline_env_replica(spline_env_replica, error)
       STOP "spline-neb"
       ! Update sline with the correct vector (forces + reparameterization)
       sline%wrk = coords%wrk - tmp_coords + sline%wrk
       DEALLOCATE(tmp_coords,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
  END SUBROUTINE reparametrize_images

END MODULE band_replica_methods
