!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Utilities for replica in BAND method
!>        Module split from other band modules.
!> \author Teodoro Laino 06.2009
! *****************************************************************************
MODULE band_replica_utils
  USE f77_blas
  USE input_constants,                 ONLY: do_band_cartesian
  USE kinds,                           ONLY: dp
  USE neb_types,                       ONLY: neb_type,&
                                             neb_var_type
  USE particle_types,                  ONLY: particle_type
  USE rmsd,                            ONLY: rmsd3
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'band_replica_utils'
  LOGICAL, PARAMETER, PRIVATE :: debug_this_module = .FALSE.

  PUBLIC :: neb_replica_distance,&
            dot_product_band

CONTAINS

! *****************************************************************************
!> \brief Computes the distance between two replica (chord length)
!> \author Teodoro Laino 09.2006
! *****************************************************************************
  SUBROUTINE neb_replica_distance(particle_set, coords, i0, i, distance, iw, rotate, error)
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particle_set
    TYPE(neb_var_type), POINTER              :: coords
    INTEGER, INTENT(IN)                      :: i0, i
    REAL(KIND=dp), INTENT(OUT)               :: distance
    INTEGER, INTENT(IN)                      :: iw
    LOGICAL, INTENT(IN), OPTIONAL            :: rotate
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'neb_replica_distance', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure, my_rotate

    failure = .FALSE.
    IF (.NOT.failure) THEN
       my_rotate = .FALSE.
       IF (PRESENT(rotate)) my_rotate=rotate
       ! The rotation of the replica is enabled exclusively when working in
       ! cartesian coordinates
       IF (my_rotate.AND.(coords%in_use==do_band_cartesian)) THEN
          CPPostcondition(PRESENT(particle_set),cp_failure_level,routineP,error,failure)
          CALL rmsd3(particle_set, coords%xyz(:,i), coords%xyz(:,i0), iw,&
                     rotate=my_rotate, error=error)
       END IF
       distance = SQRT(DOT_PRODUCT(coords%wrk(:,i)-coords%wrk(:,i0),&
                                   coords%wrk(:,i)-coords%wrk(:,i0)))
    END IF
  END SUBROUTINE neb_replica_distance

! *****************************************************************************
!> \brief  Handles the dot_product when using colvar.. in this case
!>         the scalar product needs to take into account the metric
!>         tensor
!> \author Teodoro Laino 09.2006
! *****************************************************************************
  FUNCTION dot_product_band(neb_env, ind, array1, array2, array3, error) RESULT(value)
    TYPE(neb_type), POINTER                  :: neb_env
    INTEGER, INTENT(IN)                      :: ind
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: array1, array2
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: array3
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: value

    CHARACTER(len=*), PARAMETER :: routineN = 'dot_product_band', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nsize_int
    LOGICAL                                  :: check, failure

    failure = .FALSE.
    IF (neb_env%use_colvar) THEN
       nsize_int = neb_env%nsize_int
       check =((SIZE(array1)/=SIZE(array2)).OR.&
               (SIZE(array1)/=nsize_int).OR.&
               (SIZE(array3)/=nsize_int*nsize_int))
       ! This condition should always be satisfied..
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
       value = DOT_PRODUCT(MATMUL(RESHAPE(array3,(/nsize_int,nsize_int/)),array1),array2)
    ELSE
       value = DOT_PRODUCT(array1,array2)
    END IF
  END FUNCTION dot_product_band

END MODULE band_replica_utils
