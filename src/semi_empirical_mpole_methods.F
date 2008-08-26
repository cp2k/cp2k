!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Setup and Methods for semi-empirical multipole types 
!> \author Teodoro Laino [tlaino] - 08.2008 Zurich University
! *****************************************************************************
MODULE semi_empirical_mpole_methods
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE semi_empirical_mpole_types,      ONLY: nddo_mpole_create,&
                                             nddo_mpole_release,&
                                             nddo_mpole_type,&
                                             semi_empirical_mpole_p_create,&
                                             semi_empirical_mpole_p_release,&
                                             semi_empirical_mpole_p_type
  USE semi_empirical_types,            ONLY: semi_empirical_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'semi_empirical_mpole_methods'

  PUBLIC :: semi_empirical_mpole_p_setup,&
            nddo_mpole_setup

CONTAINS

! *****************************************************************************
!> \brief Setup semi-empirical mpole type
!>        This function setup for each semi-empirical type a structure containing
!>        the multipolar expansion for all possible combination on-site of atomic
!>        orbitals ( \mu \nu |
! *****************************************************************************
  SUBROUTINE semi_empirical_mpole_p_setup(mpole, se_parameter, error)
    TYPE(semi_empirical_mpole_p_type), &
      DIMENSION(:), POINTER                  :: mpole
    TYPE(semi_empirical_type), POINTER       :: se_parameter
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'semi_empirical_mpole_p_setup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: natorb, ndim
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(mpole)) THEN
       CALL semi_empirical_mpole_p_release(mpole, error)
    END IF
    natorb = se_parameter%natorb
    ndim   = natorb*(natorb+1)/2
    CALL semi_empirical_mpole_p_create(mpole, ndim, error)

    ! Fill in information on multipole expansion due to atomic orbitals charge 
    ! distribution
    

  END SUBROUTINE semi_empirical_mpole_p_setup


! *****************************************************************************
!> \brief Setup NDDO multipole type
! *****************************************************************************
  SUBROUTINE nddo_mpole_setup(nddo_mpole, natom, error)
    TYPE(nddo_mpole_type), POINTER           :: nddo_mpole
    INTEGER, INTENT(IN)                      :: natom
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'nddo_mpole_setup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(nddo_mpole)) THEN
       CALL nddo_mpole_release(nddo_mpole, error=error)
    END IF
    CALL nddo_mpole_create(nddo_mpole, error=error)
    ! Allocate Global Arrays
    ALLOCATE(nddo_mpole%charge(natom),stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(nddo_mpole%dipole(3,natom),stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(nddo_mpole%quadrupole(3,3,natom),stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(nddo_mpole%efield0(  natom),stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(nddo_mpole%efield1(3,natom),stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(nddo_mpole%efield2(9,natom),stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE nddo_mpole_setup
  
END MODULE semi_empirical_mpole_methods
