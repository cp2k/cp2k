!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      JGH [04042007] refactoring
! *****************************************************************************
MODULE virial_types

  USE f77_blas
  USE kinds,                           ONLY: dp
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'virial_types'
  INTEGER, PRIVATE, SAVE :: last_virial_id_nr=0

  PUBLIC:: virial_type, virial_p_type, virial_set, virial_get, cp_virial,&
       virial_create, virial_retain, virial_release, zero_virial, sym_virial

! *****************************************************************************
  TYPE virial_type
     INTEGER                            :: ref_count, id_nr
     REAL (KIND=dp), DIMENSION ( 3, 3 ) :: pv_total, &
                                           pv_virial, &
                                           pv_kinetic, &
                                           pv_constraint, &
                                           pv_xc, &
                                           pv_fock_4c
     LOGICAL                            :: pv_availability,&
                                           pv_calculate,&
                                           pv_numer,&
                                           pv_diagonal
  END TYPE virial_type

! *****************************************************************************
  TYPE virial_p_type
     TYPE(virial_type),  POINTER        :: virial
  END TYPE virial_p_type

CONTAINS

! *****************************************************************************
!> \brief   copy virial_in into virial_out
!> \version 1.0
! *****************************************************************************
  SUBROUTINE cp_virial ( virial_in, virial_out )
    TYPE(virial_type), POINTER               :: virial_in, virial_out

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_virial', &
      routineP = moduleN//':'//routineN

    virial_out % pv_total = virial_in % pv_total
    virial_out % pv_kinetic = virial_in % pv_kinetic
    virial_out % pv_virial = virial_in % pv_virial
    virial_out % pv_xc = virial_in % pv_xc
    virial_out % pv_fock_4c = virial_in % pv_fock_4c
    virial_out % pv_constraint = virial_in % pv_constraint
    virial_out % pv_availability = virial_in % pv_availability
    virial_out % pv_calculate = virial_in % pv_calculate
    virial_out % pv_numer = virial_in % pv_numer
    virial_out % pv_diagonal = virial_in % pv_diagonal

  END SUBROUTINE  cp_virial

! *****************************************************************************
!> \brief   symmetrize the virial components
!> \version 1.0
! *****************************************************************************
  SUBROUTINE sym_virial (virial, error)
    TYPE(virial_type), POINTER               :: virial
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'sym_virial', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j

    DO i = 1,3
       DO j = 1, i-1
          ! Symmetrize total
          virial%pv_total(j,i) = (virial%pv_total(i,j) + virial%pv_total(j,i))*0.5_dp
          virial%pv_total(i,j) =  virial%pv_total(j,i)
          ! Symmetrize Kinetic
          virial%pv_kinetic(j,i) = (virial%pv_kinetic(i,j) + virial%pv_kinetic(j,i))*0.5_dp
          virial%pv_kinetic(i,j) =  virial%pv_kinetic(j,i)
          ! Symmetrize Virial
          virial%pv_virial(j,i) = (virial%pv_virial(i,j) + virial%pv_virial(j,i))*0.5_dp
          virial%pv_virial(i,j) =  virial%pv_virial(j,i)
          ! Symmetrize XC
          virial%pv_xc(j,i) = (virial%pv_xc(i,j) + virial%pv_xc(j,i))*0.5_dp
          virial%pv_xc(i,j) =  virial%pv_xc(j,i)
          ! Symmetrize HFX
          virial%pv_fock_4c(j,i) = (virial%pv_fock_4c(i,j) + virial%pv_fock_4c(j,i))*0.5_dp
          virial%pv_fock_4c(i,j) =  virial%pv_fock_4c(j,i)

          ! Symmetrize constraints
          virial%pv_constraint(j,i) = (virial%pv_constraint(i,j) + virial%pv_constraint(j,i))*0.5_dp
          virial%pv_constraint(i,j) =  virial%pv_constraint(j,i)
       END DO
    END DO

  END SUBROUTINE sym_virial

! *****************************************************************************
  SUBROUTINE zero_virial ( virial, reset )

    TYPE(virial_type), INTENT(OUT)           :: virial
    LOGICAL, INTENT(IN), OPTIONAL            :: reset

    CHARACTER(LEN=*), PARAMETER :: routineN = 'zero_virial', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: my_reset

    my_reset = .TRUE.
    IF (PRESENT(reset)) my_reset = reset

    virial % pv_total = 0.0_dp
    virial % pv_kinetic = 0.0_dp
    virial % pv_virial = 0.0_dp
    virial % pv_xc = 0.0_dp
    virial % pv_fock_4c = 0.0_dp
    virial % pv_constraint = 0.0_dp
    IF ( my_reset ) THEN
       virial % pv_availability = .FALSE.
       virial % pv_numer = .FALSE.
       virial % pv_calculate = .FALSE.
       virial % pv_diagonal = .FALSE.
    END IF
  END SUBROUTINE  zero_virial

! *****************************************************************************
  SUBROUTINE virial_set ( virial, pv_total, pv_virial, pv_xc, pv_fock_4c, pv_constraint, &
       pv_kinetic, pv_availability, pv_calculate, pv_numer, pv_diagonal )

    TYPE(virial_type), POINTER               :: virial
    REAL(KIND=dp), DIMENSION(3, 3), OPTIONAL :: pv_total, pv_virial, pv_xc, &
                                                pv_fock_4c, pv_constraint, &
                                                pv_kinetic
    LOGICAL, OPTIONAL                        :: pv_availability, &
                                                pv_calculate, pv_numer, &
                                                pv_diagonal

    CHARACTER(LEN=*), PARAMETER :: routineN = 'virial_set', &
      routineP = moduleN//':'//routineN

    IF ( PRESENT ( pv_total ) ) virial % pv_total = pv_total
    IF ( PRESENT ( pv_virial ) ) virial % pv_virial = pv_virial
    IF ( PRESENT ( pv_xc ) ) virial % pv_xc = pv_xc
    IF ( PRESENT ( pv_fock_4c ) ) virial % pv_fock_4c = pv_fock_4c
    IF ( PRESENT ( pv_constraint ) ) virial % pv_constraint = pv_constraint
    IF ( PRESENT ( pv_kinetic ) ) virial % pv_kinetic = pv_kinetic
    IF ( PRESENT ( pv_availability ) ) virial % pv_availability = pv_availability
    IF ( PRESENT ( pv_calculate ) ) virial % pv_calculate = pv_calculate
    IF ( PRESENT ( pv_numer ) ) virial % pv_numer = pv_numer
    IF ( PRESENT ( pv_diagonal ) ) virial % pv_diagonal = pv_diagonal

  END SUBROUTINE virial_set

! *****************************************************************************
  SUBROUTINE virial_get ( virial, pv_total, pv_virial, pv_xc, pv_fock_4c, pv_constraint, &
       pv_kinetic, pv_availability, pv_calculate, pv_numer, pv_diagonal )

    TYPE(virial_type), POINTER               :: virial
    REAL(KIND=dp), DIMENSION(3, 3), OPTIONAL :: pv_total, pv_virial, pv_xc, &
                                                pv_fock_4c, pv_constraint, &
                                                pv_kinetic
    LOGICAL, OPTIONAL                        :: pv_availability, &
                                                pv_calculate, pv_numer, &
                                                pv_diagonal

    CHARACTER(LEN=*), PARAMETER :: routineN = 'virial_get', &
      routineP = moduleN//':'//routineN

    IF ( PRESENT ( pv_total ) ) pv_total = virial % pv_total
    IF ( PRESENT ( pv_virial ) ) pv_virial = virial % pv_virial
    IF ( PRESENT ( pv_xc ) ) pv_xc = virial % pv_xc
    IF ( PRESENT ( pv_fock_4c ) ) pv_fock_4c = virial % pv_fock_4c
    IF ( PRESENT ( pv_constraint ) ) pv_constraint = virial % pv_constraint
    IF ( PRESENT ( pv_kinetic ) ) pv_kinetic = virial % pv_kinetic
    IF ( PRESENT ( pv_availability ) ) pv_availability = virial % pv_availability
    IF ( PRESENT ( pv_calculate ) ) pv_calculate = virial % pv_calculate
    IF ( PRESENT ( pv_numer ) ) pv_numer = virial % pv_numer
    IF ( PRESENT ( pv_diagonal ) ) pv_diagonal = virial % pv_diagonal

  END SUBROUTINE virial_get

! *****************************************************************************
  SUBROUTINE virial_create ( virial, error )

    TYPE(virial_type), POINTER               :: virial
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'virial_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat

    ALLOCATE ( virial, stat=istat )
    CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)

    CALL zero_virial ( virial )
    last_virial_id_nr=last_virial_id_nr+1
    virial%id_nr=last_virial_id_nr
    virial%ref_count=1
  END SUBROUTINE virial_create

! *****************************************************************************
!> \brief retains the given virial_type
!> \param virial_type the virial_type to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      see doc/ReferenceCounting.html
!> \par History
!>      04.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE virial_retain ( virial, error)
    TYPE(virial_type), POINTER               :: virial
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'virial_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(virial),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPreconditionNoFail(virial%ref_count>0,cp_failure_level,routineP,error)
       virial%ref_count=virial%ref_count+1
    END IF
  END SUBROUTINE virial_retain

! *****************************************************************************
!> \brief releases the given virial_type
!> \param virial the virial_type to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      see doc/ReferenceCounting.html
!> \par History
!>      04.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE virial_release(virial, error)
    TYPE(virial_type), POINTER               :: virial
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'virial_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat
    LOGICAL                                  :: failure

    failure=.FALSE.
    IF (ASSOCIATED(virial)) THEN
       CPPreconditionNoFail(virial%ref_count>0,cp_failure_level,routineP,error)
       virial%ref_count=virial%ref_count-1
       IF (virial%ref_count.EQ.0) THEN
          DEALLOCATE(virial,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ENDIF
       NULLIFY(virial)
    END IF
  END SUBROUTINE virial_release

END MODULE virial_types

