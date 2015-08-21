!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      - Container to hold basis sets 
!> \author JGH (09.07.2015)
! *****************************************************************************
MODULE basis_set_container_types

  USE basis_set_types,                 ONLY: deallocate_gto_basis_set,&
                                             gto_basis_set_type
  USE kinds,                           ONLY: default_string_length
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  ! Global parameters (only in this module)

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'basis_set_container_types'

! *****************************************************************************
  INTEGER, PARAMETER                       :: unknown_basis        = 100,&
                                              orbital_basis        = 101,&
                                              auxiliary_basis      = 102,&
                                              ri_aux_basis         = 103,&
                                              lri_aux_basis        = 104,&
                                              aux_fit_basis        = 105,&
                                              soft_basis           = 106,&
                                              hard_basis           = 107
! *****************************************************************************
  TYPE basis_set_container_type
     PRIVATE
     CHARACTER(LEN=default_string_length)       :: basis_type = ""
     INTEGER                                    :: basis_type_nr = 0
     TYPE(gto_basis_set_type), POINTER          :: basis_set => NULL()
  END TYPE basis_set_container_type
! *****************************************************************************

  PUBLIC :: basis_set_container_type

  PUBLIC :: remove_basis_set_container,&
            add_basis_set_to_container, get_basis_from_container,&
            remove_basis_from_container, translate_basis_type

! *****************************************************************************

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param basis ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE remove_basis_set_container(basis,error)
    TYPE(basis_set_container_type), &
      DIMENSION(:), INTENT(inout)            :: basis
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'remove_basis_set_container', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

    DO i=1,SIZE(basis)
       basis(i)%basis_type = ""
       basis(i)%basis_type_nr = 0
       IF(ASSOCIATED(basis(i)%basis_set)) THEN
          CALL deallocate_gto_basis_set(basis(i)%basis_set,error)
       END IF
    END DO

  END SUBROUTINE remove_basis_set_container

! *****************************************************************************
!> \brief ...
!> \param basis ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE init_basis_set_container(basis,error)
    TYPE(basis_set_container_type), &
      DIMENSION(:), INTENT(inout)            :: basis
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_basis_set_container', &
      routineP = moduleN//':'//routineN

    CALL remove_basis_set_container(basis,error) 

  END SUBROUTINE init_basis_set_container

! *****************************************************************************
!> \brief ...
!> \param basis_set_type ...
!> \retval basis_type_nr ...
! *****************************************************************************
  FUNCTION get_basis_type(basis_set_type) RESULT(basis_type_nr)
    CHARACTER(len=*)                         :: basis_set_type
    INTEGER                                  :: basis_type_nr

    SELECT CASE (basis_set_type)
       CASE ("ORB")
          basis_type_nr = orbital_basis
       CASE ("AUX")
          basis_type_nr = auxiliary_basis
       CASE ("RI_AUX")
          basis_type_nr = ri_aux_basis
       CASE ("LRI")
          basis_type_nr = lri_aux_basis
       CASE ("AUX_FIT")
          basis_type_nr = aux_fit_basis
       CASE ("SOFT")
          basis_type_nr = soft_basis
       CASE ("HARD")
          basis_type_nr = hard_basis
       CASE DEFAULT
          basis_type_nr = unknown_basis
    END SELECT

  END FUNCTION get_basis_type
! *****************************************************************************
!> \brief ...
!> \param basis_set_id ...
!> \retval basis_set_type ...
! *****************************************************************************
  FUNCTION translate_basis_type(basis_set_id) RESULT(basis_set_type)
    INTEGER                                  :: basis_set_id
    CHARACTER(len=default_string_length)     :: basis_set_type

    INTEGER, PARAMETER :: use_aux_basis_set = 3, use_aux_fit_basis_set = 2, &
      use_lri_basis_set = 5, use_orb_basis_set = 1, use_ri_aux_basis_set = 4

! This parameter identifies basis sets

    SELECT CASE (basis_set_id)
       CASE (use_orb_basis_set)
          basis_set_type = "ORB"
       CASE (use_ri_aux_basis_set)
          basis_set_type = "RI_AUX"
       CASE (use_aux_basis_set)
          basis_set_type = "AUX"
       CASE (use_aux_fit_basis_set)
          basis_set_type = "AUX_FIT"
       CASE (use_lri_basis_set)
          basis_set_type = "LRI"
       CASE DEFAULT
          basis_set_type = ""
    END SELECT

  END FUNCTION translate_basis_type

! *****************************************************************************
!> \brief ...
!> \param container ...
!> \param basis_set ...
!> \param basis_set_type ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE add_basis_set_to_container(container,basis_set,basis_set_type,error)
    TYPE(basis_set_container_type), &
      DIMENSION(:), INTENT(inout)            :: container
    TYPE(gto_basis_set_type), POINTER        :: basis_set
    CHARACTER(len=*)                         :: basis_set_type
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'add_basis_set_to_container', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure, success

    success = .FALSE.
    DO i=1,SIZE(container)
       IF(container(i)%basis_type_nr == 0) THEN
          container(i)%basis_type = basis_set_type
          container(i)%basis_set => basis_set
          container(i)%basis_type_nr = get_basis_type(basis_set_type)
          success = .TRUE.
          EXIT
       END IF
    END DO
    CPPostcondition(success,cp_failure_level,routineP,error,failure)

  END SUBROUTINE add_basis_set_to_container

! *****************************************************************************
!> \brief ...
!> \param container ...
!> \param inum ...
!> \param basis_type ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE remove_basis_from_container(container,inum,basis_type,error)
    TYPE(basis_set_container_type), &
      DIMENSION(:), INTENT(inout)            :: container
    INTEGER, INTENT(IN), OPTIONAL            :: inum
    CHARACTER(len=*), OPTIONAL               :: basis_type
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'remove_basis_from_container', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: basis_nr, i, ibas
    LOGICAL                                  :: failure

    IF(PRESENT(inum)) THEN
       CPPrecondition(inum <= SIZE(container),cp_failure_level,routineP,error,failure)
       CPPrecondition(inum >= 1,cp_failure_level,routineP,error,failure)
       ibas = inum
    ELSE IF(PRESENT(basis_type)) THEN
       basis_nr = get_basis_type(basis_type)
       ibas = 0
       DO i=1,SIZE(container)
          IF(container(i)%basis_type_nr == basis_nr) THEN
             ibas = i
             EXIT
          END IF
       END DO
    ELSE
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END IF
    !
    IF(ibas /= 0) THEN
       container(ibas)%basis_type = ""
       container(ibas)%basis_type_nr = 0
       IF(ASSOCIATED(container(ibas)%basis_set)) THEN
          CALL deallocate_gto_basis_set(container(ibas)%basis_set,error)
       END IF
    END IF
    ! shift other basis sets
    DO i=ibas+1,SIZE(container)
       IF(container(i)%basis_type_nr == 0) CYCLE
       container(i-1)%basis_type = container(i)%basis_type
       container(i-1)%basis_set => container(i)%basis_set
       container(i-1)%basis_type_nr = container(i)%basis_type_nr
    END DO

  END SUBROUTINE remove_basis_from_container

! *****************************************************************************
!> \brief Retrieve a basis set from the container
!> \param container ...
!> \param basis_set ...
!> \param inumbas ...
!> \param basis_type ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE get_basis_from_container(container,basis_set,inumbas,basis_type,error)
    TYPE(basis_set_container_type), &
      DIMENSION(:), INTENT(inout)            :: container
    TYPE(gto_basis_set_type), POINTER        :: basis_set
    INTEGER, OPTIONAL                        :: inumbas
    CHARACTER(len=*), OPTIONAL               :: basis_type
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_basis_from_container', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: basis_nr, i
    LOGICAL                                  :: failure

    IF(PRESENT(inumbas)) THEN
       CPPrecondition(inumbas <= SIZE(container),cp_failure_level,routineP,error,failure)
       CPPrecondition(inumbas >= 1,cp_failure_level,routineP,error,failure)
       basis_set => container(inumbas)%basis_set
       IF(PRESENT(basis_type)) THEN
          basis_type = container(inumbas)%basis_type
       END IF
    ELSE IF(PRESENT(basis_type)) THEN
       NULLIFY(basis_set)
       basis_nr = get_basis_type(basis_type)
       DO i=1,SIZE(container)
          IF(container(i)%basis_type_nr == basis_nr) THEN
             basis_set => container(i)%basis_set
             EXIT
          END IF
       END DO
    ELSE
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE get_basis_from_container
! *****************************************************************************

END MODULE basis_set_container_types
