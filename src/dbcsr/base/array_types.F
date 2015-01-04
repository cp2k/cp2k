!r-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Array objects with reference counting.
!> \author  Urban Borstnik
!> \date    2009-05-12
!> \version 1.0
!>
!> <b>Modification history:</b>
!> - Created
! *****************************************************************************
MODULE array_types

  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER :: dp = KIND(0.0d0)

  PUBLIC :: array_i1d_obj
  PUBLIC :: array_new,&
            array_hold,&
            array_release,&
            array_nullify,&
            array_exists
  PUBLIC :: array_data,&
            array_size,&
            array_equality
  PUBLIC :: array_get


  INTERFACE array_new
     MODULE PROCEDURE array_new_i1d, array_new_i1d_lb, array_new_i2d, array_new_i2d_lb,&
                      array_new_d1d, array_new_d1d_lb, array_new_d2d, array_new_d2d_lb
  END INTERFACE
  INTERFACE array_hold
     MODULE PROCEDURE array_hold_i1d, array_hold_i2d,&
                      array_hold_d1d, array_hold_d2d
  END INTERFACE
  INTERFACE array_release
     MODULE PROCEDURE array_release_i1d, array_release_i2d,&
                      array_release_d1d, array_release_d2d
  END INTERFACE
  INTERFACE array_nullify
     MODULE PROCEDURE array_nullify_i1d, array_nullify_i2d,&
                      array_nullify_d1d, array_nullify_d2d
  END INTERFACE

  INTERFACE array_exists
     MODULE PROCEDURE array_exists_i1d, array_exists_i2d,&
                      array_exists_d1d, array_exists_d2d
  END INTERFACE

  INTERFACE array_data
     MODULE PROCEDURE array_data_i1d, array_data_i2d,&
                      array_data_d1d, array_data_d2d
  END INTERFACE

  INTERFACE array_size
     MODULE PROCEDURE array_size_i1d, array_size_i2d,&
                      array_size_d1d, array_size_d2d
  END INTERFACE

  INTERFACE array_equality
     MODULE PROCEDURE array_equality_i1d, array_equality_i2d,&
                      array_equality_d1d, array_equality_d2d
  END INTERFACE

  INTERFACE array_get
     MODULE PROCEDURE array_get_i1d,&
                      array_get_d1d
  END INTERFACE


  TYPE array_i1d_type
     INTEGER, DIMENSION(:), POINTER     :: DATA
     INTEGER                            :: refcount
  END TYPE array_i1d_type
  TYPE array_i1d_obj
     TYPE(array_i1d_type), POINTER      :: low
  END TYPE array_i1d_obj

  TYPE array_i2d_type
     INTEGER, DIMENSION(:,:), POINTER   :: DATA
     INTEGER                            :: refcount
  END TYPE array_i2d_type
  TYPE array_i2d_obj
     TYPE(array_i2d_type), POINTER      :: low
  END TYPE array_i2d_obj

  TYPE array_d1d_type
     REAL(KIND=dp), DIMENSION(:), POINTER    :: DATA
     INTEGER                                 :: refcount
  END TYPE array_d1d_type
  TYPE array_d1d_obj
     TYPE(array_d1d_type), POINTER           :: low
  END TYPE array_d1d_obj

  TYPE array_d2d_type
     REAL(KIND=dp), DIMENSION(:,:), POINTER  :: DATA
     INTEGER                                 :: refcount
  END TYPE array_d2d_type
  TYPE array_d2d_obj
     TYPE(array_d2d_type), POINTER           :: low
  END TYPE array_d2d_obj

  !

CONTAINS

  !

! *****************************************************************************
!> \brief ...
!> \param array ...
!> \param DATA ...
!> \param gift ...
! *****************************************************************************
  SUBROUTINE array_new_i1d(array, DATA, gift)
    TYPE(array_i1d_obj), INTENT(OUT)         :: array
    INTEGER, DIMENSION(:), POINTER           :: DATA
    LOGICAL, INTENT(IN), OPTIONAL            :: gift

    INTEGER                                  :: lb, ub
    LOGICAL                                  :: g

    ALLOCATE (array%low)
    array%low%refcount = 1
    g = .FALSE.
    IF (PRESENT (gift)) g = gift
    IF (g) THEN
       array%low%data => DATA
    ELSE
       lb = LBOUND(DATA, 1)
       ub = UBOUND(DATA, 1)
       ALLOCATE (array%low%data(lb:ub))
       array%low%data(:) = DATA(:)
    ENDIF
  END SUBROUTINE array_new_i1d
! *****************************************************************************
!> \brief ...
!> \param array ...
!> \param DATA ...
!> \param lb ...
! *****************************************************************************
  SUBROUTINE array_new_i1d_lb(array, DATA, lb)
    TYPE(array_i1d_obj), INTENT(OUT)         :: array
    INTEGER, DIMENSION(:), INTENT(IN)        :: DATA
    INTEGER, INTENT(IN)                      :: lb

    INTEGER                                  :: ub

    ALLOCATE (array%low)
    array%low%refcount = 1
    ub = lb + SIZE(DATA) - 1
    ALLOCATE (array%low%data(lb:ub))
    array%low%data(:) = DATA(:)
  END SUBROUTINE array_new_i1d_lb
! *****************************************************************************
!> \brief ...
!> \param array ...
!> \param DATA ...
!> \param gift ...
! *****************************************************************************
  SUBROUTINE array_new_i2d(array, DATA, gift)
    TYPE(array_i2d_obj), INTENT(OUT)         :: array
    INTEGER, DIMENSION(:, :), POINTER        :: DATA
    LOGICAL, INTENT(IN), OPTIONAL            :: gift

    INTEGER, DIMENSION(2)                    :: lb, ub
    LOGICAL                                  :: g

    ALLOCATE (array%low)
    array%low%refcount = 1
    g = .FALSE.
    IF (PRESENT (gift)) g = gift
    IF (g) THEN
       array%low%data => DATA
    ELSE
       lb = LBOUND(DATA)
       ub = UBOUND(DATA)
       ALLOCATE (array%low%data(lb(1):ub(1), lb(2):ub(2)))
       array%low%data(:,:) = DATA(:,:)
    ENDIF
  END SUBROUTINE array_new_i2d
! *****************************************************************************
!> \brief ...
!> \param array ...
!> \param DATA ...
!> \param lb ...
! *****************************************************************************
  SUBROUTINE array_new_i2d_lb(array, DATA, lb)
    TYPE(array_i2d_obj), INTENT(OUT)         :: array
    INTEGER, DIMENSION(:, :), INTENT(IN)     :: DATA
    INTEGER, DIMENSION(2), INTENT(IN)        :: lb

    INTEGER, DIMENSION(2)                    :: ub

    ALLOCATE (array%low)
    array%low%refcount = 1
    ub(1) = lb(1) + SIZE(DATA,1) - 1
    ub(2) = lb(2) + SIZE(DATA,2) - 1
    ALLOCATE (array%low%data(lb(1):ub(1), lb(2):ub(2)))
    array%low%data(:,:) = DATA(:,:)
  END SUBROUTINE array_new_i2d_lb
! *****************************************************************************
!> \brief ...
!> \param array ...
!> \param DATA ...
!> \param gift ...
! *****************************************************************************
  SUBROUTINE array_new_d1d(array, DATA, gift)
    TYPE(array_d1d_obj), INTENT(OUT)         :: array
    REAL(KIND=dp), DIMENSION(:), POINTER     :: DATA
    LOGICAL, INTENT(IN), OPTIONAL            :: gift

    INTEGER                                  :: lb, ub
    LOGICAL                                  :: g

    ALLOCATE (array%low)
    array%low%refcount = 1
    g = .FALSE.
    IF (PRESENT (gift)) g = gift
    IF (g) THEN
       array%low%data => DATA
    ELSE
       lb = LBOUND(DATA, 1)
       ub = UBOUND(DATA, 1)
       ALLOCATE (array%low%data(lb:ub))
       array%low%data(:) = DATA(:)
    ENDIF
  END SUBROUTINE array_new_d1d
! *****************************************************************************
!> \brief ...
!> \param array ...
!> \param DATA ...
!> \param lb ...
! *****************************************************************************
  SUBROUTINE array_new_d1d_lb(array, DATA, lb)
    TYPE(array_d1d_obj), INTENT(OUT)         :: array
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: DATA
    INTEGER, INTENT(IN)                      :: lb

    INTEGER                                  :: ub

    ALLOCATE (array%low)
    array%low%refcount = 1
    ub = lb + SIZE(DATA) - 1
    ALLOCATE (array%low%data(lb:ub))
    array%low%data(:) = DATA(:)
  END SUBROUTINE array_new_d1d_lb
! *****************************************************************************
!> \brief ...
!> \param array ...
!> \param DATA ...
!> \param gift ...
! *****************************************************************************
  SUBROUTINE array_new_d2d(array, DATA, gift)
    TYPE(array_d2d_obj), INTENT(OUT)         :: array
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: DATA
    LOGICAL, INTENT(IN), OPTIONAL            :: gift

    INTEGER, DIMENSION(2)                    :: lb, ub
    LOGICAL                                  :: g

    ALLOCATE (array%low)
    array%low%refcount = 1
    g = .FALSE.
    IF (PRESENT (gift)) g = gift
    IF (g) THEN
       array%low%data => DATA
    ELSE
       lb = LBOUND(DATA)
       ub = UBOUND(DATA)
       ALLOCATE (array%low%data(lb(1):ub(1), lb(2):ub(2)))
       array%low%data(:,:) = DATA(:,:)
    ENDIF
  END SUBROUTINE array_new_d2d
! *****************************************************************************
!> \brief ...
!> \param array ...
!> \param DATA ...
!> \param lb ...
! *****************************************************************************
  SUBROUTINE array_new_d2d_lb(array, DATA, lb)
    TYPE(array_d2d_obj), INTENT(OUT)         :: array
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: DATA
    INTEGER, DIMENSION(2), INTENT(IN)        :: lb

    INTEGER, DIMENSION(2)                    :: ub

    ALLOCATE (array%low)
    array%low%refcount = 1
    ub(1) = lb(1) + SIZE(DATA,1) - 1
    ub(2) = lb(2) + SIZE(DATA,2) - 1
    ALLOCATE (array%low%data(lb(1):ub(1), lb(2):ub(2)))
    array%low%data(:,:) = DATA(:,:)
  END SUBROUTINE array_new_d2d_lb

  !

! *****************************************************************************
!> \brief ...
!> \param array ...
! *****************************************************************************
  PURE SUBROUTINE array_hold_i1d (array)
    TYPE(array_i1d_obj), INTENT(INOUT)       :: array

    array%low%refcount = array%low%refcount + 1
  END SUBROUTINE array_hold_i1d
! *****************************************************************************
!> \brief ...
!> \param array ...
! *****************************************************************************
  PURE SUBROUTINE array_hold_i2d (array)
    TYPE(array_i2d_obj), INTENT(INOUT)       :: array

    array%low%refcount = array%low%refcount + 1
  END SUBROUTINE array_hold_i2d
! *****************************************************************************
!> \brief ...
!> \param array ...
! *****************************************************************************
  PURE SUBROUTINE array_hold_d1d (array)
    TYPE(array_d1d_obj), INTENT(INOUT)       :: array

    array%low%refcount = array%low%refcount + 1
  END SUBROUTINE array_hold_d1d
! *****************************************************************************
!> \brief ...
!> \param array ...
! *****************************************************************************
  PURE SUBROUTINE array_hold_d2d (array)
    TYPE(array_d2d_obj), INTENT(INOUT)       :: array

    array%low%refcount = array%low%refcount + 1
  END SUBROUTINE array_hold_d2d

  !

! *****************************************************************************
!> \brief ...
!> \param array ...
! *****************************************************************************
  SUBROUTINE array_release_i1d (array)
    TYPE(array_i1d_obj), INTENT(INOUT)       :: array

    IF (ASSOCIATED (array%low)) THEN
       array%low%refcount = array%low%refcount - 1
       IF (array%low%refcount .EQ. 0) THEN
          DEALLOCATE(array%low%data)
          DEALLOCATE(array%low)
          NULLIFY (array%low)
       ENDIF
    ENDIF
  END SUBROUTINE array_release_i1d
! *****************************************************************************
!> \brief ...
!> \param array ...
! *****************************************************************************
  SUBROUTINE array_release_i2d (array)
    TYPE(array_i2d_obj), INTENT(INOUT)       :: array

    IF (ASSOCIATED (array%low)) THEN
       array%low%refcount = array%low%refcount - 1
       IF (array%low%refcount .EQ. 0 .AND. ASSOCIATED (array%low)) THEN
          DEALLOCATE(array%low%data)
          DEALLOCATE(array%low)
          NULLIFY (array%low)
       ENDIF
    ENDIF
  END SUBROUTINE array_release_i2d
! *****************************************************************************
!> \brief ...
!> \param array ...
! *****************************************************************************
  SUBROUTINE array_release_d1d (array)
    TYPE(array_d1d_obj), INTENT(INOUT)       :: array

    IF (ASSOCIATED (array%low)) THEN
       array%low%refcount = array%low%refcount - 1
       IF (array%low%refcount .EQ. 0 .AND. ASSOCIATED (array%low)) THEN
          DEALLOCATE(array%low%data)
          DEALLOCATE(array%low)
          NULLIFY (array%low)
       ENDIF
    ENDIF
  END SUBROUTINE array_release_d1d
! *****************************************************************************
!> \brief ...
!> \param array ...
! *****************************************************************************
  SUBROUTINE array_release_d2d (array)
    TYPE(array_d2d_obj), INTENT(INOUT)       :: array

    IF (ASSOCIATED (array%low)) THEN
       array%low%refcount = array%low%refcount - 1
       IF (array%low%refcount .EQ. 0 .AND. ASSOCIATED (array%low)) THEN
          DEALLOCATE(array%low%data)
          DEALLOCATE(array%low)
          NULLIFY (array%low)
       ENDIF
    ENDIF
  END SUBROUTINE array_release_d2d

  !

! *****************************************************************************
!> \brief ...
!> \param array ...
! *****************************************************************************
  PURE SUBROUTINE array_nullify_i1d (array)
    TYPE(array_i1d_obj), INTENT(INOUT)       :: array

    NULLIFY (array%low)
  END SUBROUTINE array_nullify_i1d
! *****************************************************************************
!> \brief ...
!> \param array ...
! *****************************************************************************
  PURE SUBROUTINE array_nullify_i2d (array)
    TYPE(array_i2d_obj), INTENT(INOUT)       :: array

    NULLIFY (array%low)
  END SUBROUTINE array_nullify_i2d
! *****************************************************************************
!> \brief ...
!> \param array ...
! *****************************************************************************
  PURE SUBROUTINE array_nullify_d1d (array)
    TYPE(array_d1d_obj), INTENT(INOUT)       :: array

    NULLIFY (array%low)
  END SUBROUTINE array_nullify_d1d
! *****************************************************************************
!> \brief ...
!> \param array ...
! *****************************************************************************
  PURE SUBROUTINE array_nullify_d2d (array)
    TYPE(array_d2d_obj), INTENT(INOUT)       :: array

    NULLIFY (array%low)
  END SUBROUTINE array_nullify_d2d

  !

! *****************************************************************************
!> \brief ...
!> \param array ...
!> \retval array_exists ...
! *****************************************************************************
  PURE FUNCTION array_exists_i1d (array) RESULT (array_exists)
    TYPE(array_i1d_obj), INTENT(IN)          :: array
    LOGICAL                                  :: array_exists

    array_exists = ASSOCIATED (array%low)
    IF (array_exists) array_exists = array%low%refcount .GT. 0
  END FUNCTION array_exists_i1d
! *****************************************************************************
!> \brief ...
!> \param array ...
!> \retval array_exists ...
! *****************************************************************************
  PURE FUNCTION array_exists_i2d (array) RESULT (array_exists)
    TYPE(array_i2d_obj), INTENT(IN)          :: array
    LOGICAL                                  :: array_exists

    array_exists = ASSOCIATED (array%low)
    IF (array_exists) array_exists = array%low%refcount .GT. 0
  END FUNCTION array_exists_i2d
! *****************************************************************************
!> \brief ...
!> \param array ...
!> \retval array_exists ...
! *****************************************************************************
  PURE FUNCTION array_exists_d1d (array) RESULT (array_exists)
    TYPE(array_d1d_obj), INTENT(IN)          :: array
    LOGICAL                                  :: array_exists

    array_exists = ASSOCIATED (array%low)
    IF (array_exists) array_exists = array%low%refcount .GT. 0

  END FUNCTION array_exists_d1d
! *****************************************************************************
!> \brief ...
!> \param array ...
!> \retval array_exists ...
! *****************************************************************************
  PURE FUNCTION array_exists_d2d (array) RESULT (array_exists)
    TYPE(array_d2d_obj), INTENT(IN)          :: array
    LOGICAL                                  :: array_exists

    array_exists = ASSOCIATED (array%low)
    IF (array_exists) array_exists = array%low%refcount .GT. 0

  END FUNCTION array_exists_d2d


  !

! *****************************************************************************
!> \brief ...
!> \param array ...
!> \retval DATA ...
! *****************************************************************************
  FUNCTION array_data_i1d(array) RESULT (DATA)
    TYPE(array_i1d_obj), INTENT(IN)          :: array
    INTEGER, DIMENSION(:), POINTER           :: DATA

    IF (ASSOCIATED (array%low)) THEN
       DATA => array%low%data
    ELSE
       NULLIFY (DATA)
    ENDIF
  END FUNCTION array_data_i1d
! *****************************************************************************
!> \brief ...
!> \param array ...
!> \retval DATA ...
! *****************************************************************************
  FUNCTION array_data_i2d(array) RESULT (DATA)
    TYPE(array_i2d_obj), INTENT(IN)          :: array
    INTEGER, DIMENSION(:, :), POINTER        :: DATA

    IF (ASSOCIATED (array%low)) THEN
       DATA => array%low%data
    ELSE
       NULLIFY (DATA)
    ENDIF
  END FUNCTION array_data_i2d
! *****************************************************************************
!> \brief ...
!> \param array ...
!> \retval DATA ...
! *****************************************************************************
  FUNCTION array_data_d1d(array) RESULT (DATA)
    TYPE(array_d1d_obj), INTENT(IN)          :: array
    REAL(KIND=dp), DIMENSION(:), POINTER     :: DATA

    IF (ASSOCIATED (array%low)) THEN
       DATA => array%low%data
    ELSE
       NULLIFY (DATA)
    ENDIF
  END FUNCTION array_data_d1d
! *****************************************************************************
!> \brief ...
!> \param array ...
!> \retval DATA ...
! *****************************************************************************
  FUNCTION array_data_d2d(array) RESULT (DATA)
    TYPE(array_d2d_obj), INTENT(IN)          :: array
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: DATA

    IF (ASSOCIATED (array%low)) THEN
       DATA => array%low%data
    ELSE
       NULLIFY (DATA)
    ENDIF
  END FUNCTION array_data_d2d

  !

! *****************************************************************************
!> \brief ...
!> \param array ...
!> \retval the_size ...
! *****************************************************************************
  PURE FUNCTION array_size_i1d(array) RESULT (the_size)
    TYPE(array_i1d_obj), INTENT(IN)          :: array
    INTEGER                                  :: the_size

    IF (ASSOCIATED (array%low)) THEN
       the_size = SIZE(array%low%data)
    ELSE
       the_size = 0
    ENDIF
  END FUNCTION array_size_i1d
! *****************************************************************************
!> \brief ...
!> \param array ...
!> \retval the_size ...
! *****************************************************************************
  PURE FUNCTION array_size_i2d(array) RESULT (the_size)
    TYPE(array_i2d_obj), INTENT(IN)          :: array
    INTEGER                                  :: the_size

    IF (ASSOCIATED (array%low)) THEN
       the_size = SIZE(array%low%data)
    ELSE
       the_size = 0
    ENDIF
  END FUNCTION array_size_i2d
! *****************************************************************************
!> \brief ...
!> \param array ...
!> \retval the_size ...
! *****************************************************************************
  PURE FUNCTION array_size_d1d(array) RESULT (the_size)
    TYPE(array_d1d_obj), INTENT(IN)          :: array
    INTEGER                                  :: the_size

    IF (ASSOCIATED (array%low)) THEN
       the_size = SIZE(array%low%data)
    ELSE
       the_size = 0
    ENDIF
  END FUNCTION array_size_d1d
! *****************************************************************************
!> \brief ...
!> \param array ...
!> \retval the_size ...
! *****************************************************************************
  PURE FUNCTION array_size_d2d(array) RESULT (the_size)
    TYPE(array_d2d_obj), INTENT(IN)          :: array
    INTEGER                                  :: the_size

    IF (ASSOCIATED (array%low)) THEN
       the_size = SIZE(array%low%data)
    ELSE
       the_size = 0
    ENDIF
  END FUNCTION array_size_d2d

  !

! *****************************************************************************
!> \brief ...
!> \param array1 ...
!> \param array2 ...
!> \retval are_equal ...
! *****************************************************************************
  PURE FUNCTION array_equality_i1d(array1, array2) RESULT (are_equal)
    TYPE(array_i1d_obj), INTENT(IN)          :: array1, array2
    LOGICAL                                  :: are_equal

    are_equal = .FALSE.
    IF (ASSOCIATED (array1%low) .AND. ASSOCIATED (array2%low)) THEN
       IF (SIZE(array1%low%data) .NE. SIZE(array2%low%data)) RETURN
       are_equal = ALL(array1%low%data .EQ. array2%low%data)
    ENDIF
  END FUNCTION array_equality_i1d
! *****************************************************************************
!> \brief ...
!> \param array1 ...
!> \param array2 ...
!> \retval are_equal ...
! *****************************************************************************
  PURE FUNCTION array_equality_i2d(array1, array2) RESULT (are_equal)
    TYPE(array_i2d_obj), INTENT(IN)          :: array1, array2
    LOGICAL                                  :: are_equal

    are_equal = .FALSE.
    IF (ASSOCIATED (array1%low) .AND. ASSOCIATED (array2%low)) THEN
       IF (SIZE(array1%low%data) .NE. SIZE(array2%low%data)) RETURN
       are_equal = ALL(array1%low%data .EQ. array2%low%data)
    ENDIF
  END FUNCTION array_equality_i2d
! *****************************************************************************
!> \brief ...
!> \param array1 ...
!> \param array2 ...
!> \retval are_equal ...
! *****************************************************************************
  PURE FUNCTION array_equality_d1d(array1, array2) RESULT (are_equal)
    TYPE(array_d1d_obj), INTENT(IN)          :: array1, array2
    LOGICAL                                  :: are_equal

    are_equal = .FALSE.
    IF (ASSOCIATED (array1%low) .AND. ASSOCIATED (array2%low)) THEN
       IF (SIZE(array1%low%data) .NE. SIZE(array2%low%data)) RETURN
       are_equal = ALL(array1%low%data .EQ. array2%low%data)
    ENDIF
  END FUNCTION array_equality_d1d
! *****************************************************************************
!> \brief ...
!> \param array1 ...
!> \param array2 ...
!> \retval are_equal ...
! *****************************************************************************
  PURE FUNCTION array_equality_d2d(array1, array2) RESULT (are_equal)
    TYPE(array_d2d_obj), INTENT(IN)          :: array1, array2
    LOGICAL                                  :: are_equal

    are_equal = .FALSE.
    IF (ASSOCIATED (array1%low) .AND. ASSOCIATED (array2%low)) THEN
       IF (SIZE(array1%low%data) .NE. SIZE(array2%low%data)) RETURN
       are_equal = ALL(array1%low%data .EQ. array2%low%data)
    ENDIF
  END FUNCTION array_equality_d2d

  !

! *****************************************************************************
!> \brief ...
!> \param array ...
!> \param index1 ...
!> \retval value ...
! *****************************************************************************
  PURE FUNCTION array_get_i1d(array, index1) RESULT (value)
    TYPE(array_i1d_obj), INTENT(IN)          :: array
    INTEGER, INTENT(IN)                      :: index1
    INTEGER                                  :: value

    value = array%low%data(index1)
  END FUNCTION array_get_i1d
! *****************************************************************************
!> \brief ...
!> \param array ...
!> \param index1 ...
!> \retval value ...
! *****************************************************************************
  PURE FUNCTION array_get_d1d(array, index1) RESULT (value)
    TYPE(array_d1d_obj), INTENT(IN)          :: array
    INTEGER, INTENT(IN)                      :: index1
    REAL(KIND=dp)                            :: value

    value = array%low%data(index1)
  END FUNCTION array_get_d1d

END MODULE array_types
