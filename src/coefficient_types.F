!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****s* cp2k/coefficient_types [1.0] *
!!
!!   NAME
!!     coefficient_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     apsi
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE coefficient_types
  USE ao_types,                        ONLY: ao_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             COMPLEXDATA3D,&
                                             NOSPACE,&
                                             REALDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_add,&
                                             pw_allocate,&
                                             pw_copy,&
                                             pw_fft_wrap,&
                                             pw_reduce,&
                                             pw_release,&
                                             pw_subtract,&
                                             pw_sumup,&
                                             pw_type,&
                                             pw_zero
  USE termination,                     ONLY: stop_program

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: coeff_type, coeff_allocate, coeff_deallocate, coeff_zero, coeff_init
  PUBLIC :: coeff_copy, coeff_add, coeff_sumup, coeff_subtract, coeff_reduce
  PUBLIC :: coeff_transform_space, coeff_array_type
  PUBLIC :: PLANEWAVES, ATOMICORBITALS
  PUBLIC :: PW_REALDATA1D, PW_COMPLEXDATA1D, PW_REALDATA3D, PW_COMPLEXDATA3D
  PUBLIC :: PW_NOSPACE, PW_REALSPACE, PW_RECIPROCALSPACE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN='coefficent_types'
  
  ! 'group_id' to grid?
  TYPE coeff_type
     INTEGER :: in_use
     TYPE ( pw_type ), POINTER :: pw
     TYPE ( ao_type ) :: ao  
  END TYPE coeff_type

  TYPE coeff_array_type
     TYPE(coeff_type), DIMENSION(:), POINTER :: coeff_array
  END TYPE coeff_array_type

  ! Flags for the structure member 'in_use'
  INTEGER, PARAMETER :: PLANEWAVES = 401, ATOMICORBITALS = 402, NOTHING=0

  INTERFACE coeff_allocate
     MODULE PROCEDURE &
          coeff_allocate_from_coeff, coeff_allocate_from_grid
  END INTERFACE

  INTERFACE coeff_transform_space
     MODULE PROCEDURE &
          coeff_transform_c1, coeff_transform_c1c2, coeff_transform_c1pw2, &
          coeff_transform_pw1c2
  END INTERFACE

  ! Flags for 'integral'
  PUBLIC :: SQUARE, SQUAREROOT
  INTEGER, PARAMETER :: SQUARE = 391, SQUAREROOT = 392

  INTEGER, PARAMETER :: PW_REALDATA1D = REALDATA1D
  INTEGER, PARAMETER :: PW_COMPLEXDATA1D = COMPLEXDATA1D
  INTEGER, PARAMETER :: PW_REALDATA3D = REALDATA3D
  INTEGER, PARAMETER :: PW_COMPLEXDATA3D = COMPLEXDATA3D

  INTEGER, PARAMETER :: PW_NOSPACE = NOSPACE
  INTEGER, PARAMETER :: PW_REALSPACE = REALSPACE
  INTEGER, PARAMETER :: PW_RECIPROCALSPACE = RECIPROCALSPACE

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************

!!****f* coefficient_types/coeff_init [1.0] *
!!
!!   NAME
!!     coeff_init
!!
!!   FUNCTION
!!     initializes a coefficent
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - coeff: the coefficent to initialize
!!     - pw: optional, a pw to initialize the coefficent
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE coeff_init(coeff, pw, error)
    TYPE(coeff_type), INTENT(inout)          :: coeff
    TYPE(pw_type), OPTIONAL, POINTER         :: pw
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'coeff_init', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  coeff%in_use=NOTHING
  NULLIFY(coeff%pw)
  IF (PRESENT(pw)) THEN
     IF (ASSOCIATED(pw)) THEN
        coeff%in_use=PLANEWAVES
        coeff%pw => pw
     END IF
  END IF
END SUBROUTINE coeff_init
!***************************************************************************

SUBROUTINE coeff_allocate_from_coeff ( c_out, c_model, use_basis )


    TYPE(coeff_type), INTENT(OUT)            :: c_out
    TYPE(coeff_type), INTENT(IN), TARGET     :: c_model
    INTEGER, INTENT(IN), OPTIONAL            :: use_basis

    INTEGER                                  :: stat, use_data

!------------------------------------------------------------------------------

  NULLIFY(c_out%pw)
  IF ( PRESENT ( use_basis ) ) THEN
     IF ( use_basis == PLANEWAVES ) THEN
        ALLOCATE(c_out%pw,stat=stat)
        IF (stat /=0 ) CALL stop_program(&
             "coefficient_types:coeff_allocate_from_coeff",&
             "c_out%pw allocation")
        c_out % in_use = PLANEWAVES
        CALL pw_allocate ( c_out % pw, c_model % pw % pw_grid, &
             use_data = c_model % pw % in_use )

     ELSE IF ( use_basis == PW_REALDATA1D .OR. use_basis == PW_COMPLEXDATA1D &
          .OR. use_basis == PW_REALDATA3D .OR. use_basis == PW_COMPLEXDATA3D &
          ) THEN

        c_out % in_use = PLANEWAVES

        CALL pw_allocate ( c_out % pw, c_model % pw % pw_grid, &
             use_data = use_basis )

     ELSE
        CALL stop_program ( "coeff_allocate_from_coeff", &
             "no suitable data with use_basis" )

     END IF

  ELSE
     use_data = c_model % pw % in_use
     IF ( use_data == PW_REALDATA1D .OR. use_data == PW_COMPLEXDATA1D .OR. &
          use_data == PW_REALDATA3D .OR. use_data == PW_COMPLEXDATA3D ) THEN

        c_out % in_use = PLANEWAVES
        CALL pw_allocate ( c_out % pw, c_model % pw % pw_grid, &
             use_data = c_model % pw % in_use )

     ELSE
        CALL stop_program ( "coeff_allocate_from_coeff", &
             "no suitable data without use_basis" )
     END IF
  END IF

END SUBROUTINE coeff_allocate_from_coeff

!******************************************************************************

SUBROUTINE coeff_allocate_from_grid ( c_out, pw_grid, use_data )


    TYPE(coeff_type), INTENT(OUT)            :: c_out
    TYPE(pw_grid_type), INTENT(IN), OPTIONAL :: pw_grid
    INTEGER, INTENT(IN)                      :: use_data

    INTEGER                                  :: stat

!------------------------------------------------------------------------------

  ALLOCATE(c_out%pw,stat=stat)
  IF (stat /=0 ) CALL stop_program(&
       "coefficient_types:coeff_allocate_from_grid",&
       "c_out%pw allocation")  
  IF ( PRESENT ( pw_grid ) ) THEN
     c_out % in_use = PLANEWAVES
     CALL pw_allocate ( c_out % pw, pw_grid, use_data )
  ELSE
     CALL stop_program ( "coeff_allocate_from_grid", "no suitable grid" )
  END IF

END SUBROUTINE coeff_allocate_from_grid

!******************************************************************************

SUBROUTINE coeff_deallocate ( coeff )

    TYPE(coeff_type), INTENT(INOUT)          :: coeff

!------------------------------------------------------------------------------

  IF ( coeff % in_use == PLANEWAVES ) THEN
     CALL pw_release ( coeff % pw )
     coeff%in_use=NOTHING
  ELSEIF (coeff%in_use /= NOTHING) THEN
     CALL stop_program ( "coeff_deallocate", "no valid data type" )
  END IF

END SUBROUTINE coeff_deallocate

!******************************************************************************

SUBROUTINE coeff_zero ( coeff )
  
    TYPE(coeff_type), INTENT(INOUT)          :: coeff

!------------------------------------------------------------------------------

  IF ( coeff % in_use == PLANEWAVES ) THEN
     CALL pw_zero ( coeff % pw )
  ELSE
     CALL stop_program ( "coeff_zero", "no valid data type" )
  END IF
  
END SUBROUTINE coeff_zero

!******************************************************************************

SUBROUTINE coeff_copy ( coeff1, coeff2 )
  
    TYPE(coeff_type), INTENT(IN)             :: coeff1
    TYPE(coeff_type), INTENT(INOUT)          :: coeff2

!------------------------------------------------------------------------------

  IF ( coeff1 % in_use == PLANEWAVES .AND. coeff2 % in_use == PLANEWAVES ) THEN
     CALL pw_copy ( coeff1 % pw, coeff2 % pw )
  ELSE
     CALL stop_program ( "coeff_copy", "no valid data type" )
  END IF
  
END SUBROUTINE coeff_copy

!******************************************************************************

SUBROUTINE coeff_add ( coeff1, coeff2, coeff3 )
  
    TYPE(coeff_type), INTENT(IN)             :: coeff1, coeff2
    TYPE(coeff_type), INTENT(INOUT)          :: coeff3

!------------------------------------------------------------------------------

  IF ( coeff1 % in_use == PLANEWAVES .AND. coeff2 % in_use == PLANEWAVES &
       .AND. coeff3 % in_use == PLANEWAVES ) THEN
     CALL pw_add ( coeff1 % pw, coeff2 % pw, coeff3 % pw )
  ELSE
     CALL stop_program ( "coeff_add", "no valid data type" )
  END IF
  
END SUBROUTINE coeff_add

!******************************************************************************

SUBROUTINE coeff_sumup ( coeff1, coeff2 )
  
    TYPE(coeff_type), INTENT(IN)             :: coeff1
    TYPE(coeff_type), INTENT(INOUT)          :: coeff2

!------------------------------------------------------------------------------

  IF ( coeff1 % in_use == PLANEWAVES .AND. coeff2 % in_use == PLANEWAVES ) THEN
     CALL pw_sumup ( coeff1 % pw, coeff2 % pw )
  ELSE
     CALL stop_program ( "coeff_sumup", "no valid data type" )
  END IF
  
END SUBROUTINE coeff_sumup

!******************************************************************************

SUBROUTINE coeff_subtract ( coeff1, coeff2, coeff3 )
  
    TYPE(coeff_type), INTENT(IN)             :: coeff1, coeff2
    TYPE(coeff_type), INTENT(INOUT)          :: coeff3

!------------------------------------------------------------------------------

  IF ( coeff1 % in_use == PLANEWAVES .AND. coeff2 % in_use == PLANEWAVES &
       .AND. coeff2 % in_use == PLANEWAVES ) THEN
     CALL pw_subtract ( coeff1 % pw, coeff2 % pw, coeff3 % pw )
  ELSE
     CALL stop_program ( "coeff_subtract", "no valid data type" )
  END IF
  
END SUBROUTINE coeff_subtract

!******************************************************************************

SUBROUTINE coeff_reduce ( coeff1, coeff2 )
  
    TYPE(coeff_type), INTENT(IN)             :: coeff1
    TYPE(coeff_type), INTENT(INOUT)          :: coeff2

!------------------------------------------------------------------------------

  IF ( coeff1 % in_use == PLANEWAVES .AND. coeff2 % in_use == PLANEWAVES ) THEN
     CALL pw_reduce ( coeff1 % pw, coeff2 % pw )
  ELSE
     CALL stop_program ( "coeff_reduce", "no valid data type" )
  END IF
  
END SUBROUTINE coeff_reduce

!******************************************************************************
!!****** coefficient_types/coeff_transform [1.0] *
!!
!!   NAME
!!     coeff_transform
!!
!!   FUNCTION
!!     Generic function for transformations of a coefficient_type
!!
!!   AUTHOR
!!     JGH (4-Jan-2001)
!!
!!   MODIFICATION HISTORY
!!     JGH (9-Feb-2001) added optional debug flag
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE coeff_transform_c1 ( c1, debug )
  

    TYPE(coeff_type), INTENT(INOUT), TARGET  :: c1
    LOGICAL, INTENT(IN), OPTIONAL            :: debug

!------------------------------------------------------------------------------
!..check if we are dealing with plane waves and call the correct routine

  IF ( c1 % in_use == PLANEWAVES ) THEN
    IF ( PRESENT ( debug ) ) THEN
      CALL pw_fft_wrap ( c1 % pw, debug )
    ELSE
      CALL pw_fft_wrap ( c1 % pw )
    ENDIF
  ELSE
    CALL stop_program ( "coeff_transform_c1", "not implemented")
  END IF

END SUBROUTINE coeff_transform_c1

!******************************************************************************

SUBROUTINE coeff_transform_c1c2 ( c1, c2, debug )
  

    TYPE(coeff_type), INTENT(IN)             :: c1
    TYPE(coeff_type), INTENT(INOUT)          :: c2
    LOGICAL, INTENT(IN), OPTIONAL            :: debug

!------------------------------------------------------------------------------
!..check if we are dealing with plane waves and call the correct routine

  IF ( c1 % in_use == PLANEWAVES  .AND. c2 % in_use == PLANEWAVES ) THEN
    CALL pw_fft_wrap ( c1 % pw, c2 % pw, debug )
  ELSE
    CALL stop_program ( "coeff_transform_c1c2", "not implemented")
  END IF

END SUBROUTINE coeff_transform_c1c2

!******************************************************************************

SUBROUTINE coeff_transform_c1pw2 ( c1, pw2, debug )
  

    TYPE(coeff_type), INTENT(IN)             :: c1
    TYPE(pw_type), INTENT(INOUT)             :: pw2
    LOGICAL, INTENT(IN), OPTIONAL            :: debug

!------------------------------------------------------------------------------
!..check if we are dealing with plane waves and call the correct routine

  IF ( c1 % in_use == PLANEWAVES ) THEN
    IF ( PRESENT ( debug ) ) THEN
      CALL pw_fft_wrap ( c1 % pw, pw2, debug )
    ELSE
      CALL pw_fft_wrap ( c1 % pw, pw2 )
    END IF
  ELSE
    CALL stop_program ( "coeff_transform_c1pw2", "not implemented" )
  END IF

END SUBROUTINE coeff_transform_c1pw2

!******************************************************************************

SUBROUTINE coeff_transform_pw1c2 ( pw1, c2, debug )
  

    TYPE(pw_type), INTENT(IN)                :: pw1
    TYPE(coeff_type), INTENT(INOUT)          :: c2
    LOGICAL, INTENT(IN), OPTIONAL            :: debug

!------------------------------------------------------------------------------
!..check if we are dealing with plane waves and call the correct routine

  IF ( c2 % in_use == PLANEWAVES ) THEN
    IF ( PRESENT ( debug ) ) THEN
      CALL pw_fft_wrap ( pw1, c2 % pw, debug )
    ELSE
      CALL pw_fft_wrap ( pw1, c2 % pw )
    END IF
  ELSE
    CALL stop_program ( "coeff_transform_pw1c2", "not implemented" )
  END IF

END SUBROUTINE coeff_transform_pw1c2

!!*****
!******************************************************************************

END MODULE coefficient_types

!******************************************************************************
