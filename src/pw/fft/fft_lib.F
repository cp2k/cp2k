!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2006 - 2013  CP2K developers group                          !
!-----------------------------------------------------------------------------!
MODULE fft_lib

  USE fft_kinds,                       ONLY: dp
  USE fft_plan,                        ONLY: fft_plan_type
  USE fftsg_lib,                       ONLY: fftsg1dm,&
                                             fftsg3d,&
                                             fftsg_do_cleanup,&
                                             fftsg_get_lengths
  USE fftw3_lib,                       ONLY: fftw31dm,&
                                             fftw33d,&
                                             fftw3_create_plan_1dm,&
                                             fftw3_create_plan_3d,&
                                             fftw3_destroy_plan,&
                                             fftw3_do_cleanup,&
                                             fftw3_get_lengths

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'fft_lib'

  PUBLIC :: fft_do_cleanup, fft_get_lengths, fft_create_plan_3d, fft_create_plan_1dm, fft_1dm, fft_library, fft_3d, fft_destroy_plan 

CONTAINS
! *****************************************************************************
!> \brief Interface to FFT libraries
!> \par History
!>      IAB 09-Jan-2009 : Modified to use fft_plan_type
!>                        (c) The Numerical Algorithms Group (NAG) Ltd, 2009 on behalf of the HECToR project
!> \author JGH
! *****************************************************************************
 FUNCTION fft_library ( fftlib ) RESULT ( flib )

    CHARACTER(len=*), INTENT(IN)             :: fftlib
    INTEGER                                  :: flib

  SELECT CASE ( fftlib )
    CASE DEFAULT
      flib = -1
    CASE ( "FFTSG" )
      flib = 1
    CASE ( "FFTW3" )
      flib = 3
  END SELECT

END FUNCTION fft_library

SUBROUTINE fft_do_cleanup(fft_type)
    INTEGER, INTENT(IN)                      :: fft_type

  SELECT CASE ( fft_type )
    CASE DEFAULT
      STOP "fft_do_cleanup"
    CASE ( 1 )
      CALL fftsg_do_cleanup ( )
    CASE ( 3 )
      CALL fftw3_do_cleanup ( )
  END SELECT

END SUBROUTINE

! *****************************************************************************
SUBROUTINE fft_get_lengths ( fft_type, DATA, max_length )

    INTEGER, INTENT(IN)                      :: fft_type
    INTEGER, DIMENSION(*)                    :: DATA
    INTEGER, INTENT(INOUT)                   :: max_length

  SELECT CASE ( fft_type )
    CASE DEFAULT
      STOP "fft_get_lengths"
    CASE ( 1 )
      CALL fftsg_get_lengths ( DATA, max_length )
    CASE ( 3 )
      CALL fftw3_get_lengths ( DATA, max_length )
  END SELECT

END SUBROUTINE fft_get_lengths

! *****************************************************************************

SUBROUTINE fft_create_plan_3d ( plan, fft_type, fft_in_place, fsign, n, zin, zout, plan_style)

    TYPE(fft_plan_type), INTENT(INOUT)       :: plan
    INTEGER, INTENT(IN)                      :: fft_type
    LOGICAL, INTENT(IN)                      :: fft_in_place
    INTEGER, INTENT(IN)                      :: fsign
    INTEGER, DIMENSION(3), INTENT(IN)        :: n
    COMPLEX(KIND=dp), DIMENSION(*), &
      INTENT(INOUT)                          :: zin, zout
    INTEGER, INTENT(IN)                      :: plan_style

  plan % fft_type = fft_type
  plan % fsign = fsign
  plan % fft_in_place = fft_in_place
  plan % n_3d = n
!$ plan % need_alt_plan = .FALSE.

  ! Planning only needed for FFTW3
  IF (fft_type .EQ. 3) THEN 
    CALL fftw3_create_plan_3d ( plan, zin, zout, plan_style)
  END IF

  plan % valid = .TRUE.

END SUBROUTINE fft_create_plan_3d

!
! really ugly, plan is intent out, because plan%fsign is also a status flag
! if something goes wrong, plan%fsign is set to zero, and the plan becomes invalid
!
SUBROUTINE fft_3d ( plan, scale, zin, zout, stat )

  USE fft_plan,                           ONLY: fft_plan_type
  USE fft_kinds,                          ONLY: dp
    TYPE(fft_plan_type), INTENT(IN)          :: plan
    REAL(KIND=dp), INTENT(IN)                :: scale
    COMPLEX(KIND=dp), DIMENSION(*), &
      INTENT(INOUT)                          :: zin, zout
    INTEGER, INTENT(OUT)                     :: stat

  stat=plan%fsign
  IF ( plan%n_3d(1)*plan%n_3d(2)*plan%n_3d(3) > 0 ) THEN
    SELECT CASE ( plan%fft_type )
      CASE DEFAULT
        STOP "fft_3d"
      CASE ( 1 )
        CALL fftsg3d ( plan%fft_in_place, stat, scale, plan%n_3d, zin, zout )
      CASE ( 3 )
        CALL fftw33d ( plan, scale, zin, zout ,stat )
    END SELECT
  END IF
  ! stat is set to zero on error, -1,+1 are OK
  IF (stat.EQ.0) THEN
      stat=1
  ELSE
      stat=0
  ENDIF

END SUBROUTINE fft_3d

! *****************************************************************************

SUBROUTINE fft_create_plan_1dm ( plan, fft_type, fsign, trans, n, m, zin, zout, plan_style)

  USE fft_plan,                           ONLY: fft_plan_type
  USE fft_kinds,                          ONLY: dp
    TYPE(fft_plan_type), INTENT(INOUT)       :: plan
    INTEGER, INTENT(IN)                      :: fft_type, fsign
    LOGICAL, INTENT(IN)                      :: trans
    INTEGER, INTENT(IN)                      :: n, m
    COMPLEX(KIND=dp), DIMENSION(*), &
      INTENT(IN)                             :: zin, zout
    INTEGER, INTENT(IN)                      :: plan_style

  plan % fft_type = fft_type
  plan % fsign = fsign
  plan % trans = trans
  plan % n = n
  plan % m = m
!$ plan % need_alt_plan = .FALSE.

  ! Planning only needed for FFTW3
  IF (fft_type .EQ. 3) THEN
    CALL fftw3_create_plan_1dm ( plan, zin, zout, plan_style)
  END IF
  plan % valid = .TRUE.

END SUBROUTINE fft_create_plan_1dm

SUBROUTINE fft_destroy_plan ( plan )

  USE fft_plan,                           ONLY: fft_plan_type
  USE fft_kinds,                          ONLY: dp
    TYPE(fft_plan_type), INTENT(INOUT)       :: plan

! Planning only needed for FFTW3

  IF (plan % valid .AND. plan % fft_type .EQ. 3) THEN
    CALL fftw3_destroy_plan ( plan )
    plan % valid = .FALSE.
  END IF

END SUBROUTINE

SUBROUTINE fft_1dm ( plan, zin, zout, scale, stat)


    TYPE(fft_plan_type), INTENT(IN)          :: plan
    COMPLEX(KIND=dp), DIMENSION(*), &
      INTENT(INOUT)                          :: zin, zout
    REAL(KIND=dp), INTENT(IN)                :: scale
    INTEGER, INTENT(OUT)                     :: stat

  stat=plan%fsign
  IF ( plan%n * plan%m > 0 ) THEN
    SELECT CASE ( plan%fft_type )
      CASE DEFAULT
        STOP "fft_1dm"
      CASE ( 1 )
        CALL fftsg1dm ( stat, plan%trans, plan%n, plan%m, zin, zout, scale )
      CASE ( 3 )
        CALL fftw31dm ( plan, zin, zout, scale, stat )
    END SELECT
  END IF
  ! stat is set to zero on error, -1,+1 are OK
  IF (stat.EQ.0) THEN
      stat=1
  ELSE
      stat=0
  ENDIF

END SUBROUTINE fft_1dm

END MODULE

