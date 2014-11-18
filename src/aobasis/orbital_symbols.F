!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief orbital_symbols
!> \par History
!>      none
!> \author Matthias Krack (08.06.2000)
! *****************************************************************************
MODULE orbital_symbols
  
! Index:
! FUNCTION cgf_symbol(n,lxyz) RESULT(symbol)
! FUNCTION sgf_symbol(n,l,m) RESULT(symbol)
  USE termination,                     ONLY: stop_program
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'orbital_symbols'
  PRIVATE

  CHARACTER(LEN=1), PARAMETER, DIMENSION(0:11) :: l_sym = (/"s","p","d",&
                                                            "f","g","h",&
                                                            "i","j","k",&
                                                            "l","m","n"/)

  PUBLIC :: cgf_symbol,sgf_symbol
  PUBLIC :: l_sym

CONTAINS

! *****************************************************************************
!> \brief   Build a Cartesian orbital symbol (orbital labels for printing).
!> \param n ...
!> \param lxyz ...
!> \retval symbol ...
!> \date    07.07.99
!> \author  Matthias Krack
!> \version 1.0
! *****************************************************************************
  FUNCTION cgf_symbol(n,lxyz) RESULT(symbol)
    INTEGER, INTENT(IN)                      :: n
    INTEGER, DIMENSION(3), INTENT(IN)        :: lxyz
    CHARACTER(LEN=12)                        :: symbol

    CHARACTER(len=*), PARAMETER :: routineN = 'cgf_symbol', &
      routineP = moduleN//':'//routineN
    CHARACTER(LEN=1), DIMENSION(3), &
      PARAMETER                              :: xyz = (/"x","y","z"/)

    INTEGER                                  :: i, ipos, l

    symbol = ""

    IF ((n > 0).AND.(n < 100)) THEN
      WRITE (symbol(1:2),"(I2)") n
    ELSE
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "Invalid principal quantum number specified")
    END IF

    l = SUM(lxyz(1:3))

    IF ((l >= 0).AND.(l <= 11)) THEN
      symbol(3:3) = l_sym(l)
    ELSE
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "Invalid angular momentum quantum number specified")
    END IF

    ipos = 4

    DO i=1,3
      IF (lxyz(i) > 0) THEN
        symbol(ipos:ipos) = xyz(i)
        ipos = ipos + 1
        IF (lxyz(i) > 1) THEN
          IF (lxyz(i) < 10) THEN
            WRITE (symbol(ipos:ipos),"(I1)") lxyz(i)
            ipos = ipos + 1
          ELSE IF (lxyz(i) < 100) THEN
            WRITE (symbol(ipos:ipos+1),"(I2)") lxyz(i)
            ipos = ipos + 2
          ELSE
            CALL stop_program(routineN,moduleN,__LINE__,&
                              "Invalid magnetic quantum number specified")
          END IF
        END IF
      END IF
    END DO

  END FUNCTION cgf_symbol

! *****************************************************************************
!> \brief   Build a spherical orbital symbol (orbital labels for printing).
!> \param n ...
!> \param l ...
!> \param m ...
!> \retval symbol ...
!> \date    11.03.99
!> \par Variables
!>       - l: Angular momentum quantum number l of the orbital.
!>       - m: Magnetic quantum number m of the orbital.
!>       - n: Principle quantum number n of the orbital.
!> \par History
!>  - Ignore n value for n = 0 (16.02.2009,MK)
!> \author  Matthias Krack
!> \version 1.0
! *****************************************************************************
  FUNCTION sgf_symbol(n,l,m) RESULT(symbol)
    INTEGER, INTENT(IN)                      :: n, l, m
    CHARACTER(LEN=6)                         :: symbol

    CHARACTER(len=*), PARAMETER :: routineN = 'sgf_symbol', &
      routineP = moduleN//':'//routineN
    CHARACTER(LEN=1), DIMENSION(-1:1), &
      PARAMETER                              :: yzx = (/"y","z","x"/)

    INTEGER                                  :: i

    symbol = ""

    IF (n == 0) THEN
      i = 1
    ELSE IF ((n > 0).AND.(n < 100)) THEN
      WRITE (symbol(1:2),"(I2)") n
      i = 3
    ELSE
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "Invalid principal quantum number specified")
    END IF

    IF ((l >= 0).AND.(l <= 11)) THEN
      symbol(i:i) = l_sym(l)
      i = i + 1
    ELSE
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "Invalid angular momentum quantum number specified")
    END IF

    IF (ABS(m) <= l) THEN
      IF (l == 1) THEN
        symbol(i:i) = yzx(m)
      ELSE IF (l > 1) THEN
        IF (m == 0) THEN
          WRITE (symbol(i:i),"(I1)") m
        ELSE IF (ABS(m) < 10) THEN
          WRITE (symbol(i:i+1),"(SP,I2)") m
        ELSE IF (ABS(m) < 100) THEN
          WRITE (symbol(i:i+2),"(SP,I3)") m
        END IF
      END IF
    ELSE
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "Invalid magnetic quantum number specified")
    END IF

  END FUNCTION sgf_symbol

END MODULE orbital_symbols
