!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/orbital_symbols [1.0] *
!!
!!   NAME
!!     orbital_symbols
!!
!!   FUNCTION
!!     orbital_symbols
!!
!!   AUTHOR
!!     Matthias Krack (08.06.2000)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE orbital_symbols
  USE termination, ONLY: stop_program


! *****************************************************************************

! Index:

! FUNCTION cgf_symbol(n,lxyz) RESULT(symbol)
! FUNCTION sgf_symbol(n,l,m) RESULT(symbol)

! *****************************************************************************

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=1), PARAMETER, DIMENSION(0:11) :: l_sym = (/"s","p","d",&
                                                            "f","g","h",&
                                                            "i","j","k",&
                                                            "l","m","n"/)

  PUBLIC :: cgf_symbol,sgf_symbol

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  FUNCTION cgf_symbol(n,lxyz) RESULT(symbol)

!   Purpose: Build a Cartesian orbital symbol (orbital labels for printing).

!   History: - Creation (07.07.99, Matthias Krack)

!   ***************************************************************************

    INTEGER, INTENT(IN) :: n

    INTEGER, DIMENSION(3), INTENT(IN) :: lxyz

    CHARACTER(LEN=12) :: symbol

!   *** Local parameters ***

    CHARACTER(LEN=1), PARAMETER, DIMENSION(3) :: xyz = (/"x","y","z"/)

!   *** Local variables ***

    INTEGER :: i,ipos,l

!   ---------------------------------------------------------------------------

    symbol = ""

    IF ((n > 0).AND.(n < 100)) THEN
      WRITE (symbol(1:2),"(I2)") n
    ELSE
      CALL stop_program("FUNCTION cgf_symbol (MODULE orbital_symbols)",&
                        "Invalid principal quantum number specified")
    END IF

    l = SUM(lxyz(1:3))

    IF ((l >= 0).AND.(l <= 11)) THEN
      symbol(3:3) = l_sym(l)
    ELSE
      CALL stop_program("FUNCTION cgf_symbol (MODULE orbital_symbols)",&
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
            CALL stop_program("FUNCTION cgf_symbol (MODULE orbital_symbols)",&
                              "Invalid magnetic quantum number specified")
          END IF
        END IF
      END IF
    END DO

  END FUNCTION cgf_symbol

! *****************************************************************************

  FUNCTION sgf_symbol(n,l,m) RESULT(symbol)

!   Purpose: Build a spherical orbital symbol (orbital labels for printing).

!   History: - Creation (11.03.99, Matthias Krack)

!   ***************************************************************************

!   l: Angular momentum quantum number l of the orbital.
!   m: Magnetic quantum number m of the orbital.
!   n: Principle quantum number n of the orbital.

!   ***************************************************************************

    INTEGER, INTENT(IN) :: l,m,n

    CHARACTER(LEN=6) :: symbol

!   *** Local parameters ***

    CHARACTER(LEN=1), PARAMETER, DIMENSION(-1:1) :: yzx = (/"y","z","x"/)

!   ---------------------------------------------------------------------------

    symbol = ""

    IF ((n > 0).AND.(n < 100)) THEN
      WRITE (symbol(1:2),"(I2)") n
    ELSE
      CALL stop_program("FUNCTION sgf_symbol (MODULE orbital_symbols)",&
                        "Invalid principal quantum number specified")
    END IF

    IF ((l >= 0).AND.(l <= 11)) THEN
      symbol(3:3) = l_sym(l)
    ELSE
      CALL stop_program("FUNCTION sgf_symbol (MODULE orbital_symbols)",&
                        "Invalid angular momentum quantum number specified")
    END IF

    IF (ABS(m) <= l) THEN
      IF (l == 1) THEN
        symbol(4:4) = yzx(m)
      ELSE IF (l > 1) THEN
        IF (m == 0) THEN
          WRITE (symbol(4:4),"(I1)") m
        ELSE IF (ABS(m) < 10) THEN
          WRITE (symbol(4:5),"(SP,I2)") m
        ELSE IF (ABS(m) < 100) THEN
          WRITE (symbol(4:6),"(SP,I3)") m
        END IF
      END IF
    ELSE
      CALL stop_program("FUNCTION sgf_symbol (MODULE orbital_symbols)",&
                        "Invalid magnetic quantum number specified")
    END IF

  END FUNCTION sgf_symbol

! *****************************************************************************

END MODULE orbital_symbols
