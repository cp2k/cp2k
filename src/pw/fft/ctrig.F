!-----------------------------------------------------------------------------!
!  Copyright (C) Stefan Goedecker, Lausanne, Switzerland, August 1, 1991
!  Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!  Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!  This file is distributed under the terms of the
!  GNU General Public License version 2 (or later),
!  see http://www.gnu.org/copyleft/gpl.txt .
!-----------------------------------------------------------------------------!
!  S. Goedecker: Rotating a three-dimensional array in optimal
!  positions for vector processing: Case study for a three-dimensional Fast
!  Fourier Transform, Comp. Phys. Commun. 76, 294 (1993)
! *****************************************************************************
SUBROUTINE ctrig ( n, trig, after, before, now, isign, ic )
  USE fft_kinds, ONLY                      : dp
  INTEGER, INTENT ( IN ) :: n
  INTEGER, INTENT ( IN ) :: isign
  INTEGER, INTENT ( OUT ) :: ic
  INTEGER, DIMENSION ( 7 ), INTENT ( OUT ) :: after, before, now
  REAL ( dp ) , DIMENSION ( 2, 1024 ), INTENT ( OUT ) :: trig

  INTEGER :: i, j, itt
  REAL ( dp ) :: twopi, angle
  INTEGER, PARAMETER :: nt = 82
  INTEGER, DIMENSION ( 7, nt ), PARAMETER :: idata = RESHAPE ((/ &
            3,   3, 1, 1, 1, 1, 1,       4,   4, 1, 1, 1, 1, 1,      5,   5, 1, 1, 1, 1, 1,       6,   6, 1, 1, 1, 1, 1, &
            8,   8, 1, 1, 1, 1, 1,       9,   3, 3, 1, 1, 1, 1,     12,   4, 3, 1, 1, 1, 1,      15,   5, 3, 1, 1, 1, 1, &
           16,   4, 4, 1, 1, 1, 1,      18,   6, 3, 1, 1, 1, 1,     20,   5, 4, 1, 1, 1, 1,      24,   8, 3, 1, 1, 1, 1, &
           25,   5, 5, 1, 1, 1, 1,      27,   3, 3, 3, 1, 1, 1,     30,   6, 5, 1, 1, 1, 1,      32,   8, 4, 1, 1, 1, 1, &
           36,   4, 3, 3, 1, 1, 1,      40,   8, 5, 1, 1, 1, 1,     45,   5, 3, 3, 1, 1, 1,      48,   4, 4, 3, 1, 1, 1, &
           54,   6, 3, 3, 1, 1, 1,      60,   5, 4, 3, 1, 1, 1,     64,   4, 4, 4, 1, 1, 1,      72,   8, 3, 3, 1, 1, 1, &
           75,   5, 5, 3, 1, 1, 1,      80,   5, 4, 4, 1, 1, 1,     81,   3, 3, 3, 3, 1, 1,      90,   6, 5, 3, 1, 1, 1, &
           96,   8, 4, 3, 1, 1, 1,     100,   5, 5, 4, 1, 1, 1,    108,   4, 3, 3, 3, 1, 1,     120,   8, 5, 3, 1, 1, 1, &
          125,   5, 5, 5, 1, 1, 1,     128,   8, 4, 4, 1, 1, 1,    135,   5, 3, 3, 3, 1, 1,     144,   4, 4, 3, 3, 1, 1, &
          150,   6, 5, 5, 1, 1, 1,     160,   8, 5, 4, 1, 1, 1,    162,   6, 3, 3, 3, 1, 1,     180,   5, 4, 3, 3, 1, 1, &
          192,   4, 4, 4, 3, 1, 1,     200,   8, 5, 5, 1, 1, 1,    216,   8, 3, 3, 3, 1, 1,     225,   5, 5, 3, 3, 1, 1, &
          240,   5, 4, 4, 3, 1, 1,     243,   3, 3, 3, 3, 3, 1,    256,   4, 4, 4, 4, 1, 1,     270,   6, 5, 3, 3, 1, 1, &
          288,   8, 4, 3, 3, 1, 1,     300,   5, 5, 4, 3, 1, 1,    320,   5, 4, 4, 4, 1, 1,     324,   4, 3, 3, 3, 3, 1, &
          360,   8, 5, 3, 3, 1, 1,     375,   5, 5, 5, 3, 1, 1,    384,   8, 4, 4, 3, 1, 1,     400,   5, 5, 4, 4, 1, 1, &
          405,   5, 3, 3, 3, 3, 1,     432,   4, 4, 3, 3, 3, 1,    450,   6, 5, 5, 3, 1, 1,     480,   8, 5, 4, 3, 1, 1, &
          486,   6, 3, 3, 3, 3, 1,     500,   5, 5, 5, 4, 1, 1,    512,   8, 4, 4, 4, 1, 1,     540,   5, 4, 3, 3, 3, 1, &
          576,   4, 4, 4, 3, 3, 1,     600,   8, 5, 5, 3, 1, 1,    625,   5, 5, 5, 5, 1, 1,     640,   8, 5, 4, 4, 1, 1, &
          648,   8, 3, 3, 3, 3, 1,     675,   5, 5, 3, 3, 3, 1,    720,   5, 4, 4, 3, 3, 1,     729,   3, 3, 3, 3, 3, 3, &
          750,   6, 5, 5, 5, 1, 1,     768,   4, 4, 4, 4, 3, 1,    800,   8, 5, 5, 4, 1, 1,     810,   6, 5, 3, 3, 3, 1, &
          864,   8, 4, 3, 3, 3, 1,     900,   5, 5, 4, 3, 3, 1,    960,   5, 4, 4, 4, 3, 1,     972,   4, 3, 3, 3, 3, 3, &
         1000,   8, 5, 5, 5, 1, 1,    1024,   4, 4, 4, 4, 4, 1  /),(/7,nt/))


  mloop: DO i = 1, nt
    IF ( n == idata ( 1, i ) ) THEN
      ic=0
      DO j = 1, 6
        itt = idata ( 1 + j, i )
        IF ( itt > 1 ) THEN
          ic = ic + 1
          now ( j ) = idata ( 1 + j, i )
        ELSE
          EXIT mloop
        END IF
      END DO
      EXIT mloop
    END IF
    IF ( i == nt ) THEN
      WRITE ( *, '(A,i5,A)' ) " Value of ",n, &
            " not allowed for fft, allowed values are:"
      WRITE ( *, '(15i5)' ) ( idata ( 1, j ), j = 1, nt )
      STOP 'ctrig'
    END IF
  END DO mloop

  after ( 1 ) = 1
  before ( ic ) = 1
  DO i = 2, ic
    after ( i ) = after ( i - 1 ) * now ( i - 1 )
    before ( ic - i + 1 ) = before ( ic - i + 2 ) * now ( ic - i + 2 )
  END DO

  twopi = 8._dp * ATAN ( 1._dp )
  angle = isign * twopi / REAL ( n, dp )
  trig ( 1, 1 ) = 1._dp
  trig ( 2, 1 ) = 0._dp
  DO i = 1, n - 1
    trig ( 1, i + 1 ) = COS ( REAL ( i, dp ) * angle )
    trig ( 2, i + 1 ) = SIN ( REAL ( i, dp ) * angle )
  END DO

  END SUBROUTINE ctrig
