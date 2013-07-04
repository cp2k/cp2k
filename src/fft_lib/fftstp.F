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
SUBROUTINE fftstp ( mm, nfft, m, nn, n, zin, zout,  trig, now, after, before, isign )

  USE fft_kinds, ONLY                                   : dp
  INTEGER, INTENT ( IN ) :: mm, nfft, m, nn, n, now, after, before, isign
  REAL ( dp ), DIMENSION ( 2, 1024 ), INTENT ( IN ) :: trig
  REAL ( dp ), DIMENSION ( 2, mm, m ), INTENT ( IN ) :: zin
  REAL ( dp ), DIMENSION ( 2, nn, n ), INTENT ( INOUT ) :: zout

  INTEGER :: atn, atb, ia, ib, nin1, nin2, nin3, nin4, nin5, nin6, nin7, nin8
  INTEGER :: nout1, nout2, nout3, nout4, nout5, nout6, nout7, nout8, &
             j, ias, itt, itrig
  REAL ( dp ) :: s, s1, s2, s3, s4, s5, s6, s7, s8, &
                  r, r1, r2, r3, r4, r5, r6, r7, r8, cr2, cr3, cr4, cr5, &
                  ci2, ci3, ci4, ci5, ur1, ur2, ur3, ui1, ui2, ui3, &
                  vr1, vr2, vr3, vi1, vi2, vi3, cm, cp, dm, dbl, &
                  am, ap, bm, bp, bbs, s25, s34, r34, r25, sin2, sin4
  REAL ( dp ), PARAMETER :: rt2i = 0.7071067811865475_dp  ! sqrt(0.5)
  REAL ( dp ), PARAMETER :: bb = 0.8660254037844387_dp  ! sqrt(3)/2
  REAL ( dp ), PARAMETER :: cos2 = 0.3090169943749474_dp ! cos(2*pi/5)
  REAL ( dp ), PARAMETER :: cos4 = - 0.8090169943749474_dp !  cos(4*pi/5)
  REAL ( dp ), PARAMETER :: sin2p = 0.9510565162951536_dp ! sin(2*pi/5)
  REAL ( dp ), PARAMETER :: sin4p = 0.5877852522924731_dp ! sin(4*pi/5)

!-----------------------------------------------------------------------------!

  atn = after * now
  atb = after * before

  IF ( now == 4 ) THEN
    IF ( isign == 1 ) THEN
      ia = 1
      nin1 = ia - after
      nout1 = ia - atn
      DO ib = 1, before
        nin1 = nin1 + after
        nin2 = nin1 + atb
        nin3 = nin2 + atb
        nin4 = nin3 + atb
        nout1 = nout1 + atn
        nout2 = nout1 + after
        nout3 = nout2 + after
        nout4 = nout3 + after
        DO j = 1, nfft
          r1 = zin ( 1, j, nin1 )
          s1 = zin ( 2, j, nin1 )
          r2 = zin ( 1, j, nin2 )
          s2 = zin ( 2, j, nin2 )
          r3 = zin ( 1, j, nin3 )
          s3 = zin ( 2, j, nin3 )
          r4 = zin ( 1, j, nin4 )
          s4 = zin ( 2, j, nin4 )
          r = r1 + r3
          s = r2 + r4
          zout ( 1, j, nout1 ) = r + s
          zout ( 1, j, nout3 ) = r - s
          r = r1 - r3
          s = s2 - s4
          zout ( 1, j, nout2 ) = r - s
          zout ( 1, j, nout4 ) = r + s
          r = s1 + s3
          s = s2 + s4
          zout ( 2, j, nout1 ) = r + s
          zout ( 2, j, nout3 ) = r - s
          r = s1 - s3
          s = r2 - r4
          zout ( 2, j, nout2 ) = r + s
          zout ( 2, j, nout4 ) = r - s
        END DO
      END DO
      DO ia = 2, after
        ias = ia - 1
        IF ( 2*ias == after ) THEN
          nin1 = ia - after
          nout1 = ia - atn
          DO ib = 1, before
            nin1 = nin1 + after
            nin2 = nin1 + atb
            nin3 = nin2 + atb
            nin4 = nin3 + atb
            nout1 = nout1 + atn
            nout2 = nout1 + after
            nout3 = nout2 + after
            nout4 = nout3 + after
            DO j = 1, nfft
              r1 = zin ( 1, j, nin1 )
              s1 = zin ( 2, j, nin1 )
              r = zin ( 1, j, nin2 )
              s = zin ( 2, j, nin2 )
              r2 = ( r - s ) * rt2i
              s2 = ( r + s ) * rt2i
              r3 = - zin ( 2, j, nin3 )
              s3 = zin ( 1, j, nin3 )
              r = zin ( 1, j, nin4 )
              s = zin ( 2, j, nin4 )
              r4 = - ( r + s ) * rt2i
              s4 = ( r -  s ) * rt2i
              r = r1 + r3
              s = r2 + r4
              zout ( 1, j, nout1 ) = r + s
              zout ( 1, j, nout3 ) = r - s
              r = r1 - r3
              s = s2 - s4
              zout ( 1, j, nout2 ) = r - s
              zout ( 1, j, nout4 ) = r + s
              r = s1 + s3
              s = s2 + s4
              zout ( 2, j, nout1 ) = r + s
              zout ( 2, j, nout3 ) = r - s
              r = s1 - s3
              s = r2 - r4
              zout ( 2, j, nout2 ) = r + s
              zout ( 2, j, nout4 ) = r - s
            END DO
          END DO
        ELSE
          itt = ias * before
          itrig = itt + 1
          cr2 = trig ( 1, itrig )
          ci2 = trig ( 2, itrig )
          itrig = itrig + itt
          cr3 = trig ( 1, itrig )
          ci3 = trig ( 2, itrig )
          itrig = itrig + itt
          cr4 = trig ( 1, itrig )
          ci4 = trig ( 2, itrig )
          nin1 = ia - after
          nout1 = ia - atn
          DO ib = 1, before
            nin1 = nin1 + after
            nin2 = nin1 + atb
            nin3 = nin2 + atb
            nin4 = nin3 + atb
            nout1 = nout1 + atn
            nout2 = nout1 + after
            nout3 = nout2 + after
            nout4 = nout3 + after
            DO j = 1, nfft
              r1 = zin ( 1, j, nin1 )
              s1 = zin ( 2, j, nin1 )
              r = zin ( 1, j, nin2 )
              s = zin ( 2, j, nin2 )
              r2 = r * cr2 - s * ci2
              s2 = r * ci2 + s * cr2
              r = zin ( 1, j, nin3 )
              s = zin ( 2, j, nin3 )
              r3 = r * cr3 - s * ci3
              s3 = r * ci3 + s * cr3
              r = zin ( 1, j, nin4 )
              s = zin ( 2, j, nin4 )
              r4 = r * cr4 - s * ci4
              s4 = r * ci4 + s * cr4
              r = r1 + r3
              s = r2 + r4
              zout ( 1, j, nout1 ) = r + s
              zout ( 1, j, nout3 ) = r - s
              r = r1 - r3
              s = s2 - s4
              zout ( 1, j, nout2 ) = r - s
              zout ( 1, j, nout4 ) = r + s
              r = s1 + s3
              s = s2 + s4
              zout ( 2, j, nout1 ) = r + s
              zout ( 2, j, nout3 ) = r - s
              r = s1 - s3
              s = r2 - r4
              zout ( 2, j, nout2 ) = r + s
              zout ( 2, j, nout4 ) = r - s
            END DO
          END DO
        END IF
      END DO
    ELSE
      ia = 1
      nin1 = ia - after
      nout1 = ia - atn
      DO ib = 1, before
        nin1 = nin1 + after
        nin2 = nin1 + atb
        nin3 = nin2 + atb
        nin4 = nin3 + atb
        nout1 = nout1 + atn
        nout2 = nout1 + after
        nout3 = nout2 + after
        nout4 = nout3 + after
        DO j = 1, nfft
          r1 = zin ( 1, j, nin1 )
          s1 = zin ( 2, j, nin1 )
          r2 = zin ( 1, j, nin2 )
          s2 = zin ( 2, j, nin2 )
          r3 = zin ( 1, j, nin3 )
          s3 = zin ( 2, j, nin3 )
          r4 = zin ( 1, j, nin4 )
          s4 = zin ( 2, j, nin4 )
          r = r1 + r3
          s = r2 + r4
          zout ( 1, j, nout1 ) = r + s
          zout ( 1, j, nout3 ) = r - s
          r = r1 - r3
          s = s2 - s4
          zout ( 1, j, nout2 ) = r + s
          zout ( 1, j, nout4 ) = r - s
          r = s1 + s3
          s = s2 + s4
          zout ( 2, j, nout1 ) = r + s
          zout ( 2, j, nout3 ) = r - s
          r = s1 - s3
          s = r2 - r4
          zout ( 2, j, nout2 ) = r - s
          zout ( 2, j, nout4 ) = r + s
        END DO
      END DO
      DO ia = 2, after
        ias = ia - 1
        IF ( 2 * ias == after ) THEN
          nin1 = ia - after
          nout1 = ia - atn
          DO ib = 1, before
            nin1 = nin1 + after
            nin2 = nin1 + atb
            nin3 = nin2 + atb
            nin4 = nin3 + atb
            nout1 = nout1 + atn
            nout2 = nout1 + after
            nout3 = nout2 + after
            nout4 = nout3 + after
            DO j = 1, nfft
              r1 = zin ( 1, j, nin1 )
              s1 = zin ( 2, j, nin1 )
              r = zin ( 1, j, nin2 )
              s = zin ( 2, j, nin2 )
              r2 = ( r + s ) * rt2i
              s2 = ( s - r ) * rt2i
              r3 = zin ( 2, j, nin3 )
              s3 = - zin ( 1, j, nin3 )
              r = zin ( 1, j, nin4 )
              s = zin ( 2, j, nin4 )
              r4 = ( s - r ) * rt2i
              s4 = - ( r + s ) * rt2i
              r = r1 + r3
              s = r2 + r4
              zout ( 1, j, nout1 ) = r + s
              zout ( 1, j, nout3 ) = r - s
              r = r1 - r3
              s = s2 - s4
              zout ( 1, j, nout2 ) = r + s
              zout ( 1, j, nout4 ) = r - s
              r =s1 + s3
              s =s2 + s4
              zout ( 2, j, nout1 ) = r + s
              zout ( 2, j, nout3 ) = r - s
              r = s1 - s3
              s = r2 - r4
              zout ( 2, j, nout2 ) = r - s
              zout ( 2, j, nout4 ) = r + s
            END DO
          END DO
        ELSE
          itt = ias * before
          itrig = itt + 1
          cr2 = trig ( 1, itrig )
          ci2 = trig ( 2, itrig )
          itrig = itrig + itt
          cr3 = trig ( 1, itrig )
          ci3 = trig ( 2, itrig )
          itrig = itrig + itt
          cr4 = trig ( 1, itrig )
          ci4 = trig ( 2, itrig )
          nin1 = ia - after
          nout1 = ia - atn
          DO ib = 1, before
            nin1 = nin1 + after
            nin2 = nin1 + atb
            nin3 = nin2 + atb
            nin4 = nin3 + atb
            nout1 = nout1 + atn
            nout2 = nout1 + after
            nout3 = nout2 + after
            nout4 = nout3 + after
            DO j = 1, nfft
              r1 = zin ( 1, j, nin1 )
              s1 = zin ( 2, j, nin1 )
              r = zin ( 1, j, nin2 )
              s = zin ( 2, j, nin2 )
              r2 = r * cr2 - s * ci2
              s2 = r * ci2 + s * cr2
              r = zin ( 1, j, nin3 )
              s = zin ( 2, j, nin3 )
              r3 = r * cr3 - s * ci3
              s3 = r * ci3 + s * cr3
              r = zin ( 1, j, nin4 )
              s = zin ( 2, j, nin4 )
              r4 = r * cr4 - s * ci4
              s4 = r * ci4 + s * cr4
              r = r1 + r3
              s = r2 + r4
              zout ( 1, j, nout1 ) = r + s
              zout ( 1, j, nout3 ) = r - s
              r = r1 - r3
              s = s2 - s4
              zout ( 1, j, nout2 ) = r + s
              zout ( 1, j, nout4 ) = r - s
              r = s1 + s3
              s = s2 + s4
              zout ( 2, j, nout1 ) = r + s
              zout ( 2, j, nout3 ) = r - s
              r = s1 - s3
              s = r2 - r4
              zout ( 2, j, nout2 ) = r - s
              zout ( 2, j, nout4 ) = r + s
            END DO
          END DO
        END IF
      END DO
    END IF
  ELSE IF ( now == 8 ) THEN
    IF ( isign == -1 ) THEN
      ia = 1
      nin1 = ia - after
      nout1 = ia - atn
      DO ib = 1, before
        nin1 = nin1 + after
        nin2 = nin1 + atb
        nin3 = nin2 + atb
        nin4 = nin3 + atb
        nin5 = nin4 + atb
        nin6 = nin5 + atb
        nin7 = nin6 + atb
        nin8 = nin7 + atb
        nout1 = nout1 + atn
        nout2 = nout1 + after
        nout3 = nout2 + after
        nout4 = nout3 + after
        nout5 = nout4 + after
        nout6 = nout5 + after
        nout7 = nout6 + after
        nout8 = nout7 + after
        DO j = 1, nfft
          r1 = zin ( 1, j, nin1 )
          s1 = zin ( 2, j, nin1 )
          r2 = zin ( 1, j, nin2 )
          s2 = zin ( 2, j, nin2 )
          r3 = zin ( 1, j, nin3 )
          s3 = zin ( 2, j, nin3 )
          r4 = zin ( 1, j, nin4 )
          s4 = zin ( 2, j, nin4 )
          r5 = zin ( 1, j, nin5 )
          s5 = zin ( 2, j, nin5 )
          r6 = zin ( 1, j, nin6 )
          s6 = zin ( 2, j, nin6 )
          r7 = zin ( 1, j, nin7 )
          s7 = zin ( 2, j, nin7 )
          r8 = zin ( 1, j, nin8 )
          s8 = zin ( 2, j, nin8 )
          r = r1 + r5
          s = r3 + r7
          ap = r + s
          am = r - s
          r = r2 + r6
          s = r4 + r8
          bp = r + s
          bm = r - s
          r = s1 + s5
          s = s3 + s7
          cp = r + s
          cm = r - s
          r = s2 + s6
          s = s4 + s8
          dbl = r + s
          dm = r - s
          zout ( 1, j, nout1 ) = ap + bp
          zout ( 2, j, nout1 ) = cp + dbl
          zout ( 1, j, nout5 ) = ap - bp
          zout ( 2, j, nout5 ) = cp - dbl
          zout ( 1, j, nout3 ) = am + dm
          zout ( 2, j, nout3 ) = cm - bm
          zout ( 1, j, nout7 ) = am - dm
          zout ( 2, j, nout7 ) = cm + bm
          r = r1 - r5
          s = s3 - s7
          ap = r + s
          am = r - s
          r = s1 - s5
          s = r3 - r7
          bp = r + s
          bm = r - s
          r = s4 - s8
          s = r2 - r6
          cp = r + s
          cm = r - s
          r = s2 - s6
          s = r4 - r8
          dbl = r + s
          dm = r - s
          r = ( cp + dm ) * rt2i
          s = (-cp + dm ) * rt2i
          cp = ( cm + dbl ) * rt2i
          dbl = ( cm - dbl ) * rt2i
          zout ( 1, j, nout2 ) = ap + r
          zout ( 2, j, nout2 ) = bm + s
          zout ( 1, j, nout6 ) = ap - r
          zout ( 2, j, nout6 ) = bm - s
          zout ( 1, j, nout4 ) = am + cp
          zout ( 2, j, nout4 ) = bp + dbl
          zout ( 1, j, nout8 ) = am - cp
          zout ( 2, j, nout8 ) = bp - dbl
        END DO
      END DO
    ELSE
      ia = 1
      nin1 = ia - after
      nout1 = ia - atn
      DO ib = 1, before
        nin1 = nin1 + after
        nin2 = nin1 + atb
        nin3 = nin2 + atb
        nin4 = nin3 + atb
        nin5 = nin4 + atb
        nin6 = nin5 + atb
        nin7 = nin6 + atb
        nin8 = nin7 + atb
        nout1 = nout1 + atn
        nout2 = nout1 + after
        nout3 = nout2 + after
        nout4 = nout3 + after
        nout5 = nout4 + after
        nout6 = nout5 + after
        nout7 = nout6 + after
        nout8 = nout7 + after
        DO j = 1, nfft
          r1 = zin ( 1, j, nin1 )
          s1 = zin ( 2, j, nin1 )
          r2 = zin ( 1, j, nin2 )
          s2 = zin ( 2, j, nin2 )
          r3 = zin ( 1, j, nin3 )
          s3 = zin ( 2, j, nin3 )
          r4 = zin ( 1, j, nin4 )
          s4 = zin ( 2, j, nin4 )
          r5 = zin ( 1, j, nin5 )
          s5 = zin ( 2, j, nin5 )
          r6 = zin ( 1, j, nin6 )
          s6 = zin ( 2, j, nin6 )
          r7 = zin ( 1, j, nin7 )
          s7 = zin ( 2, j, nin7 )
          r8 = zin ( 1, j, nin8 )
          s8 = zin ( 2, j, nin8 )
          r = r1 + r5
          s = r3 + r7
          ap = r + s
          am = r - s
          r = r2 + r6
          s = r4 + r8
          bp = r + s
          bm = r - s
          r = s1 + s5
          s = s3 + s7
          cp = r + s
          cm = r - s
          r = s2 + s6
          s = s4 + s8
          dbl = r + s
          dm = r - s
          zout ( 1, j, nout1 ) = ap + bp
          zout ( 2, j, nout1 ) = cp + dbl
          zout ( 1, j, nout5 ) = ap - bp
          zout ( 2, j, nout5 ) = cp - dbl
          zout ( 1, j, nout3 ) = am - dm
          zout ( 2, j, nout3 ) = cm + bm
          zout ( 1, j, nout7 ) = am + dm
          zout ( 2, j, nout7 ) = cm - bm
          r = r1 - r5
          s = -s3 + s7
          ap = r + s
          am = r - s
          r = s1 - s5
          s = r7 - r3
          bp = r + s
          bm = r - s
          r = -s4 + s8
          s = r2 - r6
          cp = r + s
          cm = r - s
          r = -s2 + s6
          s = r4 - r8
          dbl = r + s
          dm = r - s
          r = ( cp + dm ) * rt2i
          s = ( cp - dm ) * rt2i
          cp = ( cm + dbl ) * rt2i
          dbl = (-cm + dbl ) * rt2i
          zout ( 1, j, nout2 ) = ap + r
          zout ( 2, j, nout2 ) = bm + s
          zout ( 1, j, nout6 ) = ap - r
          zout ( 2, j, nout6 ) = bm - s
          zout ( 1, j, nout4 ) = am + cp
          zout ( 2, j, nout4 ) = bp + dbl
          zout ( 1, j, nout8 ) = am - cp
          zout ( 2, j, nout8 ) = bp - dbl
        END DO
      END DO
    END IF
  ELSE IF ( now == 3 ) THEN
    bbs = isign * bb
    ia = 1
    nin1 = ia - after
    nout1 = ia - atn
    DO ib = 1, before
      nin1 = nin1 + after
      nin2 = nin1 + atb
      nin3 = nin2 + atb
      nout1 = nout1 + atn
      nout2 = nout1 + after
      nout3 = nout2 + after
      DO j = 1, nfft
        r1 = zin ( 1, j, nin1 )
        s1 = zin ( 2, j, nin1 )
        r2 = zin ( 1, j, nin2 )
        s2 = zin ( 2, j, nin2 )
        r3 = zin ( 1, j, nin3 )
        s3 = zin ( 2, j, nin3 )
        r = r2 + r3
        s = s2 + s3
        zout ( 1, j, nout1 ) = r + r1
        zout ( 2, j, nout1 ) = s + s1
        r1 = r1 - 0.5_dp * r
        s1 = s1 - 0.5_dp * s
        r2 = bbs * ( r2 - r3 )
        s2 = bbs * ( s2 - s3 )
        zout ( 1, j, nout2 ) = r1 - s2
        zout ( 2, j, nout2 ) = s1 + r2
        zout ( 1, j, nout3 ) = r1 + s2
        zout ( 2, j, nout3 ) = s1 - r2
      END DO
    END DO
    DO ia = 2, after
      ias = ia - 1
      IF ( 4*ias == 3*after ) THEN
        IF ( isign == 1 ) THEN
          nin1 = ia - after
          nout1 = ia - atn
          DO ib = 1, before
            nin1 = nin1 + after
            nin2=nin1+atb
            nin3=nin2+atb
            nout1=nout1+atn
            nout2=nout1+after
            nout3=nout2+after
            DO j = 1, nfft
              r1 = zin ( 1, j, nin1 )
              s1 = zin ( 2, j, nin1 )
              r2 = -zin ( 2, j, nin2 )
              s2 = zin ( 1, j, nin2 )
              r3 = -zin ( 1, j, nin3 )
              s3 = -zin ( 2, j, nin3 )
              r = r2 + r3
              s = s2 + s3
              zout ( 1, j, nout1 ) = r + r1
              zout ( 2, j, nout1 ) = s + s1
              r1 = r1 - 0.5_dp * r
              s1 = s1 - 0.5_dp * s
              r2 = bbs*(r2-r3)
              s2 = bbs*(s2-s3)
              zout ( 1, j, nout2 ) = r1 - s2
              zout ( 2, j, nout2 ) = s1 + r2
              zout ( 1, j, nout3 ) = r1 + s2
              zout ( 2, j, nout3 ) = s1 - r2
            END DO
          END DO
        ELSE
          nin1 = ia - after
          nout1 = ia - atn
          DO ib = 1, before
            nin1 = nin1 + after
            nin2 = nin1 + atb
            nin3 = nin2 + atb
            nout1 = nout1 + atn
            nout2 = nout1 + after
            nout3 = nout2 + after
            DO j = 1, nfft
              r1 = zin ( 1, j, nin1 )
              s1 = zin ( 2, j, nin1 )
              r2 = zin ( 2, j, nin2 )
              s2 = -zin ( 1, j, nin2 )
              r3 = -zin ( 1, j, nin3 )
              s3 = -zin ( 2, j, nin3 )
              r = r2 + r3
              s = s2 + s3
              zout ( 1, j, nout1 ) = r + r1
              zout ( 2, j, nout1 ) = s + s1
              r1 = r1 - 0.5_dp * r
              s1 = s1 - 0.5_dp * s
              r2 = bbs * ( r2 - r3 )
              s2 = bbs * ( s2 - s3 )
              zout ( 1, j, nout2 ) = r1 - s2
              zout ( 2, j, nout2 ) = s1 + r2
              zout ( 1, j, nout3 ) = r1 + s2
              zout ( 2, j, nout3 ) = s1 - r2
            END DO
          END DO
        END IF
      ELSE IF ( 8 * ias == 3 * after ) THEN
        IF ( isign == 1 ) THEN
          nin1 = ia - after
          nout1 = ia - atn
          DO ib = 1, before
            nin1 = nin1 + after
            nin2 = nin1 + atb
            nin3 = nin2 + atb
            nout1 = nout1 + atn
            nout2 = nout1 + after
            nout3 = nout2 + after
            DO j = 1, nfft
              r1 = zin ( 1, j, nin1 )
              s1 = zin ( 2, j, nin1 )
              r = zin ( 1, j, nin2 )
              s = zin ( 2, j, nin2 )
              r2 = ( r - s ) * rt2i
              s2 = ( r + s ) * rt2i
              r3 = -zin ( 2, j, nin3 )
              s3 = zin ( 1, j, nin3 )
              r = r2 + r3
              s = s2 + s3
              zout ( 1, j, nout1 ) = r + r1
              zout ( 2, j, nout1 ) = s + s1
              r1 = r1 - 0.5_dp * r
              s1 = s1 - 0.5_dp * s
              r2 = bbs * ( r2 - r3 )
              s2 = bbs * ( s2 - s3 )
              zout ( 1, j, nout2 ) = r1 - s2
              zout ( 2, j, nout2 ) = s1 + r2
              zout ( 1, j, nout3 ) = r1 + s2
              zout ( 2, j, nout3 ) = s1 - r2
            END DO
          END DO
        ELSE
          nin1 = ia - after
          nout1 = ia - atn
          DO ib = 1, before
            nin1 = nin1 + after
            nin2 = nin1 + atb
            nin3 = nin2 + atb
            nout1 = nout1 + atn
            nout2 = nout1 + after
            nout3 = nout2 + after
            DO j = 1, nfft
              r1 = zin ( 1, j, nin1 )
              s1 = zin ( 2, j, nin1 )
              r = zin ( 1, j, nin2 )
              s = zin ( 2, j, nin2 )
              r2 = ( r + s ) * rt2i
              s2 = ( -r + s ) * rt2i
              r3 = zin ( 2, j, nin3 )
              s3 = -zin ( 1, j, nin3 )
              r = r2 + r3
              s = s2 + s3
              zout ( 1, j, nout1 ) = r + r1
              zout ( 2, j, nout1 ) = s + s1
              r1 = r1 - 0.5_dp * r
              s1 = s1 - 0.5_dp * s
              r2 = bbs * ( r2 - r3 )
              s2 = bbs * ( s2 - s3 )
              zout ( 1, j, nout2 ) = r1 - s2
              zout ( 2, j, nout2 ) = s1 + r2
              zout ( 1, j, nout3 ) = r1 + s2
              zout ( 2, j, nout3 ) = s1 - r2
            END DO
          END DO
        END IF
      ELSE
        itt = ias * before
        itrig = itt + 1
        cr2 = trig ( 1, itrig )
        ci2 = trig ( 2, itrig )
        itrig = itrig + itt
        cr3 = trig ( 1, itrig )
        ci3 = trig ( 2, itrig )
        nin1 = ia - after
        nout1 = ia - atn
        DO ib = 1, before
          nin1 = nin1 + after
          nin2 = nin1 + atb
          nin3 = nin2 + atb
          nout1 = nout1 + atn
          nout2 = nout1 + after
          nout3 = nout2 + after
          DO j = 1, nfft
            r1 = zin ( 1, j, nin1 )
            s1 = zin ( 2, j, nin1 )
            r = zin ( 1, j, nin2 )
            s = zin ( 2, j, nin2 )
            r2 = r * cr2 - s * ci2
            s2 = r * ci2 + s * cr2
            r = zin ( 1, j, nin3 )
            s = zin ( 2, j, nin3 )
            r3 = r * cr3 - s * ci3
            s3 = r * ci3 + s * cr3
            r = r2 + r3
            s = s2 + s3
            zout ( 1, j, nout1 ) = r + r1
            zout ( 2, j, nout1 ) = s + s1
            r1 = r1 - 0.5_dp * r
            s1 = s1 - 0.5_dp * s
            r2 = bbs * ( r2 - r3 )
            s2 = bbs * ( s2 - s3 )
            zout ( 1, j, nout2 ) = r1 - s2
            zout ( 2, j, nout2 ) = s1 + r2
            zout ( 1, j, nout3 ) = r1 + s2
            zout ( 2, j, nout3 ) = s1 - r2
          END DO
        END DO
      END IF
    END DO
  ELSE IF ( now == 5 ) THEN
    sin2 = isign * sin2p
    sin4 = isign * sin4p
    ia = 1
    nin1 = ia - after
    nout1 = ia - atn
    DO ib = 1, before
      nin1 = nin1 + after
      nin2 = nin1 + atb
      nin3 = nin2 + atb
      nin4 = nin3 + atb
      nin5 = nin4 + atb
      nout1 = nout1 + atn
      nout2 = nout1 + after
      nout3 = nout2 + after
      nout4 = nout3 + after
      nout5 = nout4 + after
      DO j = 1, nfft
        r1 = zin ( 1, j, nin1 )
        s1 = zin ( 2, j, nin1 )
        r2 = zin ( 1, j, nin2 )
        s2 = zin ( 2, j, nin2 )
        r3 = zin ( 1, j, nin3 )
        s3 = zin ( 2, j, nin3 )
        r4 = zin ( 1, j, nin4 )
        s4 = zin ( 2, j, nin4 )
        r5 = zin ( 1, j, nin5 )
        s5 = zin ( 2, j, nin5 )
        r25 = r2 + r5
        r34 = r3 + r4
        s25 = s2 - s5
        s34 = s3 - s4
        zout ( 1, j, nout1 ) = r1 + r25 + r34
        r = cos2 * r25 + cos4 * r34 + r1
        s = sin2 * s25 + sin4 * s34
        zout ( 1, j, nout2 ) = r - s
        zout ( 1, j, nout5 ) = r + s
        r = cos4 * r25 + cos2 * r34 + r1
        s = sin4 * s25 - sin2 * s34
        zout ( 1, j, nout3 ) = r - s
        zout ( 1, j, nout4 ) = r + s
        r25 = r2 - r5
        r34 = r3 - r4
        s25 = s2 + s5
        s34 = s3 + s4
        zout ( 2, j, nout1 ) = s1 + s25 + s34
        r = cos2 * s25 + cos4 * s34 + s1
        s = sin2 * r25 + sin4 * r34
        zout ( 2, j, nout2 ) = r + s
        zout ( 2, j, nout5 ) = r - s
        r = cos4 * s25 + cos2 * s34 + s1
        s = sin4 * r25 - sin2 * r34
        zout ( 2, j, nout3 ) = r + s
        zout ( 2, j, nout4 ) = r - s
      END DO
    END DO
    DO ia = 2, after
      ias = ia - 1
      IF ( 8 * ias == 5 * after ) THEN
        IF ( isign == 1 ) THEN
          nin1 = ia - after
          nout1 = ia - atn
          DO ib = 1, before
            nin1 = nin1 + after
            nin2 = nin1 + atb
            nin3 = nin2 + atb
            nin4 = nin3 + atb
            nin5 = nin4 + atb
            nout1 = nout1 + atn
            nout2 = nout1 + after
            nout3 = nout2 + after
            nout4 = nout3 + after
            nout5 = nout4 + after
            DO j = 1, nfft
              r1 = zin ( 1, j, nin1 )
              s1 = zin ( 2, j, nin1 )
              r = zin ( 1, j, nin2 )
              s = zin ( 2, j, nin2 )
              r2 = ( r - s ) * rt2i
              s2 = ( r + s ) * rt2i
              r3 = -zin ( 2, j, nin3 )
              s3 = zin ( 1, j, nin3 )
              r = zin ( 1, j, nin4 )
              s = zin ( 2, j, nin4 )
              r4 = -( r + s ) * rt2i
              s4 = ( r - s ) * rt2i
              r5 = -zin ( 1, j, nin5 )
              s5 = -zin ( 2, j, nin5 )
              r25 = r2 + r5
              r34 = r3 + r4
              s25 = s2 - s5
              s34 = s3 - s4
              zout ( 1, j, nout1 ) = r1 + r25 + r34
              r = cos2 * r25 + cos4 * r34 + r1
              s = sin2 * s25 + sin4 * s34
              zout ( 1, j, nout2 ) = r - s
              zout ( 1, j, nout5 ) = r + s
              r = cos4 * r25 + cos2 * r34 + r1
              s = sin4 * s25 - sin2 * s34
              zout ( 1, j, nout3 ) = r - s
              zout ( 1, j, nout4 ) = r + s
              r25 = r2 - r5
              r34 = r3 - r4
              s25 = s2 + s5
              s34 = s3 + s4
              zout ( 2, j, nout1 ) = s1 + s25 + s34
              r = cos2 * s25 + cos4 * s34 + s1
              s = sin2 * r25 + sin4 * r34
              zout ( 2, j, nout2 ) = r + s
              zout ( 2, j, nout5 ) = r - s
              r = cos4 * s25 + cos2 * s34 + s1
              s = sin4 * r25 - sin2 * r34
              zout ( 2, j, nout3 ) = r + s
              zout ( 2, j, nout4 ) = r - s
            END DO
          END DO
        ELSE
          nin1 = ia - after
          nout1 = ia - atn
          DO ib = 1, before
            nin1 = nin1 + after
            nin2 = nin1 + atb
            nin3 = nin2 + atb
            nin4 = nin3 + atb
            nin5 = nin4 + atb
            nout1 = nout1 + atn
            nout2 = nout1 + after
            nout3 = nout2 + after
            nout4 = nout3 + after
            nout5 = nout4 + after
            DO j = 1, nfft
              r1 = zin ( 1, j, nin1 )
              s1 = zin ( 2, j, nin1 )
              r = zin ( 1, j, nin2 )
              s = zin ( 2, j, nin2 )
              r2 = ( r + s ) * rt2i
              s2 = ( -r + s ) * rt2i
              r3 = zin ( 2, j, nin3 )
              s3 = -zin ( 1, j, nin3 )
              r = zin ( 1, j, nin4 )
              s = zin ( 2, j, nin4 )
              r4 = ( s - r ) * rt2i
              s4 = - ( r + s ) * rt2i
              r5 = -zin ( 1, j, nin5 )
              s5 = -zin ( 2, j, nin5 )
              r25 = r2 + r5
              r34 = r3 + r4
              s25 = s2 - s5
              s34 = s3 - s4
              zout ( 1, j, nout1 ) = r1 + r25 + r34
              r = cos2 * r25 + cos4 * r34 + r1
              s = sin2 * s25 + sin4 * s34
              zout ( 1, j, nout2 ) = r - s
              zout ( 1, j, nout5 ) = r + s
              r = cos4 * r25 + cos2 * r34 + r1
              s = sin4 * s25 - sin2 * s34
              zout ( 1, j, nout3 ) = r - s
              zout ( 1, j, nout4 ) = r + s
              r25 = r2 - r5
              r34 = r3 - r4
              s25 = s2 + s5
              s34 = s3 + s4
              zout ( 2, j, nout1) = s1 + s25 + s34
              r = cos2 * s25 + cos4 * s34 + s1
              s = sin2 * r25 + sin4 * r34
              zout ( 2, j, nout2 ) = r + s
              zout ( 2, j, nout5 ) = r - s
              r = cos4 * s25 + cos2 * s34 + s1
              s = sin4 * r25 - sin2 * r34
              zout ( 2, j, nout3 ) = r + s
              zout ( 2, j, nout4 ) = r - s
            END DO
          END DO
        END IF
      ELSE
        ias = ia - 1
        itt = ias * before
        itrig = itt + 1
        cr2 = trig ( 1, itrig )
        ci2 = trig ( 2, itrig )
        itrig = itrig + itt
        cr3 = trig ( 1, itrig )
        ci3 = trig ( 2, itrig )
        itrig = itrig + itt
        cr4 = trig ( 1, itrig )
        ci4 = trig ( 2, itrig )
        itrig = itrig + itt
        cr5 = trig ( 1, itrig )
        ci5 = trig ( 2, itrig )
        nin1 = ia - after
        nout1 = ia - atn
        DO ib = 1, before
          nin1 = nin1 + after
          nin2 = nin1 + atb
          nin3 = nin2 + atb
          nin4 = nin3 + atb
          nin5 = nin4 + atb
          nout1 = nout1 + atn
          nout2 = nout1 + after
          nout3 = nout2 + after
          nout4 = nout3 + after
          nout5 = nout4 + after
          DO j = 1, nfft
            r1 = zin ( 1, j, nin1 )
            s1 = zin ( 2, j, nin1 )
            r = zin ( 1, j, nin2 )
            s = zin ( 2, j, nin2 )
            r2 = r * cr2 - s * ci2
            s2 = r * ci2 + s * cr2
            r = zin ( 1, j, nin3 )
            s = zin ( 2, j, nin3 )
            r3 = r * cr3 - s * ci3
            s3 = r * ci3 + s * cr3
            r = zin ( 1, j, nin4 )
            s = zin ( 2, j, nin4 )
            r4 = r * cr4 - s * ci4
            s4 = r * ci4 + s * cr4
            r = zin ( 1, j, nin5 )
            s = zin ( 2, j, nin5 )
            r5 = r * cr5 - s * ci5
            s5 = r * ci5 + s * cr5
            r25 = r2 + r5
            r34 = r3 + r4
            s25 = s2 - s5
            s34 = s3 - s4
            zout ( 1, j, nout1 ) = r1 + r25 + r34
            r = cos2 * r25 + cos4 * r34 + r1
            s = sin2 * s25 + sin4 * s34
            zout ( 1, j, nout2 ) = r - s
            zout ( 1, j, nout5 ) = r + s
            r = cos4 * r25 + cos2 * r34 + r1
            s = sin4 * s25 - sin2 * s34
            zout ( 1, j, nout3 ) = r - s
            zout ( 1, j, nout4 ) = r + s
            r25 = r2 - r5
            r34 = r3 - r4
            s25 = s2 + s5
            s34 = s3 + s4
            zout ( 2, j, nout1 ) = s1 + s25 + s34
            r = cos2 * s25 + cos4 * s34 + s1
            s = sin2 * r25 + sin4 * r34
            zout ( 2, j, nout2 ) = r + s
            zout ( 2, j, nout5 ) = r - s
            r = cos4 * s25 + cos2 * s34 + s1
            s = sin4 * r25 - sin2 * r34
            zout ( 2, j, nout3 ) = r + s
            zout ( 2, j, nout4 ) = r - s
          END DO
        END DO
      END IF
    END DO
  ELSE IF ( now == 6 ) THEN
    bbs = isign * bb
    ia = 1
    nin1 = ia - after
    nout1 = ia - atn
    DO ib = 1, before
      nin1 = nin1 + after
      nin2 = nin1 + atb
      nin3 = nin2 + atb
      nin4 = nin3 + atb
      nin5 = nin4 + atb
      nin6 = nin5 + atb
      nout1 = nout1 + atn
      nout2 = nout1 + after
      nout3 = nout2 + after
      nout4 = nout3 + after
      nout5 = nout4 + after
      nout6 = nout5 + after
      DO j = 1, nfft
        r2 = zin ( 1, j, nin3 )
        s2 = zin ( 2, j, nin3 )
        r3 = zin ( 1, j, nin5 )
        s3 = zin ( 2, j, nin5 )
        r = r2 + r3
        s = s2 + s3
        r1 = zin ( 1, j, nin1 )
        s1 = zin ( 2, j, nin1 )
        ur1 = r + r1
        ui1 = s + s1
        r1 = r1 - 0.5_dp * r
        s1 = s1 - 0.5_dp * s
        r = r2 - r3
        s = s2 - s3
        ur2 = r1 - s * bbs
        ui2 = s1 + r * bbs
        ur3 = r1 + s * bbs
        ui3 = s1 - r * bbs

        r2 = zin ( 1, j, nin6 )
        s2 = zin ( 2, j, nin6 )
        r3 = zin ( 1, j, nin2 )
        s3 = zin ( 2, j, nin2 )
        r = r2 + r3
        s = s2 + s3
        r1 = zin ( 1, j, nin4 )
        s1 = zin ( 2, j, nin4 )
        vr1 = r + r1
        vi1 = s + s1
        r1 = r1 - 0.5_dp * r
        s1 = s1 - 0.5_dp * s
        r = r2 - r3
        s = s2 - s3
        vr2 = r1 - s * bbs
        vi2 = s1 + r * bbs
        vr3 = r1 + s * bbs
        vi3 = s1 - r * bbs

        zout ( 1, j, nout1 ) = ur1 + vr1
        zout ( 2, j, nout1 ) = ui1 + vi1
        zout ( 1, j, nout5 ) = ur2 + vr2
        zout ( 2, j, nout5 ) = ui2 + vi2
        zout ( 1, j, nout3 ) = ur3 + vr3
        zout ( 2, j, nout3 ) = ui3 + vi3
        zout ( 1, j, nout4 ) = ur1 - vr1
        zout ( 2, j, nout4 ) = ui1 - vi1
        zout ( 1, j, nout2 ) = ur2 - vr2
        zout ( 2, j, nout2 ) = ui2 - vi2
        zout ( 1, j, nout6 ) = ur3 - vr3
        zout ( 2, j, nout6 ) = ui3 - vi3
      END DO
    END DO
  ELSE
    STOP 'Error fftstp'
  END IF

!-----------------------------------------------------------------------------!

END SUBROUTINE fftstp

!-----------------------------------------------------------------------------!
