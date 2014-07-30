!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
MODULE ps_wavelet_fft3d

  USE f77_blas
  USE kinds,                           ONLY: dp
#include "../common/cp_common_uses.f90"

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ps_wavelet_fft3d'

  ! longest fft supported, must be equal to the lenght of the ctrig array
  INTEGER, PARAMETER :: ctrig_length=8192

  PUBLIC :: fourier_dim,&
            ctrig,&
            fftstp, ctrig_length

CONTAINS

! *****************************************************************************
!> \brief Give a number n_next > n compatible for the FFT
!> \param n ...
!> \param n_next ...
! *****************************************************************************
  SUBROUTINE fourier_dim(n,n_next)
    INTEGER, INTENT(in)                      :: n
    INTEGER, INTENT(out)                     :: n_next

    INTEGER, PARAMETER                       :: ndata = 149, ndata1024 = 149
    INTEGER, DIMENSION(ndata), PARAMETER :: idata = (/   3,    4,   5,     6, &
      8,    9,   12,   15,   16,   18, 20,   24,   25,   27,   30,   32,   36,&
      40,   45,   48, 54,   60,   64,   72,   75,   80,   81,   90,   96,  100&
      , 108,  120,  125,  128,  135,  144,  150,  160,  162,  180, 192,  200, &
      216,  225,  240,  243,  256,  270,  288,  300, 320,  324,  360,  375,  &
      384,  400,  405,  432,  450,  480, 486,  500,  512,  540,  576,  600,  &
      625,  640,  648,  675, 720,  729,  750,  768,  800,  810,  864,  900,  &
      960,  972, 1000, 1024, 1080, 1125, 1152, 1200, 1215, 1280, 1296, 1350,&
      1440, 1458, 1500, 1536, 1600, 1620, 1728, 1800, 1875, 1920,1944, 2000, &
      2025, 2048, 2160, 2250, 2304, 2400, 2430, 2500,2560, 2592, 2700, 2880, &
      3000, 3072, 3125, 3200, 3240, 3375,3456, 3600, 3750, 3840, 3888, 4000, &
      4050, 4096, 4320, 4500,4608, 4800, 5000, 5120, 5184, 5400, 5625, 5760, &
      6000, 6144,6400, 6480, 6750, 6912, 7200, 7500, 7680, 8000, ctrig_length &
      /)

    INTEGER                                  :: i

!Multiple of 2,3,5

    loop_data: DO i=1,ndata1024
       IF (n <= idata(i)) THEN
          n_next = idata(i)
          RETURN
       END IF
    END DO loop_data
    WRITE(unit=*,fmt=*) "fourier_dim: ",n," is bigger than ",idata(ndata1024)
    STOP
  END SUBROUTINE fourier_dim

!  Copyright (C) Stefan Goedecker, CEA Grenoble, 2002
!  This file is distributed under the terms of the
!  GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .

! --------------------------------------------------------------
!   3-dimensional complex-complex FFT routine:
!   When compared to the best vendor implementations on RISC architectures
!   it gives close to optimal performance (perhaps loosing 20 percent in speed)
!   and it is significanly faster than many not so good vendor implementations
!   as well as other portable FFT's.
!   On all vector machines tested so far (Cray, NEC, Fujitsu) is
!   was significantly faster than the vendor routines
! The theoretical background is described in :
! 1) S. Goedecker: Rotating a three-dimensional array in optimal
! positions for vector processing: Case study for a three-dimensional Fast
! Fourier Transform, Comp. Phys. Commun. \underline{76}, 294 (1993)
! Citing of this reference is greatly appreciated if the routines are used
! for scientific work.

! Presumably good compiler flags:
! IBM, serial power 2: xlf -qarch=pwr2 -O2 -qmaxmem=-1
! with OpenMP: IBM: xlf_r -qfree -O4 -qarch=pwr3 -qtune=pwr3 -qsmp=omp -qmaxmem=-1 ;
!                   a.out
! DEC: f90 -O3 -arch ev67 -pipeline
! with OpenMP: DEC: f90 -O3 -arch ev67 -pipeline -omp -lelan ;
!                   prun -N1 -c4 a.out

!-----------------------------------------------------------

! FFT PART -----------------------------------------------------------------

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param trig ...
!> \param after ...
!> \param before ...
!> \param now ...
!> \param isign ...
!> \param ic ...
! *****************************************************************************
  SUBROUTINE ctrig(n,trig,after,before,now,isign,ic)
!  Copyright (C) Stefan Goedecker, Lausanne, Switzerland, August 1, 1991
!  Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!  Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!  This file is distributed under the terms of the
!  GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .

!     Different factorizations affect the performance
!     Factoring 64 as 4*4*4 might for example be faster on some machines than 8*8.

          IMPLICIT REAL(KIND=dp) (a-h,o-z)
          IMPLICIT INTEGER (i,j,m,n)
          INTEGER after,before
          DIMENSION now(7),after(7),before(7),trig(2,ctrig_length)
          INTEGER, DIMENSION(7,149) :: idata
! The factor 6 is only allowed in the first place!
        DATA ((idata(i,j),i=1,7),j=1,76) /                        &
            3,   3, 1, 1, 1, 1, 1,       4,   4, 1, 1, 1, 1, 1,   &
            5,   5, 1, 1, 1, 1, 1,       6,   6, 1, 1, 1, 1, 1,   &
            8,   8, 1, 1, 1, 1, 1,       9,   3, 3, 1, 1, 1, 1,   &
           12,   4, 3, 1, 1, 1, 1,      15,   5, 3, 1, 1, 1, 1,   &
           16,   4, 4, 1, 1, 1, 1,      18,   6, 3, 1, 1, 1, 1,   &
           20,   5, 4, 1, 1, 1, 1,      24,   8, 3, 1, 1, 1, 1,   &
           25,   5, 5, 1, 1, 1, 1,      27,   3, 3, 3, 1, 1, 1,   &
           30,   6, 5, 1, 1, 1, 1,      32,   8, 4, 1, 1, 1, 1,   &
           36,   4, 3, 3, 1, 1, 1,      40,   8, 5, 1, 1, 1, 1,   &
           45,   5, 3, 3, 1, 1, 1,      48,   4, 4, 3, 1, 1, 1,   &
           54,   6, 3, 3, 1, 1, 1,      60,   5, 4, 3, 1, 1, 1,   &
           64,   8, 8, 1, 1, 1, 1,      72,   8, 3, 3, 1, 1, 1,   &
           75,   5, 5, 3, 1, 1, 1,      80,   5, 4, 4, 1, 1, 1,   &
           81,   3, 3, 3, 3, 1, 1,      90,   6, 5, 3, 1, 1, 1,   &
           96,   8, 4, 3, 1, 1, 1,     100,   5, 5, 4, 1, 1, 1,   &
          108,   4, 3, 3, 3, 1, 1,     120,   8, 5, 3, 1, 1, 1,   &
          125,   5, 5, 5, 1, 1, 1,     128,   8, 4, 4, 1, 1, 1,   &
          135,   5, 3, 3, 3, 1, 1,     144,   6, 8, 3, 1, 1, 1,   &
          150,   6, 5, 5, 1, 1, 1,     160,   8, 5, 4, 1, 1, 1,   &
          162,   6, 3, 3, 3, 1, 1,     180,   5, 4, 3, 3, 1, 1,   &
          192,   6, 8, 4, 1, 1, 1,     200,   8, 5, 5, 1, 1, 1,   &
          216,   8, 3, 3, 3, 1, 1,     225,   5, 5, 3, 3, 1, 1,   &
          240,   6, 8, 5, 1, 1, 1,     243,   3, 3, 3, 3, 3, 1,   &
          256,   8, 8, 4, 1, 1, 1,     270,   6, 5, 3, 3, 1, 1,   &
          288,   8, 4, 3, 3, 1, 1,     300,   5, 5, 4, 3, 1, 1,   &
          320,   5, 4, 4, 4, 1, 1,     324,   4, 3, 3, 3, 3, 1,   &
          360,   8, 5, 3, 3, 1, 1,     375,   5, 5, 5, 3, 1, 1,   &
          384,   8, 4, 4, 3, 1, 1,     400,   5, 5, 4, 4, 1, 1,   &
          405,   5, 3, 3, 3, 3, 1,     432,   4, 4, 3, 3, 3, 1,   &
          450,   6, 5, 5, 3, 1, 1,     480,   8, 5, 4, 3, 1, 1,   &
          486,   6, 3, 3, 3, 3, 1,     500,   5, 5, 5, 4, 1, 1,   &
          512,   8, 8, 8, 1, 1, 1,     540,   5, 4, 3, 3, 3, 1,   &
          576,   4, 4, 4, 3, 3, 1,     600,   8, 5, 5, 3, 1, 1,   &
          625,   5, 5, 5, 5, 1, 1,     640,   8, 5, 4, 4, 1, 1,   &
          648,   8, 3, 3, 3, 3, 1,     675,   5, 5, 3, 3, 3, 1,   &
          720,   5, 4, 4, 3, 3, 1,     729,   3, 3, 3, 3, 3, 3,   &
          750,   6, 5, 5, 5, 1, 1,     768,   4, 4, 4, 4, 3, 1,   &
          800,   8, 5, 5, 4, 1, 1,     810,   6, 5, 3, 3, 3, 1    /
        DATA ((idata(i,j),i=1,7),j=77,149) /                       &
          864,   8, 4, 3, 3, 3, 1,     900,   5, 5, 4, 3, 3, 1,   &
          960,   5, 4, 4, 4, 3, 1,     972,   4, 3, 3, 3, 3, 3,   &
         1000,   8, 5, 5, 5, 1, 1,    1024,   4, 4, 4, 4, 4, 1,   &
         1080,   6, 5, 4, 3, 3, 1,    1125,   5, 5, 5, 3, 3, 1,   &
         1152,   6, 4, 4, 4, 3, 1,    1200,   6, 8, 5, 5, 1, 1,   &
         1215,   5, 3, 3, 3, 3, 3,    1280,   8, 8, 5, 4, 1, 1,   &
         1296,   6, 8, 3, 3, 3, 1,    1350,   6, 5, 5, 3, 3, 1,   &
         1440,   6, 5, 4, 4, 3, 1,    1458,   6, 3, 3, 3, 3, 3,   &
         1500,   5, 5, 5, 4, 3, 1,    1536,   6, 8, 8, 4, 1, 1,   &
         1600,   8, 8, 5, 5, 1, 1,    1620,   5, 4, 3, 3, 3, 3,   &
         1728,   6, 8, 4, 3, 3, 1,    1800,   6, 5, 5, 4, 3, 1,   &
         1875,   5, 5, 5, 5, 3, 1,    1920,   6, 5, 4, 4, 4, 1,   &
         1944,   6, 4, 3, 3, 3, 3,    2000,   5, 5, 5, 4, 4, 1,   &
         2025,   5, 5, 3, 3, 3, 3,    2048,   8, 4, 4, 4, 4, 1,   &
         2160,   6, 8, 5, 3, 3, 1,    2250,   6, 5, 5, 5, 3, 1,   &
         2304,   6, 8, 4, 4, 3, 1,    2400,   6, 5, 5, 4, 4, 1,   &
         2430,   6, 5, 3, 3, 3, 3,    2500,   5, 5, 5, 5, 4, 1,   &
         2560,   8, 5, 4, 4, 4, 1,    2592,   6, 4, 4, 3, 3, 3,   &
         2700,   5, 5, 4, 3, 3, 3,    2880,   6, 8, 5, 4, 3, 1,   &
         3000,   6, 5, 5, 5, 4, 1,    3072,   6, 8, 4, 4, 4, 1,   &
         3125,   5, 5, 5, 5, 5, 1,    3200,   8, 5, 5, 4, 4, 1,   &
         3240,   6, 5, 4, 3, 3, 3,    3375,   5, 5, 5, 3, 3, 3,   &
         3456,   6, 4, 4, 4, 3, 3,    3600,   6, 8, 5, 5, 3, 1,   &
         3750,   6, 5, 5, 5, 5, 1,    3840,   6, 8, 5, 4, 4, 1,   &
         3888,   6, 8, 3, 3, 3, 3,    4000,   8, 5, 5, 5, 4, 1,   &
         4050,   6, 5, 5, 3, 3, 3,    4096,   8, 8, 4, 4, 4, 1,   &
         4320,   6, 5, 4, 4, 3, 3,    4500,   5, 5, 5, 4, 3, 3,   &
         4608,   6, 8, 8, 4, 3, 1,    4800,   6, 8, 5, 5, 4, 1,   &
         5000,   8, 5, 5, 5, 5, 1,    5120,   8, 8, 5, 4, 4, 1,   &
         5184,   6, 8, 4, 3, 3, 3,    5400,   6, 5, 5, 4, 3, 3,   &
         5625,   5, 5, 5, 5, 3, 3,    5760,   6, 8, 8, 5, 3, 1,   &
         6000,   6, 8, 5, 5, 5, 1,    6144,   6, 8, 8, 4, 4, 1,   &
         6400,   8, 8, 5, 5, 4, 1,    6480,   6, 8, 5, 3, 3, 3,   &
         6750,   6, 5, 5, 5, 3, 3,    6912,   6, 8, 4, 4, 3, 3,   &
         7200,   6, 5, 5, 4, 4, 3,    7500,   5, 5, 5, 5, 4, 3,   &
         7680,   6, 8, 8, 5, 4, 1,    8000,   8, 8, 5, 5, 5, 1,   &
         8192,   8, 8, 8, 4, 4, 1 /

        DO 111 i=1,149
        IF (n.eq.idata(1,i)) THEN
        ic=0
        DO 11 j=1,6
        itt=idata(1+j,i)
        IF (itt.gt.1) THEN
        ic=ic+1
        now(j)=idata(1+j,i)
        ELSE
        GOTO 1000
        ENDIF
11        CONTINUE
        GOTO 1000
        ENDIF
111        CONTINUE
        PRINT*,'VALUE OF',n,'NOT ALLOWED FOR FFT, ALLOWED VALUES ARE:'
37        FORMAT(15(i5))
        WRITE(*,37) (idata(1,i),i=1,149)
        STOP
1000        CONTINUE

        after(1)=1
        before(ic)=1
        DO i=2,ic
          after(i)=after(i-1)*now(i-1)
          before(ic-i+1)=before(ic-i+2)*now(ic-i+2)
        ENDDO

        twopi=6.283185307179586_dp
        angle=isign*twopi/n
        IF (MOD(n,2).eq.0) THEN
        nh=n/2
        trig(1,1)=1._dp
        trig(2,1)=0._dp
        trig(1,nh+1)=-1._dp
        trig(2,nh+1)=0._dp
        DO 40,i=1,nh-1
        trigc=COS(i*angle)
        trigs=SIN(i*angle)
        trig(1,i+1)=trigc
        trig(2,i+1)=trigs
        trig(1,n-i+1)=trigc
        trig(2,n-i+1)=-trigs
40      CONTINUE
        ELSE
        nh=(n-1)/2
        trig(1,1)=1._dp
        trig(2,1)=0._dp
        DO 20,i=1,nh
        trigc=COS(i*angle)
        trigs=SIN(i*angle)
        trig(1,i+1)=trigc
        trig(2,i+1)=trigs
        trig(1,n-i+1)=trigc
        trig(2,n-i+1)=-trigs
20      CONTINUE
        ENDIF

      END SUBROUTINE ctrig

!ccccccccccccccccccccccccccccccccccccccccccccccc

! *****************************************************************************
!> \brief ...
!> \param mm ...
!> \param nfft ...
!> \param m ...
!> \param nn ...
!> \param n ...
!> \param zin ...
!> \param zout ...
!> \param trig ...
!> \param after ...
!> \param now ...
!> \param before ...
!> \param isign ...
! *****************************************************************************
      SUBROUTINE fftstp(mm,nfft,m,nn,n,zin,zout,trig,after,now,before,isign)
!  Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!  Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1995, 1999
!  This file is distributed under the terms of the
!  GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .

        IMPLICIT REAL(KIND=dp) (a-h,o-z)
        IMPLICIT INTEGER (i,j,m,n)
        INTEGER after,before,atn,atb
        DIMENSION trig(2,ctrig_length),zin(2,mm,m),zout(2,nn,n)
        atn=after*now
        atb=after*before

!         sqrt(.5_dp)
        rt2i=0.7071067811865475_dp
        IF (now.eq.2) THEN
        ia=1
        nin1=ia-after
        nout1=ia-atn
        DO ib=1,before
        nin1=nin1+after
        nin2=nin1+atb
        nout1=nout1+atn
        nout2=nout1+after
        DO j=1,nfft
        r1=zin(1,j,nin1)
        s1=zin(2,j,nin1)
        r2=zin(1,j,nin2)
        s2=zin(2,j,nin2)
        zout(1,j,nout1)= r2 + r1
        zout(2,j,nout1)= s2 + s1
        zout(1,j,nout2)= r1 - r2
        zout(2,j,nout2)= s1 - s2
        ENDDO ; ENDDO
        DO 2000,ia=2,after
        ias=ia-1
        IF (2*ias.eq.after) THEN
                IF (isign.eq.1) THEN
                        nin1=ia-after
                        nout1=ia-atn
                        DO ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        DO j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r2=zin(2,j,nin2)
                        s2=zin(1,j,nin2)
                        zout(1,j,nout1)= r1 - r2
                        zout(2,j,nout1)= s2 + s1
                        zout(1,j,nout2)= r2 + r1
                        zout(2,j,nout2)= s1 - s2
                        ENDDO ; ENDDO
                ELSE
                        nin1=ia-after
                        nout1=ia-atn
                        DO ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        DO j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r2=zin(2,j,nin2)
                        s2=zin(1,j,nin2)
                        zout(1,j,nout1)= r2 + r1
                        zout(2,j,nout1)= s1 - s2
                        zout(1,j,nout2)= r1 - r2
                        zout(2,j,nout2)= s2 + s1
                        ENDDO ; ENDDO
                ENDIF
        ELSE IF (4*ias.eq.after) THEN
                IF (isign.eq.1) THEN
                        nin1=ia-after
                        nout1=ia-atn
                        DO ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        DO j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=(r - s)*rt2i
                        s2=(r + s)*rt2i
                        zout(1,j,nout1)= r2 + r1
                        zout(2,j,nout1)= s2 + s1
                        zout(1,j,nout2)= r1 - r2
                        zout(2,j,nout2)= s1 - s2
                        ENDDO ; ENDDO
                ELSE
                        nin1=ia-after
                        nout1=ia-atn
                        DO ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        DO j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=(r + s)*rt2i
                        s2=(s - r)*rt2i
                        zout(1,j,nout1)= r2 + r1
                        zout(2,j,nout1)= s2 + s1
                        zout(1,j,nout2)= r1 - r2
                        zout(2,j,nout2)= s1 - s2
                        ENDDO ; ENDDO
                ENDIF
        ELSE IF (4*ias.eq.3*after) THEN
                IF (isign.eq.1) THEN
                        nin1=ia-after
                        nout1=ia-atn
                        DO ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        DO j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=(r + s)*rt2i
                        s2=(r - s)*rt2i
                        zout(1,j,nout1)= r1 - r2
                        zout(2,j,nout1)= s2 + s1
                        zout(1,j,nout2)= r2 + r1
                        zout(2,j,nout2)= s1 - s2
                        ENDDO ; ENDDO
                ELSE
                        nin1=ia-after
                        nout1=ia-atn
                        DO ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        DO j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=(s - r)*rt2i
                        s2=(r + s)*rt2i
                        zout(1,j,nout1)= r2 + r1
                        zout(2,j,nout1)= s1 - s2
                        zout(1,j,nout2)= r1 - r2
                        zout(2,j,nout2)= s2 + s1
                        ENDDO ; ENDDO
                ENDIF
        ELSE
                itrig=ias*before+1
                cr2=trig(1,itrig)
                ci2=trig(2,itrig)
                nin1=ia-after
                nout1=ia-atn
                DO ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nout1=nout1+atn
                nout2=nout1+after
                DO j=1,nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r=zin(1,j,nin2)
                s=zin(2,j,nin2)
                r2=r*cr2 - s*ci2
                s2=r*ci2 + s*cr2
                zout(1,j,nout1)= r2 + r1
                zout(2,j,nout1)= s2 + s1
                zout(1,j,nout2)= r1 - r2
                zout(2,j,nout2)= s1 - s2
                ENDDO ; ENDDO
        ENDIF
2000        CONTINUE
        ELSE IF (now.eq.4) THEN
        IF (isign.eq.1) THEN
                ia=1
                nin1=ia-after
                nout1=ia-atn
                DO ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nin4=nin3+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                nout4=nout3+after
                DO j=1,nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r2=zin(1,j,nin2)
                s2=zin(2,j,nin2)
                r3=zin(1,j,nin3)
                s3=zin(2,j,nin3)
                r4=zin(1,j,nin4)
                s4=zin(2,j,nin4)
                r=r1 + r3
                s=r2 + r4
                zout(1,j,nout1) = r + s
                zout(1,j,nout3) = r - s
                r=r1 - r3
                s=s2 - s4
                zout(1,j,nout2) = r - s
                zout(1,j,nout4) = r + s
                r=s1 + s3
                s=s2 + s4
                zout(2,j,nout1) = r + s
                zout(2,j,nout3) = r - s
                r=s1 - s3
                s=r2 - r4
                zout(2,j,nout2) = r + s
                zout(2,j,nout4) = r - s
                ENDDO ; ENDDO
                DO 4000,ia=2,after
                ias=ia-1
                IF (2*ias.eq.after) THEN
                        nin1=ia-after
                        nout1=ia-atn
                        DO ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        DO j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=(r-s)*rt2i
                        s2=(r+s)*rt2i
                        r3=zin(2,j,nin3)
                        s3=zin(1,j,nin3)
                        r=zin(1,j,nin4)
                        s=zin(2,j,nin4)
                        r4=(r + s)*rt2i
                        s4=(r - s)*rt2i
                        r=r1 - r3
                        s=r2 - r4
                        zout(1,j,nout1) = r + s
                        zout(1,j,nout3) = r - s
                        r=r1 + r3
                        s=s2 - s4
                        zout(1,j,nout2) = r - s
                        zout(1,j,nout4) = r + s
                        r=s1 + s3
                        s=s2 + s4
                        zout(2,j,nout1) = r + s
                        zout(2,j,nout3) = r - s
                        r=s1 - s3
                        s=r2 + r4
                        zout(2,j,nout2) = r + s
                        zout(2,j,nout4) = r - s
                        ENDDO ; ENDDO
                ELSE
                        itt=ias*before
                        itrig=itt+1
                        cr2=trig(1,itrig)
                        ci2=trig(2,itrig)
                        itrig=itrig+itt
                        cr3=trig(1,itrig)
                        ci3=trig(2,itrig)
                        itrig=itrig+itt
                        cr4=trig(1,itrig)
                        ci4=trig(2,itrig)
                        nin1=ia-after
                        nout1=ia-atn
                        DO ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        DO j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=r*cr2 - s*ci2
                        s2=r*ci2 + s*cr2
                        r=zin(1,j,nin3)
                        s=zin(2,j,nin3)
                        r3=r*cr3 - s*ci3
                        s3=r*ci3 + s*cr3
                        r=zin(1,j,nin4)
                        s=zin(2,j,nin4)
                        r4=r*cr4 - s*ci4
                        s4=r*ci4 + s*cr4
                        r=r1 + r3
                        s=r2 + r4
                        zout(1,j,nout1) = r + s
                        zout(1,j,nout3) = r - s
                        r=r1 - r3
                        s=s2 - s4
                        zout(1,j,nout2) = r - s
                        zout(1,j,nout4) = r + s
                        r=s1 + s3
                        s=s2 + s4
                        zout(2,j,nout1) = r + s
                        zout(2,j,nout3) = r - s
                        r=s1 - s3
                        s=r2 - r4
                        zout(2,j,nout2) = r + s
                        zout(2,j,nout4) = r - s
                        ENDDO ; ENDDO
                ENDIF
4000                CONTINUE
        ELSE
                ia=1
                nin1=ia-after
                nout1=ia-atn
                DO ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nin4=nin3+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                nout4=nout3+after
                DO j=1,nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r2=zin(1,j,nin2)
                s2=zin(2,j,nin2)
                r3=zin(1,j,nin3)
                s3=zin(2,j,nin3)
                r4=zin(1,j,nin4)
                s4=zin(2,j,nin4)
                r=r1 + r3
                s=r2 + r4
                zout(1,j,nout1) = r + s
                zout(1,j,nout3) = r - s
                r=r1 - r3
                s=s2 - s4
                zout(1,j,nout2) = r + s
                zout(1,j,nout4) = r - s
                r=s1 + s3
                s=s2 + s4
                zout(2,j,nout1) = r + s
                zout(2,j,nout3) = r - s
                r=s1 - s3
                s=r2 - r4
                zout(2,j,nout2) = r - s
                zout(2,j,nout4) = r + s
                ENDDO ; ENDDO
                DO 4100,ia=2,after
                ias=ia-1
                IF (2*ias.eq.after) THEN
                        nin1=ia-after
                        nout1=ia-atn
                        DO ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        DO j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=(r + s)*rt2i
                        s2=(s - r)*rt2i
                        r3=zin(2,j,nin3)
                        s3=zin(1,j,nin3)
                        r=zin(1,j,nin4)
                        s=zin(2,j,nin4)
                        r4=(s - r)*rt2i
                        s4=(r + s)*rt2i
                        r=r1 + r3
                        s=r2 + r4
                        zout(1,j,nout1) = r + s
                        zout(1,j,nout3) = r - s
                        r=r1 - r3
                        s=s2 + s4
                        zout(1,j,nout2) = r + s
                        zout(1,j,nout4) = r - s
                        r=s1 - s3
                        s=s2 - s4
                        zout(2,j,nout1) = r + s
                        zout(2,j,nout3) = r - s
                        r=s1 + s3
                        s=r2 - r4
                        zout(2,j,nout2) = r - s
                        zout(2,j,nout4) = r + s
                        ENDDO ; ENDDO
                ELSE
                        itt=ias*before
                        itrig=itt+1
                        cr2=trig(1,itrig)
                        ci2=trig(2,itrig)
                        itrig=itrig+itt
                        cr3=trig(1,itrig)
                        ci3=trig(2,itrig)
                        itrig=itrig+itt
                        cr4=trig(1,itrig)
                        ci4=trig(2,itrig)
                        nin1=ia-after
                        nout1=ia-atn
                        DO ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        DO j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=r*cr2 - s*ci2
                        s2=r*ci2 + s*cr2
                        r=zin(1,j,nin3)
                        s=zin(2,j,nin3)
                        r3=r*cr3 - s*ci3
                        s3=r*ci3 + s*cr3
                        r=zin(1,j,nin4)
                        s=zin(2,j,nin4)
                        r4=r*cr4 - s*ci4
                        s4=r*ci4 + s*cr4
                        r=r1 + r3
                        s=r2 + r4
                        zout(1,j,nout1) = r + s
                        zout(1,j,nout3) = r - s
                        r=r1 - r3
                        s=s2 - s4
                        zout(1,j,nout2) = r + s
                        zout(1,j,nout4) = r - s
                        r=s1 + s3
                        s=s2 + s4
                        zout(2,j,nout1) = r + s
                        zout(2,j,nout3) = r - s
                        r=s1 - s3
                        s=r2 - r4
                        zout(2,j,nout2) = r - s
                        zout(2,j,nout4) = r + s
                        ENDDO ; ENDDO
                ENDIF
4100                CONTINUE
        ENDIF
        ELSE IF (now.eq.8) THEN
        IF (isign.eq.-1) THEN
                ia=1
                        nin1=ia-after
                        nout1=ia-atn
                        DO ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nin5=nin4+atb
                        nin6=nin5+atb
                        nin7=nin6+atb
                        nin8=nin7+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        nout5=nout4+after
                        nout6=nout5+after
                        nout7=nout6+after
                        nout8=nout7+after
                        DO j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r2=zin(1,j,nin2)
                        s2=zin(2,j,nin2)
                        r3=zin(1,j,nin3)
                        s3=zin(2,j,nin3)
                        r4=zin(1,j,nin4)
                        s4=zin(2,j,nin4)
                        r5=zin(1,j,nin5)
                        s5=zin(2,j,nin5)
                        r6=zin(1,j,nin6)
                        s6=zin(2,j,nin6)
                        r7=zin(1,j,nin7)
                        s7=zin(2,j,nin7)
                        r8=zin(1,j,nin8)
                        s8=zin(2,j,nin8)
                        r=r1 + r5
                        s=r3 + r7
                        ap=r + s
                        am=r - s
                        r=r2 + r6
                        s=r4 + r8
                        bp=r + s
                        bm=r - s
                        r=s1 + s5
                        s=s3 + s7
                        cp=r + s
                        cm=r - s
                        r=s2 + s6
                        s=s4 + s8
                        dpp=r + s
                        dm=r - s
                        zout(1,j,nout1) = ap + bp
                        zout(2,j,nout1) = cp + dpp
                        zout(1,j,nout5) = ap - bp
                        zout(2,j,nout5) = cp - dpp
                        zout(1,j,nout3) = am + dm
                        zout(2,j,nout3) = cm - bm
                        zout(1,j,nout7) = am - dm
                        zout(2,j,nout7) = cm + bm
                        r=r1 - r5
                        s=s3 - s7
                        ap=r + s
                        am=r - s
                        r=s1 - s5
                        s=r3 - r7
                        bp=r + s
                        bm=r - s
                        r=s4 - s8
                        s=r2 - r6
                        cp=r + s
                        cm=r - s
                        r=s2 - s6
                        s=r4 - r8
                        dpp=r + s
                        dm=r - s
                        r = ( cp + dm)*rt2i
                        s = ( dm - cp)*rt2i
                        cp= ( cm + dpp)*rt2i
                        dpp = ( cm - dpp)*rt2i
                        zout(1,j,nout2) = ap + r
                        zout(2,j,nout2) = bm + s
                        zout(1,j,nout6) = ap - r
                        zout(2,j,nout6) = bm - s
                        zout(1,j,nout4) = am + cp
                        zout(2,j,nout4) = bp + dpp
                        zout(1,j,nout8) = am - cp
                        zout(2,j,nout8) = bp - dpp
                        ENDDO ; ENDDO
                DO 8000,ia=2,after
                ias=ia-1
                        itt=ias*before
                        itrig=itt+1
                        cr2=trig(1,itrig)
                        ci2=trig(2,itrig)
                        itrig=itrig+itt
                        cr3=trig(1,itrig)
                        ci3=trig(2,itrig)
                        itrig=itrig+itt
                        cr4=trig(1,itrig)
                        ci4=trig(2,itrig)
                        itrig=itrig+itt
                        cr5=trig(1,itrig)
                        ci5=trig(2,itrig)
                        itrig=itrig+itt
                        cr6=trig(1,itrig)
                        ci6=trig(2,itrig)
                        itrig=itrig+itt
                        cr7=trig(1,itrig)
                        ci7=trig(2,itrig)
                        itrig=itrig+itt
                        cr8=trig(1,itrig)
                        ci8=trig(2,itrig)
                        nin1=ia-after
                        nout1=ia-atn
                        DO ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nin5=nin4+atb
                        nin6=nin5+atb
                        nin7=nin6+atb
                        nin8=nin7+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        nout5=nout4+after
                        nout6=nout5+after
                        nout7=nout6+after
                        nout8=nout7+after
                        DO j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=r*cr2 - s*ci2
                        s2=r*ci2 + s*cr2
                        r=zin(1,j,nin3)
                        s=zin(2,j,nin3)
                        r3=r*cr3 - s*ci3
                        s3=r*ci3 + s*cr3
                        r=zin(1,j,nin4)
                        s=zin(2,j,nin4)
                        r4=r*cr4 - s*ci4
                        s4=r*ci4 + s*cr4
                        r=zin(1,j,nin5)
                        s=zin(2,j,nin5)
                        r5=r*cr5 - s*ci5
                        s5=r*ci5 + s*cr5
                        r=zin(1,j,nin6)
                        s=zin(2,j,nin6)
                        r6=r*cr6 - s*ci6
                        s6=r*ci6 + s*cr6
                        r=zin(1,j,nin7)
                        s=zin(2,j,nin7)
                        r7=r*cr7 - s*ci7
                        s7=r*ci7 + s*cr7
                        r=zin(1,j,nin8)
                        s=zin(2,j,nin8)
                        r8=r*cr8 - s*ci8
                        s8=r*ci8 + s*cr8
                        r=r1 + r5
                        s=r3 + r7
                        ap=r + s
                        am=r - s
                        r=r2 + r6
                        s=r4 + r8
                        bp=r + s
                        bm=r - s
                        r=s1 + s5
                        s=s3 + s7
                        cp=r + s
                        cm=r - s
                        r=s2 + s6
                        s=s4 + s8
                        dpp=r + s
                        dm=r - s
                        zout(1,j,nout1) = ap + bp
                        zout(2,j,nout1) = cp + dpp
                        zout(1,j,nout5) = ap - bp
                        zout(2,j,nout5) = cp - dpp
                        zout(1,j,nout3) = am + dm
                        zout(2,j,nout3) = cm - bm
                        zout(1,j,nout7) = am - dm
                        zout(2,j,nout7) = cm + bm
                        r=r1 - r5
                        s=s3 - s7
                        ap=r + s
                        am=r - s
                        r=s1 - s5
                        s=r3 - r7
                        bp=r + s
                        bm=r - s
                        r=s4 - s8
                        s=r2 - r6
                        cp=r + s
                        cm=r - s
                        r=s2 - s6
                        s=r4 - r8
                        dpp=r + s
                        dm=r - s
                        r = ( cp + dm)*rt2i
                        s = ( dm - cp)*rt2i
                        cp= ( cm + dpp)*rt2i
                        dpp = ( cm - dpp)*rt2i
                        zout(1,j,nout2) = ap + r
                        zout(2,j,nout2) = bm + s
                        zout(1,j,nout6) = ap - r
                        zout(2,j,nout6) = bm - s
                        zout(1,j,nout4) = am + cp
                        zout(2,j,nout4) = bp + dpp
                        zout(1,j,nout8) = am - cp
                        zout(2,j,nout8) = bp - dpp 
                        ENDDO ; ENDDO
8000                CONTINUE

        ELSE
                ia=1
                        nin1=ia-after
                        nout1=ia-atn
                        DO ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nin5=nin4+atb
                        nin6=nin5+atb
                        nin7=nin6+atb
                        nin8=nin7+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        nout5=nout4+after
                        nout6=nout5+after
                        nout7=nout6+after
                        nout8=nout7+after
                        DO j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r2=zin(1,j,nin2)
                        s2=zin(2,j,nin2)
                        r3=zin(1,j,nin3)
                        s3=zin(2,j,nin3)
                        r4=zin(1,j,nin4)
                        s4=zin(2,j,nin4)
                        r5=zin(1,j,nin5)
                        s5=zin(2,j,nin5)
                        r6=zin(1,j,nin6)
                        s6=zin(2,j,nin6)
                        r7=zin(1,j,nin7)
                        s7=zin(2,j,nin7)
                        r8=zin(1,j,nin8)
                        s8=zin(2,j,nin8)
                        r=r1 + r5
                        s=r3 + r7
                        ap=r + s
                        am=r - s
                        r=r2 + r6
                        s=r4 + r8
                        bp=r + s
                        bm=r - s
                        r=s1 + s5
                        s=s3 + s7
                        cp=r + s
                        cm=r - s
                        r=s2 + s6
                        s=s4 + s8
                        dpp=r + s
                        dm=r - s
                        zout(1,j,nout1) = ap + bp
                        zout(2,j,nout1) = cp + dpp
                        zout(1,j,nout5) = ap - bp
                        zout(2,j,nout5) = cp - dpp
                        zout(1,j,nout3) = am - dm
                        zout(2,j,nout3) = cm + bm
                        zout(1,j,nout7) = am + dm
                        zout(2,j,nout7) = cm - bm
                        r= r1 - r5
                        s=-s3 + s7
                        ap=r + s
                        am=r - s
                        r=s1 - s5
                        s=r7 - r3
                        bp=r + s
                        bm=r - s
                        r=-s4 + s8
                        s= r2 - r6
                        cp=r + s
                        cm=r - s
                        r=-s2 + s6
                        s= r4 - r8
                        dpp=r + s
                        dm=r - s
                        r = ( cp + dm)*rt2i
                        s = ( cp - dm)*rt2i
                        cp= ( cm + dpp)*rt2i
                        dpp= ( dpp - cm)*rt2i
                        zout(1,j,nout2) = ap + r
                        zout(2,j,nout2) = bm + s
                        zout(1,j,nout6) = ap - r
                        zout(2,j,nout6) = bm - s
                        zout(1,j,nout4) = am + cp
                        zout(2,j,nout4) = bp + dpp
                        zout(1,j,nout8) = am - cp
                        zout(2,j,nout8) = bp - dpp
                        ENDDO ; ENDDO

                DO 8001,ia=2,after
                ias=ia-1
                        itt=ias*before
                        itrig=itt+1
                        cr2=trig(1,itrig)
                        ci2=trig(2,itrig)
                        itrig=itrig+itt
                        cr3=trig(1,itrig)
                        ci3=trig(2,itrig)
                        itrig=itrig+itt
                        cr4=trig(1,itrig)
                        ci4=trig(2,itrig)
                        itrig=itrig+itt
                        cr5=trig(1,itrig)
                        ci5=trig(2,itrig)
                        itrig=itrig+itt
                        cr6=trig(1,itrig)
                        ci6=trig(2,itrig)
                        itrig=itrig+itt
                        cr7=trig(1,itrig)
                        ci7=trig(2,itrig)
                        itrig=itrig+itt
                        cr8=trig(1,itrig)
                        ci8=trig(2,itrig)
                        nin1=ia-after
                        nout1=ia-atn
                        DO ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nin5=nin4+atb
                        nin6=nin5+atb
                        nin7=nin6+atb
                        nin8=nin7+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        nout5=nout4+after
                        nout6=nout5+after
                        nout7=nout6+after
                        nout8=nout7+after
                        DO j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=r*cr2 - s*ci2
                        s2=r*ci2 + s*cr2
                        r=zin(1,j,nin3)
                        s=zin(2,j,nin3)
                        r3=r*cr3 - s*ci3
                        s3=r*ci3 + s*cr3
                        r=zin(1,j,nin4)
                        s=zin(2,j,nin4)
                        r4=r*cr4 - s*ci4
                        s4=r*ci4 + s*cr4
                        r=zin(1,j,nin5)
                        s=zin(2,j,nin5)
                        r5=r*cr5 - s*ci5
                        s5=r*ci5 + s*cr5
                        r=zin(1,j,nin6)
                        s=zin(2,j,nin6)
                        r6=r*cr6 - s*ci6
                        s6=r*ci6 + s*cr6
                        r=zin(1,j,nin7)
                        s=zin(2,j,nin7)
                        r7=r*cr7 - s*ci7
                        s7=r*ci7 + s*cr7
                        r=zin(1,j,nin8)
                        s=zin(2,j,nin8)
                        r8=r*cr8 - s*ci8
                        s8=r*ci8 + s*cr8
                        r=r1 + r5
                        s=r3 + r7
                        ap=r + s
                        am=r - s
                        r=r2 + r6
                        s=r4 + r8
                        bp=r + s
                        bm=r - s
                        r=s1 + s5
                        s=s3 + s7
                        cp=r + s
                        cm=r - s
                        r=s2 + s6
                        s=s4 + s8
                        dpp=r + s
                        dm=r - s
                        zout(1,j,nout1) = ap + bp
                        zout(2,j,nout1) = cp + dpp
                        zout(1,j,nout5) = ap - bp
                        zout(2,j,nout5) = cp - dpp
                        zout(1,j,nout3) = am - dm
                        zout(2,j,nout3) = cm + bm
                        zout(1,j,nout7) = am + dm
                        zout(2,j,nout7) = cm - bm
                        r= r1 - r5
                        s=-s3 + s7
                        ap=r + s
                        am=r - s
                        r=s1 - s5
                        s=r7 - r3
                        bp=r + s
                        bm=r - s
                        r=-s4 + s8
                        s= r2 - r6
                        cp=r + s
                        cm=r - s
                        r=-s2 + s6
                        s= r4 - r8
                        dpp=r + s
                        dm=r - s
                        r = ( cp + dm)*rt2i
                        s = ( cp - dm)*rt2i
                        cp= ( cm + dpp)*rt2i
                        dpp= ( dpp - cm)*rt2i
                        zout(1,j,nout2) = ap + r
                        zout(2,j,nout2) = bm + s
                        zout(1,j,nout6) = ap - r
                        zout(2,j,nout6) = bm - s
                        zout(1,j,nout4) = am + cp
                        zout(2,j,nout4) = bp + dpp
                        zout(1,j,nout8) = am - cp
                        zout(2,j,nout8) = bp - dpp
                        ENDDO ; ENDDO
8001                CONTINUE

        ENDIF
        ELSE IF (now.eq.3) THEN
!         .5_dp*sqrt(3._dp)
        bb=isign*0.8660254037844387_dp
        ia=1
        nin1=ia-after
        nout1=ia-atn
        DO ib=1,before
        nin1=nin1+after
        nin2=nin1+atb
        nin3=nin2+atb
        nout1=nout1+atn
        nout2=nout1+after
        nout3=nout2+after
        DO j=1,nfft
        r1=zin(1,j,nin1)
        s1=zin(2,j,nin1)
        r2=zin(1,j,nin2)
        s2=zin(2,j,nin2)
        r3=zin(1,j,nin3)
        s3=zin(2,j,nin3)
        r=r2 + r3
        s=s2 + s3
        zout(1,j,nout1) = r + r1
        zout(2,j,nout1) = s + s1
        r1=r1 - .5_dp*r
        s1=s1 - .5_dp*s
        r2=bb*(r2-r3)
        s2=bb*(s2-s3)
        zout(1,j,nout2) = r1 - s2
        zout(2,j,nout2) = s1 + r2
        zout(1,j,nout3) = r1 + s2
        zout(2,j,nout3) = s1 - r2
        ENDDO ; ENDDO
        DO 3000,ia=2,after
        ias=ia-1
        IF (4*ias.eq.3*after) THEN
        IF (isign.eq.1) THEN
                nin1=ia-after
                nout1=ia-atn
                DO ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                DO j=1,nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r2=zin(2,j,nin2)
                s2=zin(1,j,nin2)
                r3=zin(1,j,nin3)
                s3=zin(2,j,nin3)
                r=r3 + r2
                s=s2 - s3
                zout(1,j,nout1) = r1 - r
                zout(2,j,nout1) = s + s1
                r1=r1 + .5_dp*r
                s1=s1 - .5_dp*s
                r2=bb*(r2-r3)
                s2=bb*(s2+s3)
                zout(1,j,nout2) = r1 - s2
                zout(2,j,nout2) = s1 - r2
                zout(1,j,nout3) = r1 + s2
                zout(2,j,nout3) = s1 + r2
                ENDDO ; ENDDO
        ELSE
                nin1=ia-after
                nout1=ia-atn
                DO ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                DO j=1,nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r2=zin(2,j,nin2)
                s2=zin(1,j,nin2)
                r3=zin(1,j,nin3)
                s3=zin(2,j,nin3)
                r=r2 - r3
                s=s2 + s3
                zout(1,j,nout1) = r + r1
                zout(2,j,nout1) = s1 - s
                r1=r1 - .5_dp*r
                s1=s1 + .5_dp*s
                r2=bb*(r2+r3)
                s2=bb*(s2-s3)
                zout(1,j,nout2) = r1 + s2
                zout(2,j,nout2) = s1 + r2
                zout(1,j,nout3) = r1 - s2
                zout(2,j,nout3) = s1 - r2
                ENDDO ; ENDDO
        ENDIF
        ELSE IF (8*ias.eq.3*after) THEN
        IF (isign.eq.1) THEN
                nin1=ia-after
                nout1=ia-atn
                DO ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                DO j=1,nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r=zin(1,j,nin2)
                s=zin(2,j,nin2)
                r2=(r - s)*rt2i
                s2=(r + s)*rt2i
                r3=zin(2,j,nin3)
                s3=zin(1,j,nin3)
                r=r2 - r3
                s=s2 + s3
                zout(1,j,nout1) = r + r1
                zout(2,j,nout1) = s + s1
                r1=r1 - .5_dp*r
                s1=s1 - .5_dp*s
                r2=bb*(r2+r3)
                s2=bb*(s2-s3)
                zout(1,j,nout2) = r1 - s2
                zout(2,j,nout2) = s1 + r2
                zout(1,j,nout3) = r1 + s2
                zout(2,j,nout3) = s1 - r2
                ENDDO ; ENDDO
        ELSE
                nin1=ia-after
                nout1=ia-atn
                DO ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                DO j=1,nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r=zin(1,j,nin2)
                s=zin(2,j,nin2)
                r2=(r + s)*rt2i
                s2=(s - r)*rt2i
                r3=zin(2,j,nin3)
                s3=zin(1,j,nin3)
                r=r2 + r3
                s=s2 - s3
                zout(1,j,nout1) = r + r1
                zout(2,j,nout1) = s + s1
                r1=r1 - .5_dp*r
                s1=s1 - .5_dp*s
                r2=bb*(r2-r3)
                s2=bb*(s2+s3)
                zout(1,j,nout2) = r1 - s2
                zout(2,j,nout2) = s1 + r2
                zout(1,j,nout3) = r1 + s2
                zout(2,j,nout3) = s1 - r2
                ENDDO ; ENDDO
        ENDIF
        ELSE
        itt=ias*before
        itrig=itt+1
        cr2=trig(1,itrig)
        ci2=trig(2,itrig)
        itrig=itrig+itt
        cr3=trig(1,itrig)
        ci3=trig(2,itrig)
        nin1=ia-after
        nout1=ia-atn
        DO ib=1,before
        nin1=nin1+after
        nin2=nin1+atb
        nin3=nin2+atb
        nout1=nout1+atn
        nout2=nout1+after
        nout3=nout2+after
        DO j=1,nfft
        r1=zin(1,j,nin1)
        s1=zin(2,j,nin1)
        r=zin(1,j,nin2)
        s=zin(2,j,nin2)
        r2=r*cr2 - s*ci2
        s2=r*ci2 + s*cr2
        r=zin(1,j,nin3)
        s=zin(2,j,nin3)
        r3=r*cr3 - s*ci3
        s3=r*ci3 + s*cr3
        r=r2 + r3
        s=s2 + s3
        zout(1,j,nout1) = r + r1
        zout(2,j,nout1) = s + s1
        r1=r1 - .5_dp*r
        s1=s1 - .5_dp*s
        r2=bb*(r2-r3)
        s2=bb*(s2-s3)
        zout(1,j,nout2) = r1 - s2
        zout(2,j,nout2) = s1 + r2
        zout(1,j,nout3) = r1 + s2
        zout(2,j,nout3) = s1 - r2
        ENDDO ; ENDDO
        ENDIF
3000        CONTINUE
        ELSE IF (now.eq.5) THEN
!         cos(2._dp*pi/5._dp)
        cos2=0.3090169943749474_dp
!         cos(4._dp*pi/5._dp)
        cos4=-0.8090169943749474_dp
!        sin(2._dp*pi/5._dp)
        sin2=isign*0.9510565162951536_dp
!         sin(4._dp*pi/5._dp)
        sin4=isign*0.5877852522924731_dp
        ia=1
        nin1=ia-after
        nout1=ia-atn
        DO ib=1,before
        nin1=nin1+after
        nin2=nin1+atb
        nin3=nin2+atb
        nin4=nin3+atb
        nin5=nin4+atb
        nout1=nout1+atn
        nout2=nout1+after
        nout3=nout2+after
        nout4=nout3+after
        nout5=nout4+after
        DO j=1,nfft
        r1=zin(1,j,nin1)
        s1=zin(2,j,nin1)
        r2=zin(1,j,nin2)
        s2=zin(2,j,nin2)
        r3=zin(1,j,nin3)
        s3=zin(2,j,nin3)
        r4=zin(1,j,nin4)
        s4=zin(2,j,nin4)
        r5=zin(1,j,nin5)
        s5=zin(2,j,nin5)
        r25 = r2 + r5
        r34 = r3 + r4
        s25 = s2 - s5
        s34 = s3 - s4
        zout(1,j,nout1) = r1 + r25 + r34
        r = r1 + cos2*r25 + cos4*r34
        s = sin2*s25 + sin4*s34
        zout(1,j,nout2) = r - s
        zout(1,j,nout5) = r + s
        r = r1 + cos4*r25 + cos2*r34
        s = sin4*s25 - sin2*s34
        zout(1,j,nout3) = r - s
        zout(1,j,nout4) = r + s
        r25 = r2 - r5
        r34 = r3 - r4
        s25 = s2 + s5
        s34 = s3 + s4
        zout(2,j,nout1) = s1 + s25 + s34
        r = s1 + cos2*s25 + cos4*s34
        s = sin2*r25 + sin4*r34
        zout(2,j,nout2) = r + s
        zout(2,j,nout5) = r - s
        r = s1 + cos4*s25 + cos2*s34
        s = sin4*r25 - sin2*r34
        zout(2,j,nout3) = r + s
        zout(2,j,nout4) = r - s
        ENDDO ; ENDDO
        DO 5000,ia=2,after
        ias=ia-1
        IF (8*ias.eq.5*after) THEN
                IF (isign.eq.1) THEN
                        nin1=ia-after
                        nout1=ia-atn
                        DO ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nin5=nin4+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        nout5=nout4+after
                        DO j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=(r - s)*rt2i
                        s2=(r + s)*rt2i
                        r3=zin(2,j,nin3)
                        s3=zin(1,j,nin3)
                        r=zin(1,j,nin4)
                        s=zin(2,j,nin4)
                        r4=(r + s)*rt2i
                        s4=(r - s)*rt2i
                        r5=zin(1,j,nin5)
                        s5=zin(2,j,nin5)
                        r25 = r2 - r5
                        r34 = r3 + r4
                        s25 = s2 + s5
                        s34 = s3 - s4
                        zout(1,j,nout1) = r1 + r25 - r34
                        r = r1 + cos2*r25 - cos4*r34
                        s = sin2*s25 + sin4*s34
                        zout(1,j,nout2) = r - s
                        zout(1,j,nout5) = r + s
                        r = r1 + cos4*r25 - cos2*r34
                        s = sin4*s25 - sin2*s34
                        zout(1,j,nout3) = r - s
                        zout(1,j,nout4) = r + s
                        r25 = r2 + r5
                        r34 = r4 - r3
                        s25 = s2 - s5
                        s34 = s3 + s4
                        zout(2,j,nout1) = s1 + s25 + s34
                        r = s1 + cos2*s25 + cos4*s34
                        s = sin2*r25 + sin4*r34
                        zout(2,j,nout2) = r + s
                        zout(2,j,nout5) = r - s
                        r = s1 + cos4*s25 + cos2*s34
                        s = sin4*r25 - sin2*r34
                        zout(2,j,nout3) = r + s
                        zout(2,j,nout4) = r - s
                        ENDDO ; ENDDO
                ELSE
                        nin1=ia-after
                        nout1=ia-atn
                        DO ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nin5=nin4+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        nout5=nout4+after
                        DO j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=(r + s)*rt2i
                        s2=(s - r)*rt2i
                        r3=zin(2,j,nin3)
                        s3=zin(1,j,nin3)
                        r=zin(1,j,nin4)
                        s=zin(2,j,nin4)
                        r4=(s - r)*rt2i
                        s4=(r + s)*rt2i
                        r5=zin(1,j,nin5)
                        s5=zin(2,j,nin5)
                        r25 = r2 - r5
                        r34 = r3 + r4
                        s25 = s2 + s5
                        s34 = s4 - s3
                        zout(1,j,nout1) = r1 + r25 + r34
                        r = r1 + cos2*r25 + cos4*r34
                        s = sin2*s25 + sin4*s34
                        zout(1,j,nout2) = r - s
                        zout(1,j,nout5) = r + s
                        r = r1 + cos4*r25 + cos2*r34
                        s = sin4*s25 - sin2*s34
                        zout(1,j,nout3) = r - s
                        zout(1,j,nout4) = r + s
                        r25 = r2 + r5
                        r34 = r3 - r4
                        s25 = s2 - s5
                        s34 = s3 + s4
                        zout(2,j,nout1) = s1 + s25 - s34
                        r = s1 + cos2*s25 - cos4*s34
                        s = sin2*r25 + sin4*r34
                        zout(2,j,nout2) = r + s
                        zout(2,j,nout5) = r - s
                        r = s1 + cos4*s25 - cos2*s34
                        s = sin4*r25 - sin2*r34
                        zout(2,j,nout3) = r + s
                        zout(2,j,nout4) = r - s
                        ENDDO ; ENDDO
                ENDIF
        ELSE
                ias=ia-1
                itt=ias*before
                itrig=itt+1
                cr2=trig(1,itrig)
                ci2=trig(2,itrig)
                itrig=itrig+itt
                cr3=trig(1,itrig)
                ci3=trig(2,itrig)
                itrig=itrig+itt
                cr4=trig(1,itrig)
                ci4=trig(2,itrig)
                itrig=itrig+itt
                cr5=trig(1,itrig)
                ci5=trig(2,itrig)
                nin1=ia-after
                nout1=ia-atn
                DO ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nin4=nin3+atb
                nin5=nin4+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                nout4=nout3+after
                nout5=nout4+after
                DO j=1,nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r=zin(1,j,nin2)
                s=zin(2,j,nin2)
                r2=r*cr2 - s*ci2
                s2=r*ci2 + s*cr2
                r=zin(1,j,nin3)
                s=zin(2,j,nin3)
                r3=r*cr3 - s*ci3
                s3=r*ci3 + s*cr3
                r=zin(1,j,nin4)
                s=zin(2,j,nin4)
                r4=r*cr4 - s*ci4
                s4=r*ci4 + s*cr4
                r=zin(1,j,nin5)
                s=zin(2,j,nin5)
                r5=r*cr5 - s*ci5
                s5=r*ci5 + s*cr5
                r25 = r2 + r5
                r34 = r3 + r4
                s25 = s2 - s5
                s34 = s3 - s4
                zout(1,j,nout1) = r1 + r25 + r34
                r = r1 + cos2*r25 + cos4*r34
                s = sin2*s25 + sin4*s34
                zout(1,j,nout2) = r - s
                zout(1,j,nout5) = r + s
                r = r1 + cos4*r25 + cos2*r34
                s = sin4*s25 - sin2*s34
                zout(1,j,nout3) = r - s
                zout(1,j,nout4) = r + s
                r25 = r2 - r5
                r34 = r3 - r4
                s25 = s2 + s5
                s34 = s3 + s4
                zout(2,j,nout1) = s1 + s25 + s34
                r = s1 + cos2*s25 + cos4*s34
                s = sin2*r25 + sin4*r34
                zout(2,j,nout2) = r + s
                zout(2,j,nout5) = r - s
                r = s1 + cos4*s25 + cos2*s34
                s = sin4*r25 - sin2*r34
                zout(2,j,nout3) = r + s
                zout(2,j,nout4) = r - s
                ENDDO ; ENDDO
        ENDIF
5000        CONTINUE
       ELSE IF (now.eq.6) THEN
!         .5_dp*sqrt(3._dp)
        bb=isign*0.8660254037844387_dp

        ia=1
        nin1=ia-after
        nout1=ia-atn
        DO ib=1,before
        nin1=nin1+after
        nin2=nin1+atb
        nin3=nin2+atb
        nin4=nin3+atb
        nin5=nin4+atb
        nin6=nin5+atb
        nout1=nout1+atn
        nout2=nout1+after
        nout3=nout2+after
        nout4=nout3+after
        nout5=nout4+after
        nout6=nout5+after
        DO j=1,nfft
        r2=zin(1,j,nin3)
        s2=zin(2,j,nin3)
        r3=zin(1,j,nin5)
        s3=zin(2,j,nin5)
        r=r2 + r3
        s=s2 + s3
        r1=zin(1,j,nin1)
        s1=zin(2,j,nin1)
        ur1 = r + r1
        ui1 = s + s1
        r1=r1 - .5_dp*r
        s1=s1 - .5_dp*s
        r=r2-r3
        s=s2-s3
        ur2 = r1 - s*bb
        ui2 = s1 + r*bb
        ur3 = r1 + s*bb
        ui3 = s1 - r*bb

        r2=zin(1,j,nin6)
        s2=zin(2,j,nin6)
        r3=zin(1,j,nin2)
        s3=zin(2,j,nin2)
        r=r2 + r3
        s=s2 + s3
        r1=zin(1,j,nin4)
        s1=zin(2,j,nin4)
        vr1 = r + r1
        vi1 = s + s1
        r1=r1 - .5_dp*r
        s1=s1 - .5_dp*s
        r=r2-r3
        s=s2-s3
        vr2 = r1 - s*bb
        vi2 = s1 + r*bb
        vr3 = r1 + s*bb
        vi3 = s1 - r*bb

        zout(1,j,nout1)=ur1+vr1
        zout(2,j,nout1)=ui1+vi1
        zout(1,j,nout5)=ur2+vr2
        zout(2,j,nout5)=ui2+vi2
        zout(1,j,nout3)=ur3+vr3
        zout(2,j,nout3)=ui3+vi3
        zout(1,j,nout4)=ur1-vr1
        zout(2,j,nout4)=ui1-vi1
        zout(1,j,nout2)=ur2-vr2
        zout(2,j,nout2)=ui2-vi2
        zout(1,j,nout6)=ur3-vr3
        zout(2,j,nout6)=ui3-vi3
        ENDDO ; ENDDO
     ELSE
        STOP 'error fftstp'
     ENDIF

   END SUBROUTINE fftstp

END MODULE ps_wavelet_fft3d
