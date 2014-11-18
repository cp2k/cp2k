!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Creates the wavelet kernel for the wavelet based poisson solver.
!> \author Florian Schiffmann (09.2007,fschiff)
! *****************************************************************************
MODULE ps_wavelet_base
  
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_alltoall
  USE ps_wavelet_fft3d,                ONLY: ctrig,&
                                             ctrig_length,&
                                             fftstp
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ps_wavelet_base'

  PUBLIC :: scramble_unpack, p_poissonsolver, s_poissonsolver, f_poissonsolver

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param n1 ...
!> \param n2 ...
!> \param n3 ...
!> \param nd1 ...
!> \param nd2 ...
!> \param nd3 ...
!> \param md1 ...
!> \param md2 ...
!> \param md3 ...
!> \param nproc ...
!> \param iproc ...
!> \param zf& ...
! *****************************************************************************
  SUBROUTINE P_PoissonSolver(n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,nproc,iproc,zf&
       ,scal,hx,hy,hz,mpi_group,error)
    INTEGER, INTENT(in)                      :: n1, n2, n3, nd1, nd2, nd3, &
                                                md1, md2, md3, nproc, iproc
    REAL(KIND=dp), &
      DIMENSION(md1, md3, md2/nproc), &
      INTENT(inout)                          :: zf
    REAL(KIND=dp), INTENT(in)                :: scal, hx, hy, hz
    INTEGER, INTENT(in)                      :: mpi_group
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER, PARAMETER                       :: ncache_optimal = 8*1024, &
                                                timing_flag = 0

    INTEGER :: i, i1, i3, i_all, i_stat, ic1, ic2, ic3, inzee, j, j2, J2stb, &
      J2stf, j3, Jp2stb, Jp2stf, lot, lzt, ma, mb, ncache, nfft
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: after1, after2, after3, &
                                                before1, before2, before3, &
                                                now1, now2, now3
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: btrig1, btrig2, btrig3, &
                                                ftrig1, ftrig2, ftrig3
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: zt, zw
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: zmpi2
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :, :)               :: zmpi1

    IF (nd1.LT.n1/2+1) STOP 'Parallel convolution:ERROR:nd1'
    IF (nd2.LT.n2/2+1) STOP 'Parallel convolution:ERROR:nd2'
    IF (nd3.LT.n3/2+1) STOP 'Parallel convolution:ERROR:nd3'
    IF (md1.LT.n1) STOP 'Parallel convolution:ERROR:md1'
    IF (md2.LT.n2) STOP 'Parallel convolution:ERROR:md2'
    IF (md3.LT.n3) STOP 'Parallel convolution:ERROR:md3'
    IF (MOD(nd3,nproc).NE.0) STOP 'Parallel convolution:ERROR:nd3'
    IF (MOD(md2,nproc).NE.0) STOP 'Parallel convolution:ERROR:md2'

    !defining work arrays dimensions
    ncache=ncache_optimal
    IF (ncache <= MAX(n1,n2,n3)*4) ncache=MAX(n1,n2,n3)*4

    lzt=n2
    IF (MOD(n2,2).EQ.0) lzt=lzt+1
    IF (MOD(n2,4).EQ.0) lzt=lzt+1 !maybe this is useless

    !Allocations
    ALLOCATE(btrig1(2,ctrig_length),stat=i_all)
    ALLOCATE(ftrig1(2,ctrig_length),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(after1(7),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(now1(7),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(before1(7),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(btrig2(2,ctrig_length),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(ftrig2(2,ctrig_length),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(after2(7),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(now2(7),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(before2(7),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(btrig3(2,ctrig_length),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(ftrig3(2,ctrig_length),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(after3(7),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(now3(7),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(before3(7),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(zw(2,ncache/4,2),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(zt(2,lzt,n1),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(zmpi2(2,n1,md2/nproc,nd3),stat=i_stat)
    i_all=i_all+i_stat
    IF (nproc.GT.1) ALLOCATE(zmpi1(2,n1,md2/nproc,nd3/nproc,nproc),stat=i_stat)
    IF (i_all+i_stat /= 0) THEN
       WRITE(*,*)' p_poissonsolver: problem of memory allocation'
       STOP
    END IF

    !calculating the FFT work arrays (beware on the HalFFT in n3 dimension)
    CALL ctrig(n3,btrig3,after3,before3,now3,1,ic3)
    CALL ctrig(n1,btrig1,after1,before1,now1,1,ic1)
    CALL ctrig(n2,btrig2,after2,before2,now2,1,ic2)
    DO  j=1,n1
       ftrig1(1,j)= btrig1(1,j)
       ftrig1(2,j)=-btrig1(2,j)
    ENDDO
    DO  j=1,n2
       ftrig2(1,j)= btrig2(1,j)
       ftrig2(2,j)=-btrig2(2,j)
    ENDDO
    DO  j=1,n3
       ftrig3(1,j)= btrig3(1,j)
       ftrig3(2,j)=-btrig3(2,j)
    ENDDO

    ! transform along z axis
    lot=ncache/(4*n3)
    IF (lot.LT.1) THEN
       WRITE(*,*) &
            'convolxc_off:ncache has to be enlarged to be able to hold at' // &
            'least one 1-d FFT of this size even though this will' // &
            'reduce the performance for shorter transform lengths'
       STOP
    ENDIF
    DO j2=1,md2/nproc
       !this condition ensures that we manage only the interesting part for the FFT
       IF (iproc*(md2/nproc)+j2.LE.n2) THEN
          DO i1=1,n1,lot
             ma=i1
             mb=MIN(i1+(lot-1),n1)
             nfft=mb-ma+1
             !inserting real data into complex array of half lenght
             CALL P_fill_upcorn(md1,md3,lot,nfft,n3,zf(i1,1,j2),zw(1,1,1))

             !performing FFT
             !input: I1,I3,J2,(Jp2)
             inzee=1
             DO i=1,ic3
                CALL fftstp(lot,nfft,n3,lot,n3,zw(1,1,inzee),zw(1,1,3-inzee), &
                     btrig3,after3(i),now3(i),before3(i),1)
                inzee=3-inzee
             ENDDO

             !output: I1,i3,J2,(Jp2)
             !exchanging components
             !input: I1,i3,J2,(Jp2)
             CALL scramble_P(i1,j2,lot,nfft,n1,n3,md2,nproc,nd3,zw(1,1,inzee),zmpi2)
             !output: I1,J2,i3,(Jp2)
          END DO
       END IF
    END DO

    !Interprocessor data transposition
    !input: I1,J2,j3,jp3,(Jp2)
    IF (nproc.GT.1) THEN
       !communication scheduling
       CALL mp_alltoall(zmpi2,zmpi1,2*n1*(md2/nproc)*(nd3/nproc),mpi_group)
    ENDIF
    !output: I1,J2,j3,Jp2,(jp3)

    !now each process perform complete convolution of its planes
    DO j3=1,nd3/nproc
       !this condition ensures that we manage only the interesting part for the FFT
       IF (iproc*(nd3/nproc)+j3.LE.n3/2+1) THEN
          Jp2stb=1
          J2stb=1
          Jp2stf=1
          J2stf=1

          ! transform along x axis
          lot=ncache/(4*n1)
          IF (lot.LT.1) THEN
             WRITE(*,*) &
                  'convolxc_off:ncache has to be enlarged to be able to hold at' // &
                  'least one 1-d FFT of this size even though this will' // &
                  'reduce the performance for shorter transform lengths'
             STOP
          ENDIF

          DO j=1,n2,lot
             ma=j
             mb=MIN(j+(lot-1),n2)
             nfft=mb-ma+1

             !reverse index ordering, leaving the planes to be transformed at the end
             !input: I1,J2,j3,Jp2,(jp3)
             IF (nproc.EQ.1) THEN
                CALL P_mpiswitch_upcorn(j3,nfft,Jp2stb,J2stb,lot,n1,md2,nd3,nproc,zmpi2,zw(1,1,1))
             ELSE
                CALL P_mpiswitch_upcorn(j3,nfft,Jp2stb,J2stb,lot,n1,md2,nd3,nproc,zmpi1,zw(1,1,1))
             ENDIF
             !output: J2,Jp2,I1,j3,(jp3)

             !performing FFT
             !input: I2,I1,j3,(jp3)
             inzee=1
             DO i=1,ic1-1
                CALL fftstp(lot,nfft,n1,lot,n1,zw(1,1,inzee),zw(1,1,3-inzee), &
                     btrig1,after1(i),now1(i),before1(i),1)
                inzee=3-inzee
             ENDDO

             !storing the last step into zt array
             i=ic1
             CALL fftstp(lot,nfft,n1,lzt,n1,zw(1,1,inzee),zt(1,j,1), &
                  btrig1,after1(i),now1(i),before1(i),1)
             !output: I2,i1,j3,(jp3)
          END DO

          !transform along y axis
          lot=ncache/(4*n2)
          IF (lot.LT.1) THEN
             WRITE(*,*) &
                  'convolxc_off:ncache has to be enlarged to be able to hold at' // &
                  'least one 1-d FFT of this size even though this will' // &
                  'reduce the performance for shorter transform lengths'
             STOP
          ENDIF

          DO j=1,n1,lot
             ma=j
             mb=MIN(j+(lot-1),n1)
             nfft=mb-ma+1

             !reverse ordering
             !input: I2,i1,j3,(jp3)
             CALL P_switch_upcorn(nfft,n2,lot,n1,lzt,zt(1,1,j),zw(1,1,1))
             !output: i1,I2,j3,(jp3)

             !performing FFT
             !input: i1,I2,j3,(jp3)
             inzee=1
             DO i=1,ic2
                CALL fftstp(lot,nfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee), &
                     btrig2,after2(i),now2(i),before2(i),1)
                inzee=3-inzee
             ENDDO
             !output: i1,i2,j3,(jp3)

             !Multiply with kernel in fourier space
             i3=iproc*(nd3/nproc)+j3
             CALL P_multkernel(n1,n2,n3,lot,nfft,j,i3,zw(1,1,inzee),hx,hy,hz)

             !TRANSFORM BACK IN REAL SPACE

             !transform along y axis
             !input: i1,i2,j3,(jp3)
             DO i=1,ic2
                CALL fftstp(lot,nfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee), &
                     ftrig2,after2(i),now2(i),before2(i),-1)
                inzee=3-inzee
             ENDDO

             !reverse ordering
             !input: i1,I2,j3,(jp3)
             CALL P_unswitch_downcorn(nfft,n2,lot,n1,lzt,zw(1,1,inzee),zt(1,1,j))
             !output: I2,i1,j3,(jp3)
          END DO

          !transform along x axis
          !input: I2,i1,j3,(jp3)
          lot=ncache/(4*n1)
          DO j=1,n2,lot
             ma=j
             mb=MIN(j+(lot-1),n2)
             nfft=mb-ma+1

             !performing FFT
             i=1
             CALL fftstp(lzt,nfft,n1,lot,n1,zt(1,j,1),zw(1,1,1), &
                  ftrig1,after1(i),now1(i),before1(i),-1)

             inzee=1
             DO i=2,ic1
                CALL fftstp(lot,nfft,n1,lot,n1,zw(1,1,inzee),zw(1,1,3-inzee), &
                     ftrig1,after1(i),now1(i),before1(i),-1)
                inzee=3-inzee
             ENDDO
             !output: I2,I1,j3,(jp3)

             !reverse ordering
             !input: J2,Jp2,I1,j3,(jp3)
             IF (nproc.EQ.1) THEN
                CALL P_unmpiswitch_downcorn(j3,nfft,Jp2stf,J2stf,lot,n1,md2,nd3,nproc,zw(1,1,inzee),zmpi2)
             ELSE
                CALL P_unmpiswitch_downcorn(j3,nfft,Jp2stf,J2stf,lot,n1,md2,nd3,nproc,zw(1,1,inzee),zmpi1)
             ENDIF
             ! output: I1,J2,j3,Jp2,(jp3)
          END DO
       ENDIF
    END DO

    !Interprocessor data transposition
    !input: I1,J2,j3,Jp2,(jp3)
    IF (nproc.GT.1) THEN
       !communication scheduling
       CALL mp_alltoall(zmpi1,zmpi2,2*n1*(md2/nproc)*(nd3/nproc),mpi_group)
    ENDIF
    !output: I1,J2,j3,jp3,(Jp2)
    !transform along z axis
    !input: I1,J2,i3,(Jp2)
    lot=ncache/(4*n3)
    DO j2=1,md2/nproc
       !this condition ensures that we manage only the interesting part for the FFT
       IF (iproc*(md2/nproc)+j2.LE.n2) THEN
          DO i1=1,n1,lot
             ma=i1
             mb=MIN(i1+(lot-1),n1)
             nfft=mb-ma+1

             !reverse ordering
             !input: I1,J2,i3,(Jp2)
             CALL unscramble_P(i1,j2,lot,nfft,n1,n3,md2,nproc,nd3,zmpi2,zw(1,1,1))
             !output: I1,i3,J2,(Jp2)

             !performing FFT
             !input: I1,i3,J2,(Jp2)
             inzee=1
             DO i=1,ic3
                CALL fftstp(lot,nfft,n3,lot,n3,zw(1,1,inzee),zw(1,1,3-inzee), &
                     ftrig3,after3(i),now3(i),before3(i),-1)
                inzee=3-inzee
             ENDDO
             !output: I1,I3,J2,(Jp2)

             !rebuild the output array
             CALL P_unfill_downcorn(md1,md3,lot,nfft,n3,zw(1,1,inzee),zf(i1,1,j2),scal)

          END DO
       ENDIF
    END DO

    !De-allocations
    DEALLOCATE(btrig1,stat=i_all)
    DEALLOCATE(ftrig1,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(after1,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(now1,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(before1,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(btrig2,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(ftrig2,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(after2,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(now2,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(before2,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(btrig3,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(ftrig3,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(after3,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(now3,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(before3,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(zmpi2,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(zw,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(zt,stat=i_stat)
    i_all=i_all+i_stat
    IF (nproc.GT.1) DEALLOCATE(zmpi1,stat=i_stat)
    IF (i_all+i_stat /= 0) THEN
       WRITE(*,*)' p_poissonsolver: problem of memory deallocation'
       STOP
    END IF
    !  call timing(iproc,'PSolv_comput  ','OF')
  END SUBROUTINE P_PoissonSolver

! *****************************************************************************
!> \brief ...
!> \param j3 ...
!> \param nfft ...
!> \param Jp2stb ...
!> \param J2stb ...
!> \param lot ...
!> \param n1 ...
!> \param md2 ...
!> \param nd3 ...
!> \param nproc ...
!> \param zmpi1 ...
!> \param zw ...
! *****************************************************************************
  SUBROUTINE P_mpiswitch_upcorn(j3,nfft,Jp2stb,J2stb,lot,n1,md2,nd3,nproc,zmpi1,zw)
    INTEGER, INTENT(in)                      :: j3, nfft
    INTEGER, INTENT(inout)                   :: Jp2stb, J2stb
    INTEGER, INTENT(in)                      :: lot, n1, md2, nd3, nproc
    REAL(KIND=dp), DIMENSION(2, n1, md2/&
      nproc, nd3/nproc, nproc), INTENT(in)   :: zmpi1
    REAL(KIND=dp), DIMENSION(2, lot, n1), &
      INTENT(inout)                          :: zw

    INTEGER                                  :: I1, J2, Jp2, mfft

    mfft=0
    DO Jp2=Jp2stb,nproc
       DO J2=J2stb,md2/nproc
          mfft=mfft+1
          IF (mfft.GT.nfft) THEN
             Jp2stb=Jp2
             J2stb=J2
             RETURN
          END IF
          DO I1=1,n1
             zw(1,mfft,I1)=zmpi1(1,I1,J2,j3,Jp2)
             zw(2,mfft,I1)=zmpi1(2,I1,J2,j3,Jp2)
          END DO
       END DO
       J2stb=1
    END DO
  END SUBROUTINE P_mpiswitch_upcorn

! *****************************************************************************
!> \brief ...
!> \param nfft ...
!> \param n2 ...
!> \param lot ...
!> \param n1 ...
!> \param lzt ...
!> \param zt ...
!> \param zw ...
! *****************************************************************************
  SUBROUTINE P_switch_upcorn(nfft,n2,lot,n1,lzt,zt,zw)
    INTEGER, INTENT(in)                      :: nfft, n2, lot, n1, lzt
    REAL(KIND=dp), DIMENSION(2, lzt, n1), &
      INTENT(in)                             :: zt
    REAL(KIND=dp), DIMENSION(2, lot, n2), &
      INTENT(inout)                          :: zw

    INTEGER                                  :: i, j

    DO j=1,nfft
       DO i=1,n2
          zw(1,j,i)=zt(1,i,j)
          zw(2,j,i)=zt(2,i,j)
       END DO
    END DO

  END SUBROUTINE P_switch_upcorn

! *****************************************************************************
!> \brief ...
!> \param nfft ...
!> \param n2 ...
!> \param lot ...
!> \param n1 ...
!> \param lzt ...
!> \param zw ...
!> \param zt ...
! *****************************************************************************
  SUBROUTINE P_unswitch_downcorn(nfft,n2,lot,n1,lzt,zw,zt)
    INTEGER, INTENT(in)                      :: nfft, n2, lot, n1, lzt
    REAL(KIND=dp), DIMENSION(2, lot, n2), &
      INTENT(in)                             :: zw
    REAL(KIND=dp), DIMENSION(2, lzt, n1), &
      INTENT(inout)                          :: zt

    INTEGER                                  :: i, j

    DO j=1,nfft
       DO i=1,n2
          zt(1,i,j)=zw(1,j,i)
          zt(2,i,j)=zw(2,j,i)
       END DO
    END DO

  END SUBROUTINE P_unswitch_downcorn

! *****************************************************************************
!> \brief ...
!> \param j3 ...
!> \param nfft ...
!> \param Jp2stf ...
!> \param J2stf ...
!> \param lot ...
!> \param n1 ...
!> \param md2 ...
!> \param nd3 ...
!> \param nproc ...
!> \param zw ...
!> \param zmpi1 ...
! *****************************************************************************
  SUBROUTINE P_unmpiswitch_downcorn(j3,nfft,Jp2stf,J2stf,lot,n1,md2,nd3,nproc,zw,zmpi1)
    INTEGER, INTENT(in)                      :: j3, nfft
    INTEGER, INTENT(inout)                   :: Jp2stf, J2stf
    INTEGER, INTENT(in)                      :: lot, n1, md2, nd3, nproc
    REAL(KIND=dp), DIMENSION(2, lot, n1), &
      INTENT(in)                             :: zw
    REAL(KIND=dp), DIMENSION(2, n1, md2/&
      nproc, nd3/nproc, nproc), &
      INTENT(inout)                          :: zmpi1

    INTEGER                                  :: I1, J2, Jp2, mfft

    mfft=0
    DO Jp2=Jp2stf,nproc
       DO J2=J2stf,md2/nproc
          mfft=mfft+1
          IF (mfft.GT.nfft) THEN
             Jp2stf=Jp2
             J2stf=J2
             RETURN
          END IF
          DO I1=1,n1
             zmpi1(1,I1,J2,j3,Jp2)=zw(1,mfft,I1)
             zmpi1(2,I1,J2,j3,Jp2)=zw(2,mfft,I1)
          END DO
       END DO
       J2stf=1
    END DO
  END SUBROUTINE P_unmpiswitch_downcorn

! *****************************************************************************
!> \brief (Based on suitable modifications of S.Goedecker routines)
!>      Restore data into output array
!> \param md1 Dimensions of the undistributed part of the real grid
!> \param md3 Dimensions of the undistributed part of the real grid
!> \param lot ...
!> \param nfft number of planes
!> \param n3 (twice the) dimension of the last FFTtransform
!> \param zw FFT work array
!> \param zf Original distributed density as well as
!>                   Distributed solution of the poisson equation (inout)
!> \param scal Needed to achieve unitarity and correct dimensions
!> \date February 2006
!> \author S. Goedecker, L. Genovese
!> \note Assuming that high frequencies are in the corners 
!>      and that n3 is multiple of 4 
!> 
!>  RESTRICTIONS on USAGE 
!>      Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994 
!>      Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999 
!>      Copyright (C) 2002 Stefan Goedecker, CEA Grenoble 
!>      This file is distributed under the terms of the 
!>      GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt . 
! *****************************************************************************
  SUBROUTINE P_unfill_downcorn(md1,md3,lot,nfft,n3,zw,zf,scal)
    INTEGER, INTENT(in)                      :: md1, md3, lot, nfft, n3
    REAL(KIND=dp), DIMENSION(2, lot, n3), &
      INTENT(in)                             :: zw
    REAL(KIND=dp), DIMENSION(md1, md3), &
      INTENT(inout)                          :: zf
    REAL(KIND=dp), INTENT(in)                :: scal

    INTEGER                                  :: i1, i3
    REAL(KIND=dp)                            :: pot1

    DO i3=1,n3
       DO i1=1,nfft
          pot1 = scal*zw(1,i1,i3)
          zf(i1,i3)= pot1
       END DO
    END DO

  END SUBROUTINE P_unfill_downcorn

! *****************************************************************************
!> \brief ...
!> \param md1 ...
!> \param md3 ...
!> \param lot ...
!> \param nfft ...
!> \param n3 ...
!> \param zf ...
!> \param zw ...
! *****************************************************************************
  SUBROUTINE P_fill_upcorn(md1,md3,lot,nfft,n3,zf,zw)
    INTEGER, INTENT(in)                      :: md1, md3, lot, nfft, n3
    REAL(KIND=dp), DIMENSION(md1, md3), &
      INTENT(in)                             :: zf
    REAL(KIND=dp), DIMENSION(2, lot, n3), &
      INTENT(inout)                          :: zw

    INTEGER                                  :: i1, i3

    DO i3=1,n3
       DO i1=1,nfft
          zw(1,i1,i3)=zf(i1,i3)
          zw(2,i1,i3)=0._dp
       END DO
    END DO

  END SUBROUTINE P_fill_upcorn

! *****************************************************************************
!> \brief (Based on suitable modifications of S.Goedecker routines)
!>      Assign the correct planes to the work array zmpi2
!>      in order to prepare for interprocessor data transposition.
!> \param i1 Starting points of the plane and number of remaining lines
!> \param j2 Starting points of the plane and number of remaining lines
!> \param lot Starting points of the plane and number of remaining lines
!> \param nfft Starting points of the plane and number of remaining lines
!> \param n1 logical dimension of the FFT transform, reference for work arrays
!> \param n3 logical dimension of the FFT transform, reference for work arrays
!> \param md2 Dimensions of real grid
!> \param nproc ...
!> \param nd3 Dimensions of the kernel
!> \param zw Work array (input)
!> \param zmpi2 Work array for multiprocessor manipulation (output)
!> \date February 2006
!> \author S. Goedecker, L. Genovese
!> \note
!>  RESTRICTIONS on USAGE
!>      Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!>      Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!>      Copyright (C) 2002 Stefan Goedecker, CEA Grenoble
!>      This file is distributed under the terms of the
!>      GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt . 
! *****************************************************************************
  SUBROUTINE scramble_P(i1,j2,lot,nfft,n1,n3,md2,nproc,nd3,zw,zmpi2)
    INTEGER, INTENT(in)                      :: i1, j2, lot, nfft, n1, n3, &
                                                md2, nproc, nd3
    REAL(KIND=dp), DIMENSION(2, lot, n3), &
      INTENT(in)                             :: zw
    REAL(KIND=dp), &
      DIMENSION(2, n1, md2/nproc, nd3), &
      INTENT(inout)                          :: zmpi2

    INTEGER                                  :: i, i3

    DO i3=1,n3/2+1
       DO i=0,nfft-1
          zmpi2(1,i1+i,j2,i3)=zw(1,i+1,i3)
          zmpi2(2,i1+i,j2,i3)=zw(2,i+1,i3)
       END DO
    END DO

  END SUBROUTINE scramble_P

! *****************************************************************************
!> \brief (Based on suitable modifications of S.Goedecker routines)
!>      Insert the correct planes of the work array zmpi2
!>      in order to prepare for backward FFT transform
!> \param i1 Starting points of the plane and number of remaining lines
!> \param j2 Starting points of the plane and number of remaining lines
!> \param lot Starting points of the plane and number of remaining lines
!> \param nfft Starting points of the plane and number of remaining lines
!> \param n1 logical dimension of the FFT transform, reference for work arrays
!> \param n3 logical dimension of the FFT transform, reference for work arrays
!> \param md2 Dimensions of real grid
!> \param nproc ...
!> \param nd3 Dimensions of the kernel
!> \param zmpi2 Work array for multiprocessor manipulation (output)
!> \param zw Work array (input)
!> \date February 2006
!> \author S. Goedecker, L. Genovese
!> \note
!>  RESTRICTIONS on USAGE
!>      Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!>      Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!>      Copyright (C) 2002 Stefan Goedecker, CEA Grenoble
!>      This file is distributed under the terms of the
!>      GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt . 
! *****************************************************************************
  SUBROUTINE unscramble_P(i1,j2,lot,nfft,n1,n3,md2,nproc,nd3,zmpi2,zw)
    INTEGER, INTENT(in)                      :: i1, j2, lot, nfft, n1, n3, &
                                                md2, nproc, nd3
    REAL(KIND=dp), &
      DIMENSION(2, n1, md2/nproc, nd3), &
      INTENT(in)                             :: zmpi2
    REAL(KIND=dp), DIMENSION(2, lot, n3), &
      INTENT(inout)                          :: zw

    INTEGER                                  :: i, i3, j3

    i3=1
    DO i=0,nfft-1
       zw(1,i+1,i3)=zmpi2(1,i1+i,j2,i3)
       zw(2,i+1,i3)=zmpi2(2,i1+i,j2,i3)
    END DO

    DO i3=2,n3/2+1
       j3=n3+2-i3
       DO i=0,nfft-1
          zw(1,i+1,j3)= zmpi2(1,i1+i,j2,i3)
          zw(2,i+1,j3)=-zmpi2(2,i1+i,j2,i3)
          zw(1,i+1,i3)= zmpi2(1,i1+i,j2,i3)
          zw(2,i+1,i3)= zmpi2(2,i1+i,j2,i3)
       END DO
    END DO

  END SUBROUTINE unscramble_P

! *****************************************************************************
!> \brief (Based on suitable modifications of S.Goedecker routines)
!>      Multiply with the kernel taking into account its symmetry
!>      Conceived to be used into convolution loops
!> \param n1 ...
!> \param n2 ...
!> \param n3 ...
!> \param lot ...
!> \param nfft ...
!> \param jS ...
!> \param i3 ...
!> \param zw Work array (input/output)
!>      n1,n2:      logical dimension of the FFT transform, reference for zw
!>      nd1,nd2:    Dimensions of POT
!>      jS,j3,nfft: starting point of the plane and number of remaining lines
!>
!>  RESTRICTIONS on USAGE
!>      Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!>      Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!>      Copyright (C) 2002 Stefan Goedecker, CEA Grenoble
!>      This file is distributed under the terms of the
!>       GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!> \param hx ...
!> \param hy ...
!> \param hz ...
!> \date February 2006
!> \author S. Goedecker, L. Genovese
! *****************************************************************************
  SUBROUTINE P_multkernel(n1,n2,n3,lot,nfft,jS,i3,zw,hx,hy,hz)
    !Argments
    INTEGER, INTENT(in)                      :: n1, n2, n3, lot, nfft, jS, i3
    REAL(KIND=dp), DIMENSION(2, lot, n2), &
      INTENT(inout)                          :: zw
    REAL(KIND=dp), INTENT(in)                :: hx, hy, hz

    INTEGER                                  :: i1, i2, j1, j2, j3
    REAL(KIND=dp)                            :: fourpi2, ker, mu3, p1, p2, pi

    pi=4._dp*ATAN(1._dp)
    fourpi2=4._dp*pi**2
    j3=i3!n3/2+1-abs(n3/2+2-i3)
    mu3=REAL(j3-1,KIND=dp)/REAL(n3,KIND=dp)
    mu3=(mu3/hy)**2 !beware of the exchanged dimension
    !Body
    !generic case
    DO i2=1,n2
       DO i1=1,nfft
          j1=i1+jS-1
          j1=j1-(j1/(n1/2+2))*n1 !n1/2+1-abs(n1/2+2-jS-i1)
          j2=i2-(i2/(n2/2+2))*n2 !n2/2+1-abs(n2/2+1-i2)
          p1=REAL(j1-1,KIND=dp)/REAL(n1,KIND=dp)
          p2=REAL(j2-1,KIND=dp)/REAL(n2,KIND=dp)
          ker=-fourpi2*((p1/hx)**2+(p2/hz)**2+mu3)!beware of the exchanged dimension
          IF (ker/=0._dp) ker=1._dp/ker
          zw(1,i1,i2)=zw(1,i1,i2)*ker
          zw(2,i1,i2)=zw(2,i1,i2)*ker
       END DO
    END DO

  END SUBROUTINE P_multkernel

! *****************************************************************************
!> \brief (Based on suitable modifications of S.Goedecker routines)
!>      Multiply with the kernel taking into account its symmetry
!>      Conceived to be used into convolution loops
!> \param nd1 ...
!> \param nd2 ...
!> \param n1 ...
!> \param n2 ...
!> \param lot ...
!> \param nfft ...
!> \param jS ...
!> \param pot Kernel, symmetric and real, half the length
!> \param zw Work array (input/output)
!>      n1,n2:    logical dimension of the FFT transform, reference for zw
!>      nd1,nd2:  Dimensions of POT
!>      jS, nfft: starting point of the plane and number of remaining lines
!>
!>  RESTRICTIONS on USAGE
!>      Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!>      Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!>      Copyright (C) 2002 Stefan Goedecker, CEA Grenoble
!>      This file is distributed under the terms of the
!>       GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!> \date February 2006
!> \author S. Goedecker, L. Genovese
! *****************************************************************************
  SUBROUTINE multkernel(nd1,nd2,n1,n2,lot,nfft,jS,pot,zw)
    !Argments
    INTEGER, INTENT(in)                      :: nd1, nd2, n1, n2, lot, nfft, &
                                                jS
    REAL(KIND=dp), DIMENSION(nd1, nd2), &
      INTENT(in)                             :: pot
    REAL(KIND=dp), DIMENSION(2, lot, n2), &
      INTENT(inout)                          :: zw

    INTEGER                                  :: i2, j, j1, j2

    DO j=1,nfft
       j1=j+jS-1
       j1=j1+(j1/(n1/2+2))*(n1+2-2*j1)
       zw(1,j,1)=zw(1,j,1)*pot(j1,1)
       zw(2,j,1)=zw(2,j,1)*pot(j1,1)
    END DO

    !generic case
    DO i2=2,n2/2
       DO j=1,nfft
          j1=j+jS-1
          j1=j1+(j1/(n1/2+2))*(n1+2-2*j1)
          j2=n2+2-i2
          zw(1,j,i2)=zw(1,j,i2)*pot(j1,i2)
          zw(2,j,i2)=zw(2,j,i2)*pot(j1,i2)
          zw(1,j,j2)=zw(1,j,j2)*pot(j1,i2)
          zw(2,j,j2)=zw(2,j,j2)*pot(j1,i2)
       END DO
    END DO

    !case i2=n2/2+1
    DO j=1,nfft
       j1=j+jS-1
       j1=j1+(j1/(n1/2+2))*(n1+2-2*j1)
       j2=n2/2+1
       zw(1,j,j2)=zw(1,j,j2)*pot(j1,j2)
       zw(2,j,j2)=zw(2,j,j2)*pot(j1,j2)
    END DO

  END SUBROUTINE multkernel

! *****************************************************************************
!> \brief !HERE POT MUST BE THE KERNEL (BEWARE THE HALF DIMENSION)
!> ****h* BigDFT/S_PoissonSolver
!>      (Based on suitable modifications of S.Goedecker routines)
!>      Applies the local FFT space Kernel to the density in Real space.
!>      Does NOT calculate the LDA exchange-correlation terms
!> \param n1 logical dimension of the transform.
!> \param n2 logical dimension of the transform.
!> \param n3 logical dimension of the transform.
!> \param nd1 Dimension of POT
!> \param nd2 Dimension of POT
!> \param nd3 Dimension of POT
!> \param md1 Dimension of ZF
!> \param md2 Dimension of ZF
!> \param md3 Dimension of ZF
!> \param nproc number of processors used as returned by MPI_COMM_SIZE
!> \param iproc [0:nproc-1] number of processor as returned by MPI_COMM_RANK
!> \param pot Kernel, only the distributed part (REAL)
!>                   POT(i1,i2,i3)
!>                   i1=1,nd1 , i2=1,nd2 , i3=1,nd3/nproc
!> \param zf Density (input/output)
!>                   ZF(i1,i3,i2)
!>                   i1=1,md1 , i2=1,md2/nproc , i3=1,md3
!> \param scal factor of renormalization of the FFT in order to acheve unitarity
!>                   and the correct dimension
!> \param hx grid spacing, used for integrating eharthree
!> \param hy grid spacing, used for integrating eharthree
!> \param hz grid spacing, used for integrating eharthree
!> \param mpi_group ...
!> \param error ...
!> \date October 2006
!> \author S. Goedecker, L. Genovese
!> \note   As transform lengths most products of the prime factors 2,3,5 are allowed.
!>                   The detailed table with allowed transform lengths can
!>                   be found in subroutine CTRIG
!> \note
!>  RESTRICTIONS on USAGE
!>      Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!>      Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!>      Copyright (C) 2002 Stefan Goedecker, CEA Grenoble
!>      This file is distributed under the terms of the
!>       GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt . 
! *****************************************************************************
  SUBROUTINE S_PoissonSolver(n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,nproc,iproc,pot,zf,&
       scal,hx,hy,hz,mpi_group,error)
    INTEGER, INTENT(in)                      :: n1, n2, n3, nd1, nd2, nd3, &
                                                md1, md2, md3, nproc, iproc
    REAL(KIND=dp), &
      DIMENSION(nd1, nd2, nd3/nproc), &
      INTENT(in)                             :: pot
    REAL(KIND=dp), &
      DIMENSION(md1, md3, md2/nproc), &
      INTENT(inout)                          :: zf
    REAL(KIND=dp), INTENT(in)                :: scal, hx, hy, hz
    INTEGER, INTENT(in)                      :: mpi_group
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'S_PoissonSolver', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: ncache_optimal = 8*1024, &
                                                timing_flag = 0

    INTEGER :: handle, i, i1, i3, i_all, i_stat, ic1, ic2, ic3, inzee, j, j2, &
      J2stb, J2stf, j3, Jp2stb, Jp2stf, lot, lzt, ma, mb, ncache, nfft
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: after1, after2, after3, &
                                                before1, before2, before3, &
                                                now1, now2, now3
    REAL(kind=dp)                            :: twopion
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: btrig1, btrig2, btrig3, &
                                                cosinarr, ftrig1, ftrig2, &
                                                ftrig3
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: zt, zw
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: zmpi2
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :, :)               :: zmpi1

    CALL timeset(routineN,handle)
    ! check input
    IF (MOD(n3,2).NE.0) STOP 'Parallel convolution:ERROR:n3'
    IF (nd1.LT.n1/2+1) STOP 'Parallel convolution:ERROR:nd1'
    IF (nd2.LT.n2/2+1) STOP 'Parallel convolution:ERROR:nd2'
    IF (nd3.LT.n3/2+1) STOP 'Parallel convolution:ERROR:nd3'
    IF (md1.LT.n1) STOP 'Parallel convolution:ERROR:md1'
    IF (md2.LT.n2) STOP 'Parallel convolution:ERROR:md2'
    IF (md3.LT.n3/2) STOP 'Parallel convolution:ERROR:md3'
    IF (MOD(nd3,nproc).NE.0) STOP 'Parallel convolution:ERROR:nd3'
    IF (MOD(md2,nproc).NE.0) STOP 'Parallel convolution:ERROR:md2'

    !defining work arrays dimensions
    ncache=ncache_optimal
    IF (ncache <= MAX(n1,n2,n3/2)*4) ncache=MAX(n1,n2,n3/2)*4

    !  if (timing_flag == 1 .and. iproc ==0) print *,'parallel ncache=',ncache

    lzt=n2
    IF (MOD(n2,2).EQ.0) lzt=lzt+1
    IF (MOD(n2,4).EQ.0) lzt=lzt+1 !maybe this is useless

    !Allocations
    ALLOCATE(btrig1(2,ctrig_length),stat=i_all)
    ALLOCATE(ftrig1(2,ctrig_length),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(after1(7),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(now1(7),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(before1(7),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(btrig2(2,ctrig_length),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(ftrig2(2,ctrig_length),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(after2(7),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(now2(7),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(before2(7),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(btrig3(2,ctrig_length),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(ftrig3(2,ctrig_length),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(after3(7),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(now3(7),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(before3(7),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(zw(2,ncache/4,2),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(zt(2,lzt,n1),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(zmpi2(2,n1,md2/nproc,nd3),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(cosinarr(2,n3/2),stat=i_stat)
    i_all=i_all+i_stat
    IF (nproc.GT.1)  THEN
       ALLOCATE(zmpi1(2,n1,md2/nproc,nd3/nproc,nproc),stat=i_stat)
       zmpi1= 0.0_dp
    END IF
    IF (i_all+i_stat /= 0) THEN
       WRITE(*,*)' s_poissonsolver: problem of memory allocation'
       STOP
    END IF
    zmpi2 = 0.0_dp
    zw = 0.0_dp
    zt = 0.0_dp

    !calculating the FFT work arrays (beware on the HalFFT in n3 dimension)
    CALL ctrig(n3/2,btrig3,after3,before3,now3,1,ic3)
    CALL ctrig(n1,btrig1,after1,before1,now1,1,ic1)
    CALL ctrig(n2,btrig2,after2,before2,now2,1,ic2)
    DO  j=1,n1
       ftrig1(1,j)= btrig1(1,j)
       ftrig1(2,j)=-btrig1(2,j)
    ENDDO
    DO  j=1,n2
       ftrig2(1,j)= btrig2(1,j)
       ftrig2(2,j)=-btrig2(2,j)
    ENDDO
    DO  j=1,n3/2
       ftrig3(1,j)= btrig3(1,j)
       ftrig3(2,j)=-btrig3(2,j)
    ENDDO

    !Calculating array of phases for HalFFT decoding
    twopion=8._dp*ATAN(1._dp)/REAL(n3,KIND=dp)
    DO i3=1,n3/2
       cosinarr(1,i3)=COS(twopion*(i3-1))
       cosinarr(2,i3)=-SIN(twopion*(i3-1))
    END DO

    !initializing integral
    !ehartree=0._dp

    ! transform along z axis
    lot=ncache/(2*n3)
    IF (lot.LT.1) THEN
       WRITE(*,*) &
            'convolxc_off:ncache has to be enlarged to be able to hold at' // &
            'least one 1-d FFT of this size even though this will' // &
            'reduce the performance for shorter transform lengths',n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,nproc,iproc
       STOP
    ENDIF

    DO j2=1,md2/nproc
       !this condition ensures that we manage only the interesting part for the FFT
       IF (iproc*(md2/nproc)+j2.LE.n2) THEN
          DO i1=1,n1,lot
             ma=i1
             mb=MIN(i1+(lot-1),n1)
             nfft=mb-ma+1

             !inserting real data into complex array of half lenght
             CALL halfill_upcorn(md1,md3,lot,nfft,n3,zf(i1,1,j2),zw(1,1,1))

             !performing FFT
             !input: I1,I3,J2,(Jp2)
             inzee=1
             DO i=1,ic3
                CALL fftstp(lot,nfft,n3/2,lot,n3/2,zw(1,1,inzee),zw(1,1,3-inzee), &
                     btrig3,after3(i),now3(i),before3(i),1)
                inzee=3-inzee
             ENDDO
             !output: I1,i3,J2,(Jp2)
             !unpacking FFT in order to restore correct result,
             !while exchanging components
             !input: I1,i3,J2,(Jp2)
             CALL scramble_unpack(i1,j2,lot,nfft,n1,n3,md2,nproc,nd3,zw(1,1,inzee),zmpi2,cosinarr)
             !output: I1,J2,i3,(Jp2)
          END DO
       END IF
    END DO
    !Interprocessor data transposition
    !input: I1,J2,j3,jp3,(Jp2)
    IF (nproc.GT.1) THEN
       CALL mp_alltoall(zmpi2,zmpi1,2*n1*(md2/nproc)*(nd3/nproc),mpi_group)
    ENDIF
    !output: I1,J2,j3,Jp2,(jp3)

    !now each process perform complete convolution of its planes
    DO j3=1,nd3/nproc
       !this condition ensures that we manage only the interesting part for the FFT
       IF (iproc*(nd3/nproc)+j3.LE.n3/2+1) THEN
          Jp2stb=1
          J2stb=1
          Jp2stf=1
          J2stf=1

          ! transform along x axis
          lot=ncache/(4*n1)
          IF (lot.LT.1) THEN
             WRITE(*,*) &
                  'convolxc_off:ncache has to be enlarged to be able to hold at' // &
                  'least one 1-d FFT of this size even though this will' // &
                  'reduce the performance for shorter transform lengths'
             STOP
          ENDIF

          DO j=1,n2,lot
             ma=j
             mb=MIN(j+(lot-1),n2)
             nfft=mb-ma+1

             !reverse index ordering, leaving the planes to be transformed at the end
             !input: I1,J2,j3,Jp2,(jp3)
             IF (nproc.EQ.1) THEN
                CALL S_mpiswitch_upcorn(j3,nfft,Jp2stb,J2stb,lot,n1,md2,nd3,nproc,zmpi2,zw(1,1,1))
             ELSE
                CALL S_mpiswitch_upcorn(j3,nfft,Jp2stb,J2stb,lot,n1,md2,nd3,nproc,zmpi1,zw(1,1,1))
             ENDIF
             !output: J2,Jp2,I1,j3,(jp3)

             !performing FFT
             !input: I2,I1,j3,(jp3)
             inzee=1
             DO i=1,ic1-1
                CALL fftstp(lot,nfft,n1,lot,n1,zw(1,1,inzee),zw(1,1,3-inzee), &
                     btrig1,after1(i),now1(i),before1(i),1)
                inzee=3-inzee
             ENDDO

             !storing the last step into zt array
             i=ic1
             CALL fftstp(lot,nfft,n1,lzt,n1,zw(1,1,inzee),zt(1,j,1), &
                  btrig1,after1(i),now1(i),before1(i),1)
             !output: I2,i1,j3,(jp3)
          END DO

          !transform along y axis
          lot=ncache/(4*n2)
          IF (lot.LT.1) THEN
             WRITE(*,*) &
                  'convolxc_off:ncache has to be enlarged to be able to hold at' // &
                  'least one 1-d FFT of this size even though this will' // &
                  'reduce the performance for shorter transform lengths'
             STOP
          ENDIF

          DO j=1,n1,lot
             ma=j
             mb=MIN(j+(lot-1),n1)
             nfft=mb-ma+1

             !reverse ordering
             !input: I2,i1,j3,(jp3)
             CALL S_switch_upcorn(nfft,n2,lot,n1,lzt,zt(1,1,j),zw(1,1,1))
             !output: i1,I2,j3,(jp3)

             !performing FFT
             !input: i1,I2,j3,(jp3)
             inzee=1
             DO i=1,ic2
                CALL fftstp(lot,nfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee), &
                     btrig2,after2(i),now2(i),before2(i),1)
                inzee=3-inzee
             ENDDO
             !output: i1,i2,j3,(jp3)

             !Multiply with kernel in fourier space
             CALL multkernel(nd1,nd2,n1,n2,lot,nfft,j,pot(1,1,j3),zw(1,1,inzee))

             !TRANSFORM BACK IN REAL SPACE

             !transform along y axis
             !input: i1,i2,j3,(jp3)
             DO i=1,ic2
                CALL fftstp(lot,nfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee), &
                     ftrig2,after2(i),now2(i),before2(i),-1)
                inzee=3-inzee
             ENDDO

             !reverse ordering
             !input: i1,I2,j3,(jp3)
             CALL S_unswitch_downcorn(nfft,n2,lot,n1,lzt,zw(1,1,inzee),zt(1,1,j))
             !output: I2,i1,j3,(jp3)
          END DO

          !transform along x axis
          !input: I2,i1,j3,(jp3)
          lot=ncache/(4*n1)
          DO j=1,n2,lot
             ma=j
             mb=MIN(j+(lot-1),n2)
             nfft=mb-ma+1

             !performing FFT
             i=1
             CALL fftstp(lzt,nfft,n1,lot,n1,zt(1,j,1),zw(1,1,1), &
                  ftrig1,after1(i),now1(i),before1(i),-1)

             inzee=1
             DO i=2,ic1
                CALL fftstp(lot,nfft,n1,lot,n1,zw(1,1,inzee),zw(1,1,3-inzee), &
                     ftrig1,after1(i),now1(i),before1(i),-1)
                inzee=3-inzee
             ENDDO
             !output: I2,I1,j3,(jp3)

             !reverse ordering
             !input: J2,Jp2,I1,j3,(jp3)
             IF (nproc.EQ.1) THEN
                CALL S_unmpiswitch_downcorn(j3,nfft,Jp2stf,J2stf,lot,n1,md2,nd3,nproc,zw(1,1,inzee),zmpi2)
             ELSE
                CALL S_unmpiswitch_downcorn(j3,nfft,Jp2stf,J2stf,lot,n1,md2,nd3,nproc,zw(1,1,inzee),zmpi1)
             ENDIF
             ! output: I1,J2,j3,Jp2,(jp3)
          END DO
       ENDIF
    END DO

    !Interprocessor data transposition
    !input: I1,J2,j3,Jp2,(jp3)
    IF (nproc.GT.1) THEN
       !communication scheduling
       CALL mp_alltoall(zmpi1,zmpi2,2*n1*(md2/nproc)*(nd3/nproc),mpi_group)
    ENDIF

    !output: I1,J2,j3,jp3,(Jp2)

    !transform along z axis
    !input: I1,J2,i3,(Jp2)
    lot=ncache/(2*n3)
    DO j2=1,md2/nproc
       !this condition ensures that we manage only the interesting part for the FFT
       IF (iproc*(md2/nproc)+j2.LE.n2) THEN
          DO i1=1,n1,lot
             ma=i1
             mb=MIN(i1+(lot-1),n1)
             nfft=mb-ma+1

             !reverse ordering and repack the FFT data in order to be backward HalFFT transformed
             !input: I1,J2,i3,(Jp2)
             CALL unscramble_pack(i1,j2,lot,nfft,n1,n3,md2,nproc,nd3,zmpi2,zw(1,1,1),cosinarr)
             !output: I1,i3,J2,(Jp2)

             !performing FFT
             !input: I1,i3,J2,(Jp2)
             inzee=1
             DO i=1,ic3
                CALL fftstp(lot,nfft,n3/2,lot,n3/2,zw(1,1,inzee),zw(1,1,3-inzee), &
                     ftrig3,after3(i),now3(i),before3(i),-1)
                inzee=3-inzee
             ENDDO
             !output: I1,I3,J2,(Jp2)

             !rebuild the output array
             CALL unfill_downcorn(md1,md3,lot,nfft,n3,zw(1,1,inzee),zf(i1,1,j2)&
                  ,scal)

             !integrate local pieces together
             !ehartree=ehartree+0.5_dp*ehartreetmp*hx*hy*hz
          END DO
       ENDIF
    END DO

    !De-allocations
    DEALLOCATE(btrig1,stat=i_all)
    DEALLOCATE(ftrig1,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(after1,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(now1,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(before1,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(btrig2,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(ftrig2,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(after2,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(now2,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(before2,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(btrig3,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(ftrig3,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(after3,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(now3,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(before3,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(zmpi2,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(zw,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(zt,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(cosinarr,stat=i_stat)
    i_all=i_all+i_stat
    IF (nproc.GT.1) DEALLOCATE(zmpi1,stat=i_stat)
    IF (i_all+i_stat /= 0) THEN
       WRITE(*,*)' s_poissonsolver: problem of memory deallocation'
       STOP
    END IF

    !  call timing(iproc,'PSolv_comput  ','OF')
    CALL timestop(handle)
  END SUBROUTINE S_PoissonSolver

! *****************************************************************************
!> \brief ...
!> \param j3 ...
!> \param nfft ...
!> \param Jp2stb ...
!> \param J2stb ...
!> \param lot ...
!> \param n1 ...
!> \param md2 ...
!> \param nd3 ...
!> \param nproc ...
!> \param zmpi1 ...
!> \param zw ...
! *****************************************************************************
  SUBROUTINE S_mpiswitch_upcorn(j3,nfft,Jp2stb,J2stb,lot,n1,md2,nd3,nproc,zmpi1,zw)
    INTEGER, INTENT(in)                      :: j3, nfft
    INTEGER, INTENT(inout)                   :: Jp2stb, J2stb
    INTEGER, INTENT(in)                      :: lot, n1, md2, nd3, nproc
    REAL(KIND=dp), DIMENSION(2, n1, md2/&
      nproc, nd3/nproc, nproc), INTENT(in)   :: zmpi1
    REAL(KIND=dp), DIMENSION(2, lot, n1), &
      INTENT(inout)                          :: zw

    INTEGER                                  :: I1, J2, Jp2, mfft

    mfft=0
    DO Jp2=Jp2stb,nproc
       DO J2=J2stb,md2/nproc
          mfft=mfft+1
          IF (mfft.GT.nfft) THEN
             Jp2stb=Jp2
             J2stb=J2
             RETURN
          END IF
          DO I1=1,n1
             zw(1,mfft,I1)=zmpi1(1,I1,J2,j3,Jp2)
             zw(2,mfft,I1)=zmpi1(2,I1,J2,j3,Jp2)
          END DO
       END DO
       J2stb=1
    END DO
  END SUBROUTINE S_mpiswitch_upcorn

! *****************************************************************************
!> \brief ...
!> \param nfft ...
!> \param n2 ...
!> \param lot ...
!> \param n1 ...
!> \param lzt ...
!> \param zt ...
!> \param zw ...
! *****************************************************************************
  SUBROUTINE S_switch_upcorn(nfft,n2,lot,n1,lzt,zt,zw)
    INTEGER, INTENT(in)                      :: nfft, n2, lot, n1, lzt
    REAL(KIND=dp), DIMENSION(2, lzt, n1), &
      INTENT(in)                             :: zt
    REAL(KIND=dp), DIMENSION(2, lot, n2), &
      INTENT(inout)                          :: zw

    INTEGER                                  :: i, j

    DO j=1,nfft
       DO i=1,n2
          zw(1,j,i)=zt(1,i,j)
          zw(2,j,i)=zt(2,i,j)
       END DO
    END DO
  END SUBROUTINE S_switch_upcorn

! *****************************************************************************
!> \brief ...
!> \param nfft ...
!> \param n2 ...
!> \param lot ...
!> \param n1 ...
!> \param lzt ...
!> \param zw ...
!> \param zt ...
! *****************************************************************************
  SUBROUTINE S_unswitch_downcorn(nfft,n2,lot,n1,lzt,zw,zt)
    INTEGER, INTENT(in)                      :: nfft, n2, lot, n1, lzt
    REAL(KIND=dp), DIMENSION(2, lot, n2), &
      INTENT(in)                             :: zw
    REAL(KIND=dp), DIMENSION(2, lzt, n1), &
      INTENT(inout)                          :: zt

    INTEGER                                  :: i, j

    DO j=1,nfft
       DO i=1,n2
          zt(1,i,j)=zw(1,j,i)
          zt(2,i,j)=zw(2,j,i)
       END DO
    END DO
  END SUBROUTINE S_unswitch_downcorn

! *****************************************************************************
!> \brief ...
!> \param j3 ...
!> \param nfft ...
!> \param Jp2stf ...
!> \param J2stf ...
!> \param lot ...
!> \param n1 ...
!> \param md2 ...
!> \param nd3 ...
!> \param nproc ...
!> \param zw ...
!> \param zmpi1 ...
! *****************************************************************************
  SUBROUTINE S_unmpiswitch_downcorn(j3,nfft,Jp2stf,J2stf,lot,n1,md2,nd3,nproc,zw,zmpi1)
    INTEGER, INTENT(in)                      :: j3, nfft
    INTEGER, INTENT(inout)                   :: Jp2stf, J2stf
    INTEGER, INTENT(in)                      :: lot, n1, md2, nd3, nproc
    REAL(KIND=dp), DIMENSION(2, lot, n1), &
      INTENT(in)                             :: zw
    REAL(KIND=dp), DIMENSION(2, n1, md2/&
      nproc, nd3/nproc, nproc), &
      INTENT(inout)                          :: zmpi1

    INTEGER                                  :: I1, J2, Jp2, mfft

    mfft=0
    DO Jp2=Jp2stf,nproc
       DO J2=J2stf,md2/nproc
          mfft=mfft+1
          IF (mfft.GT.nfft) THEN
             Jp2stf=Jp2
             J2stf=J2
             RETURN
          END IF
          DO I1=1,n1
             zmpi1(1,I1,J2,j3,Jp2)=zw(1,mfft,I1)
             zmpi1(2,I1,J2,j3,Jp2)=zw(2,mfft,I1)
          END DO
       END DO
       J2stf=1
    END DO
  END SUBROUTINE S_unmpiswitch_downcorn

! *****************************************************************************
!> \brief (Based on suitable modifications of S.Goedecker routines)
!>      Restore data into output array, calculating in the meanwhile
!>      Hartree energy of the potential
!> \param md1 Dimensions of the undistributed part of the real grid
!> \param md3 Dimensions of the undistributed part of the real grid
!> \param lot ...
!> \param nfft number of planes
!> \param n3 (twice the) dimension of the last FFTtransform.
!> \param zw FFT work array
!> \param zf Original distributed density as well as
!>                   Distributed solution of the poisson equation (inout)
!> \param scal Needed to achieve unitarity and correct dimensions
!> \date February 2006
!> \author S. Goedecker, L. Genovese
!> \note Assuming that high frequencies are in the corners 
!>      and that n3 is multiple of 4 
!> 
!>  RESTRICTIONS on USAGE 
!>      Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994 
!>      Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999 
!>      Copyright (C) 2002 Stefan Goedecker, CEA Grenoble 
!>      This file is distributed under the terms of the 
!>      GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt . 
! *****************************************************************************
  SUBROUTINE unfill_downcorn(md1,md3,lot,nfft,n3,zw,zf,scal)
    INTEGER, INTENT(in)                      :: md1, md3, lot, nfft, n3
    REAL(KIND=dp), DIMENSION(2, lot, n3/2), &
      INTENT(in)                             :: zw
    REAL(KIND=dp), DIMENSION(md1, md3), &
      INTENT(inout)                          :: zf
    REAL(KIND=dp), INTENT(in)                :: scal

    INTEGER                                  :: i1, i3
    REAL(KIND=dp)                            :: pot1

    DO i3=1,n3/4
       DO i1=1,nfft
          pot1 = scal*zw(1,i1,i3)
          !ehartreetmp =ehartreetmp + pot1* zf(i1,2*i3-1)
          zf(i1,2*i3-1)= pot1
          pot1 = scal*zw(2,i1,i3)
          !ehartreetmp =ehartreetmp + pot1* zf(i1,2*i3)
          zf(i1,2*i3)= pot1
       ENDDO
    END DO
  END SUBROUTINE unfill_downcorn

! *****************************************************************************
!> \brief ...
!> \param md1 ...
!> \param md3 ...
!> \param lot ...
!> \param nfft ...
!> \param n3 ...
!> \param zf ...
!> \param zw ...
! *****************************************************************************
  SUBROUTINE halfill_upcorn(md1,md3,lot,nfft,n3,zf,zw)
    INTEGER                                  :: md1, md3, lot, nfft, n3
    REAL(KIND=dp)                            :: zf(md1,md3), zw(2,lot,n3/2)

    INTEGER                                  :: i1, i3

    DO i3=1,n3/4
       ! WARNING: Assuming that high frequencies are in the corners
       !          and that n3 is multiple of 4
       !in principle we can relax this condition
       DO i1=1,nfft
          zw(1,i1,i3)=0._dp
          zw(2,i1,i3)=0._dp
       END DO
    END DO
    DO i3=n3/4+1,n3/2
       DO i1=1,nfft
          zw(1,i1,i3)=zf(i1,2*i3-1-n3/2)
          zw(2,i1,i3)=zf(i1,2*i3-n3/2)
       ENDDO
    END DO

  END SUBROUTINE halfill_upcorn

! *****************************************************************************
!> \brief (Based on suitable modifications of S.Goedecker routines)
!>      Assign the correct planes to the work array zmpi2
!>      in order to prepare for interprocessor data transposition.
!>      In the meanwhile, it unpacks the data of the HalFFT in order to prepare for
!>      multiplication with the kernel
!> \param i1 Starting points of the plane and number of remaining lines
!> \param j2 Starting points of the plane and number of remaining lines
!> \param lot Starting points of the plane and number of remaining lines
!> \param nfft Starting points of the plane and number of remaining lines
!> \param n1 logical dimension of the FFT transform, reference for work arrays
!> \param n3 logical dimension of the FFT transform, reference for work arrays
!> \param md2 Dimensions of real grid
!> \param nproc ...
!> \param nd3 Dimensions of the kernel
!> \param zw Work array (input)
!> \param zmpi2 Work array for multiprocessor manipulation (output)
!> \param cosinarr Array of the phases needed for unpacking
!> \date February 2006
!> \author S. Goedecker, L. Genovese
!> \note
!>  RESTRICTIONS on USAGE
!>      Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!>      Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!>      Copyright (C) 2002 Stefan Goedecker, CEA Grenoble
!>      This file is distributed under the terms of the
!>      GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt . 
! *****************************************************************************
  SUBROUTINE scramble_unpack(i1,j2,lot,nfft,n1,n3,md2,nproc,nd3,zw,zmpi2,cosinarr)
    INTEGER, INTENT(in)                      :: i1, j2, lot, nfft, n1, n3, &
                                                md2, nproc, nd3
    REAL(KIND=dp), DIMENSION(2, lot, n3/2), &
      INTENT(in)                             :: zw
    REAL(KIND=dp), &
      DIMENSION(2, n1, md2/nproc, nd3), &
      INTENT(inout)                          :: zmpi2
    REAL(KIND=dp), DIMENSION(2, n3/2), &
      INTENT(in)                             :: cosinarr

    INTEGER                                  :: i, i3, ind1, ind2
    REAL(KIND=dp)                            :: a, b, c, cp, d, feI, feR, fI, &
                                                foI, foR, fR, sp

!case i3=1 and i3=n3/2+1

    DO i=0,nfft-1
       a=zw(1,i+1,1)
       b=zw(2,i+1,1)
       zmpi2(1,i1+i,j2,1)=a+b
       zmpi2(2,i1+i,j2,1)=0._dp
       zmpi2(1,i1+i,j2,n3/2+1)=a-b
       zmpi2(2,i1+i,j2,n3/2+1)=0._dp
    END DO
    !case 2<=i3<=n3/2
    DO i3=2,n3/2
       ind1=i3
       ind2=n3/2-i3+2
       cp=cosinarr(1,i3)
       sp=cosinarr(2,i3)
       DO i=0,nfft-1
          a=zw(1,i+1,ind1)
          b=zw(2,i+1,ind1)
          c=zw(1,i+1,ind2)
          d=zw(2,i+1,ind2)
          feR=.5_dp*(a+c)
          feI=.5_dp*(b-d)
          foR=.5_dp*(a-c)
          foI=.5_dp*(b+d)
          fR=feR+cp*foI-sp*foR
          fI=feI-cp*foR-sp*foI
          zmpi2(1,i1+i,j2,ind1)=fR
          zmpi2(2,i1+i,j2,ind1)=fI
       END DO
    END DO

  END SUBROUTINE scramble_unpack

! *****************************************************************************
!> \brief (Based on suitable modifications of S.Goedecker routines)
!>      Insert the correct planes of the work array zmpi2
!>      in order to prepare for backward FFT transform
!>      In the meanwhile, it packs the data in order to be transformed with the HalFFT
!>      procedure
!> \param i1 Starting points of the plane and number of remaining lines
!> \param j2 Starting points of the plane and number of remaining lines
!> \param lot Starting points of the plane and number of remaining lines
!> \param nfft Starting points of the plane and number of remaining lines
!> \param n1 logical dimension of the FFT transform, reference for work arrays
!> \param n3 logical dimension of the FFT transform, reference for work arrays
!> \param md2 Dimensions of real grid
!> \param nproc ...
!> \param nd3 Dimensions of the kernel
!> \param zmpi2 Work array for multiprocessor manipulation (output)
!> \param zw Work array (inout)
!> \param cosinarr Array of the phases needed for packing
!> \date February 2006
!> \author S. Goedecker, L. Genovese
!> \note
!>  RESTRICTIONS on USAGE
!>      Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!>      Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!>      Copyright (C) 2002 Stefan Goedecker, CEA Grenoble
!>      This file is distributed under the terms of the
!>      GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt . 
! *****************************************************************************
  SUBROUTINE unscramble_pack(i1,j2,lot,nfft,n1,n3,md2,nproc,nd3,zmpi2,zw,cosinarr)
    INTEGER, INTENT(in)                      :: i1, j2, lot, nfft, n1, n3, &
                                                md2, nproc, nd3
    REAL(KIND=dp), &
      DIMENSION(2, n1, md2/nproc, nd3), &
      INTENT(in)                             :: zmpi2
    REAL(KIND=dp), DIMENSION(2, lot, n3/2), &
      INTENT(inout)                          :: zw
    REAL(KIND=dp), DIMENSION(2, n3/2), &
      INTENT(in)                             :: cosinarr

    INTEGER                                  :: i, i3, indA, indB
    REAL(KIND=dp)                            :: a, b, c, cp, d, ie, ih, io, &
                                                re, rh, ro, sp

    DO i3=1,n3/2
       indA=i3
       indB=n3/2+2-i3
       cp=cosinarr(1,i3)
       sp=cosinarr(2,i3)
       DO i=0,nfft-1
          a=zmpi2(1,i1+i,j2,indA)
          b=zmpi2(2,i1+i,j2,indA)
          c= zmpi2(1,i1+i,j2,indB)
          d=-zmpi2(2,i1+i,j2,indB)
          re=(a+c)
          ie=(b+d)
          ro=(a-c)*cp-(b-d)*sp
          io=(a-c)*sp+(b-d)*cp
          rh=re-io
          ih=ie+ro
          zw(1,i+1,indA)=rh
          zw(2,i+1,indA)=ih
       END DO
    END DO

  END SUBROUTINE unscramble_pack

! *****************************************************************************
!> \brief (Based on suitable modifications of S.Goedecker routines)
!>      Applies the local FFT space Kernel to the density in Real space.
!>      Calculates also the LDA exchange-correlation terms
!> \param n1 logical dimension of the transform.
!> \param n2 logical dimension of the transform.
!> \param n3 logical dimension of the transform.
!> \param nd1 Dimension of POT
!> \param nd2 Dimension of POT
!> \param nd3 Dimension of POT
!> \param md1 Dimension of ZF
!> \param md2 Dimension of ZF
!> \param md3 Dimension of ZF
!> \param nproc number of processors used as returned by MPI_COMM_SIZE
!> \param iproc [0:nproc-1] number of processor as returned by MPI_COMM_RANK
!> \param pot Kernel, only the distributed part (REAL)
!>                   POT(i1,i2,i3)
!>                   i1=1,nd1 , i2=1,nd2 , i3=1,nd3/nproc
!> \param zf Density (input/output)
!>                   ZF(i1,i3,i2)
!>                   i1=1,md1 , i2=1,md2/nproc , i3=1,md3
!> \param scal factor of renormalization of the FFT in order to acheve unitarity
!>                   and the correct dimension 
!> \param hgrid grid spacing, used for integrating quantities 
!> \param mpi_group ...
!> \param error ...
!> \param ! ...
!> \param ehartree hartree energy of the potential
!> \date February 2006
!> \author S. Goedecker, L. Genovese
!> \note As transform lengths
!>                     most products of the prime factors 2,3,5 are allowed.
!>                   The detailed table with allowed transform lengths can
!>                   be found in subroutine CTRIG
!> \note
!>  RESTRICTIONS on USAGE
!>      Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!>      Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!>      Copyright (C) 2002 Stefan Goedecker, CEA Grenoble
!>      This file is distributed under the terms of the
!>       GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt . 
! *****************************************************************************
  SUBROUTINE F_PoissonSolver(n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,nproc,iproc,pot,zf,&
       scal,hgrid,mpi_group,error)!,ehartree)
    INTEGER, INTENT(in)                      :: n1, n2, n3, nd1, nd2, nd3, &
                                                md1, md2, md3, nproc, iproc
    REAL(KIND=dp), &
      DIMENSION(nd1, nd2, nd3/nproc), &
      INTENT(in)                             :: pot
    REAL(KIND=dp), &
      DIMENSION(md1, md3, md2/nproc), &
      INTENT(inout)                          :: zf
    REAL(KIND=dp), INTENT(in)                :: scal, hgrid
    INTEGER, INTENT(in)                      :: mpi_group
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    INTEGER, PARAMETER                       :: ncache_optimal = 8*1024, &
                                                timing_flag = 0

    INTEGER :: i, i1, i3, i_all, i_stat, ic1, ic2, ic3, inzee, j, j2, J2stb, &
      J2stf, j3, Jp2stb, Jp2stf, lot, lzt, ma, mb, ncache, nfft
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: after1, after2, after3, &
                                                before1, before2, before3, &
                                                now1, now2, now3
    REAL(kind=dp)                            :: twopion
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: btrig1, btrig2, btrig3, &
                                                cosinarr, ftrig1, ftrig2, &
                                                ftrig3
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: zt, zw
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: zmpi2
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :, :)               :: zmpi1

    IF (MOD(n1,2).NE.0) STOP 'Parallel convolution:ERROR:n1'
    IF (MOD(n2,2).NE.0) STOP 'Parallel convolution:ERROR:n2'
    IF (MOD(n3,2).NE.0) STOP 'Parallel convolution:ERROR:n3'
    IF (nd1.LT.n1/2+1) STOP 'Parallel convolution:ERROR:nd1'
    IF (nd2.LT.n2/2+1) STOP 'Parallel convolution:ERROR:nd2'
    IF (nd3.LT.n3/2+1) STOP 'Parallel convolution:ERROR:nd3'
    IF (md1.LT.n1/2) STOP 'Parallel convolution:ERROR:md1'
    IF (md2.LT.n2/2) STOP 'Parallel convolution:ERROR:md2'
    IF (md3.LT.n3/2) STOP 'Parallel convolution:ERROR:md3'
    IF (MOD(nd3,nproc).NE.0) STOP 'Parallel convolution:ERROR:nd3'
    IF (MOD(md2,nproc).NE.0) STOP 'Parallel convolution:ERROR:md2'

    !defining work arrays dimensions

    ncache=ncache_optimal
    IF (ncache <= MAX(n1,n2,n3/2)*4) ncache=MAX(n1,n2,n3/2)*4
    lzt=n2/2
    IF (MOD(n2/2,2).EQ.0) lzt=lzt+1
    IF (MOD(n2/2,4).EQ.0) lzt=lzt+1

    !Allocations
    ALLOCATE(btrig1(2,ctrig_length),stat=i_all)
    ALLOCATE(ftrig1(2,ctrig_length),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(after1(7),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(now1(7),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(before1(7),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(btrig2(2,ctrig_length),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(ftrig2(2,ctrig_length),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(after2(7),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(now2(7),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(before2(7),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(btrig3(2,ctrig_length),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(ftrig3(2,ctrig_length),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(after3(7),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(now3(7),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(before3(7),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(zw(2,ncache/4,2),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(zt(2,lzt,n1),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(zmpi2(2,n1,md2/nproc,nd3),stat=i_stat)
    i_all=i_all+i_stat
    ALLOCATE(cosinarr(2,n3/2),stat=i_stat)
    i_all=i_all+i_stat
    IF (nproc.GT.1) ALLOCATE(zmpi1(2,n1,md2/nproc,nd3/nproc,nproc),stat=i_stat)
    IF (i_all+i_stat /= 0) THEN
       WRITE(*,*)' f_poissonsolver: problem of memory allocation'
       STOP
    END IF

    !calculating the FFT work arrays (beware on the HalFFT in n3 dimension)
    CALL ctrig(n3/2,btrig3,after3,before3,now3,1,ic3)
    CALL ctrig(n1,btrig1,after1,before1,now1,1,ic1)
    CALL ctrig(n2,btrig2,after2,before2,now2,1,ic2)
    DO  j=1,n1
       ftrig1(1,j)= btrig1(1,j)
       ftrig1(2,j)=-btrig1(2,j)
    ENDDO
    DO  j=1,n2
       ftrig2(1,j)= btrig2(1,j)
       ftrig2(2,j)=-btrig2(2,j)
    ENDDO
    DO  j=1,n3
       ftrig3(1,j)= btrig3(1,j)
       ftrig3(2,j)=-btrig3(2,j)
    ENDDO

    !Calculating array of phases for HalFFT decoding
    twopion=8._dp*ATAN(1._dp)/REAL(n3,KIND=dp)
    DO i3=1,n3/2
       cosinarr(1,i3)=COS(twopion*(i3-1))
       cosinarr(2,i3)=-SIN(twopion*(i3-1))
    END DO

    ! transform along z axis
    lot=ncache/(2*n3)
    IF (lot.LT.1) THEN
       WRITE(*,*) &
            'convolxc_on:ncache has to be enlarged to be able to hold at' // &
            'least one 1-d FFT of this size even though this will' // &
            'reduce the performance for shorter transform lengths'
       STOP
    ENDIF

    DO j2=1,md2/nproc
       !this condition ensures that we manage only the interesting part for the FFT
       IF (iproc*(md2/nproc)+j2.LE.n2/2) THEN
          DO i1=1,(n1/2),lot
             ma=i1
             mb=MIN(i1+(lot-1),(n1/2))
             nfft=mb-ma+1

             !inserting real data into complex array of half lenght
             CALL halfill_upcorn(md1,md3,lot,nfft,n3,zf(i1,1,j2),zw(1,1,1))

             !performing FFT
             !input: I1,I3,J2,(Jp2)
             inzee=1
             DO i=1,ic3
                CALL fftstp(lot,nfft,n3/2,lot,n3/2,zw(1,1,inzee),zw(1,1,3-inzee), &
                     btrig3,after3(i),now3(i),before3(i),1)
                inzee=3-inzee
             ENDDO
             !output: I1,i3,J2,(Jp2)

             !unpacking FFT in order to restore correct result,
             !while exchanging components
             !input: I1,i3,J2,(Jp2)
             CALL scramble_unpack(i1,j2,lot,nfft,n1/2,n3,md2,nproc,nd3,zw(1,1,inzee),zmpi2,cosinarr)
             !output: I1,J2,i3,(Jp2)
          END DO
       ENDIF
    END DO

    !Interprocessor data transposition
    !input: I1,J2,j3,jp3,(Jp2)
    IF (nproc.GT.1) THEN
       !communication scheduling
       CALL mp_alltoall(zmpi2,zmpi1,n1*(md2/nproc)*(nd3/nproc),mpi_group)
    ENDIF
    !output: I1,J2,j3,Jp2,(jp3)

    !now each process perform complete convolution of its planes
    DO j3=1,nd3/nproc
       !this condition ensures that we manage only the interesting part for the FFT
       IF (iproc*(nd3/nproc)+j3.LE.n3/2+1) THEN
          Jp2stb=1
          J2stb=1
          Jp2stf=1
          J2stf=1

          ! transform along x axis
          lot=ncache/(4*n1)
          IF (lot.LT.1) THEN
             WRITE(*,*) &
                  'convolxc_on:ncache has to be enlarged to be able to hold at' // &
                  'least one 1-d FFT of this size even though this will' // &
                  'reduce the performance for shorter transform lengths'
             STOP
          ENDIF

          DO j=1,n2/2,lot
             ma=j
             mb=MIN(j+(lot-1),n2/2)
             nfft=mb-ma+1

             !reverse index ordering, leaving the planes to be transformed at the end
             !input: I1,J2,j3,Jp2,(jp3)
             IF (nproc.EQ.1) THEN
                CALL mpiswitch_upcorn(j3,nfft,Jp2stb,J2stb,lot,n1,md2,nd3,nproc,zmpi2,zw(1,1,1))
             ELSE
                CALL mpiswitch_upcorn(j3,nfft,Jp2stb,J2stb,lot,n1,md2,nd3,nproc,zmpi1,zw(1,1,1))
             ENDIF
             !output: J2,Jp2,I1,j3,(jp3)

             !performing FFT
             !input: I2,I1,j3,(jp3)
             inzee=1
             DO i=1,ic1-1
                CALL fftstp(lot,nfft,n1,lot,n1,zw(1,1,inzee),zw(1,1,3-inzee), &
                     btrig1,after1(i),now1(i),before1(i),1)
                inzee=3-inzee
             ENDDO

             !storing the last step into zt array
             i=ic1
             CALL fftstp(lot,nfft,n1,lzt,n1,zw(1,1,inzee),zt(1,j,1), &
                  btrig1,after1(i),now1(i),before1(i),1)
             !output: I2,i1,j3,(jp3)
          END DO

          !transform along y axis
          lot=ncache/(4*n2)
          IF (lot.LT.1) THEN
             WRITE(*,*) &
                  'convolxc_on:ncache has to be enlarged to be able to hold at' // &
                  'least one 1-d FFT of this size even though this will' // &
                  'reduce the performance for shorter transform lengths'
             STOP
          ENDIF

          DO j=1,n1,lot
             ma=j
             mb=MIN(j+(lot-1),n1)
             nfft=mb-ma+1

             !reverse ordering
             !input: I2,i1,j3,(jp3)
             CALL switch_upcorn(nfft,n2,lot,n1,lzt,zt(1,1,j),zw(1,1,1))
             !output: i1,I2,j3,(jp3)

             !performing FFT
             !input: i1,I2,j3,(jp3)
             inzee=1
             DO i=1,ic2
                CALL fftstp(lot,nfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee), &
                     btrig2,after2(i),now2(i),before2(i),1)
                inzee=3-inzee
             ENDDO
             !output: i1,i2,j3,(jp3)

             !Multiply with kernel in fourier space
             CALL multkernel(nd1,nd2,n1,n2,lot,nfft,j,pot(1,1,j3),zw(1,1,inzee))

             !TRANSFORM BACK IN REAL SPACE

             !transform along y axis
             !input: i1,i2,j3,(jp3)
             DO i=1,ic2
                CALL fftstp(lot,nfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee), &
                     ftrig2,after2(i),now2(i),before2(i),-1)
                inzee=3-inzee
             ENDDO

             !reverse ordering
             !input: i1,I2,j3,(jp3)
             CALL unswitch_downcorn(nfft,n2,lot,n1,lzt,zw(1,1,inzee),zt(1,1,j))
             !output: I2,i1,j3,(jp3)
          END DO

          !transform along x axis
          !input: I2,i1,j3,(jp3)
          lot=ncache/(4*n1)
          DO j=1,n2/2,lot
             ma=j
             mb=MIN(j+(lot-1),n2/2)
             nfft=mb-ma+1

             !performing FFT
             i=1
             CALL fftstp(lzt,nfft,n1,lot,n1,zt(1,j,1),zw(1,1,1), &
                  ftrig1,after1(i),now1(i),before1(i),-1)

             inzee=1
             DO i=2,ic1
                CALL fftstp(lot,nfft,n1,lot,n1,zw(1,1,inzee),zw(1,1,3-inzee), &
                     ftrig1,after1(i),now1(i),before1(i),-1)
                inzee=3-inzee
             ENDDO
             !output: I2,I1,j3,(jp3)

             !reverse ordering
             !input: J2,Jp2,I1,j3,(jp3)
             IF (nproc.EQ.1) THEN
                CALL unmpiswitch_downcorn(j3,nfft,Jp2stf,J2stf,lot,n1,md2,nd3,nproc,zw(1,1,inzee),zmpi2)
             ELSE
                CALL unmpiswitch_downcorn(j3,nfft,Jp2stf,J2stf,lot,n1,md2,nd3,nproc,zw(1,1,inzee),zmpi1)
             ENDIF
             ! output: I1,J2,j3,Jp2,(jp3)
          END DO
       ENDIF
    END DO

    !Interprocessor data transposition
    !input: I1,J2,j3,Jp2,(jp3)
    IF (nproc.GT.1) THEN
       !communication scheduling
       CALL mp_alltoall(zmpi1,zmpi2,n1*(md2/nproc)*(nd3/nproc),mpi_group)
       !output: I1,J2,j3,jp3,(Jp2)
    ENDIF

    !transform along z axis
    !input: I1,J2,i3,(Jp2)
    lot=ncache/(2*n3)
    DO j2=1,md2/nproc
       !this condition ensures that we manage only the interesting part for the FFT
       IF (iproc*(md2/nproc)+j2.LE.n2/2) THEN
          DO i1=1,(n1/2),lot
             ma=i1
             mb=MIN(i1+(lot-1),(n1/2))
             nfft=mb-ma+1

             !reverse ordering and repack the FFT data in order to be backward HalFFT transformed
             !input: I1,J2,i3,(Jp2)
             CALL unscramble_pack(i1,j2,lot,nfft,n1/2,n3,md2,nproc,nd3,zmpi2,zw(1,1,1),cosinarr)
             !output: I1,i3,J2,(Jp2)

             !performing FFT
             !input: I1,i3,J2,(Jp2)
             inzee=1
             DO i=1,ic3
                CALL fftstp(lot,nfft,n3/2,lot,n3/2,zw(1,1,inzee),zw(1,1,3-inzee), &
                     ftrig3,after3(i),now3(i),before3(i),-1)
                inzee=3-inzee
             ENDDO
             !output: I1,I3,J2,(Jp2)

             !calculates the exchange correlation terms locally and rebuild the output array
             CALL unfill_downcorn(md1,md3,lot,nfft,n3,zw(1,1,inzee),zf(i1,1,j2)&
                  ,scal)!,ehartreetmp)

             !integrate local pieces together
             !ehartree=ehartree+0.5_dp*ehartreetmp*hgrid**3
          END DO
       ENDIF
    END DO

    !De-allocations
    DEALLOCATE(btrig1,stat=i_all)
    DEALLOCATE(ftrig1,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(after1,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(now1,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(before1,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(btrig2,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(ftrig2,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(after2,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(now2,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(before2,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(btrig3,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(ftrig3,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(after3,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(now3,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(before3,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(zmpi2,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(zw,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(zt,stat=i_stat)
    i_all=i_all+i_stat
    DEALLOCATE(cosinarr,stat=i_stat)
    i_all=i_all+i_stat
    IF (nproc.GT.1) DEALLOCATE(zmpi1,stat=i_stat)
    IF (i_all+i_stat /= 0) THEN
       WRITE(*,*)' f_poissonsolver: problem of memory deallocation'
       STOP
    END IF

  END SUBROUTINE F_PoissonSolver

! *****************************************************************************
!> \brief ...
!> \param nfft ...
!> \param n2 ...
!> \param lot ...
!> \param n1 ...
!> \param lzt ...
!> \param zt ...
!> \param zw ...
! *****************************************************************************
  SUBROUTINE switch_upcorn(nfft,n2,lot,n1,lzt,zt,zw)
    INTEGER                                  :: nfft, n2, lot, n1, lzt
    REAL(KIND=dp)                            :: zt(2,lzt,n1), zw(2,lot,n2)

    INTEGER                                  :: i, j

! WARNING: Assuming that high frequencies are in the corners
!          and that n2 is multiple of 2
! Low frequencies

    DO j=1,nfft
       DO i=n2/2+1,n2
          zw(1,j,i)=zt(1,i-n2/2,j)
          zw(2,j,i)=zt(2,i-n2/2,j)
       END DO
    END DO
    ! High frequencies
    DO i=1,n2/2
       DO j=1,nfft
          zw(1,j,i)=0._dp
          zw(2,j,i)=0._dp
       END DO
    END DO
  END SUBROUTINE switch_upcorn

! *****************************************************************************
!> \brief ...
!> \param j3 ...
!> \param nfft ...
!> \param Jp2stb ...
!> \param J2stb ...
!> \param lot ...
!> \param n1 ...
!> \param md2 ...
!> \param nd3 ...
!> \param nproc ...
!> \param zmpi1 ...
!> \param zw ...
! *****************************************************************************
  SUBROUTINE mpiswitch_upcorn(j3,nfft,Jp2stb,J2stb,lot,n1,md2,nd3,nproc,zmpi1,zw)
    INTEGER                                  :: j3, nfft, Jp2stb, J2stb, lot, &
                                                n1, md2, nd3, nproc
    REAL(KIND=dp) :: zmpi1(2,n1/2,md2/nproc,nd3/nproc,nproc), zw(2,lot,n1)

    INTEGER                                  :: i1, j2, jp2, mfft

! WARNING: Assuming that high frequencies are in the corners
!          and that n1 is multiple of 2

    mfft=0
    Main:DO Jp2=Jp2stb,nproc
       DO J2=J2stb,md2/nproc
          mfft=mfft+1
          IF (mfft.GT.nfft) THEN
             Jp2stb=Jp2
             J2stb=J2
             EXIT Main
          ENDIF
          DO I1=1,n1/2
             zw(1,mfft,I1)=0._dp
             zw(2,mfft,I1)=0._dp
          END DO
          DO I1=n1/2+1,n1
             zw(1,mfft,I1)=zmpi1(1,I1-n1/2,J2,j3,Jp2)
             zw(2,mfft,I1)=zmpi1(2,I1-n1/2,J2,j3,Jp2)
          END DO
       END DO
       J2stb=1
    END DO Main
  END SUBROUTINE mpiswitch_upcorn

! *****************************************************************************
!> \brief ...
!> \param nfft ...
!> \param n2 ...
!> \param lot ...
!> \param n1 ...
!> \param lzt ...
!> \param zw ...
!> \param zt ...
! *****************************************************************************
  SUBROUTINE unswitch_downcorn(nfft,n2,lot,n1,lzt,zw,zt)
    INTEGER                                  :: nfft, n2, lot, n1, lzt
    REAL(KIND=dp)                            :: zw(2,lot,n2), zt(2,lzt,n1)

    INTEGER                                  :: i, j

! WARNING: Assuming that high frequencies are in the corners
!          and that n2 is multiple of 2
! Low frequencies

    DO j=1,nfft
       DO i=1,n2/2
          zt(1,i,j)=zw(1,j,i)
          zt(2,i,j)=zw(2,j,i)
       END DO
    END DO
    RETURN
  END SUBROUTINE unswitch_downcorn

! *****************************************************************************
!> \brief ...
!> \param j3 ...
!> \param nfft ...
!> \param Jp2stf ...
!> \param J2stf ...
!> \param lot ...
!> \param n1 ...
!> \param md2 ...
!> \param nd3 ...
!> \param nproc ...
!> \param zw ...
!> \param zmpi1 ...
! *****************************************************************************
  SUBROUTINE unmpiswitch_downcorn(j3,nfft,Jp2stf,J2stf,lot,n1,md2,nd3,nproc,zw,zmpi1)
    INTEGER                                  :: j3, nfft, Jp2stf, J2stf, lot, &
                                                n1, md2, nd3, nproc
    REAL(KIND=dp) :: zw(2,lot,n1), zmpi1(2,n1/2,md2/nproc,nd3/nproc,nproc)

    INTEGER                                  :: i1, j2, jp2, mfft

! WARNING: Assuming that high frequencies are in the corners
!          and that n1 is multiple of 2

    mfft=0
    Main: DO Jp2=Jp2stf,nproc
       DO J2=J2stf,md2/nproc
          mfft=mfft+1
          IF (mfft.GT.nfft) THEN
             Jp2stf=Jp2
             J2stf=J2
             EXIT Main
          ENDIF
          DO I1=1,n1/2
             zmpi1(1,I1,J2,j3,Jp2)=zw(1,mfft,I1)
             zmpi1(2,I1,J2,j3,Jp2)=zw(2,mfft,I1)
          END DO
       END DO
       J2stf=1
    END DO Main
  END SUBROUTINE unmpiswitch_downcorn

! *****************************************************************************
!> \brief (Based on suitable modifications of S.Goedecker routines)
!>      Restore data into output array, calculating in the meanwhile
!>      Hartree energy of the potential
!> \param md1 Dimensions of the undistributed part of the real grid
!> \param md3 Dimensions of the undistributed part of the real grid
!> \param lot ...
!> \param nfft number of planes
!> \param n3 (twice the) dimension of the last FFTtransform.
!> \param zw FFT work array
!> \param zf Original distributed density as well as
!>                   Distributed solution of the poisson equation (inout)
!> \param scal Needed to achieve unitarity and correct dimensions
!> \param ehartreetmp Hartree energy
!> \date February 2006
!> \author S. Goedecker, L. Genovese
!> \note Assuming that high frequencies are in the corners 
!>      and that n3 is multiple of 4 
!> 
!>  RESTRICTIONS on USAGE 
!>      Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994 
!>      Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999 
!>      Copyright (C) 2002 Stefan Goedecker, CEA Grenoble 
!>      This file is distributed under the terms of the 
!>      GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt . 
! *****************************************************************************
  SUBROUTINE F_unfill_downcorn(md1,md3,lot,nfft,n3,zw,zf,scal,ehartreetmp)
    INTEGER, INTENT(in)                      :: md1, md3, lot, nfft, n3
    REAL(KIND=dp), DIMENSION(2, lot, n3/2), &
      INTENT(in)                             :: zw
    REAL(KIND=dp), DIMENSION(md1, md3), &
      INTENT(inout)                          :: zf
    REAL(KIND=dp), INTENT(in)                :: scal
    REAL(KIND=dp), INTENT(out)               :: ehartreetmp

    INTEGER                                  :: i1, i3
    REAL(KIND=dp)                            :: pot1

    ehartreetmp=0._dp
    DO i3=1,n3/4
       DO i1=1,nfft
          pot1 = scal*zw(1,i1,i3)
          ehartreetmp =ehartreetmp + pot1* zf(i1,2*i3-1)
          zf(i1,2*i3-1)= pot1
          pot1 = scal*zw(2,i1,i3)
          ehartreetmp =ehartreetmp + pot1* zf(i1,2*i3)
          zf(i1,2*i3)= pot1
       ENDDO
    END DO
  END SUBROUTINE F_unfill_downcorn

END MODULE ps_wavelet_base
