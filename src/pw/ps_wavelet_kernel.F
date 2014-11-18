!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Creates the wavelet kernel for the wavelet based poisson solver.
!> \author Florian Schiffmann (09.2007,fschiff)
! *****************************************************************************
MODULE ps_wavelet_kernel

  
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_alltoall
  USE ps_wavelet_base,                 ONLY: scramble_unpack
  USE ps_wavelet_fft3d,                ONLY: ctrig,&
                                             ctrig_length,&
                                             fftstp
  USE ps_wavelet_scaling_function,     ONLY: scaling_function,&
                                             scf_recursion
  USE ps_wavelet_util,                 ONLY: F_FFT_dimensions,&
                                             S_FFT_dimensions
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ps_wavelet_kernel'

! *** Public data types ***

  PUBLIC :: createKernel

CONTAINS

! *****************************************************************************
!> \brief Allocate a pointer which corresponds to the zero-padded FFT slice needed for
!>     calculating the convolution with the kernel expressed in the interpolating scaling
!>     function basis. The kernel pointer is unallocated on input, allocated on output.
!> \param geocode Indicates the boundary conditions (BC) of the problem:
!>             'F' free BC, isolated systems.
!>                 The program calculates the solution as if the given density is
!>                 "alone" in R^3 space.
!>             'S' surface BC, isolated in y direction, periodic in xz plane
!>                 The given density is supposed to be periodic in the xz plane,
!>                 so the dimensions in these direction mus be compatible with the FFT
!>                 Beware of the fact that the isolated direction is y!
!>             'P' periodic BC.
!>                 The density is supposed to be periodic in all the three directions,
!>                 then all the dimensions must be compatible with the FFT.
!>                 No need for setting up the kernel.
!> \param n01 dimensions of the real space grid to be hit with the Poisson Solver
!> \param n02 dimensions of the real space grid to be hit with the Poisson Solver
!> \param n03 dimensions of the real space grid to be hit with the Poisson Solver
!> \param hx  grid spacings. For the isolated BC case for the moment they are supposed to
!>                 be equal in the three directions
!> \param hy grid spacings. For the isolated BC case for the moment they are supposed to
!>                 be equal in the three directions
!> \param hz grid spacings. For the isolated BC case for the moment they are supposed to
!>                 be equal in the three directions
!> \param itype_scf order of the interpolating scaling functions used in the decomposition
!> \param iproc ,nproc: number of process, number of processes
!> \param nproc ...
!> \param kernel pointer for the kernel FFT. Unallocated on input, allocated on output.
!>                 Its dimensions are equivalent to the region of the FFT space for which the
!>                 kernel is injective. This will divide by two each direction,
!>                 since the kernel for the zero-padded convolution is real and symmetric.
!> \param mpi_group ...
!> \param error ...
!> \date February 2007
!> \author Luigi Genovese
!> \note Due to the fact that the kernel dimensions are unknown before the calling, the kernel 
!>     must be declared as pointer in input of this routine. 
!>     To avoid that, one can properly define the kernel dimensions by adding 
!>     the nd1,nd2,nd3 arguments to the PS_dim4allocation routine, then eliminating the pointer 
!>     declaration. 
! *****************************************************************************
  SUBROUTINE createKernel(geocode,n01,n02,n03,hx,hy,hz,itype_scf,iproc,nproc,kernel,mpi_group,error)

    CHARACTER(len=1), INTENT(in)             :: geocode
    INTEGER, INTENT(in)                      :: n01, n02, n03
    REAL(KIND=dp), INTENT(in)                :: hx, hy, hz
    INTEGER, INTENT(in)                      :: itype_scf, iproc, nproc
    REAL(KIND=dp), POINTER                   :: kernel(:)
    INTEGER, INTENT(in)                      :: mpi_group
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'createKernel', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i_all, m1, m2, m3, md1, md2, &
                                                md3, n1, n2, n3, nd1, nd2, &
                                                nd3, nlimd, nlimk
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: hgrid

  hgrid=MAX(hx,hy,hz)

  IF (geocode == 'P') THEN

     CALL F_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,&
                           md1,md2,md3,nd1,nd2,nd3,nproc)

     ALLOCATE(kernel(1),stat=i_all)
     CPPostcondition(i_all==0,cp_failure_level,routineP,error,failure)
     nlimd=n2
     nlimk=0

  ELSE IF (geocode == 'S') THEN

      CALL S_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,&
                            md1,md2,md3,nd1,nd2,nd3,nproc)

     ALLOCATE(kernel(nd1*nd2*nd3/nproc),stat=i_all)
     CPPostcondition(i_all==0,cp_failure_level,routineP,error,failure)

     !the kernel must be built and scattered to all the processes

     CALL Surfaces_Kernel(n1,n2,n3,m3,nd1,nd2,nd3,hx,hz,hy,&
                          itype_scf,kernel,iproc,nproc,mpi_group,error)
     !last plane calculated for the density and the kernel

     nlimd=n2
     nlimk=n3/2+1
  ELSE IF (geocode == 'F') THEN

     !Build the Kernel

     CALL F_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,&
                           md1,md2,md3,nd1,nd2,nd3,nproc)
     ALLOCATE(kernel(nd1*nd2*nd3/nproc),stat=i_all)

     CPPostcondition(i_all==0,cp_failure_level,routineP,error,failure)
     !the kernel must be built and scattered to all the processes
     CALL Free_Kernel(n01,n02,n03,n1,n2,n3,nd1,nd2,nd3,hgrid,&
                      itype_scf,iproc,nproc,kernel,mpi_group,error)

     !last plane calculated for the density and the kernel
     nlimd=n2/2
     nlimk=n3/2+1

  ELSE

     CALL cp_assert(.FALSE.,cp_assertion_failed,&
          cp_failure_level,routineP,&
          "No wavelet based poisson solver for given geometry",error,failure)

  END IF
!!!  IF (iproc==0) THEN
!!!     write(*,*)'done.'
!!!     write(*,'(1x,a,i0)') 'Allocate words for kernel ',nd1*nd2*nd3/nproc
!!!     !print the load balancing of the different dimensions on screen
!!!     write(*,'(1x,a,3(i5))')'Grid Dimensions:',n01,n02,n03
!!!     if (nproc > 1) then
!!!        write(*,'(1x,a,3(i5),a,3(i5),a,3(i5))')&
!!!             'Memory occ. per proc.  Density',md1,md3,md2/nproc,'   Kernel',nd1,nd2,nd3/nproc
!!!        write(*,'(1x,a)')'Load Balancing--------------------------------------------'
!!!        jhd=1000
!!!        jzd=1000
!!!        npd=0
!!!        load_balancing: do jproc=0,nproc-1
!!!           !print *,'jproc,jfull=',jproc,jproc*md2/nproc,(jproc+1)*md2/nproc
!!!           if ((jproc+1)*md2/nproc <= nlimd) then
!!!              jfd=jproc
!!!           else if (jproc*md2/nproc <= nlimd) then
!!!              jhd=jproc
!!!              npd=nint(real(nlimd-(jproc)*md2/nproc,KIND=dp)/real(md2/nproc,KIND=dp)*100._dp)
!!!           else
!!!              jzd=jproc
!!!              exit load_balancing
!!!           end if
!!!        end do load_balancing
!!!        write(*,'(1x,a,i3,a)')'LB_den        : processors   0  -',jfd,' work at 100%'
!!!        if (jfd < nproc-1) write(*,'(1x,a,i3,a,i3,1a)')'                processor     ',jhd,&
!!!             '    work at ',npd,'%'
!!!        if (jhd < nproc-1) write(*,'(1x,a,i3,1a,i3,a)')'                processors ',&
!!!             jzd,'  -',nproc-1,' work at   0%'
!!!        jhk=1000
!!!        jzk=1000
!!!        npk=0
!!!        if (geocode /= 'P') then
!!!           load_balancingk: do jproc=0,nproc-1
!!!              !print *,'jproc,jfull=',jproc,jproc*nd3/nproc,(jproc+1)*nd3/nproc
!!!              if ((jproc+1)*nd3/nproc <= nlimk) then
!!!                 jfk=jproc
!!!              else if (jproc*nd3/nproc <= nlimk) then
!!!                 jhk=jproc
!!!                 npk=nint(real(nlimk-(jproc)*nd3/nproc,KIND=dp)/real(nd3/nproc,KIND=dp)*100._dp)
!!!              else
!!!                 jzk=jproc
!!!                 exit load_balancingk
!!!              end if
!!!           end do load_balancingk
!!!           write(*,'(1x,a,i3,a)')'LB_ker        : processors   0  -',jfk,' work at 100%'
!!!           if (jfk < nproc-1) write(*,'(1x,a,i3,a,i3,1a)')'                processor     ',jhk,&
!!!                '    work at ',npk,'%'
!!!           if (jhk < nproc-1) write(*,'(1x,a,i3,1a,i3,a)')'                processors ',jzk,'  -',nproc-1,&
!!!                ' work at   0%'
!!!        end if
!!!        write(*,'(1x,a)')'The LB per processor is 1/3 LB_den + 2/3 LB_ker-----------'
!!!     end if
!!!
!!!  END IF
END SUBROUTINE createKernel

! *****************************************************************************
!> \brief Build the kernel of the Poisson operator with
!>     surfaces Boundary conditions
!>     in an interpolating scaling functions basis.
!>     Beware of the fact that the nonperiodic direction is y!
!> \param n1 Dimensions for the FFT
!> \param n2 Dimensions for the FFT
!> \param n3 Dimensions for the FFT
!> \param m3 Actual dimension in non-periodic direction
!> \param nker1 Dimensions of the kernel (nker3=n3/2+1) nker(1,2)=n(1,2)/2+1
!> \param nker2 Dimensions of the kernel (nker3=n3/2+1) nker(1,2)=n(1,2)/2+1
!> \param nker3 Dimensions of the kernel (nker3=n3/2+1) nker(1,2)=n(1,2)/2+1
!> \param h1 Mesh steps in the three dimensions
!> \param h2 Mesh steps in the three dimensions
!> \param h3 Mesh steps in the three dimensions
!> \param itype_scf Order of the scaling function
!> \param karray output array
!> \param iproc Number of process
!> \param nproc number of processes
!> \param mpi_group ...
!> \param error ...
!> \date October 2006
!> \author L. Genovese
! *****************************************************************************
SUBROUTINE Surfaces_Kernel(n1,n2,n3,m3,nker1,nker2,nker3,h1,h2,h3,&
                           itype_scf,karray,iproc,nproc,mpi_group,error)

    INTEGER, INTENT(in)                      :: n1, n2, n3, m3, nker1, nker2, &
                                                nker3
    REAL(KIND=dp), INTENT(in)                :: h1, h2, h3
    INTEGER, INTENT(in)                      :: itype_scf, nproc, iproc
    REAL(KIND=dp), &
      DIMENSION(nker1, nker2, nker3/nproc), &
      INTENT(out)                            :: karray
    INTEGER, INTENT(in)                      :: mpi_group
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'Surfaces_Kernel', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: n_points = 2**6, &
                                                ncache_optimal = 8*1024, &
                                                timing_flag = 0

    INTEGER :: i, i1, i2, i3, i_all, i_stat, ic, iend, imu, ind1, ind2, &
      inzee, ipolyord, ireim, istart, j2, j2nd, j2st, jnd1, jp2, jreim, &
      n_cell, n_range, n_scf, nact2, ncache, nfft, num_of_mus, shift
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: after, before, now
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: a, b, c, cp, d, diff, dx, &
                                                feI, feR, foI, foR, fR, mu1, &
                                                pi, pion, ponx, pony, sp, &
                                                value, x
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: kernel_scf, x_scf, y_scf
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: btrig, cossinarr
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: halfft_cache, kernel
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: kernel_mpi
    REAL(KIND=dp), DIMENSION(9, 8)           :: cpol

!Better if higher (1024 points are enough 10^{-14}: 2*itype_scf*n_points)
!  include "perfdata.inc"
!FFT arrays
!coefficients for the polynomial interpolation
!assign the values of the coefficients

  karray = 0.0_dp
  cpol(:,:)=1._dp

  cpol(1,2)=.25_dp

  cpol(1,3)=1._dp/3._dp

  cpol(1,4)=7._dp/12._dp
  cpol(2,4)=8._dp/3._dp

  cpol(1,5)=19._dp/50._dp
  cpol(2,5)=3._dp/2._dp

  cpol(1,6)=41._dp/272._dp
  cpol(2,6)=27._dp/34._dp
  cpol(3,6)=27._dp/272._dp

  cpol(1,7)=751._dp/2989._dp
  cpol(2,7)=73._dp/61._dp
  cpol(3,7)=27._dp/61._dp

  cpol(1,8)=-989._dp/4540._dp
  cpol(2,8)=-1472._dp/1135._dp
  cpol(3,8)=232._dp/1135._dp
  cpol(4,8)=-2624._dp/1135._dp

  !renormalize values
  cpol(1,1)=.5_dp*cpol(1,1)
  cpol(1:2,2)=2._dp/3._dp*cpol(1:2,2)
  cpol(1:2,3)=3._dp/8._dp*cpol(1:2,3)
  cpol(1:3,4)=2._dp/15._dp*cpol(1:3,4)
  cpol(1:3,5)=25._dp/144._dp*cpol(1:3,5)
  cpol(1:4,6)=34._dp/105._dp*cpol(1:4,6)
  cpol(1:4,7)=2989._dp/17280._dp*cpol(1:4,7)
  cpol(1:5,8)=-454._dp/2835._dp*cpol(1:5,8)

  !assign the complete values
  cpol(2,1)=cpol(1,1)

  cpol(3,2)=cpol(1,2)

  cpol(3,3)=cpol(2,3)
  cpol(4,3)=cpol(1,3)

  cpol(4,4)=cpol(2,4)
  cpol(5,4)=cpol(1,4)

  cpol(4,5)=cpol(3,5)
  cpol(5,5)=cpol(2,5)
  cpol(6,5)=cpol(1,5)

  cpol(5,6)=cpol(3,6)
  cpol(6,6)=cpol(2,6)
  cpol(7,6)=cpol(1,6)

  cpol(5,7)=cpol(4,7)
  cpol(6,7)=cpol(3,7)
  cpol(7,7)=cpol(2,7)
  cpol(8,7)=cpol(1,7)

  cpol(6,8)=cpol(4,8)
  cpol(7,8)=cpol(3,8)
  cpol(8,8)=cpol(2,8)
  cpol(9,8)=cpol(1,8)

  !Number of integration points : 2*itype_scf*n_points
  n_scf=2*itype_scf*n_points
  !Allocations
  i_all = 0
  ALLOCATE(x_scf(0:n_scf),stat=i_stat)
  i_all = i_all + i_stat
  ALLOCATE(y_scf(0:n_scf),stat=i_stat)
  i_all = i_all + i_stat
  CPPostcondition(i_all==0,cp_failure_level,routineP,error,failure)

  !Build the scaling function
  CALL scaling_function(itype_scf,n_scf,n_range,x_scf,y_scf)
  !Step grid for the integration
  dx = REAL(n_range,KIND=dp)/REAL(n_scf,KIND=dp)
  !Extend the range (no more calculations because fill in by 0._dp)
  n_cell = m3
  n_range = MAX(n_cell,n_range)

  !Allocations
  ncache=ncache_optimal
  !the HalFFT must be performed only in the third dimension,
  !and nker3=n3/2+1, hence
  IF (ncache <= (nker3-1)*4) ncache=nker3-1*4

  !enlarge the second dimension of the kernel to be compatible with nproc
  nact2=nker2
  enlarge_ydim: DO
     IF (nproc*(nact2/nproc) /= nact2) THEN
        nact2=nact2+1
     ELSE
        EXIT enlarge_ydim
     END IF
  END DO enlarge_ydim

  !array for the MPI procedure
  ALLOCATE(kernel(nker1,nact2/nproc,nker3),stat=i_stat)
  i_all = i_all + i_stat
  ALLOCATE(kernel_mpi(nker1,nact2/nproc,nker3/nproc,nproc),stat=i_stat)
  i_all = i_all + i_stat
  ALLOCATE(kernel_scf(n_range),stat=i_stat)
  i_all = i_all + i_stat
  ALLOCATE(halfft_cache(2,ncache/4,2),stat=i_stat)
  i_all = i_all + i_stat
  ALLOCATE(cossinarr(2,n3/2-1),stat=i_stat)
  i_all = i_all + i_stat
  ALLOCATE(btrig(2,ctrig_length),stat=i_stat)
  i_all = i_all + i_stat
  ALLOCATE(after(7),stat=i_stat)
  i_all = i_all + i_stat
  ALLOCATE(now(7),stat=i_stat)
  i_all = i_all + i_stat
  ALLOCATE(before(7),stat=i_stat)
  i_all = i_all + i_stat

  CPPostcondition(i_all==0,cp_failure_level,routineP,error,failure)
  !constants
  pi=4._dp*ATAN(1._dp)

  !arrays for the halFFT
  CALL ctrig(n3/2,btrig,after,before,now,1,ic)

  !build the phases for the HalFFT reconstruction
  pion=2._dp*pi/REAL(n3,KIND=dp)
  DO i3=2,n3/2
     x=REAL(i3-1,KIND=dp)*pion
     cossinarr(1,i3-1)= COS(x)
     cossinarr(2,i3-1)=-SIN(x)
  END DO
  !kernel=0._dp
  !kernel_mpi=0._dp

  !calculate the limits of the FFT calculations
  !that can be performed in a row remaining inside the cache
  num_of_mus=ncache/(2*n3)

  diff=0._dp
  !order of the polynomial to be used for integration (must be a power of two)
  ipolyord=8 !this part should be incorporated inside the numerical integration
  !here we have to choice the piece of the x-y grid to cover

  !let us now calculate the fraction of mu that will be considered
  j2st=iproc*(nact2/nproc)
  j2nd=MIN((iproc+1)*(nact2/nproc),n2/2+1)

  DO ind2=(n1/2+1)*j2st+1,(n1/2+1)*j2nd,num_of_mus
     istart=ind2
     iend=MIN(ind2+(num_of_mus-1),(n1/2+1)*j2nd)
     nfft=iend-istart+1
     shift=0

     !initialization of the interesting part of the cache array
     halfft_cache(:,:,:)=0._dp

     IF (istart == 1) THEN
        !i2=1
        shift=1

        CALL calculates_green_opt_muzero(n_range,n_scf,ipolyord,x_scf,y_scf,&
             cpol(1,ipolyord),dx,kernel_scf,error)

        !copy of the first zero value
        halfft_cache(1,1,1)=0._dp

        DO i3=1,m3

           value=0.5_dp*h3*kernel_scf(i3)
           !index in where to copy the value of the kernel
           CALL indices(ireim,num_of_mus,n3/2+i3,1,ind1)
           !index in where to copy the symmetric value
           CALL indices(jreim,num_of_mus,n3/2+2-i3,1,jnd1)
           halfft_cache(ireim,ind1,1) = value
           halfft_cache(jreim,jnd1,1) = value

        END DO

     END IF

     loopimpulses : DO imu=istart+shift,iend

        !here there is the value of mu associated to hgrid
        !note that we have multiplicated mu for hgrid to be comparable
        !with mu0ref

        !calculate the proper value of mu taking into account the periodic dimensions
        !corresponding value of i1 and i2
        i1=MOD(imu,n1/2+1)
        IF (i1==0) i1=n1/2+1
        i2=(imu-i1)/(n1/2+1)+1
        ponx=REAL(i1-1,KIND=dp)/REAL(n1,KIND=dp)
        pony=REAL(i2-1,KIND=dp)/REAL(n2,KIND=dp)

        mu1=2._dp*pi*SQRT((ponx/h1)**2+(pony/h2)**2)*h3

        CALL calculates_green_opt(n_range,n_scf,itype_scf,ipolyord,x_scf,y_scf,&
             cpol(1,ipolyord),mu1,dx,kernel_scf,error)

        !readjust the coefficient and define the final kernel

        !copy of the first zero value
        halfft_cache(1,imu-istart+1,1) = 0._dp
        DO i3=1,m3
           value=-0.5_dp*h3/mu1*kernel_scf(i3)
           !write(80,*)mu1,i3,kernel_scf(i03)
           !index in where to copy the value of the kernel
           CALL indices(ireim,num_of_mus,n3/2+i3,imu-istart+1,ind1)
           !index in where to copy the symmetric value
           CALL indices(jreim,num_of_mus,n3/2+2-i3,imu-istart+1,jnd1)
           halfft_cache(ireim,ind1,1)=value
           halfft_cache(jreim,jnd1,1)=value
        END DO

     END DO loopimpulses

     !now perform the FFT of the array in cache
     inzee=1
     DO i=1,ic
        CALL fftstp(num_of_mus,nfft,n3/2,num_of_mus,n3/2,&
             halfft_cache(1,1,inzee),halfft_cache(1,1,3-inzee),&
             btrig,after(i),now(i),before(i),1)
        inzee=3-inzee
     ENDDO
     !assign the values of the FFT array
     !and compare with the good results
     DO imu=istart,iend

        !corresponding value of i1 and i2
        i1=MOD(imu,n1/2+1)
        IF (i1==0) i1=n1/2+1
        i2=(imu-i1)/(n1/2+1)+1

        j2=i2-j2st

        a=halfft_cache(1,imu-istart+1,inzee)
        b=halfft_cache(2,imu-istart+1,inzee)
        kernel(i1,j2,1)=a+b
        kernel(i1,j2,n3/2+1)=a-b

        DO i3=2,n3/2
           ind1=imu-istart+1+num_of_mus*(i3-1)
           jnd1=imu-istart+1+num_of_mus*(n3/2+2-i3-1)
           cp=cossinarr(1,i3-1)
           sp=cossinarr(2,i3-1)
           a=halfft_cache(1,ind1,inzee)
           b=halfft_cache(2,ind1,inzee)
           c=halfft_cache(1,jnd1,inzee)
           d=halfft_cache(2,jnd1,inzee)
           feR=.5_dp*(a+c)
           feI=.5_dp*(b-d)
           foR=.5_dp*(a-c)
           foI=.5_dp*(b+d)
           fR=feR+cp*foI-sp*foR
           kernel(i1,j2,i3)=fR
        END DO
     END DO

  END DO

  !give to each processor a slice of the third dimension
  IF (nproc > 1) THEN
    CALL mp_alltoall(kernel,&!nker1*(nact2/nproc)*(nker3/nproc), &
         kernel_mpi,nker1*(nact2/nproc)*(nker3/nproc), &
         mpi_group)

     DO jp2=1,nproc
        DO i3=1,nker3/nproc
           DO i2=1,nact2/nproc
              j2=i2+(jp2-1)*(nact2/nproc)
              IF (j2 <= nker2) THEN
                 DO i1=1,nker1
                    karray(i1,j2,i3)=&
                         kernel_mpi(i1,i2,i3,jp2)
                 END DO
              END IF
           END DO
        END DO
     END DO

  ELSE
     karray(1:nker1,1:nker2,1:nker3)=kernel(1:nker1,1:nker2,1:nker3)
  ENDIF

  !De-allocations
  DEALLOCATE(kernel,stat=i_all)
  DEALLOCATE(kernel_mpi,stat=i_stat)
  i_all=i_all+i_stat
  DEALLOCATE(btrig,stat=i_stat)
  i_all=i_all+i_stat
  DEALLOCATE(after,stat=i_stat)
  i_all=i_all+i_stat
  DEALLOCATE(now,stat=i_stat)
  i_all=i_all+i_stat
  DEALLOCATE(before,stat=i_stat)
  i_all=i_all+i_stat
  DEALLOCATE(halfft_cache,stat=i_stat)
  i_all=i_all+i_stat
  DEALLOCATE(kernel_scf,stat=i_stat)
  i_all=i_all+i_stat
  DEALLOCATE(x_scf,stat=i_stat)
  i_all=i_all+i_stat
  DEALLOCATE(y_scf,stat=i_stat)

  CPPostcondition(i_all==0,cp_failure_level,routineP,error,failure)

END SUBROUTINE Surfaces_Kernel

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param n_scf ...
!> \param itype_scf ...
!> \param intorder ...
!> \param xval ...
!> \param yval ...
!> \param c ...
!> \param mu ...
!> \param hres ...
!> \param g_mu ...
!> \param error ...
! *****************************************************************************
SUBROUTINE calculates_green_opt(n,n_scf,itype_scf,intorder,xval,yval,c,mu,hres,g_mu,error)
    INTEGER, INTENT(in)                      :: n, n_scf, itype_scf, intorder
    REAL(KIND=dp), DIMENSION(0:n_scf), &
      INTENT(in)                             :: xval, yval
    REAL(KIND=dp), DIMENSION(intorder+1), &
      INTENT(in)                             :: c
    REAL(KIND=dp), INTENT(in)                :: mu, hres
    REAL(KIND=dp), DIMENSION(n), INTENT(out) :: g_mu
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculates_green_opt', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: mu_max = 0.2_dp

    INTEGER                                  :: i, i_all, i_stat, iend, &
                                                ikern, ivalue, izero, n_iter, &
                                                nrec
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: f, filter, fl, fr, gleft, &
                                                gltmp, gright, grtmp, mu0, &
                                                ratio, x, x0, x1
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: green, green1

  g_mu = 0.0_dp
  !We calculate the number of iterations to go from mu0 to mu0_ref
  IF (mu <= mu_max) THEN
     n_iter = 0
     mu0 = mu
  ELSE
     n_iter=1
     loop_iter: DO
        ratio=REAL(2**n_iter,KIND=dp)
        mu0=mu/ratio
        IF (mu0 <= mu_max) THEN
           EXIT loop_iter
        END IF
        n_iter=n_iter+1
     END DO loop_iter
  END IF

  !dimension needed for the correct calculation of the recursion
  nrec=2**n_iter*n

  ALLOCATE(green(-nrec:nrec),stat=i_all)
  CPPostcondition(i_all==0,cp_failure_level,routineP,error,failure)

  !initialization of the branching value
  ikern=0
  izero=0
  initialization: DO
     IF (xval(izero)>=REAL(ikern,KIND=dp) .OR. izero==n_scf) EXIT initialization
     izero=izero+1
  END DO initialization
  green=0._dp
  !now perform the interpolation in right direction
  ivalue=izero
  gright=0._dp
  loop_right: DO
     IF(ivalue >= n_scf-intorder-1) EXIT loop_right
     DO i=1,intorder+1
        x=xval(ivalue)-REAL(ikern,KIND=dp)
        f=yval(ivalue)*EXP(-mu0*x)
        filter=intorder*c(i)
        gright=gright+filter*f
        ivalue=ivalue+1
     END DO
     ivalue=ivalue-1
  END DO loop_right
  iend=n_scf-ivalue
  DO i=1,iend
     x=xval(ivalue)-REAL(ikern,KIND=dp)
     f=yval(ivalue)*EXP(-mu0*x)
     filter=intorder*c(i)
     gright=gright+filter*f
     ivalue=ivalue+1
  END DO
  gright=hres*gright

  !the scaling function is symmetric, so the same for the other direction
  gleft=gright

  green(ikern)=gleft+gright

  !now the loop until the last value
  DO ikern=1,nrec
     gltmp=0._dp
     grtmp=0._dp
     ivalue=izero
     x0=xval(izero)
     loop_integration: DO
        IF (izero==n_scf)  EXIT loop_integration
        DO i=1,intorder+1
           x=xval(ivalue)
           fl=yval(ivalue)*EXP(mu0*x)
           fr=yval(ivalue)*EXP(-mu0*x)
           filter=intorder*c(i)
           gltmp=gltmp+filter*fl
           grtmp=grtmp+filter*fr
           ivalue=ivalue+1
           IF (xval(izero)>=REAL(ikern,KIND=dp) .OR. izero==n_scf) THEN
              x1=xval(izero)
              EXIT loop_integration
           END IF
           izero=izero+1
        END DO
        ivalue=ivalue-1
        izero=izero-1
     END DO loop_integration
     gleft=EXP(-mu0)*(gleft+hres*EXP(-mu0*REAL(ikern-1,KIND=dp))*gltmp)
     IF (izero == n_scf) THEN
        gright=0._dp
     ELSE
        gright=EXP(mu0)*(gright-hres*EXP(mu0*REAL(ikern-1,KIND=dp))*grtmp)
     END IF
     green(ikern)=gleft+gright
     green(-ikern)=gleft+gright
     IF (ABS(green(ikern)) <= 1.e-20_dp) THEN
        nrec=ikern
        EXIT
     END IF
     !print *,ikern,izero,n_scf,gltmp,grtmp,gleft,gright,x0,x1,green(ikern)
  END DO
  !now we must calculate the recursion
  ALLOCATE(green1(-nrec:nrec),stat=i_all)
  CPPostcondition(i_all==0,cp_failure_level,routineP,error,failure)

  !Start the iteration to go from mu0 to mu
  CALL scf_recursion(itype_scf,n_iter,nrec,green(-nrec),green1(-nrec))

  DO i=1,MIN(n,nrec)
     g_mu(i)=green(i-1)
  END DO
  DO i=MIN(n,nrec)+1,n
     g_mu(i)=0._dp
  END DO

  DEALLOCATE(green,stat=i_all)
  DEALLOCATE(green1,stat=i_stat)
  CPPostcondition(i_stat+i_all==0,cp_failure_level,routineP,error,failure)

END SUBROUTINE calculates_green_opt

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param n_scf ...
!> \param intorder ...
!> \param xval ...
!> \param yval ...
!> \param c ...
!> \param hres ...
!> \param green ...
!> \param error ...
! *****************************************************************************
SUBROUTINE calculates_green_opt_muzero(n,n_scf,intorder,xval,yval,c,hres,green,error)
    INTEGER, INTENT(in)                      :: n, n_scf, intorder
    REAL(KIND=dp), DIMENSION(0:n_scf), &
      INTENT(in)                             :: xval, yval
    REAL(KIND=dp), DIMENSION(intorder+1), &
      INTENT(in)                             :: c
    REAL(KIND=dp), INTENT(in)                :: hres
    REAL(KIND=dp), DIMENSION(n), INTENT(out) :: green
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: i, iend, ikern, ivalue, izero
    REAL(KIND=dp)                            :: c0, c1, filter, gl0, gl1, &
                                                gr0, gr1, x, y

!initialization of the branching value

  ikern=0
  izero=0
  initialization: DO
     IF (xval(izero)>=REAL(ikern,KIND=dp) .OR. izero==n_scf) EXIT initialization
     izero=izero+1
  END DO initialization
  green=0._dp
  !first case, ikern=0
  !now perform the interpolation in right direction
  ivalue=izero
  gr1=0._dp
  loop_right: DO
     IF(ivalue >= n_scf-intorder-1) EXIT loop_right
     DO i=1,intorder+1
        x=xval(ivalue)
        y=yval(ivalue)
        filter=intorder*c(i)
        gr1=gr1+filter*x*y
        ivalue=ivalue+1
     END DO
     ivalue=ivalue-1
  END DO loop_right
  iend=n_scf-ivalue
  DO i=1,iend
     x=xval(ivalue)
     y=yval(ivalue)
     filter=intorder*c(i)
     gr1=gr1+filter*x*y
     ivalue=ivalue+1
  END DO
  gr1=hres*gr1
  !the scaling function is symmetric
  gl1=-gr1
  gl0=0.5_dp
  gr0=0.5_dp

  green(1)=2._dp*gr1

  !now the loop until the last value
  DO ikern=1,n-1
     c0=0._dp
     c1=0._dp
     ivalue=izero
     loop_integration: DO
        IF (izero==n_scf)  EXIT loop_integration
        DO i=1,intorder+1
           x=xval(ivalue)
           y=yval(ivalue)
           filter=intorder*c(i)
           c0=c0+filter*y
           c1=c1+filter*y*x
           ivalue=ivalue+1
           IF (xval(izero)>=REAL(ikern,KIND=dp) .OR. izero==n_scf) THEN
              EXIT loop_integration
           END IF
           izero=izero+1
        END DO
        ivalue=ivalue-1
        izero=izero-1
     END DO loop_integration
     c0=hres*c0
     c1=hres*c1

     gl0=gl0+c0
     gl1=gl1+c1
     gr0=gr0-c0
     gr1=gr1-c1
     !general case
     green(ikern+1)=REAL(ikern,KIND=dp)*(gl0-gr0)+gr1-gl1
     !print *,ikern,izero,n_scf,gltmp,grtmp,gleft,gright,x0,x1,green(ikern)
  END DO

END SUBROUTINE calculates_green_opt_muzero

! *****************************************************************************
!> \brief ...
!> \param var_realimag ...
!> \param nelem ...
!> \param intrn ...
!> \param extrn ...
!> \param index ...
! *****************************************************************************
SUBROUTINE indices(var_realimag,nelem,intrn,extrn,index)

    INTEGER, INTENT(out)                     :: var_realimag
    INTEGER, INTENT(in)                      :: nelem, intrn, extrn
    INTEGER, INTENT(out)                     :: index

    INTEGER                                  :: i

!real or imaginary part

  var_realimag=2-MOD(intrn,2)
!actual index over half the length

  i=(intrn+1)/2
  !check
  IF (2*(i-1)+var_realimag /= intrn) THEN
     PRINT *,'error, index=',intrn,'var_realimag=',var_realimag,'i=',i
  END IF
  !complete index to be assigned
  index=extrn+nelem*(i-1)

END SUBROUTINE indices

! *****************************************************************************
!> \brief Build the kernel of a gaussian function
!>     for interpolating scaling functions.
!>     Do the parallel HalFFT of the symmetrized function and stores into
!>     memory only 1/8 of the grid divided by the number of processes nproc
!>
!>     Build the kernel (karray) of a gaussian function
!>     for interpolating scaling functions
!>     $$ K(j) = \sum_k \omega_k \int \int \phi(x) g_k(x'-x) \delta(x'- j) dx dx' $$
!> \param n01 Mesh dimensions of the density
!> \param n02 Mesh dimensions of the density
!> \param n03 Mesh dimensions of the density
!> \param nfft1 Dimensions of the FFT grid (HalFFT in the third direction)
!> \param nfft2 Dimensions of the FFT grid (HalFFT in the third direction)
!> \param nfft3 Dimensions of the FFT grid (HalFFT in the third direction)
!> \param n1k Dimensions of the kernel FFT
!> \param n2k Dimensions of the kernel FFT
!> \param n3k Dimensions of the kernel FFT
!> \param hgrid Mesh step
!> \param itype_scf Order of the scaling function (8,14,16)
!> \param iproc ...
!> \param nproc ...
!> \param karray ...
!> \param mpi_group ...
!> \param error ...
!> \date February 2006
!> \author T. Deutsch, L. Genovese
! *****************************************************************************
SUBROUTINE Free_Kernel(n01,n02,n03,nfft1,nfft2,nfft3,n1k,n2k,n3k,&
     hgrid,itype_scf,iproc,nproc,karray,mpi_group,error)

    INTEGER, INTENT(in)                      :: n01, n02, n03, nfft1, nfft2, &
                                                nfft3, n1k, n2k, n3k
    REAL(KIND=dp), INTENT(in)                :: hgrid
    INTEGER, INTENT(in)                      :: itype_scf, iproc, nproc
    REAL(KIND=dp), &
      DIMENSION(n1k, n2k, n3k/nproc), &
      INTENT(out)                            :: karray
    INTEGER, INTENT(in)                      :: mpi_group
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'Free_Kernel', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: n_gauss = 89, n_points = 2**6
    REAL(KIND=dp), PARAMETER                 :: p0_ref = 1._dp

    INTEGER :: i, i01, i02, i03, i1, i2, i3, i_all, i_gauss, i_kern, i_stat, &
      iend, istart, istart1, n1h, n2h, n3h, n_cell, n_iter, n_range, n_scf, &
      nker1, nker2, nker3
    LOGICAL                                  :: failure
    REAL(KIND=dp) :: a1, a2, a3, a_range, absci, acc_gauss, dr_gauss, dx, &
      factor, factor2, kern, p0_cell, p0gauss, pgauss, ur_gauss
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: kern_1_scf, kernel_scf, &
                                                x_scf, y_scf
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: kp
    REAL(KIND=dp), DIMENSION(n_gauss)        :: p_gauss, w_gauss

!Do not touch !!!!
!Better if higher (1024 points are enough 10^{-14}: 2*itype_scf*n_points)
!Better p_gauss for calculation
!(the support of the exponential should be inside [-n_range/2,n_range/2])
!Number of integration points : 2*itype_scf*n_points

 n_scf=2*itype_scf*n_points
 !Set karray
 karray = 0.0_dp
 !here we must set the dimensions for the fft part, starting from the nfft
 !remember that actually nfft2 is associated to n03 and viceversa

 !dimensions that define the center of symmetry
 n1h=nfft1/2
 n2h=nfft2/2
 n3h=nfft3/2

 !Auxiliary dimensions only for building the FFT part
 nker1=nfft1
 nker2=nfft2
 nker3=nfft3/2+1

 !adjusting the last two dimensions to be multiples of nproc
 DO
    IF(MODULO(nker2,nproc) == 0) EXIT
    nker2=nker2+1
 END DO
 DO
    IF(MODULO(nker3,nproc) == 0) EXIT
    nker3=nker3+1
 END DO

 !this will be the array of the kernel in the real space
 ALLOCATE(kp(n1h+1,n3h+1,nker2/nproc),stat=i_all)
 CPPostcondition(i_all==0,cp_failure_level,routineP,error,failure)

 !defining proper extremes for the calculation of the
 !local part of the kernel

 istart=iproc*nker2/nproc+1
 iend=MIN((iproc+1)*nker2/nproc,n2h+n03)

 istart1=istart
 IF(iproc .EQ. 0) istart1=n2h-n03+2

 !Allocations
 i_all = 0
 ALLOCATE(x_scf(0:n_scf),stat=i_stat)
 i_all = i_all + i_stat
 ALLOCATE(y_scf(0:n_scf),stat=i_stat)
 i_all = i_all + i_stat
 CPPostcondition(i_all==0,cp_failure_level,routineP,error,failure)

 !Build the scaling function
 CALL scaling_function(itype_scf,n_scf,n_range,x_scf,y_scf)
 !Step grid for the integration
 dx = REAL(n_range,KIND=dp)/REAL(n_scf,KIND=dp)
 !Extend the range (no more calculations because fill in by 0._dp)
 n_cell = MAX(n01,n02,n03)
 n_range = MAX(n_cell,n_range)

 !Allocations
 ALLOCATE(kernel_scf(-n_range:n_range),stat=i_stat)
 i_all = i_all + i_stat
 ALLOCATE(kern_1_scf(-n_range:n_range),stat=i_stat)
 i_all = i_all + i_stat
 CPPostcondition(i_all==0,cp_failure_level,routineP,error,failure)

 !Lengthes of the box (use FFT dimension)
 a1 = hgrid * REAL(n01,KIND=dp)
 a2 = hgrid * REAL(n02,KIND=dp)
 a3 = hgrid * REAL(n03,KIND=dp)

 x_scf(:) = hgrid * x_scf(:)
 y_scf(:) = 1._dp/hgrid * y_scf(:)
 dx = hgrid * dx
 !To have a correct integration
 p0_cell = p0_ref/(hgrid*hgrid)

 !Initialization of the gaussian (Beylkin)
 CALL gequad(n_gauss,p_gauss,w_gauss,ur_gauss,dr_gauss,acc_gauss)
 !In order to have a range from a_range=sqrt(a1*a1+a2*a2+a3*a3)
 !(biggest length in the cube)
 !We divide the p_gauss by a_range**2 and a_gauss by a_range
 a_range = SQRT(a1*a1+a2*a2+a3*a3)
 factor = 1._dp/a_range
 !factor2 = factor*factor
 factor2 = 1._dp/(a1*a1+a2*a2+a3*a3)
 DO i_gauss=1,n_gauss
    p_gauss(i_gauss) = factor2*p_gauss(i_gauss)
 END DO
 DO i_gauss=1,n_gauss
    w_gauss(i_gauss) = factor*w_gauss(i_gauss)
 END DO

 kp(:,:,:)=0._dp
 !Use in this order (better for accuracy).
 loop_gauss: DO i_gauss=n_gauss,1,-1
    !Gaussian
    pgauss = p_gauss(i_gauss)

    !We calculate the number of iterations to go from pgauss to p0_ref
    n_iter = NINT((LOG(pgauss) - LOG(p0_cell))/LOG(4._dp))
    IF (n_iter <= 0)THEN
       n_iter = 0
       p0gauss = pgauss
    ELSE
       p0gauss = pgauss/4._dp**n_iter
    END IF

    !Stupid integration
    !Do the integration with the exponential centered in i_kern
    kernel_scf(:) = 0._dp
    DO i_kern=0,n_range
       kern = 0._dp
       DO i=0,n_scf
          absci = x_scf(i) - REAL(i_kern,KIND=dp)*hgrid
          absci = absci*absci
          kern = kern + y_scf(i)*EXP(-p0gauss*absci)*dx
       END DO
       kernel_scf(i_kern) = kern
       kernel_scf(-i_kern) = kern
       IF (ABS(kern) < 1.e-18_dp) THEN
          !Too small not useful to calculate
          EXIT
       END IF
    END DO

    !Start the iteration to go from p0gauss to pgauss
    CALL scf_recursion(itype_scf,n_iter,n_range,kernel_scf,kern_1_scf)

    !Add to the kernel (only the local part)

    DO i3=istart1,iend
       i03 = i3 - n2h -1
       DO i2=1,n02
          i02 = i2-1
          DO i1=1,n01
             i01 = i1-1
             kp(i1,i2,i3-istart+1) = kp(i1,i2,i3-istart+1) + w_gauss(i_gauss)* &
                  kernel_scf(i01)*kernel_scf(i02)*kernel_scf(i03)
          END DO
       END DO
    END DO

 END DO loop_gauss

 !De-allocations
 DEALLOCATE(kernel_scf,stat=i_all)
 DEALLOCATE(kern_1_scf,stat=i_stat)
 i_all=i_all+i_stat
 DEALLOCATE(x_scf,stat=i_stat)
 i_all=i_all+i_stat
 DEALLOCATE(y_scf,stat=i_stat)
 CPPostcondition(i_all+i_stat==0,cp_failure_level,routineP,error,failure)

!!!!END KERNEL CONSTRUCTION

!!$ if(iproc .eq. 0) print *,"Do a 3D PHalFFT for the kernel"

 CALL kernelfft(nfft1,nfft2,nfft3,nker1,nker2,nker3,n1k,n2k,n3k,nproc,iproc,&
      kp,karray,mpi_group,error)

 !De-allocations
 DEALLOCATE(kp,stat=i_all)
 CPPostcondition(i_all==0,cp_failure_level,routineP,error,failure)

END SUBROUTINE Free_Kernel

! *****************************************************************************
!> \brief ...
!> \param n1 ...
!> \param n3 ...
!> \param lot ...
!> \param nfft ...
!> \param i1 ...
!> \param zf ...
!> \param zw ...
! *****************************************************************************
SUBROUTINE inserthalf(n1,n3,lot,nfft,i1,zf,zw)
    INTEGER, INTENT(in)                      :: n1, n3, lot, nfft, i1
    REAL(KIND=dp), &
      DIMENSION(n1/2+1, n3/2+1), INTENT(in)  :: zf
    REAL(KIND=dp), DIMENSION(2, lot, n3/2), &
      INTENT(out)                            :: zw

    INTEGER                                  :: i01, i03i, i03r, i3, l1, l3

  zw = 0.0_dp
  i3=0
  DO l3=1,n3,2
     i3=i3+1
     i03r=ABS(l3-n3/2-1)+1
     i03i=ABS(l3-n3/2)+1
     DO l1=1,nfft
        i01=ABS(l1-1+i1-n1/2-1)+1
        zw(1,l1,i3)=zf(i01,i03r)
        zw(2,l1,i3)=zf(i01,i03i)
     END DO
  END DO

END SUBROUTINE inserthalf

! *****************************************************************************
!> \brief (Based on suitable modifications of S.Goedecker routines)
!>      Calculates the FFT of the distributed kernel
!> \param n1 logical dimension of the transform.
!> \param n2 logical dimension of the transform.
!> \param n3 logical dimension of the transform.
!> \param nd1 Dimensions of work arrays
!> \param nd2 Dimensions of work arrays
!> \param nd3 Dimensions of work arrays
!> \param nk1 ...
!> \param nk2 ...
!> \param nk3 ...
!> \param nproc number of processors used as returned by MPI_COMM_SIZE
!> \param iproc [0:nproc-1] number of processor as returned by MPI_COMM_RANK
!> \param zf Real kernel (input)
!>                   zf(i1,i2,i3)
!> \param zr Distributed Kernel FFT
!>                   zr(2,i1,i2,i3)
!> \param mpi_group ...
!> \param error ...
!> \date February 2006
!> \par Restrictions
!>      Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!>      Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!>      Copyright (C) 2002 Stefan Goedecker, CEA Grenoble
!>      This file is distributed under the terms of the
!>      GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!> \author S. Goedecker, L. Genovese
!> \note As transform lengths
!>                  most products of the prime factors 2,3,5 are allowed.
!>                   The detailed table with allowed transform lengths can
!>                   be found in subroutine CTRIG 
! *****************************************************************************
SUBROUTINE kernelfft(n1,n2,n3,nd1,nd2,nd3,nk1,nk2,nk3,nproc,iproc,zf,zr,mpi_group,error)

    INTEGER, INTENT(in)                      :: n1, n2, n3, nd1, nd2, nd3, &
                                                nk1, nk2, nk3, nproc, iproc
    REAL(KIND=dp), &
      DIMENSION(n1/2+1, n3/2+1, nd2/nproc), &
      INTENT(in)                             :: zf
    REAL(KIND=dp), &
      DIMENSION(nk1, nk2, nk3/nproc), &
      INTENT(inout)                          :: zr
    INTEGER, INTENT(in)                      :: mpi_group
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kernelfft', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: ncache_optimal = 8*1024, &
                                                timing_flag = 0

    INTEGER                                  :: i, i1, i3, i_all, i_stat, &
                                                ic1, ic2, ic3, inzee, j, j2, &
                                                J2st, j3, Jp2st, lot, lzt, &
                                                ma, mb, ncache, nfft
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: after1, after2, after3, &
                                                before1, before2, before3, &
                                                now1, now2, now3
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: twopion
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: cosinarr, trig1, trig2, trig3
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: zt, zw
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: zmpi2
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :, :)               :: zmpi1

!  include "perfdata.inc"
!work arrays for transpositions
!work arrays for MPI
!cache work array
!FFT work arrays
!Body
!check input

  CPPrecondition (nd1.GE.n1,cp_failure_level,routineP,error,failure)
  CPPrecondition (nd2.GE.n2 ,cp_failure_level,routineP,error,failure)
  CPPrecondition (nd3.GE.n3/2+1 ,cp_failure_level,routineP,error,failure)
  CPPrecondition (MOD(nd3,nproc).EQ.0 ,cp_failure_level,routineP,error,failure)
  CPPrecondition (MOD(nd2,nproc).EQ.0 ,cp_failure_level,routineP,error,failure)

  !defining work arrays dimensions
  ncache=ncache_optimal
  IF (ncache <= MAX(n1,n2,n3/2)*4) ncache=MAX(n1,n2,n3/2)*4
  lzt=n2
  IF (MOD(n2,2).EQ.0) lzt=lzt+1
  IF (MOD(n2,4).EQ.0) lzt=lzt+1

  !Allocations
  ALLOCATE(trig1(2,ctrig_length),stat=i_all)
  ALLOCATE(after1(7),stat=i_stat)
  i_all=i_all+i_stat
  ALLOCATE(now1(7),stat=i_stat)
  i_all=i_all+i_stat
  ALLOCATE(before1(7),stat=i_stat)
  i_all=i_all+i_stat
  ALLOCATE(trig2(2,ctrig_length),stat=i_stat)
  i_all=i_all+i_stat
  ALLOCATE(after2(7),stat=i_stat)
  i_all=i_all+i_stat
  ALLOCATE(now2(7),stat=i_stat)
  i_all=i_all+i_stat
  ALLOCATE(before2(7),stat=i_stat)
  i_all=i_all+i_stat
  ALLOCATE(trig3(2,ctrig_length),stat=i_stat)
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
  ALLOCATE(zmpi2(2,n1,nd2/nproc,nd3),stat=i_stat)
  i_all=i_all+i_stat
  ALLOCATE(cosinarr(2,n3/2),stat=i_stat)
  i_all=i_all+i_stat
  IF (nproc.GT.1) THEN
     ALLOCATE(zmpi1(2,n1,nd2/nproc,nd3/nproc,nproc),stat=i_stat)
     zmpi1 = 0.0_dp
  END IF
  CPPostcondition(i_all==0,cp_failure_level,routineP,error,failure)

  zmpi2 = 0.0_dp
  !calculating the FFT work arrays (beware on the HalFFT in n3 dimension)
  CALL ctrig(n3/2,trig3,after3,before3,now3,1,ic3)
  CALL ctrig(n1,trig1,after1,before1,now1,1,ic1)
  CALL ctrig(n2,trig2,after2,before2,now2,1,ic2)

  !Calculating array of phases for HalFFT decoding
  twopion=8._dp*ATAN(1._dp)/REAL(n3,KIND=dp)
  DO i3=1,n3/2
     cosinarr(1,i3)=COS(twopion*(i3-1))
     cosinarr(2,i3)=-SIN(twopion*(i3-1))
  END DO

  !transform along z axis

  lot=ncache/(2*n3)
  CPPostcondition(lot.GE.1,cp_failure_level,routineP,error,failure)

  DO j2=1,nd2/nproc
     !this condition ensures that we manage only the interesting part for the FFT
     IF (iproc*(nd2/nproc)+j2.LE.n2) THEN
        DO i1=1,n1,lot
           ma=i1
           mb=MIN(i1+(lot-1),n1)
           nfft=mb-ma+1

           !inserting real data into complex array of half lenght
           !input: I1,I3,J2,(Jp2)

           CALL inserthalf(n1,n3,lot,nfft,i1,zf(1,1,j2),zw(1,1,1))

           !performing FFT
           inzee=1
           DO i=1,ic3
              CALL fftstp(lot,nfft,n3/2,lot,n3/2,zw(1,1,inzee),zw(1,1,3-inzee), &
                   trig3,after3(i),now3(i),before3(i),1)
              inzee=3-inzee
           ENDDO
           !output: I1,i3,J2,(Jp2)

           !unpacking FFT in order to restore correct result,
           !while exchanging components
           !input: I1,i3,J2,(Jp2)
           CALL scramble_unpack(i1,j2,lot,nfft,n1,n3,nd2,nproc,nd3,zw(1,1,inzee),zmpi2,cosinarr)
           !output: I1,J2,i3,(Jp2)
        END DO
     ENDIF
  END DO

  !Interprocessor data transposition
  !input: I1,J2,j3,jp3,(Jp2)
  IF (nproc.GT.1) THEN
     !communication scheduling
     CALL mp_alltoall(zmpi2,&!2*n1*(nd2/nproc)*(nd3/nproc), &
          zmpi1,2*n1*(nd2/nproc)*(nd3/nproc), &
          mpi_group)
     ! output: I1,J2,j3,Jp2,(jp3)
  ENDIF

  DO j3=1,nd3/nproc
     !this condition ensures that we manage only the interesting part for the FFT
     IF (iproc*(nd3/nproc)+j3.LE.n3/2+1) THEN
        Jp2st=1
        J2st=1

        !transform along x axis
        lot=ncache/(4*n1)
        CPPostcondition(lot.GE.1,cp_failure_level,routineP,error,failure)

        DO j=1,n2,lot
           ma=j
           mb=MIN(j+(lot-1),n2)
           nfft=mb-ma+1

           !reverse ordering
           !input: I1,J2,j3,Jp2,(jp3)
           IF (nproc.EQ.1) THEN
              CALL mpiswitch(j3,nfft,Jp2st,J2st,lot,n1,nd2,nd3,nproc,zmpi2,zw(1,1,1))
           ELSE
              CALL mpiswitch(j3,nfft,Jp2st,J2st,lot,n1,nd2,nd3,nproc,zmpi1,zw(1,1,1))
           ENDIF
           !output: J2,Jp2,I1,j3,(jp3)

           !performing FFT
           !input: I2,I1,j3,(jp3)
           inzee=1
           DO i=1,ic1-1
              CALL fftstp(lot,nfft,n1,lot,n1,zw(1,1,inzee),zw(1,1,3-inzee), &
                   trig1,after1(i),now1(i),before1(i),1)
              inzee=3-inzee
           ENDDO
           !storing the last step into zt
           i=ic1
           CALL fftstp(lot,nfft,n1,lzt,n1,zw(1,1,inzee),zt(1,j,1), &
                trig1,after1(i),now1(i),before1(i),1)
           !output: I2,i1,j3,(jp3)
        END DO

        !transform along y axis, and taking only the first half
        lot=ncache/(4*n2)
        CPPostcondition(lot.GE.1,cp_failure_level,routineP,error,failure)

        DO j=1,nk1,lot
           ma=j
           mb=MIN(j+(lot-1),nk1)
           nfft=mb-ma+1

           !reverse ordering
           !input: I2,i1,j3,(jp3)
           CALL switch(nfft,n2,lot,n1,lzt,zt(1,1,j),zw(1,1,1))
           !output: i1,I2,j3,(jp3)

           !performing FFT
           !input: i1,I2,j3,(jp3)
           inzee=1
           DO i=1,ic2
              CALL fftstp(lot,nfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee), &
                   trig2,after2(i),now2(i),before2(i),1)
              inzee=3-inzee
           ENDDO

           CALL realcopy(lot,nfft,n2,nk1,nk2,zw(1,1,inzee),zr(j,1,j3))

        END DO
        !output: i1,i2,j3,(jp3)
     ENDIF
  END DO

  !De-allocations
  DEALLOCATE(trig1,stat=i_all)
  DEALLOCATE(after1,stat=i_stat)
  i_all=i_all+i_stat
  DEALLOCATE(now1,stat=i_stat)
  i_all=i_all+i_stat
  DEALLOCATE(before1,stat=i_stat)
  i_all=i_all+i_stat
  DEALLOCATE(trig2,stat=i_stat)
  i_all=i_all+i_stat
  DEALLOCATE(after2,stat=i_stat)
  i_all=i_all+i_stat
  DEALLOCATE(now2,stat=i_stat)
  i_all=i_all+i_stat
  DEALLOCATE(before2,stat=i_stat)
  i_all=i_all+i_stat
  DEALLOCATE(trig3,stat=i_stat)
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
  CPPostcondition(i_all==0,cp_failure_level,routineP,error,failure)

END SUBROUTINE kernelfft

! *****************************************************************************
!> \brief ...
!> \param lot ...
!> \param nfft ...
!> \param n2 ...
!> \param nk1 ...
!> \param nk2 ...
!> \param zin ...
!> \param zout ...
! *****************************************************************************
SUBROUTINE realcopy(lot,nfft,n2,nk1,nk2,zin,zout)
    INTEGER, INTENT(in)                      :: lot, nfft, n2, nk1, nk2
    REAL(KIND=dp), DIMENSION(2, lot, n2), &
      INTENT(in)                             :: zin
    REAL(KIND=dp), DIMENSION(nk1, nk2), &
      INTENT(inout)                          :: zout

    INTEGER                                  :: i, j

  DO i=1,nk2
     DO j=1,nfft
        zout(j,i)=zin(1,j,i)
     END DO
  END DO

END SUBROUTINE realcopy

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
SUBROUTINE switch(nfft,n2,lot,n1,lzt,zt,zw)
    INTEGER                                  :: nfft, n2, lot, n1, lzt
    REAL(KIND=dp)                            :: zt(2,lzt,n1), zw(2,lot,n2)

    INTEGER                                  :: i, j

  DO 200,j=1,nfft
     DO 100,i=1,n2
        zw(1,j,i)=zt(1,i,j)
        zw(2,j,i)=zt(2,i,j)
100     CONTINUE
200     CONTINUE
        RETURN
      END SUBROUTINE switch

! *****************************************************************************
!> \brief ...
!> \param j3 ...
!> \param nfft ...
!> \param Jp2st ...
!> \param J2st ...
!> \param lot ...
!> \param n1 ...
!> \param nd2 ...
!> \param nd3 ...
!> \param nproc ...
!> \param zmpi1 ...
!> \param zw ...
! *****************************************************************************
      SUBROUTINE mpiswitch(j3,nfft,Jp2st,J2st,lot,n1,nd2,nd3,nproc,zmpi1,zw)
    INTEGER                                  :: j3, nfft, Jp2st, J2st, lot, &
                                                n1, nd2, nd3, nproc
    REAL(KIND=dp) :: zmpi1(2,n1,nd2/nproc,nd3/nproc,nproc), zw(2,lot,n1)

    INTEGER                                  :: I1, J2, JP2, mfft

        mfft=0
        DO 300,Jp2=Jp2st,nproc
           DO 200,J2=J2st,nd2/nproc
              mfft=mfft+1
              IF (mfft.GT.nfft) THEN
                 Jp2st=Jp2
                 J2st=J2
                 RETURN
              ENDIF
        DO 100,I1=1,n1
           zw(1,mfft,I1)=zmpi1(1,I1,J2,j3,Jp2)
           zw(2,mfft,I1)=zmpi1(2,I1,J2,j3,Jp2)
100        CONTINUE
200        CONTINUE
           J2st=1
300        CONTINUE
         END SUBROUTINE mpiswitch

! *****************************************************************************
!> \brief ...
!> \param nterms ...
!> \param p ...
!> \param w ...
!> \param urange ...
!> \param drange ...
!> \param acc ...
! *****************************************************************************
         SUBROUTINE gequad(nterms,p,w,urange,drange,acc)
!
    INTEGER                                  :: nterms
    REAL(KIND=dp)                            :: p(*), w(*), urange, drange, &
                                                acc

!
!
!       range [10^(-9),1] and accuracy ~10^(-8);
!
!

           p(1)=4.96142640560223544e19_dp
           p(2)=1.37454269147978052e19_dp
           p(3)=7.58610013441204679e18_dp
           p(4)=4.42040691347806996e18_dp
           p(5)=2.61986077948367892e18_dp
           p(6)=1.56320138155496681e18_dp
           p(7)=9.35645215863028402e17_dp
           p(8)=5.60962910452691703e17_dp
           p(9)=3.3666225119686761e17_dp
           p(10)=2.0218253197947866e17_dp
           p(11)=1.21477756091902017e17_dp
           p(12)=7.3012982513608503e16_dp
           p(13)=4.38951893556421099e16_dp
           p(14)=2.63949482512262325e16_dp
           p(15)=1.58742054072786174e16_dp
           p(16)=9.54806587737665531e15_dp
           p(17)=5.74353712364571709e15_dp
           p(18)=3.455214877389445e15_dp
           p(19)=2.07871658520326804e15_dp
           p(20)=1.25064667315629928e15_dp
           p(21)=7.52469429541933745e14_dp
           p(22)=4.5274603337253175e14_dp
           p(23)=2.72414006900059548e14_dp
           p(24)=1.63912168349216752e14_dp
           p(25)=9.86275802590865738e13_dp
           p(26)=5.93457701624974985e13_dp
           p(27)=3.5709554322296296e13_dp
           p(28)=2.14872890367310454e13_dp
           p(29)=1.29294719957726902e13_dp
           p(30)=7.78003375426361016e12_dp
           p(31)=4.68148199759876704e12_dp
           p(32)=2.8169955024829868e12_dp
           p(33)=1.69507790481958464e12_dp
           p(34)=1.01998486064607581e12_dp
           p(35)=6.13759486539856459e11_dp
           p(36)=3.69320183828682544e11_dp
           p(37)=2.22232783898905102e11_dp
           p(38)=1.33725247623668682e11_dp
           p(39)=8.0467192739036288e10_dp
           p(40)=4.84199582415144143e10_dp
           p(41)=2.91360091170559564e10_dp
           p(42)=1.75321747475309216e10_dp
           p(43)=1.0549735552210995e10_dp
           p(44)=6.34815321079006586e9_dp
           p(45)=3.81991113733594231e9_dp
           p(46)=2.29857747533101109e9_dp
           p(47)=1.38313653595483694e9_dp
           p(48)=8.32282908580025358e8_dp
           p(49)=5.00814519374587467e8_dp
           p(50)=3.01358090773319025e8_dp
           p(51)=1.81337994217503535e8_dp
           p(52)=1.09117589961086823e8_dp
           p(53)=6.56599771718640323e7_dp
           p(54)=3.95099693638497164e7_dp
           p(55)=2.37745694710665991e7_dp
           p(56)=1.43060135285912813e7_dp
           p(57)=8.60844290313506695e6_dp
           p(58)=5.18000974075383424e6_dp
           p(59)=3.116998193057466e6_dp
           p(60)=1.87560993870024029e6_dp
           p(61)=1.12862197183979562e6_dp
           p(62)=679132.441326077231_dp
           p(63)=408658.421279877969_dp
           p(64)=245904.473450669789_dp
           p(65)=147969.568088321005_dp
           p(66)=89038.612357311147_dp
           p(67)=53577.7362552358895_dp
           p(68)=32239.6513926914668_dp
           p(69)=19399.7580852362791_dp
           p(70)=11673.5323603058634_dp
           p(71)=7024.38438577707758_dp
           p(72)=4226.82479307685999_dp
           p(73)=2543.43254175354295_dp
           p(74)=1530.47486269122675_dp
           p(75)=920.941785160749482_dp
           p(76)=554.163803906291646_dp
           p(77)=333.46029740785694_dp
           p(78)=200.6550575335041_dp
           p(79)=120.741366914147284_dp
           p(80)=72.6544243200329916_dp
           p(81)=43.7187810415471025_dp
           p(82)=26.3071631447061043_dp
           p(83)=15.8299486353816329_dp
           p(84)=9.52493152341244004_dp
           p(85)=5.72200417067776041_dp
           p(86)=3.36242234070940928_dp
           p(87)=1.75371394604499472_dp
           p(88)=0.64705932650658966_dp
           p(89)=0.072765905943708247_dp
           !
           w(1)=47.67445484528304247e10_dp
           w(2)=11.37485774750442175e9_dp
           w(3)=78.64340976880190239e8_dp
           w(4)=46.27335788759590498e8_dp
           w(5)=24.7380464827152951e8_dp
           w(6)=13.62904116438987719e8_dp
           w(7)=92.79560029045882433e8_dp
           w(8)=52.15931216254660251e8_dp
           w(9)=31.67018011061666244e8_dp
           w(10)=1.29291036801493046e8_dp
           w(11)=1.00139319988015862e8_dp
           w(12)=7.75892350510188341e7_dp
           w(13)=6.01333567950731271e7_dp
           w(14)=4.66141178654796875e7_dp
           w(15)=3.61398903394911448e7_dp
           w(16)=2.80225846672956389e7_dp
           w(17)=2.1730509180930247e7_dp
           w(18)=1.68524482625876965e7_dp
           w(19)=1.30701489345870338e7_dp
           w(20)=1.01371784832269282e7_dp
           w(21)=7.86264116300379329e6_dp
           w(22)=6.09861667912273717e6_dp
           w(23)=4.73045784039455683e6_dp
           w(24)=3.66928949951594161e6_dp
           w(25)=2.8462050836230259e6_dp
           w(26)=2.20777394798527011e6_dp
           w(27)=1.71256191589205524e6_dp
           w(28)=1.32843556197737076e6_dp
           w(29)=1.0304731275955989e6_dp
           w(30)=799345.206572271448_dp
           w(31)=620059.354143595343_dp
           w(32)=480986.704107449333_dp
           w(33)=373107.167700228515_dp
           w(34)=289424.08337412132_dp
           w(35)=224510.248231581788_dp
           w(36)=174155.825690028966_dp
           w(37)=135095.256919654065_dp
           w(38)=104795.442776800312_dp
           w(39)=81291.4458222430418_dp
           w(40)=63059.0493649328682_dp
           w(41)=48915.9040455329689_dp
           w(42)=37944.8484018048756_dp
           w(43)=29434.4290473253969_dp
           w(44)=22832.7622054490044_dp
           w(45)=17711.743950151233_dp
           w(46)=13739.287867104177_dp
           w(47)=10657.7895710752585_dp
           w(48)=8267.42141053961834_dp
           w(49)=6413.17397520136448_dp
           w(50)=4974.80402838654277_dp
           w(51)=3859.03698188553047_dp
           w(52)=2993.51824493299154_dp
           w(53)=2322.1211966811754_dp
           w(54)=1801.30750964719641_dp
           w(55)=1397.30379659817038_dp
           w(56)=1083.91149143250697_dp
           w(57)=840.807939169209188_dp
           w(58)=652.228524366749422_dp
           w(59)=505.944376983506128_dp
           w(60)=392.469362317941064_dp
           w(61)=304.444930257324312_dp
           w(62)=236.162932842453601_dp
           w(63)=183.195466078603525_dp
           w(64)=142.107732186551471_dp
           w(65)=110.23530215723992_dp
           w(66)=85.5113346705382257_dp
           w(67)=66.3325469806696621_dp
           w(68)=51.4552463353841373_dp
           w(69)=39.9146798429449273_dp
           w(70)=30.9624728409162095_dp
           w(71)=24.018098812215013_dp
           w(72)=18.6312338024296588_dp
           w(73)=14.4525541233150501_dp
           w(74)=11.2110836519105938_dp
           w(75)=8.69662175848497178_dp
           w(76)=6.74611236165731961_dp
           w(77)=5.23307018057529994_dp
           w(78)=4.05937850501539556_dp
           w(79)=3.14892659076635714_dp
           w(80)=2.44267408211071604_dp
           w(81)=1.89482240522855261_dp
           w(82)=1.46984505907050079_dp
           w(83)=1.14019261330527007_dp
           w(84)=0.884791217422925293_dp
           w(85)=0.692686387080616483_dp
           w(86)=0.585244576897023282_dp
           w(87)=0.576182522545327589_dp
           w(88)=0.596688817388997178_dp
           w(89)=0.607879901151108771_dp
           !
           !
           urange = 1._dp
           drange=1e-08_dp
           acc   =1e-08_dp
           !
           RETURN
         END SUBROUTINE gequad

END MODULE ps_wavelet_kernel
