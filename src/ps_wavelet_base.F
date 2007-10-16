!!****** cp2k/ ps_wavelet_base
!!
!!   NAME
!!     ps_wavelet_base
!!
!!   FUNCTION
!!     Creates the wavelet kernel for the wavelet based poisson solver.
!!
!!   AUTHOR
!!     Florian Schiffmann (09.2007,fschiff)
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************


MODULE ps_wavelet_base
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_alltoall
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE ps_wavelet_fft3d,                ONLY: fourier_dim,&
                                             ctrig,&
                                             fftstp


#include "cp_common_uses.h"

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ps_wavelet_util'


  PUBLIC :: scramble_unpack




CONTAINS


subroutine P_PoissonSolver(n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,nproc,iproc,zf&
             ,scal,hx,hy,hz,mpi_group,error)
  implicit none
  !to be preprocessed
!  include 'mpif.h'
!  include "perfdata.inc"
integer, parameter :: ncache_optimal=8*1024
integer, parameter :: timing_flag=0
integer :: count_time1,count_time2,count_rate,count_max,number_time,index_time
real(kind=dp) :: time_b0,time_b1,serial_time,parallel_time
  !Arguments
  INTEGER, INTENT(in) :: n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,nproc,iproc,mpi_group
  TYPE(cp_error_type), INTENT(inout)       :: error
  real(KIND=dp), intent(in) :: scal,hx,hy,hz
  real(KIND=dp), dimension(md1,md3,md2/nproc), intent(inout) :: zf
  !Local variables
  integer :: ncache,lzt,lot,ma,mb,nfft,ic1,ic2,ic3,Jp2stb,J2stb,Jp2stf,J2stf
  integer :: j1,j2,j3,i1,i2,i3,i,j,inzee,ierr,i_all,i_stat
  !work arrays for transpositions
  real(KIND=dp), dimension(:,:,:), allocatable :: zt
  !work arrays for MPI
  real(KIND=dp), dimension(:,:,:,:,:), allocatable :: zmpi1
  real(KIND=dp), dimension(:,:,:,:), allocatable :: zmpi2
  !cache work array
  real(KIND=dp), dimension(:,:,:), allocatable :: zw
  !FFT work arrays
  real(KIND=dp), dimension(:,:), allocatable :: btrig1,btrig2,btrig3, &
       ftrig1,ftrig2,ftrig3
  integer, dimension(:), allocatable :: after1,now1,before1, & 
       after2,now2,before2,after3,now3,before3
  
  !Body

!  call timing(iproc,'PSolv_comput  ','ON')
  ! check input
!!$  if (mod(n1,2).ne.0) stop 'Parallel convolution:ERROR:n1' !this can be avoided
!!$  if (mod(n2,2).ne.0) stop 'Parallel convolution:ERROR:n2' !this can be avoided
!!$  if (mod(n3,2).ne.0) stop 'Parallel convolution:ERROR:n3' !this can be avoided
  if (nd1.lt.n1/2+1) stop 'Parallel convolution:ERROR:nd1' 
  if (nd2.lt.n2/2+1) stop 'Parallel convolution:ERROR:nd2' 
  if (nd3.lt.n3/2+1) stop 'Parallel convolution:ERROR:nd3' 
  if (md1.lt.n1) stop 'Parallel convolution:ERROR:md1'
  if (md2.lt.n2) stop 'Parallel convolution:ERROR:md2'
  if (md3.lt.n3) stop 'Parallel convolution:ERROR:md3'
  if (mod(nd3,nproc).ne.0) stop 'Parallel convolution:ERROR:nd3'
  if (mod(md2,nproc).ne.0) stop 'Parallel convolution:ERROR:md2'
  
  !defining work arrays dimensions
  ncache=ncache_optimal
  if (ncache <= max(n1,n2,n3)*4) ncache=max(n1,n2,n3)*4

!  if (timing_flag == 1 .and. iproc ==0) print *,'parallel ncache=',ncache

  lzt=n2
  if (mod(n2,2).eq.0) lzt=lzt+1
  if (mod(n2,4).eq.0) lzt=lzt+1 !maybe this is useless
  
  !Allocations
  allocate(btrig1(2,8192),stat=i_all)
  allocate(ftrig1(2,8192),stat=i_stat)
  i_all=i_all+i_stat
  allocate(after1(7),stat=i_stat)
  i_all=i_all+i_stat
  allocate(now1(7),stat=i_stat)
  i_all=i_all+i_stat
  allocate(before1(7),stat=i_stat)
  i_all=i_all+i_stat
  allocate(btrig2(2,8192),stat=i_stat)
  i_all=i_all+i_stat
  allocate(ftrig2(2,8192),stat=i_stat)
  i_all=i_all+i_stat
  allocate(after2(7),stat=i_stat)
  i_all=i_all+i_stat
  allocate(now2(7),stat=i_stat)
  i_all=i_all+i_stat
  allocate(before2(7),stat=i_stat)
  i_all=i_all+i_stat
  allocate(btrig3(2,8192),stat=i_stat)
  i_all=i_all+i_stat
  allocate(ftrig3(2,8192),stat=i_stat)
  i_all=i_all+i_stat
  allocate(after3(7),stat=i_stat)
  i_all=i_all+i_stat
  allocate(now3(7),stat=i_stat)
  i_all=i_all+i_stat
  allocate(before3(7),stat=i_stat)
  i_all=i_all+i_stat
  allocate(zw(2,ncache/4,2),stat=i_stat)
  i_all=i_all+i_stat
  allocate(zt(2,lzt,n1),stat=i_stat)
  i_all=i_all+i_stat
  allocate(zmpi2(2,n1,md2/nproc,nd3),stat=i_stat)
  i_all=i_all+i_stat
  if (nproc.gt.1) allocate(zmpi1(2,n1,md2/nproc,nd3/nproc,nproc),stat=i_stat)
  if (i_all+i_stat /= 0) then
     write(*,*)' p_poissonsolver: problem of memory allocation'
     stop
  end if


  !calculating the FFT work arrays (beware on the HalFFT in n3 dimension)
  call ctrig(n3,btrig3,after3,before3,now3,1,ic3)
  call ctrig(n1,btrig1,after1,before1,now1,1,ic1)
  call ctrig(n2,btrig2,after2,before2,now2,1,ic2)
  do  j=1,n1
     ftrig1(1,j)= btrig1(1,j)
     ftrig1(2,j)=-btrig1(2,j)
  enddo
  do  j=1,n2
     ftrig2(1,j)= btrig2(1,j)
     ftrig2(2,j)=-btrig2(2,j)
  enddo
  do  j=1,n3
     ftrig3(1,j)= btrig3(1,j)
     ftrig3(2,j)=-btrig3(2,j)
  enddo

  ! transform along z axis
  lot=ncache/(4*n3)
  if (lot.lt.1) then  
     write(6,*) & 
          'convolxc_off:ncache has to be enlarged to be able to hold at' // &  
          'least one 1-d FFT of this size even though this will' // & 
          'reduce the performance for shorter transform lengths'
     stop
  endif
  do j2=1,md2/nproc
     !this condition ensures that we manage only the interesting part for the FFT
     if (iproc*(md2/nproc)+j2.le.n2) then
        do i1=1,n1,lot
           ma=i1
           mb=min(i1+(lot-1),n1)
           nfft=mb-ma+1
           !inserting real data into complex array of half lenght
           call P_fill_upcorn(md1,md3,lot,nfft,n3,zf(i1,1,j2),zw(1,1,1))

           !performing FFT
           !input: I1,I3,J2,(Jp2)
           inzee=1
           do i=1,ic3
              call fftstp(lot,nfft,n3,lot,n3,zw(1,1,inzee),zw(1,1,3-inzee), &
                   btrig3,after3(i),now3(i),before3(i),1)
              inzee=3-inzee
           enddo

           !output: I1,i3,J2,(Jp2)
           !exchanging components
           !input: I1,i3,J2,(Jp2)
           call scramble_P(i1,j2,lot,nfft,n1,n3,md2,nproc,nd3,zw(1,1,inzee),zmpi2)
           !output: I1,J2,i3,(Jp2)
        end do
     end if
  end do

  !Interprocessor data transposition
  !input: I1,J2,j3,jp3,(Jp2)
  if (nproc.gt.1) then
!     call timing(iproc,'PSolv_comput  ','OF')

!     call timing(iproc,'PSolv_commun  ','ON')

     !communication scheduling
     call mp_alltoall(zmpi2,& 
          zmpi1,2*n1*(md2/nproc)*(nd3/nproc), &
          mpi_group)

!     call timing(iproc,'PSolv_commun  ','OF')

!     call timing(iproc,'PSolv_comput  ','ON')
  endif
  !output: I1,J2,j3,Jp2,(jp3)

  !now each process perform complete convolution of its planes
  do j3=1,nd3/nproc
     !this condition ensures that we manage only the interesting part for the FFT
     if (iproc*(nd3/nproc)+j3.le.n3/2+1) then
	Jp2stb=1
	J2stb=1
	Jp2stf=1
	J2stf=1
        
        ! transform along x axis
        lot=ncache/(4*n1)
        if (lot.lt.1) then  
           write(6,*) & 
                'convolxc_off:ncache has to be enlarged to be able to hold at' // &  
                'least one 1-d FFT of this size even though this will' // & 
                'reduce the performance for shorter transform lengths'
           stop
        endif
        
        do j=1,n2,lot
           ma=j
           mb=min(j+(lot-1),n2)
           nfft=mb-ma+1

           !reverse index ordering, leaving the planes to be transformed at the end
           !input: I1,J2,j3,Jp2,(jp3)
           if (nproc.eq.1) then
              call P_mpiswitch_upcorn(j3,nfft,Jp2stb,J2stb,lot,n1,md2,nd3,nproc,zmpi2,zw(1,1,1))
           else
              call P_mpiswitch_upcorn(j3,nfft,Jp2stb,J2stb,lot,n1,md2,nd3,nproc,zmpi1,zw(1,1,1))
           endif
           !output: J2,Jp2,I1,j3,(jp3)

           !performing FFT
           !input: I2,I1,j3,(jp3)
           inzee=1
           do i=1,ic1-1
              call fftstp(lot,nfft,n1,lot,n1,zw(1,1,inzee),zw(1,1,3-inzee), &
                   btrig1,after1(i),now1(i),before1(i),1)
              inzee=3-inzee
           enddo

           !storing the last step into zt array
           i=ic1
           call fftstp(lot,nfft,n1,lzt,n1,zw(1,1,inzee),zt(1,j,1), & 
                btrig1,after1(i),now1(i),before1(i),1)           
           !output: I2,i1,j3,(jp3)
        end do

        !transform along y axis
        lot=ncache/(4*n2)
        if (lot.lt.1) then  
           write(6,*) & 
                'convolxc_off:ncache has to be enlarged to be able to hold at' // &  
                'least one 1-d FFT of this size even though this will' // & 
                'reduce the performance for shorter transform lengths'
           stop
        endif

        do j=1,n1,lot
           ma=j
           mb=min(j+(lot-1),n1)
           nfft=mb-ma+1

           !reverse ordering 
           !input: I2,i1,j3,(jp3)
           call P_switch_upcorn(nfft,n2,lot,n1,lzt,zt(1,1,j),zw(1,1,1))
           !output: i1,I2,j3,(jp3)
           
           !performing FFT
           !input: i1,I2,j3,(jp3)
           inzee=1
           do i=1,ic2
              call fftstp(lot,nfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee), &
                   btrig2,after2(i),now2(i),before2(i),1)
              inzee=3-inzee
           enddo
           !output: i1,i2,j3,(jp3)


           !Multiply with kernel in fourier space
           i3=iproc*(nd3/nproc)+j3
           call P_multkernel(n1,n2,n3,lot,nfft,j,i3,zw(1,1,inzee),hx,hy,hz)
           
           !TRANSFORM BACK IN REAL SPACE
           
           !transform along y axis
           !input: i1,i2,j3,(jp3)
           do i=1,ic2
              call fftstp(lot,nfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee), &
                   ftrig2,after2(i),now2(i),before2(i),-1)
              inzee=3-inzee
           enddo

           !reverse ordering
           !input: i1,I2,j3,(jp3)
           call P_unswitch_downcorn(nfft,n2,lot,n1,lzt,zw(1,1,inzee),zt(1,1,j))
           !output: I2,i1,j3,(jp3)
        end do
        
        !transform along x axis
        !input: I2,i1,j3,(jp3)
        lot=ncache/(4*n1)
        do j=1,n2,lot
           ma=j
           mb=min(j+(lot-1),n2)
           nfft=mb-ma+1

           !performing FFT
           i=1
           call fftstp(lzt,nfft,n1,lot,n1,zt(1,j,1),zw(1,1,1), &
                ftrig1,after1(i),now1(i),before1(i),-1)
           
           inzee=1
           do i=2,ic1
              call fftstp(lot,nfft,n1,lot,n1,zw(1,1,inzee),zw(1,1,3-inzee), &
                   ftrig1,after1(i),now1(i),before1(i),-1)
              inzee=3-inzee
           enddo
           !output: I2,I1,j3,(jp3)

           !reverse ordering
           !input: J2,Jp2,I1,j3,(jp3)
           if (nproc.eq.1) then
              call P_unmpiswitch_downcorn(j3,nfft,Jp2stf,J2stf,lot,n1,md2,nd3,nproc,zw(1,1,inzee),zmpi2)
           else
              call P_unmpiswitch_downcorn(j3,nfft,Jp2stf,J2stf,lot,n1,md2,nd3,nproc,zw(1,1,inzee),zmpi1)
           endif
           ! output: I1,J2,j3,Jp2,(jp3)
        end do
     endif
  end do

  !Interprocessor data transposition
  !input: I1,J2,j3,Jp2,(jp3)
  if (nproc.gt.1) then
!!     call timing(iproc,'PSolv_comput  ','OF')

!     call timing(iproc,'PSolv_commun  ','ON')
     !communication scheduling

     call mp_alltoall(zmpi1,& 
          zmpi2,2*n1*(md2/nproc)*(nd3/nproc), &
          mpi_group)
!     call timing(iproc,'PSolv_commun  ','OF')

!     call timing(iproc,'PSolv_comput  ','ON')
  endif
     !output: I1,J2,j3,jp3,(Jp2)
  !transform along z axis
  !input: I1,J2,i3,(Jp2)
  lot=ncache/(4*n3)
  do j2=1,md2/nproc
     !this condition ensures that we manage only the interesting part for the FFT
     if (iproc*(md2/nproc)+j2.le.n2) then
        do i1=1,n1,lot
           ma=i1
           mb=min(i1+(lot-1),n1)
           nfft=mb-ma+1

           !reverse ordering
           !input: I1,J2,i3,(Jp2)
           call unscramble_P(i1,j2,lot,nfft,n1,n3,md2,nproc,nd3,zmpi2,zw(1,1,1))
           !output: I1,i3,J2,(Jp2)

           !performing FFT
           !input: I1,i3,J2,(Jp2)           
           inzee=1
           do i=1,ic3
              call fftstp(lot,nfft,n3,lot,n3,zw(1,1,inzee),zw(1,1,3-inzee), &
                   ftrig3,after3(i),now3(i),before3(i),-1)
              inzee=3-inzee
           enddo
           !output: I1,I3,J2,(Jp2)

           !rebuild the output array
           call P_unfill_downcorn(md1,md3,lot,nfft,n3,zw(1,1,inzee),zf(i1,1,j2),scal)

        end do
     endif
  end do

  !De-allocations  
  deallocate(btrig1,stat=i_all)
  deallocate(ftrig1,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(after1,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(now1,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(before1,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(btrig2,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(ftrig2,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(after2,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(now2,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(before2,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(btrig3,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(ftrig3,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(after3,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(now3,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(before3,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(zmpi2,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(zw,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(zt,stat=i_stat)
  i_all=i_all+i_stat
  if (nproc.gt.1) deallocate(zmpi1,stat=i_stat)
  if (i_all+i_stat /= 0) then
     write(*,*)' p_poissonsolver: problem of memory deallocation'
     stop
  end if
!  call timing(iproc,'PSolv_comput  ','OF')
end subroutine P_PoissonSolver

subroutine P_mpiswitch_upcorn(j3,nfft,Jp2stb,J2stb,lot,n1,md2,nd3,nproc,zmpi1,zw)
  implicit none
  integer, intent(in) :: j3,nfft,lot,n1,md2,nd3,nproc
  integer, intent(inout) :: Jp2stb,J2stb
  real(KIND=dp), dimension(2,n1,md2/nproc,nd3/nproc,nproc), intent(in) :: zmpi1
  real(KIND=dp), dimension(2,lot,n1), intent(inout) :: zw
  !local variables
  integer :: I1,J2,Jp2,mfft

  mfft=0
  do Jp2=Jp2stb,nproc
     do J2=J2stb,md2/nproc
	mfft=mfft+1
	if (mfft.gt.nfft) then
           Jp2stb=Jp2
           J2stb=J2
           return
	end if
        do I1=1,n1
           zw(1,mfft,I1)=zmpi1(1,I1,J2,j3,Jp2)
           zw(2,mfft,I1)=zmpi1(2,I1,J2,j3,Jp2)
        end do
     end do
     J2stb=1
  end do

end subroutine P_mpiswitch_upcorn

subroutine P_switch_upcorn(nfft,n2,lot,n1,lzt,zt,zw)
  implicit none
  integer, intent(in) :: nfft,n2,lot,n1,lzt
  real(KIND=dp), dimension(2,lzt,n1), intent(in) :: zt
  real(KIND=dp), dimension(2,lot,n2), intent(inout) :: zw
  !local variables
  integer :: i,j

  do j=1,nfft
     do i=1,n2
	zw(1,j,i)=zt(1,i,j)
	zw(2,j,i)=zt(2,i,j)
     end do
  end do

end subroutine P_switch_upcorn

subroutine P_unswitch_downcorn(nfft,n2,lot,n1,lzt,zw,zt)
  implicit none
  integer, intent(in) :: nfft,n2,lot,n1,lzt
  real(KIND=dp), dimension(2,lot,n2), intent(in) :: zw
  real(KIND=dp), dimension(2,lzt,n1), intent(inout) :: zt
  !local variables
  integer :: i,j

  do j=1,nfft
     do i=1,n2
	zt(1,i,j)=zw(1,j,i)
	zt(2,i,j)=zw(2,j,i)
     end do
  end do

end subroutine P_unswitch_downcorn

subroutine P_unmpiswitch_downcorn(j3,nfft,Jp2stf,J2stf,lot,n1,md2,nd3,nproc,zw,zmpi1)
  implicit none
  integer, intent(in) :: j3,nfft,lot,n1,md2,nd3,nproc
  integer, intent(inout) :: Jp2stf,J2stf
  real(KIND=dp), dimension(2,lot,n1), intent(in) :: zw
  real(KIND=dp), dimension(2,n1,md2/nproc,nd3/nproc,nproc), intent(inout) :: zmpi1
  !local variables
  integer :: I1,J2,Jp2,mfft

  mfft=0
  do Jp2=Jp2stf,nproc
     do J2=J2stf,md2/nproc
	mfft=mfft+1
	if (mfft.gt.nfft) then
           Jp2stf=Jp2
           J2stf=J2
           return
	end if
        do I1=1,n1
           zmpi1(1,I1,J2,j3,Jp2)=zw(1,mfft,I1)
           zmpi1(2,I1,J2,j3,Jp2)=zw(2,mfft,I1)
        end do
     end do
     J2stf=1
  end do
end subroutine P_unmpiswitch_downcorn


!!****h* BigDFT/P_fill_downcorn
!! NAME
!!   P_fill_downcorn
!!
!! FUNCTION
!!     (Based on suitable modifications of S.Goedecker routines)
!!     Restore data into output array
!!
!! SYNOPSIS
!!     zf:          Original distributed density as well as
!!                  Distributed solution of the poisson equation (inout)
!!     zw:          FFT work array
!!     n3:          (twice the) dimension of the last FFTtransform.
!!     md1,md3:     Dimensions of the undistributed part of the real grid
!!     nfft:        number of planes
!!     scal:        Needed to achieve unitarity and correct dimensions
!!
!! WARNING
!!     Assuming that high frequencies are in the corners 
!!     and that n3 is multiple of 4   
!!
!! RESTRICTIONS on USAGE
!!     Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!!     Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!!     Copyright (C) 2002 Stefan Goedecker, CEA Grenoble
!!     This file is distributed under the terms of the
!!     GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!!
!! AUTHORS
!!    S. Goedecker, L. Genovese
!!
!! CREATION DATE
!!     February 2006
!!
!! SOURCE
!!
subroutine P_unfill_downcorn(md1,md3,lot,nfft,n3,zw,zf,scal)
  implicit none
  !Arguments
  integer, intent(in) :: md1,md3,lot,nfft,n3
  real(KIND=dp), dimension(2,lot,n3), intent(in) :: zw
  real(KIND=dp), dimension(md1,md3),intent(inout) :: zf
  real(KIND=dp), intent(in) :: scal
  !Local variables
  integer :: i3,i1
  real(KIND=dp) :: pot1

  do i3=1,n3
     do i1=1,nfft
        pot1 = scal*zw(1,i1,i3)
	zf(i1,i3)= pot1 
     end do
  end do

end subroutine P_unfill_downcorn

subroutine P_fill_upcorn(md1,md3,lot,nfft,n3,zf,zw)
  implicit none
  integer, intent(in) :: md1,md3,lot,nfft,n3
  real(KIND=dp), dimension(md1,md3), intent(in) :: zf
  real(KIND=dp), dimension(2,lot,n3), intent(inout) :: zw
  !local variables
  integer :: i1,i3

  do i3=1,n3
     do i1=1,nfft
	zw(1,i1,i3)=zf(i1,i3)
	zw(2,i1,i3)=0._dp
     end do
  end do

end subroutine P_fill_upcorn


!!****h* BigDFT/scramble_P
!! NAME
!!   scramble_P
!!
!! FUNCTION
!!     (Based on suitable modifications of S.Goedecker routines)
!!     Assign the correct planes to the work array zmpi2
!!     in order to prepare for interprocessor data transposition.
!!
!! SYNOPSIS
!!     zmpi2:          Work array for multiprocessor manipulation (output)
!!     zw:             Work array (input)
!!     n1,n3:          logical dimension of the FFT transform, reference for work arrays
!!     md2,nd3:        Dimensions of real grid and of the kernel, respectively
!!     i1,j2,lot,nfft: Starting points of the plane and number of remaining lines
!!
!! RESTRICTIONS on USAGE
!!     Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!!     Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!!     Copyright (C) 2002 Stefan Goedecker, CEA Grenoble
!!     This file is distributed under the terms of the
!!     GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!!
!! AUTHORS
!!    S. Goedecker, L. Genovese
!!
!! CREATION DATE
!!     February 2006
!!
!! SOURCE
!!
subroutine scramble_P(i1,j2,lot,nfft,n1,n3,md2,nproc,nd3,zw,zmpi2)
  implicit none
  !Arguments
  integer, intent(in) :: i1,j2,lot,nfft,n1,n3,md2,nproc,nd3
  real(KIND=dp), dimension(2,lot,n3), intent(in) :: zw
  real(KIND=dp), dimension(2,n1,md2/nproc,nd3), intent(inout) :: zmpi2
  !Local variables
  integer :: i3,i

  do i3=1,n3/2+1
     do i=0,nfft-1
        zmpi2(1,i1+i,j2,i3)=zw(1,i+1,i3)
        zmpi2(2,i1+i,j2,i3)=zw(2,i+1,i3)
     end do
  end do

end subroutine scramble_P

 
!!****h* BigDFT/unscramble_P
!! NAME
!!   unscramble_P
!!
!! FUNCTION
!!     (Based on suitable modifications of S.Goedecker routines)
!!     Insert the correct planes of the work array zmpi2
!!     in order to prepare for backward FFT transform
!!
!! SYNOPSIS
!!     zmpi2:          Work array for multiprocessor manipulation (input)
!!     zw:             Work array (output)
!!     n1,n3:          logical dimension of the FFT transform, reference for work arrays
!!     md2,nd3:        Dimensions of real grid and of the kernel, respectively
!!     i1,j2,lot,nfft: Starting points of the plane and number of remaining lines
!!
!! RESTRICTIONS on USAGE
!!     Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!!     Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!!     Copyright (C) 2002 Stefan Goedecker, CEA Grenoble
!!     This file is distributed under the terms of the
!!     GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!!
!! AUTHORS
!!    S. Goedecker, L. Genovese
!!
!! CREATION DATE
!!     February 2006
!!
!! SOURCE
!!
subroutine unscramble_P(i1,j2,lot,nfft,n1,n3,md2,nproc,nd3,zmpi2,zw)
  implicit none
  !Arguments
  integer, intent(in) :: i1,j2,lot,nfft,n1,n3,md2,nproc,nd3
  real(KIND=dp), dimension(2,lot,n3), intent(inout) :: zw
  real(KIND=dp), dimension(2,n1,md2/nproc,nd3), intent(in) :: zmpi2
  !Local variables
  integer :: i3,i,j3
  
  i3=1
  do i=0,nfft-1
     zw(1,i+1,i3)=zmpi2(1,i1+i,j2,i3)
     zw(2,i+1,i3)=zmpi2(2,i1+i,j2,i3)
  end do

  do i3=2,n3/2+1
     j3=n3+2-i3
     do i=0,nfft-1
        zw(1,i+1,j3)= zmpi2(1,i1+i,j2,i3)
        zw(2,i+1,j3)=-zmpi2(2,i1+i,j2,i3)
        zw(1,i+1,i3)= zmpi2(1,i1+i,j2,i3)
        zw(2,i+1,i3)= zmpi2(2,i1+i,j2,i3)
     end do
  end do

end subroutine unscramble_P

!!****h* BigDFT/P_multkernel
!! NAME
!!   P_multkernel
!!
!! FUNCTION
!!     (Based on suitable modifications of S.Goedecker routines)
!!     Multiply with the kernel taking into account its symmetry
!!     Conceived to be used into convolution loops
!!
!! SYNOPSIS
!!     zw:         Work array (input/output)
!!     n1,n2:      logical dimension of the FFT transform, reference for zw
!!     nd1,nd2:    Dimensions of POT
!!     jS,j3,nfft: starting point of the plane and number of remaining lines
!!
!! RESTRICTIONS on USAGE
!!     Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!!     Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!!     Copyright (C) 2002 Stefan Goedecker, CEA Grenoble
!!     This file is distributed under the terms of the
!!      GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!!
!! AUTHORS
!!    S. Goedecker, L. Genovese
!!
!! CREATION DATE
!!     February 2006
!!
!! SOURCE
!!
subroutine P_multkernel(n1,n2,n3,lot,nfft,jS,i3,zw,hx,hy,hz)
  implicit none
  !Argments
  integer, intent(in) :: n1,n2,n3,lot,nfft,jS,i3
  real(KIND=dp), intent(in) :: hx,hy,hz
  real(KIND=dp), dimension(2,lot,n2), intent(inout) :: zw
  !Local variables
  integer :: i1,j1,i2,j2,j3
  real(KIND=dp) :: ker,pi,fourpi2,mu3,p1,p2

  pi=4._dp*datan(1._dp)
  fourpi2=4._dp*pi**2
  j3=i3!n3/2+1-abs(n3/2+2-i3)
  mu3=real(j3-1,KIND=dp)/real(n3,KIND=dp)
  mu3=(mu3/hy)**2 !beware of the exchanged dimension
  !Body
  !generic case
  do i2=1,n2
     do i1=1,nfft
        j1=i1+jS-1
        j1=j1-(j1/(n1/2+2))*n1 !n1/2+1-abs(n1/2+2-jS-i1)
        j2=i2-(i2/(n2/2+2))*n2 !n2/2+1-abs(n2/2+1-i2)
        p1=real(j1-1,KIND=dp)/real(n1,KIND=dp)
        p2=real(j2-1,KIND=dp)/real(n2,KIND=dp)
        ker=-fourpi2*((p1/hx)**2+(p2/hz)**2+mu3)!beware of the exchanged dimension
        if (ker/=0._dp) ker=1._dp/ker
        zw(1,i1,i2)=zw(1,i1,i2)*ker
        zw(2,i1,i2)=zw(2,i1,i2)*ker
     end do
  end do

end subroutine P_multkernel

!!****h* BigDFT/multkernel
!! NAME
!!   multkernel
!!
!! FUNCTION
!!     (Based on suitable modifications of S.Goedecker routines)
!!     Multiply with the kernel taking into account its symmetry
!!     Conceived to be used into convolution loops
!!
!! SYNOPSIS
!!     pot:      Kernel, symmetric and real, half the length
!!     zw:       Work array (input/output)
!!     n1,n2:    logical dimension of the FFT transform, reference for zw
!!     nd1,nd2:  Dimensions of POT
!!     jS, nfft: starting point of the plane and number of remaining lines
!!
!! RESTRICTIONS on USAGE
!!     Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!!     Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!!     Copyright (C) 2002 Stefan Goedecker, CEA Grenoble
!!     This file is distributed under the terms of the
!!      GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!!
!! AUTHORS
!!    S. Goedecker, L. Genovese
!!
!! CREATION DATE
!!     February 2006
!!
!! SOURCE
!!
subroutine multkernel(nd1,nd2,n1,n2,lot,nfft,jS,pot,zw)
  implicit none
  !Argments
  integer, intent(in) :: nd1,nd2,n1,n2,lot,nfft,jS
  real(KIND=dp), dimension(nd1,nd2), intent(in) :: pot
  real(KIND=dp), dimension(2,lot,n2), intent(inout) :: zw
  !Local variables
  integer :: j,j1,i2,j2,isign

  !Body
  !case i2=1
  do j=1,nfft
     j1=j+jS-1
     !isign=(j1/(n1/2+2))
     !j1=(1-2*isign)*j1+isign*(n1+2) !n1/2+1-abs(n1/2+2-jS-i1)
     j1=j1+(j1/(n1/2+2))*(n1+2-2*j1)
!!$     j1=n1/2+1-abs(n1/2+2-jS-j)!this stands for j1=min(jS-1+j,n1+3-jS-j)
     zw(1,j,1)=zw(1,j,1)*pot(j1,1)
     zw(2,j,1)=zw(2,j,1)*pot(j1,1)
  end do

  !generic case
  do i2=2,n2/2
     do j=1,nfft
        j1=j+jS-1
        j1=j1+(j1/(n1/2+2))*(n1+2-2*j1)
!!$        j1=n1/2+1-abs(n1/2+2-jS-j)
        j2=n2+2-i2
        zw(1,j,i2)=zw(1,j,i2)*pot(j1,i2)
        zw(2,j,i2)=zw(2,j,i2)*pot(j1,i2)
        zw(1,j,j2)=zw(1,j,j2)*pot(j1,i2)
        zw(2,j,j2)=zw(2,j,j2)*pot(j1,i2)
     end do
  end do
  
  !case i2=n2/2+1
  do j=1,nfft
     j1=j+jS-1
     j1=j1+(j1/(n1/2+2))*(n1+2-2*j1)
!!$     j1=n1/2+1-abs(n1/2+2-jS-j)
     j2=n2/2+1
     zw(1,j,j2)=zw(1,j,j2)*pot(j1,j2)
     zw(2,j,j2)=zw(2,j,j2)*pot(j1,j2)
  end do

end subroutine multkernel




!!!HERE POT MUST BE THE KERNEL (BEWARE THE HALF DIMENSION)

!!****h* BigDFT/S_PoissonSolver
!! NAME
!!   PoissonSolver
!!
!! FUNCTION
!!     (Based on suitable modifications of S.Goedecker routines)
!!     Applies the local FFT space Kernel to the density in Real space.
!!     Does NOT calculate the LDA exchange-correlation terms
!!
!! SYNOPSIS
!!     zf:          Density (input/output)
!!                  ZF(i1,i3,i2)
!!                  i1=1,md1 , i2=1,md2/nproc , i3=1,md3 
!!     pot:         Kernel, only the distributed part (REAL)
!!                  POT(i1,i2,i3)
!!                  i1=1,nd1 , i2=1,nd2 , i3=1,nd3/nproc 
!!     nproc:       number of processors used as returned by MPI_COMM_SIZE
!!     iproc:       [0:nproc-1] number of processor as returned by MPI_COMM_RANK
!!     n1,n2,n3/2:  logical dimension of the transform. As transform lengths 
!!	            most products of the prime factors 2,3,5 are allowed.
!!                  The detailed table with allowed transform lengths can 
!!                  be found in subroutine CTRIG
!!     md1,md2,md3: Dimension of ZF
!!     nd1,nd2,nd3/nproc: Dimension of POT
!!     scal:        factor of renormalization of the FFT in order to acheve unitarity 
!!                  and the correct dimension
!!     hx,hy,hz:    grid spacing, used for integrating eharthree
!!     ehartree:    hartree energy of the potential
!!
!! RESTRICTIONS on USAGE
!!     Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!!     Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!!     Copyright (C) 2002 Stefan Goedecker, CEA Grenoble
!!     This file is distributed under the terms of the
!!      GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!!
!! AUTHORS
!!    S. Goedecker, L. Genovese
!!
!! CREATION DATE
!!    October 2006
!!
!! SOURCE
!!
subroutine S_PoissonSolver(n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,nproc,iproc,pot,zf&
             ,scal,hx,hy,hz,mpi_group,error)!,ehartree)
  implicit none
!  include 'mpif.h'
!  include "perfdata.inc"
integer, parameter :: ncache_optimal=8*1024
integer, parameter :: timing_flag=0
integer :: count_time1,count_time2,count_rate,count_max,number_time,index_time
real(kind=dp) :: time_b0,time_b1,serial_time,parallel_time
  !Arguments
  INTEGER, INTENT(in) :: n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,nproc,iproc,mpi_group
  real(KIND=dp), intent(in) :: scal,hx,hy,hz
  !real(KIND=dp), intent(out) :: ehartree
  real(KIND=dp), dimension(nd1,nd2,nd3/nproc), intent(in) :: pot
  real(KIND=dp), dimension(md1,md3,md2/nproc), intent(inout) :: zf
  TYPE(cp_error_type), INTENT(INOUT)       :: error
  !Local variables
    CHARACTER(len=*), PARAMETER :: routineN = 'S_PoissonSolver', &
      routineP = moduleN//':'//routineN  

  integer :: ncache,lzt,lot,ma,mb,nfft,ic1,ic2,ic3,Jp2stb,J2stb,Jp2stf,J2stf
  INTEGER :: j1,j2,j3,i1,i2,i3,i,j,inzee,ierr,i_all,i_stat,handle
  real(KIND=dp) :: twopion!,ehartreetmp
  !work arrays for transpositions
  real(KIND=dp), dimension(:,:,:), allocatable :: zt
  !work arrays for MPI
  real(KIND=dp), dimension(:,:,:,:,:), allocatable :: zmpi1
  real(KIND=dp), dimension(:,:,:,:), allocatable :: zmpi2
  !cache work array
  real(KIND=dp), dimension(:,:,:), allocatable :: zw
  !FFT work arrays
  real(KIND=dp), dimension(:,:), allocatable :: btrig1,btrig2,btrig3, &
       ftrig1,ftrig2,ftrig3,cosinarr
  integer, dimension(:), allocatable :: after1,now1,before1, & 
       after2,now2,before2,after3,now3,before3
  
  !Body

!  call timing(iproc,'PSolv_comput  ','ON')
    CALL timeset(routineN,"I","",handle)
  ! check input
  !if (mod(n1,2).ne.0) stop 'Parallel convolution:ERROR:n1' !this can be avoided
  !if (mod(n2,2).ne.0) stop 'Parallel convolution:ERROR:n2' !this can be avoided
  if (mod(n3,2).ne.0) stop 'Parallel convolution:ERROR:n3'
  if (nd1.lt.n1/2+1) stop 'Parallel convolution:ERROR:nd1' 
  if (nd2.lt.n2/2+1) stop 'Parallel convolution:ERROR:nd2' 
  if (nd3.lt.n3/2+1) stop 'Parallel convolution:ERROR:nd3'
  if (md1.lt.n1) stop 'Parallel convolution:ERROR:md1'
  if (md2.lt.n2) stop 'Parallel convolution:ERROR:md2'
  if (md3.lt.n3/2) stop 'Parallel convolution:ERROR:md3'
  if (mod(nd3,nproc).ne.0) stop 'Parallel convolution:ERROR:nd3'
  if (mod(md2,nproc).ne.0) stop 'Parallel convolution:ERROR:md2'
  
  !defining work arrays dimensions
  ncache=ncache_optimal
  if (ncache <= max(n1,n2,n3/2)*4) ncache=max(n1,n2,n3/2)*4

!  if (timing_flag == 1 .and. iproc ==0) print *,'parallel ncache=',ncache

  lzt=n2
  if (mod(n2,2).eq.0) lzt=lzt+1
  if (mod(n2,4).eq.0) lzt=lzt+1 !maybe this is useless
  
  !Allocations
  allocate(btrig1(2,8192),stat=i_all)
  allocate(ftrig1(2,8192),stat=i_stat)
  i_all=i_all+i_stat
  allocate(after1(7),stat=i_stat)
  i_all=i_all+i_stat
  allocate(now1(7),stat=i_stat)
  i_all=i_all+i_stat
  allocate(before1(7),stat=i_stat)
  i_all=i_all+i_stat
  allocate(btrig2(2,8192),stat=i_stat)
  i_all=i_all+i_stat
  allocate(ftrig2(2,8192),stat=i_stat)
  i_all=i_all+i_stat
  allocate(after2(7),stat=i_stat)
  i_all=i_all+i_stat
  allocate(now2(7),stat=i_stat)
  i_all=i_all+i_stat
  allocate(before2(7),stat=i_stat)
  i_all=i_all+i_stat
  allocate(btrig3(2,8192),stat=i_stat)
  i_all=i_all+i_stat
  allocate(ftrig3(2,8192),stat=i_stat)
  i_all=i_all+i_stat
  allocate(after3(7),stat=i_stat)
  i_all=i_all+i_stat
  allocate(now3(7),stat=i_stat)
  i_all=i_all+i_stat
  allocate(before3(7),stat=i_stat)
  i_all=i_all+i_stat
  allocate(zw(2,ncache/4,2),stat=i_stat)
  i_all=i_all+i_stat
  allocate(zt(2,lzt,n1),stat=i_stat)
  i_all=i_all+i_stat
  allocate(zmpi2(2,n1,md2/nproc,nd3),stat=i_stat)
  i_all=i_all+i_stat
  allocate(cosinarr(2,n3/2),stat=i_stat)
  i_all=i_all+i_stat
  if (nproc.gt.1)  then
      allocate(zmpi1(2,n1,md2/nproc,nd3/nproc,nproc),stat=i_stat)
      zmpi1= 0.0_dp
  end if
  if (i_all+i_stat /= 0) then
     write(*,*)' s_poissonsolver: problem of memory allocation'
     stop
  end if
  zmpi2 = 0.0_dp
  zw = 0.0_dp
  zt = 0.0_dp

  !calculating the FFT work arrays (beware on the HalFFT in n3 dimension)
  call ctrig(n3/2,btrig3,after3,before3,now3,1,ic3)
  call ctrig(n1,btrig1,after1,before1,now1,1,ic1)
  call ctrig(n2,btrig2,after2,before2,now2,1,ic2)
  do  j=1,n1
     ftrig1(1,j)= btrig1(1,j)
     ftrig1(2,j)=-btrig1(2,j)
  enddo
  do  j=1,n2
     ftrig2(1,j)= btrig2(1,j)
     ftrig2(2,j)=-btrig2(2,j)
  enddo
  do  j=1,n3/2
     ftrig3(1,j)= btrig3(1,j)
     ftrig3(2,j)=-btrig3(2,j)
  enddo

  !Calculating array of phases for HalFFT decoding
  twopion=8._dp*datan(1._dp)/real(n3,KIND=dp)
  do i3=1,n3/2
     cosinarr(1,i3)=dcos(twopion*(i3-1))
     cosinarr(2,i3)=-dsin(twopion*(i3-1))
  end do

  !initializing integral
  !ehartree=0._dp

  ! transform along z axis
  lot=ncache/(2*n3)
  if (lot.lt.1) then  
     write(6,*) & 
          'convolxc_off:ncache has to be enlarged to be able to hold at' // &  
          'least one 1-d FFT of this size even though this will' // & 
          'reduce the performance for shorter transform lengths',n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,nproc,iproc
     stop
  endif
  
  do j2=1,md2/nproc
     !this condition ensures that we manage only the interesting part for the FFT
     if (iproc*(md2/nproc)+j2.le.n2) then
        do i1=1,n1,lot
           ma=i1
           mb=min(i1+(lot-1),n1)
           nfft=mb-ma+1

           !inserting real data into complex array of half lenght
           call halfill_upcorn(md1,md3,lot,nfft,n3,zf(i1,1,j2),zw(1,1,1))

           !performing FFT
           !input: I1,I3,J2,(Jp2)
            inzee=1
           do i=1,ic3
              call fftstp(lot,nfft,n3/2,lot,n3/2,zw(1,1,inzee),zw(1,1,3-inzee), &
                   btrig3,after3(i),now3(i),before3(i),1)
              inzee=3-inzee
           enddo
           !output: I1,i3,J2,(Jp2)
            !unpacking FFT in order to restore correct result, 
           !while exchanging components
           !input: I1,i3,J2,(Jp2)
           call scramble_unpack(i1,j2,lot,nfft,n1,n3,md2,nproc,nd3,zw(1,1,inzee),zmpi2,cosinarr)
           !output: I1,J2,i3,(Jp2)
        end do
     end if
  end do
  !Interprocessor data transposition
  !input: I1,J2,j3,jp3,(Jp2)
  if (nproc.gt.1) then
!     call timing(iproc,'PSolv_comput  ','OF')

     call mp_alltoall(zmpi2,& 
          zmpi1,2*n1*(md2/nproc)*(nd3/nproc), &
          mpi_group)

!     call timing(iproc,'PSolv_commun  ','OF')

!     call timing(iproc,'PSolv_comput  ','ON')
  endif
  !output: I1,J2,j3,Jp2,(jp3)

  !now each process perform complete convolution of its planes
  do j3=1,nd3/nproc
     !this condition ensures that we manage only the interesting part for the FFT
     if (iproc*(nd3/nproc)+j3.le.n3/2+1) then
	Jp2stb=1
	J2stb=1
	Jp2stf=1
	J2stf=1
        
        ! transform along x axis
        lot=ncache/(4*n1)
        if (lot.lt.1) then  
           write(6,*) & 
                'convolxc_off:ncache has to be enlarged to be able to hold at' // &  
                'least one 1-d FFT of this size even though this will' // & 
                'reduce the performance for shorter transform lengths'
           stop
        endif
        
        do j=1,n2,lot
           ma=j
           mb=min(j+(lot-1),n2)
           nfft=mb-ma+1

           !reverse index ordering, leaving the planes to be transformed at the end
           !input: I1,J2,j3,Jp2,(jp3)
           if (nproc.eq.1) then
              call S_mpiswitch_upcorn(j3,nfft,Jp2stb,J2stb,lot,n1,md2,nd3,nproc,zmpi2,zw(1,1,1))
           else
              call S_mpiswitch_upcorn(j3,nfft,Jp2stb,J2stb,lot,n1,md2,nd3,nproc,zmpi1,zw(1,1,1))
           endif
           !output: J2,Jp2,I1,j3,(jp3)

           !performing FFT
           !input: I2,I1,j3,(jp3)
           inzee=1
           do i=1,ic1-1
              call fftstp(lot,nfft,n1,lot,n1,zw(1,1,inzee),zw(1,1,3-inzee), &
                   btrig1,after1(i),now1(i),before1(i),1)
              inzee=3-inzee
           enddo

           !storing the last step into zt array
           i=ic1
           call fftstp(lot,nfft,n1,lzt,n1,zw(1,1,inzee),zt(1,j,1), & 
                btrig1,after1(i),now1(i),before1(i),1)           
           !output: I2,i1,j3,(jp3)
        end do

        !transform along y axis
        lot=ncache/(4*n2)
        if (lot.lt.1) then  
           write(6,*) & 
                'convolxc_off:ncache has to be enlarged to be able to hold at' // &  
                'least one 1-d FFT of this size even though this will' // & 
                'reduce the performance for shorter transform lengths'
           stop
        endif

        do j=1,n1,lot
           ma=j
           mb=min(j+(lot-1),n1)
           nfft=mb-ma+1

           !reverse ordering 
           !input: I2,i1,j3,(jp3)
           call S_switch_upcorn(nfft,n2,lot,n1,lzt,zt(1,1,j),zw(1,1,1))
           !output: i1,I2,j3,(jp3)
           
           !performing FFT
           !input: i1,I2,j3,(jp3)
           inzee=1
           do i=1,ic2
              call fftstp(lot,nfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee), &
                   btrig2,after2(i),now2(i),before2(i),1)
              inzee=3-inzee
           enddo
           !output: i1,i2,j3,(jp3)
           
           !Multiply with kernel in fourier space
           call multkernel(nd1,nd2,n1,n2,lot,nfft,j,pot(1,1,j3),zw(1,1,inzee))
           
           !TRANSFORM BACK IN REAL SPACE
           
           !transform along y axis
           !input: i1,i2,j3,(jp3)
           do i=1,ic2
              call fftstp(lot,nfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee), &
                   ftrig2,after2(i),now2(i),before2(i),-1)
              inzee=3-inzee
           enddo

           !reverse ordering
           !input: i1,I2,j3,(jp3)
           call S_unswitch_downcorn(nfft,n2,lot,n1,lzt,zw(1,1,inzee),zt(1,1,j))
           !output: I2,i1,j3,(jp3)
        end do
        
        !transform along x axis
        !input: I2,i1,j3,(jp3)
        lot=ncache/(4*n1)
        do j=1,n2,lot
           ma=j
           mb=min(j+(lot-1),n2)
           nfft=mb-ma+1

           !performing FFT
           i=1
           call fftstp(lzt,nfft,n1,lot,n1,zt(1,j,1),zw(1,1,1), &
                ftrig1,after1(i),now1(i),before1(i),-1)
           
           inzee=1
           do i=2,ic1
              call fftstp(lot,nfft,n1,lot,n1,zw(1,1,inzee),zw(1,1,3-inzee), &
                   ftrig1,after1(i),now1(i),before1(i),-1)
              inzee=3-inzee
           enddo
           !output: I2,I1,j3,(jp3)

           !reverse ordering
           !input: J2,Jp2,I1,j3,(jp3)
           if (nproc.eq.1) then
              call S_unmpiswitch_downcorn(j3,nfft,Jp2stf,J2stf,lot,n1,md2,nd3,nproc,zw(1,1,inzee),zmpi2)
           else
              call S_unmpiswitch_downcorn(j3,nfft,Jp2stf,J2stf,lot,n1,md2,nd3,nproc,zw(1,1,inzee),zmpi1)
           endif
           ! output: I1,J2,j3,Jp2,(jp3)
        end do
     endif

  end do


  !Interprocessor data transposition
  !input: I1,J2,j3,Jp2,(jp3)
  if (nproc.gt.1) then
!     call timing(iproc,'PSolv_comput  ','OF')

!     call timing(iproc,'PSolv_commun  ','ON')
     !communication scheduling
     call mp_alltoall(zmpi1,& 
          zmpi2,2*n1*(md2/nproc)*(nd3/nproc), &
          mpi_group)
!     call timing(iproc,'PSolv_commun  ','OF')

!     call timing(iproc,'PSolv_comput  ','ON')
  endif

 !output: I1,J2,j3,jp3,(Jp2)

  !transform along z axis
  !input: I1,J2,i3,(Jp2)
  lot=ncache/(2*n3)
  do j2=1,md2/nproc
     !this condition ensures that we manage only the interesting part for the FFT
     if (iproc*(md2/nproc)+j2.le.n2) then
        do i1=1,n1,lot
           ma=i1
           mb=min(i1+(lot-1),n1)
           nfft=mb-ma+1

           !reverse ordering and repack the FFT data in order to be backward HalFFT transformed
           !input: I1,J2,i3,(Jp2)
           call unscramble_pack(i1,j2,lot,nfft,n1,n3,md2,nproc,nd3,zmpi2,zw(1,1,1),cosinarr)
           !output: I1,i3,J2,(Jp2)

           !performing FFT
           !input: I1,i3,J2,(Jp2)           
           inzee=1
           do i=1,ic3
              call fftstp(lot,nfft,n3/2,lot,n3/2,zw(1,1,inzee),zw(1,1,3-inzee), &
                   ftrig3,after3(i),now3(i),before3(i),-1)
              inzee=3-inzee
           enddo
           !output: I1,I3,J2,(Jp2)

           !rebuild the output array
           call unfill_downcorn(md1,md3,lot,nfft,n3,zw(1,1,inzee),zf(i1,1,j2)&
                ,scal)!,ehartreetmp)

           !integrate local pieces together
           !ehartree=ehartree+0.5_dp*ehartreetmp*hx*hy*hz
        end do
     endif
  end do

  !De-allocations  
  deallocate(btrig1,stat=i_all)
  deallocate(ftrig1,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(after1,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(now1,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(before1,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(btrig2,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(ftrig2,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(after2,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(now2,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(before2,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(btrig3,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(ftrig3,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(after3,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(now3,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(before3,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(zmpi2,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(zw,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(zt,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(cosinarr,stat=i_stat)
  i_all=i_all+i_stat
  if (nproc.gt.1) deallocate(zmpi1,stat=i_stat)
  if (i_all+i_stat /= 0) then
     write(*,*)' s_poissonsolver: problem of memory deallocation'
     stop
  end if

!  call timing(iproc,'PSolv_comput  ','OF')
    CALL timestop(0.0_dp,handle)
end subroutine S_PoissonSolver

subroutine S_mpiswitch_upcorn(j3,nfft,Jp2stb,J2stb,lot,n1,md2,nd3,nproc,zmpi1,zw)
  implicit none
  integer, intent(in) :: j3,nfft,lot,n1,md2,nd3,nproc
  integer, intent(inout) :: Jp2stb,J2stb
  real(KIND=dp), dimension(2,n1,md2/nproc,nd3/nproc,nproc), intent(in) :: zmpi1
  real(KIND=dp), dimension(2,lot,n1), intent(inout) :: zw
  !local variables
  integer :: I1,J2,Jp2,mfft

  mfft=0
  do Jp2=Jp2stb,nproc
     do J2=J2stb,md2/nproc
	mfft=mfft+1
	if (mfft.gt.nfft) then
           Jp2stb=Jp2
           J2stb=J2
           return
	end if
        do I1=1,n1
           zw(1,mfft,I1)=zmpi1(1,I1,J2,j3,Jp2)
           zw(2,mfft,I1)=zmpi1(2,I1,J2,j3,Jp2)
        end do
     end do
     J2stb=1
  end do

end subroutine S_mpiswitch_upcorn

subroutine S_switch_upcorn(nfft,n2,lot,n1,lzt,zt,zw)
  implicit none
  integer, intent(in) :: nfft,n2,lot,n1,lzt
  real(KIND=dp), dimension(2,lzt,n1), intent(in) :: zt
  real(KIND=dp), dimension(2,lot,n2), intent(inout) :: zw
  !local variables
  integer :: i,j

  do j=1,nfft
     do i=1,n2
	zw(1,j,i)=zt(1,i,j)
	zw(2,j,i)=zt(2,i,j)
     end do
  end do

end subroutine S_switch_upcorn

subroutine S_unswitch_downcorn(nfft,n2,lot,n1,lzt,zw,zt)
  implicit none
  integer, intent(in) :: nfft,n2,lot,n1,lzt
  real(KIND=dp), dimension(2,lot,n2), intent(in) :: zw
  real(KIND=dp), dimension(2,lzt,n1), intent(inout) :: zt
  !local variables
  integer :: i,j

  do j=1,nfft
     do i=1,n2
	zt(1,i,j)=zw(1,j,i)
	zt(2,i,j)=zw(2,j,i)
     end do
  end do

end subroutine S_unswitch_downcorn

subroutine S_unmpiswitch_downcorn(j3,nfft,Jp2stf,J2stf,lot,n1,md2,nd3,nproc,zw,zmpi1)
  implicit none
  integer, intent(in) :: j3,nfft,lot,n1,md2,nd3,nproc
  integer, intent(inout) :: Jp2stf,J2stf
  real(KIND=dp), dimension(2,lot,n1), intent(in) :: zw
  real(KIND=dp), dimension(2,n1,md2/nproc,nd3/nproc,nproc), intent(inout) :: zmpi1
  !local variables
  integer :: I1,J2,Jp2,mfft

  mfft=0
  do Jp2=Jp2stf,nproc
     do J2=J2stf,md2/nproc
	mfft=mfft+1
	if (mfft.gt.nfft) then
           Jp2stf=Jp2
           J2stf=J2
           return
	end if
        do I1=1,n1
           zmpi1(1,I1,J2,j3,Jp2)=zw(1,mfft,I1)
           zmpi1(2,I1,J2,j3,Jp2)=zw(2,mfft,I1)
        end do
     end do
     J2stf=1
  end do
end subroutine S_unmpiswitch_downcorn


!!****h* BigDFT/unfill_downcorn
!! NAME
!!   unfill_downcorn
!!
!! FUNCTION
!!     (Based on suitable modifications of S.Goedecker routines)
!!     Restore data into output array, calculating in the meanwhile
!!     Hartree energy of the potential 
!!
!! SYNOPSIS
!!     zf:          Original distributed density as well as
!!                  Distributed solution of the poisson equation (inout)
!!     zw:          FFT work array
!!     n3:          (twice the) dimension of the last FFTtransform.
!!     md1,md3:     Dimensions of the undistributed part of the real grid
!!     nfft:        number of planes
!!     scal:        Needed to achieve unitarity and correct dimensions
!!     ehartreetmp: Hartree energy
!!
!! WARNING
!!     Assuming that high frequencies are in the corners 
!!     and that n3 is multiple of 4   
!!
!! RESTRICTIONS on USAGE
!!     Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!!     Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!!     Copyright (C) 2002 Stefan Goedecker, CEA Grenoble
!!     This file is distributed under the terms of the
!!     GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!!
!! AUTHORS
!!    S. Goedecker, L. Genovese
!!
!! CREATION DATE
!!     February 2006
!!
!! SOURCE
!!
subroutine unfill_downcorn(md1,md3,lot,nfft,n3,zw,zf&
     ,scal)!,ehartreetmp)
  implicit none
  !Arguments
  integer, intent(in) :: md1,md3,lot,nfft,n3
  real(KIND=dp), dimension(2,lot,n3/2), intent(in) :: zw
  real(KIND=dp), dimension(md1,md3),intent(inout) :: zf
  real(KIND=dp), intent(in) :: scal
  !real(KIND=dp), intent(out) :: ehartreetmp
  !Local variables
  integer :: i3,i1
  real(KIND=dp) :: pot1

  !Body
  !ehartreetmp=0._dp
  do i3=1,n3/4
     do i1=1,nfft
        pot1 = scal*zw(1,i1,i3)
        !ehartreetmp =ehartreetmp + pot1* zf(i1,2*i3-1)
	zf(i1,2*i3-1)= pot1 
	pot1 = scal*zw(2,i1,i3)
        !ehartreetmp =ehartreetmp + pot1* zf(i1,2*i3)
	zf(i1,2*i3)= pot1 
     enddo
  end do
  
end subroutine unfill_downcorn

        subroutine halfill_upcorn(md1,md3,lot,nfft,n3,zf,zw)
	implicit real*8 (a-h,o-z)
	dimension zw(2,lot,n3/2),zf(md1,md3)
! WARNING: Assuming that high frequencies are in the corners 
!          and that n3 is multiple of 4
!in principle we can relax this condition
	
        do 90,i3=1,n3/4
	do i1=1,nfft
	zw(1,i1,i3)=0._dp
	zw(2,i1,i3)=0._dp
	enddo
90      continue
        do 100,i3=n3/4+1,n3/2
	do i1=1,nfft
	zw(1,i1,i3)=zf(i1,2*i3-1-n3/2)
	zw(2,i1,i3)=zf(i1,2*i3-n3/2)
	enddo
100     continue


	return
	end subroutine halfill_upcorn


!!****h* BigDFT/scramble_unpack
!! NAME
!!   scramble_unpack
!!
!! FUNCTION
!!     (Based on suitable modifications of S.Goedecker routines)
!!     Assign the correct planes to the work array zmpi2
!!     in order to prepare for interprocessor data transposition.
!!     In the meanwhile, it unpacks the data of the HalFFT in order to prepare for
!!     multiplication with the kernel
!!
!! SYNOPSIS
!!     zmpi2:          Work array for multiprocessor manipulation (output)
!!     zw:             Work array (input)
!!     cosinarr:      Array of the phases needed for unpacking
!!     n1,n3:          logical dimension of the FFT transform, reference for work arrays
!!     md2,nd3:        Dimensions of real grid and of the kernel, respectively
!!     i1,j2,lot,nfft: Starting points of the plane and number of remaining lines
!!
!! RESTRICTIONS on USAGE
!!     Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!!     Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!!     Copyright (C) 2002 Stefan Goedecker, CEA Grenoble
!!     This file is distributed under the terms of the
!!     GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!!
!! AUTHORS
!!    S. Goedecker, L. Genovese
!!
!! CREATION DATE
!!     February 2006
!!
!! SOURCE
!!
subroutine scramble_unpack(i1,j2,lot,nfft,n1,n3,md2,nproc,nd3,zw,zmpi2,cosinarr)
  implicit none
  !Arguments
  integer, intent(in) :: i1,j2,lot,nfft,n1,n3,md2,nproc,nd3
  real(KIND=dp), dimension(2,lot,n3/2), intent(in) :: zw
  real(KIND=dp), dimension(2,n3/2), intent(in) :: cosinarr
  real(KIND=dp), dimension(2,n1,md2/nproc,nd3), intent(inout) :: zmpi2
  !Local variables
  integer :: i3,i,ind1,ind2
  real(KIND=dp) ::  a,b,c,d,cp,sp,feR,feI,foR,foI,fR,fI
  
  !case i3=1 and i3=n3/2+1
  do i=0,nfft-1
     a=zw(1,i+1,1)
     b=zw(2,i+1,1)
     zmpi2(1,i1+i,j2,1)=a+b
     zmpi2(2,i1+i,j2,1)=0._dp
     zmpi2(1,i1+i,j2,n3/2+1)=a-b
     zmpi2(2,i1+i,j2,n3/2+1)=0._dp
  end do
  !case 2<=i3<=n3/2
  do i3=2,n3/2
     ind1=i3
     ind2=n3/2-i3+2
     cp=cosinarr(1,i3)
     sp=cosinarr(2,i3)
     do i=0,nfft-1
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
     end do
  end do

end subroutine scramble_unpack

 
!!****h* BigDFT/unscramble_pack
!! NAME
!!   unscramble_pack
!!
!! FUNCTION
!!     (Based on suitable modifications of S.Goedecker routines)
!!     Insert the correct planes of the work array zmpi2
!!     in order to prepare for backward FFT transform
!!     In the meanwhile, it packs the data in order to be transformed with the HalFFT 
!!     procedure
!!
!! SYNOPSIS
!!     zmpi2:          Work array for multiprocessor manipulation (input)
!!     zw:             Work array (output)
!!     cosinarr:       Array of the phases needed for packing
!!     n1,n3:          logical dimension of the FFT transform, reference for work arrays
!!     md2,nd3:        Dimensions of real grid and of the kernel, respectively
!!     i1,j2,lot,nfft: Starting points of the plane and number of remaining lines
!!
!! RESTRICTIONS on USAGE
!!     Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!!     Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!!     Copyright (C) 2002 Stefan Goedecker, CEA Grenoble
!!     This file is distributed under the terms of the
!!     GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!!
!! AUTHORS
!!    S. Goedecker, L. Genovese
!!
!! CREATION DATE
!!     February 2006
!!
!! SOURCE
!!
subroutine unscramble_pack(i1,j2,lot,nfft,n1,n3,md2,nproc,nd3,zmpi2,zw,cosinarr)
  implicit none
  !Arguments
  integer, intent(in) :: i1,j2,lot,nfft,n1,n3,md2,nproc,nd3
  real(KIND=dp), dimension(2,lot,n3/2), intent(inout) :: zw
  real(KIND=dp), dimension(2,n3/2), intent(in) :: cosinarr
  real(KIND=dp), dimension(2,n1,md2/nproc,nd3), intent(in) :: zmpi2
  !Local variables
  integer :: i3,i,indA,indB
  real(KIND=dp) ::  a,b,c,d,cp,sp,re,ie,ro,io,rh,ih

  !Body
  do i3=1,n3/2
     indA=i3
     indB=n3/2+2-i3
     cp=cosinarr(1,i3)
     sp=cosinarr(2,i3)
     do i=0,nfft-1
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
     end do
  end do

end subroutine unscramble_pack


!!****h* BigDFT/F_PoissonSolver
!! NAME
!!   F_PoissonSolver
!!
!! FUNCTION
!!     (Based on suitable modifications of S.Goedecker routines)
!!     Applies the local FFT space Kernel to the density in Real space.
!!     Calculates also the LDA exchange-correlation terms
!!
!! SYNOPSIS
!!     zf:          Density (input/output)
!!                  ZF(i1,i3,i2)
!!                  i1=1,md1 , i2=1,md2/nproc , i3=1,md3 
!!     pot:         Kernel, only the distributed part (REAL)
!!                  POT(i1,i2,i3)
!!                  i1=1,nd1 , i2=1,nd2 , i3=1,nd3/nproc
!!     zfpot_ion:   Distributed array of the ionization potential 
!!     nproc:       number of processors used as returned by MPI_COMM_SIZE
!!     iproc:       [0:nproc-1] number of processor as returned by MPI_COMM_RANK
!!     n1,n2,n3:    logical dimension of the transform. As transform lengths 
!!	            most products of the prime factors 2,3,5 are allowed.
!!                  The detailed table with allowed transform lengths can 
!!                  be found in subroutine CTRIG
!!     md1,md2,md3: Dimension of ZF
!!     nd1,nd2,nd3: Dimension of POT
!!     scal:        factor of renormalization of the FFT in order to acheve unitarity 
!!                  and the correct dimension
!!     hgrid:       grid spacing, used for integrating quantities
!!     ehartree:    hartree energy of the potential
!!     eexcu,vexcu: LDA exchange correlation terms   
!!
!! RESTRICTIONS on USAGE
!!     Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!!     Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!!     Copyright (C) 2002 Stefan Goedecker, CEA Grenoble
!!     This file is distributed under the terms of the
!!      GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!!
!! AUTHORS
!!    S. Goedecker, L. Genovese
!!
!! CREATION DATE
!!     February 2006
!!
!! SOURCE
!!
subroutine F_PoissonSolver(n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,nproc,iproc,pot,zf&
             ,scal,hgrid,mpi_group,error)!,ehartree)
  implicit none
integer, parameter :: ncache_optimal=8*1024
integer, parameter :: timing_flag=0
integer :: count_time1,count_time2,count_rate,count_max,number_time,index_time
real(kind=dp) :: time_b0,time_b1,serial_time,parallel_time
  !Arguments
  INTEGER, INTENT(in) :: n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,nproc,iproc,mpi_group
  real(KIND=dp), intent(in) :: scal,hgrid
  !real(KIND=dp), intent(out) :: ehartree
  real(KIND=dp), dimension(nd1,nd2,nd3/nproc), intent(in) :: pot
  real(KIND=dp), dimension(md1,md3,md2/nproc), intent(inout) :: zf
  TYPE(cp_error_type), INTENT(INOUT)       :: error
  !Local variables
  integer :: ncache,lzt,lot,ma,mb,nfft,ic1,ic2,ic3,Jp2stb,J2stb,Jp2stf,J2stf
  integer :: j1,j2,j3,i1,i2,i3,i,j,inzee,ierr,i_all,i_stat
  real(KIND=dp) :: twopion!,ehartreetmp
  !work arrays for transpositions
  real(KIND=dp), dimension(:,:,:), allocatable :: zt
  !work arrays for MPI
  real(KIND=dp), dimension(:,:,:,:,:), allocatable :: zmpi1
  real(KIND=dp), dimension(:,:,:,:), allocatable :: zmpi2
  !cache work array
  real(KIND=dp), dimension(:,:,:), allocatable :: zw
  !FFT work arrays
  real(KIND=dp), dimension(:,:), allocatable :: btrig1,btrig2,btrig3,&
       ftrig1,ftrig2,ftrig3,cosinarr
  integer, dimension(:), allocatable :: after1,now1,before1, & 
       after2,now2,before2,after3,now3,before3

!  call timing(iproc,'PSolv_comput  ','ON')
  
  !Body
  ! check input
  if (mod(n1,2).ne.0) stop 'Parallel convolution:ERROR:n1'
  if (mod(n2,2).ne.0) stop 'Parallel convolution:ERROR:n2'
  if (mod(n3,2).ne.0) stop 'Parallel convolution:ERROR:n3'
  if (nd1.lt.n1/2+1) stop 'Parallel convolution:ERROR:nd1'
  if (nd2.lt.n2/2+1) stop 'Parallel convolution:ERROR:nd2'
  if (nd3.lt.n3/2+1) stop 'Parallel convolution:ERROR:nd3'
  if (md1.lt.n1/2) stop 'Parallel convolution:ERROR:md1'
  if (md2.lt.n2/2) stop 'Parallel convolution:ERROR:md2'
  if (md3.lt.n3/2) stop 'Parallel convolution:ERROR:md3'
  if (mod(nd3,nproc).ne.0) stop 'Parallel convolution:ERROR:nd3'
  if (mod(md2,nproc).ne.0) stop 'Parallel convolution:ERROR:md2'
  
  !defining work arrays dimensions
  
  ncache=ncache_optimal
  if (ncache <= max(n1,n2,n3/2)*4) ncache=max(n1,n2,n3/2)*4
  lzt=n2/2
  if (mod(n2/2,2).eq.0) lzt=lzt+1
  if (mod(n2/2,4).eq.0) lzt=lzt+1
  
  !Allocations
  allocate(btrig1(2,8192),stat=i_all)
  allocate(ftrig1(2,8192),stat=i_stat)
  i_all=i_all+i_stat
  allocate(after1(7),stat=i_stat)
  i_all=i_all+i_stat
  allocate(now1(7),stat=i_stat)
  i_all=i_all+i_stat
  allocate(before1(7),stat=i_stat)
  i_all=i_all+i_stat
  allocate(btrig2(2,8192),stat=i_stat)
  i_all=i_all+i_stat
  allocate(ftrig2(2,8192),stat=i_stat)
  i_all=i_all+i_stat
  allocate(after2(7),stat=i_stat)
  i_all=i_all+i_stat
  allocate(now2(7),stat=i_stat)
  i_all=i_all+i_stat
  allocate(before2(7),stat=i_stat)
  i_all=i_all+i_stat
  allocate(btrig3(2,8192),stat=i_stat)
  i_all=i_all+i_stat
  allocate(ftrig3(2,8192),stat=i_stat)
  i_all=i_all+i_stat
  allocate(after3(7),stat=i_stat)
  i_all=i_all+i_stat
  allocate(now3(7),stat=i_stat)
  i_all=i_all+i_stat
  allocate(before3(7),stat=i_stat)
  i_all=i_all+i_stat
  allocate(zw(2,ncache/4,2),stat=i_stat)
  i_all=i_all+i_stat
  allocate(zt(2,lzt,n1),stat=i_stat)
  i_all=i_all+i_stat
  allocate(zmpi2(2,n1,md2/nproc,nd3),stat=i_stat)
  i_all=i_all+i_stat
  allocate(cosinarr(2,n3/2),stat=i_stat)
  i_all=i_all+i_stat
  if (nproc.gt.1) allocate(zmpi1(2,n1,md2/nproc,nd3/nproc,nproc),stat=i_stat)
  if (i_all+i_stat /= 0) then
     write(*,*)' f_poissonsolver: problem of memory allocation'
     stop
  end if


  !calculating the FFT work arrays (beware on the HalFFT in n3 dimension)
  call ctrig(n3/2,btrig3,after3,before3,now3,1,ic3)
  call ctrig(n1,btrig1,after1,before1,now1,1,ic1)
  call ctrig(n2,btrig2,after2,before2,now2,1,ic2)
  do  j=1,n1
     ftrig1(1,j)= btrig1(1,j)
     ftrig1(2,j)=-btrig1(2,j)
  enddo
  do  j=1,n2
     ftrig2(1,j)= btrig2(1,j)
     ftrig2(2,j)=-btrig2(2,j)
  enddo
  do  j=1,n3
     ftrig3(1,j)= btrig3(1,j)
     ftrig3(2,j)=-btrig3(2,j)
  enddo

  !Calculating array of phases for HalFFT decoding
  twopion=8._dp*datan(1._dp)/real(n3,KIND=dp)
  do i3=1,n3/2
     cosinarr(1,i3)=dcos(twopion*(i3-1))
     cosinarr(2,i3)=-dsin(twopion*(i3-1))
  end do

  !initializing integrals
  !ehartree=0._dp

  ! transform along z axis
  lot=ncache/(2*n3)
  if (lot.lt.1) then  
     write(6,*) & 
          'convolxc_on:ncache has to be enlarged to be able to hold at' // &  
          'least one 1-d FFT of this size even though this will' // & 
          'reduce the performance for shorter transform lengths'
     stop
  endif
  
  do j2=1,md2/nproc
     !this condition ensures that we manage only the interesting part for the FFT
     if (iproc*(md2/nproc)+j2.le.n2/2) then
        do i1=1,(n1/2),lot
           ma=i1
           mb=min(i1+(lot-1),(n1/2))
           nfft=mb-ma+1

           !inserting real data into complex array of half lenght
           call halfill_upcorn(md1,md3,lot,nfft,n3,zf(i1,1,j2),zw(1,1,1))

           !performing FFT
           !input: I1,I3,J2,(Jp2)
           inzee=1
           do i=1,ic3
              call fftstp(lot,nfft,n3/2,lot,n3/2,zw(1,1,inzee),zw(1,1,3-inzee), &
                   btrig3,after3(i),now3(i),before3(i),1)
              inzee=3-inzee
           enddo
           !output: I1,i3,J2,(Jp2)

           !unpacking FFT in order to restore correct result, 
           !while exchanging components
           !input: I1,i3,J2,(Jp2)
           call scramble_unpack(i1,j2,lot,nfft,n1/2,n3,md2,nproc,nd3,zw(1,1,inzee),zmpi2,cosinarr)
           !output: I1,J2,i3,(Jp2)
        end do
     endif
  end do

  !Interprocessor data transposition
  !input: I1,J2,j3,jp3,(Jp2)
  if (nproc.gt.1) then
!     call timing(iproc,'PSolv_comput  ','OF')

!     call timing(iproc,'PSolv_commun  ','ON')
     !communication scheduling
     call mp_alltoall(zmpi2,& !n1*(md2/nproc)*(nd3/nproc), &
!          MPI_double_precision, &
          zmpi1,n1*(md2/nproc)*(nd3/nproc), &
          mpi_group)
!     call timing(iproc,'PSolv_commun  ','OF')

!     call timing(iproc,'PSolv_comput  ','ON')
  endif
  !output: I1,J2,j3,Jp2,(jp3)

  !now each process perform complete convolution of its planes
  do j3=1,nd3/nproc
     !this condition ensures that we manage only the interesting part for the FFT
     if (iproc*(nd3/nproc)+j3.le.n3/2+1) then
	Jp2stb=1
	J2stb=1
	Jp2stf=1
	J2stf=1
        
        ! transform along x axis
        lot=ncache/(4*n1)
        if (lot.lt.1) then  
           write(6,*) & 
                'convolxc_on:ncache has to be enlarged to be able to hold at' // &  
                'least one 1-d FFT of this size even though this will' // & 
                'reduce the performance for shorter transform lengths'
           stop
        endif
        
        do j=1,n2/2,lot
           ma=j
           mb=min(j+(lot-1),n2/2)
           nfft=mb-ma+1

           !reverse index ordering, leaving the planes to be transformed at the end
           !input: I1,J2,j3,Jp2,(jp3)
           if (nproc.eq.1) then
              call mpiswitch_upcorn(j3,nfft,Jp2stb,J2stb,lot,n1,md2,nd3,nproc,zmpi2,zw(1,1,1))
           else
              call mpiswitch_upcorn(j3,nfft,Jp2stb,J2stb,lot,n1,md2,nd3,nproc,zmpi1,zw(1,1,1))
           endif
           !output: J2,Jp2,I1,j3,(jp3)
           
           !performing FFT
           !input: I2,I1,j3,(jp3)
           inzee=1
           do i=1,ic1-1
              call fftstp(lot,nfft,n1,lot,n1,zw(1,1,inzee),zw(1,1,3-inzee), &
                   btrig1,after1(i),now1(i),before1(i),1)
              inzee=3-inzee
           enddo

           !storing the last step into zt array
           i=ic1
           call fftstp(lot,nfft,n1,lzt,n1,zw(1,1,inzee),zt(1,j,1), & 
                btrig1,after1(i),now1(i),before1(i),1)           
           !output: I2,i1,j3,(jp3)
        end do

        !transform along y axis
        lot=ncache/(4*n2)
        if (lot.lt.1) then  
           write(6,*) & 
                'convolxc_on:ncache has to be enlarged to be able to hold at' // &  
                'least one 1-d FFT of this size even though this will' // & 
                'reduce the performance for shorter transform lengths'
           stop
        endif

        do j=1,n1,lot
           ma=j
           mb=min(j+(lot-1),n1)
           nfft=mb-ma+1

           !reverse ordering 
           !input: I2,i1,j3,(jp3)
           call switch_upcorn(nfft,n2,lot,n1,lzt,zt(1,1,j),zw(1,1,1))
           !output: i1,I2,j3,(jp3)
           
           !performing FFT
           !input: i1,I2,j3,(jp3)
           inzee=1
           do i=1,ic2
              call fftstp(lot,nfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee), &
                   btrig2,after2(i),now2(i),before2(i),1)
              inzee=3-inzee
           enddo
           !output: i1,i2,j3,(jp3)
           
           !Multiply with kernel in fourier space
           call multkernel(nd1,nd2,n1,n2,lot,nfft,j,pot(1,1,j3),zw(1,1,inzee))
           
           !TRANSFORM BACK IN REAL SPACE
           
           !transform along y axis
           !input: i1,i2,j3,(jp3)
           do i=1,ic2
              call fftstp(lot,nfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee), &
                   ftrig2,after2(i),now2(i),before2(i),-1)
              inzee=3-inzee
           enddo

           !reverse ordering
           !input: i1,I2,j3,(jp3)
           call unswitch_downcorn(nfft,n2,lot,n1,lzt,zw(1,1,inzee),zt(1,1,j))
           !output: I2,i1,j3,(jp3)
        end do
        
        !transform along x axis
        !input: I2,i1,j3,(jp3)
        lot=ncache/(4*n1)
        do j=1,n2/2,lot
           ma=j
           mb=min(j+(lot-1),n2/2)
           nfft=mb-ma+1

           !performing FFT
           i=1
           call fftstp(lzt,nfft,n1,lot,n1,zt(1,j,1),zw(1,1,1), &
                ftrig1,after1(i),now1(i),before1(i),-1)
           
           inzee=1
           do i=2,ic1
              call fftstp(lot,nfft,n1,lot,n1,zw(1,1,inzee),zw(1,1,3-inzee), &
                   ftrig1,after1(i),now1(i),before1(i),-1)
              inzee=3-inzee
           enddo
           !output: I2,I1,j3,(jp3)

           !reverse ordering
           !input: J2,Jp2,I1,j3,(jp3)
           if (nproc.eq.1) then
              call unmpiswitch_downcorn(j3,nfft,Jp2stf,J2stf,lot,n1,md2,nd3,nproc,zw(1,1,inzee),zmpi2)
           else
              call unmpiswitch_downcorn(j3,nfft,Jp2stf,J2stf,lot,n1,md2,nd3,nproc,zw(1,1,inzee),zmpi1)
           endif
           ! output: I1,J2,j3,Jp2,(jp3)
        end do
     endif
  end do
  
  
  !Interprocessor data transposition
  !input: I1,J2,j3,Jp2,(jp3)
  if (nproc.gt.1) then
!     call timing(iproc,'PSolv_comput  ','OF')

!     call timing(iproc,'PSolv_commun  ','ON')
     !communication scheduling
     call mp_alltoall(zmpi1,& !n1*(md2/nproc)*(nd3/nproc), &
!          MPI_double_precision, &
          zmpi2,n1*(md2/nproc)*(nd3/nproc), &
          mpi_group)
!     call timing(iproc,'PSolv_commun  ','OF')

!     call timing(iproc,'PSolv_comput  ','ON')
     !output: I1,J2,j3,jp3,(Jp2)
  endif

  !transform along z axis
  !input: I1,J2,i3,(Jp2)
  lot=ncache/(2*n3)
  do j2=1,md2/nproc
     !this condition ensures that we manage only the interesting part for the FFT
     if (iproc*(md2/nproc)+j2.le.n2/2) then
        do i1=1,(n1/2),lot
           ma=i1
           mb=min(i1+(lot-1),(n1/2))
           nfft=mb-ma+1

           !reverse ordering and repack the FFT data in order to be backward HalFFT transformed
           !input: I1,J2,i3,(Jp2)
           call unscramble_pack(i1,j2,lot,nfft,n1/2,n3,md2,nproc,nd3,zmpi2,zw(1,1,1),cosinarr)
           !output: I1,i3,J2,(Jp2)

           !performing FFT
           !input: I1,i3,J2,(Jp2)           
           inzee=1
           do i=1,ic3
              call fftstp(lot,nfft,n3/2,lot,n3/2,zw(1,1,inzee),zw(1,1,3-inzee), &
                   ftrig3,after3(i),now3(i),before3(i),-1)
              inzee=3-inzee
           enddo
           !output: I1,I3,J2,(Jp2)

           !calculates the exchange correlation terms locally and rebuild the output array
           call unfill_downcorn(md1,md3,lot,nfft,n3,zw(1,1,inzee),zf(i1,1,j2)&
                ,scal)!,ehartreetmp)

           !integrate local pieces together
           !ehartree=ehartree+0.5_dp*ehartreetmp*hgrid**3
        end do
     endif
  end do
  
  !De-allocations  
  deallocate(btrig1,stat=i_all)
  deallocate(ftrig1,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(after1,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(now1,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(before1,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(btrig2,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(ftrig2,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(after2,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(now2,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(before2,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(btrig3,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(ftrig3,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(after3,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(now3,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(before3,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(zmpi2,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(zw,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(zt,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(cosinarr,stat=i_stat)
  i_all=i_all+i_stat
  if (nproc.gt.1) deallocate(zmpi1,stat=i_stat)
  if (i_all+i_stat /= 0) then
     write(*,*)' f_poissonsolver: problem of memory deallocation'
     stop
  end if

!  call timing(iproc,'PSolv_comput  ','OF')
end subroutine F_PoissonSolver

 	subroutine switch_upcorn(nfft,n2,lot,n1,lzt,zt,zw)
	implicit real*8 (a-h,o-z)
	dimension zw(2,lot,n2),zt(2,lzt,n1)
! WARNING: Assuming that high frequencies are in the corners 
!          and that n2 is multiple of 2

! Low frequencies 
	do 100,j=1,nfft
        do 100,i=n2/2+1,n2
	zw(1,j,i)=zt(1,i-n2/2,j)
	zw(2,j,i)=zt(2,i-n2/2,j)
100	continue

! High frequencies 
        do 90,i=1,n2/2
	do 90,j=1,nfft
	zw(1,j,i)=0._dp
	zw(2,j,i)=0._dp
90	continue

	return
	end subroutine switch_upcorn

        
 	subroutine mpiswitch_upcorn(j3,nfft,Jp2stb,J2stb,lot,n1,md2,nd3,nproc,zmpi1,zw)
	implicit real*8 (a-h,o-z)
        dimension zmpi1(2,n1/2,md2/nproc,nd3/nproc,nproc),zw(2,lot,n1)
! WARNING: Assuming that high frequencies are in the corners 
!          and that n1 is multiple of 2

	mfft=0
	do 300,Jp2=Jp2stb,nproc
	do 200,J2=J2stb,md2/nproc
	mfft=mfft+1
	if (mfft.gt.nfft) then
	Jp2stb=Jp2
	J2stb=J2
	return
	endif
        do 90,I1=1,n1/2
	zw(1,mfft,I1)=0._dp
	zw(2,mfft,I1)=0._dp
90	continue
        do 100,I1=n1/2+1,n1
	zw(1,mfft,I1)=zmpi1(1,I1-n1/2,J2,j3,Jp2)
	zw(2,mfft,I1)=zmpi1(2,I1-n1/2,J2,j3,Jp2)
100	continue
200	continue
	J2stb=1
300	continue
	end subroutine mpiswitch_upcorn

 	subroutine unswitch_downcorn(nfft,n2,lot,n1,lzt,zw,zt)
	implicit real*8 (a-h,o-z)
	dimension zw(2,lot,n2),zt(2,lzt,n1)
! WARNING: Assuming that high frequencies are in the corners 
!          and that n2 is multiple of 2

! Low frequencies
	do 100,j=1,nfft
        do 100,i=1,n2/2
	zt(1,i,j)=zw(1,j,i)
	zt(2,i,j)=zw(2,j,i)
100	continue
	return
	end subroutine unswitch_downcorn

 	subroutine unmpiswitch_downcorn(j3,nfft,Jp2stf,J2stf,lot,n1,md2,nd3,nproc,zw,zmpi1)
	implicit real*8 (a-h,o-z)
        dimension zmpi1(2,n1/2,md2/nproc,nd3/nproc,nproc),zw(2,lot,n1)
! WARNING: Assuming that high frequencies are in the corners 
!          and that n1 is multiple of 2

	mfft=0
	do 300,Jp2=Jp2stf,nproc
	DO 200,J2=J2stf,md2/nproc
	mfft=mfft+1
	if (mfft.gt.nfft) then
	Jp2stf=Jp2
	J2stf=J2
	return
	endif
        do 100,I1=1,n1/2
	zmpi1(1,I1,J2,j3,Jp2)=zw(1,mfft,I1)
	zmpi1(2,I1,J2,j3,Jp2)=zw(2,mfft,I1)
100	continue
200	continue
	J2stf=1
300	continue
	end subroutine unmpiswitch_downcorn

!!****h* BigDFT/F_unfill_downcorn
!! NAME
!!   F_unfill_downcorn
!!
!! FUNCTION
!!     (Based on suitable modifications of S.Goedecker routines)
!!     Restore data into output array, calculating in the meanwhile
!!     Hartree energy of the potential 
!!
!! SYNOPSIS
!!     zf:          Original distributed density as well as
!!                  Distributed solution of the poisson equation (inout)
!!     zw:          FFT work array
!!     n3:          (twice the) dimension of the last FFTtransform.
!!     md1,md3:     Dimensions of the undistributed part of the real grid
!!     nfft:        number of planes
!!     scal:        Needed to achieve unitarity and correct dimensions
!!     ehartreetmp: Hartree energy
!!
!! WARNING
!!     Assuming that high frequencies are in the corners 
!!     and that n3 is multiple of 4   
!!
!! RESTRICTIONS on USAGE
!!     Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!!     Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!!     Copyright (C) 2002 Stefan Goedecker, CEA Grenoble
!!     This file is distributed under the terms of the
!!     GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!!
!! AUTHORS
!!    S. Goedecker, L. Genovese
!!
!! CREATION DATE
!!     February 2006
!!
!! SOURCE
!!
subroutine F_unfill_downcorn(md1,md3,lot,nfft,n3,zw,zf&
     ,scal,ehartreetmp)
  implicit none
  !Arguments
  integer, intent(in) :: md1,md3,lot,nfft,n3
  real(KIND=dp), dimension(2,lot,n3/2), intent(in) :: zw
  real(KIND=dp), dimension(md1,md3),intent(inout) :: zf
  real(KIND=dp), intent(in) :: scal
  real(KIND=dp), intent(out) :: ehartreetmp
  !Local variables
  integer :: i3,i1
  real(KIND=dp) :: pot1

  !Body
  ehartreetmp=0._dp
  do i3=1,n3/4
     do i1=1,nfft
        pot1 = scal*zw(1,i1,i3)
        ehartreetmp =ehartreetmp + pot1* zf(i1,2*i3-1)
	zf(i1,2*i3-1)= pot1 
	pot1 = scal*zw(2,i1,i3)
        ehartreetmp =ehartreetmp + pot1* zf(i1,2*i3)
	zf(i1,2*i3)= pot1 
     enddo
  end do
  
end subroutine F_unfill_downcorn
END MODULE
