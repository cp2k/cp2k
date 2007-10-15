!!****** cp2k/ ps_wavelet_scaling_function
!!
!!   NAME
!!     ps_wavelet_kernel
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


MODULE ps_wavelet_scaling_function
  USE kinds,                           ONLY: dp
  USE lazy,                            ONLY: lazy_arrays

PUBLIC :: scaling_function,&
          scf_recursion

CONTAINS

!!****h* BigDFT/scaling_function
!! NAME
!!   scaling_function
!!
!! FUNCTION
!!   Calculate the values of a scaling function in real uniform grid
!!
!! SOURCE
!!
SUBROUTINE scaling_function(itype,nd,nrange,a,x)

  IMPLICIT NONE
  !Arguments
  !Type of interpolating functions
  INTEGER, INTENT(in) :: itype
  !Number of points: must be 2**nex
  INTEGER, INTENT(in) :: nd
  INTEGER, INTENT(out) :: nrange
  REAL(KIND=dp), DIMENSION(0:nd), INTENT(out) :: a,x
  !Local variables
  REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: y
  REAL(KIND=dp),POINTER, DIMENSION(:)         ::   ch,cg,cgt,cht
  INTEGER :: i,nt,ni,i_all,m
  

  m=itype+2
  CALL lazy_arrays(itype,m,ch,cg,cgt,cht)

  ni=2*itype
  nrange = ni
  ALLOCATE(y(0:nd),stat=i_all)
  IF (i_all /= 0) THEN
     WRITE(*,*)' scaling_function: problem of memory allocation'
     STOP
  END IF
  
  ! plot scaling function
  CALL zero(nd+1,x)
  CALL zero(nd+1,y)
  nt=ni
  x(nt/2-1)=1.d0
  loop1: DO
     nt=2*nt

     CALL back_trans(nd,nt,x,y,m,ch,cg,cgt,cht)
     CALL dcopy(nt,y,1,x,1)
     IF (nt.EQ.nd) THEN
        EXIT loop1
     END IF
  END DO loop1
  
  !open (unit=1,file='scfunction',status='unknown')
  DO i=0,nd
     a(i) = 1.d0*i*ni/nd-(.5d0*ni-1.d0)
     !write(1,*) 1.d0*i*ni/nd-(.5d0*ni-1.d0),x(i)
  END DO
  !close(1)
  DEALLOCATE(ch,cg,cgt,cht)
  DEALLOCATE(y,stat=i_all)
  IF (i_all /= 0) THEN
     WRITE(*,*)' scaling_function: problem of memory deallocation'
     STOP
  END IF
END SUBROUTINE scaling_function
!!***


!!****h* BigDFT/wavelet_function
!! NAME
!!   wavelet_function
!!
!! FUNCTION
!!   Calculate the values of the wavelet function in a real uniform mesh.
!!
!! SOURCE
!!
SUBROUTINE wavelet_function(itype,nd,a,x)
  
  IMPLICIT NONE
  !Arguments
  !Type of the interpolating scaling function
  INTEGER, INTENT(in) :: itype
  !must be 2**nex
  INTEGER, INTENT(in) :: nd
  REAL(KIND=dp), DIMENSION(0:nd), INTENT(out) :: a,x
  !Local variables
  REAL(KIND=dp),POINTER, DIMENSION(:)         ::   ch,cg,cgt,cht
  REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: y
  INTEGER :: i,nt,ni,i_all,m
  
  m=itype+2
  ni=2*itype
  CALL lazy_arrays(itype,m,ch,cg,cgt,cht)
  ALLOCATE(y(0:nd),stat=i_all)
  IF (i_all /= 0) THEN
     WRITE(*,*)' wavelet_function: problem of memory allocation'
     STOP
  END IF
  
  ! plot wavelet 
  CALL zero(nd+1,x)
  CALL zero(nd+1,y)
  nt=ni
  x(nt+nt/2-1)=1.d0
  loop3: DO
     nt=2*nt
     !write(6,*) 'nd,nt',nd,nt
     CALL back_trans(nd,nt,x,y,m,ch,cg,cgt,cht)
     CALL dcopy(nd,y,1,x,1)
     IF (nt.EQ.nd) THEN
        EXIT loop3
     END IF
  END DO loop3
  
  !open (unit=1,file='wavelet',status='unknown')
  DO i=0,nd-1
     a(i) = 1.d0*i*ni/nd-(.5d0*ni-.5d0)
     !write(1,*) 1.d0*i*ni/nd-(.5d0*ni-.5d0),x(i)
  END DO
  !close(1)
  DEALLOCATE(ch,cg,cgt,cht)  
  DEALLOCATE(y,stat=i_all)
  IF (i_all /= 0) THEN
     WRITE(*,*)' wavelet_function: problem of memory deallocation'
     STOP
  END IF
  
END SUBROUTINE wavelet_function
!!***


!!****h* BigDFT/scf_recursion
!! NAME
!!   scf_recursion
!!
!! FUNCTION
!!   Do iterations to go from p0gauss to pgauss
!!   order interpolating scaling function
!!
!! SOURCE
!!
SUBROUTINE scf_recursion(itype,n_iter,n_range,kernel_scf,kern_1_scf)
  IMPLICIT NONE
  !Arguments
  INTEGER, INTENT(in) :: itype,n_iter,n_range
  real(KIND=dp), intent(inout) :: kernel_scf(-n_range:n_range)
  REAL(KIND=dp), INTENT(out) :: kern_1_scf(-n_range:n_range)
  INTEGER:: m
  REAL(KIND=dp),POINTER, DIMENSION(:)         ::   ch,cg,cgt,cht
  !Local variables
  
  m=itype+2
  CALL lazy_arrays(itype,m,ch,cg,cgt,cht)
  CALL scf_recurs(n_iter,n_range,kernel_scf,kern_1_scf,m,ch,cg,cgt,cht)
  DEALLOCATE(ch,cg,cgt,cht)  

END SUBROUTINE scf_recursion
!!***


!!****h* BigDFT/zero
!! NAME
!!   zero
!!
!! FUNCTION
!!   Set to zero an array x(n)
!!
!! SOURCE
!!
SUBROUTINE zero(n,x)
  IMPLICIT NONE
  !Arguments
  INTEGER, INTENT(in) :: n
  REAL(KIND=dp), INTENT(out) :: x(n)
  !Local variables
  INTEGER :: i
  DO i=1,n
     x(i)=0.d0
  END DO
END SUBROUTINE zero
!!***


!!****h* BigDFT/for_trans_8
!! NAME
!!   for_trans_8
!!
!! FUNCTION
!!   forward wavelet transform
!!   nd: length of data set
!!   nt length of data in data set to be transformed
!!   m filter length (m has to be even!)
!!   x input data, y output data
!!
!! SOURCE
!!
SUBROUTINE for_trans(nd,nt,x,y,m,ch,cg,cgt,cht)
  implicit none
  !Arguments
  integer, intent(in) :: nd,nt
  real(KIND=dp), intent(in) :: x(0:nd-1)
  real(KIND=dp), intent(out) :: y(0:nd-1)
  REAL(KIND=dp),POINTER, DIMENSION(:)         ::   ch,cg,cgt,cht
  !Local variables
  INTEGER :: i,j,ind,m


  do i=0,nt/2-1
     y(     i)=0.d0
     y(nt/2+i)=0.d0
     
     do j=-m+1,m
        
        ! periodically wrap index if necessary
        ind=j+2*i
        loop99: do
           if (ind.lt.0) then 
              ind=ind+nt
              cycle loop99
           end if
           if (ind.ge.nt) then 
              ind=ind-nt
              cycle loop99
           end if
           exit loop99
        end do loop99

        y(     i)=y(     i)+cht(j)*x(ind)
        y(nt/2+i)=y(nt/2+i)+cgt(j)*x(ind)
     end do
     
  end do

end subroutine for_trans
!!***


!!****h* BigDFT/back_trans
!! NAME
!!   back_trans_8
!!
!! FUNCTION
!!
!! SOURCE
!!
subroutine back_trans(nd,nt,x,y,m,ch,cg,cgt,cht)
  ! backward wavelet transform
  ! nd: length of data set
  ! nt length of data in data set to be transformed
  ! m filter length (m has to be even!)
  ! x input data, y output data
  implicit none
  !Arguments
  integer, intent(in) :: nd,nt
  real(KIND=dp), intent(in) :: x(0:nd-1)
  real(KIND=dp), intent(out) :: y(0:nd-1)
  REAL(KIND=dp),POINTER, DIMENSION(:)         ::   ch,cg,cgt,cht
  INTEGER :: i,j,ind,m


  
  do i=0,nt/2-1
     y(2*i+0)=0.d0
     y(2*i+1)=0.d0
     
     do j=-m/2,m/2-1
        
        ! periodically wrap index if necessary
        ind=i-j
        loop99: do
           if (ind.lt.0) then 
              ind=ind+nt/2
              cycle loop99
           end if
           if (ind.ge.nt/2) then 
              ind=ind-nt/2
              cycle loop99
           end if
           exit loop99
        end do loop99

        y(2*i+0)=y(2*i+0) + ch(2*j-0)*x(ind)+cg(2*j-0)*x(ind+nt/2)
        y(2*i+1)=y(2*i+1) + ch(2*j+1)*x(ind)+cg(2*j+1)*x(ind+nt/2)
     end do
        
  end do
        
end subroutine back_trans
!!***


!!****h* BigDFT/ftest
!! NAME
!!   ftest_8
!!
!! FUNCTION
!!   Tests the 4 orthogonality relations of the filters
!!
!! SOURCE
!!
subroutine ftest(m,ch,cg,cgt,cht)
  implicit none
  !Arguments
  !Local variables
  character(len=*), parameter :: fmt22 = "(a,i3,i4,4(e17.10))"
  INTEGER :: i,j,l,m
  real(KIND=dp) :: t1,t2,t3,t4,eps
  REAL(KIND=dp),POINTER, DIMENSION(:)         ::   ch,cg,cgt,cht

  
  ! do i=-m,m
  ! write(6,*) i,ch(i),cg(i)
  ! end do
  
  do i=-m,m
     do j=-m,m
        t1=0.d0
        t2=0.d0
        t3=0.d0
        t4=0.d0
        do l=-3*m,3*m
           if ( l-2*i.ge.-m .and. l-2*i.le.m  .and. &
                l-2*j.ge.-m .and. l-2*j.le.m ) then
              t1=t1+ch(l-2*i)*cht(l-2*j)
              t2=t2+cg(l-2*i)*cgt(l-2*j)
              t3=t3+ch(l-2*i)*cgt(l-2*j)
              t4=t4+cht(l-2*i)*cg(l-2*j)
           end if
        end do
        eps=1.d-10
        if (i.eq.j) then
           if (abs(t1-1.d0).gt.eps .or. abs(t2-1.d0).gt.eps .or. &
             & abs(t3).gt.eps  .or. abs(t4).gt.eps ) then 
              write(6,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           end if
        else
           if (abs(t1).gt.eps .or. abs(t2).gt.eps .or. &
             & abs(t3).gt.eps  .or. abs(t4).gt.eps ) then
              write(6,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           end if
        end if
     end do
  end do
  
  write(6,*) 'FILTER TEST PASSED'
  
end subroutine ftest
!!***


!!****h* BigDFT/scf_recursion_8
!! NAME
!!   scf_recursion_8
!!
!! FUNCTION
!!   Do iterations to go from p0gauss to pgauss
!!   8th-order interpolating scaling function
!!
!! SOURCE
!!
SUBROUTINE scf_recurs(n_iter,n_range,kernel_scf,kern_1_scf,m,ch,cg,cgt,cht)
  implicit none
  !Arguments
  integer, intent(in) :: n_iter,n_range
  real(KIND=dp), intent(inout) :: kernel_scf(-n_range:n_range)
  real(KIND=dp), intent(out) :: kern_1_scf(-n_range:n_range)
  !Local variables
  real(KIND=dp) :: kern,kern_tot
  INTEGER :: i_iter,i,j,ind,m
  REAL(KIND=dp),POINTER, DIMENSION(:)         ::   ch,cg,cgt,cht


  !Start the iteration to go from p0gauss to pgauss
  loop_iter_scf: do i_iter=1,n_iter
     kern_1_scf(:) = kernel_scf(:)
     kernel_scf(:) = 0.d0
     loop_iter_i: do i=0,n_range
        kern_tot = 0.d0
        do j=-m,m
           ind = 2*i-j
           if (abs(ind) > n_range) then
              kern = 0.d0
           else
              kern = kern_1_scf(ind)
           end if
           kern_tot = kern_tot + ch(j)*kern
        end do
        if (kern_tot == 0.d0) then
           !zero after (be sure because strictly == 0.d0)
           exit loop_iter_i
        else
           kernel_scf( i) = 0.5d0*kern_tot
           kernel_scf(-i) = kernel_scf(i)
        end if
     end do loop_iter_i
  end do loop_iter_scf
END SUBROUTINE scf_recurs
!!***



END MODULE ps_wavelet_scaling_function
