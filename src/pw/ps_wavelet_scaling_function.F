!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Creates the wavelet kernel for the wavelet based poisson solver.
!> \author Florian Schiffmann (09.2007,fschiff)
! *****************************************************************************
MODULE ps_wavelet_scaling_function
  USE kinds,                           ONLY: dp
  USE lazy,                            ONLY: lazy_arrays
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: scaling_function,&
            scf_recursion

CONTAINS

! *****************************************************************************
!> \brief Calculate the values of a scaling function in real uniform grid
!> \param itype ...
!> \param nd ...
!> \param nrange ...
!> \param a ...
!> \param x ...
! *****************************************************************************
SUBROUTINE scaling_function(itype,nd,nrange,a,x)

  !Type of interpolating functions
    INTEGER, INTENT(in)                      :: itype, nd
    INTEGER, INTENT(out)                     :: nrange
    REAL(KIND=dp), DIMENSION(0:nd), &
      INTENT(out)                            :: a, x

    INTEGER                                  :: i, i_all, m, ni, nt
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: y
    REAL(KIND=dp), DIMENSION(:), POINTER     :: cg, cgt, ch, cht

!Number of points: must be 2**nex

  a = 0.0_dp
  x = 0.0_dp
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
  x(nt/2-1)=1._dp
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
     a(i) = 1._dp*i*ni/nd-(.5_dp*ni-1._dp)
     !write(1,*) 1._dp*i*ni/nd-(.5_dp*ni-1._dp),x(i)
  END DO
  !close(1)
  DEALLOCATE(ch,cg,cgt,cht)
  DEALLOCATE(y,stat=i_all)
  IF (i_all /= 0) THEN
     WRITE(*,*)' scaling_function: problem of memory deallocation'
     STOP
  END IF
END SUBROUTINE scaling_function

! *****************************************************************************
!> \brief Calculate the values of the wavelet function in a real uniform mesh.
!> \param itype ...
!> \param nd ...
!> \param a ...
!> \param x ...
! *****************************************************************************
SUBROUTINE wavelet_function(itype,nd,a,x)

  !Type of the interpolating scaling function
    INTEGER, INTENT(in)                      :: itype, nd
    REAL(KIND=dp), DIMENSION(0:nd), &
      INTENT(out)                            :: a, x

    INTEGER                                  :: i, i_all, m, ni, nt
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: y
    REAL(KIND=dp), DIMENSION(:), POINTER     :: cg, cgt, ch, cht

!must be 2**nex

  a = 0.0_dp
  x = 0.0_dp
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
  x(nt+nt/2-1)=1._dp
  loop3: DO
     nt=2*nt
     !WRITE(*,*) 'nd,nt',nd,nt
     CALL back_trans(nd,nt,x,y,m,ch,cg,cgt,cht)
     CALL dcopy(nd,y,1,x,1)
     IF (nt.EQ.nd) THEN
        EXIT loop3
     END IF
  END DO loop3

  !open (unit=1,file='wavelet',status='unknown')
  DO i=0,nd-1
     a(i) = 1._dp*i*ni/nd-(.5_dp*ni-.5_dp)
     !write(1,*) 1._dp*i*ni/nd-(.5_dp*ni-.5_dp),x(i)
  END DO
  !close(1)
  DEALLOCATE(ch,cg,cgt,cht)
  DEALLOCATE(y,stat=i_all)
  IF (i_all /= 0) THEN
     WRITE(*,*)' wavelet_function: problem of memory deallocation'
     STOP
  END IF

END SUBROUTINE wavelet_function

! *****************************************************************************
!> \brief Do iterations to go from p0gauss to pgauss
!>    order interpolating scaling function
!> \param itype ...
!> \param n_iter ...
!> \param n_range ...
!> \param kernel_scf ...
!> \param kern_1_scf ...
! *****************************************************************************
SUBROUTINE scf_recursion(itype,n_iter,n_range,kernel_scf,kern_1_scf)
    INTEGER, INTENT(in)                      :: itype, n_iter, n_range
    REAL(KIND=dp), INTENT(inout)             :: kernel_scf(-n_range:n_range)
    REAL(KIND=dp), INTENT(out)               :: kern_1_scf(-n_range:n_range)

    INTEGER                                  :: m
    REAL(KIND=dp), DIMENSION(:), POINTER     :: cg, cgt, ch, cht

  kern_1_scf = 0.0_dp
  m=itype+2
  CALL lazy_arrays(itype,m,ch,cg,cgt,cht)
  CALL scf_recurs(n_iter,n_range,kernel_scf,kern_1_scf,m,ch,cg,cgt,cht)
  DEALLOCATE(ch,cg,cgt,cht)

END SUBROUTINE scf_recursion

! *****************************************************************************
!> \brief Set to zero an array x(n)
!> \param n ...
!> \param x ...
! *****************************************************************************
SUBROUTINE zero(n,x)
    INTEGER, INTENT(in)                      :: n
    REAL(KIND=dp), INTENT(out)               :: x(n)

    INTEGER                                  :: i

  DO i=1,n
     x(i)=0._dp
  END DO
END SUBROUTINE zero

! *****************************************************************************
!> \brief forward wavelet transform
!>    nd: length of data set
!>    nt length of data in data set to be transformed
!>    m filter length (m has to be even!)
!>    x input data, y output data
!> \param nd ...
!> \param nt ...
!> \param x ...
!> \param y ...
!> \param m ...
!> \param ch ...
!> \param cg ...
!> \param cgt ...
!> \param cht ...
! *****************************************************************************
SUBROUTINE for_trans(nd,nt,x,y,m,ch,cg,cgt,cht)
    INTEGER, INTENT(in)                      :: nd, nt
    REAL(KIND=dp), INTENT(in)                :: x(0:nd-1)
    REAL(KIND=dp), INTENT(out)               :: y(0:nd-1)
    INTEGER                                  :: m
    REAL(KIND=dp), DIMENSION(:), POINTER     :: ch, cg, cgt, cht

    INTEGER                                  :: i, ind, j

  y = 0.0_dp
  DO i=0,nt/2-1
     y(     i)=0._dp
     y(nt/2+i)=0._dp

     DO j=-m+1,m

        ! periodically wrap index if necessary
        ind=j+2*i
        loop99: DO
           IF (ind.lt.0) THEN
              ind=ind+nt
              CYCLE loop99
           END IF
           IF (ind.ge.nt) THEN
              ind=ind-nt
              CYCLE loop99
           END IF
           EXIT loop99
        END DO loop99

        y(     i)=y(     i)+cht(j)*x(ind)
        y(nt/2+i)=y(nt/2+i)+cgt(j)*x(ind)
     END DO

  END DO

END SUBROUTINE for_trans

! *****************************************************************************
!> \brief ...
!> \param nd ...
!> \param nt ...
!> \param x ...
!> \param y ...
!> \param m ...
!> \param ch ...
!> \param cg ...
!> \param cgt ...
!> \param cht ...
! *****************************************************************************
SUBROUTINE back_trans(nd,nt,x,y,m,ch,cg,cgt,cht)
  ! backward wavelet transform
  ! nd: length of data set
  ! nt length of data in data set to be transformed
  ! m filter length (m has to be even!)
  ! x input data, y output data
    INTEGER, INTENT(in)                      :: nd, nt
    REAL(KIND=dp), INTENT(in)                :: x(0:nd-1)
    REAL(KIND=dp), INTENT(out)               :: y(0:nd-1)
    INTEGER                                  :: m
    REAL(KIND=dp), DIMENSION(:), POINTER     :: ch, cg, cgt, cht

    INTEGER                                  :: i, ind, j

  y = 0.0_dp

  DO i=0,nt/2-1
     y(2*i+0)=0._dp
     y(2*i+1)=0._dp

     DO j=-m/2,m/2-1

        ! periodically wrap index if necessary
        ind=i-j
        loop99: DO
           IF (ind.lt.0) THEN
              ind=ind+nt/2
              CYCLE loop99
           END IF
           IF (ind.ge.nt/2) THEN
              ind=ind-nt/2
              CYCLE loop99
           END IF
           EXIT loop99
        END DO loop99

        y(2*i+0)=y(2*i+0) + ch(2*j-0)*x(ind)+cg(2*j-0)*x(ind+nt/2)
        y(2*i+1)=y(2*i+1) + ch(2*j+1)*x(ind)+cg(2*j+1)*x(ind+nt/2)
     END DO

  END DO

END SUBROUTINE back_trans

! *****************************************************************************
!> \brief Tests the 4 orthogonality relations of the filters
!> \param m ...
!> \param ch ...
!> \param cg ...
!> \param cgt ...
!> \param cht ...
! *****************************************************************************
SUBROUTINE ftest(m,ch,cg,cgt,cht)
    INTEGER                                  :: m
    REAL(KIND=dp), DIMENSION(:), POINTER     :: ch, cg, cgt, cht

    CHARACTER(len=*), PARAMETER              :: fmt22 = "(a,i3,i4,4(e17.10))"

    INTEGER                                  :: i, j, l
    REAL(KIND=dp)                            :: eps, t1, t2, t3, t4

! do i=-m,m
! WRITE(*,*) i,ch(i),cg(i)
! end do

  DO i=-m,m
     DO j=-m,m
        t1=0._dp
        t2=0._dp
        t3=0._dp
        t4=0._dp
        DO l=-3*m,3*m
           IF ( l-2*i.ge.-m .AND. l-2*i.le.m  .AND. &
                l-2*j.ge.-m .AND. l-2*j.le.m ) THEN
              t1=t1+ch(l-2*i)*cht(l-2*j)
              t2=t2+cg(l-2*i)*cgt(l-2*j)
              t3=t3+ch(l-2*i)*cgt(l-2*j)
              t4=t4+cht(l-2*i)*cg(l-2*j)
           END IF
        END DO
        eps=1.e-10_dp
        IF (i.eq.j) THEN
           IF (ABS(t1-1._dp).gt.eps .OR. ABS(t2-1._dp).gt.eps .OR. &
               ABS(t3).gt.eps  .OR. ABS(t4).gt.eps ) THEN
              WRITE(*,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           END IF
        ELSE
           IF (ABS(t1).gt.eps .OR. ABS(t2).gt.eps .OR. &
               ABS(t3).gt.eps  .OR. ABS(t4).gt.eps ) THEN
              WRITE(*,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           END IF
        END IF
     END DO
  END DO

  WRITE(*,*) 'FILTER TEST PASSED'

END SUBROUTINE ftest

! *****************************************************************************
!> \brief Do iterations to go from p0gauss to pgauss
!>    8th-order interpolating scaling function
!> \param n_iter ...
!> \param n_range ...
!> \param kernel_scf ...
!> \param kern_1_scf ...
!> \param m ...
!> \param ch ...
!> \param cg ...
!> \param cgt ...
!> \param cht ...
! *****************************************************************************
SUBROUTINE scf_recurs(n_iter,n_range,kernel_scf,kern_1_scf,m,ch,cg,cgt,cht)
    INTEGER, INTENT(in)                      :: n_iter, n_range
    REAL(KIND=dp), INTENT(inout)             :: kernel_scf(-n_range:n_range)
    REAL(KIND=dp), INTENT(out)               :: kern_1_scf(-n_range:n_range)
    INTEGER                                  :: m
    REAL(KIND=dp), DIMENSION(:), POINTER     :: ch, cg, cgt, cht

    INTEGER                                  :: i, i_iter, ind, j
    REAL(KIND=dp)                            :: kern, kern_tot

  kern_1_scf = 0.0_dp
  !Start the iteration to go from p0gauss to pgauss
  loop_iter_scf: DO i_iter=1,n_iter
     kern_1_scf(:) = kernel_scf(:)
     kernel_scf(:) = 0._dp
     loop_iter_i: DO i=0,n_range
        kern_tot = 0._dp
        DO j=-m,m
           ind = 2*i-j
           IF (ABS(ind) > n_range) THEN
              kern = 0._dp
           ELSE
              kern = kern_1_scf(ind)
           END IF
           kern_tot = kern_tot + ch(j)*kern
        END DO
        IF (kern_tot == 0._dp) THEN
           !zero after (be sure because strictly == 0._dp)
           EXIT loop_iter_i
        ELSE
           kernel_scf( i) = 0.5_dp*kern_tot
           kernel_scf(-i) = kernel_scf(i)
        END IF
     END DO loop_iter_i
  END DO loop_iter_scf
END SUBROUTINE scf_recurs

END MODULE ps_wavelet_scaling_function
