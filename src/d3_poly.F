!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!>   \brief
!>     Routines to efficently handle dense polynomial in 3 variables up to
!>     a given degree.
!>     Multiplication, partial evalution, affine transform (change of reference
!>     system), differentiation are efficiently implemented.
!>     some functions accept or return several polynomial together,
!>     these have to have all the same size, and are stored one after the other
!>     in an unique 1d array. This gives them an easy handling and even seem to
!>     be faster than the transposed layout.
!>   \note
!>     not all routines have been fully optimized.
!>     original available also with a BSD style license
!>   \author Fawzi Mohamed
! *****************************************************************************
MODULE d3_poly
  USE f77_blas
  USE kinds,                           ONLY: dp
#include "./common/cp_common_uses.f90"

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: poly_size1,poly_size2,poly_size3,&
        grad_size1,grad_size2,grad_size3,poly_mult1,poly_mult2,poly_mult3,&
        poly_write1,poly_write2,poly_write3,poly_random,poly_affine_t3,&
        poly_p_eval3,poly_p_eval2,poly_eval1,poly_eval2,poly_p_eval2b,&
        poly_eval3,poly_derive3,poly_p_eval3b,poly_cp2k2d3,&
        poly_padd_uneval3,poly_padd_uneval3b,poly_padd_uneval2,poly_padd_uneval2b,&
        poly_affine_t3t,poly_d32cp2k

#ifdef FD_DEBUG
#define IF_CHECK(x,y) x
#else
#define IF_CHECK(x,y) y
#endif

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'd3_poly'
! maximum grad for cached values
    INTEGER, PUBLIC, PARAMETER :: max_grad2=5
    INTEGER, PUBLIC, PARAMETER :: max_grad3=3
    INTEGER, PUBLIC, PARAMETER :: cached_dim1=max_grad2+1
    INTEGER, PUBLIC, PARAMETER :: cached_dim2=(max_grad2+1)*(max_grad2+2)/2
    INTEGER, PUBLIC, PARAMETER :: cached_dim3=(max_grad3+1)*(max_grad3+2)*(max_grad3+3)/6


! cached index -> monomial exponents
    LOGICAL, SAVE, PRIVATE :: module_initialized=.FALSE.
    INTEGER, SAVE, DIMENSION(2,cached_dim2) :: a_mono_exp2
    INTEGER, SAVE, DIMENSION(3,cached_dim3) :: a_mono_exp3
    INTEGER, SAVE, DIMENSION(cached_dim3)   :: a_reduce_idx3
    INTEGER, SAVE, DIMENSION(3,cached_dim3) :: a_deriv_idx3
    INTEGER, SAVE, DIMENSION(cached_dim2,cached_dim2) :: a_mono_mult2
    INTEGER, SAVE, DIMENSION(cached_dim3,cached_dim3) :: a_mono_mult3
    INTEGER, SAVE, DIMENSION(4,cached_dim3) :: a_mono_mult3a
!$OMP THREADPRIVATE(module_initialized,a_mono_exp2,a_mono_exp3,a_reduce_idx3,a_deriv_idx3,a_mono_mult2,a_mono_mult3,a_mono_mult3a)

CONTAINS

! *****************************************************************************
!> \brief initialization of the cache, is called by functions in this module
!> that use cached values
! *****************************************************************************
SUBROUTINE init_d3_poly_module()
    INTEGER                                  :: grad, i, ii, ij, j, subG
    INTEGER, DIMENSION(2)                    :: monoRes2
    INTEGER, DIMENSION(3)                    :: monoRes3

    IF (module_initialized) RETURN
    ii=1
    DO grad=0,max_grad2
        DO i=grad,0,-1
            a_mono_exp2(1,ii)=i
            a_mono_exp2(2,ii)=grad-i
            ii=ii+1
        END DO
    END DO
    ii=1
    DO grad=0,max_grad3
        DO i=grad,0,-1
            DO j=grad-i,0,-1
                a_mono_exp3(1,ii)=i
                a_mono_exp3(2,ii)=j
                a_mono_exp3(3,ii)=grad-i-j
                ii=ii+1
            END DO
        END DO
    END DO
    DO ii=1,cached_dim3
        subG=a_mono_exp3(2,ii)+a_mono_exp3(3,ii)
        a_reduce_idx3(ii)=subG*(subG+1)/2+a_mono_exp3(3,ii)+1
    END DO
    DO ii=1,cached_dim3
        IF (a_mono_exp3(1,ii)>0) THEN
            a_deriv_idx3(1,ii)=mono_index3(a_mono_exp3(1,ii)-1,a_mono_exp3(2,ii),a_mono_exp3(3,ii))
        ELSE
            a_deriv_idx3(1,ii)=0
        END IF
        IF (a_mono_exp3(2,ii)>0) THEN
            a_deriv_idx3(2,ii)=mono_index3(a_mono_exp3(1,ii),a_mono_exp3(2,ii)-1,a_mono_exp3(3,ii))
        ELSE
            a_deriv_idx3(2,ii)=0
        END IF
        IF (a_mono_exp3(3,ii)>0) THEN
            a_deriv_idx3(3,ii)=mono_index3(a_mono_exp3(1,ii),a_mono_exp3(2,ii),a_mono_exp3(3,ii)-1)
        ELSE
            a_deriv_idx3(3,ii)=0
        END IF
    END DO
    DO ii=1,cached_dim2
        DO ij=ii,cached_dim2
            monoRes2=a_mono_exp2(:,ii)+a_mono_exp2(:,ij)
            a_mono_mult2(ii,ij)=mono_index2(monoRes2(1),monoRes2(2))+1
            a_mono_mult2(ij,ii)=a_mono_mult2(ii,ij)
        END DO
    END DO
    DO ii=1,cached_dim3
        DO ij=ii,cached_dim3
            monoRes3=a_mono_exp3(:,ii)+a_mono_exp3(:,ij)
            a_mono_mult3(ii,ij)=mono_index3(monoRes3(1),monoRes3(2),monoRes3(3))+1
            a_mono_mult3(ij,ii)=a_mono_mult3(ii,ij)
        END DO
    END DO
    ii=1
    DO i=1,cached_dim3
       DO j=1,4
          a_mono_mult3a(j,i)=a_mono_mult3(j,i)
          ii=ii+1
       END DO
    END DO
    module_initialized=.TRUE.
END SUBROUTINE

! *****************************************************************************
!> \brief size of a polynomial in x up to the given degree
! *****************************************************************************
PURE FUNCTION poly_size1(maxgrad) RESULT(res)
    INTEGER, INTENT(in)                      :: maxgrad
    INTEGER                                  :: res

    res=maxgrad+1
END FUNCTION

! *****************************************************************************
!> \brief size of a polynomial in x,y up to the given degree
! *****************************************************************************
PURE FUNCTION poly_size2(maxgrad) RESULT(res)
    INTEGER, INTENT(in)                      :: maxgrad
    INTEGER                                  :: res

    res=(maxgrad+1)*(maxgrad+2)/2
END FUNCTION

! *****************************************************************************
!> \brief size of a polynomial in x,y,z up to the given degree
! *****************************************************************************
PURE FUNCTION poly_size3(maxgrad) RESULT(res)
    INTEGER, INTENT(in)                      :: maxgrad
    INTEGER                                  :: res

    res=(maxgrad+1)*(maxgrad+2)*(maxgrad+3)/6
END FUNCTION

! *****************************************************************************
!> \brief max grad for a polynom of the given size
! *****************************************************************************
PURE FUNCTION grad_size1(n) RESULT(res)
    INTEGER, INTENT(in)                      :: n
    INTEGER                                  :: res

    res=n-1
END FUNCTION

! *****************************************************************************
!> \brief max grad for a polynom of the given size
! *****************************************************************************
PURE FUNCTION grad_size2(n) RESULT(res)
    INTEGER, INTENT(in)                      :: n
    INTEGER                                  :: res

    res=INT(FLOOR(0.5_dp *(SQRT(1.0_dp+8.0_dp *REAL(n,dp))-1.0_dp)-2.e-6_dp))
END FUNCTION

! *****************************************************************************
!> \brief max grad for a polynom of the given size
! *****************************************************************************
PURE FUNCTION grad_size3(n) RESULT(res)
    INTEGER, INTENT(in)                      :: n
    INTEGER                                  :: res

    INTEGER                                  :: nn
    REAL(dp)                                 :: g1

    IF (n<1) THEN
        res=-1
    ELSE
        nn=n*6
        g1=(108.0_dp*nn+12.0_dp*SQRT(81.0_dp*nn*nn-12.0_dp))**(1.0_dp/3.0_dp)
        res=FLOOR(g1/6.0_dp+2.0_dp/g1-1.0_dp-2.e-6_dp)
    END IF
END FUNCTION

! *****************************************************************************
!> \brief 0-based index of monomial of the given degree
! *****************************************************************************
PURE FUNCTION mono_index1(i) RESULT(res)
    INTEGER, INTENT(in)                      :: i
    INTEGER                                  :: res

    res=i
END FUNCTION

! *****************************************************************************
!> \brief 0-based index of monomial of the given degree
! *****************************************************************************
PURE FUNCTION mono_index2(i,j) RESULT(res)
    INTEGER, INTENT(in)                      :: i, j
    INTEGER                                  :: res

    INTEGER                                  :: grad

    grad=i+j
    res=grad*(grad+1)/2+j
END FUNCTION

! *****************************************************************************
!> \brief 0-based index of monomial of the given degree
! *****************************************************************************
PURE FUNCTION mono_index3(i,j,k) RESULT(res)
    INTEGER, INTENT(in)                      :: i, j, k
    INTEGER                                  :: res

    INTEGER                                  :: grad, sgrad

    sgrad=j+k
    grad=i+sgrad
    res=grad*(grad+1)*(grad+2)/6+(sgrad)*(sgrad+1)/2+k
END FUNCTION

! *****************************************************************************
!> \brief exponents of the monomial at the given 0-based index
! *****************************************************************************
PURE FUNCTION mono_exp1(ii) RESULT(res)
    INTEGER, INTENT(in)                      :: ii
    INTEGER                                  :: res

    res=ii
END FUNCTION

! *****************************************************************************
!> \brief exponents of the monomial at the given 0-based index
! *****************************************************************************
PURE FUNCTION mono_exp2(ii) RESULT(res)
    INTEGER, INTENT(in)                      :: ii
    INTEGER, DIMENSION(2)                    :: res

    INTEGER                                  :: grad

    grad=INT(FLOOR(0.5_dp*(SQRT(9.0_dp+8.0_dp*ii)-1.0_dp)-2.e-6_dp))
    res(2)=ii-(grad)*(grad+1)/2
    res(1)=grad-res(2)
END FUNCTION

! *****************************************************************************
!> \brief exponents of the monomial at the given 0-based index
! *****************************************************************************
PURE FUNCTION mono_exp3(n) RESULT(res)
    INTEGER, INTENT(in)                      :: n
    INTEGER, DIMENSION(3)                    :: res

    INTEGER                                  :: grad, grad1, ii, nn
    REAL(dp)                                 :: g1

    nn=(n+1)*6
    g1=(108.0_dp*nn+12.0_dp*SQRT(81.0_dp*nn*nn-12.0_dp))**(1.0_dp/3.0_dp)
    grad1=INT(FLOOR(g1/6.0_dp+2.0_dp/g1-1.0_dp-2.e-6_dp))
    ii=n-grad1*(grad1+1)*(grad1+2)/6
    grad=INT(FLOOR(0.5_dp*(SQRT(9.0_dp+8.0_dp*ii)-1.0_dp)-1.e-6_dp))
    res(3)=ii-grad*(grad+1)/2
    res(2)=grad-res(3)
    res(1)=grad1-grad
END FUNCTION

! *****************************************************************************
!> \brief the index of the result of the multiplication of the two monomials
! *****************************************************************************
PURE FUNCTION mono_mult1(ii,ij) RESULT(res)
    INTEGER, INTENT(in)                      :: ii, ij
    INTEGER                                  :: res

    res=ii+ij
END FUNCTION

! *****************************************************************************
!> \brief the index of the result of the multiplication of the two monomials
! *****************************************************************************
PURE FUNCTION mono_mult2(ii,ij) RESULT(res)
    INTEGER, INTENT(in)                      :: ii, ij
    INTEGER                                  :: res

    INTEGER, DIMENSION(2)                    :: monoRes

    monoRes=mono_exp2(ii)+mono_exp2(ij)
    res=mono_index2(monoRes(1),monoRes(2))
END FUNCTION

! *****************************************************************************
!> \brief the index of the result of the multiplication of the two monomials
! *****************************************************************************
PURE FUNCTION mono_mult3(ii,ij) RESULT(res)
    INTEGER, INTENT(in)                      :: ii, ij
    INTEGER                                  :: res

    INTEGER, DIMENSION(3)                    :: monoRes

    monoRes=mono_exp3(ii)+mono_exp3(ij)
    res=mono_index3(monoRes(1),monoRes(2),monoRes(3))
END FUNCTION

! *****************************************************************************
!> \brief multiplies the polynomials p1 with p2 using pRes to store the result
! *****************************************************************************
SUBROUTINE poly_mult1(p1,p2,pRes,np1,sumUp,error)
    REAL(dp), DIMENSION(:), INTENT(in)       :: p1, p2
    REAL(dp), DIMENSION(:), INTENT(inout)    :: pRes
    INTEGER, INTENT(in), OPTIONAL            :: np1
    LOGICAL, INTENT(in), OPTIONAL            :: sumUp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'poly_mult1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ipoly, iPos, j, myNp1, &
                                                newGrad, newSize, resPos, &
                                                resShift_0, size_p1, size_p2
    LOGICAL                                  :: mySumUp

    IF (.not.module_initialized) CALL init_d3_poly_module
    mySumUp=.FALSE.
    myNp1=1
    IF (PRESENT(np1)) myNp1=np1
    IF (PRESENT(sumUp)) mySumUp=sumUp
    size_p1=SIZE(p1)/myNp1
    size_p2=SIZE(p2)
    newGrad=grad_size1(size_p1)+grad_size1(size_p2)
    newSize=SIZE(pRes)/myNp1
    CPPreconditionNoFail(newSize>=poly_size1(newGrad),cp_failure_level,routineP,error)
    IF (.not.mySumUp) pRes=0
    iPos=1
    resShift_0=0
    DO ipoly=0,myNp1-1
        DO i=1,size_p1
            resPos=resShift_0+i
            DO j=1,size_p2
                pRes(resPos)=pRes(resPos)+p1(iPos)*p2(j)
                resPos=resPos+1
            END DO
            iPos=iPos+1
        END DO
        resShift_0=resShift_0+newSize
    END DO
END SUBROUTINE

! *****************************************************************************
!> \brief multiplies p1 with p2 using pRes to store the result
! *****************************************************************************
SUBROUTINE poly_mult2(p1,p2,pRes,np1,sumUp,error)
    REAL(dp), DIMENSION(:), INTENT(in)       :: p1, p2
    REAL(dp), DIMENSION(:), INTENT(inout)    :: pRes
    INTEGER, INTENT(in), OPTIONAL            :: np1
    LOGICAL, INTENT(in), OPTIONAL            :: sumUp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'poly_mult2', &
      routineP = moduleN//':'//routineN

    INTEGER :: g1, g1g2, g2, grad1, grad2, i, ipoly, iShift, j, msize_p1, &
      myNp1, newGrad, newSize, shift1, shift2, shiftRes, shiftRes_0, size_p1, &
      size_p2, subG2
    LOGICAL                                  :: mySumUp

    IF (.not.module_initialized) CALL init_d3_poly_module
    mySumUp=.FALSE.
    IF (PRESENT(sumUp)) mySumUp=sumUp
    myNp1=1
    IF (PRESENT(np1)) myNp1=np1
    size_p1=SIZE(p1)/myNp1
    size_p2=SIZE(p2)
    grad1=grad_size2(size_p1)
    grad2=grad_size2(size_p2)
    newGrad=grad1+grad2
    newSize=SIZE(pRes)/myNp1
    CPPreconditionNoFail(newSize>=poly_size2(newGrad),cp_failure_level,routineP,error)
    IF (.not.mySumUp) pRes=0
    iShift=0
    shiftRes=0
    DO ipoly=1,myNp1
        DO i=1,MIN(size_p1,cached_dim2)
            DO j=1,MIN(size_p2,cached_dim2)
                pRes(shiftRes+a_mono_mult2(j,i))=pRes(shiftRes+a_mono_mult2(j,i))&
                    +p1(iShift+i)*p2(j)
            END DO
        END DO
        iShift=iShift+size_p1
        shiftRes=shiftRes+newSize
    END DO
    IF (grad1>max_grad2.or.grad2>max_grad2) THEN
        msize_p1=size_p1
        shiftRes_0=0
        DO ipoly=0,myNp1-1
            shift1=ipoly*size_p1
            DO g1=0,grad1
                ! shift1=g1*(g1+1)/2
                IF (g1>max_grad2) THEN
                    subG2=0
                    shift2=0
                    g1g2=shiftRes_0-1
                ELSE
                    subG2=max_grad2+1
                    shift2=cached_dim2
                    g1g2=shiftRes_0+g1*subG2-1
                END IF
                DO g2=subG2,grad2
                    ! shift2=g2*(g2+1)/2
                    shiftRes=shift1+shift2+g1g2 ! shiftRes=(g1+g2)*(g1+g2+1)/2-1+ipoly*newSize
                    DO i=1,MIN(g1+1,msize_p1-shift1)
                        DO j=1,MIN(g2+1,size_p2-shift2)
                            pRes(shiftRes+i+j)=pRes(shiftRes+i+j)+p1(shift1+i)*p2(shift2+j)
                        END DO
                    END DO
                    shift2=shift2+g2+1 !
                    g1g2=g1g2+g1 !
                END DO
                shift1=shift1+g1+1 !
            END DO
            shiftRes_0=shiftRes_0+newSize-size_p1
            msize_p1=msize_p1+size_p1
        END DO
    END IF
END SUBROUTINE

! *****************************************************************************
!> \brief multiplies p1 with p2 using pRes to store the result
! *****************************************************************************
SUBROUTINE poly_mult3(p1,p2,pRes,np1,sumUp,error)
    REAL(dp), DIMENSION(:), INTENT(in)       :: p1, p2
    REAL(dp), DIMENSION(:), INTENT(inout)    :: pRes
    INTEGER, INTENT(in), OPTIONAL            :: np1
    LOGICAL, INTENT(in), OPTIONAL            :: sumUp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'poly_mult3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: grad1, grad2, myNp1, newGrad, &
                                                newSize, size_p1, size_p2
    LOGICAL                                  :: mySumUp

    IF (.not.module_initialized) CALL init_d3_poly_module
    mySumUp=.FALSE.
    IF (PRESENT(sumUp)) mySumUp=sumUp
    myNp1=1
    IF (PRESENT(np1)) myNp1=np1
    size_p1=SIZE(p1)/myNp1
    size_p2=SIZE(p2)
    grad1=grad_size3(size_p1)
    grad2=grad_size3(size_p2)
    newGrad=grad1+grad2
    newSize=SIZE(pRes)/myNp1
    CPPreconditionNoFail(newSize>=poly_size3(newGrad),cp_failure_level,routineP,error)
    CALL poly_mult3b(p1,SIZE(p1),grad1,p2,SIZE(p2),grad2,pRes,SIZE(pRes),myNp1,mySumUp)
END SUBROUTINE

! *****************************************************************************
!> \brief low level routine of poly_mult3 without checks
! *****************************************************************************
SUBROUTINE poly_mult3b(p1,size_p1,grad1,p2,size_p2,grad2,pRes,size_pRes,np1,sumUp)
    INTEGER, INTENT(in)                      :: size_p1
    REAL(dp), &
      DIMENSION(IF_CHECK(size_p1, *)), &
      INTENT(in)                             :: p1
    INTEGER, INTENT(in)                      :: grad1, size_p2
    REAL(dp), &
      DIMENSION(IF_CHECK(size_p2, *)), &
      INTENT(in)                             :: p2
    INTEGER, INTENT(in)                      :: grad2, size_pRes
    REAL(dp), &
      DIMENSION(IF_CHECK(size_pRes, *)), &
      INTENT(inout)                          :: pRes
    INTEGER, INTENT(in)                      :: np1
    LOGICAL, INTENT(in)                      :: sumUp

    CHARACTER(len=*), PARAMETER :: routineN = 'poly_mult3b', &
      routineP = moduleN//':'//routineN

    INTEGER :: g1, g2, i, i1, i2, ipoly, iShift, j, j1, j2, msize_p1, &
      my_size_p1, newSize, shift1, shift1I, shift1J, shift2, shift2I, &
      shift2J, shiftRes, shiftRes_0, shiftResI, shiftResI_0, shiftResJ, &
      subG2, subGrad

    IF (.not.module_initialized) CALL init_d3_poly_module
    my_size_p1=size_p1/np1
    newSize=size_pRes/np1

    IF (.not.sumUp) pRes(1:size_pRes)=0.0_dp
    iShift=0
    shiftRes=0
    DO ipoly=1,np1
        DO i=1,MIN(my_size_p1,cached_dim3)
            DO j=1,MIN(size_p2,cached_dim3)
                pRes(shiftRes+a_mono_mult3(j,i))=pRes(shiftRes+a_mono_mult3(j,i))&
                    +p1(iShift+i)*p2(j)
            END DO
        END DO
        iShift=iShift+my_size_p1
        shiftRes=shiftRes+newSize
    END DO
    IF (grad1>max_grad3 .OR. grad2>max_grad3) THEN
        ! one could remove multiplications even more agressively...
        msize_p1=my_size_p1
        DO ipoly=0,np1-1
            shift1=1+ipoly*my_size_p1
            shiftRes_0=1+ipoly*newSize
            DO g1=0,grad1
                IF (g1>max_grad3) THEN
                    subG2=0
                    shift2=1
                ELSE
                    subG2=max_grad3+1
                    shift2=subG2*(subG2+1)*(subG2+2)/6+1
                END IF
                DO g2=subG2,grad2
                    shiftRes=(g1+g2)*(g1+g2+1)*(g1+g2+2)/6+shiftRes_0
                    shift1I=shift1
                    shiftResI_0=shiftRes
                    DO i1=g1,0,-1
                        IF (shift1I>msize_p1) EXIT
                        shift2I=shift2
                        shiftResI=shiftResI_0
                        subGrad=g1-i1
                        DO i2=g2,0,-1
                            !subGrad=g1+g2-i1-i2
                            !shiftResI=shiftRes+(subGrad)*(subGrad+1)/2
                            !shift2I=shift2+(g2-i2)*(g2-i2+1)/2
                            IF (shift2I>size_p2) EXIT
                            DO j1=g1-i1,0,-1
                                shift1J=shift1I+g1-i1-j1
                                IF (shift1J>msize_p1) EXIT
                                DO j2=g2-i2,0,-1
                                    shift2J=shift2I+g2-i2-j2
                                    IF (shift2J>size_p2) EXIT
                                    shiftResJ=shiftResI+(subGrad-j1-j2)
                                    ! shift1J=mono_index3(i1,j1,g1-i1-j1)+ipoly*my_size_p1+1
                                    ! shift2J=mono_index3(i2,j2,g2-i2-j2)+1
                                    ! shiftResJ=mono_index3(i1+i2,j1+j2,g1+g2-i1-i2-j1-j2)+ipoly*newSize+1
                                    pRes(shiftResJ)=pRes(shiftResJ)+p1(shift1J)*p2(shift2J)
                                END DO
                            END DO
                            subGrad=subGrad+1
                            shift2I=shift2I+(g2-i2+1)
                            shiftResI=shiftResI+subGrad
                        END DO
                        shift1I=shift1I+(g1-i1+1)
                        shiftResI_0=shiftResI_0+(g1-i1+1)
                    END DO
                    shift2=shift2+(g2+1)*(g2+2)/2
                END DO
                shift1=shift1+(g1+1)*(g1+2)/2
            END DO
            msize_p1=msize_p1+my_size_p1
        END DO
    END IF
END SUBROUTINE


! *****************************************************************************
!> \brief low level routine that multiplies with a polynomial of grad 1
! *****************************************************************************
SUBROUTINE poly_mult3ab(p1,size_p1,grad1,p2,pRes,size_pRes,np1,sumUp)
    INTEGER, INTENT(in)                      :: size_p1
    REAL(dp), &
      DIMENSION(IF_CHECK(size_p1, *)), &
      INTENT(in)                             :: p1
    INTEGER, INTENT(in)                      :: grad1
    REAL(dp), DIMENSION(IF_CHECK(4, *)), &
      INTENT(in)                             :: p2
    INTEGER, INTENT(in)                      :: size_pRes
    REAL(dp), &
      DIMENSION(IF_CHECK(size_pRes, *)), &
      INTENT(inout)                          :: pRes
    INTEGER, INTENT(in)                      :: np1
    LOGICAL, INTENT(in)                      :: sumUp

    CHARACTER(len=*), PARAMETER :: routineN = 'poly_mult3ab', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: grad2 = 1, size_p2 = 4

    INTEGER :: g1, g2, i, i1, i2, ipoly, iShift, j1, j2, msize_p1, &
      my_size_p1, newSize, shift1, shift1I, shift1J, shift2, shift2I, &
      shift2J, shiftRes, shiftRes_0, shiftResI, shiftResI_0, shiftResJ, &
      subG2, subGrad

    IF (.not.module_initialized) CALL init_d3_poly_module
    my_size_p1=size_p1/np1
    newSize=size_pRes/np1

    IF (.not.sumUp) pRes(1:size_pRes)=0.0_dp
    iShift=0
    shiftRes=0
    DO ipoly=1,np1
        DO i=1,MIN(my_size_p1,cached_dim3)
           pRes(shiftRes+a_mono_mult3a(1,i))=pRes(shiftRes+a_mono_mult3a(1,i))&
             +p1(iShift+i)*p2(1)
           pRes(shiftRes+a_mono_mult3a(2,i))=pRes(shiftRes+a_mono_mult3a(2,i))&
             +p1(iShift+i)*p2(2)
           pRes(shiftRes+a_mono_mult3a(3,i))=pRes(shiftRes+a_mono_mult3a(3,i))&
             +p1(iShift+i)*p2(3)
           pRes(shiftRes+a_mono_mult3a(4,i))=pRes(shiftRes+a_mono_mult3a(4,i))&
             +p1(iShift+i)*p2(4)
        END DO
        iShift=iShift+my_size_p1
        shiftRes=shiftRes+newSize
    END DO
    IF (grad1>max_grad3 .OR. grad2>max_grad3) THEN
        ! one could remove multiplications even more agressively...
        msize_p1=my_size_p1
        DO ipoly=0,np1-1
            shift1=1+ipoly*my_size_p1+(max_grad3+1)*(max_grad3+2)*(max_grad3+3)/6
            shiftRes_0=1+ipoly*newSize
            DO g1=max_grad3+1,grad1
                subG2=0
                shift2=1
                DO g2=subG2,grad2
                    shiftRes=(g1+g2)*(g1+g2+1)*(g1+g2+2)/6+shiftRes_0
                    shift1I=shift1
                    shiftResI_0=shiftRes
                    DO i1=g1,0,-1
                        IF (shift1I>msize_p1) EXIT
                        shift2I=shift2
                        shiftResI=shiftResI_0
                        subGrad=g1-i1
                        DO i2=g2,0,-1
                            !subGrad=g1+g2-i1-i2
                            !shiftResI=shiftRes+(subGrad)*(subGrad+1)/2
                            !shift2I=shift2+(g2-i2)*(g2-i2+1)/2
                            DO j1=g1-i1,0,-1
                                shift1J=shift1I+g1-i1-j1
                                IF (shift1J>msize_p1) EXIT
                                DO j2=g2-i2,0,-1
                                    shift2J=shift2I+g2-i2-j2
                                    shiftResJ=shiftResI+(subGrad-j1-j2)
                                    ! shift1J=mono_index3(i1,j1,g1-i1-j1)+ipoly*my_size_p1+1
                                    ! shift2J=mono_index3(i2,j2,g2-i2-j2)+1
                                    ! shiftResJ=mono_index3(i1+i2,j1+j2,g1+g2-i1-i2-j1-j2)+ipoly*newSize+1
                                    pRes(shiftResJ)=pRes(shiftResJ)+p1(shift1J)*p2(shift2J)
                                END DO
                            END DO
                            subGrad=subGrad+1
                            shift2I=shift2I+(g2-i2+1)
                            shiftResI=shiftResI+subGrad
                        END DO
                        shift1I=shift1I+(g1-i1+1)
                        shiftResI_0=shiftResI_0+(g1-i1+1)
                    END DO
                    shift2=shift2+(g2+1)*(g2+2)/2
                END DO
                shift1=shift1+(g1+1)*(g1+2)/2
            END DO
            msize_p1=msize_p1+my_size_p1
        END DO
    END IF
END SUBROUTINE

! *****************************************************************************
!> \brief writes out a 1d polynomial in a human readable form
! *****************************************************************************
SUBROUTINE poly_write1(p,out_f)
    REAL(dp), DIMENSION(:), INTENT(in)       :: p
    INTEGER, INTENT(in)                      :: out_f

    INTEGER                                  :: i
    LOGICAL                                  :: did_write

    did_write=.FALSE.
    DO i=1,SIZE(p)
        IF (p(i)/=0) THEN
            IF (p(i)>=0) WRITE(out_f,'("+")',advance='NO')
            WRITE(out_f,'(G20.10)',advance='NO') p(i)
            IF (i/=1) WRITE(out_f,'("*x^",I3)',advance='NO') i-1
            did_write=.TRUE.
        END IF
    END DO
    IF (.NOT. did_write) WRITE (out_f,'("0.0")',advance='NO')
END SUBROUTINE

! *****************************************************************************
!> \brief writes out a 2d polynomial in a human readable form
! *****************************************************************************
SUBROUTINE poly_write2(p,out_f)
    REAL(dp), DIMENSION(:), INTENT(in)       :: p
    INTEGER, INTENT(in)                      :: out_f

    INTEGER                                  :: i
    INTEGER, DIMENSION(2)                    :: mono_e
    LOGICAL                                  :: did_write

    IF (.not.module_initialized) CALL init_d3_poly_module
    did_write=.FALSE.
    DO i=1,SIZE(p)
        mono_e=mono_exp2(i-1)
        IF (p(i)/=0) THEN
            IF (p(i)>=0) WRITE(out_f,'("+")',advance='NO')
            WRITE(out_f,'(G20.10)',advance='NO') p(i)
            IF (mono_e(1)/=0) WRITE(out_f,'("*x^",I3)',advance='NO') mono_e(1)
            IF (mono_e(2)/=0) WRITE(out_f,'("*y^",I3)',advance='NO') mono_e(2)
            did_write=.TRUE.
        END IF
    END DO
    IF (.NOT. did_write) WRITE (out_f,'("0.0")',advance='NO')
END SUBROUTINE

! *****************************************************************************
!> \brief writes out the polynomial in a human readable form
! *****************************************************************************
SUBROUTINE poly_write3(p,out_f)
    REAL(dp), DIMENSION(:), INTENT(in)       :: p
    INTEGER, INTENT(in)                      :: out_f

    INTEGER                                  :: i
    INTEGER, DIMENSION(3)                    :: mono_e
    LOGICAL                                  :: did_write

    IF (.not.module_initialized) CALL init_d3_poly_module
    did_write=.FALSE.
    DO i=1,SIZE(p)
        mono_e=mono_exp3(i-1)
        IF (p(i)/=0) THEN
            IF (p(i)>=0) WRITE(out_f,'("+")',advance='NO')
            WRITE(out_f,'(G20.10)',advance='NO') p(i)
            IF (mono_e(1)/=0) WRITE(out_f,'("*x^",I3)',advance='NO') mono_e(1)
            IF (mono_e(2)/=0) WRITE(out_f,'("*y^",I3)',advance='NO') mono_e(2)
            IF (mono_e(3)/=0) WRITE(out_f,'("*z^",I3)',advance='NO') mono_e(3)
            did_write=.TRUE.
        END IF
    END DO
    IF (.NOT. did_write) WRITE (out_f,'("0.0")',advance='NO')
END SUBROUTINE

! *****************************************************************************
!> \brief random poly with coeffiecents that are easy to print exactly,
!>        of the given maximum size (for testing purposes)
! *****************************************************************************
FUNCTION poly_random(p,maxSize,minSize,error) RESULT(res)
    REAL(dp), DIMENSION(:), INTENT(out)      :: p
    INTEGER, INTENT(in)                      :: maxSize
    INTEGER, INTENT(in), OPTIONAL            :: minSize
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'poly_random', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, myMinSize, pSize
    REAL(dp)                                 :: g

    IF (.not.module_initialized) CALL init_d3_poly_module
    myMinSize=1
    IF (PRESENT(minSize)) myMinSize=minSize
    CALL RANDOM_NUMBER(g)
    pSize=MIN(maxSize,myMinSize+INT((maxSize-myMinSize+1)*g))
    CPPreconditionNoFail(SIZE(p)>=pSize,cp_failure_level,routineP,error)
    CALL RANDOM_NUMBER(p)
    DO i=1,pSize
        p(i)=REAL(INT(p(i)*200.0_dp-100.0_dp),dp)/100.0_dp
    END DO
    DO i=pSize+1,SIZE(p)
        p(i)=0.0_dp
    END DO
    res=pSize
END FUNCTION

! *****************************************************************************
!> \brief returns in the polynomials pRes the transpose of the
!> affine transformation x -> m*x+b of p
! *****************************************************************************
SUBROUTINE poly_affine_t3t(p,m,b,pRes,npoly,error)
    REAL(dp), DIMENSION(:), INTENT(in)       :: p
    REAL(dp), DIMENSION(3, 3), INTENT(in)    :: m
    REAL(dp), DIMENSION(3), INTENT(in)       :: b
    REAL(dp), DIMENSION(:), INTENT(out)      :: pRes
    INTEGER, INTENT(in), OPTIONAL            :: npoly
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'poly_affine_t3t', &
      routineP = moduleN//':'//routineN

    INTEGER :: grad, i, igrad, ii, ii1, ipoly, j, k, minResSize, monoDim1, &
      monoDim2, monoDimAtt, monoDimAtt2, monoFullDim1, monoFullDim2, &
      monoSize1, monoSize2, my_npoly, pcoeff, pIdx, pShift, rescoeff, &
      resShift, rest_size_p, size_p, size_res, start_idx1, stat
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: monoG1, monoG2
    REAL(dp), DIMENSION(4, 3)                :: basepoly

    IF (.not.module_initialized) CALL init_d3_poly_module
    my_npoly=1
    IF (PRESENT(npoly)) my_npoly=npoly
    basepoly(1,:)=b
    DO j=1,3
        DO i=1,3
            basepoly(j+1,i)=m(i,j)
        END DO
    END DO
    size_p=SIZE(pRes)/my_npoly
    size_res=SIZE(p)/my_npoly
    grad=grad_size3(size_res)
    minResSize=poly_size3(grad)
    CPPreconditionNoFail(size_res==minResSize,cp_failure_level,routineP,error)
    CPPreconditionNoFail(size_p>=minResSize,cp_failure_level,routineP,error)
    pRes=0
    IF (size_p==0) RETURN
    ii1=1
    ii=1
    DO ipoly=1,my_npoly
        pRes(ii1)=p(ii)
        ii=ii+size_res
        ii1=ii1+size_p
    END DO
    IF (size_p==1) RETURN

    ALLOCATE(monoG1((grad+1)*(grad+2)/2*minResSize),&
        monoG2((grad+1)*(grad+2)/2*minResSize),stat=stat)
    CPPostconditionNoFail(stat==0,cp_fatal_level,routineP,error)
    !monoG1=0
    !monoG2=0
    ii1=1
    DO j=1,3
        DO i=1,4
            monoG1(ii1)=basepoly(i,j)
            ii1=ii1+1
        END DO
    END DO
    ii1=2
    igrad=1
    monoDim1=4
    monoSize1=3
    monoFullDim1=monoDim1*monoSize1
    rest_size_p=size_p-1
    DO
        k=MIN(rest_size_p,monoSize1)
        !call dgemm('T','N',monoDim1,my_npoly,k,&
        !    1.0_dp,monoG1,monoDim1,p(ii1:),size_p,1.0_dp,pRes,size_res)
        resShift=0
        pShift=ii1
        DO ipoly=1,my_npoly
            pIdx=pShift
            ii=1
            DO pcoeff=1,k
                DO rescoeff=1,monoDim1
                    pRes(pIdx)=pRes(pIdx)+p(resShift+rescoeff)*monoG1(ii)
                    ii=ii+1
                END DO
                pIdx=pIdx+1
            END DO
            resShift=resShift+size_res
            pShift=pShift+size_p
        END DO

        rest_size_p=rest_size_p-k
        ii1=ii1+k
        IF (rest_size_p<=0) EXIT

        monoSize2=igrad+2+monoSize1
        monoDim2=monoDim1+monoSize2
        monoFullDim2=monoSize2*monoDim2
        monoDimAtt=monoSize1*monoDim2
        CALL poly_mult3ab(IF_CHECK(monoG1(1:monoFullDim1),monoG1(1)),monoFullDim1,igrad,&
            IF_CHECK(basepoly(:,1),basepoly(1,1)),&
            IF_CHECK(monoG2(1:monoDimAtt),monoG2(1)),monoDimAtt,monoSize1,.FALSE.)
        monoDimAtt2=monoFullDim2-monoDim2
        start_idx1=(monoSize1-igrad-1)*monoDim1
        CALL poly_mult3ab(IF_CHECK(monoG1(start_idx1+1:monoFullDim1),monoG1(start_idx1+1)),&
            monoFullDim1-start_idx1,igrad,IF_CHECK(basepoly(:,2),basepoly(1,2)),&
            IF_CHECK(monoG2(monoDimAtt+1:monoDimAtt2),monoG2(monoDimAtt+1)),&
            monoDimAtt2-monoDimAtt,igrad+1,.FALSE.)
        CALL poly_mult3ab(IF_CHECK(monoG1(monoFullDim1-monoDim1+1:monoFullDim1),monoG1(monoFullDim1-monoDim1+1)),&
            monoDim1,igrad,IF_CHECK(basepoly(:,3),basepoly(1,3)),&
            IF_CHECK(monoG2(monoDimAtt2+1:monoFullDim2),monoG2(monoDimAtt2+1)),&
            monoFullDim2-monoDimAtt2,1,.FALSE.)
        igrad=igrad+1

        ! even grads

        k=MIN(rest_size_p,monoSize2)
        !call dgemm('T','N',monoDim2,my_npoly,k,&
        !    1.0_dp,monoG2,monoDim2,p(ii1:),size_p,1.0_dp,pRes,size_res)
        resShift=0
        pShift=ii1
        DO ipoly=1,my_npoly
            pIdx=pShift
            ii=1
            DO pcoeff=1,k
                DO rescoeff=1,monoDim2
                    pRes(pIdx)=pRes(pIdx)+p(resShift+rescoeff)*monoG2(ii)
                    ii=ii+1
                END DO
                pIdx=pIdx+1
            END DO
            resShift=resShift+size_res
            pShift=pShift+size_p
        END DO

        rest_size_p=rest_size_p-k
        ii1=ii1+k
        IF (rest_size_p<=0) EXIT

        monoSize1=igrad+2+monoSize2
        monoDim1=monoDim2+monoSize1
        monoFullDim1=monoSize1*monoDim1
        monoDimAtt=monoSize2*monoDim1
        CALL poly_mult3ab(IF_CHECK(monoG2(1:monoFullDim2),monoG2(1)),monoFullDim2,igrad,&
            IF_CHECK(basepoly(:,1),basepoly(1,1)),IF_CHECK(monoG1(1:monoDimAtt),monoG1(1)),&
            monoDimAtt,monoSize2,.FALSE.)
        monoDimAtt2=monoFullDim1-monoDim1
        start_idx1=(monoSize2-igrad-1)*monoDim2
        CALL poly_mult3ab(IF_CHECK(monoG2(start_idx1+1:monoFullDim2),monoG2(start_idx1+1)),&
            monoFullDim2-start_idx1,igrad,IF_CHECK(basepoly(:,2),basepoly(1,2)),&
            IF_CHECK(monoG1(monoDimAtt+1:monoDimAtt2),monoG1(monoDimAtt+1)),monoDimAtt2-monoDimAtt,&
            igrad+1,.FALSE.)
        CALL poly_mult3ab(IF_CHECK(monoG2(monoFullDim2-monoDim2+1:monoFullDim2),monoG2(monoFullDim2-monoDim2+1)),&
            monoDim2,igrad,IF_CHECK(basepoly(:,3),basepoly(1,3)),&
            IF_CHECK(monoG1(monoDimAtt2+1:monoFullDim1),monoG1(monoDimAtt2+1)),&
            monoFullDim1-monoDimAtt2,1,.FALSE.)
        igrad=igrad+1

        ! ! alterantive to unrolling
        ! monoG1=monoG2
        ! monoSize1=monoSize2
        ! monoDim1=monoDim2
        ! monoFullDim1=monoFullDim2
    END DO
    DEALLOCATE(monoG1,monoG2,stat=stat)
    CPPostconditionNoFail(stat==0,cp_fatal_level,routineP,error)
END SUBROUTINE

! *****************************************************************************
!> \brief returns in the polynomials pRes the affine transformation x -> m*x+b of p
! *****************************************************************************
SUBROUTINE poly_affine_t3(p,m,b,pRes,npoly,error)
    REAL(dp), DIMENSION(:), INTENT(in)       :: p
    REAL(dp), DIMENSION(3, 3), INTENT(in)    :: m
    REAL(dp), DIMENSION(3), INTENT(in)       :: b
    REAL(dp), DIMENSION(:), INTENT(out)      :: pRes
    INTEGER, INTENT(in), OPTIONAL            :: npoly
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'poly_affine_t3', &
      routineP = moduleN//':'//routineN

    INTEGER :: grad, i, igrad, ii, ii1, ipoly, j, k, minResSize, monoDim1, &
      monoDim2, monoDimAtt, monoDimAtt2, monoFullDim1, monoFullDim2, &
      monoSize1, monoSize2, my_npoly, pcoeff, pIdx, pShift, rescoeff, &
      resShift, rest_size_p, size_p, size_res, start_idx1, stat
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: monoG1, monoG2
    REAL(dp), DIMENSION(4, 3)                :: basepoly

    IF (.not.module_initialized) CALL init_d3_poly_module
    my_npoly=1
    IF (PRESENT(npoly)) my_npoly=npoly
    basepoly(1,:)=b
    DO j=1,3
        DO i=1,3
            basepoly(j+1,i)=m(i,j)
        END DO
    END DO
    size_p=SIZE(p)/my_npoly
    grad=grad_size3(size_p)
    size_res=SIZE(pRes)/my_npoly
    minResSize=poly_size3(grad)
    CPPreconditionNoFail(size_res>=minResSize,cp_failure_level,routineP,error)
    pRes=0
    IF (size_p==0) RETURN
    ii1=1
    ii=1
    DO ipoly=1,my_npoly
        pRes(ii)=p(ii1)
        ii=ii+size_res
        ii1=ii1+size_p
    END DO
    IF (size_p==1) RETURN

    ALLOCATE(monoG1((grad+1)*(grad+2)/2*minResSize),&
        monoG2((grad+1)*(grad+2)/2*minResSize),stat=stat)
    CPPostconditionNoFail(stat==0,cp_fatal_level,routineP,error)
    monoG1=0
    monoG2=0
    ii1=1
    DO j=1,3
        DO i=1,4
            monoG1(ii1)=basepoly(i,j)
            ii1=ii1+1
        END DO
    END DO
    ii1=2
    igrad=1
    monoDim1=4
    monoSize1=3
    monoFullDim1=monoDim1*monoSize1
    rest_size_p=size_p-1
    DO
        k=MIN(rest_size_p,monoSize1)
        !call dgemm('T','N',monoDim1,my_npoly,k,&
        !    1.0_dp,monoG1,monoDim1,p(ii1:),size_p,1.0_dp,pRes,size_res)
        resShift=0
        pShift=ii1
        DO ipoly=1,my_npoly
            pIdx=pShift
            ii=1
            DO pcoeff=1,k
                DO rescoeff=1,monoDim1
                    pRes(resShift+rescoeff)=pRes(resShift+rescoeff)+p(pIdx)*monoG1(ii)
                    ii=ii+1
                END DO
                pIdx=pIdx+1
            END DO
            resShift=resShift+size_res
            pShift=pShift+size_p
        END DO

        rest_size_p=rest_size_p-k
        ii1=ii1+k
        IF (rest_size_p<=0) EXIT

        monoSize2=igrad+2+monoSize1
        monoDim2=monoDim1+monoSize2
        monoFullDim2=monoSize2*monoDim2
        monoDimAtt=monoSize1*monoDim2
        CALL poly_mult3ab(IF_CHECK(monoG1(1:monoFullDim1),monoG1(1)),monoFullDim1,igrad,&
            IF_CHECK(basepoly(:,1),basepoly(1,1)),&
            IF_CHECK(monoG2(1:monoDimAtt),monoG2(1)),monoDimAtt,monoSize1,.FALSE.)
        monoDimAtt2=monoFullDim2-monoDim2
        start_idx1=(monoSize1-igrad-1)*monoDim1
        CALL poly_mult3ab(IF_CHECK(monoG1(start_idx1+1:monoFullDim1),monoG1(start_idx1+1)),&
            monoFullDim1-start_idx1,igrad,IF_CHECK(basepoly(:,2),basepoly(1,2)),&
            IF_CHECK(monoG2(monoDimAtt+1:monoDimAtt2),monoG2(monoDimAtt+1)),&
            monoDimAtt2-monoDimAtt,igrad+1,.FALSE.)
        CALL poly_mult3ab(IF_CHECK(monoG1(monoFullDim1-monoDim1+1:monoFullDim1),monoG1(monoFullDim1-monoDim1+1)),&
            monoDim1,igrad,IF_CHECK(basepoly(:,3),basepoly(1,3)),&
            IF_CHECK(monoG2(monoDimAtt2+1:monoFullDim2),monoG2(monoDimAtt2+1)),&
            monoFullDim2-monoDimAtt2,1,.FALSE.)
        igrad=igrad+1

        ! even grads

        k=MIN(rest_size_p,monoSize2)
        !call dgemm('T','N',monoDim2,my_npoly,k,&
        !    1.0_dp,monoG2,monoDim2,p(ii1:),size_p,1.0_dp,pRes,size_res)
        resShift=0
        pShift=ii1
        DO ipoly=1,my_npoly
            pIdx=pShift
            ii=1
            DO pcoeff=1,k
                DO rescoeff=1,monoDim2
                    pRes(resShift+rescoeff)=pRes(resShift+rescoeff)+p(pIdx)*monoG2(ii)
                    ii=ii+1
                END DO
                pIdx=pIdx+1
            END DO
            resShift=resShift+size_res
            pShift=pShift+size_p
        END DO

        rest_size_p=rest_size_p-k
        ii1=ii1+k
        IF (rest_size_p<=0) EXIT

        monoSize1=igrad+2+monoSize2
        monoDim1=monoDim2+monoSize1
        monoFullDim1=monoSize1*monoDim1
        monoDimAtt=monoSize2*monoDim1
        CALL poly_mult3ab(IF_CHECK(monoG2(1:monoFullDim2),monoG2(1)),monoFullDim2,igrad,&
            IF_CHECK(basepoly(:,1),basepoly(1,1)),&
            IF_CHECK(monoG1(1:monoDimAtt),monoG1(1)),monoDimAtt,monoSize2,.FALSE.)
        monoDimAtt2=monoFullDim1-monoDim1
        start_idx1=(monoSize2-igrad-1)*monoDim2
        CALL poly_mult3ab(IF_CHECK(monoG2(start_idx1+1:monoFullDim2),monoG2(start_idx1+1)),&
            monoFullDim2-start_idx1,igrad,&
            IF_CHECK(basepoly(:,2),basepoly(1,2)),&
            IF_CHECK(monoG1(monoDimAtt+1:monoDimAtt2),monoG1(monoDimAtt+1)),monoDimAtt2-monoDimAtt,&
            igrad+1,.FALSE.)
        CALL poly_mult3ab(IF_CHECK(monoG2(monoFullDim2-monoDim2+1:monoFullDim2),monoG2(monoFullDim2-monoDim2+1)),&
            monoDim2,igrad,IF_CHECK(basepoly(:,3),basepoly(1,3)),&
            IF_CHECK(monoG1(monoDimAtt2+1:monoFullDim1),monoG1(monoDimAtt2+1)),&
            monoFullDim1-monoDimAtt2,1,.FALSE.)
        igrad=igrad+1

        ! ! alterantive to unrolling
        ! monoG1=monoG2
        ! monoSize1=monoSize2
        ! monoDim1=monoDim2
        ! monoFullDim1=monoFullDim2
    END DO
    DEALLOCATE(monoG1,monoG2,stat=stat)
    CPPostconditionNoFail(stat==0,cp_fatal_level,routineP,error)
END SUBROUTINE

! *****************************************************************************
!> \brief evaluates the 3d polymial at x (the result is a polynomial in two variables)
! *****************************************************************************
SUBROUTINE poly_p_eval3(p,x,pRes,npoly,error)
    REAL(dp), DIMENSION(:), INTENT(in)       :: p
    REAL(dp), INTENT(in)                     :: x
    REAL(dp), DIMENSION(:), INTENT(inout)    :: pRes
    INTEGER, INTENT(in), OPTIONAL            :: npoly
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'poly_p_eval3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: grad, my_npoly, newSize, &
                                                size_p, stat
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: xi

    IF (.not.module_initialized) CALL init_d3_poly_module
    my_npoly=1
    IF (PRESENT(npoly)) my_npoly=npoly
    size_p=SIZE(p)/my_npoly
    grad=grad_size3(size_p)
    newSize=SIZE(pRes)/my_npoly
    CPPreconditionNoFail(newSize>=poly_size2(grad),cp_failure_level,routineP,error)
    pRes=0.0
    ALLOCATE(xi(grad+1),stat=stat)
    CPPostconditionNoFail(stat==0,cp_fatal_level,routineP,error)
    CALL poly_p_eval3b(p,SIZE(p),x,pRes,SIZE(pRes),my_npoly,grad,xi)
    DEALLOCATE(xi,stat=stat)
    CPPostconditionNoFail(stat==0,cp_fatal_level,routineP,error)
END SUBROUTINE

! *****************************************************************************
!> \brief low level routine of poly_p_eval3 without checks
! *****************************************************************************
SUBROUTINE poly_p_eval3b(p,size_p,x,pRes,size_pRes,npoly,grad,xi)
    INTEGER, INTENT(in)                      :: size_p
    REAL(dp), &
      DIMENSION(IF_CHECK(size_p, *)), &
      INTENT(in)                             :: p
    REAL(dp), INTENT(in)                     :: x
    INTEGER, INTENT(in)                      :: size_pRes
    REAL(dp), &
      DIMENSION(IF_CHECK(size_pRes, *)), &
      INTENT(inout)                          :: pRes
    INTEGER, INTENT(in)                      :: npoly, grad
    REAL(dp), &
      DIMENSION(IF_CHECK(grad+1, *)), &
      INTENT(inout)                          :: xi

    CHARACTER(len=*), PARAMETER :: routineN = 'poly_p_eval3b', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, igrad, ii, ii0, inSize, &
                                                ipoly, j, msize_p, newSize, &
                                                pShift, shiftRes, shiftRes_0, &
                                                subG

    IF (.not.module_initialized) CALL init_d3_poly_module
    inSize=size_p/npoly
    newSize=size_pRes/npoly
    pRes(1:size_pRes)=0.0
    xi(1)=1.0
    DO i=1,grad
        xi(i+1)=xi(i)*x
    END DO
    shiftRes=0
    pShift=0
    DO ipoly=1,npoly
        DO ii=1,MIN(inSize,cached_dim3)
            pRes(shiftRes+a_reduce_idx3(ii))=pRes(shiftRes+a_reduce_idx3(ii))+p(pShift+ii)*xi(a_mono_exp3(1,ii)+1)
        END DO
        shiftRes=shiftRes+newSize
        pShift=pShift+inSize
    END DO
    IF (grad>max_grad3) THEN
        ii0=(max_grad3+1)*(max_grad3+2)*(max_grad3+3)/6+1
        shiftRes_0=1
        msize_p=inSize
        DO ipoly=1,npoly
            ii=ii0
grad_do:    DO igrad=max_grad3+1,grad
                !ii=igrad*(igrad+1)*(igrad+2)/6+1
                shiftRes=shiftRes_0
                subG=0
                DO i=igrad,0,-1
                    !subG=igrad-i
                    !shiftRes=subG*(subG+3)/2+1
                    DO j=subG,0,-1
                        IF (msize_p<ii) EXIT grad_do
                        pRes(shiftRes-j)=pRes(shiftRes-j)+p(ii)*xi(i+1)
                        ii=ii+1
                    END DO
                    shiftRes=shiftRes+subG+2
                    subG=subG+1
                END DO
            END DO grad_do
            ii0=ii0+inSize
            shiftRes_0=shiftRes_0+newSize
            msize_p=msize_p+inSize
        END DO
    END IF
END SUBROUTINE

! *****************************************************************************
!> \brief unevaluates a 2d polymial to a 3d polynomial at x
!>  p(a,b,c)=p(a,b,c)+sum(pRes(b,c)*(x*a)^i,i), this is *not* the inverse of poly_p_eval3
!>  adds to p
! *****************************************************************************
SUBROUTINE poly_padd_uneval3(p,x,pRes,npoly,error)
    REAL(dp), DIMENSION(:), INTENT(inout)    :: p
    REAL(dp), INTENT(in)                     :: x
    REAL(dp), DIMENSION(:), INTENT(in)       :: pRes
    INTEGER, INTENT(in), OPTIONAL            :: npoly
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'poly_padd_uneval3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: grad, my_npoly, newSize, &
                                                size_p, stat
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: xi

    IF (.not.module_initialized) CALL init_d3_poly_module
    my_npoly=1
    IF (PRESENT(npoly)) my_npoly=npoly
    size_p=SIZE(p)/my_npoly
    newSize=SIZE(pRes)/my_npoly
    grad=grad_size2(newSize)
    CPPreconditionNoFail(size_p>=poly_size3(grad),cp_failure_level,routineP,error)
    CPPreconditionNoFail(newSize==poly_size2(grad),cp_failure_level,routineP,error)
    ALLOCATE(xi(grad+1),stat=stat)
    CPPostconditionNoFail(stat==0,cp_fatal_level,routineP,error)
    CALL poly_padd_uneval3b(p,SIZE(p),x,pRes,SIZE(pRes),my_npoly,grad,xi)
    DEALLOCATE(xi,stat=stat)
    CPPostconditionNoFail(stat==0,cp_fatal_level,routineP,error)
END SUBROUTINE

! *****************************************************************************
!> \brief low level routine of poly_padd_uneval3 without checks
!> \note loop should be structured differently (more contiguous pRes access)
! *****************************************************************************
SUBROUTINE poly_padd_uneval3b(p,size_p,x,pRes,size_pRes,npoly,grad,xi)
    INTEGER, INTENT(in)                      :: size_p
    REAL(dp), &
      DIMENSION(IF_CHECK(size_p, *)), &
      INTENT(inout)                          :: p
    REAL(dp), INTENT(in)                     :: x
    INTEGER, INTENT(in)                      :: size_pRes
    REAL(dp), &
      DIMENSION(IF_CHECK(size_pRes, *)), &
      INTENT(in)                             :: pRes
    INTEGER, INTENT(in)                      :: npoly, grad
    REAL(dp), &
      DIMENSION(IF_CHECK(grad+1, *)), &
      INTENT(inout)                          :: xi

    CHARACTER(len=*), PARAMETER :: routineN = 'poly_padd_uneval3b', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, igrad, ii, ii0, inSize, &
                                                ipoly, j, msize_p, newSize, &
                                                pShift, shiftRes, shiftRes_0, &
                                                subG, upSize

    IF (.not.module_initialized) CALL init_d3_poly_module
    inSize=size_p/npoly
    newSize=size_pRes/npoly
    upSize=(grad+1)*(grad+2)*(grad+3)/6
    xi(1)=1.0
    DO i=1,grad
        xi(i+1)=xi(i)*x
    END DO
    shiftRes=0
    pShift=0
    DO ipoly=1,npoly
        DO ii=1,MIN(upSize,cached_dim3)
            p(pShift+ii)=p(pShift+ii)+pRes(shiftRes+a_reduce_idx3(ii))*xi(a_mono_exp3(1,ii)+1)
        END DO
        shiftRes=shiftRes+newSize
        pShift=pShift+inSize
    END DO
    IF (grad>max_grad3) THEN
        ii0=(max_grad3+1)*(max_grad3+2)*(max_grad3+3)/6+1
        shiftRes_0=1
        msize_p=upSize
        DO ipoly=1,npoly
            ii=ii0
grad_do:    DO igrad=max_grad3+1,grad
                !ii=igrad*(igrad+1)*(igrad+2)/6+1
                shiftRes=shiftRes_0
                subG=0
                DO i=igrad,0,-1
                    !subG=igrad-i
                    !shiftRes=subG*(subG+3)/2+1
                    DO j=subG,0,-1
                        IF (msize_p<ii) EXIT grad_do
                        p(ii)=p(ii)+pRes(shiftRes-j)*xi(i+1)
                        ii=ii+1
                    END DO
                    shiftRes=shiftRes+subG+2
                    subG=subG+1
                END DO
            END DO grad_do
            ii0=ii0+inSize
            shiftRes_0=shiftRes_0+newSize
            msize_p=msize_p+inSize
        END DO
    END IF
END SUBROUTINE

! *****************************************************************************
!> \brief evaluates the 2d polynomial at x (the result is a polynomial in one variable)
! *****************************************************************************
SUBROUTINE poly_p_eval2(p,x,pRes,npoly,error)
    REAL(dp), DIMENSION(:), INTENT(in)       :: p
    REAL(dp), INTENT(in)                     :: x
    REAL(dp), DIMENSION(:), INTENT(inout)    :: pRes
    INTEGER, INTENT(in), OPTIONAL            :: npoly
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'poly_p_eval2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: grad, my_npoly, newSize, &
                                                size_p, stat
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: xi

    IF (.not.module_initialized) CALL init_d3_poly_module
    my_npoly=1
    IF (PRESENT(npoly)) my_npoly=npoly
    size_p=SIZE(p)/my_npoly
    grad=grad_size2(size_p)
    newSize=SIZE(pRes)/my_npoly
    pRes=0.0_dp
    CPPreconditionNoFail(newSize>=poly_size1(grad),cp_failure_level,routineP,error)
    ALLOCATE(xi(grad+1),stat=stat)
    CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
    CALL poly_p_eval2b(p,SIZE(p),x,pRes,SIZE(pRes),my_npoly,grad,xi)
    DEALLOCATE(xi,stat=stat)
    CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
END SUBROUTINE

! *****************************************************************************
!> \brief low level routine of poly_p_eval2 without checks
! *****************************************************************************
SUBROUTINE poly_p_eval2b(p,size_p,x,pRes,size_pRes,npoly,grad,xi)
    INTEGER, INTENT(in)                      :: size_p
    REAL(dp), &
      DIMENSION(IF_CHECK(size_p, *)), &
      INTENT(in)                             :: p
    REAL(dp), INTENT(in)                     :: x
    INTEGER, INTENT(in)                      :: size_pRes
    REAL(dp), &
      DIMENSION(IF_CHECK(size_pRes, *)), &
      INTENT(inout)                          :: pRes
    INTEGER, INTENT(in)                      :: npoly
    INTEGER                                  :: grad
    REAL(dp), DIMENSION(IF_CHECK(grad+1, *)) :: xi

    CHARACTER(len=*), PARAMETER :: routineN = 'poly_p_eval2b', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, igrad, ii, ii0, ij, &
                                                inSize, ipoly, msize_p, &
                                                newSize, pShift, shiftRes

    IF (.not.module_initialized) CALL init_d3_poly_module
    inSize=size_p/npoly
    newSize=size_pRes/npoly
    pRes(1:size_pRes)=0.0_dp
    !CPPreconditionNoFail(newSize>grad,cp_failure_level,routineP,error)
    xi(1)=1.0_dp
    DO i=1,grad
        xi(i+1)=xi(i)*x
    END DO
    shiftRes=1
    pShift=0
    DO ipoly=1,npoly
        DO ii=1,MIN(inSize,cached_dim2)
            pRes(shiftRes+a_mono_exp2(2,ii))=pRes(shiftRes+a_mono_exp2(2,ii))+p(pShift+ii)*xi(a_mono_exp2(1,ii)+1)
        END DO
        shiftRes=shiftRes+newSize
        pShift=pShift+inSize
    END DO
    IF (grad>max_grad2) THEN
        ii0=(max_grad2+1)*(max_grad2+2)/2+1
        shiftRes=1
        msize_p=inSize
        DO ipoly=1,npoly
            ii=ii0
grad_do2:   DO igrad=max_grad2+1,grad
                !ii=igrad*(igrad+1)/2+1
                ij=shiftRes
                DO i=igrad,0,-1
                    IF (msize_p<ii) EXIT grad_do2
                    ! ij=igrad-i
                    pRes(ij)=pRes(ij)+p(ii)*xi(i+1)
                    ii=ii+1
                    ij=ij+1
                END DO
            END DO grad_do2
            msize_p=msize_p+inSize
            shiftRes=shiftRes+newSize
            ii0=ii0+inSize
        END DO
    END IF
END SUBROUTINE

! *****************************************************************************
!> \brief unevaluates a 1d polynomial to 2d at x
!>  p(a,b)=sum(pRes(b)*(x*a)^i,i), this is *not* the inverse of poly_p_eval2
!>  overwrites p
! *****************************************************************************
SUBROUTINE poly_padd_uneval2(p,x,pRes,npoly,error)
    REAL(dp), DIMENSION(:), INTENT(inout)    :: p
    REAL(dp), INTENT(in)                     :: x
    REAL(dp), DIMENSION(:), INTENT(in)       :: pRes
    INTEGER, INTENT(in), OPTIONAL            :: npoly
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'poly_padd_uneval2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: grad, my_npoly, newSize, &
                                                size_p, stat
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: xi

    IF (.not.module_initialized) CALL init_d3_poly_module
    my_npoly=1
    IF (PRESENT(npoly)) my_npoly=npoly
    size_p=SIZE(p)/my_npoly
    newSize=SIZE(pRes)/my_npoly
    grad=grad_size1(newSize)
    CPPreconditionNoFail(size_p>=poly_size2(grad),cp_failure_level,routineP,error)
    CPPreconditionNoFail(newSize==poly_size1(grad),cp_failure_level,routineP,error)
    ALLOCATE(xi(grad+1),stat=stat)
    CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
    CALL poly_padd_uneval2b(p,SIZE(p),x,pRes,SIZE(pRes),my_npoly,grad,xi)
    DEALLOCATE(xi,stat=stat)
    CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
END SUBROUTINE

! *****************************************************************************
!> \brief low level routine of poly_p_uneval2 without checks
! *****************************************************************************
SUBROUTINE poly_padd_uneval2b(p,size_p,x,pRes,size_pRes,npoly,grad,xi)
    INTEGER, INTENT(in)                      :: size_p
    REAL(dp), &
      DIMENSION(IF_CHECK(size_p, *)), &
      INTENT(inout)                          :: p
    REAL(dp), INTENT(in)                     :: x
    INTEGER, INTENT(in)                      :: size_pRes
    REAL(dp), &
      DIMENSION(IF_CHECK(size_pRes, *)), &
      INTENT(in)                             :: pRes
    INTEGER, INTENT(in)                      :: npoly
    INTEGER                                  :: grad
    REAL(dp), DIMENSION(IF_CHECK(grad+1, *)) :: xi

    CHARACTER(len=*), PARAMETER :: routineN = 'poly_padd_uneval2b', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, igrad, ii, ii0, ij, &
                                                inSize, ipoly, msize_p, &
                                                newSize, pShift, shiftRes, &
                                                upSize

    IF (.not.module_initialized) CALL init_d3_poly_module
    inSize=size_p/npoly
    upSize=(grad+1)*(grad+2)/2
    newSize=size_pRes/npoly
    !CPPreconditionNoFail(newSize>grad,cp_failure_level,routineP,error)
    xi(1)=1.0_dp
    DO i=1,grad
        xi(i+1)=xi(i)*x
    END DO
    shiftRes=1
    pShift=0
    DO ipoly=1,npoly
        DO ii=1,MIN(upSize,cached_dim2)
            p(pShift+ii)=p(pShift+ii)+pRes(shiftRes+a_mono_exp2(2,ii))*xi(a_mono_exp2(1,ii)+1)
        END DO
        shiftRes=shiftRes+newSize
        pShift=pShift+inSize
    END DO
    IF (grad>max_grad2) THEN
        ii0=(max_grad2+1)*(max_grad2+2)/2+1
        shiftRes=1
        msize_p=upSize
        DO ipoly=1,npoly
            ii=ii0
grad_do2:   DO igrad=max_grad2+1,grad
                !ii=igrad*(igrad+1)/2+1
                ij=shiftRes
                DO i=igrad,0,-1
                    IF (msize_p<ii) EXIT grad_do2
                    ! ij=igrad-i
                    p(ii)=p(ii)+pRes(ij)*xi(i+1)
                    ii=ii+1
                    ij=ij+1
                END DO
            END DO grad_do2
            msize_p=msize_p+inSize
            shiftRes=shiftRes+newSize
            ii0=ii0+inSize
        END DO
    END IF
END SUBROUTINE

! *****************************************************************************
!> \brief evaluates the 1d polynomial at the given place, results are stored contiguosly
! *****************************************************************************
SUBROUTINE poly_eval1(p,x,pRes,npoly,error)
    REAL(dp), DIMENSION(:), INTENT(in)       :: p
    REAL(dp), INTENT(in)                     :: x
    REAL(dp), DIMENSION(:), INTENT(inout)    :: pRes
    INTEGER, INTENT(in), OPTIONAL            :: npoly
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'poly_eval1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ipoly, my_npoly, pShift, &
                                                size_p
    REAL(dp)                                 :: vv, xx

    my_npoly=1
    IF (PRESENT(npoly)) my_npoly=npoly
    size_p=SIZE(p)/my_npoly
    CPPreconditionNoFail(SIZE(pRes)>=my_npoly,cp_failure_level,routineP,error)
    pShift=0
    DO ipoly=1,my_npoly
        xx=1.0_dp
        vv=0.0_dp
        DO i=1,size_p
            vv=vv+p(pShift+i)*xx
            xx=xx*x
        END DO
        pRes(ipoly)=vv
        pShift=pShift+size_p
    END DO
END SUBROUTINE

! *****************************************************************************
!> \brief evaluates the 2d polynomial at the given place, results are stored contiguosly
! *****************************************************************************
SUBROUTINE poly_eval2(p,x,y,pRes,npoly,error)
    REAL(dp), DIMENSION(:), INTENT(in)       :: p
    REAL(dp), INTENT(in)                     :: x, y
    REAL(dp), DIMENSION(:), INTENT(inout)    :: pRes
    INTEGER, INTENT(in), OPTIONAL            :: npoly
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'poly_eval2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: grad, i, igrad, ii, ipoly, j, &
                                                msize_p, my_npoly, pShift, &
                                                size_p, stat
    REAL(dp)                                 :: v
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: xi, yi

    IF (.not.module_initialized) CALL init_d3_poly_module
    my_npoly=1
    IF (PRESENT(npoly)) my_npoly=npoly
    size_p=SIZE(p)/my_npoly
    grad=grad_size2(size_p)
    CPPreconditionNoFail(SIZE(pRes)>=my_npoly,cp_failure_level,routineP,error)
    ALLOCATE(xi(grad+1),yi(grad+1),stat=stat)
    CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
    xi(1)=1.0_dp
    DO i=1,grad
        xi(i+1)=xi(i)*x
    END DO
    yi(1)=1.0_dp
    DO i=1,grad
        yi(i+1)=yi(i)*y
    END DO
    pShift=0
    DO ipoly=1,my_npoly
        v=0.0_dp
        DO ii=1,MIN(size_p,cached_dim2)
            v=v+p(pShift+ii)*xi(a_mono_exp2(1,ii)+1)*yi(a_mono_exp2(2,ii)+1)
        END DO
        pRes(ipoly)=v
        pShift=pShift+size_p
    END DO
    IF (grad>max_grad2) THEN
        pShift=(max_grad2+1)*(max_grad2+2)/2+1
        msize_p=size_p
        DO ipoly=1,my_npoly
            ii=pShift
            v=0.0_dp
grad_do4:   DO igrad=max_grad2+1,grad
                ! ii=igrad*(igrad+1)*(igrad+2)/6+1
                j=1
                DO i=igrad,0,-1
                    IF (msize_p<ii) EXIT grad_do4
                    v=v+p(ii)*xi(i+1)*yi(j)
                    j=j+1
                    ii=ii+1
                END DO
            END DO grad_do4
            pRes(ipoly)=pRes(ipoly)+v
            pShift=pShift+size_p
            msize_p=msize_p+size_p
        END DO
    END IF
END SUBROUTINE

! *****************************************************************************
!> \brief evaluates the 3d polynomial at the given place, results are stored contiguosly
! *****************************************************************************
SUBROUTINE poly_eval3(p,x,y,z,pRes,npoly,error)
    REAL(dp), DIMENSION(:), INTENT(in)       :: p
    REAL(dp), INTENT(in)                     :: x, y, z
    REAL(dp), DIMENSION(:), INTENT(inout)    :: pRes
    INTEGER, INTENT(in), OPTIONAL            :: npoly
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'poly_eval3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: grad, i, igrad, ii, ipoly, j, &
                                                k, msize_p, my_npoly, pShift, &
                                                size_p, stat
    REAL(dp)                                 :: v
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: xi, yi, zi

    IF (.not.module_initialized) CALL init_d3_poly_module
    my_npoly=1
    IF (PRESENT(npoly)) my_npoly=npoly
    size_p=SIZE(p)/my_npoly
    grad=grad_size3(size_p)
    CPPreconditionNoFail(SIZE(pRes)>=my_npoly,cp_failure_level,routineP,error)
    ALLOCATE(xi(grad+1),yi(grad+1),zi(grad+1),stat=stat)
    CPPostconditionNoFail(stat==0,cp_fatal_level,routineP,error)
    xi(1)=1.0_dp
    DO i=1,grad
        xi(i+1)=xi(i)*x
    END DO
    yi(1)=1.0_dp
    DO i=1,grad
        yi(i+1)=yi(i)*y
    END DO
    zi(1)=1.0_dp
    DO i=1,grad
        zi(i+1)=zi(i)*z
    END DO
    pShift=0
    DO ipoly=1,my_npoly
        v=0.0_dp
        DO ii=1,MIN(size_p,cached_dim3)
            v=v+p(pShift+ii)*xi(a_mono_exp3(1,ii)+1)*yi(a_mono_exp3(2,ii)+1)&
                *zi(a_mono_exp3(3,ii)+1)
        END DO
        pRes(ipoly)=v
        pShift=pShift+size_p
    END DO
    IF (grad>max_grad3) THEN
        pShift=(max_grad3+1)*(max_grad3+2)*(max_grad3+3)/6+1
        msize_p=size_p
        DO ipoly=1,my_npoly
            ii=pShift
            v=0.0_dp
grad_do3:   DO igrad=max_grad3+1,grad
                ! ii=igrad*(igrad+1)*(igrad+2)/6+1
                DO i=igrad,0,-1
                    k=1
                    DO j=igrad-i,0,-1
                        ii=(ipoly-1)*size_p+mono_index3(i,j,igrad-i-j)+1
                        IF (msize_p<ii) EXIT grad_do3
                        v=v+p(ii)*xi(i+1)*yi(j+1)*zi(k)
                        k=k+1
                        ii=ii+1
                    END DO
                END DO
            END DO grad_do3
            pRes(ipoly)=pRes(ipoly)+v
            pShift=pShift+size_p
            msize_p=msize_p+size_p
        END DO
    END IF
END SUBROUTINE

! *****************************************************************************
!> \brief returns an array with all dp/dx, the all dp/dy, and finally all dp/dz
! *****************************************************************************
SUBROUTINE poly_derive3(p,pRes,npoly,sumUp,error)
    REAL(dp), DIMENSION(:), INTENT(in)       :: p
    REAL(dp), DIMENSION(:), INTENT(inout)    :: pRes
    INTEGER, INTENT(in), OPTIONAL            :: npoly
    LOGICAL, INTENT(in), OPTIONAL            :: sumUp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'poly_derive3', &
      routineP = moduleN//':'//routineN

    INTEGER :: grad, i, igrad, ii, ii2, ipoly, j, k, msize_p, my_npoly, &
      newSize, pShift, size_p, xDerivShift, yDerivShift, yShift, zDerivShift, &
      zShift
    LOGICAL                                  :: my_sumUp

    IF (.not.module_initialized) CALL init_d3_poly_module
    my_npoly=1
    IF (PRESENT(npoly)) my_npoly=npoly
    my_sumUp=.FALSE.
    IF (PRESENT(sumUp)) my_sumUp=sumUp
    size_p=SIZE(p)/my_npoly
    newSize=SIZE(pRes)/(3*my_npoly)
    grad=grad_size3(size_p)
    CPPreconditionNoFail(newSize>=poly_size3(grad),cp_failure_level,routineP,error)
    IF (.NOT. my_sumUp) pRes=0
    xDerivShift=1
    yDerivShift=my_npoly*newSize+1
    zDerivShift=2*yDerivShift-1
    pShift=0
    DO ipoly=1,my_npoly
        ! split derivs to have less streams to follow (3 vs 5)?
        DO ii=1,MIN(cached_dim3,size_p)
            pRes(xDerivShift+a_deriv_idx3(1,ii))=pRes(xDerivShift+a_deriv_idx3(1,ii))&
                +p(pShift+ii)*a_mono_exp3(1,ii)
            pRes(yDerivShift+a_deriv_idx3(2,ii))=pRes(yDerivShift+a_deriv_idx3(2,ii))&
                +p(pShift+ii)*a_mono_exp3(2,ii)
            pRes(zDerivShift+a_deriv_idx3(3,ii))=pRes(zDerivShift+a_deriv_idx3(3,ii))&
                +p(pShift+ii)*a_mono_exp3(3,ii)
        END DO
        xDerivShift=xDerivShift+newSize
        yDerivShift=yDerivShift+newSize
        zDerivShift=zDerivShift+newSize
        pShift=pShift+size_p
    END DO
    IF (grad>max_grad3) THEN
        xDerivShift=0
        yDerivShift=my_npoly*newSize
        zDerivShift=2*yDerivShift-1
        msize_p=size_p
        xDerivShift=max_grad3*(max_grad3+1)*(max_grad3+2)/6+1
        pShift=xDerivShift+(max_grad3+1)*(max_grad3+2)/2
        DO ipoly=1,my_npoly
            ii=pShift
            ii2=xDerivShift
grad_do5:   DO igrad=max_grad3+1,grad
                yShift=yDerivShift
                zShift=zDerivShift
                DO i=igrad,0,-1
                    k=0
                    DO j=igrad-i,0,-1
                        IF (ii>msize_p) EXIT grad_do5
                        ! remove ifs?
                        IF (i>0) pRes(ii2)=pRes(ii2)+p(ii)*i
                        IF (j>0) pRes(yShift+ii2)=pRes(yShift+ii2)+p(ii)*j
                        IF (k>0) pRes(zShift+ii2)=pRes(zShift+ii2)+k*p(ii)
                        ii=ii+1
                        ii2=ii2+1
                        k=k+1
                    END DO
                    yShift=yShift-1
                    zShift=zShift-1
                END DO
                ii2=ii2-igrad-1
            END DO grad_do5
            pShift=pShift+size_p
            xDerivShift=xDerivShift+newSize
            msize_p=msize_p+size_p
        END DO
    END IF
END SUBROUTINE

! *****************************************************************************
!> \brief subroutine that converts from the cp2k poly format to the d3 poly format
! *****************************************************************************
SUBROUTINE poly_cp2k2d3(poly_cp2k,grad,poly_d3,error)
    REAL(dp), DIMENSION(:), INTENT(in)       :: poly_cp2k
    INTEGER, INTENT(in)                      :: grad
    REAL(dp), DIMENSION(:), INTENT(out)      :: poly_d3
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'poly_cp2k2d3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: cp_ii, i, j, k, mgrad, &
                                                mgrad2, sgrad, sgrad2, &
                                                sgrad2k, sgrad3, sgrad3k, &
                                                shifti, shiftj, shiftk, size_p
    LOGICAL                                  :: failure

    failure=.FALSE.
    size_p=(grad+1)*(grad+2)*(grad+3)/6
    CPPrecondition(SIZE(poly_cp2k)>=size_p,cp_failure_level,routineP,error,failure)
    CPPrecondition(SIZE(poly_d3)>=size_p,cp_failure_level,routineP,error,failure)
    cp_ii=0
    sgrad2k=0
    sgrad3k=0
    DO k=0,grad
        shiftk=k+1
        sgrad2k=sgrad2k+k
        sgrad3k=sgrad3k+sgrad2k
        sgrad2=sgrad2k
        sgrad3=sgrad3k
        DO j=0,grad-k
            sgrad=j+k
            mgrad2=sgrad2
            shiftj=mgrad2+shiftk
            mgrad=sgrad
            shifti=shiftj+sgrad3
            DO i=0,grad-j-k
                cp_ii=cp_ii+1
                poly_d3(shifti)=poly_cp2k(cp_ii)
                mgrad=mgrad+1
                mgrad2=mgrad2+mgrad
                shifti=shifti+mgrad2
            END DO
            sgrad2=sgrad2+sgrad+1
            sgrad3=sgrad3+sgrad2
        END DO
    END DO
    IF (SIZE(poly_d3)>=size_p) THEN
        poly_d3(size_p+1:)=0.0_dp
    END IF
END SUBROUTINE

! *****************************************************************************
!> \brief subroutine that converts from the d3 poly format to the cp2k poly format
! *****************************************************************************
SUBROUTINE poly_d32cp2k(poly_cp2k,grad,poly_d3,error)
    REAL(dp), DIMENSION(:), INTENT(out)      :: poly_cp2k
    INTEGER, INTENT(in)                      :: grad
    REAL(dp), DIMENSION(:), INTENT(in)       :: poly_d3
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'poly_d32cp2k', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: cp_ii, i, j, k, mgrad, &
                                                mgrad2, sgrad, sgrad2, &
                                                sgrad2k, sgrad3, sgrad3k, &
                                                shifti, shiftj, shiftk, size_p
    LOGICAL                                  :: failure

    failure=.FALSE.
    size_p=(grad+1)*(grad+2)*(grad+3)/6
    CPPrecondition(SIZE(poly_cp2k)>=size_p,cp_failure_level,routineP,error,failure)
    CPPrecondition(SIZE(poly_d3)>=size_p,cp_failure_level,routineP,error,failure)
    cp_ii=0
    sgrad2k=0
    sgrad3k=0
    DO k=0,grad
        shiftk=k+1
        sgrad2k=sgrad2k+k
        sgrad3k=sgrad3k+sgrad2k
        sgrad2=sgrad2k
        sgrad3=sgrad3k
        DO j=0,grad-k
            sgrad=j+k
            mgrad2=sgrad2
            shiftj=mgrad2+shiftk
            mgrad=sgrad
            shifti=shiftj+sgrad3
            DO i=0,grad-j-k
                cp_ii=cp_ii+1
                poly_cp2k(cp_ii)=poly_d3(shifti)
                mgrad=mgrad+1
                mgrad2=mgrad2+mgrad
                shifti=shifti+mgrad2
            END DO
            sgrad2=sgrad2+sgrad+1
            sgrad3=sgrad3+sgrad2
        END DO
    END DO
    IF (SIZE(poly_d3)>=size_p) THEN
        poly_cp2k(size_p+1:)=0.0_dp
    END IF
END SUBROUTINE

END MODULE d3_poly
