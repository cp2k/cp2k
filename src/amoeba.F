!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/amoeba [1.0] *
!!
!!   NAME
!!     amoeba, NUMERICAL RECIPIES SUBROUTINE AMOEBA ADAPTED FOR USE IN CP2K
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE amoeba
  
  USE kinds, ONLY : dbl
  USE nrutil, ONLY : assert_eq, imaxloc, iminloc, nrerror, swap
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: amoeba_evaluate
  
!!*****

CONTAINS

!!****** amoeba/amoeba_evaluate [1.0] *
!!
!!   NAME
!!     amoeba_evaluate
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
!-----------------------------------------------------------------------------!
!   NUMERICAL RECIPIES SUBROUTINE AMOEBA ADAPTED FOR USE IN CP2K              !
!-----------------------------------------------------------------------------!

SUBROUTINE amoeba_evaluate ( p, y, ftol, rtol, func, itmax )
  
  IMPLICIT NONE
  
! Arguments
  REAL (dbl), DIMENSION (:,:), INTENT (INOUT) :: p
  REAL (dbl), DIMENSION (:), INTENT (INOUT) :: y
  REAL (dbl), INTENT (IN) :: ftol
  REAL (dbl), INTENT ( OUT ) :: rtol
  INTEGER, INTENT ( IN ) :: itmax
  
! Interface
  INTERFACE
     FUNCTION func ( x )
       USE kinds, ONLY : dbl
       IMPLICIT NONE
       REAL ( dbl ), DIMENSION ( : ), INTENT ( IN ) :: x
       REAL ( dbl ) :: func
     END FUNCTION func
  END INTERFACE
  
  INTEGER :: ihi, ndim
  REAL (dbl), DIMENSION (size(p,2)) :: psum
  
  CALL amoeba_private()

!------------------------------------------------------------------------------

CONTAINS
  
  SUBROUTINE amoeba_private
    IMPLICIT NONE
    INTEGER :: i, ilo, inhi
    INTEGER :: iter
    REAL (dbl) :: ysave, ytry, ytmp
    
    ndim = assert_eq(size(p,2),size(p,1)-1,size(y)-1,'amoeba')
    iter = 0
    psum(:) = sum(p(:,:),dim=1)
    DO
       ilo = iminloc(y(:))
       ihi = imaxloc(y(:))
       ytmp = y(ihi)
       y(ihi) = y(ilo)
       inhi = imaxloc(y(:))
       y(ihi) = ytmp
       rtol = 2.0_dbl*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo)))
       IF (rtol<ftol) THEN
          CALL swap(y(1),y(ilo))
          CALL swap(p(1,:),p(ilo,:))
          RETURN
       END IF
       IF (iter>=itmax) RETURN
       ytry = amotry(-1.0_dbl)
       iter = iter + 1
       IF (ytry<=y(ilo)) THEN
          ytry = amotry(2.0_dbl)
          iter = iter + 1
       ELSE IF (ytry>=y(inhi)) THEN
          ysave = y(ihi)
          ytry = amotry(0.5_dbl)
          iter = iter + 1
          IF (ytry>=ysave) THEN
             p(:,:) = 0.5_dbl*(p(:,:)+spread(p(ilo,:),1,size(p,1)))
             DO i = 1, ndim + 1
                IF (i/=ilo) y(i) = func(p(i,:))
             END DO
             iter = iter + ndim
             psum(:) = sum(p(:,:),dim=1)
          END IF
       END IF
    END DO
  END SUBROUTINE amoeba_private
  
!BL
  FUNCTION amotry(fac)
    IMPLICIT NONE
    REAL (dbl), INTENT (IN) :: fac
    REAL (dbl) :: amotry
    REAL (dbl) :: fac1, fac2, ytry
    REAL (dbl), DIMENSION (size(p,2)) :: ptry
    
    fac1 = (1.0_dbl-fac)/ndim
    fac2 = fac1 - fac
    ptry(:) = psum(:)*fac1 - p(ihi,:)*fac2
    ytry = func(ptry)
    IF (ytry<y(ihi)) THEN
       y(ihi) = ytry
       psum(:) = psum(:) - p(ihi,:) + ptry(:)
       p(ihi,:) = ptry(:)
    END IF
    amotry = ytry
  END FUNCTION amotry
  
END SUBROUTINE amoeba_evaluate

!!*****
!******************************************************************************

END MODULE amoeba
