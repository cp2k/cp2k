
 !-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****h* cp2k/xc_atom
!!
!!   NAME
!!     xc_atom
!!     
!!
!!
!!   SOURCE
!****************************************************************************
MODULE xc_atom 

  USE cp_control_types,                 ONLY: xc_control_type
  USE cp_error_handling,               ONLY: cp_error_type
  USE kinds,                           ONLY: dp
  USE xc_functionals,                  ONLY: xc_calculate_lda,&
                                             xc_calculate_lsd,&
                                             xc_get_pot_size,&
                                             xc_get_order,&
                                             xc_set
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type

 IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PRIVATE, PARAMETER :: module_name='xc_atom'

  PUBLIC :: vxc_of_r

!******************************************************************************

CONTAINS

!******************************************************************************

  SUBROUTINE vxc_of_r(xc_control,derivative_set,nspins,na,rho,drho,ir,locpot,error)

    TYPE(xc_control_type)               :: xc_control
    TYPE(xc_derivative_set_type), &
                             POINTER    :: derivative_set
    INTEGER                             :: nspins,na
    REAL(dp), DIMENSION(:,:), &
                             POINTER    :: rho
    INTEGER, INTENT(IN), OPTIONAL       :: ir
    REAL(dp), DIMENSION(:,:,:,:), &
                     OPTIONAL, POINTER  :: drho
    REAL(dp), DIMENSION(:,:), OPTIONAL :: locpot
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                          :: error

    CHARACTER(LEN=*), PARAMETER         :: routine_name = "vxc_of_r"

    REAL(dp), DIMENSION(:), POINTER     :: rhoa,rhob,drhoa,drhob
    CHARACTER(LEN=40), DIMENSION(3)     :: funct
    LOGICAL, DIMENSION(3)               :: grad_funct
    LOGICAL                             :: docalc( 3 ),cross_term,&
                                           gradient_f,lsd
    INTEGER  :: i,ia,ifun,ind(100),istat,lp,npot,order

!******************************************************************************

   gradient_f=PRESENT(drho)
   cross_term = ANY (xc_control%crossterms) .AND. nspins>1

   lsd = (nspins==2)
   order = 1
   npot = xc_get_pot_size(lsd,gradient_f , cross_term, order)

   funct(:) = xc_control%functionals(:)
   grad_funct(:) = xc_control%gradient_functionals(:)

   docalc ( : ) = ( INDEX ( funct(:), "NONE" ) == 0 )
   CALL xc_set(xc_control%density_cut, xc_control%gradient_cut)
   IF(lsd) THEN
      ALLOCATE(rhoa(na),rhob(na),STAT=istat)
      DO ia = 1,na
        rhoa(ia) = rho(ia,1)
        rhob(ia) = rho(ia,2)
      END DO
      IF(PRESENT(drho)) THEN
        ALLOCATE(drhoa(na),drhob(na),STAT=istat)
        DO ia = 1,na
          drhoa(ia) = drho(4,ia,ir,1) 
          drhob(ia) = drho(4,ia,ir,2)
        END DO
      END IF
      DO ifun=1, 3
        IF (docalc(ifun)) THEN
          IF (gradient_f) THEN
            CALL dcopy(na*npot, 0.0_dp, 0, locpot, 1)
            lp = xc_get_pot_size(.FALSE., grad_funct(ifun), .FALSE., order)
            CALL xc_calculate_lsd(derivative_set, funct(ifun),&
                                  rhoa,rhob,locpot, order,&
                                  drhoa=drhoa,drhob=drhob,error=error)
            CALL xc_get_order(.FALSE., grad_funct(ifun), .FALSE., &
                              gradient_f, .FALSE., order, ind(1:lp))
            DO i = 1, lp
               CALL daxpy(na,1._dp, locpot(1,i), 1,& 
                          derivative_set%data(1,ind(i)-1), 1)
            END DO
          ELSE
            CALL xc_calculate_lsd(derivative_set, funct(ifun),&
                                  rhoa,rhob,&
                                  derivative_set%data, order,error=error)
          END IF
        END IF
      END DO
      DEALLOCATE(rhoa,rhob,STAT=istat)
      IF(PRESENT(drho)) THEN
        DEALLOCATE(drhoa,drhob,STAT=istat)
      END IF
   ELSE
      ALLOCATE(rhoa(na),STAT=istat)
      DO ia = 1,na
        rhoa(ia) = rho(ia,1)
      END DO
      IF(PRESENT(drho)) THEN
        ALLOCATE(drhoa(na),STAT=istat) 
        DO ia = 1,na
          drhoa(ia) = drho(4,ia,ir,1)
        END DO
      ENDIF
      DO ifun=1, 3
        IF (docalc(ifun)) THEN
          IF (gradient_f) THEN
            CALL dcopy(na*npot, 0.0_dp, 0, locpot, 1)
            lp = xc_get_pot_size(.FALSE., grad_funct(ifun), .FALSE., order)
            CALL xc_calculate_lda(derivative_set, functional=funct(ifun),&
                                 rho=rhoa, drho=drhoa, &
                                  pot=locpot, order=order, error=error)
            CALL xc_get_order(.FALSE., grad_funct(ifun), .FALSE., &
                              gradient_f, .FALSE., order, ind(1:lp))
            DO i = 1, lp
              CALL daxpy(na,1._dp, locpot(1,i), 1,& 
                         derivative_set%data(1,ind(i)-1), 1)
            END DO
          ELSE

            CALL xc_calculate_lda(derivative_set, functional=funct(ifun),rho=rhoa,&
                                  pot=derivative_set%data, order=order, error=error)

          END IF
        END IF
      END DO
      DEALLOCATE (rhoa,STAT=istat)
      IF(PRESENT(drho)) THEN
        DEALLOCATE(drhoa,STAT=istat)
      END IF
    END IF

  END SUBROUTINE vxc_of_r

!****************************************************************************
END MODULE xc_atom
