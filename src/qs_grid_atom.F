!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!
! *****************************************************************************
MODULE qs_grid_atom

  
  USE input_constants,                 ONLY: do_gapw_gcs,&
                                             do_gapw_gct,&
                                             do_gapw_log
  USE kinds,                           ONLY: dp
  USE lebedev,                         ONLY: lebedev_grid
  USE mathconstants,                   ONLY: pi
  USE memory_utilities,                ONLY: reallocate
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_grid_atom'

  TYPE grid_atom_type

     INTEGER                         :: nr, ng_sphere
     REAL(dp), DIMENSION(:), POINTER :: rad, rad2,&
          wr, wa,&
          azi, cos_azi, sin_azi, cotan_azi,&  
          pol, cos_pol, sin_pol, usin_pol
     REAL(dp), DIMENSION(:,:), &
          POINTER :: rad2l,oorad2l,weight

  END TYPE grid_atom_type

  PUBLIC :: allocate_grid_atom, &
            create_grid_atom, &
            deallocate_grid_atom

  PUBLIC :: grid_atom_type

CONTAINS

! *****************************************************************************
!> \brief   Deallocate a Gaussian-type orbital (GTO) basis set data set.
!> \param grid_atom ...
!> \param error ...
!> \date    03.11.2000
!> \author  MK
!> \author Matthias Krack (MK)
!> \version 1.0
! *****************************************************************************
  SUBROUTINE allocate_grid_atom(grid_atom,error)

    TYPE(grid_atom_type), POINTER            :: grid_atom
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_grid_atom', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.

    IF (ASSOCIATED(grid_atom)) CALL deallocate_grid_atom(grid_atom,error)

    ALLOCATE (grid_atom,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    NULLIFY (grid_atom%rad)
    NULLIFY (grid_atom%rad2)
    NULLIFY (grid_atom%wr)
    NULLIFY (grid_atom%wa)
    NULLIFY (grid_atom%weight)
    NULLIFY (grid_atom%azi)
    NULLIFY (grid_atom%cos_azi)
    NULLIFY (grid_atom%sin_azi)
    NULLIFY (grid_atom%cotan_azi)
    NULLIFY (grid_atom%pol)
    NULLIFY (grid_atom%cos_pol)
    NULLIFY (grid_atom%sin_pol)
    NULLIFY (grid_atom%usin_pol)
    NULLIFY (grid_atom%rad2l)
    NULLIFY (grid_atom%oorad2l)

  END SUBROUTINE allocate_grid_atom

! *****************************************************************************
!> \brief   Deallocate a Gaussian-type orbital (GTO) basis set data set.
!> \param grid_atom ...
!> \param error ...
!> \date    03.11.2000
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE deallocate_grid_atom(grid_atom,error)
    TYPE(grid_atom_type), POINTER            :: grid_atom
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_grid_atom', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(grid_atom)) THEN

       IF (ASSOCIATED(grid_atom%rad)) THEN
          DEALLOCATE (grid_atom%rad,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF

       IF (ASSOCIATED(grid_atom%rad2)) THEN
          DEALLOCATE (grid_atom%rad2,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF

       IF (ASSOCIATED(grid_atom%wr)) THEN
          DEALLOCATE (grid_atom%wr,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF

       IF (ASSOCIATED(grid_atom%wa)) THEN
          DEALLOCATE (grid_atom%wa,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF

       IF (ASSOCIATED(grid_atom%weight)) THEN
          DEALLOCATE (grid_atom%weight,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF

       IF (ASSOCIATED(grid_atom%azi)) THEN
          DEALLOCATE (grid_atom%azi,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF

       IF (ASSOCIATED(grid_atom%cos_azi)) THEN
          DEALLOCATE (grid_atom%cos_azi,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF

       IF (ASSOCIATED(grid_atom%sin_azi)) THEN
          DEALLOCATE (grid_atom%sin_azi,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF

       IF (ASSOCIATED(grid_atom%cotan_azi)) THEN
          DEALLOCATE (grid_atom%cotan_azi,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF

       IF (ASSOCIATED(grid_atom%pol)) THEN
          DEALLOCATE (grid_atom%pol,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF

       IF (ASSOCIATED(grid_atom%cos_pol)) THEN
          DEALLOCATE (grid_atom%cos_pol,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF

       IF (ASSOCIATED(grid_atom%sin_pol)) THEN
          DEALLOCATE (grid_atom%sin_pol,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF

       IF (ASSOCIATED(grid_atom%usin_pol)) THEN
          DEALLOCATE (grid_atom%usin_pol,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF

       IF (ASSOCIATED(grid_atom%rad2l)) THEN
          DEALLOCATE (grid_atom%rad2l,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF

       IF (ASSOCIATED(grid_atom%oorad2l)) THEN
          DEALLOCATE (grid_atom%oorad2l,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF

       DEALLOCATE (grid_atom,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ELSE
       CALL stop_program(routineN,moduleN,__LINE__,&
            "The pointer grid_atom is not associated and "//&
            "cannot be deallocated")
    END IF
  END SUBROUTINE deallocate_grid_atom

! *****************************************************************************
!> \brief ...
!> \param grid_atom ...
!> \param nr ...
!> \param na ...
!> \param llmax ...
!> \param ll ...
!> \param quadrature ...
! *****************************************************************************
  SUBROUTINE create_grid_atom(grid_atom,nr,na,llmax,ll,quadrature)

    TYPE(grid_atom_type), POINTER            :: grid_atom
    INTEGER, INTENT(IN)                      :: nr, na, llmax, ll, quadrature

    CHARACTER(len=*), PARAMETER :: routineN = 'create_grid_atom', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ia, ir, l
    REAL(dp)                                 :: cosia, pol
    REAL(dp), DIMENSION(:), POINTER          :: rad, rad2, wr

    CALL timeset(routineN,handle)

    NULLIFY (rad,rad2,wr)

    IF (ASSOCIATED(grid_atom)) THEN

       ! Allocate the radial grid arrays
       CALL reallocate(grid_atom%rad,1,nr)
       CALL reallocate(grid_atom%rad2,1,nr)
       CALL reallocate(grid_atom%wr,1,nr)
       CALL reallocate(grid_atom%wa,1,na)
       CALL reallocate(grid_atom%weight,1,na,1,nr)
       CALL reallocate(grid_atom%azi,1,na)
       CALL reallocate(grid_atom%cos_azi,1,na)
       CALL reallocate(grid_atom%sin_azi,1,na)
       CALL reallocate(grid_atom%cotan_azi,1,na)
       CALL reallocate(grid_atom%pol,1,na)
       CALL reallocate(grid_atom%cos_pol,1,na)
       CALL reallocate(grid_atom%sin_pol,1,na)
       CALL reallocate(grid_atom%usin_pol,1,na)
       CALL reallocate(grid_atom%rad2l,1,nr,0,llmax+1)
       CALL reallocate(grid_atom%oorad2l,1,nr,0,llmax+1)

       ! Calculate the radial grid for this kind
       rad => grid_atom%rad
       rad2 => grid_atom%rad2
       wr => grid_atom%wr

       CALL radial_grid(nr,rad,rad2,wr,quadrature)

       grid_atom%rad2l(:,0) = 1._dp
       grid_atom%oorad2l(:,0) = 1._dp
       DO l = 1,llmax+1
          grid_atom%rad2l(:,l)=grid_atom%rad2l(:,l-1)*rad(:)
          grid_atom%oorad2l(:,l)=grid_atom%oorad2l(:,l-1)/rad(:)
       ENDDO

       IF (ll > 0) THEN
          grid_atom%wa(1:na) = 4._dp*pi*lebedev_grid(ll)%w(1:na)
          DO ir = 1,nr
             DO ia= 1,na
                grid_atom%weight(ia,ir) = grid_atom%wr(ir)*grid_atom%wa(ia)
             END DO
          END DO

          DO ia = 1,na
             cosia = lebedev_grid(ll)%r(3,ia)
             grid_atom%cos_pol(ia) = cosia
             IF ( ABS(lebedev_grid(ll)%r(2,ia)) < EPSILON(1.0_dp) .AND. &
                  ABS(lebedev_grid(ll)%r(1,ia)) < EPSILON(1.0_dp) ) THEN
                grid_atom%azi(ia) = 0.0_dp
             ELSE
                grid_atom%azi(ia) = ATAN2( lebedev_grid(ll)%r(2,ia) , lebedev_grid(ll)%r(1,ia) )
             END IF
             grid_atom%cos_azi(ia) = COS( grid_atom%azi(ia) )
             pol =  ACOS(cosia)
             IF( grid_atom%sin_pol(ia) < 0.0_dp ) pol = - pol
             grid_atom%pol(ia) = pol
             grid_atom%sin_pol(ia)  = SIN(grid_atom%pol(ia))
             IF(ABS(grid_atom%sin_azi(ia)) > EPSILON(1.0_dp)) THEN
                grid_atom%cotan_azi(ia) = &
                     grid_atom%cos_azi(ia)/grid_atom%sin_azi(ia)
             ELSE
                grid_atom%cotan_azi(ia) = 0.0_dp
             END IF

             grid_atom%sin_pol(ia) = SIN( grid_atom%pol(ia) )
             IF(ABS(grid_atom%sin_pol(ia)) > EPSILON(1.0_dp)) THEN
                grid_atom%usin_pol(ia)  = 1.0_dp/grid_atom%sin_pol(ia)
             ELSE
                grid_atom%usin_pol(ia) = 1.0_dp
             END IF

          ENDDO

       END IF

    ELSE
       CALL stop_program(routineN,moduleN,__LINE__,&
            "The pointer grid_atom is not associated")
    END IF

    CALL timestop(handle)

  END SUBROUTINE create_grid_atom

! *****************************************************************************
!> \brief   Generate a radial grid with n points by a quadrature rule.
!> \param n ...
!> \param r ...
!> \param r2 ...
!> \param wr ...
!> \param radial_quadrature ...
!> \date    20.09.1999
!> \par Literature
!>           - A. D. Becke, J. Chem. Phys. 88, 2547 (1988)
!>           - J. M. Perez-Jorda, A. D. Becke and E. San-Fabian,
!>             J. Chem. Phys. 100, 6520 (1994)
!>           - M. Krack and A. M. Koester, J. Chem. Phys. 108, 3226 (1998)
!> \author  Matthias Krack
!> \version 1.0
! *****************************************************************************
  SUBROUTINE radial_grid(n,r,r2,wr,radial_quadrature)

    INTEGER, INTENT(IN)                      :: n
    REAL(dp), DIMENSION(:), POINTER          :: r, r2, wr
    INTEGER, INTENT(IN)                      :: radial_quadrature

    CHARACTER(len=*), PARAMETER :: routineN = 'radial_grid', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    REAL(dp)                                 :: cost, f, sint, sint2, t, w, x

    f = pi/REAL(n + 1,dp)

    SELECT CASE (radial_quadrature)

    CASE (do_gapw_gcs)

       !     *** Gauss-Chebyshev quadrature formula of the second kind ***
       !     *** u [-1,+1] -> r [0,infinity] =>  r = (1 + u)/(1 - u)   ***

       DO i=1,n
          t = REAL(i,dp)*f
          x = COS(t)
          w = f*SIN(t)**2
          r(i) = (1.0_dp + x)/(1.0_dp - x)
          r2(i) = r(i)**2
          wr(i) = w/SQRT(1.0_dp - x**2)
          wr(i) = 2.0_dp*wr(i)*r2(i)/(1.0_dp - x)**2
       END DO

    CASE (do_gapw_gct)

       !     *** Transformed Gauss-Chebyshev quadrature formula of the second kind ***
       !     *** u [-1,+1] -> r [0,infinity] => r = (1 + u)/(1 - u)                ***

       DO i=1,n
          t = REAL(i,dp)*f
          cost = COS(t)
          sint = SIN(t)
          sint2 = sint**2
          x = REAL(2*i - n - 1,dp)/REAL(n + 1,dp) -&
               2.0_dp*(1.0_dp + 2.0_dp*sint2/3.0_dp)*cost*sint/pi
          w = 16.0_dp*sint2**2/REAL(3*(n + 1),dp)
          r(n+1-i) = (1.0_dp + x)/(1.0_dp - x)
          r2(n+1-i) = r(n+1-i)**2
          wr(n+1-i) = 2.0_dp*w*r2(n+1-i)/(1.0_dp - x)**2
       END DO

    CASE (do_gapw_log)

       !     *** Transformed Gauss-Chebyshev quadrature formula of the second kind ***
       !     *** u [-1,+1] -> r [0,infinity] => r = ln(2/(1 - u))/ln(2)            ***

       DO i=1,n
          t = REAL(i,dp)*f
          cost = COS(t)
          sint = SIN(t)
          sint2 = sint**2
          x = REAL(2*i - n - 1,dp)/REAL(n + 1,dp) -&
               2.0_dp*(1.0_dp + 2.0_dp*sint2/3.0_dp)*cost*sint/pi
          w = 16.0_dp*sint2**2/REAL(3*(n + 1),dp)
          r(n+1-i) = LOG(2.0_dp/(1.0_dp - x))/LOG(2.0_dp)
          r2(n+1-i) = r(n+1-i)**2
          wr(n+1-i) = w*r2(n+1-i)/(LOG(2.0_dp)*(1.0_dp - x))
       END DO

    CASE DEFAULT

       CALL stop_program(routineN,moduleN,__LINE__,&
                         "Invalid radial quadrature type specified")

    END SELECT

  END SUBROUTINE radial_grid

END MODULE qs_grid_atom
