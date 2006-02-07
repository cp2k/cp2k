!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2002  CP2K developers group
!-----------------------------------------------------------------------------!
!!****** cp2k/cg_test [1.0] *
!!
!!   NAME
!!     cg_test
!!
!!   FUNCTION
!!     Test of Clebsch-Gordon Coefficients
!!
!!   AUTHOR
!!     JGH (28.02.2002)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*****
!******************************************************************************

MODULE cg_test

! *****************************************************************************

  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE lebedev,                         ONLY: get_number_of_lebedev_grid,&
                                             init_lebedev_grids,&
                                             lebedev_grid
  USE machine,                         ONLY: m_cputime
  USE mathconstants,                   ONLY: pi
  USE spherical_harmonics,             ONLY: clebsch_gordon,&
                                             clebsch_gordon_init,&
                                             y_lm
  USE termination,                     ONLY: stop_memory

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: clebsch_gordon_test

! *****************************************************************************

CONTAINS

!******************************************************************************
SUBROUTINE clebsch_gordon_test ( globenv )

    TYPE(global_environment_type), POINTER   :: globenv

    INTEGER, PARAMETER                       :: l = 7

    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:)                           :: a1, a2, a3
    INTEGER                                  :: ierr, il, iw, l1, l2, ll, lp, &
                                                m1, m2, mm, mp, na
    REAL(KIND=dp)                            :: ca, cga(10), cn, rga(10,21), &
                                                tend, tstart
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: b1, b2, b3, wa

  iw =globenv % scr

  IF ( globenv % ionode ) THEN

     WRITE(iw,'(/,A,/)') " Test of Clebsch-Gordon Coefficients"
     WRITE(iw,'(T40,A,T77,I4)') " Maximum l value tested:",l

     na = 500
     CALL init_lebedev_grids
     ll = get_number_of_lebedev_grid(n=na)
     na = lebedev_grid(ll)%n
     ALLOCATE (wa(na), STAT=ierr)
     IF ( ierr/= 0 ) CALL stop_memory ( "clebsch_gordon_test","wa",na) 
     ALLOCATE (a1(na),a2(na),a3(na), STAT=ierr)
     IF ( ierr/= 0 ) CALL stop_memory ( "clebsch_gordon_test","a",6*na) 
     ALLOCATE (b1(na),b2(na),b3(na), STAT=ierr)
     IF ( ierr/= 0 ) CALL stop_memory ( "clebsch_gordon_test","b",3*na) 

     wa(1:na) = 4.0_dp*pi*lebedev_grid(ll)%w(1:na)

     tstart = m_cputime ( )
     CALL clebsch_gordon_init ( l )
     tend = m_cputime ( )
     tend = tend - tstart
     WRITE(iw,'(T30,A,T71,F10.3)') " Time for Clebsch-Gordon Table [s] ",tend
     lp = (l**4 + 6*l**3 + 15*l**2 + 18*l + 8)/8
     lp = 2*lp*(l+1)
     WRITE(iw,'(T30,A,T71,I10)') "      Size of Clebsch-Gordon Table ",lp
     WRITE(iw,'(/,A)') " Start Test for Complex Spherical Harmonics "

     DO l1=0,l
       DO m1=-l1,l1
         CALL y_lm(lebedev_grid(ll)%r,a1,l1,m1)
         DO l2=0,l
           DO m2=-l2,l2
             CALL y_lm(lebedev_grid(ll)%r,a2,l2,m2)
             CALL clebsch_gordon (l1,m1,l2,m2,cga)
             DO lp=MOD(l1+l2,2),l1+l2,2
               mp = m1 + m2
               IF ( lp < ABS(mp) ) CYCLE
               CALL y_lm(lebedev_grid(ll)%r,a3,lp,mp)
               cn = SUM(a1*a2*CONJG(a3)*wa)
               il = lp/2 + 1
               ca = cga(il)
               IF ( ABS(ca-cn) > 1.e-10_dp ) THEN
                  WRITE ( 6, '(A,3I5,A,F20.12)' ) " l ",l1,l2,lp," A ",ca
                  WRITE ( 6, '(A,3I5,A,F20.12)' ) " m ",m1,m2,mp," N ",cn
                  WRITE ( 6, * )
               END IF
             END DO
           END DO
         END DO
       END DO
       WRITE(iw,'(A,i2,A)') " Test for l = ",l1," done"
     END DO

     WRITE(iw,'(/,A)') " Start Test for Real Spherical Harmonics "
     DO l1=0,l
       DO m1=-l1,l1
         CALL y_lm(lebedev_grid(ll)%r,b1,l1,m1)
         DO l2=0,l
           DO m2=-l2,l2
             CALL y_lm(lebedev_grid(ll)%r,b2,l2,m2)
             CALL clebsch_gordon (l1,m1,l2,m2,rga)
             mp = m1 + m2
             mm = m1 - m2
             IF ( m1*m2 < 0  .OR. (m1*m2==0 .AND. (m1<0 .OR. m2<0))) THEN
               mp = -ABS(mp)
               mm = -ABS(mm)
             ELSE
               mp = ABS(mp)
               mm = ABS(mm)
             END IF
             DO lp=MOD(l1+l2,2),l1+l2,2
               IF ( ABS(mp) <= lp ) THEN
                 CALL y_lm(lebedev_grid(ll)%r,b3,lp,mp)
                 cn = SUM(b1*b2*b3*wa)
                 il = lp/2 + 1
                 ca = rga(il,1)
                 IF ( ABS(ca-cn) > 1.e-10_dp ) THEN
                    WRITE ( 6, '(A,3I5,A,F20.12)' ) " l ",l1,l2,lp," A ",ca
                    WRITE ( 6, '(A,3I5,A,F20.12)' ) " m ",m1,m2,mp," N ",cn
                    WRITE ( 6, * )
                 END IF
               END IF
               IF ( mp /= mm .AND. ABS(mm) <= lp ) THEN
                 CALL y_lm(lebedev_grid(ll)%r,b3,lp,mm)
                 cn = SUM(b1*b2*b3*wa)
                 il = lp/2 + 1
                 ca = rga(il,2)
                 IF ( ABS(ca-cn) > 1.e-10_dp ) THEN
                    WRITE ( 6, '(A,3I5,A,F20.12)' ) " l ",l1,l2,lp," A ",ca
                    WRITE ( 6, '(A,3I5,A,F20.12)' ) " m ",m1,m2,mm," N ",cn
                    WRITE ( 6, * )
                 END IF
               END IF
             END DO
           END DO
         END DO
       END DO
       WRITE(iw,'(A,i2,A)') " Test for l = ",l1," done"
     END DO

     DEALLOCATE ( wa, STAT=ierr)
     IF ( ierr/= 0 ) CALL stop_memory ( "clebsch_gordon_test","wa") 
     DEALLOCATE (a1, a2, a3, STAT=ierr)
     IF ( ierr/= 0 ) CALL stop_memory ( "clebsch_gordon_test","a") 
     DEALLOCATE (b1, b2, b3, STAT=ierr)
     IF ( ierr/= 0 ) CALL stop_memory ( "clebsch_gordon_test","b") 

  END IF

END SUBROUTINE clebsch_gordon_test

! *****************************************************************************

END MODULE cg_test

! *****************************************************************************
