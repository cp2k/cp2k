!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/ai_spherical_gaussian [1.0] *
!!
!!   NAME
!!     ai_spherical_gaussian
!!
!!   FUNCTION
!!
!!   Calculates atomic integrals over normalized spherical Gaussian functions
!!
!!   X(lm,p) = N(lp) r^l exp(-p*r^2) Y_lm(theta,phi)
!!
!!   N(lp) = 2^(l+2) [(2l+1)!!]^(-1/2) (2Pi)^(-1/4) p^(l/2+3/4)
!!
!!   AUTHOR
!!     JGH : 16-06-2002
!!
!!   MODIFICATION HISTORY
!!     NONE
!!
!!*****
!******************************************************************************

MODULE ai_spherical_gaussian

  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: dfac,&
                                             fac,&
                                             gamma1,&
                                             pi,&
                                             twopi
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: sg_norm, sg_overlap, sg_kinetic, sg_nuclear, sg_erf, sg_gpot, &
            sg_proj_ol

!******************************************************************************

CONTAINS

!******************************************************************************

!------------------------------------------------------------------------------
!
!  S(l,pq) = [2 SQRT(pq)/(p+q)]^(l+1.5)
!
!------------------------------------------------------------------------------
SUBROUTINE sg_norm ( pnorm, l, p )


    REAL(KIND=dp), DIMENSION(:), INTENT(OUT)     :: pnorm
    INTEGER, INTENT(IN)                      :: l
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)      :: p

    REAL(KIND=dp)                                :: pe, pf

!------------------------------------------------------------------------------

  IF ( SIZE ( pnorm ) < SIZE ( p ) ) CALL stop_program ( "sg_norm", &
     "Incompatible array sizes" )

  pf = 2**(l+2) * twopi**(-0.25_dp) / SQRT ( dfac(2*l+1) )
  pe = 0.25_dp*(2*l+3)
  pnorm (:) = pf * p(:)**pe

END SUBROUTINE sg_norm

!------------------------------------------------------------------------------
! 
!  S(l,pq) = [2 SQRT(pq)/(p+q)]^(l+1.5) 
! 
!------------------------------------------------------------------------------
SUBROUTINE sg_overlap ( smat, l, pa, pb )


    REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)  :: smat
    INTEGER, INTENT(IN)                      :: l
    REAL(KIND=dp), DIMENSION(:), INTENT(IN), &
      TARGET                                 :: pa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN), &
      OPTIONAL, TARGET                       :: pb

    INTEGER                                  :: ip, iq, m, n, up
    LOGICAL                                  :: symmetric
    REAL(KIND=dp), DIMENSION(:), POINTER         :: p, q

!------------------------------------------------------------------------------

  p => pa

  IF ( PRESENT ( pb ) ) THEN
    q => pb
    symmetric = .FALSE.
  ELSE
    q => pa
    symmetric = .TRUE.
  END IF
  
  n = SIZE ( p ) 
  m = SIZE ( q ) 

  IF ( n > SIZE ( smat, 1 ) .OR. m > SIZE ( smat, 2 ) ) &
     CALL stop_program ( "sg_overlap", "matrix dimension too small" )

  DO ip = 1, n
    up = 1
    IF ( symmetric ) up = ip
    DO iq = up, m
      smat ( ip, iq ) = ( 2.0_dp * SQRT ( p(ip) * q(iq) ) / (p(ip) + q(iq)) ) &
                        ** ( l + 1.5_dp )
    END DO
  END DO
  IF ( symmetric ) THEN
    DO ip = 1, n
      DO iq = ip, m
        smat ( iq, ip ) = smat ( ip, iq )
      END DO
    END DO
  END IF

END SUBROUTINE sg_overlap

!------------------------------------------------------------------------------
! 
!  T(l,pq) = (2l+3) [pq/(p+q)] S(l,pq)
! 
!------------------------------------------------------------------------------
SUBROUTINE sg_kinetic ( kmat, l, pa, pb, smat )


    REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)  :: kmat
    INTEGER, INTENT(IN)                      :: l
    REAL(KIND=dp), DIMENSION(:), INTENT(IN), &
      TARGET                                 :: pa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN), &
      OPTIONAL, TARGET                       :: pb
    REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
      OPTIONAL                               :: smat

    INTEGER                                  :: ip, iq, m, n, up
    LOGICAL                                  :: symmetric
    REAL(KIND=dp)                                :: tld
    REAL(KIND=dp), DIMENSION(:), POINTER         :: p, q

!------------------------------------------------------------------------------

  p => pa

  IF ( PRESENT ( pb ) ) THEN
    q => pb
    symmetric = .FALSE.
  ELSE
    q => pa
    symmetric = .TRUE.
  END IF
  
  n = SIZE ( p ) 
  m = SIZE ( q ) 

  IF ( n > SIZE ( kmat, 1 ) .OR. m > SIZE ( kmat, 2 ) ) &
     CALL stop_program ( "sg_kinetic", "matrix dimension too small" )

  IF ( PRESENT ( smat ) ) THEN
     kmat = smat
  ELSE
     CALL sg_overlap ( kmat, l, p, q )
  END IF

  tld = 2.0_dp * REAL(l,dp) + 3.0_dp
  DO ip = 1, n
    up = 1
    IF ( symmetric ) up = ip
    DO iq = up, m
      kmat ( ip, iq ) = tld * ( p(ip) * q(iq) ) / (p(ip) + q(iq) ) * kmat ( ip, iq )
    END DO
  END DO
  IF ( symmetric ) THEN
    DO ip = 1, n
      DO iq = ip, m
        kmat ( iq, ip ) = kmat ( ip, iq )
      END DO
    END DO
  END IF

END SUBROUTINE sg_kinetic

!------------------------------------------------------------------------------
! 
!  U(l,pq) = l!/(2l+1)!! SQRT[(p+q)/2Pi] 2^(l+1.5) S(l,pq)
! 
!------------------------------------------------------------------------------
SUBROUTINE sg_nuclear ( umat, l, pa, pb, smat )


    REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)  :: umat
    INTEGER, INTENT(IN)                      :: l
    REAL(KIND=dp), DIMENSION(:), INTENT(IN), &
      TARGET                                 :: pa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN), &
      OPTIONAL, TARGET                       :: pb
    REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
      OPTIONAL                               :: smat

    INTEGER                                  :: ip, iq, m, n, up
    LOGICAL                                  :: symmetric
    REAL(KIND=dp)                                :: tld
    REAL(KIND=dp), DIMENSION(:), POINTER         :: p, q

!------------------------------------------------------------------------------

  p => pa

  IF ( PRESENT ( pb ) ) THEN
    q => pb
    symmetric = .FALSE.
  ELSE
    q => pa
    symmetric = .TRUE.
  END IF
  
  n = SIZE ( p ) 
  m = SIZE ( q ) 

  IF ( n > SIZE ( umat, 1 ) .OR. m > SIZE ( umat, 2 ) ) &
     CALL stop_program ( "sg_nuclear", "matrix dimension too small" )

  IF ( PRESENT ( smat ) ) THEN
     umat = smat
  ELSE
     CALL sg_overlap ( umat, l, p, q )
  END IF

  tld = 2.0_dp ** ( l + 1.5_dp ) * fac(l) / dfac (2*l+1) / SQRT ( twopi )
  DO ip = 1, n
    up = 1
    IF ( symmetric ) up = ip
    DO iq = up, m
      umat ( ip, iq ) = tld * SQRT ( p(ip) + q(iq) ) * umat ( ip, iq )
    END DO
  END DO
  IF ( symmetric ) THEN
    DO ip = 1, n
      DO iq = ip, m
        umat ( iq, ip ) = umat ( ip, iq )
      END DO
    END DO
  END IF

END SUBROUTINE sg_nuclear

!------------------------------------------------------------------------------
! 
!  UP(l,pq,a) = N(lp) N(lq) * [ Gamma(l+1)/(2*(p+q)^(l+1)) - 
!       Gamma(l+1.5)/SQRT(Pi)/(2l+2)/a^(2l+2) 2_F_1[l+1,l+1.5;l+2;-(p+q)/a^2] ]
! 
!  UP(l,pq,a) = SQRT[(p+q)/Pi] * 2a/(2l+1)!! * (a^2+p+q)^(-l/2-1) * F(a,p,q,l)
!
!  F(a,p,q,0) = 1 
!  F(a,p,q,1) = 2*a^2 + 3(p+q) 
!  F(a,p,q,2) = 8*a^4 + 20*a^2 (p+q) + 15(p+q)^2 
!  F(a,p,q,3) = 48*a^6 + 168*a^4 (p+q) + 210*a^2 (p+q)^2 + 105(p+q)^3
!  F(a,p,q,4) = 384*a^8 + 1728*a^6 (p+q) + 3024*a^4 (p+q)^2 + 2520*a^2 (p+q)^3
!               + 945 (p+q)^4
! 
!------------------------------------------------------------------------------
SUBROUTINE sg_erf ( upmat, l, a, pa, pb, smat )


    REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)  :: upmat
    INTEGER, INTENT(IN)                      :: l
    REAL(KIND=dp), INTENT(IN)                    :: a
    REAL(KIND=dp), DIMENSION(:), INTENT(IN), &
      TARGET                                 :: pa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN), &
      OPTIONAL, TARGET                       :: pb
    REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
      OPTIONAL                               :: smat

    INTEGER                                  :: ip, iq, m, n, up
    LOGICAL                                  :: symmetric
    REAL(KIND=dp)                                :: a2, fpol, pq, tld
    REAL(KIND=dp), DIMENSION(:), POINTER         :: p, q

!------------------------------------------------------------------------------

  p => pa

  IF ( PRESENT ( pb ) ) THEN
    q => pb
    symmetric = .FALSE.
  ELSE
    q => pa
    symmetric = .TRUE.
  END IF
  
  n = SIZE ( p ) 
  m = SIZE ( q ) 

  IF ( n > SIZE ( upmat, 1 ) .OR. m > SIZE ( upmat, 2 ) ) &
     CALL stop_program ( "sg_erf", "matrix dimension too small" )

  IF ( PRESENT ( smat ) ) THEN
     upmat = smat
  ELSE
     CALL sg_overlap ( upmat, l, p, q )
  END IF

  a2 = a*a
  tld = 2.0_dp*a / SQRT(pi) / dfac(2*l+1)
  DO ip = 1, n
    up = 1
    IF ( symmetric ) up = ip
    DO iq = up, m
      upmat ( ip, iq ) = tld * SQRT ( p(ip) + q(iq) ) * upmat ( ip, iq ) / &
                         ( a2 + p(ip) + q(iq) ) ** (0.5_dp*(l+1))
    END DO
  END DO
  IF ( symmetric ) THEN
    DO ip = 1, n
      DO iq = ip, m
        upmat ( iq, ip ) = upmat ( ip, iq )
      END DO
    END DO
  END IF

  DO ip = 1, n
    up = 1
    IF ( symmetric ) up = ip
    SELECT CASE ( l )
    CASE DEFAULT
      CALL stop_program ( "sg_erf", "l quantum number not implemented" )
    CASE ( 0 )
      ! nothing left to do
    CASE ( 1 )
      DO iq = up, m
        pq = p(ip) + q(iq)
        fpol = 2.0_dp*a2 + 3.0_dp*pq
        upmat ( ip, iq ) = upmat ( ip, iq ) * fpol
      END DO
    CASE ( 2 )
      DO iq = up, m
        pq = p(ip) + q(iq)
        fpol = 8.0_dp*a2*a2 + 20.0_dp*a2*pq + 15.0_dp*pq*pq
        upmat ( ip, iq ) = upmat ( ip, iq ) * fpol
      END DO
    CASE ( 3 )
      DO iq = up, m
        pq = p(ip) + q(iq)
        fpol = 48.0_dp*a2*a2*a2 + 168.0_dp*a2*a2*pq + 210.0_dp*a2*pq*pq &
               + 105.0_dp*pq*pq*pq
        upmat ( ip, iq ) = upmat ( ip, iq ) * fpol
      END DO
    CASE ( 4 )
      DO iq = up, m
        pq = p(ip) + q(iq)
        fpol = 384.0_dp*a2*a2*a2*a2 + 1728.0_dp*a2*a2*a2*pq + &
               3024.0_dp*a2*a2*pq*pq + 2520.0_dp*a2*pq*pq*pq + 945.0_dp*pq*pq*pq*pq
        upmat ( ip, iq ) = upmat ( ip, iq ) * fpol
      END DO
    END SELECT
  END DO
  IF ( symmetric ) THEN
    DO ip = 1, n
      DO iq = ip, m
        upmat ( iq, ip ) = upmat ( ip, iq )
      END DO
    END DO
  END IF

END SUBROUTINE sg_erf

!------------------------------------------------------------------------------
! 
!  Overlap with Projectors P(l,k,rc) for k=0,1,..
!
!  P(l,k,rc) = SQRT(2)/SQRT(Gamma[l+1.5])/rc^(l+2k+1.5) r^(l+2k) exp[-0.5(r/rc)^2]
!
!  SP(l,k,p,rc) = N(l,p) / SQRT(2) / SQRT(gamma[l+2k+1.5]) / rc^(l+2k+1.5) 
!                    * Gamma(l+k+1.5) / (p+1/2rc^2)^(l+k+1.5)
! 
!------------------------------------------------------------------------------
SUBROUTINE sg_proj_ol ( spmat, l, p, norm, k, rc )


    REAL(KIND=dp), DIMENSION(:), INTENT(OUT)     :: spmat
    INTEGER, INTENT(IN)                      :: l
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)      :: p, norm
    INTEGER, INTENT(IN)                      :: k
    REAL(KIND=dp), INTENT(IN)                    :: rc

    INTEGER                                  :: pn
    REAL(KIND=dp)                                :: orc, pe, pf

!------------------------------------------------------------------------------

  IF ( SIZE ( spmat ) < SIZE ( p ) ) CALL stop_program ( "sg_proj_ol", &
     "Array sizes incompatible" )

  pe = l + k + 1.5_dp
  pn = l + k + 1
  pf = gamma1 ( pn ) / SQRT(2.0_dp) / rc**pe / SQRT(gamma1(pn+k))
  orc = 1.0_dp / ( 2.0_dp * rc*rc )

  spmat(:) = norm(:) * pf / ( p(:) + orc )**pe

END SUBROUTINE sg_proj_ol

!------------------------------------------------------------------------------
!
!  Matrix elements for Gaussian potentials
!
!  V(k,rc) = (r/rc)^2n exp[-1/2(r/rc)^2]
!
!  VP(l,k,p,rc) = N(l,p) N(l,q) / SQRT(2) / SQRT(gamma[l+2k+1.5]) / rc^(l+2k+1.5) 
!                    * Gamma(l+k+1.5) / (p+q+1/2rc^2)^(l+k+1.5)
!
!------------------------------------------------------------------------------
SUBROUTINE sg_gpot ( vpmat, k, rc, l, pa, pb )


    REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)  :: vpmat
    INTEGER, INTENT(IN)                      :: k
    REAL(KIND=dp), INTENT(IN)                    :: rc
    INTEGER, INTENT(IN)                      :: l
    REAL(KIND=dp), DIMENSION(:), INTENT(IN), &
      TARGET                                 :: pa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN), &
      OPTIONAL, TARGET                       :: pb

    INTEGER                                  :: ierr, ip, iq, m, n, up
    LOGICAL                                  :: symmetric
    REAL(KIND=dp)                                :: orc, tld
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)     :: np, nq
    REAL(KIND=dp), DIMENSION(:), POINTER         :: p, q

!------------------------------------------------------------------------------

  p => pa

  IF ( PRESENT ( pb ) ) THEN
    q => pb
    symmetric = .FALSE.
  ELSE
    q => pa
    symmetric = .TRUE.
  END IF

  n = SIZE ( p )
  m = SIZE ( q )

  IF ( n > SIZE ( vpmat, 1 ) .OR. m > SIZE ( vpmat, 2 ) ) &
     CALL stop_program ( "sg_gpot", "matrix dimension too small" )

  ALLOCATE ( np ( n ), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "sg_gpot", "np", n )
  ALLOCATE ( nq ( m ), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "sg_gpot", "nq", m )

  CALL sg_norm ( np, l, p )
  CALL sg_norm ( nq, l, q )

  tld = gamma1 ( l + k + 1 ) * 0.5_dp
  orc = 1.0_dp / ( 2.0_dp * rc*rc )

  DO ip = 1, n
    up = 1
    IF ( symmetric ) up = ip
    DO iq = up, m
      vpmat ( ip, iq ) = tld * np ( ip ) * nq ( iq ) / &
                         ( orc + p(ip) + q(iq) ) ** (l+k+1.5_dp)
    END DO
  END DO
  IF ( symmetric ) THEN
    DO ip = 1, n
      DO iq = ip, m
        vpmat ( iq, ip ) = vpmat ( ip, iq )
      END DO
    END DO
  END IF

  DEALLOCATE ( np, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "sg_gpot", "np" )
  DEALLOCATE ( nq, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "sg_gpot", "nq" )

END SUBROUTINE sg_gpot

!******************************************************************************

END MODULE ai_spherical_gaussian

!******************************************************************************
