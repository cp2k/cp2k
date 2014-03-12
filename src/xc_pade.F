!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculate the LDA functional in the Pade approximation
!>      Literature: S. Goedecker, M. Teter and J. Hutter,
!>                  Phys. Rev. B 54, 1703 (1996)
!> \note
!>      Order of derivatives is: LDA 0; 1; 2; 3;
!>                               LSD 0; a  b; aa ab bb; aaa aab abb bbb;
!> \par History
!>      JGH (26.02.2003) : OpenMP enabled
!> \author JGH (15.02.2002)
! *****************************************************************************
MODULE xc_pade
  USE bibliography,                    ONLY: Goedecker1996,&
                                             cite_reference
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE termination,                     ONLY: stop_program
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                             xc_dset_get_derivative
  USE xc_derivative_types,             ONLY: xc_derivative_get,&
                                             xc_derivative_type
  USE xc_functionals_utilities,        ONLY: calc_fx,&
                                             calc_rs,&
                                             calc_rs_pw,&
                                             set_util,&
                                             setup_calculation
  USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
  USE xc_rho_set_types,                ONLY: xc_rho_set_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  REAL(KIND=dp), PARAMETER :: f13 = 1.0_dp/3.0_dp, &
                              f23 = 2.0_dp*f13, &
                              f43 = 4.0_dp*f13

  REAL(KIND=dp), PARAMETER :: a0 = 0.4581652932831429E+0_dp,&
                              a1 = 0.2217058676663745E+1_dp,&
                              a2 = 0.7405551735357053E+0_dp,&
                              a3 = 0.1968227878617998E-1_dp,&
                              b1 = 1.0000000000000000E+0_dp,&
                              b2 = 0.4504130959426697E+1_dp,&
                              b3 = 0.1110667363742916E+1_dp,&
                              b4 = 0.2359291751427506E-1_dp

  REAL(KIND=dp), PARAMETER :: da0 = 0.119086804055547E+0_dp,&
                              da1 = 0.6157402568883345E+0_dp,&
                              da2 = 0.1574201515892867E+0_dp,&
                              da3 = 0.3532336663397157E-2_dp,&
                              db1 = 0.0000000000000000E+0_dp,&
                              db2 = 0.2673612973836267E+0_dp,&
                              db3 = 0.2052004607777787E+0_dp,&
                              db4 = 0.4200005045691381E-2_dp

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_pade'

  PUBLIC :: pade, pade_lda_pw_eval, pade_lsd_pw_eval

  REAL(KIND=dp) :: eps_rho
  LOGICAL :: debug_flag

  INTERFACE pade
    MODULE PROCEDURE pade_lda, pade_lsd, pade_init, pade_info
  END INTERFACE

CONTAINS

! *****************************************************************************
  SUBROUTINE pade_init ( cutoff, debug )

    REAL(KIND=dp), INTENT(IN)                :: cutoff
    LOGICAL, INTENT(IN), OPTIONAL            :: debug

    eps_rho = cutoff
    CALL set_util ( cutoff )

    CALL cite_reference(Goedecker1996)

    IF ( PRESENT ( debug ) ) THEN
       debug_flag = debug
    ELSE
       debug_flag = .FALSE.
    END IF

  END SUBROUTINE pade_init

! *****************************************************************************
  SUBROUTINE pade_info ( reference, shortform, lsd, needs, max_deriv, error )

    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    LOGICAL, INTENT(IN), OPTIONAL            :: lsd
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pade_info', &
      routineP = moduleN//':'//routineN

    IF ( PRESENT ( reference ) ) THEN
       reference = "S. Goedecker, M. Teter and J. Hutter," &
                   //" Phys. Rev. B 54, 1703 (1996)"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
       shortform = "S. Goedecker et al., PRB 54, 1703 (1996)"
    END IF

    IF (PRESENT(needs)) THEN
       IF (.NOT. PRESENT(lsd)) &
            CALL stop_program(routineN,moduleN,__LINE__,"Arguments mismatch.")
       IF (lsd) THEN
          needs%rho_spin = .TRUE.
       ELSE
          needs%rho = .TRUE.
       END IF
    END IF

    IF (PRESENT(max_deriv)) max_deriv=3

  END SUBROUTINE pade_info

! *****************************************************************************
  SUBROUTINE pade_lda(derivative_set, rho, pot, order, error)

    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rho
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: pot
    INTEGER, INTENT(IN), OPTIONAL            :: order
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pade_lda', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: m, m1, m2, mp, mt, mx(0:3,2), &
                                                n, stat
    LOGICAL                                  :: calc(0:4), failure
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: rs

    failure = .FALSE.
    CALL setup_calculation(order,mx,calc,100)
    mp = SIZE ( pot, 2 )
    m = ABS(order)
    mt = mx(m,2)
    IF ( mt > mp ) CALL stop_program(routineN,moduleN,__LINE__,&
                                     "Size of POT array too small.")

    n = SIZE ( rho )
    ALLOCATE ( rs(n), STAT=stat )
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL calc_rs ( rho, rs )
    IF ( calc(0) .AND. calc(1) ) THEN
       m1=mx(0,1)
       m2=mx(1,2)
       CALL pade_lda_01(n, rho, rs, pot(:,m1:m1), pot(:,m2:m2))
    ELSE IF ( calc(0) ) THEN
       CALL pade_lda_0(n, rho, rs, pot(:,1:1))
    ELSE IF ( calc(1) ) THEN
       m1=mx(1,1)
       CALL pade_lda_1(n, rho, rs, pot(:,m1:m1))
    END IF
    IF ( calc(2) ) THEN
       m1=mx(2,1)
       CALL pade_lda_2(n, rho, rs, pot(:,m1:m1))
    END IF
    IF ( calc(3) ) THEN
       m1=mx(3,1)
       CALL pade_lda_3(n, rho, rs, pot(:,m1:m1))
    END IF

    DEALLOCATE ( rs, STAT=stat )
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE pade_lda

! *****************************************************************************
  SUBROUTINE pade_lda_pw_eval(deriv_set, rho_set, order, error)

    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    INTEGER, INTENT(IN), OPTIONAL            :: order
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pade_lda_pw_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: n, stat
    LOGICAL                                  :: calc(0:4), failure
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: rs
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: e_0, e_r, e_rr, e_rrr
    TYPE(xc_derivative_type), POINTER        :: deriv

    failure = .FALSE.
    calc = .FALSE.
    IF (order >= 0) calc(0:order) = .TRUE.
    IF (order < 0) calc(-order) = .TRUE.

    n = PRODUCT(rho_set%local_bounds(2,:)-rho_set%local_bounds(1,:)+(/1,1,1/))
    ALLOCATE ( rs(n), STAT=stat )
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL calc_rs_pw ( rho_set%rho, rs, n )
    IF ( calc(0) .AND. calc(1) ) THEN
       deriv => xc_dset_get_derivative(deriv_set, "", &
                                     allocate_deriv=.TRUE., error=error)
       CALL xc_derivative_get(deriv,deriv_data=e_0,error=error)
       deriv => xc_dset_get_derivative(deriv_set, "(rho)", &
                                     allocate_deriv=.TRUE., error=error)
       CALL xc_derivative_get(deriv,deriv_data=e_r,error=error)
       CALL pade_lda_01(n, rho_set%rho, rs, e_0, e_r)
    ELSE IF ( calc(0) ) THEN
       deriv => xc_dset_get_derivative(deriv_set, "", &
                                     allocate_deriv=.TRUE., error=error)
       CALL xc_derivative_get(deriv,deriv_data=e_0,error=error)
       CALL pade_lda_0(n, rho_set%rho, rs, e_0)
    ELSE IF ( calc(1) ) THEN
       deriv => xc_dset_get_derivative(deriv_set, "(rho)", &
                                     allocate_deriv=.TRUE., error=error)
       CALL xc_derivative_get(deriv,deriv_data=e_r,error=error)
       CALL pade_lda_1(n, rho_set%rho, rs, e_r)
    END IF
    IF ( calc(2) ) THEN
       deriv => xc_dset_get_derivative(deriv_set, "(rho)(rho)", &
                                      allocate_deriv=.TRUE., error=error)
       CALL xc_derivative_get(deriv,deriv_data=e_rr,error=error)
       CALL pade_lda_2(n, rho_set%rho, rs, e_rr)
    END IF
    IF ( calc(3) ) THEN
       deriv => xc_dset_get_derivative(deriv_set, "(rho)(rho)(rho)", &
                                       allocate_deriv=.TRUE., error=error)
       CALL xc_derivative_get(deriv,deriv_data=e_rrr,error=error)
       CALL pade_lda_3(n, rho_set%rho, rs, e_rrr)
    END IF

    DEALLOCATE ( rs, STAT=stat )
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE pade_lda_pw_eval

! *****************************************************************************
  SUBROUTINE pade_lsd ( derivative_set, rhoa, rhob, pot, order )

    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rhoa, rhob
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: pot
    INTEGER, INTENT(IN), OPTIONAL            :: order

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pade_lsd', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, m, m1, m2, mp, mt, &
                                                mx(0:3,2), n
    LOGICAL                                  :: calc(0:4)
    REAL(KIND=dp)                            :: rs
    REAL(KIND=dp), DIMENSION(4)              :: fx

    CALL setup_calculation(order,mx,calc,201)

    n  = SIZE ( rhoa )
    mp = SIZE ( pot, 2 )
    m  = ABS(order)
    mt = mx(m,2)

    IF ( mt > mp ) CALL stop_program(routineN,moduleN,__LINE__,&
                                     "Size of POT array too small.")
    IF ( n > SIZE(pot,1) ) CALL stop_program(routineN,moduleN,__LINE__,&
                                             "Size of POT array too small.")

!$omp parallel do private(i,fx,rs,m1,m2)
    DO i=1, n
       fx(1) = rhoa(i) + rhob(i)

       CALL calc_rs ( fx(1), rs )
       CALL calc_fx ( rhoa(i), rhob(i), fx, m )

       IF ( calc(0) .AND. calc(1) ) THEN
          m1=mx(0,1)
          m2=mx(1,2)
          CALL pade_lsd_01 ( rhoa(i), rhob(i), rs, fx, &
                             pot(i,m1), pot(i,m1+1), pot(i,m2) )
       ELSE IF ( calc(0) ) THEN
          CALL pade_lsd_0 ( rhoa(i), rhob(i), rs, fx, pot(i,1) )
       ELSE IF ( calc(1) ) THEN
          m1=mx(1,1)
          m2=mx(1,2)
          CALL pade_lsd_1 ( rhoa(i), rhob(i), rs, fx, &
                            pot(i,m1), pot(i,m2) )
       END IF
       IF ( calc(2) ) THEN
          m1=mx(2,1)
          m2=mx(2,2)
          CALL pade_lsd_2 ( rhoa(i), rhob(i), rs, fx, &
                            pot(i,m1), pot(i,m1+1), pot(i,m2) )
       END IF
       IF ( calc(3) ) THEN
          m1=mx(3,1)
          m2=mx(3,2)
          CALL pade_lsd_3 ( rhoa(i), rhob(i), rs, fx, &
                            pot(i,m1), pot(i,m1+1), pot(i,m1+2), pot(i,m2) )
       END IF
    END DO

  END SUBROUTINE pade_lsd

! *****************************************************************************
  SUBROUTINE pade_lsd_pw_eval ( deriv_set, rho_set, order, error )

    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    INTEGER, INTENT(IN), OPTIONAL            :: order
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    INTEGER                                  :: i, j, k
    LOGICAL                                  :: calc(0:4)
    REAL(KIND=dp)                            :: rhoa, rhob, rs
    REAL(KIND=dp), DIMENSION(4)              :: fx
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: e_0, e_ra, e_rara, e_rarara, &
                                                e_rararb, e_rarb, e_rarbrb, &
                                                e_rb, e_rbrb, e_rbrbrb
    TYPE(xc_derivative_type), POINTER        :: deriv

    calc = .FALSE.
    IF (order >= 0) calc(0:order) = .TRUE.
    IF (order < 0) calc(-order) = .TRUE.

    IF (calc(0)) THEN
       deriv => xc_dset_get_derivative(deriv_set, "", &
                                       allocate_deriv=.TRUE., error=error)
       CALL xc_derivative_get(deriv,deriv_data=e_0,error=error)
    END IF
    IF (calc(1)) THEN
       deriv => xc_dset_get_derivative(deriv_set, "(rhoa)", &
                                       allocate_deriv=.TRUE., error=error)
       CALL xc_derivative_get(deriv,deriv_data=e_ra,error=error)
       deriv => xc_dset_get_derivative(deriv_set, "(rhob)", &
                                       allocate_deriv=.TRUE., error=error)
       CALL xc_derivative_get(deriv,deriv_data=e_rb,error=error)
    END IF
    IF (calc(2)) THEN
       deriv => xc_dset_get_derivative(deriv_set, "(rhoa)(rhoa)", &
                                       allocate_deriv=.TRUE., error=error)
       CALL xc_derivative_get(deriv,deriv_data=e_rara,error=error)
       deriv => xc_dset_get_derivative(deriv_set, "(rhoa)(rhob)", &
                                       allocate_deriv=.TRUE., error=error)
       CALL xc_derivative_get(deriv,deriv_data=e_rarb,error=error)
       deriv => xc_dset_get_derivative(deriv_set, "(rhob)(rhob)", &
                                       allocate_deriv=.TRUE., error=error)
       CALL xc_derivative_get(deriv,deriv_data=e_rbrb,error=error)
    END IF
    IF (calc(3)) THEN
       deriv => xc_dset_get_derivative(deriv_set, "(rhoa)(rhoa)(rhoa)", &
                                       allocate_deriv=.TRUE., error=error)
       CALL xc_derivative_get(deriv,deriv_data=e_rarara,error=error)
       deriv => xc_dset_get_derivative(deriv_set, "(rhoa)(rhoa)(rhob)", &
                                       allocate_deriv=.TRUE., error=error)
       CALL xc_derivative_get(deriv,deriv_data=e_rararb,error=error)
       deriv => xc_dset_get_derivative(deriv_set, "(rhoa)(rhob)(rhob)", &
                                       allocate_deriv=.TRUE., error=error)
       CALL xc_derivative_get(deriv,deriv_data=e_rarbrb,error=error)
       deriv => xc_dset_get_derivative(deriv_set, "(rhob)(rhob)(rhob)", &
                                       allocate_deriv=.TRUE., error=error)
       CALL xc_derivative_get(deriv,deriv_data=e_rbrbrb,error=error)
    END IF

!$omp parallel do private(i,j,k,fx,rhoa,rhob,rs)
    DO i=rho_set%local_bounds(1,1), rho_set%local_bounds(2,1)
       DO j=rho_set%local_bounds(1,2), rho_set%local_bounds(2,2)
          DO k=rho_set%local_bounds(1,3), rho_set%local_bounds(2,3)

             rhoa = rho_set%rhoa(i,j,k)
             rhob = rho_set%rhob(i,j,k)
             fx(1) = rhoa + rhob

             CALL calc_rs ( fx(1), rs )
             CALL calc_fx ( rhoa, rhob, fx, ABS(order) )

             IF ( calc(0) .AND. calc(1) ) THEN
                CALL pade_lsd_01 ( rhoa, rhob, rs, fx, &
                                   e_0(i,j,k), e_ra(i,j,k), e_rb(i,j,k) )
             ELSE IF ( calc(0) ) THEN
                CALL pade_lsd_0 ( rhoa, rhob, rs, fx, e_0(i,j,k) )
             ELSE IF ( calc(1) ) THEN
                CALL pade_lsd_1 ( rhoa, rhob, rs, fx, &
                                  e_ra(i,j,k), e_rb(i,j,k) )
             END IF
             IF ( calc(2) ) THEN
                CALL pade_lsd_2 ( rhoa, rhob, rs, fx, &
                                  e_rara(i,j,k), e_rarb(i,j,k), e_rbrb(i,j,k) )
             END IF
             IF ( calc(3) ) THEN
                CALL pade_lsd_3 ( rhoa, rhob, rs, fx, &
                                  e_rarara(i,j,k), e_rararb(i,j,k), e_rarbrb(i,j,k), e_rbrbrb(i,j,k) )
             END IF
          END DO
       END DO
    END DO

  END SUBROUTINE pade_lsd_pw_eval

! *****************************************************************************
  SUBROUTINE pade_lda_0(n, rho, rs, pot)

    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(*), INTENT(IN)  :: rho, rs
    REAL(KIND=dp), DIMENSION(*), &
      INTENT(INOUT)                          :: pot

    INTEGER                                  :: ip
    REAL(KIND=dp)                            :: epade, p, q

!$omp parallel do private(ip,p,q,epade)

    DO ip = 1, n
      IF ( rho(ip) > eps_rho ) THEN
         p = a0 + (a1 + (a2 + a3*rs(ip))*rs(ip))*rs(ip)
         q = (b1 + (b2 + (b3 + b4*rs(ip))*rs(ip))*rs(ip))*rs(ip)
         epade = -p/q
         pot(ip) = pot(ip) + epade*rho(ip)
      END IF
    END DO

  END SUBROUTINE pade_lda_0

! *****************************************************************************
  SUBROUTINE pade_lda_1(n, rho, rs, pot)

    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(*), INTENT(IN)  :: rho, rs
    REAL(KIND=dp), DIMENSION(*), &
      INTENT(INOUT)                          :: pot

    INTEGER                                  :: ip
    REAL(KIND=dp)                            :: depade, dpv, dq, epade, p, q

!$omp parallel do private(ip,p,q,epade,dpv,dq,depade)

    DO ip = 1, n
      IF ( rho(ip) > eps_rho ) THEN

         p = a0 + (a1 + (a2 + a3*rs(ip))*rs(ip))*rs(ip)
         q = (b1 + (b2 + (b3 + b4*rs(ip))*rs(ip))*rs(ip))*rs(ip)
         epade = -p/q

         dpv = a1 + (2.0_dp*a2 + 3.0_dp*a3*rs(ip))*rs(ip)
         dq = b1 + (2.0_dp*b2 + (3.0_dp*b3 + 4.0_dp*b4*rs(ip))*rs(ip))*rs(ip)
         depade = f13 * rs(ip) * (dpv*q - p*dq) / (q*q)

         pot(ip) = pot(ip) + epade + depade

      END IF
    END DO

  END SUBROUTINE pade_lda_1

! *****************************************************************************
  SUBROUTINE pade_lda_01(n, rho, rs, pot0, pot1)

    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(*), INTENT(IN)  :: rho, rs
    REAL(KIND=dp), DIMENSION(*), &
      INTENT(INOUT)                          :: pot0, pot1

    INTEGER                                  :: ip
    REAL(KIND=dp)                            :: depade, dpv, dq, epade, p, q

!$omp parallel do private(ip,p,q,epade,dpv,dq,depade)

    DO ip = 1, n
      IF ( rho(ip) > eps_rho ) THEN

         p = a0 + (a1 + (a2 + a3*rs(ip))*rs(ip))*rs(ip)
         q = (b1 + (b2 + (b3 + b4*rs(ip))*rs(ip))*rs(ip))*rs(ip)
         epade = -p/q

         dpv = a1 + (2.0_dp*a2 + 3.0_dp*a3*rs(ip))*rs(ip)
         dq = b1 + (2.0_dp*b2 + (3.0_dp*b3 + 4.0_dp*b4*rs(ip))*rs(ip))*rs(ip)
         depade = f13 * rs(ip) * (dpv*q - p*dq) / (q*q)

         pot0(ip) = pot0(ip) + epade * rho(ip)
         pot1(ip) = pot1(ip) + epade + depade

      END IF
    END DO

  END SUBROUTINE pade_lda_01

! *****************************************************************************
  SUBROUTINE pade_lda_2(n, rho, rs, pot)

    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(*), INTENT(IN)  :: rho, rs
    REAL(KIND=dp), DIMENSION(*), &
      INTENT(INOUT)                          :: pot

    INTEGER                                  :: ip
    REAL(KIND=dp)                            :: d2p, d2q, dpv, dq, p, q, rsr, &
                                                t1, t2, t3

!$omp parallel do private(ip,p,q,dpv,dq,d2p,d2q,rsr,t1,t2,t3)

    DO ip = 1, n
      IF ( rho(ip) > eps_rho ) THEN

         p = a0 + (a1 + (a2 + a3*rs(ip))*rs(ip))*rs(ip)
         q = (b1 + (b2 + (b3 + b4*rs(ip))*rs(ip))*rs(ip))*rs(ip)

         dpv = a1 + (2.0_dp*a2 + 3.0_dp*a3*rs(ip))*rs(ip)
         dq = b1 + (2.0_dp*b2 + (3.0_dp*b3 + 4.0_dp*b4*rs(ip))*rs(ip))*rs(ip)

         d2p = 2.0_dp*a2 + 6.0_dp*a3*rs(ip)
         d2q = 2.0_dp*b2 + (6.0_dp*b3 + 12.0_dp*b4*rs(ip))*rs(ip)

         rsr = rs(ip)/rho(ip)
         t1 = (p*dq - dpv*q)/(q*q)
         t2 = (d2p*q - p*d2q)/(q*q)
         t3 = (p*dq*dq - dpv*q*dq)/(q*q*q)

         pot(ip) = pot(ip) - f13 * ( f23*t1 + f13*t2*rs(ip) + f23*t3*rs(ip) ) * rsr

      END IF
    END DO

  END SUBROUTINE pade_lda_2

! *****************************************************************************
  SUBROUTINE pade_lda_3(n, rho, rs, pot)

    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(*), INTENT(IN)  :: rho, rs
    REAL(KIND=dp), DIMENSION(*), &
      INTENT(INOUT)                          :: pot

    INTEGER                                  :: ip
    REAL(KIND=dp)                            :: ab1, ab2, ab3, d2p, d2q, d3p, &
                                                d3q, dpv, dq, p, q, rsr1, &
                                                rsr2, rsr3

!$omp parallel do private(ip,p,q,dpv,dq,d2p,d2q,d3p,d3q,ab1,ab2,ab3,rsr1,rsr2,rsr3)

    DO ip = 1, n
      IF ( rho(ip) > eps_rho ) THEN

         p = a0 + (a1 + (a2 + a3*rs(ip))*rs(ip))*rs(ip)
         q = (b1 + (b2 + (b3 + b4*rs(ip))*rs(ip))*rs(ip))*rs(ip)

         dpv = a1 + (2.0_dp*a2 + 3.0_dp*a3*rs(ip))*rs(ip)
         dq = b1 + (2.0_dp*b2 + (3.0_dp*b3 + 4.0_dp*b4*rs(ip))*rs(ip))*rs(ip)

         d2p = 2.0_dp*a2 + 6.0_dp*a3*rs(ip)
         d2q = 2.0_dp*b2 + (6.0_dp*b3 + 12.0_dp*b4*rs(ip))*rs(ip)

         d3p = 6.0_dp*a3
         d3q = 6.0_dp*b3 + 24.0_dp*b4*rs(ip)

         ab1 = (dpv*q - p*dq)/(q*q)
         ab2 = (d2p*q*q - p*q*d2q - 2.0_dp*dpv*q*dq + 2.0_dp*p*dq*dq)/(q*q*q)
         ab3 = (d3p*q*q - p*q*d3q - 3.0_dp*dpv*q*d2q + 3.0_dp*p*dq*d2q)/(q*q*q)
         ab3 = ab3 - 3.0_dp*ab2*dq/q
         rsr1 = rs(ip) / (rho(ip)*rho(ip))
         rsr2 = f13*f13*rs(ip) * rsr1
         rsr3 = f13*rs(ip) * rsr2
         rsr1 = -f23*f23*f23*rsr1
         pot(ip) = pot(ip) + rsr1*ab1 + rsr2*ab2 + rsr3*ab3

      END IF
    END DO

  END SUBROUTINE pade_lda_3

! *****************************************************************************
  SUBROUTINE pade_lsd_0(rhoa, rhob, rs, fx, pot0 )

    REAL(KIND=dp), INTENT(IN)                :: rhoa, rhob, rs
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: fx
    REAL(KIND=dp), INTENT(INOUT)             :: pot0

    REAL(KIND=dp)                            :: fa0, fa1, fa2, fa3, fb1, fb2, &
                                                fb3, fb4, p, q, rhoab

      rhoab = rhoa + rhob

      IF ( rhoab > eps_rho ) THEN

         fa0 = a0 + fx(1)*da0
         fa1 = a1 + fx(1)*da1
         fa2 = a2 + fx(1)*da2
         fa3 = a3 + fx(1)*da3
         fb1 = b1 + fx(1)*db1
         fb2 = b2 + fx(1)*db2
         fb3 = b3 + fx(1)*db3
         fb4 = b4 + fx(1)*db4

         p = fa0 + (fa1 + (fa2 + fa3*rs)*rs)*rs
         q = (fb1 + (fb2 + (fb3 + fb4*rs)*rs)*rs)*rs

         pot0 = pot0 - p/q*rhoab

      END IF

  END SUBROUTINE pade_lsd_0

! *****************************************************************************
  SUBROUTINE pade_lsd_1 ( rhoa, rhob, rs, fx, pota, potb )

    REAL(KIND=dp), INTENT(IN)                :: rhoa, rhob, rs
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: fx
    REAL(KIND=dp), INTENT(INOUT)             :: pota, potb

    REAL(KIND=dp)                            :: dc, dpv, dq, dr, dx, fa0, &
                                                fa1, fa2, fa3, fb1, fb2, fb3, &
                                                fb4, p, q, rhoab, xp, xq

      rhoab = rhoa + rhob

      IF ( rhoab > eps_rho ) THEN

         fa0 = a0 + fx(1)*da0
         fa1 = a1 + fx(1)*da1
         fa2 = a2 + fx(1)*da2
         fa3 = a3 + fx(1)*da3
         fb1 = b1 + fx(1)*db1
         fb2 = b2 + fx(1)*db2
         fb3 = b3 + fx(1)*db3
         fb4 = b4 + fx(1)*db4

         p = fa0 + (fa1 + (fa2 + fa3*rs)*rs)*rs
         q = (fb1 + (fb2 + (fb3 + fb4*rs)*rs)*rs)*rs
         dpv = fa1 + (2.0_dp*fa2 + 3.0_dp*fa3*rs)*rs
         dq = fb1 + (2.0_dp*fb2 + (3.0_dp*fb3 + &
                     4.0_dp*fb4*rs)*rs)*rs
         xp = da0 + (da1 + (da2 + da3*rs)*rs)*rs
         xq = (db1 + (db2 + (db3 + db4*rs)*rs)*rs)*rs

         dr = (dpv*q - p*dq)/(q*q)
         dx = 2.0_dp*(xp*q - p*xq)/(q*q)*fx(2)/rhoab
         dc = f13*rs*dr - p/q

         pota = pota + dc - dx*rhob
         potb = potb + dc + dx*rhoa

      END IF

  END SUBROUTINE pade_lsd_1

! *****************************************************************************
  SUBROUTINE pade_lsd_01 ( rhoa, rhob, rs, fx, pot0, pota, potb )

    REAL(KIND=dp), INTENT(IN)                :: rhoa, rhob, rs
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: fx
    REAL(KIND=dp), INTENT(INOUT)             :: pot0, pota, potb

    REAL(KIND=dp)                            :: dc, dpv, dq, dr, dx, fa0, &
                                                fa1, fa2, fa3, fb1, fb2, fb3, &
                                                fb4, p, q, rhoab, xp, xq

      rhoab = rhoa + rhob

      IF ( rhoab > eps_rho ) THEN

         fa0 = a0 + fx(1)*da0
         fa1 = a1 + fx(1)*da1
         fa2 = a2 + fx(1)*da2
         fa3 = a3 + fx(1)*da3
         fb1 = b1 + fx(1)*db1
         fb2 = b2 + fx(1)*db2
         fb3 = b3 + fx(1)*db3
         fb4 = b4 + fx(1)*db4

         p = fa0 + (fa1 + (fa2 + fa3*rs)*rs)*rs
         q = (fb1 + (fb2 + (fb3 + fb4*rs)*rs)*rs)*rs
         dpv = fa1 + (2.0_dp*fa2 + 3.0_dp*fa3*rs)*rs
         dq = fb1 + (2.0_dp*fb2 + (3.0_dp*fb3 + &
                     4.0_dp*fb4*rs)*rs)*rs
         xp = da0 + (da1 + (da2 + da3*rs)*rs)*rs
         xq = (db1 + (db2 + (db3 + db4*rs)*rs)*rs)*rs

         dr = (dpv*q - p*dq)/(q*q)
         dx = 2.0_dp*(xp*q - p*xq)/(q*q)*fx(2)/rhoab
         dc = f13*rs*dr - p/q

         pot0 = pot0 - p/q*rhoab
         pota = pota + dc - dx*rhob
         potb = potb + dc + dx*rhoa

      END IF

  END SUBROUTINE pade_lsd_01

! *****************************************************************************
  SUBROUTINE pade_lsd_2 (rhoa, rhob, rs, fx, potaa, potab, potbb )

    REAL(KIND=dp), INTENT(IN)                :: rhoa, rhob, rs
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: fx
    REAL(KIND=dp), INTENT(INOUT)             :: potaa, potab, potbb

    REAL(KIND=dp)                            :: d2p, d2q, dpv, dq, dr, drr, &
                                                dx, dxp, dxq, dxr, dxx, fa0, &
                                                fa1, fa2, fa3, fb1, fb2, fb3, &
                                                fb4, or, p, q, rhoab, xp, xq, &
                                                xt, yt

      rhoab = rhoa + rhob

      IF ( rhoab > eps_rho ) THEN

         fa0 = a0 + fx(1)*da0
         fa1 = a1 + fx(1)*da1
         fa2 = a2 + fx(1)*da2
         fa3 = a3 + fx(1)*da3
         fb1 = b1 + fx(1)*db1
         fb2 = b2 + fx(1)*db2
         fb3 = b3 + fx(1)*db3
         fb4 = b4 + fx(1)*db4

         p = fa0 + (fa1 + (fa2 + fa3*rs)*rs)*rs
         q = (fb1 + (fb2 + (fb3 + fb4*rs)*rs)*rs)*rs

         dpv = fa1 + (2.0_dp*fa2 + 3.0_dp*fa3*rs)*rs
         dq = fb1 + (2.0_dp*fb2 + (3.0_dp*fb3 + &
                     4.0_dp*fb4*rs)*rs)*rs

         d2p = 2.0_dp*fa2 + 6.0_dp*fa3*rs
         d2q = 2.0_dp*fb2 + (6.0_dp*fb3 + 12.0_dp*fb4*rs)*rs

         xp = da0 + (da1 + (da2 + da3*rs)*rs)*rs
         xq = (db1 + (db2 + (db3 + db4*rs)*rs)*rs)*rs

         dxp = da1 + (2.0_dp*da2 + 3.0_dp*da3*rs)*rs
         dxq = db1 + (2.0_dp*db2 + (3.0_dp*db3 + &
                     4.0_dp*db4*rs)*rs)*rs

         dr = (dpv*q - p*dq)/(q*q)
         drr = (d2p*q*q-p*q*d2q-2.0_dp*dpv*q*dq+2.0_dp*p*dq*dq)/(q*q*q)
         dx = (xp*q - p*xq)/(q*q)
         dxx = 2.0_dp*xq*(p*xq - xp*q)/(q*q*q)
         dxr = (dxp*q*q+dpv*xq*q-xp*dq*q-p*dxq*q-2.0_dp*dpv*q*xq+2.0_dp*p*dq*xq)/(q*q*q)

         or = 1.0_dp/rhoab
         yt = rhob*or
         xt = rhoa*or

         potaa = potaa + f23*f13*dr*rs*or - f13*f13*drr*rs*rs*or &
                     + f43*rs*fx(2)*dxr*yt*or &
                     - 4.0_dp*fx(2)*fx(2)*dxx*yt*yt*or &
                     - 4.0_dp*dx*fx(3)*yt*yt*or
         potab = potab + f23*f13*dr*rs*or - f13*f13*drr*rs*rs*or &
                     + f23*rs*fx(2)*dxr*(yt-xt)*or &
                     + 4.0_dp*fx(2)*fx(2)*dxx*xt*yt*or &
                     + 4.0_dp*dx*fx(3)*xt*yt*or
         potbb = potbb + f23*f13*dr*rs*or - f13*f13*drr*rs*rs*or &
                     - f43*rs*fx(2)*dxr*xt*or &
                     - 4.0_dp*fx(2)*fx(2)*dxx*xt*xt*or &
                     - 4.0_dp*dx*fx(3)*xt*xt*or

      END IF

  END SUBROUTINE pade_lsd_2

! *****************************************************************************
  SUBROUTINE pade_lsd_3 ( rhoa, rhob, rs, fx, potaaa, potaab, potabb, potbbb )

    REAL(KIND=dp), INTENT(IN)                :: rhoa, rhob, rs
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: fx
    REAL(KIND=dp), INTENT(INOUT)             :: potaaa, potaab, potabb, potbbb

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pade_lsd_3', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp) :: d2p, d2q, d2xp, d2xq, d3p, d3q, dpv, dq, dr, drr, drrr, &
      dx, dxp, dxq, dxr, dxrr, dxx, dxxr, dxxx, fa0, fa1, fa2, fa3, fb1, fb2, &
      fb3, fb4, or, p, q, rhoab, xp, xq, xt, yt

    IF ( .NOT. debug_flag ) CALL stop_program(routineN,moduleN,__LINE__,&
                                              "Routine not tested")

      rhoab = rhoa + rhob

      IF ( rhoab > eps_rho ) THEN

         fa0 = a0 + fx(1)*da0
         fa1 = a1 + fx(1)*da1
         fa2 = a2 + fx(1)*da2
         fa3 = a3 + fx(1)*da3
         fb1 = b1 + fx(1)*db1
         fb2 = b2 + fx(1)*db2
         fb3 = b3 + fx(1)*db3
         fb4 = b4 + fx(1)*db4

         p = fa0 + (fa1 + (fa2 + fa3*rs)*rs)*rs
         q = (fb1 + (fb2 + (fb3 + fb4*rs)*rs)*rs)*rs

         dpv = fa1 + (2.0_dp*fa2 + 3.0_dp*fa3*rs)*rs
         dq = fb1 + (2.0_dp*fb2 + (3.0_dp*fb3 + &
                     4.0_dp*fb4*rs)*rs)*rs

         d2p = 2.0_dp*fa2 + 6.0_dp*fa3*rs
         d2q = 2.0_dp*fb2 + (6.0_dp*fb3 + 12.0_dp*fb4*rs)*rs

         d3p = 6.0_dp*fa3
         d3q = 6.0_dp*fb3 + 24.0_dp*fb4*rs

         xp = da0 + (da1 + (da2 + da3*rs)*rs)*rs
         xq = (db1 + (db2 + (db3 + db4*rs)*rs)*rs)*rs

         dxp = da1 + (2.0_dp*da2 + 3.0_dp*da3*rs)*rs
         dxq = db1 + (2.0_dp*db2 + (3.0_dp*db3 + &
                     4.0_dp*db4*rs)*rs)*rs

         d2xp = 2.0_dp*da2 + 6.0_dp*da3*rs
         d2xq = 2.0_dp*db2 + (6.0_dp*db3 + 12.0_dp*db4*rs)*rs

         dr = (dpv*q - p*dq)/(q*q)
         drr = (d2p*q*q-p*q*d2q-2.0_dp*dpv*q*dq+2.0_dp*p*dq*dq)/(q*q*q)
         drrr = (d3p*q*q*q-3.0_dp*d2p*dq*q*q+6.0_dp*dpv*dq*dq*q-3.0_dp*dpv*d2q*q*q - &
                6.0_dp*p*dq*dq*dq+6.0_dp*p*dq*d2q*q-p*d3q*q*q)/(q*q*q*q)
         dx = (xp*q - p*xq)/(q*q)
         dxx = 2.0_dp*xq*(p*xq - xp*q)/(q*q*q)
         dxxx = 6.0_dp*xq*(q*xp*xq-p*xq*xq)/(q*q*q*q)
         dxr = (dxp*q*q+dpv*xq*q-xp*dq*q-p*dxq*q-2.0_dp*dpv*q*xq+2.0_dp*p*dq*xq)/(q*q*q)
         dxxr = 2.0_dp*(2.0_dp*dxq*q*p*xq-dxq*q*q*xp+xq*xq*q*dpv-xq*q*q*dxp + &
                2.0_dp*xq*q*xp*dq-3.0_dp*xq*xq*dq*p)/(q*q*q*q)
         dxrr = (q*q*q*d2xp-2.0_dp*q*q*dxp*dq-q*q*xp*d2q-q*q*d2p*xq - &
                2.0_dp*q*q*dpv*dxq-q*q*p*d2xq+4.0_dp*dq*q*dpv*xq+4.0_dp*dq*q*p*dxq + &
                2.0_dp*dq*dq*q*xp-6.0_dp*dq*dq*p*xq+2.0_dp*d2q*q*p*xq)/(q*q*q*q)

         or = 1.0_dp/rhoab
         yt = rhob*or
         xt = rhoa*or

         potaaa = potaaa + 8.0_dp/27.0_dp*dr*rs*or*or + &
                     1.0_dp/9.0_dp*drr*rs*rs*or*or + &
                     1.0_dp/27.0_dp*drrr*rs**3*or*or + &
                     dxr*or*or*yt*rs*(-8.0_dp/3.0_dp*fx(2)+4.0_dp*fx(3)*yt)
         potaab = potaab + 0.0_dp
         potabb = potabb + 0.0_dp
         potbbb = potbbb + 0.0_dp

      END IF

  END SUBROUTINE pade_lsd_3

END MODULE xc_pade

