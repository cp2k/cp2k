!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****h* cp2k/atom_util *
!!
!!   NAME
!!     atom_util
!!
!!   FUNCTION
!!     
!!
!!   SYNOPSIS
!!     
!!
!!   MODIFICATION HISTORY
!!     , TCH, working version
!!
!!   SEE ALSO
!!     
!!   SOURCE
!!
MODULE atom_util

  USE atom_types,                      ONLY: atom_calculation_type,&
                                             atom_type
  USE xc_xalpha_new,                   ONLY: xalpha_lda_0
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: PI
  USE radial_poisson,                  ONLY: rg_poisson
  USE radial_util,                     ONLY: radial_integrate
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: calculate_energy, &
            calculate_potential, &
            calculate_rho, &
            count_nodes, &
            DGESS, &
            find_matchp, &
            signum

CONTAINS  

! ****************************************************************************

  FUNCTION calculate_energy ( atom, ac, rho ) RESULT (res)


    TYPE(atom_type), INTENT(IN)              :: atom
    TYPE(atom_calculation_type), INTENT(IN)  :: ac
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)      :: rho
    REAL(KIND=dp)                                :: res

    INTEGER                                  :: errno, i, l, n, np
    REAL(KIND=dp)                                :: res_int
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)     :: t, y0

! FUNCTION
! VARIABLES
! ------------------------------------------------------------

    np = ac%rg%n

    DO n=1, atom%ns
       DO l=0, atom%shells(n)%no - 1
          res = res + REAL(atom%shells(n)%orbits(l)%ne, KIND=dp) * atom%shells(n)%orbits(l)%ev
       END DO
    END DO
    
    SELECT CASE ( ac%method )
       
    CASE DEFAULT
       
    CASE ( 'SIMPLE')
       !! nothing to do

    CASE ( 'XA' )
       ALLOCATE ( t(0:np+1), y0(0:np+1) )
       
       t(0) = 0.0_dp; t(np+1) = 0.0_dp
       y0 = 0.0_dp

       DO i=1, np
          t(i) = rho(i) * ac%rg%r(i)**2
       END DO
       CALL DSCAL(np, 4.0_dp*pi, t(1:np), 1)
       CALL rg_poisson(t, y0, 0, ac%rg, 'NUMEROV')
       
       DO i=1, np
          t(i) = rho(i) * y0(i) * ac%rg%r(i)
       END DO
       CALL radial_integrate(np, t(1:np), ac%rg, res_int, errno)
       res = res - 0.5_dp*4.0_dp*pi*res_int

       DO i=1, np
          t(i) = rho(i)**(4.0_dp/3.0_dp) * ac%rg%r(i)**2
       END DO
       CALL radial_integrate(np, t(1:np), ac%rg, res_int, errno)
       res = res + 0.25_dp*4.0_dp*pi*res_int*ac%xapara%c
       
       DEALLOCATE(t, y0)
                               
    END SELECT

  END FUNCTION calculate_energy

! ****************************************************************************
      
  FUNCTION calculate_rho( atom, np )
!!! ***********************************************************
!!! Calculate the electronic density.
!!! Returns the full density (with angular function) and
!!! assumes that the wavefunction is the full wavefunction.
!!! ***********************************************************
    TYPE(ATOM_TYPE), INTENT(INOUT)           :: atom
    INTEGER, INTENT(IN)                      :: np
    REAL(KIND=dp), DIMENSION(1:np)               :: calculate_rho

    INTEGER                                  :: k, l, n

! cannot use DSCAL here

    calculate_rho = 0.0_dp

    DO n=1, atom%ns
       DO l=0, atom%shells(n)%no-1
          DO k = 1, np
             calculate_rho(k) = calculate_rho(k) + &
                  REAL(atom%shells(n)%orbits(l)%ne, KIND=dp) * &
                  (atom%shells(n)%orbits(l)%wf(k))**2
          END DO
       END DO
    END DO
    
  END FUNCTION calculate_rho

! ****************************************************************************
  
  FUNCTION find_matchp(pot, ev)
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)      :: pot
    REAL(KIND=dp), INTENT(IN)                    :: ev
    INTEGER                                  :: find_matchp

    INTEGER                                  :: k

    find_matchp = SIZE(pot) / 3

    ! to find best matching point: find point with v - e = 0
    DO k=2, SIZE(pot)
       IF ( signum(pot(k) - ev) /= signum(pot(k-1) - ev) ) THEN
          find_matchp = k
          EXIT
       ENDIF
    ENDDO

  END FUNCTION find_matchp

! ****************************************************************************
  
  SUBROUTINE calculate_potential(pot, n, l, rho, atom, ac)
!!! arguments
    REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)   :: pot
    INTEGER, INTENT(IN)                      :: n, l
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)      :: rho
    TYPE(ATOM_TYPE), INTENT(IN)              :: atom
    TYPE(ATOM_CALCULATION_TYPE), INTENT(IN)  :: ac

    INTEGER                                  :: i, np
    REAL(KIND=dp)                                :: ex, vx, z_eff
    REAL(KIND=dp), DIMENSION(0:ac%rg%n+1)        :: rc, y0
    REAL(KIND=dp), DIMENSION(size(rho))          :: r13

!!! variables

    np = ac%rg%n

!!! zentrifugal term
    pot(1:np) = REAL(l*(l+1),KIND=dp) / (2.0_dp*ac%rg%r(1:np)**2)

!!!!!! method specific
    SELECT CASE (ac%method)

    CASE('XA')
       rc = 0.0_dp; y0 = 0.0_dp
       rc(1:np) = rho(1:np) * ac%rg%r(1:np)**2
       CALL DSCAL(np, 4.0_dp*pi, rc(1:np), 1)
       CALL rg_poisson(rc, y0, 0, ac%rg, 'NUMEROV')
       z_eff = REAL(atom%z, KIND=dp)
       DO i=1, np
          r13(i)=rho(i)**(1.0_dp/3.0_dp)
       END DO       
       call xalpha_lda_0(np, rho, r13, pot)
       DO i=1, np
          pot(i) = pot(i) + (y0(i)-z_eff)/ac%rg%r(i)
       END DO       

    CASE ('SIMPLE')
       z_eff = REAL(atom%z, KIND=dp)
       pot(1:np) = pot(1:np) - z_eff / ac%rg%r(1:np)

    CASE DEFAULT
       WRITE (6,*) "Unkown method"

    END SELECT


  END SUBROUTINE calculate_potential

! **************************************************************************

  FUNCTION signum(a)
    

    REAL(KIND=dp), INTENT(IN)                    :: a
    INTEGER                                  :: signum

    IF (a < 0.0_dp) THEN 
       signum = -1
    ELSE
       signum = 1
    ENDIF    

  END FUNCTION signum

! **************************************************************************

  SUBROUTINE count_nodes(nodes, f)


    INTEGER, INTENT(OUT)                     :: nodes
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)      :: f

    INTEGER                                  :: i, n, s

    n = SIZE(f)
    IF ((n < 4) .OR. (f(1) < 0.0_dp)) THEN ! not enough or bad function data
       nodes = -1
       RETURN
    END IF
    
    nodes = 0
    s = signum(f(1))
    DO i=2, n-1

       IF (signum(f(i)) /= s) THEN
          s = signum(f(i))
          IF ( ABS(f(i+1)) > ABS(f(i)) .AND. &
               signum(f(i+1)) == s ) nodes = nodes + 1
       ENDIF

    ENDDO
    
  END SUBROUTINE count_nodes

! ****************************************************************************
      
  SUBROUTINE DGESS(ALPHA, X)
    REAL(KIND=dp), INTENT(IN)                    :: ALPHA
    REAL(KIND=dp), INTENT(INOUT)                 :: X

    X = ALPHA * X
  END SUBROUTINE DGESS
 
! **************************************************************************

END MODULE atom_util
!!*****
