MODULE atom_calculations

  USE atom_type, ONLY : atom_info, atom_init
  USE kinds, ONLY : dbl
  USE termination, ONLY : stop_program
  USE radial_util, ONLY : radial_integrate
  USE radial_grids, ONLY : radial_grid_type

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: &
       calculateElectronDensity, &
       getGroundStateEnergy, &
       calculatePotential, &
       countNodes, &
       initialEigenvalues, &
       rgNormalize, &
       signum

CONTAINS

! **************************************************************************

  SUBROUTINE calculatePotential(v, n, l, a, rho)

    USE functionals, ONLY : slater
    USE mathconstants, ONLY : pi
    USE radial_poisson, ONLY : rg_poisson
    
    IMPLICIT NONE

    ! ARGUMENTS
    REAL(DBL), DIMENSION(:), INTENT(OUT) :: v  ! v to compute
    INTEGER, INTENT(IN) :: n, l                ! main and angular quantum number
    TYPE(atom_info), INTENT(IN) :: a           ! info on the atom to be computed
    REAL(DBL), DIMENSION(:), INTENT(IN) :: rho ! the electron density

    ! VARIABLES
    INTEGER :: errno, np, k

    ! needed for XA
    REAL(DBL), DIMENSION(0:a%rg%n+1) :: y0  ! r times coulomb repulsion potential
    REAL(DBL), DIMENSION(0:a%rg%n+1) :: rho_calc ! electron density (for calc)
    REAL(DBL) :: ex, vx

    ! ----------------------------------------------------------
    np = a%rg%n

    ! zentrifugal term, common to all potentials
    v(1:np) = DBLE(l*(l+1)) / (2.0_dbl * (a%rg%r(1:np))**2)

    SELECT CASE (a%method)

    CASE DEFAULT

       CALL stop_program(&
            "atom_calculations/calculatePotential","Method not implemented")

    CASE ('XA')

       ! calculate the electrostatic electron potential
       rho_calc = 0.0_dbl
       rho_calc(1:np) = rho(1:np) * (4.0_dbl*pi * (a%rg%r(1:np))**2)
       CALL rg_poisson(rho_calc, y0, 0, a%rg, 'NUMEROV')

       v(1:np) = v(1:np) &                 
            - DBLE(a%z) / (a%rg%r(1:np)) &  ! el-core interaction
            + y0(1:np) / (a%rg%r(1:np))     ! electrostatic el-el repulsion

       ! calculate and add the exchange part
       DO k=1, np
          CALL slater(rho(k), ex, vx)
          v(k) = v(k) + vx
       ENDDO


    CASE ('SIMPLE')
       v(1:np) = v(1:np) - DBLE(a%z) / (a%rg%r(1:np))

    END SELECT

  END SUBROUTINE calculatePotential

! **************************************************************************
  
  SUBROUTINE calculateElectronDensity(n, phi, c, w)

    ! Calculates the electron density (n), given a set of wave functions (phi).
    ! c specifies the number of electrons in an orbit and w is a possible weight
    ! (in r-space)

    IMPLICIT NONE

    ! arguments
    REAL(DBL), DIMENSION(1:), INTENT(OUT) :: n          ! the density to be calculated
    REAL(DBL), DIMENSION(1:,1:), INTENT(IN) :: phi      ! the wavefunctions
    REAL(DBL), DIMENSION(1:), INTENT(IN) :: c           ! number of electrons in the orbit
    REAL(DBL), DIMENSION(1:), OPTIONAL, INTENT(IN) :: w ! weights for the wavefunctions
    
    ! variables
    INTEGER :: i                                        ! loop variable
    REAL(DBL), DIMENSION(1:SIZE(n),1:SIZE(c)) :: psi    ! the wavefunctions to work with
    
    IF (PRESENT(w)) THEN
       DO i=1, size(w)
          psi(i,:) = phi(i,:)*w(i)
       ENDDO
    ELSE
       psi = phi
    ENDIF

    n = 0.0_dbl
    DO i=1, SIZE(c)
       IF (c(i) /= 0) n = n + c(i)*psi(:,i)**2
    ENDDO          

  END SUBROUTINE calculateElectronDensity

! **************************************************************************

  SUBROUTINE rgNormalize(f, rg)

    USE mathconstants, ONLY : pi

    IMPLICIT NONE

    REAL(DBL), DIMENSION(1:), INTENT(INOUT) :: f
    TYPE(radial_grid_type), INTENT(IN) :: rg

    REAL(DBL) :: res
    INTEGER :: errno

    CALL radial_integrate(rg%n, 4.0_dbl*pi * f**2, rg, res, errno)
    f = f/SQRT(res)

  END SUBROUTINE rgNormalize

! **************************************************************************

  SUBROUTINE countNodes(nodes, f)

    IMPLICIT NONE

    INTEGER, INTENT(OUT) :: nodes
    REAL(DBL), DIMENSION(:), INTENT(IN) :: f
    
    INTEGER :: n, i
    INTEGER :: s

    n = SIZE(f)
    IF ((n < 4) .OR. (f(1) < 0.0_dbl)) THEN ! not enough or bad function data
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
    
  END SUBROUTINE countNodes

! **************************************************************************

  SUBROUTINE initialEigenvalues(ew, a)
    
    IMPLICIT NONE
    REAL(DBL), DIMENSION(:), INTENT(OUT) :: ew ! the vector with the eigenvalues
    TYPE(atom_info), INTENT(IN) :: a        ! 

    ! variables
    INTEGER :: n, l, orbit                     ! main and angular quantum numbers
    INTEGER :: counter

    counter = 0
    DO n=1, 10
       DO l=0, n-1

          orbit = (n*(n-1))/2 + 1 + l

          IF (a%elec_conf(orbit) == 0) CYCLE

          IF (a%method == 'SIMPLE') THEN
             ew(orbit) = -0.5_dbl * (DBLE(a%z)/DBLE(n))**2
          ELSE
             ew(orbit) = -0.5_dbl * (DBLE(a%z - counter)/DBLE(n))**2
          END IF
          ew(orbit) = ew(orbit) * 1.05_dbl ! shift the guesses by 5%

       END DO
       counter = counter + SUM(a%elec_conf((n*(n-1)/2 + 1):(n*(n-1)/2 + n)))
    END DO
  END SUBROUTINE initialEigenvalues

! **************************************************************************

  SUBROUTINE getGroundStateEnergy(a, rho)

    USE mathconstants, ONLY : pi
    USE radial_poisson, ONLY : rg_poisson

    IMPLICIT NONE

    ! ARGUMENTS
    TYPE(atom_info), INTENT(INOUT) :: a
    REAL(DBL), DIMENSION(:), INTENT(IN) :: rho

    ! VARIABLES
    REAL(DBL) :: res
    REAL(DBL), DIMENSION(0:a%rg%n+1) :: y0, rho_calc
    INTEGER :: errno, np

    ! ------------------------------------------------------------
    np = a%rg%n

    a%ground_state_energy = SUM( DBLE(a%elec_conf) * a%eigenvalues )

    SELECT CASE (a%method)
       
    CASE DEFAULT
       CALL stop_program(&
            "atom_calculations/getGroundStateEnergy","Method not implemented")
       
    CASE ('SIMPLE')
       !! nothing to do

    CASE ('XA')
       rho_calc = 0.0_dbl

       rho_calc(1:np) = rho(1:np) * 4.0_dbl*pi * (a%rg%r(1:np))**2
       CALL rg_poisson(rho_calc, y0, 0, a%rg, 'NUMEROV')
       
       rho_calc(1:np) = rho * y0(1:np) * 4.0_dbl*pi * (a%rg%r(1:np))
       CALL radial_integrate(np, rho_calc(1:np), a%rg, res, errno)
       a%ground_state_energy = a%ground_state_energy - 0.5_dbl*res

       rho_calc(1:np) = rho**(4._dbl/3._dbl) * (4.0_dbl*pi*a%rg%r(1:np)**2)
       CALL radial_integrate(np, rho_calc(1:np), a%rg, res, errno)
       a%ground_state_energy = a%ground_state_energy + 0.25_dbl * a%mpara(2) * res

    END SELECT

  END SUBROUTINE getGroundStateEnergy

! **************************************************************************

  FUNCTION signum(a)
    
    IMPLICIT NONE

    INTEGER :: signum    
    REAL(DBL), INTENT(IN) :: a

    IF (a < 0.0_dbl) THEN 
       signum = -1
    ELSE
       signum = 1
    ENDIF    

  END FUNCTION signum

! **************************************************************************

END MODULE atom_calculations
