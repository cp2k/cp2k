!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****m* cp2k/atom_calculations *
!!
!!   NAME
!!     atom_calculations
!!
!!   FUNCTION
!!     Contains functions and subroutines needed by atom_module
!!
!!   SYNOPSIS
!!     calculate_electron_density ( n, phi, c, w )
!!     calculate_gs_energy ( atom, parameters, rho )
!!     calculate_potential ( atom, parameters, v, n, l, rho )
!!     count_nodes ( nodes, f )
!!     initialize_eigenvalues ( atom, parameters, eigenvalues )
!!     rg_normalize ( f, rg )
!!     signum ( a )
!!
!!   MODIFICATION HISTORY
!!     14-FEB-2002, TCH, working version
!!
!!   SEE ALSO
!!     atom_module, atom_parameters_type, atom_type
!!
!!****
MODULE atom_calculations

  USE atom_parameters_type, ONLY : AM_PARAMETERS_INFO, &
       am_get_method, am_get_method_parameters
  USE atom_type, ONLY : AM_ATOM_INFO, am_get_elec_conf
  USE functionals, ONLY : slater
  USE kinds, ONLY : dbl
  USE mathconstants, ONLY : pi
  USE radial_poisson, ONLY : rg_poisson
  USE radial_util, ONLY : radial_integrate
  USE termination, ONLY : stop_program

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: &
       calculate_electron_density, &
       calculate_gs_energy, &
       calculate_potential, &
       count_nodes, &
       initialize_eigenvalues, &
       rg_normalize, &
       signum

CONTAINS

! **************************************************************************

!!****f* cp2k/atom_calculations/calculate_potential *
!!
!!   NAME
!!     calculate_potential
!!
!!   FUNCTION
!!     Calculate the potential term of the Kohn-Sham equations.
!!     Currently works for the most simple (independend, no pseudopotential)
!!     and the X-alpha method of Slater. Uses the functionals defined in 
!!     functionals.F.
!!
!!   ARGUMENTS
!!     v - REAL ( DBL ), DIMENSION - the potential to be calculated
!!     atom - TYPE ( AM_ATOM_INFO ) - the atom itself
!!     parameters - TYPE ( AM_PARAMETERS_INFO ) - parameters about the method
!!     n - INTEGER - main quantum number
!!     l - INTEGER - angular quantum number
!!     rho - REAL ( DBL ), DIMENSION - electron density
!!
!!   MODIFICATION HISTORY
!!     14-FEB-2002, TCH, working version
!!
!!   SEE ALSO
!!     atom_parameters_type, atom_type, functionals.F
!!
!!****
  SUBROUTINE calculate_potential ( v, atom, parameters, n, l, rho )

    IMPLICIT NONE

    ! ARGUMENTS
    REAL ( DBL ), DIMENSION ( : ), INTENT ( OUT ) :: v
    TYPE ( AM_ATOM_INFO ), INTENT ( IN ) :: atom
    TYPE ( AM_PARAMETERS_INFO ), INTENT ( IN ) :: parameters
    INTEGER, INTENT ( IN ) :: n, l                
    REAL ( DBL ), DIMENSION ( : ), INTENT ( IN ) :: rho

    ! VARIABLES
    INTEGER :: errno, np, k

    ! needed for XA
    REAL ( DBL ), DIMENSION ( 0:atom%rg%n+1 ) :: y0  ! r times coulomb repulsion potential
    REAL ( DBL ), DIMENSION ( 0:atom%rg%n+1 ) :: rho_calc ! electron density (for calc)
    REAL ( DBL ) :: ex, vx

    ! ----------------------------------------------------------
    np = atom%rg%n

    ! zentrifugal term, common to all potentials
    v(1:np) = DBLE(l*(l+1)) / (2.0_dbl * ( atom%rg%r(1:np))**2)

    SELECT CASE ( am_get_method ( parameters ) )

    CASE DEFAULT

       CALL stop_program(&
            "atom_calculations/calculate_potential","Method not implemented")

    CASE ('XA')

       ! calculate the electrostatic electron potential
       rho_calc = 0.0_dbl
       rho_calc(1:np) = rho(1:np) * (4.0_dbl*pi * (atom%rg%r(1:np))**2)
       CALL rg_poisson(rho_calc, y0, 0, atom%rg, 'NUMEROV')

       v(1:np) = v(1:np) &                 
            - DBLE(atom%z) / (atom%rg%r(1:np)) &  ! el-core interaction
            + y0(1:np) / (atom%rg%r(1:np))     ! electrostatic el-el repulsion

       ! calculate and add the exchange part
       DO k=1, np
          CALL slater(rho(k), ex, vx)
          v(k) = v(k) + vx
       ENDDO


    CASE ('SIMPLE')
       v(1:np) = v(1:np) - DBLE(atom%z) / (atom%rg%r(1:np))

    END SELECT

  END SUBROUTINE calculate_potential

! **************************************************************************
  
!!****f* cp2k/atom_calculations/calculate_electron_density *
!!
!!   NAME
!!     calculate_electron_density
!!
!!   FUNCTION
!!     Calculates the electron density given a set of wavefunctions (WF) with 
!!     coefficients and weights. The weights are optional and will be applied
!!     componentwise to each WF. The coefficients give the number of electrons
!!     in the specific orbit.
!!
!!   ARGUMENTS
!!     n - REAL(DBL), DIMENSION(:) - the array to store the density in
!!     phi - REAL(DBL), DIMENSION(:,:) - the wavefunctions
!!     c - INTEGER, DIMENSION(:) - coefficients
!!     w - REAL(DBL), DIMENSION(:), OPTIONAL - weights
!!
!!   MODIFICATION HISTORY
!!     15-FEB-2002, TCH, working version
!!
!!   SEE ALSO
!!     
!!
!!****
  SUBROUTINE calculate_electron_density(n, phi, c, w)

    ! Calculates the electron density (n), given a set of wave functions (phi).
    ! c specifies the number of electrons in an orbit and w is a possible weight
    ! (in r-space)

    IMPLICIT NONE

    ! arguments
    REAL(DBL), DIMENSION(1:), INTENT(OUT) :: n          ! the density to be calculated
    REAL(DBL), DIMENSION(1:,1:), INTENT(IN) :: phi      ! the wavefunctions
    INTEGER, DIMENSION(1:), INTENT(IN) :: c             ! number of electrons in the orbit
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
       IF (c(i) /= 0) n = n + DBLE(c(i))*psi(:,i)**2
    ENDDO          

  END SUBROUTINE calculate_electron_density

! **************************************************************************

  SUBROUTINE rg_normalize(f, rg)

    USE mathconstants, ONLY : pi
    USE radial_grids, ONLY : RADIAL_GRID_TYPE
    USE radial_util, ONLY : radial_integrate

    IMPLICIT NONE

    REAL ( DBL ), DIMENSION ( 1: ), INTENT ( INOUT ) :: f
    TYPE ( RADIAL_GRID_TYPE ), INTENT ( IN ) :: rg

    REAL ( DBL ) :: res
    INTEGER :: errno

    CALL radial_integrate ( rg%n, 4.0_dbl*pi * f**2, rg, res, errno )
    f = f / SQRT ( res )

  END SUBROUTINE rg_normalize

! **************************************************************************

  SUBROUTINE count_nodes(nodes, f)

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
    
  END SUBROUTINE count_nodes

! **************************************************************************

  SUBROUTINE initialize_eigenvalues ( atom, parameters, eigenvalues )
    
    IMPLICIT NONE

    TYPE ( AM_ATOM_INFO ), INTENT ( IN ) :: atom
    TYPE ( AM_PARAMETERS_INFO ), INTENT ( IN ) :: parameters
    REAL ( DBL ), DIMENSION ( : ), INTENT ( OUT ) :: eigenvalues

    ! variables
    INTEGER :: n, l, orbit
    INTEGER :: counter

    counter = 0
    DO n=1, 10
       DO l=0, n-1

          orbit = (n*(n-1))/2 + 1 + l

          IF ( am_get_elec_conf( atom, orbit ) == 0 ) THEN 
             eigenvalues( orbit ) = 0.0_dbl
             CYCLE
          END IF

          IF ( am_get_method ( parameters ) == 'SIMPLE' ) THEN
             eigenvalues( orbit ) = -0.5_dbl * ( DBLE ( atom%z ) / DBLE ( n ) )**2
          ELSE
             eigenvalues( orbit) = -0.5_dbl * ( DBLE ( atom%z - counter) / DBLE ( n ) )**2
          END IF
          eigenvalues( orbit ) = eigenvalues( orbit ) * 1.05_dbl
          
          counter = counter + am_get_elec_conf ( atom, orbit )
       END DO
    END DO
  END SUBROUTINE initialize_eigenvalues

! **************************************************************************

  SUBROUTINE calculate_gs_energy ( atom, parameters, rho )

    IMPLICIT NONE

    ! ARGUMENTS
    TYPE ( AM_ATOM_INFO ), INTENT ( INOUT ) :: atom
    TYPE ( AM_PARAMETERS_INFO ), INTENT ( IN ) :: parameters
    REAL ( DBL ), DIMENSION ( : ), INTENT ( IN ) :: rho

    ! VARIABLES
    REAL ( DBL ) :: res
    REAL ( DBL ), DIMENSION(0:atom%rg%n+1) :: y0, rho_calc
    INTEGER :: errno, np, orbit

    ! ------------------------------------------------------------
    np = atom%rg%n
       
    atom%ground_state_energy = SUM( DBLE ( am_get_elec_conf ( atom ) ) * atom%eigenvalues )

    SELECT CASE ( am_get_method ( parameters ) )
       
    CASE DEFAULT
       CALL stop_program(&
            "atom_calculations/calculate_gs_energy","Method not implemented")
       
    CASE ( 'SIMPLE')
       !! nothing to do

    CASE ( 'XA' )
       rho_calc = 0.0_dbl

       rho_calc( 1:np ) = rho( 1:np ) * 4.0_dbl*pi * ( atom%rg%r( 1:np ) )**2
       CALL rg_poisson ( rho_calc, y0, 0, atom%rg, 'NUMEROV' )
       
       rho_calc( 1:np ) = rho * y0( 1:np ) * 4.0_dbl*pi * ( atom%rg%r( 1:np ) )
       CALL radial_integrate ( np, rho_calc( 1:np ), atom%rg, res, errno )
       atom%ground_state_energy = atom%ground_state_energy - 0.5_dbl*res

       rho_calc( 1:np ) = rho**( 4._dbl/3._dbl ) * ( 4.0_dbl*pi*atom%rg%r( 1:np )**2 )
       CALL radial_integrate ( np, rho_calc( 1:np ), atom%rg, res, errno )
       atom%ground_state_energy = atom%ground_state_energy &
            + 0.25_dbl * res * am_get_method_parameters ( parameters, 2 )
                               
    END SELECT

  END SUBROUTINE calculate_gs_energy

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
