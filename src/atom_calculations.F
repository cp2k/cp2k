MODULE atom_calculations

  USE atom_type, ONLY : atom_info, atom_init
  USE kinds, ONLY : dbl
  USE termination, ONLY : stop_program
  USE radial_util, ONLY : radial_integrate
  USE radial_grids, ONLY : radial_grid_type

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: bfsdm, bff, testfsdm, normalize, count_nodes, &
       calc_electron_density, guess_ew, compute_gse

CONTAINS

  !! ===========================================================================
  !! Build the Full Second Derivative Matrix
  !!
  !!
  !! ===========================================================================
  SUBROUTINE bfsdm(a, h, n)
    
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n                    ! size of the matrix
    REAL(DBL), DIMENSION(n,n), INTENT(OUT) :: a ! the matrix in which to store the second derivative
    REAL(DBL), INTENT(IN) :: h                  ! stepsize

    INTEGER :: i                                ! loop variables

    DO i=3, n-2
       a(i,i-2) = -5.0_dbl
       a(i,i-1) = 80.0_dbl
       a(i,  i) = -150.0_dbl
       a(i,i+1) = 80.0_dbl
       a(i,i+2) = -5.0_dbl
    END DO
    a(2,1) = 80.0_dbl; a(2,2) = -150.0_dbl; a(2,3) = 80.0_dbl; a(2,4) = -5.0_dbl
    a(n-1,n) = 80.0_dbl; a(n-1,n-1) = -150.0_dbl; a(n-1,n-2) = 80.0_dbl; a(n-1,n-3) = -5.0_dbl
    a(1,1) = -75._dbl; a(1,2) = -20._dbl; a(1,3) = 70._dbl; a(1,4) = -30._dbl; a(1,5) = 5._dbl
    a(n,n) = a(1,1); a(n,n-1) = a(1,2); a(n,n-2) = a(1,3); a(n,n-3) = a(1,4); a(n,n-4) = a(1,5)
    
    a = a / (h*h*60.0_dbl)

  END SUBROUTINE bfsdm
  

  !! ===========================================================================
  !! Build the F-Function
  !! ===========================================================================
  SUBROUTINE bff(f, n, l, atom, rho)

    USE mathconstants, ONLY : pi
    USE radial_poisson, ONLY : rg_poisson
    USE density_types, ONLY : density_type
    USE potential_types, ONLY : potential_type
    
    IMPLICIT NONE

    ! ARGUMENTS
    REAL(DBL), DIMENSION(:), INTENT(OUT) :: f  ! f to compute
    INTEGER, INTENT(IN) :: n, l                ! main and angular quantum number
    TYPE(atom_info), INTENT(IN) :: atom        ! info on the atom to be computed
    REAL(DBL), DIMENSION(:), INTENT(IN) :: rho ! the electron density

    ! VARIABLES
    INTEGER :: errno
    INTEGER :: i,j
    INTEGER :: k, ngp
    REAL(DBL), DIMENSION(atom%rg%n) :: v       ! potential
    REAL(DBL) :: r
    REAL(DBL) :: w                             ! weights of radial transformation

    ! needed for XA
    REAL(DBL), DIMENSION(0:atom%rg%n+1) :: y0  ! r times coulomb repulsion potential
    REAL(DBL), DIMENSION(0:atom%rg%n+1) :: rho_calc ! electron density (for calc)

    ! needed for IPM
    INTEGER :: totnel                          ! total number of electrons
    REAL(DBL) :: omega, etha, xi


    ! ----------------------------------------------------------
    ngp = atom%rg%n
    totnel = SUM(atom%elconf)

    ! zentrifugal term
    v(1:ngp) = DBLE(l)*(DBLE(l)+1.0_dbl)/(2.0_dbl*atom%rg%r(1:ngp)**2)

    SELECT CASE (atom%method)

    CASE DEFAULT

       CALL stop_program(&
            "atom_calculations/bff","Method not implemented")

    CASE ('XA')

       ! solve for the electrostatic potential of the electron density
       rho_calc = 0.0_dbl
       rho_calc(1:ngp) = rho(1:ngp)
       CALL rg_poisson(rho_calc, y0, 0, atom%rg, 'DIFF6P')

       DO k=1, ngp
          r = atom%rg%r(k)
          v(k) = v(k) - DBLE(atom%z)/r     ! electron-core interaction
          v(k) = v(k) + y0(k)/r            ! electron-electron interaction
          v(k) = v(k) - atom%mpara(2) * &  ! exchange potential
               (rho(k)/(4.0_dbl*pi*r**2))**(1.0_dbl/3.0_dbl)
       END DO
       
    CASE ('HF')
       
    CASE ('IPM')
       CALL stop_program(&
            "atom_calculations/bff","Method not implemented")
       DO k=1, ngp
          r = atom%rg%r(k)

          ! N = 2 hardcoded for now
          ! read from a table later
          etha = 1.770_dbl + 1.1402_dbl * DBLE(atom%z - totnel)
          xi = 2.625_dbl + 1.2996_dbl * DBLE(atom%z -totnel)

          omega = 1.0_dbl/((etha/xi)*(EXP(xi*r)-1.0_dbl) + 1.0_dbl)
          v(k) = v(k) + 2.0_dbl*(DBLE(totnel - 1)*(1.0_dbl - omega) - DBLE(atom%z))/r
       END DO

    CASE ('SIMPLE')
       DO k=1, ngp
          r = atom%rg%r(k)
          v(k) = v(k) - DBLE(atom%z)/r
       END DO

    END SELECT

    ! transform the function for s-space
    DO k=1, ngp
       w = atom%rg%w(k)
       f(k) = -(atom%rg%wdd(k))/w - (2.0_dbl/(w**4))*v(k)
    END DO
    
  END SUBROUTINE bff

  
  !! ===========================================================================
  !! calculate the electron density
  !! actually this calculates 4*pi*r^2 * the density
  !! ===========================================================================
  SUBROUTINE calc_electron_density(rho, wf, atom)
    IMPLICIT NONE
    REAL(DBL), DIMENSION(:), INTENT(OUT) :: rho
    REAL(DBL), DIMENSION(:,0:), INTENT(IN) :: wf
    TYPE(atom_info), INTENT(IN) :: atom

    INTEGER :: n, l ! main and angular quantum numbers
    INTEGER :: k    ! loop index
    INTEGER :: orbit
    
    rho = 0.0_dbl
    DO k=1, atom%rg%n
       DO n=1, 10
          DO l=0, n-1
             IF (atom%elconf(n,l) /= 0) THEN
                orbit = (n*(n-1))/2+1+l
                rho(k) = rho(k) + atom%elconf(n,l) * &
                     (wf(k, orbit))**2
             END IF
          END DO
       END DO
    END DO

  END SUBROUTINE calc_electron_density


  !! ===========================================================================
  !! Test the Second Derivative Matrix (for debugging purposes only)
  !!
  !!
  !! ===========================================================================
  SUBROUTINE testfsdm(derivative, rg)
    IMPLICIT NONE
    TYPE(radial_grid_type), INTENT(IN) :: rg
    REAL(DBL), DIMENSION(:,:), INTENT(IN) :: derivative
    
    REAL(DBL), DIMENSION(rg%n,rg%n) :: a
    REAL(DBL), DIMENSION(rg%n) :: tf, res, cres
    REAL(DBL) :: r, w
    INTEGER :: k, i, j

    a = derivative
    DO k=1, rg%n
       r = rg%r(k)
       w = rg%w(k)
       tf(k) = w*r*EXP(-r)
       cres(k) = (-2.0*exp(-r) + r*exp(-r))/(w**3)
       a(k,k) = a(k,k) - rg%wdd(k)/w
    END DO
    
    DO i=1, rg%n
       r = rg%r(i)
       DO j=1, rg%n
          res(i) = res(i) + a(i,j)*tf(j)
       END DO
       WRITE (6,'(1X,I4,3F14.7,ES12.3)') i, r, res(i), cres(i), abs((res(i) - cres(i))/cres(i))
    END DO
    
  END SUBROUTINE testfsdm


  !! ===========================================================================
  !! Normalize a function f
  !!
  !! f := f/sqrt(int(|f|^2))
  !! ===========================================================================
  SUBROUTINE normalize(f, rg)
    IMPLICIT NONE
    REAL(DBL), DIMENSION(:), INTENT(INOUT) :: f
    TYPE(radial_grid_type), INTENT(IN) :: rg
    REAL(DBL) :: res
    INTEGER :: info
    CALL radial_integrate(rg%n, f*f, rg, res, info)
    f = f/SQRT(res)
  END SUBROUTINE normalize


  !!===============================================================
  !! counts the nuber of nodes in a function f
  !! f must start positively                  
  !!
  !!===============================================================
  SUBROUTINE count_nodes(nodes, f)

    IMPLICIT NONE

    INTEGER, INTENT(OUT) :: nodes
    REAL(DBL), DIMENSION(:), INTENT(IN) :: f
    
    INTEGER :: n, i
    INTEGER :: signum

    n = SIZE(f)
    IF ((n < 4) .OR. (f(2) < 0.0_dbl)) THEN ! not enough or bad function data
       nodes = -1
       RETURN
    END IF
    
    nodes = 0
    signum = 1
    DO i=3, n-1
       IF (INT(f(i)/ABS(f(i))) /= signum) THEN
          signum = INT(f(i)/ABS(f(i)))
          IF (ABS(f(i+1)) > ABS(f(i))) THEN
             nodes = nodes + 1
          END IF
       END IF
    END DO
        
  END SUBROUTINE count_nodes


  !!===============================================================
  !!
  !!
  !!
  !!===============================================================
  SUBROUTINE guess_ew(ew, atom)
    
    IMPLICIT NONE
    REAL(DBL), DIMENSION(:), INTENT(OUT) :: ew ! the vector with the eigenvalues
    TYPE(atom_info), INTENT(IN) :: atom        ! 

    INTEGER :: n, l                            ! main and angular quantum numbers
    INTEGER :: counter, orbit

    counter = 0
    DO n=1, 10 ! arbitrary limit
       DO l=0, n-1

          IF (atom%elconf(n,l) == 0) CYCLE

          orbit = (n*(n-1))/2 + 1 + l

          IF (atom%method == 'SIMPLE') THEN
             ew(orbit) = -0.5_dbl * (dble(atom%z)/DBLE(n))**2
          ELSE
             ew(orbit) = -0.5_dbl * (DBLE(atom%z - counter)/DBLE(n))**2
          END IF
          ew(orbit) = ew(orbit) * 1.05_dbl ! shift the guesses by 5%

          ! calculate the running total of electrons
          counter = counter + atom%elconf(n,l)

       END DO
    END DO
  END SUBROUTINE guess_ew


  !!===============================================================
  !!
  !!
  !!
  !!===============================================================
  SUBROUTINE compute_gse(atom, rho)

    USE mathconstants, ONLY : pi
    USE radial_poisson, ONLY : rg_poisson

    IMPLICIT NONE

    ! ARGUMENTS
    TYPE(atom_info), INTENT(INOUT) :: atom
    REAL(DBL), DIMENSION(:), INTENT(IN) :: rho

    ! VARIABLES
    REAL(DBL) :: gse
    REAL(DBL) :: res, alpha, cx
    REAL(DBL), DIMENSION(0:atom%rg%n+1) :: y0, rho_calc
    INTEGER :: k, n, l, errno, ngp

    ! ------------------------------------------------------------
    ngp = atom%rg%n

    gse = 0.0_dbl
    DO n=1, 10
       DO l=0, n-1
          IF (atom%elconf(n,l) /= 0) &
               gse = gse + atom%elconf(n,l) * atom%oew((n*(n-1))/2 + l + 1)
       END DO
    END DO

    SELECT CASE (atom%method)
       
    CASE DEFAULT
       CALL stop_program(&
            "atom_calculations/compute_gse","Method not implemented")
       
    CASE ('SIMPLE')
       !! nothing to do

    CASE ('XA')
       rho_calc = 0.0_dbl
       rho_calc(1:ngp) = rho(1:ngp)
       CALL rg_poisson(rho_calc, y0, 0, atom%rg, 'NUMEROV')
       CALL radial_integrate(ngp, &
            rho*y0/(atom%rg%r(1:ngp)), &
            atom%rg, res, errno)
       gse = gse - 0.5_dbl*res

       rho_calc = 0.0_dbl
       rho_calc(1:ngp) = (rho(1:ngp) / &
            (4.0_dbl*pi*atom%rg%r(1:ngp)**2))**(4.0_dbl/3.0_dbl)
       rho_calc(1:ngp) = rho_calc(1:ngp)*(4.0_dbl*pi*atom%rg%r(1:ngp)**2)
       CALL radial_integrate(ngp, rho_calc(1:ngp), atom%rg, res, errno)
       gse = gse + 0.25_dbl * atom%mpara(2) * res

    END SELECT

    atom%gse = gse
  END SUBROUTINE compute_gse

END MODULE atom_calculations
