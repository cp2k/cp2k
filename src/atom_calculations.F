MODULE atom_calculations

  USE atom_type, ONLY : atom_info, atom_init
  USE kinds, ONLY : dbl
  USE termination, ONLY : stop_program
  USE radial_util, ONLY : radial_integrate
  USE radial_grids, ONLY : radial_grid_type
  
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: bfsdm, bff, testfsdm, normalize, reportf, count_nodes, &
       calc_electron_density, guess_ew
  
  CONTAINS

  !! ===========================================================================
  !! Build the Full Second Derivative Matrix
  !! ===========================================================================
  SUBROUTINE bfsdm(a,h,n)

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
!!!    a(1,1) = 720.0_dbl; a(1,2) = -1440.0_dbl; a(1,3) = 720.0_dbl
!!!    a(n,n) = 720.0_dbl; a(n,n-1) = -1440.0_dbl; a(n,n-2) = 720.0_dbl
    a(1,1) = -75.0_dbl; a(1,2) = -20.0_dbl; a(1,3) = 70.0_dbl; a(1,4) = -30.0_dbl; a(1,5) = 5.0_dbl
    a(n,n) = a(1,1); a(n,n-1) = a(1,2); a(n,n-2) = a(1,3); a(n,n-3) = a(1,4); a(n,n-4) = a(1,5)
    
    a = a / (h*h*60.0_dbl)

  END SUBROUTINE bfsdm
  

  !! ===========================================================================
  !! Build the F-Function
  !! ===========================================================================
  SUBROUTINE bff(f, n, l, atom, rho)

    USE mathconstants, ONLY : pi
    
    IMPLICIT NONE
    
    REAL(DBL), DIMENSION(:), INTENT(OUT) :: f  ! f to compute
    INTEGER, INTENT(IN) :: n, l                ! main and angular quantum number
    TYPE(atom_info), INTENT(IN) :: atom        ! info on the atom to be computed
    REAL(DBL), DIMENSION(:), INTENT(IN) :: rho ! the electron density


    ! generally needed
    INTEGER :: i,j
    INTEGER :: k, ngp
    REAL(DBL), DIMENSION(atom%rg%n) :: v
    REAL(DBL) :: r, w

    ! needed for DFT
    REAL(DBL), DIMENSION(0:atom%rg%n+1) :: y0   ! r times coulomb repulsion potential
    REAL(DBL), DIMENSION(0:atom%rg%n+1) :: rho_calc ! for the calculations in radial_poission
                                                ! we need a rho with dimension ngp+2

    ! needed for IPM
    INTEGER :: totnel                           ! total number of electrons
    REAL(DBL) :: omega, etha, xi

    ngp = atom%rg%n

    DO k=1, ngp
       r = atom%rg%r(k)
       v(k) = DBLE(l)*(DBLE(l)+1.0_dbl)/(2.0_dbl*r*r)          ! zentrifugal term
    END DO

    SELECT CASE (atom%method)
    CASE DEFAULT

    CASE ('DFT')
       
       rho_calc = 0.0_dbl
       rho_calc(1:ngp) = rho(1:ngp)

       CALL radial_poisson(rho_calc, y0, 0, atom%rg, 'NUMEROV')

!!       do k=1, ngp, ngp/20
!!          r = atom%rg%r(k)
!!          write (6,'(1x,3F14.8)') r, y0(k)/r, rho_calc(k)
!!       end do
!!       write (6,*) 
!!       write (6,*)

       DO k=1, ngp
          r = atom%rg%r(k)
          
          v(k) = v(k) - DBLE(atom%z)/r ! electron-ion interaction
          v(k) = v(k) + y0(k)/r        ! electron-electron interaction
          v(k) = v(k)
       END DO
       
    CASE ('HF')
       
    CASE ('IPM')
       DO k=1, ngp
          r = atom%rg%r(k)
          totnel = SUM(atom%nel)

          ! N = 2 hardcoded for now
          ! read from a table later
          etha = 1.770_dbl + 1.1402_dbl * DBLE(atom%z - totnel)
          xi = 2.625_dbl + 1.2996_dbl * DBLE(atom%z -totnel)

          omega = 1.0_dbl/((etha/xi)*(EXP(xi*r)-1.0_dbl) + 1.0_dbl)
          v(k) = v(k) + 2.0_dbl*(DBLE(totnel - 1)*(1.0_dbl - omega) - DBLE(atom%z))/r
       END DO

    CASE ('SIMPLE')
       totnel = 0
       DO i=1, n-1
          totnel = totnel + SUM(atom%elconf(i,:))
       END DO
       IF (l>0) totnel = totnel + SUM(atom%elconf(n,0:l-1))
       !! testing
       totnel = 0
       DO k=1, ngp
          r = atom%rg%r(k)
          v(k) = v(k) - DBLE(atom%z - totnel)/r
       END DO

    END SELECT

    ! transform the function for s-space
    DO k=1, ngp
       w = atom%rg%w(k)
       f(k) = -(atom%rg%wdd(k))/w - (2.0_dbl/(w**4))*v(k)    ! build the f-function (s-space)
    END DO
    
  END SUBROUTINE bff

  
  !! ===========================================================================
  !! calculate the electron density
  !! ===========================================================================
  SUBROUTINE calc_electron_density(rho, wf, atom)
    IMPLICIT NONE
    REAL(DBL), DIMENSION(:), INTENT(OUT) :: rho
    REAL(DBL), DIMENSION(:,:), INTENT(IN) :: wf
    TYPE(atom_info), INTENT(IN) :: atom

    INTEGER :: n, l ! main and angular quantum numbers
    INTEGER :: k    ! loop index
    
    rho = 0.0_dbl
    DO k=1, atom%rg%n
       DO n=1, 10
          DO l=0, n-1
             IF (atom%elconf(n,l) /= 0) &
                  rho(k) = rho(k) + DBLE(2*(2*l+1)) * wf(k, n*(n-1)/2+1+l)**2
          END DO
       END DO
    END DO

  END SUBROUTINE calc_electron_density


  !! ===========================================================================
  !! Test the Second Derivative Matrix (for debugging purposes only)
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
  !! Diagonalize the matrix
  !! ===========================================================================
  SUBROUTINE diagonalize(n, a, eigenvalues, eigenvectors)
    IMPLICIT NONE
    ! ARGUMENTS
    INTEGER, INTENT(IN) :: n                               ! size of involved matrizes
    REAL(DBL), DIMENSION(n,n), INTENT(INOUT) :: a          ! matrix to diagonalize
    REAL(DBL), DIMENSION(n,2), INTENT(OUT) :: eigenvalues  ! complex
    REAL(DBL), DIMENSION(n,n), INTENT(OUT) :: eigenvectors ! the right eigenvectors
    ! VARIABLES
    INTEGER :: info                                        !
    INTEGER, PARAMETER :: lwork = 10                       ! variables needed for lapack call
    REAL(DBL), DIMENSION(:), ALLOCATABLE :: VR             !
    REAL(DBL), DIMENSION(lwork*n) :: work                  !
    ! EXECUTION SECTION
    CALL dgeev('N', 'V', n, &
         a, n, &
         eigenvalues(:,1), eigenvalues(:,2), &
         vr, 1, &
         eigenvectors, n, &
         work, lwork*n, &
         info)
  END SUBROUTINE diagonalize
  
  
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


  !! ===========================================================================
  !! Report a function
  !! ===========================================================================
  SUBROUTINE reportf(f, p, rg)
    IMPLICIT NONE

    REAL(DBL), DIMENSION(:), INTENT(IN) :: f, p
    TYPE(radial_grid_type), INTENT(IN) :: rg

    INTEGER :: k, n
    REAL(DBL) :: r

    n = SIZE(f)

    WRITE (6,'(1X,A3,2X,A12,2X,A10,2X,A10,2X,A10)') "k", "r(k)", "f(k)", "p(k)", "|f-p|/p"
    WRITE (6,*)
    DO k=1, n, n/20
       r = rg%r(k)
       WRITE (6,'(1X,I3,(A),F12.8,(a),F10.8,(a),F10.8,(a),F10.8)') k, "  ", r, "  ", &
            f(k), "  ", p(k), "  ", ABS(f(k)-p(k))/p(k)
    END DO
    WRITE (6,*)
  END SUBROUTINE reportf

  
  ! counts the nuber of nodes in a function f
  ! f must start positively
  SUBROUTINE count_nodes(f, nodes)

    IMPLICIT NONE

    REAL(DBL), DIMENSION(:), INTENT(IN) :: f
    INTEGER, INTENT(OUT) :: nodes
    
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
       IF (int(f(i)/abs(f(i))) /= signum) THEN
          nodes = nodes + 1
          signum = int(f(i)/abs(f(i)))
       END IF
    END DO
          
  END SUBROUTINE count_nodes

  SUBROUTINE guess_ew(ew, atom)
    
    IMPLICIT NONE
    REAL(DBL), DIMENSION(:), INTENT(OUT) :: ew ! the vector with the eigenvalues
    TYPE(atom_info), INTENT(IN) :: atom        ! 

    INTEGER :: n, l                            ! main and angular quantum numbers
    INTEGER :: counter, index

    counter = 0
    DO n=1, 10 ! arbitrary limit
       DO l=0, n-1

          IF (atom%elconf(n,l) == 0) CYCLE

          index = n*(n-1)/2 + 1 + l

          IF (atom%method == 'SIMPLE') THEN
             ew(index) = -0.5_dbl * (dble(atom%z)/dble(n))**2
          ELSE
             ew(index) = -0.5_dbl * (dble(atom%z - counter)/dble(n))**2
          END IF
          ew(index) = ew(index) * (1.01_dbl + 0.0_dbl*atom%precision)

          ! calculate the running total of electrons
          counter = counter + atom%elconf(n,l)

       END DO
    END DO
  END SUBROUTINE guess_ew

END MODULE atom_calculations
