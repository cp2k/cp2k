!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/atom_module [1.0] *
!!
!!   FUNCTION
!!
!!   MODIFICATION HISTORY
!!     JGH, 17-NOV-2000, began work
!!
!!   SOURCE
!******************************************************************************


MODULE atom_module

  USE atom_type, ONLY : atom_info, atom_init, atom_destroy
  USE global_types, ONLY : global_environment_type
  USE header, ONLY : atom_header
  USE kinds, ONLY : dbl
  USE mathconstants, ONLY : zero
  USE timings, ONLY : timeset, timestop


  IMPLICIT NONE

  PRIVATE
  PUBLIC :: atomic_code

  !!*****
  !-----------------------------------------------------------------------------!

CONTAINS

  !-----------------------------------------------------------------------------!
  !!****** atom_module/atomic_code [1.0] *
  !!
  !!   FUNCTION
  !!
  !!
  !!   MODIFICATION HISTORY
  !!     JGH, 17-NOV-2000, began work
  !!
  !!   SOURCE
  !******************************************************************************

  SUBROUTINE atomic_code(globenv)

    IMPLICIT NONE

    !Arguments
    TYPE(global_environment_type), INTENT(INOUT) :: globenv

    !Locals
    TYPE(atom_info) :: this     ! the atom itself

    CHARACTER(20) :: tempstring

    INTEGER :: handle           ! a handle for the timer
    INTEGER :: i                ! loop variable

    !------------------------------------------------------------------------------
    CALL timeset ( 'ATOMIC_CODE', 'I', ' ', handle )

    IF ( globenv % ionode ) CALL atom_header ( globenv % scr )

    ! read input file and setup calculations
    CALL atom_init(this, globenv, .true.)

    ! do calculations
    CALL atom_do_calc(this)

    CALL timestop ( zero, handle )

    CALL atom_destroy(this)

  END SUBROUTINE atomic_code

  !-----------------------------------------------------------------------------!
  !!****** atom_module/atom_do_calc [1.0] *
  !!
  !!   FUNCTION
  !!     Does the actual computation
  !!     (no interactions)
  !!
  !!   HISTORY
  !!     TCH, 28-OCT-2001, began work
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE atom_do_calc(this)

    USE atom_calculations, ONLY : bfsdm, testfsdm, bff, normalize, reportf, count_nodes, &
         calc_electron_density, guess_ew
    USE radial_util, ONLY : radial_integrate
    USE termination, ONLY : stop_memory
   
    IMPLICIT NONE
    
    ! arguments
    TYPE(atom_info), INTENT(INOUT) :: this ! the atom_info record

    ! variables
    INTEGER :: errorcode      ! to store an errorcode
    INTEGER :: g, i, j, k     ! loop variables
    INTEGER :: ngp              ! holds the number of interior grid points
    
    REAL(DBL), DIMENSION(this%rg%n, this%rg%n) :: sdm          ! second derivative matrix 
    REAL(DBL), DIMENSION(this%rg%n, this%rg%n) :: hm           ! h operator matrix
    REAL(DBL), DIMENSION(this%rg%n, this%rg%n) :: hcalcm       ! h operator matrix for lapack calculations
    REAL(DBL), DIMENSION(this%rg%n) :: f
    REAL(DBL), DIMENSION(this%rg%n) :: rho                     ! electron density
    REAL(DBL), DIMENSION(55) :: sigma                          ! guess for the eigenvalues

    REAL(DBL), DIMENSION(0:this%rg%n+1) :: testwf

    REAL(DBL), DIMENSION(this%rg%n) :: ipiv, y                 ! helpers for the calculations
 
    REAL(DBL) :: corr         ! the percentage of sigma to be corrected
    REAL(DBL) :: corr_tresh   ! the threshold for the corrections to be applied
    
    REAL(DBL) :: old_ew       ! store the old eigenwert for sc criterium
    REAL(DBL) :: w            ! cache for the weights of the radial transformation
    REAL(DBL) :: fac          ! a shorthand for a factor
    REAL(DBL) :: r ! DEBUGGING ONLY

         
    INTEGER :: l              ! angular quantum number
    INTEGER :: n              ! main quantum number
    INTEGER :: orbitnum, counter ! to keep track which orbital we are computing
    INTEGER :: nodes

    CHARACTER(50) :: fmt

    ! parameters
    corr = 25.0_dbl
    corr_tresh = 100.0_dbl

    ! execution section
    ngp = this%rg%n

    sdm = 0.0_dbl
    CALL bfsdm(sdm, this%rg%h, ngp)

    ! construct the starting functions
    DO k=1, 55
       CALL construct_initial_wf(this%orbitals(:,k), this)
    END DO

    ! approximations for the eigenvalues
    CALL guess_ew(sigma, this)

    DO 

       IF (this%method == 'DFT') CALL calc_electron_density(rho, this%orbitals, this)
      
       DO n = 1, 10
          DO l=0, n-1

             IF (this%elconf(n,l) == 0) CYCLE

             orbitnum = n*(n-1)/2 + 1 + l
             WRITE (6,'(1X,3(A,I2))') "n = ", n, "  l = ", l, "  orbitnum = ", orbitnum

             old_ew = sigma(orbitnum)

             ! construct the model dependend part of the operator
             CALL bff(f, n, l, this, rho)          

             DO ! until convergance

                hm = sdm
                DO j=1, ngp
                   hm(j,j) = hm(j,j) + f(j)
                   w = this%rg%w(j)
                   DO i=1, ngp
                      hm(i,j) = hm(i,j) * w
                   END DO
                END DO

                ! build the operator for the eigenvector calculation
                hcalcm = hm
                DO k=1, ngp
                   fac = 2.0_dbl / (this%rg%w(k)**3)
                   hcalcm(k,k) = hcalcm(k,k) + sigma(orbitnum) * fac
                   this%orbitals(k,orbitnum) = this%orbitals(k,orbitnum) * fac
                END DO

                ! solve the linear equations
                CALL dgesv(ngp, 1, hcalcm, ngp, ipiv, this%orbitals(:,orbitnum), ngp, errorcode) 

                ! normalize and take the positive solution
                CALL normalize(this%orbitals(:,orbitnum), this%rg)
                ! choose the solution that starts positively
                IF (this%orbitals(1,orbitnum) < 0.0_dbl) this%orbitals(:,orbitnum) = -1.0_dbl * this%orbitals(:,orbitnum) 

                ! now check the number of nodes and adjust the approximation
                ! to the eigenvalue accordingly
                ! there should be n-(l+1) nodes in the solution
                CALL count_nodes(this%orbitals(:,orbitnum), nodes)
                IF (nodes < n-(l+1)) THEN
!!                   WRITE (6,'(1X,A,I2,A,F14.8,A,F14.8)') &
!!                        "Too FEW nodes  : ", nodes, " ", sigma(orbitnum), " -> " , &
!!                        (1.0_dbl-corr/100.0_dbl) * sigma(orbitnum)
                   sigma(orbitnum) = (1.0_dbl - corr/100.0_dbl) * sigma(orbitnum)
                   CALL construct_initial_wf(this%orbitals(:,orbitnum), this)
                   CYCLE
                ELSE IF (nodes > n-(l+1)) THEN
!!                   WRITE (6,'(1X,A,I2,A,F14.8,A,F14.8)') &
!!                        "Too MANY nodes : ", nodes, " ", sigma(orbitnum), " -> " , &
!!                        (1.0_dbl+corr/100.0_dbl) * sigma(orbitnum)
                   sigma(orbitnum) = (1.0_dbl + corr/100.0_dbl) * sigma(orbitnum)
                   CALL construct_initial_wf(this%orbitals(:,orbitnum), this)
                   CYCLE
                END IF

                ! compute the new approximation to the eigenvalue
                y = 0.0_dbl
                DO i=1, ngp
                   DO j=1, ngp
                      y(i)=y(i) + hm(i,j)*this%orbitals(j,orbitnum)
                   END DO
                   y(i) = -1.0_dbl * y(i) * (this%rg%w(i)**3) / 2.0_dbl
                END DO
                y = y * this%orbitals(:,orbitnum)
                do k=1, ngp
                   write (6,*) k, y(k)
                end do
                CALL radial_integrate(ngp, this%orbitals(:,orbitnum), this%rg, this%oew(orbitnum), errorcode)   
                IF (ABS((this%oew(orbitnum) - old_ew)/this%oew(orbitnum)) < this%precision) EXIT
                old_ew = this%oew(orbitnum)

                ! Correct sigma to improve convergence
                IF ((sigma(orbitnum) - this%oew(orbitnum))/sigma(orbitnum) > 1.0_dbl) THEN
!!                   WRITE (6,'(1X,A,F14.8,A,F14.8)') &
!!                        "Sigma to BIG  : ", sigma(orbitnum), " -> " , &
!!                        (1.0_dbl - corr/400.0_dbl) * sigma(orbitnum)
                   sigma(orbitnum) = (1.0_dbl - corr/400.0_dbl) * sigma(orbitnum)
                ELSE IF ((sigma(orbitnum) - this%oew(orbitnum))/sigma(orbitnum) < -1.0_dbl) THEN
!!                   WRITE (6,'(1X,A,F14.8,A,F14.8)') &
!!                        "Sigma to SMALL: ", sigma(orbitnum), " -> " , &
!!                        (1.0_dbl + corr/400.0_dbl) * sigma(orbitnum)
                   sigma(orbitnum) = (1.0_dbl + corr/400.0_dbl) * sigma(orbitnum)
                END IF

             END DO

             WRITE (6, '(1X,F16.9)') this%oew(orbitnum)
             WRITE (6,*)

          END DO ! l loop
       END DO ! n loop

       ! these methods dont need the self consistency loop
       IF (this%method == 'SIMPLE' .OR. this%method == 'IPM') EXIT

    END DO ! self-consistency loop

  END SUBROUTINE atom_do_calc


  !! ===============================================================
  !! Fill the initial wavefunctions with a reasonable function.
  !! The function must be 0 at the inner and outer borders.
  !! ===============================================================
  SUBROUTINE construct_initial_wf(f, this)

    USE atom_calculations, ONLY : normalize
    
    IMPLICIT NONE

    REAL(DBL), DIMENSION(:), INTENT(OUT) :: f
    TYPE(atom_info), INTENT(IN) :: this

    INTEGER :: k

    f = exp(-this%rg%r)*this%rg%r**2
    CALL normalize(f, this%rg)
    
  END SUBROUTINE construct_initial_wf

END MODULE atom_module

!******************************************************************************
