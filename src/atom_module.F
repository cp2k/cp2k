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

  USE global_types, ONLY : global_environment_type
  USE header, ONLY : atom_header
  USE kinds, ONLY : dbl
  USE mathconstants, ONLY : zero, pi
  USE timings, ONLY : timeset, timestop
  USE termination

  USE atom_calculations
  USE atom_type, ONLY : atom_info, atom_init

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: atomic_code, atom_do_calc

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

  END SUBROUTINE atomic_code

  !-----------------------------------------------------------------------------!
  !!****** atom_module/atom_do_calc [1.0] *
  !!
  !!   FUNCTION
  !!     Does the actual computation using the independend particle model
  !!     (no interactions)
  !!
  !!   HISTORY
  !!     TCH, 28-OCT-2001, began work
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE atom_do_calc(this)

    USE radial_util, ONLY : radial_integrate
    
    IMPLICIT NONE
    
    ! arguments
    TYPE(atom_info), INTENT(IN) :: this ! the atom_info record

    ! variables
    INTEGER :: errorcode      ! to store an errorcode
    INTEGER :: g, i, j, k     ! loop variables
    INTEGER :: n              ! holds the number of interior grid points
    
    REAL(DBL), DIMENSION(this%rg%n, this%rg%n) :: sdm          ! second derivative matrix 
    REAL(DBL), DIMENSION(this%rg%n, this%rg%n) :: hm           ! h operator matrix
    REAL(DBL), DIMENSION(this%rg%n, this%rg%n) :: hcalcm       ! h operator matrix for lapack calculations
    REAL(DBL), DIMENSION(:,:), ALLOCATABLE :: wf               ! the wavefunctions

    REAL(DBL), DIMENSION(this%rg%n) :: f, ipiv, y
    REAL(DBL), DIMENSION(:), ALLOCATABLE :: ew                 ! vector of eigenwerte
 
    REAL(DBL) :: sigma        ! guess for the eigenwert
    REAL(DBL) :: old_ew       ! store the old eigenwert for sc criterium
    REAL(DBL) :: w            ! cache for the weights of the radial transformation
    REAL(DBL) :: r ! DEBUGGING ONLY
         
    INTEGER :: l              ! angular quantum number
    INTEGER :: index, counter ! to keep track which wf we are computing

    CHARACTER(50) :: fmt

    ! execution section
    
    n = this%rg%n

    ALLOCATE(ew(SUM(this%nel)), STAT=errorcode)
    IF (errorcode /= 0) CALL stop_memory("atom_module/atom_ipm", "ew", SUM(this%nel))
    ALLOCATE(wf(1:n, 0:SUM(this%nel)), STAT=errorcode)
    IF (errorcode /= 0) CALL stop_memory("atom_module/atom_ipm", "wf", n*(SUM(this%nel)+1))

    sdm = 0.0_dbl
    CALL bfsdm(sdm, this%rg%h, n)

    ! construct the starting functions
    wf = 1.0_dbl
    CALL normalize(wf(:,1), this%rg)
    DO k=2, SUM(this%nel)
       wf(:,k) = wf(:,1)
    END DO

    counter = 0
    DO l=0, 9
       
       IF (this%nel(l) > 0) THEN 
          
          WRITE (6,'(1X,(A),I1)') "l = ", l

          DO index=1+counter, this%nel(l)+counter

             sigma = dble(this%z)**2*0.5_dbl/(DBLE(index-counter+l)**2) + 1000.0_dbl*this%precision
             old_ew = sigma

             DO                        ! self-consistent loop
                CALL bff(f, l, this)

                hm = sdm
                DO j=1, n
                   hm(j,j) = hm(j,j) + f(j)
                   w = this%rg%w(j)
                   DO i=1, n
                      hm(i,j) = hm(i,j) * (w**4) / 2.0_dbl
                   END DO
                END DO

                ! build the operator for the eigenvector calculation
                hcalcm = hm
                DO k=1, n
                   hcalcm(k,k) = hcalcm(k,k) - sigma
                END DO

                CALL dgesv(n, 1, hcalcm, n, ipiv, wf(:,index), n, errorcode) ! solve the linear equations

                ! normalize and take the positive solution
                CALL normalize(wf(:,index), this%rg)
                ! choose the solution that starts positively
                IF (wf(1,index) < 0.0_dbl) wf(:,index) = -1.0_dbl * wf(:,index) 

                ! compute the approximation to the eigenvalue
                y = 0.0_dbl
                DO i=1, n
                   DO j=1, n
                      y(i)=y(i) + hm(i,j)*wf(j,index)
                   END DO
                END DO
                CALL radial_integrate(this%rg%n, wf(:,index)*y, this%rg, ew(index), errorcode)   
!!                write (6,"1X,F12.9") -ew(index)
                IF (ABS((ew(index) - old_ew)/ew(index)) < this%precision) EXIT
                old_ew = ew(index)

             END DO
             WRITE (6, '(1X,F12.9)') -ew(index)

             DO k=1, n
                r = this%rg%r(k)
                w = this%rg%w(k)
                wf(k,index) = wf(k,index) * (w**3) / 2.0_dbl

                wf(k,0) = 2.0 * (real(this%z)/2.0)**1.5 * (1.0 - real(this%z)*r/2.0) * exp(-0.5*dble(this%z)*r) * r
             END DO
             CALL normalize(wf(:,index), this%rg)
             IF (wf(1,index) < 0.0_dbl) wf(:,index) = -1.0_dbl * wf(:,index) 
             

             CALL reportf(wf(:,index), wf(:,0), this%rg)

          END DO ! index loop
          counter = counter + this%nel(l)

          
          WRITE (6,*)

       END IF
    END DO ! l loop

    DEALLOCATE(ew, wf)

  END SUBROUTINE atom_do_calc

END MODULE atom_module

!******************************************************************************
