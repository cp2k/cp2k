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

    INTEGER :: handle           ! a handle for the timer

    !------------------------------------------------------------------------------
    CALL timeset ( 'ATOMIC_CODE', 'I', ' ', handle )

    IF ( globenv % ionode ) CALL atom_header ( globenv % scr )

    CALL atom_init(this, globenv, .false.)
    CALL atom_do_calc(this)
    CALL timestop ( zero, handle )
    CALL atom_destroy(this)

  END SUBROUTINE atomic_code

  !-----------------------------------------------------------------------------!
  !!****** atom_module/atom_do_calc [1.0] *
  !!
  !!   FUNCTION
  !!     Does the actual computation
  !!
  !!   HISTORY
  !!     TCH, 28-OCT-2001, began work
  !!     TCH, 06-DEC-2001, V0.5: calculates Xa method
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE atom_do_calc(this)

    USE atom_calculations, ONLY : &
         bfsdm, bff, normalize, count_nodes, &
         calc_electron_density, guess_ew, compute_gse
    USE radial_util, ONLY : radial_integrate
    USE termination, ONLY : stop_memory, stop_program
   
    IMPLICIT NONE
    
    ! arguments
    TYPE(atom_info), INTENT(INOUT) :: this ! the atom_info record

    ! variables
    INTEGER :: errorcode      ! to store an errorcode
    INTEGER :: g, i, j, k     ! loop variables
    INTEGER :: ngp            ! holds the number of interior grid points
    
    REAL(DBL), DIMENSION(1:this%rg%n) :: f
    REAL(DBL), DIMENSION(0:this%rg%n+1) :: s, t1, t2
    REAL(DBL), DIMENSION(1:this%rg%n) :: rho  ! 4*pi*r^2 * electron density

    REAL(DBL) :: h2           ! h**2
    REAL(DBL) :: scale        ! scale factor to match wf in the middle
    REAL(DBL) :: etot         ! total computed energy from last iteration
    REAL(DBL) :: ldd          ! logarhytmic derivative difference
    REAL(DBL) :: ldd1, ldd2
    REAL(DBL) :: c1, c2       ! parameters for eigenvalue corrections
         
    INTEGER :: iteration
    INTEGER :: l              ! angular quantum number
    INTEGER :: n              ! main quantum number
    INTEGER :: orbit          ! to keep track which orbital we are computing
    INTEGER :: nodes          ! number of nodes in the wavefunction

    INTEGER :: mp             ! matching point

    ! ------------------------------------------------

    ! initialisations
    c1 = 0.1_dbl
    c2 = 0.01_dbl
    h2 = this%rg%h**2
    ngp = this%rg%n
    s = 0.0_dbl
    t1 = 0.0_dbl; t2 = 0.0_dbl
    mp = 3*ngp/4
    iteration = 0

    ! construct the starting functions
!!!    DO k = 1, ngp
!!!       this%orbitals(k,0) = exp(-this%rg%r(k))*this%rg%r(k)**2
!!!    END DO
!!!    this%orbitals(:,0) = EXP(-1.0_dbl*this%rg%r(1:ngp))*this%rg%r(1:ngp)**2
    this%orbitals(:,0) = 1.0_dbl
    CALL normalize(this%orbitals(:,0), this%rg)
    DO k=1, 55
       this%orbitals(:,k) = this%orbitals(:,0)
    END DO

    ! approximations for the eigenvalues
    CALL guess_ew(this%oew, this)

    etot = 0.0_dbl

    DO 
       
       IF (this%method == 'XA') THEN
          CALL calc_electron_density(rho, this%orbitals, this)
       END IF

       n_loop: DO n = 1, 10
          l_loop: DO l=0, n-1

             IF (this%elconf(n,l) == 0) CYCLE

             orbit = (n*(n-1))/2 + 1 + l

             ! construct the model dependend part of the operator
             CALL bff(f, n, l, this, rho)          

             shooting_loop: DO

                iteration = iteration + 1
                if (iteration > 100) exit
                
                DO k=1, ngp
                   f(k) = -1.0_dbl*(f(k) + 2.0_dbl*this%oew(orbit)/(this%rg%w(k)**4))
                END DO
                   
                t1(0) = 0.0_dbl
                t1(1) = 1.0E-7_dbl
                t2(ngp+1) = 0.0_dbl
                t2(ngp) = 0.5E-7_dbl
!!!                t2(ngp+1) = 0.5-9_dbl
!!!                t2(ngp) = t2(ngp+1) * EXP(-this%rg%h * SQRT(f(ngp)))
                
                ! numerov iteration
                DO k=2, (mp+1)
                   t1(k) = t1(k-1) &
                        * (2.0_dbl+f(k-1)*h2/(1.0_dbl-f(k-1)*h2/12.0_dbl)) &
                        - t1(k-2)
                END DO
                DO k=ngp-1, (mp-1), -1
                   t2(k) = t2(k+1) &
                        * (2.0_dbl+f(k+1)*h2/(1.0_dbl-f(k+1)*h2/12.0_dbl)) &
                        - t2(k+2)
                END DO

                ! get wP
                DO k=1, mp+1
                   t1(k) = t1(k) / (1.0_dbl-h2*f(k)/12.0_dbl)
                END DO
                DO k=ngp, (mp-1), -1
                   t2(k) = t2(k) / (1.0_dbl-h2*f(k)/12.0_dbl)
                END DO
                
                scale = t1(mp) / t2(mp)
                t2(mp-1:ngp) = scale * t2(mp-1:ngp)
                s(1:mp) = t1(1:mp)
                s(mp+1:ngp) = t2(mp+1:ngp)

                this%orbitals(1:ngp,orbit) = s(1:ngp) / this%rg%w(1:ngp)
                CALL normalize(this%orbitals(:,orbit), this%rg)

!!!                DO k=mp-10, mp+10
!!!                   WRITE (6,*) this%rg%r(k), this%orbitals(k,orbit)
!!!                END DO

!!!                ! number of nodes in wavefunction should be n-(l+1)
!!!                CALL count_nodes(nodes, this%orbitals(:,orbit))
!!!                IF (nodes < n-(l+1)) THEN
!!!                   this%oew(orbit) = (1.0_dbl-c1) * this%oew(orbit)
!!!                   CYCLE shooting_loop
!!!                ELSE IF (nodes > n-(l+1)) THEN
!!!                   this%oew(orbit) = (1.0_dbl+c1) * this%oew(orbit)
!!!                   CYCLE shooting_loop
!!!                END IF

!!!                ldd = (t1(mp+1) - t1(mp-1) - t2(mp+1) + t2(mp-1)) &
!!!                     / (2.0_dbl * this%rg%h * s(mp))
!!!                ldd = ((this%orbitals(mp,orbit) - this%orbitals(mp-1,orbit)) &
!!!                     / (this%rg%r(mp) - this%rg%r(mp-1)) &
!!!                     - (this%orbitals(mp+1, orbit) - this%orbitals(mp,orbit)) &
!!!                     / (this%rg%r(mp+1) - this%rg%r(mp))) &
!!!                     / (this%rg%h * this%orbitals(mp,orbit))
                ldd1 = (t1(mp+1) - t1(mp-1))/(2.0_dbl*this%rg%h)
                ldd2 = (t2(mp+1) - t2(mp-1))/(2.0_dbl*this%rg%h)
                ldd = (ldd1 - ldd2) / s(mp)
!!!                WRITE (6,'(1X,4ES13.4)') ldd1, ldd2, ldd1-ldd2, (ldd1-ldd2)/s(mp)
                WRITE (6,'(1X,ES13.4,A,ES13.4)') this%oew(orbit), "  ", ldd
                IF (ABS(ldd) < this%precision) EXIT shooting_loop
                this%oew(orbit) = this%oew(orbit) - c2 * ldd
!!!                EXIT
             END DO shooting_loop

!!!             CALL normalize(this%orbitals(1:ngp,orbit), this%rg)

!!!             DO k=1, ngp
!!!                WRITE (6,'(2(1X,F15.7))') this%rg%r(k), this%orbitals(k,orbit)
!!!             END DO

          END DO l_loop
!!!!!!!!!!!!!!!!!
          EXIT
       END DO n_loop

       CALL compute_gse(this, rho)

       ! these methods dont need the self consistency loop
       IF (this%method == 'SIMPLE' .OR. this%method == 'IPM') EXIT

       ! self consistency condition
       IF (ABS((this%gse - etot)/this%gse) < this%precision) EXIT
       etot = this%gse

    END DO ! self-consistency loop

    WRITE (6,'(1X,A,F10.5)') "GSE: ", this%gse

  END SUBROUTINE atom_do_calc
  
END MODULE atom_module

!******************************************************************************
