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
!!     TCH, 21-DEC-2001, version 0.1 (Xa works)
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

    CALL atom_init(this, globenv,.true.)
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
         calculateElectronDensity, &
         getGroundStateEnergy, calculatePotential, &
         countNodes, initialEigenvalues, rgNormalize, signum
    USE radial_util, ONLY : radial_integrate
    USE termination, ONLY : stop_memory, stop_program
    USE mathconstants, ONLY : pi

    IMPLICIT NONE
    
    ! arguments
    TYPE(atom_info), INTENT(INOUT) :: this ! the atom_info record

    ! variables
    INTEGER :: np            ! holds the number of interior grid points
    
    REAL(DBL), DIMENSION(1:this%rg%n) :: f, fprime
    REAL(DBL), DIMENSION(0:this%rg%n+1) :: t1, t2
    REAL(DBL), DIMENSION(1:this%rg%n) :: rho, new_rho  ! 4*pi*r^2 * electron density

    REAL(DBL) :: a            ! parameter to squish the ldd to a useable number
    REAL(DBL) :: h2           ! h**2 (square of the stepwidth in s-space)
    REAL(DBL) :: etot         ! total computed energy from last iteration
    REAL(DBL) :: ldd          ! logarhytmic derivative difference
         
    INTEGER :: l, n           ! angular and main quantum numbers
    INTEGER :: orbit          ! orbital number orbit = (n*(n-1))/2 + l + 1
    INTEGER :: nodes          ! number of nodes in the wavefunction

    
    INTEGER :: errno
    INTEGER :: k              ! loop variables
    INTEGER :: iteration      ! keep track of the number of iterations done
    INTEGER :: matchp         ! matching point (the coordinate where the two shoots should meet)
    INTEGER, DIMENSION(55) :: shoot_loop_counter, node_correction_counter

    LOGICAL :: itermax        ! to tell if the maximum number of iterations were reached

    ! ------------------------------------------------

    IF (this%log_wf) THEN
       OPEN (UNIT=10, FILE=this%wf_logfile, STATUS='REPLACE', ACTION='WRITE', IOSTAT=errno)
    ENDIF

    ! initialisations
    a = this%lddrsv * (2.0_dbl/pi)
    etot = 0.0_dbl
    h2 = this%rg%h**2
    iteration = 0
    itermax = .FALSE.
    np = this%rg%n
    node_correction_counter = 0
    shoot_loop_counter = 0
    t1 = 0.0_dbl; t2 = 0.0_dbl

    ! construct the initial wavefunctions
    this%orbitals = 0.0_dbl
    DO orbit=1, 55
       IF (this%elec_conf(orbit) /= 0) THEN
          this%orbitals(:,orbit) = &
               EXP(-1.0_dbl*(1.5_dbl*this%rg%r(1:np)/DBLE(orbit))**2) &
               * (2.25_dbl/SQRT(4.0_dbl*pi)) * this%rg%r(1:np)**2 &
               / DBLE(orbit)**2
          CALL rgNormalize(this%orbitals(:,orbit), this%rg)
       ENDIF
    ENDDO
    
    IF (this%method == 'XA') THEN
       CALL calculateElectronDensity( &
            rho, &
            this%orbitals(1:np,1:55), &
            DBLE(this%elec_conf), &
            1.0_dbl / this%rg%r(1:np) )
    ENDIF
    
    ! approximations for the eigenvalues
    CALL initialEigenvalues(this%eigenvalues, this)

    ! self-consistency loop
    sc_loop: DO 
       
       IF (iteration == this%maxiter) THEN
          itermax = .TRUE.
          EXIT sc_loop
       END IF
       iteration = iteration + 1

       WRITE (6,*)
       WRITE (6,"(2X,'Iteration ',I3)") iteration

       SELECT CASE (this%method)

       CASE DEFAULT
          WRITE (6,'(1X,3A)') "Method (", TRIM(this%method), ") is not available."
          EXIT sc_loop

       CASE ('IPM')
          WRITE (6,'(1X,3A)') "Method (", TRIM(this%method), ") not yet implemented."
          EXIT sc_loop

       CASE ('SIMPLE')
          ! nothing to do here

       CASE ('XA')
          CALL calculateElectronDensity( &
               new_rho, &
               this%orbitals(1:np,1:55), &
               DBLE(this%elec_conf), &
               1.0_dbl / this%rg%r(1:np) )
          rho = this%dma*new_rho + (1.0_dbl - this%dma)*rho

       ENDSELECT
       
       n_loop: DO n = 1, 10
          l_loop: DO l=0, n-1

             orbit = (n*(n-1))/2 + l + 1

             IF (this%elec_conf(orbit) == 0) CYCLE l_loop

             ! construct the model dependend part of the potential operator
             CALL calculatePotential(f, n, l, this, rho)          

             ! to find best matching point: find point with
             ! v - e = 0
             DO k=2, np
                IF ( signum(f(k)-this%eigenvalues(orbit)) /= &
                     signum(f(k-1)-this%eigenvalues(orbit)) ) THEN
                   matchp = k
                   EXIT
                ENDIF
             ENDDO
             
             ! transform the potential to the function in s-space
             f(1:np) = -(this%rg%wdd(1:np)) / this%rg%w(1:np) &
                  - (2.0_dbl / (this%rg%w(1:np))**4) * f(1:np)
    
             shoot_loop_counter(orbit) = 0
             shooting_loop: DO

                shoot_loop_counter(orbit) = shoot_loop_counter(orbit) + 1
                IF (shoot_loop_counter(orbit) > this%maxshoot) EXIT shooting_loop

                fprime = -1.0_dbl *(f + 2.0_dbl*this%eigenvalues(orbit)/(this%rg%w(1:np)**4))
                   
                t1(0) = 0.0_dbl ; t1(1) = 1.0E-4_dbl
                t2(np+1) = 0.0_dbl ; t2(np) = 0.5E-4_dbl
                
                ! numerov iteration
                DO k=2, (matchp+1)
                   t1(k) = t1(k-1) &
                        * (2.0_dbl+fprime(k-1)*h2/(1.0_dbl-fprime(k-1)*h2/12.0_dbl)) &
                        - t1(k-2)
                ENDDO
                DO k=np-1, (matchp-1), -1
                   t2(k) = t2(k+1) &
                        * (2.0_dbl+fprime(k+1)*h2/(1.0_dbl-fprime(k+1)*h2/12.0_dbl)) &
                        - t2(k+2)
                ENDDO

                ! get w*P
                DO k=1, matchp+1
                   t1(k) = t1(k) / (1.0_dbl-h2*fprime(k)/12.0_dbl)
                ENDDO
                DO k=np, (matchp-1), -1
                   t2(k) = t2(k) / (1.0_dbl-h2*fprime(k)/12.0_dbl)
                ENDDO
                
                ! rescale t2 so that they have the same value at the matching point
                t2 = (t1(matchp) / t2(matchp)) * t2

                ! extract P (t contains wP)
                t1(1:np) = t1(1:np) / this%rg%w(1:np)
                t2(1:np) = t2(1:np) / this%rg%w(1:np)

                ! compute the logarhitmic derivative difference
                ldd = (t1(matchp+1) - t1(matchp-1)) / ABS(t1(matchp)) - &
                     (t2(matchp+1) - t2(matchp-1)) / ABS(t2(matchp))
                ldd = ldd / ((this%rg%r(matchp+1)-this%rg%r(matchp-1))*2.0_dbl)

                this%orbitals(1:matchp,orbit) = t1(1:matchp)
                this%orbitals(matchp+1:np,orbit) = t2(matchp+1:np)
                CALL rgNormalize(this%orbitals(:,orbit), this%rg)

                IF (this%orbitals(1,orbit) < 0.0_dbl) &
                     this%orbitals(:,orbit) = this%orbitals(:,orbit) * (-1.0_dbl)

                ! number of nodes in wavefunction should be n-(l+1)
                CALL countNodes(nodes, this%orbitals(:,orbit))
                IF (nodes /= n-(l+1)) THEN
                   IF (nodes > n-(l+1)) THEN
                      this%eigenvalues(orbit) = this%eigenvalues(orbit) * (1.0_dbl + this%ncv)
                   ELSE
                      this%eigenvalues(orbit) = this%eigenvalues(orbit) * (1.0_dbl - this%ncv)
                   END IF
                   node_correction_counter(orbit) = node_correction_counter(orbit) + 1
                   CYCLE shooting_loop
                END IF

                ! exit condition
                IF (ABS(ldd) < this%precision**1.5) THEN
                   EXIT shooting_loop
                ENDIF
                
                !! "crush" the logarhitmic derivative difference
                ldd = a * ATAN(ldd/a)

                this%eigenvalues(orbit) = this%eigenvalues(orbit) + &
                     ldd * signum(this%orbitals(matchp,orbit)) * &
                     ABS(this%eigenvalues(orbit)) * &
                     (DBLE(orbit) / DBLE(this%z))**1.5 * & 
                     0.3_dbl

             END DO shooting_loop
             
             WRITE (6,"(2X,'n=',I1,' l=',I1,'  ',F15.9,' mp=',I4)") &
                  n, l, this%eigenvalues(orbit), matchp

          END DO l_loop
       END DO n_loop

       CALL getGroundStateEnergy(this, rho)

       ! these methods dont need the self consistency loop
       IF (this%method == 'SIMPLE' .OR. this%method == 'IPM') EXIT

       IF (iteration == 1) etot = 2.0_dbl * this%ground_state_energy

       ! self consistency condition
       IF ( ABS((this%ground_state_energy - etot) / this%ground_state_energy) < &
            this%precision) EXIT sc_loop
       etot = this%ground_state_energy

    END DO sc_loop

    WRITE (6,*)
    WRITE (6,*)
    DO n=1, 9
       DO l=0, n-1
          orbit = (n*(n-1))/2 + 1 + l
          IF (this%elec_conf(orbit) /= 0) THEN
             WRITE (6, "(' n = ',I1,' l = ',I1,' : ',F18.10,' lc = ',I5,' nc = ',I6)") &
                  n, l, this%eigenvalues(orbit), &
                  shoot_loop_counter(orbit), node_correction_counter(orbit)
          ENDIF
       ENDDO
    ENDDO
    WRITE (6,*)
    
    WRITE (6,'(1X,A,F18.10)') "GSE: ", this%ground_state_energy
    WRITE (6,*)
    WRITE (6,'(1X,A,I3)') "Number of iterations needed: ", iteration
    IF (itermax) THEN
       WRITE (6,'(3X,A)') "The simulation stopped because the maximum number of "
       WRITE (6,'(3X,A)') "iterations were reached!"
    ENDIF
    WRITE (6,*)

    IF (this%log_wf) THEN
       DO k=1, np
          WRITE (10,'(5(2X,F14.8))') this%rg%r(k), &
               this%orbitals(k,3), this%orbitals(k,0), this%orbitals(k,0), this%orbitals(k,0)
       ENDDO
       CLOSE (10)
   ENDIF

  END SUBROUTINE atom_do_calc
  
END MODULE atom_module

!******************************************************************************
