!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****m* cp2k/atom_module [1.0] *
!!
!!   NAME
!!     atom_module
!!
!!   FUNCTION
!!     A module to compute the ground state energy of atoms. The atom
!!     to be calculated must be specified in the input file. This
!!     module uses the functionals defined in functionals.F.
!!
!!   SYNOPSIS
!!     atomic_code ( globenv )
!!
!!   MODIFICATION HISTORY
!!     JGH, 17-NOV-2000, began work
!!     TCH, 21-DEC-2001, version 0.1 (Xa works)
!!     TCH, 12-FEB-2002, working version
!!
!!   SEE ALSO
!!     atom_calculations, atom_type, atom_parameters_type
!!
!!****
MODULE atom_module

  USE atom_calculations, ONLY : &
       calculate_electron_density, &
       calculate_gs_energy, calculate_potential, &
       count_nodes, initialize_eigenvalues, rg_normalize, signum
  USE atom_parameters_type, ONLY : AM_PARAMETERS_INFO, am_init, am_dealloc_ref, &
       am_get_dma, am_get_method, am_get_maxshoot, am_get_ncv, &
       am_get_precision, am_get_maxiter, am_get_ldd_rescale_value
  USE atom_type, ONLY : AM_ATOM_INFO, am_init, am_dealloc_ref, &
       am_get_elec_conf, am_get_z
  USE global_types, ONLY : global_environment_type
  USE header, ONLY : atom_header
  USE kinds, ONLY : dbl
  USE mathconstants, ONLY : zero, pi
  USE cp_error_handling, ONLY : cp_assert, cp_error_type, cp_init, cp_dealloc_ref, &
       cp_note_level, cp_warning_level, cp_failure_level, cp_fatal_level
  USE cp_log_handling, ONLY : cp_init_default_loggers, cp_to_string
  USE radial_util, ONLY : radial_integrate
  USE termination, ONLY : stop_memory, stop_program
  USE timings, ONLY : timeset, timestop

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: atomic_code

  !-----------------------------------------------------------------------------!

CONTAINS

!!****f* cp2k/atom_module/atomic_code [1.0] *
!!
!!   NAME
!!     atomic_code
!!
!!   FUNCTION
!!     This is the public subroutine of the atom_module module. It computes the
!!     ground state energy and writes it on the screen. The main part of the work
!!     is delegated to atom_do_calc.
!!
!!   ARGUMENTS
!!     globenv - TYPE ( GLOBAL_ENVIRONMENT_TYPE ) - contains the environment
!!
!!   MODIFICATION HISTORY
!!     JGH, 17-NOV-2000, began work
!!     TCH, 02-FEB-2002, working version
!!
!!   SEE ALSO
!!     atom_do_calc
!!
!!****
  SUBROUTINE atomic_code(globenv)

    IMPLICIT NONE

    ! Arguments
    TYPE (GLOBAL_ENVIRONMENT_TYPE), INTENT(INOUT) :: globenv

    ! Locals
    TYPE (AM_ATOM_INFO) :: atom
    TYPE (CP_ERROR_TYPE) :: error
    TYPE (AM_PARAMETERS_INFO) :: parameters

    INTEGER :: handle ! for the timer
    LOGICAL :: condition, failure

    !------------------------------------------------------------------------------
    CALL timeset('ATOMIC_CODE', 'I', ' ', handle)

    IF (globenv%ionode) CALL atom_header(globenv%scr)

    CALL am_init(parameters, globenv)
    CALL am_init(atom, globenv)
    CALL cp_init_default_loggers
    CALL cp_init(error)
    
    CALL atom_do_calc ( atom, parameters )

    CALL timestop ( zero, handle )

    CALL cp_dealloc_ref ( error )
    CALL am_dealloc_ref ( atom )
    CALL am_dealloc_ref ( parameters )

  END SUBROUTINE atomic_code

  !-----------------------------------------------------------------------------!

!!****f* cp2k/atom_module/atom_do_calc [1.0] *
!!
!!   NAME
!!     atom_do_calc
!!
!!   FUNCTION
!!     Calculates the ground state electron density and energy.
!!
!!   ARGUMENTS
!!     atom - TYPE ( AM_ATOM_INFO ) - the atom to compute
!!     parameters - TYPE ( AM_PARAMETERS_INFO ) - an object holding program parameters
!!
!!   MODIFICATION HISTORY
!!     08-FEB-2002, TCH, working version
!!
!!   SEE ALSO
!!     AM_ATOM_INFO, AM_PARAMETERS_INFO
!!
!!****
  SUBROUTINE atom_do_calc ( atom, parameters )
    
    IMPLICIT NONE
    
    ! Arguments
    TYPE (AM_ATOM_INFO), INTENT (INOUT) :: atom
    TYPE (AM_PARAMETERS_INFO), INTENT (IN) :: parameters

    ! Variables    
    INTEGER :: errno        
    INTEGER :: np
    INTEGER :: l, n, orbit  ! orbit = (n*(n-1))/2 + l + 1
    INTEGER :: nodes
    INTEGER :: iteration
    INTEGER :: k
    INTEGER :: matchp       ! matching point (the coordinate where the two shoots should meet)
    INTEGER, DIMENSION(55) :: shoot_loop_counter, node_correction_counter

    LOGICAL :: itermax      ! to tell if the maximum number of iterations were reached

    REAL(DBL) :: a          ! parameter to squish the ldd to a useable number
    REAL(DBL) :: h2         ! h**2 (square of the stepwidth in s-space)
    REAL(DBL) :: gse_li     ! ground state energy from last iteration
    REAL(DBL) :: ldd        ! logarhytmic derivative difference

    REAL(DBL), DIMENSION(1:atom%rg%n) :: f, fprime
    REAL(DBL), DIMENSION(0:atom%rg%n+1) :: t1, t2
    REAL(DBL), DIMENSION(1:atom%rg%n) :: rho, new_rho ! electron density
    REAL(DBL), DIMENSION(1:atom%rg%n,55) :: wf

    ! ------------------------------------------------

    IF (atom%log_wf) THEN
       OPEN (UNIT=10, FILE=atom%wf_logfile, STATUS='REPLACE', ACTION='WRITE', IOSTAT=errno)
    ENDIF

    ! variable initialisations
    np = atom%rg%n
    gse_li = 0.0_dbl
    t1 = 0.0_dbl; t2 = 0.0_dbl
    a = am_get_ldd_rescale_value ( parameters ) * (2.0_dbl/pi)
    h2 = atom%rg%h**2
    iteration = 0 ; itermax = .FALSE.
    node_correction_counter = 0
    shoot_loop_counter = 0

    IF (am_get_method(parameters) == 'XA') THEN
       DO orbit=1, 55
          IF (am_get_elec_conf(atom, orbit) /= 0) &
               wf(1:np,orbit) = atom%orbitals(1:np,orbit) / atom%rg%r(1:np)
       END DO
       CALL calculate_electron_density(rho, wf, am_get_elec_conf(atom))
    END IF

    ! approximations for the eigenvalues
    CALL initialize_eigenvalues ( atom, parameters, atom%eigenvalues )

    ! self-consistency loop
    sc_loop: DO 
       
       IF ( iteration == am_get_maxiter ( parameters ) ) THEN
          itermax = .TRUE.
          EXIT sc_loop
       END IF
       iteration = iteration + 1

       WRITE (6,*)
       WRITE (6,"(2X,'Iteration ',I3)") iteration

       SELECT CASE (am_get_method(parameters))

       CASE DEFAULT
          WRITE (6,'(1X,3A)') "Method (", TRIM(am_get_method(parameters)), ") is not available."
          EXIT sc_loop

       CASE ('IPM')
          WRITE (6,'(1X,3A)') "Method (", TRIM(am_get_method(parameters)), ") not yet implemented."
          EXIT sc_loop

       CASE ('SIMPLE')
          ! nothing to do here

       CASE ('XA')
          DO orbit=1, 55
             IF (am_get_elec_conf(atom,orbit) /= 0) &
                  wf(1:np,orbit) = atom%orbitals(1:np,orbit) / atom%rg%r(1:np)
          END DO
          CALL calculate_electron_density(new_rho, wf, am_get_elec_conf(atom))
          rho = am_get_dma(parameters) * new_rho + (1.0_dbl - am_get_dma(parameters))*rho

       ENDSELECT
       
       n_loop: DO n = 1, 10
          l_loop: DO l=0, n-1

             orbit = (n*(n-1))/2 + l + 1

             IF (am_get_elec_conf(atom, orbit) == 0) CYCLE l_loop

             ! construct the model dependend part of the potential operator
             CALL calculate_potential(f, atom, parameters, n, l, rho)          

             ! to find best matching point: find point with
             ! v - e = 0
             DO k=2, np
                IF ( signum(f(k) - atom%eigenvalues(orbit)) /= &
                     signum(f(k-1) - atom%eigenvalues(orbit))) THEN
                   matchp = k
                   EXIT
                ENDIF
             ENDDO
             
             ! transform the potential to the function in s-space
             f(1:np) = - atom%rg%wdd(1:np) / atom%rg%w(1:np) &
                  - (2.0_dbl / atom%rg%w(1:np)**4) * f(1:np)
    
             shoot_loop_counter(orbit) = 0
             shooting_loop: DO

                shoot_loop_counter(orbit) = shoot_loop_counter(orbit) + 1
                IF (shoot_loop_counter(orbit) > am_get_maxshoot(parameters)) EXIT shooting_loop

                fprime = -1.0_dbl *(f + 2.0_dbl*atom%eigenvalues(orbit)/(atom%rg%w(1:np)**4))
                   
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
                t1(1:np) = t1(1:np) / atom%rg%w(1:np)
                t2(1:np) = t2(1:np) / atom%rg%w(1:np)

                ! compute the logarhitmic derivative difference
                ldd = (t1(matchp+1) - t1(matchp-1)) / ABS(t1(matchp)) - &
                     (t2(matchp+1) - t2(matchp-1)) / ABS(t2(matchp))
                ldd = ldd / ((atom%rg%r(matchp+1)-atom%rg%r(matchp-1))*2.0_dbl)

                atom%orbitals(1:matchp,orbit) = t1(1:matchp)
                atom%orbitals(matchp+1:np,orbit) = t2(matchp+1:np)
                CALL rg_normalize(atom%orbitals(:,orbit), atom%rg)

                IF (atom%orbitals(1,orbit) < 0.0_dbl) &
                     atom%orbitals(:,orbit) = atom%orbitals(:,orbit) * (-1.0_dbl)

                ! number of nodes in wavefunction should be n-(l+1)
                CALL count_nodes(nodes, atom%orbitals(:,orbit))
                IF (nodes /= n-(l+1)) THEN
                   IF (nodes > n-(l+1)) THEN
                      atom%eigenvalues(orbit) = atom%eigenvalues(orbit) &
                           * (1.0_dbl + am_get_ncv(parameters))
                   ELSE
                      atom%eigenvalues(orbit) = atom%eigenvalues(orbit) &
                           * (1.0_dbl - am_get_ncv(parameters))
                   END IF
                   node_correction_counter(orbit) = node_correction_counter(orbit) + 1
                   CYCLE shooting_loop
                END IF

                ! exit condition
                IF ( ABS ( ldd ) < am_get_precision ( parameters )**1.5 ) THEN
                   EXIT shooting_loop
                ENDIF
                
                !! "crush" the logarhitmic derivative difference
                ldd = a * ATAN(ldd/a)

                atom%eigenvalues(orbit) = atom%eigenvalues(orbit) + &
                     ldd * signum(atom%orbitals(matchp,orbit)) * &
                     ABS(atom%eigenvalues(orbit)) * &
                     (DBLE(orbit) / DBLE(atom%z))**1.5 * & 
                     0.3_dbl

             END DO shooting_loop
             
             WRITE (6,"(2X,'n=',I1,' l=',I1,'  ',F15.9,' mp=',I4)") &
                  n, l, atom%eigenvalues(orbit), matchp

          END DO l_loop
       END DO n_loop

       CALL calculate_gs_energy ( atom, parameters, rho )

       ! these methods dont need the self consistency loop
       IF ( am_get_method ( parameters )  == 'SIMPLE' .OR. am_get_method ( parameters ) == 'IPM') EXIT

       IF (iteration == 1) gse_li = 2.0_dbl * atom%ground_state_energy

       ! self consistency condition
       IF ( ABS((atom%ground_state_energy - gse_li) / atom%ground_state_energy) < &
            am_get_precision ( parameters ) ) EXIT sc_loop
       gse_li = atom%ground_state_energy

    END DO sc_loop

    WRITE (6,*)
    WRITE (6,*)
    DO n=1, 9
       DO l=0, n-1
          orbit = (n*(n-1))/2 + 1 + l
          IF ( am_get_elec_conf ( atom, orbit ) /= 0) THEN
             WRITE (6, "(' n = ',I1,' l = ',I1,' : ',F18.10,' lc = ',I5,' nc = ',I6)") &
                  n, l, atom%eigenvalues(orbit), &
                  shoot_loop_counter(orbit), node_correction_counter(orbit)
          ENDIF
       ENDDO
    ENDDO
    WRITE (6,*)
    
    WRITE (6,'(1X,A,F18.10)') "GSE: ", atom%ground_state_energy
    WRITE (6,*)
    WRITE (6,'(1X,A,I3)') "Number of iterations needed: ", iteration
    IF (itermax) THEN
       WRITE (6,'(3X,A)') "The simulation stopped because the maximum number of "
       WRITE (6,'(3X,A)') "iterations were reached!"
    ENDIF
    WRITE (6,*)

    IF (atom%log_wf) THEN
       DO k=1, np
          WRITE (10,'(5(2X,F14.8))') atom%rg%r(k), &
               atom%orbitals(k,3), atom%orbitals(k,0), atom%orbitals(k,0), atom%orbitals(k,0)
       ENDDO
       CLOSE (10)
   ENDIF

  END SUBROUTINE atom_do_calc
  
END MODULE atom_module

