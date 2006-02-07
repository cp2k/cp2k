!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****h* cp2k/atom_module [1.0] *
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
!!   SOURCE
!!
MODULE atom_module

  USE atom_types,                      ONLY: ATOM_CALCULATION_TYPE,&
                                             ATOM_TYPE,&
                                             init
  USE atom_util
  USE global_types,                    ONLY: global_environment_type
  USE header,                          ONLY: atom_header
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi,&
                                             zero
  USE radial_grids,                    ONLY: radial_grid_type
  USE radial_util,                     ONLY: radial_integrate,&
                                             radial_normalize
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

! global variables
  INTEGER                              :: np
  REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: wdd_w, spot, tl, tu
  
  PRIVATE
  PUBLIC :: atomic_code

! ****************************************************************************
  
CONTAINS

!!*****
!!****f* atom_module/atomic_code [1.0] *
!!
!!   NAME
!!     atomic_code
!!
!!   SYNOPSIS
!!     Subroutine atomic_code(globenv)
!!       Type(global_environment_type), pointer:: globenv
!!     End Subroutine atomic_code
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
!!     TCH, 20-JUN-2002, working version
!!
!!   SEE ALSO
!!
!!   SOURCE
!!
  SUBROUTINE atomic_code(globenv)


    TYPE(GLOBAL_ENVIRONMENT_TYPE), POINTER   :: globenv

    INTEGER                                  :: handle
    TYPE(ATOM_CALCULATION_TYPE)              :: ac
    TYPE(ATOM_TYPE)                          :: atom

! for the timer

    CALL timeset('ATOMIC_CODE', 'I', ' ', handle)
    
    IF (globenv%ionode) CALL atom_header(globenv%scr)
    
    CALL init(atom, ac, globenv)
    CALL init_globals(ac)
    CALL calculate(atom, ac)
    
    CALL timestop ( zero, handle )
    
  END SUBROUTINE atomic_code
!!*****
  
! ****************************************************************************
  
  SUBROUTINE calculate ( atom, ac )

    TYPE(ATOM_TYPE), INTENT(INOUT)           :: atom
    TYPE(ATOM_CALCULATION_TYPE), INTENT(IN)  :: ac

    INTEGER                                  :: err, it, l, matchp, max_iter, &
                                                max_shoot, n, nodes, si
    LOGICAL                                  :: shoot_converged
    REAL(KIND=dp)                            :: a, eg, eg_li, ev, ldd, &
                                                mix_ang, ncv, PRECISION, res
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: pot, rho, wf

! variables
! ground state enegy (and of last iteration)

    a = 5.0_dp
    ncv = 0.1_dp  ! node correction value
    max_iter = ac%maxiter
    max_shoot = ac%maxshoot ! maximum number of shooting iterations
    mix_ang = 0.4_dp
    PRECISION = ac%precision

    ALLOCATE (pot(1:np), rho(1:np), wf(1:np))

    ! initialize the electronic density
    IF (ac%method == 'XA') rho = calculate_rho(atom, np)

    it = 1
    self_consistent : DO

       IF (it > max_iter) EXIT self_consistent
       WRITE (6,*)
       WRITE (6,'(1X,A,I3)') "Iteration ", it

       ! rho = (1-ma)*rho + ma*new_rho
       IF (ac%method == 'XA') THEN
          CALL DSCAL(np, 1.0_dp - mix_ang, rho(1:np), 1)
          CALL DAXPY(np, mix_ang, calculate_rho(atom, np), 1, rho(1:np), 1)
! DEBUG
          CALL radial_integrate(np,rho*ac%rg%r(1:np)**2,ac%rg,res,err)
          IF (ABS(4.0_dp*pi*res - REAL(atom%ne, KIND=dp)) > 1.E-6) &
               WRITE (6,*) "!! Integral(rho) = ", 4.0_dp*pi*res
       END IF

       shoot_converged = .TRUE.
       DO n=1, atom%ns
          DO l=0, atom%shells(n)%no-1

             CALL calculate_potential(pot, n, l, rho, atom, ac)

             si = 1
             shooting : DO

                IF (si > max_shoot) THEN
                   shoot_converged = .FALSE.
                   WRITE (6,'(1X,3(A,I3),A)') "I ", it, "  N ", n, "  L ", l, " : didn't converge"
                   EXIT shooting
                END IF

                ev = atom%shells(n)%orbits(l)%ev
                matchp = find_matchp(pot, ev)
                ldd = shoot(wf, ev, pot, matchp, ac%rg)

                ! check number of nodes
                CALL count_nodes(nodes, wf)
                IF (nodes /= n-(l+1)) THEN
                   IF (nodes > n-(l+1)) THEN
                      CALL DGESS(1.0_dp + ncv, atom%shells(n)%orbits(l)%ev)
                   ELSE                                                 
                      CALL DGESS(1.0_dp - ncv, atom%shells(n)%orbits(l)%ev)
                   END IF
                   si = si + 1
                   CYCLE shooting
                END IF

                ! we only take the new wavefunction 
                ! if it has the correct number of nodes
                ! but still does this to often
                atom%shells(n)%orbits(l)%wf(1:np) = wf(1:np) / ac%rg%r(1:np)
                CALL radial_normalize(atom%shells(n)%orbits(l)%wf(1:np), ac%rg)

                ! exit condition
                IF ( ABS(ldd) < 1.E-5_dp ) EXIT shooting
                
                ! otherwise calculate new ev and reshoot
                ! "crush" the logarhitmic derivative difference
                ldd = a * ATAN(ldd/a)
                CALL DGESS(1.0_dp - 0.3_dp*ldd*signum(wf(matchp)) * &
                     (REAL((n*(n-1))/2 + l + 1,KIND=dp) / REAL(atom%z,KIND=dp))**1.5, &
                     atom%shells(n)%orbits(l)%ev)
                                     
                si = si + 1

             END DO shooting
             
          END DO
       END DO

       IF (.NOT.shoot_converged) THEN
          WRITE (6,*) "WARNING: At least one of the shooting loops didn't converge!"
          ! maybe do something here
       END IF

       eg = calculate_energy( atom, ac, rho )
       WRITE (6,'(1X,A,F14.6)') "ground state energy : ", eg
       IF (ABS(eg - eg_li) < 1.E-4_dp) THEN
          EXIT self_consistent
       END IF
       eg_li = eg

       it = it + 1
       
    END DO self_consistent

    IF (it > max_iter) THEN
       WRITE (6,*) "WARNING: Self consistency not reached withing specified iterations!"
       ! maybe do something here
    END IF

    ! write solution
    WRITE (6,'(1X,A,F14.6)') "ground state energy : ", eg

    DEALLOCATE(wf)

  END SUBROUTINE calculate

! ****************************************************************************
  
  FUNCTION shoot(wf, ev, pot, matchp, rg)
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: wf
    REAL(KIND=dp), INTENT(IN)                :: ev
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: pot
    INTEGER, INTENT(IN)                      :: matchp
    TYPE(radial_grid_type), INTENT(IN)       :: rg
    REAL(KIND=dp)                            :: shoot

    INTEGER                                  :: k
    REAL(KIND=dp)                            :: h2, h2_12, scale

! function
! variables

    h2 = rg%h * rg%h
    h2_12 = h2 / 12.0_dp

    tl(0) = 0.0_dp    ; tl(1) = 1.0E-4_dp
    tu(np+1) = 0.0_dp ; tu(np) = 0.5E-4_dp

    DO k=1, np
       spot(k) = wdd_w(k) - 2.0_dp*(ev - pot(k))/rg%w(k)**4
    END DO

    ! numerov iteration
    DO k=1, matchp
       ! calculate spot at k
       tl(k+1) = tl(k) * (2.0_dp+spot(k)*h2/(1.0_dp-spot(k)*h2_12)) - tl(k-1)
    ENDDO
    DO k=np, matchp, -1
       ! calculate spot at k+1
       tu(k-1) = tu(k) * (2.0_dp+spot(k)*h2/(1.0_dp-spot(k)*h2_12)) - tu(k+1)
    ENDDO

    ! get w*P
    DO k=1, matchp+1
       tl(k) = tl(k) / (1.0_dp-h2_12*spot(k))
    ENDDO
    DO k=np, (matchp-1), -1
       tu(k) = tu(k) / (1.0_dp-h2_12*spot(k))
    ENDDO
                
    scale = tl(matchp) / tu(matchp)
    CALL DSCAL(np+1-(matchp-1)+1, scale, tu(matchp-1:np+1), 1)

    tl(1:matchp+1) = tl(1:matchp+1) / rg%w(1:matchp+1)
    tu(matchp-1:np) = tu(matchp-1:np) / rg%w(matchp-1:np)
    
    ! compute the logarhitmic derivative difference
    shoot = (tl(matchp+1) - tl(matchp-1)) / ABS(tl(matchp)) - &
          (tu(matchp+1) - tu(matchp-1)) / ABS(tu(matchp))
    shoot = shoot / ((rg%r(matchp+1)-rg%r(matchp-1))*2.0_dp)
    
    CALL DCOPY(matchp, tl(1:matchp), 1, wf(1:matchp), 1)
    CALL DCOPY(np-(matchp+1)+1, tu(matchp+1:np), 1, wf(matchp+1:np), 1)

    ! make the wf start positively
    IF (signum(wf(1)) == -1) CALL DSCAL(np, -1.0_dp, wf, 1)

  END FUNCTION shoot

! ****************************************************************************
  
  SUBROUTINE init_globals(ac)
    TYPE(atom_calculation_type), INTENT(IN)  :: ac

    INTEGER                                  :: err

    np = ac%rg%n
    ALLOCATE (wdd_w(1:np), spot(1:np), tl(0:np+1), tu(0:np+1), stat=err)
    IF (err /= 0) THEN
       CALL stop_program("atom_module/init_globals","allocate somehow didn't work")
    END IF
    wdd_w(1:np) = ac%rg%wdd(1:np)/ac%rg%w(1:np)
  END SUBROUTINE init_globals

! ****************************************************************************
  
END MODULE atom_module
