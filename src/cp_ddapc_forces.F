!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

MODULE cp_ddapc_forces
!!****h* cp2k/cp_ddapc_forces *
!!
!!   NAME
!!    cp_ddapc_forces
!!
!!   FUNCTION
!!     Density Derived atomic point charges from a QM calculation
!!     (see J. Chem. Phys. Vol. 103 pp. 7422-7428)
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     08.2005 created [tlaino]
!!
!!   SOURCE
!****************************************************************************  
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE cell_types,                      ONLY: cell_type
  USE cp2k_units,                      ONLY: get_cp2k_units,&
                                             unit_convert_type
  USE erf_fn,                          ONLY: erfc
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi,&
                                             rootpi
  USE particle_types,                  ONLY: particle_type
  USE pw_spline_utils,                 ONLY: Eval_d_Interp_Spl3_pbc
  USE pw_types,                        ONLY: pw_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_ddapc_forces'
  PUBLIC :: ewald_ddapc_force, reset_ch_pulay
  
CONTAINS

!!****f* cp_ddapc_forces/ewald_ddapc_force [1.0] *
!!
!!   NAME
!!     ewald_ddapc_force
!!
!!   FUNCTION
!!     Evaluates the Ewald term E2 and E3 energy term for the decoupling/coupling
!!     of periodic images
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     08.2005 created [tlaino]
!!
!!*** **********************************************************************
  RECURSIVE SUBROUTINE ewald_ddapc_force(qs_env, coeff, qmmm_periodic, multipole_section, &
       particle_set, radii, dq, charges, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(pw_type), POINTER                   :: coeff
    LOGICAL, INTENT(IN)                      :: qmmm_periodic
    TYPE(section_vals_type), POINTER         :: multipole_section
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), DIMENSION(:), POINTER     :: radii
    REAL(KIND=dp), DIMENSION(:, :, :), &
      OPTIONAL, POINTER                      :: dq
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: charges
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ewald_ddapc_force', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, ip1, ip2, iparticle1, iparticle2, k1, k2, k3, n_rep, &
      nmax1, nmax2, nmax3, r1, r2, r3, rmax1, rmax2, rmax3, sfact, stat
    LOGICAL                                  :: analyt, failure
    REAL(KIND=dp)                            :: alpha, dfact, eps, eps0, fac, &
                                                fac2, fac3, fs, galpha, gsq, &
                                                gsqi, q1t, q2t, r, rcut, t1, &
                                                t2, tol, tol1
    REAL(KIND=dp), DIMENSION(3)              :: drvec, fvec, gvec, ra, rb, &
                                                rvec
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: d_el, M
    TYPE(cell_type), POINTER                 :: cell
    TYPE(unit_convert_type), POINTER         :: units

    failure = .FALSE.
    NULLIFY(d_el, units, M)
    CALL timeset(routineN,"I","",handle)
    IF (.NOT.failure) THEN
       CPPostcondition(PRESENT(charges),cp_failure_level,routineP,error,failure)
       CPPostcondition(ASSOCIATED(radii),cp_failure_level,routineP,error,failure)
       CALL get_cp2k_units ( units )
       sfact =  1.0_dp
       IF (qmmm_periodic) THEN 
          sfact = -1.0_dp
          CALL get_qs_env(qs_env, super_cell=cell)
       ELSE
          CALL get_qs_env(qs_env, cell=cell)
       END IF
       rcut = MIN(cell%hmat(1,1),cell%hmat(2,2),cell%hmat(3,3))/2.0_dp
       CALL section_vals_val_get(multipole_section,"RCUT",n_rep_val=n_rep,error=error)
       IF (n_rep==1) CALL section_vals_val_get(multipole_section,"RCUT",r_val=rcut,error=error)
       CALL section_vals_val_get(multipole_section,"EWALD_PRECISION",r_val=eps,error=error)
       CALL section_vals_val_get(multipole_section,"ANALYTICAL_GTERM",l_val=analyt,error=error)
       !
       ! Setting-up parameters for Ewald summation
       !
       eps0   = units%eps0
       eps    = MIN(ABS(eps),0.5_dp)
       tol    = SQRT(ABS(LOG(eps*rcut)))
       alpha  = SQRT(ABS(LOG(eps*rcut*tol)))/rcut
       galpha = 1.0_dp/(4.0_dp*alpha*alpha)
       tol1   = SQRT(-LOG(eps*rcut*(2.0_dp*tol*alpha)**2))
       nmax1  = NINT(0.25_dp + cell%hmat(1,1)*alpha*tol1/pi)
       nmax2  = NINT(0.25_dp + cell%hmat(2,2)*alpha*tol1/pi)
       nmax3  = NINT(0.25_dp + cell%hmat(3,3)*alpha*tol1/pi)

       rmax1  = CEILING(rcut/cell%hmat(1,1))
       rmax2  = CEILING(rcut/cell%hmat(2,2))
       rmax3  = CEILING(rcut/cell%hmat(3,3))

       ALLOCATE(d_el(3,SIZE(particle_set)), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       d_el = 0.0_dp
       fac    = 1.d0/cell%deth
       fac2   = 1.0_dp / ( 4.0_dp * pi * eps0 )
       fac3   = fac/8.0_dp
       fvec   = 2.0_dp * pi / (/cell%hmat(1,1),cell%hmat(2,2),cell%hmat(3,3)/)
       !
       DO iparticle1 = 1, SIZE(particle_set)
          ip1 = (iparticle1-1)*SIZE(radii)          
          q1t   = SUM(charges(ip1+1:ip1+SIZE(radii)))
          DO iparticle2 = 1, iparticle1
             ip2 = (iparticle2-1)*SIZE(radii)
             q2t   = SUM(charges(ip2+1:ip2+SIZE(radii)))
             !
             ! Real-Space Contribution
             !
             rvec   = particle_set(iparticle1)%r - particle_set(iparticle2)%r
             dfact = 1.0_dp
             IF (iparticle1 == iparticle2) dfact = 0.5_dp
             IF (iparticle1 /= iparticle2) THEN
                ra = rvec 
                r  = SQRT(DOT_PRODUCT(ra,ra))
                IF (r<=rcut) THEN 
                   t1 = erfc(alpha*r) / r
                   drvec = ra / r * q1t * q2t * fac2 * sfact
                   t2 = -2.0_dp * alpha * EXP(-alpha*alpha*r*r) / (r*rootpi) - t1/r
                   d_el(1:3,iparticle1) = d_el(1:3,iparticle1) - t2 * drvec
                   d_el(1:3,iparticle2) = d_el(1:3,iparticle2) + t2 * drvec
                END IF
             END IF
             DO r1 = -rmax1, rmax1
                DO r2 = -rmax2, rmax2
                   DO r3 = -rmax3, rmax3
                      IF ((r1==0).AND.(r2==0).AND.(r3==0)) CYCLE
                      rb = (/REAL(r1,KIND=dp),REAL(r2,KIND=dp),REAL(r3,KIND=dp)/)
                      ra = rvec + MATMUL(cell%hmat,rb) 
                      r  = SQRT(DOT_PRODUCT(ra,ra))
                      IF (r<=rcut) THEN 
                         t1 = erfc(alpha*r) / r
                         drvec = ra / r * q1t * q2t * fac2 * dfact * sfact
                         t2 = -2.0_dp * alpha * EXP(-alpha*alpha*r*r) / (r*rootpi) - t1/r
                         d_el(1:3,iparticle1) = d_el(1:3,iparticle1) - t2 * drvec
                         d_el(1:3,iparticle2) = d_el(1:3,iparticle2) + t2 * drvec
                      END IF
                   END DO
                END DO
             END DO
             !
             ! G-space Contribution
             !
             IF (analyt) THEN
                DO k1 = 0, nmax1
                   DO k2 = -nmax2, nmax2
                      DO k3 = -nmax3, nmax3
                         IF (k1 == 0.AND.k2 == 0.AND.k3 == 0) CYCLE
                         fs = 2.0_dp; IF (k1==0) fs = 1.0_dp
                         gvec = fvec * (/REAL(k1,KIND=dp),REAL(k2,KIND=dp),REAL(k3,KIND=dp)/)
                         gsq  = DOT_PRODUCT(gvec,gvec)
                         gsqi = fs/gsq
                         t1   = fac * gsqi * EXP(-galpha*gsq)
                         t2 = -SIN(DOT_PRODUCT(gvec,rvec)) * t1 * q1t * q2t / eps0 * dfact * sfact
                         d_el(1:3,iparticle1) = d_el(1:3,iparticle1) - t2 * gvec
                         d_el(1:3,iparticle2) = d_el(1:3,iparticle2) + t2 * gvec
                      END DO
                   END DO
                END DO
             ELSE
                gvec = Eval_d_Interp_Spl3_pbc( rvec, coeff, error ) * q1t * q2t / eps0 * dfact * sfact
                d_el(1:3,iparticle1) = d_el(1:3,iparticle1) -  gvec
                d_el(1:3,iparticle2) = d_el(1:3,iparticle2) +  gvec                
             END IF
             IF (iparticle1/=iparticle2) THEN 
                ra = rvec 
                r  = SQRT(DOT_PRODUCT(ra,ra))
                t2    = - 1.0_dp /(r*r) * sfact
                drvec = ra / r * q1t * q2t 
                d_el(1:3,iparticle1) = d_el(1:3,iparticle1) + t2 * drvec
                d_el(1:3,iparticle2) = d_el(1:3,iparticle2) - t2 * drvec
             END IF
          END DO
       END DO
       M => qs_env%cp_ddapc_env%Md
       IF (qmmm_periodic) M => qs_env%cp_ddapc_env%Mr
       CALL cp_decpl_ddapc_forces(qs_env, M, charges, dq, d_el, particle_set, error)
       DEALLOCATE(d_el, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)          
    END IF
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE ewald_ddapc_force


!!****f* cp_ddapc_forces/cp_decpl_ddapc_forces [1.0] *
!!
!!   NAME
!!     cp_decpl_ddapc_forces
!!
!!   FUNCTION
!!     Evaluation of the pulay forces due to the fitted charge density
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     08.2005 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE cp_decpl_ddapc_forces(qs_env, M, charges, dq, d_el, particle_set, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: M
    REAL(KIND=dp), DIMENSION(:), POINTER     :: charges
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: dq
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: d_el
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_decpl_ddapc_forces', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iatom, ikind, j, &
                                                k, natom, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind, kind_of
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: uv
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: chf
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force

    failure = .FALSE.
    CALL timeset(routineN,"I","",handle)
    IF (.NOT.failure) THEN
       natom = SIZE(particle_set)
       CALL get_qs_env(qs_env=qs_env,&
                       atomic_kind_set=atomic_kind_set,&
                       force=force)
       ALLOCATE (atom_of_kind(natom),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE (kind_of(natom),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE (chf(3,natom),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                                atom_of_kind=atom_of_kind,&
                                kind_of=kind_of)

       ALLOCATE(uv(SIZE(M,1)), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       uv   = MATMUL(M,charges)
       DO k = 1, natom
          DO j = 1, 3
             chf(j,k) =  DOT_PRODUCT(uv,dq(:,k,j))
          END DO
       END DO
       DO iatom=1,natom
          ikind = kind_of(iatom)
          i = atom_of_kind(iatom)
          force(ikind)%ch_pulay(1:3,i) = force(ikind)%ch_pulay(1:3,i) + chf(1:3,iatom) + d_el(1:3,iatom) 
       END DO
       DEALLOCATE(atom_of_kind, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(kind_of, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(chf, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(uv, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL timestop(0.0_dp, handle)
  END SUBROUTINE cp_decpl_ddapc_forces

!!****f* cp_ddapc_forces/reset_ch_pulay [1.0] *
!!
!!   NAME
!!     reset_ch_pulay
!!
!!   FUNCTION
!!     Evaluation of the pulay forces due to the fitted charge density
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     08.2005 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE reset_ch_pulay(qs_env, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'reset_ch_pulay', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ind
    LOGICAL                                  :: failure
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force

    failure = .FALSE.

    CALL timeset(routineN,"I","",handle)
    IF (.NOT.failure) THEN
       CALL get_qs_env(qs_env=qs_env,&
                       force=force)
       DO ind = 1, SIZE(force)
          force(ind)%ch_pulay = 0.0_dp
       END DO
    END IF
    CALL timestop(0.0_dp, handle)
  END SUBROUTINE reset_ch_pulay

END MODULE cp_ddapc_forces
