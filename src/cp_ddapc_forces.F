!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Density Derived atomic point charges from a QM calculation
!>      (see J. Chem. Phys. Vol. 103 pp. 7422-7428)
!> \par History
!>      08.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
MODULE cp_ddapc_forces
  
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE cell_types,                      ONLY: cell_type
  USE cp_control_types,                ONLY: ddapc_restraint_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE erf_fn,                          ONLY: erfc
  USE input_constants,                 ONLY: do_ddapc_constraint,&
                                             do_ddapc_restraint,&
                                             weight_type_mass,&
                                             weight_type_unit
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: fourpi,&
                                             pi,&
                                             rootpi,&
                                             twopi
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE pw_spline_utils,                 ONLY: Eval_d_Interp_Spl3_pbc
  USE pw_types,                        ONLY: pw_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE spherical_harmonics,             ONLY: dlegendre,&
                                             legendre
  USE timings,                         ONLY: timeset,&
                                             timestop
!NB for reducing results of calculations that use dq, which is now spread over nodes
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_ddapc_forces'
  PUBLIC :: ewald_ddapc_force,&
            reset_ch_pulay,&
            evaluate_restraint_functional,&
            restraint_functional_force,&
            solvation_ddapc_force

CONTAINS

! *****************************************************************************
!> \brief Evaluates the Ewald term E2 and E3 energy term for the decoupling/coupling
!>      of periodic images
!> \param qs_env ...
!> \param coeff ...
!> \param apply_qmmm_periodic ...
!> \param factor ...
!> \param multipole_section ...
!> \param cell ...
!> \param particle_set ...
!> \param radii ...
!> \param dq ...
!> \param charges ...
!> \par History
!>      08.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  RECURSIVE SUBROUTINE ewald_ddapc_force(qs_env, coeff, apply_qmmm_periodic,&
       factor, multipole_section, cell, particle_set, radii, dq, charges)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(pw_type), POINTER                   :: coeff
    LOGICAL, INTENT(IN)                      :: apply_qmmm_periodic
    REAL(KIND=dp), INTENT(IN)                :: factor
    TYPE(section_vals_type), POINTER         :: multipole_section
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), DIMENSION(:), POINTER     :: radii
    REAL(KIND=dp), DIMENSION(:, :, :), &
      OPTIONAL, POINTER                      :: dq
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: charges

    CHARACTER(len=*), PARAMETER :: routineN = 'ewald_ddapc_force', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, ip1, ip2, iparticle1, iparticle2, k1, k2, k3, n_rep, &
      nmax1, nmax2, nmax3, r1, r2, r3, rmax1, rmax2, rmax3, stat
    LOGICAL                                  :: analyt, failure
    REAL(KIND=dp)                            :: alpha, eps, fac, fac3, fs, &
                                                galpha, gsq, gsqi, ij_fac, &
                                                q1t, q2t, r, r2tmp, rcut, &
                                                rcut2, t1, t2, tol, tol1
    REAL(KIND=dp), DIMENSION(3)              :: drvec, fvec, gvec, ra, rvec
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: d_el, M
    TYPE(cp_para_env_type), POINTER          :: para_env

    failure = .FALSE.
    NULLIFY(d_el, M, para_env)
    CALL timeset(routineN,handle)
    CALL get_qs_env(qs_env, para_env=para_env)
    CPASSERT(PRESENT(charges))
    CPASSERT(ASSOCIATED(radii))
    CPASSERT(cell%orthorhombic)
    rcut = MIN(cell%hmat(1,1),cell%hmat(2,2),cell%hmat(3,3))/2.0_dp
    CALL section_vals_val_get(multipole_section,"RCUT",n_rep_val=n_rep)
    IF (n_rep==1) CALL section_vals_val_get(multipole_section,"RCUT",r_val=rcut)
    CALL section_vals_val_get(multipole_section,"EWALD_PRECISION",r_val=eps)
    CALL section_vals_val_get(multipole_section,"ANALYTICAL_GTERM",l_val=analyt)
    rcut2=rcut**2
    !
    ! Setting-up parameters for Ewald summation
    !
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
    CPASSERT(stat==0)
    d_el = 0.0_dp
    fac    = 1.e0_dp/cell%deth
    fac3   = fac/8.0_dp
    fvec   = twopi / (/cell%hmat(1,1),cell%hmat(2,2),cell%hmat(3,3)/)
    !
    DO iparticle1 = 1, SIZE(particle_set)
       !NB parallelization
       IF (MOD(iparticle1, para_env%num_pe) /= para_env%mepos) CYCLE
       ip1 = (iparticle1-1)*SIZE(radii)
       q1t   = SUM(charges(ip1+1:ip1+SIZE(radii)))
       DO iparticle2 = 1, iparticle1
          ij_fac = 1.0_dp
          IF (iparticle1==iparticle2) ij_fac=0.5_dp

          ip2 = (iparticle2-1)*SIZE(radii)
          q2t   = SUM(charges(ip2+1:ip2+SIZE(radii)))
          !
          ! Real-Space Contribution
          !
          rvec   = particle_set(iparticle1)%r - particle_set(iparticle2)%r
          IF (iparticle1 /= iparticle2) THEN
             ra = rvec
             r2tmp  = DOT_PRODUCT(ra,ra)
             IF (r2tmp<=rcut2) THEN
                r  = SQRT(r2tmp)
                t1 = erfc(alpha*r) / r
                drvec = ra / r * q1t * q2t * factor
                t2 = -2.0_dp * alpha * EXP(-alpha*alpha*r*r) / (r*rootpi) - t1/r
                d_el(1:3,iparticle1) = d_el(1:3,iparticle1) - t2 * drvec
                d_el(1:3,iparticle2) = d_el(1:3,iparticle2) + t2 * drvec
             END IF
          END IF
          DO r1 = -rmax1, rmax1
             DO r2 = -rmax2, rmax2
                DO r3 = -rmax3, rmax3
                   IF ((r1==0).AND.(r2==0).AND.(r3==0)) CYCLE
                   ra(1) = rvec(1) + cell%hmat(1,1)*r1
                   ra(2) = rvec(2) + cell%hmat(2,2)*r2
                   ra(3) = rvec(3) + cell%hmat(3,3)*r3
                   r2tmp  = DOT_PRODUCT(ra,ra)
                   IF (r2tmp<=rcut2) THEN
                      r  = SQRT(r2tmp)
                      t1 = erfc(alpha*r) / r
                      drvec = ra / r * q1t * q2t * factor * ij_fac
                      t2 = -2.0_dp * alpha * EXP(-alpha*alpha*r*r) / (r*rootpi) - t1/r
                      d_el(1,iparticle1) = d_el(1,iparticle1) - t2 * drvec(1)
                      d_el(2,iparticle1) = d_el(2,iparticle1) - t2 * drvec(2)
                      d_el(3,iparticle1) = d_el(3,iparticle1) - t2 * drvec(3)
                      d_el(1,iparticle2) = d_el(1,iparticle2) + t2 * drvec(1)
                      d_el(2,iparticle2) = d_el(2,iparticle2) + t2 * drvec(2)
                      d_el(3,iparticle2) = d_el(3,iparticle2) + t2 * drvec(3)
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
                      t2 = -SIN(DOT_PRODUCT(gvec,rvec)) * t1 * q1t * q2t * factor * fourpi
                      d_el(1:3,iparticle1) = d_el(1:3,iparticle1) - t2 * gvec
                      d_el(1:3,iparticle2) = d_el(1:3,iparticle2) + t2 * gvec
                   END DO
                END DO
             END DO
          ELSE
             gvec = Eval_d_Interp_Spl3_pbc( rvec, coeff) * q1t * q2t * factor * fourpi
             d_el(1:3,iparticle1) = d_el(1:3,iparticle1) -  gvec
             d_el(1:3,iparticle2) = d_el(1:3,iparticle2) +  gvec
          END IF
          IF (iparticle1/=iparticle2) THEN
             ra = rvec
             r  = SQRT(DOT_PRODUCT(ra,ra))
             t2    = - 1.0_dp /(r*r) * factor
             drvec = ra / r * q1t * q2t
             d_el(1:3,iparticle1) = d_el(1:3,iparticle1) + t2 * drvec
             d_el(1:3,iparticle2) = d_el(1:3,iparticle2) - t2 * drvec
          END IF
       END DO ! iparticle2
    END DO ! iparticle1
    !NB parallelization
    CALL mp_sum(d_el, para_env%group)
    M => qs_env%cp_ddapc_env%Md
    IF (apply_qmmm_periodic) M => qs_env%cp_ddapc_env%Mr
    CALL cp_decpl_ddapc_forces(qs_env, M, charges, dq, d_el, particle_set)
    DEALLOCATE(d_el, stat=stat)
    CPASSERT(stat==0)
    CALL timestop(handle)
  END SUBROUTINE ewald_ddapc_force

! *****************************************************************************
!> \brief Evaluation of the pulay forces due to the fitted charge density
!> \param qs_env ...
!> \param M ...
!> \param charges ...
!> \param dq ...
!> \param d_el ...
!> \param particle_set ...
!> \par History
!>      08.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE cp_decpl_ddapc_forces(qs_env, M, charges, dq, d_el, particle_set)
    TYPE(qs_environment_type), POINTER       :: qs_env
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: M
    REAL(KIND=dp), DIMENSION(:), POINTER     :: charges
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: dq
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: d_el
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

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
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force

    failure = .FALSE.
    NULLIFY(para_env)
    CALL timeset(routineN,handle)
    natom = SIZE(particle_set)
    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    para_env=para_env,&
                    force=force)
    ALLOCATE (atom_of_kind(natom),stat=stat)
    CPASSERT(stat==0)
    ALLOCATE (kind_of(natom),stat=stat)
    CPASSERT(stat==0)
    ALLOCATE (chf(3,natom),stat=stat)
    CPASSERT(stat==0)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             atom_of_kind=atom_of_kind,&
                             kind_of=kind_of)

    ALLOCATE(uv(SIZE(M,1)), stat=stat)
    CPASSERT(stat==0)
    uv(:) = MATMUL(M,charges)
    DO k = 1, natom
       DO j = 1, 3
          chf(j,k) =  DOT_PRODUCT(uv,dq(:,k,j))
       END DO
    END DO
    !NB now that get_ddapc returns dq that's spread over nodes, must reduce chf here
    CALL mp_sum(chf, para_env%group)
    DO iatom=1,natom
       ikind = kind_of(iatom)
       i = atom_of_kind(iatom)
       force(ikind)%ch_pulay(1:3,i) = force(ikind)%ch_pulay(1:3,i) + chf(1:3,iatom) + d_el(1:3,iatom)
    END DO
    DEALLOCATE(atom_of_kind, stat=stat)
    CPASSERT(stat==0)
    DEALLOCATE(kind_of, stat=stat)
    CPASSERT(stat==0)
    DEALLOCATE(chf, stat=stat)
    CPASSERT(stat==0)
    DEALLOCATE(uv, stat=stat)
    CPASSERT(stat==0)
    CALL timestop(handle)
  END SUBROUTINE cp_decpl_ddapc_forces

! *****************************************************************************
!> \brief Evaluation of the pulay forces due to the fitted charge density
!> \param qs_env ...
!> \par History
!>      08.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE reset_ch_pulay(qs_env)
    TYPE(qs_environment_type), POINTER       :: qs_env

    CHARACTER(len=*), PARAMETER :: routineN = 'reset_ch_pulay', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ind
    LOGICAL                                  :: failure
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force

    failure = .FALSE.

    CALL timeset(routineN,handle)
    CALL get_qs_env(qs_env=qs_env,&
                    force=force)
    DO ind = 1, SIZE(force)
       force(ind)%ch_pulay = 0.0_dp
    END DO
    CALL timestop(handle)
  END SUBROUTINE reset_ch_pulay

! *****************************************************************************
!> \brief computes energy and derivatives given a set of charges
!> \param ddapc_restraint_control ...
!> \param n_gauss ...
!> \param uv derivate of energy wrt the corresponding charge entry
!> \param charges current value of the charges (one number for each gaussian used)
!>
!> \param energy_res energy due to the restraint
!> \par History
!>      02.2006 [Joost VandeVondele]
!>               modified [Teo]
!> \note
!>       should be easy to adapt for other specialized cases
! *****************************************************************************
  SUBROUTINE evaluate_restraint_functional(ddapc_restraint_control, n_gauss, uv,&
       charges, energy_res)
    TYPE(ddapc_restraint_type), &
      INTENT(INOUT)                          :: ddapc_restraint_control
    INTEGER, INTENT(in)                      :: n_gauss
    REAL(KIND=dp), DIMENSION(:)              :: uv
    REAL(KIND=dp), DIMENSION(:), POINTER     :: charges
    REAL(KIND=dp), INTENT(INOUT)             :: energy_res

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'evaluate_restraint_functional', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: I, ind
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dE, order_p

! order parameter (i.e. the sum of the charges of the selected atoms)

    failure = .FALSE.
    order_p = 0.0_dp
    DO I=1,ddapc_restraint_control%natoms
       ind     = (ddapc_restraint_control%atoms(I)-1)*n_gauss
       order_p = order_p + ddapc_restraint_control%coeff(I)*SUM(charges(ind+1:ind+n_gauss))
    ENDDO
    ddapc_restraint_control%ddapc_order_p = order_p

    SELECT CASE(ddapc_restraint_control%functional_form)
    CASE(do_ddapc_restraint)
      ! the restraint energy
      energy_res = ddapc_restraint_control%strength*(order_p-ddapc_restraint_control%target)**2.0_dp

      ! derivative of the energy
      dE=2.0_dp*ddapc_restraint_control%strength*(order_p-ddapc_restraint_control%target)
      DO I=1,ddapc_restraint_control%natoms
         ind     = (ddapc_restraint_control%atoms(I)-1)*n_gauss
         uv(ind+1:ind+n_gauss)  = dE*ddapc_restraint_control%coeff(I)
      ENDDO
    CASE(do_ddapc_constraint)
      energy_res = ddapc_restraint_control%strength*(order_p-ddapc_restraint_control%target)

      ! derivative of the energy
      DO I=1,ddapc_restraint_control%natoms
         ind     = (ddapc_restraint_control%atoms(I)-1)*n_gauss
         uv(ind+1:ind+n_gauss)  = ddapc_restraint_control%strength*ddapc_restraint_control%coeff(I)
      ENDDO

    CASE DEFAULT
       CPASSERT(.FALSE.)
    END SELECT

  END SUBROUTINE evaluate_restraint_functional

! *****************************************************************************
!> \brief computes derivatives for DDAPC restraint
!> \param qs_env ...
!> \param ddapc_restraint_control ...
!> \param dq ...
!> \param charges ...
!> \param n_gauss ...
!> \param particle_set ...
!> \par History
!>      02.2006 [Joost VandeVondele]
!>              modified [Teo]
!> \note
!>       should be easy to adapt for other specialized cases
! *****************************************************************************
  SUBROUTINE restraint_functional_force(qs_env, ddapc_restraint_control, dq, charges,&
       n_gauss, particle_set)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(ddapc_restraint_type), &
      INTENT(INOUT)                          :: ddapc_restraint_control
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: dq
    REAL(KIND=dp), DIMENSION(:), POINTER     :: charges
    INTEGER, INTENT(in)                      :: n_gauss
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    CHARACTER(len=*), PARAMETER :: routineN = 'restraint_functional_force', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iatom, ikind, j, &
                                                k, natom, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind, kind_of
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: dum
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: uv
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: chf
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force

    failure = .FALSE.
    NULLIFY(para_env)
    CALL timeset(routineN,handle)
    natom = SIZE(particle_set)
    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    para_env=para_env,&
                    force=force)
    ALLOCATE (atom_of_kind(natom),stat=stat)
    CPASSERT(stat==0)
    ALLOCATE (kind_of(natom),stat=stat)
    CPASSERT(stat==0)
    ALLOCATE (chf(3,natom),stat=stat)
    CPASSERT(stat==0)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             atom_of_kind=atom_of_kind,&
                             kind_of=kind_of)

    ALLOCATE(uv(SIZE(dq,1)), stat=stat)
    CPASSERT(stat==0)
    uv = 0.0_dp
    CALL evaluate_restraint_functional(ddapc_restraint_control, n_gauss, uv,&
         charges, dum)
    DO k = 1, natom
       DO j = 1, 3
          chf(j,k) =  DOT_PRODUCT(uv,dq(:,k,j))
       END DO
    END DO
    !NB now that get_ddapc returns dq that's spread over nodes, must reduce chf here
    CALL mp_sum(chf, para_env%group)
    DO iatom=1,natom
       ikind = kind_of(iatom)
       i = atom_of_kind(iatom)
       force(ikind)%ch_pulay(1:3,i) = force(ikind)%ch_pulay(1:3,i) + chf(1:3,iatom)
    END DO
    DEALLOCATE(atom_of_kind, stat=stat)
    CPASSERT(stat==0)
    DEALLOCATE(kind_of, stat=stat)
    CPASSERT(stat==0)
    DEALLOCATE(chf, stat=stat)
    CPASSERT(stat==0)
    DEALLOCATE(uv, stat=stat)
    CPASSERT(stat==0)
    CALL timestop(handle)

  END SUBROUTINE restraint_functional_force

! *****************************************************************************
!> \brief Evaluates the electrostatic potential due to a simple solvation model
!>      Spherical cavity in a dieletric medium
!> \param qs_env ...
!> \param solvation_section ...
!> \param particle_set ...
!> \param radii ...
!> \param dq ...
!> \param charges ...
!> \par History
!>      08.2006 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE solvation_ddapc_force(qs_env, solvation_section, particle_set,&
       radii, dq, charges)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type), POINTER         :: solvation_section
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), DIMENSION(:), POINTER     :: radii
    REAL(KIND=dp), DIMENSION(:, :, :), &
      OPTIONAL, POINTER                      :: dq
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: charges

    CHARACTER(len=*), PARAMETER :: routineN = 'solvation_ddapc_force', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ip1, ip2, iparticle1, &
                                                iparticle2, l, lmax, n_rep1, &
                                                n_rep2, stat, weight
    INTEGER, DIMENSION(:), POINTER           :: list
    LOGICAL                                  :: failure, fixed_center
    REAL(KIND=dp) :: center(3), dcos1(3), dcos2(3), dpos1(3), dpos2(3), &
      eps_in, eps_out, factor1(3), factor2(3), lr, mass, mycos, pos1, pos1i, &
      pos2, pos2i, ptcos, q1t, q2t, r1(3), r1s, r2(3), r2s, Rs, rvec(3)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: pos, R0
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: d_el, LocP, M

    failure = .FALSE.
    fixed_center = .FALSE.
    NULLIFY(d_el, M)
    eps_in = 1.0_dp
    CALL section_vals_val_get(solvation_section,"EPS_OUT",r_val=eps_out)
    CALL section_vals_val_get(solvation_section,"LMAX",i_val=lmax)
    CALL section_vals_val_get(solvation_section,"SPHERE%RADIUS",r_val=Rs)
    CALL section_vals_val_get(solvation_section,"SPHERE%CENTER%XYZ",n_rep_val=n_rep1)
    IF (n_rep1/=0) THEN
       CALL section_vals_val_get(solvation_section,"SPHERE%CENTER%XYZ",r_vals=R0)
       center = R0
    ELSE
       CALL section_vals_val_get(solvation_section,"SPHERE%CENTER%ATOM_LIST",&
            n_rep_val=n_rep2)
       IF (n_rep2/=0) THEN
          CALL section_vals_val_get(solvation_section,"SPHERE%CENTER%ATOM_LIST",i_vals=list)
          CALL section_vals_val_get(solvation_section,"SPHERE%CENTER%WEIGHT_TYPE",i_val=weight)
          ALLOCATE(R0(SIZE(list)),stat=stat)
          CPASSERT(stat==0)
          SELECT CASE(weight)
          CASE (weight_type_unit)
             R0 = 0.0_dp
             DO i = 1, SIZE(list)
                R0 = R0 + particle_set(list(i))%r
             END DO
             R0 = R0 / REAL(SIZE(list),KIND=dp)
          CASE (weight_type_mass)
             R0 = 0.0_dp
             mass = 0.0_dp
             DO i = 1, SIZE(list)
                R0 = R0 + particle_set(list(i))%r * particle_set(list(i))%atomic_kind%mass
                mass = mass + particle_set(list(i))%atomic_kind%mass
             END DO
             R0 = R0 / mass
          END SELECT
          center = R0
          DEALLOCATE(R0, stat=stat)
          CPASSERT(stat==0)
       END IF
    END IF
    CPASSERT(n_rep1/=0.OR.n_rep2/=0)
    ! Potential calculation
    ALLOCATE(LocP(0:lmax,SIZE(particle_set)),stat=stat)
    CPASSERT(stat==0)
    ALLOCATE(pos(SIZE(particle_set)),stat=stat)
    CPASSERT(stat==0)
    ALLOCATE(d_el(3,SIZE(particle_set)), stat=stat)
    CPASSERT(stat==0)
    d_el = 0.0_dp
    ! Determining the single atomic contribution to the dielectric dipole
    DO i = 1, SIZE(particle_set)
       rvec = particle_set(i)%r-center
       r2s  = DOT_PRODUCT(rvec,rvec)
       r1s  = SQRT(r2s)
       LocP(:,i) = 0.0_dp
       IF (r1s/=0.0_dp) THEN
          DO l = 0, lmax
             LocP(l,i) = ( r1s**l * REAL(l+1,KIND=dp) * (eps_in-eps_out))/&
                  (Rs**(2*l+1)*eps_in*(REAL(l,KIND=dp)*eps_in+REAL(l+1,KIND=dp)*eps_out))
          END DO
       END IF
       pos(i) =  r1s
    END DO
    ! Computes the full derivatives of the interaction energy
    DO iparticle1 = 1, SIZE(particle_set)
       ip1 = (iparticle1-1)*SIZE(radii)
       q1t = SUM(charges(ip1+1:ip1+SIZE(radii)))
       DO iparticle2 = 1, iparticle1
          ip2 = (iparticle2-1)*SIZE(radii)
          q2t = SUM(charges(ip2+1:ip2+SIZE(radii)))
          !
          r1 = particle_set(iparticle1)%r-center
          r2 = particle_set(iparticle2)%r-center
          pos1  = pos(iparticle1)
          pos2  = pos(iparticle2)
          factor1 = 0.0_dp
          factor2 = 0.0_dp
          IF (pos1*pos2/=0.0_dp) THEN
             pos1i = 1.0_dp/pos1
             pos2i = 1.0_dp/pos1
             dpos1 = pos1i*r1
             dpos2 = pos2i*r2
             ptcos = DOT_PRODUCT(r1,r2)
             mycos = ptcos/(pos1*pos2)
             IF (ABS(mycos) > 1.0_dp) mycos = SIGN(1.0_dp, mycos)
             dcos1 = (r2*(pos1*pos2)-pos2*dpos1*ptcos)/(pos1*pos2)**2
             dcos2 = (r1*(pos1*pos2)-pos1*dpos2*ptcos)/(pos1*pos2)**2

             DO l = 1, lmax
                lr =  REAL(l,KIND=dp)
                factor1 = factor1 + lr*LocP(l,iparticle2)*pos1**(l-1)*legendre(mycos,l,0)*dpos1&
                     + LocP(l,iparticle2)*pos1**l*dlegendre(mycos,l,0)*dcos1
                factor2 = factor2 + lr*LocP(l,iparticle1)*pos2**(l-1)*legendre(mycos,l,0)*dpos2&
                     + LocP(l,iparticle1)*pos2**l*dlegendre(mycos,l,0)*dcos2
             END DO
          END IF
          factor1 = factor1 * q1t * q2t
          factor2 = factor2 * q1t * q2t
          d_el(1:3,iparticle1) = d_el(1:3,iparticle1) + 0.5_dp * factor1
          d_el(1:3,iparticle2) = d_el(1:3,iparticle2) + 0.5_dp * factor2
       END DO
    END DO
    DEALLOCATE(pos,stat=stat)
    CPASSERT(stat==0)
    DEALLOCATE(LocP,stat=stat)
    CPASSERT(stat==0)
    M => qs_env%cp_ddapc_env%Ms
    CALL cp_decpl_ddapc_forces(qs_env, M, charges, dq, d_el, particle_set)
    DEALLOCATE(d_el, stat=stat)
    CPASSERT(stat==0)
  END SUBROUTINE solvation_ddapc_force

END MODULE cp_ddapc_forces
