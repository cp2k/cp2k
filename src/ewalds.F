!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      JGH (15-Mar-2001) : New routine ewald_setup (former pme_setup)
!>      JGH (23-Mar-2001) : Get rid of global variable ewald_grp
! *****************************************************************************
MODULE ewalds

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE bibliography,                    ONLY: Ewald1921,&
                                             cite_reference
  USE cell_types,                      ONLY: cell_type
  USE dg_rho0_types,                   ONLY: dg_rho0_type
  USE dg_types,                        ONLY: dg_get,&
                                             dg_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE ewald_environment_types,         ONLY: ewald_env_get,&
                                             ewald_environment_type
  USE ewald_pw_types,                  ONLY: ewald_pw_get,&
                                             ewald_pw_type
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: fourpi,&
                                             oorootpi,&
                                             pi
  USE mathlib,                         ONLY: matvec_3x3
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_poisson_types,                ONLY: do_ewald_none
  USE pw_pool_types,                   ONLY: pw_pool_type
  USE shell_potential_types,           ONLY: get_shell,&
                                             shell_kind_type
  USE structure_factor_types,          ONLY: structure_factor_type
  USE structure_factors,               ONLY: structure_factor_allocate,&
                                             structure_factor_deallocate,&
                                             structure_factor_evaluate
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ewalds'

  PRIVATE
  PUBLIC :: ewald_evaluate, ewald_self, ewald_self_atom, ewald_print

CONTAINS

! *****************************************************************************
!> \brief computes the potential and the force from the g-space part of
!>      the 1/r potential
!>      Ref.: J.-P. Hansen, Enrico Fermi School, 1985
!>      Note: Only the positive G-vectors are used in the sum.
!> \param ewald_env ...
!> \param ewald_pw ...
!> \param cell ...
!> \param atomic_kind_set ...
!> \param particle_set ...
!> \param local_particles ...
!> \param fg_coulomb ...
!> \param vg_coulomb ...
!> \param pv_g ...
!> \param use_virial ...
!> \param charges ...
!> \param e_coulomb ...
!> \param pv_coulomb ...
!> \param error ...
!> \par History
!>      JGH (21-Feb-2001) : changed name
!> \author CJM
! *****************************************************************************
  SUBROUTINE ewald_evaluate (ewald_env, ewald_pw, cell, atomic_kind_set, particle_set,&
       local_particles, fg_coulomb, vg_coulomb, pv_g, use_virial, charges, e_coulomb,&
       pv_coulomb, error )
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(cell_type), POINTER                 :: cell
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_particles
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: fg_coulomb
    REAL(KIND=dp), INTENT(OUT)               :: vg_coulomb
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: pv_g
    LOGICAL, INTENT(IN)                      :: use_virial
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: charges, e_coulomb
    REAL(KIND=dp), DIMENSION(:, :, :), &
      OPTIONAL, POINTER                      :: pv_coulomb
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ewald_evaluate', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp)                         :: snode
    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:)                           :: summe
    INTEGER :: gpt, group, handle, iparticle, iparticle_kind, &
      iparticle_local, isos, lp, mp, nnodes, node, np, nparticle_kind, &
      nparticle_local
    INTEGER, DIMENSION(:, :), POINTER        :: bds
    LOGICAL                                  :: atenergy, atstress, failure, &
                                                use_charge_array
    REAL(KIND=dp)                            :: alpha, denom, e_igdotr, &
                                                factor, four_alpha_sq, gauss, &
                                                pref, q
    REAL(KIND=dp), DIMENSION(3)              :: vec
    REAL(KIND=dp), DIMENSION(:), POINTER     :: charge
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: rho0
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(dg_rho0_type), POINTER              :: dg_rho0
    TYPE(dg_type), POINTER                   :: dg
    TYPE(pw_grid_type), POINTER              :: pw_grid
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(structure_factor_type)              :: exp_igr

    CALL timeset(routineN,handle)
    CALL cite_reference(Ewald1921)
    use_charge_array=PRESENT(charges)
    IF (use_charge_array) use_charge_array=ASSOCIATED(charges)
    atenergy=PRESENT(e_coulomb)
    IF (atenergy) atenergy=ASSOCIATED(e_coulomb)
    IF (atenergy) e_coulomb=0._dp

    atstress=PRESENT(pv_coulomb)
    IF (atstress) atstress=ASSOCIATED(pv_coulomb)
    IF (atstress) pv_coulomb=0._dp

    ! pointing
    CALL ewald_env_get (ewald_env, alpha=alpha, group = group ,error=error)
    CALL ewald_pw_get (ewald_pw, pw_big_pool=pw_pool, dg = dg )
    CALL dg_get ( dg, dg_rho0=dg_rho0 )
    rho0    => dg_rho0 % density % pw % cr3d
    pw_grid => pw_pool % pw_grid
    bds     =>  pw_grid % bounds

    ! allocating
    nparticle_kind = SIZE ( atomic_kind_set )
    nnodes = 0
    DO iparticle_kind = 1, nparticle_kind
       nnodes = nnodes + local_particles%n_el(iparticle_kind)
    ENDDO

    CALL structure_factor_allocate ( pw_grid % bounds, nnodes, exp_igr)

    ALLOCATE (summe(1:pw_grid%ngpts_cut),STAT=isos)
    CPPostcondition(isos==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (charge(1:nnodes),STAT=isos)
    CPPostcondition(isos==0,cp_failure_level,routineP,error,failure)

    ! Initializing vg_coulomb and fg_coulomb
    vg_coulomb = 0.0_dp
    fg_coulomb = 0.0_dp
    IF (use_virial) pv_g = 0.0_dp
    ! defining four_alpha_sq
    four_alpha_sq = 4.0_dp * alpha ** 2
    ! zero node count
    node = 0
    DO iparticle_kind = 1, nparticle_kind
       nparticle_local = local_particles%n_el(iparticle_kind)
       IF (use_charge_array) THEN
          DO iparticle_local=1,nparticle_local
             node = node + 1
             iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
             charge(node) = charges(iparticle)
             CALL matvec_3x3 (vec, cell % h_inv, particle_set ( iparticle ) % r )
             CALL structure_factor_evaluate ( vec, pw_grid % npts, exp_igr % lb, &
                  exp_igr % ex(:,node), exp_igr % ey(:,node), exp_igr % ez(:,node) )
          END DO
       ELSE
          atomic_kind => atomic_kind_set(iparticle_kind)
          CALL get_atomic_kind(atomic_kind=atomic_kind, qeff=q)
          DO iparticle_local=1,nparticle_local
             node = node + 1
             iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
             charge(node) = q
             CALL matvec_3x3 (vec, cell % h_inv, particle_set ( iparticle ) % r )
             CALL structure_factor_evaluate ( vec, pw_grid % npts, exp_igr % lb, &
                  exp_igr % ex(:,node), exp_igr % ey(:,node), exp_igr % ez(:,node) )
          END DO
       END IF
    END DO

    summe ( : ) = CMPLX ( 0.0_dp, 0.0_dp,KIND=dp)
    ! looping over the positive g-vectors
    DO gpt = 1, pw_grid % ngpts_cut_local

       lp = pw_grid % mapl % pos ( pw_grid % g_hat ( 1, gpt ) )
       mp = pw_grid % mapm % pos ( pw_grid % g_hat ( 2, gpt ) )
       np = pw_grid % mapn % pos ( pw_grid % g_hat ( 3, gpt ) )

       lp = lp + bds ( 1, 1 )
       mp = mp + bds ( 1, 2 )
       np = np + bds ( 1, 3 )

       ! initializing sum to be used in the energy and force
       DO node = 1, nnodes
          summe ( gpt ) = summe ( gpt ) + charge ( node ) * &
               ( exp_igr % ex ( lp, node ) &
               * exp_igr % ey ( mp, node ) &
               * exp_igr % ez ( np, node ) )
       END DO
    END DO
    CALL mp_sum ( summe, group )

    pref = fourpi / pw_grid % vol

    ! looping over the positive g-vectors
    DO gpt = 1, pw_grid % ngpts_cut_local
       ! computing the potential energy
       lp = pw_grid % mapl % pos ( pw_grid % g_hat ( 1, gpt ) )
       mp = pw_grid % mapm % pos ( pw_grid % g_hat ( 2, gpt ) )
       np = pw_grid % mapn % pos ( pw_grid % g_hat ( 3, gpt ) )

       lp = lp + bds ( 1, 1 )
       mp = mp + bds ( 1, 2 )
       np = np + bds ( 1, 3 )

       IF ( pw_grid % gsq ( gpt ) <= 1.0E-10_dp ) CYCLE

       gauss = ( rho0 ( lp, mp, np ) * pw_grid % vol ) ** 2 / pw_grid % gsq ( gpt )
       factor = gauss * REAL ( summe ( gpt ) * CONJG ( summe ( gpt ) ), KIND=dp )
       vg_coulomb = vg_coulomb + factor

       ! atomic energies
       IF (atenergy) THEN
          DO node = 1, nnodes
             snode = CONJG ( exp_igr % ex ( lp, node ) &
                           * exp_igr % ey ( mp, node ) &
                           * exp_igr % ez ( np, node ) )
             e_coulomb(node) = e_coulomb(node) + gauss*charge(node)*REAL(summe(gpt)*snode, KIND=dp)
          END DO
       END IF

       ! computing the force
       node = 0
       DO node = 1, nnodes
          e_igdotr = AIMAG ( summe(gpt) * CONJG &
               ( exp_igr % ex ( lp, node ) &
               * exp_igr % ey ( mp, node ) &
               * exp_igr % ez ( np, node ) ) )
          fg_coulomb ( :, node ) = fg_coulomb ( :, node ) &
               + charge ( node ) * gauss * e_igdotr * pw_grid % g ( :, gpt )
       END DO

       ! compute the virial P*V
       denom = 1.0_dp / four_alpha_sq + 1.0_dp / pw_grid % gsq ( gpt )
       IF (use_virial) THEN
          pv_g(1,1) = pv_g(1,1) + factor*(1.0_dp - 2.0_dp*pw_grid%g(1,gpt)*pw_grid%g(1,gpt)*denom)
          pv_g(1,2) = pv_g(1,2) - factor*(         2.0_dp*pw_grid%g(1,gpt)*pw_grid%g(2,gpt)*denom)
          pv_g(1,3) = pv_g(1,3) - factor*(         2.0_dp*pw_grid%g(1,gpt)*pw_grid%g(3,gpt)*denom)
          pv_g(2,1) = pv_g(2,1) - factor*(         2.0_dp*pw_grid%g(2,gpt)*pw_grid%g(1,gpt)*denom)
          pv_g(2,2) = pv_g(2,2) + factor*(1.0_dp - 2.0_dp*pw_grid%g(2,gpt)*pw_grid%g(2,gpt)*denom)
          pv_g(2,3) = pv_g(2,3) - factor*(         2.0_dp*pw_grid%g(2,gpt)*pw_grid%g(3,gpt)*denom)
          pv_g(3,1) = pv_g(3,1) - factor*(         2.0_dp*pw_grid%g(3,gpt)*pw_grid%g(1,gpt)*denom)
          pv_g(3,2) = pv_g(3,2) - factor*(         2.0_dp*pw_grid%g(3,gpt)*pw_grid%g(2,gpt)*denom)
          pv_g(3,3) = pv_g(3,3) + factor*(1.0_dp - 2.0_dp*pw_grid%g(3,gpt)*pw_grid%g(3,gpt)*denom)
       END IF
       IF (atstress) THEN
          DO node = 1, nnodes
             snode = CONJG ( exp_igr % ex ( lp, node ) &
                           * exp_igr % ey ( mp, node ) &
                           * exp_igr % ez ( np, node ) )
             factor = gauss * charge ( node ) * REAL ( summe(gpt) * snode, KIND=dp )
             pv_coulomb(1,1,node) = pv_coulomb(1,1,node) + factor*(1.0_dp - 2.0_dp*pw_grid%g(1,gpt)*pw_grid%g(1,gpt)*denom)
             pv_coulomb(1,2,node) = pv_coulomb(1,2,node) - factor*(         2.0_dp*pw_grid%g(1,gpt)*pw_grid%g(2,gpt)*denom)
             pv_coulomb(1,3,node) = pv_coulomb(1,3,node) - factor*(         2.0_dp*pw_grid%g(1,gpt)*pw_grid%g(3,gpt)*denom)
             pv_coulomb(2,1,node) = pv_coulomb(2,1,node) - factor*(         2.0_dp*pw_grid%g(2,gpt)*pw_grid%g(1,gpt)*denom)
             pv_coulomb(2,2,node) = pv_coulomb(2,2,node) + factor*(1.0_dp - 2.0_dp*pw_grid%g(2,gpt)*pw_grid%g(2,gpt)*denom)
             pv_coulomb(2,3,node) = pv_coulomb(2,3,node) - factor*(         2.0_dp*pw_grid%g(2,gpt)*pw_grid%g(3,gpt)*denom)
             pv_coulomb(3,1,node) = pv_coulomb(3,1,node) - factor*(         2.0_dp*pw_grid%g(3,gpt)*pw_grid%g(1,gpt)*denom)
             pv_coulomb(3,2,node) = pv_coulomb(3,2,node) - factor*(         2.0_dp*pw_grid%g(3,gpt)*pw_grid%g(2,gpt)*denom)
             pv_coulomb(3,3,node) = pv_coulomb(3,3,node) + factor*(1.0_dp - 2.0_dp*pw_grid%g(3,gpt)*pw_grid%g(3,gpt)*denom)
          END DO
       END IF
    END DO

    vg_coulomb = vg_coulomb * pref
    IF (use_virial) pv_g = pv_g * pref
    IF (atenergy) e_coulomb = e_coulomb * pref
    IF (atstress) pv_coulomb = pv_coulomb * pref

    fg_coulomb = fg_coulomb * (2.0_dp*pref)

    CALL structure_factor_deallocate ( exp_igr )

    DEALLOCATE ( charge, summe, STAT = isos )
    CPPostcondition(isos==0,cp_failure_level,routineP,error,failure)
    NULLIFY ( charge )

    CALL timestop(handle)

  END SUBROUTINE ewald_evaluate

! *****************************************************************************
!> \brief Computes the self interaction from g-space
!>      and the neutralizing background
!> \param ewald_env ...
!> \param cell ...
!> \param atomic_kind_set ...
!> \param local_particles ...
!> \param e_self ...
!> \param e_neut ...
!> \param charges ...
!> \param error ...
!> \par History
!>      none
!> \author CJM
! *****************************************************************************
  SUBROUTINE ewald_self ( ewald_env, cell, atomic_kind_set, local_particles, e_self,&
       e_neut, charges, error)

    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(cell_type), POINTER                 :: cell
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_particles
    REAL(KIND=dp), INTENT(OUT)               :: e_self, e_neut
    REAL(KIND=dp), DIMENSION(:), POINTER     :: charges
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ewald_self', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ewald_type, group, &
                                                iparticle_kind, &
                                                nparticle_kind, &
                                                nparticle_local
    LOGICAL                                  :: is_shell
    REAL(KIND=dp)                            :: alpha, mm_radius, q, q_neutg, &
                                                q_self, q_sum, qcore, qshell
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(shell_kind_type), POINTER           :: shell

    CALL ewald_env_get ( ewald_env, ewald_type = ewald_type, &
         alpha = alpha, group = group ,error=error)
    q_neutg = 0.0_dp
    q_self = 0.0_dp
    q_sum = 0.0_dp
    nparticle_kind = SIZE ( atomic_kind_set )
    IF (ASSOCIATED(charges)) THEN
       q_self = DOT_PRODUCT(charges,charges)
       q_sum  = SUM(charges)
       ! check and abort..
       DO iparticle_kind=1,nparticle_kind
          atomic_kind => atomic_kind_set(iparticle_kind)
          CALL get_atomic_kind(atomic_kind=atomic_kind, mm_radius=mm_radius)
          IF(mm_radius>0.0_dp)THEN
             CALL cp_unimplemented_error(fromWhere=routineP, &
                       message="Array of charges not implemented for mm_radius>0.0 !!",&
                       error=error, error_level=cp_failure_level)
          END IF
       END DO
    ELSE
       DO iparticle_kind=1,nparticle_kind
          atomic_kind => atomic_kind_set(iparticle_kind)
          CALL get_atomic_kind(atomic_kind=atomic_kind, mm_radius=mm_radius, &
               qeff=q, shell_active=is_shell, shell=shell)
          nparticle_local = local_particles%n_el(iparticle_kind)
          IF (is_shell) THEN
             CALL get_shell(shell=shell, charge_core=qcore, charge_shell=qshell, error=error)
             ! MI: the core-shell ES interaction, when core and shell belong to the same ion, is excluded
             !     in the nonbond correction term. Therefore, here the self interaction is computed entirely
             q_self = q_self + qcore*qcore*nparticle_local + qshell*qshell*nparticle_local
             q_sum  = q_sum + qcore*nparticle_local + qshell*nparticle_local
             IF (mm_radius > 0) THEN
                ! the core is always a point charge
                q_neutg = q_neutg + 2.0_dp*qshell*mm_radius**2
             END IF
          ELSE
             q_self = q_self + q*q*nparticle_local
             q_sum  = q_sum + q*nparticle_local
             IF (mm_radius > 0) THEN
                q_neutg = q_neutg + 2.0_dp*q*mm_radius**2
             END IF
          END IF
       END DO

       CALL mp_sum ( q_self, group )
       CALL mp_sum ( q_sum, group )
    END IF

    e_neut = 0.0_dp
    e_self = 0.0_dp
    IF ( ewald_type /= do_ewald_none) THEN
       e_self = -q_self *  alpha * oorootpi
       e_neut = -q_sum*pi/(2.0_dp*cell%deth)*(q_sum/alpha**2 - q_neutg)
    END IF

  END SUBROUTINE ewald_self

! *****************************************************************************
!> \brief Computes the self interaction per atom
!> \param ewald_env ...
!> \param cell ...
!> \param atomic_kind_set ...
!> \param local_particles ...
!> \param e_self ...
!> \param charges ...
!> \param error ...
!> \par History
!>      none
!> \author JHU from ewald_self
! *****************************************************************************
  SUBROUTINE ewald_self_atom ( ewald_env, cell, atomic_kind_set, local_particles, e_self,&
       charges, error)

    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(cell_type), POINTER                 :: cell
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_particles
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: e_self(:)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: charges
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ewald_self_atom', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ewald_type, ii, &
                                                iparticle_kind, &
                                                iparticle_local, &
                                                nparticle_kind, &
                                                nparticle_local
    LOGICAL                                  :: is_shell
    REAL(KIND=dp)                            :: alpha, fself, q, qcore, qshell
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(shell_kind_type), POINTER           :: shell

    CALL ewald_env_get(ewald_env,ewald_type=ewald_type,alpha=alpha,error=error)

    fself = alpha*oorootpi

    IF (ewald_type /= do_ewald_none) THEN
       nparticle_kind = SIZE(atomic_kind_set)
       IF (ASSOCIATED(charges)) THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="Atomic energy not implemented for charges",&
               error=error, error_level=cp_failure_level)
       ELSE
          DO iparticle_kind=1,nparticle_kind
             atomic_kind => atomic_kind_set(iparticle_kind)
             nparticle_local = local_particles%n_el(iparticle_kind)
             CALL get_atomic_kind(atomic_kind=atomic_kind,qeff=q,&
                                  shell_active=is_shell,shell=shell)
             IF (is_shell) THEN
                CALL get_shell(shell=shell,charge_core=qcore,charge_shell=qshell,error=error)
                DO iparticle_local=1,nparticle_local
                   ii = local_particles%list(iparticle_kind)%array(iparticle_local)
                   e_self(ii) = e_self(ii) - (qcore*qcore + qshell*qshell)*fself
                END DO
             ELSE
                DO iparticle_local=1,nparticle_local
                   ii = local_particles%list(iparticle_kind)%array(iparticle_local)
                   e_self(ii) = e_self(ii) - q*q*fself
                END DO
             END IF
          END DO
       END IF
    END IF

  END SUBROUTINE ewald_self_atom

! *****************************************************************************
!> \brief ...
!> \param iw ...
!> \param pot_nonbond ...
!> \param e_gspace ...
!> \param e_self ...
!> \param e_neut ...
!> \param e_bonded ...
!> \par History
!>      none
!> \author CJM
! *****************************************************************************
  SUBROUTINE ewald_print ( iw, pot_nonbond, e_gspace, e_self, e_neut, e_bonded )

    INTEGER, INTENT(IN)                      :: iw
    REAL(KIND=dp), INTENT(IN)                :: pot_nonbond, e_gspace, &
                                                e_self, e_neut, e_bonded

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ewald_print', &
      routineP = moduleN//':'//routineN

    IF (iw>0) THEN
       WRITE ( iw, '( A, A )' ) ' *********************************', &
            '**********************************************'
       WRITE ( iw, '( A, A, T35, A, T56, E25.15 )' ) ' INITIAL GSPACE ENERGY', &
            '[hartree]', '= ', e_gspace
       WRITE ( iw, '( A, A, T35, A, T56, E25.15 )' ) ' INITIAL NONBONDED ENERGY', &
            '[hartree]', '= ', pot_nonbond
       WRITE ( iw, '( A, A, T35, A, T56, E25.15 )' ) ' SELF ENERGY CORRECTION', &
            '[hartree]', '= ', e_self
       WRITE ( iw, '( A, A, T35, A, T56, E25.15 )' ) ' NEUT. BACKGROUND',&
            '[hartree]', '= ', e_neut
       WRITE ( iw, '( A, A, T35, A, T56, E25.15 )' ) ' BONDED CORRECTION',&
            '[hartree]','= ', e_bonded
       WRITE ( iw, '( A, A )' ) ' *********************************', &
            '**********************************************'
    END IF
  END SUBROUTINE ewald_print

END MODULE ewalds
