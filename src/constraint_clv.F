!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****** cp2k/constraint_clv [1.0] *
!!
!!   NAME
!!     constraint_clv
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
MODULE constraint_clv

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_p_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE kinds,                           ONLY: dp
  USE md,                              ONLY: simulation_parameters_type
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE molecule_kind_types,             ONLY: dist_constraint_type,&
                                             fixd_constraint_type,&
                                             get_molecule_kind,&
                                             molecule_kind_type
  USE molecule_types_new,              ONLY: get_molecule,&
                                             local_dist_constraint_type,&
                                             molecule_type
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: dotprod_3d,&
                                             matvec_3x3
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: shake,&
            rattle,&
            shake_roll,&
            rattle_roll,&
            fix_atom_control

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'constraint_clv'

CONTAINS
!!****** constraint_clv/fix_atom_control [1.0] *
!!
!!   NAME
!!     fix_atom_control
!!
!!   FUNCTION
!!     allows for fix atom constraints
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE fix_atom_control( force_env, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fix_atom_control', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ikind, iparticle, &
                                                iparticle_local, isubsys, &
                                                nfixed_atoms, nkind, &
                                                nparticle_local
    LOGICAL                                  :: failure
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    failure = .FALSE.
    CALL timeset(routineN,'I','',handle)
    IF (.NOT. failure) THEN
       CALL force_env_get(force_env=force_env, subsys=subsys)

       DO isubsys = 1, SIZE(subsys)
          CALL cp_subsys_get(subsys=subsys(isubsys)%subsys,&
               particles=particles, &
               local_particles=local_particles,&
               molecule_kinds_new=molecule_kinds )

          nkind             =  molecule_kinds % n_els
          molecule_kind_set => molecule_kinds % els
          particle_set      => particles%els

          Mol_kind:  DO ikind = 1, nkind
             molecule_kind => molecule_kind_set(ikind)
             CALL get_molecule_kind ( molecule_kind, nfixd=nfixed_atoms, fixd_list=fixd_list )

             IF (nfixed_atoms ==0) CYCLE Mol_Kind
             DO i = 1, SIZE(local_particles%n_el)
                nparticle_local   = local_particles%n_el(i)
                DO iparticle_local=1,nparticle_local
                   iparticle = local_particles%list(i)%array(iparticle_local)
                   IF (ANY(fixd_list(:)%fixd == iparticle)) particle_set(iparticle)%f(:) = 0.0_dp
                END DO
             END DO
          END DO Mol_kind

       END DO
    END IF
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE fix_atom_control

!!****** constraint_clv/shake [1.0] *
!!
!!   NAME
!!     shake
!!
!!   FUNCTION
!!     shake algorithm for collective variables constraints
!!     updates the multiplier one molecule type at a time
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE shake( molecule, particle_set, pos, vel, simpar, max_sigma, ishake )

    TYPE(molecule_type), POINTER             :: molecule
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: pos( :, : ), vel( :, : )
    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar
    REAL(KIND=dp), INTENT(INOUT)             :: max_sigma
    INTEGER, INTENT(INOUT)                   :: ishake

    INTEGER                                  :: first_atom, iconst, index_a, &
                                                index_b, ndist
    REAL(KIND=dp)                            :: del_lam, dtby2, dtsqby2, &
                                                fdotf, fdotf_sum, imass1, &
                                                imass2, mass
    REAL(KIND=dp), DIMENSION(3)              :: fc_t, r12
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(dist_constraint_type), POINTER      :: dist_list( : )
    TYPE(local_dist_constraint_type), &
      POINTER                                :: ldist( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

    dtsqby2 = simpar%dt*simpar%dt*.5_dp
    dtby2 = simpar%dt*.5_dp
    molecule_kind => molecule % molecule_kind
    CALL get_molecule_kind ( molecule_kind, ndist = ndist, dist_list = dist_list )
    CALL get_molecule ( molecule, first_atom = first_atom, ldist=ldist )
    IF (ishake==1) THEN
       DO iconst = 1, ndist
          index_a = dist_list ( iconst ) % a + first_atom -1
          index_b = dist_list ( iconst ) % b + first_atom -1
          atomic_kind=>particle_set ( index_a ) % atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
          imass1 = 1.0_dp/mass
          atomic_kind=>particle_set ( index_b ) % atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
          imass2 = 1.0_dp/mass
          ldist ( iconst ) % f = -2.0_dp * ( ldist ( iconst ) % ra_old - &
               ldist ( iconst ) % rb_old )

          ! update positions with full multiplier
          pos ( :, index_a ) = pos ( :, index_a ) + &
               imass1*dtsqby2*ldist ( iconst ) % lambda * ldist ( iconst ) % f ( : )
          pos ( :, index_b ) = pos ( :, index_b ) - &
               imass2*dtsqby2*ldist ( iconst ) % lambda * ldist ( iconst ) % f ( : )

          ! update velocites with full multiplier
          vel ( :, index_a ) = vel ( :, index_a ) + &
               imass1*dtby2*ldist ( iconst ) % lambda * ldist ( iconst ) % f ( : ) 
          vel ( :, index_b ) = vel ( :, index_b ) - &
               imass2*dtby2*ldist ( iconst ) % lambda * ldist ( iconst ) % f ( : ) 
       END DO

    ELSE

       DO iconst = 1, ndist
          ! computing the constraint
          index_a = dist_list ( iconst ) % a + first_atom -1
          index_b = dist_list ( iconst ) % b + first_atom -1
          r12 = pos ( :, index_a )  -  pos ( :, index_b )
          ldist ( iconst ) % sigma = DOT_PRODUCT(r12,r12) - dist_list ( iconst ) % dab * &
               dist_list ( iconst ) % dab
          ! computing fc(dt)
          fc_t = -2.0_dp*r12
          atomic_kind=>particle_set ( index_a ) % atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
          imass1 = 1.0_dp/mass
          atomic_kind=>particle_set ( index_b ) % atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
          imass2 = 1.0_dp/mass
          fdotf = DOT_PRODUCT( ldist ( iconst ) % f , fc_t )
          fdotf_sum = fdotf*imass1 + fdotf*imass2
          del_lam = 2.0_dp*ldist ( iconst ) % sigma/(simpar%dt*simpar%dt*fdotf_sum)
          ldist ( iconst ) % lambda = ldist ( iconst ) % lambda + del_lam

          ! update positions by delta_lam
          pos ( :, index_a ) = pos ( :, index_a ) + &
               imass1*dtsqby2*del_lam * ldist ( iconst ) % f ( : )
          pos ( :, index_b ) = pos ( :, index_b ) - &
               imass2*dtsqby2*del_lam * ldist ( iconst ) % f ( : )

          ! update velocites by delta_lam
          vel ( :, index_a ) = vel ( :, index_a ) + &
               imass1*dtby2*del_lam * ldist ( iconst ) % f ( : ) 
          vel ( :, index_b ) = vel ( :, index_b ) - &
               imass2*dtby2*del_lam * ldist ( iconst ) % f ( : ) 
       END DO
    END IF

    ! computing the constraint and value of tolerance
    DO iconst = 1, ndist
       index_a = dist_list ( iconst ) % a + first_atom -1
       index_b = dist_list ( iconst ) % b + first_atom -1
       r12 ( : ) = pos ( :, index_a ) -  pos ( :, index_b )
       ldist ( iconst ) % sigma = DOT_PRODUCT(r12,r12) - dist_list ( iconst ) % dab * &
            dist_list ( iconst ) % dab
       max_sigma = MAX(ABS(ldist ( iconst ) % sigma),max_sigma)
    END DO

  END SUBROUTINE shake

!!****** constraint_clv/rattle [1.0] *
!!
!!   NAME
!!     rattle
!!
!!   FUNCTION
!!     rattle algorithm for collective variables constraints
!!     updates the multiplier one molecule type at a time
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE rattle ( molecule, particle_set, vel, simpar, max_sigma, irattle )

    TYPE(molecule_type), POINTER             :: molecule
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: vel( :, : )
    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar
    REAL(KIND=dp), INTENT(INOUT)             :: max_sigma
    INTEGER, INTENT(INOUT)                   :: irattle

    INTEGER                                  :: first_atom, iconst, index_a, &
                                                index_b, ndist
    REAL(KIND=dp)                            :: del_lam, dtby2, fdotf, &
                                                fdotf_sum, imass1, imass2, &
                                                mass
    REAL(KIND=dp), DIMENSION(3)              :: r12, v12
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(dist_constraint_type), POINTER      :: dist_list( : )
    TYPE(local_dist_constraint_type), &
      POINTER                                :: ldist( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

    dtby2 = simpar%dt*.5_dp
    molecule_kind => molecule % molecule_kind
    CALL get_molecule_kind ( molecule_kind, ndist = ndist, dist_list = dist_list )
    CALL get_molecule ( molecule, first_atom = first_atom, ldist=ldist )
    IF (irattle==1) THEN
       DO iconst = 1, ndist
          index_a = dist_list ( iconst ) % a + first_atom -1
          index_b = dist_list ( iconst ) % b + first_atom -1
          atomic_kind=>particle_set ( index_a ) % atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
          imass1 = 1.0_dp/mass
          atomic_kind=>particle_set ( index_b ) % atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
          imass2 = 1.0_dp/mass
          ldist ( iconst ) % f = -2.0_dp * ( particle_set ( index_a ) % r - &
               particle_set ( index_b ) % r )
          ! update velocites with full multiplier
          vel ( :, index_a ) = vel ( :, index_a ) + &
               imass1*dtby2*ldist ( iconst ) % lambda * ldist ( iconst ) % f ( : )
          vel ( :, index_b ) = vel ( :, index_b ) - &
               imass2*dtby2*ldist ( iconst ) % lambda * ldist ( iconst ) % f ( : )
       END DO

    ELSE

       DO iconst = 1, ndist

          ! computing sigma
          index_a = dist_list ( iconst ) % a + first_atom -1
          index_b = dist_list ( iconst ) % b + first_atom -1
          r12 ( : ) = particle_set ( index_a ) % r ( : ) -  particle_set ( index_b ) % r ( : )
          v12 ( : ) = vel ( :, index_a ) -  vel ( :, index_b ) 
          ldist ( iconst ) % sigma = 2.0_dp*DOT_PRODUCT(r12,v12)

          ! computing f(dt)
          ldist ( iconst ) % f = -2.0_dp * r12
          atomic_kind=>particle_set ( index_a ) % atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
          imass1 = 1.0_dp/mass
          atomic_kind=>particle_set ( index_b ) % atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
          imass2 = 1.0_dp/mass
          fdotf = DOT_PRODUCT( ldist ( iconst ) % f , ldist ( iconst ) % f )
          fdotf_sum = fdotf*imass1 + fdotf*imass2
          del_lam = 2.0_dp*ldist ( iconst ) % sigma/(simpar%dt*fdotf_sum)
          ldist ( iconst ) % lambda = ldist ( iconst ) % lambda + del_lam

          !update velocites with del_lam
          vel ( :, index_a ) = vel ( :, index_a ) + &
               imass1*dtby2*del_lam * ldist ( iconst ) % f ( : )
          vel ( :, index_b ) = vel ( :, index_b ) - &
               imass2*dtby2*del_lam * ldist ( iconst ) % f ( : )
       END DO
    END IF

    ! computing the constraint and value of the tolerance
    DO iconst = 1, ndist
       index_a = dist_list ( iconst ) % a + first_atom -1
       index_b = dist_list ( iconst ) % b + first_atom -1
       r12 = particle_set ( index_a ) % r -  particle_set ( index_b ) % r
       v12 = vel ( :, index_a ) -  vel ( :, index_b )
       ldist ( iconst ) % sigma = DOT_PRODUCT(r12,v12)
       max_sigma = MAX ( ABS ( ldist ( iconst ) % sigma ), max_sigma )
    END DO

  END SUBROUTINE rattle

!!****** constraint_clv/rattle_roll [1.0] *
!!
!!   NAME
!!     rattle_roll
!!
!!   FUNCTION
!!     rattle algorithm (box allowed to change) for collective variables constraints
!!     updates the multiplier one molecule type at a time
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE rattle_roll ( molecule, particle_set, vel, r_rattle,  &
       simpar, veps, max_sigma, irattle )

    TYPE(molecule_type), POINTER             :: molecule
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: vel( :, : )
    REAL(KIND=dp), INTENT(IN)                :: r_rattle( :, : )
    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar
    REAL(KIND=dp), INTENT(IN)                :: veps( :, : )
    REAL(KIND=dp), INTENT(INOUT)             :: max_sigma
    INTEGER, INTENT(INOUT)                   :: irattle

    INTEGER                                  :: first_atom, iconst, index_a, &
                                                index_b, ndist
    REAL(KIND=dp)                            :: del_lam, dtby2, fdotf, &
                                                fdotf_sum, imass1, imass2, &
                                                mass
    REAL(KIND=dp), DIMENSION(3)              :: f_roll, r1, r12, r2, v12
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(dist_constraint_type), POINTER      :: dist_list( : )
    TYPE(local_dist_constraint_type), &
      POINTER                                :: ldist( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

    dtby2 = simpar%dt*.5_dp
    molecule_kind => molecule % molecule_kind
    CALL get_molecule_kind ( molecule_kind, ndist = ndist, dist_list = dist_list )
    CALL get_molecule ( molecule, first_atom = first_atom, ldist=ldist )
    IF (irattle==1) THEN
       DO iconst = 1, SIZE ( dist_list )
          index_a = dist_list ( iconst ) % a + first_atom -1
          index_b = dist_list ( iconst ) % b + first_atom -1
          atomic_kind=>particle_set ( index_a ) % atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
          imass1 = 1.0_dp/mass
          atomic_kind=>particle_set ( index_b ) % atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
          imass2 = 1.0_dp/mass
          ldist ( iconst ) % f = -2.0_dp * ( particle_set ( index_a ) % r - &
               particle_set ( index_b ) % r )
          ! update velocites with full multiplier
          vel ( :, index_a ) = vel ( :, index_a ) + &
               imass1*dtby2*ldist ( iconst ) % lambda * ldist ( iconst ) % f ( : )
          vel ( :, index_b ) = vel ( :, index_b ) - &
               imass2*dtby2*ldist ( iconst ) % lambda * ldist ( iconst ) % f ( : )
       END DO
       max_sigma = -1.E+10_dp
    ELSE
       DO iconst = 1, SIZE ( dist_list )
          index_a = dist_list ( iconst ) % a + first_atom -1
          index_b = dist_list ( iconst ) % b + first_atom -1
          r1 ( : ) = particle_set ( index_a ) % r ( : ) 
          r2 ( : ) = particle_set ( index_b ) % r ( : ) 
          r12 ( : ) = r1 ( : ) - r2 ( : )
          v12 ( : ) = vel ( :, index_a ) -  vel ( :, index_b ) 

          !  rolling f(dt)
          ldist ( iconst ) % f = -2.0_dp * r12
          f_roll = MATVEC_3x3 ( r_rattle, ldist ( iconst ) % f )

          ! computing sigma
          v12 ( : ) = vel ( :, index_a ) + MATVEC_3X3 ( veps, r1 ( : ) ) &
               - vel ( :, index_b ) - MATVEC_3X3 ( veps, r2 ( : ) )
          ldist ( iconst ) % sigma = 2.0_dp * DOTPROD_3D ( r12, v12 )

          atomic_kind=>particle_set ( index_a ) % atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
          imass1 = 1.0_dp/mass
          atomic_kind=>particle_set ( index_b ) % atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
          imass2 = 1.0_dp/mass

          fdotf = DOT_PRODUCT( f_roll , ldist ( iconst ) % f )
          fdotf_sum = fdotf*imass1 + fdotf*imass2
          del_lam = 2.0_dp*ldist ( iconst ) % sigma/(simpar%dt*fdotf_sum)
          ldist ( iconst ) % lambda = ldist ( iconst ) % lambda + del_lam

          !update velocites with del_lam
          vel ( :, index_a ) = vel ( :, index_a ) + &
               Imass1*dtby2*del_lam * f_roll ( : )
          vel ( :, index_b ) = vel ( :, index_b ) - &
               imass2*dtby2*del_lam * f_roll ( : )
       END DO
       ! computing the constraint and value of the tolerance
       DO iconst = 1, SIZE ( dist_list )
          index_a = dist_list ( iconst ) % a + first_atom -1
          index_b = dist_list ( iconst ) % b + first_atom -1
          r1 ( : ) = particle_set ( index_a ) % r ( : ) 
          r2 ( : ) = particle_set ( index_b ) % r ( : ) 
          r12 ( : ) = r1 ( : ) - r2 ( : )
          v12 ( : ) = vel ( :, index_a ) + MATVEC_3X3 ( veps, r1 ( : ) ) &
               - vel ( :, index_b ) - MATVEC_3X3 ( veps, r2 ( : ) )
          ldist ( iconst ) % sigma = 2.0_dp * DOTPROD_3D ( r12, v12 )
          max_sigma = MAX ( ABS ( ldist ( iconst ) % sigma ), max_sigma )
       END DO
    END IF
  END SUBROUTINE rattle_roll

!!****** constraint_clv/shake_roll [1.0] *
!!
!!   NAME
!!     shake_roll
!!
!!   FUNCTION
!!     shake algorithm (box allowed to change) for collective variables constraints
!!     updates the multiplier one molecule type at a time
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE shake_roll( molecule, particle_set, pos, vel, r_shake, v_shake, &
       simpar, max_sigma, ishake )

    TYPE(molecule_type), POINTER             :: molecule
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: pos( :, : ), vel( :, : )
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: r_shake, v_shake
    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar
    REAL(KIND=dp), INTENT(INOUT)             :: max_sigma
    INTEGER, INTENT(INOUT)                   :: ishake

    INTEGER                                  :: first_atom, iconst, index_a, &
                                                index_b, ndist
    REAL(KIND=dp)                            :: del_lam, dtby2, dtsqby2, &
                                                fdotf, fdotf_sum, imass1, &
                                                imass2, mass
    REAL(KIND=dp), DIMENSION(3)              :: f_roll_r, f_roll_v, fc_t, r12
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(dist_constraint_type), POINTER      :: dist_list( : )
    TYPE(local_dist_constraint_type), &
      POINTER                                :: ldist( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

    dtsqby2 = simpar%dt*simpar%dt*.5_dp
    dtby2 = simpar%dt*.5_dp
    molecule_kind => molecule % molecule_kind
    CALL get_molecule_kind ( molecule_kind, ndist = ndist, dist_list = dist_list )
    CALL get_molecule ( molecule, first_atom = first_atom, ldist=ldist )
    IF (ishake==1) THEN
       DO iconst = 1, ndist
          index_a = dist_list ( iconst ) % a + first_atom -1
          index_b = dist_list ( iconst ) % b + first_atom -1
          atomic_kind=>particle_set ( index_a ) % atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
          imass1 = 1.0_dp/mass
          atomic_kind=>particle_set ( index_b ) % atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
          imass2 = 1.0_dp/mass
          ldist ( iconst ) % f = -2.0_dp * ( ldist ( iconst ) % ra_old - &
               ldist ( iconst ) % rb_old )

          ! rotate fconst:
          f_roll_r = MATVEC_3X3 ( r_shake, ldist ( iconst ) % f )
          f_roll_v = MATVEC_3X3 ( v_shake, ldist ( iconst ) % f )
          ! update positions with full multiplier
          pos ( :, index_a ) = pos ( :, index_a ) + &
               imass1*dtsqby2*ldist ( iconst ) % lambda * f_roll_r
          pos ( :, index_b ) = pos ( :, index_b ) - &
               imass2*dtsqby2*ldist ( iconst ) % lambda * f_roll_r

          ! update velocites with full multiplier
          vel ( :, index_a ) = vel ( :, index_a ) + &
               imass1*dtby2*ldist ( iconst ) % lambda * f_roll_v
          vel ( :, index_b ) = vel ( :, index_b ) - &
               imass2*dtby2*ldist ( iconst ) % lambda * f_roll_v
       END DO

    ELSE

       DO iconst = 1, ndist
          ! computing the constraint
          index_a = dist_list ( iconst ) % a + first_atom -1
          index_b = dist_list ( iconst ) % b + first_atom -1
          r12 = pos ( :, index_a ) -  pos ( :, index_b )
          ldist ( iconst ) % sigma = DOT_PRODUCT(r12,r12) - dist_list ( iconst ) % dab * &
               dist_list ( iconst ) % dab
          ! computing fc(dt)
          fc_t = -2.0_dp*r12
          atomic_kind=>particle_set ( index_a ) % atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
          imass1 = 1.0_dp/mass
          atomic_kind=>particle_set ( index_b ) % atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
          imass2 = 1.0_dp/mass
          ! rotate fc:
          f_roll_r = MATVEC_3X3 ( r_shake, ldist ( iconst ) % f )
          f_roll_v = MATVEC_3X3 ( v_shake, ldist ( iconst ) % f )
          fdotf = DOTPROD_3D ( f_roll_r, fc_t )
          fdotf_sum = fdotf * imass1 + fdotf * imass2
          del_lam = 2.0_dp*ldist ( iconst ) % sigma/(simpar%dt*simpar%dt*fdotf_sum)
          ldist ( iconst ) % lambda = ldist ( iconst ) % lambda + del_lam

          ! update positions by delta_lam
          pos ( :, index_a ) = pos ( :, index_a ) + &
               imass1*dtsqby2*del_lam * f_roll_r
          pos ( :, index_b ) = pos ( :, index_b ) - &
               imass2*dtsqby2*del_lam * f_roll_r

          ! update velocites by delta_lam
          vel ( :, index_a ) = vel ( :, index_a ) + &
               imass1*dtby2*del_lam * f_roll_v
          vel ( :, index_b ) = vel ( :, index_b ) - &
               imass2*dtby2*del_lam * f_roll_v
       END DO
    END IF

    ! computing the constraint and value of tolerance
    DO iconst = 1, ndist
       index_a = dist_list ( iconst ) % a + first_atom -1
       index_b = dist_list ( iconst ) % b + first_atom -1
       r12 ( : ) = pos ( :, index_a ) -  pos ( :, index_b )
       ldist ( iconst ) % sigma = DOT_PRODUCT(r12,r12) - dist_list ( iconst ) % dab * &
            dist_list ( iconst ) % dab
       max_sigma = MAX(ABS(ldist ( iconst ) % sigma),max_sigma)
    END DO

  END SUBROUTINE shake_roll

END MODULE constraint_clv
