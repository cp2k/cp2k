!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      none
! *****************************************************************************
MODULE constraint_3x3

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE constraint_fxd,                  ONLY: check_fixed_atom_cns_g3x3
  USE kinds,                           ONLY: dp
  USE linear_systems,                  ONLY: solve_system
  USE mathlib,                         ONLY: dotprod_3d,&
                                             matvec_3x3
  USE molecule_kind_types,             ONLY: fixd_constraint_type,&
                                             g3x3_constraint_type,&
                                             get_molecule_kind,&
                                             molecule_kind_type
  USE molecule_types_new,              ONLY: get_molecule,&
                                             global_constraint_type,&
                                             local_g3x3_constraint_type,&
                                             molecule_type
  USE particle_types,                  ONLY: particle_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: shake_3x3_int,&
            rattle_3x3_int,&
            shake_roll_3x3_int,&
            rattle_roll_3x3_int,&
            shake_3x3_ext,&
            rattle_3x3_ext,&
            shake_roll_3x3_ext,&
            rattle_roll_3x3_ext

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'constraint_3x3'

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param molecule ...
!> \param particle_set ...
!> \param pos ...
!> \param vel ...
!> \param dt ...
!> \param ishake ...
!> \param max_sigma ...
!> \par History
!>      none
! *****************************************************************************
  SUBROUTINE shake_3x3_int( molecule, particle_set, pos, vel, dt, ishake,&
       max_sigma)

    TYPE(molecule_type), POINTER             :: molecule
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: pos( :, : ), vel( :, : )
    REAL(kind=dp), INTENT(in)                :: dt
    INTEGER, INTENT(IN)                      :: ishake
    REAL(KIND=dp), INTENT(INOUT)             :: max_sigma

    INTEGER                                  :: first_atom, ng3x3
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(g3x3_constraint_type), POINTER      :: g3x3_list( : )
    TYPE(local_g3x3_constraint_type), &
      POINTER                                :: lg3x3( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

    NULLIFY(fixd_list)
    molecule_kind => molecule % molecule_kind
    CALL get_molecule_kind ( molecule_kind, ng3x3 = ng3x3,&
         g3x3_list = g3x3_list, fixd_list=fixd_list )
    CALL get_molecule ( molecule, first_atom = first_atom, lg3x3=lg3x3 )
    ! Real Shake
    CALL shake_3x3_low( fixd_list, g3x3_list, lg3x3, first_atom, ng3x3, &
       particle_set, pos, vel, dt, ishake, max_sigma)

  END SUBROUTINE shake_3x3_int

! *****************************************************************************
!> \brief ...
!> \param molecule ...
!> \param particle_set ...
!> \param pos ...
!> \param vel ...
!> \param r_shake ...
!> \param v_shake ...
!> \param dt ...
!> \param ishake ...
!> \param max_sigma ...
!> \par History
!>      none
! *****************************************************************************
  SUBROUTINE shake_roll_3x3_int( molecule, particle_set, pos, vel, r_shake,  &
       v_shake, dt, ishake, max_sigma)

    TYPE(molecule_type), POINTER             :: molecule
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: pos( :, : ), vel( :, : )
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: r_shake, v_shake
    REAL(kind=dp), INTENT(in)                :: dt
    INTEGER, INTENT(IN)                      :: ishake
    REAL(KIND=dp), INTENT(INOUT)             :: max_sigma

    INTEGER                                  :: first_atom, ng3x3
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(g3x3_constraint_type), POINTER      :: g3x3_list( : )
    TYPE(local_g3x3_constraint_type), &
      POINTER                                :: lg3x3( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

    NULLIFY(fixd_list)
    molecule_kind => molecule % molecule_kind
    CALL get_molecule_kind ( molecule_kind, ng3x3 = ng3x3,&
         g3x3_list = g3x3_list, fixd_list=fixd_list )
    CALL get_molecule ( molecule, first_atom = first_atom, lg3x3=lg3x3 )
    ! Real Shake
    CALL shake_roll_3x3_low( fixd_list, g3x3_list, lg3x3, first_atom, ng3x3, &
       particle_set, pos, vel, r_shake, v_shake, dt, ishake, max_sigma)

  END SUBROUTINE shake_roll_3x3_int

! *****************************************************************************
!> \brief ...
!> \param molecule ...
!> \param particle_set ...
!> \param vel ...
!> \param r_rattle ...
!> \param dt ...
!> \param veps ...
!> \par History
!>      none
! *****************************************************************************
  SUBROUTINE rattle_roll_3x3_int( molecule, particle_set, vel, r_rattle, dt, veps)

    TYPE(molecule_type), POINTER             :: molecule
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: vel( :, : )
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: r_rattle
    REAL(kind=dp), INTENT(in)                :: dt
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: veps

    INTEGER                                  :: first_atom
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(g3x3_constraint_type), POINTER      :: g3x3_list( : )
    TYPE(local_g3x3_constraint_type), &
      POINTER                                :: lg3x3( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

    NULLIFY(fixd_list)
    molecule_kind => molecule % molecule_kind
    CALL get_molecule_kind ( molecule_kind, &
         g3x3_list = g3x3_list, fixd_list=fixd_list )
    CALL get_molecule ( molecule, first_atom = first_atom, lg3x3=lg3x3 )
    ! Real Rattle
    CALL rattle_roll_3x3_low( fixd_list, g3x3_list, lg3x3, first_atom, &
       particle_set, vel, r_rattle, dt, veps)

  END SUBROUTINE rattle_roll_3x3_int

! *****************************************************************************
!> \brief ...
!> \param molecule ...
!> \param particle_set ...
!> \param vel ...
!> \param dt ...
!> \par History
!>      none
! *****************************************************************************
  SUBROUTINE rattle_3x3_int( molecule, particle_set, vel, dt)

    TYPE(molecule_type), POINTER             :: molecule
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: vel( :, : )
    REAL(kind=dp), INTENT(in)                :: dt

    INTEGER                                  :: first_atom
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(g3x3_constraint_type), POINTER      :: g3x3_list( : )
    TYPE(local_g3x3_constraint_type), &
      POINTER                                :: lg3x3( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

    NULLIFY(fixd_list)
    molecule_kind => molecule % molecule_kind
    CALL get_molecule_kind ( molecule_kind, &
         g3x3_list = g3x3_list, fixd_list=fixd_list )
    CALL get_molecule ( molecule, first_atom = first_atom, lg3x3=lg3x3 )
    ! Real Rattle
    CALL rattle_3x3_low( fixd_list, g3x3_list, lg3x3, first_atom, &
       particle_set, vel, dt)

  END SUBROUTINE rattle_3x3_int

! *****************************************************************************
!> \brief ...
!> \param gci ...
!> \param particle_set ...
!> \param pos ...
!> \param vel ...
!> \param dt ...
!> \param ishake ...
!> \param max_sigma ...
!> \par History
!>      none
! *****************************************************************************
  SUBROUTINE shake_3x3_ext( gci, particle_set, pos, vel, dt, ishake,&
       max_sigma)

    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: pos( :, : ), vel( :, : )
    REAL(kind=dp), INTENT(in)                :: dt
    INTEGER, INTENT(IN)                      :: ishake
    REAL(KIND=dp), INTENT(INOUT)             :: max_sigma

    INTEGER                                  :: first_atom, ng3x3
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(g3x3_constraint_type), POINTER      :: g3x3_list( : )
    TYPE(local_g3x3_constraint_type), &
      POINTER                                :: lg3x3( : )

    first_atom = 1
    ng3x3 = gci%ng3x3
    g3x3_list => gci%g3x3_list
    fixd_list => gci%fixd_list
    lg3x3     => gci%lg3x3
    ! Real Shake
    CALL shake_3x3_low( fixd_list, g3x3_list, lg3x3, first_atom, ng3x3, &
       particle_set, pos, vel, dt, ishake, max_sigma)

  END SUBROUTINE shake_3x3_ext

! *****************************************************************************
!> \brief ...
!> \param gci ...
!> \param particle_set ...
!> \param pos ...
!> \param vel ...
!> \param r_shake ...
!> \param v_shake ...
!> \param dt ...
!> \param ishake ...
!> \param max_sigma ...
!> \par History
!>      none
! *****************************************************************************
  SUBROUTINE shake_roll_3x3_ext(gci, particle_set, pos, vel, r_shake,  &
       v_shake, dt, ishake, max_sigma)

    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: pos( :, : ), vel( :, : )
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: r_shake, v_shake
    REAL(kind=dp), INTENT(in)                :: dt
    INTEGER, INTENT(IN)                      :: ishake
    REAL(KIND=dp), INTENT(INOUT)             :: max_sigma

    INTEGER                                  :: first_atom, ng3x3
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(g3x3_constraint_type), POINTER      :: g3x3_list( : )
    TYPE(local_g3x3_constraint_type), &
      POINTER                                :: lg3x3( : )

    first_atom = 1
    ng3x3 = gci%ng3x3
    g3x3_list => gci%g3x3_list
    fixd_list => gci%fixd_list
    lg3x3     => gci%lg3x3
    ! Real Shake
    CALL shake_roll_3x3_low( fixd_list, g3x3_list, lg3x3, first_atom, ng3x3, &
       particle_set, pos, vel, r_shake, v_shake, dt, ishake, max_sigma)

  END SUBROUTINE shake_roll_3x3_ext

! *****************************************************************************
!> \brief ...
!> \param gci ...
!> \param particle_set ...
!> \param vel ...
!> \param r_rattle ...
!> \param dt ...
!> \param veps ...
!> \par History
!>      none
! *****************************************************************************
  SUBROUTINE rattle_roll_3x3_ext( gci, particle_set, vel, r_rattle, dt, veps)

    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: vel( :, : )
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: r_rattle
    REAL(kind=dp), INTENT(in)                :: dt
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: veps

    INTEGER                                  :: first_atom
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(g3x3_constraint_type), POINTER      :: g3x3_list( : )
    TYPE(local_g3x3_constraint_type), &
      POINTER                                :: lg3x3( : )

    first_atom = 1
    g3x3_list => gci%g3x3_list
    fixd_list => gci%fixd_list
    lg3x3     => gci%lg3x3
    ! Real Rattle
    CALL rattle_roll_3x3_low( fixd_list, g3x3_list, lg3x3, first_atom, &
       particle_set, vel, r_rattle, dt, veps)

  END SUBROUTINE rattle_roll_3x3_ext

! *****************************************************************************
!> \brief ...
!> \param gci ...
!> \param particle_set ...
!> \param vel ...
!> \param dt ...
!> \par History
!>      none
! *****************************************************************************
  SUBROUTINE rattle_3x3_ext( gci, particle_set, vel, dt)

    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: vel( :, : )
    REAL(kind=dp), INTENT(in)                :: dt

    INTEGER                                  :: first_atom
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(g3x3_constraint_type), POINTER      :: g3x3_list( : )
    TYPE(local_g3x3_constraint_type), &
      POINTER                                :: lg3x3( : )

    first_atom = 1
    g3x3_list => gci%g3x3_list
    fixd_list => gci%fixd_list
    lg3x3     => gci%lg3x3
    ! Real Rattle
    CALL rattle_3x3_low( fixd_list, g3x3_list, lg3x3, first_atom, &
       particle_set, vel, dt)

  END SUBROUTINE rattle_3x3_ext

! *****************************************************************************
!> \brief ...
!> \param fixd_list ...
!> \param g3x3_list ...
!> \param lg3x3 ...
!> \param first_atom ...
!> \param ng3x3 ...
!> \param particle_set ...
!> \param pos ...
!> \param vel ...
!> \param dt ...
!> \param ishake ...
!> \param max_sigma ...
!> \par History
!>      none
! *****************************************************************************
  SUBROUTINE shake_3x3_low( fixd_list, g3x3_list, lg3x3, first_atom, ng3x3, &
       particle_set, pos, vel, dt, ishake, max_sigma)

    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(g3x3_constraint_type), POINTER      :: g3x3_list( : )
    TYPE(local_g3x3_constraint_type), &
      POINTER                                :: lg3x3( : )
    INTEGER, INTENT(IN)                      :: first_atom, ng3x3
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: pos( :, : ), vel( :, : )
    REAL(kind=dp), INTENT(in)                :: dt
    INTEGER, INTENT(IN)                      :: ishake
    REAL(KIND=dp), INTENT(INOUT)             :: max_sigma

    INTEGER                                  :: iconst, index_a, index_b, &
                                                index_c
    REAL(KIND=dp)                            :: dtby2, dtsqby2, idtsq, &
                                                imass1, imass2, imass3, sigma
    REAL(KIND=dp), DIMENSION(3)              :: fc1, fc2, fc3, r0_12, r0_13, &
                                                r0_23, r1, r12, r13, r2, r23, &
                                                r3, v1, v2, v3, vec
    REAL(KIND=dp), DIMENSION(3, 1)           :: bvec
    REAL(KIND=dp), DIMENSION(3, 3)           :: amat, atemp
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

    dtsqby2 = dt*dt*.5_dp
    idtsq = 1.0_dp/dt/dt
    dtby2 = dt*.5_dp
    DO iconst = 1, ng3x3
       IF (g3x3_list(iconst)%restraint%active) CYCLE
       index_a = g3x3_list ( iconst ) % a + first_atom -1
       index_b = g3x3_list ( iconst ) % b + first_atom -1
       index_c = g3x3_list ( iconst ) % c + first_atom -1
       IF (ishake==1) THEN
          r0_12 ( : ) = pos ( :, index_a ) - pos ( :,index_b )
          r0_13 ( : ) = pos ( :, index_a ) - pos ( :,index_c )
          r0_23 ( : ) = pos ( :, index_b ) - pos ( :,index_c )
          atomic_kind => particle_set ( index_a ) % atomic_kind
          imass1 = 1.0_dp/atomic_kind % mass
          atomic_kind => particle_set ( index_b ) % atomic_kind
          imass2 = 1.0_dp/atomic_kind % mass
          atomic_kind => particle_set ( index_c ) % atomic_kind
          imass3 = 1.0_dp/atomic_kind % mass
          lg3x3 ( iconst ) % fa = -2.0_dp*( lg3x3 ( iconst ) % ra_old - &
               lg3x3 ( iconst ) % rb_old )
          lg3x3 ( iconst ) % fb = -2.0_dp*( lg3x3 ( iconst ) % ra_old - &
               lg3x3 ( iconst ) % rc_old )
          lg3x3 ( iconst ) % fc = -2.0_dp*( lg3x3 ( iconst ) % rb_old - &
               lg3x3 ( iconst ) % rc_old )
          ! Check for fixed atom constraints
          CALL check_fixed_atom_cns_g3x3(imass1, imass2, imass3,&
               index_a, index_b, index_c, fixd_list, lg3x3(iconst))
          ! construct matrix
          amat(1,1)=(imass1+imass2)*DOTPROD_3D(r0_12,lg3x3 ( iconst ) % fa )
          amat(1,2)=imass1*DOTPROD_3D(r0_12, lg3x3 ( iconst) % fb )
          amat(1,3)=-imass2*DOTPROD_3D(r0_12, lg3x3 ( iconst ) % fc )
          amat(2,1)=imass1*DOTPROD_3D(r0_13, lg3x3 ( iconst ) % fa )
          amat(2,2)=(imass1+imass3)*DOTPROD_3D(r0_13,lg3x3 ( iconst ) % fb )
          amat(2,3)=imass3*DOTPROD_3D(r0_13,lg3x3 ( iconst ) % fc )
          amat(3,1)=-imass2*DOTPROD_3D(r0_23, lg3x3 ( iconst ) % fa )
          amat(3,2)=imass3*DOTPROD_3D(r0_23, lg3x3 ( iconst ) % fb )
          amat(3,3)=(imass3+imass2)*DOTPROD_3D(r0_23,lg3x3 ( iconst ) % fc )
          ! Store values
          lg3x3 ( iconst ) % r0_12  = r0_12
          lg3x3 ( iconst ) % r0_13  = r0_13
          lg3x3 ( iconst ) % r0_23  = r0_23
          lg3x3 ( iconst ) % amat   = amat
          lg3x3 ( iconst ) % lambda_old = 0.0_dp
          lg3x3 ( iconst ) % del_lambda = 0.0_dp
       ELSE
          ! Retrieve values
          r0_12  = lg3x3 ( iconst ) % r0_12
          r0_13  = lg3x3 ( iconst ) % r0_13
          r0_23  = lg3x3 ( iconst ) % r0_23
          amat   = lg3x3 ( iconst ) % amat
          imass1 = lg3x3 ( iconst ) % imass1
          imass2 = lg3x3 ( iconst ) % imass2
          imass3 = lg3x3 ( iconst ) % imass3
       END IF
       ! Iterate until convergence
       vec= lg3x3 ( iconst ) % lambda ( 1 ) * lg3x3 ( iconst ) % fa *(imass1+imass2) + &
            lg3x3 ( iconst ) % lambda ( 2 ) * imass1*lg3x3 ( iconst ) % fb - &
            lg3x3 ( iconst ) % lambda ( 3 ) * imass2*lg3x3 ( iconst ) % fc
       bvec(1,1)=g3x3_list ( iconst ) % dab * g3x3_list ( iconst) % dab &
            -dtsqby2*dtsqby2*DOTPROD_3D(vec,vec)-DOTPROD_3D(r0_12,r0_12)
       vec= lg3x3 ( iconst ) % lambda ( 1 ) * lg3x3 ( iconst ) % fa * imass1 + &
            lg3x3 ( iconst ) % lambda ( 2 ) * ( imass1 + imass3 ) * lg3x3 ( iconst ) % fb + &
            lg3x3 ( iconst ) % lambda ( 3 ) * imass3 * lg3x3 ( iconst ) % fc
       bvec(2,1)=g3x3_list ( iconst ) % dac * g3x3_list ( iconst) % dac &
            -dtsqby2*dtsqby2*DOTPROD_3D(vec,vec)-DOTPROD_3D(r0_13,r0_13)
       vec= -lg3x3 ( iconst ) % lambda ( 1 ) * lg3x3 ( iconst ) % fa * imass2 + &
            lg3x3 ( iconst ) % lambda ( 2 ) * imass3 * lg3x3 ( iconst ) % fb + &
            lg3x3 ( iconst ) % lambda ( 3 ) * ( imass2 + imass3 ) * lg3x3 ( iconst ) % fc
       bvec(3,1)=g3x3_list ( iconst ) % dbc * g3x3_list ( iconst) % dbc &
            -dtsqby2*dtsqby2*DOTPROD_3D(vec,vec)-DOTPROD_3D(r0_23,r0_23)
       bvec = bvec * idtsq
       ! get lambda
       atemp = amat
       CALL solve_system ( atemp, 3, bvec )
       lg3x3 ( iconst ) % lambda ( : )    = bvec ( :, 1 )
       lg3x3 ( iconst ) % del_lambda( : ) = lg3x3 ( iconst ) % lambda ( : ) -&
                                            lg3x3 ( iconst ) % lambda_old ( : )
       lg3x3 ( iconst ) % lambda_old ( : )= lg3x3 ( iconst ) % lambda ( : )
       fc1= lg3x3 ( iconst ) % del_lambda ( 1 ) * lg3x3 ( iconst ) % fa + &
            lg3x3 ( iconst ) % del_lambda ( 2 ) * lg3x3 ( iconst ) % fb
       fc2=-lg3x3 ( iconst ) % del_lambda ( 1 ) * lg3x3 ( iconst ) % fa + &
            lg3x3 ( iconst ) % del_lambda ( 3 ) * lg3x3 ( iconst ) % fc
       fc3=-lg3x3 ( iconst ) % del_lambda ( 2 ) * lg3x3 ( iconst ) % fb - &
            lg3x3 ( iconst ) % del_lambda ( 3 ) * lg3x3 ( iconst ) % fc
       r1 ( : ) = pos ( :, index_a ) + imass1*dtsqby2*fc1 ( : )
       r2 ( : ) = pos ( :, index_b ) + imass2*dtsqby2*fc2 ( : )
       r3 ( : ) = pos ( :, index_c ) + imass3*dtsqby2*fc3 ( : )
       v1 ( : ) = vel ( :, index_a ) + imass1*dtby2*fc1 ( : )
       v2 ( : ) = vel ( :, index_b ) + imass2*dtby2*fc2 ( : )
       v3 ( : ) = vel ( :, index_c ) + imass3*dtby2*fc3 ( : )
       r12=r1-r2
       r13=r1-r3
       r23=r2-r3

       ! compute the tolerance:
       sigma = DOT_PRODUCT(r12,r12) - g3x3_list ( iconst ) % dab *  &
            g3x3_list ( iconst ) % dab
       max_sigma=MAX(max_sigma,ABS(sigma))
       sigma = DOT_PRODUCT(r13,r13) -  g3x3_list ( iconst ) % dac * &
            g3x3_list ( iconst ) % dac
       max_sigma=MAX(max_sigma,ABS(sigma))
       sigma = DOT_PRODUCT(r23,r23) - g3x3_list ( iconst ) % dbc *  &
            g3x3_list ( iconst ) % dbc
       max_sigma=MAX(max_sigma,ABS(sigma))
       ! update positions with full multiplier
       pos ( :, index_a ) = r1 ( : )
       pos ( :, index_b ) = r2 ( : )
       pos ( :, index_c ) = r3 ( : )

       ! update velocites with full multiplier
       vel ( :, index_a ) = v1 ( : )
       vel ( :, index_b ) = v2 ( : )
       vel ( :, index_c ) = v3 ( : )
    END DO

  END SUBROUTINE shake_3x3_low

! *****************************************************************************
!> \brief ...
!> \param fixd_list ...
!> \param g3x3_list ...
!> \param lg3x3 ...
!> \param first_atom ...
!> \param ng3x3 ...
!> \param particle_set ...
!> \param pos ...
!> \param vel ...
!> \param r_shake ...
!> \param v_shake ...
!> \param dt ...
!> \param ishake ...
!> \param max_sigma ...
!> \par History
!>      none
! *****************************************************************************
  SUBROUTINE shake_roll_3x3_low( fixd_list, g3x3_list, lg3x3, first_atom, ng3x3, &
       particle_set, pos, vel, r_shake, v_shake, dt, ishake, max_sigma)

    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(g3x3_constraint_type), POINTER      :: g3x3_list( : )
    TYPE(local_g3x3_constraint_type), &
      POINTER                                :: lg3x3( : )
    INTEGER, INTENT(IN)                      :: first_atom, ng3x3
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: pos( :, : ), vel( :, : )
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: r_shake, v_shake
    REAL(kind=dp), INTENT(in)                :: dt
    INTEGER, INTENT(IN)                      :: ishake
    REAL(KIND=dp), INTENT(INOUT)             :: max_sigma

    INTEGER                                  :: iconst, index_a, index_b, &
                                                index_c
    REAL(KIND=dp)                            :: dtby2, dtsqby2, idtsq, &
                                                imass1, imass2, imass3, sigma
    REAL(KIND=dp), DIMENSION(3)              :: f_roll1, f_roll2, f_roll3, &
                                                fc1, fc2, fc3, r0_12, r0_13, &
                                                r0_23, r1, r12, r13, r2, r23, &
                                                r3, v1, v2, v3, vec
    REAL(KIND=dp), DIMENSION(3, 1)           :: bvec
    REAL(KIND=dp), DIMENSION(3, 3)           :: amat, atemp
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

    dtsqby2 = dt*dt*.5_dp
    idtsq = 1.0_dp/dt/dt
    dtby2 = dt*.5_dp
    DO iconst = 1, ng3x3
       IF (g3x3_list(iconst)%restraint%active) CYCLE
       index_a = g3x3_list ( iconst ) % a + first_atom -1
       index_b = g3x3_list ( iconst ) % b + first_atom -1
       index_c = g3x3_list ( iconst ) % c + first_atom -1
       IF (ishake==1) THEN
          r0_12 ( : ) = pos ( :, index_a ) - pos ( :,index_b )
          r0_13 ( : ) = pos ( :, index_a ) - pos ( :,index_c )
          r0_23 ( : ) = pos ( :, index_b ) - pos ( :,index_c )
          atomic_kind => particle_set ( index_a ) % atomic_kind
          imass1 = 1.0_dp/atomic_kind % mass
          atomic_kind => particle_set ( index_b ) % atomic_kind
          imass2 = 1.0_dp/atomic_kind % mass
          atomic_kind => particle_set ( index_c ) % atomic_kind
          imass3 = 1.0_dp/atomic_kind % mass
          lg3x3 ( iconst ) % fa = -2.0_dp*( lg3x3 ( iconst ) % ra_old - &
               lg3x3 ( iconst ) % rb_old )
          lg3x3 ( iconst ) % fb = -2.0_dp*( lg3x3 ( iconst ) % ra_old - &
               lg3x3 ( iconst ) % rc_old )
          lg3x3 ( iconst ) % fc = -2.0_dp*( lg3x3 ( iconst ) % rb_old - &
               lg3x3 ( iconst ) % rc_old )
          ! Check for fixed atom constraints
          CALL check_fixed_atom_cns_g3x3(imass1, imass2, imass3,&
               index_a, index_b, index_c, fixd_list, lg3x3(iconst))
          ! rotate fconst:
          CALL matvec_3x3 (f_roll1, r_shake, lg3x3 ( iconst ) % fa )
          CALL matvec_3x3 (f_roll2, r_shake, lg3x3 ( iconst ) % fb )
          CALL matvec_3x3 (f_roll3, r_shake, lg3x3 ( iconst ) % fc )
          ! construct matrix
          amat(1,1)=(imass1+imass2)*DOTPROD_3D(r0_12,f_roll1 )
          amat(1,2)=imass1*DOTPROD_3D(r0_12, f_roll2 )
          amat(1,3)=-imass2*DOTPROD_3D(r0_12, f_roll3 )
          amat(2,1)=imass1*DOTPROD_3D(r0_13, f_roll1 )
          amat(2,2)=(imass1+imass3)*DOTPROD_3D(r0_13,f_roll2 )
          amat(2,3)=imass3*DOTPROD_3D(r0_13,f_roll3)
          amat(3,1)=-imass2*DOTPROD_3D(r0_23, f_roll1 )
          amat(3,2)=imass3*DOTPROD_3D(r0_23, f_roll2 )
          amat(3,3)=(imass3+imass2)*DOTPROD_3D(r0_23,f_roll3)
          ! Store values
          lg3x3 ( iconst ) % r0_12  = r0_12
          lg3x3 ( iconst ) % r0_13  = r0_13
          lg3x3 ( iconst ) % r0_23  = r0_23
          lg3x3 ( iconst ) % f_roll1 = f_roll1
          lg3x3 ( iconst ) % f_roll2 = f_roll2
          lg3x3 ( iconst ) % f_roll3 = f_roll3
          lg3x3 ( iconst ) % amat   = amat
          lg3x3 ( iconst ) % lambda_old = 0.0_dp
          lg3x3 ( iconst ) % del_lambda = 0.0_dp
       ELSE
          ! Retrieve values
          r0_12  = lg3x3 ( iconst ) % r0_12
          r0_13  = lg3x3 ( iconst ) % r0_13
          r0_23  = lg3x3 ( iconst ) % r0_23
          f_roll1 = lg3x3 ( iconst ) % f_roll1
          f_roll2 = lg3x3 ( iconst ) % f_roll2
          f_roll3 = lg3x3 ( iconst ) % f_roll3
          amat   = lg3x3 ( iconst ) % amat
          imass1 = lg3x3 ( iconst ) % imass1
          imass2 = lg3x3 ( iconst ) % imass2
          imass3 = lg3x3 ( iconst ) % imass3
       END IF
       ! Iterate until convergence
       vec= lg3x3 ( iconst ) % lambda ( 1 ) * f_roll1 *(imass1+imass2) + &
            lg3x3 ( iconst ) % lambda ( 2 ) * imass1* f_roll2 - &
            lg3x3 ( iconst ) % lambda ( 3 ) * imass2* f_roll3
       bvec(1,1)=g3x3_list ( iconst ) % dab * g3x3_list ( iconst) % dab &
            -dtsqby2*dtsqby2*DOTPROD_3D(vec,vec)-DOTPROD_3D(r0_12,r0_12)
       vec= lg3x3 ( iconst ) % lambda ( 1 ) * f_roll1 * imass1 + &
            lg3x3 ( iconst ) % lambda ( 2 ) * ( imass1 + imass3 ) * f_roll2 + &
            lg3x3 ( iconst ) % lambda ( 3 ) * imass3 * f_roll3
       bvec(2,1)=g3x3_list ( iconst ) % dac * g3x3_list ( iconst) % dac &
            -dtsqby2*dtsqby2*DOTPROD_3D(vec,vec)-DOTPROD_3D(r0_13,r0_13)
       vec= -lg3x3 ( iconst ) % lambda ( 1 ) * f_roll1 * imass2 + &
            lg3x3 ( iconst ) % lambda ( 2 ) * imass3 * f_roll2 + &
            lg3x3 ( iconst ) % lambda ( 3 ) * ( imass2 + imass3 ) * f_roll3
       bvec(3,1)=g3x3_list ( iconst ) % dbc * g3x3_list ( iconst) % dbc &
            -dtsqby2*dtsqby2*DOTPROD_3D(vec,vec)-DOTPROD_3D(r0_23,r0_23)
       bvec = bvec * idtsq

       ! get lambda
       atemp = amat
       CALL solve_system ( atemp, 3, bvec )
       lg3x3 ( iconst ) % lambda ( : )    = bvec ( :, 1 )
       lg3x3 ( iconst ) % del_lambda( : ) = lg3x3 ( iconst ) % lambda ( : ) -&
                                            lg3x3 ( iconst ) % lambda_old ( : )
       lg3x3 ( iconst ) % lambda_old ( : )= lg3x3 ( iconst ) % lambda ( : )

       fc1= lg3x3 ( iconst ) % del_lambda ( 1 ) * lg3x3 ( iconst ) % fa + &
            lg3x3 ( iconst ) % del_lambda ( 2 ) * lg3x3 ( iconst ) % fb
       fc2=-lg3x3 ( iconst ) % del_lambda ( 1 ) * lg3x3 ( iconst ) % fa + &
            lg3x3 ( iconst ) % del_lambda ( 3 ) * lg3x3 ( iconst ) % fc
       fc3=-lg3x3 ( iconst ) % del_lambda ( 2 ) * lg3x3 ( iconst ) % fb - &
            lg3x3 ( iconst ) % del_lambda ( 3 ) * lg3x3 ( iconst ) % fc
       CALL MATVEC_3x3 (vec, r_shake, fc1 )
       r1 ( : ) = pos ( :, index_a ) + imass1*dtsqby2*vec
       CALL MATVEC_3x3 (vec, r_shake, fc2 )
       r2 ( : ) = pos ( :, index_b ) + imass2*dtsqby2*vec
       CALL MATVEC_3x3 (vec, r_shake, fc3 )
       r3 ( : ) = pos ( :, index_c ) + imass3*dtsqby2*vec
       CALL MATVEC_3x3 (vec, v_shake, fc1 )
       v1 ( : ) = vel ( :, index_a ) + imass1*dtby2*vec
       CALL MATVEC_3x3 (vec, v_shake, fc2 )
       v2 ( : ) = vel ( :, index_b ) + imass2*dtby2*vec
       CALL MATVEC_3x3 (vec, v_shake, fc3 )
       v3 ( : ) = vel ( :, index_c ) + imass3*dtby2*vec
       r12=r1-r2
       r13=r1-r3
       r23=r2-r3

       ! compute the tolerance:
       sigma = DOT_PRODUCT(r12,r12) - g3x3_list ( iconst ) % dab *  &
            g3x3_list ( iconst ) % dab
       max_sigma=MAX(max_sigma,ABS(sigma))
       sigma = DOT_PRODUCT(r13,r13) -  g3x3_list ( iconst ) % dac * &
            g3x3_list ( iconst ) % dac
       max_sigma=MAX(max_sigma,ABS(sigma))
       sigma = DOT_PRODUCT(r23,r23) - g3x3_list ( iconst ) % dbc *  &
            g3x3_list ( iconst ) % dbc
       max_sigma=MAX(max_sigma,ABS(sigma))

       ! update positions with full multiplier
       pos ( :, index_a ) = r1 ( : )
       pos ( :, index_b ) = r2 ( : )
       pos ( :, index_c ) = r3 ( : )

       ! update velocites with full multiplier
       vel ( :, index_a ) = v1 ( : )
       vel ( :, index_b ) = v2 ( : )
       vel ( :, index_c ) = v3 ( : )
    END DO
  END SUBROUTINE shake_roll_3x3_low

! *****************************************************************************
!> \brief ...
!> \param fixd_list ...
!> \param g3x3_list ...
!> \param lg3x3 ...
!> \param first_atom ...
!> \param particle_set ...
!> \param vel ...
!> \param r_rattle ...
!> \param dt ...
!> \param veps ...
!> \par History
!>      none
! *****************************************************************************
  SUBROUTINE rattle_roll_3x3_low( fixd_list, g3x3_list, lg3x3, first_atom, &
       particle_set, vel, r_rattle, dt, veps)
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(g3x3_constraint_type), POINTER      :: g3x3_list( : )
    TYPE(local_g3x3_constraint_type), &
      POINTER                                :: lg3x3( : )
    INTEGER, INTENT(IN)                      :: first_atom
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: vel( :, : )
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: r_rattle
    REAL(kind=dp), INTENT(in)                :: dt
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: veps

    INTEGER                                  :: iconst, index_a, index_b, &
                                                index_c
    REAL(KIND=dp)                            :: dtby2, idt, imass1, imass2, &
                                                imass3, mass
    REAL(KIND=dp), DIMENSION(3)              :: f_roll1, f_roll2, f_roll3, &
                                                fc1, fc2, fc3, lambda, r12, &
                                                r13, r23, v12, v13, v23, vec
    REAL(KIND=dp), DIMENSION(3, 1)           :: bvec
    REAL(KIND=dp), DIMENSION(3, 3)           :: amat
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

    idt = 1.0_dp/dt
    dtby2 = dt*.5_dp
    DO iconst = 1, SIZE ( g3x3_list )
       IF (g3x3_list(iconst)%restraint%active) CYCLE
       index_a = g3x3_list ( iconst ) % a + first_atom - 1
       index_b = g3x3_list ( iconst ) % b + first_atom - 1
       index_c = g3x3_list ( iconst ) % c + first_atom - 1
       v12 ( : ) = vel ( :, index_a ) - vel ( :, index_b )
       v13 ( : ) = vel ( :, index_a ) - vel ( :, index_c )
       v23 ( : ) = vel ( :, index_b ) - vel ( :, index_c )
       r12 ( : ) = particle_set ( index_a ) % r ( : ) - particle_set ( index_b ) % r ( : )
       r13 ( : ) = particle_set ( index_a ) % r ( : ) - particle_set ( index_c ) % r ( : )
       r23 ( : ) = particle_set ( index_b ) % r ( : ) - particle_set ( index_c ) % r ( : )
       atomic_kind=>particle_set ( index_a ) % atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
       imass1 = 1.0_dp/mass
       atomic_kind=>particle_set ( index_b ) % atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
       imass2 = 1.0_dp/mass
       atomic_kind=>particle_set ( index_c ) % atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
       imass3 = 1.0_dp/mass
       lg3x3 ( iconst ) % fa = -2.0_dp* r12
       lg3x3 ( iconst ) % fb = -2.0_dp* r13
       lg3x3 ( iconst ) % fc = -2.0_dp* r23
       ! Check for fixed atom constraints
       CALL check_fixed_atom_cns_g3x3(imass1, imass2, imass3,&
            index_a, index_b, index_c, fixd_list, lg3x3(iconst))
       ! roll the fc
       CALL MATVEC_3x3 (f_roll1, r_rattle, lg3x3 ( iconst ) % fa )
       CALL MATVEC_3x3 (f_roll2, r_rattle, lg3x3 ( iconst ) % fb )
       CALL MATVEC_3x3 (f_roll3, r_rattle, lg3x3 ( iconst ) % fc )

       ! construct matrix
       amat(1,1)=(imass1+imass2)*DOTPROD_3D(r12, f_roll1 )
       amat(1,2)=imass1*DOTPROD_3D(r12, f_roll2 )
       amat(1,3)=-imass2*DOTPROD_3D(r12,f_roll3 )
       amat(2,1)=imass1*DOTPROD_3D(r13, f_roll1 )
       amat(2,2)=(imass1+imass3)*DOTPROD_3D(r13, f_roll2 )
       amat(2,3)=imass3*DOTPROD_3D(r13, f_roll3 )
       amat(3,1)=-imass2*DOTPROD_3D(r23,f_roll1 )
       amat(3,2)=imass3*DOTPROD_3D(r23, f_roll2 )
       amat(3,3)=(imass2+imass3)*DOTPROD_3D(r23, f_roll3 )

       ! construct solution vector
       CALL matvec_3x3 (vec, veps, r12 )
       bvec ( 1, 1 ) = DOTPROD_3D ( r12, v12 + vec )
       CALL matvec_3x3 (vec, veps, r13 )
       bvec ( 2, 1 ) = DOTPROD_3D ( r13, v13 + vec )
       CALL matvec_3x3 (vec, veps, r23 )
       bvec ( 3, 1 ) = DOTPROD_3D ( r23, v23 + vec )
       bvec = -bvec * 2.0_dp * idt

       ! get lambda
       CALL solve_system ( amat, 3, bvec )
       lambda ( : ) = bvec ( :, 1 )
       lg3x3 ( iconst ) % lambda ( : )= lambda

       fc1= lambda ( 1 ) * f_roll1 + &
            lambda ( 2 ) * f_roll2
       fc2=-lambda ( 1 ) * f_roll1 + &
            lambda ( 3 ) * f_roll3
       fc3=-lambda ( 2 ) * f_roll2 - &
            lambda ( 3 ) * f_roll3
       vel ( :, index_a ) = vel ( :, index_a ) + imass1 * dtby2 * fc1 ( : )
       vel ( :, index_b ) = vel ( :, index_b ) + imass2 * dtby2 * fc2 ( : )
       vel ( :, index_c ) = vel ( :, index_c ) + imass3 * dtby2 * fc3 ( : )
    END DO
  END SUBROUTINE rattle_roll_3x3_low

! *****************************************************************************
!> \brief ...
!> \param fixd_list ...
!> \param g3x3_list ...
!> \param lg3x3 ...
!> \param first_atom ...
!> \param particle_set ...
!> \param vel ...
!> \param dt ...
!> \par History
!>      none
! *****************************************************************************
  SUBROUTINE rattle_3x3_low( fixd_list, g3x3_list, lg3x3, first_atom, &
       particle_set, vel, dt)

    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(g3x3_constraint_type), POINTER      :: g3x3_list( : )
    TYPE(local_g3x3_constraint_type), &
      POINTER                                :: lg3x3( : )
    INTEGER, INTENT(IN)                      :: first_atom
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: vel( :, : )
    REAL(kind=dp), INTENT(in)                :: dt

    INTEGER                                  :: iconst, index_a, index_b, &
                                                index_c
    REAL(KIND=dp)                            :: dtby2, idt, imass1, imass2, &
                                                imass3, mass
    REAL(KIND=dp), DIMENSION(3)              :: fc1, fc2, fc3, r12, r13, r23, &
                                                v12, v13, v23
    REAL(KIND=dp), DIMENSION(3, 1)           :: bvec
    REAL(KIND=dp), DIMENSION(3, 3)           :: amat
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

    idt = 1.0_dp/dt
    dtby2 = dt*.5_dp
    DO iconst = 1, SIZE ( g3x3_list )
       IF (g3x3_list(iconst)%restraint%active) CYCLE
       index_a = g3x3_list ( iconst ) % a + first_atom - 1
       index_b = g3x3_list ( iconst ) % b + first_atom - 1
       index_c = g3x3_list ( iconst ) % c + first_atom - 1
       v12 ( : ) = vel ( :, index_a ) - vel ( :, index_b )
       v13 ( : ) = vel ( :, index_a ) - vel ( :, index_c )
       v23 ( : ) = vel ( :, index_b ) - vel ( :, index_c )
       r12 ( : ) = particle_set ( index_a ) % r ( : ) - particle_set ( index_b ) % r ( : )
       r13 ( : ) = particle_set ( index_a ) % r ( : ) - particle_set ( index_c ) % r ( : )
       r23 ( : ) = particle_set ( index_b ) % r ( : ) - particle_set ( index_c ) % r ( : )
       atomic_kind=>particle_set ( index_a ) % atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
       imass1 = 1.0_dp/mass
       atomic_kind=>particle_set ( index_b ) % atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
       imass2 = 1.0_dp/mass
       atomic_kind=>particle_set ( index_c ) % atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
       imass3 = 1.0_dp/mass
       lg3x3 ( iconst ) % fa = -2.0_dp* r12
       lg3x3 ( iconst ) % fb = -2.0_dp* r13
       lg3x3 ( iconst ) % fc = -2.0_dp* r23
       ! Check for fixed atom constraints
       CALL check_fixed_atom_cns_g3x3(imass1, imass2, imass3,&
            index_a, index_b, index_c, fixd_list, lg3x3(iconst))
       ! construct matrix
       amat(1,1)=(imass1+imass2)*DOTPROD_3D(r12,lg3x3 ( iconst ) % fa )
       amat(1,2)=imass1*DOTPROD_3D(r12,lg3x3 ( iconst ) % fb )
       amat(1,3)=-imass2*DOTPROD_3D(r12,lg3x3 ( iconst ) % fc )
       amat(2,1)=imass1*DOTPROD_3D(r13,lg3x3 ( iconst ) % fa )
       amat(2,2)=(imass1+imass3)*DOTPROD_3D(r13,lg3x3 ( iconst ) % fb )
       amat(2,3)=imass3*DOTPROD_3D(r13,lg3x3 ( iconst ) % fc )
       amat(3,1)=-imass2*DOTPROD_3D(r23,lg3x3 ( iconst ) % fa )
       amat(3,2)=imass3*DOTPROD_3D(r23,lg3x3 ( iconst ) % fb )
       amat(3,3)=(imass2+imass3)*DOTPROD_3D(r23,lg3x3 ( iconst ) % fc )

       ! construct solution vector
       bvec(1,1)=DOTPROD_3D(r12,v12)
       bvec(2,1)=DOTPROD_3D(r13,v13)
       bvec(3,1)=DOTPROD_3D(r23,v23)
       bvec=-bvec*2.0_dp*idt

       ! get lambda
       CALL solve_system ( amat, 3, bvec )
       lg3x3 ( iconst ) % lambda ( : )=bvec ( :, 1 )

       fc1= lg3x3 ( iconst ) % lambda ( 1 ) * lg3x3 ( iconst ) % fa + &
            lg3x3 ( iconst ) % lambda ( 2 ) * lg3x3 ( iconst ) % fb
       fc2=-lg3x3 ( iconst ) % lambda ( 1 ) * lg3x3 ( iconst ) % fa + &
            lg3x3 ( iconst ) % lambda ( 3 ) * lg3x3 ( iconst ) % fc
       fc3=-lg3x3 ( iconst ) % lambda ( 2 ) * lg3x3 ( iconst ) % fb - &
            lg3x3 ( iconst ) % lambda ( 3 ) * lg3x3 ( iconst ) % fc
       vel ( :, index_a ) = vel ( :, index_a ) + imass1 * dtby2 * fc1 ( : )
       vel ( :, index_b ) = vel ( :, index_b ) + imass2 * dtby2 * fc2 ( : )
       vel ( :, index_c ) = vel ( :, index_c ) + imass3 * dtby2 * fc3 ( : )
    END DO
  END SUBROUTINE rattle_3x3_low

END MODULE constraint_3x3
