!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****** cp2k/constraint_4x6 [1.0] *
!!
!!   NAME
!!     constraint_4x6
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
MODULE constraint_4x6

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE kinds,                           ONLY: dp
  USE linear_systems,                  ONLY: solve_system
  USE md,                              ONLY: simulation_parameters_type
  USE molecule_kind_types,             ONLY: g4x6_constraint_type,&
                                             get_molecule_kind,&
                                             molecule_kind_type
  USE molecule_types_new,              ONLY: get_molecule,&
                                             local_g4x6_constraint_type,&
                                             molecule_type
  USE particle_types,                  ONLY: particle_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: dotprod_3d,&
                                             matvec_3x3
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: shake_4x6,&
            rattle_4x6,&
            shake_roll_4x6,&
            rattle_roll_4x6

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'constraint_4x6'

CONTAINS

!!****** constraint_4x6/shake_4x6 [1.0] *
!!
!!   NAME
!!     shake_4x6
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
  SUBROUTINE shake_4x6( molecule, particle_set, pos, vel, simpar )
    TYPE(molecule_type), POINTER             :: molecule
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: pos( :, : ), vel( :, : )
    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar

    INTEGER                                  :: first_atom, iconst, index_a, &
                                                index_b, index_c, index_d, &
                                                ng4x6
    REAL(KIND=dp)                            :: dtby2, dtsqby2, idtsq, &
                                                imass1, imass2, imass3, &
                                                imass4, mass, max_sigma, sigma
    REAL(KIND=dp), DIMENSION(3) :: fc1, fc2, fc3, fc4, r0_12, r0_13, r0_14, &
      r0_23, r0_24, r0_34, r1, r12, r13, r14, r2, r23, r24, r3, r34, r4, v1, &
      v2, v3, v4, vec
    REAL(KIND=dp), DIMENSION(6, 1)           :: bvec
    REAL(KIND=dp), DIMENSION(6, 6)           :: amat, atemp
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(g4x6_constraint_type), POINTER      :: g4x6_list( : )
    TYPE(local_g4x6_constraint_type), &
      POINTER                                :: lg4x6( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    INTEGER                                  :: ishake

    dtsqby2 = simpar%dt*simpar%dt*.5_dp
    idtsq = 1.0_dp/simpar%dt/simpar%dt
    dtby2 = simpar%dt*.5_dp
    molecule_kind => molecule % molecule_kind
    CALL get_molecule_kind ( molecule_kind, ng4x6 = ng4x6, g4x6_list = g4x6_list )
    CALL get_molecule ( molecule, first_atom = first_atom, lg4x6=lg4x6 )
    DO iconst = 1, ng4x6
       ishake=0
       index_a = g4x6_list ( iconst ) % a + first_atom -1
       index_b = g4x6_list ( iconst ) % b + first_atom -1
       index_c = g4x6_list ( iconst ) % c + first_atom -1
       index_d = g4x6_list ( iconst ) % d + first_atom -1
       r0_12 ( : ) = pos ( :, index_a ) - pos ( :,index_b ) 
       r0_13 ( : ) = pos ( :, index_a ) - pos ( :,index_c ) 
       r0_14 ( : ) = pos ( :, index_a ) - pos ( :,index_d ) 
       r0_23 ( : ) = pos ( :, index_b ) - pos ( :,index_c ) 
       r0_24 ( : ) = pos ( :, index_b ) - pos ( :,index_d ) 
       r0_34 ( : ) = pos ( :, index_c ) - pos ( :,index_d ) 
       atomic_kind => particle_set ( index_a ) % atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
       imass1 = 1.0_dp/mass
       atomic_kind => particle_set ( index_b ) % atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
       imass2 = 1.0_dp/mass
       atomic_kind => particle_set ( index_c ) % atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
       imass3 = 1.0_dp/mass
       atomic_kind => particle_set ( index_d ) % atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
       imass4 = 1.0_dp/mass
       lg4x6 ( iconst ) % fa = -2.0_dp*( lg4x6 ( iconst ) % ra_old - &
            lg4x6 ( iconst ) % rb_old )
       lg4x6 ( iconst ) % fb = -2.0_dp*( lg4x6 ( iconst ) % ra_old - &
            lg4x6 ( iconst ) % rc_old )
       lg4x6 ( iconst ) % fc = -2.0_dp*( lg4x6 ( iconst ) % ra_old - &
            lg4x6 ( iconst ) % rd_old )
       lg4x6 ( iconst ) % fd = -2.0_dp*( lg4x6 ( iconst ) % rb_old - &
            lg4x6 ( iconst ) % rc_old )
       lg4x6 ( iconst ) % fe = -2.0_dp*( lg4x6 ( iconst ) % rb_old - &
            lg4x6 ( iconst ) % rd_old )
       lg4x6 ( iconst ) % ff = -2.0_dp*( lg4x6 ( iconst ) % rc_old - &
            lg4x6 ( iconst ) % rd_old )
       ! construct matrix
       amat(1,1)=(imass1+imass2)*DOTPROD_3D(r0_12,lg4x6 ( iconst ) % fa )
       amat(1,2)=imass1*DOTPROD_3D(r0_12, lg4x6 ( iconst) % fb )
       amat(1,3)=imass1*DOTPROD_3D(r0_12, lg4x6 ( iconst ) % fc )
       amat(1,4)=-imass2*DOTPROD_3D(r0_12, lg4x6 ( iconst ) % fd )
       amat(1,5)=-imass2*DOTPROD_3D(r0_12, lg4x6 ( iconst ) % fe )
       amat(1,6)=0.0_dp

       amat(2,1)=imass1*DOTPROD_3D(r0_13,lg4x6 ( iconst ) % fa )
       amat(2,2)=(imass1+imass3)*DOTPROD_3D(r0_13, lg4x6 ( iconst) % fb )
       amat(2,3)=imass1*DOTPROD_3D(r0_13, lg4x6 ( iconst ) % fc )
       amat(2,4)=imass3*DOTPROD_3D(r0_13, lg4x6 ( iconst ) % fd )
       amat(2,5)=0.0_dp
       amat(2,6)=-imass3*DOTPROD_3D(r0_13, lg4x6 ( iconst ) % ff )

       amat(3,1)=imass1*DOTPROD_3D(r0_14,lg4x6 ( iconst ) % fa )
       amat(3,2)=imass1*DOTPROD_3D(r0_14, lg4x6 ( iconst) % fb )
       amat(3,3)=(imass1+imass4)*DOTPROD_3D(r0_14, lg4x6 ( iconst ) % fc )
       amat(3,4)=0.0_dp
       amat(3,5)=imass4*DOTPROD_3D(r0_14, lg4x6 ( iconst ) % fe )
       amat(3,6)=imass4*DOTPROD_3D(r0_14, lg4x6 ( iconst ) % ff )

       amat(4,1)=-imass2*DOTPROD_3D(r0_23,lg4x6 ( iconst ) % fa )
       amat(4,2)=imass3*DOTPROD_3D(r0_23, lg4x6 ( iconst) % fb )
       amat(4,3)=0.0_dp
       amat(4,4)=(imass3+imass2)*DOTPROD_3D(r0_23, lg4x6 ( iconst ) % fd )
       amat(4,5)=imass2*DOTPROD_3D(r0_23, lg4x6 ( iconst ) % fe )
       amat(4,6)=-imass3*DOTPROD_3D(r0_23, lg4x6 ( iconst ) % ff )

       amat(5,1)=-imass2*DOTPROD_3D(r0_24,lg4x6 ( iconst ) % fa )
       amat(5,2)=0.0_dp
       amat(5,3)=imass4*DOTPROD_3D(r0_24, lg4x6 ( iconst) % fc )
       amat(5,4)=imass2*DOTPROD_3D(r0_24, lg4x6 ( iconst ) % fd )
       amat(5,5)=(imass4+imass2)*DOTPROD_3D(r0_24, lg4x6 ( iconst ) % fe )
       amat(5,6)=imass4*DOTPROD_3D(r0_24, lg4x6 ( iconst ) % ff )

       amat(6,1)=0.0_dp
       amat(6,2)=-imass3*DOTPROD_3D(r0_34,lg4x6 ( iconst ) % fb )
       amat(6,3)=imass4*DOTPROD_3D(r0_34, lg4x6 ( iconst) % fc )
       amat(6,4)=-imass3*DOTPROD_3D(r0_34, lg4x6 ( iconst ) % fd )
       amat(6,5)=imass4*DOTPROD_3D(r0_34, lg4x6 ( iconst ) % fe )
       amat(6,6)=(imass3+imass4)*DOTPROD_3D(r0_34, lg4x6 ( iconst ) % ff )

       ! Iterate until convergence:
       max_sigma = -1.0E+10_dp
       CONS: DO WHILE (ABS(max_sigma)>=simpar%shake_tol)
          max_sigma=0.0_dp
          ishake=ishake+1
          vec= lg4x6 ( iconst ) % lambda ( 1 ) * lg4x6 ( iconst ) % fa *(imass1+imass2) + &
               lg4x6 ( iconst ) % lambda ( 2 ) * imass1*lg4x6 ( iconst ) % fb + &
               lg4x6 ( iconst ) % lambda ( 3 ) * imass1*lg4x6 ( iconst ) % fc - &
               lg4x6 ( iconst ) % lambda ( 4 ) * imass2*lg4x6 ( iconst ) % fd - &
               lg4x6 ( iconst ) % lambda ( 5 ) * imass2*lg4x6 ( iconst ) % fe
          bvec(1,1)=g4x6_list ( iconst ) % dab * g4x6_list ( iconst) % dab &
               -dtsqby2*dtsqby2*DOTPROD_3D(vec,vec)-DOTPROD_3D(r0_12,r0_12)
          vec= lg4x6 ( iconst ) % lambda ( 2 ) * lg4x6 ( iconst ) % fb *(imass1+imass3) + &
               lg4x6 ( iconst ) % lambda ( 1 ) * imass1*lg4x6 ( iconst ) % fa + &
               lg4x6 ( iconst ) % lambda ( 3 ) * imass1*lg4x6 ( iconst ) % fc + &
               lg4x6 ( iconst ) % lambda ( 4 ) * imass3*lg4x6 ( iconst ) % fd - &
               lg4x6 ( iconst ) % lambda ( 6 ) * imass3*lg4x6 ( iconst ) % ff
          bvec(2,1)=g4x6_list ( iconst ) % dac * g4x6_list ( iconst) % dac &
               -dtsqby2*dtsqby2*DOTPROD_3D(vec,vec)-DOTPROD_3D(r0_13,r0_13)
          vec= lg4x6 ( iconst ) % lambda ( 3 ) * lg4x6 ( iconst ) % fc *(imass1+imass4) + &
               lg4x6 ( iconst ) % lambda ( 1 ) * imass1*lg4x6 ( iconst ) % fa + &
               lg4x6 ( iconst ) % lambda ( 2 ) * imass1*lg4x6 ( iconst ) % fb + &
               lg4x6 ( iconst ) % lambda ( 5 ) * imass4*lg4x6 ( iconst ) % fe + &
               lg4x6 ( iconst ) % lambda ( 6 ) * imass4*lg4x6 ( iconst ) % ff
          bvec(3,1)=g4x6_list ( iconst ) % dad * g4x6_list ( iconst) % dad &
               -dtsqby2*dtsqby2*DOTPROD_3D(vec,vec)-DOTPROD_3D(r0_14,r0_14)
          vec= lg4x6 ( iconst ) % lambda ( 4 ) * lg4x6 ( iconst ) % fd *(imass2+imass3) - &
               lg4x6 ( iconst ) % lambda ( 1 ) * imass2*lg4x6 ( iconst ) % fa + &
               lg4x6 ( iconst ) % lambda ( 2 ) * imass3*lg4x6 ( iconst ) % fb + &
               lg4x6 ( iconst ) % lambda ( 5 ) * imass2*lg4x6 ( iconst ) % fe - &
               lg4x6 ( iconst ) % lambda ( 6 ) * imass3*lg4x6 ( iconst ) % ff
          bvec(4,1)=g4x6_list ( iconst ) % dbc * g4x6_list ( iconst) % dbc &
               -dtsqby2*dtsqby2*DOTPROD_3D(vec,vec)-DOTPROD_3D(r0_23,r0_23)
          vec= lg4x6 ( iconst ) % lambda ( 5 ) * lg4x6 ( iconst ) % fe *(imass2+imass4) - &
               lg4x6 ( iconst ) % lambda ( 1 ) * imass2*lg4x6 ( iconst ) % fa + &
               lg4x6 ( iconst ) % lambda ( 3 ) * imass4*lg4x6 ( iconst ) % fc + &
               lg4x6 ( iconst ) % lambda ( 4 ) * imass2*lg4x6 ( iconst ) % fd + &
               lg4x6 ( iconst ) % lambda ( 6 ) * imass4*lg4x6 ( iconst ) % ff
          bvec(5,1)=g4x6_list ( iconst ) % dbd * g4x6_list ( iconst) % dbd &
               -dtsqby2*dtsqby2*DOTPROD_3D(vec,vec)-DOTPROD_3D(r0_24,r0_24)
          vec= lg4x6 ( iconst ) % lambda ( 6 ) * lg4x6 ( iconst ) % ff *(imass3+imass4) - &
               lg4x6 ( iconst ) % lambda ( 2 ) * imass3*lg4x6 ( iconst ) % fb + &
               lg4x6 ( iconst ) % lambda ( 3 ) * imass4*lg4x6 ( iconst ) % fc - &
               lg4x6 ( iconst ) % lambda ( 4 ) * imass3*lg4x6 ( iconst ) % fd + &
               lg4x6 ( iconst ) % lambda ( 5 ) * imass4*lg4x6 ( iconst ) % fe
          bvec(6,1)=g4x6_list ( iconst ) % dcd * g4x6_list ( iconst) % dcd &
               -dtsqby2*dtsqby2*DOTPROD_3D(vec,vec)-DOTPROD_3D(r0_34,r0_34)
          bvec = bvec * idtsq

          ! get lambda
          atemp = amat
          CALL solve_system ( atemp, 6, bvec )
          lg4x6 ( iconst ) % lambda ( : ) = bvec ( :, 1 )

          fc1=lg4x6 ( iconst ) % lambda ( 1 ) * lg4x6 ( iconst ) % fa + &
               lg4x6 ( iconst ) % lambda ( 2 ) * lg4x6 ( iconst ) % fb + &
               lg4x6 ( iconst ) % lambda ( 3 ) * lg4x6 ( iconst ) % fc 
          fc2=-lg4x6 ( iconst ) % lambda ( 1 ) * lg4x6 ( iconst ) % fa + &
               lg4x6 ( iconst ) % lambda ( 4 ) * lg4x6 ( iconst ) % fd + &
               lg4x6 ( iconst ) % lambda ( 5 ) * lg4x6 ( iconst ) % fe 
          fc3=-lg4x6 ( iconst ) % lambda ( 2 ) * lg4x6 ( iconst ) % fb - &
               lg4x6 ( iconst ) % lambda ( 4 ) * lg4x6 ( iconst ) % fd + &
               lg4x6 ( iconst ) % lambda ( 6 ) * lg4x6 ( iconst ) % ff 
          fc4=-lg4x6 ( iconst ) % lambda ( 3 ) * lg4x6 ( iconst ) % fc - &
               lg4x6 ( iconst ) % lambda ( 5 ) * lg4x6 ( iconst ) % fe - &
               lg4x6 ( iconst ) % lambda ( 6 ) * lg4x6 ( iconst ) % ff 
          r1 ( : ) = pos ( :, index_a ) + imass1*dtsqby2*fc1 ( : )
          r2 ( : ) = pos ( :, index_b ) + imass2*dtsqby2*fc2 ( : )
          r3 ( : ) = pos ( :, index_c ) + imass3*dtsqby2*fc3 ( : )
          r4 ( : ) = pos ( :, index_d ) + imass4*dtsqby2*fc4 ( : )
          v1 ( : ) = vel ( :, index_a ) + imass1*dtby2*fc1 ( : )
          v2 ( : ) = vel ( :, index_b ) + imass2*dtby2*fc2 ( : )
          v3 ( : ) = vel ( :, index_c ) + imass3*dtby2*fc3 ( : )
          v4 ( : ) = vel ( :, index_d ) + imass4*dtby2*fc4 ( : )
          r12=r1-r2
          r13=r1-r3
          r14=r1-r4
          r23=r2-r3
          r24=r2-r4
          r34=r3-r4
          ! compute the tolerance:
          sigma = DOT_PRODUCT(r12,r12) - g4x6_list ( iconst ) % dab *  &
               g4x6_list ( iconst ) % dab
          max_sigma=MAX(max_sigma,ABS(sigma))
          sigma = DOT_PRODUCT(r13,r13) - g4x6_list ( iconst ) % dac *  &
               g4x6_list ( iconst ) % dac
          max_sigma=MAX(max_sigma,ABS(sigma))
          sigma = DOT_PRODUCT(r14,r14) - g4x6_list ( iconst ) % dad *  &
               g4x6_list ( iconst ) % dad
          max_sigma=MAX(max_sigma,ABS(sigma))
          sigma = DOT_PRODUCT(r23,r23) - g4x6_list ( iconst ) % dbc *  &
               g4x6_list ( iconst ) % dbc
          max_sigma=MAX(max_sigma,ABS(sigma))
          sigma = DOT_PRODUCT(r24,r24) - g4x6_list ( iconst ) % dbd *  &
               g4x6_list ( iconst ) % dbd
          max_sigma=MAX(max_sigma,ABS(sigma))
          sigma = DOT_PRODUCT(r34,r34) - g4x6_list ( iconst ) % dcd *  &
               g4x6_list ( iconst ) % dcd
          max_sigma=MAX(max_sigma,ABS(sigma))
       END DO CONS

       ! update positions with full multiplier
       pos ( :, index_a ) = r1 ( : )
       pos ( :, index_b ) = r2 ( : )
       pos ( :, index_c ) = r3 ( : )
       pos ( :, index_d ) = r4 ( : )

       ! update velocites with full multiplier
       vel ( :, index_a ) = v1 ( : )
       vel ( :, index_b ) = v2 ( : )
       vel ( :, index_c ) = v3 ( : )
       vel ( :, index_d ) = v4 ( : )
    END DO

  END SUBROUTINE shake_4x6

!!****** constraint_4x6/shake_roll_4x6 [1.0] *
!!
!!   NAME
!!     shake_roll_4x6
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
  SUBROUTINE shake_roll_4x6( molecule, particle_set, pos, vel, r_shake,  &
       v_shake, simpar )
    TYPE(molecule_type), POINTER             :: molecule
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: pos( :, : ), vel( :, : )
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: r_shake, v_shake
    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar

    INTEGER                                  :: first_atom, iconst, index_a, &
                                                index_b, index_c, index_d, &
                                                ng4x6, ishake
    REAL(KIND=dp)                            :: dtby2, dtsqby2, idtsq, &
                                                imass1, imass2, imass3, &
                                                imass4, mass, max_sigma, sigma
    REAL(KIND=dp), DIMENSION(3) :: f_roll1, f_roll2, f_roll3, f_roll4, &
      f_roll5, f_roll6, fc1, fc2, fc3, fc4, r0_12, r0_13, r0_14, r0_23, &
      r0_24, r0_34, r1, r12, r13, r14, r2, r23, r24, r3, r34, r4, v1, v2, v3, &
      v4, vec
    REAL(KIND=dp), DIMENSION(6, 1)           :: bvec
    REAL(KIND=dp), DIMENSION(6, 6)           :: amat, atemp
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(g4x6_constraint_type), POINTER      :: g4x6_list( : )
    TYPE(local_g4x6_constraint_type), &
      POINTER                                :: lg4x6( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

    dtsqby2 = simpar%dt*simpar%dt*.5_dp
    idtsq = 1.0_dp/simpar%dt/simpar%dt
    dtby2 = simpar%dt*.5_dp
    molecule_kind => molecule % molecule_kind
    CALL get_molecule_kind ( molecule_kind, ng4x6 = ng4x6, g4x6_list = g4x6_list )
    CALL get_molecule ( molecule, first_atom = first_atom, lg4x6=lg4x6 )
    DO iconst = 1, ng4x6
       ishake=0
       index_a = g4x6_list ( iconst ) % a + first_atom -1
       index_b = g4x6_list ( iconst ) % b + first_atom -1
       index_c = g4x6_list ( iconst ) % c + first_atom -1
       index_d = g4x6_list ( iconst ) % d + first_atom -1
       r0_12 ( : ) = pos ( :, index_a ) - pos ( :,index_b ) 
       r0_13 ( : ) = pos ( :, index_a ) - pos ( :,index_c ) 
       r0_23 ( : ) = pos ( :, index_b ) - pos ( :,index_c ) 
       r0_14 ( : ) = pos ( :, index_a ) - pos ( :,index_d ) 
       r0_24 ( : ) = pos ( :, index_b ) - pos ( :,index_d ) 
       r0_34 ( : ) = pos ( :, index_c ) - pos ( :,index_d ) 
       atomic_kind => particle_set ( index_a ) % atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
       imass1 = 1.0_dp/mass
       atomic_kind => particle_set ( index_b ) % atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
       imass2 = 1.0_dp/mass
       atomic_kind => particle_set ( index_c ) % atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
       imass3 = 1.0_dp/mass
       atomic_kind => particle_set ( index_d ) % atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
       imass4 = 1.0_dp/mass
       lg4x6 ( iconst ) % fa = -2.0_dp*( lg4x6 ( iconst ) % ra_old - &
            lg4x6 ( iconst ) % rb_old )
       lg4x6 ( iconst ) % fb = -2.0_dp*( lg4x6 ( iconst ) % ra_old - &
            lg4x6 ( iconst ) % rc_old )
       lg4x6 ( iconst ) % fc = -2.0_dp*( lg4x6 ( iconst ) % ra_old - &
            lg4x6 ( iconst ) % rd_old )
       lg4x6 ( iconst ) % fd = -2.0_dp*( lg4x6 ( iconst ) % rb_old - &
            lg4x6 ( iconst ) % rc_old )
       lg4x6 ( iconst ) % fe = -2.0_dp*( lg4x6 ( iconst ) % rb_old - &
            lg4x6 ( iconst ) % rd_old )
       lg4x6 ( iconst ) % ff = -2.0_dp*( lg4x6 ( iconst ) % rc_old - &
            lg4x6 ( iconst ) % rd_old )
       ! rotate fconst:
       f_roll1 = matvec_3x3 ( r_shake, lg4x6 ( iconst ) % fa )
       f_roll2 = matvec_3x3 ( r_shake, lg4x6 ( iconst ) % fb )
       f_roll3 = matvec_3x3 ( r_shake, lg4x6 ( iconst ) % fc )
       f_roll4 = matvec_3x3 ( r_shake, lg4x6 ( iconst ) % fd )
       f_roll5 = matvec_3x3 ( r_shake, lg4x6 ( iconst ) % fe )
       f_roll6 = matvec_3x3 ( r_shake, lg4x6 ( iconst ) % ff )

       ! construct matrix
       amat(1,1)=(imass1+imass2)*DOTPROD_3D(r0_12, f_roll1)
       amat(1,2)=imass1*DOTPROD_3D(r0_12, f_roll2 )
       amat(1,3)=imass1*DOTPROD_3D(r0_12, f_roll3 )
       amat(1,4)=-imass2*DOTPROD_3D(r0_12, f_roll4 )
       amat(1,5)=-imass2*DOTPROD_3D(r0_12, f_roll5 )
       amat(1,6)=0.0_dp

       amat(2,1)=imass1*DOTPROD_3D(r0_13, f_roll1 )
       amat(2,2)=(imass1+imass3)*DOTPROD_3D(r0_13, f_roll2 )
       amat(2,3)=imass1*DOTPROD_3D(r0_13, f_roll3 )
       amat(2,4)=imass3*DOTPROD_3D(r0_13, f_roll4 )
       amat(2,5)=0.0_dp
       amat(2,6)=-imass3*DOTPROD_3D(r0_13, f_roll6 )

       amat(3,1)=imass1*DOTPROD_3D(r0_14, f_roll1 )
       amat(3,2)=imass1*DOTPROD_3D(r0_14, f_roll2 )
       amat(3,3)=(imass1+imass4)*DOTPROD_3D(r0_14, f_roll3 )
       amat(3,4)=0.0_dp
       amat(3,5)=imass4*DOTPROD_3D(r0_14, f_roll5 )
       amat(3,6)=imass4*DOTPROD_3D(r0_14, f_roll6 )

       amat(4,1)=-imass2*DOTPROD_3D(r0_23, f_roll1 )
       amat(4,2)=imass3*DOTPROD_3D(r0_23, f_roll2)
       amat(4,3)=0.0_dp
       amat(4,4)=(imass3+imass2)*DOTPROD_3D(r0_23, f_roll4 )
       amat(4,5)=imass2*DOTPROD_3D(r0_23, f_roll5 )
       amat(4,6)=-imass3*DOTPROD_3D(r0_23, f_roll6 )

       amat(5,1)=-imass2*DOTPROD_3D(r0_24, f_roll1 )
       amat(5,2)=0.0_dp
       amat(5,3)=imass4*DOTPROD_3D(r0_24, f_roll3 )
       amat(5,4)=imass2*DOTPROD_3D(r0_24, f_roll4 )
       amat(5,5)=(imass4+imass2)*DOTPROD_3D(r0_24, f_roll5 )
       amat(5,6)=imass4*DOTPROD_3D(r0_24, f_roll6 )

       amat(6,1)=0.0_dp
       amat(6,2)=-imass3*DOTPROD_3D(r0_34, f_roll2 )
       amat(6,3)=imass4*DOTPROD_3D(r0_34, f_roll3 )
       amat(6,4)=-imass3*DOTPROD_3D(r0_34, f_roll4 )
       amat(6,5)=imass4*DOTPROD_3D(r0_34, f_roll5 )
       amat(6,6)=(imass3+imass4)*DOTPROD_3D(r0_34, f_roll6 )

       ! Iterate until convergence:
       max_sigma = -1.0E+10_dp
       CONS: DO WHILE (ABS(max_sigma)>=simpar%shake_tol)
          max_sigma=0.0_dp
          ishake=ishake+1
          vec= lg4x6 ( iconst ) % lambda ( 1 ) * f_roll1*(imass1+imass2) + &
               lg4x6 ( iconst ) % lambda ( 2 ) * imass1*f_roll2 + &
               lg4x6 ( iconst ) % lambda ( 3 ) * imass1*f_roll3 - &
               lg4x6 ( iconst ) % lambda ( 4 ) * imass2*f_roll4 - &
               lg4x6 ( iconst ) % lambda ( 5 ) * imass2*f_roll5
          bvec(1,1)=g4x6_list ( iconst ) % dab * g4x6_list ( iconst) % dab &
               -dtsqby2*dtsqby2*DOTPROD_3D(vec,vec)-DOTPROD_3D(r0_12,r0_12)
          vec= lg4x6 ( iconst ) % lambda ( 2 ) * f_roll2 *(imass1+imass3) + &
               lg4x6 ( iconst ) % lambda ( 1 ) * imass1*f_roll1 + &
               lg4x6 ( iconst ) % lambda ( 3 ) * imass1*f_roll3 + &
               lg4x6 ( iconst ) % lambda ( 4 ) * imass3*f_roll4 - &
               lg4x6 ( iconst ) % lambda ( 6 ) * imass3*f_roll6
          bvec(2,1)=g4x6_list ( iconst ) % dac * g4x6_list ( iconst) % dac &
               -dtsqby2*dtsqby2*DOTPROD_3D(vec,vec)-DOTPROD_3D(r0_13,r0_13)
          vec= lg4x6 ( iconst ) % lambda ( 3 ) * f_roll3 *(imass1+imass4) + &
               lg4x6 ( iconst ) % lambda ( 1 ) * imass1*f_roll1 + &
               lg4x6 ( iconst ) % lambda ( 2 ) * imass1*f_roll2 + &
               lg4x6 ( iconst ) % lambda ( 5 ) * imass4*f_roll5 + &
               lg4x6 ( iconst ) % lambda ( 6 ) * imass4*f_roll6
          bvec(3,1)=g4x6_list ( iconst ) % dad * g4x6_list ( iconst) % dad &
               -dtsqby2*dtsqby2*DOTPROD_3D(vec,vec)-DOTPROD_3D(r0_14,r0_14)
          vec= lg4x6 ( iconst ) % lambda ( 4 ) * f_roll4 *(imass2+imass3) - &
               lg4x6 ( iconst ) % lambda ( 1 ) * imass2*f_roll1 + &
               lg4x6 ( iconst ) % lambda ( 2 ) * imass3*f_roll2 + &
               lg4x6 ( iconst ) % lambda ( 5 ) * imass2*f_roll5 - &
               lg4x6 ( iconst ) % lambda ( 6 ) * imass3*f_roll6
          bvec(4,1)=g4x6_list ( iconst ) % dbc * g4x6_list ( iconst) % dbc &
               -dtsqby2*dtsqby2*DOTPROD_3D(vec,vec)-DOTPROD_3D(r0_23,r0_23)
          vec= lg4x6 ( iconst ) % lambda ( 5 ) * f_roll5 *(imass2+imass4) - &
               lg4x6 ( iconst ) % lambda ( 1 ) * imass2*f_roll1 + &
               lg4x6 ( iconst ) % lambda ( 3 ) * imass4*f_roll3 + &
               lg4x6 ( iconst ) % lambda ( 4 ) * imass2*f_roll4 + &
               lg4x6 ( iconst ) % lambda ( 6 ) * imass4*f_roll6
          bvec(5,1)=g4x6_list ( iconst ) % dbd * g4x6_list ( iconst) % dbd &
               -dtsqby2*dtsqby2*DOTPROD_3D(vec,vec)-DOTPROD_3D(r0_24,r0_24)
          vec= lg4x6 ( iconst ) % lambda ( 6 ) * f_roll6 *(imass3+imass4) - &
               lg4x6 ( iconst ) % lambda ( 2 ) * imass3*f_roll2 + &
               lg4x6 ( iconst ) % lambda ( 3 ) * imass4*f_roll3 - &
               lg4x6 ( iconst ) % lambda ( 4 ) * imass3*f_roll4 + &
               lg4x6 ( iconst ) % lambda ( 5 ) * imass4*f_roll5
          bvec(6,1)=g4x6_list ( iconst ) % dcd * g4x6_list ( iconst) % dcd &
               -dtsqby2*dtsqby2*DOTPROD_3D(vec,vec)-DOTPROD_3D(r0_34,r0_34)
          bvec = bvec * idtsq

          ! get lambda
          atemp = amat
          CALL solve_system ( atemp, 6, bvec )
          lg4x6 ( iconst ) % lambda ( : ) = bvec ( :, 1 )

          fc1=lg4x6 ( iconst ) % lambda ( 1 ) * lg4x6 ( iconst ) % fa + &
               lg4x6 ( iconst ) % lambda ( 2 ) * lg4x6 ( iconst ) % fb + &
               lg4x6 ( iconst ) % lambda ( 3 ) * lg4x6 ( iconst ) % fc 
          fc2=-lg4x6 ( iconst ) % lambda ( 1 ) * lg4x6 ( iconst ) % fa + &
               lg4x6 ( iconst ) % lambda ( 4 ) * lg4x6 ( iconst ) % fd + &
               lg4x6 ( iconst ) % lambda ( 5 ) * lg4x6 ( iconst ) % fe 
          fc3=-lg4x6 ( iconst ) % lambda ( 2 ) * lg4x6 ( iconst ) % fb - &
               lg4x6 ( iconst ) % lambda ( 4 ) * lg4x6 ( iconst ) % fd + &
               lg4x6 ( iconst ) % lambda ( 6 ) * lg4x6 ( iconst ) % ff 
          fc4=-lg4x6 ( iconst ) % lambda ( 3 ) * lg4x6 ( iconst ) % fc - &
               lg4x6 ( iconst ) % lambda ( 5 ) * lg4x6 ( iconst ) % fe - &
               lg4x6 ( iconst ) % lambda ( 6 ) * lg4x6 ( iconst ) % ff 
          r1 ( : ) = pos ( :, index_a ) + imass1*dtsqby2*MATVEC_3x3(r_shake,fc1)
          r2 ( : ) = pos ( :, index_b ) + imass2*dtsqby2*MATVEC_3x3(r_shake,fc2)
          r3 ( : ) = pos ( :, index_c ) + imass3*dtsqby2*MATVEC_3x3(r_shake,fc3)
          r4 ( : ) = pos ( :, index_d ) + imass4*dtsqby2*MATVEC_3x3(r_shake,fc4)
          v1 ( : ) = vel ( :, index_a ) + imass1*dtby2*MATVEC_3x3(r_shake,fc1)
          v2 ( : ) = vel ( :, index_b ) + imass2*dtby2*MATVEC_3x3(r_shake,fc2)
          v3 ( : ) = vel ( :, index_c ) + imass3*dtby2*MATVEC_3x3(r_shake,fc3)
          v4 ( : ) = vel ( :, index_d ) + imass4*dtby2*MATVEC_3x3(r_shake,fc4)

          r12=r1-r2
          r13=r1-r3
          r23=r2-r3
          r14=r1-r4
          r24=r2-r4
          r34=r3-r4
          ! compute the tolerance:
          sigma = DOT_PRODUCT(r12,r12) - g4x6_list ( iconst ) % dab *  &
               g4x6_list ( iconst ) % dab
          max_sigma=MAX(max_sigma,ABS(sigma))
          sigma = DOT_PRODUCT(r13,r13) - g4x6_list ( iconst ) % dac *  &
               g4x6_list ( iconst ) % dac
          max_sigma=MAX(max_sigma,ABS(sigma))
          sigma = DOT_PRODUCT(r14,r14) - g4x6_list ( iconst ) % dad *  &
               g4x6_list ( iconst ) % dad
          max_sigma=MAX(max_sigma,ABS(sigma))
          sigma = DOT_PRODUCT(r23,r23) - g4x6_list ( iconst ) % dbc *  &
               g4x6_list ( iconst ) % dbc
          max_sigma=MAX(max_sigma,ABS(sigma))
          sigma = DOT_PRODUCT(r24,r24) - g4x6_list ( iconst ) % dbd *  &
               g4x6_list ( iconst ) % dbd
          max_sigma=MAX(max_sigma,ABS(sigma))
          sigma = DOT_PRODUCT(r34,r34) - g4x6_list ( iconst ) % dcd *  &
               g4x6_list ( iconst ) % dcd
          max_sigma=MAX(max_sigma,ABS(sigma))
       END DO CONS

       ! update positions with full multiplier
       pos ( :, index_a ) = r1 ( : )
       pos ( :, index_b ) = r2 ( : )
       pos ( :, index_c ) = r3 ( : )
       pos ( :, index_d ) = r4 ( : )

       ! update velocites with full multiplier
       vel ( :, index_a ) = v1 ( : )
       vel ( :, index_b ) = v2 ( : )
       vel ( :, index_c ) = v3 ( : )
       vel ( :, index_d ) = v4 ( : )
    END DO

  END SUBROUTINE shake_roll_4x6

!!****** constraint_4x6/rattle_4x6 [1.0] *
!!
!!   NAME
!!     rattle_4x6
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
  SUBROUTINE rattle_4x6( molecule, particle_set, vel, simpar)
    TYPE(molecule_type), POINTER             :: molecule
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: vel( :, : )
    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar

    INTEGER                                  :: first_atom, iconst, index_a, &
                                                index_b, index_c, index_d, &
                                                ng4x6
    REAL(KIND=dp)                            :: dtby2, idt, imass1, imass2, &
                                                imass3, imass4, mass
    REAL(KIND=dp), DIMENSION(3)              :: fc1, fc2, fc3, fc4, r12, r13, &
                                                r14, r23, r24, r34, v12, v13, &
                                                v14, v23, v24, v34
    REAL(KIND=dp), DIMENSION(6, 1)           :: bvec
    REAL(KIND=dp), DIMENSION(6, 6)           :: amat
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(g4x6_constraint_type), POINTER      :: g4x6_list( : )
    TYPE(local_g4x6_constraint_type), &
      POINTER                                :: lg4x6( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

    idt = 1.0_dp/simpar%dt
    dtby2 = simpar%dt*.5_dp
    molecule_kind => molecule % molecule_kind
    CALL get_molecule_kind ( molecule_kind, ng4x6 = ng4x6, g4x6_list = g4x6_list )
    CALL get_molecule ( molecule, first_atom = first_atom, lg4x6=lg4x6 )
    DO iconst = 1, SIZE ( g4x6_list )
       index_a = g4x6_list ( iconst ) % a + first_atom -1
       index_b = g4x6_list ( iconst ) % b + first_atom -1
       index_c = g4x6_list ( iconst ) % c + first_atom -1
       index_d = g4x6_list ( iconst ) % d + first_atom -1
       v12 ( : ) = vel ( :, index_a ) - vel ( :, index_b ) 
       v13 ( : ) = vel ( :, index_a ) - vel ( :, index_c ) 
       v14 ( : ) = vel ( :, index_a ) - vel ( :, index_d ) 
       v23 ( : ) = vel ( :, index_b ) - vel ( :, index_c ) 
       v24 ( : ) = vel ( :, index_b ) - vel ( :, index_d ) 
       v34 ( : ) = vel ( :, index_c ) - vel ( :, index_d ) 

       r12 ( : ) = particle_set ( index_a ) % r ( : ) - particle_set ( index_b ) % r ( : )
       r13 ( : ) = particle_set ( index_a ) % r ( : ) - particle_set ( index_c ) % r ( : )
       r14 ( : ) = particle_set ( index_a ) % r ( : ) - particle_set ( index_d ) % r ( : )
       r23 ( : ) = particle_set ( index_b ) % r ( : ) - particle_set ( index_c ) % r ( : )
       r24 ( : ) = particle_set ( index_b ) % r ( : ) - particle_set ( index_d ) % r ( : )
       r34 ( : ) = particle_set ( index_c ) % r ( : ) - particle_set ( index_d ) % r ( : )
       atomic_kind => particle_set ( index_a ) % atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
       imass1 = 1.0_dp/mass
       atomic_kind => particle_set ( index_b ) % atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
       imass2 = 1.0_dp/mass
       atomic_kind => particle_set ( index_c ) % atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
       imass3 = 1.0_dp/mass
       atomic_kind => particle_set ( index_d ) % atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
       imass4 = 1.0_dp/mass
       lg4x6 ( iconst ) % fa = -2.0_dp* r12
       lg4x6 ( iconst ) % fb = -2.0_dp* r13
       lg4x6 ( iconst ) % fc = -2.0_dp* r14
       lg4x6 ( iconst ) % fd = -2.0_dp* r23
       lg4x6 ( iconst ) % fe = -2.0_dp* r24
       lg4x6 ( iconst ) % ff = -2.0_dp* r34
       !
       ! construct matrix
       amat(1,1)=(imass1+imass2)*DOTPROD_3D(r12,lg4x6 ( iconst ) % fa )
       amat(1,2)=imass1*DOTPROD_3D(r12, lg4x6 ( iconst) % fb )
       amat(1,3)=imass1*DOTPROD_3D(r12, lg4x6 ( iconst ) % fc )
       amat(1,4)=-imass2*DOTPROD_3D(r12, lg4x6 ( iconst ) % fd )
       amat(1,5)=-imass2*DOTPROD_3D(r12, lg4x6 ( iconst ) % fe )
       amat(1,6)=0.0_dp

       amat(2,1)=imass1*DOTPROD_3D(r13,lg4x6 ( iconst ) % fa )
       amat(2,2)=(imass1+imass3)*DOTPROD_3D(r13, lg4x6 ( iconst) % fb )
       amat(2,3)=imass1*DOTPROD_3D(r13, lg4x6 ( iconst ) % fc )
       amat(2,4)=imass3*DOTPROD_3D(r13, lg4x6 ( iconst ) % fd )
       amat(2,5)=0.0_dp
       amat(2,6)=-imass3*DOTPROD_3D(r13, lg4x6 ( iconst ) % ff )

       amat(3,1)=imass1*DOTPROD_3D(r14,lg4x6 ( iconst ) % fa )
       amat(3,2)=imass1*DOTPROD_3D(r14, lg4x6 ( iconst) % fb )
       amat(3,3)=(imass1+imass4)*DOTPROD_3D(r14, lg4x6 ( iconst ) % fc )
       amat(3,4)=0.0_dp
       amat(3,5)=imass4*DOTPROD_3D(r14, lg4x6 ( iconst ) % fe )
       amat(3,6)=imass4*DOTPROD_3D(r14, lg4x6 ( iconst ) % ff )

       amat(4,1)=-imass2*DOTPROD_3D(r23,lg4x6 ( iconst ) % fa )
       amat(4,2)=imass3*DOTPROD_3D(r23, lg4x6 ( iconst) % fb )
       amat(4,3)=0.0_dp
       amat(4,4)=(imass3+imass2)*DOTPROD_3D(r23, lg4x6 ( iconst ) % fd )
       amat(4,5)=imass2*DOTPROD_3D(r23, lg4x6 ( iconst ) % fe )
       amat(4,6)=-imass3*DOTPROD_3D(r23, lg4x6 ( iconst ) % ff )

       amat(5,1)=-imass2*DOTPROD_3D(r24,lg4x6 ( iconst ) % fa )
       amat(5,2)=0.0_dp
       amat(5,3)=imass4*DOTPROD_3D(r24, lg4x6 ( iconst) % fc )
       amat(5,4)=imass2*DOTPROD_3D(r24, lg4x6 ( iconst ) % fd )
       amat(5,5)=(imass4+imass2)*DOTPROD_3D(r24, lg4x6 ( iconst ) % fe )
       amat(5,6)=imass4*DOTPROD_3D(r24, lg4x6 ( iconst ) % ff )

       amat(6,1)=0.0_dp
       amat(6,2)=-imass3*DOTPROD_3D(r34,lg4x6 ( iconst ) % fb )
       amat(6,3)=imass4*DOTPROD_3D(r34, lg4x6 ( iconst) % fc )
       amat(6,4)=-imass3*DOTPROD_3D(r34, lg4x6 ( iconst ) % fd )
       amat(6,5)=imass4*DOTPROD_3D(r34, lg4x6 ( iconst ) % fe )
       amat(6,6)=(imass3+imass4)*DOTPROD_3D(r34, lg4x6 ( iconst ) % ff )

       ! construct solution vector
       bvec(1,1)=DOTPROD_3D(r12,v12)
       bvec(2,1)=DOTPROD_3D(r13,v13)
       bvec(3,1)=DOTPROD_3D(r14,v14)
       bvec(4,1)=DOTPROD_3D(r23,v23)
       bvec(5,1)=DOTPROD_3D(r24,v24)
       bvec(6,1)=DOTPROD_3D(r34,v34)
       bvec=-bvec*2.0_dp*idt

       ! get lambda
       CALL solve_system ( amat, 6, bvec )
       lg4x6 ( iconst ) % lambda ( : )=bvec ( :, 1 )

       fc1=lg4x6 ( iconst ) % lambda ( 1 ) * lg4x6 ( iconst ) % fa + &
            lg4x6 ( iconst ) % lambda ( 2 ) * lg4x6 ( iconst ) % fb + &
            lg4x6 ( iconst ) % lambda ( 3 ) * lg4x6 ( iconst ) % fc 
       fc2=-lg4x6 ( iconst ) % lambda ( 1 ) * lg4x6 ( iconst ) % fa + &
            lg4x6 ( iconst ) % lambda ( 4 ) * lg4x6 ( iconst ) % fd + &
            lg4x6 ( iconst ) % lambda ( 5 ) * lg4x6 ( iconst ) % fe 
       fc3=-lg4x6 ( iconst ) % lambda ( 2 ) * lg4x6 ( iconst ) % fb - &
            lg4x6 ( iconst ) % lambda ( 4 ) * lg4x6 ( iconst ) % fd + &
            lg4x6 ( iconst ) % lambda ( 6 ) * lg4x6 ( iconst ) % ff 
       fc4=-lg4x6 ( iconst ) % lambda ( 3 ) * lg4x6 ( iconst ) % fc - &
            lg4x6 ( iconst ) % lambda ( 5 ) * lg4x6 ( iconst ) % fe - &
            lg4x6 ( iconst ) % lambda ( 6 ) * lg4x6 ( iconst ) % ff 
       vel ( :, index_a ) = vel ( :, index_a ) + imass1 * dtby2 * fc1 ( : )
       vel ( :, index_b ) = vel ( :, index_b ) + imass2 * dtby2 * fc2 ( : )
       vel ( :, index_c ) = vel ( :, index_c ) + imass3 * dtby2 * fc3 ( : )
       vel ( :, index_d ) = vel ( :, index_d ) + imass4 * dtby2 * fc4 ( : )
    END DO

  END SUBROUTINE rattle_4x6

!!****** constraint_4x6/rattle_roll_4x6 [1.0] *
!!
!!   NAME
!!     rattle_roll_4x6
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
  SUBROUTINE rattle_roll_4x6( molecule, particle_set, vel, r_rattle, simpar, &
       veps )
    TYPE(molecule_type), POINTER             :: molecule
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: vel( :, : )
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: r_rattle
    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: veps

    INTEGER                                  :: first_atom, iconst, index_a, &
                                                index_b, index_c, index_d, &
                                                ng4x6
    REAL(KIND=dp)                            :: dtby2, idt, imass1, imass2, &
                                                imass3, imass4, mass
    REAL(KIND=dp), DIMENSION(3) :: f_roll1, f_roll2, f_roll3, f_roll4, &
      f_roll5, f_roll6, fc1, fc2, fc3, fc4, r12, r13, r14, r23, r24, r34, &
      v12, v13, v14, v23, v24, v34
    REAL(KIND=dp), DIMENSION(6)              :: lambda
    REAL(KIND=dp), DIMENSION(6, 1)           :: bvec
    REAL(KIND=dp), DIMENSION(6, 6)           :: amat
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(g4x6_constraint_type), POINTER      :: g4x6_list( : )
    TYPE(local_g4x6_constraint_type), &
      POINTER                                :: lg4x6( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

    idt = 1.0_dp/simpar%dt
    dtby2 = simpar%dt*.5_dp
    molecule_kind => molecule % molecule_kind
    CALL get_molecule_kind ( molecule_kind, ng4x6 = ng4x6, g4x6_list = g4x6_list )
    CALL get_molecule ( molecule, first_atom = first_atom, lg4x6=lg4x6 )
    DO iconst = 1, SIZE ( g4x6_list )
       index_a = g4x6_list ( iconst ) % a + first_atom -1
       index_b = g4x6_list ( iconst ) % b + first_atom -1
       index_c = g4x6_list ( iconst ) % c + first_atom -1
       index_d = g4x6_list ( iconst ) % d + first_atom -1
       v12 ( : ) = vel ( :, index_a ) - vel ( :, index_b ) 
       v13 ( : ) = vel ( :, index_a ) - vel ( :, index_c ) 
       v14 ( : ) = vel ( :, index_a ) - vel ( :, index_d ) 
       v23 ( : ) = vel ( :, index_b ) - vel ( :, index_c ) 
       v24 ( : ) = vel ( :, index_b ) - vel ( :, index_d ) 
       v34 ( : ) = vel ( :, index_c ) - vel ( :, index_d ) 

       r12 ( : ) = particle_set ( index_a ) % r ( : ) - particle_set ( index_b ) % r ( : )
       r13 ( : ) = particle_set ( index_a ) % r ( : ) - particle_set ( index_c ) % r ( : )
       r14 ( : ) = particle_set ( index_a ) % r ( : ) - particle_set ( index_d ) % r ( : )
       r23 ( : ) = particle_set ( index_b ) % r ( : ) - particle_set ( index_c ) % r ( : )
       r24 ( : ) = particle_set ( index_b ) % r ( : ) - particle_set ( index_d ) % r ( : )
       r34 ( : ) = particle_set ( index_c ) % r ( : ) - particle_set ( index_d ) % r ( : )
       atomic_kind => particle_set ( index_a ) % atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
       imass1 = 1.0_dp/mass
       atomic_kind => particle_set ( index_b ) % atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
       imass2 = 1.0_dp/mass
       atomic_kind => particle_set ( index_c ) % atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
       imass3 = 1.0_dp/mass
       atomic_kind => particle_set ( index_d ) % atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
       imass4 = 1.0_dp/mass
       lg4x6 ( iconst ) % fa = -2.0_dp* r12
       lg4x6 ( iconst ) % fb = -2.0_dp* r13
       lg4x6 ( iconst ) % fc = -2.0_dp* r14
       lg4x6 ( iconst ) % fd = -2.0_dp* r23
       lg4x6 ( iconst ) % fe = -2.0_dp* r24
       lg4x6 ( iconst ) % ff = -2.0_dp* r34
       ! roll the fc
       f_roll1 = MATVEC_3x3 ( r_rattle, lg4x6 ( iconst ) % fa )
       f_roll2 = MATVEC_3x3 ( r_rattle, lg4x6 ( iconst ) % fb )
       f_roll3 = MATVEC_3x3 ( r_rattle, lg4x6 ( iconst ) % fc )
       f_roll4 = MATVEC_3x3 ( r_rattle, lg4x6 ( iconst ) % fd )
       f_roll5 = MATVEC_3x3 ( r_rattle, lg4x6 ( iconst ) % fe )
       f_roll6 = MATVEC_3x3 ( r_rattle, lg4x6 ( iconst ) % ff )
       ! construct matrix
       amat(1,1)=(imass1+imass2)*DOTPROD_3D(r12,f_roll1 )
       amat(1,2)=imass1*DOTPROD_3D(r12, f_roll2 )
       amat(1,3)=imass1*DOTPROD_3D(r12, f_roll3 )
       amat(1,4)=-imass2*DOTPROD_3D(r12, f_roll4 )
       amat(1,5)=-imass2*DOTPROD_3D(r12, f_roll5 )
       amat(1,6)=0.0_dp

       amat(2,1)=imass1*DOTPROD_3D(r13,f_roll1 )
       amat(2,2)=(imass1+imass3)*DOTPROD_3D(r13, f_roll2)
       amat(2,3)=imass1*DOTPROD_3D(r13, f_roll3 )
       amat(2,4)=imass3*DOTPROD_3D(r13, f_roll4 )
       amat(2,5)=0.0_dp
       amat(2,6)=-imass3*DOTPROD_3D(r13, f_roll6 )

       amat(3,1)=imass1*DOTPROD_3D(r14,f_roll1 )
       amat(3,2)=imass1*DOTPROD_3D(r14, f_roll2 )
       amat(3,3)=(imass1+imass4)*DOTPROD_3D(r14, f_roll3 )
       amat(3,4)=0.0_dp
       amat(3,5)=imass4*DOTPROD_3D(r14, f_roll5 )
       amat(3,6)=imass4*DOTPROD_3D(r14, f_roll6 )

       amat(4,1)=-imass2*DOTPROD_3D(r23,f_roll1 )
       amat(4,2)=imass3*DOTPROD_3D(r23, f_roll2 )
       amat(4,3)=0.0_dp
       amat(4,4)=(imass3+imass2)*DOTPROD_3D(r23, f_roll4 )
       amat(4,5)=imass2*DOTPROD_3D(r23, f_roll5 )
       amat(4,6)=-imass3*DOTPROD_3D(r23, f_roll6 )

       amat(5,1)=-imass2*DOTPROD_3D(r24,f_roll1 )
       amat(5,2)=0.0_dp
       amat(5,3)=imass4*DOTPROD_3D(r24, f_roll3)
       amat(5,4)=imass2*DOTPROD_3D(r24, f_roll4 )
       amat(5,5)=(imass4+imass2)*DOTPROD_3D(r24, f_roll5 )
       amat(5,6)=imass4*DOTPROD_3D(r24, f_roll6 )

       amat(6,1)=0.0_dp
       amat(6,2)=-imass3*DOTPROD_3D(r34,f_roll2 )
       amat(6,3)=imass4*DOTPROD_3D(r34, f_roll3 )
       amat(6,4)=-imass3*DOTPROD_3D(r34, f_roll4 )
       amat(6,5)=imass4*DOTPROD_3D(r34, f_roll5 )
       amat(6,6)=(imass3+imass4)*DOTPROD_3D(r34, f_roll6 )

       ! construct solution vector
       bvec ( 1, 1 ) = DOTPROD_3D ( r12, v12 + MATVEC_3X3 ( veps, r12 ) )
       bvec ( 2, 1 ) = DOTPROD_3D ( r13, v13 + MATVEC_3X3 ( veps, r13 ) )
       bvec ( 3, 1 ) = DOTPROD_3D ( r14, v14 + MATVEC_3X3 ( veps, r14 ) )
       bvec ( 4, 1 ) = DOTPROD_3D ( r23, v23 + MATVEC_3X3 ( veps, r23 ) )
       bvec ( 5, 1 ) = DOTPROD_3D ( r24, v24 + MATVEC_3X3 ( veps, r24 ) )
       bvec ( 6, 1 ) = DOTPROD_3D ( r34, v34 + MATVEC_3X3 ( veps, r34 ) )
       bvec=-bvec*2.0_dp*idt


       ! get lambda
       CALL solve_system ( amat, 6, bvec )
       lambda ( : ) = bvec ( :, 1 )
       lg4x6 ( iconst ) % lambda ( : )= lambda 

       fc1=lambda ( 1 ) * f_roll1 + &
            lambda ( 2 ) * f_roll2 + &
            lambda ( 3 ) * f_roll3 
       fc2=-lambda ( 1 ) * f_roll1 + &
            lambda ( 4 ) * f_roll4 + &
            lambda ( 5 ) * f_roll5 
       fc3=-lambda ( 2 ) * f_roll2 - &
            lambda ( 4 ) * f_roll4 + &
            lambda ( 6 ) * f_roll6 
       fc4=-lambda ( 3 ) * f_roll3 - &
            lambda ( 5 ) * f_roll5 - &
            lambda ( 6 ) * f_roll6 
       vel ( :, index_a ) = vel ( :, index_a ) + imass1 * dtby2 * fc1 ( : )
       vel ( :, index_b ) = vel ( :, index_b ) + imass2 * dtby2 * fc2 ( : )
       vel ( :, index_c ) = vel ( :, index_c ) + imass3 * dtby2 * fc3 ( : )
       vel ( :, index_d ) = vel ( :, index_d ) + imass4 * dtby2 * fc4 ( : )
    END DO

  END SUBROUTINE rattle_roll_4x6

END MODULE constraint_4x6
