!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****** cp2k/constraint [1.0] *
!!
!!   NAME
!!     constraint
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
MODULE constraint
  USE cell_types,                      ONLY: cell_type
  USE constraint_3x3,                  ONLY: rattle_3x3,&
                                             rattle_roll_3x3,&
                                             shake_3x3,&
                                             shake_roll_3x3
  USE constraint_4x6,                  ONLY: rattle_4x6,&
                                             rattle_roll_4x6,&
                                             shake_4x6,&
                                             shake_roll_4x6
  USE constraint_clv,                  ONLY: rattle,&
                                             rattle_roll,&
                                             shake,&
                                             shake_roll
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE input_constants,                 ONLY: npt_f_ensemble,&
                                             npt_i_ensemble
  USE kinds,                           ONLY: dp
  USE md,                              ONLY: simulation_parameters_type
  USE message_passing,                 ONLY: mp_sum
  USE molecule_kind_types,             ONLY: colvar_constraint_type,&
                                             g3x3_constraint_type,&
                                             g4x6_constraint_type,&
                                             get_molecule_kind,&
                                             molecule_kind_type,&
                                             colvar_counters
  USE molecule_types_new,              ONLY: get_molecule,&
                                             local_colvar_constraint_type,&
                                             local_g3x3_constraint_type,&
                                             local_g4x6_constraint_type,&
                                             molecule_type
  USE particle_types,                  ONLY: particle_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: matmul_3x3,&
                                             transpose_3d
  USE virial_types,                    ONLY: virial_type
  USE colvar_types, ONLY: dist_colvar_id
  USE colvar_methods, ONLY: colvar_eval_mol_f
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: shake_control, rattle_control,  getold,  pv_constraint
  PUBLIC :: shake_roll_control, rattle_roll_control
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'constraint'

CONTAINS
!!****** constraint/shake_control [1.0] *
!!
!!   NAME
!!     shake_control
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
  SUBROUTINE shake_control( local_molecules, molecule_set, molecule_kind_set, &
       particle_set, pos, vel, simpar )

    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: pos( :, : ), vel( :, : )
    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar

    CHARACTER(LEN=*), PARAMETER              :: routineN = 'shake_control', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ikind, imol, &
                                                ishake, n3x3con, n4x6con, &
                                                nkind, nmol_per_kind
    TYPE(colvar_counters)                    :: ncolv
    REAL(KIND=dp)                            :: max_sigma
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule

    CALL timeset ( routineN, 'I', ' ', handle )
    nkind = SIZE ( molecule_kind_set )

    MOL:  DO ikind = 1, nkind
       nmol_per_kind = local_molecules % n_el ( ikind )
       DO imol = 1, nmol_per_kind
          i = local_molecules % list ( ikind ) % array ( imol )
          molecule => molecule_set ( i ) 
          molecule_kind => molecule % molecule_kind
          CALL get_molecule_kind ( molecule_kind, ncolv = ncolv,  &
               ng3x3 = n3x3con, ng4x6=n4x6con )

          ! 3x3
          IF ( n3x3con /= 0 ) CALL shake_3x3( molecule, particle_set,  &
               pos, vel, simpar, ishake )

          ! 4x6
          IF ( n4x6con /= 0 ) CALL shake_4x6( molecule, particle_set, &
               pos, vel, simpar, ishake )

          !
          ! collective variables
          !
          ishake = 0       
          IF ( ncolv%ndist /= 0 ) THEN
             max_sigma = -1.0E+10_dp
             CONS: DO WHILE (ABS(max_sigma)>=simpar%shake_tol)
                max_sigma = 0.0_dp
                ishake = ishake + 1
                CALL shake( molecule, particle_set, pos, vel, simpar, &
                     max_sigma,ishake )
             END DO CONS
          END IF
       END DO
    END DO MOL

    CALL timestop ( 0.0_dp, handle )

  END SUBROUTINE shake_control

!!****** constraint/rattle_control [1.0] *
!!
!!   NAME
!!     rattle_control
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
  SUBROUTINE rattle_control( local_molecules, molecule_set, molecule_kind_set,  &
       particle_set, vel, simpar )

    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: vel( :, : )
    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar

    CHARACTER(LEN=*), PARAMETER              :: routineN = 'rattle_control', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ikind, imol, &
                                                irattle, n3x3con, n4x6con, &
                                                 nkind, nmol_per_kind
    TYPE(colvar_counters)                    :: ncolv
    REAL(KIND=dp)                            :: max_sigma
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule

    CALL timeset ( routineN, 'I', ' ', handle )
    nkind = SIZE ( molecule_kind_set )
    MOL: DO ikind = 1, nkind
       nmol_per_kind = local_molecules % n_el ( ikind )
       DO imol = 1, nmol_per_kind
          i = local_molecules % list ( ikind ) % array ( imol )
          molecule => molecule_set ( i ) 
          molecule_kind => molecule % molecule_kind
          CALL get_molecule_kind ( molecule_kind, ncolv = ncolv, ng3x3 = n3x3con, ng4x6 = n4x6con )

          ! 3x3
          IF ( n3x3con /= 0 ) CALL rattle_3x3 ( molecule, particle_set, vel, simpar )

          ! 4x6
          IF ( n4x6con /= 0 ) CALL rattle_4x6 ( molecule, particle_set, vel, simpar )

          irattle = 0
          ! dist
          IF ( ncolv%ndist /= 0 ) THEN
             max_sigma = -1.0E+10_dp
             CONS: DO WHILE ( ABS ( max_sigma ) >= simpar % shake_tol )
                max_sigma = 0.0_dp
                irattle = irattle + 1
                CALL rattle ( molecule, particle_set, vel, simpar, max_sigma, irattle )
             END DO CONS
          END IF
       END DO
    END DO MOL
    CALL timestop ( 0.0_dp, handle )

  END SUBROUTINE rattle_control

!!****** constraint/shake_roll_control [1.0] *
!!
!!   NAME
!!     shake_roll_control
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
  SUBROUTINE shake_roll_control( local_molecules, molecule_set, molecule_kind_set, &
       particle_set, pos, vel, simpar, roll_tol, iroll,  &
       vector_r, vector_v, group, u )

    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: pos( :, : ), vel( :, : )
    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar
    REAL(KIND=dp), INTENT(OUT)               :: roll_tol
    INTEGER, INTENT(INOUT)                   :: iroll
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: vector_r, vector_v
    INTEGER, INTENT(IN)                      :: group
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN), OPTIONAL                   :: u

    CHARACTER(LEN=*), PARAMETER :: routineN = 'shake_roll_control', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ikind, imol, &
                                                ishake, n3x3con, n4x6con, &
                                                nkind, nmol_per_kind
    TYPE(colvar_counters)                    :: ncolv
    REAL(KIND=dp)                            :: max_sigma
    REAL(KIND=dp), DIMENSION(3, 3)           :: r_shake, v_shake
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule

    CALL timeset ( routineN, 'I', ' ', handle )
    nkind = SIZE ( molecule_kind_set )
    ! setting up for roll
    IF ( simpar % ensemble == npt_i_ensemble ) THEN
       CALL get_roll_matrix ('SHAKE', r_shake, v_shake, vector_r, vector_v )
    ELSE IF ( simpar % ensemble == npt_f_ensemble ) THEN
       CALL get_roll_matrix ('SHAKE', r_shake, v_shake, vector_r, vector_v, u )
    END IF

    MOL:  DO ikind = 1, nkind
       nmol_per_kind = local_molecules % n_el ( ikind )
       DO imol = 1, nmol_per_kind
          i = local_molecules % list ( ikind ) % array ( imol )
          molecule => molecule_set ( i ) 
          molecule_kind => molecule % molecule_kind
          CALL get_molecule_kind ( molecule_kind, ncolv = ncolv,  &
               ng3x3 = n3x3con, ng4x6=n4x6con )

          ! 3x3
          IF ( n3x3con /= 0 ) CALL shake_roll_3x3( molecule, particle_set,  &
               pos, vel, r_shake, v_shake, &
               simpar, ishake )

          ! 4x6
          IF ( n4x6con /= 0 ) CALL shake_roll_4x6( molecule, particle_set, &
               pos, vel, r_shake, v_shake, &
               simpar, ishake )

          ishake = 0

          ! dist
          IF ( ncolv%ndist /= 0 ) THEN
             max_sigma = -1.0E+10_dp
             CONS: DO WHILE (ABS(max_sigma)>=simpar%shake_tol)
                max_sigma = 0.0_dp
                ishake = ishake + 1
                CALL shake_roll( molecule, particle_set, pos, vel, &
                     r_shake, v_shake, simpar, max_sigma,ishake )
             END DO CONS
          END IF
       END DO
    END DO MOL

    CALL check_tol ( roll_tol, iroll, 'SHAKE', r_shake )

    CALL timestop ( 0.0_dp, handle )

  END SUBROUTINE shake_roll_control

!!****** constraint/rattle_roll_control [1.0] *
!!
!!   NAME
!!     rattle_roll_control
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
  SUBROUTINE rattle_roll_control( local_molecules, molecule_set, molecule_kind_set, &
       particle_set, vel, simpar, virial, vector,   &
       veps, mass, roll_tol, iroll, para_env, u )

    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: vel( :, : )
    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar
    TYPE(virial_type), POINTER               :: virial
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: vector
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: veps
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: mass
    REAL(KIND=dp), INTENT(OUT)               :: roll_tol
    INTEGER, INTENT(INOUT)                   :: iroll
    TYPE(cp_para_env_type), INTENT(IN)       :: para_env
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN), OPTIONAL                   :: u

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rattle_roll_control', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ikind, imol, &
                                                irattle, n3x3con, n4x6con, &
                                                nkind, nmol_per_kind
    TYPE(colvar_counters)                    :: ncolv
    REAL(KIND=dp)                            :: max_sigma
    REAL(KIND=dp), DIMENSION(3, 3)           :: r_rattle
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule

    CALL timeset ( routineN, 'I', ' ', handle )
    ! initialize locals
    nkind = SIZE ( molecule_kind_set )
    ! setting up for roll
    IF ( simpar % ensemble == npt_i_ensemble ) THEN
       CALL get_roll_matrix ('RATTLE', v_shake=r_rattle, vector_v=vector )
    ELSE IF ( simpar % ensemble == npt_f_ensemble ) THEN
       CALL get_roll_matrix ('RATTLE', v_shake=r_rattle, vector_v=vector, u=u )
    END IF
    MOL: DO ikind = 1, nkind
       nmol_per_kind = local_molecules % n_el ( ikind )
       DO imol = 1, nmol_per_kind
          i = local_molecules % list ( ikind ) % array ( imol )
          molecule => molecule_set ( i ) 
          molecule_kind => molecule % molecule_kind
          CALL get_molecule_kind ( molecule_kind, ncolv = ncolv,  &
               ng3x3 = n3x3con, ng4x6 = n4x6con )
          ! 3x3
          IF ( n3x3con /= 0 ) CALL rattle_roll_3x3 ( molecule, particle_set, vel, r_rattle, &
               simpar, veps )

          ! 4x6
          IF ( n4x6con /= 0 ) CALL rattle_roll_4x6 ( molecule, particle_set, vel, r_rattle, &
               simpar, veps )
          irattle = 0 
          ! dist
          IF ( ncolv%ndist /= 0 )  THEN
             max_sigma = -1.0E+10_dp
             CONS: DO WHILE (ABS(max_sigma)>=simpar%shake_tol)
                max_sigma = 0.0_dp
                irattle = irattle + 1
                CALL rattle_roll ( molecule, particle_set, vel, r_rattle, &
                     simpar, veps, max_sigma, irattle )
             END DO CONS
          ENDIF
       END DO
    END DO MOL


    CALL check_tol ( roll_tol, iroll, 'RATTLE', veps = veps )

    CALL timestop ( 0.0_dp, handle )

  END SUBROUTINE rattle_roll_control

!!****** constraint/getold [1.0] *
!!
!!   NAME
!!     getold
!!
!!   FUNCTION
!!     saves all of the old variables
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE getold( local_molecules, molecule_set, molecule_kind_set, particle_set ,&
       cell, error)

    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'getold', &
      routineP = moduleN//':'//routineN
    INTEGER                                  :: first_atom, i, iconst, ikind, &
                                                imol, index, n3x3con, &
                                                n4x6con, nkind, &
                                                nmol_per_kind
    TYPE(colvar_constraint_type), POINTER      :: colv_list( : )
    TYPE(g3x3_constraint_type), POINTER      :: g3x3_list( : )
    TYPE(g4x6_constraint_type), POINTER      :: g4x6_list( : )
    TYPE(local_colvar_constraint_type), &
      POINTER                                :: lcolv( : )
    TYPE(local_g3x3_constraint_type), &
      POINTER                                :: lg3x3( : )
    TYPE(local_g4x6_constraint_type), &
      POINTER                                :: lg4x6( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule

    nkind = SIZE ( molecule_kind_set )
    DO ikind = 1, nkind
       nmol_per_kind = local_molecules % n_el ( ikind )
       DO imol = 1, nmol_per_kind
          i = local_molecules % list ( ikind ) % array ( imol )
          molecule => molecule_set ( i )
          molecule_kind => molecule % molecule_kind
          CALL get_molecule_kind ( molecule_kind, ng3x3=n3x3con, ng4x6=n4x6con, & 
               colv_list=colv_list, g3x3_list=g3x3_list, g4x6_list=g4x6_list )
          CALL get_molecule ( molecule, first_atom=first_atom, &
               lcolv=lcolv, lg3x3=lg3x3, lg4x6=lg4x6 )
          ! distance constraints
          IF (ASSOCIATED(colv_list)) THEN
             DO iconst = 1, SIZE(colv_list)
                CALL colvar_eval_mol_f(lcolv ( iconst ) % colvar, cell,&
                     particles=particle_set, error=error)
                IF ( colv_list ( iconst ) % type_id == dist_colvar_id ) THEN
                   index = colv_list ( iconst ) % i_atoms(1) + first_atom - 1
                   lcolv ( iconst ) % ra_old = particle_set ( index ) % r 
                   index = colv_list ( iconst ) % i_atoms(2) + first_atom - 1
                   lcolv ( iconst ) % rb_old = particle_set ( index ) % r 
                END IF
             ENDDO
          END IF
          ! 3x3 constraints
          DO iconst = 1, n3x3con
             index = g3x3_list ( iconst ) % a + first_atom - 1
             lg3x3 ( iconst ) % ra_old = particle_set ( index ) % r 
             index = g3x3_list ( iconst ) % b + first_atom - 1
             lg3x3 ( iconst ) % rb_old = particle_set ( index ) % r 
             index = g3x3_list ( iconst ) % c + first_atom - 1
             lg3x3 ( iconst ) % rc_old = particle_set ( index ) % r 
          ENDDO
          ! 4x6 constraints
          DO iconst = 1, n4x6con
             index = g4x6_list ( iconst ) % a + first_atom - 1
             lg4x6 ( iconst ) % ra_old = particle_set ( index ) % r 
             index = g4x6_list ( iconst ) % b + first_atom - 1
             lg4x6 ( iconst ) % rb_old = particle_set ( index ) % r 
             index = g4x6_list ( iconst ) % c + first_atom - 1
             lg4x6 ( iconst ) % rc_old = particle_set ( index ) % r 
             index = g4x6_list ( iconst ) % d + first_atom - 1
             lg4x6 ( iconst ) % rd_old = particle_set ( index ) % r 
          ENDDO
       END DO
    END DO

  END SUBROUTINE getold

!!****** constraint/pv_constraint [1.0] *
!!
!!   NAME
!!     pv_constraint
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
  SUBROUTINE pv_constraint ( local_molecules, molecule_set, molecule_kind_set, &
       particle_set, virial, group )

    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(virial_type), INTENT(INOUT)         :: virial
    INTEGER, INTENT(IN)                      :: group

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pv_constraint', &
      routineP = moduleN//':'//routineN
    INTEGER :: first_atom, i, iconst, ikind, imol, index_a, index_b, index_c, &
      index_d, ng3x3, ng4x6, nkind, nmol_per_kind
    REAL(KIND=dp) :: fc1( 3 ), fc2( 3 ), fc3( 3 ), fc4( 3 ), lambda, &
      lambda_3x3( 3 ), lambda_4x6( 6 ), pv( 3, 3 )
    TYPE(colvar_constraint_type), POINTER      :: colv_list( : )
    TYPE(g3x3_constraint_type), POINTER      :: g3x3_list( : )
    TYPE(g4x6_constraint_type), POINTER      :: g4x6_list( : )
    TYPE(local_colvar_constraint_type), &
      POINTER                                :: lcolv( : )
    TYPE(local_g3x3_constraint_type), &
      POINTER                                :: lg3x3( : )
    TYPE(local_g4x6_constraint_type), &
      POINTER                                :: lg4x6( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule

    pv = 0.0_dp
    nkind = SIZE ( molecule_kind_set )
    MOL: DO ikind = 1, nkind
       nmol_per_kind = local_molecules % n_el ( ikind )
       DO imol = 1, nmol_per_kind
          i = local_molecules % list ( ikind ) % array ( imol )
          molecule => molecule_set ( i ) 
          molecule_kind => molecule % molecule_kind
          CALL get_molecule_kind ( molecule_kind, ng3x3 = ng3x3, &
               ng4x6 = ng4x6, g3x3_list=g3x3_list, g4x6_list=g4x6_list, &
               colv_list = colv_list )
          CALL get_molecule ( molecule, first_atom = first_atom, lg3x3=lg3x3, &
               lg4x6 = lg4x6, lcolv = lcolv )
          ! dist
          IF (ASSOCIATED(colv_list)) THEN
             DO iconst = 1, SIZE(colv_list)
                IF ( colv_list ( iconst ) % type_id == dist_colvar_id ) THEN
                   !  pv gets updated with FULL multiplier
                   
                   lambda = lcolv ( iconst ) % lambda
                   index_a = colv_list ( iconst ) % i_atoms(1) + first_atom -1
                   index_b = colv_list ( iconst ) % i_atoms(2) + first_atom -1
                   !pv(1,1)
                   pv(1,1) = pv(1,1) + lambda*lcolv ( iconst ) %f(1)*particle_set ( index_a ) % r ( 1 )
                   pv(1,1) = pv(1,1) - lambda*lcolv ( iconst ) %f(1)*particle_set ( index_b ) % r ( 1 )
                   !pv(1,2)
                   pv(1,2) = pv(1,2) + lambda*lcolv ( iconst ) %f(1)*particle_set ( index_a ) % r ( 2 )
                   pv(1,2) = pv(1,2) - lambda*lcolv ( iconst ) %f(1)*particle_set ( index_b ) % r ( 2 )
                   !pv(1,3)
                   pv(1,3) = pv(1,3) + lambda*lcolv ( iconst ) %f(1)*particle_set ( index_a ) % r ( 3 )
                   pv(1,3) = pv(1,3) - lambda*lcolv ( iconst ) %f(1)*particle_set ( index_b ) % r ( 3 )
                   !pv(2,1)
                   pv(2,1) = pv(2,1) + lambda*lcolv ( iconst ) %f(2)*particle_set ( index_a ) % r ( 1 )
                   pv(2,1) = pv(2,1) - lambda*lcolv ( iconst ) %f(2)*particle_set ( index_b ) % r ( 1 )
                   !pv(2,2)
                   pv(2,2) = pv(2,2) + lambda*lcolv ( iconst ) %f(2)*particle_set ( index_a ) % r ( 2 )
                   pv(2,2) = pv(2,2) - lambda*lcolv ( iconst ) %f(2)*particle_set ( index_b ) % r ( 2 )
                   !pv(2,3)
                   pv(2,3) = pv(2,3) + lambda*lcolv ( iconst ) %f(2)*particle_set ( index_a ) % r ( 3 )
                   pv(2,3) = pv(2,3) - lambda*lcolv ( iconst ) %f(2)*particle_set ( index_b ) % r ( 3 )
                   !pv(3,1)
                   pv(3,1) = pv(3,1) + lambda*lcolv ( iconst ) %f(3)*particle_set ( index_a ) % r ( 1 )
                   pv(3,1) = pv(3,1) - lambda*lcolv ( iconst ) %f(3)*particle_set ( index_b ) % r ( 1 )
                   !pv(3,2)
                   pv(3,2) = pv(3,2) + lambda*lcolv ( iconst ) %f(3)*particle_set ( index_a ) % r ( 2 )
                   pv(3,2) = pv(3,2) - lambda*lcolv ( iconst ) %f(3)*particle_set ( index_b ) % r ( 2 )
                   !pv(3,3)
                   pv(3,3) = pv(3,3) + lambda*lcolv ( iconst ) %f(3)*particle_set ( index_a ) % r ( 3 )
                   pv(3,3) = pv(3,3) - lambda*lcolv ( iconst ) %f(3)*particle_set ( index_b ) % r ( 3 )
                END IF
             END DO
          END IF
          ! 3x3
          DO iconst = 1, ng3x3

             !  pv gets updated with FULL multiplier
             lambda_3x3 = lg3x3 ( iconst ) % lambda


             fc1=lambda_3x3(1)*lg3x3 ( iconst ) %fa+ &
                  lambda_3x3(2)*lg3x3 ( iconst ) %fb
             fc2=-lambda_3x3(1)*lg3x3 ( iconst ) %fa+ &
                  lambda_3x3(3)*lg3x3 ( iconst ) %fc
             fc3=-lambda_3x3(2)*lg3x3 ( iconst ) %fb- &
                  lambda_3x3(3)*lg3x3 ( iconst ) %fc
             index_a = g3x3_list ( iconst ) % a + first_atom -1
             index_b = g3x3_list ( iconst ) % b + first_atom -1
             index_c = g3x3_list ( iconst ) % c + first_atom -1

             !pv(1,1)
             pv(1,1) = pv(1,1) + fc1(1)*particle_set ( index_a ) % r ( 1 )
             pv(1,1) = pv(1,1) + fc2(1)*particle_set ( index_b ) % r ( 1 )
             pv(1,1) = pv(1,1) + fc3(1)*particle_set ( index_c ) % r ( 1 )
             !pv(1,2)
             pv(1,2) = pv(1,2) + fc1(1)*particle_set ( index_a ) % r ( 2 )
             pv(1,2) = pv(1,2) + fc2(1)*particle_set ( index_b ) % r ( 2 )
             pv(1,2) = pv(1,2) + fc3(1)*particle_set ( index_c ) % r ( 2 )
             !pv(1,3)
             pv(1,3) = pv(1,3) + fc1(1)*particle_set ( index_a ) % r ( 3 )
             pv(1,3) = pv(1,3) + fc2(1)*particle_set ( index_b ) % r ( 3 )
             pv(1,3) = pv(1,3) + fc3(1)*particle_set ( index_c ) % r ( 3 )
             !pv(2,1)
             pv(2,1) = pv(2,1) + fc1(2)*particle_set ( index_a ) % r ( 1 )
             pv(2,1) = pv(2,1) + fc2(2)*particle_set ( index_b ) % r ( 1 )
             pv(2,1) = pv(2,1) + fc3(2)*particle_set ( index_c ) % r ( 1 )
             !pv(2,2)
             pv(2,2) = pv(2,2) + fc1(2)*particle_set ( index_a ) % r ( 2 )
             pv(2,2) = pv(2,2) + fc2(2)*particle_set ( index_b ) % r ( 2 )
             pv(2,2) = pv(2,2) + fc3(2)*particle_set ( index_c ) % r ( 2 )
             !pv(2,3)
             pv(2,3) = pv(2,3) + fc1(2)*particle_set ( index_a ) % r ( 3 )
             pv(2,3) = pv(2,3) + fc2(2)*particle_set ( index_b ) % r ( 3 )
             pv(2,3) = pv(2,3) + fc3(2)*particle_set ( index_c ) % r ( 3 )
             !pv(3,1)
             pv(3,1) = pv(3,1) + fc1(3)*particle_set ( index_a ) % r ( 1 )
             pv(3,1) = pv(3,1) + fc2(3)*particle_set ( index_b ) % r ( 1 )
             pv(3,1) = pv(3,1) + fc3(3)*particle_set ( index_c ) % r ( 1 )
             !pv(3,2)
             pv(3,2) = pv(3,2) + fc1(3)*particle_set ( index_a ) % r ( 2 )
             pv(3,2) = pv(3,2) + fc2(3)*particle_set ( index_b ) % r ( 2 )
             pv(3,2) = pv(3,2) + fc3(3)*particle_set ( index_c ) % r ( 2 )
             !pv(3,3)
             pv(3,3) = pv(3,3) + fc1(3)*particle_set ( index_a ) % r ( 3 )
             pv(3,3) = pv(3,3) + fc2(3)*particle_set ( index_b ) % r ( 3 )
             pv(3,3) = pv(3,3) + fc3(3)*particle_set ( index_c ) % r ( 3 )
          END DO

          ! 4x6
          DO iconst = 1, ng4x6

             !  pv gets updated with FULL multiplier
             lambda_4x6 = lg4x6 ( iconst ) % lambda

             fc1=lambda_4x6(1)*lg4x6 ( iconst ) %fa+ &
                  lambda_4x6(2)*lg4x6 ( iconst ) %fb+ &
                  lambda_4x6(3)*lg4x6 ( iconst ) %fc
             fc2=-lambda_4x6(1)*lg4x6 ( iconst ) %fa+ &
                  lambda_4x6(4)*lg4x6 ( iconst ) %fd+ &
                  lambda_4x6(5)*lg4x6 ( iconst ) %fe
             fc3=-lambda_4x6(2)*lg4x6 ( iconst ) %fb- &
                  lambda_4x6(4)*lg4x6 ( iconst ) %fd+ &
                  lambda_4x6(6)*lg4x6 ( iconst ) %ff
             fc4=-lambda_4x6(3)*lg4x6 ( iconst ) %fc- &
                  lambda_4x6(5)*lg4x6 ( iconst ) %fe- &
                  lambda_4x6(6)*lg4x6 ( iconst ) %ff
             index_a = g4x6_list ( iconst ) % a + first_atom -1
             index_b = g4x6_list ( iconst ) % b + first_atom -1
             index_c = g4x6_list ( iconst ) % c + first_atom -1
             index_d = g4x6_list ( iconst ) % d + first_atom -1

             !pv(1,1)
             pv(1,1) = pv(1,1) + fc1(1)*particle_set(index_a)%r(1)
             pv(1,1) = pv(1,1) + fc2(1)*particle_set(index_b)%r(1)
             pv(1,1) = pv(1,1) + fc3(1)*particle_set(index_c)%r(1)
             pv(1,1) = pv(1,1) + fc4(1)*particle_set(index_d)%r(1)
             !pv(1,2)
             pv(1,2) = pv(1,2) + fc1(1)*particle_set(index_a)%r(2)
             pv(1,2) = pv(1,2) + fc2(1)*particle_set(index_b)%r(2)
             pv(1,2) = pv(1,2) + fc3(1)*particle_set(index_c)%r(2)
             pv(1,2) = pv(1,2) + fc4(1)*particle_set(index_d)%r(2)
             !pv(1,3)
             pv(1,3) = pv(1,3) + fc1(1)*particle_set(index_a)%r(3)
             pv(1,3) = pv(1,3) + fc2(1)*particle_set(index_b)%r(3)
             pv(1,3) = pv(1,3) + fc3(1)*particle_set(index_c)%r(3)
             pv(1,3) = pv(1,3) + fc4(1)*particle_set(index_d)%r(3)
             !pv(2,1)
             pv(2,1) = pv(2,1) + fc1(2)*particle_set(index_a)%r(1)
             pv(2,1) = pv(2,1) + fc2(2)*particle_set(index_b)%r(1)
             pv(2,1) = pv(2,1) + fc3(2)*particle_set(index_c)%r(1)
             pv(2,1) = pv(2,1) + fc4(2)*particle_set(index_d)%r(1)
             !pv(2,2)
             pv(2,2) = pv(2,2) + fc1(2)*particle_set(index_a)%r(2)
             pv(2,2) = pv(2,2) + fc2(2)*particle_set(index_b)%r(2)
             pv(2,2) = pv(2,2) + fc3(2)*particle_set(index_c)%r(2)
             pv(2,2) = pv(2,2) + fc4(2)*particle_set(index_d)%r(2)
             !pv(2,3)
             pv(2,3) = pv(2,3) + fc1(2)*particle_set(index_a)%r(3)
             pv(2,3) = pv(2,3) + fc2(2)*particle_set(index_b)%r(3)
             pv(2,3) = pv(2,3) + fc3(2)*particle_set(index_c)%r(3)
             pv(2,3) = pv(2,3) + fc4(2)*particle_set(index_d)%r(3)
             !pv(3,1)
             pv(3,1) = pv(3,1) + fc1(3)*particle_set(index_a)%r(1)
             pv(3,1) = pv(3,1) + fc2(3)*particle_set(index_b)%r(1)
             pv(3,1) = pv(3,1) + fc3(3)*particle_set(index_c)%r(1)
             pv(3,1) = pv(3,1) + fc4(3)*particle_set(index_d)%r(1)
             !pv(3,2)
             pv(3,2) = pv(3,2) + fc1(3)*particle_set(index_a)%r(2)
             pv(3,2) = pv(3,2) + fc2(3)*particle_set(index_b)%r(2)
             pv(3,2) = pv(3,2) + fc3(3)*particle_set(index_c)%r(2)
             pv(3,2) = pv(3,2) + fc4(3)*particle_set(index_d)%r(2)
             !pv(3,3)
             pv(3,3) = pv(3,3) + fc1(3)*particle_set(index_a)%r(3)
             pv(3,3) = pv(3,3) + fc2(3)*particle_set(index_b)%r(3)
             pv(3,3) = pv(3,3) + fc3(3)*particle_set(index_c)%r(3)
             pv(3,3) = pv(3,3) + fc4(3)*particle_set(index_d)%r(3)

          END DO
       END DO
    END DO MOL
    CALL mp_sum ( pv, group )
    virial % pv_constraint = pv

  END SUBROUTINE pv_constraint

!!****** constraint/check_tol [1.0] *
!!
!!   NAME
!!     check_tol
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
  SUBROUTINE check_tol(roll_tol,iroll,char,matrix,veps)

    REAL(KIND=dp), INTENT(OUT)               :: roll_tol
    INTEGER, INTENT(INOUT)                   :: iroll
    CHARACTER(LEN=*), INTENT(IN)             :: char
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN), OPTIONAL                   :: matrix, veps

    CHARACTER(LEN=*), PARAMETER :: routineN = 'check_tol', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp)                            :: local_tol
    REAL(KIND=dp), DIMENSION(3, 3)           :: diff_rattle, diff_shake
    REAL(KIND=dp), DIMENSION(3, 3), SAVE     :: matrix_old, veps_old

    SELECT CASE ( char )
    CASE ( 'SHAKE' )
       IF (iroll==1) THEN
          matrix_old = matrix
          roll_tol = -1.E10_dp
       ELSE
          roll_tol = 0.0_dp
          diff_shake = ABS ( matrix_old - matrix )
          local_tol = MAXVAL ( diff_shake )
          roll_tol = MAX ( roll_tol, local_tol )
          matrix_old = matrix
       END IF
       iroll = iroll + 1
    CASE ( 'RATTLE' )
       IF (iroll==1) THEN
          veps_old = veps
          roll_tol = -1.E+10_dp
       ELSE
          roll_tol = 0.0_dp
          ! compute tolerance on veps
          diff_rattle = ABS ( veps - veps_old )
          local_tol = MAXVAL ( diff_rattle )
          roll_tol = MAX ( roll_tol, local_tol )
          veps_old = veps
       END IF
       iroll = iroll + 1
    END SELECT

  END SUBROUTINE check_tol

!!****** constraint/get_roll_matrix [1.0] *
!!
!!   NAME
!!     get_roll_matrix
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
  SUBROUTINE get_roll_matrix(char,r_shake,v_shake,vector_r,vector_v,u)


    CHARACTER(len=*), INTENT(IN)             :: char
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT), OPTIONAL                  :: r_shake, v_shake
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(IN), OPTIONAL                   :: vector_r, vector_v
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN), OPTIONAL                   :: u

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_roll_matrix', &
      routineP = moduleN//':'//routineN
    INTEGER                                  :: i
    REAL(KIND=dp), DIMENSION(3, 3)           :: diag

    IF ( PRESENT ( r_shake ) ) r_shake = 0.0_dp
    IF ( PRESENT ( v_shake ) ) v_shake = 0.0_dp
    diag = 0.0_dp

    SELECT CASE ( char )
    CASE ( 'SHAKE' )
       IF ( PRESENT ( u ) .AND. PRESENT ( vector_v ) .AND.  &
            PRESENT ( vector_r ) ) THEN
          diag ( 1, 1 ) = vector_r ( 1 )
          diag ( 2, 2 ) = vector_r ( 2 )
          diag ( 3, 3 ) = vector_r ( 3 )
          r_shake = MATMUL_3X3 ( MATMUL_3X3 ( u, diag ),TRANSPOSE_3D ( u ) )
          diag ( 1, 1 ) = vector_v ( 1 )
          diag ( 2, 2 ) = vector_v ( 2 )
          diag ( 3, 3 ) = vector_v ( 3 )
          v_shake = MATMUL_3X3 ( MATMUL_3X3 ( u, diag ), TRANSPOSE_3D ( u ) )
          diag = MATMUL_3x3 ( r_shake, v_shake )
          r_shake = diag
       ELSEIF ( .NOT. PRESENT ( u ) .AND. PRESENT ( vector_v ) .AND.  &
            PRESENT ( vector_r ) ) THEN
          DO i = 1, 3
             r_shake ( i, i ) = vector_r ( i ) * vector_v ( i )
             v_shake ( i, i ) = vector_v ( i )
          END DO
       ELSE
          CALL stop_program ( "get_roll_matrix", "not sufficient parameters" )
       END IF
    CASE ( 'RATTLE' )
       IF ( PRESENT ( u ) .AND. PRESENT ( vector_v ) ) THEN
          diag ( 1, 1 ) = vector_v ( 1 )
          diag ( 2, 2 ) = vector_v ( 2 )
          diag ( 3, 3 ) = vector_v ( 3 )
          v_shake = MATMUL_3x3 ( MATMUL_3X3 ( u, diag ), TRANSPOSE_3D ( u ) )
       ELSEIF ( .NOT. PRESENT ( u ) .AND. PRESENT ( vector_v ) ) THEN
          DO i = 1, 3
             v_shake ( i, i ) = vector_v ( i )
          END DO
       ELSE
          CALL stop_program ( "get_roll_matrix", "not sufficient parameters" )
       END IF
    END SELECT

  END SUBROUTINE get_roll_matrix

END MODULE constraint
