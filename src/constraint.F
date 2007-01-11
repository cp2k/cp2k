!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
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
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_types,                      ONLY: cell_type
  USE colvar_types,                    ONLY: colvar_counters
  USE constraint_3x3,                  ONLY: rattle_3x3_int,&
                                             rattle_roll_3x3_int,&
                                             shake_3x3_int,&
                                             shake_roll_3x3_int
  USE constraint_4x6,                  ONLY: rattle_4x6_int,&
                                             rattle_roll_4x6_int,&
                                             shake_4x6_int,&
                                             shake_roll_4x6_int
  USE constraint_clv,                  ONLY: rattle_colv_int,&
                                             rattle_roll_colv_int,&
                                             shake_colv_int,&
                                             shake_roll_colv_int
  USE constraint_util,                 ONLY: check_tol,&
                                             get_roll_matrix
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE input_constants,                 ONLY: npt_f_ensemble,&
                                             npt_i_ensemble
  USE kinds,                           ONLY: dp
  USE md,                              ONLY: simulation_parameters_type
  USE molecule_kind_types,             ONLY: get_molecule_kind,&
                                             molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE particle_types,                  ONLY: particle_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_types,                    ONLY: virial_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: shake_control,&
            rattle_control,&
            shake_roll_control,&
            rattle_roll_control

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
       particle_set, pos, vel, dt, shake_tol, log_unit, cell, error )

    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: pos( :, : ), vel( :, : )
    REAL(kind=dp), INTENT(in)                :: dt, shake_tol
    INTEGER, INTENT(in)                      :: log_unit
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'shake_control', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ikind, imol, &
                                                ishake, k, n3x3con, n4x6con, &
                                                nkind, nmol_per_kind
    REAL(KIND=dp)                            :: mass, max_sigma
    REAL(KIND=dp), DIMENSION(SIZE(pos, 2))   :: imass
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(colvar_counters)                    :: ncolv
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule

    CALL timeset ( routineN, 'I', ' ', handle )
    nkind = SIZE ( molecule_kind_set )
    DO k = 1, SIZE(pos,2)
       atomic_kind=>particle_set (k) % atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
       imass(k) = 1.0_dp/mass
    END DO
    MOL:  DO ikind = 1, nkind
       nmol_per_kind = local_molecules % n_el ( ikind )
       DO imol = 1, nmol_per_kind
          i = local_molecules % list ( ikind ) % array ( imol )
          molecule => molecule_set ( i ) 
          molecule_kind => molecule % molecule_kind
          CALL get_molecule_kind ( molecule_kind, ncolv = ncolv,  &
               ng3x3 = n3x3con, ng4x6=n4x6con )
          ishake = 0
          max_sigma = -1.0E+10_dp
          Shake_Intra_Loop: DO WHILE (ABS(max_sigma)>=shake_tol)
             max_sigma = 0.0_dp
             ishake  = ishake + 1
             ! 3x3
             IF ( n3x3con /= 0 ) &
                  CALL shake_3x3_int( molecule, particle_set, pos, vel, dt, ishake,&
                  log_unit, max_sigma, error)
             ! 4x6
             IF ( n4x6con /= 0 ) &
                  CALL shake_4x6_int( molecule, particle_set, pos, vel, dt, ishake,&
                  log_unit, max_sigma, error)
             ! Collective Variables
             IF ( ncolv%ntot /= 0 ) &
                  CALL shake_colv_int(molecule, particle_set, pos, vel, dt, ishake,&
                  log_unit, cell, imass, max_sigma, error )
          END DO Shake_Intra_Loop
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
       particle_set, vel, dt, rattle_tol, log_unit, cell, error )

    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: vel( :, : )
    REAL(kind=dp), INTENT(in)                :: dt, rattle_tol
    INTEGER, INTENT(in)                      :: log_unit
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rattle_control', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ikind, imol, &
                                                irattle, k, n3x3con, n4x6con, &
                                                nkind, nmol_per_kind
    REAL(KIND=dp)                            :: mass, max_sigma
    REAL(KIND=dp), DIMENSION(SIZE(vel, 2))   :: imass
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(colvar_counters)                    :: ncolv
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule

    CALL timeset ( routineN, 'I', ' ', handle )
    nkind = SIZE ( molecule_kind_set )
    DO k = 1, SIZE(vel,2)
       atomic_kind=>particle_set (k) % atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
       imass(k) = 1.0_dp/mass
    END DO
    MOL: DO ikind = 1, nkind
       nmol_per_kind = local_molecules % n_el ( ikind )
       DO imol = 1, nmol_per_kind
          i = local_molecules % list ( ikind ) % array ( imol )
          molecule => molecule_set ( i ) 
          molecule_kind => molecule % molecule_kind
          CALL get_molecule_kind ( molecule_kind, ncolv = ncolv, ng3x3 = n3x3con,&
               ng4x6 = n4x6con )
          irattle   = 0
          max_sigma = -1.0E+10_dp
          Rattle_Intra_Loop: DO WHILE ( ABS ( max_sigma ) >=  rattle_tol )
             max_sigma = 0.0_dp
             irattle = irattle + 1         
             ! 3x3
             IF ( n3x3con /= 0 ) &
                  CALL rattle_3x3_int( molecule, particle_set, vel, dt, error)
             ! 4x6
             IF ( n4x6con /= 0 ) &
                  CALL rattle_4x6_int( molecule, particle_set, vel, dt, error)
             ! Collective Variables
             IF ( ncolv%ntot /= 0 ) &
                  CALL rattle_colv_int( molecule, particle_set, vel, dt, irattle, log_unit,&
                  cell, imass, max_sigma, error )
          END DO Rattle_Intra_Loop
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
       vector_r, vector_v, group, u, cell, error )

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
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'shake_roll_control', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ikind, imol, &
                                                ishake, k, log_unit, n3x3con, &
                                                n4x6con, nkind, nmol_per_kind
    REAL(KIND=dp)                            :: dt, mass, max_sigma, shake_tol
    REAL(KIND=dp), DIMENSION(3, 3)           :: r_shake, v_shake
    REAL(KIND=dp), DIMENSION(SIZE(pos, 2))   :: imass
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(colvar_counters)                    :: ncolv
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule

    CALL timeset ( routineN, 'I', ' ', handle )
    nkind = SIZE ( molecule_kind_set )
    shake_tol = simpar%shake_tol
    dt = simpar%dt
    log_unit = simpar%info_constraint
    ! setting up for roll
    IF ( simpar % ensemble == npt_i_ensemble ) THEN
       CALL get_roll_matrix ('SHAKE', r_shake, v_shake, vector_r, vector_v )
    ELSE IF ( simpar % ensemble == npt_f_ensemble ) THEN
       CALL get_roll_matrix ('SHAKE', r_shake, v_shake, vector_r, vector_v, u )
    END IF
    DO k = 1, SIZE(pos,2)
       atomic_kind=>particle_set (k) % atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
       imass(k) = 1.0_dp/mass
    END DO
    MOL:  DO ikind = 1, nkind
       nmol_per_kind = local_molecules % n_el ( ikind )
       DO imol = 1, nmol_per_kind
          i = local_molecules % list ( ikind ) % array ( imol )
          molecule => molecule_set ( i ) 
          molecule_kind => molecule % molecule_kind
          CALL get_molecule_kind ( molecule_kind, ncolv = ncolv,  &
               ng3x3 = n3x3con, ng4x6=n4x6con )
          ishake = 0
          max_sigma = -1.0E+10_dp
          Shake_Roll_Intra_Loop: DO WHILE (ABS(max_sigma)>=shake_tol)
             max_sigma = 0.0_dp
             ishake  = ishake + 1
             ! 3x3
             IF ( n3x3con /= 0 ) &
                  CALL shake_roll_3x3_int( molecule, particle_set, pos, vel, r_shake,&
                  v_shake, dt, ishake, log_unit, max_sigma, error)
             ! 4x6
             IF ( n4x6con /= 0 ) &
                  CALL shake_roll_4x6_int( molecule, particle_set, pos, vel, r_shake,&
                  v_shake, dt, ishake, log_unit, max_sigma, error)
             ! Collective Variables
             IF ( ncolv%ntot /= 0 ) &
                  CALL shake_roll_colv_int( molecule, particle_set, pos, vel, r_shake,&
                  v_shake, dt, ishake, log_unit, cell, imass, max_sigma, error )
          END DO Shake_Roll_Intra_Loop
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
       veps, mass, roll_tol, iroll, para_env, u, cell, error )

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
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rattle_roll_control', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ikind, imol, &
                                                irattle, k, log_unit, &
                                                n3x3con, n4x6con, nkind, &
                                                nmol_per_kind
    REAL(KIND=dp)                            :: amass, dt, max_sigma, &
                                                rattle_tol
    REAL(KIND=dp), DIMENSION(3, 3)           :: r_rattle
    REAL(KIND=dp), DIMENSION(SIZE(vel, 2))   :: imass
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(colvar_counters)                    :: ncolv
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule

    CALL timeset ( routineN, 'I', ' ', handle )
    ! initialize locals
    nkind = SIZE ( molecule_kind_set )
    rattle_tol = simpar%shake_tol
    dt         = simpar%dt
    log_unit   = simpar%info_constraint
    ! setting up for roll
    IF ( simpar % ensemble == npt_i_ensemble ) THEN
       CALL get_roll_matrix ('RATTLE', v_shake=r_rattle, vector_v=vector )
    ELSE IF ( simpar % ensemble == npt_f_ensemble ) THEN
       CALL get_roll_matrix ('RATTLE', v_shake=r_rattle, vector_v=vector, u=u )
    END IF
    DO k = 1, SIZE(vel,2)
       atomic_kind=>particle_set (k) % atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=amass)
       imass(k) = 1.0_dp/amass
    END DO
    MOL: DO ikind = 1, nkind
       nmol_per_kind = local_molecules % n_el ( ikind )
       DO imol = 1, nmol_per_kind
          i = local_molecules % list ( ikind ) % array ( imol )
          molecule => molecule_set ( i ) 
          molecule_kind => molecule % molecule_kind
          CALL get_molecule_kind ( molecule_kind, ncolv = ncolv,  &
               ng3x3 = n3x3con, ng4x6 = n4x6con )
          max_sigma = -1.0E+10_dp
          irattle = 0
          Rattle_Roll_Intramolecular: DO WHILE (ABS(max_sigma)>=rattle_tol)
             max_sigma = 0.0_dp
             irattle = irattle + 1
             ! 3x3
             IF ( n3x3con /= 0 ) &
                  CALL rattle_roll_3x3_int( molecule, particle_set, vel, r_rattle, dt,&
                  veps, error )
             ! 4x6
             IF ( n4x6con /= 0 ) &
                  CALL rattle_roll_4x6_int( molecule, particle_set, vel, r_rattle, dt,&
                  veps, error )
             ! Collective Variables
             IF ( ncolv%ntot /= 0 ) &
                  CALL rattle_roll_colv_int( molecule, particle_set, vel, r_rattle, dt,&
                  irattle, log_unit, veps, cell, imass, max_sigma, error  )
          END DO Rattle_Roll_Intramolecular
       END DO
    END DO MOL
    CALL check_tol ( roll_tol, iroll, 'RATTLE', veps = veps )
    CALL timestop ( 0.0_dp, handle )
  END SUBROUTINE rattle_roll_control

END MODULE constraint
