!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      Teodoro Laino [tlaino] 2007 - Extension to Intermolecular constraints
! *****************************************************************************
MODULE constraint
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_types,                      ONLY: cell_type
  USE colvar_types,                    ONLY: colvar_counters
  USE constraint_3x3,                  ONLY: rattle_3x3_ext,&
                                             rattle_3x3_int,&
                                             rattle_roll_3x3_ext,&
                                             rattle_roll_3x3_int,&
                                             shake_3x3_ext,&
                                             shake_3x3_int,&
                                             shake_roll_3x3_ext,&
                                             shake_roll_3x3_int
  USE constraint_4x6,                  ONLY: rattle_4x6_ext,&
                                             rattle_4x6_int,&
                                             rattle_roll_4x6_ext,&
                                             rattle_roll_4x6_int,&
                                             shake_4x6_ext,&
                                             shake_4x6_int,&
                                             shake_roll_4x6_ext,&
                                             shake_roll_4x6_int
  USE constraint_clv,                  ONLY: &
       rattle_colv_ext, rattle_colv_int, rattle_roll_colv_ext, &
       rattle_roll_colv_int, shake_colv_ext, shake_colv_int, &
       shake_roll_colv_ext, shake_roll_colv_int, shake_update_colv_ext, &
       shake_update_colv_int
  USE constraint_util,                 ONLY: check_tol,&
                                             get_roll_matrix,&
                                             restore_temporary_set,&
                                             update_temporary_set
  USE constraint_vsite,                ONLY: shake_vsite_ext,&
                                             shake_vsite_int
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE input_constants,                 ONLY: npt_f_ensemble,&
                                             npt_i_ensemble
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_sum
  USE molecule_kind_types,             ONLY: get_molecule_kind,&
                                             get_molecule_kind_set,&
                                             molecule_kind_type
  USE molecule_types_new,              ONLY: global_constraint_type,&
                                             molecule_type
  USE particle_types,                  ONLY: particle_type
  USE simpar_types,                    ONLY: simpar_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: shake_control,&
            rattle_control,&
            shake_roll_control,&
            rattle_roll_control,&
            shake_update_targets

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'constraint'
  INTEGER, PARAMETER, PRIVATE :: Max_Shake_Iter = 1000

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param gci ...
!> \param local_molecules ...
!> \param molecule_set ...
!> \param molecule_kind_set ...
!> \param particle_set ...
!> \param pos ...
!> \param vel ...
!> \param dt ...
!> \param shake_tol ...
!> \param log_unit ...
!> \param lagrange_mult ...
!> \param dump_lm ...
!> \param cell ...
!> \param group ...
!> \param local_particles ...
!> \par History
!>      Teodoro Laino [tlaino] 2007 - Extension to Intermolecular constraints
! *****************************************************************************
  SUBROUTINE shake_control( gci, local_molecules, molecule_set, molecule_kind_set, &
       particle_set, pos, vel, dt, shake_tol, log_unit, lagrange_mult, dump_lm, &
       cell, group, local_particles)

    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: pos( :, : ), vel( :, : )
    REAL(kind=dp), INTENT(in)                :: dt, shake_tol
    INTEGER, INTENT(in)                      :: log_unit, lagrange_mult
    LOGICAL, INTENT(IN)                      :: dump_lm
    TYPE(cell_type), POINTER                 :: cell
    INTEGER, INTENT(in)                      :: group
    TYPE(distribution_1d_type), POINTER      :: local_particles

    CHARACTER(LEN=*), PARAMETER :: routineN = 'shake_control', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, ikind, imol, ishake_ext, ishake_int, k, n3x3con, &
      n4x6con, nconstraint, nkind, nmol_per_kind, nvsitecon
    LOGICAL                                  :: do_ext_constraint
    REAL(KIND=dp)                            :: int_max_sigma, mass, max_sigma
    REAL(KIND=dp), DIMENSION(SIZE(pos, 2))   :: imass
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(colvar_counters)                    :: ncolv
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule

    CALL timeset(routineN,handle)
    nkind = SIZE ( molecule_kind_set )
    DO k = 1, SIZE(pos,2)
       atomic_kind=>particle_set (k) % atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
       imass(k) = 1.0_dp/MAX(mass,EPSILON(0.0_dp))
    END DO
    do_ext_constraint = (gci%ntot/=0)
    ishake_ext = 0
    max_sigma  = -1.0E+10_dp
    Shake_Inter_Loop: DO WHILE ((ABS(max_sigma)>=shake_tol).AND.(ishake_ext<=Max_Shake_Iter))
       max_sigma = 0.0_dp
       ishake_ext  = ishake_ext + 1
       ! Intramolecular Constraints
       MOL:  DO ikind = 1, nkind
          nmol_per_kind = local_molecules % n_el ( ikind )
          DO imol = 1, nmol_per_kind
             i = local_molecules % list ( ikind ) % array ( imol )
             molecule => molecule_set ( i )
             molecule_kind => molecule % molecule_kind
             CALL get_molecule_kind ( molecule_kind, ncolv = ncolv,  &
                                      ng3x3 = n3x3con, ng4x6=n4x6con,&
                                      nconstraint=nconstraint,nvsite=nvsitecon)
             IF (nconstraint==0) CYCLE
             ishake_int = 0
             int_max_sigma  = -1.0E+10_dp
             Shake_Intra_Loop: DO WHILE ((ABS(int_max_sigma)>=shake_tol).AND.(ishake_int<=Max_Shake_Iter))
                int_max_sigma   = 0.0_dp
                ishake_int  = ishake_int + 1
                ! 3x3
                IF ( n3x3con /= 0 ) &
                     CALL shake_3x3_int( molecule, particle_set, pos, vel, dt, ishake_int,&
                     int_max_sigma)
                ! 4x6
                IF ( n4x6con /= 0 ) &
                     CALL shake_4x6_int( molecule, particle_set, pos, vel, dt, ishake_int,&
                     int_max_sigma)
                ! Collective Variables
                IF ( ncolv%ntot /= 0 ) &
                     CALL shake_colv_int(molecule, particle_set, pos, vel, dt, ishake_int,&
                     cell, imass, int_max_sigma)
             END DO Shake_Intra_Loop
             max_sigma = MAX(max_sigma, int_max_sigma)
             CALL shake_int_info(log_unit,i,ishake_int,max_sigma)
             ! Virtual Site
             IF ( nvsitecon /= 0 ) &
                CALL shake_vsite_int( molecule, pos)
          END DO
       END DO MOL
       ! Intermolecular constraints
       IF (do_ext_constraint) THEN
          CALL update_temporary_set ( group, pos = pos, vel = vel )
          ! 3x3
          IF ( gci%ng3x3 /= 0 ) &
               CALL shake_3x3_ext( gci, particle_set, pos, vel, dt, ishake_ext,&
               max_sigma)
          ! 4x6
          IF ( gci%ng4x6 /= 0 ) &
               CALL shake_4x6_ext( gci, particle_set, pos, vel, dt, ishake_ext,&
               max_sigma)
          ! Collective Variables
          IF ( gci%ncolv%ntot /= 0 ) &
               CALL shake_colv_ext(gci, particle_set, pos, vel, dt, ishake_ext,&
               cell, imass, max_sigma)
          ! Virtual Site
          IF ( gci%nvsite /= 0 ) &
               CALL shake_vsite_ext( gci, pos)
          CALL restore_temporary_set( particle_set, local_particles, pos=pos, vel=vel )
       END IF
       CALL shake_ext_info(log_unit,ishake_ext,max_sigma)
    END DO Shake_Inter_Loop
    CALL dump_lagrange_mult(dump_lm, lagrange_mult,local_molecules,molecule_set,gci,&
         molecule_kind_set,group,"S")

    CALL timestop(handle)
  END SUBROUTINE shake_control

! *****************************************************************************
!> \brief ...
!> \param gci ...
!> \param local_molecules ...
!> \param molecule_set ...
!> \param molecule_kind_set ...
!> \param particle_set ...
!> \param vel ...
!> \param dt ...
!> \param rattle_tol ...
!> \param log_unit ...
!> \param lagrange_mult ...
!> \param dump_lm ...
!> \param cell ...
!> \param group ...
!> \param local_particles ...
!> \par History
!>      Teodoro Laino [tlaino] 2007 - Extension to Intermolecular constraints
! *****************************************************************************
  SUBROUTINE rattle_control( gci, local_molecules, molecule_set, molecule_kind_set,  &
       particle_set, vel, dt, rattle_tol, log_unit, lagrange_mult, dump_lm, cell, group,&
       local_particles)

    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: vel( :, : )
    REAL(kind=dp), INTENT(in)                :: dt, rattle_tol
    INTEGER, INTENT(in)                      :: log_unit, lagrange_mult
    LOGICAL, INTENT(IN)                      :: dump_lm
    TYPE(cell_type), POINTER                 :: cell
    INTEGER, INTENT(in)                      :: group
    TYPE(distribution_1d_type), POINTER      :: local_particles

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rattle_control', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, ikind, imol, irattle_ext, irattle_int, k, n3x3con, &
      n4x6con, nconstraint, nkind, nmol_per_kind
    LOGICAL                                  :: do_ext_constraint
    REAL(KIND=dp)                            :: int_max_sigma, mass, max_sigma
    REAL(KIND=dp), DIMENSION(SIZE(vel, 2))   :: imass
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(colvar_counters)                    :: ncolv
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule

    CALL timeset(routineN,handle)
    nkind = SIZE ( molecule_kind_set )
    DO k = 1, SIZE(vel,2)
       atomic_kind=>particle_set (k) % atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
       imass(k) = 1.0_dp/MAX(mass,EPSILON(0.0_dp))
    END DO
    do_ext_constraint = (gci%ntot/=0)
    irattle_ext = 0
    max_sigma  = -1.0E+10_dp
    Rattle_Inter_Loop: DO WHILE (ABS(max_sigma)>=rattle_tol)
       max_sigma = 0.0_dp
       irattle_ext  = irattle_ext + 1
       ! Intramolecular Constraints
       MOL: DO ikind = 1, nkind
          nmol_per_kind = local_molecules % n_el ( ikind )
          DO imol = 1, nmol_per_kind
             i = local_molecules % list ( ikind ) % array ( imol )
             molecule => molecule_set ( i )
             molecule_kind => molecule % molecule_kind
             CALL get_molecule_kind ( molecule_kind, ncolv = ncolv, ng3x3 = n3x3con,&
                                      ng4x6 = n4x6con ,nconstraint=nconstraint )
             IF (nconstraint==0) CYCLE
             irattle_int   = 0
             int_max_sigma = -1.0E+10_dp
             Rattle_Intra_Loop: DO WHILE ( ABS ( int_max_sigma ) >=  rattle_tol )
                int_max_sigma = 0.0_dp
                irattle_int = irattle_int + 1
                ! 3x3
                IF ( n3x3con /= 0 ) &
                     CALL rattle_3x3_int( molecule, particle_set, vel, dt)
                ! 4x6
                IF ( n4x6con /= 0 ) &
                     CALL rattle_4x6_int( molecule, particle_set, vel, dt)
                ! Collective Variables
                IF ( ncolv%ntot /= 0 ) &
                     CALL rattle_colv_int( molecule, particle_set, vel, dt, &
                     irattle_int, cell, imass, int_max_sigma)
             END DO Rattle_Intra_Loop
             max_sigma = MAX(max_sigma, int_max_sigma)
             CALL rattle_int_info(log_unit,i,irattle_int,max_sigma)
          END DO
       END DO MOL
       ! Intermolecular Constraints
       IF (do_ext_constraint) THEN
          CALL update_temporary_set ( group, vel = vel )
          ! 3x3
          IF ( gci%ng3x3 /= 0 ) &
               CALL rattle_3x3_ext( gci, particle_set, vel, dt)
          ! 4x6
          IF ( gci%ng4x6 /= 0 ) &
               CALL rattle_4x6_ext( gci, particle_set, vel, dt)
          ! Collective Variables
          IF ( gci%ncolv%ntot /= 0 ) &
               CALL rattle_colv_ext( gci, particle_set, vel, dt, &
               irattle_ext, cell, imass, max_sigma)
          CALL restore_temporary_set ( particle_set, local_particles, vel = vel )
       END IF
       CALL rattle_ext_info(log_unit,irattle_ext,max_sigma)
    END DO Rattle_Inter_Loop
    CALL dump_lagrange_mult(dump_lm,lagrange_mult,local_molecules,molecule_set,gci,&
         molecule_kind_set,group,"R")
    CALL timestop(handle)

  END SUBROUTINE rattle_control

! *****************************************************************************
!> \brief ...
!> \param gci ...
!> \param local_molecules ...
!> \param molecule_set ...
!> \param molecule_kind_set ...
!> \param particle_set ...
!> \param pos ...
!> \param vel ...
!> \param dt ...
!> \param simpar ...
!> \param roll_tol ...
!> \param iroll ...
!> \param vector_r ...
!> \param vector_v ...
!> \param group ...
!> \param u ...
!> \param cell ...
!> \param local_particles ...
!> \par History
!>      Teodoro Laino [tlaino] 2007 - Extension to Intermolecular constraints
! *****************************************************************************
  SUBROUTINE shake_roll_control( gci, local_molecules, molecule_set,&
       molecule_kind_set, particle_set, pos, vel, dt, simpar, roll_tol, iroll,&
       vector_r, vector_v, group, u, cell, local_particles)

    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: pos( :, : ), vel( :, : )
    REAL(KIND=dp), INTENT(IN)                :: dt
    TYPE(simpar_type), INTENT(IN)            :: simpar
    REAL(KIND=dp), INTENT(OUT)               :: roll_tol
    INTEGER, INTENT(INOUT)                   :: iroll
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: vector_r, vector_v
    INTEGER, INTENT(IN)                      :: group
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN), OPTIONAL                   :: u
    TYPE(cell_type), POINTER                 :: cell
    TYPE(distribution_1d_type), POINTER      :: local_particles

    CHARACTER(LEN=*), PARAMETER :: routineN = 'shake_roll_control', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, ikind, imol, ishake_ext, ishake_int, k, &
      lagrange_mult, log_unit, n3x3con, n4x6con, nconstraint, nkind, &
      nmol_per_kind, nvsitecon
    LOGICAL                                  :: do_ext_constraint, dump_lm
    REAL(KIND=dp)                            :: int_max_sigma, mass, &
                                                max_sigma, shake_tol
    REAL(KIND=dp), DIMENSION(3, 3)           :: r_shake, v_shake
    REAL(KIND=dp), DIMENSION(SIZE(pos, 2))   :: imass
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(colvar_counters)                    :: ncolv
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule

    CALL timeset(routineN,handle)
    nkind = SIZE ( molecule_kind_set )
    shake_tol = simpar%shake_tol
    log_unit = simpar%info_constraint
    lagrange_mult = simpar%lagrange_multipliers
    dump_lm = simpar%dump_lm
    ! setting up for roll
    IF ( simpar % ensemble == npt_i_ensemble ) THEN
       CALL get_roll_matrix ('SHAKE', r_shake, v_shake, vector_r, vector_v )
    ELSE IF ( simpar % ensemble == npt_f_ensemble ) THEN
       CALL get_roll_matrix ('SHAKE', r_shake, v_shake, vector_r, vector_v, u )
    END IF
    DO k = 1, SIZE(pos,2)
       atomic_kind=>particle_set (k) % atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
       imass(k) = 1.0_dp/MAX(mass,EPSILON(0.0_dp))
    END DO
    do_ext_constraint = (gci%ntot/=0)
    ishake_ext = 0
    max_sigma  = -1.0E+10_dp
    Shake_Inter_Loop: DO WHILE (ABS(max_sigma)>=shake_tol)
       max_sigma = 0.0_dp
       ishake_ext  = ishake_ext + 1
       ! Intramolecular Constraints
       MOL:  DO ikind = 1, nkind
          nmol_per_kind = local_molecules % n_el ( ikind )
          DO imol = 1, nmol_per_kind
             i = local_molecules % list ( ikind ) % array ( imol )
             molecule => molecule_set ( i )
             molecule_kind => molecule % molecule_kind
             CALL get_molecule_kind ( molecule_kind, ncolv = ncolv,  &
                                      ng3x3 = n3x3con, ng4x6=n4x6con,&
                                      nconstraint=nconstraint,nvsite=nvsitecon)
             IF (nconstraint==0) CYCLE
             ishake_int = 0
             int_max_sigma = -1.0E+10_dp
             Shake_Roll_Intra_Loop: DO WHILE (ABS(int_max_sigma)>=shake_tol)
                int_max_sigma = 0.0_dp
                ishake_int  = ishake_int + 1
                ! 3x3
                IF ( n3x3con /= 0 ) &
                     CALL shake_roll_3x3_int( molecule, particle_set, pos, vel, r_shake,&
                     v_shake, dt, ishake_int, int_max_sigma)
                ! 4x6
                IF ( n4x6con /= 0 ) &
                     CALL shake_roll_4x6_int( molecule, particle_set, pos, vel, r_shake,&
                     dt, ishake_int, int_max_sigma)
                ! Collective Variables
                IF ( ncolv%ntot /= 0 ) &
                     CALL shake_roll_colv_int( molecule, particle_set, pos, vel, r_shake,&
                     v_shake, dt, ishake_int, cell, imass, int_max_sigma)
             END DO Shake_Roll_Intra_Loop
             max_sigma = MAX(max_sigma, int_max_sigma)
             CALL shake_int_info(log_unit,i,ishake_int,max_sigma)
             ! Virtual Site
             IF ( nvsitecon /= 0 ) THEN
                CALL cp_unimplemented_error(fromWhere=routineP, &
                     message="Virtual Site Constraint/Restraint not implemented for SHAKE_ROLL!", &
                     error_level=cp_failure_level)
             END IF
          END DO
       END DO MOL
       ! Intermolecular constraints
       IF (do_ext_constraint) THEN
          CALL update_temporary_set ( group, pos = pos, vel = vel )
          ! 3x3
          IF ( gci%ng3x3 /= 0 ) &
               CALL shake_roll_3x3_ext( gci, particle_set, pos, vel, r_shake,&
               v_shake, dt, ishake_ext, max_sigma)
          ! 4x6
          IF ( gci%ng4x6 /= 0 ) &
               CALL shake_roll_4x6_ext( gci, particle_set, pos, vel, r_shake,&
               dt, ishake_ext, max_sigma)
          ! Collective Variables
          IF ( gci%ncolv%ntot /= 0 ) &
               CALL shake_roll_colv_ext( gci, particle_set, pos, vel, r_shake,&
               v_shake, dt, ishake_ext, cell, imass, max_sigma)
          ! Virtual Site
          IF ( gci%nvsite /= 0 ) &
               CALL cp_unimplemented_error(fromWhere=routineP, &
               message="Virtual Site Constraint/Restraint not implemented for SHAKE_ROLL!", &
               error_level=cp_failure_level)
          CALL restore_temporary_set ( particle_set, local_particles, pos = pos, vel = vel )
       END IF
       CALL shake_ext_info(log_unit,ishake_ext,max_sigma)
    END DO Shake_Inter_Loop
    CALL dump_lagrange_mult(dump_lm,lagrange_mult,local_molecules,molecule_set,gci,&
         molecule_kind_set,group,"S")
    CALL check_tol ( roll_tol, iroll, 'SHAKE', r_shake )
    CALL timestop(handle)

  END SUBROUTINE shake_roll_control

! *****************************************************************************
!> \brief ...
!> \param gci ...
!> \param local_molecules ...
!> \param molecule_set ...
!> \param molecule_kind_set ...
!> \param particle_set ...
!> \param vel ...
!> \param dt ...
!> \param simpar ...
!> \param vector ...
!> \param veps ...
!> \param roll_tol ...
!> \param iroll ...
!> \param para_env ...
!> \param u ...
!> \param cell ...
!> \param local_particles ...
!> \par History
!>      Teodoro Laino [tlaino] 2007 - Extension to Intermolecular constraints
! *****************************************************************************
  SUBROUTINE rattle_roll_control( gci, local_molecules, molecule_set, &
       molecule_kind_set, particle_set, vel, dt, simpar, vector,&
       veps, roll_tol, iroll, para_env, u, cell, local_particles)

    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: vel( :, : )
    REAL(KIND=dp), INTENT(IN)                :: dt
    TYPE(simpar_type), INTENT(IN)            :: simpar
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: vector
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: veps
    REAL(KIND=dp), INTENT(OUT)               :: roll_tol
    INTEGER, INTENT(INOUT)                   :: iroll
    TYPE(cp_para_env_type), INTENT(IN)       :: para_env
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN), OPTIONAL                   :: u
    TYPE(cell_type), POINTER                 :: cell
    TYPE(distribution_1d_type), POINTER      :: local_particles

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rattle_roll_control', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, ikind, imol, irattle_ext, irattle_int, k, &
      lagrange_mult, log_unit, n3x3con, n4x6con, nconstraint, nkind, &
      nmol_per_kind
    LOGICAL                                  :: do_ext_constraint, dump_lm
    REAL(KIND=dp)                            :: int_max_sigma, mass, &
                                                max_sigma, rattle_tol
    REAL(KIND=dp), DIMENSION(3, 3)           :: r_rattle
    REAL(KIND=dp), DIMENSION(SIZE(vel, 2))   :: imass
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(colvar_counters)                    :: ncolv
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule

    CALL timeset(routineN,handle)
    ! initialize locals
    nkind = SIZE ( molecule_kind_set )
    rattle_tol = simpar%shake_tol
    log_unit   = simpar%info_constraint
    lagrange_mult = simpar%lagrange_multipliers
    dump_lm = simpar%dump_lm
    ! setting up for roll
    IF ( simpar % ensemble == npt_i_ensemble ) THEN
       CALL get_roll_matrix ('RATTLE', v_shake=r_rattle, vector_v=vector )
    ELSE IF ( simpar % ensemble == npt_f_ensemble ) THEN
       CALL get_roll_matrix ('RATTLE', v_shake=r_rattle, vector_v=vector, u=u )
    END IF
    DO k = 1, SIZE(vel,2)
       atomic_kind=>particle_set (k) % atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
       imass(k) = 1.0_dp/MAX(mass,EPSILON(0.0_dp))
    END DO
    do_ext_constraint = (gci%ntot/=0)
    irattle_ext = 0
    max_sigma  = -1.0E+10_dp
    Rattle_Inter_Loop: DO WHILE (ABS(max_sigma)>=rattle_tol)
       max_sigma = 0.0_dp
       irattle_ext  = irattle_ext + 1
       ! Intramolecular Constraints
       MOL: DO ikind = 1, nkind
          nmol_per_kind = local_molecules % n_el ( ikind )
          DO imol = 1, nmol_per_kind
             i = local_molecules % list ( ikind ) % array ( imol )
             molecule => molecule_set ( i )
             molecule_kind => molecule % molecule_kind
             CALL get_molecule_kind ( molecule_kind, ncolv = ncolv,  &
                                      ng3x3 = n3x3con, ng4x6 = n4x6con,&
                                      nconstraint=nconstraint )
             IF (nconstraint==0) CYCLE
             int_max_sigma = -1.0E+10_dp
             irattle_int = 0
             Rattle_Roll_Intramolecular: DO WHILE (ABS(int_max_sigma)>=rattle_tol)
                int_max_sigma = 0.0_dp
                irattle_int = irattle_int + 1
                ! 3x3
                IF ( n3x3con /= 0 ) &
                     CALL rattle_roll_3x3_int( molecule, particle_set, vel, r_rattle, dt,&
                     veps)
                ! 4x6
                IF ( n4x6con /= 0 ) &
                     CALL rattle_roll_4x6_int( molecule, particle_set, vel, r_rattle, dt,&
                     veps)
                ! Collective Variables
                IF ( ncolv%ntot /= 0 ) &
                     CALL rattle_roll_colv_int( molecule, particle_set, vel, r_rattle, dt,&
                     irattle_int, veps, cell, imass, int_max_sigma)
             END DO Rattle_Roll_Intramolecular
             max_sigma = MAX(max_sigma, int_max_sigma)
             CALL rattle_int_info(log_unit,i,irattle_int,max_sigma)
          END DO
       END DO MOL
       ! Intermolecular Constraints
       IF (do_ext_constraint) THEN
          CALL update_temporary_set ( para_env%group, vel = vel )
          ! 3x3
          IF ( gci%ng3x3 /= 0 )&
               CALL rattle_roll_3x3_ext( gci, particle_set, vel, r_rattle, dt,&
               veps)
          ! 4x6
          IF ( gci%ng4x6 /= 0 ) &
               CALL rattle_roll_4x6_ext( gci, particle_set, vel, r_rattle, dt,&
               veps)
          ! Collective Variables
          IF ( gci%ncolv%ntot /= 0 ) &
               CALL rattle_roll_colv_ext( gci, particle_set, vel, r_rattle, dt,&
               irattle_ext, veps, cell, imass, max_sigma)
          CALL restore_temporary_set ( particle_set, local_particles, vel = vel )
       END IF
       CALL rattle_ext_info(log_unit,irattle_ext,max_sigma)
    END DO Rattle_Inter_Loop
    CALL dump_lagrange_mult(dump_lm,lagrange_mult,local_molecules,molecule_set,gci,&
         molecule_kind_set,para_env%group,"R")
    CALL check_tol ( roll_tol, iroll, 'RATTLE', veps = veps )
    CALL timestop(handle)
  END SUBROUTINE rattle_roll_control

! *****************************************************************************
!> \brief ...
!> \param dump_lm ...
!> \param log_unit ...
!> \param local_molecules ...
!> \param molecule_set ...
!> \param gci ...
!> \param molecule_kind_set ...
!> \param group ...
!> \param id_type ...
!> \par History
!>      Teodoro Laino [tlaino] 2007 - Dumps lagrange multipliers
! *****************************************************************************
  SUBROUTINE dump_lagrange_mult(dump_lm,log_unit,local_molecules,molecule_set,gci,&
       molecule_kind_set, group, id_type)
    LOGICAL, INTENT(IN)                      :: dump_lm
    INTEGER, INTENT(IN)                      :: log_unit
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    INTEGER, INTENT(IN)                      :: group
    CHARACTER(LEN=1), INTENT(IN)             :: id_type

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dump_lagrange_mult', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: label
    INTEGER                                  :: handle, i, ikind, imol, j, &
                                                my_index, n3x3con, n4x6con, &
                                                nconstraint, nkind
    LOGICAL                                  :: do_ext_constraint, &
                                                do_int_constraint, failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: lagr
    TYPE(colvar_counters)                    :: ncolv
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule

    failure = .FALSE.
    CALL timeset(routineN,handle)
    ! Total number of intramolecular constraints (distributed)
    CALL get_molecule_kind_set(molecule_kind_set=molecule_kind_set,&
         nconstraint=nconstraint)
    do_int_constraint = (nconstraint>0)
    do_ext_constraint = (gci%ntot>0)
    IF (dump_lm.AND.(do_int_constraint.OR.do_ext_constraint)) THEN
       nkind = SIZE (  molecule_kind_set )
       ALLOCATE(lagr(nconstraint))
       lagr = 0.0_dp
       ! Dump lagrange multipliers for Intramolecular Constraints
       my_index = 0
       IF (do_int_constraint) THEN
          MOL:  DO ikind = 1, nkind
             molecule_kind => molecule_kind_set(ikind)
             CALL get_molecule_kind ( molecule_kind,&
                                      ncolv = ncolv,&
                                      ng3x3 = n3x3con,&
                                      ng4x6 = n4x6con )
             DO imol = 1, molecule_kind%nmolecule
                i = molecule_kind%molecule_list(imol)
                IF (ANY(local_molecules % list ( ikind ) % array == i)) THEN
                   molecule => molecule_set ( i )
                   ! Collective Variables
                   DO j = 1, ncolv%ntot
                      lagr(my_index+1) = molecule%lci%lcolv(j)%lambda
                      my_index = my_index + 1
                   END DO
                   ! 3x3
                   DO j = 1, n3x3con
                      lagr(my_index+1:my_index+3) = molecule%lci%lg3x3(j)%lambda(:)
                      my_index = my_index + 3
                   END DO
                   ! 4x6
                   DO j = 1, n4x6con
                      lagr(my_index+1:my_index+6) = molecule%lci%lg4x6(j)%lambda(:)
                      my_index = my_index + 6
                   END DO
                ELSE
                   my_index = my_index + ncolv%ntot + 3*n3x3con + 6*n4x6con
                END IF
             END DO
          END DO MOL
          CALL mp_sum(lagr,group)
       END IF
       ! Intermolecular constraints
       IF (do_ext_constraint) THEN
          CALL reallocate(lagr,1,SIZE(lagr)+gci%ntot)
          ! Collective Variables
          DO j = 1, gci%ncolv%ntot
             lagr(my_index+1) = gci%lcolv(j)%lambda
             my_index = my_index + 1
          END DO
          ! 3x3
          DO j = 1, gci%ng3x3
             lagr(my_index+1:my_index+3) = gci%lg3x3(j)%lambda(:)
             my_index = my_index + 3
          END DO
          ! 4x6
          DO j = 1, gci%ng4x6
             lagr(my_index+1:my_index+6) = gci%lg4x6(j)%lambda(:)
             my_index = my_index + 6
          END DO
       END IF
       IF (log_unit>0) THEN
          IF (id_type=="S") THEN
             label = "Shake  Lagrangian Multipliers:"
          ELSEIF (id_type=="R") THEN
             label = "Rattle Lagrangian Multipliers:"
          ELSE
             CPASSERT(.FALSE.)
          END IF
          WRITE(log_unit,FMT='(A,T40,4F15.9)')TRIM(label),lagr(1:MIN(4,SIZE(lagr)))
          DO j = 5, SIZE(lagr), 4
             WRITE(log_unit,FMT='(T40,4F15.9)')lagr(j:MIN(j+3,SIZE(lagr)))
          END DO
       END IF
       DEALLOCATE(lagr)
    END IF
    CALL timestop(handle)

  END SUBROUTINE dump_lagrange_mult

! *****************************************************************************
!> \brief Dumps convergence info about shake - intramolecular constraint loop
!> \param log_unit ...
!> \param i ...
!> \param ishake_int ...
!> \param max_sigma ...
!> \par History
!>      Teodoro Laino [tlaino] 2007 - University of Zurich
! *****************************************************************************
  SUBROUTINE shake_int_info(log_unit,i,ishake_int,max_sigma)
    INTEGER, INTENT(IN)                      :: log_unit, i, ishake_int
    REAL(KIND=dp), INTENT(IN)                :: max_sigma

    CHARACTER(LEN=*), PARAMETER :: routineN = 'shake_int_info', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (log_unit>0) THEN
       ! Dump info if requested
       WRITE(log_unit,'("SHAKE_INFO|",2X,2(A,I6),A,F15.9)')&
            "Molecule Nr.:",i," Nr. Iterations:",ishake_int," Max. Err.:",max_sigma
    END IF
    ! Notify a not converged SHAKE
    CALL cp_assert(ishake_int<=Max_Shake_Iter,cp_warning_level,cp_assertion_failed,routineP,&
         "Shake NOT converged in "//cp_to_string(Max_Shake_Iter)//" iterations in the "//&
         " intramolecular constraint loop for Molecule nr. "//cp_to_string(i)//&
         ". CP2K continues but results could be meaningless. "//&
CPSourceFileRef)
  END SUBROUTINE shake_int_info

! *****************************************************************************
!> \brief Dumps convergence info about shake - intermolecular constraint loop
!> \param log_unit ...
!> \param ishake_ext ...
!> \param max_sigma ...
!> \par History
!>      Teodoro Laino [tlaino] 2007 - University of Zurich
! *****************************************************************************
  SUBROUTINE shake_ext_info(log_unit,ishake_ext,max_sigma)
    INTEGER, INTENT(IN)                      :: log_unit, ishake_ext
    REAL(KIND=dp), INTENT(IN)                :: max_sigma

    CHARACTER(LEN=*), PARAMETER :: routineN = 'shake_ext_info', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (log_unit>0) THEN
       ! Dump info if requested
       WRITE(log_unit,'("SHAKE_INFO|",2X,A,I6,A,F15.9)')&
            "External Shake      Nr. Iterations:",ishake_ext,&
            " Max. Err.:",max_sigma
    END IF
    ! Notify a not converged SHAKE
    CALL cp_assert(ishake_ext<=Max_Shake_Iter,cp_warning_level,cp_assertion_failed,routineP,&
         "Shake NOT converged in "//cp_to_string(Max_Shake_Iter)//" iterations in the "//&
         " intermolecular constraint. CP2K continues but results could be meaningless. "//&
CPSourceFileRef)
  END SUBROUTINE shake_ext_info

! *****************************************************************************
!> \brief Dumps convergence info about rattle - intramolecular constraint loop
!> \param log_unit ...
!> \param i ...
!> \param irattle_int ...
!> \param max_sigma ...
!> \par History
!>      Teodoro Laino [tlaino] 2007 - University of Zurich
! *****************************************************************************
  SUBROUTINE rattle_int_info(log_unit,i,irattle_int,max_sigma)
    INTEGER, INTENT(IN)                      :: log_unit, i, irattle_int
    REAL(KIND=dp), INTENT(IN)                :: max_sigma

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rattle_int_info', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (log_unit>0) THEN
       ! Dump info if requested
       WRITE(log_unit,'("RATTLE_INFO|",1X,2(A,I6),A,F15.9)')&
            "Molecule Nr.:",i," Nr. Iterations:",irattle_int," Max. Err.:",max_sigma
    END IF
    ! Notify a not converged RATTLE
    CALL cp_assert(irattle_int<=Max_shake_Iter,cp_warning_level,cp_assertion_failed,routineP,&
         "Rattle NOT converged in "//cp_to_string(Max_Shake_Iter)//" iterations in the "//&
         " intramolecular constraint loop for Molecule nr. "//cp_to_string(i)//&
         ". CP2K continues but results could be meaningless. "//&
CPSourceFileRef)
  END SUBROUTINE rattle_int_info

! *****************************************************************************
!> \brief Dumps convergence info about rattle - intermolecular constraint loop
!> \param log_unit ...
!> \param irattle_ext ...
!> \param max_sigma ...
!> \par History
!>      Teodoro Laino [tlaino] 2007 - University of Zurich
! *****************************************************************************
  SUBROUTINE rattle_ext_info(log_unit,irattle_ext,max_sigma)
    INTEGER, INTENT(IN)                      :: log_unit, irattle_ext
    REAL(KIND=dp), INTENT(IN)                :: max_sigma

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rattle_ext_info', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (log_unit>0) THEN
       ! Dump info if requested
       WRITE(log_unit,'("RATTLE_INFO|",1X,A,I6,A,F15.9)')&
            "External Rattle     Nr. Iterations:",irattle_ext,&
            " Max. Err.:",max_sigma
    END IF
    ! Notify a not converged RATTLE
    CALL cp_assert(irattle_ext<=Max_shake_Iter,cp_warning_level,cp_assertion_failed,routineP,&
         "Rattle NOT converged in "//cp_to_string(Max_Shake_Iter)//" iterations in the "//&
         " intermolecular constraint. CP2K continues but results could be meaningless. "//&
CPSourceFileRef)
  END SUBROUTINE rattle_ext_info

! *****************************************************************************
!> \brief Updates the TARGET of the COLLECTIVE constraints if the growth speed
!>        is different from zero.
!> \param gci ...
!> \param local_molecules ...
!> \param molecule_set ...
!> \param molecule_kind_set ...
!> \param dt ...
!> \param root_section ...
!> \date    02.2008
!> \author  Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE shake_update_targets( gci, local_molecules, molecule_set, &
       molecule_kind_set, dt, root_section)

    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    REAL(kind=dp), INTENT(in)                :: dt
    TYPE(section_vals_type), POINTER         :: root_section

    CHARACTER(LEN=*), PARAMETER :: routineN = 'shake_update_targets', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ikind, imol, &
                                                nkind, nmol_per_kind
    LOGICAL                                  :: do_ext_constraint
    TYPE(colvar_counters)                    :: ncolv
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule
    TYPE(section_vals_type), POINTER         :: motion_section

    CALL timeset(routineN,handle)
    motion_section => section_vals_get_subs_vals(root_section,"MOTION")
    nkind = SIZE ( molecule_kind_set )
    do_ext_constraint = (gci%ntot/=0)
    ! Intramolecular Constraints
    MOL:  DO ikind = 1, nkind
       nmol_per_kind = local_molecules % n_el ( ikind )
       DO imol = 1, nmol_per_kind
          i = local_molecules % list ( ikind ) % array ( imol )
          molecule => molecule_set ( i )
          molecule_kind => molecule % molecule_kind
          CALL get_molecule_kind ( molecule_kind, ncolv=ncolv )

          ! Updating TARGETS for Collective Variables only
          IF ( ncolv%ntot /= 0 ) CALL shake_update_colv_int(molecule, dt, motion_section)
       END DO
    END DO MOL
    ! Intermolecular constraints
    IF (do_ext_constraint) THEN
       ! Collective Variables
       IF (gci%ncolv%ntot/= 0) CALL shake_update_colv_ext(gci, dt, motion_section)
    END IF
    CALL timestop(handle)
  END SUBROUTINE shake_update_targets

END MODULE constraint
