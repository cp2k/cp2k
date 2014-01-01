!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Module that handles the COLLECTIVE constraints
!> \par History
!>      none
!> \author Teodoro Laino [tlaino]
! *****************************************************************************
MODULE constraint_clv
  USE cell_types,                      ONLY: cell_type
  USE colvar_methods,                  ONLY: colvar_eval_mol_f
  USE colvar_types,                    ONLY: colvar_type,&
                                             diff_colvar
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set
  USE kinds,                           ONLY: dp
  USE mathlib,                         ONLY: matvec_3x3
  USE molecule_kind_types,             ONLY: colvar_constraint_type,&
                                             fixd_constraint_type,&
                                             get_molecule_kind,&
                                             molecule_kind_type
  USE molecule_types_new,              ONLY: get_molecule,&
                                             global_constraint_type,&
                                             local_colvar_constraint_type,&
                                             molecule_type
  USE particle_types,                  ONLY: particle_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: shake_roll_colv_int,&
            rattle_roll_colv_int,&
            shake_colv_int,&
            rattle_colv_int,&
            shake_roll_colv_ext,&
            rattle_roll_colv_ext,&
            shake_colv_ext,&
            rattle_colv_ext,&
            shake_update_colv_int,&
            shake_update_colv_ext

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'constraint_clv'

CONTAINS

! *****************************************************************************
!> \brief Intramolecular subroutine
!>      shake_colv algorithm for collective variables constraints
!>      updates the multiplier one molecule type at a time
!> \par History
!>      none
!> \author Teodoro Laino [tlaino]
! *****************************************************************************
  SUBROUTINE shake_colv_int( molecule, particle_set, pos, vel, dt, ishake,&
       cell,  imass, max_sigma, error )

    TYPE(molecule_type), POINTER             :: molecule
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: pos( :, : ), vel( :, : )
    REAL(kind=dp), INTENT(in)                :: dt
    INTEGER, INTENT(IN)                      :: ishake
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(:)              :: imass
    REAL(KIND=dp), INTENT(INOUT)             :: max_sigma
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'shake_colv_int', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(colvar_constraint_type), POINTER    :: colv_list( : )
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(local_colvar_constraint_type), &
      POINTER                                :: lcolv( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

    NULLIFY(fixd_list)
    failure = .FALSE.
    molecule_kind => molecule % molecule_kind
    CALL get_molecule_kind ( molecule_kind, colv_list=colv_list, fixd_list=fixd_list )
    CALL get_molecule ( molecule, lcolv=lcolv )
    ! Real Shake
    CALL shake_colv_low( fixd_list, colv_list, lcolv, &
       particle_set, pos, vel, dt, ishake, cell,  imass, max_sigma,&
       error )

  END SUBROUTINE shake_colv_int

! *****************************************************************************
!> \brief Intramolecular subroutine for updating the TARGET value of collective
!>        constraints
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE shake_update_colv_int( molecule, dt, motion_section, error)

    TYPE(molecule_type), POINTER             :: molecule
    REAL(kind=dp), INTENT(in)                :: dt
    TYPE(section_vals_type), POINTER         :: motion_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'shake_update_colv_int', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(colvar_constraint_type), POINTER    :: colv_list( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

    failure = .FALSE.
    molecule_kind => molecule % molecule_kind
    CALL get_molecule_kind ( molecule_kind, colv_list=colv_list)
    ! Real update of the Shake target
    CALL shake_update_colv_low( colv_list, dt, motion_section, error)

  END SUBROUTINE shake_update_colv_int

! *****************************************************************************
!> \brief Intramolecular subroutine
!>      rattle algorithm for collective variables constraints
!>      updates the multiplier one molecule type at a time
!> \par History
!>      none
!> \author Teodoro Laino [tlaino]
! *****************************************************************************
  SUBROUTINE rattle_colv_int( molecule, particle_set, vel, dt, irattle,&
       cell, imass, max_sigma, error )

    TYPE(molecule_type), POINTER             :: molecule
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: vel( :, : )
    REAL(kind=dp), INTENT(in)                :: dt
    INTEGER, INTENT(IN)                      :: irattle
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(:)              :: imass
    REAL(KIND=dp), INTENT(INOUT)             :: max_sigma
    TYPE(cp_error_type), INTENT(inout)       :: error

    TYPE(colvar_constraint_type), POINTER    :: colv_list( : )
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(local_colvar_constraint_type), &
      POINTER                                :: lcolv( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

    NULLIFY(fixd_list)
    molecule_kind => molecule % molecule_kind
    CALL get_molecule_kind ( molecule_kind, colv_list = colv_list, fixd_list=fixd_list )
    CALL get_molecule ( molecule, lcolv=lcolv )
    ! Real Rattle
    CALL rattle_colv_low ( fixd_list, colv_list, lcolv, &
       particle_set, vel, dt, irattle, cell, imass, max_sigma, error )

  END SUBROUTINE rattle_colv_int

! *****************************************************************************
!> \brief Intramolecular subroutine
!>      shake algorithm (box allowed to change) for collective variables constraints
!>      updates the multiplier one molecule type at a time
!> \par History
!>      none
!> \author Teodoro Laino [tlaino]
! *****************************************************************************
  SUBROUTINE shake_roll_colv_int( molecule, particle_set, pos, vel, r_shake, v_shake, &
       dt, ishake, cell, imass, max_sigma, error )

    TYPE(molecule_type), POINTER             :: molecule
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: pos( :, : ), vel( :, : )
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: r_shake, v_shake
    REAL(kind=dp), INTENT(in)                :: dt
    INTEGER, INTENT(in)                      :: ishake
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(:)              :: imass
    REAL(KIND=dp), INTENT(INOUT)             :: max_sigma
    TYPE(cp_error_type), INTENT(inout)       :: error

    TYPE(colvar_constraint_type), POINTER    :: colv_list( : )
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(local_colvar_constraint_type), &
      POINTER                                :: lcolv( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

    NULLIFY(fixd_list)
    molecule_kind => molecule % molecule_kind
    CALL get_molecule_kind ( molecule_kind, colv_list = colv_list, fixd_list=fixd_list )
    CALL get_molecule ( molecule, lcolv=lcolv )
    ! Real Shake
    CALL shake_roll_colv_low( fixd_list, colv_list, lcolv, &
       particle_set, pos, vel, r_shake, v_shake, dt, ishake, cell,&
       imass, max_sigma, error )

  END SUBROUTINE shake_roll_colv_int

! *****************************************************************************
!> \brief Intramolecular subroutine
!>      rattle algorithm (box allowed to change) for collective variables constraints
!>      updates the multiplier one molecule type at a time
!> \par History
!>      none
!> \author Teodoro Laino [tlaino]
! *****************************************************************************
  SUBROUTINE rattle_roll_colv_int ( molecule, particle_set, vel, r_rattle,  &
       dt, irattle, veps,  cell, imass, max_sigma, error )

    TYPE(molecule_type), POINTER             :: molecule
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: vel( :, : )
    REAL(KIND=dp), INTENT(IN)                :: r_rattle( :, : ), dt
    INTEGER, INTENT(in)                      :: irattle
    REAL(KIND=dp), INTENT(IN)                :: veps( :, : )
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(:)              :: imass
    REAL(KIND=dp), INTENT(INOUT)             :: max_sigma
    TYPE(cp_error_type), INTENT(inout)       :: error

    TYPE(colvar_constraint_type), POINTER    :: colv_list( : )
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(local_colvar_constraint_type), &
      POINTER                                :: lcolv( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

    NULLIFY(fixd_list)
    molecule_kind => molecule % molecule_kind
    CALL get_molecule_kind ( molecule_kind, colv_list = colv_list, fixd_list=fixd_list )
    CALL get_molecule ( molecule, lcolv=lcolv )
    ! Real Rattle
    CALL rattle_roll_colv_low (fixd_list, colv_list, lcolv, &
       particle_set, vel, r_rattle, dt, irattle, veps,  cell,&
       imass, max_sigma, error )

  END SUBROUTINE rattle_roll_colv_int

! *****************************************************************************
!> \brief Intermolecular subroutine
!>      shake_colv algorithm for collective variables constraints
!>      updates the multiplier one molecule type at a time
!> \par History
!>      none
!> \author Teodoro Laino [tlaino]
! *****************************************************************************
  SUBROUTINE shake_colv_ext( gci, particle_set, pos, vel, dt, ishake,&
       cell,  imass, max_sigma, error )

    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: pos( :, : ), vel( :, : )
    REAL(kind=dp), INTENT(in)                :: dt
    INTEGER, INTENT(IN)                      :: ishake
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(:)              :: imass
    REAL(KIND=dp), INTENT(INOUT)             :: max_sigma
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'shake_colv_ext', &
      routineP = moduleN//':'//routineN

    TYPE(colvar_constraint_type), POINTER    :: colv_list( : )
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(local_colvar_constraint_type), &
      POINTER                                :: lcolv( : )

    colv_list => gci%colv_list
    fixd_list => gci%fixd_list
    lcolv => gci%lcolv
    ! Real Shake
    CALL shake_colv_low( fixd_list, colv_list, lcolv, &
       particle_set, pos, vel, dt, ishake, cell,  imass, max_sigma,&
       error )

  END SUBROUTINE shake_colv_ext

! *****************************************************************************
!> \brief Intermolecular subroutine for updating the TARGET value for collective
!>        constraints
!> \author Teodoro Laino [tlaino]
! *****************************************************************************
  SUBROUTINE shake_update_colv_ext( gci, dt, motion_section, error)

    TYPE(global_constraint_type), POINTER    :: gci
    REAL(kind=dp), INTENT(in)                :: dt
    TYPE(section_vals_type), POINTER         :: motion_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'shake_update_colv_ext', &
      routineP = moduleN//':'//routineN

    TYPE(colvar_constraint_type), POINTER    :: colv_list( : )

    colv_list => gci%colv_list
    ! Real update of the Shake target
    CALL shake_update_colv_low( colv_list, dt, motion_section, error)

  END SUBROUTINE shake_update_colv_ext

! *****************************************************************************
!> \brief Intermolecular subroutine
!>      rattle algorithm for collective variables constraints
!>      updates the multiplier one molecule type at a time
!> \par History
!>      none
!> \author Teodoro Laino [tlaino]
! *****************************************************************************
  SUBROUTINE rattle_colv_ext( gci, particle_set, vel, dt, irattle,&
       cell, imass, max_sigma, error )

    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: vel( :, : )
    REAL(kind=dp), INTENT(in)                :: dt
    INTEGER, INTENT(IN)                      :: irattle
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(:)              :: imass
    REAL(KIND=dp), INTENT(INOUT)             :: max_sigma
    TYPE(cp_error_type), INTENT(inout)       :: error

    TYPE(colvar_constraint_type), POINTER    :: colv_list( : )
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(local_colvar_constraint_type), &
      POINTER                                :: lcolv( : )

    colv_list => gci%colv_list
    fixd_list => gci%fixd_list
    lcolv => gci%lcolv
    ! Real Rattle
    CALL rattle_colv_low ( fixd_list, colv_list, lcolv, &
       particle_set, vel, dt, irattle, cell, imass, max_sigma, error )

  END SUBROUTINE rattle_colv_ext

! *****************************************************************************
!> \brief Intermolecular subroutine
!>      shake algorithm (box allowed to change) for collective variables constraints
!>      updates the multiplier one molecule type at a time
!> \par History
!>      none
!> \author Teodoro Laino [tlaino]
! *****************************************************************************
  SUBROUTINE shake_roll_colv_ext( gci, particle_set, pos, vel, r_shake, v_shake, &
       dt, ishake, cell, imass, max_sigma, error )

    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: pos( :, : ), vel( :, : )
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: r_shake, v_shake
    REAL(kind=dp), INTENT(in)                :: dt
    INTEGER, INTENT(in)                      :: ishake
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(:)              :: imass
    REAL(KIND=dp), INTENT(INOUT)             :: max_sigma
    TYPE(cp_error_type), INTENT(inout)       :: error

    TYPE(colvar_constraint_type), POINTER    :: colv_list( : )
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(local_colvar_constraint_type), &
      POINTER                                :: lcolv( : )

    colv_list => gci%colv_list
    fixd_list => gci%fixd_list
    lcolv => gci%lcolv
    ! Real Shake
    CALL shake_roll_colv_low( fixd_list, colv_list, lcolv, &
       particle_set, pos, vel, r_shake, v_shake, dt, ishake, cell,&
       imass, max_sigma, error )

  END SUBROUTINE shake_roll_colv_ext

! *****************************************************************************
!> \brief Intermolecular subroutine
!>      rattle algorithm (box allowed to change) for collective variables constraints
!>      updates the multiplier one molecule type at a time
!> \par History
!>      none
!> \author Teodoro Laino [tlaino]
! *****************************************************************************
  SUBROUTINE rattle_roll_colv_ext ( gci, particle_set, vel, r_rattle,  &
       dt, irattle, veps,  cell, imass, max_sigma, error )

    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: vel( :, : )
    REAL(KIND=dp), INTENT(IN)                :: r_rattle( :, : ), dt
    INTEGER, INTENT(in)                      :: irattle
    REAL(KIND=dp), INTENT(IN)                :: veps( :, : )
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(:)              :: imass
    REAL(KIND=dp), INTENT(INOUT)             :: max_sigma
    TYPE(cp_error_type), INTENT(inout)       :: error

    TYPE(colvar_constraint_type), POINTER    :: colv_list( : )
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(local_colvar_constraint_type), &
      POINTER                                :: lcolv( : )

    colv_list => gci%colv_list
    fixd_list => gci%fixd_list
    lcolv => gci%lcolv
    ! Real Rattle
    CALL rattle_roll_colv_low (fixd_list, colv_list, lcolv, &
       particle_set, vel, r_rattle, dt, irattle, veps,  cell,&
       imass, max_sigma, error )

  END SUBROUTINE rattle_roll_colv_ext

! *****************************************************************************
!> \brief Real Shake subroutine - Low Level
!>      shake_colv algorithm for collective variables constraints
!>      updates the multiplier one molecule type at a time
!> \par History
!>      none
!> \author Teodoro Laino [tlaino]
! *****************************************************************************
  SUBROUTINE shake_colv_low( fixd_list, colv_list, lcolv, &
       particle_set, pos, vel, dt, ishake, cell,  imass, max_sigma,&
       error )
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(colvar_constraint_type), POINTER    :: colv_list( : )
    TYPE(local_colvar_constraint_type), &
      POINTER                                :: lcolv( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: pos( :, : ), vel( :, : )
    REAL(kind=dp), INTENT(in)                :: dt
    INTEGER, INTENT(IN)                      :: ishake
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(:)              :: imass
    REAL(KIND=dp), INTENT(INOUT)             :: max_sigma
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'shake_colv_low', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iconst
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: del_lam, dtby2, dtsqby2, &
                                                fdotf_sum

    failure = .FALSE.
    dtsqby2 = dt*dt*.5_dp
    dtby2   = dt*.5_dp
    IF (ishake==1) THEN
       DO iconst = 1, SIZE(colv_list)
          IF (colv_list(iconst)%restraint%active) CYCLE
          ! Update positions
          CALL update_con_colv(pos, dtsqby2, lcolv(iconst), &
               lambda=lcolv(iconst)%lambda,&
               imass=imass, error=error)
          ! Update velocities
          CALL update_con_colv(vel, dtby2, lcolv(iconst), &
               lambda=lcolv(iconst)%lambda,&
               imass=imass, error=error)
       END DO
    ELSE
       DO iconst = 1, SIZE(colv_list)
          IF (colv_list(iconst)%restraint%active) CYCLE
          ! Update colvar
          CALL colvar_eval_mol_f( lcolv ( iconst ) % colvar, cell, particles=particle_set,&
               pos=pos, fixd_list=fixd_list, error=error)
          lcolv ( iconst ) % sigma = diff_colvar(lcolv(iconst)%colvar,&
               colv_list(iconst)%expected_value)
          fdotf_sum = eval_Jac_colvar(lcolv ( iconst ) % colvar,&
               lcolv ( iconst ) % colvar_old, imass=imass, error=error)
          del_lam = 2.0_dp*lcolv ( iconst ) % sigma/(dt*dt*fdotf_sum)
          lcolv ( iconst ) % lambda = lcolv ( iconst ) % lambda + del_lam

          ! Update positions
          CALL update_con_colv(pos, dtsqby2, lcolv(iconst), &
               lambda=del_lam,&
               imass=imass, error=error)
          ! Update velocities
          CALL update_con_colv(vel, dtby2, lcolv(iconst), &
               lambda=del_lam,&
               imass=imass, error=error)
       END DO
    END IF
    ! computing the constraint and value of tolerance
    DO iconst = 1, SIZE(colv_list)
       IF (colv_list(iconst)%restraint%active) CYCLE
       CALL colvar_eval_mol_f( lcolv ( iconst ) % colvar, cell, particles=particle_set,&
            pos=pos, fixd_list=fixd_list, error=error)
       lcolv ( iconst ) % sigma = diff_colvar(lcolv ( iconst ) % colvar,&
            colv_list ( iconst ) % expected_value)
       max_sigma = MAX(ABS(lcolv ( iconst ) % sigma),max_sigma)
    END DO
  END SUBROUTINE shake_colv_low

! *****************************************************************************
!> \brief Real Shake subroutine - Low Level - for updating the TARGET value
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date 02.2008
! *****************************************************************************
  SUBROUTINE shake_update_colv_low(colv_list, dt, motion_section, error )
    TYPE(colvar_constraint_type), POINTER    :: colv_list( : )
    REAL(kind=dp), INTENT(in)                :: dt
    TYPE(section_vals_type), POINTER         :: motion_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'shake_update_colv_low', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iconst, irep, n_rep
    LOGICAL                                  :: do_update_colvar, explicit, &
                                                failure
    REAL(KIND=dp)                            :: clv_target, limit, &
                                                new_clv_target, value
    TYPE(section_vals_type), POINTER         :: collective_sections

    failure = .FALSE.
    ! Update globally for restart
    collective_sections => section_vals_get_subs_vals(motion_section,"CONSTRAINT%COLLECTIVE",error=error)
    CALL section_vals_get(collective_sections, n_repetition=n_rep, error=error)
    IF (n_rep/=0) THEN
       DO irep = 1, n_rep
          CALL section_vals_val_get(collective_sections, "TARGET_GROWTH", r_val=value,&
             i_rep_section=irep, error=error)
          IF ( value /= 0.0_dp ) THEN
             CALL section_vals_val_get(collective_sections, "TARGET", r_val=clv_target,&
                i_rep_section=irep, error=error)
             new_clv_target = clv_target + value * dt
             ! Check limits..
             CALL section_vals_val_get(collective_sections, "TARGET_LIMIT", explicit=explicit,&
                i_rep_section=irep, error=error)
             do_update_colvar = .TRUE.
             IF (explicit) THEN
                CALL section_vals_val_get(collective_sections, "TARGET_LIMIT", r_val=limit,&
                   i_rep_section=irep, error=error)
                IF (value>0.0_dp) THEN
                   IF (clv_target==limit) THEN
                      do_update_colvar = .FALSE.
                   ELSE IF (new_clv_target>=limit) THEN
                      new_clv_target = limit
                   END IF
                ELSE
                   IF (clv_target==limit) THEN
                      do_update_colvar = .FALSE.
                   ELSE IF (new_clv_target<=limit) THEN
                      new_clv_target = limit
                   END IF
                ENDIF
             END IF
             IF (do_update_colvar) THEN
                CALL section_vals_val_set(collective_sections, "TARGET", r_val=new_clv_target,&
                   i_rep_section=irep, error=error)
             END IF
          END IF
       END DO
    END IF

    ! Update locally the value to each processor
    DO iconst = 1, SIZE(colv_list)
       ! Update local to each processor
       IF (colv_list(iconst)%expected_value_growth_speed == 0.0_dp) CYCLE
       CALL section_vals_val_get(collective_sections, "TARGET", &
             r_val=colv_list(iconst)%expected_value,&
             i_rep_section=colv_list(iconst)%inp_seq_num, error=error)
    END DO

  END SUBROUTINE shake_update_colv_low

! *****************************************************************************
!> \brief Real Rattle - Low Level
!>      rattle algorithm for collective variables constraints
!>      updates the multiplier one molecule type at a time
!> \par History
!>      none
!> \author Teodoro Laino [tlaino]
! *****************************************************************************
  SUBROUTINE rattle_colv_low ( fixd_list, colv_list, lcolv, &
       particle_set, vel, dt, irattle, cell, imass, max_sigma, error )

    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(colvar_constraint_type), POINTER    :: colv_list( : )
    TYPE(local_colvar_constraint_type), &
      POINTER                                :: lcolv( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: vel( :, : )
    REAL(kind=dp), INTENT(in)                :: dt
    INTEGER, INTENT(IN)                      :: irattle
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(:)              :: imass
    REAL(KIND=dp), INTENT(INOUT)             :: max_sigma
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: iconst
    REAL(KIND=dp)                            :: del_lam, dtby2, fdotf_sum

    dtby2     = dt*.5_dp
    IF (irattle==1) THEN
       DO iconst = 1, SIZE(colv_list)
          IF (colv_list(iconst)%restraint%active) CYCLE
          ! Update colvar_old
          CALL colvar_eval_mol_f(lcolv ( iconst ) % colvar_old, cell,&
               particles=particle_set, fixd_list=fixd_list, error=error)
          ! Update velocities
          CALL update_con_colv(vel, dtby2, lcolv(iconst), &
               lambda=lcolv(iconst)%lambda,&
               imass=imass,error=error)
       END DO
    ELSE
       DO iconst = 1, SIZE(colv_list)
          IF (colv_list(iconst)%restraint%active) CYCLE
          lcolv ( iconst ) % sigma = rattle_con_eval(lcolv ( iconst ) % colvar_old, vel,error=error)
          fdotf_sum = eval_Jac_colvar(lcolv ( iconst ) % colvar_old,&
               lcolv ( iconst ) % colvar_old, imass=imass,error=error)
          del_lam = 2.0_dp*lcolv ( iconst ) % sigma/(dt*fdotf_sum)
          lcolv ( iconst ) % lambda = lcolv ( iconst ) % lambda + del_lam

          ! Update velocities
          CALL update_con_colv(vel, dtby2, lcolv(iconst), &
               lambda=del_lam,&
               imass=imass,error=error)
       END DO
    END IF

    DO iconst = 1, SIZE(colv_list)
       IF (colv_list(iconst)%restraint%active) CYCLE
       lcolv ( iconst ) % sigma = rattle_con_eval(lcolv ( iconst ) % colvar_old, vel,error=error)
       max_sigma = MAX(ABS(lcolv ( iconst ) % sigma),max_sigma)
    END DO

  END SUBROUTINE rattle_colv_low

! *****************************************************************************
!> \brief Real shake_roll - Low Level
!>      shake algorithm (box allowed to change) for collective variables constraints
!>      updates the multiplier one molecule type at a time
!> \par History
!>      none
!> \author Teodoro Laino [tlaino]
! *****************************************************************************
  SUBROUTINE shake_roll_colv_low( fixd_list, colv_list, lcolv, &
       particle_set, pos, vel, r_shake, v_shake, dt, ishake, cell,&
       imass, max_sigma, error )

    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(colvar_constraint_type), POINTER    :: colv_list( : )
    TYPE(local_colvar_constraint_type), &
      POINTER                                :: lcolv( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: pos( :, : ), vel( :, : )
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: r_shake, v_shake
    REAL(kind=dp), INTENT(in)                :: dt
    INTEGER, INTENT(in)                      :: ishake
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(:)              :: imass
    REAL(KIND=dp), INTENT(INOUT)             :: max_sigma
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: iconst
    REAL(KIND=dp)                            :: del_lam, dtby2, dtsqby2, &
                                                fdotf_sum

    dtsqby2 = dt*dt*.5_dp
    dtby2 = dt*.5_dp
    IF (ishake==1) THEN
       DO iconst = 1, SIZE(colv_list)
          IF (colv_list(iconst)%restraint%active) CYCLE
          ! Update positions
          CALL update_con_colv(pos, dtsqby2, lcolv(iconst), &
               lambda=lcolv(iconst)%lambda,&
               roll=.TRUE.,rmat=r_shake,imass=imass,error=error)
          ! Update velocities
          CALL update_con_colv(vel, dtby2, lcolv(iconst), &
               lambda=lcolv(iconst)%lambda,&
               roll=.TRUE.,rmat=v_shake,imass=imass,error=error)
       END DO
    ELSE
       DO iconst = 1, SIZE(colv_list)
          IF (colv_list(iconst)%restraint%active) CYCLE
          ! Update colvar
          CALL colvar_eval_mol_f( lcolv ( iconst ) % colvar, cell, particles=particle_set,&
               pos=pos, fixd_list=fixd_list, error=error)
          lcolv ( iconst ) % sigma = diff_colvar(lcolv ( iconst ) % colvar,&
               colv_list ( iconst ) % expected_value)
          fdotf_sum = eval_Jac_colvar(lcolv ( iconst ) % colvar,&
               lcolv ( iconst ) % colvar_old, roll=.TRUE., rmat=r_shake,&
               imass=imass, error=error)
          del_lam = 2.0_dp*lcolv ( iconst ) % sigma/(dt*dt*fdotf_sum)
          lcolv ( iconst ) % lambda = lcolv ( iconst ) % lambda + del_lam

          ! Update positions
          CALL update_con_colv(pos, dtsqby2, lcolv(iconst), &
               lambda=del_lam,&
               roll=.TRUE., rmat=r_shake, imass=imass,error=error)
          ! Update velocities
          CALL update_con_colv(vel, dtby2, lcolv(iconst), &
               lambda=del_lam,&
               roll=.TRUE., rmat=v_shake, imass=imass,error=error)
       END DO
    END IF
    ! computing the constraint and value of tolerance
    DO iconst = 1, SIZE(colv_list)
       IF (colv_list(iconst)%restraint%active) CYCLE
       CALL colvar_eval_mol_f( lcolv ( iconst ) % colvar, cell, particles=particle_set,&
            pos=pos, fixd_list=fixd_list, error=error)
       lcolv ( iconst ) % sigma = diff_colvar(lcolv ( iconst ) % colvar,&
            colv_list ( iconst ) % expected_value)
       max_sigma = MAX(ABS(lcolv ( iconst ) % sigma),max_sigma)
    END DO

  END SUBROUTINE shake_roll_colv_low

! *****************************************************************************
!> \brief Real Rattle_roll - Low Level
!>      rattle algorithm (box allowed to change) for collective variables constraints
!>      updates the multiplier one molecule type at a time
!> \par History
!>      none
!> \author Teodoro Laino [tlaino]
! *****************************************************************************
  SUBROUTINE rattle_roll_colv_low (fixd_list, colv_list, lcolv, &
       particle_set, vel, r_rattle, dt, irattle, veps,  cell,&
       imass, max_sigma, error )

    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(colvar_constraint_type), POINTER    :: colv_list( : )
    TYPE(local_colvar_constraint_type), &
      POINTER                                :: lcolv( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: vel( :, : )
    REAL(KIND=dp), INTENT(IN)                :: r_rattle( :, : ), dt
    INTEGER, INTENT(in)                      :: irattle
    REAL(KIND=dp), INTENT(IN)                :: veps( :, : )
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(:)              :: imass
    REAL(KIND=dp), INTENT(INOUT)             :: max_sigma
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: iconst
    REAL(KIND=dp)                            :: del_lam, dtby2, fdotf_sum

    dtby2 = dt*.5_dp
    IF (irattle==1) THEN
       DO iconst = 1, SIZE ( colv_list )
          IF (colv_list(iconst)%restraint%active) CYCLE
          ! Update colvar_old
          CALL colvar_eval_mol_f(lcolv ( iconst ) % colvar_old, cell,&
               particles=particle_set, fixd_list=fixd_list, error=error)
          ! Update velocities
          CALL update_con_colv(vel, dtby2, lcolv(iconst), &
               lambda=lcolv(iconst)%lambda,&
               imass=imass,error=error)
       END DO
    ELSE
       DO iconst = 1, SIZE ( colv_list )
          IF (colv_list(iconst)%restraint%active) CYCLE
          lcolv ( iconst ) % sigma = rattle_con_eval(lcolv ( iconst ) % colvar_old, vel,&
               roll=.TRUE., veps=veps, rmat=r_rattle, particles=particle_set,&
               error=error)
          fdotf_sum = eval_Jac_colvar(lcolv ( iconst ) % colvar_old,&
               lcolv ( iconst ) % colvar_old, roll=.TRUE.,&
               rmat=r_rattle,imass=imass, error=error)
          del_lam = 2.0_dp*lcolv ( iconst ) % sigma/(dt*fdotf_sum)
          lcolv ( iconst ) % lambda = lcolv ( iconst ) % lambda + del_lam
          ! Update velocities
          CALL update_con_colv(vel, dtby2, lcolv(iconst), &
               lambda=del_lam,&
               roll=.TRUE.,rmat=r_rattle,imass=imass,error=error)
       END DO
    END IF
    ! computing the constraint and value of the tolerance
    DO iconst = 1, SIZE(colv_list)
       IF (colv_list(iconst)%restraint%active) CYCLE
       lcolv ( iconst ) % sigma = rattle_con_eval(lcolv ( iconst ) % colvar_old, vel,&
            roll=.TRUE., veps=veps, rmat=r_rattle, particles=particle_set,&
            error=error)
       max_sigma = MAX(ABS(lcolv ( iconst ) % sigma),max_sigma)
    END DO

  END SUBROUTINE rattle_roll_colv_low

! *****************************************************************************
!> \brief Update position/velocities
!> \par History
!>      Teodoro Laino [teo] created 04.2006
!> \author Teodoro Laino [tlaino]
! *****************************************************************************
  SUBROUTINE  update_con_colv(wrk, fac, lcolv, lambda, roll, rmat, imass, error)
    REAL(KIND=dp), INTENT(INOUT)             :: wrk( :, : )
    REAL(KIND=dp), INTENT(IN)                :: fac
    TYPE(local_colvar_constraint_type), &
      INTENT(IN)                             :: lcolv
    REAL(KIND=dp), INTENT(IN)                :: lambda
    LOGICAL, INTENT(in), OPTIONAL            :: roll
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN), OPTIONAL                   :: rmat
    REAL(KIND=dp), DIMENSION(:)              :: imass
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'update_con_colv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iatm, ind
    LOGICAL                                  :: failure, my_roll
    REAL(KIND=dp), DIMENSION(3)              :: f_roll

    failure = .FALSE.
    my_roll = .FALSE.
    IF (PRESENT(roll)) THEN
       my_roll = roll
       IF (my_roll) THEN
          CPPostcondition(PRESENT(rmat),cp_failure_level,routineP,error,failure)
       END IF
    END IF
    IF (.NOT.failure) THEN
       DO iatm = 1, SIZE( lcolv % colvar_old % i_atom )
          ind = lcolv % colvar_old % i_atom (iatm)
          !
          IF (my_roll) THEN
             ! If ROLL rotate forces
             CALL matvec_3x3 (f_roll, rmat, lcolv % colvar_old % dsdr(:,iatm)  )
          ELSE
             f_roll = lcolv % colvar_old % dsdr(:,iatm)
          END IF
          wrk ( :, ind ) = wrk ( :, ind ) - imass(ind) * fac * lambda * f_roll
       END DO
    END IF
  END SUBROUTINE update_con_colv

! *****************************************************************************
!> \brief Evaluates the Jacobian of the collective variables constraints
!> \par History
!>      Teodoro Laino [teo] created 04.2006
!> \author Teodoro Laino [tlaino]
! *****************************************************************************
  FUNCTION eval_Jac_colvar(colvar, colvar_old, roll, rmat, imass, error) RESULT(res)
    TYPE(colvar_type), POINTER               :: colvar, colvar_old
    LOGICAL, INTENT(IN), OPTIONAL            :: roll
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN), OPTIONAL                   :: rmat
    REAL(KIND=dp), DIMENSION(:)              :: imass
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'eval_Jac_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iatom
    LOGICAL                                  :: failure, my_roll
    REAL(KIND=dp), DIMENSION(3)              :: tmp1, tmp2, tmp3

    failure = .FALSE.
    my_roll = .FALSE.
    IF (PRESENT(roll)) THEN
       my_roll = roll
       IF (my_roll) THEN
          CPPostcondition(PRESENT(rmat),cp_failure_level,routineP,error,failure)
       END IF
    END IF

    res = 0.0_dp
    DO i = 1, SIZE(colvar%i_atom)
       iatom = colvar%i_atom(i)
       tmp1 = colvar%dsdr(1:3,i)
       IF (my_roll) THEN
          tmp3 = colvar_old%dsdr(1:3,i)
          CALL  matvec_3x3 (tmp2, rmat, tmp3 )
       ELSE
          tmp2 = colvar_old%dsdr(1:3,i)
       END IF
       res = res + DOT_PRODUCT(tmp1,tmp2)*imass(iatom)
    END DO

  END FUNCTION eval_Jac_colvar

! *****************************************************************************
!> \brief Evaluates the constraint for the rattle scheme
!> \par History
!>      Teodoro Laino [teo] created 04.2006
!> \author Teodoro Laino [tlaino]
! *****************************************************************************
  FUNCTION rattle_con_eval(colvar, vel, roll, veps, rmat, particles, error) RESULT(res)
    TYPE(colvar_type), POINTER               :: colvar
    REAL(KIND=dp), INTENT(INOUT)             :: vel( :, : )
    LOGICAL, INTENT(IN), OPTIONAL            :: roll
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN), OPTIONAL                   :: veps, rmat
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particles
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'rattle_con_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iatm, ind
    LOGICAL                                  :: failure, my_roll
    REAL(KIND=dp)                            :: tmp(3)
    REAL(KIND=dp), DIMENSION(3)              :: f_roll, pos, v_roll

    failure = .FALSE.
    my_roll = .FALSE.
    IF (PRESENT(roll)) THEN
       my_roll = roll
       IF (my_roll) THEN
          CPPostcondition(PRESENT(rmat),cp_failure_level,routineP,error,failure)
          CPPostcondition(PRESENT(veps),cp_failure_level,routineP,error,failure)
          CPPostcondition(PRESENT(particles),cp_failure_level,routineP,error,failure)
       END IF
    END IF
    res = 0.0_dp
    DO iatm = 1, SIZE( colvar % i_atom )
       ind = colvar % i_atom (iatm)
       IF (my_roll) THEN
          pos    = particles ( ind ) % r
          CALL matvec_3x3 (f_roll, rmat, colvar % dsdr(:,iatm) )
          CALL matvec_3x3 (tmp, veps, pos )
          v_roll = vel(:,ind) + tmp
       ELSE
          f_roll = colvar % dsdr(:,iatm)
          v_roll = vel(:,ind)
       END IF
       res = res + DOT_PRODUCT(f_roll, v_roll)
    END DO

  END FUNCTION rattle_con_eval

END MODULE constraint_clv
