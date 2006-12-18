!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2006  CP2K developers group                          !
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
  USE cell_types,                      ONLY: cell_type
  USE colvar_methods,                  ONLY: colvar_eval_mol_f
  USE colvar_types,                    ONLY: colvar_counters,&
                                             colvar_type
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_p_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE kinds,                           ONLY: dp
  USE mathlib,                         ONLY: matvec_3x3
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE molecule_kind_types,             ONLY: colvar_constraint_type,&
                                             fixd_constraint_type,&
                                             get_molecule_kind,&
                                             molecule_kind_type
  USE molecule_types_new,              ONLY: get_molecule,&
                                             local_colvar_constraint_type,&
                                             molecule_type
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type,&
                                             update_particle_set
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: shake_roll_colv,&
            rattle_roll_colv,&
            fix_atom_control,&
            shake_colv,&
            rattle_colv

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
!!     - optionally apply fix atom constraint to random forces (Langevin)
!!       (04.10.206,MK)
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE fix_atom_control( force_env, error, w)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout)  &
                                             :: error
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL                               :: w

    CHARACTER(len=*), PARAMETER :: routineN = 'fix_atom_control', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, ii, ikind, iparticle, iparticle_local, isubsys, &
      my_atm_fixed, nfixed_atoms, nkind, nparticle_local, stat
    LOGICAL                                  :: failure, zero_force
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: force
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
       CALL force_env_get(force_env=force_env, subsys=subsys,error=error)
       DO isubsys = 1, SIZE(subsys)
          CALL cp_subsys_get(subsys=subsys(isubsys)%subsys,&
               particles=particles, &
               local_particles=local_particles,&
               molecule_kinds_new=molecule_kinds ,error=error)

          nkind             =  molecule_kinds % n_els
          molecule_kind_set => molecule_kinds % els
          particle_set      => particles%els
          my_atm_fixed      =  0
          DO ikind = 1, nkind
             molecule_kind => molecule_kind_set(ikind)
             CALL get_molecule_kind ( molecule_kind, nfixd=nfixed_atoms )
             my_atm_fixed = my_atm_fixed + nfixed_atoms
          END DO
          IF (my_atm_fixed /=0) THEN
             IF (.NOT.PRESENT(w)) THEN
                ! Allocate scratch array
                ALLOCATE(force(3,SIZE(particles%els)),stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                force = 0.0_dp
                DO i = 1, SIZE(local_particles%n_el)
                   nparticle_local   = local_particles%n_el(i)
                   DO iparticle_local=1,nparticle_local
                      iparticle = local_particles%list(i)%array(iparticle_local)
                      force(:,iparticle)  = particle_set(iparticle)%f(:)
                   END DO
                END DO
             END IF

             Mol_kind:  DO ikind = 1, nkind
                molecule_kind => molecule_kind_set(ikind)
                CALL get_molecule_kind ( molecule_kind, nfixd=nfixed_atoms, fixd_list=fixd_list )

                IF (nfixed_atoms ==0) CYCLE Mol_Kind
                DO i = 1, SIZE(local_particles%n_el)
                   nparticle_local   = local_particles%n_el(i)
                   DO iparticle_local=1,nparticle_local
                      iparticle = local_particles%list(i)%array(iparticle_local)
                      DO ii = 1, SIZE(fixd_list)
                         IF (fixd_list(ii)%fixd == iparticle) THEN
                            zero_force = ((fixd_list(ii)%fixd == iparticle).AND.&
                                 (.NOT.fixd_list(ii)%restraint%active))
                            IF (zero_force) THEN
                               ! Constraint
                               IF (PRESENT(w)) THEN
                                  w(:,iparticle) = 0.0_dp
                               ELSE
                                  force(:,iparticle) = 0.0_dp
                               END IF
                               EXIT
                            END IF
                         END IF
                      END DO
                   END DO
                END DO
             END DO Mol_kind
             IF (.NOT.PRESENT(w)) THEN
                CALL update_particle_set(particle_set,force_env%para_env%group,for=force)
                DEALLOCATE(force,stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
          END IF
       END DO
    END IF
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE fix_atom_control


!!****** constraint_clv/shake_colv [1.0] *
!!
!!   NAME
!!     shake_colv
!!
!!   FUNCTION
!!     shake_colv algorithm for collective variables constraints
!!     updates the multiplier one molecule type at a time
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE shake_colv( molecule, particle_set, pos, vel, dt, shake_tol,&
       log_unit, cell, ncolv, imass, error )

    TYPE(molecule_type), POINTER             :: molecule
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: pos( :, : ), vel( :, : )
    REAL(kind=dp), INTENT(in)                :: dt, shake_tol
    INTEGER, INTENT(in)                      :: log_unit
    TYPE(cell_type), POINTER                 :: cell
    TYPE(colvar_counters), INTENT(IN)        :: ncolv
    REAL(KIND=dp), DIMENSION(:)              :: imass
    TYPE(cp_error_type), INTENT(inout)  &
                                             :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'shake_colv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: first_atom, iconst, ishake
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: del_lam, dtby2, dtsqby2, &
                                                fdotf_sum, max_sigma
    TYPE(colvar_constraint_type), POINTER    :: colv_list( : )
    TYPE(local_colvar_constraint_type), &
      POINTER                                :: lcolv( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    
    NULLIFY(fixd_list)
    failure = .FALSE.
    ishake = 0
    max_sigma = -1.0E+10_dp
    dtsqby2 = dt*dt*.5_dp
    dtby2   = dt*.5_dp
    CONS: DO WHILE (ABS(max_sigma)>=shake_tol)
       max_sigma = 0.0_dp
       ishake  = ishake + 1
       molecule_kind => molecule % molecule_kind
       CALL get_molecule_kind ( molecule_kind, colv_list=colv_list, fixd_list=fixd_list )
       CALL get_molecule ( molecule, first_atom = first_atom, lcolv=lcolv )

       IF (ishake==1) THEN
          DO iconst = 1, SIZE(colv_list)
             IF (colv_list(iconst)%restraint%active) CYCLE
             ! Update positions
             CALL update_con_colv(pos, dtsqby2, lcolv(iconst), &
                  particles=particle_set, lambda=lcolv(iconst)%lambda,&
                  imass=imass, error=error)
             ! Update velocities
             CALL update_con_colv(vel, dtby2, lcolv(iconst), &
                  particles=particle_set, lambda=lcolv(iconst)%lambda,&
                  imass=imass, error=error)
          END DO
       ELSE
          DO iconst = 1, SIZE(colv_list)
             IF (colv_list(iconst)%restraint%active) CYCLE
             ! Update colvar
             CALL colvar_eval_mol_f( lcolv ( iconst ) % colvar, cell, particles=particle_set,&
                  pos=pos, fixd_list=fixd_list, error=error)
             lcolv ( iconst ) % sigma = lcolv ( iconst ) % colvar % ss - &
                  colv_list ( iconst ) % expected_value
             fdotf_sum = eval_Jac_colvar(lcolv ( iconst ) % colvar,&
                  lcolv ( iconst ) % colvar_old, particle_set, imass=imass, error=error)
             del_lam = 2.0_dp*lcolv ( iconst ) % sigma/(dt*dt*fdotf_sum)
             lcolv ( iconst ) % lambda = lcolv ( iconst ) % lambda + del_lam

             ! Update positions
             CALL update_con_colv(pos, dtsqby2, lcolv(iconst), &
                  particles=particle_set, lambda=del_lam,&
                  imass=imass, error=error)
             ! Update velocities
             CALL update_con_colv(vel, dtby2, lcolv(iconst), &
                  particles=particle_set, lambda=del_lam,&
                  imass=imass, error=error)
          END DO
       END IF
       ! computing the constraint and value of tolerance
       DO iconst = 1, SIZE(colv_list)
          IF (colv_list(iconst)%restraint%active) CYCLE
          CALL colvar_eval_mol_f( lcolv ( iconst ) % colvar, cell, particles=particle_set,&
               pos=pos, fixd_list=fixd_list, error=error)
          lcolv ( iconst ) % sigma = lcolv ( iconst ) % colvar % ss - &
               colv_list ( iconst ) % expected_value
          max_sigma = MAX(ABS(lcolv ( iconst ) % sigma),max_sigma)
       END DO
       IF (log_unit>0) &
            WRITE(log_unit,'(T2,"SHAKE : ",I5,2F15.9)')ishake,ABS(max_sigma),&
            shake_tol
    END DO CONS

  END SUBROUTINE shake_colv

!!****** constraint_clv/rattle_colv [1.0] *
!!
!!   NAME
!!     rattle_colv
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
  SUBROUTINE rattle_colv ( molecule, particle_set, vel, dt, shake_tol,&
       log_unit, ncolv, cell, imass, error )

    TYPE(molecule_type), POINTER             :: molecule
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: vel( :, : )
    REAL(kind=dp), INTENT(in)                :: dt, shake_tol
    INTEGER, INTENT(in)                      :: log_unit
    TYPE(colvar_counters), INTENT(IN)        :: ncolv
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(:)              :: imass
    TYPE(cp_error_type), INTENT(inout)  &
                                             :: error

    INTEGER                                  :: first_atom, iconst, irattle
    REAL(KIND=dp)                            :: del_lam, dtby2, fdotf_sum, &
                                                max_sigma
    TYPE(colvar_constraint_type), POINTER    :: colv_list( : )
    TYPE(local_colvar_constraint_type), &
      POINTER                                :: lcolv( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    
    NULLIFY(fixd_list)
    irattle   = 0
    max_sigma = -1.0E+10_dp
    dtby2     = dt*.5_dp
    CONS: DO WHILE ( ABS ( max_sigma ) >=  shake_tol )
       max_sigma = 0.0_dp
       irattle = irattle + 1
       molecule_kind => molecule % molecule_kind
       CALL get_molecule_kind ( molecule_kind, colv_list = colv_list, fixd_list=fixd_list )
       CALL get_molecule ( molecule, first_atom = first_atom, lcolv=lcolv )
       IF (irattle==1) THEN
          DO iconst = 1, SIZE(colv_list)
             IF (colv_list(iconst)%restraint%active) CYCLE
             ! Update colvar_old
             CALL colvar_eval_mol_f(lcolv ( iconst ) % colvar_old, cell,&
                  particles=particle_set, fixd_list=fixd_list, error=error)
             ! Update velocities
             CALL update_con_colv(vel, dtby2, lcolv(iconst), &
                  particles=particle_set, lambda=lcolv(iconst)%lambda,&
                  imass=imass,error=error)
          END DO
       ELSE
          DO iconst = 1, SIZE(colv_list)
             IF (colv_list(iconst)%restraint%active) CYCLE
             lcolv ( iconst ) % sigma = rattle_con_eval(lcolv ( iconst ) % colvar_old, vel,error=error)
             fdotf_sum = eval_Jac_colvar(lcolv ( iconst ) % colvar_old,&
                  lcolv ( iconst ) % colvar_old, particle_set,imass=imass,error=error)
             del_lam = 2.0_dp*lcolv ( iconst ) % sigma/(dt*fdotf_sum)
             lcolv ( iconst ) % lambda = lcolv ( iconst ) % lambda + del_lam

             ! Update velocities
             CALL update_con_colv(vel, dtby2, lcolv(iconst), &
                  particles=particle_set, lambda=del_lam,&
                  imass=imass,error=error)
          END DO
       END IF

       DO iconst = 1, SIZE(colv_list)
          IF (colv_list(iconst)%restraint%active) CYCLE
          lcolv ( iconst ) % sigma = rattle_con_eval(lcolv ( iconst ) % colvar_old, vel,error=error)
          max_sigma = MAX(ABS(lcolv ( iconst ) % sigma),max_sigma)
       END DO
       IF (log_unit>0) &
            WRITE(log_unit,'(T2,"RATTLE: ",I5,2F15.9)')irattle,ABS(max_sigma),&
            shake_tol
    END DO CONS

  END SUBROUTINE rattle_colv

!!****** constraint_clv/shake_roll_colv [1.0] *
!!
!!   NAME
!!     shake_roll_colv
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
  SUBROUTINE shake_roll_colv( molecule, particle_set, pos, vel, r_shake, v_shake, &
       dt, shake_tol, log_unit, ncolv, cell, imass, error )

    TYPE(molecule_type), POINTER             :: molecule
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: pos( :, : ), vel( :, : )
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: r_shake, v_shake
    REAL(kind=dp), INTENT(in)                :: dt, shake_tol
    INTEGER, INTENT(in)                      :: log_unit
    TYPE(colvar_counters), INTENT(IN)        :: ncolv
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(:)              :: imass
    TYPE(cp_error_type), INTENT(inout)  &
                                             :: error

    INTEGER                                  :: first_atom, iconst, ishake
    REAL(KIND=dp)                            :: del_lam, dtby2, dtsqby2, &
                                                fdotf_sum, max_sigma
    TYPE(colvar_constraint_type), POINTER    :: colv_list( : )
    TYPE(local_colvar_constraint_type), &
      POINTER                                :: lcolv( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    
    NULLIFY(fixd_list)
    dtsqby2 = dt*dt*.5_dp
    dtby2 = dt*.5_dp
    max_sigma = -1.0E+10_dp
    ishake = 0
    CONS: DO WHILE (ABS(max_sigma)>=shake_tol)
       max_sigma = 0.0_dp
       ishake = ishake + 1
       molecule_kind => molecule % molecule_kind
       CALL get_molecule_kind ( molecule_kind, colv_list = colv_list, fixd_list=fixd_list )
       CALL get_molecule ( molecule, first_atom = first_atom, lcolv=lcolv )
       IF (ishake==1) THEN
          DO iconst = 1, SIZE(colv_list)
             IF (colv_list(iconst)%restraint%active) CYCLE
             ! Update positions
             CALL update_con_colv(pos, dtsqby2, lcolv(iconst), &
                  particles=particle_set, lambda=lcolv(iconst)%lambda,&
                  roll=.TRUE.,rmat=r_shake,imass=imass,error=error)
             ! Update velocities
             CALL update_con_colv(vel, dtby2, lcolv(iconst), &
                  particles=particle_set, lambda=lcolv(iconst)%lambda,&
                  roll=.TRUE.,rmat=v_shake,imass=imass,error=error)
          END DO
       ELSE
          DO iconst = 1, SIZE(colv_list)
             IF (colv_list(iconst)%restraint%active) CYCLE
             ! Update colvar
             CALL colvar_eval_mol_f( lcolv ( iconst ) % colvar, cell, particles=particle_set,&
                  pos=pos, fixd_list=fixd_list, error=error)
             lcolv ( iconst ) % sigma = lcolv ( iconst ) % colvar % ss - &
                  colv_list ( iconst ) % expected_value
             fdotf_sum = eval_Jac_colvar(lcolv ( iconst ) % colvar,&
                  lcolv ( iconst ) % colvar_old, particle_set, roll=.TRUE., rmat=r_shake,&
                  imass=imass, error=error)
             del_lam = 2.0_dp*lcolv ( iconst ) % sigma/(dt*dt*fdotf_sum)
             lcolv ( iconst ) % lambda = lcolv ( iconst ) % lambda + del_lam

             ! Update positions
             CALL update_con_colv(pos, dtsqby2, lcolv(iconst), &
                  particles=particle_set, lambda=del_lam,&
                  roll=.TRUE., rmat=r_shake, imass=imass,error=error)
             ! Update velocities
             CALL update_con_colv(vel, dtby2, lcolv(iconst), &
                  particles=particle_set, lambda=del_lam,&
                  roll=.TRUE., rmat=v_shake, imass=imass,error=error)
          END DO
       END IF

       ! computing the constraint and value of tolerance
       DO iconst = 1, SIZE(colv_list)
          IF (colv_list(iconst)%restraint%active) CYCLE
          CALL colvar_eval_mol_f( lcolv ( iconst ) % colvar, cell, particles=particle_set,&
               pos=pos, fixd_list=fixd_list, error=error)
          lcolv ( iconst ) % sigma = lcolv ( iconst ) % colvar % ss - &
               colv_list ( iconst ) % expected_value
          max_sigma = MAX(ABS(lcolv ( iconst ) % sigma),max_sigma)
       END DO
       IF (log_unit>0) &
            WRITE(log_unit,'(T2,"SHAKE : ",I5,2F15.9)')ishake,ABS(max_sigma),&
            shake_tol
    END DO CONS

  END SUBROUTINE shake_roll_colv

!!****** constraint_clv/rattle_roll_colv [1.0] *
!!
!!   NAME
!!     rattle_roll_colv
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
  SUBROUTINE rattle_roll_colv ( molecule, particle_set, vel, r_rattle,  &
       dt, shake_tol, log_unit, veps, ncolv, cell, imass, error )

    TYPE(molecule_type), POINTER             :: molecule
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: vel( :, : )
    REAL(KIND=dp), INTENT(IN)                :: r_rattle( :, : ), dt, &
                                                shake_tol
    INTEGER, INTENT(in)                      :: log_unit
    REAL(KIND=dp), INTENT(IN)                :: veps( :, : )
    TYPE(colvar_counters), INTENT(IN)        :: ncolv
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(:)              :: imass
    TYPE(cp_error_type), INTENT(inout)  &
                                             :: error

    INTEGER                                  :: first_atom, iconst, irattle
    REAL(KIND=dp)                            :: del_lam, dtby2, fdotf_sum, &
                                                max_sigma
    TYPE(colvar_constraint_type), POINTER    :: colv_list( : )
    TYPE(local_colvar_constraint_type), &
      POINTER                                :: lcolv( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    
    NULLIFY(fixd_list)
    dtby2 = dt*.5_dp
    max_sigma = -1.0E+10_dp
    irattle = 0
    CONS: DO WHILE (ABS(max_sigma)>=shake_tol)
       max_sigma = 0.0_dp
       irattle = irattle + 1
       molecule_kind => molecule % molecule_kind
       CALL get_molecule_kind ( molecule_kind, colv_list = colv_list, fixd_list=fixd_list )
       CALL get_molecule ( molecule, first_atom = first_atom, lcolv=lcolv )
       IF (irattle==1) THEN
          DO iconst = 1, SIZE ( colv_list )
             IF (colv_list(iconst)%restraint%active) CYCLE
             ! Update colvar_old
             CALL colvar_eval_mol_f(lcolv ( iconst ) % colvar_old, cell,&
                  particles=particle_set, fixd_list=fixd_list, error=error)
             ! Update velocities
             CALL update_con_colv(vel, dtby2, lcolv(iconst), &
                  particles=particle_set, lambda=lcolv(iconst)%lambda,&
                  imass=imass,error=error)
          END DO
       ELSE
          DO iconst = 1, SIZE ( colv_list )
             IF (colv_list(iconst)%restraint%active) CYCLE
             lcolv ( iconst ) % sigma = rattle_con_eval(lcolv ( iconst ) % colvar_old, vel,&
                  roll=.TRUE., veps=veps, rmat=r_rattle, particles=particle_set,&
                  error=error)
             fdotf_sum = eval_Jac_colvar(lcolv ( iconst ) % colvar_old,&
                  lcolv ( iconst ) % colvar_old, particle_set,roll=.TRUE.,&
                  rmat=r_rattle,imass=imass, error=error)
             del_lam = 2.0_dp*lcolv ( iconst ) % sigma/(dt*fdotf_sum)
             lcolv ( iconst ) % lambda = lcolv ( iconst ) % lambda + del_lam

             ! Update velocities
             CALL update_con_colv(vel, dtby2, lcolv(iconst), &
                  particles=particle_set, lambda=del_lam,&
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
       IF (log_unit>0) &
            WRITE(log_unit,'(T2,"RATTLE: ",I5,2F15.9)')irattle,ABS(max_sigma),&
            shake_tol
    END DO CONS
  END SUBROUTINE rattle_roll_colv

!!****** constraint_clv/update_con_colv [1.0] *
!!
!!   NAME
!!     update_con_colv
!!
!!   FUNCTION
!!     Update position/velocities
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     Teodoro Laino [teo] created 04.2006
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE  update_con_colv(wrk, fac, lcolv, particles, lambda, roll, rmat, imass,&
       error)
    REAL(KIND=dp), INTENT(INOUT)             :: wrk( :, : )
    REAL(KIND=dp), INTENT(IN)                :: fac
    TYPE(local_colvar_constraint_type), &
      INTENT(IN)                             :: lcolv
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles
    REAL(KIND=dp), INTENT(IN)                :: lambda
    LOGICAL, INTENT(in), OPTIONAL            :: roll
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN), OPTIONAL                   :: rmat
    REAL(KIND=dp), DIMENSION(:)              :: imass
    TYPE(cp_error_type), INTENT(inout)  &
                                             :: error

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

!!****** constraint_clv/eval_Jac_colvar [1.0] *
!!
!!   NAME
!!     eval_Jac_colvar
!!
!!   FUNCTION
!!     Evaluates the Jacobian of the collective variables constraints
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     Teodoro Laino [teo] created 04.2006
!!
!!   SOURCE
!******************************************************************************
  FUNCTION eval_Jac_colvar(colvar, colvar_old, particles, roll, rmat, imass,&
       error) RESULT(res)
    TYPE(colvar_type), POINTER               :: colvar, colvar_old
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles
    LOGICAL, INTENT(IN), OPTIONAL            :: roll
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN), OPTIONAL                   :: rmat
    REAL(KIND=dp), DIMENSION(:)              :: imass
    TYPE(cp_error_type), INTENT(inout)  &
                                             :: error
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

!!****** constraint_clv/rattle_con_eval [1.0] *
!!
!!   NAME
!!     rattle_con_eval
!!
!!   FUNCTION
!!     Evaluates the constraint for the rattle scheme
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     Teodoro Laino [teo] created 04.2006
!!
!!   SOURCE
!******************************************************************************
  FUNCTION rattle_con_eval(colvar, vel, roll, veps, rmat, particles, error) RESULT(res)
    TYPE(colvar_type), POINTER               :: colvar
    REAL(KIND=dp), INTENT(INOUT)             :: vel( :, : )
    LOGICAL, INTENT(IN), OPTIONAL            :: roll
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN), OPTIONAL                   :: veps, rmat
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particles
    TYPE(cp_error_type), INTENT(inout)  &
                                             :: error
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
          CALL MATVEC_3x3 (f_roll, rmat, colvar % dsdr(:,iatm) )
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
