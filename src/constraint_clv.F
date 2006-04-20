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
  USE cell_types,                      ONLY: cell_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE colvar_types,                    ONLY: dist_colvar_id,&
                                             colvar_type,&
                                             colvar_counters
  USE colvar_methods,                  ONLY: colvar_eval_mol_f
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_p_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE kinds,                           ONLY: dp
  USE mathlib,                         ONLY: invert_matrix
  USE md,                              ONLY: simulation_parameters_type
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE molecule_kind_types,             ONLY: colvar_constraint_type,&
                                             fixd_constraint_type,&
                                             get_molecule_kind,&
                                             molecule_kind_type
  USE molecule_types_new,              ONLY: get_molecule,&
                                             local_colvar_constraint_type,&
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
  SUBROUTINE shake_colv( molecule, particle_set, pos, vel, simpar, cell, ncolv,&
       imass, error )

    TYPE(molecule_type), POINTER             :: molecule
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: pos( :, : ), vel( :, : )
    TYPE(cell_type), POINTER                 :: cell
    TYPE(colvar_counters), INTENT(IN)        :: ncolv
    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    REAL(KIND=dp),DIMENSION(:)               :: imass

    CHARACTER(len=*), PARAMETER :: routineN = 'shake_colv', &
         routineP = moduleN//':'//routineN
    INTEGER                                  :: first_atom, iconst, index_a, &
                                                index_b
    REAL(KIND=dp)                            :: del_lam, dtby2, dtsqby2, dtsqi,&
                                                fdotf, fdotf_sum, imass1, &
                                                imass2, mass
    REAL(KIND=dp), DIMENSION(3)              :: fc_t, r12
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(colvar_constraint_type), POINTER    :: colv_list( : )
    TYPE(local_colvar_constraint_type), &
      POINTER                                :: lcolv( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    REAL(KIND=dp)                            :: max_sigma
    INTEGER                                  :: ishake
    LOGICAL                                  :: failure
    INTEGER                                  :: stat, jconst

    failure = .FALSE.
    ishake = 0       
    max_sigma = -1.0E+10_dp
    dtsqi   = 1.0_dp/(simpar%dt*simpar%dt)
    dtsqby2 = simpar%dt*simpar%dt*.5_dp
    dtby2   = simpar%dt*.5_dp
    CONS: DO WHILE (ABS(max_sigma)>=simpar%shake_tol)
       max_sigma = 0.0_dp
       ishake  = ishake + 1
       molecule_kind => molecule % molecule_kind
       CALL get_molecule_kind ( molecule_kind, colv_list = colv_list )
       CALL get_molecule ( molecule, first_atom = first_atom, lcolv=lcolv )

       IF (ishake==1) THEN
          DO iconst = 1, SIZE(colv_list)             
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
             ! Update colvar
             CALL colvar_eval_mol_f( lcolv ( iconst ) % colvar, cell, particles=particle_set,&
                  pos=pos, error=error)
             lcolv ( iconst ) % sigma = lcolv ( iconst ) % colvar % ss - &
                  colv_list ( iconst ) % expected_value
             fdotf_sum = eval_Jac_colvar(lcolv ( iconst ) % colvar,&
                  lcolv ( iconst ) % colvar_old, particle_set, imass=imass, error=error)
             del_lam = 2.0_dp*lcolv ( iconst ) % sigma/(simpar%dt*simpar%dt*fdotf_sum)
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
          CALL colvar_eval_mol_f( lcolv ( iconst ) % colvar, cell, particles=particle_set,&
               pos=pos, error=error)
          lcolv ( iconst ) % sigma = lcolv ( iconst ) % colvar % ss - &
               colv_list ( iconst ) % expected_value
          max_sigma = MAX(ABS(lcolv ( iconst ) % sigma),max_sigma)
       END DO
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
  SUBROUTINE rattle_colv ( molecule, particle_set, vel, simpar, ncolv, cell, imass, error )

    TYPE(molecule_type), POINTER             :: molecule
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: vel( :, : )
    TYPE(colvar_counters), INTENT(IN)        :: ncolv
    TYPE(cell_type), POINTER                 :: cell
    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar
    REAL(KIND=dp), DIMENSION(:)              :: imass
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error 

    INTEGER                                  :: first_atom, iconst, index_a, &
                                                index_b
    REAL(KIND=dp)                            :: del_lam, dtby2, fdotf, &
                                                fdotf_sum, imass1, imass2, &
                                                mass
    REAL(KIND=dp), DIMENSION(3)              :: r12, v12
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(colvar_constraint_type), POINTER    :: colv_list( : )
    TYPE(local_colvar_constraint_type), &
      POINTER                                :: lcolv( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    REAL(KIND=dp)                            :: max_sigma
    INTEGER                                  :: irattle
    
    irattle   = 0
    max_sigma = -1.0E+10_dp
    dtby2     = simpar%dt*.5_dp
    CONS: DO WHILE ( ABS ( max_sigma ) >= simpar % shake_tol )
       max_sigma = 0.0_dp
       irattle = irattle + 1
       molecule_kind => molecule % molecule_kind
       CALL get_molecule_kind ( molecule_kind, colv_list = colv_list )
       CALL get_molecule ( molecule, first_atom = first_atom, lcolv=lcolv )
       IF (irattle==1) THEN
          DO iconst = 1, SIZE(colv_list)
             ! Update colvar_old
             CALL colvar_eval_mol_f(lcolv ( iconst ) % colvar_old, cell,&
                  particles=particle_set, error=error)                
             ! Update velocities 
             CALL update_con_colv(vel, dtby2, lcolv(iconst), &
                  particles=particle_set, lambda=lcolv(iconst)%lambda,&
                  imass=imass,error=error)                
          END DO
       ELSE             
          DO iconst = 1, SIZE(colv_list)
             lcolv ( iconst ) % sigma = rattle_con_eval(lcolv ( iconst ) % colvar_old, vel)
             fdotf_sum = eval_Jac_colvar(lcolv ( iconst ) % colvar_old,&
                  lcolv ( iconst ) % colvar_old, particle_set,imass=imass,error=error)
             del_lam = 2.0_dp*lcolv ( iconst ) % sigma/(simpar%dt*fdotf_sum)
             lcolv ( iconst ) % lambda = lcolv ( iconst ) % lambda + del_lam
             
             ! Update velocities 
             CALL update_con_colv(vel, dtby2, lcolv(iconst), &
                  particles=particle_set, lambda=del_lam,&
                  imass=imass,error=error)  
          END DO
       END IF
       
       DO iconst = 1, SIZE(colv_list)
          lcolv ( iconst ) % sigma = rattle_con_eval(lcolv ( iconst ) % colvar_old, vel)
          max_sigma = MAX(ABS(lcolv ( iconst ) % sigma),max_sigma)
       END DO
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
       simpar, ncolv, cell, imass, error )

    TYPE(molecule_type), POINTER             :: molecule
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: pos( :, : ), vel( :, : )
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: r_shake, v_shake
    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar
    REAL(KIND=dp), DIMENSION(:)    :: imass
    TYPE(cell_type), POINTER                 :: cell
    TYPE(colvar_counters), INTENT(IN)        :: ncolv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error    

    INTEGER                                  :: first_atom, iconst, index_a, &
                                                index_b,k
    REAL(KIND=dp)                            :: del_lam, dtby2, dtsqby2, &
                                                fdotf, fdotf_sum, imass1, &
                                                imass2, mass
    REAL(KIND=dp), DIMENSION(3)              :: f_roll_r, f_roll_v, fc_t, r12
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(colvar_constraint_type), POINTER    :: colv_list( : )
    TYPE(local_colvar_constraint_type), &
      POINTER                                :: lcolv( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    REAL(KIND=dp)                            :: max_sigma
    INTEGER                                  :: ishake

    dtsqby2 = simpar%dt*simpar%dt*.5_dp
    dtby2 = simpar%dt*.5_dp
    max_sigma = -1.0E+10_dp
    ishake = 0
    CONS: DO WHILE (ABS(max_sigma)>=simpar%shake_tol)
       max_sigma = 0.0_dp
       ishake = ishake + 1
       molecule_kind => molecule % molecule_kind
       CALL get_molecule_kind ( molecule_kind, colv_list = colv_list )
       CALL get_molecule ( molecule, first_atom = first_atom, lcolv=lcolv )
       IF (ishake==1) THEN
          DO iconst = 1, SIZE(colv_list)
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
             ! Update colvar
             CALL colvar_eval_mol_f( lcolv ( iconst ) % colvar, cell, particles=particle_set,&
                  pos=pos, error=error)
             lcolv ( iconst ) % sigma = lcolv ( iconst ) % colvar % ss - &
                  colv_list ( iconst ) % expected_value
             fdotf_sum = eval_Jac_colvar(lcolv ( iconst ) % colvar,&
                  lcolv ( iconst ) % colvar_old, particle_set, roll=.TRUE., rmat=r_shake,&
                  imass=imass, error=error)
             del_lam = 2.0_dp*lcolv ( iconst ) % sigma/(simpar%dt*simpar%dt*fdotf_sum)
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
          CALL colvar_eval_mol_f( lcolv ( iconst ) % colvar, cell, particles=particle_set,&
               pos=pos, error=error)
          lcolv ( iconst ) % sigma = lcolv ( iconst ) % colvar % ss - &
               colv_list ( iconst ) % expected_value
          max_sigma = MAX(ABS(lcolv ( iconst ) % sigma),max_sigma)
       END DO
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
       simpar, veps, ncolv, cell, imass, error )

    TYPE(molecule_type), POINTER             :: molecule
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: vel( :, : )
    REAL(KIND=dp), INTENT(IN)                :: r_rattle( :, : )
    TYPE(cell_type), POINTER                 :: cell
    TYPE(colvar_counters), INTENT(IN)        :: ncolv
    REAL(KIND=dp), DIMENSION(:)    :: imass
    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar
    REAL(KIND=dp), INTENT(IN)                :: veps( :, : )
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    INTEGER                                  :: first_atom, iconst, index_a, &
                                                index_b,k
    REAL(KIND=dp)                            :: del_lam, dtby2, fdotf, &
                                                fdotf_sum, imass1, imass2, &
                                                mass
    REAL(KIND=dp), DIMENSION(3)              :: f_roll, r1, r12, r2, v12
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(colvar_constraint_type), POINTER    :: colv_list( : )
    TYPE(local_colvar_constraint_type), &
      POINTER                                :: lcolv( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    REAL(KIND=dp)                            :: max_sigma
    INTEGER                                  :: irattle

    dtby2 = simpar%dt*.5_dp
    max_sigma = -1.0E+10_dp
    irattle = 0
    CONS: DO WHILE (ABS(max_sigma)>=simpar%shake_tol)
       max_sigma = 0.0_dp
       irattle = irattle + 1
       molecule_kind => molecule % molecule_kind
       CALL get_molecule_kind ( molecule_kind, colv_list = colv_list )
       CALL get_molecule ( molecule, first_atom = first_atom, lcolv=lcolv )
       IF (irattle==1) THEN
          DO iconst = 1, SIZE ( colv_list )
             ! Update colvar_old
             CALL colvar_eval_mol_f(lcolv ( iconst ) % colvar_old, cell,&
                  particles=particle_set, error=error)                
             ! Update velocities 
             CALL update_con_colv(vel, dtby2, lcolv(iconst), &
                  particles=particle_set, lambda=lcolv(iconst)%lambda,&
                  imass=imass,error=error)
          END DO
       ELSE
          DO iconst = 1, SIZE ( colv_list )
             lcolv ( iconst ) % sigma = rattle_con_eval(lcolv ( iconst ) % colvar_old, vel,&
                  roll=.TRUE., veps=veps, rmat=r_rattle, particles=particle_set,&
                  error=error)
             fdotf_sum = eval_Jac_colvar(lcolv ( iconst ) % colvar_old,&
                  lcolv ( iconst ) % colvar_old, particle_set,roll=.TRUE.,&
                  rmat=r_rattle,imass=imass, error=error)
             del_lam = 2.0_dp*lcolv ( iconst ) % sigma/(simpar%dt*fdotf_sum)
             lcolv ( iconst ) % lambda = lcolv ( iconst ) % lambda + del_lam
             
             ! Update velocities 
             CALL update_con_colv(vel, dtby2, lcolv(iconst), &
                  particles=particle_set, lambda=del_lam,&
                  roll=.TRUE.,rmat=r_rattle,imass=imass,error=error)  
          END DO
       END IF
       ! computing the constraint and value of the tolerance
       DO iconst = 1, SIZE(colv_list)
          lcolv ( iconst ) % sigma = rattle_con_eval(lcolv ( iconst ) % colvar_old, vel,&
               roll=.TRUE., veps=veps, rmat=r_rattle, particles=particle_set,&
               error=error)
          max_sigma = MAX(ABS(lcolv ( iconst ) % sigma),max_sigma)
       END DO
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
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(INOUT)             :: wrk( :, : )
    REAL(KIND=dp), INTENT(IN)                :: fac, lambda
    REAL(KIND=dp), DIMENSION(:)              :: imass
    TYPE(local_colvar_constraint_type),&
         INTENT(IN)                          :: lcolv
    TYPE(particle_type), POINTER, DIMENSION(:) :: particles
    LOGICAL, INTENT(in), OPTIONAL            :: roll
    REAL(KIND=dp), DIMENSION(:, :), OPTIONAL, &
      INTENT(IN)                             :: rmat
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'update_con_colv', &
         routineP = moduleN//':'//routineN
    LOGICAL :: failure, my_roll
    REAL(KIND=dp)  :: del_lam, mass
    INTEGER :: stat, iatm, ind
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
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
             f_roll = MATVEC_3X3 ( rmat, lcolv % colvar_old % dsdr(:,iatm)  )
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
    IMPLICIT NONE
    TYPE(colvar_type), POINTER :: colvar, colvar_old
    TYPE(particle_type), POINTER, DIMENSION(:) :: particles
    REAL(KIND=dp)              :: res
    REAL(KIND=dp), DIMENSION(:) :: imass
    LOGICAL, INTENT(IN), OPTIONAL :: roll
    REAL(KIND=dp), DIMENSION(:, :), OPTIONAL, &
      INTENT(IN)                             :: rmat    
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'eval_Jac_colvar', &
         routineP = moduleN//':'//routineN
    REAL(KIND=dp), DIMENSION(:), POINTER  :: wrk1, wrk2
    REAL(KIND=dp), DIMENSION(3) :: tmp1, tmp2
    INTEGER :: stat, i, ind, iatom
    LOGICAL :: failure, my_roll
    REAL(KIND=dp)                            :: mass
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

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
       tmp2 = colvar_old%dsdr(1:3,i)
       IF (my_roll) THEN          
          tmp2 =  MATVEC_3X3 ( rmat, tmp2 )
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
    IMPLICIT NONE
    TYPE(colvar_type), POINTER               :: colvar
    REAL(KIND=dp), INTENT(INOUT)             :: vel( :, : )
    LOGICAL, INTENT(IN), OPTIONAL            :: roll
    REAL(KIND=dp), INTENT(IN), DIMENSION(:,:),&
         OPTIONAL               :: rmat, veps
    TYPE(particle_type), POINTER, DIMENSION(:),&
         OPTIONAL               :: particles
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    REAL(KIND=dp)               :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'rattle_con_eval', &
         routineP = moduleN//':'//routineN
    INTEGER :: iatm, ind
    LOGICAL :: my_roll, failure
    REAL(KIND=dp), DIMENSION(3) :: v_roll, f_roll, pos

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
          f_roll = MATVEC_3x3 ( rmat, colvar % dsdr(:,iatm) )
          v_roll = vel(:,ind) + MATVEC_3X3 ( veps, pos )
       ELSE          
          f_roll = colvar % dsdr(:,iatm)
          v_roll = vel(:,ind)
       END IF
       res = res + DOT_PRODUCT(f_roll, v_roll)
    END DO

  END FUNCTION rattle_con_eval

END MODULE constraint_clv
