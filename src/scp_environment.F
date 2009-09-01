!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief initialize scp environment
!> \author CJM
! *****************************************************************************
MODULE scp_environment
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE cp_dbcsr_operations,             ONLY: dbcsr_add_block_node,&
                                             dbcsr_allocate_matrix,&
                                             dbcsr_deallocate_matrix,&
                                             sm_from_dbcsr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE dbcsr_methods,                   ONLY: dbcsr_get_info
  USE dbcsr_operations,                ONLY: dbcsr_copy,&
                                             dbcsr_get_diag,&
                                             dbcsr_set
  USE dbcsr_types,                     ONLY: dbcsr_obj,&
                                             dbcsr_repl_full
  USE dbcsr_util,                      ONLY: dbcsr_finalize
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE distribution_2d_types,           ONLY: distribution_2d_type
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: get_particle_set,&
                                             particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE scp_coeff_types,                 ONLY: aux_coeff_set_create,&
                                             aux_coeff_set_initialize,&
                                             aux_coeff_set_release,&
                                             aux_coeff_set_type
  USE scp_energy_types,                ONLY: allocate_scp_energy,&
                                             scp_energy_type
  USE scp_environment_methods,         ONLY: scp_qs_env_setup
  USE scp_environment_types,           ONLY: scp_environment_type,&
                                             set_scp_env
  USE scp_force_types,                 ONLY: allocate_scp_force,&
                                             init_scp_force,&
                                             scp_force_type
  USE semi_empirical_types,            ONLY: get_se_param,&
                                             semi_empirical_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: scp_dft_init, scp_nddo_init

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'scp_environment'


CONTAINS

! *****************************************************************************
!> \brief Initializes the SCP environment using the QS environment
! *****************************************************************************
  SUBROUTINE scp_dft_init ( scp_env, qs_env, error )

    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(qs_environment_type), OPTIONAL, &
      POINTER                                :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scp_dft_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nkind, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: natom_of_kind
    LOGICAL                                  :: failure
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(aux_coeff_set_type), POINTER        :: aux_coeff_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(scp_energy_type), POINTER           :: energy
    TYPE(scp_force_type), POINTER            :: force( : )
    TYPE(section_vals_type), POINTER         :: dft_section, input

    CALL timeset(routineN,handle)

    failure = .FALSE.
    ! Nullifying  pointers
    NULLIFY ( aux_coeff_set )
    NULLIFY ( energy, force )
    NULLIFY ( para_env )
    NULLIFY ( input )
    NULLIFY ( atomic_kind_set )
    NULLIFY ( local_particles )

    IF ( PRESENT ( qs_env ) ) THEN
      CALL get_qs_env (qs_env, para_env = para_env, input = input, &
                       atomic_kind_set = atomic_kind_set, &
                       local_particles = local_particles, error = error )

      dft_section =>  section_vals_get_subs_vals( input, "DFT", error = error )
      ! Set the input section of the SCP
      CALL set_scp_env(scp_env,input=dft_section,&
           error=error)

      ! Allocate and initialize coefs pointer
      CALL aux_coeff_set_create(aux_coeff_set=aux_coeff_set,atomic_kind_set=atomic_kind_set,&
           distribution=local_particles, error=error)
      CALL aux_coeff_set_initialize(aux_coeff_set,para_env,error)
      CALL set_scp_env (scp_env=scp_env,aux_coeff_set = aux_coeff_set, error=error)
      CALL aux_coeff_set_release(aux_coeff_set, error=error)

      ! Allocate the data structure for energies ***
      CALL allocate_scp_energy (energy)

      ! Allocate the data structure for energies ***
      nkind = SIZE(atomic_kind_set)
      ALLOCATE (natom_of_kind(nkind),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             natom_of_kind=natom_of_kind)
      CALL allocate_scp_force(force,natom_of_kind)
      DEALLOCATE (natom_of_kind,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      CALL init_scp_force(force)

      ! set the scp_env
      CALL set_scp_env (scp_env=scp_env,energy=energy,force=force,error=error)

      !Sets up pw, grids, densities
      CALL scp_qs_env_setup(scp_env,qs_env,error)
    END IF
    CALL timestop(handle)

  END SUBROUTINE scp_dft_init
! *****************************************************************************
!> \brief Initializes the SCP environment using the QS environment
! *****************************************************************************
  SUBROUTINE scp_nddo_init ( scp_env, qs_env, error )

    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scp_nddo_init', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, iatom, ii, ikind, ilist, jj, natom, natorb, &
      neighbor_list_id, nkind, nrow, nsgf_a, nsgf_b, stat
    INTEGER, ALLOCATABLE                     :: first_sgf_a( : ), &
                                                first_sgf_b( : ), &
                                                last_sgf_a( : ), &
                                                last_sgf_b( : )
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: app, asp, ass
    REAL(KIND=dp), DIMENSION(9, 9)           :: a_scp
    REAL(KIND=dp), POINTER                   :: ascp_diag( : ), block( :, : )
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dbcsr_obj), POINTER                 :: ks_scp, pscp
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(scp_energy_type), POINTER           :: energy
    TYPE(semi_empirical_type), POINTER       :: se_kind

    CALL timeset(routineN,handle)

    failure = .FALSE.
    ! Nullifying  pointers
    NULLIFY ( energy )
    NULLIFY ( particle_set )
    NULLIFY ( para_env )
    NULLIFY ( atomic_kind_set, atomic_kind )
    NULLIFY ( distribution_2d )
    NULLIFY ( pscp, ks_scp, block, ascp_diag, se_kind, ascp_diag, atom_list )

    CALL get_qs_env (qs_env, para_env = para_env, &
                     atomic_kind_set = atomic_kind_set, &
                     particle_set = particle_set, &
                     distribution_2d = distribution_2d, &
                     neighbor_list_id=neighbor_list_id, &
                     error = error )

    ! Allocate the data structure for energies ***
    CALL allocate_scp_energy (energy)
    ! set the scp_env
    CALL set_scp_env (scp_env=scp_env,energy=energy,error=error)
    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,nsgf=nsgf_a)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,nsgf=nsgf_b)

    ALLOCATE (first_sgf_a(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (last_sgf_a(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (first_sgf_b(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (last_sgf_b(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL get_particle_set(particle_set=particle_set,&
                          first_sgf=first_sgf_a,&
                          last_sgf=last_sgf_a,&
                          error=error)
    CALL get_particle_set(particle_set=particle_set,&
                          first_sgf=first_sgf_b,&
                          last_sgf=last_sgf_b,&
                          error=error)

    CALL dbcsr_allocate_matrix(matrix=pscp,&
                               nblock_row=natom,&
                               nblock_col=natom,&
                               nrow=nsgf_a,&
                               ncol=nsgf_b,&
                               first_row=first_sgf_a,&
                               last_row=last_sgf_a,&
                               first_col=first_sgf_b,&
                               last_col=last_sgf_b,&
                               matrix_name="PSCP",&
                               matrix_symmetry="symmetric",&
                               sparsity_id=neighbor_list_id, &
                               distribution_2d=distribution_2d,&
                               replication_type=dbcsr_repl_full,&
                               error=error)

    DO ikind=1,nkind
       atomic_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,&
                            se_parameter=se_kind, &
                            atom_list=atom_list)
       CALL get_se_param(se_kind, ass=ass,app=app,asp=asp,natorb=natorb)
       IF (natorb>4) THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="SCP-NDDO not implemened for d-orbitals!", &
               error=error, error_level=cp_failure_level)
       END IF
       a_scp( 1  ,  1) = ass
       a_scp( 1  ,2:4) = asp
       a_scp( 2:4,  1) = asp
       a_scp( 2:4,2:4) = app

       ! Loop over the replicated list of particles, filling the sparse matrix structure to exploit
       DO ilist=1,SIZE ( atom_list )
          iatom = atom_list ( ilist )

          NULLIFY(block)
          CALL dbcsr_add_block_node(matrix=pscp,&
                                    block_row=iatom,&
                                    block_col=iatom,&
                                    BLOCK=BLOCK)
          DO ii=1,SIZE(block,1)
             DO jj=1,SIZE(block,2)
                block(ii,jj) = a_scp (ii,jj) 
             END DO
          END DO
       END DO ! ilist
    END DO ! ikind 

    CALL dbcsr_finalize(pscp,error=error)

    CALL dbcsr_get_info ( pscp, nfullrows_total=nrow )
    ALLOCATE (ascp_diag(nrow),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL dbcsr_get_diag( pscp, ascp_diag, error=error )
    CALL set_scp_env ( scp_env=scp_env, ascp_diag=ascp_diag, error = error )

    CALL dbcsr_set(pscp,0.0_dp,error=error)
    ALLOCATE(ks_scp)!sm->dbcsr
    CALL dbcsr_copy(ks_scp,pscp,'KS_SCP',error=error)

    CALL sm_from_dbcsr ( scp_env%pscp, pscp, distribution_2d, error=error )!sm->dbcsr
    CALL sm_from_dbcsr ( scp_env%ks_scp, ks_scp, distribution_2d, error=error )!sm->dbcsr
    CALL dbcsr_deallocate_matrix ( ks_scp, error )!sm->dbcsr
    CALL dbcsr_deallocate_matrix ( pscp, error )!sm->dbcsr

    DEALLOCATE (first_sgf_a,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (last_sgf_a,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (first_sgf_b,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (last_sgf_b,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE scp_nddo_init

END MODULE scp_environment

