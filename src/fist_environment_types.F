!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/fist_environment_types [1.0] *
!!
!!   NAME
!!     fist_environment_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM SEPT-01-02 
!!
!!   MODIFICATION HISTORY
!!   gt sept-23-02 added atomic_kind_set to replica_environment_type
!!   to allow use of kind_based neighbor list
!!
!!   SOURCE
!******************************************************************************

MODULE fist_environment_types
!  USE cp_para_types,                   ONLY: cp_para_env_type
!  USE ewald_parameters_types,          ONLY: ewald_parameters_type
!  USE kinds,                           ONLY: dbl
!  USE linklist_types,                  ONLY: linklist_internal_data_type
!  USE md,                              ONLY: virial_type
!  USE pair_potential_types,            ONLY: pair_potential_type
!  USE simulation_cell,                 ONLY: cell_type
!  PRIVATE
!
!  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='fist_environment_types'
!
!  TYPE pimd_environment_type
!     REAL ( dbl ), DIMENSION ( : , : ), POINTER :: u
!     REAL ( dbl ) :: lbd
!     REAL ( dbl ) :: slbd
!     REAL ( dbl ) :: meff
!     REAL ( dbl ) :: ismeff
!     REAL ( dbl ), DIMENSION ( 2 , 2 ) :: harmonic_prop
!  END TYPE pimd_environment_type
!
!  TYPE fist_thermodynamic_type
!     TYPE ( virial_type ), POINTER ::  virial
!     REAL ( dbl ) :: kin, pot, e_gspace, e_self, e_neut, e_bonded
!  END TYPE fist_thermodynamic_type
!
!  TYPE replica_environment_type
!    TYPE ( linklist_internal_data_type), POINTER :: ll_data ( : )
!    TYPE ( pimd_environment_type ), POINTER :: pimd_env
!  END TYPE replica_environment_type
!
!  TYPE fist_environment_type
!    !PRIVATE
!    INTEGER :: id_nr, ref_count 
!    TYPE ( cell_type ), POINTER                    :: cell
!    TYPE ( cell_type ), POINTER                    :: cell_ref
!    TYPE ( ewald_parameters_type ), POINTER        :: ewald_param
!    TYPE ( pair_potential_type ), POINTER          :: potparm ( :, : )
!    TYPE ( replica_environment_type ), POINTER     :: rep_env ( : )
!    TYPE ( fist_thermodynamic_type ), POINTER      :: thermo
!    TYPE ( cp_para_env_type ), POINTER             :: para_env
!    TYPE ( cp_subsystem_type ), POINTER            :: subsys
!  END TYPE fist_environment_type
!
!! *****************************************************************************
!
!!!****s* fist_environment_types/fist_environment_p_type [1.0] *
!!!
!!!   NAME
!!!     fist_environment_p_type
!!!
!!!   FUNCTION
!!!     to build arrays of pointers
!!!
!!!   ATTRIBUTES
!!!     - fist_env: the pointer to the fist_env
!!!
!!!   AUTHOR
!!!     CJM
!!!
!!!   MODIFICATION HISTORY
!!!     11/03
!!!
!!!   SOURCE
!  !***************************************************************************
!  TYPE fist_environment_p_type
!     TYPE(fist_environment_type), POINTER :: fist_env
!  END TYPE fist_environment_p_type
!!!***
!  !****************************************************************************
!
!! *** Public data types ***

!  PUBLIC :: fist_environment_type,        & 
!            fist_environment_p_type,      & 
!            replica_environment_type,     & 
!            fist_thermodynamic_type,      &
!            pimd_environment_type
!
!! *** Public subroutines ***
!
!  PUBLIC :: get_fist_env,&
!            init_fist_env,&
!            set_fist_env,&
!            replica_env_init
!
!  CHARACTER(LEN=*), PRIVATE, PARAMETER :: moduleN='fist_environment_types'
!  INTEGER, PRIVATE, SAVE :: last_fist_env_id_nr=0
!
!!!***
!! *****************************************************************************
!
!CONTAINS
!
!! *****************************************************************************
!
!  SUBROUTINE get_fist_env( fist_env, atomic_kind_set, molecule, particle_set, &
!                           local_particles, local_molecules, molecule_kind_set, &
!                           molecule_set, cell, cell_ref, rep_env, ewald_param,  &
!                           potparm, thermo, para_env, subsys, id_nr, error )
!    IMPLICIT NONE
!
!!   Purpose: Get the FIST environment.
!
!!   ***************************************************************************
!
!    TYPE ( fist_environment_type ), INTENT(IN)              :: fist_env
!    INTEGER, INTENT(out), OPTIONAL                          :: id_nr
!    TYPE ( cell_type ), OPTIONAL, POINTER                   :: cell, cell_ref
!    TYPE ( ewald_parameters_type ), OPTIONAL, POINTER       :: ewald_param
!    TYPE ( pair_potential_type ), OPTIONAL, POINTER         :: potparm ( :, : )
!    TYPE ( fist_thermodynamic_type ), OPTIONAL, POINTER     :: thermo
!    TYPE ( replica_environment_type ), OPTIONAL, POINTER    :: rep_env ( : )
!    TYPE(atomic_kind_type), OPTIONAL, POINTER               :: atomic_kind_set
!    TYPE(molecule_structure_type), OPTIONAL, POINTER        :: molecule ( : )
!    TYPE(particle_type), OPTIONAL, POINTER                  :: particle_set ( : )
!    TYPE(cp_para_env_type), OPTIONAL, POINTER               :: para_env
!    INTEGER, INTENT(out), OPTIONAL                          :: id_nr
!    TYPE(distribution_1d_type),OPTIONAL, POINTER            :: local_particles,&
!                                                               local_molecules
!    TYPE(cp_error_type), OPTIONAL,INTENT(INOUT)             :: error
!    TYPE(molecule_kind_type),  OPTIONAL, POINTER            :: molecule_kind_set ( : )
!    TYPE(molecule_type),  OPTIONAL, POINTER                 :: molecule_set ( : )
!    TYPE(cp_subsystem_type), POINTER, OPTIONAL              :: subsys
!
!    LOGICAL :: failure
!    CHARACTER(len=*), PARAMETER :: routineN='get_fist_env',&
!                                   routineP=moduleN//':'//routineN
!    TYPE (atomic_kind_list_type), POINTER :: atomic_kinds
!    TYPE (particle_list_type), POINTER :: particles
!    TYPE (mol_struct_list_type), POINTER :: local_molecules_old
!    TYPE (mol_new_list_type), POINTER :: molecules_new
!    TYPE (mol_kind_new_list_type), POINTER :: molecule_kinds_new
!
!!   ---------------------------------------------------------------------------
!    failure=.FALSE.
!    NULLIFY(atomic_kinds, particles,local_molecules_old, molecules_new,&
!         molecule_kinds_new)
!!    CPPrecondition(fist_env%ref_count>0,cp_failure_level,routineP,error,failure)
!    CPPrecondition(ASSOCIATED(fist_env%subsys),cp_failure_level,routineP,error,failure)
!
!    IF ( PRESENT ( cell ) ) cell => fist_env % cell
!    IF ( PRESENT ( cell_ref ) ) cell => fist_env % cell_ref
!    IF ( PRESENT ( ewald_param ) ) ewald_param => fist_env % ewald_param
!    IF ( PRESENT ( potparm ) ) potparm => fist_env % potparm
!    IF ( PRESENT ( thermo ) ) thermo => fist_env % thermo
!    IF ( PRESENT ( rep_env ) ) rep_env => fist_env % rep_env
!    IF ( PRESENT ( para_env ) ) para_env => fist_env % para_env
!    CALL cp_subsys_get(fist_env%subsys,&
!                      atomic_kinds=atomic_kinds,&
!                      local_molecules_new=local_molecules,&
!                      local_molecules=local_molecules_old,&
!                      local_particles=local_particles,&
!                      particles=particles,&
!                      molecule_kinds_new=molecule_kinds_new,&
!                      molecules_new=molecules_new,&
!                      error=error)
!    IF (PRESENT(atomic_kind_set)) atomic_kind_set => atomic_kinds%els
!    IF (PRESENT(molecule)) THEN
!      IF (ASSOCIATED(local_molecules_old)) THEN
!        molecule => local_molecules_old%els
!      END IF
!    END IF
!    IF (PRESENT(particle_set)) particle_set => particles%els
!    IF (PRESENT(molecule_kind_set)) molecule_kind_set => molecule_kinds_new%els
!    IF (PRESENT(molecule_set)) molecule_set => molecules_new%els
!
!  END SUBROUTINE get_fist_env
!
!! *****************************************************************************
!
!  SUBROUTINE init_fist_env ( fist_env, para_env, error )
!
!!   Purpose: Initialise the FIST environment.
!
!!   ***************************************************************************
!
!    TYPE(fist_environment_type), INTENT(OUT)     :: fist_env
!    TYPE(cp_para_env_type), POINTER              :: para_env
!    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error
!
!!   ---------------------------------------------------------------------------
!
!    NULLIFY (fist_env%cell)
!    NULLIFY (fist_env%cell_ref)
!    NULLIFY (fist_env%ewald_param)
!    NULLIFY (fist_env%potparm)
!    NULLIFY (fist_env%thermo)
!    NULLIFY (fist_env%rep_env)
!    CALL cp_subsys_create(fist_env%subsys, para_env=para_env, error=error)
!    CALL cp_para_env_retain(para_env,error=error)
!    fist_env%para_env => para_env
!    fist_env%ref_count=1
!    last_fist_env_id_nr=last_fist_env_id_nr+1
!    fist_env%id_nr=last_fist_env_id_nr
!
!  END SUBROUTINE init_fist_env
!
!! *****************************************************************************
!
!  SUBROUTINE set_fist_env( fist_env, atomic_kind_set, molecule, particle_set, &
!                           local_particles, local_molecules, molecule_kind_set, &
!                           molecule_set, cell, cell_ref, rep_env, ewald_param,  &
!                           potparm, thermo, subsys, error )
!
!    IMPLICIT NONE
!!   Purpose: Set the FIST environment.
!
!    TYPE ( fist_environment_type ), INTENT(out)             :: fist_env
!    TYPE ( cell_type ), OPTIONAL, POINTER                   :: cell
!    TYPE ( cell_type ), OPTIONAL, POINTER                   :: cell_ref
!    TYPE ( ewald_parameters_type ), OPTIONAL, POINTER       :: ewald_param
!    TYPE ( pair_potential_type ), OPTIONAL, POINTER         :: potparm ( :, : )
!    TYPE ( fist_thermodynamic_type ), OPTIONAL, POINTER     :: thermo
!    TYPE ( replica_environment_type ), OPTIONAL, POINTER    :: rep_env ( : )
!! locals
!    CHARACTER(len=*), PARAMETER :: routineN="set_fist_env",&
!         routineP=moduleN//':'//routineN
!    LOGICAL :: failure
!    TYPE (atomic_kind_list_type), POINTER :: atomic_kinds
!    TYPE (particle_list_type), POINTER :: particles
!    TYPE (mol_struct_list_type), POINTER :: local_molecules_old
!    TYPE (mol_new_list_type), POINTER :: molecules_new
!    TYPE (mol_kind_new_list_type), POINTER :: molecule_kinds_new
!
!    failure=.FALSE.
!    CPPrecondition(fist_env%ref_count>0,cp_failure_level,routineP,error,failure)
!  
!    IF (PRESENT(cell)) THEN
!       CALL cell_retain(cell, error=error)
!       CALL cell_release(fist_env%cell,error=error)
!       fist_env%cell => cell
!    END IF
!    IF (PRESENT(cell_ref)) THEN
!       CALL cell_retain(cell_ref, error=error)
!       CALL cell_release(fist_env%cell_ref,error=error)
!       fist_env%cell_ref => cell_ref
!    END IF
!    IF ( PRESENT ( ewald_param ) ) fist_env % ewald_param => ewald_param
!    IF ( PRESENT ( potparm ) ) fist_env % potparm => potparm
!    IF ( PRESENT ( thermo ) ) fist_env % thermo => thermo
!    IF ( PRESENT ( rep_env ) ) fist_env % rep_env => rep_env
!    IF (PRESENT(subsys)) THEN
!       CALL cp_subsys_retain(subsys,error=error)
!       CALL cp_subsys_release(fist_env%subsys,error=error)
!       fist_env%subsys => subsys
!    END IF
!    IF (PRESENT(atomic_kind_set)) THEN
!      CALL atomic_kind_list_create(atomic_kinds,&
!                                   els_ptr=atomic_kind_set,&
!                                   error=error)
!      CALL cp_subsys_set(fist_env%subsys,&
!                        atomic_kinds=atomic_kinds,&
!                        error=error)
!      CALL atomic_kind_list_release(atomic_kinds,error=error)
!    END IF
!    IF (PRESENT(molecule)) THEN
!      CALL mol_struct_list_create(local_molecules_old,&
!                                  els_ptr=molecule,&
!                                  error=error)
!      CALL cp_subsys_set(fist_env%subsys,&
!                        local_molecules=local_molecules_old,&
!                        error=error)
!      CALL mol_struct_list_release(local_molecules_old,error=error)
!    END IF
!    IF (PRESENT(particle_set)) THEN
!      CALL particle_list_create(particles,&
!                                els_ptr=particle_set,&
!                                error=error)
!      CALL cp_subsys_set(fist_env%subsys,&
!                        particles=particles,&
!                        error=error)
!      CALL particle_list_release(particles,error=error)
!    END IF
!    IF (PRESENT(local_particles)) THEN
!      CALL cp_subsys_set(fist_env%subsys,&
!                        local_particles=local_particles,&
!                        error=error)
!    END IF
!    IF (PRESENT(local_molecules)) THEN
!      CALL cp_subsys_set(fist_env%subsys,&
!                        local_molecules_new=local_molecules,&
!                        error=error)
!    END IF
!    IF (PRESENT(molecule_kind_set)) THEN
!      CALL mol_kind_new_list_create(molecule_kinds_new,&
!                                    els_ptr=molecule_kind_set,&
!                                    error=error)
!      CALL cp_subsys_set(fist_env%subsys,&
!                        molecule_kinds_new=molecule_kinds_new,&
!                        error=error)
!      CALL mol_kind_new_list_release(molecule_kinds_new,error=error)
!    END IF
!    IF (PRESENT(molecule_set)) THEN
!      CALL mol_new_list_create(molecules_new,&
!                               els_ptr=molecule_set,&
!                               error=error)
!      CALL cp_subsys_set(fist_env%subsys,&
!                        molecules_new=molecules_new,&
!                        error=error)
!      CALL mol_new_list_release(molecules_new,error=error)
!    END IF
!
!!***************************************************************************
!  END SUBROUTINE set_fist_env
!!***************************************************************************
!!
!!****f* fist_environment_type/replica_env_init [1.0] *
!!
!!   NAME
!!     replica_env_init
!!
!!   SYNOPSIS
!!     Subroutine replica_env_init(replica_env)
!!       Type(replica_environment_type), Intent (INOUT):: replica_env
!!     End Subroutine replica_env_init
!!
!!   FUNCTION
!!     initializes a replica environment
!!
!!   MODIFICATION HISTORY
!!     06.2003 created [fawzi]
!!
!!*** **********************************************************************
!SUBROUTINE replica_env_init(replica_env)
!  TYPE(replica_environment_type), INTENT(inout) :: replica_env
!  
!  NULLIFY( replica_env%ll_data, replica_env%pimd_env)
!
!END SUBROUTINE replica_env_init
! *****************************************************************************

!!****f* fist_environment_types/fist_env_create [1.0] *
!!
!!   NAME
!!     fist_env_create
!!
!!   SYNOPSIS
!!     Subroutine fist_env_create(fist_env, para_env, error)
!!       Type(fist_environment_type), Pointer:: fist_env
!!       Type(cp_para_env_type), Pointer:: para_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine fist_env_create
!!
!!   FUNCTION
!!     allocates and intitializes a fist_env
!!
!!   ARGUMENTS
!!     - fist_env: the object to create
!!     - para_env: the parallel environement for the qs_env
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     12.2002 created [fawzi]
!!
!!*** **********************************************************************
!SUBROUTINE fist_env_create(fist_env,para_env,error)
!  TYPE(fist_environment_type), POINTER :: fist_env
!  TYPE(cp_para_env_type), POINTER :: para_env
!  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
!  
!  LOGICAL :: failure
!  CHARACTER(len=*), PARAMETER :: routineN='fist_env_create',&
!        routineP=moduleN//':'//routineN
!  INTEGER :: stat
!
!  failure=.FALSE.
!  
!  ALLOCATE(fist_env, stat=stat)
!  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!  IF (.NOT. failure) THEN
!     CALL init_fist_env(fist_env,para_env=para_env, error=error)
!  END IF
!END SUBROUTINE fist_env_create
!!***************************************************************************

!!****f* fist_environment_types/fist_env_retain [1.0] *
!!
!!   NAME
!!     fist_env_retain
!!
!!   SYNOPSIS
!!     Subroutine fist_env_retain(fist_env, error)
!!       Type(fist_environment_type), Pointer:: fist_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine fist_env_retain
!!
!!   FUNCTION
!!     retains the given fist_env (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - fist_env: the object to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     12.2002 created [fawzi]
!!
!!*** **********************************************************************
!SUBROUTINE fist_env_retain(fist_env,error)
!  TYPE(fist_environment_type), POINTER :: fist_env
!  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
!  LOGICAL :: failure
!  CHARACTER(len=*), PARAMETER :: routineN='fist_env_retain',&
!        routineP=moduleN//':'//routineN
!  failure=.FALSE.
!!!  
!  CPPrecondition(ASSOCIATED(fist_env),cp_failure_level,routineP,error,failure)
!  IF (.NOT. failure) THEN
!     CPPrecondition(fist_env%ref_count>0,cp_failure_level,routineP,error,failure)
!     fist_env%ref_count=fist_env%ref_count+1
!  END IF
!END SUBROUTINE fist_env_retain
!***************************************************************************

!!****f* fist_environment_types/fist_env_release [1.0] *
!!
!!   NAME
!!     fist_env_release
!!
!!   SYNOPSIS
!!     Subroutine fist_env_release(fist_env, error)
!!       Type(fist_environment_type), Pointer:: fist_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine fist_env_release
!!
!!   FUNCTION
!!     releases the given fist_env (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - fist_env: the object to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     12.2002 created [fawzi]
!!
!!*** **********************************************************************
!SUBROUTINE fist_env_release(fist_env,error)
!  TYPE(fist_environment_type), POINTER :: fist_env
!  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
!
!  LOGICAL :: failure
!  CHARACTER(len=*), PARAMETER :: routineN='fist_env_release',&
!       routineP=moduleN//':'//routineN
!  INTEGER :: i, stat
!  failure=.FALSE.
!
!  IF (ASSOCIATED(fist_env)) THEN
!     CPPrecondition(fist_env%ref_count>0,cp_failure_level,routineP,error,failure)
!     fist_env%ref_count=fist_env%ref_count-1
!     IF (qs_env%ref_count<1) THEN
!        CALL cell_release(fist_env%cell)
!        CALL cell_release(fist_env%cell_ref)
!        CALL cp_para_env_release(fist_env%para_env,error=error)
!        DEALLOCATE(fist_env, stat=stat)
!        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
!     END IF
!  END IF
!  NULLIFY(fist_env)
!END SUBROUTINE fist_env_release
!***************************************************************************
END MODULE fist_environment_types
