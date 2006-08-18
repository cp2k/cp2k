!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


!!****h* cp2k/ep_types [1.0] *
!!
!!   NAME
!!     ep_types
!!
!!   FUNCTION
!!     the types needed for the calculation of the enrgy with perturbation
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     11.2002 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE ep_types
  USE cp_fm_types,                     ONLY: cp_fm_p_type
  USE cp_fm_vect,                      ONLY: cp_fm_vect_dealloc
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE global_types,                    ONLY: global_environment_type,&
                                             globenv_release,&
                                             globenv_retain
  USE input_section_types,             ONLY: section_vals_release,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE preconditioner,                  ONLY: destroy_preconditioner,&
                                             preconditioner_p_type
  USE qs_environment_types,            ONLY: qs_env_release,&
                                             qs_environment_type
  USE qs_p_env_types,                  ONLY: p_env_release,&
                                             qs_p_env_type
  USE qs_p_sparse_psi,                 ONLY: p_proj_release,&
                                             qs_p_projection_p_type
  USE replica_types,                   ONLY: rep_env_release,&
                                             replica_env_type
  USE sparse_matrix_types,             ONLY: deallocate_matrix,&
                                             real_matrix_p_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ep_types'
  INTEGER, SAVE :: last_ep_env_id=0, last_force_id=0
  LOGICAL, SAVE, PRIVATE :: module_initialized=.FALSE.

  PUBLIC :: ep_env_type, ep_env_p_type, ep_energy_type, ep_force_type
  PUBLIC :: ep_env_retain, ep_env_release, ep_env_get
  PUBLIC :: ep_envs_get_ep_env, ep_env_create, ep_env_calc_e_f
  PUBLIC :: ep_energy_zero, ep_dealloc_local_sm, ep_force_create,&
       ep_force_retain, ep_force_release, ep_force_zero
  !!***
  !****************************************************************************


  !!****s* ep_types/ep_energy_type *
  !!
  !!   NAME
  !!     ep_energy_type
  !!
  !!   FUNCTION
  !!     stores various energies calculated during an ep calculation
  !!
  !!   NOTES
  !!     -
  !!
  !!   ATTRIBUTES
  !!     - e_no_int: energy of the single fragments (no interaction)
  !!     - e0: the energy of the orthogonalized frozen densities
  !!       (energy of psi0, E(psi0))
  !!     - e1: the correction energy (E(psi0+psi1)-E(psi0))
  !!     - e_tot: total energy
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** **********************************************************************
  TYPE ep_energy_type
     REAL(dp) :: e_no_int,e0,e1,e_tot
  END TYPE ep_energy_type
  !***************************************************************************

  !!****s* ep_types/ep_force_type *
  !!
  !!   NAME
  !!     ep_force_type
  !!
  !!   FUNCTION
  !!     stores various forces calculated during an ep calculation
  !!     (only those not stored in main_qs_env)
  !!
  !!   NOTES
  !!     -
  !!
  !!   ATTRIBUTES
  !!     - f0_internal: intermolecular forces due to e0
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** **********************************************************************
  TYPE ep_force_type
     INTEGER :: id_nr, ref_count
     REAL(dp), DIMENSION(:,:), POINTER :: f0_internal
  END TYPE ep_force_type
  !***************************************************************************

  !!****s* ep_types/ep_env_type [1.0] *
  !!
  !!   NAME
  !!     ep_env_type
  !!
  !!   FUNCTION
  !!     environement for an ep calculation
  !!
  !!   NOTES
  !!     -
  !!
  !!   ATTRIBUTES
  !!     - id_nr: identification number (unique)
  !!     - ref_count: reference count (see doc/ReferenceCounting.html)
  !!     - mol_envs: qs_envs that represent the pieces of the full system
  !!       handled by a replica_env
  !!     - sub_proj: projection from the full system to the pieces
  !!     - main_qs_env: qs_env of the full system
  !!     - main_p_env: perturbation environment for the full system
  !!     - sub_p_env: perturbation environement for the molecular subsystem
  !!     - m_pi_Hrho_psi0d: the rhs of the linear system to solve for
  !!       the perturbation ( - pi Hrho psi0d )
  !!     - nspins: number of spins
  !!     - nmol: number of molecules replicas in mol_envs
  !!     - nat: number of atoms (in main_qs_env)
  !!     - nat_per_mol: number of atoms in each mol
  !!     - sub_nao: number of atom orbitals in the mol_env
  !!     - sub_nmo: number of molecular orbitals in the mol_env
  !!     - full_nao: number of atom orbitals in main_qs_env
  !!     - full_nmo: number of molecular orbitals in main_qs_env
  !!     - at2sub: which molecule owns each atom
  !!     - para_env: the global parallel environment for the ep calculation
  !!     - input: the input for the ep calculation
  !!     - globenv: global environment for the main_qs
  !!     - energy: the energy of the last configuration evaluated
  !!     - local_H_full: diagonal elements of the KohnSham matrix
  !!     - local_matrix_w: matrix_w for the intra forces
  !!
  !!   AUTHOR
  !!     Fawzi Mohamed
  !!
  !!   MODIFICATION HISTORY
  !!     11.2002 created [fawzi]
  !!
  !!   SOURCE
  !***************************************************************************
  TYPE ep_env_type
     INTEGER :: id_nr, ref_count, f_env_id,nspins,nat,nmol,nat_per_mol
     TYPE(global_environment_type), POINTER :: globenv
     TYPE(cp_para_env_type), POINTER :: para_env
     TYPE(section_vals_type), POINTER :: input
     TYPE(replica_env_type), POINTER :: mol_envs
     TYPE(qs_p_projection_p_type), &
          DIMENSION(:), POINTER        :: sub_proj
     TYPE(qs_environment_type), POINTER :: main_qs_env
     TYPE(qs_p_env_type), POINTER :: main_p_env, sub_p_env
     TYPE(cp_fm_p_type), DIMENSION(:), POINTER :: m_pi_Hrho_psi0d
     TYPE(cp_fm_p_type), DIMENSION(:), POINTER :: psi1
     TYPE(preconditioner_p_type), DIMENSION(:), POINTER :: precond
     INTEGER, DIMENSION(:), POINTER :: sub_nmo,sub_nao,full_nmo,full_nao,at2sub
     TYPE(ep_energy_type) :: energy
     TYPE(ep_force_type), POINTER :: force
  END TYPE ep_env_type
  !!***
  !****************************************************************************

  !!****s* ep_types/ep_env_p_type [1.0] *
  !!
  !!   NAME
  !!     ep_env_p_type
  !!
  !!   FUNCTION
  !!     type to have array of pointers
  !!
  !!   NOTES
  !!     -
  !!
  !!   ATTRIBUTES
  !!     - ep_env: the pointer to the ep env
  !!
  !!   AUTHOR
  !!     Fawzi Mohamed
  !!
  !!   MODIFICATION HISTORY
  !!     12.2002 created [fawzi]
  !!
  !!   SOURCE
  !***************************************************************************
  TYPE ep_env_p_type
     TYPE(ep_env_type), POINTER :: ep_env
  END TYPE ep_env_p_type
  !!***
  !****************************************************************************

  TYPE(ep_env_p_type), DIMENSION(:), POINTER, PRIVATE :: ep_envs

CONTAINS


  !!****f* ep_types/ep_energy_zero *
  !!
  !!   NAME
  !!     ep_energy_zero
  !!
  !!   FUNCTION
  !!     zeroes out the energies
  !!
  !!   NOTES
  !!     -
  !!
  !!   ARGUMENTS
  !!     - energy: the energy object to zero out
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** **********************************************************************
  SUBROUTINE ep_energy_zero(ep_energy,error)
    TYPE(ep_energy_type), INTENT(out)        :: ep_energy
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_energy_zero', &
      routineP = moduleN//':'//routineN

    ep_energy%e_no_int=0._dp
    ep_energy%e0=0._dp
    ep_energy%e1=0._dp
    ep_energy%e_tot=0._dp
  END SUBROUTINE ep_energy_zero
  !***************************************************************************

  !!****f* ep_types/ep_env_retain [1.0] *
  !!
  !!   NAME
  !!     ep_env_retain
  !!
  !!   FUNCTION
  !!     retains an ep environment (see doc/ReferenceCounting.html)
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     - ep_env: the environment to retain
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     Fawzi Mohamed
  !!
  !!   MODIFICATION HISTORY
  !!     11.2002 created [fawzi]
  !!
  !!*** **********************************************************************
  SUBROUTINE ep_env_retain(ep_env, error)
    TYPE(ep_env_type), POINTER               :: ep_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(ep_env),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPrecondition(ep_env%ref_count>0,cp_failure_level,routineP,error,failure)
       ep_env%ref_count=ep_env%ref_count+1
    END IF
  END SUBROUTINE ep_env_retain
  !***************************************************************************

  !!****f* ep_types/ep_env_release [1.0] *
  !!
  !!   NAME
  !!     ep_env_release
  !!
  !!   FUNCTION
  !!     releases the given ep_env (see doc/ReferenceCounting.html)
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     - ep_env: the ep_env to release
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     Fawzi Mohamed
  !!
  !!   MODIFICATION HISTORY
  !!     11.2002 created [fawzi]
  !!
  !!*** **********************************************************************
  SUBROUTINE ep_env_release(ep_env, error)
    TYPE(ep_env_type), POINTER               :: ep_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ierr, stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(ep_env)) THEN
       CPPrecondition(ep_env%ref_count>0,cp_failure_level,routineP,error,failure)
       ep_env%ref_count=ep_env%ref_count-1
       IF (ep_env%ref_count<1) THEN
          CALL rep_env_release(ep_env%mol_envs)
          IF (ASSOCIATED(ep_env%sub_proj)) THEN
             DO i=1,SIZE(ep_env%sub_proj)
                CALL p_proj_release(ep_env%sub_proj(i)%projection,error=error)
             END DO
             DEALLOCATE(ep_env%sub_proj,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          CALL qs_env_release(ep_env%main_qs_env)
          CALL p_env_release(ep_env%main_p_env,error=error)
          CALL p_env_release(ep_env%sub_p_env,error=error)
          CALL cp_fm_vect_dealloc(ep_env%m_pi_Hrho_psi0d,error=error)
          CALL cp_fm_vect_dealloc(ep_env%psi1,error=error)
          CALL cp_destroy_fenv(ep_env%f_env_id,ierr)
          CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
          IF (ASSOCIATED(ep_env%precond)) THEN
             DO i=1,1 !SIZE(ep_env%precond) ! same precond for all spins
                CALL destroy_preconditioner(ep_env%precond(i)%preconditioner)
             END DO
             DEALLOCATE(ep_env%precond,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF (ASSOCIATED(ep_env%sub_nmo)) THEN
             DEALLOCATE(ep_env%sub_nmo,ep_env%sub_nao,ep_env%full_nmo,ep_env%full_nao,&
                  stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(ep_env%at2sub)) THEN
             DEALLOCATE(ep_env%at2sub,stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
          CALL cp_para_env_release(ep_env%para_env,error=error)
          CALL section_vals_release(ep_env%input,error=error)
          CALL globenv_release(ep_env%globenv,error=error)
          CALL ep_envs_rm_ep_env(ep_env,error=error)
          CALL ep_force_release(ep_env%force,error=error)
          DEALLOCATE(ep_env, stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF
    NULLIFY(ep_env)
  END SUBROUTINE ep_env_release
  !***************************************************************************

  !!****f* ep_methods/ep_env_get [1.0] *
  !!
  !!   NAME
  !!     ep_env_get
  !!
  !!   FUNCTION
  !!     returns the various attributes of the ep_env
  !!
  !!   NOTES
  !!     -
  !!
  !!   ARGUMENTS
  !!     - ep_env: the ep_env to get the info from
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling
  !!
  !!     For the other arguments see the attributes of ep_env_type
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!   MODIFICATION HISTORY
  !!     12.2002 created [fawzi]
  !!
  !!*** **********************************************************************
  SUBROUTINE ep_env_get(ep_env,id_nr, mol_envs, sub_proj, main_qs_env,&
       main_p_env,sub_p_env, preconditioner, m_pi_Hrho_psi0d, psi1, f_env_id,&
       globenv,at2sub,error)
    TYPE(ep_env_type), POINTER               :: ep_env
    INTEGER, INTENT(out), OPTIONAL           :: id_nr
    TYPE(replica_env_type), OPTIONAL, &
      POINTER                                :: mol_envs
    TYPE(qs_p_projection_p_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: sub_proj
    TYPE(qs_environment_type), OPTIONAL, &
      POINTER                                :: main_qs_env
    TYPE(qs_p_env_type), OPTIONAL, POINTER   :: main_p_env, sub_p_env
    TYPE(preconditioner_p_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: preconditioner
    TYPE(cp_fm_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: m_pi_Hrho_psi0d, psi1
    INTEGER, INTENT(out), OPTIONAL           :: f_env_id
    TYPE(global_environment_type), &
      OPTIONAL, POINTER                      :: globenv
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: at2sub
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_env_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(ep_env),cp_failure_level,routineP,error,failure)
    IF (.not.failure) THEN
       CPPrecondition(ep_env%ref_count>0,cp_failure_level,routineP,error,failure)
    END IF

    IF (.NOT. failure) THEN
       IF (PRESENT(id_nr)) id_nr=ep_env%id_nr
       IF (PRESENT(mol_envs)) mol_envs => ep_env%mol_envs
       IF (PRESENT(sub_proj)) sub_proj => ep_env%sub_proj
       IF (PRESENT(main_qs_env)) main_qs_env => ep_env%main_qs_env
       IF (PRESENT(main_p_env)) main_p_env => ep_env%main_p_env
       IF (PRESENT(sub_p_env)) sub_p_env => ep_env%sub_p_env
       IF (PRESENT(m_pi_Hrho_psi0d)) m_pi_Hrho_psi0d => ep_env%m_pi_Hrho_psi0d
       IF (PRESENT(psi1)) psi1 => ep_env%psi1
!FM     IF (PRESENT(nspins)) THEN
!FM        CPPrecondition(ASSOCIATED(ep_env%main_qs_env),cp_failure_level,routineP,error,failure)
!FM        IF (.not.failure) THEN
!FM           nspins=ep_env%main_qs_env%dft_control%nspins
!FM        END IF
!FM     END IF
       IF (PRESENT(preconditioner)) THEN
          preconditioner => ep_env%precond
       END IF
       IF (PRESENT(globenv)) globenv => ep_env%globenv
       IF (PRESENT(f_env_id)) f_env_id=ep_env%f_env_id
       IF (PRESENT(at2sub)) at2sub => ep_env%at2sub
    END IF
  END SUBROUTINE ep_env_get
  !***************************************************************************

  !!****f* ep_methods/ep_env_create [1.0] *
  !!
  !!   NAME
  !!     ep_env_create
  !!
  !!   FUNCTION
  !!     creates the environement to perform an ep perturbation
  !!
  !!   NOTES
  !!     Could take more advantage of what is calculated in main_qs_env
  !!     (for example the S matrix)
  !!     
  !!
  !!   INPUTS
  !!     - ep_env: the ep_env env to create
  !!     - main_qs_env: the qs_env of the full system
  !!     - projections: the projections into the subsystems
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     Fawzi Mohamed
  !!
  !!   MODIFICATION HISTORY
  !!     11.2002 created [fawzi]
  !!
  !!*** **********************************************************************
  SUBROUTINE ep_env_create(ep_env,globenv, error)
    TYPE(ep_env_type), POINTER               :: ep_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    CALL timeset(routineN,handle)
    CPPrecondition(.NOT.ASSOCIATED(ep_env),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       ALLOCATE(ep_env,stat=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)

       last_ep_env_id=last_ep_env_id+1
       ep_env%id_nr=last_ep_env_id
       ep_env%ref_count=1
       ep_env%f_env_id=-1
       ep_env%nspins=-1
       ep_env%nat=-1
       ep_env%nat_per_mol=-1
       ep_env%nmol=-1
       NULLIFY(ep_env%mol_envs,ep_env%sub_proj,ep_env%main_qs_env,&
            ep_env%main_p_env,ep_env%sub_p_env,ep_env%m_pi_Hrho_psi0d,&
            ep_env%psi1, ep_env%precond,ep_env%sub_nmo,ep_env%sub_nao,&
            ep_env%full_nmo,ep_env%full_nao,ep_env%input,ep_env%at2sub,&
            ep_env%force)
       CALL ep_energy_zero(ep_env%energy)
       ep_env%globenv => globenv
       CALL globenv_retain(globenv,error=error)
       ep_env%para_env=>globenv%para_env
       CALL cp_para_env_retain(ep_env%para_env,error=error)

       CALL ep_envs_add_ep_env(ep_env,error=error)
       CALL cp_ep_init(ep_env%id_nr,ierr)

       CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL timestop(handle)
  END SUBROUTINE ep_env_create

  !!****f* ep_types/ep_envs_get_ep_env *
  !!
  !!   NAME
  !!     ep_envs_get_ep_env
  !!
  !!   FUNCTION
  !!     returns the replica environment with the given id_nr
  !!
  !!   NOTES
  !!     -
  !!
  !!   ARGUMENTS
  !!     - id_nr: the id_nr of the requested ep_envs
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** **********************************************************************
  FUNCTION ep_envs_get_ep_env(id_nr) RESULT(res)
    INTEGER, INTENT(in)                      :: id_nr
    TYPE(ep_env_type), POINTER               :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_envs_get_ep_env', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

    NULLIFY(res)
    IF (module_initialized) THEN
       IF (ASSOCIATED(ep_envs)) THEN
          DO i=1,SIZE(ep_envs)
             IF (ep_envs(i)%ep_env%id_nr==id_nr) THEN
                res => ep_envs(i)%ep_env
                EXIT
             END IF
          END DO
       END IF
    END IF
  END FUNCTION ep_envs_get_ep_env
  !***************************************************************************

  !!****f* ep_types/ep_envs_add_ep_env *
  !!
  !!   NAME
  !!     ep_envs_add_ep_env
  !!
  !!   FUNCTION
  !!     adds the given ep_env to the list of controlled ep_envs.
  !!
  !!   NOTES
  !!     -
  !!
  !!   ARGUMENTS
  !!     - ep_env: the ep_env to add
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** **********************************************************************
  SUBROUTINE ep_envs_add_ep_env(ep_env,error)
    TYPE(ep_env_type), POINTER               :: ep_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_envs_add_ep_env', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure
    TYPE(ep_env_p_type), DIMENSION(:), &
      POINTER                                :: new_ep_envs
    TYPE(ep_env_type), POINTER               :: ep_env2

    failure=.FALSE.

    IF (ASSOCIATED(ep_env)) THEN
       ep_env2 => ep_envs_get_ep_env(ep_env%id_nr)
       IF (.NOT.ASSOCIATED(ep_env2)) THEN
          IF (module_initialized) THEN
             IF (ASSOCIATED(ep_envs)) THEN
                ALLOCATE(new_ep_envs(SIZE(ep_envs)+1),stat=stat)
                CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
                DO i=1,SIZE(ep_envs)
                   new_ep_envs(i)%ep_env => ep_envs(i)%ep_env
                END DO
                DEALLOCATE(ep_envs,stat=stat)
                CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
                ep_envs => new_ep_envs
             ELSE
                ALLOCATE(ep_envs(1),stat=stat)
                CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
             END IF
          ELSE
             ALLOCATE(ep_envs(1),stat=stat)
             CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
          END IF
          ep_envs(SIZE(ep_envs))%ep_env => ep_env
          module_initialized=.TRUE.
       END IF
    END IF
  END SUBROUTINE ep_envs_add_ep_env
  !***************************************************************************

  !!****f* ep_types/ep_envs_rm_ep_env *
  !!
  !!   NAME
  !!     ep_envs_rm_ep_env
  !!
  !!   FUNCTION
  !!     removes the given ep_env to the list of controlled ep_envs.
  !!
  !!   NOTES
  !!     -
  !!
  !!   ARGUMENTS
  !!     - ep_env: the ep_env to remove
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** **********************************************************************
  SUBROUTINE ep_envs_rm_ep_env(ep_env,error)
    TYPE(ep_env_type), POINTER               :: ep_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_envs_rm_ep_env', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ii, stat
    LOGICAL                                  :: failure
    TYPE(ep_env_p_type), DIMENSION(:), &
      POINTER                                :: new_ep_envs

    failure=.FALSE.

    IF (ASSOCIATED(ep_env)) THEN
       CPPrecondition(module_initialized,cp_failure_level,routineP,error,failure)
       ALLOCATE(new_ep_envs(SIZE(ep_envs)-1),stat=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
       ii=0
       DO i=1,SIZE(ep_envs)
          IF (ep_envs(i)%ep_env%id_nr/=ep_env%id_nr) THEN
             ii=ii+1
             new_ep_envs(ii)%ep_env => ep_envs(i)%ep_env
          END IF
       END DO
       CPPostcondition(ii==SIZE(new_ep_envs),cp_failure_level,routineP,error,failure)
       DEALLOCATE(ep_envs,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ep_envs => new_ep_envs
       IF (SIZE(ep_envs)==0) THEN
          DEALLOCATE(ep_envs,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF
  END SUBROUTINE ep_envs_rm_ep_env
  !***************************************************************************

  !!****f* ep_types/ep_env_calc_e_f *
  !!
  !!   NAME
  !!     ep_env_calc_e_f
  !!
  !!   FUNCTION
  !!     updates forces and energy in the given ep_env
  !!
  !!   NOTES
  !!     -
  !!
  !!   ARGUMENTS
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** **********************************************************************
  SUBROUTINE ep_env_calc_e_f(ep_env,calc_f,error)
    TYPE(ep_env_type), POINTER               :: ep_env
    LOGICAL, INTENT(in)                      :: calc_f
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_env_calc_e_f', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, my_calc_f
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (.NOT.failure) THEN
       IF (calc_f) THEN
          my_calc_f=1
       ELSE
          my_calc_f=0
       END IF
       CALL cp_ep_calc_e_f(ep_env%id_nr,my_calc_f,ierr)
       CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
  END SUBROUTINE ep_env_calc_e_f
  !***************************************************************************

  !!****f* ep_types/ep_dealloc_local_sm *
  !!
  !!   NAME
  !!     ep_dealloc_local_sm
  !!
  !!   FUNCTION
  !!     deallocates a local spare matrix
  !!
  !!   NOTES
  !!     -
  !!
  !!   ARGUMENTS
  !!     - local_sm: the matrix to deallocate
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** **********************************************************************
  SUBROUTINE ep_dealloc_local_sm(local_sm,error)
    TYPE(real_matrix_p_type), &
      DIMENSION(:, :), POINTER               :: local_sm
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_dealloc_local_sm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ispin, stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(local_sm),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       DO i=1,SIZE(local_sm,2)
          DO ispin=1,SIZE(local_sm,1)
             IF (ASSOCIATED(local_sm(ispin,i)%matrix)) THEN
                CALL deallocate_matrix(local_sm(ispin,i)%matrix)
             END IF
          END DO
       END DO
       DEALLOCATE(local_sm,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF
  END SUBROUTINE ep_dealloc_local_sm
  !***************************************************************************

  !!****f* ep_types/ep_force_create *
  !!
  !!   NAME
  !!     ep_force_create
  !!
  !!   FUNCTION
  !!     creates an ep_force_type
  !!
  !!   NOTES
  !!     -
  !!
  !!   ARGUMENTS
  !!     - force: the type to create
  !!     - nat: the number of atoms
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** **********************************************************************
  SUBROUTINE ep_force_create(force,nat,error)
    TYPE(ep_force_type), POINTER             :: force
    INTEGER, INTENT(in)                      :: nat
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_force_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(.NOT.ASSOCIATED(force),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       ALLOCATE(force,stat=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
    END IF
    IF (.NOT. failure) THEN
       last_force_id=last_force_id+1
       force%id_nr=last_force_id
       force%ref_count=1
       ALLOCATE(force%f0_internal(3,nat),stat=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
    END IF
  END SUBROUTINE ep_force_create
  !***************************************************************************

  !!****f* ep_types/ep_force_zero *
  !!
  !!   NAME
  !!     ep_force_zero
  !!
  !!   FUNCTION
  !!     zeroes the forces
  !!
  !!   NOTES
  !!     -
  !!
  !!   ARGUMENTS
  !!     - force: the force type to zero out
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** **********************************************************************
  SUBROUTINE ep_force_zero(force,error)
    TYPE(ep_force_type), POINTER             :: force
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_force_zero', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(force),cp_failure_level,routineP,error,failure)
    CPPrecondition(force%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       force%f0_internal=0._dp
    END IF
  END SUBROUTINE ep_force_zero
  !***************************************************************************

  !!****f* ep_types/ep_force_retain *
  !!
  !!   NAME
  !!     ep_force_retain
  !!
  !!   FUNCTION
  !!     retains the given force type
  !!
  !!   NOTES
  !!     -
  !!
  !!   ARGUMENTS
  !!     - force: the force type to retain
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** **********************************************************************
  SUBROUTINE ep_force_retain(force,error)
    TYPE(ep_force_type), POINTER             :: force
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_force_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(force),cp_failure_level,routineP,error,failure)
    CPPrecondition(force%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       force%ref_count=force%ref_count+1
    END IF
  END SUBROUTINE ep_force_retain
  !***************************************************************************

  !!****f* ep_types/ep_force_release *
  !!
  !!   NAME
  !!     ep_force_release
  !!
  !!   FUNCTION
  !!     retains the given force type
  !!
  !!   NOTES
  !!     -
  !!
  !!   ARGUMENTS
  !!     - force: the force type to retain
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** **********************************************************************
  SUBROUTINE ep_force_release(force,error)
    TYPE(ep_force_type), POINTER             :: force
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_force_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    IF (ASSOCIATED(force)) THEN
       CPPrecondition(force%ref_count>0,cp_failure_level,routineP,error,failure)
       force%ref_count=force%ref_count-1
       IF (force%ref_count==0) THEN
          DEALLOCATE(force%f0_internal,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          DEALLOCATE(force,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF
    NULLIFY(force)
  END SUBROUTINE ep_force_release
  !***************************************************************************

END MODULE ep_types
