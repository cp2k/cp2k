!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/pao_qs_env_methods [1.0] *
!!
!!   NAME
!!     pao_qs_env_methods
!!
!!   FUNCTION
!!     pao related methods that work on the qs_env methods
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE pao_qs_env_methods
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE cp_array_utils,                  ONLY: cp_2d_r_output
  USE cp_b_matrix_structure,           ONLY: cp_b_dims_type,&
                                             cp_b_matrix_struct_type,&
                                             cp_bmstruct_create,&
                                             cp_bmstruct_get,&
                                             cp_bmstruct_release
  USE cp_block_matrix,                 ONLY: cp_block_matrix_p_type,&
                                             cp_block_matrix_type,&
                                             cp_get_matrix,&
                                             cp_sp_create,&
                                             cp_sp_release
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_debug,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_matrix_utils,                 ONLY: cp_sm_set_dist_sparsity_from_sm
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE mathlib,                         ONLY: power_matrix,&
                                             symmetrize_matrix
  USE pao_bas_projection_types,        ONLY: pao_bas_projection_type
  USE pao_glob_angles_methods,         ONLY: pao_g_ang_did_change
  USE pao_types,                       ONLY: pao_env_get,&
                                             pao_env_type,&
                                             pao_glob_angles_type
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE sparse_matrix_types,             ONLY: add_block_node,&
                                             get_block_node,&
                                             real_matrix_p_type
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='pao_qs_env_methods'

  PUBLIC :: cp_get_b_projection, cp_did_change
  ! underlying routines
  PUBLIC :: pao_qs_env_get_b_projection, pao_qs_env_get, &
       pao_qs_env_did_change, pao_qs_env_get_atomic_ortho
  PUBLIC :: pao_bas_projection_type ! for the intel compiler
  INTERFACE cp_get_b_projection
     MODULE PROCEDURE pao_qs_env_get_b_projection
  END INTERFACE
  INTERFACE cp_did_change
     MODULE PROCEDURE pao_qs_env_did_change
  END INTERFACE
!***
!****************************************************************************
CONTAINS

!!****f* pao_qs_env_methods/pao_qs_env_get_b_projection [1.0] *
!!
!!   NAME
!!     pao_qs_env_get_b_projection
!!
!!   SYNOPSIS
!!     Function pao_qs_env_get_b_projection(qs_env, atom_nr, error)&
!!         Result(res)
!!       Type(pao_bas_projection_type), Pointer:: res
!!       Type(qs_environment_type), Intent (IN):: qs_env
!!       Integer, Intent (IN):: atom_nr
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function pao_qs_env_get_b_projection
!!
!!   FUNCTION
!!     returns the base projection projection on atom_nr
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     qs_env: the qs environement you get the data from
!!     atom_nr: the number of the atom
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     - 1.2002 created [fawzi]
!!     - 19.04.02 used new projection
!!     - 23.04.02 moved to qs
!!
!!*** **********************************************************************
  FUNCTION pao_qs_env_get_b_projection(qs_env, atom_nr, error) RESULT(res)
    TYPE(qs_environment_type), INTENT(in)    :: qs_env
    INTEGER, INTENT(in)                      :: atom_nr
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    TYPE(pao_bas_projection_type), POINTER   :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_qs_env_get_b_projection', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    failure=.FALSE.

    CALL get_qs_env(qs_env, particle_set=particle_set)
    CALL get_atomic_kind(particle_set(atom_nr)%atomic_kind,pao_bas_proj=res)
  END FUNCTION pao_qs_env_get_b_projection
!***************************************************************************

!!****f* pao_qs_env_methods/pao_qs_env_get [1.0] *
!!
!!   NAME
!!     pao_qs_env_get
!!
!!   SYNOPSIS
!!     Subroutine pao_qs_env_get(qs_env, global_env, pao_env,&
!!         pao_max_min_bas, pao_max_full_bas, pao_max_excl_bas,&
!!         pao_max_l_angles, pao_tot_min_bas, pao_tot_full_bas,&
!!         pao_tot_excl_bas, pao_min_density_m, pao_min_bas_m_struct,&
!!         pao_full_bas_m_struct, pao_min_full_bas_m_struct,&
!!         pao_full_min_bas_m_struct, pao_full_bas_diag_m_struct,&
!!         pao_angles_att, pao_w_subspace_pos_precision,&
!!         pao_w_in_subspace_precision, pao_atomic_ortho,&
!!         pao_atomic_ortho_valid, n_atoms, n_electrons, error)
!!       Type(qs_environment_type), Intent (IN), Target:: qs_env
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Type(pao_env_type), Pointer, Optional:: pao_env
!!       Integer, Intent (OUT), Optional:: pao_max_min_bas,&
!!         pao_max_full_bas, pao_max_excl_bas, pao_max_l_angles,&
!!         pao_tot_min_bas, pao_tot_full_bas, pao_tot_excl_bas, n_atoms,&
!!         n_electrons
!!       Real(KIND = dp), Optional, Intent (OUT)::&
!!         pao_w_subspace_pos_precision
!!       Real(KIND = dp), Optional, Intent (OUT)::&
!!         pao_w_in_subspace_precision
!!       Type(pao_glob_angles_type), Optional, Pointer:: pao_angles_att
!!       Type(cp_block_matrix_p_type), Dimension(:), Optional, Pointer::&
!!         pao_min_density_m
!!       Type(cp_b_matrix_struct_type), Pointer, Optional::&
!!         pao_min_bas_m_struct, pao_full_bas_m_struct,&
!!         pao_min_full_bas_m_struct, pao_full_min_bas_m_struct,&
!!         pao_full_bas_diag_m_struct
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!       Logical, Dimension(:), Pointer, Optional:: pao_atomic_ortho_valid
!!       Type(cp_block_matrix_type), Pointer, Optional:: pao_atomic_ortho
!!     End Subroutine pao_qs_env_get
!!
!!   FUNCTION
!!     returns pao related info from the qs_environment_type
!!
!!   NOTES
!!     use this function rather than get_atomic_kind_set: some values
!!     might be cached in pao_env in the future
!!
!!   INPUTS
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE pao_qs_env_get(qs_env, global_env, pao_env, pao_max_min_bas,&
     pao_max_full_bas, pao_max_excl_bas, pao_max_l_angles,&
     pao_tot_min_bas, pao_tot_full_bas, pao_tot_excl_bas,&
     pao_min_density_m, pao_min_bas_m_struct, pao_full_bas_m_struct,&
     pao_min_full_bas_m_struct, pao_full_min_bas_m_struct,&
     pao_full_bas_diag_m_struct,&
     pao_angles_att, pao_w_subspace_pos_precision,&
     pao_w_in_subspace_precision,pao_atomic_ortho,&
     pao_atomic_ortho_valid, n_atoms, n_electrons, error)
    TYPE(qs_environment_type), INTENT(in), &
      TARGET                                 :: qs_env
    TYPE(global_environment_type), &
      INTENT(in), TARGET                     :: global_env
    TYPE(pao_env_type), OPTIONAL, POINTER    :: pao_env
    INTEGER, INTENT(out), OPTIONAL :: pao_max_min_bas, pao_max_full_bas, &
      pao_max_excl_bas, pao_max_l_angles, pao_tot_min_bas, pao_tot_full_bas, &
      pao_tot_excl_bas
    TYPE(cp_block_matrix_p_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: pao_min_density_m
    TYPE(cp_b_matrix_struct_type), OPTIONAL, POINTER :: pao_min_bas_m_struct, &
      pao_full_bas_m_struct, pao_min_full_bas_m_struct, &
      pao_full_min_bas_m_struct, pao_full_bas_diag_m_struct
    TYPE(pao_glob_angles_type), OPTIONAL, &
      POINTER                                :: pao_angles_att
    REAL(KIND=dp), INTENT(out), OPTIONAL     :: pao_w_subspace_pos_precision, &
                                                pao_w_in_subspace_precision
    TYPE(cp_block_matrix_type), OPTIONAL, &
      POINTER                                :: pao_atomic_ortho
    LOGICAL, DIMENSION(:), OPTIONAL, POINTER :: pao_atomic_ortho_valid
    INTEGER, INTENT(out), OPTIONAL           :: n_atoms, n_electrons
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_qs_env_get', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: natoms
    LOGICAL                                  :: failure
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_subsystem_type), POINTER         :: subsys
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pao_env_type), POINTER              :: p_env
    TYPE(particle_list_type), POINTER        :: particles

  NULLIFY(atomic_kind_set, p_env, dft_control, subsys, particles)
  CALL get_qs_env(qs_env, subsys=subsys, error=error)
  CALL cp_subsys_get(subsys,particles=particles, error=error)
  natoms = particles%n_els
  failure=.FALSE.
  
  CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, pao_env=p_env)
  CPPrecondition(ASSOCIATED(p_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(p_env%initialized,cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     IF (PRESENT(pao_max_min_bas)) pao_max_min_bas=p_env%max_min_bas
     IF (PRESENT(pao_max_full_bas)) pao_max_full_bas=p_env%max_full_bas
     IF (PRESENT(pao_max_excl_bas)) pao_max_excl_bas=p_env%max_excl_bas
     IF (PRESENT(pao_tot_min_bas)) pao_tot_min_bas=p_env%tot_min_bas
     IF (PRESENT(pao_tot_full_bas)) pao_tot_full_bas=p_env%tot_full_bas
     IF (PRESENT(pao_tot_excl_bas)) pao_tot_excl_bas=p_env%tot_excl_bas
     IF (PRESENT(pao_max_l_angles)) pao_max_l_angles=p_env%max_l_angles
     CALL pao_env_get(p_env,&
          wanted_subspace_pos_precision=pao_w_subspace_pos_precision,&
          wanted_in_subspace_precision=pao_w_in_subspace_precision,&
          angles_att=pao_angles_att,min_density_m=pao_min_density_m,&
          min_bas_m_struct=pao_min_bas_m_struct,&
          min_full_bas_m_struct=pao_min_full_bas_m_struct,&
          full_min_bas_m_struct=pao_full_min_bas_m_struct,&
          full_bas_diag_m_struct=pao_full_bas_diag_m_struct,&
          full_bas_m_struct=pao_full_bas_m_struct, error=error)
     IF (PRESENT(pao_atomic_ortho_valid)) THEN
        IF (.NOT.ASSOCIATED(pao_env%atomic_ortho_valid)) THEN
           ALLOCATE(p_env%atomic_ortho_valid(natoms))
           p_env%atomic_ortho_valid=.FALSE.
        END IF
        pao_atomic_ortho_valid => p_env%atomic_ortho_valid
     END IF
     IF (PRESENT(pao_atomic_ortho)) THEN
        IF (.NOT.ASSOCIATED(p_env%atomic_ortho_valid)) THEN
           ALLOCATE(p_env%atomic_ortho_valid(natoms))
           p_env%atomic_ortho_valid=.FALSE.
        END IF
        IF (.NOT.ASSOCIATED(p_env%atomic_ortho)) THEN
           CALL cp_sp_create(p_env%atomic_ortho,&
                matrix_struct=p_env%full_min_bas_m_struct,&
                error=error)
        END IF
        pao_atomic_ortho => p_env%atomic_ortho
     END IF
     IF (PRESENT(n_atoms)) n_atoms = natoms
     IF (PRESENT(n_electrons)) THEN
        CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set,&
             dft_control=dft_control)
        CPPrecondition(ASSOCIATED(atomic_kind_set),cp_failure_level,routineP,error,failure)
        CPPrecondition(ASSOCIATED(dft_control),cp_failure_level,routineP,error,failure)
        IF (.NOT.failure) THEN
           CALL get_atomic_kind_set(atomic_kind_set,nelectron=n_electrons)
           n_electrons=n_electrons-dft_control%charge
        ELSE
           n_electrons=-1
        END IF
     END IF
  END IF
  IF (PRESENT(pao_env)) pao_env => p_env
END SUBROUTINE pao_qs_env_get
!***************************************************************************

!!****f* pao_qs_env_methods/pao_qs_env_did_change [1.0] *
!!
!!   NAME
!!     pao_qs_env_did_change
!!
!!   SYNOPSIS
!!     Subroutine pao_qs_env_did_change(qs_env, global_env, h_changed,&
!!         s_changed, angles_changed, full_reset, error)
!!       Type(qs_environment_type), Intent (INOUT):: qs_env
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Logical, Intent (IN), Optional:: h_changed, s_changed,&
!!         angles_changed, full_reset
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_qs_env_did_change
!!
!!   FUNCTION
!!     invalidates the cache when some data changes
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - qs_env: the environement that has changed
!!     - s_changed: if the overlap has changed (implies h_changed
!!     - angles_changed: if the angles that describe tha actual subspace 
!!       have changed
!!     - full_reset: if the structures should be regenerated from scratch
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
SUBROUTINE pao_qs_env_did_change(qs_env,global_env,&
     h_changed, s_changed,&
     angles_changed, full_reset, error)
    TYPE(qs_environment_type), INTENT(inout) :: qs_env
    TYPE(global_environment_type), &
      INTENT(in), TARGET                     :: global_env
    LOGICAL, INTENT(in), OPTIONAL            :: h_changed, s_changed, &
                                                angles_changed, full_reset
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_qs_env_did_change', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure, full_r, s_c
    TYPE(cp_b_dims_type), POINTER            :: b_dims
    TYPE(cp_b_matrix_struct_type), POINTER   :: old_m_struct
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(pao_env_type), POINTER              :: pao_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: overlap_rmpv

  failure=.FALSE.
  s_c=.FALSE.; full_r=.FALSE.
  
  CALL get_qs_env(qs_env,pao_env=pao_env)
  CPPrecondition(ASSOCIATED(pao_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(pao_env%initialized,cp_failure_level,routineP,error,failure)
  END IF
  IF (.NOT. failure) THEN
     IF (PRESENT(s_changed)) s_c=s_changed
     IF (PRESENT(full_reset)) full_r=full_reset
     IF (s_c.OR.full_r) THEN
        CALL get_qs_env(qs_env, s=overlap_rmpv)
        ! invalidate min_density_m
        pao_env%min_d_m_valid=.FALSE.
        DO i=1,SIZE(pao_env%min_density_m)
           CALL cp_sp_release(pao_env%min_density_m(i)%matrix,error=error)
        END DO
        ! invalidate matrix structure
        CALL cp_bmstruct_get(pao_env%min_bas_m_struct,&
            b_dims=b_dims, error=error)
        old_m_struct => pao_env%min_bas_m_struct
        CALL cp_bmstruct_create(pao_env%min_bas_m_struct,&
               symmetric=.TRUE.,&
               para_env=global_env%para_env,&
               b_dims=b_dims,&
               error=error)
        CALL cp_sm_set_dist_sparsity_from_sm(pao_env%min_bas_m_struct,&
             para_env=global_env%para_env,&
             sample_matrix=overlap_rmpv(1)%matrix,&
             error=error)
        CALL cp_bmstruct_release(old_m_struct,error=error)

        CALL cp_bmstruct_get(pao_env%full_bas_m_struct,&
             para_env=para_env, b_dims=b_dims, error=error)
        old_m_struct => pao_env%full_bas_m_struct
        CALL cp_bmstruct_create(pao_env%full_bas_m_struct,&
             symmetric=.FALSE.,&
             para_env=para_env,&
             b_dims=b_dims,&
             d_struct=pao_env%min_bas_m_struct%distribution,&
             sparsity=pao_env%min_bas_m_struct%sparsity,&
             error=error)
        CALL cp_bmstruct_release(old_m_struct,error=error)

        CALL cp_bmstruct_get(pao_env%full_min_bas_m_struct,&
             para_env=para_env, b_dims=b_dims, error=error)
        old_m_struct => pao_env%full_min_bas_m_struct
        CALL cp_bmstruct_create(pao_env%full_min_bas_m_struct,&
             symmetric=.FALSE.,&
             para_env=para_env,&
             b_dims=b_dims,&
             d_struct=pao_env%min_bas_m_struct%distribution,&
             sparsity=pao_env%min_bas_m_struct%sparsity,&
             error=error)
        CALL cp_bmstruct_release(old_m_struct,error=error)

        CALL cp_bmstruct_get(pao_env%min_full_bas_m_struct,&
             para_env=para_env, b_dims=b_dims, error=error)
        old_m_struct => pao_env%min_full_bas_m_struct
        CALL cp_bmstruct_create(pao_env%min_full_bas_m_struct,&
             symmetric=.FALSE.,&
             para_env=para_env,&
             b_dims=b_dims,&
             d_struct=pao_env%min_bas_m_struct%distribution,&
             sparsity=pao_env%min_bas_m_struct%sparsity,&
             error=error)
        CALL cp_bmstruct_release(old_m_struct,error=error)
        
        DO i=1,SIZE(pao_env%min_density_m)
           CALL cp_sp_create(pao_env%min_density_m(i)%matrix,&
                matrix_struct= pao_env%min_bas_m_struct,&
                error=error)
        END DO
        ! invalidate atomic orthogonalisation
        IF (ASSOCIATED(pao_env%atomic_ortho_valid)) &
             pao_env%atomic_ortho_valid=.FALSE.
        ! in general the atom distribution changes
        CALL cp_sp_release(pao_env%atomic_ortho,error=error)
        CALL pao_qs_env_preinit_atom_ortho(qs_env, global_env, error=error)
     END IF
     CALL pao_g_ang_did_change(pao_env%angles_att,qs_env=qs_env,&
          global_env=global_env, h_changed=h_changed,&
          s_changed=s_changed, angles_changed=angles_changed,&
          full_reset=full_reset, error=error)
  END IF
END SUBROUTINE pao_qs_env_did_change
!***************************************************************************

!!****f* pao_qs_env_methods/pao_qs_env_get_atomic_ortho [1.0] *
!!
!!   NAME
!!     pao_qs_env_get_atomic_ortho
!!
!!   SYNOPSIS
!!     Function pao_qs_env_get_atomic_ortho(qs_env, global_env, atom_nr,&
!!         error) Result(res)
!!       Type(qs_environment_type), Intent (IN), Target:: qs_env
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Integer, Intent (IN):: atom_nr
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!       Real(KIND = dp), Dimension(:,:), Pointer:: res
!!     End Function pao_qs_env_get_atomic_ortho
!!
!!   FUNCTION
!!     returns the atomic orthogonalization matrix, calculating it
!!     if necessary.
!!
!!   NOTES
!!     could cache tmp_m (but submatrix indexing on most compilers is
!!     not optimized)
!!
!!   INPUTS
!!     - qs_env: the qs_env that contains the overlap matrix
!!     - global_env: the global environment (for the communication)
!!     - atom_nr: the atom you want the atomic orthogonalization from
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     6.2002 created [fawzi]
!!
!!*** **********************************************************************
FUNCTION pao_qs_env_get_atomic_ortho(qs_env,global_env,atom_nr,error)&
     RESULT(res)
    TYPE(qs_environment_type), INTENT(in), &
      TARGET                                 :: qs_env
    TYPE(global_environment_type), &
      INTENT(in), TARGET                     :: global_env
    INTEGER, INTENT(in)                      :: atom_nr
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_qs_env_get_atomic_ortho', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: n_atoms, stat
    INTEGER, SAVE                            :: iter = 0
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: eps_eigval
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: tmp_m
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block_val
    TYPE(cp_block_matrix_type), POINTER      :: atomic_ortho
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(pao_env_type), POINTER              :: pao_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: s

  failure=.FALSE.
  NULLIFY(res, atomic_ortho, block_val, pao_env, s)

  CALL pao_qs_env_get(qs_env,global_env,pao_env=pao_env,n_atoms=n_atoms)
  CPPrecondition(ASSOCIATED(pao_env),cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     IF (.NOT.ASSOCIATED(pao_env%atomic_ortho_valid)) THEN
        ALLOCATE(pao_env%atomic_ortho_valid(n_atoms),stat=stat)
        CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
        IF (.NOT.failure) THEN
           pao_env%atomic_ortho_valid=.FALSE.
        END IF
     END IF
     IF (.NOT.ASSOCIATED(pao_env%atomic_ortho)) THEN
        CALL cp_sp_create(pao_env%atomic_ortho,&
             matrix_struct=pao_env%full_bas_diag_m_struct,&
             error=error)
     END IF
     atomic_ortho => pao_env%atomic_ortho
  END IF
  CPPrecondition(ASSOCIATED(atomic_ortho),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(pao_env%atomic_ortho_valid),cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     IF (pao_env%atomic_ortho_valid(atom_nr)) THEN
        CALL get_block_node(cp_get_matrix(atomic_ortho),&
             BLOCK=res,block_row=atom_nr,&
             block_col=atom_nr)
        CPPostcondition(ASSOCIATED(res),cp_failure_level,routineP,error,failure)
     ELSE
        CALL get_qs_env(qs_env,s=s)
        CALL get_block_node(cp_get_matrix(atomic_ortho),&
             BLOCK=res,block_row=atom_nr,&
             block_col=atom_nr)
        IF (.NOT.ASSOCIATED(res)) THEN
           CALL add_block_node(cp_get_matrix(atomic_ortho),&
                block_row=atom_nr,block_col=atom_nr,BLOCK=res)
        END IF
        CPPostcondition(ASSOCIATED(res),cp_failure_level,routineP,error,failure)
        IF (.NOT.failure) THEN
           CPPrecondition(SIZE(res,1)==SIZE(res,2),cp_failure_level,routineP,error,failure)
        END IF
        IF (.NOT.failure) THEN
           ALLOCATE(tmp_m(SIZE(res,1),SIZE(res,2)),stat=stat)
           CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
        END IF
        IF (.NOT.failure) THEN
           CALL get_block_node(s(1)%matrix,block_row=atom_nr,&
                block_col=atom_nr,BLOCK=block_val)
           CPPrecondition(SIZE(block_val,1)==SIZE(res,1),cp_failure_level,routineP,error,failure)
           CPPrecondition(SIZE(block_val,2)==SIZE(res,2),cp_failure_level,routineP,error,failure)
           IF (.NOT.failure) THEN
              tmp_m=block_val
              eps_eigval= 10 * EPSILON(0.0_dp) ! to do: use the scf value
              CALL power_matrix(tmp_m,res,-0.5_dp,threshold=eps_eigval)
              CALL symmetrize_matrix(res,"upper_to_lower")
              IF (cp_debug.AND.debug_this_module) THEN
                 iter=iter+1
                 logger => cp_error_get_logger(error)
                 CALL cp_2d_r_output(logger,&
                      outputName="atomic_ortho",&
                      fromWhere=routineP, iter=iter,array=res,&
                      comment="atom_nr="//cp_to_string(atom_nr),&
                      error=error)
              END IF
           END IF
           DEALLOCATE(tmp_m,stat=stat) ! just to be sure...
           CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
           pao_env%atomic_ortho_valid(atom_nr)=.TRUE.
        END IF
     END IF
  END IF
  IF (failure) NULLIFY(res)
END FUNCTION pao_qs_env_get_atomic_ortho
!***************************************************************************

!!****f* pao_qs_env_methods/pao_qs_env_preinit_atom_ortho [1.0] *
!!
!!   NAME
!!     pao_qs_env_preinit_atom_ortho
!!
!!   FUNCTION
!!     If pao is run in parallel shares the ortho
!!
!!   NOTES
!!     could cache tmp_m (but submatrix indexing on most compilers is
!!     not optimized)
!!
!!   INPUTS
!!     - qs_env: the qs_env that contains the overlap matrix
!!     - global_env: the global environment (for the communication)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     6.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE pao_qs_env_preinit_atom_ortho(qs_env,global_env,error)
    TYPE(qs_environment_type), INTENT(in), &
      TARGET                                 :: qs_env
    TYPE(global_environment_type), &
      INTENT(in), TARGET                     :: global_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'pao_qs_env_preinit_atom_ortho', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: n_atoms, stat
    LOGICAL                                  :: failure
    TYPE(cp_block_matrix_type), POINTER      :: atomic_ortho
    TYPE(pao_env_type), POINTER              :: pao_env

  failure=.FALSE.
  NULLIFY(atomic_ortho, pao_env)

  CALL pao_qs_env_get(qs_env,global_env,pao_env=pao_env,n_atoms=n_atoms)
  CPPrecondition(ASSOCIATED(pao_env),cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     IF (.NOT.ASSOCIATED(pao_env%atomic_ortho_valid)) THEN
        ALLOCATE(pao_env%atomic_ortho_valid(n_atoms),stat=stat)
        CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
        IF (.NOT.failure) THEN
           pao_env%atomic_ortho_valid=.FALSE.
        END IF
     END IF
     IF (.NOT.ASSOCIATED(pao_env%atomic_ortho)) THEN
        CALL cp_sp_create(pao_env%atomic_ortho,&
             matrix_struct=pao_env%full_bas_diag_m_struct,&
             error=error)
     END IF
     atomic_ortho => pao_env%atomic_ortho
  END IF
  CPPrecondition(ASSOCIATED(atomic_ortho),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(pao_env%atomic_ortho_valid),cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
#ifdef __parallel
!!FM     ! stupid parallelization, every processor has everything
!!FM     do i=1,n_atoms
!!FM        if (atomic_ortho%matrix_struct%distribution%block_owner(i,i) ==&
!!FM             atomic_ortho%global_env%mepos) then
!!FM           call pao_qs_env_get_atomic_ortho(qs_env,&
!!FM                global_env=atomic_ortho%global_env,&
!!FM                atom_nr=i,error=error)
!!FM        end if
!!FM     end do
!!FM     call cp_sp_block_handle_init(block_handle,error=error)
!!FM     do i=1,nblock_rows
!!FM        call cp_sp_get_local_block(matrix=atomic_ortho,&
!!FM             block_row=i, block_col=i,&
!!FM             block_handle=block_handle,alloc_block=.true.,error=error,&
!!FM             owner=owner)
!!FM        call mp_bcast(block_handle%block_val,owner,atomic_ortho%global_env%group)
!!FM        call cp_sp_block_handle_release(block_handle,error=error)
!!FM     end do
#endif
  END IF
END SUBROUTINE pao_qs_env_preinit_atom_ortho
!***************************************************************************

END MODULE pao_qs_env_methods
