!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/pao_obj_function [1.0] *
!!
!!   NAME
!!     pao_obj_function
!!
!!   FUNCTION
!!     contain the objective function to be optimized from the pao procedure
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!*****************************************************************************
MODULE pao_obj_function
!  use pao_env_methods
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE cp_array_utils,                  ONLY: cp_1d_r_output
  USE cp_b_matrix_structure,           ONLY: cp_b_matrix_struct_type
  USE cp_block_matrix,                 ONLY: cp_block_matrix_p_type,&
                                             cp_block_matrix_type,&
                                             cp_get_matrix,&
                                             cp_sp_alloc_nonsparse_blocks,&
                                             cp_sp_create,&
                                             cp_sp_release
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm,&
                                             cp_fm_symm,&
                                             cp_fm_upper_to_full
  USE cp_fm_diag,                      ONLY: cp_fm_power,&
                                             cp_fm_syevx
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_release,&
                                             cp_fm_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_get_default_logger,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_matrix_utils,                 ONLY: cp_d_struct_create,&
                                             cp_dealloc_ref,&
                                             cp_get,&
                                             cp_init,&
                                             cp_init_diagonal_blocks,&
                                             cp_matrix_block_iterator,&
                                             cp_next,&
                                             cp_set_to,&
                                             cp_sm_b_i_dealloc_ref,&
                                             cp_sm_b_i_init,&
                                             cp_sm_output
  USE cp_output_handling,              ONLY: cp_write_output
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: copy_fm_to_sm,&
                                             copy_sm_to_fm
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE pao_glob_angles_methods,         ONLY: pao_g_ang_did_change,&
                                             pao_g_ang_set
  USE pao_l_angles_methods,            ONLY: cp_get_g_index,&
                                             cp_set,&
                                             pao_l_ang_dealloc_ref,&
                                             pao_l_ang_init,&
                                             pao_local_angles_type
  USE pao_proj_methods,                ONLY: cp_m2f_with_B_M_pi,&
                                             pao_g_ang_get_min_m
  USE pao_qs_env_methods,              ONLY: pao_qs_env_get,&
                                             pao_qs_env_get_atomic_ortho
  USE pao_types,                       ONLY: cp_get,&
                                             pao_env_type,&
                                             pao_glob_angles_type
  USE pao_unitary_t,                   ONLY: cp_dealloc_ref,&
                                             cp_init,&
                                             cp_next,&
                                             pao_calc_u_du_get,&
                                             pao_calc_u_du_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_matrix_pools,                 ONLY: qs_matrix_pools_type
  USE qs_mo_methods,                   ONLY: calculate_density_matrix
  USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                             deallocate_mo_set,&
                                             get_mo_set,&
                                             init_mo_set,&
                                             mo_set_p_type,&
                                             mo_set_type
  USE sparse_matrix_types,             ONLY: get_block_node,&
                                             get_matrix_info,&
                                             real_matrix_p_type,&
                                             real_matrix_type,&
                                             symmetrise_diagonal_blocks
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE
  PRIVATE
  
  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN="pao_obj_function"
  INTEGER, PRIVATE, SAVE :: last_gf_env_id=0
  
  ! type
  PUBLIC :: pao_gerd_functional_type

  ! generic functions
  PUBLIC :: cp_create, cp_retain, cp_release, cp_eval_at, cp_valid

  ! underlying routines
  PUBLIC :: pao_gerd_f_eval_at, pao_gerd_f_valid,&
       pao_gerd_f_create, pao_gerd_f_retain, pao_gerd_f_release
  
  INTERFACE cp_create
     MODULE PROCEDURE pao_gerd_f_create
  END INTERFACE
  INTERFACE cp_release
     MODULE PROCEDURE pao_gerd_f_release
  END INTERFACE
  INTERFACE cp_retain
     MODULE PROCEDURE pao_gerd_f_retain
  END INTERFACE
  INTERFACE cp_valid
     MODULE PROCEDURE pao_gerd_f_valid
  END INTERFACE
  INTERFACE cp_eval_at
     MODULE PROCEDURE pao_gerd_f_eval_at
  END INTERFACE

!!***
!****************************************************************************

!!****s* pao_obj_function/pao_gerd_functional_type [1.0] *
!!
!!   NAME
!!     pao_gerd_functional_type
!!
!!   FUNCTION
!!     the environement for the gert functional, basically all the
!!     temporary space, this way it is not allocated everytime,
!!     and it is not a save attribute (so it is multithread safe;)
!!
!!   NOTES
!!     the functional is:
!!       f = Tr[P (N U i)^* H (N U i)] + Tr[V(P - P (N U i)^* S (N U i) P)]
!!     where P is the density matrix in the minimal basis, 
!!     N an atomic ortogonalization term ( S_{at}^{-1/2} )
!!     in the full basis, U a rotation between the full basis and the
!!     min basis, i an injection into the min basis,
!!     H the hamiltonian in the full basis, V langrange multiplicators
!!     for the idempotency of P, S is the the overlap matrix in the 
!!     full basis, and Tr[A] is the trace of a with respect to the
!!     standard scalar product (Tr[A]=sum(A_{i,i})).
!!     This functional was proposed by Gerd Breghold in his dissertation
!!     at the Max-Planck-Institut, and is derived from the functional used 
!!     by Lee and Head-Gordon in Chem. Phys. 107,p 9085, (1997)
!!
!!   ATTRIBUTES
!!     - ref_cont: reference count (see doc/ReferenceCounting.html)
!!     - id_nr: identification number (unique)
!!     - pao_env: a pointer to the actual pao_env. It is assumed that the
!!       size of the min and full basis, and the pao env do not change
!!       between sucessive evaluations
!!     - min_hamiltonian_m: cp_fm_type matrix with the hamiltonian in the 
!!       min basis
!!     - min_overlap_m: cp_fm_type matrix with the overlap in the min basis
!!     - density_m: cp_fm_type matrix with the density matrix
!!     - work_m: cp_fm_type matrix in the small basis
!!     - alpha_mo: molecular orbitals
!!     - scf_work1, scf_work2, scf_work3: cp_fm_type matrices in the big basis
!!     - diag_matrix: diagonal sparse matrix
!!     - full_work_m: sparse matrix in the full basis
!!     - min_work_m: sparse matrix in the minimal basis
!!     - d_struct: the distribution structure of the diagonal blocks between
!!       the processors (d_struct(i) gives on which processor is blocks i)
!!     - n_eval: the number of function evaluations
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  TYPE pao_gerd_functional_type
     LOGICAL :: initialized, should_dealloc_mo, rebuild_rho
     INTEGER :: n_eval, ref_count, id_nr
     TYPE(qs_environment_type), POINTER :: qs_env
     TYPE(cp_fm_type), POINTER :: min_blacs_1, min_blacs_2,&
          min_blacs_3
     TYPE(cp_fm_type), POINTER :: scf_work1, scf_work2, scf_work3
     TYPE(mo_set_p_type), DIMENSION(:), POINTER :: mo_set
     TYPE(cp_block_matrix_type), POINTER :: diag_matrix, full_work_m,&
          min_work_m
     INTEGER, DIMENSION(:), POINTER :: d_struct
     TYPE(global_environment_type), POINTER :: global_env
  END TYPE pao_gerd_functional_type
!!***
!****************************************************************************
CONTAINS

!!****f* pao_obj_function/pao_gerd_f_create [1.0] *
!!
!!   NAME
!!     pao_gerd_f_create
!!
!!   SYNOPSIS
!!     Subroutine pao_gerd_f_create(gf_env, qs_env, global_env, scf_work1,&
!!         scf_work2, scf_work3, mo_orbital, should_dealloc_mo, error)
!!       Type(pao_gerd_functional_type), Pointer:: gf_env
!!       Type(qs_environment_type), Intent (IN), Target:: qs_env
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Type(cp_fm_type), Intent (IN), Target:: scf_work1, scf_work2,&
!!         scf_work3
!!       Type(mo_set_p_type), Dimension(:), Target, Intent (IN), Optional::&
!!         mo_orbital
!!       Logical, Intent (IN), Optional:: should_dealloc_mo
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_gerd_f_create
!!
!!   FUNCTION
!!     initializes gerd s functional
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     gf_env: the object to be initialized
!!     pao_env: the pao environement of this functional
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling.
!!
!!     see pao_gerd_functional_type for the other variables
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE pao_gerd_f_create(gf_env,qs_env,global_env, scf_work1,scf_work2,&
       scf_work3, mo_orbital, should_dealloc_mo, error)
    TYPE(pao_gerd_functional_type), POINTER  :: gf_env
    TYPE(qs_environment_type), INTENT(in), &
      TARGET                                 :: qs_env
    TYPE(global_environment_type), &
      INTENT(in), TARGET                     :: global_env
    TYPE(cp_fm_type), INTENT(in), TARGET     :: scf_work1, scf_work2, &
                                                scf_work3
    TYPE(mo_set_p_type), DIMENSION(:), &
      INTENT(in), OPTIONAL, TARGET           :: mo_orbital
    LOGICAL, INTENT(in), OPTIONAL            :: should_dealloc_mo
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'pao_obj_function:pao_gerd_f_create'

    INTEGER                                  :: i, n_atoms, nelectron, &
                                                ntot_full_bas, ntot_min_bas, &
                                                stat
    LOGICAL                                  :: failure
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_b_matrix_struct_type), POINTER   :: full_bas_m_struct, &
                                                min_bas_m_struct
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(pao_env_type), POINTER              :: pao_env
    TYPE(qs_matrix_pools_type), POINTER      :: mpools
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: full_overlap_m
    TYPE(real_matrix_type), POINTER          :: diagMatrix

    NULLIFY(diagMatrix, min_bas_m_struct, full_bas_m_struct, pao_env,&
         atomic_kind_set,mpools,fm_struct_tmp)

    failure=.FALSE.

    ALLOCATE(gf_env,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       gf_env%ref_count=1
       last_gf_env_id=last_gf_env_id+1
       gf_env%id_nr=last_gf_env_id
       gf_env%qs_env => qs_env
       gf_env%global_env => global_env
       gf_env%n_eval = 0
    END IF
    IF (.NOT.failure) THEN
       CALL get_qs_env(gf_env%qs_env, s=full_overlap_m, pao_env=pao_env,&
            atomic_kind_set=atomic_kind_set, mpools=mpools)
       CALL pao_qs_env_get(gf_env%qs_env, global_env=global_env,&
            pao_tot_full_bas=ntot_full_bas,&
            pao_tot_min_bas=ntot_min_bas,&
            pao_min_bas_m_struct=min_bas_m_struct,&
            pao_full_bas_diag_m_struct=full_bas_m_struct, n_atoms=n_atoms,&
            error=error)
       CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ntot_min_bas, &
                                ncol_global=ntot_min_bas,  &
                                para_env=scf_work1%matrix_struct%para_env, &
                                context=scf_work1%matrix_struct%context)
       CALL cp_fm_create(gf_env%min_blacs_1,fm_struct_tmp,name='min blacs 1')
       CALL cp_fm_create(gf_env%min_blacs_3,fm_struct_tmp,name='min blacs 3')
       CALL cp_fm_create(gf_env%min_blacs_3,fm_struct_tmp,name='min blacs 3')
       CALL cp_fm_struct_release(fm_struct_tmp)
       gf_env%scf_work1 => scf_work1
       gf_env%scf_work2 => scf_work2
       gf_env%scf_work3 => scf_work3
       IF (PRESENT(mo_orbital)) THEN
          ALLOCATE(gf_env%mo_set(SIZE(mo_orbital)),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO i=1,SIZE(mo_orbital)
             gf_env%mo_set(i)%mo_set => mo_orbital(i)%mo_set
          END DO
          gf_env%should_dealloc_mo=.FALSE.
       ELSE
          CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
               nelectron=nelectron)
          ALLOCATE(gf_env%mo_set(1),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          CALL allocate_mo_set(gf_env%mo_set(1)%mo_set,nao=ntot_min_bas,&
               nmo=nelectron/2,&
               nelectron=nelectron,maxocc=2.0_dp)
          CALL init_mo_set(gf_env%mo_set(1)%mo_set,&
               mpools%ao_mo_fm_pools(1)%pool,&
               "pao_min_bas_mo")
          gf_env%should_dealloc_mo=.TRUE.
       END IF
       IF (PRESENT(should_dealloc_mo)) &
            gf_env%should_dealloc_mo=should_dealloc_mo
       CALL cp_sp_create(gf_env%diag_matrix,&
            matrix_struct=full_bas_m_struct)
       CALL cp_d_struct_create(gf_env%d_struct,&
            matrix=full_overlap_m(1)%matrix,&
            para_env=global_env%para_env,error=error)
       diagMatrix => cp_get_matrix(gf_env%diag_matrix) ! for xlf...
       CALL cp_init_diagonal_blocks(diagMatrix,&
            para_env=global_env%para_env,&
            default_value=0.0_dp,&
            d_struct=gf_env%d_struct)
       CALL cp_sp_create(gf_env%min_work_m,&
            matrix_struct=min_bas_m_struct)
       CALL cp_sp_alloc_nonsparse_blocks(gf_env%min_work_m,error=error)
       CALL cp_sp_create(gf_env%full_work_m,&
            matrix_struct=full_bas_m_struct)
       gf_env%rebuild_rho=.TRUE.
       gf_env%initialized=.TRUE.
    END IF
  END SUBROUTINE pao_gerd_f_create
!***************************************************************************

!!****f* pao_obj_function/pao_gerd_f_retain [1.0] *
!!
!!   NAME
!!     pao_gerd_f_retain
!!
!!   SYNOPSIS
!!     Subroutine pao_gerd_f_retain(gf_env, error)
!!       Type(pao_gerd_functional_type), Pointer:: gf_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_gerd_f_retain
!!
!!   FUNCTION
!!     retains the given object
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - gf_env: the object to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE pao_gerd_f_retain(gf_env,error)
    TYPE(pao_gerd_functional_type), POINTER  :: gf_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_gerd_f_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(gf_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPreconditionNoFail(gf_env%ref_count>0,cp_failure_level,routineP,error)
     gf_env%ref_count=gf_env%ref_count+1
  END IF
END SUBROUTINE pao_gerd_f_retain
!***************************************************************************

!!****f* pao_obj_function/pao_gerd_f_release [1.0] *
!!
!!   NAME
!!     pao_gerd_f_release
!!
!!   SYNOPSIS
!!     Subroutine pao_gerd_f_release(gf_env, error)
!!       Type(pao_gerd_functional_type), Pointer:: gf_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_gerd_f_release
!!
!!   FUNCTION
!!     releases the functional object
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     gf_env: the functional to be released
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE pao_gerd_f_release(gf_env, error)
    TYPE(pao_gerd_functional_type), POINTER  :: gf_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'pao_obj_function:pao_gerd_f_release'

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(gf_env)) THEN
       CPPreconditionNoFail(gf_env%ref_count>0,cp_failure_level,routineP,error)
       gf_env%ref_count=gf_env%ref_count-1
       IF (gf_env%ref_count==0) THEN
          CALL cp_fm_release(gf_env%min_blacs_1)
          CALL cp_fm_release(gf_env%min_blacs_2)
          CALL cp_fm_release(gf_env%min_blacs_3)
          NULLIFY(gf_env%scf_work1, gf_env%scf_work2, gf_env%scf_work3)
          CALL cp_sp_release(gf_env%min_work_m, error=error)
          CALL cp_sp_release(gf_env%full_work_m, error=error)
          CALL cp_sp_release(gf_env%diag_matrix, error=error)
          DEALLOCATE(gf_env%d_struct, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          IF (ASSOCIATED(gf_env%mo_set)) THEN
             IF (gf_env%should_dealloc_mo) THEN
                DO i=1,SIZE(gf_env%mo_set)
                   CALL deallocate_mo_set(gf_env%mo_set(i)%mo_set)
                END DO
             END IF
             DEALLOCATE(gf_env%mo_set,stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
          NULLIFY(gf_env%global_env)
          DEALLOCATE(gf_env,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF
  END SUBROUTINE pao_gerd_f_release
!***************************************************************************

!!****f* pao_obj_function/pao_gerd_f_valid [1.0] *
!!
!!   NAME
!!     pao_gerd_f_valid
!!
!!   SYNOPSIS
!!     Function pao_gerd_f_valid(gf_env, error) Result(res)
!!       Logical:: res
!!       Type(pao_gerd_functional_type), Pointer:: gf_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function pao_gerd_f_valid
!!
!!   FUNCTION
!!     returns true if the functional object is valid (minimal testing)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     gf_env: the functional to test
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  FUNCTION pao_gerd_f_valid(gf_env, error) RESULT(res)
    TYPE(pao_gerd_functional_type), POINTER  :: gf_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'pao_obj_function:pao_gerd_f_valid'

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(gf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(gf_env%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT.(ASSOCIATED(gf_env%scf_work1).AND.&
         ASSOCIATED(gf_env%scf_work2).AND.&
         ASSOCIATED(gf_env%scf_work3).AND.&
         ASSOCIATED(gf_env%qs_env).AND.&
         ASSOCIATED(gf_env%diag_matrix).AND.&
         ASSOCIATED(gf_env%full_work_m).AND.&
         ASSOCIATED(gf_env%min_work_m).AND.&
         ASSOCIATED(gf_env%min_blacs_1).AND.&
         ASSOCIATED(gf_env%min_blacs_2).AND.&
         ASSOCIATED(gf_env%min_blacs_3))) THEN
       CPPrecondition(ASSOCIATED(gf_env%scf_work1),cp_warning_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(gf_env%scf_work2),cp_warning_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(gf_env%scf_work3),cp_warning_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(gf_env%qs_env),cp_warning_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(gf_env%diag_matrix),cp_warning_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(gf_env%full_work_m),cp_warning_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(gf_env%min_work_m),cp_warning_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(gf_env%min_blacs_1),cp_warning_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(gf_env%min_blacs_2),cp_warning_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(gf_env%min_blacs_3),cp_warning_level,routineP,error,failure)
    END IF
    res=.NOT.failure
  END FUNCTION pao_gerd_f_valid
!***************************************************************************

!!****f* pao_obj_function/pao_gerd_f_eval_at [1.0] *
!!
!!   NAME
!!     pao_gerd_f_eval_at
!!
!!   SYNOPSIS
!!     Subroutine pao_gerd_f_eval_at(gf_env, x, master, para_env, f,&
!!         gradient, error)
!!       Type(pao_gerd_functional_type), Pointer:: gf_env
!!       Real(KIND = dp), Intent (IN), Dimension(:):: x
!!       Real(KIND = dp), Intent (OUT), Optional:: f
!!       Real(KIND = dp), Optional, Pointer, Dimension(:):: gradient
!!       Integer, Intent (IN):: master
!!       Type(cp_para_env_type), Pointer:: para_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_gerd_f_eval_at
!!
!!   FUNCTION
!!     evaluates the functional at one point (i.e. for the given angles
!!     and unitary transformation
!!
!!   NOTES
!!     min_p, and the density are twice the P used for the pao functional
!!     (closed shell electrons)
!!
!!   INPUTS
!!     gf_env: the gerd functional environement
!!     angles: angles defining the point where the functional is
!!            evaluated. If not given uses the ones in pao env
!!     f: the value of the functional
!!     df_du: the value of the derivative with respect to u
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE pao_gerd_f_eval_at(gf_env,x, master, para_env, f, gradient, error)
    TYPE(pao_gerd_functional_type), POINTER  :: gf_env
    REAL(KIND=dp), DIMENSION(:), INTENT(in)  :: x
    INTEGER, INTENT(in)                      :: master
    TYPE(cp_para_env_type), POINTER          :: para_env
    REAL(KIND=dp), INTENT(out), OPTIONAL     :: f
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: gradient
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'pao_obj_function:pao_gerd_f_eval_at'

    INTEGER :: block_col, block_row, excl_bas, full_bas_size, min_bas, &
      min_bas_size, size_null_space, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mo_eigenvalues, &
                                                mo_occupation_numbers
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block_val, du
    TYPE(cp_block_matrix_p_type), &
      DIMENSION(:), POINTER                  :: h_matrixes, min_d_smpv, &
                                                s_matrixes
    TYPE(cp_block_matrix_type), POINTER      :: h_matrix, s_matrix
    TYPE(cp_error_type)                      :: iError
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_matrix_block_iterator)           :: block_iterator
    TYPE(global_environment_type), POINTER   :: g_env
    TYPE(pao_calc_u_du_type)                 :: calc_u_du
    TYPE(pao_env_type), POINTER              :: pao_env
    TYPE(pao_glob_angles_type), POINTER      :: angles
    TYPE(pao_local_angles_type)              :: local_angles
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: full_h_rmpv, full_s_rmpv
    TYPE(real_matrix_type), POINTER          :: diagMatrix, min_density_m

!         full_hamiltonian_m
!type(cp_block_matrix_type) :: full_pao_m, min_pao_m

    CALL cp_error_init(iError, template_error=error)
    NULLIFY(diagMatrix, mo_eigenvalues, mo_occupation_numbers,min_density_m,&
         min_d_smpv,full_h_rmpv,full_s_rmpv, h_matrixes, s_matrixes, h_matrix,&
         s_matrix, angles, pao_env)
    g_env => gf_env%global_env
    logger => cp_error_get_logger(error)
    failure=.FALSE.
    
    CPAssert(cp_valid(gf_env,error=error),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       gf_env%n_eval = gf_env%n_eval+1
       CALL get_qs_env(gf_env%qs_env,pao_env=pao_env,&
            k=full_h_rmpv, s=full_s_rmpv)
       CPPrecondition(ASSOCIATED(pao_env),cp_failure_level,routineP,error,failure)
       IF (cp_debug .AND. debug_this_module) THEN
          CALL cp_1d_r_output(logger=logger,&
               outputName='obj_f_x',fromWhere=routineP, iter=gf_env%n_eval,&
               array=x,error=error)
       END IF
    END IF
    IF (.NOT.failure) THEN
       CALL pao_qs_env_get(gf_env%qs_env,global_env=g_env,pao_angles_att=angles,&
            pao_min_density_m=min_d_smpv,&
            pao_tot_full_bas=full_bas_size,&
            pao_tot_min_bas=min_bas_size, error=error)
       CPPostcondition(ASSOCIATED(angles),cp_failure_level,routineP,error,failure)
       CPPrecondition(SIZE(min_d_smpv)==1,cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT.failure) THEN
       min_density_m => min_d_smpv(1)%matrix%sm
       CALL pao_g_ang_set(angles,qs_env=gf_env%qs_env,global_env=g_env,&
            angles=x,error=error)
       CALL pao_g_ang_did_change(angles, qs_env=gf_env%qs_env,&
            global_env=g_env,&
            angles_changed=.TRUE., s_changed=.FALSE.,&
            h_changed=.FALSE., full_reset=.FALSE., error=error)
       CALL pao_g_ang_get_min_m(angles, min_hamiltonian_m=h_matrixes,&
            min_overlap_m=s_matrixes, qs_env=gf_env%qs_env, global_env=g_env,&
            error=error)
       h_matrix => h_matrixes(1)%matrix
       s_matrix => s_matrixes(1)%matrix
       IF (cp_debug .AND. debug_this_module) THEN
          CALL cp_sm_output(logger=logger, &
               outputName="obj_f_min_h", &
               fromWhere=routineP, iter=gf_env%n_eval, matrix=h_matrix%sm,&
               para_env=g_env%para_env,error=error)
          CALL cp_sm_output(logger=logger, &
               outputName="obj_f_min_s", &
               fromWhere=routineP, iter=gf_env%n_eval, matrix=s_matrix%sm,&
               para_env=g_env%para_env,error=error)
          CALL cp_sm_output(logger=logger, &
               outputName="obj_f_nui", &
               fromWhere=routineP, iter=gf_env%n_eval, &
               matrix=pao_env%angles_att%NUi_injection%sm,&
               para_env=g_env%para_env,error=error)
       END IF


       ! calculates the new density matrix with diagonalization
       ! do not do it every time?
       CALL copy_sm_to_fm(real_matrix=cp_get_matrix(s_matrix),&
            fm=gf_env%min_blacs_1)
       CALL calculate_ortho_matrix(ortho=gf_env%min_blacs_1,&
            work=gf_env%min_blacs_2,ndep=size_null_space, globenv=g_env)
       CALL copy_sm_to_fm(real_matrix=cp_get_matrix(h_matrix),&
            fm=gf_env%min_blacs_2)
       CALL orthogonalise_matrix(ortho=gf_env%min_blacs_1,&
            ks=gf_env%min_blacs_2,&
            work=gf_env%min_blacs_3,&
            globenv=g_env)
       
       ! min_blacs1 = S^{-1/2}, min_blacs2= H_{KS}
       CALL eigensolver(ks=gf_env%min_blacs_2,&
            mo_set=gf_env%mo_set(1)%mo_set,&
            ortho=gf_env%min_blacs_1, work=gf_env%min_blacs_3,&
            globenv=g_env)

       CALL calculate_density_matrix(mo_set=gf_env%mo_set(1)%mo_set,&
            density_matrix=gf_env%min_blacs_2)
       
       CALL copy_fm_to_sm(gf_env%min_blacs_2,&
            min_density_m)

       CALL symmetrise_diagonal_blocks(min_density_m)

       CALL get_mo_set(mo_set=gf_env%mo_set(1)%mo_set,&
                       eigenvalues=mo_eigenvalues,&
                       occupation_numbers=mo_occupation_numbers)

       IF (PRESENT(f)) THEN
          f=DOT_PRODUCT(mo_occupation_numbers,mo_eigenvalues)/2.0_dp
          IF (cp_debug .AND. debug_this_module) THEN
             CALL cp_sm_output(logger=logger, &
                  outputName="obj_f_min_p", &
                  fromWhere=routineP, iter=gf_env%n_eval,&
                  matrix=min_density_m,&
                  para_env=g_env%para_env,error=error)
             CALL cp_write_output(logger=logger,&
                  fromWhere=routineP,&
                  outputName="obj_f_f", iter=gf_env%n_eval,&
                  message="f ="//cp_to_string(f), error=error)
          END IF
       END IF

       IF (PRESENT(gradient)) THEN
          ! uses the filtred P matrix
          ! should change it soon
          
          diagMatrix => cp_get_matrix(gf_env%diag_matrix) ! for xlf...
          ! ** gf_env%diag_matrix = df/du_{i,j}
                    
          ! 2 H B P i^*
          CALL copy_sm_to_fm(full_h_rmpv(1)%matrix,&
               gf_env%scf_work1)
          CALL cp_m2f_with_B_M_pi(min_density_m,cp_get_matrix(gf_env%full_work_m),&
               glob_angles=angles, qs_env=gf_env%qs_env, global_env=g_env,&
               error=error)
          CALL copy_sm_to_fm(cp_get_matrix(gf_env%full_work_m),&
               gf_env%scf_work2)
          CALL cp_fm_symm(side="L",uplo="U",m=full_bas_size,n=full_bas_size,&
               alpha=1.0_dp,matrix_a=gf_env%scf_work1,&
               matrix_b=gf_env%scf_work2,&
               beta=0.0_dp,matrix_c=gf_env%scf_work3)
          CALL cp_set_to(diagMatrix,value=0.0_dp)
          CALL add_blacs_to_block_diag_sm(gf_env%scf_work3,&
               diagMatrix, &
               d_struct=gf_env%d_struct, error=error)
          IF (cp_debug .AND. debug_this_module) THEN
             CALL cp_sm_output(logger=logger, &
               outputName="obj_f_BPpi", &
               fromWhere=routineP, iter=gf_env%n_eval, &
               matrix=cp_get_matrix(gf_env%full_work_m),&
               para_env=g_env%para_env,error=error)
             CALL cp_sm_output(logger=logger, &
               outputName="obj_f_2HBPpi", &
               fromWhere=routineP, iter=gf_env%n_eval, matrix=diagMatrix,&
               para_env=g_env%para_env,error=error)
          END IF
          ! 2 S B P minH P i^*
          CALL copy_sm_to_fm(h_matrix%sm,&
               gf_env%min_blacs_1)
          CALL cp_fm_upper_to_full(gf_env%min_blacs_1,gf_env%min_blacs_2)
          CALL copy_sm_to_fm(min_density_m,&
               gf_env%min_blacs_2)
          CALL cp_fm_symm('R','U',min_bas_size,min_bas_size,&
               1.0_dp,gf_env%min_blacs_2,gf_env%min_blacs_1, 0.0_dp,&
               gf_env%min_blacs_3)
          CALL cp_fm_symm('L','U',min_bas_size,min_bas_size,&
               -0.5_dp,gf_env%min_blacs_2,gf_env%min_blacs_3, 0.0_dp,&
               gf_env%min_blacs_1)
          CALL copy_fm_to_sm(gf_env%min_blacs_1,&
               cp_get_matrix(gf_env%min_work_m))
          CALL cp_m2f_with_B_M_pi(cp_get_matrix(gf_env%min_work_m),&
               cp_get_matrix(gf_env%full_work_m),glob_angles=angles,&
               qs_env=gf_env%qs_env, global_env=g_env,&
               error=error)
          CALL copy_sm_to_fm(cp_get_matrix(gf_env%full_work_m),&
               gf_env%scf_work2)
          CALL copy_sm_to_fm(full_s_rmpv(1)%matrix,&
               gf_env%scf_work1)
          CALL cp_fm_symm('L','U',full_bas_size,full_bas_size,&
               1.0_dp,gf_env%scf_work1,gf_env%scf_work2,0.0_dp,&
               gf_env%scf_work3)
          CALL add_blacs_to_block_diag_sm(gf_env%scf_work3,&
               diagMatrix,&
               d_struct=gf_env%d_struct, error=error)
          ! diag_matrix = N diag_matrix = 2 N H B P i^* - 2 N S B P minH P i^*
          CALL cp_l_multiply_atomic_ortho(diagMatrix,qs_env=gf_env%qs_env,&
               global_env=g_env,error=error)
          
          IF (cp_debug .AND. debug_this_module) THEN
             CALL cp_sm_output(logger=logger, &
                  outputName="obj_f_df_du", &
                  fromWhere=routineP, iter=gf_env%n_eval, matrix=diagMatrix,&
                  para_env=g_env%para_env,error=error)
          END IF
          
          ! ** gradient = df/du_{i,j} du_{i,j}/dtheta_k
          CPAssert(para_env%num_pe==1,cp_failure_level,routineP,error,failure)
          ! not thought in parallel
          gradient=0.0_dp
          CALL cp_init(block_iterator,matrix=cp_get_matrix(gf_env%diag_matrix))
          CALL pao_l_ang_init(local_angles, qs_env=gf_env%qs_env, &
               global_env=g_env, atom_nr=1,&
               glob_angles=angles, error=error)
          DO
             IF (.NOT.cp_next(block_iterator)) EXIT
             CALL cp_get(block_iterator, block_row=block_row,&
                  block_col=block_col, block_val=block_val)
             IF (block_col==block_row) THEN
                CALL cp_set(local_angles, atom_nr=block_col)
                CALL cp_init(calc_u_du,local_angles)
                ALLOCATE(du(SIZE(block_val,1),SIZE(block_val,2)),stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                IF (.NOT.failure) THEN
                   DO
                      IF (.NOT.cp_next(calc_u_du)) EXIT
                      CALL pao_calc_u_du_get(calc_u_du, du=du, &
                           min_bas=min_bas, excl_bas=excl_bas)
                      gradient(cp_get_g_index(local_angles,&
                           min_bas=min_bas, excl_bas=excl_bas))=&
                           SUM(du*block_val)
                   END DO
                   CALL cp_dealloc_ref(calc_u_du)
                   DEALLOCATE(du,stat=stat)
                   CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
                END IF
             END IF
          END DO
          CALL pao_l_ang_dealloc_ref(local_angles,error=error)
          CALL cp_sm_b_i_dealloc_ref(block_iterator,error=error)
          CALL mp_sum(gradient,g_env%group)
          IF (cp_debug.AND.debug_this_module) THEN
             CALL cp_1d_r_output(logger,"obj_f_grad",&
                  routineP,gf_env%n_eval,gradient,error=error)
          END IF
       END IF
    END IF
  END SUBROUTINE pao_gerd_f_eval_at
!***************************************************************************

!=========== qs_scf rip off (to avoid circular ref) to do =========

  SUBROUTINE calculate_ortho_matrix(ortho,work,ndep,globenv)

!   Purpose: Calculate the orthogonalization matrix (S**(-1/2))

!   History: - Creation (01.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(cp_fm_type), POINTER                :: ortho, work
    INTEGER, INTENT(OUT)                     :: ndep
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    REAL(KIND=dp), PARAMETER                 :: eps_eigval = 1.0E-5_dp, &
                                                work_syevx = 1.0_dp

    INTEGER                                  :: handle

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_ortho_matrix","I","",handle)

    CALL cp_fm_power(ortho,work,-0.5_dp,eps_eigval,ndep,work_syevx)
    CALL cp_fm_upper_to_full(ortho,work)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE calculate_ortho_matrix


  SUBROUTINE orthogonalise_matrix(ortho,ks,work,globenv)

!   Purpose: Orthogonalisation matrix (S**(-1/2))

!   History: - Creation (01.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(cp_fm_type), POINTER                :: ortho, ks, work
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    INTEGER                                  :: handle, nao

!   ---------------------------------------------------------------------------

    CALL timeset("orthogonalise_matrix","I","",handle)

    CALL cp_fm_get_info(matrix=ks,nrow_global=nao)
    CALL cp_fm_symm("L","U",nao,nao,1.0_dp,ks,ortho,0.0_dp,work)
    CALL cp_fm_gemm("T","N",nao,nao,nao,1.0_dp,ortho,work,0.0_dp,ks)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE orthogonalise_matrix

! *****************************************************************************

  SUBROUTINE eigensolver(ks,mo_set,ortho,work,globenv)

!   Purpose: Diagonalise the Kohn-Sham matrix to get a new set of MO eigen-
!            vectors and MO eigenvalues.

!   History: - Creation (01.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(cp_fm_type), POINTER                :: ks
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_fm_type), POINTER                :: ortho, work
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(len=*), PARAMETER :: routineN = "eigensolver", &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nao, nmo
    INTEGER, SAVE                            :: iter = 0
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mo_eigenvalues
    TYPE(cp_fm_type), POINTER                :: mo_eigenvectors

!   ---------------------------------------------------------------------------

    CALL timeset("eigensolver","I","",handle)

    CALL get_mo_set(mo_set=mo_set,&
                    nao=nao,&
                    nmo=nmo,&
                    eigenvalues=mo_eigenvalues,&
                    eigenvectors=mo_eigenvectors)

!   *** Diagonalise the Kohn-Sham matrix ***

    CALL cp_fm_syevx(ks,work,mo_eigenvalues,nmo,1.0_dp)

    CALL cp_fm_gemm("N","N",nao,nmo,nao,1.0_dp,ortho,work,0.0_dp,&
                    mo_eigenvectors)

    IF (cp_debug .AND. debug_this_module) THEN
       iter=iter+1
       CALL cp_1d_r_output(logger=cp_get_default_logger(),&
            outputName="obj_f_mo_eval", &
            fromWhere=routineP, iter=iter, array=mo_eigenvalues)
!       call cp_2d_r_output(logger=cp_get_default_logger(),&
!            outputName="obj_f_mo_evect", &
!            fromWhere=routineP, iter=iter, array=mo_eigenvectors)
    END IF
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE eigensolver

!!****f* pao_obj_function/cp_l_compose_atomic_ortho [1.0] *
!!
!!   NAME
!!     cp_l_compose_atomic_ortho
!!
!!   FUNCTION
!!     left multiplies with the atomic orthogonalization
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - matrix: the matrix to be multiplied, on output contains
!!       N matrix (where N is the atomic orthogonalization)
!!     - qs_env: the qs_environment that defines the atomic ortho
!!     - global_env: global environement (to communicate)
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
SUBROUTINE cp_l_multiply_atomic_ortho(matrix,qs_env,global_env,error)
    TYPE(real_matrix_type), POINTER          :: matrix
    TYPE(qs_environment_type), INTENT(in), &
      TARGET                                 :: qs_env
    TYPE(global_environment_type), &
      INTENT(in), TARGET                     :: global_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_l_multiply_atomic_ortho', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: block_row, max_full_bas, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: t1_m
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a_ortho, block_val
    TYPE(cp_matrix_block_iterator)           :: iterator

! pao_env might be modified

  failure=.FALSE.
  NULLIFY(block_val,a_ortho)
  
  CPPrecondition(ASSOCIATED(matrix),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL pao_qs_env_get(qs_env,global_env=global_env,&
          pao_max_full_bas=max_full_bas,error=error)
     ALLOCATE(t1_m(max_full_bas,max_full_bas),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF (.NOT.failure) THEN
     CALL cp_sm_b_i_init(iterator,matrix=matrix,error=error)
     DO
        IF (.NOT. cp_next(iterator,block_row=block_row,&
             block_val=block_val,error=error)) EXIT
        a_ortho => pao_qs_env_get_atomic_ortho(qs_env=qs_env,&
             global_env=global_env,atom_nr=block_row,error=error)
        CPPrecondition(ASSOCIATED(block_val),cp_failure_level,routineP,error,failure)
        CPPrecondition(ASSOCIATED(a_ortho),cp_failure_level,routineP,error,failure)
        IF (.NOT.failure) THEN
           CALL cp_assert(SIZE(block_val,1)==SIZE(a_ortho,2).AND.&
                SIZE(a_ortho,1)==SIZE(a_ortho,2),cp_failure_level,&
                cp_assertion_failed,routineP,"size mismach in "//&
                CPSourceFileRef,&
                error=error,failure=failure)
           CALL cp_assert(SIZE(t1_m,1)>=SIZE(a_ortho,1).AND.&
                SIZE(t1_m,2)>=SIZE(block_val,2),cp_failure_level,&
                cp_assertion_failed,routineP,"t1_m too small in "//&
                CPSourceFileRef,&
                error=error,failure=failure)           
        END IF
        IF (.NOT.failure) THEN
           CALL dgemm('N','N',SIZE(a_ortho,1),SIZE(block_val,2),SIZE(a_ortho,2),&
                1.0_dp,a_ortho,SIZE(a_ortho,1),block_val,SIZE(block_val,1),&
                0.0_dp,t1_m,SIZE(t1_m,1))
           block_val=t1_m(1:SIZE(block_val,1),1:SIZE(block_val,2))
        END IF
     END DO
     CALL cp_sm_b_i_dealloc_ref(iterator,error=error)
     DEALLOCATE(t1_m,stat=stat)
     CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
  END IF
END SUBROUTINE cp_l_multiply_atomic_ortho
!***************************************************************************

!!****f* blacs/add_blacs_to_block_diag_sm [1.0] *
!!
!!   NAME
!!     add_blacs_to_block_diag_sm
!!
!!   SYNOPSIS
!!     Subroutine add_blacs_to_block_diag_sm(source_m, target_m, para_env,&
!!         d_struct, error)
!!       Type(cp_fm_type), Pointer:: source_m
!!       Type(real_matrix_type), Intent (INOUT), Target:: target_m
!!       Type(cp_para_env_type), Pointer:: para_env
!!       Integer, Dimension(:), Pointer, Optional:: d_struct
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine add_blacs_to_block_diag_sm
!!
!!   FUNCTION
!!     add the diagonal blocs of a blacs matrix to a sparse matrix
!!
!!   NOTES
!!     - moved here from qs_blacs, presumably needs a serious rework (parallel)
!!
!!   INPUTS
!!     source_m: the blacs matrix whose diagonal blocks are copied
!!     target_m: the sparse matrix that will contain the diagonal
!!               blocks (the blocks should already be allocated)
!!     d_struct: the distribution of the the blocks among the processors:
!!               the number of the processor that has each block.
!!               If an associated pointer is given its contents should
!!               be rightly initialized, if not associated then
!!               after the call it will contain the actual distribution
!!               and you are responsible of deallocating it.
!!     error: variable to control error logging, stopping,...
!!            see module cp_error_handling
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     JVDV : see notes
!!
!!*** **********************************************************************
  SUBROUTINE add_blacs_to_block_diag_sm(source_m, target_m, &
       d_struct, error)
    TYPE(cp_fm_type), POINTER                :: source_m
    TYPE(real_matrix_type), INTENT(inout), &
      TARGET                                 :: target_m
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: d_struct
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'add_blacs_to_block_diag_sm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, nblock_col, &
                                                ncols, nrows, sm_block_nr, &
                                                stat
    INTEGER, DIMENSION(:), POINTER           :: first_col, first_row, &
                                                last_col, last_row, &
                                                my_d_struct
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block_val
    TYPE(cp_matrix_block_iterator)           :: iterator
    TYPE(real_matrix_type), POINTER          :: target_ptr

    failure=.FALSE.
    target_ptr => target_m
    NULLIFY(my_d_struct)

    CALL timeset(routineN//','//moduleN,'I',"",handle)
    CPPrecondition(ASSOCIATED(target_ptr),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CALL get_matrix_info(target_ptr, nrow=nrows, ncol=ncols,&
            nblock_col=nblock_col)
       CPAssert(ncols==nrows,cp_warning_level,routineP,error,failure)
       ncols=MIN(nrows,ncols)
       IF (PRESENT(d_struct)) my_d_struct => d_struct
       IF (.NOT.ASSOCIATED(my_d_struct)) THEN
! build distribution structure
! assume every diagonal block is on one (and just one) proc
          ALLOCATE(my_d_struct(ncols),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          my_d_struct=0
          CALL cp_init(iterator, matrix=target_m, error=error)
          DO
             IF (.NOT.cp_next(iterator)) EXIT
             CALL cp_get(iterator,block_row=i,block_col=j,error=error)
             IF (i==j) THEN
                my_d_struct(i)=source_m%matrix_struct%para_env%mepos
             ELSE
                CALL cp_error_message(cp_warning_level,&
                     routineP,"there is an off diagonal block! in "//&
                     CPSourceFileRef,error)
! just as info (only on the local proc)
             END IF
          END DO
          CALL cp_dealloc_ref(iterator,error=error)
          CALL mp_sum(my_d_struct,source_m%matrix_struct%para_env%group)
       END IF
    END IF
    CPAssert(ASSOCIATED(my_d_struct),cp_failure_level,routineP,error,failure)
    CPAssert(SIZE(my_d_struct)==nblock_col,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPAssert(ALL(my_d_struct<source_m%matrix_struct%para_env%num_pe),cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT.failure) THEN
!group = global_env%group
!source = global_env%source

       CALL get_matrix_info(matrix=target_ptr,&
            first_row=first_row,&
            first_col=first_col,&
            last_row=last_row,&
            last_col=last_col)
#if (__SCALAPACK)
!MK       call blacs_pinfo(blacs_proc_nr, n_blacs_proc)
!MK       CPPrecondition(blacs_proc_nr==global_env%mepos,cp_failure_level,routineP,error,failure)
!MK       call blacs_gridinfo(source_m%context,nprow,npcol,myprow,mypcol)
!MK       start_proc_row=source_m%matrix_struct%descriptor(7)
!MK       start_proc_col=source_m%matrix_struct%descriptor(8)
!MK!    if ( (myprow >= n_proc_row) .or. (mypcol >= n_proc_col) ) then
!MK!       has_b_matrix=.false.
!MK!    else
!MK!       has_b_matrix=.true.
!MK!    end if
!MK
!MK       do sm_block_nr=1,ncols
!MK          CPPrecondition(first_row(sm_block_nr)==first_col(sm_block_nr),cp_warning_level,routineP,error,failure)
!MK! the processor that has the sm_block
!MK          sm_block_proc=d_struct(sm_block_nr)
!MK
!MK! 0 based indexing
!MK          b_block_row_start= (first_row(sm_block_nr)-1)/source_m%nrow_block
!MK          b_block_row_stop= (last_row(sm_block_nr)-1)/source_m%nrow_block
!MK          b_block_col_start= (first_col(sm_block_nr)-1)/source_m%ncol_block
!MK          b_block_col_stop= (last_col(sm_block_nr)-1)/source_m%ncol_block
!MK
!MK          do b_block_row= b_block_row_start,b_block_row_stop
!MK             do b_block_col= b_block_col_start,b_block_col_stop
!MK                blacs_block_proc_row= mod(b_block_row+start_proc_row,nprow)
!MK                blacs_block_proc_col= mod(b_block_col+start_proc_col,npcol)
!MK                blacs_proc_nr= blacs_pnum(source_matrix%context,&
!MK                     blacs_block_proc_row, blacs_block_proc_col)
!MK! the blacs proc nr should be just the mpi nr
!MK                if (blacs_block_proc_nr == sm_block_nr) then
!MK                   if (blacs_block_proc_nr == global_env%mepos) then
!MK                      start_row_of_b= first_col(sm_block_nr)-b_block_row*&
!MK                           source_m%nrow_block
!MK                      end_row_of_b= min(source_m%nrow_block,&
!MK                           last_col(sm_block_nr)-b_block_row*&
!MK                           source_m%nrow_block)
!MK                      start_row_of_sm=max(0,-start_row_of_sm)
!MK                      start_row_of_b=max(0,start_row_of_sm)
!MK                      end_row_of_sm=start_row_of_sm+end_row_of_b-&
!MK                           start_row_of_b
!MK! to do
!MK!source_m%p(blacs_block_proc_row,blacs_block_proc_col)&
!MK!     %block(start_row_of_b:end_row_of_b,&
!MK!     start_col_of_b,end_col_of_b)
!MK
!MK                   end if
!MK                end if
!MK             end do
!MK          end do
!MK       end do
!MK
#else
       DO sm_block_nr=1,nblock_col
          CALL get_block_node(target_ptr,sm_block_nr,sm_block_nr,&
               BLOCK=block_val)
          block_val=block_val+source_m%local_data &
               (first_row(sm_block_nr):last_row(sm_block_nr),&
               first_col(sm_block_nr):last_col(sm_block_nr))
       END DO
#endif
    END IF
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE add_blacs_to_block_diag_sm
!***************************************************************************
END MODULE pao_obj_function

