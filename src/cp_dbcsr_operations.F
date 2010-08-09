!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   DBCSR operations in CP2K
!> \author  Urban Borstnik
!> \date    2009-05-12
!> \version 0.8
!>
!> <b>Modification history:</b>
!> - Created 2009-05-12
! *****************************************************************************
MODULE cp_dbcsr_operations
  USE array_types,                     ONLY: array_data,&
                                             array_i1d_obj,&
                                             array_new,&
                                             array_nullify,&
                                             array_release
  USE cp_blacs_env,                    ONLY: cp_blacs_env_type,&
                                             cp_blacs_env_write,&
                                             get_blacs_info
  USE cp_cfm_types,                    ONLY: cp_cfm_type
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_checksum, cp_dbcsr_copy, cp_dbcsr_create, &
       cp_dbcsr_distribute, cp_dbcsr_filter, cp_dbcsr_finalize, &
       cp_dbcsr_init, cp_dbcsr_multiply, cp_dbcsr_norm, cp_dbcsr_release, &
       cp_dbcsr_replicate_all, cp_dbcsr_scale, cp_dbcsr_sum_replicated, &
       cp_dbcsr_trace, cp_dbcsr_transposed, cp_dbcsr_verify_matrix, &
       cp_dbcsr_work_create
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_obj_p_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_frobenius_norm,&
                                             cp_fm_gemm
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_release,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply,&
                                             cp_sm_plus_fm_fm_t
  USE dbcsr_block_access,              ONLY: dbcsr_access_flush,&
                                             dbcsr_get_block,&
                                             dbcsr_get_block_p,&
                                             dbcsr_get_tmp_block2d,&
                                             dbcsr_put_block,&
                                             dbcsr_remove_block,&
                                             dbcsr_reserve_block2d
  USE dbcsr_block_operations,          ONLY: block_set
  USE dbcsr_data_methods,              ONLY: dbcsr_get_data,&
                                             dbcsr_get_data_p
  USE dbcsr_dist_operations
  USE dbcsr_error_handling,            ONLY: dbcsr_error_set,&
                                             dbcsr_error_stop,&
                                             dbcsr_error_type
  USE dbcsr_io
  USE dbcsr_iterator_operations,       ONLY: dbcsr_iterator_blocks_left,&
                                             dbcsr_iterator_next_block,&
                                             dbcsr_iterator_start,&
                                             dbcsr_iterator_stop
  USE dbcsr_methods,                   ONLY: &
       dbcsr_col_block_sizes, dbcsr_distribution, &
       dbcsr_distribution_col_dist, dbcsr_distribution_init, &
       dbcsr_distribution_local_cols, dbcsr_distribution_local_rows, &
       dbcsr_distribution_mp, dbcsr_distribution_ncols, &
       dbcsr_distribution_new, dbcsr_distribution_nlocal_cols, &
       dbcsr_distribution_nlocal_rows, dbcsr_distribution_nrows, &
       dbcsr_distribution_release, dbcsr_distribution_row_dist, &
       dbcsr_get_data_size, dbcsr_get_data_type, dbcsr_get_info, &
       dbcsr_get_matrix_type, dbcsr_mp_get_coordinates, dbcsr_mp_group, &
       dbcsr_mp_hold, dbcsr_mp_mynode, dbcsr_mp_mypcol, dbcsr_mp_myprow, &
       dbcsr_mp_new, dbcsr_mp_npcols, dbcsr_mp_nprows, dbcsr_mp_numnodes, &
       dbcsr_mp_pgrid, dbcsr_mp_release, dbcsr_nblkcols_total, &
       dbcsr_nblkrows_total, dbcsr_nfullcols_total, dbcsr_nfullrows_total, &
       dbcsr_row_block_sizes, dbcsr_valid_index, dbcsr_wm_use_mutable
  USE dbcsr_operations,                ONLY: dbcsr_replace_blocks
  USE dbcsr_transformations,           ONLY: dbcsr_complete_redistribute
  USE dbcsr_types,                     ONLY: &
       dbcsr_distribution_obj, dbcsr_iterator, dbcsr_mp_obj, &
       dbcsr_norm_frobenius, dbcsr_norm_gershgorin, dbcsr_norm_maxabsnorm, &
       dbcsr_obj, dbcsr_type_antisymmetric, dbcsr_type_complex_8, &
       dbcsr_type_normal, dbcsr_type_real_8, dbcsr_type_symmetric
  USE dbcsr_util,                      ONLY: convert_sizes_to_offsets,&
                                             dbcsr_set_debug,&
                                             swap
  USE dbcsr_work_operations,           ONLY: add_work_coordinate
  USE distribution_2d_types,           ONLY: distribution_2d_get,&
                                             distribution_2d_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             int_8,&
                                             sp
  USE machine,                         ONLY: m_walltime
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_sum
  USE qs_neighbor_list_types,          ONLY: &
       first_list, first_node, get_neighbor_list, get_neighbor_list_set, &
       get_neighbor_node, neighbor_list_set_p_type, neighbor_list_type, &
       neighbor_node_type, next
  USE sparse_matrix_types,             ONLY: &
       add_block_node, allocate_matrix, checksum_matrix, deallocate_matrix, &
       first_block_node, get_block_node, next_block_node, &
       real_block_node_type, real_matrix_type, replicate_matrix
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  !$ USE OMP_LIB
#include "cp_common_uses.h"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dbcsr_operations'

  PUBLIC :: cp_dbcsr_check, cp_dbcsr_multiply_local
  PUBLIC :: cp_dbcsr_multiply_vec, cp_dbcsr_pack_vec,&
       packed_vec_scale, cp_dbcsr_unpack_vec,&
       cp_dbcsr_mult_pack_vec_local, packed_vec_bif_tech, &
       packed_vec_bif_tech2, cp_dbcsr_copy_vec, packed_vec_build_u, packed_vec_bcast, &
       packed_vec_ini

  ! CP2K API emulation
  PUBLIC :: cp_dbcsr_add_block_node,&
            cp_dbcsr_deallocate_matrix,&
            cp_dbcsr_allocate_matrix_set, cp_dbcsr_deallocate_matrix_set,&
            test_sm,&
            cp_dbcsr_from_fm, cp_dbcsr_from_sm, copy_fm_to_dbcsr, copy_dbcsr_to_fm,&
            copy_dbcsr_to_cfm, copy_cfm_to_dbcsr, &
            sm_from_dbcsr, cp_dbcsr_type_from_sm_symmetry, &
            cp_dbcsr_sm_fm_multiply, cp_dbcsr_plus_fm_fm_t,&
            cp_dbcsr_get_id_nr,&
            cp_dbcsr_alloc_block_from_nbl


  ! distribution_2d_type compatibility
  PUBLIC :: cp_dbcsr_dist2d_to_dist

  PUBLIC :: cp_dbcsr_copy_columns_hack

  INTERFACE cp_dbcsr_allocate_matrix_set
     MODULE PROCEDURE allocate_dbcsr_matrix_set, allocate_dbcsr_matrix_set_2d
  END INTERFACE

  INTERFACE cp_dbcsr_deallocate_matrix_set
     MODULE PROCEDURE deallocate_dbcsr_matrix_set,&
                      deallocate_dbcsr_matrix_set_2d
  END INTERFACE

  INTERFACE cp_dbcsr_plus_fm_fm_t
     MODULE PROCEDURE cp_dbcsr_plus_fm_fm_t_native
  END INTERFACE

  PRIVATE

  INTEGER, SAVE, PRIVATE :: last_matrix_id=0

#define DBG IF (dbg) WRITE(*,*)routineN//" ",
#define DBGV IF (bcsr_verbose) WRITE(*,*)routineN//" ",
#define DBGI IF (info) WRITE(*,*)routineN//" ",
#define DEBUG_HEADER        LOGICAL :: dbg, info
#define DEBUG_BODY        dbg = .FALSE. ; CALL dbcsr_set_debug(dbg, info=info)
#define DEBUG_BODYY       dbg = .TRUE. ; CALL dbcsr_set_debug(dbg, info=info)

CONTAINS

! *****************************************************************************
!> \brief hack for dbcsr_copy_columns
!> \author vw
! *****************************************************************************
  SUBROUTINE cp_dbcsr_copy_columns_hack(matrix_b, matrix_a,&
       ncol, source_start, target_start, para_env, blacs_env, error)

    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_b
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a
    INTEGER, INTENT(IN)                      :: ncol, source_start, &
                                                target_start
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_copy_columns_hack', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nfullcols_total, &
                                                nfullrows_total
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct
    TYPE(cp_fm_type), POINTER                :: fm_matrix_a, fm_matrix_b

    NULLIFY(fm_matrix_a, fm_matrix_b, fm_struct)
    CALL dbcsr_get_info(matrix_a,nfullrows_total=nfullrows_total,nfullcols_total=nfullcols_total)
    CALL cp_fm_struct_create(fm_struct,context=blacs_env,nrow_global=nfullrows_total,&
         ncol_global=nfullcols_total,para_env=para_env,error=error)
    CALL cp_fm_create(fm_matrix_a,fm_struct,name="fm_matrix_a",error=error)
    CALL cp_fm_struct_release(fm_struct,error=error)

    CALL dbcsr_get_info(matrix_b,nfullrows_total=nfullrows_total,nfullcols_total=nfullcols_total)
    CALL cp_fm_struct_create(fm_struct,context=blacs_env,nrow_global=nfullrows_total,&
         ncol_global=nfullcols_total,para_env=para_env,error=error)
    CALL cp_fm_create(fm_matrix_b,fm_struct,name="fm_matrix_b",error=error)
    CALL cp_fm_struct_release(fm_struct,error=error)

    CALL copy_dbcsr_to_fm(matrix_a, fm_matrix_a, error=error)
    CALL copy_dbcsr_to_fm(matrix_b, fm_matrix_b, error=error)

    CALL cp_fm_to_fm(fm_matrix_a, fm_matrix_b, ncol, source_start, target_start)

    CALL copy_fm_to_dbcsr(fm_matrix_b, matrix_b, error=error)

    CALL cp_fm_release(fm_matrix_a, error=error)
    CALL cp_fm_release(fm_matrix_b, error=error)

  END SUBROUTINE cp_dbcsr_copy_columns_hack


! *****************************************************************************
!> \brief Creates a DBCSR distribution from a distribution_2d
!> \param[in] dist2d          distribution_2d
!> \param[out] dist           DBCSR distribution
!> \param[in,out] error       cp2k error
!> \par History
!>    move form dbcsr_operation 01.2010
! *****************************************************************************
  SUBROUTINE cp_dbcsr_dist2d_to_dist(dist2d, dist, error, mp_obj)
    TYPE(distribution_2d_type), INTENT(IN), &
      TARGET                                 :: dist2d
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT)                            :: dist
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    TYPE(dbcsr_mp_obj), INTENT(IN), OPTIONAL :: mp_obj

    INTEGER                                  :: mypcol, myproc, myprow, &
                                                numproc
    INTEGER, DIMENSION(:), POINTER           :: col_dist_data, row_dist_data
    INTEGER, DIMENSION(:, :), POINTER        :: pgrid
    TYPE(array_i1d_obj)                      :: cd, rd
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dbcsr_mp_obj)                       :: mp_env
    TYPE(distribution_2d_type), POINTER      :: dist2d_p

!
!   ---------------------------------------------------------------------------

    dist2d_p => dist2d
    CALL distribution_2d_get(dist2d_p, error=error,&
         row_distribution=row_dist_data, col_distribution=col_dist_data,&
         blacs_env=blacs_env)
    CALL get_blacs_info(blacs_env, para_env=para_env,&
         my_process_row=myprow, my_process_column=mypcol,&
         blacs2mpi=pgrid)
    myproc = para_env%mepos
    numproc = para_env%num_pe
    IF (PRESENT (mp_obj)) THEN
       mp_env = mp_obj
       CALL dbcsr_mp_hold (mp_env)
    ELSE
       CALL dbcsr_mp_new(mp_env, pgrid, para_env%group, myproc, numproc,&
            myprow, mypcol)
    ENDIF
    CALL array_nullify (rd)
    CALL array_nullify (cd)
    CALL array_new(rd, row_dist_data)
    CALL array_new(cd, col_dist_data)
    CALL dbcsr_distribution_new(dist, mp_env, rd, cd)
    CALL dbcsr_mp_release (mp_env)
    CALL array_release (rd)
    CALL array_release (cd)
  END SUBROUTINE cp_dbcsr_dist2d_to_dist

! *****************************************************************************
!> \brief allocate the blocks of a dbcsr based on the neighbor list
!> \note
!> \param matrix        the matrix
!> \param sab_orb       the corresponding neighbor list
!> \param error
!> \par History
!>      11.2009 created vw
!> \author vw
! *****************************************************************************

  SUBROUTINE cp_dbcsr_alloc_block_from_nbl(matrix,sab_orb,error)

    TYPE(dbcsr_obj), POINTER                 :: matrix
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'cp_dbcsr_alloc_block_from_nbl', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=1)                         :: symmetry
    INTEGER                                  :: iab, iatom, icol, ilist, &
                                                inode, irow, jatom, &
                                                last_jatom, nlist, nnode
    LOGICAL                                  :: failure, new_atom_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block
    TYPE(neighbor_list_type), POINTER        :: sab_orb_neighbor_list, &
                                                sab_orb_neighbor_list_local
    TYPE(neighbor_node_type), POINTER        :: sab_orb_neighbor_node

    failure = .FALSE.

    symmetry = dbcsr_get_matrix_type(matrix)

    CPPrecondition(ASSOCIATED(sab_orb),cp_failure_level,routineP,error,failure)

    !nkind = SIZE(atomic_kind_set,1)
    !DO ikind=1,nkind
    !DO jkind=1,nkind
    DO iab = 1, SIZE(sab_orb)
       !iab = ikind + nkind*(jkind - 1)

       IF (.NOT.ASSOCIATED(sab_orb(iab)%neighbor_list_set)) CYCLE

       CALL get_neighbor_list_set(neighbor_list_set=&
            sab_orb(iab)%neighbor_list_set,&
            nlist=nlist)

       NULLIFY ( sab_orb_neighbor_list )

       DO ilist=1,nlist

          IF ( .NOT. ASSOCIATED(sab_orb_neighbor_list) ) THEN
             sab_orb_neighbor_list => first_list(sab_orb(iab)%neighbor_list_set)
          ELSE
             sab_orb_neighbor_list => next(sab_orb_neighbor_list)
          END IF
          sab_orb_neighbor_list_local => sab_orb_neighbor_list

          CALL get_neighbor_list(neighbor_list=sab_orb_neighbor_list_local,&
               atom=iatom,&
               nnode=nnode)

          last_jatom = 0

          sab_orb_neighbor_node => first_node(sab_orb_neighbor_list_local)

          DO inode=1,nnode

             CALL get_neighbor_node(neighbor_node=sab_orb_neighbor_node,&
                  neighbor=jatom)

             IF (jatom /= last_jatom) THEN
                new_atom_b = .TRUE.
                last_jatom = jatom
             ELSE
                new_atom_b = .FALSE.
                sab_orb_neighbor_node => next(sab_orb_neighbor_node)
                CYCLE
             END IF
             NULLIFY(block)
             IF(symmetry==dbcsr_type_normal) THEN
                CALL cp_dbcsr_add_block_node(matrix=matrix,&
                     block_row=iatom,&
                     block_col=jatom,&
                     block=block,error=error)
             ELSE
                IF(iatom<=jatom) THEN
                   irow = iatom
                   icol = jatom
                ELSE
                   irow = jatom
                   icol = iatom
                END IF
                CALL cp_dbcsr_add_block_node(matrix=matrix,&
                     block_row=irow,&
                     block_col=icol,&
                     block=block,error=error)
             END IF

             block(:,:) = 0.0_dp

             sab_orb_neighbor_node => next(sab_orb_neighbor_node)

          END DO  ! inode

       END DO  ! ilist

       !END DO  ! jkind

    END DO  ! ikind

    CALL cp_dbcsr_finalize( matrix, error=error )

  END SUBROUTINE cp_dbcsr_alloc_block_from_nbl



! *****************************************************************************
!> \brief   Checks the new bcsr routines, to find bugs early
!> \author  Urban
!> \date    2008-10-07
! *****************************************************************************
  SUBROUTINE cp_dbcsr_check(matrix_s, matrix_x, hint, error)

    TYPE(real_matrix_type), POINTER          :: matrix_s
    TYPE(cp_fm_type), POINTER                :: matrix_x
    INTEGER, INTENT(IN)                      :: hint
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_check', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk_nze, col, col_size, &
                                                colsize, handle, k, &
                                                output_unit, row, row_size, &
                                                rowsize
    INTEGER(KIND=int_8)                      :: flop
    INTEGER, DIMENSION(:), POINTER           :: cbs, rbs
    LOGICAL                                  :: failure, found, ldbg, tr
    REAL(kind=dp)                            :: checksum, newtime, oldtime, &
                                                trace, walltime1, walltime2
    REAL(KIND=dp), DIMENSION(:), POINTER     :: block1d
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block2d, block2d_2
    TYPE(array_i1d_obj)                      :: col_blk_size, &
                                                col_blk_size_right, &
                                                row_blk_size
    TYPE(cp_fm_type), POINTER                :: fmout, matrix_x2
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dbcsr_distribution_obj)             :: dist_right, product_dist, &
                                                weird_dist
    TYPE(dbcsr_error_type)                   :: dbcsr_error
    TYPE(dbcsr_iterator)                     :: iter
    TYPE(dbcsr_obj)                          :: matrix_a, matrix_b, &
                                                matrix_bt, matrix_c
    TYPE(dbcsr_obj), POINTER                 :: matrix_a2
    TYPE(distribution_2d_type), POINTER      :: dist2d
    TYPE(real_matrix_type), POINTER          :: matrix_s2

    CALL timeset(routineN,handle)

    failure=.FALSE.
    ldbg = .TRUE.
    logger => cp_error_get_logger (error)
    output_unit = cp_logger_get_default_io_unit (logger)
    !Do the new multiplication
    CALL cp_dbcsr_init (matrix_a, error)
    CALL cp_dbcsr_from_sm(matrix_a, matrix_s, error)
    CALL cp_dbcsr_verify_matrix(matrix_a, error)
!write(*,*)'-###MATRIX A'
!call dbcsr_print_dist(matrix_a%m%dist)
!call dbcsr_print(matrix_a)
    !WRITE(*,*)'S checksum=',&
    !     checksum_matrix(matrix_s, matrix_s%distribution_2d%blacs_env%para_env)
    IF (ldbg) THEN
       checksum = checksum_matrix(&
            matrix_s, matrix_s%distribution_2d%blacs_env%para_env)
       IF (ldbg.AND.(output_unit>0)) THEN
          WRITE(output_unit, '(T2,A,F12.4)')"   SM S  checksum", checksum
       ENDIF
    ENDIF
    checksum = 0.0_dp
    checksum = cp_dbcsr_checksum(matrix_a, error=error)
    !WRITE(*,*)'BCSR S checksum',checksum
    IF (ldbg.AND.(output_unit>0)) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR S  checksum", checksum
    ENDIF
    !
    ! Convert back to see if we get the same thing.
    dist2d => matrix_s%distribution_2d
    NULLIFY (matrix_s2)
    CALL sm_from_dbcsr (matrix_s2, matrix_a, dist2d, error)
    checksum = checksum_matrix(matrix_s2,&
                               matrix_s2%distribution_2d%blacs_env%para_env)
    !WRITE(*,*)'BCSR S checksum',checksum
    IF (ldbg.AND.(output_unit>0)) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"  SM2 S  checksum", checksum
    ENDIF
    ! Check that pointer reuse works.
    CALL sm_from_dbcsr (matrix_s2, matrix_a, dist2d, error)
    checksum = checksum_matrix(matrix_s2,&
                               matrix_s2%distribution_2d%blacs_env%para_env)
    !WRITE(*,*)'BCSR S checksum',checksum
    IF (ldbg.AND.(output_unit>0)) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"  SM2 S2 checksum", checksum
    ENDIF
    
    CALL deallocate_matrix (matrix_s2, error)
    !
    ! test the scale and trace routines on a square matrix
    CALL cp_dbcsr_scale(matrix_a,alpha_scalar=2.0_dp,error=error)
    checksum = 0.0_dp
    checksum = cp_dbcsr_checksum(matrix_a, error=error)
    CALL cp_dbcsr_trace(matrix_a,trace,error=error)
    IF (ldbg.AND.(output_unit>0)) THEN
       WRITE(output_unit, '(T2,2(A,F12.4))')"DBCSR 2*S  checksum", checksum, " trace(2*S)",trace
    ENDIF
    CALL cp_dbcsr_scale(matrix_a,alpha_scalar=0.5_dp,error=error)
    !
    ! test the copy and add
    ALLOCATE (matrix_a2)
    CALL cp_dbcsr_init (matrix_a2, error)
    CALL cp_dbcsr_copy(matrix_a2,matrix_a,error=error)
    checksum = 0.0_dp
    checksum = cp_dbcsr_checksum(matrix_a2, error=error)
    IF (ldbg.AND.(output_unit>0)) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR Sbis  checksum", checksum
    ENDIF
    CALL cp_dbcsr_add(matrix_a2,matrix_a,alpha_scalar=-2.0_dp,beta_scalar=2.0_dp,error=error)
    checksum = 0.0_dp
    checksum = cp_dbcsr_checksum(matrix_a2, error=error)
    IF (ldbg.AND.(output_unit>0)) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR copy and add: 2*S-2*S  checksum", checksum
    ENDIF
    !CALL cp_dbcsr_release(matrix_a2)
    !
    ! Replace blocks
    CALL dbcsr_replace_blocks(matrix_a2, matrix_a, error=dbcsr_error)
    CALL cp_dbcsr_verify_matrix(matrix_a2, error)
    checksum = cp_dbcsr_checksum(matrix_a2, error=error)
    IF (ldbg.AND.(output_unit>0)) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR replace=S checksum", checksum
    ENDIF
    CALL cp_dbcsr_release(matrix_a2, error=error)
    !
    ! filter
    CALL cp_dbcsr_copy(matrix_a2,matrix_a,error=error)
    CALL cp_dbcsr_filter(matrix_a2,1.0_dp,error=error)
    checksum = cp_dbcsr_checksum(matrix_a2, error=error)
    IF (ldbg.AND.(output_unit>0)) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR filter(S,1.0)  checksum", checksum
    ENDIF
    CALL cp_dbcsr_release(matrix_a2, error=error)
    !
    ! some norms here
    CALL cp_dbcsr_norm(matrix_a, which_norm=dbcsr_norm_maxabsnorm, &
         norm_scalar=checksum, error=error)
    IF (ldbg.AND.(output_unit>0)) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR maxabsval(S)", checksum
    ENDIF
    CALL cp_dbcsr_norm(matrix_a, which_norm=dbcsr_norm_frobenius, &
         norm_scalar=checksum, error=error)
    IF (ldbg.AND.(output_unit>0)) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR Frobenius(S)", checksum
    ENDIF
    CALL cp_dbcsr_norm(matrix_a, which_norm=dbcsr_norm_gershgorin, &
         norm_scalar=checksum, error=error)
    IF (ldbg.AND.(output_unit>0)) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR Gershgorin(S)", checksum
    ENDIF
    !
    ! Fill a matrix with all zeros.
    CALL cp_dbcsr_init (matrix_a2, error)
    CALL cp_dbcsr_from_sm(matrix_a2, matrix_s, error)
    ! Now delete all the blocks.
    !WRITE(*,*)'Deleting all blocks!'
    CALL dbcsr_iterator_start (iter, matrix_a2)
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col, k, &
            row_size=row_size, col_size=col_size)
       blk_nze = row_size * col_size
       CALL dbcsr_remove_block (matrix_a2, row, col, k)
    ENDDO
    CALL dbcsr_iterator_stop (iter)
    !WRITE(*,*)"blk_p=", matrix_a2%blk_p
    CALL cp_dbcsr_finalize (matrix_a2, error=error)
    !WRITE(*,*)"blk_p=", matrix_a2%nblks
    IF (ldbg.AND.(output_unit>0)) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR 0 checksum ", checksum
    ENDIF
    !write(*,*)'zeroing...'
    !call dbcsr_work_create(matrix_a2, dbcsr_get_num_blocks (matrix_a2),&
    !     dbcsr_nfullrows_total(matrix_a2)*dbcsr_nfullcols_total(matrix_a2))
    rbs => array_data (dbcsr_row_block_sizes (matrix_a2))
    cbs => array_data (dbcsr_col_block_sizes (matrix_a2))
    ALLOCATE (block2d(MAXVAL (rbs), MAXVAL (cbs)))
    block2d(:,:) = 0
    ALLOCATE (block1d(MAXVAL(rbs)*MAXVAL(cbs)))
    block1d(:) = 0
    k = 0
    DO row = 1, MIN(10, dbcsr_nblkrows_total (matrix_a2))
       rowsize = rbs (row)
       DO col = 1, MIN(10, dbcsr_nblkcols_total (matrix_a2))
          colsize = cbs (col)
          CALL dbcsr_put_block (matrix_a2, row, col,&
               block2d(1:rowsize, 1:colsize))
       ENDDO
    ENDDO
    CALL cp_dbcsr_finalize (matrix_a2, error=error)
    checksum = cp_dbcsr_checksum(matrix_a2, error=error)
    IF (ldbg.AND.(output_unit>0)) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR 0 checksum2", checksum
    ENDIF
    DEALLOCATE (block2d)
    DEALLOCATE (block1d)
    CALL cp_dbcsr_release (matrix_a2, error=error)
    !
    ! Tests (undefined) block addition then subsequently setting the
    ! data area.
    CALL cp_dbcsr_from_sm(matrix_a2, matrix_s, error)
    !write(*,*)'zeroing..., second time around'
    rbs => array_data (dbcsr_row_block_sizes (matrix_a2))
    cbs => array_data (dbcsr_col_block_sizes (matrix_a2))
    DO row = 1, dbcsr_nblkrows_total (matrix_a2)
       rowsize = rbs (row)
       DO col = row, dbcsr_nblkcols_total (matrix_a2)
          colsize = cbs (col)
          NULLIFY (block2d)
          !CALL dbcsr_reserve_block2d (matrix_a2, row, col,&
          !     block2d)
          CALL cp_dbcsr_add_block_node (matrix_a2, row, col, block2d, error)
!write(*,*)routineN//" setting block at virt. coor", row, col, rowsize, colsize          
          !block2d(1:rowsize, 1:colsize) = 1.0_dp
          block2d(:,:) = 1.0_dp
       ENDDO
    ENDDO
    DO row = 1, dbcsr_nblkrows_total (matrix_a2)
       DO col = 1, dbcsr_nblkcols_total (matrix_a2)
          CALL dbcsr_get_tmp_block2d (matrix_a2, row, col, block2d_2, tr, found)
          IF (found) THEN
             block2d_2(:,:) = 0.0_dp
          ENDIF
          !write(*,*)'size1:',size(block2d_2), rowsize, colsize
          !write(*,'(10(1X,F5.1))')block2d_2
       ENDDO
    ENDDO
    DO row = 1, dbcsr_nblkrows_total (matrix_a2)
       DO col = 1, dbcsr_nblkcols_total (matrix_a2)
          CALL dbcsr_get_tmp_block2d (matrix_a2, row, col, block2d, tr, found)
          !write(*,*)'size3:',size(block2d), rowsize, colsize
          !write(*,'(10(1X,F5.1))')block2d
          !IF (.NOT. found) WRITE(*,*)"Oops, missing block."
       ENDDO
    ENDDO
    CALL cp_dbcsr_finalize (matrix_a2, error=error)
    CALL cp_dbcsr_verify_matrix(matrix_a2, error)
    checksum = cp_dbcsr_checksum(matrix_a2, error=error)
    IF (ldbg.AND.(output_unit>0)) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR 0 checksum4", checksum
    ENDIF
    !
    ! Tests block getting
    DO row = 1, dbcsr_nblkrows_total (matrix_a2)
       DO col = 1, dbcsr_nblkcols_total (matrix_a2)
          CALL dbcsr_get_block_p (matrix_a2, row, col, block=block2d,&
               tr=tr,found=found)
          !WRITE(*,*)'getting block_node:',row, col, ASSOCIATED(block2d)
       ENDDO
    ENDDO
    !
    !CALL cp_dbcsr_release (matrix_a2)
    CALL cp_dbcsr_deallocate_matrix (matrix_a2, error)
    !
    checksum = 0.0_dp
    CALL cp_fm_frobenius_norm(matrix_x, checksum, error)
    !WRITE(*,*)'x checksum',checksum*checksum
    IF (ldbg.AND.(output_unit>0)) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"   FM x  checksum", checksum*checksum
    ENDIF
    CALL cp_fm_get_info(matrix_x, ncol_global=k, error=error)
    IF (MOD(hint,2) .EQ. 0) THEN
       CALL dbcsr_create_dist_r_unrot (dist_right, matrix_a%m%dist, k,&
            col_blk_size_right)
    ELSE
       CALL dbcsr_create_dist_r_rot (dist_right, matrix_a%m%dist, k,&
            col_blk_size_right)
    ENDIF
    CALL cp_dbcsr_init (matrix_b, error)
    CALL cp_dbcsr_from_fm(matrix_b, matrix_x, 0.0_dp, dist_right,&
         row_blk_size=matrix_a%m%row_blk_size,&
         col_blk_size=col_blk_size_right,&
         error=error)
    CALL cp_dbcsr_verify_matrix(matrix_b, error)
    !CALL cp_dbcsr_print (matrix_b, error)

    checksum = cp_dbcsr_checksum(matrix_b, error=error)
    !WRITE(*,*)'DBCSR x checksum',checksum
    IF (ldbg.AND.(output_unit>0)) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR x  checksum", checksum
    ENDIF
    !
    ! Test replacing a block.
    !CALL cp_dbcsr_print (matrix_b, error)
    CALL dbcsr_iterator_start (iter, matrix_b)
    IF (dbcsr_iterator_blocks_left (iter)) THEN
       checksum = cp_dbcsr_checksum(matrix_b, local=.TRUE., error=error)
       IF (ldbg.AND.(output_unit>0)) THEN
          WRITE(output_unit, '(T2,A,F12.4)')"DBCSR x checksumL", checksum
       ENDIF
       CALL dbcsr_iterator_next_block (iter, row, col, block2d, tr, &
            row_size=row_size, col_size=col_size)
       blk_nze = row_size * col_size
       !write(*,*)'got block',row,col,tr
       CALL dbcsr_put_block (matrix_b, row, col, block2d, tr)
       CALL cp_dbcsr_verify_matrix (matrix_b, error, local=.TRUE.)
       !CALL cp_dbcsr_print (matrix_b, error)
       checksum = cp_dbcsr_checksum(matrix_b, error=error, local=.TRUE.)
       IF (ldbg.AND.(output_unit>0)) THEN
          WRITE(output_unit, '(T2,A,F12.4)')"DBCSR x checksum2", checksum
       ENDIF
       ALLOCATE (block2d_2(SIZE(block2d,1),SIZE(block2d,2)))
       block2d_2(:,:) = block2d(:,:)
       !write(*,'(8(1X,F7.2))')block2d_2
       block2d(:,:) = 0.0_dp
       CALL dbcsr_put_block (matrix_b, row, col, block2d, tr)
       CALL cp_dbcsr_finalize (matrix_b, error=error)
       CALL cp_dbcsr_verify_matrix (matrix_b, error, local=.TRUE.)
       checksum = cp_dbcsr_checksum(matrix_b, error=error, local=.TRUE.)
       IF (ldbg.AND.(output_unit>0)) THEN
          WRITE(output_unit, '(T2,A,F12.4)')"DBCSR x wrng cksm", checksum
       ENDIF
       !CALL cp_dbcsr_print (matrix_b, error)
       CALL dbcsr_put_block (matrix_b, row, col, block2d_2, tr)
       CALL cp_dbcsr_finalize (matrix_b, error=error)
       CALL cp_dbcsr_verify_matrix (matrix_b, error, local=.TRUE.)
       checksum = cp_dbcsr_checksum(matrix_b, error=error, local=.TRUE.)
       IF (ldbg.AND.(output_unit>0)) THEN
          WRITE(output_unit, '(T2,A,F12.4)')"DBCSR x checksum3", checksum
       ENDIF
       !write(*,*)'removing block'
       CALL dbcsr_remove_block (matrix_b, row, col, blk_nze)
       CALL cp_dbcsr_finalize (matrix_b, error=error)
       CALL cp_dbcsr_verify_matrix (matrix_b, error, local=.TRUE.)
       checksum = cp_dbcsr_checksum(matrix_b, error=error, local=.TRUE.)
       IF (ldbg.AND.(output_unit>0)) THEN
          WRITE(output_unit, '(T2,A,F12.4)')"DBCSR x wrng cks4", checksum
       ENDIF
       !CALL cp_dbcsr_print (matrix_b, error)
       !write(*,*)'use block',row,col,tr
       !write(*,'(8(1X,F7.2))')block2d_2
       CALL dbcsr_put_block (matrix_b, row, col, block2d_2, tr)
       CALL cp_dbcsr_finalize (matrix_b, error=error)
       CALL cp_dbcsr_verify_matrix (matrix_b, error, local=.TRUE.)
       !CALL cp_dbcsr_print (matrix_b, error)
       checksum = cp_dbcsr_checksum(matrix_b, error=error, local=.TRUE.)
       IF (ldbg.AND.(output_unit>0)) THEN
          WRITE(output_unit, '(T2,A,F12.4)')"DBCSR x checksm4L", checksum
       ENDIF
       DEALLOCATE (block2d_2)
    ENDIF
    CALL dbcsr_iterator_stop (iter)
    CALL cp_dbcsr_verify_matrix (matrix_b, error)
    checksum = cp_dbcsr_checksum(matrix_b, error=error)
    IF (ldbg.AND.(output_unit>0)) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR x checksm4g", checksum
    ENDIF
    CALL cp_dbcsr_finalize (matrix_b, error=error)
    CALL cp_dbcsr_verify_matrix (matrix_b, error)
    checksum = cp_dbcsr_checksum(matrix_b, error=error)
    IF (ldbg.AND.(output_unit>0)) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR x checksum4", checksum
    ENDIF
    !
    ! test the scale routine on a general matrix
    CALL cp_dbcsr_scale(matrix_b,alpha_scalar=3.0_dp,error=error)
    checksum = cp_dbcsr_checksum(matrix_b, error=error)
    IF (ldbg.AND.(output_unit>0)) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR 3*x  checksum", checksum
    ENDIF
    CALL cp_dbcsr_scale(matrix_b,alpha_scalar=1.0_dp/3.0_dp,error=error)
    checksum = cp_dbcsr_checksum(matrix_b, error=error)
    IF (ldbg.AND.(output_unit>0)) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR 3*x/3  checksum", checksum
    ENDIF
    !
    ! Replication
    CALL cp_dbcsr_init (matrix_bt, error)
    CALL cp_dbcsr_copy (matrix_bt, matrix_b, error=error)
    CALL cp_dbcsr_replicate_all (matrix_bt, error=error)
    CALL cp_dbcsr_verify_matrix (matrix_bt, error)
    checksum = cp_dbcsr_checksum(matrix_bt, error=error)
    IF (ldbg.AND.(output_unit>0)) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR x.*p checksum", checksum
    ENDIF
    ! Sum the replicates.
    CALL cp_dbcsr_sum_replicated (matrix_bt, error)
    CALL cp_dbcsr_verify_matrix (matrix_bt, error)
    checksum = cp_dbcsr_checksum(matrix_bt, error=error)
    IF (ldbg.AND.(output_unit>0)) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSRpx.*p checksum", checksum
    ENDIF
    ! Again, to test unity of distribution of replication
    !CALL cp_dbcsr_release(matrix_bt)
    CALL cp_dbcsr_copy (matrix_bt, matrix_b, error=error)
    CALL cp_dbcsr_replicate_all (matrix_bt, error=error)
    CALL cp_dbcsr_verify_matrix (matrix_bt, error)
    CALL cp_dbcsr_distribute (matrix_bt, error=error)
    CALL cp_dbcsr_verify_matrix (matrix_bt, error)
    checksum = cp_dbcsr_checksum(matrix_bt, error=error)
    IF (ldbg.AND.(output_unit>0)) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR x s  checksum", checksum
    ENDIF
    !CALL cp_dbcsr_release(matrix_bt)
    !
    CALL cp_dbcsr_copy (matrix_bt, matrix_b, error=error)
    CALL cp_dbcsr_replicate_all (matrix_bt, error=error)
    CALL cp_dbcsr_verify_matrix (matrix_bt, error)
    CALL cp_dbcsr_distribute (matrix_bt, fast=.TRUE., error=error)
    CALL cp_dbcsr_verify_matrix (matrix_bt, error)
    !checksum = cp_dbcsr_checksum(matrix_bt, error=error)
    IF (ldbg.AND.(output_unit>0)) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR x d  checksum", checksum
    ENDIF

    !CALL cp_dbcsr_release (matrix_bt)
    ! Transpose
    CALL cp_dbcsr_transposed (matrix_bt, matrix_b, error=error)
    CALL cp_dbcsr_verify_matrix (matrix_bt, error)
    !
    CALL dbcsr_distribution_new (product_dist,&
         dbcsr_distribution_mp (matrix_a%m%dist),&
         dbcsr_distribution_row_dist (matrix_a%m%dist),&
         dbcsr_distribution_col_dist (matrix_b%m%dist))
    !
    walltime1 = m_walltime()
    flop = 0
    CALL cp_dbcsr_init (matrix_c, error)
    CALL cp_dbcsr_create (matrix_c, "C", product_dist, dbcsr_type_normal,&
         matrix_a%m%row_blk_size, matrix_b%m%col_blk_size,&
         error=error)
    CALL cp_dbcsr_multiply("N", "N",&
         1.0_dp, matrix_a, matrix_b,&
         0.0_dp, matrix_c,&
         error=error, flop=flop)
    !CALL dbcsr_multiply("N", "N", matrix_bt, matrix_a, matrix_c,&
    !     dist_result, hint, error)
    !CALL dbcsr_distribution_release (dist_result)
    walltime2 = m_walltime()
    newtime = walltime2-walltime1
    IF (newtime .LT. 0.1_dp) newtime = 0.1_dp
    !
    CALL cp_dbcsr_verify_matrix(matrix_c, error)
    !
    checksum = cp_dbcsr_checksum(matrix_c, error=error)
    !WRITE(*,*)'BCSR Sx checksum',checksum,'time=',newtime
    CALL mp_sum (flop, dbcsr_mp_group(dbcsr_distribution_mp(&
         dbcsr_distribution(matrix_c))))
    IF (ldbg.AND.(output_unit>0)) THEN
       WRITE(output_unit, '(T2,A,F12.4,A,F12.4,";",F12.4,A)')&
            "DBCSR Sx checksum", checksum,&
            " time ", newtime,&
            REAL(flop, dp)/1000000000.0_dp/newtime, " GFLOP/s"
       !WRITE(output_unit,'(A,I30,5X,EN12.4)')'FLOP=',flop, REAL (flop, dp)
    ENDIF
    !
    ! Verify copy_dbcsr_to_fm
    CALL cp_dbcsr_copy(matrix_bt, matrix_b, error=error)
    CALL dbcsr_iterator_start(iter,matrix_bt)
    DO WHILE(dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter,row,col,k, &
            row_size=row_size, col_size=col_size)
       blk_nze = row_size * col_size
       !write(*,*)'i have',row,col
       IF (MOD(row+col,3) .EQ. 0) &
            CALL dbcsr_remove_block(matrix_bt,row,col,blk_nze,k)
    ENDDO
    CALL dbcsr_iterator_stop(iter)
    !call dbcsr_remove_block(matrix_bt, 1, 1)
    !call dbcsr_finalize(matrix_bt)
    CALL cp_dbcsr_finalize(matrix_bt, error=error)
    CALL cp_dbcsr_verify_matrix(matrix_bt, error)
    CALL copy_fm_to_dbcsr(matrix_x, matrix_bt, error=error)
    CALL cp_dbcsr_verify_matrix (matrix_bt, error)
    checksum = cp_dbcsr_checksum(matrix_bt, error=error)
    IF (ldbg.AND.(output_unit>0)) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"copy fm2dbcsr chksm", checksum
    ENDIF

    CALL cp_dbcsr_release (matrix_bt, error=error)
    ! check trace(c'*s*c)
    !call dbcsr_replicate_all (matrix_b)
    !call dbcsr_replicate_all (matrix_c)
    CALL cp_dbcsr_trace(matrix_b, matrix_c, trace, error=error)
    IF (ldbg.AND.(output_unit>0)) THEN
       WRITE(output_unit, '(T2,A,F12.4,A,I4)')"DBCSR trace(C'*S*C)", trace," should be equal to ",k
    ENDIF

    CALL cp_dbcsr_release (matrix_c, error=error)

    ! Test redistribution
    CALL dbcsr_create_dist_elemental (weird_dist,&
         dbcsr_nfullrows_total(matrix_b),&
         dbcsr_nfullcols_total(matrix_b),&
         dbcsr_distribution_mp(dist_right),&
         row_blk_size, col_blk_size)
    CALL cp_dbcsr_create (matrix_bt, "Elemental", weird_dist, dbcsr_type_normal,&
         row_blk_size, col_blk_size, 0, 0,error=error)
    CALL dbcsr_distribution_release (weird_dist)
    CALL array_release (row_blk_size)
    CALL array_release (col_blk_size)

    CALL dbcsr_complete_redistribute (matrix_b, matrix_bt, error=dbcsr_error)
    CALL cp_dbcsr_verify_matrix (matrix_bt, error)
    checksum = cp_dbcsr_checksum(matrix_bt, error=error)
    IF (ldbg.AND.(output_unit>0)) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"Elemental B chcksum", checksum
    ENDIF
    !
    CALL dbcsr_create_dist_fullest (weird_dist,&
         dbcsr_nfullrows_total(matrix_b),&
         dbcsr_nfullcols_total(matrix_b),&
         dbcsr_distribution_mp(dist_right),&
         row_blk_size, col_blk_size)
    CALL cp_dbcsr_create (matrix_bt, "Fullest", weird_dist, dbcsr_type_normal,&
         row_blk_size, col_blk_size, 0, 0,error=error)
    CALL dbcsr_distribution_release (weird_dist)
    CALL array_release (row_blk_size)
    CALL array_release (col_blk_size)
    CALL dbcsr_complete_redistribute (matrix_b, matrix_bt, error=dbcsr_error)
    CALL cp_dbcsr_verify_matrix (matrix_bt, error)
    checksum = cp_dbcsr_checksum(matrix_bt, error=error)
    IF (ldbg.AND.(output_unit>0)) THEN
       WRITE(output_unit, '(T2,A,F12.4)')" Fullest B checksum", checksum
    ENDIF

    ! Test dbcsr->fm conversion
    CALL cp_fm_get_info(matrix_x, ncol_global=k, error=error)
    CALL cp_fm_create(matrix_x2, matrix_x%matrix_struct,&
         "x reconversion", error=error)
    CALL cp_fm_create(fmout, matrix_x%matrix_struct, "FMout", error=error)
    CALL copy_dbcsr_to_fm (matrix_b, matrix_x2, error=error)
    CALL cp_sm_fm_multiply(matrix_s, matrix_x2, fmout, k, error=error)
    CALL cp_fm_frobenius_norm(fmout, checksum, error)
    IF (ldbg.AND.(output_unit>0)) THEN
       WRITE(output_unit, '(T2,A,F12.4,A,F12.4)')" Conv Sx checksum",&
            checksum*checksum
    ENDIF
    CALL cp_fm_release (fmout, error)
    CALL cp_fm_release (matrix_x2, error)
    
    

    CALL cp_dbcsr_release (matrix_b, error=error)
    CALL cp_dbcsr_init (matrix_b, error)
    CALL cp_dbcsr_release (matrix_bt, error=error)
    CALL cp_dbcsr_init (matrix_bt, error)

    CALL timestop(handle)

    CALL timeset(routineN//'std', handle)

    CALL cp_fm_get_info(matrix_x, ncol_global=k, error=error)
    CALL cp_fm_create(fmout, matrix_x%matrix_struct, "FMout", error=error)
    walltime1 = m_walltime()
    CALL cp_sm_fm_multiply(matrix_s, matrix_x, fmout, k, error=error)
    walltime2 = m_walltime()
    oldtime = walltime2-walltime1
    CALL cp_fm_frobenius_norm(fmout, checksum, error)
    !WRITE(*,*)'Sx checksum',checksum*checksum,'time=',oldtime
    IF (ldbg.AND.(output_unit>0)) THEN
       WRITE(output_unit, '(T2,A,F12.4,A,F12.4)')" SMFM Sx checksum",&
            checksum*checksum, " time ", oldtime
    ENDIF
    CALL cp_dbcsr_from_fm(matrix_b, fmout, 0.0_dp, dist_right,&
         row_blk_size=matrix_a%m%row_blk_size,&
         col_blk_size=col_blk_size_right,&
         error=error)

    !CALL cp_fm_release (fmout, error)
    !CALL cp_fm_create(fmout, matrix_x%matrix_struct, "FMout", error=error)
    CALL cp_dbcsr_sm_fm_multiply(matrix_a, matrix_x, fmout, k, beta=0.0_dp, error=error)
    CALL cp_fm_frobenius_norm(fmout, checksum, error)
    IF (ldbg.AND.(output_unit>0)) THEN
       WRITE(output_unit, '(T2,A,F12.4)')" dbcsr_SM_FM Sx checksum",&
            checksum*checksum
    ENDIF

    ALLOCATE (matrix_a2)
    CALL cp_dbcsr_init (matrix_a2, error)
    CALL cp_dbcsr_copy(matrix_a2, matrix_a, error=error)
    CALL cp_dbcsr_plus_fm_fm_t(matrix_a2,matrix_x,ncol=k,error=error)
    checksum = cp_dbcsr_checksum(matrix_a2, error=error)
    IF (ldbg.AND.(output_unit>0)) THEN
       WRITE(output_unit, '(T2,A,F12.4)')" S = S+C*C^t Sx1checksum",&
            checksum
    ENDIF
    CALL cp_dbcsr_release (matrix_a2, error=error)
    DEALLOCATE (matrix_a2)

    ALLOCATE (matrix_a2)
    CALL cp_dbcsr_init (matrix_a2, error)
    CALL cp_dbcsr_copy(matrix_a2, matrix_a, error=error)
    CALL cp_dbcsr_plus_fm_fm_t(matrix_a2,matrix_x,matrix_x,ncol=k,error=error)
    checksum = cp_dbcsr_checksum(matrix_a2, error=error)
    IF (ldbg.AND.(output_unit>0)) THEN
       WRITE(output_unit, '(T2,A,F12.4)')" S = S+C*C^t Sx2checksum",&
            checksum
    ENDIF
    CALL cp_dbcsr_release (matrix_a2, error=error)
    DEALLOCATE (matrix_a2)



    NULLIFY (matrix_s2)
    CALL replicate_matrix(matrix_s, matrix_s2, "copy of S",&
         allocate_blocks=.TRUE., error=error)
    CALL cp_sm_plus_fm_fm_t(matrix_s2,matrix_x,&
         ncol=dbcsr_nfullcols_total(matrix_b),&
         error=error)
    !checksum = cp_dbcsr_checksum(matrix_a2, error=error)
    checksum = checksum_matrix(&
         matrix_s2, matrix_s%distribution_2d%blacs_env%para_env)
    IF (ldbg.AND.(output_unit>0)) THEN
       WRITE(output_unit, '(T2,A,F12.4)')" S = S+C*C^t Sx3checksum",&
            checksum
    ENDIF
    CALL deallocate_matrix (matrix_s2, error=error)
    
    CALL array_release (col_blk_size_right)
    CALL cp_dbcsr_release (matrix_a, error=error)
    CALL cp_dbcsr_release (matrix_b, error=error)
    CALL cp_fm_release (fmout, error)
    CALL dbcsr_distribution_release (dist_right)
    CALL dbcsr_distribution_release (product_dist)

    CALL timestop(handle)


  END SUBROUTINE cp_dbcsr_check

! *****************************************************************************
!> \brief multiply a dbcsr with a column vector
!>        c (packed column) = A (dbscr) * B_b (dbcsr)
!> \param[in]  matrix_a  DBCSR matrix
!> \param[in]  matrix_b  DBCSR matrix
!> \param[in]  b_col     the column of the matrix B
!> \param[out] pkd       the result in a packed format
!> 
! *****************************************************************************
  SUBROUTINE cp_dbcsr_multiply_vec(matrix_a, matrix_b, a_row_beg, a_row_end, b_col, pkd, error)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a, matrix_b
    INTEGER, INTENT(in)                      :: a_row_beg, a_row_end, b_col
    REAL(dp), DIMENSION(:), POINTER          :: pkd
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_multiply_vec', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: col, nblkcols_total, node, &
                                                pcol, prow, row, timing_handle
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: tr
    REAL(dp), DIMENSION(:), POINTER          :: pkd_b
    TYPE(dbcsr_mp_obj)                       :: mp_obj

    CALL timeset(routineN, timing_handle)

    nblkcols_total = dbcsr_nblkcols_total(matrix_a)
    col_blk_size => array_data (matrix_a%m%col_blk_size)
    row_blk_size => array_data (matrix_a%m%row_blk_size)

    mp_obj = dbcsr_distribution_mp (dbcsr_distribution (matrix_b))


    ALLOCATE(pkd_b(SIZE(pkd)))
    !
    ! packed the b_col
    CALL cp_dbcsr_pack_vec(matrix_b, b_col, pkd_b, 'column', error)
    !
    ! send the packed col to the right guy
    row = 1
    col = b_col
    !write(*,*) 'cp_dbcsr_multiply_vec: row, col',row,col
    tr = .FALSE.
    CALL dbcsr_get_stored_coordinates (matrix_b, row, col, tr, node)
    !write(*,*) 'cp_dbcsr_multiply_vec: node',node
    CALL dbcsr_mp_get_coordinates(mp_obj, node, prow, pcol)
    !write(*,*) 'cp_dbcsr_multiply_vec: prow, pcol',prow, pcol
    CALL packed_vec_alltoall(pkd_b, pcol, 'column', nblkcols_total, &
         row_blk_size*col_blk_size(b_col), mp_obj, error)
    !
    ! local multiply
    CALL cp_dbcsr_mult_pack_vec_local(matrix_a, pkd_b, a_row_beg, a_row_end, b_col, pkd, error)
    !
    ! sum the local products
    CALL packed_vec_alltoall(pkd, -1, 'all', nblkcols_total, &
         row_blk_size*col_blk_size(b_col), mp_obj, error)
    DEALLOCATE(pkd_b)

    CALL timestop(timing_handle)

  END SUBROUTINE cp_dbcsr_multiply_vec

! *****************************************************************************
!> \brief multiply a dbcsr with a packed vector as
!>        c (packed column) = A (dbscr) * b (packed column)
!> \param[in]  matrix      DBCSR matrix
!> \param[in]  pkd_b       packed vector b
!> \param[in]  b_col       the column of the packed b (this should be changed in b_col_size)
!> \param[out] pkd_c       packed vector c
!> 
! *****************************************************************************
  SUBROUTINE cp_dbcsr_mult_pack_vec_local(matrix_a, pkd_b, a_row_beg, a_row_end, b_col, pkd_c, error)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a
    REAL(dp), DIMENSION(:), INTENT(in)       :: pkd_b
    INTEGER, INTENT(in)                      :: a_row_beg, a_row_end, b_col
    REAL(dp), DIMENSION(:), INTENT(inout)    :: pkd_c
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_mult_pack_vec_local', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: b_col_size, blk, c_offset, &
                                                col, col_size, &
                                                nblkrows_total, row, &
                                                row_size, timing_handle
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size
    LOGICAL                                  :: tr
    REAL(dp), DIMENSION(:), POINTER          :: data_p
    TYPE(dbcsr_iterator)                     :: iter

    CALL timeset(routineN, timing_handle)

    nblkrows_total = dbcsr_nblkrows_total(matrix_a)
    col_blk_size => array_data (matrix_a%m%col_blk_size)

    b_col_size = col_blk_size(b_col)

    IF(pkd_b(nblkrows_total+1).LE.0) CALL stop_program(routineN,"1pkd_b(n+1).LE.0")

    CALL packed_vec_ini(pkd_c, nblkrows_total, error)
    !pkd_c(:) = 0.0_dp
    c_offset = nblkrows_total + 2
    !pkd_c(nblkrows_total+1) = REAL(c_offset,dp)

    CALL dbcsr_iterator_start(iter, matrix_a)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, data_p, tr, blk, &
            row_size=row_size, col_size=col_size)
       IF(row.LT.a_row_beg.OR.row.GT.a_row_end) CYCLE
       IF(INT(pkd_b(col)).LE.0) CYCLE
       !
       IF(INT(pkd_c(row)).LE.0) THEN
          pkd_c(row) = REAL(c_offset,dp)
          pkd_c( INT(pkd_c(row)):INT(pkd_c(row))+row_size*b_col_size-1 ) = 0.0_dp
          c_offset = c_offset + row_size*b_col_size
          pkd_c(nblkrows_total+1) = REAL(c_offset,dp)
       ENDIF
       CALL dgemm('N','N',row_size,b_col_size,col_size,&
            &     1.0_dp,data_p(1),row_size,&
            &            pkd_b(INT(pkd_b(col))),col_size,&
            &     1.0_dp,pkd_c(INT(pkd_c(row))),row_size)
    ENDDO
    CALL dbcsr_iterator_stop(iter)

    IF(pkd_b(nblkrows_total+1).LE.0) CALL stop_program(routineN,"2pkd_b(n+1).LE.0")
    IF(pkd_c(nblkrows_total+1).LE.0) CALL stop_program(routineN,"2pkd_c(n+1).LE.0")

    CALL timestop(timing_handle)

  END SUBROUTINE cp_dbcsr_mult_pack_vec_local

! *****************************************************************************
!> \brief multiply a dbcsr with a replicated array
!>        c = alpha_scalar * A (dbscr) * b + c 
!> \param[in]  matrix       DBCSR matrix
!> \param[in]  vec_b        vectors b
!> \param[inout] vec_c      vectors c
!> \param[in]  ncol         nbr of columns
!> \param[in]  alpha        alpha
!> 
! *****************************************************************************
  SUBROUTINE cp_dbcsr_multiply_local(matrix_a, vec_b, vec_c, ncol, alpha, error)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: vec_b
    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: vec_c
    INTEGER, INTENT(in), OPTIONAL            :: ncol
    REAL(dp), INTENT(IN), OPTIONAL           :: alpha
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_multiply_local', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, my_ncol, row, &
                                                timing_handle
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: offset
    INTEGER, DIMENSION(:), POINTER           :: row_blk_size
    LOGICAL                                  :: tr
    REAL(dp)                                 :: my_alpha, my_alpha2
    REAL(dp), DIMENSION(:, :), POINTER       :: data_d
    TYPE(dbcsr_iterator)                     :: iter

    CALL timeset(routineN, timing_handle)
    
    my_alpha = 1.0_dp
    IF (PRESENT(alpha)) my_alpha = alpha

    my_ncol = SIZE(vec_b,2)
    IF(PRESENT(ncol)) my_ncol = ncol

    my_alpha2 = 0.0_dp
    IF(dbcsr_get_matrix_type(matrix_a).EQ.dbcsr_type_symmetric) my_alpha2 = my_alpha
    IF(dbcsr_get_matrix_type(matrix_a).EQ.dbcsr_type_antisymmetric) my_alpha2 = -my_alpha

    !> this should be precomputed somewhere else
    row_blk_size => array_data (matrix_a%m%row_blk_size)
    ALLOCATE(offset(matrix_a%m%nblkrows_total))
    offset(1) = 1
    DO row = 2,matrix_a%m%nblkrows_total
       offset(row) = offset(row-1) + row_blk_size(row-1)
    ENDDO
    !<

    CALL dbcsr_iterator_start(iter, matrix_a)

    DO WHILE (dbcsr_iterator_blocks_left(iter))

       CALL dbcsr_iterator_next_block(iter, row, col, data_d, tr, blk)

       CALL dgemm('N','N',SIZE(data_d,1),my_ncol,SIZE(data_d,2),&
                  my_alpha,data_d(1,1),SIZE(data_d,1),&
                  vec_b(offset(col),1),SIZE(vec_b,1),&
                  1.0_dp,vec_c(offset(row),1),SIZE(vec_c,1))

       IF((dbcsr_get_matrix_type(matrix_a).EQ.dbcsr_type_symmetric.OR.&
          dbcsr_get_matrix_type(matrix_a).EQ.dbcsr_type_antisymmetric)) THEN
          IF(row.NE.col) THEN
             CALL dgemm('T','N',SIZE(data_d,1),my_ncol,SIZE(data_d,2),&
                  my_alpha2,data_d(1,1),SIZE(data_d,1),&
                  vec_b(offset(row),1),SIZE(vec_b,1),&
                  1.0_dp,vec_c(offset(col),1),SIZE(vec_c,1))
          ENDIF
       ENDIF
    ENDDO

    CALL dbcsr_iterator_stop(iter)

    DEALLOCATE(offset)

    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_multiply_local


! *****************************************************************************
!> \brief multiply a dbcsr with a fm matrix
!> \param[in]  matrix         DBCSR matrix
!> \param[]    fm_in, fm_out  fm matrices
!> \param[in]  ncol           nbr of columns
!> \param[in]  alpha          alpha
!> \param[in]  beta           beta
!> 
! *****************************************************************************
  SUBROUTINE cp_dbcsr_sm_fm_multiply(matrix, fm_in, fm_out, ncol, alpha, beta, error)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    TYPE(cp_fm_type), POINTER                :: fm_in, fm_out
    INTEGER, INTENT(IN)                      :: ncol
    REAL(dp), INTENT(IN), OPTIONAL           :: alpha, beta
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_sm_fm_multiply', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: k_in, k_out, timing_handle, timing_handle_mult
    TYPE(dbcsr_obj)                          :: in, out, out_copy
    TYPE(dbcsr_distribution_obj)             :: dist_right_in, &
                                                dist_right_out, product_dist
    TYPE(array_i1d_obj)                      :: col_blk_size_right_in, &
                                                col_blk_size_right_out
    REAL(dp)                                 :: my_alpha, my_beta, tmp
    INTEGER, DIMENSION(:), POINTER           :: in_col_blk_sizes,&
                                                out_col_blk_sizes

    CALL timeset(routineN, timing_handle)

    my_alpha=1.0_dp
    my_beta=0.0_dp
    IF (PRESENT(alpha)) my_alpha=alpha
    IF (PRESENT(beta)) my_beta=beta

    CALL cp_fm_get_info(fm_out, ncol_global=k_out, error=error)

    CALL cp_fm_get_info(fm_in, ncol_global=k_in, error=error)
    !write(*,*)routineN//" -----------------------------------"
    !IF (k_in .NE. k_out) &
    !   WRITE(*,'(3(A,I5,1X),2(A,F5.2,1X))')&
    !   routineN//" ncol", ncol,'k_in',k_in,'k_out',k_out,&
    !   'alpha',my_alpha,'beta',my_beta

    IF (ncol.GT.0.AND.k_out.GT.0.AND.k_in.GT.0) THEN

       CALL dbcsr_create_dist_r_unrot (dist_right_in, matrix%m%dist, k_in, &
            col_blk_size_right_in)
       CALL cp_dbcsr_init(in, error)
       CALL cp_dbcsr_create(in, "D", dist_right_in, dbcsr_type_normal, &
            matrix%m%row_blk_size, col_blk_size_right_in,&
            0, 0, error=error)
       CALL cp_dbcsr_finalize(in, error=error)

       CALL cp_dbcsr_init(out, error)
       CALL dbcsr_distribution_new (product_dist,&
            dbcsr_distribution_mp (matrix%m%dist),&
            dbcsr_distribution_row_dist (matrix%m%dist),&
            dbcsr_distribution_col_dist (dist_right_in))
       in_col_blk_sizes => array_data (col_blk_size_right_in)
       CALL array_nullify (col_blk_size_right_out)
       CALL array_new (col_blk_size_right_out, in_col_blk_sizes, lb=1)
       out_col_blk_sizes => array_data (col_blk_size_right_out)
       CALL match_col_sizes (out_col_blk_sizes, in_col_blk_sizes, k_out)

       !if (k_in .ne. k_out) then
       !   write(*,*)routineN//" in cs", in_col_blk_sizes
       !   write(*,*)routineN//" out cs", out_col_blk_sizes
       !endif

       CALL cp_dbcsr_create(out, "D", product_dist, dbcsr_type_normal, &
            matrix%m%row_blk_size, col_blk_size_right_out,&
            0, 0, error=error)
       CALL cp_dbcsr_finalize(out, error=error)

       CALL copy_fm_to_dbcsr(fm_in, in, error=error)
       IF(ncol.NE.k_out.OR.my_beta.NE.0.0_dp) &
            CALL copy_fm_to_dbcsr(fm_out, out, error=error)

       CALL timeset(routineN//' (core)', timing_handle_mult)
       CALL cp_dbcsr_multiply("N", "N", my_alpha, matrix, in,&
            my_beta, out, last_column=ncol, error=error)
       CALL timestop(timing_handle_mult)

       CALL copy_dbcsr_to_fm(out, fm_out,error)

       CALL cp_dbcsr_release(in, error=error)
       CALL cp_dbcsr_release(out, error=error)
       CALL array_release(col_blk_size_right_in)
       CALL array_release(col_blk_size_right_out)
       CALL dbcsr_distribution_release(dist_right_in)
       CALL dbcsr_distribution_release(product_dist)

    ENDIF

    CALL timestop(timing_handle)

  CONTAINS
    SUBROUTINE match_col_sizes (sizes1, sizes2, full_num)
    INTEGER, DIMENSION(:), INTENT(INOUT)     :: sizes1
    INTEGER, DIMENSION(:), INTENT(IN)        :: sizes2
    INTEGER, INTENT(IN)                      :: full_num

    INTEGER                                  :: left, n1, n2, p, rm, used

      n1 = SIZE(sizes1)
      n2 = SIZE(sizes2)
      CALL cp_assert (n1 .EQ. n2, cp_fatal_level, cp_caller_error,&
           "match_col_sizes", "distributions must be equal!")
      sizes1(1:n1) = sizes2(1:n1)
      used = SUM (sizes1(1:n1))
      ! If sizes1 does not cover everything, then we increase the
      ! size of the last block; otherwise we reduce the blocks
      ! (from the end) until it is small enough.
      IF (used .LT. full_num) THEN
         sizes1(n1) = sizes1(n1) + full_num-used
      ELSE
         left = used - full_num
         p = n1
         DO WHILE (left .GT. 0 .AND. p .GT. 0)
            rm = MIN(left, sizes1(p))
            sizes1(p) = sizes1(p) - rm
            left = left - rm
            p = p-1
         ENDDO
      ENDIF
    END SUBROUTINE match_col_sizes
  END SUBROUTINE cp_dbcsr_sm_fm_multiply

  SUBROUTINE cp_dbcsr_plus_fm_fm_t_native(sparse_matrix,matrix_v,matrix_g,ncol,&
       alpha,error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: sparse_matrix
    TYPE(cp_fm_type), POINTER                :: matrix_v
    TYPE(cp_fm_type), OPTIONAL, POINTER      :: matrix_g
    INTEGER, INTENT(IN)                      :: ncol
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: alpha
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_plus_fm_fm_t_native', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: k, nao, timing_handle
    LOGICAL                                  :: check_product
    REAL(KIND=dp)                            :: my_alpha, norm
    TYPE(array_i1d_obj)                      :: col_blk_size_left, &
                                                col_dist_left
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: fm_matrix
    TYPE(dbcsr_distribution_obj)             :: dist_left
    TYPE(dbcsr_mp_obj)                       :: mp
    TYPE(dbcsr_obj)                          :: mat_g, mat_v, sparse_matrix2, &
                                                sparse_matrix3

    check_product = .FALSE.

    CALL timeset(routineN, timing_handle)
    IF (ncol .GT. 0) THEN
       CALL cp_assert (dbcsr_valid_index (sparse_matrix), cp_fatal_level,&
            cp_caller_error, routineN, "sparse_matrix must pre-exist", error)
       !
       ! Setup matrix_v
       CALL cp_fm_get_info(matrix_v, ncol_global=k, error=error)
       !WRITE(*,*)routineN//'truncated mult k, ncol',k,ncol,' PRESENT (matrix_g)',PRESENT (matrix_g)
       mp = dbcsr_distribution_mp (dbcsr_distribution(sparse_matrix))
       CALL create_bl_distribution (col_dist_left, col_blk_size_left,&
            k, dbcsr_mp_npcols (mp))
       CALL dbcsr_distribution_new (dist_left, mp,&
            dbcsr_distribution_row_dist (dbcsr_distribution(sparse_matrix)),&
            col_dist_left)
       CALL array_release (col_dist_left)
       CALL cp_dbcsr_init (mat_v, error)
       CALL cp_dbcsr_create(mat_v, "DBCSR matrix_v", dist_left, dbcsr_type_normal,&
            dbcsr_row_block_sizes (sparse_matrix), col_blk_size_left, 0, 0,&
            dbcsr_get_data_type (sparse_matrix), special=.TRUE.,error=error)
       CALL cp_dbcsr_finalize (mat_v, error=error)
       CALL copy_fm_to_dbcsr(matrix_v, mat_v, error=error)
       CALL cp_dbcsr_verify_matrix(mat_v, error)
       !
       ! Setup matrix_g
       IF(PRESENT (matrix_g)) THEN
          CALL cp_dbcsr_init(mat_g, error)
          CALL cp_dbcsr_create(mat_g, "DBCSR matrix_g", dist_left,&
               dbcsr_type_normal,&
               dbcsr_row_block_sizes (sparse_matrix),&
               dbcsr_col_block_sizes (mat_v),&
               data_type=dbcsr_get_data_type (sparse_matrix), special=.TRUE.,error=error)
          CALL cp_dbcsr_finalize (mat_g, error=error)
          CALL copy_fm_to_dbcsr(matrix_g, mat_g, error=error)
       ENDIF
       !
       CALL array_release (col_blk_size_left)
       CALL dbcsr_distribution_release (dist_left)
       !
       !
       IF(check_product) THEN
          NULLIFY(fm_matrix)
          CALL cp_fm_get_info(matrix_v,nrow_global=nao,error=error)
          CALL cp_fm_struct_create(fm_struct_tmp,context=matrix_v%matrix_struct%context,nrow_global=nao,&
               ncol_global=nao,para_env=matrix_v%matrix_struct%para_env,error=error)
          CALL cp_fm_create(fm_matrix,fm_struct_tmp,name="fm matrix",error=error)
          CALL cp_fm_struct_release(fm_struct_tmp,error=error)
          CALL copy_dbcsr_to_fm(sparse_matrix,fm_matrix, error=error)
          CALL cp_dbcsr_init(sparse_matrix3, error)
          CALL cp_dbcsr_copy(sparse_matrix3,sparse_matrix,error=error)
       ENDIF
       !
       my_alpha = 1.0_dp
       IF(PRESENT (alpha)) my_alpha = alpha
       IF(PRESENT (matrix_g)) THEN
          CALL cp_dbcsr_multiply("N", "T", my_alpha, mat_v, mat_g,&
               1.0_dp, sparse_matrix,&
               retain_sparsity=.TRUE.,&
               last_k = ncol,&
               error=error)
       ELSE
          CALL cp_dbcsr_multiply("N", "T", my_alpha, mat_v, mat_v,&
               1.0_dp, sparse_matrix,&
               retain_sparsity=.TRUE.,&
               last_k = ncol,&
               error=error)
       ENDIF

       IF(check_product) THEN
          IF(PRESENT (matrix_g)) THEN
             CALL cp_fm_gemm("N","T",nao,nao,ncol,my_alpha,matrix_v,matrix_g,&
                  1.0_dp,fm_matrix,error=error)
          ELSE
             CALL cp_fm_gemm("N","T",nao,nao,ncol,my_alpha,matrix_v,matrix_v,&
                  1.0_dp,fm_matrix,error=error)
          ENDIF

          CALL cp_dbcsr_init(sparse_matrix2, error)
          CALL cp_dbcsr_copy(sparse_matrix2,sparse_matrix,error=error)
          CALL cp_dbcsr_scale(sparse_matrix2,alpha_scalar=0.0_dp,error=error)
          CALL copy_fm_to_dbcsr(fm_matrix,sparse_matrix2,keep_sparsity=.TRUE., error=error)
          CALL cp_dbcsr_add(sparse_matrix2,sparse_matrix,alpha_scalar=1.0_dp,&
               beta_scalar=-1.0_dp,error=error)
          CALL cp_dbcsr_norm(sparse_matrix2,which_norm=dbcsr_norm_frobenius,&
               norm_scalar=norm,error=error)
          WRITE(*,*) 'nao=',nao,' k=',k,' ncol=',ncol,' my_alpha=',my_alpha
          WRITE(*,*) 'PRESENT (matrix_g)',PRESENT (matrix_g)
          WRITE(*,*) 'matrix_type=',dbcsr_get_matrix_type(sparse_matrix)
          WRITE(*,*) 'norm(sm+alpha*v*g^t - fm+alpha*v*g^t)/n=',norm/REAL(nao,dp)
          IF(norm/REAL(nao,dp).GT.1e-12_dp) THEN
             !WRITE(*,*) 'fm_matrix'
             !DO j=1,SIZE(fm_matrix%local_data,2)
             !   DO i=1,SIZE(fm_matrix%local_data,1)
             !      WRITE(*,'(A,I3,A,I3,A,E26.16,A)') 'a(',i,',',j,')=',fm_matrix%local_data(i,j),';'
             !   ENDDO
             !ENDDO
             !WRITE(*,*) 'mat_v'
             !CALL cp_dbcsr_print(mat_v,matlab_format=.TRUE.)
             !WRITE(*,*) 'mat_g'
             !CALL cp_dbcsr_print(mat_g,matlab_format=.TRUE.)
             !WRITE(*,*) 'sparse_matrix'
             !CALL cp_dbcsr_print(sparse_matrix,matlab_format=.TRUE.)
             !WRITE(*,*) 'sparse_matrix2 (-sm + sparse(fm))'
             !CALL cp_dbcsr_print(sparse_matrix2,matlab_format=.TRUE.)
             !WRITE(*,*) 'sparse_matrix3 (copy of sm input)'
             !CALL cp_dbcsr_print(sparse_matrix3,matlab_format=.TRUE.)
             !stop
          ENDIF
          CALL cp_dbcsr_release(sparse_matrix2, error=error)
          CALL cp_dbcsr_release(sparse_matrix3, error=error)
          CALL cp_fm_release(fm_matrix,error=error)
       ENDIF
       CALL cp_dbcsr_release (mat_v, error=error)
       IF(PRESENT (matrix_g)) CALL cp_dbcsr_release (mat_g, error=error)
    ENDIF
    CALL timestop(timing_handle)

  END SUBROUTINE cp_dbcsr_plus_fm_fm_t_native

! *****************************************************************************
!> \brief pack a row or column of a dbcsr
!> \param[in]  matrix      DBCSR matrix
!> \param[in]  ivec        the vector to pack
!> \param[out] pkd_vec     packed vector
!> \param[in]  what        what to pack (row or column)
!> 
! *****************************************************************************
  SUBROUTINE cp_dbcsr_pack_vec(matrix, ivec, pkd_vec, what, error)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER, INTENT(IN)                      :: ivec
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: pkd_vec
    CHARACTER(LEN=*), INTENT(IN)             :: what
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_pack_vec', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: col, col_size, &
                                                nblkrows_total, nze, row, &
                                                row_size, timing_handle
    LOGICAL                                  :: tr
    REAL(KIND=dp), DIMENSION(:), POINTER     :: data_p
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)
    
    nblkrows_total = dbcsr_nblkrows_total(matrix)
    pkd_vec(1:nblkrows_total+1) = 0.0_dp ! should be big enough to hold all the data
    !
    ! let's go
    SELECT CASE(what)
    CASE('column')
       ! we need nblkrows_total+1 to store rows 
       ! and 1 extra for adding new data
       pkd_vec(nblkrows_total+1) = REAL(nblkrows_total+2,dp)
       CALL dbcsr_iterator_start(iter, matrix)
       DO WHILE (dbcsr_iterator_blocks_left (iter))
          CALL dbcsr_iterator_next_block (iter, row, col, data_p, tr, &
               row_size=row_size, col_size=col_size)
          IF(col.ne.ivec) CYCLE
          nze = row_size*col_size
          !write(*,*) 'row',row,' col',col,' nze',nze
          IF(nze.le.0) CYCLE
          !
          ! let's copy the block
          pkd_vec(row) = pkd_vec(nblkrows_total+1)
          pkd_vec(nblkrows_total+1) = pkd_vec(nblkrows_total+1) + REAL(nze,dp)
          IF(SIZE(pkd_vec).LT.INT(pkd_vec(nblkrows_total+1))-1) THEN
             WRITE(*,*) 'SIZE(pkd_vec)',SIZE(pkd_vec)
             WRITE(*,*) 'pkd_vec(nblkrows_total+1)',pkd_vec(nblkrows_total+1)
             CALL stop_program(routineN,"col: SIZE(pkd_vec).LT.pkd_vec(nblkrows_total+1)")
          ENDIF
          CALL dcopy(nze,data_p(1),1,pkd_vec(INT(pkd_vec(row))),1)
       ENDDO
       CALL dbcsr_iterator_stop(iter)
    CASE('row')
       ! we need nblkrows_total+1 to store rows 
       ! and 1 extra for adding new data
       pkd_vec(nblkrows_total+1) = REAL(nblkrows_total+2,dp)
       CALL dbcsr_iterator_start(iter, matrix)
       DO WHILE (dbcsr_iterator_blocks_left (iter))
          CALL dbcsr_iterator_next_block (iter, row, col, data_p, tr, &
               row_size=row_size, col_size=col_size)
          IF(row.ne.ivec) CYCLE
          nze = row_size*col_size
          !write(*,*) 'row',row,' col',col,' nze',nze
          IF(nze.le.0) CYCLE
          !
          ! let's copy the block
          pkd_vec(col) = pkd_vec(nblkrows_total+1)
          pkd_vec(nblkrows_total+1) = pkd_vec(nblkrows_total+1) + nze
          IF(SIZE(pkd_vec).LT.INT(pkd_vec(nblkrows_total+1))-1) &
               CALL stop_program(routineN,"SIZE(pkd_vec).LT.pkd_vec(nblkrows_total+1)")
          CALL dcopy(nze,data_p(1),1,pkd_vec(INT(pkd_vec(col))),1)
       ENDDO
       CALL dbcsr_iterator_stop(iter)
    CASE DEFAULT
       CALL stop_program(routineN,"pack what?")       
    END SELECT

    CALL timestop(timing_handle)

  END SUBROUTINE cp_dbcsr_pack_vec

! *****************************************************************************
!> \brief unpack a row or column of a dbcsr
!> \param[inout]  matrix      DBCSR matrix
!> \param[in]     ivec        the vector to pack
!> \param[in]     pkd_vec     packed vector
!> \param[in]     what        what to unpack (row or column)
!> 
! *****************************************************************************
  SUBROUTINE cp_dbcsr_unpack_vec(matrix, ivec, pkd_vec, what, do_sum, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: ivec
    REAL(dp), DIMENSION(:), INTENT(IN)       :: pkd_vec
    CHARACTER(LEN=*), INTENT(IN)             :: what
    LOGICAL, INTENT(in), OPTIONAL            :: do_sum
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_unpack_vec', &
      routineP = moduleN//':'//routineN

    INTEGER :: col, col_size, mynode, nblkcols_total, nblkrows_total, node, &
      nze, row, row_size, timing_handle
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: my_do_sum, tr

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)
    
    my_do_sum = .FALSE.
    IF(PRESENT(do_sum))my_do_sum=do_sum

    row_blk_size => array_data (dbcsr_row_block_sizes (matrix))
    col_blk_size => array_data (dbcsr_col_block_sizes (matrix))
    nblkrows_total = dbcsr_nblkrows_total(matrix)
    nblkcols_total = dbcsr_nblkcols_total(matrix)
    mynode = dbcsr_mp_mynode (dbcsr_distribution_mp (dbcsr_distribution (matrix)))

    !
    ! let's go
    SELECT CASE(what)
    CASE('column')
       IF(INT(pkd_vec(nblkrows_total+1)).LE.0) &
            CALL stop_program(routineN,"1pkd_vec(n+1).LE.0")
       col = ivec
       col_size = col_blk_size(col)
       DO row = 1, nblkrows_total
          IF(INT(pkd_vec(row)).LE.0) CYCLE
          row_size = row_blk_size(row)
          nze = row_size*col_size
          ! add block only on the correct node
          tr = .FALSE.
          CALL dbcsr_get_stored_coordinates (matrix, row, col, tr, node)
          IF(node.EQ.mynode) CALL dbcsr_put_block(matrix, row, col, &
               pkd_vec(INT(pkd_vec(row)):INT(pkd_vec(row)+nze-1)), &
               summation=my_do_sum)
       ENDDO ! row
       CALL cp_dbcsr_finalize(matrix, error=error)
    CASE('row')
       IF(INT(pkd_vec(nblkrows_total+1)).LE.0) &
            CALL stop_program(routineN,"2pkd_vec(n+1).LE.0")
       row = ivec
       row_size = row_blk_size(row)
       DO col = 1,nblkcols_total
          IF(INT(pkd_vec(col)).LE.0) CYCLE
          col_size = col_blk_size(col)
          nze = row_size*col_size
          ! add block only on the correct node
          tr = .FALSE.
          CALL dbcsr_get_stored_coordinates (matrix, row, col, tr, node)
          IF(node.EQ.mynode) CALL dbcsr_put_block(matrix, row, col, &
               pkd_vec(INT(pkd_vec(col)):INT(pkd_vec(col)+nze-1)),&
               summation=my_do_sum)
       ENDDO ! col
       CALL cp_dbcsr_finalize(matrix, error=error)
    CASE DEFAULT
       CALL stop_program(routineN,"unpack what?")
    END SELECT

    CALL timestop(timing_handle)

  END SUBROUTINE cp_dbcsr_unpack_vec


! *****************************************************************************
!> \brief Does one of the operations:
!>      left)  pkd_ji = alpha * block_jj (block) * pkd_ji for j = 1,..., nvec
!>      right) pkd_ij = alpha * pkd_ij * block_jj (block) for j = 1,..., nvec
!>
! *****************************************************************************
  SUBROUTINE packed_vec_scale(alpha, block, pkd, nvec, ivec, vec_blk_size, side, error)
    REAL(dp), INTENT(IN)                     :: alpha
    REAL(dp), DIMENSION(:), INTENT(in)       :: block
    REAL(dp), DIMENSION(:), INTENT(inout)    :: pkd
    INTEGER, INTENT(in)                      :: nvec, ivec
    INTEGER, DIMENSION(:), INTENT(in)        :: vec_blk_size
    CHARACTER(LEN=*), INTENT(IN)             :: side
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'packed_vec_scale', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iblk, k, m, n, offset
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: buff

!   ---------------------------------------------------------------------------

    ALLOCATE(buff(SIZE(pkd)))!too long buff allocated here...
    buff=0.0_dp
    SELECT CASE(side)
    CASE('left')
       iblk = 1
       DO i = 1,nvec
          offset = INT(pkd(i))
          IF(offset.GT.0) THEN
             m = vec_blk_size(i)
             n = vec_blk_size(ivec)
             k = vec_blk_size(i)
             CALL dgemm('N','N',m,n,k,alpha,block(iblk),m,pkd(offset),k,0.0_dp,buff(1),m)
             CALL dcopy(m*n,buff(1),1,pkd(offset),1)
          ENDIF
          iblk = iblk + vec_blk_size(i)**2
       ENDDO
    CASE('right')
       iblk = 1
       DO i = 1,nvec
          offset = INT(pkd(i))
          IF(offset.GT.0) THEN
             m = vec_blk_size(ivec)
             n = vec_blk_size(i)
             k = vec_blk_size(i)
             CALL dgemm('N','N',m,n,k,alpha,pkd(offset),m,block(iblk),k,0.0_dp,buff(1),m)
             CALL dcopy(m*n,buff(1),1,pkd(offset),1)
          ENDIF
          iblk = iblk + vec_blk_size(i)**2
       ENDDO
    CASE DEFAULT
       CALL stop_program(routineN,"side?")
    END SELECT
    DEALLOCATE(buff)

  END SUBROUTINE packed_vec_scale

! *****************************************************************************
!> \brief Does the operations:
!>      v(k,:) = v(k,:) - uAk(i,k)' * inv(d(i,i)) * v(i,:) / s
!>      u(k,:) = u(k,:) -   v(i,k)' * inv(d(i,i)) * u(i,:) / s
!>
! *****************************************************************************
  SUBROUTINE packed_vec_bif_tech(mat_v, mat_u, pkd_v_fac, pkd_u_fac, ivec, pkd_v, pkd_u, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: mat_v, mat_u
    REAL(dp), DIMENSION(:), INTENT(IN)       :: pkd_v_fac, pkd_u_fac
    INTEGER, INTENT(IN)                      :: ivec
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: pkd_v, pkd_u
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'packed_vec_bif_tech', &
      routineP = moduleN//':'//routineN

    INTEGER :: blk, col, col_size, k_row, k_row_size, offset, row, row_size, &
      timing_handle, u_offset, u_offset_last, ufac_offset, v_offset, &
      v_offset_last, vfac_offset
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    REAL(KIND=dp), DIMENSION(:), POINTER     :: data_u, data_v

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)
    
    row_blk_size => array_data (dbcsr_row_block_sizes (mat_u))
    col_blk_size => array_data (dbcsr_col_block_sizes (mat_u))

    IF(pkd_u(mat_u%m%nblkrows_total+1).LE.0) CALL stop_program(routineN,"1pkd_u(n+1).LE.0")
    IF(pkd_v(mat_v%m%nblkrows_total+1).LE.0) CALL stop_program(routineN,"1pkd_v(n+1).LE.0")
    CALL dbcsr_get_data (mat_v, data_v)
    CALL dbcsr_get_data (mat_u, data_u)
    !
    !v(k,:) (pkd) = v(k,:) (pkd) - (inv(d(i,i)) * uAk(i,k))' / s (pkd) * v(i,:) (bcsr) i = 1,ivec-1
    k_row = ivec
    k_row_size = row_blk_size(k_row)
    v_offset_last = INT(pkd_v(mat_v%m%nblkrows_total+1))
    DO row = 1, mat_v%m%nblkrows_total
       IF(row.GE.k_row) EXIT !i = 1,ivec-1
       row_size = row_blk_size(row)
       ufac_offset = INT(pkd_v_fac(row))
       IF(ufac_offset.LE.0) CYCLE
       DO blk = mat_v%m%row_p(row)+1,mat_v%m%row_p(row+1)
          col = mat_v%m%col_i(blk)
          col_size = col_blk_size(col)
          offset = ABS(mat_v%m%blk_p(blk))
          v_offset = INT(pkd_v(col))
          IF(v_offset.LE.0) THEN
             v_offset = v_offset_last
             v_offset_last = v_offset_last + k_row_size * col_size
             pkd_v(v_offset:v_offset+k_row_size * col_size-1) = 0.0_dp
             pkd_v(mat_v%m%nblkrows_total+1) = REAL(v_offset_last,dp)
          ENDIF
          !
          ! let's multiply and add
          CALL dgemm('T', 'N', k_row_size, col_size, row_size,&
               &    -1.0_dp, pkd_v_fac(ufac_offset), row_size,&
               &             data_v(offset), row_size,&
                  &     1.0_dp, pkd_v(v_offset), k_row_size)
       ENDDO
    ENDDO ! row
    !
    !u(k,:) (pkd) = u(k,:) (pkd) - (inv(d(i,i)) * v(i,k))' / s (pkd) * u(i,:) (bcsr) i = 1,ivec-1
    u_offset_last = INT(pkd_u(mat_u%m%nblkrows_total+1))
    DO row = 1, mat_u%m%nblkrows_total
       IF(row.GE.k_row) EXIT !i = 1,ivec-1
       row_size = row_blk_size(row)
       vfac_offset = INT(pkd_u_fac(row))
       IF(vfac_offset.LE.0) CYCLE
       DO blk = mat_u%m%row_p(row)+1,mat_u%m%row_p(row+1)
          col = mat_u%m%col_i(blk)
          col_size = col_blk_size(col)
          offset = ABS(mat_u%m%blk_p(blk))
          u_offset = INT(pkd_u(col))
          IF(u_offset.LE.0) THEN
             u_offset = u_offset_last
             u_offset_last = u_offset_last + k_row_size * col_size
             pkd_u(u_offset:u_offset+k_row_size * col_size-1) = 0.0_dp
             pkd_u(mat_u%m%nblkrows_total+1) = REAL(u_offset_last,dp)
          ENDIF
          !
          ! let's multiply and add
          CALL dgemm('T', 'N', k_row_size, col_size, row_size,&
               &    -1.0_dp, pkd_u_fac(vfac_offset), row_size,&
               &             data_u(offset), row_size,&
               &     1.0_dp, pkd_u(u_offset), k_row_size)
       ENDDO
    ENDDO ! row

    IF(pkd_u(mat_u%m%nblkrows_total+1).LE.0) CALL stop_program(routineN,"2pkd_u(n+1).LE.0")
    IF(pkd_v(mat_v%m%nblkrows_total+1).LE.0) CALL stop_program(routineN,"2pkd_v(n+1).LE.0")

    CALL timestop(timing_handle)

  END SUBROUTINE packed_vec_bif_tech


  SUBROUTINE packed_vec_build_u(pkd_u, pkd_v, k, n, s, vec_blk_size, error)
    REAL(dp), DIMENSION(:), INTENT(OUT)      :: pkd_u
    REAL(dp), DIMENSION(:), INTENT(IN)       :: pkd_v
    INTEGER, INTENT(IN)                      :: k, n
    REAL(dp), INTENT(IN)                     :: s
    INTEGER, DIMENSION(:), INTENT(IN)        :: vec_blk_size
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'packed_vec_build_u', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, i_blk_size, k_blk_size, &
                                                u_offset

    pkd_u(1:n+1) = 0.0_dp
    u_offset = n+2
    !
    ! u(k,:) = [ v(k,1:k-1)/s -1 0 ... 0]
    k_blk_size = vec_blk_size(k)
    DO i = 1,k-1
       IF(INT(pkd_v(i)).GT.0) THEN
          i_blk_size = vec_blk_size(i)
          CALL dcopy(k_blk_size*i_blk_size, pkd_v(INT(pkd_v(i))),1,pkd_u(u_offset),1)
          CALL dscal(k_blk_size*i_blk_size, -1.0_dp/s, pkd_u(u_offset),1)
          pkd_u(i) = u_offset
          u_offset = u_offset + k_blk_size*i_blk_size
       ENDIF
    ENDDO
    !
    CALL block_set(k_blk_size, k_blk_size, pkd_u(u_offset:u_offset+ k_blk_size**2-1), &
         &         1.0_dp, 0.0_dp)
    pkd_u(k) = u_offset
    u_offset = u_offset + k_blk_size**2
    pkd_u(n+1) = u_offset

  END SUBROUTINE packed_vec_build_u

! *****************************************************************************
!> \brief Does the operations:
!>      v(k,:) = v(k,:) - uAk(i,k)' * inv(d(i,i)) * v(i,:) / s
!>
! *****************************************************************************
  SUBROUTINE packed_vec_bif_tech2(mat_v, pkd_v_fac, ivec, pkd_v, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: mat_v
    REAL(dp), DIMENSION(:), INTENT(IN)       :: pkd_v_fac
    INTEGER, INTENT(IN)                      :: ivec
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: pkd_v
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'packed_vec_bif_tech2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: col, col_size, k_row, &
                                                k_row_size, nblkrows_total, &
                                                row, row_size, timing_handle
    INTEGER, DIMENSION(:), POINTER           :: row_blk_size
    LOGICAL                                  :: tr
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: buff
    REAL(KIND=dp), DIMENSION(:), POINTER     :: data_v
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)
    
    nblkrows_total = dbcsr_nblkrows_total(mat_v)
    row_blk_size => array_data (dbcsr_row_block_sizes(mat_v))

    IF(pkd_v(nblkrows_total+1).LE.0) CALL stop_program(routineN,"1pkd_v(n+1).LE.0")

    ALLOCATE(buff(SIZE(pkd_v)))

    CALL packed_vec_ini(buff, nblkrows_total, error)
    !
    !v(k,:) (pkd) = v(k,:) (pkd) - (inv(d(i,i)) * uAk(i,k))' / s (pkd) * v(i,:) (bcsr) i = 1,ivec-1
    k_row = ivec
    k_row_size = row_blk_size(k_row)

    CALL dbcsr_iterator_start(iter, mat_v)
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col, data_v, tr, &
            row_size=row_size, col_size=col_size)
       IF(row.GE.k_row) EXIT !i = 1,ivec-1
       IF(INT(pkd_v_fac(row)).LE.0) CYCLE
       IF(INT(pkd_v(col)).LE.0) THEN
          pkd_v(col) = pkd_v(nblkrows_total+1)
          pkd_v(INT(pkd_v(col)):INT(pkd_v(col)) + k_row_size * col_size-1) = 0.0_dp
          !v_offset = v_offset + k_row_size * col_size
          pkd_v(nblkrows_total+1) = pkd_v(nblkrows_total+1) + REAL(k_row_size * col_size,dp)
       ENDIF
       !
       ! let's multiply and add
       CALL dgemm('T', 'N', k_row_size, col_size, row_size,&
            &    -1.0_dp, pkd_v_fac(INT(pkd_v_fac(row))), row_size,&
            &             data_v(1), row_size,&
            &     1.0_dp, pkd_v(INT(pkd_v(col))), k_row_size)
    ENDDO
    CALL dbcsr_iterator_stop(iter)
    !
    IF(pkd_v(nblkrows_total+1).LE.0) CALL stop_program(routineN,"2pkd_v(n+1).LE.0")
    DEALLOCATE(buff)

    CALL timestop(timing_handle)

  END SUBROUTINE packed_vec_bif_tech2

  SUBROUTINE packed_vec_ini(pkd_vec, n, error)
    REAL(dp), DIMENSION(:), INTENT(OUT)      :: pkd_vec
    INTEGER, INTENT(IN)                      :: n
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    pkd_vec(1:n) = 0.0_dp
    pkd_vec(n+1) = n + 2
  END SUBROUTINE packed_vec_ini

  SUBROUTINE packed_vec_alltoall(pkd_vec, source, scope, n, vec_blk_size, mp_obj, error)
    REAL(dp), DIMENSION(:), POINTER          :: pkd_vec
    INTEGER, INTENT(IN)                      :: source
    CHARACTER(LEN=*), INTENT(IN)             :: scope
    INTEGER, INTENT(IN)                      :: n
    INTEGER, DIMENSION(:), INTENT(in)        :: vec_blk_size
    TYPE(dbcsr_mp_obj)                       :: mp_obj
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'packed_vec_alltoall', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: inode, irow, mp_group, &
                                                mynode, mypcol, myprow, &
                                                npcols, nprows, numnodes, &
                                                src, timing_handle
    INTEGER, DIMENSION(:, :), POINTER        :: blacs2mpi
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: buff, buff2

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)

    numnodes = dbcsr_mp_numnodes (mp_obj)
    mynode = dbcsr_mp_mynode (mp_obj)
    myprow = dbcsr_mp_myprow (mp_obj)
    mypcol = dbcsr_mp_mypcol (mp_obj)
    npcols = dbcsr_mp_npcols (mp_obj)
    nprows = dbcsr_mp_nprows (mp_obj)
    blacs2mpi => dbcsr_mp_pgrid (mp_obj)
    mp_group = dbcsr_mp_group (mp_obj)

    ALLOCATE(buff(SIZE(pkd_vec)))
    ALLOCATE(buff2(10*SIZE(pkd_vec)))! that sucks !

    !buff2(1:int(pkd_vec(n+1))-1) = pkd_vec(1:int(pkd_vec(n+1))-1)
    CALL dcopy(INT(pkd_vec(n+1))-1, pkd_vec(1), 1, buff2(1), 1)

    SELECT CASE(scope)
    CASE('row')
       CALL stop_program(routineN,"NYI")
    CASE('column')
       !
       ! simple hack
       DO irow = 0,nprows-1
          !buff = pkd_vec
          CALL dcopy(INT(pkd_vec(n+1))-1, pkd_vec(1), 1, buff(1), 1)
          src = blacs2mpi(irow,source)
          CALL mp_bcast(buff,src,mp_group)
          IF(src.NE.mynode) THEN
             CALL packed_vec_add(buff2,buff,vec_blk_size,n,error=error)
          ENDIF
       ENDDO
    CASE('all')
       !
       ! simple hack
       DO inode = 0,numnodes-1
          IF(inode.eq.mynode) THEN
             !buff(1:int(pkd_vec(n+1))-1) = pkd_vec(1:int(pkd_vec(n+1))-1)
             CALL dcopy(INT(pkd_vec(n+1))-1, pkd_vec(1), 1, buff(1), 1)
             !else
             !buff(1:int(pkd_vec(n+1))-1) = huge(0.0_dp)
          ENDIF
          !IF(numnodes.ne.1) then
          CALL mp_bcast(buff,inode,mp_group)
          IF(inode.ne.mynode) THEN
             CALL packed_vec_add(buff2,buff,vec_blk_size,n,error=error)
          ENDIF
          !endif
       ENDDO
    CASE DEFAULT
       CALL stop_program(routineN,"how do you wanna bcast")
    END SELECT
    
    IF(SIZE(pkd_vec).lt.INT(buff2(n+1))-1) THEN
       DEALLOCATE(pkd_vec)
       ALLOCATE(pkd_vec(INT(buff2(n+1))))
    ENDIF

    pkd_vec(1:INT(buff2(n+1)) - 1) = buff2(1:INT(buff2(n+1)) - 1)

    DEALLOCATE(buff)
    DEALLOCATE(buff2)
       
    CALL timestop(timing_handle)

  END SUBROUTINE packed_vec_alltoall

  SUBROUTINE packed_vec_bcast(pkd_vec, source, scope, do_summation, n, vec_blk_size, mp_obj, error)
    REAL(dp), DIMENSION(:), POINTER          :: pkd_vec
    INTEGER, INTENT(IN)                      :: source
    CHARACTER(LEN=*), INTENT(IN)             :: scope
    LOGICAL, INTENT(in)                      :: do_summation
    INTEGER, INTENT(IN)                      :: n
    INTEGER, DIMENSION(:), INTENT(in)        :: vec_blk_size
    TYPE(dbcsr_mp_obj)                       :: mp_obj
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'packed_vec_bcast', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: icol, inode, irow, mp_group, &
                                                mynode, mypcol, myprow, &
                                                npcols, nprows, numnodes, &
                                                src, timing_handle
    INTEGER, DIMENSION(:, :), POINTER        :: blacs2mpi
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: buff, buff2

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)

    numnodes = dbcsr_mp_numnodes (mp_obj)
    mynode = dbcsr_mp_mynode (mp_obj)
    myprow = dbcsr_mp_myprow (mp_obj)
    mypcol = dbcsr_mp_mypcol (mp_obj)
    npcols = dbcsr_mp_npcols (mp_obj)
    nprows = dbcsr_mp_nprows (mp_obj)
    blacs2mpi => dbcsr_mp_pgrid (mp_obj)
    mp_group = dbcsr_mp_group (mp_obj)

    ALLOCATE(buff(SIZE(pkd_vec)))
    ALLOCATE(buff2(10*SIZE(pkd_vec)))! that sucks !

    SELECT CASE(scope)
    CASE('rowise')
       !
       ! simple hack
       DO icol = 0,npcols-1
          buff = pkd_vec
          src = blacs2mpi(source,icol)
          CALL mp_bcast(buff,src,mp_group)
          IF(mypcol.EQ.icol) THEN
             pkd_vec = buff
          ENDIF
       ENDDO
    CASE('columnwise')
       !
       ! simple hack
       DO irow = 0,nprows-1
          buff = pkd_vec
          src = blacs2mpi(irow,source)
          CALL mp_bcast(buff,src,mp_group)
          IF(myprow.EQ.irow) THEN
             pkd_vec = buff
          ENDIF
       ENDDO
    CASE('all')
       !
       ! simple hack
       !buff2(1:int(pkd_vec(n+1))-1) = pkd_vec(1:int(pkd_vec(n+1))-1)
       CALL dcopy( INT(pkd_vec(n+1))-1, pkd_vec(1), 1, buff2(1), 1)
       DO inode = 0,numnodes-1
          IF(inode.eq.mynode) THEN
             !buff = pkd_vec
             CALL dcopy( INT(pkd_vec(n+1))-1, pkd_vec(1), 1, buff(1), 1)
          ELSE
             !buff = huge(0.0_dp)
          ENDIF
          CALL mp_bcast(buff,inode,mp_group)
          IF(inode.ne.mynode) THEN
             CALL packed_vec_add(buff2,buff,vec_blk_size,n,do_summation,error=error)
          ENDIF
       ENDDO
       IF(SIZE(pkd_vec).lt.INT(buff2(n+1))-1) THEN
          DEALLOCATE(pkd_vec)
          ALLOCATE(pkd_vec(INT(buff2(n+1))))
       ENDIF
       !pkd_vec(1:int(buff2(n+1)) - 1) = buff2(1:int(buff2(n+1)) - 1)
       CALL dcopy( INT(buff2(n+1)) - 1, buff2(1), 1, pkd_vec(1), 1)

    CASE DEFAULT
       CALL stop_program(routineN,"how do you wanna bcast")
    END SELECT

    DEALLOCATE(buff)
    DEALLOCATE(buff2)

    CALL timestop(timing_handle)

  END SUBROUTINE packed_vec_bcast

  SUBROUTINE packed_vec_reduce(pkd_vec, to, scope, vec_blk_size, n, mp_obj, error)
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: pkd_vec
    INTEGER, INTENT(IN)                      :: to
    CHARACTER(LEN=*), INTENT(IN)             :: scope
    INTEGER, DIMENSION(:), INTENT(IN)        :: vec_blk_size
    INTEGER, INTENT(IN)                      :: n
    TYPE(dbcsr_mp_obj)                       :: mp_obj
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'packed_vec_reduce', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: icol, irow, mp_group, mynode, &
                                                mypcol, myprow, npcols, &
                                                nprows, numnodes, src, &
                                                timing_handle
    INTEGER, DIMENSION(:, :), POINTER        :: blacs2mpi
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: buff, sums

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)

    numnodes = dbcsr_mp_numnodes (mp_obj)
    mynode = dbcsr_mp_mynode (mp_obj)
    myprow = dbcsr_mp_myprow (mp_obj)
    mypcol = dbcsr_mp_mypcol (mp_obj)
    npcols = dbcsr_mp_npcols (mp_obj)
    nprows = dbcsr_mp_nprows (mp_obj)
    blacs2mpi => dbcsr_mp_pgrid (mp_obj)
    mp_group = dbcsr_mp_group (mp_obj)

    IF(pkd_vec(n+1).LE.0) CALL stop_program(routineN,"a_offset_last.LE.0!!!!!!!")

    ALLOCATE(buff(SIZE(pkd_vec)),sums(SIZE(pkd_vec)))
    sums = 0.0_dp
    SELECT CASE(scope)
    CASE('rowise')
       CALL stop_program(routineN,"more work here")
       !
       ! simple hack
       DO icol = 0,npcols-1
          buff = pkd_vec
          src = blacs2mpi(to,icol)
          CALL mp_bcast(buff,src,mp_group)
          IF(mypcol.EQ.icol) THEN
             pkd_vec = buff
             CALL packed_vec_add(sums,buff,vec_blk_size,n,error=error)
          ENDIF
       ENDDO
    CASE('columnwise')
       !
       ! simple hack
       DO irow = 0,nprows-1
          buff = pkd_vec
          src = blacs2mpi(irow,to)
          CALL mp_bcast(buff,src,mp_group)
          IF(myprow.EQ.irow) THEN
             CALL packed_vec_add(sums,buff,vec_blk_size,n,error=error)
          ENDIF
       ENDDO
    CASE('all')
       CALL stop_program(routineN,"more work here")
    CASE DEFAULT
       CALL stop_program(routineN,"how do you wanna reduce")
    END SELECT
    pkd_vec = sums
    DEALLOCATE(buff,sums)

    IF(pkd_vec(n+1).LE.0) CALL stop_program(routineN,"a_offset_last.LE.0!!!!!!!")

    CALL timestop(timing_handle)

  END SUBROUTINE packed_vec_reduce

  SUBROUTINE packed_vec_add(pkd_vec_a,pkd_vec_b,vec_blk_size,n,do_summation,error)
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: pkd_vec_a
    REAL(dp), DIMENSION(:), INTENT(IN)       :: pkd_vec_b
    INTEGER, DIMENSION(:), INTENT(IN)        :: vec_blk_size
    INTEGER, INTENT(IN)                      :: n
    LOGICAL, INTENT(in), OPTIONAL            :: do_summation
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'packed_vec_add', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: a_offset, a_offset_last, &
                                                abeg, aend, b_offset, bbeg, &
                                                bend, i, timing_handle
    LOGICAL                                  :: my_do_summation

!   ---------------------------------------------------------------------------
! this points to the last empty entry

    CALL timeset(routineN, timing_handle)

    my_do_summation = .TRUE.
    IF(PRESENT(do_summation)) my_do_summation=do_summation

    a_offset_last = INT(pkd_vec_a(n+1)) 
    IF(a_offset_last.LE.0) CALL stop_program(routineN,"a_offset_last.LT.0!!!!!!!")
    DO i = 1,n
       b_offset = INT(pkd_vec_b(i))
       a_offset = INT(pkd_vec_a(i))
       IF(b_offset.LE.0) CYCLE
       IF(a_offset.GT.0) THEN
          ! the block exsits in a, just add
          !CALL daxpy(vec_blk_size(i),1.0_dp,pkd_vec_b(b_offset),1,pkd_vec_a(a_offset),1)
          abeg = a_offset
          aend = a_offset + vec_blk_size(i) - 1
          bbeg = b_offset
          bend = b_offset + vec_blk_size(i) - 1
          IF(my_do_summation) THEN
             pkd_vec_a(abeg:aend) = pkd_vec_a(abeg:aend) + pkd_vec_b(bbeg:bend)
          ELSE
             pkd_vec_a(abeg:aend) = pkd_vec_b(bbeg:bend)
          ENDIF
       ELSE
          ! the block doesnt exsit in a, copy b at the end
          !CALL dcopy(vec_blk_size(i),pkd_vec_b(b_offset),1,pkd_vec_a(a_offset_last),1)
          abeg = a_offset_last
          aend = a_offset_last + vec_blk_size(i) - 1
          bbeg = b_offset
          bend = b_offset + vec_blk_size(i) - 1
          pkd_vec_a(i) = a_offset_last
          pkd_vec_a(abeg:aend) = pkd_vec_b(bbeg:bend)
          a_offset_last = a_offset_last + vec_blk_size(i)
       ENDIF
    ENDDO
    ! reset the last empty entry if needed
    pkd_vec_a(n+1) = REAL(a_offset_last,dp)
    CALL timestop(timing_handle)
  END SUBROUTINE packed_vec_add

! *****************************************************************************
!> \brief copy a vector (row or column) from a matrix to another matrix
!> \param[inout] matrix_b  matrix
!> \param[in]    matrix_a  matrix
!> \param[in]    what      row or column
!> \param[in]    ivec      the column/row index
!> 
! *****************************************************************************
  SUBROUTINE cp_dbcsr_copy_vec(matrix_a, matrix_b, what, ivec, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_b
    CHARACTER(len=*), INTENT(in)             :: what
    INTEGER, INTENT(IN)                      :: ivec
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_copy_vec', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: col, row, timing_handle
    LOGICAL                                  :: copy_column, tr
    REAL(KIND=dp), DIMENSION(:), POINTER     :: data_d
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)

    SELECT CASE(what)
    CASE('column')
       copy_column = .TRUE.
    CASE('row')
       copy_column = .FALSE.
    CASE DEFAULT
       CALL stop_program(routineN,"copy what?")
    END SELECT

    CALL dbcsr_iterator_start(iter, matrix_b)
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col, data_d, tr)
       IF( (copy_column.AND.col.EQ.ivec) .OR. (.NOT.copy_column.AND.row.EQ.ivec) ) &
            CALL dbcsr_put_block(matrix_a, row, col, data_d, tr, summation=.FALSE.)
       !IF(.NOT.copy_column.AND.row.EQ.ivec) &
       !     CALL dbcsr_put_block(matrix_a, row, col, data_d, tr,&
       !     summation=.FALSE.)
    ENDDO
    CALL dbcsr_iterator_stop(iter)
    CALL cp_dbcsr_finalize(matrix_a, error=error)

    CALL timestop(timing_handle)

  END SUBROUTINE cp_dbcsr_copy_vec


! *****************************************************************************
!> \brief Emulation of sparse_matrix_types/add_block_node mapped
!>        to add_real_matrix_block.
!>
!> It adds a block to the dbcsr matrix and returns a rank-2 pointer to the
!> block. Currently it only and always uses the mutable data.
!> \param[in,out] matrix      DBCSR matrix
!> \param[in]  row            the row
!> \param[in]  col            the column
!> \param[in]  block          the block to put
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE cp_dbcsr_add_block_node (matrix, block_row, block_col, block, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: block_row, block_col
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_add_block_node', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: c, ithread, p, r
    LOGICAL                                  :: dbg, existed, is_there, tr

!   ---------------------------------------------------------------------------
!    CALL cp_assert (matrix%m%work_mutable, cp_warning_level, cp_caller_error,&
!         routineN, "Mutable not defined upon DBCSR creation, forcing its use.")

    dbg = .FALSE.

    ithread = 0
!$  ithread = omp_get_thread_num()
    IF (.NOT. ASSOCIATED (matrix%m%wms)) THEN
       CALL cp_dbcsr_work_create (matrix, work_mutable=.TRUE., error=error)
       matrix%m%valid = .FALSE.
    ENDIF
!$  CALL cp_assert (SIZE (matrix%m%wms) .GE. omp_get_num_threads(),&
!$       cp_fatal_level, cp_wrong_args_error, routineN,&
!$       "Too few threads.", error=error)
    CALL cp_assert (dbcsr_wm_use_mutable (matrix%m%wms(ithread+1)),&
         cp_warning_level,&
         cp_unimplemented_error_nr, routineN,&
         "Data loss due to no conversion of appendable to mutable data")
    is_there = ASSOCIATED(block)
    !r = row ; c = col ; tr = .FALSE.
    !CALL dbcsr_get_stored_coordinates (matrix, r, c, tr)
    !CALL dbcsr_reserve_block2d (matrix, row, col, block)
    !write(*,*) 'add_block_node: block_row',block_row,' block_col',block_col
    CALL dbcsr_reserve_block2d (matrix, block_row, block_col, block,&
         existed=existed)
!
    IF (dbg) THEN
       r = block_row ; c = block_col ; tr = .FALSE.
       CALL dbcsr_get_stored_coordinates (matrix, r, c, tr, p)
       CALL cp_assert (p .EQ. dbcsr_mp_mynode (dbcsr_distribution_mp (&
            matrix%m%dist)),&
            cp_warning_level, cp_internal_error, routineN,&
            "Adding non-local element", error=error)
    ENDIF
    CALL cp_assert (.NOT.existed, cp_warning_level, cp_wrong_args_error,&
         routineN, "You should not add existing blocks according to old API.")
    IF(.NOT.is_there) block(:,:) = 0.0_dp
  END SUBROUTINE cp_dbcsr_add_block_node

! *****************************************************************************
!> \brief returns the id of the given matrix (from cp_sm_get_id_nr)
!> \param matrix the matrix you want info about
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
FUNCTION cp_dbcsr_get_id_nr(matrix,error) RESULT(res)
    TYPE(dbcsr_obj), POINTER                 :: matrix
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_get_id_nr', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
!  CPPrecondition(associated(matrix),cp_failure_level,routineP,error,failure)
  CPPrecondition(matrix%m%refcount>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     res=matrix%m%id_nr
  ELSE
     res=0
  END IF
END FUNCTION cp_dbcsr_get_id_nr

! *****************************************************************************
!> \brief Deallocates a DBCSR matrix for compatibility with CP2K
!> \param[in,out] matrix      DBCSR matrix
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE cp_dbcsr_deallocate_matrix(matrix, error)
    TYPE(dbcsr_obj), POINTER                 :: matrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_deallocate_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle

    CALL timeset(routineN, timing_handle)

    CALL cp_dbcsr_release (matrix, error=error)
    CALL cp_assert (.NOT. dbcsr_valid_index(matrix), cp_warning_level,&
         cp_caller_error, routineN,&
         'You should not "deallocate" a referenced matrix. '//&
         'Avoid pointers to DBCSR matrices.')
    DEALLOCATE (matrix)
    NULLIFY (matrix)

    CALL timestop(timing_handle)

  END SUBROUTINE cp_dbcsr_deallocate_matrix

! *****************************************************************************
! CP2k-compatible matrix sets
! *****************************************************************************


! *****************************************************************************
!> \brief   Allocate and initialize a real matrix 1-dimensional set.
!> \param[in,out] matrix_set  Set containing the DBCSR matrices
!> \param[in] nmatrix         Size of set
!> \param[in,out] error       cp2k error
!> \par History
!>      2009-08-17 Adapted from sparse_matrix_type for DBCSR
! *****************************************************************************
  SUBROUTINE allocate_dbcsr_matrix_set(matrix_set, nmatrix, error)
    TYPE(cp_dbcsr_obj_p_type), &
      DIMENSION(:), POINTER                  :: matrix_set
    INTEGER, INTENT(IN)                      :: nmatrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_dbcsr_matrix_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: imatrix, istat

    IF (ASSOCIATED(matrix_set)) CALL cp_dbcsr_deallocate_matrix_set(matrix_set,error=error)
    ALLOCATE (matrix_set(nmatrix),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"matrix_set",0)
    DO imatrix=1,nmatrix
       NULLIFY (matrix_set(imatrix)%matrix)
       !ALLOCATE (matrix_set(imatrix)%matrix)
       !CALL cp_dbcsr_init (matrix_set(imatrix)%matrix)
    END DO
  END SUBROUTINE allocate_dbcsr_matrix_set

! *****************************************************************************
!> \brief   Allocate and initialize a real matrix 2-dimensional set.
!> \param[in,out] matrix_set  Set containing the DBCSR matrix pointer type
!> \param[in] nmatrix         Size of set
!> \param[in,out] error       cp2k error
!> \par History
!>      2009-08-17 Adapted from sparse_matrix_type for DBCSR
! *****************************************************************************
  SUBROUTINE allocate_dbcsr_matrix_set_2d(matrix_set,nmatrix,mmatrix,error)
    TYPE(cp_dbcsr_obj_p_type), &
      DIMENSION(:, :), POINTER               :: matrix_set
    INTEGER, INTENT(IN)                      :: nmatrix, mmatrix
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_dbcsr_matrix_set_2d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: imatrix, istat, jmatrix

    IF (ASSOCIATED(matrix_set)) CALL cp_dbcsr_deallocate_matrix_set(matrix_set,error=error)
    ALLOCATE (matrix_set(nmatrix,mmatrix),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"matrix_set",0)
    DO jmatrix=1,mmatrix
      DO imatrix=1,nmatrix
         NULLIFY (matrix_set(imatrix,jmatrix)%matrix)
         !ALLOCATE (matrix_set(imatrix,jmatrix)%matrix)
         !CALL cp_dbcsr_init (matrix_set(imatrix,jmatrix)%matrix)
      END DO
    END DO
  END SUBROUTINE allocate_dbcsr_matrix_set_2d


! *****************************************************************************
!> \brief Deallocate a real matrix set and release all of the member matrices.
!> \param[in,out] matrix_set  Set containing the DBCSR matrix pointer type
!> \param[in,out] error       cp2k error
!> \par History
!>      2009-08-17 Adapted from sparse_matrix_type for DBCSR
! *****************************************************************************
  SUBROUTINE deallocate_dbcsr_matrix_set(matrix_set,error)

    TYPE(cp_dbcsr_obj_p_type), &
      DIMENSION(:), POINTER                  :: matrix_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_dbcsr_matrix_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: imatrix, istat, timing_handle

    CALL timeset(routineN, timing_handle)

    IF (ASSOCIATED(matrix_set)) THEN
      DO imatrix=1,SIZE(matrix_set)
        CALL cp_dbcsr_deallocate_matrix(matrix_set(imatrix)%matrix,error=error)
      END DO
      DEALLOCATE (matrix_set,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineP,"matrix_set")
    END IF

    CALL timestop(timing_handle)

  END SUBROUTINE deallocate_dbcsr_matrix_set

! *****************************************************************************
!> \brief Deallocate a real matrix set and release all of the member matrices.
!> \param[in,out] matrix_set  Set containing the DBCSR matrix pointer type
!> \param[in,out] error       cp2k error
!> \par History
!>      2009-08-17 Adapted from sparse_matrix_type for DBCSR
! *****************************************************************************
  SUBROUTINE deallocate_dbcsr_matrix_set_2d(matrix_set,error)

    TYPE(cp_dbcsr_obj_p_type), &
      DIMENSION(:, :), POINTER               :: matrix_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'deallocate_dbcsr_matrix_set_2d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: imatrix, istat, jmatrix

    IF (ASSOCIATED(matrix_set)) THEN
      DO jmatrix=1,SIZE(matrix_set,2)
        DO imatrix=1,SIZE(matrix_set,1)
          CALL cp_dbcsr_deallocate_matrix(matrix_set(imatrix,jmatrix)%matrix,error=error)
        END DO
      END DO
      DEALLOCATE (matrix_set,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineP,"matrix_set")
    END IF
  END SUBROUTINE deallocate_dbcsr_matrix_set_2d

! *****************************************************************************
!> \brief Maps sparse_matrix_type symmetry to DBCSR matrix type
!> \param[out] dbcsr_matrix_type        matrix_type of the DBCSR matrix
!> \param[in] sm_symmetry               sparse_matrix_symmetry
! *****************************************************************************
  SUBROUTINE cp_dbcsr_type_from_sm_symmetry (dbcsr_matrix_type, sm_symmetry)
    CHARACTER, INTENT(OUT)                   :: dbcsr_matrix_type
    CHARACTER(LEN=*), INTENT(IN)             :: sm_symmetry

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'cp_dbcsr_type_from_sm_symmetry', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    SELECT CASE (sm_symmetry)
    CASE ("symmetric")
       dbcsr_matrix_type = dbcsr_type_symmetric
    CASE ("antisymmetric")
       dbcsr_matrix_type = dbcsr_type_antisymmetric
    CASE ("none", "no symmetry")
       dbcsr_matrix_type = dbcsr_type_normal
    CASE default
       CALL cp_assert (.FALSE., cp_warning_level, cp_caller_error, routineN,&
            "Unknown matrix type "//sm_symmetry)
    END SELECT
  END SUBROUTINE cp_dbcsr_type_from_sm_symmetry


! *****************************************************************************
!> \brief Converts a cp2k full matrix into a DBCSR matrix.
!> \param[out] matrix         the created BCSR matrix
!> \param[in] fm              the cpk full matrix
!> \param[in] threshold       the threshold for determining sparsity
!> \param[in] row_blk_sizes   sizes of row blocks
!> \param[in] col_blk_sizes   sizes of column blocks
!> \param[in] distribution    the distribution to use for the new matrix
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE cp_dbcsr_from_fm(matrix, fm, threshold, distribution, row_blk_size,&
       col_blk_size, error)
    TYPE(dbcsr_obj), INTENT(OUT)             :: matrix
    TYPE(cp_fm_type), POINTER                :: fm
    REAL(KIND=dp), INTENT(IN)                :: threshold
    TYPE(dbcsr_distribution_obj)             :: distribution
    TYPE(array_i1d_obj), INTENT(IN)          :: row_blk_size, col_blk_size
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_from_fm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle

    CALL timeset(routineN, timing_handle)

    !CALL cp_dbcsr_init (matrix, error)! the matrix should already be initialized
    CALL cp_dbcsr_create(matrix, fm%name, distribution, dbcsr_type_normal,&
         row_blk_size, col_blk_size,&
         0, 0, dbcsr_type_real_8, error=error)
    CALL cp_dbcsr_finalize(matrix, error=error)
    CALL copy_fm_to_dbcsr(fm, matrix, error=error)
    CALL cp_dbcsr_verify_matrix(matrix, error)
    CALL timestop(timing_handle)

  END SUBROUTINE cp_dbcsr_from_fm


! *****************************************************************************
!> \brief   Copy a BLACS matrix to a dbcsr matrix.
!>
!>          real_matrix=beta*real_matrix+alpha*fm
!>          beta defaults to 0, alpha to 1
!> \author  Urban Borstnik
!> \date    2009-10-13
!> \par History
!>          2009-10-13 rewritten based on copy_dbcsr_to_fm
!> \version 2.0
!> \param[in] fm              full matrix
!> \param[out] matrix         DBCSR matrix
!> \param[in] alpha           (optional) scaling of FM
!> \param[in] beta            (optional) scaling of existing SM
!> \param[in] keep_sparsity   (optional) retains the sparsity of the input
!>                            matrix
! *****************************************************************************
  SUBROUTINE copy_fm_to_dbcsr(fm,matrix,alpha,beta,keep_sparsity,error)
    TYPE(cp_fm_type), POINTER                :: fm
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    REAL(kind=dp), INTENT(IN), OPTIONAL      :: alpha, beta
    LOGICAL, INTENT(IN), OPTIONAL            :: keep_sparsity
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'copy_fm_to_dbcsr', &
      routineP = moduleN//':'//routineN

    INTEGER :: blk_p, col, col_l, col_size, error_handler, group, handle, &
      nblkcols_local, nblkcols_total, nblkrows_local, nblkrows_total, &
      ncol_block, ncol_global, nfullcols_local, nfullcols_total, &
      nfullrows_local, nfullrows_total, nrow_block, nrow_global, nze, row, &
      row_l, row_size
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_col, first_row, &
                                                last_col, last_row, &
                                                local_col_sizes, &
                                                local_row_sizes
    INTEGER, DIMENSION(:), POINTER           :: cbs, local_cols, local_rows, &
                                                rbs
    REAL(kind=dp)                            :: my_beta
    REAL(KIND=dp), DIMENSION(:), POINTER     :: blk_1d_dp
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: fm_block
    REAL(KIND=sp), DIMENSION(:), POINTER     :: blk_1d_sp
    REAL(kind=sp), DIMENSION(:, :), POINTER  :: fm_block_sp
    TYPE(array_i1d_obj)                      :: col_blk_size, row_blk_size
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(dbcsr_distribution_obj)             :: bc_dist
    TYPE(dbcsr_error_type)                   :: dbcsr_error
    TYPE(dbcsr_obj)                          :: bc_mat

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL dbcsr_error_set(routineN, error_handler, dbcsr_error)
    CALL dbcsr_access_flush (matrix, error=dbcsr_error)
    CALL timeset(routineN,handle)

    my_beta=0._dp
    IF (PRESENT(beta)) THEN
       CALL cp_assert (beta .EQ. my_beta, cp_fatal_level,&
            cp_unimplemented_error_nr, routineN,&
            "beta not supported, use matrix addition instead")
       my_beta=beta
    ENDIF
    group = fm%matrix_struct%para_env%group
    context => fm%matrix_struct%context
    nrow_block = fm%matrix_struct%nrow_block
    ncol_block = fm%matrix_struct%ncol_block
    nrow_global = fm%matrix_struct%nrow_global
    ncol_global = fm%matrix_struct%ncol_global
    
    CALL dbcsr_get_info(matrix,&
         nfullrows_total=nfullrows_total,&
         nfullcols_total=nfullcols_total)

    CALL cp_assert (nrow_global.eq.nfullrows_total,&
         cp_fatal_level, cp_caller_error, routineN,&
         "FM and DBCSR matrix sizes do not match in rows")
    CALL cp_assert (ncol_global.eq.nfullcols_total,&
         cp_fatal_level, cp_caller_error, routineN,&
         "FM and DBCSR matrix sizes do not match in columns")

    ! Create a block-cyclic distribution compatible with the FM matrix.
    CALL dbcsr_distribution_init (bc_dist)
    CALL dbcsr_create_dist_block_cyclic (bc_dist,&
         nfullrows_total, nfullcols_total,& ! Actual full matrix size
         nrow_block, ncol_block,&           ! BLACS parameters
         dbcsr_distribution_mp (dbcsr_distribution (matrix)),&
         row_blk_size, col_blk_size)        ! block-cyclic row/col sizes

    ! Create the block-cyclic DBCSR matrix
    CALL cp_dbcsr_init (bc_mat, error)
    CALL cp_dbcsr_create (bc_mat, "Block-cyclic "//matrix%m%name, bc_dist,&
         dbcsr_get_matrix_type(matrix), row_blk_size, col_blk_size, 0, 0,&
         data_type=dbcsr_get_data_type(matrix),error=error)

    !call dbcsr_finalize (bc_mat)
    CALL dbcsr_distribution_release (bc_dist)
    CALL array_release (row_blk_size)
    CALL array_release (col_blk_size)

    CALL dbcsr_get_info(bc_mat,&
         nblkrows_total=nblkrows_total,&
         nblkcols_total=nblkcols_total,&
         nblkrows_local=nblkrows_local,&
         nblkcols_local=nblkcols_local,&
         nfullrows_local=nfullrows_local,&
         nfullcols_local=nfullcols_local,&
         nfullrows_total=nfullrows_total,&
         nfullcols_total=nfullcols_total,&
         local_rows=local_rows,&
         local_cols=local_cols,&
         row_blk_size=row_blk_size,&
         col_blk_size=col_blk_size)

    !WRITE(*,*)routineN//" sizes",nblkrows_total,&
    !     nblkcols_total,&
    !     nblkrows_local,&
    !     nblkcols_local,&
    !     nfullrows_local,&
    !     nfullcols_local,&
    !     nfullrows_total,&
    !     nfullcols_total


    rbs => array_data (row_blk_size)
    cbs => array_data (col_blk_size)
    ALLOCATE (local_row_sizes (nblkrows_total))
    ALLOCATE (local_col_sizes (nblkcols_total))
    local_row_sizes(:) = 0
    IF (nblkrows_local .GE. 1) THEN
       FORALL (row = 1 : nblkrows_local)
          local_row_sizes(local_rows(row)) = rbs(local_rows(row))
       END FORALL
    ENDIF
    IF (dbg) WRITE(*,*)routineN//" local row sizes",local_row_sizes
    local_col_sizes(:) = 0
    IF (nblkcols_local .GE. 1) THEN
       FORALL (col = 1 : nblkcols_local)
          local_col_sizes(local_cols(col)) = cbs(local_cols(col))
       END FORALL
    ENDIF
    IF (dbg) WRITE(*,*)routineN//" local col sizes",local_col_sizes

    ALLOCATE (first_row(nblkrows_total),last_row(nblkrows_total))
    ALLOCATE (first_col(nblkcols_total),last_col(nblkcols_total))
    CALL convert_sizes_to_offsets (local_row_sizes, first_row, last_row)
    CALL convert_sizes_to_offsets (local_col_sizes, first_col, last_col)

    ! Copy the FM data to the block-cyclic DBCSR matrix.  This step
    ! could be skipped with appropriate DBCSR index manipulation.
    fm_block => fm%local_data
    fm_block_sp => fm%local_data_sp

    CALL cp_dbcsr_work_create (bc_mat, nblks_guess=nblkrows_local*nblkcols_local,&
         sizedata_guess=nfullrows_local*nfullcols_local, work_mutable=.FALSE.,&
         n=1, error=error)
    blk_p = 1
    bc_rows: DO row_l = 1, nblkrows_local
       row = local_rows (row_l)
       row_size = rbs(row)
       bc_cols: DO col_l = 1, nblkcols_local
          col = local_cols (col_l)
          col_size = cbs(col)
          nze = row_size*col_size
          !WRITE(*,*)routineN//" Adding block",row,col,"size",nze
          CALL add_work_coordinate(bc_mat%m%wms(1), row, col, blk_p, error=dbcsr_error)
          IF (fm%use_sp) THEN
             !blk_1d_sp => bc_mat%m%wms(1)%data_area%d%r_sp(blk_p:blk_p+nze-1)
             blk_1d_sp => dbcsr_get_data_p (bc_mat%m%wms(1)%data_area,&
                  coersion=REAL(0.0, KIND=sp), lb=blk_p, ub=blk_p+nze-1)
          ELSE
             !blk_1d_dp => bc_mat%m%wms(1)%data_area%d%r_dp(blk_p:blk_p+nze-1)
             blk_1d_dp => dbcsr_get_data_p (bc_mat%m%wms(1)%data_area,&
                  coersion=REAL(0.0, KIND=dp), lb=blk_p, ub=blk_p+nze-1)
          ENDIF
          CALL cp_assert (nze .EQ. (last_row(row)-first_row(row)+1)*(last_col(col)-first_col(col)+1),&
               cp_fatal_level, cp_internal_error, routineN,&
               "Block size does not match block row/col sizes")
          IF (fm%use_sp) THEN
             blk_1d_sp(1:nze) = RESHAPE(&
                  fm_block_sp(&
                     first_row(row):last_row(row),first_col(col):last_col(col)&
                  ), (/ nze /))
          ELSE
             blk_1d_dp(1:nze) = RESHAPE(&
                  fm_block(&
                     first_row(row):last_row(row),first_col(col):last_col(col)&
                  ), (/ nze /))
          ENDIF
          blk_p = blk_p + nze
       ENDDO bc_cols
    ENDDO bc_rows
    bc_mat%m%wms(1)%datasize = blk_p - 1
    CALL cp_dbcsr_finalize (bc_mat, reshuffle=.FALSE., error=error)

    ! Now convert to the desired matrix distribution
    IF (PRESENT (alpha)) THEN
       CALL stop_program(routineN,'no more alpha... cleanme')
    ELSE
       CALL dbcsr_complete_redistribute (bc_mat, matrix,&
            keep_sparsity=keep_sparsity, error=dbcsr_error)
    ENDIF
    CALL cp_dbcsr_release (bc_mat, error=error)

    CALL timestop(handle)
    CALL dbcsr_error_stop(error_handler, dbcsr_error)
  END SUBROUTINE copy_fm_to_dbcsr

! *****************************************************************************
!> \brief   Copy a BLACS matrix to a dbcsr matrix.
!>
!>          real_matrix=beta*real_matrix+alpha*fm
!>          beta defaults to 0, alpha to 1
!> \author  VW
!> \date    2010
!> \par History
!>          2010  copied from copy_dbcsr_to_fm
!> \version 2.0
!> \param[in] fm              full matrix
!> \param[out] matrix         DBCSR matrix
!> \param[in] alpha           (optional) scaling of FM
!> \param[in] beta            (optional) scaling of existing SM
!> \param[in] keep_sparsity   (optional) retains the sparsity of the input
!>                            matrix
! *****************************************************************************
  SUBROUTINE copy_cfm_to_dbcsr(fm,matrix,keep_sparsity,error)
    TYPE(cp_cfm_type), POINTER               :: fm
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    LOGICAL, INTENT(IN), OPTIONAL            :: keep_sparsity
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'copy_cfm_to_dbcsr', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: blk_1d
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: fm_block
    INTEGER :: blk_p, col, col_l, col_size, error_handler, group, handle, &
      nblkcols_local, nblkcols_total, nblkrows_local, nblkrows_total, &
      ncol_block, ncol_global, nfullcols_local, nfullcols_total, &
      nfullrows_local, nfullrows_total, nrow_block, nrow_global, nze, row, &
      row_l, row_size
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_col, first_row, &
                                                last_col, last_row, &
                                                local_col_sizes, &
                                                local_row_sizes
    INTEGER, DIMENSION(:), POINTER           :: cbs, local_cols, local_rows, &
                                                rbs
    TYPE(array_i1d_obj)                      :: col_blk_size, row_blk_size
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(dbcsr_distribution_obj)             :: bc_dist
    TYPE(dbcsr_error_type)                   :: dbcsr_error
    TYPE(dbcsr_obj)                          :: bc_mat

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL dbcsr_error_set(routineN, error_handler, dbcsr_error)
    CALL dbcsr_access_flush (matrix, error=dbcsr_error)
    CALL timeset(routineN,handle)

    group = fm%matrix_struct%para_env%group
    context => fm%matrix_struct%context
    nrow_block = fm%matrix_struct%nrow_block
    ncol_block = fm%matrix_struct%ncol_block
    nrow_global = fm%matrix_struct%nrow_global
    ncol_global = fm%matrix_struct%ncol_global
    
    CALL dbcsr_get_info(matrix,&
         nfullrows_total=nfullrows_total,&
         nfullcols_total=nfullcols_total)

    CALL cp_assert (nrow_global.eq.nfullrows_total,&
         cp_fatal_level, cp_caller_error, routineN,&
         "FM and DBCSR matrix sizes do not match in rows")
    CALL cp_assert (ncol_global.eq.nfullcols_total,&
         cp_fatal_level, cp_caller_error, routineN,&
         "FM and DBCSR matrix sizes do not match in columns")

    ! Create a block-cyclic distribution compatible with the FM matrix.
    CALL dbcsr_distribution_init (bc_dist)
    CALL dbcsr_create_dist_block_cyclic (bc_dist,&
         nfullrows_total, nfullcols_total,& ! Actual full matrix size
         nrow_block, ncol_block,&           ! BLACS parameters
         dbcsr_distribution_mp (dbcsr_distribution (matrix)),&
         row_blk_size, col_blk_size)        ! block-cyclic row/col sizes

    ! Create the block-cyclic DBCSR matrix
    CALL cp_dbcsr_init (bc_mat, error)
    CALL cp_dbcsr_create (bc_mat, "Block-cyclic "//matrix%m%name, bc_dist,&
         dbcsr_get_matrix_type(matrix), row_blk_size, col_blk_size, 0, 0,&
         dbcsr_type_complex_8,error=error) ! type hard coded !
    !call dbcsr_finalize (bc_mat)
    CALL dbcsr_distribution_release (bc_dist)
    CALL array_release (row_blk_size)
    CALL array_release (col_blk_size)

    CALL dbcsr_get_info(bc_mat,&
         nblkrows_total=nblkrows_total,&
         nblkcols_total=nblkcols_total,&
         nblkrows_local=nblkrows_local,&
         nblkcols_local=nblkcols_local,&
         nfullrows_local=nfullrows_local,&
         nfullcols_local=nfullcols_local,&
         nfullrows_total=nfullrows_total,&
         nfullcols_total=nfullcols_total,&
         local_rows=local_rows,&
         local_cols=local_cols,&
         row_blk_size=row_blk_size,&
         col_blk_size=col_blk_size)

    rbs => array_data (row_blk_size)
    cbs => array_data (col_blk_size)
    ALLOCATE (local_row_sizes (nblkrows_total))
    ALLOCATE (local_col_sizes (nblkcols_total))
    local_row_sizes(:) = 0
    IF (nblkrows_local .GE. 1) THEN
       FORALL (row = 1 : nblkrows_local)
          local_row_sizes(local_rows(row)) = rbs(local_rows(row))
       END FORALL
    ENDIF
    IF (dbg) WRITE(*,*)routineN//" local row sizes",local_row_sizes
    local_col_sizes(:) = 0
    IF (nblkcols_local .GE. 1) THEN
       FORALL (col = 1 : nblkcols_local)
          local_col_sizes(local_cols(col)) = cbs(local_cols(col))
       END FORALL
    ENDIF
    IF (dbg) WRITE(*,*)routineN//" local col sizes",local_col_sizes

    ALLOCATE (first_row(nblkrows_total),last_row(nblkrows_total))
    ALLOCATE (first_col(nblkcols_total),last_col(nblkcols_total))
    CALL convert_sizes_to_offsets (local_row_sizes, first_row, last_row)
    CALL convert_sizes_to_offsets (local_col_sizes, first_col, last_col)

    ! Copy the FM data to the block-cyclic DBCSR matrix.  This step
    ! could be skipped with appropriate DBCSR index manipulation.
    fm_block => fm%local_data
    CALL cp_dbcsr_work_create (bc_mat, nblks_guess=nblkrows_local*nblkcols_local,&
         sizedata_guess=nfullrows_local*nfullcols_local, work_mutable=.FALSE.,&
         n=1, error=error)
    blk_p = 1
    bc_rows: DO row_l = 1, nblkrows_local
       row = local_rows (row_l)
       row_size = rbs(row)
       bc_cols: DO col_l = 1, nblkcols_local
          col = local_cols (col_l)
          col_size = cbs(col)
          nze = row_size*col_size
          !WRITE(*,*)routineN//" Adding block",row,col,"size",nze
          CALL add_work_coordinate(bc_mat%m%wms(1), row, col, blk_p, error=dbcsr_error)
          !blk_1d => bc_mat%m%wms(1)%data_area%d%c_dp(blk_p:blk_p+nze-1)
          blk_1d => dbcsr_get_data_p (bc_mat%m%wms(1)%data_area,&
               coersion=CMPLX(0.0, KIND=dp), lb=blk_p, ub=blk_p+nze-1)
          CALL cp_assert (nze .EQ. (last_row(row)-first_row(row)+1)*(last_col(col)-first_col(col)+1),&
               cp_fatal_level, cp_internal_error, routineN,&
               "Block size does not match block row/col sizes")
          blk_1d(1:nze) = RESHAPE(&
               fm_block(&
               first_row(row):last_row(row),first_col(col):last_col(col)&
               ), (/ nze /))
          blk_p = blk_p + nze
       ENDDO bc_cols
    ENDDO bc_rows
    bc_mat%m%wms(1)%datasize = blk_p - 1
    CALL cp_dbcsr_finalize (bc_mat, reshuffle=.FALSE., error=error)

    ! Now convert to the desired matrix distribution
    CALL dbcsr_complete_redistribute (bc_mat, matrix, keep_sparsity=keep_sparsity,&
         error=dbcsr_error)
    CALL cp_dbcsr_release (bc_mat, error=error)

    CALL timestop(handle)
    CALL dbcsr_error_stop(error_handler, dbcsr_error)
  END SUBROUTINE copy_cfm_to_dbcsr

! *****************************************************************************
!> \brief   Copy a DBCSR matrix to a BLACS matrix
!> \param[in] matrix          DBCSR matrix
!> \param[out] fm             full matrix
! *****************************************************************************
  SUBROUTINE copy_dbcsr_to_fm(matrix, fm, error)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    TYPE(cp_fm_type), POINTER                :: fm
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'copy_dbcsr_to_fm', &
      routineP = moduleN//':'//routineN

    INTEGER :: col, error_handle, group, handle, mypcol, mype, myprow, &
      nblkcols_local, nblkcols_total, nblkrows_local, nblkrows_total, &
      ncol_block, ncol_global, nfullcols_total, nfullrows_total, npcol, npe, &
      nprow, nrow_block, nrow_global, row
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_col, first_row, &
                                                last_col, last_row, &
                                                local_col_sizes, &
                                                local_row_sizes
    INTEGER, DIMENSION(:), POINTER           :: col_blk_sizes, local_cols, &
                                                local_rows, row_blk_sizes
    LOGICAL                                  :: tr
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: blk_2d, fm_block
    REAL(KIND=sp), DIMENSION(:, :), POINTER  :: fm_block_sp
    TYPE(array_i1d_obj)                      :: cbs, rbs
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(dbcsr_distribution_obj)             :: bc_dist
    TYPE(dbcsr_error_type)                   :: dbcsr_error
    TYPE(dbcsr_iterator)                     :: iter
    TYPE(dbcsr_obj)                          :: bc_mat

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY

    CALL timeset(routineN,handle)
    CALL dbcsr_error_set (routineN, error_handle, error=dbcsr_error)
    CALL dbcsr_access_flush (matrix, error=dbcsr_error)

    ! info about the full matrix
    group = fm%matrix_struct%para_env%group
    context => fm%matrix_struct%context
    mype=context%my_pid
    npe=context%n_pid
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)
    nrow_block = fm%matrix_struct%nrow_block
    ncol_block = fm%matrix_struct%ncol_block
    nrow_global = fm%matrix_struct%nrow_global
    ncol_global = fm%matrix_struct%ncol_global

    CALL dbcsr_get_info(matrix,&
         nfullrows_total=nfullrows_total,&
         nfullcols_total=nfullcols_total)

    ! Convert DBCSR to a block-cyclic one
    CALL dbcsr_distribution_init (bc_dist)
    CALL dbcsr_create_dist_block_cyclic (bc_dist,&
         nfullrows_total, nfullcols_total,&
         nrow_block, ncol_block,&
         dbcsr_distribution_mp (dbcsr_distribution(matrix)),&
         rbs, cbs)

    CALL cp_dbcsr_init (bc_mat, error)
    CALL cp_dbcsr_create (bc_mat, "Block-cyclic"//matrix%m%name, bc_dist,&
         dbcsr_type_normal, rbs, cbs, 0, 0, error=error)
    CALL cp_dbcsr_finalize (bc_mat, error=error)
    CALL dbcsr_distribution_release (bc_dist)
    CALL array_release (rbs)
    CALL array_release (cbs)
    CALL dbcsr_complete_redistribute (matrix, bc_mat, error=dbcsr_error)

    ! Find the local extents of the local blocked rows so that index lookups
    ! into the FM matrix work correctly.
    row_blk_sizes => array_data (rbs)
    col_blk_sizes => array_data (cbs)
    local_rows => array_data (dbcsr_distribution_local_rows (bc_dist))
    local_cols => array_data (dbcsr_distribution_local_cols (bc_dist))
    ALLOCATE (local_row_sizes (dbcsr_distribution_nrows (bc_dist)))
    ALLOCATE (local_col_sizes (dbcsr_distribution_ncols (bc_dist)))
    nblkrows_local = dbcsr_distribution_nlocal_rows (bc_dist)
    nblkcols_local = dbcsr_distribution_nlocal_cols (bc_dist)
    local_row_sizes(:) = 0
    IF (nblkrows_local .GE. 1) THEN
       FORALL (row = 1 : nblkrows_local)
          local_row_sizes(local_rows(row)) = row_blk_sizes(local_rows(row))
       END FORALL
    ENDIF
    IF (dbg) WRITE(*,*)routineN//" local row sizes",local_row_sizes
    local_col_sizes(:) = 0
    IF (nblkcols_local .GE. 1) THEN
       FORALL (col = 1 : nblkcols_local)
          local_col_sizes(local_cols(col)) = col_blk_sizes(local_cols(col))
       END FORALL
    ENDIF
    IF (dbg) WRITE(*,*)routineN//" local col sizes",local_col_sizes
    nblkrows_total = dbcsr_distribution_nrows (bc_dist)
    nblkcols_total = dbcsr_distribution_ncols (bc_dist)
    ALLOCATE (first_row(nblkrows_total),last_row(nblkrows_total))
    ALLOCATE (first_col(nblkcols_total),last_col(nblkcols_total))
    CALL convert_sizes_to_offsets (local_row_sizes, first_row, last_row)
    CALL convert_sizes_to_offsets (local_col_sizes, first_col, last_col)
    !
    ! Now copy data to the FM matrix
    fm_block => fm%local_data
    fm_block_sp => fm%local_data_sp
    IF(fm%use_sp) THEN
       fm_block_sp=0.0_sp
    ELSE
       fm_block=0.0_dp
    ENDIF

    IF(dbg) THEN
       WRITE(*,*)routineN//" FM data size is", UBOUND(fm_block)
       WRITE(*,*)routineN//" dbcsr data size is", dbcsr_get_data_size(bc_mat)
       WRITE(*,*)routineN//" FM block sizes are",nrow_block,'/',nfullrows_total
       WRITE(*,*)routineN//" FM block sizes are",ncol_block,'/',nfullcols_total
       WRITE(*,*)routineN//" dbcsr row sizes are",bc_mat%m%row_blk_size%low%data
       WRITE(*,*)routineN//" dbcsr col sizes are",bc_mat%m%col_blk_size%low%data
    ENDIF
    !
    CALL dbcsr_iterator_start(iter, bc_mat)
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col, blk_2d, tr)
       ! Convert absolute coordinates to FM-local coordinates
       IF (.NOT. tr) THEN

          IF(fm%use_sp) THEN
             fm_block_sp(first_row(row):last_row(row),first_col(col):last_col(col))&
                  = REAL(blk_2d(:,:),sp)
          ELSE
             IF (dbg) THEN
                WRITE(*,*)routineN//" blk2d size",UBOUND(blk_2d)
                WRITE(*,*)routineN//" want to set coor.",row,col
                WRITE(*,*)routineN//" local extents",&
                     first_row(row),last_row(row),first_col(col),last_col(col)
             ENDIF
             fm_block(first_row(row):last_row(row),first_col(col):last_col(col))&
                  = blk_2d(:,:)
          ENDIF
       ELSE
          IF(fm%use_sp) THEN
             fm_block_sp(first_col(col):last_col(col),first_row(row):last_row(row))&
                  = TRANSPOSE(REAL(blk_2d(:,:),sp))
          ELSE
             WRITE(*,*)'transposing...'

             fm_block(first_col(col):last_col(col),first_row(row):last_row(row))&
                  = TRANSPOSE(blk_2d(:,:))
          ENDIF
       ENDIF
    ENDDO
    CALL dbcsr_iterator_stop(iter)

    CALL cp_dbcsr_release (bc_mat, error=error)
    CALL dbcsr_error_stop (error_handle, error=dbcsr_error)
    CALL timestop(handle)
  END SUBROUTINE copy_dbcsr_to_fm

! *****************************************************************************
!> \brief   Copy a DBCSR matrix to a BLACS matrix
!> \param[in] matrix          DBCSR matrix
!> \param[out] fm             full matrix
! *****************************************************************************
  SUBROUTINE copy_dbcsr_to_cfm(matrix, fm, error)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    TYPE(cp_cfm_type), POINTER               :: fm
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'copy_dbcsr_to_cfm', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: blk_2d, fm_block
    INTEGER :: col, error_handle, group, handle, mypcol, mype, myprow, &
      nblkcols_local, nblkcols_total, nblkrows_local, nblkrows_total, &
      ncol_block, ncol_global, nfullcols_total, nfullrows_total, npcol, npe, &
      nprow, nrow_block, nrow_global, row
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_col, first_row, &
                                                last_col, last_row, &
                                                local_col_sizes, &
                                                local_row_sizes
    INTEGER, DIMENSION(:), POINTER           :: col_blk_sizes, local_cols, &
                                                local_rows, row_blk_sizes
    LOGICAL                                  :: tr
    TYPE(array_i1d_obj)                      :: cbs, rbs
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(dbcsr_distribution_obj)             :: bc_dist
    TYPE(dbcsr_error_type)                   :: dbcsr_error
    TYPE(dbcsr_iterator)                     :: iter
    TYPE(dbcsr_obj)                          :: bc_mat

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY

    CALL timeset(routineN,handle)
    CALL dbcsr_error_set (routineN, error_handle, error=dbcsr_error)
    CALL dbcsr_access_flush (matrix, error=dbcsr_error)

    ! info about the full matrix
    group = fm%matrix_struct%para_env%group
    context => fm%matrix_struct%context
    mype=context%my_pid
    npe=context%n_pid
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)
    nrow_block = fm%matrix_struct%nrow_block
    ncol_block = fm%matrix_struct%ncol_block
    nrow_global = fm%matrix_struct%nrow_global
    ncol_global = fm%matrix_struct%ncol_global

    CALL dbcsr_get_info(matrix,&
         nfullrows_total=nfullrows_total,&
         nfullcols_total=nfullcols_total)

    ! Convert DBCSR to a block-cyclic one
    CALL dbcsr_distribution_init (bc_dist)
    CALL dbcsr_create_dist_block_cyclic (bc_dist,&
         nfullrows_total, nfullcols_total,&
         nrow_block, ncol_block,&
         dbcsr_distribution_mp (dbcsr_distribution(matrix)),&
         rbs, cbs)

    CALL cp_dbcsr_init (bc_mat, error)
    CALL cp_dbcsr_create (bc_mat, "Block-cyclic"//matrix%m%name, bc_dist,&
         dbcsr_type_normal, rbs, cbs, 0, 0, dbcsr_get_data_type(matrix),&
         error=error)
    CALL cp_dbcsr_finalize (bc_mat, error=error)
    CALL dbcsr_distribution_release (bc_dist)
    CALL array_release (rbs)
    CALL array_release (cbs)
    CALL dbcsr_complete_redistribute (matrix, bc_mat, error=dbcsr_error)

    ! Find the local extents of the local blocked rows so that index lookups
    ! into the FM matrix work correctly.
    row_blk_sizes => array_data (rbs)
    col_blk_sizes => array_data (cbs)
    local_rows => array_data (dbcsr_distribution_local_rows (bc_dist))
    local_cols => array_data (dbcsr_distribution_local_cols (bc_dist))
    ALLOCATE (local_row_sizes (dbcsr_distribution_nrows (bc_dist)))
    ALLOCATE (local_col_sizes (dbcsr_distribution_ncols (bc_dist)))
    nblkrows_local = dbcsr_distribution_nlocal_rows (bc_dist)
    nblkcols_local = dbcsr_distribution_nlocal_cols (bc_dist)
    local_row_sizes(:) = 0
    IF (nblkrows_local .GE. 1) THEN
       FORALL (row = 1 : nblkrows_local)
          local_row_sizes(local_rows(row)) = row_blk_sizes(local_rows(row))
       END FORALL
    ENDIF
    IF (dbg) WRITE(*,*)routineN//" local row sizes",local_row_sizes
    local_col_sizes(:) = 0
    IF (nblkcols_local .GE. 1) THEN
       FORALL (col = 1 : nblkcols_local)
          local_col_sizes(local_cols(col)) = col_blk_sizes(local_cols(col))
       END FORALL
    ENDIF
    IF (dbg) WRITE(*,*)routineN//" local col sizes",local_col_sizes
    nblkrows_total = dbcsr_distribution_nrows (bc_dist)
    nblkcols_total = dbcsr_distribution_ncols (bc_dist)
    ALLOCATE (first_row(nblkrows_total),last_row(nblkrows_total))
    ALLOCATE (first_col(nblkcols_total),last_col(nblkcols_total))
    CALL convert_sizes_to_offsets (local_row_sizes, first_row, last_row)
    CALL convert_sizes_to_offsets (local_col_sizes, first_col, last_col)
    !
    ! Now copy data to the FM matrix
    fm_block => fm%local_data
    fm_block=(0.0_dp,0.0_dp)

    IF(dbg) THEN
       WRITE(*,*)routineN//" FM data size is", UBOUND(fm_block)
       WRITE(*,*)routineN//" dbcsr data size is", dbcsr_get_data_size(bc_mat)
       WRITE(*,*)routineN//" FM block sizes are",nrow_block,'/',nfullrows_total
       WRITE(*,*)routineN//" FM block sizes are",ncol_block,'/',nfullcols_total
       WRITE(*,*)routineN//" dbcsr row sizes are",bc_mat%m%row_blk_size%low%data
       WRITE(*,*)routineN//" dbcsr col sizes are",bc_mat%m%col_blk_size%low%data
    ENDIF
    !
    CALL dbcsr_iterator_start(iter, bc_mat)
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col, blk_2d, tr)
       ! Convert absolute coordinates to FM-local coordinates
       IF (.NOT. tr) THEN
          IF (dbg) THEN
             WRITE(*,*)routineN//" blk2d size",UBOUND(blk_2d)
             WRITE(*,*)routineN//" want to set coor.",row,col
             WRITE(*,*)routineN//" local extents",&
                  first_row(row),last_row(row),first_col(col),last_col(col)
          ENDIF
          fm_block(first_row(row):last_row(row),first_col(col):last_col(col))&
               = blk_2d(:,:)
       ELSE
          WRITE(*,*)'transposing...'
          fm_block(first_col(col):last_col(col),first_row(row):last_row(row))&
               = TRANSPOSE(CMPLX(blk_2d(:,:)))
       ENDIF
    ENDDO
    CALL dbcsr_iterator_stop(iter)

    CALL cp_dbcsr_release (bc_mat, error=error)
    CALL dbcsr_error_stop (error_handle, error=dbcsr_error)
    CALL timestop(handle)
  END SUBROUTINE copy_dbcsr_to_cfm

! *****************************************************************************
!> \brief Create a DBCSR matrix from a real_matrix_type
!> \par Limitations
!>      Data is never moved among processes, so it is more efficient
!>      and gracefully handles any misuse of distribution_2d.
!> \param[out] matrix         bcsr matrix
!> \param[in] sm    cp2k sparse matrix
!> \param[in,out] error     cp2k error
! *****************************************************************************
  SUBROUTINE cp_dbcsr_from_sm(matrix, sm, error, distribution_2d, mp_obj)
    TYPE(dbcsr_obj), INTENT(OUT)             :: matrix
    TYPE(real_matrix_type), POINTER          :: sm
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    TYPE(distribution_2d_type), OPTIONAL, &
      POINTER                                :: distribution_2d
    TYPE(dbcsr_mp_obj), INTENT(IN), OPTIONAL :: mp_obj

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_from_sm', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: data_exchange_info_cnt = 3, &
                                                metalen = 2

    INTEGER :: absrow, blk, block_col, mypcol, myprow, nblkcols_local, &
      nblkcols_total, nblkrows_local, nblkrows_total, nfullcols_local, &
      nfullcols_total, nfullrows_local, nfullrows_total, npcols, nprows, nze, &
      row_l, stat, timing_handle
    INTEGER, DIMENSION(:), POINTER :: blkcols_local, blkrows_local, &
      col_blk_begin, col_blk_end, col_blk_size, pcol_dist, prow_dist, &
      row_blk_begin, row_blk_end, row_blk_size
    INTEGER, DIMENSION(:, :), POINTER        :: blacs2mpi
    LOGICAL                                  :: symmetric
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: blkdata_2d
    TYPE(array_i1d_obj)                      :: cbs_array, rbs_array
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dbcsr_distribution_obj)             :: dist
    TYPE(distribution_2d_type), POINTER      :: dist2d
    TYPE(real_block_node_type), POINTER      :: sm_block_node

!hack for c**py kg*.F

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
    IF(PRESENT(distribution_2d)) THEN
       dist2d => distribution_2d
    ELSE
       dist2d => sm%distribution_2d
    ENDIF
    CALL dbcsr_distribution_init (dist)
    CALL cp_dbcsr_dist2d_to_dist (dist2d, dist, mp_obj=mp_obj, error=error)
    symmetric = (sm%symmetry .EQ. "symmetric")
    nblkrows_total = sm%nblock_row
    nblkcols_total = sm%nblock_col
    row_blk_begin => sm%first_row
    row_blk_end => sm%last_row
    col_blk_begin => sm%first_col
    col_blk_end => sm%last_col
    CALL distribution_2d_get(dist2d, row_distribution=prow_dist,&
         col_distribution=pcol_dist, error=error, blacs_env=blacs_env,&
         n_flat_local_rows=nblkrows_local, n_flat_local_cols=nblkcols_local,&
         n_row_distribution=nblkrows_total, n_col_distribution=nblkcols_total,&
         flat_local_rows=blkrows_local, flat_local_cols=blkcols_local)
    nfullrows_total = MAXVAL(row_blk_end)
    nfullcols_total = MAXVAL(col_blk_end)
    nfullrows_local = 0
    nfullcols_local = 0
    DO row_l = 1, nblkrows_local
       absrow = blkrows_local(row_l)
       nfullrows_local = nfullrows_local +&
            row_blk_end(absrow) - row_blk_begin(absrow) + 1
    ENDDO
    DO row_l = 1, nblkcols_local
       absrow = blkcols_local(row_l)
       nfullcols_local = nfullcols_local +&
            col_blk_end(absrow) - col_blk_begin(absrow) + 1
    ENDDO
    CALL get_blacs_info(blacs_env, para_env=para_env, blacs2mpi=blacs2mpi,&
         my_process_row=myprow, my_process_column=mypcol,&
         number_of_process_rows=nprows, number_of_process_columns=npcols)
    IF (dbg) CALL cp_blacs_env_write(blacs_env, 6, error)
    IF (dbg) WRITE(*,*)routineP//' myprow, mypcol',myprow,'/',nprows,';',&
         mypcol,'/',npcols
    CALL cp_assert(SIZE(prow_dist) .EQ. SIZE(pcol_dist), cp_warning_level,&
         cp_wrong_args_error, routineN,&
         'Unequal row and column distributions for symmetric matrix.', error)
    CALL cp_assert(MAXVAL(prow_dist).LE.UBOUND(blacs2mpi,1), cp_failure_level,&
         cp_wrong_args_error, routineN,&
         'Row distribution references unexistent processor rows', error)
    IF (dbg) &
         CALL cp_assert(MAXVAL(prow_dist).EQ.UBOUND(blacs2mpi,1), cp_warning_level,&
         cp_wrong_args_error, routineN,&
         'Range of row distribution not equal to processor rows', error)
    CALL cp_assert(MAXVAL(pcol_dist).LE.UBOUND(blacs2mpi,2), cp_failure_level,&
         cp_wrong_args_error, routineN,&
         'Col distribution references unexistent processor cols', error)
    IF (dbg) &
         CALL cp_assert(MAXVAL(pcol_dist).EQ.UBOUND(blacs2mpi,2), cp_warning_level,&
         cp_wrong_args_error, routineN,&
         'Range of col distribution not equal to processor cols', error)
    ALLOCATE (row_blk_size (nblkrows_total), stat=stat)
    IF (stat /= 0) CALL stop_memory (routineP, 'row_blk_size', nblkrows_total)
    ALLOCATE (col_blk_size (nblkcols_total), stat=stat)
    IF (stat /= 0) CALL stop_memory (routineP, 'col_blk_size', nblkcols_total)
    row_blk_size = row_blk_end(:)-row_blk_begin(:)+1
    col_blk_size = col_blk_end(:)-col_blk_begin(:)+1
    !
    CALL array_new (rbs_array, row_blk_size, gift=.TRUE.)
    CALL array_new (cbs_array, col_blk_size, gift=.TRUE.)
    !IF (ASSOCIATED (matrix)) CALL dbcsr_deallocate (matrix, error)
    !ALLOCATE (matrix)
    CALL cp_dbcsr_init (matrix, error)
    CALL cp_dbcsr_create(matrix, sm%name, dist, sm%symmetry,&
         rbs_array, cbs_array, nblks=0, nze=0, data_type=dbcsr_type_real_8,&
         error=error)
    ! Count
    nze = 0
    blk = 0
    DO absrow = 1, nblkrows_total
       sm_block_node => first_block_node(sm, absrow)
       DO WHILE (ASSOCIATED(sm_block_node))
          CALL get_block_node(block_node=sm_block_node, block_col=block_col,&
               block=blkdata_2d)
          ! Only add data if the block is associated
          IF (ASSOCIATED (blkdata_2d)) THEN
             IF(SIZE(blkdata_2d).GT.0) THEN
                nze = nze + SIZE(blkdata_2d)
                blk = blk + 1
             ENDIF
          ENDIF
          sm_block_node => next_block_node(sm_block_node)
       END DO
    ENDDO
    ! Allocate
    CALL cp_dbcsr_work_create(matrix,&
            nblks_guess=blk,&
            sizedata_guess=nze, error=error)
    ! Now copy the data
    DO absrow = 1, nblkrows_total
       sm_block_node => first_block_node(sm, absrow)
       DO WHILE (ASSOCIATED(sm_block_node))
          CALL get_block_node(block_node=sm_block_node, block_col=block_col,&
               block=blkdata_2d)
          ! Only add data if the block is associated
          IF (ASSOCIATED (blkdata_2d)) THEN
             IF(SIZE(blkdata_2d).GT.0) THEN
                CALL dbcsr_put_block(matrix, absrow, block_col, blkdata_2d)
             ENDIF
          ENDIF
          sm_block_node => next_block_node(sm_block_node)
       END DO
    ENDDO
    matrix%m%sparsity_id = sm%sparsity_id
    !
    CALL cp_dbcsr_finalize (matrix, error=error)
    CALL dbcsr_distribution_release (dist)
    CALL array_release (rbs_array)
    CALL array_release (cbs_array)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_from_sm


! *****************************************************************************
!> \brief Create a DBCSR matrix from a real_matrix_type
!> \par Limitations
!>      Data is never moved among processes, so it is more efficient
!>      and gracefully handles any misuse of distribution_2d.
!> \param[out] sm             cp2k sparse matrix
!> \param[in] matrix          DBCSR matrix
!> \param[dist] distribution_2d  distribution_2d for the sparse matrix
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE sm_from_dbcsr (sm, matrix, distribution_2d, error)
    TYPE(real_matrix_type), POINTER          :: sm
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    TYPE(distribution_2d_type), OPTIONAL, &
      POINTER                                :: distribution_2d
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'sm_from_dbcsr', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: data_exchange_info_cnt = 3, &
                                                metalen = 2

    CHARACTER                                :: dbcsr_matrix_type
    CHARACTER(LEN=40)                        :: sm_symmetry
    CHARACTER(LEN=default_string_length)     :: name
    INTEGER :: blk, block_col, error_handle, nblkcols_total, nblkrows_total, &
      nfullcols_total, nfullrows_total, stored_col, stored_row, timing_handle
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: sm_first_col, sm_first_row, &
                                                sm_last_col, sm_last_row
    INTEGER, DIMENSION(4)                    :: ver
    LOGICAL                                  :: sm_filled, symmetric, tr
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: blkdata_2d, buff2d
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(dbcsr_error_type)                   :: dbcsr_error
    TYPE(dbcsr_iterator)                     :: iter
    TYPE(distribution_2d_type), POINTER      :: dist2d
    TYPE(real_block_node_type), POINTER      :: sm_block_node

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
    CALL dbcsr_error_set (routineN, error_handle, error=dbcsr_error)
    CALL dbcsr_access_flush (matrix, error=dbcsr_error)
    ! Sanity checks.
    CALL cp_assert (dbcsr_valid_index (matrix), cp_warning_level,&
         cp_caller_error, routineN,&
         "DBCSR Matrix is declared invalid; expect problems.", error=error)
    CALL cp_assert (dbcsr_get_data_type (matrix) .EQ. dbcsr_type_real_8,&
         cp_failure_level, cp_unimplemented_error_nr, routineN,&
         "Non double-precision conversions not possible.", error=error)
    ! Get the dbcsr info and dimensions
    CALL dbcsr_get_info (matrix, nblkrows_total=nblkrows_total,&
         nblkcols_total=nblkcols_total, nfullrows_total=nfullrows_total,&
         nfullcols_total=nfullcols_total, name=name,&
         matrix_type=dbcsr_matrix_type)
    ! Handle dist2d
    IF (PRESENT (distribution_2d)) THEN
       dist2d => distribution_2d
    ELSE
       CALL cp_assert (.FALSE., cp_fatal_level,&
            cp_caller_error, routineN,&
            "distribution_2d is not anymore in DBCSR matrix",&
            error=error)
    ENDIF
    ! Get info and verify.
    CALL distribution_2d_get(dist2d, error=error, blacs_env=blacs_env,&
         n_row_distribution=ver(3), n_col_distribution=ver(4))
    CALL cp_assert (nblkrows_total.EQ.ver(3).AND.nblkcols_total.EQ.ver(4),&
         cp_warning_level, cp_internal_error, routineN,&
         "Specified matrix dimensions/distributions seem incompatible.")
    ! Convert type
    symmetric = .FALSE.

    SELECT CASE (dbcsr_matrix_type)
    CASE (dbcsr_type_normal)
       sm_symmetry = "none"
    CASE (dbcsr_type_symmetric)
       sm_symmetry = "symmetric"
       symmetric = .TRUE.
    CASE (dbcsr_type_antisymmetric)
       sm_symmetry = "antisymmetric"
       symmetric = .TRUE.
    CASE default
       CALL cp_unimplemented_error (routineN,&
            "Can not convert this matrix type or invalid matrix type")
    END SELECT
    ! Get row/column offsets.
    ALLOCATE (sm_first_row(nblkrows_total), sm_last_row(nblkrows_total),&
         sm_first_col(nblkcols_total), sm_last_col(nblkcols_total))
    CALL convert_sizes_to_offsets (array_data (dbcsr_row_block_sizes (matrix)),&
         sm_first_row, sm_last_row)
    CALL convert_sizes_to_offsets (array_data (dbcsr_col_block_sizes (matrix)),&
         sm_first_col, sm_last_col)
    !
    ! Only fill an existing matrix if it seems to have some
    ! structure. Currently "has structure" == "row 1 has some block".
    sm_filled = .FALSE.
    IF (ASSOCIATED(sm) .AND. nblkrows_total .GE. 1) THEN
       sm_block_node => first_block_node (sm, 1)
       sm_filled = ASSOCIATED (sm_block_node)
    ENDIF
    sm_struct_exists: IF (sm_filled) THEN
       CALL dbcsr_fill_sm_with_dbcsr(sm, matrix, error)
    ELSE
       CALL allocate_matrix(matrix=sm,&
            nblock_row=nblkrows_total,&
            nblock_col=nblkcols_total,&
            nrow=nfullrows_total,&
            ncol=nfullcols_total,&
            first_row=sm_first_row,&
            last_row=sm_last_row,&
            first_col=sm_first_col,&
            last_col=sm_last_col,&
            matrix_name=TRIM(name),&
            matrix_symmetry=sm_symmetry,&
            sparsity_id=matrix%m%sparsity_id,&
            distribution_2d=dist2d, error=error,&
            keep_pointer=.TRUE.)
       DEALLOCATE (sm_first_row, sm_last_row, sm_first_col, sm_last_col)
       !
       NULLIFY (buff2d)
       CALL dbcsr_iterator_start(iter, matrix)
       DO WHILE (dbcsr_iterator_blocks_left (iter))
          NULLIFY(blkdata_2d)
          CALL dbcsr_iterator_next_block(iter, stored_row, stored_col,&
               blkdata_2d, tr, blk)
          !
          IF(SIZE(blkdata_2d).EQ.0) CYCLE
          ! We need to get the upper triangular form for the
          ! sparse_matrix_type.
          IF (symmetric .AND. stored_row .GT. stored_col) THEN
             CALL swap (stored_row, stored_col)
             tr = .NOT. tr
          ENDIF
          IF (.NOT. tr) THEN
             CALL add_block_node(matrix=sm,&
                  block_row=stored_row, block_col=stored_col,&
                  block=blkdata_2d, error=error)
          ELSE
             ! This should never be the case when converting a
             ! sparse_matrix_type-compatible matrix back to the
             ! sparse_matrix_type, but can be the case with other
             ! matrices.
             !CALL cp_assert (.FALSE., cp_warning_level, cp_internal_error,&
             !     routineN, "Should not need to transpose")
             ALLOCATE (buff2d(SIZE(blkdata_2d,1),SIZE(blkdata_2d,2)))
             buff2d = TRANSPOSE(blkdata_2d)
             CALL add_block_node(matrix=sm,&
                  block_row=stored_row, block_col=stored_col,&
                  block=buff2d, error=error)
             DEALLOCATE (buff2d)
          ENDIF
       ENDDO
       CALL dbcsr_iterator_stop(iter)
    ENDIF sm_struct_exists
    CALL dbcsr_error_stop (error_handle, error=dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE sm_from_dbcsr


! *****************************************************************************
!> \brief Fills an existing real_matrix_type with data from a DBCSR matrix.
!>
!>        This routine is for compatibility: when a SM (real_matrix_type) is
!>        passed to a routine, converted to a DBCSR, which is transformed
!>        and we want to only update the passed SM.
!> \param[in,out] sm          cp2k sparse matrix
!> \param[in] matrix          DBCSR matrix
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_fill_sm_with_dbcsr (sm, matrix, error)
    TYPE(real_matrix_type), POINTER          :: sm
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_fill_sm_with_dbcsr', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: data_exchange_info_cnt = 3, &
                                                metalen = 2

    INTEGER                                  :: block_col, error_handle, &
                                                sm_col, sm_row, timing_handle
    LOGICAL                                  :: found
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: blkdata_2d
    TYPE(dbcsr_error_type)                   :: dbcsr_error
    TYPE(real_block_node_type), POINTER      :: sm_block_node

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
    CALL dbcsr_error_set (routineN, error_handle, error=dbcsr_error)
    CALL dbcsr_access_flush (matrix, error=dbcsr_error)
    !
    DO sm_row = 1, dbcsr_nblkrows_total(matrix)
       sm_block_node => first_block_node(sm, sm_row)
       DO WHILE (ASSOCIATED (sm_block_node))
          NULLIFY (blkdata_2d)
          CALL get_block_node (block_node=sm_block_node, block_col=sm_col,&
               block=blkdata_2d)
          ! There are several options now:
          ! * data exists in both: replace sm with dbcsr
          ! * data exists in SM but not DBCSR: clear existing block
          !
          !write(*,*)routineN//" row,col,assoc",&
          !     sm_row, sm_col, associated(blkdata_2d)
          IF (ASSOCIATED (blkdata_2d)) THEN
             CALL dbcsr_get_block (matrix, sm_row, sm_col, blkdata_2d,&
                  .FALSE., found)
             IF (.NOT. found) blkdata_2d(:,:) = 0.0_dp
          ELSEIF (found) THEN
             ! * data exists in DBCSR but not SM: add it
             !CALL add_block_node (matrix=sm, block_row=stored_row,&
             !     block_col=stored_col, block=buff2d, error=error)
             CALL cp_assert (.FALSE., cp_failure_level, cp_internal_error,&
                  routineP, "Block found in DBCSR but not SM.")
          ENDIF
          ! * no data exists in either: do nothing
          sm_block_node => next_block_node(sm_block_node)
       ENDDO
    ENDDO
    CALL dbcsr_error_stop (error_handle, error=dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_fill_sm_with_dbcsr

  SUBROUTINE test_sm(sm, error)
    TYPE(real_matrix_type), POINTER          :: sm
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'test_sm', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: data_exchange_info_cnt = 3, &
                                                metalen = 2

    INTEGER                                  :: blk, block_col, icol, irow, &
                                                proc, timing_handle
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: blkdata_2d
    TYPE(real_block_node_type), POINTER      :: block_node

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)

    proc = sm%distribution_2d%blacs_env%para_env%mepos

    DO irow=1,sm%nblock_row
       WRITE(*,*) '1proc',proc,' irow',irow
       block_node => first_block_node(sm,irow)
       DO WHILE (ASSOCIATED(block_node))
          CALL get_block_node(block_node=block_node,block_col=icol,BLOCK=blkdata_2d)

          WRITE(*,*) '1sm_test: proc',proc,' irow',irow,' icol',icol,'blk',blkdata_2d(1,1),&
               SIZE(blkdata_2d,1),SIZE(blkdata_2d,2)
          block_node => next_block_node(block_node)
       END DO
    END DO

    DO irow = 1,sm%nblock_row
       DO icol = irow,sm%nblock_row
          WRITE(*,*) '2proc',proc,' irow',irow,' icol',icol
          CALL get_block_node(sm,block_row=irow,block_col=icol,BLOCK=blkdata_2d)
          IF(ASSOCIATED(blkdata_2d)) THEN
             WRITE(*,*) '2sm_test: proc',proc,' irow',irow,' icol',icol,'blk',blkdata_2d(1,1)
          ENDIF
       ENDDO
    ENDDO

    CALL timestop(timing_handle)
  END SUBROUTINE test_sm

END MODULE cp_dbcsr_operations
