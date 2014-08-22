!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
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
  USE cp_blacs_env,                    ONLY: cp_blacs_env_type,&
                                             get_blacs_info
  USE cp_cfm_types,                    ONLY: cp_cfm_type
  USE cp_dbcsr_interface,              ONLY: &
       add_work_coordinate, array_data, array_i1d_obj, array_new, &
       array_nullify, array_release, block_set, convert_sizes_to_offsets, &
       cp_dbcsr_add, cp_dbcsr_col_block_sizes, &
       cp_dbcsr_complete_redistribute, cp_dbcsr_copy, cp_dbcsr_create, &
       cp_dbcsr_distribution, cp_dbcsr_finalize, cp_dbcsr_get_data_size, &
       cp_dbcsr_get_data_type, cp_dbcsr_get_info, cp_dbcsr_get_matrix_type, &
       cp_dbcsr_get_stored_coordinates, cp_dbcsr_init, cp_dbcsr_iterator, &
       cp_dbcsr_iterator_blocks_left, cp_dbcsr_iterator_next_block, &
       cp_dbcsr_iterator_start, cp_dbcsr_iterator_stop, cp_dbcsr_mp_new, &
       cp_dbcsr_mp_release, cp_dbcsr_multiply, cp_dbcsr_nblkcols_total, &
       cp_dbcsr_nblkrows_total, cp_dbcsr_norm, cp_dbcsr_put_block, &
       cp_dbcsr_release, cp_dbcsr_row_block_sizes, cp_dbcsr_scale, &
       cp_dbcsr_type, cp_dbcsr_valid_index, cp_dbcsr_verify_matrix, &
       cp_dbcsr_work_create, create_bl_distribution, &
       dbcsr_create_dist_block_cyclic, dbcsr_create_dist_r_unrot, &
       dbcsr_distribution_col_dist, dbcsr_distribution_init, &
       dbcsr_distribution_local_cols, dbcsr_distribution_local_rows, &
       dbcsr_distribution_mp, dbcsr_distribution_ncols, &
       dbcsr_distribution_new, dbcsr_distribution_nlocal_cols, &
       dbcsr_distribution_nlocal_rows, dbcsr_distribution_nrows, &
       dbcsr_distribution_obj, dbcsr_distribution_release, &
       dbcsr_distribution_row_dist, dbcsr_error_set, dbcsr_error_stop, &
       dbcsr_error_type, dbcsr_get_data_p, dbcsr_mp_get_coordinates, &
       dbcsr_mp_group, dbcsr_mp_hold, dbcsr_mp_mynode, dbcsr_mp_mypcol, &
       dbcsr_mp_myprow, dbcsr_mp_npcols, dbcsr_mp_nprows, dbcsr_mp_numnodes, &
       dbcsr_mp_obj, dbcsr_mp_pgrid, dbcsr_norm_frobenius, &
       dbcsr_type_antisymmetric, dbcsr_type_complex_8, &
       dbcsr_type_no_symmetry, dbcsr_type_real_4, dbcsr_type_real_8, &
       dbcsr_type_symmetric
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_release,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_2d_types,           ONLY: distribution_2d_get,&
                                             distribution_2d_type
  USE kinds,                           ONLY: dp,&
                                             real_4,&
                                             sp
  USE message_passing,                 ONLY: mp_bcast
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  !$ USE OMP_LIB
#include "common/cp_common_uses.f90"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dbcsr_operations'


  PUBLIC :: cp_dbcsr_multiply_local, cp_dbcsr_multiply_vec, cp_dbcsr_pack_vec,&
            packed_vec_scale, cp_dbcsr_unpack_vec,&
            packed_vec_bif_tech2, cp_dbcsr_copy_vec, packed_vec_build_u, packed_vec_bcast,&
            packed_vec_ini

  INTERFACE cp_dbcsr_multiply_local
     MODULE PROCEDURE cp_dbcsr_multiply_local_d,&
                      cp_dbcsr_multiply_local_s
  END INTERFACE

  ! CP2K API emulation
  PUBLIC :: cp_dbcsr_from_fm, copy_fm_to_dbcsr, copy_dbcsr_to_fm,&
            copy_dbcsr_to_cfm, copy_cfm_to_dbcsr,&
            cp_dbcsr_sm_fm_multiply, cp_dbcsr_plus_fm_fm_t,&
            copy_dbcsr_to_fm_bc, copy_fm_to_dbcsr_bc

  ! distribution_2d_type compatibility
  PUBLIC :: cp_dbcsr_dist2d_to_dist

  PUBLIC :: cp_dbcsr_copy_columns_hack

  INTERFACE cp_dbcsr_plus_fm_fm_t
     MODULE PROCEDURE cp_dbcsr_plus_fm_fm_t_native
  END INTERFACE

  PRIVATE

CONTAINS

! *****************************************************************************
!> \brief hack for dbcsr_copy_columns
!> \param matrix_b ...
!> \param matrix_a ...
!> \param ncol ...
!> \param source_start ...
!> \param target_start ...
!> \param para_env ...
!> \param blacs_env ...
!> \param error ...
!> \author vw
! *****************************************************************************
  SUBROUTINE cp_dbcsr_copy_columns_hack(matrix_b, matrix_a,&
       ncol, source_start, target_start, para_env, blacs_env, error)

    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_b
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_a
    INTEGER, INTENT(IN)                      :: ncol, source_start, &
                                                target_start
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_dbcsr_copy_columns_hack', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nfullcols_total, &
                                                nfullrows_total
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct
    TYPE(cp_fm_type), POINTER                :: fm_matrix_a, fm_matrix_b

    NULLIFY(fm_matrix_a, fm_matrix_b, fm_struct)
    CALL cp_dbcsr_get_info(matrix_a,nfullrows_total=nfullrows_total,nfullcols_total=nfullcols_total)
    CALL cp_fm_struct_create(fm_struct,context=blacs_env,nrow_global=nfullrows_total,&
         ncol_global=nfullcols_total,para_env=para_env,error=error)
    CALL cp_fm_create(fm_matrix_a,fm_struct,name="fm_matrix_a",error=error)
    CALL cp_fm_struct_release(fm_struct,error=error)

    CALL cp_dbcsr_get_info(matrix_b,nfullrows_total=nfullrows_total,nfullcols_total=nfullcols_total)
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
!> \param mp_obj ...
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
       CALL cp_dbcsr_mp_new(mp_env, pgrid, para_env%group, myproc, numproc,&
            myprow, mypcol)
    ENDIF
    CALL array_nullify (rd)
    CALL array_nullify (cd)
    CALL array_new(rd, row_dist_data)
    CALL array_new(cd, col_dist_data)
    CALL dbcsr_distribution_new(dist, mp_env, rd, cd)
    CALL cp_dbcsr_mp_release (mp_env)
    CALL array_release (rd)
    CALL array_release (cd)
  END SUBROUTINE cp_dbcsr_dist2d_to_dist

! *****************************************************************************
!> \brief multiply a dbcsr with a column vector
!>        c (packed column) = A (dbscr) * B_b (dbcsr)
!> \param[in]  matrix_a  DBCSR matrix
!> \param[in]  matrix_b  DBCSR matrix
!> \param a_row_beg ...
!> \param a_row_end ...
!> \param[in]  b_col     the column of the matrix B
!> \param[out] pkd       the result in a packed format
!>
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_multiply_vec(matrix_a, matrix_b, a_row_beg, a_row_end, b_col, pkd, error)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_a, matrix_b
    INTEGER, INTENT(in)                      :: a_row_beg, a_row_end, b_col
    REAL(dp), DIMENSION(:), POINTER          :: pkd
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_dbcsr_multiply_vec', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: col, nblkcols_total, node, &
                                                pcol, prow, row, timing_handle
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: tr
    REAL(dp), DIMENSION(:), POINTER          :: pkd_b
    TYPE(dbcsr_mp_obj)                       :: mp_obj

    CALL timeset(routineN, timing_handle)

    nblkcols_total = cp_dbcsr_nblkcols_total(matrix_a)
    col_blk_size => array_data (cp_dbcsr_col_block_sizes(matrix_a))
    row_blk_size => array_data (cp_dbcsr_row_block_sizes(matrix_a))

    mp_obj = dbcsr_distribution_mp (cp_dbcsr_distribution (matrix_b))


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
    CALL cp_dbcsr_get_stored_coordinates (matrix_b, row, col, tr, node)
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
!> \param[in] matrix_a DBCSR matrix
!> \param[in]  pkd_b       packed vector b
!> \param a_row_beg ...
!> \param a_row_end ...
!> \param[in]  b_col       the column of the packed b (this should be changed in b_col_size)
!> \param[out] pkd_c       packed vector c
!>
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_mult_pack_vec_local(matrix_a, pkd_b, a_row_beg, a_row_end, b_col, pkd_c, error)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_a
    REAL(dp), DIMENSION(:), INTENT(in)       :: pkd_b
    INTEGER, INTENT(in)                      :: a_row_beg, a_row_end, b_col
    REAL(dp), DIMENSION(:), INTENT(inout)    :: pkd_c
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_dbcsr_mult_pack_vec_local', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: b_col_size, blk, c_offset, &
                                                col, col_size, &
                                                nblkrows_total, row, &
                                                row_size, timing_handle
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size
    REAL(dp), DIMENSION(:), POINTER          :: data_p
    TYPE(cp_dbcsr_iterator)                  :: iter

    CALL timeset(routineN, timing_handle)

    nblkrows_total = cp_dbcsr_nblkrows_total(matrix_a)
    col_blk_size =>  array_data (cp_dbcsr_col_block_sizes(matrix_a))

    b_col_size = col_blk_size(b_col)

    IF (pkd_b(nblkrows_total+1).LE.0) CALL stop_program(routineN,moduleN,__LINE__,&
                                                        "1pkd_b(n+1).LE.0")

    CALL packed_vec_ini(pkd_c, nblkrows_total, error)
    !pkd_c(:) = 0.0_dp
    c_offset = nblkrows_total + 2
    !pkd_c(nblkrows_total+1) = REAL(c_offset,dp)

    CALL cp_dbcsr_iterator_start(iter, matrix_a)
    DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
       CALL cp_dbcsr_iterator_next_block(iter, row, col, data_p, blk, &
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
                  1.0_dp,data_p(1),row_size,&
                         pkd_b(INT(pkd_b(col))),col_size,&
                  1.0_dp,pkd_c(INT(pkd_c(row))),row_size)
    ENDDO
    CALL cp_dbcsr_iterator_stop(iter)

    IF (pkd_b(nblkrows_total+1).LE.0) CALL stop_program(routineN,moduleN,__LINE__,&
                                                        "2pkd_b(n+1).LE.0")
    IF (pkd_c(nblkrows_total+1).LE.0) CALL stop_program(routineN,moduleN,__LINE__,&
                                                        "2pkd_c(n+1).LE.0")

    CALL timestop(timing_handle)

  END SUBROUTINE cp_dbcsr_mult_pack_vec_local

! *****************************************************************************
!> \brief multiply a dbcsr with a replicated array
!>        c = alpha_scalar * A (dbscr) * b + c
!> \param[in] matrix_a DBSCR matrxx
!> \param[in]  vec_b        vectors b
!> \param[inout] vec_c      vectors c
!> \param[in]  ncol         nbr of columns
!> \param[in]  alpha        alpha
!>
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_multiply_local_d(matrix_a, vec_b, vec_c, ncol, alpha, error)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_a
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: vec_b
    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: vec_c
    INTEGER, INTENT(in), OPTIONAL            :: ncol
    REAL(dp), INTENT(IN), OPTIONAL           :: alpha
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_dbcsr_multiply_local_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, coloff, my_ncol, &
                                                row, rowoff, timing_handle
    LOGICAL                                  :: has_symm
    REAL(dp)                                 :: my_alpha, my_alpha2
    REAL(dp), DIMENSION(:, :), POINTER       :: data_d
    TYPE(cp_dbcsr_iterator)                  :: iter

    CALL timeset(routineN, timing_handle)


    my_alpha = 1.0_dp
    IF (PRESENT(alpha)) my_alpha = alpha

    my_ncol = SIZE(vec_b,2)
    IF(PRESENT(ncol)) my_ncol = ncol

    my_alpha2 = 0.0_dp
    IF(cp_dbcsr_get_matrix_type(matrix_a).EQ.dbcsr_type_symmetric) my_alpha2 = my_alpha
    IF(cp_dbcsr_get_matrix_type(matrix_a).EQ.dbcsr_type_antisymmetric) my_alpha2 = -my_alpha

    has_symm=(cp_dbcsr_get_matrix_type(matrix_a).EQ.dbcsr_type_symmetric.OR.&
          cp_dbcsr_get_matrix_type(matrix_a).EQ.dbcsr_type_antisymmetric)

    !$OMP PARALLEL DEFAULT(NONE) SHARED(matrix_a,vec_b,vec_c,ncol,my_alpha2,my_alpha,my_ncol,has_symm) &
    !$OMP          PRIVATE(iter,row,col,data_d,blk,rowoff,coloff)
    CALL cp_dbcsr_iterator_start(iter, matrix_a, read_only=.TRUE., dynamic = .TRUE., dynamic_byrows=.TRUE.)
    DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
       CALL cp_dbcsr_iterator_next_block(iter, row, col, data_d, blk, row_offset=rowoff, col_offset=coloff)
       IF (my_ncol.ne.1) THEN
          CALL dgemm('N','N',&
               SIZE(data_d,1),my_ncol,SIZE(data_d,2),&
               my_alpha, data_d(1,1),          SIZE(data_d,1),&
               vec_b(coloff,1), SIZE(vec_b,1), &
               1.0_dp,   vec_c(rowoff,1), SIZE(vec_c,1))
       ELSE
          CALL dgemv('N',SIZE(data_d,1),SIZE(data_d,2),&
               my_alpha, data_d(1,1),          SIZE(data_d,1),&
               vec_b(coloff,1),1,&
               1.0_dp,   vec_c(rowoff,1),1)
       ENDIF
    ENDDO
    CALL cp_dbcsr_iterator_stop(iter)
    !$OMP END PARALLEL

    ! FIXME ... in the symmetric case, the writes to vec_c depend on the column, not the row. This makes OMP-ing more difficult
    ! needs e.g. a buffer for vec_c and a reduction of that buffer.
    IF (has_symm) THEN
       CALL cp_dbcsr_iterator_start(iter, matrix_a)
       DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
          CALL cp_dbcsr_iterator_next_block(iter, row, col, data_d, blk, row_offset=rowoff, col_offset=coloff)
          IF(row.NE.col) THEN
             IF (my_ncol.ne.1) THEN
                CALL dgemm('T','N', &
                     SIZE(data_d,2), my_ncol, SIZE(data_d,1),&
                     my_alpha2, data_d(1,1),          SIZE(data_d,1), &
                     vec_b(rowoff,1), SIZE(vec_b,1),  &
                     1.0_dp,    vec_c(coloff,1), SIZE(vec_c,1))
             ELSE
                CALL dgemv('T',SIZE(data_d,1),SIZE(data_d,2),&
                     my_alpha2, data_d(1,1),          SIZE(data_d,1), &
                     vec_b(rowoff,1),1,  &
                     1.0_dp,    vec_c(coloff,1), 1)
             ENDIF
          ENDIF
       ENDDO
       CALL cp_dbcsr_iterator_stop(iter)
    ENDIF

    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_multiply_local_d

! *****************************************************************************
!> \brief ...
!> \param matrix_a ...
!> \param vec_b ...
!> \param vec_c ...
!> \param ncol ...
!> \param alpha ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_multiply_local_s(matrix_a, vec_b, vec_c, ncol, alpha, error)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_a
    REAL(KIND=real_4), DIMENSION(:, :), &
      INTENT(IN)                             :: vec_b
    REAL(KIND=real_4), DIMENSION(:, :), &
      INTENT(INOUT)                          :: vec_c
    INTEGER, INTENT(in), OPTIONAL            :: ncol
    REAL(KIND=real_4), INTENT(IN), OPTIONAL  :: alpha
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_dbcsr_multiply_local_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, coloff, my_ncol, &
                                                row, rowoff, timing_handle
    REAL(KIND=real_4)                        :: my_alpha, my_alpha2
    REAL(KIND=real_4), DIMENSION(:, :), &
      POINTER                                :: data_d
    TYPE(cp_dbcsr_iterator)                  :: iter

    CALL timeset(routineN, timing_handle)

    my_alpha = 1.0
    IF (PRESENT(alpha)) my_alpha = alpha

    my_ncol = SIZE(vec_b,2)
    IF(PRESENT(ncol)) my_ncol = ncol

    my_alpha2 = 0.0
    IF(cp_dbcsr_get_matrix_type(matrix_a).EQ.dbcsr_type_symmetric) my_alpha2 = my_alpha
    IF(cp_dbcsr_get_matrix_type(matrix_a).EQ.dbcsr_type_antisymmetric) my_alpha2 = -my_alpha

    CALL cp_dbcsr_iterator_start(iter, matrix_a)

    DO WHILE (cp_dbcsr_iterator_blocks_left(iter))

       CALL cp_dbcsr_iterator_next_block(iter, row, col, data_d, blk, row_offset=rowoff, col_offset=coloff)

       CALL sgemm('N','N',&
                  SIZE(data_d,1),my_ncol,SIZE(data_d,2),&
                  my_alpha, data_d(1,1),          SIZE(data_d,1),&
                            vec_b(coloff,1), SIZE(vec_b,1), &
                  1.0,      vec_c(rowoff,1), SIZE(vec_c,1))

       IF((cp_dbcsr_get_matrix_type(matrix_a).EQ.dbcsr_type_symmetric.OR.&
          cp_dbcsr_get_matrix_type(matrix_a).EQ.dbcsr_type_antisymmetric)) THEN
          IF(row.NE.col) THEN
             CALL sgemm('T','N', &
                        SIZE(data_d,2), my_ncol, SIZE(data_d,1),&
                        my_alpha2, data_d(1,1),          SIZE(data_d,1), &
                                   vec_b(rowoff,1), SIZE(vec_b,1),  &
                        1.0,       vec_c(coloff,1), SIZE(vec_c,1))
          ENDIF
       ENDIF
    ENDDO

    CALL cp_dbcsr_iterator_stop(iter)

    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_multiply_local_s



! *****************************************************************************
!> \brief multiply a dbcsr with a fm matrix
!> \param[in]  matrix         DBCSR matrix
!> \param fm_in full matrix
!> \param fm_out full matrix
!> \param[in]  ncol           nbr of columns
!> \param[in]  alpha          alpha
!> \param[in]  beta           beta
!>
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_sm_fm_multiply(matrix, fm_in, fm_out, ncol, alpha, beta, error)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    TYPE(cp_fm_type), POINTER                :: fm_in, fm_out
    INTEGER, INTENT(IN)                      :: ncol
    REAL(dp), INTENT(IN), OPTIONAL           :: alpha, beta
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_dbcsr_sm_fm_multiply', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: k_in, k_out, timing_handle, &
                                                timing_handle_mult
    INTEGER, DIMENSION(:), POINTER           :: in_col_blk_sizes, &
                                                out_col_blk_sizes
    TYPE(cp_dbcsr_type)                      :: in, out
    TYPE(array_i1d_obj)                      :: col_blk_size_right_in, &
                                                col_blk_size_right_out
    REAL(dp)                                 :: my_alpha, my_beta
    TYPE(dbcsr_distribution_obj)             :: dist_right_in, product_dist

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

       CALL dbcsr_create_dist_r_unrot (dist_right_in, matrix%matrix%m%dist, k_in, &
            col_blk_size_right_in)
       CALL cp_dbcsr_init(in, error)
       CALL cp_dbcsr_create(in, "D", dist_right_in, dbcsr_type_no_symmetry, &
            cp_dbcsr_row_block_sizes(matrix), col_blk_size_right_in,&
            0, 0, error=error)

       CALL cp_dbcsr_init(out, error)
       CALL dbcsr_distribution_new (product_dist,&
            dbcsr_distribution_mp (cp_dbcsr_distribution(matrix)),&
            dbcsr_distribution_row_dist (cp_dbcsr_distribution(matrix)),&
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

       CALL cp_dbcsr_create(out, "D", product_dist, dbcsr_type_no_symmetry, &
            cp_dbcsr_row_block_sizes(matrix), col_blk_size_right_out,&
            0, 0, error=error)

       CALL copy_fm_to_dbcsr(fm_in, in, error=error)
       IF(ncol.NE.k_out.OR.my_beta.NE.0.0_dp) &
            CALL copy_fm_to_dbcsr(fm_out, out, error=error)

       CALL timeset(routineN//'_core', timing_handle_mult)
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

  END SUBROUTINE cp_dbcsr_sm_fm_multiply

! *****************************************************************************
!> \brief ...
!> \param sizes1 ...
!> \param sizes2 ...
!> \param full_num ...
! *****************************************************************************
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
! *****************************************************************************
!> \brief performs the multiplication sparse_matrix*dense_mat*dens_mat^T
!>        if matrix_g is not explicitly given, matrix_v^T will be used
!>        this can be important to save the necessary redistribute for a 
!>        different matrix_g and increase performance.
!> \param sparse_matrix ...
!> \param matrix_v ...
!> \param matrix_g ...
!> \param ncol ...
!> \param alpha ...
!> \param keep_sparsity ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_plus_fm_fm_t_native(sparse_matrix,matrix_v,matrix_g,ncol,&
       alpha,keep_sparsity,error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: sparse_matrix
    TYPE(cp_fm_type), POINTER                :: matrix_v
    TYPE(cp_fm_type), OPTIONAL, POINTER      :: matrix_g
    INTEGER, INTENT(IN)                      :: ncol
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: alpha
    LOGICAL, INTENT(IN), OPTIONAL            :: keep_sparsity
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_dbcsr_plus_fm_fm_t_native', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: k, nao, timing_handle
    LOGICAL                                  :: check_product, &
                                                my_keep_sparsity
    REAL(KIND=dp)                            :: my_alpha, norm
    TYPE(array_i1d_obj)                      :: col_blk_size_left, &
                                                col_dist_left
    TYPE(cp_dbcsr_type)                      :: mat_g, mat_v, sparse_matrix2, &
                                                sparse_matrix3
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: fm_matrix
    TYPE(dbcsr_distribution_obj)             :: dist_left
    TYPE(dbcsr_mp_obj)                       :: mp

    check_product = .FALSE.

    CALL timeset(routineN, timing_handle)

    my_keep_sparsity=.TRUE.
    IF(PRESENT(keep_sparsity)) my_keep_sparsity=keep_sparsity

    IF (ncol .GT. 0) THEN
       CALL cp_assert (cp_dbcsr_valid_index (sparse_matrix), cp_fatal_level,&
            cp_caller_error, routineN, "sparse_matrix must pre-exist", error)
       !
       ! Setup matrix_v
       CALL cp_fm_get_info(matrix_v, ncol_global=k, error=error)
       !WRITE(*,*)routineN//'truncated mult k, ncol',k,ncol,' PRESENT (matrix_g)',PRESENT (matrix_g)
       mp = dbcsr_distribution_mp (cp_dbcsr_distribution(sparse_matrix))
       CALL create_bl_distribution (col_dist_left, col_blk_size_left,&
            k, dbcsr_mp_npcols (mp))
       CALL dbcsr_distribution_new (dist_left, mp,&
            dbcsr_distribution_row_dist (cp_dbcsr_distribution(sparse_matrix)),&
            col_dist_left)
       CALL array_release (col_dist_left)
       CALL cp_dbcsr_init (mat_v, error)
       CALL cp_dbcsr_create(mat_v, "DBCSR matrix_v", dist_left, dbcsr_type_no_symmetry,&
            cp_dbcsr_row_block_sizes (sparse_matrix), col_blk_size_left, 0, 0,&
            cp_dbcsr_get_data_type (sparse_matrix), error=error)
       CALL copy_fm_to_dbcsr(matrix_v, mat_v, error=error)
       CALL cp_dbcsr_verify_matrix(mat_v, error)
       !
       ! Setup matrix_g
       IF(PRESENT (matrix_g)) THEN
          CALL cp_dbcsr_init(mat_g, error)
          CALL cp_dbcsr_create(mat_g, "DBCSR matrix_g", dist_left,&
               dbcsr_type_no_symmetry,&
               cp_dbcsr_row_block_sizes (sparse_matrix),&
               cp_dbcsr_col_block_sizes (mat_v),&
               data_type=cp_dbcsr_get_data_type (sparse_matrix), error=error)
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
               retain_sparsity=my_keep_sparsity,&
               last_k = ncol,&
               error=error)
       ELSE
          CALL cp_dbcsr_multiply("N", "T", my_alpha, mat_v, mat_v,&
               1.0_dp, sparse_matrix,&
               retain_sparsity=my_keep_sparsity,&
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
          CALL copy_fm_to_dbcsr(fm_matrix,sparse_matrix2,keep_sparsity=my_keep_sparsity, error=error)
          CALL cp_dbcsr_add(sparse_matrix2,sparse_matrix,alpha_scalar=1.0_dp,&
               beta_scalar=-1.0_dp,error=error)
          CALL cp_dbcsr_norm(sparse_matrix2,which_norm=dbcsr_norm_frobenius,&
               norm_scalar=norm,error=error)
          WRITE(*,*) 'nao=',nao,' k=',k,' ncol=',ncol,' my_alpha=',my_alpha
          WRITE(*,*) 'PRESENT (matrix_g)',PRESENT (matrix_g)
          WRITE(*,*) 'matrix_type=',cp_dbcsr_get_matrix_type(sparse_matrix)
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
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_pack_vec(matrix, ivec, pkd_vec, what, error)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER, INTENT(IN)                      :: ivec
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: pkd_vec
    CHARACTER(LEN=*), INTENT(IN)             :: what
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_dbcsr_pack_vec', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: col, col_size, &
                                                nblkrows_total, nze, row, &
                                                row_size, timing_handle
    REAL(KIND=dp), DIMENSION(:), POINTER     :: data_p
    TYPE(cp_dbcsr_iterator)                  :: iter

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)

    nblkrows_total = cp_dbcsr_nblkrows_total(matrix)
    pkd_vec(1:nblkrows_total+1) = 0.0_dp ! should be big enough to hold all the data
    !
    ! let's go
    SELECT CASE(what)
    CASE('column')
       ! we need nblkrows_total+1 to store rows
       ! and 1 extra for adding new data
       pkd_vec(nblkrows_total+1) = REAL(nblkrows_total+2,dp)
       CALL cp_dbcsr_iterator_start(iter, matrix)
       DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
          CALL cp_dbcsr_iterator_next_block (iter, row, col, data_p, &
               row_size=row_size, col_size=col_size)
          IF(col.ne.ivec) CYCLE
          nze = row_size*col_size
          !write(*,*) 'row',row,' col',col,' nze',nze
          IF(nze.le.0) CYCLE
          !
          ! let's copy the block
          pkd_vec(row) = pkd_vec(nblkrows_total+1)
          pkd_vec(nblkrows_total+1) = pkd_vec(nblkrows_total+1) + REAL(nze,dp)
          IF (SIZE(pkd_vec).LT.INT(pkd_vec(nblkrows_total+1))-1) THEN
             WRITE(*,*) 'SIZE(pkd_vec)',SIZE(pkd_vec)
             WRITE(*,*) 'pkd_vec(nblkrows_total+1)',pkd_vec(nblkrows_total+1)
             CALL stop_program(routineN,moduleN,__LINE__,&
                               "col: SIZE(pkd_vec).LT.pkd_vec(nblkrows_total+1)")
          ENDIF
          CALL dcopy(nze,data_p(1),1,pkd_vec(INT(pkd_vec(row))),1)
       ENDDO
       CALL cp_dbcsr_iterator_stop(iter)
    CASE('row')
       ! we need nblkrows_total+1 to store rows
       ! and 1 extra for adding new data
       pkd_vec(nblkrows_total+1) = REAL(nblkrows_total+2,dp)
       CALL cp_dbcsr_iterator_start(iter, matrix)
       DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
          CALL cp_dbcsr_iterator_next_block (iter, row, col, data_p, &
               row_size=row_size, col_size=col_size)
          IF(row.ne.ivec) CYCLE
          nze = row_size*col_size
          !write(*,*) 'row',row,' col',col,' nze',nze
          IF(nze.le.0) CYCLE
          !
          ! let's copy the block
          pkd_vec(col) = pkd_vec(nblkrows_total+1)
          pkd_vec(nblkrows_total+1) = pkd_vec(nblkrows_total+1) + nze
          IF (SIZE(pkd_vec).LT.INT(pkd_vec(nblkrows_total+1))-1) &
             CALL stop_program(routineN,moduleN,__LINE__,&
                               "SIZE(pkd_vec).LT.pkd_vec(nblkrows_total+1)")
          CALL dcopy(nze,data_p(1),1,pkd_vec(INT(pkd_vec(col))),1)
       ENDDO
       CALL cp_dbcsr_iterator_stop(iter)
    CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,"pack what?")
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
!> \param do_sum ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_unpack_vec(matrix, ivec, pkd_vec, what, do_sum, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    INTEGER, INTENT(IN)                      :: ivec
    REAL(dp), DIMENSION(:), INTENT(IN)       :: pkd_vec
    CHARACTER(LEN=*), INTENT(IN)             :: what
    LOGICAL, INTENT(in), OPTIONAL            :: do_sum
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_dbcsr_unpack_vec', &
      routineP = moduleN//':'//routineN

    INTEGER :: col, col_s, col_size, mynode, nblkcols_total, nblkrows_total, &
      node, nze, row, row_s, row_size, timing_handle
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: my_do_sum, tr

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)

    my_do_sum = .FALSE.
    IF(PRESENT(do_sum))my_do_sum=do_sum

    row_blk_size => array_data (cp_dbcsr_row_block_sizes (matrix))
    col_blk_size => array_data (cp_dbcsr_col_block_sizes (matrix))
    nblkrows_total = cp_dbcsr_nblkrows_total(matrix)
    nblkcols_total = cp_dbcsr_nblkcols_total(matrix)
    mynode = dbcsr_mp_mynode (dbcsr_distribution_mp (cp_dbcsr_distribution (matrix)))

    !
    ! let's go
    SELECT CASE(what)
    CASE('column')
       IF (INT(pkd_vec(nblkrows_total+1)).LE.0) &
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "1pkd_vec(n+1).LE.0")
       col = ivec
       col_size = col_blk_size(col)
       DO row = 1, nblkrows_total
          IF(INT(pkd_vec(row)).LE.0) CYCLE
          row_size = row_blk_size(row)
          nze = row_size*col_size
          ! add block only on the correct node
          tr = .FALSE.
          row_s=row ; col_s=col
          CALL cp_dbcsr_get_stored_coordinates (matrix, row_s, col_s, tr, node)
          IF(node.EQ.mynode) CALL cp_dbcsr_put_block(matrix, row, col, &
               pkd_vec(INT(pkd_vec(row)):INT(pkd_vec(row)+nze-1)), &
               summation=my_do_sum)
       ENDDO ! row
       CALL cp_dbcsr_finalize(matrix, error=error)
    CASE('row')
       IF (INT(pkd_vec(nblkrows_total+1)).LE.0) &
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "2pkd_vec(n+1).LE.0")
       row = ivec
       row_size = row_blk_size(row)
       DO col = 1,nblkcols_total
          IF(INT(pkd_vec(col)).LE.0) CYCLE
          col_size = col_blk_size(col)
          nze = row_size*col_size
          ! add block only on the correct node
          tr = .FALSE.
          row_s=row ; col_s=col
          CALL cp_dbcsr_get_stored_coordinates (matrix, row_s, col_s, tr, node)
          IF(node.EQ.mynode) CALL cp_dbcsr_put_block(matrix, row, col, &
               pkd_vec(INT(pkd_vec(col)):INT(pkd_vec(col)+nze-1)),&
               summation=my_do_sum)
       ENDDO ! col
       CALL cp_dbcsr_finalize(matrix, error=error)
    CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,"unpack what?")
    END SELECT

    CALL timestop(timing_handle)

  END SUBROUTINE cp_dbcsr_unpack_vec

! *****************************************************************************
!> \brief Does one of the operations:
!>      left)  pkd_ji = alpha * block_jj (block) * pkd_ji for j = 1,..., nvec
!>      right) pkd_ij = alpha * pkd_ij * block_jj (block) for j = 1,..., nvec
!>
!> \param alpha ...
!> \param block ...
!> \param pkd ...
!> \param nvec ...
!> \param ivec ...
!> \param vec_blk_size ...
!> \param side ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE packed_vec_scale(alpha, block, pkd, nvec, ivec, vec_blk_size, side, error)
    REAL(dp), INTENT(IN)                     :: alpha
    REAL(dp), DIMENSION(:), INTENT(in)       :: block
    REAL(dp), DIMENSION(:), INTENT(inout)    :: pkd
    INTEGER, INTENT(in)                      :: nvec, ivec
    INTEGER, DIMENSION(:), INTENT(in)        :: vec_blk_size
    CHARACTER(LEN=*), INTENT(IN)             :: side
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'packed_vec_scale', &
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
       CALL stop_program(routineN,moduleN,__LINE__,"side?")
    END SELECT
    DEALLOCATE(buff)

  END SUBROUTINE packed_vec_scale

! *****************************************************************************
!> \brief ...
!> \param pkd_u ...
!> \param pkd_v ...
!> \param k ...
!> \param n ...
!> \param s ...
!> \param vec_blk_size ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE packed_vec_build_u(pkd_u, pkd_v, k, n, s, vec_blk_size, error)
    REAL(dp), DIMENSION(:), INTENT(OUT)      :: pkd_u
    REAL(dp), DIMENSION(:), INTENT(IN)       :: pkd_v
    INTEGER, INTENT(IN)                      :: k, n
    REAL(dp), INTENT(IN)                     :: s
    INTEGER, DIMENSION(:), INTENT(IN)        :: vec_blk_size
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'packed_vec_build_u', &
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
                   1.0_dp, 0.0_dp)
    pkd_u(k) = u_offset
    u_offset = u_offset + k_blk_size**2
    pkd_u(n+1) = u_offset

  END SUBROUTINE packed_vec_build_u

! *****************************************************************************
!> \brief Does the operations:
!>      v(k,:) = v(k,:) - uAk(i,k)' * inv(d(i,i)) * v(i,:) / s
!>
!> \param mat_v ...
!> \param pkd_v_fac ...
!> \param ivec ...
!> \param pkd_v ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE packed_vec_bif_tech2(mat_v, pkd_v_fac, ivec, pkd_v, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: mat_v
    REAL(dp), DIMENSION(:), INTENT(IN)       :: pkd_v_fac
    INTEGER, INTENT(IN)                      :: ivec
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: pkd_v
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'packed_vec_bif_tech2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: col, col_size, k_row, &
                                                k_row_size, nblkrows_total, &
                                                row, row_size, timing_handle
    INTEGER, DIMENSION(:), POINTER           :: row_blk_size
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: buff
    REAL(KIND=dp), DIMENSION(:), POINTER     :: data_v
    TYPE(cp_dbcsr_iterator)                  :: iter

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)

    nblkrows_total = cp_dbcsr_nblkrows_total(mat_v)
    row_blk_size => array_data (cp_dbcsr_row_block_sizes(mat_v))

    IF (pkd_v(nblkrows_total+1).LE.0) CALL stop_program(routineN,moduleN,__LINE__,&
                                                        "1pkd_v(n+1).LE.0")

    ALLOCATE(buff(SIZE(pkd_v)))

    CALL packed_vec_ini(buff, nblkrows_total, error)
    !
    !v(k,:) (pkd) = v(k,:) (pkd) - (inv(d(i,i)) * uAk(i,k))' / s (pkd) * v(i,:) (bcsr) i = 1,ivec-1
    k_row = ivec
    k_row_size = row_blk_size(k_row)

    CALL cp_dbcsr_iterator_start(iter, mat_v)
    DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
       CALL cp_dbcsr_iterator_next_block (iter, row, col, data_v, &
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
                          data_v(1), row_size,&
                  1.0_dp, pkd_v(INT(pkd_v(col))), k_row_size)
    ENDDO
    CALL cp_dbcsr_iterator_stop(iter)
    !
    IF (pkd_v(nblkrows_total+1).LE.0) CALL stop_program(routineN,moduleN,__LINE__,&
                                                        "2pkd_v(n+1).LE.0")
    DEALLOCATE(buff)

    CALL timestop(timing_handle)

  END SUBROUTINE packed_vec_bif_tech2

! *****************************************************************************
!> \brief ...
!> \param pkd_vec ...
!> \param n ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE packed_vec_ini(pkd_vec, n, error)
    REAL(dp), DIMENSION(:), INTENT(OUT)      :: pkd_vec
    INTEGER, INTENT(IN)                      :: n
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    pkd_vec(1:n) = 0.0_dp
    pkd_vec(n+1) = n + 2
  END SUBROUTINE packed_vec_ini

! *****************************************************************************
!> \brief ...
!> \param pkd_vec ...
!> \param source ...
!> \param scope ...
!> \param n ...
!> \param vec_blk_size ...
!> \param mp_obj ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE packed_vec_alltoall(pkd_vec, source, scope, n, vec_blk_size, mp_obj, error)
    REAL(dp), DIMENSION(:), POINTER          :: pkd_vec
    INTEGER, INTENT(IN)                      :: source
    CHARACTER(LEN=*), INTENT(IN)             :: scope
    INTEGER, INTENT(IN)                      :: n
    INTEGER, DIMENSION(:), INTENT(in)        :: vec_blk_size
    TYPE(dbcsr_mp_obj)                       :: mp_obj
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'packed_vec_alltoall', &
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
       CALL stop_program(routineN,moduleN,__LINE__,"NYI")
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
       CALL stop_program(routineN,moduleN,__LINE__,"how do you wanna bcast")
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

! *****************************************************************************
!> \brief ...
!> \param pkd_vec ...
!> \param source ...
!> \param scope ...
!> \param do_summation ...
!> \param n ...
!> \param vec_blk_size ...
!> \param mp_obj ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE packed_vec_bcast(pkd_vec, source, scope, do_summation, n, vec_blk_size, mp_obj, error)
    REAL(dp), DIMENSION(:), POINTER          :: pkd_vec
    INTEGER, INTENT(IN)                      :: source
    CHARACTER(LEN=*), INTENT(IN)             :: scope
    LOGICAL, INTENT(in)                      :: do_summation
    INTEGER, INTENT(IN)                      :: n
    INTEGER, DIMENSION(:), INTENT(in)        :: vec_blk_size
    TYPE(dbcsr_mp_obj)                       :: mp_obj
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'packed_vec_bcast', &
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
          buff(:) = pkd_vec
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
          buff(:) = pkd_vec
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
       CALL stop_program(routineN,moduleN,__LINE__,"how do you wanna bcast")
    END SELECT

    DEALLOCATE(buff)
    DEALLOCATE(buff2)

    CALL timestop(timing_handle)

  END SUBROUTINE packed_vec_bcast

! *****************************************************************************
!> \brief ...
!> \param pkd_vec_a ...
!> \param pkd_vec_b ...
!> \param vec_blk_size ...
!> \param n ...
!> \param do_summation ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE packed_vec_add(pkd_vec_a,pkd_vec_b,vec_blk_size,n,do_summation,error)
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: pkd_vec_a
    REAL(dp), DIMENSION(:), INTENT(IN)       :: pkd_vec_b
    INTEGER, DIMENSION(:), INTENT(IN)        :: vec_blk_size
    INTEGER, INTENT(IN)                      :: n
    LOGICAL, INTENT(in), OPTIONAL            :: do_summation
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'packed_vec_add', &
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
    IF (a_offset_last.LE.0) CALL stop_program(routineN,moduleN,__LINE__,&
                                              "a_offset_last.LT.0!!!!!!!")
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
!> \param[in]    matrix_a  matrix
!> \param[inout] matrix_b  matrix
!> \param[in]    what      row or column
!> \param[in]    ivec      the column/row index
!>
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_copy_vec(matrix_a, matrix_b, what, ivec, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_a
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_b
    CHARACTER(LEN=*), INTENT(in)             :: what
    INTEGER, INTENT(IN)                      :: ivec
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_dbcsr_copy_vec', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: col, row, timing_handle
    LOGICAL                                  :: copy_column
    REAL(KIND=dp), DIMENSION(:), POINTER     :: data_d
    TYPE(cp_dbcsr_iterator)                  :: iter

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)

    SELECT CASE(what)
    CASE('column')
       copy_column = .TRUE.
    CASE('row')
       copy_column = .FALSE.
    CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,"copy what?")
    END SELECT

    CALL cp_dbcsr_iterator_start(iter, matrix_b)
    DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
       CALL cp_dbcsr_iterator_next_block (iter, row, col, data_d)
       IF( (copy_column.AND.col.EQ.ivec) .OR. (.NOT.copy_column.AND.row.EQ.ivec) ) &
            CALL cp_dbcsr_put_block(matrix_a, row, col, data_d, summation=.FALSE.)
       !IF(.NOT.copy_column.AND.row.EQ.ivec) &
       !     CALL dbcsr_put_block(matrix_a, row, col, data_d, tr,&
       !     summation=.FALSE.)
    ENDDO
    CALL cp_dbcsr_iterator_stop(iter)
    CALL cp_dbcsr_finalize(matrix_a, error=error)

    CALL timestop(timing_handle)

  END SUBROUTINE cp_dbcsr_copy_vec

! *****************************************************************************
!> \brief Converts a cp2k full matrix into a DBCSR matrix.
!> \param[out] matrix         the created BCSR matrix
!> \param[in] fm              the cpk full matrix
!> \param[in] threshold       the threshold for determining sparsity
!> \param[in] distribution    the distribution to use for the new matrix
!> \param[in] row_blk_size sizes of row blocks
!> \param[in] col_blk_size sizes of column blocks
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE cp_dbcsr_from_fm(matrix, fm, threshold, distribution, row_blk_size,&
       col_blk_size, error)
    TYPE(cp_dbcsr_type), INTENT(OUT)         :: matrix
    TYPE(cp_fm_type), POINTER                :: fm
    REAL(KIND=dp), INTENT(IN)                :: threshold
    TYPE(dbcsr_distribution_obj)             :: distribution
    TYPE(array_i1d_obj), INTENT(IN)          :: row_blk_size, col_blk_size
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_dbcsr_from_fm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle

    CALL timeset(routineN, timing_handle)

    !CALL cp_dbcsr_init (matrix, error)! the matrix should already be initialized
    CALL cp_dbcsr_create(matrix, fm%name, distribution, dbcsr_type_no_symmetry,&
         row_blk_size, col_blk_size,&
         0, 0, dbcsr_type_real_8, error=error)
    CALL copy_fm_to_dbcsr(fm, matrix, error=error)
    CALL cp_dbcsr_verify_matrix(matrix, error)
    CALL timestop(timing_handle)

  END SUBROUTINE cp_dbcsr_from_fm


! *****************************************************************************
!> \brief   Copy a BLACS matrix to a dbcsr matrix.
!>
!>          real_matrix=beta*real_matrix+alpha*fm
!>          beta defaults to 0, alpha to 1
!> \param[in] fm              full matrix
!> \param[out] matrix         DBCSR matrix
!> \param[in] alpha           (optional) scaling of FM
!> \param[in] beta            (optional) scaling of existing SM
!> \param[in] keep_sparsity   (optional) retains the sparsity of the input
!>                            matrix
!> \param error ...
!> \date    2009-10-13
!> \par History
!>          2009-10-13 rewritten based on copy_dbcsr_to_fm
!> \author  Urban Borstnik
!> \version 2.0
! *****************************************************************************
  SUBROUTINE copy_fm_to_dbcsr(fm,matrix,alpha,beta,keep_sparsity,error)
    TYPE(cp_fm_type), POINTER                :: fm
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
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
    TYPE(cp_dbcsr_type)                      :: bc_mat
    TYPE(dbcsr_distribution_obj)             :: bc_dist
    TYPE(dbcsr_error_type)                   :: dbcsr_error

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, dbcsr_error)
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

    CALL cp_dbcsr_get_info(matrix,&
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
         dbcsr_distribution_mp (cp_dbcsr_distribution (matrix)),&
         row_blk_size, col_blk_size)        ! block-cyclic row/col sizes

    ! Create the block-cyclic DBCSR matrix
    CALL cp_dbcsr_init (bc_mat, error)
    CALL cp_dbcsr_create (bc_mat, "Block-cyclic "//matrix%matrix%m%name, bc_dist,&
         cp_dbcsr_get_matrix_type(matrix), row_blk_size, col_blk_size, 0, 0,&
         data_type=cp_dbcsr_get_data_type(matrix),error=error)

    !call dbcsr_finalize (bc_mat)
    CALL dbcsr_distribution_release (bc_dist)
    CALL array_release (row_blk_size)
    CALL array_release (col_blk_size)

    CALL cp_dbcsr_get_info(bc_mat,&
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
    local_col_sizes(:) = 0
    IF (nblkcols_local .GE. 1) THEN
       FORALL (col = 1 : nblkcols_local)
          local_col_sizes(local_cols(col)) = cbs(local_cols(col))
       END FORALL
    ENDIF

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
          CALL add_work_coordinate(bc_mat%matrix%m%wms(1), row, col, blk_p, error=dbcsr_error)
          IF (fm%use_sp) THEN
             !blk_1d_sp => bc_mat%m%wms(1)%data_area%d%r_sp(blk_p:blk_p+nze-1)
             blk_1d_sp => dbcsr_get_data_p (bc_mat%matrix%m%wms(1)%data_area,&
                  coersion=REAL(0.0, KIND=sp), lb=blk_p, ub=blk_p+nze-1)
          ELSE
             !blk_1d_dp => bc_mat%m%wms(1)%data_area%d%r_dp(blk_p:blk_p+nze-1)
             blk_1d_dp => dbcsr_get_data_p (bc_mat%matrix%m%wms(1)%data_area,&
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
    bc_mat%matrix%m%wms(1)%datasize = blk_p - 1
    CALL cp_dbcsr_finalize (bc_mat, reshuffle=.FALSE., error=error)

    ! Now convert to the desired matrix distribution
    IF (PRESENT (alpha)) THEN
       CALL stop_program(routineN,moduleN,__LINE__,'no more alpha... clean me')
    ELSE
       CALL cp_dbcsr_complete_redistribute (bc_mat, matrix,&
            keep_sparsity=keep_sparsity, error=error)
    ENDIF
    CALL cp_dbcsr_release (bc_mat, error=error)

    CALL timestop(handle)
    CALL dbcsr_error_stop(error_handler, dbcsr_error)
  END SUBROUTINE copy_fm_to_dbcsr

! *****************************************************************************
!> \brief ...
!> \param fm ...
!> \param bc_mat ...
!> \param alpha ...
!> \param beta ...
!> \param keep_sparsity ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE copy_fm_to_dbcsr_bc(fm,bc_mat,alpha,beta,keep_sparsity,error)
    TYPE(cp_fm_type), POINTER                :: fm
    TYPE(cp_dbcsr_type)                      :: bc_mat
    REAL(kind=dp), INTENT(IN), OPTIONAL      :: alpha, beta
    LOGICAL, INTENT(IN), OPTIONAL            :: keep_sparsity
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'copy_fm_to_dbcsr_bc', &
      routineP = moduleN//':'//routineN

    INTEGER :: blk_p, col, col_l, col_size, data_type, error_handler, group, &
      handle, my_blk, mypcol, myproc, myprow, nblkcols_local, nblkcols_total, &
      nblkrows_local, nblkrows_total, ncol_block, ncol_global, &
      nfullcols_local, nfullcols_total, nfullrows_local, nfullrows_total, &
      nrow_block, nrow_global, numproc, nze, row, row_l, row_size
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_col, first_row, &
                                                last_col, last_row, &
                                                local_col_sizes, &
                                                local_row_sizes
    INTEGER, DIMENSION(:), POINTER           :: blk_map, cbs, local_cols, &
                                                local_rows, rbs
    INTEGER, DIMENSION(:, :), POINTER        :: pgrid
    REAL(kind=dp)                            :: my_beta
    REAL(KIND=dp), DIMENSION(:), POINTER     :: blk_1d_dp
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: fm_block
    REAL(KIND=sp), DIMENSION(:), POINTER     :: blk_1d_sp
    REAL(kind=sp), DIMENSION(:, :), POINTER  :: fm_block_sp
    TYPE(array_i1d_obj)                      :: col_blk_size, row_blk_size
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(dbcsr_distribution_obj)             :: bc_dist
    TYPE(dbcsr_error_type)                   :: dbcsr_error
    TYPE(dbcsr_mp_obj)                       :: mp_env

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, dbcsr_error)
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
    myproc = fm%matrix_struct%para_env%mepos
    numproc = fm%matrix_struct%para_env%num_pe
    myprow = context%mepos(1)
    mypcol = context%mepos(2)
    pgrid => context%blacs2mpi

    nfullrows_total=nrow_global
    nfullcols_total=ncol_global
    CALL cp_dbcsr_mp_new(mp_env, pgrid, group, myproc, numproc,&
                      myprow, mypcol)

    ! Create a block-cyclic distribution compatible with the FM matrix.
    CALL dbcsr_distribution_init (bc_dist)
    CALL dbcsr_create_dist_block_cyclic (bc_dist,&
         nfullrows_total, nfullcols_total,& ! Actual full matrix size
         nrow_block, ncol_block,&           ! BLACS parameters
         mp_env,&
         row_blk_size, col_blk_size)        ! block-cyclic row/col sizes

    CALL cp_dbcsr_mp_release (mp_env)

    ! Create the block-cyclic DBCSR matrix
    data_type=dbcsr_type_real_8
    IF(fm%use_sp)data_type=dbcsr_type_real_4
    CALL cp_dbcsr_init (bc_mat, error)
    CALL cp_dbcsr_create (bc_mat, "Block-cyclic ", bc_dist,&
         dbcsr_type_no_symmetry, row_blk_size, col_blk_size, 0, 0,&
         data_type=data_type,error=error)

    !call dbcsr_finalize (bc_mat)
    CALL dbcsr_distribution_release (bc_dist)
    CALL array_release (row_blk_size)
    CALL array_release (col_blk_size)

    CALL cp_dbcsr_get_info(bc_mat,&
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
    local_col_sizes(:) = 0
    IF (nblkcols_local .GE. 1) THEN
       FORALL (col = 1 : nblkcols_local)
          local_col_sizes(local_cols(col)) = cbs(local_cols(col))
       END FORALL
    ENDIF

    ALLOCATE (first_row(nblkrows_total),last_row(nblkrows_total))
    ALLOCATE (first_col(nblkcols_total),last_col(nblkcols_total))
    ALLOCATE (blk_map(0:nblkrows_total*nblkcols_total))
    blk_map(0)=1
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
          my_blk=(row_l-1)*nblkcols_local+col_l-1
          nze = row_size*col_size
          CALL add_work_coordinate(bc_mat%matrix%m%wms(1), row, col, blk_map(my_blk), error=dbcsr_error)
          blk_map(my_blk+1)=blk_map(my_blk)+nze
       ENDDO bc_cols
    ENDDO bc_rows

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(my_blk,nze,row_l,col_l,row,col,blk_1d_dp,blk_1d_sp) &
!$OMP SHARED(local_rows, nblkrows_local, nblkcols_local, local_cols, blk_map, fm, bc_mat,&
!$OMP        last_row, first_row, last_col, first_col, fm_block, fm_block_sp)
    bc_rows1: DO row_l = 1, nblkrows_local
       row = local_rows (row_l)
       bc_cols1: DO col_l = 1, nblkcols_local
          col = local_cols (col_l)
          my_blk=(row_l-1)*nblkcols_local+col_l-1
          nze=blk_map(my_blk+1)-blk_map(my_blk)
          IF (fm%use_sp) THEN
             !blk_1d_sp => bc_mat%m%wms(1)%data_area%d%r_sp(blk_p:blk_p+nze-1)
             blk_1d_sp => dbcsr_get_data_p (bc_mat%matrix%m%wms(1)%data_area,&
                  coersion=REAL(0.0, KIND=sp), lb=blk_map(my_blk), ub=blk_map(my_blk)+nze-1)
          ELSE
             !blk_1d_dp => bc_mat%m%wms(1)%data_area%d%r_dp(blk_p:blk_p+nze-1)
             blk_1d_dp => dbcsr_get_data_p (bc_mat%matrix%m%wms(1)%data_area,&
                  coersion=REAL(0.0, KIND=dp), lb=blk_map(my_blk), ub=blk_map(my_blk)+nze-1)
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
       ENDDO bc_cols1
    ENDDO bc_rows1
!$OMP END PARALLEL DO
    bc_mat%matrix%m%wms(1)%datasize = blk_map(nblkrows_local*nblkcols_local) - 1
    CALL cp_dbcsr_finalize (bc_mat, reshuffle=.FALSE., error=error)
    DEALLOCATE(blk_map)

    CALL timestop(handle)
    CALL dbcsr_error_stop(error_handler, dbcsr_error)
  END SUBROUTINE copy_fm_to_dbcsr_bc
! *****************************************************************************
!> \brief   Copy a BLACS matrix to a dbcsr matrix.
!>
!>          real_matrix=beta*real_matrix+alpha*fm
!>          beta defaults to 0, alpha to 1
!> \param[in] fm              full matrix
!> \param[out] matrix         DBCSR matrix
!> \param[in] keep_sparsity   (optional) retains the sparsity of the input
!>                            matrix
!> \param error ...
!> \date    2010
!> \par History
!>          2010  copied from copy_dbcsr_to_fm
!> \author  VW
!> \version 2.0
! *****************************************************************************
  SUBROUTINE copy_cfm_to_dbcsr(fm,matrix,keep_sparsity,error)
    TYPE(cp_cfm_type), POINTER               :: fm
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
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
    TYPE(cp_dbcsr_type)                      :: bc_mat
    TYPE(dbcsr_distribution_obj)             :: bc_dist
    TYPE(dbcsr_error_type)                   :: dbcsr_error

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, dbcsr_error)
    CALL timeset(routineN,handle)

    group = fm%matrix_struct%para_env%group
    context => fm%matrix_struct%context
    nrow_block = fm%matrix_struct%nrow_block
    ncol_block = fm%matrix_struct%ncol_block
    nrow_global = fm%matrix_struct%nrow_global
    ncol_global = fm%matrix_struct%ncol_global

    CALL cp_dbcsr_get_info(matrix,&
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
         dbcsr_distribution_mp (cp_dbcsr_distribution (matrix)),&
         row_blk_size, col_blk_size)        ! block-cyclic row/col sizes

    ! Create the block-cyclic DBCSR matrix
    CALL cp_dbcsr_init (bc_mat, error)
    CALL cp_dbcsr_create (bc_mat, "Block-cyclic "//matrix%matrix%m%name, bc_dist,&
         cp_dbcsr_get_matrix_type(matrix), row_blk_size, col_blk_size, 0, 0,&
         dbcsr_type_complex_8,error=error) ! type hard coded !
    !call dbcsr_finalize (bc_mat)
    CALL dbcsr_distribution_release (bc_dist)
    CALL array_release (row_blk_size)
    CALL array_release (col_blk_size)

    CALL cp_dbcsr_get_info(bc_mat,&
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
    local_col_sizes(:) = 0
    IF (nblkcols_local .GE. 1) THEN
       FORALL (col = 1 : nblkcols_local)
          local_col_sizes(local_cols(col)) = cbs(local_cols(col))
       END FORALL
    ENDIF

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
          CALL add_work_coordinate(bc_mat%matrix%m%wms(1), row, col, blk_p, error=dbcsr_error)
          !blk_1d => bc_mat%m%wms(1)%data_area%d%c_dp(blk_p:blk_p+nze-1)
          blk_1d => dbcsr_get_data_p (bc_mat%matrix%m%wms(1)%data_area,&
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
    bc_mat%matrix%m%wms(1)%datasize = blk_p - 1
    CALL cp_dbcsr_finalize (bc_mat, reshuffle=.FALSE., error=error)

    ! Now convert to the desired matrix distribution
    CALL cp_dbcsr_complete_redistribute (bc_mat, matrix, keep_sparsity=keep_sparsity,&
         error=error)
    CALL cp_dbcsr_release (bc_mat, error=error)

    CALL timestop(handle)
    CALL dbcsr_error_stop(error_handler, dbcsr_error)
  END SUBROUTINE copy_cfm_to_dbcsr
! *****************************************************************************
!> \brief   Copy a DBCSR matrix to a BLACS matrix
!> \param[in] matrix          DBCSR matrix
!> \param[out] fm             full matrix
!> \param error ...
! *****************************************************************************
  SUBROUTINE copy_dbcsr_to_fm(matrix, fm, error)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    TYPE(cp_fm_type), POINTER                :: fm
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'copy_dbcsr_to_fm', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

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
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: blk_2d, fm_block
    REAL(KIND=sp), DIMENSION(:, :), POINTER  :: fm_block_sp
    TYPE(array_i1d_obj)                      :: cbs, rbs
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_type)                      :: bc_mat
    TYPE(dbcsr_distribution_obj)             :: bc_dist
    TYPE(dbcsr_error_type)                   :: dbcsr_error

!   ---------------------------------------------------------------------------

    CALL timeset(routineN,handle)
    CALL dbcsr_error_set (routineN, error_handle, error=dbcsr_error)

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

    CALL cp_dbcsr_get_info(matrix,&
         nfullrows_total=nfullrows_total,&
         nfullcols_total=nfullcols_total)

    ! Convert DBCSR to a block-cyclic one
    CALL dbcsr_distribution_init (bc_dist)
    CALL dbcsr_create_dist_block_cyclic (bc_dist,&
         nfullrows_total, nfullcols_total,&
         nrow_block, ncol_block,&
         dbcsr_distribution_mp (cp_dbcsr_distribution(matrix)),&
         rbs, cbs)

    CALL cp_dbcsr_init (bc_mat, error)
    CALL cp_dbcsr_create (bc_mat, "Block-cyclic"//matrix%matrix%m%name, bc_dist,&
         dbcsr_type_no_symmetry, rbs, cbs, 0, 0, error=error)
    CALL dbcsr_distribution_release (bc_dist)
    CALL array_release (rbs)
    CALL array_release (cbs)
    CALL cp_dbcsr_complete_redistribute (matrix, bc_mat, error=error)

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
    local_col_sizes(:) = 0
    IF (nblkcols_local .GE. 1) THEN
       FORALL (col = 1 : nblkcols_local)
          local_col_sizes(local_cols(col)) = col_blk_sizes(local_cols(col))
       END FORALL
    ENDIF
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
       WRITE(*,*)routineN//" dbcsr data size is", cp_dbcsr_get_data_size(bc_mat)
       WRITE(*,*)routineN//" FM block sizes are",nrow_block,'/',nfullrows_total
       WRITE(*,*)routineN//" FM block sizes are",ncol_block,'/',nfullcols_total
       WRITE(*,*)routineN//" dbcsr row sizes are",bc_mat%matrix%m%row_blk_size%low%data
       WRITE(*,*)routineN//" dbcsr col sizes are",bc_mat%matrix%m%col_blk_size%low%data
    ENDIF
    !
    CALL cp_dbcsr_iterator_start(iter, bc_mat)
    DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
       CALL cp_dbcsr_iterator_next_block (iter, row, col, blk_2d)
       ! Convert absolute coordinates to FM-local coordinates
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
    ENDDO
    CALL cp_dbcsr_iterator_stop(iter)

    CALL cp_dbcsr_release (bc_mat, error=error)
    CALL dbcsr_error_stop (error_handle, error=dbcsr_error)
    CALL timestop(handle)
  END SUBROUTINE copy_dbcsr_to_fm

! *****************************************************************************
!> \brief   Copy a DBCSR matrix to a BLACS matrix
!> \param[in] matrix          DBCSR matrix
!> \param[out] fm             full matrix
!> \param error ...
! *****************************************************************************
  SUBROUTINE copy_dbcsr_to_cfm(matrix, fm, error)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    TYPE(cp_cfm_type), POINTER               :: fm
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'copy_dbcsr_to_cfm', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

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
    TYPE(array_i1d_obj)                      :: cbs, rbs
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_type)                      :: bc_mat
    TYPE(dbcsr_distribution_obj)             :: bc_dist
    TYPE(dbcsr_error_type)                   :: dbcsr_error

!   ---------------------------------------------------------------------------

    CALL timeset(routineN,handle)
    CALL dbcsr_error_set (routineN, error_handle, error=dbcsr_error)

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

    CALL cp_dbcsr_get_info(matrix,&
         nfullrows_total=nfullrows_total,&
         nfullcols_total=nfullcols_total)

    ! Convert DBCSR to a block-cyclic one
    CALL dbcsr_distribution_init (bc_dist)
    CALL dbcsr_create_dist_block_cyclic (bc_dist,&
         nfullrows_total, nfullcols_total,&
         nrow_block, ncol_block,&
         dbcsr_distribution_mp (cp_dbcsr_distribution(matrix)),&
         rbs, cbs)

    CALL cp_dbcsr_init (bc_mat, error)
    CALL cp_dbcsr_create (bc_mat, "Block-cyclic"//matrix%matrix%m%name, bc_dist,&
         dbcsr_type_no_symmetry, rbs, cbs, 0, 0, cp_dbcsr_get_data_type(matrix),&
         error=error)
    CALL dbcsr_distribution_release (bc_dist)
    CALL array_release (rbs)
    CALL array_release (cbs)
    CALL cp_dbcsr_complete_redistribute (matrix, bc_mat, error=error)

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
    local_col_sizes(:) = 0
    IF (nblkcols_local .GE. 1) THEN
       FORALL (col = 1 : nblkcols_local)
          local_col_sizes(local_cols(col)) = col_blk_sizes(local_cols(col))
       END FORALL
    ENDIF
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
       WRITE(*,*)routineN//" dbcsr data size is", cp_dbcsr_get_data_size(bc_mat)
       WRITE(*,*)routineN//" FM block sizes are",nrow_block,'/',nfullrows_total
       WRITE(*,*)routineN//" FM block sizes are",ncol_block,'/',nfullcols_total
       WRITE(*,*)routineN//" dbcsr row sizes are",bc_mat%matrix%m%row_blk_size%low%data
       WRITE(*,*)routineN//" dbcsr col sizes are",bc_mat%matrix%m%col_blk_size%low%data
    ENDIF
    !
    CALL cp_dbcsr_iterator_start(iter, bc_mat)
    DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
       CALL cp_dbcsr_iterator_next_block (iter, row, col, blk_2d)
       ! Convert absolute coordinates to FM-local coordinates
       IF (dbg) THEN
          WRITE(*,*)routineN//" blk2d size",UBOUND(blk_2d)
          WRITE(*,*)routineN//" want to set coor.",row,col
          WRITE(*,*)routineN//" local extents",&
               first_row(row),last_row(row),first_col(col),last_col(col)
       ENDIF
       fm_block(first_row(row):last_row(row),first_col(col):last_col(col))&
            = blk_2d(:,:)
    ENDDO
    CALL cp_dbcsr_iterator_stop(iter)

    CALL cp_dbcsr_release (bc_mat, error=error)
    CALL dbcsr_error_stop (error_handle, error=dbcsr_error)
    CALL timestop(handle)
  END SUBROUTINE copy_dbcsr_to_cfm
! *****************************************************************************
!> \brief   Copy a DBCSR_BLACS matrix to a BLACS matrix 
!> \param bc_mat DBCSR matrix
!> \param[out] fm             full matrix
!> \param error ...
! *****************************************************************************
  SUBROUTINE copy_dbcsr_to_fm_bc(bc_mat, fm, error)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: bc_mat
    TYPE(cp_fm_type), POINTER                :: fm
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'copy_dbcsr_to_fm_bc', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

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
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: blk_2d, fm_block
    REAL(KIND=sp), DIMENSION(:, :), POINTER  :: blk_2d_sp, fm_block_sp
    TYPE(array_i1d_obj)                      :: cbs, rbs
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(dbcsr_distribution_obj)             :: bc_dist
    TYPE(dbcsr_error_type)                   :: dbcsr_error

!   ---------------------------------------------------------------------------

    CALL timeset(routineN,handle)
    CALL dbcsr_error_set (routineN, error_handle, error=dbcsr_error)

    ! info about the full matrix
    group = fm%matrix_struct%para_env%group
    context => fm%matrix_struct%context
    mype=context%my_pid
    npe=context%n_pid
    myprow=context%mepos(1); mypcol=context%mepos(2)
    nprow=context%num_pe(1); npcol=context%num_pe(2)
    nrow_block = fm%matrix_struct%nrow_block; ncol_block = fm%matrix_struct%ncol_block
    nrow_global = fm%matrix_struct%nrow_global; ncol_global = fm%matrix_struct%ncol_global
    nfullrows_total = nrow_global; nfullcols_total=ncol_global

    ! Convert DBCSR to a block-cyclic one
    CALL dbcsr_distribution_init (bc_dist)
    CALL dbcsr_create_dist_block_cyclic (bc_dist,&
         nfullrows_total, nfullcols_total,&
         nrow_block, ncol_block,&
         dbcsr_distribution_mp (cp_dbcsr_distribution(bc_mat)),&
         rbs, cbs)

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
    local_col_sizes(:) = 0
    IF (nblkcols_local .GE. 1) THEN
       FORALL (col = 1 : nblkcols_local)
          local_col_sizes(local_cols(col)) = col_blk_sizes(local_cols(col))
       END FORALL
    ENDIF
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
       WRITE(*,*)routineN//" FM data size is", UBOUND(fm_block_sp), cp_dbcsr_get_data_type(bc_mat)
       WRITE(*,*)routineN//" dbcsr data size is", cp_dbcsr_get_data_size(bc_mat)
       WRITE(*,*)routineN//" FM block sizes are",nrow_block,'/',nfullrows_total
       WRITE(*,*)routineN//" FM block sizes are",ncol_block,'/',nfullcols_total
       WRITE(*,*)routineN//" dbcsr row sizes are",bc_mat%matrix%m%row_blk_size%low%data
       WRITE(*,*)routineN//" dbcsr col sizes are",bc_mat%matrix%m%col_blk_size%low%data
    ENDIF
    !
    CALL cp_dbcsr_iterator_start(iter, bc_mat)
    DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
       IF(cp_dbcsr_get_data_type (bc_mat)== dbcsr_type_real_4)THEN
          CALL cp_dbcsr_iterator_next_block (iter, row, col, blk_2d_sp)
         ! Convert absolute coordinates to FM-local coordinates
          IF(fm%use_sp) THEN
            fm_block_sp(first_row(row):last_row(row),first_col(col):last_col(col))&
                 = blk_2d_sp(:,:)
          ELSE
            fm_block(first_row(row):last_row(row),first_col(col):last_col(col))&
                 =REAL(blk_2d_sp(:,:),dp)
          END IF
       ELSE
          CALL cp_dbcsr_iterator_next_block (iter, row, col, blk_2d)
          IF(fm%use_sp) THEN
            fm_block_sp(first_row(row):last_row(row),first_col(col):last_col(col))&
                 = REAL(blk_2d_sp(:,:),sp)
          ELSE
            fm_block(first_row(row):last_row(row),first_col(col):last_col(col))&
                 = blk_2d(:,:)
          END IF
       ENDIF
    ENDDO
    CALL cp_dbcsr_iterator_stop(iter)
    CALL dbcsr_distribution_release (bc_dist)
    CALL array_release (rbs)
    CALL array_release (cbs)

    CALL dbcsr_error_stop (error_handle, error=dbcsr_error)
    CALL timestop(handle)
  END SUBROUTINE copy_dbcsr_to_fm_bc

END MODULE cp_dbcsr_operations
