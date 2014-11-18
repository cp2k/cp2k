!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Testing DBCSR to CSR conversion with random matrices
!> \par History
!>       2014.10 created [Patrick Seewald]
!> \author Patrick Seewald
! *****************************************************************************
PROGRAM dbcsr_test_csr_conversions
  USE dbcsr_api,                       ONLY: &
       array_data, array_i1d_obj, array_new, array_release, csr_destroy, &
       csr_type, csr_write, dbcsr_col_block_sizes, dbcsr_create, &
       dbcsr_distribution, dbcsr_distribution_mp, dbcsr_distribution_new, &
       dbcsr_distribution_obj, dbcsr_distribution_release, dbcsr_error_type, &
       dbcsr_finalize, dbcsr_finalize_lib, dbcsr_get_stored_coordinates, &
       dbcsr_init, dbcsr_init_lib, dbcsr_max_col_size, dbcsr_max_row_size, &
       dbcsr_mp_mynode, dbcsr_mp_new, dbcsr_mp_obj, dbcsr_mp_release, &
       dbcsr_nblkcols_total, dbcsr_nblkrows_total, dbcsr_obj, dbcsr_print, &
       dbcsr_put_block, dbcsr_release, dbcsr_row_block_sizes, &
       dbcsr_set_default_config, dbcsr_type_no_symmetry, dbcsr_type_real_8
  USE dbcsr_csr_conversions,           ONLY: csr_conversion_test
  USE kinds,                           ONLY: real_8
  USE machine,                         ONLY: default_output_unit
  USE message_passing,                 ONLY: mp_cart_create,&
                                             mp_cart_rank,&
                                             mp_comm_free,&
                                             mp_environ,&
                                             mp_world_finalize,&
                                             mp_world_init

  IMPLICIT NONE

  TYPE(dbcsr_obj)              :: matrix_a
  TYPE(csr_type)               :: matrix_b
  TYPE(dbcsr_error_type)       :: error

  TYPE(array_i1d_obj)          :: col_blk_sizes, row_blk_sizes
  INTEGER                      :: nblkrows_total, nblkcols_total

  TYPE(array_i1d_obj)          :: col_dist, row_dist

  INTEGER                      :: mp_comm, group, numnodes, mynode, &
                                  prow, pcol, io_unit
  INTEGER, DIMENSION(:,:), POINTER         :: pgrid

  INTEGER, DIMENSION(:), POINTER           :: rbs, cbs, rdist, cdist
  INTEGER, DIMENSION(2)                    :: npdims, myploc, coord
  TYPE(dbcsr_mp_obj)           :: mp_env

  INTEGER                      :: max_blks_total, max_blk_size, k, seedsz
  INTEGER, ALLOCATABLE, DIMENSION(:)        ::seed
  
  REAL                         :: rn
  REAL, ALLOCATABLE, DIMENSION(:)        :: rn_array

  REAL(KIND=real_8)            :: norm, sparsity
 
  CHARACTER(LEN=10)            :: k_str,mynode_str 

  ! Set up everything as in the dbcsr example codes
  CALL mp_world_init(mp_comm)
  npdims(:) = 0
  CALL mp_cart_create (mp_comm, 2, npdims, myploc, group)
  CALL mp_environ (numnodes, mynode, group)
  ALLOCATE (pgrid(0:npdims(1)-1, 0:npdims(2)-1))
  DO prow = 0, npdims(1)-1
     DO pcol = 0, npdims(2)-1
        coord = (/ prow, pcol /)
        CALL mp_cart_rank (group, coord, pgrid(prow, pcol))
     ENDDO
  ENDDO
  CALL dbcsr_mp_new (mp_env, pgrid, group, mynode, numnodes,&
       myprow=myploc(1), mypcol=myploc(2))
  DEALLOCATE(pgrid)

  CALL dbcsr_init_lib (group=mp_comm, error=error)
  CALL dbcsr_set_default_config (error=error)

  io_unit = 0
  IF (mynode.EQ.0) io_unit = default_output_unit

  ! Set seed for random number generator
  CALL RANDOM_SEED(size = seedsz)
  ALLOCATE(seed(seedsz))
  seed = 434358235

  ! Maximum number of blocks and maximum block sizes (in 1 dimension)
  max_blks_total = 50
  max_blk_size = 10

  DO k=1,100 ! test 100 matrices

    ! set seed to synchronize random numbers on different nodes
    ! (for creating the global matrix data)

    CALL RANDOM_SEED(put = seed + k*100)
    
    CALL RANDOM_NUMBER(rn) 
    nblkrows_total = FLOOR(rn*(max_blks_total)) + 1

    CALL RANDOM_NUMBER(rn)
    nblkcols_total = FLOOR(rn*(max_blks_total)) + 1

    ALLOCATE(rn_array(MAX(nblkcols_total,nblkrows_total)))
    ALLOCATE(cbs(nblkcols_total))
    ALLOCATE(rbs(nblkrows_total))
    ALLOCATE(rdist(nblkrows_total))
    ALLOCATE(cdist(nblkcols_total))

    CALL RANDOM_NUMBER(rn_array)
    cbs = FLOOR(rn_array(1:nblkcols_total)*(max_blk_size)) + 1

    CALL RANDOM_NUMBER(rn_array)
    rbs = FLOOR(rn_array(1:nblkrows_total)*(max_blk_size)) + 1

    CALL RANDOM_NUMBER(rn)
    sparsity = rn

    CALL array_new (row_blk_sizes, rbs, gift=.TRUE.)
    CALL array_new (col_blk_sizes, cbs, gift=.TRUE.)

    CALL RANDOM_NUMBER(rn_array)
    rdist = FLOOR(rn_array(1:nblkrows_total) * npdims(1))
    CALL RANDOM_NUMBER(rn_array)
    cdist = FLOOR(rn_array(1:nblkcols_total) * npdims(2))

    CALL array_new (row_dist, rdist, gift=.TRUE.)
    CALL array_new (col_dist, cdist, gift=.TRUE.)

    CALL make_random_dbcsr_matrix(matrix_a, mp_env,nblkcols_total,nblkrows_total,&
      col_blk_sizes, row_blk_sizes, col_dist, row_dist, sparsity, error)

    CALL csr_conversion_test(matrix_a, matrix_b, norm, error=error)
    WRITE (UNIT=k_str,FMT='(I0)') k
    WRITE (UNIT=mynode_str, FMT='(I0)') mynode

    !OPEN (UNIT=1000,FILE="test_mat"//TRIM(k_str),ACTION="write",STATUS="replace")
    !CALL csr_write(matrix_b,1000,error=error)
    !CLOSE(1000)

    !IF (k .EQ. 1) THEN
    !  OPEN (UNIT=200+mynode,FILE="dbcsr_mat"//TRIM(mynode_str),ACTION="write",STATUS="replace")
    !  CALL dbcsr_print(matrix_a, matlab_format = .TRUE.,unit_nr = 200+mynode, error=error)
    !  CLOSE(200+mynode)
    !ENDIF

    CALL csr_destroy(matrix_b, error)

    IF (norm > EPSILON(norm)) THEN
      IF(io_unit > 0) WRITE(io_unit,*) "Conversion error > 0 for matrix no.", k_str
    ELSE       
      IF(io_unit > 0) WRITE(io_unit,*) "Conversion OK!"
    ENDIF

    CALL dbcsr_release (matrix_a)
    CALL array_release (row_dist)
    CALL array_release (col_dist)
    CALL array_release (row_blk_sizes)
    CALL array_release (col_blk_sizes)
    DEALLOCATE(rn_array)

  ENDDO

  DEALLOCATE(seed)

  CALL dbcsr_mp_release (mp_env)
  CALL dbcsr_finalize_lib (mp_comm, io_unit, error)
  CALL mp_comm_free(group)
  CALL mp_world_finalize()

CONTAINS

! *****************************************************************************
!> \brief Create a DBCSR matrix with random values and random blocks
!> \param matrix_a ...
!> \param mp_env ...
!> \param nblkcols_total ...
!> \param nblkrows_total ...
!> \param col_blk_sizes ...
!> \param row_blk_sizes ...
!> \param col_dist ...
!> \param row_dist ...
!> \param sparsity ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE make_random_dbcsr_matrix(matrix_a, mp_env, nblkcols_total, nblkrows_total,&
                        col_blk_sizes, row_blk_sizes, col_dist, row_dist, sparsity, error) 
    TYPE(dbcsr_obj), INTENT(OUT)             :: matrix_a
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    INTEGER, INTENT(IN)                      :: nblkcols_total, nblkrows_total
    TYPE(array_i1d_obj), INTENT(IN)          :: col_blk_sizes, row_blk_sizes, &
                                                col_dist, row_dist
    REAL(real_8), INTENT(IN)                 :: sparsity
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    INTEGER                                  :: col, col_s, max_nze, &
                                                node_holds_blk, nze, row, &
                                                row_s
    LOGICAL                                  :: tr
    REAL(real_8)                             :: rn
    REAL(real_8), ALLOCATABLE, DIMENSION(:)  :: values
    TYPE(dbcsr_distribution_obj)             :: dist

    CALL dbcsr_distribution_new (dist, mp_env, row_dist, col_dist)
  
    CALL dbcsr_init (matrix_a)
  
    CALL dbcsr_create (matrix=matrix_a,&
         name="this is my matrix a",&
         dist=dist, &
         matrix_type=dbcsr_type_no_symmetry,&
         row_blk_size=row_blk_sizes,&
         col_blk_size=col_blk_sizes,&
         data_type=dbcsr_type_real_8,&
         error=error)
  
    rbs => array_data (dbcsr_row_block_sizes (matrix_a))
    cbs => array_data (dbcsr_col_block_sizes (matrix_a))
    mynode = dbcsr_mp_mynode (dbcsr_distribution_mp (dbcsr_distribution (matrix_a)))
    max_nze = dbcsr_max_row_size (matrix_a) * dbcsr_max_col_size (matrix_a)
  
    ALLOCATE (values (max_nze))
  
    DO row = 1, dbcsr_nblkrows_total (matrix_a)
      DO col = 1, dbcsr_nblkcols_total (matrix_a)
        CALL RANDOM_NUMBER(rn)
        IF (rn .GT. sparsity) THEN
          tr = .FALSE.
          row_s=row ; col_s=col
          CALL dbcsr_get_stored_coordinates (matrix_a, row_s, col_s, tr, node_holds_blk)
          IF(node_holds_blk .EQ. mynode) THEN
            nze = rbs(row_s) * cbs(col_s)
            CALL RANDOM_NUMBER (values (1:nze))
            CALL dbcsr_put_block (matrix_a, row_s, col_s, values(1:nze))
          ENDIF
        ENDIF
      ENDDO
    ENDDO
    DEALLOCATE (values)
  
    CALL dbcsr_finalize(matrix_a, error=error)
    CALL dbcsr_distribution_release (dist)

  END SUBROUTINE make_random_dbcsr_matrix

END PROGRAM dbcsr_test_csr_conversions
