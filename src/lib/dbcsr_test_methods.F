!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Tests for CP2K DBCSR operations
!> \author  Urban Borstnik
!> \date    2010-02-08
!> \version 1.0
!>
!> <b>Modification history:</b>
!> - Created 2010-02-08
! *****************************************************************************
MODULE dbcsr_test_methods
  USE array_types,                     ONLY: array_data,&
                                             array_i1d_obj,&
                                             array_new,&
                                             array_nullify,&
                                             array_release,&
                                             array_size
  USE dbcsr_blas_operations,           ONLY: dbcsr_lapack_larnv
  USE dbcsr_block_access
  USE dbcsr_data_methods
  USE dbcsr_error_handling,            ONLY: dbcsr_assert,&
                                             dbcsr_error_type,&
                                             dbcsr_fatal_level,&
                                             dbcsr_internal_error
  USE dbcsr_kinds,                     ONLY: default_string_length,&
                                             real_8
  USE dbcsr_message_passing,           ONLY: mp_comm_null,&
                                             mp_environ
  USE dbcsr_methods
  USE dbcsr_operations
  USE dbcsr_ptr_util
  USE dbcsr_transformations
  USE dbcsr_types
  USE dbcsr_util
  USE dbcsr_work_operations

  !$ USE OMP_LIB

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: dbcsr_random_dist, dbcsr_make_random_matrix, dbcsr_make_null_mp, &
       dbcsr_make_random_block_sizes, dbcsr_make_null_dist, dbcsr_test_read_args, &
       atoi, atol, ator

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_test_methods'

CONTAINS

  SUBROUTINE dbcsr_test_read_args(narg, args)
    INTEGER, INTENT(out)                     :: narg
    CHARACTER(len=*), DIMENSION(:), &
      INTENT(out)                            :: args

    CHARACTER(len=default_string_length)     :: line
    INTEGER                                  :: istat

    narg = 0
    DO
       READ(5,*,IOSTAT=istat) line
       IF(istat.NE.0) EXIT
       IF(line(1:1).EQ.'#') CYCLE
       narg = narg + 1
       args(narg) = line
       !WRITE(*,*) 'we read <'//TRIM(args(narg))//'>'
    ENDDO

  END SUBROUTINE dbcsr_test_read_args



  SUBROUTINE dbcsr_random_dist (dist_array, dist_size, nbins)
    TYPE(array_i1d_obj), INTENT(out)         :: dist_array
    INTEGER, INTENT(in)                      :: dist_size, nbins

    INTEGER                                  :: i
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: grid_dist

!REAL, ALLOCATABLE, DIMENSION(:)          :: grid_dist
!

    ALLOCATE (grid_dist(dist_size))
    CALL array_nullify (dist_array)
    !CALL RANDOM_NUMBER (grid_dist)
    FORALL (i = 1 : dist_size)
       grid_dist(i) = MODULO (nbins-i, nbins)
    END FORALL
    !CALL array_new (dist_array, INT(grid_dist*REAL(nbins)), lb=1)
    CALL array_new (dist_array, grid_dist, lb=1)
  END SUBROUTINE dbcsr_random_dist


!> \brief Creates a random matrix.
!>        Only the master process gets the data.
  SUBROUTINE dbcsr_make_random_matrix (matrix, row_blk_sizes, col_blk_sizes,&
       name, sparsity, mp_group, data_type, symmetry, error)
    TYPE(dbcsr_obj), INTENT(out)             :: matrix
    TYPE(array_i1d_obj), INTENT(in)          :: row_blk_sizes, col_blk_sizes
    CHARACTER(len=*), INTENT(in)             :: name
    REAL(kind=real_8), INTENT(in)            :: sparsity
    INTEGER, INTENT(in)                      :: mp_group
    CHARACTER, INTENT(in), OPTIONAL          :: data_type, symmetry
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_make_random_matrix', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: my_data_type, my_symmetry
    INTEGER                                  :: col, max_nze, nze, row
    INTEGER, DIMENSION(4)                    :: iseed = (/4,3,2,1/)
    INTEGER, DIMENSION(:), POINTER           :: cbs, rbs
    REAL(kind=real_8)                        :: my_sparsity
    REAL(kind=real_8), DIMENSION(1)          :: value
    TYPE(dbcsr_data_obj)                     :: data_values
    TYPE(dbcsr_distribution_obj)             :: dist

!   ---------------------------------------------------------------------------
! Create the matrix

    CALL dbcsr_make_null_dist (dist, array_size (row_blk_sizes),&
         array_size (col_blk_sizes), group=mp_group)
    CALL dbcsr_init (matrix)
    my_data_type = dbcsr_type_real_default
    IF(PRESENT(data_type)) my_data_type = data_type
    my_symmetry = dbcsr_type_normal
    IF(PRESENT(symmetry)) my_symmetry = symmetry
    CALL dbcsr_create (matrix, name,&
         dist, my_symmetry,&
         row_blk_sizes,&
         col_blk_sizes,&
         data_type=my_data_type,&
         error=error)
    rbs => array_data (row_blk_sizes)
    cbs => array_data (col_blk_sizes)
    !
    IF (sparsity .GT. 1) THEN
       my_sparsity = sparsity / 100.0
    ELSE
       my_sparsity = sparsity
    ENDIF
    IF (dbcsr_mp_mynode (dbcsr_distribution_mp (dist)) .EQ. 0) THEN
       max_nze = dbcsr_max_row_size (matrix) * dbcsr_max_col_size (matrix)
       CALL dbcsr_work_create (matrix,&
            nblks_guess=INT(dbcsr_nblkrows_total(matrix)*dbcsr_nblkcols_total(matrix)*(1.0-sparsity)),&
            sizedata_guess=INT(dbcsr_nfullrows_total(matrix)*dbcsr_nfullcols_total(matrix)*(1.0-sparsity)),&
            work_mutable = .FALSE.,error=error)

       CALL dbcsr_data_init (data_values)
       CALL dbcsr_data_new (data_values, my_data_type, data_size=max_nze)

       DO row = 1, dbcsr_nblkrows_total (matrix)
          DO col = 1, dbcsr_nblkcols_total (matrix)
             !
             ! build the upper matrix if some symmetry
             IF(my_symmetry.NE.dbcsr_type_normal.AND.col.LT.row) CYCLE

             CALL dlarnv( 1, iseed, 1, value )
             IF (value(1) .LT. my_sparsity) CYCLE
             nze = rbs(row) * cbs(col)
             CALL dbcsr_lapack_larnv( 1, iseed, nze, data_values, error )
             CALL dbcsr_put_block (matrix, row, col, data_values)
          ENDDO
       ENDDO

       CALL dbcsr_data_release (data_values)

    ENDIF
    CALL dbcsr_distribution_release (dist)
    CALL dbcsr_finalize (matrix, error=error)
    CALL dbcsr_verify_matrix (matrix, error=error)

  END SUBROUTINE dbcsr_make_random_matrix

  SUBROUTINE dbcsr_make_random_block_sizes(block_sizes, size_sum, size_mix)
    TYPE(array_i1d_obj), INTENT(out)         :: block_sizes
    INTEGER, INTENT(in)                      :: size_sum
    INTEGER, DIMENSION(:), INTENT(in)        :: size_mix

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'dbcsr_make_random_block_sizes', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: block_size, current_sum, &
                                                istat, nblocks, nsize_mix, &
                                                selector
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: mixer
    INTEGER, DIMENSION(:), POINTER           :: sizes
    TYPE(dbcsr_error_type)                   :: error

!

    NULLIFY (sizes)
    nsize_mix = SIZE (size_mix) / 2
    ALLOCATE (mixer (3, nsize_mix))
    mixer(1, :) = size_mix(1:nsize_mix*2-1:2)
    mixer(2, :) = size_mix(2:nsize_mix*2:2)
    mixer(3, :) = 1
    nblocks = 0
    current_sum = 0
    CALL ensure_array_size (sizes, lb=1, ub=1, error=error)

    selector = 1
    !
    DO WHILE (current_sum .LT. size_sum)
       nblocks = nblocks+1
       !CALL RANDOM_NUMBER(value)
       !block_size = MIN (INT (value(1) * size_max),&
       !                  size_sum - current_sum)
       block_size = MIN (mixer(2, selector),&
            size_sum - current_sum)
       sizes(nblocks) = block_size
       current_sum = current_sum + block_size
       CALL ensure_array_size (sizes, ub=nblocks+1, factor=2.0, error=error)
       mixer(3, selector) = mixer(3, selector) + 1
       IF (mixer(3, selector) .GT. mixer(1, selector)) THEN
          mixer(3, selector) = 1
          selector = MOD (selector, nsize_mix)+1
       ENDIF
    ENDDO
    CALL array_new (block_sizes, sizes(1:nblocks), lb=1)
    current_sum = SUM (array_data (block_sizes))
    CALL dbcsr_assert (current_sum, "EQ", size_sum, dbcsr_fatal_level,&
         dbcsr_internal_error, routineN, "Incorrect block sizes",__LINE__,error)
    DEALLOCATE(mixer, sizes, STAT=istat)
    CALL dbcsr_assert (istat==0, dbcsr_fatal_level, dbcsr_internal_error, &
         routineN, "deallocation error",__LINE__,error)
        
  END SUBROUTINE dbcsr_make_random_block_sizes


  SUBROUTINE dbcsr_make_null_mp (mp_env, group)
    TYPE(dbcsr_mp_obj), INTENT(out)          :: mp_env
    INTEGER, INTENT(in), OPTIONAL            :: group

    INTEGER                                  :: mynode, numnodes

    IF (PRESENT (group)) THEN
       CALL mp_environ (numnodes, mynode, group)
       CALL dbcsr_mp_new (mp_env, RESHAPE( (/ 1 /), (/1,1/)),&
            group, mynode, numnodes,&
            myprow=0, mypcol=0)
    ELSE
       CALL dbcsr_mp_new (mp_env, RESHAPE( (/ 1 /), (/1,1/)),&
            MP_COMM_NULL, 0, 1,&
            myprow=0, mypcol=0)
    ENDIF
  END SUBROUTINE dbcsr_make_null_mp
  !
  SUBROUTINE dbcsr_make_null_dist (distribution, nblkrows, nblkcols, group)
    TYPE(dbcsr_distribution_obj), &
      INTENT(out)                            :: distribution
    INTEGER, INTENT(in)                      :: nblkrows, nblkcols
    INTEGER, INTENT(in), OPTIONAL            :: group

    INTEGER                                  :: i
    TYPE(array_i1d_obj)                      :: col_dist, row_dist
    TYPE(dbcsr_mp_obj)                       :: mp_env

    CALL dbcsr_make_null_mp (mp_env, group=group)
    CALL array_new (row_dist, (/ (0, i = 1, nblkrows) /), lb=1)
    CALL array_new (col_dist, (/ (0, i = 1, nblkcols) /), lb=1)
    CALL dbcsr_distribution_new (distribution, mp_env,&
         row_dist, col_dist)
    CALL array_release (row_dist)
    CALL array_release (col_dist)
    CALL dbcsr_mp_release (mp_env)
  END SUBROUTINE dbcsr_make_null_dist

  FUNCTION atoi(a)
    CHARACTER(len=*), INTENT(in)             :: a
    INTEGER                                  :: atoi

    READ(a,'(I6)') atoi
  END FUNCTION atoi

  FUNCTION atol(a)
    CHARACTER(len=*), INTENT(in)             :: a
    LOGICAL                                  :: atol

    READ(a,'(L1)') atol
  END FUNCTION atol

  FUNCTION ator(a)
    CHARACTER(len=*), INTENT(in)             :: a
    REAL(real_8)                             :: ator

    READ(a,'(E26.15)') ator
  END FUNCTION ator

END MODULE dbcsr_test_methods
