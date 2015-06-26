!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   DBCSR operations on distributions
!> \author  Urban Borstnik
!> \date    2010-02-18
!> \version 0.9
!>
!> <b>Modification history:</b>
!> - 2010-02-18 Moved from dbcsr_util
! *****************************************************************************
MODULE dbcsr_dist_operations

  USE array_types,                     ONLY: &
       array_data, array_equality, array_exists, array_hold, array_i1d_obj, &
       array_new, array_nullify, array_release, array_size
  USE dbcsr_config,                    ONLY: max_elements_per_block
  USE dbcsr_dist_methods,              ONLY: &
       dbcsr_distribution_col_clusters, dbcsr_distribution_col_dist, &
       dbcsr_distribution_has_threads, dbcsr_distribution_hold, &
       dbcsr_distribution_init, dbcsr_distribution_local_cols, &
       dbcsr_distribution_local_rows, dbcsr_distribution_make_threads, &
       dbcsr_distribution_mp, dbcsr_distribution_ncols, &
       dbcsr_distribution_new, dbcsr_distribution_no_threads, &
       dbcsr_distribution_nrows, dbcsr_distribution_processor, &
       dbcsr_distribution_release, dbcsr_distribution_row_clusters, &
       dbcsr_distribution_row_dist, dbcsr_distribution_thread_dist
  USE dbcsr_error_handling,            ONLY: &
       dbcsr_assert, dbcsr_caller_error, dbcsr_error_set, dbcsr_error_stop, &
       dbcsr_error_type, dbcsr_fatal_level, dbcsr_internal_error, &
       dbcsr_warning_level, dbcsr_wrong_args_error
  USE dbcsr_methods,                   ONLY: dbcsr_distribution,&
                                             dbcsr_is_initialized,&
                                             dbcsr_release_locals
  USE dbcsr_mp_methods,                ONLY: dbcsr_mp_mypcol,&
                                             dbcsr_mp_myprow,&
                                             dbcsr_mp_new_transposed,&
                                             dbcsr_mp_npcols,&
                                             dbcsr_mp_nprows,&
                                             dbcsr_mp_release
  USE dbcsr_toollib,                   ONLY: gcd,&
                                             lcm
  USE dbcsr_types,                     ONLY: &
       dbcsr_distribution_obj, dbcsr_imagedistribution_obj, dbcsr_mp_obj, &
       dbcsr_obj, dbcsr_slot_home_pcol, dbcsr_slot_home_prow, &
       dbcsr_slot_home_vpcol, dbcsr_slot_home_vprow, &
       dbcsr_slot_nblkcols_local, dbcsr_slot_nblkrows_local, dbcsr_type
  USE kinds,                           ONLY: sp
  USE min_heap,                        ONLY: heap_fill,&
                                             heap_get_first,&
                                             heap_new,&
                                             heap_release,&
                                             heap_reset_first,&
                                             heap_t

  !$ USE OMP_LIB

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_dist_operations'

  INTEGER :: idid = 0


  ! Queries
  PUBLIC :: dbcsr_get_block_index, checker_square_proc,&
            dbcsr_get_stored_coordinates, dbcsr_get_stored_block_info,&
            checker_tr, get_stored_canonical, dbcsr_find_column

  ! New/transformed distributions
  PUBLIC :: dbcsr_transpose_distribution,&
            dbcsr_transpose_dims,&
            dbcsr_reblocking_targets,&
            create_bl_distribution,&
            dbcsr_create_dist_r_unrot,&
            dbcsr_create_dist_block_cyclic

  PUBLIC :: dbcsr_create_image_dist

  ! Helper routines
  PUBLIC :: dbcsr_dist_bin,&
            make_sizes_dense,&
            image_calculator

  PUBLIC :: dbcsr_get_local_rows, dbcsr_get_local_cols,&
            dbcsr_get_local_vrows, dbcsr_get_local_vcols,&
            dbcsr_reset_locals, dbcsr_reset_vlocals

  ! For dense distributions
  PUBLIC :: dbcsr_make_dists_dense


  INTERFACE dbcsr_get_stored_block_info
     MODULE PROCEDURE get_stored_block_info_type, get_stored_block_info_obj
  END INTERFACE

  INTERFACE dbcsr_get_stored_coordinates
     MODULE PROCEDURE get_stored_coordinates_type, get_stored_coordinates_obj
  END INTERFACE

  INTERFACE dbcsr_get_block_index
     MODULE PROCEDURE get_block_index_type, get_block_index_obj
  END INTERFACE



  LOGICAL, PARAMETER :: careful_mod = .FALSE.
  LOGICAL, PARAMETER :: debug_mod = .FALSE.

CONTAINS



! *****************************************************************************
!> \brief Determines whether a transpose must be applied
!> \param[in] row   The absolute matrix row.
!> \param[in] column          The absolute matrix column.
!> \retval transpose ...
! *****************************************************************************
  ELEMENTAL FUNCTION checker_tr(row, column) RESULT(transpose)
    INTEGER, INTENT(IN)                      :: row, column
    LOGICAL                                  :: transpose

    transpose = BTEST(column+row, 0) .EQV. column.GE.row

  END FUNCTION checker_tr

! *****************************************************************************
!> \brief Determines the home process for a given logical matrix element.
!> \param[in] row             logical row
!> \param[in] col             logical column
!> \param[in] pgrid           process grid
!> \param[in] row_dist        row distribution
!> \param[in] col_dist        column distribution
!> \retval process            home process of the given element
!> \note This routine is a more low-level version of
!>       dbcsr_get_stored_coordinate without the setting the row and column
!>       to the stored position.
!> \note It assumes a square matrix.
! *****************************************************************************
  PURE FUNCTION checker_square_proc(row, col, pgrid,&
       row_dist, col_dist) RESULT(process)
    INTEGER, INTENT(IN)                      :: row, col
    INTEGER, DIMENSION(0:, 0:), INTENT(IN)   :: pgrid
    INTEGER, DIMENSION(:), INTENT(IN)        :: row_dist, col_dist
    INTEGER                                  :: process

    IF (.NOT. checker_tr (row, col)) THEN
       process = pgrid (row_dist (row), col_dist(col))
    ELSE
       process = pgrid (row_dist (col), col_dist(row))
    ENDIF
  END FUNCTION checker_square_proc


! *****************************************************************************
!> \brief Sets the correct source matrix, row, column and possible data
!>        transposition for a given matrix and matrix logical row and
!>        column.
!> \param[in] matrix          DBCSR matrix
!> \param[in,out] row         input is logical row, output is lookup row
!> \param[in,out] column      input is logical column, output is lookup
!>                            column
!> \param[in,out] transpose   input: current block state; output: whether the
!>                            data should be considered transposed
!> \param[out] processor      (optional) returns the processor on which this
!>                            block resides
! *****************************************************************************
  PURE SUBROUTINE get_stored_coordinates_type(matrix, row, column, transpose,&
       processor)
    TYPE(dbcsr_type), INTENT(IN)             :: matrix
    INTEGER, INTENT(INOUT)                   :: row, column
    LOGICAL, INTENT(INOUT)                   :: transpose
    INTEGER, INTENT(OUT), OPTIONAL           :: processor

!   ---------------------------------------------------------------------------
! SM-compatible way

    IF (PRESENT (processor)) THEN
       IF (matrix%symmetry .AND. checker_tr(row, column)) THEN
          ! The transposed way.
          processor = dbcsr_distribution_processor (matrix%dist, column, row)
       ELSE
          ! The default way.
          processor = dbcsr_distribution_processor (matrix%dist, row, column)
       ENDIF
    ENDIF
  END SUBROUTINE get_stored_coordinates_type

! *****************************************************************************
!> \brief Canonical logic
!> \param matrix ...
!> \param row ...
!> \param column ...
!> \param transpose ...
!> \param processor ...
! *****************************************************************************
  PURE SUBROUTINE get_stored_canonical(matrix, row, column,&
       transpose, processor)
    TYPE(dbcsr_type), INTENT(IN)             :: matrix
    INTEGER, INTENT(INOUT)                   :: row, column
    LOGICAL, INTENT(INOUT)                   :: transpose
    INTEGER, INTENT(OUT), OPTIONAL           :: processor

    INTEGER                                  :: tmp
    LOGICAL                                  :: straight

! The old way
!straight = matrix%transpose .OR. matrix%symmetry

    straight = matrix%symmetry
    straight = (.NOT. matrix%symmetry)&
         .OR.&
         (straight .EQV. .NOT.checker_tr(row, column))
    !transpose = .NOT. straight! .NEQV. transpose
    transpose =  straight .EQV. transpose
    IF (.NOT. straight) THEN
       tmp = row
       row = column
       column = tmp
    ENDIF
    IF (PRESENT (processor)) THEN
       processor = dbcsr_distribution_processor (matrix%dist, row, column)
    END IF
  END SUBROUTINE get_stored_canonical

! *****************************************************************************
!> \brief Sets the correct source matrix, row, column and possible data
!>        transposition for a given matrix and matrix logical row and
!>        column.
!> \param[in] matrix          DBCSR matrix
!> \param[in,out] row         input is logical row, output is lookup row
!> \param[in,out] column      input is logical column, output is lookup
!>                            column
!> \param[out] transpose      whether the data must be transposed
!> \param[out] processor      (optional) returns the processor on which this
!>                            block resides
! *****************************************************************************
  SUBROUTINE get_stored_coordinates_obj(matrix, row, column, transpose, processor)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER, INTENT(INOUT)                   :: row, column
    LOGICAL, INTENT(INOUT)                   :: transpose
    INTEGER, INTENT(OUT), OPTIONAL           :: processor

    CHARACTER(len=*), PARAMETER :: routineN = 'get_stored_coordinates_obj', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    IF (debug_mod) THEN
       CALL dbcsr_assert (dbcsr_is_initialized(matrix),&
            dbcsr_fatal_level, dbcsr_caller_error, routineN,&
            "matrix%m not associated",__LINE__,error)
    ENDIF
    IF (PRESENT (processor)) THEN
       CALL dbcsr_get_stored_coordinates (matrix%m, row, column, transpose, processor)
    ELSE
       CALL dbcsr_get_stored_coordinates (matrix%m, row, column, transpose)
    ENDIF
  END SUBROUTINE get_stored_coordinates_obj


! *****************************************************************************
!> \brief Looks up a block's index given logical coordinates.
!> \param[in] matrix          DBCSR matrix
!> \param[in] row             logical row
!> \param[in] column          logical column
!> \param[out] stored_row     row where block is actually stored
!> \param[out] stored_column  column where block is actually stored
!> \param[out] transpose      whether the data must be transposed
!> \param[out] found          whether the block was found
!> \param[out] block_number   returns the block number of the row
!>                            and column
!> \param[out] data_offset    (optional) data offset for the block; 0 if
!>                            nonexistent
! *****************************************************************************
  PURE SUBROUTINE get_block_index_type(matrix, row, column, stored_row,&
       stored_column, transpose, found, block_number, data_offset)
    TYPE(dbcsr_type), INTENT(IN)             :: matrix
    INTEGER, INTENT(IN)                      :: row, column
    INTEGER, INTENT(OUT)                     :: stored_row, stored_column
    LOGICAL, INTENT(OUT)                     :: transpose, found
    INTEGER, INTENT(OUT)                     :: block_number
    INTEGER, INTENT(OUT), OPTIONAL           :: data_offset

!   ---------------------------------------------------------------------------

    stored_row = row
    stored_column = column
    transpose = .FALSE.
    CALL dbcsr_get_stored_coordinates (matrix, stored_row, stored_column,&
         transpose)
    CALL get_stored_block_info_type(matrix, stored_row, stored_column,&
         found, block_number, data_offset=data_offset, transposed=transpose)
  END SUBROUTINE get_block_index_type

! *****************************************************************************
!> \brief Sets the correct source matrix, row, column and possible data
!>        transposition for a given matrix and matrix logical row and
!>        column.
!> \param[in] matrix          DBCSR matrix
!> \param[in] row             logical row
!> \param[in] column          logical column
!> \param[out] stored_row     row where block is actually stored
!> \param[out] stored_column  column where block is actually stored
!> \param[out] transpose      whether the data must be transposed
!> \param[out] found          whether the block was found
!> \param[out] block_number   returns the block number of the row
!>                            and column
!> \param[out] data_offset    (optional) data offset for the block; 0 if
!>                            nonexistent
! *****************************************************************************
  PURE SUBROUTINE get_block_index_obj(matrix, row, column, stored_row,&
       stored_column, transpose, found, block_number, data_offset)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER, INTENT(IN)                      :: row, column
    INTEGER, INTENT(OUT)                     :: stored_row, stored_column
    LOGICAL, INTENT(OUT)                     :: transpose, found
    INTEGER, INTENT(OUT)                     :: block_number
    INTEGER, INTENT(OUT), OPTIONAL           :: data_offset

!   ---------------------------------------------------------------------------

    IF (PRESENT (data_offset)) THEN
       CALL dbcsr_get_block_index (matrix%m, row, column, stored_row,&
            stored_column, transpose, found, block_number,&
            data_offset=data_offset)
    ELSE
       CALL dbcsr_get_block_index (matrix%m, row, column, stored_row,&
            stored_column, transpose, found, block_number)
    ENDIF
  END SUBROUTINE get_block_index_obj

! *****************************************************************************
!> \brief Returns the index to a queried block, given a real (stored) row and
!>        column
!> \param[in] matrix          bcsr matrix
!> \param[in] row             input is logical row, output is lookup row
!> \param[in] column          input is logical column, output is lookup
!>                            column
!> \param[out] found          whether the block was found
!> \param[out] block_number   returns the block number of the row
!>                            and column
!> \param[out] data_offset    (optional) data offset for the block; 0 if
!>                            nonexistent
!> \param[out] transposed     (optional) whether the block is stored transposed
!>                            according to its position
! *****************************************************************************
  PURE SUBROUTINE get_stored_block_info_type(matrix, row, column,&
       found, block_number, data_offset, transposed)
    TYPE(dbcsr_type), INTENT(IN)             :: matrix
    INTEGER, INTENT(IN)                      :: row, column
    LOGICAL, INTENT(OUT)                     :: found
    INTEGER, INTENT(OUT)                     :: block_number
    INTEGER, INTENT(OUT), OPTIONAL           :: data_offset
    LOGICAL, INTENT(OUT), OPTIONAL           :: transposed

    INTEGER                                  :: blk_last, offset

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED (matrix%row_p)) THEN
       blk_last = matrix%row_p(row+1)
       IF (blk_last .GT. 0) THEN
          CALL dbcsr_find_column(column, matrix%row_p(row)+1, blk_last,&
               matrix%col_i, matrix%blk_p,&
               block_number, found)
       ELSE
          found = .FALSE.
       ENDIF
    ELSE
       found = .FALSE.
    ENDIF
    IF (found) THEN
       IF (PRESENT (data_offset) .OR. PRESENT (transposed)) THEN
          offset = matrix%blk_p(block_number)
       ENDIF
       IF (PRESENT (data_offset)) THEN
          data_offset = ABS (offset)
       ENDIF
       IF (PRESENT (transposed)) THEN
          transposed = offset .LT. 0
       ENDIF
    ELSE
       IF (PRESENT (data_offset)) THEN
          data_offset = 0
       ENDIF
       IF (PRESENT (transposed)) THEN
          transposed = .FALSE.
       ENDIF
    ENDIF
  END SUBROUTINE get_stored_block_info_type

! *****************************************************************************
!> \brief Returns the index to a queried block, given a real (stored) row and
!>        column
!> \param[in] matrix          bcsr matrix
!> \param[in] row             input is logical row, output is lookup row
!> \param[in] column          input is logical column, output is lookup
!>                            column
!> \param[out] found          whether the block was found
!> \param[out] block_number   returns the block number of the row
!>                            and column; 0 if nonexistent
!> \param[out] data_offset    (optional) data offset for the block; 0 if
!>                            nonexistent
!> \param[out] transposed     (optional) whether the block is stored transposed
!>                            according to its position
! *****************************************************************************
  PURE SUBROUTINE get_stored_block_info_obj(matrix, row, column,&
       found, block_number, data_offset, transposed)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER, INTENT(IN)                      :: row, column
    LOGICAL, INTENT(OUT)                     :: found
    INTEGER, INTENT(OUT)                     :: block_number
    INTEGER, INTENT(OUT), OPTIONAL           :: data_offset
    LOGICAL, INTENT(OUT), OPTIONAL           :: transposed

!   ---------------------------------------------------------------------------

    IF (PRESENT (data_offset)) THEN
       IF (PRESENT (transposed)) THEN
          CALL dbcsr_get_stored_block_info (matrix%m, row, column, found,&
               block_number, data_offset, transposed)
       ELSE
          CALL dbcsr_get_stored_block_info (matrix%m, row, column, found,&
               block_number, data_offset)
       ENDIF
    ELSE
       IF (PRESENT (transposed)) THEN
          CALL dbcsr_get_stored_block_info (matrix%m, row, column, found,&
               block_number, transposed=transposed)
       ELSE
          CALL dbcsr_get_stored_block_info (matrix%m, row, column, found,&
               block_number)
       ENDIF
    ENDIF
  END SUBROUTINE get_stored_block_info_obj





! *****************************************************************************
! Distribution transformations
! *****************************************************************************



! *****************************************************************************
!> \brief Prints a distribution
!> \param[in] dist  distribution to print
! *****************************************************************************
  SUBROUTINE dbcsr_print_dist(dist)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist

    TYPE(dbcsr_mp_obj)                       :: mp

    mp = dbcsr_distribution_mp (dist)
    WRITE(*,'(1X,A,I5,"x",I5,1X,I5)')"Distribution process grid size, threads",&
         dbcsr_mp_nprows (mp), dbcsr_mp_npcols (mp), dist%d%num_threads
    IF (array_size(dist%d%row_dist_block).gt.0) &
         WRITE(*,*)"Distribution: Row dist", array_data(dist%d%row_dist_block)
    IF (array_size(dist%d%col_dist_block).gt.0) &
         WRITE(*,*)"Distribution: Col dist", array_data(dist%d%col_dist_block)
    IF (array_size(dist%d%local_rows).gt.0) &
         WRITE(*,*)"Distribution: Local rows", array_data(dist%d%local_rows)
    IF (array_size(dist%d%local_cols).gt.0) &
         WRITE(*,*)"Distribution: Local cols", array_data(dist%d%local_cols)
    IF (dist%d%has_thread_dist) &
         WRITE(*,*)"Distribution: Thread", array_data(dist%d%thread_dist)
  END SUBROUTINE dbcsr_print_dist


! *****************************************************************************
!> \brief Transposes a distribution
!> \param[out] dist_tr        transposed distribution
!> \param[in] dist_normal current distribution
! *****************************************************************************
  SUBROUTINE dbcsr_transpose_distribution(dist_tr, dist_normal)
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT)                            :: dist_tr
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist_normal

    TYPE(dbcsr_mp_obj)                       :: mp_env_tr

!   ---------------------------------------------------------------------------

    CALL dbcsr_mp_new_transposed (mp_env_tr, dbcsr_distribution_mp (&
         dist_normal))
    CALL dbcsr_distribution_new(dist_tr, mp_env_tr, &
         dbcsr_distribution_col_dist (dist_normal), &
         dbcsr_distribution_row_dist (dist_normal), &
         dbcsr_distribution_local_cols (dist_normal), &
         dbcsr_distribution_local_rows (dist_normal), &
         dbcsr_distribution_col_clusters (dist_normal), &
         dbcsr_distribution_row_clusters (dist_normal))
    CALL dbcsr_mp_release (mp_env_tr)
  END SUBROUTINE dbcsr_transpose_distribution

! *****************************************************************************
!> \brief Transposes a distribution but keeps the same mp_env
!> \param[out] dist_tr        transposed distribution
!> \param[in] dist_normal current distribution
! *****************************************************************************
  SUBROUTINE dbcsr_transpose_dims(dist_tr, dist_normal)
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT)                            :: dist_tr
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist_normal

    INTEGER, DIMENSION(:), POINTER           :: cda, rda
    TYPE(dbcsr_mp_obj)                       :: mp

!   ---------------------------------------------------------------------------

    mp = dbcsr_distribution_mp (dist_normal)
    ALLOCATE (rda (dbcsr_distribution_ncols (dist_normal)))
    CALL evenly_bin (rda, dbcsr_mp_nprows (mp))

    ALLOCATE (cda (dbcsr_distribution_nrows (dist_normal)))
    CALL evenly_bin (cda, dbcsr_mp_npcols (mp))

    CALL dbcsr_distribution_new(dist_tr,mp,rda,cda,reuse_arrays=.TRUE.)
  END SUBROUTINE dbcsr_transpose_dims


! *****************************************************************************
!> \brief Evenly distributes elements into bins
!> \param[out] array          elements to bin
!> \param[in] nbins           number of bins into which elemnts are distributed
! *****************************************************************************
  PURE SUBROUTINE evenly_bin (array, nbins)
    INTEGER, DIMENSION(:), INTENT(OUT)       :: array
    INTEGER, INTENT(IN)                      :: nbins

    INTEGER                                  :: bin, first, last, maxbinsize

!   ---------------------------------------------------------------------------

    maxbinsize = CEILING (REAL(SIZE (array)) / REAL(nbins))
    DO bin = 1, nbins
       first = 1 + (bin-1)*maxbinsize
       last = MIN (1 + (bin)*maxbinsize-1, SIZE (array))
       array (first : last) = bin - 1
    ENDDO
  END SUBROUTINE evenly_bin


! *****************************************************************************
!> \brief Makes new distribution with decimation and multiplicity
!> \param[out] new_bins      new real distribution
!> \param[out] images        new image distribution
!> \param[in] source_bins    Basis for the new distribution and images
!> \param[in] nbins          number of bins in the new real distribution
!> \param[in] multiplicity   multiplicity
!> \param[in] nimages        number of images in the new distribution
!> \par Definition of multiplicity and nimages
!>      Multiplicity and decimation (number of images) are used to
!>      match process grid coordinates on non-square process
!>      grids. Given source_nbins and target_nbins, their relation is
!>                source_nbins * target_multiplicity
!>              = target_nbins * target_nimages.
!>      It is best when both multiplicity and nimages are small. To
!>      get these two factors, then, one can use the following formulas:
!>          nimages      = lcm(source_nbins, target_nbins) / target_nbins
!>          multiplicity = target_nbins / gcd(source_nbins, target_nbins)
!>      from the target's point of view (nimages = target_nimages).
!> \par Mapping
!>      The new distribution comprises of real bins and images within
!>      bins. These can be view as target_nbins*nimages virtual
!>      columns. These same virtual columns are also
!>      source_nbins*multiplicity in number. Therefore these virtual
!>      columns are mapped from source_nbins*multiplicity onto
!>      target_bins*nimages (each target bin has nimages images):
!>      Source 4: |1 2 3|4 5 6|7 8 9|A B C| (4*3)
!>      Target 6: |1 2|3 4|5 6|7 8|9 A|B C| (6*2)
!>      multiplicity=3, nimages=2, 12 virtual columns (1-C).
!>      Source bin elements are evenly mapped into one of multiplicity
!>      virtual columns. Other (non-even, block-size aware) mappings
!>      could be better.
! *****************************************************************************
  SUBROUTINE rebin_distribution (new_bins, images, source_bins,&
       nbins, multiplicity, nimages)
    INTEGER, DIMENSION(:), INTENT(OUT)       :: new_bins, images
    INTEGER, DIMENSION(:), INTENT(IN)        :: source_bins
    INTEGER, INTENT(IN)                      :: nbins, multiplicity, nimages

    CHARACTER(len=*), PARAMETER :: routineN = 'rebin_distribution', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: bin, i, old_nbins, virtual_bin
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: bin_multiplier
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    CALL dbcsr_assert (MOD(nbins*nimages,multiplicity), "EQ", 0,&
         dbcsr_warning_level, dbcsr_internal_error,&
         routineN, "mulitplicity is not divisor of new process grid coordinate",__LINE__,error)
    old_nbins = (nbins * nimages) / multiplicity
    ALLOCATE (bin_multiplier(0:old_nbins-1))
    bin_multiplier(:) = 0
    DO i = 1, SIZE (new_bins)
       IF (i .LE. SIZE (source_bins)) THEN
          bin = source_bins (i)
       ELSE
          ! Fill remainder with a cyclic distribution
          bin = MOD (i, old_nbins)
       ENDIF
       virtual_bin = bin*multiplicity+bin_multiplier(bin)
       new_bins(i) = virtual_bin / nimages
       images(i) = 1 + MOD (virtual_bin, nimages)
       bin_multiplier(bin) = bin_multiplier(bin)+1
       IF (bin_multiplier(bin) .GE. multiplicity) THEN
          bin_multiplier(bin) = 0
       ENDIF
    ENDDO
  END SUBROUTINE rebin_distribution

! *****************************************************************************
!> \brief Makes new distribution with decimation and multiplicity
!>        Multiplicity is beeing ignored, maybe this is a bug
!> \param[out] images        new image distribution
!> \param[in] my_bins        Basis for the new images
!> \param[in] nbins          number of bins in the new real distribution
!> \param[in] multiplicity   multiplicity
!> \param[in] nimages        number of images in the new distribution
!> \par Definition of multiplicity and nimages
!>      Multiplicity and decimation (number of images) are used to
!>      match process grid coordinates on non-square process
!>      grids. Given source_nbins and target_nbins, their relation is
!>                source_nbins * target_multiplicity
!>              = target_nbins * target_nimages.
!>      It is best when both multiplicity and nimages are small. To
!>      get these two factors, then, one can use the following formulas:
!>          nimages      = lcm(source_nbins, target_nbins) / target_nbins
!>          multiplicity = target_nbins / gcd(source_nbins, target_nbins)
!>      from the target's point of view (nimages = target_nimages).
!> \par Mapping
!>      The new distribution comprises of real bins and images within
!>      bins. These can be view as target_nbins*nimages virtual
!>      columns. These same virtual columns are also
!>      source_nbins*multiplicity in number. Therefore these virtual
!>      columns are mapped from source_nbins*multiplicity onto
!>      target_bins*nimages (each target bin has nimages images):
!>      Source 4: |1 2 3|4 5 6|7 8 9|A B C| (4*3)
!>      Target 6: |1 2|3 4|5 6|7 8|9 A|B C| (6*2)
!>      multiplicity=3, nimages=2, 12 virtual columns (1-C).
!>      Source bin elements are evenly mapped into one of multiplicity
!>      virtual columns. Other (non-even, block-size aware) mappings
!>      could be better.
! *****************************************************************************
  SUBROUTINE reimage_distribution (images, my_bins,&
       nbins, multiplicity, nimages)
    INTEGER, DIMENSION(:), INTENT(OUT)       :: images
    INTEGER, DIMENSION(:), INTENT(IN)        :: my_bins
    INTEGER, INTENT(IN)                      :: nbins, multiplicity, nimages

    CHARACTER(len=*), PARAMETER :: routineN = 'reimage_distribution', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: bin, i
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: bin_multiplier

!   ---------------------------------------------------------------------------

    ALLOCATE (bin_multiplier(0:nbins-1))
    bin_multiplier(:) = 0
    DO i = 1, SIZE (my_bins)
       bin = my_bins (i)
       images(i) = 1 + bin_multiplier(bin)
       bin_multiplier(bin) = bin_multiplier(bin)+1
       IF (bin_multiplier(bin) .GE. nimages) THEN
          bin_multiplier(bin) = 0
       ENDIF
    ENDDO
  END SUBROUTINE reimage_distribution

! *****************************************************************************
!> \brief Makes new virtual distribution of rows/columns.
!> \param[in] nelements      number of elements
!> \param[out] vbins         virtual bins
!> \param[in] bins           bins to which elements belong
!> \param[in] decimation     matching between bins
!> \param[in] images         images to which element belong
! *****************************************************************************
  PURE SUBROUTINE make_vdistribution (nelements, vbins, bins, decimation, images)
    INTEGER, INTENT(IN)                      :: nelements
    INTEGER, DIMENSION(nelements), &
      INTENT(OUT)                            :: vbins
    INTEGER, DIMENSION(nelements), &
      INTENT(IN)                             :: bins
    INTEGER, INTENT(IN)                      :: decimation
    INTEGER, DIMENSION(nelements), &
      INTENT(IN)                             :: images

    CHARACTER(len=*), PARAMETER :: routineN = 'make_vdistribution', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: el

!   ---------------------------------------------------------------------------

    DO el = 1 , nelements
       vbins(el) = bins(el) * decimation + images(el)-1
    END DO
  END SUBROUTINE make_vdistribution


! *****************************************************************************
!> \brief Makes new distribution with multiplicity
!> \param[out] new_bins      new real distribution
!> \param[out] images        new image distribution
!> \param[in] source_bins    Basis for the new distribution and images
!> \param[in] source_images  Basis for the new distribution and images
!> \param[in] nbins          number of bins in the new real distribution
!> \param[in] multiplicity   multiplicity
!> \param[in] nimages        number of images in the new distribution
!> \par Definition of multiplicity and nimages
!>      Multiplicity and number of images are used to match process
!>      grid coordinates on non-square process grids. Given
!>      source_nbins and target_nbins, their relation is
!>          source_nbins * multiplicity = target_nbins * nimages.
!>      It is best when both multiplicity and nimages are small. To
!>      get these two factors, then, one can use the following formulas:
!>          nimages      = lcm(source_nbins, target_nbins) / target_nbins
!>          multiplicity = target_nbins / gcd(source_nbins, target_nbins)
!> \par Mapping
!>      The new distribution comprises of real bins and images within
!>      bins. These can be view as target_nbins*nimages virtual
!>      columns. These same virtual columns are also
!>      source_nbins*multiplicity in number. Therefore these virtual
!>      columns are mapped from source_nbins*multiplicity onto
!>      target_bins*nimages (each target bin has nimages images):
!>      Source 4: |1 2 3|4 5 6|7 8 9|A B C| (4*3)
!>      Target 6: |1 2|3 4|5 6|7 8|9 A|B C| (6*2)
!>      multiplicity=3, nimages=2, 12 virtual columns (1-C).
!>      Source bin elements are evenly mapped into one of multiplicity
!>      virtual columns. Other (non-even, block-size aware) mappings
!>      could be better.
! *****************************************************************************
  SUBROUTINE rebin_imaged_distribution (new_bins, images,&
       source_bins, source_images, nbins, multiplicity, nimages)
    INTEGER, DIMENSION(:), INTENT(OUT)       :: new_bins, images
    INTEGER, DIMENSION(:), INTENT(IN)        :: source_bins, source_images
    INTEGER, INTENT(IN)                      :: nbins, multiplicity, nimages

    CHARACTER(len=*), PARAMETER :: routineN = 'rebin_imaged_distribution', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, virtual_bin
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    DO i = 1, SIZE (new_bins)
       IF (i .LE. SIZE (source_bins)) THEN
          virtual_bin = source_bins(i)*multiplicity + source_images(i)-1
       ELSE
          ! Fill remainder with a cyclic distribution
          virtual_bin = MOD (i, nbins*nimages)
       ENDIF
       new_bins(i) = virtual_bin / nimages
       images(i) = 1 + MOD (virtual_bin, nimages)
       CALL dbcsr_assert (new_bins(i).lt.nbins, dbcsr_fatal_level, dbcsr_internal_error,&
            routineN, "Wrong bin calculation",__LINE__,error)
       CALL dbcsr_assert (images(i).le.nimages, dbcsr_fatal_level, dbcsr_internal_error,&
            routineN, "Wrong image calculation",__LINE__,error)
    ENDDO
  END SUBROUTINE rebin_imaged_distribution


! *****************************************************************************
!> \brief ...
!> \param imgdist ...
!> \param dist ...
!> \param template ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_new_image_dist (imgdist, dist,&
       template, error)
    TYPE(dbcsr_imagedistribution_obj), &
      INTENT(OUT)                            :: imgdist
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    TYPE(dbcsr_imagedistribution_obj), &
      INTENT(IN)                             :: template
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_new_image_dist', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    idid = idid+1
    ALLOCATE (imgdist%i)
    imgdist%i%refcount = 1
    imgdist%i%id = idid
    imgdist%i%row_decimation   = template%i%row_decimation
    imgdist%i%row_multiplicity = template%i%row_multiplicity
    imgdist%i%col_decimation   = template%i%col_decimation
    imgdist%i%col_multiplicity = template%i%col_multiplicity
    !
    NULLIFY (imgdist%i%other_vl_rows)
    NULLIFY (imgdist%i%other_vl_cols)
    CALL array_nullify (imgdist%i%global_vrow_map)
    CALL array_nullify (imgdist%i%global_vcol_map)
    imgdist%i%has_other_vl_rows = .FALSE.
    imgdist%i%has_other_vl_cols = .FALSE.
    imgdist%i%has_global_vrow_map = .FALSE.
    imgdist%i%has_global_vcol_map = .FALSE.
    !
    imgdist%i%main = dist
    CALL dbcsr_distribution_hold (imgdist%i%main)
    !
  END SUBROUTINE dbcsr_new_image_dist

! *****************************************************************************
!> \brief Creates an image distribution given the other compatibility images
!> \param[out] imgdist                  distribution repetition
!> \param[in] dist                      distribution for which to form the
!>                                      image distribution
!> \param[in] match_row_pdist           (optional) match the new row
!>                                      distribution to this row distribution
!> \param[in] match_row_idist           (optional) match the row
!>                                      distribution to these row images
!> \param[in] match_row_nbins           number of bins in the distribution to
!>                                      match the local rows
!> \param[in] match_col_pdist           (optional) match the new column
!>                                      distribution to this column
!>                                      distribution
!> \param[in] match_col_idist           (optional) match the new column
!>                                      distribution to these column images
!> \param[in] match_col_nbins           number of bins in the distribution to
!>                                      match the local columns
!> \param nclusters_rows ...
!> \param nclusters_cols ...
! *****************************************************************************
  SUBROUTINE dbcsr_create_image_dist(imgdist, dist,&
       match_row_pdist, match_row_idist, match_row_nbins,&
       match_col_pdist, match_col_idist, match_col_nbins,&
       nclusters_rows, nclusters_cols)
    TYPE(dbcsr_imagedistribution_obj), &
      INTENT(OUT)                            :: imgdist
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: match_row_pdist, &
                                                match_row_idist
    INTEGER, INTENT(IN)                      :: match_row_nbins
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: match_col_pdist, &
                                                match_col_idist
    INTEGER, INTENT(IN)                      :: match_col_nbins, &
                                                nclusters_rows, nclusters_cols

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_create_image_dist', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ncols, npcols, nprows, nrows
    INTEGER, DIMENSION(:), POINTER           :: col_dist_data, col_img_data, &
                                                col_vdist, row_dist_data, &
                                                row_img_data, row_vdist
    LOGICAL                                  :: new_col_dist, new_row_dist
    TYPE(dbcsr_distribution_obj)             :: new_dist
    TYPE(dbcsr_mp_obj)                       :: mp_env

!   ---------------------------------------------------------------------------

    idid = idid+1
    ALLOCATE (imgdist%i)
    imgdist%i%refcount = 1
    imgdist%i%id = idid
    CALL dbcsr_distribution_init (new_dist)
    CALL dbcsr_distribution_init (imgdist%i%main)
    mp_env = dbcsr_distribution_mp (dist)
    ! Determine the factors.
    nrows = dbcsr_distribution_nrows (dist)
    ncols = dbcsr_distribution_ncols (dist)
    nprows = dbcsr_mp_nprows (mp_env)
    npcols = dbcsr_mp_npcols (mp_env)
    IF (debug_mod) WRITE (*,'(1X,A,I5,"x",I5)')routineN//"pgrid",&
         nprows, npcols
    !
    ! Create the new row distribution and row image distribution
    imgdist%i%row_decimation = nclusters_rows / nprows
    imgdist%i%row_multiplicity = nclusters_rows / gcd (nclusters_rows, match_row_nbins)
    new_row_dist = .FALSE.
    !
    IF (debug_mod) WRITE (*,*)routineN//'row decimation, multiplicity',&
         imgdist%i%row_decimation, imgdist%i%row_multiplicity
    IF (debug_mod) WRITE(*,*)routineN//" nprows, match prows", nprows, match_row_nbins
    ALLOCATE (row_img_data (nrows))
    IF (imgdist%i%row_decimation.EQ.1 .AND. imgdist%i%row_multiplicity.EQ.1 .AND. &
         .NOT. PRESENT (match_row_pdist)) THEN
       row_dist_data => dbcsr_distribution_row_dist (dist)
       row_img_data(:) = 1
    ELSE
       IF (PRESENT (match_row_pdist)) THEN
          ALLOCATE (row_dist_data (nrows))
          new_row_dist = .TRUE.
          IF (PRESENT (match_row_idist)) THEN
             CALL rebin_imaged_distribution (row_dist_data, row_img_data,&
                  match_row_pdist, match_row_idist,&
                  nprows,&
                  imgdist%i%row_multiplicity, imgdist%i%row_decimation)
          ELSE
             CALL rebin_distribution (row_dist_data, row_img_data,&
                  match_row_pdist, &
                  nprows,&
                  imgdist%i%row_multiplicity, imgdist%i%row_decimation)
          ENDIF
       ELSE
          row_dist_data => dbcsr_distribution_row_dist (dist)
          CALL reimage_distribution (row_img_data,&
               row_dist_data,&
               nprows,&
               imgdist%i%row_multiplicity, imgdist%i%row_decimation) !the multiplicity argument is beeing ignored
       ENDIF
    ENDIF
    !
    ALLOCATE (row_vdist (nrows))
    CALL make_vdistribution (nrows, row_vdist, row_dist_data,&
         imgdist%i%row_decimation, row_img_data)
    CALL array_new (imgdist%i%vrow_dist, row_vdist, gift=.TRUE.)
    !
    ! Create the new column distribution and column image distribution
    imgdist%i%col_decimation = nclusters_cols / npcols
    imgdist%i%col_multiplicity = nclusters_cols / gcd (nclusters_cols, match_col_nbins)
    new_col_dist = .FALSE.
    !
    IF (debug_mod) WRITE (*,*)routineN//'col decimation, multiplicity',&
         imgdist%i%col_decimation, imgdist%i%col_multiplicity
    IF (debug_mod) WRITE(*,*)routineN//" npcols, match pcols", npcols, match_col_nbins
    ALLOCATE (col_img_data (ncols))
    IF (imgdist%i%col_decimation.EQ.1 .AND. imgdist%i%col_multiplicity.EQ.1 .AND. &
         .NOT. PRESENT (match_col_pdist)) THEN
       col_dist_data => dbcsr_distribution_col_dist (dist)
       col_img_data(:) = 1
    ELSE
       IF (PRESENT (match_col_pdist)) THEN
          ALLOCATE (col_dist_data (ncols))
          new_col_dist = .TRUE.
          IF (PRESENT (match_col_idist)) THEN
             CALL rebin_imaged_distribution (col_dist_data, col_img_data,&
                  match_col_pdist, match_col_idist,&
                  npcols,&
                  imgdist%i%col_multiplicity, imgdist%i%col_decimation)
          ELSE
             CALL rebin_distribution (col_dist_data, col_img_data,&
                  match_col_pdist,&
                  npcols,&
                  imgdist%i%col_multiplicity, imgdist%i%col_decimation)
          ENDIF
       ELSE
          col_dist_data => dbcsr_distribution_col_dist (dist)
          CALL reimage_distribution (col_img_data,&
               col_dist_data,&
               npcols,&
               imgdist%i%col_multiplicity, imgdist%i%col_decimation) !the multiplicity argument is beeing ignored
       ENDIF
    ENDIF
    !
    ALLOCATE (col_vdist (ncols))
    CALL make_vdistribution (ncols, col_vdist, col_dist_data,&
         imgdist%i%col_decimation, col_img_data)
    CALL array_new (imgdist%i%vcol_dist, col_vdist, gift=.TRUE.)
    !
    ! Copy the row & column distribution from old distribution
    IF (new_row_dist.AND.new_col_dist) THEN
       CALL dbcsr_distribution_new(new_dist,&
            mp_env,&
            row_dist_data, col_dist_data, reuse_arrays=.TRUE.)
    ELSE
       CALL dbcsr_distribution_new(new_dist,&
            mp_env,&
            row_dist_data, col_dist_data)
       IF (new_row_dist) DEALLOCATE(row_dist_data)
       IF (new_col_dist) DEALLOCATE(col_dist_data)
    ENDIF
    ! Now setup the list of distributions, ... .
    ! Now finish the distribution image.
    imgdist%i%main = new_dist
    CALL array_new (imgdist%i%col_image, col_img_data, gift=.TRUE.)
    CALL array_new (imgdist%i%row_image, row_img_data, gift=.TRUE.)
    !
    imgdist%i%has_other_vl_rows = .FALSE.
    imgdist%i%has_other_vl_cols = .FALSE.
    imgdist%i%has_global_vrow_map = .FALSE.
    imgdist%i%has_global_vcol_map = .FALSE.
    !
  END SUBROUTINE dbcsr_create_image_dist



! *****************************************************************************
!> \brief Calculates the intersections of blocks
!> \param[out] ints           intersections of blocks
!> \param[in] numints        maximum number of expected intersections 
!> \param[out] n_src_dsts     offset and number intersections belonging
!>                            to source blocks
!> \param[in] src_sizes       sizes of source blocks
!> \param[in] dst_sizes       sizes of target blocks
!> \param[in,out] error       error
!> \par ints output format
!>      The ints array should be up to twice as large as the number of
!>      intersecting blocks. Each entry is comprised of the target
!>      block and the common length along with the offsets of the
!>      intersection in the old and new blocks.
!> \par n_src_dsts format
!>      This arrays stored the number of intersecting blocks in common
!>      (position 2) and the offset of the first common intersecting
!>      block (position 1).
!> \par Interpretation (Mapping old blocks into new blocks)
!>      The old element belongs to block B. Lookup row B in
!>      n_src_dsts.  The first element F tells you which is the first
!>      new block to map into and the second element tells you into
!>      how many new blocks N you have to map.  You then look up rows
!>      F to F+N-1 in ints. The first block tells you into which block
!>      it is mapped and the secend element tells you how many
!>      elements they have in common. The third element specifies the
!>      offset of the intesection in the old block while the fourth
!>      specifies the offset of the intersection in the new block.
!> \par History
!>      UB [2009] Created
!>      UB [2010-03-22] Offsets of the intersections in the old/new blocks
!>                      are calculated & saved. Zero-length blocks are
!>                      skipped.
!> \note This routine is used in the counting and sending loops in
!>       dbcsr_complete_redistribute
! *****************************************************************************
  SUBROUTINE dbcsr_reblocking_targets (ints, numints, n_src_dsts,&
       src_sizes, dst_sizes, error)
    INTEGER, INTENT(INOUT)                   :: numints
    INTEGER, DIMENSION(4, numints), &
      INTENT(OUT)                            :: ints
    INTEGER, DIMENSION(:, :), INTENT(OUT)    :: n_src_dsts
    INTEGER, DIMENSION(:), INTENT(IN)        :: src_sizes, dst_sizes
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_reblocking_targets', &
      routineP = moduleN//':'//routineN

    INTEGER :: common_extent, current_dst, current_int, current_src, dst_off, &
      n_dst, n_src, s_dst, s_src, src_off

!   ---------------------------------------------------------------------------

    n_src = SIZE(src_sizes)
    n_dst = SIZE(dst_sizes)
    current_int = 0
    current_src = 0
    s_src = 0 ! HUGE(0)
    DO WHILE (s_src .EQ. 0 .AND. current_src .LE. n_src)
       current_src = current_src + 1
       src_off = 1
       IF (current_src .LE. n_src) THEN
          s_src = src_sizes (current_src)
          n_src_dsts(:, current_src) = (/ current_int+1, 0 /)
       ENDIF
    END DO
    current_dst = 0
    s_dst = 0 ! HUGE(0)
    DO WHILE (s_dst .EQ. 0 .AND. current_dst .LE. n_dst)
       current_dst = current_dst + 1
       dst_off = 1
       IF (current_dst .LE. n_dst) s_dst = dst_sizes (current_dst)
    END DO
    current_int = current_int + 1
    DO WHILE (current_src .LE. n_src .AND. current_dst .LE. n_dst)
       CALL dbcsr_assert (current_int, "LE", numints, &
            dbcsr_fatal_level, dbcsr_internal_error, routineN,&
            "Ran out of space.", __LINE__, error=error)
       ! Calculate how many elements do the current blocks have in
       ! common and record these as going to the current target block.
       common_extent = MIN (s_src, s_dst)
       ints(1, current_int) = current_dst ! target block
       ints(2, current_int) = common_extent
       ints(3, current_int) = src_off
       ints(4, current_int) = dst_off
       ! We've used up the common extents.
       s_src = s_src - common_extent
       s_dst = s_dst - common_extent
       src_off = src_off + common_extent
       dst_off = dst_off + common_extent
       ! We've used up another block.
       n_src_dsts(2,current_src) = n_src_dsts(2,current_src) + 1
       ! Check if we've used up the whole source block.
       DO WHILE (s_src .EQ. 0 .AND. current_src .LE. n_src)
          current_src = current_src + 1
          src_off = 1
          IF (current_src .LE. n_src) THEN
             s_src = src_sizes (current_src)
             n_src_dsts(:, current_src) = (/ current_int+1, 0 /)
          ENDIF
       END DO
       DO WHILE (s_dst .EQ. 0 .AND. current_dst .LE. n_dst)
          current_dst = current_dst + 1
          dst_off = 1
          IF (current_dst .LE. n_dst) s_dst = dst_sizes (current_dst)
       END DO
       current_int = current_int + 1
    ENDDO
    IF (current_src .LT. n_src) &
         n_src_dsts(:, current_src+1:n_src) = -7
    numints = current_int - 1
  END SUBROUTINE dbcsr_reblocking_targets



! *****************************************************************************
!> \brief Distributes elements into blocks and into bins
!>
!> \param[out] block_distribution       block distribution to bins
!> \param[out] block_size       sizes of blocks
!> \param[in] nelements number of elements to bin
!> \param[in] nbins             number of bins
!> \par Term clarification
!>      An example: blocks are atom blocks and bins are process rows/colums.
! *****************************************************************************
  SUBROUTINE create_bl_distribution (block_distribution,&
       block_size, nelements, nbins)
    INTEGER, DIMENSION(:), INTENT(OUT), &
      POINTER                                :: block_distribution, block_size
    INTEGER, INTENT(IN)                      :: nelements, nbins

    CHARACTER(len=*), PARAMETER :: routineN = 'create_bl_distribution', &
      routineP = moduleN//':'//routineN

    INTEGER :: bin, blk_layer, element_stack, els, estimated_blocks, &
      max_blocks_per_bin, nblks, nblocks, stat
    INTEGER, DIMENSION(:), POINTER           :: blk_dist, blk_sizes
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    NULLIFY (block_distribution)
    NULLIFY (block_size)
    ! Define the sizes on which we build the distribution.
    IF (nelements .GT. 0) THEN

       nblocks = CEILING(REAL(nelements)/REAL(max_elements_per_block))
       max_blocks_per_bin = CEILING(REAL(nblocks)/REAL(nbins))

       IF (debug_mod) THEN
          WRITE(*,'(1X,A,1X,A,I7,A,I7,A)')routineP,"For",nelements,&
               " elements and",nbins," bins"
          WRITE(*,'(1X,A,1X,A,I7,A)')routineP,"There are",&
               max_elements_per_block," max elements per block"
          WRITE(*,'(1X,A,1X,A,I7,A)')routineP,"There are",&
               nblocks," blocks"
          WRITE(*,'(1X,A,1X,A,I7,A)')routineP,"There are",&
               max_blocks_per_bin," max blocks/bin"
       ENDIF

       estimated_blocks = max_blocks_per_bin * nbins
       ALLOCATE (blk_dist (estimated_blocks), stat=stat)
       CALL dbcsr_assert (stat == 0, dbcsr_fatal_level, &
            dbcsr_caller_error, routineN,&
            "blk_dist",__LINE__,error)
       ALLOCATE (blk_sizes (estimated_blocks), stat=stat)
       CALL dbcsr_assert (stat == 0,dbcsr_fatal_level, &
            dbcsr_caller_error, routineN,&
            "blk_sizes",__LINE__,error)
       element_stack = 0
       nblks = 0
       DO blk_layer = 1, max_blocks_per_bin
          DO bin = 0, nbins-1
             els = MIN (max_elements_per_block, nelements-element_stack)
             IF (els .GT. 0) THEN
                element_stack = element_stack + els
                nblks = nblks + 1
                blk_dist(nblks) = bin
                blk_sizes(nblks) = els
                IF (debug_mod) WRITE(*,'(1X,A,I5,A,I5,A,I5)')routineP//" Assigning",&
                     els," elements as block",nblks," to bin",bin
             ENDIF
          ENDDO
       ENDDO
       ! Create the output arrays.
       IF (nblks .EQ. estimated_blocks) THEN
          block_distribution => blk_dist
          block_size => blk_sizes
       ELSE
          ALLOCATE (block_distribution (nblks), stat=stat)
          CALL dbcsr_assert (stat == 0, dbcsr_fatal_level, &
               dbcsr_caller_error, routineN,&
               "blk_dist",__LINE__,error)
          block_distribution(:)=blk_dist(1:nblks)
          DEALLOCATE (blk_dist)
          ALLOCATE (block_size (nblks), stat=stat)
          CALL dbcsr_assert (stat == 0,dbcsr_fatal_level, &
               dbcsr_caller_error, routineN,&
               "blk_sizes",__LINE__,error)
          block_size(:)=blk_sizes(1:nblks)
          DEALLOCATE(blk_sizes)
       ENDIF
    ELSE
       ALLOCATE (block_distribution (0), stat=stat)
       CALL dbcsr_assert (stat == 0, dbcsr_fatal_level,&
            dbcsr_caller_error, routineN,&
            "blk_dist",__LINE__,error)
       ALLOCATE (block_size (0), stat=stat)
       CALL dbcsr_assert (stat == 0, dbcsr_fatal_level,&
            dbcsr_caller_error, routineN,&
            "blk_sizes",__LINE__,error)
    ENDIF
1579 FORMAT (I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5)
    IF (debug_mod) THEN
       WRITE(*,'(1X,A,A)')routineP//" Distribution"
       WRITE(*,1579) block_distribution(:)
       WRITE(*,'(1X,A,A)')routineP//" Sizes"
       WRITE(*,1579) block_size(:)
    ENDIF
  END SUBROUTINE create_bl_distribution

! *****************************************************************************
!> \brief Creates a new distribution for the right matrix in a matrix
!>        multiplication with unrotated grid.
!> \param[out] dist_right     new distribution for the right matrix
!> \param[in] dist_left       the distribution of the left matrix
!> \param[in] ncolumns        number of columns in right matrix
!> \param[out] right_col_blk_sizes      sizes of blocks in the created column
!> \par The new row distribution for the right matrix is the same as the row
!>      distribution of the left matrix, while the column distribution is
!>      created so that it is appropriate to the parallel environment.
! *****************************************************************************
  SUBROUTINE dbcsr_create_dist_r_unrot(dist_right, dist_left, ncolumns,&
       right_col_blk_sizes)
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT)                            :: dist_right
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist_left
    INTEGER, INTENT(IN)                      :: ncolumns
    INTEGER, DIMENSION(:), INTENT(OUT), &
      POINTER                                :: right_col_blk_sizes

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_create_dist_r_unrot', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: multiplicity, new_npcols, &
                                                new_nprows, nimages
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: tmp_images
    INTEGER, DIMENSION(:), POINTER           :: right_col_dist, right_row_dist
    TYPE(dbcsr_mp_obj)                       :: right_mp_env

!   ---------------------------------------------------------------------------

    right_mp_env = dbcsr_distribution_mp (dist_left)
    ! Create the column distribution
    CALL create_bl_distribution (right_col_dist,&
         right_col_blk_sizes, ncolumns, dbcsr_mp_npcols (right_mp_env))
    ! Create an even row distribution.
    ALLOCATE (right_row_dist(dbcsr_distribution_ncols (dist_left)))
    ALLOCATE (tmp_images(dbcsr_distribution_ncols (dist_left)))
    new_nprows = dbcsr_mp_nprows (right_mp_env)
    new_npcols = dbcsr_mp_npcols (right_mp_env)
    nimages = lcm (new_nprows, new_npcols) / new_nprows
    multiplicity = new_nprows / gcd (new_nprows, new_npcols)
    CALL rebin_distribution (right_row_dist, tmp_images,& 
         dbcsr_distribution_col_dist (dist_left),&
         new_nprows, multiplicity, nimages)
    CALL dbcsr_distribution_new (dist_right, right_mp_env,&
         right_row_dist, right_col_dist, &
         reuse_arrays=.TRUE.)
    DEALLOCATE(tmp_images)
  END SUBROUTINE dbcsr_create_dist_r_unrot

! *****************************************************************************
!> \brief Creates a block-cyclic compatible distribution
!>
!>        All blocks in a dimension, except for possibly the last
!>        block, have the same size.
!> \param[out] dist           the elemental distribution
!> \param[in] nrows           number of full rows
!> \param[in] ncolumns        number of full columns
!> \param[in] nrow_block      size of row blocks
!> \param[in] ncol_block      size of column blocks
!> \param[in] mp_env          multiprocess environment
!> \param[out] row_blk_sizes  row block sizes
!> \param[out] col_blk_sizes  column block sizes
! *****************************************************************************
  SUBROUTINE dbcsr_create_dist_block_cyclic (dist, nrows, ncolumns,&
       nrow_block, ncol_block, mp_env, row_blk_sizes, col_blk_sizes)
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT)                            :: dist
    INTEGER, INTENT(IN)                      :: nrows, ncolumns, nrow_block, &
                                                ncol_block
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    INTEGER, DIMENSION(:), INTENT(OUT), &
      POINTER                                :: row_blk_sizes, col_blk_sizes

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'dbcsr_create_dist_block_cyclic', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nblkcols, nblkrows, npcols, &
                                                nprows, pdim, sz
    INTEGER, DIMENSION(:), POINTER           :: cd_data, rd_data

!   ---------------------------------------------------------------------------
! Row sizes

    IF(nrow_block.EQ.0) THEN
       nblkrows = 0
       sz = 0
    ELSE
       nblkrows = nrows / nrow_block
       sz = MOD(nrows, nrow_block)
    ENDIF
    IF (sz .GT. 0) nblkrows = nblkrows + 1
    ALLOCATE (row_blk_sizes(nblkrows), rd_data(nblkrows))
    row_blk_sizes = nrow_block
    IF (sz .NE. 0) row_blk_sizes(nblkrows) = sz

    ! Column sizes
    IF(ncol_block.EQ.0) THEN
       nblkcols = 0
       sz = 0
    ELSE
       nblkcols = ncolumns / ncol_block
       sz = MOD(ncolumns, ncol_block)
    ENDIF
    IF (sz .GT. 0) nblkcols = nblkcols + 1
    ALLOCATE (col_blk_sizes(nblkcols), cd_data(nblkcols))
    col_blk_sizes = ncol_block
    IF (sz .NE. 0) col_blk_sizes(nblkcols) = sz
    !
    IF(debug_mod) THEN
       WRITE(*,*)routineN//" nrows,nrow_block,nblkrows=",&
            nrows,nrow_block,nblkrows
       WRITE(*,*)routineN//" ncols,ncol_block,nblkcols=",&
            ncolumns,ncol_block,nblkcols
    ENDIF
    ! Calculate process row distribution
    nprows = dbcsr_mp_nprows (mp_env)
    DO pdim = 0 , MIN(nprows-1,nblkrows-1)
       rd_data(1+pdim:nblkrows:nprows) = pdim
    END DO
    ! Calculate process column distribution
    npcols = dbcsr_mp_npcols (mp_env)
    DO pdim = 0 , MIN(npcols-1, nblkcols-1)
       cd_data(1+pdim:nblkcols:npcols) = pdim
    END DO
    !
    IF(debug_mod) THEN
       WRITE(*,*)routineN//" row_dist",&
            rd_data
       WRITE(*,*)routineN//" col_dist",&
            cd_data
    ENDIF
    !
    CALL dbcsr_distribution_new (dist, mp_env,&
         rd_data, cd_data, reuse_arrays=.TRUE.)
    !
  END SUBROUTINE dbcsr_create_dist_block_cyclic

! *****************************************************************************
!> \brief Matches row/block sizes and offsets to a given distribution
!> \param old_sizes ...
!> \param mapping ...
!> \param nel_new ...
!> \param new_sizes ...
!> \param error ...
!> \note Used for making matrices dense/undense
! *****************************************************************************
  SUBROUTINE make_sizes_dense (old_sizes, mapping, nel_new, new_sizes, error)
    TYPE(array_i1d_obj), INTENT(IN)          :: old_sizes, mapping
    INTEGER, INTENT(IN)                      :: nel_new
    TYPE(array_i1d_obj), INTENT(OUT)         :: new_sizes
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_sizes_dense', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: el, nel_old
    INTEGER, DIMENSION(:), POINTER           :: map, new_s, old_s

!   ---------------------------------------------------------------------------

    map => array_data (mapping)
    old_s => array_data (old_sizes)
    nel_old = array_size (old_sizes)
    ALLOCATE (new_s (nel_new))
    new_s(:) = 0
    DO el = 1, nel_old
       new_s(map(el)) = new_s(map(el)) + old_s(el)
    ENDDO
    CALL array_new (new_sizes, new_s, gift=.TRUE.)
  END SUBROUTINE make_sizes_dense


! *****************************************************************************
!> \brief Makes a 1-D distribution dense.
!> \param[in] nbins Number of bins in the main distribution
!> \param[in] old_dist         Current distribution
!> \param[out] dense_dist      Dense distribution
!> \param[out] dist_map        Map from current to dense distribution
!> \param[in] nsubdist         Number of bins in the subdistribution
!> \param[in] subdist          Subdistribution
! *****************************************************************************
  SUBROUTINE dbcsr_make_1dist_dense (nbins, old_dist, dense_dist, dist_map, nsubdist, subdist)
    INTEGER, INTENT(IN)                      :: nbins
    TYPE(array_i1d_obj), INTENT(IN)          :: old_dist
    TYPE(array_i1d_obj), INTENT(OUT)         :: dense_dist, dist_map
    INTEGER, INTENT(IN)                      :: nsubdist
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: subdist

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_make_1dist_dense', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: b, i, n_new_bins
    INTEGER, DIMENSION(:), POINTER           :: dense, map, old_d

!   ---------------------------------------------------------------------------

    IF (nsubdist .EQ. 0) THEN
       n_new_bins = nbins
    ELSE
       n_new_bins = nbins * nsubdist
    ENDIF
    old_d => array_data (old_dist)
    ALLOCATE (dense (n_new_bins))
    ALLOCATE (map (array_size (old_dist)))
    !
    IF (nsubdist .EQ. 0) THEN
       dense(:) = (/ (b, b = 0, n_new_bins-1) /)
       map(:) = old_d(:)+1
    ELSE
       DO i = 1, nbins
          dense((i-1)*nsubdist+1 : (i)*nsubdist) = i-1
       ENDDO
       map(:) = old_d(:)*nsubdist + subdist(:) + 1
    ENDIF
    !
    CALL array_new (dense_dist, dense, gift=.TRUE.)
    CALL array_new (dist_map, map, gift=.TRUE.)
  END SUBROUTINE dbcsr_make_1dist_dense


! *****************************************************************************
!> \brief Prepares distributions for making dense matrices.
!> \param product_dist ...
!> \param left_rdist ...
!> \param right_rdist ...
!> \param dense_product_dist ...
!> \param dense_left_rdist ...
!> \param dense_right_rdist ...
!> \param partial ...
!> \param m_map ...
!> \param k_vmap ...
!> \param n_map ...
!> \param old_m_sizes ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_make_dists_dense (product_dist, left_rdist, right_rdist,&
       dense_product_dist, dense_left_rdist, dense_right_rdist,&
       partial,&
       m_map, k_vmap, n_map,&
       old_m_sizes,&
       error)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: product_dist
    TYPE(dbcsr_imagedistribution_obj), &
      INTENT(IN)                             :: left_rdist, right_rdist
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT)                            :: dense_product_dist
    TYPE(dbcsr_imagedistribution_obj), &
      INTENT(OUT)                            :: dense_left_rdist, &
                                                dense_right_rdist
    LOGICAL, INTENT(IN)                      :: partial
    TYPE(array_i1d_obj), INTENT(OUT)         :: m_map, k_vmap, n_map
    TYPE(array_i1d_obj), INTENT(IN)          :: old_m_sizes
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_make_dists_dense', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle, i, j, k_nbins, &
                                                m_nbins, n_nbins, nthreads
    INTEGER, DIMENSION(:), POINTER           :: tdist
    TYPE(array_i1d_obj) :: new_k_idist, new_k_pdist, new_k_vdist, new_m_dist, &
      new_m_sizes, new_n_dist, old_k_vdist, old_m_dist, old_n_dist
    TYPE(dbcsr_distribution_obj)             :: dense_left_dist, &
                                                dense_right_dist

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handle, error)
    !
    CALL dbcsr_assert (dbcsr_distribution_has_threads (product_dist), &
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Product distribution must have threads.", __LINE__, error=error)
    tdist => array_data (dbcsr_distribution_thread_dist (product_dist))
    old_m_dist = product_dist%d%row_dist_block
    old_n_dist = product_dist%d%col_dist_block
    old_k_vdist = right_rdist%i%vrow_dist
    m_nbins = dbcsr_mp_nprows (product_dist%d%mp_env)
    n_nbins = dbcsr_mp_npcols (product_dist%d%mp_env)
    k_nbins = dbcsr_mp_nprows (right_rdist%i%main%d%mp_env) * right_rdist%i%row_decimation
    CALL dbcsr_assert (array_equality (old_k_vdist, left_rdist%i%vcol_dist),&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "k distribution mismtach", __LINE__, error=error)
    nthreads = product_dist%d%num_threads
    !
    IF (partial) THEN
       new_m_dist = old_m_dist
       CALL array_hold (new_m_dist)
       new_n_dist = old_n_dist
       CALL array_hold (new_n_dist)
       dense_product_dist = product_dist
       CALL dbcsr_distribution_hold (dense_product_dist)
       CALL array_new (m_map, (/ (i, i=1, array_size(new_m_dist)) /), lb=1)
       CALL array_new (n_map, (/ (i, i=1, array_size(new_n_dist)) /), lb=1)
    ELSE
       CALL dbcsr_make_1dist_dense (m_nbins, old_m_dist, new_m_dist, m_map,&
            nthreads, tdist)
       CALL dbcsr_make_1dist_dense (n_nbins, old_n_dist, new_n_dist, n_map, 0)
       CALL dbcsr_distribution_new (dense_product_dist, product_dist%d%mp_env,&
            array_data(new_m_dist), array_data(new_n_dist))
       CALL make_sizes_dense (old_m_sizes, m_map, array_size(new_m_dist), new_m_sizes, error=error)
       CALL dbcsr_distribution_make_threads (dense_product_dist,&
            array_data(new_m_sizes))
       CALL array_release (new_m_sizes)
       tdist => array_data (dbcsr_distribution_thread_dist (dense_product_dist))
       ! Resets the thread distribution to be in-order.
       DO i = 1, m_nbins
          tdist((i-1)*nthreads+1 : (i)*nthreads) = (/ (j, j=0,nthreads-1) /)
       ENDDO
    ENDIF
    !
    CALL dbcsr_make_1dist_dense (k_nbins, old_k_vdist, new_k_vdist, k_vmap, 0)
    CALL v_to_p_i_dist_o (new_k_vdist,&
         left_rdist%i%col_decimation, new_k_pdist, new_k_idist)
    ! Left
    CALL dbcsr_distribution_new (dense_left_dist, left_rdist%i%main%d%mp_env,&
         array_data(new_m_dist), array_data(new_k_pdist))
    CALL dbcsr_distribution_no_threads (dense_left_dist)
    dense_left_dist%d%thread_dist = dbcsr_distribution_thread_dist (dense_product_dist)
    CALL array_hold (dense_left_dist%d%thread_dist)
    dense_left_dist%d%has_thread_dist = .TRUE.
    CALL dbcsr_new_image_dist (dense_left_rdist, dense_left_dist, left_rdist, error)
    CALL dbcsr_distribution_release (dense_left_dist)
    CALL array_new (dense_left_rdist%i%row_image,&
         (/ (1, i=1, array_size (new_m_dist)) /), lb=1)
    dense_left_rdist%i%col_image = new_k_idist
    CALL array_hold (new_k_idist)
    dense_left_rdist%i%vrow_dist = new_m_dist
    CALL array_hold (new_m_dist)
    dense_left_rdist%i%vcol_dist = new_k_vdist
    CALL array_hold (new_k_vdist)
    !
    CALL array_release (new_k_pdist)
    CALL array_release (new_k_idist)
    ! Right
    CALL v_to_p_i_dist_o (new_k_vdist,&
         right_rdist%i%row_decimation, new_k_pdist, new_k_idist)
    CALL dbcsr_distribution_new (dense_right_dist, right_rdist%i%main%d%mp_env,&
         array_data(new_k_pdist), array_data(new_n_dist))
    CALL dbcsr_new_image_dist (dense_right_rdist, dense_right_dist, right_rdist, error)
    CALL dbcsr_distribution_release (dense_right_dist)
    CALL array_new (dense_right_rdist%i%col_image,&
         (/ (1, i=1, array_size (new_n_dist)) /), lb=1)
    dense_right_rdist%i%row_image = new_k_idist
    CALL array_hold (new_k_idist)
    dense_right_rdist%i%vrow_dist = new_k_vdist
    CALL array_hold (new_k_vdist)
    dense_right_rdist%i%vcol_dist = new_n_dist
    CALL array_hold (new_n_dist)
    !
    CALL array_release (new_k_idist)
    CALL array_release (new_k_pdist)
    CALL array_release (new_m_dist)
    CALL array_release (new_n_dist)
    CALL array_release (new_k_vdist)
    !
    CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_make_dists_dense


! *****************************************************************************
!> \brief Finds the block that has the given column.
!>
!>        If the block having the queried column is found, the blk parameter
!>        is set to this block number and the found parameter is true.
!>        Otherwise found is false and the block number is invalid.
!> \param[in] find_col        column to find
!> \param[in] frst_blk        first block number in row
!> \param[in] last_blk        last block number in row
!> \param[in] col_i           col indices
!> \param[in] blk_p           block pointers
!> \param[out] blk            block number with searched-for column
!> \param[out] found          flag specified whether a block that has the
!>                            correct column was found
!> \par Index validity
!>      The blk_p array of block pointers is a required parameter to enable
!>      the detection of deleted blocks.
! *****************************************************************************
  PURE SUBROUTINE dbcsr_find_column(find_col,frst_blk,last_blk,col_i,blk_p,&
       blk,found)
    INTEGER, INTENT(IN)                      :: find_col, frst_blk, last_blk
    INTEGER, DIMENSION(:), INTENT(IN)        :: col_i, blk_p
    INTEGER, INTENT(OUT)                     :: blk
    LOGICAL, INTENT(OUT)                     :: found

    INTEGER                                  :: col, high, low

    found = .FALSE.
    low = frst_blk
    high = last_blk
    blk = (low + high) / 2
    DO WHILE (blk .GE. low .AND. blk .LE. high)
       col = col_i(blk)
       IF(col.EQ.find_col) THEN
          found =  blk_p(blk) .NE. 0
          EXIT
       ELSEIF (col .LT. find_col) THEN
          low = blk+1
       ELSE
          high = blk-1
       ENDIF
       blk = (low + high) / 2
    ENDDO
  END SUBROUTINE dbcsr_find_column

! *****************************************************************************
!> \brief Finds the local virtual elements
!>
!> All elements are mapped at once.  Therefore an entry in the
!> resulting array can be used as a lookup index for any of the local
!> element arrays.  The distribution itself tells into which array to
!> look.
!> \param[out] local_elements local virtual elements
!> \param[in] bin_distribution           distribution of elements to bins
!> \param[in] nbins                      number of bins
! *****************************************************************************
  SUBROUTINE find_all_local_elements (local_elements, &
       bin_distribution, nbins)
    INTEGER, INTENT(IN)                      :: nbins
    INTEGER, DIMENSION(:), INTENT(IN)        :: bin_distribution
    TYPE(array_i1d_obj), &
      DIMENSION(0:nbins-1), INTENT(INOUT)    :: local_elements

    INTEGER                                  :: bin, el
    INTEGER, DIMENSION(0:nbins-1)            :: nlve

    nlve(:) = 0
    DO el = 1, SIZE (bin_distribution)
       bin = bin_distribution (el)
       nlve(bin) = nlve(bin) + 1
       local_elements(bin)%low%data(nlve(bin)) = el
    ENDDO
  END SUBROUTINE find_all_local_elements

! *****************************************************************************
!> \brief Points to local virtual elements.
!>
!> All elements are mapped at once.  Therefore an entry in the
!> resulting array can be used as a lookup index for any of the local
!> element arrays.  The distribution itself tells into which array to
!> look.
!> \param[out] global_elements   enumerated local elements
!> \param[in] bin_distribution           distribution of elements to bins
!> \param[in] nbins                      number of bins
! *****************************************************************************
  PURE SUBROUTINE map_all_local_elements (global_elements, &
       bin_distribution, nbins)
    INTEGER, DIMENSION(:), INTENT(OUT)       :: global_elements
    INTEGER, DIMENSION(:), INTENT(IN)        :: bin_distribution
    INTEGER, INTENT(IN)                      :: nbins

    INTEGER                                  :: bin, el
    INTEGER, DIMENSION(0:nbins-1)            :: nlve

    nlve(:) = 0
    DO el = 1, SIZE (bin_distribution)
       bin = bin_distribution (el)
       nlve(bin) = nlve(bin) + 1
       global_elements(el) = nlve(bin)
    ENDDO
  END SUBROUTINE map_all_local_elements


! *****************************************************************************
!> \brief Determines mapping from local to global rows
!> \param[in] dist              mapping for this distribution
!> \param[in] local_rows        local elements for specified row
!> \param[in] local_prow find local elemnts for this local row
!> \param[in,out] error         error
! *****************************************************************************
  SUBROUTINE dbcsr_get_local_rows (dist, local_rows, local_prow, error)
    TYPE(dbcsr_distribution_obj), &
      INTENT(INOUT)                          :: dist
    TYPE(array_i1d_obj), INTENT(OUT)         :: local_rows
    INTEGER, INTENT(IN)                      :: local_prow
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_local_rows', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: el, error_handle, nprows, prow
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: nle
    INTEGER, DIMENSION(:), POINTER           :: itmp, row_dist

    CALL dbcsr_error_set(routineN, error_handle, error)
    ! If the current local row mappings do not exist, create them.
    IF (.NOT. dist%d%has_other_l_rows) THEN
       dist%d%has_other_l_rows = .TRUE.
       nprows = dbcsr_mp_nprows (dbcsr_distribution_mp (dist))
       ALLOCATE (dist%d%other_l_rows(0:dbcsr_mp_nprows (dist%d%mp_env)-1))
       ALLOCATE (nle (0:nprows-1))
       row_dist => dbcsr_distribution_row_dist (dist)
       ! Count the number of local elements per row.
       nle(:) = 0
       DO el = 1, SIZE (row_dist)
          prow = row_dist(el)
          nle(prow) = nle(prow) + 1
       ENDDO
       DO prow = 0, nprows-1
          ALLOCATE(itmp(nle(prow)))
          itmp=0
          CALL array_new (dist%d%other_l_rows(prow),&
               itmp, lb=1)
          DEALLOCATE(itmp)
       ENDDO
       DEALLOCATE (nle)
       CALL find_all_local_elements (dist%d%other_l_rows, row_dist, nprows)
    ELSE
       IF (careful_mod) THEN
          CALL dbcsr_assert (ASSOCIATED (dist%d%other_l_rows),&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "Local rows mapping does not exist.", __LINE__, error=error)
       ENDIF
    ENDIF
    local_rows = dist%d%other_l_rows(local_prow)
    CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_get_local_rows

! *****************************************************************************
!> \brief Determines mapping from local to global columnss
!> \param[in] dist              mapping for this distribution
!> \param[in] local_cols        local elements for specified column
!> \param[in] local_pcol         find local elemnts for this local column 
!> \param[in,out] error         error
! *****************************************************************************
  SUBROUTINE dbcsr_get_local_cols (dist, local_cols, local_pcol, error)
    TYPE(dbcsr_distribution_obj), &
      INTENT(INOUT)                          :: dist
    TYPE(array_i1d_obj), INTENT(OUT)         :: local_cols
    INTEGER, INTENT(IN)                      :: local_pcol
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_local_cols', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: el, error_handle, npcols, pcol
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: nle
    INTEGER, DIMENSION(:), POINTER           :: col_dist, itmp

    CALL dbcsr_error_set(routineN, error_handle, error)
    ! If the current local col mappings do not exist, create them.
    IF (.NOT. dist%d%has_other_l_cols) THEN
       dist%d%has_other_l_cols = .TRUE.
       npcols = dbcsr_mp_npcols (dbcsr_distribution_mp (dist))
       ALLOCATE (dist%d%other_l_cols(0:dbcsr_mp_npcols (dist%d%mp_env)-1))
       ALLOCATE (nle (0:npcols-1))
       col_dist => dbcsr_distribution_col_dist (dist)
       ! Count the number of local elements per col.
       nle(:) = 0
       DO el = 1, SIZE (col_dist)
          pcol = col_dist(el)
          nle(pcol) = nle(pcol) + 1
       ENDDO
       DO pcol = 0, npcols-1
          ALLOCATE(itmp(nle(pcol)))
          itmp=0
          CALL array_new (dist%d%other_l_cols(pcol),&
               itmp, lb=1)
          DEALLOCATE(itmp)
       ENDDO
       DEALLOCATE (nle)
       CALL find_all_local_elements (dist%d%other_l_cols, col_dist, npcols)
    ELSE
       IF (careful_mod) THEN
          CALL dbcsr_assert (ASSOCIATED (dist%d%other_l_cols),&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "Local columns mapping does not exist.", __LINE__, error=error)
       ENDIF
    ENDIF
    local_cols = dist%d%other_l_cols(local_pcol)
    CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_get_local_cols


! *****************************************************************************
!> \brief Determines mapping from global to local rows
!> \param[in] dist              mapping for this distribution
!> \param[in] row_map           mapping to local rows
!> \param[in,out] error         error
! *****************************************************************************
  SUBROUTINE dbcsr_get_global_row_map (dist, row_map, error)
    TYPE(dbcsr_distribution_obj), &
      INTENT(INOUT)                          :: dist
    TYPE(array_i1d_obj), INTENT(OUT)         :: row_map
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_global_row_map', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle, nprows
    INTEGER, DIMENSION(:), POINTER           :: rmap, row_dist

    CALL dbcsr_error_set(routineN, error_handle, error)
    ! If the current local row mappings do not exist, create them.
    IF (.NOT. dist%d%has_global_row_map) THEN
       row_dist => dbcsr_distribution_row_dist (dist)
       ALLOCATE (rmap (SIZE (row_dist)))
       nprows = dbcsr_mp_nprows (dbcsr_distribution_mp (dist))
       CALL map_all_local_elements (rmap, row_dist, nprows)
       CALL array_new (dist%d%global_row_map, rmap, gift=.TRUE.)
       dist%d%has_global_row_map = .TRUE.
    ELSE
       IF (careful_mod) THEN
          CALL dbcsr_assert (array_exists (dist%d%global_row_map),&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "Row map does not exist.", __LINE__, error=error)
       ENDIF
    ENDIF
    row_map = dist%d%global_row_map
    CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_get_global_row_map


! *****************************************************************************
!> \brief Determines mapping from global to local columns
!> \param[in] dist              mapping for this distribution
!> \param[in] col_map           mapping to local columns
!> \param[in,out] error         error
! *****************************************************************************
  SUBROUTINE dbcsr_get_global_col_map (dist, col_map, error)
    TYPE(dbcsr_distribution_obj), &
      INTENT(INOUT)                          :: dist
    TYPE(array_i1d_obj), INTENT(OUT)         :: col_map
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_global_col_map', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle, npcols
    INTEGER, DIMENSION(:), POINTER           :: cmap, col_dist

    CALL dbcsr_error_set(routineN, error_handle, error)
    ! If the current local col mappings do not exist, create them.
    IF (.NOT. dist%d%has_global_col_map) THEN
       col_dist => dbcsr_distribution_col_dist (dist)
       ALLOCATE (cmap (SIZE (col_dist)))
       npcols = dbcsr_mp_npcols (dbcsr_distribution_mp (dist))
       CALL map_all_local_elements (cmap, col_dist, npcols)
       CALL array_new (dist%d%global_col_map, cmap, gift=.TRUE.)
       dist%d%has_global_col_map = .TRUE.
    ELSE
       IF (careful_mod) THEN
          CALL dbcsr_assert (array_exists (dist%d%global_col_map),&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "Column map does not exist.", __LINE__, error=error)
       ENDIF
    ENDIF
    col_map = dist%d%global_col_map
    CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_get_global_col_map


! *****************************************************************************
!> \brief Resets local rows, columns to the correct arrays and values.
!> \param matrix ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_reset_locals (matrix, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    TYPE(dbcsr_distribution_obj)             :: dist

!   ---------------------------------------------------------------------------

    dist = dbcsr_distribution (matrix)
    CALL dbcsr_release_locals (matrix, error=error)
    ! Rows
    IF (dbg) &
         WRITE(*,*)"reset local rows for ", TRIM(matrix%m%name),&
         matrix%m%nblkrows_local, "prow", matrix%m%index(dbcsr_slot_home_prow),&
         dbcsr_mp_myprow (dbcsr_distribution_mp (matrix%m%dist))
    CALL dbcsr_get_local_rows (dist, matrix%m%local_rows,&
         matrix%m%index(dbcsr_slot_home_prow), error=error)
    CALL array_hold (matrix%m%local_rows)
    IF (dbg) WRITE(*,*)"local rows", matrix%m%local_rows%low%data
    matrix%m%nblkrows_local = array_size (matrix%m%local_rows)
    CALL dbcsr_get_global_row_map (dist, matrix%m%global_rows, error=error)
    CALL array_hold (matrix%m%global_rows)
    matrix%m%has_local_rows = .TRUE.
    matrix%m%has_global_rows = .TRUE.
    ! Columns
    IF (dbg) &
         WRITE(*,*)"reset local cols for ", TRIM(matrix%m%name),&
         matrix%m%nblkcols_local, "pcol", matrix%m%index(dbcsr_slot_home_pcol),&
         dbcsr_mp_mypcol (dbcsr_distribution_mp (matrix%m%dist))
    CALL dbcsr_get_local_cols (dist, matrix%m%local_cols,&
         matrix%m%index(dbcsr_slot_home_pcol), error=error)
    CALL array_hold (matrix%m%local_cols)
    IF (dbg) WRITE(*,*)"local cols", matrix%m%local_cols%low%data
    matrix%m%nblkcols_local = array_size (matrix%m%local_cols)
    CALL dbcsr_get_global_col_map (dist, matrix%m%global_cols, error=error)
    CALL array_hold (matrix%m%global_cols)
    matrix%m%has_local_cols = .TRUE.
    matrix%m%has_global_cols = .TRUE.
    !
  END SUBROUTINE dbcsr_reset_locals

! *****************************************************************************
!> \brief Determines mapping from local to global virtual process rows
!> \param[in] imgdist           image distribution
!> \param[out] local_vrows       local rows 
!> \param[in] local_vprow       the local virtual process row
!> \param[in,out] error         error
! *****************************************************************************
  SUBROUTINE dbcsr_get_local_vrows (imgdist, local_vrows, local_vprow, error)
    TYPE(dbcsr_imagedistribution_obj), &
      INTENT(INOUT)                          :: imgdist
    TYPE(array_i1d_obj), INTENT(OUT)         :: local_vrows
    INTEGER, INTENT(IN)                      :: local_vprow
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_local_vrows', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: el, error_handle, nvprows, &
                                                vprow
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: itmp, nle
    INTEGER, DIMENSION(:), POINTER           :: vrow_dist

    IF (careful_mod) CALL dbcsr_error_set(routineN, error_handle, error)
    ! If the current local row mappings do not exist, create them.
    IF (.NOT. imgdist%i%has_other_vl_rows) THEN
       imgdist%i%has_other_vl_rows = .TRUE.
       nvprows = dbcsr_mp_nprows (dbcsr_distribution_mp (imgdist%i%main))&
               * imgdist%i%row_decimation
       ALLOCATE (imgdist%i%other_vl_rows(0:nvprows-1))
       ALLOCATE (nle (0:nvprows-1))
       vrow_dist => array_data (imgdist%i%vrow_dist)
       ! Count the number of local elements per row.
       nle(:) = 0
       DO el = 1, SIZE (vrow_dist)
          vprow = vrow_dist(el)
          nle(vprow) = nle(vprow) + 1
       ENDDO
       DO vprow = 0, nvprows-1
          ALLOCATE(itmp(nle(vprow)))
          itmp=0
          CALL array_new (imgdist%i%other_vl_rows(vprow),&
               itmp, lb=1)
          DEALLOCATE(itmp)
       ENDDO
       DEALLOCATE (nle)
       CALL find_all_local_elements (imgdist%i%other_vl_rows, vrow_dist, nvprows)
    ELSE
       IF (careful_mod) THEN
          CALL dbcsr_assert (ASSOCIATED (imgdist%i%other_vl_rows),&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "Local rows mapping does not exist.", __LINE__, error=error)
       ENDIF
    ENDIF
    local_vrows = imgdist%i%other_vl_rows(local_vprow)
    IF (careful_mod) CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_get_local_vrows

! *****************************************************************************
!> \brief Determines mapping from local to global virtual process columnss
!> \param[in] imgdist           image distribution
!> \param[out] local_vcols       local columnss 
!> \param[in] local_vpcol       the local virtual process column
!> \param[in,out] error         error
! *****************************************************************************
  SUBROUTINE dbcsr_get_local_vcols (imgdist, local_vcols, local_vpcol, error)
    TYPE(dbcsr_imagedistribution_obj), &
      INTENT(INOUT)                          :: imgdist
    TYPE(array_i1d_obj), INTENT(OUT)         :: local_vcols
    INTEGER, INTENT(IN)                      :: local_vpcol
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_local_vcols', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: el, error_handle, nvpcols, &
                                                vpcol
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: nle
    INTEGER, DIMENSION(:), POINTER           :: itmp, vcol_dist

    IF (careful_mod) CALL dbcsr_error_set(routineN, error_handle, error)
    ! If the current local col mappings do not exist, create them.
    IF (.NOT. imgdist%i%has_other_vl_cols) THEN
       imgdist%i%has_other_vl_cols = .TRUE.
       nvpcols = dbcsr_mp_npcols (dbcsr_distribution_mp (imgdist%i%main))&
               * imgdist%i%col_decimation
       ALLOCATE (imgdist%i%other_vl_cols(0:nvpcols-1))
       ALLOCATE (nle (0:nvpcols-1))
       vcol_dist => array_data (imgdist%i%vcol_dist)
       ! Count the number of local elements per col.
       nle(:) = 0
       DO el = 1, SIZE (vcol_dist)
          vpcol = vcol_dist(el)
          nle(vpcol) = nle(vpcol) + 1
       ENDDO
       DO vpcol = 0, nvpcols-1
          ALLOCATE(itmp(nle(vpcol)))
          itmp=0
          CALL array_new (imgdist%i%other_vl_cols(vpcol),&
               itmp, lb=1)
          DEALLOCATE(itmp)
       ENDDO
       DEALLOCATE (nle)
       CALL find_all_local_elements (imgdist%i%other_vl_cols, vcol_dist, nvpcols)
    ELSE
       IF (careful_mod) THEN
          CALL dbcsr_assert (ASSOCIATED (imgdist%i%other_vl_cols),&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "Local cols mapping does not exist.", __LINE__, error=error)
       ENDIF
    ENDIF
    local_vcols = imgdist%i%other_vl_cols(local_vpcol)
    IF (careful_mod) CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_get_local_vcols

! *****************************************************************************
!> \brief Determines mapping from global to virtual local rows
!> \param[in] imgdist           mapping for this image distribution
!> \param[in] vrow_map           mapping to local rows
!> \param[in,out] error         error
! *****************************************************************************
  SUBROUTINE dbcsr_get_global_vrow_map (imgdist, vrow_map, error)
    TYPE(dbcsr_imagedistribution_obj), &
      INTENT(INOUT)                          :: imgdist
    TYPE(array_i1d_obj), INTENT(OUT)         :: vrow_map
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_global_vrow_map', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle, nvprows
    INTEGER, DIMENSION(:), POINTER           :: rmap, vrow_dist

    IF (careful_mod) CALL dbcsr_error_set(routineN, error_handle, error)
    ! If the current local row mappings do not exist, create them.
    IF (.NOT. imgdist%i%has_global_vrow_map) THEN
       vrow_dist => array_data (imgdist%i%vrow_dist)
       ALLOCATE (rmap (SIZE (vrow_dist)))
       nvprows = dbcsr_mp_nprows (dbcsr_distribution_mp (imgdist%i%main))&
               * imgdist%i%row_decimation
       CALL map_all_local_elements (rmap, vrow_dist, nvprows)
       CALL array_new (imgdist%i%global_vrow_map, rmap, gift=.TRUE.)
       imgdist%i%has_global_vrow_map = .TRUE.
    ELSE
       IF (careful_mod) THEN
          CALL dbcsr_assert (array_exists (imgdist%i%global_vrow_map),&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "Row map does not exist.", __LINE__, error=error)
       ENDIF
    ENDIF
    vrow_map = imgdist%i%global_vrow_map
    IF (careful_mod) CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_get_global_vrow_map

! *****************************************************************************
!> \brief Determines mapping from global to virtual local columnss
!> \param[in] imgdist           mapping for this image distribution
!> \param[in] vcol_map           mapping to local columnss 
!> \param[in,out] error         error
! *****************************************************************************
  SUBROUTINE dbcsr_get_global_vcol_map (imgdist, vcol_map, error)
    TYPE(dbcsr_imagedistribution_obj), &
      INTENT(INOUT)                          :: imgdist
    TYPE(array_i1d_obj), INTENT(OUT)         :: vcol_map
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_global_vcol_map', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle, nvpcols
    INTEGER, DIMENSION(:), POINTER           :: rmap, vcol_dist

    IF (careful_mod) CALL dbcsr_error_set(routineN, error_handle, error)
    ! If the current local col mappings do not exist, create them.
    IF (.NOT. imgdist%i%has_global_vcol_map) THEN
       vcol_dist => array_data (imgdist%i%vcol_dist)
       ALLOCATE (rmap (SIZE (vcol_dist)))
       nvpcols = dbcsr_mp_npcols (dbcsr_distribution_mp (imgdist%i%main))&
               * imgdist%i%col_decimation
       CALL map_all_local_elements (rmap, vcol_dist, nvpcols)
       CALL array_new (imgdist%i%global_vcol_map, rmap, gift=.TRUE.)
       imgdist%i%has_global_vcol_map = .TRUE.
    ELSE
       IF (careful_mod) THEN
          CALL dbcsr_assert (array_exists (imgdist%i%global_vcol_map),&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "Col map does not exist.", __LINE__, error=error)
       ENDIF
    ENDIF
    vcol_map = imgdist%i%global_vcol_map
    IF (careful_mod) CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_get_global_vcol_map



! *****************************************************************************
!> \brief Resets local rows, columns to the correct arrays and values
!>         for images.
!> \param matrix ...
!> \param imgdist ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_reset_vlocals (matrix, imgdist, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    TYPE(dbcsr_imagedistribution_obj), &
      INTENT(INOUT)                          :: imgdist
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    TYPE(dbcsr_distribution_obj)             :: dist

!   ---------------------------------------------------------------------------

    dist = dbcsr_distribution (matrix)
    ! Rows
    CALL dbcsr_release_locals (matrix, error)
    CALL dbcsr_get_local_vrows (imgdist, matrix%m%local_rows,&
         matrix%m%index(dbcsr_slot_home_vprow), error=error)
    CALL array_hold (matrix%m%local_rows)
    matrix%m%has_local_rows = .TRUE.
    !write(*,'(1X,A,7(1X,I7))')"reset_vlocals my homevprow",&
    !     matrix%m%index(dbcsr_slot_home_vprow), &
    !     matrix%m%index(dbcsr_slot_home_prow), matrix%m%index(dbcsr_slot_home_rowi),&
    !     matrix%m%dist%d%mp_env%mp%myprow, matrix%m%dist%d%mp_env%mp%mypcol,&
    !     SIZE(matrix%m%dist%d%mp_env%mp%pgrid,1), SIZE(matrix%m%dist%d%mp_env%mp%pgrid,2)
    matrix%m%nblkrows_local = array_size (matrix%m%local_rows)
    matrix%m%index(dbcsr_slot_nblkrows_local) = array_size (matrix%m%local_rows)
    CALL dbcsr_get_global_vrow_map (imgdist, matrix%m%global_rows, error=error)
    CALL array_hold (matrix%m%global_rows)
    matrix%m%has_global_rows = .TRUE.
    ! Columns
    CALL dbcsr_get_local_vcols (imgdist, matrix%m%local_cols,&
         matrix%m%index(dbcsr_slot_home_vpcol), error=error)
    CALL array_hold (matrix%m%local_cols)
    matrix%m%has_local_cols = .TRUE.
    matrix%m%nblkcols_local = array_size (matrix%m%local_cols)
    matrix%m%index(dbcsr_slot_nblkcols_local) = array_size (matrix%m%local_cols)
    CALL dbcsr_get_global_vcol_map (imgdist, matrix%m%global_cols, error=error)
    CALL array_hold (matrix%m%global_cols)
    matrix%m%has_global_cols = .TRUE.
  END SUBROUTINE dbcsr_reset_vlocals

! *****************************************************************************
!> \brief Converts virtual 1-D distribution to process and image
!> \param nel ...
!> \param vdist ...
!> \param nim ...
!> \param pdist ...
!> \param idist ...
! *****************************************************************************
  pure SUBROUTINE v_to_p_i_dist (nel, vdist, nim, pdist, idist)
    INTEGER, INTENT(in)                      :: nel
    INTEGER, DIMENSION(1:nel), INTENT(in)    :: vdist
    INTEGER, INTENT(in)                      :: nim
    INTEGER, DIMENSION(1:nel), INTENT(out)   :: pdist, idist

    INTEGER                                  :: i

    DO i = 1 , nel
       pdist(i) = vdist(i) / nim
       idist(i) = MOD (vdist(i), nim) + 1
    END DO
  END SUBROUTINE v_to_p_i_dist
! *****************************************************************************
!> \brief ...
!> \param vdist ...
!> \param nim ...
!> \param pdist ...
!> \param idist ...
! *****************************************************************************
  SUBROUTINE v_to_p_i_dist_o (vdist, nim, pdist, idist)
    TYPE(array_i1d_obj), INTENT(in)          :: vdist
    INTEGER, INTENT(in)                      :: nim
    TYPE(array_i1d_obj), INTENT(out)         :: pdist, idist

    INTEGER                                  :: nel
    INTEGER, DIMENSION(:), POINTER           :: id, pd, vd

    nel = array_size (vdist)
    vd => array_data (vdist)
    ALLOCATE (pd(nel), id(nel))
    CALL v_to_p_i_dist (nel, vd, nim, pd, id)
    CALL array_new (pdist, pd, gift=.TRUE.)
    CALL array_new (idist, id, gift=.TRUE.)
  END SUBROUTINE v_to_p_i_dist_o


! *****************************************************************************
!> \brief Determines whether most blocks are stored transposed in normally.
!> \param[in] blk_p           Pointers to blocks
!> \retval mostly_non_transposed ...
!> \note Tries to be quick and not necessarily accurate.
! *****************************************************************************
  PURE FUNCTION mostly_non_transposed (blk_p)
    INTEGER, DIMENSION(:), INTENT(IN) :: blk_p
    LOGICAL                           :: mostly_non_transposed

    INTEGER            :: n, str, sntr
    INTEGER, PARAMETER :: allcheck_cutoff = 8

    n = SIZE (blk_p)
    IF (n .EQ. 0) THEN
       mostly_non_transposed = .TRUE.
       RETURN
    ENDIF
    str = 0
    sntr = 0
    CALL check_range (blk_p, 1, allcheck_cutoff, sntr, str)
    IF (n .GT. 4*allcheck_cutoff) THEN
       CALL check_range (blk_p, (n-allcheck_cutoff)/2, (n+allcheck_cutoff)/2,&
            sntr, str)
       CALL check_range (blk_p, n-allcheck_cutoff, n, sntr, str)
    ENDIF
    IF (str .EQ. 0) THEN
       mostly_non_transposed = .TRUE.
       RETURN
    ELSE
       ! Bias towards .TRUE.
       mostly_non_transposed = ((2*str) / (1+str+sntr)) .EQ. 0
    ENDIF
    RETURN
  CONTAINS
! *****************************************************************************
!> \brief ...
!> \param blk_p ...
!> \param lb ...
!> \param ub ...
!> \param sntr ...
!> \param str ...
! *****************************************************************************
    PURE SUBROUTINE check_range (blk_p, lb, ub, sntr, str)
    INTEGER, DIMENSION(:), INTENT(IN)        :: blk_p
    INTEGER, INTENT(IN)                      :: lb, ub
    INTEGER, INTENT(INOUT)                   :: sntr, str

    INTEGER                                  :: b1, b2

      b1 = MAX (1, lb)
      b2 = MIN (SIZE(blk_p), ub)
      sntr = sntr + COUNT (blk_p(b1:b2) .GT. 0)
      str = str + COUNT (blk_p(b1:b2) .LT. 0)
    END SUBROUTINE check_range
  END FUNCTION mostly_non_transposed

! *****************************************************************************
!> \brief Transform between virtual process rows/columns and actual process rows/columns and images therein.
!> \param image_dist ...
!> \param prow ...
!> \param rowi ...
!> \param pcol ...
!> \param coli ...
!> \param vprow ...
!> \param vpcol ...
!> \param myprow ...
!> \param mypcol ...
!> \param myrowi ...
!> \param mycoli ...
!> \param myvprow ...
!> \param myvpcol ...
!> \param vprow_shift ...
!> \param vpcol_shift ...
!> \param shifting ...
!> \param error ...
!> \par Shifting
!>  (L)eft and (R)ight shifting are "shifts from", (l)eft and (r)ight
!>  are "shifts to".  A caller (or the my* specifications) would use
!>  L/R to see which data he has (i.e., from where his data was
!>  shifted).  To see where the caller's data goes to, use l/r.
! *****************************************************************************
  SUBROUTINE image_calculator (image_dist,&
       prow, rowi, pcol, coli, vprow, vpcol,&
       myprow, mypcol, myrowi, mycoli, myvprow, myvpcol,&
       vprow_shift, vpcol_shift,&
       shifting, error)
    TYPE(dbcsr_imagedistribution_obj), &
      INTENT(IN)                             :: image_dist
    INTEGER, INTENT(OUT), OPTIONAL           :: prow, rowi, pcol, coli, &
                                                vprow, vpcol
    INTEGER, INTENT(IN), OPTIONAL            :: myprow, mypcol, myrowi, &
                                                mycoli, myvprow, myvpcol, &
                                                vprow_shift, vpcol_shift
    CHARACTER, INTENT(IN), OPTIONAL          :: shifting
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'image_calculator', &
      routineP = moduleN//':'//routineN

    INTEGER :: col_mult, my_pcol, my_prow, ncol_images, npcols, nprows, &
      nrow_images, nvpcols, nvprows, row_mult, vcol, vrow
    TYPE(dbcsr_mp_obj)                       :: mp

!   ---------------------------------------------------------------------------

    IF (careful_mod) THEN
       CALL dbcsr_assert (.not.PRESENT (myvprow), "IMP",&
            PRESENT(mycoli),&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Must specify either (process row and row image) or (virtual process row)",&
            __LINE__, error=error)
       CALL dbcsr_assert (.not.PRESENT (myvpcol), "IMP",&
            PRESENT(mycoli),&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Must specify either (process col and col image) or (virtual process col)",&
            __LINE__, error=error)
    ENDIF
    !
    mp = image_dist%i%main%d%mp_env
    nprows = SIZE (mp%mp%pgrid, 1)
    npcols = SIZE (mp%mp%pgrid, 2)
    nrow_images = image_dist%i%row_decimation
    ncol_images = image_dist%i%col_decimation
    row_mult = image_dist%i%row_multiplicity
    col_mult = image_dist%i%col_multiplicity
    nvprows = nprows * nrow_images
    nvpcols = npcols * ncol_images
    !
    IF (PRESENT (myprow)) THEN
       my_prow = myprow
    ELSE
       my_prow = mp%mp%myprow
    ENDIF
    IF (PRESENT (mypcol)) THEN
       my_pcol = mypcol
    ELSE
       my_pcol = mp%mp%mypcol
    ENDIF
    !
    IF (.NOT. PRESENT (myvprow)) THEN
       vrow = my_prow*nrow_images+myrowi-1
    ELSE
       vrow = myvprow
    ENDIF
    IF (.NOT.PRESENT (myvpcol)) THEN
       vcol = my_pcol*ncol_images+mycoli-1
    ELSE
       vcol = myvpcol
    ENDIF
    !
    IF (PRESENT (vprow_shift)) vrow = vrow + vprow_shift
    IF (PRESENT (vpcol_shift)) vcol = vcol + vpcol_shift
    IF (PRESENT (shifting)) THEN
       SELECT CASE (shifting)
       CASE ('R')
          vrow = vrow + my_pcol*row_mult
       CASE ('L')
          vcol = vcol + my_prow*col_mult
       CASE ('r')
          vrow = vrow - my_pcol*row_mult
       CASE ('l')
          vcol = vcol - my_prow*col_mult
       END SELECT
    ENDIF
    vrow = MODULO (vrow, nvprows)
    vcol = MODULO (vcol, nvpcols)
    IF (PRESENT (prow)) prow = vrow / nrow_images
    IF (PRESENT (rowi)) rowi = MODULO (vrow, nrow_images)+1
    IF (PRESENT (pcol)) pcol = vcol / ncol_images
    IF (PRESENT (coli)) coli = MODULO (vcol, ncol_images)+1
    IF (PRESENT (vprow)) vprow = vrow
    IF (PRESENT (vpcol)) vpcol = vcol
  END SUBROUTINE image_calculator

! *****************************************************************************
!> \brief Creates a sane 1-D distribution
!> \param[out] bin_dist    Distribution of elements to bins
!> \param[in] nelements    Number of elements
!> \param[in] nbins        Number of bins
!> \param[in] element_sizes  (optional) sizes of elements
!> \param[in,out] error    error
! *****************************************************************************
  SUBROUTINE dbcsr_dist_bin (bin_dist, nelements, nbins, element_sizes,&
       error)
    INTEGER, DIMENSION(:), INTENT(OUT), &
      POINTER                                :: bin_dist
    INTEGER, INTENT(IN)                      :: nelements, nbins
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: element_sizes
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_dist_bin', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: bin, bin_cnt, error_handle, i
    LOGICAL                                  :: found, heap_error
    REAL(kind=sp)                            :: rn
    TYPE(heap_t)                             :: bin_heap

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handle, error)
    ALLOCATE (bin_dist (nelements))
    IF (PRESENT (element_sizes)) THEN
       CALL dbcsr_assert (SIZE(element_sizes), "EQ", nelements,&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Array of element sizes does not match the number of elements.",&
            __LINE__, error=error)
       CALL heap_new (bin_heap, nbins, heap_error)
       CALL heap_fill (bin_heap, (/ (bin, bin = 0, nbins-1) /),&
                                 (/ (0,   bin = 0, nbins-1) /), heap_error)
       DO i = 1, nelements
          CALL heap_get_first (bin_heap, bin, bin_cnt, found, heap_error)
          bin_dist(i) = bin
          bin_cnt = bin_cnt + element_sizes(i)
          CALL heap_reset_first (bin_heap, bin_cnt, heap_error)
       ENDDO
       CALL heap_release (bin_heap, heap_error)
    ELSE
       DO i = 1, nelements
          CALL RANDOM_NUMBER(rn)
          bin_dist(i) = MOD (INT (rn*REAL (nbins, kind=sp)), nbins)
       ENDDO
    ENDIF
    CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE dbcsr_dist_bin

END MODULE dbcsr_dist_operations
