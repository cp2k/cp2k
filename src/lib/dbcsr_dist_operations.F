!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
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

  USE array_types,                     ONLY: array_data,&
                                             array_hold,&
                                             array_i1d_obj,&
                                             array_new,&
                                             array_nullify,&
                                             array_release,&
                                             array_size
  USE dbcsr_error_handling
  USE dbcsr_methods,                   ONLY: &
       dbcsr_distribution_col_dist, dbcsr_distribution_init, &
       dbcsr_distribution_local_cols, dbcsr_distribution_local_rows, &
       dbcsr_distribution_mp, dbcsr_distribution_ncols, &
       dbcsr_distribution_new, dbcsr_distribution_nrows, &
       dbcsr_distribution_processor, dbcsr_distribution_row_dist, &
       dbcsr_mp_new_transposed, dbcsr_mp_npcols, dbcsr_mp_nprows, &
       dbcsr_mp_release
  USE dbcsr_types,                     ONLY: dbcsr_distribution_obj,&
                                             dbcsr_imagedistribution_type,&
                                             dbcsr_magic_number,&
                                             dbcsr_mp_obj,&
                                             dbcsr_obj,&
                                             dbcsr_type
  USE dbcsr_util,                      ONLY: convert_sizes_to_offsets,&
                                             dbcsr_set_debug,&
                                             gcd,&
                                             lcm

  !$ USE OMP_LIB

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_dist_operations'

  REAL, PARAMETER                      :: default_resize_factor = 1.618034


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
            dbcsr_create_dist_r_rot,&
            dbcsr_create_dist_l_unrot,&
            dbcsr_create_dist_l_rot,&
            dbcsr_create_dist_elemental,&
            dbcsr_create_dist_fullest,&
            dbcsr_create_dist_block_cyclic

  PUBLIC :: dbcsr_create_image_dist,&
            make_image_distribution_dense

  ! Helper routines
  PUBLIC ::  evenly_bin,&
             rebin_distribution,&
             reimage_distribution,&
             rebin_imaged_distribution,&
             test_create_bl_distribution,&
             match_sizes_to_dist, make_sizes_dense,&
             find_local_virtual_elements,&
             mostly_non_transposed,&
             dbcsr_print_dist


  INTERFACE dbcsr_get_stored_block_info
     MODULE PROCEDURE get_stored_block_info_type, get_stored_block_info_obj
  END INTERFACE

  INTERFACE dbcsr_get_stored_coordinates
     MODULE PROCEDURE get_stored_coordinates_type, get_stored_coordinates_obj
  END INTERFACE

  INTERFACE dbcsr_get_block_index
     MODULE PROCEDURE get_block_index_type, get_block_index_obj
  END INTERFACE



#define DBG IF (dbg) WRITE(*,*)routineN//" ",
#define DBGV IF (bcsr_verbose) WRITE(*,*)routineN//" ",
#define DBGI IF (info) WRITE(*,*)routineN//" ",
#define DEBUG_HEADER        LOGICAL :: dbg, info
#define DEBUG_BODY        dbg = .FALSE. ; CALL dbcsr_set_debug(dbg, info=info)
#define DEBUG_BODYY       dbg = .TRUE. ; CALL dbcsr_set_debug(dbg, info=info)

  LOGICAL, PARAMETER :: bcsr_debug =   .TRUE.
  LOGICAL, PARAMETER :: bcsr_info =    .FALSE.
  LOGICAL, PARAMETER :: bcsr_verbose = .FALSE.

CONTAINS



! *****************************************************************************
!> \brief Determines whether a transpose must be applied
!> \param[in] row   The absolute matrix row.
!> \param[in] column          The absolute matrix column.
! *****************************************************************************
  ELEMENTAL FUNCTION checker_tr(row, column) RESULT(transpose)
    INTEGER, INTENT(IN)                      :: row, column
    LOGICAL                                  :: transpose

    transpose = BTEST(column+row, 0) .EQV. column.GT.row

  END FUNCTION checker_tr

! *****************************************************************************
!> \brief Determines the home process for a given logical matrix element.
!> \note This routine is a more low-level version of
!>       dbcsr_get_stored_coordinate without the setting the row and column
!>       to the stored position.
!> \note It assumes a square matrix.
!> \param[in] row             logical row
!> \param[in] col             logical column
!> \param[in] pgrid           process grid
!> \param[in] row_dist        row distribution
!> \param[in] col_dist        column distribution
!> \retval process            home process of the given element
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
       IF ((matrix%symmetry .AND. checker_tr(row, column)) .NEQV. transpose) THEN
          ! The transposed way.
          processor = dbcsr_distribution_processor (matrix%dist, column, row)
       ELSE
          ! The default way.
          processor = dbcsr_distribution_processor (matrix%dist, row, column)
       ENDIF
    ENDIF
  END SUBROUTINE get_stored_coordinates_type

! \brief Canonical logic
  PURE SUBROUTINE get_stored_canonical(matrix, row, column,&
       transpose, processor)
    TYPE(dbcsr_type), INTENT(IN)             :: matrix
    INTEGER, INTENT(INOUT)                   :: row, column
    LOGICAL, INTENT(INOUT)                   :: transpose
    INTEGER, INTENT(OUT), OPTIONAL           :: processor

    INTEGER                                  :: tmp
    LOGICAL                                  :: straight

! The old way

    straight = matrix%transpose .OR. matrix%symmetry
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

    CALL dbcsr_assert (matrix%m%initialized,'EQ',dbcsr_magic_number,&
         dbcsr_fatal_level, dbcsr_caller_error, routineN,&
         "matrix%m not associated",__LINE__,error)
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
!> \param[out] processor      (optional) returns the processor on which this
!>                            block resides
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
    IF (PRESENT (data_offset)) THEN
       CALL get_stored_block_info_type(matrix, stored_row, stored_column,&
            found, block_number, data_offset)
    ELSE
       CALL get_stored_block_info_type(matrix, stored_row, stored_column,&
            found, block_number)
    ENDIF
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
!> \param[out] processor      (optional) returns the processor on which this
!>                            block resides
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
               matrix%col_i(1:blk_last), matrix%blk_p(1:blk_last),&
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
    WRITE(*,'(1X,A,I5,"x",I5)')"Distribution process grid size",&
         dbcsr_mp_nprows (mp), dbcsr_mp_npcols (mp)
    IF (array_size(dist%d%row_dist).gt.0) &
         WRITE(*,*)"Distribution: Row dist", array_data(dist%d%row_dist)
    IF (array_size(dist%d%col_dist).gt.0) &
         WRITE(*,*)"Distribution: Col dist", array_data(dist%d%col_dist)
    IF (array_size(dist%d%local_rows).gt.0) &
         WRITE(*,*)"Distribution: Local rows", array_data(dist%d%local_rows)
    IF (array_size(dist%d%local_cols).gt.0) &
    WRITE(*,*)"Distribution: Local cols", array_data(dist%d%local_cols)
  END SUBROUTINE dbcsr_print_dist


! *****************************************************************************
!> \brief Transposes a distribution
!> \param[out] dist_tr        transposed distribution
!> \param[in] dist2d          current distribution
! *****************************************************************************
  SUBROUTINE dbcsr_transpose_distribution(dist_tr, dist_normal)
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT)                            :: dist_tr
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist_normal

    TYPE(array_i1d_obj)                      :: cd, lc, lr, rd
    TYPE(dbcsr_mp_obj)                       :: mp_env_tr

!   ---------------------------------------------------------------------------

    rd = dbcsr_distribution_col_dist (dist_normal)
    cd = dbcsr_distribution_row_dist (dist_normal)
    lr = dbcsr_distribution_local_cols (dist_normal)
    lc = dbcsr_distribution_local_rows (dist_normal)
    CALL dbcsr_mp_new_transposed (mp_env_tr, dbcsr_distribution_mp (&
         dist_normal))
    CALL dbcsr_distribution_new(dist_tr, mp_env_tr, rd, cd, lr, lc)
    CALL dbcsr_mp_release (mp_env_tr)
  END SUBROUTINE dbcsr_transpose_distribution

! *****************************************************************************
!> \brief Transposes a distribution but keeps the same mp_env
!> \param[out] dist_tr        transposed distribution
!> \param[in] dist2d          current distribution
! *****************************************************************************
  SUBROUTINE dbcsr_transpose_dims(dist_tr, dist_normal)
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT)                            :: dist_tr
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist_normal

    INTEGER, DIMENSION(:), POINTER           :: cda, rda
    TYPE(array_i1d_obj)                      :: cd, new_cd, new_rd, rd
    TYPE(dbcsr_mp_obj)                       :: mp

!   ---------------------------------------------------------------------------

    mp = dbcsr_distribution_mp (dist_normal)
    rd = dbcsr_distribution_row_dist (dist_normal)
    cd = dbcsr_distribution_col_dist (dist_normal)
    ALLOCATE (rda (array_size (cd)))
    CALL evenly_bin (rda, dbcsr_mp_nprows (mp))
    CALL array_new (new_rd, rda, gift=.TRUE.)

    ALLOCATE (cda (array_size (rd)))
    CALL evenly_bin (cda, dbcsr_mp_npcols (mp))
    CALL array_new (new_cd, cda, gift=.TRUE.)

    CALL dbcsr_distribution_new(dist_tr,&
         mp, new_rd, new_cd)
    CALL array_release (new_rd)
    CALL array_release (new_cd)
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
!> \param[out] new_bins      new real distribution
!> \param[out] images        new image distribution
!> \param[in] source_bins    Basis for the new distribution and images
!> \param[in] nbins          number of bins in the new real distribution
!> \param[in] multiplicity   multiplicity
!> \param[in] nimages        number of images in the new distribution
! *****************************************************************************
  SUBROUTINE rebin_distribution (new_bins, images, source_bins,&
       nbins, multiplicity, nimages)
    INTEGER, DIMENSION(:), INTENT(OUT)       :: new_bins, images
    INTEGER, DIMENSION(:), INTENT(IN)        :: source_bins
    INTEGER, INTENT(IN)                      :: nbins, multiplicity, nimages

    CHARACTER(len=*), PARAMETER :: routineN = 'rebin_distribution', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: debug_level = 2

    INTEGER                                  :: bin, i, old_nbins, virtual_bin
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: bin_multiplier
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------
!CALL dbcsr_assert (SIZE(source_bins), "EQ", SIZE(new_bins),&
!     dbcsr_warning_level, dbcsr_internal_error,&
!     routineN, "source/new bin array sizes not equal",__LINE__,error)
!CALL dbcsr_assert (SIZE(source_bins), "EQ", SIZE(images),&
!     dbcsr_warning_level, dbcsr_internal_error,&
!     routineN, "source/image array sizes not equal",__LINE__,error)

    CALL dbcsr_assert (MOD (nbins, multiplicity), "EQ", 0,&
         dbcsr_warning_level, dbcsr_internal_error,&
         routineN, "mulitplicity is not divisor of new process grid coordinate",__LINE__,error)
    old_nbins = nbins * nimages / multiplicity
    ALLOCATE (bin_multiplier(0:old_nbins-1))
    bin_multiplier(:) = 0
    DO i = 1, SIZE (new_bins)
       IF (i .LE. SIZE (source_bins)) THEN
          bin = source_bins (i)
       ELSE
          ! Fill remainder with a cyclic distribution
          bin = MOD (i, old_nbins)
       ENDIF
       bin_multiplier(bin) = bin_multiplier(bin)+1
       IF (bin_multiplier(bin) .GE. multiplicity) THEN
          bin_multiplier(bin) = 0
       ENDIF
       virtual_bin = bin*multiplicity+bin_multiplier(bin)
       new_bins(i) = virtual_bin / nimages
       images(i) = 1 + MOD (virtual_bin, nimages)
    ENDDO
  END SUBROUTINE rebin_distribution

! *****************************************************************************
!> \brief Makes new distribution with decimation and multiplicity
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
!> \param[out] images        new image distribution
!> \param[in] my_bins        Basis for the new images
!> \param[in] nbins          number of bins in the new real distribution
!> \param[in] multiplicity   multiplicity
!> \param[in] nimages        number of images in the new distribution
! *****************************************************************************
  SUBROUTINE reimage_distribution (images, my_bins,&
       nbins, multiplicity, nimages)
    INTEGER, DIMENSION(:), INTENT(OUT)       :: images
    INTEGER, DIMENSION(:), INTENT(IN)        :: my_bins
    INTEGER, INTENT(IN)                      :: nbins, multiplicity, nimages

    CHARACTER(len=*), PARAMETER :: routineN = 'reimage_distribution', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: debug_level = 2

    INTEGER                                  :: bin, i
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: bin_multiplier

!   ---------------------------------------------------------------------------
!CALL bcsr_assert (debug_level, SIZE(my_bins).EQ.SIZE(images),&
!     routineN, "array sizes not equal")

    ALLOCATE (bin_multiplier(0:nbins-1))
    bin_multiplier(:) = 0
    DO i = 1, SIZE (my_bins)
       bin = my_bins (i)
       bin_multiplier(bin) = bin_multiplier(bin)+1
       IF (bin_multiplier(bin) .GE. nimages) THEN
          bin_multiplier(bin) = 0
       ENDIF
       images(i) = 1 + bin_multiplier(bin)
    ENDDO
  END SUBROUTINE reimage_distribution


! *****************************************************************************
!> \brief Makes new distribution with multiplicity
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
!> \param[out] new_bins      new real distribution
!> \param[out] images        new image distribution
!> \param[in] source_bins    Basis for the new distribution and images
!> \param[in] source_images  Basis for the new distribution and images
!> \param[in] nbins          number of bins in the new real distribution
!> \param[in] multiplicity   multiplicity
!> \param[in] nimages        number of images in the new distribution
! *****************************************************************************
  SUBROUTINE rebin_imaged_distribution (new_bins, images,&
       source_bins, source_images, nbins, multiplicity, nimages)
    INTEGER, DIMENSION(:), INTENT(OUT)       :: new_bins, images
    INTEGER, DIMENSION(:), INTENT(IN)        :: source_bins, source_images
    INTEGER, INTENT(IN)                      :: nbins, multiplicity, nimages

    CHARACTER(len=*), PARAMETER :: routineN = 'rebin_imaged_distribution', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: debug_level = 2

    INTEGER                                  :: bin, i, virtual_bin
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------
!CALL dbcsr_assert (SIZE(source_bins), "EQ", SIZE(new_bins),&
!     dbcsr_warning_level, dbcsr_internal_error,&
!     routineN, "source/new array sizes not equal",__LINE__,error)
!CALL dbcsr_assert (SIZE(source_bins), "EQ", SIZE(images),&
!     dbcsr_warning_level, dbcsr_internal_error,&
!     routineN, "source/image array sizes not equal",__LINE__,error)

    CALL dbcsr_assert (MOD (nbins, multiplicity), "EQ", 0,&
         dbcsr_warning_level, dbcsr_internal_error,&
         routineN, "mulitplicity is not divisor of new process grid coordinate",__LINE__,error)
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
!> \brief Creates an image distribution with column images.
!>
!> The column images make the column distribution compatible with the
!> row distribution of the supplied other_dist distribution.
!>
!> This routine is used to make a left matrix (with new_row_dist)
!> column-comptable with a given right matrix (other_dist) for
!> multiplication with non-transposed processor grids.
!>
!> \param[out] imgdist        image distribution with columns
!> \param[in] new_row_dist    use this row distribution
!> \param[in] other_dist      makes the distribution compatible with and based
!>                            on this one
! *****************************************************************************
  SUBROUTINE dbcsr_create_col_image_dist(imgdist, new_row_dist, other_dist)
    TYPE(dbcsr_imagedistribution_type), &
      INTENT(OUT)                            :: imgdist
    TYPE(array_i1d_obj), INTENT(IN)          :: new_row_dist
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: other_dist

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_create_col_image_dist', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: debug_level, new_ncols, &
                                                npcols, nprows
    INTEGER, DIMENSION(:), POINTER           :: new_col_dist_data, &
                                                new_col_img_data, rd
    TYPE(array_i1d_obj)                      :: new_col_dist, new_col_image
    TYPE(dbcsr_distribution_obj)             :: new_dist
    TYPE(dbcsr_mp_obj)                       :: other_mp_env

    DEBUG_HEADER
!
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    IF (dbg) debug_level = 5
    CALL dbcsr_distribution_init (new_dist)
    CALL dbcsr_distribution_init (imgdist%main)
    CALL array_nullify (imgdist%col_image)
    CALL array_nullify (imgdist%row_image)
    other_mp_env = dbcsr_distribution_mp (other_dist)
    ! Determine the factors.
    new_ncols = dbcsr_distribution_nrows (other_dist) ! They must be compatible
    nprows = dbcsr_mp_nprows (other_mp_env)
    npcols = dbcsr_mp_npcols (other_mp_env)
    imgdist%col_decimation = lcm (nprows, npcols) / npcols
    imgdist%col_multiplicity = npcols / gcd (nprows, npcols)
    DBG 'image col decimation, multiplicity',&
         imgdist%col_decimation, imgdist%col_multiplicity
    imgdist%row_decimation = 1
    imgdist%row_multiplicity = 1
    ! Create the new column distribution and column image distribution
    ALLOCATE (new_col_dist_data (new_ncols))
    ALLOCATE (new_col_img_data (new_ncols))
    rd => array_data (dbcsr_distribution_row_dist (other_dist))
    CALL rebin_distribution (new_col_dist_data, new_col_img_data, rd,&
         npcols, imgdist%col_multiplicity, imgdist%col_decimation)
    CALL array_new (new_col_dist, new_col_dist_data, gift=.TRUE.)
    CALL array_new (new_col_image, new_col_img_data, gift=.TRUE.)
    ! Copy the row distribution from old distribution
    CALL dbcsr_distribution_new(new_dist,&
         other_mp_env,&
         new_row_dist, new_col_dist)
    CALL array_release (new_col_dist)
    ! Now finish the distribution image.
    imgdist%main = new_dist
    imgdist%col_image = new_col_image
  END SUBROUTINE dbcsr_create_col_image_dist


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
! *****************************************************************************
  SUBROUTINE dbcsr_create_image_dist(imgdist, dist,&
       match_row_pdist, match_row_idist, match_row_nbins,&
       match_col_pdist, match_col_idist, match_col_nbins)
    TYPE(dbcsr_imagedistribution_type), &
      INTENT(OUT)                            :: imgdist
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL, TARGET                       :: match_row_pdist, &
                                                match_row_idist
    INTEGER, INTENT(IN)                      :: match_row_nbins
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: match_col_pdist, &
                                                match_col_idist
    INTEGER, INTENT(IN)                      :: match_col_nbins

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_create_image_dist', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ncols, npcols, nprows, nrows
    INTEGER, DIMENSION(:), POINTER           :: col_dist_data, col_img_data, &
                                                row_dist_data, row_img_data
    LOGICAL                                  :: dbg
    TYPE(array_i1d_obj)                      :: col_dist, col_image, &
                                                row_dist, row_image
    TYPE(dbcsr_distribution_obj)             :: new_dist
    TYPE(dbcsr_mp_obj)                       :: mp_env

!   ---------------------------------------------------------------------------

    dbg = .FALSE.

    CALL dbcsr_distribution_init (new_dist)
    CALL dbcsr_distribution_init (imgdist%main)
    CALL array_nullify (imgdist%col_image)
    CALL array_nullify (imgdist%row_image)
    mp_env = dbcsr_distribution_mp (dist)
    ! Determine the factors.
    nrows = dbcsr_distribution_nrows (dist)
    ncols = dbcsr_distribution_ncols (dist)
    nprows = dbcsr_mp_nprows (mp_env)
    npcols = dbcsr_mp_npcols (mp_env)
    ! Create the new row distribution and row image distribution
    imgdist%row_decimation = lcm (nprows, match_row_nbins) / nprows
    imgdist%row_multiplicity = nprows / gcd (nprows, match_row_nbins)
    IF (dbg) WRITE (*,*)routineN//'row decimation, multiplicity',&
         imgdist%row_decimation, imgdist%row_multiplicity
    IF (dbg) WRITE(*,*)routineN//" nprows, match prows", nprows, match_row_nbins
    IF (match_row_nbins .EQ. nprows .AND. .NOT. PRESENT (match_row_pdist)) THEN
       row_dist = dbcsr_distribution_row_dist (dist)
       CALL array_hold (row_dist)
       ALLOCATE (row_img_data (nrows))
       row_img_data(:) = 1
       CALL array_new (row_image, row_img_data, gift=.TRUE.)
    ELSE
       ALLOCATE (row_img_data (nrows))
       IF (PRESENT (match_row_pdist)) THEN
          ALLOCATE (row_dist_data (nrows))
          IF (PRESENT (match_row_idist)) THEN
             CALL rebin_imaged_distribution (row_dist_data, row_img_data,&
                  match_row_pdist, match_row_idist,&
                  nprows,&
                  imgdist%row_multiplicity, imgdist%row_decimation)
          ELSE
             CALL rebin_distribution (row_dist_data, row_img_data,&
                  match_row_pdist, &
                  nprows,&
                  imgdist%row_multiplicity, imgdist%row_decimation)
          ENDIF
          CALL array_new (row_dist, row_dist_data, gift=.TRUE.)
       ELSE
          row_dist = dbcsr_distribution_row_dist (dist)
          CALL array_hold (row_dist)
          CALL reimage_distribution (row_img_data,&
               array_data (row_dist),&
               nprows,&
               imgdist%row_multiplicity, imgdist%row_decimation)
       ENDIF
       CALL array_new (row_image, row_img_data, gift=.TRUE.)
    ENDIF
    ! Create the new column distribution and column image distribution
    imgdist%col_decimation = lcm (npcols, match_col_nbins) / npcols
    imgdist%col_multiplicity = npcols / gcd (npcols, match_col_nbins)
    IF (dbg) WRITE (*,*)routineN//'col decimation, multiplicity',&
         imgdist%col_decimation, imgdist%col_multiplicity
    IF (dbg) WRITE(*,*)routineN//" npcols, match pcols", npcols, match_col_nbins
    IF (match_col_nbins .EQ. npcols .AND. .NOT. PRESENT (match_col_pdist)) THEN
       col_dist = dbcsr_distribution_col_dist (dist)
       CALL array_hold (col_dist)
       ALLOCATE (col_img_data (ncols))
       col_img_data(:) = 1
       CALL array_new (col_image, col_img_data, gift=.TRUE.)
    ELSE
       ALLOCATE (col_img_data (ncols))
       IF (PRESENT (match_col_pdist)) THEN
          ALLOCATE (col_dist_data (ncols))
          IF (PRESENT (match_col_idist)) THEN
             CALL rebin_imaged_distribution (col_dist_data, col_img_data,&
                  match_col_pdist, match_col_idist,&
                  npcols,&
                  imgdist%col_multiplicity, imgdist%col_decimation)
          ELSE
             CALL rebin_distribution (col_dist_data, col_img_data,&
                  match_col_pdist,&
                  npcols,&
                  imgdist%col_multiplicity, imgdist%col_decimation)
          ENDIF
          CALL array_new (col_dist, col_dist_data, gift=.TRUE.)
       ELSE
          col_dist = dbcsr_distribution_col_dist (dist)
          CALL array_hold (col_dist)
          CALL reimage_distribution (col_img_data,&
               array_data (col_dist),&
               npcols,&
               imgdist%col_multiplicity, imgdist%col_decimation)
       ENDIF
       CALL array_new (col_image, col_img_data, gift=.TRUE.)
    ENDIF
    ! Copy the row distribution from old distribution
    CALL dbcsr_distribution_new(new_dist,&
         mp_env,&
         row_dist, col_dist)
    CALL array_release (row_dist)
    CALL array_release (col_dist)
    ! Now finish the distribution image.
    imgdist%main = new_dist
    imgdist%col_image = col_image
    imgdist%row_image = row_image
  END SUBROUTINE dbcsr_create_image_dist



! *****************************************************************************
!> \brief Calculates the intersections of blocks
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
!> \note This routine is used in the counting and sending loops in
!>       dbcsr_complete_redistribute
!> \param[out] ints           intersections of blocks
!> \param[in] num_ints        maximum number of expected intersections
!> \param[out] n_src_dsts     offset and number intersections belonging
!>                            to source blocks
!> \param[in] src_sizes       sizes of source blocks
!> \param[in] dst_sizes       sizes of target blocks
!> \param[in,out] error       error
!> \par History
!>      UB [2009] Created
!>      UB [2010-03-22] Offsets of the intersections in the old/new blocks
!>                      are calculated & saved. Zero-length blocks are
!>                      skipped.
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
!> \brief Test routine for create_bl_distribution.
! *****************************************************************************
  SUBROUTINE test_create_bl_distribution()

    INTEGER, PARAMETER                       :: rpt = 10

    INTEGER                                  :: bins, elements, i
    REAL                                     :: rn
    TYPE(array_i1d_obj)                      :: block_distribution, block_size

    DO i = 1, rpt
       CALL RANDOM_NUMBER (rn)
       elements = INT(rn * 500)
       IF (i .EQ. 1) elements = 0
       CALL RANDOM_NUMBER (rn)
       bins = MAX(1, INT(rn*8))
       CALL create_bl_distribution (block_distribution,&
            block_size, elements, bins)
       CALL array_release (block_distribution)
       CALL array_release (block_size)
    ENDDO
  END SUBROUTINE test_create_bl_distribution

! *****************************************************************************
!> \brief Distributes elements into blocks and into bins
!>
!> \par Term clarification
!>      An example: blocks are atom blocks and bins are process rows/colums.
!> \param[out] block_distribution       block distribution to bins
!> \param[out] block_size       sizes of blocks
!> \param[in] elements          number of elements to bin
!> \param[in] nbins             number of bins
! *****************************************************************************
  SUBROUTINE create_bl_distribution (block_distribution,&
       block_size, nelements, nbins)
    TYPE(array_i1d_obj), INTENT(OUT)         :: block_distribution, block_size
    INTEGER, INTENT(IN)                      :: nelements, nbins

    CHARACTER(len=*), PARAMETER :: routineN = 'create_bl_distribution', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: max_elements_per_block = 16

    INTEGER :: bin, blk_layer, element_stack, els, estimated_blocks, &
      max_blocks_per_bin, nblks, nblocks, stat
    INTEGER, DIMENSION(:), POINTER           :: blk_dist, blk_sizes
    TYPE(dbcsr_error_type)                   :: error

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL array_nullify (block_distribution)
    CALL array_nullify (block_size)
    ! Define the sizes on which we build the distribution.
    IF (nelements .GT. 0) THEN

       nblocks = CEILING(REAL(nelements)/REAL(max_elements_per_block))
       max_blocks_per_bin = CEILING(REAL(nblocks)/REAL(nbins))

       IF (dbg) THEN
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
       CALL dbcsr_assert (stat == 0,&
            dbcsr_fatal_level, dbcsr_caller_error, routineN,&
            "blk_dist",__LINE__,error)
       ALLOCATE (blk_sizes (estimated_blocks), stat=stat)
       CALL dbcsr_assert (stat == 0,&
            dbcsr_fatal_level, dbcsr_caller_error, routineN,&
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
                IF (dbg) WRITE(*,'(1X,A,I5,A,I5,A,I5)')routineP//" Assigning",&
                     els," elements as block",nblks," to bin",bin
             ENDIF
          ENDDO
       ENDDO
       ! Create the output arrays.
       DBG "Wound up with",nblks,'blocks instead of',estimated_blocks
       IF (nblks .EQ. estimated_blocks) THEN
          CALL array_new (block_distribution, blk_dist, gift=.TRUE.)
          CALL array_new (block_size, blk_sizes, gift=.TRUE.)
       ELSE
          CALL array_new (block_distribution, blk_dist(1:nblks), lb=1)
          CALL array_new (block_size, blk_sizes(1:nblks), lb=1)
          DEALLOCATE (blk_dist)
          DEALLOCATE (blk_sizes)
       ENDIF
    ELSE
       ALLOCATE (blk_dist (0), stat=stat)
       CALL dbcsr_assert (stat == 0,&
            dbcsr_fatal_level, dbcsr_caller_error, routineN,&
            "blk_dist",__LINE__,error)
       ALLOCATE (blk_sizes (0), stat=stat)
       CALL dbcsr_assert (stat == 0,&
            dbcsr_fatal_level, dbcsr_caller_error, routineN,&
            "blk_sizes",__LINE__,error)
       CALL array_new (block_distribution, blk_dist, gift=.TRUE.)
       CALL array_new (block_size, blk_sizes, gift=.TRUE.)
    ENDIF
1579 FORMAT (I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5)
    IF (dbg) THEN
       WRITE(*,'(1X,A,A)')routineP//" Distribution"
       WRITE(*,1579)array_data (block_distribution)
       WRITE(*,'(1X,A,A)')routineP//" Sizes"
       WRITE(*,1579)array_data (block_size)
    ENDIF
  END SUBROUTINE create_bl_distribution

! *****************************************************************************
!> \brief Creates a new distribution for the right matrix in a matrix
!>        multiplication with unrotated grid.
!> \par The new row distribution for the right matrix is the same as the row
!>      distribution of the left matrix, while the column distribution is
!>      created so that it is appropriate to the parallel environment.
!> \param[out] dist_right     new distribution for the right matrix
!> \param[in] dist_left       the distribution of the left matrix
!> \param[in] ncolumns        number of columns in right matrix
!> \param[out] right_col_blk_sizes      sizes of blocks in the created column
! *****************************************************************************
  SUBROUTINE dbcsr_create_dist_r_unrot(dist_right, dist_left, ncolumns,&
       right_col_blk_sizes)
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT)                            :: dist_right
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist_left
    INTEGER, INTENT(IN)                      :: ncolumns
    TYPE(array_i1d_obj), INTENT(OUT)         :: right_col_blk_sizes

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_create_dist_r_unrot', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: multiplicity, new_npcols, &
                                                new_nprows, nimages
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: tmp_images
    INTEGER, DIMENSION(:), POINTER           :: right_rd
    TYPE(array_i1d_obj)                      :: left_col_dist, &
                                                right_col_dist, right_row_dist
    TYPE(dbcsr_mp_obj)                       :: right_mp_env

!   ---------------------------------------------------------------------------

    right_mp_env = dbcsr_distribution_mp (dist_left)
    ! Create the column distribution
    CALL create_bl_distribution (right_col_dist,&
         right_col_blk_sizes, ncolumns, dbcsr_mp_npcols (right_mp_env))
    ! Create an even row distribution.
    left_col_dist = dbcsr_distribution_col_dist (dist_left)
    CALL array_nullify (right_row_dist)
    ALLOCATE (right_rd(array_size (left_col_dist)))
    ALLOCATE (tmp_images(array_size (left_col_dist)))
    new_nprows = dbcsr_mp_nprows (right_mp_env)
    new_npcols = dbcsr_mp_npcols (right_mp_env)
    nimages = lcm (new_nprows, new_npcols) / new_nprows
    multiplicity = new_nprows / gcd (new_nprows, new_npcols)
    CALL rebin_distribution (right_rd, tmp_images, array_data (left_col_dist),&
         new_nprows, multiplicity, nimages)
    CALL array_new (right_row_dist, right_rd, gift=.TRUE.)
    CALL dbcsr_distribution_new (dist_right, right_mp_env,&
         right_row_dist, right_col_dist)
    CALL array_release (right_col_dist)
    CALL array_release (right_row_dist)
  END SUBROUTINE dbcsr_create_dist_r_unrot

! *****************************************************************************
!> \brief Creates a new distribution for the right matrix in a matrix
!>        multiplication with rotated grid.
!> \par The new row distribution for the right matrix is the same as the
!>      column distribution of the left matrix, while the column distribution
!>      is created so that it is appropriate to the parallel environment.
!> \param[out] dist_right     new distribution for the right matrix
!> \param[in] dist_left       the distribution of the left matrix
!> \param[in] ncolumns        number of columns in right matrix
!> \param[out] right_col_blk_sizes      sizes of blocks in the created column
! *****************************************************************************
  SUBROUTINE dbcsr_create_dist_r_rot(dist_right, dist_left, ncolumns,&
       right_col_blk_sizes)
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT)                            :: dist_right
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist_left
    INTEGER, INTENT(IN)                      :: ncolumns
    TYPE(array_i1d_obj), INTENT(OUT)         :: right_col_blk_sizes

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_create_dist_r_rot', &
      routineP = moduleN//':'//routineN

    TYPE(array_i1d_obj)                      :: right_col_dist
    TYPE(dbcsr_mp_obj)                       :: right_mp_env

!   ---------------------------------------------------------------------------
!CALL test_create_bl_distribution (error)

    CALL dbcsr_mp_new_transposed (right_mp_env,&
         dbcsr_distribution_mp (dist_left))
    CALL create_bl_distribution (right_col_dist,&
         right_col_blk_sizes, ncolumns, dbcsr_mp_npcols (right_mp_env))
    CALL dbcsr_distribution_new (dist_right, right_mp_env,&
         dbcsr_distribution_col_dist (dist_left), right_col_dist,&
         local_rows=dbcsr_distribution_local_cols (dist_left))
    CALL dbcsr_mp_release (right_mp_env)
    CALL array_release (right_col_dist)
  END SUBROUTINE dbcsr_create_dist_r_rot

! *****************************************************************************
!> \brief Creates a new distribution for the left matrix in a matrix
!>        multiplication with unrotated grid.
!> \par The new column distribution for the left matrix is the same as the
!>      column distribution of the right matrix, while the row distribution is
!>      created so that it is appropriate to the parallel environment.
!> \param[out] dist_left      new distribution for the left matrix
!> \param[in] dist_right      old distribution of the right matrix
!> \param[in] nrows           number of rows in the left matrix
!> \param[out] left_row_blk_sizes       sizes of blocks in the created rows
! *****************************************************************************
  SUBROUTINE dbcsr_create_dist_l_unrot(dist_left, dist_right, nrows,&
       left_row_blk_sizes)
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT)                            :: dist_left
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist_right
    INTEGER, INTENT(IN)                      :: nrows
    TYPE(array_i1d_obj), INTENT(OUT)         :: left_row_blk_sizes

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_create_dist_l_unrot', &
      routineP = moduleN//':'//routineN

    TYPE(array_i1d_obj)                      :: left_row_dist
    TYPE(dbcsr_mp_obj)                       :: left_mp_env

!   ---------------------------------------------------------------------------

    left_mp_env = dbcsr_distribution_mp (dist_left)
    CALL create_bl_distribution (left_row_dist,&
         left_row_blk_sizes, nrows, dbcsr_mp_nprows (left_mp_env))
    CALL dbcsr_distribution_new (dist_left, left_mp_env,&
         left_row_dist, dbcsr_distribution_col_dist (dist_right),&
         local_cols=dbcsr_distribution_local_cols (dist_right))
    CALL array_release (left_row_dist)
  END SUBROUTINE dbcsr_create_dist_l_unrot

! *****************************************************************************
!> \brief Creates a new distribution for the left matrix in a matrix
!>        multiplication with rotated grid.
!> \par The new column distribution for the left matrix is the same as the
!>      row distribution of the right matrix, while the row distribution is
!>      created so that it is appropriate to the parallel environment.
!> \param[out] dist_left      new distribution for the left matrix
!> \param[in] dist_right      the distribution of the right matrix
!> \param[in] nrows           (optional) number of rows in left matrix
!> \param[out] left_row_blk_sizes      sizes of blocks in the created rows
! *****************************************************************************
  SUBROUTINE dbcsr_create_dist_l_rot(dist_right, dist_left, nrows,&
       left_row_blk_sizes)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist_right
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT)                            :: dist_left
    INTEGER, INTENT(IN)                      :: nrows
    TYPE(array_i1d_obj), INTENT(OUT)         :: left_row_blk_sizes

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_create_dist_l_rot', &
      routineP = moduleN//':'//routineN

    TYPE(array_i1d_obj)                      :: left_row_dist
    TYPE(dbcsr_mp_obj)                       :: left_mp_env

!   ---------------------------------------------------------------------------

    CALL dbcsr_mp_new_transposed (left_mp_env,&
         dbcsr_distribution_mp (dist_right))
    CALL create_bl_distribution (left_row_dist,&
         left_row_blk_sizes, nrows, dbcsr_mp_nprows (left_mp_env))
    CALL dbcsr_distribution_new (dist_left, left_mp_env,&
         left_row_dist, dbcsr_distribution_row_dist (dist_right),&
         local_cols=dbcsr_distribution_local_rows (dist_right))
    CALL dbcsr_mp_release (left_mp_env)
    CALL array_release (left_row_dist)
  END SUBROUTINE dbcsr_create_dist_l_rot

! *****************************************************************************
!> \brief Creates a distribution where all block sizes are 1.
!> \param[out] dist           the elemental distribution
!> \param[in] nrows           number of full rows
!> \param[in] ncolumns        number of full columns
!> \param[in] mp_env          mp environment
!> \param[out] row_blk_sizes  (optional) row block sizes (all 1)
!> \param[out] col_blk_sizes  (optional) column block sizes (all 1)
! *****************************************************************************
  SUBROUTINE dbcsr_create_dist_elemental(dist, nrows, ncolumns,&
       mp_env, row_blk_sizes, col_blk_sizes)
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT)                            :: dist
    INTEGER, INTENT(IN)                      :: nrows, ncolumns
    TYPE(dbcsr_mp_obj)                       :: mp_env
    TYPE(array_i1d_obj), INTENT(OUT), &
      OPTIONAL                               :: row_blk_sizes, col_blk_sizes

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_create_dist_elemental', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: npdims, pdim, sz
    INTEGER, DIMENSION(:), POINTER           :: cbs_data, cd_data, rbs_data, &
                                                rd_data
    TYPE(array_i1d_obj)                      :: cbs, col_dist, rbs, row_dist

!   ---------------------------------------------------------------------------

    ALLOCATE (rbs_data(nrows), cbs_data(ncolumns))
    rbs_data(:) = 1
    cbs_data(:) = 1
    CALL array_new (rbs, rbs_data, gift=.TRUE.)
    CALL array_new (cbs, cbs_data, gift=.TRUE.)
    ALLOCATE (rd_data(nrows), cd_data(ncolumns))
    ! Calclate process row distribution
    npdims = dbcsr_mp_nprows (mp_env)
    sz = CEILING (REAL(nrows) / REAL(npdims))
    FORALL (pdim = 0 : npdims-1)
       rd_data(pdim*sz+1 : MIN((pdim+1)*sz, nrows)) = pdim
    END FORALL
    CALL array_new (row_dist, rd_data, gift=.TRUE.)
    ! Calclate process column distribution
    npdims = dbcsr_mp_npcols (mp_env)
    sz = CEILING (REAL(ncolumns) / REAL(npdims))
    FORALL (pdim = 0 : npdims-1)
       cd_data(pdim*sz+1 : MIN((pdim+1)*sz, ncolumns)) = pdim
    END FORALL
    CALL array_new (col_dist, cd_data, gift=.TRUE.)
    !
    CALL dbcsr_distribution_new (dist, mp_env,&
         row_dist, col_dist)
    CALL array_release (row_dist)
    CALL array_release (col_dist)
    IF (PRESENT (row_blk_sizes)) THEN
       row_blk_sizes = rbs
    ELSE
       CALL array_release (rbs)
    ENDIF
    IF (PRESENT (col_blk_sizes)) THEN
       col_blk_sizes = cbs
    ELSE
       CALL array_release (cbs)
    ENDIF
  END SUBROUTINE dbcsr_create_dist_elemental


! *****************************************************************************
!> \brief Creates a distribution where all block sizes are as large as
!>        possible for the given distribution.
!> \param[out] dist           the elemental distribution
!> \param[in] nrows           number of full rows
!> \param[in] ncolumns        number of full columns
!> \param[in] mp_env          mp environment
!> \param[out] row_blk_sizes  row block sizes
!> \param[out] col_blk_sizes  column block sizes
! *****************************************************************************
  SUBROUTINE dbcsr_create_dist_fullest(dist, nrows, ncolumns,&
       mp_env, row_blk_sizes, col_blk_sizes)
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT)                            :: dist
    INTEGER, INTENT(IN)                      :: nrows, ncolumns
    TYPE(dbcsr_mp_obj)                       :: mp_env
    TYPE(array_i1d_obj), INTENT(OUT)         :: row_blk_sizes, col_blk_sizes

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_create_dist_fullest', &
      routineP = moduleN//':'//routineN

    TYPE(array_i1d_obj)                      :: col_dist, row_dist

!   ---------------------------------------------------------------------------

    CALL create_bl_distribution(row_dist, row_blk_sizes,&
         nrows, dbcsr_mp_nprows (mp_env))
    CALL create_bl_distribution(col_dist, col_blk_sizes,&
         ncolumns, dbcsr_mp_npcols (mp_env))
    !
    CALL dbcsr_distribution_new (dist, mp_env,&
         row_dist, col_dist)
    !
    CALL array_release (row_dist)
    CALL array_release (col_dist)
  END SUBROUTINE dbcsr_create_dist_fullest

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
    TYPE(array_i1d_obj), INTENT(OUT)         :: row_blk_sizes, col_blk_sizes

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'dbcsr_create_dist_block_cyclic', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nblkcols, nblkrows, npcols, &
                                                nprows, pdim, sz
    INTEGER, DIMENSION(:), POINTER           :: cbs_data, cd_data, rbs_data, &
                                                rd_data
    TYPE(array_i1d_obj)                      :: col_dist, row_dist

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY

    ! Row sizes
    IF(nrow_block.EQ.0) THEN
       nblkrows = 0
       sz = 0
    ELSE
       nblkrows = nrows / nrow_block
       sz = MOD(nrows, nrow_block)
    ENDIF
    IF (sz .GT. 0) nblkrows = nblkrows + 1
    ALLOCATE (rbs_data(nblkrows), rd_data(nblkrows))
    rbs_data(1:nblkrows) = nrow_block
    IF (sz .NE. 0) rbs_data(nblkrows) = sz

    ! Column sizes
    IF(ncol_block.EQ.0) THEN
       nblkcols = 0
       sz = 0
    ELSE
       nblkcols = ncolumns / ncol_block
       sz = MOD(ncolumns, ncol_block)
    ENDIF
    IF (sz .GT. 0) nblkcols = nblkcols + 1
    ALLOCATE (cbs_data(nblkcols), cd_data(nblkcols))
    cbs_data(1:nblkcols) = ncol_block
    IF (sz .NE. 0) cbs_data(nblkcols) = sz
    CALL array_nullify (row_blk_sizes)
    CALL array_nullify (col_blk_sizes)
    CALL array_new (row_blk_sizes, rbs_data, gift=.TRUE.)
    CALL array_new (col_blk_sizes, cbs_data, gift=.TRUE.)
    !
    IF(dbg) THEN
       WRITE(*,*)routineN//" nrows,nrow_block,nblkrows=",&
            nrows,nrow_block,nblkrows
       WRITE(*,*)routineN//" ncols,ncol_block,nblkcols=",&
            ncolumns,ncol_block,nblkcols
    ENDIF
    ! Calculate process row distribution
    nprows = dbcsr_mp_nprows (mp_env)
    FORALL (pdim = 0 : MIN(nprows-1,nblkrows-1))
       rd_data(1+pdim:nblkrows:nprows) = pdim
    END FORALL
    CALL array_new (row_dist, rd_data, gift=.TRUE.)
    ! Calculate process column distribution
    npcols = dbcsr_mp_npcols (mp_env)
    FORALL (pdim = 0 : MIN(npcols-1, nblkcols-1))
       cd_data(1+pdim:nblkcols:npcols) = pdim
    END FORALL
    CALL array_new (col_dist, cd_data, gift=.TRUE.)
    !
    IF(dbg) THEN
       WRITE(*,*)routineN//" row_dist",&
            rd_data
       WRITE(*,*)routineN//" col_dist",&
            cd_data
    ENDIF
    !
    CALL dbcsr_distribution_new (dist, mp_env,&
         row_dist, col_dist)
    !
    CALL array_release (row_dist)
    CALL array_release (col_dist)
  END SUBROUTINE dbcsr_create_dist_block_cyclic


! *****************************************************************************
!> \brief Make the image distribution reflect the dense matrix blocks.
!> \note Used for making matrices dense/undense
! *****************************************************************************
  SUBROUTINE make_image_distribution_dense (image_dist, join_rows, join_cols)
    TYPE(dbcsr_imagedistribution_type), &
      INTENT(INOUT)                          :: image_dist
    LOGICAL, INTENT(IN), OPTIONAL            :: join_rows, join_cols

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'make_image_distribution_dense', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ncol_images, npcols, &
                                                nprows, nrow_images, nvpcols, &
                                                nvprows
    INTEGER, DIMENSION(:), POINTER           :: new_cid, new_rid
    LOGICAL                                  :: jcols, jrows

!   ---------------------------------------------------------------------------

    jcols = .TRUE.
    IF (PRESENT (join_cols)) jcols = join_cols
    jrows = .TRUE.
    IF (PRESENT (join_rows)) jrows = join_rows
    nrow_images = image_dist%row_decimation
    ncol_images = image_dist%col_decimation
    nprows = dbcsr_mp_nprows (dbcsr_distribution_mp (image_dist%main))
    npcols = dbcsr_mp_npcols (dbcsr_distribution_mp (image_dist%main))
    nvprows = nprows * nrow_images
    nvpcols = npcols * ncol_images
    IF (jrows) THEN
       CALL array_release (image_dist%row_image)
       CALL array_nullify (image_dist%row_image)
       ALLOCATE (new_rid (nvprows))
       FORALL (i = 1 : nvprows)
          new_rid(i) = MOD (i-1, nrow_images) + 1
       END FORALL
       CALL array_new (image_dist%row_image, new_rid, gift=.TRUE.)
    ENDIF
    IF (jcols) THEN
       CALL array_release (image_dist%col_image)
       CALL array_nullify (image_dist%col_image)
       ALLOCATE (new_cid (nvpcols))
       FORALL (i = 1 : nvpcols)
          new_cid(i) = MOD (i-1, ncol_images) + 1
       END FORALL
       CALL array_new (image_dist%col_image, new_cid, gift=.TRUE.)
    ENDIF
  END SUBROUTINE make_image_distribution_dense


! *****************************************************************************
!> \brief Matches row/block sizes and offsets to a given distribution, inplace.
!> \note Used for making matrices dense/undense
! *****************************************************************************
  SUBROUTINE match_sizes_to_dist (new_distribution, old_distribution,&
       row_blk_sizes, col_blk_sizes, row_blk_offsets, col_blk_offsets,&
       old_image_dist, old_col_idist, old_col_nimages)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: new_distribution, &
                                                old_distribution
    TYPE(array_i1d_obj), INTENT(INOUT)       :: row_blk_sizes, col_blk_sizes, &
                                                row_blk_offsets, &
                                                col_blk_offsets
    TYPE(dbcsr_imagedistribution_type), &
      INTENT(IN), OPTIONAL                   :: old_image_dist, &
                                                old_col_idist, old_col_nimages

    CHARACTER(len=*), PARAMETER :: routineN = 'match_sizes_to_dist', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nvpcols, nvprows, old_npcols, &
                                                old_nprows
    INTEGER, DIMENSION(:), POINTER :: cbo, cbs, new_pcol_dist, new_prow_dist, &
      old_cbo, old_icol_dist, old_irow_dist, old_pcol_dist, old_prow_dist, &
      old_rbo, rbo, rbs
    LOGICAL                                  :: images
    TYPE(array_i1d_obj)                      :: old_col_blk_offsets, &
                                                old_col_blk_sizes, &
                                                old_row_blk_offsets, &
                                                old_row_blk_sizes
    TYPE(dbcsr_mp_obj)                       :: mp, old_mp

!   ---------------------------------------------------------------------------

    images = PRESENT (old_image_dist)
    mp = dbcsr_distribution_mp (new_distribution)
    old_mp = dbcsr_distribution_mp (old_distribution)
    ! Save old data
    old_prow_dist => array_data (dbcsr_distribution_row_dist (old_distribution))
    old_pcol_dist => array_data (dbcsr_distribution_col_dist (old_distribution))
    old_row_blk_sizes = row_blk_sizes
    old_col_blk_sizes = col_blk_sizes
    old_row_blk_offsets = row_blk_offsets
    old_col_blk_offsets = col_blk_offsets
    old_rbo => array_data (row_blk_offsets)
    old_cbo => array_data (col_blk_offsets)
    old_nprows = dbcsr_mp_nprows (old_mp)
    old_npcols = dbcsr_mp_npcols (old_mp)
    CALL array_nullify (row_blk_sizes)
    CALL array_nullify (col_blk_sizes)
    CALL array_nullify (row_blk_offsets)
    CALL array_nullify (col_blk_offsets)
    new_prow_dist => array_data (dbcsr_distribution_row_dist (new_distribution))
    new_pcol_dist => array_data (dbcsr_distribution_col_dist (new_distribution))
    IF (images) THEN
       nvprows = dbcsr_mp_nprows(mp)*old_image_dist%row_decimation
       nvpcols = dbcsr_mp_npcols(mp)*old_image_dist%col_decimation
       ALLOCATE (rbs(nvprows), cbs(nvpcols), rbo(nvprows+1), cbo(nvpcols+1))
       old_irow_dist => array_data (old_image_dist%row_image)
       old_icol_dist => array_data (old_image_dist%col_image)
       CALL match_sizes_to_el_dist(rbs, old_rbo, new_prow_dist,&
            old_prow_dist, old_irow_dist, old_image_dist%row_decimation)
       CALL match_sizes_to_el_dist(cbs, old_cbo, new_pcol_dist,&
            old_pcol_dist, old_icol_dist, old_image_dist%col_decimation)
    ELSE
       nvprows = dbcsr_mp_nprows (mp)
       nvpcols = dbcsr_mp_npcols (mp)
       ALLOCATE (rbs(nvprows), cbs(nvpcols), rbo(nvprows+1), cbo(nvpcols+1))
       CALL match_sizes_to_el_dist (rbs, old_rbo, new_prow_dist,&
            old_prow_dist)
       CALL match_sizes_to_el_dist (cbs, old_cbo, new_pcol_dist,&
            old_pcol_dist)
    ENDIF
    CALL convert_sizes_to_offsets (rbs, rbo)
    CALL convert_sizes_to_offsets (cbs, cbo)
    CALL array_new (row_blk_sizes, rbs, gift=.TRUE.)
    CALL array_new (col_blk_sizes, cbs, gift=.TRUE.)
    CALL array_new (row_blk_offsets, rbo, gift=.TRUE.)
    CALL array_new (col_blk_offsets, cbo, gift=.TRUE.)
    CALL array_release (old_row_blk_sizes)
    CALL array_release (old_col_blk_sizes)
    CALL array_release (old_row_blk_offsets)
    CALL array_release (old_col_blk_offsets)
  END SUBROUTINE match_sizes_to_dist

! *****************************************************************************
!> \brief Matches row/block sizes and offsets to a given distribution, inplace.
!> \note Used for making matrices dense/undense
! *****************************************************************************
  SUBROUTINE make_sizes_dense (new_distribution, old_distribution,&
       blk_sizes, blk_offsets, nbins,&
       old_image_dist, nimages)
    TYPE(array_i1d_obj), INTENT(IN)          :: new_distribution, &
                                                old_distribution
    TYPE(array_i1d_obj), INTENT(INOUT)       :: blk_sizes, blk_offsets
    INTEGER, INTENT(IN)                      :: nbins
    INTEGER, DIMENSION(:), INTENT(in), &
      OPTIONAL                               :: old_image_dist
    INTEGER, INTENT(IN), OPTIONAL            :: nimages

    CHARACTER(len=*), PARAMETER :: routineN = 'make_sizes_dense', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nvbins
    INTEGER, DIMENSION(:), POINTER           :: bo, bs, new_bin_dist, &
                                                old_bin_dist, old_bo
    LOGICAL                                  :: images
    TYPE(array_i1d_obj)                      :: old_blk_offsets, old_blk_sizes

!   ---------------------------------------------------------------------------

    images = PRESENT (old_image_dist) .AND. PRESENT (nimages)
    ! Save old data
    old_bin_dist => array_data (old_distribution)
    old_blk_sizes = blk_sizes
    old_blk_offsets = blk_offsets
    old_bo => array_data (blk_offsets)
    CALL array_nullify (blk_sizes)
    CALL array_nullify (blk_offsets)
    new_bin_dist => array_data (new_distribution)
    IF (images) THEN
       nvbins = nbins * nimages
       ALLOCATE (bs(nvbins), bo(nvbins+1))
       CALL match_sizes_to_el_dist(bs, old_bo, new_bin_dist,&
            old_bin_dist, old_image_dist, nimages)
    ELSE
       nvbins = nbins
       ALLOCATE (bs(nvbins), bo(nvbins+1))
       CALL match_sizes_to_el_dist (bs, old_bo, new_bin_dist,&
            old_bin_dist)
    ENDIF
    CALL convert_sizes_to_offsets (bs, bo)
    CALL array_new (blk_sizes, bs, gift=.TRUE.)
    CALL array_new (blk_offsets, bo, gift=.TRUE.)
    CALL array_release (old_blk_sizes)
    CALL array_release (old_blk_offsets)
  END SUBROUTINE make_sizes_dense


! *****************************************************************************
!> \brief Determines block sizes from old offsets and old and new
!>        distributions.
!> \note Used for making matrices dense/undense
! *****************************************************************************
  PURE SUBROUTINE match_sizes_to_el_dist(new_sizes, old_offsets,&
       new_dist, old_dist, old_idist, old_nimages)
    INTEGER, DIMENSION(:), INTENT(OUT)       :: new_sizes
    INTEGER, DIMENSION(:), INTENT(IN)        :: old_offsets, new_dist, &
                                                old_dist
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: old_idist
    INTEGER, INTENT(IN), OPTIONAL            :: old_nimages

    INTEGER                                  :: new_bin, old_bin, old_el, &
                                                old_size
    LOGICAL                                  :: images

    images = PRESENT (old_idist) .AND. PRESENT (old_nimages)
    new_sizes(:) = 0
    DO old_el = 1, SIZE(old_dist)
       old_bin = old_dist(old_el)
       old_size = old_offsets(old_el+1)-old_offsets(old_el)
       IF (images) THEN
          new_bin = old_bin*old_nimages + old_idist(old_el)
       ELSE
          new_bin = 1 + old_bin
       ENDIF
       new_sizes(new_bin) = new_sizes(new_bin) + old_size
    ENDDO
  END SUBROUTINE match_sizes_to_el_dist


! *****************************************************************************
!> \brief Finds the block that has the given column.
!>
!>        If the block having the queried column is found, the blk parameter
!>        is set to this block number and the found parameter is true.
!>        Otherwise found is false and the block number is invalid.
!> \par Index validity
!>      The blk_p array of block pointers is a required parameter to enable
!>      the detection of deleted blocks.
!> \param[in] find_col        column to find
!> \param[in] frst_blk        first block number in row
!> \param[in] last_blk        last block number in row
!> \param[in] col_i           col indices
!> \param[in] blk_p           block pointers
!> \param[out] blk            block number with searched-for column
!> \param[out] found          flag specified whether a block that has the
!>                            correct column was found
! *****************************************************************************
  PURE SUBROUTINE dbcsr_find_column(find_col,frst_blk,last_blk,col_i,blk_p,&
       blk,found)
    INTEGER, INTENT(IN)                      :: find_col, frst_blk, last_blk
    INTEGER, DIMENSION(last_blk), INTENT(IN) :: col_i, blk_p
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
!> \brief Determines local virtual elements
!>        Analogous to the local_rows
!> \param[out] local_virtual_elements    enumerated local virtual elements
!> \param[out] nvle                      number of local virtual elements
!> \param[in] grid_distribution          distribution of elements to bins
!> \param[in] image_distribution         distribution into images of all
!>                                       elements
!> \param[in] local_grid                 local grid
!> \param[in] local_image                local image
!> \note virtual_local_elements should be large enough to accomodate
!>       all possible data
!> \note Used for making matrices dense/undense
! *****************************************************************************
  PURE SUBROUTINE find_local_virtual_elements (local_virtual_elements, nlve,&
       grid_distribution, image_distribution, local_grid, local_image)
    INTEGER, DIMENSION(:), INTENT(OUT)       :: local_virtual_elements
    INTEGER, INTENT(OUT)                     :: nlve
    INTEGER, DIMENSION(:), INTENT(IN)        :: grid_distribution
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: image_distribution
    INTEGER, INTENT(IN)                      :: local_grid
    INTEGER, INTENT(IN), OPTIONAL            :: local_image

    INTEGER                                  :: el

    nlve = 0
    IF (PRESENT (image_distribution) .AND. PRESENT (local_image)) THEN
       DO el = 1, SIZE (grid_distribution)
          IF (grid_distribution(el) .EQ. local_grid&
               .AND. image_distribution(el) .EQ. local_image) THEN
             nlve = nlve + 1
             local_virtual_elements(nlve) = el
          ENDIF
       ENDDO
    ELSE
       DO el = 1, SIZE (grid_distribution)
          IF (grid_distribution(el) .EQ. local_grid) THEN
             nlve = nlve + 1
             local_virtual_elements(nlve) = el
          ENDIF
       ENDDO
    ENDIF
  END SUBROUTINE find_local_virtual_elements

! *****************************************************************************
!> \brief Determines whether most blocks are stored transposed in normally.
!> \note Tries to be quick and not necessarily accurate.
!> \param[in] blk_p           Pointers to blocks
! *****************************************************************************
  PURE FUNCTION mostly_non_transposed (blk_p)
    INTEGER, DIMENSION(:), INTENT(IN) :: blk_p
    LOGICAL                           :: mostly_non_transposed

    INTEGER            :: n, lb, ub, str, sntr
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
    PURE SUBROUTINE check_range (blk_p, lb, ub, sntr, str)
    INTEGER, DIMENSION(:), INTENT(IN)        :: blk_p
    INTEGER, INTENT(IN)                      :: lb, ub
    INTEGER, INTENT(INOUT)                   :: sntr, str

    INTEGER                                  :: b1, b2

      b1 = MAX (1, lb)
      b2 = MIN (SIZE(blk_p), b2)
      sntr = sntr + COUNT (blk_p(b1:b2) .GT. 0)
      str = str + COUNT (blk_p(b1:b2) .LT. 0)
    END SUBROUTINE check_range
  END FUNCTION mostly_non_transposed


END MODULE dbcsr_dist_operations
