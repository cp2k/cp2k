!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   DBCSR input/output
!> \author  Urban Borstnik
!> \date    2010-02-18
!> \version 0.9
!>
!> <b>Modification history:</b>
!> - 2010-02-18 Moved from dbcsr_util
! *****************************************************************************
MODULE dbcsr_io

  USE array_types,                     ONLY: array_data,&
                                             array_i1d_obj,&
                                             array_new,&
                                             array_release
  USE dbcsr_data_methods,              ONLY: dbcsr_data_clear_pointer,&
                                             dbcsr_data_init,&
                                             dbcsr_data_new,&
                                             dbcsr_get_data
  USE dbcsr_error_handling
  USE dbcsr_methods,                   ONLY: &
       dbcsr_data_area, dbcsr_distribution, dbcsr_distribution_mp, &
       dbcsr_get_data_size, dbcsr_get_data_size_referenced, &
       dbcsr_get_data_type, dbcsr_get_matrix_type, dbcsr_get_num_blocks, &
       dbcsr_mp_group, dbcsr_mp_numnodes, dbcsr_mp_pgrid, dbcsr_name, &
       dbcsr_nblkcols_total, dbcsr_nblkrows_total, dbcsr_valid_index
  USE dbcsr_mp_operations,             ONLY: dbcsr_mp_type_from_anytype
  USE dbcsr_transformations,           ONLY: dbcsr_datablock_redistribute
  USE dbcsr_types,                     ONLY: dbcsr_data_obj,&
                                             dbcsr_distribution_obj,&
                                             dbcsr_mp_obj,&
                                             dbcsr_obj,&
                                             dbcsr_type_complex_4,&
                                             dbcsr_type_complex_8,&
                                             dbcsr_type_real_4,&
                                             dbcsr_type_real_8
  USE dbcsr_work_operations,           ONLY: dbcsr_create
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             real_4,&
                                             real_4_size,&
                                             real_8,&
                                             real_8_size,&
                                             sp
  USE message_passing,                 ONLY: &
       file_amode_create, file_amode_rdonly, file_amode_wronly, file_offset, &
       mp_allgather, mp_environ, mp_file_close, mp_file_get_size, &
       mp_file_open, mp_file_read_at_all, mp_file_write_at, &
       mp_file_write_at_all, mp_sum, mp_type_descriptor_type, mp_type_size, &
       mpi_character_size, mpi_integer_size

  !$ USE OMP_LIB

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_io'


  ! Main
  PUBLIC :: dbcsr_print
  ! Low-level printing
! PUBLIC :: dbcsr_printmat, dbcsr_print2dmat
  ! Utility printing
  PUBLIC :: all_print
  PUBLIC :: dbcsr_binary_write
  PUBLIC :: dbcsr_binary_read
  PUBLIC :: dbcsr_print_block_sum

  LOGICAL, PARAMETER :: bcsr_debug =   .TRUE.
  LOGICAL, PARAMETER :: bcsr_info =    .FALSE.
  LOGICAL, PARAMETER :: bcsr_verbose = .FALSE.

  INTERFACE dbcsr_printmat
     MODULE PROCEDURE printmat_s, printmat_d, printmat_c, printmat_z
  END INTERFACE

CONTAINS


! *****************************************************************************
!> \brief Prints a BCSR matrix (block-style, not full)
!> \param[in] matrix          matrix
!> \param[in,out] error       cp2k error
!> \param[in] nodata          (optional) don't print actual data
! *****************************************************************************
  SUBROUTINE dbcsr_print(matrix, nodata, matlab_format, variable_name, unit_nr, error)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    LOGICAL, INTENT(IN), OPTIONAL            :: nodata, matlab_format
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: variable_name
    INTEGER, INTENT(IN), OPTIONAL            :: unit_nr
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_print', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=real_4), DIMENSION(:), &
      POINTER                                :: c_sp
    COMPLEX(KIND=real_8), DIMENSION(:), &
      POINTER                                :: c_dp
    INTEGER                                  :: bc, blk, blk_p, br, &
                                                error_handler, iunit, m, mn, n
    INTEGER, DIMENSION(:), POINTER           :: col_blk_offset, col_blk_size, &
                                                row_blk_offset, row_blk_size
    LOGICAL                                  :: my_matlab_format, tr, yesprint
    REAL(KIND=dp)                            :: blk_cs
    REAL(KIND=real_4), DIMENSION(:), POINTER :: r_sp
    REAL(KIND=real_8), DIMENSION(:), POINTER :: r_dp

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL dbcsr_assert (dbcsr_valid_index (matrix),&
         dbcsr_warning_level, dbcsr_caller_error,&
         routineN, "Can not print invalid matrix.",__LINE__,error)

    iunit=6
    IF (PRESENT(unit_nr)) iunit=unit_nr

    my_matlab_format = .FALSE.
    IF(PRESENT(matlab_format)) my_matlab_format = matlab_format
    yesprint = .TRUE.
    IF (PRESENT (nodata)) yesprint = .NOT. nodata
    WRITE(iunit,*)routineP//' Contents of matrix named ', matrix%m%name
    WRITE(iunit,*)routineP//' Flags ', matrix%m%symmetry,&
         matrix%m%negate_real, matrix%m%negate_imaginary, "type",&
         dbcsr_get_data_type (matrix), "serial", matrix%m%serial_number
    WRITE(iunit,'(1X,A,3(1X,I9,1X,A))')routineP, matrix%m%nblks, "blocks", &
         matrix%m%nze, "nzes,", dbcsr_get_data_size (matrix), "data els", &
         dbcsr_get_data_size_referenced (matrix), "used"
    WRITE(iunit,'(1X,A,I5,A,I5)')routineP//" Full size",&
         matrix%m%nfullrows_total, "x",matrix%m%nfullcols_total
    WRITE(iunit,'(1X,A,I5,A,I5)')routineP//" Blocked size",&
         matrix%m%nblkrows_total, "x",matrix%m%nblkcols_total
    SELECT CASE (matrix%m%data_type)
    CASE (dbcsr_type_real_8)
       CALL dbcsr_get_data (matrix%m%data_area, r_dp)
    CASE (dbcsr_type_real_4)
       CALL dbcsr_get_data (matrix%m%data_area, r_sp)
    CASE (dbcsr_type_complex_8)
       CALL dbcsr_get_data (matrix%m%data_area, c_dp)
    CASE (dbcsr_type_complex_4)
       CALL dbcsr_get_data (matrix%m%data_area, c_sp)
    END SELECT
    row_blk_size => array_data (matrix%m%row_blk_size)
    col_blk_size => array_data (matrix%m%col_blk_size)
    row_blk_offset => array_data (matrix%m%row_blk_offset)
    col_blk_offset => array_data (matrix%m%col_blk_offset)

    IF (matrix%m%nblks .GT. 0) THEN
       DO br = 1, matrix%m%nblkrows_total
          m = row_blk_size(br)
          DO blk = matrix%m%row_p(br)+1, matrix%m%row_p(br+1)
             bc = matrix%m%col_i(blk)
             n = col_blk_size(bc)
             mn = m*n
             blk_p = ABS(matrix%m%blk_p(blk))
             tr = matrix%m%blk_p(blk) .LT. 0
             block_exists: IF (blk_p .NE. 0) THEN
                IF (mn .GT. 0) THEN
                   SELECT CASE (matrix%m%data_type)
                   CASE (dbcsr_type_real_8)
                      blk_cs = REAL (DOT_PRODUCT(r_dp(blk_p:blk_p+mn-1),&
                           &                     r_dp(blk_p:blk_p+mn-1)),KIND=dp)
                      !CALL &
                      !     dbcsr_printmat(r_dp(blk_p:blk_p+mn-1),m,n, tr=tr)
                   CASE (dbcsr_type_real_4)
                      blk_cs = REAL (DOT_PRODUCT(r_sp(blk_p:blk_p+mn-1),&
                           &                     r_sp(blk_p:blk_p+mn-1)),KIND=dp)
                      !CALL &
                      !     dbcsr_printmat(r_sp(blk_p:blk_p+mn-1),m,n, tr=tr)
                   CASE (dbcsr_type_complex_8)
                      blk_cs = REAL (DOT_PRODUCT(c_dp(blk_p:blk_p+mn-1),&
                           &                     c_dp(blk_p:blk_p+mn-1)),KIND=dp)
                      !CALL &
                      !     dbcsr_printmat(c_dp(blk_p:blk_p+mn-1),m,n, tr=tr)
                   CASE (dbcsr_type_complex_4)
                      blk_cs = REAL (DOT_PRODUCT(c_sp(blk_p:blk_p+mn-1),&
                           &                     c_sp(blk_p:blk_p+mn-1)),KIND=dp)
                      !CALL &
                      !     dbcsr_printmat(c_sp(blk_p:blk_p+mn-1),m,n, tr=tr)
                   END SELECT
                ELSE
                   blk_cs = 0.0_dp
                ENDIF
                !WRITE(iunit,*)routineP//' chksum for (',br,',',bc,') at',&
                !     blk_p,'l',mn,'= ', blk_cs,'size',m,n
                IF(.not.my_matlab_format) WRITE(iunit,'(A,I6,",",I6,A,I7,A,I6,I6,"=",I7,A,E12.3)')&
                     !" Checksum for (",br,bc,") at ",blk_p," size ",m,n,mn,&
                     " Checksum for (",br,bc,") at ",matrix%m%blk_p(blk)," size ",m,n,mn,&
                     " checksum=",blk_cs
                IF (yesprint .AND. blk_p .NE. 0) THEN
                   IF (mn .GT. 0) THEN
                      SELECT CASE (matrix%m%data_type)
                      CASE (dbcsr_type_real_8)
                         !WRITE(iunit,'(10(1X,F7.2))')r_dp(blk_p:blk_p+mn-1)
                         IF(my_matlab_format) THEN
                            CALL dbcsr_printmat_matlab_d(r_dp(blk_p:blk_p+mn-1),m,n, &
                                 row_blk_offset(br), col_blk_offset(bc),iunit,tr=tr, &
                                 variable_name=variable_name)
                         ELSE
                            CALL dbcsr_printmat(r_dp(blk_p:blk_p+mn-1),m,n,iunit=iunit, tr=tr)
                         ENDIF
                      CASE (dbcsr_type_real_4)
                         IF(my_matlab_format) THEN
                            CALL dbcsr_printmat_matlab_s(r_sp(blk_p:blk_p+mn-1),m,n, &
                                 row_blk_offset(br), col_blk_offset(bc),iunit,tr=tr, &
                                 variable_name=variable_name)
                         ELSE
                            CALL dbcsr_printmat(r_sp(blk_p:blk_p+mn-1),m,n,iunit=iunit, tr=tr)
                         ENDIF
                      CASE (dbcsr_type_complex_8)
                         IF(my_matlab_format) THEN
                            CALL dbcsr_printmat_matlab_z(c_dp(blk_p:blk_p+mn-1),m,n, &
                                 row_blk_offset(br), col_blk_offset(bc),iunit,tr=tr, &
                                 variable_name=variable_name)
                         ELSE
                            CALL dbcsr_printmat(c_dp(blk_p:blk_p+mn-1),m,n,iunit=iunit, tr=tr)
                         ENDIF
                      CASE (dbcsr_type_complex_4)
                         IF(my_matlab_format) THEN
                            CALL dbcsr_printmat_matlab_c(c_sp(blk_p:blk_p+mn-1),m,n, &
                                 row_blk_offset(br), col_blk_offset(bc),iunit,tr=tr, &
                                 variable_name=variable_name)
                         ELSE
                            CALL dbcsr_printmat(c_sp(blk_p:blk_p+mn-1),m,n,iunit=iunit, tr=tr)
                         ENDIF
                      END SELECT
                   ENDIF
                ENDIF
             ENDIF block_exists
          ENDDO
       ENDDO
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_print


! *****************************************************************************
!> \brief Prints the elements of a matrix.
!> \param[in] rows,columns    is the logical (possibly detransposed) matrix
!>                            size, not the stored size.
!> \param[in] tr    (optional) specifies whether the elements are stored
!>                  transposed
! *****************************************************************************

  SUBROUTINE dbcsr_printmat_matlab_d(matrix, rows, cols, r_offset, c_offset, iunit, tr, variable_name)
    REAL(KIND=real_8), DIMENSION(:), &
      INTENT(IN)                             :: matrix
    INTEGER, INTENT(IN)                      :: rows, cols, r_offset, &
                                                c_offset, iunit
    LOGICAL, INTENT(IN), OPTIONAL            :: tr
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: variable_name

    INTEGER                                  :: c, c_off, m, n, r, r_off
    LOGICAL                                  :: t

!   ---------------------------------------------------------------------------

    m = rows
    n = cols
    r_off = r_offset
    c_off = c_offset
    t = .FALSE.
    IF(PRESENT(tr)) THEN
       IF(tr) THEN
          t = .TRUE.
          m = cols
          n = rows
          r_off = c_offset
          c_off = r_offset
       ENDIF
    ENDIF

    DO c = 1,cols
    DO r = 1,rows
       IF(.NOT.t) THEN
          IF(ABS(matrix(r+(c-1)*rows)).GT.1e-18_real_8) THEN
          IF(PRESENT(variable_name)) THEN
              WRITE(iunit,'(A,I4,A,I4,A,E23.16,A)') variable_name//'(',r+r_offset-1,',',c+c_offset-1,')=',matrix(r+(c-1)*rows),';'
          ELSE
             WRITE(iunit,'(A,I4,A,I4,A,E23.16,A)') 'a(',r+r_offset-1,',',c+c_offset-1,')=',matrix(r+(c-1)*rows),';'
          ENDIF
          ENDIF
       ELSE
          IF(ABS(matrix((r-1)*cols+c)).GT.1e-18_real_8) THEN
          IF(PRESENT(variable_name)) THEN
             WRITE(iunit,'(A,I4,A,I4,A,E23.16,A)') variable_name//'(',r+r_offset-1,',',c+c_offset-1,')=',matrix((r-1)*cols+c),';'
          ELSE
             WRITE(iunit,'(A,I4,A,I4,A,E23.16,A)') 'a(',r+r_offset-1,',',c+c_offset-1,')=',matrix((r-1)*cols+c),';'
          ENDIF
          ENDIF
       ENDIF
    ENDDO
    ENDDO
  END SUBROUTINE dbcsr_printmat_matlab_d

  SUBROUTINE dbcsr_printmat_matlab_s(matrix, rows, cols, r_offset, c_offset, iunit, tr, variable_name)
    REAL(KIND=real_4), DIMENSION(:), &
      INTENT(IN)                             :: matrix
    INTEGER, INTENT(IN)                      :: rows, cols, r_offset, &
                                                c_offset, iunit
    LOGICAL, INTENT(IN), OPTIONAL            :: tr
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: variable_name

    INTEGER                                  :: c, c_off, m, n, r, r_off
    LOGICAL                                  :: t

!   ---------------------------------------------------------------------------

    m = rows
    n = cols
    r_off = r_offset
    c_off = c_offset
    t = .FALSE.
    IF(PRESENT(tr)) THEN
       IF(tr) THEN
          t = .TRUE.
          m = cols
          n = rows
          r_off = c_offset
          c_off = r_offset
       ENDIF
    ENDIF

    DO c = 1,cols
    DO r = 1,rows
       IF(.NOT.t) THEN
          IF(ABS(matrix(r+(c-1)*rows)).GT.1e-8_real_4) THEN
          IF(PRESENT(variable_name)) THEN
              WRITE(iunit,'(A,I4,A,I4,A,E15.7,A)') variable_name//'(',r+r_offset-1,',',c+c_offset-1,')=',matrix(r+(c-1)*rows),';'
          ELSE
             WRITE(iunit,'(A,I4,A,I4,A,E15.7,A)') 'a(',r+r_offset-1,',',c+c_offset-1,')=',matrix(r+(c-1)*rows),';'
          ENDIF
          ENDIF
       ELSE
          IF(ABS(matrix((r-1)*cols+c)).GT.1e-8_real_4) THEN
          IF(PRESENT(variable_name)) THEN
             WRITE(iunit,'(A,I4,A,I4,A,E15.7,A)') variable_name//'(',r+r_offset-1,',',c+c_offset-1,')=',matrix((r-1)*cols+c),';'
          ELSE
             WRITE(iunit,'(A,I4,A,I4,A,E15.7,A)') 'a(',r+r_offset-1,',',c+c_offset-1,')=',matrix((r-1)*cols+c),';'
          ENDIF
          ENDIF
       ENDIF
    ENDDO
    ENDDO
  END SUBROUTINE dbcsr_printmat_matlab_s

  SUBROUTINE dbcsr_printmat_matlab_z(matrix, rows, cols, r_offset, c_offset, iunit, tr, variable_name)
    COMPLEX(KIND=real_8), DIMENSION(:), &
      INTENT(IN)                             :: matrix
    INTEGER, INTENT(IN)                      :: rows, cols, r_offset, &
                                                c_offset, iunit
    LOGICAL, INTENT(IN), OPTIONAL            :: tr
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: variable_name

    INTEGER                                  :: c, c_off, m, n, r, r_off
    LOGICAL                                  :: t

!   ---------------------------------------------------------------------------

    m = rows
    n = cols
    r_off = r_offset
    c_off = c_offset
    t = .FALSE.
    IF(PRESENT(tr)) THEN
       IF(tr) THEN
          t = .TRUE.
          m = cols
          n = rows
          r_off = c_offset
          c_off = r_offset
       ENDIF
    ENDIF

    DO c = 1,cols
       DO r = 1,rows
          IF(.NOT.t) THEN
             IF(PRESENT(variable_name)) THEN
                WRITE(iunit,'(A,I3,A,I3,A,E23.16,A,E23.16,A)') variable_name//'(',r+r_offset-1,',',c+c_offset-1,')=',&
                     REAL(matrix(r+(c-1)*rows)),'+',AIMAG(matrix(r+(c-1)*rows)),'i;'
             ELSE
                WRITE(iunit,'(A,I3,A,I3,A,E23.16,A,E23.16,A)') 'a(',r+r_offset-1,',',c+c_offset-1,')=',&
                     REAL(matrix(r+(c-1)*rows)),'+',AIMAG(matrix(r+(c-1)*rows)),'i;'
             ENDIF
          ELSE
             IF(PRESENT(variable_name)) THEN
                WRITE(iunit,'(A,I3,A,I3,A,E23.16,A,E23.16,A)') variable_name//'(',r+r_offset-1,',',c+c_offset-1,')=',&
                     REAL(matrix((r-1)*cols+c)),'+',AIMAG(matrix((r-1)*cols+c)),'i;'
             ELSE
                WRITE(iunit,'(A,I3,A,I3,A,E23.16,A,E23.16,A)') 'a(',r+r_offset-1,',',c+c_offset-1,')=',&
                     REAL(matrix((r-1)*cols+c)),'+',AIMAG(matrix((r-1)*cols+c)),'i;'
             ENDIF
          ENDIF
       ENDDO
    ENDDO
  END SUBROUTINE dbcsr_printmat_matlab_z

  SUBROUTINE dbcsr_printmat_matlab_c(matrix, rows, cols, r_offset, c_offset, iunit, tr, variable_name)
    COMPLEX(KIND=real_4), DIMENSION(:), &
      INTENT(IN)                             :: matrix
    INTEGER, INTENT(IN)                      :: rows, cols, r_offset, &
                                                c_offset, iunit
    LOGICAL, INTENT(IN), OPTIONAL            :: tr
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: variable_name

    INTEGER                                  :: c, c_off, m, n, r, r_off
    LOGICAL                                  :: t

!   ---------------------------------------------------------------------------

    m = rows
    n = cols
    r_off = r_offset
    c_off = c_offset
    t = .FALSE.
    IF(PRESENT(tr)) THEN
       IF(tr) THEN
          t = .TRUE.
          m = cols
          n = rows
          r_off = c_offset
          c_off = r_offset
       ENDIF
    ENDIF

    DO c = 1,cols
       DO r = 1,rows
          IF(.NOT.t) THEN
             IF(PRESENT(variable_name)) THEN
                WRITE(iunit,'(A,I3,A,I3,A,E15.7,A,E15.7,A)') variable_name//'(',r+r_offset-1,',',c+c_offset-1,')=',&
                     REAL(matrix(r+(c-1)*rows)),'+',AIMAG(matrix(r+(c-1)*rows)),'i;'
             ELSE
                WRITE(iunit,'(A,I3,A,I3,A,E15.7,A,E15.7,A)') 'a(',r+r_offset-1,',',c+c_offset-1,')=',&
                     REAL(matrix(r+(c-1)*rows)),'+',AIMAG(matrix(r+(c-1)*rows)),'i;'
             ENDIF
          ELSE
             IF(PRESENT(variable_name)) THEN
                WRITE(iunit,'(A,I3,A,I3,A,E15.7,A,E15.7,A)') variable_name//'(',r+r_offset-1,',',c+c_offset-1,')=',&
                     REAL(matrix((r-1)*cols+c)),'+',AIMAG(matrix((r-1)*cols+c)),'i;'
             ELSE
                WRITE(iunit,'(A,I3,A,I3,A,E15.7,A,E15.7,A)') 'a(',r+r_offset-1,',',c+c_offset-1,')=',&
                     REAL(matrix((r-1)*cols+c)),'+',AIMAG(matrix((r-1)*cols+c)),'i;'
             ENDIF
          ENDIF
       ENDDO
    ENDDO
  END SUBROUTINE dbcsr_printmat_matlab_c


! *****************************************************************************
!> \brief Prints the elements of a matrix.
!> \param[in] rows,columns    is the logical (possibly detransposed) matrix
!>                            size, not the stored size.
!> \param[in] tr    (optional) specifies whether the elements are stored
!>                  transposed
! *****************************************************************************
  SUBROUTINE printmat_s(matrix, rows, cols, iunit, title, tr)
    REAL(KIND=real_4), DIMENSION(:), &
      INTENT(IN)                             :: matrix
    INTEGER, INTENT(IN)                      :: rows, cols, iunit
    CHARACTER(*), INTENT(IN), OPTIONAL       :: title
    LOGICAL, INTENT(IN), OPTIONAL            :: tr

    CHARACTER(30)                            :: f
    INTEGER                                  :: m, n, r
    LOGICAL                                  :: t
    REAL(KIND=dp)                            :: bit_bucket

!   ---------------------------------------------------------------------------

    m = rows
    n = cols
    t = .FALSE.
    IF(PRESENT(title)) WRITE(iunit,*)title
    IF(PRESENT(tr)) THEN
       IF(tr) THEN
          t = .TRUE.
          m = cols
          n = rows
       ENDIF
    ENDIF
    DO r = LBOUND(matrix,1), UBOUND(matrix,1)
       bit_bucket = matrix(r)
    ENDDO
    bit_bucket = 0.0_dp
    DO r = LBOUND(matrix,1), UBOUND(matrix,1)
       bit_bucket = bit_bucket + matrix(r)
    ENDDO
    IF(m.GT.10000) m=0
    IF(n.GT.10000) n=0
    IF(m*n.LT.1 .OR. m*n.GT.SIZE(matrix)) RETURN
    WRITE(f,FMT="('(',I4,'(F9.4))')")cols
    DO r = 1, rows
       IF(.NOT.t) THEN
          WRITE(iunit,FMT=f)matrix(r:r+(cols-1)*rows:rows)
       ELSE
          WRITE(iunit,FMT=f)matrix((r-1)*cols+1:r*cols)
       ENDIF
    ENDDO
  END SUBROUTINE printmat_s

  SUBROUTINE printmat_d(matrix, rows, cols, iunit, title, tr)
    REAL(KIND=real_8), DIMENSION(:), &
      INTENT(IN)                             :: matrix
    INTEGER, INTENT(IN)                      :: rows, cols, iunit
    CHARACTER(*), INTENT(IN), OPTIONAL       :: title
    LOGICAL, INTENT(IN), OPTIONAL            :: tr

    IF (PRESENT (title)) THEN
       IF (PRESENT (tr)) THEN
          CALL printmat_s(REAL(matrix, KIND=sp), rows, cols, iunit, title, tr)
       ELSE
          CALL printmat_s(REAL(matrix, KIND=sp), rows, cols, iunit, title)
       ENDIF
    ELSE
       IF (PRESENT (tr)) THEN
          CALL printmat_s(REAL(matrix, KIND=sp), rows, cols, iunit, tr=tr)
       ELSE
          CALL printmat_s(REAL(matrix, KIND=sp), rows, cols, iunit)
       ENDIF
    ENDIF
  END SUBROUTINE printmat_d

  SUBROUTINE printmat_c(matrix, rows, cols, iunit, title, tr)
    COMPLEX(KIND=real_4), DIMENSION(:), &
      INTENT(IN)                             :: matrix
    INTEGER, INTENT(IN)                      :: rows, cols, iunit
    CHARACTER(*), INTENT(IN), OPTIONAL       :: title
    LOGICAL, INTENT(IN), OPTIONAL            :: tr

    IF (PRESENT (title)) THEN
       IF (PRESENT (tr)) THEN
          CALL printmat_s(REAL(matrix, KIND=sp), rows, cols, iunit, title, tr)
       ELSE
          CALL printmat_s(REAL(matrix, KIND=sp), rows, cols, iunit, title)
       ENDIF
    ELSE
       IF (PRESENT (tr)) THEN
          CALL printmat_s(REAL(matrix, KIND=sp), rows, cols, iunit, tr=tr)
       ELSE
          CALL printmat_s(REAL(matrix, KIND=sp), rows, cols, iunit)
       ENDIF
    ENDIF
  END SUBROUTINE printmat_c

  SUBROUTINE printmat_z(matrix, rows, cols, iunit, title, tr)
    COMPLEX(KIND=real_8), DIMENSION(:), &
      INTENT(IN)                             :: matrix
    INTEGER, INTENT(IN)                      :: rows, cols, iunit
    CHARACTER(*), INTENT(IN), OPTIONAL       :: title
    LOGICAL, INTENT(IN), OPTIONAL            :: tr

    IF (PRESENT (title)) THEN
       IF (PRESENT (tr)) THEN
          CALL printmat_s(REAL(matrix, KIND=sp), rows, cols, iunit, title, tr)
       ELSE
          CALL printmat_s(REAL(matrix, KIND=sp), rows, cols, iunit, title)
       ENDIF
    ELSE
       IF (PRESENT (tr)) THEN
          CALL printmat_s(REAL(matrix, KIND=sp), rows, cols, iunit, tr=tr)
       ELSE
          CALL printmat_s(REAL(matrix, KIND=sp), rows, cols, iunit)
       ENDIF
    ENDIF
  END SUBROUTINE printmat_z

! *****************************************************************************
!> \brief Prints the values from all processors
!> \param[in] val   local value to print
!> \param[in] desc  description of value
!> \param[in] mp    message-passing environment
! *****************************************************************************
  SUBROUTINE all_print (val, desc, mp)
    INTEGER, INTENT(IN)                      :: val
    CHARACTER(LEN=*), INTENT(IN)             :: desc
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp

    CHARACTER(LEN=80)                        :: fmt
    INTEGER                                  :: n
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: infos

!

    n = dbcsr_mp_numnodes (mp)
    WRITE(fmt,'("(1X,A,",I5,"(1X,I4))")')n
    ALLOCATE (infos(n))
    CALL mp_allgather(val, infos, dbcsr_mp_group(mp))
    WRITE(*,fmt)desc,infos
  END SUBROUTINE all_print
  
! *****************************************************************************
!> \brief Writes a DBCSR matrix in a file
!>   file's header consists of 3 sub-headers:
!>   sub-header1 contains: 
!>     1 string: (of length version_len) the current version of this routine, 
!>     1 string: (of length default_string_length) matrix_name, 
!>     1 character: matrix_type, 
!>     4 integers: numnodes, data_type, nblkrows_total, nblkcols_total, 
!>     2 vectors:  row_blk_size (length = nblkrows_total), 
!>                 col_blk_size (length = nblkcols_total), 
!>   sub-header2 contains:
!>     2 integers: nblks, data_area_size, 
!>   sub-header3 contains: 
!>     3 vectors:  row_p (length = nblkrows_total+1),  
!>                 col_i (length = nblks), 
!>                 blk_p (length = nblks); 
!>   and the file's body contains the block data
!>
!> \param[in] matrix    DBCSR matrix
!> \param[in] filepath  path to the file
!>
!> \par History
!>      11.2012 created [Hossein Bani-Hashemian]
!> \author Hossein Bani-Hashemian
!> \version DBCSRv_1.0
! *****************************************************************************
  SUBROUTINE dbcsr_binary_write(matrix, filepath, error)

    IMPLICIT NONE
   
    TYPE(dbcsr_obj), INTENT(IN)           :: matrix
    CHARACTER (len=*), INTENT(IN)         :: filepath
    TYPE(dbcsr_error_type), INTENT(INOUT) :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_binary_write', &
                                   routineP = moduleN//':'//routineN

    INTEGER                               :: nblkrows_total, nblkcols_total, &
                                             nblks, size_of_pgrid, &
                                             thefile, i, sendbuf, data_area_size, &
                                             data_type, type_size, &
                                             mp_group, mynode, numnodes, & 
                                             ginfo_size, linfo_size, error_handler                
    INTEGER, DIMENSION(:), POINTER        :: row_p, col_i, blk_p,&
                                             row_blk_size, col_blk_size
    INTEGER, DIMENSION(:,:), POINTER      :: pgrid                              
    TYPE(mp_type_descriptor_type)         :: mp_type
    TYPE(dbcsr_mp_obj)                    :: mp_env
    TYPE(dbcsr_distribution_obj)          :: distribution
    TYPE(dbcsr_data_obj)                  :: data_area
    COMPLEX(sp),DIMENSION(:),POINTER      :: c_sp
    COMPLEX(dp),DIMENSION(:),POINTER      :: c_dp
    REAL(sp),DIMENSION(:),POINTER         :: r_sp
    REAL(dp),DIMENSION(:),POINTER         :: r_dp
    CHARACTER                             :: matrix_type
    CHARACTER(LEN=80)                     :: matrix_name_v_1_0
    CHARACTER(LEN=default_string_length)  :: matrix_name
    INTEGER, PARAMETER                    :: version_len=10
    CHARACTER(LEN=version_len), PARAMETER :: version = "DBCSRv_1.0"   
    INTEGER, ALLOCATABLE, DIMENSION(:) :: linfo_sizes, da_sizes
    INTEGER(kind=file_offset), ALLOCATABLE, DIMENSION(:) :: bdata_disps, bdata_offsets, & 
                                                          subh2_disps, subh2_offsets, & 
                                                          subh3_disps, subh3_offsets
    INTEGER(kind=file_offset), PARAMETER                 :: BOF = 0
    INTEGER, PARAMETER          :: char_count = version_len+default_string_length+1 !version, matrix_name, matrix_type

   CALL dbcsr_error_set(routineN, error_handler, error)

    CALL dbcsr_assert (default_string_length==80,&
         dbcsr_warning_level, dbcsr_caller_error,&
         routineN, "Changing the default string length affects the format of the writen matrix. Version needs to be adjusted",&
         __LINE__,error)

    nblkrows_total = dbcsr_nblkrows_total(matrix)
    nblkcols_total = dbcsr_nblkcols_total(matrix)
    distribution = dbcsr_distribution(matrix) 
    matrix_name =  dbcsr_name(matrix)
    data_area = dbcsr_data_area(matrix%m)
    matrix_type = dbcsr_get_matrix_type(matrix)
    data_type = dbcsr_get_data_type(matrix)
    mp_env = dbcsr_distribution_mp(distribution)
    mp_group = dbcsr_mp_group(mp_env)
    nblks = dbcsr_get_num_blocks(matrix)
    row_p => matrix%m%row_p
    col_i => matrix%m%col_i
    blk_p => matrix%m%blk_p
    row_blk_size => array_data (matrix%m%row_blk_size)
    col_blk_size => array_data (matrix%m%col_blk_size)
    pgrid => dbcsr_mp_pgrid(mp_env)
    size_of_pgrid = SIZE(pgrid)
 
    CALL mp_environ(numnodes,mynode,mp_group)

    ALLOCATE(linfo_sizes(numnodes), da_sizes(numnodes), &
             subh2_disps(numnodes), subh2_offsets(numnodes),&
             subh3_disps(numnodes), subh3_offsets(numnodes),&
             bdata_disps(numnodes), bdata_offsets(numnodes))
    subh2_disps(:) = (/((i-1)*2,i=1,numnodes)/)
    subh3_disps = BOF
    bdata_disps = BOF
    linfo_sizes = BOF
    subh2_offsets = BOF
    subh3_offsets = BOF
    bdata_offsets = BOF
    da_sizes = BOF

    ginfo_size = char_count+4+nblkrows_total+nblkcols_total
    linfo_size = 1+nblkrows_total+2*nblks

    sendbuf = linfo_size
    CALL mp_allgather(sendbuf,linfo_sizes,mp_group)
    CALL cumsum_l(INT(linfo_sizes,kind=file_offset),subh3_disps)
    subh3_disps(:) = CSHIFT(subh3_disps,shift=-1)+ginfo_size+2*numnodes
    subh3_disps(1)=ginfo_size+2*numnodes

    data_area_size = dbcsr_get_data_size_referenced(matrix)
    sendbuf = data_area_size
    CALL mp_allgather(sendbuf,da_sizes,mp_group) 
    CALL cumsum_l(INT(da_sizes,kind=file_offset),bdata_disps)
    bdata_disps(:) = CSHIFT(bdata_disps,shift=-1)+SUM(INT(linfo_sizes,KIND=file_offset))+&
                    ginfo_size+numnodes*2
    bdata_disps(1)=SUM(INT(linfo_sizes,KIND=file_offset))+ginfo_size+numnodes*2

    CALL mp_file_open(mp_group, thefile, filepath, file_amode_create+file_amode_wronly)

    IF (mynode .EQ. 0) THEN
       CALL mp_file_write_at(thefile, BOF, version)
       matrix_name_v_1_0 = matrix_name
       CALL mp_file_write_at(thefile, BOF+version_len*mpi_character_size, matrix_name_v_1_0 ) 
       CALL mp_file_write_at(thefile, BOF+(version_len+default_string_length)*mpi_character_size, matrix_type) 
       CALL mp_file_write_at(thefile, BOF+char_count*mpi_character_size, &
                             (/size_of_pgrid, data_type, &
                              nblkrows_total, nblkcols_total, &
                              row_blk_size, col_blk_size/))             
    END IF
!> write sub-header2 
    subh2_disps = subh2_disps+ginfo_size
    subh2_offsets(:) = BOF + (subh2_disps-char_count)*mpi_integer_size +&
                    char_count*mpi_character_size
    CALL mp_file_write_at_all(thefile, subh2_offsets(mynode+1), (/nblks, data_area_size/))
!> write sub-header3 
    subh3_offsets(:) = BOF + (subh3_disps-char_count)*mpi_integer_size +&
                    char_count*mpi_character_size
    CALL mp_file_write_at_all(thefile, subh3_offsets(mynode+1), (/row_p, col_i, blk_p/))
!> write block data 
    mp_type = dbcsr_mp_type_from_anytype(data_area)
    CALL mp_type_size(mp_type,type_size)
    bdata_offsets(:) = BOF + (/((bdata_disps(i)-bdata_disps(1))*type_size,i=1,numnodes)/)+&
                     (bdata_disps(1)-char_count)*mpi_integer_size+ &
                     char_count*mpi_character_size
    SELECT CASE (data_type)
    CASE (dbcsr_type_real_4)
       r_sp => data_area%d%r_sp
       CALL mp_file_write_at_all(thefile, bdata_offsets(mynode+1), r_sp, msglen = data_area_size) 
    CASE (dbcsr_type_real_8)
       r_dp => data_area%d%r_dp
       CALL mp_file_write_at_all(thefile, bdata_offsets(mynode+1), r_dp, msglen = data_area_size)
    CASE (dbcsr_type_complex_4)
       c_sp => data_area%d%c_sp
       CALL mp_file_write_at_all(thefile, bdata_offsets(mynode+1), c_sp, msglen = data_area_size) 
    CASE (dbcsr_type_complex_8)
       c_dp => data_area%d%c_dp
       CALL mp_file_write_at_all(thefile, bdata_offsets(mynode+1), c_dp, msglen = data_area_size) 
    END SELECT
    CALL mp_file_close(thefile)

    DEALLOCATE(linfo_sizes,da_sizes)
    DEALLOCATE(subh2_disps,subh2_offsets,subh3_disps,subh3_offsets)
    DEALLOCATE(bdata_disps,bdata_offsets)

    CALL dbcsr_error_stop(error_handler, error)
    
    CONTAINS
      SUBROUTINE cumsum_l(arr,cumsum)
    INTEGER(kind=file_offset), &
      DIMENSION(:), INTENT(IN)               :: arr
    INTEGER(kind=file_offset), &
      DIMENSION(SIZE(arr)), INTENT(OUT)      :: cumsum

    INTEGER                                  :: i

        cumsum(1)=arr(1)
        DO i=2,SIZE(arr) 
          cumsum(i)=cumsum(i-1)+arr(i)
        END DO
      END SUBROUTINE cumsum_l
  END SUBROUTINE dbcsr_binary_write

! *****************************************************************************
!> \brief Reads a DBCSR matrix from a file 
!> \param[in] filepath             path to the file
!> \param[in] distribution         row and column distribution
!> \param[in](optional) groupid    message passing environment identifier
!> \param[out] matrix_new          DBCSR matrix
!>
!> \par History
!>      11.2012 created [Hossein Bani-Hashemian]
!> \author Hossein Bani-Hashemian
!> \version DBCSRv_1.0
! *****************************************************************************
  SUBROUTINE dbcsr_binary_read(filepath, distribution, groupid, matrix_new, error)

    IMPLICIT NONE

    CHARACTER (len=*), INTENT(IN)                :: filepath
    TYPE(dbcsr_distribution_obj), INTENT(IN)     :: distribution
    INTEGER, INTENT(IN), OPTIONAL                :: groupid
    TYPE(dbcsr_obj), INTENT(INOUT)               :: matrix_new
    TYPE(dbcsr_error_type), INTENT(INOUT)        :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_binary_read', &
                                   routineP = moduleN//':'//routineN

    INTEGER                               :: nblkrows_total, nblkcols_total, &
                                             nblks, darea_size, data_type, type_size, &
                                             globalinfo_size, &
                                             size_of_pgrid, &
                                             thefile, i, j, &
                                             nblocks, &
                                             share_size, order, cur_blks,  &
                                             job_count, start_index, end_index,&
                                             localinfo_length, blockdata_length,&
                                             group_id, worker_id, group_list_size, error_handler, linfo_length
    CHARACTER                             :: matrix_type
    CHARACTER(LEN=default_string_length)  :: matrix_name
    INTEGER, PARAMETER                    :: version_len=10
    CHARACTER(LEN=version_len)            :: version   
    CHARACTER(LEN=80)                     :: matrix_name_v_1_0
    CHARACTER(LEN=version_len), PARAMETER :: version_v_1_0="DBCSRv_1.0"

    INTEGER, DIMENSION(:), POINTER        :: row_p, col_i, blk_p, ginfo,&
                                             proc_nblks, proc_darea_sizes
    INTEGER, DIMENSION(4)                 :: values
    INTEGER, ALLOCATABLE, DIMENSION(:)    :: linfo_lens, bdata_lens
    INTEGER, ALLOCATABLE, DIMENSION(:), TARGET     :: ginfo_vec, linfo_vec, &
                                                      rowp, coli, blkp
    INTEGER, ALLOCATABLE, DIMENSION(:,:), TARGET   :: val_data
    TYPE(array_i1d_obj)                            :: row_blk_size, col_blk_size
    TYPE(dbcsr_mp_obj)                             :: mp_env
    TYPE(dbcsr_data_obj)                           :: dblk
    REAL(sp)                                       :: rsp_dummy(1)
    REAL(dp)                                       :: rdp_dummy(1)
    COMPLEX(sp)                                    :: csp_dummy(1)
    COMPLEX(dp)                                    :: cdp_dummy(1)
    REAL(sp), ALLOCATABLE, DIMENSION(:),TARGET     :: rsp
    REAL(dp), ALLOCATABLE, DIMENSION(:),TARGET     :: rdp
    COMPLEX(sp), ALLOCATABLE, DIMENSION(:),TARGET  :: csp
    COMPLEX(dp), ALLOCATABLE, DIMENSION(:),TARGET  :: cdp
    INTEGER(kind=file_offset), ALLOCATABLE, DIMENSION(:)   :: subh2_offsets, &
                                                            subh3_disps, subh3_offsets, &
                                                            bdata_disps, bdata_offsets
    INTEGER(kind=file_offset), PARAMETER    :: BOF = 0
    INTEGER(kind=file_offset)               :: offset, subh2_start, subh3_start, bdata_start, file_size,&
                                               localinfo_offset, blockdata_offset, sum_nblks, subh3_length, data_area_size
    INTEGER, PARAMETER                      :: char_count = 1+version_len+default_string_length

    CALL dbcsr_error_set(routineN, error_handler, error)


    IF (PRESENT(groupid)) THEN
       group_id = groupid
    ELSE
       mp_env = dbcsr_distribution_mp(distribution)
       group_id = dbcsr_mp_group(mp_env)
    END IF

    CALL mp_environ(group_list_size,worker_id,group_id)
    CALL mp_file_open(group_id, thefile, filepath, file_amode_rdonly)
    
!> read version, matrix name and matrix type
    CALL mp_file_read_at_all(thefile, BOF, version)

    CALL dbcsr_assert (version==version_v_1_0,&
         dbcsr_warning_level, dbcsr_caller_error,&
         routineN, "Trying to read an unknown version of the matrix data file. Good luck!",&
         __LINE__,error)

    CALL mp_file_read_at_all(thefile, BOF+version_len*mpi_character_size, matrix_name_v_1_0) 
    matrix_name=matrix_name_v_1_0

    CALL mp_file_read_at_all(thefile, BOF+(version_len+default_string_length)*mpi_character_size, matrix_type) 
!> read 4 integer values form sub-header1
    CALL mp_file_read_at_all(thefile, BOF+char_count*mpi_character_size, values) 
    size_of_pgrid = values(1) 
    data_type = values(2)
    nblkrows_total = values(3)  
    nblkcols_total = values(4)
!> read 2 vectors, row_blk_size and col_blk_size, from sub-header1   
    globalinfo_size = nblkrows_total+nblkcols_total
    ALLOCATE(ginfo_vec(globalinfo_size))
    CALL mp_file_read_at_all(thefile, BOF+char_count*mpi_character_size+4*mpi_integer_size, ginfo_vec)
    ginfo => ginfo_vec(1:nblkrows_total)
    CALL array_new (row_blk_size, ginfo)
    ginfo => ginfo_vec(nblkrows_total+1:globalinfo_size)
    CALL array_new (col_blk_size, ginfo)

!> compute the offsets where sub-header2 and sub-header3 start 
    subh2_start = (4+globalinfo_size)*mpi_integer_size+char_count*mpi_character_size
    subh3_start = subh2_start+2*size_of_pgrid*mpi_integer_size
 
!> compute the offsets in sub-header2 and read 2 integers nblocks, data_area_size
    ! number of data chunks from sub-header 2 and 3 to be read by every node rounded up 
    ! to the next integer to make it even for all the nodes in the specified mpi group 
    share_size = CEILING(REAL(size_of_pgrid)/group_list_size) 

    ALLOCATE(subh2_offsets(share_size))
    subh2_offsets = BOF
    DO i = 1, share_size
      offset = subh2_start+mpi_integer_size*2*(worker_id+(i-1)*group_list_size)
      IF (offset .GE. subh3_start) EXIT
      subh2_offsets(i) = offset
    END DO

    ALLOCATE(val_data(3,share_size))
    val_data(:,:) = 0
    DO i = 1, share_size
       CALL mp_file_read_at_all(thefile, subh2_offsets(i), values, msglen=2)
       nblocks = values(1)    
       data_area_size = values(2)
       IF (subh2_offsets(i) .EQ. 0 ) EXIT
       val_data(1,i) = nblocks  
       val_data(2,i) = data_area_size
       val_data(3,i) = worker_id+(i-1)*group_list_size+1   ! order 
       ! order = indices of an array of length size_of_pgrid to be accessed by the current node 
    END DO
    nblks = SUM(val_data(1,:))
    darea_size = SUM(val_data(2,:))
    proc_nblks => val_data(1,:)           ! to be passed to dbcsr_datablock_redistribute
    proc_darea_sizes => val_data(2,:)     ! to be passed to dbcsr_datablock_redistribute

!> compute the offsets in sub-header3 and read 3 vectors row_p, col_i, blk_p
    ! actual number of chunks to be read by the current node 
    job_count = COUNT(val_data(3,:) .NE. 0)
    CALL mp_file_get_size(thefile, file_size)

    ALLOCATE(linfo_lens(size_of_pgrid))
    ALLOCATE(subh3_disps(size_of_pgrid))
    ALLOCATE(subh3_offsets(size_of_pgrid))
    linfo_lens = 0; subh3_disps = 0
    DO i = 1, size_of_pgrid
       DO j = 1, share_size
          order = val_data(3,j)
          IF (i .EQ. order) linfo_lens(order) = &
                            1+nblkrows_total+2*val_data(1,j)
       END DO
    END DO
    CALL mp_sum(linfo_lens,group_id)
    CALL cumsum_l(INT(linfo_lens,kind=file_offset),subh3_disps)
    subh3_disps(:) = CSHIFT(subh3_disps,shift=-1)
    subh3_disps(1) = BOF
    subh3_offsets(:) = subh3_start+subh3_disps*mpi_integer_size

    sum_nblks = INT(nblks,kind=file_offset)  
    CALL mp_sum(sum_nblks,group_id)
    subh3_length = size_of_pgrid*INT(1+nblkrows_total,KIND=file_offset) + 2*sum_nblks

    linfo_length = nblkrows_total+1 + 2*MAXVAL(val_data(1,:))

    ALLOCATE(linfo_vec(linfo_length))
    ALLOCATE(rowp((nblkrows_total+1)*job_count))
    ALLOCATE(coli(nblks))
    ALLOCATE(blkp(nblks))
    DO i = 1, share_size
       order = val_data(3,i)
       cur_blks = val_data(1,i)
       IF (order .EQ. 0) THEN
          localinfo_offset = file_size
          localinfo_length = 0
       ELSE
          localinfo_offset = subh3_offsets(order)
          localinfo_length = linfo_lens(order)
       END IF
       CALL mp_file_read_at_all(thefile, localinfo_offset, linfo_vec, msglen=localinfo_length)
       IF (localinfo_length .EQ. 0) EXIT

       rowp((i-1)*(nblkrows_total+1)+1:i*(nblkrows_total+1)) = linfo_vec(1:nblkrows_total+1)
       start_index = SUM(val_data(1,1:i-1))+1
       end_index = SUM(val_data(1,1:i))
       coli(start_index:end_index) = &
                linfo_vec(nblkrows_total+2:cur_blks+nblkrows_total+1)
       blkp(start_index:end_index) = &
                linfo_vec(cur_blks+nblkrows_total+2:2*cur_blks+nblkrows_total+1)
    END DO
    row_p => rowp
    col_i => coli
    blk_p => blkp 

!> compute the offsets and read block data
    ALLOCATE(bdata_lens(size_of_pgrid))
    ALLOCATE(bdata_disps(size_of_pgrid))
    ALLOCATE(bdata_offsets(size_of_pgrid))
    bdata_lens = 0
    DO i = 1, size_of_pgrid
       DO j = 1, share_size
          order = val_data(3,j)
          IF (i .EQ. order) bdata_lens(order) = val_data(2,j)
       END DO
    END DO
    CALL mp_sum(bdata_lens,group_id)
    CALL cumsum_l(INT(bdata_lens,kind=file_offset),bdata_disps)
    bdata_disps(:) = CSHIFT(bdata_disps,shift=-1)
    bdata_disps(1) = BOF

    bdata_start = subh3_start+subh3_length*mpi_integer_size
    SELECT CASE (data_type)
    CASE (dbcsr_type_real_4)
         type_size = real_4_size
    CASE (dbcsr_type_real_8)
         type_size = real_8_size
    CASE (dbcsr_type_complex_4)
         type_size = 2*real_4_size
    CASE (dbcsr_type_complex_8)
         type_size = 2*real_8_size
    END SELECT
    bdata_offsets(:) = bdata_start+bdata_disps*type_size  

    SELECT CASE (data_type)
    CASE (dbcsr_type_real_4)
         ALLOCATE(rsp(darea_size))
         DO i = 1, share_size
            order = val_data(3,i)
            !> use dummy one-sized data array as buffer in place of empty array 
            !> when nothing is supposed to be read (order = 0)  
            IF (order .EQ. 0) THEN
               blockdata_offset = file_size
               CALL mp_file_read_at_all(thefile, blockdata_offset, rsp_dummy)
            ELSE
               start_index = SUM(val_data(2,1:i-1))+1
               end_index = SUM(val_data(2,1:i))
               blockdata_length = bdata_lens(order)
               blockdata_offset = bdata_offsets(order)
               CALL mp_file_read_at_all(thefile, blockdata_offset, rsp(start_index:end_index), &
                                        msglen = blockdata_length)
            END IF
         END DO
    CASE (dbcsr_type_real_8)
         ALLOCATE(rdp(darea_size))
         DO i = 1, share_size
            order = val_data(3,i)
            IF (order .EQ. 0)  THEN
               blockdata_offset = file_size
               CALL mp_file_read_at_all(thefile, blockdata_offset, rdp_dummy)
            ELSE
               start_index = SUM(val_data(2,1:i-1))+1
               end_index = SUM(val_data(2,1:i))
               blockdata_length = bdata_lens(order)
               blockdata_offset = bdata_offsets(order)
               CALL mp_file_read_at_all(thefile, blockdata_offset, rdp(start_index:end_index), &
                                        msglen = blockdata_length)
            END IF
         END DO
    CASE (dbcsr_type_complex_4)
         ALLOCATE(csp(darea_size))
         DO i = 1, share_size
            order = val_data(3,i)
            IF (order .EQ. 0) THEN
               blockdata_offset = file_size
               CALL mp_file_read_at_all(thefile, blockdata_offset, csp_dummy)
            ELSE
               start_index = SUM(val_data(2,1:i-1))+1
               end_index = SUM(val_data(2,1:i))
               blockdata_length = bdata_lens(order)
               blockdata_offset = bdata_offsets(order)
               CALL mp_file_read_at_all(thefile, blockdata_offset, csp(start_index:end_index), &
                                        msglen = blockdata_length)
            END IF
         END DO
    CASE (dbcsr_type_complex_8)
         ALLOCATE(cdp(darea_size))
         DO i = 1, share_size
            order = val_data(3,i)
            IF (order .EQ. 0) THEN
               blockdata_offset = file_size
               CALL mp_file_read_at_all(thefile, blockdata_offset, cdp_dummy)
            ELSE
               start_index = SUM(val_data(2,1:i-1))+1
               end_index = SUM(val_data(2,1:i))
               blockdata_length = bdata_lens(order)
               blockdata_offset = bdata_offsets(order)
               CALL mp_file_read_at_all(thefile, blockdata_offset, cdp(start_index:end_index), &
                                        msglen = blockdata_length)
            END IF
         END DO
    END SELECT
    CALL dbcsr_data_init (dblk)
    CALL dbcsr_data_new (dblk, data_type)
    dblk%d%r_dp => rdp
    dblk%d%r_sp => rsp
    dblk%d%c_dp => cdp
    dblk%d%c_sp => csp

    CALL mp_file_close(thefile)

    CALL dbcsr_create (matrix_new, matrix_name, distribution, matrix_type, &
                       row_blk_size, col_blk_size, nblks=nblks, nze=darea_size, & 
                       data_type=data_type, error=error)
    CALL dbcsr_datablock_redistribute (dblk, row_p, col_i, blk_p, proc_nblks, proc_darea_sizes, matrix_new, error)

    DEALLOCATE(subh2_offsets, subh3_offsets, bdata_offsets)
    DEALLOCATE(subh3_disps, bdata_disps)
    DEALLOCATE(linfo_lens, bdata_lens)
    DEALLOCATE(val_data, ginfo_vec, linfo_vec)
    DEALLOCATE(rowp, coli, blkp)
    IF (ALLOCATED(rdp)) DEALLOCATE(rdp)
    IF (ALLOCATED(rsp)) DEALLOCATE(rsp)
    IF (ALLOCATED(cdp)) DEALLOCATE(cdp)
    IF (ALLOCATED(csp)) DEALLOCATE(csp)
    CALL dbcsr_data_clear_pointer (dblk) 
    DEALLOCATE (dblk%d)
    CALL array_release(row_blk_size)
    CALL array_release(col_blk_size)

    CALL dbcsr_error_stop(error_handler, error)
    CONTAINS
      SUBROUTINE cumsum_l(arr,cumsum)
    INTEGER(kind=file_offset), &
      DIMENSION(:), INTENT(IN)               :: arr
    INTEGER(kind=file_offset), &
      DIMENSION(SIZE(arr)), INTENT(OUT)      :: cumsum

    INTEGER                                  :: i

        cumsum(1)=arr(1)
        DO i=2,SIZE(arr) 
          cumsum(i)=cumsum(i-1)+arr(i)
        END DO
      END SUBROUTINE cumsum_l

  END SUBROUTINE dbcsr_binary_read

! *****************************************************************************
!> \brief Prints the sum of elements for each block
!> \param[in] matrix          matrix
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_print_block_sum(matrix, unit_nr, error)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER, INTENT(IN), OPTIONAL            :: unit_nr
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_print_block_sum', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=real_4)                     :: blk_sum_c_sp
    COMPLEX(KIND=real_4), DIMENSION(:), &
      POINTER                                :: c_sp
    COMPLEX(KIND=real_8)                     :: blk_sum_c_dp
    COMPLEX(KIND=real_8), DIMENSION(:), &
      POINTER                                :: c_dp
    INTEGER                                  :: bc, blk, blk_p, br, &
                                                error_handler, iunit, m, mn, n
    INTEGER, DIMENSION(:), POINTER           :: col_blk_offset, col_blk_size, &
                                                row_blk_offset, row_blk_size
    REAL(KIND=real_4)                        :: blk_sum_r_sp
    REAL(KIND=real_4), DIMENSION(:), POINTER :: r_sp
    REAL(KIND=real_8)                        :: blk_sum_r_dp
    REAL(KIND=real_8), DIMENSION(:), POINTER :: r_dp

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL dbcsr_assert (dbcsr_valid_index (matrix),&
         dbcsr_warning_level, dbcsr_caller_error,&
         routineN, "Can not print invalid matrix.",__LINE__,error)

    iunit=6
    IF (PRESENT(unit_nr)) iunit=unit_nr

    SELECT CASE (matrix%m%data_type)
    CASE (dbcsr_type_real_8)
       CALL dbcsr_get_data (matrix%m%data_area, r_dp)
    CASE (dbcsr_type_real_4)
       CALL dbcsr_get_data (matrix%m%data_area, r_sp)
    CASE (dbcsr_type_complex_8)
       CALL dbcsr_get_data (matrix%m%data_area, c_dp)
    CASE (dbcsr_type_complex_4)
       CALL dbcsr_get_data (matrix%m%data_area, c_sp)
    END SELECT
    row_blk_size => array_data (matrix%m%row_blk_size)
    col_blk_size => array_data (matrix%m%col_blk_size)
    row_blk_offset => array_data (matrix%m%row_blk_offset)
    col_blk_offset => array_data (matrix%m%col_blk_offset)

    IF (matrix%m%nblks .GT. 0) THEN
       DO br = 1, matrix%m%nblkrows_total
          m = row_blk_size(br)
          DO blk = matrix%m%row_p(br)+1, matrix%m%row_p(br+1)
             bc = matrix%m%col_i(blk)
             n = col_blk_size(bc)
             mn = m*n
             blk_p = ABS(matrix%m%blk_p(blk))
             block_exists: IF (blk_p .NE. 0) THEN
                IF (mn .GT. 0) THEN
                   SELECT CASE (matrix%m%data_type)
                   CASE (dbcsr_type_real_8)
                      blk_sum_r_dp = SUM(r_dp(blk_p:blk_p+mn-1))
                      WRITE(iunit,'("Block (",I6,",",I6,") sum = ",ES18.9)')&
                         br,bc,blk_sum_r_dp
                   CASE (dbcsr_type_real_4)
                      blk_sum_r_sp = SUM(r_sp(blk_p:blk_p+mn-1))
                      WRITE(iunit,'(" (",I6,",",I6,") sum = ",ES18.9)')&
                         br,bc,blk_sum_r_sp
                   CASE (dbcsr_type_complex_8)
                      blk_sum_c_dp = SUM(c_dp(blk_p:blk_p+mn-1))
                      WRITE(iunit,'(" (",I6,",",I6,") sum = ",ES18.9," I*",ES18.9)')&
                         br,bc,REAL(blk_sum_c_dp),AIMAG(blk_sum_c_dp)
                   CASE (dbcsr_type_complex_4)
                      blk_sum_c_sp = SUM(c_sp(blk_p:blk_p+mn-1))
                      WRITE(iunit,'(" (",I6,",",I6,") sum = ",ES18.9," I*",ES18.9)')&
                         br,bc,REAL(blk_sum_c_sp),AIMAG(blk_sum_c_sp)
                   END SELECT
                ELSE
                   blk_sum_r_dp = 0.0_dp
                   WRITE(iunit,'(" (",I6,",",I6,") sum = ",ES18.9)')&
                      br,bc,blk_sum_r_dp
                ENDIF
             ENDIF block_exists
          ENDDO
       ENDDO
    ENDIF

    CALL dbcsr_error_stop(error_handler, error)
  
  END SUBROUTINE dbcsr_print_block_sum

END MODULE dbcsr_io
