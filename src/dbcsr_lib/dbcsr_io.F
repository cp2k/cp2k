!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2011  CP2K developers group                          !
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

  USE array_types,                     ONLY: array_data
  USE dbcsr_data_methods,              ONLY: dbcsr_get_data
  USE dbcsr_error_handling
  USE dbcsr_kinds,                     ONLY: dp,&
                                             real_4,&
                                             real_8,&
                                             sp
  USE dbcsr_message_passing,           ONLY: mp_allgather
  USE dbcsr_methods,                   ONLY: dbcsr_get_data_size,&
                                             dbcsr_get_data_size_used,&
                                             dbcsr_get_data_type,&
                                             dbcsr_mp_group,&
                                             dbcsr_mp_numnodes,&
                                             dbcsr_valid_index
  USE dbcsr_types,                     ONLY: dbcsr_mp_obj,&
                                             dbcsr_obj,&
                                             dbcsr_type_complex_4,&
                                             dbcsr_type_complex_8,&
                                             dbcsr_type_real_4,&
                                             dbcsr_type_real_8

  !$ USE OMP_LIB

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_io'

  REAL, PARAMETER                      :: default_resize_factor = 1.618034



  INTERFACE dbcsr_printmat
     MODULE PROCEDURE printmat_s, printmat_d, printmat_c, printmat_z
     MODULE PROCEDURE printmat_s_2d, printmat_d_2d, printmat_c_2d, printmat_z_2d
  END INTERFACE

  ! Main
  PUBLIC :: dbcsr_print
  ! Low-level printing
! PUBLIC :: dbcsr_printmat, dbcsr_print2dmat
  ! Utility printing
  PUBLIC :: all_print
  PUBLIC :: print_xfer_timings


  LOGICAL, PARAMETER :: bcsr_debug =   .TRUE.
  LOGICAL, PARAMETER :: bcsr_info =    .FALSE.
  LOGICAL, PARAMETER :: bcsr_verbose = .FALSE.

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
         dbcsr_get_data_type (matrix)
    WRITE(iunit,'(1X,A,3(1X,I9,1X,A))')routineP, matrix%m%nblks, "blocks", &
         matrix%m%nze, "nzes,", dbcsr_get_data_size (matrix), "data els", &
         dbcsr_get_data_size_used (matrix, error=error), "used"
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

  SUBROUTINE printmat_s_2d(matrix, title, iunit, tr)
    REAL(KIND=real_4), DIMENSION(:, :), &
      INTENT(IN)                             :: matrix
    CHARACTER(*), INTENT(IN), OPTIONAL       :: title
    INTEGER                                  :: iunit
    LOGICAL, INTENT(IN), OPTIONAL            :: tr

!   ---------------------------------------------------------------------------

    CALL dbcsr_printmat(RESHAPE(matrix, (/ SIZE(matrix) /)),&
         SIZE(matrix, 1), SIZE(matrix, 2), iunit=iunit, title=title, tr=tr)
  END SUBROUTINE printmat_s_2d


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

  SUBROUTINE printmat_d_2d(matrix, title, iunit, tr)
    REAL(KIND=real_8), DIMENSION(:, :), &
      INTENT(IN)                             :: matrix
    CHARACTER(*), INTENT(IN), OPTIONAL       :: title
    INTEGER                                  :: iunit
    LOGICAL, INTENT(IN), OPTIONAL            :: tr

!   ---------------------------------------------------------------------------

    CALL dbcsr_printmat(RESHAPE(matrix, (/ SIZE(matrix) /)),&
         SIZE(matrix, 1), SIZE(matrix, 2), iunit=iunit, title=title, tr=tr)
  END SUBROUTINE printmat_d_2d


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

  SUBROUTINE printmat_c_2d(matrix, title, iunit, tr)
    COMPLEX(KIND=real_4), DIMENSION(:, :), &
      INTENT(IN)                             :: matrix
    CHARACTER(*), INTENT(IN), OPTIONAL       :: title
    INTEGER                                  :: iunit
    LOGICAL, INTENT(IN), OPTIONAL            :: tr

!   ---------------------------------------------------------------------------

    CALL dbcsr_printmat(RESHAPE(matrix, (/ SIZE(matrix) /)),&
         SIZE(matrix, 1), SIZE(matrix, 2), iunit=iunit, title=title, tr=tr)
  END SUBROUTINE printmat_c_2d

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

  SUBROUTINE printmat_z_2d(matrix, title, iunit, tr)
    COMPLEX(KIND=real_8), DIMENSION(:, :), &
      INTENT(IN)                             :: matrix
    CHARACTER(*), INTENT(IN), OPTIONAL       :: title
    INTEGER                                  :: iunit
    LOGICAL, INTENT(IN), OPTIONAL            :: tr

!   ---------------------------------------------------------------------------

    CALL dbcsr_printmat(RESHAPE(matrix, (/ SIZE(matrix) /)),&
         SIZE(matrix, 1), SIZE(matrix, 2), iunit=iunit, title=title, tr=tr)
  END SUBROUTINE printmat_z_2d


! *****************************************************************************
!> \brief Prints the elements of a matrix.
!> rows and columns is the matrix size. tr specifies whether
!> the elements are stored transposed
! *****************************************************************************
  SUBROUTINE dbcsr_print2dmat(matrix, title, iunit, tr)
    REAL(KIND=real_8), DIMENSION(:, :), &
      INTENT(IN)                             :: matrix
    CHARACTER(*), INTENT(IN), OPTIONAL       :: title
    INTEGER                                  :: iunit
    LOGICAL, INTENT(IN), OPTIONAL            :: tr

    CHARACTER(30)                            :: f
    INTEGER                                  :: cols, m, n, r, rows
    LOGICAL                                  :: t

!   ---------------------------------------------------------------------------

    rows = SIZE(matrix,1)
    cols = SIZE(matrix,2)
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
    WRITE(iunit,*)'Print Matrix 2D:',rows,cols,t,'range',&
         LBOUND(matrix),':',UBOUND(matrix)
    IF(m.GT.10000) m=0
    IF(n.GT.10000) n=0
    IF(m*n.LT.1 .OR. m*n.GT.SIZE(matrix)) RETURN
    WRITE(iunit,FMT="('(',I4,'(F9.4))')")cols
    DO r = 1, rows
       !WRITE(iunit,FMT=f)matrix(r,1:(n-1)*m+r:n)
       IF(.NOT.t) THEN
          WRITE(iunit,FMT=f)matrix(r,:)
       ELSE
          WRITE(iunit,FMT=f)matrix(:,r)
       ENDIF
    ENDDO
  END SUBROUTINE dbcsr_print2dmat


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

  SUBROUTINE print_xfer_timings (timings)
    REAL(kind=dp), DIMENSION(:, :, :, :), &
      INTENT(in)                             :: timings

    INTEGER, PARAMETER                       :: id_bytes = 3, id_recv = 2, &
                                                id_send = 1, id_time = 1, &
                                                id_waittime = 2

    INTEGER                                  :: tick

    DO tick = 1, SIZE(timings,4)-1
       WRITE(*,*)"Tick", tick
       WRITE(*,'(A12,3(1X,EN12.3),1X,3(1X,EN12.3))')"Left recv",&
            timings(1,id_send,:,tick),&
            timings(2,id_send,:,tick)
       WRITE(*,'(A12,3(1X,EN12.3),1X,3(1X,EN12.3))')"Left send",&
            timings(1,id_recv,:,tick),&
            timings(2,id_recv,:,tick)
       WRITE(*,'(A12,3(1X,EN12.3),1X,3(1X,EN12.3))')"Right recv",&
            timings(3,id_send,:,tick),&
            timings(4,id_send,:,tick)
       WRITE(*,'(A12,3(1X,EN12.3),1X,3(1X,EN12.3))')"Right send",&
            timings(3,id_recv,:,tick),&
            timings(4,id_recv,:,tick)
    ENDDO
  END SUBROUTINE print_xfer_timings


END MODULE dbcsr_io
