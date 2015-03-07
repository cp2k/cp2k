!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief represent a full matrix distribued on many processors
!> \par History
!>      3) separated structure object, removed globenv, renamed to full matrix
!>         many changes (fawzi 08.2002)
!> \author Matthias Krack (22.05.2001)
! *****************************************************************************
MODULE cp_fm_types
  USE cp_array_utils_r,                ONLY: cp_2d_r_write
  USE cp_blacs_calls,                  ONLY: cp_blacs_dgebr2d,&
                                             cp_blacs_dgebs2d,&
                                             cp_blacs_gridexit,&
                                             cp_blacs_gridinfo,&
                                             cp_blacs_gridinit
  USE cp_blacs_env,                    ONLY: cp_blacs_env_type
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_fm_struct,                    ONLY: cp_fm_struct_equivalent,&
                                             cp_fm_struct_get,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_retain,&
                                             cp_fm_struct_type,&
                                             cp_fm_struct_write
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: default_path_length,&
                                             dp,&
                                             sp
  USE message_passing,                 ONLY: cp2k_is_parallel,&
                                             mp_bcast,&
                                             mp_max,&
                                             mp_recv,&
                                             mp_send,&
                                             mp_sum
  USE parallel_rng_types,              ONLY: UNIFORM,&
                                             create_rng_stream,&
                                             delete_rng_stream,&
                                             get_rng_stream,&
                                             random_numbers,&
                                             reset_to_next_rng_substream,&
                                             rng_stream_type
  USE string_utilities,                ONLY: compress
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_fm_types'
  LOGICAL, PARAMETER          :: debug_this_module = .TRUE.

  INTEGER, SAVE :: last_fm_id_nr = 0
  INTEGER, PRIVATE :: mm_type=1

  PUBLIC :: cp_fm_type,&
            cp_fm_p_type

  PUBLIC :: cp_fm_add_to_element,&
            cp_fm_create,&
            cp_fm_retain,&
            cp_fm_release,&
            cp_fm_get_info,&
            cp_fm_set_element,&
            cp_fm_get_element,&
            cp_fm_get_diag,&  ! get diagonal
            cp_fm_set_all,&  ! set all elements and diagonal
            cp_fm_set_submatrix,&  ! set a submatrix to given values
            cp_fm_get_submatrix,&  ! get a submatrix of given values
            cp_fm_init_random,&
            cp_fm_maxabsval,&  ! find the maximum absolute value
            cp_fm_maxabsrownorm,&  ! find the maximum of the sum of the abs of the elements of a row
            cp_fm_write,&  ! writes lots of information and data (might need further development)
            cp_fm_to_fm,&  ! copy (parts of) a fm to a fm
            cp_fm_vectorsnorm,&  ! compute the norm of the column-vectors
            cp_fm_to_fm_submat,&  ! copy (parts of) a fm to a fm
            cp_fm_to_fm_triangular,&
            cp_fm_copy_general,&
            cp_fm_write_unformatted,&  ! writes a full matrix to an open unit
            cp_fm_read_unformatted,&  ! reads a full matrix from an open unit
            cp_fm_setup,&  ! allows to set flags for fms
            cp_fm_get_mm_type

  PUBLIC :: cp_fm_indxg2p,&
            cp_fm_indxg2l,&
            cp_fm_indxl2g

  INTERFACE cp_fm_to_fm
    MODULE PROCEDURE cp_fm_to_fm_matrix, &  ! a full matrix
                     cp_fm_to_fm_columns    ! just a number of columns
  END INTERFACE

! *****************************************************************************
!> \brief represent a full matrix
!> \param name the name of the matrix, used for printing
!> \param id_nr identification number that makes this matrix unique
!> \param ref_cont reference count (see doc/ReferenceCounting.html)
!> \param print_count sequential number to make the printed files unique
!> \param matrix_struct structure of this matrix
!> \param local_data array with the data of the matrix (its contents
!>        depend on the matrix type used: in parallel runs it will be
!>        in scalapack format, in sequential, it will simply contain
!>        the matrix)
!> \par History
!>      08.2002 created [fawzi]
!> \author fawzi
! *****************************************************************************
  TYPE cp_fm_type
!    PRIVATE
     CHARACTER(LEN=60) :: name
     INTEGER :: id_nr, ref_count, print_count
     LOGICAL :: use_sp
     TYPE(cp_fm_struct_type), POINTER :: matrix_struct
     REAL(KIND = dp), DIMENSION(:,:), POINTER :: local_data
     REAL(KIND = sp), DIMENSION(:,:), POINTER :: local_data_sp
  END TYPE cp_fm_type

! *****************************************************************************
!> \brief just to build arrays of pointers to matrixes
!> \param matrix the pointer to the matrix
!> \par History
!>      08.2002 created [fawzi]
!> \author fawzi
! *****************************************************************************
  TYPE cp_fm_p_type
     TYPE(cp_fm_type), POINTER :: matrix
  END TYPE cp_fm_p_type

CONTAINS

! *****************************************************************************
!> \brief creates a new full matrix with the given structure
!> \param matrix the matrix to be created
!> \param matrix_struct the structure of matrix
!> \param name ...
!> \param use_sp ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      preferred allocation routine
! *****************************************************************************
  SUBROUTINE cp_fm_create(matrix,matrix_struct,name,use_sp,error)
    TYPE(cp_fm_type), POINTER                :: matrix
    TYPE(cp_fm_struct_type), POINTER         :: matrix_struct
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: name
    LOGICAL, INTENT(in), OPTIONAL            :: use_sp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ncol_local, npcol, &
                                                nprow, nrow_local, stat
    LOGICAL                                  :: failure
    TYPE(cp_blacs_env_type), POINTER         :: context

    failure=.FALSE.

    CALL timeset(routineN,handle)

#if defined(__parallel) && ! defined(__SCALAPACK)
     CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
     routineP,"full matrixes need scalapack for parallel runs "//&
CPSourceFileRef,&
     error)
#endif

    CPPrecondition(ASSOCIATED(matrix_struct),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       ALLOCATE(matrix,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT. failure) THEN
       context => matrix_struct%context
       matrix%matrix_struct => matrix_struct
       CALL cp_fm_struct_retain(matrix%matrix_struct,error=error)
       last_fm_id_nr=last_fm_id_nr+1
       matrix%id_nr=last_fm_id_nr
       matrix%ref_count=1
       matrix%print_count=0

       matrix%use_sp = .FALSE.
       IF(PRESENT(use_sp)) matrix%use_sp = use_sp

       nprow=context%num_pe(1)
       npcol=context%num_pe(2)
       NULLIFY(matrix%local_data)
       NULLIFY(matrix%local_data_sp)

       ! OK, we allocate here at least a 1 x 1 matrix
       ! this must (and is) compatible with the descinit call
       ! in cp_fm_struct
       nrow_local=matrix_struct%local_leading_dimension
       ncol_local=MAX(1,matrix_struct%ncol_locals(context%mepos(2)))
       IF(matrix%use_sp) THEN
          ALLOCATE(matrix%local_data_sp(nrow_local,ncol_local),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ELSE
          ALLOCATE(matrix%local_data(nrow_local,ncol_local),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ENDIF

       ! JVDV we should remove this, as it is up to the user to zero afterwards
       IF(matrix%use_sp) THEN
          CALL scopy(nrow_local*ncol_local,0.0_sp,0,matrix%local_data_sp,1)
       ELSE
          CALL dcopy(nrow_local*ncol_local,0.0_dp,0,matrix%local_data,1)
       ENDIF

       IF (PRESENT(name)) THEN
          matrix%name=name
       ELSE
          matrix%name='full matrix'//cp_to_string(matrix%id_nr)
       END IF
    END IF

    CALL timestop(handle)

  END SUBROUTINE cp_fm_create

! *****************************************************************************
!> \brief retains a full matrix
!> \param matrix the matrix to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE cp_fm_retain(matrix,error)
    TYPE(cp_fm_type), POINTER                :: matrix
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(matrix),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPrecondition(matrix%ref_count>0,cp_failure_level,routineP,error,failure)
       matrix%ref_count=matrix%ref_count+1
    END IF

  END SUBROUTINE cp_fm_retain

! *****************************************************************************
!> \brief releases a full matrix
!> \param matrix the matrix to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE cp_fm_release(matrix,error)
    TYPE(cp_fm_type), POINTER                :: matrix
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CALL timeset(routineN,handle)

    IF (ASSOCIATED(matrix)) THEN
       CPPrecondition(matrix%ref_count>0,cp_failure_level,routineP,error,failure)
       matrix%ref_count=matrix%ref_count-1
       IF (matrix%ref_count<1) THEN
          IF (ASSOCIATED(matrix%local_data)) THEN
             DEALLOCATE(matrix%local_data,stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(matrix%local_data_sp)) THEN
             DEALLOCATE(matrix%local_data_sp,stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
          matrix%name=""
          CALL cp_fm_struct_release(matrix%matrix_struct,error=error)
          DEALLOCATE(matrix,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
    END IF
    NULLIFY(matrix)

    CALL timestop(handle)

  END SUBROUTINE cp_fm_release

! *****************************************************************************
!> \brief fills a matrix with random numbers
!> \param matrix : to be initialized
!> \param ncol : numbers of cols to fill
!> \param start_col : starting at coll number
!> \param error ...
!> \author Joost VandeVondele
!> \note
!>      the value of a_ij is independent of the number of cpus
! *****************************************************************************
  SUBROUTINE cp_fm_init_random(matrix,ncol,start_col,error)
    TYPE(cp_fm_type), POINTER                :: matrix
    INTEGER, INTENT(IN), OPTIONAL            :: ncol, start_col
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_init_random', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, icol_global, icol_local, irow_local, my_ncol, &
      my_start_col, ncol_global, ncol_local, nrow_global, nrow_local
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    LOGICAL                                  :: failure
    LOGICAL, SAVE                            :: FIRST = .TRUE.
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: buff
    REAL(KIND=dp), DIMENSION(3, 2), SAVE     :: seed
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: local_data
    TYPE(rng_stream_type), POINTER           :: rng

    failure = .FALSE.

    CALL timeset(routineN,handle)

    ! first setup a parallel rng
    IF (FIRST) THEN
       seed(:,:)=RESHAPE((/1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp/),(/3,2/))
       FIRST=.FALSE.
    ENDIF
    ! guarantee same seed over all tasks
    CALL mp_bcast(seed,0,matrix%matrix_struct%para_env%group)

    NULLIFY(rng)
    CALL create_rng_stream(rng,"cp_fm_init_random_stream",distribution_type=UNIFORM, &
                           extended_precision=.TRUE.,seed=seed,error=error)

    CPPrecondition(.NOT.matrix%use_sp,cp_failure_level,routineP,error,failure)

    CALL cp_fm_get_info(matrix,nrow_global=nrow_global,ncol_global=ncol_global, &
                               nrow_local=nrow_local,ncol_local=ncol_local,&
                               local_data=local_data,&
                               row_indices=row_indices, col_indices=col_indices, error=error)

    my_start_col = 1
    IF (PRESENT(start_col)) my_start_col=start_col
    my_ncol = matrix%matrix_struct%ncol_global
    IF (PRESENT(ncol)) my_ncol=ncol

    CALL cp_assert(ncol_global>=(my_start_col+my_ncol-1),&
         cp_failure_level,cp_assertion_failed,routineP,&
         "ncol_global>=(my_start_col+my_ncol-1)",&
         error,failure)

    ALLOCATE(buff(nrow_global))

    ! each global row has its own substream, in order to reach the stream for the local col,
    ! we just reset to the next substream
    ! following this, we fill the full buff with random numbers, and pick those we need
    icol_global=0
    DO icol_local=1,ncol_local
       CPPrecondition(col_indices(icol_local)>icol_global,cp_failure_level,routineP,error,failure)
       DO
         CALL reset_to_next_rng_substream(rng,error=error)
         icol_global=icol_global+1
         IF (icol_global==col_indices(icol_local)) EXIT
       ENDDO
       CALL random_numbers(buff,rng,error)
       DO irow_local=1,nrow_local
          local_data(irow_local,icol_local)=buff(row_indices(irow_local))
       ENDDO
    ENDDO

    DEALLOCATE(buff)

    ! store seed before deletion (unclear if this is the proper seed)
    CALL get_rng_stream(rng,ig=seed,error=error)
    CALL delete_rng_stream(rng,error)

    CALL timestop(handle)

  END SUBROUTINE cp_fm_init_random

! *****************************************************************************
!> \brief set all elements of a matrix to the same value,
!>      and optionally the diagonal to a different one
!> \param matrix input matrix
!> \param alpha scalar used to set all elements of the matrix
!> \param beta scalar used to set diagonal of the matrix
!> \param error ...
!> \note
!>      can be used to zero a matrix
!>      can be used to create a unit matrix (I-matrix) alpha=0.0_dp beta=1.0_dp
! *****************************************************************************
  SUBROUTINE cp_fm_set_all(matrix,alpha,beta,error)

    TYPE(cp_fm_type), POINTER                :: matrix
    REAL(KIND=dp), INTENT(IN)                :: alpha
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: beta
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_set_all', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, matrix_size, n
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)

    failure = .FALSE.

    IF(matrix%use_sp) THEN
       matrix_size = SIZE(matrix%local_data_sp,1)*SIZE(matrix%local_data_sp,2)
       CALL scopy (matrix_size, REAL(alpha,sp), 0, matrix%local_data_sp, 1)
    ELSE
       matrix_size = SIZE(matrix%local_data,1)*SIZE(matrix%local_data,2)
       CALL dcopy (matrix_size, alpha, 0, matrix%local_data, 1)
    ENDIF

    IF (PRESENT(beta)) THEN
       CPPrecondition(.NOT.matrix%use_sp,cp_failure_level,routineP,error,failure)
      n = MIN(matrix%matrix_struct%nrow_global,matrix%matrix_struct%ncol_global)
      DO i=1,n
         CALL cp_fm_set_element(matrix,i,i,beta,error=error)
      END DO
    END IF

    CALL timestop(handle)

  END SUBROUTINE cp_fm_set_all

! *****************************************************************************
!> \brief returns the diagonal elements of a fm
!> \param matrix ...
!> \param diag ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_fm_get_diag(matrix,diag,error)

    IMPLICIT NONE

    ! arguments
    TYPE(cp_fm_type), POINTER                  :: matrix
    REAL(KIND = dp), DIMENSION(:), INTENT(OUT) :: diag
    TYPE(cp_error_type), INTENT(inout)         :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_get_diag', &
      routineP = moduleN//':'//routineN

    ! locals
    INTEGER :: i,nrow_global
    LOGICAL :: failure

#if defined(__SCALAPACK)
    INTEGER, DIMENSION(9) :: desca
    TYPE(cp_blacs_env_type), POINTER :: context
    INTEGER :: icol_local,ipcol,iprow,irow_local,mypcol,myprow,npcol,&
               nprow
    REAL(KIND = dp), DIMENSION(:,:), POINTER :: a
    REAL(KIND = sp), DIMENSION(:,:), POINTER :: a_sp
#endif

    failure = .FALSE.

    CALL cp_fm_get_info(matrix,nrow_global=nrow_global,error=error)

#if defined(__SCALAPACK)
    diag=0.0_dp
    context => matrix%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)

    a => matrix%local_data
    a_sp => matrix%local_data_sp
    desca(:) = matrix%matrix_struct%descriptor(:)

    DO i=1,nrow_global
       CALL infog2l(i,i,desca,nprow,npcol,myprow,mypcol,&
                    irow_local,icol_local,iprow,ipcol)
       IF ((iprow == myprow).AND.(ipcol == mypcol)) THEN
          IF(matrix%use_sp) THEN
             diag(i) = REAL(a_sp(irow_local,icol_local),dp)
          ELSE
             diag(i) = a(irow_local,icol_local)
          ENDIF
       ENDIF
    END DO
#else
    DO i=1, nrow_global
       IF(matrix%use_sp) THEN
          diag(i)=REAL(matrix%local_data_sp(i,i),dp)
       ELSE
          diag(i)=matrix%local_data(i,i)
       ENDIF
    ENDDO
#endif
    CALL mp_sum(diag,matrix%matrix_struct%para_env%group)

  END SUBROUTINE cp_fm_get_diag

! *****************************************************************************
!> \brief returns an element of a fm
!>      this value is valid on every cpu
!>      using this call is expensive
!> \param matrix the matrix to read
!> \param irow_global the row
!> \param icol_global the col
!> \param alpha the value of matrix(irow_global, icol_global)
!> \param local true if the element is on this cpu, false otherwise
!> \note
!>      - modified semantics. now this function always returns the value
!>        previously the value was zero on cpus that didn't own the relevant
!>        part of the matrix (Joost VandeVondele, May 2003)
!>      - usage of the function should be avoided, as it is likely to rather slow
!>        using row_indices/col_indices/local_data + some smart scheme normally
!>        yields a real parallel code
! *****************************************************************************
  SUBROUTINE cp_fm_get_element(matrix,irow_global,icol_global,alpha,local)

    IMPLICIT NONE

    ! arguments
    TYPE(cp_fm_type), POINTER          :: matrix
    REAL(KIND = dp), INTENT(OUT)                     :: alpha
    INTEGER, INTENT(IN)                       :: icol_global,&
                                                 irow_global
    LOGICAL, INTENT(OUT), OPTIONAL            :: local
    LOGICAL :: failure

    ! locals
#if defined(__SCALAPACK)
    INTEGER, DIMENSION(9) :: desca
    TYPE(cp_blacs_env_type), POINTER :: context
    INTEGER :: icol_local,ipcol,iprow,irow_local,mypcol,myprow,npcol,&
               nprow
    REAL(KIND = dp), DIMENSION(:,:), POINTER :: a
#endif

    failure = .FALSE.

#if defined(__SCALAPACK)
    context => matrix%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)

    a => matrix%local_data
    desca(:) = matrix%matrix_struct%descriptor(:)

    CALL infog2l(irow_global,icol_global,desca,nprow,npcol,myprow,mypcol,&
                 irow_local,icol_local,iprow,ipcol)

    IF ((iprow == myprow).AND.(ipcol == mypcol)) THEN
      alpha = a(irow_local,icol_local)
      CALL cp_blacs_dgebs2d(context%group, 'All', ' ', 1, 1, alpha , 1 )
      IF (PRESENT(local)) local=.TRUE.
    ELSE
      CALL cp_blacs_dgebr2d(context%group, 'All', ' ', 1, 1, alpha , 1 , iprow, ipcol )
      IF (PRESENT(local)) local=.FALSE.
    END IF

#else
    IF (PRESENT(local)) local=.TRUE.
    alpha = matrix%local_data(irow_global,icol_global)
#endif

  END SUBROUTINE cp_fm_get_element

! *****************************************************************************
!> \brief sets an element of a matrix
!> \param matrix ...
!> \param irow_global ...
!> \param icol_global ...
!> \param alpha ...
!> \param error ...
!> \note
!>      we expect all cpus to have the same arguments in the call to this function
!>      (otherwise one should use local_data tricks)
! *****************************************************************************
  SUBROUTINE cp_fm_set_element(matrix,irow_global,icol_global,alpha,error)
    TYPE(cp_fm_type), POINTER                :: matrix
    INTEGER, INTENT(IN)                      :: irow_global, icol_global
    REAL(KIND=dp), INTENT(IN)                :: alpha
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_set_element', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: mypcol, myprow, npcol, nprow
    LOGICAL                                  :: failure
    TYPE(cp_blacs_env_type), POINTER         :: context
#if defined(__SCALAPACK)
    INTEGER                                  :: icol_local, ipcol, iprow, &
                                                irow_local
    INTEGER, DIMENSION(9)                    :: desca
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a
#endif

    failure = .FALSE.

    context => matrix%matrix_struct%context
    myprow = context%mepos(1)
    mypcol = context%mepos(2)
    nprow = context%num_pe(1)
    npcol = context%num_pe(2)

    CPPrecondition(.NOT.matrix%use_sp,cp_failure_level,routineP,error,failure)

#if defined(__SCALAPACK)

    a => matrix%local_data

    desca(:) = matrix%matrix_struct%descriptor(:)

    CALL infog2l(irow_global,icol_global,desca,nprow,npcol,myprow,mypcol,&
         irow_local,icol_local,iprow,ipcol)

    IF ((iprow == myprow).AND.(ipcol == mypcol)) THEN
       a(irow_local,icol_local) = alpha
    END IF

#else

    matrix%local_data(irow_global,icol_global) = alpha

#endif
  END SUBROUTINE cp_fm_set_element

! *****************************************************************************
!> \brief sets a submatrix of a full matrix
!>       fm(start_row:start_row+n_rows,start_col:start_col+n_cols)
!>       = alpha*op(new_values)(1:n_rows,1:n_cols)+ beta
!>       * fm(start_row:start_row+n_rows,start_col:start_col+n_cols)
!> \param fm the full to change
!> \param new_values a replicated full matrix with the new values
!> \param start_row the starting row of b_matrix (defaults to 1)
!> \param start_col the starting col of b_matrix (defaults to 1)
!> \param n_rows the number of row to change in b (defaults to
!>        size(op(new_values),1))
!> \param n_cols the number of columns to change in b (defaults to
!>        size(op(new_values),2))
!> \param alpha rescaling factor for the new values (defaults to 1.0)
!> \param beta rescaling factor for the old values (defaults to 0.0)
!> \param transpose if new_values should be transposed: if true
!>        op(new_values)=new_values^T, else op(new_values)=new_values
!>        (defaults to false)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      07.2002 created borrowing from Joost's blacs_replicated_copy [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      optimized for full column updates and alpha=1.0, beta=0.0
!>      the new_values need to be valid on all cpus
! *****************************************************************************
  SUBROUTINE cp_fm_set_submatrix(fm,new_values,start_row,&
       start_col, n_rows, n_cols, alpha, beta, transpose, error)
    TYPE(cp_fm_type), POINTER                :: fm
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(in)                             :: new_values
    INTEGER, INTENT(in), OPTIONAL            :: start_row, start_col, n_rows, &
                                                n_cols
    REAL(KIND=dp), INTENT(in), OPTIONAL      :: alpha, beta
    LOGICAL, INTENT(in), OPTIONAL            :: transpose
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_set_submatrix', &
      routineP = moduleN//':'//routineN

    INTEGER :: i, i0, j, j0, ncol, ncol_block, ncol_global, ncol_local, nrow, &
      nrow_block, nrow_global, nrow_local, this_col, this_row
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    LOGICAL                                  :: failure, tr_a
    REAL(KIND=dp)                            :: al, be
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: full_block

    al=1.0_dp; be=0.0_dp; i0=1; j0=1; tr_a=.FALSE.
    ! can be called too many times, making it a bit useless
    ! CALL timeset(routineN//','//moduleN,handle)

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(fm),cp_failure_level,routineP,error,failure)
    IF (.not.failure) THEN
       CPPrecondition(fm%ref_count>0,cp_failure_level,routineP,error,failure)
    END IF

    CPPrecondition(.NOT.fm%use_sp,cp_failure_level,routineP,error,failure)

    IF (.not.failure) THEN
       IF (PRESENT(alpha)) al=alpha
       IF (PRESENT(beta)) be=beta
       IF (PRESENT(start_row)) i0=start_row
       IF (PRESENT(start_col)) j0=start_col
       IF (PRESENT(transpose)) tr_a=transpose
       IF (tr_a) THEN
          nrow=SIZE(new_values,2)
          ncol=SIZE(new_values,1)
       ELSE
          nrow=SIZE(new_values,1)
          ncol=SIZE(new_values,2)
       END IF
       IF (PRESENT(n_rows)) nrow=n_rows
       IF (PRESENT(n_cols)) ncol=n_cols

       full_block => fm%local_data

       CALL cp_fm_get_info(matrix=fm,&
            nrow_global=nrow_global,ncol_global=ncol_global,&
            nrow_block =nrow_block ,ncol_block =ncol_block ,&
            nrow_local =nrow_local ,ncol_local =ncol_local ,&
            row_indices=row_indices,col_indices=col_indices,error=error)

       IF (al==1.0.AND.be==0.0) THEN
          DO j=1,ncol_local
             this_col=col_indices(j)-j0+1
             IF (this_col.GE.1 .AND. this_col.LE.ncol) THEN
                IF (tr_a) THEN
                   IF (i0==1.AND.nrow_global==nrow) THEN
                      DO i=1,nrow_local
                         full_block(i,j)=new_values(this_col,row_indices(i))
                      END DO
                   ELSE
                      DO i=1,nrow_local
                         this_row=row_indices(i)-i0+1
                         IF (this_row>=1 .AND. this_row<=nrow) THEN
                            full_block(i,j)=new_values(this_col,this_row)
                         END IF
                      END DO
                   END IF
                ELSE
                   IF (i0==1.AND.nrow_global==nrow) THEN
                      DO i=1,nrow_local
                         full_block(i,j)=new_values(row_indices(i),this_col)
                      END DO
                   ELSE
                      DO i=1,nrow_local
                         this_row=row_indices(i)-i0+1
                         IF (this_row>=1 .AND. this_row<=nrow) THEN
                            full_block(i,j)=new_values(this_row,this_col)
                         END IF
                      END DO
                   END IF
                END IF
             END IF
          END DO
       ELSE
          DO j=1,ncol_local
             this_col=col_indices(j)-j0+1
             IF (this_col.GE.1 .AND. this_col.LE.ncol) THEN
                IF (tr_a) THEN
                   DO i=1,nrow_local
                      this_row=row_indices(i)-i0+1
                      IF (this_row>=1 .AND. this_row<=nrow) THEN
                         full_block(i,j)=al*new_values(this_col,this_row)+&
                              be*full_block(i,j)
                      END IF
                   END DO
                ELSE
                   DO i=1,nrow_local
                      this_row=row_indices(i)-i0+1
                      IF (this_row>=1 .AND. this_row<=nrow) THEN
                         full_block(i,j)=al*new_values(this_row,this_col)+&
                              be*full_block(i,j)
                      END IF
                   END DO
                END IF
             END IF
          END DO
       END IF
    END IF

    ! CALL timestop(handle)

  END SUBROUTINE cp_fm_set_submatrix

! *****************************************************************************
!> \brief gets a submatrix of a full matrix
!>       op(target_m)(1:n_rows,1:n_cols)
!>       =fm(start_row:start_row+n_rows,start_col:start_col+n_cols)
!>      target_m is replicated on all cpus
!>      using this call is expensive
!> \param fm the full you want to get the info from
!> \param target_m a replicated full matrix that will contain the result
!> \param start_row the starting row of b_matrix (defaults to 1)
!> \param start_col the starting col of b_matrix (defaults to 1)
!> \param n_rows the number of row to change in b (defaults to
!>        size(op(new_values),1))
!> \param n_cols the number of columns to change in b (defaults to
!>        size(op(new_values),2))
!> \param transpose if target_m should be transposed: if true
!>        op(target_m)=target_m^T, else op(target_m)=target_m
!>        (defaults to false)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      07.2002 created borrowing from Joost's blacs_replicated_copy [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      optimized for full column updates. Zeros out a little too much
!>      of target_m
!>      the target_m is replicated and valid on all cpus
! *****************************************************************************
  SUBROUTINE cp_fm_get_submatrix(fm,target_m, start_row,&
       start_col, n_rows, n_cols, transpose, error)
    TYPE(cp_fm_type), POINTER                :: fm
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(out)                            :: target_m
    INTEGER, INTENT(in), OPTIONAL            :: start_row, start_col, n_rows, &
                                                n_cols
    LOGICAL, INTENT(in), OPTIONAL            :: transpose
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_get_submatrix', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, i0, j, j0, ncol, ncol_global, ncol_local, nrow, &
      nrow_global, nrow_local, this_col, this_row
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    LOGICAL                                  :: failure, tr_a
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: full_block
    TYPE(cp_para_env_type), POINTER          :: para_env

    CALL timeset(routineN,handle)

    i0=1; j0=1; tr_a=.FALSE.
    failure=.FALSE.

    CPPrecondition(ASSOCIATED(fm),cp_failure_level,routineP,error,failure)
    IF (.not.failure) THEN
       CPPrecondition(fm%ref_count>0,cp_failure_level,routineP,error,failure)
    END IF

    CPPrecondition(.NOT.fm%use_sp,cp_failure_level,routineP,error,failure)

    IF (.not.failure) THEN
       IF (PRESENT(start_row)) i0=start_row
       IF (PRESENT(start_col)) j0=start_col
       IF (PRESENT(transpose)) tr_a=transpose
       IF (tr_a) THEN
          nrow=SIZE(target_m,2)
          ncol=SIZE(target_m,1)
       ELSE
          nrow=SIZE(target_m,1)
          ncol=SIZE(target_m,2)
       END IF
       IF (PRESENT(n_rows)) nrow=n_rows
       IF (PRESENT(n_cols)) ncol=n_cols

       para_env => fm%matrix_struct%para_env

       full_block => fm%local_data
#if defined(__SCALAPACK)
       ! zero-out whole target_m
       IF (SIZE(target_m,1)*SIZE(target_m,2)/=0) THEN
          CALL dcopy(SIZE(target_m,1)*SIZE(target_m,2),0.0_dp,0,target_m(1,1),1)
       END IF
#endif

       CALL cp_fm_get_info(matrix=fm, &
            nrow_global=nrow_global,ncol_global=ncol_global, &
            nrow_local =nrow_local ,ncol_local =ncol_local , &
            row_indices=row_indices,col_indices=col_indices,error=error)

       DO j=1,ncol_local
          this_col=col_indices(j)-j0+1
          IF (this_col.GE.1 .AND. this_col.LE.ncol) THEN
             IF (tr_a) THEN
                IF (i0==1.AND.nrow_global==nrow) THEN
                   DO i=1,nrow_local
                      target_m(this_col,row_indices(i))=full_block(i,j)
                   END DO
                ELSE
                   DO i=1,nrow_local
                      this_row=row_indices(i)-i0+1
                      IF (this_row>=1 .AND. this_row<=nrow) THEN
                         target_m(this_col,this_row)=full_block(i,j)
                      END IF
                   END DO
                END IF
             ELSE
                IF (i0==1.AND.nrow_global==nrow) THEN
                   DO i=1,nrow_local
                      target_m(row_indices(i),this_col)=full_block(i,j)
                   END DO
                ELSE
                   DO i=1,nrow_local
                      this_row=row_indices(i)-i0+1
                      IF (this_row>=1 .AND. this_row<=nrow) THEN
                         target_m(this_row,this_col)=full_block(i,j)
                      END IF
                   END DO
                END IF
             END IF
          END IF
       END DO

       CALL mp_sum(target_m,para_env%group)
    END IF

    CALL timestop(handle)

  END SUBROUTINE cp_fm_get_submatrix

! *****************************************************************************
!> \brief returns all kind of information about the full matrix
!> \param matrix ...
!> \param name ...
!> \param nrow_global ...
!> \param ncol_global ...
!> \param nrow_block ...
!> \param ncol_block ...
!> \param nrow_local ...
!> \param ncol_local ...
!> \param row_indices ...
!> \param col_indices ...
!> \param local_data ...
!> \param context ...
!> \param nrow_locals ...
!> \param ncol_locals ...
!> \param matrix_struct ...
!> \param para_env ...
!> \param error ...
!> \note
!>       see also cp_fm_struct for explaination
!>       - nrow_local, ncol_local, row_indices, col_indices, local_data are hooks for efficient
!>         access to the local blacs block
! *****************************************************************************
  SUBROUTINE cp_fm_get_info(matrix,name,nrow_global,ncol_global,&
       nrow_block,ncol_block,nrow_local,ncol_local,&
       row_indices,col_indices,local_data,context,&
       nrow_locals, ncol_locals, matrix_struct,para_env,error)

    TYPE(cp_fm_type), POINTER                  :: matrix
    CHARACTER(LEN=*), OPTIONAL, INTENT(OUT)    :: name
    INTEGER, OPTIONAL, INTENT(OUT)             :: ncol_block,ncol_global,&
                                                  nrow_block,nrow_global,&
                                                  nrow_local,ncol_local
    INTEGER, OPTIONAL, DIMENSION(:), POINTER   :: row_indices,col_indices,&
                                                  nrow_locals, ncol_locals
    TYPE(cp_para_env_type), POINTER, OPTIONAL  :: para_env
    TYPE(cp_blacs_env_type), POINTER, OPTIONAL :: context
    TYPE(cp_fm_struct_type),POINTER,OPTIONAL   :: matrix_struct
    TYPE(cp_error_type),INTENT(inout):: error
    REAL(KIND = dp), DIMENSION(:,:),POINTER, OPTIONAL :: local_data

    CHARACTER(len=*), PARAMETER :: routineN='cp_fm_get_info',&
         routineP=moduleN//':'//routineN
#if defined(__SCALAPACK)
    INTEGER , EXTERNAL :: indxl2g
#endif

    IF (PRESENT(name)) name = matrix%name
    IF (PRESENT(context)) context => matrix%matrix_struct%context
    IF (PRESENT(matrix_struct)) matrix_struct => matrix%matrix_struct
    IF (PRESENT(local_data)) local_data => matrix%local_data ! not hiding things anymore :-(
    IF (PRESENT(row_indices)) THEN
       CALL cp_fm_struct_get(matrix%matrix_struct, row_indices=row_indices,&
            error=error)
    ENDIF
    IF (PRESENT(col_indices)) THEN
       CALL cp_fm_struct_get(matrix%matrix_struct, col_indices=col_indices,&
            error=error)
    ENDIF
    IF (PRESENT(nrow_locals)) THEN
       nrow_locals => matrix%matrix_struct%nrow_locals
    END IF
    IF (PRESENT(ncol_locals)) THEN
       ncol_locals => matrix%matrix_struct%ncol_locals
    END IF

    CALL cp_fm_struct_get(matrix%matrix_struct, nrow_local=nrow_local,&
          ncol_local=ncol_local, nrow_global=nrow_global,&
          ncol_global=ncol_global, nrow_block=nrow_block, &
          ncol_block=ncol_block, error=error)

    IF (PRESENT(para_env)) para_env => matrix%matrix_struct%para_env

  END SUBROUTINE cp_fm_get_info

! *****************************************************************************
!> \brief find the maximum absolute value of the matrix element
!>      maxval(abs(matrix))
!> \param matrix ...
!> \param a_max ...
!> \param ir_max ...
!> \param ic_max ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_fm_maxabsval(matrix,a_max,ir_max,ic_max,error)
    TYPE(cp_fm_type), POINTER                :: matrix
    REAL(KIND=dp), INTENT(OUT)               :: a_max
    INTEGER, INTENT(OUT), OPTIONAL           :: ir_max, ic_max
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_maxabsval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ic_max_local, &
                                                ir_max_local, istat, j, &
                                                mepos, ncol_local, &
                                                nrow_local, num_pe
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ic_max_vec, ir_max_vec
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    LOGICAL                                  :: failure
    REAL(dp)                                 :: my_max
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: a_max_vec
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: my_block
    REAL(KIND=sp), DIMENSION(:, :), POINTER  :: my_block_sp

    CALL timeset(routineN,handle)

    failure = .FALSE.
    my_block => matrix%local_data
    my_block_sp => matrix%local_data_sp

    CALL cp_fm_get_info(matrix, nrow_local=nrow_local, ncol_local=ncol_local,&
                        row_indices=row_indices,col_indices=col_indices,error=error)

    IF(matrix%use_sp) THEN
       a_max = REAL(MAXVAL(ABS(my_block_sp(1:nrow_local,1:ncol_local))),dp)
    ELSE
       a_max = MAXVAL(ABS(my_block(1:nrow_local,1:ncol_local)))
    ENDIF

    IF(PRESENT(ir_max)) THEN
      num_pe = matrix%matrix_struct%para_env%num_pe
      mepos = matrix%matrix_struct%para_env%mepos
      ALLOCATE(ir_max_vec(0:num_pe-1),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ir_max_vec(0:num_pe-1) = 0
      ALLOCATE(ic_max_vec(0:num_pe-1),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ic_max_vec(0:num_pe-1) = 0
      ALLOCATE(a_max_vec(0:num_pe-1),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      a_max_vec(0:num_pe-1) = 0.0_dp
      my_max = 0.0_dp

      IF((ncol_local > 0) .AND. (nrow_local > 0)) THEN
        DO i = 1,ncol_local
          DO j = 1,nrow_local
             IF(matrix%use_sp) THEN
                IF( ABS(my_block_sp(j,i)) > my_max) THEN
                   my_max = REAL(my_block_sp(j,i),dp)
                   ir_max_local = j
                   ic_max_local = i
                ENDIF
             ELSE
                IF( ABS(my_block(j,i)) > my_max) THEN
                   my_max =  my_block(j,i)
                   ir_max_local = j
                   ic_max_local = i
                ENDIF
            END IF
          END DO
        END DO

        a_max_vec(mepos) = my_max
        ir_max_vec(mepos) = row_indices(ir_max_local)
        ic_max_vec(mepos) = col_indices(ic_max_local)

      END IF

      CALL mp_sum(a_max_vec,matrix%matrix_struct%para_env%group)
      CALL mp_sum(ir_max_vec,matrix%matrix_struct%para_env%group)
      CALL mp_sum(ic_max_vec,matrix%matrix_struct%para_env%group)

      my_max = 0.0_dp
      DO i = 0,num_pe-1
        IF(a_max_vec(i)>my_max) THEN
           ir_max = ir_max_vec(i)
           ic_max = ic_max_vec(i)
        END IF
      END DO

      DEALLOCATE(ir_max_vec,ic_max_vec,a_max_vec, STAT = istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      CPPostcondition(ic_max>0,cp_failure_level,routineP,error,failure)
      CPPostcondition(ir_max>0,cp_failure_level,routineP,error,failure)

    END IF

    CALL mp_max(a_max,matrix%matrix_struct%para_env%group)

    CALL timestop(handle)

  END SUBROUTINE cp_fm_maxabsval

! *****************************************************************************
!> \brief find the maximum over the rows of the sum of the absolute values of the elements of a given row
!>      = || A ||_infinity
!> \param matrix ...
!> \param a_max ...
!> \param error ...
!> \note
!>      for a real symmetric matrix it holds that || A ||_2 = |lambda_max| < || A ||_infinity
!>      Hence this can be used to estimate an upper bound for the eigenvalues of a matrix
!>      http://mathworld.wolfram.com/MatrixNorm.html
!>      (but the bound is not so tight in the general case)
! *****************************************************************************
  SUBROUTINE cp_fm_maxabsrownorm(matrix, a_max,error)
    TYPE(cp_fm_type), POINTER                :: matrix
    REAL(KIND=dp), INTENT(OUT)               :: a_max
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_maxabsrownorm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, ncol_local, &
                                                nrow_global, nrow_local
    INTEGER, DIMENSION(:), POINTER           :: row_indices
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: values
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: my_block

    CALL timeset(routineN,handle)

    my_block => matrix%local_data

    CPPrecondition(.NOT.matrix%use_sp,cp_failure_level,routineP,error,failure)

    CALL cp_fm_get_info(matrix, row_indices=row_indices, nrow_global=nrow_global, &
                                nrow_local=nrow_local, ncol_local=ncol_local,error=error)

    ! the efficiency could be improved by making use of the row-col distribution of scalapack
    ALLOCATE(values(nrow_global))
    values=0.0_dp
    DO j=1,ncol_local
       DO i=1,nrow_local
          values(row_indices(i))=values(row_indices(i))+ABS(my_block(i,j))
       ENDDO
    ENDDO
    CALL mp_sum(values,matrix%matrix_struct%para_env%group)
    a_max = MAXVAL(values)
    DEALLOCATE(values)

    CALL timestop(handle)
  END SUBROUTINE

! *****************************************************************************
!> \brief find the inorm of each column norm_{j}= sqrt( \sum_{i} A_{ij}*A_{ij} )
!> \param matrix ...
!> \param norm_array ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_fm_vectorsnorm(matrix, norm_array, error)
    TYPE(cp_fm_type), POINTER                :: matrix
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: norm_array
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_vectorsnorm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, ncol_global, &
                                                ncol_local, nrow_local
    INTEGER, DIMENSION(:), POINTER           :: col_indices
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: my_block

    failure = .FALSE.

    CALL timeset(routineN,handle)

    my_block => matrix%local_data

    CPPrecondition(.NOT.matrix%use_sp,cp_failure_level,routineP,error,failure)

    CALL cp_fm_get_info(matrix, col_indices=col_indices, ncol_global=ncol_global, &
                                nrow_local=nrow_local, ncol_local=ncol_local, error=error)

    ! the efficiency could be improved by making use of the row-col distribution of scalapack
    norm_array=0.0_dp
    DO j=1,ncol_local
       DO i=1,nrow_local
          norm_array(col_indices(j))=norm_array(col_indices(j))+my_block(i,j)*my_block(i,j)
       ENDDO
    ENDDO
    CALL mp_sum(norm_array,matrix%matrix_struct%para_env%group)
    norm_array = SQRT(norm_array)

    CALL timestop(handle)
  END SUBROUTINE

! *****************************************************************************
!> \brief writes out the matrix (if long description, might need some refactoring
!> \param matrix the matrix to output
!> \param unit_nr the unit where to output
!> \param long_description if a long description should be written (is
!>        needed to see the content of the matrix)
!>        (defaults to false)
!> \param local if the unit is a local unit or a global unit
!>        (defaults to false, i.e. global)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      uniquing_nr should be smartened: don't add if no file exist, else
!>      increase nr until no file exists
!>      routine to help debugging
! *****************************************************************************
   SUBROUTINE cp_fm_write(matrix, unit_nr, long_description, local, error)
    TYPE(cp_fm_type), POINTER                :: matrix
    INTEGER, INTENT(in)                      :: unit_nr
    LOGICAL, INTENT(in), OPTIONAL            :: long_description, local
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_write', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_path_length)       :: filename
    INTEGER                                  :: data_unit, desc_unit, iostat
    LOGICAL                                  :: exists, failure, loc, &
                                                my_long_description
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env

     failure=.FALSE.; loc=.FALSE.; my_long_description=.FALSE.
     filename=' '

     IF (PRESENT(local)) loc=local
     IF (PRESENT(long_description)) my_long_description=long_description
     logger=>cp_error_get_logger(error)
     para_env=>logger%para_env
     IF (ASSOCIATED(matrix)) THEN
        data_unit=unit_nr
        desc_unit=unit_nr
        IF (.not.loc.and.my_long_description) THEN
           matrix%print_count=matrix%print_count+1
           IF (para_env%mepos==para_env%source) THEN
              filename=TRIM(matrix%name)//'-'//&
                   ADJUSTL(cp_to_string(matrix%print_count))
              CALL compress(filename,full=.TRUE.)
              WRITE (unit=desc_unit,&
                   fmt="(' <cp_fm_types>:{ id_nr=',i10,' written to files ',a,'*')",&
                   iostat=iostat) matrix%id_nr, TRIM(filename)
           END IF
           CALL cp_logger_generate_filename(logger=logger,res=filename,&
                root=TRIM(matrix%name)//"-"//&
                TRIM(ADJUSTL(cp_to_string(matrix%print_count))),&
                postfix=".desc",local=.TRUE.)
           INQUIRE (FILE=TRIM(filename),EXIST=exists)
           CALL cp_assert(.NOT.exists,cp_warning_level,cp_assertion_failed,&
                routineP,"file "//TRIM(filename)//" exists, overwriting",error)
           CALL open_file(TRIM(filename),file_status="unknown",&
                file_action="WRITE",&
                unit_number=desc_unit)

           CALL cp_logger_generate_filename(logger=logger,res=filename,&
             root=TRIM(matrix%name)//"-"//&
             TRIM(ADJUSTL(cp_to_string(matrix%print_count))),&
             postfix=".dat",local=.TRUE.)
           INQUIRE (FILE=TRIM(filename),EXIST=exists)
           CALL cp_assert(.NOT.exists,cp_warning_level,cp_assertion_failed,&
                routineP,"file "//TRIM(filename)//" exists, overwriting",error)
           CALL open_file(TRIM(filename),file_status="unknown",&
                file_action="WRITE",&
                unit_number=data_unit)
        END IF

        IF (loc .OR. para_env%mepos==para_env%source.or.long_description) THEN
           WRITE (unit=desc_unit,&
                fmt="(' <cp_fm_types>:{ id_nr=',i10,' ref_count=',i10,',')",&
                iostat=iostat) matrix%id_nr, matrix%ref_count
           WRITE (unit=desc_unit,&
                fmt="(' name=',a,',')",&
                iostat=iostat) matrix%name
           WRITE(unit=desc_unit,fmt="(a)",iostat=iostat) " matrix_structure="
           CALL cp_fm_struct_write(matrix%matrix_struct,unit_nr=desc_unit,&
                long_description=my_long_description,error=error)
           WRITE(unit=desc_unit,fmt="(a)",iostat=iostat) " ,"
           IF (my_long_description) THEN
              WRITE(unit=desc_unit,fmt="(a)",iostat=iostat) " local_data=("
!BEG:ORIG:LT:2015/02/03
              ! CALL cp_2d_r_write(matrix%local_data,unit_nr=data_unit,error=error)
!END:ORIG:LT:2015/02/03
!BEG:DEBUG:LT:2015/02/03
              CALL cp_2d_r_write(matrix%local_data, unit_nr=data_unit, &
                                 el_format="f9.4", &
                                 error=error)
!END:DEBUG:LT:2015/02/03
              IF (desc_unit/=data_unit) THEN
                 WRITE(unit=desc_unit,fmt="(a)",iostat=iostat) &
                      "  written to file "//TRIM(filename)
              END IF
           ELSE IF (loc) THEN
              IF (ASSOCIATED(matrix%local_data)) THEN
                 WRITE (unit=desc_unit,&
                      fmt="(' local_data=(REAL(wp, DIMENSION(',i8,', ',i8,'))')",&
                      iostat=iostat) SIZE(matrix%local_data,1), SIZE(matrix%local_data,2)
              END IF
           END IF
           WRITE (unit=desc_unit,fmt="(a)",iostat=iostat) " }"
        END IF

        IF (.not.loc.and.my_long_description) THEN
           CALL close_file(desc_unit)
           CALL close_file(data_unit)
        END IF
     ELSE
        IF (loc .OR. para_env%mepos==para_env%source) THEN
           WRITE (unit=unit_nr,fmt="(a)",iostat=iostat) " <cp_fm_types>:*null*"
        END IF
     END IF
   END SUBROUTINE cp_fm_write

! *****************************************************************************
!> \brief copy one identically sized matrix in the other
!> \param source ...
!> \param destination ...
!> \param error ...
!> \note
!>      see also cp_fm_to_fm_columns
! *****************************************************************************
  SUBROUTINE cp_fm_to_fm_matrix(source,destination,error)

    TYPE(cp_fm_type), POINTER                :: source, destination
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_fm_to_fm_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, npcol, nprow
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)

    failure=.TRUE.

    nprow = source%matrix_struct%context%num_pe(1)
    npcol = source%matrix_struct%context%num_pe(2)

    IF ((.NOT. cp2k_is_parallel).OR.&
         cp_fm_struct_equivalent(source%matrix_struct,&
         destination%matrix_struct,error=error)) THEN
       IF(source%use_sp.AND.destination%use_sp) THEN
          CALL cp_assert(SIZE(source%local_data_sp,1)==SIZE(destination%local_data_sp,1).AND.&
               SIZE(source%local_data_sp,2)==SIZE(destination%local_data_sp,2),cp_failure_level,&
               cp_assertion_failed, routineP,&
               "Cannot copy full matrix <"//TRIM(source%name)//&
               "> to full matrix <"//TRIM(destination%name)//&
               ">. The local_data blocks have different sizes."//&
CPSourceFileRef,&
               error=error)
          CALL scopy(SIZE(source%local_data_sp,1)*SIZE(source%local_data_sp,2),&
               source%local_data_sp(1,1),1,destination%local_data_sp(1,1),1)
       ELSEIF(source%use_sp.AND..NOT.destination%use_sp) THEN
          CALL cp_assert(SIZE(source%local_data_sp,1)==SIZE(destination%local_data,1).AND.&
               SIZE(source%local_data_sp,2)==SIZE(destination%local_data,2),cp_failure_level,&
               cp_assertion_failed, routineP,&
               "Cannot copy full matrix <"//TRIM(source%name)//&
               "> to full matrix <"//TRIM(destination%name)//&
               ">. The local_data blocks have different sizes."//&
CPSourceFileRef,&
               error=error)
          destination%local_data=REAL(source%local_data_sp,dp)
       ELSEIF(.NOT.source%use_sp.AND.destination%use_sp) THEN
          CALL cp_assert(SIZE(source%local_data,1)==SIZE(destination%local_data_sp,1).AND.&
               SIZE(source%local_data,2)==SIZE(destination%local_data_sp,2),cp_failure_level,&
               cp_assertion_failed, routineP,&
               "Cannot copy full matrix <"//TRIM(source%name)//&
               "> to full matrix <"//TRIM(destination%name)//&
               ">. The local_data blocks have different sizes."//&
CPSourceFileRef,&
               error=error)
          destination%local_data_sp=REAL(source%local_data,sp)
       ELSE
          CALL cp_assert(SIZE(source%local_data,1)==SIZE(destination%local_data,1).AND.&
               SIZE(source%local_data,2)==SIZE(destination%local_data,2),cp_failure_level,&
               cp_assertion_failed, routineP,&
               "Cannot copy full matrix <"//TRIM(source%name)//&
               "> to full matrix <"//TRIM(destination%name)//&
               ">. The local_data blocks have different sizes."//&
CPSourceFileRef,&
               error=error)
          CALL dcopy(SIZE(source%local_data,1)*SIZE(source%local_data,2),&
               source%local_data(1,1),1,destination%local_data(1,1),1)
       ENDIF
    ELSE
      CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
    END IF

    CALL timestop(handle)

  END SUBROUTINE cp_fm_to_fm_matrix

! *****************************************************************************
!> \brief copy just a subset of columns of a fm to a fm
!> \param msource ...
!> \param mtarget ...
!> \param ncol ...
!> \param source_start ...
!> \param target_start ...
! *****************************************************************************
  SUBROUTINE cp_fm_to_fm_columns(msource,mtarget,ncol,source_start,&
                                 target_start)

    TYPE(cp_fm_type), POINTER                :: msource, mtarget
    INTEGER, INTENT(IN)                      :: ncol
    INTEGER, INTENT(IN), OPTIONAL            :: source_start, target_start

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_fm_to_fm_columns', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, n, ss, ts
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a, b
#if defined(__SCALAPACK)
    INTEGER                                  :: i
    INTEGER, DIMENSION(9)                    :: desca, descb
#endif

    CALL timeset(routineN,handle)

    ss = 1
    ts = 1

    IF (PRESENT(source_start)) ss = source_start
    IF (PRESENT(target_start)) ts = target_start

    n = msource%matrix_struct%nrow_global

    a => msource%local_data
    b => mtarget%local_data

#if defined(__SCALAPACK)
    desca(:) = msource%matrix_struct%descriptor(:)
    descb(:) = mtarget%matrix_struct%descriptor(:)
    DO i=0,ncol-1
      CALL pdcopy(n,a(1,1),1,ss+i,desca,1,b(1,1),1,ts+i,descb,1)
    END DO
#else
    CALL dcopy(ncol*n,a(1,ss),1,b(1,ts),1)
#endif

    CALL timestop(handle)

  END SUBROUTINE cp_fm_to_fm_columns

! *****************************************************************************
!> \brief copy just a triangular matrix
!> \param msource ...
!> \param mtarget ...
!> \param uplo ...
!> \param error ...
! *****************************************************************************

  SUBROUTINE cp_fm_to_fm_triangular (msource, mtarget, uplo, error)

    TYPE(cp_fm_type), POINTER                :: msource, mtarget
    CHARACTER(LEN=*), INTENT(IN)             :: uplo
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_fm_to_fm_triangular', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ncol, nrow
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a, b
#if defined(__SCALAPACK)
    INTEGER, DIMENSION(9)                    :: desca, descb
#endif

    CALL timeset(routineN,handle)

    failure=.TRUE.

    nrow = msource%matrix_struct%nrow_global
    ncol = msource%matrix_struct%ncol_global

    a => msource%local_data
    b => mtarget%local_data

#if defined(__SCALAPACK)
    desca(:) = msource%matrix_struct%descriptor(:)
    descb(:) = mtarget%matrix_struct%descriptor(:)
    CALL pdlacpy(uplo,nrow, ncol,a(1,1),1,1,desca,b(1,1),1,1,descb)
#else
    CALL dlacpy(uplo,nrow, ncol,a(1,1),nrow,b(1,1),nrow)
#endif

    CALL timestop(handle)

  END SUBROUTINE cp_fm_to_fm_triangular

! *****************************************************************************
!> \brief copy just a part ot the matrix
!> \param msource ...
!> \param mtarget ...
!> \param nrow ...
!> \param ncol ...
!> \param s_firstrow ...
!> \param s_firstcol ...
!> \param t_firstrow ...
!> \param t_firstcol ...
!> \param error ...
! *****************************************************************************

  SUBROUTINE cp_fm_to_fm_submat(msource, mtarget, nrow, ncol, s_firstrow, s_firstcol, t_firstrow, t_firstcol, error)

    TYPE(cp_fm_type), POINTER                :: msource, mtarget
    INTEGER, INTENT(IN)                      :: nrow, ncol, s_firstrow, &
                                                s_firstcol, t_firstrow, &
                                                t_firstcol
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_fm_to_fm_submat', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, na, nb, ss, ts
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a, b
#if defined(__SCALAPACK)
    INTEGER, DIMENSION(9)                    :: desca, descb
#endif

    CALL timeset(routineN,handle)

    failure = .FALSE.

    a => msource%local_data
    b => mtarget%local_data


    na = msource%matrix_struct%nrow_global
    nb = mtarget%matrix_struct%nrow_global
!    nrow must be <= na and nb
    CALL cp_assert(nrow<=na, cp_failure_level,&
           cp_assertion_failed, routineP, &
           "cannot copy because nrow > number of rows of source matrix"//&
CPSourceFileRef,&
           error=error)
    CALL cp_assert(nrow<=nb, cp_failure_level,&
           cp_assertion_failed, routineP, &
           "cannot copy because nrow > number of rows of target matrix"//&
CPSourceFileRef,&
           error=error)
    na = msource%matrix_struct%ncol_global
    nb = mtarget%matrix_struct%ncol_global
!    ncol must be <= na_col and nb_col
    CALL cp_assert(ncol<=na, cp_failure_level,&
           cp_assertion_failed, routineP, &
           "cannot copy because nrow > number of rows of source matrix"//&
CPSourceFileRef,&
           error=error)
    CALL cp_assert(ncol<=nb, cp_failure_level,&
           cp_assertion_failed, routineP, &
           "cannot copy because nrow > number of rows of target matrix"//&
CPSourceFileRef,&
           error=error)

#if defined(__SCALAPACK)
    desca(:) = msource%matrix_struct%descriptor(:)
    descb(:) = mtarget%matrix_struct%descriptor(:)
    DO i=0,ncol-1
      ss = s_firstcol+i
      ts = t_firstcol+i
      CALL pdcopy(nrow,a(1,1),s_firstrow,ss,desca,1,b(1,1),t_firstrow,ts,descb,1)
    END DO
#else
    DO i = 0,ncol-1
      ss = s_firstcol+i
      ts = t_firstcol+i
      CALL dcopy(nrow,a(s_firstrow,ss),1,b(t_firstrow,ts),1)
    END DO
#endif


    CALL timestop(handle)
  END SUBROUTINE cp_fm_to_fm_submat

! *****************************************************************************
!> \brief General copy of a fm matrix to another fm matrix.
!>        The two matrices can have different contexts.
!>        Matrices have to have the same global shape.
!>
!>        Summary of distribution routines for dense matrices
!>        The following will copy A(iA:iA+M-1,jA:jA+N-1) to B(iB:iB+M-1,jB:jB+N-1):
!>
!>        call pdgemr2d(M,N,Aloc,iA,jA,descA,Bloc,iB,jB,descB,context)
!>
!>        Here A and B are global matrices, Aloc, and Bloc, are local matrices, 
!>        and context is a process grid that covers all parts of either A or B. 
!>        If a whole global matrix is to be distributed, use the call:
!>
!>        call pdgemr2d(M,N,A,1,1,descA,Bloc,1,1,descB,context)
!>
!>        A process that is not a part of the context of A should set descA(2) 
!>        to -1, and similarly for B.
!>
!> \param source          input fm matrix
!> \param destination     output fm matrix
!> \param global_context  process grid that covers all parts of either A or B.
!> \param error           CP2K error handling type
! *****************************************************************************
  SUBROUTINE cp_fm_copy_general(source,destination,global_context,error)

    TYPE(cp_fm_type), POINTER                :: source, destination
    INTEGER, INTENT(IN)                      :: global_context
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_fm_copy_general', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure, debug
    INTEGER                                  :: handle
#if defined(__SCALAPACK)
    INTEGER                                  :: sdat, ddat, m, n
    INTEGER, DIMENSION(9)                    :: desca, descb
    REAL(KIND=dp), DIMENSION(1,1), TARGET    :: dummy
    REAL(KIND=dp), DIMENSION(:,:), POINTER   :: smat, dmat
#endif

    CALL timeset(routineN,handle)

    failure=.FALSE.
    debug = debug_this_module

    IF (.NOT. cp2k_is_parallel) THEN
       CALL cp_fm_to_fm(source,destination,error)
    ELSE
#ifdef __SCALAPACK
       m = 1
       n = 1
       sdat = 0
       ddat = 0
       NULLIFY(smat,dmat)
       IF(ASSOCIATED(source)) THEN
          sdat = 1
          desca = source%matrix_struct%descriptor
          CALL cp_assert(.NOT.source%use_sp,&
               cp_failure_level, cp_assertion_failed, routineP, &
               "only DP kind implemnted"//&
CPSourceFileRef,&
               error=error)
          m = source%matrix_struct%nrow_global
          n = source%matrix_struct%ncol_global
          smat => source%local_data
       ELSE
          desca = -1
          smat => dummy
       END IF
       IF(ASSOCIATED(destination)) THEN
          ddat = 1
          descb = destination%matrix_struct%descriptor
          CALL cp_assert(.NOT.destination%use_sp,&
               cp_failure_level, cp_assertion_failed, routineP, &
               "only DP kind implemnted"//&
CPSourceFileRef,&
               error=error)
          m = destination%matrix_struct%nrow_global
          n = destination%matrix_struct%ncol_global
          dmat => destination%local_data
       ELSE
          descb = -1
          dmat => dummy
       END IF

       IF(sdat*ddat == 1) THEN
          CALL cp_assert(source%matrix_struct%nrow_global==&
               destination%matrix_struct%nrow_global,&
               cp_failure_level,&
               cp_assertion_failed, routineP, &
               "cannot copy between full matrixes of different sizes"//&
CPSourceFileRef,&
               error=error)
          CALL cp_assert(source%matrix_struct%ncol_global==&
               destination%matrix_struct%ncol_global,&
               cp_failure_level,&
               cp_assertion_failed, routineP, &
               "cannot copy between full matrixes of differen sizes"//&
CPSourceFileRef,&
               error=error)
       END IF

       CALL pdgemr2d(m,n,smat,1,1,desca,dmat,1,1,descb,global_context)

#else
       CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
#endif
    END IF

    CALL timestop(handle)

  END SUBROUTINE cp_fm_copy_general

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param irow_global ...
!> \param icol_global ...
!> \param alpha ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_fm_add_to_element(matrix,irow_global,icol_global,alpha,error)

    ! Add alpha to the matrix element specified by the global indices
    ! irow_global and icol_global

    ! - Creation (05.05.06,MK)

    TYPE(cp_fm_type), POINTER                :: matrix
    INTEGER, INTENT(IN)                      :: irow_global, icol_global
    REAL(KIND=dp), INTENT(IN)                :: alpha
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_fm_add_to_element', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: mypcol, myprow,npcol, nprow
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a
    TYPE(cp_blacs_env_type), POINTER         :: context
#if defined(__SCALAPACK)
    INTEGER                                  :: icol_local, ipcol, iprow, &
                                                irow_local
    INTEGER, DIMENSION(9)                    :: desca
#endif

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(matrix),cp_failure_level,routineP,error,failure)

    context => matrix%matrix_struct%context

    myprow = context%mepos(1)
    mypcol = context%mepos(2)

    nprow = context%num_pe(1)
    npcol = context%num_pe(2)

    a => matrix%local_data

#if defined(__SCALAPACK)

    desca(:) = matrix%matrix_struct%descriptor(:)

    CALL infog2l(irow_global,icol_global,desca,nprow,npcol,myprow,mypcol,&
                 irow_local,icol_local,iprow,ipcol)

    IF ((iprow == myprow).AND.(ipcol == mypcol)) THEN
      a(irow_local,icol_local) = a(irow_local,icol_local) + alpha
    END IF

#else

    a(irow_global,icol_global) = a(irow_global,icol_global) + alpha

#endif

  END SUBROUTINE cp_fm_add_to_element

! *****************************************************************************
!> \brief ...
!> \param fm ...
!> \param unit ...
!> \param error ...
! *****************************************************************************
 SUBROUTINE cp_fm_write_unformatted(fm,unit,error)
    TYPE(cp_fm_type), POINTER                :: fm
    INTEGER                                  :: unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_write_unformatted', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, j, max_block,&
                                                ncol_global, nrow_global
    LOGICAL                                  :: failure
    TYPE(cp_para_env_type), POINTER          :: para_env
#if defined(__SCALAPACK)
    INTEGER                                  :: i, i_block, icol_local,&
                                                ictxt_loc, in, info, ipcol, &
                                                iprow, irow_local, istat, &
                                                mepos, mypcol, myprow, &
                                                npcol, nprow, num_pe, rb, tag
    INTEGER, DIMENSION(9)                    :: desc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: vecbuf
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: newdat
    INTEGER, EXTERNAL                        :: numroc
#endif

   failure = .FALSE.
   CALL timeset(routineN,handle)
   CALL cp_fm_get_info(fm, nrow_global=nrow_global,ncol_global=ncol_global,ncol_block=max_block,&
        para_env=para_env,error=error)

#if defined(__SCALAPACK)
    num_pe=para_env%num_pe
    mepos =para_env%mepos
    rb=nrow_global
    tag =0
    ! get a new context
    ictxt_loc=para_env%group
    CALL cp_blacs_gridinit(ictxt_loc,'R',1,num_pe)
    CALL cp_blacs_gridinfo(ictxt_loc,nprow,npcol,myprow,mypcol)
    CALL descinit(desc,nrow_global,ncol_global,rb,max_block,0,0,ictxt_loc,nrow_global,info)
    CPPostcondition(info==0,cp_failure_level,routineP,error,failure)
    in=numroc(ncol_global,max_block,mypcol,0,npcol)

    ALLOCATE(newdat(nrow_global,MAX(1,in)))

            ! do the actual scalapack to cols reordering
    CALL pdgemr2d(nrow_global,ncol_global,fm%local_data(1,1),1,1,&
         fm%matrix_struct%descriptor, &
         newdat(1,1),1,1,desc,ictxt_loc)


   ALLOCATE(vecbuf(nrow_global*max_block),STAT=istat)
   CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
   DO i=1,ncol_global,MAX(max_block,1)
      i_block=MIN(max_block,ncol_global-i+1)
      CALL infog2l(1,i,desc,nprow,npcol,myprow,mypcol,&
               irow_local,icol_local,iprow,ipcol)
      IF (ipcol==mypcol) THEN
         DO j=1,i_block
            vecbuf((j-1)*nrow_global+1:nrow_global*j)=newdat(:,icol_local+j-1)
         END DO
      END IF

      IF (ipcol==0) THEN
         ! do nothing
      ELSE
         IF (ipcol==mypcol) THEN
            CALL mp_send(vecbuf(:),0,tag,para_env%group)
         ENDIF
         IF (mypcol==0) THEN
            CALL mp_recv(vecbuf(:),ipcol,tag,para_env%group)
         ENDIF
      ENDIF

      IF (unit>0)THEN
         DO j=1,i_block
            WRITE (unit) vecbuf((j-1)*nrow_global+1:nrow_global*j)
         END DO
      END IF

   END DO
   DEALLOCATE(vecbuf,STAT=istat)
   CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

   CALL cp_blacs_gridexit(ictxt_loc)

   DEALLOCATE(newdat,STAT=istat)
   CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

#else

   IF (unit>0)THEN
      DO j=1,ncol_global
         WRITE (unit) fm%local_data(:,j)
      END DO
   END IF

#endif
   CALL timestop(handle)

 END SUBROUTINE cp_fm_write_unformatted

! *****************************************************************************
!> \brief ...
!> \param fm ...
!> \param unit ...
!> \param error ...
! *****************************************************************************
 SUBROUTINE cp_fm_read_unformatted(fm,unit,error)
    TYPE(cp_fm_type), POINTER                :: fm
    INTEGER                                  :: unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_read_unformatted', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, j, max_block, &
                                                ncol_global, nrow_global
    LOGICAL                                  :: failure
    TYPE(cp_para_env_type), POINTER          :: para_env
#if defined(__SCALAPACK)
    INTEGER                                  :: k, n_cols, stat
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: vecbuf
#endif

   CALL timeset(routineN,handle)
   failure=.FALSE.

   CALL cp_fm_get_info(fm, nrow_global=nrow_global,ncol_global=ncol_global,ncol_block=max_block,&
        para_env=para_env,error=error)

#if defined(__SCALAPACK)

   ! the parallel case could be made more efficient (see cp_fm_write_unformatted)

   ALLOCATE(vecbuf(nrow_global,max_block),stat=stat)
   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

   DO j=1,ncol_global,max_block

      n_cols=MIN(max_block,ncol_global-j+1)
      IF (para_env%mepos==0) THEN
         DO k=1,n_cols
            READ(unit) vecbuf(:,k)
         ENDDO
      ENDIF
      CALL mp_bcast(vecbuf,0,para_env%group)
      CALL cp_fm_set_submatrix(fm,vecbuf,start_row=1,start_col=j,n_cols=n_cols,error=error)

   ENDDO

   DEALLOCATE(vecbuf,stat=stat)
   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

#else

   DO j=1,ncol_global
      READ (unit) fm%local_data(:,j)
   END DO

#endif

    CALL timestop(handle)

 END SUBROUTINE cp_fm_read_unformatted

! *****************************************************************************
!> \brief wrapper to scalapack function INDXG2P that computes the process 
!>         coordinate which posseses the entry of a distributed matrix specified 
!>         by a global index INDXGLOB.
!>
!>  Arguments
!>  =========
!>
!>  INDXGLOB  (global input) INTEGER
!>            The global index of the element.
!>
!>  NB        (global input) INTEGER
!>            Block size, size of the blocks the distributed matrix is
!>            split into.
!>
!>  IPROC     (local dummy) INTEGER
!>            Dummy argument in this case in order to unify the calling
!>            sequence of the tool-routines.
!>
!>  ISRCPROC  (global input) INTEGER
!>            The coordinate of the process that possesses the first
!>            row/column of the distributed matrix.
!>
!>  NPROCS    (global input) INTEGER
!>            The total number processes over which the matrix is
!>            distributed.         
!>
!> \param INDXGLOB ...
!> \param NB ...
!> \param IPROC ...
!> \param ISRCPROC ...
!> \param NPROCS ...
!> \retval G2P ...
!> \author Mauro Del Ben [MDB] - 12.2012
! *****************************************************************************
 FUNCTION  cp_fm_indxg2p( INDXGLOB, NB, IPROC, ISRCPROC, NPROCS ) RESULT(G2P)
   INTEGER, INTENT(IN)                      :: INDXGLOB, NB, IPROC, ISRCPROC, NPROCS
   INTEGER                                  :: G2P

   CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_fm_indxg2p', &
     routineP = moduleN//':'//routineN

#if defined(__SCALAPACK)
   INTEGER, EXTERNAL :: indxg2p
#endif

#if defined(__SCALAPACK)

   G2P = indxg2p( INDXGLOB, NB, IPROC, ISRCPROC, NPROCS )

#else

   G2P = 0

#endif

 END FUNCTION  cp_fm_indxg2p

! *****************************************************************************
!> \brief wrapper to scalapack function INDXG2L that computes the local index 
!>         of a distributed matrix entry pointed to by the global index INDXGLOB.
!>
!>  Arguments
!>  =========
!>
!>  INDXGLOB  (global input) INTEGER
!>            The global index of the distributed matrix entry.
!>
!>  NB        (global input) INTEGER
!>            Block size, size of the blocks the distributed matrix is
!>            split into.
!>
!>  IPROC     (local dummy) INTEGER
!>            Dummy argument in this case in order to unify the calling
!>            sequence of the tool-routines.
!>
!>  ISRCPROC  (local dummy) INTEGER
!>            Dummy argument in this case in order to unify the calling
!>            sequence of the tool-routines.
!>
!>  NPROCS    (global input) INTEGER
!>            The total number processes over which the distributed
!>            matrix is distributed.
!>
!> \param INDXGLOB ...
!> \param NB ...
!> \param IPROC ...
!> \param ISRCPROC ...
!> \param NPROCS ...
!> \retval G2L ...
!> \author Mauro Del Ben [MDB] - 12.2012
! *****************************************************************************
 FUNCTION  cp_fm_indxg2l( INDXGLOB, NB, IPROC, ISRCPROC, NPROCS ) RESULT(G2L)
   INTEGER, INTENT(IN)                      :: INDXGLOB, NB, IPROC, ISRCPROC, NPROCS
   INTEGER                                  :: G2L

   CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_fm_indxg2l', &
     routineP = moduleN//':'//routineN

#if defined(__SCALAPACK)
   INTEGER, EXTERNAL :: indxg2l
#endif

#if defined(__SCALAPACK)

   G2L = indxg2l( INDXGLOB, NB, IPROC, ISRCPROC, NPROCS )

#else

   G2L = INDXGLOB

#endif

 END FUNCTION  cp_fm_indxg2l

! *****************************************************************************
!> \brief wrapper to scalapack function INDXL2G that computes the global index
!>         of a distributed matrix entry pointed to by the local index INDXLOC 
!>         of the process indicated by IPROC.
!>
!>  Arguments
!>  =========
!>
!>  INDXLOC   (global input) INTEGER
!>            The local index of the distributed matrix entry.
!>
!>  NB        (global input) INTEGER
!>            Block size, size of the blocks the distributed matrix is
!>            split into.
!>
!>  IPROC     (local input) INTEGER
!>            The coordinate of the process whose local array row or
!>            column is to be determined.
!>
!>  ISRCPROC  (global input) INTEGER
!>            The coordinate of the process that possesses the first
!>            row/column of the distributed matrix.
!>
!>  NPROCS    (global input) INTEGER
!>            The total number processes over which the distributed
!>            matrix is distributed.
!>
!> \param INDXLOC ...
!> \param NB ...
!> \param IPROC ...
!> \param ISRCPROC ...
!> \param NPROCS ...
!> \retval L2G ...
!> \author Mauro Del Ben [MDB] - 12.2012
! *****************************************************************************
 FUNCTION  cp_fm_indxl2g( INDXLOC, NB, IPROC, ISRCPROC, NPROCS ) RESULT(L2G)
   INTEGER, INTENT(IN)                      :: INDXLOC, NB, IPROC, ISRCPROC, NPROCS 
   INTEGER                                  :: L2G

   CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_fm_indxl2g', &
     routineP = moduleN//':'//routineN

#if defined(__SCALAPACK)
   INTEGER, EXTERNAL :: indxl2g
#endif

#if defined(__SCALAPACK)

   L2G = indxl2g( INDXLOC, NB, IPROC, ISRCPROC, NPROCS )

#else

   L2G = INDXLOC

#endif

 END FUNCTION  cp_fm_indxl2g

! *****************************************************************************
!> \brief ...
!> \param mult_type ...
! *****************************************************************************
 SUBROUTINE cp_fm_setup(mult_type)
    INTEGER                                  :: mult_type

   mm_type=mult_type
 END SUBROUTINE cp_fm_setup

! *****************************************************************************
!> \brief ...
!> \retval res ...
! *****************************************************************************
 FUNCTION cp_fm_get_mm_type()  RESULT(res)
    INTEGER                                  :: res

   res=mm_type
 END FUNCTION

END MODULE cp_fm_types
