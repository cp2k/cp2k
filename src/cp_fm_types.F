!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2002  CP2K developers group
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_fm_types [1.0] *
!!
!!   NAME
!!     cp_fm_types
!!
!!   FUNCTION
!!     represent a full matrix distribued on many processors
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Matthias Krack (22.05.2001)
!!
!!   MODIFICATION HISTORY
!!     1) Changed the size work array for syevx (Joost VandeVondele 28.02.02)
!!     1) Introduced some cholesky guys         (Joost VandeVondele    05.02)
!!     2) Introduced a symmetric matrix x distributed vector multiply ( JVdV)
!!     3) separated structure object, removed globenv, renamed to full matrix
!!        many changes (fawzi 08.2002)
!!
!!   SOURCE
!****************************************************************************
MODULE cp_fm_types
  USE cp_array_utils,                  ONLY: cp_2d_r_output,&
                                             cp_2d_r_write
  USE cp_blacs_env,                    ONLY: cp_blacs_env_from_globenv,&
                                             cp_blacs_env_release,&
                                             cp_blacs_env_retain
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error,&
                                             cp_unimplemented_error
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_equivalent,&
                                             cp_fm_struct_get,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_retain,&
                                             cp_fm_struct_type,&
                                             cp_fm_struct_write
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_logger_generate_filename,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_env,                     ONLY: cp_para_env_from_globenv,&
                                             cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: wp=>dp
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_max,&
                                             mp_sum,&
                                             mp_sync
  USE qs_parser,                       ONLY: close_file,&
                                             open_file
  USE string_utilities,                ONLY: compress_all
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_fm_types'

  INTEGER,SAVE,PRIVATE :: last_fm_id_nr=0


  PUBLIC :: cp_full_matrix_type, cp_full_matrix_p_type
  PUBLIC :: cp_fm_create, cp_fm_create2, cp_fm_retain, cp_fm_release,&
       cp_fm_init_random, cp_fm_set_element, cp_fm_set_all,&
       cp_fm_set_submatrix, cp_fm_get_submatrix, cp_fm_get_info,&
       cp_fm_maxval, cp_fm_write, cp_fm_get_element
!***
!****************************************************************************

!!****s* cp_fm_types/cp_full_matrix_type [1.0] *
!!
!!   NAME
!!     cp_full_matrix_type
!!
!!   FUNCTION
!!     represent a full matrix
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - name: the name of the matrix, used for printing
!!     - id_nr: identification number that makes this matrix unique
!!     - ref_cont: reference count (see doc/ReferenceCounting.html)
!!     - print_count: sequential number to make the printed files unique
!!     - matrix_struct: structure of this matrix
!!     - local_data: array with the data of the matrix (its contents
!!       depend on the matrix type used: in parallel runs it will be
!!       in scalapack format, in sequential, it will simply contain
!!       the matrix)
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!   SOURCE
  !***************************************************************************
  TYPE cp_full_matrix_type
!    PRIVATE
     CHARACTER(LEN=60) :: name
     INTEGER :: id_nr, ref_count, print_count
     TYPE(cp_fm_struct_type), POINTER :: matrix_struct
     REAL(kind=wp), DIMENSION(:,:), POINTER :: local_data
  END TYPE cp_full_matrix_type
!!***
  !****************************************************************************

!!****s* cp_fm_types/cp_full_matrix_p_type [1.0] *
!!
!!   NAME
!!     cp_full_matrix_p_type
!!
!!   FUNCTION
!!     just to build arrays of pointers to matrixes
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - matrix: the pointer to the matrix
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!   SOURCE
  !***************************************************************************
  TYPE cp_full_matrix_p_type
     TYPE(cp_full_matrix_type), POINTER :: matrix
  END TYPE cp_full_matrix_p_type
!!***
  !****************************************************************************

CONTAINS

  SUBROUTINE cp_fm_create2(new_matrix,nrow_global,ncol_global,&
       nrow_block,ncol_block,name,globenv,&
       context,para_env,error)

!   Purpose: Allocate a new distributed BLACS matrix.

!   History: - Creation (23.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(cp_full_matrix_type), POINTER        :: new_matrix
    CHARACTER(LEN=*), INTENT(IN)              :: name
    INTEGER, INTENT(IN)                       :: ncol_block,&
         ncol_global,nrow_block,&
         nrow_global
    TYPE(global_environment_type), INTENT(IN), OPTIONAL :: globenv
    TYPE(cp_blacs_env_type), POINTER, OPTIONAL          :: context
    TYPE(cp_para_env_type), POINTER, OPTIONAL           :: para_env
    TYPE(cp_error_type), INTENT(inout),OPTIONAL         :: error

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routineN ="cp_fm_create2",&
         routineP=moduleN//':'//routineN

!   *** Local variables ***

    CHARACTER(LEN=40) :: message
    TYPE(cp_fm_struct_type), POINTER :: matrix_struct
    TYPE(cp_blacs_env_type), POINTER       :: ctxt
    TYPE(cp_para_env_type), POINTER        :: para_e
    LOGICAL :: failure
    failure=.FALSE.

    NULLIFY(ctxt,para_e)
    IF (PRESENT(globenv)) THEN
       ctxt => cp_blacs_env_from_globenv(globenv)
       para_e => cp_para_env_from_globenv(globenv)
    END IF
    IF (PRESENT(context)) ctxt => context
    IF (PRESENT(para_env)) para_e => para_env
    CPPrecondition(ASSOCIATED(ctxt),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(para_e),cp_failure_level,routineP,error,failure)
    CALL cp_fm_struct_create(matrix_struct,nrow_global=nrow_global,&
         ncol_global=ncol_global,nrow_block=nrow_block,&
         ncol_block=ncol_block, context=ctxt,&
         para_env=para_e)
    CALL cp_fm_create(new_matrix,matrix_struct=matrix_struct)
    CALL cp_fm_struct_release(matrix_struct)
  END SUBROUTINE cp_fm_create2

!!****f* cp_fm_types/cp_fm_create [1.0] *
!!
!!   NAME
!!     cp_fm_create
!!
!!   SYNOPSIS
!!     Subroutine cp_fm_create(matrix, matrix_struct, name, error)
!!       Type(cp_full_matrix_type), Pointer:: matrix
!!       Type(cp_fm_struct_type), Pointer:: matrix_struct
!!       Character(Len=*), Intent (IN), Optional:: name
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_fm_create
!!
!!   FUNCTION
!!     creates a new full matrix with the given structure
!!
!!   NOTES
!!     preferred allocation routine
!!
!!   INPUTS
!!     - matrix: the matrix to be created
!!     - matrix_struct: the structure of matrix
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE cp_fm_create(matrix,matrix_struct,name,error)
    TYPE(cp_full_matrix_type), POINTER :: matrix
    TYPE(cp_fm_struct_type), POINTER :: matrix_struct
    CHARACTER(len=*), INTENT(in), OPTIONAL :: name
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='cp_fm_create',&
         routineP=moduleN//':'//routineN
    INTEGER :: stat, nprow, npcol, iprow, ipcol, nrow_local, ncol_local
    TYPE(cp_blacs_env_type), POINTER :: context
    failure=.FALSE.

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

       nprow=context%num_pe(1)
       npcol=context%num_pe(2)
       NULLIFY(matrix%local_data)

       ! OK, we allocate here at least a 1 x 1 matrix 
       ! this must (and is) compatible with the descinit call
       ! in cp_fm_struct
       nrow_local=matrix_struct%local_leading_dimension
       ncol_local=MAX(1,matrix_struct%ncol_locals(context%mepos(2)))
       ALLOCATE(matrix%local_data(nrow_local,ncol_local),stat=stat)

       CALL dcopy(nrow_local*ncol_local,0.0_wp,0,matrix%local_data,1)
! matrix%local_data(:,:) = 0.0_wp

       IF (PRESENT(name)) THEN
          matrix%name=name
       ELSE
          matrix%name='full matrix'//cp_to_string(matrix%id_nr)
       END IF
    END IF
  END SUBROUTINE cp_fm_create
!***************************************************************************

!!****f* cp_fm_types/cp_fm_retain [1.0] *
!!
!!   NAME
!!     cp_fm_retain
!!
!!   SYNOPSIS
!!     Subroutine cp_fm_retain(matrix, error)
!!       Type(cp_full_matrix_type), Pointer:: matrix
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_fm_retain
!!
!!   FUNCTION
!!     retains a full matrix
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - matrix: the matrix to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE cp_fm_retain(matrix,error)
    TYPE(cp_full_matrix_type), POINTER :: matrix
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='cp_fm_retain',&
         routineP=moduleN//':'//routineN

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(matrix),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPrecondition(matrix%ref_count>0,cp_failure_level,routineP,error,failure)
       matrix%ref_count=matrix%ref_count+1
    END IF

  END SUBROUTINE cp_fm_retain
!***************************************************************************

!!****f* cp_fm_types/cp_fm_release [1.0] *
!!
!!   NAME
!!     cp_fm_release
!!
!!   SYNOPSIS
!!     Subroutine cp_fm_release(matrix, error)
!!       Type(cp_full_matrix_type), Pointer:: matrix
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_fm_release
!!
!!   FUNCTION
!!     releases a full matrix
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - matrix: the matrix to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE cp_fm_release(matrix,error)
    TYPE(cp_full_matrix_type), POINTER :: matrix
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='cp_fm_release',&
         routineP=moduleN//':'//routineN
    INTEGER :: iprow,ipcol, stat
    TYPE(cp_blacs_env_type), POINTER :: context

    failure=.FALSE.

    IF (ASSOCIATED(matrix)) THEN
       CPPrecondition(matrix%ref_count>0,cp_failure_level,routineP,error,failure)
       matrix%ref_count=matrix%ref_count-1
       context=> matrix%matrix_struct%context
       IF (matrix%ref_count<1) THEN
          IF (ASSOCIATED(matrix%local_data)) THEN
             DEALLOCATE(matrix%local_data,stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
          matrix%name=""
          CALL cp_fm_struct_release(matrix%matrix_struct,error=error)
          DEALLOCATE(matrix,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
    END IF
    NULLIFY(matrix)

  END SUBROUTINE cp_fm_release
!***************************************************************************

! *****************************************************************************
! init ncol vectors of a full matrix with random numbers
! *****************************************************************************
  SUBROUTINE cp_fm_init_random(matrix,ncol,start_col,error)

    TYPE(cp_full_matrix_type), POINTER           :: matrix
    INTEGER, INTENT(IN)                          :: ncol
    INTEGER, INTENT(IN), OPTIONAL                :: start_col
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

    INTEGER :: handle,n,i,my_start_col
    REAL(wp), DIMENSION(:,:), POINTER :: a,buff
    INTEGER :: nprow,npcol,myprow,mypcol,nrow_block,ncol_block
    INTEGER, DIMENSION(9) :: desca
    INTEGER, DIMENSION(4),SAVE :: ISEED
    LOGICAL, SAVE :: FIRST = .TRUE.
    LOGICAL :: failure
    INTEGER :: icol
    TYPE(cp_blacs_env_type), POINTER :: context
    CHARACTER(len=*), PARAMETER :: routineN='cp_fm_init_random',&
                                   routineP=moduleN//':'//routineN

    CALL timeset("cp_fm_init_random","I","",handle)
    IF (FIRST) THEN
       ISEED(1)=4
       ISEED(2)=3
       ISEED(3)=2
       ISEED(4)=1 ! last one has to be odd
       FIRST=.FALSE.
    ENDIF

    my_start_col = 1
    IF (PRESENT(start_col)) my_start_col=start_col

    n = matrix%matrix_struct%nrow_global
    CPPostcondition(matrix%matrix_struct%ncol_global>=(my_start_col+ncol-1),cp_failure_level,routineP,error,failure)
    ALLOCATE(buff(n,1))
    DO icol=my_start_col, my_start_col+ncol-1
       IF (matrix%matrix_struct%para_env%mepos.eq.0) THEN
          CALL dlarnv(1,iseed,n,buff(1,1))
       ENDIF
       CALL mp_bcast(buff,0,matrix%matrix_struct%para_env%group) 
       ! write(6,*) buff
       CALL cp_fm_set_submatrix(matrix,buff,start_col=icol)

    ENDDO
    ! write(6,*) matrix%local_data
    DEALLOCATE(buff)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE cp_fm_init_random
! *****************************************************************************

  SUBROUTINE cp_fm_set_all(matrix,alpha,beta,error)

!   Purpose: Set the FULL matrix elements to alpha.
!   if beta is present the diagonal gets the value beta

!   ***************************************************************************

    TYPE(cp_full_matrix_type), POINTER           :: matrix
    REAL(wp), INTENT(IN)                         :: alpha
    REAL(wp), INTENT(IN), optional               :: beta
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   *** Local variables ***

    INTEGER :: matrix_size,i

!   ---------------------------------------------------------------------------

    matrix_size = size(matrix%local_data,1)*size(matrix%local_data,2)
    CALL dcopy (matrix_size, alpha, 0, matrix%local_data, 1)
    IF (PRESENT(beta)) THEN
      DO i=1,matrix%matrix_struct%nrow_global
         CALL cp_fm_set_element(matrix,i,i,beta)
      ENDDO
    ENDIF

  END SUBROUTINE cp_fm_set_all

! *****************************************************************************

!!****f* cp2k/cp_fm_types/cp_fm_get_element *
!!
!!   NAME
!!     cp_fm_get_element
!!
!!   FUNCTION
!!     return an element of a full_matrix
!!
!!   ARGUMENTS
!!     - matrix: the matrix to read
!!     - irow_global: the row
!!     - icol_global: the col
!!     - alpha: the value of matrix(irow_global, icol_global)
!!     - local: true if the element is on this cpu, false otherwise
!!
!!   NOTES
!!     - modified semantics. now this function always returns the value 
!!       previously the value was zero on cpus that didn't own the relevant 
!!       part of the matrix (Joost VandeVondele, May 2003)
!!     - usage of the function should be avoided, as it is likely to rather slow
!!       using row_indices/col_indices/local_data + some smart scheme normally 
!!       yields a real parallel code
!!
!!****
  SUBROUTINE cp_fm_get_element(matrix,irow_global,icol_global,alpha,local)
    
    IMPLICIT NONE

    ! arguments
    TYPE(cp_full_matrix_type), POINTER          :: matrix
    REAL(wp), INTENT(OUT)                     :: alpha
    INTEGER, INTENT(IN)                       :: icol_global,&
                                                 irow_global
    LOGICAL, INTENT(OUT), OPTIONAL            :: local

    ! locals
    INTEGER :: icol_local,ipcol,iprow,irow_local,mypcol,myprow,npcol,&
               nprow
    INTEGER, DIMENSION(9) :: desca
    REAL(wp), DIMENSION(:,:), POINTER :: a
    TYPE(cp_blacs_env_type), POINTER :: context


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
      CALL DGEBS2D(context%group, 'All', ' ', 1, 1, alpha , 1 )
      IF (PRESENT(local)) local=.TRUE.
    ELSE
      CALL DGEBR2D(context%group, 'All', ' ', 1, 1, alpha , 1 , iprow, ipcol )
      IF (PRESENT(local)) local=.FALSE.
    END IF

#else 
    IF (PRESENT(local)) local=.TRUE.
    alpha = matrix%local_data(irow_global,icol_global)
#endif

  END SUBROUTINE cp_fm_get_element

! *****************************************************************************

  SUBROUTINE cp_fm_set_element(matrix,irow_global,icol_global,alpha,error)

!   Purpose: Set the FULL matrix element (irow_global,icol_global) to alpha.

!   History: - Creation (08.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(cp_full_matrix_type), POINTER          :: matrix
    REAL(wp), INTENT(IN)                      :: alpha
    INTEGER, INTENT(IN)                       :: icol_global,&
         irow_global
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   *** Local variables ***

    INTEGER :: icol_local,ipcol,iprow,irow_local,mypcol,myprow,npcol,&
         nprow

    INTEGER, DIMENSION(9) :: desca

    REAL(wp), DIMENSION(:,:), POINTER :: a
    TYPE(cp_blacs_env_type), POINTER :: context

!   ---------------------------------------------------------------------------

    context => matrix%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)

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

!!****f* cp_fm_types/cp_fm_set_submatrix [1.0] *
!!
!!   NAME
!!     cp_fm_set_submatrix
!!
!!   SYNOPSIS
!!     Subroutine cp_fm_set_submatrix(full_matrix, new_values, start_row,&
!!         start_col, n_rows, n_cols, alpha, beta, transpose_new_val,&
!!         error)
!!       Type(cp_full_matrix_type), Pointer:: full_matrix
!!       Real(Kind=wp), Dimension(:,:), Intent (IN):: new_values
!!       Integer, Intent (IN), Optional:: start_row, start_col, n_rows,&
!!         n_cols
!!       Real(Kind=wp), Intent (IN), Optional:: alpha, beta
!!       Logical, Intent (IN), Optional:: transpose_new_val
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine cp_fm_set_submatrix
!!
!!   FUNCTION
!!     sets a submatrix of a full matrix
!!      full_matrix(start_row:start_row+n_rows,start_col:start_col+n_cols)
!!      = alpha*op(new_values)(1:n_rows,1:n_cols)+ beta 
!!      * full_matrix(start_row:start_row+n_rows,start_col:start_col+n_cols)
!!
!!   NOTES
!!     optimized for full column updates and alpha=1.0, beta=0.0
!!
!!   INPUTS
!!     - new_values: a replicated full matrix with the new values
!!     - full_matrix: the full to change
!!     - start_row: the starting row of b_matrix (defaults to 1)
!!     - start_col: the starting col of b_matrix (defaults to 1)
!!     - n_rows: the number of row to change in b (defaults to
!!       size(op(new_values),1))
!!     - n_cols: the number of columns to change in b (defaults to
!!       size(op(new_values),2))
!!     - alpha: rescaling factor for the new values (defaults to 1.0)
!!     - beta: rescaling factor for the old values (defaults to 0.0)
!!     - transpose_a: if new_values should be transposed: if true
!!       op(new_values)=new_values^T, else op(new_values)=new_values
!!       (defaults to false)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2002 created borrowing from Joost's blacs_replicated_copy [fawzi]
!!
!!*** *********************************************************************
  SUBROUTINE cp_fm_set_submatrix(full_matrix,new_values,start_row,&
       start_col, n_rows, n_cols, alpha, beta, transpose_new_val, error)
    TYPE(cp_full_matrix_type) , POINTER :: full_matrix
    REAL(wp), DIMENSION(:,:), INTENT(in) :: new_values
    INTEGER, INTENT(in), OPTIONAL :: start_row, start_col,&
         n_rows, n_cols
    REAL(kind=wp), INTENT(in), OPTIONAL :: alpha, beta
    LOGICAL, INTENT(in), OPTIONAL :: transpose_new_val
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

    REAL(kind=wp), DIMENSION(:,:), POINTER :: full_block
    REAL(kind=wp) :: al,be
    INTEGER :: i0,j0,ncol,nrow,i,j,nrow_global,ncol_global,ncol_block,&
         this_col,this_row,nrow_local,ncol_local,nrow_block,handle
    INTEGER, DIMENSION(:), POINTER :: row_indices,col_indices
    CHARACTER(len=*), PARAMETER :: routineN='cp_fm_set_submatrix',&
         routineP=moduleN//':'//routineN
    LOGICAL :: tr_a, failure
    TYPE(cp_para_env_type), POINTER :: para_env
    TYPE(cp_blacs_env_type), POINTER :: context
    al=1.0_wp; be=0.0_wp; i0=1; j0=1; tr_a=.FALSE.
    ! can be called too many times, making it a bit useless
    ! CALL timeset(routineN//','//moduleN,"I","",handle)

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(full_matrix),cp_failure_level,routineP,error,failure)
    IF (.not.failure) THEN
       CPPrecondition(full_matrix%ref_count>0,cp_failure_level,routineP,error,failure)
    END IF
    IF (.not.failure) THEN
       IF (PRESENT(alpha)) al=alpha
       IF (PRESENT(beta)) be=beta
       IF (PRESENT(start_row)) i0=start_row
       IF (PRESENT(start_col)) j0=start_col
       IF (PRESENT(transpose_new_val)) tr_a=transpose_new_val
       IF (tr_a) THEN
          nrow=SIZE(new_values,2)
          ncol=SIZE(new_values,1)
       ELSE
          nrow=SIZE(new_values,1)
          ncol=SIZE(new_values,2)
       END IF
       IF (PRESENT(n_rows)) nrow=n_rows
       IF (PRESENT(n_cols)) ncol=n_cols

       para_env => full_matrix%matrix_struct%para_env
       context => full_matrix%matrix_struct%context

       full_block => full_matrix%local_data

       CALL cp_fm_get_info(matrix=full_matrix,&
            nrow_global=nrow_global,ncol_global=ncol_global,&
            nrow_block =nrow_block ,ncol_block =ncol_block ,&
            nrow_local =nrow_local ,ncol_local =ncol_local ,&
            row_indices=row_indices,col_indices=col_indices)

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

    ! CALL timestop(0.0_wp,handle)

  END SUBROUTINE cp_fm_set_submatrix

!!****f* cp_fm_types/cp_fm_get_submatrix [1.0] *
!!
!!   NAME
!!     cp_fm_get_submatrix
!!
!!   SYNOPSIS
!!     Subroutine cp_fm_get_submatrix(full_matrix, target_m, start_row,&
!!         start_col, n_rows, n_cols, transpose_target_m, error)
!!       Type(cp_full_matrix_type), Pointer:: full_matrix
!!       Real(Kind=wp), Dimension(:,:), Intent (OUT):: target_m
!!       Integer, Intent (IN), Optional:: start_row, start_col, n_rows,&
!!         n_cols
!!       Logical, Intent (IN), Optional:: transpose_target_m
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine cp_fm_get_submatrix
!!
!!   FUNCTION
!!     gets a submatrix of a full matrix
!!      op(target_m)(1:n_rows,1:n_cols)
!!      =full_matrix(start_row:start_row+n_rows,start_col:start_col+n_cols)
!!
!!   NOTES
!!     optimized for full column updates. Zeros out a little too much
!!     of target_m
!!
!!   INPUTS
!!     - target_m: a replicated full matrix that will contain the result
!!     - full_matrix: the full you want to get the info from
!!     - start_row: the starting row of b_matrix (defaults to 1)
!!     - start_col: the starting col of b_matrix (defaults to 1)
!!     - n_rows: the number of row to change in b (defaults to
!!       size(op(new_values),1))
!!     - n_cols: the number of columns to change in b (defaults to
!!       size(op(new_values),2))
!!     - transpose_target_m: if target_m should be transposed: if true
!!       op(target_m)=target_m^T, else op(target_m)=target_m
!!       (defaults to false)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2002 created borrowing from Joost's blacs_replicated_copy [fawzi]
!!
!!*** *********************************************************************
  SUBROUTINE cp_fm_get_submatrix(full_matrix,target_m, start_row,&
       start_col, n_rows, n_cols, transpose_target_m, error)
    TYPE(cp_full_matrix_type) , POINTER :: full_matrix
    REAL(wp), DIMENSION(:,:), INTENT(out) :: target_m
    INTEGER, INTENT(in), OPTIONAL :: start_row, start_col,&
         n_rows, n_cols
    LOGICAL, INTENT(in), OPTIONAL :: transpose_target_m
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

    REAL(kind=wp), DIMENSION(:,:), POINTER :: full_block
    INTEGER :: i0,j0,ncol,nrow,i,j,nrow_global,ncol_global,ncol_block,&
         this_col,this_row,nrow_local,ncol_local,nrow_block,handle
    INTEGER, DIMENSION(:), POINTER :: row_indices,col_indices
    CHARACTER(len=*), PARAMETER :: routineN='cp_fm_get_submatrix',&
         routineP=moduleN//':'//routineN
    LOGICAL :: tr_a, failure
    TYPE(cp_para_env_type), POINTER :: para_env
    TYPE(cp_blacs_env_type), POINTER :: context

    CALL timeset(routineN//','//moduleN,"I","",handle)

    i0=1; j0=1; tr_a=.FALSE.
    failure=.FALSE.

    CPPrecondition(ASSOCIATED(full_matrix),cp_failure_level,routineP,error,failure)
    IF (.not.failure) THEN
       CPPrecondition(full_matrix%ref_count>0,cp_failure_level,routineP,error,failure)
    END IF
    IF (.not.failure) THEN
       IF (PRESENT(start_row)) i0=start_row
       IF (PRESENT(start_col)) j0=start_col
       IF (PRESENT(transpose_target_m)) tr_a=transpose_target_m
       IF (tr_a) THEN
          nrow=SIZE(target_m,2)
          ncol=SIZE(target_m,1)
       ELSE
          nrow=SIZE(target_m,1)
          ncol=SIZE(target_m,2)
       END IF
       IF (PRESENT(n_rows)) nrow=n_rows
       IF (PRESENT(n_cols)) ncol=n_cols

       para_env => full_matrix%matrix_struct%para_env
       context => full_matrix%matrix_struct%context

       full_block => full_matrix%local_data
#if defined(__SCALAPACK)
       ! zero-out whole target_m 
       CALL dcopy(SIZE(target_m,1)*SIZE(target_m,2),0.0_wp,0,target_m(1,1),1)
#endif

       CALL cp_fm_get_info(matrix=full_matrix, &
            nrow_global=nrow_global,ncol_global=ncol_global, &
            nrow_block =nrow_block ,ncol_block =ncol_block , &
            nrow_local =nrow_local ,ncol_local =ncol_local , &
            row_indices=row_indices,col_indices=col_indices)


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

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE cp_fm_get_submatrix

  SUBROUTINE cp_fm_get_info(matrix,name,nrow_global,ncol_global,&
       nrow_block,ncol_block,nrow_local,ncol_local,&
       row_indices,col_indices,local_data,context,&
       nrow_locals, ncol_locals, matrix_struct,para_env,error)

!   Purpose: Return informations about the specified FULL matrix.
!   the indices vectors are cached.
!   indices is an array that maps the local indices to the global ones
!   local data is a little hole to access the blacs block, needed in some rare cases
!   for efficient use

!   History: - Creation (08.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(cp_full_matrix_type), POINTER        :: matrix
    CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: name
    INTEGER, OPTIONAL, INTENT(OUT)          :: ncol_block,ncol_global,&
         nrow_block,nrow_global,&
         nrow_local,ncol_local
    INTEGER, OPTIONAL, DIMENSION(:), POINTER   :: row_indices,col_indices,&
         nrow_locals, ncol_locals
    TYPE(cp_para_env_type), POINTER, OPTIONAL :: para_env
    TYPE(cp_blacs_env_type), POINTER, OPTIONAL :: context
    TYPE(cp_fm_struct_type),POINTER,OPTIONAL :: matrix_struct
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error
    REAL(wp), DIMENSION(:,:), POINTER, OPTIONAL :: local_data

    CHARACTER(len=*), PARAMETER :: routineN='cp_fm_get_info',&
         routineP=moduleN//':'//routineN
    INTEGER i,nprow,npcol,myprow,mypcol, stat
    TYPE(cp_blacs_env_type), POINTER :: ctxt
    LOGICAL :: failure
#if defined(__SCALAPACK)
    INTEGER , EXTERNAL :: indxl2g
#endif

!   ---------------------------------------------------------------------------
    failure=.FALSE.
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
       nrow_locals => matrix%matrix_struct%ncol_locals
    END IF

    CALL cp_fm_struct_get(matrix%matrix_struct, nrow_local=nrow_local,&
          ncol_local=ncol_local, nrow_global=nrow_global,&
          ncol_global=ncol_global, nrow_block=nrow_block, &
          ncol_block=ncol_block, error=error)

    IF (PRESENT(para_env)) para_env => matrix%matrix_struct%para_env

  END SUBROUTINE cp_fm_get_info

! *****************************************************************************

  SUBROUTINE cp_fm_maxval(matrix,a_max,error)

!   Purpose: Get the maximum absolute element of a BLACS matrix.

!   History: - Creation (11.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(cp_full_matrix_type), POINTER        :: matrix
    REAL(wp), INTENT(OUT)                     :: a_max
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

!   *** Local variables ***

    INTEGER :: handle

    REAL(wp), DIMENSION(:,:), POINTER :: my_block

!   ---------------------------------------------------------------------------

    CALL timeset("cp_fm_maxval","I","",handle)

    my_block => matrix%local_data

    a_max = MAXVAL(ABS(my_block))

    CALL mp_max(a_max,matrix%matrix_struct%para_env%group)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE cp_fm_maxval

!!****f* cp_fm_types/cp_fm_write [1.0] *
!!
!!   NAME
!!     cp_fm_write
!!
!!   SYNOPSIS
!!     Subroutine cp_fm_write(matrix, unit_nr, long_description, local,&
!!         error)
!!       Type(cp_full_matrix_type), Pointer:: matrix
!!       Integer, Intent (IN):: unit_nr
!!       Logical, Optional, Intent (IN):: long_description, local
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_fm_write
!!
!!   FUNCTION
!!     writes out the matrix
!!
!!   NOTES
!!     uniquing_nr should be smartened: don't add if no file exist, else
!!     increase nr until no file exists.
!!
!!   INPUTS
!!     - matrix: the matrix to output
!!     - unit_nr: the unit where to output
!!     - long_description: if a long description should be written
!!       (defaults to false)
!!     - local: if the unit is a local unit or a global unit 
!!       (defaults to false, i.e. global)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE cp_fm_write(matrix, unit_nr, long_description, local, error)
    TYPE(cp_full_matrix_type), POINTER :: matrix
    INTEGER, INTENT(in) :: unit_nr
    LOGICAL, OPTIONAL, INTENT(in) :: long_description, local
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure, loc, my_long_description, exists
    CHARACTER(len=*), PARAMETER :: routineN='cp_fm_write',&
         routineP=moduleN//':'//routineN
    INTEGER :: iostat, desc_unit, data_unit
    TYPE(cp_logger_type), POINTER :: logger
    TYPE(cp_para_env_type), POINTER :: para_env
    CHARACTER(len=60):: filename

    failure=.FALSE.; loc=.FALSE.; my_long_description=.FALSE.
    filename=' '

    IF (PRESENT(local)) loc=local
    IF (PRESENT(long_description)) my_long_description=long_description
    logger=>cp_error_get_logger(error) 
    para_env=>logger%para_env
    IF (ASSOCIATED(matrix)) THEN
       data_unit=unit_nr
       desc_unit=unit_nr
       IF (.not.loc.and.long_description) THEN
          matrix%print_count=matrix%print_count+1
          IF (para_env%mepos==para_env%source) THEN
             filename=TRIM(matrix%name)//'-'//&
                  ADJUSTL(cp_to_string(matrix%print_count))
             CALL compress_all(filename)
             WRITE (unit=desc_unit,&
                  fmt="(' <cp_fm_types>:{ id_nr=',i10,' written to files ',a,'*')",&
                  iostat=iostat) matrix%id_nr, TRIM(filename)
          END IF
          CALL cp_logger_generate_filename(logger=logger,res=filename,&
               root=TRIM(matrix%name)//"-"//&
               TRIM(ADJUSTL(cp_to_string(matrix%print_count))),&
               postfix=".desc",local=.TRUE.)
          INQUIRE (FILE=TRIM(filename),EXIST=exists)
          CPPreconditionNoFail(.not.exists,cp_warning_level,routineP,error)
          CALL open_file(TRIM(filename),file_status="unknown",&
               file_action="WRITE",&
               unit_number=desc_unit)

          CALL cp_logger_generate_filename(logger=logger,res=filename,&
            root=TRIM(matrix%name)//"-"//&
            TRIM(ADJUSTL(cp_to_string(matrix%print_count))),&
            postfix=".dat",local=.TRUE.)
          INQUIRE (FILE=TRIM(filename),EXIST=exists)
          CPPreconditionNoFail(.not.exists,cp_warning_level,routineP,error)
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
             CALL cp_2d_r_write(matrix%local_data,unit_nr=data_unit,error=error)
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

       IF (.not.loc.and.long_description) THEN
          CALL close_file(desc_unit)
          CALL close_file(data_unit)
       END IF
    ELSE
       IF (loc .OR. para_env%mepos==para_env%source) THEN
          WRITE (unit=unit_nr,fmt="(a)",iostat=iostat) " <cp_fm_types>:*null*"
       END IF
    END IF
  END SUBROUTINE cp_fm_write

END MODULE cp_fm_types
