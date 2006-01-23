!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001,2002,2003  CP2K developers group
!-----------------------------------------------------------------------------!


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
!!     3) separated structure object, removed globenv, renamed to full matrix
!!        many changes (fawzi 08.2002)
!!
!!   SOURCE
!****************************************************************************
MODULE cp_fm_types
  USE cp_array_utils,                  ONLY: cp_2d_r_write
  USE cp_blacs_calls,                  ONLY: cp_blacs_dgebr2d,&
                                             cp_blacs_dgebs2d
  USE cp_fm_struct,                    ONLY: cp_fm_struct_equivalent,&
                                             cp_fm_struct_get,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_retain,&
                                             cp_fm_struct_type,&
                                             cp_fm_struct_write
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: cp2k_is_parallel,&
                                             mp_bcast,&
                                             mp_max,&
                                             mp_sum
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE string_utilities,                ONLY: compress
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_fm_types'

  INTEGER,SAVE,PRIVATE :: last_fm_id_nr=0



  PUBLIC :: cp_fm_type, &
            cp_fm_p_type

  PUBLIC :: cp_fm_create, &
            cp_fm_retain, &
            cp_fm_release,&
            cp_fm_get_info, &
            cp_fm_set_element, &
            cp_fm_get_element, &
            cp_fm_set_all, &           ! set all elements and diagonal
            cp_fm_set_submatrix, &     ! set a submatrix to given values
            cp_fm_get_submatrix, &     ! get a submatrix of given values
            cp_fm_init_random, &         
            cp_fm_maxabsval, &         ! find the maximum absolute value
            cp_fm_maxabsrownorm, &     ! find the maximum of the sum of the abs of the elements of a row
            cp_fm_to_fm,&             ! copy (parts of) a fm to a fm
            cp_fm_vectorsnorm          ! compute the norm of the column-vectors      


  INTERFACE cp_fm_to_fm
    MODULE PROCEDURE cp_fm_to_fm_matrix, &  ! a full matrix
                     cp_fm_to_fm_columns    ! just a number of columns
  END INTERFACE 

!****************************************************************************

!!****s* cp_fm_types/cp_fm_type [1.0] *
!!
!!   NAME
!!     cp_fm_type
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
  TYPE cp_fm_type
!    PRIVATE
     CHARACTER(LEN=60) :: name
     INTEGER :: id_nr, ref_count, print_count
     TYPE(cp_fm_struct_type), POINTER :: matrix_struct
     REAL(KIND = dp), DIMENSION(:,:), POINTER :: local_data
  END TYPE cp_fm_type
!!***
  !****************************************************************************

!!****s* cp_fm_types/cp_fm_p_type [1.0] *
!!
!!   NAME
!!     cp_fm_p_type
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
  TYPE cp_fm_p_type
     TYPE(cp_fm_type), POINTER :: matrix
  END TYPE cp_fm_p_type
!!***
  !****************************************************************************

CONTAINS

!!****f* cp_fm_types/cp_fm_create [1.0] *
!!
!!   NAME
!!     cp_fm_create
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
    TYPE(cp_fm_type), POINTER                :: matrix
    TYPE(cp_fm_struct_type), POINTER         :: matrix_struct
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: name
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ncol_local, npcol, nprow, &
                                                nrow_local, stat
    LOGICAL                                  :: failure
    TYPE(cp_blacs_env_type), POINTER         :: context

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
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       ! JVDV we should remove this, as it is up to the user to zero afterwards
       CALL dcopy(nrow_local*ncol_local,0.0_dp,0,matrix%local_data,1)

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
    TYPE(cp_fm_type), POINTER                :: matrix
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
!***************************************************************************

!!****f* cp_fm_types/cp_fm_release [1.0] *
!!
!!   NAME
!!     cp_fm_release
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
    TYPE(cp_fm_type), POINTER                :: matrix
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(matrix)) THEN
       CPPrecondition(matrix%ref_count>0,cp_failure_level,routineP,error,failure)
       matrix%ref_count=matrix%ref_count-1
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

!!****f* cp_fm_types/cp_fm_init_random [1.0] *
!!
!!   NAME
!!     cp_fm_init_random
!!
!!   FUNCTION
!!     fills a matrix with random numbers
!!
!!   NOTES
!!     the value of a_ij is independent of the number of cpus 
!!
!!   INPUTS
!!     matrix :: to be initialized
!!     ncol :: numbers of cols to fill
!!     start_col :: starting at coll number
!!
!!   AUTHOR
!!     Joost VandeVondele
!!
!!   MODIFICATION HISTORY
!!
!!*** *********************************************************************
  SUBROUTINE cp_fm_init_random(matrix,ncol,start_col,error)
    TYPE(cp_fm_type), POINTER                :: matrix
    INTEGER, INTENT(IN), OPTIONAL            :: ncol, start_col
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_init_random', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, icol, my_ncol, &
                                                my_start_col, n
    INTEGER, DIMENSION(4), SAVE              :: ISEED
    LOGICAL                                  :: failure
    LOGICAL, SAVE                            :: FIRST = .TRUE.
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: buff

    CALL timeset("cp_fm_init_random",handle)
    IF (FIRST) THEN
       ISEED(1)=4
       ISEED(2)=3
       ISEED(3)=2
       ISEED(4)=1 ! last one has to be odd
       FIRST=.FALSE.
    ENDIF

    my_start_col = 1
    IF (PRESENT(start_col)) my_start_col=start_col
    my_ncol = matrix%matrix_struct%ncol_global
    IF (PRESENT(ncol)) my_ncol=ncol
    n = matrix%matrix_struct%nrow_global
    call cp_assert(matrix%matrix_struct%ncol_global>=(my_start_col+my_ncol-1),&
         cp_failure_level,cp_assertion_failed,routineP,&
         "matrix%matrix_struct%ncol_global>=(my_start_col+my_ncol-1)",&
         error,failure)
    ALLOCATE(buff(n,1))
    DO icol=my_start_col, my_start_col+my_ncol-1
       IF (matrix%matrix_struct%para_env%mepos.eq.0) THEN
          CALL dlarnv(1,iseed,n,buff(1,1))
       ENDIF
       CALL mp_bcast(buff,0,matrix%matrix_struct%para_env%group) 
       ! write(6,*) buff
       CALL cp_fm_set_submatrix(matrix,buff,start_col=icol)

    ENDDO
    ! write(6,*) matrix%local_data
    DEALLOCATE(buff)

    CALL timestop(handle)

  END SUBROUTINE cp_fm_init_random

!!****f* cp_fm_types/cp_fm_set_all [1.0] *
!!
!!   NAME
!!     cp_fm_set_all
!!
!!   FUNCTION
!!     set all elements of a matrix to the same value, 
!!     and optionally the diagonal to a different one
!!
!!   NOTES
!!     can be used to zero a matrix
!!     can be used to create a unit matrix (I-matrix) alpha=0.0_dp beta=1.0_dp
!!
!!   INPUTS
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!*** *********************************************************************

  SUBROUTINE cp_fm_set_all(matrix,alpha,beta,error)
    TYPE(cp_fm_type), POINTER                :: matrix
    REAL(KIND=dp), INTENT(IN)                :: alpha
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: beta
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    INTEGER                                  :: i, matrix_size

!   ---------------------------------------------------------------------------

    matrix_size = SIZE(matrix%local_data,1)*SIZE(matrix%local_data,2)
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
!!     returns an element of a fm 
!!     this value is valid on every cpu
!!     using this call is expensive
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
    TYPE(cp_fm_type), POINTER          :: matrix
    REAL(KIND = dp), INTENT(OUT)                     :: alpha
    INTEGER, INTENT(IN)                       :: icol_global,&
                                                 irow_global
    LOGICAL, INTENT(OUT), OPTIONAL            :: local

    ! locals
#if defined(__SCALAPACK)
    INTEGER, DIMENSION(9) :: desca
    TYPE(cp_blacs_env_type), POINTER :: context
    INTEGER :: icol_local,ipcol,iprow,irow_local,mypcol,myprow,npcol,&
               nprow
    REAL(KIND = dp), DIMENSION(:,:), POINTER :: a
#endif

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

!!****f* cp_fm_types/cp_fm_set_element [1.0] *
!!
!!   NAME
!!     cp_fm_set_element
!!
!!   FUNCTION
!!     sets an element of a matrix
!!
!!   NOTES
!!     we expect all cpus to have the same arguments in the call to this function
!!     (otherwise one should use local_data tricks)
!!
!!   INPUTS
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!*** *********************************************************************
  SUBROUTINE cp_fm_set_element(matrix,irow_global,icol_global,alpha,error)
    TYPE(cp_fm_type), POINTER                :: matrix
    INTEGER, INTENT(IN)                      :: irow_global, icol_global
    REAL(KIND=dp), INTENT(IN)                :: alpha
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    INTEGER                                  :: icol_local, ipcol, iprow, &
                                                irow_local, mypcol, myprow, &
                                                npcol, nprow
    INTEGER, DIMENSION(9)                    :: desca
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a
    TYPE(cp_blacs_env_type), POINTER         :: context

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
!!   FUNCTION
!!     sets a submatrix of a full matrix
!!      fm(start_row:start_row+n_rows,start_col:start_col+n_cols)
!!      = alpha*op(new_values)(1:n_rows,1:n_cols)+ beta 
!!      * fm(start_row:start_row+n_rows,start_col:start_col+n_cols)
!!
!!   NOTES
!!     optimized for full column updates and alpha=1.0, beta=0.0
!!     the new_values need to be valid on all cpus
!!
!!   INPUTS
!!     - new_values: a replicated full matrix with the new values
!!     - fm: the full to change
!!     - start_row: the starting row of b_matrix (defaults to 1)
!!     - start_col: the starting col of b_matrix (defaults to 1)
!!     - n_rows: the number of row to change in b (defaults to
!!       size(op(new_values),1))
!!     - n_cols: the number of columns to change in b (defaults to
!!       size(op(new_values),2))
!!     - alpha: rescaling factor for the new values (defaults to 1.0)
!!     - beta: rescaling factor for the old values (defaults to 0.0)
!!     - transpose: if new_values should be transposed: if true
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
  SUBROUTINE cp_fm_set_submatrix(fm,new_values,start_row,&
       start_col, n_rows, n_cols, alpha, beta, transpose, error)
    TYPE(cp_fm_type), POINTER                :: fm
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(in)                             :: new_values
    INTEGER, INTENT(in), OPTIONAL            :: start_row, start_col, n_rows, &
                                                n_cols
    REAL(KIND=dp), INTENT(in), OPTIONAL      :: alpha, beta
    LOGICAL, INTENT(in), OPTIONAL            :: transpose
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
    ! CALL timeset(routineN//','//moduleN,"I","",handle)

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(fm),cp_failure_level,routineP,error,failure)
    IF (.not.failure) THEN
       CPPrecondition(fm%ref_count>0,cp_failure_level,routineP,error,failure)
    END IF
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

    ! CALL timestop(0.0_dp,handle)

  END SUBROUTINE cp_fm_set_submatrix

!!****f* cp_fm_types/cp_fm_get_submatrix [1.0] *
!!
!!   NAME
!!     cp_fm_get_submatrix
!!
!!   FUNCTION
!!     gets a submatrix of a full matrix
!!      op(target_m)(1:n_rows,1:n_cols)
!!      =fm(start_row:start_row+n_rows,start_col:start_col+n_cols)
!!     target_m is replicated on all cpus
!!     using this call is expensive
!!
!!   NOTES
!!     optimized for full column updates. Zeros out a little too much
!!     of target_m
!!     the target_m is replicated and valid on all cpus
!!
!!   INPUTS
!!     - target_m: a replicated full matrix that will contain the result
!!     - fm: the full you want to get the info from
!!     - start_row: the starting row of b_matrix (defaults to 1)
!!     - start_col: the starting col of b_matrix (defaults to 1)
!!     - n_rows: the number of row to change in b (defaults to
!!       size(op(new_values),1))
!!     - n_cols: the number of columns to change in b (defaults to
!!       size(op(new_values),2))
!!     - transpose: if target_m should be transposed: if true
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
  SUBROUTINE cp_fm_get_submatrix(fm,target_m, start_row,&
       start_col, n_rows, n_cols, transpose, error)
    TYPE(cp_fm_type), POINTER                :: fm
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(out)                            :: target_m
    INTEGER, INTENT(in), OPTIONAL            :: start_row, start_col, n_rows, &
                                                n_cols
    LOGICAL, INTENT(in), OPTIONAL            :: transpose
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_get_submatrix', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, i0, j, j0, ncol, ncol_global, &
      ncol_local, nrow, nrow_global, nrow_local, this_col, &
      this_row
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
       CALL dcopy(SIZE(target_m,1)*SIZE(target_m,2),0.0_dp,0,target_m(1,1),1)
#endif

       CALL cp_fm_get_info(matrix=fm, &
            nrow_global=nrow_global,ncol_global=ncol_global, &
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

    CALL timestop(handle)

  END SUBROUTINE cp_fm_get_submatrix


!!****f* cp_fm_types/cp_fm_get_info [1.0] *
!!
!!   NAME
!!     cp_fm_get_info
!!
!!   FUNCTION
!!      returns all kind of information about the full matrix
!!
!!   NOTES
!!      see also cp_fm_struct for explaination
!!      - nrow_local, ncol_local, row_indices, col_indices, local_data are hooks for efficient 
!!        access to the local blacs block 
!!
!!   INPUTS
!!      matrix_a
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!*** *********************************************************************
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
    TYPE(cp_error_type),INTENT(inout),OPTIONAL :: error
    REAL(KIND = dp), DIMENSION(:,:),POINTER, OPTIONAL :: local_data

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


!!****f* cp_fm_types/cp_fm_maxabsval [1.0] *
!!
!!   NAME
!!     cp_fm_maxabsval
!!
!!   FUNCTION
!!     find the maximum absolute value of the matrix element
!!     maxval(abs(matrix))
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!*** *********************************************************************
  SUBROUTINE cp_fm_maxabsval(matrix,a_max,ir_max,ic_max,error)
    TYPE(cp_fm_type), POINTER                :: matrix
    REAL(KIND=dp), INTENT(OUT)               :: a_max
    INTEGER, INTENT(OUT), OPTIONAL           :: ir_max,ic_max
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_maxabsval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ic_max_local, &
                                                ir_max_local, istat, j, mepos, &
                                                ncol_local, nrow_local, num_pe
    INTEGER, DIMENSION(:), ALLOCATABLE       :: ir_max_vec, ic_max_vec
    INTEGER, DIMENSION(:), POINTER           :: row_indices, col_indices
    LOGICAL :: failure
    REAL(dp)                                 :: my_max
    REAL(dp), DIMENSION(:), ALLOCATABLE      :: a_max_vec
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: my_block

!   ---------------------------------------------------------------------------

    CALL timeset("cp_fm_maxabsval",handle)

    failure = .FALSE.
    my_block => matrix%local_data

    CALL cp_fm_get_info(matrix, nrow_local=nrow_local, ncol_local=ncol_local,&
                        row_indices=row_indices,col_indices=col_indices,error=error)
    a_max = MAXVAL(ABS(my_block(1:nrow_local,1:ncol_local)))

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
            IF( ABS(my_block(j,i)) > my_max) THEN
              my_max =  my_block(j,i)
              ir_max_local = j
              ic_max_local = i
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

!!****f* cp_fm_types/cp_fm_maxabsrownorm [1.0] *
!!
!!   NAME
!!     cp_fm_maxabsrownorm
!!
!!   FUNCTION
!!     find the maximum over the rows of the sum of the absolute values of the elements of a given row 
!!     = || A ||_infinity
!!
!!   NOTES
!!     for a real symmetric matrix it holds that || A ||_2 = |lambda_max| < || A ||_infinity
!!     Hence this can be used to estimate an upper bound for the eigenvalues of a matrix
!!     http://mathworld.wolfram.com/MatrixNorm.html
!!     (but the bound is not so tight in the general case)
!!
!!   INPUTS
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!*** *********************************************************************
  SUBROUTINE cp_fm_maxabsrownorm(matrix, a_max)
    TYPE(cp_fm_type), POINTER                :: matrix
    REAL(KIND=dp), INTENT(OUT)               :: a_max

    INTEGER                                  :: handle, i, j, ncol_local, &
                                                nrow_global, nrow_local
    INTEGER, DIMENSION(:), POINTER           :: row_indices
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: values
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: my_block

!   ---------------------------------------------------------------------------

    CALL timeset("cp_fm_maxabsrownorm",handle)

    my_block => matrix%local_data
    CALL cp_fm_get_info(matrix, row_indices=row_indices, nrow_global=nrow_global, &
                                nrow_local=nrow_local, ncol_local=ncol_local) 

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

!!****f* cp_fm_types/cp_fm_vectorsnorm [1.0] *
!!
!!   NAME
!!     cp_fm_vectorsnorm
!!
!!   FUNCTION
!!     find the inorm of each column norm_{j}= sqrt( \sum_{i} A_{ij}*A_{ij} )
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!*** *********************************************************************
  SUBROUTINE cp_fm_vectorsnorm(matrix, norm_array)
    TYPE(cp_fm_type), POINTER                :: matrix
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: norm_array

    INTEGER                                  :: handle, i, j, ncol_local, &
                                                ncol_global, nrow_local
    INTEGER, DIMENSION(:), POINTER           :: col_indices
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: values
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: my_block

!   ---------------------------------------------------------------------------

    CALL timeset("cp_fm_maxabsrownorm",handle)

    my_block => matrix%local_data
    CALL cp_fm_get_info(matrix, col_indices=col_indices, ncol_global=ncol_global, &
                                nrow_local=nrow_local, ncol_local=ncol_local)

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


!!****f* cp_fm_types/cp_fm_to_fm_matrix [1.0] *
!!
!!   NAME
!!     cp_fm_to_fm_matrix  ! interface cp_fm_to_fm
!!
!!   SYNOPSIS
!!     SUBROUTINE cp_fm_to_fm_matrix(source,destination,error)
!!       TYPE(cp_fm_type), POINTER :: source,destination
!!       TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error
!!     END SUBROUTINE
!!
!!   FUNCTION
!!     copy one identically sized matrix in the other
!!
!!   NOTES
!!     see also cp_fm_to_fm_columns
!!
!!   INPUTS
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!*** *********************************************************************
  SUBROUTINE cp_fm_to_fm_matrix(source,destination,error)
    TYPE(cp_fm_type), POINTER                :: source, destination
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE cp_fm_to_fm (MODULE cp_fm_types)"

    INTEGER                                  :: npcol, nprow
    LOGICAL                                  :: failure

!   ---------------------------------------------------------------------------

    failure=.TRUE.
    nprow = source%matrix_struct%context%num_pe(1)
    npcol = source%matrix_struct%context%num_pe(2)

    IF(.NOT. cp2k_is_parallel .OR.&
         cp_fm_struct_equivalent(source%matrix_struct,&
         destination%matrix_struct,error=error)) THEN
      CALL cp_assert(SIZE(source%local_data,1)==SIZE(destination%local_data,1).AND.&
         SIZE(source%local_data,2)==SIZE(destination%local_data,2),cp_failure_level,&
         cp_assertion_failed, routine, &
         "internal error, local_data has different sizes"//&
         CPSourceFileRef,&
         error=error)
      CALL dcopy(SIZE(source%local_data,1)*SIZE(source%local_data,2),&
           source%local_data(1,1),1,destination%local_data(1,1),1)
   ELSE
      CALL cp_assert(source%matrix_struct%nrow_global==&
           destination%matrix_struct%nrow_global,&
           cp_failure_level,&
           cp_assertion_failed, routine, &
           "cannot copy between full matrixes of differen sizes"//&
           CPSourceFileRef,&
           error=error)
      CALL cp_assert(source%matrix_struct%ncol_global==&
           destination%matrix_struct%ncol_global,&
           cp_failure_level,&
           cp_assertion_failed, routine, &
           "cannot copy between full matrixes of differen sizes"//&
           CPSourceFileRef,&
           error=error)
#ifdef __SCALAPACK
      CALL pdcopy(source%matrix_struct%nrow_global*&
           source%matrix_struct%ncol_global,&
           source%local_data(1,1),1,1,source%matrix_struct%descriptor,1,&
           destination%local_data(1,1),1,1,destination%matrix_struct%descriptor,1)
#else
      CPAssert(.FALSE.,cp_failure_level,routine,error,failure)
#endif
   END IF
  END SUBROUTINE cp_fm_to_fm_matrix

!!****f* cp_fm_types/cp_fm_to_fm_columns [1.0] *
!!
!!   NAME
!!     cp_fm_to_fm_columns ! interface cp_fm_to_fm
!!
!!   SYNOPSIS
!!     SUBROUTINE cp_fm_to_fm_columns(msource,mtarget,ncol,source_start,target_start)
!!       TYPE(cp_fm_type), POINTER                 :: msource,mtarget
!!       INTEGER, INTENT(IN)                       :: ncol
!!       INTEGER, INTENT(IN), OPTIONAL             :: source_start,target_start
!!     END SUBROUTINE
!!
!!   FUNCTION
!!     copy just a subset of columns of a fm to a fm
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!*** *********************************************************************
  SUBROUTINE cp_fm_to_fm_columns(msource,mtarget,ncol,source_start,target_start)
    TYPE(cp_fm_type), POINTER                :: msource, mtarget
    INTEGER, INTENT(IN)                      :: ncol
    INTEGER, INTENT(IN), OPTIONAL            :: source_start, target_start

    INTEGER                                  :: i, mypcol, myprow, n, ss, ts
    INTEGER, DIMENSION(9)                    :: desca, descb
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a, b
    TYPE(cp_blacs_env_type), POINTER         :: context

    ss=1
    ts=1
    IF (PRESENT(source_start)) ss=source_start
    IF (PRESENT(target_start)) ts=target_start

    context => msource%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    n = msource%matrix_struct%nrow_global

    a => msource%local_data
    b => mtarget%local_data

#if defined(__SCALAPACK)
    desca(:) = msource%matrix_struct%descriptor(:)
    descb(:) = mtarget%matrix_struct%descriptor(:)
    DO i=0,ncol-1
      CALL pdcopy(n,a(1,1),1,ss+i,desca,1,b(1,1),1,ts+i,descb,1)
    ENDDO
#else
    CALL dcopy(ncol*n,a(1,ss),1,b(1,ts),1)
#endif

  END SUBROUTINE cp_fm_to_fm_columns

END MODULE cp_fm_types
