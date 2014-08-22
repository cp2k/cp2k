!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief represent a complex full matrix distribued on many processors
!> \author Joost VandeVondele, based on Fawzi's cp_fm_* routines
! *****************************************************************************
MODULE cp_cfm_types

  USE cp_blacs_calls,                  ONLY: cp_blacs_zgebr2d,&
                                             cp_blacs_zgebs2d
  USE cp_blacs_env,                    ONLY: cp_blacs_env_type
  USE cp_fm_struct,                    ONLY: cp_fm_struct_equivalent,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_retain,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: cp2k_is_parallel
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_cfm_types'

  INTEGER,SAVE,PRIVATE :: last_cfm_id_nr=0

  PUBLIC :: cp_cfm_type, cp_cfm_p_type
  PUBLIC :: cp_cfm_create, cp_cfm_release
  PUBLIC :: cp_cfm_set_all, cp_cfm_get_info, cp_cfm_get_element
  PUBLIC :: cp_cfm_to_cfm, cp_cfm_to_fm

  INTERFACE cp_cfm_to_cfm
    MODULE PROCEDURE cp_cfm_to_cfm_matrix, &  ! a full matrix
                     cp_cfm_to_cfm_columns    ! just a number of columns
  END INTERFACE

! *****************************************************************************
!> \brief represent a complex full matrix
!> \param name the name of the matrix, used for printing
!> \param id_nr identification number that makes this matrix unique
!> \param ref_cont reference count (see doc/ReferenceCounting.html)
!> \param print_count sequential number to make the printed files unique
!> \param matrix_struct structure of this matrix
!> \param local_data array with the data of the matrix (its contents
!>        depend on the matrix type used: in parallel runs it will be
!>        in scalapack format, in sequential, it will simply contain
!>        the matrix)
! *****************************************************************************
  TYPE cp_cfm_type
!    PRIVATE
     CHARACTER(LEN=60) :: name
     INTEGER :: id_nr, ref_count, print_count
     TYPE(cp_fm_struct_type), POINTER :: matrix_struct
     COMPLEX(KIND = dp), DIMENSION(:,:), POINTER :: local_data
  END TYPE cp_cfm_type

! *****************************************************************************
!> \brief just to build arrays of pointers to matrixes
!> \param matrix the pointer to the matrix
! *****************************************************************************
  TYPE cp_cfm_p_type
     TYPE(cp_cfm_type), POINTER :: matrix
  END TYPE cp_cfm_p_type

CONTAINS

! *****************************************************************************
!> \brief creates a new full matrix with the given structure
!> \param matrix the matrix to be created
!> \param matrix_struct the structure of matrix
!> \param name ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      preferred allocation routine
! *****************************************************************************
  SUBROUTINE cp_cfm_create(matrix,matrix_struct,name,error)
    TYPE(cp_cfm_type), POINTER               :: matrix
    TYPE(cp_fm_struct_type), POINTER         :: matrix_struct
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: name
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_create', &
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
       last_cfm_id_nr=last_cfm_id_nr+1
       matrix%id_nr=last_cfm_id_nr
       matrix%ref_count=1
       matrix%print_count=0

       nprow=context%num_pe(1)
       npcol=context%num_pe(2)
       NULLIFY(matrix%local_data)

       nrow_local=matrix_struct%local_leading_dimension
       ncol_local=MAX(1,matrix_struct%ncol_locals(context%mepos(2)))
       ALLOCATE(matrix%local_data(nrow_local,ncol_local),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       CALL dcopy(2*nrow_local*ncol_local,0.0_dp,0,matrix%local_data,1)
! matrix%local_data(:,:) = 0.0_dp

       IF (PRESENT(name)) THEN
          matrix%name=name
       ELSE
          matrix%name='full matrix'//cp_to_string(matrix%id_nr)
       END IF
    END IF
  END SUBROUTINE cp_cfm_create

! *****************************************************************************
!> \brief retains a full matrix
!> \param matrix the matrix to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE cp_cfm_retain(matrix,error)
    TYPE(cp_cfm_type), POINTER               :: matrix
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(matrix),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPrecondition(matrix%ref_count>0,cp_failure_level,routineP,error,failure)
       matrix%ref_count=matrix%ref_count+1
    END IF
  END SUBROUTINE cp_cfm_retain

! *****************************************************************************
!> \brief releases a full matrix
!> \param matrix the matrix to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE cp_cfm_release(matrix,error)
    TYPE(cp_cfm_type), POINTER               :: matrix
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_release', &
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
  END SUBROUTINE cp_cfm_release

! init ncol vectors of a full matrix with random numbers
! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param ncol ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_cfm_init_random(matrix,ncol,error)

    TYPE(cp_cfm_type), POINTER               :: matrix
    INTEGER, INTENT(IN)                      :: ncol
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_init_random', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: a
    INTEGER                                  :: handle, mypcol, myprow, n, &
                                                npcol, nprow
    INTEGER, DIMENSION(4), SAVE              :: ISEED
    LOGICAL, SAVE                            :: FIRST = .TRUE.
    TYPE(cp_blacs_env_type), POINTER         :: context
#if defined(__SCALAPACK)
    INTEGER                                  :: i, icol_local, ipcol, &
                                                iprow, irow_local
    INTEGER, DIMENSION(9)                    :: desca
#endif

    CALL timeset(routineN,handle)

    context => matrix%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)
    n = matrix%matrix_struct%nrow_global

    a => matrix%local_data

#if defined(__SCALAPACK)
    IF (FIRST) THEN
       ISEED(1)=myprow
       ISEED(2)=mypcol
       ISEED(3)=2
       ISEED(4)=1 ! last one has to be odd
       FIRST=.FALSE.
    ENDIF
    desca(:)=matrix%matrix_struct%descriptor(:)
    DO i=1,ncol
       CALL infog2l(1,i,desca,nprow,npcol,myprow,mypcol,&
            irow_local,icol_local,iprow,ipcol)
       IF (ipcol == mypcol) THEN
          CALL dlarnv(1,iseed,2*SIZE(a,1),a(1,icol_local))
       ENDIF
    ENDDO
#else
    IF (FIRST) THEN
       ISEED(1)=4
       ISEED(2)=3
       ISEED(3)=2
       ISEED(4)=1 ! last one has to be odd
       FIRST=.FALSE.
    ENDIF
    CALL dlarnv(1,iseed,2*n*ncol,a(1,1))
#endif

    CALL timestop(handle)

  END SUBROUTINE CP_cFM_INIT_RANDOM

! *****************************************************************************
!> \brief   Set the FULL matrix elements to alpha.
!> \param matrix ...
!> \param alpha ...
!> \param beta ...
!> \param error ...
!> \date    12.06.2001
!> \author  Matthias Krack
!> \version 1.0
! *****************************************************************************
  SUBROUTINE cp_cfm_set_all(matrix,alpha,beta,error)
    TYPE(cp_cfm_type), POINTER               :: matrix
    COMPLEX(KIND=dp), INTENT(IN)             :: alpha
    COMPLEX(KIND=dp), INTENT(IN), OPTIONAL   :: beta
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: i

    matrix%local_data(:,:) = alpha
    IF (PRESENT(beta)) THEN
      DO i=1,matrix%matrix_struct%nrow_global
         CALL cp_cfm_set_element(matrix,i,i,beta,error=error)
      ENDDO
    ENDIF

  END SUBROUTINE cp_cfm_set_all

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param irow_global ...
!> \param icol_global ...
!> \param alpha ...
!> \par History
!>      , TCH, created
!>      always return the answer
! *****************************************************************************
  SUBROUTINE cp_cfm_get_element(matrix,irow_global,icol_global,alpha)

    TYPE(cp_cfm_type), POINTER               :: matrix
    INTEGER, INTENT(IN)                      :: irow_global, icol_global
    COMPLEX(KIND=dp), INTENT(OUT)            :: alpha

    INTEGER                                  :: mypcol, myprow, npcol, nprow
    TYPE(cp_blacs_env_type), POINTER         :: context
#if defined(__SCALAPACK)
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: a
    INTEGER                                  :: icol_local, ipcol, iprow, &
                                                irow_local
    INTEGER, DIMENSION(9)                    :: desca
#endif

    context => matrix%matrix_struct%context
    myprow = context%mepos(1)
    mypcol = context%mepos(2)
    nprow = context%num_pe(1)
    npcol = context%num_pe(2)

#if defined(__SCALAPACK)
    a => matrix%local_data
    desca(:) = matrix%matrix_struct%descriptor(:)

    CALL infog2l(irow_global,icol_global,desca,nprow,npcol,myprow,mypcol,&
                 irow_local,icol_local,iprow,ipcol)

    IF ((iprow == myprow).AND.(ipcol == mypcol)) THEN
      alpha = a(irow_local,icol_local)
      CALL cp_blacs_ZGEBS2D(context%group, 'All', ' ', 1, 1, alpha , 1 )
    ELSE
      CALL cp_blacs_ZGEBR2D(context%group, 'All', ' ', 1, 1, alpha , 1 , iprow, ipcol )
    END IF
#else
    alpha = matrix%local_data(irow_global,icol_global)
#endif

  END SUBROUTINE cp_cfm_get_element

! *****************************************************************************
!> \brief   Set the FULL matrix element (irow_global,icol_global) to alpha.
!> \param matrix ...
!> \param irow_global ...
!> \param icol_global ...
!> \param alpha ...
!> \param error ...
!> \date    12.06.2001
!> \author  Matthias Krack
!> \version 1.0
! *****************************************************************************
  SUBROUTINE cp_cfm_set_element(matrix,irow_global,icol_global,alpha,error)

    TYPE(cp_cfm_type), POINTER               :: matrix
    INTEGER, INTENT(IN)                      :: irow_global, icol_global
    COMPLEX(KIND=dp), INTENT(IN)             :: alpha
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: mypcol, myprow, npcol, nprow
    TYPE(cp_blacs_env_type), POINTER         :: context
#if defined(__SCALAPACK)
    INTEGER                                  :: icol_local, ipcol, iprow, &
                                                irow_local
    INTEGER, DIMENSION(9)                    :: desca
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: a
#endif

    context => matrix%matrix_struct%context
    myprow = context%mepos(1)
    mypcol = context%mepos(2)
    nprow = context%num_pe(1)
    npcol = context%num_pe(2)

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
  END SUBROUTINE cp_cfm_set_element

! *****************************************************************************
!> \brief   Return informations about the specified FULL matrix.
!>          the indices vectors are cached.
!>          indices is an array that maps the local indices to the global ones
!>          local data is a little hole to access the blacs block, needed in some rare cases
!>          for efficient use
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
!> \param matrix_struct ...
!> \param para_env ...
!> \param error ...
!> \date    12.06.2001
!> \author  Matthias Krack
!> \version 1.0
! *****************************************************************************
  SUBROUTINE cp_cfm_get_info(matrix,name,nrow_global,ncol_global,&
       nrow_block,ncol_block,nrow_local,ncol_local,&
       row_indices,col_indices,local_data,context,&
       matrix_struct,para_env,error)
    TYPE(cp_cfm_type), POINTER        :: matrix
    CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: name
    INTEGER, OPTIONAL, INTENT(OUT)          :: ncol_block,ncol_global,&
         nrow_block,nrow_global,&
         nrow_local,ncol_local
    INTEGER, OPTIONAL, DIMENSION(:), POINTER   :: row_indices,col_indices
    TYPE(cp_para_env_type), POINTER, OPTIONAL :: para_env
    TYPE(cp_blacs_env_type), POINTER, OPTIONAL :: context
    TYPE(cp_fm_struct_type),POINTER,OPTIONAL :: matrix_struct
    TYPE(cp_error_type), INTENT(inout)  :: error
    COMPLEX(KIND = dp), DIMENSION(:,:), POINTER, OPTIONAL :: local_data

    CHARACTER(len=*), PARAMETER :: routineN='cp_cfm_get_info',&
         routineP=moduleN//':'//routineN
    INTEGER i,nprow,npcol,myprow,mypcol, stat
    TYPE(cp_blacs_env_type), POINTER :: ctxt
    LOGICAL :: failure
#if defined(__SCALAPACK)
    INTEGER , EXTERNAL :: indxl2g
#endif

    failure=.FALSE.
    IF (PRESENT(name)) name = matrix%name
    IF (PRESENT(nrow_global)) nrow_global = matrix%matrix_struct%nrow_global
    IF (PRESENT(ncol_global)) ncol_global = matrix%matrix_struct%ncol_global
    IF (PRESENT(nrow_block)) nrow_block = matrix%matrix_struct%nrow_block
    IF (PRESENT(ncol_block)) ncol_block = matrix%matrix_struct%ncol_block
    IF (PRESENT(context)) context => matrix%matrix_struct%context
    IF (PRESENT(matrix_struct)) matrix_struct => matrix%matrix_struct
    IF (PRESENT(local_data)) local_data => matrix%local_data ! not hiding things anymore :-(

    ctxt => matrix%matrix_struct%context
    myprow=ctxt%mepos(1)
    mypcol=ctxt%mepos(2)
    nprow=ctxt%num_pe(1)
    npcol=ctxt%num_pe(2)

    IF (PRESENT(nrow_local)) THEN
       nrow_local = matrix%matrix_struct%nrow_locals(myprow)
    ENDIF

    IF (PRESENT(ncol_local)) THEN
       ncol_local = matrix%matrix_struct%ncol_locals(mypcol)
    ENDIF

    IF (PRESENT(row_indices)) THEN
       row_indices => matrix%matrix_struct%row_indices
       IF (.NOT.ASSOCIATED(row_indices)) THEN
          ALLOCATE(matrix%matrix_struct%row_indices &
               (matrix%matrix_struct%nrow_locals(myprow)),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          row_indices => matrix%matrix_struct%row_indices
#ifdef __SCALAPACK
          DO i=1,SIZE(row_indices)
             matrix%matrix_struct%row_indices(i)=&
                  indxl2g(i,matrix%matrix_struct%nrow_block,myprow,&
                  matrix%matrix_struct%first_p_pos(1),nprow)
          ENDDO
#else
          DO i=1,SIZE(row_indices)
             row_indices(i)=i
          ENDDO
#endif
       END IF
    ENDIF

    IF (PRESENT(col_indices)) THEN
       col_indices => matrix%matrix_struct%col_indices
       IF (.NOT.ASSOCIATED(col_indices)) THEN
          ALLOCATE(matrix%matrix_struct%col_indices &
               (matrix%matrix_struct%ncol_locals(mypcol)),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          col_indices => matrix%matrix_struct%col_indices
#ifdef __SCALAPACK
          DO i=1,SIZE(col_indices)
             col_indices(i)=&
                  indxl2g(i,matrix%matrix_struct%ncol_block,mypcol,&
                  matrix%matrix_struct%first_p_pos(2),npcol)
          ENDDO
#else
          DO i=1,SIZE(col_indices)
             col_indices(i)=i
          ENDDO
#endif
       END IF
    ENDIF

    IF (PRESENT(para_env)) para_env => matrix%matrix_struct%para_env

  END SUBROUTINE cp_cfm_get_info

! *****************************************************************************
!> \brief copy one identically sized matrix in the other
!> \param source ...
!> \param destination ...
!> \param error ...
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE cp_cfm_to_cfm_matrix(source,destination,error)
    TYPE(cp_cfm_type), POINTER               :: source, destination
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_to_cfm_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: npcol, nprow
    LOGICAL                                  :: failure

    failure=.TRUE.
    nprow = source%matrix_struct%context%num_pe(1)
    npcol = source%matrix_struct%context%num_pe(2)

    IF(.NOT. cp2k_is_parallel .OR.&
         cp_fm_struct_equivalent(source%matrix_struct,&
         destination%matrix_struct,error=error)) THEN
      CALL cp_assert(SIZE(source%local_data,1)==SIZE(destination%local_data,1).AND.&
         SIZE(source%local_data,2)==SIZE(destination%local_data,2),cp_failure_level,&
         cp_assertion_failed, routineP, &
         "internal error, local_data has different sizes"//&
CPSourceFileRef,&
         error=error)
      CALL dcopy(SIZE(source%local_data,1)*SIZE(source%local_data,2)*2,&
           source%local_data(1,1),1,destination%local_data(1,1),1)
   ELSE
      CALL cp_assert(source%matrix_struct%nrow_global==&
           destination%matrix_struct%nrow_global,&
           cp_failure_level,&
           cp_assertion_failed, routineP, &
           "cannot copy between full matrixes of differen sizes"//&
CPSourceFileRef,&
           error=error)
      CALL cp_assert(source%matrix_struct%ncol_global==&
           destination%matrix_struct%ncol_global,&
           cp_failure_level,&
           cp_assertion_failed, routineP, &
           "cannot copy between full matrixes of differen sizes"//&
CPSourceFileRef,&
           error=error)
#ifdef __SCALAPACK
      CALL pzcopy(source%matrix_struct%nrow_global*&
           source%matrix_struct%ncol_global,&
           source%local_data(1,1),1,1,source%matrix_struct%descriptor,1,&
           destination%local_data(1,1),1,1,destination%matrix_struct%descriptor,1)
#else
      CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
#endif
   END IF
  END SUBROUTINE cp_cfm_to_cfm_matrix

! *****************************************************************************
!> \brief copy just a subset of columns of a cfm to a cfm
!> \param msource ...
!> \param mtarget ...
!> \param ncol ...
!> \param source_start ...
!> \param target_start ...
! *****************************************************************************
  SUBROUTINE cp_cfm_to_cfm_columns(msource,mtarget,ncol,source_start,&
                                   target_start)

    TYPE(cp_cfm_type), POINTER               :: msource, mtarget
    INTEGER, INTENT(IN)                      :: ncol
    INTEGER, INTENT(IN), OPTIONAL            :: source_start, target_start

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_cfm_to_cfm_columns', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, n, ss, ts
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: a, b
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
      CALL pzcopy(n,a(1,1),1,ss+i,desca,1,b(1,1),1,ts+i,descb,1)
    END DO
#else
    CALL zcopy(ncol*n,a(1,ss),1,b(1,ts),1)
#endif

    CALL timestop(handle)

  END SUBROUTINE cp_cfm_to_cfm_columns

! *****************************************************************************
!> \brief copy real and imaginary parts of a cfm to fm
!>        works only for absolutely equivalent structures
!> \param msource ...
!> \param mtargetr ...
!> \param mtargeti ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_cfm_to_fm(msource,mtargetr,mtargeti,error)

    TYPE(cp_cfm_type), POINTER               :: msource
    TYPE(cp_fm_type), POINTER                :: mtargetr, mtargeti
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_cfm_to_fm', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: a
    INTEGER                                  :: handle
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: bi, br

    CALL timeset(routineN,handle)

    a => msource%local_data
    br => mtargetr%local_data
    bi => mtargeti%local_data

    ! This routine is only for equivalent structures
    CALL cp_assert(cp_fm_struct_equivalent(msource%matrix_struct,&
         mtargetr%matrix_struct,error=error),cp_failure_level,&
         cp_assertion_failed, routineP, &
         "internal error, local_data has different sizes"//&
CPSourceFileRef,&
         error=error)
    CALL cp_assert(cp_fm_struct_equivalent(msource%matrix_struct,&
         mtargeti%matrix_struct,error=error),cp_failure_level,&
         cp_assertion_failed, routineP, &
         "internal error, local_data has different sizes"//&
CPSourceFileRef,&
         error=error)
    CALL cp_assert(SIZE(a,1)==SIZE(br,1).AND.SIZE(a,2)==SIZE(br,2)&
         .AND.SIZE(a,1)==SIZE(bi,1).AND.SIZE(a,2)==SIZE(bi,2),cp_failure_level,&
         cp_assertion_failed, routineP, &
         "internal error, local_data has different sizes"//&
CPSourceFileRef,&
         error=error)

    br = REAL(a)
    bi = AIMAG(a)

    CALL timestop(handle)

  END SUBROUTINE cp_cfm_to_fm

END MODULE cp_cfm_types
