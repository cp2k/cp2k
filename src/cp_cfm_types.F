!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2002  CP2K developers group
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_cfm_types [1.0] *
!!
!!   NAME
!!     cp_cfm_types
!!
!!   FUNCTION
!!     represent a complex full matrix distribued on many processors
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Joost VandeVondele, based on Fawzi's cp_fm_* routines
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!****************************************************************************
MODULE cp_cfm_types
  USE cp_array_utils,                  ONLY: cp_2d_r_output,&
                                             cp_2d_r_write
  USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                             cp_blacs_env_from_globenv,&
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
  USE sparse_matrix_types,             ONLY: first_block_node,&
                                             get_block_node,&
                                             get_matrix_info,&
                                             next_block_node,&
                                             real_block_node_type,&
                                             real_matrix_type,&
                                             sparse_plus_rep_rept,&
                                             sparse_times_replicated
  USE string_utilities,                ONLY: compress_all
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_cfm_types'

  INTEGER,SAVE,PRIVATE :: last_cfm_id_nr=0


  PUBLIC :: cp_complex_full_matrix_type, cp_complex_full_matrix_p_type
  PUBLIC :: cp_cfm_create, cp_cfm_create2, cp_cfm_retain, cp_cfm_release,&
       cp_cfm_init_random, cp_cfm_set_element, cp_cfm_set_all,&
       cp_cfm_get_info, cp_cfm_get_element
!***
!****************************************************************************

!!****s* cp_fm_types/cp_complex_full_matrix_type [1.0] *
!!
!!   NAME
!!     cp_complex_full_matrix_type
!!
!!   FUNCTION
!!     represent a complex full matrix
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
!!
!!   SOURCE
  !***************************************************************************
  TYPE cp_complex_full_matrix_type
!    PRIVATE
     CHARACTER(LEN=60) :: name
     INTEGER :: id_nr, ref_count, print_count
     TYPE(cp_fm_struct_type), POINTER :: matrix_struct
     COMPLEX(kind=wp), DIMENSION(:,:), POINTER :: local_data
  END TYPE cp_complex_full_matrix_type
!!***
  !****************************************************************************

!!****s* cp_cfm_types/cp_complex_full_matrix_p_type [1.0] *
!!
!!   NAME
!!     cp_complex_full_matrix_p_type
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
!!   SOURCE
  !***************************************************************************
  TYPE cp_complex_full_matrix_p_type
     TYPE(cp_complex_full_matrix_type), POINTER :: matrix
  END TYPE cp_complex_full_matrix_p_type
!!***
  !****************************************************************************

CONTAINS

  SUBROUTINE cp_cfm_create2(new_matrix,nrow_global,ncol_global,&
       nrow_block,ncol_block,name,globenv,&
       context,para_env,error)

!   ***************************************************************************

    TYPE(cp_complex_full_matrix_type), POINTER        :: new_matrix
    CHARACTER(LEN=*), INTENT(IN)              :: name
    INTEGER, INTENT(IN)                       :: ncol_block,&
         ncol_global,nrow_block,&
         nrow_global
    TYPE(global_environment_type), INTENT(IN), OPTIONAL :: globenv
    TYPE(cp_blacs_env_type), POINTER, OPTIONAL          :: context
    TYPE(cp_para_env_type), POINTER, OPTIONAL           :: para_env
    TYPE(cp_error_type), INTENT(inout),OPTIONAL         :: error

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routineN ="cp_cfm_create2",&
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
    CALL cp_cfm_create(new_matrix,matrix_struct=matrix_struct)
    CALL cp_fm_struct_release(matrix_struct)
  END SUBROUTINE cp_cfm_create2

!!****f* cp_fm_types/cp_cfm_create [1.0] *
!!
!!   NAME
!!     cp_cfm_create
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
!!
!!*** **********************************************************************
  SUBROUTINE cp_cfm_create(matrix,matrix_struct,name,error)
    TYPE(cp_complex_full_matrix_type), POINTER :: matrix
    TYPE(cp_fm_struct_type), POINTER :: matrix_struct
    CHARACTER(len=*), INTENT(in), OPTIONAL :: name
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='cp_cfm_create',&
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
       last_cfm_id_nr=last_cfm_id_nr+1
       matrix%id_nr=last_cfm_id_nr
       matrix%ref_count=1
       matrix%print_count=0

       nprow=context%num_pe(1)
       npcol=context%num_pe(2)
       NULLIFY(matrix%local_data)

       nrow_local=matrix_struct%nrow_local(context%mepos(1))
       ncol_local=matrix_struct%ncol_local(context%mepos(2))
       ALLOCATE(matrix%local_data(nrow_local,ncol_local),stat=stat)

       CALL dcopy(2*nrow_local*ncol_local,0.0_wp,0,matrix%local_data,1)
! matrix%local_data(:,:) = 0.0_wp

       IF (PRESENT(name)) THEN
          matrix%name=name
       ELSE
          matrix%name='full matrix'//cp_to_string(matrix%id_nr)
       END IF
    END IF
  END SUBROUTINE cp_cfm_create
!***************************************************************************

!!****f* cp_fm_types/cp_cfm_retain [1.0] *
!!
!!   NAME
!!     cp_cfm_retain
!!
!!   SYNOPSIS
!!     Subroutine cp_cfm_retain(matrix, error)
!!       Type(cp_complex_full_matrix_type), Pointer:: matrix
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_cfm_retain
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
!!*** **********************************************************************
  SUBROUTINE cp_cfm_retain(matrix,error)
    TYPE(cp_complex_full_matrix_type), POINTER :: matrix
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='cp_cfm_retain',&
         routineP=moduleN//':'//routineN
    failure=.FALSE.

    CPPrecondition(ASSOCIATED(matrix),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPrecondition(matrix%ref_count>0,cp_failure_level,routineP,error,failure)
       matrix%ref_count=matrix%ref_count+1
    END IF
  END SUBROUTINE cp_cfm_retain
!***************************************************************************

!!****f* cp_cfm_types/cp_cfm_release [1.0] *
!!
!!   NAME
!!     cp_fm_release
!!
!!   SYNOPSIS
!!     Subroutine cp_cfm_release(matrix, error)
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
!!
!!*** **********************************************************************
  SUBROUTINE cp_cfm_release(matrix,error)
    TYPE(cp_complex_full_matrix_type), POINTER :: matrix
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='cp_cfm_release',&
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
  END SUBROUTINE cp_cfm_release
!***************************************************************************

! *****************************************************************************
! init ncol vectors of a full matrix with random numbers
! *****************************************************************************
  SUBROUTINE cp_cfm_init_random(matrix,ncol,error)

    TYPE(cp_complex_full_matrix_type), POINTER           :: matrix
    INTEGER, INTENT(IN)                          :: ncol
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

    INTEGER :: handle,n,i
    COMPLEX(wp), DIMENSION(:,:), POINTER :: a
    INTEGER :: nprow,npcol,myprow,mypcol,nrow_block,ncol_block
    INTEGER, DIMENSION(9) :: desca
    INTEGER, DIMENSION(4),SAVE :: ISEED
    LOGICAL, SAVE :: FIRST = .TRUE.
    INTEGER :: irow_global,icol_global,irow_local,icol_local,iprow,ipcol
    TYPE(cp_blacs_env_type), POINTER :: context

    CALL timeset("cp_cfm_init_random","I","",handle)

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

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE CP_cFM_INIT_RANDOM
! *****************************************************************************

  SUBROUTINE cp_cfm_set_all(matrix,alpha,error)

!   Purpose: Set the FULL matrix elements to alpha.

!   History: - Creation (12.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(cp_complex_full_matrix_type), POINTER           :: matrix
    COMPLEX(wp), INTENT(IN)                         :: alpha
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   *** Local variables ***

    INTEGER :: mypcol,myprow
    TYPE(cp_blacs_env_type), POINTER :: context

!   ---------------------------------------------------------------------------

    context => matrix%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)

    matrix%local_data(:,:) = alpha

  END SUBROUTINE cp_cfm_set_all

! *****************************************************************************

!!****f* cp2k/m/s *
!!
!!   NAME
!!     s
!!
!!   FUNCTION
!!     -
!!
!!   ARGUMENTS
!!     -
!!
!!   NOTES
!!     -
!!
!!   MODIFICATION HISTORY
!!     , TCH, created
!!
!!   SEE ALSO
!!     
!!
!!****
  SUBROUTINE cp_cfm_get_element(matrix,irow_global,icol_global,alpha)
    
    IMPLICIT NONE

    ! arguments
    TYPE(cp_complex_full_matrix_type), POINTER          :: matrix
    COMPLEX(wp), INTENT(OUT)                     :: alpha
    INTEGER, INTENT(IN)                       :: icol_global,&
                                                 irow_global

    ! locals
    INTEGER :: icol_local,ipcol,iprow,irow_local,mypcol,myprow,npcol,&
               nprow
    INTEGER, DIMENSION(9) :: desca
    COMPLEX(wp), DIMENSION(:,:), POINTER :: a
    TYPE(cp_blacs_env_type), POINTER :: context

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
      alpha = a(irow_local,icol_local)
    ELSE
      alpha = CMPLX(0.0_wp,0.0_wp)
    END IF

#else
    alpha = matrix%local_data(irow_global,icol_global)
#endif

  END SUBROUTINE cp_cfm_get_element

! *****************************************************************************

  SUBROUTINE cp_cfm_set_element(matrix,irow_global,icol_global,alpha,error)

!   Purpose: Set the FULL matrix element (irow_global,icol_global) to alpha.


!   ***************************************************************************

    TYPE(cp_complex_full_matrix_type), POINTER          :: matrix
    COMPLEX(wp), INTENT(IN)                      :: alpha
    INTEGER, INTENT(IN)                       :: icol_global,&
         irow_global
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   *** Local variables ***

    INTEGER :: icol_local,ipcol,iprow,irow_local,mypcol,myprow,npcol,&
         nprow

    INTEGER, DIMENSION(9) :: desca

    COMPLEX(wp), DIMENSION(:,:), POINTER :: a
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
  END SUBROUTINE cp_cfm_set_element

  SUBROUTINE cp_cfm_get_info(matrix,name,nrow_global,ncol_global,&
       nrow_block,ncol_block,nrow_local,ncol_local,&
       row_indices,col_indices,local_data,context,&
       matrix_struct,para_env,error)

!   Purpose: Return informations about the specified FULL matrix.
!   the indices vectors are cached.
!   indices is an array that maps the local indices to the global ones
!   local data is a little hole to access the blacs block, needed in some rare cases
!   for efficient use
!   ***************************************************************************

    TYPE(cp_complex_full_matrix_type), POINTER        :: matrix
    CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: name
    INTEGER, OPTIONAL, INTENT(OUT)          :: ncol_block,ncol_global,&
         nrow_block,nrow_global,&
         nrow_local,ncol_local
    INTEGER, OPTIONAL, DIMENSION(:), POINTER   :: row_indices,col_indices
    TYPE(cp_para_env_type), POINTER, OPTIONAL :: para_env
    TYPE(cp_blacs_env_type), POINTER, OPTIONAL :: context
    TYPE(cp_fm_struct_type),POINTER,OPTIONAL :: matrix_struct
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error
    COMPLEX(wp), DIMENSION(:,:), POINTER, OPTIONAL :: local_data

    CHARACTER(len=*), PARAMETER :: routineN='cp_cfm_get_info',&
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
       nrow_local = matrix%matrix_struct%nrow_local(myprow)
    ENDIF

    IF (PRESENT(ncol_local)) THEN
       ncol_local = matrix%matrix_struct%ncol_local(mypcol)
    ENDIF

    IF (PRESENT(row_indices)) THEN
       row_indices => matrix%matrix_struct%row_indices
       IF (.NOT.ASSOCIATED(row_indices)) THEN
          ALLOCATE(matrix%matrix_struct%row_indices &
               (matrix%matrix_struct%nrow_local(myprow)),stat=stat)
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
               (matrix%matrix_struct%ncol_local(mypcol)),stat=stat)
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
END MODULE cp_cfm_types
