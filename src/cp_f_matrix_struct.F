!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_f_matrix_struct [1.0] *
!!
!!   NAME
!!     cp_f_matrix_struct
!!
!!   FUNCTION
!!     represent the structure of a full matrix
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE cp_f_matrix_struct
  USE cp_log_handling, ONLY: cp_fatal_level, cp_failure_level,&
       cp_warning_level, cp_note_level, cp_to_string, cp_log,&
       cp_logger_get_default_unit_nr, cp_logger_type
  USE cp_error_handling, ONLY: cp_debug, cp_error_type, cp_error_init,&
       cp_error_dealloc_ref, cp_error_message, cp_assert,&
       cp_assertion_failed, cp_internal_error, cp_unimplemented_error,&
       cp_error_get_logger
  USE kinds, ONLY: wp=>dp
  USE cp_para_types, ONLY: cp_para_env_type, cp_blacs_env_type
  USE cp_para_env, ONLY: cp_para_env_retain, cp_para_env_release, cp_para_env_write
  USE cp_blacs_env, ONLY: cp_blacs_env_retain, cp_blacs_env_release,&
       cp_blacs_env_write
  USE timings, ONLY: timeset, timestop
  USE message_passing,     ONLY: mp_bcast,&
                                 mp_max,&
                                 mp_sum,&
                                 mp_sync
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_f_matrix_struct'
  INTEGER,SAVE,PRIVATE :: last_fmstruct_id_nr=0

! the default blacs block sizes
! consider using #ifdefs to give them the optimal values
  INTEGER, PARAMETER :: optimal_blacs_col_block_size=32
  INTEGER, PARAMETER :: optimal_blacs_row_block_size=32

  PUBLIC :: cp_f_matrix_struct_type
  PUBLIC :: optimal_blacs_row_block_size, optimal_blacs_col_block_size
  PUBLIC :: cp_fmstruct_create,cp_fmstruct_retain, cp_fmstruct_release,&
       cp_fmstruct_equivalent, cp_fmstruct_write
!***
!****************************************************************************

!!****s* qs_blacs/cp_f_matrix_struct_type [1.0] *
!!
!!   NAME
!!     cp_f_matrix_struct_type
!!
!!   FUNCTION
!!     keeps the information about the structure of a full matrix
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - para_env:
!!     - context:
!!     - descriptor:
!!     - ncol_block:
!!     - nrow_block:
!!     - nrow_global:
!!     - ncol_block
!!     - ncol_global
!!     - first_p_pos
!!     - row_indices
!!     - col_indices
!!     - nrow_local
!!     - ncol_local
!!     - id_nr
!!     - ref_count
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!   SOURCE
!***************************************************************************
  TYPE cp_f_matrix_struct_type
     TYPE(cp_para_env_type), POINTER :: para_env
     TYPE(cp_blacs_env_type), POINTER :: context
     INTEGER, DIMENSION(9) :: descriptor
     ! when using blacs
     ! ncol_block=descriptor(6)
     ! ncol_global=descriptor(4)
     ! nrow_block=descriptor(5)
     ! nrow_global=descriptor(3)
     INTEGER :: ncol_block, nrow_block, nrow_global, ncol_global
     INTEGER, DIMENSION(2) :: first_p_pos
     INTEGER, DIMENSION(:), POINTER :: row_indices, col_indices,&
          nrow_local, ncol_local
     INTEGER :: id_nr, ref_count
  END TYPE cp_f_matrix_struct_type
!!***
!****************************************************************************

CONTAINS

!!****f* qs_blacs/cp_fmstruct_create [1.0] *
!!
!!   NAME
!!     cp_fmstruct_create
!!
!!   SYNOPSIS
!!     Subroutine cp_fmstruct_create(fmstruct, para_env, nrow_global,&
!!         ncol_global, nrow_block, ncol_block, context, descriptor,&
!!         first_p_pos, error)
!!       Type(cp_f_matrix_struct_type), Pointer:: fmstruct
!!       Type(cp_para_env_type), Pointer:: para_env
!!       Integer, Intent (IN):: nrow_global, ncol_global, nrow_block,&
!!         ncol_block
!!       Type(cp_blacs_env_type), Pointer:: context
!!       Integer, Dimension(9), Intent (IN), Optional:: descriptor
!!       Integer, Optional, Dimension(2):: first_p_pos
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_fmstruct_create
!!
!!   FUNCTION
!!     allocates and initializes a full matrix structure
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - fmstruct: the pointer that will point to the new structure
!!     - para_env: the parallel environment
!!     - nrow_global: the number of row of the full matrix
!!     - ncol_global: the number of colums of the full matrix
!!     - nrow_block: the number of rows of a block of the matrix
!!     - ncol_block: the number of colums of a block of the matrix
!!     - context: the blacs context of this matrix
!!     - descriptor: the scalapack descriptor of the matrix (if not given
!!       a new one is allocated
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
SUBROUTINE cp_fmstruct_create(fmstruct,para_env,context, nrow_global,&
     ncol_global, nrow_block, ncol_block, descriptor, first_p_pos, error)
  TYPE(cp_f_matrix_struct_type), POINTER :: fmstruct
  TYPE(cp_para_env_type), POINTER        :: para_env
  INTEGER, INTENT(in)                    :: nrow_global,ncol_global
  INTEGER, INTENT(in), OPTIONAL          :: nrow_block,ncol_block
  TYPE(cp_blacs_env_type), POINTER       :: context
  INTEGER, DIMENSION(9), INTENT(in), OPTIONAL  :: descriptor
  INTEGER,OPTIONAL, DIMENSION(2)               :: first_p_pos
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
#if defined(__SCALAPACK)
    INTEGER, EXTERNAL :: numroc
#endif

  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='cp_fmstruct_create',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat
  failure=.FALSE.
  
#if defined(__parallel) && ! defined(__SCALAPACK)
  CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
       routineP,"full matrixes need scalapack for parallel runs "//&
       CPSourceFileRef,&
       error)
#endif

  ALLOCATE(fmstruct,stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     fmstruct%context => context
     CALL cp_blacs_env_retain(fmstruct%context,error=error)
     fmstruct%para_env => para_env
     CALL cp_para_env_retain(para_env,error=error)
     
     fmstruct%nrow_global=nrow_global
     fmstruct%ncol_global=ncol_global
     ! checks for more processes than matrix elements
     CPAssert(context%num_pe(1)<=nrow_global,cp_warning_level,routineP,error,failure)
     CPAssert(context%num_pe(2)<=ncol_global,cp_warning_level,routineP,error,failure)
     fmstruct%nrow_block=optimal_blacs_row_block_size
     IF (PRESENT(nrow_block)) fmstruct%nrow_block=nrow_block
     fmstruct%ncol_block=optimal_blacs_col_block_size
     IF (PRESENT(ncol_block)) fmstruct%ncol_block=ncol_block
     fmstruct%nrow_block=MAX(1,&
          MIN(fmstruct%nrow_block,nrow_global/context%num_pe(1)))
     fmstruct%ncol_block=MAX(1,&
          MIN(fmstruct%ncol_block,ncol_global/context%num_pe(2)))
     IF (fmstruct%nrow_global==fmstruct%ncol_global) THEN
        ! square matrix -> square blocks (otherwise some op fail)
        fmstruct%nrow_block=MIN(fmstruct%nrow_block,fmstruct%ncol_block)
        fmstruct%ncol_block=fmstruct%nrow_block
     END IF
     ALLOCATE(fmstruct%nrow_local(0:(context%num_pe(1)-1)),&
          fmstruct%ncol_local(0:(context%num_pe(2)-1)),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF (.NOT.failure) THEN
     fmstruct%first_p_pos=context%source
     IF (PRESENT(first_p_pos)) fmstruct%first_p_pos=first_p_pos

     fmstruct%nrow_local=0
     fmstruct%ncol_local=0
#ifdef __SCALAPACK
     fmstruct%nrow_local(context%mepos(1))= &
          numroc(fmstruct%nrow_global,fmstruct%nrow_block,&
          fmstruct%context%mepos(1),fmstruct%first_p_pos(1),&
          fmstruct%context%num_pe(1))
     fmstruct%ncol_local(context%mepos(2))= &
          numroc(fmstruct%ncol_global,fmstruct%ncol_block,&
          fmstruct%context%mepos(2),fmstruct%first_p_pos(2),&
          fmstruct%context%num_pe(2))
     CALL mp_sum(fmstruct%nrow_local,fmstruct%para_env%group)
     CALL mp_sum(fmstruct%ncol_local,fmstruct%para_env%group)
     fmstruct%nrow_local(:)=fmstruct%nrow_local(:)/fmstruct%context%num_pe(2)
     fmstruct%ncol_local(:)=fmstruct%ncol_local(:)/fmstruct%context%num_pe(1)
#else
     ! block = full matrix
     fmstruct%nrow_block=fmstruct%nrow_global
     fmstruct%ncol_block=fmstruct%ncol_global
     fmstruct%nrow_local(context%mepos(1))=fmstruct%nrow_global
     fmstruct%ncol_local(context%mepos(2))=fmstruct%ncol_global
#endif

     NULLIFY(fmstruct%row_indices, fmstruct%col_indices)
     last_fmstruct_id_nr=last_fmstruct_id_nr+1
     fmstruct%id_nr=last_fmstruct_id_nr
     fmstruct%ref_count=1

     IF (PRESENT(descriptor)) THEN
        fmstruct%descriptor=descriptor
     ELSE
        fmstruct%descriptor=0
#if defined(__SCALAPACK)
        CALL descinit(fmstruct%descriptor,fmstruct%nrow_global,&
             fmstruct%ncol_global,fmstruct%nrow_block,&
             fmstruct%ncol_block,fmstruct%first_p_pos(1),&
             fmstruct%first_p_pos(2),fmstruct%context%group,&
             fmstruct%nrow_local(fmstruct%context%mepos(1)),stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
#endif
     END IF
  END IF
END SUBROUTINE cp_fmstruct_create
!***************************************************************************

!!****f* qs_blacs/cp_fmstruct_retain [1.0] *
!!
!!   NAME
!!     cp_fmstruct_retain
!!
!!   SYNOPSIS
!!     Subroutine cp_fmstruct_retain(fmstruct, error)
!!       Type(cp_f_matrix_struct_type), Pointer:: fmstruct
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_fmstruct_retain
!!
!!   FUNCTION
!!     retains a full matrix structure
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - fmstruct: the structure to retain
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
SUBROUTINE cp_fmstruct_retain(fmstruct,error)
  TYPE(cp_f_matrix_struct_type), POINTER :: fmstruct
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='cp_fmstruct_retain',&
        routineP=moduleN//':'//routineN
  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(fmstruct),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(fmstruct%ref_count>0,cp_failure_level,routineP,error,failure)
     fmstruct%ref_count=fmstruct%ref_count+1
  END IF
END SUBROUTINE cp_fmstruct_retain
!***************************************************************************

!!****f* qs_blacs/cp_fmstruct_release [1.0] *
!!
!!   NAME
!!     cp_fmstruct_release
!!
!!   SYNOPSIS
!!     Subroutine cp_fmstruct_release(fmstruct, error)
!!       Type(cp_f_matrix_struct_type), Pointer:: fmstruct
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_fmstruct_release
!!
!!   FUNCTION
!!     releases a full matrix structure
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - fmstruct: the structure to release
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
SUBROUTINE cp_fmstruct_release(fmstruct,error)
  TYPE(cp_f_matrix_struct_type), POINTER :: fmstruct
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='cp_fmstruct_release',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat
  failure=.FALSE.
  
  IF (ASSOCIATED(fmstruct)) THEN
     CPPrecondition(fmstruct%ref_count>0,cp_failure_level,routineP,error,failure)
     fmstruct%ref_count=fmstruct%ref_count-1
     IF (fmstruct%ref_count<1) THEN
        CALL cp_blacs_env_release(fmstruct%context,error=error)
        CALL cp_para_env_release(fmstruct%para_env,error=error)
        IF (ASSOCIATED(fmstruct%row_indices)) THEN
           DEALLOCATE(fmstruct%row_indices,stat=stat)
           CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
        ENDIF
        IF (ASSOCIATED(fmstruct%col_indices)) THEN
           DEALLOCATE(fmstruct%col_indices,stat=stat)
           CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
        END IF
        DEALLOCATE(fmstruct, stat=stat)
        CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
     END IF
  END IF
  NULLIFY(fmstruct)
END SUBROUTINE cp_fmstruct_release
!***************************************************************************

!!****f* qs_blacs/cp_fmstruct_equivalent [1.0] *
!!
!!   NAME
!!     cp_fmstruct_equivalent
!!
!!   SYNOPSIS
!!     Function cp_fmstruct_equivalent(fmstruct1, fmstruct2, error)&
!!         Result(res)
!!       Type(cp_f_matrix_struct_type), Pointer:: fmstruct1, fmstruct2
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!       Logical:: res
!!     End Function cp_fmstruct_equivalent
!!
!!   FUNCTION
!!     returns true if the two matrix structures are equivalent, false
!!     otherwise.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - fmstruct1: one of the full matrix structures to compare
!!     - fmstruct2: the second of the full matrix structures to compare
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
FUNCTION cp_fmstruct_equivalent(fmstruct1,fmstruct2,error) RESULT(res)
  TYPE(cp_f_matrix_struct_type), POINTER :: fmstruct1, fmstruct2
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  LOGICAL :: res
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='cp_fmstruct_equivalent',&
        routineP=moduleN//':'//routineN
  INTEGER :: i
  failure=.FALSE.

  CPPrecondition(ASSOCIATED(fmstruct1),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(fmstruct2),cp_failure_level,routineP,error,failure)
  CPPrecondition(fmstruct1%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(fmstruct2%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     IF (fmstruct1%id_nr==fmstruct2%id_nr) THEN
        res=.TRUE.
     ELSE
        res= (fmstruct1%context%group==fmstruct2%context%group).AND.&
             (fmstruct1%nrow_global==fmstruct2%nrow_global).AND.&
             (fmstruct1%ncol_global==fmstruct2%ncol_global)
        DO i=1,9
           res=res.AND.(fmstruct1%descriptor(i)==fmstruct1%descriptor(i))
        END DO
     END IF
  END IF
END FUNCTION cp_fmstruct_equivalent
!***************************************************************************

!!****f* cp_f_matrix_struct/cp_fmstruct_write [1.0] *
!!
!!   NAME
!!     cp_fmstruct_write
!!
!!   FUNCTION
!!     writes a description of the matrix structure to the given unit
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - fmstruct: the structure to print
!!     - unit_nr: the number of the unit where to write the description
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
SUBROUTINE cp_fmstruct_write(fmstruct, unit_nr, long_description,error)
  TYPE(cp_f_matrix_struct_type), POINTER :: fmstruct
  INTEGER, INTENT(in) :: unit_nr
  LOGICAL, INTENT(in), OPTIONAL :: long_description
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure,l_desc
  INTEGER :: iostat
  CHARACTER (len=80) :: fmtstring
  CHARACTER(len=*), PARAMETER :: routineN='cp_fmstruct_write',&
        routineP=moduleN//':'//routineN
  failure=.FALSE.; l_desc=.FALSE.
  
  IF (PRESENT(long_description)) l_desc=long_description
  IF (ASSOCIATED(fmstruct)) THEN
     WRITE (unit=unit_nr,&
          fmt="('<cp_f_matrix_struct>:{ id_nr=',i10,', ref_count=',i10,',')",&
          iostat=iostat) fmstruct%id_nr, fmstruct%ref_count
     CPPostcondition(iostat==0,cp_warning_level,routineP,error,failure)
     WRITE (unit=unit_nr,fmt="(' nrow_global=',i8,', ncol_global=',i8,',')",&
          iostat=iostat) fmstruct%nrow_global, fmstruct%ncol_global
     WRITE (unit=unit_nr,fmt="(' nrow_block=',i8,', ncol_block=',i8,',')",&
          iostat=iostat) fmstruct%nrow_global, fmstruct%ncol_global
     IF (ASSOCIATED(fmstruct%nrow_local)) THEN
        WRITE (fmtstring,'(A)') "(' nrow_local=(',"//cp_to_string(SIZE(fmstruct%nrow_local))//&
             "(i8),'),')"
        WRITE (unit=unit_nr,fmt=fmtstring,iostat=iostat) fmstruct%nrow_local
     ELSE
        WRITE (unit=unit_nr,fmt="(a)") " nrow_local=*null*,"
     END IF
     IF (ASSOCIATED(fmstruct%ncol_local)) THEN
        WRITE (fmtstring,'(A)') "(' ncol_local=(',"//cp_to_string(SIZE(fmstruct%ncol_local))//&
             "(i8),'),')"
        WRITE (unit=unit_nr,fmt=fmtstring,iostat=iostat) fmstruct%ncol_local
     ELSE
        WRITE (unit=unit_nr,fmt="(a)") " ncol_local=*null*,"
     END IF
     WRITE (unit=unit_nr,fmt="(' first_p_pos=(',i8,', ',i8,'),')",&
          iostat=iostat) fmstruct%first_p_pos(1), fmstruct%first_p_pos(2)
     IF (ASSOCIATED(fmstruct%row_indices)) THEN
        IF (l_desc) THEN
           WRITE (unit=unit_nr,fmt="(a)",iostat=iostat) " row_indices=(,"
           WRITE (unit=unit_nr,fmt="(i1)",iostat=iostat)
           WRITE (unit=unit_nr,fmt="(a)",iostat=iostat) "),"
        ELSE
           WRITE (unit=unit_nr,&
                fmt="(' row_indices=<integer,dimension(',i8,'),')",&
                iostat=iostat) SIZE(fmstruct%row_indices)
        END IF
     ELSE
        WRITE (unit=unit_nr,fmt="(a)") " row_indices=*null*,"
     END IF
     IF (ASSOCIATED(fmstruct%col_indices)) THEN
        IF (l_desc) THEN
           WRITE (unit=unit_nr,fmt="(a)",iostat=iostat) " col_indices=("
           WRITE (unit=unit_nr,fmt="(i1)",iostat=iostat)
           WRITE (unit=unit_nr,fmt="(a)",iostat=iostat) "),"
        ELSE
           WRITE (unit=unit_nr,&
                fmt="(' col_indices=<integer,dimension(',i8,'),')",&
                iostat=iostat) SIZE(fmstruct%col_indices)
        END IF
     ELSE
        WRITE (unit=unit_nr,fmt="(a)") " col_indices=*null*,"
     END IF
     IF (ASSOCIATED(fmstruct%context)) THEN
        IF (l_desc) THEN
           WRITE (unit=unit_nr,fmt="(a)",iostat=iostat) " context="
           CALL cp_blacs_env_write(fmstruct%context,unit_nr,error=error)
        ELSE
           WRITE (unit=unit_nr,fmt="(a,i10)",iostat=iostat) " context%group=",&
                fmstruct%context%group
        END IF
     ELSE
        WRITE (unit=unit_nr,fmt="(a)",iostat=iostat) " context=*null*"
     END IF
     IF (ASSOCIATED(fmstruct%para_env)) THEN
        IF (l_desc) THEN
           WRITE (unit=unit_nr,fmt="(a)",iostat=iostat) " para_env="
           CALL cp_para_env_write(fmstruct%para_env,unit_nr,error=error)
        ELSE
           WRITE (unit=unit_nr,fmt="(a,i10)",iostat=iostat) &
                " para_env%group=",fmstruct%para_env%group
        END IF
     ELSE
        WRITE (unit=unit_nr,fmt="(a)",iostat=iostat) " para_env=*null*"
     END IF
     WRITE (unit=unit_nr, fmt="(' descriptor=(',9i8,') }')",iostat=iostat)&
          fmstruct%descriptor
  ELSE
     WRITE (unit=unit_nr, fmt="(a)",iostat=iostat)"<cp_f_matrix_struct>:*null*"
  END IF
END SUBROUTINE cp_fmstruct_write
!***************************************************************************

END MODULE cp_f_matrix_struct
