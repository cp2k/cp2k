!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_fm_struct [1.0] *
!!
!!   NAME
!!     cp_fm_struct
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
MODULE cp_fm_struct
  USE cp_blacs_env,                    ONLY: cp_blacs_env_release,&
                                             cp_blacs_env_retain,&
                                             cp_blacs_env_write
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
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain,&
                                             cp_para_env_write
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE kinds,                           ONLY: wp=>dp
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_max,&
                                             mp_sum,&
                                             mp_sync
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_fm_struct'
  INTEGER,SAVE,PRIVATE :: last_fmstruct_id_nr=0

! the default blacs block sizes
! consider using #ifdefs to give them the optimal values
  INTEGER, PARAMETER :: optimal_blacs_col_block_size=32
  INTEGER, PARAMETER :: optimal_blacs_row_block_size=32

  PUBLIC :: cp_fm_struct_type
  PUBLIC :: optimal_blacs_row_block_size, optimal_blacs_col_block_size
  PUBLIC :: cp_fm_struct_create,cp_fm_struct_retain, cp_fm_struct_release,&
       cp_fm_struct_equivalent, cp_fm_struct_write, cp_fm_struct_get
!***
!****************************************************************************

!!****s* cp_fm_struct/cp_fm_struct_type [1.0] *
!!
!!   NAME
!!     cp_fm_struct_type
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
  TYPE cp_fm_struct_type
     TYPE(cp_para_env_type), POINTER :: para_env
     TYPE(cp_blacs_env_type), POINTER :: context
     INTEGER, DIMENSION(9) :: descriptor
     ! when using blacs
     ! ncol_block=descriptor(6)
     ! ncol_global=descriptor(4)
     ! nrow_block=descriptor(5)
     ! nrow_global=descriptor(3)
     INTEGER :: nrow_block, ncol_block, nrow_global, ncol_global
     INTEGER, DIMENSION(2) :: first_p_pos
     INTEGER, DIMENSION(:), POINTER :: row_indices, col_indices,&
          nrow_local, ncol_local
     INTEGER :: id_nr, ref_count
  END TYPE cp_fm_struct_type
!!***
!****************************************************************************

CONTAINS

!!****f* cp_fm_struct/cp_fm_struct_create [1.0] *
!!
!!   NAME
!!     cp_fm_struct_create
!!
!!   SYNOPSIS
!!     Subroutine cp_fm_struct_create(fmstruct, para_env, context,&
!!         nrow_global, ncol_global, nrow_block, ncol_block, descriptor,&
!!         first_p_pos, template_fmstruct, error)
!!       Type(cp_fm_struct_type), Pointer:: fmstruct
!!       Type(cp_para_env_type), Pointer, Optional:: para_env
!!       Integer, Intent (IN), Optional:: nrow_global, ncol_global
!!       Integer, Intent (IN), Optional:: nrow_block, ncol_block
!!       Type(cp_blacs_env_type), Pointer, Optional:: context
!!       Integer, Dimension(9), Intent (IN), Optional:: descriptor
!!       Integer, Optional, Dimension(2):: first_p_pos
!!       Type(cp_fm_struct_type), Pointer, Optional:: template_fmstruct
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_fm_struct_create
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
!!     - template_fmstruct: a matrix structure where to take the default values
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
SUBROUTINE cp_fm_struct_create(fmstruct,para_env,context, nrow_global,&
     ncol_global, nrow_block, ncol_block, descriptor, first_p_pos, &
     template_fmstruct, error)
  TYPE(cp_fm_struct_type), POINTER       :: fmstruct
  TYPE(cp_para_env_type), POINTER, OPTIONAL    :: para_env
  INTEGER, INTENT(in), OPTIONAL                :: nrow_global,ncol_global
  INTEGER, INTENT(in), OPTIONAL                :: nrow_block,ncol_block
  TYPE(cp_blacs_env_type), POINTER, OPTIONAL   :: context
  INTEGER, DIMENSION(9), INTENT(in), OPTIONAL  :: descriptor
  INTEGER,OPTIONAL, DIMENSION(2)               :: first_p_pos
  TYPE(cp_fm_struct_type), POINTER, OPTIONAL :: template_fmstruct
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
#if defined(__SCALAPACK)
    INTEGER, EXTERNAL :: numroc
#endif

  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='cp_fm_struct_create',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat, i
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
     fmstruct%nrow_block=optimal_blacs_row_block_size
     fmstruct%ncol_block=optimal_blacs_col_block_size

     IF (.NOT.PRESENT(template_fmstruct)) THEN
        CPPrecondition(PRESENT(para_env),cp_failure_level,routineP,error,failure)
        CPPrecondition(PRESENT(context),cp_failure_level,routineP,error,failure)
        CPPrecondition(PRESENT(nrow_global),cp_failure_level,routineP,error,failure)
        CPPrecondition(PRESENT(ncol_global),cp_failure_level,routineP,error,failure)
     ELSE
        fmstruct%context => template_fmstruct%context
        fmstruct%para_env => template_fmstruct%para_env
        fmstruct%descriptor = template_fmstruct %descriptor
        fmstruct %nrow_block = template_fmstruct %nrow_block
        fmstruct %nrow_global = template_fmstruct %nrow_global
        fmstruct %ncol_block = template_fmstruct %ncol_block
        fmstruct %first_p_pos = template_fmstruct %first_p_pos
     END IF

  END IF
  IF (.NOT.failure) THEN
     
     IF (PRESENT(context)) fmstruct%context => context
     CALL cp_blacs_env_retain(fmstruct%context,error=error)
     IF (PRESENT(para_env)) fmstruct%para_env => para_env
     CALL cp_para_env_retain(fmstruct%para_env,error=error)
     
     IF (PRESENT(nrow_global)) fmstruct%nrow_global=nrow_global
     IF (PRESENT(ncol_global)) fmstruct%ncol_global=ncol_global
     ! checks for more processes than matrix elements
     CALL cp_assert(fmstruct%context%num_pe(1)<=fmstruct%nrow_global,&
          cp_warning_level,cp_assertion_failed,routineP,&
          "context%num_pe(1)<=nrow_global"//&
          CPSourceFileRef,&
          error)
     CALL cp_assert(fmstruct%context%num_pe(2)<=fmstruct%ncol_global,&
          cp_warning_level,cp_assertion_failed,routineP,&
          "context%num_pe(2)<=ncol_global"//&
          CPSourceFileRef,&
          error)
     IF (PRESENT(nrow_block)) fmstruct%nrow_block=nrow_block
     IF (PRESENT(ncol_block)) fmstruct%ncol_block=ncol_block
     fmstruct%nrow_block=MAX(1,&
          MIN(fmstruct%nrow_block,&
          fmstruct%nrow_global/fmstruct%context%num_pe(1)))
     fmstruct%ncol_block=MAX(1,&
          MIN(fmstruct%ncol_block,&
          fmstruct%ncol_global/fmstruct%context%num_pe(2)))
     IF (fmstruct%nrow_global==fmstruct%ncol_global) THEN
        ! square matrix -> square blocks (otherwise some op fail)
        fmstruct%nrow_block=MIN(fmstruct%nrow_block,fmstruct%ncol_block)
        fmstruct%ncol_block=fmstruct%nrow_block
     END IF
     ALLOCATE(fmstruct%nrow_local(0:(fmstruct%context%num_pe(1)-1)),&
          fmstruct%ncol_local(0:(fmstruct%context%num_pe(2)-1)),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF (.NOT.failure) THEN
     IF (.NOT.PRESENT(template_fmstruct)) &
          fmstruct%first_p_pos=fmstruct%context%source
     IF (PRESENT(first_p_pos)) fmstruct%first_p_pos=first_p_pos

     fmstruct%nrow_local=0
     fmstruct%ncol_local=0
#ifdef __SCALAPACK
     fmstruct%nrow_local(fmstruct%context%mepos(1))= &
          numroc(fmstruct%nrow_global,fmstruct%nrow_block,&
          fmstruct%context%mepos(1),fmstruct%first_p_pos(1),&
          fmstruct%context%num_pe(1))
     fmstruct%ncol_local(fmstruct%context%mepos(2))= &
          numroc(fmstruct%ncol_global,fmstruct%ncol_block,&
          fmstruct%context%mepos(2),fmstruct%first_p_pos(2),&
          fmstruct%context%num_pe(2))
     CALL mp_sum(fmstruct%nrow_local,fmstruct%para_env%group)
     CALL mp_sum(fmstruct%ncol_local,fmstruct%para_env%group)
     fmstruct%nrow_local(:)=fmstruct%nrow_local(:)/fmstruct%context%num_pe(2)
     fmstruct%ncol_local(:)=fmstruct%ncol_local(:)/fmstruct%context%num_pe(1)
     do i=0,(size(fmstruct%nrow_local)-1)
        fmstruct%nrow_local(i)=max(1,fmstruct%nrow_local(i))
     end do
     do i=0,(size(fmstruct%ncol_local)-1)
        fmstruct%ncol_local(i)=max(1,fmstruct%ncol_local(i))
     end do
#else
     ! block = full matrix
     fmstruct%nrow_block=fmstruct%nrow_global
     fmstruct%ncol_block=fmstruct%ncol_global
     fmstruct%nrow_local(fmstruct%context%mepos(1))=fmstruct%nrow_global
     fmstruct%ncol_local(fmstruct%context%mepos(2))=fmstruct%ncol_global
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
END SUBROUTINE cp_fm_struct_create
!***************************************************************************

!!****f* cp_fm_struct/cp_fm_struct_retain [1.0] *
!!
!!   NAME
!!     cp_fm_struct_retain
!!
!!   SYNOPSIS
!!     Subroutine cp_fm_struct_retain(fmstruct, error)
!!       Type(cp_fm_struct_type), Pointer:: fmstruct
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_fm_struct_retain
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
SUBROUTINE cp_fm_struct_retain(fmstruct,error)
  TYPE(cp_fm_struct_type), POINTER :: fmstruct
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='cp_fm_struct_retain',&
        routineP=moduleN//':'//routineN
  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(fmstruct),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(fmstruct%ref_count>0,cp_failure_level,routineP,error,failure)
     fmstruct%ref_count=fmstruct%ref_count+1
  END IF
END SUBROUTINE cp_fm_struct_retain
!***************************************************************************

!!****f* cp_fm_struct/cp_fm_struct_release [1.0] *
!!
!!   NAME
!!     cp_fm_struct_release
!!
!!   SYNOPSIS
!!     Subroutine cp_fm_struct_release(fmstruct, error)
!!       Type(cp_fm_struct_type), Pointer:: fmstruct
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_fm_struct_release
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
SUBROUTINE cp_fm_struct_release(fmstruct,error)
  TYPE(cp_fm_struct_type), POINTER :: fmstruct
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='cp_fm_struct_release',&
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
END SUBROUTINE cp_fm_struct_release
!***************************************************************************

!!****f* cp_fm_struct/cp_fm_struct_equivalent [1.0] *
!!
!!   NAME
!!     cp_fm_struct_equivalent
!!
!!   SYNOPSIS
!!     Function cp_fm_struct_equivalent(fmstruct1, fmstruct2, error)&
!!         Result(res)
!!       Type(cp_fm_struct_type), Pointer:: fmstruct1, fmstruct2
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!       Logical:: res
!!     End Function cp_fm_struct_equivalent
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
FUNCTION cp_fm_struct_equivalent(fmstruct1,fmstruct2,error) RESULT(res)
  TYPE(cp_fm_struct_type), POINTER :: fmstruct1, fmstruct2
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  LOGICAL :: res
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='cp_fm_struct_equivalent',&
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
END FUNCTION cp_fm_struct_equivalent
!***************************************************************************

!!****f* cp_fm_struct/cp_fm_struct_write [1.0] *
!!
!!   NAME
!!     cp_fm_struct_write
!!
!!   SYNOPSIS
!!     Subroutine cp_fm_struct_write(fmstruct, unit_nr, long_description,&
!!         error)
!!       Type(cp_fm_struct_type), Pointer:: fmstruct
!!       Integer, Intent (IN):: unit_nr
!!       Logical, Intent (IN), Optional:: long_description
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_fm_struct_write
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
SUBROUTINE cp_fm_struct_write(fmstruct, unit_nr, long_description,error)
  TYPE(cp_fm_struct_type), POINTER :: fmstruct
  INTEGER, INTENT(in) :: unit_nr
  LOGICAL, INTENT(in), OPTIONAL :: long_description
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure,l_desc
  INTEGER :: iostat
  CHARACTER (len=80) :: fmtstring
  CHARACTER(len=*), PARAMETER :: routineN='cp_fm_struct_write',&
        routineP=moduleN//':'//routineN
  failure=.FALSE.; l_desc=.FALSE.
  
  IF (PRESENT(long_description)) l_desc=long_description
  IF (ASSOCIATED(fmstruct)) THEN
     WRITE (unit=unit_nr,&
          fmt="('<cp_fm_struct>:{ id_nr=',i10,', ref_count=',i10,',')",&
          iostat=iostat) fmstruct%id_nr, fmstruct%ref_count
     CPPostcondition(iostat==0,cp_warning_level,routineP,error,failure)
     WRITE (unit=unit_nr,fmt="(' nrow_global=',i8,', ncol_global=',i8,',')",&
          iostat=iostat) fmstruct%nrow_global, fmstruct%ncol_global
     WRITE (unit=unit_nr,fmt="(' nrow_block=',i8,', ncol_block=',i8,',')",&
          iostat=iostat) fmstruct%nrow_block, fmstruct%ncol_block
     IF (ASSOCIATED(fmstruct%nrow_local)) THEN
        WRITE (fmtstring,'(A,i6,A)') "(' nrow_local=(',",&
             SIZE(fmstruct%nrow_local),"(i8),'),')"
        WRITE (unit=unit_nr,fmt=fmtstring,iostat=iostat) fmstruct%nrow_local
     ELSE
        WRITE (unit=unit_nr,fmt="(a)") " nrow_local=*null*,"
     END IF
     IF (ASSOCIATED(fmstruct%ncol_local)) THEN
        WRITE (fmtstring,'(A,i6,A)') "(' ncol_local=(',",SIZE(fmstruct%ncol_local),&
             "(i8),'),')"
        WRITE (unit=unit_nr,fmt=fmtstring,iostat=iostat) fmstruct%ncol_local
     ELSE
        WRITE (unit=unit_nr,fmt="(a)") " ncol_local=*null*,"
     END IF
     WRITE (unit=unit_nr,fmt="(' first_p_pos=(',i8,', ',i8,'),')",&
          iostat=iostat) fmstruct%first_p_pos(1), fmstruct%first_p_pos(2)
     IF (ASSOCIATED(fmstruct%row_indices)) THEN
        IF (l_desc) THEN
           WRITE (unit=unit_nr,fmt="(a)",iostat=iostat,advance='no')&
                " row_indices=(,"
           WRITE (unit=unit_nr,fmt="(i8)",iostat=iostat, advance='no')&
                fmstruct%row_indices
           WRITE (unit=unit_nr,fmt="(a)",iostat=iostat) "),"
        ELSE
           WRITE (unit=unit_nr,&
                fmt="(' row_indices=<INTEGER,DIMENSION(',i8,'),')",&
                iostat=iostat) SIZE(fmstruct%row_indices)
        END IF
     ELSE
        WRITE (unit=unit_nr,fmt="(a)") " row_indices=*null*,"
     END IF
     IF (ASSOCIATED(fmstruct%col_indices)) THEN
        IF (l_desc) THEN
           WRITE (unit=unit_nr,fmt="(a)",iostat=iostat, advance='no')&
                " col_indices=("
           WRITE (unit=unit_nr,fmt="(i8)",iostat=iostat, advance='no')&
                fmstruct%col_indices
           WRITE (unit=unit_nr,fmt="(a)",iostat=iostat) "),"
        ELSE
           WRITE (unit=unit_nr,&
                fmt="(' col_indices=<INTEGER,DIMENSION(',i8,'),')",&
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
     WRITE (unit=unit_nr, fmt="(a)",iostat=iostat)"<cp_fm_struct>:*null*"
  END IF
END SUBROUTINE cp_fm_struct_write
!***************************************************************************

!!****f* cp_fm_struct/cp_fm_struct_get [1.0] *
!!
!!   NAME
!!     cp_fm_struct_get
!!
!!   SYNOPSIS
!!     Subroutine cp_fm_struct_get(fmstruct, para_env, context, descriptor,&
!!         ncol_block, nrow_block, nrow_global, ncol_global, id_nr,&
!!         ref_count, first_p_pos, row_indices, col_indices, nrow_local,&
!!         ncol_local, error)
!!       Type(cp_fm_struct_type), Pointer:: fmstruct
!!       Type(cp_para_env_type), Pointer, Optional:: para_env
!!       Type(cp_blacs_env_type), Pointer, Optional:: context
!!       Integer, Dimension(9), Intent (OUT), Optional:: descriptor
!!       Integer, Intent (OUT), Optional:: ncol_block, nrow_block,&
!!         nrow_global, ncol_global, id_nr, ref_count
!!       Integer, Dimension(2), Intent (OUT), Optional:: first_p_pos
!!       Integer, Dimension(:), Pointer, Optional:: row_indices,&
!!         col_indices, nrow_local, ncol_local
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine cp_fm_struct_get
!!
!!   FUNCTION
!!     returns the values of various attributes of the matrix structure
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - fmstruct: the stru
!!     - unit_nr: the number of the unit where to write the descr
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
SUBROUTINE cp_fm_struct_get(fmstruct, para_env, context,&
     descriptor, ncol_block, nrow_block, nrow_global,&
     ncol_global, id_nr, ref_count, first_p_pos, row_indices,&
     col_indices, nrow_local, ncol_local, error)
  TYPE(cp_fm_struct_type), POINTER :: fmstruct
  TYPE(cp_para_env_type), POINTER, OPTIONAL :: para_env
  TYPE(cp_blacs_env_type), POINTER, OPTIONAL :: context
  INTEGER, DIMENSION(9), INTENT(OUT), OPTIONAL :: descriptor
  INTEGER, INTENT(out), OPTIONAL :: ncol_block, nrow_block, nrow_global,&
       ncol_global, id_nr, ref_count
  INTEGER, DIMENSION(2), INTENT(out), OPTIONAL :: first_p_pos
  INTEGER, DIMENSION(:), POINTER, OPTIONAL :: row_indices, col_indices,&
       nrow_local, ncol_local
  TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='cp_fm_struct_get',&
       routineP=moduleN//':'//routineN
  failure=.FALSE.

  CPPrecondition(ASSOCIATED(fmstruct),cp_failure_level,routineP,error,failure)
  CPPrecondition(fmstruct%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     IF (PRESENT(para_env)) para_env => fmstruct%para_env
     IF (PRESENT(context)) context => fmstruct%context
     IF (PRESENT(descriptor)) descriptor = fmstruct%descriptor
     IF (PRESENT(ncol_block)) ncol_block=fmstruct%ncol_block
     IF (PRESENT(nrow_block)) nrow_block=fmstruct%nrow_block
     IF (PRESENT(nrow_global)) nrow_global=fmstruct%nrow_global
     IF (PRESENT(ncol_global)) ncol_global=fmstruct%ncol_global
     IF (PRESENT(id_nr)) id_nr=fmstruct%id_nr
     IF (PRESENT(ref_count)) ref_count=fmstruct%ref_count
     IF (PRESENT(first_p_pos)) first_p_pos=fmstruct%first_p_pos
     IF (PRESENT(row_indices)) row_indices => fmstruct%row_indices
     IF (PRESENT(col_indices)) col_indices => fmstruct%col_indices
     IF (PRESENT(nrow_local)) nrow_local => nrow_local
     IF (PRESENT(ncol_local)) ncol_local => ncol_local
  END IF
END SUBROUTINE cp_fm_struct_get
!***************************************************************************

END MODULE cp_fm_struct
