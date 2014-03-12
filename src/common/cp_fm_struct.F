!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief represent the structure of a full matrix
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE cp_fm_struct
  USE cp_blacs_env,                    ONLY: cp_blacs_env_release,&
                                             cp_blacs_env_retain,&
                                             cp_blacs_env_write
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain,&
                                             cp_para_env_write
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE f77_blas
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE machine,                         ONLY: m_flush
  USE message_passing,                 ONLY: mp_sum
#include "cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_fm_struct'
  INTEGER,SAVE,PRIVATE :: last_fmstruct_id_nr=0

! the default blacs block sizes
! consider using #ifdefs to give them the optimal values
! these can be changed using scf_control
! *** these are used by default
  INTEGER :: optimal_blacs_col_block_size=32
  INTEGER :: optimal_blacs_row_block_size=32

  PUBLIC :: cp_fm_struct_type,cp_fm_struct_p_type
  PUBLIC :: optimal_blacs_row_block_size, optimal_blacs_col_block_size
  PUBLIC :: cp_fm_struct_create,cp_fm_struct_retain, cp_fm_struct_release,&
       cp_fm_struct_equivalent, cp_fm_struct_compatible, cp_fm_struct_write,&
       cp_fm_struct_get,cp_fm_struct_double
!***

! *****************************************************************************
!> \brief keeps the information about the structure of a full matrix
!> \param para_env the parallel environment of the matrixes with this structure
!> \param context the blacs context (parallel environment for scalapack),
!>        should be compatible wih para_env
!> \param descriptor the scalapack descriptor of the matrixes, when using
!>        scalapack (ncol_block=descriptor(6), ncol_global=descriptor(4),
!>        nrow_block=descriptor(5), nrow_global=descriptor(3))
!> \param ncol_block number of columns of a scalapack block
!> \param nrow_block number of rows of a scalapack block
!> \param nrow_global number of rows of the matrix
!> \param ncol_global number of rows
!> \param first_p_pos position of the first processor (for scalapack)
!> \param row_indices real (global) indexes of the rows (defined only for
!>        the local rows really used)
!> \param col_indices real (global) indexes of the cols (defined only for
!>        the local cols really used)
!> \param nrow_locals nrow_locals(i) number of local rows of the matrix really
!>        used on the processors with context%mepos(1)==i
!> \param ncol_locals ncol_locals(i) number of local rows of the matrix really
!>        used on the processors with context%mepos(2)==i
!> \param id_nr an identification number (guaranteed to be unique)
!> \param ref_count reference count (see doc/ReferenceCounting.html)
!> \param local_leading_dimension leading dimension of the data that is
!>        stored on this processor
!>
!>      readonly attributes:
!> \param nrow_local number of local rows really used on the actual processor
!> \param ncol_local number of local cols really used on the actual processor
!> \note
!>      use cp_fm_struct_get to extract information from this structure
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE cp_fm_struct_type
     TYPE(cp_para_env_type), POINTER :: para_env
     TYPE(cp_blacs_env_type), POINTER :: context
     INTEGER, DIMENSION(9) :: descriptor
     INTEGER :: nrow_block, ncol_block, nrow_global, ncol_global
     INTEGER, DIMENSION(2) :: first_p_pos
     INTEGER, DIMENSION(:), POINTER :: row_indices, col_indices,&
          nrow_locals, ncol_locals
     INTEGER :: id_nr, ref_count, local_leading_dimension
  END TYPE cp_fm_struct_type
! *****************************************************************************
  TYPE cp_fm_struct_p_type
     TYPE(cp_fm_struct_type), POINTER :: struct
  END TYPE cp_fm_struct_p_type

CONTAINS

! *****************************************************************************
!> \brief allocates and initializes a full matrix structure
!> \param fmstruct the pointer that will point to the new structure
!> \param para_env the parallel environment
!> \param nrow_global the number of row of the full matrix
!> \param ncol_global the number of colums of the full matrix
!> \param nrow_block the number of rows of a block of the matrix ! useful defaults are in, do not specify
!> \param ncol_block the number of colums of a block of the matrix ! useful defaults are in, do not specify
!> \param context the blacs context of this matrix
!> \param local_leading_dimension the leading dimension of the locally stored
!>        data block
!> \param descriptor the scalapack descriptor of the matrix (if not given
!>        a new one is allocated
!> \param template_fmstruct a matrix structure where to take the default values
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE cp_fm_struct_create(fmstruct,para_env,context, nrow_global,&
     ncol_global, nrow_block, ncol_block, descriptor, first_p_pos, &
     local_leading_dimension, template_fmstruct, square_blocks, force_block, error)
  TYPE(cp_fm_struct_type), POINTER       :: fmstruct
  TYPE(cp_para_env_type), POINTER, OPTIONAL    :: para_env
  INTEGER, INTENT(in), OPTIONAL                :: nrow_global,ncol_global
  INTEGER, INTENT(in), OPTIONAL                :: nrow_block,ncol_block
  INTEGER, INTENT(in), OPTIONAL                :: local_leading_dimension
  TYPE(cp_blacs_env_type), POINTER, OPTIONAL   :: context
  INTEGER, DIMENSION(9), INTENT(in), OPTIONAL  :: descriptor
  INTEGER,OPTIONAL, DIMENSION(2)               :: first_p_pos
  TYPE(cp_fm_struct_type), POINTER, OPTIONAL   :: template_fmstruct
  LOGICAL, OPTIONAL, INTENT(in)                :: square_blocks
  LOGICAL, OPTIONAL, INTENT(in)                :: force_block
  TYPE(cp_error_type), INTENT(inout) :: error

  INTEGER ::  dumblock

  INTEGER :: iunit
  TYPE(cp_logger_type), POINTER :: logger
#if defined(__SCALAPACK)
    INTEGER, EXTERNAL :: numroc
#endif

  LOGICAL :: failure, my_square_blocks, my_force_block
  CHARACTER(len=*), PARAMETER :: routineN='cp_fm_struct_create',&
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
     fmstruct%nrow_block=optimal_blacs_row_block_size
     fmstruct%ncol_block=optimal_blacs_col_block_size

     IF (.NOT.PRESENT(template_fmstruct)) THEN
        CPPrecondition(PRESENT(context),cp_failure_level,routineP,error,failure)
        CPPrecondition(PRESENT(nrow_global),cp_failure_level,routineP,error,failure)
        CPPrecondition(PRESENT(ncol_global),cp_failure_level,routineP,error,failure)
        fmstruct %local_leading_dimension = 1
     ELSE
        fmstruct%context => template_fmstruct%context
        fmstruct%para_env => template_fmstruct%para_env
        fmstruct%descriptor = template_fmstruct %descriptor
        fmstruct %nrow_block = template_fmstruct %nrow_block
        fmstruct %nrow_global = template_fmstruct %nrow_global
        fmstruct %ncol_block = template_fmstruct %ncol_block
        fmstruct %ncol_global = template_fmstruct %ncol_global
        fmstruct %first_p_pos = template_fmstruct %first_p_pos
        fmstruct %local_leading_dimension = &
             template_fmstruct %local_leading_dimension
     END IF

  END IF
  my_force_block=.FALSE.
  IF (PRESENT(force_block)) my_force_block=force_block
  IF (.NOT.failure) THEN

     IF (PRESENT(context)) THEN
        fmstruct%context => context
        fmstruct%para_env => context%para_env
     END IF
     IF (PRESENT(para_env)) fmstruct%para_env => para_env
     CALL cp_blacs_env_retain(fmstruct%context,error=error)
     CALL cp_para_env_retain(fmstruct%para_env,error=error)

     IF (PRESENT(nrow_global)) fmstruct%nrow_global=nrow_global
     IF (PRESENT(ncol_global)) fmstruct%ncol_global=ncol_global

     ! try to avoid small left-over blocks (anyway naive)
     IF (PRESENT(nrow_block)) THEN
         fmstruct%nrow_block=nrow_block
     END IF
     IF(.NOT. my_force_block) THEN
       dumblock = CEILING(REAL(fmstruct%nrow_global,KIND=dp)/ &
                REAL(fmstruct%context%num_pe(1),KIND=dp))
       fmstruct%nrow_block=MAX(1, MIN(fmstruct%nrow_block, dumblock))
     END IF
     IF (PRESENT(ncol_block)) THEN
       fmstruct%ncol_block=ncol_block
     END IF
     IF(.NOT. my_force_block) THEN
       dumblock = CEILING(REAL(fmstruct%ncol_global,KIND=dp)/ &
                  REAL(fmstruct%context%num_pe(2),KIND=dp))
       fmstruct%ncol_block=MAX(1, MIN(fmstruct%ncol_block, dumblock))
     END IF

     ! square matrix -> square blocks (otherwise some op fail)
     my_square_blocks = fmstruct%nrow_global==fmstruct%ncol_global
     IF (PRESENT(square_blocks)) my_square_blocks=square_blocks
     IF (my_square_blocks) THEN
        fmstruct%nrow_block=MIN(fmstruct%nrow_block,fmstruct%ncol_block)
        fmstruct%ncol_block=fmstruct%nrow_block
     END IF

     ALLOCATE(fmstruct%nrow_locals(0:(fmstruct%context%num_pe(1)-1)),&
          fmstruct%ncol_locals(0:(fmstruct%context%num_pe(2)-1)),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF (.NOT.failure) THEN
     IF (.NOT.PRESENT(template_fmstruct)) &
          fmstruct%first_p_pos=(/0,0/)
     IF (PRESENT(first_p_pos)) fmstruct%first_p_pos=first_p_pos

     fmstruct%nrow_locals=0
     fmstruct%ncol_locals=0
#ifdef __SCALAPACK
     fmstruct%nrow_locals(fmstruct%context%mepos(1))= &
          numroc(fmstruct%nrow_global,fmstruct%nrow_block,&
          fmstruct%context%mepos(1),fmstruct%first_p_pos(1),&
          fmstruct%context%num_pe(1))
     fmstruct%ncol_locals(fmstruct%context%mepos(2))= &
          numroc(fmstruct%ncol_global,fmstruct%ncol_block,&
          fmstruct%context%mepos(2),fmstruct%first_p_pos(2),&
          fmstruct%context%num_pe(2))
     CALL mp_sum(fmstruct%nrow_locals,fmstruct%para_env%group)
     CALL mp_sum(fmstruct%ncol_locals,fmstruct%para_env%group)
     fmstruct%nrow_locals(:)=fmstruct%nrow_locals(:)/fmstruct%context%num_pe(2)
     fmstruct%ncol_locals(:)=fmstruct%ncol_locals(:)/fmstruct%context%num_pe(1)

     IF (SUM(fmstruct%ncol_locals).NE.fmstruct%ncol_global .OR. SUM(fmstruct%nrow_locals).NE.fmstruct%nrow_global) THEN
        ! try to collect some output if this is going to happen again
        ! this seems to trigger on blanc, but should really never happen
        logger => cp_error_get_logger(error)
        iunit=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
        WRITE(iunit,*) "mepos",fmstruct%context%mepos(1:2),"numpe",fmstruct%context%num_pe(1:2)
        WRITE(iunit,*) "ncol_global",fmstruct%ncol_global
        WRITE(iunit,*) "nrow_global",fmstruct%nrow_global
        WRITE(iunit,*) "ncol_locals",fmstruct%ncol_locals
        WRITE(iunit,*) "nrow_locals",fmstruct%nrow_locals
        CALL m_flush(iunit)
     ENDIF

     CALL cp_assert(SUM(fmstruct%ncol_locals).EQ.fmstruct%ncol_global,&
          cp_failure_level,cp_assertion_failed,routineP,&
          "sum of local cols not equal global cols"//&
CPSourceFileRef,&
          error)
     CALL cp_assert(SUM(fmstruct%nrow_locals).EQ.fmstruct%nrow_global,&
          cp_failure_level,cp_assertion_failed,routineP,&
          "sum of local row not equal global rows"//&
CPSourceFileRef,&
          error)
#else
     ! block = full matrix
     fmstruct%nrow_block=fmstruct%nrow_global
     fmstruct%ncol_block=fmstruct%ncol_global
     fmstruct%nrow_locals(fmstruct%context%mepos(1))=fmstruct%nrow_global
     fmstruct%ncol_locals(fmstruct%context%mepos(2))=fmstruct%ncol_global
#endif

     fmstruct%local_leading_dimension = MAX(fmstruct%local_leading_dimension,&
          fmstruct%nrow_locals(fmstruct%context%mepos(1)))
     IF (PRESENT(local_leading_dimension)) THEN
        CALL cp_assert(MAX(1,fmstruct%nrow_locals(fmstruct%context%mepos(1)))&
             <= local_leading_dimension,cp_failure_level,cp_assertion_failed,&
             routineP,"local_leading_dimension too small ("//&
             cp_to_string(local_leading_dimension)//"<"//&
             cp_to_string(fmstruct%local_leading_dimension)//")"//&
CPSourceFileRef,&
             error=error,failure=failure)
        fmstruct%local_leading_dimension=local_leading_dimension
     END IF

     NULLIFY(fmstruct%row_indices, fmstruct%col_indices)
     last_fmstruct_id_nr=last_fmstruct_id_nr+1
     fmstruct%id_nr=last_fmstruct_id_nr
     fmstruct%ref_count=1

     IF (PRESENT(descriptor)) THEN
        fmstruct%descriptor=descriptor
     ELSE
        fmstruct%descriptor=0
#if defined(__SCALAPACK)
        ! local leading dimension needs to be at least 1
        CALL descinit(fmstruct%descriptor,fmstruct%nrow_global,&
             fmstruct%ncol_global,fmstruct%nrow_block,&
             fmstruct%ncol_block,fmstruct%first_p_pos(1),&
             fmstruct%first_p_pos(2),fmstruct%context%group,&
             fmstruct%local_leading_dimension,stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
#endif
     END IF
  END IF
END SUBROUTINE cp_fm_struct_create

! *****************************************************************************
!> \brief retains a full matrix structure
!> \param fmstruct the structure to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE cp_fm_struct_retain(fmstruct,error)
    TYPE(cp_fm_struct_type), POINTER         :: fmstruct
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_struct_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(fmstruct),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(fmstruct%ref_count>0,cp_failure_level,routineP,error,failure)
     fmstruct%ref_count=fmstruct%ref_count+1
  END IF
END SUBROUTINE cp_fm_struct_retain

! *****************************************************************************
!> \brief releases a full matrix structure
!> \param fmstruct the structure to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE cp_fm_struct_release(fmstruct,error)
    TYPE(cp_fm_struct_type), POINTER         :: fmstruct
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_struct_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

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
        IF (ASSOCIATED(fmstruct%nrow_locals)) THEN
           DEALLOCATE(fmstruct%nrow_locals,stat=stat)
           CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
        END IF
        IF (ASSOCIATED(fmstruct%ncol_locals)) THEN
           DEALLOCATE(fmstruct%ncol_locals,stat=stat)
           CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
        END IF
        DEALLOCATE(fmstruct, stat=stat)
        CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
     END IF
  END IF
  NULLIFY(fmstruct)
END SUBROUTINE cp_fm_struct_release

! *****************************************************************************
!> \brief returns true if the two matrix structures are equivalent, false
!>      otherwise.
!> \param fmstruct 1: one of the full matrix structures to compare
!> \param fmstruct 2: the second of the full matrix structures to compare
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
FUNCTION cp_fm_struct_equivalent(fmstruct1,fmstruct2,error) RESULT(res)
    TYPE(cp_fm_struct_type), POINTER         :: fmstruct1, fmstruct2
    TYPE(cp_error_type), INTENT(inout)       :: error
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_struct_equivalent', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure

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
             (fmstruct1%ncol_global==fmstruct2%ncol_global).AND.&
             (fmstruct1%local_leading_dimension==&
             fmstruct2%local_leading_dimension)
        DO i=1,9
           res=res.AND.(fmstruct1%descriptor(i)==fmstruct1%descriptor(i))
        END DO
     END IF
  END IF
END FUNCTION cp_fm_struct_equivalent

! *****************************************************************************
!> \brief returns true if the two matrix structures are compatible,
!>      (i.e. in the same parallel environment) false otherwise.
!> \param fmstruct 1: one of the full matrix structures to compare
!> \param fmstruct 2: the second of the full matrix structures to compare
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
FUNCTION cp_fm_struct_compatible(fmstruct1,fmstruct2,error) RESULT(res)
    TYPE(cp_fm_struct_type), POINTER         :: fmstruct1, fmstruct2
    TYPE(cp_error_type), INTENT(inout)       :: error
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_struct_compatible', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(fmstruct1),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(fmstruct2),cp_failure_level,routineP,error,failure)
  CPPrecondition(fmstruct1%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(fmstruct2%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     IF (fmstruct1%id_nr==fmstruct2%id_nr) THEN
        res=.TRUE.
     ELSE
        res= (fmstruct1%context%group==fmstruct2%context%group)
     END IF
  END IF
END FUNCTION cp_fm_struct_compatible

! *****************************************************************************
!> \brief writes a description of the matrix structure to the given unit
!> \param fmstruct the structure to print
!> \param unit_nr the number of the unit where to write the description
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE cp_fm_struct_write(fmstruct, unit_nr, long_description,error)
    TYPE(cp_fm_struct_type), POINTER         :: fmstruct
    INTEGER, INTENT(in)                      :: unit_nr
    LOGICAL, INTENT(in), OPTIONAL            :: long_description
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_struct_write', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: fmtstring
    INTEGER                                  :: iostat
    LOGICAL                                  :: failure, l_desc

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
     WRITE (unit=unit_nr,fmt="(' local_leading_dimension=',i8,',')",&
          iostat=iostat) fmstruct%local_leading_dimension
     IF (ASSOCIATED(fmstruct%nrow_locals)) THEN
        WRITE (fmtstring,'(A,i6,A)') "(' nrow_locals=(',",&
             SIZE(fmstruct%nrow_locals),"(i8),'),')"
        WRITE (unit=unit_nr,fmt=fmtstring,iostat=iostat) fmstruct%nrow_locals
     ELSE
        WRITE (unit=unit_nr,fmt="(a)") " nrow_locals=*null*,"
     END IF
     IF (ASSOCIATED(fmstruct%ncol_locals)) THEN
        WRITE (fmtstring,'(A,i6,A)') "(' ncol_locals=(',",SIZE(fmstruct%ncol_locals),&
             "(i8),'),')"
        WRITE (unit=unit_nr,fmt=fmtstring,iostat=iostat) fmstruct%ncol_locals
     ELSE
        WRITE (unit=unit_nr,fmt="(a)") " ncol_locals=*null*,"
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

! *****************************************************************************
!> \brief returns the values of various attributes of the matrix structure
!> \param fmstruct the structure you whant info about
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!>      see cp_fm_struct_type attributes for the other arguments
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE cp_fm_struct_get(fmstruct, para_env, context,&
     descriptor, ncol_block, nrow_block, nrow_global,&
     ncol_global, id_nr, ref_count, first_p_pos, row_indices,&
     col_indices, nrow_local, ncol_local, nrow_locals, ncol_locals,&
     local_leading_dimension, error)
  TYPE(cp_fm_struct_type), POINTER :: fmstruct
  TYPE(cp_para_env_type), POINTER, OPTIONAL :: para_env
  TYPE(cp_blacs_env_type), POINTER, OPTIONAL :: context
  INTEGER, DIMENSION(9), INTENT(OUT), OPTIONAL :: descriptor
  INTEGER, INTENT(out), OPTIONAL :: ncol_block, nrow_block, nrow_global,&
       ncol_global, id_nr, ref_count, nrow_local, ncol_local,&
       local_leading_dimension
  INTEGER, DIMENSION(2), INTENT(out), OPTIONAL :: first_p_pos
  INTEGER, DIMENSION(:), POINTER, OPTIONAL :: row_indices, col_indices,&
       nrow_locals, ncol_locals
  TYPE(cp_error_type), INTENT(inout)  :: error

  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='cp_fm_struct_get',&
       routineP=moduleN//':'//routineN
  INTEGER i,nprow,npcol,myprow,mypcol, stat
#if defined(__SCALAPACK)
  INTEGER , EXTERNAL :: indxl2g
#endif

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
     IF (PRESENT(nrow_locals)) nrow_locals => fmstruct%nrow_locals
     IF (PRESENT(ncol_locals)) ncol_locals => fmstruct%ncol_locals
     IF (PRESENT(local_leading_dimension)) local_leading_dimension =&
          fmstruct%local_leading_dimension

    myprow=fmstruct%context%mepos(1)
    mypcol=fmstruct%context%mepos(2)
    nprow=fmstruct%context%num_pe(1)
    npcol=fmstruct%context%num_pe(2)

    IF (PRESENT(nrow_local)) nrow_local=fmstruct%nrow_locals(myprow)
    IF (PRESENT(ncol_local)) ncol_local=fmstruct%ncol_locals(mypcol)

    IF (PRESENT(row_indices)) THEN
       row_indices => fmstruct%row_indices
       IF (.NOT.ASSOCIATED(row_indices)) THEN
          ! the max should go away
          ALLOCATE(fmstruct%row_indices &
               (MAX(fmstruct%nrow_locals(myprow),1)),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          row_indices => fmstruct%row_indices
#ifdef __SCALAPACK
          DO i=1,SIZE(row_indices)
          ! the max should go away
             fmstruct%row_indices(i)=&
                  indxl2g(i,fmstruct%nrow_block,myprow,&
                  fmstruct%first_p_pos(1),nprow)
          ENDDO
#else
          DO i=1,SIZE(row_indices)
             row_indices(i)=i
          ENDDO
#endif
       END IF
    ENDIF

    IF (PRESENT(col_indices)) THEN
       col_indices => fmstruct%col_indices
       IF (.NOT.ASSOCIATED(col_indices)) THEN
          ALLOCATE(fmstruct%col_indices &
               (MAX(fmstruct%ncol_locals(mypcol),1)),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          col_indices => fmstruct%col_indices
#ifdef __SCALAPACK
          DO i=1,SIZE(col_indices)
             col_indices(i)=&
                  indxl2g(i,fmstruct%ncol_block,mypcol,&
                  fmstruct%first_p_pos(2),npcol)
          ENDDO
#else
          DO i=1,SIZE(col_indices)
             col_indices(i)=i
          ENDDO
#endif
       END IF

    ENDIF
  END IF
END SUBROUTINE cp_fm_struct_get

! *****************************************************************************
!> \brief creates a struct with twice the number of blocks on each core.
!>        If matrix A has to be multiplied with B anc C, a
!>        significant speedup of pdgemm can be acchieved by joining the matrices
!>        in a new one with this structure (see arnoldi in rt_matrix_exp)
!> \param fmstruct the struct to create
!> \param struct struct of either A or B
!> \param col, row  in which direction the matrix should be enlarged
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      06.2009 created [fschiff]
!> \author Florian Schiffmann
! *****************************************************************************
SUBROUTINE cp_fm_struct_double(fmstruct,struct,context,col,row,error)
    TYPE(cp_fm_struct_type), POINTER         :: fmstruct, struct
    TYPE(cp_blacs_env_type), POINTER         :: context
    LOGICAL                                  :: col, row
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_struct_double', &
      routineP = moduleN//':'//routineN

    INTEGER :: nblocks, ncol_block, ncol_global, nempty, newdim_col, &
      newdim_row, nfilled, nprocs_col, nprocs_row, nrow_block, nrow_global
    TYPE(cp_para_env_type), POINTER          :: para_env

    CALL cp_fm_struct_get(struct, nrow_global=nrow_global,&
          ncol_global=ncol_global, nrow_block=nrow_block, &
          ncol_block=ncol_block, error=error)
    newdim_row=nrow_global
    newdim_col=ncol_global
    nprocs_row=context%num_pe(1)
    nprocs_col=context%num_pe(2)
    para_env => struct%para_env

    IF(col)THEN
       IF(ncol_global==0)THEN
          newdim_col=0
       ELSE
          nblocks=CEILING(ncol_global/REAL(ncol_block,dp))
          nfilled=MOD(nblocks,nprocs_col)
          nempty=MOD(nprocs_col-nfilled,nprocs_col)
          newdim_col=2*ncol_global+nempty*ncol_block+2*MOD(ncol_block-MOD(ncol_global,ncol_block),ncol_block)
       END IF
    END IF

    IF(row)THEN
       IF(ncol_global==0)THEN
          newdim_row=0
       ELSE
          nblocks=CEILING(nrow_global/REAL(nrow_block,dp))
          nfilled=MOD(nblocks,nprocs_row)
          nempty=MOD(nprocs_row-nfilled,nprocs_row)
          newdim_row=2*nrow_global+nempty*nrow_block+2*MOD(nrow_block-MOD(nrow_global,nrow_block),nrow_block)
       END IF
    END IF


    CALL cp_fm_struct_create(fmstruct=fmstruct,para_env=para_env,&
         context=context,&
         nrow_global=newdim_row,&
         ncol_global=newdim_col,&
         ncol_block=ncol_block,&
         nrow_block=nrow_block,&
         error=error)

  END SUBROUTINE cp_fm_struct_double


END MODULE cp_fm_struct
