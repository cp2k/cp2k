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
module cp_f_matrix_struct
  use cp_log_handling, only: cp_fatal_level, cp_failure_level,&
       cp_warning_level, cp_note_level, cp_to_string, cp_log,&
       cp_logger_get_default_unit_nr
  use cp_error_handling, only: cp_debug, cp_error_type, cp_error_init,&
       cp_error_dealloc_ref, cp_error_message, cp_assert,&
       cp_assertion_failed, cp_internal_error, cp_unimplemented_error,&
       cp_error_get_logger
  use kinds, only: wp=>dp
  use cp_para_types, only: cp_para_env_type, cp_blacs_env_type
  use cp_para_env, only: cp_para_env_retain, cp_para_env_release, cp_para_env_write
  use cp_blacs_env, only: cp_blacs_env_retain, cp_blacs_env_release,&
       cp_blacs_env_write
  use timings, only: timeset, timestop
  USE message_passing,     ONLY: mp_bcast,&
                                 mp_max,&
                                 mp_sum,&
                                 mp_sync
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='cp_f_matrix_struct'
  integer,save,private :: last_fmstruct_id_nr=0

  public :: cp_f_matrix_struct_type
  public :: cp_fmstruct_create,cp_fmstruct_retain, cp_fmstruct_release,&
       cp_fmstruct_equivalent
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
  type cp_f_matrix_struct_type
     type(cp_para_env_type), pointer :: para_env
     type(cp_blacs_env_type), pointer :: context
     integer, dimension(9) :: descriptor
     ! when using blacs
     ! ncol_block=descriptor(6)
     ! ncol_global=descriptor(4)
     ! nrow_block=descriptor(5)
     ! nrow_global=descriptor(3)
     integer :: ncol_block, nrow_block, nrow_global, ncol_global
     integer, dimension(2) :: first_p_pos
     integer, dimension(:), pointer :: row_indices, col_indices,&
          nrow_local, ncol_local
     integer :: id_nr, ref_count
  end type cp_f_matrix_struct_type
  !!***
  !****************************************************************************

contains

!!****f* qs_blacs/cp_fmstruct_create [1.0] *
!!
!!   NAME
!!     cp_fmstruct_create
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
subroutine cp_fmstruct_create(fmstruct,para_env,nrow_global,ncol_global,&
     nrow_block, ncol_block, context, descriptor, first_p_pos, error)
  type(cp_f_matrix_struct_type), pointer :: fmstruct
  type(cp_para_env_type), pointer :: para_env
  integer, intent(in) :: nrow_global,ncol_global,nrow_block,ncol_block
  type(cp_blacs_env_type), pointer :: context
  integer, dimension(9), intent(in), optional :: descriptor
  integer,optional, dimension(2) :: first_p_pos
  type(cp_error_type), optional, intent(inout) :: error
  
#if defined(__SCALAPACK)
    INTEGER, EXTERNAL :: numroc
#endif

  logical :: failure
  character(len=*), parameter :: routineN='cp_fmstruct_create',&
        routineP=moduleN//':'//routineN
  integer :: stat
  failure=.false.
  
  allocate(fmstruct,stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  if (.not. failure) then
     fmstruct%context => context
     call cp_blacs_env_retain(fmstruct%context,error=error)
     fmstruct%para_env => para_env
     call cp_para_env_retain(para_env,error=error)
     
     fmstruct%nrow_global=nrow_global
     fmstruct%ncol_global=ncol_global
     fmstruct%nrow_block=min(nrow_block,nrow_global/context%num_pe(1))
     fmstruct%ncol_block=min(ncol_block,ncol_global/context%num_pe(2))
     ! checks for more processes than matrix elements
     CPAssert(fmstruct%nrow_block>0,cp_failure_level,routineP,error,failure)
     CPAssert(fmstruct%ncol_block>0,cp_failure_level,routineP,error,failure)
     if (fmstruct%nrow_global==fmstruct%ncol_global) then
        ! square matrix -> square blocks (otherwise some op fail)
        fmstruct%nrow_block=min(fmstruct%nrow_block,fmstruct%ncol_block)
        fmstruct%ncol_block=fmstruct%nrow_block
     end if
     allocate(fmstruct%nrow_local(0:(context%num_pe(1)-1)),&
          fmstruct%ncol_local(0:(context%num_pe(2)-1)),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  end if
  if (.not.failure) then
     fmstruct%first_p_pos=context%source
     if (present(first_p_pos)) fmstruct%first_p_pos=first_p_pos

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
     call mp_sum(fmstruct%nrow_local,fmstruct%para_env%group)
     call mp_sum(fmstruct%ncol_local,fmstruct%para_env%group)
#else
     ! block = full matrix
     fmstruct%nrow_block=fmstruct%nrow_global
     fmstruct%ncol_block=fmstruct%ncol_global
     fmstruct%nrow_local(context%mepos(1))=fmstruct%nrow_global
     fmstruct%ncol_local(context%mepos(2))=fmstruct%ncol_global
#endif

     nullify(fmstruct%row_indices, fmstruct%col_indices)
     last_fmstruct_id_nr=last_fmstruct_id_nr+1
     fmstruct%id_nr=last_fmstruct_id_nr
     fmstruct%ref_count=1

     if (present(descriptor)) then
        fmstruct%descriptor=descriptor
     else
        fmstruct%descriptor=0
#if defined(__SCALAPACK)
        CALL descinit(fmstruct%descriptor,fmstruct%nrow_global,&
             fmstruct%ncol_global,fmstruct%nrow_block,&
             fmstruct%ncol_block,fmstruct%first_p_pos(1),&
             fmstruct%first_p_pos(2),fmstruct%context%group,&
             fmstruct%nrow_local(fmstruct%context%mepos(1)),stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
#endif
     end if
  end if
end subroutine cp_fmstruct_create
!***************************************************************************

!!****f* qs_blacs/cp_fmstruct_retain [1.0] *
!!
!!   NAME
!!     cp_fmstruct_retain
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
subroutine cp_fmstruct_retain(fmstruct,error)
  type(cp_f_matrix_struct_type), pointer :: fmstruct
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_fmstruct_retain',&
        routineP=moduleN//':'//routineN
  failure=.false.
  
  CPPrecondition(associated(fmstruct),cp_failure_level,routineP,error,failure)
  if (.not. failure) then
     CPPrecondition(fmstruct%ref_count>0,cp_failure_level,routineP,error,failure)
     fmstruct%ref_count=fmstruct%ref_count+1
  end if
end subroutine cp_fmstruct_retain
!***************************************************************************

!!****f* qs_blacs/cp_fmstruct_release [1.0] *
!!
!!   NAME
!!     cp_fmstruct_release
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
subroutine cp_fmstruct_release(fmstruct,error)
  type(cp_f_matrix_struct_type), pointer :: fmstruct
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_fmstruct_release',&
        routineP=moduleN//':'//routineN
  integer :: stat
  failure=.false.
  
  if (associated(fmstruct)) then
     CPPrecondition(fmstruct%ref_count>0,cp_failure_level,routineP,error,failure)
     fmstruct%ref_count=fmstruct%ref_count-1
     if (fmstruct%ref_count<1) then
        call cp_blacs_env_release(fmstruct%context,error=error)
        call cp_para_env_release(fmstruct%para_env,error=error)
        if (associated(fmstruct%row_indices)) then
           deallocate(fmstruct%row_indices,stat=stat)
           CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
        endif
        if (associated(fmstruct%col_indices)) then
           deallocate(fmstruct%col_indices,stat=stat)
           CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
        end if
        deallocate(fmstruct, stat=stat)
        CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
     end if
  end if
  nullify(fmstruct)
end subroutine cp_fmstruct_release
!***************************************************************************

!!****f* qs_blacs/cp_fmstruct_equivalent [1.0] *
!!
!!   NAME
!!     cp_fmstruct_equivalent
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
function cp_fmstruct_equivalent(fmstruct1,fmstruct2,error) result(res)
  type(cp_f_matrix_struct_type), pointer :: fmstruct1, fmstruct2
  type(cp_error_type), optional, intent(inout) :: error
  logical :: res
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_fmstruct_equivalent',&
        routineP=moduleN//':'//routineN
  integer :: i
  failure=.false.

  CPPrecondition(associated(fmstruct1),cp_failure_level,routineP,error,failure)
  CPPrecondition(associated(fmstruct2),cp_failure_level,routineP,error,failure)
  CPPrecondition(fmstruct1%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(fmstruct2%ref_count>0,cp_failure_level,routineP,error,failure)
  if (.not. failure) then
     if (fmstruct1%id_nr==fmstruct2%id_nr) then
        res=.true.
     else
        res= (fmstruct1%context%group==fmstruct2%context%group).and.&
             (fmstruct1%nrow_global==fmstruct2%nrow_global).and.&
             (fmstruct1%ncol_global==fmstruct2%ncol_global)
        do i=1,9
           res=res.and.(fmstruct1%descriptor(i)==fmstruct1%descriptor(i))
        end do
     end if
  end if
end function cp_fmstruct_equivalent
!***************************************************************************

end module cp_f_matrix_struct
