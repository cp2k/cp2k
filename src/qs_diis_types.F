!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/qs_diis_types [1.0] *
!!
!!   NAME
!!     qs_diis_types
!!
!!   FUNCTION
!!     buffer for the diis of the scf
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Matthias Krack
!!
!!   MODIFICATION HISTORY
!!     02.2003 rewamped [fawzi]
!!
!!   SOURCE
!****************************************************************************
module qs_diis_types
  use cp_log_handling, only: cp_fatal_level, cp_failure_level,&
       cp_warning_level, cp_note_level, cp_to_string, cp_log
  use cp_error_handling, only: cp_debug, cp_error_type, cp_error_init,&
       cp_error_dealloc_ref, cp_error_message, cp_assert,&
       cp_assertion_failed, cp_internal_error, cp_unimplemented_error,&
       cp_error_get_logger, cp_a_l, cp_error_check
  USE cp_fm_types, ONLY: cp_fm_p_type, cp_fm_type,&
       cp_fm_release, cp_fm_retain
  use kinds, only: dp
  use cp_para_types, only: cp_para_env_type
  use timings, only: timeset, timestop
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='qs_diis_types'

  public :: qs_diis_buffer_type
  public :: qs_diis_b_retain, qs_diis_b_release
!!***
!****************************************************************************

  !!****s* qs_diis_types/qs_diis_buffer_type [1.0] *
  !!
  !!   NAME
  !!     qs_diis_buffer_type
  !!
  !!   FUNCTION
  !!     keeps a buffer with the previous values of s,p,k
  !!
  !!   NOTES
  !!     -
  !!
  !!   ATTRIBUTES
  !!     -
  !!
  !!   AUTHOR
  !!     Matthias Krack
  !!
  !!   MODIFICATION HISTORY
  !!     02.2003 rewamped [fawzi]
  !!
  !!   SOURCE
  !***************************************************************************
  TYPE qs_diis_buffer_type
    INTEGER                                          :: nbuffer,ncall,&
         id_nr,ref_count
    TYPE(cp_fm_p_type), DIMENSION(:,:), POINTER :: error,PARAMETER
    REAL(KIND = dp), DIMENSION(:,:), POINTER                :: b_matrix
  END TYPE qs_diis_buffer_type
  !!***
  !****************************************************************************

  !!****s* qs_diis_types/qs_diis_buffer_p_type [1.0] *
  !!
  !!   NAME
  !!     qs_diis_buffer_p_type
  !!
  !!   FUNCTION
  !!     build array of pointers to diis buffers
  !!
  !!   NOTES
  !!     -
  !!
  !!   ATTRIBUTES
  !!     - diis_buffer: the diis buffer pointer
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!   MODIFICATION HISTORY
  !!     02.2003 created [fawzi]
  !!
  !!   SOURCE
  !***************************************************************************
  TYPE qs_diis_buffer_p_type
     TYPE(qs_diis_buffer_type), POINTER :: diis_buffer
  END TYPE qs_diis_buffer_p_type
  !!***
  !****************************************************************************

CONTAINS

!!****f* qs_diis_types/qs_diis_b_retain [1.0] *
!!
!!   NAME
!!     qs_diis_b_retain
!!
!!   FUNCTION
!!     retains a diis buffer (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - diis_buffer: the buffer to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE qs_diis_b_retain(diis_buffer,error)
  type(qs_diis_buffer_type), pointer :: diis_buffer
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='qs_diis_b_retain',&
        routineP=moduleN//':'//routineN

  failure=.false.
  
  CPPrecondition(associated(diis_buffer),cp_failure_level,routineP,error,failure)
  if (.not. failure) then
     CPPreconditionNoFail(diis_buffer%ref_count>0,cp_failure_level,routineP,error)
     diis_buffer%ref_count=diis_buffer%ref_count+1
  end if
END SUBROUTINE qs_diis_b_retain
!***************************************************************************

!!****f* qs_diis_types/qs_diis_b_release [1.0] *
!!
!!   NAME
!!     qs_diis_b_release
!!
!!   FUNCTION
!!     releases the given diis buffer (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - diis_buffer: the buffer to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!*** **********************************************************************
subroutine qs_diis_b_release(diis_buffer,error)
  type(qs_diis_buffer_type), pointer :: diis_buffer  
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='qs_diis_b_release',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat,i,j

  failure=.false.
  
  if (associated(diis_buffer)) then
     CPPreconditionNoFail(diis_buffer%ref_count>0,cp_failure_level,routineP,error)
     diis_buffer%ref_count=diis_buffer%ref_count-1
     if (diis_buffer%ref_count<1) then
        if (associated(diis_buffer%b_matrix)) then
           deallocate(diis_buffer%b_matrix,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        end if
        if (associated(diis_buffer%error)) then
           do j=1,size(diis_buffer%error,2)
              do i=1,size(diis_buffer%error,1)
                 call cp_fm_release(diis_buffer%error(i,j)%matrix,error=error)
              end do
           end do
           deallocate(diis_buffer%error,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        end if
        if (associated(diis_buffer%parameter)) then
           do j=1,size(diis_buffer%parameter,2)
              do i=1,size(diis_buffer%parameter,1)
                 call cp_fm_release(diis_buffer%parameter(i,j)%matrix,&
                      error=error)
              end do
           end do
           deallocate(diis_buffer%parameter,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        end if
        deallocate(diis_buffer,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     end if
  end if
end subroutine qs_diis_b_release
!***************************************************************************

end module qs_diis_types
