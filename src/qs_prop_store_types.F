!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****h* cp2k/qs_prop_store_types *
!!
!!   NAME
!!     qs_prop_store_types
!!
!!   FUNCTION
!!     types to control calcualtion ,ouput and in memory storage of the 
!!     properties that can be calculated by qs
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     08.2005 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE qs_prop_store_types
  use kinds, only: dp
  use cp_para_types, only: cp_para_env_type
  use timings, only: timeset, timestop
  USE cp_output_handling, only: cp_out_default
#include "cp_common_uses.h"

  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='qs_prop_store_types'
  INTEGER, PRIVATE, SAVE :: last_pstore_id_nr=0
!!***
!****************************************************************************


  !!****s* qs_prop_store_types/qs_prop_store_type *
  !!
  !!   NAME
  !!     qs_prop_store_type
  !!
  !!   FUNCTION
  !!     stores the properties of the qs_env
  !!
  !!   NOTES
  !!     -
  !!
  !!   ATTRIBUTES
  !!     -
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!   SOURCE
  !***************************************************************************
  TYPE qs_prop_store_type
     INTEGER :: id_nr, ref_count
     INTEGER, DIMENSION(max_iter_r_level) :: last_wannier_centers, &
          last_wannier_wfn, last_atom_pos, last_dipole
     REAL(kind=dp), DIMENSION(:,:), POINTER :: wannier_centers, &
          wannier_spreads
     REAL(kind=dp), DIMENSION(3) :: total_dipole,wrapped_total_dipole
  END TYPE qs_prop_store_type
  !!***
  !****************************************************************************

contains


!!****f* qs_prop_store_types/qs_pstore_create *
!!
!!   NAME
!!     qs_pstore_create
!!
!!   FUNCTION
!!     allocates and initializes a property store
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - pstore: the property store to create
!!     - control_flags: flags to use (defaults to all cp_out_default)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE qs_pstore_create(pstore,error)
  TYPE(qs_prop_store_type), pointer :: pstore
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='qs_pstore_create',&
        routineP=moduleN//':'//routineN
  integer :: stat

  failure=.false.
  
  CPPrecondition(.NOT.ASSOCIATED(pstore),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     ALLOCATE(pstore,stat=stat)
     CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
  END IF
  IF(.NOT.failure) THEN
     last_pstore_id_nr=last_pstore_id_nr+1
     pstore%id_nr=last_pstore_id_nr
     pstore%ref_count=1

     pstore%last_wannier_centers=cp_iter_undefined
     pstore%last_wannier_wfn=cp_iter_undefined
     pstore%last_atom_pos=cp_iter_undefined
     pstore%last_dipole=cp_iter_undefined

     nullify(pstore%wannier_centers, pstore%wannier_spreads)
     pstore%total_dipole=0._dp
     pstore%wrapped_total_dipole=0._dp
  END IF
END SUBROUTINE qs_pstore_create
!***************************************************************************

!!****f* qs_out_flags_types/qs_pstore_retain *
!!
!!   NAME
!!     qs_pstore_retain
!!
!!   FUNCTION
!!     retains the given pstore (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - pstore: the pstore to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE qs_pstore_retain(pstore,error)
  TYPE(qs_prop_store_type), pointer :: pstore
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='qs_pstore_retain',&
        routineP=moduleN//':'//routineN

  failure=.false.
  
  CPPrecondition(ASSOCIATED(pstore),cp_failure_level,routineP,error,failure)
  if (.not. failure) then
     CPPrecondition(pstore%ref_count>0,cp_failure_level,routineP,error,failure)
     pstore%ref_count=pstore%ref_count+1
  end if
END SUBROUTINE qs_pstore_retain
!***************************************************************************

!!****f* qs_out_flags_types/qs_pstore_release *
!!
!!   NAME
!!     qs_pstore_release
!!
!!   FUNCTION
!!     releases the given pstore (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - pstore: the pstore to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE qs_pstore_release(pstore,error)
  TYPE(qs_prop_store_type), pointer :: pstore
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='qs_pstore_release',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat

  failure=.false.
  
  IF (ASSOCIATED(pstore)) THEN
     CPPrecondition(pstore%ref_count>0,cp_failure_level,routineP,error,failure)
     pstore%ref_count=pstore%ref_count-1
     IF (pstore%ref_count==0) THEN
        IF (ASSOCIATED(pstore%wannier_centers)) THEN
           DEALLOCATE(pstore%wannier_centers,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(pstore%wannier_spreads)) THEN
           DEALLOCATE(pstore%wannier_spreads,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        DEALLOCATE(pstore,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
END SUBROUTINE qs_pstore_release
!***************************************************************************

END MODULE qs_prop_store_types
