!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_random_types [1.0] *
!!
!!   NAME
!!     cp_random_types
!!
!!   FUNCTION
!!     flexible random number generator
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     07.2004 created
!!
!!   SOURCE
!****************************************************************************
MODULE cp_random_types
  use cp_log_handling, only: cp_fatal_level, cp_failure_level,&
       cp_warning_level, cp_note_level, cp_to_string, cp_log
  use cp_error_handling, only: cp_debug, cp_error_type, cp_error_init,&
       cp_error_dealloc_ref, cp_error_message, cp_assert,&
       cp_assertion_failed, cp_internal_error, cp_unimplemented_error,&
       cp_error_get_logger, cp_a_l
  USE kinds, ONLY: dp, default_string_length
  use cp_para_types, only: cp_para_env_type
  use timings, only: timeset, timestop
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='cp_random_types'
  INTEGER, SAVE, private :: last_rg_id_number=0

  INTEGER, PUBLIC, PARAMETER :: rg_luxury_0=1, rg_luxury_1=2,rg_luxury_2=3,&
       rg_luxury_3=4,rg_luxury_4=5, rg_fortran_default=6

  public cp_random_generator_type

!***
!****************************************************************************
  
  !!****s* cp_random_generators/cp_random_generator_type [1.0] *
  !!
  !!   NAME
  !!     cp_random_generator_type
  !!
  !!   FUNCTION
  !!     represents a random generator (with its state and method used)
  !!
  !!   NOTES
  !!     -
  !!
  !!   ATTRIBUTES
  !!     - id_nr: identification number (different in each instance)
  !!     - ref_count: reference count (see doc/ReferenceCounting.html)
  !!     - seed_size: size of the seed (number of integers)
  !!     - method_id: method used to generate the random numbers
  !!     - i1_status: integer array that stores the actual status of
  !!       the random number generator
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!   SOURCE
  !***************************************************************************
  type cp_random_generator_type
     INTEGER :: id_nr,ref_count,seed_size,method_size
     INTEGER, DIMENSION(:), POINTER :: i1_status
  end type cp_random_generator_type
  !!***
  !****************************************************************************

contains

!!****f* cp_random_types/rg_method_id_2_label *
!!
!!   NAME
!!     rg_method_id_2_label
!!
!!   FUNCTION
!!     returns the label of the method id
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - method_id: the id of the method
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
  FUNCTION rg_method_id_2_label(method_id,error) RESULT(res)
    INTEGER, intent(in) :: method_id
    type(cp_error_type), optional, intent(inout) :: error
    CHARACTER(len=default_string_length) :: res

    logical :: failure
    character(len=*), parameter :: routineN='rg_method_id_2_label',&
         routineP=moduleN//':'//routineN

    SELECT CASE(method_id)
    case (rg_luxury_0)
       res="LUXURY_L0"
    CASE (rg_luxury_1)
       res="LUXURY_L1"
    CASE(rg_luxury_2)
       res="LUXURY_L2"
    CASE(rg_luxury_3)
       res="LUXURY_L3"
    CASE(rg_luxury_4)
       res="LUXURY_L4"
    CASE(rg_fortran_default)
       res="FORTRAN_DEFAULT"
    CASE default
       res="UNKNOWN_"//TRIM(adjustl(cp_to_string(method_id)))
       CPAssertNoFail(.FALSE.,cp_failure_level,routineP,error)
    END SELECT
  END FUNCTION rg_method_id_2_label
!***************************************************************************

! rg_get(rg,method_id,id_nr,ref_count,seed_size,method_size)
! rg_retain/release

END MODULE cp_random_types
