!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_random_methods [1.0] *
!!
!!   NAME
!!     cp_random_methods
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
MODULE cp_random_methods
  use cp_log_handling, only: cp_fatal_level, cp_failure_level,&
       cp_warning_level, cp_note_level, cp_to_string, cp_log
  use cp_error_handling, only: cp_debug, cp_error_type, cp_error_init,&
       cp_error_dealloc_ref, cp_error_message, cp_assert,&
       cp_assertion_failed, cp_internal_error, cp_unimplemented_error,&
       cp_error_get_logger, cp_a_l
  use kinds, only: dp
  use cp_para_types, only: cp_para_env_type
  use timings, only: timeset, timestop
  USE cp_random_types, ONLY: cp_random_generator_type, rg_luxury_0, &
       rg_luxury_1,rg_luxury_3,rg_luxury_4,rg_fortran_default
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='cp_random_methods'
  INTEGER, SAVE, private :: last_rg_id_number=0

  LOGICAL, SAVE, private :: default_rg_initialized=.false.
  TYPE(cp_random_generator_type), POINTER, SAVE, PRIVATE :: default_generator
!***
!****************************************************************************
contains

!!****f* cp_random_methods/rg_create *
!!
!!   NAME
!!     rg_create
!!
!!   FUNCTION
!!     creates a random generator
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - rg: will contain the random number generator
!!     - method_id: the id of the method to use (defaults to rg_luxury_l4)
!!     - seed: the seed to initialize the random generator
!!       (the default generator always accept also an array of size 1)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE rg_create(rg,method_id,seed,error)
  TYPE(cp_random_generator_type), pointer :: rg
  INTEGER, INTENT(in), optional :: method_id
  INTEGER, DIMENSION(:), INTENT(in), optional :: seed
  
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='rg_create',&
        routineP=moduleN//':'//routineN
  failure=.false.
  
  if (.not. failure) then
     
  end if
END SUBROUTINE rg_create
!***************************************************************************

!rg_create_copy(orig_rg,copy_rg) ! to store the status
!rg_get_state(rg,old_state) ! old_state is a cp_random_generator_type
!
!rg_seed(rg,seed,old_state)
!rg_seed_random(rg)
!rg_skip(rg,n)
!
!rg_get_int(rg,lower_bound,upper_bound) ! returns int in [lower_bound,upper_bound]
!rg_get_free_int(rg,lower_bound,upper_bound) ! lower_bound, upper_nound are optional *OUTPUT* values
!rg_get_real(rg)! returns a real in (0,1)
!! 1d array version
!rg_reals_get(rg,reals)
!rg_ints_get(rg,ints,lower_bound,upper_bound)
!rg_free_ints_get(rg,ints,lower_bound,upper_bound)
!
!! io functions
!rg_read_state(rg,section)
!rg_write_state(rg,unit_nr) ! write section_vals?
!
!! default rg
!rg_get_default_rg()
!get_int(lower_bound,upper_bound)
!get_free_int(lower_bound,upper_bound) ! lower_bound, upper_nound are optional *OUTPUT* values
!get_real()! returns a real in (0,1), unless one explictly chooses a generator with a non uniform prob.dist.
!! 1d array version
!reals_get(reals)
!ints_get(ints,lower_bound,upper_bound)
!free_ints_get(ints,lower_bound,upper_bound)

END MODULE cp_random_methods
