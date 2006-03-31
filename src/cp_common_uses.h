! Common use statements and preprocessor macros
! should be included in the use statements

  USE cp_log_handling
      !  , only: cp_logger_type, cp_iteration_info_type,&
      !  cp_note_level, cp_warning_level, cp_failure_level,&
      !  cp_fatal_level,cp_get_default_logger, cp_add_default_logger, &
      !  cp_rm_default_logger, cp_logger_create, cp_logger_retain, &
      !  cp_logger_release, cp_logger_would_log, cp_logger_set,&
      !  cp_logger_get_default_unit_nr, cp_logger_get_unit_nr, &
      !  cp_logger_set_log_level, cp_logger_generate_filename,&
      !  cp_to_string, cp_log, cp_iteration_info_create, cp_iteration_info_retain, &
      !  cp_iteration_info_release, cp_iter_undefined, max_iter_r_level
  USE cp_error_handling
      ! , only: cp_error_type, cp_debug, cp_no_error, &
      !  cp_caller_error, cp_wrong_args_error,&
      !  cp_precondition_failed, cp_internal_error, cp_postcondition_failed,&
      !  cp_invariant_failed, cp_assertion_failed, cp_unimplemented_error_nr, &
      !  cp_assert, cp_a_l, cp_simple_assert, cp_unimplemented_error, &
      !  cp_error_init, cp_error_dealloc_ref, cp_error_set,&
      !  cp_error_get, cp_error_reset, cp_error_get_level,&
      !  cp_error_get_print_level, cp_error_get_nr,&
      !  cp_error_get_logger, cp_error_get_stop_level,&
      !  cp_error_common_stop, cp_error_handle_error,&
      !  cp_error_message, cp_error_propagate_error,&
      !  cp_error_check, cp_error_synchronize_error
#ifndef CP_PREP_GLOBALS_H
#define CP_PREP_GLOBALS_H

! The following macros are here to facilitate the use of error handling 
! proposed in cp_error_handling.
! they assume at least
! 'use cp_error_handling, only: cp_assert, cp_a_l, cp_error_type'
! and 'use cp_log_handling, only: cp_to_string'
! They ere useful because they give a reference to the file and line
! number in the error message.


! this macro expands to a string that contains the filename.
! if the path is long the filename could make that some lines of code
! become too long and overlow (fortran compilers have a maximum line length)
! in this case substitute __FILE__ with "file" down here.
! obviously then the error messages will not give the filename.
! (otherwise make the file reference in the makefile relative vs. absolute)
#define CPSourceFileRef __FILE__//' line '//cp_to_string(__LINE__)

! if the following macro is defined the longest form of macro 
! expansions is used (and the error messages are more meaningful)

#ifdef FD_LONG_EXPANSIONS

! inlines the test but does not write the file name, but that should
! be easily recovered from the routineP variable (that should contain
! the module name).

! (cp_precondition_failed = 200)
#define CPPrecondition(cond,level,routineP,error,failure) \
IF(.NOT.(cond)) CALL cp_assert(.FALSE.,level,200,routineP,\
'PRECONDITION(cond) failed at line '//cp_to_string(__LINE__),error,failure)

! (cp_postcondition_failed = -200)
#define CPPostcondition(cond,level,routineP,error,failure) \
IF(.NOT.(cond)) CALL cp_assert(.FALSE.,level,-200,routineP,\
'POSTCONDITION(cond) failed at line '//cp_to_string(__LINE__),error,failure)

! (cp_invariant_failed = -100)
#define CPInvariant(cond,level,routineP,error,failure) \
IF(.NOT.(cond)) CALL cp_assert(.FALSE., level, -100,routineP,\
'INVARIANT (cond) failed at line '//cp_to_string(__LINE__),error,failure)

! (cp_assert_failed = -300)
#define CPAssert(cond,level,routineP,error,failure) \
IF(.NOT.(cond)) CALL cp_assert(.FALSE.,level,-300,routineP,\
'ASSERTION (cond) failed at line '//cp_to_string(__LINE__),error,failure)

#define CPErrorMessage(level,routineP,msg,error) \
CALL cp_error_message(level,routineP,msg//' at line '//cp_to_string(__LINE__),error)

#define CPPreconditionNoFail(cond,level,routineP,error) \
IF(.NOT.(cond)) CALL cp_assert(.FALSE.,level,200,routineP,\
'PRECONDITION(cond) failed at line '//cp_to_string(__LINE__),error)

! (cp_postcondition_failed = -200)
#define CPPostconditionNoFail(cond,level,routineP,error) \
IF(.NOT.(cond)) CALL cp_assert(.FALSE.,level,-200,routineP,\
'POSTCONDITION(cond) failed at line '//cp_to_string(__LINE__),error)

! (cp_invariant_failed = -100)
#define CPInvariantNoFail(cond,level,routineP,error) \
IF(.NOT.(cond)) CALL cp_assert(.FALSE., level, -100,routineP,\
'INVARIANT (cond) failed at line '//cp_to_string(__LINE__),error)

! (cp_assert_failed = -300)
#define CPAssertNoFail(cond,level,routineP,error) \
IF(.NOT.(cond)) CALL cp_assert(.FALSE.,level,-300,routineP,\
'ASSERTION (cond) failed at line '//cp_to_string(__LINE__),error)


#define CPPreconditionNoErr(cond, level, routineN) \
IF(.NOT.(cond)) CALL cp_assert(.FALSE., level, 200, routineN ,\
"PRECONDITION failed on line " //cp_to_string(__LINE__))

#define CPPostconditionNoErr(cond, level, routineN) \
IF(.NOT.(cond)) CALL cp_assert(.FALSE., level, -200, routineN ,\
"POSTCONDITION failed on line " //cp_to_string(__LINE__))

#define CPInvariantNoErr(cond,level,routineP) \
IF(.NOT.(cond)) CALL cp_assert(.FALSE., level, -100,routineP,\
'INVARIANT (cond) failed at line '//cp_to_string(__LINE__))

! (cp_assert_failed = -300)
#define CPAssertNoErr(cond,level,routineP) \
IF(.NOT.(cond)) CALL cp_assert(.FALSE.,level,-300,routineP,\
'ASSERTION (cond) failed at line '//cp_to_string(__LINE__))

#else
! same as before but trying to use a small amount of characters
! the test is not inlined (you have always a function call)

#define CPPrecondition(cond,level,routineP,error,failure) \
IF(.NOT.(cond))CALL cp_a_l(0==1,level,routineP,__LINE__,error,failure)
#define CPPostcondition(cond,level,routineP,error,failure) \
IF(.NOT.(cond))CALL cp_a_l(0==1,level,routineP,__LINE__,error,failure)
#define CPInvariant(cond,level,routineP,error,failure) \
IF(.NOT.(cond))CALL cp_a_l(0==1,level,routineP,__LINE__,error,failure)
#define CPAssert(cond,level,routineP,error,failure) \
IF(.NOT.(cond))CALL cp_a_l(0==1,level,routineP,__LINE__,error,failure)
#define CPErrorMessage(level,routineP,msg,error) \
CALL cp_error_message(level,routineP,msg,error)
#define CPPreconditionNoFail(cond,level,routineP,error) \
IF(.NOT.(cond))CALL cp_a_l(0==1,level,routineP,__LINE__,error)
#define CPPostconditionNoFail(cond,level,routineP,error) \
IF(.NOT.(cond))CALL cp_a_l(0==1,level,routineP,__LINE__,error)
#define CPInvariantNoFail(cond,level,routineP,error) \
IF(.NOT.(cond))CALL cp_a_l(0==1,level,routineP,__LINE__,error)
#define CPAssertNoFail(cond,level,routineP,error) \
IF(.NOT.(cond))CALL cp_a_l(0==1,level,routineP,__LINE__,error)
#define CPPreconditionNoErr(cond, level, routineN) \
IF(.NOT.(cond))CALL cp_a_l(0==1,level,routineN,__LINE__)
#define CPPostconditionNoErr(cond, level, routineN) \
IF(.NOT.(cond))CALL cp_a_l(0==1,level,routineN ,__LINE__)
#define CPInvariantNoErr(cond,level,routineP) \
IF(.NOT.(cond))CALL cp_a_l(0==1, level, routineP,__LINE__)
#define CPAssertNoErr(cond,level,routineP) \
IF(.NOT.(cond))CALL cp_a_l(0==1,level,routineP,__LINE__)

#endif

#define CPAAssert(condition) \
CALL cp_simple_assert(condition, routineP, __LINE__)

#define CPTest WRITE(*,*) "#LINE"

#endif
