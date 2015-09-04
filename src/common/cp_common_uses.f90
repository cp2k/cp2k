! Common use statements and preprocessor macros
! should be included in the use statements

  USE cp_log_handling,                 ONLY: cp_logger_type,&
                                             cp_logger_p_type,&
                                             cp_logger_create,&
                                             cp_logger_release,&
                                             cp_logger_set,&
                                             cp_log,&
                                             cp_warning_level,&
                                             cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_logger_get_default_io_unit,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_default_logger_stack_size,&
                                             cp_logger_get_unit_nr,&
                                             cp_logger_generate_filename,&
                                             cp_get_default_logger,&
                                             cp_logger_would_log,&
                                             cp_add_default_logger,&
                                             cp_rm_default_logger,&
                                             cp_to_string
  USE cp_error_handling,               ONLY: cp_error_type,&
                                             cp_debug,&
                                             cp_assertion_failed,&
                                             cp_internal_error,&
                                             cp_assert,&
                                             cp_unimplemented_error,&
                                             cp_error_get_logger,&
                                             cp_precondition_failed,&
                                             cp_error_init,&
                                             cp_caller_error,&
                                             cp_wrong_args_error,&
                                             cp_unimplemented_error_nr,&
                                             cp_a_l,&
                                             cp_error_dealloc_ref,&
                                             cp_error_message

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
#ifdef __SHORT_FILE__
#define CPSourceFileRef __SHORT_FILE__//' line '//TRIM(ADJUSTL(cp_to_string(__LINE__)))
#else
#define CPSourceFileRef __FILE__//' line '//TRIM(ADJUSTL(cp_to_string(__LINE__)))
#endif

! if the following macro is defined the longest form of macro
! expansions is used (and the error messages are more meaningful)

! inlines the test but does not write the file name, but that should
! be easily recovered from the routineP variable (that should contain
! the module name).
!
! We are trying to use a small amount of characters
! the test is not inlined (you have always a function call)

#define CPPrecondition(cond,level,routineP,error,failure) \
IF(.NOT.(cond))CALL cp_a_l(level,routineP,__LINE__,error,failure)
#define CPPostcondition(cond,level,routineP,error,failure) \
IF(.NOT.(cond))CALL cp_a_l(level,routineP,__LINE__,error,failure)
#define CPInvariant(cond,level,routineP,error,failure) \
IF(.NOT.(cond))CALL cp_a_l(level,routineP,__LINE__,error,failure)
#define CPAssert(cond,level,routineP,error,failure) \
IF(.NOT.(cond))CALL cp_a_l(level,routineP,__LINE__,error,failure)
#define CPErrorMessage(level,routineP,msg,error) \
CALL cp_error_message(level,routineP,msg,error)
#define CPPreconditionNoFail(cond,level,routineP,error) \
IF(.NOT.(cond))CALL cp_a_l(level,routineP,__LINE__,error)
#define CPPostconditionNoFail(cond,level,routineP,error) \
IF(.NOT.(cond))CALL cp_a_l(level,routineP,__LINE__,error)
#define CPInvariantNoFail(cond,level,routineP,error) \
IF(.NOT.(cond))CALL cp_a_l(level,routineP,__LINE__,error)
#define CPAssertNoFail(cond,level,routineP,error) \
IF(.NOT.(cond))CALL cp_a_l(level,routineP,__LINE__,error)
#define CPPreconditionNoErr(cond, level, routineN) \
IF(.NOT.(cond))CALL cp_a_l(level,routineN,__LINE__)
#define CPPostconditionNoErr(cond, level, routineN) \
IF(.NOT.(cond))CALL cp_a_l(level,routineN ,__LINE__)
#define CPInvariantNoErr(cond,level,routineP) \
IF(.NOT.(cond))CALL cp_a_l(level,routineP,__LINE__)
#define CPAssertNoErr(cond,level,routineP) \
IF(.NOT.(cond))CALL cp_a_l(level,routineP,__LINE__)

! new assert macro, meant to replace all of the above
#define ASSERT(cond) \
IF(.NOT.(cond))CALL cp_a_l(CP_FAILURE_LEVEL,routineP,__LINE__)

