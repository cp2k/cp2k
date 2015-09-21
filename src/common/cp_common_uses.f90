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
  USE cp_error_handling,               ONLY: cp_debug,&
                                             cp_abort,&
                                             cp_warn,&
                                             cp_assertion_failed,&
                                             cp_internal_error,&
                                             cp_assert,&
                                             cp_unimplemented_error,&
                                             cp_precondition_failed,&
                                             cp_caller_error,&
                                             cp_wrong_args_error,&
                                             cp_unimplemented_error_nr,&
                                             cp__a,cp__b,cp__w

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


#define CPWARN(msg) CALL cp__w(routineP,__LINE__,msg)
#define CPABORT(msg) CALL cp__b(routineP,__LINE__,msg)
#define CPASSERT(cond) IF(.NOT.(cond))CALL cp__a(routineP,__LINE__)

