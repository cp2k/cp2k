! Common use statements and preprocessor macros
! should be included in the use statements

  use cp_log_handling, only: cp_logger_type, cp_iteration_info_type,&
        cp_note_level, cp_warning_level, cp_failure_level,&
        cp_fatal_level,cp_get_default_logger, cp_add_default_logger, &
        cp_rm_default_logger, cp_logger_create, cp_logger_retain, &
        cp_logger_release, cp_logger_would_log, cp_logger_set,&
        cp_logger_get_default_unit_nr, cp_logger_get_unit_nr, &
        cp_logger_set_log_level, cp_logger_generate_filename,&
        cp_to_string, cp_log, cp_iteration_info_create, cp_iteration_info_retain, &
        cp_iteration_info_release
  use cp_error_handling, only: cp_error_type, cp_debug, cp_no_error, &
        cp_caller_error, cp_wrong_args_error,&
        cp_precondition_failed, cp_internal_error, cp_postcondition_failed,&
        cp_invariant_failed, cp_assertion_failed, cp_unimplemented_error_nr, &
        cp_assert, cp_a_l, cp_simple_assert, cp_unimplemented_error, &
        cp_error_init, cp_error_dealloc_ref, cp_error_set,&
        cp_error_get, cp_error_reset, cp_error_get_level,&
        cp_error_get_print_level, cp_error_get_nr,&
        cp_error_get_logger, cp_error_get_stop_level,&
        cp_error_common_stop, cp_error_handle_error,&
        cp_error_message, cp_error_propagate_error,&
        cp_error_check, cp_error_synchronize_error

#include "cp_prep_globals.h"
