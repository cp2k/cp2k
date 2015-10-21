! Common use statements and preprocessor macros
! should be included in the use statements
  USE cp_log_handling,                 ONLY: cp_logger_type,&
                                             cp_logger_p_type,&
                                             cp_logger_create,&
                                             cp_logger_release,&
                                             cp_logger_set,&
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
