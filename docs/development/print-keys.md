# Printkeys

When writing output in cp2k, to the screen or a file, one should use a printkey.

A printkey allows the user to control the output such as turning if off, redirecting it to a
different file, or outputting it after a certain number of iteration.

To each printkey belongs an input-section, which is always structured in a similar way. For example,
the printkey used to write trajectories is controlled by
[this input section](#CP2K_INPUT.MOTION.PRINT.TRAJECTORY).

## Creating a Printkey Input-section

```
 USE cp_output_handling,              ONLY: cp_print_key_section_create
 !...
 CALL section_create(subsection,name="print",&
      description="Controls the printing properties during an MD run",&
      n_keywords=0, n_subsections=1, repeats=.TRUE., required=.FALSE.,error=error)
 CALL cp_print_key_section_create(print_key,"TRAJECTORY",&
      description="Controls the output of the trajectory",&
      print_level=low_print_level, common_iter_levels=1,&
      filename="",unit_str="angstrom",error=error)
 CALL section_add_subsection(subsection,print_key,error=error)
 CALL section_release(print_key,error=error)
 !...
```

## Using a Printkey for output

```
 USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                            cp_print_key_unit_nr
 !...
 TYPE(cp_logger_type), POINTER :: logger
 INTEGER                       :: traj_unit
 !...
 logger => cp_error_get_logger(error)
 traj_unit = cp_print_key_unit_nr(logger,root_section,"MOTION%PRINT%TRAJECTORY",error=error)
 IF (traj_unit > 0) THEN
   WRITE(traj_unit, *) "here goes the output"
 ENDIF
 CALL cp_print_key_finished_output(traj_unit,logger,root_section,"MOTION%PRINT%TRAJECTORY",error=error)
```
