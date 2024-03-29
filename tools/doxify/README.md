# DOXYGEN

This file describes the scripts used to add missing comments and parse existing comments for
subroutines and functions in CP2K. The objective of the scripts is to ensure that all subroutines
and functions in CP2K have a standard DOXYGEN comment block and to flag any missing comments such
that the CP2K developers can add/update these when time is available.

The driver script `doxyify.sh` will process CP2K source code (`*.F` files) processing one file at a
time. The following steps are carried out for each `*.F` file:

1. Run `remove_double_ampersands.pl` - removes any any double ampersand characters
1. Run `fixcomments.pl` - does the addition of missing comments or updating of existing comments
1. Run `remove_extra_comments.pl` - removes any extra comment lines which arise from application of
   `fixcomments.pl`
1. Finally, overwrite the original `*.F` file with the updated version

To run the script, you should execute it with the name of the source file you want to process:

```shell
./doxyify.sh full_path_to_cp2k_src/myfile.F
```

or with a list of files, e.g.

```shell
./doxyify.sh full_path_to_cp2k_src/myfile.F full_path_to_cp2k_src/myfile2.F full_path_to_cp2k_src/myfile2.F
```

In debug mode, the script will produce a lot of output if the whole source tree is processed so you
may wish to redirect the stdout to file, e.g.

```shell
./doxyify.sh > comments.txt
```

## `remove_double_ampersands.pl`

Removes any double ampersand characters that occur at the end of one line and the start of the next
line, e.g:

```fortran
SUBROUTINE fred(a,b, &
   & c,d)
```

becomes

```fortran
SUBROUTINE fred(a,b, &
     c,d)
```

This is performed throughout the source code and is necessary as the ampersand is used by
`fixcomments.pl` to determine whether a subroutine/function definition extends across multiple
lines.

## `fixcomments.pl`

The script adds in missing comments and parses existing comments checking whether any data is
missing. For procedures with no existing comments, e.g. a subroutine with 3 arguments `a`, `b` and
`c`:

```fortran
SUBROUTINE fred (a,b,c)
```

will be updated to (the ... are to indicate comments are required)

```fortran
! *****************************************************************************
!> \brief ...
!> \param a ...
!> \param b ...
!> \param c ...
! *****************************************************************************
SUBROUTINE fred (a,b,c)
```

The scripts works as follows:

- We read through the source file until we see a comment block e.g. something beginning with `!>`

- We then loop over the comment block looking for entries for `\brief`, `\param`, `\date`, `\par`,
  `\version`, `\note`, and `\author`. If any of these items exist we store the details in separate
  variables. We also keep a copy of the entire comment block in the `oldheader` variable as we need
  to print the comment block out unchanged for MODULE and TYPE. When parsing the header we ensure
  that we don't throw any text away.

- Once the comment block has been read in we continue reading in the code line by line. If we don't
  encounter any procedures then the output file will be identical to the input file.

- However, if we encounter a procedure we then read in the line(s) containing the definition, e.g.

  ```fortran
  SUBROUTINE fred(a,b,c)
  ```

  to an array called `@string`. We do this by splitting over space, comma, brackets etc. We also
  take into account of the fact a subroutine / function definition can extend over multiple lines
  via the ampersand (`&`) character and ampersand variable. The `@string` array then contains the
  type of procedure, its name and the arguments, e.g. for the example above `@string` would contain
  5 elements: `SUBROUTINE`, `fred`, `a`, `b`, `c`.

- Once we know the name and arguments of the procedure we loop over the arguments checking to see
  whether the comment block contained a match. If it is does we use the existing text, if no match
  is found then the standard text detailed above is output. We also check whether any text exists
  for the `\brief` entry and re-use this if available, otherwise the `\brief ...` text is added as
  detailed above. The script can also handle `\date`, `\par`, `\author`, `\note` and `\version`
  entries if these should be required in future. At present these entries are copied across
  unchanged if they exist and otherwise they are not added in.

- Finally anything else that was found in the comment block is output and suitably annotated. For
  example unmatched procedure arguments are flagged up, lines which begin with `!>`, `!`, or
  `!> \something_or_other` are also annotated.

- The annotations used for marking up the comment block are as follows:

  - `...` Added to `\brief` and `\param` only if no existing information is available. Other
    commonly occurring entries such as `\par`, `\author`, `\version`, `\note` and `\date` if they
    exist are passed through unchanged.

  - `UNMATCHED_PROCEDURE_ARGUMENT` - gets appended on to any procedure arguments in the comment
    block that don't match those in the procedure definition.

  - `UNKNOWN_DOXYGEN_COMMENT` - for parts of the comment block which look like a standard Doxygen
    block (e.g. lines that look like: `!> \something or other`) but not part of the standard header
    above.

  - `UNKNOWN_COMMENT` - lines that begin `!>` or `!` inside the comment block. These get placed into
    a `!> \note with UNKNOWN_COMMENT` added at the start of the line.

The script also includes checks to ensure that the annotations are not added to a line multiple
times.

## `remove_extra_comments.pl`

Post-processing script to remove any double lines which begin `! ******` e.g.

```fortran
! *****************************************************************************
! *****************************************************************************
```

gets changed to

```fortran
! *****************************************************************************
```

This is necessary because occasionally a procedure may be missing the start or end `! ****` line or
may have a duplicate one.
