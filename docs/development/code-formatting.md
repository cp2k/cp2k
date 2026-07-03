# Code Formatting

Uniform formatting of CP2K source code is enabled by a
[fprettify](https://github.com/fortran-lang/fprettify) that is an almost complete auto-formatter for
Fortran 90 source code. As a rule of thumb, developers should not worry about the format of their
code and just let prettify do its magic by running `./make_pretty.sh`.

## Formatting conventions

The following formatting conventions are automatically enforced by the `./make_pretty.sh` command:

- Sorting and alignment of variable declarations and `USE` statements, removal of unused list
  entries.
- Indentation with a relative width of 3 characters.
- Line continuations are aligned with the previous opening delimiter `(`, `[` or `(/` or with an
  assignment operator `=` or `=>`. If none of the above is present, a default hanging indent of 3
  characters is applied.
- All operators are surrounded by exactly one whitespace character, except for arithmetic operators.
- Removal of extraneous whitespace and consecutive blank lines.
- Uppercase notation for all Fortran and OpenMP keywords.

## Manual formatting

The following formatting decisions are still manual and are never changed by prettify:

- Positions of line breaks (except for variable declarations and `USE` statements).
- No indentation of subsequent `DO` / `IF` statements that are aligned with each other. There may be
  cases where manual alignment is preferred over the automatic formatting conventions. The following
  options for manual formatting are provided:
- No automatic realignment of line continuations that are prefixed with an `&`.
- No auto-formatting of lines to which a comment starting with `!&` is attached.
- No auto-formatting of code blocks enclosed between two comment lines starting with `!&<` and
  `!&>`.

## Examples

A few examples to illustrate how to deal with cases where auto-formatting produces unsatisfying
results:

- Reduce hanging indent by inserting linebreaks directly after assignment operator and opening
  delimiter:

```
! No:
long_result_var_name = long_function_name(arg_1, arg_2, &
                                          arg_3, arg_4, arg_5)+ &
                       foo
! Yes:
long_result_var_name = &
   long_function_name( &
      arg_1, arg_2, &
      arg_3, arg_4, arg_5)+ &
   foo
```

- Avoid linebreaks in deeply nested expressions:

```
! No:
bessj0 = (r1+y*(r2+y*(r3+y*(r4+y* &
                            (r5+y*r6)))))/(s1+y*(s2+y*(s3+y* &
                                                       (s4+y*(s5+y*s6)))))
! Yes:
bessj0 = (r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/ &
         (s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
```

- Alignment by explicit bracketing:

```
! No:
foo = bar+foobar(x1, y1, z1)* &
      foobar(x2, y2, z2)* &
      foobar(x3, y3, z3)
! Yes:
foo = bar+(foobar(x1, y1, z1)* &
           foobar(x2, y2, z2)* &
           foobar(x3, y3, z3))
```

- Special vertical alignment may require manual formatting:

```
! Auto-formatting:
align_me = [-1, 10, 0, &
            0, 1000, 0, &
            0, -1, 1]
! Manual alignment (!& disables whitespace formatting):
align_me = [-1,   10, 0, & !&
             0, 1000, 0, & !&
             0,   -1, 1] !&
! Alternatively:
!&<
align_me = [-1,   10, 0, &
             0, 1000, 0, &
             0,   -1, 1]
!&>
```
