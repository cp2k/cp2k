# Error Handling

In CP2K there exist three convenient macros for error handling. They are defined in
`src/base/base_uses.f90`, and are therefore always available. However, these macros can only be used
for **short messages** (approximately less than 7 words), because of Fortran's line length
limitation.

```
   CPWARN("A short warning message")
   CPABORT("A short error message")
   CPASSERT(2 > 1) ! For sanity checks only, not for user errors.
```

For longer messages one can use these two routines:

```
   CALL cp_warn(__LOCATION__,"A long warning message "//&
                "which can span multiple lines and contains "//&
                "some valuable advice for the user.")

   CALL cp_abort(__LOCATION__,"A long error message "//&
                 "with even more valuable advice.")      
```

```{note}
Messages from `CPWARN` and `cp_warn` are only printed on the first MPI-rank.
Warnings issued on other ranks are ignored.
```

## Cheat Sheet

If you want to...

- print a **short warning** message, then you should use: `CPWARN("your short message")`
- print a **longer warning** message, then you should use:

```
  CALL cp_warn(__LOCATION__,"your lengthy message, "//&
              "which can even span multiple lines")
```

- **stop** the program with a **short** error message, then you should use:
  `CPABORT("your short message")`
- **stop** the program with a **longer** error message, then you should use:

```
   CALL cp_abort(__LOCATION__,"your lengthy message, "//&
               "which can even span multiple lines")
```

- have a **simple assertion** without a custom error message, then you should use:
  `CPASSERT(your_logical_expression)`
- have a **simple assertion** with a **short** error message, then you should use:

```
IF(.NOT.your_logical_expression) CPABORT("your short message")
```

- have a **complex assertion** with a **short** error message, then you should use:

```
  IF(.NOT.your_super_complicated_logical_expression_that_takes_alot_of_space)&
    CPABORT("your short message")
```

- have a **lengthy assertion** with a **longer** error message, then you should use:

```
   IF(.NOT.your_super_long_logical_expression_that_takes_alot_of_space)&
      CALL cp_abort(__LOCATION__,"your lengthy message"//&
                    "which again can even span multiple lines")
```
