# Quickstep - Orbital Transformation Method (Regression Test)

**Authors:** V. Weber and U. Borstnik

The purpose of this benchmark is to analyse the time spent in the following
linear algebra routines as displayed at the end of the cp2k output.

**Main Routines:**

- `make_preconditioner`
- `cp_dbcsr_multiply_d`
- `apply_preconditioner`

**Secondary Routines** (separate timings for different multiplication options):

- `dbcsr_mult_NSS_NRN`
- `dbcsr_mult_NRN_TRN`
- `dbcsr_mult_NSN_NRN`
- `dbcsr_mult_TRN_NRN`
- `dbcsr_mult_NRN_NSN`

## Results Archive

On rosa the total run times are approximately:

| Input file    | Number cores | Runtime               |
| ------------- | ------------:| ---------------------:|
| H2O-256.inp   | 576 cores    | 5 min                 |
| H2O-1024.inp  | 576 cores    | 40 min                |
| H2O-4096.inp  | 2304 cores   | 15 min                |
| H2O-65536.inp |              | doesn't run so far... |

Runs performed with cp2k from 30.01.09 (no libdbcsr routines used), compiled with gfortran and linked to the default libraries.

