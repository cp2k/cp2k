
MODULE kpoints

  USE kinds,                           ONLY: dbl

  IMPLICIT NONE

  TYPE kpts_type
     REAL ( dbl ), DIMENSION ( :, : ), POINTER :: xkpt
     REAL ( dbl ), DIMENSION ( : ), POINTER :: wkpt
  END TYPE kpts_type

  TYPE ( kpts_type ) :: kset

END MODULE kpoints
