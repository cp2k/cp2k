
MODULE kpoints

  USE kinds,                           ONLY: dp

  IMPLICIT NONE

  TYPE kpts_type
     REAL (KIND=dp), DIMENSION ( :, : ), POINTER :: xkpt
     REAL (KIND=dp), DIMENSION ( : ), POINTER :: wkpt
  END TYPE kpts_type

  TYPE ( kpts_type ) :: kset

END MODULE kpoints
