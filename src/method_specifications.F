MODULE method_specifications

! Purpose: Calculation method specification.

! History: - Creation (13.06.2000, Matthias Krack)

! *****************************************************************************

  USE functionals,                     ONLY: xc_type

  IMPLICIT NONE

  PRIVATE

  TYPE(xc_type)    :: xc_info
  CHARACTER(LEN=8) :: method = "GAPW"
  INTEGER          :: maxder = 0
  LOGICAL          :: allchem = .FALSE.,&
                      gapw = .TRUE.,&
                      gpw = .TRUE.,&
                      oks = .FALSE.,&
                      periodic = .TRUE.

! *** Public variables ***

  PUBLIC :: allchem,&
            gapw,&
            gpw,&
            maxder,&
            method,&
            oks,&
            periodic,&
            xc_info

! *** Public subroutines ***

  PUBLIC :: get_maxder

! *****************************************************************************

CONTAINS

! *****************************************************************************
  SUBROUTINE get_maxder(maxder_index)
    INTEGER, INTENT(in), OPTIONAL            :: maxder_index

! maxder

    IF (PRESENT(maxder_index)) maxder = maxder_index
  END SUBROUTINE get_maxder
! *****************************************************************************

END MODULE method_specifications
