!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \par History
!> \author JGH
! **************************************************************************************************
MODULE fist_efield_types
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'fist_efield_types'

! **************************************************************************************************
   TYPE fist_efield_type
      LOGICAL                                :: apply_field = .FALSE.
      LOGICAL                                :: displacement = .FALSE.
      REAL(KIND=dp)                          :: strength
      REAL(KIND=dp), DIMENSION(3)            :: polarisation
      REAL(KIND=dp), DIMENSION(3)            :: dfilter
   END TYPE fist_efield_type
! **************************************************************************************************

   PRIVATE

   PUBLIC :: fist_efield_type
   PUBLIC :: read_efield_section

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Read input section PERIODIC_EFIELD
!> \param input_section ...
!> \param efield ...
!> \par History
!> \author JGH
! **************************************************************************************************
   SUBROUTINE read_efield_section(input_section, efield)
      TYPE(section_vals_type), POINTER                   :: input_section
      TYPE(fist_efield_type), POINTER                    :: efield

      CHARACTER(len=*), PARAMETER :: routineN = 'read_efield_section', &
         routineP = moduleN//':'//routineN

      REAL(KIND=dp), DIMENSION(:), POINTER               :: pp
      TYPE(section_vals_type), POINTER                   :: tmp_section

      IF (.NOT. ASSOCIATED(efield)) ALLOCATE (efield)

      ! Read the finite field input section for periodic fields
      tmp_section => section_vals_get_subs_vals(input_section, "PERIODIC_EFIELD")
      CALL section_vals_get(tmp_section, explicit=efield%apply_field)
      IF (efield%apply_field) THEN
         CALL section_vals_val_get(tmp_section, "POLARISATION", r_vals=pp)
         efield%polarisation(1:3) = pp(1:3)
         CALL section_vals_val_get(tmp_section, "D_FILTER", r_vals=pp)
         efield%dfilter(1:3) = pp(1:3)
         CALL section_vals_val_get(tmp_section, "INTENSITY", r_val=efield%strength)
         CALL section_vals_val_get(tmp_section, "DISPLACEMENT_FIELD", l_val=efield%displacement)
      END IF

   END SUBROUTINE read_efield_section

! **************************************************************************************************

END MODULE fist_efield_types
