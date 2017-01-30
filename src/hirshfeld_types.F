!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief The types needed for the calculation of Hirshfeld charges and
!>        related functions
!> \par History
!>      11.2014 created [JGH]
!> \author JGH
! **************************************************************************************************
MODULE hirshfeld_types

   USE input_constants,                 ONLY: radius_default,&
                                              shape_function_gaussian
   USE kinds,                           ONLY: dp
   USE pw_types,                        ONLY: pw_p_type,&
                                              pw_release
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hirshfeld_types'

   PUBLIC :: hirshfeld_type
   PUBLIC :: create_hirshfeld_type, release_hirshfeld_type
   PUBLIC :: get_hirshfeld_info, set_hirshfeld_info

! **************************************************************************************************
!> \brief quantities needed for a Hirshfeld based partitioning of real space
!> \author JGH
! **************************************************************************************************
   TYPE hirshfeld_type
      LOGICAL                       :: iterative, &
                                       use_bohr
      INTEGER                       :: shape_function_type
      INTEGER                       :: ref_charge, &
                                       radius_type
      TYPE(shape_fn), DIMENSION(:), &
         POINTER                    :: kind_shape_fn
      REAL(KIND=dp), DIMENSION(:), &
         POINTER                    :: charges
      TYPE(pw_p_type), POINTER      :: fnorm
   END TYPE hirshfeld_type

   TYPE shape_fn
      INTEGER                       :: numexp
      REAL(KIND=dp), DIMENSION(:), &
         POINTER                    :: zet
      REAL(KIND=dp), DIMENSION(:), &
         POINTER                    :: coef
   END TYPE shape_fn

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param hirshfeld_env ...
! **************************************************************************************************
   SUBROUTINE create_hirshfeld_type(hirshfeld_env)
      TYPE(hirshfeld_type), POINTER                      :: hirshfeld_env

      CHARACTER(len=*), PARAMETER :: routineN = 'create_hirshfeld_type', &
         routineP = moduleN//':'//routineN

      IF (ASSOCIATED(hirshfeld_env)) THEN
         CALL release_hirshfeld_type(hirshfeld_env)
      END IF

      ALLOCATE (hirshfeld_env)

      hirshfeld_env%iterative = .FALSE.
      hirshfeld_env%use_bohr = .FALSE.
      hirshfeld_env%shape_function_type = shape_function_gaussian
      hirshfeld_env%radius_type = radius_default
      NULLIFY (hirshfeld_env%kind_shape_fn)
      NULLIFY (hirshfeld_env%charges)
      NULLIFY (hirshfeld_env%fnorm)

   END SUBROUTINE create_hirshfeld_type

! **************************************************************************************************
!> \brief ...
!> \param hirshfeld_env ...
! **************************************************************************************************
   SUBROUTINE release_hirshfeld_type(hirshfeld_env)
      TYPE(hirshfeld_type), POINTER                      :: hirshfeld_env

      CHARACTER(len=*), PARAMETER :: routineN = 'release_hirshfeld_type', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ikind
      TYPE(shape_fn), DIMENSION(:), POINTER              :: kind_shape

      IF (ASSOCIATED(hirshfeld_env)) THEN

         IF (ASSOCIATED(hirshfeld_env%kind_shape_fn)) THEN
            kind_shape => hirshfeld_env%kind_shape_fn
            DO ikind = 1, SIZE(kind_shape)
               IF (ASSOCIATED(hirshfeld_env%kind_shape_fn(ikind)%zet)) THEN
                  DEALLOCATE (kind_shape(ikind)%zet)
               END IF
               IF (ASSOCIATED(hirshfeld_env%kind_shape_fn(ikind)%coef)) THEN
                  DEALLOCATE (kind_shape(ikind)%coef)
               END IF
            END DO
            DEALLOCATE (kind_shape)
         END IF

         IF (ASSOCIATED(hirshfeld_env%charges)) THEN
            DEALLOCATE (hirshfeld_env%charges)
         END IF

         IF (ASSOCIATED(hirshfeld_env%fnorm)) THEN
            CALL pw_release(hirshfeld_env%fnorm%pw)
            DEALLOCATE (hirshfeld_env%fnorm)
         ENDIF

         DEALLOCATE (hirshfeld_env)

      END IF

   END SUBROUTINE release_hirshfeld_type

! **************************************************************************************************
!> \brief Get information from a Hirshfeld env
!> \param hirshfeld_env the env that holds the information
!> \param shape_function_type the type of shape function used
!> \param iterative logical which determins if iterative Hirshfeld charges should be computed
!> \param ref_charge the reference charge type (core charge or mulliken)
!> \param fnorm normalization of the shape function
!> \param radius_type the type of radius used for building the shape functions
!> \param use_bohr logical which determines if angstrom or bohr units are used to build the
!>                 shape functions
! **************************************************************************************************
   SUBROUTINE get_hirshfeld_info(hirshfeld_env, shape_function_type, iterative, &
                                 ref_charge, fnorm, radius_type, use_bohr)
      TYPE(hirshfeld_type), POINTER                      :: hirshfeld_env
      INTEGER, INTENT(OUT), OPTIONAL                     :: shape_function_type
      LOGICAL, INTENT(OUT), OPTIONAL                     :: iterative
      INTEGER, INTENT(OUT), OPTIONAL                     :: ref_charge
      TYPE(pw_p_type), OPTIONAL, POINTER                 :: fnorm
      INTEGER, INTENT(OUT), OPTIONAL                     :: radius_type
      LOGICAL, INTENT(OUT), OPTIONAL                     :: use_bohr

      CHARACTER(len=*), PARAMETER :: routineN = 'get_hirshfeld_info', &
         routineP = moduleN//':'//routineN

      CPASSERT(ASSOCIATED(hirshfeld_env))

      IF (PRESENT(shape_function_type)) THEN
         shape_function_type = hirshfeld_env%shape_function_type
      END IF
      IF (PRESENT(iterative)) THEN
         iterative = hirshfeld_env%iterative
      END IF
      IF (PRESENT(use_bohr)) THEN
         use_bohr = hirshfeld_env%use_bohr
      END IF
      IF (PRESENT(radius_type)) THEN
         radius_type = hirshfeld_env%radius_type
      END IF
      IF (PRESENT(ref_charge)) THEN
         ref_charge = hirshfeld_env%ref_charge
      END IF
      IF (PRESENT(fnorm)) THEN
         fnorm => hirshfeld_env%fnorm
      END IF

   END SUBROUTINE get_hirshfeld_info

! **************************************************************************************************
!> \brief Set values of a Hirshfeld env
!> \param hirshfeld_env the env that holds the information
!> \param shape_function_type the type of shape function used
!> \param iterative logical which determins if iterative Hirshfeld charges should be computed
!> \param ref_charge the reference charge type (core charge or mulliken)
!> \param fnorm normalization of the shape function
!> \param radius_type the type of radius used for building the shape functions
!> \param use_bohr logical which determines if angstrom or bohr units are used to build the
!>                 shape functions
! **************************************************************************************************
   SUBROUTINE set_hirshfeld_info(hirshfeld_env, shape_function_type, iterative, &
                                 ref_charge, fnorm, radius_type, use_bohr)
      TYPE(hirshfeld_type), POINTER                      :: hirshfeld_env
      INTEGER, INTENT(IN), OPTIONAL                      :: shape_function_type
      LOGICAL, INTENT(IN), OPTIONAL                      :: iterative
      INTEGER, INTENT(IN), OPTIONAL                      :: ref_charge
      TYPE(pw_p_type), OPTIONAL, POINTER                 :: fnorm
      INTEGER, INTENT(IN), OPTIONAL                      :: radius_type
      LOGICAL, INTENT(IN), OPTIONAL                      :: use_bohr

      CHARACTER(len=*), PARAMETER :: routineN = 'set_hirshfeld_info', &
         routineP = moduleN//':'//routineN

      CPASSERT(ASSOCIATED(hirshfeld_env))

      IF (PRESENT(shape_function_type)) THEN
         hirshfeld_env%shape_function_type = shape_function_type
      END IF
      IF (PRESENT(iterative)) THEN
         hirshfeld_env%iterative = iterative
      END IF
      IF (PRESENT(use_bohr)) THEN
         hirshfeld_env%use_bohr = use_bohr
      END IF
      IF (PRESENT(radius_type)) THEN
         hirshfeld_env%radius_type = radius_type
      END IF
      IF (PRESENT(ref_charge)) THEN
         hirshfeld_env%ref_charge = ref_charge
      END IF
      IF (PRESENT(fnorm)) THEN
         hirshfeld_env%fnorm => fnorm
      END IF

   END SUBROUTINE set_hirshfeld_info
! **************************************************************************************************

END MODULE hirshfeld_types
