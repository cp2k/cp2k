!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Wrapper to allocate and optionally assign array from source array.
!>        Will be redundant once compilers have implemented Fortran 2008 standard :-/
!> \author Patrick Seewald
! **************************************************************************************************

MODULE allocate_wrap_c
   USE kinds,                           ONLY: real_4

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: allocate_any

   INTERFACE allocate_any
      MODULE PROCEDURE allocate_1d_c
      MODULE PROCEDURE allocate_2d_c
      MODULE PROCEDURE allocate_3d_c
      MODULE PROCEDURE allocate_4d_c
      MODULE PROCEDURE allocate_5d_c
      MODULE PROCEDURE allocate_6d_c
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param sizes ...
!> \param array ...
!> \param source ...
! **************************************************************************************************
   PURE SUBROUTINE allocate_1d_c(sizes, array, source)
      INTEGER, DIMENSION(1), INTENT(IN)                  :: sizes
      COMPLEX(kind=real_4), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: array
      COMPLEX(kind=real_4), DIMENSION(:), INTENT(IN), &
         OPTIONAL                                        :: source

      ALLOCATE (array(sizes(1)))
      IF (PRESENT(source)) THEN
         array(:) = source(:)
      ENDIF
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param sizes ...
!> \param array ...
!> \param source ...
! **************************************************************************************************
   PURE SUBROUTINE allocate_2d_c(sizes, array, source)
      INTEGER, DIMENSION(2), INTENT(IN)                  :: sizes
      COMPLEX(kind=real_4), ALLOCATABLE, &
         DIMENSION(:, :), INTENT(OUT)                    :: array
      COMPLEX(kind=real_4), DIMENSION(:, :), &
         INTENT(IN), OPTIONAL                            :: source

      ALLOCATE (array(sizes(1), sizes(2)))
      IF (PRESENT(source)) THEN
         array(:, :) = source(:, :)
      ENDIF
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param sizes ...
!> \param array ...
!> \param source ...
! **************************************************************************************************
   PURE SUBROUTINE allocate_3d_c(sizes, array, source)
      INTEGER, DIMENSION(3), INTENT(IN)                  :: sizes
      COMPLEX(kind=real_4), ALLOCATABLE, &
         DIMENSION(:, :, :), INTENT(OUT)                 :: array
      COMPLEX(kind=real_4), DIMENSION(:, :, :), &
         INTENT(IN), OPTIONAL                            :: source

      ALLOCATE (array(sizes(1), sizes(2), sizes(3)))
      IF (PRESENT(source)) THEN
         array(:, :, :) = source(:, :, :)
      ENDIF
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param sizes ...
!> \param array ...
!> \param source ...
! **************************************************************************************************
   PURE SUBROUTINE allocate_4d_c(sizes, array, source)
      INTEGER, DIMENSION(4), INTENT(IN)                  :: sizes
      COMPLEX(kind=real_4), ALLOCATABLE, &
         DIMENSION(:, :, :, :), INTENT(OUT)              :: array
      COMPLEX(kind=real_4), DIMENSION(:, :, :, :), &
         INTENT(IN), OPTIONAL                            :: source

      ALLOCATE (array(sizes(1), sizes(2), sizes(3), sizes(4)))
      IF (PRESENT(source)) THEN
         array(:, :, :, :) = source(:, :, :, :)
      ENDIF
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param sizes ...
!> \param array ...
!> \param source ...
! **************************************************************************************************
   PURE SUBROUTINE allocate_5d_c(sizes, array, source)
      INTEGER, DIMENSION(5), INTENT(IN)                  :: sizes
      COMPLEX(kind=real_4), ALLOCATABLE, &
         DIMENSION(:, :, :, :, :), INTENT(OUT)           :: array
      COMPLEX(kind=real_4), DIMENSION(:, :, :, :, :), &
         INTENT(IN), OPTIONAL                            :: source

      ALLOCATE (array(sizes(1), sizes(2), sizes(3), sizes(4), sizes(5)))
      IF (PRESENT(source)) THEN
         array(:, :, :, :, :) = source(:, :, :, :, :)
      ENDIF
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param sizes ...
!> \param array ...
!> \param source ...
! **************************************************************************************************
   PURE SUBROUTINE allocate_6d_c(sizes, array, source)
      INTEGER, DIMENSION(6), INTENT(IN)                  :: sizes
      COMPLEX(kind=real_4), ALLOCATABLE, &
         DIMENSION(:, :, :, :, :, :), INTENT(OUT)        :: array
      COMPLEX(kind=real_4), &
         DIMENSION(:, :, :, :, :, :), INTENT(IN), &
         OPTIONAL                                        :: source

      ALLOCATE (array(sizes(1), sizes(2), sizes(3), sizes(4), sizes(5), sizes(6)))
      IF (PRESENT(source)) THEN
         array(:, :, :, :, :, :) = source(:, :, :, :, :, :)
      ENDIF
   END SUBROUTINE

END MODULE
