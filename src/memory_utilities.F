!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!****** cp2k/memory_utilities [1.0] *
!!
!!   NAME
!!     memory_utilities
!!
!!   FUNCTION
!!     Utility routines for the memory handling.
!!
!!   AUTHOR
!!     Matthias Krack (25.06.1999)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
MODULE memory_utilities

  USE f77_blas
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             dp_size,&
                                             int_8,&
                                             int_size
  USE termination,                     ONLY: stop_memory
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'memory_utilities'

  PUBLIC :: reallocate

  INTERFACE reallocate
    MODULE PROCEDURE reallocate_c1,reallocate_c2,reallocate_c3,reallocate_c4,&
                     reallocate_i1,reallocate_i2,reallocate_i3,reallocate_i4,&
                     reallocate_r1,reallocate_r2,reallocate_r3,reallocate_r4,&
                     reallocate_r5,reallocate_s1,reallocate_l1,reallocate_8i1,&
                     reallocate_8i2
  END INTERFACE

!!***
! *****************************************************************************

CONTAINS

!!****** memory_utilities/reallocate_c1 [1.0] *
!!
!!   NAME
!!     reallocate_c1
!!
!!   FUNCTION
!!     (Re)Allocate a complex vector with a new dimension 
!!
!!   AUTHOR
!!     Matthias Krack (28.11.2005,MK)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE reallocate_c1(p,lb1_new,ub1_new)
    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new

    CHARACTER(len=*), PARAMETER :: routineN = 'reallocate_c1', &
         routineP = moduleN//':'//routineN
    
    COMPLEX(KIND=dp), PARAMETER              :: zero=(0.0_dp,0.0_dp)
    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:)                           :: work
#include "reallocate_1.f90"

  END SUBROUTINE reallocate_c1

!!****** memory_utilities/reallocate_c2 [1.0] *
!!
!!   NAME
!!     reallocate_c2
!!
!!   FUNCTION
!!     (Re)Allocate a complex array 2D with a new dimension 
!!
!!   AUTHOR
!!     Matthias Krack (28.11.2005,MK)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE reallocate_c2(p,lb1_new,ub1_new,lb2_new,ub2_new)
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new, lb2_new, &
                                                ub2_new

    CHARACTER(len=*), PARAMETER :: routineN = 'reallocate_c2', &
         routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), PARAMETER              :: zero=(0.0_dp,0.0_dp)
    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: work
#include "reallocate_2.f90"

  END SUBROUTINE reallocate_c2

!!****** memory_utilities/reallocate_c3 [1.0] *
!!
!!   NAME
!!     reallocate_c3
!!
!!   FUNCTION
!!     (Re)Allocate a complex array 3D with a new dimension 
!!
!!   AUTHOR
!!     Matthias Krack (28.11.2005,MK)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE reallocate_c3(p,lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new)
    COMPLEX(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new, lb2_new, &
                                                ub2_new, lb3_new, ub3_new

    CHARACTER(len=*), PARAMETER :: routineN = 'reallocate_c3', &
         routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), PARAMETER              :: zero=(0.0_dp,0.0_dp)
    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: work
#include "reallocate_3.f90"

  END SUBROUTINE reallocate_c3

!!****** memory_utilities/reallocate_c4 [1.0] *
!!
!!   NAME
!!     reallocate_c4
!!
!!   FUNCTION
!!     (Re)Allocate a complex array 4D with a new dimension 
!!
!!   AUTHOR
!!     Matthias Krack (28.11.2005,MK)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE reallocate_c4(p,lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new,&
       lb4_new,ub4_new)
    COMPLEX(KIND=dp), &
      DIMENSION(:, :, :, :), POINTER         :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new, lb2_new, &
                                                ub2_new, lb3_new, ub3_new, &
                                                lb4_new, ub4_new

    CHARACTER(len=*), PARAMETER :: routineN = 'reallocate_c4', &
         routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), PARAMETER              :: zero=(0.0_dp,0.0_dp)
    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: work
#include "reallocate_4.f90"

  END SUBROUTINE reallocate_c4

!!****** memory_utilities/reallocate_i1 [1.0] *
!!
!!   NAME
!!     reallocate_i1
!!
!!   FUNCTION
!!     (Re)Allocate an integer vector with a new dimension.
!!
!!   AUTHOR
!!     Matthias Krack (18.07.2002,MK)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE reallocate_i1(p,lb1_new,ub1_new)
    INTEGER, DIMENSION(:), POINTER           :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new

    CHARACTER(len=*), PARAMETER :: routineN = 'reallocate_i1', &
          routineP = moduleN//':'//routineN

    INTEGER, PARAMETER                       :: zero=0
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: work
#include "reallocate_1.f90"

  END SUBROUTINE reallocate_i1

!!****** memory_utilities/reallocate_8i1 [1.0] *
!!
!!   NAME
!!     reallocate_8i1
!!
!!   FUNCTION
!!     (Re)Allocate an integer (int_8) vector with a new dimension.
!!
!!   AUTHOR
!!     Matthias Krack (18.07.2002,MK)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE reallocate_8i1(p,lb1_new,ub1_new)
    INTEGER(kind=int_8), DIMENSION(:), &
      POINTER                                :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new

    CHARACTER(len=*), PARAMETER :: routineN = 'reallocate_8i1', &
          routineP = moduleN//':'//routineN

    INTEGER(kind=int_8), PARAMETER           :: zero=0
    INTEGER(kind=int_8), ALLOCATABLE, &
      DIMENSION(:)                           :: work
#include "reallocate_1.f90"

  END SUBROUTINE reallocate_8i1

!!****** memory_utilities/reallocate_8i2 [1.0] *
!!
!!   NAME
!!     reallocate_8i2
!!
!!   FUNCTION
!!     (Re)Allocate an integer (int_8) array with a new dimension.
!!
!!   AUTHOR
!!     Matthias Krack (18.07.2002,MK)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE reallocate_8i2(p,lb1_new,ub1_new,lb2_new,ub2_new)
    INTEGER(kind=int_8), DIMENSION(:, :), &
      POINTER                                :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new, lb2_new, &
                                                ub2_new

    CHARACTER(len=*), PARAMETER :: routineN = 'reallocate_8i2', &
         routineP = moduleN//':'//routineN

    INTEGER(kind=int_8), PARAMETER           :: zero=0
    INTEGER(kind=int_8), ALLOCATABLE, &
      DIMENSION(:, :)                        :: work
#include "reallocate_2.f90"

  END SUBROUTINE reallocate_8i2

!!****** memory_utilities/reallocate_i2 [1.0] *
!!
!!   NAME
!!     reallocate_i2
!!
!!   FUNCTION
!!     (Re)Allocate an integer array with a new dimension.
!!
!!   AUTHOR
!!     Matthias Krack (18.07.2002,MK)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE reallocate_i2(p,lb1_new,ub1_new,lb2_new,ub2_new)
    INTEGER, DIMENSION(:, :), POINTER        :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new, lb2_new, &
                                                ub2_new

    CHARACTER(len=*), PARAMETER :: routineN = 'reallocate_i2', &
         routineP = moduleN//':'//routineN

    INTEGER, PARAMETER                       :: zero=0
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: work
#include "reallocate_2.f90"

  END SUBROUTINE reallocate_i2

!!****** memory_utilities/reallocate_i3 [1.0] *
!!
!!   NAME
!!     reallocate_i3
!!
!!   FUNCTION
!!     (Re)Allocate an integer array 3D with a new dimension.
!!
!!   AUTHOR
!!     Matthias Krack (18.07.2002,MK)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE reallocate_i3(p,lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new)
    INTEGER, DIMENSION(:, :, :), POINTER     :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new, lb2_new, &
                                                ub2_new, lb3_new, ub3_new

    CHARACTER(len=*), PARAMETER :: routineN = 'reallocate_i3', &
         routineP = moduleN//':'//routineN

    INTEGER, PARAMETER                       :: zero=0
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: work
#include "reallocate_3.f90"

  END SUBROUTINE reallocate_i3

!!****** memory_utilities/reallocate_i3 [1.0] *
!!
!!   NAME
!!     reallocate_i3
!!
!!   FUNCTION
!!     (Re)Allocate an integer array 3D with a new dimension.
!!
!!   AUTHOR
!!     Matthias Krack (04.10.2002,MK)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE reallocate_i4(p,lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new,&
       lb4_new,ub4_new)
    INTEGER, DIMENSION(:, :, :, :), POINTER  :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new, lb2_new, &
                                                ub2_new, lb3_new, ub3_new, &
                                                lb4_new, ub4_new

    CHARACTER(len=*), PARAMETER :: routineN = 'reallocate_i4', &
         routineP = moduleN//':'//routineN

    INTEGER, PARAMETER                       :: zero=0
    INTEGER, ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: work
#include "reallocate_4.f90"

  END SUBROUTINE reallocate_i4

!!****** memory_utilities/reallocate_r1 [1.0] *
!!
!!   NAME
!!     reallocate_r1
!!
!!   FUNCTION
!!     (Re)Allocate an real vector with a new dimension.
!!
!!   AUTHOR
!!     Matthias Krack (18.07.2002,MK)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE reallocate_r1(p,lb1_new,ub1_new)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new

    CHARACTER(len=*), PARAMETER :: routineN = 'reallocate_r1', &
         routineP = moduleN//':'//routineN

    REAL(KIND=dp), PARAMETER                 :: zero=0.0_dp
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: work
#include "reallocate_1.f90"

  END SUBROUTINE reallocate_r1

!!****** memory_utilities/reallocate_r2 [1.0] *
!!
!!   NAME
!!     reallocate_r2
!!
!!   FUNCTION
!!     (Re)Allocate an real array with new dimensions.
!!
!!   AUTHOR
!!     Matthias Krack (18.07.2002,MK)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE reallocate_r2(p,lb1_new,ub1_new,lb2_new,ub2_new)
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new, lb2_new, &
                                                ub2_new

    CHARACTER(len=*), PARAMETER :: routineN = 'reallocate_r2', &
         routineP = moduleN//':'//routineN

    REAL(KIND=dp), PARAMETER                 :: zero=0.0_dp
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: work
#include "reallocate_2.f90"

  END SUBROUTINE reallocate_r2

!!****** memory_utilities/reallocate_r3 [1.0] *
!!
!!   NAME
!!     reallocate_r3
!!
!!   FUNCTION
!!     (Re)Allocate an real array 3D with new dimensions.
!!
!!   AUTHOR
!!     Matthias Krack (18.07.2002,MK)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE reallocate_r3(p,lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new)
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new, lb2_new, &
                                                ub2_new, lb3_new, ub3_new

    CHARACTER(len=*), PARAMETER :: routineN = 'reallocate_r3', &
         routineP = moduleN//':'//routineN

    REAL(KIND=dp), PARAMETER                 :: zero=0.0_dp
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: work
#include "reallocate_3.f90"

  END SUBROUTINE reallocate_r3

!!****** memory_utilities/reallocate_r4 [1.0] *
!!
!!   NAME
!!     reallocate_r4
!!
!!   FUNCTION
!!     (Re)Allocate an real array 4D with new dimensions.
!!
!!   AUTHOR
!!     Matthias Krack (04.10.2002,MK)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE reallocate_r4(p,lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new,&
                             lb4_new,ub4_new)
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new, lb2_new, &
                                                ub2_new, lb3_new, ub3_new, &
                                                lb4_new, ub4_new

    CHARACTER(len=*), PARAMETER :: routineN = 'reallocate_r4', &
         routineP = moduleN//':'//routineN

    REAL(KIND=dp), PARAMETER                 :: zero=0.0_dp
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: work
#include "reallocate_4.f90"

  END SUBROUTINE reallocate_r4

!!****** memory_utilities/reallocate_r5 [1.0] *
!!
!!   NAME
!!     reallocate_r5
!!
!!   FUNCTION
!!     (Re)Allocate an real array 5D with new dimensions.
!!
!!   AUTHOR
!!     Matthias Krack (04.10.2002,MK)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE reallocate_r5(p,lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new,&
                             lb4_new,ub4_new,lb5_new,ub5_new)
    REAL(KIND=dp), &
      DIMENSION(:, :, :, :, :), POINTER      :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new, lb2_new, &
                                                ub2_new, lb3_new, ub3_new, &
                                                lb4_new, ub4_new, lb5_new, &
                                                ub5_new

    CHARACTER(len=*), PARAMETER :: routineN = 'reallocate_r5', &
         routineP = moduleN//':'//routineN

    REAL(KIND=dp), PARAMETER                 :: zero=0.0_dp
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :, :)               :: work
#include "reallocate_5.f90"

  END SUBROUTINE reallocate_r5

!!****** memory_utilities/reallocate_l1 [1.0] *
!!
!!   NAME
!!     reallocate_l1
!!
!!   FUNCTION
!!     (Re)Allocate a logical vector with a new dimension.
!!
!!   AUTHOR
!!     Matthias Krack (18.07.2002,MK)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE reallocate_l1(p,lb1_new,ub1_new)
    LOGICAL, DIMENSION(:), POINTER           :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new

    CHARACTER(len=*), PARAMETER :: routineN = 'reallocate_l1', &
         routineP = moduleN//':'//routineN

    LOGICAL, PARAMETER                       :: zero=.FALSE.
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: work
#include "reallocate_1.f90"

  END SUBROUTINE reallocate_l1

!!****** memory_utilities/reallocate_s1 [1.0] *
!!
!!   NAME
!!     reallocate_s1
!!
!!   FUNCTION
!!     (Re)Allocate a vector of string variables of default string length.
!!
!!   AUTHOR
!!     Matthias Krack (18.07.2002,MK)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   NOTES
!!     (Maybe outdated) This routine doesnt work on SUN/Solaris!!! 
!!                      It should probably not be used. (Thomas Chassaing)
!!
!******************************************************************************
  SUBROUTINE reallocate_s1(p,lb1_new,ub1_new)
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new

    CHARACTER(len=*), PARAMETER :: routineN = 'reallocate_s1', &
         routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length),&
         PARAMETER                           :: zero=""
    CHARACTER(LEN=default_string_length), &
      ALLOCATABLE, DIMENSION(:)              :: work
#include "reallocate_1.f90"

  END SUBROUTINE reallocate_s1

END MODULE memory_utilities
