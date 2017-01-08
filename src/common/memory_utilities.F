!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utility routines for the memory handling.
!>        Uses stop_memory to handle exceptions
!> \par History
!>      none
!> \author Matthias Krack (25.06.1999)
! **************************************************************************************************
MODULE memory_utilities

   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp,&
                                              dp_size,&
                                              int_8,&
                                              int_8_size,&
                                              int_size
   USE termination,                     ONLY: stop_memory
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'memory_utilities'

   PUBLIC :: reallocate

   INTERFACE reallocate
      MODULE PROCEDURE reallocate_c1, reallocate_c2, reallocate_c3, reallocate_c4, &
         reallocate_i1, reallocate_i2, reallocate_i3, reallocate_i4, &
         reallocate_r1, reallocate_r2, reallocate_r3, reallocate_r4, &
         reallocate_r5, reallocate_s1, reallocate_l1, reallocate_8i1, &
         reallocate_8i2
   END INTERFACE

CONTAINS

#! *************************************************************************************************
#!> \brief Fypp macro for common subroutine body
#!> \author Ole Schuett
#! *************************************************************************************************
#:def reallocate(rank)

      INTEGER :: istat
#:for i in range(1, rank+1)
      INTEGER :: lb${i}$, lb${i}$_old, ub${i}$, ub${i}$_old
#:endfor

    IF (ASSOCIATED(p)) THEN

#:for i in range(1, rank+1)
       lb${i}$_old = LBOUND(p,${i}$)
       ub${i}$_old = UBOUND(p,${i}$)
       lb${i}$ = MAX(lb${i}$_new,lb${i}$_old)
       ub${i}$ = MIN(ub${i}$_new,ub${i}$_old)
#:endfor

#:set old_bounds = ','.join(['lb{0}:ub{0}'.format(i) for i in range(1, rank+1)])
#:set old_size  = '*'.join(['(ub{0}-lb{0}+1)'.format(i) for i in range(1, rank+1)])

       ALLOCATE (work( ${old_bounds}$ ), STAT=istat)
       IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,__LINE__,&
                           "work",t_size*${old_size}$)
       END IF
       work(${old_bounds}$) = p(${old_bounds}$)
       DEALLOCATE (p)
    END IF

#:set new_bounds = ','.join(['lb{0}_new:ub{0}_new'.format(i) for i in range(1, rank+1)])
#:set new_size  = '*'.join(['(ub{0}_new-lb{0}_new+1)'.format(i) for i in range(1, rank+1)])

    ALLOCATE (p(${new_bounds}$),STAT=istat)
    IF (istat /= 0) THEN
       CALL stop_memory(routineN,moduleN,__LINE__,&
                        "p",t_size*${new_size}$)
    END IF
    p = zero

    IF (ASSOCIATED(p).AND.ALLOCATED(work)) THEN
       p(${old_bounds}$) = work(${old_bounds}$)
       DEALLOCATE (work)
    END IF

#:enddef

! **************************************************************************************************
!> \brief (Re)Allocate a complex vector with a new dimension
!> \param p ...
!> \param lb1_new ...
!> \param ub1_new ...
!> \par History
!>      none
!> \author Matthias Krack (28.11.2005,MK)
! **************************************************************************************************
   SUBROUTINE reallocate_c1(p, lb1_new, ub1_new)

      COMPLEX(KIND=dp), DIMENSION(:), POINTER            :: p
      INTEGER, INTENT(IN)                                :: lb1_new, ub1_new

      CHARACTER(LEN=*), PARAMETER :: routineN = 'reallocate_c1', routineP = moduleN//':'//routineN
      COMPLEX(KIND=dp), PARAMETER                        :: zero = (0.0_dp, 0.0_dp)
      INTEGER, PARAMETER                                 :: t_size = 2*dp_size

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: work

$: reallocate(rank=1)

   END SUBROUTINE reallocate_c1

! **************************************************************************************************
!> \brief (Re)Allocate a complex array 2D with a new dimension
!> \param p ...
!> \param lb1_new ...
!> \param ub1_new ...
!> \param lb2_new ...
!> \param ub2_new ...
!> \par History
!>      none
!> \author Matthias Krack (28.11.2005,MK)
! **************************************************************************************************
   SUBROUTINE reallocate_c2(p, lb1_new, ub1_new, lb2_new, ub2_new)

      COMPLEX(KIND=dp), DIMENSION(:, :), POINTER         :: p
      INTEGER, INTENT(IN)                                :: lb1_new, ub1_new, lb2_new, ub2_new

      CHARACTER(LEN=*), PARAMETER :: routineN = 'reallocate_c2', routineP = moduleN//':'//routineN
      COMPLEX(KIND=dp), PARAMETER                        :: zero = (0.0_dp, 0.0_dp)
      INTEGER, PARAMETER                                 :: t_size = 2*dp_size

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: work

$: reallocate(rank=2)

   END SUBROUTINE reallocate_c2

! **************************************************************************************************
!> \brief (Re)Allocate a complex array 3D with a new dimension
!> \param p ...
!> \param lb1_new ...
!> \param ub1_new ...
!> \param lb2_new ...
!> \param ub2_new ...
!> \param lb3_new ...
!> \param ub3_new ...
!> \par History
!>      none
!> \author Matthias Krack (28.11.2005,MK)
! **************************************************************************************************
   SUBROUTINE reallocate_c3(p, lb1_new, ub1_new, lb2_new, ub2_new, lb3_new, ub3_new)

      COMPLEX(KIND=dp), DIMENSION(:, :, :), POINTER      :: p
      INTEGER, INTENT(IN)                                :: lb1_new, ub1_new, lb2_new, ub2_new, &
                                                            lb3_new, ub3_new

      CHARACTER(LEN=*), PARAMETER :: routineN = 'reallocate_c3', routineP = moduleN//':'//routineN
      COMPLEX(KIND=dp), PARAMETER                        :: zero = (0.0_dp, 0.0_dp)
      INTEGER, PARAMETER                                 :: t_size = 2*dp_size

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: work

$: reallocate(rank=3)

   END SUBROUTINE reallocate_c3

! **************************************************************************************************
!> \brief (Re)Allocate a complex array 4D with a new dimension
!> \param p ...
!> \param lb1_new ...
!> \param ub1_new ...
!> \param lb2_new ...
!> \param ub2_new ...
!> \param lb3_new ...
!> \param ub3_new ...
!> \param lb4_new ...
!> \param ub4_new ...
!> \par History
!>      none
!> \author Matthias Krack (28.11.2005,MK)
! **************************************************************************************************
   SUBROUTINE reallocate_c4(p, lb1_new, ub1_new, lb2_new, ub2_new, lb3_new, ub3_new, &
                            lb4_new, ub4_new)

      COMPLEX(KIND=dp), DIMENSION(:, :, :, :), POINTER   :: p
      INTEGER, INTENT(IN)                                :: lb1_new, ub1_new, lb2_new, ub2_new, &
                                                            lb3_new, ub3_new, lb4_new, ub4_new

      CHARACTER(LEN=*), PARAMETER :: routineN = 'reallocate_c4', routineP = moduleN//':'//routineN
      COMPLEX(KIND=dp), PARAMETER                        :: zero = (0.0_dp, 0.0_dp)
      INTEGER, PARAMETER                                 :: t_size = 2*dp_size

      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :)                           :: work

$: reallocate(rank=4)

   END SUBROUTINE reallocate_c4

! **************************************************************************************************
!> \brief (Re)Allocate an integer vector with a new dimension.
!> \param p ...
!> \param lb1_new ...
!> \param ub1_new ...
!> \par History
!>      none
!> \author Matthias Krack (18.07.2002,MK)
! **************************************************************************************************
   SUBROUTINE reallocate_i1(p, lb1_new, ub1_new)

      INTEGER, DIMENSION(:), POINTER                     :: p
      INTEGER, INTENT(IN)                                :: lb1_new, ub1_new

      CHARACTER(LEN=*), PARAMETER :: routineN = 'reallocate_i1', routineP = moduleN//':'//routineN
      INTEGER, PARAMETER                                 :: t_size = int_size, zero = 0

      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: work

$: reallocate(rank=1)

   END SUBROUTINE reallocate_i1

! **************************************************************************************************
!> \brief (Re)Allocate an integer array with a new dimension.
!> \param p ...
!> \param lb1_new ...
!> \param ub1_new ...
!> \param lb2_new ...
!> \param ub2_new ...
!> \par History
!>      none
!> \author Matthias Krack (18.07.2002,MK)
! **************************************************************************************************
   SUBROUTINE reallocate_i2(p, lb1_new, ub1_new, lb2_new, ub2_new)

      INTEGER, DIMENSION(:, :), POINTER                  :: p
      INTEGER, INTENT(IN)                                :: lb1_new, ub1_new, lb2_new, ub2_new

      CHARACTER(LEN=*), PARAMETER :: routineN = 'reallocate_i2', routineP = moduleN//':'//routineN
      INTEGER, PARAMETER                                 :: t_size = int_size, zero = 0

      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: work

$: reallocate(rank=2)

   END SUBROUTINE reallocate_i2

! **************************************************************************************************
!> \brief (Re)Allocate an integer array 3D with a new dimension.
!> \param p ...
!> \param lb1_new ...
!> \param ub1_new ...
!> \param lb2_new ...
!> \param ub2_new ...
!> \param lb3_new ...
!> \param ub3_new ...
!> \par History
!>      none
!> \author Matthias Krack (18.07.2002,MK)
! **************************************************************************************************
   SUBROUTINE reallocate_i3(p, lb1_new, ub1_new, lb2_new, ub2_new, lb3_new, ub3_new)

      INTEGER, DIMENSION(:, :, :), POINTER               :: p
      INTEGER, INTENT(IN)                                :: lb1_new, ub1_new, lb2_new, ub2_new, &
                                                            lb3_new, ub3_new

      CHARACTER(LEN=*), PARAMETER :: routineN = 'reallocate_i3', routineP = moduleN//':'//routineN
      INTEGER, PARAMETER                                 :: t_size = int_size, zero = 0

      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: work

$: reallocate(rank=3)

   END SUBROUTINE reallocate_i3

! **************************************************************************************************
!> \brief (Re)Allocate an integer array 3D with a new dimension.
!> \param p ...
!> \param lb1_new ...
!> \param ub1_new ...
!> \param lb2_new ...
!> \param ub2_new ...
!> \param lb3_new ...
!> \param ub3_new ...
!> \param lb4_new ...
!> \param ub4_new ...
!> \par History
!>      none
!> \author Matthias Krack (04.10.2002,MK)
! **************************************************************************************************
   SUBROUTINE reallocate_i4(p, lb1_new, ub1_new, lb2_new, ub2_new, lb3_new, ub3_new, &
                            lb4_new, ub4_new)

      INTEGER, DIMENSION(:, :, :, :), POINTER            :: p
      INTEGER, INTENT(IN)                                :: lb1_new, ub1_new, lb2_new, ub2_new, &
                                                            lb3_new, ub3_new, lb4_new, ub4_new

      CHARACTER(LEN=*), PARAMETER :: routineN = 'reallocate_i4', routineP = moduleN//':'//routineN
      INTEGER, PARAMETER                                 :: t_size = int_size, zero = 0

      INTEGER, ALLOCATABLE, DIMENSION(:, :, :, :)        :: work

$: reallocate(rank=4)

   END SUBROUTINE reallocate_i4

! **************************************************************************************************
!> \brief (Re)Allocate an integer (int_8) vector with a new dimension.
!> \param p ...
!> \param lb1_new ...
!> \param ub1_new ...
!> \par History
!>      none
!> \author Matthias Krack (18.07.2002,MK)
! **************************************************************************************************
   SUBROUTINE reallocate_8i1(p, lb1_new, ub1_new)

      INTEGER(KIND=int_8), DIMENSION(:), POINTER         :: p
      INTEGER, INTENT(IN)                                :: lb1_new, ub1_new

      CHARACTER(LEN=*), PARAMETER :: routineN = 'reallocate_8i1', routineP = moduleN//':'//routineN
      INTEGER(KIND=int_8), PARAMETER                     :: zero = 0
      INTEGER, PARAMETER                                 :: t_size = int_8_size

      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:)     :: work

$: reallocate(rank=1)

   END SUBROUTINE reallocate_8i1

! **************************************************************************************************
!> \brief (Re)Allocate an integer (int_8) array with a new dimension.
!> \param p ...
!> \param lb1_new ...
!> \param ub1_new ...
!> \param lb2_new ...
!> \param ub2_new ...
!> \par History
!>      none
!> \author Matthias Krack (18.07.2002,MK)
! **************************************************************************************************
   SUBROUTINE reallocate_8i2(p, lb1_new, ub1_new, lb2_new, ub2_new)

      INTEGER(kind=int_8), DIMENSION(:, :), POINTER      :: p
      INTEGER, INTENT(IN)                                :: lb1_new, ub1_new, lb2_new, ub2_new

      CHARACTER(LEN=*), PARAMETER :: routineN = 'reallocate_8i2', routineP = moduleN//':'//routineN
      INTEGER(KIND=int_8), PARAMETER                     :: zero = 0
      INTEGER, PARAMETER                                 :: t_size = int_8_size

      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:, :)  :: work

$: reallocate(rank=2)

   END SUBROUTINE reallocate_8i2

! **************************************************************************************************
!> \brief (Re)Allocate an real vector with a new dimension.
!> \param p ...
!> \param lb1_new ...
!> \param ub1_new ...
!> \par History
!>      none
!> \author Matthias Krack (18.07.2002,MK)
! **************************************************************************************************
   SUBROUTINE reallocate_r1(p, lb1_new, ub1_new)

      REAL(KIND=dp), DIMENSION(:), POINTER               :: p
      INTEGER, INTENT(IN)                                :: lb1_new, ub1_new

      CHARACTER(LEN=*), PARAMETER :: routineN = 'reallocate_r1', routineP = moduleN//':'//routineN
      INTEGER, PARAMETER                                 :: t_size = dp_size
      REAL(KIND=dp), PARAMETER                           :: zero = 0.0_dp

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: work

$: reallocate(rank=1)

   END SUBROUTINE reallocate_r1

! **************************************************************************************************
!> \brief (Re)Allocate an real array with new dimensions.
!> \param p ...
!> \param lb1_new ...
!> \param ub1_new ...
!> \param lb2_new ...
!> \param ub2_new ...
!> \par History
!>      none
!> \author Matthias Krack (18.07.2002,MK)
! **************************************************************************************************
   SUBROUTINE reallocate_r2(p, lb1_new, ub1_new, lb2_new, ub2_new)

      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: p
      INTEGER, INTENT(IN)                                :: lb1_new, ub1_new, lb2_new, ub2_new

      CHARACTER(LEN=*), PARAMETER :: routineN = 'reallocate_r2', routineP = moduleN//':'//routineN
      INTEGER, PARAMETER                                 :: t_size = dp_size
      REAL(KIND=dp), PARAMETER                           :: zero = 0.0_dp

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: work

$: reallocate(rank=2)

   END SUBROUTINE reallocate_r2

! **************************************************************************************************
!> \brief (Re)Allocate an real array 3D with new dimensions.
!> \param p ...
!> \param lb1_new ...
!> \param ub1_new ...
!> \param lb2_new ...
!> \param ub2_new ...
!> \param lb3_new ...
!> \param ub3_new ...
!> \par History
!>      none
!> \author Matthias Krack (18.07.2002,MK)
! **************************************************************************************************
   SUBROUTINE reallocate_r3(p, lb1_new, ub1_new, lb2_new, ub2_new, lb3_new, ub3_new)

      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: p
      INTEGER, INTENT(IN)                                :: lb1_new, ub1_new, lb2_new, ub2_new, &
                                                            lb3_new, ub3_new

      CHARACTER(LEN=*), PARAMETER :: routineN = 'reallocate_r3', routineP = moduleN//':'//routineN
      INTEGER, PARAMETER                                 :: t_size = dp_size
      REAL(KIND=dp), PARAMETER                           :: zero = 0.0_dp

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: work

$: reallocate(rank=3)

   END SUBROUTINE reallocate_r3

! **************************************************************************************************
!> \brief (Re)Allocate an real array 4D with new dimensions.
!> \param p ...
!> \param lb1_new ...
!> \param ub1_new ...
!> \param lb2_new ...
!> \param ub2_new ...
!> \param lb3_new ...
!> \param ub3_new ...
!> \param lb4_new ...
!> \param ub4_new ...
!> \par History
!>      none
!> \author Matthias Krack (04.10.2002,MK)
! **************************************************************************************************
   SUBROUTINE reallocate_r4(p, lb1_new, ub1_new, lb2_new, ub2_new, lb3_new, ub3_new, &
                            lb4_new, ub4_new)

      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER      :: p
      INTEGER, INTENT(IN)                                :: lb1_new, ub1_new, lb2_new, ub2_new, &
                                                            lb3_new, ub3_new, lb4_new, ub4_new

      CHARACTER(LEN=*), PARAMETER :: routineN = 'reallocate_r4', routineP = moduleN//':'//routineN
      INTEGER, PARAMETER                                 :: t_size = dp_size
      REAL(KIND=dp), PARAMETER                           :: zero = 0.0_dp

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: work

$: reallocate(rank=4)

   END SUBROUTINE reallocate_r4

! **************************************************************************************************
!> \brief (Re)Allocate an real array 5D with new dimensions.
!> \param p ...
!> \param lb1_new ...
!> \param ub1_new ...
!> \param lb2_new ...
!> \param ub2_new ...
!> \param lb3_new ...
!> \param ub3_new ...
!> \param lb4_new ...
!> \param ub4_new ...
!> \param lb5_new ...
!> \param ub5_new ...
!> \par History
!>      none
!> \author Matthias Krack (04.10.2002,MK)
! **************************************************************************************************
   SUBROUTINE reallocate_r5(p, lb1_new, ub1_new, lb2_new, ub2_new, lb3_new, ub3_new, &
                            lb4_new, ub4_new, lb5_new, ub5_new)

      REAL(KIND=dp), DIMENSION(:, :, :, :, :), POINTER   :: p
      INTEGER, INTENT(IN)                                :: lb1_new, ub1_new, lb2_new, ub2_new, &
                                                            lb3_new, ub3_new, lb4_new, ub4_new, &
                                                            lb5_new, ub5_new

      CHARACTER(LEN=*), PARAMETER :: routineN = 'reallocate_r5', routineP = moduleN//':'//routineN
      INTEGER, PARAMETER                                 :: t_size = dp_size
      REAL(KIND=dp), PARAMETER                           :: zero = 0.0_dp

      REAL(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :, :)                        :: work

$: reallocate(rank=5)

   END SUBROUTINE reallocate_r5

! **************************************************************************************************
!> \brief (Re)Allocate a logical vector with a new dimension.
!> \param p ...
!> \param lb1_new ...
!> \param ub1_new ...
!> \par History
!>      none
!> \author Matthias Krack (18.07.2002,MK)
! **************************************************************************************************
   SUBROUTINE reallocate_l1(p, lb1_new, ub1_new)

      LOGICAL, DIMENSION(:), POINTER                     :: p
      INTEGER, INTENT(IN)                                :: lb1_new, ub1_new

      CHARACTER(LEN=*), PARAMETER :: routineN = 'reallocate_l1', routineP = moduleN//':'//routineN
      INTEGER, PARAMETER                                 :: t_size = int_size
      LOGICAL, PARAMETER                                 :: zero = .FALSE.

      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: work

$: reallocate(rank=1)

   END SUBROUTINE reallocate_l1

! **************************************************************************************************
!> \brief (Re)Allocate a vector of string variables
!> \param p_short ...
!> \param lb_new ...
!> \param ub_new ...
!> \param p_long ...
!> \par History
!>      none
!> \author Matthias Krack (18.07.2002,MK)
!> \note
!>      (Maybe outdated) This routine doesnt work on SUN/Solaris!!!
!>                       It should probably not be used. (Thomas Chassaing)
! **************************************************************************************************
   SUBROUTINE reallocate_s1(p_short, lb_new, ub_new, p_long)

      CHARACTER(LEN=default_string_length), &
         DIMENSION(:), OPTIONAL, POINTER                 :: p_short
      INTEGER, INTENT(IN)                                :: lb_new, ub_new
      CHARACTER(LEN=default_path_length), DIMENSION(:), &
         OPTIONAL, POINTER                               :: p_long

      CHARACTER(LEN=*), PARAMETER :: routineN = 'reallocate_s1', routineP = moduleN//':'//routineN

      IF (PRESENT(p_short)) THEN
         CALL reallocate_ss1(p_short, lb_new, ub_new)
      ELSE IF (PRESENT(p_long)) THEN
         CALL reallocate_ls1(p_long, lb_new, ub_new)
      ELSE
         CPABORT("At least one of the two optional arguments is required")
      END IF

   END SUBROUTINE reallocate_s1

! **************************************************************************************************
!> \brief (Re)Allocate a vector of string variables of default string length.
!> \param p ...
!> \param lb1_new ...
!> \param ub1_new ...
!> \par History
!>      none
!> \author Thomas Chassaing
!> \author Matthias Krack (18.07.2002,MK)
!> \note
!>      (Maybe outdated) This routine doesnt work on SUN/Solaris!!!
!>                       It should probably not be used. (Thomas Chassaing)
! **************************************************************************************************
   SUBROUTINE reallocate_ss1(p, lb1_new, ub1_new)

      CHARACTER(LEN=default_string_length), &
         DIMENSION(:), POINTER                           :: p
      INTEGER, INTENT(IN)                                :: lb1_new, ub1_new

      CHARACTER(LEN=*), PARAMETER :: routineN = 'reallocate_ss1', routineP = moduleN//':'//routineN
      CHARACTER(LEN=default_string_length), PARAMETER    :: zero = ""
      INTEGER, PARAMETER :: t_size = int_size*default_string_length

      CHARACTER(LEN=default_string_length), &
         ALLOCATABLE, DIMENSION(:)                       :: work

$: reallocate(rank=1)

   END SUBROUTINE reallocate_ss1

! **************************************************************************************************
!> \brief (Re)Allocate a vector of string variables of default path length.
!> \param p ...
!> \param lb1_new ...
!> \param ub1_new ...
!> \par History
!>      none
!> \author Matthias Krack (18.07.2002,MK)
!> \note
!>      (Maybe outdated) This routine doesnt work on SUN/Solaris!!!
!>                       It should probably not be used. (Thomas Chassaing)
! **************************************************************************************************
   SUBROUTINE reallocate_ls1(p, lb1_new, ub1_new)

      CHARACTER(LEN=default_path_length), DIMENSION(:), &
         POINTER                                         :: p
      INTEGER, INTENT(IN)                                :: lb1_new, ub1_new

      CHARACTER(LEN=*), PARAMETER :: routineN = 'reallocate_ls1', routineP = moduleN//':'//routineN
      CHARACTER(LEN=default_path_length), PARAMETER      :: zero = ""
      INTEGER, PARAMETER :: t_size = int_size*default_path_length

      CHARACTER(LEN=default_path_length), ALLOCATABLE, &
         DIMENSION(:)                                    :: work

$: reallocate(rank=1)

   END SUBROUTINE reallocate_ls1

END MODULE memory_utilities
