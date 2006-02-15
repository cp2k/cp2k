!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002 CP2K developers group                                  !
!-----------------------------------------------------------------------------!
!!****** cp2k/parallel_include [1.0] *
!!
!!   NAME
!!     parallel_include
!!
!!   FUNCTION
!!     Loads external information for the Message Passing Systems.
!!
!!   AUTHOR
!!     JGH,MK
!!
!!   MODIFICATION HISTORY
!!     - Changed from fixed to free format (13.06.2002,MK)
!!
!!*** *************************************************************************

MODULE parallel_include

  IMPLICIT NONE

#if defined(__parallel)
  INCLUDE "mpif.h"
#if defined (__SGL)
  INTEGER     :: MP_STD_REAL = MPI_REAL
  INTEGER     :: MP_STD_COMPLEX = MPI_COMPLEX
#else
  INTEGER     :: MP_STD_REAL = MPI_DOUBLE_PRECISION
  INTEGER     :: MP_STD_COMPLEX = MPI_DOUBLE_COMPLEX
#endif
#endif
  PUBLIC

  INTEGER :: dummy

END MODULE parallel_include
