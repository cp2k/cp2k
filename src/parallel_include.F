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
!!   SOURCE
!******************************************************************************

MODULE parallel_include

  IMPLICIT NONE

#if defined(__parallel)
  INCLUDE "mpif.h"

#endif
  PUBLIC

  INTEGER :: dummy

END MODULE parallel_include
