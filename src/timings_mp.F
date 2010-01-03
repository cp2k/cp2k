!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief wrappers for timings for MPI calls: way to fix the circular dependency
!> \author Teodoro Laino [tlaino] - University of Zurich - 08.2008
! *****************************************************************************

SUBROUTINE timeset_mp(name, handle)
  USE timings,                         ONLY:  timeset
  CHARACTER(LEN=*), INTENT(IN)             :: name
  INTEGER                                  :: handle

  CALL timeset(name, handle)
END SUBROUTINE timeset_mp

SUBROUTINE timestop_mp(handle)
  USE timings,                         ONLY:  timestop
  INTEGER                                  :: handle
  
  CALL timestop(handle)
END SUBROUTINE timestop_mp
