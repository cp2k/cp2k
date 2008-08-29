!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief wrappers for timings for MPI calls: way to fix the circular dependency
!> \author Teodoro Laino [tlaino] - University of Zurich - 08.2008
! *****************************************************************************

INTERFACE
   SUBROUTINE timeset_mp(name, handle)
     CHARACTER(LEN=*), INTENT(IN)             :: name
     INTEGER                                  :: handle

   END SUBROUTINE timeset_mp

   SUBROUTINE timestop_mp(handle)
     INTEGER                                  :: handle

   END SUBROUTINE timestop_mp
END INTERFACE
