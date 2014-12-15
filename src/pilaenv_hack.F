#if defined(__PILAENV_BLOCKSIZE)
!
! This function overwrites the one present in Scalapack.....
! most libraries do not provide a good default, but an old value of 32, which leads to poor pdgemm performance
! set -D__PILAENV_BLOCKSIZE=1024 or some number (depending on the hardware, and application), might lead to better performance
!
! *****************************************************************************
!> \brief ...
!> \param ICTXT ...
!> \param PREC ...
!> \retval PILAENV ...
! *****************************************************************************
INTEGER FUNCTION PILAENV (ICTXT, PREC) 
  INTEGER :: ICTXT
  CHARACTER(LEN=1) :: PREC

  PILAENV = __PILAENV_BLOCKSIZE
END FUNCTION PILAENV
#endif

! *****************************************************************************
!> \brief ...
! *****************************************************************************
SUBROUTINE NAG_dummy()

END SUBROUTINE NAG_dummy
