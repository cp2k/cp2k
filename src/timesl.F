!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/timesl [1.0] *
!!
!!   NAME
!!     timesl
!!
!!   FUNCTION
!!     Provides cputime, wallclock and date information
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE timesl
  
  USE kinds, ONLY : dbl, sgl
  USE machine, ONLY : m_cputime, m_walltime, m_datum
  
  PRIVATE
  PUBLIC :: walltime, cputime, datum

!******************************************************************************
  
CONTAINS

!******************************************************************************

FUNCTION cputime() RESULT ( ct )
  
  IMPLICIT NONE
  
  REAL ( dbl ) :: ct
  
  ct = m_cputime()
END FUNCTION cputime

!******************************************************************************

FUNCTION walltime() RESULT ( wt )
  
  IMPLICIT NONE
  
  REAL ( dbl ) :: wt
  
  wt = m_walltime()
  
END FUNCTION walltime

!******************************************************************************

SUBROUTINE datum ( cal_date )
  
  IMPLICIT NONE
  
  CHARACTER ( LEN = * ) :: cal_date
  
  CALL m_datum ( cal_date )
  
END SUBROUTINE datum

!******************************************************************************

END MODULE timesl

!******************************************************************************
!!*****
