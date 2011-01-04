!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2011  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Defines the basic variable types
!> \author Teodoro Laino - 01.2011
! *****************************************************************************
MODULE lib_kinds
  IMPLICIT NONE
  
  PRIVATE
#if defined(__SGL)
  INTEGER, PARAMETER, PUBLIC :: wp = KIND(0.0)
  INTEGER, PARAMETER, PUBLIC :: dp = KIND(0.0)
  INTEGER, PARAMETER, PUBLIC :: sp = KIND(0.0)
#else
  INTEGER, PARAMETER, PUBLIC :: wp = SELECTED_REAL_KIND ( 14, 200 )
  INTEGER, PARAMETER, PUBLIC :: dp = SELECTED_REAL_KIND ( 14, 200 )
  INTEGER, PARAMETER, PUBLIC :: sp = KIND(0.0)
#endif

END MODULE lib_kinds

