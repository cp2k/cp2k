!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!
! *****************************************************************************
!> \brief an exclusion type
! *****************************************************************************
MODULE exclusion_types
  IMPLICIT NONE
  PRIVATE

! *****************************************************************************
!> \brief A type used to store lists of exclusions and onfos
!> \par History
!>      12.2010 created [Joost VandeVondele]
!> \author fawzi
! *****************************************************************************
  TYPE exclusion_type
     INTEGER, POINTER, DIMENSION(:)        :: list_exclude_vdw
     INTEGER, POINTER, DIMENSION(:)        :: list_exclude_ei
     INTEGER, POINTER, DIMENSION(:)        :: list_onfo
  END TYPE
  
  PUBLIC :: exclusion_type

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'exclusion_types'

END MODULE exclusion_types
