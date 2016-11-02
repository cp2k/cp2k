!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Wrapper to allocate and optionally assign array from source array.
!>        Will be redundant once compilers implement Fortran 2008 standard :-/
!> \author Patrick Seewald
! **************************************************************************************************
MODULE allocate_wrap
   USE allocate_wrap_c,                 ONLY: allocate_any
   USE allocate_wrap_d,                 ONLY: allocate_any
   USE allocate_wrap_i,                 ONLY: allocate_any
   USE allocate_wrap_s,                 ONLY: allocate_any
   USE allocate_wrap_z,                 ONLY: allocate_any

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: allocate_any

END MODULE
