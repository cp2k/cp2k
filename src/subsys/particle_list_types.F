!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief represent a simple array based list of the given type
!> \par History
!>      08.2003 created [fawzi]
!>      01.2017 ported to Fypp [Ole Schuett]
!> \author Fawzi Mohamed
! **************************************************************************************************
MODULE particle_list_types
   USE particle_types,                  ONLY: deallocate_particle_set,&
                                              particle_type

#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

#:include "subsys_list.fypp"

$:subsys_list_macro(shortname='particle', longname='particle')

END MODULE particle_list_types
