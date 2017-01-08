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
MODULE mol_new_list_types
   USE molecule_types_new,              ONLY: deallocate_molecule_set,&
                                              molecule_type

#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

#:include "subsys_list.fypp"

$:subsys_list_macro(shortname='mol_new', longname='molecule')

END MODULE mol_new_list_types
