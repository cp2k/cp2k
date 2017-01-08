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
MODULE atomic_kind_list_types
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              deallocate_atomic_kind_set
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

#:include "subsys_list.fypp"

$:subsys_list_macro(shortname='atomic_kind', longname='atomic_kind')

END MODULE atomic_kind_list_types
