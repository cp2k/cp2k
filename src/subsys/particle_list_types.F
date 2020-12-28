!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief represent a simple array based list of the given type
!> \par History
!>      08.2003 created [fawzi]
!>      01.2017 ported to Fypp [Ole Schuett]
!> \author Fawzi Mohamed
! **************************************************************************************************
MODULE particle_list_types
   USE particle_types, ONLY: deallocate_particle_set, &
                             particle_type

#include "base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

#:include "subsys_list.fypp"

$:subsys_list_macro(name='particle')

END MODULE particle_list_types
