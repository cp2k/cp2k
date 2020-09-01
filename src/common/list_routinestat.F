!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief An array-based list which grows on demand.
!>        When the internal array is full, a new array of twice the size will be
!>        allocated and the items are copied over.
!>
!>        This list can also be used as a stack.
!>        Have look at list_push(), list_pop() and list_peek().
!> \par History
!>      12.2012 first version [ole]
!> \author Ole Schuett
! **************************************************************************************************

#:include 'list.fypp'
#:set n_inst = 1
MODULE list_routinestat
   USE timings_base_type, ONLY: ${uselist(usetimingsbase)}$

$:list_body(valuetype[n_inst], valuetype_in[n_inst], valuetype_out[n_inst], value_assign[n_inst])

END MODULE list_routinestat
