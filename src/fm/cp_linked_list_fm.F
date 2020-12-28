!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

#:include '../common/cp_linked_list.fypp'

MODULE cp_linked_list_fm
   USE cp_fm_types, ONLY: cp_fm_type, cp_fm_p_type
#include "base/base_uses.f90"

#:mute
#:set nametype1 = ['fm']
#:set type1 = ['type(cp_fm_type),pointer']
#:set type1in = type1
#:set type1out = type1

#:set eq = ['=>']
#:set arrayeq = eq
#:set type1arrayEl = ['type(cp_fm_p_type)']
#:set arrayEl = ['%matrix']
#:set private_routines = ''
#:endmute

$:inst(nametype1, type1, type1in, type1out, eq, arrayeq, type1arrayEl, arrayEl, private_routines)
END MODULE
