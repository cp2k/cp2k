!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

#:include '../common/cp_linked_list.fypp'

MODULE cp_linked_list_input
   USE kinds, only: dp, default_string_length
   USE input_val_types, only: val_type, val_p_type
#include "../base/base_uses.f90"

#:mute
#:set nametype1 = ['int', 'real', 'logical', 'char', 'val']
#:set type1 = ['integer', 'REAL(kind=dp)','logical', 'character(len=default_string_length)', 'type(val_type),pointer']
#:set type1in = [_ + ', intent(in)' for _ in type1]
#:set type1out = [_ + ', intent(out)' for _ in type1]

#:set eq = ['=', '=', '=', '=', '=>']
#:set arrayeq = eq
#:set type1arrayEl = ['integer', 'REAL(kind=dp)', 'logical', 'character(len=default_string_length)','type(val_p_type)']
#:set arrayEl = ['', '', '', '', '%val']
#:set private_routines = ''
#:endmute

$:inst(nametype1, type1, type1in, type1out, eq, arrayeq, type1arrayEl, arrayEl, private_routines)
END MODULE
