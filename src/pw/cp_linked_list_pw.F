!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

#:include '../common/cp_linked_list.fypp'
#:include 'pw_types.fypp'

MODULE cp_linked_list_pw
   USE kinds, ONLY: dp
   #:for kind in pw_kinds
      USE pw_types, ONLY: pw_${kind}$_type, pw_${kind}$_p_type
   #:endfor
   USE pw_types, ONLY: pw_type, pw_p_type
   #:for kind in pw_kinds
      USE cp_array_utils, ONLY: cp_${kind[1:]}$_${kind[0]}$_cp_type
   #:endfor
   USE realspace_grid_types, ONLY: realspace_grid_type, realspace_grid_p_type
#include "../base/base_uses.f90"

   #:mute
      #:set nametype1 = ['pw', 'rs']
      #:set type1 = ['type(pw_type),pointer', 'type(realspace_grid_type),pointer']
      #:for kind, type in zip(pw_kinds, pw_types)
         $:nametype1.append(kind[1:]+'_'+kind[0])
         $:nametype1.append('pw_'+kind)
         $:type1.append(type+',pointer,contiguous')
         $:type1.append('type(pw_'+kind+'_type),pointer')
      #:endfor
      #:set type1in = type1
      #:set type1out = type1

      #:set eq = ['=>','=>']
      #:set arrayeq = eq
      #:set type1arrayEl = ['type(pw_p_type)','type(realspace_grid_p_type)']
      #:for kind, type in zip(pw_kinds, pw_types)
         $:eq.append('=>')
         $:eq.append('=>')
         $:type1arrayEl.append('type(cp_'+kind[1:]+'_'+kind[0]+'_cp_type)')
         $:type1arrayEl.append('type(pw_'+kind+'_p_type)')
      #:endfor

      #:set arrayEl = ['%pw','%rs_grid']

      #:set private_routines = ''

      #:set default_init = [' => NULL()', ' => NULL()']

      #:for kind, type in zip(pw_kinds, pw_types)
         $:arrayEl.append('%array')
         $:arrayEl.append('%pw')
         $:default_init.append(' => NULL()')
         $:default_init.append(' => NULL()')
      #:endfor

      #:set common_dir = '../common'
   #:endmute

   $:inst(nametype1, type1, type1in, type1out, eq, arrayeq, type1arrayEl, arrayEl, private_routines, default_init)
END MODULE
