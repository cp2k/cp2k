!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

#:include '../common/cp_linked_list.fypp'

MODULE cp_linked_list_pw
   USE kinds, ONLY: dp
   USE pw_types, ONLY: pw_type, pw_p_type
   USE cp_array_utils, ONLY: cp_3d_r_p_type
   USE realspace_grid_types, ONLY: realspace_grid_type, realspace_grid_p_type
#include "../base/base_uses.f90"

   #:mute
      #:set nametype1 = ['3d_r', 'pw', 'rs']
      #:set type1 = ['REAL(kind=dp),dimension(:,:,:),pointer', 'type(pw_type),pointer', 'type(realspace_grid_type),pointer']
      #:set type1in = type1
      #:set type1out = type1

      #:set eq = ['=>','=>','=>','=>']
      #:set arrayeq = eq
      #:set type1arrayEl = ['type(cp_3d_r_p_type)','type(pw_p_type)','type(realspace_grid_p_type)']

      #:set arrayEl = ['%array','%pw','%rs_grid']

      #:set private_routines = ''

      #:set common_dir = '../common'
   #:endmute

   $:inst(nametype1, type1, type1in, type1out, eq, arrayeq, type1arrayEl, arrayEl, private_routines)
END MODULE
