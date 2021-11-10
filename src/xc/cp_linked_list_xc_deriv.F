!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

#:include '../common/cp_linked_list.fypp'

MODULE cp_linked_list_xc_deriv
   USE xc_derivative_types, ONLY: xc_derivative_type, xc_derivative_p_type
#include "../base/base_uses.f90"

   #:mute
      #:set nametype1 = ['xc_deriv']
      #:set type1 = ['TYPE(xc_derivative_type), POINTER']
      #:set type1in = type1
      #:set type1out = type1

      #:set eq = ['=>']
      #:set arrayeq = eq
      #:set type1arrayEl = ['type(xc_derivative_p_type)']
      #:set arrayEl = ['%deriv']

      #:def get_private_routines()
! **************************************************************************************************
!> \brief private compare function
!> \param el1 ...
!> \param el2 ...
!> \return ...
! **************************************************************************************************
         function cp_sll_deriv_less_q(el1, el2) result(res)
            type(xc_derivative_type), INTENT(IN) :: el1, el2
            logical :: res

            CPASSERT(associated(el1%split_desc))
            CPASSERT(associated(el2%split_desc))
     res = size(el1%split_desc) < size(el2%split_desc) .or. (size(el1%split_desc) == size(el2%split_desc) .and. el1%desc < el2%desc)
         end function
      #:enddef

      #:set private_routines = get_private_routines()
   #:endmute

   $:inst(nametype1, type1, type1in, type1out, eq, arrayeq, type1arrayEl, arrayEl, private_routines)
END MODULE
