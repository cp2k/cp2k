!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
MODULE openPMD_api
   USE ISO_C_BINDING, ONLY: C_ASSOCIATED, &
                            C_BOOL, &
                            C_CHAR, &
                            C_FLOAT, &
                            C_DOUBLE, &
                            C_F_POINTER, &
                            C_INT, &
                            C_NULL_CHAR, &
                            C_NULL_PTR, &
                            C_PTR, &
                            C_INT64_T

   USE kinds, ONLY: sp, int_8, dp, default_string_length

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

CONTAINS


! **************************************************************************************************
!> \brief Releases a Torch dictionary and all its ressources.
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE torch_dict_release(dict)
      TYPE(torch_dict_type), INTENT(INOUT)               :: dict

      CPABORT("CP2K was compiled without Torch library.")
      MARK_USED(dict)
   END SUBROUTINE torch_dict_release

#if defined(__LIBTORCH)
      INTERFACE
         SUBROUTINE torch_c_model_eval(model, inputs, outputs) BIND(C, name="torch_c_model_eval")
            IMPORT :: C_PTR
            TYPE(C_PTR), VALUE                        :: model
            TYPE(C_PTR), VALUE                        :: inputs
            TYPE(C_PTR), VALUE                        :: outputs
         END SUBROUTINE torch_c_model_eval
      END INTERFACE

      CPASSERT(C_ASSOCIATED(model%c_ptr))
      CPASSERT(C_ASSOCIATED(inputs%c_ptr))
      CPASSERT(C_ASSOCIATED(outputs%c_ptr))
      CALL torch_c_model_eval(model=model%c_ptr, &
                              inputs=inputs%c_ptr, &
                              outputs=outputs%c_ptr)
#else
      CPABORT("CP2K was compiled without Torch library.")
      MARK_USED(model)
      MARK_USED(inputs)
      MARK_USED(outputs)
#endif
   END SUBROUTINE torch_model_eval

END MODULE openPMD_api
