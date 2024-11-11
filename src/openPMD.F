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


END MODULE openPMD_api
