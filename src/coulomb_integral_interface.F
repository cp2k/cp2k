!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Common interface for two- and three-center Coulomb integrals.
! **************************************************************************************************
MODULE coulomb_integral_interface
   USE integral_library_types,          ONLY: coulomb_operator_type
   USE libint_2c_3c,                    ONLY: compute_coulomb_2c,&
                                              compute_coulomb_2c_derivs,&
                                              compute_coulomb_3c,&
                                              compute_coulomb_3c_derivs,&
                                              coulomb_integral_cleanup,&
                                              coulomb_integral_context_type,&
                                              coulomb_integral_init,&
                                              coulomb_integral_max_am

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: coulomb_operator_type, coulomb_integral_context_type, &
             coulomb_integral_init, coulomb_integral_cleanup, coulomb_integral_max_am, &
             compute_coulomb_2c, compute_coulomb_3c, &
             compute_coulomb_2c_derivs, compute_coulomb_3c_derivs

END MODULE coulomb_integral_interface
