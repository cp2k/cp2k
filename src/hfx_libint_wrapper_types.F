!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Provides types that are used to interface with the libint library
!> \par History
!>      10.2009 created [Manuel Guidon]
!> \author Manuel Guidon
! **************************************************************************************************
MODULE hfx_libint_wrapper_types

   USE libint_f,                        ONLY: libint2_max_am_eri,&
                                              libint_t

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: libint_t, &
             prim_data_f_size, &
             build_eri_size, &
             build_deriv1_eri_size, &
             libint_max_am, &
             libderiv_max_am1

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_libint_wrapper_types'

   INTEGER, PARAMETER :: libint_max_am = libint2_max_am_eri
   INTEGER, PARAMETER :: libderiv_max_am1 = libint_max_am
   INTEGER, PARAMETER :: prim_data_f_size = 4*(libint_max_am)+1
   INTEGER, PARAMETER :: libint_vrr_classes_size = 2*(libint_max_am)+1
   INTEGER, PARAMETER :: libint_dvrr_classes_size = 2*(libderiv_max_am1)+1
   INTEGER, PARAMETER :: build_eri_size = libint_max_am
   INTEGER, PARAMETER :: build_deriv1_eri_size = libderiv_max_am1

END MODULE hfx_libint_wrapper_types
