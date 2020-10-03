!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Integrate single or product functions over a potential on a RS grid
!> \par History
!>      Refactored from earlier versions by Joost VandeVondele (2002,2003,2007)
!> \author JGH [04.2014]
! **************************************************************************************************
!
! This module acts as a common container for the routines from the low level
! modules
!           grid_api
!           qs_integrate_potential_product
!           qs_integrate_potential_single
!
! **************************************************************************************************
MODULE qs_integrate_potential
   USE grid_api,                        ONLY: integrate_pgf_product
   USE qs_integrate_potential_product,  ONLY: integrate_v_rspace
   USE qs_integrate_potential_single,   ONLY: integrate_ppl_rspace,&
                                              integrate_rho_nlcc,&
                                              integrate_v_core_rspace,&
                                              integrate_v_rspace_diagonal,&
                                              integrate_v_rspace_one_center
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_integrate_potential'

! *** Public subroutines ***

   ! included from qs_integrate_potential_product
   PUBLIC :: integrate_v_rspace

   ! included from qs_integrate_potential_single
   PUBLIC :: integrate_v_rspace_one_center, &
             integrate_v_rspace_diagonal, &
             integrate_v_core_rspace, &
             integrate_ppl_rspace, &
             integrate_rho_nlcc

   ! included from grid_api
   PUBLIC :: integrate_pgf_product

END MODULE qs_integrate_potential
