!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief contains the structure
!> \par History
!>      11.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
MODULE xc_rho_cflags_types

#include "../base/base_uses.f90"
   IMPLICIT NONE
   PRIVATE
   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .FALSE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_rho_cflags_types'

   PUBLIC :: xc_rho_cflags_type
   PUBLIC :: xc_rho_cflags_setall, &
             xc_rho_cflags_equal

! **************************************************************************************************
!> \brief contains a flag for each component of xc_rho_set, so that you can
!>      use it to tell which components you need, which ones you need,....
!> \param rho flags for rho (the total density)
!> \param rho_spin flag for rhoa and rhob (the spin density with LSD)
!> \param drho flag for drho (the gradient of rho)
!> \param drho_spin flag for drhoa and drhob (the gradient of the spin
!>        density)
!> \param norm_drho flag for norm_drho (the norm of the gradient of rho)
!> \param norm_drho_spin flag for norm_drhoa, norm_drhob (the norm of the
!>        gradient of the spin density)
!> \param drhoa_drhob flag for drhoa_drhob (the scalar product of the
!>        gradient of the two spin densities)
!> \param rho_ 1_3: flag for rho**(1.0_dp/3.0_dp)
!> \param rho_spin_ 1_3: flag for rhoa**(1.0_dp/3.0_dp) and rhob**(1.0_dp/3.0_dp)
!> \param tau flags for the kinetic (KS) part of rho
!> \param tau_spin flags for the kinetic (KS) part of rhoa and rhob
!> \note
!>      low_level type without retain/release
!> \par History
!>      11.2003 created [fawzi]
!>      12.2008 added laplace parts [mguidon]
!> \author fawzi
! **************************************************************************************************
   TYPE xc_rho_cflags_type
      LOGICAL :: rho, rho_spin, drho, drho_spin, &
                 norm_drho, norm_drho_spin, drhoa_drhob, &
                 rho_1_3, rho_spin_1_3, tau, tau_spin, laplace_rho, laplace_rho_spin
   END TYPE xc_rho_cflags_type

CONTAINS

! **************************************************************************************************
!> \brief sets all the flags to the given value
!> \param cflags the flags to set
!> \param value the value to set
! **************************************************************************************************
   SUBROUTINE xc_rho_cflags_setall(cflags, value)
      TYPE(xc_rho_cflags_type), INTENT(out)              :: cflags
      LOGICAL, INTENT(in)                                :: value

      CHARACTER(len=*), PARAMETER :: routineN = 'xc_rho_cflags_setall', &
         routineP = moduleN//':'//routineN

      cflags%rho = value
      cflags%rho_spin = value
      cflags%drho = value
      cflags%drho_spin = value
      cflags%norm_drho = value
      cflags%norm_drho_spin = value
      cflags%drhoa_drhob = value
      cflags%rho_1_3 = value
      cflags%rho_spin_1_3 = value
      cflags%tau = value
      cflags%tau_spin = value
      cflags%laplace_rho = value
      cflags%laplace_rho_spin = value
   END SUBROUTINE xc_rho_cflags_setall

! **************************************************************************************************
!> \brief return true if the two cflags are equal
!> \param cflags1 the flags to compare
!> \param cflags2 the flags to compare
!> \retval equal ...
! **************************************************************************************************
   FUNCTION xc_rho_cflags_equal(cflags1, cflags2) RESULT(equal)
      TYPE(xc_rho_cflags_type), INTENT(inout)            :: cflags1
      TYPE(xc_rho_cflags_type), INTENT(in)               :: cflags2
      LOGICAL                                            :: equal

      CHARACTER(len=*), PARAMETER :: routineN = 'xc_rho_cflags_equal', &
         routineP = moduleN//':'//routineN

      equal = ((cflags1%rho .EQV. cflags2%rho) .AND. &
               (cflags1%rho_spin .EQV. cflags2%rho_spin) .AND. &
               (cflags1%drho .EQV. cflags2%drho) .AND. &
               (cflags1%drho_spin .EQV. cflags2%drho_spin) .AND. &
               (cflags1%norm_drho .EQV. cflags2%norm_drho) .AND. &
               (cflags1%norm_drho_spin .EQV. cflags2%norm_drho_spin) .AND. &
               (cflags1%drhoa_drhob .EQV. cflags2%drhoa_drhob) .AND. &
               (cflags1%rho_1_3 .EQV. cflags2%rho_1_3) .AND. &
               (cflags1%rho_spin_1_3 .EQV. cflags2%rho_spin_1_3) .AND. &
               (cflags1%tau .EQV. cflags2%tau) .AND. &
               (cflags1%tau_spin .EQV. cflags2%tau_spin) .AND. &
               (cflags1%laplace_rho .EQV. cflags2%laplace_rho) .AND. &
               (cflags1%laplace_rho_spin .EQV. cflags2%laplace_rho_spin))

   END FUNCTION xc_rho_cflags_equal

END MODULE xc_rho_cflags_types
