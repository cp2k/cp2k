! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the Lesser GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! Lesser GNU General Public License for more details.
!
! You should have received a copy of the Lesser GNU General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

!> Implementation of the spherical-wave expanded (Koide) dispersion damping function.
module dftd4_damping_koide
   use dftd4_damping_type, only : damping_twobody
   use dftd4_param_type, only : param_type
   use mctc_env, only : error_type, fatal_error, wp
   implicit none
   private

   public :: koide_damping_twobody, new_koide_damping_twobody

   !> Spherical-wave expanded (Koide) dispersion damping for two-body terms
   type, extends(damping_twobody) :: koide_damping_twobody
   contains
      !> Evaluate two-body damping factor
      procedure :: get_2b_damp
      !> Evaluate two-body damping factor with derivatives
      procedure :: get_2b_derivs
      !> Check the availability of the parameters for two-body damping
      procedure :: check_2b_params
   end type koide_damping_twobody

   real(wp), parameter :: fact(2:6) = [2.0_wp, 6.0_wp, 24.0_wp, 120.0_wp, 720.0_wp]

   character(len=*), parameter :: koide_label = "Spherical-wave expanded (Koide) damping"
   character(len=*), parameter :: koide_label_short = "Koide"

contains


!> Create a new instance of the Koide damping for two-body terms
subroutine new_koide_damping_twobody(self)
   !> Instance of the two-body Koide damping function
   type(koide_damping_twobody), intent(out) :: self

   self%label = koide_label
   self%label_short = koide_label_short
end subroutine new_koide_damping_twobody


!> Evaluate two-body damping factor
pure subroutine get_2b_damp(self, param, r2, rdamp, d6, d8)
   !> Instance of the two-body damping function
   class(koide_damping_twobody), intent(in) :: self
   !> Damping parameters
   type(param_type), intent(in) :: param
   !> Square of interatomic distance
   real(wp), intent(in) :: r2
   !> Damping radius
   real(wp), intent(in) :: rdamp
   !> Damping factor for C6/R^6 term
   real(wp), intent(out) :: d6
   !> Damping factor for C8/R^8 term
   real(wp), intent(out) :: d8

   real(wp) :: r1, r6, r8, rpoly, R_scaled6, R_scaled8, exp_R
   real(wp) :: P0, P1, P2, P3, phi2_11, phi0_11, chi_11, phi3_12, phi1_12, chi_12
   
   r1 = sqrt(r2)
   r6 = r2 * r2 * r2
   r8 = r6 * r2
   
   ! First-order polynomial scaling of the damping radius
   rpoly = param%a1 * rdamp + param%a2
   
   ! Scaled distance depending on the element damping radius
   R_scaled6 = param%rs6 * r1 / rpoly 
   exp_R = exp(-R_scaled6)

   ! Dipole-dipole interaction
   P2 = (1.0_wp + R_scaled6 * (1.0_wp + R_scaled6 * (1.0_wp/fact(2) &
      & + R_scaled6 * (1.0_wp/fact(3) + R_scaled6 * (1.0_wp/fact(4) &
      & + R_scaled6 * (1699.0_wp/207360.0_wp &
      & + R_scaled6 * (259.0_wp/207360.0_wp &
      & + R_scaled6 * (197.0_wp/1451520.0_wp &
      & + R_scaled6 * (19.0_wp/2177280.0_wp &
      & + R_scaled6 * (1.0_wp/4354560.0_wp))))))))))

   P0 = (R_scaled6**3) * (89.0_wp/13824.0_wp &
      & + R_scaled6 * (89.0_wp/13824.0_wp &
      & + R_scaled6 * (119.0_wp/41472.0_wp &
      & + R_scaled6 * (5.0_wp/6912.0_wp &
      & + R_scaled6 * (11.0_wp/103680.0_wp &
      & + R_scaled6 * (1.0_wp/124416.0_wp &
      & + R_scaled6 * (1.0_wp/4354560.0_wp)))))))

   phi2_11 = 1.0_wp - exp_R * P2
   phi0_11 = exp_R * P0

   chi_11 = phi2_11**2 + 0.5_wp * phi0_11**2

   ! Scaled distance depending on the element damping radius
   R_scaled8 = param%rs8 * r1 / rpoly 
   exp_R = exp(-R_scaled8)

   ! Dipole-quadrupole interaction
   P3 = (1.0_wp + R_scaled8 * (1.0_wp + R_scaled8 * (1.0_wp/fact(2) &
      & + R_scaled8 * (1.0_wp/fact(3) + R_scaled8 * (1.0_wp/fact(4) &
      & + R_scaled8 * (1.0_wp/fact(5) + R_scaled8 * (1.0_wp/fact(6) &
      & + R_scaled8 * (109.0_wp/552960.0_wp &
      & + R_scaled8 * (13.0_wp/552960.0_wp &
      & + R_scaled8 * (211.0_wp/96768000.0_wp &
      & + R_scaled8 * (57.0_wp/435456000.0_wp &
      & + R_scaled8 * (1.0_wp/290304000.0_wp))))))))))))

   P1 = (R_scaled8**5) * (47.0_wp/552960.0_wp &
      & + R_scaled8 * (47.0_wp/552960.0_wp &
      & + R_scaled8 * (7.0_wp/184320.0_wp &
      & + R_scaled8 * (1.0_wp/103680.0_wp &
      & + R_scaled8 * (209.0_wp/145152000.0_wp &
      & + R_scaled8 * (11.0_wp/96768000.0_wp &
      & + R_scaled8 * (1.0_wp/290304000.0_wp)))))))

   phi3_12 = 1.0_wp - exp_R * P3
   phi1_12 = exp_R * P1

   chi_12 = phi3_12**2 + (2.0_wp/3.0_wp) * phi1_12**2

   ! Assemble damping factors
   d6 = param%s6 * chi_11 / r6
   d8 = param%s8 * chi_12 / r8

end subroutine get_2b_damp


!> Evaluate two-body damping factor with derivarives
pure subroutine get_2b_derivs(self, param, r2, rdamp, d6, d8, d6dr, d8dr)
   !> Instance of the two-body damping function
   class(koide_damping_twobody), intent(in) :: self
   !> Damping parameters
   type(param_type), intent(in) :: param
   !> Square of interatomic distance
   real(wp), intent(in) :: r2
   !> Damping radius
   real(wp), intent(in) :: rdamp
   !> Damping factor for C6/R^6 term
   real(wp), intent(out) :: d6
   !> Damping factor for C8/R^8 term
   real(wp), intent(out) :: d8
   !> Derivative of damping factor for C6/R^6 w.r.t the interatomic distance
   real(wp), intent(out) :: d6dr
   !> Derivative of damping factor for C8/R^8 w.r.t the interatomic distance
   real(wp), intent(out) :: d8dr
   
   real(wp) :: r1, r6, r8, rpoly, R_scaled6, R_scaled8, exp_R
   real(wp) :: P2, dP2, phi2, dphi2
   real(wp) :: P0, dP0, phi0, dphi0
   real(wp) :: P3, dP3, phi3, dphi3
   real(wp) :: P1, dP1, phi1, dphi1
   real(wp) :: chi_11, dchi_11, chi_12, dchi_12

   r1 = sqrt(r2)
   r6 = r2 * r2 * r2
   r8 = r6 * r2

   ! First-order polynomial scaling of the damping radius
   rpoly = param%a1 * rdamp + param%a2

   ! Scaled distance depending on the element damping radius
   R_scaled6 = param%rs6 * r1 / rpoly
   exp_R = exp(-R_scaled6)

   ! Dipole-dipole interaction
   P2 = 1.0_wp + R_scaled6 * (1.0_wp + R_scaled6 * (1.0_wp/fact(2) &
      & + R_scaled6 * (1.0_wp/fact(3) + R_scaled6 * (1.0_wp/fact(4) &
      & + R_scaled6 * (1699.0_wp/207360.0_wp &
      & + R_scaled6 * (259.0_wp/207360.0_wp &
      & + R_scaled6 * (197.0_wp/1451520.0_wp &
      & + R_scaled6 * (19.0_wp/2177280.0_wp &
      & + R_scaled6 * (1.0_wp/4354560.0_wp)))))))))

   dP2 = 1.0_wp + R_scaled6 * (1.0_wp + R_scaled6 * (1.0_wp/fact(2) &
      & + R_scaled6 * (1.0_wp/fact(3) &
      & + R_scaled6 * (5.0_wp * 1699.0_wp/207360.0_wp &
      & + R_scaled6 * (6.0_wp * 259.0_wp/207360.0_wp &
      & + R_scaled6 * (7.0_wp * 197.0_wp/1451520.0_wp &
      & + R_scaled6 * (8.0_wp * 19.0_wp/2177280.0_wp &
      & + R_scaled6 * (9.0_wp * 1.0_wp/4354560.0_wp))))))))

   phi2 = 1.0_wp - exp_R * P2
   dphi2 = exp_R * (P2 - dP2)

   P0 = R_scaled6**3 * (89.0_wp/13824.0_wp &
      & + R_scaled6 * (89.0_wp/13824.0_wp &
      & + R_scaled6 * (119.0_wp/41472.0_wp &
      & + R_scaled6 * (5.0_wp/6912.0_wp &
      & + R_scaled6 * (11.0_wp/103680.0_wp &
      & + R_scaled6 * (1.0_wp/124416.0_wp &
      & + R_scaled6 * (1.0_wp/4354560.0_wp)))))))

   dP0 = R_scaled6**2 * (3.0_wp * 89.0_wp/13824.0_wp &
      & + R_scaled6 * (4.0_wp * 89.0_wp/13824.0_wp &
      & + R_scaled6 * (5.0_wp * 119.0_wp/41472.0_wp &
      & + R_scaled6 * (6.0_wp * 5.0_wp/6912.0_wp &
      & + R_scaled6 * (7.0_wp * 11.0_wp/103680.0_wp &
      & + R_scaled6 * (8.0_wp * 1.0_wp/124416.0_wp &
      & + R_scaled6 * (9.0_wp * 1.0_wp/4354560.0_wp)))))))

   phi0 = exp_R * P0
   dphi0 = exp_R * (dP0 - P0)

   chi_11 = phi2**2 + 0.5_wp * phi0**2
   dchi_11 = 2.0_wp * phi2 * dphi2 + phi0 * dphi0

   ! Scaled distance depending on the element damping radius
   R_scaled8 = param%rs8 * r1 / rpoly
   exp_R = exp(-R_scaled8)

   ! Dipole-quadrupole interaction
   P3 = 1.0_wp + R_scaled8 * (1.0_wp + R_scaled8 * (1.0_wp/fact(2) &
      & + R_scaled8 * (1.0_wp/fact(3) + R_scaled8 * (1.0_wp/fact(4) &
      & + R_scaled8 * (1.0_wp/fact(5) + R_scaled8 * (1.0_wp/fact(6) &
      & + R_scaled8 * (109.0_wp/552960.0_wp &
      & + R_scaled8 * (13.0_wp/552960.0_wp &
      & + R_scaled8 * (211.0_wp/96768000.0_wp &
      & + R_scaled8 * (57.0_wp/435456000.0_wp &
      & + R_scaled8 * (1.0_wp/290304000.0_wp)))))))))))

   dP3 = 1.0_wp + R_scaled8 * (1.0_wp + R_scaled8 * (1.0_wp/fact(2) &
      & + R_scaled8 * (1.0_wp/fact(3) + R_scaled8 * (1.0_wp/fact(4) &
      & + R_scaled8 * (1.0_wp/fact(5) &
      & + R_scaled8 * (7.0_wp * 109.0_wp/552960.0_wp &
      & + R_scaled8 * (8.0_wp * 13.0_wp/552960.0_wp &
      & + R_scaled8 * (9.0_wp * 211.0_wp/96768000.0_wp &
      & + R_scaled8 * (10.0_wp * 57.0_wp/435456000.0_wp &
      & + R_scaled8 * (11.0_wp * 1.0_wp/290304000.0_wp))))))))))

   phi3 = 1.0_wp - exp_R * P3
   dphi3 = exp_R * (P3 - dP3)

   P1 = R_scaled8**5 * (47.0_wp/552960.0_wp &
      & + R_scaled8 * (47.0_wp/552960.0_wp &
      & + R_scaled8 * (7.0_wp/184320.0_wp &
      & + R_scaled8 * (1.0_wp/103680.0_wp &
      & + R_scaled8 * (209.0_wp/145152000.0_wp &
      & + R_scaled8 * (11.0_wp/96768000.0_wp &
      & + R_scaled8 * (1.0_wp/290304000.0_wp)))))))

   dP1 = R_scaled8**4 * (5.0_wp * 47.0_wp/552960.0_wp &
      & + R_scaled8 * (6.0_wp * 47.0_wp/552960.0_wp &
      & + R_scaled8 * (7.0_wp * 7.0_wp/184320.0_wp &
      & + R_scaled8 * (8.0_wp * 1.0_wp/103680.0_wp &
      & + R_scaled8 * (9.0_wp * 209.0_wp/145152000.0_wp &
      & + R_scaled8 * (10.0_wp * 11.0_wp/96768000.0_wp &
      & + R_scaled8 * (11.0_wp * 1.0_wp/290304000.0_wp)))))))

   phi1 = exp_R * P1
   dphi1 = exp_R * (dP1 - P1)

   chi_12 = phi3**2 + (2.0_wp/3.0_wp) * phi1**2
   dchi_12 = 2.0_wp * phi3 * dphi3 + (4.0_wp/3.0_wp) * phi1 * dphi1

   ! Assemble damping factors and their derivatives
   d6 = param%s6 * chi_11 / r6
   d8 = param%s8 * chi_12 / r8

   d6dr = ( -6.0_wp * d6 + param%s6 * R_scaled6 * dchi_11 / r6 ) / r2
   d8dr = ( -8.0_wp * d8 + param%s8 * R_scaled8 * dchi_12 / r8 ) / r2

end subroutine get_2b_derivs


!> Check the availability of the parameters for two-body damping
subroutine check_2b_params(self, error, param)
   !> Instance of the two-body damping function
   class(koide_damping_twobody), intent(in) :: self
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Damping parameters
   type(param_type), intent(in) :: param

   if (.not. allocated(param%s6)) then
      call fatal_error(error, "Missing parameter in koide damping: s6")
   end if
   if (.not. allocated(param%s8)) then
      call fatal_error(error, "Missing parameter in koide damping: s8")
   end if
   if (.not. allocated(param%rs6)) then
      call fatal_error(error, "Missing parameter in koide damping: rs6")
   end if
   if (.not. allocated(param%rs8)) then
      call fatal_error(error, "Missing parameter in koide damping: rs8")
   end if

end subroutine check_2b_params

end module dftd4_damping_koide