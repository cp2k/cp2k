! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later
!
! tblite is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! tblite is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with tblite.  If not, see <https://www.gnu.org/licenses/>.

module tblite_features
   implicit none
   private

   public :: get_tblite_feature


contains

function get_tblite_feature(feature) result(has_feature)
   !> Feature name
   character(len=*), intent(in) :: feature
   !> Whether the feature is enabled
   logical :: has_feature

   select case(feature)
   case("color")
      has_feature = color_support()
   case default
      has_feature = .false.
   end select
end function get_tblite_feature

!> Check output terminal for color support
function color_support() result(color)
   use, intrinsic :: iso_fortran_env, only : output_unit
#if defined __INTEL_COMPILER
   use ifport, only : isatty
#endif
   !> Whether color can be supported
   logical :: color

   color = .false.
#if defined __GFORTRAN__ || defined __INTEL_COMPILER
   color = isatty(output_unit)
#endif
end function color_support

end module tblite_features
