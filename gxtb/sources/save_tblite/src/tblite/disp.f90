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

!> @dir tblite/disp
!> Contains implemenations for the evaluation of London-dispersion interactions.

!> @file tblite/disp.f90
!> Reexports of the available dispersion interactions.

!> Proxy module for dispersion interactions.
module tblite_disp
   use mctc_env, only : wp
   use tblite_disp_cache, only : dispersion_cache
   use tblite_disp_d3, only : d3_dispersion, new_d3_dispersion
   use tblite_disp_d4, only : d4_dispersion, new_d4_dispersion, &
      & new_d4s_dispersion, new_d4srev_dispersion, twobody_damping_function, &
      & threebody_damping_function, get_damping_function_id
   use tblite_disp_type, only : dispersion_type
   implicit none
   public

end module tblite_disp
