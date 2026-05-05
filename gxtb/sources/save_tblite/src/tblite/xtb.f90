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

!> @dir tblite/xtb
!> Routines related to the extend tight-binding Hamiltonian

!> @file tblite/xtb.f90
!> Proxy module for extended tight-binding Hamiltonian related settings

!> This module contains reexports for the #tblite_xtb_calculator::xtb_calculator class
!> and the possible parametrizations via #tblite_xtb_gfn1, #tblite_xtb_gfn2, 
!> #tblite_xtb_gxtb, and #tblite_xtb_ipea1. The main entry point for performing 
!> calculations is provided with #tblite_xtb_singlepoint::xtb_singlepoint.
module tblite_xtb
   use tblite_xtb_calculator, only : xtb_calculator, new_xtb_calculator, param_h0spec
   use tblite_xtb_gfn1, only : new_gfn1_calculator, gfn1_h0spec, export_gfn1_param
   use tblite_xtb_gfn2, only : new_gfn2_calculator, gfn2_h0spec, export_gfn2_param
   use tblite_xtb_gxtb, only : new_gxtb_calculator, gxtb_h0spec, export_gxtb_param
   use tblite_xtb_ipea1, only : new_ipea1_calculator, ipea1_h0spec, export_ipea1_param
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   implicit none

   public :: xtb_calculator, new_xtb_calculator, param_h0spec
   public :: new_gfn1_calculator, gfn1_h0spec, export_gfn1_param
   public :: new_gfn2_calculator, gfn2_h0spec, export_gfn2_param
   public :: new_gxtb_calculator, gxtb_h0spec, export_gxtb_param
   public :: new_ipea1_calculator, ipea1_h0spec, export_ipea1_param
   public :: xtb_singlepoint
end module tblite_xtb
