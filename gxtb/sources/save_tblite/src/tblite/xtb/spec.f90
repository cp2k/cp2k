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

!> @file tblite/xtb/spec.f90
!> Provides generator base class for Hamiltonian parameters.

!> Definition of an abstract base class of the generator for the Hamiltonian.
module tblite_xtb_spec
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_data_atomicrad, only : get_atomic_rad
   use mctc_data_vdwrad, only : get_vdw_rad
   use tblite_basis_type, only : basis_type
   implicit none
   private


   !> Specification of the effective Hamiltonian.
   !>
   !> An instance of the specification is consumed by the constructor of the
   !> Hamiltonian to generate the respective entries in the derived type.
   type, public, abstract :: tb_h0spec
   contains
      !> Generator for the self energy / atomic levels of the Hamiltonian
      procedure(get_selfenergy), deferred :: get_selfenergy
      !> Generator for the coordination number dependent shift of the self energy
      procedure :: get_cnshift
      !> Generator for the electronegativity scaled coordination number dependent shift of the self energy
      procedure :: get_cnenshift
      !> Generator for the linear partial charge dependent shift of the self energy
      procedure :: get_q1shift
      !> Generator for the quadratic partial charge dependent shift of the self energy
      procedure :: get_q2shift
      !> Generator for the enhancement factor to for scaling Hamiltonian elements
      procedure(get_hscale), deferred :: get_hscale
      !> Generator for the atomic radii used in the distant dependent scaling
      procedure :: get_rad
      !> Generator for the polynomial parameters for the square-root distant dependent scaling
      procedure :: get_shpoly
      !> Generator for the polynomial parameters for the linear distant dependent scaling
      procedure :: get_shpoly2
      !> Generator for the polynomial parameters for the square distant dependent scaling
      procedure :: get_shpoly4
      !> Generator for the coefficients of the anisotropic Hamiltonian
      procedure :: get_anisotropy
      !> Generator for the scaled van der Waals radii
      procedure :: get_rvdw
      !> Generator for the reference occupation numbers of the atoms
      procedure(get_reference_occ), deferred :: get_reference_occ
      !> Generator for the diatomic frame scaling factors
      procedure :: get_diat_scale
   end type tb_h0spec

   abstract interface
      !> Generator for the enhancement factor to for scaling Hamiltonian elements
      subroutine get_hscale(self, mol, bas, hscale)
         import :: tb_h0spec, structure_type, basis_type, wp
         !> Instance of the Hamiltonian specification
         class(tb_h0spec), intent(in) :: self
         !> Molecular structure data
         type(structure_type), intent(in) :: mol
         !> Basis set information
         class(basis_type), intent(in) :: bas
         !> Scaling parameters for the Hamiltonian elements
         real(wp), intent(out) :: hscale(:, :, :, :)
      end subroutine get_hscale

      !> Generator for the self energy / atomic levels of the Hamiltonian
      subroutine get_selfenergy(self, mol, bas, selfenergy)
         import :: tb_h0spec, structure_type, basis_type, wp
         !> Instance of the Hamiltonian specification
         class(tb_h0spec), intent(in) :: self
         !> Molecular structure data
         type(structure_type), intent(in) :: mol
         !> Basis set information
         class(basis_type), intent(in) :: bas
         !> Self energy / atomic levels
         real(wp), intent(out) :: selfenergy(:, :)
      end subroutine get_selfenergy

      !> Generator for the reference occupation numbers of the atoms
      subroutine get_reference_occ(self, mol, bas, refocc)
         import :: tb_h0spec, structure_type, basis_type, wp
         !> Instance of the Hamiltonian specification
         class(tb_h0spec), intent(in) :: self
         !> Molecular structure data
         type(structure_type), intent(in) :: mol
         !> Basis set information
         class(basis_type), intent(in) :: bas
         !> Reference occupation numbers
         real(wp), intent(out) :: refocc(:, :)
      end subroutine get_reference_occ
   end interface

contains


!> Stub implementation of the coordination number dependent shift generator
subroutine get_cnshift(self, mol, bas, kcn)
   !> Instance of the Hamiltonian specification
   class(tb_h0spec), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Coordination number dependent shift
   real(wp), intent(out) :: kcn(:, :)

   kcn(:, :) = 0.0_wp
end subroutine get_cnshift


!> Stub implementation of the electronegativity scaled coordination number dependent shift generator
subroutine get_cnenshift(self, mol, bas, kcn_en)
   !> Instance of the Hamiltonian specification
   class(tb_h0spec), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Coordination number dependent shift
   real(wp), intent(out) :: kcn_en(:, :)

   kcn_en(:, :) = 0.0_wp
end subroutine get_cnenshift


!> Stub implementation of the linear partial charge dependent shift generator
subroutine get_q1shift(self, mol, bas, kq1)
   !> Instance of the Hamiltonian specification
   class(tb_h0spec), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Linear partial charge dependent shift
   real(wp), intent(out) :: kq1(:, :)

   kq1(:, :) = 0.0_wp
end subroutine get_q1shift


!> Stub implementation of the quadratic partial charge dependent shift generator
subroutine get_q2shift(self, mol, bas, kq2)
   !> Instance of the Hamiltonian specification
   class(tb_h0spec), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Quadratic partial charge dependent shift
   real(wp), intent(out) :: kq2(:, :)

   kq2(:, :) = 0.0_wp
end subroutine get_q2shift


!> Stub implementation of the atomic radii generator
subroutine get_rad(self, mol, bas, rad)
   !> Instance of the Hamiltonian specification
   class(tb_h0spec), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Atomic radii
   real(wp), intent(out) :: rad(:)

   rad(:) = get_atomic_rad(mol%num)
end subroutine get_rad


!> Stub implementation of the shell polynomial for the square-root distant dependent scaling
subroutine get_shpoly(self, mol, bas, shpoly)
   !> Instance of the Hamiltonian specification
   class(tb_h0spec), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Parameters for square-root distance-dependent scaling
   real(wp), intent(out) :: shpoly(:, :)

   shpoly(:, :) = 0.0_wp
end subroutine get_shpoly


!> Stub implementation of the shell polynomial for the linear distant dependent scaling
subroutine get_shpoly2(self, mol, bas, shpoly2)
   !> Instance of the Hamiltonian specification
   class(tb_h0spec), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Parameters for linear distance-dependent scaling
   real(wp), intent(out) :: shpoly2(:, :)

   shpoly2(:, :) = 0.0_wp
end subroutine get_shpoly2


!> Stub implementation of the shell polynomial for the square distant dependent scaling
subroutine get_shpoly4(self, mol, bas, shpoly4)
   !> Instance of the Hamiltonian specification
   class(tb_h0spec), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Parameters for square distance-dependent scaling
   real(wp), intent(out) :: shpoly4(:, :)

   shpoly4(:, :) = 0.0_wp
end subroutine get_shpoly4


!> Stub implementation of the coefficients for the anisotropic Hamiltonian
subroutine get_anisotropy(self, mol, bas, dip_scale, aniso_exp)
   !> Instance of the Hamiltonian specification
   class(tb_h0spec), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Parameters for dipole Hamiltonian correction
   real(wp), intent(out) :: dip_scale(:, :)
   !> Exponent for the anisotropic counting function
   real(wp), intent(out) :: aniso_exp

   dip_scale(:, :) = 0.0_wp
   aniso_exp = 0.0_wp
end subroutine get_anisotropy


!> Stub implementation of scaled van der Waals radii generator
subroutine get_rvdw(self, mol, bas, rvdw)
   !> Instance of the Hamiltonian specification
   class(tb_h0spec), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Pairwise scaled van der Waals radii
   real(wp), intent(out) :: rvdw(:, :)

   integer :: isp, izp, jsp, jzp

   do isp = 1, mol%nid
      izp = mol%num(isp)
      do jsp = 1, isp
         jzp = mol%num(jsp)
         rvdw(jsp, isp) = get_vdw_rad(jzp, izp) 
         rvdw(isp, jsp) = rvdw(jsp, isp)
      end do
   end do
end subroutine get_rvdw


!> Stub implementation of the diatomic frame scaling factor generator
subroutine get_diat_scale(self, mol, bas, ksig, kpi, kdel, do_diat_scale)
   !> Instance of the Hamiltonian specification
   class(tb_h0spec), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Diatomic frame scaling of sigma bonding contribution
   real(wp), intent(out) :: ksig(:, :)
   !> Diatomic frame scaling of pi bonding contribution
   real(wp), intent(out) :: kpi(:, :)
   !> Diatomic frame scaling of delta bonding contribution
   real(wp), intent(out) :: kdel(:, :)
   !> Flag indicating if the diatomic frame scaling is applied
   logical, intent(out) :: do_diat_scale

   ksig(:, :) = 0.0_wp
   kpi(:, :) = 0.0_wp
   kdel(:, :) = 0.0_wp
   do_diat_scale = .false.
end subroutine get_diat_scale


end module tblite_xtb_spec
