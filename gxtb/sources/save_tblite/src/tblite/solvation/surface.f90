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

!> @file tblite/solvation/surface.f90
!> Provides a surface integrator

!> Surface integrator for solvent accessible surface area
module tblite_solvation_surface
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use mctc_io_convert, only : aatoau
   use tblite_mesh_lebedev, only : get_angular_grid, grid_size, list_bisection
   use tblite_adjlist, only : adjacency_list, new_adjacency_list
   implicit none
   private

   public :: new_surface_integrator

   type, public :: surface_integrator
      !> Number of angular grid points
      integer :: nang
      !> Angular grid coordinates
      real(wp), allocatable :: ang_grid(:, :)
      !> Angular grid weights
      real(wp), allocatable :: ang_weight(:)
      !> Cut-off radius for the nearest neighbour list
      real(wp) :: srcut
      !> Atom specific surface data
      real(wp), allocatable :: vdwsa(:)
      real(wp), allocatable :: wrp(:)
      real(wp), allocatable :: trj2(:, :)
      real(wp) :: ah0, ah1, ah3
   contains
      procedure :: get_surface
   end type surface_integrator

   !> real space cut-offs
   real(wp), parameter :: tolsesp = 1.e-6_wp

   real(wp), parameter :: default_offset = 2.0_wp*aatoau
   real(wp), parameter :: default_smoothing = 0.3_wp*aatoau

contains

!> Initialize data straucture
subroutine new_surface_integrator(self, num, rad, probe, nang, offset, smoothing)
   !> Instance of the surface integrator
   type(surface_integrator), intent(out) :: self
   !> Atomic numbers
   integer, intent(in) :: num(:)
   !> Van-der-Waals Radii
   real(wp), intent(in) :: rad(:)
   !> Probe radius of the solvent
   real(wp), intent(in) :: probe
   !> Number of angular grid points for integration
   integer, intent(in) :: nang
   !> Offset for surface integration cutoff
   real(wp), intent(in), optional :: offset
   !> Smooting function parameter
   real(wp), intent(in), optional :: smoothing

   integer :: iat, izp, iang, ierr, nat
   real(wp) :: r, w, w3

   nat = size(num)

   allocate (self%vdwsa(nat))
   allocate (self%trj2(2, nat))
   allocate (self%wrp(nat))
   if (present(smoothing)) then
      w = smoothing
   else
      w = default_smoothing
   end if
   w3 = w*w*w
   self%ah0 = 0.5_wp
   self%ah1 = 3._wp/(4.0_wp*w)
   self%ah3 = -1._wp/(4.0_wp*w3)
   do iat = 1, nat
      izp = num(iat)
      self%vdwsa(iat) = rad(izp) + probe
      self%trj2(1, iat) = (self%vdwsa(iat) - w)**2
      self%trj2(2, iat) = (self%vdwsa(iat) + w)**2
      r = self%vdwsa(iat) + w
      self%wrp(iat) = (0.25_wp/w + &
         &            3.0_wp*self%ah3*(0.2_wp*r*r - 0.5_wp*r*self%vdwsa(iat) + &
         &            self%vdwsa(iat)*self%vdwsa(iat)/3.0_wp))*r*r*r
      r = self%vdwsa(iat) - w
      self%wrp(iat) = self%wrp(iat) - (0.25/w + &
         &    3.0_wp*self%ah3*(0.2_wp*r*r - 0.5_wp*r*self%vdwsa(iat) + &
         &            self%vdwsa(iat)*self%vdwsa(iat)/3.0_wp))*r*r*r
   end do
   self%srcut = 2*(w + maxval(self%vdwsa))
   if (present(offset)) then
      self%srcut = self%srcut + offset
   else
      self%srcut = self%srcut + default_offset
   end if

   iang = list_bisection(grid_size, nang)
   allocate (self%ang_grid(3, grid_size(iang)))
   allocate (self%ang_weight(grid_size(iang)))
   call get_angular_grid(iang, self%ang_grid, self%ang_weight, ierr)
   self%ang_weight(:) = self%ang_weight*4*pi

end subroutine new_surface_integrator

subroutine get_surface(self, mol, surface, dsdr)
   !> Instance of the surface integrator
   class(surface_integrator), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Surface area for each sphere
   real(wp), intent(out) :: surface(:)
   !> Derivative of surface area w.r.t. to coordinate displacements
   real(wp), intent(out), optional :: dsdr(:, :, :)

   type(adjacency_list) :: list
   real(wp), parameter :: trans(3, 1) = 0.0_wp

   call new_adjacency_list(list, mol, trans, self%srcut, .true.)

   ! compute solvent accessible surface and its derivatives
   if (present(dsdr)) then
      call compute_numsa(mol%nat, mol%xyz, list, self%vdwsa, &
         & self%wrp, self%trj2, self%ah0, self%ah1, self%ah3, self%ang_weight, self%ang_grid, &
         & surface, dsdr)
   else
      call compute_surface(mol%nat, mol%xyz, list, self%vdwsa, &
         & self%wrp, self%trj2, self%ah0, self%ah1, self%ah3, self%ang_weight, self%ang_grid, &
         & surface)
   end if

end subroutine get_surface

subroutine compute_surface(nat, xyz, list, vdwsa, &
      & wrp, trj2, ah0, ah1, ah3, ang_weight, ang_grid, surface)
   !> Number of atoms
   integer, intent(in) :: nat
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)
   !> Neighbourlist
   type(adjacency_list), intent(in) :: list
   !> Van-der-Waals radii including probe radius of solvent
   real(wp), intent(in) :: vdwsa(:)
   !> Radial weights for the numerical integration
   real(wp), intent(in) :: wrp(:)
   !> Radial smoothing function
   real(wp), intent(in) :: trj2(:, :)
   real(wp), intent(in)  :: ah0, ah1, ah3
   !> Angular weights for the numerical integration
   real(wp), intent(in) :: ang_weight(:)
   !> Angular grid for each atom
   real(wp), intent(in) :: ang_grid(:, :)
   !> Surface area for each atom, including surface tension
   real(wp), intent(out) :: surface(:)

   integer :: iat, ip, nno, ino
   real(wp) :: rsas, sasai, xyza(3), xyzp(3), sasap, wr, wsa

   surface(:) = 0.0_wp

   !$omp parallel do default(none) shared(surface, ah0, ah1, ah3) &
   !$omp shared(nat, vdwsa, list, xyz, wrp, ang_grid, ang_weight, trj2) &
   !$omp private(iat, rsas, nno, sasai, xyza, wr, ip, xyzp, wsa, sasap, ino)
   do iat = 1, nat

      rsas = vdwsa(iat)
      ino = list%inl(iat)
      nno = list%nnl(iat)

      ! initialize storage
      sasai = 0.0_wp

      ! atomic position
      xyza(:) = xyz(:, iat)
      ! radial atomic weight
      wr = wrp(iat)

      ! loop over grid points
      do ip = 1, size(ang_grid, 2)
         ! grid point position
         xyzp(:) = xyza(:) + rsas*ang_grid(:, ip)
         ! atomic surface function at the grid point
         call compute_w_sp(nat, list%nlat(ino+1:ino+nno), trj2, vdwsa, xyz, xyzp, &
            & ah0, ah1, ah3, sasap)

         if (sasap > tolsesp) then
            ! numerical quadrature weight
            wsa = ang_weight(ip)*wr*sasap
            ! accumulate the surface area
            sasai = sasai + wsa
         end if
      end do

      surface(iat) = sasai
   end do

end subroutine compute_surface

subroutine compute_numsa(nat, xyz, list, vdwsa, &
      & wrp, trj2, ah0, ah1, ah3, ang_weight, ang_grid, surface, dsdrt)
   !> Number of atoms
   integer, intent(in) :: nat
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)
   !> Neighbourlist
   type(adjacency_list), intent(in) :: list
   !> Van-der-Waals radii including probe radius of solvent
   real(wp), intent(in) :: vdwsa(:)
   !> Radial weights for the numerical integration
   real(wp), intent(in) :: wrp(:)
   !> Radial smoothing function
   real(wp), intent(in) :: trj2(:, :)
   real(wp), intent(in)  :: ah0, ah1, ah3
   !> Angular weights for the numerical integration
   real(wp), intent(in) :: ang_weight(:)
   !> Angular grid for each atom
   real(wp), intent(in) :: ang_grid(:, :)
   !> Surface area for each atom, including surface tension
   real(wp), intent(out) :: surface(:)
   !> Derivative of surface area w.r.t. cartesian coordinates
   real(wp), intent(out) :: dsdrt(:, :, :)

   integer :: iat, jat, ip, jj, nni, nno, ino
   real(wp) :: rsas, sasai, xyza(3), xyzp(3), sasap, wr, wsa, drjj(3)
   real(wp), allocatable :: grds(:, :), grads(:, :)
   integer, allocatable :: grdi(:)

   surface(:) = 0.0_wp
   dsdrt(:, :, :) = 0.0_wp

   ! allocate space for the gradient storage
   allocate (grads(3, nat), source=0.0_wp)
   allocate (grds(3, maxval(list%nnl)))
   allocate (grdi(maxval(list%nnl)))

   !$omp parallel do default(none) shared(surface, dsdrt, ah0, ah1, ah3) &
   !$omp shared(nat, vdwsa, list, xyz, wrp, ang_grid, ang_weight, trj2) &
   !$omp private(iat, jat, rsas, nno, grads, sasai, xyza, wr, ip, xyzp, wsa, &
   !$omp& sasap, jj, nni, grdi, grds, drjj, ino)
   do iat = 1, nat

      rsas = vdwsa(iat)
      ino = list%inl(iat)
      nno = list%nnl(iat)

      ! initialize storage
      grads = 0.0_wp
      sasai = 0.0_wp

      ! atomic position
      xyza(:) = xyz(:, iat)
      ! radial atomic weight
      wr = wrp(iat)

      ! loop over grid points
      do ip = 1, size(ang_grid, 2)
         ! grid point position
         xyzp(:) = xyza(:) + rsas*ang_grid(:, ip)
         ! atomic surface function at the grid point
         call compute_w_spg(nat, list%nlat(ino+1:ino+nno), trj2, vdwsa, xyz, xyzp, &
            & ah0, ah1, ah3, sasap, grds, nni, grdi)

         if (sasap > tolsesp) then
            ! numerical quadrature weight
            wsa = ang_weight(ip)*wr*sasap
            ! accumulate the surface area
            sasai = sasai + wsa
            ! accumulate the surface gradient
            do jj = 1, nni
               jat = grdi(jj)
               drjj(:) = wsa*grds(:, jj)
               grads(:, iat) = grads(:, iat) + drjj(:)
               grads(:, jat) = grads(:, jat) - drjj(:)
            end do
         end if
      end do

      surface(iat) = sasai
      dsdrt(:, :, iat) = grads
   end do

end subroutine compute_numsa

pure subroutine compute_w_sp(nat, nnlists, trj2, vdwsa, xyza, xyzp, ah0, ah1, ah3, &
      & sasap)
   integer, intent(in)  :: nat
   integer, intent(in)  :: nnlists(:)
   real(wp), intent(in)  :: xyza(:, :)
   real(wp), intent(in)  :: xyzp(:)
   real(wp), intent(in)  :: ah0, ah1, ah3
   real(wp), intent(out) :: sasap
   real(wp), intent(in)  :: trj2(:, :)
   real(wp), intent(in)  :: vdwsa(:)

   integer  :: img, jat
   real(wp) :: tj(3), tj2, sqtj
   real(wp) :: uj, ah3uj2
   real(wp) :: sasaij, dsasaij

   sasap = 1.0_wp
   do img = 1, size(nnlists)
      jat = nnlists(img)
      ! compute the distance to the atom
      tj(:) = xyzp(:) - xyza(:, jat)
      tj2 = dot_product(tj, tj)
      ! if within the outer cut-off compute
      if (tj2 < trj2(2, jat)) then
         if (tj2 <= trj2(1, jat)) then
            sasap = 0.0_wp
            return
         else
            sqtj = sqrt(tj2)
            uj = sqtj - vdwsa(jat)
            ah3uj2 = ah3*uj*uj
            dsasaij = ah1 + 3.0_wp*ah3uj2
            sasaij = ah0 + (ah1 + ah3uj2)*uj

            ! accumulate the molecular surface
            sasap = sasap*sasaij
         end if
      end if
   end do

end subroutine compute_w_sp

pure subroutine compute_w_spg(nat, nnlists, trj2, vdwsa, xyza, xyzp, ah0, ah1, ah3, &
      & sasap, grds, nni, grdi)
   integer, intent(in)  :: nat
   integer, intent(in)  :: nnlists(:)
   integer, intent(out) :: nni
   real(wp), intent(in)  :: xyza(:, :)
   real(wp), intent(in)  :: xyzp(:)
   real(wp), intent(in)  :: ah0, ah1, ah3
   real(wp), intent(out) :: sasap
   real(wp), intent(out) :: grds(:, :)
   integer, intent(out) :: grdi(:)
   real(wp), intent(in)  :: trj2(:, :)
   real(wp), intent(in)  :: vdwsa(:)

   integer  :: img, jat
   real(wp) :: tj(3), tj2, sqtj
   real(wp) :: uj, ah3uj2
   real(wp) :: sasaij, dsasaij

   nni = 0
   sasap = 1.0_wp
   do img = 1, size(nnlists)
      jat = nnlists(img)
      ! compute the distance to the atom
      tj(:) = xyzp(:) - xyza(:, jat)
      tj2 = dot_product(tj, tj)
      ! if within the outer cut-off compute
      if (tj2 < trj2(2, jat)) then
         if (tj2 <= trj2(1, jat)) then
            sasap = 0.0_wp
            return
         else
            sqtj = sqrt(tj2)
            uj = sqtj - vdwsa(jat)
            ah3uj2 = ah3*uj*uj
            dsasaij = ah1 + 3.0_wp*ah3uj2
            sasaij = ah0 + (ah1 + ah3uj2)*uj

            ! accumulate the molecular surface
            sasap = sasap*sasaij
            ! compute the gradient wrt the neighbor
            dsasaij = dsasaij/(sasaij*sqtj)
            nni = nni + 1
            grdi(nni) = jat
            grds(:, nni) = dsasaij*tj(:)
         end if
      end if
   end do

end subroutine compute_w_spg

end module tblite_solvation_surface
