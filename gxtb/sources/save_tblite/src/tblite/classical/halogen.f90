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

!> @dir tblite/classical
!> Contains classical correction potentials

!> @file tblite/classical/halogen.f90
!> Provides a classical halogen bonding correction

!> Classical correction term for halogen bonding contributions
module tblite_classical_halogen
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_data_atomicrad, only : get_atomic_rad
   use tblite_container, only : container_cache
   use tblite_cutoff, only : get_lattice_points
   use tblite_classical_type, only : classical_type
   implicit none
   private

   public :: new_halogen_correction


   !> Container for evaluating halogen bonding energy terms by a classical potential
   type, public, extends(classical_type) :: halogen_correction
      !> Interaction strength of the halogen bond
      real(wp), allocatable :: bond_strength(:)
      !> Atomic radii of all atoms
      real(wp), allocatable :: rad(:)
      !> Damping for the interaction potential
      real(wp) :: damping
      !> Suitable acceptors for halogen bonding interactions
      logical, allocatable :: acceptor(:)
      !> Suitable donors for halogen bonding interactions
      logical, allocatable :: halogen(:)
      !> Real-space cutoff
      real(wp) :: cutoff = 20.0_wp
   contains
      !> Entry point for evaluation of energy and gradient
      procedure :: get_engrad
   end type halogen_correction

   real(wp), parameter :: alp = 6.0_wp, lj = 12.0_wp, lj2 = lj * 0.5_wp
   character(len=*), parameter :: label = "halogen-bond correction"

contains


!> Construct new halogen bonding correction
subroutine new_halogen_correction(self, mol, damping, rad_scale, bond_strength, rad, cutoff)
   !> Instance of the halogen bond correction
   type(halogen_correction), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Damping factor
   real(wp), intent(in) :: damping
   !> Scaling of the atomic radii
   real(wp), intent(in) :: rad_scale
   !> Strength of the halogen bond
   real(wp), intent(in) :: bond_strength(:)
   !> Atomic radii for each species
   real(wp), intent(in), optional :: rad(:)
   !> Real-space cutoff
   real(wp), intent(in), optional :: cutoff

   self%label = label
   allocate(self%rad(mol%nid), self%bond_strength(mol%nid), self%halogen(mol%nid), &
      & self%acceptor(mol%nid))
   if (present(rad)) then
      self%rad(:) = rad * rad_scale
   else
      self%rad(:) = get_atomic_rad(mol%num) * rad_scale
   end if
   if (present(cutoff)) self%cutoff = cutoff
   self%bond_strength(:) = bond_strength
   self%damping = damping

   self%halogen(:) = is_halogen(mol%num)
   self%acceptor(:) = is_acceptor(mol%num)

end subroutine new_halogen_correction


!> Check whether an atom qualifies as halogen bond doner
elemental function is_halogen(num) result(halogen)
   integer, intent(in) :: num
   logical :: halogen
   halogen = any(num == [17, 35, 53, 85])
end function is_halogen

!> Check whether an atom qualifies as halogen bond acceptor
elemental function is_acceptor(num) result(acceptor)
   integer, intent(in) :: num
   logical :: acceptor
   acceptor = any(num == [7, 8, 15, 16])
end function is_acceptor


!> Evaluate classical interaction for energy and derivatives
subroutine get_engrad(self, mol, cache, energies, gradient, sigma)
   !> Instance of the halogen bond correction
   class(halogen_correction), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different runs
   type(container_cache), intent(inout) :: cache
   !> Halogen-bonding energy
   real(wp), intent(inout) :: energies(:)
   !> Molecular gradient of the halogen-bonding energy
   real(wp), contiguous, intent(inout), optional :: gradient(:, :)
   !> Strain derivatives of the halogen-bonding energy
   real(wp), contiguous, intent(inout), optional :: sigma(:, :)

   integer, allocatable :: list(:, :)
   real(wp), allocatable :: trans(:, :)

   if (count(self%halogen) * count(self%acceptor) == 0) return

   call get_lattice_points(mol%periodic, mol%lattice, self%cutoff, trans)

   call get_xbond_list(mol, trans, self%cutoff, self%halogen, self%acceptor, list)
   if (.not.allocated(list)) return

   if (present(gradient) .and. present(sigma)) then
      call get_xbond_derivs(mol, trans, list, self%damping, self%bond_strength, self%rad, &
         & energies, gradient, sigma)
   else
      call get_xbond_energy(mol, trans, list, self%damping, self%bond_strength, self%rad, &
         & energies)
   end if

end subroutine get_engrad


!> Collect all triples for halogen bonding interactions
subroutine get_xbond_list(mol, trans, cutoff, halogen, acceptor, list)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Lattice points
   real(wp), intent(in) :: trans(:, :)
   !> Real space cutoff
   real(wp), intent(in) :: cutoff
   !> Species can act as halogen bond donor
   logical, intent(in) :: halogen(:)
   !> Species can act as halogen bond acceptor
   logical, intent(in) :: acceptor(:)
   !> List of all halogen bonds
   integer, allocatable, intent(out) :: list(:, :)

   integer :: iat, jat, kat, izp, jzp, jtr, ktr, nxb, ij
   real(wp) :: vec(3), r1, dist

   nxb = 0
   call resize(list)

   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, mol%nat
         jzp = mol%id(jat)
         if (.not.(halogen(izp).and.acceptor(jzp))) cycle
         do jtr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, jtr)
            r1 = norm2(vec)
            if (r1 > cutoff) cycle
            nxb = nxb + 1
            if (nxb > size(list, 2)) call resize(list)
            list([1, 2, 4], nxb) = [iat, jat, jtr]
         end do
      end do
   end do

   if (nxb == 0) then
      deallocate(list)
      return
   end if

   call resize(list, nxb)

   do ij = 1, size(list, 2)
      iat = list(1, ij)
      jat = list(2, ij)
      dist = huge(1.0_wp)
      do kat = 1, mol%nat
         do ktr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - mol%xyz(:, kat) - trans(:, ktr)
            r1 = norm2(vec)
            if (r1 < dist .and. r1 > 0.0_wp) then
               list([3, 5], ij) = [kat, ktr]
               dist = r1
            end if
         end do
      end do
   end do
end subroutine get_xbond_list


!> Reallocate list of integers
pure subroutine resize(var, n)
   !> Instance of the array to be resized
   integer, allocatable, intent(inout) :: var(:, :)
   !> Dimension of the final array size
   integer, intent(in), optional :: n

   integer, parameter :: nitem = 5, initial_size = 20
   integer, allocatable :: tmp(:, :)
   integer :: this_size, new_size

   if (allocated(var)) then
      this_size = size(var, 2)
      call move_alloc(var, tmp)
   else
      this_size = initial_size
   end if

   if (present(n)) then
      new_size = n
   else
      new_size = this_size + this_size/2 + 1
   end if

   allocate(var(nitem, new_size))

   if (allocated(tmp)) then
      this_size = min(size(tmp, 2), size(var, 2))
      var(:, :this_size) = tmp(:, :this_size)
      deallocate(tmp)
   end if

end subroutine resize


!> Get energy contributions from halogen bonding interactions
subroutine get_xbond_energy(mol, trans, list, damping, bond_strength, rad, &
      & energies)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Lattice points
   real(wp), intent(in) :: trans(:, :)
   !> List of halogen bonded atoms
   integer, intent(in) :: list(:, :)
   !> Damping factor
   real(wp), intent(in) :: damping
   !> Strength of the halogen bond
   real(wp), intent(in) :: bond_strength(:)
   !> Atomic radii for each species
   real(wp), intent(in) :: rad(:)
   !> Repulsion energy
   real(wp), intent(inout) :: energies(:)

   integer :: ijk, jat, kat, xat, xzp, jzp, jtr, ktr
   real(wp) :: cc, r0jx, t13, t14, d2jx, rjx, term, aterm
   real(wp) :: xy, d2kx, d2jk, dxj(3), dxk(3), dkj(3)

   do ijk = 1, size(list, 2)
      xat = list(1, ijk)
      jat = list(2, ijk)
      kat = list(3, ijk)
      jtr = list(4, ijk)
      ktr = list(5, ijk)
      xzp = mol%id(xat)
      jzp = mol%id(jat)
      cc = bond_strength(xzp)
      r0jx = rad(xzp) + rad(jzp)
      dxj = mol%xyz(:, jat) - mol%xyz(:, xat)
      dxk = mol%xyz(:, kat) + trans(:, ktr) - mol%xyz(:, xat)
      dkj = mol%xyz(:, jat) + trans(:, jtr) - mol%xyz(:, kat) - trans(:, ktr)
      d2jx = sum(dxj*dxj)
      d2kx = sum(dxk*dxk)
      d2jk = sum(dkj*dkj)
      rjx = sqrt(d2jx)
      ! angle part. term = cos angle kat-xat-jat
      xy = sqrt(d2kx*d2jx)
      term = (d2kx+d2jx-d2jk) / xy
      aterm = (0.5_wp-0.25_wp*term)**alp
      t13 = r0jx/rjx
      t14 = t13**lj
      energies(xat) = energies(xat) + aterm*cc*(t14-damping*t13**lj2) / (1.0_wp+t14)
   end do

end subroutine get_xbond_energy


!> Get energy and its derivatives for halogen bonding interactions
subroutine get_xbond_derivs(mol, trans, list, damping, bond_strength, rad, &
      & energies, gradient, sigma)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Lattice points
   real(wp), intent(in) :: trans(:, :)
   !> Damping factor
   real(wp), intent(in) :: damping
   !> List of halogen bonded atoms
   integer, intent(in) :: list(:, :)
   !> Strength of the halogen bond
   real(wp), intent(in) :: bond_strength(:)
   !> Atomic radii for each species
   real(wp), intent(in) :: rad(:)
   !> Repulsion energy
   real(wp), intent(inout) :: energies(:)
   !> Molecular gradient of the repulsion energy
   real(wp), intent(inout) :: gradient(:, :)
   !> Strain derivatives of the repulsion energy
   real(wp), intent(inout) :: sigma(:, :)

   integer :: ijk, jat, kat, xat, xzp, jzp, jtr, ktr
   real(wp) :: cc, r0jx, t13, t14
   real(wp) :: d2jx, rjx, rkx, term, aterm, xy, d2kx, d2jk
   real(wp) :: dxj(3), dxk(3), dkj(3), dcosterm
   real(wp) :: dtermlj, termlj, prefactor, numerator, denominator

   do ijk = 1, size(list, 2)
      xat = list(1, ijk)
      jat = list(2, ijk)
      kat = list(3, ijk)
      jtr = list(4, ijk)
      ktr = list(5, ijk)
      xzp = mol%id(xat)
      jzp = mol%id(jat)
      cc = bond_strength(xzp)
      r0jx = rad(xzp)+rad(jzp)

      dxj = mol%xyz(:, jat) - mol%xyz(:, xat)
      dxk = mol%xyz(:, kat) + trans(:, ktr) - mol%xyz(:, xat)
      dkj = mol%xyz(:, jat) + trans(:, jtr) - mol%xyz(:, kat) - trans(:, ktr)

      d2jx = sum(dxj*dxj)
      d2kx = sum(dxk*dxk)
      d2jk = sum(dkj*dkj)
      rjx = sqrt(d2jx)+1.0e-18_wp
      rkx = sqrt(d2kx)+1.0e-18_wp

      xy = sqrt(d2kx*d2jx)
      term = (d2kx+d2jx-d2jk) / xy
      ! now compute angular damping function
      aterm = (0.5_wp-0.25_wp*term)**alp

      t13 = r0jx/rjx
      t14 = t13**lj
      energies(xat) = energies(xat) + aterm*cc*(t14-damping*t13**lj2) / (1.0_wp+t14)

      ! set up weighted inverted distance and compute the modified Lennard-Jones potential
      t14 = (r0jx/rjx)**lj2 ! (rov/r)^lj2 ; lj2 = 6 in GFN1
      numerator = (t14*t14 - damping*t14)
      denominator = (1.0_wp + t14*t14)
      termLJ = numerator/denominator

      ! LJ derivative
      ! denominator part
      dtermlj = 2.0_wp*lj2*numerator*t14*t14/(rjx*denominator*denominator)
      ! numerator part
      dtermlj = dtermlj+lj2*t14*(damping - 2.0_wp*t14)/(rjx*denominator)
      ! scale w/ angular damping term
      dtermlj = dtermlj*aterm*cc/rjx
      ! gradient for the acceptor
      gradient(:, jat) = gradient(:, jat)+dtermlj*dxj(:)
      ! halogen gradient
      gradient(:, xat) = gradient(:, xat)-dtermlj*dxj(:)
      ! strain contribution
      sigma(:, :) = sigma + spread(dxj, 1, 3) * spread(dxj, 2, 3) * dtermlj

      ! cosine term derivative
      prefactor = -0.250_wp*alp*(0.5_wp-0.25_wp*term)**(alp-1.0_wp)
      prefactor = prefactor*cc*termlj
      ! AX part
      dcosterm = 2.0_wp/rkx - term/rjx
      dcosterm = dcosterm*prefactor/rjx
      ! gradient for the acceptor
      gradient(:, jat) = gradient(:, jat)+dcosterm*dxj(:)
      ! halogen gradient
      gradient(:, xat) = gradient(:, xat)-dcosterm*dxj(:)
      ! strain contribution
      sigma(:, :) = sigma + spread(dxj, 1, 3) * spread(dxj, 2, 3) * dcosterm
      ! KX part
      dcosterm = 2.0_wp/rjx - term/rkx
      dcosterm = dcosterm*prefactor/rkx
      ! gradient for the acceptor
      gradient(:, kat) = gradient(:, kat)+dcosterm*dxk(:)
      ! halogen gradient
      gradient(:, xat) = gradient(:, xat)-dcosterm*dxk(:)
      ! strain contribution
      sigma(:, :) = sigma + spread(dxk, 1, 3) * spread(dxk, 2, 3) * dcosterm
      ! JK part
      t13 = 2.0_wp*prefactor/xy
      ! acceptor
      gradient(:, jat) = gradient(:, jat)-t13*dkj(:)
      ! neighbor
      gradient(:, kat) = gradient(:, kat)+t13*dkj(:)
      ! strain contribution
      sigma(:, :) = sigma - spread(dkj, 1, 3) * spread(dkj, 2, 3) * t13
   end do

end subroutine get_xbond_derivs


end module tblite_classical_halogen
