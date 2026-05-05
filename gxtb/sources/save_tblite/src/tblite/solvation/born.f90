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

!> @file tblite/solvation/born.f90
!> Provides a Born radii integrator

!> Integrator for Born radii based on the Onufriev-Bashford-Case model
module tblite_solvation_born
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use mctc_io_convert, only : aatoau
   use tblite_adjlist, only : adjacency_list, new_adjacency_list
   implicit none
   private

   public :: new_born_integrator

   !> Implementation of GBOBC integrator
   type, public :: born_integrator
      !> van der Waals radii of the particles
      real(wp), allocatable :: vdwr(:)
      !> pair descreening approximation radii
      real(wp), allocatable :: rho(:)
      !> offset van der Waals radii
      real(wp), allocatable :: svdw(:)
      !> cut-off radius for the Born radius NN list
      real(wp) :: lrcut
      !> Scaling factor for Born radii
      real(wp) :: born_scale
      !> Volume polynome correction, default parameters correspond to GBOBCII
      real(wp) :: obc(3)
   contains
      !> Calculate Born radii for a given geometry
      procedure :: get_rad
   end type born_integrator

   real(wp), parameter :: lrcut_default = 35.0_wp * aatoau
   real(wp), parameter :: born_scale_default = 1.0_wp
   real(wp), parameter :: born_offset_default = 0.0_wp
   real(wp), parameter :: descreening_default = 0.8_wp
   real(wp), parameter :: obc_default(3) = [1.0_wp, 0.8_wp, 4.85_wp]

contains

!> Create new Born radii integrator
subroutine new_born_integrator(self, mol, vdwrad, descreening, born_scale, born_offset, &
      & obc, rcutoff)
   !> Instance of the Born integrator
   type(born_integrator), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Van-der-Waals Radii
   real(wp), intent(in) :: vdwRad(:)
   !> Dielectric descreening parameter
   real(wp), intent(in), optional :: descreening(:)
   !> Scaling factor for Born radii
   real(wp), intent(in), optional :: born_scale
   !> Offset parameter for Born radii integration
   real(wp), intent(in), optional :: born_offset
   !> GBOBC integrator parameters
   real(wp), intent(in), optional :: obc(3)
   !> Real-space cutoff for Born radii integration
   real(wp), intent(in), optional :: rCutoff

   self%lrcut = lrcut_default
   if (present(rCutoff)) then
      self%lrcut = rCutoff
   end if

   self%born_scale = born_scale_default
   if (present(born_scale)) then
      self%born_scale = born_scale
   end if

   self%obc = obc_default
   if (present(obc)) then
      self%obc = obc
   end if

   self%vdwr = vdwRad(mol%id)

   if (present(descreening)) then
      self%rho = self%vdwr * descreening(mol%id)
   else
      self%rho = self%vdwr * descreening_default
   end if

   if (present(born_offset)) then
      self%svdw = self%vdwr - born_offset
   else
      self%svdw = self%vdwr - born_offset_default
   end if
end subroutine new_born_integrator

!> Calculate Born radii
subroutine get_rad(self, mol, rad, draddr)
   !> Instance of the Born integrator
   class(born_integrator), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Born radii
   real(wp), intent(out) :: rad(:)
   !> Derivative of Born radii w.r.t. cartesian displacements
   real(wp), intent(out), optional :: draddr(:, :, :)

   type(adjacency_list) :: list
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp), allocatable :: brdr(:, :, :)

   call new_adjacency_list(list, mol, trans, self%lrcut)
   allocate(brdr(3, mol%nat, mol%nat))

   call compute_bornr(mol%nat, mol%xyz, list, &
      & self%vdwr, self%rho, self%svdw, self%born_scale, self%obc, rad, brdr)

   if (present(draddr)) then
      draddr(:, :, :) = brdr
   end if
end subroutine get_rad

subroutine compute_bornr(nat, xyz, list, vdwr, rho, svdw, c1, obc, &
      & brad, brdr)
   !> Number of atoms
   integer, intent(in) :: nat
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)
   !> Neighbourlist
   type(adjacency_list), intent(in) :: list
   !> Van-der-Waals radii
   real(wp), intent(in) :: vdwr(:)
   !> Descreened van-der-Waals radii
   real(wp), intent(in) :: rho(:)
   !> van-der-Waals radii with offset
   real(wp), intent(in) :: svdw(:)
   !> Scaling factor for the Born radii
   real(wp), intent(in) :: c1
   !> Volume polynome correction
   real(wp), intent(in) :: obc(3)
   !> Born radii
   real(wp), intent(out) :: brad(:)

   !> Derivative of Born radii w.r.t. cartesian coordinates
   real(wp), intent(out) :: brdr(:, :, :)

   integer :: iat
   real(wp) :: br, dpsi, svdwi, vdwri, s1, v1, s2, arg, arg2
   real(wp) :: th, ch

   call compute_psi(nat, xyz, list, vdwr, rho, brad, brdr)

   do iat = 1, nat

      br = brad(iat)

      svdwi = svdw(iat)
      vdwri = vdwr(iat)
      s1 = 1.0_wp/svdwi
      v1 = 1.0_wp/vdwri
      s2 = 0.5_wp*svdwi

      br = br*s2

      arg2 = br*(obc(3)*br-obc(2))
      arg = br*(obc(1)+arg2)
      arg2 = 2.0_wp*arg2+obc(1)+obc(3)*br*br

      th = tanh(arg)
      ch = cosh(arg)

      br = 1.0_wp/(s1-v1*th)
      ! Include GBMV2-like scaling
      br = c1*br

      dpsi = ch*(s1-v1*th)
      dpsi = s2*v1*arg2/(dpsi*dpsi)
      dpsi = c1*dpsi

      brad(iat) = br
      brdr(:, :, iat) = brdr(:, :, iat) * dpsi

   end do

end subroutine compute_bornr


pure subroutine compute_psi(nat, xyz, list, vdwr, rho, psi, dpsidr)
   !> Number of atoms
   integer, intent(in) :: nat
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)
   !> Neighbourlist
   type(adjacency_list), intent(in) :: list
   !> Van-der-Waals radii
   real(wp), intent(in) :: vdwr(:)
   !> Descreened van-der-Waals radii
   real(wp), intent(in) :: rho(:)
   !> Integrated value of Psi
   real(wp), intent(out) :: psi(:)
   !> Derivative of Psi w.r.t. cartesian coordinates
   real(wp), intent(out) :: dpsidr(:, :, :)
   real(wp), allocatable :: dpsitr(:, :)

   integer  :: iat, jat, img, inl
   real(wp) :: vec(3), r, rhoi, rhoj
   real(wp) :: gi, gj, ap, am, lnab, rhab, ab, dgi, dgj
   real(wp) :: drjj(3)
   real(wp) :: rh1, rhr1, r24, r1, aprh1, r12
   real(wp) :: rvdwi, rvdwj
   logical :: ijov, jiov

   allocate(dpsitr(3, nat))
   psi(:) = 0.0_wp
   dpsidr(:, :, :) = 0.0_wp
   dpsitr(:, :) = 0.0_wp

   do iat = 1, nat
      inl = list%inl(iat)
      do img = 1, list%nnl(iat)
         jat = list%nlat(inl+img)

         vec(:) = xyz(:, iat) - xyz(:, jat)
         r = norm2(vec)

         rhoi = rho(iat)
         rhoj = rho(jat)
         rvdwi = vdwr(iat)
         rvdwj = vdwr(jat)

         ijov = r < (rvdwi+rhoj)
         jiov = r < (rhoi+rvdwj)

         if (.not.(ijov .or. jiov)) then
            ! nonoverlaping spheres
            if(abs(rhoi-rhoj) < 1.e-8_wp) then
               ! equal reduced radii
               r1 = 1.0_wp/r
               ap = r+rhoj
               am = r-rhoj
               ab = ap*am
               rhab = rhoj/ab
               lnab = 0.5_wp*log(am/ap)*r1
               gi = rhab+lnab
               dgi = -2.0_wp*rhab/ab+(rhab-lnab)*r1*r1
               ! accumulate psi
               psi(iat) = psi(iat)+gi
               psi(jat) = psi(jat)+gi
               ! accumulate psi gradient
               drjj(:) = dgi*vec(:)
               dpsitr(:, iat) = dpsitr(:, iat)+drjj(:)
               dpsidr(:, jat, iat) = dpsidr(:, jat, iat)-drjj(:)
               dpsitr(:, jat) = dpsitr(:, jat)-drjj(:)
               dpsidr(:, iat, jat) = dpsidr(:, iat, jat)+drjj(:)
            else
               ! unequal reduced radii
               ! ij contribution
               r1 = 1.0_wp/r
               ap = r+rhoj
               am = r-rhoj
               ab = ap*am
               rhab = rhoj/ab
               lnab = 0.5_wp*log(am/ap)*r1
               gi = rhab+lnab
               dgi = -2.0_wp*rhab/ab+(rhab-lnab)*r1*r1
               ! ji contribution
               ap = r+rhoi
               am = r-rhoi
               ab = ap*am
               rhab = rhoi/ab
               lnab = 0.5_wp*log(am/ap)*r1
               gj = rhab+lnab
               dgj = -2.0_wp*rhab/ab+(rhab-lnab)*r1*r1
               ! accumulate psi
               psi(iat) = psi(iat)+gi
               psi(jat) = psi(jat)+gj
               ! accumulate psi gradient
               drjj(:) = dgi*vec(:)
               dpsitr(:, iat) = dpsitr(:, iat)+drjj(:)
               dpsidr(:, jat, iat) = dpsidr(:, jat, iat)-drjj(:)

               drjj(:) = dgj*vec(:)
               dpsitr(:, jat) = dpsitr(:, jat)-drjj(:)
               dpsidr(:, iat, jat) = dpsidr(:, iat, jat)+drjj(:)
            end if

         else if (.not.ijov .and. jiov) then

            ! ij contribution
            r1 = 1.0_wp/r
            ap = r+rhoj
            am = r-rhoj
            ab = ap*am
            rhab = rhoj/ab
            lnab = 0.5_wp*log(am/ap)*r1
            gi = rhab+lnab
            dgi = -2.0_wp*rhab/ab+(rhab-lnab)*r1*r1
            ! accumulate psi
            psi(iat) = psi(iat)+gi
            ! accumulate psi gradient
            drjj(:) = dgi*vec(:)
            dpsitr(:, iat) = dpsitr(:, iat)+drjj(:)
            dpsidr(:, jat, iat) = dpsidr(:, jat, iat)-drjj(:)

            if((r+rhoi) > rvdwj) then
               ! ji contribution
               r1 = 1.0_wp/r
               r12 = 0.5_wp*r1
               r24 = r12*r12

               ap = r+rhoi
               am = r-rhoi
               rh1 = 1.0_wp/rvdwj
               rhr1 = 1.0_wp/ap
               aprh1 = ap*rh1
               lnab = log(aprh1)

               gj = rh1-rhr1+r12*(0.5_wp*am*(rhr1-rh1*aprh1)-lnab)

               dgj = rhr1*rhr1*(1.0_wp-0.25_wp*am*r1*(1.0_wp+aprh1*aprh1))+ &
                  &         rhoi*r24*(rhr1-rh1*aprh1)+ &
                  &         r12*(r1*lnab-rhr1)
               dgj = dgj*r1
               ! accumulate psi
               psi(jat) = psi(jat)+gj
               ! accumulate psi gradient
               drjj(:) = dgj*vec(:)
               dpsitr(:, jat) = dpsitr(:, jat)-drjj(:)
               dpsidr(:, iat, jat) = dpsidr(:, iat, jat)+drjj(:)
            end if

         else if (ijov .and. .not.jiov) then

            if((r+rhoj) > rvdwi) then
               ! ij contribution
               r1 = 1.0_wp/r
               r12 = 0.5_wp*r1
               r24 = r12*r12

               ap = r+rhoj
               am = r-rhoj
               rh1 = 1.0_wp/rvdwi
               rhr1 = 1.0_wp/ap
               aprh1 = ap*rh1
               lnab = log(aprh1)

               gi = rh1-rhr1+r12*(0.5_wp*am*(rhr1-rh1*aprh1)-lnab)

               dgi = rhr1*rhr1*(1.0_wp-0.25_wp*am*r1*(1.0_wp+aprh1*aprh1))+ &
                  &         rhoj*r24*(rhr1-rh1*aprh1)+ &
                  &         r12*(r1*lnab-rhr1)
               dgi = dgi*r1
               ! accumulate psi
               psi(iat) = psi(iat)+gi
               ! accumulate psi gradient
               drjj(:) = dgi*vec(:)
               dpsitr(:, iat) = dpsitr(:, iat)+drjj(:)
               dpsidr(:, jat, iat) = dpsidr(:, jat, iat)-drjj(:)
            end if

            ! ji contribution
            ap = r+rhoi
            am = r-rhoi
            ab = ap*am
            rhab = rhoi/ab
            lnab = 0.5_wp*log(am/ap)*r1
            gj = rhab+lnab
            dgj = -2.0_wp*rhab/ab+(rhab-lnab)*r1*r1
            ! accumulate psi
            psi(jat) = psi(jat)+gj
            ! accumulate psi gradient
            drjj(:) = dgj*vec(:)
            dpsitr(:, jat) = dpsitr(:, jat)-drjj(:)
            dpsidr(:, iat, jat) = dpsidr(:, iat, jat)+drjj(:)

         else if (ijov .and. jiov) then
            ! overlaping spheres
            if((r+rhoj) > rvdwi) then
               ! ij contribution
               r1 = 1.0_wp/r
               r12 = 0.5_wp*r1
               r24 = r12*r12

               ap = r+rhoj
               am = r-rhoj
               rh1 = 1.0_wp/rvdwi
               rhr1 = 1.0_wp/ap
               aprh1 = ap*rh1
               lnab = log(aprh1)

               gi = rh1-rhr1+r12*(0.5_wp*am*(rhr1-rh1*aprh1)-lnab)

               dgi = rhr1*rhr1*(1.0_wp-0.25_wp*am*r1*(1.0_wp+aprh1*aprh1))+ &
                  &         rhoj*r24*(rhr1-rh1*aprh1)+ &
                  &         r12*(r1*lnab-rhr1)
               dgi = dgi*r1
               ! accumulate psi
               psi(iat) = psi(iat)+gi
               ! accumulate psi gradient
               drjj(:) = dgi*vec(:)
               dpsitr(:, iat) = dpsitr(:, iat)+drjj(:)
               dpsidr(:, jat, iat) = dpsidr(:, jat, iat)-drjj(:)
            end if

            if((r+rhoi) > rvdwj) then
               ! ji contribution
               r1 = 1.0_wp/r
               r12 = 0.5_wp*r1
               r24 = r12*r12

               ap = r+rhoi
               am = r-rhoi
               rh1 = 1.0_wp/rvdwj
               rhr1 = 1.0_wp/ap
               aprh1 = ap*rh1
               lnab = log(aprh1)

               gj = rh1-rhr1+r12*(0.5_wp*am*(rhr1-rh1*aprh1)-lnab)

               dgj = rhr1*rhr1*(1.0_wp-0.25_wp*am*r1*(1.0_wp+aprh1*aprh1))+ &
                  &         rhoi*r24*(rhr1-rh1*aprh1)+ &
                  &         r12*(r1*lnab-rhr1)
               dgj = dgj*r1
               ! accumulate psi
               psi(jat) = psi(jat)+gj
               ! accumulate psi gradient
               drjj(:) = dgj*vec(:)
               dpsitr(:, jat) = dpsitr(:, jat)-drjj(:)
               dpsidr(:, iat, jat) = dpsidr(:, iat, jat)+drjj(:)
            end if

         end if

      end do
   end do

   ! save one-center terms
   do iat = 1, nat
      dpsidr(:, iat, iat) = dpsitr(:, iat)
   end do

end subroutine compute_psi

end module tblite_solvation_born
