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

!> @file tblite/coulomb/multipole/gfn2.f90
!> Provides an implemenation of the zero-damped GFN2 multipole electrostatic

!> Zero-damped anisotropic second-order electrostatics using a damped multipole expansion
module tblite_coulomb_multipole_gfn2
   use mctc_env, only : error_type, wp
   use mctc_io, only : structure_type
   use mctc_io_math, only : matdet_3x3, matinv_3x3
   use mctc_io_constants, only : pi
   use tblite_blas, only : dot, gemv, symv, gemm
   use tblite_container_cache, only : container_cache
   use tblite_coulomb_cache, only : coulomb_cache
   use tblite_coulomb_multipole_type, only : damped_multipole
   use tblite_scf_potential, only : potential_type
   use tblite_utils_average, only : average_type
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   public :: new_gfn2_multipole


   !> Container to handle zero-damped GFN2 multipole electrostatics
   type, public, extends(damped_multipole) :: gfn2_multipole
      !> Shift for the generation of the multipolar damping radii
      real(wp) :: shift = 0.0_wp
      !> Exponent for the generation of the multipolar damping radii
      real(wp) :: kradexp = 0.0_wp
      !> Maximum radius for the multipolar damping radii
      real(wp) :: rmax = 0.0_wp
      !> Base radii for the multipolar damping radii
      real(wp), allocatable :: rad(:)
      !> Averager for pairwise multipolar damping radii
      type(average_type), allocatable :: average
      !> Valence coordination number
      real(wp), allocatable :: valence_cn(:)
   contains
      !> Update cache from container
      procedure :: update
      !> Evaluate multipole damping radius
      procedure :: get_mrad_pair
      !> Evaluate multipole damping radius derivatives
      procedure :: get_mrad_derivs
      !> Evaluate pairwise damping factors
      procedure :: get_damping_pair
      !> Evaluate derivatives of damping factors
      procedure :: get_damping_derivs
   end type gfn2_multipole

   character(len=*), parameter :: label = "zero-damped anisotropic electrostatics"

contains


!> Create a new zero-damped GFN2 anisotropic electrostatics container
subroutine new_gfn2_multipole(self, mol, dkernel, qkernel, shift, kradexp, rmax, &
   & rad, average, vcn, kdmp3, kdmp5, kdmp7, kdmp9, kexp3, kexp5, kexp7, kexp9)
   !> Instance of the multipole container
   type(gfn2_multipole), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Kernel for on-site dipole exchange-correlation
   real(wp), intent(in) :: dkernel(:)
   !> Kernel for on-site quadrupolar exchange-correlation
   real(wp), intent(in) :: qkernel(:)
   !> Shift for the generation of the multipolar damping radii
   real(wp), intent(in) :: shift
   !> Exponent for the generation of the multipolar damping radii
   real(wp), intent(in) :: kradexp
   !> Maximum radius for the multipolar damping radii
   real(wp), intent(in) :: rmax
   !> Base radii for the multipolar damping radii
   real(wp), intent(in) :: rad(:)
   !> Averager for pairwise multipolar damping radii
   type(average_type), intent(in) :: average
   !> Valence coordination number
   real(wp), intent(in) :: vcn(:)
   !> Damping prefactor for inverse quadratic contributions
   real(wp), intent(in), optional :: kdmp3
   !> Damping prefactor for inverse cubic contributions
   real(wp), intent(in), optional :: kdmp5
   !> Damping prefactor for inverse quartic contributions
   real(wp), intent(in), optional :: kdmp7
   !> Damping prefactor for inverse quintic contributions
   real(wp), intent(in), optional :: kdmp9
   !> Damping function exponent for inverse quadratic interactions
   real(wp), intent(in), optional :: kexp3
   !> Damping function exponent for inverse cubic interactions
   real(wp), intent(in), optional :: kexp5
   !> Damping function exponent for inverse quartic interactions
   real(wp), intent(in), optional :: kexp7
   !> Damping function exponent for inverse quintic interactions
   real(wp), intent(in), optional :: kexp9

   self%label = label
   ! Prefactors and exponents for the zero-damping function
   if (present(kdmp3) .and. present(kexp3)) then
      self%kdmp3 = kdmp3
      self%kexp3 = kexp3
   else
      self%kdmp3 = 0.0_wp
      self%kexp3 = 0.0_wp
   end if
   if (present(kdmp5) .and. present(kexp5)) then
      self%kdmp5 = kdmp5
      self%kexp5 = kexp5
   else
      self%kdmp5 = 0.0_wp
      self%kexp5 = 0.0_wp
   end if
   if (present(kdmp7) .and. present(kexp7)) then
      self%kdmp7 = kdmp7
      self%kexp7 = kexp7
   else
      self%kdmp7 = 0.0_wp
      self%kexp7 = 0.0_wp
   end if
   if (present(kdmp9) .and. present(kexp9)) then
      self%kdmp9 = kdmp9
      self%kexp9 = kexp9
   else
      self%kdmp9 = 0.0_wp
      self%kexp9 = 0.0_wp
   end if
   ! Kernels for on-site exchange-correlation contributions
   self%dkernel = dkernel
   self%qkernel = qkernel
   ! Parameter for the multipolar damping radii
   self%shift = shift
   self%kradexp = kradexp
   self%rmax = rmax
   self%rad = rad
   self%average = average
   self%valence_cn = vcn

   allocate(self%dip_scale(mol%nid), source = 0.0_wp)

end subroutine new_gfn2_multipole


!> Update cache from container
subroutine update(self, mol, cache)
   !> Instance of the multipole container
   class(gfn2_multipole), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache

   type(coulomb_cache), pointer :: ptr

   call taint(cache, ptr)

   ! Multipole interaction matrices
   if (.not.allocated(ptr%amat_sd)) then
      allocate(ptr%amat_sd(3, mol%nat, mol%nat))
   end if
   if (.not.allocated(ptr%amat_dd)) then
      allocate(ptr%amat_dd(3, mol%nat, 3, mol%nat))
   end if
   if (.not.allocated(ptr%amat_sq)) then
      allocate(ptr%amat_sq(6, mol%nat, mol%nat))
   end if
   if (.not.allocated(ptr%amat_dq)) then
      allocate(ptr%amat_dq(3, mol%nat, 6, mol%nat))
   end if
   if (.not.allocated(ptr%amat_qq)) then
      allocate(ptr%amat_qq(6, mol%nat, 6, mol%nat))
   end if

   ! Multipole damping radii and their derivatives
   if (.not.allocated(ptr%mrad)) then
      allocate(ptr%mrad(mol%nat))
   end if
   if (.not.allocated(ptr%dmrdcn)) then
      allocate(ptr%dmrdcn(mol%nat))
   end if
   call get_mrad(mol, self%shift, self%kradexp, self%rmax, self%rad, &
      & self%valence_cn, ptr%cn, ptr%mrad, ptr%dmrdcn)

   call self%get_multipole_matrix(mol, ptr, ptr%amat_sd, ptr%amat_dd, &
      & ptr%amat_sq, ptr%amat_dq, ptr%amat_qq)

end subroutine update


!> Calculate multipole damping radii
subroutine get_mrad(mol, shift, kradexp, rmax, rad, valence_cn, cn, mrad, dmrdcn)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Shift for the generation of the multipolar damping radii
   real(wp), intent(in) :: shift
   !> Exponent for the generation of the multipolar damping radii
   real(wp), intent(in) :: kradexp
   !> Maximum radius for the multipolar damping radii
   real(wp), intent(in) :: rmax
   !> Base radii for the multipolar damping radii
   real(wp), intent(in) :: rad(:)
   !> Valence coordination number
   real(wp), intent(in) :: valence_cn(:)
   !> Coordination numbers for all atoms
   real(wp), intent(in) :: cn(:)
   !> Multipole damping radii for all atoms
   real(wp), intent(out) :: mrad(:)
   !> Derivative of multipole damping radii with repect to the coordination numbers
   real(wp), intent(out) :: dmrdcn(:)

   integer :: iat, izp
   real(wp) :: arg, exponential, delta_rad

   do iat = 1, mol%nat
      izp = mol%id(iat)
      arg = cn(iat) - valence_cn(izp) - shift
      exponential = exp(-kradexp*arg)
      delta_rad = (rmax - rad(izp)) / (1.0_wp + exponential)
      mrad(iat) = rad(izp) + delta_rad
      dmrdcn(iat) = delta_rad * kradexp * exponential / (1.0_wp + exponential)
   end do
end subroutine get_mrad


!> Evaluate multipole damping radius for a diatomic pair
pure subroutine get_mrad_pair(self, cache, iat, jat, isp, jsp, mrad)
   !> Instance of the multipole containeris
   class(gfn2_multipole), intent(in) :: self
   !> Reusable data container
   type(coulomb_cache), intent(in) :: cache
   !> Index of atom i
   integer, intent(in) :: iat
   !> Index of atom j
   integer, intent(in) :: jat
   !> Species of atom i
   integer, intent(in) :: isp
   !> Species of atom j
   integer, intent(in) :: jsp
   !> Multipole damping radius for the pair
   real(wp), intent(out) :: mrad

   mrad = self%average%value(cache%mrad(iat), cache%mrad(jat))

end subroutine get_mrad_pair


!> Evaluate multipole damping radius derivatives for a diatomic pair
pure subroutine get_mrad_derivs(self, cache, iat, jat, isp, jsp, mrad, &
   & dmraddcni, dmraddcnj)
   !> Instance of the multipole containeris
   class(gfn2_multipole), intent(in) :: self
   !> Reusable data container
   type(coulomb_cache), intent(in) :: cache
   !> Index of atom i
   integer, intent(in) :: iat
   !> Index of atom j
   integer, intent(in) :: jat
   !> Species of atom i
   integer, intent(in) :: isp
   !> Species of atom j
   integer, intent(in) :: jsp
   !> Multipole damping radius for the pair
   real(wp), intent(out) :: mrad
   !> Derivative of the multipole damping radius w.r.t. the CN of atom i
   real(wp), intent(out) :: dmraddcni
   !> Derivative of the multipole damping radius w.r.t. the CN of atom j
   real(wp), intent(out) :: dmraddcnj

   mrad = self%average%value(cache%mrad(iat), cache%mrad(jat))
   dmraddcni = self%average%deriv(cache%mrad(iat), cache%mrad(jat)) * cache%dmrdcn(iat)
   dmraddcnj = self%average%deriv(cache%mrad(jat), cache%mrad(iat)) * cache%dmrdcn(jat)

end subroutine get_mrad_derivs


!> Evaluate damping function for a diatomic pair
pure subroutine get_damping_pair(self, r, rinv, mrad, fdmp3, fdmp5, &
   & fdmp7, fdmp9)
   !> Instance of the multipole containeris
   class(gfn2_multipole), intent(in) :: self
   !> Interatomic distance
   real(wp), intent(in) :: r
   !> Inverse interatomic distance
   real(wp), intent(in) :: rinv
   !> Pairwise reference radius
   real(wp), intent(in) :: mrad
   !> Damping factor for inverse quadratic interactions
   real(wp), intent(out) :: fdmp3
   !> Damping factor for inverse cubic interactions
   real(wp), intent(out) :: fdmp5
   !> Damping factor for inverse quartic interactions
   real(wp), intent(out) :: fdmp7
   !> Damping factor for inverse quintic interactions
   real(wp), intent(out) :: fdmp9

   real(wp) :: rr

   rr = mrad * rinv
   fdmp3 = self%kdmp3 / (1.0_wp + 6.0_wp * rr**self%kexp3)
   fdmp5 = self%kdmp5 / (1.0_wp + 6.0_wp * rr**self%kexp5)
   fdmp7 = self%kdmp7 / (1.0_wp + 6.0_wp * rr**self%kexp7)
   fdmp9 = self%kdmp9 / (1.0_wp + 6.0_wp * rr**self%kexp9)

end subroutine get_damping_pair


!> Evaluate derivatives of the damping function for a diatomic pair
pure subroutine get_damping_derivs(self, r, rinv, mrad, fdmp3, fdmp5, &
   & fdmp7, fdmp9, dfdmp3dr, dfdmp5dr, dfdmp7dr, dfdmp9dr, &
   & dfdmp3dmrad, dfdmp5dmrad, dfdmp7dmrad, dfdmp9dmrad)
   !> Instance of the multipole container
   class(gfn2_multipole), intent(in) :: self
   !> Interatomic distance
   real(wp), intent(in) :: r
   !> Inverse interatomic distance
   real(wp), intent(in) :: rinv
   !> Averaged pairwise reference radius
   real(wp), intent(in) :: mrad
   !> Damping factor for inverse quadratic interactions
   real(wp), intent(out) :: fdmp3
   !> Damping factor for inverse cubic interactions
   real(wp), intent(out) :: fdmp5
   !> Damping factor for inverse quartic interactions
   real(wp), intent(out) :: fdmp7
   !> Damping factor for inverse quintic interactions
   real(wp), intent(out) :: fdmp9
   !> Derivative of damping factor for inverse quadratic interactions w.r.t. the distance
   real(wp), intent(out) :: dfdmp3dr
   !> Derivative of damping factor for inverse cubic interactions w.r.t. the distance
   real(wp), intent(out) :: dfdmp5dr
   !> Derivative of damping factor for inverse quartic interactions w.r.t. the distance
   real(wp), intent(out) :: dfdmp7dr
   !> Derivative of damping factor for inverse quintic interactions w.r.t. the distance
   real(wp), intent(out) :: dfdmp9dr
   !> Derivative of damping factor for inverse quadratic interactions w.r.t. the radius
   real(wp), intent(out) :: dfdmp3dmrad
   ! Derivative of damping factor for inverse cubic interactions w.r.t. the radius
   real(wp), intent(out) :: dfdmp5dmrad
   !> Derivative of damping factor for inverse quartic interactions w.r.t. the radius
   real(wp), intent(out) :: dfdmp7dmrad
   !> Derivative of damping factor for inverse quintic interactions w.r.t. the radius
   real(wp), intent(out) :: dfdmp9dmrad

   real(wp) :: rr, tmp

   rr = mrad * rinv
   fdmp3 = 1.0_wp / (1.0_wp + 6.0_wp * rr**self%kexp3)
   fdmp5 = 1.0_wp / (1.0_wp + 6.0_wp * rr**self%kexp5)
   fdmp7 = 1.0_wp / (1.0_wp + 6.0_wp * rr**self%kexp7)
   fdmp9 = 1.0_wp / (1.0_wp + 6.0_wp * rr**self%kexp9)

   ! Derivatives w.r.t. the interatomic distance and the damping radius
   tmp = -6.0_wp * self%kdmp3 * self%kexp3 * rr**self%kexp3 * fdmp3**2
   dfdmp3dr = -tmp * rinv
   dfdmp3dmrad = tmp / mrad
   tmp = -6.0_wp * self%kdmp5 * self%kexp5 * rr**self%kexp5 * fdmp5**2
   dfdmp5dr = -tmp * rinv
   dfdmp5dmrad = tmp / mrad
   tmp = -6.0_wp * self%kdmp7 * self%kexp7 * rr**self%kexp7 * fdmp7**2
   dfdmp7dr = -tmp * rinv
   dfdmp7dmrad = tmp / mrad
   tmp = -6.0_wp * self%kdmp9 * self%kexp9 * rr**self%kexp9 * fdmp9**2
   dfdmp9dr = -tmp * rinv
   dfdmp9dmrad = tmp / mrad

   ! Scale the damping function by the prefactors
   fdmp3 = self%kdmp3 * fdmp3
   fdmp5 = self%kdmp5 * fdmp5
   fdmp7 = self%kdmp7 * fdmp7
   fdmp9 = self%kdmp9 * fdmp9

end subroutine get_damping_derivs

!> Inspect container cache and reallocate it in case of type mismatch
subroutine taint(cache, ptr)
   !> Instance of the container cache
   type(container_cache), target, intent(inout) :: cache
   !> Reference to the container cache
   type(coulomb_cache), pointer, intent(out) :: ptr

   if (allocated(cache%raw)) then
      call view(cache, ptr)
      if (associated(ptr)) return
      deallocate(cache%raw)
   end if

   if (.not.allocated(cache%raw)) then
      block
         type(coulomb_cache), allocatable :: tmp
         allocate(tmp)
         call move_alloc(tmp, cache%raw)
      end block
   end if

   call view(cache, ptr)
end subroutine taint

!> Return reference to container cache after resolving its type
subroutine view(cache, ptr)
   !> Instance of the container cache
   type(container_cache), target, intent(inout) :: cache
   !> Reference to the container cache
   type(coulomb_cache), pointer, intent(out) :: ptr
   nullify(ptr)
   select type(target => cache%raw)
   type is(coulomb_cache)
      ptr => target
   end select
end subroutine view

end module tblite_coulomb_multipole_gfn2
