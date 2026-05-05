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

!> @file tblite/coulomb/multipole/gxtb.f90
!> Provides an implemenation of the error-function damped g-xTB multipole second-order electrostatic

!> Error-function damped anisotropic second-order electrostatics using a multipole expansion
module tblite_coulomb_multipole_gxtb
   use mctc_env, only : error_type, wp
   use mctc_io, only : structure_type
   use mctc_io_math, only : matdet_3x3, matinv_3x3
   use mctc_io_constants, only : pi
   use tblite_blas, only : dot, gemv, symv, gemm
   use tblite_container_cache, only : container_cache
   use tblite_coulomb_cache, only : coulomb_cache
   use tblite_coulomb_multipole_type, only : damped_multipole
   use tblite_scf_potential, only : potential_type
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   public :: new_gxtb_multipole


   !> Container to handle the error-function damped g-xTB multipole electrostatics
   type, public, extends(damped_multipole) :: gxtb_multipole
      !> Pairwise multipole damping radii
      real(wp), allocatable :: mrad(:, :)
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
   end type gxtb_multipole

   real(wp), parameter :: sqrtpi = sqrt(pi)
   character(len=*), parameter :: label = "erf-damped anisotropic electrostatics"

contains


!> Create a new error-function-damped g-xTB anisotropic electrostatics container
subroutine new_gxtb_multipole(self, mol, mrad, dip_scale, kdmp3, kdmp5, &
   & kdmp7, kdmp9, kexp3, kexp5, kexp7, kexp9)
   !> Instance of the multipole container
   type(gxtb_multipole), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Pairwise multipole damping radii
   real(wp), intent(in) :: mrad(:, :)
   !> Scaling factors for the dipole moments
   real(wp), intent(in) :: dip_scale(:)
   !> Damping prefactor for inverse quadratic interactions
   real(wp), intent(in), optional :: kdmp3
   !> Damping prefactor for inverse cubic interactions
   real(wp), intent(in), optional :: kdmp5
   !> Damping prefactor for inverse quartic interactions
   real(wp), intent(in), optional :: kdmp7
   !> Damping prefactor for inverse quintic interactions
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

   self%mrad = mrad
   self%dip_scale = dip_scale

   ! Prefactors and exponents for the damping function
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

end subroutine new_gxtb_multipole


!> Update cache from container
subroutine update(self, mol, cache)
   !> Instance of the multipole container
   class(gxtb_multipole), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   
   type(coulomb_cache), pointer :: ptr
   integer :: iat, jat, isp, jsp

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

   call self%get_multipole_matrix(mol, ptr, ptr%amat_sd, ptr%amat_dd, &
      & ptr%amat_sq, ptr%amat_dq, ptr%amat_qq)
end subroutine update


!> Evaluate multipole damping radius for a diatomic pair
pure subroutine get_mrad_pair(self, cache, iat, jat, isp, jsp, mrad)
   !> Instance of the multipole containeris
   class(gxtb_multipole), intent(in) :: self
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

   mrad = self%mrad(isp, jsp)

end subroutine get_mrad_pair


!> Evaluate multipole damping radius derivatives for a diatomic pair
pure subroutine get_mrad_derivs(self, cache, iat, jat, isp, jsp, mrad, &
   & dmraddcni, dmraddcnj)
   !> Instance of the multipole containeris
   class(gxtb_multipole), intent(in) :: self
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

   mrad = self%mrad(isp, jsp)
   dmraddcni = 0.0_wp
   dmraddcnj = 0.0_wp

end subroutine get_mrad_derivs


!> Evaluate damping function for a diatomic pair
pure subroutine get_damping_pair(self, r, rinv, mrad, fdmp3, fdmp5, &
   & fdmp7, fdmp9)
   !> Instance of the multipole containeris
   class(gxtb_multipole), intent(in) :: self
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

   rr = r - mrad
   fdmp3 = 0.5_wp * self%kdmp3 * (1.0_wp - erf(-self%kexp3 * rr))
   fdmp5 = 0.5_wp * self%kdmp5 * (1.0_wp - erf(-self%kexp5 * rr))
   fdmp7 = 0.5_wp * self%kdmp7 * (1.0_wp - erf(-self%kexp7 * rr))
   fdmp9 = 0.5_wp * self%kdmp9 * (1.0_wp - erf(-self%kexp9 * rr))

end subroutine get_damping_pair


!> Evaluate derivatives of the damping function for a diatomic pair
pure subroutine get_damping_derivs(self, r, rinv, mrad, fdmp3, fdmp5, &
   & fdmp7, fdmp9, dfdmp3dr, dfdmp5dr, dfdmp7dr, dfdmp9dr, &
   & dfdmp3dmrad, dfdmp5dmrad, dfdmp7dmrad, dfdmp9dmrad)
   !> Instance of the multipole container
   class(gxtb_multipole), intent(in) :: self
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

   real(wp) :: rr, rr2

   rr = r - mrad
   rr2 = rr**2
   fdmp3 = 0.5_wp * self%kdmp3 * (1.0_wp - erf(-self%kexp3 * rr))
   fdmp5 = 0.5_wp * self%kdmp5 * (1.0_wp - erf(-self%kexp5 * rr))
   fdmp7 = 0.5_wp * self%kdmp7 * (1.0_wp - erf(-self%kexp7 * rr))
   fdmp9 = 0.5_wp * self%kdmp9 * (1.0_wp - erf(-self%kexp9 * rr))

   ! Derivatives w.r.t. the interatomic distance
   dfdmp3dr = self%kdmp3 * self%kexp3 * exp(-self%kexp3**2 * rr2) / sqrtpi
   dfdmp5dr = self%kdmp5 * self%kexp5 * exp(-self%kexp5**2 * rr2) / sqrtpi
   dfdmp7dr = self%kdmp7 * self%kexp7 * exp(-self%kexp7**2 * rr2) / sqrtpi
   dfdmp9dr = self%kdmp9 * self%kexp9 * exp(-self%kexp9**2 * rr2) / sqrtpi

   ! Derivatives w.r.t. the damping radius
   dfdmp3dmrad = -dfdmp3dr
   dfdmp5dmrad = -dfdmp5dr
   dfdmp7dmrad = -dfdmp7dr
   dfdmp9dmrad = -dfdmp9dr

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

end module tblite_coulomb_multipole_gxtb
