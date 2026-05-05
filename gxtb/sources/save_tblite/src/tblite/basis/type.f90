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

!> @file tblite/basis/type.f90
!> Provides data types for managing basis set information

!> Gaussian type basis set data
module tblite_basis_type
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use tblite_basis_cache, only : basis_cache, cgto_cache
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   public :: new_basis, new_cgto
   public :: get_cutoff, integral_cutoff, basis_set, maxg

   !> Maximum contraction length of basis functions.
   !> The limit is chosen as twice the maximum size returned by the STO-NG expansion
   integer, parameter :: maxg = 12

   !> Default accuracy for basis set creation
   real(wp), parameter :: default_accuracy = 1.0_wp

   !> Two over pi
   real(wp), parameter :: top = 2.0_wp / pi
   !> Double factorial, see OEIS A001147
   real(wp), parameter :: dfactorial(8) = &
      & [1.0_wp,1.0_wp,3.0_wp,15.0_wp,105.0_wp,945.0_wp,10395.0_wp,135135.0_wp]

   !> Wrapped cgto_type for creation of lists with variable CGTOs
   type, public :: cgto_container
      !> Actual CGTO
      class(cgto_type), allocatable :: raw
   end type cgto_container

   !> Contracted Gaussian type basis function
   type, public :: cgto_type
      !> Angular momentum of this basis function
      integer :: ang = -1
      !> Contraction length of this basis function
      integer :: nprim = 0
      !> Exponent of the primitive Gaussian functions
      real(wp) :: alpha(maxg) = 0.0_wp
      !> Contraction coefficients of the primitive Gaussian functions,
      !> might contain normalization
      real(wp) :: coeff(maxg) = 0.0_wp
   contains 
      !> Update CGTO cache
      procedure :: update => cgto_update
      !> Get normalization factor and its derivatives
      procedure :: get_normalization
      !> Get (scaled) coefficient of the CGTO
      procedure :: get_coeffs
      !> Get coefficient derivatives of the CGTO
      procedure :: get_coeff_derivs
      !> Get effective atomic charge
      procedure :: get_qeff
   end type cgto_type

   !> Collection of information regarding the basis set of a system
   type, public :: basis_type
      !> Maximum angular momentum of all basis functions,
      !> used to determine scratch size in integral calculation
      integer :: maxl = 0
      !> Number of shells in this basis set
      integer :: nsh = 0
      !> Number of spherical atomic orbitals in this basis set
      integer :: nao = 0
      !> Number of cartesian atomic orbitals in this basis set
      integer :: nao_cart = 0
      !> Integral cutoff as maximum exponent of Gaussian product theoreom to consider
      real(wp) :: intcut = 0.0_wp
      !> Smallest primitive exponent in the basis set
      real(wp) :: min_alpha = huge(0.0_wp)
      !> Number of shells for each species
      integer, allocatable :: nsh_id(:)
      !> Number of shells for each atom
      integer, allocatable :: nsh_at(:)
      !> Number of spherical atomic orbitals for each shell
      integer, allocatable :: nao_sh(:)
      !> Number of cartesian atomic orbitals for each shell
      integer, allocatable :: nao_cart_sh(:)
      !> Index offset for each shell in the atomic orbital space
      integer, allocatable :: iao_sh(:)
      !> Index offset for each shell in the cartesian atomic orbital space
      integer, allocatable :: iao_cart_sh(:)
      !> Index offset for each atom in the shell space
      integer, allocatable :: ish_at(:)
      !> Mapping from spherical atomic orbitals to the respective atom
      integer, allocatable :: ao2at(:)
      !> Mapping from spherical atomic orbitals to the respective shell
      integer, allocatable :: ao2sh(:)
      !> Mapping from shells to the respective atom
      integer, allocatable :: sh2at(:)
      !> Contracted Gaussian basis functions forming the basis set
      type(cgto_container), allocatable :: cgto(:, :)
      !> Optional scaled contracted Gaussian basis functions for the H0 construction
      type(cgto_container), allocatable :: cgto_h0(:, :)
      !> Logical flag indicating a charge-dependent basis set
      logical :: charge_dependent = .false.
   contains
      !> Update basis cache
      procedure :: update => basis_update
      !> Get basis set gradient contributions
      procedure :: get_basis_gradient
   end type basis_type

   !> Get optimal real space cutoff for integral evaluation
   interface get_cutoff
      module procedure :: get_cutoff
   end interface get_cutoff

   !> Possible basis sets
   type :: enum_basis_set
      !> Stewart type Slater basis set expansion 
      integer :: sto_ng = 1
      !> Charge-dependent valence single-zeta q-vSZP basis set
      integer :: q_vszp = 2
      !> Custom basis set
      integer :: custom = 3
   end type enum_basis_set

   !> Actual enumerator for possible repulsion kernels
   type(enum_basis_set), parameter :: basis_set = enum_basis_set()

contains

!> Create a new basis set
subroutine new_basis(self, mol, nshell, cgto, accuracy, cgto_h0, charge_dependent)
   !> Instance of the basis set data
   type(basis_type), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Number of shells per species
   integer, intent(in) :: nshell(:)
   !> Contracted Gaussian basis functions for each shell and species
   type(cgto_container), intent(in) :: cgto(:, :)
   !> Calculation accuracy
   real(wp), intent(in), optional :: accuracy
   !> Optional scaled contracted Gaussian basis functions for the H0 construction
   type(cgto_container), intent(in), optional :: cgto_h0(:, :)
   !> Logical flag indicating a charge-dependent basis set
   logical, intent(in), optional :: charge_dependent

   integer :: iat, isp, ish, iao, ii, iicart
   real(wp) :: min_alpha, acc

   self%nsh_id = nshell

   self%cgto = cgto
   if (present(cgto_h0)) then
      self%cgto_h0 = cgto_h0
   end if

   ! Flag to decide if the basis set should be treated as charge-dependent
   if (present(charge_dependent)) then
      self%charge_dependent = charge_dependent
   else
      self%charge_dependent = .false.
   end if

   ! Integral cutoffs depending on the accuracy
   if (present(accuracy)) then
      acc = accuracy
   else
      acc = default_accuracy
   end if
   self%intcut = integral_cutoff(acc)

   ! Make count of shells for each atom
   self%nsh_at = nshell(mol%id)

   ! Create mapping between atoms and shells
   self%nsh = sum(self%nsh_at)
   allocate(self%ish_at(mol%nat), self%sh2at(self%nsh))
   ii = 0
   do iat = 1, mol%nat
      self%ish_at(iat) = ii
      do ish = 1, self%nsh_at(iat)
         self%sh2at(ii+ish) = iat
      end do
      ii = ii + self%nsh_at(iat)
   end do

   ! Make count of spherical orbitals for each shell
   allocate(self%nao_sh(self%nsh), self%nao_cart_sh(self%nsh))
   do iat = 1, mol%nat
      isp = mol%id(iat)
      ii = self%ish_at(iat)
      do ish = 1, self%nsh_at(iat)
         associate(p_cgto => self%cgto(ish, isp)%raw)
            self%nao_sh(ii+ish) = 2*p_cgto%ang + 1
            self%nao_cart_sh(ii+ish) = (p_cgto%ang + 1)*(p_cgto%ang + 2)/2
         end associate
      end do
   end do

   ! Create mapping between shells and spherical or cartesian orbitals,
   ! also map directly back to atoms
   self%nao = sum(self%nao_sh)
   self%nao_cart = sum(self%nao_cart_sh)
   allocate(self%iao_sh(self%nsh), self%iao_cart_sh(self%nsh), &
      & self%ao2sh(self%nao), self%ao2at(self%nao))
   ii = 0
   iicart = 0
   do ish = 1, self%nsh
      self%iao_sh(ish) = ii
      self%iao_cart_sh(ish) = iicart
      do iao = 1, self%nao_sh(ish)
         self%ao2sh(ii+iao) = ish
         self%ao2at(ii+iao) = self%sh2at(ish)
      end do
      ii = ii + self%nao_sh(ish)
      iicart = iicart + self%nao_cart_sh(ish)
   end do

   ii = 0
   do iat = 1, mol%nat
      isp = mol%id(iat)
      do ish = 1, nshell(isp)
         self%iao_sh(ish+self%ish_at(iat)) = ii
         associate(p_cgto => self%cgto(ish, isp)%raw)
            ii = ii + 2*p_cgto%ang + 1
         end associate
      end do
   end do

   min_alpha = huge(acc)
   do isp = 1, size(nshell)
      do ish = 1, nshell(isp)
         associate(p_cgto => self%cgto(ish, isp)%raw)
            self%maxl = max(self%maxl, p_cgto%ang)
            min_alpha = min(min_alpha, minval(p_cgto%alpha(:p_cgto%nprim)))
         end associate
      end do
   end do

   self%min_alpha = min_alpha

end subroutine new_basis

!> Update basis cache
subroutine basis_update(self, mol, cache, grad, wfn_aux)
   !> Instance of the basis type
   class(basis_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different runs
   type(basis_cache), intent(inout) :: cache
   !> Flag to indicate if gradients are required
   logical, intent(in) :: grad
   !> Optional auxiliary wavefunction data
   type(wavefunction_type), intent(in), optional :: wfn_aux

   integer :: max_nsh, ish, iat

   max_nsh = maxval(self%nsh_at)
   if (.not. allocated(cache%cgto)) allocate(cache%cgto(max_nsh, mol%nat))
   do iat = 1, mol%nat
      do ish = 1, self%nsh_at(iat)
         ! Update the cache without charge dependence
         associate(p_cgto => self%cgto(ish, mol%id(iat))%raw)
            call p_cgto%update(cache%cgto(ish, iat), grad)
         end associate
      end do
   end do 

   if (allocated(self%cgto_h0)) then
      if (.not. allocated(cache%cgto_h0)) allocate(cache%cgto_h0(max_nsh, mol%nat))
      do iat = 1, mol%nat
         do ish = 1, self%nsh_at(iat)
            ! Update the cache without charge dependence
            associate(p_cgto_h0 => self%cgto_h0(ish, mol%id(iat))%raw)
               call p_cgto_h0%update(cache%cgto_h0(ish, iat), grad)
            end associate
         end do
      end do 
   end if 
end subroutine basis_update


!> Contract gradients w.r.t. the basis set specific charge and CN
subroutine get_basis_gradient(self, mol, dEdcnbas, dEdqbas, dEdq, gradient, sigma)
   !> Instance of the basis type
   class(basis_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Derivative of the energy w.r.t. the basis set coordination number
   real(wp), intent(in), optional :: dEdcnbas(:)
   !> Derivative of the energy w.r.t. the basis set atomic charges
   real(wp), intent(in), optional :: dEdqbas(:)
   !> Derivative of the energy w.r.t. the atomic partial charges from an external model
   real(wp), intent(inout), optional :: dEdq(:)
   !> Derivative of the energy w.r.t. coordinate displacements
   real(wp), intent(inout) :: gradient(:, :)
   !> Derivative of the energy w.r.t. strain deformations
   real(wp), intent(inout) :: sigma(:, :)

end subroutine get_basis_gradient


!> Setup CGTO with coefficients and exponents
subroutine new_cgto(self, ng, l, expos, coeffs, norm)
   !> Instance of the cgto data
   type(cgto_type), intent(out) :: self
   !> Number of primitive Gaussian functions
   integer, intent(in) :: ng
   !> Azimudal quantum number of shell
   integer, intent(in) :: l
   !> Exponent of Gaussian function
   real(wp), intent(in) :: expos(:)
   !> Coefficients of the Gaussian function
   real(wp), intent(in) :: coeffs(:)
   !> Include normalization in contraction coefficients
   logical, intent(in) :: norm

   real(wp) :: normalizer(maxg)

   self%ang = l
   self%nprim = ng
   self%alpha(1:ng) = expos(1:ng)

   normalizer = 1.0_wp
   if (norm) then
      normalizer(1:ng) = (top*self%alpha(1:ng))**0.75_wp &
         & * sqrt(4.0_wp*self%alpha(1:ng))**l / sqrt(dfactorial(l+1))
   endif
   self%coeff(1:ng) = coeffs(1:ng) * normalizer(1:ng)

end subroutine new_cgto

!> Update CGTO cache
subroutine cgto_update(self, cache, grad, cn, q)
   !> Instance of the basis type
   class(cgto_type), intent(in) :: self
   !> Cached data between different runs
   type(cgto_cache), intent(inout) :: cache
   !> Flag to get gradient information
   logical, intent(in) :: grad
   !> Coordination number
   real(wp), intent(in), optional :: cn
   !> Effective charge
   real(wp), intent(in), optional :: q

   cache%qeff = 0.0_wp

   call self%get_normalization(cache, grad)
end subroutine cgto_update


!> Get normalization constant for the CGTOs
pure subroutine get_normalization(self, cache, grad)
   !> Instance of the basis set data
   class(cgto_type), intent(in) :: self
   !> Cached data between different runs
   type(cgto_cache), intent(inout) :: cache
   !> Flag to get gradient of normalization
   logical, intent(in) :: grad

   cache%norm = 1.0_wp
end subroutine get_normalization

!> Get (scaled) coefficients of the CGTOs
pure subroutine get_coeffs(self, cache, coeff)
   !> Instance of the basis set data
   class(cgto_type), intent(in) :: self
   !> Cached data between different runs
   type(cgto_cache), intent(in) :: cache
   !> Scaled coefficients of the CGTO
   real(wp), intent(out) :: coeff(:)
   
   coeff = self%coeff(1:self%nprim) * cache%norm
end subroutine get_coeffs

!> Get coefficient derivatives of the CGTOs
pure subroutine get_coeff_derivs(self, cache, dcoeff)
   !> Instance of the basis set data
   class(cgto_type), intent(in) :: self
   !> Cached data between different runs
   type(cgto_cache), intent(in) :: cache
   !> Derivative of coefficients of the CGTO w.r.t. effective charge
   real(wp), intent(out) :: dcoeff(:)
   
   dcoeff = 0.0_wp

end subroutine get_coeff_derivs

!> Update the effective environment-dependent charge
subroutine get_qeff(cgto, qat, cn, qeff, dqeffdcn, dqeffdq)
   !> Instance of the cgto data
   class(cgto_type), intent(in) :: cgto
   !> Atom-resolved charge
   real(wp), intent(in) :: qat
   !> Atom-resolved coordination number
   real(wp), intent(in) :: cn
   !> Scaled effective nuclear charge
   real(wp), intent(out) :: qeff
   !> Derivative of effective environment-dependent charge w.r.t. coordination number
   real(wp), intent(out), optional :: dqeffdcn
   !> Derivative of effective environment-dependent charge w.r.t. charge
   real(wp), intent(out), optional :: dqeffdq

   qeff = 0.0_wp
   if (present(dqeffdcn)) dqeffdcn = 0.0_wp
   if (present(dqeffdq)) dqeffdq = 0.0_wp

end subroutine get_qeff


!> Determine required real space cutoff for the basis set
pure function get_cutoff(self, acc) result(cutoff)
   !> Instance of the basis set data
   class(basis_type), intent(in) :: self
   !> Accuracy for the integral cutoff
   real(wp), intent(in), optional :: acc
   !> Required realspace cutoff
   real(wp) :: cutoff

   real(wp) :: intcut
   real(wp), parameter :: max_cutoff = 40.0_wp

   if (present(acc)) then
      intcut = integral_cutoff(acc)
   else
      intcut = self%intcut
   end if
   ! ai * aj * cutoff2 / (ai + aj) == intcut
   cutoff = min(sqrt(2.0_wp*intcut/self%min_alpha), max_cutoff)

end function get_cutoff


!> Create integral cutoff from accuracy value
pure function integral_cutoff(acc) result(intcut)
   !> Accuracy for the integral cutoff
   real(wp), intent(in) :: acc
   !> Integral cutoff
   real(wp) :: intcut

   real(wp), parameter :: min_intcut = 5.0_wp, max_intcut = 25.0_wp, &
      & max_acc = 1.0e-4_wp, min_acc = 1.0e+3_wp

   intcut = clip(max_intcut - 10*log10(clip(acc, min_acc, max_acc)), min_intcut, max_intcut)
end function integral_cutoff


pure function clip(val, min_val, max_val) result(res)
   real(wp), intent(in) :: val, min_val, max_val
   real(wp) :: res
   res = min(max(val, min_val), max_val)
end function clip

end module tblite_basis_type
