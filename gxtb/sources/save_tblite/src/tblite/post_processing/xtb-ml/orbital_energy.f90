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

!> @file tblite/post-processing/xtb-ml/orbital_energy.f90
!> Orbital energy based xtbml features
module tblite_xtbml_orbital_energy
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_convert, only : autoev
   use tblite_container, only : container_cache
   use tblite_integral_type, only : integral_type
   use tblite_output_format, only : format_string
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtbml_atomic_frontier, only : atomic_frontier_orbitals
   use tblite_xtbml_convolution, only : xtbml_convolution_type
   use tblite_xtbml_feature_type, only : xtbml_feature_type
   implicit none
   private
   character(len=*), parameter :: label = "orbital energy-based features"
   type, public, extends(xtbml_feature_type) :: xtbml_orbital_features_type

      real(wp),allocatable ::  response(:)
      real(wp),allocatable ::  egap(:, :)
      real(wp),allocatable ::  chempot(:, :)
      real(wp),allocatable ::  ehoao(:)
      real(wp),allocatable ::  eluao(:)
      real(wp),allocatable ::  ext_chempot(:,:)
      real(wp),allocatable ::  ext_egap(:,:)
      real(wp),allocatable ::  ext_eluao(:,:)
      real(wp),allocatable ::  ext_ehoao(:,:)
   contains
      procedure :: compute_features
      procedure :: compute_extended
      procedure, private :: allocate
      procedure, private :: allocate_extended
      procedure :: setup
   end type xtbml_orbital_features_type

contains

!> Setup the container
subroutine setup(self)
   class(xtbml_orbital_features_type) :: self
   self%label = label
   if (allocated(self%dict)) deallocate(self%dict)
   allocate(self%dict)
   if (allocated(self%dict_ext)) deallocate(self%dict_ext)
   allocate(self%dict_ext)
end subroutine setup

!> Allocate memory for the features
subroutine allocate(self, nat, nspin)
   !> Instance of feature container
   class(xtbml_orbital_features_type), intent(inout) :: self
   !> Number of atoms
   integer, intent(in) :: nat
   !> Number of spins
   integer, intent(in) :: nspin

   allocate(self%response(nat), source=0.0_wp)
   allocate(self%egap(nat, nspin), source=0.0_wp)
   allocate(self%chempot(nat, nspin), source=0.0_wp)
   allocate(self%ehoao(nat), source=0.0_wp)
   allocate(self%eluao(nat), source=0.0_wp)

end subroutine allocate

!> Allocate memory for the extended features
subroutine allocate_extended(self, nat, n_a)
   !> Instance of feature container
   class(xtbml_orbital_features_type), intent(inout) :: self
   !> Number of atoms
   integer, intent(in) :: nat
   !> Number of convolution kernels
   integer, intent(in) :: n_a

   allocate(self%ext_chempot(nat, n_a), source=0.0_wp)
   allocate(self%ext_egap(nat, n_a), source=0.0_wp)
   allocate(self%ext_eluao(nat, n_a), source=0.0_wp)
   allocate(self%ext_ehoao(nat, n_a), source=0.0_wp)

end subroutine allocate_extended

subroutine compute_features(self, mol, wfn, integrals, calc, cache_list)
   !> Instance of feature container
   class(xtbml_orbital_features_type), intent(inout) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wavefunction strcuture data
   type(wavefunction_type), intent(in) :: wfn
   !> Integral container
   type(integral_type), intent(in) :: integrals
   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc
   !> Cache list for storing caches of various interactions
   type(container_cache), intent(inout) :: cache_list(:)

   integer :: i, j, nspin, spin
   real(wp) :: nel_
   character(len=6), allocatable :: spin_label(:)

   if (wfn%nuhf > 0 .or. wfn%nspin > 1) then
      nspin = 2
   else
      nspin = 1
   end if

   if (nspin > 1) then
      spin_label = [character(len=6) :: "_alpha", "_beta"]
   else
      spin_label = [""]
   end if
   call self%allocate(mol%nat, nspin)
   self%label = label

   do spin = 1, nspin
      if (wfn%nspin > 1) then
         call atomic_frontier_orbitals(wfn%focc(:, spin), wfn%emo(:, spin)*autoev, &
            calc%bas%ao2at, wfn%coeff(:, :, spin), integrals%overlap(:, :), &
            self%response, self%egap(:, spin), self%chempot(:, spin), self%ehoao, &
            self%eluao)
      else
         call atomic_frontier_orbitals(wfn%focc(:, spin), wfn%emo(:, 1)*autoev, &
            calc%bas%ao2at, wfn%coeff(:, :, 1), integrals%overlap(:, :), &
            self%response, self%egap(:, spin), self%chempot(:, spin), self%ehoao, &
            self%eluao)
      endif
      associate(dict => self%dict)
         call dict%add_entry(trim("response"//spin_label(spin)), self%response)
         call dict%add_entry(trim("gap"//spin_label(spin)), self%egap(:, spin))
         call dict%add_entry(trim("chem_pot"//spin_label(spin)), self%chempot(:, spin))
         call dict%add_entry(trim("HOAO"//spin_label(spin)), self%ehoao)
         call dict%add_entry(trim("LUAO"//spin_label(spin)), self%eluao)
      end associate
   end do
end subroutine compute_features

subroutine compute_extended(self, mol, wfn, integrals, calc, cache_list, convolution)
   !> Instance of feature container
   class(xtbml_orbital_features_type), intent(inout) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wavefunction strcuture data
   type(wavefunction_type), intent(in) :: wfn
   !> Integral container
   type(integral_type), intent(in) :: integrals
   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc
   !> Cache list for storing caches of various interactions
   type(container_cache), intent(inout) :: cache_list(:)
   !> Convolution container
   type(xtbml_convolution_type), intent(in) :: convolution

   integer :: n, i, nspin, spin
   character(len=:), allocatable :: a_label
   real(wp), allocatable :: beta(:, :, :)
   character(len=6), allocatable :: spin_label(:)
   nspin = size(self%chempot, 2)
   call self%allocate_extended(mol%nat, convolution%n_a)
   if (nspin > 1) then
      spin_label = [character(len=6) :: "_alpha", "_beta"]
   else
      spin_label = [""]
   end if

   n=convolution%n_a
   allocate(beta(mol%nat, mol%nat, n))
   call get_beta(convolution%kernel, convolution%cn, beta)

   do spin = 1, nspin
      call get_chem_pot_ext(beta, self%chempot(:, spin), self%ext_chempot)

      call get_e_gap_ext(beta, self%egap(:, spin), self%ext_egap)

      call get_ehoao_ext(self%ext_chempot, self%ext_egap, self%ext_ehoao)

      call get_eluao_ext(self%ext_chempot, self%ext_egap, self%ext_eluao)

      associate( dict => self%dict_ext)
         do i = 1, n
            a_label = "_"//trim(adjustl(format_string(convolution%a(i), '(f12.2)')) )
            if (a_label .eq. "_1.00") a_label = ""
            call dict%add_entry(trim("ext_gap"//spin_label(spin))//a_label, self%ext_egap(:, i))
            call dict%add_entry(trim("ext_chem_pot"//spin_label(spin))//a_label, self%ext_chempot(:, i))
            call dict%add_entry(trim("ext_HOAO"//spin_label(spin))//a_label, self%ext_ehoao(:, i))
            call dict%add_entry(trim("ext_LUAO"//spin_label(spin))//a_label, self%ext_eluao(:, i))
         end do
      end associate
   end do
end subroutine compute_extended

subroutine get_beta(kernel, cn, beta)
   intrinsic :: sum
   !> Convolution kernel
   real(wp), intent(in) :: kernel(:, :, :)
   !> Coordination number
   real(wp), intent(in) :: cn(:, :)
   !> Convolution coefficients
   real(wp), intent(out) :: beta(:, :, :)
   real(wp) :: sigma_tot
   real(wp) :: damp_func
   integer :: nat, n_a
   integer :: a, b, k
   real(wp), parameter :: eps = 1.0e-12_wp 

   nat = size(kernel, 1)
   n_a = size(kernel, 3)
   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(kernel, cn, beta, n_a, nat) private(a, b, k)
   do k = 1, n_a
      do a = 1, nat
         do b = 1, nat
            beta(a, b, k) = 1.0_wp/(kernel(a, b, k) *(cn(a, k) +1.0_wp)) 
         end do
      end do
   end do
   
end subroutine get_beta

!> Get the extended chemical potential
subroutine get_chem_pot_ext(beta, chempot, chempot_ext)
   !> beta matrix
   real(wp), intent(in) :: beta(:, :, :)
   !> Chemical potential
   real(wp), intent(in) :: chempot(:)
   !> Extended chemical potential
   real(wp), intent(out) :: chempot_ext(:, :)
   integer :: a, b, k
   !$omp parallel do default(none) schedule(runtime) reduction(+:chempot_ext) &
   !$omp shared(beta, chempot) private(a, b, k)
   do k = 1, size(beta, 3)
      do a = 1, size(chempot, 1)
         do b = 1, size(chempot, 1)
            chempot_ext(a, k) = chempot_ext(a, k) + beta(a, b, k) * chempot(b)
         end do
      end do
   end do
end subroutine get_chem_pot_ext

!> Get the extended energy gap
subroutine get_e_gap_ext(beta, e_gap, e_gap_ext)
   !> beta amtrix
   real(wp), intent(in) :: beta(:, :, :)
   !> Energy gap
   real(wp), intent(in) :: e_gap(:)
   !> Extended energy gap
   real(wp), intent(out) :: e_gap_ext(:, :)
   integer :: a, b, k
   !$omp parallel do default(none) schedule(runtime) reduction(+:e_gap_ext) &
   !$omp shared(beta, e_gap) private(a, b, k)
   do k = 1, size(beta, 3)
      do a = 1, size(e_gap)
         do b = 1, size(e_gap)
            e_gap_ext(a, k) = e_gap_ext(a, k) + beta(a, b, k) * e_gap(b)
         end do
      end do
   end do
   
end subroutine get_e_gap_ext

!> Compute the extended HOAO 
subroutine get_ehoao_ext(chempot_ext, e_gap_ext, ehoao_ext)
   !> extended chemical potential
   real(wp), intent(in) :: chempot_ext(:, :)
   !> extended gap 
   real(wp), intent(in) :: e_gap_ext(:, :)
   !> extended HOAO
   real(wp), intent(out) :: ehoao_ext(:, :)
   integer :: a, k
   do k = 1, size(chempot_ext, 2)
      do a = 1, size(chempot_ext, 1)
         ehoao_ext(a, k) = chempot_ext(a, k) - e_gap_ext(a, k) / 2.0_wp
      end do
   end do
end subroutine get_ehoao_ext

!> Compute the extended LUAO
subroutine get_eluao_ext(chempot_ext, e_gap_ext, eluao_ext)
   !> extended chemical potential
   real(wp), intent(in) :: chempot_ext(:, :)
   !> extended gap 
   real(wp), intent(in) :: e_gap_ext(:, :)
   !> extended LUAO
   real(wp), intent(out) :: eluao_ext(:, :)
   integer :: a, k
   do k = 1, size(chempot_ext, 2)
      do a = 1, size(chempot_ext, 1)
         eluao_ext(a, k) = chempot_ext(a, k) + e_gap_ext(a, k) / 2.0_wp
      end do
   end do
end subroutine get_eluao_ext

end module tblite_xtbml_orbital_energy
