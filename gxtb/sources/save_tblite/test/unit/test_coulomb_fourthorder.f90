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

module test_coulomb_fourthorder
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type, new
   use mctc_ncoord, only : new_ncoord, ncoord_type, cn_count
   use mstore, only : get_structure
   use tblite_basis_slater, only : slater_to_gauss
   use tblite_basis_type
   use tblite_ceh_ceh, only : get_effective_qat
   use tblite_container_cache, only : container_cache
   use tblite_coulomb_fourthorder, only : onsite_fourthorder, new_onsite_fourthorder
   use tblite_coulomb_type, only : coulomb_type
   use tblite_cutoff, only : get_lattice_points
   use tblite_scf, only : new_potential, potential_type   
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   implicit none
   private

   public :: collect_coulomb_fourthorder

   real(wp), parameter :: cutoff = 25.0_wp
   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr1 = 1e5*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

   abstract interface
      subroutine coulomb_maker(coulomb, mol, shell)
         import :: coulomb_type, structure_type
         class(coulomb_type), allocatable, intent(out) :: coulomb
         type(structure_type), intent(in) :: mol
         logical, intent(in) :: shell
      end subroutine coulomb_maker
   end interface

   abstract interface
      subroutine charge_maker(wfn, mol, nshell, error)
         import :: wavefunction_type, structure_type, error_type
         class(wavefunction_type), intent(inout) :: wfn
         type(structure_type), intent(in) :: mol
         integer, optional, intent(in) :: nshell(:)
         type(error_type), allocatable, intent(out) :: error
      end subroutine charge_maker
   end interface

contains


!> Collect all exported unit tests
subroutine collect_coulomb_fourthorder(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("energy-atom-gxtb-m01", test_e_gxtb_m01_atom), &
      new_unittest("energy-shell-gxtb-m01", test_e_gxtb_m01_shell), &
      new_unittest("energy-atom-gxtb-m02", test_e_gxtb_m02_atom), &
      new_unittest("energy-shell-gxtb-m02", test_e_gxtb_m02_shell), &
      new_unittest("energy-atom-pbc-gxtb", test_e_gxtb_oxacb_atom), &
      new_unittest("energy-shell-pbc-gxtb", test_e_gxtb_oxacb_shell), &
      new_unittest("energy-atom-pbcsc-gxtb", test_e_gxtb_oxacb_sc), &
      new_unittest("potential-atom-gxtb-m01", test_p_gxtb_m01_atom), &
      new_unittest("potential-shell-gxtb-m01", test_p_gxtb_m01_shell), &
      new_unittest("potential-atom-gxtb-m02", test_p_gxtb_m02_atom), &
      new_unittest("potential-shell-gxtb-m02", test_p_gxtb_m02_shell), &
      new_unittest("potential-gradient-atom-gxtb-effceh-lih", test_pg_ceh_lih_atom), &
      new_unittest("potential-gradient-atom-gxtb-effceh-m15", test_pg_ceh_m15_atom), &
      new_unittest("potential-gradient-atom-gxtb-effceh-m16", test_pg_ceh_m16_shell), &
      new_unittest("potential-sigma-atom-gxtb-effceh-lih", test_ps_ceh_lih_atom), &
      new_unittest("potential-sigma-atom-gxtb-effceh-co2", test_ps_ceh_co2_atom), &
      new_unittest("potential-sigma-atom-gxtb-effceh-m05", test_ps_ceh_m05_atom), &
      new_unittest("potential-sigma-shell-gxtb-effceh-m17", test_ps_ceh_m17_shell) &
      ]

end subroutine collect_coulomb_fourthorder


!> Factory to setup the CEH basis set for testing of the potential gradient
subroutine make_basis(bas, mol, ng)
   type(basis_type), intent(out) :: bas
   type(structure_type), intent(in) :: mol
   integer, intent(in) :: ng

   integer, parameter :: nsh(20) = [ &
      & 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 2, 3]
   integer, parameter :: lsh(3, 20) = reshape([ &
      & 0, 0, 0,  0, 0, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0, &
      & 0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 2, &
      & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2, &
      & 0, 1, 0,  0, 1, 2], shape(lsh))
   integer, parameter :: pqn(3, 20) = reshape([ &
      & 1, 0, 0,  1, 0, 0,  2, 2, 0,  2, 2, 0,  2, 2, 0,  2, 2, 0, &
      & 2, 2, 0,  2, 2, 0,  2, 2, 0,  2, 2, 0,  3, 3, 0,  3, 3, 3, &
      & 3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3, &
      & 4, 4, 0,  4, 4, 3], shape(pqn))
   real(wp), parameter :: zeta(3, 20) = reshape([ &
      & 1.23363166_wp, 0.00000000_wp, 0.00000000_wp, 2.27004605_wp, 0.00000000_wp, 0.00000000_wp, &
      & 0.86185456_wp, 1.42017184_wp, 0.00000000_wp, 1.76817995_wp, 1.44095844_wp, 0.00000000_wp, &
      & 2.06339837_wp, 1.52051807_wp, 0.00000000_wp, 2.56058582_wp, 1.86484737_wp, 0.00000000_wp, &
      & 2.71233631_wp, 2.19848968_wp, 0.00000000_wp, 3.21585650_wp, 2.41309737_wp, 0.00000000_wp, &
      & 3.82146807_wp, 2.63063636_wp, 0.00000000_wp, 4.62721228_wp, 2.53599954_wp, 0.00000000_wp, &
      & 0.93221172_wp, 1.55333839_wp, 0.00000000_wp, 1.77220557_wp, 1.59942632_wp, 2.98596647_wp, &
      & 2.26040231_wp, 1.78718151_wp, 2.00990188_wp, 1.85259089_wp, 1.81733349_wp, 1.65269988_wp, &
      & 2.65701241_wp, 2.03189759_wp, 2.03883661_wp, 2.60609998_wp, 2.16530440_wp, 2.41888232_wp, &
      & 2.78818934_wp, 2.24732894_wp, 1.99081182_wp, 2.55424399_wp, 2.20946190_wp, 1.93619550_wp, &
      & 1.73713827_wp, 1.33788617_wp, 0.00000000_wp, 2.47982574_wp, 1.07250770_wp, 2.11920764_wp], &
      shape(zeta))

   integer :: isp, izp, ish, stat
   integer, allocatable :: nshell(:)
   type(cgto_container), allocatable :: cgto(:, :)

   nshell = nsh(mol%num)
   allocate(cgto(maxval(nshell), mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do ish = 1, nshell(isp)
         allocate(cgto_type :: cgto(ish, isp)%raw)
         call slater_to_gauss(ng, pqn(ish, izp), lsh(ish, izp), zeta(ish, izp), &
            cgto(ish, isp)%raw, .true., stat)
      end do
   end do

   call new_basis(bas, mol, nshell, cgto, accuracy=1.0_wp)

end subroutine make_basis


!> Factory to create onsite fourth-order objects using the g-xTB parameter
subroutine make_coulomb_gxtb(coulomb, mol, shell)

   !> New coulomb object
   class(coulomb_type), allocatable, intent(out) :: coulomb

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Return a shell resolved object
   logical, intent(in) :: shell
   
   !> Second derivative of Hubbard parameters or chemical hardness 
   !> for onsite fourth-order tight-binding
   real(wp), parameter :: p_hubbard_second_deriv = 0.0360000000_wp
   !> Global Parameter: Shell-dependence of the fourth-order tight-binding
   real(wp), parameter :: p_kshell(0:3) = &
      & [1.0000000000_wp, 1.1500000000_wp, 1.3000000000_wp, 1.4500000000_wp]

   integer, parameter :: shell_count(20) = [ &
      & 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 2, 3]

   real(wp), allocatable :: hubbard_second_derivs(:, :)
   integer :: isp, izp, ish
   integer, allocatable :: nshell(:)
   type(onsite_fourthorder), allocatable :: tmp

   allocate(tmp)
   if (shell) then
      allocate(hubbard_second_derivs(maxval(shell_count), mol%nid), source=0.0_wp)
      do isp = 1, mol%nid
         izp = mol%num(isp)
         do ish = 1, shell_count(izp)
            hubbard_second_derivs(ish, isp) = p_kshell(ish-1) * p_hubbard_second_deriv
         end do
      end do
      nshell = shell_count(mol%num)
      call new_onsite_fourthorder(tmp, mol, hubbard_second_derivs, nshell)
   else
      allocate(hubbard_second_derivs(1, mol%nid))
      do isp = 1, mol%nid
         izp = mol%num(isp)
         hubbard_second_derivs(1, isp) = p_hubbard_second_deriv
      end do
      call new_onsite_fourthorder(tmp, mol, hubbard_second_derivs)
   end if
   call move_alloc(tmp, coulomb)

end subroutine make_coulomb_gxtb


!> Procedure to create CN based effective charges and gradients from CEH
!> only for testing purposes (CEH does not contain a fourth-order term)
subroutine get_charges_effceh(wfn, mol, nshell, error)

   !> New wavefunction object
   class(wavefunction_type), intent(inout) :: wfn

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Return shell-resolved charges
   integer, optional, intent(in) :: nshell(:)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: ceh_cov_radii(20) = 0.5 * [ &
      &  2.4040551903_wp,  1.8947380542_wp,  3.4227634078_wp,  3.5225408137_wp, &
      &  3.6150631704_wp,  2.8649682108_wp,  2.4695867541_wp,  2.3533691180_wp, &
      &  2.4992147462_wp,  3.3442607441_wp,  4.4665909451_wp,  4.3877250907_wp, &
      &  4.6647077385_wp,  4.2086223530_wp,  4.4750280107_wp,  4.2847281423_wp, &
      &  3.8560304959_wp,  3.9017061017_wp,  5.2392192639_wp,  5.1872031383_wp]
   real(wp), parameter :: pauling_en_ceh(20) = (1.0_wp/3.98_wp) * [ &
      &  1.9435211923_wp,  3.6116085622_wp,  2.4630915335_wp,  2.0658837656_wp, &
      &  2.3619778807_wp,  2.9484294262_wp,  3.8753937411_wp,  4.6235054741_wp, &
      &  3.9800000000_wp,  3.5124865506_wp,  2.3578254072_wp,  2.4225832022_wp, &
      &  2.1120078826_wp,  2.4607564741_wp,  2.7410779326_wp,  3.3517034720_wp, &
      &  4.1093492601_wp,  3.7979559518_wp,  2.4147937668_wp,  2.1974781961_wp]

   real(wp), allocatable :: lattr(:, :)
   real(wp), allocatable :: cn_en(:), dcn_endr(:, :, :), dcn_endL(:, :, :)
   class(ncoord_type), allocatable :: ncoord_en
   integer :: iat, ii, ish

   allocate(cn_en(mol%nat), dcn_endr(3, mol%nat, mol%nat), dcn_endL(3, 3, mol%nat))

   call new_ncoord(ncoord_en, mol, cn_count%erf_en, error, &
      rcov=ceh_cov_radii(mol%num), en=pauling_en_ceh(mol%num))
   if (allocated(error)) return
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
   call ncoord_en%get_coordination_number(mol, lattr, cn_en, dcn_endr, dcn_endL)

   call get_effective_qat(mol, cn_en, wfn%qat, &
      dcn_endr, dcn_endL, wfn%dqatdr, wfn%dqatdL)

   if (present(nshell)) then
      ii = 0
      do iat = 1, mol%nat
         do ish = 1, nshell(iat)
            wfn%qsh(ii+ish, :) = wfn%qat(iat, :) / real(nshell(iat), wp)
            wfn%dqshdr(:, :, ii+ish, :) = wfn%dqatdr(:, :, iat, :) &
               & / real(nshell(iat), wp)
            wfn%dqshdL(:, :, ii+ish, :) = wfn%dqatdL(:, :, iat, :) &
               & / real(nshell(iat), wp)
         end do
         ii = ii + nshell(iat)
      end do
   end if

end subroutine get_charges_effceh


subroutine test_generic(error, mol, qat, qsh, make_coulomb, ref, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Atomic partial charges for this structure
   real(wp), intent(in) :: qat(:)

   !> Shell-resolved partial charges for this structure
   real(wp), intent(in), optional :: qsh(:)

   !> Factory to create new onsite fourth-order objects
   procedure(coulomb_maker) :: make_coulomb

   !> Reference value to check against
   real(wp), intent(in) :: ref

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   class(coulomb_type), allocatable :: coulomb
   type(container_cache) :: cache
   type(wavefunction_type) :: wfn
   type(potential_type) :: pot
   type(basis_type) :: bas
   real(wp) :: energy(mol%nat)
   real(wp) :: thr_

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   call make_basis(bas, mol, 6)
   call new_wavefunction(wfn, mol%nat, bas%nsh, bas%nao, 1, 0.0_wp, .true.)
   call new_potential(pot, mol, bas, 1, .true.)

   energy = 0.0_wp
   if (present(qsh)) then
      wfn%qsh = reshape(qsh, [size(qsh), 1])
   else
      wfn%qsh = reshape(qat, [size(qat), 1])
   end if
   wfn%qat = reshape(qat, [size(qat), 1])

   call make_coulomb(coulomb, mol, present(qsh))
   call coulomb%update(mol, cache)
   call coulomb%get_potential(mol, cache, wfn, pot)
   call coulomb%get_energy(mol, cache, wfn, energy)

   call check(error, sum(energy), ref, thr=thr_)
   if (allocated(error)) then
      print*, ref, sum(energy), ref - sum(energy)
   end if

end subroutine test_generic


subroutine test_numpot(error, mol, qat, qsh, make_coulomb, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   !> Shell-resolved partial charges
   real(wp), intent(in), optional :: qsh(:)

   !> Factory to create new coulomb objects
   procedure(coulomb_maker) :: make_coulomb

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   integer :: iat, ish
   type(basis_type) :: bas
   class(coulomb_type), allocatable :: coulomb
   type(container_cache) :: cache
   type(wavefunction_type) :: wfn
   type(potential_type) :: pot
   real(wp), allocatable :: numvat(:), numvsh(:)
   real(wp) :: er(mol%nat), el(mol%nat)
   real(wp), parameter :: step = 1.0e-5_wp
   real(wp) :: thr_

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   if (present(qsh)) then
      wfn%qsh = reshape(qsh, [size(qsh), 1])
      bas%nsh = size(qsh)
      bas%nao = size(qsh)
   else
      wfn%qsh = reshape(qat, [size(qat), 1])
      bas%nsh = size(qat)
      bas%nao = size(qat)
   end if
   wfn%qat = reshape(qat, [size(qat), 1])

   ! Setup potential with dummy basis set
   call new_potential(pot, mol, bas, 1, .true.)

   ! Setup coulomb object
   call make_coulomb(coulomb, mol, present(qsh))
   call coulomb%update(mol, cache)

   ! Numerical atomic potential
   allocate(numvat(mol%nat), source=0.0_wp)
   do iat = 1, mol%nat
      er = 0.0_wp
      el = 0.0_wp
      ! Right hand side
      wfn%qat(iat, 1) = wfn%qat(iat, 1) + step
      call coulomb%get_energy(mol, cache, wfn, er)
      ! Left hand side
      wfn%qat(iat, 1) = wfn%qat(iat, 1) - 2*step
      call coulomb%get_energy(mol, cache, wfn, el)

      wfn%qat(iat, 1) = wfn%qat(iat, 1) + step
      numvat(iat) = 0.5_wp*(sum(er) - sum(el))/step
   end do

   if (present(qsh)) then
      ! Numerical shell potential
      allocate(numvsh(bas%nsh), source=0.0_wp)   
      do ish = 1, bas%nsh
         er = 0.0_wp
         el = 0.0_wp
         ! Right hand side
         wfn%qsh(ish, 1) = wfn%qsh(ish, 1) + step
         call coulomb%get_energy(mol, cache, wfn, er)
         ! Left hand side
         wfn%qsh(ish, 1) = wfn%qsh(ish, 1) - 2*step
         call coulomb%get_energy(mol, cache, wfn, el)

         wfn%qsh(ish, 1) = wfn%qsh(ish, 1) + step
         numvsh(ish) = 0.5_wp*(sum(er) - sum(el))/step
      end do
   end if

   ! Analytic potentials
   call pot%reset()
   call coulomb%get_potential(mol, cache, wfn, pot)

   if (any(abs(pot%vat(:, 1) - numvat) > thr_)) then
      call test_failed(error, "Atom-resolved potential does not match")
      write(*,*) "numerical potential:"
      print'(3es21.14)', numvat
      write(*,*) "analytical potential:"
      print'(3es21.14)', pot%vat(: ,1)
      write(*,*) "difference:"
      print'(3es21.14)', pot%vat(: ,1) - numvat
   end if

   if (present(qsh)) then
      if (any(abs(pot%vsh(:,1) - numvsh) > thr_)) then
         call test_failed(error, "Shell-resolved potential does not match")
         write(*,*) "numerical potential:"
         print'(3es21.14)', numvsh
         write(*,*) "analytical potential:"
         print'(3es21.14)', pot%vsh(: ,1)
         write(*,*) "difference:"
         print'(3es21.14)', pot%vsh(: ,1) - numvsh
      end if
   end if

end subroutine test_numpot


subroutine test_numpotgrad(error, mol, get_charges, make_coulomb, shell, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Procedure to provide atomic charges (and their gradients)
   procedure(charge_maker) :: get_charges

   !> Factory to create new onsite fourth-order objects
   procedure(coulomb_maker) :: make_coulomb

   !> Test shell-resolved
   logical, intent(in) :: shell

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   integer :: iat, ic
   type(basis_type) :: bas
   type(potential_type) :: potl, potr
   class(coulomb_type), allocatable :: coulomb
   type(container_cache) :: cache
   type(wavefunction_type) :: wfn
   real(wp), allocatable :: numpotgrad(:, :, :)
   integer, allocatable :: nshell(:)
   real(wp), parameter :: step = 1.0e-5_wp
   real(wp) :: thr_

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   call make_basis(bas, mol, 6)
   call new_wavefunction(wfn, mol%nat, bas%nsh, bas%nao, 1, 0.0_wp, .true.)
   call new_potential(potr, mol, bas, 1, .true.)
   call new_potential(potl, mol, bas, 1, .true.)
   call make_coulomb(coulomb, mol, shell)

   if (shell) then
      allocate(numpotgrad(3, mol%nat, bas%nsh), source=0.0_wp)
      allocate(nshell(mol%nat))
      nshell = bas%nsh_at
   else
      allocate(numpotgrad(3, mol%nat, mol%nat), source=0.0_wp)
   end if

   do iat = 1, mol%nat
      do ic = 1, 3
         call potr%reset
         call potl%reset
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call get_charges(wfn, mol, nshell, error)
         if (allocated(error)) return
         call coulomb%update(mol, cache)
         call coulomb%get_potential(mol, cache, wfn, potr)

         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call get_charges(wfn, mol, nshell, error)
         if (allocated(error)) return
         call coulomb%update(mol, cache)
         call coulomb%get_potential(mol, cache, wfn, potl)

         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         if (shell) then
            numpotgrad(ic, iat, :) = 0.5_wp*(potr%vsh(:, 1) - potl%vsh(:, 1))/step
         else
            numpotgrad(ic, iat, :) = 0.5_wp*(potr%vat(:, 1) - potl%vat(:, 1))/step
         end if
      end do
   end do

   call get_charges(wfn, mol, nshell, error)
   if (allocated(error)) return
   call coulomb%update(mol, cache)
   call coulomb%get_potential(mol, cache, wfn, potl)
   call coulomb%get_potential_gradient(mol, cache, wfn, potl)

   if (shell) then
      if (any(abs(potl%dvshdr(:, :, :, 1) - numpotgrad) > thr_)) then
         call test_failed(error, "Gradient of shell-resolved potential does not match")
         print '(3es21.14)', potl%dvshdr(:, :, :, 1) - numpotgrad
      end if
   else
      if (any(abs(potl%dvatdr(:, :, :, 1) - numpotgrad) > thr_)) then
         call test_failed(error, "Gradient of atom-resolved potential does not match")
         print '(3es21.14)', potl%dvatdr(:, :, :, 1) - numpotgrad
      end if
   end if

end subroutine test_numpotgrad


subroutine test_numpotsigma(error, mol, get_charges, make_coulomb, shell, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Procedure to provide atomic charges (and their gradients)
   procedure(charge_maker) :: get_charges

   !> Factory to create new electrostatic objects
   procedure(coulomb_maker) :: make_coulomb

   !> Test shell-resolved
   logical, intent(in) :: shell

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   integer :: ic, jc
   type(basis_type) :: bas
   type(potential_type) :: potl, potr
   class(coulomb_type), allocatable :: coulomb
   type(container_cache) :: cache
   type(wavefunction_type) :: wfn
   real(wp) :: eps(3, 3)
   real(wp), allocatable :: numpotsigma(:, :, :)
   real(wp), allocatable :: xyz(:, :), lattice(:, :)
   integer, allocatable :: nshell(:)
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))
   real(wp), parameter :: step = 5.0e-5_wp
   real(wp) :: thr_

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   allocate(xyz(3, mol%nat), source=0.0_wp)

   call make_basis(bas, mol, 6)
   call new_wavefunction(wfn, mol%nat, bas%nsh, bas%nao, 1, 0.0_wp, .true.)
   call new_potential(potr, mol, bas, 1, .true.)
   call new_potential(potl, mol, bas, 1, .true.)
   call make_coulomb(coulomb, mol, shell)

   if (shell) then
      allocate(numpotsigma(3, 3, bas%nsh), source=0.0_wp)
      allocate(nshell(mol%nat))
      nshell = bas%nsh_at
   else
      allocate(numpotsigma(3, 3, mol%nat), source=0.0_wp)
   end if

   eps(:, :) = unity
   xyz(:, :) = mol%xyz
   if (any(mol%periodic)) lattice = mol%lattice
   do ic = 1, 3
      do jc = 1, 3
         call potr%reset
         call potl%reset
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         if (allocated(lattice)) mol%lattice(:, :) = matmul(eps, lattice)
         call get_charges(wfn, mol, nshell, error)
         if (allocated(error)) return
         call coulomb%update(mol, cache)
         call coulomb%get_potential(mol, cache, wfn, potr)

         eps(jc, ic) = eps(jc, ic) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         if (allocated(lattice)) mol%lattice(:, :) = matmul(eps, lattice)
         call get_charges(wfn, mol, nshell, error)
         if (allocated(error)) return
         call coulomb%update(mol, cache)
         call coulomb%get_potential(mol, cache, wfn, potl)

         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         if (allocated(lattice)) mol%lattice = lattice
         if (shell) then
            numpotsigma(jc, ic, :) = 0.5_wp*(potr%vsh(:, 1) - potl%vsh(:, 1))/step
         else
            numpotsigma(jc, ic, :) = 0.5_wp*(potr%vat(:, 1) - potl%vat(:, 1))/step
         end if
      end do
   end do

   call get_charges(wfn, mol, nshell, error)
   if (allocated(error)) return
   call coulomb%update(mol, cache)
   call coulomb%get_potential(mol, cache, wfn, potl)
   call coulomb%get_potential_gradient(mol, cache, wfn, potl)

   if (shell) then
      if (any(abs(potl%dvshdL(:, :, :, 1) - numpotsigma) > thr_)) then
         call test_failed(error, "Sigma of shell-resolved potential does not match")
         print '(3es21.14)', potl%dvshdL(:, :, :, 1) - numpotsigma
      end if
   else
      if (any(abs(potl%dvatdL(:, :, :, 1) - numpotsigma) > thr_)) then
         call test_failed(error, "Sigma of atom-resolved potential does not match")
         print '(3es21.14)', potl%dvatdL(:, :, :, 1) - numpotsigma
      end if
   end if

end subroutine test_numpotsigma


subroutine test_e_gxtb_m01_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [ &
      &  9.115222032717878E-1_wp, -7.820978364181586E-2_wp, -8.085043301444819E-1_wp, &
      & -1.370008894958090E-1_wp, -4.522670977383796E-1_wp,  2.366794200051806E-1_wp, &
      & -1.056480024872493E-1_wp, -7.542974658029267E-1_wp, -6.084417930861796E-1_wp, &
      &  3.252603520948062E-1_wp,  2.590209057852253E-1_wp, -3.245713596482476E-1_wp, &
      &  9.588528411034325E-2_wp,  7.104890489542338E-1_wp, -3.468331265256113E-1_wp, &
      &  1.076916633846319E+0_wp]

   call get_structure(mol, "MB16-43", "01")
   call test_generic(error, mol, qat, make_coulomb=make_coulomb_gxtb, &
      & ref=4.8976351076135145E-3_wp)

end subroutine test_e_gxtb_m01_atom

subroutine test_e_gxtb_m01_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [ &
      &  9.115222032717878E-1_wp, -7.820978364181586E-2_wp, -8.085043301444819E-1_wp, &
      & -1.370008894958090E-1_wp, -4.522670977383796E-1_wp,  2.366794200051806E-1_wp, &
      & -1.056480024872493E-1_wp, -7.542974658029267E-1_wp, -6.084417930861796E-1_wp, &
      &  3.252603520948062E-1_wp,  2.590209057852253E-1_wp, -3.245713596482476E-1_wp, &
      &  9.588528411034325E-2_wp,  7.104890489542338E-1_wp, -3.468331265256113E-1_wp, &
      &  1.076916633846319E+0_wp]
   real(wp), parameter :: qsh(*) = [ &
      &  4.452001114394941E-1_wp,  4.663220918322937E-1_wp, -7.820978364181586E-2_wp, &
      & -1.925657752689702E-1_wp, -6.159385548755116E-1_wp, -1.370008894958090E-1_wp, &
      & -9.716756357792056E-2_wp, -3.550995341604590E-1_wp,  2.366794200051806E-1_wp, &
      & -1.056480024872493E-1_wp, -1.692417149620185E-1_wp, -5.850557508409082E-1_wp, &
      & -1.817505622195756E-1_wp, -4.266912308666040E-1_wp,  3.252603520948062E-1_wp, &
      &  2.590209057852253E-1_wp, -7.093896295382240E-2_wp, -2.932293779939634E-1_wp, &
      &  3.959698129953824E-2_wp,  1.473378889672090E-1_wp, -5.145260485686576E-2_wp, &
      &  3.544877620709267E-1_wp,  3.560012868833071E-1_wp, -6.786424039218564E-2_wp, &
      & -2.789688861334256E-1_wp,  9.044625085911038E-1_wp,  5.497811147662567E-1_wp, &
      & -3.773269895110409E-1_wp]

   call get_structure(mol, "MB16-43", "01")
   call test_generic(error, mol, qat, qsh, make_coulomb_gxtb, &
      & ref=1.985899050208160E-3_wp)

end subroutine test_e_gxtb_m01_shell

subroutine test_e_gxtb_m02_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [ &
      & -1.097288775394156E-1_wp, -4.735023909155955E-1_wp, -1.707389583649019E-1_wp, &
      & -6.433491759460934E-1_wp,  8.965254778844174E-1_wp,  3.351267282573080E-1_wp, &
      & -6.657210518576395E-2_wp, -5.047599279669956E-2_wp,  5.146699702222350E-1_wp, &
      &  3.455270064324552E-1_wp,  8.365308656920378E-2_wp,  7.840093894994569E-1_wp, &
      & -6.236807887791791E-1_wp, -7.429581513991801E-2_wp, -6.465871716055771E-2_wp, &
      & -6.825088374693604E-1_wp]

   call get_structure(mol, "MB16-43", "02")
   call test_generic(error, mol, qat, make_coulomb=make_coulomb_gxtb, &
      & ref=2.5677988121880326E-3_wp)

end subroutine test_e_gxtb_m02_atom

subroutine test_e_gxtb_m02_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [ &
      & -1.097288775394156E-1_wp, -4.735023909155955E-1_wp, -1.707389583649019E-1_wp, &
      & -6.433491759460934E-1_wp,  8.965254778844174E-1_wp,  3.351267282573080E-1_wp, &
      & -6.657210518576395E-2_wp, -5.047599279669956E-2_wp,  5.146699702222350E-1_wp, &
      &  3.455270064324552E-1_wp,  8.365308656920378E-2_wp,  7.840093894994569E-1_wp, &
      & -6.236807887791791E-1_wp, -7.429581513991801E-2_wp, -6.465871716055771E-2_wp, &
      & -6.825088374693604E-1_wp]
   real(wp), parameter :: qsh(*) = [ &
      & -1.097288775394156E-1_wp, -1.649240032771995E-1_wp, -4.263046913300617E-1_wp, &
      &  1.177263036916656E-1_wp,  8.265314281056635E-2_wp, -2.533921011754683E-1_wp, &
      & -1.316164635096304E-1_wp, -5.117327124364630E-1_wp,  6.670371585258130E-1_wp, &
      &  2.689765878292192E-1_wp, -3.948826847061487E-2_wp,  3.351267282573080E-1_wp, &
      & -6.657210518576395E-2_wp, -5.047599279669956E-2_wp,  1.895960251239481E-1_wp, &
      &  2.395294429005155E-1_wp,  8.554450219777143E-2_wp,  3.455270064324552E-1_wp, &
      &  1.534569746334520E-1_wp, -6.980388806424820E-2_wp,  3.986364504419448E-1_wp, &
      &  3.853729390575120E-1_wp, -1.606712300918751E-1_wp, -4.630095586873040E-1_wp, &
      & -7.429581513991801E-2_wp, -6.465871716055771E-2_wp, -2.013953379594855E-1_wp, &
      & -6.189894955055433E-1_wp,  1.378759959956683E-1_wp]

   call get_structure(mol, "MB16-43", "02")
   call test_generic(error, mol, qat, qsh, make_coulomb=make_coulomb_gxtb, &
      & ref=9.522043568493133E-4_wp)

end subroutine test_e_gxtb_m02_shell

subroutine test_e_gxtb_oxacb_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   ! Currently a regression test based on GFN2 charges
   real(wp), parameter :: qat(*) = [ &
      & 0.43142759318744994_wp,  0.43150967285357533_wp,  0.43146165224958366_wp, &
      & 0.43147407117646708_wp,  0.90726023139894296_wp,  0.90728336242768870_wp, &
      & 0.90737069516079072_wp,  0.90712544694131525_wp, -0.77945631396614101_wp, &
      &-0.78097308159910250_wp, -0.78002089184600454_wp, -0.77924355344945684_wp, &
      &-0.55876993965122157_wp, -0.55823164158343652_wp, -0.55927229556741787_wp, &
      &-0.55894500773303069_wp]

   call get_structure(mol, "X23", "oxacb")
   call test_generic(error, mol, qat, make_coulomb=make_coulomb_gxtb, &
      & ref=7.0782005419069314E-3_wp)

end subroutine test_e_gxtb_oxacb_atom

subroutine test_e_gxtb_oxacb_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   ! Currently a regression test based on GFN2 charges
   real(wp), parameter :: qat(*) = [ &
      & 0.43142759318744994_wp,  0.43150967285357533_wp,  0.43146165224958366_wp, &
      & 0.43147407117646708_wp,  0.90726023139894296_wp,  0.90728336242768870_wp, &
      & 0.90737069516079072_wp,  0.90712544694131525_wp, -0.77945631396614101_wp, &
      &-0.78097308159910250_wp, -0.78002089184600454_wp, -0.77924355344945684_wp, &
      &-0.55876993965122157_wp, -0.55823164158343652_wp, -0.55927229556741787_wp, &
      &-0.55894500773303069_wp]
   real(wp), parameter :: qsh(*) = [ &
      & 0.43142759318744994_wp,  0.43150967285357533_wp,  0.43146165224958366_wp, &
      & 0.43147407117646708_wp,  0.12631625958529813_wp,  0.78094397181364483_wp, &
      & 0.12627921822540056_wp,  0.78100414420228814_wp,  0.12624063356686954_wp, &
      & 0.78113006159392118_wp,  0.12626556398329924_wp,  0.78085988295801601_wp, &
      & 0.19373236130795646_wp, -0.97318867527409747_wp,  0.19324264343503161_wp, &
      &-0.97421572503413412_wp,  0.19355130642179486_wp, -0.97357219826779939_wp, &
      & 0.19382295575165664_wp, -0.97306650920111348_wp,  0.23748355725433878_wp, &
      &-0.79625349690556035_wp,  0.23769535856264201_wp, -0.79592700014607853_wp, &
      & 0.23731240854425084_wp, -0.79658470411166871_wp,  0.23739675443019159_wp, &
      &-0.79634176216322228_wp]

   call get_structure(mol, "X23", "oxacb")
   call test_generic(error, mol, qat, qsh, make_coulomb_gxtb, &
      & ref=1.1775374572754271E-2_wp)

end subroutine test_e_gxtb_oxacb_shell

subroutine test_e_gxtb_oxacb_sc(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   ! Currently a regression test based on GFN2 charges
   real(wp), parameter :: qat1(*) = [ &
      & 3.41731844312030E-1_wp,  3.41716020106239E-1_wp,  3.41730526585671E-1_wp, &
      & 3.41714427217954E-1_wp,  3.80996046757999E-1_wp,  3.80989821246195E-1_wp, &
      & 3.81000747720282E-1_wp,  3.80990494183703E-1_wp, -3.70406587264474E-1_wp, &
      &-3.70407565207006E-1_wp, -3.70417590212352E-1_wp, -3.70399716470705E-1_wp, &
      &-3.52322260586075E-1_wp, -3.52304269439196E-1_wp, -3.52313440903261E-1_wp, &
      &-3.52298498047004E-1_wp]
   integer, parameter :: supercell(*) = [2, 2, 2]
   real(wp), parameter :: qat(*) = [spread(qat1, 2, product(supercell))]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "X23", "oxacb")
   call make_supercell(mol, supercell)
   call test_generic(error, mol, qat, qsh, make_coulomb_gxtb, &
      ref=3.3089982933808496E-3_wp, thr_in=1.0e-7_wp)

end subroutine test_e_gxtb_oxacb_sc

subroutine make_supercell(mol, rep)
   type(structure_type), intent(inout) :: mol
   integer, intent(in) :: rep(3)

   real(wp), allocatable :: xyz(:, :), lattice(:, :)
   integer, allocatable :: num(:)
   integer :: i, j, k, c

   num = reshape(spread([mol%num(mol%id)], 2, product(rep)), [product(rep)*mol%nat])
   lattice = reshape(&
      & [rep(1)*mol%lattice(:, 1), rep(2)*mol%lattice(:, 2), rep(3)*mol%lattice(:, 3)], &
      & shape(mol%lattice))
   allocate(xyz(3, product(rep)*mol%nat))
   c = 0
   do i = 0, rep(1) - 1
      do j = 0, rep(2) - 1
         do k = 0, rep(3) - 1
            xyz(:, c+1:c+mol%nat) = mol%xyz &
               & + spread(matmul(mol%lattice, [real(wp):: i, j, k]), 2, mol%nat)
            c = c + mol%nat
         end do
      end do
   end do

   call new(mol, num, xyz, lattice=lattice)
end subroutine make_supercell


subroutine test_p_gxtb_m01_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [ &
      &  9.115222032717878E-1_wp, -7.820978364181586E-2_wp, -8.085043301444819E-1_wp, &
      & -1.370008894958090E-1_wp, -4.522670977383796E-1_wp,  2.366794200051806E-1_wp, &
      & -1.056480024872493E-1_wp, -7.542974658029267E-1_wp, -6.084417930861796E-1_wp, &
      &  3.252603520948062E-1_wp,  2.590209057852253E-1_wp, -3.245713596482476E-1_wp, &
      &  9.588528411034325E-2_wp,  7.104890489542338E-1_wp, -3.468331265256113E-1_wp, &
      &  1.076916633846319E+0_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "MB16-43", "01")
   call test_numpot(error, mol, qat, qsh, make_coulomb_gxtb, thr_in=thr1)

end subroutine test_p_gxtb_m01_atom

subroutine test_p_gxtb_m01_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [ &
      &  9.115222032717878E-1_wp, -7.820978364181586E-2_wp, -8.085043301444819E-1_wp, &
      & -1.370008894958090E-1_wp, -4.522670977383796E-1_wp,  2.366794200051806E-1_wp, &
      & -1.056480024872493E-1_wp, -7.542974658029267E-1_wp, -6.084417930861796E-1_wp, &
      &  3.252603520948062E-1_wp,  2.590209057852253E-1_wp, -3.245713596482476E-1_wp, &
      &  9.588528411034325E-2_wp,  7.104890489542338E-1_wp, -3.468331265256113E-1_wp, &
      &  1.076916633846319E+0_wp]
   real(wp), parameter :: qsh(*) = [ &
      &  4.452001114394941E-1_wp,  4.663220918322937E-1_wp, -7.820978364181586E-2_wp, &
      & -1.925657752689702E-1_wp, -6.159385548755116E-1_wp, -1.370008894958090E-1_wp, &
      & -9.716756357792056E-2_wp, -3.550995341604590E-1_wp,  2.366794200051806E-1_wp, &
      & -1.056480024872493E-1_wp, -1.692417149620185E-1_wp, -5.850557508409082E-1_wp, &
      & -1.817505622195756E-1_wp, -4.266912308666040E-1_wp,  3.252603520948062E-1_wp, &
      &  2.590209057852253E-1_wp, -7.093896295382240E-2_wp, -2.932293779939634E-1_wp, &
      &  3.959698129953824E-2_wp,  1.473378889672090E-1_wp, -5.145260485686576E-2_wp, &
      &  3.544877620709267E-1_wp,  3.560012868833071E-1_wp, -6.786424039218564E-2_wp, &
      & -2.789688861334256E-1_wp,  9.044625085911038E-1_wp,  5.497811147662567E-1_wp, &
      & -3.773269895110409E-1_wp]

   call get_structure(mol, "MB16-43", "01")
   call test_numpot(error, mol, qat, qsh, make_coulomb_gxtb, thr_in=thr1)

end subroutine test_p_gxtb_m01_shell

subroutine test_p_gxtb_m02_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [ &
      & -1.097288775394156E-1_wp, -4.735023909155955E-1_wp, -1.707389583649019E-1_wp, &
      & -6.433491759460934E-1_wp,  8.965254778844174E-1_wp,  3.351267282573080E-1_wp, &
      & -6.657210518576395E-2_wp, -5.047599279669956E-2_wp,  5.146699702222350E-1_wp, &
      &  3.455270064324552E-1_wp,  8.365308656920378E-2_wp,  7.840093894994569E-1_wp, &
      & -6.236807887791791E-1_wp, -7.429581513991801E-2_wp, -6.465871716055771E-2_wp, &
      & -6.825088374693604E-1_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "MB16-43", "02")
   call test_numpot(error, mol, qat, qsh, make_coulomb_gxtb, thr_in=thr1)

end subroutine test_p_gxtb_m02_atom

subroutine test_p_gxtb_m02_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [ &
      & -1.097288775394156E-1_wp, -4.735023909155955E-1_wp, -1.707389583649019E-1_wp, &
      & -6.433491759460934E-1_wp,  8.965254778844174E-1_wp,  3.351267282573080E-1_wp, &
      & -6.657210518576395E-2_wp, -5.047599279669956E-2_wp,  5.146699702222350E-1_wp, &
      &  3.455270064324552E-1_wp,  8.365308656920378E-2_wp,  7.840093894994569E-1_wp, &
      & -6.236807887791791E-1_wp, -7.429581513991801E-2_wp, -6.465871716055771E-2_wp, &
      & -6.825088374693604E-1_wp]
   real(wp), parameter :: qsh(*) = [ &
      & -1.097288775394156E-1_wp, -1.649240032771995E-1_wp, -4.263046913300617E-1_wp, &
      &  1.177263036916656E-1_wp,  8.265314281056635E-2_wp, -2.533921011754683E-1_wp, &
      & -1.316164635096304E-1_wp, -5.117327124364630E-1_wp,  6.670371585258130E-1_wp, &
      &  2.689765878292192E-1_wp, -3.948826847061487E-2_wp,  3.351267282573080E-1_wp, &
      & -6.657210518576395E-2_wp, -5.047599279669956E-2_wp,  1.895960251239481E-1_wp, &
      &  2.395294429005155E-1_wp,  8.554450219777143E-2_wp,  3.455270064324552E-1_wp, &
      &  1.534569746334520E-1_wp, -6.980388806424820E-2_wp,  3.986364504419448E-1_wp, &
      &  3.853729390575120E-1_wp, -1.606712300918751E-1_wp, -4.630095586873040E-1_wp, &
      & -7.429581513991801E-2_wp, -6.465871716055771E-2_wp, -2.013953379594855E-1_wp, &
      & -6.189894955055433E-1_wp,  1.378759959956683E-1_wp]

   call get_structure(mol, "MB16-43", "02")
   call test_numpot(error, mol, qat, qsh, make_coulomb_gxtb, thr_in=thr1)

end subroutine test_p_gxtb_m02_shell


subroutine test_pg_ceh_lih_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "LiH")
   call test_numpotgrad(error, mol, get_charges_effceh, make_coulomb_gxtb, &
      & .false., thr_in=thr1)

end subroutine test_pg_ceh_lih_atom

subroutine test_pg_ceh_m15_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "15")
   call test_numpotgrad(error, mol, get_charges_effceh, make_coulomb_gxtb, &
      & .false., thr_in=thr1)

end subroutine test_pg_ceh_m15_atom

subroutine test_pg_ceh_m16_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "16")
   call test_numpotgrad(error, mol, get_charges_effceh, make_coulomb_gxtb, &
      & .true., thr_in=thr1)

end subroutine test_pg_ceh_m16_shell

subroutine test_ps_ceh_lih_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "LiH")
   call test_numpotsigma(error, mol, get_charges_effceh, make_coulomb_gxtb, &
      & .false., thr_in=thr1)

end subroutine test_ps_ceh_lih_atom

subroutine test_ps_ceh_co2_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "CO2")
   call test_numpotsigma(error, mol, get_charges_effceh, make_coulomb_gxtb, &
      & .false., thr_in=thr1)

end subroutine test_ps_ceh_co2_atom

subroutine test_ps_ceh_m05_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "05")
   call test_numpotsigma(error, mol, get_charges_effceh, make_coulomb_gxtb, &
      & .false., thr_in=thr1)

end subroutine test_ps_ceh_m05_atom

subroutine test_ps_ceh_m17_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "17")
   call test_numpotsigma(error, mol, get_charges_effceh, make_coulomb_gxtb, &
      & .true., thr_in=thr1)

end subroutine test_ps_ceh_m17_shell

end module test_coulomb_fourthorder
