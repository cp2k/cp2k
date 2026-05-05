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

module test_solvation_cds
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mctc_io_convert, only : kcaltoau
   use mstore, only : get_structure
   use tblite_container, only : container_cache
   use tblite_scf_potential, only : potential_type
   use tblite_solvation_cds, only : cds_solvation, cds_input
   use tblite_solvation_data_cds, only : get_cds_param
   use tblite_solvation_data, only : solvent_data, get_solvent_data, & 
      & get_vdw_rad_cosmo, get_vdw_rad_bondi, get_vdw_rad_d3
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   public :: collect_solvation_cds

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

   real(wp), parameter :: tension_water(20) = 1.0e-5_wp * [&
      &-0.08499967_wp, 0.46780225_wp,-2.87013596_wp,-3.95935069_wp,-0.29783987_wp, &
      &-0.48323273_wp, 0.00133622_wp, 0.20448945_wp, 0.20150600_wp, 0.36379863_wp, &
      &-3.47082133_wp,-0.93451053_wp,-1.46342018_wp,-0.32774697_wp,-0.38015204_wp, &
      &-0.35311116_wp,-0.19972593_wp,-0.12891363_wp,-1.19450558_wp,-1.61289300_wp]
   real(wp), parameter :: hbond_water(20) = -kcaltoau * [&
      & 6.70894947_wp, 0.00000000_wp, 0.00000000_wp, 0.00000000_wp, 0.00000000_wp, &
      & 1.26459036_wp, 3.52206160_wp, 2.30440543_wp, 1.98829409_wp, 0.00000000_wp, &
      & 0.00000000_wp, 0.00000000_wp, 0.00000000_wp, 0.00000000_wp, 2.68116653_wp, &
      & 0.38262428_wp, 1.02948365_wp, 0.00000000_wp, 0.00000000_wp, 0.00000000_wp]**2

contains


!> Collect all exported unit tests
subroutine collect_solvation_cds(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("sasa-e", test_e_sasa), &
      new_unittest("sasa-g", test_g_sasa), &
      new_unittest("sasa-p", test_p_sasa), &
      new_unittest("energy-neutral-alpb-gfn1", test_e_alpb_gfn1_all_solvents), &
      new_unittest("energy-neutral-alpb-gfn2", test_e_alpb_gfn2_all_solvents), &
      new_unittest("energy-neutral-gbsa-gfn1", test_e_gbsa_gfn1_all_solvents), &
      new_unittest("energy-neutral-gbsa-gfn2", test_e_gbsa_gfn2_all_solvents), &
      new_unittest("energy-charged-alpb-gfn1", test_e_charged_alpb_gfn1), &
      new_unittest("energy-charged-alpb-gfn2", test_e_charged_alpb_gfn2), &
      new_unittest("energy-charged-gbsa-gfn1", test_e_charged_gbsa_gfn1), &
      new_unittest("energy-charged-gbsa-gfn2", test_e_charged_gbsa_gfn2), &
      new_unittest("gradient-scf", test_g_cds), &
      new_unittest("gradient-nonscf", test_g_cds_nonscf), &
      new_unittest("potential", test_p_cds), &
      new_unittest("unsupported-solvent", test_unsupported_solvent, should_fail=.true.) &
      ]

end subroutine collect_solvation_cds


subroutine test_e(error, mol, input, qat, ref, method)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Solvation model input
   type(cds_input), intent(in) :: input

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   !> Reference energy
   real(wp), intent(in) :: ref

   !> Method for parameter selection
   character(len=*), intent(in), optional :: method

   type(cds_solvation) :: solv
   type(cds_input), allocatable :: scratch_input
   type(wavefunction_type) :: wfn
   type(potential_type) :: pot
   type(container_cache) :: cache
   real(wp) :: energy(mol%nat)

   energy = 0.0_wp
   wfn%qat = reshape(qat, [size(qat), 1])
   allocate(pot%vat(size(qat, 1), 1))

   scratch_input = input

   if (allocated(input%solvent) .and. present(method)) then
      call get_cds_param(scratch_input, mol, method, error)
      if(allocated(error)) then
         call test_failed(error, "No CDS parameters found for the method/solvent")
         return
      end if
   end if

   solv = cds_solvation(mol, scratch_input, method)

   call solv%update(mol, cache)
   call solv%get_potential(mol, cache, wfn, pot)
   call solv%get_energy(mol, cache, wfn, energy)
   !> CDS has some non-selfconsistent part as well:
   call solv%get_engrad(mol, cache, energy)

   if (abs(sum(energy) - ref) > thr) then
      call test_failed(error, "Energy does not match reference")
      print '(es21.14)', sum(energy)
      print '(es21.14)', ref
      print '(a)', "---"
   end if
end subroutine test_e


subroutine test_g(error, mol, input, qat, method)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Solvation model input
   type(cds_input), intent(in) :: input

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   !> Method for parameter selection
   character(len=*), intent(in), optional :: method

   type(cds_solvation) :: solv
   type(cds_input), allocatable :: scratch_input
   type(wavefunction_type) :: wfn
   type(potential_type) :: pot
   type(container_cache) :: cache
   real(wp), parameter :: step = 1.0e-4_wp
   real(wp), allocatable :: gradient(:, :), numg(:, :)
   real(wp) :: energy(mol%nat), er(mol%nat), el(mol%nat), sigma(3, 3)
   integer :: ii, ic

   wfn%qat = reshape(qat, [size(qat), 1])
   allocate(pot%vat(size(qat, 1), 1))

   scratch_input = input

   if (allocated(input%solvent) .and. present(method)) then   
      call get_cds_param(scratch_input, mol, method, error)
      if(allocated(error)) then
         call test_failed(error, "No CDS parameters found for the method/solvent")
         return
      end if
   end if

   solv = cds_solvation(mol, scratch_input, method)

   allocate(numg(3, mol%nat), gradient(3, mol%nat))
   do ii = 1, mol%nat
      do ic = 1, 3
         er = 0.0_wp
         el = 0.0_wp
         mol%xyz(ic, ii) = mol%xyz(ic, ii) + step
         call solv%update(mol, cache)
         call solv%get_potential(mol, cache, wfn, pot)
         call solv%get_energy(mol, cache, wfn, er)

         mol%xyz(ic, ii) = mol%xyz(ic, ii) - 2*step
         call solv%update(mol, cache)
         call solv%get_potential(mol, cache, wfn, pot)
         call solv%get_energy(mol, cache, wfn, el)

         mol%xyz(ic, ii) = mol%xyz(ic, ii) + step
         numg(ic, ii) = 0.5_wp*(sum(er) - sum(el))/step
      end do
   end do

   energy = 0.0_wp
   gradient(:, :) = 0.0_wp

   call solv%update(mol, cache)
   call solv%get_potential(mol, cache, wfn, pot)
   call solv%get_energy(mol, cache, wfn, energy)
   call solv%get_gradient(mol, cache, wfn, gradient, sigma)

   if (any(abs(gradient - numg) > thr2)) then
      call test_failed(error, "Gradient does not match")
      print '(3es20.13)', gradient
      print '(a)', "---"
      print '(3es20.13)', numg
      print '(a)', "---"
      print '(3es20.13)', gradient - numg
   end if
end subroutine test_g


subroutine test_g_nonscf(error, mol, input, method)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Solvation model input
   type(cds_input), intent(in) :: input

   !> Method for parameter selection
   character(len=*), intent(in), optional :: method

   type(cds_solvation) :: solv
   type(cds_input), allocatable :: scratch_input
   type(container_cache) :: cache
   real(wp), parameter :: step = 1.0e-4_wp
   real(wp), allocatable :: gradient(:, :), numg(:, :)
   real(wp) :: energy(mol%nat), er(mol%nat), el(mol%nat)
   integer :: ii, ic

   scratch_input = input

   if (allocated(input%solvent) .and. present(method)) then   
      call get_cds_param(scratch_input, mol, method, error)
      if(allocated(error)) then
         call test_failed(error, "No CDS parameters found for the method/solvent")
         return
      end if
   end if

   solv = cds_solvation(mol, scratch_input, method)

   allocate(numg(3, mol%nat), gradient(3, mol%nat))
   do ii = 1, mol%nat
      do ic = 1, 3
         er = 0.0_wp
         el = 0.0_wp
         mol%xyz(ic, ii) = mol%xyz(ic, ii) + step
         call solv%update(mol, cache)
         call solv%get_engrad(mol, cache, er)

         mol%xyz(ic, ii) = mol%xyz(ic, ii) - 2*step
         call solv%update(mol, cache)
         call solv%get_engrad(mol, cache, el)

         mol%xyz(ic, ii) = mol%xyz(ic, ii) + step
         numg(ic, ii) = 0.5_wp*(sum(er) - sum(el))/step
      end do
   end do

   energy = 0.0_wp
   gradient(:, :) = 0.0_wp

   call solv%update(mol, cache)
   call solv%get_engrad(mol, cache, energy, gradient)

   if (any(abs(gradient - numg) > thr2)) then
      call test_failed(error, "Gradient does not match")
      print '(3es20.13)', gradient
      print '(a)', "---"
      print '(3es20.13)', numg
      print '(a)', "---"
      print '(3es20.13)', gradient - numg
   end if
end subroutine test_g_nonscf

subroutine test_p(error, mol, input, qat, method)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Solvation model input
   type(cds_input), intent(in) :: input

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   !> Method for parameter selection
   character(len=*), intent(in), optional :: method

   type(cds_solvation) :: solv
   type(cds_input), allocatable :: scratch_input
   type(wavefunction_type) :: wfn
   type(potential_type) :: pot
   type(container_cache) :: cache
   real(wp), parameter :: step = 1.0e-4_wp
   real(wp), allocatable :: vat(:)
   real(wp) :: energy(mol%nat), er(mol%nat), el(mol%nat)
   integer :: ii

   wfn%qat = reshape(qat, [size(qat), 1])
   allocate(pot%vat(size(qat, 1), 1))

   scratch_input = input

   if (allocated(input%solvent) .and. present(method)) then   
      call get_cds_param(scratch_input, mol, method, error)
      if(allocated(error)) then
         call test_failed(error, "No CDS parameters found for the method/solvent")
         return
      end if
   end if

   solv = cds_solvation(mol, scratch_input, method)

   call solv%update(mol, cache)

   allocate(vat(mol%nat))
   do ii = 1, mol%nat
      er = 0.0_wp
      el = 0.0_wp
      wfn%qat(ii, 1) = wfn%qat(ii, 1) + step
      call solv%get_potential(mol, cache, wfn, pot)
      call solv%get_energy(mol, cache, wfn, er)

      wfn%qat(ii, 1) = wfn%qat(ii, 1) - 2*step
      call solv%get_potential(mol, cache, wfn, pot)
      call solv%get_energy(mol, cache, wfn, el)

      wfn%qat(ii, 1) = wfn%qat(ii, 1) + step
      vat(ii) = 0.5_wp*(sum(er) - sum(el))/step
   end do

   energy = 0.0_wp
   pot%vat(:, :) = 0.0_wp
   call solv%get_potential(mol, cache, wfn, pot)
   call solv%get_energy(mol, cache, wfn, energy)

   if (any(abs([pot%vat] - vat) > thr)) then
      call test_failed(error, "Potential does not match")
      print '(3es20.13)', pot%vat
      print '(a)', "---"
      print '(3es20.13)', vat
      print '(a)', "---"
      print '(3es20.13)', [pot%vat] - vat
   end if
end subroutine test_p


subroutine test_e_sasa(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(cds_input) :: input
   real(wp), parameter :: qat(*) = [&
      &-8.99890404486076E-2_wp, 9.42168087556583E-2_wp,-1.49387631509499E-1_wp, &
      &-2.99114121895542E-1_wp, 4.85527734875224E-1_wp,-6.83156326406137E-2_wp, &
      & 1.50011293889337E-2_wp, 2.79368544459465E-1_wp,-1.24072452878322E-1_wp, &
      &-9.36760994051244E-2_wp,-2.19062123031622E-1_wp, 2.14538817685587E-1_wp, &
      & 3.06156726072831E-1_wp,-3.86105514712244E-1_wp,-1.51265171389388E-3_wp, &
      & 3.64255069977693E-2_wp]
   real(wp), allocatable :: rad(:), tension(:), hbond(:)

   call get_structure(mol, "MB16-43", "04")
   rad = get_vdw_rad_cosmo(mol%num)
   tension = tension_water(mol%num)
   hbond = hbond_water(mol%num)

   input = cds_input(probe=0.3_wp, nang=110, rad=rad, tension=tension)
   call test_e(error, mol, input, qat, -2.0260929722303264E-3_wp)

   if (allocated(error)) return

   input = cds_input(probe=0.3_wp, nang=110, rad=rad, tension=tension, hbond=hbond)
   call test_e(error, mol, input, qat, -3.89595604489100E-3_wp)

end subroutine test_e_sasa


subroutine test_g_sasa(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(cds_input) :: input
   real(wp), parameter :: qat(*) = [&
      & 2.08159387594211E-1_wp,-3.78010519998818E-1_wp, 3.36498247356244E-2_wp, &
      &-4.11556158912895E-1_wp, 8.14928196660512E-2_wp,-2.00886649303053E-1_wp, &
      & 2.44756994282684E-1_wp, 2.54580499189089E-2_wp, 2.59835128092562E-1_wp, &
      & 4.21683321877209E-1_wp, 1.37097163086023E-1_wp, 4.06951664942900E-2_wp, &
      &-1.10955378625897E-1_wp,-6.44033540918074E-2_wp,-1.91525919028143E-1_wp, &
      &-9.54898757869102E-2_wp]
   real(wp), allocatable :: rad(:), tension(:), hbond(:)

   call get_structure(mol, "MB16-43", "06")
   rad = get_vdw_rad_bondi(mol%num)
   tension = tension_water(mol%num)
   hbond = hbond_water(mol%num)

   input = cds_input(probe=0.3_wp, nang=110, rad=rad, tension=tension, hbond=hbond)
   call test_g(error, mol, input, qat)

   if (allocated(error)) return

   call test_g_nonscf(error, mol, input)

end subroutine test_g_sasa


subroutine test_p_sasa(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(cds_input) :: input
   real(wp), parameter :: qat(*) = [&
      &-2.05668345919710E-1_wp,-3.99553123071811E-1_wp, 3.29242774348191E-1_wp, &
      &-3.11737933844111E-1_wp, 3.58851882478133E-2_wp, 3.21886835736497E-1_wp, &
      & 4.14743455841314E-2_wp, 2.95727359478547E-2_wp,-5.06347224522431E-1_wp, &
      & 3.43067182413129E-1_wp, 6.88373767679515E-1_wp, 7.03357390141253E-2_wp, &
      &-9.62424552888750E-2_wp,-1.32209348056625E-1_wp, 9.78998441832186E-2_wp, &
      &-3.05979982450903E-1_wp]
   real(wp), allocatable :: rad(:), tension(:), hbond(:)

   call get_structure(mol, "MB16-43", "08")
   rad = get_vdw_rad_d3(mol%num)
   tension = tension_water(mol%num)
   hbond = hbond_water(mol%num)

   input = cds_input(probe=2.2_wp, nang=110, rad=rad, tension=tension, hbond=hbond)
   call test_p(error, mol, input, qat)

end subroutine test_p_sasa


subroutine test_e_alpb_gfn1_all_solvents(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(solvent_data) :: solvent
   real(wp), parameter :: qat(*) = [&
      &-8.99890404486076E-2_wp, 9.42168087556583E-2_wp,-1.49387631509499E-1_wp, &
      &-2.99114121895542E-1_wp, 4.85527734875224E-1_wp,-6.83156326406137E-2_wp, &
      & 1.50011293889337E-2_wp, 2.79368544459465E-1_wp,-1.24072452878322E-1_wp, &
      &-9.36760994051244E-2_wp,-2.19062123031622E-1_wp, 2.14538817685587E-1_wp, &
      & 3.06156726072831E-1_wp,-3.86105514712244E-1_wp,-1.51265171389388E-3_wp, &
      & 3.64255069977693E-2_wp]
   type(cds_input) :: input
   integer, parameter :: nsolvents = 25
   character(len=*), parameter :: solvents(*) = [character(len=nsolvents):: &
      & "acetone", "acetonitrile", "aniline", "benzaldehyde", "benzene", &
      & "ch2cl2", "chcl3", "cs2", "dioxane", "dmf", "dmso", "ethanol", &
      & "ether", "ethylacetate", "furane", "hexadecane", "hexane", &
      & "nitromethane", "methanol", "octanol", "phenol", "thf", "toluene", &
      & "water", "woctanol"] 
   real(wp), parameter :: refs(*) = [&
      &-1.87575467176401E-2_wp, -2.19756935070354E-2_wp, -7.59897506913909E-3_wp, &
      &-1.28403478505020E-2_wp, -1.39686149805164E-2_wp, -1.08664447124936E-2_wp, &
      &-1.04964421005823E-2_wp, -1.14057753202729E-2_wp, -7.06543816921921E-3_wp, &
      &-1.44521134955677E-2_wp, -1.44911028376612E-2_wp, -8.97886030389678E-3_wp, &
      &-1.95763354166587E-2_wp, -1.14400229683406E-2_wp, -8.60633510983436E-3_wp, &
      &-1.12789113389271E-2_wp, -1.11505748030504E-2_wp, -8.81015006824253E-3_wp, &
      &-9.43366007096800E-3_wp, -9.76784539429431E-3_wp, -5.54432225159284E-3_wp, &
      &-1.20189216613106E-2_wp, -1.09141276818639E-2_wp, -4.39565407359322E-3_wp, &
      &-1.03868393814325E-2_wp]
   integer :: i

   call get_structure(mol, "MB16-43", "04")

   ! Check GFN1/ALPB for all available solvents
   do i = 1, nsolvents
      solvent = get_solvent_data(solvents(i))
      input = cds_input(solvent=solvent%solvent, alpb=.true.)
      call test_e(error, mol, input, qat, refs(i), method='gfn1') 
      if (allocated(error)) return
   end do 

end subroutine test_e_alpb_gfn1_all_solvents

subroutine test_e_alpb_gfn2_all_solvents(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(solvent_data) :: solvent
   real(wp), parameter :: qat(*) = [&
      & 2.29208822115185E-1_wp, 4.70816658242009E-2_wp,-3.52834459119718E-2_wp, &
      & 3.09012847396269E-2_wp, 3.28898468854920E-1_wp,-2.01747405019535E-1_wp, &
      & 5.46554362391008E-2_wp,-1.09681283064574E-1_wp,-3.47340505091849E-1_wp, &
      & 1.84567817865267E-1_wp,-2.07552337277185E-1_wp, 4.67140380802351E-1_wp, &
      &-1.84261319200178E-2_wp,-1.05015595324833E-1_wp, 6.52511545312054E-2_wp, &
      &-3.82658324740237E-1_wp]
   type(cds_input) :: input
   integer, parameter :: nsolvents = 25
   character(len=*), parameter :: solvents(*) = [character(len=nsolvents):: &
      & "acetone", "acetonitrile", "aniline", "benzaldehyde", "benzene", &
      & "ch2cl2", "chcl3", "cs2", "dioxane", "dmf", "dmso", "ethanol", &
      & "ether", "ethylacetate", "furane", "hexadecane", "hexane", &
      & "nitromethane", "methanol", "octanol", "phenol", "thf", "toluene", &
      & "water", "woctanol"] 
   real(wp), parameter :: refs(*) = [&
      &-2.05399633856027E-2_wp, -1.76286670133248E-2_wp, -1.52489808627813E-2_wp, &
      &-1.61094315706730E-2_wp, -2.00524539538480E-2_wp, -2.04518654284223E-2_wp, &
      &-1.89495890022200E-2_wp, -1.83634224905815E-2_wp, -1.35338999058532E-2_wp, &
      &-2.02119670921433E-2_wp, -2.06860521914816E-2_wp, -1.55826272873075E-2_wp, &
      &-1.91024778367152E-2_wp, -1.64456391051830E-2_wp, -1.51588817057209E-2_wp, &
      &-7.85038297075149E-3_wp, -1.80332280208901E-2_wp, -1.27398470763579E-2_wp, &
      &-1.43153922624655E-2_wp,  1.17718049638770E-3_wp, -1.25432724254740E-2_wp, &
      &-2.00228656875909E-2_wp, -1.75195809653790E-2_wp,  1.88693640907118E-3_wp, &
      &-1.59803235935055E-2_wp]
   integer :: i

   call get_structure(mol, "MB16-43", "05")

   ! Check GFN2/ALPB for all available solvents
   do i = 1, nsolvents
      solvent = get_solvent_data(solvents(i))
      input = cds_input(solvent=solvent%solvent, alpb=.true.)
      call test_e(error, mol, input, qat, refs(i), method='gfn2') 
      if (allocated(error)) return
   end do 

end subroutine test_e_alpb_gfn2_all_solvents

subroutine test_e_gbsa_gfn1_all_solvents(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(solvent_data) :: solvent
   real(wp), parameter :: qat(*) = [&
      & 2.08159387594211E-1_wp,-3.78010519998818E-1_wp, 3.36498247356244E-2_wp, &
      &-4.11556158912895E-1_wp, 8.14928196660512E-2_wp,-2.00886649303053E-1_wp, &
      & 2.44756994282684E-1_wp, 2.54580499189089E-2_wp, 2.59835128092562E-1_wp, &
      & 4.21683321877209E-1_wp, 1.37097163086023E-1_wp, 4.06951664942900E-2_wp, &
      &-1.10955378625897E-1_wp,-6.44033540918074E-2_wp,-1.91525919028143E-1_wp, &
      &-9.54898757869102E-2_wp]
   type(cds_input) :: input
   integer, parameter :: nsolvents = 12
   character(len=*), parameter :: solvents(*) = [character(len=nsolvents):: &
      & "acetone", "acetonitrile", "benzene", "ch2cl2", "chcl3", "cs2", &
      & "dmso", "ether", "methanol", "thf", "toluene", "water"] 
   real(wp), parameter :: refs(*) = [&
      &-2.16720081179071E-3_wp, -3.28426923492598E-2_wp, -2.52277980102607E-2_wp, &
      &-2.41431963258825E-2_wp, -2.41431963258825E-2_wp, -4.54197252450546E-3_wp, &
      &-3.28426923492598E-2_wp, -2.62339810005497E-2_wp, -2.55870388212558E-2_wp, &
      &-2.95720285475705E-2_wp, -2.52277980102607E-2_wp, -1.87472790863901E-2_wp]
   integer :: i

   call get_structure(mol, "MB16-43", "06")

   ! Check GFN1/GBSA for all available solvents
   do i = 1, nsolvents
      solvent = get_solvent_data(solvents(i))
      input = cds_input(solvent=solvent%solvent, alpb=.false.)
      call test_e(error, mol, input, qat, refs(i), method='gfn1') 
      if (allocated(error)) return
   end do 

end subroutine test_e_gbsa_gfn1_all_solvents

subroutine test_e_gbsa_gfn2_all_solvents(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(solvent_data) :: solvent
   real(wp), parameter :: qat(*) = [&
      &-1.57321098180703E-1_wp, 1.65233008998668E-1_wp, 3.22320267782066E-1_wp, &
      & 3.63564544135336E-2_wp, 4.85639267214320E-2_wp,-3.59203277893926E-1_wp, &
      &-1.93841260011383E-1_wp,-3.86495230324447E-1_wp, 3.10104147485353E-1_wp, &
      & 8.34907519580185E-2_wp,-3.62672063405622E-1_wp, 3.64143595819311E-1_wp, &
      & 3.34640678947868E-1_wp,-4.69881543486815E-1_wp,-1.89222615863620E-1_wp, &
      & 4.53784257040286E-1_wp]
   type(cds_input) :: input
   integer, parameter :: nsolvents = 14
   character(len=*), parameter :: solvents(*) = [character(len=nsolvents):: &
      & "acetone", "acetonitrile", "benzene", "ch2cl2", "chcl3", "cs2", &
      & "dmf", "dmso", "ether", "hexane", "methanol", "thf", "toluene", &
      & "water"]
   real(wp), parameter :: refs(*) = [&
      &-2.21027923191889E-2_wp, -1.67263825829947E-2_wp, -1.54769936490801E-2_wp, &
      &-1.61275987597537E-2_wp, -1.46478122238807E-2_wp, -1.64141563713530E-2_wp, &
      &-2.22619239865685E-2_wp, -2.50545882957244E-2_wp, -1.93615578636444E-2_wp, &
      &-1.49825471291675E-2_wp, -1.73962890782651E-2_wp, -2.06110330162246E-2_wp, &
      &-1.54769936490801E-2_wp, -8.02197393839715E-3_wp]
   integer :: i

   call get_structure(mol, "MB16-43", "07")

   ! Check GFN2/GBSA for all available solvents
   do i = 1, nsolvents
      solvent = get_solvent_data(solvents(i))
      input = cds_input(solvent=solvent%solvent, alpb=.false.)
      call test_e(error, mol, input, qat, refs(i), method='gfn2') 
      if (allocated(error)) return
   end do 

end subroutine test_e_gbsa_gfn2_all_solvents


subroutine test_e_charged_alpb_gfn1(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(solvent_data) :: solvent
   type(cds_input) :: input
   real(wp), parameter :: qat(*) = [&
      & 2.91737167991831E-1_wp,-4.44921281646692E-1_wp, 7.03886898899541E-2_wp, &
      & 4.32458130745259E-2_wp, 3.90393832017139E-2_wp, 9.56766032864156E-2_wp, &
      & 6.01866392558221E-2_wp,-3.31751380042108E-1_wp, 1.56707781479786E-1_wp, &
      & 8.89443731456522E-2_wp,-6.03526873354187E-2_wp, 3.83058968439732E-2_wp, &
      & 9.51289165540499E-2_wp,-9.20860756307965E-2_wp, 7.12122392663975E-2_wp, &
      & 2.76742553014572E-1_wp,-4.38083248713717E-1_wp,-1.76798298145349E-1_wp, &
      & 2.14382387079637E-1_wp, 3.25215078814299E-1_wp,-4.30846325150764E-1_wp, &
      & 1.12908046155543E-1_wp, 6.61648434849400E-2_wp, 9.89171412174742E-2_wp, &
      & 5.80323177221752E-2_wp,-4.43388471029542E-1_wp, 3.38493755177292E-1_wp, &
      &-3.52969753592199E-1_wp, 4.59552645375000E-1_wp,-5.33442813398395E-1_wp, &
      &-5.23460418980519E-1_wp,-3.14158345734236E-1_wp, 8.63777453496141E-2_wp, &
      & 5.51896920278464E-2_wp, 5.60421076086023E-2_wp, 1.02875052304244E-1_wp, &
      & 5.43513415521492E-2_wp,-3.13791816210054E-1_wp, 1.68862983166011E-1_wp, &
      & 8.93759201291419E-2_wp,-8.90121909290432E-2_wp, 4.15150867426933E-2_wp, &
      & 1.22221651251480E-1_wp,-8.26447904844349E-2_wp, 1.00154406589010E-1_wp, &
      & 2.77572579099978E-1_wp,-4.29147369583175E-1_wp,-1.78581481555413E-1_wp, &
      & 2.09487890121871E-1_wp, 3.17109649407645E-1_wp,-4.62948099575570E-1_wp, &
      & 9.81620738022878E-2_wp, 5.14984224707990E-2_wp, 9.63222020737258E-2_wp, &
      & 3.80443799704811E-2_wp,-4.41189291092377E-1_wp, 3.13549746324888E-1_wp, &
      &-4.41335746902051E-1_wp, 3.01219329594079E-1_wp]

   call get_structure(mol, "UPU23", "0a")
   solvent = get_solvent_data("water")
   input = cds_input(solvent=solvent%solvent, alpb=.true.)
   call test_e(error, mol, input, qat, -1.74474048674407E-02_wp, method='gfn1')

end subroutine test_e_charged_alpb_gfn1

subroutine test_e_charged_alpb_gfn2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(solvent_data) :: solvent
   type(cds_input) :: input
   real(wp), parameter :: qat(*) = [&
      & 2.91737167991831E-1_wp,-4.44921281646692E-1_wp, 7.03886898899541E-2_wp, &
      & 4.32458130745259E-2_wp, 3.90393832017139E-2_wp, 9.56766032864156E-2_wp, &
      & 6.01866392558221E-2_wp,-3.31751380042108E-1_wp, 1.56707781479786E-1_wp, &
      & 8.89443731456522E-2_wp,-6.03526873354187E-2_wp, 3.83058968439732E-2_wp, &
      & 9.51289165540499E-2_wp,-9.20860756307965E-2_wp, 7.12122392663975E-2_wp, &
      & 2.76742553014572E-1_wp,-4.38083248713717E-1_wp,-1.76798298145349E-1_wp, &
      & 2.14382387079637E-1_wp, 3.25215078814299E-1_wp,-4.30846325150764E-1_wp, &
      & 1.12908046155543E-1_wp, 6.61648434849400E-2_wp, 9.89171412174742E-2_wp, &
      & 5.80323177221752E-2_wp,-4.43388471029542E-1_wp, 3.38493755177292E-1_wp, &
      &-3.52969753592199E-1_wp, 4.59552645375000E-1_wp,-5.33442813398395E-1_wp, &
      &-5.23460418980519E-1_wp,-3.14158345734236E-1_wp, 8.63777453496141E-2_wp, &
      & 5.51896920278464E-2_wp, 5.60421076086023E-2_wp, 1.02875052304244E-1_wp, &
      & 5.43513415521492E-2_wp,-3.13791816210054E-1_wp, 1.68862983166011E-1_wp, &
      & 8.93759201291419E-2_wp,-8.90121909290432E-2_wp, 4.15150867426933E-2_wp, &
      & 1.22221651251480E-1_wp,-8.26447904844349E-2_wp, 1.00154406589010E-1_wp, &
      & 2.77572579099978E-1_wp,-4.29147369583175E-1_wp,-1.78581481555413E-1_wp, &
      & 2.09487890121871E-1_wp, 3.17109649407645E-1_wp,-4.62948099575570E-1_wp, &
      & 9.81620738022878E-2_wp, 5.14984224707990E-2_wp, 9.63222020737258E-2_wp, &
      & 3.80443799704811E-2_wp,-4.41189291092377E-1_wp, 3.13549746324888E-1_wp, &
      &-4.41335746902051E-1_wp, 3.01219329594079E-1_wp]

   call get_structure(mol, "UPU23", "0a")
   solvent = get_solvent_data("water")
   input = cds_input(solvent=solvent%solvent, alpb=.true.)
   call test_e(error, mol, input, qat, -7.40555251354138E-03_wp, method='gfn2')

end subroutine test_e_charged_alpb_gfn2

subroutine test_e_charged_gbsa_gfn1(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(solvent_data) :: solvent
   type(cds_input) :: input
   real(wp), parameter :: qat(*) = [&
      & 2.91737167991831E-1_wp,-4.44921281646692E-1_wp, 7.03886898899541E-2_wp, &
      & 4.32458130745259E-2_wp, 3.90393832017139E-2_wp, 9.56766032864156E-2_wp, &
      & 6.01866392558221E-2_wp,-3.31751380042108E-1_wp, 1.56707781479786E-1_wp, &
      & 8.89443731456522E-2_wp,-6.03526873354187E-2_wp, 3.83058968439732E-2_wp, &
      & 9.51289165540499E-2_wp,-9.20860756307965E-2_wp, 7.12122392663975E-2_wp, &
      & 2.76742553014572E-1_wp,-4.38083248713717E-1_wp,-1.76798298145349E-1_wp, &
      & 2.14382387079637E-1_wp, 3.25215078814299E-1_wp,-4.30846325150764E-1_wp, &
      & 1.12908046155543E-1_wp, 6.61648434849400E-2_wp, 9.89171412174742E-2_wp, &
      & 5.80323177221752E-2_wp,-4.43388471029542E-1_wp, 3.38493755177292E-1_wp, &
      &-3.52969753592199E-1_wp, 4.59552645375000E-1_wp,-5.33442813398395E-1_wp, &
      &-5.23460418980519E-1_wp,-3.14158345734236E-1_wp, 8.63777453496141E-2_wp, &
      & 5.51896920278464E-2_wp, 5.60421076086023E-2_wp, 1.02875052304244E-1_wp, &
      & 5.43513415521492E-2_wp,-3.13791816210054E-1_wp, 1.68862983166011E-1_wp, &
      & 8.93759201291419E-2_wp,-8.90121909290432E-2_wp, 4.15150867426933E-2_wp, &
      & 1.22221651251480E-1_wp,-8.26447904844349E-2_wp, 1.00154406589010E-1_wp, &
      & 2.77572579099978E-1_wp,-4.29147369583175E-1_wp,-1.78581481555413E-1_wp, &
      & 2.09487890121871E-1_wp, 3.17109649407645E-1_wp,-4.62948099575570E-1_wp, &
      & 9.81620738022878E-2_wp, 5.14984224707990E-2_wp, 9.63222020737258E-2_wp, &
      & 3.80443799704811E-2_wp,-4.41189291092377E-1_wp, 3.13549746324888E-1_wp, &
      &-4.41335746902051E-1_wp, 3.01219329594079E-1_wp]

   call get_structure(mol, "UPU23", "0a")
   solvent = get_solvent_data("water")
   input = cds_input(solvent=solvent%solvent, alpb=.false.)
   call test_e(error, mol, input, qat, 8.78720929142896E-03_wp, method='gfn1')

end subroutine test_e_charged_gbsa_gfn1

subroutine test_e_charged_gbsa_gfn2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(solvent_data) :: solvent
   type(cds_input) :: input
   real(wp), parameter :: qat(*) = [&
      & 2.91737167991831E-1_wp,-4.44921281646692E-1_wp, 7.03886898899541E-2_wp, &
      & 4.32458130745259E-2_wp, 3.90393832017139E-2_wp, 9.56766032864156E-2_wp, &
      & 6.01866392558221E-2_wp,-3.31751380042108E-1_wp, 1.56707781479786E-1_wp, &
      & 8.89443731456522E-2_wp,-6.03526873354187E-2_wp, 3.83058968439732E-2_wp, &
      & 9.51289165540499E-2_wp,-9.20860756307965E-2_wp, 7.12122392663975E-2_wp, &
      & 2.76742553014572E-1_wp,-4.38083248713717E-1_wp,-1.76798298145349E-1_wp, &
      & 2.14382387079637E-1_wp, 3.25215078814299E-1_wp,-4.30846325150764E-1_wp, &
      & 1.12908046155543E-1_wp, 6.61648434849400E-2_wp, 9.89171412174742E-2_wp, &
      & 5.80323177221752E-2_wp,-4.43388471029542E-1_wp, 3.38493755177292E-1_wp, &
      &-3.52969753592199E-1_wp, 4.59552645375000E-1_wp,-5.33442813398395E-1_wp, &
      &-5.23460418980519E-1_wp,-3.14158345734236E-1_wp, 8.63777453496141E-2_wp, &
      & 5.51896920278464E-2_wp, 5.60421076086023E-2_wp, 1.02875052304244E-1_wp, &
      & 5.43513415521492E-2_wp,-3.13791816210054E-1_wp, 1.68862983166011E-1_wp, &
      & 8.93759201291419E-2_wp,-8.90121909290432E-2_wp, 4.15150867426933E-2_wp, &
      & 1.22221651251480E-1_wp,-8.26447904844349E-2_wp, 1.00154406589010E-1_wp, &
      & 2.77572579099978E-1_wp,-4.29147369583175E-1_wp,-1.78581481555413E-1_wp, &
      & 2.09487890121871E-1_wp, 3.17109649407645E-1_wp,-4.62948099575570E-1_wp, &
      & 9.81620738022878E-2_wp, 5.14984224707990E-2_wp, 9.63222020737258E-2_wp, &
      & 3.80443799704811E-2_wp,-4.41189291092377E-1_wp, 3.13549746324888E-1_wp, &
      &-4.41335746902051E-1_wp, 3.01219329594079E-1_wp]

   call get_structure(mol, "UPU23", "0a")
   solvent = get_solvent_data("water")
   input = cds_input(solvent=solvent%solvent, alpb=.false.)
   call test_e(error, mol, input, qat, -1.68750629525887E-02_wp, method='gfn2')

end subroutine test_e_charged_gbsa_gfn2


subroutine test_g_cds(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(solvent_data) :: solvent
   type(cds_input) :: input
   real(wp), parameter :: qat(*) = [&
      & 2.08159387594211E-1_wp,-3.78010519998818E-1_wp, 3.36498247356244E-2_wp, &
      &-4.11556158912895E-1_wp, 8.14928196660512E-2_wp,-2.00886649303053E-1_wp, &
      & 2.44756994282684E-1_wp, 2.54580499189089E-2_wp, 2.59835128092562E-1_wp, &
      & 4.21683321877209E-1_wp, 1.37097163086023E-1_wp, 4.06951664942900E-2_wp, &
      &-1.10955378625897E-1_wp,-6.44033540918074E-2_wp,-1.91525919028143E-1_wp, &
      &-9.54898757869102E-2_wp]

   call get_structure(mol, "MB16-43", "06")
   solvent = get_solvent_data("water")
   input = cds_input(solvent=solvent%solvent, alpb=.true.)
   call test_g(error, mol, input, qat, method='gfn2')

end subroutine test_g_cds


subroutine test_g_cds_nonscf(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(solvent_data) :: solvent
   type(cds_input) :: input

   call get_structure(mol, "MB16-43", "07")
   solvent = get_solvent_data("water")
   input = cds_input(solvent=solvent%solvent, alpb=.true.)
   call test_g_nonscf(error, mol, input, method='gfn2')

end subroutine test_g_cds_nonscf


subroutine test_p_cds(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(solvent_data) :: solvent
   type(cds_input) :: input
   real(wp), parameter :: qat(*) = [&
      &-2.05668345919710E-1_wp,-3.99553123071811E-1_wp, 3.29242774348191E-1_wp, &
      &-3.11737933844111E-1_wp, 3.58851882478133E-2_wp, 3.21886835736497E-1_wp, &
      & 4.14743455841314E-2_wp, 2.95727359478547E-2_wp,-5.06347224522431E-1_wp, &
      & 3.43067182413129E-1_wp, 6.88373767679515E-1_wp, 7.03357390141253E-2_wp, &
      &-9.62424552888750E-2_wp,-1.32209348056625E-1_wp, 9.78998441832186E-2_wp, &
      &-3.05979982450903E-1_wp]

   call get_structure(mol, "MB16-43", "08")
   solvent = get_solvent_data("water")
   input = cds_input(solvent=solvent%solvent, alpb=.true.)
   call test_p(error, mol, input, qat, method='gfn2')

end subroutine test_p_cds


subroutine test_unsupported_solvent(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(solvent_data) :: solvent
   real(wp), parameter :: qat(*) = [&
      &-8.99890404486076E-2_wp, 9.42168087556583E-2_wp,-1.49387631509499E-1_wp, &
      &-2.99114121895542E-1_wp, 4.85527734875224E-1_wp,-6.83156326406137E-2_wp, &
      & 1.50011293889337E-2_wp, 2.79368544459465E-1_wp,-1.24072452878322E-1_wp, &
      &-9.36760994051244E-2_wp,-2.19062123031622E-1_wp, 2.14538817685587E-1_wp, &
      & 3.06156726072831E-1_wp,-3.86105514712244E-1_wp,-1.51265171389388E-3_wp, &
      & 3.64255069977693E-2_wp]
   type(cds_input) :: input

   call get_structure(mol, "MB16-43", "04")

   ! Check GFN1/GBSA unavailable solvent
   solvent = get_solvent_data("aniline")
   input = cds_input(solvent=solvent%solvent, alpb=.false.)
   call get_cds_param(input, mol, 'gfn1', error)

end subroutine test_unsupported_solvent


end module test_solvation_cds
