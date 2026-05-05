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

module test_solvation_shift
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   use tblite_container, only : container_cache
   use tblite_scf_potential, only : potential_type
   use tblite_solvation_shift, only : shift_solvation, shift_input, solution_state
   use tblite_solvation_data_shift, only : get_shift_param
   use tblite_solvation_data, only : solvent_data, get_solvent_data
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   public :: collect_solvation_shift


   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

contains


!> Collect all exported unit tests
subroutine collect_solvation_shift(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("solvation-shift-alpb-gfn1", test_e_alpb_gfn1_all_solvents), &
      new_unittest("solvation-shift-alpb-gfn2", test_e_alpb_gfn2_all_solvents), &
      new_unittest("solvation-shift-gbsa-gfn1", test_e_gbsa_gfn1_all_solvents), &
      new_unittest("solvation-shift-gbsa-gfn2", test_e_gbsa_gfn2_all_solvents), &
      new_unittest("unsupported-solvent", test_unsupported_solvent, should_fail=.true.) &
      ]

end subroutine collect_solvation_shift


subroutine test_e(error, mol, input, ref, method)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Solvation model input
   type(shift_input), intent(in) :: input

   !> Reference energy
   real(wp), intent(in) :: ref

   !> Method for parameter selection
   character(len=*), intent(in), optional :: method

   type(shift_solvation) :: solv
   type(shift_input), allocatable :: scratch_input
   type(wavefunction_type) :: wfn
   type(potential_type) :: pot
   type(container_cache) :: cache
   real(wp) :: energy(mol%nat)

   energy = 0.0_wp

   scratch_input = input

   if (allocated(input%solvent) .and. present(method)) then 
      call get_shift_param(scratch_input, method, error)
      if(allocated(error)) then
         call test_failed(error, "No shift parameters found for the method/solvent")
         return
      end if
   end if

   solv = shift_solvation(scratch_input)

   call solv%update(mol, cache)
   call solv%get_potential(mol, cache, wfn, pot)
   call solv%get_energy(mol, cache, wfn, energy)
   call solv%get_engrad(mol, cache, energy)

   if (abs(sum(energy) - ref) > thr) then
      call test_failed(error, "Energy does not match reference")
      print *, sum(energy),'reference:',ref
   end if
end subroutine test_e

subroutine test_e_alpb_gfn1_all_solvents(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(solvent_data) :: solvent
   type(shift_input) :: input
   integer, parameter :: nsolvents = 25
   character(len=*), parameter :: solvents(*) = [character(len=nsolvents):: &
      & "acetone", "acetonitrile", "aniline", "benzaldehyde", "benzene", &
      & "ch2cl2", "chcl3", "cs2", "dioxane", "dmf", "dmso", "ethanol", &
      & "ether", "ethylacetate", "furane", "hexadecane", "hexane", &
      & "nitromethane", "methanol", "octanol", "phenol", "thf", "toluene", &
      & "water", "woctanol"] 
   real(wp), parameter :: refs_bar1mol(*) = [&
      & 1.63261835085536E-2_wp, 9.31280584349187E-3_wp, 6.10306008943574E-3_wp, &
      & 9.67917852709751E-3_wp, 8.08739100629020E-3_wp, 6.37136322581741E-3_wp, &
      & 5.92108557794216E-3_wp, 8.53454208782486E-3_wp, 5.66153228555004E-3_wp, &
      & 1.09663650850496E-2_wp, 1.06555118341220E-2_wp, 7.00077557148834E-3_wp, &
      & 1.66506474225119E-2_wp, 7.86585894064636E-3_wp, 5.92595471558278E-3_wp, &
      & 7.25557204326916E-3_wp, 7.19072153234925E-3_wp, 8.58706969316386E-3_wp, &
      & 7.40662420091461E-3_wp, 6.21424037098337E-3_wp, 6.70655886619314E-3_wp, &
      & 9.26552024665788E-3_wp, 6.58621292223285E-3_wp, 7.27683949906328E-3_wp, &
      & 6.41785372708611E-3_wp]
   real(wp), parameter :: refs_reference(*) = [&
      & 1.87907057789662E-2_wp, 1.21001962664228E-2_wp, 8.36302467503810E-3_wp, &
      & 1.18341556383070E-2_wp, 1.03599692855620E-2_wp, 8.96891666720903E-3_wp, &
      & 8.30441563658662E-3_wp, 1.11842017076946E-2_wp, 7.98302610150582E-3_wp, &
      & 1.33878533380648E-2_wp, 1.31525893890045E-2_wp, 9.68319014349151E-3_wp, &
      & 1.87880893670823E-2_wp, 1.02312181766696E-2_wp, 8.40073773949523E-3_wp, &
      & 9.17774520460835E-3_wp, 9.11289469368843E-3_wp, 1.13509387695088E-2_wp, &
      & 1.04351691764993E-2_wp, 7.96298743088067E-3_wp, 9.00182470817862E-3_wp, &
      & 1.16308794826811E-2_wp, 8.70282029304879E-3_wp, 1.10699922271961E-2_wp, &
      & 8.16660078698342E-3_wp]
   integer :: i   

   call get_structure(mol, "MB16-43", "01")

   do i = 1, nsolvents
      ! bar1mol state
      solvent = get_solvent_data(solvents(i))
      input = shift_input(state=solution_state%bar1mol, solvent=solvent%solvent, &
         & alpb=.true.)
      call test_e(error, mol, input, refs_bar1mol(i),  method='gfn1') 
      if(allocated(error)) return  
      ! reference state
      input%state = solution_state%reference
      call test_e(error, mol, input, refs_reference(i), method='gfn1') 
      if(allocated(error)) return  
   end do 

end subroutine test_e_alpb_gfn1_all_solvents

subroutine test_e_alpb_gfn2_all_solvents(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(solvent_data) :: solvent
   type(shift_input) :: input
   integer, parameter :: nsolvents = 25
   character(len=*), parameter :: solvents(*) = [character(len=nsolvents):: &
      & "acetone", "acetonitrile", "aniline", "benzaldehyde", "benzene", &
      & "ch2cl2", "chcl3", "cs2", "dioxane", "dmf", "dmso", "ethanol", &
      & "ether", "ethylacetate", "furane", "hexadecane", "hexane", &
      & "nitromethane", "methanol", "octanol", "phenol", "thf", "toluene", &
      & "water", "woctanol"] 
   real(wp), parameter :: refs_bar1mol(*) = [&
      & 8.95537409262417E-3_wp, 8.22988757996136E-3_wp, 6.00365174170799E-3_wp, &
      & 5.29205499543163E-3_wp, 7.70568145464064E-3_wp, 7.17993018930200E-3_wp, &
      & 6.92225207178554E-3_wp, 8.01170085026346E-3_wp, 5.80234080502634E-3_wp, &
      & 9.03982927966215E-3_wp, 1.02259268778103E-2_wp, 6.44249646259719E-3_wp, &
      & 7.98467040578235E-3_wp, 5.45002953461799E-3_wp, 5.01783774913864E-3_wp, &
      & 7.42815990763873E-3_wp, 7.84458864824169E-3_wp, 5.44750394725157E-3_wp, &
      & 6.97615264444313E-3_wp, 5.37090220339446E-3_wp, 5.45892065137498E-3_wp, &
      & 8.35604827278334E-3_wp, 5.23937596608403E-3_wp, 4.11199252528792E-3_wp, &
      & 5.54333830909911E-3_wp]
   real(wp), parameter :: refs_reference(*) = [&
      & 1.14198963630367E-2_wp, 1.10172780028923E-2_wp, 8.26361632731035E-3_wp, &
      & 7.44703210664113E-3_wp, 9.97825973391251E-3_wp, 9.77748363069363E-3_wp, &
      & 9.30558213043000E-3_wp, 1.06613604701332E-2_wp, 8.12383462098212E-3_wp, &
      & 1.14613175326774E-2_wp, 1.27230044326928E-2_wp, 9.12491103460036E-3_wp, &
      & 1.01221123503528E-2_wp, 7.81538877064125E-3_wp, 7.49262077305109E-3_wp, &
      & 9.35033306897792E-3_wp, 9.76676180958088E-3_wp, 8.21137302359659E-3_wp, &
      & 1.00046976200278E-2_wp, 7.11964926329176E-3_wp, 7.75418649336047E-3_wp, &
      & 1.07214075088065E-2_wp, 7.51195424535590E-3_wp, 7.90514525342083E-3_wp, &
      & 7.29208536899641E-3_wp]
   integer :: i   

   call get_structure(mol, "MB16-43", "01")

   do i = 1, nsolvents
      ! bar1mol state
      solvent = get_solvent_data(solvents(i))
      input = shift_input(state=solution_state%bar1mol, solvent=solvent%solvent, &
         & alpb=.true.)
      call test_e(error, mol, input, refs_bar1mol(i), method='gfn2') 
      if(allocated(error)) return
      ! reference state
      input%state = solution_state%reference
      call test_e(error, mol, input, refs_reference(i), method='gfn2') 
      if(allocated(error)) return
   end do 

end subroutine test_e_alpb_gfn2_all_solvents


subroutine test_e_gbsa_gfn1_all_solvents(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(solvent_data) :: solvent
   type(shift_input) :: input
   integer, parameter :: nsolvents = 12
   character(len=*), parameter :: solvents(*) = [character(len=nsolvents):: &
      & "acetone", "acetonitrile", "benzene", "ch2cl2", "chcl3", "cs2", &
      & "dmso", "ether", "methanol", "thf", "toluene", "water"]
   real(wp), parameter :: refs_bar1mol(*) = [&
      & 5.21213510246909E-3_wp, 7.25142009001955E-3_wp, 7.26000582493303E-3_wp, &
      & 6.57710945808429E-3_wp, 6.57710945808429E-3_wp, 8.41667542117321E-3_wp, &
      & 7.25142009001955E-3_wp, 6.49277863570250E-3_wp, 5.40945961328245E-3_wp, &
      & 5.87994173698486E-3_wp, 7.26000582493303E-3_wp, 5.10670346719552E-3_wp]
   real(wp), parameter :: refs_reference(*) = [&
      & 7.67665737288168E-3_wp, 1.00388105129505E-2_wp, 9.53258410420490E-3_wp, &
      & 9.17466289947592E-3_wp, 8.96043951672875E-3_wp, 1.10709445028420E-2_wp, &
      & 9.74849764490206E-3_wp, 8.63022058027298E-3_wp, 8.43800458886716E-3_wp, &
      & 8.24530097300811E-3_wp, 9.37661319574898E-3_wp, 8.89796593470179E-3_wp]
   integer :: i   

   call get_structure(mol, "MB16-43", "01")

   do i = 1, nsolvents
      ! bar1mol state
      solvent = get_solvent_data(solvents(i))
      input = shift_input(state=solution_state%bar1mol, solvent=solvent%solvent, &
         & alpb=.false.)
      call test_e(error, mol, input, refs_bar1mol(i), method='gfn1') 
      if(allocated(error)) return
      ! reference state
      input%state = solution_state%reference
      call test_e(error, mol, input, refs_reference(i), method='gfn1') 
      if(allocated(error)) return  
   end do 

end subroutine test_e_gbsa_gfn1_all_solvents

subroutine test_e_gbsa_gfn2_all_solvents(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(solvent_data) :: solvent
   type(shift_input) :: input
   integer, parameter :: nsolvents = 14
   character(len=*), parameter :: solvents(*) = [character(len=nsolvents):: &
      & "acetone", "acetonitrile", "benzene", "ch2cl2", "chcl3", "cs2", &
      & "dmf", "dmso", "ether", "hexane", "methanol", "thf", "toluene", &
      & "water"]
   real(wp), parameter :: refs_bar1mol(*) = [&
      & 7.04678383302601E-3_wp, 5.07850817974356E-3_wp, 7.55511270005529E-3_wp, &
      & 7.47206200177234E-3_wp, 7.13802760282821E-3_wp, 7.51040424359422E-3_wp, &
      & 6.38830191828239E-3_wp, 8.19297097398559E-3_wp, 6.76591459154802E-3_wp, &
      & 7.72252129600845E-3_wp, 6.38816235066849E-3_wp, 6.43059788527107E-3_wp, &
      & 7.55511270005529E-3_wp, 4.88867595421726E-3_wp]
   real(wp), parameter :: refs_reference(*) = [&
      & 9.51130610343859E-3_wp, 7.86589860267454E-3_wp, 9.82769097932716E-3_wp, &
      & 1.00696154431639E-2_wp, 9.52135766147267E-3_wp, 1.01600638634640E-2_wp, &
      & 8.80979017129767E-3_wp, 1.06900485288681E-2_wp, 8.90335653611850E-3_wp, &
      & 9.64469445734764E-3_wp, 9.41670732625320E-3_wp, 8.79595712129432E-3_wp, &
      & 9.67172007087124E-3_wp, 8.68182868235016E-3_wp]
   integer :: i   

   call get_structure(mol, "MB16-43", "01")

   do i = 1, nsolvents
      ! bar1mol state
      solvent = get_solvent_data(solvents(i))
      input = shift_input(state=solution_state%bar1mol, solvent=solvent%solvent, &
         & alpb=.false.)
      call test_e(error, mol, input, refs_bar1mol(i), method='gfn2') 
      if(allocated(error)) return 
      ! reference state
      input%state = solution_state%reference
      call test_e(error, mol, input, refs_reference(i), method='gfn2') 
      if(allocated(error)) return 
   end do 

end subroutine test_e_gbsa_gfn2_all_solvents


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
   type(shift_input) :: input

   call get_structure(mol, "MB16-43", "04")

   ! Check GFN1/GBSA unavailable solvent
   solvent = get_solvent_data("aniline")
   input = shift_input(state=solution_state%bar1mol, solvent=solvent%solvent, &
      & alpb=.false.)
   call get_shift_param(input, 'gfn1', error)

end subroutine test_unsupported_solvent


end module test_solvation_shift
