! This file is part of s-dftd3.
! SPDX-Identifier: LGPL-3.0-or-later
!
! s-dftd3 is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! s-dftd3 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with s-dftd3.  If not, see <https://www.gnu.org/licenses/>.

module test_periodic
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type, new
   use mstore, only : get_structure
   use dftd3
   implicit none
   private

   public :: collect_periodic

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))
   real(wp), parameter :: thr3 = 1000*sqrt(epsilon(1.0_wp))
   type(realspace_cutoff), parameter :: cutoff = &
      & realspace_cutoff(cn=30_wp, disp2=60.0_wp, disp3=15.0_wp)


contains


!> Collect all exported unit tests
subroutine collect_periodic(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("PBE-D3(BJ)", test_pbed3bj_acetic), &
      & new_unittest("PBEsol-D3(BJ)", test_pbesold3bj_adaman), &
      & new_unittest("TPSS-D3(BJ)", test_tpssd3bj_ammonia), &
      & new_unittest("HSE06-D3(BJ)", test_hse06d3bj_anthracene), &
      & new_unittest("BLYP-D3(0)", test_blypd3zero_benzene), &
      & new_unittest("M06L-D3(0)", test_m06ld3zero_cyanamide), &
      & new_unittest("rPW86PBE-D3(0)", test_rpw86pbed3zero_co2), &
      & new_unittest("revSSB-D3(0)", test_revssbd3zero_cytosine), &
      & new_unittest("HSEsol-D3(BJ)-ATM", test_hsesold3bjatm_oxacb), &
      ! new_unittest("PWGGA-D3(BJ)-ATM", test_pwggad3bjatm_pyrazine), &
      & new_unittest("B3PW91-D3(0)-ATM", test_b3pw91d3zeroatm_urea), &
      & new_unittest("RPBE-D3(0)-ATM", test_rpbed3zeroatm_hexamine), &
      & new_unittest("revPBE-D3(BJ)", test_revpbed3bj_1d), &
      & new_unittest("BP86-D3(0)", test_bp86d3zero_2d) &
      & ]

end subroutine collect_periodic


subroutine test_dftd3_gen(error, mol, param, ref)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Damping parameters
   class(damping_param), intent(in) :: param

   !> Expected dispersion energy
   real(wp), intent(in) :: ref

   type(d3_model) :: d3
   real(wp) :: energy

   call new_d3_model(d3, mol)
   call get_dispersion(mol, d3, param, cutoff, energy)

   call check(error, energy, ref, thr=thr)
   if (allocated(error)) then
      print*,energy
   end if

end subroutine test_dftd3_gen


subroutine test_numgrad(error, mol, param, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Damping parameters
   class(damping_param), intent(in) :: param

   !> Gradient threshold
   real(wp), intent(in), optional :: thr_in

   integer :: iat, ic
   type(d3_model) :: d3
   real(wp) :: energy, er, el, sigma(3, 3), thr_grad, numerr
   real(wp), allocatable :: gradient(:, :), numgrad(:, :)
   real(wp), parameter :: step = 1.0e-6_wp
   logical :: tight

   allocate(gradient(3, mol%nat), numgrad(3, mol%nat))
   call new_d3_model(d3, mol)
   thr_grad = thr2
   if (present(thr_in)) thr_grad = thr_in
   tight = present(thr_in)

   do iat = 1, mol%nat
      do ic = 1, 3
         if (tight) then
            call ridders_cartesian(ic, iat, numgrad(ic, iat), numerr)
         else
            mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
            call get_dispersion(mol, d3, param, cutoff, er)
            mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
            call get_dispersion(mol, d3, param, cutoff, el)
            mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
            numgrad(ic, iat) = 0.5_wp*(er - el)/step
         end if
      end do
   end do

   call get_dispersion(mol, d3, param, cutoff, energy, gradient, sigma)

   if (any(abs(gradient - numgrad) > thr_grad)) then
      call test_failed(error, "Gradient of dispersion energy does not match")
      print'(3es21.14)', gradient-numgrad
   end if

contains

subroutine ridders_cartesian(ic, iat, deriv, err)

   integer, intent(in) :: ic, iat
   real(wp), intent(out) :: deriv, err

   integer, parameter :: max_tab = 10
   real(wp), parameter :: initial_step = 1.0e-4_wp
   real(wp), parameter :: step_div = sqrt(2.0_wp)
   real(wp), parameter :: step_div_2 = step_div**2
   integer :: iter, order
   real(wp) :: table(max_tab, max_tab), factor, err_est, prev_deriv, prev_err
   real(wp) :: h, x0

   x0 = mol%xyz(ic, iat)
   h = initial_step
   err = huge(1.0_wp)
   prev_err = err

   mol%xyz(ic, iat) = x0 + h
   call get_dispersion(mol, d3, param, cutoff, er)
   mol%xyz(ic, iat) = x0 - h
   call get_dispersion(mol, d3, param, cutoff, el)
   table(1, 1) = 0.5_wp*(er - el)/h
   deriv = table(1, 1)
   prev_deriv = deriv

   do iter = 2, max_tab
      h = h / step_div
      mol%xyz(ic, iat) = x0 + h
      call get_dispersion(mol, d3, param, cutoff, er)
      mol%xyz(ic, iat) = x0 - h
      call get_dispersion(mol, d3, param, cutoff, el)
      table(iter, 1) = 0.5_wp*(er - el)/h

      factor = step_div_2
      do order = 1, iter - 1
         factor = factor * step_div_2
         table(iter, order + 1) = (factor*table(iter, order) &
            & - table(iter - 1, order))/(factor - 1.0_wp)
         err_est = max(abs(table(iter, order + 1) - table(iter, order)), &
            & abs(table(iter, order + 1) - table(iter - 1, order)))
         if (err_est <= err) then
            err = err_est
            deriv = table(iter, order + 1)
         end if
      end do

      if (abs(table(iter, iter) - table(iter - 1, iter - 1)) >= 2*err &
         & .and. iter > 2) then
         deriv = prev_deriv
         err = prev_err
         exit
      end if

      prev_deriv = deriv
      prev_err = err
   end do

   mol%xyz(ic, iat) = x0

end subroutine ridders_cartesian

end subroutine test_numgrad


subroutine test_numsigma(error, mol, param, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Damping parameters
   class(damping_param), intent(in) :: param

   !> Strain derivative threshold
   real(wp), intent(in), optional :: thr_in

   integer :: ic, jc
   type(d3_model) :: d3
   real(wp) :: energy, er, el, sigma(3, 3), eps(3, 3), numsigma(3, 3), lattice(3, 3)
   real(wp) :: thr_sigma, numerr
   real(wp), allocatable :: gradient(:, :), xyz(:, :)
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))
   real(wp), parameter :: step = 1.0e-7_wp
   logical :: tight

   allocate(gradient(3, mol%nat), xyz(3, mol%nat))
   call new_d3_model(d3, mol)
   thr_sigma = thr3
   if (present(thr_in)) thr_sigma = thr_in
   tight = present(thr_in)

   eps(:, :) = unity
   xyz(:, :) = mol%xyz
   lattice(:, :) = mol%lattice
   do ic = 1, 3
      do jc = 1, 3
         if (tight) then
            call ridders_strain(ic, jc, numsigma(jc, ic), numerr)
         else
            eps(jc, ic) = eps(jc, ic) + step
            mol%xyz(:, :) = matmul(eps, xyz)
            mol%lattice(:, :) = matmul(eps, lattice)
            call get_dispersion(mol, d3, param, cutoff, er)
            eps(jc, ic) = eps(jc, ic) - 2*step
            mol%xyz(:, :) = matmul(eps, xyz)
            mol%lattice(:, :) = matmul(eps, lattice)
            call get_dispersion(mol, d3, param, cutoff, el)
            eps(jc, ic) = eps(jc, ic) + step
            mol%xyz(:, :) = xyz
            mol%lattice(:, :) = lattice
            numsigma(jc, ic) = 0.5_wp*(er - el)/step
         end if
      end do
   end do

   call get_dispersion(mol, d3, param, cutoff, energy, gradient, sigma)

   if (any(abs(sigma - numsigma) > thr_sigma)) then
      call test_failed(error, "Strain derivatives do not match")
      print'(3es21.14)', sigma-numsigma
   end if

contains

subroutine ridders_strain(ic, jc, deriv, err)

   integer, intent(in) :: ic, jc
   real(wp), intent(out) :: deriv, err

   integer, parameter :: max_tab = 10
   real(wp), parameter :: initial_step = 5.0e-8_wp
   real(wp), parameter :: step_div = sqrt(2.0_wp)
   real(wp), parameter :: step_div_2 = step_div**2
   integer :: iter, order
   real(wp) :: table(max_tab, max_tab), factor, err_est, prev_deriv, prev_err
   real(wp) :: h

   h = initial_step
   err = huge(1.0_wp)
   prev_err = err

   call strain_energy(ic, jc, h, er)
   call strain_energy(ic, jc, -h, el)
   table(1, 1) = 0.5_wp*(er - el)/h
   deriv = table(1, 1)
   prev_deriv = deriv

   do iter = 2, max_tab
      h = h / step_div
      call strain_energy(ic, jc, h, er)
      call strain_energy(ic, jc, -h, el)
      table(iter, 1) = 0.5_wp*(er - el)/h

      factor = step_div_2
      do order = 1, iter - 1
         factor = factor * step_div_2
         table(iter, order + 1) = (factor*table(iter, order) &
            & - table(iter - 1, order))/(factor - 1.0_wp)
         err_est = max(abs(table(iter, order + 1) - table(iter, order)), &
            & abs(table(iter, order + 1) - table(iter - 1, order)))
         if (err_est <= err) then
            err = err_est
            deriv = table(iter, order + 1)
         end if
      end do

      if (abs(table(iter, iter) - table(iter - 1, iter - 1)) >= 2*err &
         & .and. iter > 2) then
         deriv = prev_deriv
         err = prev_err
         exit
      end if

      prev_deriv = deriv
      prev_err = err
   end do

   mol%xyz(:, :) = xyz
   mol%lattice(:, :) = lattice

end subroutine ridders_strain

subroutine strain_energy(ic, jc, step, energy)

   integer, intent(in) :: ic, jc
   real(wp), intent(in) :: step
   real(wp), intent(out) :: energy
   real(wp) :: strain(3, 3)

   strain(:, :) = unity
   strain(jc, ic) = strain(jc, ic) + step
   mol%xyz(:, :) = matmul(strain, xyz)
   mol%lattice(:, :) = matmul(strain, lattice)
   call get_dispersion(mol, d3, param, cutoff, energy)

end subroutine strain_energy

end subroutine test_numsigma


subroutine test_pbed3bj_acetic(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, &
      & a1 = 0.4289_wp, s8 = 0.7875_wp, a2 = 4.4407_wp)

   call get_structure(mol, "X23", "acetic")
   call new_rational_damping(param, inp)
   call test_dftd3_gen(error, mol, param, -6.6732836815486210E-2_wp)

end subroutine test_pbed3bj_acetic


subroutine test_pbesold3bj_adaman(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, &
      & a1 = 0.4466_wp, s8 = 2.9491_wp, a2 = 6.1742_wp)

   call get_structure(mol, "X23", "adaman")
   call new_rational_damping(param, inp)
   call test_dftd3_gen(error, mol, param, -7.9313514714713124E-2_wp)

end subroutine test_pbesold3bj_adaman


subroutine test_tpssd3bj_ammonia(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, &
      & a1 = 0.4535_wp, s8 = 1.9435_wp, a2 = 4.4752_wp)

   call get_structure(mol, "X23", "ammonia")
   call new_rational_damping(param, inp)
   call test_numgrad(error, mol, param)

end subroutine test_tpssd3bj_ammonia


subroutine test_hse06d3bj_anthracene(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, &
      & a1 = 0.383_wp, s8 = 2.310_wp, a2 = 5.685_wp)

   call get_structure(mol, "X23", "anthracene")
   call new_rational_damping(param, inp)
   call test_numsigma(error, mol, param)

end subroutine test_hse06d3bj_anthracene


subroutine test_blypd3zero_benzene(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(zero_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, rs8 = 1.0_wp, &
      & rs6 = 1.094_wp, s8 = 1.682_wp)

   call get_structure(mol, "X23", "benzene")
   call new_zero_damping(param, inp)
   call test_dftd3_gen(error, mol, param, -0.16647358662571451_wp)

end subroutine test_blypd3zero_benzene


subroutine test_m06ld3zero_cyanamide(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(zero_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, rs8 = 1.0_wp, &
      & rs6 = 1.581_wp, s8 = 0.000_wp)

   call get_structure(mol, "X23", "cyanamide")
   call new_zero_damping(param, inp)
   call test_dftd3_gen(error, mol, param, -2.3222556050310025E-2_wp)

end subroutine test_m06ld3zero_cyanamide


subroutine test_rpw86pbed3zero_co2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(zero_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, rs8 = 1.0_wp, &
      & rs6 = 1.224_wp, s8 = 0.901_wp)

   call get_structure(mol, "X23", "CO2")
   call new_zero_damping(param, inp)
   call test_numgrad(error, mol, param)

end subroutine test_rpw86pbed3zero_co2


subroutine test_revssbd3zero_cytosine(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(zero_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, rs8 = 1.0_wp, &
      & rs6 = 1.221_wp, s8 = 0.560_wp)

   call get_structure(mol, "X23", "cytosine")
   call new_zero_damping(param, inp)
   call test_numsigma(error, mol, param)

end subroutine test_revssbd3zero_cytosine


subroutine test_hsesold3bjatm_oxacb(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 14.0_wp, &
      & a1 = 0.4650_wp, s8 = 2.9215_wp, a2 = 6.2003_wp)

   call get_structure(mol, "X23", "oxacb")
   call new_rational_damping(param, inp)
   call test_numgrad(error, mol, param, 1.0e-11_wp)

end subroutine test_hsesold3bjatm_oxacb


subroutine test_pwggad3bjatm_pyrazine(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 14.0_wp, &
      & a1 = 0.2211_wp, s8 = 2.6910_wp, a2 = 6.7278_wp)

   call get_structure(mol, "X23", "pyrazine")
   call new_rational_damping(param, inp)
   call test_numsigma(error, mol, param, 1.0e-8_wp)

end subroutine test_pwggad3bjatm_pyrazine


subroutine test_b3pw91d3zeroatm_urea(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(zero_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 14.0_wp, rs8 = 1.0_wp, &
      & rs6 = 1.176_wp, s8 = 1.775_wp)

   call get_structure(mol, "X23", "urea")
   call new_zero_damping(param, inp)
   call test_numgrad(error, mol, param, 1.0e-11_wp)

end subroutine test_b3pw91d3zeroatm_urea


subroutine test_rpbed3zeroatm_hexamine(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(zero_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 14.0_wp, rs8 = 1.0_wp, &
      & rs6 = 0.872_wp, s8 = 0.514_wp)

   call get_structure(mol, "X23", "hexamine")
   call new_zero_damping(param, inp)
   call test_numsigma(error, mol, param, 1.0e-8_wp)

end subroutine test_rpbed3zeroatm_hexamine


subroutine test_revpbed3bj_1d(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nat = 32
   integer, parameter :: num(nat) = 6
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & 1.36794785746435_wp, 13.45808943446053_wp,  8.83754983226359_wp, &
      & 3.69183290816438_wp, 13.13552229161569_wp, 10.16652201690950_wp, &
      & 1.36792668081267_wp, 10.38660504434782_wp, 13.04411926632965_wp, &
      & 3.69180534781206_wp, 11.55414582295511_wp, 12.33193380846742_wp, &
      & 1.36791549262702_wp,  3.53066844289674_wp, 10.38660588677206_wp, &
      & 1.36792046664920_wp,  7.73723910626293_wp, 13.45809224934817_wp, &
      & 3.69181279359489_wp,  6.40826717723392_wp, 13.13552570942280_wp, &
      & 1.36792009865062_wp,  3.11669712338516_wp,  7.73723850632628_wp, &
      & 3.69181515738094_wp,  3.43926499914873_wp,  6.40826580885474_wp, &
      & 3.69178443989294_wp,  4.24285720771059_wp, 11.55415026712869_wp, &
      & 1.36790824853106_wp,  6.18818490375705_wp,  3.53066863732142_wp, &
      & 3.69178194163078_wp,  5.02063901427657_wp,  4.24285736953327_wp, &
      & 1.36794124909207_wp, 13.04411858182861_wp,  6.18818324080182_wp, &
      & 1.36792249732236_wp,  8.83755133592807_wp,  3.11669686076913_wp, &
      & 3.69182456413952_wp, 10.16652118921143_wp,  3.43926084011816_wp, &
      & 3.69181444966104_wp, 12.33193631088573_wp,  5.02063847821044_wp, &
      & 6.01572566324028_wp, 13.45790756713123_wp,  8.83752222635545_wp, &
      & 8.33965926123256_wp, 13.13576644753615_wp, 10.16660228658307_wp, &
      & 6.01574747573805_wp, 10.38654070512969_wp, 13.04391961251944_wp, &
      & 8.33964066450677_wp, 11.55427002850905_wp, 12.33211653730939_wp, &
      & 6.01574728097580_wp,  3.53087013230607_wp, 10.38654217813321_wp, &
      & 6.01568913853645_wp,  7.73726406411719_wp, 13.45790864082374_wp, &
      & 8.33963586549168_wp,  6.40818371470975_wp, 13.13576911116618_wp, &
      & 6.01568179676984_wp,  3.11688332536281_wp,  7.73726611148835_wp, &
      & 8.33963704688671_wp,  3.43902559351770_wp,  6.40818390180453_wp, &
      & 8.33962496288127_wp,  4.24267007149867_wp, 11.55427031066552_wp, &
      & 6.01573464280675_wp,  6.18824653544318_wp,  3.53086861480278_wp, &
      & 8.33961857277245_wp,  5.02052001792996_wp,  4.24267413625204_wp, &
      & 6.01575677304189_wp, 13.04392044501564_wp,  6.18824448603611_wp, &
      & 6.01568344836224_wp,  8.83752193432504_wp,  3.11688171781516_wp, &
      & 8.33964228963694_wp, 10.16660428027860_wp,  3.43902155668011_wp, &
      & 8.33965118613331_wp, 12.33211762632282_wp,  5.02051902430387_wp],&
      & shape(xyz))
   logical, parameter :: periodic(3) = [.true., .false., .false.]
   real(wp), parameter :: lattice(3, 3) = reshape([&
      & 9.2955628527586_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], &
      & shape(lattice))
   real(wp), parameter :: lattice_vac(3, 3) = reshape([&
      & 9.2955628527586_wp, 0.0_wp, 0.0_wp, 0.0_wp, 100.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 100.0_wp], &
      & shape(lattice))
   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, &
      & a1 = 0.5238_wp, s8 = 2.3550_wp, a2 = 3.5016_wp)

   call new(mol, num, xyz, periodic=periodic, lattice=lattice)
   call new_rational_damping(param, inp)
   call test_dftd3_gen(error, mol, param, -0.25436101196186600_wp)

   call new(mol, num, xyz, periodic=[.true.], lattice=lattice_vac)
   call test_dftd3_gen(error, mol, param, -0.25436101196186600_wp)

end subroutine test_revpbed3bj_1d


subroutine test_bp86d3zero_2d(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nat = 4
   integer, parameter :: num(nat) = 6
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & -0.12918412100093_wp,  0.06210659750976_wp, -2.13384498734326_wp, &
      &  0.12856915667443_wp, -0.07403227791901_wp,  4.02358027265954_wp, &
      & -0.12317720857511_wp,  2.75170732207802_wp, -2.13345350602279_wp, &
      &  2.44816466162280_wp,  1.28612566399214_wp,  4.02317048854901_wp],&
      & shape(xyz))
   logical, parameter :: periodic(3) = [.true., .true., .false.]
   real(wp), parameter :: lattice(3, 3) = reshape([&
      &  4.68837849314507_wp, 0.00000000000000_wp, 0.00000000000000_wp, &
      & -2.36282788044783_wp, 4.04978545156612_wp, 0.00000000000000_wp, &
      &  0.00000000000000_wp, 0.00000000000000_wp, 1.00000000000000_wp],&
      & shape(lattice))
   real(wp), parameter :: lattice_vac(3, 3) = reshape([&
      &  4.68837849314507_wp, 0.00000000000000_wp, 0.00000000000000_wp, &
      & -2.36282788044783_wp, 4.04978545156612_wp, 0.00000000000000_wp, &
      &  0.00000000000000_wp, 0.00000000000000_wp, 100.000000000000_wp],&
      & shape(lattice))
   type(structure_type) :: mol
   type(zero_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, &
      & rs6 = 1.139_wp, s8 = 1.683_wp)

   call new(mol, num, xyz, periodic=periodic, lattice=lattice)
   call new_zero_damping(param, inp)
   call test_dftd3_gen(error, mol, param, -1.5983949078995030E-2_wp)

   call new(mol, num, xyz, periodic=[.true.], lattice=lattice_vac)
   call test_dftd3_gen(error, mol, param, -1.5983949078995030E-2_wp)

end subroutine test_bp86d3zero_2d


end module test_periodic
