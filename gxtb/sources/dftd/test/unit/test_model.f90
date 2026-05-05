! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the Lesser GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! Lesser GNU General Public License for more details.
!
! You should have received a copy of the Lesser GNU General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

module test_model
   use dftd4_cache, only : dispersion_cache
   use dftd4_cutoff, only : get_lattice_points
   use dftd4_data, only : get_covalent_rad
   use dftd4_model, only : dispersion_model, new_dispersion_model, d4_qmod, &
      & get_dispersion_model_id, dftd_models
   use dftd4_model_d4, only : d4_model, new_d4_model
   use dftd4_model_d4s, only : d4s_model, new_d4s_model
   use dftd4_model_d4srev, only : d4srev_model, new_d4srev_model
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, &
      & check, test_failed
   use mctc_io_structure, only : new, structure_type
   use mctc_ncoord, only : new_ncoord, ncoord_type, cn_count
   use mstore, only : get_structure
   use multicharge, only : get_charges
   implicit none
   private

   public :: collect_model

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))


contains


!> Collect all exported unit tests
subroutine collect_model(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("c6-D4-mb01", test_c6_d4_mb01), &
      & new_unittest("c6-D4S-mb01", test_c6_d4s_mb01), &
      & new_unittest("c6-D4-mb02", test_c6_d4_mb02), &
      & new_unittest("c6-D4-EEQBC-mb02", test_c6_d4_eeqbc_mb02), &
      & new_unittest("c6-D4-mb03", test_c6_d4_mb03), &
      & new_unittest("dc6-D4-mb04", test_dc6_d4_mb04), &
      & new_unittest("dc6-D4S-mb04", test_dc6_d4s_mb04), &
      & new_unittest("dc6-D4-mb05", test_dc6_d4_mb05), &
      & new_unittest("dc6-D4S-mb05", test_dc6_d4s_mb05), &
      & new_unittest("dc6-D4S-EEQBC-mb05", test_dc6_d4s_eeqbc_mb05), &
      & new_unittest("dc6-D4-mb06", test_dc6_d4_mb06), &
      & new_unittest("dc6-D4S-mb06", test_dc6_d4s_mb06), &
      & new_unittest("c6-D4-gfn2", test_c6_d4_mb07), &
      & new_unittest("dc6-D4-gfn2", test_dc6_d4_mb08), &
      & new_unittest("dc6-D4S-gfn2", test_dc6_d4s_mb08), &
      & new_unittest("pol-D4-mb09", test_pol_d4_mb09), &
      & new_unittest("pol-D4s-mb09", test_pol_d4s_mb09), &
      & new_unittest("dpol-D4-mb10", test_dpol_d4_mb10), &
      & new_unittest("dpol-D4S-mb10", test_dpol_d4s_mb10), &
      & new_unittest("model-D4-error", test_d4_model_error, should_fail=.true.), &
      & new_unittest("model-D4S-error", test_d4s_model_error, should_fail=.true.), &
      & new_unittest("model-wrapper", test_model_wrapper), &
      & new_unittest("model-wrapper-fail", test_model_wrapper_fail, should_fail=.true.), &
      & new_unittest("id-from-name", test_id_from_name), &
      & new_unittest("id-from-name-empty", test_id_from_name_empty, should_fail=.true.), &
      & new_unittest("id-from-name-wrong", test_id_from_name_wrong, should_fail=.true.), &
      & new_unittest("id-from-model", test_id_from_model) &
      & ]

end subroutine collect_model


subroutine test_c6_gen(error, mol, d4, ref_c6, with_cn, with_q, qat)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type) :: mol

   !> Dispersion model
   class(dispersion_model), intent(in) :: d4

   !> Reference C6 coefficients
   real(wp), intent(in) :: ref_c6(:, :)

   !> Calculate coordination number
   logical, intent(in) :: with_cn

   !> Calculate atomic charges
   logical, intent(in) :: with_q

   !> Atomic charges
   real(wp), optional, intent(in) :: qat(:)

   real(wp), allocatable :: cn(:), q(:)
   type(dispersion_cache) :: cache
   real(wp), parameter :: cutoff = 30.0_wp
   real(wp), allocatable :: lattr(:, :)
   class(ncoord_type), allocatable :: ncoord

   allocate(cn(mol%nat), q(mol%nat))
   cn(:) = 0.0_wp
   q(:) = 0.0_wp

   if (with_cn) then
      call new_ncoord(ncoord, mol, cn_count%dftd4, &
         & cutoff=cutoff, rcov=d4%rcov, en=d4%en, error=error)
      if (allocated(error)) return
      call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
      call ncoord%get_coordination_number(mol, lattr, cn)
   end if
   if (with_q) then
      if(present(qat)) then
         q(:) = qat
      else
         call get_charges(d4%mchrg, mol, error, q)
         if (allocated(error)) return
      end if
   end if

   call d4%update(mol, cache, cn, q, grad=.false.)
   if (any(abs(cache%c6 - ref_c6) > thr2)) then
      call test_failed(error, "C6 coefficients do not match")
      where(abs(cache%c6) < thr) cache%c6 = 0.0_wp
      print'(3(es20.13,"_wp,"), " &")', cache%c6
   end if

end subroutine test_c6_gen


subroutine test_dc6_gen(error, mol, d4, with_cn, with_q, qat)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type) :: mol

   !> Dispersion model
   class(dispersion_model), intent(in) :: d4

   !> Calculate coordination number
   logical, intent(in) :: with_cn

   !> Calculate atomic charges
   logical, intent(in) :: with_q

   !> Atomic charges
   real(wp), optional, intent(in) :: qat(:)

   integer :: iat, jat
   real(wp), allocatable :: cn(:), q(:)
   real(wp), allocatable :: c6(:, :), c6dcn(:, :), c6dq(:, :)
   real(wp), allocatable :: c6r(:, :), c6l(:, :), numdcn(:, :), numdq(:, :)
   type(dispersion_cache) :: cache
   real(wp), parameter :: cutoff = 30.0_wp, lattr(3, 1) = 0.0_wp
   real(wp), parameter :: step = 5.0e-6_wp
   class(ncoord_type), allocatable :: ncoord

   allocate(cn(mol%nat), q(mol%nat), c6(mol%nat, mol%nat), &
      & c6dcn(mol%nat, mol%nat), c6dq(mol%nat, mol%nat), &
      & c6r(mol%nat, mol%nat), c6l(mol%nat, mol%nat), &
      & numdcn(mol%nat, mol%nat), numdq(mol%nat, mol%nat))
   cn(:) = 0.0_wp
   q(:) = 0.0_wp

   if (with_cn) then
      call new_ncoord(ncoord, mol, cn_count%dftd4, &
         & cutoff=cutoff, rcov=d4%rcov, en=d4%en, error=error)
      if (allocated(error)) return
      call ncoord%get_coordination_number(mol, lattr, cn)
   end if
   if (with_q) then
      if(present(qat)) then
         q(:) = qat
      else
         call get_charges(d4%mchrg, mol, error, q)
         if (allocated(error)) return
      end if
   end if

   numdcn(:, :) = 0.0_wp
   if (with_cn) then
      do iat = 1, mol%nat
         cn(iat) = cn(iat) + step
         call d4%update(mol, cache, cn, q, grad=.false.)
         c6r = cache%c6
         cn(iat) = cn(iat) - 2*step
         call d4%update(mol, cache, cn, q, grad=.false.)
         c6l = cache%c6
         cn(iat) = cn(iat) + step
         do jat = 1, mol%nat
            if (iat .ne. jat) then
               numdcn(iat, jat) = numdcn(iat, jat) + 0.5_wp*(c6r(iat, jat) - c6l(iat, jat))/step
            else
               numdcn(iat, jat) = numdcn(iat, jat) + 0.25_wp*(c6r(iat, jat) - c6l(iat, jat))/step
            end if
         end do
      end do
   end if

   numdq(:, :) = 0.0_wp
   if (with_q) then
      do iat = 1, mol%nat
         q(iat) = q(iat) + step
         call d4%update(mol, cache, cn, q, grad=.false.)
         c6r = cache%c6
         q(iat) = q(iat) - 2*step
         call d4%update(mol, cache, cn, q, grad=.false.)
         c6l = cache%c6
         q(iat) = q(iat) + step
         do jat = 1, mol%nat
            if (iat .ne. jat) then
               numdq(iat, jat) = numdq(iat, jat) + 0.5_wp*(c6r(iat, jat) - c6l(iat, jat))/step
            else
               numdq(iat, jat) = numdq(iat, jat) + 0.25_wp*(c6r(iat, jat) - c6l(iat, jat))/step
            end if
         end do
      end do
   end if

   call d4%update(mol, cache, cn, q, grad=.true.)
   c6 = cache%c6
   c6dcn = cache%dc6dcn
   c6dq = cache%dc6dq

   if (with_cn .and. any(abs(c6dcn - numdcn) > 2*thr2)) then
     call test_failed(error, "C6 CN derivatives do not match")
     print'(3es21.14)', c6dcn
     print'("---")'
     print'(3es21.14)', numdcn
     print'("---")'
     print'(3es21.14)', c6dcn - numdcn
   end if

   if (with_q .and. any(abs(c6dq - numdq) > 2*thr2)) then
      call test_failed(error, "C6 Q derivatives do not match")
      print'(3es21.14)', c6dq
      print'("---")'
      print'(3es21.14)', numdq
      print'("---")'
      print'(3es21.14)', c6dq - numdq
   end if

end subroutine test_dc6_gen


subroutine test_pol_gen(error, mol, d4, ref, ref_qq, with_cn, &
   & with_q, qat)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type) :: mol

   !> Dispersion model
   class(dispersion_model), intent(in) :: d4

   !> Reference dipole-dipole polarizabilities
   real(wp), intent(in) :: ref(:)

   !> Reference quadrupole-quadrupole polarizabilities
   real(wp), intent(in) :: ref_qq(:)

   !> Calculate coordination number
   logical, intent(in) :: with_cn

   !> Calculate atomic charges
   logical, intent(in) :: with_q

   !> Atomic charges
   real(wp), optional, intent(in) :: qat(:)

   real(wp), allocatable :: cn(:), q(:), alpha(:), alphaqq(:)
   type(dispersion_cache) :: cache
   real(wp), parameter :: cutoff = 30.0_wp
   real(wp), allocatable :: lattr(:, :)
   class(ncoord_type), allocatable :: ncoord

   allocate(cn(mol%nat), q(mol%nat), alpha(mol%nat), &
      & alphaqq(mol%nat))
   cn(:) = 0.0_wp
   q(:) = 0.0_wp

   if (with_cn) then
      call new_ncoord(ncoord, mol, cn_count%dftd4, &
         & cutoff=cutoff, rcov=d4%rcov, en=d4%en, error=error)
      if (allocated(error)) return
      call get_lattice_points(mol%periodic, mol%lattice, &
         & cutoff, lattr)
      call ncoord%get_coordination_number(mol, lattr, cn)
   end if

   if (with_q) then
      if(present(qat)) then
         q(:) = qat
      else
         call get_charges(d4%mchrg, mol, error, q)
         if (allocated(error)) return
      end if
   end if

   call d4%update(mol, cache, cn, q, grad=.false.)

   call d4%get_polarizabilities(cache, alpha, alphaqq)

   if (any(abs(alpha - ref) > thr2)) then
      call test_failed(error, "Polarizabilities do not match")
      where(abs(alpha) < thr) alpha = 0.0_wp
      print'(3(es20.13,"_wp,"), " &")', alpha
   end if
   
   if (any(abs(alphaqq - ref_qq) > thr2)) then
      call test_failed(error, &
         & "Quadrupole polarizabilities do not match")
      where(abs(alphaqq) < thr) alphaqq = 0.0_wp
      print'(3(es20.13,"_wp,"), " &")', alphaqq
   end if

end subroutine test_pol_gen


subroutine test_dpol_gen(error, mol, d4, with_cn, with_q, qat)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type) :: mol

   !> Dispersion model
   class(dispersion_model), intent(in) :: d4

   !> Calculate coordination number
   logical, intent(in) :: with_cn

   !> Calculate atomic charges
   logical, intent(in) :: with_q

   !> Atomic charges
   real(wp), optional, intent(in) :: qat(:)

   integer :: iat
   real(wp), allocatable :: cn(:), q(:)
   real(wp), allocatable :: alpha(:), alphadcn(:), alphadq(:)
   real(wp), allocatable :: alphaqq(:), alphaqqcn(:), alphaqqq(:)
   real(wp), allocatable :: alphar(:), alphal(:), numdcn(:), numdq(:)
   real(wp), allocatable :: qqr(:), qql(:), numqqcn(:), numqqq(:)
   type(dispersion_cache) :: cache
   real(wp), parameter :: cutoff = 30.0_wp, lattr(3, 1) = 0.0_wp
   real(wp), parameter :: step = 1.0e-5_wp
   class(ncoord_type), allocatable :: ncoord

   allocate(cn(mol%nat), q(mol%nat), &
      & alpha(mol%nat), alphadcn(mol%nat), alphadq(mol%nat), &
      & alphaqq(mol%nat), alphaqqcn(mol%nat), alphaqqq(mol%nat), &
      & alphar(mol%nat), alphal(mol%nat), numdcn(mol%nat), &
      & numdq(mol%nat), qqr(mol%nat), qql(mol%nat), &
      & numqqcn(mol%nat), numqqq(mol%nat))
   cn(:) = 0.0_wp
   q(:) = 0.0_wp

   if (with_cn) then
      call new_ncoord(ncoord, mol, cn_count%dftd4, &
         & cutoff=cutoff, rcov=d4%rcov, en=d4%en, error=error)
      if (allocated(error)) return
      call ncoord%get_coordination_number(mol, lattr, cn)
   end if

   if (with_q) then
      if(present(qat)) then
         q(:) = qat
      else
         call get_charges(d4%mchrg, mol, error, q)
         if (allocated(error)) return
      end if
   end if

   if (with_cn) then
      do iat = 1, mol%nat
         cn(iat) = cn(iat) + step
         call d4%update(mol, cache, cn, q, grad=.false.)
         call d4%get_polarizabilities(cache, alphar, qqr)
         cn(iat) = cn(iat) - 2*step
         call d4%update(mol, cache, cn, q, grad=.false.)
         call d4%get_polarizabilities(cache, alphal, qql)
         cn(iat) = cn(iat) + step
         numdcn(iat) = 0.5_wp*(alphar(iat) - alphal(iat))/step
         numqqcn(iat) = 0.5_wp*(qqr(iat) - qql(iat))/step
      end do
   end if

   if (with_q) then
      do iat = 1, mol%nat
         q(iat) = q(iat) + step
         call d4%update(mol, cache, cn, q, grad=.false.)
         call d4%get_polarizabilities(cache, alphar, qqr)
         q(iat) = q(iat) - 2*step
         call d4%update(mol, cache, cn, q, grad=.false.)
         call d4%get_polarizabilities(cache, alphal, qql)
         q(iat) = q(iat) + step
         numdq(iat) = 0.5_wp*(alphar(iat) - alphal(iat))/step
         numqqq(iat) = 0.5_wp*(qqr(iat) - qql(iat))/step
      end do
   end if

   call d4%update(mol, cache, cn, q, grad=.true.)
   call d4%get_polarizabilities(cache, alpha, alphaqq, &
      & dadcn=alphadcn, dadq=alphadq, daqqdcn=alphaqqcn, &
      & daqqdq=alphaqqq)


   if (with_cn .and. any(abs(alphadcn - numdcn) > thr2)) then
     call test_failed(error, "Polarizability CN derivatives do not match")
     print'(3es21.14)', alphadcn
     print'("---")'
     print'(3es21.14)', numdcn
     print'("---")'
     print'(3es21.14)', alphadcn - numdcn
   end if

   if (with_q .and. any(abs(alphadq - numdq) > thr2)) then
      call test_failed(error, "Polarizability Q derivatives do not match")
      print'(3es21.14)', alphadq
      print'("---")'
      print'(3es21.14)', numdq
      print'("---")'
      print'(3es21.14)', alphadq - numdq
   end if

   if (with_cn .and. any(abs(alphaqqcn - numqqcn) > thr2)) then
     call test_failed(error, &
        & "Quadrupole polarizability CN derivatives do not match")
     print'(3es21.14)', alphaqqcn
     print'("---")'
     print'(3es21.14)', numqqcn
     print'("---")'
     print'(3es21.14)', alphaqqcn - numqqcn
   end if

   if (with_q .and. any(abs(alphaqqq - numqqq) > thr2)) then
      call test_failed(error, &
         & "Quadrupole polarizability Q derivatives do not match")
      print'(3es21.14)', alphaqqq
      print'("---")'
      print'(3es21.14)', numqqq
      print'("---")'
      print'(3es21.14)', alphaqqq - numqqq
   end if

end subroutine test_dpol_gen


subroutine test_c6_d4_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   real(wp), parameter :: ref(16, 16) = reshape([&
      & 1.0673751057011E+02_wp, 1.3341905517525E+01_wp, 2.5653153287381E+01_wp, &
      & 1.3286256032550E+01_wp, 1.8172837654743E+01_wp, 1.3407379807008E+01_wp, &
      & 1.3286079405123E+01_wp, 2.5572059731044E+01_wp, 3.2413191615239E+01_wp, &
      & 1.3488257086895E+01_wp, 1.3407544045253E+01_wp, 7.3738556418539E+01_wp, &
      & 4.4740018054061E+01_wp, 4.7377208046689E+01_wp, 3.2412969775949E+01_wp, &
      & 1.0847087588239E+02_wp, 1.3341905517525E+01_wp, 3.0564966218613E+00_wp, &
      & 6.5553944080817E+00_wp, 3.0470665530425E+00_wp, 5.1002826817488E+00_wp, &
      & 3.0675915503404E+00_wp, 3.0470366226884E+00_wp, 6.5291578413266E+00_wp, &
      & 7.7571265678988E+00_wp, 3.0812965885013E+00_wp, 3.0676193812890E+00_wp, &
      & 1.7231850034693E+01_wp, 1.0121118362811E+01_wp, 1.0554487753558E+01_wp, &
      & 7.7570715810400E+00_wp, 2.0891158034439E+01_wp, 2.5653153287381E+01_wp, &
      & 6.5553944080817E+00_wp, 1.4581513138206E+01_wp, 6.5365717713675E+00_wp, &
      & 1.1626599584016E+01_wp, 6.5775401434039E+00_wp, 6.5365120296846E+00_wp, &
      & 1.4518124333712E+01_wp, 1.6957726525989E+01_wp, 6.6048957183646E+00_wp, &
      & 6.5775956946247E+00_wp, 3.7235277639546E+01_wp, 2.1680017664472E+01_wp, &
      & 2.2530197721674E+01_wp, 1.6957605864802E+01_wp, 4.3281495531031E+01_wp, &
      & 1.3286256032550E+01_wp, 3.0470665530425E+00_wp, 6.5365717713675E+00_wp, &
      & 3.0376708717485E+00_wp, 5.0864296817319E+00_wp, 3.0581210229485E+00_wp, &
      & 3.0376410505379E+00_wp, 6.5103995227002E+00_wp, 7.7339832042333E+00_wp, &
      & 3.0717760845559E+00_wp, 3.0581487524093E+00_wp, 1.7179411037999E+01_wp, &
      & 1.0089722744091E+01_wp, 1.0521470167813E+01_wp, 7.7339283801084E+00_wp, &
      & 2.0820557350817E+01_wp, 1.8172837654743E+01_wp, 5.1002826817488E+00_wp, &
      & 1.1626599584016E+01_wp, 5.0864296817319E+00_wp, 9.4188339482987E+00_wp, &
      & 5.1165814010224E+00_wp, 5.0863857133152E+00_wp, 1.1573411106197E+01_wp, &
      & 1.3364030360454E+01_wp, 5.1367144347811E+00_wp, 5.1166222853615E+00_wp, &
      & 2.9120203955066E+01_wp, 1.6851140998078E+01_wp, 1.7468396982731E+01_wp, &
      & 1.3363934948055E+01_wp, 3.2809103726717E+01_wp, 1.3407379807008E+01_wp, &
      & 3.0675915503404E+00_wp, 6.5775401434039E+00_wp, 3.0581210229485E+00_wp, &
      & 5.1165814010224E+00_wp, 3.0787340802721E+00_wp, 3.0580909641819E+00_wp, &
      & 6.5512279033320E+00_wp, 7.7843558458146E+00_wp, 3.0924979182465E+00_wp, &
      & 3.0787620306260E+00_wp, 1.7293547030655E+01_wp, 1.0158056814792E+01_wp, &
      & 1.0593334527578E+01_wp, 7.7843006674916E+00_wp, 2.0974223126583E+01_wp, &
      & 1.3286079405123E+01_wp, 3.0470366226884E+00_wp, 6.5365120296846E+00_wp, &
      & 3.0376410505379E+00_wp, 5.0863857133152E+00_wp, 3.0580909641819E+00_wp, &
      & 3.0376112296738E+00_wp, 6.5103399851583E+00_wp, 7.7339097488767E+00_wp, &
      & 3.0717458671673E+00_wp, 3.0581186933205E+00_wp, 1.7179244600430E+01_wp, &
      & 1.0089623096676E+01_wp, 1.0521365372392E+01_wp, 7.7338549252683E+00_wp, &
      & 2.0820333269371E+01_wp, 2.5572059731044E+01_wp, 6.5291578413266E+00_wp, &
      & 1.4518124333712E+01_wp, 6.5103995227002E+00_wp, 1.1573411106197E+01_wp, &
      & 6.5512279033320E+00_wp, 6.5103399851583E+00_wp, 1.4455069876704E+01_wp, &
      & 1.6886749979232E+01_wp, 6.5784900026390E+00_wp, 6.5512832647310E+00_wp, &
      & 3.7083479708919E+01_wp, 2.1593323739087E+01_wp, 2.2440726387215E+01_wp, &
      & 1.6886629834681E+01_wp, 4.3119210273214E+01_wp, 3.2413191615239E+01_wp, &
      & 7.7571265678988E+00_wp, 1.6957726525989E+01_wp, 7.7339832042333E+00_wp, &
      & 1.3364030360454E+01_wp, 7.7843558458146E+00_wp, 7.7339097488767E+00_wp, &
      & 1.6886749979232E+01_wp, 1.9888865465619E+01_wp, 7.8179908794822E+00_wp, &
      & 7.7844241487907E+00_wp, 4.3903940006890E+01_wp, 2.5674131507899E+01_wp, &
      & 2.6728387684242E+01_wp, 1.9888724314305E+01_wp, 5.2175766976938E+01_wp, &
      & 1.3488257086895E+01_wp, 3.0812965885013E+00_wp, 6.6048957183646E+00_wp, &
      & 3.0717760845559E+00_wp, 5.1367144347811E+00_wp, 3.0924979182465E+00_wp, &
      & 3.0717458671673E+00_wp, 6.5784900026390E+00_wp, 7.8179908794822E+00_wp, &
      & 3.1063343888377E+00_wp, 3.0925260160962E+00_wp, 1.7369758399270E+01_wp, &
      & 1.0203685130015E+01_wp, 1.0641320102372E+01_wp, 7.8179354646525E+00_wp, &
      & 2.1076829489466E+01_wp, 1.3407544045253E+01_wp, 3.0676193812890E+00_wp, &
      & 6.5775956946247E+00_wp, 3.0581487524093E+00_wp, 5.1166222853615E+00_wp, &
      & 3.0787620306260E+00_wp, 3.0581186933205E+00_wp, 6.5512832647310E+00_wp, &
      & 7.7844241487907E+00_wp, 3.0925260160962E+00_wp, 3.0787899812795E+00_wp, &
      & 1.7293701793791E+01_wp, 1.0158149472637E+01_wp, 1.0593431972331E+01_wp, &
      & 7.7843689699874E+00_wp, 2.0974431490280E+01_wp, 7.3738556418539E+01_wp, &
      & 1.7231850034693E+01_wp, 3.7235277639546E+01_wp, 1.7179411037999E+01_wp, &
      & 2.9120203955066E+01_wp, 1.7293547030655E+01_wp, 1.7179244600430E+01_wp, &
      & 3.7083479708919E+01_wp, 4.3903940006890E+01_wp, 1.7369758399270E+01_wp, &
      & 1.7293701793791E+01_wp, 9.7300368597444E+01_wp, 5.7048187086676E+01_wp, &
      & 5.9450286555495E+01_wp, 4.3903628523784E+01_wp, 1.1700185746899E+02_wp, &
      & 4.4740018054061E+01_wp, 1.0121118362811E+01_wp, 2.1680017664472E+01_wp, &
      & 1.0089722744091E+01_wp, 1.6851140998078E+01_wp, 1.0158056814792E+01_wp, &
      & 1.0089623096676E+01_wp, 2.1593323739087E+01_wp, 2.5674131507899E+01_wp, &
      & 1.0203685130015E+01_wp, 1.0158149472637E+01_wp, 5.7048187086676E+01_wp, &
      & 3.3522232984596E+01_wp, 3.4965993468169E+01_wp, 2.5673949524761E+01_wp, &
      & 6.9388660920841E+01_wp, 4.7377208046689E+01_wp, 1.0554487753558E+01_wp, &
      & 2.2530197721674E+01_wp, 1.0521470167813E+01_wp, 1.7468396982731E+01_wp, &
      & 1.0593334527578E+01_wp, 1.0521365372392E+01_wp, 2.2440726387215E+01_wp, &
      & 2.6728387684242E+01_wp, 1.0641320102372E+01_wp, 1.0593431972331E+01_wp, &
      & 5.9450286555495E+01_wp, 3.4965993468169E+01_wp, 3.6486691424878E+01_wp, &
      & 2.6728198279961E+01_wp, 7.2679869961361E+01_wp, 3.2412969775949E+01_wp, &
      & 7.7570715810400E+00_wp, 1.6957605864802E+01_wp, 7.7339283801084E+00_wp, &
      & 1.3363934948055E+01_wp, 7.7843006674916E+00_wp, 7.7338549252683E+00_wp, &
      & 1.6886629834681E+01_wp, 1.9888724314305E+01_wp, 7.8179354646525E+00_wp, &
      & 7.7843689699874E+00_wp, 4.3903628523784E+01_wp, 2.5673949524761E+01_wp, &
      & 2.6728198279961E+01_wp, 1.9888583164001E+01_wp, 5.2175398249584E+01_wp, &
      & 1.0847087588239E+02_wp, 2.0891158034439E+01_wp, 4.3281495531031E+01_wp, &
      & 2.0820557350817E+01_wp, 3.2809103726717E+01_wp, 2.0974223126583E+01_wp, &
      & 2.0820333269371E+01_wp, 4.3119210273214E+01_wp, 5.2175766976938E+01_wp, &
      & 2.1076829489466E+01_wp, 2.0974431490280E+01_wp, 1.1700185746899E+02_wp, &
      & 6.9388660920841E+01_wp, 7.2679869961361E+01_wp, 5.2175398249584E+01_wp, &
      & 1.5002113554136E+02_wp], [16, 16])

   call get_structure(mol, "MB16-43", "01")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) then 
      call test_failed(error, "D4 model could not be created")
      return
   end if
   call test_c6_gen(error, mol, d4, ref, with_cn=.true., with_q=.false.)

end subroutine test_c6_d4_mb01

subroutine test_c6_d4s_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4s_model) :: d4s
   real(wp), parameter :: ref(16, 16) = reshape([&
      & 2.4365306222463E+02_wp, 2.0332338587140E+01_wp, 3.0502619153371E+01_wp, &
      & 1.8298596791412E+01_wp, 2.0742419322199E+01_wp, 2.0716929307831E+01_wp, &
      & 1.8190372677731E+01_wp, 3.0488441450911E+01_wp, 4.0457627599127E+01_wp, &
      & 2.0967433681204E+01_wp, 2.0717596735073E+01_wp, 9.9580698415550E+01_wp, &
      & 6.1846685567225E+01_wp, 7.1021885273123E+01_wp, 4.0423599454699E+01_wp, &
      & 1.8347543119346E+02_wp, 2.0332338587140E+01_wp, 3.1822696668068E+00_wp, &
      & 6.6304723933976E+00_wp, 3.1140882425890E+00_wp, 5.1231960123924E+00_wp, &
      & 3.2264928164160E+00_wp, 3.1133031099925E+00_wp, 6.6045195385247E+00_wp, &
      & 7.9287997139231E+00_wp, 3.2692094735311E+00_wp, 3.2265895526420E+00_wp, &
      & 1.8012665879686E+01_wp, 1.0906617349170E+01_wp, 1.1497231041362E+01_wp, &
      & 7.9271978729409E+00_wp, 2.3992975724594E+01_wp, 3.0502619153371E+01_wp, &
      & 6.6304723933976E+00_wp, 1.4581513138206E+01_wp, 6.5569791757550E+00_wp, &
      & 1.1626599584016E+01_wp, 6.6886761622036E+00_wp, 6.5564027372405E+00_wp, &
      & 1.4518124333712E+01_wp, 1.6973077732097E+01_wp, 6.7494612462405E+00_wp, &
      & 6.6888091930108E+00_wp, 3.7235277639546E+01_wp, 2.1774193034488E+01_wp, &
      & 2.2669077958442E+01_wp, 1.6972811391289E+01_wp, 4.3476182851475E+01_wp, &
      & 1.8298596791412E+01_wp, 3.1140882425890E+00_wp, 6.5569791757550E+00_wp, &
      & 3.0476210856053E+00_wp, 5.0902564954875E+00_wp, 3.1571995014269E+00_wp, &
      & 3.0468556933868E+00_wp, 6.5312738981987E+00_wp, 7.7505874377698E+00_wp, &
      & 3.1988421451067E+00_wp, 3.1572938054399E+00_wp, 1.7286653842521E+01_wp, &
      & 1.0253209303737E+01_wp, 1.0801256909450E+01_wp, 7.7490213169745E+00_wp, &
      & 2.1815226920390E+01_wp, 2.0742419322199E+01_wp, 5.1231960123924E+00_wp, &
      & 1.1626599584016E+01_wp, 5.0902564954875E+00_wp, 9.4188339482987E+00_wp, &
      & 5.1535678750056E+00_wp, 5.0900696911622E+00_wp, 1.1573411106197E+01_wp, &
      & 1.3364030360454E+01_wp, 5.1873221228506E+00_wp, 5.1536397442244E+00_wp, &
      & 2.9120203955066E+01_wp, 1.6877692257805E+01_wp, 1.7511681863102E+01_wp, &
      & 1.3363934948055E+01_wp, 3.2885793756667E+01_wp, 2.0716929307831E+01_wp, &
      & 3.2264928164160E+00_wp, 6.6886761622036E+00_wp, 3.1571995014269E+00_wp, &
      & 5.1535678750056E+00_wp, 3.2714371494340E+00_wp, 3.1564015650249E+00_wp, &
      & 6.6625272357290E+00_wp, 8.0408657246332E+00_wp, 3.3148504223411E+00_wp, &
      & 3.2715354632174E+00_wp, 1.8375751583069E+01_wp, 1.1176495233369E+01_wp, &
      & 1.1784690045716E+01_wp, 8.0392414215684E+00_wp, 2.4688315344031E+01_wp, &
      & 1.8190372677731E+01_wp, 3.1133031099925E+00_wp, 6.5564027372405E+00_wp, &
      & 3.0468556933868E+00_wp, 5.0900696911622E+00_wp, 3.1564015650249E+00_wp, &
      & 3.0460905284855E+00_wp, 6.5306994015382E+00_wp, 7.7484183050892E+00_wp, &
      & 3.1980318410700E+00_wp, 3.1564958410302E+00_wp, 1.7273782041288E+01_wp, &
      & 1.0238064084132E+01_wp, 1.0785125058701E+01_wp, 7.7468526190666E+00_wp, &
      & 2.1747178438015E+01_wp, 3.0488441450911E+01_wp, 6.6045195385247E+00_wp, &
      & 1.4518124333712E+01_wp, 6.5312738981987E+00_wp, 1.1573411106197E+01_wp, &
      & 6.6625272357290E+00_wp, 6.5306994015382E+00_wp, 1.4455069876704E+01_wp, &
      & 1.6901253376081E+01_wp, 6.7231075524626E+00_wp, 6.6626598183941E+00_wp, &
      & 3.7083479708919E+01_wp, 2.1697223691340E+01_wp, 2.2589483264101E+01_wp, &
      & 1.6900988187810E+01_wp, 4.3337080992823E+01_wp, 4.0457627599127E+01_wp, &
      & 7.9287997139231E+00_wp, 1.6973077732097E+01_wp, 7.7505874377698E+00_wp, &
      & 1.3364030360454E+01_wp, 8.0408657246332E+00_wp, 7.7484183050892E+00_wp, &
      & 1.6901253376081E+01_wp, 1.9903312598032E+01_wp, 8.1478434370286E+00_wp, &
      & 8.0411092346947E+00_wp, 4.3913859005520E+01_wp, 2.5705612403693E+01_wp, &
      & 2.6926314672038E+01_wp, 1.9897989085535E+01_wp, 5.2296415572987E+01_wp, &
      & 2.0967433681204E+01_wp, 3.2692094735311E+00_wp, 6.7494612462405E+00_wp, &
      & 3.1988421451067E+00_wp, 5.1873221228506E+00_wp, 3.3148504223411E+00_wp, &
      & 3.1980318410700E+00_wp, 6.7231075524626E+00_wp, 8.1478434370286E+00_wp, &
      & 3.3589365803364E+00_wp, 3.3149502599414E+00_wp, 1.8693867510742E+01_wp, &
      & 1.1396899314516E+01_wp, 1.2019452289797E+01_wp, 8.1461976917604E+00_wp, &
      & 2.5207079835493E+01_wp, 2.0717596735073E+01_wp, 3.2265895526420E+00_wp, &
      & 6.6888091930108E+00_wp, 3.1572938054399E+00_wp, 5.1536397442244E+00_wp, &
      & 3.2715354632174E+00_wp, 3.1564958410302E+00_wp, 6.6626598183941E+00_wp, &
      & 8.0411092346947E+00_wp, 3.3149502599414E+00_wp, 3.2716337804517E+00_wp, &
      & 1.8376503071173E+01_wp, 1.1177031770379E+01_wp, 1.1785261535282E+01_wp, &
      & 8.0394848828217E+00_wp, 2.4689626643585E+01_wp, 9.9580698415550E+01_wp, &
      & 1.8012665879686E+01_wp, 3.7235277639546E+01_wp, 1.7286653842521E+01_wp, &
      & 2.9120203955066E+01_wp, 1.8375751583069E+01_wp, 1.7273782041288E+01_wp, &
      & 3.7083479708919E+01_wp, 4.3913859005520E+01_wp, 1.8693867510742E+01_wp, &
      & 1.8376503071173E+01_wp, 9.7300368597444E+01_wp, 5.7050336368043E+01_wp, &
      & 5.9923299874942E+01_wp, 4.3905285101770E+01_wp, 1.1697722659238E+02_wp, &
      & 6.1846685567225E+01_wp, 1.0906617349170E+01_wp, 2.1774193034488E+01_wp, &
      & 1.0253209303737E+01_wp, 1.6877692257805E+01_wp, 1.1176495233369E+01_wp, &
      & 1.0238064084132E+01_wp, 2.1697223691340E+01_wp, 2.5705612403693E+01_wp, &
      & 1.1396899314516E+01_wp, 1.1177031770379E+01_wp, 5.7050336368043E+01_wp, &
      & 3.3659827509442E+01_wp, 3.5934002503412E+01_wp, 2.5688699926209E+01_wp, &
      & 7.0354767087578E+01_wp, 7.1021885273123E+01_wp, 1.1497231041362E+01_wp, &
      & 2.2669077958442E+01_wp, 1.0801256909450E+01_wp, 1.7511681863102E+01_wp, &
      & 1.1784690045716E+01_wp, 1.0785125058701E+01_wp, 2.2589483264101E+01_wp, &
      & 2.6926314672038E+01_wp, 1.2019452289797E+01_wp, 1.1785261535282E+01_wp, &
      & 5.9923299874942E+01_wp, 3.5934002503412E+01_wp, 3.8409602775833E+01_wp, &
      & 2.6908606755692E+01_wp, 7.6809703882464E+01_wp, 4.0423599454699E+01_wp, &
      & 7.9271978729409E+00_wp, 1.6972811391289E+01_wp, 7.7490213169745E+00_wp, &
      & 1.3363934948055E+01_wp, 8.0392414215684E+00_wp, 7.7468526190666E+00_wp, &
      & 1.6900988187810E+01_wp, 1.9897989085535E+01_wp, 8.1461976917604E+00_wp, &
      & 8.0394848828217E+00_wp, 4.3905285101770E+01_wp, 2.5688699926209E+01_wp, &
      & 2.6908606755692E+01_wp, 1.9892667009015E+01_wp, 5.2233988011324E+01_wp, &
      & 1.8347543119346E+02_wp, 2.3992975724594E+01_wp, 4.3476182851475E+01_wp, &
      & 2.1815226920390E+01_wp, 3.2885793756667E+01_wp, 2.4688315344031E+01_wp, &
      & 2.1747178438015E+01_wp, 4.3337080992823E+01_wp, 5.2296415572987E+01_wp, &
      & 2.5207079835493E+01_wp, 2.4689626643585E+01_wp, 1.1697722659238E+02_wp, &
      & 7.0354767087578E+01_wp, 7.6809703882464E+01_wp, 5.2233988011324E+01_wp, &
      & 1.5010894587960E+02_wp], [16, 16])

   call get_structure(mol, "MB16-43", "01")
   call new_d4s_model(error, d4s, mol)
   if (allocated(error)) then 
      call test_failed(error, "D4S model could not be created")
      return
   end if
   call test_c6_gen(error, mol, d4s, ref, with_cn=.true., with_q=.false.)

end subroutine test_c6_d4s_mb01

subroutine test_c6_d4_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   real(wp), parameter :: ref(16, 16) = reshape([&
      & 3.0589395270838E+00_wp, 2.0548666062381E+01_wp, 1.0486006240981E+01_wp, &
      & 6.5056567344045E+00_wp, 1.4992082319852E+01_wp, 3.0801316359883E+00_wp, &
      & 3.0579827712864E+00_wp, 3.0577825167623E+00_wp, 2.1537656255928E+01_wp, &
      & 3.0803284996092E+00_wp, 1.0172957112767E+01_wp, 5.2220400277510E+00_wp, &
      & 5.0976789500531E+00_wp, 3.0600561665054E+00_wp, 3.0481646087415E+00_wp, &
      & 2.0549746985773E+01_wp, 2.0548666062381E+01_wp, 1.3851681015795E+02_wp, &
      & 7.0565183559524E+01_wp, 4.3366575795645E+01_wp, 1.0382773857053E+02_wp, &
      & 2.0694474784198E+01_wp, 2.0542083265900E+01_wp, 2.0540705448509E+01_wp, &
      & 1.4602962434315E+02_wp, 2.0695829271052E+01_wp, 6.8403707744864E+01_wp, &
      & 3.5378883322712E+01_wp, 3.3779518127935E+01_wp, 2.0556348911102E+01_wp, &
      & 2.0474531058723E+01_wp, 1.3852545717271E+02_wp, 1.0486006240981E+01_wp, &
      & 7.0565183559524E+01_wp, 3.5978462430120E+01_wp, 2.2210612243768E+01_wp, &
      & 5.2203292609600E+01_wp, 1.0559559576149E+01_wp, 1.0482685543632E+01_wp, &
      & 1.0481990502464E+01_wp, 7.4184741679162E+01_wp, 1.0560242848209E+01_wp, &
      & 3.4889842830742E+01_wp, 1.7979649047208E+01_wp, 1.7349677314061E+01_wp, &
      & 1.0489881860628E+01_wp, 1.0448608774602E+01_wp, 7.0569262061317E+01_wp, &
      & 6.5056567344045E+00_wp, 4.3366575795645E+01_wp, 2.2210612243768E+01_wp, &
      & 1.4330344748644E+01_wp, 2.9585795587556E+01_wp, 6.5476501757141E+00_wp, &
      & 6.5037608651161E+00_wp, 6.5033640487047E+00_wp, 4.4690011014072E+01_wp, &
      & 6.5480402728463E+00_wp, 2.1600731051093E+01_wp, 1.0877146737095E+01_wp, &
      & 1.1511272048731E+01_wp, 6.5078694227322E+00_wp, 6.4843055841691E+00_wp, &
      & 4.3367198644475E+01_wp, 1.4992082319852E+01_wp, 1.0382773857053E+02_wp, &
      & 5.2203292609600E+01_wp, 2.9585795587556E+01_wp, 9.4034692488024E+01_wp, &
      & 1.5118842906512E+01_wp, 1.4986359485625E+01_wp, 1.4985161663310E+01_wp, &
      & 1.1441184252246E+02_wp, 1.5120020446142E+01_wp, 5.0280957257019E+01_wp, &
      & 2.7544654796568E+01_wp, 2.1792156965418E+01_wp, 1.4998761497880E+01_wp, &
      & 1.4927632152271E+01_wp, 1.0384223027605E+02_wp, 3.0801316359883E+00_wp, &
      & 2.0694474784198E+01_wp, 1.0559559576149E+01_wp, 6.5476501757141E+00_wp, &
      & 1.5118842906512E+01_wp, 3.1014972488871E+00_wp, 3.0791670470416E+00_wp, &
      & 3.0789651529940E+00_wp, 2.1697037089362E+01_wp, 3.1016957242694E+00_wp, &
      & 1.0243880281179E+01_wp, 5.2604271712471E+00_wp, 5.1287635630518E+00_wp, &
      & 3.0812574175580E+00_wp, 3.0692685012547E+00_wp, 2.0695574926403E+01_wp, &
      & 3.0579827712864E+00_wp, 2.0542083265900E+01_wp, 1.0482685543632E+01_wp, &
      & 6.5037608651161E+00_wp, 1.4986359485625E+01_wp, 3.0791670470416E+00_wp, &
      & 3.0570263691305E+00_wp, 3.0568261886256E+00_wp, 2.1530460722105E+01_wp, &
      & 3.0793638378966E+00_wp, 1.0169755159019E+01_wp, 5.2203069712509E+00_wp, &
      & 5.0962755794007E+00_wp, 3.0590989979692E+00_wp, 3.0472118356314E+00_wp, &
      & 2.0543163321624E+01_wp, 3.0577825167623E+00_wp, 2.0540705448509E+01_wp, &
      & 1.0481990502464E+01_wp, 6.5033640487047E+00_wp, 1.4985161663310E+01_wp, &
      & 3.0789651529940E+00_wp, 3.0568261886256E+00_wp, 3.0566260236134E+00_wp, &
      & 2.1528954655244E+01_wp, 3.0791619286187E+00_wp, 1.0169084971573E+01_wp, &
      & 5.2199442324981E+00_wp, 5.0959818458034E+00_wp, 3.0588986570565E+00_wp, &
      & 3.0470124147068E+00_wp, 2.0541785322625E+01_wp, 2.1537656255928E+01_wp, &
      & 1.4602962434315E+02_wp, 7.4184741679162E+01_wp, 4.4690011014072E+01_wp, &
      & 1.1441184252246E+02_wp, 2.1697037089362E+01_wp, 2.1530460722105E+01_wp, &
      & 2.1528954655244E+01_wp, 1.5556225735539E+02_wp, 2.1698517654043E+01_wp, &
      & 7.1805295866219E+01_wp, 3.7626982953065E+01_wp, 3.4358625569838E+01_wp, &
      & 2.1546054236631E+01_wp, 2.1456620646243E+01_wp, 1.4604162267121E+02_wp, &
      & 3.0803284996092E+00_wp, 2.0695829271052E+01_wp, 1.0560242848209E+01_wp, &
      & 6.5480402728463E+00_wp, 1.5120020446142E+01_wp, 3.1016957242694E+00_wp, &
      & 3.0793638378966E+00_wp, 3.0791619286187E+00_wp, 2.1698517654043E+01_wp, &
      & 3.1018942146242E+00_wp, 1.0244539120363E+01_wp, 5.2607837677566E+00_wp, &
      & 5.1290523228678E+00_wp, 3.0814543661047E+00_wp, 3.0694645453914E+00_wp, &
      & 2.0696929591790E+01_wp, 1.0172957112767E+01_wp, 6.8403707744864E+01_wp, &
      & 3.4889842830742E+01_wp, 2.1600731051093E+01_wp, 5.0280957257019E+01_wp, &
      & 1.0243880281179E+01_wp, 1.0169755159019E+01_wp, 1.0169084971573E+01_wp, &
      & 7.1805295866219E+01_wp, 1.0244539120363E+01_wp, 3.3841320400868E+01_wp, &
      & 1.7407271135443E+01_wp, 1.6904472228624E+01_wp, 1.0176694145552E+01_wp, &
      & 1.0136896928684E+01_wp, 6.8407464235121E+01_wp, 5.2220400277510E+00_wp, &
      & 3.5378883322712E+01_wp, 1.7979649047208E+01_wp, 1.0877146737095E+01_wp, &
      & 2.7544654796568E+01_wp, 5.2604271712471E+00_wp, 5.2203069712509E+00_wp, &
      & 5.2199442324981E+00_wp, 3.7626982953065E+01_wp, 5.2607837677566E+00_wp, &
      & 1.7407271135443E+01_wp, 9.1050543972585E+00_wp, 8.3876487926031E+00_wp, &
      & 5.2240626956176E+00_wp, 5.2025224639999E+00_wp, 3.5381660815106E+01_wp, &
      & 5.0976789500531E+00_wp, 3.3779518127935E+01_wp, 1.7349677314061E+01_wp, &
      & 1.1511272048731E+01_wp, 2.1792156965418E+01_wp, 5.1287635630518E+00_wp, &
      & 5.0962755794007E+00_wp, 5.0959818458034E+00_wp, 3.4358625569838E+01_wp, &
      & 5.1290523228678E+00_wp, 1.6904472228624E+01_wp, 8.3876487926031E+00_wp, &
      & 9.4043400097466E+00_wp, 5.0993168382139E+00_wp, 5.0818742857664E+00_wp, &
      & 3.3779053209122E+01_wp, 3.0600561665054E+00_wp, 2.0556348911102E+01_wp, &
      & 1.0489881860628E+01_wp, 6.5078694227322E+00_wp, 1.4998761497880E+01_wp, &
      & 3.0812574175580E+00_wp, 3.0590989979692E+00_wp, 3.0588986570565E+00_wp, &
      & 2.1546054236631E+01_wp, 3.0814543661047E+00_wp, 1.0176694145552E+01_wp, &
      & 5.2240626956176E+00_wp, 5.0993168382139E+00_wp, 3.0611732876384E+00_wp, &
      & 3.0492765999282E+00_wp, 2.0557430847157E+01_wp, 3.0481646087415E+00_wp, &
      & 2.0474531058723E+01_wp, 1.0448608774602E+01_wp, 6.4843055841691E+00_wp, &
      & 1.4927632152271E+01_wp, 3.0692685012547E+00_wp, 3.0472118356314E+00_wp, &
      & 3.0470124147068E+00_wp, 2.1456620646243E+01_wp, 3.0694645453914E+00_wp, &
      & 1.0136896928684E+01_wp, 5.2025224639999E+00_wp, 5.0818742857664E+00_wp, &
      & 3.0492765999282E+00_wp, 3.0374345431516E+00_wp, 2.0475602210498E+01_wp, &
      & 2.0549746985773E+01_wp, 1.3852545717271E+02_wp, 7.0569262061317E+01_wp, &
      & 4.3367198644475E+01_wp, 1.0384223027605E+02_wp, 2.0695574926403E+01_wp, &
      & 2.0543163321624E+01_wp, 2.0541785322625E+01_wp, 1.4604162267121E+02_wp, &
      & 2.0696929591790E+01_wp, 6.8407464235121E+01_wp, 3.5381660815106E+01_wp, &
      & 3.3779053209122E+01_wp, 2.0557430847157E+01_wp, 2.0475602210498E+01_wp, &
      & 1.3853411065714E+02_wp], [16, 16])

   call get_structure(mol, "MB16-43", "02")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) then 
      call test_failed(error, "D4 model could not be created")
      return
   end if
   call test_c6_gen(error, mol, d4, ref, with_cn=.true., with_q=.false.)

end subroutine test_c6_d4_mb02

subroutine test_c6_d4_eeqbc_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   real(wp), parameter :: ref(16, 16) = reshape([&
      1.0570679100062E+01_wp, 3.8815647635771E+01_wp, 3.3698374113141E+01_wp, &
      1.3727374704448E+01_wp, 7.4074408768038E+01_wp, 6.6257190381866E+00_wp, &
      9.3685106345191E+00_wp, 9.8654303943211E+00_wp, 5.5985549225446E+01_wp, &
      6.8968671962915E+00_wp, 3.3119557687207E+01_wp, 7.5961669722395E+01_wp, &
      1.0453095328420E+01_wp, 9.4909328258469E+00_wp, 1.0152147444406E+01_wp, &
      3.9321967045592E+01_wp, 3.8815647635771E+01_wp, 1.4365097800474E+02_wp, &
      1.2243745731009E+02_wp, 5.1950152544985E+01_wp, 2.6302723286397E+02_wp, &
      2.4329711751288E+01_wp, 3.4401272067689E+01_wp, 3.6225966783814E+01_wp, &
      2.0299287006316E+02_wp, 2.5325370711556E+01_wp, 1.2033442371883E+02_wp, &
      2.6459567906515E+02_wp, 3.9952908598325E+01_wp, 3.4850807674286E+01_wp, &
      3.7278794883305E+01_wp, 1.4552479134462E+02_wp, 3.3698374113141E+01_wp, &
      1.2243745731009E+02_wp, 1.0936170613305E+02_wp, 4.2120145790738E+01_wp, &
      2.5058609065255E+02_wp, 2.1122196294471E+01_wp, 2.9865969182918E+01_wp, &
      3.1450104678046E+01_wp, 1.8230632688512E+02_wp, 2.1986592232084E+01_wp, &
      1.0748326683492E+02_wp, 2.6767061141066E+02_wp, 3.1720015995405E+01_wp, &
      3.0256240116702E+01_wp, 3.2364132842832E+01_wp, 1.2403455706090E+02_wp, &
      1.3727374704448E+01_wp, 5.1950152544985E+01_wp, 4.2120145790738E+01_wp, &
      2.0044123974403E+01_wp, 8.5203266890613E+01_wp, 8.6043410326444E+00_wp, &
      1.2166205660514E+01_wp, 1.2811519332064E+01_wp, 6.9434919615275E+01_wp, &
      8.9564614907052E+00_wp, 4.1396673676707E+01_wp, 8.1494900263643E+01_wp, &
      1.5860100108505E+01_wp, 1.2325186486304E+01_wp, 1.3183857981588E+01_wp, &
      5.2627801175610E+01_wp, 7.4074408768038E+01_wp, 2.6302723286397E+02_wp, &
      2.5058609065255E+02_wp, 8.5203266890613E+01_wp, 6.3119494197519E+02_wp, &
      4.6429961194634E+01_wp, 6.5650170601150E+01_wp, 6.9132353445228E+01_wp, &
      4.2082340967655E+02_wp, 4.8330041531010E+01_wp, 2.4628192650508E+02_wp, &
      7.4319763495018E+02_wp, 6.2620932975588E+01_wp, 6.6508048449570E+01_wp, &
      7.1141533344429E+01_wp, 2.6645821499576E+02_wp, 6.6257190381866E+00_wp, &
      2.4329711751288E+01_wp, 2.1122196294471E+01_wp, 8.6043410326444E+00_wp, &
      4.6429961194634E+01_wp, 4.1530115858622E+00_wp, 5.8721978676113E+00_wp, &
      6.1836679899946E+00_wp, 3.5091834295130E+01_wp, 4.3229676971317E+00_wp, &
      2.0759393207118E+01_wp, 4.7612899463495E+01_wp, 6.5520173368126E+00_wp, &
      5.9489322983985E+00_wp, 6.3633827083525E+00_wp, 2.4647073589756E+01_wp, &
      9.3685106345191E+00_wp, 3.4401272067689E+01_wp, 2.9865969182918E+01_wp, &
      1.2166205660514E+01_wp, 6.5650170601150E+01_wp, 5.8721978676113E+00_wp, &
      8.3030608230821E+00_wp, 8.7434675377445E+00_wp, 4.9618497386313E+01_wp, &
      6.1125092401059E+00_wp, 2.9352979639818E+01_wp, 6.7322799592498E+01_wp, &
      9.2642992773630E+00_wp, 8.4115603424129E+00_wp, 8.9975772034897E+00_wp, &
      3.4850009441179E+01_wp, 9.8654303943211E+00_wp, 3.6225966783814E+01_wp, &
      3.1450104678046E+01_wp, 1.2811519332064E+01_wp, 6.9132353445228E+01_wp, &
      6.1836679899946E+00_wp, 8.7434675377445E+00_wp, 9.2072340806018E+00_wp, &
      5.2250336401588E+01_wp, 6.4367258356648E+00_wp, 3.0909905405409E+01_wp, &
      7.0893701169900E+01_wp, 9.7556915115436E+00_wp, 8.8577220332062E+00_wp, &
      9.4748221015511E+00_wp, 3.6698505856052E+01_wp, 5.5985549225446E+01_wp, &
      2.0299287006316E+02_wp, 1.8230632688512E+02_wp, 6.9434919615275E+01_wp, &
      4.2082340967655E+02_wp, 3.5091834295130E+01_wp, 4.9618497386313E+01_wp, &
      5.2250336401588E+01_wp, 3.0411027771033E+02_wp, 3.6527917862636E+01_wp, &
      1.7917496235998E+02_wp, 4.5248422961037E+02_wp, 5.2171208291926E+01_wp, &
      5.0266882750583E+01_wp, 5.3768877582277E+01_wp, 2.0564075143356E+02_wp, &
      6.8968671962915E+00_wp, 2.5325370711556E+01_wp, 2.1986592232084E+01_wp, &
      8.9564614907052E+00_wp, 4.8330041531010E+01_wp, 4.3229676971317E+00_wp, &
      6.1125092401059E+00_wp, 6.4367258356648E+00_wp, 3.6527917862636E+01_wp, &
      4.4998790213018E+00_wp, 2.1608941942741E+01_wp, 4.9561389871427E+01_wp, &
      6.8201493572785E+00_wp, 6.1923839186838E+00_wp, 6.6237951240831E+00_wp, &
      2.5655720133327E+01_wp, 3.3119557687207E+01_wp, 1.2033442371883E+02_wp, &
      1.0748326683492E+02_wp, 4.1396673676707E+01_wp, 2.4628192650508E+02_wp, &
      2.0759393207118E+01_wp, 2.9352979639818E+01_wp, 3.0909905405409E+01_wp, &
      1.7917496235998E+02_wp, 2.1608941942741E+01_wp, 1.0563709234246E+02_wp, &
      2.6307299707896E+02_wp, 3.1175180569910E+01_wp, 2.9736547127724E+01_wp, &
      3.1808233865699E+01_wp, 1.2190409106049E+02_wp, 7.5961669722395E+01_wp, &
      2.6459567906515E+02_wp, 2.6767061141066E+02_wp, 8.1494900263643E+01_wp, &
      7.4319763495018E+02_wp, 4.7612899463495E+01_wp, 6.7322799592498E+01_wp, &
      7.0893701169900E+01_wp, 4.5248422961037E+02_wp, 4.9561389871427E+01_wp, &
      2.6307299707896E+02_wp, 1.0050045877384E+03_wp, 5.8656870507959E+01_wp, &
      6.8202534373614E+01_wp, 7.2954070769255E+01_wp, 2.6804712033204E+02_wp, &
      1.0453095328420E+01_wp, 3.9952908598325E+01_wp, 3.1720015995405E+01_wp, &
      1.5860100108505E+01_wp, 6.2620932975588E+01_wp, 6.5520173368126E+00_wp, &
      9.2642992773630E+00_wp, 9.7556915115436E+00_wp, 5.2171208291926E+01_wp, &
      6.8201493572785E+00_wp, 3.1175180569910E+01_wp, 5.8656870507959E+01_wp, &
      1.2706815830838E+01_wp, 9.3853596959180E+00_wp, 1.0039219242208E+01_wp, &
      4.0474062674538E+01_wp, 9.4909328258469E+00_wp, 3.4850807674286E+01_wp, &
      3.0256240116702E+01_wp, 1.2325186486304E+01_wp, 6.6508048449570E+01_wp, &
      5.9489322983985E+00_wp, 8.4115603424129E+00_wp, 8.8577220332062E+00_wp, &
      5.0266882750583E+01_wp, 6.1923839186838E+00_wp, 2.9736547127724E+01_wp, &
      6.8202534373614E+01_wp, 9.3853596959180E+00_wp, 8.5214776696997E+00_wp, &
      9.1151522547293E+00_wp, 3.5305408884060E+01_wp, 1.0152147444406E+01_wp, &
      3.7278794883305E+01_wp, 3.2364132842832E+01_wp, 1.3183857981588E+01_wp, &
      7.1141533344429E+01_wp, 6.3633827083525E+00_wp, 8.9975772034897E+00_wp, &
      9.4748221015511E+00_wp, 5.3768877582277E+01_wp, 6.6237951240831E+00_wp, &
      3.1808233865699E+01_wp, 7.2954070769255E+01_wp, 1.0039219242208E+01_wp, &
      9.1151522547293E+00_wp, 9.7501869801678E+00_wp, 3.7765067265033E+01_wp, &
      3.9321967045592E+01_wp, 1.4552479134462E+02_wp, 1.2403455706090E+02_wp, &
      5.2627801175610E+01_wp, 2.6645821499576E+02_wp, 2.4647073589756E+01_wp, &
      3.4850009441179E+01_wp, 3.6698505856052E+01_wp, 2.0564075143356E+02_wp, &
      2.5655720133327E+01_wp, 1.2190409106049E+02_wp, 2.6804712033204E+02_wp, &
      4.0474062674538E+01_wp, 3.5305408884060E+01_wp, 3.7765067265033E+01_wp, &
      1.4742304709681E+02_wp], [16, 16])

   call get_structure(mol, "MB16-43", "02")
   call new_d4_model(error, d4, mol)
   call new_d4_model(error, d4, mol, qmod=d4_qmod%eeqbc)
   if (allocated(error)) then 
      call test_failed(error, "D4 model could not be created")
      return
   end if
   call test_c6_gen(error, mol, d4, ref, with_cn=.false., with_q=.true.)

end subroutine test_c6_d4_eeqbc_mb02

subroutine test_c6_d4_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   real(wp), parameter :: ref(16, 16) = reshape([&
      & 2.6969645804632E+01_wp, 2.8729803918609E+01_wp, 8.0783557426273E+00_wp, &
      & 1.2572598284646E+01_wp, 3.9383726731882E+01_wp, 6.0100611314024E+01_wp, &
      & 2.8877446590256E+01_wp, 7.9285493229103E+00_wp, 7.5097624231199E+00_wp, &
      & 7.1772575509757E+00_wp, 1.9756649427858E+01_wp, 6.3114791359146E+01_wp, &
      & 2.9607901769799E+01_wp, 5.3279952989770E+00_wp, 3.4651507628618E+01_wp, &
      & 8.1388552398630E+00_wp, 2.8729803918609E+01_wp, 3.0858718804005E+01_wp, &
      & 8.5314526128999E+00_wp, 1.3158527195522E+01_wp, 4.0242276058964E+01_wp, &
      & 6.2518208475630E+01_wp, 3.0769885934381E+01_wp, 8.3738757317567E+00_wp, &
      & 7.9309682052198E+00_wp, 7.5798206446138E+00_wp, 2.1471517618889E+01_wp, &
      & 6.6436238381645E+01_wp, 3.1548229972467E+01_wp, 5.6258043424413E+00_wp, &
      & 3.4923486847229E+01_wp, 8.5961975589899E+00_wp, 8.0783557426273E+00_wp, &
      & 8.5314526128999E+00_wp, 2.4429152644641E+00_wp, 3.8470627309562E+00_wp, &
      & 1.2418251750013E+01_wp, 1.8522898371215E+01_wp, 8.6470455821784E+00_wp, &
      & 2.3973918803547E+00_wp, 2.2709711896641E+00_wp, 2.1704185011555E+00_wp, &
      & 5.7906292579429E+00_wp, 1.9170869616128E+01_wp, 8.8657656755301E+00_wp, &
      & 1.6115591315647E+00_wp, 1.1136757306622E+01_wp, 2.4609116715906E+00_wp, &
      & 1.2572598284646E+01_wp, 1.3158527195522E+01_wp, 3.8470627309562E+00_wp, &
      & 6.1873811861492E+00_wp, 2.1045559351132E+01_wp, 3.0123412184969E+01_wp, &
      & 1.3450271401073E+01_wp, 3.7748067440421E+00_wp, 3.5762880002650E+00_wp, &
      & 3.4179332169982E+00_wp, 8.7894143888807E+00_wp, 3.0447639291352E+01_wp, &
      & 1.3790466532352E+01_wp, 2.5387769376831E+00_wp, 1.9680032604563E+01_wp, &
      & 3.8746391243000E+00_wp, 3.9383726731882E+01_wp, 4.0242276058964E+01_wp, &
      & 1.2418251750013E+01_wp, 2.1045559351132E+01_wp, 8.0792669914257E+01_wp, &
      & 1.0524757517715E+02_wp, 4.2068874980499E+01_wp, 1.2180230990667E+01_wp, &
      & 1.1544195337928E+01_wp, 1.1032977672474E+01_wp, 2.5675572843387E+01_wp, &
      & 1.0047223990865E+02_wp, 4.3132752418890E+01_wp, 8.2028913059178E+00_wp, &
      & 8.3004740679372E+01_wp, 1.2500821517587E+01_wp, 6.0100611314024E+01_wp, &
      & 6.2518208475630E+01_wp, 1.8522898371215E+01_wp, 3.0123412184969E+01_wp, &
      & 1.0524757517715E+02_wp, 1.4754968408298E+02_wp, 6.4275900488174E+01_wp, &
      & 1.8173467936826E+01_wp, 1.7219167833650E+01_wp, 1.6456702330152E+01_wp, &
      & 4.1319144188437E+01_wp, 1.4726172884744E+02_wp, 6.5901567272212E+01_wp, &
      & 1.2226230276259E+01_wp, 1.0042992019072E+02_wp, 1.8653608543109E+01_wp, &
      & 2.8877446590256E+01_wp, 3.0769885934381E+01_wp, 8.6470455821784E+00_wp, &
      & 1.3450271401073E+01_wp, 4.2068874980499E+01_wp, 6.4275900488174E+01_wp, &
      & 3.0920640741304E+01_wp, 8.4867269265489E+00_wp, 8.0384251534349E+00_wp, &
      & 7.6825133853359E+00_wp, 2.1168508106272E+01_wp, 6.7543056296917E+01_wp, &
      & 3.1702779636090E+01_wp, 5.7030141029870E+00_wp, 3.6962950298764E+01_wp, &
      & 8.7118494250596E+00_wp, 7.9285493229103E+00_wp, 8.3738757317567E+00_wp, &
      & 2.3973918803547E+00_wp, 3.7748067440421E+00_wp, 1.2180230990667E+01_wp, &
      & 1.8173467936826E+01_wp, 8.4867269265489E+00_wp, 2.3527194049706E+00_wp, &
      & 2.2286519594791E+00_wp, 2.1299730839455E+00_wp, 5.6843960534634E+00_wp, &
      & 1.8812501010282E+01_wp, 8.7013919729876E+00_wp, 1.5815237509997E+00_wp, &
      & 1.0919772339795E+01_wp, 2.4150564158072E+00_wp, 7.5097624231199E+00_wp, &
      & 7.9309682052198E+00_wp, 2.2709711896641E+00_wp, 3.5762880002650E+00_wp, &
      & 1.1544195337928E+01_wp, 1.7219167833650E+01_wp, 8.0384251534349E+00_wp, &
      & 2.2286519594791E+00_wp, 2.1111293622441E+00_wp, 2.0176540534750E+00_wp, &
      & 5.3830571387835E+00_wp, 1.7821531942660E+01_wp, 8.2417506804364E+00_wp, &
      & 1.4981298823014E+00_wp, 1.0352898443990E+01_wp, 2.2877009236808E+00_wp, &
      & 7.1772575509757E+00_wp, 7.5798206446138E+00_wp, 2.1704185011555E+00_wp, &
      & 3.4179332169982E+00_wp, 1.1032977672474E+01_wp, 1.6456702330152E+01_wp, &
      & 7.6825133853359E+00_wp, 2.1299730839455E+00_wp, 2.0176540534750E+00_wp, &
      & 1.9283175881762E+00_wp, 5.1447273217295E+00_wp, 1.7032428915733E+01_wp, &
      & 7.8768364100925E+00_wp, 1.4317965486738E+00_wp, 9.8943978385702E+00_wp, &
      & 2.1864075234916E+00_wp, 1.9756649427858E+01_wp, 2.1471517618889E+01_wp, &
      & 5.7906292579429E+00_wp, 8.7894143888807E+00_wp, 2.5675572843387E+01_wp, &
      & 4.1319144188437E+01_wp, 2.1168508106272E+01_wp, 5.6843960534634E+00_wp, &
      & 5.3830571387835E+00_wp, 5.1447273217295E+00_wp, 1.5194417412576E+01_wp, &
      & 4.4825188609896E+01_wp, 2.1704001718253E+01_wp, 3.8172825287704E+00_wp, &
      & 2.1539554554163E+01_wp, 5.8355457277839E+00_wp, 6.3114791359146E+01_wp, &
      & 6.6436238381645E+01_wp, 1.9170869616128E+01_wp, 3.0447639291352E+01_wp, &
      & 1.0047223990865E+02_wp, 1.4726172884744E+02_wp, 6.7543056296917E+01_wp, &
      & 1.8812501010282E+01_wp, 1.7821531942660E+01_wp, 1.7032428915733E+01_wp, &
      & 4.4825188609896E+01_wp, 1.5096768893659E+02_wp, 6.9251466330288E+01_wp, &
      & 1.2648593650528E+01_wp, 9.1879170429485E+01_wp, 1.9310583681271E+01_wp, &
      & 2.9607901769799E+01_wp, 3.1548229972467E+01_wp, 8.8657656755301E+00_wp, &
      & 1.3790466532352E+01_wp, 4.3132752418890E+01_wp, 6.5901567272212E+01_wp, &
      & 3.1702779636090E+01_wp, 8.7013919729876E+00_wp, 8.2417506804364E+00_wp, &
      & 7.8768364100925E+00_wp, 2.1704001718253E+01_wp, 6.9251466330288E+01_wp, &
      & 3.2504702769580E+01_wp, 5.8472671394681E+00_wp, 3.7897574429130E+01_wp, &
      & 8.9322087944904E+00_wp, 5.3279952989770E+00_wp, 5.6258043424413E+00_wp, &
      & 1.6115591315647E+00_wp, 2.5387769376831E+00_wp, 8.2028913059178E+00_wp, &
      & 1.2226230276259E+01_wp, 5.7030141029870E+00_wp, 1.5815237509997E+00_wp, &
      & 1.4981298823014E+00_wp, 1.4317965486738E+00_wp, 3.8172825287704E+00_wp, &
      & 1.2648593650528E+01_wp, 5.8472671394681E+00_wp, 1.0631312285249E+00_wp, &
      & 7.3621267055024E+00_wp, 1.6234254583398E+00_wp, 3.4651507628618E+01_wp, &
      & 3.4923486847229E+01_wp, 1.1136757306622E+01_wp, 1.9680032604563E+01_wp, &
      & 8.3004740679372E+01_wp, 1.0042992019072E+02_wp, 3.6962950298764E+01_wp, &
      & 1.0919772339795E+01_wp, 1.0352898443990E+01_wp, 9.8943978385702E+00_wp, &
      & 2.1539554554163E+01_wp, 9.1879170429485E+01_wp, 3.7897574429130E+01_wp, &
      & 7.3621267055024E+00_wp, 9.3031746769708E+01_wp, 1.1206049984339E+01_wp, &
      & 8.1388552398630E+00_wp, 8.5961975589899E+00_wp, 2.4609116715906E+00_wp, &
      & 3.8746391243000E+00_wp, 1.2500821517587E+01_wp, 1.8653608543109E+01_wp, &
      & 8.7118494250596E+00_wp, 2.4150564158072E+00_wp, 2.2877009236808E+00_wp, &
      & 2.1864075234916E+00_wp, 5.8355457277839E+00_wp, 1.9310583681271E+01_wp, &
      & 8.9322087944904E+00_wp, 1.6234254583398E+00_wp, 1.1206049984339E+01_wp, &
      & 2.4790453596084E+00_wp], [16, 16])

   call get_structure(mol, "MB16-43", "03")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) then 
      call test_failed(error, "D4 model could not be created")
      return
   end if
   call test_c6_gen(error, mol, d4, ref, with_cn=.true., with_q=.true.)

end subroutine test_c6_d4_mb03

subroutine test_dc6_d4_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4

   call get_structure(mol, "MB16-43", "04")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) then 
      call test_failed(error, "D4 model could not be created")
      return
   end if
   call test_dc6_gen(error, mol, d4, with_cn=.true., with_q=.false.)

end subroutine test_dc6_d4_mb04

subroutine test_dc6_d4s_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4s_model) :: d4s

   call get_structure(mol, "MB16-43", "04")
   call new_d4s_model(error, d4s, mol)
   if (allocated(error)) then 
      call test_failed(error, "D4S model could not be created")
      return
   end if
   call test_dc6_gen(error, mol, d4s, with_cn=.true., with_q=.false.)

end subroutine test_dc6_d4s_mb04

subroutine test_dc6_d4_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4

   call get_structure(mol, "MB16-43", "05")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) then 
      call test_failed(error, "D4 model could not be created")
      return
   end if
   call test_dc6_gen(error, mol, d4, with_cn=.false., with_q=.true.)

end subroutine test_dc6_d4_mb05

subroutine test_dc6_d4s_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4s_model) :: d4s

   call get_structure(mol, "MB16-43", "05")
   call new_d4s_model(error, d4s, mol)
   if (allocated(error)) then 
      call test_failed(error, "D4S model could not be created")
      return
   end if
   call test_dc6_gen(error, mol, d4s, with_cn=.false., with_q=.true.)

end subroutine test_dc6_d4s_mb05

subroutine test_dc6_d4s_eeqbc_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4s_model) :: d4s

   call get_structure(mol, "MB16-43", "05")
   call new_d4s_model(error, d4s, mol, qmod=d4_qmod%eeqbc)
   if (allocated(error)) then 
      call test_failed(error, "D4S model could not be created")
      return
   end if
   call test_dc6_gen(error, mol, d4s, with_cn=.false., with_q=.true.)

end subroutine test_dc6_d4s_eeqbc_mb05

subroutine test_dc6_d4_mb06(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4

   call get_structure(mol, "MB16-43", "06")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) then 
      call test_failed(error, "D4 model could not be created")
      return
   end if
   call test_dc6_gen(error, mol, d4, with_cn=.true., with_q=.true.)

end subroutine test_dc6_d4_mb06

subroutine test_dc6_d4s_mb06(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4s_model) :: d4s

   call get_structure(mol, "MB16-43", "06")
   call new_d4s_model(error, d4s, mol)
   if (allocated(error)) then 
      call test_failed(error, "D4S model could not be created")
      return
   end if
   call test_dc6_gen(error, mol, d4s, with_cn=.true., with_q=.true.)

end subroutine test_dc6_d4s_mb06

subroutine test_c6_d4_mb07(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: qat(16) = [&
      &-1.57324183192355E-1_wp, 1.65228028672395E-1_wp, 3.22320366812437E-1_wp, &
      & 3.63579860576008E-2_wp, 4.85677294229959E-2_wp,-3.59193331069774E-1_wp, &
      &-1.93844259127416E-1_wp,-3.86492630088218E-1_wp, 3.10104713332147E-1_wp, &
      & 8.34803229863654E-2_wp,-3.62667644019899E-1_wp, 3.64142434058147E-1_wp, &
      & 3.34644499696670E-1_wp,-4.69889877462762E-1_wp,-1.89224201365947E-1_wp, &
      & 4.53790045287620E-1_wp]
   real(wp), parameter :: ref(16, 16) = reshape([&
      & 2.5899776752056E+01_wp, 2.1718975287903E+01_wp, 4.0896559404825E+00_wp, &
      & 1.9068415125987E+01_wp, 2.3793938067778E+01_wp, 3.0452162987514E+01_wp, &
      & 5.0890472676690E+01_wp, 3.4184801377206E+01_wp, 2.9794164793817E+01_wp, &
      & 5.0038725736932E+01_wp, 3.1096288315400E+01_wp, 3.7613777716825E+00_wp, &
      & 3.9882410743378E+00_wp, 4.8634793174522E+01_wp, 2.4997216709159E+01_wp, &
      & 3.1670568222466E+00_wp, 2.1718975287903E+01_wp, 1.8511876415478E+01_wp, &
      & 3.4348560517258E+00_wp, 1.6416126383445E+01_wp, 1.9955717709287E+01_wp, &
      & 2.5580529446023E+01_wp, 4.1460242401148E+01_wp, 2.8711434912189E+01_wp, &
      & 2.4534366589169E+01_wp, 4.0232199445588E+01_wp, 2.6113343013771E+01_wp, &
      & 3.1591367921103E+00_wp, 3.3496939128841E+00_wp, 4.0854373962435E+01_wp, &
      & 2.1287944173968E+01_wp, 2.6604014464326E+00_wp, 4.0896559404825E+00_wp, &
      & 3.4348560517258E+00_wp, 6.4613960289572E-01_wp, 3.0193840734568E+00_wp, &
      & 3.7572241617445E+00_wp, 4.8113918184410E+00_wp, 7.9897749151934E+00_wp, &
      & 5.4009811476696E+00_wp, 4.6856517803282E+00_wp, 7.8256165209907E+00_wp, &
      & 4.9128723609097E+00_wp, 5.9427366784289E-01_wp, 6.3011722860587E-01_wp, &
      & 7.6842179141693E+00_wp, 3.9552290087112E+00_wp, 5.0038977112109E-01_wp, &
      & 1.9068415125987E+01_wp, 1.6416126383445E+01_wp, 3.0193840734568E+00_wp, &
      & 1.4653130435190E+01_wp, 1.7521926999029E+01_wp, 2.2488882972412E+01_wp, &
      & 3.5652597710167E+01_wp, 2.5238571136210E+01_wp, 2.1242519514252E+01_wp, &
      & 3.4196881867623E+01_wp, 2.2952234544006E+01_wp, 2.7770138595030E+00_wp, &
      & 2.9445321791501E+00_wp, 3.5916746334507E+01_wp, 1.8876211586780E+01_wp, &
      & 2.3388667861345E+00_wp, 2.3793938067778E+01_wp, 1.9955717709287E+01_wp, &
      & 3.7572241617445E+00_wp, 1.7521926999029E+01_wp, 2.1859348620277E+01_wp, &
      & 2.7976875373198E+01_wp, 4.6738350872009E+01_wp, 3.1406055614332E+01_wp, &
      & 2.7366044676824E+01_wp, 4.5948203012455E+01_wp, 2.8568548657463E+01_wp, &
      & 3.4556303884227E+00_wp, 3.6640530770944E+00_wp, 4.4681540401234E+01_wp, &
      & 2.2967987820257E+01_wp, 2.9096243152442E+00_wp, 3.0452162987514E+01_wp, &
      & 2.5580529446023E+01_wp, 4.8113918184410E+00_wp, 2.2488882972412E+01_wp, &
      & 2.7976875373198E+01_wp, 3.5827458959970E+01_wp, 5.9469834950904E+01_wp, &
      & 4.0217681198927E+01_wp, 3.4880771776963E+01_wp, 5.8234629326369E+01_wp, &
      & 3.6582969286005E+01_wp, 4.4251790562564E+00_wp, 4.6920836370245E+00_wp, &
      & 5.7219618191610E+01_wp, 2.9456279627371E+01_wp, 3.7260931204162E+00_wp, &
      & 5.0890472676690E+01_wp, 4.1460242401148E+01_wp, 7.9897749151934E+00_wp, &
      & 3.5652597710167E+01_wp, 4.6738350872009E+01_wp, 5.9469834950904E+01_wp, &
      & 1.0713531598689E+02_wp, 6.6785267065801E+01_wp, 6.1391355088361E+01_wp, &
      & 1.0946639476652E+02_wp, 6.0774463131515E+01_wp, 7.3484447367073E+00_wp, &
      & 7.7915601482861E+00_wp, 9.4978522631689E+01_wp, 4.7586719998488E+01_wp, &
      & 6.1849272739308E+00_wp, 3.4184801377206E+01_wp, 2.8711434912189E+01_wp, &
      & 5.4009811476696E+00_wp, 2.5238571136210E+01_wp, 3.1406055614332E+01_wp, &
      & 4.0217681198927E+01_wp, 6.6785267065801E+01_wp, 4.5145967259725E+01_wp, &
      & 3.9166628211263E+01_wp, 6.5413078285301E+01_wp, 4.1065941214127E+01_wp, &
      & 4.9674418069468E+00_wp, 5.2670525958601E+00_wp, 6.4231190889839E+01_wp, &
      & 3.3061153040629E+01_wp, 4.1826808228299E+00_wp, 2.9794164793817E+01_wp, &
      & 2.4534366589169E+01_wp, 4.6856517803282E+00_wp, 2.1242519514252E+01_wp, &
      & 2.7366044676824E+01_wp, 3.4880771776963E+01_wp, 6.1391355088361E+01_wp, &
      & 3.9166628211263E+01_wp, 3.5483366231146E+01_wp, 6.2152639920027E+01_wp, &
      & 3.5637244661212E+01_wp, 4.3095375883142E+00_wp, 4.5694234361968E+00_wp, &
      & 5.5707650512993E+01_wp, 2.8168797141138E+01_wp, 3.6276376496133E+00_wp, &
      & 5.0038725736932E+01_wp, 4.0232199445588E+01_wp, 7.8256165209907E+00_wp, &
      & 3.4196881867623E+01_wp, 4.5948203012455E+01_wp, 5.8234629326369E+01_wp, &
      & 1.0946639476652E+02_wp, 6.5413078285301E+01_wp, 6.2152639920027E+01_wp, &
      & 1.1451912162654E+02_wp, 5.9539056863648E+01_wp, 7.1974700323621E+00_wp, &
      & 7.6314252828682E+00_wp, 9.3005759723011E+01_wp, 4.6024084470642E+01_wp, &
      & 6.0564658155562E+00_wp, 3.1096288315400E+01_wp, 2.6113343013771E+01_wp, &
      & 4.9128723609097E+00_wp, 2.2952234544006E+01_wp, 2.8568548657463E+01_wp, &
      & 3.6582969286005E+01_wp, 6.0774463131515E+01_wp, 4.1065941214127E+01_wp, &
      & 3.5637244661212E+01_wp, 5.9539056863648E+01_wp, 3.7354718892385E+01_wp, &
      & 4.5185137849098E+00_wp, 4.7910472039631E+00_wp, 5.8426234588868E+01_wp, &
      & 3.0069099854223E+01_wp, 3.8046669972155E+00_wp, 3.7613777716825E+00_wp, &
      & 3.1591367921103E+00_wp, 5.9427366784289E-01_wp, 2.7770138595030E+00_wp, &
      & 3.4556303884227E+00_wp, 4.4251790562564E+00_wp, 7.3484447367073E+00_wp, &
      & 4.9674418069468E+00_wp, 4.3095375883142E+00_wp, 7.1974700323621E+00_wp, &
      & 4.5185137849098E+00_wp, 5.4657103623079E-01_wp, 5.7953741705933E-01_wp, &
      & 7.0674020034531E+00_wp, 3.6377387176181E+00_wp, 4.6022324073896E-01_wp, &
      & 3.9882410743378E+00_wp, 3.3496939128841E+00_wp, 6.3011722860587E-01_wp, &
      & 2.9445321791501E+00_wp, 3.6640530770944E+00_wp, 4.6920836370245E+00_wp, &
      & 7.7915601482861E+00_wp, 5.2670525958601E+00_wp, 4.5694234361968E+00_wp, &
      & 7.6314252828682E+00_wp, 4.7910472039631E+00_wp, 5.7953741705933E-01_wp, &
      & 6.1449216340646E-01_wp, 7.4936722059056E+00_wp, 3.8571663116461E+00_wp, &
      & 4.8798159569671E-01_wp, 4.8634793174522E+01_wp, 4.0854373962435E+01_wp, &
      & 7.6842179141693E+00_wp, 3.5916746334507E+01_wp, 4.4681540401234E+01_wp, &
      & 5.7219618191610E+01_wp, 9.4978522631689E+01_wp, 6.4231190889839E+01_wp, &
      & 5.5707650512993E+01_wp, 9.3005759723011E+01_wp, 5.8426234588868E+01_wp, &
      & 7.0674020034531E+00_wp, 7.4936722059056E+00_wp, 9.1384787005601E+01_wp, &
      & 4.7044291543750E+01_wp, 5.9509000047757E+00_wp, 2.4997216709159E+01_wp, &
      & 2.1287944173968E+01_wp, 3.9552290087112E+00_wp, 1.8876211586780E+01_wp, &
      & 2.2967987820257E+01_wp, 2.9456279627371E+01_wp, 4.7586719998488E+01_wp, &
      & 3.3061153040629E+01_wp, 2.8168797141138E+01_wp, 4.6024084470642E+01_wp, &
      & 3.0069099854223E+01_wp, 3.6377387176181E+00_wp, 3.8571663116461E+00_wp, &
      & 4.7044291543750E+01_wp, 2.4500571380658E+01_wp, 3.0634834288830E+00_wp, &
      & 3.1670568222466E+00_wp, 2.6604014464326E+00_wp, 5.0038977112109E-01_wp, &
      & 2.3388667861345E+00_wp, 2.9096243152442E+00_wp, 3.7260931204162E+00_wp, &
      & 6.1849272739308E+00_wp, 4.1826808228299E+00_wp, 3.6276376496133E+00_wp, &
      & 6.0564658155562E+00_wp, 3.8046669972155E+00_wp, 4.6022324073896E-01_wp, &
      & 4.8798159569671E-01_wp, 5.9509000047757E+00_wp, 3.0634834288830E+00_wp, &
      & 3.8751757297523E-01_wp], [16, 16])

   type(structure_type) :: mol
   type(d4_model) :: d4

   call get_structure(mol, "MB16-43", "06")
   call new_d4_model(error, d4, mol, qmod=d4_qmod%gfn2)
   if (allocated(error)) then 
      call test_failed(error, "D4 model could not be created")
      return
   end if
   call test_c6_gen(error, mol, d4, ref, with_cn=.true., with_q=.true., qat=qat)

end subroutine test_c6_d4_mb07

subroutine test_dc6_d4_mb08(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: qat(16) = [&
      &-2.05667914412001E-1_wp,-3.99554326663093E-1_wp, 3.29243862111419E-1_wp, &
      &-3.11738256025803E-1_wp, 3.58849862618190E-2_wp, 3.21889825709581E-1_wp, &
      & 4.14746199807777E-2_wp, 2.95730046338104E-2_wp,-5.06348926564523E-1_wp, &
      & 3.43067357217139E-1_wp, 6.88375720293390E-1_wp, 7.03350634832100E-2_wp, &
      &-9.62426987152087E-2_wp,-1.32210876939567E-1_wp, 9.79003738112971E-2_wp, &
      &-3.05981814182260E-1_wp]

   type(structure_type) :: mol
   type(d4_model) :: d4

   call get_structure(mol, "MB16-43", "08")
   call new_d4_model(error, d4, mol, qmod=d4_qmod%gfn2)
   if (allocated(error)) then 
      call test_failed(error, "D4 model could not be created")
      return
   end if
   call test_dc6_gen(error, mol, d4, with_cn=.true., with_q=.true., qat=qat)

end subroutine test_dc6_d4_mb08

subroutine test_dc6_d4s_mb08(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: qat(16) = [&
      &-2.05667914412001E-1_wp,-3.99554326663093E-1_wp, 3.29243862111419E-1_wp, &
      &-3.11738256025803E-1_wp, 3.58849862618190E-2_wp, 3.21889825709581E-1_wp, &
      & 4.14746199807777E-2_wp, 2.95730046338104E-2_wp,-5.06348926564523E-1_wp, &
      & 3.43067357217139E-1_wp, 6.88375720293390E-1_wp, 7.03350634832100E-2_wp, &
      &-9.62426987152087E-2_wp,-1.32210876939567E-1_wp, 9.79003738112971E-2_wp, &
      &-3.05981814182260E-1_wp]

   type(structure_type) :: mol
   type(d4s_model) :: d4s

   call get_structure(mol, "MB16-43", "08")
   call new_d4s_model(error, d4s, mol, qmod=d4_qmod%gfn2)
   if (allocated(error)) then 
      call test_failed(error, "D4S model could not be created")
      return
   end if
   call test_dc6_gen(error, mol, d4s, with_cn=.true., with_q=.true., qat=qat)

end subroutine test_dc6_d4s_mb08


subroutine test_pol_d4_mb09(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   real(wp), parameter :: ref(16) = [&
      & 2.7207807469636E+00_wp, 2.7207697669106E+00_wp, 2.7340295240471E+00_wp, &
      & 2.7373201241818E+00_wp, 2.0051462957754E+01_wp, 2.7348043044018E+00_wp, &
      & 6.9923962318719E+00_wp, 9.2345547677329E+00_wp, 2.7338175995796E+00_wp, &
      & 2.7347959316428E+00_wp, 2.3018401095949E+01_wp, 2.7347926176379E+00_wp, &
      & 1.5079046616252E+01_wp, 3.5315073387745E+00_wp, 2.7340548152913E+00_wp, &
      & 9.2442543133947E+00_wp]
   real(wp), parameter :: ref_qq(16) = [&
      & 5.4615565063498E+00_wp, 5.4615344655514E+00_wp, 5.4881514257539E+00_wp, &
      & 5.4947568086372E+00_wp, 1.0065572291228E+02_wp, 5.4897066803227E+00_wp, &
      & 2.1710888390047E+01_wp, 3.3655020315897E+01_wp, 5.4877260193864E+00_wp, &
      & 5.4896898732732E+00_wp, 1.1239537836978E+02_wp, 5.4896832209087E+00_wp, &
      & 5.6234642849808E+01_wp, 8.4341312278194E+00_wp, 5.4882021941076E+00_wp, &
      & 3.3690369979688E+01_wp]

   call get_structure(mol, "MB16-43", "09")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) then 
      call test_failed(error, "D4 model could not be created")
      return
   end if
   call test_pol_gen(error, mol, d4, ref, ref_qq, &
      & with_cn=.true., with_q=.false.)

end subroutine test_pol_d4_mb09

subroutine test_pol_d4s_mb09(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4s_model) :: d4s
   real(wp), parameter :: ref(16) = [&
      & 2.0796638857659E+00_wp, 2.0317560362159E+00_wp, 2.1362479991008E+00_wp, &
      & 2.4840420351444E+00_wp, 3.5144617579624E+01_wp, 2.1798399760895E+00_wp, &
      & 8.7285804463012E+00_wp, 1.0754720091309E+01_wp, 2.1068749047750E+00_wp, &
      & 2.2167964014782E+00_wp, 2.4223007635619E+01_wp, 2.1626725547918E+00_wp, &
      & 1.5342832501948E+01_wp, 4.0831999012216E+00_wp, 2.2319948533110E+00_wp, &
      & 1.0723200212288E+01_wp]
   real(wp), parameter :: ref_qq(16) = [&
      & 4.1746112173873E+00_wp, 4.0784434435940E+00_wp, 4.2881952805960E+00_wp, &
      & 4.9863392903782E+00_wp, 1.7642138613056E+02_wp, 4.3756995919275E+00_wp, &
      & 2.7101615753611E+01_wp, 3.9195211059822E+01_wp, 4.2292332291312E+00_wp, &
      & 4.4498840354034E+00_wp, 1.1827729029096E+02_wp, 4.3412385859865E+00_wp, &
      & 5.7218385751360E+01_wp, 9.7517123688813E+00_wp, 4.4803926324619E+00_wp, &
      & 3.9080337934321E+01_wp]

   call get_structure(mol, "MB16-43", "09")
   call new_d4s_model(error, d4s, mol)
   if (allocated(error)) then 
      call test_failed(error, "D4S model could not be created")
      return
   end if
   call test_pol_gen(error, mol, d4s, ref, ref_qq, &
      & with_cn=.true., with_q=.true.)

end subroutine test_pol_d4s_mb09

subroutine test_dpol_d4_mb10(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4

   call get_structure(mol, "MB16-43", "10")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) then 
      call test_failed(error, "D4 model could not be created")
      return
   end if
   call test_dpol_gen(error, mol, d4, with_cn=.true., with_q=.true.)

end subroutine test_dpol_d4_mb10

subroutine test_dpol_d4s_mb10(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4s_model) :: d4s

   call get_structure(mol, "MB16-43", "10")
   call new_d4s_model(error, d4s, mol)
   if (allocated(error)) then 
      call test_failed(error, "D4S model could not be created")
      return
   end if
   call test_dpol_gen(error, mol, d4s, with_cn=.true., with_q=.true.)

end subroutine test_dpol_d4s_mb10


subroutine test_d4_model_error(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4

   integer, parameter :: nat = 2
   character(len=*), parameter :: sym(nat) = [character(len=4) :: "Rf", "Db"]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp],&
      & [3, nat])

   call new(mol, sym, xyz)
   call new_d4_model(error, d4, mol)
  
end subroutine test_d4_model_error

subroutine test_d4s_model_error(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4s_model) :: d4s

   integer, parameter :: nat = 2
   character(len=*), parameter :: sym(nat) = [character(len=4) :: "H", "Db"]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp],&
      & [3, nat])

   call new(mol, sym, xyz)
   call new_d4s_model(error, d4s, mol)
  
end subroutine test_d4s_model_error

subroutine test_model_wrapper(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(dispersion_model), allocatable :: disp

   integer, parameter :: nat = 2
   character(len=*), parameter :: sym(nat) = [character(len=4) :: "H", "H"]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp],&
      & [3, nat])

   call new(mol, sym, xyz)

   call new_dispersion_model(error, disp, mol, "d4")
   if (allocated(error)) then
      call test_failed(error, "D4 model could not be created")
      return
   end if

   ! check if the model is of type d4_model
   select type (disp)
      type is (d4_model)
      class default
         call test_failed(error, "Expected d4_model but got something else")
         return
   end select

   deallocate(disp)
   call new_dispersion_model(error, disp, mol, "D4S")
   if (allocated(error)) then
      call test_failed(error, "D4S model could not be created")
      return
   end if

   ! check if the model is of type d4s_model
   select type (disp)
      type is (d4s_model)
      class default
         call test_failed(error, "Expected d4s_model but got something else")
         return
   end select

end subroutine test_model_wrapper

subroutine test_model_wrapper_fail(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(dispersion_model), allocatable :: disp

   integer, parameter :: nat = 2
   character(len=*), parameter :: sym(nat) = [character(len=4) :: "H", "H"]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp],&
      & [3, nat])

   call new(mol, sym, xyz)

   call new_dispersion_model(error, disp, mol, "wrong")

   if (.not. allocated(error)) then
      call test_failed(error, "Expected an error for key 'wrong'")
      return
   end if

   if (allocated(disp)) then
      call test_failed(error, "Model should not be allocated with invalid key")
      return
   end if

end subroutine test_model_wrapper_fail


subroutine test_id_from_name(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: id

   call get_dispersion_model_id(error, "d4", id)
   if (allocated(error)) return
   call check(error, id == dftd_models%d4, "d4 name ID mismatch")
   if (allocated(error)) return

   ! Test case-insensitivity
   call get_dispersion_model_id(error, "D4", id)
   if (allocated(error)) return
   call check(error, id == dftd_models%d4, "D4 case-insensitive ID mismatch")
   if (allocated(error)) return

   call get_dispersion_model_id(error, "d4s", id)
   if (allocated(error)) return
   call check(error, id == dftd_models%d4s, "d4s name ID mismatch")
   if (allocated(error)) return

   call get_dispersion_model_id(error, "d4srev", id)
   if (allocated(error)) return
   call check(error, id == dftd_models%d4srev, "d4srev name ID mismatch")

end subroutine test_id_from_name


subroutine test_id_from_name_empty(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: id

   call get_dispersion_model_id(error, "", id)

end subroutine test_id_from_name_empty


subroutine test_id_from_name_wrong(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: id

   call get_dispersion_model_id(error, "d5", id)

end subroutine test_id_from_name_wrong


subroutine test_id_from_model(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(d4s_model) :: d4s
   type(d4srev_model) :: d4srev
   integer :: id

   integer, parameter :: nat = 2
   character(len=*), parameter :: sym(nat) = [character(len=4) :: "H", "H"]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.4_wp], [3, nat])

   call new(mol, sym, xyz)

   call new_d4_model(error, d4, mol)
   if (allocated(error)) return

   call get_dispersion_model_id(error, d4, id)
   if (allocated(error)) return
   call check(error, id == dftd_models%d4, "d4 model ID mismatch")
   if (allocated(error)) return

   call new_d4s_model(error, d4s, mol)
   if (allocated(error)) return

   call get_dispersion_model_id(error, d4s, id)
   if (allocated(error)) return
   call check(error, id == dftd_models%d4s, "d4s model ID mismatch")
   if (allocated(error)) return

   call new_d4srev_model(error, d4srev, mol)
   if (allocated(error)) return

   call get_dispersion_model_id(error, d4srev, id)
   if (allocated(error)) return
   call check(error, id == dftd_models%d4srev, "d4srev model ID mismatch")

end subroutine test_id_from_model


end module test_model
