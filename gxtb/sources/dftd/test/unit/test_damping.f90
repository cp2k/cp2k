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

module test_damping
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type, new
   use mstore, only : get_structure
   use dftd4
   implicit none
   private

   public :: collect_damping

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))


contains


!> Collect all exported unit tests
subroutine collect_damping(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("damp-rational-2b-m01", test_damp_rational_2b_mb01), &
      & new_unittest("damp-rational-2b-m02", test_damp_rational_2b_mb02), &
      & new_unittest("grad-rational-2b-m03", test_grad_rational_2b_mb03), &
      & new_unittest("grad-rational-2b-m04", test_grad_rational_2b_mb04), &
      & new_unittest("damp-screened-2b-m01", test_damp_screened_2b_mb01), &
      & new_unittest("damp-screened-2b-m02", test_damp_screened_2b_mb02), &
      & new_unittest("grad-screened-2b-m03", test_grad_screened_2b_mb03), &
      & new_unittest("grad-screened-2b-m04", test_grad_screened_2b_mb04), &
      & new_unittest("damp-zero-2b-m01", test_damp_zero_2b_mb01), &
      & new_unittest("damp-zero-2b-m02", test_damp_zero_2b_mb02), &
      & new_unittest("grad-zero-2b-m03", test_grad_zero_2b_mb03), &
      & new_unittest("grad-zero-2b-m04", test_grad_zero_2b_mb04), &
      & new_unittest("damp-mzero-2b-m01", test_damp_mzero_2b_mb01), &
      & new_unittest("damp-mzero-2b-m02", test_damp_mzero_2b_mb02), &
      & new_unittest("grad-mzero-2b-m03", test_grad_mzero_2b_mb03), &
      & new_unittest("grad-mzero-2b-m04", test_grad_mzero_2b_mb04), &
      & new_unittest("damp-optpower-2b-m01", test_damp_optpower_2b_mb01), &
      & new_unittest("damp-optpower-2b-m02", test_damp_optpower_2b_mb02), &
      & new_unittest("grad-optpower-2b-m03", test_grad_optpower_2b_mb03), &
      & new_unittest("grad-optpower-2b-m04", test_grad_optpower_2b_mb04), &
      & new_unittest("damp-cso_2b_m01", test_damp_cso_2b_mb01), &
      & new_unittest("damp-cso_2b_m02", test_damp_cso_2b_mb02), &
      & new_unittest("grad-cso_2b_m03", test_grad_cso_2b_mb03), &
      & new_unittest("grad-cso_2b_m04", test_grad_cso_2b_mb04), &
      & new_unittest("damp-koide-2b-m01", test_damp_koide_2b_mb01), &
      & new_unittest("damp-koide-2b-m02", test_damp_koide_2b_mb02), &
      & new_unittest("grad-koide-2b-m03", test_grad_koide_2b_mb03), &
      & new_unittest("grad-koide-2b-m04", test_grad_koide_2b_mb04), &
      & new_unittest("damp-rational-3b-m01", test_damp_rational_3b_mb01), &
      & new_unittest("damp-rational-3b-m02", test_damp_rational_3b_mb02), &
      & new_unittest("grad-rational-3b-m03", test_grad_rational_3b_mb03), &
      & new_unittest("grad-rational-3b-m04", test_grad_rational_3b_mb04), &
      & new_unittest("damp-screened-3b-m01", test_damp_screened_3b_mb01), &
      & new_unittest("damp-screened-3b-m02", test_damp_screened_3b_mb02), &
      & new_unittest("grad-screened-3b-m03", test_grad_screened_3b_mb03), &
      & new_unittest("grad-screened-3b-m04", test_grad_screened_3b_mb04), &
      & new_unittest("damp-zero-3b-m01", test_damp_zero_3b_mb01), &
      & new_unittest("damp-zero-3b-m02", test_damp_zero_3b_mb02), &
      & new_unittest("grad-zero-3b-m03", test_grad_zero_3b_mb03), &
      & new_unittest("grad-zero-3b-m04", test_grad_zero_3b_mb04), &
      & new_unittest("damp-zero-avg-3b-m01", test_damp_zero_avg_3b_mb01), &
      & new_unittest("damp-zero-avg-3b-m02", test_damp_zero_avg_3b_mb02), &
      & new_unittest("grad-zero-avg-3b-m03", test_grad_zero_avg_3b_mb03), &
      & new_unittest("grad-zero-avg-3b-m04", test_grad_zero_avg_3b_mb04), &
      & new_unittest("damping-empty", test_damping_empty, should_fail=.true.), &
      & new_unittest("damping-3b-empty", test_damping_3b_empty), &
      & new_unittest("params-rational", test_params_rational), &
      & new_unittest("params-rational-empty", test_params_rational_empty, should_fail=.true.), &
      & new_unittest("params-screened", test_params_screened), &
      & new_unittest("params-screened-empty", test_params_screened_empty, should_fail=.true.), &
      & new_unittest("params-zero", test_params_zero), &
      & new_unittest("params-zero-empty", test_params_zero_empty, should_fail=.true.), &
      & new_unittest("params-mzero-zero-avg", test_params_mzero_zero_avg), &
      & new_unittest("params-mzero-zero-avg-empty", test_params_mzero_zero_avg_empty, should_fail=.true.), &
      & new_unittest("params-optpower", test_params_optpower), &
      & new_unittest("params-optpower-empty", test_params_optpower_empty, should_fail=.true.), &
      & new_unittest("params-cso", test_params_cso), &
      & new_unittest("params-cso-empty", test_params_cso_empty, should_fail=.true.), &
      & new_unittest("params-cso-noa1", test_params_cso_noa1, should_fail=.true.), &
      & new_unittest("params-koide", test_params_koide), &
      & new_unittest("params-koide-empty", test_params_koide_empty, should_fail=.true.), &
      & new_unittest("id-from-name", test_id_from_name), &
      & new_unittest("id-from-name-empty", test_id_from_name_empty, should_fail=.true.), &
      & new_unittest("id-from-name-wrong", test_id_from_name_wrong, should_fail=.true.), &
      & new_unittest("id-from-damping", test_id_from_damping) &
      & ]

end subroutine collect_damping


subroutine test_damping_2b_gen(error, mol, d4, damp, param, ref6, ref8)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Dispersion model
   class(dispersion_model), intent(in) :: d4

   !> Damping function
   type(damping_type), intent(in) :: damp
   
   !> Damping parameters
   type(param_type), intent(in) :: param

   !> Expected damping function for C6 term
   real(wp), intent(in) :: ref6(:, :)

   !> Expected damping function for C8 term
   real(wp), intent(in), optional :: ref8(:, :)

   real(wp), allocatable :: fdmp6(:, :), fdmp8(:, :)
   integer :: iat, jat, jtr, izp, jzp
   real(wp) :: vec(3), r2, rdamp, d6, d8

   allocate(fdmp6(mol%nat, mol%nat), fdmp8(mol%nat, mol%nat), source=0.0_wp)

   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, mol%nat
         jzp = mol%id(jat)
         if (jat == iat) cycle
         vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat)
         r2 = vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3)
         ! Get damping radius
         rdamp = d4%get_2b_rdamp(izp, jzp)
         ! Get damping function value
         call damp%damping_2b%get_2b_damp(param, r2, rdamp, d6, d8)

         fdmp6(iat, jat) = d6
         fdmp8(iat, jat) = d8
      end do
   end do

   if (any(abs(fdmp6 - ref6) > thr)) then
      call test_failed(error, "Damping function for C6 does not match")
      write(*,*) 'fdmp6'
      print'(3es21.14)', fdmp6
      write(*,*) 'ref'
      print'(3es21.14)', ref6
      write(*,*) "diff:"
      print'(3es21.14)', abs(fdmp6 - ref6)
   end if

   if (present(ref8)) then
      if (any(abs(fdmp8 - ref8) > thr)) then
         call test_failed(error, "Damping function for C8 does not match")
         write(*,*) 'fdmp8'
         print'(3es21.14)', fdmp8
         write(*,*) 'ref'
         print'(3es21.14)', ref8
         write(*,*) "diff:"
         print'(3es21.14)', abs(fdmp8 - ref8)
      end if
   end if

   ! write(*,*) "fdmp6:"
   ! print '(*(6x,"&", 3(es21.14e2, "_wp":, ","), "&", /))', fdmp6

   ! write(*,*) "fdmp8:"
   ! print '(*(6x,"&", 3(es21.14e2, "_wp":, ","), "&", /))', fdmp8

end subroutine test_damping_2b_gen

subroutine test_damping_2b_numgrad(error, mol, d4, damp, param)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Dispersion model
   class(dispersion_model), intent(in) :: d4

   !> Damping function
   type(damping_type), intent(in) :: damp
   
   !> Damping parameters
   type(param_type), intent(in) :: param

   integer :: iat, jat, izp, jzp
   real(wp) :: vec(3), r2, r, rdamp
   real(wp) :: d6, d8, d6dr, d8dr
   real(wp) :: d6_p, d8_p, d6_m, d8_m, r_p, r2_p, r_m, r2_m
   real(wp) :: num_d6dr, num_d8dr
   real(wp), parameter :: step = 1.0e-5_wp

   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, mol%nat
         if (jat == iat) cycle
         jzp = mol%id(jat)

         vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat)
         r2 = vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3)
         r = sqrt(r2)

         ! Get damping radius
         rdamp = d4%get_2b_rdamp(izp, jzp)

         ! Get analytical damping derivatives
         call damp%damping_2b%get_2b_derivs(param, r2, rdamp, d6, d8, d6dr, d8dr)

         ! Forward step
         r_p = r + step
         r2_p = r_p * r_p
         call damp%damping_2b%get_2b_damp(param, r2_p, rdamp, d6_p, d8_p)

         ! Backward step
         r_m = r - step
         r2_m = r_m * r_m
         call damp%damping_2b%get_2b_damp(param, r2_m, rdamp, d6_m, d8_m)

         ! Calculate the derivative devided by r
         num_d6dr = (d6_p - d6_m) / (2.0_wp * step * r)
         num_d8dr = (d8_p - d8_m) / (2.0_wp * step * r)

         if (abs(d6dr - num_d6dr) > thr) then
            call test_failed(error, "Gradient of C6 damping does not match.")
            print'(3es21.14)', d6dr, num_d6dr, abs(d6dr - num_d6dr)
            return
         end if

         if (abs(d8dr - num_d8dr) > thr) then
            call test_failed(error, "Gradient of C8 damping does not match.")
            print'(3es21.14)', d8dr, num_d8dr, abs(d8dr - num_d8dr)
            return
         end if
      end do
   end do

end subroutine test_damping_2b_numgrad


subroutine test_damping_3b_gen(error, mol, d4, damp, param, ref9)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Dispersion model
   class(dispersion_model), intent(in) :: d4

   !> Damping function
   type(damping_type), intent(in) :: damp
   
   !> Damping parameters
   type(param_type), intent(in) :: param

   !> Expected damping function for C9 term (slice for reference atom 1)
   real(wp), intent(in) :: ref9(:, :)

   integer :: iat, jat, kat, izp, jzp, kzp
   real(wp) :: vij(3), vjk(3), vik(3), r2ij, r2jk, r2ik, r1, r2
   real(wp) :: rdamp, rdampij, rdampik, rdampjk, d9
   
   real(wp), allocatable :: fdmp9(:, :)

   allocate(fdmp9(mol%nat, mol%nat), source=0.0_wp)

   ! Fix the first atom to compute a 2D slice
   iat = 1
   izp = mol%id(iat)

   do jat = 1, mol%nat
      if (jat == iat) cycle
      jzp = mol%id(jat)
      
      rdampij = d4%get_3b_rdamp(izp, jzp)
      
      vij(:) = mol%xyz(:, jat) - mol%xyz(:, iat)
      r2ij = vij(1)*vij(1) + vij(2)*vij(2) + vij(3)*vij(3)

      do kat = 1, mol%nat
         if (kat == iat .or. kat == jat) cycle
         kzp = mol%id(kat)
         
         rdampik = d4%get_3b_rdamp(izp, kzp)
         rdampjk = d4%get_3b_rdamp(jzp, kzp)
         rdamp = rdampij * rdampik * rdampjk
         
         vik(:) = mol%xyz(:, kat) - mol%xyz(:, iat)
         r2ik = vik(1)*vik(1) + vik(2)*vik(2) + vik(3)*vik(3)
         
         vjk(:) = mol%xyz(:, kat) - mol%xyz(:, jat)
         r2jk = vjk(1)*vjk(1) + vjk(2)*vjk(2) + vjk(3)*vjk(3)
         
         r2 = r2ij * r2ik * r2jk
         r1 = sqrt(r2)

         call damp%damping_3b%get_3b_damp(param, r1, r2ij, r2ik, r2jk, &
               & rdamp, rdampij, rdampik, rdampjk, d9)
         
         fdmp9(jat, kat) = d9
      end do
   end do

   ! Check Damping values against reference
   if (any(abs(fdmp9 - ref9) > thr)) then
      call test_failed(error, "Damping function for C9 does not match reference.")
      write(*,*) 'fdmp9 (slice for reference atom 1)'
      print'(3es21.14)', fdmp9
      write(*,*) 'ref9'
      print'(3es21.14)', ref9
   end if

   ! write(*,*) "fdmp9:"
   ! print '(*(6x,"&", 3(es21.14e2, "_wp":, ","), "&", /))', fdmp9

end subroutine test_damping_3b_gen


subroutine test_damping_3b_numgrad(error, mol, d4, damp, param)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Dispersion model
   class(dispersion_model), intent(in) :: d4

   !> Damping function
   type(damping_type), intent(in) :: damp
   
   !> Damping parameters
   type(param_type), intent(in) :: param

   integer :: iat, jat, kat, izp, jzp, kzp
   real(wp) :: vij(3), vjk(3), vik(3), r2ij, r2jk, r2ik, r1, r2
   real(wp) :: rdamp, rdampij, rdampik, rdampjk
   real(wp) :: d9, d9drij, d9drik, d9drjk
   real(wp) :: d9_p, d9_m, r_p, r2ij_p, r_m, r2ij_m, num_d9drij
   real(wp), parameter :: step = 1.0e-5_wp

   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, mol%nat
         if (jat == iat) cycle
         jzp = mol%id(jat)
         
         rdampij = d4%get_3b_rdamp(izp, jzp)
         
         vij(:) = mol%xyz(:, jat) - mol%xyz(:, iat)
         r2ij = vij(1)*vij(1) + vij(2)*vij(2) + vij(3)*vij(3)

         do kat = 1, mol%nat
            if (kat == iat .or. kat == jat) cycle
            kzp = mol%id(kat)
            
            rdampik = d4%get_3b_rdamp(izp, kzp)
            rdampjk = d4%get_3b_rdamp(jzp, kzp)
            rdamp = rdampij * rdampik * rdampjk
            
            vik(:) = mol%xyz(:, kat) - mol%xyz(:, iat)
            r2ik = vik(1)*vik(1) + vik(2)*vik(2) + vik(3)*vik(3)
            
            vjk(:) = mol%xyz(:, kat) - mol%xyz(:, jat)
            r2jk = vjk(1)*vjk(1) + vjk(2)*vjk(2) + vjk(3)*vjk(3)
            
            r2 = r2ij * r2ik * r2jk
            r1 = sqrt(r2)

            ! Calculate analytical derivatives
            call damp%damping_3b%get_3b_derivs(param, r1, r2ij, r2ik, r2jk, rdamp, &
                  & rdampij, rdampik, rdampjk, d9, d9drij, d9drik, d9drjk)
            
            ! Forward step along r_ij
            r_p = sqrt(r2ij) + step
            r2ij_p = r_p * r_p
            r1 = sqrt(r2ij_p * r2ik * r2jk)
            call damp%damping_3b%get_3b_damp(param, r1, r2ij_p, r2ik, r2jk, &
                  & rdamp, rdampij, rdampik, rdampjk, d9_p)

            ! Backward step along r_ij
            r_m = sqrt(r2ij) - step
            r2ij_m = r_m * r_m
            r1 = sqrt(r2ij_m * r2ik * r2jk)
            call damp%damping_3b%get_3b_damp(param, r1, r2ij_m, r2ik, r2jk, &
                  & rdamp, rdampij, rdampik, rdampjk, d9_m)

            ! Evaluate numerical gradient derived by r_ij
            num_d9drij = (d9_p - d9_m) / (2.0_wp * step * sqrt(r2ij))

            if (abs(d9drij - num_d9drij) > thr2) then
               call test_failed(error, "Numerical gradient of C9 damping does not match")
               print'(3es21.14)', d9drij, num_d9drij, abs(d9drij - num_d9drij)
               return
            end if

         end do
      end do
   end do

end subroutine test_damping_3b_numgrad


subroutine test_damp_rational_2b_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   real(wp), parameter :: ref6(16, 16) = reshape([&
      & 0.00000000000000E+00_wp, 2.92771306317172E-06_wp, 5.20224505553292E-06_wp,&
      & 6.66208228183840E-06_wp, 2.79102642376119E-06_wp, 2.47928100514870E-06_wp,&
      & 6.07839019056014E-06_wp, 1.25282067093488E-06_wp, 4.79094216702847E-06_wp,&
      & 9.47285853578846E-07_wp, 1.74117672073466E-06_wp, 2.18921026846665E-06_wp,&
      & 3.22820642251887E-06_wp, 3.40454116220084E-06_wp, 2.68857834138824E-06_wp,&
      & 1.77077884309694E-06_wp, 2.92771306317172E-06_wp, 0.00000000000000E+00_wp,&
      & 6.18254144126491E-07_wp, 1.74464020773302E-05_wp, 2.09127083160798E-07_wp,&
      & 6.41063915401796E-07_wp, 1.75525036460926E-05_wp, 4.08958088304507E-06_wp,&
      & 2.11314690739555E-06_wp, 5.92353719933680E-06_wp, 4.35667233844064E-07_wp,&
      & 9.99362993741571E-06_wp, 1.19269749452743E-05_wp, 6.33360308104738E-07_wp,&
      & 1.00495276607781E-06_wp, 5.02172301884075E-06_wp, 5.20224505553292E-06_wp,&
      & 6.18254144126491E-07_wp, 0.00000000000000E+00_wp, 2.89742578663825E-06_wp,&
      & 1.21509524605811E-05_wp, 3.11306144348237E-06_wp, 4.33102082478293E-06_wp,&
      & 1.14167907809834E-06_wp, 1.10769440892605E-05_wp, 6.96789563794283E-07_wp,&
      & 4.12551856586243E-06_wp, 8.75400857570605E-07_wp, 1.81779989565535E-06_wp,&
      & 9.56606736612782E-06_wp, 4.57796398178965E-06_wp, 3.05871145927592E-06_wp,&
      & 6.66208228183840E-06_wp, 1.74464020773302E-05_wp, 2.89742578663825E-06_wp,&
      & 0.00000000000000E+00_wp, 5.96376683225238E-07_wp, 1.63081390636728E-06_wp,&
      & 1.74855011426249E-05_wp, 2.15313091371526E-06_wp, 6.38411204942689E-06_wp,&
      & 1.99254976331599E-06_wp, 8.11172473535117E-07_wp, 1.00280658649558E-05_wp,&
      & 1.19252130284801E-05_wp, 2.33422573276254E-06_wp, 2.11530157272343E-06_wp,&
      & 5.54611652433839E-06_wp, 2.79102642376119E-06_wp, 2.09127083160798E-07_wp,&
      & 1.21509524605811E-05_wp, 5.96376683225238E-07_wp, 0.00000000000000E+00_wp,&
      & 9.55241777658837E-06_wp, 1.35678146890306E-06_wp, 1.51618281511852E-06_wp,&
      & 1.18380380107070E-05_wp, 6.79671446246605E-07_wp, 1.45630005001257E-05_wp,&
      & 6.10156661278750E-07_wp, 5.72848898554305E-07_wp, 1.02783366074826E-05_wp,&
      & 1.03433947343193E-05_wp, 3.17225466253770E-06_wp, 2.47928100514870E-06_wp,&
      & 6.41063915401796E-07_wp, 3.11306144348237E-06_wp, 1.63081390636728E-06_wp,&
      & 9.55241777658837E-06_wp, 0.00000000000000E+00_wp, 2.78574381838554E-06_wp,&
      & 4.37338264875123E-06_wp, 1.41629255330471E-05_wp, 1.59917548591407E-06_wp,&
      & 1.81976175363918E-05_wp, 5.85554985190050E-06_wp, 1.88376384621664E-06_wp,&
      & 8.18054633609327E-06_wp, 1.50167759215813E-05_wp, 7.08187190862342E-06_wp,&
      & 6.07839019056014E-06_wp, 1.75525036460926E-05_wp, 4.33102082478293E-06_wp,&
      & 1.74855011426249E-05_wp, 1.35678146890306E-06_wp, 2.78574381838554E-06_wp,&
      & 0.00000000000000E+00_wp, 1.20064929143903E-05_wp, 1.12719644747794E-05_wp,&
      & 1.37536709126639E-05_wp, 2.63775076024445E-06_wp, 1.01089760862117E-05_wp,&
      & 1.19216735841076E-05_wp, 4.56658691989427E-06_wp, 5.45603770389709E-06_wp,&
      & 7.94049169668054E-06_wp, 1.25282067093488E-06_wp, 4.08958088304507E-06_wp,&
      & 1.14167907809834E-06_wp, 2.15313091371526E-06_wp, 1.51618281511852E-06_wp,&
      & 4.37338264875123E-06_wp, 1.20064929143903E-05_wp, 0.00000000000000E+00_wp,&
      & 7.50454779771964E-06_wp, 1.55152023983390E-05_wp, 8.04587745132172E-06_wp,&
      & 5.87046301707544E-06_wp, 5.11674027455015E-06_wp, 2.47871145178959E-06_wp,&
      & 8.55886331397517E-06_wp, 6.22219249067286E-06_wp, 4.79094216702847E-06_wp,&
      & 2.11314690739555E-06_wp, 1.10769440892605E-05_wp, 6.38411204942689E-06_wp,&
      & 1.18380380107070E-05_wp, 1.41629255330471E-05_wp, 1.12719644747794E-05_wp,&
      & 7.50454779771964E-06_wp, 0.00000000000000E+00_wp, 4.18351499507724E-06_wp,&
      & 1.41856357526654E-05_wp, 5.72371125289642E-06_wp, 5.44945220826229E-06_wp,&
      & 9.16506887753609E-06_wp, 1.18114199614298E-05_wp, 5.92422590841259E-06_wp,&
      & 9.47285853578846E-07_wp, 5.92353719933680E-06_wp, 6.96789563794283E-07_wp,&
      & 1.99254976331599E-06_wp, 6.79671446246605E-07_wp, 1.59917548591407E-06_wp,&
      & 1.37536709126639E-05_wp, 1.55152023983390E-05_wp, 4.18351499507724E-06_wp,&
      & 0.00000000000000E+00_wp, 3.11233699262440E-06_wp, 5.06384084573986E-06_wp,&
      & 5.80735004636926E-06_wp, 1.28585153889435E-06_wp, 4.16315490672691E-06_wp,&
      & 7.59384703209654E-06_wp, 1.74117672073466E-06_wp, 4.35667233844064E-07_wp,&
      & 4.12551856586243E-06_wp, 8.11172473535117E-07_wp, 1.45630005001257E-05_wp,&
      & 1.81976175363918E-05_wp, 2.63775076024445E-06_wp, 8.04587745132172E-06_wp,&
      & 1.41856357526654E-05_wp, 3.11233699262440E-06_wp, 0.00000000000000E+00_wp,&
      & 2.54278819289747E-06_wp, 1.16275908765845E-06_wp, 9.93229141146896E-06_wp,&
      & 1.50166132680335E-05_wp, 7.18486068808024E-06_wp, 2.18921026846665E-06_wp,&
      & 9.99362993741571E-06_wp, 8.75400857570605E-07_wp, 1.00280658649558E-05_wp,&
      & 6.10156661278750E-07_wp, 5.85554985190050E-06_wp, 1.01089760862117E-05_wp,&
      & 5.87046301707544E-06_wp, 5.72371125289642E-06_wp, 5.06384084573986E-06_wp,&
      & 2.54278819289747E-06_wp, 0.00000000000000E+00_wp, 6.63746713188362E-06_wp,&
      & 1.44331589514881E-06_wp, 5.57988036309748E-06_wp, 4.09082624773462E-06_wp,&
      & 3.22820642251887E-06_wp, 1.19269749452743E-05_wp, 1.81779989565535E-06_wp,&
      & 1.19252130284801E-05_wp, 5.72848898554305E-07_wp, 1.88376384621664E-06_wp,&
      & 1.19216735841076E-05_wp, 5.11674027455015E-06_wp, 5.44945220826229E-06_wp,&
      & 5.80735004636926E-06_wp, 1.16275908765845E-06_wp, 6.63746713188362E-06_wp,&
      & 0.00000000000000E+00_wp, 1.83647878871814E-06_wp, 2.73838745156716E-06_wp,&
      & 4.07341455541469E-06_wp, 3.40454116220084E-06_wp, 6.33360308104738E-07_wp,&
      & 9.56606736612782E-06_wp, 2.33422573276254E-06_wp, 1.02783366074826E-05_wp,&
      & 8.18054633609327E-06_wp, 4.56658691989427E-06_wp, 2.47871145178959E-06_wp,&
      & 9.16506887753609E-06_wp, 1.28585153889435E-06_wp, 9.93229141146896E-06_wp,&
      & 1.44331589514881E-06_wp, 1.83647878871814E-06_wp, 0.00000000000000E+00_wp,&
      & 8.12412187352909E-06_wp, 3.59658326858472E-06_wp, 2.68857834138824E-06_wp,&
      & 1.00495276607781E-06_wp, 4.57796398178965E-06_wp, 2.11530157272343E-06_wp,&
      & 1.03433947343193E-05_wp, 1.50167759215813E-05_wp, 5.45603770389709E-06_wp,&
      & 8.55886331397517E-06_wp, 1.18114199614298E-05_wp, 4.16315490672691E-06_wp,&
      & 1.50166132680335E-05_wp, 5.57988036309748E-06_wp, 2.73838745156716E-06_wp,&
      & 8.12412187352909E-06_wp, 0.00000000000000E+00_wp, 5.89134391594683E-06_wp,&
      & 1.77077884309694E-06_wp, 5.02172301884075E-06_wp, 3.05871145927592E-06_wp,&
      & 5.54611652433839E-06_wp, 3.17225466253770E-06_wp, 7.08187190862342E-06_wp,&
      & 7.94049169668054E-06_wp, 6.22219249067286E-06_wp, 5.92422590841259E-06_wp,&
      & 7.59384703209654E-06_wp, 7.18486068808024E-06_wp, 4.09082624773462E-06_wp,&
      & 4.07341455541469E-06_wp, 3.59658326858472E-06_wp, 5.89134391594683E-06_wp,&
      & 0.00000000000000E+00_wp], shape(ref6))

   real(wp), parameter :: ref8(16, 16) = reshape([&
      & 0.00000000000000E+00_wp, 5.03464173832655E-08_wp, 8.85630798816101E-08_wp,&
      & 1.25589211938184E-07_wp, 4.74868842272753E-08_wp, 3.98242609707400E-08_wp,&
      & 1.18525618886366E-07_wp, 1.53446017906119E-08_wp, 8.16776989937909E-08_wp,&
      & 9.99771902538170E-09_wp, 2.39574654741861E-08_wp, 3.38775726570832E-08_wp,&
      & 5.00666783855794E-08_wp, 5.17111125726352E-08_wp, 4.51557963362175E-08_wp,&
      & 2.29281184655165E-08_wp, 5.03464173832655E-08_wp, 0.00000000000000E+00_wp,&
      & 5.25937938863503E-09_wp, 4.58923543197251E-07_wp, 1.20784865259198E-09_wp,&
      & 5.49374005833298E-09_wp, 4.60413144915037E-07_wp, 7.52327065924785E-08_wp,&
      & 2.92280423952520E-08_wp, 1.25889825751569E-07_wp, 3.24851298297493E-09_wp,&
      & 2.33025636296754E-07_wp, 2.61954208746454E-07_wp, 5.49430732989573E-09_wp,&
      & 1.02857510308749E-08_wp, 1.02624713797033E-07_wp, 8.85630798816101E-08_wp,&
      & 5.25937938863503E-09_wp, 0.00000000000000E+00_wp, 4.57706086075924E-08_wp,&
      & 2.94692823667863E-07_wp, 5.07503561042260E-08_wp, 8.17256533344214E-08_wp,&
      & 1.24087001773201E-08_wp, 2.59191814116142E-07_wp, 6.19622972587363E-09_wp,&
      & 7.61889722755355E-08_wp, 8.73413584799701E-09_wp, 2.46234976017273E-08_wp,&
      & 1.95339731966989E-07_wp, 9.05471086032081E-08_wp, 5.36673912021823E-08_wp,&
      & 1.25589211938184E-07_wp, 4.58923543197251E-07_wp, 4.57706086075924E-08_wp,&
      & 0.00000000000000E+00_wp, 4.99647294011708E-09_wp, 1.99433747784876E-08_wp,&
      & 4.59479810863950E-07_wp, 2.99306581346282E-08_wp, 1.42450868042157E-07_wp,&
      & 2.64267152803500E-08_wp, 7.58108825911776E-09_wp, 2.33660739390052E-07_wp,&
      & 2.61948144277068E-07_wp, 3.44761371499872E-08_wp, 2.92704923836813E-08_wp,&
      & 1.15117270153361E-07_wp, 4.74868842272753E-08_wp, 1.20784865259198E-09_wp,&
      & 2.94692823667863E-07_wp, 4.99647294011708E-09_wp, 0.00000000000000E+00_wp,&
      & 2.43618825603933E-07_wp, 1.55380314433007E-08_wp, 1.84050996564182E-08_wp,&
      & 2.83057918901539E-07_wp, 5.97508022905764E-09_wp, 3.76879857913355E-07_wp,&
      & 5.25995156093282E-09_wp, 4.81609482223786E-09_wp, 2.15268295349482E-07_wp,&
      & 2.53474457066021E-07_wp, 5.62861345749027E-08_wp, 3.98242609707400E-08_wp,&
      & 5.49374005833298E-09_wp, 5.07503561042260E-08_wp, 1.99433747784876E-08_wp,&
      & 2.43618825603933E-07_wp, 0.00000000000000E+00_wp, 4.25136988771421E-08_wp,&
      & 8.28811198976246E-08_wp, 3.47473864704498E-07_wp, 1.94039661290948E-08_wp,&
      & 4.67896483659639E-07_wp, 1.27586037198218E-07_wp, 2.53469419801484E-08_wp,&
      & 1.93313198634048E-07_wp, 3.55811379624215E-07_wp, 1.45078023243346E-07_wp,&
      & 1.18525618886366E-07_wp, 4.60413144915037E-07_wp, 8.17256533344214E-08_wp,&
      & 4.59479810863950E-07_wp, 1.55380314433007E-08_wp, 4.25136988771421E-08_wp,&
      & 0.00000000000000E+00_wp, 3.11047412983800E-07_wp, 2.89053137025698E-07_wp,&
      & 3.78079410373412E-07_wp, 3.93295940346492E-08_wp, 2.35129533493642E-07_wp,&
      & 2.61935423556151E-07_wp, 9.04672568649688E-08_wp, 1.14196160880360E-07_wp,&
      & 1.55108574863938E-07_wp, 1.53446017906119E-08_wp, 7.52327065924785E-08_wp,&
      & 1.24087001773201E-08_wp, 2.99306581346282E-08_wp, 1.84050996564182E-08_wp,&
      & 8.28811198976246E-08_wp, 3.11047412983800E-07_wp, 0.00000000000000E+00_wp,&
      & 1.75769861391565E-07_wp, 3.71551957504019E-07_wp, 1.94841632302204E-07_wp,&
      & 1.26249438340652E-07_wp, 1.06479113080395E-07_wp, 3.85100702500764E-08_wp,&
      & 2.04700095828997E-07_wp, 1.10677721683708E-07_wp, 8.16776989937909E-08_wp,&
      & 2.92280423952520E-08_wp, 2.59191814116142E-07_wp, 1.42450868042157E-07_wp,&
      & 2.83057918901539E-07_wp, 3.47473864704498E-07_wp, 2.89053137025698E-07_wp,&
      & 1.75769861391565E-07_wp, 0.00000000000000E+00_wp, 7.80168905804317E-08_wp,&
      & 3.47775418273532E-07_wp, 1.21787390455946E-07_wp, 1.15061186295724E-07_wp,&
      & 1.85023145226633E-07_wp, 2.60171954283083E-07_wp, 1.03924251638814E-07_wp,&
      & 9.99771902538170E-09_wp, 1.25889825751569E-07_wp, 6.19622972587363E-09_wp,&
      & 2.64267152803500E-08_wp, 5.97508022905764E-09_wp, 1.94039661290948E-08_wp,&
      & 3.78079410373412E-07_wp, 3.71551957504019E-07_wp, 7.80168905804317E-08_wp,&
      & 0.00000000000000E+00_wp, 4.98180708823302E-08_wp, 1.04689518066842E-07_wp,&
      & 1.26166201047386E-07_wp, 1.47377125812002E-08_wp, 7.74697254860845E-08_wp,&
      & 1.51876405386020E-07_wp, 2.39574654741861E-08_wp, 3.24851298297493E-09_wp,&
      & 7.61889722755355E-08_wp, 7.58108825911776E-09_wp, 3.76879857913355E-07_wp,&
      & 4.67896483659639E-07_wp, 3.93295940346492E-08_wp, 1.94841632302204E-07_wp,&
      & 3.47775418273532E-07_wp, 4.98180708823302E-08_wp, 0.00000000000000E+00_wp,&
      & 3.90797979426098E-08_wp, 1.27914156336358E-08_wp, 2.34126769317346E-07_wp,&
      & 3.55810907381323E-07_wp, 1.46604951268558E-07_wp, 3.38775726570832E-08_wp,&
      & 2.33025636296754E-07_wp, 8.73413584799701E-09_wp, 2.33660739390052E-07_wp,&
      & 5.25995156093282E-09_wp, 1.27586037198218E-07_wp, 2.35129533493642E-07_wp,&
      & 1.26249438340652E-07_wp, 1.21787390455946E-07_wp, 1.04689518066842E-07_wp,&
      & 3.90797979426098E-08_wp, 0.00000000000000E+00_wp, 1.22267842505619E-07_wp,&
      & 1.83306992354658E-08_wp, 1.18147655143357E-07_wp, 6.44265026055596E-08_wp,&
      & 5.00666783855794E-08_wp, 2.61954208746454E-07_wp, 2.46234976017273E-08_wp,&
      & 2.61948144277068E-07_wp, 4.81609482223786E-09_wp, 2.53469419801484E-08_wp,&
      & 2.61935423556151E-07_wp, 1.06479113080395E-07_wp, 1.15061186295724E-07_wp,&
      & 1.26166201047386E-07_wp, 1.27914156336358E-08_wp, 1.22267842505619E-07_wp,&
      & 0.00000000000000E+00_wp, 2.58871007846658E-08_wp, 4.46578643524190E-08_wp,&
      & 6.58004481541880E-08_wp, 5.17111125726352E-08_wp, 5.49430732989573E-09_wp,&
      & 1.95339731966989E-07_wp, 3.44761371499872E-08_wp, 2.15268295349482E-07_wp,&
      & 1.93313198634048E-07_wp, 9.04672568649688E-08_wp, 3.85100702500764E-08_wp,&
      & 1.85023145226633E-07_wp, 1.47377125812002E-08_wp, 2.34126769317346E-07_wp,&
      & 1.83306992354658E-08_wp, 2.58871007846658E-08_wp, 0.00000000000000E+00_wp,&
      & 1.73354790258678E-07_wp, 6.03103883856771E-08_wp, 4.51557963362175E-08_wp,&
      & 1.02857510308749E-08_wp, 9.05471086032081E-08_wp, 2.92704923836813E-08_wp,&
      & 2.53474457066021E-07_wp, 3.55811379624215E-07_wp, 1.14196160880360E-07_wp,&
      & 2.04700095828997E-07_wp, 2.60171954283083E-07_wp, 7.74697254860845E-08_wp,&
      & 3.55810907381323E-07_wp, 1.18147655143357E-07_wp, 4.46578643524190E-08_wp,&
      & 1.73354790258678E-07_wp, 0.00000000000000E+00_wp, 1.03746059757864E-07_wp,&
      & 2.29281184655165E-08_wp, 1.02624713797033E-07_wp, 5.36673912021823E-08_wp,&
      & 1.15117270153361E-07_wp, 5.62861345749027E-08_wp, 1.45078023243346E-07_wp,&
      & 1.55108574863938E-07_wp, 1.10677721683708E-07_wp, 1.03924251638814E-07_wp,&
      & 1.51876405386020E-07_wp, 1.46604951268558E-07_wp, 6.44265026055596E-08_wp,&
      & 6.58004481541880E-08_wp, 6.03103883856771E-08_wp, 1.03746059757864E-07_wp,&
      & 0.00000000000000E+00_wp], shape(ref8))

   ! PBE-D4(BJ)
   param = param_type(&
      & s6 = 1.0_wp, s8 = 0.95948085_wp, a1 = 0.38574991_wp, a2 = 4.80688534_wp)

   call get_structure(mol, "MB16-43", "01")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%rational, -1)
   if (allocated(error)) return
   call test_damping_2b_gen(error, mol, d4, damp, param, ref6, ref8)

end subroutine test_damp_rational_2b_mb01

subroutine test_damp_rational_2b_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   real(wp), parameter :: ref6(16, 16) = reshape([&
      & 0.00000000000000E+00_wp, 9.24455324085766E-06_wp, 4.17389082745223E-06_wp,&
      & 1.22906563552554E-05_wp, 2.24949391527265E-06_wp, 1.80531872454847E-05_wp,&
      & 1.74431677764355E-05_wp, 1.41160069860077E-06_wp, 5.41854475837807E-06_wp,&
      & 1.12025816041255E-05_wp, 1.19248475290049E-05_wp, 5.53682099605787E-06_wp,&
      & 1.62666196361523E-06_wp, 3.04985669941946E-06_wp, 1.77231662164648E-05_wp,&
      & 2.82937867236546E-06_wp, 9.24455324085766E-06_wp, 0.00000000000000E+00_wp,&
      & 6.15069103896672E-06_wp, 5.36510942388375E-06_wp, 3.67444270722425E-06_wp,&
      & 7.30193280547424E-06_wp, 9.44415986607326E-06_wp, 9.07173850747533E-06_wp,&
      & 3.96673281944419E-06_wp, 5.34283497753362E-06_wp, 6.14448799614771E-06_wp,&
      & 4.13001542025359E-06_wp, 5.97482812926627E-06_wp, 2.75633023546466E-06_wp,&
      & 9.22230546749978E-06_wp, 2.74175094680576E-06_wp, 4.17389082745223E-06_wp,&
      & 6.15069103896672E-06_wp, 0.00000000000000E+00_wp, 2.12021610982332E-06_wp,&
      & 4.31335636547904E-06_wp, 2.61527483953398E-06_wp, 7.69930373222560E-06_wp,&
      & 1.19268429880358E-05_wp, 5.01517761254365E-06_wp, 2.94486531995931E-06_wp,&
      & 6.00409504356252E-06_wp, 2.74626172255652E-06_wp, 4.10627284144087E-06_wp,&
      & 7.88035422070603E-06_wp, 9.52234658694059E-06_wp, 4.31019422608816E-06_wp,&
      & 1.22906563552554E-05_wp, 5.36510942388375E-06_wp, 2.12021610982332E-06_wp,&
      & 0.00000000000000E+00_wp, 4.44050786288524E-06_wp, 1.55143053576223E-05_wp,&
      & 2.46982937985267E-06_wp, 6.72254567451801E-07_wp, 3.08274731057595E-06_wp,&
      & 1.55138233806164E-05_wp, 6.35291421638953E-06_wp, 6.95901367087207E-06_wp,&
      & 8.71105221985713E-06_wp, 1.05531327118784E-06_wp, 8.51367515276888E-06_wp,&
      & 5.57331223717342E-06_wp, 2.24949391527265E-06_wp, 3.67444270722425E-06_wp,&
      & 4.31335636547904E-06_wp, 4.44050786288524E-06_wp, 0.00000000000000E+00_wp,&
      & 3.89487654205009E-06_wp, 1.54374187139748E-06_wp, 6.00437018378008E-06_wp,&
      & 2.94248689928276E-06_wp, 6.61279924194428E-06_wp, 2.70810665746670E-06_wp,&
      & 2.88291770572169E-06_wp, 6.96863072734967E-06_wp, 2.63570252815228E-06_wp,&
      & 4.66843859787461E-06_wp, 3.84233811283023E-06_wp, 1.80531872454847E-05_wp,&
      & 7.30193280547424E-06_wp, 2.61527483953398E-06_wp, 1.55143053576223E-05_wp,&
      & 3.89487654205009E-06_wp, 0.00000000000000E+00_wp, 6.22417494678132E-06_wp,&
      & 7.77704111783744E-07_wp, 4.26091010532200E-06_wp, 1.83140911950459E-05_wp,&
      & 1.04709904661254E-05_wp, 8.35292491208758E-06_wp, 5.06615934671203E-06_wp,&
      & 1.62392734951728E-06_wp, 1.43592811510725E-05_wp, 5.39469957751509E-06_wp,&
      & 1.74431677764355E-05_wp, 9.44415986607326E-06_wp, 7.69930373222560E-06_wp,&
      & 2.46982937985267E-06_wp, 1.54374187139748E-06_wp, 6.22417494678132E-06_wp,&
      & 0.00000000000000E+00_wp, 4.62904384828798E-06_wp, 6.82984296252412E-06_wp,&
      & 2.24982717190279E-06_wp, 1.19264510667878E-05_wp, 1.77813158843865E-06_wp,&
      & 6.38782298599451E-07_wp, 7.05260704712990E-06_wp, 1.75820269738059E-05_wp,&
      & 1.61080400653750E-06_wp, 1.41160069860077E-06_wp, 9.07173850747533E-06_wp,&
      & 1.19268429880358E-05_wp, 6.72254567451801E-07_wp, 6.00437018378008E-06_wp,&
      & 7.77704111783744E-07_wp, 4.62904384828798E-06_wp, 0.00000000000000E+00_wp,&
      & 7.22705803001736E-06_wp, 8.72251794675711E-07_wp, 4.52477835911936E-06_wp,&
      & 1.47281893976449E-06_wp, 2.15016273853699E-06_wp, 5.17462582896099E-06_wp,&
      & 4.80354211894679E-06_wp, 2.37030393013641E-06_wp, 5.41854475837807E-06_wp,&
      & 3.96673281944419E-06_wp, 5.01517761254365E-06_wp, 3.08274731057595E-06_wp,&
      & 2.94248689928276E-06_wp, 4.26091010532200E-06_wp, 6.82984296252412E-06_wp,&
      & 7.22705803001736E-06_wp, 0.00000000000000E+00_wp, 5.44778570910793E-06_wp,&
      & 4.75036083500940E-06_wp, 1.71906411002963E-06_wp, 2.12307897718637E-06_wp,&
      & 9.23495995360032E-06_wp, 9.20318153323740E-06_wp, 4.44529382472647E-06_wp,&
      & 1.12025816041255E-05_wp, 5.34283497753362E-06_wp, 2.94486531995931E-06_wp,&
      & 1.55138233806164E-05_wp, 6.61279924194428E-06_wp, 1.83140911950459E-05_wp,&
      & 2.24982717190279E-06_wp, 8.72251794675711E-07_wp, 5.44778570910793E-06_wp,&
      & 0.00000000000000E+00_wp, 6.70114339676236E-06_wp, 8.35804958406384E-06_wp,&
      & 1.00407761831985E-05_wp, 2.12343221407206E-06_wp, 1.17013915256564E-05_wp,&
      & 9.90429754947780E-06_wp, 1.19248475290049E-05_wp, 6.14448799614771E-06_wp,&
      & 6.00409504356252E-06_wp, 6.35291421638953E-06_wp, 2.70810665746670E-06_wp,&
      & 1.04709904661254E-05_wp, 1.19264510667878E-05_wp, 4.52477835911936E-06_wp,&
      & 4.75036083500940E-06_wp, 6.70114339676236E-06_wp, 0.00000000000000E+00_wp,&
      & 3.43483680929703E-06_wp, 1.86919004742358E-06_wp, 6.22505286850088E-06_wp,&
      & 1.19087771605109E-05_wp, 3.16119553209217E-06_wp, 5.53682099605787E-06_wp,&
      & 4.13001542025359E-06_wp, 2.74626172255652E-06_wp, 6.95901367087207E-06_wp,&
      & 2.88291770572169E-06_wp, 8.35292491208758E-06_wp, 1.77813158843865E-06_wp,&
      & 1.47281893976449E-06_wp, 1.71906411002963E-06_wp, 8.35804958406384E-06_wp,&
      & 3.43483680929703E-06_wp, 0.00000000000000E+00_wp, 7.58883857928991E-06_wp,&
      & 5.68911747983374E-07_wp, 3.56966016651432E-06_wp, 2.42295965196860E-06_wp,&
      & 1.62666196361523E-06_wp, 5.97482812926627E-06_wp, 4.10627284144087E-06_wp,&
      & 8.71105221985713E-06_wp, 6.96863072734967E-06_wp, 5.06615934671203E-06_wp,&
      & 6.38782298599451E-07_wp, 2.15016273853699E-06_wp, 2.12307897718637E-06_wp,&
      & 1.00407761831985E-05_wp, 1.86919004742358E-06_wp, 7.58883857928991E-06_wp,&
      & 0.00000000000000E+00_wp, 5.08814270555239E-07_wp, 1.83326716256904E-06_wp,&
      & 4.43789947063537E-06_wp, 3.04985669941946E-06_wp, 2.75633023546466E-06_wp,&
      & 7.88035422070603E-06_wp, 1.05531327118784E-06_wp, 2.63570252815228E-06_wp,&
      & 1.62392734951728E-06_wp, 7.05260704712990E-06_wp, 5.17462582896099E-06_wp,&
      & 9.23495995360032E-06_wp, 2.12343221407206E-06_wp, 6.22505286850088E-06_wp,&
      & 5.68911747983374E-07_wp, 5.08814270555239E-07_wp, 0.00000000000000E+00_wp,&
      & 1.62313402225260E-05_wp, 8.01337021340942E-06_wp, 1.77231662164648E-05_wp,&
      & 9.22230546749978E-06_wp, 9.52234658694059E-06_wp, 8.51367515276888E-06_wp,&
      & 4.66843859787461E-06_wp, 1.43592811510725E-05_wp, 1.75820269738059E-05_wp,&
      & 4.80354211894679E-06_wp, 9.20318153323740E-06_wp, 1.17013915256564E-05_wp,&
      & 1.19087771605109E-05_wp, 3.56966016651432E-06_wp, 1.83326716256904E-06_wp,&
      & 1.62313402225260E-05_wp, 0.00000000000000E+00_wp, 8.53853842081169E-06_wp,&
      & 2.82937867236546E-06_wp, 2.74175094680576E-06_wp, 4.31019422608816E-06_wp,&
      & 5.57331223717342E-06_wp, 3.84233811283023E-06_wp, 5.39469957751509E-06_wp,&
      & 1.61080400653750E-06_wp, 2.37030393013641E-06_wp, 4.44529382472647E-06_wp,&
      & 9.90429754947780E-06_wp, 3.16119553209217E-06_wp, 2.42295965196860E-06_wp,&
      & 4.43789947063537E-06_wp, 8.01337021340942E-06_wp, 8.53853842081169E-06_wp,&
      & 0.00000000000000E+00_wp], shape(ref6))

   real(wp), parameter :: ref8(16, 16) = reshape([&
      & 0.00000000000000E+00_wp, 2.11210483590081E-07_wp, 7.96208340706285E-08_wp,&
      & 3.17971089881000E-07_wp, 3.39857912171813E-08_wp, 4.66494380021760E-07_wp,&
      & 4.58877150877435E-07_wp, 1.63019363789754E-08_wp, 1.14330554623730E-07_wp,&
      & 3.00066852243986E-07_wp, 2.61946863193696E-07_wp, 1.17096857293223E-07_wp,&
      & 2.00354462984439E-08_wp, 4.83917894620377E-08_wp, 4.62671574716581E-07_wp,&
      & 4.59332369570651E-08_wp, 2.11210483590081E-07_wp, 0.00000000000000E+00_wp,&
      & 1.09751691011256E-07_wp, 1.12212264466722E-07_wp, 5.83168200179793E-08_wp,&
      & 1.67095448854070E-07_wp, 2.14824809077459E-07_wp, 2.07918655779563E-07_wp,&
      & 6.69355302266826E-08_wp, 1.12897000972169E-07_wp, 1.09715512206864E-07_wp,&
      & 6.72326743215182E-08_wp, 1.28957386307585E-07_wp, 4.42306943392841E-08_wp,&
      & 2.10794989373786E-07_wp, 4.63856206301636E-08_wp, 7.96208340706285E-08_wp,&
      & 1.09751691011256E-07_wp, 0.00000000000000E+00_wp, 3.07346005256722E-08_wp,&
      & 6.99646903805197E-08_wp, 4.06127419800053E-08_wp, 1.80409970973696E-07_wp,&
      & 2.61953761071953E-07_wp, 8.36450653756568E-08_wp, 4.82024288591218E-08_wp,&
      & 1.17110224244693E-07_wp, 4.63578730652473E-08_wp, 7.87608883925409E-08_wp,&
      & 1.85329402168888E-07_wp, 2.25620679957650E-07_wp, 8.20886351541257E-08_wp,&
      & 3.17971089881000E-07_wp, 1.12212264466722E-07_wp, 3.07346005256722E-08_wp,&
      & 0.00000000000000E+00_wp, 8.56162728123089E-08_wp, 3.71549549073827E-07_wp,&
      & 3.64041324432481E-08_wp, 5.89883564520814E-09_wp, 5.39593746371130E-08_wp,&
      & 3.71548223303573E-07_wp, 1.38933827612862E-07_wp, 1.29729346806790E-07_wp,&
      & 2.12997196191633E-07_wp, 1.09885490189436E-08_wp, 2.09715618545862E-07_wp,&
      & 1.17422158299238E-07_wp, 3.39857912171813E-08_wp, 5.83168200179793E-08_wp,&
      & 6.99646903805197E-08_wp, 8.56162728123089E-08_wp, 0.00000000000000E+00_wp,&
      & 7.40041461443801E-08_wp, 1.97479990387609E-08_wp, 1.26471751577837E-07_wp,&
      & 4.30914877036311E-08_wp, 1.39194467995283E-07_wp, 4.51928418509007E-08_wp,&
      & 4.14646642230336E-08_wp, 1.29514974551665E-07_wp, 4.27002378052075E-08_wp,&
      & 9.40351189902575E-08_wp, 5.97098709349609E-08_wp, 4.66494380021760E-07_wp,&
      & 1.67095448854070E-07_wp, 4.06127419800053E-08_wp, 3.71549549073827E-07_wp,&
      & 7.40041461443801E-08_wp, 0.00000000000000E+00_wp, 1.35147214938322E-07_wp,&
      & 7.15548549533804E-09_wp, 8.33828121848827E-08_wp, 4.68872233547726E-07_wp,&
      & 2.44118264463100E-07_wp, 1.74522068717940E-07_wp, 1.01792096333465E-07_wp,&
      & 1.98255745461404E-08_wp, 3.94470205748802E-07_wp, 1.14370274235708E-07_wp,&
      & 4.58877150877435E-07_wp, 2.14824809077459E-07_wp, 1.80409970973696E-07_wp,&
      & 3.64041324432481E-08_wp, 1.97479990387609E-08_wp, 1.35147214938322E-07_wp,&
      & 0.00000000000000E+00_wp, 8.82349888493819E-08_wp, 1.49089751879610E-07_wp,&
      & 3.13746247687834E-08_wp, 2.61952425139452E-07_wp, 2.40005583246208E-08_wp,&
      & 5.48846389894480E-09_wp, 1.61450627079454E-07_wp, 4.60816246581173E-07_wp,&
      & 2.04356947185081E-08_wp, 1.63019363789754E-08_wp, 2.07918655779563E-07_wp,&
      & 2.61953761071953E-07_wp, 5.89883564520814E-09_wp, 1.26471751577837E-07_wp,&
      & 7.15548549533804E-09_wp, 8.82349888493819E-08_wp, 0.00000000000000E+00_wp,&
      & 1.57650405773702E-07_wp, 8.37575348623444E-09_wp, 8.92983255273805E-08_wp,&
      & 1.83136275314924E-08_wp, 2.97174138743047E-08_wp, 1.03630866891566E-07_wp,&
      & 9.30777542299611E-08_wp, 3.55745864766774E-08_wp, 1.14330554623730E-07_wp,&
      & 6.69355302266826E-08_wp, 8.36450653756568E-08_wp, 5.39593746371130E-08_wp,&
      & 4.30914877036311E-08_wp, 8.33828121848827E-08_wp, 1.49089751879610E-07_wp,&
      & 1.57650405773702E-07_wp, 0.00000000000000E+00_wp, 1.15098734396820E-07_wp,&
      & 8.15463115770116E-08_wp, 2.48925891831949E-08_wp, 3.14665877813885E-08_wp,&
      & 1.86891975321464E-07_wp, 1.86733067499200E-07_wp, 7.15154037373456E-08_wp,&
      & 3.00066852243986E-07_wp, 1.12897000972169E-07_wp, 4.82024288591218E-08_wp,&
      & 3.71548223303573E-07_wp, 1.39194467995283E-07_wp, 4.68872233547726E-07_wp,&
      & 3.13746247687834E-08_wp, 8.37575348623444E-09_wp, 1.15098734396820E-07_wp,&
      & 0.00000000000000E+00_wp, 1.52221605989081E-07_wp, 1.74587023731427E-07_wp,&
      & 2.59064050163494E-07_wp, 2.89100513249767E-08_wp, 3.16216514908993E-07_wp,&
      & 2.22322838013796E-07_wp, 2.61946863193696E-07_wp, 1.09715512206864E-07_wp,&
      & 1.17110224244693E-07_wp, 1.38933827612862E-07_wp, 4.51928418509007E-08_wp,&
      & 2.44118264463100E-07_wp, 2.61952425139452E-07_wp, 8.92983255273805E-08_wp,&
      & 8.15463115770116E-08_wp, 1.52221605989081E-07_wp, 0.00000000000000E+00_wp,&
      & 6.05426185253787E-08_wp, 2.54365641919303E-08_wp, 1.38377264182846E-07_wp,&
      & 2.61883891570842E-07_wp, 5.60827444820214E-08_wp, 1.17096857293223E-07_wp,&
      & 6.72326743215182E-08_wp, 4.63578730652473E-08_wp, 1.29729346806790E-07_wp,&
      & 4.14646642230336E-08_wp, 1.74522068717940E-07_wp, 2.40005583246208E-08_wp,&
      & 1.83136275314924E-08_wp, 2.48925891831949E-08_wp, 1.74587023731427E-07_wp,&
      & 6.05426185253787E-08_wp, 0.00000000000000E+00_wp, 1.44865542878192E-07_wp,&
      & 4.80129545913076E-09_wp, 6.53030425997779E-08_wp, 3.92444925908198E-08_wp,&
      & 2.00354462984439E-08_wp, 1.28957386307585E-07_wp, 7.87608883925409E-08_wp,&
      & 2.12997196191633E-07_wp, 1.29514974551665E-07_wp, 1.01792096333465E-07_wp,&
      & 5.48846389894480E-09_wp, 2.97174138743047E-08_wp, 3.14665877813885E-08_wp,&
      & 2.59064050163494E-07_wp, 2.54365641919303E-08_wp, 1.44865542878192E-07_wp,&
      & 0.00000000000000E+00_wp, 4.02338600616182E-09_wp, 2.37113371142762E-08_wp,&
      & 8.80995909483606E-08_wp, 4.83917894620377E-08_wp, 4.42306943392841E-08_wp,&
      & 1.85329402168888E-07_wp, 1.09885490189436E-08_wp, 4.27002378052075E-08_wp,&
      & 1.98255745461404E-08_wp, 1.61450627079454E-07_wp, 1.03630866891566E-07_wp,&
      & 1.86891975321464E-07_wp, 2.89100513249767E-08_wp, 1.38377264182846E-07_wp,&
      & 4.80129545913076E-09_wp, 4.02338600616182E-09_wp, 0.00000000000000E+00_wp,&
      & 4.37947692641701E-07_wp, 1.84872232299874E-07_wp, 4.62671574716581E-07_wp,&
      & 2.10794989373786E-07_wp, 2.25620679957650E-07_wp, 2.09715618545862E-07_wp,&
      & 9.40351189902575E-08_wp, 3.94470205748802E-07_wp, 4.60816246581173E-07_wp,&
      & 9.30777542299611E-08_wp, 1.86733067499200E-07_wp, 3.16216514908993E-07_wp,&
      & 2.61883891570842E-07_wp, 6.53030425997779E-08_wp, 2.37113371142762E-08_wp,&
      & 4.37947692641701E-07_wp, 0.00000000000000E+00_wp, 1.96887119769173E-07_wp,&
      & 4.59332369570651E-08_wp, 4.63856206301636E-08_wp, 8.20886351541257E-08_wp,&
      & 1.17422158299238E-07_wp, 5.97098709349609E-08_wp, 1.14370274235708E-07_wp,&
      & 2.04356947185081E-08_wp, 3.55745864766774E-08_wp, 7.15154037373456E-08_wp,&
      & 2.22322838013796E-07_wp, 5.60827444820214E-08_wp, 3.92444925908198E-08_wp,&
      & 8.80995909483606E-08_wp, 1.84872232299874E-07_wp, 1.96887119769173E-07_wp,&
      & 0.00000000000000E+00_wp], shape(ref8))

   ! PBE-D4(BJ)
   param = param_type(&
      & s6 = 1.0_wp, s8 = 0.95948085_wp, a1 = 0.38574991_wp, a2 = 4.80688534_wp)

   call get_structure(mol, "MB16-43", "02")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%rational, -1)
   if (allocated(error)) return
   call test_damping_2b_gen(error, mol, d4, damp, param, ref6, ref8)

end subroutine test_damp_rational_2b_mb02

subroutine test_grad_rational_2b_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   ! PBE-D4(BJ)
   param = param_type(&
      & s6 = 1.0_wp, s8 = 0.95948085_wp, a1 = 0.38574991_wp, a2 = 4.80688534_wp)

   call get_structure(mol, "MB16-43", "03")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%rational, -1)
   if (allocated(error)) return
   call test_damping_2b_numgrad(error, mol, d4, damp, param)

end subroutine test_grad_rational_2b_mb03

subroutine test_grad_rational_2b_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   ! PBE-D4(BJ)
   param = param_type(&
      & s6 = 1.0_wp, s8 = 0.95948085_wp, a1 = 0.38574991_wp, a2 = 4.80688534_wp)

   call get_structure(mol, "MB16-43", "04")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%rational, -1)
   if (allocated(error)) return
   call test_damping_2b_numgrad(error, mol, d4, damp, param)

end subroutine test_grad_rational_2b_mb04


subroutine test_damp_screened_2b_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   real(wp), parameter :: ref6(16, 16) = reshape([&
      & 0.00000000000000E+00_wp, 5.03534916438524E-06_wp, 2.68346738443249E-06_wp,&
      & 1.15352946921458E-05_wp, 5.30906721979157E-06_wp, 3.84581103887485E-06_wp,&
      & 2.10853362264733E-05_wp, 1.63847652035871E-06_wp, 4.91919637075494E-06_wp,&
      & 1.09639740962485E-06_wp, 2.32133123207574E-06_wp, 4.96877775437257E-06_wp,&
      & 3.14045139299545E-06_wp, 1.53600250014087E-06_wp, 5.57531113851225E-06_wp,&
      & 7.35843877543454E-07_wp, 5.03534916438524E-06_wp, 0.00000000000000E+00_wp,&
      & 6.43897886169665E-07_wp, 1.65683545423019E-04_wp, 2.11817960065162E-07_wp,&
      & 6.64053889636435E-07_wp, 1.70690982679544E-04_wp, 5.55206703698342E-06_wp,&
      & 2.45889349372963E-06_wp, 8.70902858496695E-06_wp, 4.46164686373305E-07_wp,&
      & 5.08046170386433E-05_wp, 1.03992116931341E-05_wp, 6.68822178660848E-07_wp,&
      & 1.07697039492201E-06_wp, 1.27732455059780E-05_wp, 2.68346738443249E-06_wp,&
      & 6.43897886169665E-07_wp, 0.00000000000000E+00_wp, 3.56229117700954E-06_wp,&
      & 7.56160765484231E-05_wp, 3.89390074484321E-06_wp, 6.00658936551557E-06_wp,&
      & 1.25364832069415E-06_wp, 6.55876590160280E-05_wp, 7.29534554267928E-07_wp,&
      & 5.61850407128601E-06_wp, 9.65408208881092E-07_wp, 2.24314054532269E-06_wp,&
      & 5.64801113398239E-06_wp, 7.27650717340789E-06_wp, 5.91488634504335E-06_wp,&
      & 1.15352946921458E-05_wp, 1.65683545423019E-04_wp, 3.56229117700954E-06_wp,&
      & 0.00000000000000E+00_wp, 6.18794208364236E-07_wp, 1.78831446403148E-06_wp,&
      & 1.67580536738972E-04_wp, 2.49985296053724E-06_wp, 1.10898654161329E-05_wp,&
      & 2.23281790544317E-06_wp, 8.48335822099510E-07_wp, 5.14243444663048E-05_wp,&
      & 1.07276082796066E-05_wp, 2.90111451467537E-06_wp, 2.46181140089688E-06_wp,&
      & 1.63521964640480E-05_wp, 5.30906721979157E-06_wp, 2.11817960065162E-07_wp,&
      & 7.56160765484231E-05_wp, 6.18794208364236E-07_wp, 0.00000000000000E+00_wp,&
      & 2.26242977220302E-05_wp, 1.47865136616739E-06_wp, 1.70586425169611E-06_wp,&
      & 7.37021891884309E-05_wp, 7.08941969235615E-07_wp, 9.72897136251241E-05_wp,&
      & 6.49359069737587E-07_wp, 6.06541151511531E-07_wp, 8.44197962215214E-06_wp,&
      & 4.35215543746399E-05_wp, 5.89360411831741E-06_wp, 3.84581103887485E-06_wp,&
      & 6.64053889636435E-07_wp, 3.89390074484321E-06_wp, 1.78831446403148E-06_wp,&
      & 2.26242977220302E-05_wp, 0.00000000000000E+00_wp, 3.27905713233834E-06_wp,&
      & 6.08836290160338E-06_wp, 1.08840306946438E-04_wp, 1.75034094994188E-06_wp,&
      & 1.61892019481197E-04_wp, 1.16559601230590E-05_wp, 2.23644504054718E-06_wp,&
      & 2.50113959025705E-05_wp, 1.84576380140450E-05_wp, 2.86522572157955E-05_wp,&
      & 2.10853362264733E-05_wp, 1.70690982679544E-04_wp, 6.00658936551557E-06_wp,&
      & 1.67580536738972E-04_wp, 1.47865136616739E-06_wp, 3.27905713233834E-06_wp,&
      & 0.00000000000000E+00_wp, 4.96150371576745E-05_wp, 4.27831067302917E-05_wp,&
      & 5.17254813710344E-05_wp, 3.07591991079842E-06_wp, 5.28820604833307E-05_wp,&
      & 1.13757099970310E-05_wp, 7.38824537822471E-06_wp, 8.56338638305307E-06_wp,&
      & 9.95573787135679E-06_wp, 1.63847652035871E-06_wp, 5.55206703698342E-06_wp,&
      & 1.25364832069415E-06_wp, 2.49985296053724E-06_wp, 1.70586425169611E-06_wp,&
      & 6.08836290160338E-06_wp, 4.96150371576745E-05_wp, 0.00000000000000E+00_wp,&
      & 1.88953759391406E-05_wp, 1.89357089632387E-05_wp, 1.66528513095078E-05_wp,&
      & 1.52235547478993E-05_wp, 1.08780037788119E-05_wp, 3.34293320357185E-06_wp,&
      & 2.68949564728784E-05_wp, 2.34664804968978E-06_wp, 4.91919637075494E-06_wp,&
      & 2.45889349372963E-06_wp, 6.55876590160280E-05_wp, 1.10898654161329E-05_wp,&
      & 7.37021891884309E-05_wp, 1.08840306946438E-04_wp, 4.27831067302917E-05_wp,&
      & 1.88953759391406E-05_wp, 0.00000000000000E+00_wp, 5.79707020682059E-06_wp,&
      & 1.09118040289588E-04_wp, 1.51460718063339E-05_wp, 1.30901984402942E-05_wp,&
      & 6.64109822278374E-06_wp, 2.08635433861639E-05_wp, 2.30299581645720E-06_wp,&
      & 1.09639740962485E-06_wp, 8.70902858496695E-06_wp, 7.29534554267928E-07_wp,&
      & 2.23281790544317E-06_wp, 7.08941969235615E-07_wp, 1.75034094994188E-06_wp,&
      & 5.17254813710344E-05_wp, 1.89357089632387E-05_wp, 5.79707020682059E-06_wp,&
      & 0.00000000000000E+00_wp, 3.74115304527744E-06_wp, 8.90278822001552E-06_wp,&
      & 1.12669210136344E-05_wp, 1.44096289474786E-06_wp, 5.75805642969075E-06_wp,&
      & 2.34758650472405E-05_wp, 2.32133123207574E-06_wp, 4.46164686373305E-07_wp,&
      & 5.61850407128601E-06_wp, 8.48335822099510E-07_wp, 9.72897136251241E-05_wp,&
      & 1.61892019481197E-04_wp, 3.07591991079842E-06_wp, 1.66528513095078E-05_wp,&
      & 1.09118040289588E-04_wp, 3.74115304527744E-06_wp, 0.00000000000000E+00_wp,&
      & 3.24720289690559E-06_wp, 1.28814685728107E-06_wp, 4.76403436825604E-05_wp,&
      & 1.85366239164814E-05_wp, 2.86352079094568E-05_wp, 4.96877775437256E-06_wp,&
      & 5.08046170386433E-05_wp, 9.65408208881092E-07_wp, 5.14243444663048E-05_wp,&
      & 6.49359069737587E-07_wp, 1.16559601230590E-05_wp, 5.28820604833307E-05_wp,&
      & 1.52235547478993E-05_wp, 1.51460718063339E-05_wp, 8.90278822001552E-06_wp,&
      & 3.24720289690559E-06_wp, 0.00000000000000E+00_wp, 5.56302714229750E-06_wp,&
      & 1.83414982505854E-06_wp, 1.42313763547916E-05_wp, 1.35474767921758E-06_wp,&
      & 3.14045139299545E-06_wp, 1.03992116931341E-05_wp, 2.24314054532269E-06_wp,&
      & 1.07276082796066E-05_wp, 6.06541151511531E-07_wp, 2.23644504054718E-06_wp,&
      & 1.13757099970310E-05_wp, 1.08780037788119E-05_wp, 1.30901984402942E-05_wp,&
      & 1.12669210136344E-05_wp, 1.28814685728107E-06_wp, 5.56302714229750E-06_wp,&
      & 0.00000000000000E+00_wp, 2.49800874011755E-06_wp, 3.89687935545383E-06_wp,&
      & 2.89113736324840E-06_wp, 1.53600250014087E-06_wp, 6.68822178660848E-07_wp,&
      & 5.64801113398239E-06_wp, 2.90111451467537E-06_wp, 8.44197962215214E-06_wp,&
      & 2.50113959025705E-05_wp, 7.38824537822471E-06_wp, 3.34293320357185E-06_wp,&
      & 6.64109822278374E-06_wp, 1.44096289474786E-06_wp, 4.76403436825604E-05_wp,&
      & 1.83414982505854E-06_wp, 2.49800874011755E-06_wp, 0.00000000000000E+00_wp,&
      & 3.70277903206410E-05_wp, 7.87359974019103E-06_wp, 5.57531113851225E-06_wp,&
      & 1.07697039492201E-06_wp, 7.27650717340789E-06_wp, 2.46181140089688E-06_wp,&
      & 4.35215543746399E-05_wp, 1.84576380140450E-05_wp, 8.56338638305307E-06_wp,&
      & 2.68949564728784E-05_wp, 2.08635433861639E-05_wp, 5.75805642969075E-06_wp,&
      & 1.85366239164814E-05_wp, 1.42313763547916E-05_wp, 3.89687935545383E-06_wp,&
      & 3.70277903206410E-05_wp, 0.00000000000000E+00_wp, 2.95794775698417E-06_wp,&
      & 7.35843877543454E-07_wp, 1.27732455059780E-05_wp, 5.91488634504335E-06_wp,&
      & 1.63521964640480E-05_wp, 5.89360411831741E-06_wp, 2.86522572157954E-05_wp,&
      & 9.95573787135676E-06_wp, 2.34664804968978E-06_wp, 2.30299581645719E-06_wp,&
      & 2.34758650472404E-05_wp, 2.86352079094567E-05_wp, 1.35474767921758E-06_wp,&
      & 2.89113736324840E-06_wp, 7.87359974019103E-06_wp, 2.95794775698416E-06_wp,&
      & 0.00000000000000E+00_wp], shape(ref6))

   real(wp), parameter :: ref8(16, 16) = reshape([&
      & 0.00000000000000E+00_wp, 8.28086617060915E-08_wp, 3.57792642799852E-08_wp,&
      & 2.50077471589419E-07_wp, 8.88642902739787E-08_wp, 5.78122981082288E-08_wp,&
      & 5.58912323371465E-07_wp, 1.85334643734401E-08_wp, 8.02715885043798E-08_wp,&
      & 1.08474330641728E-08_wp, 2.94907486282681E-08_wp, 8.13521572303127E-08_wp,&
      & 4.41257817307720E-08_wp, 1.70043045502458E-08_wp, 9.48553410038026E-08_wp,&
      & 6.37407572831104E-09_wp, 8.28086617060915E-08_wp, 0.00000000000000E+00_wp,&
      & 5.33489214409222E-09_wp, 8.73123110603287E-06_wp, 1.21149051824341E-09_wp,&
      & 5.55871077294144E-09_wp, 9.08483581069493E-06_wp, 9.43284242281929E-08_wp,&
      & 3.18436323843386E-08_wp, 1.71921133506486E-07_wp, 3.27112069616806E-09_wp,&
      & 1.80540186703872E-06_wp, 2.17789465958412E-07_wp, 5.61199401868333E-09_wp,&
      & 1.05919193336971E-08_wp, 2.86486826906670E-07_wp, 3.57792642799852E-08_wp,&
      & 5.33489214409222E-09_wp, 0.00000000000000E+00_wp, 5.22006044058733E-08_wp,&
      & 3.06799864877664E-06_wp, 5.87781817541937E-08_wp, 1.04762752805643E-07_wp,&
      & 1.29699193497879E-08_wp, 2.53784923205909E-06_wp, 6.30130934390886E-09_wp,&
      & 9.58364190274570E-08_wp, 9.15484426760728E-09_wp, 2.81737695822022E-08_wp,&
      & 9.65080861967912E-08_wp, 1.35290383223564E-07_wp, 1.02635637280416E-07_wp,&
      & 2.50077471589419E-07_wp, 8.73123110603287E-06_wp, 5.22006044058733E-08_wp,&
      & 0.00000000000000E+00_wp, 5.05938774969931E-09_wp, 2.08270693635819E-08_wp,&
      & 8.86477561659876E-06_wp, 3.25528435220703E-08_wp, 2.37285579341720E-07_wp,&
      & 2.80010329601741E-08_wp, 7.70536801679902E-09_wp, 1.83482506171271E-06_wp,&
      & 2.27007490088747E-07_wp, 3.96999055675365E-08_wp, 3.18940263985877E-08_wp,&
      & 3.98233483719282E-07_wp, 8.88642902739787E-08_wp, 1.21149051824341E-09_wp,&
      & 3.06799864877664E-06_wp, 5.05938774969931E-09_wp, 0.00000000000000E+00_wp,&
      & 6.13954877811801E-07_wp, 1.61630765462122E-08_wp, 1.95567012770084E-08_wp,&
      & 2.96490076287824E-06_wp, 6.06527589798159E-09_wp, 4.29330963668247E-06_wp,&
      & 5.39530738035960E-09_wp, 4.92625255487705E-09_wp, 1.64928373932152E-07_wp,&
      & 1.46884344099167E-06_wp, 1.02143544654681E-07_wp, 5.78122981082288E-08_wp,&
      & 5.55871077294144E-09_wp, 5.87781817541937E-08_wp, 2.08270693635819E-08_wp,&
      & 6.13954877811801E-07_wp, 0.00000000000000E+00_wp, 4.67413899309928E-08_wp,&
      & 1.06668704158924E-07_wp, 4.98604264193123E-06_wp, 2.02395032461751E-08_wp,&
      & 8.46584342537644E-06_wp, 2.53571466203672E-07_wp, 2.80616983147280E-08_wp,&
      & 7.01811077260374E-07_wp, 4.68027254363654E-07_wp, 8.41230027213864E-07_wp,&
      & 5.58912323371465E-07_wp, 9.08483581069493E-06_wp, 1.04762752805643E-07_wp,&
      & 8.86477561659876E-06_wp, 1.61630765462122E-08_wp, 4.67413899309928E-08_wp,&
      & 0.00000000000000E+00_wp, 1.74925882038971E-06_wp, 1.43570782883665E-06_wp,&
      & 1.84916514676949E-06_wp, 4.29209845613365E-08_wp, 1.90449923521437E-06_wp,&
      & 2.45475221461323E-07_wp, 1.38067478106305E-07_wp, 1.68098445714244E-07_wp,&
      & 2.05494837413711E-07_wp, 1.85334643734401E-08_wp, 9.43284242281929E-08_wp,&
      & 1.29699193497879E-08_wp, 3.25528435220703E-08_wp, 1.95567012770084E-08_wp,&
      & 1.06668704158924E-07_wp, 1.74925882038971E-06_wp, 0.00000000000000E+00_wp,&
      & 4.82884976492883E-07_wp, 4.84259784592653E-07_wp, 4.08025949077726E-07_wp,&
      & 3.62013202229254E-07_wp, 2.31260743908459E-07_wp, 4.79593449889770E-08_wp,&
      & 7.73150614742160E-07_wp, 2.99203674789581E-08_wp, 8.02715885043798E-08_wp,&
      & 3.18436323843386E-08_wp, 2.53784923205909E-06_wp, 2.37285579341720E-07_wp,&
      & 2.96490076287824E-06_wp, 4.98604264193123E-06_wp, 1.43570782883665E-06_wp,&
      & 4.82884976492883E-07_wp, 0.00000000000000E+00_wp, 9.99189184746263E-08_wp,&
      & 5.00301403714627E-06_wp, 3.59558582644546E-07_wp, 2.96004253338168E-07_wp,&
      & 1.19772242962055E-07_wp, 5.51087301186098E-07_wp, 2.91805748580956E-08_wp,&
      & 1.08474330641728E-08_wp, 1.71921133506486E-07_wp, 6.30130934390886E-09_wp,&
      & 2.80010329601741E-08_wp, 6.06527589798159E-09_wp, 2.02395032461751E-08_wp,&
      & 1.84916514676949E-06_wp, 4.84259784592653E-07_wp, 9.99189184746263E-08_wp,&
      & 0.00000000000000E+00_wp, 5.57241703386082E-08_wp, 1.77039850632016E-07_wp,&
      & 2.42350160026607E-07_wp, 1.56161285304105E-08_wp, 9.90233315017990E-08_wp,&
      & 6.44958519021423E-07_wp, 2.94907486282681E-08_wp, 3.27112069616806E-09_wp,&
      & 9.58364190274570E-08_wp, 7.70536801679902E-09_wp, 4.29330963668247E-06_wp,&
      & 8.46584342537644E-06_wp, 4.29209845613365E-08_wp, 4.08025949077726E-07_wp,&
      & 5.00301403714627E-06_wp, 5.57241703386082E-08_wp, 0.00000000000000E+00_wp,&
      & 4.61369497742773E-08_wp, 1.34479726009382E-08_wp, 1.65705208311747E-06_wp,&
      & 4.70699600335174E-07_wp, 8.40562670095430E-07_wp, 8.13521572303126E-08_wp,&
      & 1.80540186703872E-06_wp, 9.15484426760728E-09_wp, 1.83482506171271E-06_wp,&
      & 5.39530738035960E-09_wp, 2.53571466203672E-07_wp, 1.90449923521437E-06_wp,&
      & 3.62013202229254E-07_wp, 3.59558582644546E-07_wp, 1.77039850632016E-07_wp,&
      & 4.61369497742773E-08_wp, 0.00000000000000E+00_wp, 9.45767856532812E-08_wp,&
      & 2.15418363766496E-08_wp, 3.30901586776019E-07_wp, 1.43829354758462E-08_wp,&
      & 4.41257817307720E-08_wp, 2.17789465958412E-07_wp, 2.81737695822022E-08_wp,&
      & 2.27007490088747E-07_wp, 4.92625255487705E-09_wp, 2.80616983147280E-08_wp,&
      & 2.45475221461323E-07_wp, 2.31260743908459E-07_wp, 2.96004253338168E-07_wp,&
      & 2.42350160026607E-07_wp, 1.34479726009382E-08_wp, 9.45767856532812E-08_wp,&
      & 0.00000000000000E+00_wp, 3.25208271131549E-08_wp, 5.88381386437705E-08_wp,&
      & 3.95179687006130E-08_wp, 1.70043045502458E-08_wp, 5.61199401868333E-09_wp,&
      & 9.65080861967912E-08_wp, 3.96999055675365E-08_wp, 1.64928373932152E-07_wp,&
      & 7.01811077260374E-07_wp, 1.38067478106305E-07_wp, 4.79593449889770E-08_wp,&
      & 1.19772242962055E-07_wp, 1.56161285304105E-08_wp, 1.65705208311747E-06_wp,&
      & 2.15418363766496E-08_wp, 3.25208271131549E-08_wp, 0.00000000000000E+00_wp,&
      & 1.18414996245397E-06_wp, 1.50291385807284E-07_wp, 9.48553410038026E-08_wp,&
      & 1.05919193336971E-08_wp, 1.35290383223564E-07_wp, 3.18940263985877E-08_wp,&
      & 1.46884344099167E-06_wp, 4.68027254363654E-07_wp, 1.68098445714244E-07_wp,&
      & 7.73150614742160E-07_wp, 5.51087301186098E-07_wp, 9.90233315017990E-08_wp,&
      & 4.70699600335174E-07_wp, 3.30901586776019E-07_wp, 5.88381386437705E-08_wp,&
      & 1.18414996245397E-06_wp, 0.00000000000000E+00_wp, 4.07402456707778E-08_wp,&
      & 6.37407572831104E-09_wp, 2.86486826906670E-07_wp, 1.02635637280416E-07_wp,&
      & 3.98233483719282E-07_wp, 1.02143544654681E-07_wp, 8.41230027213861E-07_wp,&
      & 2.05494837413710E-07_wp, 2.99203674789581E-08_wp, 2.91805748580955E-08_wp,&
      & 6.44958519021421E-07_wp, 8.40562670095428E-07_wp, 1.43829354758462E-08_wp,&
      & 3.95179687006130E-08_wp, 1.50291385807284E-07_wp, 4.07402456707777E-08_wp,&
      & 0.00000000000000E+00_wp], shape(ref8))

   ! PBE-D4(sc)
   param = param_type(&
      & s6 = 1.0_wp, s8 = 0.95948085_wp, a1 = 0.38574991_wp, a2 = 4.80688534_wp, &
      & a3 = 0.6_wp, a4 = 0.6_wp)

   call get_structure(mol, "MB16-43", "01")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%screened, -1)
   if (allocated(error)) return
   call test_damping_2b_gen(error, mol, d4, damp, param, ref6, ref8)

end subroutine test_damp_screened_2b_mb01

subroutine test_damp_screened_2b_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   real(wp), parameter :: ref6(16, 16) = reshape([&
      & 0.00000000000000E+00_wp, 4.43960420211676E-05_wp, 6.41345687794813E-06_wp,&
      & 5.43473835735561E-05_wp, 3.07735482891268E-06_wp, 1.76208295472050E-04_wp,&
      & 1.65524465375002E-04_wp, 1.52809249035624E-06_wp, 1.27989503910954E-05_wp,&
      & 2.81777903592627E-05_wp, 1.07952212351110E-05_wp, 1.39068520804979E-05_wp,&
      & 1.80502334861211E-06_wp, 3.65124041443327E-06_wp, 1.77278270694809E-04_wp,&
      & 3.81404195252165E-06_wp, 4.43960420211676E-05_wp, 0.00000000000000E+00_wp,&
      & 3.22587114414972E-06_wp, 1.35424832633020E-05_wp, 3.08967516925818E-06_wp,&
      & 2.11341568295526E-05_wp, 4.75331651220247E-05_wp, 4.16806906846245E-05_wp,&
      & 7.42022184995262E-06_wp, 1.03881822200816E-05_wp, 3.37080308468720E-06_wp,&
      & 3.27649702770743E-06_wp, 1.58581444678564E-05_wp, 3.68250830960063E-06_wp,&
      & 4.40441651007290E-05_wp, 5.31255408039084E-06_wp, 6.41345687794813E-06_wp,&
      & 3.22587114414972E-06_wp, 0.00000000000000E+00_wp, 2.72225619572626E-06_wp,&
      & 2.30698218889007E-06_wp, 3.34830646392932E-06_wp, 2.11651381832761E-05_wp,&
      & 1.04239540141317E-05_wp, 1.73947962166423E-06_wp, 3.90827402667924E-06_wp,&
      & 2.09809452987787E-05_wp, 5.95495320361910E-06_wp, 6.81598110793212E-06_wp,&
      & 2.25257526850053E-05_wp, 4.08542230443641E-05_wp, 1.23064663678396E-05_wp,&
      & 5.43473835735561E-05_wp, 1.35424832633020E-05_wp, 2.72225619572626E-06_wp,&
      & 0.00000000000000E+00_wp, 1.26111425541010E-05_wp, 1.94600300333417E-05_wp,&
      & 2.93711610809132E-06_wp, 7.02683804956986E-07_wp, 5.33587257384541E-06_wp,&
      & 1.97364759832113E-05_wp, 1.80487463224943E-05_wp, 5.52680809312983E-06_wp,&
      & 2.37142882973005E-05_wp, 1.13228562245393E-06_wp, 1.87699871395873E-05_wp,&
      & 1.48574156690666E-05_wp, 3.07735482891268E-06_wp, 3.08967516925818E-06_wp,&
      & 2.30698218889007E-06_wp, 1.26111425541010E-05_wp, 0.00000000000000E+00_wp,&
      & 7.26041977240245E-06_wp, 1.89337338088248E-06_wp, 1.92297909671841E-05_wp,&
      & 1.29922099711315E-06_wp, 2.50009487016315E-05_wp, 6.29445817878738E-06_wp,&
      & 9.29317435867104E-07_wp, 4.67434047786745E-06_wp, 3.84819780622928E-06_wp,&
      & 1.04095371099382E-05_wp, 1.39999293649882E-06_wp, 1.76208295472050E-04_wp,&
      & 2.11341568295526E-05_wp, 3.34830646392932E-06_wp, 1.94600300333417E-05_wp,&
      & 7.26041977240245E-06_wp, 0.00000000000000E+00_wp, 9.37452032635273E-06_wp,&
      & 8.11799611675010E-07_wp, 7.85416625069978E-06_wp, 1.37034156254103E-04_wp,&
      & 5.75597121686933E-05_wp, 3.45448181164043E-05_wp, 7.31779402304623E-06_wp,&
      & 1.78003687226040E-06_wp, 6.08824043136650E-05_wp, 1.05841336125743E-05_wp,&
      & 1.65524465375002E-04_wp, 4.75331651220247E-05_wp, 2.11651381832761E-05_wp,&
      & 2.93711610809132E-06_wp, 1.89337338088248E-06_wp, 9.37452032635273E-06_wp,&
      & 0.00000000000000E+00_wp, 6.17188068559074E-06_wp, 2.33057007328833E-05_wp,&
      & 2.56099234625799E-06_wp, 1.04972959073836E-05_wp, 2.21318205962249E-06_wp,&
      & 6.64570107635276E-07_wp, 1.13880637761096E-05_wp, 1.71989313556295E-04_wp,&
      & 1.88839438152465E-06_wp, 1.52809249035624E-06_wp, 4.16806906846245E-05_wp,&
      & 1.04239540141317E-05_wp, 7.02683804956986E-07_wp, 1.92297909671841E-05_wp,&
      & 8.11799611675010E-07_wp, 6.17188068559074E-06_wp, 0.00000000000000E+00_wp,&
      & 2.75280879215100E-05_wp, 9.15371218863842E-07_wp, 7.27962499335901E-06_wp,&
      & 1.75927064477405E-06_wp, 2.47319846843358E-06_wp, 7.18128644275564E-06_wp,&
      & 6.48600270412218E-06_wp, 3.02449891044266E-06_wp, 1.27989503910954E-05_wp,&
      & 7.42022184995262E-06_wp, 1.73947962166423E-06_wp, 5.33587257384541E-06_wp,&
      & 1.29922099711315E-06_wp, 7.85416625069978E-06_wp, 2.33057007328833E-05_wp,&
      & 2.75280879215100E-05_wp, 0.00000000000000E+00_wp, 1.29571991633019E-05_wp,&
      & 6.01035798501628E-06_wp, 3.27341836266361E-06_wp, 2.90536771276773E-06_wp,&
      & 6.76338635126860E-06_wp, 9.07718578357070E-06_wp, 1.41263643948278E-06_wp,&
      & 2.81777903592627E-05_wp, 1.03881822200816E-05_wp, 3.90827402667924E-06_wp,&
      & 1.97364759832113E-05_wp, 2.50009487016315E-05_wp, 1.37034156254103E-04_wp,&
      & 2.56099234625799E-06_wp, 9.15371218863842E-07_wp, 1.29571991633019E-05_wp,&
      & 0.00000000000000E+00_wp, 1.51367974055984E-05_wp, 3.44740740284382E-05_wp,&
      & 2.55235222921964E-05_wp, 2.39847999208443E-06_wp, 3.15086959214516E-05_wp,&
      & 5.34113327585157E-05_wp, 1.07952212351110E-05_wp, 3.37080308468720E-06_wp,&
      & 2.09809452987787E-05_wp, 1.80487463224943E-05_wp, 6.29445817878738E-06_wp,&
      & 5.75597121686933E-05_wp, 1.04972959073836E-05_wp, 7.27962499335901E-06_wp,&
      & 6.01035798501628E-06_wp, 1.51367974055984E-05_wp, 0.00000000000000E+00_wp,&
      & 9.27773272749164E-06_wp, 2.28297744645850E-06_wp, 1.29346623625015E-05_wp,&
      & 1.36397933510193E-05_wp, 6.32588326503760E-06_wp, 1.39068520804979E-05_wp,&
      & 3.27649702770743E-06_wp, 5.95495320361910E-06_wp, 5.52680809312983E-06_wp,&
      & 9.29317435867104E-07_wp, 3.45448181164043E-05_wp, 2.21318205962249E-06_wp,&
      & 1.75927064477405E-06_wp, 3.27341836266361E-06_wp, 3.44740740284382E-05_wp,&
      & 9.27773272749164E-06_wp, 0.00000000000000E+00_wp, 5.70978262207076E-06_wp,&
      & 6.07095222292657E-07_wp, 5.88947764915525E-06_wp, 5.16801654325163E-06_wp,&
      & 1.80502334861211E-06_wp, 1.58581444678564E-05_wp, 6.81598110793212E-06_wp,&
      & 2.37142882973005E-05_wp, 4.67434047786745E-06_wp, 7.31779402304623E-06_wp,&
      & 6.64570107635276E-07_wp, 2.47319846843358E-06_wp, 2.90536771276773E-06_wp,&
      & 2.55235222921964E-05_wp, 2.28297744645850E-06_wp, 5.70978262207076E-06_wp,&
      & 0.00000000000000E+00_wp, 5.25042630057703E-07_wp, 2.06301340198709E-06_wp,&
      & 8.36789382074104E-06_wp, 3.65124041443327E-06_wp, 3.68250830960063E-06_wp,&
      & 2.25257526850053E-05_wp, 1.13228562245393E-06_wp, 3.84819780622928E-06_wp,&
      & 1.78003687226040E-06_wp, 1.13880637761096E-05_wp, 7.18128644275564E-06_wp,&
      & 6.76338635126860E-06_wp, 2.39847999208443E-06_wp, 1.29346623625015E-05_wp,&
      & 6.07095222292657E-07_wp, 5.25042630057703E-07_wp, 0.00000000000000E+00_wp,&
      & 1.08799926141897E-04_wp, 2.76889210706195E-05_wp, 1.77278270694809E-04_wp,&
      & 4.40441651007290E-05_wp, 4.08542230443641E-05_wp, 1.87699871395873E-05_wp,&
      & 1.04095371099382E-05_wp, 6.08824043136650E-05_wp, 1.71989313556295E-04_wp,&
      & 6.48600270412218E-06_wp, 9.07718578357070E-06_wp, 3.15086959214516E-05_wp,&
      & 1.36397933510193E-05_wp, 5.88947764915525E-06_wp, 2.06301340198709E-06_wp,&
      & 1.08799926141897E-04_wp, 0.00000000000000E+00_wp, 3.39594405563334E-05_wp,&
      & 3.81404195252165E-06_wp, 5.31255408039084E-06_wp, 1.23064663678396E-05_wp,&
      & 1.48574156690666E-05_wp, 1.39999293649882E-06_wp, 1.05841336125743E-05_wp,&
      & 1.88839438152465E-06_wp, 3.02449891044266E-06_wp, 1.41263643948278E-06_wp,&
      & 5.34113327585157E-05_wp, 6.32588326503760E-06_wp, 5.16801654325163E-06_wp,&
      & 8.36789382074104E-06_wp, 2.76889210706195E-05_wp, 3.39594405563334E-05_wp,&
      & 0.00000000000000E+00_wp], shape(ref6))

   real(wp), parameter :: ref8(16, 16) = reshape([&
      & 0.00000000000000E+00_wp, 1.50832634355079E-06_wp, 1.14329740188187E-07_wp,&
      & 1.97518524340416E-06_wp, 4.29476835166620E-08_wp, 9.47846719570597E-06_wp,&
      & 8.72005524387999E-06_wp, 1.68876479350852E-08_wp, 2.87255786243287E-07_wp,&
      & 8.22707688782856E-07_wp, 2.28917173018236E-07_wp, 3.20879095642001E-07_wp,&
      & 2.10869326000276E-08_wp, 5.39457063886465E-08_wp, 9.55528518358186E-06_wp,&
      & 5.71764164336513E-08_wp, 1.50832634355079E-06_wp, 0.00000000000000E+00_wp,&
      & 4.57332781016569E-08_wp, 3.09718671616739E-07_wp, 4.31770938859225E-08_wp,&
      & 5.60638449417400E-07_wp, 1.65208335794686E-06_wp, 1.38659465100004E-06_wp,&
      & 1.38864797300686E-07_wp, 2.17481535139261E-07_wp, 4.84931975301740E-08_wp,&
      & 4.66927388420669E-08_wp, 3.82272277927690E-07_wp, 5.45625457506452E-08_wp,&
      & 1.49240772067572E-06_wp, 8.89421172047510E-08_wp, 1.14329740188187E-07_wp,&
      & 4.57332781016569E-08_wp, 0.00000000000000E+00_wp, 3.64704933109258E-08_wp,&
      & 2.92479411318643E-08_wp, 4.80621557211404E-08_wp, 5.61734531656741E-07_wp,&
      & 2.18480640471488E-07_wp, 2.00722213793549E-08_wp, 5.90676446000717E-08_wp,&
      & 5.55225895434974E-07_wp, 1.03563673889473E-07_wp, 1.23995944738482E-07_wp,&
      & 6.10391850570317E-07_wp, 1.35005745704094E-06_wp, 2.72613569988802E-07_wp,&
      & 1.97518524340416E-06_wp, 3.09718671616739E-07_wp, 3.64704933109258E-08_wp,&
      & 0.00000000000000E+00_wp, 2.81649439353691E-07_wp, 5.02220365133747E-07_wp,&
      & 4.03581389583065E-08_wp, 5.99399307581999E-09_wp, 8.94630247023868E-08_wp,&
      & 5.11755428845807E-07_wp, 4.54254284690135E-07_wp, 9.37566665014987E-08_wp,&
      & 6.53706954072466E-07_wp, 1.13234204885091E-08_wp, 4.78617173781924E-07_wp,&
      & 3.50451031847996E-07_wp, 4.29476835166620E-08_wp, 4.31770938859225E-08_wp,&
      & 2.92479411318643E-08_wp, 2.81649439353691E-07_wp, 0.00000000000000E+00_wp,&
      & 1.34891717908075E-07_wp, 2.24742220386185E-08_wp, 4.94313413755495E-07_wp,&
      & 1.36023418590515E-08_wp, 7.01420244702676E-07_wp, 1.11510076658348E-07_wp,&
      & 8.70138580851840E-09_wp, 7.49888629864903E-08_wp, 5.78601419483570E-08_wp,&
      & 2.18077838978822E-07_wp, 1.50269475494392E-08_wp, 9.47846719570597E-06_wp,&
      & 5.60638449417400E-07_wp, 4.80621557211404E-08_wp, 5.02220365133747E-07_wp,&
      & 1.34891717908075E-07_wp, 0.00000000000000E+00_wp, 1.89656785158349E-07_wp,&
      & 7.26610099746616E-09_wp, 1.49796993739532E-07_wp, 6.77862201906525E-06_wp,&
      & 2.13236295616368E-06_wp, 1.07947732721414E-06_wp, 1.36314865464553E-07_wp,&
      & 2.06986319642856E-08_wp, 2.29804592935127E-06_wp, 2.22968440408455E-07_wp,&
      & 8.72005524387999E-06_wp, 1.65208335794686E-06_wp, 5.61734531656741E-07_wp,&
      & 4.03581389583065E-08_wp, 2.24742220386185E-08_wp, 1.89656785158349E-07_wp,&
      & 0.00000000000000E+00_wp, 1.08624137270821E-07_wp, 6.38732770475943E-07_wp,&
      & 3.36186821605443E-08_wp, 2.20532650650447E-07_wp, 2.76731862470966E-08_wp,&
      & 5.56447311304009E-09_wp, 2.45830726937907E-07_wp, 9.17708868688325E-06_wp,&
      & 2.23954560533109E-08_wp, 1.68876479350852E-08_wp, 1.38659465100004E-06_wp,&
      & 2.18480640471488E-07_wp, 5.99399307581999E-09_wp, 4.94313413755495E-07_wp,&
      & 7.26610099746616E-09_wp, 1.08624137270821E-07_wp, 0.00000000000000E+00_wp,&
      & 7.97512881963463E-07_wp, 8.52771443433429E-09_wp, 1.35367680598225E-07_wp,&
      & 2.03772943717781E-08_wp, 3.20908782642877E-08_wp, 1.32934991192183E-07_wp,&
      & 1.16057304180595E-07_wp, 4.19669641971747E-08_wp, 2.87255786243287E-07_wp,&
      & 1.38864797300686E-07_wp, 2.00722213793549E-08_wp, 8.94630247023868E-08_wp,&
      & 1.36023418590515E-08_wp, 1.49796993739532E-07_wp, 6.38732770475943E-07_wp,&
      & 7.97512881963463E-07_wp, 0.00000000000000E+00_wp, 2.92001101868296E-07_wp,&
      & 1.04850401487980E-07_wp, 4.66342500413747E-08_wp, 3.97775276014124E-08_wp,&
      & 1.22721853766814E-07_wp, 1.81678949579366E-07_wp, 1.52081662348863E-08_wp,&
      & 8.22707688782856E-07_wp, 2.17481535139261E-07_wp, 5.90676446000717E-08_wp,&
      & 5.11755428845807E-07_wp, 7.01420244702676E-07_wp, 6.77862201906525E-06_wp,&
      & 3.36186821605443E-08_wp, 8.52771443433429E-09_wp, 2.92001101868296E-07_wp,&
      & 0.00000000000000E+00_wp, 3.59265054390042E-07_wp, 1.07653079509855E-06_wp,&
      & 7.21036286492258E-07_wp, 3.08047563310233E-08_wp, 9.54868463008901E-07_wp,&
      & 1.92995654764448E-06_wp, 2.28917173018236E-07_wp, 4.84931975301740E-08_wp,&
      & 5.55225895434974E-07_wp, 4.54254284690135E-07_wp, 1.11510076658348E-07_wp,&
      & 2.13236295616368E-06_wp, 2.20532650650447E-07_wp, 1.35367680598225E-07_wp,&
      & 1.04850401487980E-07_wp, 3.59265054390042E-07_wp, 0.00000000000000E+00_wp,&
      & 1.87050462988074E-07_wp, 2.88428700109008E-08_wp, 2.91324117883242E-07_wp,&
      & 3.12689548467284E-07_wp, 1.12252979117870E-07_wp, 3.20879095642001E-07_wp,&
      & 4.66927388420669E-08_wp, 1.03563673889473E-07_wp, 9.37566665014987E-08_wp,&
      & 8.70138580851840E-09_wp, 1.07947732721414E-06_wp, 2.76731862470966E-08_wp,&
      & 2.03772943717781E-08_wp, 4.66342500413747E-08_wp, 1.07653079509855E-06_wp,&
      & 1.87050462988074E-07_wp, 0.00000000000000E+00_wp, 9.79179720837000E-08_wp,&
      & 4.93225359442576E-09_wp, 1.02048199941860E-07_wp, 8.57303976466581E-08_wp,&
      & 2.10869326000276E-08_wp, 3.82272277927690E-07_wp, 1.23995944738482E-07_wp,&
      & 6.53706954072466E-07_wp, 7.49888629864903E-08_wp, 1.36314865464553E-07_wp,&
      & 5.56447311304009E-09_wp, 3.20908782642877E-08_wp, 3.97775276014124E-08_wp,&
      & 7.21036286492258E-07_wp, 2.88428700109008E-08_wp, 9.79179720837000E-08_wp,&
      & 0.00000000000000E+00_wp, 4.06408142130737E-09_wp, 2.51983658413843E-08_wp,&
      & 1.63001346629557E-07_wp, 5.39457063886465E-08_wp, 5.45625457506452E-08_wp,&
      & 6.10391850570317E-07_wp, 1.13234204885091E-08_wp, 5.78601419483570E-08_wp,&
      & 2.06986319642856E-08_wp, 2.45830726937907E-07_wp, 1.32934991192183E-07_wp,&
      & 1.22721853766814E-07_wp, 3.08047563310233E-08_wp, 2.91324117883242E-07_wp,&
      & 4.93225359442576E-09_wp, 4.06408142130737E-09_wp, 0.00000000000000E+00_wp,&
      & 4.98357630119030E-06_wp, 8.03731560456559E-07_wp, 9.55528518358186E-06_wp,&
      & 1.49240772067572E-06_wp, 1.35005745704094E-06_wp, 4.78617173781924E-07_wp,&
      & 2.18077838978822E-07_wp, 2.29804592935127E-06_wp, 9.17708868688325E-06_wp,&
      & 1.16057304180595E-07_wp, 1.81678949579366E-07_wp, 9.54868463008901E-07_wp,&
      & 3.12689548467284E-07_wp, 1.02048199941860E-07_wp, 2.51983658413843E-08_wp,&
      & 4.98357630119030E-06_wp, 0.00000000000000E+00_wp, 1.05515682893766E-06_wp,&
      & 5.71764164336513E-08_wp, 8.89421172047510E-08_wp, 2.72613569988802E-07_wp,&
      & 3.50451031847996E-07_wp, 1.50269475494392E-08_wp, 2.22968440408455E-07_wp,&
      & 2.23954560533109E-08_wp, 4.19669641971747E-08_wp, 1.52081662348863E-08_wp,&
      & 1.92995654764448E-06_wp, 1.12252979117870E-07_wp, 8.57303976466581E-08_wp,&
      & 1.63001346629557E-07_wp, 8.03731560456559E-07_wp, 1.05515682893766E-06_wp,&
      & 0.00000000000000E+00_wp], shape(ref8))

   ! PBE-D4(sc)
   param = param_type(&
      & s6 = 1.0_wp, s8 = 0.95948085_wp, a1 = 0.38574991_wp, a2 = 4.80688534_wp, &
      & a3 = 0.6_wp, a4 = 0.6_wp)

   call get_structure(mol, "MB16-43", "02")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%screened, -1)
   if (allocated(error)) return
   call test_damping_2b_gen(error, mol, d4, damp, param, ref6, ref8)

end subroutine test_damp_screened_2b_mb02

subroutine test_grad_screened_2b_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   ! PBE-D4(sc)
   param = param_type(&
      & s6 = 1.0_wp, s8 = 0.95948085_wp, a1 = 0.38574991_wp, a2 = 4.80688534_wp, &
      & a3 = 0.6_wp, a4 = 0.6_wp)

   call get_structure(mol, "MB16-43", "03")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%screened, -1)
   if (allocated(error)) return
   call test_damping_2b_numgrad(error, mol, d4, damp, param)

end subroutine test_grad_screened_2b_mb03

subroutine test_grad_screened_2b_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   ! PBE-D4(sc)
   param = param_type(&
      & s6 = 1.0_wp, s8 = 0.95948085_wp, a1 = 0.38574991_wp, a2 = 4.80688534_wp, &
      & a3 = 0.6_wp, a4 = 0.6_wp)

   call get_structure(mol, "MB16-43", "04")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%screened, -1)
   if (allocated(error)) return
   call test_damping_2b_numgrad(error, mol, d4, damp, param)

end subroutine test_grad_screened_2b_mb04


subroutine test_damp_zero_2b_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   real(wp), parameter :: ref6(16, 16) = reshape([&
      & 0.00000000000000E+00_wp, 6.88413960858574E-07_wp, 8.04183806038227E-10_wp,&
      & 8.43518801588689E-09_wp, 2.09424879273975E-07_wp, 8.82411451324456E-07_wp,&
      & 3.98496864081583E-08_wp, 4.36371289848057E-07_wp, 2.14518494955439E-09_wp,&
      & 9.29574558248900E-07_wp, 1.14180923555798E-06_wp, 8.53883088553101E-09_wp,&
      & 1.05463511200401E-09_wp, 4.07332114713338E-10_wp, 8.09086078865742E-08_wp,&
      & 1.51003831897163E-10_wp, 6.88413960858574E-07_wp, 0.00000000000000E+00_wp,&
      & 6.43848856464688E-07_wp, 1.33436770586981E-05_wp, 2.11817283674042E-07_wp,&
      & 6.64044850498431E-07_wp, 1.16333084408331E-05_wp, 5.48851965067575E-06_wp,&
      & 2.45307844764535E-06_wp, 8.66197349003322E-06_wp, 4.46162285126318E-07_wp,&
      & 1.31217743785190E-06_wp, 2.83785009638503E-09_wp, 6.68220759727773E-07_wp,&
      & 1.07659853971432E-06_wp, 6.13693024005148E-07_wp, 8.04183806038227E-10_wp,&
      & 6.43848856464688E-07_wp, 0.00000000000000E+00_wp, 3.54767667257219E-06_wp,&
      & 2.53399064953528E-06_wp, 3.87426187592109E-06_wp, 5.92422419231073E-06_wp,&
      & 1.25093775276718E-06_wp, 1.09323626089675E-06_wp, 7.29460216251989E-07_wp,&
      & 5.55241339125867E-06_wp, 9.51171703775904E-07_wp, 2.05586510099161E-06_wp,&
      & 9.68925308596002E-10_wp, 6.17265534491058E-06_wp, 2.93736641126022E-07_wp,&
      & 8.43518801588689E-09_wp, 1.33436770586981E-05_wp, 3.54767667257219E-06_wp,&
      & 0.00000000000000E+00_wp, 6.18770102249371E-07_wp, 1.78806888614815E-06_wp,&
      & 1.27093605301791E-05_wp, 2.49535209485431E-06_wp, 1.02784900314725E-05_wp,&
      & 2.23230330544886E-06_wp, 8.48315373167665E-07_wp, 1.27289217052490E-06_wp,&
      & 3.20642208487729E-09_wp, 2.82316118656716E-06_wp, 2.45597336375221E-06_wp,&
      & 4.28516469660418E-07_wp, 2.09424879273975E-07_wp, 2.11817283674042E-07_wp,&
      & 2.53399064953528E-06_wp, 6.18770102249371E-07_wp, 0.00000000000000E+00_wp,&
      & 1.93649622389691E-05_wp, 1.47821174861062E-06_wp, 1.70161779392116E-06_wp,&
      & 1.75693268797298E-06_wp, 7.08904037529611E-07_wp, 1.19702267060771E-05_wp,&
      & 6.47203580303879E-07_wp, 6.05078236021314E-07_wp, 2.86309273617930E-09_wp,&
      & 5.19204484542973E-06_wp, 5.08290063587853E-07_wp, 8.82411451324456E-07_wp,&
      & 6.64044850498431E-07_wp, 3.87426187592109E-06_wp, 1.78806888614815E-06_wp,&
      & 1.93649622389691E-05_wp, 0.00000000000000E+00_wp, 3.27720638843620E-06_wp,&
      & 6.00224790066316E-06_wp, 2.21360133151320E-06_wp, 1.75011232500897E-06_wp,&
      & 2.62277182123362E-06_wp, 6.35870529619291E-06_wp, 2.20329053439577E-06_wp,&
      & 4.64381968077780E-06_wp, 6.54025035819969E-09_wp, 9.23774740118466E-08_wp,&
      & 3.98496864081583E-08_wp, 1.16333084408331E-05_wp, 5.92422419231073E-06_wp,&
      & 1.27093605301791E-05_wp, 1.47821174861062E-06_wp, 3.27720638843620E-06_wp,&
      & 0.00000000000000E+00_wp, 1.63384075042767E-05_wp, 1.45330584141231E-05_wp,&
      & 3.87147601989630E-05_wp, 3.07442407392868E-06_wp, 1.18243245093900E-06_wp,&
      & 3.97958713368447E-09_wp, 5.93818739345938E-06_wp, 8.20788811013322E-06_wp,&
      & 6.39102137199099E-09_wp, 4.36371289848057E-07_wp, 5.48851965067575E-06_wp,&
      & 1.25093775276718E-06_wp, 2.49535209485431E-06_wp, 1.70161779392116E-06_wp,&
      & 6.00224790066316E-06_wp, 1.63384075042767E-05_wp, 0.00000000000000E+00_wp,&
      & 7.05465842100888E-06_wp, 5.66695680559431E-09_wp, 1.45028314421774E-05_wp,&
      & 1.42723353454195E-06_wp, 2.33415777586855E-06_wp, 2.71548133737837E-06_wp,&
      & 5.43531377659134E-06_wp, 3.84349823958801E-10_wp, 2.14518494955439E-09_wp,&
      & 2.45307844764535E-06_wp, 1.09323626089675E-06_wp, 1.02784900314725E-05_wp,&
      & 1.75693268797298E-06_wp, 2.21360133151320E-06_wp, 1.45330584141231E-05_wp,&
      & 7.05465842100888E-06_wp, 0.00000000000000E+00_wp, 5.69741872504042E-06_wp,&
      & 2.13315145420947E-06_wp, 1.07007466708843E-06_wp, 1.48531873770816E-06_wp,&
      & 2.11957055102939E-09_wp, 2.09271892515983E-08_wp, 4.41097584633871E-10_wp,&
      & 9.29574558248900E-07_wp, 8.66197349003322E-06_wp, 7.29460216251989E-07_wp,&
      & 2.23230330544886E-06_wp, 7.08904037529611E-07_wp, 1.75011232500897E-06_wp,&
      & 3.87147601989630E-05_wp, 5.66695680559431E-09_wp, 5.69741872504042E-06_wp,&
      & 0.00000000000000E+00_wp, 3.73828312471502E-06_wp, 6.16779078100834E-06_wp,&
      & 6.81039589022605E-06_wp, 1.43322954631533E-06_wp, 5.66059229851410E-06_wp,&
      & 3.28960821435868E-08_wp, 1.14180923555798E-06_wp, 4.46162285126318E-07_wp,&
      & 5.55241339125867E-06_wp, 8.48315373167665E-07_wp, 1.19702267060771E-05_wp,&
      & 2.62277182123362E-06_wp, 3.07442407392868E-06_wp, 1.45028314421774E-05_wp,&
      & 2.13315145420947E-06_wp, 3.73828312471502E-06_wp, 0.00000000000000E+00_wp,&
      & 3.11574043359528E-06_wp, 1.28281823813947E-06_wp, 1.83600547334721E-06_wp,&
      & 6.66344268554393E-09_wp, 7.85587366521644E-08_wp, 8.53883088553101E-09_wp,&
      & 1.31217743785190E-06_wp, 9.51171703775904E-07_wp, 1.27289217052490E-06_wp,&
      & 6.47203580303879E-07_wp, 6.35870529619291E-06_wp, 1.18243245093900E-06_wp,&
      & 1.42723353454195E-06_wp, 1.07007466708843E-06_wp, 6.16779078100834E-06_wp,&
      & 3.11574043359528E-06_wp, 0.00000000000000E+00_wp, 2.47215437865693E-09_wp,&
      & 1.06403813982370E-06_wp, 1.15832501506677E-06_wp, 2.90461298013617E-10_wp,&
      & 1.05463511200401E-09_wp, 2.83785009638503E-09_wp, 2.05586510099161E-06_wp,&
      & 3.20642208487729E-09_wp, 6.05078236021314E-07_wp, 2.20329053439577E-06_wp,&
      & 3.97958713368447E-09_wp, 2.33415777586855E-06_wp, 1.48531873770816E-06_wp,&
      & 6.81039589022605E-06_wp, 1.28281823813947E-06_wp, 2.47215437865693E-09_wp,&
      & 0.00000000000000E+00_wp, 1.10202356104784E-06_wp, 2.68497141215432E-06_wp,&
      & 9.85563379409269E-10_wp, 4.07332114713338E-10_wp, 6.68220759727773E-07_wp,&
      & 9.68925308596002E-10_wp, 2.82316118656716E-06_wp, 2.86309273617930E-09_wp,&
      & 4.64381968077780E-06_wp, 5.93818739345938E-06_wp, 2.71548133737837E-06_wp,&
      & 2.11957055102939E-09_wp, 1.43322954631533E-06_wp, 1.83600547334721E-06_wp,&
      & 1.06403813982370E-06_wp, 1.10202356104784E-06_wp, 0.00000000000000E+00_wp,&
      & 1.86383755719533E-07_wp, 5.09571309160253E-09_wp, 8.09086078865742E-08_wp,&
      & 1.07659853971432E-06_wp, 6.17265534491058E-06_wp, 2.45597336375221E-06_wp,&
      & 5.19204484542973E-06_wp, 6.54025035819969E-09_wp, 8.20788811013322E-06_wp,&
      & 5.43531377659134E-06_wp, 2.09271892515983E-08_wp, 5.66059229851410E-06_wp,&
      & 6.66344268554393E-09_wp, 1.15832501506677E-06_wp, 2.68497141215432E-06_wp,&
      & 1.86383755719533E-07_wp, 0.00000000000000E+00_wp, 8.30063071743697E-10_wp,&
      & 1.51003831897163E-10_wp, 6.13693024005148E-07_wp, 2.93736641126022E-07_wp,&
      & 4.28516469660418E-07_wp, 5.08290063587853E-07_wp, 9.23774740118466E-08_wp,&
      & 6.39102137199099E-09_wp, 3.84349823958801E-10_wp, 4.41097584633871E-10_wp,&
      & 3.28960821435868E-08_wp, 7.85587366521644E-08_wp, 2.90461298013617E-10_wp,&
      & 9.85563379409269E-10_wp, 5.09571309160253E-09_wp, 8.30063071743697E-10_wp,&
      & 0.00000000000000E+00_wp], shape(ref6))

   real(wp), parameter :: ref8(16, 16) = reshape([&
      & 0.00000000000000E+00_wp, 4.90402033666035E-08_wp, 1.77139492327965E-10_wp,&
      & 2.39988987970123E-09_wp, 2.94736969326410E-08_wp, 3.84225697386835E-08_wp,&
      & 1.12489427597280E-08_wp, 1.26056307128132E-08_wp, 4.51884863897259E-10_wp,&
      & 8.12445426848118E-09_wp, 2.14493777713098E-08_wp, 1.28790433525962E-09_wp,&
      & 1.65293527161265E-10_wp, 6.38520132825057E-11_wp, 1.40154933328535E-08_wp,&
      & 1.52623531381544E-11_wp, 4.90402033666035E-08_wp, 0.00000000000000E+00_wp,&
      & 4.01445177110866E-09_wp, 6.89054296688666E-06_wp, 9.11634806317231E-10_wp,&
      & 4.18287536383341E-09_wp, 6.76002232886062E-06_wp, 7.09696915430405E-08_wp,&
      & 2.39612334853694E-08_wp, 1.29373088993629E-07_wp, 2.46148640437487E-09_wp,&
      & 5.03053354379216E-07_wp, 1.45954098147194E-09_wp, 4.22292421128987E-09_wp,&
      & 7.97028672658607E-09_wp, 1.09753433980225E-07_wp, 1.77139492327965E-10_wp,&
      & 4.01445177110866E-09_wp, 0.00000000000000E+00_wp, 3.92781249161240E-08_wp,&
      & 1.12021753596273E-06_wp, 4.42267643251914E-08_wp, 7.88179808611140E-08_wp,&
      & 9.75944383227772E-09_wp, 5.10261731435047E-07_wp, 4.74166959832738E-09_wp,&
      & 7.21039717355144E-08_wp, 6.88704042118722E-09_wp, 2.11548974304374E-08_wp,&
      & 3.85686475293189E-10_wp, 1.01393228904868E-07_wp, 3.93962598676176E-08_wp,&
      & 2.39988987970123E-09_wp, 6.89054296688666E-06_wp, 3.92781249161240E-08_wp,&
      & 0.00000000000000E+00_wp, 3.80713887066749E-09_wp, 1.56721468896394E-08_wp,&
      & 6.86042323044523E-06_wp, 2.44951085059198E-08_wp, 1.78420302277857E-07_wp,&
      & 2.10704565842232E-08_wp, 5.79821340528338E-09_wp, 4.95088671808398E-07_wp,&
      & 1.64910160080065E-09_wp, 2.98575052772871E-08_wp, 2.39991507883536E-08_wp,&
      & 1.01151179374868E-07_wp, 2.94736969326410E-08_wp, 9.11634806317231E-10_wp,&
      & 1.12021753596273E-06_wp, 3.80713887066749E-09_wp, 0.00000000000000E+00_wp,&
      & 4.63542088161735E-07_wp, 1.21625199168623E-08_wp, 1.47157068519785E-08_wp,&
      & 8.28379953083375E-07_wp, 4.56405907040336E-09_wp, 3.17425379498889E-06_wp,&
      & 4.05971618094543E-09_wp, 3.70683038934871E-09_wp, 1.23766772047411E-09_wp,&
      & 9.20865537879056E-07_wp, 5.16681081995315E-08_wp, 3.84225697386835E-08_wp,&
      & 4.18287536383341E-09_wp, 4.42267643251914E-08_wp, 1.56721468896394E-08_wp,&
      & 4.63542088161735E-07_wp, 0.00000000000000E+00_wp, 3.51722345779119E-08_wp,&
      & 8.02515181499602E-08_wp, 1.35581331759839E-06_wp, 1.52300107445634E-08_wp,&
      & 2.37986835239811E-06_wp, 1.86679393968125E-07_wp, 2.11103532295114E-08_wp,&
      & 4.65976619208918E-07_wp, 4.52069406934493E-09_wp, 3.00969808665323E-08_wp,&
      & 1.12489427597280E-08_wp, 6.76002232886062E-06_wp, 7.88179808611140E-08_wp,&
      & 6.86042323044523E-06_wp, 1.21625199168623E-08_wp, 3.51722345779119E-08_wp,&
      & 0.00000000000000E+00_wp, 1.32542516372584E-06_wp, 1.07489403199117E-06_wp,&
      & 1.43835151087666E-06_wp, 3.22974596214541E-08_wp, 4.74951117161584E-07_wp,&
      & 2.04674828564456E-09_wp, 1.03287138870540E-07_wp, 1.26423795645105E-07_wp,&
      & 2.11925206334422E-09_wp, 1.26056307128132E-08_wp, 7.09696915430405E-08_wp,&
      & 9.75944383227772E-09_wp, 2.44951085059198E-08_wp, 1.47157068519785E-08_wp,&
      & 8.02515181499602E-08_wp, 1.32542516372584E-06_wp, 0.00000000000000E+00_wp,&
      & 3.48678836808107E-07_wp, 4.09548210649839E-09_wp, 3.07071067264859E-07_wp,&
      & 1.93383140810334E-07_wp, 1.53217490580092E-07_wp, 3.58669539635157E-08_wp,&
      & 5.23818063059040E-07_wp, 9.86482672127482E-11_wp, 4.51884863897259E-10_wp,&
      & 2.39612334853694E-08_wp, 5.10261731435047E-07_wp, 1.78420302277857E-07_wp,&
      & 8.28379953083375E-07_wp, 1.35581331759839E-06_wp, 1.07489403199117E-06_wp,&
      & 3.48678836808107E-07_wp, 0.00000000000000E+00_wp, 7.51682112823424E-08_wp,&
      & 1.31741780151134E-06_wp, 1.70339635718447E-07_wp, 1.68224264035909E-07_wp,&
      & 8.06949770485052E-10_wp, 1.07070076212749E-08_wp, 1.08281148697716E-10_wp,&
      & 8.12445426848118E-09_wp, 1.29373088993629E-07_wp, 4.74166959832738E-09_wp,&
      & 2.10704565842232E-08_wp, 4.56405907040336E-09_wp, 1.52300107445634E-08_wp,&
      & 1.43835151087666E-06_wp, 4.09548210649839E-09_wp, 7.51682112823424E-08_wp,&
      & 0.00000000000000E+00_wp, 4.19315685770362E-08_wp, 1.31728286336873E-07_wp,&
      & 1.79392910543320E-07_wp, 1.17499801440150E-08_wp, 7.44947463007276E-08_wp,&
      & 1.08851755216956E-08_wp, 2.14493777713098E-08_wp, 2.46148640437487E-09_wp,&
      & 7.21039717355144E-08_wp, 5.79821340528338E-09_wp, 3.17425379498889E-06_wp,&
      & 2.37986835239811E-06_wp, 3.22974596214541E-08_wp, 3.07071067264859E-07_wp,&
      & 1.31741780151134E-06_wp, 4.19315685770362E-08_wp, 0.00000000000000E+00_wp,&
      & 3.46870362672762E-08_wp, 1.01188267116996E-08_wp, 6.13537889065857E-07_wp,&
      & 4.60584582424786E-09_wp, 2.57210893249569E-08_wp, 1.28790433525962E-09_wp,&
      & 5.03053354379216E-07_wp, 6.88704042118722E-09_wp, 4.95088671808398E-07_wp,&
      & 4.05971618094543E-09_wp, 1.86679393968125E-07_wp, 4.74951117161584E-07_wp,&
      & 1.93383140810334E-07_wp, 1.70339635718447E-07_wp, 1.31728286336873E-07_wp,&
      & 3.46870362672762E-08_wp, 0.00000000000000E+00_wp, 6.84354326087749E-10_wp,&
      & 1.58444154697615E-08_wp, 1.66574803970585E-07_wp, 5.18471189323666E-11_wp,&
      & 1.65293527161265E-10_wp, 1.45954098147194E-09_wp, 2.11548974304374E-08_wp,&
      & 1.64910160080065E-09_wp, 3.70683038934871E-09_wp, 2.11103532295114E-08_wp,&
      & 2.04674828564456E-09_wp, 1.53217490580092E-07_wp, 1.68224264035909E-07_wp,&
      & 1.79392910543320E-07_wp, 1.01188267116996E-08_wp, 6.84354326087749E-10_wp,&
      & 0.00000000000000E+00_wp, 2.34480342171035E-08_wp, 4.36962395643603E-08_wp,&
      & 1.80007676676439E-10_wp, 6.38520132825057E-11_wp, 4.22292421128987E-09_wp,&
      & 3.85686475293189E-10_wp, 2.98575052772871E-08_wp, 1.23766772047411E-09_wp,&
      & 4.65976619208918E-07_wp, 1.03287138870540E-07_wp, 3.58669539635157E-08_wp,&
      & 8.06949770485052E-10_wp, 1.17499801440150E-08_wp, 6.13537889065857E-07_wp,&
      & 1.58444154697615E-08_wp, 2.34480342171035E-08_wp, 0.00000000000000E+00_wp,&
      & 6.87589070399636E-08_wp, 9.29011951398606E-10_wp, 1.40154933328535E-08_wp,&
      & 7.97028672658607E-09_wp, 1.01393228904868E-07_wp, 2.39991507883536E-08_wp,&
      & 9.20865537879056E-07_wp, 4.52069406934493E-09_wp, 1.26423795645105E-07_wp,&
      & 5.23818063059040E-07_wp, 1.07070076212749E-08_wp, 7.44947463007276E-08_wp,&
      & 4.60584582424786E-09_wp, 1.66574803970585E-07_wp, 4.36962395643603E-08_wp,&
      & 6.87589070399636E-08_wp, 0.00000000000000E+00_wp, 2.03763735310248E-10_wp,&
      & 1.52623531381544E-11_wp, 1.09753433980225E-07_wp, 3.93962598676176E-08_wp,&
      & 1.01151179374868E-07_wp, 5.16681081995315E-08_wp, 3.00969808665323E-08_wp,&
      & 2.11925206334422E-09_wp, 9.86482672127482E-11_wp, 1.08281148697716E-10_wp,&
      & 1.08851755216956E-08_wp, 2.57210893249569E-08_wp, 5.18471189323666E-11_wp,&
      & 1.80007676676439E-10_wp, 9.29011951398606E-10_wp, 2.03763735310248E-10_wp,&
      & 0.00000000000000E+00_wp], shape(ref8))

   ! PBE-D3(0) (using D4 r4/r2 radii instead of D3 vdW radii)
   param = param_type(&
      & s6 = 1.0_wp, s8 = 0.722_wp, a1 = 1.0_wp, a2 = 0.0_wp, &
      & rs6 = 1.217_wp, rs8 = 1.0_wp, alp = 14.0_wp)

   call get_structure(mol, "MB16-43", "01")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%zero, -1)
   if (allocated(error)) return
   call test_damping_2b_gen(error, mol, d4, damp, param, ref6, ref8)

end subroutine test_damp_zero_2b_mb01

subroutine test_damp_zero_2b_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   real(wp), parameter :: ref6(16, 16) = reshape([&
      & 0.00000000000000E+00_wp, 9.01705849264917E-07_wp, 5.45561439306869E-06_wp,&
      & 1.51946975105081E-05_wp, 1.99975128279883E-06_wp, 4.34489747120837E-06_wp,&
      & 1.33963323128116E-05_wp, 1.52794708791284E-06_wp, 1.60432360805069E-06_wp,&
      & 2.60958239754062E-05_wp, 3.28428291428594E-09_wp, 1.22038974098681E-06_wp,&
      & 1.80416885064417E-06_wp, 3.64859359329591E-06_wp, 8.97509255934878E-06_wp,&
      & 3.44312903091752E-06_wp, 9.01705849264917E-07_wp, 0.00000000000000E+00_wp,&
      & 9.30310315072224E-10_wp, 9.72801144012400E-07_wp, 1.06419841772418E-09_wp,&
      & 2.97656362010527E-06_wp, 7.46242260907044E-07_wp, 1.04617002193055E-06_wp,&
      & 4.11497988679676E-09_wp, 4.90405805872025E-06_wp, 1.02381089813213E-09_wp,&
      & 1.17638119754792E-09_wp, 1.37693495675094E-06_wp, 3.34999221798228E-06_wp,&
      & 9.19792252998876E-07_wp, 1.61806219757397E-07_wp, 5.45561439306869E-06_wp,&
      & 9.30310315072224E-10_wp, 0.00000000000000E+00_wp, 2.38147221678861E-06_wp,&
      & 7.01988052385574E-10_wp, 3.22393893136350E-06_wp, 5.40298210888355E-06_wp,&
      & 2.86505721882952E-09_wp, 3.40154362497530E-10_wp, 3.70337380076441E-06_wp,&
      & 4.19271847501332E-08_wp, 6.18155501580136E-08_wp, 4.04805092715763E-06_wp,&
      & 5.12707170946965E-06_wp, 2.43353063964647E-06_wp, 9.58218767201640E-08_wp,&
      & 1.51946975105081E-05_wp, 9.72801144012400E-07_wp, 2.38147221678861E-06_wp,&
      & 0.00000000000000E+00_wp, 1.24197394522869E-07_wp, 6.45393876041118E-09_wp,&
      & 2.92942083080631E-06_wp, 7.02618201266191E-07_wp, 8.40186543771942E-07_wp,&
      & 6.88715357186242E-09_wp, 1.34014614738065E-06_wp, 2.35505184164687E-09_wp,&
      & 1.09314255139469E-05_wp, 1.13196388765754E-06_wp, 1.56950366755713E-05_wp,&
      & 8.61738339115528E-07_wp, 1.99975128279883E-06_wp, 1.06419841772418E-09_wp,&
      & 7.01988052385574E-10_wp, 1.24197394522869E-07_wp, 0.00000000000000E+00_wp,&
      & 1.44845636917032E-06_wp, 1.61344230592836E-06_wp, 4.24001811501497E-07_wp,&
      & 3.27699589459142E-10_wp, 2.53318718705978E-07_wp, 3.01180075322275E-08_wp,&
      & 2.05512219678339E-10_wp, 1.68341647509018E-09_wp, 2.01694273550196E-06_wp,&
      & 9.97581551461384E-07_wp, 3.32053223246824E-10_wp, 4.34489747120837E-06_wp,&
      & 2.97656362010527E-06_wp, 3.22393893136350E-06_wp, 6.45393876041118E-09_wp,&
      & 1.44845636917032E-06_wp, 0.00000000000000E+00_wp, 9.31459573232734E-06_wp,&
      & 8.11781953663496E-07_wp, 2.46006043874364E-06_wp, 1.42199810287113E-06_wp,&
      & 1.15019379323360E-06_wp, 8.85520127358723E-08_wp, 7.22866025327590E-06_wp,&
      & 1.77979506235292E-06_wp, 4.05095892511959E-05_wp, 4.88101215383324E-06_wp,&
      & 1.33963323128116E-05_wp, 7.46242260907044E-07_wp, 5.40298210888355E-06_wp,&
      & 2.92942083080631E-06_wp, 1.61344230592836E-06_wp, 9.31459573232734E-06_wp,&
      & 0.00000000000000E+00_wp, 6.15676613573339E-06_wp, 7.12122291714440E-07_wp,&
      & 2.56017972323474E-06_wp, 2.94624818524831E-09_wp, 1.94457467193816E-06_wp,&
      & 6.64539527503119E-07_wp, 1.12747911970943E-05_wp, 1.11642484439595E-05_wp,&
      & 1.84973723603671E-06_wp, 1.52794708791284E-06_wp, 1.04617002193055E-06_wp,&
      & 2.86505721882952E-09_wp, 7.02618201266191E-07_wp, 4.24001811501497E-07_wp,&
      & 8.11781953663496E-07_wp, 6.15676613573339E-06_wp, 0.00000000000000E+00_wp,&
      & 5.27610612468442E-07_wp, 9.15344869711176E-07_wp, 5.89013360878893E-06_wp,&
      & 1.62766520767394E-06_wp, 2.47075870779500E-06_wp, 7.15635454034322E-06_wp,&
      & 6.46819148114272E-06_wp, 2.84601678442178E-06_wp, 1.60432360805069E-06_wp,&
      & 4.11497988679676E-09_wp, 3.40154362497530E-10_wp, 8.40186543771942E-07_wp,&
      & 3.27699589459142E-10_wp, 2.46006043874364E-06_wp, 7.12122291714440E-07_wp,&
      & 5.27610612468442E-07_wp, 0.00000000000000E+00_wp, 1.58259611541755E-06_wp,&
      & 2.88648585176153E-09_wp, 1.80076422091063E-08_wp, 1.68458847871656E-06_wp,&
      & 2.17756374040625E-09_wp, 4.62789837511975E-09_wp, 2.73455132586791E-10_wp,&
      & 2.60958239754062E-05_wp, 4.90405805872025E-06_wp, 3.70337380076441E-06_wp,&
      & 6.88715357186242E-09_wp, 2.53318718705978E-07_wp, 1.42199810287113E-06_wp,&
      & 2.56017972323474E-06_wp, 9.15344869711176E-07_wp, 1.58259611541755E-06_wp,&
      & 0.00000000000000E+00_wp, 6.55188367592859E-06_wp, 8.76089280655923E-08_wp,&
      & 2.08671032693074E-05_wp, 2.39782681007998E-06_wp, 2.85564290014006E-05_wp,&
      & 4.34947803753967E-07_wp, 3.28428291428594E-09_wp, 1.02381089813213E-09_wp,&
      & 4.19271847501332E-08_wp, 1.34014614738065E-06_wp, 3.01180075322275E-08_wp,&
      & 1.15019379323360E-06_wp, 2.94624818524831E-09_wp, 5.89013360878893E-06_wp,&
      & 2.88648585176153E-09_wp, 6.55188367592859E-06_wp, 0.00000000000000E+00_wp,&
      & 2.79089675221936E-08_wp, 2.16749904608163E-06_wp, 6.79119587976615E-06_wp,&
      & 7.11246257059381E-09_wp, 2.61186201304170E-07_wp, 1.22038974098681E-06_wp,&
      & 1.17638119754792E-09_wp, 6.18155501580136E-08_wp, 2.35505184164687E-09_wp,&
      & 2.05512219678339E-10_wp, 8.85520127358723E-08_wp, 1.94457467193816E-06_wp,&
      & 1.62766520767394E-06_wp, 1.80076422091063E-08_wp, 8.76089280655923E-08_wp,&
      & 2.79089675221936E-08_wp, 0.00000000000000E+00_wp, 2.26658129830572E-09_wp,&
      & 6.03022595120772E-07_wp, 2.49917765003255E-06_wp, 3.60302393384392E-08_wp,&
      & 1.80416885064417E-06_wp, 1.37693495675094E-06_wp, 4.04805092715763E-06_wp,&
      & 1.09314255139469E-05_wp, 1.68341647509018E-09_wp, 7.22866025327590E-06_wp,&
      & 6.64539527503119E-07_wp, 2.47075870779500E-06_wp, 1.68458847871656E-06_wp,&
      & 2.08671032693074E-05_wp, 2.16749904608163E-06_wp, 2.26658129830572E-09_wp,&
      & 0.00000000000000E+00_wp, 5.25028689064895E-07_wp, 2.06167982840535E-06_wp,&
      & 2.54928574960063E-06_wp, 3.64859359329591E-06_wp, 3.34999221798228E-06_wp,&
      & 5.12707170946965E-06_wp, 1.13196388765754E-06_wp, 2.01694273550196E-06_wp,&
      & 1.77979506235292E-06_wp, 1.12747911970943E-05_wp, 7.15635454034322E-06_wp,&
      & 2.17756374040625E-09_wp, 2.39782681007998E-06_wp, 6.79119587976615E-06_wp,&
      & 6.03022595120772E-07_wp, 5.25028689064895E-07_wp, 0.00000000000000E+00_wp,&
      & 3.19835532404969E-05_wp, 2.11860255652576E-06_wp, 8.97509255934878E-06_wp,&
      & 9.19792252998876E-07_wp, 2.43353063964647E-06_wp, 1.56950366755713E-05_wp,&
      & 9.97581551461384E-07_wp, 4.05095892511959E-05_wp, 1.11642484439595E-05_wp,&
      & 6.46819148114272E-06_wp, 4.62789837511975E-09_wp, 2.85564290014006E-05_wp,&
      & 7.11246257059381E-09_wp, 2.49917765003255E-06_wp, 2.06167982840535E-06_wp,&
      & 3.19835532404969E-05_wp, 0.00000000000000E+00_wp, 1.54811161907370E-06_wp,&
      & 3.44312903091752E-06_wp, 1.61806219757397E-07_wp, 9.58218767201640E-08_wp,&
      & 8.61738339115528E-07_wp, 3.32053223246824E-10_wp, 4.88101215383324E-06_wp,&
      & 1.84973723603671E-06_wp, 2.84601678442178E-06_wp, 2.73455132586791E-10_wp,&
      & 4.34947803753967E-07_wp, 2.61186201304170E-07_wp, 3.60302393384392E-08_wp,&
      & 2.54928574960063E-06_wp, 2.11860255652576E-06_wp, 1.54811161907370E-06_wp,&
      & 0.00000000000000E+00_wp], shape(ref6))

   real(wp), parameter :: ref8(16, 16) = reshape([&
      & 0.00000000000000E+00_wp, 3.38313968908567E-07_wp, 8.56769878673858E-08_wp,&
      & 1.49145901861886E-06_wp, 3.17973306994815E-08_wp, 3.75351688839210E-06_wp,&
      & 6.89216156351992E-06_wp, 1.27077799442787E-08_wp, 1.68292584279487E-07_wp,&
      & 6.23165940790036E-07_wp, 1.68914625862650E-09_wp, 1.66478562972431E-07_wp,&
      & 1.58676277110931E-08_wp, 4.05933177288144E-08_wp, 6.18459294734525E-06_wp,&
      & 4.29149424583112E-08_wp, 3.38313968908567E-07_wp, 0.00000000000000E+00_wp,&
      & 2.37671089648647E-10_wp, 1.47226543172709E-07_wp, 1.81354697395772E-10_wp,&
      & 3.46227613256310E-07_wp, 2.98239498118331E-07_wp, 3.67533186583854E-07_wp,&
      & 7.83920334826650E-10_wp, 1.58394408056639E-07_wp, 2.61557733816910E-10_wp,&
      & 2.18174369986450E-10_wp, 1.98501182885193E-07_wp, 4.09621766709891E-08_wp,&
      & 3.42388519157805E-07_wp, 2.45709629704660E-08_wp, 8.56769878673858E-08_wp,&
      & 2.37671089648647E-10_wp, 0.00000000000000E+00_wp, 2.73452698683199E-08_wp,&
      & 1.32660083859589E-10_wp, 3.61376378314511E-08_wp, 3.90013111261496E-07_wp,&
      & 1.47353390862414E-09_wp, 7.19203172702896E-11_wp, 4.43952670386430E-08_wp,&
      & 1.17666189614926E-08_wp, 1.11058112494132E-08_wp, 9.14492580640533E-08_wp,&
      & 4.17770946971026E-07_wp, 6.34747537200202E-07_wp, 2.23768815580111E-08_wp,&
      & 1.49145901861886E-06_wp, 1.47226543172709E-07_wp, 2.73452698683199E-08_wp,&
      & 0.00000000000000E+00_wp, 2.92960005234579E-08_wp, 4.66422943404331E-09_wp,&
      & 3.03680000229031E-08_wp, 4.51041770333478E-09_wp, 5.50278626897909E-08_wp,&
      & 4.97731116571373E-09_wp, 2.21972557453627E-07_wp, 6.80575189849368E-10_wp,&
      & 4.83152238587731E-07_wp, 8.52073848355318E-09_wp, 3.60355795378375E-07_wp,&
      & 1.50257415068768E-07_wp, 3.17973306994815E-08_wp, 1.81354697395772E-10_wp,&
      & 1.32660083859589E-10_wp, 2.92960005234579E-08_wp, 0.00000000000000E+00_wp,&
      & 8.76703254079051E-08_wp, 1.68357212231353E-08_wp, 1.10117479853480E-07_wp,&
      & 4.62214880046112E-11_wp, 7.81436006761070E-08_wp, 5.38858594576079E-09_wp,&
      & 2.81964610594669E-11_wp, 4.85469799113998E-10_wp, 4.22896400231628E-08_wp,&
      & 1.16082234904963E-07_wp, 5.65927289491721E-11_wp, 3.75351688839210E-06_wp,&
      & 3.46227613256310E-07_wp, 3.61376378314511E-08_wp, 4.66422943404331E-09_wp,&
      & 8.76703254079051E-08_wp, 0.00000000000000E+00_wp, 1.42723249248880E-07_wp,&
      & 5.46766941423818E-09_wp, 1.04595067954914E-07_wp, 1.31673911590776E-06_wp,&
      & 4.90415406421207E-07_wp, 3.28033982379406E-08_wp, 1.02564326893702E-07_wp,&
      & 1.55754993629445E-08_wp, 1.81406347649791E-06_wp, 1.62135552527443E-07_wp,&
      & 6.89216156351992E-06_wp, 2.98239498118331E-07_wp, 3.90013111261496E-07_wp,&
      & 3.03680000229031E-08_wp, 1.68357212231353E-08_wp, 1.42723249248880E-07_wp,&
      & 0.00000000000000E+00_wp, 8.17376337889341E-08_wp, 1.86365349672466E-07_wp,&
      & 2.52976491338952E-08_wp, 1.51529132508726E-09_wp, 2.07518035917784E-08_wp,&
      & 4.18721030592290E-09_wp, 1.85021332250425E-07_wp, 6.69461297437840E-06_wp,&
      & 1.68456042858389E-08_wp, 1.27077799442787E-08_wp, 3.67533186583854E-07_wp,&
      & 1.47353390862414E-09_wp, 4.51041770333478E-09_wp, 1.10117479853480E-07_wp,&
      & 5.46766941423818E-09_wp, 8.17376337889341E-08_wp, 0.00000000000000E+00_wp,&
      & 1.61759036766789E-07_wp, 6.41702071604589E-09_wp, 1.01288663217432E-07_wp,&
      & 1.53048687570244E-08_wp, 2.41477824955660E-08_wp, 1.00031928657220E-07_wp,&
      & 8.73310932937407E-08_wp, 3.15357634809007E-08_wp, 1.68292584279487E-07_wp,&
      & 7.83920334826650E-10_wp, 7.19203172702896E-11_wp, 5.50278626897909E-08_wp,&
      & 4.62214880046112E-11_wp, 1.04595067954914E-07_wp, 1.86365349672466E-07_wp,&
      & 1.61759036766789E-07_wp, 0.00000000000000E+00_wp, 1.69798625308576E-07_wp,&
      & 6.10155982391626E-10_wp, 2.58278731704038E-09_wp, 2.92592995296342E-08_wp,&
      & 8.35905612996506E-10_wp, 1.77649745005639E-09_wp, 5.21430993516171E-11_wp,&
      & 6.23165940790036E-07_wp, 1.58394408056639E-07_wp, 4.43952670386430E-08_wp,&
      & 4.97731116571373E-09_wp, 7.81436006761070E-08_wp, 1.31673911590776E-06_wp,&
      & 2.52976491338952E-08_wp, 6.41702071604589E-09_wp, 1.69798625308576E-07_wp,&
      & 0.00000000000000E+00_wp, 2.61277221738345E-07_wp, 3.24595512527009E-08_wp,&
      & 5.45305656796274E-07_wp, 2.31802156555730E-08_wp, 7.25128991028232E-07_wp,&
      & 1.91471666715141E-07_wp, 1.68914625862650E-09_wp, 2.61557733816910E-10_wp,&
      & 1.17666189614926E-08_wp, 2.21972557453627E-07_wp, 5.38858594576079E-09_wp,&
      & 4.90415406421207E-07_wp, 1.51529132508726E-09_wp, 1.01288663217432E-07_wp,&
      & 6.10155982391626E-10_wp, 2.61277221738345E-07_wp, 0.00000000000000E+00_wp,&
      & 5.57900625399690E-09_wp, 2.16786873630531E-08_wp, 2.14218133798898E-07_wp,&
      & 3.65800626368950E-09_wp, 3.84571244260274E-08_wp, 1.66478562972431E-07_wp,&
      & 2.18174369986450E-10_wp, 1.11058112494132E-08_wp, 6.80575189849368E-10_wp,&
      & 2.81964610594669E-11_wp, 3.28033982379406E-08_wp, 2.07518035917784E-08_wp,&
      & 1.53048687570244E-08_wp, 2.58278731704038E-09_wp, 3.24595512527009E-08_wp,&
      & 5.57900625399690E-09_wp, 0.00000000000000E+00_wp, 7.11339694955262E-10_wp,&
      & 3.71106246499744E-09_wp, 7.34295167143090E-08_wp, 6.12705588066961E-09_wp,&
      & 1.58676277110931E-08_wp, 1.98501182885193E-07_wp, 9.14492580640533E-08_wp,&
      & 4.83152238587731E-07_wp, 4.85469799113998E-10_wp, 1.02564326893702E-07_wp,&
      & 4.18721030592290E-09_wp, 2.41477824955660E-08_wp, 2.92592995296342E-08_wp,&
      & 5.45305656796274E-07_wp, 2.16786873630531E-08_wp, 7.11339694955262E-10_wp,&
      & 0.00000000000000E+00_wp, 3.05818110695778E-09_wp, 1.89613825892962E-08_wp,&
      & 1.13487043936559E-07_wp, 4.05933177288144E-08_wp, 4.09621766709891E-08_wp,&
      & 4.17770946971026E-07_wp, 8.52073848355318E-09_wp, 4.22896400231628E-08_wp,&
      & 1.55754993629445E-08_wp, 1.85021332250425E-07_wp, 1.00031928657220E-07_wp,&
      & 8.35905612996506E-10_wp, 2.31802156555730E-08_wp, 2.14218133798898E-07_wp,&
      & 3.71106246499744E-09_wp, 3.05818110695778E-09_wp, 0.00000000000000E+00_wp,&
      & 4.30198176883569E-06_wp, 4.09322436310400E-07_wp, 6.18459294734525E-06_wp,&
      & 3.42388519157805E-07_wp, 6.34747537200202E-07_wp, 3.60355795378375E-07_wp,&
      & 1.16082234904963E-07_wp, 1.81406347649791E-06_wp, 6.69461297437840E-06_wp,&
      & 8.73310932937407E-08_wp, 1.77649745005639E-09_wp, 7.25128991028232E-07_wp,&
      & 3.65800626368950E-09_wp, 7.34295167143090E-08_wp, 1.89613825892962E-08_wp,&
      & 4.30198176883569E-06_wp, 0.00000000000000E+00_wp, 4.17370943281672E-07_wp,&
      & 4.29149424583112E-08_wp, 2.45709629704660E-08_wp, 2.23768815580111E-08_wp,&
      & 1.50257415068768E-07_wp, 5.65927289491721E-11_wp, 1.62135552527443E-07_wp,&
      & 1.68456042858389E-08_wp, 3.15357634809007E-08_wp, 5.21430993516171E-11_wp,&
      & 1.91471666715141E-07_wp, 3.84571244260274E-08_wp, 6.12705588066961E-09_wp,&
      & 1.13487043936559E-07_wp, 4.09322436310400E-07_wp, 4.17370943281672E-07_wp,&
      & 0.00000000000000E+00_wp], shape(ref8))

   ! PBE-D3(0) (using D4 r4/r2 radii instead of D3 vdW radii)
   param = param_type(&
      & s6 = 1.0_wp, s8 = 0.722_wp, a1 = 1.0_wp, a2 = 0.0_wp, &
      & rs6 = 1.217_wp, rs8 = 1.0_wp, alp = 14.0_wp)

   call get_structure(mol, "MB16-43", "02")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%zero, -1)
   if (allocated(error)) return
   call test_damping_2b_gen(error, mol, d4, damp, param, ref6, ref8)

end subroutine test_damp_zero_2b_mb02

subroutine test_grad_zero_2b_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   ! PBE-D3(0) (using D4 r4/r2 radii instead of D3 vdW radii)
   param = param_type(&
      & s6 = 1.0_wp, s8 = 0.722_wp, a1 = 1.0_wp, a2 = 0.0_wp, &
      & rs6 = 1.217_wp, rs8 = 1.0_wp, alp = 14.0_wp)

   call get_structure(mol, "MB16-43", "03")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%zero, -1)
   if (allocated(error)) return
   call test_damping_2b_numgrad(error, mol, d4, damp, param)

end subroutine test_grad_zero_2b_mb03

subroutine test_grad_zero_2b_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   ! PBE-D3(0) (using D4 r4/r2 radii instead of D3 vdW radii)
   param = param_type(&
      & s6 = 1.0_wp, s8 = 0.722_wp, a1 = 1.0_wp, a2 = 0.0_wp, &
      & rs6 = 1.217_wp, rs8 = 1.0_wp, alp = 14.0_wp)

   call get_structure(mol, "MB16-43", "04")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%zero, -1)
   if (allocated(error)) return
   call test_damping_2b_numgrad(error, mol, d4, damp, param)

end subroutine test_grad_zero_2b_mb04


subroutine test_damp_mzero_2b_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   real(wp), parameter :: ref6(16, 16) = reshape([&
      & 0.00000000000000E+00_wp, 6.71086519398558E-07_wp, 4.22276969437660E-09_wp,&
      & 2.05335563715780E-08_wp, 2.58154905915366E-07_wp, 8.26277614304922E-07_wp,&
      & 6.71830743580388E-08_wp, 4.55322168704039E-07_wp, 8.34110186393240E-09_wp,&
      & 8.80669061621277E-07_wp, 1.03721025158510E-06_wp, 2.43397372948269E-08_wp,&
      & 5.66760122908048E-09_wp, 3.06941111116173E-09_wp, 1.23891696483496E-07_wp,&
      & 2.01497470458597E-09_wp, 6.71086519398558E-07_wp, 0.00000000000000E+00_wp,&
      & 6.43770699293235E-07_wp, 7.98933409379490E-06_wp, 2.11816045239509E-07_wp,&
      & 6.64028423290081E-07_wp, 7.03548821971774E-06_wp, 5.41316312563907E-06_wp,&
      & 2.44546344001824E-06_wp, 8.59520361855053E-06_wp, 4.46157794163729E-07_wp,&
      & 1.09096404914722E-06_wp, 7.87618549789826E-09_wp, 6.67484027681618E-07_wp,&
      & 1.07605447635691E-06_wp, 6.02119810374154E-07_wp, 4.22276969437660E-09_wp,&
      & 6.43770699293235E-07_wp, 0.00000000000000E+00_wp, 3.52888242659943E-06_wp,&
      & 1.89033555350764E-06_wp, 3.84937678555132E-06_wp, 5.82818391671496E-06_wp,&
      & 1.24769150155529E-06_wp, 9.23713174051740E-07_wp, 7.29343218458289E-07_wp,&
      & 5.47423546954281E-06_wp, 9.39968114336764E-07_wp, 1.94941007214675E-06_wp,&
      & 3.95635006166417E-09_wp, 5.51443851342009E-06_wp, 3.36570917151925E-07_wp,&
      & 2.05335563715780E-08_wp, 7.98933409379490E-06_wp, 3.52888242659943E-06_wp,&
      & 0.00000000000000E+00_wp, 6.18729663531687E-07_wp, 1.78765771106392E-06_wp,&
      & 7.63563458734344E-06_wp, 2.48924450266700E-06_wp, 9.55703345434010E-06_wp,&
      & 2.23145936080057E-06_wp, 8.48278902524467E-07_wp, 1.06261906370145E-06_wp,&
      & 8.61195771370002E-09_wp, 2.75423980549028E-06_wp, 2.44832963106352E-06_wp,&
      & 4.46006654454675E-07_wp, 2.58154905915366E-07_wp, 2.11816045239509E-07_wp,&
      & 1.89033555350764E-06_wp, 6.18729663531687E-07_wp, 0.00000000000000E+00_wp,&
      & 1.67862766180292E-05_wp, 1.47753540639153E-06_wp, 1.69638400456181E-06_wp,&
      & 1.37871128407520E-06_wp, 7.08841195475576E-07_wp, 7.47850047659377E-06_wp,&
      & 6.45068321580347E-07_wp, 6.03564689059911E-07_wp, 8.31173261756385E-09_wp,&
      & 3.64334359915934E-06_wp, 5.22860634507504E-07_wp, 8.26277614304922E-07_wp,&
      & 6.64028423290081E-07_wp, 3.84937678555132E-06_wp, 1.78765771106392E-06_wp,&
      & 1.67862766180292E-05_wp, 0.00000000000000E+00_wp, 3.27428429822684E-06_wp,&
      & 5.90213536129803E-06_wp, 1.64951153141424E-06_wp, 1.74972878595017E-06_wp,&
      & 1.84476739439197E-06_wp, 4.99474698415791E-06_wp, 2.17156035291770E-06_wp,&
      & 3.38844399643015E-06_wp, 1.36098559988374E-08_wp, 1.26969478755155E-07_wp,&
      & 6.71830743580388E-08_wp, 7.03548821971774E-06_wp, 5.82818391671496E-06_wp,&
      & 7.63563458734344E-06_wp, 1.47753540639153E-06_wp, 3.27428429822684E-06_wp,&
      & 0.00000000000000E+00_wp, 1.07799049665037E-05_wp, 9.75674939918651E-06_wp,&
      & 2.98988265058451E-05_wp, 3.07204693214464E-06_wp, 9.96994548555582E-07_wp,&
      & 1.00963674505342E-08_wp, 5.19248509003585E-06_wp, 7.85528670509833E-06_wp,&
      & 1.60672498563368E-08_wp, 4.55322168704039E-07_wp, 5.41316312563907E-06_wp,&
      & 1.24769150155529E-06_wp, 2.48924450266700E-06_wp, 1.69638400456181E-06_wp,&
      & 5.90213536129803E-06_wp, 1.07799049665037E-05_wp, 0.00000000000000E+00_wp,&
      & 5.18260108514214E-06_wp, 1.20765602708635E-08_wp, 1.28271938197434E-05_wp,&
      & 1.21622081628146E-06_wp, 1.88796762166374E-06_wp, 2.45484533995458E-06_wp,&
      & 3.90287162348491E-06_wp, 2.47184538238848E-09_wp, 8.34110186393240E-09_wp,&
      & 2.44546344001824E-06_wp, 9.23713174051740E-07_wp, 9.55703345434010E-06_wp,&
      & 1.37871128407520E-06_wp, 1.64951153141424E-06_wp, 9.75674939918651E-06_wp,&
      & 5.18260108514214E-06_wp, 0.00000000000000E+00_wp, 5.58646096084686E-06_wp,&
      & 1.59674849132040E-06_wp, 9.52697896187369E-07_wp, 1.26342101528530E-06_wp,&
      & 6.94019415571521E-09_wp, 3.57520646439345E-08_wp, 2.74452746958803E-09_wp,&
      & 8.80669061621277E-07_wp, 8.59520361855053E-06_wp, 7.29343218458289E-07_wp,&
      & 2.23145936080057E-06_wp, 7.08841195475576E-07_wp, 1.74972878595017E-06_wp,&
      & 2.98988265058451E-05_wp, 1.20765602708635E-08_wp, 5.58646096084686E-06_wp,&
      & 0.00000000000000E+00_wp, 3.73381313587490E-06_wp, 5.12743153387544E-06_wp,&
      & 5.43560408623701E-06_wp, 1.42500202176917E-06_wp, 5.55189038350986E-06_wp,&
      & 5.59587822941605E-08_wp, 1.03721025158510E-06_wp, 4.46157794163729E-07_wp,&
      & 5.47423546954281E-06_wp, 8.48278902524467E-07_wp, 7.47850047659377E-06_wp,&
      & 1.84476739439197E-06_wp, 3.07204693214464E-06_wp, 1.28271938197434E-05_wp,&
      & 1.59674849132040E-06_wp, 3.73381313587490E-06_wp, 0.00000000000000E+00_wp,&
      & 3.00876951215399E-06_wp, 1.27701559681588E-06_wp, 1.45721013990716E-06_wp,&
      & 1.38044125204005E-08_wp, 1.11480364461904E-07_wp, 2.43397372948269E-08_wp,&
      & 1.09096404914722E-06_wp, 9.39968114336764E-07_wp, 1.06261906370145E-06_wp,&
      & 6.45068321580347E-07_wp, 4.99474698415791E-06_wp, 9.96994548555582E-07_wp,&
      & 1.21622081628146E-06_wp, 9.52697896187369E-07_wp, 5.12743153387544E-06_wp,&
      & 3.00876951215399E-06_wp, 0.00000000000000E+00_wp, 8.49096894522788E-09_wp,&
      & 9.76405201581498E-07_wp, 1.02012869903792E-06_wp, 2.36716156559370E-09_wp,&
      & 5.66760122908048E-09_wp, 7.87618549789826E-09_wp, 1.94941007214675E-06_wp,&
      & 8.61195771370002E-09_wp, 6.03564689059911E-07_wp, 2.17156035291770E-06_wp,&
      & 1.00963674505342E-08_wp, 1.88796762166374E-06_wp, 1.26342101528530E-06_wp,&
      & 5.43560408623701E-06_wp, 1.27701559681588E-06_wp, 8.49096894522788E-09_wp,&
      & 0.00000000000000E+00_wp, 1.00325140833800E-06_wp, 2.34140976718037E-06_wp,&
      & 5.14702380525962E-09_wp, 3.06941111116173E-09_wp, 6.67484027681618E-07_wp,&
      & 3.95635006166417E-09_wp, 2.75423980549028E-06_wp, 8.31173261756385E-09_wp,&
      & 3.38844399643015E-06_wp, 5.19248509003585E-06_wp, 2.45484533995458E-06_wp,&
      & 6.94019415571521E-09_wp, 1.42500202176917E-06_wp, 1.45721013990716E-06_wp,&
      & 9.76405201581498E-07_wp, 1.00325140833800E-06_wp, 0.00000000000000E+00_wp,&
      & 2.19215384049659E-07_wp, 1.60095365459624E-08_wp, 1.23891696483496E-07_wp,&
      & 1.07605447635691E-06_wp, 5.51443851342009E-06_wp, 2.44832963106352E-06_wp,&
      & 3.64334359915934E-06_wp, 1.36098559988374E-08_wp, 7.85528670509833E-06_wp,&
      & 3.90287162348491E-06_wp, 3.57520646439345E-08_wp, 5.55189038350986E-06_wp,&
      & 1.38044125204005E-08_wp, 1.02012869903792E-06_wp, 2.34140976718037E-06_wp,&
      & 2.19215384049659E-07_wp, 0.00000000000000E+00_wp, 4.15962323857527E-09_wp,&
      & 2.01497470458597E-09_wp, 6.02119810374154E-07_wp, 3.36570917151925E-07_wp,&
      & 4.46006654454675E-07_wp, 5.22860634507504E-07_wp, 1.26969478755155E-07_wp,&
      & 1.60672498563368E-08_wp, 2.47184538238848E-09_wp, 2.74452746958803E-09_wp,&
      & 5.59587822941605E-08_wp, 1.11480364461904E-07_wp, 2.36716156559370E-09_wp,&
      & 5.14702380525962E-09_wp, 1.60095365459624E-08_wp, 4.15962323857527E-09_wp,&
      & 0.00000000000000E+00_wp], shape(ref6))

   real(wp), parameter :: ref8(16, 16) = reshape([&
      & 0.00000000000000E+00_wp, 1.21898469951458E-07_wp, 5.15179714964540E-09_wp,&
      & 3.50753641756649E-08_wp, 1.07389251983842E-07_wp, 8.84482782580646E-08_wp,&
      & 1.15018407640645E-07_wp, 2.86934049621058E-08_wp, 1.00294101748756E-08_wp,&
      & 1.72973792366646E-08_wp, 4.65228974688172E-08_wp, 1.94725221847150E-08_wp,&
      & 4.89543872022831E-09_wp, 2.56953046152930E-09_wp, 8.14934828100275E-08_wp,&
      & 1.02963846629603E-09_wp, 1.21898469951458E-07_wp, 0.00000000000000E+00_wp,&
      & 8.52365568187441E-09_wp, 1.98331537320159E-05_wp, 1.93562166533445E-09_wp,&
      & 8.88125843310596E-09_wp, 2.05616687095849E-05_wp, 1.50697238780745E-07_wp,&
      & 5.08760886272661E-08_wp, 2.74698111417406E-07_wp, 5.22633229192801E-09_wp,&
      & 2.00928479016612E-06_wp, 2.40214937008432E-08_wp, 8.96632774445039E-09_wp,&
      & 1.69228697501624E-08_wp, 3.57554656038717E-07_wp, 5.15179714964540E-09_wp,&
      & 8.52365568187441E-09_wp, 0.00000000000000E+00_wp, 8.33987999227048E-08_wp,&
      & 4.08985032321254E-06_wp, 9.39067699312449E-08_wp, 1.67365428784939E-07_wp,&
      & 2.07218960836287E-08_wp, 2.38559061186967E-06_wp, 1.00677164359899E-08_wp,&
      & 1.53106163571095E-07_wp, 1.46247387003844E-08_wp, 4.49663189814936E-08_wp,&
      & 8.95566079922316E-09_wp, 2.15760997067542E-07_wp, 1.30266574432536E-07_wp,&
      & 3.50753641756649E-08_wp, 1.98331537320159E-05_wp, 8.33987999227048E-08_wp,&
      & 0.00000000000000E+00_wp, 8.08347923498787E-09_wp, 3.32757774658101E-08_wp,&
      & 2.01338569895564E-05_wp, 5.20094737523046E-08_wp, 3.79121915197992E-07_wp,&
      & 4.47377139395814E-08_wp, 1.23110129426922E-08_wp, 2.00910301580584E-06_wp,&
      & 2.63616767266899E-08_wp, 6.34100668395331E-08_wp, 5.09565991120235E-08_wp,&
      & 4.14558783507467E-07_wp, 1.07389251983842E-07_wp, 1.93562166533445E-09_wp,&
      & 4.08985032321254E-06_wp, 8.08347923498787E-09_wp, 0.00000000000000E+00_wp,&
      & 9.86064822336263E-07_wp, 2.58240009583992E-08_wp, 3.12454125382257E-08_wp,&
      & 3.41637262804322E-06_wp, 9.69060475930760E-09_wp, 7.97362841298503E-06_wp,&
      & 8.61994046302747E-09_wp, 7.87060987724363E-09_wp, 2.12003785105482E-08_wp,&
      & 2.30310503980809E-06_wp, 1.42561733899929E-07_wp, 8.84482782580646E-08_wp,&
      & 8.88125843310596E-09_wp, 9.39067699312449E-08_wp, 3.32757774658101E-08_wp,&
      & 9.86064822336263E-07_wp, 0.00000000000000E+00_wp, 7.46793141639161E-08_wp,&
      & 1.70410080348642E-07_wp, 5.96497358031085E-06_wp, 3.23370143308418E-08_wp,&
      & 1.12591334187055E-05_wp, 4.01970743709781E-07_wp, 4.48276493119789E-08_wp,&
      & 1.08710846068754E-06_wp, 5.76356587971203E-08_wp, 2.51102505682088E-07_wp,&
      & 1.15018407640645E-07_wp, 2.05616687095849E-05_wp, 1.67365428784939E-07_wp,&
      & 2.01338569895564E-05_wp, 2.58240009583992E-08_wp, 7.46793141639161E-08_wp,&
      & 0.00000000000000E+00_wp, 2.92430710365252E-06_wp, 2.36650695726872E-06_wp,&
      & 3.06755765430034E-06_wp, 6.85754387301605E-08_wp, 2.00156147501618E-06_wp,&
      & 3.11034951408092E-08_wp, 2.20019479406905E-07_wp, 2.68528934880020E-07_wp,&
      & 3.19290216170841E-08_wp, 2.86934049621058E-08_wp, 1.50697238780745E-07_wp,&
      & 2.07218960836287E-08_wp, 5.20094737523046E-08_wp, 3.12454125382257E-08_wp,&
      & 1.70410080348642E-07_wp, 2.92430710365252E-06_wp, 0.00000000000000E+00_wp,&
      & 7.63955635227348E-07_wp, 5.33256598714920E-08_wp, 6.53051391502828E-07_wp,&
      & 5.12997090921948E-07_wp, 3.52857705577122E-07_wp, 7.64071400565058E-08_wp,&
      & 1.20899073481871E-06_wp, 3.44329313524202E-09_wp, 1.00294101748756E-08_wp,&
      & 5.08760886272661E-08_wp, 2.38559061186967E-06_wp, 3.79121915197992E-07_wp,&
      & 3.41637262804322E-06_wp, 5.96497358031085E-06_wp, 2.36650695726872E-06_wp,&
      & 7.63955635227348E-07_wp, 0.00000000000000E+00_wp, 1.59620435334065E-07_wp,&
      & 5.87009780750477E-06_wp, 4.87231506500954E-07_wp, 4.28388454773203E-07_wp,&
      & 1.53855434042880E-08_wp, 1.14647008377049E-07_wp, 3.66861888334553E-09_wp,&
      & 1.72973792366646E-08_wp, 2.74698111417406E-07_wp, 1.00677164359899E-08_wp,&
      & 4.47377139395814E-08_wp, 9.69060475930760E-09_wp, 3.23370143308418E-08_wp,&
      & 3.06755765430034E-06_wp, 5.33256598714920E-08_wp, 1.59620435334065E-07_wp,&
      & 0.00000000000000E+00_wp, 8.90311165943919E-08_wp, 2.81561793026435E-07_wp,&
      & 3.84913376127282E-07_wp, 2.49489132391654E-08_wp, 1.58189953385520E-07_wp,&
      & 1.15026278099364E-07_wp, 4.65228974688172E-08_wp, 5.22633229192801E-09_wp,&
      & 1.53106163571095E-07_wp, 1.23110129426922E-08_wp, 7.97362841298503E-06_wp,&
      & 1.12591334187055E-05_wp, 6.85754387301605E-08_wp, 6.53051391502828E-07_wp,&
      & 5.87009780750477E-06_wp, 8.90311165943919E-08_wp, 0.00000000000000E+00_wp,&
      & 7.36788136347899E-08_wp, 2.14852616493390E-08_wp, 2.11403967277541E-06_wp,&
      & 5.84896959781819E-08_wp, 2.23839220557527E-07_wp, 1.94725221847150E-08_wp,&
      & 2.00928479016612E-06_wp, 1.46247387003844E-08_wp, 2.00910301580584E-06_wp,&
      & 8.61994046302747E-09_wp, 4.01970743709781E-07_wp, 2.00156147501618E-06_wp,&
      & 5.12997090921948E-07_wp, 4.87231506500954E-07_wp, 2.81561793026435E-07_wp,&
      & 7.36788136347899E-08_wp, 0.00000000000000E+00_wp, 1.36120369798859E-08_wp,&
      & 3.41227017988999E-08_wp, 4.58520889626603E-07_wp, 2.24040285056455E-09_wp,&
      & 4.89543872022831E-09_wp, 2.40214937008432E-08_wp, 4.49663189814936E-08_wp,&
      & 2.63616767266899E-08_wp, 7.87060987724363E-09_wp, 4.48276493119789E-08_wp,&
      & 3.11034951408092E-08_wp, 3.52857705577122E-07_wp, 4.28388454773203E-07_wp,&
      & 3.84913376127282E-07_wp, 2.14852616493390E-08_wp, 1.36120369798859E-08_wp,&
      & 0.00000000000000E+00_wp, 5.11482079950630E-08_wp, 9.34686001962333E-08_wp,&
      & 5.20573786149825E-09_wp, 2.56953046152930E-09_wp, 8.96632774445039E-09_wp,&
      & 8.95566079922316E-09_wp, 6.34100668395331E-08_wp, 2.12003785105482E-08_wp,&
      & 1.08710846068754E-06_wp, 2.20019479406905E-07_wp, 7.64071400565058E-08_wp,&
      & 1.53855434042880E-08_wp, 2.49489132391654E-08_wp, 2.11403967277541E-06_wp,&
      & 3.41227017988999E-08_wp, 5.11482079950630E-08_wp, 0.00000000000000E+00_wp,&
      & 4.87457409729713E-07_wp, 1.68081456482769E-08_wp, 8.14934828100275E-08_wp,&
      & 1.69228697501624E-08_wp, 2.15760997067542E-07_wp, 5.09565991120235E-08_wp,&
      & 2.30310503980809E-06_wp, 5.76356587971203E-08_wp, 2.68528934880020E-07_wp,&
      & 1.20899073481871E-06_wp, 1.14647008377049E-07_wp, 1.58189953385520E-07_wp,&
      & 5.84896959781819E-08_wp, 4.58520889626603E-07_wp, 9.34686001962333E-08_wp,&
      & 4.87457409729713E-07_wp, 0.00000000000000E+00_wp, 5.68325753769145E-09_wp,&
      & 1.02963846629603E-09_wp, 3.57554656038717E-07_wp, 1.30266574432536E-07_wp,&
      & 4.14558783507467E-07_wp, 1.42561733899929E-07_wp, 2.51102505682088E-07_wp,&
      & 3.19290216170841E-08_wp, 3.44329313524202E-09_wp, 3.66861888334553E-09_wp,&
      & 1.15026278099364E-07_wp, 2.23839220557527E-07_wp, 2.24040285056455E-09_wp,&
      & 5.20573786149825E-09_wp, 1.68081456482769E-08_wp, 5.68325753769145E-09_wp,&
      & 0.00000000000000E+00_wp], shape(ref8))

   ! B3LYP-D3(0M) (using D4 r4/r2 radii instead of D3 vdW radii)
   param = param_type(&
      & s6 = 1.0_wp, s8 = 1.532981_wp, a1 = 1.0_wp, a2 = 0.0_wp, &
      & rs6 = 1.338153_wp, rs8 = 1.0_wp, alp = 14.0_wp, bet=0.013988_wp)

   call get_structure(mol, "MB16-43", "01")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%mzero, -1)
   if (allocated(error)) return
   call test_damping_2b_gen(error, mol, d4, damp, param, ref6, ref8)

end subroutine test_damp_mzero_2b_mb01

subroutine test_damp_mzero_2b_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   real(wp), parameter :: ref6(16, 16) = reshape([&
      & 0.00000000000000E+00_wp, 7.98504500860863E-07_wp, 4.89798365614754E-06_wp,&
      & 9.89201743120253E-06_wp, 1.76216996939039E-06_wp, 2.88711043934147E-06_wp,&
      & 8.01869938049789E-06_wp, 1.52770017648918E-06_wp, 1.35089965751532E-06_wp,&
      & 2.38228564774253E-05_wp, 8.76484715399140E-09_wp, 1.06683313253425E-06_wp,&
      & 1.80288320472036E-06_wp, 3.64446063956308E-06_wp, 5.54703877986797E-06_wp,&
      & 3.21215570673873E-06_wp, 7.98504500860863E-07_wp, 0.00000000000000E+00_wp,&
      & 4.43181839685148E-09_wp, 8.81656218710672E-07_wp, 5.54349147813544E-09_wp,&
      & 2.28418794511425E-06_wp, 6.79248326705810E-07_wp, 9.07165778912931E-07_wp,&
      & 1.35785798029972E-08_wp, 3.85063792981731E-06_wp, 4.72786386320211E-09_wp,&
      & 5.76300448064152E-09_wp, 1.17828307506268E-06_wp, 3.13829921416598E-06_wp,&
      & 8.12211799918111E-07_wp, 2.11109612433207E-07_wp, 4.89798365614754E-06_wp,&
      & 4.43181839685148E-09_wp, 0.00000000000000E+00_wp, 2.21219440840645E-06_wp,&
      & 4.06824573479002E-09_wp, 3.11953076888727E-06_wp, 3.94662969367110E-06_wp,&
      & 7.93122047272150E-09_wp, 2.45289047207329E-09_wp, 3.54212378551807E-06_wp,&
      & 6.99785759899758E-08_wp, 1.00895179628025E-07_wp, 3.34495981511373E-06_wp,&
      & 3.73918982866345E-06_wp, 1.86856676015171E-06_wp, 1.37173253002558E-07_wp,&
      & 9.89201743120253E-06_wp, 8.81656218710672E-07_wp, 2.21219440840645E-06_wp,&
      & 0.00000000000000E+00_wp, 1.67629087056597E-07_wp, 1.33302082632489E-08_wp,&
      & 2.91922161324575E-06_wp, 7.02514549893264E-07_wp, 7.91682960879880E-07_wp,&
      & 1.40067605907471E-08_wp, 1.14751154820593E-06_wp, 8.10195708989375E-09_wp,&
      & 7.93851616879955E-06_wp, 1.13148115817095E-06_wp, 1.35078313514848E-05_wp,&
      & 7.94912225837276E-07_wp, 1.76216996939039E-06_wp, 5.54349147813544E-09_wp,&
      & 4.06824573479002E-09_wp, 1.67629087056597E-07_wp, 0.00000000000000E+00_wp,&
      & 1.25410200244622E-06_wp, 1.50371589535908E-06_wp, 4.40057856038095E-07_wp,&
      & 2.79028948605787E-09_wp, 2.87122986859715E-07_wp, 5.90358572023611E-08_wp,&
      & 2.12572934446326E-09_wp, 6.40497358677365E-09_wp, 1.73953042199488E-06_wp,&
      & 9.05547660565330E-07_wp, 2.61790546927051E-09_wp, 2.88711043934147E-06_wp,&
      & 2.28418794511425E-06_wp, 3.11953076888727E-06_wp, 1.33302082632489E-08_wp,&
      & 1.25410200244622E-06_wp, 0.00000000000000E+00_wp, 9.23031862953183E-06_wp,&
      & 8.11750351938870E-07_wp, 2.01009584217281E-06_wp, 1.07748771107559E-06_wp,&
      & 9.70377552596310E-07_wp, 1.19985080540792E-07_wp, 7.11910867426684E-06_wp,&
      & 1.77939002474727E-06_wp, 2.93903102160927E-05_wp, 3.82236144245935E-06_wp,&
      & 8.01869938049789E-06_wp, 6.79248326705810E-07_wp, 3.94662969367110E-06_wp,&
      & 2.91922161324575E-06_wp, 1.50371589535908E-06_wp, 9.23031862953183E-06_wp,&
      & 0.00000000000000E+00_wp, 6.13445846313322E-06_wp, 6.70220029623945E-07_wp,&
      & 2.55886456220848E-06_wp, 8.09474384883401E-09_wp, 1.81835238062180E-06_wp,&
      & 6.64488561230309E-07_wp, 1.11194474902869E-05_wp, 6.77364894879456E-06_wp,&
      & 1.81661024640718E-06_wp, 1.52770017648918E-06_wp, 9.07165778912931E-07_wp,&
      & 7.93122047272150E-09_wp, 7.02514549893264E-07_wp, 4.40057856038095E-07_wp,&
      & 8.11750351938870E-07_wp, 6.13445846313322E-06_wp, 0.00000000000000E+00_wp,&
      & 5.20231835134893E-07_wp, 9.15298158173242E-07_wp, 5.16551935292340E-06_wp,&
      & 1.55519591656686E-06_wp, 2.46722319575462E-06_wp, 7.12017618864599E-06_wp,&
      & 6.44204834837955E-06_wp, 2.71928456389130E-06_wp, 1.35089965751532E-06_wp,&
      & 1.35785798029972E-08_wp, 2.45289047207329E-09_wp, 7.91682960879880E-07_wp,&
      & 2.79028948605787E-09_wp, 2.01009584217281E-06_wp, 6.70220029623945E-07_wp,&
      & 5.20231835134893E-07_wp, 0.00000000000000E+00_wp, 1.33465901535016E-06_wp,&
      & 1.02532049092283E-08_wp, 4.19362893653177E-08_wp, 1.48739521734509E-06_wp,&
      & 7.05825600154848E-09_wp, 1.21862263870599E-08_wp, 2.22376661334240E-09_wp,&
      & 2.38228564774253E-05_wp, 3.85063792981731E-06_wp, 3.54212378551807E-06_wp,&
      & 1.40067605907471E-08_wp, 2.87122986859715E-07_wp, 1.07748771107559E-06_wp,&
      & 2.55886456220848E-06_wp, 9.15298158173242E-07_wp, 1.33465901535016E-06_wp,&
      & 0.00000000000000E+00_wp, 4.94418505686422E-06_wp, 1.18953190905801E-07_wp,&
      & 1.75313897964021E-05_wp, 2.39676292694007E-06_wp, 2.55125763296910E-05_wp,&
      & 4.29862683683663E-07_wp, 8.76484715399140E-09_wp, 4.72786386320211E-09_wp,&
      & 6.99785759899758E-08_wp, 1.14751154820593E-06_wp, 5.90358572023611E-08_wp,&
      & 9.70377552596310E-07_wp, 8.09474384883401E-09_wp, 5.16551935292340E-06_wp,&
      & 1.02532049092283E-08_wp, 4.94418505686422E-06_wp, 0.00000000000000E+00_wp,&
      & 5.47256609383407E-08_wp, 2.08737107116351E-06_wp, 5.26804858716149E-06_wp,&
      & 1.55799226943174E-08_wp, 3.06309949346805E-07_wp, 1.06683313253425E-06_wp,&
      & 5.76300448064152E-09_wp, 1.00895179628025E-07_wp, 8.10195708989375E-09_wp,&
      & 2.12572934446326E-09_wp, 1.19985080540792E-07_wp, 1.81835238062180E-06_wp,&
      & 1.55519591656686E-06_wp, 4.19362893653177E-08_wp, 1.18953190905801E-07_wp,&
      & 5.47256609383407E-08_wp, 0.00000000000000E+00_wp, 7.69757638271475E-09_wp,&
      & 5.99582673556361E-07_wp, 2.07669729594416E-06_wp, 6.78973404610769E-08_wp,&
      & 1.80288320472036E-06_wp, 1.17828307506268E-06_wp, 3.34495981511373E-06_wp,&
      & 7.93851616879955E-06_wp, 6.40497358677365E-09_wp, 7.11910867426684E-06_wp,&
      & 6.64488561230309E-07_wp, 2.46722319575462E-06_wp, 1.48739521734509E-06_wp,&
      & 1.75313897964021E-05_wp, 2.08737107116351E-06_wp, 7.69757638271475E-09_wp,&
      & 0.00000000000000E+00_wp, 5.25004957917854E-07_wp, 2.05970434761062E-06_wp,&
      & 2.07061941927127E-06_wp, 3.64446063956308E-06_wp, 3.13829921416598E-06_wp,&
      & 3.73918982866345E-06_wp, 1.13148115817095E-06_wp, 1.73953042199488E-06_wp,&
      & 1.77939002474727E-06_wp, 1.11194474902869E-05_wp, 7.12017618864599E-06_wp,&
      & 7.05825600154848E-09_wp, 2.39676292694007E-06_wp, 5.26804858716149E-06_wp,&
      & 5.99582673556361E-07_wp, 5.25004957917854E-07_wp, 0.00000000000000E+00_wp,&
      & 1.91878986278909E-05_wp, 1.67888520650329E-06_wp, 5.54703877986797E-06_wp,&
      & 8.12211799918111E-07_wp, 1.86856676015171E-06_wp, 1.35078313514848E-05_wp,&
      & 9.05547660565330E-07_wp, 2.93903102160927E-05_wp, 6.77364894879456E-06_wp,&
      & 6.44204834837955E-06_wp, 1.21862263870599E-08_wp, 2.55125763296910E-05_wp,&
      & 1.55799226943174E-08_wp, 2.07669729594416E-06_wp, 2.05970434761062E-06_wp,&
      & 1.91878986278909E-05_wp, 0.00000000000000E+00_wp, 1.27376160391766E-06_wp,&
      & 3.21215570673873E-06_wp, 2.11109612433207E-07_wp, 1.37173253002558E-07_wp,&
      & 7.94912225837276E-07_wp, 2.61790546927051E-09_wp, 3.82236144245935E-06_wp,&
      & 1.81661024640718E-06_wp, 2.71928456389130E-06_wp, 2.22376661334240E-09_wp,&
      & 4.29862683683663E-07_wp, 3.06309949346805E-07_wp, 6.78973404610769E-08_wp,&
      & 2.07061941927127E-06_wp, 1.67888520650329E-06_wp, 1.27376160391766E-06_wp,&
      & 0.00000000000000E+00_wp], shape(ref6))

   real(wp), parameter :: ref8(16, 16) = reshape([&
      & 0.00000000000000E+00_wp, 1.48148396464278E-06_wp, 1.82313714057822E-07_wp,&
      & 3.33169730869128E-06_wp, 6.81535152081391E-08_wp, 1.60151831897279E-05_wp,&
      & 1.98068973782486E-05_wp, 2.69817026804452E-08_wp, 4.20179743468692E-07_wp,&
      & 1.32413898600715E-06_wp, 2.68488482750308E-08_wp, 4.49533563848694E-07_wp,&
      & 3.36908726078567E-08_wp, 8.61896633537605E-08_wp, 2.08835878809300E-05_wp,&
      & 9.12348921377013E-08_wp, 1.48148396464278E-06_wp, 0.00000000000000E+00_wp,&
      & 6.33480696392203E-09_wp, 4.20023131213960E-07_wp, 5.22884325308065E-09_wp,&
      & 8.42385301794741E-07_wp, 1.43417822133790E-06_wp, 1.48648269993552E-06_wp,&
      & 1.49262203943442E-08_wp, 3.43196544289534E-07_wp, 6.77971812901394E-09_wp,&
      & 5.95646350009990E-09_wp, 5.36330380679720E-07_wp, 8.70731209292253E-08_wp,&
      & 1.48387782639824E-06_wp, 1.00318405409452E-07_wp, 1.82313714057822E-07_wp,&
      & 6.33480696392203E-09_wp, 0.00000000000000E+00_wp, 5.81697723888975E-08_wp,&
      & 4.21468087712717E-09_wp, 7.67566585403119E-08_wp, 8.79077098427500E-07_wp,&
      & 2.41962808932839E-08_wp, 2.78301051902363E-09_wp, 9.43142504439601E-08_wp,&
      & 1.18742210131302E-07_wp, 7.54262663467248E-08_wp, 1.96440241649525E-07_wp,&
      & 9.52153218328438E-07_wp, 1.88934597019596E-06_wp, 1.57365635663033E-07_wp,&
      & 3.33169730869128E-06_wp, 4.20023131213960E-07_wp, 5.81697723888975E-08_wp,&
      & 0.00000000000000E+00_wp, 1.86375494560774E-07_wp, 5.90737175330945E-08_wp,&
      & 6.44794447878536E-08_wp, 9.57671224030740E-09_wp, 1.33247956108049E-07_wp,&
      & 6.21828540130331E-08_wp, 6.24036718510598E-07_wp, 1.35578896890497E-08_wp,&
      & 1.04664158580279E-06_wp, 1.80916089999918E-08_wp, 7.66873474567603E-07_wp,&
      & 4.56978291584491E-07_wp, 6.81535152081391E-08_wp, 5.22884325308065E-09_wp,&
      & 4.21468087712717E-09_wp, 1.86375494560774E-07_wp, 0.00000000000000E+00_wp,&
      & 2.04339537108286E-07_wp, 3.58349514415628E-08_wp, 4.81280782636885E-07_wp,&
      & 2.07669083747706E-09_wp, 4.59031184303239E-07_wp, 5.14654599274830E-08_wp,&
      & 1.51143399017427E-09_wp, 1.05816872819898E-08_wp, 9.13642228998111E-08_wp,&
      & 3.08174479200771E-07_wp, 2.37310597782308E-09_wp, 1.60151831897279E-05_wp,&
      & 8.42385301794741E-07_wp, 7.67566585403119E-08_wp, 5.90737175330945E-08_wp,&
      & 2.04339537108286E-07_wp, 0.00000000000000E+00_wp, 3.03046474115155E-07_wp,&
      & 1.16091878484505E-08_wp, 2.32614331144122E-07_wp, 6.82545901779947E-06_wp,&
      & 2.15268629735516E-06_wp, 2.79011164760119E-07_wp, 2.17786185194961E-07_wp,&
      & 3.30705711322258E-08_wp, 3.87984035407848E-06_wp, 3.51676515764725E-07_wp,&
      & 1.98068973782486E-05_wp, 1.43417822133790E-06_wp, 8.79077098427500E-07_wp,&
      & 6.44794447878536E-08_wp, 3.58349514415628E-08_wp, 3.03046474115155E-07_wp,&
      & 0.00000000000000E+00_wp, 1.73550657791492E-07_wp, 7.07746648553221E-07_wp,&
      & 5.37130893456718E-08_wp, 2.47158695919099E-08_wp, 4.41422120486552E-08_wp,&
      & 8.89046323114757E-09_wp, 3.92868444530406E-07_wp, 2.07065259512601E-05_wp,&
      & 3.57737439119365E-08_wp, 2.69817026804452E-08_wp, 1.48648269993552E-06_wp,&
      & 2.41962808932839E-08_wp, 9.57671224030740E-09_wp, 4.81280782636885E-07_wp,&
      & 1.16091878484505E-08_wp, 1.73550657791492E-07_wp, 0.00000000000000E+00_wp,&
      & 7.36026258240751E-07_wp, 1.36248910363141E-08_wp, 2.15734067747382E-07_wp,&
      & 3.25275929692621E-08_wp, 5.12717900599789E-08_wp, 2.12395491402786E-07_wp,&
      & 1.85427349715018E-07_wp, 6.70028259582930E-08_wp, 4.20179743468692E-07_wp,&
      & 1.49262203943442E-08_wp, 2.78301051902363E-09_wp, 1.33247956108049E-07_wp,&
      & 2.07669083747706E-09_wp, 2.32614331144122E-07_wp, 7.07746648553221E-07_wp,&
      & 7.36026258240751E-07_wp, 0.00000000000000E+00_wp, 4.26026553875884E-07_wp,&
      & 1.24805552643245E-08_wp, 2.61341463745820E-08_wp, 6.29752086092855E-08_wp,&
      & 1.57944936738505E-08_wp, 2.79001461079696E-08_wp, 2.24893657882281E-09_wp,&
      & 1.32413898600715E-06_wp, 3.43196544289534E-07_wp, 9.43142504439601E-08_wp,&
      & 6.21828540130331E-08_wp, 4.59031184303239E-07_wp, 6.82545901779947E-06_wp,&
      & 5.37130893456718E-08_wp, 1.36248910363141E-08_wp, 4.26026553875884E-07_wp,&
      & 0.00000000000000E+00_wp, 5.68085091288365E-07_wp, 2.76724904443632E-07_wp,&
      & 1.16088620538689E-06_wp, 4.92172533049722E-08_wp, 1.54123574628677E-06_wp,&
      & 1.13607520217337E-06_wp, 2.68488482750308E-08_wp, 6.77971812901394E-09_wp,&
      & 1.18742210131302E-07_wp, 6.24036718510598E-07_wp, 5.14654599274830E-08_wp,&
      & 2.15268629735516E-06_wp, 2.47158695919099E-08_wp, 2.15734067747382E-07_wp,&
      & 1.24805552643245E-08_wp, 5.68085091288365E-07_wp, 0.00000000000000E+00_wp,&
      & 5.91545098568960E-08_wp, 4.60550650897138E-08_wp, 4.61885181551847E-07_wp,&
      & 4.88021292124180E-08_wp, 1.36382057915442E-07_wp, 4.49533563848694E-07_wp,&
      & 5.95646350009990E-09_wp, 7.54262663467248E-08_wp, 1.35578896890497E-08_wp,&
      & 1.51143399017427E-09_wp, 2.79011164760119E-07_wp, 4.41422120486552E-08_wp,&
      & 3.25275929692621E-08_wp, 2.61341463745820E-08_wp, 2.76724904443632E-07_wp,&
      & 5.91545098568960E-08_wp, 0.00000000000000E+00_wp, 1.40099968546971E-08_wp,&
      & 7.87988331717984E-09_wp, 1.60171886071864E-07_wp, 5.12239719541201E-08_wp,&
      & 3.36908726078567E-08_wp, 5.36330380679720E-07_wp, 1.96440241649525E-07_wp,&
      & 1.04664158580279E-06_wp, 1.05816872819898E-08_wp, 2.17786185194961E-07_wp,&
      & 8.89046323114757E-09_wp, 5.12717900599789E-08_wp, 6.29752086092855E-08_wp,&
      & 1.16088620538689E-06_wp, 4.60550650897138E-08_wp, 1.40099968546971E-08_wp,&
      & 0.00000000000000E+00_wp, 6.49326005773950E-09_wp, 4.02597031574862E-08_wp,&
      & 2.52881581666097E-07_wp, 8.61896633537605E-08_wp, 8.70731209292253E-08_wp,&
      & 9.52153218328438E-07_wp, 1.80916089999918E-08_wp, 9.13642228998111E-08_wp,&
      & 3.30705711322258E-08_wp, 3.92868444530406E-07_wp, 2.12395491402786E-07_wp,&
      & 1.57944936738505E-08_wp, 4.92172533049722E-08_wp, 4.61885181551847E-07_wp,&
      & 7.87988331717984E-09_wp, 6.49326005773950E-09_wp, 0.00000000000000E+00_wp,&
      & 9.60590459620985E-06_wp, 1.13532159447882E-06_wp, 2.08835878809300E-05_wp,&
      & 1.48387782639824E-06_wp, 1.88934597019596E-06_wp, 7.66873474567603E-07_wp,&
      & 3.08174479200771E-07_wp, 3.87984035407848E-06_wp, 2.07065259512601E-05_wp,&
      & 1.85427349715018E-07_wp, 2.79001461079696E-08_wp, 1.54123574628677E-06_wp,&
      & 4.88021292124180E-08_wp, 1.60171886071864E-07_wp, 4.02597031574862E-08_wp,&
      & 9.60590459620985E-06_wp, 0.00000000000000E+00_wp, 1.35660927565807E-06_wp,&
      & 9.12348921377013E-08_wp, 1.00318405409452E-07_wp, 1.57365635663033E-07_wp,&
      & 4.56978291584491E-07_wp, 2.37310597782308E-09_wp, 3.51676515764725E-07_wp,&
      & 3.57737439119365E-08_wp, 6.70028259582930E-08_wp, 2.24893657882281E-09_wp,&
      & 1.13607520217337E-06_wp, 1.36382057915442E-07_wp, 5.12239719541201E-08_wp,&
      & 2.52881581666097E-07_wp, 1.13532159447882E-06_wp, 1.35660927565807E-06_wp,&
      & 0.00000000000000E+00_wp], shape(ref8))

   ! B3LYP-D3(0M) (using D4 r4/r2 radii instead of D3 vdW radii)
   param = param_type(&
      & s6 = 1.0_wp, s8 = 1.532981_wp, a1 = 1.0_wp, a2 = 0.0_wp, &
      & rs6 = 1.338153_wp, rs8 = 1.0_wp, alp = 14.0_wp, bet=0.013988_wp)

   call get_structure(mol, "MB16-43", "02")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%mzero, -1)
   if (allocated(error)) return
   call test_damping_2b_gen(error, mol, d4, damp, param, ref6, ref8)

end subroutine test_damp_mzero_2b_mb02

subroutine test_grad_mzero_2b_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   ! B3LYP-D3(0M) (using D4 r4/r2 radii instead of D3 vdW radii)
   param = param_type(&
      & s6 = 1.0_wp, s8 = 1.532981_wp, a1 = 1.0_wp, a2 = 0.0_wp, &
      & rs6 = 1.338153_wp, rs8 = 1.0_wp, alp = 14.0_wp, bet=0.013988_wp)

   call get_structure(mol, "MB16-43", "03")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%mzero, -1)
   if (allocated(error)) return
   call test_damping_2b_numgrad(error, mol, d4, damp, param)

end subroutine test_grad_mzero_2b_mb03

subroutine test_grad_mzero_2b_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   ! B3LYP-D3(0M) (using D4 r4/r2 radii instead of D3 vdW radii)
   param = param_type(&
      & s6 = 1.0_wp, s8 = 1.532981_wp, a1 = 1.0_wp, a2 = 0.0_wp, &
      & rs6 = 1.338153_wp, rs8 = 1.0_wp, alp = 14.0_wp, bet=0.013988_wp)

   call get_structure(mol, "MB16-43", "04")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%mzero, -1)
   if (allocated(error)) return
   call test_damping_2b_numgrad(error, mol, d4, damp, param)

end subroutine test_grad_mzero_2b_mb04


subroutine test_damp_optpower_2b_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   real(wp), parameter :: ref6(16, 16) = reshape([&
      & 0.00000000000000E+00_wp, 4.53751409234784E-06_wp, 2.28093498914889E-06_wp,&
      & 4.46054764473265E-06_wp, 4.60625110886985E-06_wp, 3.59123092721622E-06_wp,&
      & 8.26125258189966E-06_wp, 1.59683147364618E-06_wp, 3.89766734997624E-06_wp,&
      & 1.08676126194321E-06_wp, 2.25179397528647E-06_wp, 4.13444626661782E-06_wp,&
      & 4.12781795452078E-06_wp, 2.81995812537323E-06_wp, 4.66829574082281E-06_wp,&
      & 2.76075671447151E-06_wp, 4.53751409234784E-06_wp, 0.00000000000000E+00_wp,&
      & 6.43202538738371E-07_wp, 1.23494486372834E-05_wp, 2.11785109252475E-07_wp,&
      & 6.63475052932561E-07_wp, 1.15522728413166E-05_wp, 5.34281355479756E-06_wp,&
      & 2.43312100685905E-06_wp, 8.18867068141421E-06_wp, 4.45964153983154E-07_wp,&
      & 1.38841816469217E-05_wp, 7.63353764290331E-07_wp, 6.67679095548494E-07_wp,&
      & 1.07409637427841E-06_wp, 9.13939335312676E-06_wp, 2.28093498914889E-06_wp,&
      & 6.43202538738371E-07_wp, 0.00000000000000E+00_wp, 3.49687924104030E-06_wp,&
      & 1.38013028635696E-05_wp, 3.81121391418054E-06_wp, 5.74989463449390E-06_wp,&
      & 1.24815321334895E-06_wp, 1.22244858667022E-05_wp, 7.28564706260068E-07_wp,&
      & 5.40267613371258E-06_wp, 9.61036003902942E-07_wp, 2.20356060170158E-06_wp,&
      & 7.82998482790917E-07_wp, 6.69614205232480E-06_wp, 5.08180846878463E-06_wp,&
      & 4.46054764473265E-06_wp, 1.23494486372834E-05_wp, 3.49687924104030E-06_wp,&
      & 0.00000000000000E+00_wp, 6.18221719047843E-07_wp, 1.78021844003738E-06_wp,&
      & 1.20602979497379E-05_wp, 2.47420186655528E-06_wp, 9.81702091680599E-06_wp,&
      & 2.21821354346384E-06_wp, 8.47224080418201E-07_wp, 1.37763670223736E-05_wp,&
      & 8.11398927815234E-07_wp, 2.84493836386526E-06_wp, 2.43595781833980E-06_wp,&
      & 1.02342761125968E-05_wp, 4.60625110886985E-06_wp, 2.11785109252475E-07_wp,&
      & 1.38013028635696E-05_wp, 6.18221719047843E-07_wp, 0.00000000000000E+00_wp,&
      & 1.65317572364393E-05_wp, 1.47282621995612E-06_wp, 1.69455967212618E-06_wp,&
      & 1.28898984835748E-05_wp, 7.08119403968648E-07_wp, 1.67265768449121E-05_wp,&
      & 6.47998009788295E-07_wp, 6.05440384891609E-07_wp, 1.12659962616122E-06_wp,&
      & 1.69562037016116E-05_wp, 5.14553341484802E-06_wp, 3.59123092721622E-06_wp,&
      & 6.63475052932561E-07_wp, 3.81121391418054E-06_wp, 1.78021844003738E-06_wp,&
      & 1.65317572364393E-05_wp, 0.00000000000000E+00_wp, 3.23860120158441E-06_wp,&
      & 5.82251350451159E-06_wp, 1.03207398235837E-05_wp, 1.74269406905433E-06_wp,&
      & 5.55856430438795E-06_wp, 9.68841551739198E-06_wp, 2.20817972005972E-06_wp,&
      & 1.47343966347136E-05_wp, 5.82907742989221E-07_wp, 9.35515001292760E-06_wp,&
      & 8.26125258189966E-06_wp, 1.15522728413166E-05_wp, 5.74989463449390E-06_wp,&
      & 1.20602979497379E-05_wp, 1.47282621995612E-06_wp, 3.23860120158441E-06_wp,&
      & 0.00000000000000E+00_wp, 1.97484412088618E-05_wp, 1.91812500142715E-05_wp,&
      & 2.31479899758357E-05_wp, 3.04176410060418E-06_wp, 1.35055601560270E-05_wp,&
      & 9.03915170758513E-07_wp, 6.75815983178248E-06_wp, 7.89628231055961E-06_wp,&
      & 2.92434946074405E-06_wp, 1.59683147364618E-06_wp, 5.34281355479756E-06_wp,&
      & 1.24815321334895E-06_wp, 2.47420186655528E-06_wp, 1.69455967212618E-06_wp,&
      & 5.82251350451159E-06_wp, 1.97484412088618E-05_wp, 0.00000000000000E+00_wp,&
      & 1.33312218130915E-05_wp, 4.87595139505589E-07_wp, 1.34078645827083E-05_wp,&
      & 1.06333098673331E-05_wp, 8.75606917678213E-06_wp, 3.23026649814091E-06_wp,&
      & 1.53524792678599E-05_wp, 1.19605588122354E-06_wp, 3.89766734997624E-06_wp,&
      & 2.43312100685905E-06_wp, 1.22244858667022E-05_wp, 9.81702091680599E-06_wp,&
      & 1.28898984835748E-05_wp, 1.03207398235837E-05_wp, 1.91812500142715E-05_wp,&
      & 1.33312218130915E-05_wp, 0.00000000000000E+00_wp, 5.55174188428394E-06_wp,&
      & 1.01481034474557E-05_wp, 1.04234405114430E-05_wp, 9.72099700526711E-06_wp,&
      & 1.27233980411381E-06_wp, 2.09484923527703E-06_wp, 1.39290508450814E-06_wp,&
      & 1.08676126194321E-06_wp, 8.18867068141421E-06_wp, 7.28564706260068E-07_wp,&
      & 2.21821354346384E-06_wp, 7.08119403968648E-07_wp, 1.74269406905433E-06_wp,&
      & 2.31479899758357E-05_wp, 4.87595139505589E-07_wp, 5.55174188428394E-06_wp,&
      & 0.00000000000000E+00_wp, 3.68382880599655E-06_wp, 7.87566423051389E-06_wp,&
      & 9.49331045143998E-06_wp, 1.43215106912185E-06_wp, 5.51698672275792E-06_wp,&
      & 6.33777898250981E-06_wp, 2.25179397528647E-06_wp, 4.45964153983154E-07_wp,&
      & 5.40267613371258E-06_wp, 8.47224080418201E-07_wp, 1.67265768449121E-05_wp,&
      & 5.55856430438795E-06_wp, 3.04176410060418E-06_wp, 1.34078645827083E-05_wp,&
      & 1.01481034474557E-05_wp, 3.68382880599655E-06_wp, 0.00000000000000E+00_wp,&
      & 3.16950225440771E-06_wp, 1.28160520804594E-06_wp, 1.47805052629972E-05_wp,&
      & 5.88371553776396E-07_wp, 8.89521085261751E-06_wp, 4.13444626661782E-06_wp,&
      & 1.38841816469217E-05_wp, 9.61036003902942E-07_wp, 1.37763670223736E-05_wp,&
      & 6.47998009788295E-07_wp, 9.68841551739198E-06_wp, 1.35055601560270E-05_wp,&
      & 1.06333098673331E-05_wp, 1.04234405114430E-05_wp, 7.87566423051389E-06_wp,&
      & 3.16950225440771E-06_wp, 0.00000000000000E+00_wp, 2.60493976005574E-06_wp,&
      & 1.79502616833638E-06_wp, 1.00997951828103E-05_wp, 1.98141154443650E-06_wp,&
      & 4.12781795452078E-06_wp, 7.63353764290331E-07_wp, 2.20356060170158E-06_wp,&
      & 8.11398927815234E-07_wp, 6.05440384891609E-07_wp, 2.20817972005972E-06_wp,&
      & 9.03915170758513E-07_wp, 8.75606917678213E-06_wp, 9.72099700526711E-06_wp,&
      & 9.49331045143998E-06_wp, 1.28160520804594E-06_wp, 2.60493976005574E-06_wp,&
      & 0.00000000000000E+00_wp, 2.41325789855384E-06_wp, 3.71911609471677E-06_wp,&
      & 3.39030156058252E-06_wp, 2.81995812537323E-06_wp, 6.67679095548494E-07_wp,&
      & 7.82998482790917E-07_wp, 2.84493836386526E-06_wp, 1.12659962616122E-06_wp,&
      & 1.47343966347136E-05_wp, 6.75815983178248E-06_wp, 3.23026649814091E-06_wp,&
      & 1.27233980411381E-06_wp, 1.43215106912185E-06_wp, 1.47805052629972E-05_wp,&
      & 1.79502616833638E-06_wp, 2.41325789855384E-06_wp, 0.00000000000000E+00_wp,&
      & 1.01904663735165E-05_wp, 5.89075317366339E-06_wp, 4.66829574082281E-06_wp,&
      & 1.07409637427841E-06_wp, 6.69614205232480E-06_wp, 2.43595781833980E-06_wp,&
      & 1.69562037016116E-05_wp, 5.82907742989221E-07_wp, 7.89628231055961E-06_wp,&
      & 1.53524792678599E-05_wp, 2.09484923527703E-06_wp, 5.51698672275792E-06_wp,&
      & 5.88371553776396E-07_wp, 1.00997951828103E-05_wp, 3.71911609471677E-06_wp,&
      & 1.01904663735165E-05_wp, 0.00000000000000E+00_wp, 1.90581359721660E-06_wp,&
      & 2.76075671447151E-06_wp, 9.13939335312676E-06_wp, 5.08180846878463E-06_wp,&
      & 1.02342761125968E-05_wp, 5.14553341484802E-06_wp, 9.35515001292760E-06_wp,&
      & 2.92434946074405E-06_wp, 1.19605588122354E-06_wp, 1.39290508450814E-06_wp,&
      & 6.33777898250981E-06_wp, 8.89521085261751E-06_wp, 1.98141154443650E-06_wp,&
      & 3.39030156058252E-06_wp, 5.89075317366339E-06_wp, 1.90581359721660E-06_wp,&
      & 0.00000000000000E+00_wp], shape(ref6))

   real(wp), parameter :: ref8(16, 16) = reshape([&
      & 0.00000000000000E+00_wp, 6.32364743943135E-08_wp, 4.39021921758097E-08_wp,&
      & 9.40589772034651E-08_wp, 6.59736678614997E-08_wp, 4.53205348056845E-08_wp,&
      & 1.80013397541053E-07_wp, 1.49396172241417E-08_wp, 7.49844291080504E-08_wp,&
      & 8.82306787112145E-09_wp, 2.37059681541874E-08_wp, 6.12747562839808E-08_wp,&
      & 7.29958941433630E-08_wp, 4.87966774042347E-08_wp, 6.86735378863763E-08_wp,&
      & 4.12691996335100E-08_wp, 6.32364743943135E-08_wp, 0.00000000000000E+00_wp,&
      & 4.35303673239313E-09_wp, 3.51927867466300E-07_wp, 9.88768954170322E-10_wp,&
      & 4.53594662226795E-09_wp, 3.28390362221155E-07_wp, 7.54451230007061E-08_wp,&
      & 2.58797784038785E-08_wp, 1.35369511464771E-07_wp, 2.66956944379010E-09_wp,&
      & 3.55289290623288E-07_wp, 1.86918759335787E-08_wp, 4.57821031920428E-09_wp,&
      & 8.63785414325859E-09_wp, 1.77428995681375E-07_wp, 4.39021921758097E-08_wp,&
      & 4.35303673239313E-09_wp, 0.00000000000000E+00_wp, 4.22486320510843E-08_wp,&
      & 3.66419993922388E-07_wp, 4.74948545474834E-08_wp, 8.35062804674934E-08_wp,&
      & 1.05700785517978E-08_wp, 3.14193378797574E-07_wp, 5.14118838792299E-09_wp,&
      & 7.66143654558456E-08_wp, 7.46046666054275E-09_wp, 2.28116126506564E-08_wp,&
      & 1.79349134387897E-08_wp, 1.04877830954855E-07_wp, 7.55371839809857E-08_wp,&
      & 9.40589772034651E-08_wp, 3.51927867466300E-07_wp, 4.22486320510843E-08_wp,&
      & 0.00000000000000E+00_wp, 4.12843028897922E-09_wp, 1.69724092052482E-08_wp,&
      & 3.43371190172409E-07_wp, 2.64590213320963E-08_wp, 1.78377853882320E-07_wp,&
      & 2.27989451516078E-08_wp, 6.28678671489163E-09_wp, 3.52412296397130E-07_wp,&
      & 1.98685955072911E-08_wp, 3.21132080524731E-08_wp, 2.59204729699534E-08_wp,&
      & 2.14227967858935E-07_wp, 6.59736678614997E-08_wp, 9.88768954170322E-10_wp,&
      & 3.66419993922388E-07_wp, 4.12843028897922E-09_wp, 0.00000000000000E+00_wp,&
      & 3.85596061464473E-07_wp, 1.31747636109422E-08_wp, 1.59228086816871E-08_wp,&
      & 3.37796522631020E-07_wp, 4.94887476507319E-09_wp, 4.70794113654290E-07_wp,&
      & 4.40085374851028E-09_wp, 4.01864516294056E-09_wp, 2.63884314526765E-08_wp,&
      & 4.45662162153426E-07_wp, 7.60976294643870E-08_wp, 4.53205348056845E-08_wp,&
      & 4.53594662226795E-09_wp, 4.74948545474834E-08_wp, 1.69724092052482E-08_wp,&
      & 3.85596061464473E-07_wp, 0.00000000000000E+00_wp, 3.79521206621985E-08_wp,&
      & 8.49712626448811E-08_wp, 2.76272833873640E-07_wp, 1.64946596070098E-08_wp,&
      & 1.55905428373636E-07_wp, 1.80625461367961E-07_wp, 2.27814765726215E-08_wp,&
      & 3.49652670379961E-07_wp, 1.52995685434314E-08_wp, 2.13232358607083E-07_wp,&
      & 1.80013397541053E-07_wp, 3.28390362221155E-07_wp, 8.35062804674934E-08_wp,&
      & 3.43371190172409E-07_wp, 1.31747636109422E-08_wp, 3.79521206621985E-08_wp,&
      & 0.00000000000000E+00_wp, 5.43459819117352E-07_wp, 5.16261585859265E-07_wp,&
      & 6.62091889101715E-07_wp, 3.48717553862799E-08_wp, 3.45136538142096E-07_wp,&
      & 2.21346554607825E-08_wp, 1.06546169923086E-07_wp, 1.30496398197804E-07_wp,&
      & 6.38877037223032E-08_wp, 1.49396172241417E-08_wp, 7.54451230007061E-08_wp,&
      & 1.05700785517978E-08_wp, 2.64590213320963E-08_wp, 1.59228086816871E-08_wp,&
      & 8.49712626448811E-08_wp, 5.43459819117352E-07_wp, 0.00000000000000E+00_wp,&
      & 2.92966817184002E-07_wp, 1.29238077547227E-08_wp, 2.82012138483036E-07_wp,&
      & 2.18085711228978E-07_wp, 1.60151860485155E-07_wp, 3.84595489287261E-08_wp,&
      & 3.71377150941087E-07_wp, 2.40823389807691E-08_wp, 7.49844291080504E-08_wp,&
      & 2.58797784038785E-08_wp, 3.14193378797574E-07_wp, 1.78377853882320E-07_wp,&
      & 3.37796522631020E-07_wp, 2.76272833873640E-07_wp, 5.16261585859265E-07_wp,&
      & 2.92966817184002E-07_wp, 0.00000000000000E+00_wp, 7.96680975186651E-08_wp,&
      & 2.71469067010499E-07_wp, 2.13450724940440E-07_wp, 1.89825647585690E-07_wp,&
      & 2.87982455244389E-08_wp, 5.12724911006471E-08_wp, 2.76627485075032E-08_wp,&
      & 8.82306787112145E-09_wp, 1.35369511464771E-07_wp, 5.14118838792299E-09_wp,&
      & 2.27989451516078E-08_wp, 4.94887476507319E-09_wp, 1.64946596070098E-08_wp,&
      & 6.62091889101715E-07_wp, 1.29238077547227E-08_wp, 7.96680975186651E-08_wp,&
      & 0.00000000000000E+00_wp, 4.51753091581298E-08_wp, 1.33054757320603E-07_wp,&
      & 1.74727471346673E-07_wp, 1.27173197322637E-08_wp, 7.89783992595800E-08_wp,&
      & 1.41062594167417E-07_wp, 2.37059681541874E-08_wp, 2.66956944379010E-09_wp,&
      & 7.66143654558456E-08_wp, 6.28678671489163E-09_wp, 4.70794113654290E-07_wp,&
      & 1.55905428373636E-07_wp, 3.48717553862799E-08_wp, 2.82012138483036E-07_wp,&
      & 2.71469067010499E-07_wp, 4.51753091581298E-08_wp, 0.00000000000000E+00_wp,&
      & 3.72218233931707E-08_wp, 1.09565171477322E-08_wp, 3.80659980914974E-07_wp,&
      & 1.54429863272843E-08_wp, 2.02064353975638E-07_wp, 6.12747562839808E-08_wp,&
      & 3.55289290623288E-07_wp, 7.46046666054275E-09_wp, 3.52412296397130E-07_wp,&
      & 4.40085374851028E-09_wp, 1.80625461367961E-07_wp, 3.45136538142096E-07_wp,&
      & 2.18085711228978E-07_wp, 2.13450724940440E-07_wp, 1.33054757320603E-07_wp,&
      & 3.72218233931707E-08_wp, 0.00000000000000E+00_wp, 5.39238546642001E-08_wp,&
      & 1.74058636244167E-08_wp, 2.02787721137275E-07_wp, 3.54889840728167E-08_wp,&
      & 7.29958941433630E-08_wp, 1.86918759335787E-08_wp, 2.28116126506564E-08_wp,&
      & 1.98685955072911E-08_wp, 4.01864516294056E-09_wp, 2.27814765726215E-08_wp,&
      & 2.21346554607825E-08_wp, 1.60151860485155E-07_wp, 1.89825647585690E-07_wp,&
      & 1.74727471346673E-07_wp, 1.09565171477322E-08_wp, 5.39238546642001E-08_wp,&
      & 0.00000000000000E+00_wp, 2.60764104108733E-08_wp, 4.68100908797005E-08_wp,&
      & 6.20824516993795E-08_wp, 4.87966774042347E-08_wp, 4.57821031920428E-09_wp,&
      & 1.79349134387897E-08_wp, 3.21132080524731E-08_wp, 2.63884314526765E-08_wp,&
      & 3.49652670379961E-07_wp, 1.06546169923086E-07_wp, 3.84595489287261E-08_wp,&
      & 2.87982455244389E-08_wp, 1.27173197322637E-08_wp, 3.80659980914974E-07_wp,&
      & 1.74058636244167E-08_wp, 2.60764104108733E-08_wp, 0.00000000000000E+00_wp,&
      & 2.41063212681091E-07_wp, 1.11040237184179E-07_wp, 6.86735378863763E-08_wp,&
      & 8.63785414325859E-09_wp, 1.04877830954855E-07_wp, 2.59204729699534E-08_wp,&
      & 4.45662162153426E-07_wp, 1.52995685434314E-08_wp, 1.30496398197804E-07_wp,&
      & 3.71377150941087E-07_wp, 5.12724911006471E-08_wp, 7.89783992595800E-08_wp,&
      & 1.54429863272843E-08_wp, 2.02787721137275E-07_wp, 4.68100908797005E-08_wp,&
      & 2.41063212681091E-07_wp, 0.00000000000000E+00_wp, 3.79094747460024E-08_wp,&
      & 4.12691996335100E-08_wp, 1.77428995681375E-07_wp, 7.55371839809857E-08_wp,&
      & 2.14227967858935E-07_wp, 7.60976294643870E-08_wp, 2.13232358607083E-07_wp,&
      & 6.38877037223032E-08_wp, 2.40823389807691E-08_wp, 2.76627485075032E-08_wp,&
      & 1.41062594167417E-07_wp, 2.02064353975638E-07_wp, 3.54889840728167E-08_wp,&
      & 6.20824516993795E-08_wp, 1.11040237184179E-07_wp, 3.79094747460024E-08_wp,&
      & 0.00000000000000E+00_wp], shape(ref8))

   ! B3LYP-D3(op)
   param = param_type(&
      & s6 = 1.0_wp, s8 = 0.78311_wp, a1 = 0.30_wp, a2 = 4.25_wp, bet=4.0_wp)

   call get_structure(mol, "MB16-43", "01")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%optpower, -1)
   if (allocated(error)) return
   call test_damping_2b_gen(error, mol, d4, damp, param, ref6, ref8)

end subroutine test_damp_optpower_2b_mb01

subroutine test_damp_optpower_2b_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   real(wp), parameter :: ref6(16, 16) = reshape([&
      & 0.00000000000000E+00_wp, 1.33667526140109E-05_wp, 5.97277144531062E-06_wp,&
      & 1.95793403772385E-05_wp, 2.96643415763569E-06_wp, 7.13376863752912E-06_wp,&
      & 1.23731325958912E-05_wp, 1.52276395677972E-06_wp, 9.61873253209367E-06_wp,&
      & 1.94884570294028E-05_wp, 8.21188483317195E-07_wp, 9.98963849163361E-06_wp,&
      & 1.79512391176352E-06_wp, 3.59748374962813E-06_wp, 1.01810552941586E-05_wp,&
      & 3.68337997693672E-06_wp, 1.33667526140109E-05_wp, 0.00000000000000E+00_wp,&
      & 1.87203460117897E-06_wp, 9.70872471421717E-06_wp, 3.82629438073943E-06_wp,&
      & 1.32377540011231E-05_wp, 1.27531438696494E-05_wp, 1.37738296723562E-05_wp,&
      & 5.56662922870609E-06_wp, 8.76297972968094E-06_wp, 1.96282729493118E-06_wp,&
      & 3.58785572423950E-06_wp, 1.08588497087850E-05_wp, 3.56326999546552E-06_wp,&
      & 1.34254367026028E-05_wp, 4.57667130535908E-06_wp, 5.97277144531062E-06_wp,&
      & 1.87203460117897E-06_wp, 0.00000000000000E+00_wp, 2.65640259500625E-06_wp,&
      & 2.75867188380896E-06_wp, 3.26641849815590E-06_wp, 1.38551836023078E-05_wp,&
      & 7.67003354788988E-07_wp, 1.60620511561977E-06_wp, 3.78552118682441E-06_wp,&
      & 8.40362072018497E-06_wp, 4.88006283620467E-06_wp, 6.18771106061767E-06_wp,&
      & 1.42094212323316E-05_wp, 1.54183389647785E-05_wp, 8.07326463170093E-06_wp,&
      & 1.95793403772385E-05_wp, 9.70872471421717E-06_wp, 2.65640259500625E-06_wp,&
      & 0.00000000000000E+00_wp, 8.29977060216204E-06_wp, 5.20349915582215E-07_wp,&
      & 2.89781505511052E-06_wp, 7.01806178852054E-07_wp, 4.79276910993979E-06_wp,&
      & 5.37529361786253E-07_wp, 1.16030784034296E-05_wp, 2.34297476713772E-06_wp,&
      & 1.55457914690998E-05_wp, 1.12915829158447E-06_wp, 1.44946239389790E-05_wp,&
      & 1.01792480344766E-05_wp, 2.96643415763569E-06_wp, 3.82629438073943E-06_wp,&
      & 2.75867188380896E-06_wp, 8.29977060216204E-06_wp, 0.00000000000000E+00_wp,&
      & 6.29808733010180E-06_wp, 1.86231196089148E-06_wp, 1.09670581401908E-05_wp,&
      & 2.92123454078739E-06_wp, 1.12100133586650E-05_wp, 4.99981772942420E-06_wp,&
      & 2.46882247846634E-06_wp, 1.99268695033221E-06_wp, 3.65062394726259E-06_wp,&
      & 8.17534798311016E-06_wp, 2.27348009716206E-06_wp, 7.13376863752912E-06_wp,&
      & 1.32377540011231E-05_wp, 3.26641849815590E-06_wp, 5.20349915582215E-07_wp,&
      & 6.29808733010180E-06_wp, 0.00000000000000E+00_wp, 8.74633687973410E-06_wp,&
      & 8.10810921443359E-07_wp, 6.83260722891744E-06_wp, 4.10112649348709E-06_wp,&
      & 1.30874699035609E-05_wp, 7.95127627092669E-06_wp, 6.92472751600628E-06_wp,&
      & 1.77204011513183E-06_wp, 2.31266533967621E-05_wp, 8.88506801640934E-06_wp,&
      & 1.23731325958912E-05_wp, 1.27531438696494E-05_wp, 1.38551836023078E-05_wp,&
      & 2.89781505511052E-06_wp, 1.86231196089148E-06_wp, 8.74633687973410E-06_wp,&
      & 0.00000000000000E+00_wp, 5.95840603711244E-06_wp, 1.22509413043608E-05_wp,&
      & 2.53997564266140E-06_wp, 7.77792482422420E-07_wp, 2.17152210731337E-06_wp,&
      & 6.63877688907866E-07_wp, 1.03599964077264E-05_wp, 1.13232236348158E-05_wp,&
      & 1.86783720321871E-06_wp, 1.52276395677972E-06_wp, 1.37738296723562E-05_wp,&
      & 7.67003354788988E-07_wp, 7.01806178852054E-07_wp, 1.09670581401908E-05_wp,&
      & 8.10810921443359E-07_wp, 5.95840603711244E-06_wp, 0.00000000000000E+00_wp,&
      & 1.23721289624316E-05_wp, 9.14009756921576E-07_wp, 6.67255892388518E-06_wp,&
      & 1.73653931321344E-06_wp, 2.45035830999960E-06_wp, 6.86479944476029E-06_wp,&
      & 6.24305693088261E-06_wp, 2.95327057668660E-06_wp, 9.61873253209367E-06_wp,&
      & 5.56662922870609E-06_wp, 1.60620511561977E-06_wp, 4.79276910993979E-06_wp,&
      & 2.92123454078739E-06_wp, 6.83260722891744E-06_wp, 1.22509413043608E-05_wp,&
      & 1.23721289624316E-05_wp, 0.00000000000000E+00_wp, 9.68900039192982E-06_wp,&
      & 4.41164660210862E-06_wp, 2.89555056350540E-06_wp, 2.80130751305417E-06_wp,&
      & 1.26732521137762E-06_wp, 1.84574387206583E-06_wp, 1.72456754906960E-06_wp,&
      & 1.94884570294028E-05_wp, 8.76297972968094E-06_wp, 3.78552118682441E-06_wp,&
      & 5.37529361786253E-07_wp, 1.12100133586650E-05_wp, 4.10112649348709E-06_wp,&
      & 2.53997564266140E-06_wp, 9.14009756921576E-07_wp, 9.68900039192982E-06_wp,&
      & 0.00000000000000E+00_wp, 1.16131522016546E-05_wp, 7.91644273728343E-06_wp,&
      & 1.76012112140726E-05_wp, 2.38081923839281E-06_wp, 2.05014580972116E-05_wp,&
      & 1.06907801597797E-05_wp, 8.21188483317195E-07_wp, 1.96282729493118E-06_wp,&
      & 8.40362072018497E-06_wp, 1.16030784034296E-05_wp, 4.99981772942420E-06_wp,&
      & 1.30874699035609E-05_wp, 7.77792482422420E-07_wp, 6.67255892388518E-06_wp,&
      & 4.41164660210862E-06_wp, 1.16131522016546E-05_wp, 0.00000000000000E+00_wp,&
      & 6.51807094260904E-06_wp, 2.24579846541067E-06_wp, 1.04787576225949E-05_wp,&
      & 1.20822895095446E-06_wp, 5.34609792331310E-06_wp, 9.98963849163361E-06_wp,&
      & 3.58785572423950E-06_wp, 4.88006283620467E-06_wp, 2.34297476713772E-06_wp,&
      & 2.46882247846634E-06_wp, 7.95127627092669E-06_wp, 2.17152210731337E-06_wp,&
      & 1.73653931321344E-06_wp, 2.89555056350540E-06_wp, 7.91644273728343E-06_wp,&
      & 6.51807094260904E-06_wp, 0.00000000000000E+00_wp, 1.95540925316546E-06_wp,&
      & 6.05748834433472E-07_wp, 5.36907197330584E-06_wp, 4.30371451837279E-06_wp,&
      & 1.79512391176352E-06_wp, 1.08588497087850E-05_wp, 6.18771106061767E-06_wp,&
      & 1.55457914690998E-05_wp, 1.99268695033221E-06_wp, 6.92472751600628E-06_wp,&
      & 6.63877688907866E-07_wp, 2.45035830999960E-06_wp, 2.80130751305417E-06_wp,&
      & 1.76012112140726E-05_wp, 2.24579846541067E-06_wp, 1.95540925316546E-06_wp,&
      & 0.00000000000000E+00_wp, 5.24673148716667E-07_wp, 2.04889670218354E-06_wp,&
      & 7.19852145732336E-06_wp, 3.59748374962813E-06_wp, 3.56326999546552E-06_wp,&
      & 1.42094212323316E-05_wp, 1.12915829158447E-06_wp, 3.65062394726259E-06_wp,&
      & 1.77204011513183E-06_wp, 1.03599964077264E-05_wp, 6.86479944476029E-06_wp,&
      & 1.26732521137762E-06_wp, 2.38081923839281E-06_wp, 1.04787576225949E-05_wp,&
      & 6.05748834433472E-07_wp, 5.24673148716667E-07_wp, 0.00000000000000E+00_wp,&
      & 1.91421856762549E-05_wp, 1.42408812471664E-05_wp, 1.01810552941586E-05_wp,&
      & 1.34254367026028E-05_wp, 1.54183389647785E-05_wp, 1.44946239389790E-05_wp,&
      & 8.17534798311016E-06_wp, 2.31266533967621E-05_wp, 1.13232236348158E-05_wp,&
      & 6.24305693088261E-06_wp, 1.84574387206583E-06_wp, 2.05014580972116E-05_wp,&
      & 1.20822895095446E-06_wp, 5.36907197330584E-06_wp, 2.04889670218354E-06_wp,&
      & 1.91421856762549E-05_wp, 0.00000000000000E+00_wp, 1.43913617204881E-05_wp,&
      & 3.68337997693672E-06_wp, 4.57667130535908E-06_wp, 8.07326463170093E-06_wp,&
      & 1.01792480344766E-05_wp, 2.27348009716206E-06_wp, 8.88506801640934E-06_wp,&
      & 1.86783720321871E-06_wp, 2.95327057668660E-06_wp, 1.72456754906960E-06_wp,&
      & 1.06907801597797E-05_wp, 5.34609792331310E-06_wp, 4.30371451837279E-06_wp,&
      & 7.19852145732336E-06_wp, 1.42408812471664E-05_wp, 1.43913617204881E-05_wp,&
      & 0.00000000000000E+00_wp], shape(ref6))

   real(wp), parameter :: ref8(16, 16) = reshape([&
      & 0.00000000000000E+00_wp, 3.35364136986283E-07_wp, 8.94153926747017E-08_wp,&
      & 5.42462497778662E-07_wp, 3.43890625843630E-08_wp, 2.00596870163115E-07_wp,&
      & 3.52629715880317E-07_wp, 1.37678300958219E-08_wp, 1.86378147160415E-07_wp,&
      & 4.87627320929193E-07_wp, 2.01083673063426E-08_wp, 1.99075749394554E-07_wp,&
      & 1.71772740596519E-08_wp, 4.37474756286207E-08_wp, 2.88276773926833E-07_wp,&
      & 4.58340584426447E-08_wp, 3.35364136986283E-07_wp, 0.00000000000000E+00_wp,&
      & 3.77043691737067E-08_wp, 1.91893100311171E-07_wp, 6.90215830097287E-08_wp,&
      & 2.99636981334192E-07_wp, 3.19386768085265E-07_wp, 3.45472834875674E-07_wp,&
      & 1.06305293827155E-07_wp, 1.57001131894376E-07_wp, 3.95453180462193E-08_wp,&
      & 6.61775629646020E-08_wp, 2.25557916815608E-07_wp, 4.37913775317971E-08_wp,&
      & 3.36857208876853E-07_wp, 6.57000734700838E-08_wp, 8.94153926747017E-08_wp,&
      & 3.77043691737067E-08_wp, 0.00000000000000E+00_wp, 2.94187981506990E-08_wp,&
      & 5.07072197096132E-08_wp, 3.87628595842265E-08_wp, 3.14750318644833E-07_wp,&
      & 1.87812596974155E-08_wp, 3.04002342368421E-08_wp, 4.74361201386913E-08_wp,&
      & 1.83076014057780E-07_wp, 7.38082644078717E-08_wp, 9.51465264528266E-08_wp,&
      & 3.28381579503881E-07_wp, 3.94793280108555E-07_wp, 1.56265658955887E-07_wp,&
      & 5.42462497778662E-07_wp, 1.91893100311171E-07_wp, 2.94187981506990E-08_wp,&
      & 0.00000000000000E+00_wp, 1.61535026218437E-07_wp, 1.37920169878076E-08_wp,&
      & 3.27516219732684E-08_wp, 4.89057655149497E-09_wp, 6.81374953495221E-08_wp,&
      & 1.42473845401921E-08_wp, 2.50434734689608E-07_wp, 4.90815941136755E-08_wp,&
      & 3.66358250599633E-07_wp, 9.23408435954459E-09_wp, 3.17800261558626E-07_wp,&
      & 2.07222501432201E-07_wp, 3.43890625843630E-08_wp, 6.90215830097287E-08_wp,&
      & 5.07072197096132E-08_wp, 1.61535026218437E-07_wp, 0.00000000000000E+00_wp,&
      & 9.98472602559524E-08_wp, 1.82091829839405E-08_wp, 2.39464129011740E-07_wp,&
      & 4.89200266018652E-08_wp, 2.56534054124561E-07_wp, 7.81468457217415E-08_wp,&
      & 4.06320391680674E-08_wp, 4.16723095983350E-08_wp, 4.58488459695493E-08_wp,&
      & 1.47807876355687E-07_wp, 4.01860913703209E-08_wp, 2.00596870163115E-07_wp,&
      & 2.99636981334192E-07_wp, 3.87628595842265E-08_wp, 1.37920169878076E-08_wp,&
      & 9.98472602559524E-08_wp, 0.00000000000000E+00_wp, 1.48507051721338E-07_wp,&
      & 5.92856390254893E-09_wp, 1.11096080368343E-07_wp, 1.14831039105291E-07_wp,&
      & 3.35362931847121E-07_wp, 1.84264443332695E-07_wp, 1.07815647630029E-07_wp,&
      & 1.68679834757750E-08_wp, 6.72565821944595E-07_wp, 1.60274834970302E-07_wp,&
      & 3.52629715880317E-07_wp, 3.19386768085265E-07_wp, 3.14750318644833E-07_wp,&
      & 3.27516219732684E-08_wp, 1.82091829839405E-08_wp, 1.48507051721338E-07_wp,&
      & 0.00000000000000E+00_wp, 8.70553309221111E-08_wp, 2.80608794391833E-07_wp,&
      & 2.73521794759115E-08_wp, 1.90455027867258E-08_wp, 2.23916592584356E-08_wp,&
      & 4.54041741051787E-09_wp, 1.88866194407788E-07_wp, 3.21657677421330E-07_wp,&
      & 1.81974796108419E-08_wp, 1.37678300958219E-08_wp, 3.45472834875674E-07_wp,&
      & 1.87812596974155E-08_wp, 4.89057655149497E-09_wp, 2.39464129011740E-07_wp,&
      & 5.92856390254893E-09_wp, 8.70553309221111E-08_wp, 0.00000000000000E+00_wp,&
      & 2.90945255201667E-07_wp, 6.95733905363337E-09_wp, 1.04626848335819E-07_wp,&
      & 1.65406200761680E-08_wp, 2.60964805097942E-08_wp, 1.05864003160136E-07_wp,&
      & 9.28378494274742E-08_wp, 3.38652134037087E-08_wp, 1.86378147160415E-07_wp,&
      & 1.06305293827155E-07_wp, 3.04002342368421E-08_wp, 6.81374953495221E-08_wp,&
      & 4.89200266018652E-08_wp, 1.11096080368343E-07_wp, 2.80608794391833E-07_wp,&
      & 2.90945255201667E-07_wp, 0.00000000000000E+00_wp, 1.88524595958996E-07_wp,&
      & 8.54555959234894E-08_wp, 3.55444957988440E-08_wp, 3.18556494307875E-08_wp,&
      & 2.87501575908618E-08_wp, 4.19003320011497E-08_wp, 3.15402370452959E-08_wp,&
      & 4.87627320929193E-07_wp, 1.57001131894376E-07_wp, 4.74361201386913E-08_wp,&
      & 1.42473845401921E-08_wp, 2.56534054124561E-07_wp, 1.14831039105291E-07_wp,&
      & 2.73521794759115E-08_wp, 6.95733905363337E-09_wp, 1.88524595958996E-07_wp,&
      & 0.00000000000000E+00_wp, 2.37322023515433E-07_wp, 1.83411277135384E-07_wp,&
      & 4.26249385009395E-07_wp, 2.50725079357949E-08_wp, 5.29761448021107E-07_wp,&
      & 2.64614944120299E-07_wp, 2.01083673063426E-08_wp, 3.95453180462193E-08_wp,&
      & 1.83076014057780E-07_wp, 2.50434734689608E-07_wp, 7.81468457217415E-08_wp,&
      & 3.35362931847121E-07_wp, 1.90455027867258E-08_wp, 1.04626848335819E-07_wp,&
      & 8.54555959234894E-08_wp, 2.37322023515433E-07_wp, 0.00000000000000E+00_wp,&
      & 1.17556313967583E-07_wp, 2.33707071054951E-08_wp, 2.02693506371111E-07_wp,&
      & 2.95902335606418E-08_wp, 8.15137931060279E-08_wp, 1.99075749394554E-07_wp,&
      & 6.61775629646020E-08_wp, 7.38082644078717E-08_wp, 4.90815941136755E-08_wp,&
      & 4.06320391680674E-08_wp, 1.84264443332695E-07_wp, 2.23916592584356E-08_wp,&
      & 1.65406200761680E-08_wp, 3.55444957988440E-08_wp, 1.83411277135384E-07_wp,&
      & 1.17556313967583E-07_wp, 0.00000000000000E+00_wp, 4.19302249374329E-08_wp,&
      & 4.02297634700321E-09_wp, 7.85715226303798E-08_wp, 6.22079473325066E-08_wp,&
      & 1.71772740596519E-08_wp, 2.25557916815608E-07_wp, 9.51465264528266E-08_wp,&
      & 3.66358250599633E-07_wp, 4.16723095983350E-08_wp, 1.07815647630029E-07_wp,&
      & 4.54041741051787E-09_wp, 2.60964805097942E-08_wp, 3.18556494307875E-08_wp,&
      & 4.26249385009395E-07_wp, 2.33707071054951E-08_wp, 4.19302249374329E-08_wp,&
      & 0.00000000000000E+00_wp, 3.31647895335784E-09_wp, 2.05141984417781E-08_wp,&
      & 1.19751135659404E-07_wp, 4.37474756286207E-08_wp, 4.37913775317971E-08_wp,&
      & 3.28381579503881E-07_wp, 9.23408435954459E-09_wp, 4.58488459695493E-08_wp,&
      & 1.68679834757750E-08_wp, 1.88866194407788E-07_wp, 1.05864003160136E-07_wp,&
      & 2.87501575908618E-08_wp, 2.50725079357949E-08_wp, 2.02693506371111E-07_wp,&
      & 4.02297634700321E-09_wp, 3.31647895335784E-09_wp, 0.00000000000000E+00_wp,&
      & 5.59330632431777E-07_wp, 3.42586885498081E-07_wp, 2.88276773926833E-07_wp,&
      & 3.36857208876853E-07_wp, 3.94793280108555E-07_wp, 3.17800261558626E-07_wp,&
      & 1.47807876355687E-07_wp, 6.72565821944595E-07_wp, 3.21657677421330E-07_wp,&
      & 9.28378494274742E-08_wp, 4.19003320011497E-08_wp, 5.29761448021107E-07_wp,&
      & 2.95902335606418E-08_wp, 7.85715226303798E-08_wp, 2.05141984417781E-08_wp,&
      & 5.59330632431777E-07_wp, 0.00000000000000E+00_wp, 3.56436862923370E-07_wp,&
      & 4.58340584426447E-08_wp, 6.57000734700838E-08_wp, 1.56265658955887E-07_wp,&
      & 2.07222501432201E-07_wp, 4.01860913703209E-08_wp, 1.60274834970302E-07_wp,&
      & 1.81974796108419E-08_wp, 3.38652134037087E-08_wp, 3.15402370452959E-08_wp,&
      & 2.64614944120299E-07_wp, 8.15137931060279E-08_wp, 6.22079473325066E-08_wp,&
      & 1.19751135659404E-07_wp, 3.42586885498081E-07_wp, 3.56436862923370E-07_wp,&
      & 0.00000000000000E+00_wp], shape(ref8))

   ! B3LYP-D3(op)
   param = param_type(&
      & s6 = 1.0_wp, s8 = 0.78311_wp, a1 = 0.30_wp, a2 = 4.25_wp, bet=4.0_wp)

   call get_structure(mol, "MB16-43", "02")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%optpower, -1)
   if (allocated(error)) return
   call test_damping_2b_gen(error, mol, d4, damp, param, ref6, ref8)

end subroutine test_damp_optpower_2b_mb02

subroutine test_grad_optpower_2b_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   ! B3LYP-D3(op)
   param = param_type(&
      & s6 = 1.0_wp, s8 = 0.78311_wp, a1 = 0.30_wp, a2 = 4.25_wp, bet=4.0_wp)

   call get_structure(mol, "MB16-43", "03")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%optpower, -1)
   if (allocated(error)) return
   call test_damping_2b_numgrad(error, mol, d4, damp, param)

end subroutine test_grad_optpower_2b_mb03

subroutine test_grad_optpower_2b_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   ! B3LYP-D3(op)
   param = param_type(&
      & s6 = 1.0_wp, s8 = 0.78311_wp, a1 = 0.30_wp, a2 = 4.25_wp, bet=4.0_wp)

   call get_structure(mol, "MB16-43", "04")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%optpower, -1)
   if (allocated(error)) return
   call test_damping_2b_numgrad(error, mol, d4, damp, param)

end subroutine test_grad_optpower_2b_mb04


subroutine test_damp_cso_2b_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   real(wp), parameter :: ref6(16, 16) = reshape([&
      & 0.00000000000000E+00_wp, 7.21954610218185E-06_wp, 2.91044100485606E-05_wp,&
      & 2.81229494547426E-05_wp, 7.56084971394266E-06_wp, 5.82267984152246E-06_wp,&
      & 2.30904232658061E-05_wp, 2.77639185258663E-06_wp, 2.62145073703324E-05_wp,&
      & 1.91177935002093E-06_wp, 3.79162624788233E-06_wp, 8.08159278273943E-06_wp,&
      & 2.04263192062399E-05_wp, 2.47951107755410E-05_wp, 7.94288755041460E-06_wp,&
      & 1.39943877261144E-05_wp, 7.21954610218185E-06_wp, 0.00000000000000E+00_wp,&
      & 7.76325314166202E-07_wp, 2.94540577905322E-05_wp, 2.14599905287923E-07_wp,&
      & 7.03447541235027E-07_wp, 2.96244981817746E-05_wp, 7.45128957143892E-06_wp,&
      & 3.65114363242784E-06_wp, 9.91553282060012E-06_wp, 4.57113708900886E-07_wp,&
      & 2.50041606059446E-05_wp, 3.11370275434631E-05_wp, 1.04956748917678E-06_wp,&
      & 1.49603183734453E-06_wp, 1.37574419685040E-05_wp, 2.91044100485606E-05_wp,&
      & 7.76325314166202E-07_wp, 0.00000000000000E+00_wp, 5.10336404013263E-06_wp,&
      & 2.71027343114045E-05_wp, 5.52703611528623E-06_wp, 7.92863658615115E-06_wp,&
      & 2.00120870585132E-06_wp, 2.70192092932689E-05_wp, 9.04078679981771E-07_wp,&
      & 7.52229805602183E-06_wp, 1.67279741693846E-06_wp, 3.66281325418424E-06_wp,&
      & 3.10898138145489E-05_wp, 9.38473026954073E-06_wp, 8.19584681715848E-06_wp,&
      & 2.81229494547426E-05_wp, 2.94540577905322E-05_wp, 5.10336404013263E-06_wp,&
      & 0.00000000000000E+00_wp, 7.02284814789231E-07_wp, 2.17895904670661E-06_wp,&
      & 2.95168683785348E-05_wp, 3.64488254429234E-06_wp, 1.22403940157566E-05_wp,&
      & 2.79465276440766E-06_wp, 9.25540186928027E-07_wp, 2.51204772044108E-05_wp,&
      & 3.11305323498021E-05_wp, 4.52766180390250E-06_wp, 3.65539157324152E-06_wp,&
      & 1.59843208476159E-05_wp, 7.56084971394266E-06_wp, 2.14599905287923E-07_wp,&
      & 2.71027343114045E-05_wp, 7.02284814789231E-07_wp, 0.00000000000000E+00_wp,&
      & 1.77293893040822E-05_wp, 1.97590968969825E-06_wp, 2.66164091110503E-06_wp,&
      & 2.72713665872278E-05_wp, 8.25295126356291E-07_wp, 2.74613220531758E-05_wp,&
      & 1.10803334336718E-06_wp, 1.02363874632347E-06_wp, 3.10364946215907E-05_wp,&
      & 2.31120660349012E-05_wp, 8.14965030704975E-06_wp, 5.82267984152246E-06_wp,&
      & 7.03447541235027E-07_wp, 5.52703611528623E-06_wp, 2.17895904670661E-06_wp,&
      & 1.77293893040822E-05_wp, 0.00000000000000E+00_wp, 4.20451197138566E-06_wp,&
      & 8.01246143826856E-06_wp, 2.91828966984039E-05_wp, 2.12641119216651E-06_wp,&
      & 3.06628952121091E-05_wp, 1.27894581863013E-05_wp, 3.58747574885183E-06_wp,&
      & 1.89302296010045E-05_wp, 3.11747959862119E-05_wp, 2.40765019428028E-05_wp,&
      & 2.30904232658061E-05_wp, 2.96244981817746E-05_wp, 7.92863658615115E-06_wp,&
      & 2.95168683785348E-05_wp, 1.97590968969825E-06_wp, 4.20451197138566E-06_wp,&
      & 0.00000000000000E+00_wp, 2.36045983242636E-05_wp, 2.26584938431090E-05_wp,&
      & 2.34428477259037E-05_wp, 3.93761424766232E-06_wp, 2.53948686629791E-05_wp,&
      & 3.11174911860238E-05_wp, 9.49844439188331E-06_wp, 1.03395174797729E-05_wp,&
      & 3.00083626294379E-05_wp, 2.77639185258663E-06_wp, 7.45128957143892E-06_wp,&
      & 2.00120870585132E-06_wp, 3.64488254429234E-06_wp, 2.66164091110503E-06_wp,&
      & 8.01246143826856E-06_wp, 2.36045983242636E-05_wp, 0.00000000000000E+00_wp,&
      & 1.66037689832874E-05_wp, 3.11820333182898E-05_wp, 1.53813173307863E-05_wp,&
      & 1.50617783116314E-05_wp, 1.23324844145714E-05_wp, 5.17099304704220E-06_wp,&
      & 1.94745262654924E-05_wp, 3.06358545098983E-05_wp, 2.62145073703324E-05_wp,&
      & 3.65114363242784E-06_wp, 2.70192092932689E-05_wp, 1.22403940157566E-05_wp,&
      & 2.72713665872278E-05_wp, 2.91828966984039E-05_wp, 2.26584938431090E-05_wp,&
      & 1.66037689832874E-05_wp, 0.00000000000000E+00_wp, 7.76747971199152E-06_wp,&
      & 2.92353547626740E-05_wp, 1.50693723950491E-05_wp, 1.38235463197693E-05_wp,&
      & 3.09440040842352E-05_wp, 3.08956941099971E-05_wp, 3.04132365183011E-05_wp,&
      & 1.91177935002093E-06_wp, 9.91553282060012E-06_wp, 9.04078679981771E-07_wp,&
      & 2.79465276440766E-06_wp, 8.25295126356291E-07_wp, 2.12641119216651E-06_wp,&
      & 2.34428477259037E-05_wp, 3.11820333182898E-05_wp, 7.76747971199152E-06_wp,&
      & 0.00000000000000E+00_wp, 4.79696948110858E-06_wp, 1.07916837058346E-05_wp,&
      & 1.25228020929737E-05_wp, 2.36445998630318E-06_wp, 7.72662295573573E-06_wp,&
      & 2.74610615401257E-05_wp, 3.79162624788233E-06_wp, 4.57113708900886E-07_wp,&
      & 7.52229805602183E-06_wp, 9.25540186928027E-07_wp, 2.74613220531758E-05_wp,&
      & 3.06628952121091E-05_wp, 3.93761424766232E-06_wp, 1.53813173307863E-05_wp,&
      & 2.92353547626740E-05_wp, 4.79696948110858E-06_wp, 0.00000000000000E+00_wp,&
      & 4.99997618858190E-06_wp, 2.11536232967112E-06_wp, 2.42762586017987E-05_wp,&
      & 3.11744007737288E-05_wp, 2.47243903216222E-05_wp, 8.08159278273943E-06_wp,&
      & 2.50041606059446E-05_wp, 1.67279741693846E-06_wp, 2.51204772044108E-05_wp,&
      & 1.10803334336718E-06_wp, 1.27894581863013E-05_wp, 2.53948686629791E-05_wp,&
      & 1.50617783116314E-05_wp, 1.50693723950491E-05_wp, 1.07916837058346E-05_wp,&
      & 4.99997618858190E-06_wp, 0.00000000000000E+00_wp, 2.97056709614169E-05_wp,&
      & 3.07398415370347E-06_wp, 1.45384687352259E-05_wp, 2.83312666113635E-05_wp,&
      & 2.04263192062399E-05_wp, 3.11370275434631E-05_wp, 3.66281325418424E-06_wp,&
      & 3.11305323498021E-05_wp, 1.02363874632347E-06_wp, 3.58747574885183E-06_wp,&
      & 3.11174911860238E-05_wp, 1.23324844145714E-05_wp, 1.38235463197693E-05_wp,&
      & 1.25228020929737E-05_wp, 2.11536232967112E-06_wp, 2.97056709614169E-05_wp,&
      & 0.00000000000000E+00_wp, 4.04316821553614E-06_wp, 5.87317256913425E-06_wp,&
      & 2.54786827682883E-05_wp, 2.47951107755410E-05_wp, 1.04956748917678E-06_wp,&
      & 3.10898138145489E-05_wp, 4.52766180390250E-06_wp, 3.10364946215907E-05_wp,&
      & 1.89302296010045E-05_wp, 9.49844439188331E-06_wp, 5.17099304704220E-06_wp,&
      & 3.09440040842352E-05_wp, 2.36445998630318E-06_wp, 2.42762586017987E-05_wp,&
      & 3.07398415370347E-06_wp, 4.04316821553614E-06_wp, 0.00000000000000E+00_wp,&
      & 2.51031961902857E-05_wp, 1.76220425152720E-05_wp, 7.94288755041460E-06_wp,&
      & 1.49603183734453E-06_wp, 9.38473026954073E-06_wp, 3.65539157324152E-06_wp,&
      & 2.31120660349012E-05_wp, 3.11747959862119E-05_wp, 1.03395174797729E-05_wp,&
      & 1.94745262654924E-05_wp, 3.08956941099971E-05_wp, 7.72662295573573E-06_wp,&
      & 3.11744007737288E-05_wp, 1.45384687352259E-05_wp, 5.87317256913425E-06_wp,&
      & 2.51031961902857E-05_wp, 0.00000000000000E+00_wp, 2.99518223527900E-05_wp,&
      & 1.39943877261144E-05_wp, 1.37574419685040E-05_wp, 8.19584681715848E-06_wp,&
      & 1.59843208476159E-05_wp, 8.14965030704975E-06_wp, 2.40765019428028E-05_wp,&
      & 3.00083626294379E-05_wp, 3.06358545098983E-05_wp, 3.04132365183011E-05_wp,&
      & 2.74610615401257E-05_wp, 2.47243903216222E-05_wp, 2.83312666113635E-05_wp,&
      & 2.54786827682883E-05_wp, 1.76220425152720E-05_wp, 2.99518223527900E-05_wp,&
      & 0.00000000000000E+00_wp], shape(ref6))

   ! B3LYP-D4(cso)
   param = param_type(&
      & s6 = 1.0_wp, a1 = 1.0_wp, a2 = 6.25_wp, a3 = 0.86_wp, a4 = 2.5_wp)

   call get_structure(mol, "MB16-43", "01")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%cso, -1)
   if (allocated(error)) return
   call test_damping_2b_gen(error, mol, d4, damp, param, ref6)

end subroutine test_damp_cso_2b_mb01

subroutine test_damp_cso_2b_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   real(wp), parameter :: ref6(16, 16) = reshape([&
      & 0.00000000000000E+00_wp, 2.43000010913675E-05_wp, 8.58252575241756E-06_wp,&
      & 2.42070045019778E-05_wp, 4.83183715528291E-06_wp, 3.04295712852862E-05_wp,&
      & 2.94488618815561E-05_wp, 1.82028926156310E-06_wp, 1.36344568045398E-05_wp,&
      & 1.91626376395400E-05_wp, 3.11291852550130E-05_wp, 1.43408554830367E-05_wp,&
      & 2.47493878975118E-06_wp, 4.68333670530977E-06_wp, 2.98986579568937E-05_wp,&
      & 5.74568091202493E-06_wp, 2.43000010913675E-05_wp, 0.00000000000000E+00_wp,&
      & 3.00928327659778E-05_wp, 1.41579679185762E-05_wp, 2.32444343607691E-05_wp,&
      & 1.76501755664030E-05_wp, 2.50493332657260E-05_wp, 2.36619744131097E-05_wp,&
      & 2.04775961125126E-05_wp, 1.19393438927170E-05_wp, 3.00131316156748E-05_wp,&
      & 2.52302904135600E-05_wp, 1.54068163419404E-05_wp, 5.58156441268021E-06_wp,&
      & 2.42173102278172E-05_wp, 7.58012899578817E-06_wp, 8.58252575241756E-06_wp,&
      & 3.00928327659778E-05_wp, 0.00000000000000E+00_wp, 4.34093916544585E-06_wp,&
      & 2.72360799115303E-05_wp, 5.12383700695206E-06_wp, 1.75584591275585E-05_wp,&
      & 3.11365410066731E-05_wp, 2.98061958931816E-05_wp, 5.83307678687040E-06_wp,&
      & 2.27079186648955E-05_wp, 8.41676920082710E-06_wp, 9.00510039449592E-06_wp,&
      & 1.80699417174436E-05_wp, 2.29727075129656E-05_wp, 1.41729685382475E-05_wp,&
      & 2.42070045019778E-05_wp, 1.41579679185762E-05_wp, 4.34093916544585E-06_wp,&
      & 0.00000000000000E+00_wp, 1.42273909842792E-05_wp, 3.11799471668091E-05_wp,&
      & 4.26459826801520E-06_wp, 8.63623067448581E-07_wp, 7.54558713010900E-06_wp,&
      & 3.11788290712811E-05_wp, 1.64984841227649E-05_wp, 3.00384759698485E-05_wp,&
      & 1.83410134697032E-05_wp, 1.53619977718647E-06_wp, 1.63365449222837E-05_wp,&
      & 1.49509138424825E-05_wp, 4.83183715528291E-06_wp, 2.32444343607691E-05_wp,&
      & 2.72360799115303E-05_wp, 1.42273909842792E-05_wp, 0.00000000000000E+00_wp,&
      & 9.45116958466236E-06_wp, 3.15763551508444E-06_wp, 1.74582463787922E-05_wp,&
      & 2.25750309958254E-05_wp, 2.03921625911410E-05_wp, 9.04159319998638E-06_wp,&
      & 2.37928703060361E-05_wp, 3.02808430428247E-05_wp, 5.81894492581201E-06_wp,&
      & 1.20604464038713E-05_wp, 2.73018563882278E-05_wp, 3.04295712852862E-05_wp,&
      & 1.76501755664030E-05_wp, 5.12383700695206E-06_wp, 3.11799471668091E-05_wp,&
      & 9.45116958466236E-06_wp, 0.00000000000000E+00_wp, 1.04568279385174E-05_wp,&
      & 8.80584409120106E-07_wp, 9.96582885898119E-06_wp, 3.08522000876677E-05_wp,&
      & 2.60416291776344E-05_wp, 2.70456225002849E-05_wp, 9.08006566280255E-06_wp,&
      & 2.16750062864409E-06_wp, 2.44424787458908E-05_wp, 1.20795569230188E-05_wp,&
      & 2.94488618815561E-05_wp, 2.50493332657260E-05_wp, 1.75584591275585E-05_wp,&
      & 4.26459826801520E-06_wp, 3.15763551508444E-06_wp, 1.04568279385174E-05_wp,&
      & 0.00000000000000E+00_wp, 7.56701953908077E-06_wp, 1.89316290662374E-05_wp,&
      & 3.24530008106241E-06_wp, 3.11350960449229E-05_wp, 3.62538145140410E-06_wp,&
      & 7.64190910855114E-07_wp, 1.19394050952368E-05_wp, 2.96719222172290E-05_wp,&
      & 3.10701786501833E-06_wp, 1.82028926156310E-06_wp, 2.36619744131097E-05_wp,&
      & 3.11365410066731E-05_wp, 8.63623067448581E-07_wp, 1.74582463787922E-05_wp,&
      & 8.80584409120106E-07_wp, 7.56701953908077E-06_wp, 0.00000000000000E+00_wp,&
      & 2.06218247618925E-05_wp, 1.00911370185354E-06_wp, 9.40002324135158E-06_wp,&
      & 2.94863761021160E-06_wp, 3.46636846011453E-06_wp, 8.56006056685999E-06_wp,&
      & 7.88507232216160E-06_wp, 4.72740476266026E-06_wp, 1.36344568045398E-05_wp,&
      & 2.04775961125126E-05_wp, 2.98061958931816E-05_wp, 7.54558713010900E-06_wp,&
      & 2.25750309958254E-05_wp, 9.96582885898119E-06_wp, 1.89316290662374E-05_wp,&
      & 2.06218247618925E-05_wp, 0.00000000000000E+00_wp, 1.37343034439179E-05_wp,&
      & 2.52996638354403E-05_wp, 5.22787887061988E-06_wp, 4.60324407500863E-06_wp,&
      & 3.09499154926972E-05_wp, 3.07584466879696E-05_wp, 2.92033594539849E-05_wp,&
      & 1.91626376395400E-05_wp, 1.19393438927170E-05_wp, 5.83307678687040E-06_wp,&
      & 3.11788290712811E-05_wp, 2.03921625911410E-05_wp, 3.08522000876677E-05_wp,&
      & 3.24530008106241E-06_wp, 1.00911370185354E-06_wp, 1.37343034439179E-05_wp,&
      & 0.00000000000000E+00_wp, 1.48338309792304E-05_wp, 2.70745298171097E-05_wp,&
      & 1.86735826480592E-05_wp, 3.02284926087310E-06_wp, 2.00086749939017E-05_wp,&
      & 2.68292690471312E-05_wp, 3.11291852550130E-05_wp, 3.00131316156748E-05_wp,&
      & 2.27079186648955E-05_wp, 1.64984841227649E-05_wp, 9.04159319998638E-06_wp,&
      & 2.60416291776344E-05_wp, 3.11350960449229E-05_wp, 9.40002324135158E-06_wp,&
      & 2.52996638354403E-05_wp, 1.48338309792304E-05_wp, 0.00000000000000E+00_wp,&
      & 1.25684414958557E-05_wp, 3.70843935092366E-06_wp, 1.35883469517094E-05_wp,&
      & 3.10700347068988E-05_wp, 8.62826298956265E-06_wp, 1.43408554830367E-05_wp,&
      & 2.52302904135600E-05_wp, 8.41676920082710E-06_wp, 3.00384759698485E-05_wp,&
      & 2.37928703060361E-05_wp, 2.70456225002849E-05_wp, 3.62538145140410E-06_wp,&
      & 2.94863761021160E-06_wp, 5.22787887061988E-06_wp, 2.70745298171097E-05_wp,&
      & 1.25684414958557E-05_wp, 0.00000000000000E+00_wp, 3.04597596753017E-05_wp,&
      & 1.06293315990166E-06_wp, 8.10871381970167E-06_wp, 7.61327318626073E-06_wp,&
      & 2.47493878975118E-06_wp, 1.54068163419404E-05_wp, 9.00510039449592E-06_wp,&
      & 1.83410134697032E-05_wp, 3.02808430428247E-05_wp, 9.08006566280255E-06_wp,&
      & 7.64190910855114E-07_wp, 3.46636846011453E-06_wp, 4.60324407500863E-06_wp,&
      & 1.86735826480592E-05_wp, 3.70843935092366E-06_wp, 3.04597596753017E-05_wp,&
      & 0.00000000000000E+00_wp, 5.79494973182554E-07_wp, 2.86365325782075E-06_wp,&
      & 1.04039576804406E-05_wp, 4.68333670530977E-06_wp, 5.58156441268021E-06_wp,&
      & 1.80699417174436E-05_wp, 1.53619977718647E-06_wp, 5.81894492581201E-06_wp,&
      & 2.16750062864409E-06_wp, 1.19394050952368E-05_wp, 8.56006056685999E-06_wp,&
      & 3.09499154926972E-05_wp, 3.02284926087310E-06_wp, 1.35883469517094E-05_wp,&
      & 1.06293315990166E-06_wp, 5.79494973182554E-07_wp, 0.00000000000000E+00_wp,&
      & 2.74963046831576E-05_wp, 1.99584543456620E-05_wp, 2.98986579568937E-05_wp,&
      & 2.42173102278172E-05_wp, 2.29727075129656E-05_wp, 1.63365449222837E-05_wp,&
      & 1.20604464038713E-05_wp, 2.44424787458908E-05_wp, 2.96719222172290E-05_wp,&
      & 7.88507232216160E-06_wp, 3.07584466879696E-05_wp, 2.00086749939017E-05_wp,&
      & 3.10700347068988E-05_wp, 8.10871381970167E-06_wp, 2.86365325782075E-06_wp,&
      & 2.74963046831576E-05_wp, 0.00000000000000E+00_wp, 2.17536835420974E-05_wp,&
      & 5.74568091202493E-06_wp, 7.58012899578817E-06_wp, 1.41729685382475E-05_wp,&
      & 1.49509138424825E-05_wp, 2.73018563882278E-05_wp, 1.20795569230188E-05_wp,&
      & 3.10701786501833E-06_wp, 4.72740476266026E-06_wp, 2.92033594539849E-05_wp,&
      & 2.68292690471312E-05_wp, 8.62826298956265E-06_wp, 7.61327318626073E-06_wp,&
      & 1.04039576804406E-05_wp, 1.99584543456620E-05_wp, 2.17536835420974E-05_wp,&
      & 0.00000000000000E+00_wp], shape(ref6))

   ! B3LYP-D4(cso)
   param = param_type(&
      & s6 = 1.0_wp, a1 = 1.0_wp, a2 = 6.25_wp, a3 = 0.86_wp, a4 = 2.5_wp)

   call get_structure(mol, "MB16-43", "02")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%cso, -1)
   if (allocated(error)) return
   call test_damping_2b_gen(error, mol, d4, damp, param, ref6)

end subroutine test_damp_cso_2b_mb02

subroutine test_grad_cso_2b_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   ! B3LYP-D4(cso)
   param = param_type(&
      & s6 = 1.0_wp, a1 = 1.0_wp, a2 = 6.25_wp, a3 = 0.86_wp, a4 = 2.5_wp)

   call get_structure(mol, "MB16-43", "03")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%cso, -1)
   if (allocated(error)) return
   call test_damping_2b_numgrad(error, mol, d4, damp, param)

end subroutine test_grad_cso_2b_mb03

subroutine test_grad_cso_2b_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   ! B3LYP-D4(cso)
   param = param_type(&
      & s6 = 1.0_wp, a1 = 1.0_wp, a2 = 6.25_wp, a3 = 0.86_wp, a4 = 2.5_wp)

   call get_structure(mol, "MB16-43", "04")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%cso, -1)
   if (allocated(error)) return
   call test_damping_2b_numgrad(error, mol, d4, damp, param)

end subroutine test_grad_cso_2b_mb04


subroutine test_damp_koide_2b_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   real(wp), parameter :: ref6(16, 16) = reshape([&
      & 0.00000000000000E+00_wp, 5.90548829293452E-09_wp, 6.25412185525528E-09_wp,&
      & 7.94154869196974E-09_wp, 5.09266977787729E-09_wp, 5.68059539613240E-09_wp,&
      & 7.41079232430687E-09_wp, 3.95051398712273E-09_wp, 5.70801747300151E-09_wp,&
      & 4.55465299974758E-09_wp, 5.24168330531565E-09_wp, 3.28263620693583E-09_wp,&
      & 3.88992627439869E-09_wp, 4.05611303397122E-09_wp, 4.54964824678742E-09_wp,&
      & 2.16169368803934E-09_wp, 5.90548829293452E-09_wp, 0.00000000000000E+00_wp,&
      & 7.38498518552886E-09_wp, 2.07894649699503E-08_wp, 5.52906061024305E-09_wp,&
      & 8.46300631706588E-09_wp, 2.09117534391124E-08_wp, 1.18221120249956E-08_wp,&
      & 9.90513149886915E-09_wp, 1.47367092527110E-08_wp, 7.54799324087481E-09_wp,&
      & 1.23182242085933E-08_wp, 1.50742955605302E-08_wp, 6.15180823167874E-09_wp,&
      & 8.23381234411572E-09_wp, 7.58294130818892E-09_wp, 6.25412185525528E-09_wp,&
      & 7.38498518552886E-09_wp, 0.00000000000000E+00_wp, 1.09298552909373E-08_wp,&
      & 1.46861833522848E-08_wp, 1.11112878637269E-08_wp, 1.19767266857234E-08_wp,&
      & 7.54339261292298E-09_wp, 1.33496475811201E-08_wp, 7.63700293055344E-09_wp,&
      & 1.18455620254499E-08_wp, 5.60086400495995E-09_wp, 6.77404356579160E-09_wp,&
      & 1.20305417727104E-08_wp, 1.01571587419510E-08_wp, 5.50861108840061E-09_wp,&
      & 7.94154869196974E-09_wp, 2.07894649699503E-08_wp, 1.09298552909373E-08_wp,&
      & 0.00000000000000E+00_wp, 7.62269330938668E-09_wp, 1.08888420474014E-08_wp,&
      & 2.08335500650099E-08_wp, 1.01986668259122E-08_wp, 1.27738958871549E-08_wp,&
      & 1.14439123342815E-08_wp, 9.04786080669313E-09_wp, 1.23398237082115E-08_wp,&
      & 1.50512409160924E-08_wp, 8.49938711303868E-09_wp, 9.90751948646844E-09_wp,&
      & 7.81358375552647E-09_wp, 5.09266977787729E-09_wp, 5.52906061024305E-09_wp,&
      & 1.46861833522848E-08_wp, 7.62269330938668E-09_wp, 0.00000000000000E+00_wp,&
      & 1.51374830105515E-08_wp, 9.51936321638474E-09_wp, 8.49485540647830E-09_wp,&
      & 1.42628562323390E-08_wp, 7.91126235900418E-09_wp, 1.76524169688548E-08_wp,&
      & 5.38978831692585E-09_wp, 5.37826574973416E-09_wp, 1.28271174433308E-08_wp,&
      & 1.33434107842202E-08_wp, 5.92505478150488E-09_wp, 5.68059539613240E-09_wp,&
      & 8.46300631706588E-09_wp, 1.11112878637269E-08_wp, 1.08888420474014E-08_wp,&
      & 1.51374830105515E-08_wp, 0.00000000000000E+00_wp, 1.24003284734753E-08_wp,&
      & 1.20031745551227E-08_wp, 1.68765559123843E-08_wp, 1.08351978268207E-08_wp,&
      & 2.19833234904247E-08_wp, 1.03649522763642E-08_wp, 8.08827272083890E-09_wp,&
      & 1.15301681771158E-08_wp, 1.91600879363687E-08_wp, 8.62398908126134E-09_wp,&
      & 7.41079232430687E-09_wp, 2.09117534391124E-08_wp, 1.19767266857234E-08_wp,&
      & 2.08335500650099E-08_wp, 9.51936321638474E-09_wp, 1.24003284734753E-08_wp,&
      & 0.00000000000000E+00_wp, 1.56328342277459E-08_wp, 1.49605830017022E-08_wp,&
      & 1.83685157373003E-08_wp, 1.22419597732037E-08_wp, 1.23915622937806E-08_wp,&
      & 1.50088137018247E-08_wp, 9.90459432279881E-09_wp, 1.23125224768497E-08_wp,&
      & 9.57571547166595E-09_wp, 3.95051398712273E-09_wp, 1.18221120249956E-08_wp,&
      & 7.54339261292298E-09_wp, 1.01986668259122E-08_wp, 8.49485540647830E-09_wp,&
      & 1.20031745551227E-08_wp, 1.56328342277459E-08_wp, 0.00000000000000E+00_wp,&
      & 1.14921757058533E-08_wp, 1.98682763234573E-08_wp, 1.38562524185172E-08_wp,&
      & 8.81668545779098E-09_wp, 8.62153079348212E-09_wp, 7.26961374828290E-09_wp,&
      & 1.19560236944380E-08_wp, 7.65762342571572E-09_wp, 5.70801747300151E-09_wp,&
      & 9.90513149886915E-09_wp, 1.33496475811201E-08_wp, 1.27738958871549E-08_wp,&
      & 1.42628562323390E-08_wp, 1.68765559123843E-08_wp, 1.49605830017022E-08_wp,&
      & 1.14921757058533E-08_wp, 0.00000000000000E+00_wp, 1.15885544519407E-08_wp,&
      & 1.69021382451719E-08_wp, 8.50339640981889E-09_wp, 8.51283172059904E-09_wp,&
      & 1.13688357004018E-08_wp, 1.45313076476434E-08_wp, 7.24825910759652E-09_wp,&
      & 4.55465299974758E-09_wp, 1.47367092527110E-08_wp, 7.63700293055344E-09_wp,&
      & 1.14439123342815E-08_wp, 7.91126235900418E-09_wp, 1.08351978268207E-08_wp,&
      & 1.83685157373003E-08_wp, 1.98682763234573E-08_wp, 1.15885544519407E-08_wp,&
      & 0.00000000000000E+00_wp, 1.27253555639217E-08_wp, 1.00030093903294E-08_wp,&
      & 1.04916132829602E-08_wp, 7.38319650489934E-09_wp, 1.15757513617802E-08_wp,&
      & 9.05218825865469E-09_wp, 5.24168330531565E-09_wp, 7.54799324087481E-09_wp,&
      & 1.18455620254499E-08_wp, 9.04786080669313E-09_wp, 1.76524169688548E-08_wp,&
      & 2.19833234904247E-08_wp, 1.22419597732037E-08_wp, 1.38562524185172E-08_wp,&
      & 1.69021382451719E-08_wp, 1.27253555639217E-08_wp, 0.00000000000000E+00_wp,&
      & 8.54136867702080E-09_wp, 7.20236977672459E-09_wp, 1.24020779304092E-08_wp,&
      & 1.91565679449265E-08_wp, 8.69680807511640E-09_wp, 3.28263620693582E-09_wp,&
      & 1.23182242085933E-08_wp, 5.60086400495995E-09_wp, 1.23398237082115E-08_wp,&
      & 5.38978831692585E-09_wp, 1.03649522763642E-08_wp, 1.23915622937806E-08_wp,&
      & 8.81668545779098E-09_wp, 8.50339640981889E-09_wp, 1.00030093903294E-08_wp,&
      & 8.54136867702080E-09_wp, 0.00000000000000E+00_wp, 7.99592711495649E-09_wp,&
      & 4.91538331272441E-09_wp, 8.44088286698059E-09_wp, 4.90962227072692E-09_wp,&
      & 3.88992627439869E-09_wp, 1.50742955605302E-08_wp, 6.77404356579160E-09_wp,&
      & 1.50512409160924E-08_wp, 5.37826574973416E-09_wp, 8.08827272083890E-09_wp,&
      & 1.50088137018247E-08_wp, 8.62153079348212E-09_wp, 8.51283172059904E-09_wp,&
      & 1.04916132829602E-08_wp, 7.20236977672459E-09_wp, 7.99592711495649E-09_wp,&
      & 0.00000000000000E+00_wp, 5.28601527349756E-09_wp, 7.20735394414811E-09_wp,&
      & 4.85317066687809E-09_wp, 4.05611303397122E-09_wp, 6.15180823167874E-09_wp,&
      & 1.20305417727104E-08_wp, 8.49938711303868E-09_wp, 1.28271174433308E-08_wp,&
      & 1.15301681771158E-08_wp, 9.90459432279881E-09_wp, 7.26961374828290E-09_wp,&
      & 1.13688357004018E-08_wp, 7.38319650489934E-09_wp, 1.24020779304092E-08_wp,&
      & 4.91538331272441E-09_wp, 5.28601527349756E-09_wp, 0.00000000000000E+00_wp,&
      & 9.85820882850090E-09_wp, 4.48151542548455E-09_wp, 4.54964824678742E-09_wp,&
      & 8.23381234411572E-09_wp, 1.01571587419510E-08_wp, 9.90751948646844E-09_wp,&
      & 1.33434107842202E-08_wp, 1.91600879363687E-08_wp, 1.23125224768497E-08_wp,&
      & 1.19560236944380E-08_wp, 1.45313076476434E-08_wp, 1.15757513617802E-08_wp,&
      & 1.91565679449265E-08_wp, 8.44088286698059E-09_wp, 7.20735394414811E-09_wp,&
      & 9.85820882850090E-09_wp, 0.00000000000000E+00_wp, 7.13934006967499E-09_wp,&
      & 2.16169368803934E-09_wp, 7.58294130818892E-09_wp, 5.50861108840061E-09_wp,&
      & 7.81358375552646E-09_wp, 5.92505478150488E-09_wp, 8.62398908126134E-09_wp,&
      & 9.57571547166595E-09_wp, 7.65762342571572E-09_wp, 7.24825910759654E-09_wp,&
      & 9.05218825865470E-09_wp, 8.69680807511637E-09_wp, 4.90962227072692E-09_wp,&
      & 4.85317066687809E-09_wp, 4.48151542548455E-09_wp, 7.13934006967500E-09_wp,&
      & 0.00000000000000E+00_wp], shape(ref6))

   real(wp), parameter :: ref8(16, 16) = reshape([&
      & 0.00000000000000E+00_wp, 1.46572100071570E-10_wp, 8.07128937462069E-11_wp,&
      & 1.29974737294804E-10_wp, 1.18066112737297E-10_wp, 1.41700202592567E-10_wp,&
      & 1.49454425633320E-10_wp, 8.97351083058245E-11_wp, 8.80111251991700E-11_wp,&
      & 1.10129498247758E-10_wp, 1.30459271708433E-10_wp, 6.19678408978143E-11_wp,&
      & 6.04075118535431E-11_wp, 5.53641129434336E-11_wp, 9.97591835747222E-11_wp,&
      & 2.90254768950882E-11_wp, 1.46572100071570E-10_wp, 0.00000000000000E+00_wp,&
      & 2.09034766807522E-10_wp, 4.99985134487769E-10_wp, 1.37270247585517E-10_wp,&
      & 2.50132258980799E-10_wp, 4.92279951277773E-10_wp, 3.84803856992653E-10_wp,&
      & 3.10151224122767E-10_wp, 5.09696066108902E-10_wp, 2.12987361162512E-10_wp,&
      & 3.03054095820550E-10_wp, 1.26113447096229E-10_wp, 1.64407241237300E-10_wp,&
      & 2.42948259844199E-10_wp, 1.89787942124071E-10_wp, 8.07128937462069E-11_wp,&
      & 2.09034766807522E-10_wp, 0.00000000000000E+00_wp, 3.51513963193868E-10_wp,&
      & 3.60296040138773E-10_wp, 3.58512036198655E-10_wp, 3.90226543996074E-10_wp,&
      & 2.16363240549659E-10_wp, 3.12959613130810E-10_wp, 2.18949848089819E-10_wp,&
      & 3.85634111896265E-10_wp, 1.45464037399196E-10_wp, 1.85659866284909E-10_wp,&
      & 1.02425819771858E-10_wp, 3.06284030349269E-10_wp, 1.30720738543319E-10_wp,&
      & 1.29974737294804E-10_wp, 4.99985134487769E-10_wp, 3.51513963193868E-10_wp,&
      & 0.00000000000000E+00_wp, 2.17740371272162E-10_wp, 3.52968981202561E-10_wp,&
      & 4.97258602476481E-10_wp, 3.22530225824178E-10_wp, 4.09474555473078E-10_wp,&
      & 3.76981356218739E-10_wp, 2.74394585939252E-10_wp, 3.02545698969840E-10_wp,&
      & 1.29265048260295E-10_wp, 2.51016048115937E-10_wp, 3.10246257249180E-10_wp,&
      & 1.90908499659943E-10_wp, 1.18066112737297E-10_wp, 1.37270247585517E-10_wp,&
      & 3.60296040138773E-10_wp, 2.17740371272162E-10_wp, 0.00000000000000E+00_wp,&
      & 4.85752294199320E-10_wp, 2.95037022602114E-10_wp, 2.53326656457074E-10_wp,&
      & 3.41324536307942E-10_wp, 2.29245227156888E-10_wp, 4.65961397550333E-10_wp,&
      & 1.37987336832920E-10_wp, 1.37498304537020E-10_wp, 1.26908828954402E-10_wp,&
      & 3.62788661732625E-10_wp, 1.45176400521889E-10_wp, 1.41700202592567E-10_wp,&
      & 2.50132258980799E-10_wp, 3.58512036198655E-10_wp, 3.52968981202561E-10_wp,&
      & 4.85752294199320E-10_wp, 0.00000000000000E+00_wp, 4.18018215144115E-10_wp,&
      & 3.91141617773393E-10_wp, 3.78301544195624E-10_wp, 3.50649303031036E-10_wp,&
      & 4.04253608918124E-10_wp, 3.02455044439563E-10_wp, 2.36251392634850E-10_wp,&
      & 3.20812677998192E-10_wp, 1.41242232943822E-10_wp, 1.82244605530237E-10_wp,&
      & 1.49454425633320E-10_wp, 4.92279951277773E-10_wp, 3.90226543996074E-10_wp,&
      & 4.97258602476481E-10_wp, 2.95037022602114E-10_wp, 4.18018215144115E-10_wp,&
      & 0.00000000000000E+00_wp, 4.52684634235575E-10_wp, 4.34897172487277E-10_wp,&
      & 5.74854519920367E-10_wp, 4.11285070506367E-10_wp, 3.01270019980432E-10_wp,&
      & 1.34976867615242E-10_wp, 2.95274732974381E-10_wp, 3.96732200874991E-10_wp,&
      & 1.36555829762004E-10_wp, 8.97351083058245E-11_wp, 3.84803856992653E-10_wp,&
      & 2.16363240549659E-10_wp, 3.22530225824178E-10_wp, 2.53326656457074E-10_wp,&
      & 3.91141617773393E-10_wp, 4.52684634235575E-10_wp, 0.00000000000000E+00_wp,&
      & 3.32253014837556E-10_wp, 1.35038389365013E-10_wp, 4.42430141026438E-10_wp,&
      & 2.31511944993679E-10_wp, 2.33661797902900E-10_wp, 2.01419463772972E-10_wp,&
      & 3.34905448430024E-10_wp, 7.80205268001701E-11_wp, 8.80111251991700E-11_wp,&
      & 3.10151224122767E-10_wp, 3.12959613130810E-10_wp, 4.09474555473078E-10_wp,&
      & 3.41324536307942E-10_wp, 3.78301544195624E-10_wp, 4.34897172487277E-10_wp,&
      & 3.32253014837556E-10_wp, 0.00000000000000E+00_wp, 3.73505138427143E-10_wp,&
      & 3.76826889665747E-10_wp, 2.19569599347592E-10_wp, 2.24353952986883E-10_wp,&
      & 1.18380535476782E-10_wp, 1.84135122386066E-10_wp, 7.83099083524109E-11_wp,&
      & 1.10129498247758E-10_wp, 5.09696066108902E-10_wp, 2.18949848089819E-10_wp,&
      & 3.76981356218739E-10_wp, 2.29245227156888E-10_wp, 3.50649303031036E-10_wp,&
      & 5.74854519920367E-10_wp, 1.35038389365013E-10_wp, 3.73505138427143E-10_wp,&
      & 0.00000000000000E+00_wp, 4.31716037009045E-10_wp, 2.94908672116325E-10_wp,&
      & 3.08825762767101E-10_wp, 2.10157863725409E-10_wp, 3.73065178594643E-10_wp,&
      & 1.67413533725333E-10_wp, 1.30459271708433E-10_wp, 2.12987361162512E-10_wp,&
      & 3.85634111896265E-10_wp, 2.74394585939252E-10_wp, 4.65961397550333E-10_wp,&
      & 4.04253608918124E-10_wp, 4.11285070506367E-10_wp, 4.42430141026438E-10_wp,&
      & 3.76826889665747E-10_wp, 4.31716037009045E-10_wp, 0.00000000000000E+00_wp,&
      & 2.51740717966354E-10_wp, 2.03399534515101E-10_wp, 3.14014351315624E-10_wp,&
      & 1.41802650315821E-10_wp, 1.80297545292712E-10_wp, 6.19678408978145E-11_wp,&
      & 3.03054095820550E-10_wp, 1.45464037399196E-10_wp, 3.02545698969840E-10_wp,&
      & 1.37987336832920E-10_wp, 3.02455044439563E-10_wp, 3.01270019980432E-10_wp,&
      & 2.31511944993679E-10_wp, 2.19569599347592E-10_wp, 2.94908672116325E-10_wp,&
      & 2.51740717966354E-10_wp, 0.00000000000000E+00_wp, 1.08894028666122E-10_wp,&
      & 1.20583988077873E-10_wp, 2.19107060793089E-10_wp, 5.98940546991756E-11_wp,&
      & 6.04075118535431E-11_wp, 1.26113447096229E-10_wp, 1.85659866284909E-10_wp,&
      & 1.29265048260295E-10_wp, 1.37498304537020E-10_wp, 2.36251392634850E-10_wp,&
      & 1.34976867615242E-10_wp, 2.33661797902900E-10_wp, 2.24353952986883E-10_wp,&
      & 3.08825762767101E-10_wp, 2.03399534515101E-10_wp, 1.08894028666122E-10_wp,&
      & 0.00000000000000E+00_wp, 1.31527363876826E-10_wp, 1.97443536534411E-10_wp,&
      & 7.09505698650581E-11_wp, 5.53641129434336E-11_wp, 1.64407241237300E-10_wp,&
      & 1.02425819771858E-10_wp, 2.51016048115937E-10_wp, 1.26908828954402E-10_wp,&
      & 3.20812677998192E-10_wp, 2.95274732974381E-10_wp, 2.01419463772972E-10_wp,&
      & 1.18380535476782E-10_wp, 2.10157863725409E-10_wp, 3.14014351315624E-10_wp,&
      & 1.20583988077873E-10_wp, 1.31527363876826E-10_wp, 0.00000000000000E+00_wp,&
      & 2.15167970434677E-10_wp, 8.05096835379294E-11_wp, 9.97591835747222E-11_wp,&
      & 2.42948259844199E-10_wp, 3.06284030349269E-10_wp, 3.10246257249180E-10_wp,&
      & 3.62788661732625E-10_wp, 1.41242232943822E-10_wp, 3.96732200874991E-10_wp,&
      & 3.34905448430024E-10_wp, 1.84135122386066E-10_wp, 3.73065178594643E-10_wp,&
      & 1.41802650315821E-10_wp, 2.19107060793089E-10_wp, 1.97443536534411E-10_wp,&
      & 2.15167970434677E-10_wp, 0.00000000000000E+00_wp, 8.66753794405566E-11_wp,&
      & 2.90254768950882E-11_wp, 1.89787942124070E-10_wp, 1.30720738543319E-10_wp,&
      & 1.90908499659943E-10_wp, 1.45176400521889E-10_wp, 1.82244605530238E-10_wp,&
      & 1.36555829762005E-10_wp, 7.80205268001701E-11_wp, 7.83099083524098E-11_wp,&
      & 1.67413533725331E-10_wp, 1.80297545292712E-10_wp, 5.98940546991756E-11_wp,&
      & 7.09505698650581E-11_wp, 8.05096835379294E-11_wp, 8.66753794405580E-11_wp,&
      & 0.00000000000000E+00_wp], shape(ref8))

   ! PBE-D4(koide)
   param = param_type(&
      & s6 = 1.0_wp, s8 = 0.95948085_wp, a1 = 0.38574991_wp, a2 = 4.80688534_wp, &
      & rs6 = 2.0_wp, rs8 = 4.0_wp)

   call get_structure(mol, "MB16-43", "01")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%koide, -1)
   if (allocated(error)) return
   call test_damping_2b_gen(error, mol, d4, damp, param, ref6, ref8)

end subroutine test_damp_koide_2b_mb01

subroutine test_damp_koide_2b_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   real(wp), parameter :: ref6(16, 16) = reshape([&
      & 0.00000000000000E+00_wp, 1.14557698260524E-08_wp, 9.70075048806017E-09_wp,&
      & 1.57785037426664E-08_wp, 6.39690303222029E-09_wp, 2.16623201473160E-08_wp,&
      & 2.07858663083318E-08_wp, 1.04962099477579E-08_wp, 8.54598614291202E-09_wp,&
      & 1.72249792805384E-08_wp, 1.50466333279507E-08_wp, 8.44219793559008E-09_wp,&
      & 9.96260910612204E-09_wp, 1.26655429024987E-08_wp, 2.11289212936867E-08_wp,&
      & 8.30707164163668E-09_wp, 1.14557698260524E-08_wp, 0.00000000000000E+00_wp,&
      & 7.46578990490820E-09_wp, 8.15185738880176E-09_wp, 4.39026921660827E-09_wp,&
      & 1.04869283424384E-08_wp, 1.15804218583985E-08_wp, 1.13540881254013E-08_wp,&
      & 4.82786047440901E-09_wp, 9.63132581557590E-09_wp, 7.44761922255097E-09_wp,&
      & 4.92158010713194E-09_wp, 8.87459479091038E-09_wp, 8.25799441348217E-09_wp,&
      & 1.14423820389406E-08_wp, 4.91903519577719E-09_wp, 9.70075048806017E-09_wp,&
      & 7.46578990490820E-09_wp, 0.00000000000000E+00_wp, 7.01647625943397E-09_wp,&
      & 5.14780330738237E-09_wp, 8.72290919393981E-09_wp, 1.13172935999127E-08_wp,&
      & 1.50725176464964E-08_wp, 6.08322202128651E-09_wp, 8.96118548095991E-09_wp,&
      & 7.34539597165189E-09_wp, 4.48526083857271E-09_wp, 8.63258022951499E-09_wp,&
      & 1.13967763887773E-08_wp, 1.21753809581686E-08_wp, 6.03653161962977E-09_wp,&
      & 1.57785037426664E-08_wp, 8.15185738880176E-09_wp, 7.01647625943397E-09_wp,&
      & 0.00000000000000E+00_wp, 6.26287373213277E-09_wp, 1.98449922981651E-08_wp,&
      & 1.05330597287599E-08_wp, 7.56102912269640E-09_wp, 6.17737201959698E-09_wp,&
      & 1.98330869644676E-08_wp, 9.15881690426689E-09_wp, 8.41778314010628E-09_wp,&
      & 1.28846439913220E-08_wp, 8.54208015316975E-09_wp, 1.40607708344622E-08_wp,&
      & 8.24236265951067E-09_wp, 6.39690303222031E-09_wp, 4.39026921660827E-09_wp,&
      & 5.14780330738237E-09_wp, 6.26287373213277E-09_wp, 0.00000000000000E+00_wp,&
      & 7.26569385442179E-09_wp, 5.87133801598873E-09_wp, 8.19632783116026E-09_wp,&
      & 3.51016700477798E-09_wp, 8.49374769539539E-09_wp, 4.16731915813239E-09_wp,&
      & 3.43469709827163E-09_wp, 8.47005616180196E-09_wp, 6.63111007765357E-09_wp,&
      & 7.60813294984684E-09_wp, 4.59252842060985E-09_wp, 2.16623201473160E-08_wp,&
      & 1.04869283424384E-08_wp, 8.72290919393981E-09_wp, 1.98449922981651E-08_wp,&
      & 7.26569385442177E-09_wp, 0.00000000000000E+00_wp, 1.49043138755038E-08_wp,&
      & 8.94173503248868E-09_wp, 8.03495279628172E-09_wp, 2.23271840430803E-08_wp,&
      & 1.27413445746448E-08_wp, 9.98382161768481E-09_wp, 1.29862870922552E-08_wp,&
      & 1.08772452284555E-08_wp, 1.86658679555380E-08_wp, 9.65454790070431E-09_wp,&
      & 2.07858663083318E-08_wp, 1.15804218583985E-08_wp, 1.13172935999127E-08_wp,&
      & 1.05330597287599E-08_wp, 5.87133801598871E-09_wp, 1.49043138755038E-08_wp,&
      & 0.00000000000000E+00_wp, 1.39349354822391E-08_wp, 9.17382391501234E-09_wp,&
      & 1.17862941212354E-08_wp, 1.50672878483470E-08_wp, 6.44478606071434E-09_wp,&
      & 7.77364371534729E-09_wp, 1.53401142023903E-08_wp, 2.09473643110181E-08_wp,&
      & 7.29846154986097E-09_wp, 1.04962099477579E-08_wp, 1.13540881254013E-08_wp,&
      & 1.50725176464964E-08_wp, 7.56102912269640E-09_wp, 8.19632783116025E-09_wp,&
      & 8.94173503248868E-09_wp, 1.39349354822391E-08_wp, 0.00000000000000E+00_wp,&
      & 9.36847038303384E-09_wp, 9.23211669488791E-09_wp, 9.88342047710461E-09_wp,&
      & 6.17053434225762E-09_wp, 1.06615368178969E-08_wp, 1.42911111337430E-08_wp,&
      & 1.40522371412051E-08_wp, 7.97989559175407E-09_wp, 8.54598614291202E-09_wp,&
      & 4.82786047440901E-09_wp, 6.08322202128651E-09_wp, 6.17737201959698E-09_wp,&
      & 3.51016700477798E-09_wp, 8.03495279628172E-09_wp, 9.17382391501234E-09_wp,&
      & 9.36847038303384E-09_wp, 0.00000000000000E+00_wp, 8.55862465803331E-09_wp,&
      & 5.66596183846194E-09_wp, 3.08336985162779E-09_wp, 6.03313367940650E-09_wp,&
      & 1.14589185958898E-08_wp, 1.12932362907123E-08_wp, 5.36464812937315E-09_wp,&
      & 1.72249792805384E-08_wp, 9.63132581557590E-09_wp, 8.96118548095991E-09_wp,&
      & 1.98330869644676E-08_wp, 8.49374769539538E-09_wp, 2.23271840430803E-08_wp,&
      & 1.17862941212354E-08_wp, 9.23211669488791E-09_wp, 8.55862465803331E-09_wp,&
      & 0.00000000000000E+00_wp, 1.08846085609738E-08_wp, 9.98858446482475E-09_wp,&
      & 1.53486420438740E-08_wp, 1.16227065225159E-08_wp, 1.74413363711312E-08_wp,&
      & 1.19099461488647E-08_wp, 1.50466333279507E-08_wp, 7.44761922255097E-09_wp,&
      & 7.34539597165189E-09_wp, 9.15881690426689E-09_wp, 4.16731915813239E-09_wp,&
      & 1.27413445746448E-08_wp, 1.50672878483470E-08_wp, 9.88342047710461E-09_wp,&
      & 5.66596183846194E-09_wp, 1.08846085609738E-08_wp, 0.00000000000000E+00_wp,&
      & 4.78648756225713E-09_wp, 7.21193594816660E-09_wp, 1.06770201157759E-08_wp,&
      & 1.48836665753766E-08_wp, 5.53428623624055E-09_wp, 8.44219793559008E-09_wp,&
      & 4.92158010713194E-09_wp, 4.48526083857271E-09_wp, 8.41778314010628E-09_wp,&
      & 3.43469709827163E-09_wp, 9.98382161768481E-09_wp, 6.44478606071434E-09_wp,&
      & 6.17053434225762E-09_wp, 3.08336985162779E-09_wp, 9.98858446482475E-09_wp,&
      & 4.78648756225713E-09_wp, 0.00000000000000E+00_wp, 9.24956842672786E-09_wp,&
      & 4.87647804434493E-09_wp, 7.55592728560088E-09_wp, 3.98004523579369E-09_wp,&
      & 9.96260910612204E-09_wp, 8.87459479091038E-09_wp, 8.63258022951499E-09_wp,&
      & 1.28846439913220E-08_wp, 8.47005616180196E-09_wp, 1.29862870922552E-08_wp,&
      & 7.77364371534729E-09_wp, 1.06615368178969E-08_wp, 6.03313367940650E-09_wp,&
      & 1.53486420438740E-08_wp, 7.21193594816660E-09_wp, 9.24956842672786E-09_wp,&
      & 0.00000000000000E+00_wp, 7.27960859052478E-09_wp, 1.02595060014340E-08_wp,&
      & 8.20181641529113E-09_wp, 1.26655429024987E-08_wp, 8.25799441348217E-09_wp,&
      & 1.13967763887773E-08_wp, 8.54208015316975E-09_wp, 6.63111007765355E-09_wp,&
      & 1.08772452284555E-08_wp, 1.53401142023903E-08_wp, 1.42911111337430E-08_wp,&
      & 1.14589185958898E-08_wp, 1.16227065225159E-08_wp, 1.06770201157759E-08_wp,&
      & 4.87647804434493E-09_wp, 7.27960859052478E-09_wp, 0.00000000000000E+00_wp,&
      & 1.97506658721727E-08_wp, 1.08104973218505E-08_wp, 2.11289212936867E-08_wp,&
      & 1.14423820389406E-08_wp, 1.21753809581686E-08_wp, 1.40607708344622E-08_wp,&
      & 7.60813294984684E-09_wp, 1.86658679555380E-08_wp, 2.09473643110181E-08_wp,&
      & 1.40522371412051E-08_wp, 1.12932362907123E-08_wp, 1.74413363711312E-08_wp,&
      & 1.48836665753766E-08_wp, 7.55592728560088E-09_wp, 1.02595060014340E-08_wp,&
      & 1.97506658721727E-08_wp, 0.00000000000000E+00_wp, 1.10669251531723E-08_wp,&
      & 8.30707164163668E-09_wp, 4.91903519577719E-09_wp, 6.03653161962977E-09_wp,&
      & 8.24236265951067E-09_wp, 4.59252842060985E-09_wp, 9.65454790070431E-09_wp,&
      & 7.29846154986097E-09_wp, 7.97989559175407E-09_wp, 5.36464812937315E-09_wp,&
      & 1.19099461488647E-08_wp, 5.53428623624055E-09_wp, 3.98004523579369E-09_wp,&
      & 8.20181641529113E-09_wp, 1.08104973218505E-08_wp, 1.10669251531723E-08_wp,&
      & 0.00000000000000E+00_wp], shape(ref6))

   real(wp), parameter :: ref8(16, 16) = reshape([&
      & 0.00000000000000E+00_wp, 2.78344428280108E-10_wp, 2.89726383032657E-10_wp,&
      & 4.51020348885246E-10_wp, 1.69467730525696E-10_wp, 4.34708231461063E-10_wp,&
      & 5.00205174253056E-10_wp, 3.36012428939871E-10_wp, 2.26283092397584E-10_wp,&
      & 5.70811012721970E-10_wp, 1.29890858690800E-10_wp, 2.19853096367033E-10_wp,&
      & 3.13470752045924E-10_wp, 4.29209114091478E-10_wp, 4.77480231297244E-10_wp,&
      & 2.40648301985472E-10_wp, 2.78344428280108E-10_wp, 0.00000000000000E+00_wp,&
      & 9.01155856417038E-11_wp, 2.09536874810803E-10_wp, 6.65698410229814E-11_wp,&
      & 2.85667667198120E-10_wp, 2.75506601342228E-10_wp, 2.80303090928473E-10_wp,&
      & 8.39332823527320E-11_wp, 2.75398308855002E-10_wp, 9.14650145895225E-11_wp,&
      & 7.33598571882231E-11_wp, 2.32419659044975E-10_wp, 2.39080167932265E-10_wp,&
      & 2.78620296537302E-10_wp, 1.12311549065646E-10_wp, 2.89726383032657E-10_wp,&
      & 9.01155856417038E-11_wp, 0.00000000000000E+00_wp, 1.93527757565139E-10_wp,&
      & 7.04316066634236E-11_wp, 2.58806874249903E-10_wp, 3.19884727074113E-10_wp,&
      & 1.26357696187031E-10_wp, 6.91941906073247E-11_wp, 2.66871186674355E-10_wp,&
      & 1.49212204809146E-10_wp, 9.68899445398839E-11_wp, 2.44736581134180E-10_wp,&
      & 3.20335298054312E-10_wp, 3.17602145311148E-10_wp, 1.34857440790375E-10_wp,&
      & 4.51020348885246E-10_wp, 2.09536874810803E-10_wp, 1.93527757565139E-10_wp,&
      & 0.00000000000000E+00_wp, 1.42233998084950E-10_wp, 1.38873021028707E-10_wp,&
      & 3.35914572507722E-10_wp, 2.15952143297426E-10_wp, 1.55468216752505E-10_wp,&
      & 1.40821826334697E-10_wp, 2.38953246383922E-10_wp, 1.10515058914088E-10_wp,&
      & 3.81374082605421E-10_wp, 2.55175217391481E-10_wp, 4.46021188896738E-10_wp,&
      & 2.10197735585964E-10_wp, 1.69467730525696E-10_wp, 6.65698410229814E-11_wp,&
      & 7.04316066634236E-11_wp, 1.42233998084950E-10_wp, 0.00000000000000E+00_wp,&
      & 1.90502201382386E-10_wp, 1.53546047031894E-10_wp, 1.99365115504022E-10_wp,&
      & 4.82106026674320E-11_wp, 1.97627495922231E-10_wp, 8.60776718299949E-11_wp,&
      & 4.45558283750831E-11_wp, 1.04821988398972E-10_wp, 1.75943507510684E-10_wp,&
      & 1.95874812942583E-10_wp, 5.86045621078018E-11_wp, 4.34708231461063E-10_wp,&
      & 2.85667667198120E-10_wp, 2.58806874249903E-10_wp, 1.38873021028707E-10_wp,&
      & 1.90502201382386E-10_wp, 0.00000000000000E+00_wp, 5.15379309109629E-10_wp,&
      & 2.69962146443269E-10_wp, 2.18327839027037E-10_wp, 3.67392027615542E-10_wp,&
      & 3.05891985396770E-10_wp, 1.99949932100771E-10_wp, 4.31863671575587E-10_wp,&
      & 3.52467457942993E-10_wp, 5.72433653748293E-10_wp, 2.75833881592673E-10_wp,&
      & 5.00205174253056E-10_wp, 2.75506601342228E-10_wp, 3.19884727074113E-10_wp,&
      & 3.35914572507722E-10_wp, 1.53546047031894E-10_wp, 5.15379309109629E-10_wp,&
      & 0.00000000000000E+00_wp, 4.80455239428834E-10_wp, 2.28860702168893E-10_wp,&
      & 3.91755494565176E-10_wp, 1.27074993704338E-10_wp, 1.73532046039781E-10_wp,&
      & 2.23745229551497E-10_wp, 5.29319245227162E-10_wp, 4.89952044265243E-10_wp,&
      & 2.06268156173653E-10_wp, 3.36012428939871E-10_wp, 2.80303090928473E-10_wp,&
      & 1.26357696187031E-10_wp, 2.15952143297426E-10_wp, 1.99365115504021E-10_wp,&
      & 2.69962146443269E-10_wp, 4.80455239428834E-10_wp, 0.00000000000000E+00_wp,&
      & 2.27634361600183E-10_wp, 2.82120608400972E-10_wp, 2.94716467604735E-10_wp,&
      & 1.64579674329182E-10_wp, 3.42460444981896E-10_wp, 4.93828683394407E-10_wp,&
      & 4.84920582708860E-10_wp, 2.29947570146549E-10_wp, 2.26283092397584E-10_wp,&
      & 8.39332823527320E-11_wp, 6.91941906073247E-11_wp, 1.55468216752505E-10_wp,&
      & 4.82106026674320E-11_wp, 2.18327839027037E-10_wp, 2.28860702168893E-10_wp,&
      & 2.27634361600183E-10_wp, 0.00000000000000E+00_wp, 2.26420654502789E-10_wp,&
      & 9.09999821736905E-11_wp, 6.02455942843927E-11_wp, 1.56733467068586E-10_wp,&
      & 1.19119310306544E-10_wp, 1.36144453894978E-10_wp, 6.24640417814111E-11_wp,&
      & 5.70811012721970E-10_wp, 2.75398308855002E-10_wp, 2.66871186674355E-10_wp,&
      & 1.40821826334697E-10_wp, 1.97627495922231E-10_wp, 3.67392027615542E-10_wp,&
      & 3.91755494565176E-10_wp, 2.82120608400972E-10_wp, 2.26420654502789E-10_wp,&
      & 0.00000000000000E+00_wp, 3.15437060171633E-10_wp, 1.99757265450206E-10_wp,&
      & 4.88340304027799E-10_wp, 3.84703736226700E-10_wp, 5.73042061408769E-10_wp,&
      & 2.65523961585653E-10_wp, 1.29890858690800E-10_wp, 9.14650145895225E-11_wp,&
      & 1.49212204809146E-10_wp, 2.38953246383922E-10_wp, 8.60776718299949E-11_wp,&
      & 3.05891985396770E-10_wp, 1.27074993704338E-10_wp, 2.94716467604735E-10_wp,&
      & 9.09999821736905E-11_wp, 3.15437060171633E-10_wp, 0.00000000000000E+00_wp,&
      & 9.86667576341274E-11_wp, 2.02092497233031E-10_wp, 3.12221788917944E-10_wp,&
      & 1.51172047834594E-10_wp, 1.30670204211734E-10_wp, 2.19853096367033E-10_wp,&
      & 7.33598571882231E-11_wp, 9.68899445398839E-11_wp, 1.10515058914088E-10_wp,&
      & 4.45558283750831E-11_wp, 1.99949932100771E-10_wp, 1.73532046039781E-10_wp,&
      & 1.64579674329182E-10_wp, 6.02455942843927E-11_wp, 1.99757265450206E-10_wp,&
      & 9.86667576341274E-11_wp, 0.00000000000000E+00_wp, 1.13801092316488E-10_wp,&
      & 1.20785787282092E-10_wp, 2.05115086500754E-10_wp, 8.27621982768348E-11_wp,&
      & 3.13470752045924E-10_wp, 2.32419659044975E-10_wp, 2.44736581134180E-10_wp,&
      & 3.81374082605421E-10_wp, 1.04821988398972E-10_wp, 4.31863671575587E-10_wp,&
      & 2.23745229551497E-10_wp, 3.42460444981896E-10_wp, 1.56733467068586E-10_wp,&
      & 4.88340304027799E-10_wp, 2.02092497233031E-10_wp, 1.13801092316488E-10_wp,&
      & 0.00000000000000E+00_wp, 2.04202221386513E-10_wp, 3.25814715217996E-10_wp,&
      & 2.23485290188131E-10_wp, 4.29209114091478E-10_wp, 2.39080167932265E-10_wp,&
      & 3.20335298054312E-10_wp, 2.55175217391481E-10_wp, 1.75943507510684E-10_wp,&
      & 3.52467457942993E-10_wp, 5.29319245227162E-10_wp, 4.93828683394407E-10_wp,&
      & 1.19119310306544E-10_wp, 3.84703736226700E-10_wp, 3.12221788917944E-10_wp,&
      & 1.20785787282092E-10_wp, 2.04202221386513E-10_wp, 0.00000000000000E+00_wp,&
      & 5.48769063081561E-10_wp, 2.85843658976709E-10_wp, 4.77480231297244E-10_wp,&
      & 2.78620296537302E-10_wp, 3.17602145311148E-10_wp, 4.46021188896738E-10_wp,&
      & 1.95874812942584E-10_wp, 5.72433653748293E-10_wp, 4.89952044265243E-10_wp,&
      & 4.84920582708860E-10_wp, 1.36144453894978E-10_wp, 5.73042061408769E-10_wp,&
      & 1.51172047834594E-10_wp, 2.05115086500754E-10_wp, 3.25814715217996E-10_wp,&
      & 5.48769063081561E-10_wp, 0.00000000000000E+00_wp, 2.84213803587800E-10_wp,&
      & 2.40648301985472E-10_wp, 1.12311549065646E-10_wp, 1.34857440790375E-10_wp,&
      & 2.10197735585964E-10_wp, 5.86045621078018E-11_wp, 2.75833881592673E-10_wp,&
      & 2.06268156173653E-10_wp, 2.29947570146549E-10_wp, 6.24640417814111E-11_wp,&
      & 2.65523961585653E-10_wp, 1.30670204211734E-10_wp, 8.27621982768348E-11_wp,&
      & 2.23485290188131E-10_wp, 2.85843658976709E-10_wp, 2.84213803587800E-10_wp,&
      & 0.00000000000000E+00_wp], shape(ref8))

   ! PBE-D4(koide)
   param = param_type(&
      & s6 = 1.0_wp, s8 = 0.95948085_wp, a1 = 0.38574991_wp, a2 = 4.80688534_wp, &
      & rs6 = 2.0_wp, rs8 = 4.0_wp)

   call get_structure(mol, "MB16-43", "02")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%koide, -1)
   if (allocated(error)) return
   call test_damping_2b_gen(error, mol, d4, damp, param, ref6, ref8)

end subroutine test_damp_koide_2b_mb02

subroutine test_grad_koide_2b_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   ! PBE-D4(koide)
   param = param_type(&
      & s6 = 1.0_wp, s8 = 0.95948085_wp, a1 = 0.38574991_wp, a2 = 4.80688534_wp, &
      & rs6 = 2.0_wp, rs8 = 4.0_wp)

   call get_structure(mol, "MB16-43", "03")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%koide, -1)
   if (allocated(error)) return
   call test_damping_2b_numgrad(error, mol, d4, damp, param)

end subroutine test_grad_koide_2b_mb03

subroutine test_grad_koide_2b_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   ! PBE-D4(koide)
   param = param_type(&
      & s6 = 1.0_wp, s8 = 0.95948085_wp, a1 = 0.38574991_wp, a2 = 4.80688534_wp, &
      & rs6 = 2.0_wp, rs8 = 4.0_wp)

   call get_structure(mol, "MB16-43", "04")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%koide, -1)
   if (allocated(error)) return
   call test_damping_2b_numgrad(error, mol, d4, damp, param)

end subroutine test_grad_koide_2b_mb04


subroutine test_damp_rational_3b_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   real(wp), parameter :: ref9(16, 16) = reshape([&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 5.89753787094091E-02_wp, 1.88509729040082E-02_wp, 2.47440763557434E-01_wp,&
      & 2.60462781340476E-01_wp, 2.83454922685748E-02_wp, 2.17355957264545E-01_wp,&
      & 7.43733114593417E-02_wp, 2.29367617020142E-01_wp, 2.96400734393158E-01_wp,&
      & 6.90954191410616E-02_wp, 5.13583677975919E-03_wp, 8.31934571212963E-02_wp,&
      & 2.06480742624636E-01_wp, 6.60390065726004E-02_wp, 0.00000000000000E+00_wp,&
      & 5.89753787094091E-02_wp, 0.00000000000000E+00_wp, 1.56237856384890E-02_wp,&
      & 1.73654140684522E-02_wp, 5.02312554691199E-02_wp, 2.23963455475711E-02_wp,&
      & 6.43698151016757E-02_wp, 6.40628231580018E-03_wp, 7.73352767274106E-02_wp,&
      & 5.20272950560371E-02_wp, 4.33551672924808E-02_wp, 2.22874714754533E-02_wp,&
      & 1.10282722990614E-03_wp, 3.60240869579873E-02_wp, 1.85262353488588E-02_wp,&
      & 0.00000000000000E+00_wp, 1.88509729040082E-02_wp, 1.56237856384890E-02_wp,&
      & 0.00000000000000E+00_wp, 7.51214015600882E-02_wp, 7.69354070105746E-02_wp,&
      & 9.49555432938421E-03_wp, 8.06903787867501E-02_wp, 1.82876290009549E-02_wp,&
      & 9.34303882735737E-02_wp, 9.18123423411253E-02_wp, 2.23057117042037E-02_wp,&
      & 1.74710250243169E-03_wp, 2.24341249006189E-02_wp, 6.06582283391375E-02_wp,&
      & 1.97872549301519E-02_wp, 0.00000000000000E+00_wp, 2.47440763557434E-01_wp,&
      & 1.73654140684522E-02_wp, 7.51214015600882E-02_wp, 0.00000000000000E+00_wp,&
      & 1.34476756534431E-01_wp, 1.08481741471873E-01_wp, 2.42348562825782E-01_wp,&
      & 2.47909831594960E-02_wp, 3.02451383868491E-01_wp, 8.57991583092925E-02_wp,&
      & 1.77286072762406E-01_wp, 1.03293574187133E-01_wp, 5.30771830581002E-03_wp,&
      & 8.48865044319875E-02_wp, 7.35437444332664E-02_wp, 0.00000000000000E+00_wp,&
      & 2.60462781340476E-01_wp, 5.02312554691199E-02_wp, 7.69354070105746E-02_wp,&
      & 1.34476756534431E-01_wp, 0.00000000000000E+00_wp, 1.11581426451101E-01_wp,&
      & 2.26829506886783E-01_wp, 2.19871127024973E-02_wp, 3.14089331854819E-01_wp,&
      & 4.25286380228427E-02_wp, 1.24990097518577E-01_wp, 1.00728286022353E-01_wp,&
      & 4.41616567642880E-02_wp, 7.50139357546020E-03_wp, 4.38304178752185E-02_wp,&
      & 0.00000000000000E+00_wp, 2.83454922685748E-02_wp, 2.23963455475711E-02_wp,&
      & 9.49555432938421E-03_wp, 1.08481741471873E-01_wp, 1.11581426451101E-01_wp,&
      & 0.00000000000000E+00_wp, 6.24247464538160E-02_wp, 1.95649829628904E-02_wp,&
      & 7.33106401469414E-02_wp, 1.24491950518668E-01_wp, 3.43624005731382E-02_wp,&
      & 2.96846399516899E-03_wp, 2.94674692779246E-02_wp, 7.63292842426705E-02_wp,&
      & 9.38110546827499E-03_wp, 0.00000000000000E+00_wp, 2.17355957264545E-01_wp,&
      & 6.43698151016757E-02_wp, 8.06903787867501E-02_wp, 2.42348562825782E-01_wp,&
      & 2.26829506886783E-01_wp, 6.24247464538160E-02_wp, 0.00000000000000E+00_wp,&
      & 5.53677493483528E-02_wp, 1.06880044881108E-02_wp, 2.00136163695598E-01_wp,&
      & 1.22257769057687E-01_wp, 7.81704383519018E-02_wp, 7.70210248289603E-02_wp,&
      & 1.24395744336467E-01_wp, 1.37841234063769E-02_wp, 0.00000000000000E+00_wp,&
      & 7.43733114593417E-02_wp, 6.40628231580018E-03_wp, 1.82876290009549E-02_wp,&
      & 2.47909831594960E-02_wp, 2.19871127024973E-02_wp, 1.95649829628904E-02_wp,&
      & 5.53677493483528E-02_wp, 0.00000000000000E+00_wp, 8.54169320890958E-02_wp,&
      & 2.40289468420661E-02_wp, 3.63379070057053E-02_wp, 2.20765216755316E-02_wp,&
      & 2.34822201004519E-03_wp, 7.26755306384071E-03_wp, 4.73611397864343E-03_wp,&
      & 0.00000000000000E+00_wp, 2.29367617020142E-01_wp, 7.73352767274106E-02_wp,&
      & 9.34303882735737E-02_wp, 3.02451383868491E-01_wp, 3.14089331854819E-01_wp,&
      & 7.33106401469414E-02_wp, 1.06880044881108E-02_wp, 8.54169320890958E-02_wp,&
      & 0.00000000000000E+00_wp, 3.13145650133690E-01_wp, 1.66595360615635E-01_wp,&
      & 9.13123504133823E-02_wp, 1.01221231389186E-01_wp, 2.14320265897356E-01_wp,&
      & 4.07149514408583E-02_wp, 0.00000000000000E+00_wp, 2.96400734393158E-01_wp,&
      & 5.20272950560371E-02_wp, 9.18123423411253E-02_wp, 8.57991583092925E-02_wp,&
      & 4.25286380228427E-02_wp, 1.24491950518668E-01_wp, 2.00136163695598E-01_wp,&
      & 2.40289468420661E-02_wp, 3.13145650133690E-01_wp, 0.00000000000000E+00_wp,&
      & 1.80945595568419E-01_wp, 1.20077883488168E-01_wp, 3.74902347224943E-02_wp,&
      & 8.34797576626734E-03_wp, 4.63339935874524E-02_wp, 0.00000000000000E+00_wp,&
      & 6.90954191410616E-02_wp, 4.33551672924808E-02_wp, 2.23057117042037E-02_wp,&
      & 1.77286072762406E-01_wp, 1.24990097518577E-01_wp, 3.43624005731382E-02_wp,&
      & 1.22257769057687E-01_wp, 3.63379070057053E-02_wp, 1.66595360615635E-01_wp,&
      & 1.80945595568419E-01_wp, 0.00000000000000E+00_wp, 1.36930631831915E-02_wp,&
      & 5.45785858789404E-02_wp, 9.28282463222766E-02_wp, 1.66751298668331E-02_wp,&
      & 0.00000000000000E+00_wp, 5.13583677975919E-03_wp, 2.22874714754533E-02_wp,&
      & 1.74710250243169E-03_wp, 1.03293574187133E-01_wp, 1.00728286022353E-01_wp,&
      & 2.96846399516899E-03_wp, 7.81704383519018E-02_wp, 2.20765216755316E-02_wp,&
      & 9.13123504133823E-02_wp, 1.20077883488168E-01_wp, 1.36930631831915E-02_wp,&
      & 0.00000000000000E+00_wp, 2.99475552989371E-02_wp, 7.38738448274877E-02_wp,&
      & 1.35794448957664E-02_wp, 0.00000000000000E+00_wp, 8.31934571212963E-02_wp,&
      & 1.10282722990614E-03_wp, 2.24341249006189E-02_wp, 5.30771830581002E-03_wp,&
      & 4.41616567642880E-02_wp, 2.94674692779246E-02_wp, 7.70210248289603E-02_wp,&
      & 2.34822201004519E-03_wp, 1.01221231389186E-01_wp, 3.74902347224943E-02_wp,&
      & 5.45785858789404E-02_wp, 2.99475552989371E-02_wp, 0.00000000000000E+00_wp,&
      & 2.46950229009068E-02_wp, 1.63421353384773E-02_wp, 0.00000000000000E+00_wp,&
      & 2.06480742624636E-01_wp, 3.60240869579873E-02_wp, 6.06582283391375E-02_wp,&
      & 8.48865044319875E-02_wp, 7.50139357546020E-03_wp, 7.63292842426705E-02_wp,&
      & 1.24395744336467E-01_wp, 7.26755306384071E-03_wp, 2.14320265897356E-01_wp,&
      & 8.34797576626734E-03_wp, 9.28282463222766E-02_wp, 7.38738448274877E-02_wp,&
      & 2.46950229009068E-02_wp, 0.00000000000000E+00_wp, 1.46977909201259E-02_wp,&
      & 0.00000000000000E+00_wp, 6.60390065726004E-02_wp, 1.85262353488588E-02_wp,&
      & 1.97872549301519E-02_wp, 7.35437444332664E-02_wp, 4.38304178752185E-02_wp,&
      & 9.38110546827499E-03_wp, 1.37841234063769E-02_wp, 4.73611397864343E-03_wp,&
      & 4.07149514408583E-02_wp, 4.63339935874524E-02_wp, 1.66751298668331E-02_wp,&
      & 1.35794448957664E-02_wp, 1.63421353384773E-02_wp, 1.46977909201259E-02_wp,&
      & 0.00000000000000E+00_wp], shape(ref9))

   ! PBE-D4-ATM(BJ)
   param = param_type(&
      & s9=1.0, a1 = 0.38574991_wp, a2 = 4.80688534_wp)

   call get_structure(mol, "MB16-43", "01")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%rational, &
      & threebody_damping_function%rational)
   if (allocated(error)) return
   call test_damping_3b_gen(error, mol, d4, damp, param, ref9)

end subroutine test_damp_rational_3b_mb01

subroutine test_damp_rational_3b_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   real(wp), parameter :: ref9(16, 16) = reshape([&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 1.82231767176402E-02_wp, 4.50616293785077E-02_wp, 4.15159296277611E-02_wp,&
      & 1.72245217251960E-02_wp, 1.71031106459543E-02_wp, 7.32179196354087E-02_wp,&
      & 3.76364759089327E-02_wp, 6.80844562142755E-02_wp, 1.29906017154505E-03_wp,&
      & 2.65762034385782E-02_wp, 9.74853310766595E-02_wp, 1.31926206810714E-01_wp,&
      & 1.59033306259119E-02_wp, 9.55225367259856E-02_wp, 0.00000000000000E+00_wp,&
      & 1.82231767176402E-02_wp, 0.00000000000000E+00_wp, 1.27629169809255E-01_wp,&
      & 5.91862194326416E-02_wp, 5.21144850650301E-02_wp, 4.88845916339636E-02_wp,&
      & 1.70020318717953E-02_wp, 2.81516180274396E-02_wp, 1.64060683213572E-01_wp,&
      & 6.48876769036187E-03_wp, 1.20598615127684E-01_wp, 2.39068528713025E-01_wp,&
      & 1.67044841586299E-01_wp, 3.37947429211108E-02_wp, 1.45392390504375E-01_wp,&
      & 0.00000000000000E+00_wp, 4.50616293785077E-02_wp, 1.27629169809255E-01_wp,&
      & 0.00000000000000E+00_wp, 8.53611294399357E-02_wp, 1.13991693285535E-03_wp,&
      & 4.69625269971565E-02_wp, 2.17160986330350E-01_wp, 8.35096812244028E-02_wp,&
      & 3.77712081424424E-03_wp, 5.62362430308576E-03_wp, 1.70486229902428E-02_wp,&
      & 1.09306154233064E-01_wp, 1.84847261000541E-01_wp, 2.81686165882178E-02_wp,&
      & 9.14961134964628E-02_wp, 0.00000000000000E+00_wp, 4.15159296277611E-02_wp,&
      & 5.91862194326416E-02_wp, 8.53611294399357E-02_wp, 0.00000000000000E+00_wp,&
      & 4.44473927655016E-02_wp, 8.38222514837560E-02_wp, 1.86193983984827E-01_wp,&
      & 6.02168452299122E-02_wp, 9.44258005425578E-02_wp, 1.10897025025922E-02_wp,&
      & 5.30000496257445E-02_wp, 4.75993611830862E-02_wp, 2.56758112718652E-01_wp,&
      & 5.11622397712132E-02_wp, 6.08302785089858E-02_wp, 0.00000000000000E+00_wp,&
      & 1.72245217251960E-02_wp, 5.21144850650301E-02_wp, 1.13991693285535E-03_wp,&
      & 4.44473927655016E-02_wp, 0.00000000000000E+00_wp, 1.60408557909856E-02_wp,&
      & 8.87194909980438E-02_wp, 3.28944158645641E-02_wp, 5.87602055043863E-03_wp,&
      & 1.50139277293597E-03_wp, 1.36853141706912E-02_wp, 6.22527946700777E-02_wp,&
      & 7.30097996865720E-02_wp, 8.43921307246712E-03_wp, 4.37681556660063E-02_wp,&
      & 0.00000000000000E+00_wp, 1.71031106459543E-02_wp, 4.88845916339636E-02_wp,&
      & 4.69625269971565E-02_wp, 8.38222514837559E-02_wp, 1.60408557909856E-02_wp,&
      & 0.00000000000000E+00_wp, 9.78983053352314E-02_wp, 3.40649555387963E-02_wp,&
      & 6.47543842796438E-02_wp, 3.03381626842943E-04_wp, 5.89452045420041E-02_wp,&
      & 1.24346344569446E-01_wp, 7.71451770262709E-02_wp, 6.50517777861042E-03_wp,&
      & 8.83449833764136E-02_wp, 0.00000000000000E+00_wp, 7.32179196354087E-02_wp,&
      & 1.70020318717953E-02_wp, 2.17160986330350E-01_wp, 1.86193983984827E-01_wp,&
      & 8.87194909980438E-02_wp, 9.78983053352314E-02_wp, 0.00000000000000E+00_wp,&
      & 1.23550547426184E-01_wp, 2.84026568950313E-01_wp, 1.73782675818426E-02_wp,&
      & 2.38934868118472E-01_wp, 4.20574303803425E-01_wp, 3.31504776306283E-01_wp,&
      & 8.52327810038826E-02_wp, 3.20304349276087E-01_wp, 0.00000000000000E+00_wp,&
      & 3.76364759089327E-02_wp, 2.81516180274396E-02_wp, 8.35096812244028E-02_wp,&
      & 6.02168452299122E-02_wp, 3.28944158645641E-02_wp, 3.40649555387963E-02_wp,&
      & 1.23550547426184E-01_wp, 0.00000000000000E+00_wp, 9.32769096333806E-02_wp,&
      & 3.83129130638591E-03_wp, 1.02350004748365E-01_wp, 2.13928208114226E-01_wp,&
      & 2.00536822535515E-02_wp, 6.57430399490490E-03_wp, 3.44649979295415E-02_wp,&
      & 0.00000000000000E+00_wp, 6.80844562142755E-02_wp, 1.64060683213572E-01_wp,&
      & 3.77712081424424E-03_wp, 9.44258005425577E-02_wp, 5.87602055043863E-03_wp,&
      & 6.47543842796438E-02_wp, 2.84026568950313E-01_wp, 9.32769096333806E-02_wp,&
      & 0.00000000000000E+00_wp, 8.35733134464843E-03_wp, 4.41434262504296E-02_wp,&
      & 1.49184212688693E-01_wp, 2.27577268496875E-01_wp, 3.37854338793665E-02_wp,&
      & 6.91301943625666E-02_wp, 0.00000000000000E+00_wp, 1.29906017154505E-03_wp,&
      & 6.48876769036187E-03_wp, 5.62362430308576E-03_wp, 1.10897025025922E-02_wp,&
      & 1.50139277293597E-03_wp, 3.03381626842943E-04_wp, 1.73782675818426E-02_wp,&
      & 3.83129130638591E-03_wp, 8.35733134464843E-03_wp, 0.00000000000000E+00_wp,&
      & 7.02201025513863E-03_wp, 2.03526214132487E-02_wp, 1.35025979369222E-02_wp,&
      & 3.65478437144764E-04_wp, 1.24405639169578E-02_wp, 0.00000000000000E+00_wp,&
      & 2.65762034385782E-02_wp, 1.20598615127684E-01_wp, 1.70486229902428E-02_wp,&
      & 5.30000496257445E-02_wp, 1.36853141706912E-02_wp, 5.89452045420041E-02_wp,&
      & 2.38934868118472E-01_wp, 1.02350004748365E-01_wp, 4.41434262504296E-02_wp,&
      & 7.02201025513863E-03_wp, 0.00000000000000E+00_wp, 3.18629558322787E-02_wp,&
      & 2.43772483299753E-01_wp, 4.28292098380329E-02_wp, 1.32045445605594E-01_wp,&
      & 0.00000000000000E+00_wp, 9.74853310766595E-02_wp, 2.39068528713025E-01_wp,&
      & 1.09306154233064E-01_wp, 4.75993611830862E-02_wp, 6.22527946700777E-02_wp,&
      & 1.24346344569446E-01_wp, 4.20574303803425E-01_wp, 2.13928208114226E-01_wp,&
      & 1.49184212688693E-01_wp, 2.03526214132487E-02_wp, 3.18629558322787E-02_wp,&
      & 0.00000000000000E+00_wp, 4.41399639656738E-01_wp, 9.69027712477043E-02_wp,&
      & 2.42931379540518E-01_wp, 0.00000000000000E+00_wp, 1.31926206810714E-01_wp,&
      & 1.67044841586299E-01_wp, 1.84847261000541E-01_wp, 2.56758112718652E-01_wp,&
      & 7.30097996865720E-02_wp, 7.71451770262709E-02_wp, 3.31504776306283E-01_wp,&
      & 2.00536822535515E-02_wp, 2.27577268496875E-01_wp, 1.35025979369222E-02_wp,&
      & 2.43772483299753E-01_wp, 4.41399639656738E-01_wp, 0.00000000000000E+00_wp,&
      & 3.30009197087160E-02_wp, 1.64373907317028E-01_wp, 0.00000000000000E+00_wp,&
      & 1.59033306259119E-02_wp, 3.37947429211108E-02_wp, 2.81686165882178E-02_wp,&
      & 5.11622397712132E-02_wp, 8.43921307246712E-03_wp, 6.50517777861042E-03_wp,&
      & 8.52327810038826E-02_wp, 6.57430399490490E-03_wp, 3.37854338793665E-02_wp,&
      & 3.65478437144764E-04_wp, 4.28292098380329E-02_wp, 9.69027712477043E-02_wp,&
      & 3.30009197087160E-02_wp, 0.00000000000000E+00_wp, 3.81595405263840E-02_wp,&
      & 0.00000000000000E+00_wp, 9.55225367259856E-02_wp, 1.45392390504375E-01_wp,&
      & 9.14961134964628E-02_wp, 6.08302785089858E-02_wp, 4.37681556660063E-02_wp,&
      & 8.83449833764136E-02_wp, 3.20304349276087E-01_wp, 3.44649979295415E-02_wp,&
      & 6.91301943625666E-02_wp, 1.24405639169578E-02_wp, 1.32045445605594E-01_wp,&
      & 2.42931379540518E-01_wp, 1.64373907317028E-01_wp, 3.81595405263840E-02_wp,&
      & 0.00000000000000E+00_wp], shape(ref9))

   ! PBE-D4-ATM(BJ)
   param = param_type(&
      & s9=1.0, a1 = 0.38574991_wp, a2 = 4.80688534_wp)

   call get_structure(mol, "MB16-43", "02")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%rational, &
      & threebody_damping_function%rational)
   if (allocated(error)) return
   call test_damping_3b_gen(error, mol, d4, damp, param, ref9)

end subroutine test_damp_rational_3b_mb02

subroutine test_grad_rational_3b_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   ! PBE-D4-ATM(BJ)
   param = param_type(&
      & s9=1.0, a1 = 0.38574991_wp, a2 = 4.80688534_wp)

   call get_structure(mol, "MB16-43", "03")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%rational, &
      & threebody_damping_function%rational)
   if (allocated(error)) return
   call test_damping_3b_numgrad(error, mol, d4, damp, param)

end subroutine test_grad_rational_3b_mb03

subroutine test_grad_rational_3b_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   ! PBE-D4-ATM(BJ)
   param = param_type(&
      & s9=1.0, a1 = 0.38574991_wp, a2 = 4.80688534_wp)

   call get_structure(mol, "MB16-43", "04")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%rational, &
      & threebody_damping_function%rational)
   if (allocated(error)) return
   call test_damping_3b_numgrad(error, mol, d4, damp, param)

end subroutine test_grad_rational_3b_mb04


subroutine test_damp_screened_3b_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   real(wp), parameter :: ref9(16, 16) = reshape([&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 1.07295754928311E-01_wp, 2.03093057580183E-01_wp, 9.93259734077058E-01_wp,&
      & 9.98000976303815E-01_wp, 4.72230123717093E-01_wp, 9.98428501706680E-01_wp,&
      & 2.35914372342270E-01_wp, 9.98436931536776E-01_wp, 9.98443649648898E-01_wp,&
      & 7.95626385241856E-01_wp, 1.15013142555723E-02_wp, 1.53617081588838E-01_wp,&
      & 9.85038800019169E-01_wp, 2.27847186738620E-01_wp, 0.00000000000000E+00_wp,&
      & 1.07295754928311E-01_wp, 0.00000000000000E+00_wp, 2.94997392958996E-02_wp,&
      & 8.79889450403619E-02_wp, 1.07406850802627E-01_wp, 7.14143053919836E-02_wp,&
      & 1.07454848212686E-01_wp, 1.97261318391755E-02_wp, 1.07458845957096E-01_wp,&
      & 1.07452912296645E-01_wp, 9.89214084862321E-02_wp, 3.37735968437893E-02_wp,&
      & 5.84743660205558E-04_wp, 1.05986994164686E-01_wp, 2.48291806338567E-02_wp,&
      & 0.00000000000000E+00_wp, 2.03093057580183E-01_wp, 2.94997392958996E-02_wp,&
      & 0.00000000000000E+00_wp, 2.73085778907253E-01_wp, 2.74389330185878E-01_wp,&
      & 1.33294980023688E-01_wp, 2.74511253154613E-01_wp, 6.48352308618775E-02_wp,&
      & 2.74521501436158E-01_wp, 2.74511039161304E-01_wp, 2.17499916387597E-01_wp,&
      & 3.36218570924160E-03_wp, 4.22351996474504E-02_wp, 2.70825477080916E-01_wp,&
      & 6.14183891411338E-02_wp, 0.00000000000000E+00_wp, 9.93259734077058E-01_wp,&
      & 8.79889450403619E-02_wp, 2.73085778907253E-01_wp, 0.00000000000000E+00_wp,&
      & 9.91348784124517E-01_wp, 6.61113206177342E-01_wp, 9.94732471458364E-01_wp,&
      & 1.86592035401233E-01_wp, 9.94769502586781E-01_wp, 8.73216932392309E-01_wp,&
      & 9.15736625773958E-01_wp, 3.12650337447421E-01_wp, 8.00078486341653E-03_wp,&
      & 9.32780801528847E-01_wp, 2.30307981511207E-01_wp, 0.00000000000000E+00_wp,&
      & 9.98000976303815E-01_wp, 1.07406850802627E-01_wp, 2.74389330185878E-01_wp,&
      & 9.91348784124517E-01_wp, 0.00000000000000E+00_wp, 6.64268820130590E-01_wp,&
      & 9.99455387597621E-01_wp, 1.57100112111396E-01_wp, 9.99517949005366E-01_wp,&
      & 3.91401758406345E-01_wp, 9.18334297759345E-01_wp, 3.14142577324156E-01_wp,&
      & 1.50952602019618E-01_wp, 3.03062215415006E-02_wp, 1.65066881092263E-01_wp,&
      & 0.00000000000000E+00_wp, 4.72230123717093E-01_wp, 7.14143053919836E-02_wp,&
      & 1.33294980023688E-01_wp, 6.61113206177342E-01_wp, 6.64268820130590E-01_wp,&
      & 0.00000000000000E+00_wp, 6.43057061985552E-01_wp, 1.52943913567527E-01_wp,&
      & 6.53671615244436E-01_wp, 6.64563516073725E-01_wp, 5.18952486183161E-01_wp,&
      & 9.08903510379644E-03_wp, 1.02216252484974E-01_wp, 6.55545751298823E-01_wp,&
      & 2.37475121939602E-02_wp, 0.00000000000000E+00_wp, 9.98428501706680E-01_wp,&
      & 1.07454848212686E-01_wp, 2.74511253154613E-01_wp, 9.94732471458364E-01_wp,&
      & 9.99455387597621E-01_wp, 6.43057061985552E-01_wp, 0.00000000000000E+00_wp,&
      & 2.34596119701894E-01_wp, 2.62448450742322E-02_wp, 9.98432171714139E-01_wp,&
      & 9.07490586506203E-01_wp, 3.12900933280019E-01_wp, 1.53840294951054E-01_wp,&
      & 9.67506163485578E-01_wp, 1.18449923773657E-02_wp, 0.00000000000000E+00_wp,&
      & 2.35914372342270E-01_wp, 1.97261318391755E-02_wp, 6.48352308618775E-02_wp,&
      & 1.86592035401233E-01_wp, 1.57100112111396E-01_wp, 1.52943913567527E-01_wp,&
      & 2.34596119701894E-01_wp, 0.00000000000000E+00_wp, 2.36267186406483E-01_wp,&
      & 1.55167223686867E-01_wp, 2.13792783647549E-01_wp, 7.35290753580661E-02_wp,&
      & 2.10170834642848E-03_wp, 2.59464659939155E-02_wp, 3.28166703538808E-03_wp,&
      & 0.00000000000000E+00_wp, 9.98436931536776E-01_wp, 1.07458845957096E-01_wp,&
      & 2.74521501436158E-01_wp, 9.94769502586781E-01_wp, 9.99517949005366E-01_wp,&
      & 6.53671615244436E-01_wp, 2.62448450742322E-02_wp, 2.36267186406483E-01_wp,&
      & 0.00000000000000E+00_wp, 9.99960763373456E-01_wp, 9.19889619078767E-01_wp,&
      & 3.13808195513043E-01_wp, 1.53850578179327E-01_wp, 9.86512565382484E-01_wp,&
      & 1.01442538124418E-01_wp, 0.00000000000000E+00_wp, 9.98443649648898E-01_wp,&
      & 1.07452912296645E-01_wp, 2.74511039161304E-01_wp, 8.73216932392309E-01_wp,&
      & 3.91401758406345E-01_wp, 6.64563516073725E-01_wp, 9.98432171714139E-01_wp,&
      & 1.55167223686867E-01_wp, 9.99960763373456E-01_wp, 0.00000000000000E+00_wp,&
      & 9.20510243653362E-01_wp, 3.14282085938683E-01_wp, 1.38317302129852E-01_wp,&
      & 3.05978390299168E-02_wp, 1.55334612233348E-01_wp, 0.00000000000000E+00_wp,&
      & 7.95626385241856E-01_wp, 9.89214084862321E-02_wp, 2.17499916387597E-01_wp,&
      & 9.15736625773958E-01_wp, 9.18334297759345E-01_wp, 5.18952486183161E-01_wp,&
      & 9.07490586506203E-01_wp, 2.13792783647549E-01_wp, 9.19889619078767E-01_wp,&
      & 9.20510243653362E-01_wp, 0.00000000000000E+00_wp, 3.74358659963616E-02_wp,&
      & 1.41625695206204E-01_wp, 8.95320272897329E-01_wp, 1.93515715827938E-02_wp,&
      & 0.00000000000000E+00_wp, 1.15013142555723E-02_wp, 3.37735968437893E-02_wp,&
      & 3.36218570924160E-03_wp, 3.12650337447421E-01_wp, 3.14142577324156E-01_wp,&
      & 9.08903510379644E-03_wp, 3.12900933280019E-01_wp, 7.35290753580661E-02_wp,&
      & 3.13808195513043E-01_wp, 3.14282085938683E-01_wp, 3.74358659963616E-02_wp,&
      & 0.00000000000000E+00_wp, 4.83514848134555E-02_wp, 3.10036191315462E-01_wp,&
      & 1.43663108670396E-02_wp, 0.00000000000000E+00_wp, 1.53617081588838E-01_wp,&
      & 5.84743660205558E-04_wp, 4.22351996474504E-02_wp, 8.00078486341653E-03_wp,&
      & 1.50952602019618E-01_wp, 1.02216252484974E-01_wp, 1.53840294951054E-01_wp,&
      & 2.10170834642848E-03_wp, 1.53850578179327E-01_wp, 1.38317302129852E-01_wp,&
      & 1.41625695206204E-01_wp, 4.83514848134555E-02_wp, 0.00000000000000E+00_wp,&
      & 1.11151516336064E-01_wp, 2.14914962671584E-02_wp, 0.00000000000000E+00_wp,&
      & 9.85038800019169E-01_wp, 1.05986994164686E-01_wp, 2.70825477080916E-01_wp,&
      & 9.32780801528847E-01_wp, 3.03062215415006E-02_wp, 6.55545751298823E-01_wp,&
      & 9.67506163485578E-01_wp, 2.59464659939155E-02_wp, 9.86512565382484E-01_wp,&
      & 3.05978390299168E-02_wp, 8.95320272897329E-01_wp, 3.10036191315462E-01_wp,&
      & 1.11151516336064E-01_wp, 0.00000000000000E+00_wp, 1.96837768768324E-02_wp,&
      & 0.00000000000000E+00_wp, 2.27847186738620E-01_wp, 2.48291806338567E-02_wp,&
      & 6.14183891411338E-02_wp, 2.30307981511207E-01_wp, 1.65066881092263E-01_wp,&
      & 2.37475121939602E-02_wp, 1.18449923773657E-02_wp, 3.28166703538808E-03_wp,&
      & 1.01442538124418E-01_wp, 1.55334612233347E-01_wp, 1.93515715827938E-02_wp,&
      & 1.43663108670396E-02_wp, 2.14914962671584E-02_wp, 1.96837768768323E-02_wp,&
      & 0.00000000000000E+00_wp], shape(ref9))

   ! PBE-D4-ATM(sc)
   param = param_type(&
      & s9 = 1.0_wp, a1 = 0.38574991_wp, a2 = 4.80688534_wp, a3 = 0.6_wp, a4 = 0.6_wp)

   call get_structure(mol, "MB16-43", "01")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%screened, &
      & threebody_damping_function%screened)
   if (allocated(error)) return
   call test_damping_3b_gen(error, mol, d4, damp, param, ref9)

end subroutine test_damp_screened_3b_mb01

subroutine test_damp_screened_3b_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   real(wp), parameter :: ref9(16, 16) = reshape([&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 7.30555435687545E-02_wp, 8.19679222733779E-01_wp, 2.17619523462881E-01_wp,&
      & 4.20961042640968E-01_wp, 5.35973473048937E-01_wp, 7.70630208827885E-01_wp,&
      & 4.13441423593889E-01_wp, 8.62159309477586E-01_wp, 3.05384077892867E-03_wp,&
      & 1.83946990088554E-01_wp, 8.52832594486557E-01_wp, 8.66531949087696E-01_wp,&
      & 4.93480856127545E-01_wp, 8.60836230705940E-01_wp, 0.00000000000000E+00_wp,&
      & 7.30555435687545E-02_wp, 0.00000000000000E+00_wp, 9.59505618982751E-01_wp,&
      & 1.41535735758051E-01_wp, 4.94265631529614E-01_wp, 7.33113881529456E-01_wp,&
      & 3.68187082989642E-02_wp, 6.91432812594211E-02_wp, 9.96590394338660E-01_wp,&
      & 2.69781771676147E-02_wp, 9.67449147879111E-01_wp, 9.99313815935570E-01_wp,&
      & 9.86111022568024E-01_wp, 6.10717695001026E-01_wp, 9.38702024361094E-01_wp,&
      & 0.00000000000000E+00_wp, 8.19679222733779E-01_wp, 9.59505618982751E-01_wp,&
      & 0.00000000000000E+00_wp, 9.08853579513940E-01_wp, 1.32534532134884E-02_wp,&
      & 7.11850756709475E-01_wp, 9.59674599466436E-01_wp, 9.49791581755763E-01_wp,&
      & 2.75760930435853E-02_wp, 3.70551389760579E-02_wp, 1.07162508759830E-01_wp,&
      & 9.51561640109324E-01_wp, 9.59674158113458E-01_wp, 6.26976034295376E-01_wp,&
      & 9.41438896457459E-01_wp, 0.00000000000000E+00_wp, 2.17619523462881E-01_wp,&
      & 1.41535735758051E-01_wp, 9.08853579513940E-01_wp, 0.00000000000000E+00_wp,&
      & 4.93235578475587E-01_wp, 7.41724377748335E-01_wp, 9.49862708359327E-01_wp,&
      & 1.70536900341949E-01_wp, 8.85922574796557E-01_wp, 3.78250599225937E-02_wp,&
      & 1.29658740625132E-01_wp, 9.22350213369940E-02_wp, 9.99799704811262E-01_wp,&
      & 6.49660574849514E-01_wp, 1.09223718094233E-01_wp, 0.00000000000000E+00_wp,&
      & 4.20961042640968E-01_wp, 4.94265631529614E-01_wp, 1.32534532134884E-02_wp,&
      & 4.93235578475587E-01_wp, 0.00000000000000E+00_wp, 3.66668609754361E-01_wp,&
      & 4.94351153366459E-01_wp, 4.89178427164992E-01_wp, 1.41031744259146E-01_wp,&
      & 1.60615484902270E-02_wp, 2.74550077327247E-01_wp, 4.94330438241347E-01_wp,&
      & 4.94350924706956E-01_wp, 3.15817721178343E-01_wp, 4.93365583137589E-01_wp,&
      & 0.00000000000000E+00_wp, 5.35973473048937E-01_wp, 7.33113881529456E-01_wp,&
      & 7.11850756709475E-01_wp, 7.41724377748335E-01_wp, 3.66668609754361E-01_wp,&
      & 0.00000000000000E+00_wp, 7.41756505946967E-01_wp, 6.97551854070343E-01_wp,&
      & 7.39368089129866E-01_wp, 1.09259744161709E-03_wp, 7.32146549363313E-01_wp,&
      & 7.41762945567340E-01_wp, 7.41658346618077E-01_wp, 3.41303373562005E-01_wp,&
      & 7.41745009167991E-01_wp, 0.00000000000000E+00_wp, 7.70630208827885E-01_wp,&
      & 3.68187082989642E-02_wp, 9.59674599466436E-01_wp, 9.49862708359327E-01_wp,&
      & 4.94351153366459E-01_wp, 7.41756505946967E-01_wp, 0.00000000000000E+00_wp,&
      & 9.09357071707728E-01_wp, 9.96771434264449E-01_wp, 3.94326579405436E-02_wp,&
      & 9.87038993057594E-01_wp, 9.99999902285071E-01_wp, 9.99981186795075E-01_wp,&
      & 6.54759785416111E-01_wp, 9.99969483371687E-01_wp, 0.00000000000000E+00_wp,&
      & 4.13441423593889E-01_wp, 6.91432812594211E-02_wp, 9.49791581755763E-01_wp,&
      & 1.70536900341949E-01_wp, 4.89178427164992E-01_wp, 6.97551854070343E-01_wp,&
      & 9.09357071707728E-01_wp, 0.00000000000000E+00_wp, 9.78920746240139E-01_wp,&
      & 1.13060237039229E-02_wp, 9.63308355700973E-01_wp, 9.91125036424010E-01_wp,&
      & 5.71677675296173E-02_wp, 5.75323505536952E-02_wp, 7.53079183274958E-02_wp,&
      & 0.00000000000000E+00_wp, 8.62159309477586E-01_wp, 9.96590394338660E-01_wp,&
      & 2.75760930435853E-02_wp, 8.85922574796557E-01_wp, 1.41031744259146E-01_wp,&
      & 7.39368089129866E-01_wp, 9.96771434264449E-01_wp, 9.78920746240139E-01_wp,&
      & 0.00000000000000E+00_wp, 3.91525341107812E-02_wp, 5.50795290759323E-01_wp,&
      & 9.92538160160780E-01_wp, 9.96770948092317E-01_wp, 6.49743681024023E-01_wp,&
      & 7.17751265569413E-01_wp, 0.00000000000000E+00_wp, 3.05384077892867E-03_wp,&
      & 2.69781771676147E-02_wp, 3.70551389760579E-02_wp, 3.78250599225937E-02_wp,&
      & 1.60615484902270E-02_wp, 1.09259744161709E-03_wp, 3.94326579405436E-02_wp,&
      & 1.13060237039229E-02_wp, 3.91525341107812E-02_wp, 0.00000000000000E+00_wp,&
      & 3.52553115891594E-02_wp, 3.94438977265268E-02_wp, 3.93456640505740E-02_wp,&
      & 1.52995547542866E-03_wp, 3.91771545521930E-02_wp, 0.00000000000000E+00_wp,&
      & 1.83946990088554E-01_wp, 9.67449147879111E-01_wp, 1.07162508759830E-01_wp,&
      & 1.29658740625132E-01_wp, 2.74550077327247E-01_wp, 7.32146549363313E-01_wp,&
      & 9.87038993057594E-01_wp, 9.63308355700973E-01_wp, 5.50795290759323E-01_wp,&
      & 3.52553115891594E-02_wp, 0.00000000000000E+00_wp, 9.00975627701351E-02_wp,&
      & 9.87039496654248E-01_wp, 6.45882548285014E-01_wp, 9.64329388355316E-01_wp,&
      & 0.00000000000000E+00_wp, 8.52832594486557E-01_wp, 9.99313815935570E-01_wp,&
      & 9.51561640109324E-01_wp, 9.22350213369940E-02_wp, 4.94330438241347E-01_wp,&
      & 7.41762945567340E-01_wp, 9.99999902285071E-01_wp, 9.91125036424010E-01_wp,&
      & 9.92538160160780E-01_wp, 3.94438977265268E-02_wp, 9.00975627701351E-02_wp,&
      & 0.00000000000000E+00_wp, 9.99999530777005E-01_wp, 6.54767067165967E-01_wp,&
      & 9.98037664068364E-01_wp, 0.00000000000000E+00_wp, 8.66531949087696E-01_wp,&
      & 9.86111022568024E-01_wp, 9.59674158113458E-01_wp, 9.99799704811262E-01_wp,&
      & 4.94350924706956E-01_wp, 7.41658346618077E-01_wp, 9.99981186795075E-01_wp,&
      & 5.71677675296173E-02_wp, 9.96770948092317E-01_wp, 3.93456640505740E-02_wp,&
      & 9.87039496654248E-01_wp, 9.99999530777005E-01_wp, 0.00000000000000E+00_wp,&
      & 5.95563966968360E-01_wp, 9.63566771372407E-01_wp, 0.00000000000000E+00_wp,&
      & 4.93480856127545E-01_wp, 6.10717695001026E-01_wp, 6.26976034295376E-01_wp,&
      & 6.49660574849514E-01_wp, 3.15817721178343E-01_wp, 3.41303373562005E-01_wp,&
      & 6.54759785416111E-01_wp, 5.75323505536952E-02_wp, 6.49743681024023E-01_wp,&
      & 1.52995547542866E-03_wp, 6.45882548285014E-01_wp, 6.54767067165967E-01_wp,&
      & 5.95563966968360E-01_wp, 0.00000000000000E+00_wp, 6.13557316378678E-01_wp,&
      & 0.00000000000000E+00_wp, 8.60836230705940E-01_wp, 9.38702024361094E-01_wp,&
      & 9.41438896457459E-01_wp, 1.09223718094233E-01_wp, 4.93365583137589E-01_wp,&
      & 7.41745009167991E-01_wp, 9.99969483371687E-01_wp, 7.53079183274958E-02_wp,&
      & 7.17751265569413E-01_wp, 3.91771545521930E-02_wp, 9.64329388355316E-01_wp,&
      & 9.98037664068364E-01_wp, 9.63566771372407E-01_wp, 6.13557316378678E-01_wp,&
      & 0.00000000000000E+00_wp], shape(ref9))

   ! PBE-D4-ATM(sc)
   param = param_type(&
      & s9 = 1.0_wp, a1 = 0.38574991_wp, a2 = 4.80688534_wp, a3 = 0.6_wp, a4 = 0.6_wp)

   call get_structure(mol, "MB16-43", "02")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%screened, &
      & threebody_damping_function%screened)
   if (allocated(error)) return
   call test_damping_3b_gen(error, mol, d4, damp, param, ref9)

end subroutine test_damp_screened_3b_mb02

subroutine test_grad_screened_3b_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   ! PBE-D4-ATM(sc)
   param = param_type(&
      & s9 = 1.0_wp, a1 = 0.38574991_wp, a2 = 4.80688534_wp, a3 = 0.6_wp, a4 = 0.6_wp)

   call get_structure(mol, "MB16-43", "03")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%screened, &
      & threebody_damping_function%screened)
   if (allocated(error)) return
   call test_damping_3b_numgrad(error, mol, d4, damp, param)

end subroutine test_grad_screened_3b_mb03

subroutine test_grad_screened_3b_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   ! PBE-D4-ATM(sc)
   param = param_type(&
      & s9 = 1.0_wp, a1 = 0.38574991_wp, a2 = 4.80688534_wp, a3 = 0.6_wp, a4 = 0.6_wp)

   call get_structure(mol, "MB16-43", "04")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%screened, &
      & threebody_damping_function%screened)
   if (allocated(error)) return
   call test_damping_3b_numgrad(error, mol, d4, damp, param)

end subroutine test_grad_screened_3b_mb04


subroutine test_damp_zero_3b_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   real(wp), parameter :: ref9(16, 16) = reshape([&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 3.84958274411613E-05_wp, 2.57090925826337E-04_wp, 2.76380681220778E-01_wp,&
      & 5.85628057053494E-01_wp, 3.28919025437052E-03_wp, 6.04524615775837E-01_wp,&
      & 2.70690284904394E-04_wp, 7.03259830076992E-01_wp, 6.67472917223392E-01_wp,&
      & 3.74012570939608E-03_wp, 2.11379465983135E-09_wp, 6.98207552116749E-05_wp,&
      & 1.30209247400614E-01_wp, 5.31828313468121E-05_wp, 0.00000000000000E+00_wp,&
      & 3.84958274411613E-05_wp, 0.00000000000000E+00_wp, 4.65599243789285E-08_wp,&
      & 5.59844972700630E-06_wp, 4.45087044251475E-05_wp, 6.96864423493808E-07_wp,&
      & 4.59875042901148E-05_wp, 2.82214255935219E-09_wp, 5.34847522317304E-05_wp,&
      & 5.07067683130700E-05_wp, 1.20483895701871E-06_wp, 2.78808803954404E-08_wp,&
      & 1.78107648037276E-14_wp, 9.78723947575291E-06_wp, 4.18376149231598E-09_wp,&
      & 0.00000000000000E+00_wp, 2.57090925826337E-04_wp, 4.65599243789285E-08_wp,&
      & 0.00000000000000E+00_wp, 3.34362589356344E-04_wp, 7.08482991190459E-04_wp,&
      & 4.41125287944397E-06_wp, 7.31808225252275E-04_wp, 3.25864836785096E-07_wp,&
      & 8.51086656710509E-04_wp, 8.07502952859222E-04_wp, 4.33733804150941E-06_wp,&
      & 3.16648861205810E-12_wp, 8.43245468252637E-08_wp, 1.57505736723684E-04_wp,&
      & 4.17664371582014E-08_wp, 0.00000000000000E+00_wp, 2.76380681220778E-01_wp,&
      & 5.59844972700630E-06_wp, 3.34362589356344E-04_wp, 0.00000000000000E+00_wp,&
      & 3.16108776632186E-01_wp, 5.00748215422631E-03_wp, 3.30158670816269E-01_wp,&
      & 2.84537292268091E-05_wp, 3.83993020029454E-01_wp, 2.26195773903379E-01_wp,&
      & 8.65653462403015E-03_wp, 2.01305034987913E-04_wp, 5.52255386584737E-10_wp,&
      & 4.64677167650136E-02_wp, 3.99462018873432E-05_wp, 0.00000000000000E+00_wp,&
      & 5.85628057053494E-01_wp, 4.45087044251475E-05_wp, 7.08482991190459E-04_wp,&
      & 3.16108776632186E-01_wp, 0.00000000000000E+00_wp, 1.06102704267787E-02_wp,&
      & 6.99047542339666E-01_wp, 3.89309403165595E-05_wp, 8.13646328904477E-01_wp,&
      & 2.89330991693689E-02_wp, 1.74107640152280E-02_wp, 4.26204897924299E-04_wp,&
      & 6.24304095724368E-05_wp, 7.88064372214096E-07_wp, 3.54463453250717E-06_wp,&
      & 0.00000000000000E+00_wp, 3.28919025437052E-03_wp, 6.96864423493808E-07_wp,&
      & 4.41125287944397E-06_wp, 5.00748215422631E-03_wp, 1.06102704267787E-02_wp,&
      & 0.00000000000000E+00_wp, 9.58574334405471E-03_wp, 4.32361896424984E-06_wp,&
      & 1.24431511829500E-02_wp, 1.20931704540142E-02_wp, 5.85150685510802E-05_wp,&
      & 6.92095307915803E-11_wp, 1.24559739114794E-06_wp, 2.35267208329314E-03_wp,&
      & 5.29839228969435E-10_wp, 0.00000000000000E+00_wp, 6.04524615775837E-01_wp,&
      & 4.59875042901148E-05_wp, 7.31808225252275E-04_wp, 3.30158670816269E-01_wp,&
      & 6.99047542339666E-01_wp, 9.58574334405471E-03_wp, 0.00000000000000E+00_wp,&
      & 2.91458631509662E-04_wp, 2.70900249799994E-06_wp, 7.89813913769793E-01_wp,&
      & 1.15702772892420E-02_wp, 3.56338361202660E-04_wp, 8.22091258099520E-05_wp,&
      & 1.22974430992778E-01_wp, 9.79843271185354E-10_wp, 0.00000000000000E+00_wp,&
      & 2.70690284904394E-04_wp, 2.82214255935219E-09_wp, 3.25864836785096E-07_wp,&
      & 2.84537292268091E-05_wp, 3.89309403165594E-05_wp, 4.32361896424983E-06_wp,&
      & 2.91458631509662E-04_wp, 0.00000000000000E+00_wp, 3.75723688947003E-04_wp,&
      & 4.18780573306241E-05_wp, 4.52806618089108E-06_wp, 1.30523183258156E-07_wp,&
      & 6.22706077053355E-13_wp, 1.35254435586739E-08_wp, 7.04982479747265E-13_wp,&
      & 0.00000000000000E+00_wp, 7.03259830076992E-01_wp, 5.34847522317304E-05_wp,&
      & 8.51086656710509E-04_wp, 3.83993020029454E-01_wp, 8.13646328904477E-01_wp,&
      & 1.24431511829500E-02_wp, 2.70900249799994E-06_wp, 3.75723688947003E-04_wp,&
      & 0.00000000000000E+00_wp, 9.27319703769622E-01_wp, 2.14286730145667E-02_wp,&
      & 4.91743825788054E-04_wp, 9.69787817653094E-05_wp, 1.80712871829786E-01_wp,&
      & 7.12161741604516E-07_wp, 0.00000000000000E+00_wp, 6.67472917223392E-01_wp,&
      & 5.07067683130700E-05_wp, 8.07502952859222E-04_wp, 2.26195773903379E-01_wp,&
      & 2.89330991693689E-02_wp, 1.20931704540142E-02_wp, 7.89813913769793E-01_wp,&
      & 4.18780573306242E-05_wp, 9.27319703769622E-01_wp, 0.00000000000000E+00_wp,&
      & 2.08541141801709E-02_wp, 4.86107639668330E-04_wp, 3.07908365977363E-05_wp,&
      & 9.28016472311441E-07_wp, 3.05920956752519E-06_wp, 0.00000000000000E+00_wp,&
      & 3.74012570939608E-03_wp, 1.20483895701871E-06_wp, 4.33733804150941E-06_wp,&
      & 8.65653462403015E-03_wp, 1.74107640152280E-02_wp, 5.85150685510802E-05_wp,&
      & 1.15702772892420E-02_wp, 4.52806618089108E-06_wp, 2.14286730145667E-02_wp,&
      & 2.08541141801709E-02_wp, 0.00000000000000E+00_wp, 1.34396310211122E-09_wp,&
      & 2.09065971920767E-06_wp, 2.33821552869626E-03_wp, 1.05920054446400E-10_wp,&
      & 0.00000000000000E+00_wp, 2.11379465983135E-09_wp, 2.78808803954404E-08_wp,&
      & 3.16648861205810E-12_wp, 2.01305034987913E-04_wp, 4.26204897924299E-04_wp,&
      & 6.92095307915803E-11_wp, 3.56338361202660E-04_wp, 1.30523183258156E-07_wp,&
      & 4.91743825788054E-04_wp, 4.86107639668330E-04_wp, 1.34396310211122E-09_wp,&
      & 0.00000000000000E+00_wp, 4.70522836888680E-08_wp, 9.21929768293452E-05_wp,&
      & 1.85122003794545E-11_wp, 0.00000000000000E+00_wp, 6.98207552116749E-05_wp,&
      & 1.78107648037276E-14_wp, 8.43245468252637E-08_wp, 5.52255386584737E-10_wp,&
      & 6.24304095724368E-05_wp, 1.24559739114794E-06_wp, 8.22091258099520E-05_wp,&
      & 6.22706077053355E-13_wp, 9.69787817653094E-05_wp, 3.07908365977363E-05_wp,&
      & 2.09065971920767E-06_wp, 4.70522836888680E-08_wp, 0.00000000000000E+00_wp,&
      & 7.29015621614462E-07_wp, 6.19344949195053E-11_wp, 0.00000000000000E+00_wp,&
      & 1.30209247400614E-01_wp, 9.78723947575291E-06_wp, 1.57505736723684E-04_wp,&
      & 4.64677167650136E-02_wp, 7.88064372214093E-07_wp, 2.35267208329314E-03_wp,&
      & 1.22974430992778E-01_wp, 1.35254435586739E-08_wp, 1.80712871829786E-01_wp,&
      & 9.28016472311439E-07_wp, 2.33821552869626E-03_wp, 9.21929768293452E-05_wp,&
      & 7.29015621614462E-07_wp, 0.00000000000000E+00_wp, 1.02515985188826E-09_wp,&
      & 0.00000000000000E+00_wp, 5.31828313468121E-05_wp, 4.18376149231597E-09_wp,&
      & 4.17664371582014E-08_wp, 3.99462018873431E-05_wp, 3.54463453250716E-06_wp,&
      & 5.29839228969434E-10_wp, 9.79843271185354E-10_wp, 7.04982479747264E-13_wp,&
      & 7.12161741604514E-07_wp, 3.05920956752518E-06_wp, 1.05920054446400E-10_wp,&
      & 1.85122003794545E-11_wp, 6.19344949195053E-11_wp, 1.02515985188826E-09_wp,&
      & 0.00000000000000E+00_wp], shape(ref9))

   ! PBE-D3-ATM(0)
   param = param_type(&
      & s9 = 1.0_wp, a1 = 1.0_wp, a2 = 0.0_wp, rs9 = 1.0_wp, alp = 14.0_wp)

   call get_structure(mol, "MB16-43", "01")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%zero, &
      & threebody_damping_function%zero)
   if (allocated(error)) return
   call test_damping_3b_gen(error, mol, d4, damp, param, ref9)

end subroutine test_damp_zero_3b_mb01

subroutine test_damp_zero_3b_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   real(wp), parameter :: ref9(16, 16) = reshape([&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 6.17905651622535E-06_wp, 8.88362089493197E-02_wp, 6.39813024378638E-05_wp,&
      & 1.20000665904744E-02_wp, 1.20905692719995E-02_wp, 4.68776333725664E-02_wp,&
      & 2.68550946770069E-04_wp, 1.80830425762267E-01_wp, 5.46766518703445E-11_wp,&
      & 3.00414410565887E-05_wp, 1.14892732809327E-01_wp, 1.93692487996373E-01_wp,&
      & 1.00703571000559E-02_wp, 6.31843732893026E-02_wp, 0.00000000000000E+00_wp,&
      & 6.17905651622534E-06_wp, 0.00000000000000E+00_wp, 8.27305421015712E-01_wp,&
      & 9.11281468565087E-05_wp, 8.54002270321799E-02_wp, 3.49541192037622E-01_wp,&
      & 5.76215508887912E-06_wp, 1.01111919667210E-05_wp, 9.79941014288468E-01_wp,&
      & 1.05540475092633E-07_wp, 7.99073781764198E-02_wp, 9.47302641914448E-01_wp,&
      & 8.07246311695763E-01_wp, 1.17361261678259E-01_wp, 9.57654645451599E-02_wp,&
      & 0.00000000000000E+00_wp, 8.88362089493197E-02_wp, 8.27305421015712E-01_wp,&
      & 0.00000000000000E+00_wp, 9.98098884396490E-02_wp, 2.95774980929268E-07_wp,&
      & 3.55807518265145E-01_wp, 8.44265418684717E-01_wp, 4.31777414478916E-01_wp,&
      & 3.80664818855444E-06_wp, 3.40474800728992E-06_wp, 4.27385760985077E-05_wp,&
      & 7.83764527086980E-01_wp, 8.44221028351088E-01_wp, 2.14665810683177E-01_wp,&
      & 4.02749695696973E-01_wp, 0.00000000000000E+00_wp, 6.39813024378638E-05_wp,&
      & 9.11281468565087E-05_wp, 9.98098884396490E-02_wp, 0.00000000000000E+00_wp,&
      & 6.65216121731541E-02_wp, 4.02985592909640E-01_wp, 2.32876024939611E-01_wp,&
      & 7.75303029562302E-05_wp, 1.07693132296473E-01_wp, 4.62367202084943E-07_wp,&
      & 3.42437331611572E-05_wp, 4.63041550171567E-05_wp, 9.13539427845179E-01_wp,&
      & 1.54248864119296E-01_wp, 4.24741258445952E-05_wp, 0.00000000000000E+00_wp,&
      & 1.20000665904744E-02_wp, 8.54002270321799E-02_wp, 2.95774980929268E-07_wp,&
      & 6.65216121731541E-02_wp, 0.00000000000000E+00_wp, 3.64764147951233E-02_wp,&
      & 8.65739377401451E-02_wp, 5.21357260582560E-02_wp, 1.13012608419811E-03_wp,&
      & 1.13314273509487E-07_wp, 6.44550906480921E-04_wp, 8.65032672005331E-02_wp,&
      & 8.65698107201649E-02_wp, 2.14999320925298E-02_wp, 7.99619135913060E-02_wp,&
      & 0.00000000000000E+00_wp, 1.20905692719995E-02_wp, 3.49541192037622E-01_wp,&
      & 3.55807518265145E-01_wp, 4.02985592909640E-01_wp, 3.64764147951233E-02_wp,&
      & 0.00000000000000E+00_wp, 4.21436896176385E-01_wp, 8.87937058546489E-02_wp,&
      & 4.19170628896386E-01_wp, 1.90864672541352E-11_wp, 2.48079746175605E-01_wp,&
      & 4.21492098404702E-01_wp, 4.21208201290204E-01_wp, 3.70303822448184E-02_wp,&
      & 4.18061328402200E-01_wp, 0.00000000000000E+00_wp, 4.68776333725664E-02_wp,&
      & 5.76215508887912E-06_wp, 8.44265418684717E-01_wp, 2.32876024939611E-01_wp,&
      & 8.65739377401451E-02_wp, 4.21436896176385E-01_wp, 0.00000000000000E+00_wp,&
      & 1.40192304007176E-01_wp, 9.94471514774982E-01_wp, 7.28969252644294E-06_wp,&
      & 5.90697340645504E-01_wp, 9.99900452420775E-01_wp, 9.99722320535322E-01_wp,&
      & 2.57485528838593E-01_wp, 9.89178452846266E-01_wp, 0.00000000000000E+00_wp,&
      & 2.68550946770069E-04_wp, 1.01111919667210E-05_wp, 4.31777414478916E-01_wp,&
      & 7.75303029562302E-05_wp, 5.21357260582560E-02_wp, 8.87937058546489E-02_wp,&
      & 1.40192304007176E-01_wp, 0.00000000000000E+00_wp, 4.64987660795181E-01_wp,&
      & 3.19062368690783E-09_wp, 3.15539248886049E-02_wp, 6.56610071886996E-01_wp,&
      & 1.15040390083971E-05_wp, 1.10830180587577E-05_wp, 1.19216030637952E-05_wp,&
      & 0.00000000000000E+00_wp, 1.80830425762267E-01_wp, 9.79941014288468E-01_wp,&
      & 3.80664818855444E-06_wp, 1.07693132296472E-01_wp, 1.13012608419811E-03_wp,&
      & 4.19170628896386E-01_wp, 9.94471514774982E-01_wp, 4.64987660795181E-01_wp,&
      & 0.00000000000000E+00_wp, 6.78235476171432E-06_wp, 7.26799711338704E-03_wp,&
      & 9.79808955535939E-01_wp, 9.94415867201521E-01_wp, 2.54268587942120E-01_wp,&
      & 6.14074937513014E-02_wp, 0.00000000000000E+00_wp, 5.46766518703444E-11_wp,&
      & 1.05540475092633E-07_wp, 3.40474800728992E-06_wp, 4.62367202084943E-07_wp,&
      & 1.13314273509487E-07_wp, 1.90864672541352E-11_wp, 7.28969252644294E-06_wp,&
      & 3.19062368690783E-09_wp, 6.78235476171432E-06_wp, 0.00000000000000E+00_wp,&
      & 1.63542639484609E-07_wp, 7.37470548062281E-06_wp, 6.99094343027975E-06_wp,&
      & 5.45204269563152E-11_wp, 2.93234752830162E-06_wp, 0.00000000000000E+00_wp,&
      & 3.00414410565887E-05_wp, 7.99073781764198E-02_wp, 4.27385760985077E-05_wp,&
      & 3.42437331611572E-05_wp, 6.44550906480921E-04_wp, 2.48079746175605E-01_wp,&
      & 5.90697340645504E-01_wp, 3.15539248886049E-02_wp, 7.26799711338704E-03_wp,&
      & 1.63542639484609E-07_wp, 0.00000000000000E+00_wp, 3.07008417219397E-05_wp,&
      & 5.93472036654064E-01_wp, 1.40678682186211E-01_wp, 5.59084480660837E-02_wp,&
      & 0.00000000000000E+00_wp, 1.14892732809327E-01_wp, 9.47302641914448E-01_wp,&
      & 7.83764527086980E-01_wp, 4.63041550171567E-05_wp, 8.65032672005331E-02_wp,&
      & 4.21492098404702E-01_wp, 9.99900452420775E-01_wp, 6.56610071886996E-01_wp,&
      & 9.79808955535939E-01_wp, 7.37470548062281E-06_wp, 3.07008417219397E-05_wp,&
      & 0.00000000000000E+00_wp, 9.99921550108017E-01_wp, 2.57514359624277E-01_wp,&
      & 8.65987944492338E-01_wp, 0.00000000000000E+00_wp, 1.93692487996373E-01_wp,&
      & 8.07246311695763E-01_wp, 8.44221028351088E-01_wp, 9.13539427845179E-01_wp,&
      & 8.65698107201649E-02_wp, 4.21208201290204E-01_wp, 9.99722320535322E-01_wp,&
      & 1.15040390083971E-05_wp, 9.94415867201521E-01_wp, 6.99094343027975E-06_wp,&
      & 5.93472036654064E-01_wp, 9.99921550108017E-01_wp, 0.00000000000000E+00_wp,&
      & 2.14773743193850E-01_wp, 5.40798875695556E-01_wp, 0.00000000000000E+00_wp,&
      & 1.00703571000559E-02_wp, 1.17361261678259E-01_wp, 2.14665810683177E-01_wp,&
      & 1.54248864119296E-01_wp, 2.14999320925298E-02_wp, 3.70303822448184E-02_wp,&
      & 2.57485528838593E-01_wp, 1.10830180587577E-05_wp, 2.54268587942120E-01_wp,&
      & 5.45204269563152E-11_wp, 1.40678682186211E-01_wp, 2.57514359624277E-01_wp,&
      & 2.14773743193850E-01_wp, 0.00000000000000E+00_wp, 1.00940093339268E-01_wp,&
      & 0.00000000000000E+00_wp, 6.31843732893026E-02_wp, 9.57654645451602E-02_wp,&
      & 4.02749695696973E-01_wp, 4.24741258445952E-05_wp, 7.99619135913060E-02_wp,&
      & 4.18061328402200E-01_wp, 9.89178452846266E-01_wp, 1.19216030637952E-05_wp,&
      & 6.14074937513014E-02_wp, 2.93234752830163E-06_wp, 5.59084480660837E-02_wp,&
      & 8.65987944492338E-01_wp, 5.40798875695556E-01_wp, 1.00940093339268E-01_wp,&
      & 0.00000000000000E+00_wp], shape(ref9))

   ! PBE-D3-ATM(0)
   param = param_type(&
      & s9 = 1.0_wp, a1 = 1.0_wp, a2 = 0.0_wp, rs9 = 1.0_wp, alp = 14.0_wp)

   call get_structure(mol, "MB16-43", "02")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%zero, &
      & threebody_damping_function%zero)
   if (allocated(error)) return
   call test_damping_3b_gen(error, mol, d4, damp, param, ref9)

end subroutine test_damp_zero_3b_mb02

subroutine test_grad_zero_3b_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   ! PBE-D3-ATM(0)
   param = param_type(&
      & s9 = 1.0_wp, a1 = 1.0_wp, a2 = 0.0_wp, rs9 = 1.0_wp, alp = 14.0_wp)

   call get_structure(mol, "MB16-43", "03")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%zero, &
      & threebody_damping_function%zero)
   if (allocated(error)) return
   call test_damping_3b_numgrad(error, mol, d4, damp, param)

end subroutine test_grad_zero_3b_mb03

subroutine test_grad_zero_3b_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   ! PBE-D3-ATM(0)
   param = param_type(&
      & s9 = 1.0_wp, a1 = 1.0_wp, a2 = 0.0_wp, rs9 = 1.0_wp, alp = 14.0_wp)

   call get_structure(mol, "MB16-43", "04")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%zero, &
      & threebody_damping_function%zero)
   if (allocated(error)) return
   call test_damping_3b_numgrad(error, mol, d4, damp, param)

end subroutine test_grad_zero_3b_mb04


subroutine test_damp_zero_avg_3b_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   real(wp), parameter :: ref9(16, 16) = reshape([&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 7.50926885874549E-01_wp, 1.03472039034664E-01_wp, 9.94955104253121E-01_wp,&
      & 9.95783862217283E-01_wp, 2.07685274526683E-01_wp, 9.63786991051598E-01_wp,&
      & 6.47430868672408E-01_wp, 9.88336498142147E-01_wp, 9.97907413576063E-01_wp,&
      & 2.05661058759505E-01_wp, 1.93901292997906E-03_wp, 6.17468626841848E-01_wp,&
      & 9.66936598439209E-01_wp, 6.43096184768207E-02_wp, 0.00000000000000E+00_wp,&
      & 7.50926885874549E-01_wp, 0.00000000000000E+00_wp, 5.31417187373872E-02_wp,&
      & 2.26692816667183E-02_wp, 4.78780732941612E-01_wp, 8.46813969073427E-02_wp,&
      & 5.65889202458667E-01_wp, 1.48238232765329E-03_wp, 8.99753781672112E-01_wp,&
      & 5.05804669593623E-01_wp, 9.81373822527633E-02_wp, 1.66062111671122E-02_wp,&
      & 2.61159224817707E-05_wp, 9.24227812067843E-02_wp, 1.95935785077942E-03_wp,&
      & 0.00000000000000E+00_wp, 1.03472039034664E-01_wp, 5.31417187373872E-02_wp,&
      & 0.00000000000000E+00_wp, 8.57812895031348E-01_wp, 8.85017538360905E-01_wp,&
      & 1.91338260784361E-02_wp, 7.77097294779149E-01_wp, 3.85040983965269E-02_wp,&
      & 9.45059304522389E-01_wp, 9.53211680378844E-01_wp, 1.75864828322975E-02_wp,&
      & 1.46891368546297E-04_wp, 3.50361495832334E-02_wp, 5.19843605244662E-01_wp,&
      & 3.85817118751364E-03_wp, 0.00000000000000E+00_wp, 9.94955104253121E-01_wp,&
      & 2.26692816667183E-02_wp, 8.57812895031348E-01_wp, 0.00000000000000E+00_wp,&
      & 8.65051285529819E-01_wp, 8.83490528779475E-01_wp, 9.65826955568633E-01_wp,&
      & 3.71918282220539E-02_wp, 9.96060591075288E-01_wp, 7.14768453737651E-01_wp,&
      & 8.03191374538036E-01_wp, 5.62671506986008E-01_wp, 9.65653698150175E-04_wp,&
      & 3.92179632480836E-01_wp, 5.15999712260301E-02_wp, 0.00000000000000E+00_wp,&
      & 9.95783862217283E-01_wp, 4.78780732941612E-01_wp, 8.85017538360905E-01_wp,&
      & 8.65051285529819E-01_wp, 0.00000000000000E+00_wp, 9.22404374343703E-01_wp,&
      & 9.68352675232523E-01_wp, 5.93565516508122E-02_wp, 9.97264863258457E-01_wp,&
      & 5.83201378499281E-01_wp, 5.56375745185779E-01_wp, 5.75937547611680E-01_wp,&
      & 1.03812253018588E-01_wp, 1.72930912671396E-02_wp, 2.66563492897457E-02_wp,&
      & 0.00000000000000E+00_wp, 2.07685274526683E-01_wp, 8.46813969073427E-02_wp,&
      & 1.91338260784361E-02_wp, 8.83490528779475E-01_wp, 9.22404374343703E-01_wp,&
      & 0.00000000000000E+00_wp, 4.45172724672839E-01_wp, 3.22352122208847E-02_wp,&
      & 7.82639417088819E-01_wp, 9.48756294945225E-01_wp, 4.06592767183819E-02_wp,&
      & 4.12253805677341E-04_wp, 4.16047631287577E-02_wp, 5.03952725435618E-01_wp,&
      & 8.12144253527105E-04_wp, 0.00000000000000E+00_wp, 9.63786991051598E-01_wp,&
      & 5.65889202458667E-01_wp, 7.77097294779149E-01_wp, 9.65826955568633E-01_wp,&
      & 9.68352675232523E-01_wp, 4.45172724672839E-01_wp, 0.00000000000000E+00_wp,&
      & 2.12561164612196E-01_wp, 1.04581750532789E-01_wp, 9.53884956521610E-01_wp,&
      & 3.70082242680622E-01_wp, 1.88079247180095E-01_wp, 2.51055559426600E-01_wp,&
      & 6.27729953884485E-01_wp, 1.86670540556215E-03_wp, 0.00000000000000E+00_wp,&
      & 6.47430868672407E-01_wp, 1.48238232765329E-03_wp, 3.85040983965268E-02_wp,&
      & 3.71918282220539E-02_wp, 5.93565516508121E-02_wp, 3.22352122208846E-02_wp,&
      & 2.12561164612196E-01_wp, 0.00000000000000E+00_wp, 7.55602867513857E-01_wp,&
      & 8.38166230635594E-02_wp, 2.10483952021063E-02_wp, 7.23028581722678E-03_wp,&
      & 8.54008747836572E-05_wp, 2.54264337107546E-03_wp, 8.90093861708039E-05_wp,&
      & 0.00000000000000E+00_wp, 9.88336498142147E-01_wp, 8.99753781672112E-01_wp,&
      & 9.45059304522389E-01_wp, 9.96060591075288E-01_wp, 9.97264863258457E-01_wp,&
      & 7.82639417088819E-01_wp, 1.04581750532789E-01_wp, 7.55602867513857E-01_wp,&
      & 0.00000000000000E+00_wp, 9.96669819983155E-01_wp, 8.04513882890817E-01_wp,&
      & 5.05720734448749E-01_wp, 7.44574806608069E-01_wp, 9.63024142162321E-01_wp,&
      & 3.82728725659678E-02_wp, 0.00000000000000E+00_wp, 9.97907413576063E-01_wp,&
      & 5.05804669593623E-01_wp, 9.53211680378844E-01_wp, 7.14768453737651E-01_wp,&
      & 5.83201378499281E-01_wp, 9.48756294945225E-01_wp, 9.53884956521610E-01_wp,&
      & 8.38166230635594E-02_wp, 9.96669819983155E-01_wp, 0.00000000000000E+00_wp,&
      & 8.34356377509535E-01_wp, 7.55571590504508E-01_wp, 8.31689337900791E-02_wp,&
      & 2.56863806968870E-02_wp, 3.56013742404578E-02_wp, 0.00000000000000E+00_wp,&
      & 2.05661058759505E-01_wp, 9.81373822527633E-02_wp, 1.75864828322975E-02_wp,&
      & 8.03191374538036E-01_wp, 5.56375745185779E-01_wp, 4.06592767183819E-02_wp,&
      & 3.70082242680622E-01_wp, 2.10483952021063E-02_wp, 8.04513882890817E-01_wp,&
      & 8.34356377509535E-01_wp, 0.00000000000000E+00_wp, 1.11088375699252E-03_wp,&
      & 3.51453696769664E-02_wp, 1.59659741941246E-01_wp, 4.76517017456913E-04_wp,&
      & 0.00000000000000E+00_wp, 1.93901292997906E-03_wp, 1.66062111671122E-02_wp,&
      & 1.46891368546297E-04_wp, 5.62671506986008E-01_wp, 5.75937547611680E-01_wp,&
      & 4.12253805677341E-04_wp, 1.88079247180095E-01_wp, 7.23028581722678E-03_wp,&
      & 5.05720734448749E-01_wp, 7.55571590504508E-01_wp, 1.11088375699252E-03_wp,&
      & 0.00000000000000E+00_wp, 8.49094044309582E-03_wp, 1.37202817117622E-01_wp,&
      & 2.64545940626275E-04_wp, 0.00000000000000E+00_wp, 6.17468626841848E-01_wp,&
      & 2.61159224817707E-05_wp, 3.50361495832334E-02_wp, 9.65653698150175E-04_wp,&
      & 1.03812253018588E-01_wp, 4.16047631287577E-02_wp, 2.51055559426600E-01_wp,&
      & 8.54008747836572E-05_wp, 7.44574806608069E-01_wp, 8.31689337900791E-02_wp,&
      & 3.51453696769664E-02_wp, 8.49094044309582E-03_wp, 0.00000000000000E+00_wp,&
      & 9.66715134613245E-03_wp, 3.96010896254589E-04_wp, 0.00000000000000E+00_wp,&
      & 9.66936598439209E-01_wp, 9.24227812067843E-02_wp, 5.19843605244662E-01_wp,&
      & 3.92179632480836E-01_wp, 1.72930912671395E-02_wp, 5.03952725435618E-01_wp,&
      & 6.27729953884485E-01_wp, 2.54264337107546E-03_wp, 9.63024142162321E-01_wp,&
      & 2.56863806968870E-02_wp, 1.59659741941246E-01_wp, 1.37202817117622E-01_wp,&
      & 9.66715134613245E-03_wp, 0.00000000000000E+00_wp, 1.07749713470057E-03_wp,&
      & 0.00000000000000E+00_wp, 6.43096184768206E-02_wp, 1.95935785077942E-03_wp,&
      & 3.85817118751364E-03_wp, 5.15999712260300E-02_wp, 2.66563492897457E-02_wp,&
      & 8.12144253527105E-04_wp, 1.86670540556215E-03_wp, 8.90093861708037E-05_wp,&
      & 3.82728725659678E-02_wp, 3.56013742404578E-02_wp, 4.76517017456913E-04_wp,&
      & 2.64545940626275E-04_wp, 3.96010896254589E-04_wp, 1.07749713470057E-03_wp,&
      & 0.00000000000000E+00_wp], shape(ref9))

   ! PBE-D3-ATM(0avg)
   param = param_type(&
      & s9 = 1.0_wp, a1 = 1.0_wp, a2 = 0.0_wp, rs9 = 1.0_wp, alp = 14.0_wp)

   call get_structure(mol, "MB16-43", "01")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%zero, &
      & threebody_damping_function%zero_avg)
   if (allocated(error)) return
   call test_damping_3b_gen(error, mol, d4, damp, param, ref9)

end subroutine test_damp_zero_avg_3b_mb01

subroutine test_damp_zero_avg_3b_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   real(wp), parameter :: ref9(16, 16) = reshape([&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,&
      & 8.11776046693404E-02_wp, 5.35904955325969E-01_wp, 1.17865230240150E-01_wp,&
      & 2.77285422826114E-01_wp, 2.37941080385233E-01_wp, 9.58790999777654E-01_wp,&
      & 9.27412053956574E-02_wp, 8.94294193948170E-01_wp, 4.07836244660862E-04_wp,&
      & 4.31754559159966E-02_wp, 9.57544988721795E-01_wp, 9.89431173054700E-01_wp,&
      & 2.16436700913049E-01_wp, 7.20054551091144E-01_wp, 0.00000000000000E+00_wp,&
      & 8.11776046693404E-02_wp, 0.00000000000000E+00_wp, 9.74007558567710E-01_wp,&
      & 3.85114870671393E-01_wp, 9.37710618848747E-01_wp, 8.74285904785547E-01_wp,&
      & 8.14738331909834E-01_wp, 1.24846658522385E-01_wp, 9.93988212645279E-01_wp,&
      & 2.08260643596266E-02_wp, 7.32258182090425E-01_wp, 9.97542504172901E-01_wp,&
      & 9.95124260814792E-01_wp, 7.48453950304181E-01_wp, 9.17809784473359E-01_wp,&
      & 0.00000000000000E+00_wp, 5.35904955325969E-01_wp, 9.74007558567710E-01_wp,&
      & 0.00000000000000E+00_wp, 7.36719880966844E-01_wp, 1.26028761812331E-02_wp,&
      & 9.66266709193354E-01_wp, 9.99811465789211E-01_wp, 7.65278647674620E-01_wp,&
      & 1.41048705649098E-01_wp, 3.50794269410874E-02_wp, 8.06695202659924E-02_wp,&
      & 9.92501555660204E-01_wp, 9.99462204218298E-01_wp, 8.40433660571830E-01_wp,&
      & 8.99912852414793E-01_wp, 0.00000000000000E+00_wp, 1.17865230240150E-01_wp,&
      & 3.85114870671393E-01_wp, 7.36719880966844E-01_wp, 0.00000000000000E+00_wp,&
      & 6.87527450092233E-01_wp, 9.25344328933978E-01_wp, 9.91363027990398E-01_wp,&
      & 1.63301706142477E-01_wp, 8.96937542133691E-01_wp, 2.39804037706329E-02_wp,&
      & 1.19911579765600E-01_wp, 7.81634156895843E-01_wp, 9.95488406154786E-01_wp,&
      & 7.17511936508872E-01_wp, 3.63452091923689E-01_wp, 0.00000000000000E+00_wp,&
      & 2.77285422826114E-01_wp, 9.37710618848747E-01_wp, 1.26028761812331E-02_wp,&
      & 6.87527450092232E-01_wp, 0.00000000000000E+00_wp, 8.45694375501546E-01_wp,&
      & 9.99553132042622E-01_wp, 5.32435882657587E-01_wp, 3.78921068288422E-01_wp,&
      & 5.29351380327714E-03_wp, 1.07712003721318E-01_wp, 9.93708150359169E-01_wp,&
      & 9.98381535631330E-01_wp, 4.90113700029164E-01_wp, 8.50152891812622E-01_wp,&
      & 0.00000000000000E+00_wp, 2.37941080385233E-01_wp, 8.74285904785547E-01_wp,&
      & 9.66266709193354E-01_wp, 9.25344328933978E-01_wp, 8.45694375501546E-01_wp,&
      & 0.00000000000000E+00_wp, 9.98903930235544E-01_wp, 4.71184326832522E-01_wp,&
      & 9.94660167865794E-01_wp, 3.20630583005129E-04_wp, 8.31637964870210E-01_wp,&
      & 9.99503758014888E-01_wp, 9.96532347923936E-01_wp, 3.36731542980010E-01_wp,&
      & 9.77231266863896E-01_wp, 0.00000000000000E+00_wp, 9.58790999777654E-01_wp,&
      & 8.14738331909834E-01_wp, 9.99811465789211E-01_wp, 9.91363027990398E-01_wp,&
      & 9.99553132042622E-01_wp, 9.98903930235544E-01_wp, 0.00000000000000E+00_wp,&
      & 9.78364375030667E-01_wp, 9.99960376579845E-01_wp, 8.11844641651628E-01_wp,&
      & 9.97225373124671E-01_wp, 9.99977331269260E-01_wp, 9.99960066989907E-01_wp,&
      & 9.98541519002600E-01_wp, 9.99448119762207E-01_wp, 0.00000000000000E+00_wp,&
      & 9.27412053956574E-02_wp, 1.24846658522385E-01_wp, 7.65278647674620E-01_wp,&
      & 1.63301706142477E-01_wp, 5.32435882657587E-01_wp, 4.71184326832522E-01_wp,&
      & 9.78364375030667E-01_wp, 0.00000000000000E+00_wp, 9.04258186307076E-01_wp,&
      & 2.16419240538481E-03_wp, 3.92299924138832E-01_wp, 9.91457018935562E-01_wp,&
      & 4.80505908738086E-01_wp, 3.49985670394344E-02_wp, 1.50529636675080E-01_wp,&
      & 0.00000000000000E+00_wp, 8.94294193948170E-01_wp, 9.93988212645279E-01_wp,&
      & 1.41048705649098E-01_wp, 8.96937542133691E-01_wp, 3.78921068288422E-01_wp,&
      & 9.94660167865794E-01_wp, 9.99960376579845E-01_wp, 9.04258186307076E-01_wp,&
      & 0.00000000000000E+00_wp, 2.00164794481779E-01_wp, 5.97740698297142E-01_wp,&
      & 9.98642300125744E-01_wp, 9.99834994409711E-01_wp, 9.53501647447812E-01_wp,&
      & 9.23136263541009E-01_wp, 0.00000000000000E+00_wp, 4.07836244660862E-04_wp,&
      & 2.08260643596266E-02_wp, 3.50794269410874E-02_wp, 2.39804037706329E-02_wp,&
      & 5.29351380327714E-03_wp, 3.20630583005129E-04_wp, 8.11844641651628E-01_wp,&
      & 2.16419240538481E-03_wp, 2.00164794481779E-01_wp, 0.00000000000000E+00_wp,&
      & 7.42229238910919E-03_wp, 8.05970016851937E-01_wp, 5.82740485522375E-01_wp,&
      & 4.18585371938164E-04_wp, 8.19968074227248E-02_wp, 0.00000000000000E+00_wp,&
      & 4.31754559159966E-02_wp, 7.32258182090425E-01_wp, 8.06695202659924E-02_wp,&
      & 1.19911579765600E-01_wp, 1.07712003721318E-01_wp, 8.31637964870210E-01_wp,&
      & 9.97225373124671E-01_wp, 3.92299924138832E-01_wp, 5.97740698297142E-01_wp,&
      & 7.42229238910919E-03_wp, 0.00000000000000E+00_wp, 5.75596408032770E-01_wp,&
      & 9.97611055693031E-01_wp, 6.42831619659488E-01_wp, 7.37585078013968E-01_wp,&
      & 0.00000000000000E+00_wp, 9.57544988721795E-01_wp, 9.97542504172901E-01_wp,&
      & 9.92501555660204E-01_wp, 7.81634156895843E-01_wp, 9.93708150359169E-01_wp,&
      & 9.99503758014888E-01_wp, 9.99977331269260E-01_wp, 9.91457018935562E-01_wp,&
      & 9.98642300125744E-01_wp, 8.05970016851937E-01_wp, 5.75596408032770E-01_wp,&
      & 0.00000000000000E+00_wp, 9.99986629023464E-01_wp, 9.98467750161668E-01_wp,&
      & 9.96879845936036E-01_wp, 0.00000000000000E+00_wp, 9.89431173054700E-01_wp,&
      & 9.95124260814792E-01_wp, 9.99462204218298E-01_wp, 9.95488406154786E-01_wp,&
      & 9.98381535631330E-01_wp, 9.96532347923936E-01_wp, 9.99960066989907E-01_wp,&
      & 4.80505908738086E-01_wp, 9.99834994409711E-01_wp, 5.82740485522375E-01_wp,&
      & 9.97611055693031E-01_wp, 9.99986629023464E-01_wp, 0.00000000000000E+00_wp,&
      & 9.70993162840575E-01_wp, 9.93594924411051E-01_wp, 0.00000000000000E+00_wp,&
      & 2.16436700913049E-01_wp, 7.48453950304181E-01_wp, 8.40433660571830E-01_wp,&
      & 7.17511936508872E-01_wp, 4.90113700029164E-01_wp, 3.36731542980010E-01_wp,&
      & 9.98541519002600E-01_wp, 3.49985670394344E-02_wp, 9.53501647447812E-01_wp,&
      & 4.18585371938164E-04_wp, 6.42831619659488E-01_wp, 9.98467750161668E-01_wp,&
      & 9.70993162840575E-01_wp, 0.00000000000000E+00_wp, 7.61967675436271E-01_wp,&
      & 0.00000000000000E+00_wp, 7.20054551091144E-01_wp, 9.17809784473359E-01_wp,&
      & 8.99912852414793E-01_wp, 3.63452091923689E-01_wp, 8.50152891812622E-01_wp,&
      & 9.77231266863896E-01_wp, 9.99448119762207E-01_wp, 1.50529636675080E-01_wp,&
      & 9.23136263541009E-01_wp, 8.19968074227249E-02_wp, 7.37585078013968E-01_wp,&
      & 9.96879845936036E-01_wp, 9.93594924411051E-01_wp, 7.61967675436271E-01_wp,&
      & 0.00000000000000E+00_wp], shape(ref9))

   ! PBE-D3-ATM(0avg)
   param = param_type(&
      & s9 = 1.0_wp, a1 = 1.0_wp, a2 = 0.0_wp, rs9 = 1.0_wp, alp = 14.0_wp)

   call get_structure(mol, "MB16-43", "02")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%zero, &
      & threebody_damping_function%zero_avg)
   if (allocated(error)) return
   call test_damping_3b_gen(error, mol, d4, damp, param, ref9)

end subroutine test_damp_zero_avg_3b_mb02

subroutine test_grad_zero_avg_3b_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   ! PBE-D3-ATM(0avg)
   param = param_type(&
      & s9 = 1.0_wp, a1 = 1.0_wp, a2 = 0.0_wp, rs9 = 1.0_wp, alp = 14.0_wp)

   call get_structure(mol, "MB16-43", "03")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%zero, &
      & threebody_damping_function%zero_avg)
   if (allocated(error)) return
   call test_damping_3b_numgrad(error, mol, d4, damp, param)

end subroutine test_grad_zero_avg_3b_mb03

subroutine test_grad_zero_avg_3b_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   ! PBE-D3-ATM(0avg)
   param = param_type(&
      & s9 = 1.0_wp, a1 = 1.0_wp, a2 = 0.0_wp, rs9 = 1.0_wp, alp = 14.0_wp)

   call get_structure(mol, "MB16-43", "04")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, twobody_damping_function%zero, &
      & threebody_damping_function%zero_avg)
   if (allocated(error)) return
   call test_damping_3b_numgrad(error, mol, d4, damp, param)

end subroutine test_grad_zero_avg_3b_mb04


subroutine test_damping_empty(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(damping_type) :: damp

   call new_damping(error, damp, -1, -1)
  
end subroutine test_damping_empty

subroutine test_damping_3b_empty(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(damping_type) :: damp

   call new_damping(error, damp, twobody_damping_function%rational, -1)
   if (allocated(error)) return
   if (allocated(damp%damping_3b)) then
      call test_failed(error, "3-body damping function should not be allocated")
   end if
  
end subroutine test_damping_3b_empty


subroutine test_params_rational(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(param_type) :: param
   type(damping_type) :: damp

   ! PBE-D4(BJ)
   param = param_type(&
      & s6 = 1.0_wp, s8 = 0.95948085_wp, a1 = 0.38574991_wp, a2 = 4.80688534_wp, &
      & s9=1.0)

   call new_damping(error, damp, twobody_damping_function%rational, &
      & threebody_damping_function%rational)
   if (allocated(error)) return
   call damp%check_params(error, param)
  
end subroutine test_params_rational

subroutine test_params_rational_empty(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(param_type) :: param
   type(damping_type) :: damp

   call new_damping(error, damp, twobody_damping_function%rational, &
      & threebody_damping_function%rational)
   if (allocated(error)) return
   call damp%check_params(error, param)
  
end subroutine test_params_rational_empty


subroutine test_params_screened(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(param_type) :: param
   type(damping_type) :: damp

   ! PBE-D4(sc)
   param = param_type(&
      & s6 = 1.0_wp, s8 = 0.95948085_wp, a1 = 0.38574991_wp, a2 = 4.80688534_wp, &
      & a3 = 0.6_wp, a4 = 0.6_wp, s9=1.0)

   call new_damping(error, damp, twobody_damping_function%screened, &
      & threebody_damping_function%screened)
   if (allocated(error)) return
   call damp%check_params(error, param)
  
end subroutine test_params_screened

subroutine test_params_screened_empty(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(param_type) :: param
   type(damping_type) :: damp

   call new_damping(error, damp, twobody_damping_function%screened, &
      & threebody_damping_function%screened)
   if (allocated(error)) return
   call damp%check_params(error, param)
  
end subroutine test_params_screened_empty


subroutine test_params_zero(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(param_type) :: param
   type(damping_type) :: damp

   ! PBE-D4(0)
   param = param_type(&
      & s6 = 1.0_wp, s8 = 0.722_wp, s9 = 1.0_wp, a1 = 1.0_wp, a2 = 0.0_wp, &
      & rs6 = 1.217_wp, rs8 = 1.0_wp, rs9 = 1.0_wp, alp = 14.0_wp)

   call new_damping(error, damp, twobody_damping_function%zero, &
      & threebody_damping_function%zero)
   if (allocated(error)) return
   call damp%check_params(error, param)
  
end subroutine test_params_zero

subroutine test_params_zero_empty(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(param_type) :: param
   type(damping_type) :: damp

   call new_damping(error, damp, twobody_damping_function%zero, &
      & threebody_damping_function%zero)
   if (allocated(error)) return
   call damp%check_params(error, param)
  
end subroutine test_params_zero_empty


subroutine test_params_mzero_zero_avg(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(param_type) :: param
   type(damping_type) :: damp

   ! B3LYP-D3(0M)
   param = param_type(&
      & s6 = 1.0_wp, s8 = 1.532981_wp, s9 = 1.0_wp, a1 = 1.0_wp, a2 = 0.0_wp, &
      & rs6 = 1.338153_wp, rs8 = 1.0_wp, rs9 = 1.0_wp, alp = 14.0_wp, bet=0.013988_wp)

   call new_damping(error, damp, twobody_damping_function%mzero, &
      & threebody_damping_function%zero_avg)
   if (allocated(error)) return
   call damp%check_params(error, param)
  
end subroutine test_params_mzero_zero_avg

subroutine test_params_mzero_zero_avg_empty(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(param_type) :: param
   type(damping_type) :: damp

   call new_damping(error, damp, twobody_damping_function%mzero, &
      & threebody_damping_function%zero_avg)
   if (allocated(error)) return
   call damp%check_params(error, param)
  
end subroutine test_params_mzero_zero_avg_empty


subroutine test_params_optpower(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(param_type) :: param
   type(damping_type) :: damp

   ! B3LYP-D3(op)
   param = param_type(&
      & s6 = 1.0_wp, s8 = 0.78311_wp, a1 = 0.30_wp, a2 = 4.25_wp, bet=4.0_wp)

   call new_damping(error, damp, twobody_damping_function%optpower, -1)
   if (allocated(error)) return
   call damp%check_params(error, param)
  
end subroutine test_params_optpower

subroutine test_params_optpower_empty(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(param_type) :: param
   type(damping_type) :: damp

   call new_damping(error, damp, twobody_damping_function%optpower, -1)
   if (allocated(error)) return
   call damp%check_params(error, param)
  
end subroutine test_params_optpower_empty


subroutine test_params_cso(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(param_type) :: param
   type(damping_type) :: damp

   ! B3LYP-D4(cso)
   param = param_type(&
      & s6 = 1.0_wp, a1 = 1.0_wp, a2 = 6.25_wp, a3 = 0.86_wp, a4 = 2.5_wp)

   call new_damping(error, damp, twobody_damping_function%cso, -1)
   if (allocated(error)) return
   call damp%check_params(error, param)
  
end subroutine test_params_cso

subroutine test_params_cso_empty(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(param_type) :: param
   type(damping_type) :: damp

   call new_damping(error, damp, twobody_damping_function%cso, -1)
   if (allocated(error)) return
   call damp%check_params(error, param)
  
end subroutine test_params_cso_empty


subroutine test_params_cso_noa1(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(param_type) :: param
   type(damping_type) :: damp

   ! B3LYP-D4(cso)
   param = param_type(&
      & s6 = 1.0_wp, a1 = 0.0_wp, a2 = 6.25_wp, a3 = 0.86_wp, a4 = 2.5_wp)

   call new_damping(error, damp, twobody_damping_function%cso, -1)
   if (allocated(error)) return
   call damp%check_params(error, param)
  
end subroutine test_params_cso_noa1


subroutine test_params_koide(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(param_type) :: param
   type(damping_type) :: damp

   ! PBE-D4(koide)
   param = param_type(&
      & s6 = 1.0_wp, s8 = 0.95948085_wp, a1 = 0.38574991_wp, a2 = 4.80688534_wp, &
      & rs6 = 2.0_wp, rs8 = 4.0_wp)

   call new_damping(error, damp, twobody_damping_function%koide, -1)
   if (allocated(error)) return
   call damp%check_params(error, param)
  
end subroutine test_params_koide

subroutine test_params_koide_empty(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(param_type) :: param
   type(damping_type) :: damp

   call new_damping(error, damp, twobody_damping_function%koide, -1)
   if (allocated(error)) return
   call damp%check_params(error, param)
  
end subroutine test_params_koide_empty


subroutine test_id_from_name(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: id_2b, id_3b

   call get_damping_function_id(error, damping_2b_name="rational", damping_2b_id=id_2b)
   if (allocated(error)) return
   call check(error, id_2b == twobody_damping_function%rational, &
      & "rational 2b name ID mismatch")
   if (allocated(error)) return

   call get_damping_function_id(error, damping_2b_name="screened", damping_2b_id=id_2b)
   if (allocated(error)) return
   call check(error, id_2b == twobody_damping_function%screened, &
      & "screened 2b name ID mismatch")
   if (allocated(error)) return

   call get_damping_function_id(error, damping_2b_name="zero", damping_2b_id=id_2b)
   if (allocated(error)) return
   call check(error, id_2b == twobody_damping_function%zero, &
      & "zero 2b name ID mismatch")
   if (allocated(error)) return

   call get_damping_function_id(error, damping_2b_name="mzero", damping_2b_id=id_2b)
   if (allocated(error)) return
   call check(error, id_2b == twobody_damping_function%mzero, &
      & "mzero 2b name ID mismatch")
   if (allocated(error)) return

   call get_damping_function_id(error, damping_2b_name="optpower", damping_2b_id=id_2b)
   if (allocated(error)) return
   call check(error, id_2b == twobody_damping_function%optpower, &
      & "optpower 2b name ID mismatch")
   if (allocated(error)) return

   call get_damping_function_id(error, damping_2b_name="cso", damping_2b_id=id_2b)
   if (allocated(error)) return
   call check(error, id_2b == twobody_damping_function%cso, &
      & "cso 2b name ID mismatch")
   if (allocated(error)) return

   call get_damping_function_id(error, damping_2b_name="koide", damping_2b_id=id_2b)
   if (allocated(error)) return
   call check(error, id_2b == twobody_damping_function%koide, &
      & "koide 2b name ID mismatch")
   if (allocated(error)) return

   ! Test case-insensitivity
   call get_damping_function_id(error, damping_2b_name="RATIONAL", damping_2b_id=id_2b)
   if (allocated(error)) return
   call check(error, id_2b == twobody_damping_function%rational, &
      & "RATIONAL 2b case-insensitivity mismatch")
   if (allocated(error)) return

   call get_damping_function_id(error, damping_3b_name="rational", damping_3b_id=id_3b)
   if (allocated(error)) return
   call check(error, id_3b == threebody_damping_function%rational, &
      & "rational 3b name ID mismatch")
   if (allocated(error)) return

   call get_damping_function_id(error, damping_3b_name="screened", damping_3b_id=id_3b)
   if (allocated(error)) return
   call check(error, id_3b == threebody_damping_function%screened, &
      & "screened 3b name ID mismatch")
   if (allocated(error)) return

   call get_damping_function_id(error, damping_3b_name="zero", damping_3b_id=id_3b)
   if (allocated(error)) return
   call check(error, id_3b == threebody_damping_function%zero, &
      & "zero 3b name ID mismatch")
   if (allocated(error)) return

   call get_damping_function_id(error, damping_3b_name="zero_avg", damping_3b_id=id_3b)
   if (allocated(error)) return
   call check(error, id_3b == threebody_damping_function%zero_avg, &
      & "zero_avg 3b name ID mismatch")

end subroutine test_id_from_name


subroutine test_id_from_name_empty(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: id_2b

   call get_damping_function_id(error, damping_2b_name="", &
      & damping_2b_id=id_2b)

end subroutine test_id_from_name_empty


subroutine test_id_from_name_wrong(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: id_2b

   call get_damping_function_id(error, damping_2b_name="unnamed_damping", &
      & damping_2b_id=id_2b)

end subroutine test_id_from_name_wrong


subroutine test_id_from_damping(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(damping_type) :: damp
   integer :: id_2b, id_3b

   call new_damping(error, damp, twobody_damping_function%rational, -1)
   if (allocated(error)) return
   call get_damping_function_id(error, damp, damping_2b_id=id_2b)
   if (allocated(error)) return
   call check(error, id_2b == twobody_damping_function%rational, &
      & "rational 2b damping ID mismatch")
   if (allocated(error)) return

   call new_damping(error, damp, twobody_damping_function%screened, -1)
   if (allocated(error)) return
   call get_damping_function_id(error, damp, damping_2b_id=id_2b)
   if (allocated(error)) return
   call check(error, id_2b == twobody_damping_function%screened, &
      & "screened 2b damping ID mismatch")
   if (allocated(error)) return

   call new_damping(error, damp, twobody_damping_function%zero, -1)
   if (allocated(error)) return
   call get_damping_function_id(error, damp, damping_2b_id=id_2b)
   if (allocated(error)) return
   call check(error, id_2b == twobody_damping_function%zero, &
      & "zero 2b damping ID mismatch")
   if (allocated(error)) return

   call new_damping(error, damp, twobody_damping_function%mzero, -1)
   if (allocated(error)) return
   call get_damping_function_id(error, damp, damping_2b_id=id_2b)
   if (allocated(error)) return
   call check(error, id_2b == twobody_damping_function%mzero, &
      & "mzero 2b damping ID mismatch")
   if (allocated(error)) return

   call new_damping(error, damp, twobody_damping_function%optpower, -1)
   if (allocated(error)) return
   call get_damping_function_id(error, damp, damping_2b_id=id_2b)
   if (allocated(error)) return
   call check(error, id_2b == twobody_damping_function%optpower, &
      & "optpower 2b damping ID mismatch")
   if (allocated(error)) return

   call new_damping(error, damp, twobody_damping_function%cso, -1)
   if (allocated(error)) return
   call get_damping_function_id(error, damp, damping_2b_id=id_2b)
   if (allocated(error)) return
   call check(error, id_2b == twobody_damping_function%cso, &
      & "cso 2b damping ID mismatch")
   if (allocated(error)) return

   call new_damping(error, damp, twobody_damping_function%koide, -1)
   if (allocated(error)) return
   call get_damping_function_id(error, damp, damping_2b_id=id_2b)
   if (allocated(error)) return
   call check(error, id_2b == twobody_damping_function%koide, &
      & "koide 2b damping ID mismatch")
   if (allocated(error)) return

   call new_damping(error, damp, twobody_damping_function%rational, &
      & threebody_damping_function%rational)
   if (allocated(error)) return
   call get_damping_function_id(error, damp, damping_3b_id=id_3b)
   if (allocated(error)) return
   call check(error, id_3b == threebody_damping_function%rational, &
      & "rational 3b damping ID mismatch")
   if (allocated(error)) return

   call new_damping(error, damp, twobody_damping_function%rational, &
      & threebody_damping_function%screened)
   if (allocated(error)) return
   call get_damping_function_id(error, damp, damping_3b_id=id_3b)
   if (allocated(error)) return
   call check(error, id_3b == threebody_damping_function%screened, &
      & "screened 3b damping ID mismatch")
   if (allocated(error)) return

   call new_damping(error, damp, twobody_damping_function%rational, &
      & threebody_damping_function%zero)
   if (allocated(error)) return
   call get_damping_function_id(error, damp, damping_3b_id=id_3b)
   if (allocated(error)) return
   call check(error, id_3b == threebody_damping_function%zero, &
      & "zero 3b damping ID mismatch")
   if (allocated(error)) return

   call new_damping(error, damp, twobody_damping_function%rational, &
      & threebody_damping_function%zero_avg)
   if (allocated(error)) return
   call get_damping_function_id(error, damp, damping_3b_id=id_3b)
   if (allocated(error)) return
   call check(error, id_3b == threebody_damping_function%zero_avg, &
      & "zero_avg 3b damping ID mismatch")

end subroutine test_id_from_damping


end module test_damping
