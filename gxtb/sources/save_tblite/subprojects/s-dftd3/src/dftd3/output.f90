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

module dftd3_output
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_convert, only : autoaa, autokcal, autoev
   use mctc_io_math, only : matinv_3x3
   use dftd3_damping, only : damping_param
   use dftd3_damping_mzero, only : mzero_damping_param
   use dftd3_damping_optimizedpower, only : optimizedpower_damping_param
   use dftd3_damping_rational, only : rational_damping_param
   use dftd3_damping_zero, only : zero_damping_param
   use dftd3_gcp, only : gcp_param
   use dftd3_model, only : d3_model
   use dftd3_version, only : get_dftd3_version
   implicit none
   private

   public :: ascii_atomic_radii, ascii_atomic_references, ascii_system_properties
   public :: ascii_energy_atom
   public :: ascii_results, ascii_damping_param, ascii_pairwise, ascii_gcp_param
   public :: turbomole_gradient, turbomole_gradlatt
   public :: json_results, tagged_result


contains


subroutine ascii_atomic_radii(unit, mol, disp)

   !> Unit for output
   integer, intent(in) :: unit

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Dispersion model
   class(d3_model), intent(in) :: disp

   integer :: isp

   write(unit, '(a,":")') "Atomic radii (in Angstrom)"
   write(unit, '(43("-"))')
   write(unit, '(a4,5x,*(1x,a10))') "Z", "R(cov)", "R(vdw)", "r4/r2"
   write(unit, '(43("-"))')
   do isp = 1, mol%nid
      write(unit, '(i4, 1x, a4, *(1x,f10.4))') &
         & mol%num(isp), mol%sym(isp), &
         & disp%rcov(isp)*autoaa, disp%rvdw(isp, isp)*autoaa/2, &
         & disp%r4r2(isp)*autoaa
   end do
   write(unit, '(43("-"))')

end subroutine ascii_atomic_radii


subroutine ascii_atomic_references(unit, mol, disp)

   !> Unit for output
   integer, intent(in) :: unit

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Dispersion model
   class(d3_model), intent(in) :: disp

   integer :: isp, iref, mref

   mref = maxval(disp%ref)
   write(unit, '(a,":")') "Atomic reference systems (in atomic units)"
   write(unit, '(76("-"))')
   write(unit, '(a4, 5x)', advance='no') "Z"
   do iref = 1, 3
      write(unit, '(a4, 1x, a7, 1x, a9)', advance='no') "#", "CN", "C6(AA)"
   end do
   write(unit, '(a)')
   write(unit, '(76("-"))')
   do isp = 1, mol%nid
      write(unit, '(i4, 1x, a4)', advance='no') &
         & mol%num(isp), mol%sym(isp)
      do iref = 1, disp%ref(isp)
         write(unit, '(i4, 1x, f7.4, 1x, f9.4)', advance='no') &
            iref, disp%cn(iref, isp), disp%c6(iref, iref, isp, isp)
         if ((iref == 3 .and. disp%ref(isp) > 3) .or. &
             (iref == 6 .and. disp%ref(isp) > 6)) then
            write(unit, '(/,9x)', advance='no')
         end if
      end do
      write(unit, '(a)')
   end do
   write(unit, '(76("-"))')

end subroutine ascii_atomic_references


subroutine ascii_system_properties(unit, mol, disp, cn, c6)

   !> Unit for output
   integer, intent(in) :: unit

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Dispersion model
   class(d3_model), intent(in) :: disp

   !> Coordination numbers
   real(wp), intent(in) :: cn(:)

   !> Atomic dispersion coefficients
   real(wp), intent(in) :: c6(:, :)

   integer :: iat, isp

   write(unit, '(a,":")') "Dispersion properties (in atomic units)"
   write(unit, '(56("-"))')
   write(unit, '(a6,1x,a4,5x,*(1x,a12))') "#", "Z", "CN", "C6(AA)", "C8(AA)"
   write(unit, '(56("-"))')
   do iat = 1, mol%nat
      isp = mol%id(iat)
      write(unit, '(i6,1x,i4,1x,a4,*(1x,f12.4))') &
         & iat, mol%num(isp), mol%sym(isp), cn(iat), c6(iat, iat), &
         & c6(iat, iat)*3*disp%r4r2(isp)**2
   end do
   write(unit, '(56("-"))')

end subroutine ascii_system_properties


!> Print atom-resolved dispersion energies
subroutine ascii_energy_atom(unit, mol, energies, label)

   !> Unit for output
   integer, intent(in) :: unit

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Atom-resolved dispersion energies
   real(wp), allocatable, intent(in) :: energies(:)

   !> Label for the output
   character(len=*), intent(in), optional :: label

   integer :: iat, isp
   character(len=:), allocatable :: label_

   label_ = "dispersion"
   if (present(label)) label_ = label

   write(unit, '(a,":")') "Atom-resolved "//label_//" energies"
   write(unit, '(50("-"))')
   write(unit, '(a6,1x,a4,1x,4x,a15,1x,a15)') "#", "Z", "[Hartree]", "[kcal/mol]"
   write(unit, '(50("-"))')
   do iat = 1, mol%nat
      isp = mol%id(iat)
      write(unit, '(i6,1x,i4,1x,a4,e15.8,1x,f15.8)') &
         & iat, mol%num(isp), mol%sym(isp), energies(iat), energies(iat)*autokcal
   end do
   write(unit, '(50("-"))')
   write(unit, '(a)')

end subroutine ascii_energy_atom

   
subroutine ascii_results(unit, mol, energy, gradient, sigma, label)

   !> Unit for output
   integer, intent(in) :: unit

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   real(wp), intent(in) :: energy
   real(wp), intent(in), optional :: gradient(:, :)
   real(wp), intent(in), optional :: sigma(:, :)

   !> Label for the output
   character(len=*), intent(in), optional :: label

   integer :: iat, isp
   logical :: grad
   character(len=1), parameter :: comp(3) = ["x", "y", "z"]
   character(len=:), allocatable :: label_

   label_ = "Dispersion"
   if (present(label)) label_ = label

   grad = present(gradient) .and. present(sigma)

   write(unit, '(a,":", t25, es20.13, 1x, a)') &
      & label_//" energy", energy, "Eh"
   write(unit, '(a)')
   if (grad) then
      write(unit, '(a,":", t25, es20.13, 1x, a)') &
         & "Gradient norm", norm2(gradient), "Eh/a0"
      write(unit, '(50("-"))')
      write(unit, '(a6,1x,a4,5x,*(1x,a10))') "#", "Z", "dE/dx", "dE/dy", "dE/dz"
      write(unit, '(50("-"))')
      do iat = 1, mol%nat
         isp = mol%id(iat)
         write(unit, '(i6,1x,i4,1x,a4,*(es11.3))') &
            & iat, mol%num(isp), mol%sym(isp), gradient(:, iat)
      end do
      write(unit, '(50("-"))')
      write(unit, '(a)')

      write(unit, '(a,":")') &
         & "Virial"
      write(unit, '(50("-"))')
      write(unit, '(a15,1x,*(1x,a10))') "component", "x", "y", "z"
      write(unit, '(50("-"))')
      do iat = 1, 3
         write(unit, '(2x,4x,1x,a4,1x,4x,*(es11.3))') &
            & comp(iat), sigma(:, iat)
      end do
      write(unit, '(50("-"))')
      write(unit, '(a)')
   end if

end subroutine ascii_results


subroutine ascii_pairwise(unit, mol, pair_disp2, pair_disp3)

   !> Unit for output
   integer, intent(in) :: unit

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   real(wp), intent(in) :: pair_disp2(:, :)
   real(wp), intent(in) :: pair_disp3(:, :)

   integer :: iat, jat, isp, jsp
   real(wp) :: disp, e2, e3

   e2 = 0.0_wp
   e3 = 0.0_wp

   write(unit, '(a,":")') "Pairwise representation of dispersion (in kcal/mol)"
   write(unit, '(82("-"))')
   write(unit, '(2(a6,1x,a4,5x),*(1x,a10:,1x,a7))') &
      "#", "Z", "#", "Z", "additive", "(rel.)", "non-add.", "(rel.)", "total"
   write(unit, '(82("-"))')
   do iat = 1, mol%nat
      isp = mol%id(iat)
      do jat = 1, mol%nat
         jsp = mol%id(jat)
         e2 = e2 + pair_disp2(jat, iat)
         e3 = e3 + pair_disp3(jat, iat)
         disp = pair_disp2(jat, iat) + pair_disp3(jat, iat)
         if (abs(disp) < epsilon(disp)) cycle
         write(unit, '(2(i6,1x,i4,1x,a4),*(1x,es10.2:,1x,"(",i4,"%)"))') &
            & iat, mol%num(isp), mol%sym(isp), &
            & jat, mol%num(jsp), mol%sym(jsp), &
            & pair_disp2(jat, iat) * autokcal, nint(pair_disp2(jat, iat)/disp*100), &
            & pair_disp3(jat, iat) * autokcal, nint(pair_disp3(jat, iat)/disp*100), &
            & disp * autokcal
      end do
   end do
   write(unit, '(82("-"))')
   disp = e2 + e3
   write(unit, '(1x, a, t33,*(1x,es10.2:,1x,"(",i4,"%)"))') &
      & "total dispersion energy", &
      & e2 * autokcal, nint(e2/disp*100), &
      & e3 * autokcal, nint(e3/disp*100), &
      & disp * autokcal
   write(unit, '(82("-"))')
   write(unit, '(a)')

end subroutine ascii_pairwise


subroutine ascii_damping_param(unit, param, method)

   !> Unit for output
   integer, intent(in) :: unit

   !> Damping parameters
   class(damping_param), intent(in) :: param

   !> Method name
   character(len=*), intent(in), optional :: method

   select type(param)
   type is (zero_damping_param)
      write(unit, '(a, ":", 1x)', advance="no") "Zero (Chai-Head-Gordon) damping"
      if (present(method)) then
         write(unit, '(a, "-")', advance="no") method
      end if
      if (abs(param%s9) > 0) then
         write(unit, '(a)') "D3(0)-ATM"
      else
         write(unit, '(a)') "D3(0)"
      end if
      write(unit, '(20("-"))')
      write(unit, '(a4, t10, f10.4)') &
         & "s6", param%s6, &
         & "s8", param%s8, &
         & "s9", param%s9, &
         & "rs6", param%rs6, &
         & "rs8", param%rs8, &
         & "alp", param%alp
      write(unit, '(20("-"))')
      write(unit, '(a)')
   type is (mzero_damping_param)
      write(unit, '(a, ":", 1x)', advance="no") "Modified zero damping"
      if (present(method)) then
         write(unit, '(a, "-")', advance="no") method
      end if
      if (abs(param%s9) > 0) then
         write(unit, '(a)') "D3(0M)-ATM"
      else
         write(unit, '(a)') "D3(0M)"
      end if
      write(unit, '(20("-"))')
      write(unit, '(a5, t10, f10.4)') &
         & "s6", param%s6, &
         & "s8", param%s8, &
         & "s9", param%s9, &
         & "rs6", param%rs6, &
         & "rs8", param%rs8, &
         & "alp", param%alp, &
         & "beta", param%bet
      write(unit, '(20("-"))')
      write(unit, '(a)')
   type is (optimizedpower_damping_param)
      write(unit, '(a, ":", 1x)', advance="no") "Optimized power damping"
      if (present(method)) then
         write(unit, '(a, "-")', advance="no") method
      end if
      if (abs(param%s9) > 0) then
         write(unit, '(a)') "D3(op)-ATM"
      else
         write(unit, '(a)') "D3(op)"
      end if
      write(unit, '(20("-"))')
      write(unit, '(a5, t10, f10.4)') &
         & "s6", param%s6, &
         & "s8", param%s8, &
         & "s9", param%s9, &
         & "a1", param%a1, &
         & "a1", param%a2, &
         & "alp", param%alp, &
         & "beta", param%bet
      write(unit, '(20("-"))')
      write(unit, '(a)')
   type is (rational_damping_param)
      write(unit, '(a, ":", 1x)', advance="no") "Rational (Becke-Johnson) damping"
      if (present(method)) then
         write(unit, '(a, "-")', advance="no") method
      end if
      if (abs(param%s9) > 0) then
         write(unit, '(a)') "D3(BJ)-ATM"
      else
         write(unit, '(a)') "D3(BJ)"
      end if
      write(unit, '(21("-"))')
      write(unit, '(a4, t10, f10.4)') &
         & "s6", param%s6, &
         & "s8", param%s8, &
         & "s9", param%s9, &
         & "a1", param%a1, &
         & "a2", param%a2, &
         & "alp", param%alp
      write(unit, '(20("-"))')
      write(unit, '(a)')
   end select

end subroutine ascii_damping_param


subroutine ascii_gcp_param(unit, mol, param, method)

   !> Unit for output
   integer, intent(in) :: unit

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Counter-poise parameters
   type(gcp_param), intent(in) :: param

   !> Method name
   character(len=*), intent(in), optional :: method

   integer :: isp

   write(unit, '(a,":")') "Global counter-poise parameters"
   write(unit, '(20("-"))')
   if (param%sigma > 0.0_wp .and. param%alpha > 0.0_wp .and. param%beta > 0.0_wp) then
      write(unit, '(a6, t10, f10.4)') &
         & "sigma", param%sigma, &
         & "alpha", param%alpha, &
         & "beta", param%beta
   end if
   if (param%damp) then
      write(unit, '(a6, t10, f10.4)') &
         & "dscal", param%dmp_scal, &
         & "dexpo", param%dmp_exp
   end if
   if (param%srb) then
      write(unit, '(a6, t10, f10.4)') &
         & "rscal", param%rscal, &
         & "qscal", param%qscal
   end if
   write(unit, '(20("-"))')
   write(unit, '(a)')

   if (allocated(param%emiss) .and. allocated(param%xv) &
      & .and. allocated(param%slater)) then
      write(unit, '(a,":")') "Atomic counter-poise parameters"
      write(unit, '(47("-"))')
      write(unit, '(a4,5x,a4,*(1x,a10))') "Z", "Zeff", "Emiss[Eh]", "Virtual", "Slater"
      write(unit, '(47("-"))')
      do isp = 1, mol%nid
         write(unit, '(i4, 1x, a4, i4, *(1x,f10.4))') &
            & mol%num(isp), mol%sym(isp), param%zeff(isp), &
            & param%emiss(isp), param%xv(isp), param%slater(isp)
      end do
      write(unit, '(47("-"))')
      write(unit, '(a)')
   end if

end subroutine ascii_gcp_param


subroutine turbomole_gradlatt(mol, fname, energy, sigma, stat)
   type(structure_type),intent(in) :: mol
   character(len=*),intent(in) :: fname
   real(wp),intent(in) :: energy
   real(wp),intent(in) :: sigma(3,3)
   integer, intent(out) :: stat
   character(len=:),allocatable :: line
   integer  :: i,j,icycle,line_number
   integer  :: err
   integer  :: igrad ! file handle
   logical  :: exist
   real(wp) :: escf
   real(wp) :: glat(3,3), inv_lat(3,3), gradlatt(3, 3)
   real(wp) :: dlat(3,3)
   stat = 0

   inv_lat = matinv_3x3(mol%lattice)

   do i = 1, 3
      do j = 1, 3
         gradlatt(i,j) = sigma(i,1)*inv_lat(j,1) &
            & + sigma(i,2)*inv_lat(j,2) &
            & + sigma(i,3)*inv_lat(j,3)
      enddo
   enddo

   icycle = 1
   i = 0
   escf = 0.0_wp
   line_number = 0

   inquire(file=fname,exist=exist)
   if (exist) then
      open(newunit=igrad,file=fname)
      read_file: do
         call getline(igrad,line,iostat=err)
         if (err.ne.0) exit read_file
         i=i+1
         if (index(line,'cycle') > 0) line_number = i
      enddo read_file
      if (line_number < 2) then
         stat = 1
         return
      endif

      rewind(igrad)
      skip_lines: do i = 1, line_number-1
         read(igrad,'(a)')
      enddo skip_lines
      call getline(igrad,line)
      read(line(10:17),*,iostat=err) icycle
      read(line(33:51),*,iostat=err) escf

      do i = 1, 3
         call getline(igrad,line)
         read(line,*,iostat=err) dlat(1,i),dlat(2,i),dlat(3,i)
      enddo
      if (any(abs(dlat-mol%lattice) > 1.0e-8_wp)) then
         stat = 1
         return
      endif
      do i = 1, 3
         call getline(igrad,line)
         read(line,*,iostat=err) glat(1,i),glat(2,i),glat(3,i)
      enddo
      do i = 1, 3
         backspace(igrad)
         backspace(igrad)
      enddo
      backspace(igrad)
   else
      open(newunit=igrad,file=fname)
      write(igrad,'("$gradlatt")')
   endif

   write(igrad,'(2x,"cycle =",1x,i6,4x,"SCF energy =",f18.11,3x,'//&
                   '"|dE/dlatt| =",f10.6)') &
      icycle, energy+escf, norm2(gradlatt+glat)
   do i = 1, 3
      write(igrad,'(3(F20.14,2x))') mol%lattice(1,i),mol%lattice(2,i),mol%lattice(3,i)
   enddo
   do i = 1, 3
      write(igrad,'(3D22.13)') gradlatt(1,i)+glat(1,i),gradlatt(2,i)+glat(2,i),gradlatt(3,i)+glat(3,i)
   enddo
   write(igrad,'("$end")')
   close(igrad)

end subroutine turbomole_gradlatt


subroutine turbomole_gradient(mol, fname, energy, gradient, stat)
   type(structure_type),intent(in) :: mol
   character(len=*),intent(in) :: fname
   real(wp),intent(in) :: energy
   real(wp),intent(in) :: gradient(:, :)
   integer, intent(out) :: stat
   character(len=:),allocatable :: line
   integer  :: i,icycle,line_number
   integer  :: err
   integer  :: igrad ! file handle
   logical  :: exist
   real(wp) :: escf
   real(wp),allocatable :: gscf(:,:)
   real(wp),allocatable :: xyz (:,:)
   allocate( gscf(3,mol%nat), source = 0.0_wp )
   stat = 0
   icycle = 1
   i = 0
   escf = 0.0_wp
   line_number = 0

   inquire(file=fname,exist=exist)
   if (exist) then
      open(newunit=igrad,file=fname)
      read_file: do
         call getline(igrad,line,iostat=err)
         if (err.ne.0) exit read_file
         i=i+1
         if (index(line,'cycle') > 0) line_number = i
      enddo read_file
      if (line_number < 2) then
         stat = 1
         return
      endif

      rewind(igrad)
      skip_lines: do i = 1, line_number-1
         read(igrad,'(a)')
      enddo skip_lines
      call getline(igrad,line)
      read(line(10:17),*,iostat=err) icycle
      read(line(33:51),*,iostat=err) escf

      allocate(xyz(3,mol%nat))
      do i = 1, mol%nat
         call getline(igrad,line)
         read(line,*,iostat=err) xyz(1,i),xyz(2,i),xyz(3,i)
      enddo
      if (any(abs(xyz-mol%xyz) > 1.0e-8_wp)) then
         stat = 1
         return
      endif
      do i = 1, mol%nat
         call getline(igrad,line)
         read(line,*,iostat=err) gscf(1,i),gscf(2,i),gscf(3,i)
      enddo
      do i = 1, mol%nat
         backspace(igrad)
         backspace(igrad)
      enddo
      backspace(igrad)
   else
      open(newunit=igrad,file=fname)
      write(igrad,'("$grad")')
   endif

   write(igrad,'(2x,"cycle =",1x,i6,4x,"SCF energy =",f18.11,3x,'//&
                   '"|dE/dxyz| =",f10.6)') &
      icycle, energy+escf, norm2(gradient+gscf)
   do i = 1, mol%nat
      write(igrad,'(3(F20.14,2x),4x,a2)') mol%xyz(1,i),mol%xyz(2,i),mol%xyz(3,i),mol%sym(i)
   enddo
   do i = 1, mol%nat
      write(igrad,'(3D22.13)') gradient(1,i)+gscf(1,i),gradient(2,i)+gscf(2,i),gradient(3,i)+gscf(3,i)
   enddo
   write(igrad,'("$end")')
   close(igrad)

end subroutine turbomole_gradient


!> reads a line from unit into an allocatable character
subroutine getline(unit,line,iostat)
   integer,intent(in) :: unit
   character(len=:),allocatable,intent(out) :: line
   integer,intent(out),optional :: iostat

   integer,parameter  :: buffersize=256
   character(len=buffersize) :: buffer
   integer :: size
   integer :: stat

   line = ''
   do
      read(unit,'(a)',advance='no',iostat=stat,size=size)  &
      &    buffer
      if (stat.gt.0) then
         if (present(iostat)) iostat=stat
         return ! an error occurred
      endif
      line = line // buffer(:size)
      if (stat.lt.0) then
         if (is_iostat_eor(stat)) stat = 0
         if (present(iostat)) iostat=stat
         return
      endif
   enddo

end subroutine getline


subroutine json_results(unit, indentation, energy, gradient, sigma, cn, c6, &
      & pairwise_energy2, pairwise_energy3, param)
   integer, intent(in) :: unit
   character(len=*), intent(in), optional :: indentation
   real(wp), intent(in), optional :: energy
   real(wp), intent(in), optional :: gradient(:, :)
   real(wp), intent(in), optional :: sigma(:, :)
   real(wp), intent(in), optional :: cn(:)
   real(wp), intent(in), optional :: c6(:, :)
   real(wp), intent(in), optional :: pairwise_energy2(:, :)
   real(wp), intent(in), optional :: pairwise_energy3(:, :)
   class(damping_param), intent(in), optional :: param
   character(len=:), allocatable :: indent, version_string
   character(len=*), parameter :: jsonkey = "('""',a,'"":',1x)"
   real(wp), allocatable :: array(:)

   call get_dftd3_version(string=version_string)

   if (present(indentation)) then
      indent = indentation
   end if

   write(unit, '("{")', advance='no')
   if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
   write(unit, jsonkey, advance='no') 'version'
   write(unit, '(1x,a)', advance='no') '"'//version_string//'"'
   if (present(energy)) then
      write(unit, '(",")', advance='no')
      if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
      write(unit, jsonkey, advance='no') 'energy'
      write(unit, '(1x,es25.16)', advance='no') energy
   end if
   if (present(sigma)) then
      write(unit, '(",")', advance='no')
      if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
      write(unit, jsonkey, advance='no') 'virial'
      array = reshape(sigma, [product(shape(sigma))])
      call write_json_array(unit, array, indent)
   end if
   if (present(gradient)) then
      write(unit, '(",")', advance='no')
      if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
      write(unit, jsonkey, advance='no') 'gradient'
      array = reshape(gradient, [product(shape(gradient))])
      call write_json_array(unit, array, indent)
   end if
   if (present(cn)) then
      write(unit, '(",")', advance='no')
      if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
      write(unit, jsonkey, advance='no') 'coordination numbers'
      call write_json_array(unit, cn, indent)
   end if
   if (present(c6)) then
      write(unit, '(",")', advance='no')
      if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
      write(unit, jsonkey, advance='no') 'c6 coefficients'
      array = reshape(c6, [product(shape(c6))])
      call write_json_array(unit, array, indent)
   end if
   if (present(pairwise_energy2)) then
      write(unit, '(",")', advance='no')
      if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
      write(unit, jsonkey, advance='no') 'additive pairwise energy'
      array = reshape(pairwise_energy2, [size(pairwise_energy2)])
      call write_json_array(unit, array, indent)
   end if
   if (present(pairwise_energy3)) then
      write(unit, '(",")', advance='no')
      if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
      write(unit, jsonkey, advance='no') 'non-additive pairwise energy'
      array = reshape(pairwise_energy3, [size(pairwise_energy3)])
      call write_json_array(unit, array, indent)
   end if
   if (present(param)) then
      write(unit, '(",")', advance='no')
      if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
      write(unit, jsonkey, advance='no') 'damping parameters'
      select type(param)
      type is(rational_damping_param)
         call write_json_param(unit, "rational", s6=param%s6, s8=param%s8, &
            & s9=param%s9, a1=param%a1, a2=param%a2, alp=param%alp, indent=indent)
      type is(zero_damping_param)
         call write_json_param(unit, "zero", s6=param%s6, s8=param%s8, &
            & s9=param%s9, rs6=param%rs6, rs8=param%rs8, alp=param%alp, indent=indent)
      type is(mzero_damping_param)
         call write_json_param(unit, "mzero", s6=param%s6, s8=param%s8, s9=param%s9, &
            & rs6=param%rs6, rs8=param%rs8, alp=param%alp, bet=param%bet, indent=indent)
      type is(optimizedpower_damping_param)
         call write_json_param(unit, "optimizedpower", s6=param%s6, s8=param%s8, &
            & s9=param%s9, a1=param%a1, a2=param%a2, alp=param%alp, bet=param%bet, &
            & indent=indent)
      class default
         call write_json_param(unit, "unknown", indent=indent)
      end select
   end if
   if (allocated(indent)) write(unit, '(/)', advance='no')
   write(unit, '("}")')

end subroutine json_results

subroutine write_json_param(unit, damping, s6, s8, s9, a1, a2, rs6, rs8, alp, bet, indent)
   integer, intent(in) :: unit
   character(len=*), intent(in) :: damping
   real(wp), intent(in), optional :: s6, s8, s9, a1, a2, rs6, rs8, alp, bet
   character(len=:), allocatable, intent(in) :: indent
   character(len=*), parameter :: jsonkeyval = "('""',a,'"":',1x,'""',a,'""')"

   write(unit, '("{")', advance='no')

   if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 2)
   write(unit, jsonkeyval, advance='no') 'damping', damping

   if (present(s6)) then
      call write_json_keyval(unit, 's6', s6, indent)
   end if

   if (present(s8)) then
      call write_json_keyval(unit, 's8', s8, indent)
   end if

   if (present(s9)) then
      call write_json_keyval(unit, 's9', s9, indent)
   end if

   if (present(a1)) then
      call write_json_keyval(unit, 'a1', a1, indent)
   end if

   if (present(a2)) then
      call write_json_keyval(unit, 'a2', a2, indent)
   end if

   if (present(rs6)) then
      call write_json_keyval(unit, 'rs6', rs6, indent)
   end if

   if (present(rs8)) then
      call write_json_keyval(unit, 'rs8', rs8, indent)
   end if

   if (present(alp)) then
      call write_json_keyval(unit, 'alp', alp, indent)
   end if

   if (present(bet)) then
      call write_json_keyval(unit, 'bet', bet, indent)
   end if

   if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
   write(unit, '("}")', advance='no')
end subroutine write_json_param

subroutine write_json_keyval(unit, key, val, indent)
   integer, intent(in) :: unit
   character(len=*), intent(in) :: key
   real(wp), intent(in) :: val
   character(len=:), allocatable, intent(in) :: indent
   character(len=*), parameter :: jsonkeyval = "('""',a,'"":',1x,es23.16)"

   write(unit, '(",")', advance='no')
   if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 2)
   write(unit, jsonkeyval, advance='no') key, val

end subroutine write_json_keyval

subroutine write_json_array(unit, array, indent)
   integer, intent(in) :: unit
   real(wp), intent(in) :: array(:)
   character(len=:), allocatable, intent(in) :: indent
   integer :: i
   write(unit, '("[")', advance='no')
   do i = 1, size(array)
      if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 2)
      write(unit, '(es23.16)', advance='no') array(i)
      if (i /= size(array)) write(unit, '(",")', advance='no')
   end do
   if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
   write(unit, '("]")', advance='no')
end subroutine write_json_array


subroutine tagged_result(unit, energy, gradient, sigma)
   integer, intent(in) :: unit
   real(wp), intent(in), optional :: energy
   real(wp), intent(in), optional :: gradient(:, :)
   real(wp), intent(in), optional :: sigma(:, :)
   character(len=*), parameter :: tag_header = &
      & '(a,t20,":",a,":",i0,":",*(i0:,","))'

   if (present(energy)) then
      write(unit, tag_header) "energy", "real", 0
      write(unit, '(3es24.16)') energy
   end if
   if (present(gradient)) then
      write(unit, tag_header) "gradient", "real", 2, 3, size(gradient, 2)
      write(unit, '(3es24.16)') gradient
   end if
   if (present(sigma)) then
      write(unit, tag_header) "virial", "real", 2, 3, 3
      write(unit, '(3es24.16)') sigma
   end if

end subroutine tagged_result

end module dftd3_output
