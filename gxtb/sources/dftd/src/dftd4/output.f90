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

module dftd4_output
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_convert, only : autoaa, autokcal, autoev
   use mctc_io_math, only : matinv_3x3
   use dftd4_damping, only : damping_type
   use dftd4_model, only : dispersion_model
   use dftd4_param, only : param_type
   use dftd4_version, only : get_dftd4_version
   implicit none
   private

   public :: ascii_atomic_radii, ascii_atomic_references, ascii_system_properties
   public :: ascii_results, ascii_damping_param, ascii_pairwise
   public :: turbomole_gradient, turbomole_gradlatt
   public :: json_results, tagged_result


contains


subroutine ascii_atomic_radii(unit, mol, disp)
   !DEC$ ATTRIBUTES DLLEXPORT :: ascii_atomic_radii

   !> Unit for output
   integer, intent(in) :: unit

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Dispersion model
   class(dispersion_model), intent(in) :: disp

   integer :: isp

   write(unit, '(a,":")') "Atomic data, radii in Ångström"
   write(unit, '(54("-"))')
   write(unit, '(a4,5x,*(1x,a10))') &
      "Z", "R(cov)", "r4/r2", "hardness", "EN"
   write(unit, '(54("-"))')
   do isp = 1, mol%nid
      write(unit, '(i4, 1x, a4, *(1x,f10.4))') &
         & mol%num(isp), mol%sym(isp), &
         & disp%rcov(isp)*autoaa, &
         & disp%r4r2(isp)*autoaa, &
         & disp%eta(isp), &
         & disp%en(isp)
   end do
   write(unit, '(54("-"))')
   write(unit, '(a)')

end subroutine ascii_atomic_radii


subroutine ascii_atomic_references(unit, mol, disp)
   !DEC$ ATTRIBUTES DLLEXPORT :: ascii_atomic_references

   !> Unit for output
   integer, intent(in) :: unit

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Dispersion model
   class(dispersion_model), intent(in) :: disp

   integer :: isp, iref, mref

   mref = maxval(disp%ref)
   write(unit, '(a,":")') "Atomic reference systems (in atomic units)"
   write(unit, '(76("-"))')
   write(unit, '(a4, 5x)', advance='no') "Z"
   do iref = 1, 2
      write(unit, '(a4, 2(1x, a7), 1x, a12)', advance='no') &
         "#", "CN", "q+Z", "C6(AA)"
   end do
   write(unit, '(a)')
   write(unit, '(76("-"))')
   do isp = 1, mol%nid
      write(unit, '(i4, 1x, a4)', advance='no') &
         & mol%num(isp), mol%sym(isp)
      do iref = 1, disp%ref(isp)
         write(unit, '(i4, 2(1x, f7.4), 1x, f12.4)', advance='no') &
            iref, disp%cn(iref, isp), disp%q(iref, isp) + disp%zeff(isp), &
            disp%c6(iref, iref, isp, isp)
         if (iref == 2 .and. disp%ref(isp) > 2) then
            write(unit, '(/,9x)', advance='no')
         end if
         if (iref == 4 .and. disp%ref(isp) > 4) then
            write(unit, '(/,9x)', advance='no')
         end if
         if (iref == 6 .and. disp%ref(isp) > 6) then
            write(unit, '(/,9x)', advance='no')
         end if
      end do
      write(unit, '(a)')
   end do
   write(unit, '(76("-"))')
   write(unit, '(a)')

end subroutine ascii_atomic_references


subroutine ascii_system_properties(unit, mol, disp, cn, q, c6, alpha, alphaqq)
   !DEC$ ATTRIBUTES DLLEXPORT :: ascii_system_properties

   !> Unit for output
   integer, intent(in) :: unit

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Dispersion model
   class(dispersion_model), intent(in) :: disp

   !> Coordination numbers
   real(wp), intent(in) :: cn(:)

   !> Atomic partial charges
   real(wp), intent(in) :: q(:)

   !> Atomic dispersion coefficients
   real(wp), intent(in) :: c6(:, :)

   !> Atomic static dipole-dipole polarizabilities
   real(wp), intent(in) :: alpha(:)

   !> Atomic static quadrupole-quadrupole polarizabilities
   real(wp), intent(in) :: alphaqq(:)

   integer :: iat, isp, jat
   real(wp) :: sum_c8

   sum_c8 = 0.0_wp

   write(unit, '(a,":")') "Atomic properties (in atomic units)"
   write(unit, '(87("-"))')
   write(unit, '(a6,1x,a4,5x,*(1x,a10))') "#", "Z", "CN", "q", "C6(AA)", &
      & "C8(AA)", "alpha(0)", "alphaQQ(0)"
   write(unit, '(87("-"))')
   do iat = 1, mol%nat
      isp = mol%id(iat)
      write(unit, '(i6,1x,i4,1x,a4, 2(1x,f10.4), 1x,f11.4, 1x,f13.4, 1x,f10.4)') &
         & iat, mol%num(isp), mol%sym(isp), cn(iat), q(iat), c6(iat, iat), &
         & c6(iat, iat)*3*disp%r4r2(isp)**2, alpha(iat), alphaqq(iat)
      do jat = 1, mol%nat
         sum_c8 = sum_c8 + 3*c6(jat, iat)*disp%r4r2(mol%id(jat))*disp%r4r2(isp)
      end do
   end do
   write(unit, '(87("-"))')
   write(unit, '(a)')

   write(unit, '(a,":")') "Molecular properties (in atomic units)"
   write(unit, '(40("-"))')
   write(unit, '(1x, a, t20, f19.4)') &
      "molecular C6",  sum(c6), &
      "molecular C8",  sum_c8
   write(unit, '(40("-"))')
   write(unit, '(a)')

end subroutine ascii_system_properties


subroutine ascii_results(unit, mol, energy, gradient, sigma)
   !DEC$ ATTRIBUTES DLLEXPORT :: ascii_results

   !> Unit for output
   integer, intent(in) :: unit

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   real(wp), intent(in) :: energy
   real(wp), intent(in), optional :: gradient(:, :)
   real(wp), intent(in), optional :: sigma(:, :)

   integer :: iat, isp
   logical :: grad
   character(len=1), parameter :: comp(3) = ["x", "y", "z"]

   grad = present(gradient) .and. present(sigma)

   write(unit, '(a,":", t25, es20.13, 1x, a)') &
      & "Dispersion energy", energy, "Eh"
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
   !DEC$ ATTRIBUTES DLLEXPORT :: ascii_pairwise

   !> Unit for output
   integer, intent(in) :: unit

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

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


subroutine ascii_damping_param(unit, verbosity, model, damp, param, method, modified)
   !DEC$ ATTRIBUTES DLLEXPORT :: ascii_damping_param
   
   !> Unit for output
   integer, intent(in) :: unit

   !> Level of verbosity
   integer, intent(in) :: verbosity
   
   !> Dispersion model
   class(dispersion_model), intent(in) :: model

   !> Damping function
   type(damping_type), intent(in) :: damp

   !> Damping parameters
   type(param_type), intent(in) :: param

   !> Method name
   character(len=*), intent(in), optional :: method

   !> Flag indicating whether the parameters were modified via the CLI
   logical, intent(in), optional :: modified

   logical :: atm

   ! Check for the ATM term
   atm = .false.
   if (allocated(param%s9)) then
      if (abs(param%s9) > 0) then
         atm = .true.
      end if
   end if

   ! Output the full method string including dispersion
   if (present(method)) then
      write(unit, '(a, ":", 1x)', advance="no") "Method"
      write(unit, '(a, "-")', advance="no") method
   end if

   write(unit, '(a, "(")', advance="no") model%label
   write(unit, '(a, ")")', advance="no") damp%damping_2b%label_short

   if (atm) then
      write(unit, '("-", a, "(")', advance="no") "ATM"
      write(unit, '(a, ")")') damp%damping_3b%label_short
   else
      write(unit, '(a)') ""
   end if

   ! Specify two and three-body damping functions
   if (verbosity > 2) then
      if (allocated(damp%damping_2b)) then
         write(unit, '(a, ":", 1x)', advance="no") "Two-body damping: "
         write(unit, '(a, 1x)', advance="yes") damp%damping_2b%label
      end if

      if (allocated(damp%damping_3b)) then
         write(unit, '(a, ":", 1x)', advance="no") "Three-body damping: "
         write(unit, '(a, 1x)', advance="yes") damp%damping_3b%label
      end if
   end if

   ! Output the damping parameters
   write(unit, '(21("-"))')

   if (allocated(param%s6)) then
      write(unit, '(a4, t10, f10.4)') "s6", param%s6
   end if
   
   if (allocated(param%s8)) then
      write(unit, '(a4, t10, f10.4)') "s8", param%s8
   end if
   
   if (allocated(param%s9)) then
      write(unit, '(a4, t10, f10.4)') "s9", param%s9
   end if

   write(unit, '(a4, t10, f10.4)') "a1", param%a1
   write(unit, '(a4, t10, f10.4)') "a2", param%a2

   if (allocated(param%a3)) then
      write(unit, '(a4, t10, f10.4)') "a3", param%a3
   end if

   if (allocated(param%a4)) then
      write(unit, '(a4, t10, f10.4)') "a4", param%a4
   end if

   if (allocated(param%rs6)) then
      write(unit, '(a4, t10, f10.4)') "rs6", param%rs6
   end if

   if (allocated(param%rs8)) then
      write(unit, '(a4, t10, f10.4)') "rs8", param%rs8
   end if

   if (allocated(param%rs9)) then
      write(unit, '(a4, t10, f10.4)') "rs9", param%rs9
   end if

   if (allocated(param%alp)) then
      write(unit, '(a4, t10, f10.4)') "alp", param%alp
   end if

   if (allocated(param%bet)) then
      write(unit, '(a4, t10, f10.4)') "bet", param%bet
   end if

   write(unit, '(20("-"))')

   ! Add warning if the default parameters were modified via the CLI
   if (verbosity > 1) then
      if (present(modified) .and. present(method)) then
         if (modified) then
            write(unit, '(a)') "Warning: The default "//method//" damping parameters were modified."
         end if
      end if
   end if
   write(unit, '(a)')


end subroutine ascii_damping_param


subroutine turbomole_gradlatt(mol, fname, energy, sigma, stat)
   !DEC$ ATTRIBUTES DLLEXPORT :: turbomole_gradlatt
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
   !DEC$ ATTRIBUTES DLLEXPORT :: turbomole_gradient
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


subroutine json_results(unit, indentation, energy, gradient, sigma, hessian, &
      & cn, q, c6, alpha, pairwise_energy2, pairwise_energy3)
   !DEC$ ATTRIBUTES DLLEXPORT :: json_results
   integer, intent(in) :: unit
   character(len=*), intent(in), optional :: indentation
   real(wp), intent(in), optional :: energy
   real(wp), intent(in), optional :: gradient(:, :)
   real(wp), intent(in), optional :: sigma(:, :)
   real(wp), intent(in), optional :: hessian(:, :, :, :)
   real(wp), intent(in), optional :: cn(:)
   real(wp), intent(in), optional :: q(:)
   real(wp), intent(in), optional :: c6(:, :)
   real(wp), intent(in), optional :: alpha(:)
   real(wp), intent(in), optional :: pairwise_energy2(:, :)
   real(wp), intent(in), optional :: pairwise_energy3(:, :)
   character(len=:), allocatable :: indent, version_string
   character(len=*), parameter :: jsonkey = "('""',a,'"":',1x)"
   real(wp), allocatable :: array(:)

   call get_dftd4_version(string=version_string)

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
      array = reshape(sigma, [size(sigma)])
      call write_json_array(unit, array, indent)
   end if
   if (present(gradient)) then
      write(unit, '(",")', advance='no')
      if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
      write(unit, jsonkey, advance='no') 'gradient'
      array = reshape(gradient, [size(gradient)])
      call write_json_array(unit, array, indent)
   end if
   if (present(hessian)) then
      write(unit, '(",")', advance='no')
      if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
      write(unit, jsonkey, advance='no') 'hessian'
      array = reshape(hessian, [size(hessian)])
      call write_json_array(unit, array, indent)
   end if
   if (present(cn)) then
      write(unit, '(",")', advance='no')
      if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
      write(unit, jsonkey, advance='no') 'coordination numbers'
      call write_json_array(unit, cn, indent)
   end if
   if (present(q)) then
      write(unit, '(",")', advance='no')
      if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
      write(unit, jsonkey, advance='no') 'partial charges'
      call write_json_array(unit, q, indent)
   end if
   if (present(c6)) then
      write(unit, '(",")', advance='no')
      if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
      write(unit, jsonkey, advance='no') 'c6 coefficients'
      array = reshape(c6, [size(c6)])
      call write_json_array(unit, array, indent)
   end if
   if (present(alpha)) then
      write(unit, '(",")', advance='no')
      if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
      write(unit, jsonkey, advance='no') 'polarizabilities'
      call write_json_array(unit, alpha, indent)
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
   if (allocated(indent)) write(unit, '(/)', advance='no')
   write(unit, '("}")')

end subroutine json_results


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


subroutine tagged_result(unit, energy, gradient, sigma, hessian)
   !DEC$ ATTRIBUTES DLLEXPORT :: tagged_result
   integer, intent(in) :: unit
   real(wp), intent(in), optional :: energy
   real(wp), intent(in), optional :: gradient(:, :)
   real(wp), intent(in), optional :: sigma(:, :)
   real(wp), intent(in), optional :: hessian(:, :, :, :)
   character(len=*), parameter :: tag_header = &
      & '(a,t20,":",a,":",i0,":",*(i0:,","))'

   if (present(energy)) then
      write(unit, tag_header) "energy", "real", 0
      write(unit, '(3es24.16)') energy
   end if
   if (present(gradient)) then
      write(unit, tag_header) "gradient", "real", 2, shape(gradient)
      write(unit, '(3es24.16)') gradient
   end if
   if (present(sigma)) then
      write(unit, tag_header) "virial", "real", 2, shape(sigma)
      write(unit, '(3es24.16)') sigma
   end if
   if (present(hessian)) then
      write(unit, tag_header) "hessian", "real", 4, shape(hessian)
      write(unit, '(3es24.16)') hessian
   end if

end subroutine tagged_result

end module dftd4_output
